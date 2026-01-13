/**
 * @file NearFieldExplosion.cpp
 * @brief Implementation of CRAM3D-style near-field explosion physics
 * 
 * See NearFieldExplosion.hpp for detailed documentation.
 */

#include "NearFieldExplosion.hpp"
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

namespace FSRM {

// =============================================================================
// NearFieldExplosionSolver Implementation
// =============================================================================

NearFieldExplosionSolver::NearFieldExplosionSolver() {
    // Default initialization
    source = UndergroundExplosionSource();
    eos = MieGruneisenEOS();
    strength = PressureDependentStrength();
    damage_model = DamageEvolutionModel();
}

void NearFieldExplosionSolver::setSource(const UndergroundExplosionSource& src) {
    source = src;
    source.computeRiseTime();
    shock_model.setFromSource(source);
    
    // Initialize seismic source parameters
    scalar_moment = source.scalarMoment();
    corner_frequency = 2.5 * std::pow(source.yield_kt, -1.0/3.0);
}

void NearFieldExplosionSolver::setEOS(const MieGruneisenEOS& e) {
    eos = e;
}

void NearFieldExplosionSolver::setStrengthModel(const PressureDependentStrength& s) {
    strength = s;
}

void NearFieldExplosionSolver::setDamageModel(const DamageEvolutionModel& d) {
    damage_model = d;
}

void NearFieldExplosionSolver::configure(const std::map<std::string, std::string>& config) {
    // Parse source parameters
    auto getDouble = [&config](const std::string& key, double def) {
        auto it = config.find(key);
        if (it != config.end() && !it->second.empty()) {
            try { return std::stod(it->second); } catch (...) {}
        }
        return def;
    };
    
    source.yield_kt = getDouble("yield_kt", 100.0);
    source.depth = getDouble("depth_of_burial", 1000.0);
    source.location[0] = getDouble("source_x", 0.0);
    source.location[1] = getDouble("source_y", 0.0);
    source.location[2] = getDouble("source_z", -source.depth);
    
    source.host_density = getDouble("host_density", 2650.0);
    source.host_vp = getDouble("p_wave_velocity", 5500.0);
    source.host_vs = getDouble("s_wave_velocity", 3200.0);
    source.host_porosity = getDouble("porosity", 0.01);
    
    // Parse EOS parameters
    eos.rho0 = source.host_density;
    eos.c0 = getDouble("bulk_sound_speed", 4500.0);
    eos.s = getDouble("hugoniot_slope", 1.4);
    eos.gamma0 = getDouble("gruneisen_gamma", 1.5);
    
    // Parse strength model
    strength.Y0 = getDouble("cohesion", 50.0e6);
    strength.A = getDouble("pressure_coefficient", 0.5);
    strength.n = getDouble("pressure_exponent", 0.5);
    strength.Y_max = getDouble("max_strength", 2.0e9);
    strength.Y_min = getDouble("residual_strength", 0.1e6);
    strength.softening_strain = getDouble("softening_strain", 0.05);
    
    // Parse damage model
    damage_model.strain_threshold = getDouble("damage_strain_threshold", 0.001);
    damage_model.strain_saturation = getDouble("damage_strain_saturation", 0.10);
    damage_model.pressure_threshold = getDouble("damage_pressure_threshold", 1.0e9);
    damage_model.pressure_saturation = getDouble("damage_pressure_saturation", 10.0e9);
    damage_model.damage_rate = getDouble("damage_rate", 10.0);
    damage_model.max_damage = getDouble("max_damage", 0.99);
    
    source.computeRiseTime();
    shock_model.setFromSource(source);
    scalar_moment = source.scalarMoment();
    corner_frequency = 2.5 * std::pow(source.yield_kt, -1.0/3.0);
}

void NearFieldExplosionSolver::initialize() {
    current_time = 0.0;
    cavity_radius_current = source.cavityRadius();
    cavity_velocity = 0.0;
    cavity_acceleration = 0.0;
    
    // Initialize energy
    energy.total = source.energyJoules();
    energy.kinetic = 0.0;
    energy.internal = energy.total;
    energy.plastic = 0.0;
    energy.seismic = 0.0;
    energy.fracture = 0.0;
    energy.cavity_pv = 0.0;
    
    initialized = true;
}

void NearFieldExplosionSolver::prestep(double time, double /*dt*/) {
    current_time = time;
}

void NearFieldExplosionSolver::step(double dt) {
    if (!initialized) {
        initialize();
    }
    
    // Update cavity expansion
    computeCavityExpansion(dt);
    
    current_time += dt;
}

void NearFieldExplosionSolver::poststep(double /*time*/, double /*dt*/) {
    // Post-step updates (e.g., energy accounting)
}

void NearFieldExplosionSolver::computeCavityExpansion(double dt) {
    // Simple cavity expansion model
    // Cavity radius approaches equilibrium exponentially
    
    double R_eq = source.cavityRadius();  // Equilibrium radius
    double tau = source.rise_time * 10.0;  // Expansion time constant
    
    if (current_time < tau) {
        // Expansion phase
        double expansion_rate = R_eq / tau * std::exp(-current_time / tau);
        cavity_velocity = expansion_rate;
        cavity_radius_current = R_eq * (1.0 - std::exp(-current_time / tau));
    } else {
        // Equilibrium
        cavity_velocity = 0.0;
        cavity_radius_current = R_eq;
    }
    
    // P*dV work
    double dR = cavity_velocity * dt;
    if (dR > 0) {
        double P_cavity = source.initialCavityPressure() * 
                         std::pow(R_eq / cavity_radius_current, 3);
        double dV = 4.0 * M_PI * cavity_radius_current * cavity_radius_current * dR;
        energy.cavity_pv += P_cavity * dV;
    }
}

void NearFieldExplosionSolver::computeShockPosition(double /*dt*/) {
    // Shock position tracking (for wave propagation coupling)
}

double NearFieldExplosionSolver::distanceToSource(const std::array<double, 3>& point) const {
    double dx = point[0] - source.location[0];
    double dy = point[1] - source.location[1];
    double dz = point[2] - source.location[2];
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

void NearFieldExplosionSolver::updateMaterialPoint(NearFieldMaterialState& state,
                                                   double dt, double time) {
    // Get distance from source
    std::array<double, 3> point = {0, 0, 0};  // Would be passed in real implementation
    double r = std::sqrt(state.strain[0] * state.strain[0] +
                        state.strain[1] * state.strain[1] +
                        state.strain[2] * state.strain[2]);
    if (r < 1e-10) r = 1.0;  // Avoid division by zero
    
    // Update pressure from EOS
    double mu = state.density / eos.rho0 - 1.0;
    state.pressure = eos.pressure(state.density, state.internal_energy);
    
    // Update peak pressure
    if (state.pressure > state.peak_pressure) {
        state.peak_pressure = state.pressure;
    }
    
    // Compute yield strength
    double Y = strength.yieldStrength(state.pressure, state.damage, 
                                      state.equivalent_plastic_strain);
    
    // Check for yielding
    double sigma_eq = state.vonMisesStress();
    
    if (sigma_eq > Y) {
        // Plastic flow
        double dgamma = (sigma_eq - Y) / (3.0 * eos.rho0 * eos.c0 * eos.c0);
        state.equivalent_plastic_strain += dgamma;
        
        // Update plastic strain tensor (radial return)
        double scale = Y / sigma_eq;
        for (int i = 0; i < 6; ++i) {
            double dev = (i < 3) ? state.stress[i] - state.pressure : state.stress[i];
            state.plastic_strain[i] += dev * (1.0 - scale) / (2.0 * eos.rho0 * eos.c0 * eos.c0);
        }
        
        // Plastic work
        double plastic_work = sigma_eq * dgamma;
        energy.plastic += plastic_work * dt;
    }
    
    // Update damage
    double strain_eq = state.deviatoricStrainMagnitude();
    double dD = damage_model.computeDamageIncrement(state.damage, strain_eq,
                                                     state.pressure, dt);
    state.damage += dD;
    state.damage = std::min(damage_model.max_damage, std::max(0.0, state.damage));
    
    // Update damage state
    state.damage_state = damage_model.getDamageState(state.damage, state.temperature);
    
    // Check failure mode
    if (state.damage > 0.1) {
        if (state.pressure > ExplosionPhysics::CRUSHING_PRESSURE) {
            state.failure_mode = FailureMode::COMPRESSIVE;
        } else if (state.pressure < -strength.tensileStrength(0)) {
            state.failure_mode = FailureMode::TENSILE;
            state.is_spalled = true;
        } else if (sigma_eq > Y * 0.9) {
            state.failure_mode = FailureMode::SHEAR;
        } else {
            state.failure_mode = FailureMode::MIXED;
        }
    }
    
    // Update temperature from internal energy
    // Simplified: T = T0 + E / (ρ * c_v)
    double cv = 800.0;  // Specific heat capacity (J/kg·K)
    state.temperature = 300.0 + state.internal_energy / cv;
    if (state.temperature > state.peak_temperature) {
        state.peak_temperature = state.temperature;
    }
    
    // Fracture energy
    if (dD > 0) {
        double G_c = 100.0;  // Fracture energy (J/m²)
        energy.fracture += G_c * dD * 1.0;  // Approximate
    }
}

void NearFieldExplosionSolver::computeSourceFunction(double time,
                                                     std::array<double, 6>& moment_tensor,
                                                     double& moment_rate) {
    // Use RDP model
    double tau = 1.0 / (2.0 * M_PI * corner_frequency);
    double t_norm = time / tau;
    
    // Brune source with overshoot
    double psi, psi_dot;
    if (time <= 0) {
        psi = 0.0;
        psi_dot = 0.0;
    } else {
        double overshoot = source.overshoot;
        psi = (1.0 - (1.0 + t_norm) * std::exp(-t_norm)) * overshoot;
        psi_dot = t_norm * std::exp(-t_norm) * overshoot / tau;
    }
    
    // Steady-state RDP
    double psi_inf = 1.5e5 * std::pow(source.yield_kt, 0.8);  // m³
    
    // Scale by medium properties
    double K = source.host_density * source.host_vp * source.host_vp;
    double M0 = 4.0 * M_PI * K * psi_inf * psi;
    double M0_dot = 4.0 * M_PI * K * psi_inf * psi_dot;
    
    // Moment tensor (predominantly isotropic)
    double iso_frac = 0.70;    // 70% isotropic
    double clvd_frac = 0.25;   // 25% CLVD
    // double dc_frac = 0.05;  // 5% double-couple (unused directly)
    
    double M_iso = iso_frac * M0 / 3.0;
    double M_clvd = clvd_frac * M0;
    
    moment_tensor[0] = M_iso + M_clvd / 3.0;       // Mxx
    moment_tensor[1] = M_iso + M_clvd / 3.0;       // Myy
    moment_tensor[2] = M_iso - 2.0 * M_clvd / 3.0; // Mzz
    moment_tensor[3] = 0.0;  // Mxy
    moment_tensor[4] = 0.0;  // Mxz
    moment_tensor[5] = 0.0;  // Myz
    
    moment_rate = M0_dot;
}

bool NearFieldExplosionSolver::isInCavity(const std::array<double, 3>& point, 
                                          double /*time*/) const {
    double r = distanceToSource(point);
    return r < cavity_radius_current;
}

bool NearFieldExplosionSolver::isInCrushedZone(const std::array<double, 3>& point) const {
    double r = distanceToSource(point);
    return r >= cavity_radius_current && r < source.crushedZoneRadius();
}

bool NearFieldExplosionSolver::isInFracturedZone(const std::array<double, 3>& point) const {
    double r = distanceToSource(point);
    return r >= source.crushedZoneRadius() && r < source.fracturedZoneRadius();
}

bool NearFieldExplosionSolver::isSpalled(const std::array<double, 3>& /*point*/) const {
    // Would check spall map in real implementation
    return false;
}

double NearFieldExplosionSolver::getPeakPressure(const std::array<double, 3>& point) const {
    double r = distanceToSource(point);
    return shock_model.peakPressure(r);
}

double NearFieldExplosionSolver::getPeakTemperature(const std::array<double, 3>& point) const {
    // Approximate temperature from peak pressure
    double P = getPeakPressure(point);
    double T0 = 300.0;
    double gamma = eos.gamma0;
    // Rough estimate: T/T0 = (P/P0)^(γ-1)/γ
    double T = T0 * std::pow(1.0 + P / (eos.rho0 * eos.c0 * eos.c0), gamma);
    return T;
}

double NearFieldExplosionSolver::getDamage(const std::array<double, 3>& point) const {
    double r = distanceToSource(point);
    double R_c = source.cavityRadius();
    
    // Simple radial damage profile
    if (r < R_c) return 1.0;  // Cavity
    
    double R_cr = source.crushedZoneRadius();
    double R_fr = source.fracturedZoneRadius();
    double R_dm = source.damagedZoneRadius();
    
    if (r < R_cr) {
        // Crushed zone: high damage
        return 0.9 + 0.09 * (1.0 - (r - R_c) / (R_cr - R_c));
    } else if (r < R_fr) {
        // Fractured zone
        return 0.4 + 0.5 * (1.0 - (r - R_cr) / (R_fr - R_cr));
    } else if (r < R_dm) {
        // Damaged zone
        return 0.1 + 0.3 * (1.0 - (r - R_fr) / (R_dm - R_fr));
    }
    
    return 0.0;  // Undamaged
}

DamageState NearFieldExplosionSolver::getDamageState(const std::array<double, 3>& point) const {
    double D = getDamage(point);
    double T = getPeakTemperature(point);
    return damage_model.getDamageState(D, T);
}

double NearFieldExplosionSolver::getCavityRadius(double /*time*/) const {
    return cavity_radius_current;
}

double NearFieldExplosionSolver::getCrushedZoneRadius() const {
    return source.crushedZoneRadius();
}

double NearFieldExplosionSolver::getFracturedZoneRadius() const {
    return source.fracturedZoneRadius();
}

double NearFieldExplosionSolver::getScalarMoment() const {
    return scalar_moment;
}

double NearFieldExplosionSolver::getCornerFrequency() const {
    return corner_frequency;
}

void NearFieldExplosionSolver::getMomentTensor(double time, 
                                               std::array<double, 6>& M) const {
    double moment_rate;
    std::array<double, 6> M_temp;
    const_cast<NearFieldExplosionSolver*>(this)->computeSourceFunction(time, M_temp, moment_rate);
    M = M_temp;
}

NearFieldExplosionSolver::EnergyPartition 
NearFieldExplosionSolver::getEnergyPartition() const {
    EnergyPartition e = energy;
    
    // Estimate seismic energy from moment
    // E_s ≈ M0 / (2 * μ) where μ is shear modulus
    double mu = source.host_density * source.host_vs * source.host_vs;
    e.seismic = scalar_moment / (2.0 * mu) * 0.0001;  // ~0.01% seismic efficiency
    
    // Kinetic energy (approximate)
    e.kinetic = 0.1 * energy.total * std::exp(-current_time / source.rise_time);
    
    // Internal = total - kinetic - plastic - seismic - fracture
    e.internal = energy.total - e.kinetic - e.plastic - e.seismic - e.fracture;
    if (e.internal < 0) e.internal = 0;
    
    return e;
}

void NearFieldExplosionSolver::writeOutput(const std::string& filename, double time) {
    std::ofstream out(filename);
    
    out << "# Near-field explosion output at t = " << time << " s\n";
    out << "# Yield: " << source.yield_kt << " kt\n";
    out << "# Depth: " << source.depth << " m\n";
    out << "#\n";
    
    // Zone radii
    out << "# Zone Radii:\n";
    out << "#   Cavity: " << cavity_radius_current << " m\n";
    out << "#   Crushed zone: " << source.crushedZoneRadius() << " m\n";
    out << "#   Fractured zone: " << source.fracturedZoneRadius() << " m\n";
    out << "#   Damaged zone: " << source.damagedZoneRadius() << " m\n";
    out << "#\n";
    
    // Seismic source
    out << "# Seismic Source:\n";
    out << "#   Scalar moment: " << scalar_moment << " N·m\n";
    out << "#   Moment magnitude: " << ((std::log10(scalar_moment) - 9.1) / 1.5) << "\n";
    out << "#   Body wave magnitude: " << source.bodyWaveMagnitude() << "\n";
    out << "#   Corner frequency: " << corner_frequency << " Hz\n";
    out << "#\n";
    
    // Energy
    auto e = getEnergyPartition();
    out << "# Energy Partition:\n";
    out << "#   Total: " << e.total << " J\n";
    out << "#   Kinetic: " << e.kinetic << " J (" << 100*e.kinetic/e.total << "%)\n";
    out << "#   Internal: " << e.internal << " J (" << 100*e.internal/e.total << "%)\n";
    out << "#   Plastic: " << e.plastic << " J (" << 100*e.plastic/e.total << "%)\n";
    out << "#   Seismic: " << e.seismic << " J (" << 100*e.seismic/e.total << "%)\n";
    out << "#   Fracture: " << e.fracture << " J (" << 100*e.fracture/e.total << "%)\n";
    out << "#\n";
    
    // Radial profiles
    out << "# Radial profiles (radius, pressure, temperature, damage)\n";
    out << "# r(m)  P_peak(Pa)  T_peak(K)  damage\n";
    
    double R_max = 5.0 * source.damagedZoneRadius();
    int n_points = 200;
    double dr = R_max / n_points;
    
    for (int i = 1; i <= n_points; ++i) {
        double r = i * dr;
        std::array<double, 3> point = {source.location[0] + r, 
                                       source.location[1],
                                       source.location[2]};
        
        double P = getPeakPressure(point);
        double T = getPeakTemperature(point);
        double D = getDamage(point);
        
        out << std::scientific << std::setprecision(6)
            << r << "  " << P << "  " << T << "  " << D << "\n";
    }
}

void NearFieldExplosionSolver::writeZoneFile(const std::string& filename) {
    std::ofstream out(filename);
    
    out << "# Near-field damage zone definitions\n";
    out << "# All dimensions in meters, centered at source location\n";
    out << "#\n";
    out << "source_x = " << source.location[0] << "\n";
    out << "source_y = " << source.location[1] << "\n";
    out << "source_z = " << source.location[2] << "\n";
    out << "depth = " << source.depth << "\n";
    out << "yield_kt = " << source.yield_kt << "\n";
    out << "\n";
    out << "# Zone radii\n";
    out << "cavity_radius = " << source.cavityRadius() << "\n";
    out << "crushed_zone_radius = " << source.crushedZoneRadius() << "\n";
    out << "fractured_zone_radius = " << source.fracturedZoneRadius() << "\n";
    out << "damaged_zone_radius = " << source.damagedZoneRadius() << "\n";
    out << "\n";
    out << "# Zone volumes (spherical)\n";
    double R_c = source.cavityRadius();
    double R_cr = source.crushedZoneRadius();
    double R_fr = source.fracturedZoneRadius();
    double R_dm = source.damagedZoneRadius();
    
    out << "cavity_volume = " << (4.0/3.0 * M_PI * R_c * R_c * R_c) << "\n";
    out << "crushed_volume = " << (4.0/3.0 * M_PI * (R_cr*R_cr*R_cr - R_c*R_c*R_c)) << "\n";
    out << "fractured_volume = " << (4.0/3.0 * M_PI * (R_fr*R_fr*R_fr - R_cr*R_cr*R_cr)) << "\n";
    out << "damaged_volume = " << (4.0/3.0 * M_PI * (R_dm*R_dm*R_dm - R_fr*R_fr*R_fr)) << "\n";
}

// =============================================================================
// Analysis Functions
// =============================================================================

NearFieldAnalysisResults analyzeNearField(const NearFieldExplosionSolver& solver) {
    NearFieldAnalysisResults results;
    
    const auto& src = solver.getSource();
    
    // Zone dimensions
    results.cavity_radius = src.cavityRadius();
    results.crushed_zone_radius = src.crushedZoneRadius();
    results.fractured_zone_radius = src.fracturedZoneRadius();
    results.damaged_zone_radius = src.damagedZoneRadius();
    
    // Zone volumes
    double R_c = results.cavity_radius;
    double R_cr = results.crushed_zone_radius;
    double R_fr = results.fractured_zone_radius;
    double R_dm = results.damaged_zone_radius;
    
    results.cavity_volume = 4.0/3.0 * M_PI * R_c * R_c * R_c;
    results.crushed_volume = 4.0/3.0 * M_PI * (R_cr*R_cr*R_cr - R_c*R_c*R_c);
    results.fractured_volume = 4.0/3.0 * M_PI * (R_fr*R_fr*R_fr - R_cr*R_cr*R_cr);
    results.damaged_volume = 4.0/3.0 * M_PI * (R_dm*R_dm*R_dm - R_fr*R_fr*R_fr);
    results.total_damage_volume = results.crushed_volume + results.fractured_volume + 
                                  results.damaged_volume;
    
    // Energy
    auto energy = solver.getEnergyPartition();
    results.total_energy = energy.total;
    results.seismic_energy = energy.seismic;
    results.seismic_efficiency = energy.seismic / energy.total;
    
    // Seismic source
    results.scalar_moment = solver.getScalarMoment();
    results.moment_magnitude = (std::log10(results.scalar_moment) - 9.1) / 1.5;
    results.body_wave_magnitude = src.bodyWaveMagnitude();
    results.corner_frequency = solver.getCornerFrequency();
    
    // Radial profiles
    int n_radii = 50;
    double R_max = 5.0 * R_dm;
    
    results.radii.resize(n_radii);
    results.peak_pressures.resize(n_radii);
    results.peak_velocities.resize(n_radii);
    results.peak_accelerations.resize(n_radii);
    results.damage_values.resize(n_radii);
    
    const auto& shock = solver.getShockModel();
    
    for (int i = 0; i < n_radii; ++i) {
        double r = (i + 1) * R_max / n_radii;
        results.radii[i] = r;
        
        std::array<double, 3> point = {src.location[0] + r, 
                                       src.location[1], 
                                       src.location[2]};
        
        results.peak_pressures[i] = shock.peakPressure(r);
        results.peak_velocities[i] = shock.peakParticleVelocity(r, src.host_density, 
                                                                 src.host_vp);
        results.damage_values[i] = solver.getDamage(point);
        
        // Estimate peak acceleration
        double tau_pos = shock.positivePhaseDuration(r);
        if (tau_pos > 0) {
            results.peak_accelerations[i] = results.peak_velocities[i] / tau_pos;
        }
    }
    
    // Average damage in damaged zone
    double sum_damage = 0.0;
    int count = 0;
    for (int i = 0; i < n_radii; ++i) {
        if (results.radii[i] < R_dm) {
            sum_damage += results.damage_values[i];
            ++count;
        }
    }
    results.average_damage = count > 0 ? sum_damage / count : 0.0;
    
    return results;
}

void writeNearFieldPlotScript(const std::string& base_name,
                              const NearFieldAnalysisResults& results) {
    std::ofstream out(base_name + "_plot.gp");
    
    out << R"(# Gnuplot script for near-field explosion visualization
set terminal pngcairo enhanced size 1600,1200 font 'Arial,12'
set output ')" << base_name << R"(_results.png'

set multiplot layout 2,3 title 'Underground Explosion Near-Field Effects' font ',16'

# Color palette for zones
set palette defined (0 'green', 0.3 'yellow', 0.6 'orange', 0.9 'red', 1.0 'black')

# 1. Peak Pressure vs Distance
set xlabel 'Distance from Source (m)'
set ylabel 'Peak Pressure (GPa)'
set title 'Shock Pressure Attenuation'
set logscale xy
set grid
set key top right
plot ')" << base_name << R"(_radial.dat' u 1:($2/1e9) w l lw 2 title 'Peak Pressure', \
     )" << results.crushed_zone_radius << R"( notitle w impulse lc 'red' dt 2, \
     )" << results.fractured_zone_radius << R"( notitle w impulse lc 'orange' dt 2, \
     )" << results.damaged_zone_radius << R"( notitle w impulse lc 'yellow' dt 2
unset logscale

# 2. Peak Velocity vs Distance
set ylabel 'Peak Particle Velocity (m/s)'
set title 'Ground Motion Attenuation'
set logscale xy
plot ')" << base_name << R"(_radial.dat' u 1:3 w l lw 2 lc 'blue' title 'Peak Velocity'
unset logscale

# 3. Damage Profile
set ylabel 'Damage [-]'
set title 'Radial Damage Distribution'
set logscale x
set yrange [0:1]
plot ')" << base_name << R"(_radial.dat' u 1:4 w filledcurves x1 lc rgb '#FF6666' title 'Damage'
unset logscale
set yrange [*:*]

# 4. Zone Diagram (cross-section)
set title 'Damage Zone Cross-Section'
set xlabel 'X Distance (m)'
set ylabel 'Z Depth (m)'
set size ratio -1
unset key
set xrange [-)" << 2*results.damaged_zone_radius << R"(:)" << 2*results.damaged_zone_radius << R"(]
set yrange [-)" << 2*results.damaged_zone_radius << R"(:)" << 0.5*results.damaged_zone_radius << R"(]

# Draw zones as circles
set object 1 circle at 0,0 size )" << results.cavity_radius << R"( fc rgb 'black' fs solid front
set object 2 circle at 0,0 size )" << results.crushed_zone_radius << R"( fc rgb 'red' fs solid back
set object 3 circle at 0,0 size )" << results.fractured_zone_radius << R"( fc rgb 'orange' fs solid back
set object 4 circle at 0,0 size )" << results.damaged_zone_radius << R"( fc rgb 'yellow' fs solid back

# Surface line
set arrow from -)" << 2*results.damaged_zone_radius << R"(,0 to )" << 2*results.damaged_zone_radius << R"(,0 nohead lw 2 lc 'brown'

plot 1/0 notitle

unset object 1
unset object 2
unset object 3
unset object 4
unset arrow

# 5. Source Time Function
set title 'Seismic Source Time Function'
set xlabel 'Time (s)'
set ylabel 'Normalized Moment Rate'
set xrange [0:5]
set size noratio
set key top right

# Approximate Brune source
corner_freq = )" << results.corner_frequency << R"(
tau = 1.0 / (2.0 * pi * corner_freq)
brune(t) = (t > 0) ? t/tau * exp(-t/tau) : 0
plot brune(x) w l lw 2 title sprintf('f_c = %.2f Hz', corner_freq)

# 6. Summary text
unset border
unset tics
unset xlabel
unset ylabel
set title 'Summary'
set label 1 sprintf('Yield: %.1f kt', )" << results.body_wave_magnitude << R"() at screen 0.72, 0.35
set label 2 sprintf('Cavity Radius: %.1f m', )" << results.cavity_radius << R"() at screen 0.72, 0.32
set label 3 sprintf('Crushed Zone: %.1f m', )" << results.crushed_zone_radius << R"() at screen 0.72, 0.29
set label 4 sprintf('Fractured Zone: %.1f m', )" << results.fractured_zone_radius << R"() at screen 0.72, 0.26
set label 5 sprintf('M_w = %.2f', )" << results.moment_magnitude << R"() at screen 0.72, 0.23
set label 6 sprintf('m_b = %.2f', )" << results.body_wave_magnitude << R"() at screen 0.72, 0.20
set label 7 sprintf('Seismic Efficiency: %.4f%%', )" << 100*results.seismic_efficiency << R"() at screen 0.72, 0.17
plot 1/0 notitle

unset multiplot
print 'Plot saved to: )" << base_name << R"(_results.png'
)";

    // Write radial data file
    std::ofstream dat(base_name + "_radial.dat");
    dat << "# Radial profiles from near-field explosion\n";
    dat << "# radius(m) pressure(Pa) velocity(m/s) damage\n";
    
    for (size_t i = 0; i < results.radii.size(); ++i) {
        dat << std::scientific << std::setprecision(6)
            << results.radii[i] << " "
            << results.peak_pressures[i] << " "
            << results.peak_velocities[i] << " "
            << results.damage_values[i] << "\n";
    }
}

// =============================================================================
// Configuration Parser
// =============================================================================

bool NearFieldExplosionConfig::parse(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open config file: " << filename << std::endl;
        return false;
    }
    
    std::map<std::string, std::string> config;
    std::string line;
    
    while (std::getline(file, line)) {
        // Skip comments and empty lines
        if (line.empty() || line[0] == '#' || line[0] == ';') continue;
        
        // Remove inline comments
        size_t comment_pos = line.find('#');
        if (comment_pos != std::string::npos) {
            line = line.substr(0, comment_pos);
        }
        
        // Find key=value
        size_t eq_pos = line.find('=');
        if (eq_pos == std::string::npos) continue;
        
        std::string key = line.substr(0, eq_pos);
        std::string value = line.substr(eq_pos + 1);
        
        // Trim whitespace
        auto trim = [](std::string& s) {
            while (!s.empty() && (s.front() == ' ' || s.front() == '\t')) s.erase(0, 1);
            while (!s.empty() && (s.back() == ' ' || s.back() == '\t')) s.pop_back();
        };
        trim(key);
        trim(value);
        
        config[key] = value;
    }
    
    return parse(config);
}

bool NearFieldExplosionConfig::parse(const std::map<std::string, std::string>& config) {
    auto getDouble = [&config](const std::string& key, double def) {
        auto it = config.find(key);
        if (it != config.end() && !it->second.empty()) {
            try { return std::stod(it->second); } catch (...) {}
        }
        return def;
    };
    
    // Source parameters
    source.yield_kt = getDouble("yield_kt", 100.0);
    source.depth = getDouble("depth_of_burial", 1000.0);
    source.location[0] = getDouble("source_x", 0.0);
    source.location[1] = getDouble("source_y", 0.0);
    source.location[2] = -source.depth;
    
    source.host_density = getDouble("host_density", 2650.0);
    source.host_vp = getDouble("p_wave_velocity", 5500.0);
    source.host_vs = getDouble("s_wave_velocity", 3200.0);
    source.host_porosity = getDouble("porosity", 0.01);
    source.overshoot = getDouble("overshoot_factor", 1.2);
    
    source.computeRiseTime();
    
    // EOS parameters
    eos.rho0 = source.host_density;
    eos.c0 = getDouble("bulk_sound_speed", 4500.0);
    eos.s = getDouble("hugoniot_slope", 1.4);
    eos.gamma0 = getDouble("gruneisen_gamma", 1.5);
    eos.a = getDouble("gruneisen_a", 0.5);
    
    // Strength parameters
    strength.Y0 = getDouble("cohesion", 50.0e6);
    strength.A = getDouble("pressure_coefficient", 0.5);
    strength.n = getDouble("pressure_exponent", 0.5);
    strength.B = getDouble("tension_coefficient", 5.0);
    strength.Y_max = getDouble("max_strength", 2.0e9);
    strength.Y_min = getDouble("residual_strength", 0.1e6);
    strength.softening_strain = getDouble("softening_strain", 0.05);
    strength.residual_strength_ratio = getDouble("residual_strength_ratio", 0.1);
    
    // Damage parameters
    damage.strain_threshold = getDouble("damage_strain_threshold", 0.001);
    damage.strain_saturation = getDouble("damage_strain_saturation", 0.10);
    damage.pressure_threshold = getDouble("damage_pressure_threshold", 1.0e9);
    damage.pressure_saturation = getDouble("damage_pressure_saturation", 10.0e9);
    damage.tensile_threshold = getDouble("tensile_damage_threshold", 5.0e6);
    damage.damage_rate = getDouble("damage_rate", 10.0);
    damage.max_damage = getDouble("max_damage", 0.99);
    damage.healing_enabled = (getDouble("healing_enabled", 0.0) > 0.5);
    damage.healing_rate = getDouble("healing_rate", 0.01);
    
    configured = true;
    return true;
}

std::unique_ptr<NearFieldExplosionSolver> 
NearFieldExplosionConfig::createSolver() const {
    auto solver = std::make_unique<NearFieldExplosionSolver>();
    
    solver->setSource(source);
    solver->setEOS(eos);
    solver->setStrengthModel(strength);
    solver->setDamageModel(damage);
    solver->initialize();
    
    return solver;
}

} // namespace FSRM
