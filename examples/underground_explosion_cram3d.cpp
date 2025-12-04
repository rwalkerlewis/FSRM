/**
 * @file underground_explosion_cram3d.cpp
 * @brief Comprehensive underground explosion example with CRAM3D-style near-field effects
 * 
 * This example demonstrates near-field deformation modeling for underground explosions,
 * implementing physics from the CRAM3D (Computational Rock Analysis Model - 3D) code:
 * https://apps.dtic.mil/sti/tr/pdf/AD1064154.pdf
 * 
 * Features:
 * 1. Near-field damage zones (crushed, fractured, disturbed)
 * 2. Cavity formation and growth
 * 3. Shock wave propagation and attenuation
 * 4. Pressure-dependent strength model
 * 5. Damage evolution and material degradation
 * 6. Spall (tensile failure) modeling
 * 7. Seismic source generation (RDP model)
 * 8. Energy partitioning
 * 
 * Based on Nevada Test Site empirical relations and theoretical models.
 * 
 * Usage:
 *   ./underground_explosion_cram3d [config_file] [yield_kt]
 * 
 * Examples:
 *   ./underground_explosion_cram3d                           # Default 100 kt
 *   ./underground_explosion_cram3d config.txt                # With config
 *   ./underground_explosion_cram3d "" 150                    # 150 kt yield
 */

#include "NearFieldExplosion.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>
#include <chrono>
#include <memory>

using namespace FSRM;

// =============================================================================
// Helper Functions
// =============================================================================

// Print a separator line
void printSeparator(char c = '=', int length = 78) {
    std::cout << std::string(length, c) << std::endl;
}

// Format number with SI prefix
std::string formatWithPrefix(double value, const std::string& unit) {
    const char* prefixes[] = {"", "k", "M", "G", "T", "P"};
    int idx = 0;
    while (std::abs(value) >= 1000.0 && idx < 5) {
        value /= 1000.0;
        ++idx;
    }
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2) << value << " " << prefixes[idx] << unit;
    return oss.str();
}

// Get damage zone color for visualization
std::string getDamageZoneColor(DamageState state) {
    switch (state) {
        case DamageState::VAPORIZED: return "magenta";
        case DamageState::MELTED: return "red";
        case DamageState::CRUSHED: return "orange";
        case DamageState::COMMINUTED: return "yellow";
        case DamageState::FRACTURED: return "green";
        case DamageState::MICRO_FRACTURED: return "cyan";
        case DamageState::INTACT: return "blue";
        default: return "white";
    }
}

// Get damage state name
std::string getDamageStateName(DamageState state) {
    switch (state) {
        case DamageState::VAPORIZED: return "Vaporized";
        case DamageState::MELTED: return "Melted";
        case DamageState::CRUSHED: return "Crushed";
        case DamageState::COMMINUTED: return "Comminuted";
        case DamageState::FRACTURED: return "Fractured";
        case DamageState::MICRO_FRACTURED: return "Micro-fractured";
        case DamageState::INTACT: return "Intact";
        default: return "Unknown";
    }
}

// =============================================================================
// Output Functions
// =============================================================================

/**
 * @brief Write comprehensive time-history output
 */
void writeTimeHistory(const std::string& filename,
                     const std::vector<double>& times,
                     const std::vector<double>& cavity_radii,
                     const std::vector<double>& cavity_pressures,
                     const std::vector<std::array<double, 6>>& moment_tensors,
                     const std::vector<double>& moment_rates) {
    std::ofstream out(filename);
    
    out << "# Time history of underground explosion\n";
    out << "# time(s) cavity_r(m) cavity_P(Pa) Mxx(N-m) Myy(N-m) Mzz(N-m) Mdot(N-m/s)\n";
    
    for (size_t i = 0; i < times.size(); ++i) {
        out << std::scientific << std::setprecision(6)
            << times[i] << " "
            << cavity_radii[i] << " "
            << cavity_pressures[i] << " "
            << moment_tensors[i][0] << " "
            << moment_tensors[i][1] << " "
            << moment_tensors[i][2] << " "
            << moment_rates[i] << "\n";
    }
}

/**
 * @brief Write ground motion at various distances
 */
void writeGroundMotion(const std::string& filename,
                       const ShockAttenuationModel& shock,
                       const UndergroundExplosionSource& source,
                       double max_distance) {
    std::ofstream out(filename);
    
    out << "# Ground motion vs distance for underground explosion\n";
    out << "# Yield: " << source.yield_kt << " kt\n";
    out << "# Depth: " << source.depth << " m\n";
    out << "#\n";
    out << "# distance(m) peak_P(Pa) peak_v(m/s) peak_a(m/s2) arrival(s) duration(s)\n";
    
    int n_points = 100;
    double dr = max_distance / n_points;
    
    for (int i = 1; i <= n_points; ++i) {
        double r = i * dr;
        double P = shock.peakPressure(r);
        double v = shock.peakParticleVelocity(r, source.host_density, source.host_vp);
        double tau = shock.positivePhaseDuration(r);
        double a = (tau > 0) ? v / tau : 0;
        double t_arr = shock.arrivalTime(r);
        
        out << std::scientific << std::setprecision(6)
            << r << " " << P << " " << v << " " << a << " "
            << t_arr << " " << tau << "\n";
    }
}

/**
 * @brief Write damage zone visualization data (2D cross-section)
 */
void writeDamageZoneMap(const std::string& filename,
                        const NearFieldExplosionSolver& solver,
                        double extent) {
    std::ofstream out(filename);
    
    const auto& src = solver.getSource();
    
    out << "# Damage zone map (XZ cross-section at Y=0)\n";
    out << "# x(m) z(m) damage pressure(Pa) temperature(K) damage_state\n";
    
    int nx = 200;
    int nz = 200;
    double dx = 2.0 * extent / nx;
    double dz = 2.0 * extent / nz;
    
    for (int j = 0; j < nz; ++j) {
        double z = src.location[2] - extent + j * dz;
        for (int i = 0; i < nx; ++i) {
            double x = src.location[0] - extent + i * dx;
            
            std::array<double, 3> point = {x, src.location[1], z};
            
            double D = solver.getDamage(point);
            double P = solver.getPeakPressure(point);
            double T = solver.getPeakTemperature(point);
            int state = static_cast<int>(solver.getDamageState(point));
            
            out << std::scientific << std::setprecision(6)
                << x << " " << z << " " << D << " " << P << " " << T << " "
                << state << "\n";
        }
        out << "\n";  // Blank line for gnuplot pm3d
    }
}

/**
 * @brief Write seismic source time function
 */
void writeSourceTimeFunction(const std::string& filename,
                             const RDPSeismicSource& rdp,
                             double duration) {
    std::ofstream out(filename);
    
    out << "# Seismic source time function (Reduced Displacement Potential)\n";
    out << "# Yield: " << rdp.yield_kt << " kt\n";
    out << "# Corner frequency: " << rdp.corner_frequency << " Hz\n";
    out << "# Scalar moment: " << rdp.scalarMoment() << " N-m\n";
    out << "#\n";
    out << "# time(s) psi(m3) psi_dot(m3/s) M0(N-m) M0_dot(N-m/s)\n";
    
    double K = rdp.density * rdp.p_velocity * rdp.p_velocity;
    
    int n_points = 500;
    double dt = duration / n_points;
    
    for (int i = 0; i <= n_points; ++i) {
        double t = i * dt;
        double psi = rdp.psi(t);
        double psi_dot = rdp.psiDot(t);
        double M0 = 4.0 * M_PI * K * psi;
        double M0_dot = 4.0 * M_PI * K * psi_dot;
        
        out << std::scientific << std::setprecision(6)
            << t << " " << psi << " " << psi_dot << " " << M0 << " " << M0_dot << "\n";
    }
}

/**
 * @brief Write comprehensive gnuplot visualization script
 */
void writeVisualizationScript(const std::string& base_name,
                              const NearFieldAnalysisResults& results,
                              const UndergroundExplosionSource& source) {
    std::ofstream out(base_name + "_visualize.gp");
    
    out << R"(#!/usr/bin/gnuplot
# Visualization script for underground explosion CRAM3D simulation
# Generated automatically - edit as needed

set terminal pngcairo enhanced size 2400,1800 font 'Arial,14'
set output ')" << base_name << R"(_complete_results.png'

# Set up multiplot
set multiplot layout 3,3 title 'Underground Explosion Near-Field Analysis\n)" 
    << source.yield_kt << " kt at " << source.depth << R"( m depth' font ',18'

#------------------------------------------------------------------------------
# Plot 1: Shock Pressure Attenuation (log-log)
#------------------------------------------------------------------------------
set title 'Shock Pressure Attenuation'
set xlabel 'Distance from Source (m)'
set ylabel 'Peak Pressure (GPa)'
set logscale xy
set grid
set key top right

# Mark zone boundaries
R_cavity = )" << results.cavity_radius << R"(
R_crushed = )" << results.crushed_zone_radius << R"(
R_fractured = )" << results.fractured_zone_radius << R"(
R_damaged = )" << results.damaged_zone_radius << R"(

set arrow from R_cavity, graph 0 to R_cavity, graph 1 nohead dt 2 lc 'black' lw 2
set arrow from R_crushed, graph 0 to R_crushed, graph 1 nohead dt 3 lc 'red' lw 2
set arrow from R_fractured, graph 0 to R_fractured, graph 1 nohead dt 3 lc 'orange' lw 2
set arrow from R_damaged, graph 0 to R_damaged, graph 1 nohead dt 3 lc 'yellow' lw 2

plot ')" << base_name << R"(_ground_motion.dat' u 1:($2/1e9) w l lw 3 lc rgb '#0000AA' title 'Peak Pressure', \
     5 w l dt 2 lc 'gray' title 'Crushing (5 GPa)', \
     0.5 w l dt 3 lc 'gray' title 'Fracture (0.5 GPa)'

unset arrow
unset logscale

#------------------------------------------------------------------------------
# Plot 2: Particle Velocity
#------------------------------------------------------------------------------
set title 'Peak Particle Velocity'
set ylabel 'Velocity (m/s)'
set logscale xy

plot ')" << base_name << R"(_ground_motion.dat' u 1:3 w l lw 3 lc rgb '#00AA00' title 'Peak Velocity', \
     1.0 w l dt 2 lc 'gray' title 'Spall threshold (1 m/s)'

unset logscale

#------------------------------------------------------------------------------
# Plot 3: Damage Distribution
#------------------------------------------------------------------------------
set title 'Radial Damage Distribution'
set ylabel 'Damage Parameter D [-]'
set yrange [0:1.1]
set logscale x

# Fill regions
set style fill solid 0.3 border

plot ')" << base_name << R"(_radial.dat' u 1:4 w filledcurves x1 lc rgb '#FF6666' title 'Damage', \
     ')" << base_name << R"(_radial.dat' u 1:4 w l lw 2 lc 'red' notitle

unset logscale
set yrange [*:*]

#------------------------------------------------------------------------------
# Plot 4: Damage Zone Cross-Section (Color map)
#------------------------------------------------------------------------------
set title 'Damage Zone Cross-Section (XZ plane)'
set xlabel 'X Distance (m)'
set ylabel 'Z Depth (m)'
set size ratio -1

# Color palette for damage
set palette defined (0 '#0000FF', 0.2 '#00FFFF', 0.4 '#00FF00', \
                     0.6 '#FFFF00', 0.8 '#FF8800', 1.0 '#FF0000')
set cbrange [0:1]
set cblabel 'Damage'

# Draw surface
set arrow from graph 0, first 0 to graph 1, first 0 nohead lw 3 lc 'brown'

plot ')" << base_name << R"(_damage_map.dat' u 1:2:3 w image notitle

unset arrow
set size noratio

#------------------------------------------------------------------------------
# Plot 5: Source Time Function
#------------------------------------------------------------------------------
set title 'Seismic Source Time Function'
set xlabel 'Time (s)'
set ylabel 'Moment Rate (N·m/s)'
unset logscale

set y2label 'Cumulative Moment (N·m)'
set y2tics
set ytics nomirror

plot ')" << base_name << R"(_source.dat' u 1:5 w l lw 3 lc 'blue' axes x1y1 title 'Moment Rate', \
     ')" << base_name << R"(_source.dat' u 1:4 w l lw 2 lc 'red' axes x1y2 title 'Moment'

unset y2label
unset y2tics

#------------------------------------------------------------------------------
# Plot 6: Cavity and Zone Evolution 
#------------------------------------------------------------------------------
set title 'Cavity and Zone Radii'
set xlabel 'Time (s)'
set ylabel 'Radius (m)'

R_c_eq = )" << results.cavity_radius << R"(

plot ')" << base_name << R"(_history.dat' u 1:2 w l lw 3 lc 'black' title 'Cavity', \
     R_c_eq w l dt 2 lc 'gray' title 'Equilibrium'

#------------------------------------------------------------------------------
# Plot 7: Energy Partitioning (pie chart approximation as bar)
#------------------------------------------------------------------------------
set title 'Energy Partitioning'
set xlabel ''
set ylabel 'Energy Fraction (%)'
set style data histogram
set style histogram cluster gap 1
set style fill solid 0.7 border -1
set boxwidth 0.8
set xtics ('Kinetic' 0, 'Internal' 1, 'Plastic' 2, 'Seismic' 3, 'Fracture' 4)

$energy << EOD
)" << 10.0 << R"(
)" << 80.0 << R"(
)" << 9.0 << R"(
)" << 0.01 << R"(
)" << 0.99 << R"(
EOD

plot $energy u 1 w boxes lc rgb '#4477AA' notitle

#------------------------------------------------------------------------------
# Plot 8: Arrival Time vs Distance
#------------------------------------------------------------------------------
set title 'Shock Wave Arrival Time'
set xlabel 'Distance (m)'
set ylabel 'Arrival Time (s)'
set xtics auto
set logscale x

plot ')" << base_name << R"(_ground_motion.dat' u 1:5 w l lw 3 lc 'purple' title 'Shock Arrival'

unset logscale

#------------------------------------------------------------------------------
# Plot 9: Summary Information
#------------------------------------------------------------------------------
unset border
unset tics
unset xlabel
unset ylabel
set title 'Summary'

set label 1 sprintf("Yield: %.1f kt TNT", )" << source.yield_kt << R"() at screen 0.72, 0.30 left
set label 2 sprintf("Depth of Burial: %.0f m", )" << source.depth << R"() at screen 0.72, 0.27 left
set label 3 sprintf("Cavity Radius: %.1f m", )" << results.cavity_radius << R"() at screen 0.72, 0.24 left
set label 4 sprintf("Crushed Zone: %.1f m", )" << results.crushed_zone_radius << R"() at screen 0.72, 0.21 left
set label 5 sprintf("Fractured Zone: %.1f m", )" << results.fractured_zone_radius << R"() at screen 0.72, 0.18 left
set label 6 sprintf("Damaged Zone: %.1f m", )" << results.damaged_zone_radius << R"() at screen 0.72, 0.15 left
set label 7 sprintf("Moment Magnitude: M_w %.2f", )" << results.moment_magnitude << R"() at screen 0.72, 0.12 left
set label 8 sprintf("Body Wave Magnitude: m_b %.2f", )" << results.body_wave_magnitude << R"() at screen 0.72, 0.09 left
set label 9 sprintf("Scalar Moment: %.2e N-m", )" << results.scalar_moment << R"() at screen 0.72, 0.06 left
set label 10 sprintf("Corner Frequency: %.2f Hz", )" << results.corner_frequency << R"() at screen 0.72, 0.03 left

plot 1/0 notitle

# Clean up
unset multiplot
unset label 1
unset label 2
unset label 3
unset label 4
unset label 5
unset label 6
unset label 7
unset label 8
unset label 9
unset label 10

print sprintf("Visualization saved to: )" << base_name << R"(_complete_results.png")
)";
}

// =============================================================================
// Main Simulation
// =============================================================================

int main(int argc, char* argv[]) {
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Print header
    printSeparator('=');
    std::cout << "  CRAM3D-Style Underground Explosion Near-Field Simulation\n";
    std::cout << "  =========================================================\n";
    std::cout << "  Near-field damage modeling with:\n";
    std::cout << "    - Cavity formation and expansion\n";
    std::cout << "    - Crushed, fractured, and damaged zones\n";
    std::cout << "    - Shock wave attenuation\n";
    std::cout << "    - Material damage evolution\n";
    std::cout << "    - Seismic source generation (RDP model)\n";
    printSeparator('=');
    std::cout << std::endl;
    
    // Parse command line arguments
    std::string config_file = "";
    double yield_kt = 100.0;  // Default 100 kt
    
    if (argc > 1 && std::string(argv[1]) != "") {
        config_file = argv[1];
    }
    if (argc > 2) {
        try {
            yield_kt = std::stod(argv[2]);
        } catch (...) {
            std::cerr << "Warning: Invalid yield value, using default 100 kt\n";
        }
    }
    
    // ==========================================================================
    // Setup
    // ==========================================================================
    std::cout << "Setting up simulation...\n\n";
    
    // Create source parameters
    UndergroundExplosionSource source;
    source.yield_kt = yield_kt;
    source.depth = 1000.0;  // 1 km depth
    source.location = {0.0, 0.0, -source.depth};
    
    // Host rock properties (granite)
    source.host_density = 2650.0;     // kg/m³
    source.host_vp = 5500.0;          // m/s
    source.host_vs = 3200.0;          // m/s
    source.host_porosity = 0.01;
    source.overshoot = 1.2;
    source.computeRiseTime();
    
    // Create EOS
    MieGruneisenEOS eos;
    eos.rho0 = source.host_density;
    eos.c0 = 4500.0;                  // Bulk sound speed
    eos.s = 1.4;                      // Hugoniot slope
    eos.gamma0 = 1.5;                 // Grüneisen parameter
    eos.a = 0.5;
    
    // Create strength model
    PressureDependentStrength strength;
    strength.Y0 = 50.0e6;             // 50 MPa cohesion
    strength.A = 0.5;                 // Pressure coefficient
    strength.n = 0.5;                 // Pressure exponent
    strength.B = 5.0;                 // Tension coefficient
    strength.Y_max = 2.0e9;           // 2 GPa maximum
    strength.Y_min = 0.1e6;           // 0.1 MPa residual
    strength.softening_strain = 0.05;
    
    // Create damage model
    DamageEvolutionModel damage;
    damage.strain_threshold = 0.001;
    damage.strain_saturation = 0.10;
    damage.pressure_threshold = 1.0e9;
    damage.pressure_saturation = 10.0e9;
    damage.tensile_threshold = 5.0e6;
    damage.damage_rate = 10.0;
    damage.max_damage = 0.99;
    
    // Create solver
    NearFieldExplosionSolver solver;
    solver.setSource(source);
    solver.setEOS(eos);
    solver.setStrengthModel(strength);
    solver.setDamageModel(damage);
    
    // If config file provided, override with file settings
    if (!config_file.empty()) {
        NearFieldExplosionConfig config;
        if (config.parse(config_file)) {
            std::cout << "Loaded configuration from: " << config_file << "\n\n";
            solver = *config.createSolver();
        } else {
            std::cerr << "Warning: Could not parse config file, using defaults\n\n";
        }
    }
    
    solver.initialize();
    
    // Create seismic source model
    RDPSeismicSource rdp;
    rdp.initialize(source);
    
    // ==========================================================================
    // Print Configuration
    // ==========================================================================
    std::cout << "Explosion Parameters:\n";
    printSeparator('-', 40);
    std::cout << "  Yield:             " << source.yield_kt << " kt TNT\n";
    std::cout << "  Energy:            " << formatWithPrefix(source.energyJoules(), "J") << "\n";
    std::cout << "  Depth of Burial:   " << source.depth << " m\n";
    std::cout << "  Rise Time:         " << source.rise_time * 1000.0 << " ms\n";
    std::cout << std::endl;
    
    std::cout << "Host Rock Properties:\n";
    printSeparator('-', 40);
    std::cout << "  Density:           " << source.host_density << " kg/m³\n";
    std::cout << "  P-wave velocity:   " << source.host_vp << " m/s\n";
    std::cout << "  S-wave velocity:   " << source.host_vs << " m/s\n";
    std::cout << "  Porosity:          " << source.host_porosity * 100 << " %\n";
    std::cout << std::endl;
    
    std::cout << "Computed Zone Radii:\n";
    printSeparator('-', 40);
    std::cout << "  Cavity radius:     " << std::fixed << std::setprecision(1) 
              << source.cavityRadius() << " m\n";
    std::cout << "  Crushed zone:      " << source.crushedZoneRadius() << " m\n";
    std::cout << "  Fractured zone:    " << source.fracturedZoneRadius() << " m\n";
    std::cout << "  Damaged zone:      " << source.damagedZoneRadius() << " m\n";
    std::cout << std::endl;
    
    std::cout << "Seismic Source:\n";
    printSeparator('-', 40);
    std::cout << "  Scalar moment:     " << std::scientific << std::setprecision(2) 
              << rdp.scalarMoment() << " N·m\n";
    std::cout << "  Moment magnitude:  " << std::fixed << std::setprecision(2) 
              << rdp.momentMagnitude() << "\n";
    std::cout << "  Body wave mag:     " << source.bodyWaveMagnitude() << "\n";
    std::cout << "  Corner frequency:  " << rdp.corner_frequency << " Hz\n";
    std::cout << std::endl;
    
    // ==========================================================================
    // Time Integration
    // ==========================================================================
    std::cout << "Running time integration...\n";
    printSeparator('-', 40);
    
    double dt = 0.001;  // 1 ms time step
    double t_max = 5.0;  // 5 seconds
    int n_steps = static_cast<int>(t_max / dt);
    int output_interval = 100;
    
    // Storage for time history
    std::vector<double> times;
    std::vector<double> cavity_radii;
    std::vector<double> cavity_pressures;
    std::vector<std::array<double, 6>> moment_tensors;
    std::vector<double> moment_rates;
    
    // Time stepping loop
    double time = 0.0;
    for (int step = 0; step <= n_steps; ++step) {
        solver.prestep(time, dt);
        solver.step(dt);
        solver.poststep(time, dt);
        
        // Store time history
        if (step % output_interval == 0) {
            times.push_back(time);
            cavity_radii.push_back(solver.getCavityRadius(time));
            
            // Cavity pressure (decaying)
            double R_c = source.cavityRadius();
            double R_current = solver.getCavityRadius(time);
            double P_cavity = source.initialCavityPressure() * 
                             std::pow(R_c / R_current, 3);
            cavity_pressures.push_back(P_cavity);
            
            // Moment tensor
            std::array<double, 6> M;
            double M_dot;
            solver.computeSourceFunction(time, M, M_dot);
            moment_tensors.push_back(M);
            moment_rates.push_back(M_dot);
        }
        
        // Progress output
        if (step % (n_steps / 10) == 0) {
            std::cout << "  t = " << std::fixed << std::setprecision(2) 
                      << time << " s (" << (100 * step / n_steps) << "%)\n";
        }
        
        time += dt;
    }
    
    std::cout << "  Done!\n\n";
    
    // ==========================================================================
    // Analysis
    // ==========================================================================
    std::cout << "Performing analysis...\n";
    printSeparator('-', 40);
    
    NearFieldAnalysisResults results = analyzeNearField(solver);
    
    // Print results
    std::cout << "\nDamage Zone Volumes:\n";
    std::cout << "  Cavity volume:     " << formatWithPrefix(results.cavity_volume, "m³") << "\n";
    std::cout << "  Crushed volume:    " << formatWithPrefix(results.crushed_volume, "m³") << "\n";
    std::cout << "  Fractured volume:  " << formatWithPrefix(results.fractured_volume, "m³") << "\n";
    std::cout << "  Damaged volume:    " << formatWithPrefix(results.damaged_volume, "m³") << "\n";
    std::cout << "  Total damage vol:  " << formatWithPrefix(results.total_damage_volume, "m³") << "\n";
    std::cout << std::endl;
    
    std::cout << "Energy Partitioning:\n";
    auto energy = solver.getEnergyPartition();
    std::cout << "  Total energy:      " << formatWithPrefix(energy.total, "J") << "\n";
    std::cout << "  Seismic energy:    " << formatWithPrefix(energy.seismic, "J") 
              << " (" << std::scientific << std::setprecision(4) 
              << 100.0 * energy.seismic / energy.total << " %)\n";
    std::cout << "  Seismic efficiency: " << 100.0 * results.seismic_efficiency << " %\n";
    std::cout << std::endl;
    
    // ==========================================================================
    // Output Files
    // ==========================================================================
    std::string base_name = "explosion_cram3d";
    
    std::cout << "Writing output files...\n";
    printSeparator('-', 40);
    
    // Time history
    writeTimeHistory(base_name + "_history.dat", times, cavity_radii, 
                    cavity_pressures, moment_tensors, moment_rates);
    std::cout << "  " << base_name << "_history.dat\n";
    
    // Ground motion
    double max_dist = 5.0 * source.damagedZoneRadius();
    writeGroundMotion(base_name + "_ground_motion.dat", solver.getShockModel(),
                     source, max_dist);
    std::cout << "  " << base_name << "_ground_motion.dat\n";
    
    // Damage zone map
    double map_extent = 1.5 * source.damagedZoneRadius();
    writeDamageZoneMap(base_name + "_damage_map.dat", solver, map_extent);
    std::cout << "  " << base_name << "_damage_map.dat\n";
    
    // Source time function
    double source_duration = 10.0 / rdp.corner_frequency;  // 10 x corner period
    writeSourceTimeFunction(base_name + "_source.dat", rdp, source_duration);
    std::cout << "  " << base_name << "_source.dat\n";
    
    // Radial profiles (for plotting)
    solver.writeOutput(base_name + "_radial.dat", t_max);
    std::cout << "  " << base_name << "_radial.dat\n";
    
    // Zone definitions
    solver.writeZoneFile(base_name + "_zones.dat");
    std::cout << "  " << base_name << "_zones.dat\n";
    
    // Visualization script
    writeVisualizationScript(base_name, results, source);
    std::cout << "  " << base_name << "_visualize.gp\n";
    
    // Near-field analysis summary
    writeNearFieldPlotScript(base_name + "_analysis", results);
    std::cout << "  " << base_name << "_analysis_plot.gp\n";
    
    std::cout << std::endl;
    
    // ==========================================================================
    // Summary
    // ==========================================================================
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        end_time - start_time).count();
    
    printSeparator('=');
    std::cout << "  Simulation Complete!\n";
    printSeparator('=');
    std::cout << "\n";
    
    std::cout << "Results Summary:\n";
    printSeparator('-', 60);
    std::cout << std::left << std::setw(30) << "  Parameter" 
              << std::right << std::setw(25) << "Value" << "\n";
    printSeparator('-', 60);
    
    std::cout << std::left << std::setw(30) << "  Yield" 
              << std::right << std::setw(20) << source.yield_kt << " kt\n";
    std::cout << std::left << std::setw(30) << "  Depth of Burial" 
              << std::right << std::setw(20) << source.depth << " m\n";
    std::cout << std::left << std::setw(30) << "  Cavity Radius" 
              << std::right << std::setw(20) << std::fixed << std::setprecision(1) 
              << results.cavity_radius << " m\n";
    std::cout << std::left << std::setw(30) << "  Crushed Zone Radius" 
              << std::right << std::setw(20) << results.crushed_zone_radius << " m\n";
    std::cout << std::left << std::setw(30) << "  Fractured Zone Radius" 
              << std::right << std::setw(20) << results.fractured_zone_radius << " m\n";
    std::cout << std::left << std::setw(30) << "  Damaged Zone Radius" 
              << std::right << std::setw(20) << results.damaged_zone_radius << " m\n";
    std::cout << std::left << std::setw(30) << "  Moment Magnitude" 
              << std::right << std::setw(20) << std::setprecision(2) 
              << results.moment_magnitude << "\n";
    std::cout << std::left << std::setw(30) << "  Body Wave Magnitude" 
              << std::right << std::setw(20) << results.body_wave_magnitude << "\n";
    std::cout << std::left << std::setw(30) << "  Corner Frequency" 
              << std::right << std::setw(20) << results.corner_frequency << " Hz\n";
    
    printSeparator('-', 60);
    std::cout << "\n";
    
    std::cout << "Execution time: " << duration / 1000.0 << " s\n\n";
    
    std::cout << "To visualize results, run:\n";
    std::cout << "  gnuplot " << base_name << "_visualize.gp\n";
    std::cout << "  gnuplot " << base_name << "_analysis_plot.gp\n\n";
    
    std::cout << "Output files are in the current directory.\n\n";
    
    return 0;
}
