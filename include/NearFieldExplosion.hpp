/**
 * @file NearFieldExplosion.hpp
 * @brief Near-field explosion physics for underground detonations
 * 
 * Implements comprehensive near-field effects for underground explosions including
 * cavity formation, damage zones, shock wave propagation, and seismic source generation.
 * 
 * Physical basis and references:
 * - CRAM3D methodology: https://apps.dtic.mil/sti/tr/pdf/AD1064154.pdf
 * - Rimer & Cherry (1982) - Nuclear test phenomenology
 * - Terhune et al. (1977) - Near-field ground motion from underground explosions
 * - Mueller & Murphy (1971) - Seismic characteristics of underground nuclear detonations
 * - Rodean (1971) - Nuclear explosion seismology
 * - Glasstone & Dolan (1977) - Effects of nuclear weapons
 * - NTS empirical scaling relations
 * 
 * Key Physics:
 * 1. Shock wave propagation with Mie-Grüneisen EOS
 * 2. Material strength degradation under pressure
 * 3. Cavity formation (elastic rebound + vaporization)
 * 4. Damage zones (crushed, fractured, disturbed)
 * 5. Spall and tensile failure
 * 6. Seismic source generation (RDP model)
 * 7. Energy partitioning
 * 
 * @author FSRM Development Team
 */

#ifndef NEAR_FIELD_EXPLOSION_HPP
#define NEAR_FIELD_EXPLOSION_HPP

#include <vector>
#include <array>
#include <string>
#include <memory>
#include <functional>
#include <cmath>
#include <map>

namespace FSRM {

// =============================================================================
// Constants for Underground Explosion Physics
// =============================================================================
namespace ExplosionPhysics {
    // TNT equivalence
    constexpr double JOULES_PER_KT = 4.184e12;        // J/kiloton
    constexpr double JOULES_PER_MT = 4.184e15;        // J/megaton
    
    // Physical constants
    constexpr double BOLTZMANN = 1.381e-23;           // J/K
    constexpr double SPECIFIC_GAS_CONSTANT = 287.0;   // J/(kg·K) for air
    
    // Empirical scaling parameters (NTS data)
    constexpr double CAVITY_SCALING_COEFF = 55.0;     // meters/(kt^0.3)
    constexpr double CAVITY_SCALING_EXP = 0.295;      // cube-root scaling exponent
    
    // Damage zone ratios (relative to cavity radius)
    constexpr double CRUSHED_ZONE_RATIO = 2.5;        // ~2-3x cavity radius
    constexpr double FRACTURED_ZONE_RATIO = 5.0;      // ~4-6x cavity radius
    constexpr double DAMAGED_ZONE_RATIO = 10.0;       // ~8-12x cavity radius
    
    // Pressure thresholds
    constexpr double VAPORIZATION_PRESSURE = 100.0e9; // ~100 GPa
    constexpr double MELTING_PRESSURE = 30.0e9;       // ~30 GPa
    constexpr double CRUSHING_PRESSURE = 5.0e9;       // ~5 GPa (high porosity collapse)
    constexpr double FRACTURE_PRESSURE = 0.5e9;       // ~500 MPa
    constexpr double DAMAGE_PRESSURE = 0.1e9;         // ~100 MPa
}

// =============================================================================
// Enumerations
// =============================================================================

/**
 * @brief Equation of state type for material under shock
 */
enum class EOSType {
    LINEAR_US_UP,         // Linear Us-Up Hugoniot
    MIE_GRUNEISEN,        // Mie-Grüneisen EOS
    TILLOTSON,            // Tillotson (for melting/vaporization)
    JWL,                  // Jones-Wilkins-Lee (for explosives)
    SESAME                // Tabular EOS
};

/**
 * @brief Material damage state
 */
enum class DamageState {
    INTACT,               // Undamaged
    MICRO_FRACTURED,      // Microfractures (damage < 0.2)
    FRACTURED,            // Macroscopic fractures (damage 0.2-0.6)
    COMMINUTED,           // Highly fractured (damage 0.6-0.9)
    CRUSHED,              // Completely crushed (damage > 0.9)
    MELTED,               // Melted material
    VAPORIZED             // Vaporized (cavity interior)
};

/**
 * @brief Failure mode
 */
enum class FailureMode {
    NONE,                 // No failure
    COMPRESSIVE,          // Crushing under compression
    SHEAR,                // Shear failure
    TENSILE,              // Tensile/spall failure
    MIXED                 // Mixed-mode failure
};

// =============================================================================
// Material Models
// =============================================================================

/**
 * @brief Mie-Grüneisen Equation of State parameters
 * 
 * P = P_H(V) + γ/V * (E - E_H)
 * 
 * where P_H, E_H are Hugoniot reference curve values
 * and γ is the Grüneisen parameter
 */
struct MieGruneisenEOS {
    double rho0 = 2650.0;         // Reference density (kg/m³)
    double c0 = 4500.0;           // Bulk sound speed (m/s)
    double s = 1.4;               // Linear Us-Up slope
    double gamma0 = 1.5;          // Grüneisen parameter at reference
    double a = 0.5;               // First-order volume correction
    
    // Optional higher-order terms
    double s2 = 0.0;              // Quadratic Us-Up term
    double s3 = 0.0;              // Cubic Us-Up term
    
    // Compute Hugoniot pressure from compression μ = ρ/ρ₀ - 1
    double hugoniotPressure(double mu) const {
        if (mu >= 0) {
            // Compression
            double eta = mu / (1.0 + mu);
            double denom = 1.0 - s * eta;
            if (s2 != 0 || s3 != 0) {
                denom = 1.0 - (s * eta + s2 * eta * eta + s3 * eta * eta * eta);
            }
            if (denom < 0.01) denom = 0.01;
            return rho0 * c0 * c0 * mu / (denom * denom);
        } else {
            // Tension (expanded state)
            return rho0 * c0 * c0 * mu;
        }
    }
    
    // Compute pressure with internal energy
    double pressure(double rho, double E) const {
        double mu = rho / rho0 - 1.0;
        double P_H = hugoniotPressure(mu);
        double gamma = gamma0 * (1.0 - a * mu);
        double V = 1.0 / rho;
        double V0 = 1.0 / rho0;
        double E_H = P_H * (V0 - V) / 2.0;  // Hugoniot energy
        return P_H + gamma / V * (E - E_H);
    }
    
    // Sound speed
    double soundSpeed(double rho, double E) const {
        double mu = rho / rho0 - 1.0;
        double gamma = gamma0 * (1.0 - a * mu);
        double P = pressure(rho, E);
        // Isentropic bulk modulus approximation
        double K = rho * c0 * c0 * std::pow(rho / rho0, gamma);
        return std::sqrt(K / rho + gamma * P / rho);
    }
};

/**
 * @brief Pressure-dependent strength model (CRAM3D style)
 * 
 * Y = Y₀ + A*P^n for P ≥ 0
 * Y = max(Y_min, Y₀ + B*P) for P < 0
 * 
 * With damage: Y_eff = Y * (1 - D)
 */
struct PressureDependentStrength {
    double Y0 = 50.0e6;           // Cohesive strength at P=0 (Pa)
    double A = 0.5;               // Pressure coefficient
    double n = 0.5;               // Pressure exponent
    double B = 5.0;               // Tension coefficient
    double Y_max = 2.0e9;         // Maximum yield strength (Pa)
    double Y_min = 0.1e6;         // Minimum (residual) strength (Pa)
    
    // Softening parameters
    double softening_strain = 0.05;    // Strain for full softening
    double residual_strength_ratio = 0.1;
    
    // Compute yield strength
    double yieldStrength(double P, double damage = 0.0, 
                         double plastic_strain = 0.0) const {
        double Y;
        if (P >= 0) {
            Y = Y0 + A * std::pow(P, n);
        } else {
            Y = Y0 + B * P;
        }
        
        Y = std::max(Y_min, std::min(Y, Y_max));
        
        // Damage softening
        Y *= (1.0 - damage);
        
        // Strain softening
        if (plastic_strain > 0 && softening_strain > 0) {
            double softening = 1.0 - (1.0 - residual_strength_ratio) *
                              std::min(1.0, plastic_strain / softening_strain);
            Y *= softening;
        }
        
        return std::max(Y_min, Y);
    }
    
    // Tensile strength (cutoff)
    double tensileStrength(double damage = 0.0) const {
        double T0 = Y0 / 10.0;  // Typical rock: tensile ~ 10% of cohesion
        return T0 * (1.0 - damage);
    }
};

/**
 * @brief Damage evolution model
 * 
 * Based on strain accumulation and pressure history
 */
struct DamageEvolutionModel {
    // Damage thresholds
    double strain_threshold = 0.001;   // Strain for damage initiation
    double strain_saturation = 0.10;   // Strain for full damage
    
    // Pressure-based damage
    double pressure_threshold = 1.0e9;  // Pressure for pressure damage
    double pressure_saturation = 10.0e9;
    
    // Tensile damage
    double tensile_threshold = 5.0e6;   // Tensile stress for damage
    
    // Damage rate parameters
    double damage_rate = 10.0;          // Rate coefficient
    double max_damage = 0.99;           // Maximum damage
    
    // Healing (for long-term recovery)
    bool healing_enabled = false;
    double healing_rate = 0.01;         // 1/s
    
    // Compute damage increment
    double computeDamageIncrement(double current_damage, double strain,
                                  double pressure, double dt) const {
        if (current_damage >= max_damage) return 0.0;
        
        double d_dot = 0.0;
        
        // Strain-based damage
        if (strain > strain_threshold) {
            double strain_factor = (strain - strain_threshold) / 
                                   (strain_saturation - strain_threshold);
            strain_factor = std::min(1.0, strain_factor);
            d_dot += damage_rate * strain_factor;
        }
        
        // Pressure-based damage (crushing)
        if (pressure > pressure_threshold) {
            double press_factor = (pressure - pressure_threshold) /
                                  (pressure_saturation - pressure_threshold);
            press_factor = std::min(1.0, press_factor);
            d_dot += damage_rate * press_factor;
        }
        
        // Healing (optional)
        if (healing_enabled && d_dot == 0.0 && current_damage > 0.0) {
            d_dot = -healing_rate * current_damage;
        }
        
        double dD = d_dot * dt;
        double new_damage = current_damage + dD;
        
        return std::max(0.0, std::min(max_damage, new_damage)) - current_damage;
    }
    
    // Get damage state enum
    DamageState getDamageState(double damage, double temperature = 300.0) const {
        if (temperature > 2000.0) return DamageState::VAPORIZED;
        if (temperature > 1500.0) return DamageState::MELTED;
        if (damage > 0.9) return DamageState::CRUSHED;
        if (damage > 0.6) return DamageState::COMMINUTED;
        if (damage > 0.2) return DamageState::FRACTURED;
        if (damage > 0.01) return DamageState::MICRO_FRACTURED;
        return DamageState::INTACT;
    }
};

// =============================================================================
// Near-Field State Variables
// =============================================================================

/**
 * @brief Complete material state for near-field explosion physics
 */
struct NearFieldMaterialState {
    // Thermodynamic state
    double density = 2650.0;          // kg/m³
    double pressure = 0.0;            // Pa
    double internal_energy = 0.0;     // J/kg
    double temperature = 300.0;       // K
    double entropy = 0.0;             // J/(kg·K)
    
    // Stress tensor (Voigt: xx, yy, zz, xy, xz, yz)
    std::array<double, 6> stress = {0, 0, 0, 0, 0, 0};
    
    // Strain tensor
    std::array<double, 6> strain = {0, 0, 0, 0, 0, 0};
    std::array<double, 6> strain_rate = {0, 0, 0, 0, 0, 0};
    std::array<double, 6> plastic_strain = {0, 0, 0, 0, 0, 0};
    double equivalent_plastic_strain = 0.0;
    
    // Velocity
    std::array<double, 3> velocity = {0, 0, 0};
    
    // Damage and failure
    double damage = 0.0;              // Scalar damage [0, 1]
    DamageState damage_state = DamageState::INTACT;
    FailureMode failure_mode = FailureMode::NONE;
    double peak_pressure = 0.0;       // Maximum pressure experienced
    double peak_temperature = 0.0;    // Maximum temperature
    
    // Zone classification
    bool in_cavity = false;
    bool in_crushed_zone = false;
    bool in_fractured_zone = false;
    bool is_spalled = false;
    
    // Helper functions
    double meanStress() const {
        return (stress[0] + stress[1] + stress[2]) / 3.0;
    }
    
    double vonMisesStress() const {
        double s1 = stress[0] - pressure;
        double s2 = stress[1] - pressure;
        double s3 = stress[2] - pressure;
        double J2 = 0.5 * (s1*s1 + s2*s2 + s3*s3) +
                   stress[3]*stress[3] + stress[4]*stress[4] + stress[5]*stress[5];
        return std::sqrt(3.0 * J2);
    }
    
    double volumetricStrain() const {
        return strain[0] + strain[1] + strain[2];
    }
    
    double deviatoricStrainMagnitude() const {
        double e_mean = volumetricStrain() / 3.0;
        double e1 = strain[0] - e_mean;
        double e2 = strain[1] - e_mean;
        double e3 = strain[2] - e_mean;
        return std::sqrt(2.0/3.0 * (e1*e1 + e2*e2 + e3*e3 + 
                        2.0*(strain[3]*strain[3] + strain[4]*strain[4] + 
                             strain[5]*strain[5])));
    }
};

// =============================================================================
// Explosion Source Model
// =============================================================================

/**
 * @brief Underground explosion source parameters
 */
struct UndergroundExplosionSource {
    // Basic parameters
    double yield_kt = 100.0;          // Yield in kilotons TNT
    double depth = 1000.0;            // Depth of burial (m)
    std::array<double, 3> location = {0, 0, -1000};  // Source location (m)
    
    // Medium properties (for scaling)
    double host_density = 2650.0;     // kg/m³
    double host_vp = 5500.0;          // P-wave velocity (m/s)
    double host_vs = 3200.0;          // S-wave velocity (m/s)
    double host_porosity = 0.01;      // Fraction
    
    // Overburden
    double overburden_stress = 0.0;   // Lithostatic stress at depth (Pa)
    
    // Source time function parameters
    double rise_time = 0.001;         // Rise time (s) - yield dependent
    double overshoot = 1.2;           // Source overshoot factor
    
    // Compute energy in Joules
    double energyJoules() const {
        return yield_kt * ExplosionPhysics::JOULES_PER_KT;
    }
    
    // Compute cavity radius using NTS empirical scaling
    // R_c = 55 * W^0.295 * (ρ/2.65)^(-1/3.4) meters
    double cavityRadius() const {
        double rho_ratio = host_density / 2650.0;
        return ExplosionPhysics::CAVITY_SCALING_COEFF * 
               std::pow(yield_kt, ExplosionPhysics::CAVITY_SCALING_EXP) *
               std::pow(rho_ratio, -1.0/3.4);
    }
    
    // Compute zone radii
    double crushedZoneRadius() const {
        return cavityRadius() * ExplosionPhysics::CRUSHED_ZONE_RATIO;
    }
    
    double fracturedZoneRadius() const {
        return cavityRadius() * ExplosionPhysics::FRACTURED_ZONE_RATIO;
    }
    
    double damagedZoneRadius() const {
        return cavityRadius() * ExplosionPhysics::DAMAGED_ZONE_RATIO;
    }
    
    // Initial cavity pressure (for energy release)
    double initialCavityPressure() const {
        double R_c = cavityRadius();
        double V_c = 4.0/3.0 * M_PI * R_c * R_c * R_c;
        // P = γ * E / V (ideal gas approximation)
        // γ ≈ 1.25 for detonation products
        return 1.25 * energyJoules() / V_c;
    }
    
    // Compute scaled distance (Hopkinson scaling)
    double scaledDistance(double r) const {
        double W_kg = yield_kt * 1000.0 * 4.184e6 / 4.184e9;  // kg TNT equivalent
        return r / std::pow(W_kg, 1.0/3.0);
    }
    
    // Auto-compute rise time from yield
    void computeRiseTime() {
        // Empirical: τ ≈ 0.001 * W^0.4 seconds
        rise_time = 0.001 * std::pow(yield_kt, 0.4);
    }
    
    // Seismic moment (empirical relation)
    // log10(M0) ≈ 17.0 + 1.0 * log10(W_kt) for contained tests
    double scalarMoment() const {
        return std::pow(10.0, 17.0 + std::log10(yield_kt));
    }
    
    // Body wave magnitude mb
    // mb ≈ 4.0 + 0.75 * log10(W_kt) (approximate)
    double bodyWaveMagnitude() const {
        return 4.0 + 0.75 * std::log10(yield_kt);
    }
};

// =============================================================================
// Shock Wave Attenuation Model
// =============================================================================

/**
 * @brief Shock wave propagation and attenuation
 * 
 * Peak pressure decay: P(r) = P_0 * (r_0/r)^n
 * where n depends on material and regime
 */
class ShockAttenuationModel {
public:
    // Peak pressure at cavity wall
    double P0 = 1.0e11;               // Pa (~100 GPa)
    double r0 = 100.0;                // Reference radius (cavity) (m)
    
    // Decay exponents for different regimes
    double n_strong = 3.0;            // Strong shock (near cavity)
    double n_weak = 1.87;             // Weak shock (far field)
    double transition_radius = 500.0;  // Strong-weak transition (m)
    
    // Arrival time parameters
    double c0 = 5000.0;               // Ambient sound speed (m/s)
    
    // Set from explosion source
    void setFromSource(const UndergroundExplosionSource& src) {
        r0 = src.cavityRadius();
        P0 = src.initialCavityPressure();
        c0 = src.host_vp;
        
        // Transition at ~5x cavity radius
        transition_radius = 5.0 * r0;
    }
    
    // Peak pressure at radius r
    double peakPressure(double r) const {
        if (r < r0) return P0;
        
        double n;
        if (r < transition_radius) {
            // Strong shock regime
            n = n_strong;
        } else {
            // Weak shock regime with smooth transition
            double f = (r - transition_radius) / (2.0 * transition_radius);
            f = std::min(1.0, std::max(0.0, f));
            n = n_strong + (n_weak - n_strong) * f;
        }
        
        return P0 * std::pow(r0 / r, n);
    }
    
    // Peak particle velocity at radius r
    // u ≈ P / (ρ * c) in weak shock limit
    double peakParticleVelocity(double r, double rho, double c) const {
        double P = peakPressure(r);
        return P / (rho * c);
    }
    
    // Shock arrival time at radius r
    double arrivalTime(double r) const {
        if (r <= r0) return 0.0;
        
        // In strong shock, velocity > c0
        // Approximate with average velocity
        double P_avg = (P0 + peakPressure(r)) / 2.0;
        double rho = 2650.0;  // Approximate
        double v_shock = c0 + P_avg / (rho * c0) * 0.5;
        
        return (r - r0) / v_shock;
    }
    
    // Positive phase duration
    double positivePhaseDuration(double r) const {
        // Empirical: τ+ ≈ τ_0 * (r/r_0)^α
        double tau0 = 0.001 * std::pow(P0 / 1e9, 0.3);  // Base duration
        return tau0 * std::pow(r / r0, 0.8);
    }
};

// =============================================================================
// CRAM3D-Style Near-Field Solver
// =============================================================================

/**
 * @brief Main near-field explosion solver
 * 
 * Integrates all near-field physics:
 * - Shock propagation
 * - Material response (EOS + strength)
 * - Damage evolution
 * - Cavity formation
 * - Zone delineation
 */
class NearFieldExplosionSolver {
public:
    NearFieldExplosionSolver();
    ~NearFieldExplosionSolver() = default;
    
    // Configuration
    void setSource(const UndergroundExplosionSource& src);
    void setEOS(const MieGruneisenEOS& eos);
    void setStrengthModel(const PressureDependentStrength& strength);
    void setDamageModel(const DamageEvolutionModel& damage);
    void configure(const std::map<std::string, std::string>& config);
    
    // Initialize
    void initialize();
    
    // Time stepping
    void prestep(double time, double dt);
    void step(double dt);
    void poststep(double time, double dt);
    
    // Material point update (for coupling with FE solver)
    void updateMaterialPoint(NearFieldMaterialState& state,
                            double dt, double time);
    
    // Source function (for coupling with wave propagation)
    void computeSourceFunction(double time,
                              std::array<double, 6>& moment_tensor,
                              double& moment_rate);
    
    // Zone queries
    bool isInCavity(const std::array<double, 3>& point, double time) const;
    bool isInCrushedZone(const std::array<double, 3>& point) const;
    bool isInFracturedZone(const std::array<double, 3>& point) const;
    bool isSpalled(const std::array<double, 3>& point) const;
    
    // Field queries
    double getPeakPressure(const std::array<double, 3>& point) const;
    double getPeakTemperature(const std::array<double, 3>& point) const;
    double getDamage(const std::array<double, 3>& point) const;
    DamageState getDamageState(const std::array<double, 3>& point) const;
    
    // Zone radii (time-dependent for cavity)
    double getCavityRadius(double time) const;
    double getCrushedZoneRadius() const;
    double getFracturedZoneRadius() const;
    
    // Seismic source parameters
    double getScalarMoment() const;
    double getCornerFrequency() const;
    void getMomentTensor(double time, std::array<double, 6>& M) const;
    
    // Energy partitioning
    struct EnergyPartition {
        double total = 0.0;           // Total energy (J)
        double kinetic = 0.0;         // Kinetic energy
        double internal = 0.0;        // Internal (thermal)
        double plastic = 0.0;         // Plastic work (dissipated)
        double seismic = 0.0;         // Radiated seismic
        double fracture = 0.0;        // Fracture energy
        double cavity_pv = 0.0;       // P*dV work for cavity
    };
    EnergyPartition getEnergyPartition() const;
    
    // Output
    void writeOutput(const std::string& filename, double time);
    void writeZoneFile(const std::string& filename);
    
    // Accessors
    const UndergroundExplosionSource& getSource() const { return source; }
    const ShockAttenuationModel& getShockModel() const { return shock_model; }
    
private:
    // Models
    UndergroundExplosionSource source;
    MieGruneisenEOS eos;
    PressureDependentStrength strength;
    DamageEvolutionModel damage_model;
    ShockAttenuationModel shock_model;
    
    // State
    double current_time = 0.0;
    double cavity_radius_current = 0.0;
    bool initialized = false;
    
    // Cavity expansion model
    double cavity_velocity = 0.0;
    double cavity_acceleration = 0.0;
    
    // Energy tracking
    EnergyPartition energy;
    
    // Seismic source
    double scalar_moment = 0.0;
    double corner_frequency = 1.0;
    
    // Internal methods
    void computeCavityExpansion(double dt);
    void computeShockPosition(double dt);
    double distanceToSource(const std::array<double, 3>& point) const;
    
    // Reduced Displacement Potential (Mueller-Murphy)
    double computeRDP(double t) const;
    double computeRDPDerivative(double t) const;
};

// =============================================================================
// Spall Model
// =============================================================================

/**
 * @brief Spall (tensile failure) model
 * 
 * Spall occurs when tensile stress exceeds dynamic tensile strength.
 * Critical for surface ground motion in underground explosions.
 */
class SpallModel {
public:
    // Dynamic tensile strength (strain-rate dependent)
    double static_tensile_strength = 10.0e6;  // Pa
    double strain_rate_exponent = 0.1;        // Rate sensitivity
    double reference_strain_rate = 1.0;       // 1/s
    
    // Spall parameters
    double spall_velocity_threshold = 1.0;    // m/s (velocity for spall)
    double spall_energy_per_area = 1000.0;    // J/m² (fracture energy)
    
    // Compute dynamic tensile strength
    double dynamicTensileStrength(double strain_rate) const {
        double rate_factor = std::pow(std::abs(strain_rate) / reference_strain_rate,
                                     strain_rate_exponent);
        return static_tensile_strength * (1.0 + rate_factor);
    }
    
    // Check for spall
    bool checkSpall(double tensile_stress, double strain_rate) const {
        return tensile_stress > dynamicTensileStrength(strain_rate);
    }
    
    // Spall layer thickness (empirical)
    double spallThickness(double yield_kt, double depth) const {
        // Approximate spall depth ~ cavity radius * factor
        double R_c = ExplosionPhysics::CAVITY_SCALING_COEFF * 
                     std::pow(yield_kt, ExplosionPhysics::CAVITY_SCALING_EXP);
        double factor = 2.0 * std::exp(-depth / (10.0 * R_c));
        return R_c * factor;
    }
    
    // Maximum spall velocity
    double maxSpallVelocity(double peak_stress, double rho, double c) const {
        return peak_stress / (rho * c);
    }
};

// =============================================================================
// Seismic Source Generation
// =============================================================================

/**
 * @brief Reduced Displacement Potential (RDP) seismic source
 * 
 * Based on Mueller-Murphy (1971) and later refinements.
 * 
 * ψ(ω) = ψ_∞ / (1 + (ω/ωc)²)
 * 
 * where ψ_∞ = K * W^α for yield W and constants K, α
 */
class RDPSeismicSource {
public:
    // Source parameters (from explosion)
    double yield_kt = 100.0;
    std::array<double, 3> location = {0, 0, -1000};
    
    // Medium parameters
    double density = 2650.0;
    double p_velocity = 5500.0;
    double s_velocity = 3200.0;
    
    // RDP parameters (Mueller-Murphy)
    double psi_infinity = 0.0;        // Steady-state RDP (computed)
    double corner_frequency = 1.0;    // Hz
    double overshoot = 1.2;           // Overshoot factor
    
    // Source time function type
    enum class STFType {
        STEP,                 // Heaviside
        EXPONENTIAL,          // exp(-t/τ)
        HASKELL,              // Linear rise + constant
        BRUNE                 // Brune omega-squared
    };
    STFType stf_type = STFType::BRUNE;
    
    // Initialize from source parameters
    void initialize(const UndergroundExplosionSource& src) {
        yield_kt = src.yield_kt;
        location = src.location;
        density = src.host_density;
        p_velocity = src.host_vp;
        s_velocity = src.host_vs;
        
        // Compute steady-state RDP (Mueller-Murphy relation)
        // ψ_∞ ≈ 1.5e5 * W^0.8 m³ (for hard rock)
        psi_infinity = 1.5e5 * std::pow(yield_kt, 0.8);
        
        // Corner frequency (Patton relation)
        // fc ≈ 2.5 * W^(-1/3) Hz
        corner_frequency = 2.5 * std::pow(yield_kt, -1.0/3.0);
    }
    
    // Compute RDP at time t
    double psi(double t) const {
        double tau = 1.0 / (2.0 * M_PI * corner_frequency);
        
        switch (stf_type) {
            case STFType::STEP:
                return t > 0 ? psi_infinity : 0.0;
                
            case STFType::EXPONENTIAL:
                return t > 0 ? psi_infinity * (1.0 - std::exp(-t / tau)) : 0.0;
                
            case STFType::HASKELL: {
                double rise_time = 2.0 * tau;
                if (t <= 0) return 0.0;
                if (t < rise_time) return psi_infinity * t / rise_time;
                return psi_infinity;
            }
            
            case STFType::BRUNE:
            default: {
                // Brune source with overshoot
                if (t <= 0) return 0.0;
                double x = t / tau;
                return psi_infinity * (1.0 - (1.0 + x) * std::exp(-x)) * overshoot;
            }
        }
    }
    
    // Compute RDP derivative (proportional to moment rate)
    double psiDot(double t) const {
        double tau = 1.0 / (2.0 * M_PI * corner_frequency);
        
        switch (stf_type) {
            case STFType::STEP:
                return 0.0;  // Delta function
                
            case STFType::EXPONENTIAL:
                return t > 0 ? (psi_infinity / tau) * std::exp(-t / tau) : 0.0;
                
            case STFType::HASKELL: {
                double rise_time = 2.0 * tau;
                if (t <= 0 || t > rise_time) return 0.0;
                return psi_infinity / rise_time;
            }
            
            case STFType::BRUNE:
            default: {
                if (t <= 0) return 0.0;
                double x = t / tau;
                return (psi_infinity * overshoot / tau) * x * std::exp(-x);
            }
        }
    }
    
    // Scalar seismic moment
    double scalarMoment() const {
        double K = density * p_velocity * p_velocity;
        return 4.0 * M_PI * K * psi_infinity;
    }
    
    // Moment magnitude
    double momentMagnitude() const {
        double M0 = scalarMoment();
        return (std::log10(M0) - 9.1) / 1.5;
    }
    
    // Moment tensor (predominantly isotropic for explosions)
    void momentTensor(double t, std::array<double, 6>& M,
                     double iso_fraction = 0.7,
                     double clvd_fraction = 0.25) const {
        double psi_t = psi(t);
        double K = density * p_velocity * p_velocity;
        double M0 = 4.0 * M_PI * K * psi_t;
        
        // Isotropic component
        double M_iso = iso_fraction * M0 / 3.0;
        
        // CLVD component (vertical axis)
        double M_clvd = clvd_fraction * M0;
        
        // DC component (remaining)
        double M_dc = (1.0 - iso_fraction - clvd_fraction) * M0;
        
        // Mxx, Myy, Mzz, Mxy, Mxz, Myz
        M[0] = M_iso + M_clvd / 3.0;  // Mxx
        M[1] = M_iso + M_clvd / 3.0;  // Myy
        M[2] = M_iso - 2.0 * M_clvd / 3.0 + M_dc;  // Mzz
        M[3] = 0.0;  // Mxy
        M[4] = 0.0;  // Mxz
        M[5] = 0.0;  // Myz
    }
    
    // Moment rate tensor
    void momentRateTensor(double t, std::array<double, 6>& Mdot,
                         double iso_fraction = 0.7,
                         double clvd_fraction = 0.25) const {
        double psi_dot = psiDot(t);
        double K = density * p_velocity * p_velocity;
        double M0_dot = 4.0 * M_PI * K * psi_dot;
        
        double M_iso_dot = iso_fraction * M0_dot / 3.0;
        double M_clvd_dot = clvd_fraction * M0_dot;
        double M_dc_dot = (1.0 - iso_fraction - clvd_fraction) * M0_dot;
        
        Mdot[0] = M_iso_dot + M_clvd_dot / 3.0;
        Mdot[1] = M_iso_dot + M_clvd_dot / 3.0;
        Mdot[2] = M_iso_dot - 2.0 * M_clvd_dot / 3.0 + M_dc_dot;
        Mdot[3] = 0.0;
        Mdot[4] = 0.0;
        Mdot[5] = 0.0;
    }
};

// =============================================================================
// Analysis and Output
// =============================================================================

/**
 * @brief Near-field analysis results
 */
struct NearFieldAnalysisResults {
    // Zone dimensions
    double cavity_radius = 0.0;
    double crushed_zone_radius = 0.0;
    double fractured_zone_radius = 0.0;
    double damaged_zone_radius = 0.0;
    
    // Cavity properties
    double cavity_volume = 0.0;
    double cavity_pressure = 0.0;
    double cavity_temperature = 0.0;
    
    // Damage statistics
    double crushed_volume = 0.0;
    double fractured_volume = 0.0;
    double damaged_volume = 0.0;
    double total_damage_volume = 0.0;
    double average_damage = 0.0;
    
    // Energy
    double total_energy = 0.0;
    double seismic_energy = 0.0;
    double seismic_efficiency = 0.0;
    
    // Seismic source
    double scalar_moment = 0.0;
    double moment_magnitude = 0.0;
    double body_wave_magnitude = 0.0;
    double corner_frequency = 0.0;
    
    // Spall
    bool spall_occurred = false;
    double spall_velocity = 0.0;
    double spall_area = 0.0;
    
    // Peak values at various distances
    std::vector<double> radii;
    std::vector<double> peak_pressures;
    std::vector<double> peak_velocities;
    std::vector<double> peak_accelerations;
    std::vector<double> damage_values;
};

/**
 * @brief Perform comprehensive near-field analysis
 */
NearFieldAnalysisResults analyzeNearField(const NearFieldExplosionSolver& solver);

/**
 * @brief Write Gnuplot visualization script
 */
void writeNearFieldPlotScript(const std::string& base_name,
                              const NearFieldAnalysisResults& results);

// =============================================================================
// Configuration Parser
// =============================================================================

/**
 * @brief Parse near-field explosion configuration
 */
class NearFieldExplosionConfig {
public:
    NearFieldExplosionConfig() = default;
    
    bool parse(const std::string& filename);
    bool parse(const std::map<std::string, std::string>& config);
    
    // Create solver with configuration
    std::unique_ptr<NearFieldExplosionSolver> createSolver() const;
    
    // Access individual components
    UndergroundExplosionSource getSource() const { return source; }
    MieGruneisenEOS getEOS() const { return eos; }
    PressureDependentStrength getStrength() const { return strength; }
    DamageEvolutionModel getDamageModel() const { return damage; }
    
private:
    UndergroundExplosionSource source;
    MieGruneisenEOS eos;
    PressureDependentStrength strength;
    DamageEvolutionModel damage;
    
    bool configured = false;
};

} // namespace FSRM

#endif // NEAR_FIELD_EXPLOSION_HPP
