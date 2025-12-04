/**
 * @file ExplosionImpactPhysics.hpp
 * @brief Physics models for explosions, nuclear detonations, and hypervelocity impacts
 * 
 * Implements source models and coupled physics for:
 * - Underground nuclear tests (cavity, damage zones, seismic coupling)
 * - Atmospheric nuclear detonations (fireball, blast, radiation, EMP)
 * - Surface explosions and cratering
 * - Hypervelocity impact events (asteroids, meteorites)
 * - Radiation transport and electromagnetic pulse generation
 * - Atmospheric coupling and shock propagation
 */

#ifndef EXPLOSION_IMPACT_PHYSICS_HPP
#define EXPLOSION_IMPACT_PHYSICS_HPP

#include <vector>
#include <array>
#include <string>
#include <memory>
#include <functional>
#include <cmath>
#include <map>

namespace FSRM {

// Forward declarations
class PhysicsKernel;

// =============================================================================
// Constants
// =============================================================================
namespace ExplosionConstants {
    constexpr double JOULES_PER_KILOTON = 4.184e12;     // J/kt TNT
    constexpr double JOULES_PER_MEGATON = 4.184e15;     // J/Mt TNT
    constexpr double SPEED_OF_LIGHT = 2.998e8;          // m/s
    constexpr double ELECTRON_CHARGE = 1.602e-19;       // Coulombs
    constexpr double ELECTRON_MASS = 9.109e-31;         // kg
    constexpr double BOLTZMANN = 1.381e-23;             // J/K
    constexpr double STEFAN_BOLTZMANN = 5.670e-8;       // W/m²·K⁴
    constexpr double EARTH_MAGNETIC_FIELD = 5.0e-5;     // Tesla
    constexpr double ATMOSPHERE_SCALE_HEIGHT = 8500.0;  // m
}

// =============================================================================
// Enumerations
// =============================================================================

/**
 * @brief Types of explosive sources
 */
enum class ExplosionType {
    CHEMICAL_HE,              // High explosive (conventional)
    NUCLEAR_FISSION,          // Pure fission weapon
    NUCLEAR_THERMONUCLEAR,    // Fusion-boosted or staged thermonuclear
    NUCLEAR_ENHANCED_RAD,     // Enhanced radiation (neutron bomb)
    VOLCANIC,                 // Volcanic explosion
    INDUSTRIAL                // Industrial explosion (e.g., fertilizer)
};

/**
 * @brief Types of impact events
 */
enum class ImpactorType {
    STONY_ASTEROID,           // Ordinary chondrite
    CARBONACEOUS_ASTEROID,    // C-type asteroid
    IRON_ASTEROID,            // Iron-nickel meteorite
    COMETARY,                 // Icy comet nucleus
    ICY_BODY,                 // Outer solar system
    ARTIFICIAL                // Kinetic impactor (DART-like)
};

/**
 * @brief Burst configurations for nuclear weapons
 */
enum class BurstType {
    UNDERGROUND_CONTAINED,    // Fully contained, no surface effects
    UNDERGROUND_VENTING,      // Partial containment, some venting
    SURFACE,                  // Contact burst
    LOW_ALTITUDE,             // Below fireball max radius from ground
    OPTIMUM_ALTITUDE,         // Max blast effect at surface
    HIGH_ALTITUDE,            // Above significant blast
    EXOATMOSPHERIC            // Space detonation (EMP optimized)
};

/**
 * @brief EMP component types
 */
enum class EMPComponent {
    E1_FAST,                  // Gamma-induced Compton current (ns rise)
    E2_INTERMEDIATE,          // Scattered gamma (μs-ms)
    E3_SLOW,                  // MHD-EMP from fireball (seconds)
    SGEMP                     // System-generated EMP (internal)
};

/**
 * @brief Crater formation stages
 */
enum class CraterStage {
    CONTACT,                  // Impactor contacts target
    COMPRESSION,              // Shock propagation through impactor
    EXCAVATION,               // Material ejection phase
    MODIFICATION,             // Collapse and final form
    RELAXATION                // Long-term isostatic adjustment
};

// =============================================================================
// Explosion Source Models
// =============================================================================

/**
 * @brief Parameters for nuclear explosion source
 */
struct NuclearSourceParameters {
    double yield_kt = 100.0;              // Yield in kilotons TNT
    double depth_of_burial = 0.0;         // meters (positive = underground)
    double height_of_burst = 0.0;         // meters AGL (for air bursts)
    double fission_fraction = 0.5;        // Fraction of yield from fission
    
    // Location
    double x = 0.0, y = 0.0, z = 0.0;
    
    // Derived quantities
    double total_energy_joules() const { 
        return yield_kt * ExplosionConstants::JOULES_PER_KILOTON; 
    }
    
    // Scaling laws for cavity and damage zones
    double cavity_radius(double rock_density = 2650.0) const;
    double crushed_zone_radius() const;
    double fractured_zone_radius() const;
    
    // Seismic source parameters
    double scalar_moment() const;        // N·m
    double body_wave_magnitude() const;  // m_b
    double surface_wave_magnitude() const; // M_s
};

/**
 * @brief Mueller-Murphy source model for underground nuclear explosions
 * 
 * Reduced Displacement Potential (RDP) model for seismic source:
 * ψ(ω) = K * W^n / [1 + (ω/ω_c)²]
 */
class MuellerMurphySource {
public:
    MuellerMurphySource();
    
    void setParameters(const NuclearSourceParameters& params);
    void setMediumProperties(double rho, double vp, double vs);
    
    // Source time function (moment rate)
    double momentRate(double t) const;
    
    // Reduced displacement potential in frequency domain
    std::complex<double> rdp(double omega) const;
    
    // Corner frequency
    double getCornerFrequency() const;
    
    // Overshoot parameter
    double getOvershoot() const;
    
    // Moment tensor (isotropic + CLVD)
    void getMomentTensor(double t, double* M) const;
    
private:
    NuclearSourceParameters source_params;
    double density;
    double p_velocity;
    double s_velocity;
    
    // Derived
    double corner_frequency;
    double scalar_moment;
    double overshoot;
    double rise_time;
    
    void computeDerivedQuantities();
};

/**
 * @brief Spherical cavity source (equivalent to explosion)
 */
class SphericalCavitySource {
public:
    SphericalCavitySource();
    
    void setCavityRadius(double R);
    void setMediumProperties(double rho, double K, double G);
    void setOverpressure(double P0);
    void setRiseTime(double tau);
    
    // Radial displacement at distance r and time t
    double displacement(double r, double t) const;
    
    // Radial velocity at distance r and time t
    double velocity(double r, double t) const;
    
    // Stress components
    void stress(double r, double t, double& sigma_rr, double& sigma_tt) const;
    
    // Far-field equivalent source strength
    double equivalentMoment() const;
    
private:
    double cavity_radius;
    double density;
    double bulk_modulus;
    double shear_modulus;
    double initial_pressure;
    double rise_time;
    
    double p_velocity;
    double s_velocity;
};

// =============================================================================
// Atmospheric Explosion Models
// =============================================================================

/**
 * @brief Fireball model for atmospheric nuclear detonations
 * 
 * Models the evolution of the fireball including:
 * - Initial X-ray emission and absorption
 * - Hydrodynamic expansion
 * - Thermal radiation pulse
 * - Rise and stabilization
 */
class FireballModel {
public:
    enum class Phase {
        INITIAL_FLASH,      // X-ray heating
        MINIMUM,            // First minimum
        SECOND_MAXIMUM,     // Shock breakout
        COOLING             // Gradual cooling
    };
    
    FireballModel();
    
    void setYield(double yield_kt);
    void setBurstHeight(double height_m);
    void setAtmosphericDensity(double rho);
    
    // Fireball radius as function of time
    double radius(double t) const;
    
    // Fireball temperature
    double temperature(double t) const;
    
    // Thermal radiation power
    double thermalPower(double t) const;
    
    // Thermal fluence at distance r at time t
    double thermalFluence(double r, double t) const;
    
    // Maximum fireball radius
    double maxRadius() const;
    
    // Rise dynamics (mushroom cloud)
    double cloudTopHeight(double t) const;
    double cloudStabilizationHeight() const;
    double cloudStabilizationTime() const;
    
    // Current phase
    Phase getCurrentPhase(double t) const;
    
private:
    double yield_kt;
    double burst_height;
    double ambient_density;
    
    // Time markers
    double t_first_minimum;
    double t_second_maximum;
    double t_stabilization;
    
    void computeTimeMarkers();
};

/**
 * @brief Atmospheric blast wave model
 * 
 * Sedov-Taylor self-similar solution with modifications for:
 * - Non-ideal atmosphere
 * - Ground reflection (Mach stem)
 * - Thermal precursor
 */
class AtmosphericBlastWave {
public:
    AtmosphericBlastWave();
    
    void setYield(double yield_kt);
    void setBurstHeight(double height_m);
    void setAtmosphericProfile(std::function<double(double)> rho_z);
    void setGroundReflection(bool enable);
    
    // Shock front position at time t
    double shockRadius(double t) const;
    
    // Shock velocity at time t
    double shockVelocity(double t) const;
    
    // Peak overpressure at distance r
    double peakOverpressure(double r) const;
    
    // Dynamic pressure at distance r
    double dynamicPressure(double r) const;
    
    // Positive phase duration at distance r
    double positivePhaseDuration(double r) const;
    
    // Impulse at distance r
    double positiveImpulse(double r) const;
    double negativeImpulse(double r) const;
    
    // Arrival time at distance r
    double arrivalTime(double r) const;
    
    // Mach stem height at distance r (for surface bursts)
    double machStemHeight(double r) const;
    
    // Is point in Mach reflection region?
    bool inMachRegion(double r, double z) const;
    
    // Full pressure waveform at distance r
    void pressureWaveform(double r, double* times, double* pressures, int n_points) const;
    
private:
    double yield_kt;
    double burst_height;
    bool ground_reflection;
    std::function<double(double)> density_profile;
    
    // Glasstone-Dolan scaling parameters
    double sedovTaylorRadius(double t) const;
};

// =============================================================================
// Radiation Transport
// =============================================================================

/**
 * @brief Prompt radiation model (gamma and neutron)
 */
class PromptRadiationModel {
public:
    PromptRadiationModel();
    
    void setYield(double yield_kt);
    void setFissionFraction(double fraction);
    void setBurstHeight(double height_m);
    
    // Gamma dose at distance r (rad or Gy)
    double gammaDose(double r, double shielding_thickness = 0.0) const;
    
    // Neutron dose at distance r
    double neutronDose(double r, double shielding_thickness = 0.0) const;
    
    // Total prompt radiation dose
    double totalPromptDose(double r) const;
    
    // Dose rate vs time at fixed distance
    double doseRate(double r, double t) const;
    
    // Gamma spectrum (energy distribution)
    double gammaSpectrum(double energy_MeV) const;
    
    // Neutron spectrum  
    double neutronSpectrum(double energy_MeV) const;
    
private:
    double yield_kt;
    double fission_fraction;
    double burst_height;
    
    // Attenuation coefficients
    double gammaAttenuation(double r) const;
    double neutronAttenuation(double r) const;
};

/**
 * @brief Fallout model for nuclear debris
 */
class FalloutModel {
public:
    FalloutModel();
    
    void setYield(double yield_kt, double fission_yield_kt);
    void setDetonationParameters(double height_m, double x, double y);
    void setWindProfile(std::function<void(double z, double& vx, double& vy)> wind);
    
    // Particle size distribution (log-normal)
    double particleSizeDistribution(double diameter_um) const;
    
    // Cloud rise parameters
    double cloudTopHeight() const;
    double cloudBottomHeight() const;
    double cloudRadius() const;
    
    // Ground deposition at point (x, y) at time t (Bq/m²)
    double depositionRate(double x, double y, double t) const;
    
    // Total accumulated activity at point (x, y) after time t
    double accumulatedActivity(double x, double y, double t) const;
    
    // Dose rate at point (x, y) at time t (from ground shine)
    double groundDoseRate(double x, double y, double t) const;
    
    // Generate fallout contours (returns polygons for dose rate levels)
    void generateContours(double t, const std::vector<double>& levels,
                         std::vector<std::vector<std::pair<double, double>>>& contours) const;
    
private:
    double total_yield_kt;
    double fission_yield_kt;
    double burst_height;
    double burst_x, burst_y;
    std::function<void(double, double&, double&)> wind_profile;
    
    // Particle settling velocity
    double settlingVelocity(double diameter_um) const;
    
    // Activity vs time (Way-Wigner)
    double activityDecay(double t) const;
};

// =============================================================================
// Electromagnetic Pulse (EMP) Models
// =============================================================================

/**
 * @brief High-Altitude EMP (HEMP) model
 * 
 * Models all three EMP components:
 * - E1: Fast (gamma-induced Compton current)
 * - E2: Intermediate (scattered gamma, neutrons)
 * - E3: Slow (MHD-EMP from fireball distortion of Earth's field)
 */
class EMPModel {
public:
    EMPModel();
    
    void setYield(double yield_kt);
    void setBurstAltitude(double altitude_m);
    void setGeomagnetic(double field_T, double dip_angle_deg);
    
    // E1 component (fast)
    double E1_electricField(double r, double t) const;     // V/m
    double E1_peakField() const;
    double E1_riseTime() const;                            // seconds
    double E1_decayTime() const;
    
    // E2 component (intermediate)
    double E2_electricField(double r, double t) const;
    double E2_duration() const;
    
    // E3 component (slow, MHD-EMP)
    double E3_electricField(double r, double t) const;
    double E3_duration() const;
    
    // Total electric field at point and time
    void totalElectricField(double x, double y, double z, double t,
                           double& Ex, double& Ey, double& Ez) const;
    
    // Magnetic field
    void totalMagneticField(double x, double y, double z, double t,
                           double& Bx, double& By, double& Bz) const;
    
    // Induced current in a conductor of given length and orientation
    double inducedCurrent(double length, double orientation_angle, double t) const;
    
    // Coupled energy to power line
    double powerLineCoupling(double line_length, double line_height, double t) const;
    
    // Source region geometry
    void getSourceRegion(double& min_alt, double& max_alt, double& radius) const;
    
private:
    double yield_kt;
    double burst_altitude;
    double geomagnetic_field;
    double geomagnetic_dip;
    
    // Compton electron parameters
    double comptonElectronFlux(double altitude) const;
    double atmosphericConductivity(double altitude) const;
    
    // E1 source current
    double comptonCurrent(double altitude, double t) const;
};

/**
 * @brief Surface-burst EMP (SBEMP) model
 */
class SurfaceEMPModel {
public:
    SurfaceEMPModel();
    
    void setYield(double yield_kt);
    void setGroundConductivity(double sigma);  // S/m
    
    // Electric field from ground wave
    double groundWaveField(double r, double t) const;
    
    // Magnetic field disturbance
    double magneticDisturbance(double r, double t) const;
    
private:
    double yield_kt;
    double ground_conductivity;
};

// =============================================================================
// Impact Physics Models
// =============================================================================

/**
 * @brief Impactor parameters
 */
struct ImpactorParameters {
    ImpactorType type = ImpactorType::STONY_ASTEROID;
    double diameter = 1000.0;             // meters
    double density = 3000.0;              // kg/m³
    double velocity = 20000.0;            // m/s
    double impact_angle = 45.0;           // degrees from horizontal
    double azimuth = 0.0;                 // degrees
    double x = 0.0, y = 0.0, z = 0.0;     // Impact point
    
    // Derived quantities
    double mass() const { return (4.0/3.0) * M_PI * pow(diameter/2.0, 3) * density; }
    double kineticEnergy() const { return 0.5 * mass() * velocity * velocity; }
    double kineticEnergy_Mt() const { 
        return kineticEnergy() / ExplosionConstants::JOULES_PER_MEGATON; 
    }
};

/**
 * @brief Target material properties for impact
 */
struct ImpactTargetProperties {
    double density = 2700.0;              // kg/m³
    double porosity = 0.0;                // Volume fraction
    double strength = 100.0e6;            // Pa
    double sound_speed = 5000.0;          // m/s
    double gravity = 9.81;                // m/s² (planetary surface gravity)
    
    // Tillotson EOS parameters
    double tillotson_A = 18.0e9;
    double tillotson_B = 18.0e9;
    double tillotson_E0 = 16.0e6;
    double tillotson_alpha = 5.0;
    double tillotson_beta = 5.0;
};

/**
 * @brief Pi-group crater scaling model
 * 
 * Uses Holsapple-Schmidt scaling laws:
 * D = K * (ρ_i/ρ_t)^ν * (gR/v²)^-μ * R
 */
class CraterScalingModel {
public:
    enum class ScalingRegime {
        STRENGTH,             // Small craters, strength-dominated
        GRAVITY               // Large craters, gravity-dominated
    };
    
    CraterScalingModel();
    
    void setImpactor(const ImpactorParameters& imp);
    void setTarget(const ImpactTargetProperties& tgt);
    
    // Scaling regime
    ScalingRegime getScalingRegime() const;
    
    // Transient crater dimensions
    double transientCraterDiameter() const;
    double transientCraterDepth() const;
    double transientCraterVolume() const;
    
    // Final crater dimensions (after modification)
    double finalCraterDiameter() const;
    double finalCraterDepth() const;
    bool isComplexCrater() const;
    
    // Central peak (for complex craters)
    double centralPeakHeight() const;
    
    // Rim
    double rimHeight() const;
    double rimWidth() const;
    
    // Excavation depth
    double maxExcavationDepth() const;
    
    // Ejecta
    double ejectaMass() const;
    double ejectaVelocityAtRim() const;
    double ejectaBlanketThickness(double r) const;
    
    // Formation time
    double craterFormationTime() const;
    
private:
    ImpactorParameters impactor;
    ImpactTargetProperties target;
    
    // Pi-group scaling exponents
    double mu;           // Velocity exponent
    double nu;           // Density ratio exponent
    double K1, K2;       // Scaling constants
    
    void computeScalingParameters();
    double piV() const;  // Volume pi-group
    double pi2() const;  // Gravity pi-group
    double pi3() const;  // Strength pi-group
};

/**
 * @brief Z-model for crater excavation flow
 * 
 * Maxwell Z-model describes the flow field during excavation:
 * v_r = A/r² * (1 - z/z_c)^Z
 */
class ZModelExcavation {
public:
    ZModelExcavation();
    
    void setCraterParameters(double transient_radius, double transient_depth);
    void setZExponent(double Z);
    void setFlowCenter(double depth);
    
    // Velocity field during excavation
    void excavationVelocity(double r, double z, double t,
                           double& vr, double& vz) const;
    
    // Particle trajectory (Lagrangian)
    void particleTrajectory(double r0, double z0, 
                           std::vector<double>& r, std::vector<double>& z,
                           double t_max, double dt) const;
    
    // Ejecta launch position and velocity
    void ejectaLaunchConditions(double r_final, 
                               double& v_launch, double& angle_launch) const;
    
    // Landing distance for particle launched from rim
    double ejectaLandingDistance(double v_launch, double angle) const;
    
private:
    double crater_radius;
    double crater_depth;
    double Z_exponent;
    double flow_center_depth;
    
    double excavation_time;
};

/**
 * @brief Shock wave propagation from impact
 */
class ImpactShockWave {
public:
    ImpactShockWave();
    
    void setImpactor(const ImpactorParameters& imp);
    void setTarget(const ImpactTargetProperties& tgt);
    
    // Peak shock pressure at distance r
    double peakPressure(double r) const;
    
    // Shock velocity at distance r
    double shockVelocity(double r) const;
    
    // Particle velocity at distance r
    double particleVelocity(double r) const;
    
    // Hugoniot relations
    void hugoniotState(double pressure, 
                      double& density, double& energy, double& temperature) const;
    
    // Shock attenuation
    double pressureDecayExponent() const;
    
    // Isobaric core radius
    double isobaricCoreRadius() const;
    
    // Shock metamorphism boundaries
    double vaporizationRadius() const;    // Complete vaporization
    double meltRadius() const;            // Complete melting
    double incipientMeltRadius() const;   // Partial melting
    double shockedQuartzRadius() const;   // PDF formation
    double fractureRadius() const;        // Pervasive fracturing
    
private:
    ImpactorParameters impactor;
    ImpactTargetProperties target;
    
    double initial_pressure;  // Peak pressure at impact
    double isobaric_radius;
    double decay_exponent;
    
    void computeShockParameters();
};

/**
 * @brief Seismic source from impact
 */
class ImpactSeismicSource {
public:
    ImpactSeismicSource();
    
    void setImpactor(const ImpactorParameters& imp);
    void setTarget(const ImpactTargetProperties& tgt);
    void setSeismicEfficiency(double eta);  // Typically 10^-4 to 10^-5
    
    // Seismic moment
    double scalarMoment() const;
    
    // Equivalent earthquake magnitude
    double momentMagnitude() const;
    double bodyWaveMagnitude() const;
    double surfaceWaveMagnitude() const;
    
    // Source time function (moment rate)
    double momentRate(double t) const;
    
    // Equivalent point source parameters
    void equivalentMomentTensor(double* M) const;  // Predominantly isotropic
    
    // Corner frequency
    double cornerFrequency() const;
    
    // Source duration
    double sourceDuration() const;
    
private:
    ImpactorParameters impactor;
    ImpactTargetProperties target;
    double seismic_efficiency;
    
    double seismic_energy;
    double scalar_moment;
    double source_radius;
    double source_duration;
    
    void computeSourceParameters();
};

// =============================================================================
// Coupled Physics Manager
// =============================================================================

/**
 * @brief Manager for coupled explosion/impact physics
 * 
 * Coordinates multiple physics models for comprehensive simulation.
 */
class ExplosionImpactManager {
public:
    ExplosionImpactManager();
    
    // Configuration
    void configureNuclearTest(const NuclearSourceParameters& params, BurstType type);
    void configureImpact(const ImpactorParameters& imp, const ImpactTargetProperties& tgt);
    
    // Enable/disable physics
    void enableRadiation(bool enable);
    void enableEMP(bool enable);
    void enableAtmosphere(bool enable);
    void enableThermal(bool enable);
    void enableSeismic(bool enable);
    
    // Time stepping
    double suggestTimeStep(double current_time) const;
    
    // Source terms for different physics kernels
    void getMomentSource(double t, double* M) const;
    void getHeatSource(double x, double y, double z, double t, double& Q) const;
    void getRadiationSource(double x, double y, double z, double t,
                           double& gamma_flux, double& neutron_flux) const;
    void getEMPFields(double x, double y, double z, double t,
                     double& Ex, double& Ey, double& Ez,
                     double& Bx, double& By, double& Bz) const;
    
    // State queries
    double getCavityRadius(double t) const;
    double getFireballRadius(double t) const;
    double getShockRadius(double t) const;
    double getCraterRadius(double t) const;
    
    // Damage assessment
    double getDamageLevel(double x, double y, double z) const;
    double getPeakPressure(double r) const;
    double getPeakTemperature(double r) const;
    
private:
    // Source models
    std::unique_ptr<MuellerMurphySource> nuclear_source;
    std::unique_ptr<FireballModel> fireball;
    std::unique_ptr<AtmosphericBlastWave> blast_wave;
    std::unique_ptr<PromptRadiationModel> radiation;
    std::unique_ptr<FalloutModel> fallout;
    std::unique_ptr<EMPModel> emp;
    
    std::unique_ptr<CraterScalingModel> crater;
    std::unique_ptr<ImpactShockWave> impact_shock;
    std::unique_ptr<ImpactSeismicSource> impact_seismic;
    
    // Configuration flags
    bool is_nuclear;
    bool is_impact;
    bool radiation_enabled;
    bool emp_enabled;
    bool atmosphere_enabled;
    bool thermal_enabled;
    bool seismic_enabled;
};

// =============================================================================
// Configuration Parser
// =============================================================================

/**
 * @brief Parse explosion/impact configuration from file
 */
class ExplosionImpactConfig {
public:
    ExplosionImpactConfig();
    
    bool parse(const std::string& filename);
    bool parse(const std::map<std::string, std::string>& config);
    
    // Create simulation manager
    std::unique_ptr<ExplosionImpactManager> createManager() const;
    
    // Query configuration
    bool isNuclearTest() const { return is_nuclear; }
    bool isImpact() const { return is_impact; }
    double getYield() const { return yield_kt; }
    double getImpactorDiameter() const { return impactor_diameter; }
    
private:
    // Nuclear parameters
    bool is_nuclear = false;
    double yield_kt = 0.0;
    double burst_height = 0.0;
    double depth_of_burial = 0.0;
    BurstType burst_type;
    double fission_fraction = 0.5;
    
    // Impact parameters
    bool is_impact = false;
    double impactor_diameter = 0.0;
    double impactor_density = 3000.0;
    double impact_velocity = 20000.0;
    double impact_angle = 45.0;
    ImpactorType impactor_type;
    
    // Physics enables
    bool enable_radiation = true;
    bool enable_emp = true;
    bool enable_atmosphere = true;
    bool enable_thermal = true;
    bool enable_seismic = true;
    
    // Parsing helpers
    void parseNuclearSection(const std::map<std::string, std::string>& config);
    void parseImpactSection(const std::map<std::string, std::string>& config);
    void parsePhysicsSection(const std::map<std::string, std::string>& config);
};

} // namespace FSRM

#endif // EXPLOSION_IMPACT_PHYSICS_HPP
