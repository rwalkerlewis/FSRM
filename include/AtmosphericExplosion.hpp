/**
 * @file AtmosphericExplosion.hpp
 * @brief Comprehensive atmospheric explosion physics with instabilities and cloud dynamics
 * 
 * Implements advanced physics models for atmospheric explosions including:
 * - Rayleigh-Taylor and Richtmyer-Meshkov instabilities
 * - Kelvin-Helmholtz instabilities in shear layers
 * - Mushroom cloud formation and dynamics
 * - Buoyant thermal plume physics
 * - Enhanced atmospheric stratification effects
 * - Multi-phase flow for debris and particulates
 * - Turbulent mixing and entrainment
 * 
 * Physical basis and references:
 * - Rayleigh (1883) - Original RT instability analysis
 * - Taylor (1950) - RT instability growth rates
 * - Richtmyer (1960), Meshkov (1969) - Shock-driven instabilities
 * - Morton, Taylor & Turner (1956) - Buoyant plume theory
 * - Glasstone & Dolan (1977) - Nuclear weapons effects
 * - Brode (1955, 1959) - Fireball and blast physics
 * - Penner et al. (1986) - Atmospheric nuclear test phenomenology
 * - Carrier et al. (1985) - Convective column models
 * 
 * @author FSRM Development Team
 */

#ifndef ATMOSPHERIC_EXPLOSION_HPP
#define ATMOSPHERIC_EXPLOSION_HPP

#include <vector>
#include <array>
#include <string>
#include <memory>
#include <functional>
#include <cmath>
#include <complex>
#include <map>
#include <random>

namespace FSRM {

// Forward declarations
class AtmosphericProfile;
class TopographyModel;

// =============================================================================
// Physical Constants
// =============================================================================
namespace AtmosphericExplosionConstants {
    // Thermodynamic constants
    constexpr double R_AIR = 287.05;                      // J/(kg·K) gas constant
    constexpr double CP_AIR = 1005.0;                     // J/(kg·K) specific heat at constant pressure
    constexpr double CV_AIR = 718.0;                      // J/(kg·K) specific heat at constant volume
    constexpr double GAMMA_AIR = 1.4;                     // Ratio of specific heats
    constexpr double STEFAN_BOLTZMANN = 5.670374e-8;      // W/(m²·K⁴)
    constexpr double GRAVITY = 9.80665;                   // m/s²
    
    // Standard atmosphere
    constexpr double SEA_LEVEL_PRESSURE = 101325.0;       // Pa
    constexpr double SEA_LEVEL_TEMPERATURE = 288.15;      // K
    constexpr double SEA_LEVEL_DENSITY = 1.225;           // kg/m³
    constexpr double SCALE_HEIGHT = 8500.0;               // m
    constexpr double TROPOPAUSE_HEIGHT = 11000.0;         // m
    constexpr double STRATOPAUSE_HEIGHT = 50000.0;        // m
    
    // Nuclear explosion
    constexpr double JOULES_PER_KT = 4.184e12;            // J/kt TNT
    
    // Instability parameters
    constexpr double ATWOOD_NUMBER_THRESHOLD = 0.01;      // Minimum Atwood for RT
    constexpr double RT_AMPLITUDE_SEED = 0.001;           // Initial perturbation amplitude (fraction)
    constexpr double TURBULENT_PRANDTL = 0.9;             // Turbulent Prandtl number
    constexpr double ENTRAINMENT_COEFFICIENT = 0.12;      // Morton-Taylor-Turner
}

// =============================================================================
// Enumerations
// =============================================================================

/**
 * @brief Instability type enumeration
 */
enum class InstabilityType {
    RAYLEIGH_TAYLOR,              // Density-driven (heavy over light)
    RICHTMYER_MESHKOV,            // Shock-accelerated interface
    KELVIN_HELMHOLTZ,             // Shear-driven
    RAYLEIGH_BENARD,              // Thermal convection
    COMBINED                       // Multiple instability modes
};

/**
 * @brief Cloud evolution phase
 */
enum class CloudPhase {
    FIREBALL,                     // Initial fireball formation
    RISE,                         // Buoyant rise phase
    TOROIDAL_FORMATION,           // Vortex ring formation
    STEM_DEVELOPMENT,             // Dust stem entrainment
    STABILIZATION,                // Neutral buoyancy reached
    DISPERSION                    // Horizontal spreading
};

/**
 * @brief Entrainment model type
 */
enum class EntrainmentModel {
    MORTON_TAYLOR_TURNER,         // Classic MTT model
    PRIESTLEY_BALL,               // Priestley-Ball for strong plumes
    BRIGGS,                       // Briggs bent-over plume
    NUMERICAL_LES,                // Large eddy simulation
    HYBRID                         // Combined analytical/numerical
};

// =============================================================================
// Rayleigh-Taylor Instability Model
// =============================================================================

/**
 * @brief Rayleigh-Taylor instability parameters and evolution
 * 
 * Models the growth of perturbations at a density interface under acceleration.
 * 
 * Linear growth rate: σ = √(A·g·k) where:
 *   A = Atwood number = (ρ_heavy - ρ_light)/(ρ_heavy + ρ_light)
 *   g = effective acceleration
 *   k = wavenumber
 * 
 * Nonlinear growth: h = α·A·g·t² (bubble/spike penetration)
 */
struct RayleighTaylorInstability {
    // Interface properties
    double density_heavy = 2.0;           // kg/m³ (denser fluid)
    double density_light = 1.0;           // kg/m³ (lighter fluid)
    double interface_position = 0.0;      // m (initial interface height)
    
    // Perturbation parameters
    double initial_amplitude = 0.001;     // m (initial perturbation amplitude)
    double dominant_wavelength = 100.0;   // m (most unstable mode)
    double spectral_width = 0.5;          // Spectral bandwidth parameter
    bool multimode = true;                // Multiple wavelength modes
    
    // Growth parameters (empirically calibrated)
    double alpha_bubble = 0.05;           // Bubble penetration coefficient
    double alpha_spike = 0.06;            // Spike penetration coefficient
    
    // Turbulent mixing zone
    double mixing_efficiency = 0.5;       // Fraction of PE converted to mixing
    
    // Computed quantities
    double atwood_number() const {
        return (density_heavy - density_light) / (density_heavy + density_light);
    }
    
    double linear_growth_rate(double wavenumber, double acceleration) const {
        double A = atwood_number();
        if (A <= 0) return 0.0;
        return std::sqrt(A * acceleration * wavenumber);
    }
    
    double most_unstable_wavenumber(double thickness) const {
        // Most unstable mode at wavelength ~ interface thickness
        return 2.0 * M_PI / (4.0 * thickness);
    }
    
    double bubble_penetration(double time, double acceleration) const {
        double A = atwood_number();
        if (A <= 0) return 0.0;
        return alpha_bubble * A * acceleration * time * time;
    }
    
    double spike_penetration(double time, double acceleration) const {
        double A = atwood_number();
        if (A <= 0) return 0.0;
        return alpha_spike * A * acceleration * time * time;
    }
    
    double mixing_zone_width(double time, double acceleration) const {
        return bubble_penetration(time, acceleration) + spike_penetration(time, acceleration);
    }
};

/**
 * @brief Richtmyer-Meshkov instability for shock-accelerated interfaces
 * 
 * Impulsive analog of RT - interface accelerated by shock passage.
 * Growth rate: dh/dt = Δv · A · k · h₀ (linear phase)
 * Late time: h ~ t^θ where θ ≈ 0.2-0.4
 */
struct RichtmyerMeshkovInstability {
    // Shock parameters
    double shock_mach_number = 5.0;       // Incident shock Mach number
    double shock_velocity = 0.0;          // m/s (computed from Mach)
    double post_shock_velocity = 0.0;     // m/s (interface velocity after shock)
    
    // Interface parameters
    double pre_shock_atwood = 0.5;        // Atwood number before shock
    double post_shock_atwood = 0.0;       // Atwood number after (computed)
    double initial_amplitude = 0.01;      // m
    double wavenumber = 0.1;              // 1/m
    
    // Growth exponent
    double theta = 0.3;                   // Late-time growth exponent
    
    // Compute post-shock conditions
    void computePostShockState(double gamma = 1.4) {
        double M = shock_mach_number;
        double M2 = M * M;
        
        // Density ratio across shock (normal shock relations)
        double density_ratio = (gamma + 1.0) * M2 / ((gamma - 1.0) * M2 + 2.0);
        
        // Post-shock Atwood (compressed heavy fluid)
        post_shock_atwood = pre_shock_atwood * density_ratio / 
                           (1.0 + pre_shock_atwood * (density_ratio - 1.0));
    }
    
    double initial_growth_rate() const {
        // Impulsive growth rate from Richtmyer's theory
        return post_shock_velocity * post_shock_atwood * wavenumber * initial_amplitude;
    }
    
    double amplitude(double time) const {
        if (time <= 0) return initial_amplitude;
        
        // Linear phase for early time
        double t_nonlinear = 1.0 / (post_shock_atwood * wavenumber * post_shock_velocity);
        
        if (time < t_nonlinear) {
            return initial_amplitude + initial_growth_rate() * time;
        } else {
            // Late-time power law growth
            double h_transition = initial_amplitude + initial_growth_rate() * t_nonlinear;
            return h_transition * std::pow(time / t_nonlinear, theta);
        }
    }
};

/**
 * @brief Kelvin-Helmholtz instability for shear layers
 * 
 * Occurs at interfaces with velocity shear.
 * Growth rate: σ = k · Δu / 2 (inviscid, equal density)
 * With density contrast: modified by Richardson number
 */
struct KelvinHelmholtzInstability {
    // Shear parameters
    double velocity_upper = 100.0;        // m/s
    double velocity_lower = 0.0;          // m/s
    double shear_layer_thickness = 10.0;  // m
    
    // Density parameters
    double density_upper = 1.0;           // kg/m³
    double density_lower = 1.2;           // kg/m³
    
    // Stratification
    double buoyancy_frequency = 0.01;     // rad/s (N = sqrt(g/θ dθ/dz))
    
    double velocity_difference() const {
        return velocity_upper - velocity_lower;
    }
    
    double gradient_richardson_number() const {
        // Ri = N² / (du/dz)²
        double shear = velocity_difference() / shear_layer_thickness;
        if (std::abs(shear) < 1e-10) return 1e10;
        return buoyancy_frequency * buoyancy_frequency / (shear * shear);
    }
    
    bool is_unstable() const {
        // KH unstable when Ri < 0.25
        return gradient_richardson_number() < 0.25;
    }
    
    double growth_rate(double wavenumber) const {
        double Ri = gradient_richardson_number();
        double du = velocity_difference();
        
        if (Ri >= 0.25) return 0.0;  // Stable
        
        // Growth rate with stratification
        return 0.5 * wavenumber * du * std::sqrt(1.0 - 4.0 * Ri);
    }
    
    double most_unstable_wavelength() const {
        // Most unstable at wavelength ~ 7-10 times shear layer thickness
        return 7.0 * shear_layer_thickness;
    }
};

// =============================================================================
// Enhanced Atmospheric Model
// =============================================================================

/**
 * @brief Extended atmospheric structure model
 * 
 * Includes detailed vertical profiles for:
 * - Temperature, pressure, density
 * - Wind velocity and shear
 * - Humidity and condensation
 * - Turbulence intensity
 */
class ExtendedAtmosphericModel {
public:
    ExtendedAtmosphericModel();
    
    // Initialization
    void setStandardAtmosphere();
    void loadFromProfile(const std::string& filename);
    void setCustomProfile(const std::vector<double>& altitudes,
                         const std::vector<double>& temperatures,
                         const std::vector<double>& pressures);
    
    // Wind profile
    void setWindProfile(std::function<void(double z, double& u, double& v)> wind_func);
    void setLogWindProfile(double u_ref, double z_ref, double z0);
    void setJetStream(double jet_altitude, double jet_speed, double jet_width);
    
    // Basic properties
    double temperature(double z) const;
    double pressure(double z) const;
    double density(double z) const;
    double soundSpeed(double z) const;
    
    // Derived quantities
    double potentialTemperature(double z) const;
    double virtualPotentialTemperature(double z, double humidity) const;
    double buoyancyFrequency(double z) const;     // Brunt-Väisälä frequency
    double absoluteVorticity(double z, double latitude) const;
    
    // Stability
    double richardsonNumber(double z) const;
    bool isStaticallyStable(double z) const;
    bool isDynamicallyStable(double z) const;
    double convectiveAvailablePotentialEnergy(double z_start, double z_end) const;
    
    // Wind
    void wind(double z, double& u, double& v) const;
    double windShear(double z) const;
    double windDirection(double z) const;
    
    // Turbulence
    double turbulentKineticEnergy(double z) const;
    double turbulentDissipationRate(double z) const;
    double eddyViscosity(double z) const;
    double eddyDiffusivity(double z) const;
    
    // Layer boundaries
    double tropopauseHeight() const;
    double stratopauseHeight() const;
    double boundaryLayerHeight() const;
    
    // Humidity and condensation
    void setHumidityProfile(std::function<double(double)> humidity_func);
    double relativeHumidity(double z) const;
    double saturationVaporPressure(double T) const;
    double mixingRatio(double z) const;
    double liftingCondensationLevel(double z_surface, double T_surface, double RH_surface) const;
    
private:
    // Altitude breakpoints for standard atmosphere
    std::vector<double> altitude_levels_;
    std::vector<double> temperature_levels_;
    std::vector<double> lapse_rates_;
    
    // Wind profile function
    std::function<void(double, double&, double&)> wind_func_;
    
    // Humidity function
    std::function<double(double)> humidity_func_;
    
    // Jet stream parameters
    bool has_jet_stream_ = false;
    double jet_altitude_ = 10000.0;
    double jet_speed_ = 50.0;
    double jet_width_ = 3000.0;
    
    // Interpolation helpers
    double interpolateTemperature(double z) const;
    double computePressure(double z) const;
};

// =============================================================================
// Mushroom Cloud Physics
// =============================================================================

/**
 * @brief Comprehensive mushroom cloud model
 * 
 * Models the complete evolution of a nuclear mushroom cloud including:
 * - Initial fireball rise
 * - Buoyant thermal with vortex ring
 * - Dust stem entrainment
 * - Cap evolution and stabilization
 * - Rayleigh-Taylor instabilities at cloud boundaries
 */
class MushroomCloudModel {
public:
    MushroomCloudModel();
    
    // Configuration
    void setYield(double yield_kt);
    void setBurstHeight(double height_m);
    void setAtmosphere(std::shared_ptr<ExtendedAtmosphericModel> atm);
    void setGroundProperties(double dust_availability, double soil_density);
    
    // Enable physics
    void enableRayleighTaylor(bool enable);
    void enableKelvinHelmholtz(bool enable);
    void enableEntrainment(bool enable);
    void enableCondensation(bool enable);
    
    // Time evolution
    void initialize();
    void step(double dt);
    double getCurrentTime() const { return current_time_; }
    
    // Cloud state queries
    CloudPhase getCurrentPhase() const;
    
    // Geometry
    double cloudTopHeight() const;          // m
    double cloudBottomHeight() const;       // m
    double cloudCenterHeight() const;       // m
    double cloudRadius() const;             // m (cap radius)
    double stemRadius() const;              // m
    double stemHeight() const;              // m
    
    // Dynamics
    double riseVelocity() const;            // m/s
    double rotationalVelocity() const;      // m/s (toroidal)
    double entrainmentRate() const;         // kg/s
    double buoyancyFlux() const;            // m⁴/s³
    
    // Thermodynamics
    double cloudTemperature() const;        // K
    double cloudDensity() const;            // kg/m³
    double thermalExcess() const;           // K above ambient
    double buoyancy() const;                // m/s²
    
    // Instability characteristics
    double rtMixingZoneWidth() const;       // m (RT mixing layer)
    double khRollupScale() const;           // m (KH vortex size)
    double turbulentIntensity() const;      // m/s (rms velocity)
    
    // Mass and composition
    double totalMass() const;               // kg
    double debrisMass() const;              // kg (entrained material)
    double fissionProductMass() const;      // kg
    double waterVaporMass() const;          // kg (condensed)
    
    // Time markers
    double fireballMaxTime() const;         // s
    double riseStartTime() const;           // s
    double stabilizationTime() const;       // s
    double stabilizationHeight() const;     // m
    
    // Detailed profile
    void getVerticalProfile(std::vector<double>& z,
                           std::vector<double>& radius,
                           std::vector<double>& temperature,
                           std::vector<double>& density,
                           std::vector<double>& velocity) const;
    
    // Cross-section at height
    void getCrossSection(double z, double& radius, double& temperature,
                        double& w_velocity, double& mixing_ratio) const;
    
    // Instability state
    const RayleighTaylorInstability& getRTState() const { return rt_state_; }
    const KelvinHelmholtzInstability& getKHState() const { return kh_state_; }
    
private:
    // Configuration
    double yield_kt_;
    double burst_height_;
    std::shared_ptr<ExtendedAtmosphericModel> atmosphere_;
    
    // Physics enables
    bool enable_rt_ = true;
    bool enable_kh_ = true;
    bool enable_entrainment_ = true;
    bool enable_condensation_ = true;
    
    // Ground properties
    double dust_availability_ = 0.1;      // Fraction of surface material available
    double soil_density_ = 1500.0;        // kg/m³
    
    // State variables
    double current_time_ = 0.0;
    CloudPhase current_phase_ = CloudPhase::FIREBALL;
    
    // Cloud state
    double cloud_center_z_ = 0.0;
    double cloud_radius_ = 0.0;
    double cloud_temperature_ = 0.0;
    double cloud_density_ = 0.0;
    double rise_velocity_ = 0.0;
    double cloud_mass_ = 0.0;
    double debris_mass_ = 0.0;
    
    // Stem state
    double stem_radius_ = 0.0;
    double stem_top_z_ = 0.0;
    
    // Toroidal vortex
    double vortex_strength_ = 0.0;        // Circulation (m²/s)
    double vortex_core_radius_ = 0.0;     // m
    
    // Instability states
    RayleighTaylorInstability rt_state_;
    KelvinHelmholtzInstability kh_state_;
    
    // Internal methods
    void updateFireballPhase(double dt);
    void updateRisePhase(double dt);
    void updateStabilizationPhase(double dt);
    void updateDispersionPhase(double dt);
    
    void computeEntrainment(double dt);
    void computeInstabilities(double dt);
    void computeCondensation(double dt);
    void updateStem(double dt);
    
    double computeBuoyancy() const;
    double computeDragForce() const;
    double effectiveGravity() const;
    
    // Scaling relations
    double initialFireballRadius() const;
    double maxFireballRadius() const;
    double initialCloudTemperature() const;
    double initialCloudMass() const;
};

// =============================================================================
// Buoyant Plume Model (Morton-Taylor-Turner)
// =============================================================================

/**
 * @brief General buoyant plume model
 * 
 * Extends MTT theory with:
 * - Variable stratification
 * - Cross-flow effects
 * - Moisture and condensation
 * - Turbulent structure
 */
class BuoyantPlumeModel {
public:
    BuoyantPlumeModel();
    
    // Source configuration
    void setSourceStrength(double buoyancy_flux);    // m⁴/s³
    void setSourceMomentum(double momentum_flux);    // m⁴/s²
    void setSourceRadius(double radius);             // m
    void setSourceHeight(double height);             // m
    void setSourceTemperature(double temperature);   // K
    
    void setAtmosphere(std::shared_ptr<ExtendedAtmosphericModel> atm);
    void setEntrainmentModel(EntrainmentModel model);
    
    // Solve for plume trajectory and properties
    void solve();
    
    // Results at height z
    double plumeRadius(double z) const;
    double plumeVelocity(double z) const;
    double plumeTemperature(double z) const;
    double plumeDensity(double z) const;
    double plumeMassFlux(double z) const;
    double plumeBuoyancyFlux(double z) const;
    
    // Key heights
    double neutralBuoyancyHeight() const;
    double maximumRiseHeight() const;
    double spreadingHeight() const;
    
    // Bent-over plume (with crosswind)
    void setCrosswind(double speed, double direction);
    double plumeTrajectoryX(double z) const;
    double plumeTrajectoryY(double z) const;
    
private:
    // Source parameters
    double source_buoyancy_flux_;
    double source_momentum_flux_;
    double source_radius_;
    double source_height_;
    double source_temperature_;
    
    // Atmosphere
    std::shared_ptr<ExtendedAtmosphericModel> atmosphere_;
    EntrainmentModel entrainment_model_;
    
    // Crosswind
    double crosswind_speed_ = 0.0;
    double crosswind_direction_ = 0.0;
    
    // Solution storage
    std::vector<double> z_profile_;
    std::vector<double> radius_profile_;
    std::vector<double> velocity_profile_;
    std::vector<double> temperature_profile_;
    std::vector<double> x_trajectory_;
    std::vector<double> y_trajectory_;
    
    double neutral_height_;
    double max_height_;
    
    // Entrainment coefficient
    double entrainmentCoefficient(double Ri_local) const;
    
    // ODE system for plume integration
    void plumeEquations(double z, const std::vector<double>& state,
                       std::vector<double>& dstate_dz) const;
};

// =============================================================================
// Multi-Phase Debris Model
// =============================================================================

/**
 * @brief Multi-phase model for particulates and debris
 * 
 * Tracks multiple size classes of particles:
 * - Large debris (ballistic)
 * - Medium particles (settling with entrainment)
 * - Fine particulates (suspended in cloud)
 * - Condensed water droplets
 */
class DebrisTransportModel {
public:
    /**
     * @brief Particle size class
     */
    struct ParticleSizeClass {
        double diameter_min;           // m
        double diameter_max;           // m
        double median_diameter;        // m
        double density;                // kg/m³
        double mass_fraction;          // Fraction of total debris
        double settling_velocity;      // m/s (terminal)
        
        // Radioactive properties
        double specific_activity;      // Bq/kg
        double gamma_energy;           // MeV (mean)
    };
    
    DebrisTransportModel();
    
    // Configuration
    void setInitialDebrisMass(double mass);
    void setParticleSizeDistribution(const std::vector<ParticleSizeClass>& classes);
    void setDefaultNuclearDistribution(double yield_kt, bool surface_burst);
    
    void setAtmosphere(std::shared_ptr<ExtendedAtmosphericModel> atm);
    void setCloudModel(std::shared_ptr<MushroomCloudModel> cloud);
    
    // Time evolution
    void initialize();
    void step(double dt);
    
    // Particle tracking
    void getParticleDistribution(double z, std::vector<double>& mass_density) const;
    double totalSuspendedMass() const;
    double massSettledBelow(double z) const;
    
    // Fallout from particles
    double falloutRate(double x, double y, double t) const;  // kg/(m²·s)
    double accumulatedFallout(double x, double y, double t) const;  // kg/m²
    
    // Activity
    double activityConcentration(double x, double y, double z, double t) const;  // Bq/m³
    double groundActivity(double x, double y, double t) const;  // Bq/m²
    
private:
    std::vector<ParticleSizeClass> size_classes_;
    double total_debris_mass_;
    
    std::shared_ptr<ExtendedAtmosphericModel> atmosphere_;
    std::shared_ptr<MushroomCloudModel> cloud_;
    
    // Particle state (mass in each size class at each height)
    std::vector<std::vector<double>> particle_mass_;  // [class][height_bin]
    std::vector<double> height_bins_;
    
    // Ground deposition
    std::vector<std::vector<double>> ground_deposition_;  // [x][y]
    std::vector<double> x_grid_, y_grid_;
    
    double current_time_ = 0.0;
    
    // Physics
    double computeSettlingVelocity(const ParticleSizeClass& pc, double z) const;
    double computeEntrainmentFlux(double z, double cloud_velocity) const;
};

// =============================================================================
// Turbulent Mixing Model
// =============================================================================

/**
 * @brief Turbulent mixing and entrainment model
 * 
 * Provides subgrid-scale models for:
 * - Turbulent viscosity and diffusivity
 * - Entrainment at cloud boundaries
 * - Mixing efficiency in unstable layers
 */
class TurbulentMixingModel {
public:
    TurbulentMixingModel();
    
    // Mixing length models
    double mixingLength(double z, double boundary_layer_height) const;
    double smagorinskyViscosity(double strain_rate, double grid_size) const;
    
    // Turbulent Prandtl number
    double turbulentPrandtl(double Ri) const;
    
    // Entrainment velocity
    double entrainmentVelocity(double w, double Ri_bulk) const;
    double detrainmentVelocity(double w, double Ri_bulk) const;
    
    // Mixing in unstable layers
    double mixingEfficiency(double Ri, double Pe) const;
    double molecularMixednessFraction(double time, double thickness, double diffusivity) const;
    
    // Turbulent kinetic energy budget
    double shearProduction(double viscosity, double shear) const;
    double buoyantProduction(double diffusivity, double N2, double w_rms) const;
    double dissipationRate(double k, double length_scale) const;
    
    // Turbulent spectra
    double kolmogorovWavenumber(double epsilon, double nu) const;
    double batchelorWavenumber(double epsilon, double nu, double D) const;
    
private:
    double karman_constant_ = 0.41;
    double smagorinsky_constant_ = 0.18;
};

// =============================================================================
// Atmospheric Explosion Solver
// =============================================================================

/**
 * @brief Complete solver for atmospheric explosion physics
 * 
 * Integrates:
 * - Fireball and thermal radiation
 * - Blast wave propagation
 * - Mushroom cloud dynamics
 * - Instability evolution
 * - Debris transport
 * - Fallout generation
 */
class AtmosphericExplosionSolver {
public:
    AtmosphericExplosionSolver();
    ~AtmosphericExplosionSolver() = default;
    
    // Configuration
    void setYield(double yield_kt);
    void setBurstHeight(double height_m);
    void setBurstLocation(double x, double y, double z);
    void setFissionFraction(double fraction);
    
    void setAtmosphere(std::shared_ptr<ExtendedAtmosphericModel> atm);
    void setTopography(std::shared_ptr<TopographyModel> topo);
    
    // Enable physics modules
    void enableMushroomCloud(bool enable);
    void enableInstabilities(bool enable);
    void enableDebrisTransport(bool enable);
    void enableTurbulentMixing(bool enable);
    void enableRadiation(bool enable);
    
    // Initialize and run
    void initialize();
    void step(double dt);
    void run(double end_time, double dt);
    
    double getCurrentTime() const { return current_time_; }
    
    // Access sub-models
    std::shared_ptr<MushroomCloudModel> getCloudModel() { return cloud_model_; }
    std::shared_ptr<DebrisTransportModel> getDebrisModel() { return debris_model_; }
    std::shared_ptr<TurbulentMixingModel> getMixingModel() { return mixing_model_; }
    
    // Output
    void writeOutput(const std::string& filename) const;
    void writeCloudProfile(const std::string& filename) const;
    void writeFalloutMap(const std::string& filename) const;
    
    // Instability diagnostics
    void getInstabilityGrowthHistory(std::vector<double>& times,
                                     std::vector<double>& rt_amplitude,
                                     std::vector<double>& kh_amplitude) const;
    
private:
    // Configuration
    double yield_kt_ = 100.0;
    double burst_height_ = 0.0;
    double burst_x_ = 0.0, burst_y_ = 0.0, burst_z_ = 0.0;
    double fission_fraction_ = 0.5;
    
    // Models
    std::shared_ptr<ExtendedAtmosphericModel> atmosphere_;
    std::shared_ptr<TopographyModel> topography_;
    std::shared_ptr<MushroomCloudModel> cloud_model_;
    std::shared_ptr<DebrisTransportModel> debris_model_;
    std::shared_ptr<TurbulentMixingModel> mixing_model_;
    
    // Physics enables
    bool enable_cloud_ = true;
    bool enable_instabilities_ = true;
    bool enable_debris_ = true;
    bool enable_mixing_ = true;
    bool enable_radiation_ = true;
    
    // State
    double current_time_ = 0.0;
    bool initialized_ = false;
    
    // History tracking
    std::vector<double> time_history_;
    std::vector<double> rt_amplitude_history_;
    std::vector<double> kh_amplitude_history_;
    std::vector<double> cloud_height_history_;
    
    // Internal methods
    void updateCoupledPhysics(double dt);
    void recordHistory();
};

// =============================================================================
// Configuration Parser
// =============================================================================

/**
 * @brief Parse atmospheric explosion configuration
 */
class AtmosphericExplosionConfig {
public:
    AtmosphericExplosionConfig() = default;
    
    bool parse(const std::string& filename);
    bool parse(const std::map<std::string, std::string>& config);
    
    std::unique_ptr<AtmosphericExplosionSolver> createSolver() const;
    
    // Access configuration
    double getYield() const { return yield_kt_; }
    double getBurstHeight() const { return burst_height_; }
    bool hasInstabilities() const { return enable_instabilities_; }
    
private:
    double yield_kt_ = 100.0;
    double burst_height_ = 0.0;
    double burst_x_ = 0.0, burst_y_ = 0.0, burst_z_ = 0.0;
    double fission_fraction_ = 0.5;
    
    bool enable_cloud_ = true;
    bool enable_instabilities_ = true;
    bool enable_debris_ = true;
    bool enable_mixing_ = true;
    bool enable_radiation_ = true;
    
    std::string atmosphere_model_ = "US_STANDARD_1976";
    
    bool configured_ = false;
};

} // namespace FSRM

#endif // ATMOSPHERIC_EXPLOSION_HPP
