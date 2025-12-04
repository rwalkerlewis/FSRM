/**
 * @file AtmosphericInfrasound.hpp
 * @brief Atmospheric physics and infrasound modeling for FSRM
 * 
 * Implements physics models for:
 * - Atmospheric acoustic wave propagation
 * - Infrasound (sub-audible sound < 20 Hz)
 * - Topographic effects on wave propagation
 * - Wind and temperature stratification effects
 * - Atmospheric attenuation and absorption
 * 
 * All calculations use SI base units (m, kg, s) internally.
 * User can specify input/output in any supported units.
 * 
 * @note Infrasound applications include:
 * - Nuclear explosion monitoring (CTBT/IMS)
 * - Volcanic eruption detection
 * - Meteorite/bolide tracking
 * - Severe weather monitoring
 * - Industrial accident detection
 */

#ifndef ATMOSPHERIC_INFRASOUND_HPP
#define ATMOSPHERIC_INFRASOUND_HPP

#include "PhysicsKernel.hpp"
#include "UnitSystem.hpp"
#include <vector>
#include <array>
#include <complex>
#include <functional>
#include <memory>
#include <cmath>

namespace FSRM {

// =============================================================================
// Physical Constants (SI base units)
// =============================================================================

namespace AtmosphericConstants {
    constexpr double GAMMA_AIR = 1.4;                    // Ratio of specific heats
    constexpr double R_AIR = 287.05;                     // Gas constant for air (J/kg/K)
    constexpr double SEA_LEVEL_PRESSURE = 101325.0;      // Pa
    constexpr double SEA_LEVEL_TEMPERATURE = 288.15;     // K (15°C)
    constexpr double SEA_LEVEL_DENSITY = 1.225;          // kg/m³
    constexpr double SCALE_HEIGHT = 8500.0;              // m
    constexpr double EARTH_RADIUS = 6371000.0;           // m
    constexpr double GRAVITY = 9.80665;                  // m/s²
    constexpr double SOUND_SPEED_0 = 343.0;              // m/s at 20°C
    
    // Molecular absorption constants
    constexpr double OXYGEN_RELAXATION_FREQ = 24.0;      // Hz at 20°C, 1 atm
    constexpr double NITROGEN_RELAXATION_FREQ = 9.0;     // Hz at 20°C, 1 atm
}

// =============================================================================
// Atmospheric Profile Models
// =============================================================================

/**
 * @brief Standard atmosphere types
 */
enum class AtmosphereModel {
    US_STANDARD_1976,        ///< US Standard Atmosphere 1976
    ICAO_STANDARD,           ///< ICAO International Standard
    TROPICAL,                ///< Tropical atmosphere
    MIDLATITUDE_SUMMER,      ///< Mid-latitude summer
    MIDLATITUDE_WINTER,      ///< Mid-latitude winter
    SUBARCTIC_SUMMER,        ///< Sub-arctic summer
    SUBARCTIC_WINTER,        ///< Sub-arctic winter
    CUSTOM                   ///< User-defined profile
};

/**
 * @brief Atmospheric state at a given altitude
 */
struct AtmosphericState {
    double altitude;         ///< Altitude above sea level (m)
    double temperature;      ///< Temperature (K)
    double pressure;         ///< Pressure (Pa)
    double density;          ///< Density (kg/m³)
    double sound_speed;      ///< Sound speed (m/s)
    double wind_east;        ///< Eastward wind component (m/s)
    double wind_north;       ///< Northward wind component (m/s)
    double humidity;         ///< Relative humidity (0-1)
    
    /// Compute sound speed from temperature
    static double soundSpeedFromTemp(double T) {
        return std::sqrt(AtmosphericConstants::GAMMA_AIR * 
                        AtmosphericConstants::R_AIR * T);
    }
    
    /// Effective sound speed along propagation direction
    double effectiveSoundSpeed(double azimuth_rad) const {
        double wind_along = wind_east * std::sin(azimuth_rad) + 
                           wind_north * std::cos(azimuth_rad);
        return sound_speed + wind_along;
    }
};

/**
 * @brief Atmospheric profile manager
 * 
 * Provides atmospheric properties as a function of altitude,
 * including temperature, pressure, wind, and derived quantities.
 */
class AtmosphericProfile {
public:
    AtmosphericProfile(AtmosphereModel model = AtmosphereModel::US_STANDARD_1976);
    
    /// Set standard atmosphere model
    void setModel(AtmosphereModel model);
    
    /// Set custom profile from data arrays
    void setCustomProfile(const std::vector<double>& altitudes,
                         const std::vector<double>& temperatures,
                         const std::vector<double>& pressures);
    
    /// Set wind profile
    void setWindProfile(std::function<void(double z, double& u, double& v)> wind_func);
    
    /// Set humidity profile
    void setHumidityProfile(std::function<double(double z)> humidity_func);
    
    /// Get atmospheric state at altitude z (m)
    AtmosphericState getState(double z) const;
    
    /// Get temperature at altitude (K)
    double temperature(double z) const;
    
    /// Get pressure at altitude (Pa)
    double pressure(double z) const;
    
    /// Get density at altitude (kg/m³)
    double density(double z) const;
    
    /// Get sound speed at altitude (m/s)
    double soundSpeed(double z) const;
    
    /// Get wind velocity at altitude
    void wind(double z, double& u_east, double& v_north) const;
    
    /// Get effective sound speed along azimuth direction
    double effectiveSoundSpeed(double z, double azimuth_rad) const;
    
    /// Check for acoustic duct (waveguide)
    bool hasAcousticDuct(double z_min, double z_max, double azimuth_rad) const;
    
    /// Find turning point altitude for ray at given launch angle
    double turningPointAltitude(double z0, double launch_angle, double azimuth) const;
    
private:
    AtmosphereModel model_;
    std::vector<double> alt_data_;
    std::vector<double> temp_data_;
    std::vector<double> pres_data_;
    std::function<void(double, double&, double&)> wind_func_;
    std::function<double(double)> humidity_func_;
    
    // Standard atmosphere calculations
    double usStandard1976Temperature(double z) const;
    double usStandard1976Pressure(double z) const;
};

// =============================================================================
// Topography Model
// =============================================================================

/**
 * @brief Digital elevation model interface for topographic effects
 */
class TopographyModel {
public:
    TopographyModel();
    virtual ~TopographyModel() = default;
    
    /// Load topography from file (GeoTIFF, NetCDF, etc.)
    virtual bool loadFromFile(const std::string& filename);
    
    /// Set uniform flat terrain at given elevation
    void setFlatTerrain(double elevation_m);
    
    /// Set analytical terrain function
    void setTerrainFunction(std::function<double(double, double)> elev_func);
    
    /// Get terrain elevation at (x, y) in meters
    virtual double elevation(double x, double y) const;
    
    /// Get terrain slope (gradient) at (x, y)
    virtual void slope(double x, double y, double& dz_dx, double& dz_dy) const;
    
    /// Get terrain normal vector at (x, y)
    virtual void normal(double x, double y, double& nx, double& ny, double& nz) const;
    
    /// Check if point is above terrain
    bool isAboveTerrain(double x, double y, double z) const {
        return z > elevation(x, y);
    }
    
    /// Get terrain roughness parameter for acoustic scattering
    double roughness(double x, double y) const;
    
    /// Set roughness function
    void setRoughnessFunction(std::function<double(double, double)> rough_func);
    
protected:
    std::function<double(double, double)> elevation_func_;
    std::function<double(double, double)> roughness_func_;
    double default_elevation_ = 0.0;
    double default_roughness_ = 0.0;
};

// =============================================================================
// Infrasound Source Models
// =============================================================================

/**
 * @brief Infrasound source types
 */
enum class InfrasoundSourceType {
    POINT_MONOPOLE,          ///< Idealized point source
    EXPLOSION,               ///< Chemical/nuclear explosion
    VOLCANIC,                ///< Volcanic eruption
    BOLIDE,                  ///< Meteor/asteroid entry
    EARTHQUAKE,              ///< Seismo-acoustic coupling
    OCEAN_SWELL,             ///< Microbaroms from ocean waves
    SEVERE_WEATHER,          ///< Thunderstorm, tornado
    INDUSTRIAL,              ///< Industrial explosion/accident
    ROCKET_LAUNCH            ///< Rocket or missile launch
};

/**
 * @brief Infrasound source parameters
 */
struct InfrasoundSourceParams {
    InfrasoundSourceType type = InfrasoundSourceType::POINT_MONOPOLE;
    
    // Location (SI base units: meters)
    double x = 0.0;          ///< X coordinate (m)
    double y = 0.0;          ///< Y coordinate (m)
    double z = 0.0;          ///< Z coordinate (m) - altitude above sea level
    
    // Source strength
    double yield_kt = 0.0;           ///< Equivalent TNT yield (kt) for explosions
    double peak_overpressure = 0.0;  ///< Peak overpressure at source (Pa)
    double duration = 1.0;           ///< Source duration (s)
    double dominant_frequency = 1.0; ///< Dominant frequency (Hz)
    
    // Time
    double onset_time = 0.0;         ///< Event start time (s)
    
    // Compute source strength from yield using Reed's relation
    double computePeakOverpressure(double distance_m) const;
    
    // Compute dominant period from yield (empirical)
    double computeDominantPeriod() const;
};

/**
 * @brief Infrasound source model
 * 
 * Generates acoustic source terms for infrasound propagation.
 */
class InfrasoundSource {
public:
    InfrasoundSource();
    explicit InfrasoundSource(const InfrasoundSourceParams& params);
    
    void setParameters(const InfrasoundSourceParams& params);
    const InfrasoundSourceParams& getParameters() const { return params_; }
    
    /// Pressure source term at location (x, y, z) and time t
    double pressureSource(double x, double y, double z, double t) const;
    
    /// Velocity source term
    void velocitySource(double x, double y, double z, double t,
                       double& vx, double& vy, double& vz) const;
    
    /// Frequency spectrum (amplitude vs frequency)
    double spectrum(double frequency) const;
    
    /// Source time function (pressure rate vs time)
    double sourceTimeFunction(double t) const;
    
    /// Directivity pattern (for non-isotropic sources)
    double directivity(double theta, double phi) const;
    
    /// Get frequency range of significant energy
    void getFrequencyRange(double& f_min, double& f_max) const;
    
private:
    InfrasoundSourceParams params_;
    
    // Source time function types
    double explosionSTF(double t) const;
    double volcanicSTF(double t) const;
    double bolideSTF(double t) const;
};

// =============================================================================
// Atmospheric Attenuation
// =============================================================================

/**
 * @brief Atmospheric absorption model for infrasound
 * 
 * Computes frequency-dependent attenuation due to:
 * - Classical absorption (viscosity, thermal conductivity)
 * - Molecular relaxation (O2, N2)
 * - Humidity effects
 */
class AtmosphericAttenuation {
public:
    AtmosphericAttenuation();
    
    /// Set atmospheric conditions
    void setConditions(double temperature_K, double pressure_Pa, 
                      double humidity_fraction);
    
    /// Set conditions from atmospheric state
    void setConditions(const AtmosphericState& state);
    
    /// Compute absorption coefficient (Np/m) at frequency f
    double absorptionCoefficient(double frequency_Hz) const;
    
    /// Compute absorption in dB/km at frequency f
    double absorptionDB(double frequency_Hz) const;
    
    /// Compute total attenuation over distance
    double totalAttenuation(double frequency_Hz, double distance_m) const;
    
    /// Apply attenuation to complex amplitude
    std::complex<double> attenuate(std::complex<double> amplitude,
                                   double frequency_Hz, 
                                   double distance_m) const;
    
private:
    double temperature_;
    double pressure_;
    double humidity_;
    
    // Relaxation frequencies
    double f_ro_;  // O2 relaxation frequency
    double f_rn_;  // N2 relaxation frequency
    
    void computeRelaxationFrequencies();
    double classicalAbsorption(double f) const;
    double molecularAbsorption(double f) const;
};

// =============================================================================
// Ray Tracing for Infrasound
// =============================================================================

/**
 * @brief Ray path for geometric acoustics
 */
struct RayPath {
    std::vector<std::array<double, 3>> points;  ///< (x, y, z) along ray
    std::vector<double> times;                   ///< Travel time at each point
    std::vector<double> amplitudes;              ///< Amplitude at each point
    
    double total_distance = 0.0;
    double travel_time = 0.0;
    double final_amplitude = 1.0;
    int num_bounces = 0;
    bool reached_receiver = false;
};

/**
 * @brief Ray tracing for infrasound propagation
 * 
 * Uses geometric ray theory with atmospheric refraction.
 */
class InfrasoundRayTracer {
public:
    InfrasoundRayTracer();
    
    /// Set atmospheric profile
    void setAtmosphere(std::shared_ptr<AtmosphericProfile> atm);
    
    /// Set topography
    void setTopography(std::shared_ptr<TopographyModel> topo);
    
    /// Set source location
    void setSource(double x, double y, double z);
    
    /// Set receiver location
    void setReceiver(double x, double y, double z);
    
    /// Trace single ray with given launch angles
    RayPath traceRay(double elevation_angle, double azimuth) const;
    
    /// Find eigenrays connecting source to receiver
    std::vector<RayPath> findEigenrays(int max_bounces = 5, 
                                       double angle_tolerance = 0.1) const;
    
    /// Compute transmission loss to receiver
    double transmissionLoss(double frequency_Hz) const;
    
    /// Set maximum propagation distance
    void setMaxDistance(double d) { max_distance_ = d; }
    
    /// Set step size for ray integration
    void setStepSize(double ds) { step_size_ = ds; }
    
private:
    std::shared_ptr<AtmosphericProfile> atmosphere_;
    std::shared_ptr<TopographyModel> topography_;
    std::array<double, 3> source_;
    std::array<double, 3> receiver_;
    double max_distance_ = 1000000.0;  // 1000 km
    double step_size_ = 100.0;         // 100 m
    
    // Ray equations (Snell's law in 2D/3D)
    void rayStep(double& x, double& y, double& z,
                double& px, double& py, double& pz,
                double ds) const;
};

// =============================================================================
// Parabolic Equation Method
// =============================================================================

/**
 * @brief Parabolic equation (PE) solver for infrasound
 * 
 * Solves the one-way wave equation using split-step Fourier method.
 * More accurate than ray tracing for low frequencies and complex media.
 */
class ParabolicEquationSolver {
public:
    ParabolicEquationSolver();
    
    /// Set atmospheric profile
    void setAtmosphere(std::shared_ptr<AtmosphericProfile> atm);
    
    /// Set topography
    void setTopography(std::shared_ptr<TopographyModel> topo);
    
    /// Set source
    void setSource(const InfrasoundSource& src);
    
    /// Set frequency (Hz)
    void setFrequency(double f);
    
    /// Set grid parameters
    void setGrid(double dz, double dr, double z_max, double r_max);
    
    /// Run propagation calculation
    void propagate();
    
    /// Get pressure field at range r
    std::vector<std::complex<double>> getPressure(double r) const;
    
    /// Get transmission loss at (r, z)
    double transmissionLoss(double r, double z) const;
    
    /// Get pressure time series at receiver
    void getTimeSeries(double r, double z, 
                      std::vector<double>& times,
                      std::vector<double>& pressures) const;
    
private:
    std::shared_ptr<AtmosphericProfile> atmosphere_;
    std::shared_ptr<TopographyModel> topography_;
    InfrasoundSource source_;
    double frequency_ = 1.0;
    double dz_ = 50.0;        // Vertical step (m)
    double dr_ = 500.0;       // Range step (m)
    double z_max_ = 120000.0; // Max altitude (m)
    double r_max_ = 1000000.0; // Max range (m)
    
    std::vector<std::vector<std::complex<double>>> field_;
};

// =============================================================================
// Infrasound Physics Kernel
// =============================================================================

/**
 * @brief Physics kernel for infrasound propagation
 * 
 * Implements the linearized acoustic wave equation in a stratified,
 * moving atmosphere with topographic boundary conditions.
 * 
 * Governing equations:
 *   ∂p/∂t + ρc² ∇·v + v·∇p₀ = S
 *   ρ ∂v/∂t + ∇p = 0
 * 
 * where p is acoustic pressure perturbation, v is velocity perturbation,
 * ρ is density, c is sound speed, and S is source term.
 * 
 * All quantities are in SI base units.
 */
class InfrasoundKernel : public PhysicsKernel {
public:
    InfrasoundKernel();
    ~InfrasoundKernel() override = default;
    
    // PhysicsKernel interface
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    int getNumFields() const override { return 2; }  // pressure, velocity
    int getNumComponents(int field) const override { 
        return field == 0 ? 1 : 3; 
    }
    
    // Configuration
    void setAtmosphere(std::shared_ptr<AtmosphericProfile> atm);
    void setTopography(std::shared_ptr<TopographyModel> topo);
    void addSource(const InfrasoundSource& src);
    void setAttenuation(bool enable);
    
    /// Set frequency range for simulation
    void setFrequencyRange(double f_min, double f_max);
    
    /// Get atmospheric properties at location (for coupling)
    AtmosphericState getAtmosphericState(double x, double y, double z) const;
    
private:
    std::shared_ptr<AtmosphericProfile> atmosphere_;
    std::shared_ptr<TopographyModel> topography_;
    std::vector<InfrasoundSource> sources_;
    AtmosphericAttenuation attenuation_;
    bool use_attenuation_ = true;
    double f_min_ = 0.01;
    double f_max_ = 20.0;
    
    // Cached atmospheric properties at quadrature points
    mutable double cached_sound_speed_;
    mutable double cached_density_;
    mutable std::array<double, 3> cached_wind_;
};

// =============================================================================
// Atmospheric Blast Kernel
// =============================================================================

/**
 * @brief Physics kernel for atmospheric blast wave propagation
 * 
 * Solves the Euler equations for compressible flow with strong shocks.
 * Includes spherical divergence and atmospheric stratification.
 * 
 * Conservation equations:
 *   ∂ρ/∂t + ∇·(ρv) = 0
 *   ∂(ρv)/∂t + ∇·(ρv⊗v) + ∇p = ρg
 *   ∂E/∂t + ∇·((E+p)v) = ρv·g
 * 
 * where E = ρe + ½ρv² is total energy density.
 */
class AtmosphericBlastKernel : public PhysicsKernel {
public:
    AtmosphericBlastKernel();
    ~AtmosphericBlastKernel() override = default;
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    int getNumFields() const override { return 3; }  // ρ, ρv, E
    int getNumComponents(int field) const override;
    
    void setAtmosphere(std::shared_ptr<AtmosphericProfile> atm);
    void setExplosionSource(double yield_kt, double x, double y, double z);
    void setBurstHeight(double height_m);
    
    /// Get peak overpressure at distance
    double peakOverpressure(double distance) const;
    
    /// Get shock arrival time at distance
    double arrivalTime(double distance) const;
    
private:
    std::shared_ptr<AtmosphericProfile> atmosphere_;
    double yield_kt_ = 1.0;
    double burst_x_ = 0.0, burst_y_ = 0.0, burst_z_ = 0.0;
    double burst_height_ = 0.0;
    double gamma_ = 1.4;
};

// =============================================================================
// Coupled Physics Manager for Atmospheric Simulations
// =============================================================================

/**
 * @brief Manager for coupled atmospheric physics simulations
 * 
 * Coordinates:
 * - Source (explosion, volcanic, seismic)
 * - Near-field blast wave
 * - Far-field infrasound propagation
 * - Ground coupling
 * - Atmospheric refraction
 */
class AtmosphericPhysicsManager {
public:
    AtmosphericPhysicsManager();
    
    void setAtmosphere(std::shared_ptr<AtmosphericProfile> atm);
    void setTopography(std::shared_ptr<TopographyModel> topo);
    
    // Source configuration
    void configureExplosionSource(double yield_kt, double x, double y, double z);
    void configureVolcanicSource(double x, double y, double z, double vei);
    void configureSeismicSource(double x, double y, double z, double Mw);
    
    // Solver configuration
    void enableBlastWave(bool enable);
    void enableInfrasound(bool enable);
    void enableGroundCoupling(bool enable);
    
    // Get kernels for simulation
    std::shared_ptr<AtmosphericBlastKernel> getBlastKernel();
    std::shared_ptr<InfrasoundKernel> getInfrasoundKernel();
    
    // Compute infrasound signal at receiver
    void computeInfrasoundSignal(double rx, double ry, double rz,
                                std::vector<double>& times,
                                std::vector<double>& pressures);
    
    // Compute transmission loss map
    void computeTransmissionLossMap(double z_receiver,
                                   std::vector<double>& x_grid,
                                   std::vector<double>& y_grid,
                                   std::vector<std::vector<double>>& TL);
    
private:
    std::shared_ptr<AtmosphericProfile> atmosphere_;
    std::shared_ptr<TopographyModel> topography_;
    std::shared_ptr<AtmosphericBlastKernel> blast_kernel_;
    std::shared_ptr<InfrasoundKernel> infrasound_kernel_;
    
    bool enable_blast_ = true;
    bool enable_infrasound_ = true;
    bool enable_ground_coupling_ = true;
};

} // namespace FSRM

#endif // ATMOSPHERIC_INFRASOUND_HPP
