/**
 * @file HydrodynamicModel.hpp
 * @brief Comprehensive hydrodynamic modeling for ocean and coastal simulations
 *
 * Implements hydrodynamic solvers for various ocean and coastal applications:
 * - 3D primitive equations (ocean circulation)
 * - 2D shallow water equations (coastal hydrodynamics)
 * - Storm surge modeling with wind and pressure forcing
 * - Wave-current interaction
 * - Sediment transport and morphodynamics
 * - Estuarine and river-ocean coupling
 *
 * Numerical methods:
 * - Finite volume and finite difference discretizations
 * - Sigma/terrain-following coordinates for bathymetry
 * - Split-explicit time integration
 * - Wetting/drying for inundation
 *
 * Physical processes:
 * - Coriolis acceleration (f-plane, beta-plane)
 * - Wind and atmospheric pressure forcing
 * - Tidal forcing (boundary and body tides)
 * - Baroclinic pressure gradients
 * - Bottom friction (various formulations)
 * - Horizontal and vertical mixing
 *
 * @author FSRM Development Team
 */

#ifndef HYDRODYNAMIC_MODEL_HPP
#define HYDRODYNAMIC_MODEL_HPP

#include "FSRM.hpp"
#include "TsunamiModel.hpp"
#include "CoordinateSystem.hpp"
#include <vector>
#include <array>
#include <string>
#include <map>
#include <memory>
#include <functional>
#include <complex>

namespace FSRM {

// Forward declarations
class OceanPhysics;
class WaveModel;

// =============================================================================
// Enumerations
// =============================================================================

/**
 * @brief Hydrodynamic model type
 */
enum class HydrodynamicModelType {
    SHALLOW_WATER_2D,      // 2D depth-averaged shallow water equations
    PRIMITIVE_3D,          // 3D primitive equations
    BOUSSINESQ,            // Boussinesq equations with dispersion
    NON_HYDROSTATIC,       // Full 3D non-hydrostatic
    QUASI_GEOSTROPHIC      // Quasi-geostrophic approximation
};

/**
 * @brief Vertical coordinate system
 */
enum class VerticalCoordinateType {
    Z_LEVEL,               // Fixed z-levels
    SIGMA,                 // Terrain-following sigma coordinates
    HYBRID_Z_SIGMA,        // Hybrid z-sigma coordinates
    ISOPYCNAL,             // Density (isopycnal) coordinates
    ADAPTIVE               // Adaptive vertical coordinates
};

/**
 * @brief Turbulence closure model
 */
enum class TurbulenceClosureType {
    CONSTANT,              // Constant eddy viscosity
    SMAGORINSKY,           // Smagorinsky horizontal mixing
    KPP,                   // K-Profile Parameterization
    MELLOR_YAMADA_25,      // Mellor-Yamada level 2.5
    K_EPSILON,             // k-epsilon model
    GLS,                   // Generic Length Scale
    LES                    // Large Eddy Simulation
};

/**
 * @brief Advection scheme
 */
enum class AdvectionScheme {
    UPWIND_FIRST,          // First-order upwind
    CENTERED,              // Second-order centered
    QUICK,                 // QUICK scheme
    MUSCL,                 // MUSCL with limiters
    WENO,                  // Weighted ENO
    PPM                    // Piecewise Parabolic Method
};

/**
 * @brief Time stepping scheme
 */
enum class TimeSteppingScheme {
    FORWARD_EULER,         // Explicit Euler (for testing)
    LEAPFROG,              // Leapfrog with Robert filter
    ADAMS_BASHFORTH_2,     // 2nd order Adams-Bashforth
    ADAMS_BASHFORTH_3,     // 3rd order Adams-Bashforth
    RK3,                   // 3rd order Runge-Kutta
    RK4,                   // 4th order Runge-Kutta
    SPLIT_EXPLICIT         // Split barotropic/baroclinic
};

/**
 * @brief Bottom friction formulation
 */
enum class BottomFrictionType {
    LINEAR,                // Linear friction τ = r·u
    QUADRATIC,             // Quadratic friction τ = Cd·|u|·u
    MANNING,               // Manning formulation
    CHEZY,                 // Chézy formulation
    LOG_LAYER,             // Logarithmic boundary layer
    WAVE_CURRENT           // Combined wave-current friction
};

/**
 * @brief Tidal constituent type
 */
enum class TidalConstituent {
    M2,                    // Principal lunar semidiurnal
    S2,                    // Principal solar semidiurnal
    N2,                    // Larger lunar elliptic
    K2,                    // Lunisolar semidiurnal
    K1,                    // Lunar diurnal
    O1,                    // Principal lunar diurnal
    P1,                    // Principal solar diurnal
    Q1,                    // Larger lunar elliptic diurnal
    M4,                    // Shallow water overtide
    MS4,                   // Shallow water compound
    MN4,                   // Shallow water compound
    M6                     // Shallow water overtide
};

// =============================================================================
// Data Structures
// =============================================================================

/**
 * @brief Wind forcing data
 */
struct WindForcing {
    // Wind components at 10m height
    std::vector<double> u10;           // East component (m/s)
    std::vector<double> v10;           // North component (m/s)
    
    // Time series for time-varying wind
    std::vector<double> time;          // Time points (s)
    std::vector<std::vector<double>> u10_series;
    std::vector<std::vector<double>> v10_series;
    
    // Wind stress formulation
    bool use_bulk_formula = true;
    double drag_coefficient = 0.0013;   // Default Cd
    double air_density = 1.225;         // kg/m³
    
    // Wind stress at current time
    void computeWindStress(double t, const std::vector<double>& x,
                          const std::vector<double>& y,
                          std::vector<double>& tau_x,
                          std::vector<double>& tau_y) const;
    
    // Large and Pond (1981) drag coefficient
    static double dragCoefficientLargePond(double wind_speed);
    
    // COARE 3.0 bulk flux algorithm
    static double dragCoefficientCOARE(double wind_speed, double sst, 
                                       double air_temp, double humidity);
};

/**
 * @brief Atmospheric pressure forcing
 */
struct AtmosphericPressureForcing {
    // Sea level pressure field (Pa)
    std::vector<double> pressure;
    
    // Time series
    std::vector<double> time;
    std::vector<std::vector<double>> pressure_series;
    
    // Reference pressure
    double reference_pressure = 101325.0;  // Pa
    
    // Inverse barometer effect
    bool enable_inverse_barometer = true;
    double inverse_barometer_factor = 1.0;  // cm/hPa
    
    // Compute pressure-driven sea level change
    double computeInverseBarometer(double p, double rho_water = 1025.0) const;
};

/**
 * @brief Tidal forcing specification
 */
struct TidalForcing {
    struct Constituent {
        TidalConstituent type;
        double amplitude;              // meters
        double phase;                  // degrees
        double frequency;              // rad/s
        double nodal_factor;           // f
        double equilibrium_argument;   // u+V
    };
    
    std::vector<Constituent> constituents;
    
    // Boundary or body tide
    enum ForcingType { BOUNDARY, BODY, BOTH };
    ForcingType forcing_type = BOUNDARY;
    
    // Get tidal elevation at time t
    double getTidalElevation(double t, double lon = 0.0, double lat = 0.0) const;
    
    // Get tidal current at time t
    void getTidalCurrent(double t, double lon, double lat,
                        double& u, double& v) const;
    
    // Load tidal constituents from data file
    bool loadFromFile(const std::string& filename);
    
    // Set up common constituent sets
    void setupM2Only();
    void setupPrimaryConstituents();    // M2, S2, N2, K1, O1
    void setupFullConstituents();       // All major constituents
};

/**
 * @brief River discharge input
 */
struct RiverDischarge {
    std::string name;
    
    // Location
    double longitude;
    double latitude;
    int cell_i, cell_j;              // Grid cell indices
    
    // Discharge specification
    double discharge;                 // m³/s
    double temperature;               // °C
    double salinity;                  // PSU
    
    // Time series
    std::vector<double> time;
    std::vector<double> discharge_series;
    std::vector<double> temperature_series;
    std::vector<double> salinity_series;
    
    // Get discharge at time t
    double getDischarge(double t) const;
};

/**
 * @brief Open boundary condition specification
 */
struct OpenBoundaryCondition {
    enum Type {
        CLAMPED,                      // Fixed elevation/velocity
        RADIATION,                    // Orlanski radiation
        FLATHER,                      // Flather (radiation + external data)
        CHAPMAN,                      // Chapman free surface
        SPONGE,                       // Sponge layer relaxation
        NESTED                        // Nested model data
    };
    
    enum Location { WEST, EAST, SOUTH, NORTH };
    Location location;
    Type type;
    
    // External data for clamped/Flather boundaries
    std::vector<double> elevation;    // Sea level (m)
    std::vector<double> u_normal;     // Normal velocity (m/s)
    std::vector<double> u_tangent;    // Tangential velocity (m/s)
    
    // Time series
    std::vector<double> time;
    
    // Sponge layer parameters
    double sponge_width = 10000.0;    // meters
    double sponge_strength = 0.001;   // 1/s
    
    // Add tidal components
    TidalForcing tidal;
};

/**
 * @brief Hydrodynamic model configuration
 */
struct HydrodynamicConfig {
    // Model type
    HydrodynamicModelType model_type = HydrodynamicModelType::SHALLOW_WATER_2D;
    
    // Grid
    int nx, ny, nz;
    double dx, dy;
    double lon_min, lon_max;
    double lat_min, lat_max;
    
    // Vertical configuration (for 3D)
    VerticalCoordinateType vertical_type = VerticalCoordinateType::SIGMA;
    int num_sigma_levels = 20;
    std::vector<double> sigma_levels;
    double theta_s = 5.0;             // Surface stretching parameter
    double theta_b = 0.4;             // Bottom stretching parameter
    double h_c = 10.0;                // Critical depth for stretching
    
    // Time stepping
    TimeSteppingScheme time_scheme = TimeSteppingScheme::SPLIT_EXPLICIT;
    double dt;                        // Baroclinic time step (s)
    double dt_barotropic;             // Barotropic time step (s)
    int barotropic_steps = 30;        // Barotropic steps per baroclinic
    double end_time;
    double cfl_number = 0.8;
    bool adaptive_timestep = true;
    
    // Physical parameters
    double g = 9.81;                  // Gravity (m/s²)
    double rho_0 = 1025.0;            // Reference density (kg/m³)
    double f0 = 1e-4;                 // Coriolis parameter (rad/s)
    double beta = 1.6e-11;            // Beta parameter (1/m/s)
    bool use_coriolis = true;
    bool use_beta_plane = false;
    double latitude_reference = 45.0;
    
    // Advection
    AdvectionScheme advection_scheme = AdvectionScheme::MUSCL;
    
    // Turbulence
    TurbulenceClosureType turbulence_closure = TurbulenceClosureType::SMAGORINSKY;
    double horizontal_viscosity = 10.0;    // m²/s
    double horizontal_diffusivity = 1.0;   // m²/s
    double vertical_viscosity = 1e-4;      // m²/s
    double vertical_diffusivity = 1e-5;    // m²/s
    double smagorinsky_coefficient = 0.1;
    
    // Bottom friction
    BottomFrictionType friction_type = BottomFrictionType::QUADRATIC;
    double bottom_drag_coefficient = 0.0025;
    double manning_n = 0.025;
    double roughness_length = 0.01;   // m (for log-layer)
    
    // Wetting/drying
    bool enable_wetting_drying = true;
    double dry_depth = 0.1;           // meters
    double min_depth = 0.5;           // meters
    
    // Baroclinic (3D)
    bool use_baroclinic = true;
    bool use_nonlinear_eos = true;
    
    // Output
    int output_interval = 3600;       // seconds
    bool output_mean_fields = true;
    bool output_hourly = true;
    bool output_daily = true;
    std::string output_format = "NETCDF";
    std::string output_dir = "output";
    
    HydrodynamicConfig() = default;
};

// =============================================================================
// Main Hydrodynamic Model Class
// =============================================================================

/**
 * @brief General hydrodynamic model for ocean and coastal simulations
 *
 * Supports 2D (depth-averaged) and 3D (primitive equations) modes.
 * Can be coupled with wave models, sediment transport, and ecosystem models.
 */
class HydrodynamicModel {
public:
    HydrodynamicModel();
    ~HydrodynamicModel();
    
    // =========================================================================
    // Initialization
    // =========================================================================
    
    /**
     * @brief Initialize model with configuration
     */
    void initialize(const HydrodynamicConfig& config);
    
    /**
     * @brief Load bathymetry
     */
    void setBathymetry(const BathymetryGrid& bathymetry);
    void setBathymetry(const std::vector<double>& depth);
    
    /**
     * @brief Set initial conditions
     */
    void setInitialElevation(const std::vector<double>& eta);
    void setInitialVelocity(const std::vector<double>& u,
                           const std::vector<double>& v);
    void setInitialTemperature(const std::vector<double>& T);
    void setInitialSalinity(const std::vector<double>& S);
    void setInitialFromFile(const std::string& filename);
    
    /**
     * @brief Set up boundary conditions
     */
    void addOpenBoundary(const OpenBoundaryCondition& bc);
    void setTidalForcing(const TidalForcing& tidal);
    void addRiverInput(const RiverDischarge& river);
    
    /**
     * @brief Set forcing fields
     */
    void setWindForcing(const WindForcing& wind);
    void setAtmosphericPressure(const AtmosphericPressureForcing& pressure);
    
    // =========================================================================
    // Time Integration
    // =========================================================================
    
    /**
     * @brief Advance model by one time step
     * @return Actual time step taken
     */
    double step();
    
    /**
     * @brief Run simulation to specified end time
     */
    void run(double end_time);
    
    /**
     * @brief Run with callback at each output time
     */
    void run(double end_time, 
             std::function<void(double, const HydrodynamicModel&)> callback);
    
    // =========================================================================
    // Barotropic Mode (Fast, 2D)
    // =========================================================================
    
    /**
     * @brief Compute barotropic (depth-averaged) mode
     */
    void stepBarotropic();
    
    /**
     * @brief Compute barotropic pressure gradient
     */
    void computeBarotropicPressureGradient();
    
    /**
     * @brief Apply boundary conditions for barotropic mode
     */
    void applyBarotropicBoundaryConditions();
    
    // =========================================================================
    // Baroclinic Mode (Slow, 3D)
    // =========================================================================
    
    /**
     * @brief Compute baroclinic (3D) mode
     */
    void stepBaroclinic();
    
    /**
     * @brief Compute baroclinic pressure gradient
     */
    void computeBaroclinicPressureGradient();
    
    /**
     * @brief Solve vertical diffusion implicitly
     */
    void solveVerticalDiffusion();
    
    // =========================================================================
    // Physics Components
    // =========================================================================
    
    /**
     * @brief Compute advection terms
     */
    void computeAdvection();
    
    /**
     * @brief Compute Coriolis acceleration
     */
    void computeCoriolis();
    
    /**
     * @brief Compute horizontal mixing
     */
    void computeHorizontalMixing();
    
    /**
     * @brief Compute vertical mixing
     */
    void computeVerticalMixing();
    
    /**
     * @brief Compute bottom friction
     */
    void computeBottomFriction();
    
    /**
     * @brief Compute surface forcing (wind stress)
     */
    void computeSurfaceForcing();
    
    /**
     * @brief Apply wetting/drying
     */
    void applyWettingDrying();
    
    // =========================================================================
    // Equation of State
    // =========================================================================
    
    /**
     * @brief Compute density from temperature and salinity
     * @param T Temperature (°C)
     * @param S Salinity (PSU)
     * @param p Pressure (dbar)
     * @return Density (kg/m³)
     */
    double computeDensity(double T, double S, double p = 0.0) const;
    
    /**
     * @brief UNESCO equation of state
     */
    double densityUNESCO(double T, double S, double p) const;
    
    /**
     * @brief TEOS-10 equation of state (simplified)
     */
    double densityTEOS10(double T, double S, double p) const;
    
    /**
     * @brief Linear equation of state
     */
    double densityLinear(double T, double S) const;
    
    // =========================================================================
    // Solution Access
    // =========================================================================
    
    /**
     * @brief Get sea surface elevation
     */
    const std::vector<double>& getElevation() const { return eta; }
    double getElevation(double lon, double lat) const;
    
    /**
     * @brief Get depth-averaged velocity
     */
    void getDepthAveragedVelocity(std::vector<double>& u_bar,
                                  std::vector<double>& v_bar) const;
    double getDepthAveragedU(double lon, double lat) const;
    double getDepthAveragedV(double lon, double lat) const;
    
    /**
     * @brief Get 3D velocity field (for 3D models)
     */
    const std::vector<double>& getVelocityU() const { return u; }
    const std::vector<double>& getVelocityV() const { return v; }
    const std::vector<double>& getVelocityW() const { return w; }
    
    /**
     * @brief Get tracer fields
     */
    const std::vector<double>& getTemperature() const { return temperature; }
    const std::vector<double>& getSalinity() const { return salinity; }
    const std::vector<double>& getDensity() const { return density; }
    
    /**
     * @brief Get current simulation time
     */
    double getCurrentTime() const { return current_time; }
    
    /**
     * @brief Get total water depth (bathymetry + elevation)
     */
    void getTotalDepth(std::vector<double>& H) const;
    
    // =========================================================================
    // Derived Quantities
    // =========================================================================
    
    /**
     * @brief Compute transport streamfunction
     */
    void computeStreamfunction(std::vector<double>& psi) const;
    
    /**
     * @brief Compute vorticity
     */
    void computeVorticity(std::vector<double>& vort) const;
    
    /**
     * @brief Compute kinetic energy
     */
    double computeKineticEnergy() const;
    
    /**
     * @brief Compute volume transport
     */
    double computeVolumeTransport(int i_start, int j_start,
                                  int i_end, int j_end) const;
    
    /**
     * @brief Compute heat transport
     */
    double computeHeatTransport(int i_start, int j_start,
                               int i_end, int j_end) const;
    
    /**
     * @brief Compute mixed layer depth
     */
    void computeMixedLayerDepth(std::vector<double>& mld,
                                double delta_T = 0.5) const;
    
    // =========================================================================
    // I/O and Diagnostics
    // =========================================================================
    
    /**
     * @brief Write output at current time
     */
    void writeOutput(const std::string& filename) const;
    
    /**
     * @brief Write restart file
     */
    void writeRestart(const std::string& filename) const;
    
    /**
     * @brief Read restart file
     */
    void readRestart(const std::string& filename);
    
    /**
     * @brief Compute and print diagnostics
     */
    void printDiagnostics() const;
    
    // =========================================================================
    // Model Coupling
    // =========================================================================
    
    /**
     * @brief Couple with wave model for wave-current interaction
     */
    void setWaveModel(std::shared_ptr<WaveModel> wave);
    
    /**
     * @brief Get wave-enhanced bottom stress
     */
    void getWaveEnhancedFriction(std::vector<double>& tau_b) const;
    
    /**
     * @brief Provide currents to wave model
     */
    void provideCurrentsToWaveModel();
    
private:
    // Configuration
    HydrodynamicConfig config;
    
    // Grid
    int nx, ny, nz;
    double dx, dy;
    double lon_min, lat_min;
    std::vector<double> sigma;        // Sigma levels (0 at surface, -1 at bottom)
    
    // Bathymetry
    std::vector<double> h;            // Still water depth (positive)
    std::vector<double> H;            // Total depth (h + eta)
    std::vector<int> mask;            // Land/sea mask (0=land, 1=sea)
    
    // State variables - 2D (barotropic)
    std::vector<double> eta;          // Sea surface elevation
    std::vector<double> u_bar;        // Depth-averaged u velocity
    std::vector<double> v_bar;        // Depth-averaged v velocity
    
    // State variables - 3D (baroclinic)
    std::vector<double> u;            // u velocity (nx × ny × nz)
    std::vector<double> v;            // v velocity
    std::vector<double> w;            // Vertical velocity (omega in sigma)
    std::vector<double> temperature;  // Temperature
    std::vector<double> salinity;     // Salinity
    std::vector<double> density;      // In-situ density
    
    // Tendencies
    std::vector<double> du_dt, dv_dt, deta_dt;
    std::vector<double> dT_dt, dS_dt;
    
    // Pressure gradient
    std::vector<double> dP_dx, dP_dy; // Baroclinic pressure gradient
    
    // Mixing coefficients
    std::vector<double> Kv;           // Vertical diffusivity
    std::vector<double> Av;           // Vertical viscosity
    std::vector<double> Kh;           // Horizontal diffusivity
    std::vector<double> Ah;           // Horizontal viscosity
    
    // Forcing
    WindForcing wind_forcing;
    AtmosphericPressureForcing pressure_forcing;
    TidalForcing tidal_forcing;
    std::vector<RiverDischarge> rivers;
    std::vector<OpenBoundaryCondition> open_boundaries;
    
    // Current forcing fields
    std::vector<double> tau_x, tau_y; // Wind stress
    std::vector<double> p_atm;        // Atmospheric pressure
    
    // Coupled wave model
    std::shared_ptr<WaveModel> wave_model;
    
    // Time stepping
    double current_time;
    double dt, dt_baro;
    int step_count;
    
    // Leapfrog storage
    std::vector<double> eta_old, u_bar_old, v_bar_old;
    std::vector<double> u_old, v_old;
    
    // Index helpers
    int idx2d(int i, int j) const { return i + j * nx; }
    int idx3d(int i, int j, int k) const { return i + j * nx + k * nx * ny; }
    
    // Grid metrics
    double f_coriolis(int j) const;   // Coriolis at cell j
    double cell_area(int i, int j) const;
    
    // Sigma coordinate functions
    double z_at_sigma(int i, int j, int k) const;
    double dz_at_sigma(int i, int j, int k) const;
    void computeSigmaCoordinates();
    
    // Numerical methods
    void computeFlux(const std::vector<double>& q, const std::vector<double>& u,
                    std::vector<double>& flux, char direction);
    double minmod(double a, double b) const;
    double muscl_reconstruct(double qL, double qC, double qR, char side) const;
    
    // Turbulence closure
    void computeTurbulentMixing();
    void computeKPP();
    void computeMellorYamada();
    void computeSmagorinsky();
};

// =============================================================================
// Coastal Hydrodynamics Model
// =============================================================================

/**
 * @brief Specialized coastal hydrodynamics model
 *
 * Optimized for nearshore and coastal applications including:
 * - Storm surge
 * - Tidal flats
 * - Estuaries
 * - Coastal flooding
 */
class CoastalHydrodynamicsModel : public HydrodynamicModel {
public:
    CoastalHydrodynamicsModel();
    
    // Storm surge specific methods
    void initializeStormSurge(const HydrodynamicConfig& config);
    void setHurricaneTrack(const std::vector<double>& time,
                          const std::vector<double>& lon,
                          const std::vector<double>& lat,
                          const std::vector<double>& pressure,
                          const std::vector<double>& max_wind,
                          const std::vector<double>& rmw);
    
    /**
     * @brief Compute Holland (1980) hurricane wind profile
     */
    void computeHollandWinds(double t, std::vector<double>& u10,
                            std::vector<double>& v10);
    
    /**
     * @brief Compute parametric hurricane pressure field
     */
    void computeHollandPressure(double t, std::vector<double>& p);
    
    // Inundation mapping
    void computeMaxInundation(std::vector<double>& max_eta,
                             std::vector<double>& max_depth,
                             std::vector<double>& max_velocity);
    
    // High-water marks
    void extractHighWaterMarks(std::vector<double>& hwm_lon,
                              std::vector<double>& hwm_lat,
                              std::vector<double>& hwm_elevation);
    
private:
    // Hurricane track data
    std::vector<double> track_time;
    std::vector<double> track_lon, track_lat;
    std::vector<double> central_pressure;
    std::vector<double> max_wind_speed;
    std::vector<double> radius_max_wind;
    
    // Holland parameters
    double holland_B = 1.5;
    double ambient_pressure = 101325.0;
    
    // Maximum value tracking
    std::vector<double> max_eta;
    std::vector<double> max_depth;
    std::vector<double> max_speed;
};

// =============================================================================
// Estuarine Model
// =============================================================================

/**
 * @brief Estuarine dynamics model
 *
 * Models river-ocean interaction including:
 * - Salt intrusion
 * - Estuarine circulation (gravitational)
 * - Tidal pumping
 * - Stratification and mixing
 */
class EstuarineModel : public HydrodynamicModel {
public:
    EstuarineModel();
    
    /**
     * @brief Initialize for estuary modeling
     */
    void initializeEstuary(const HydrodynamicConfig& config,
                          double river_width,
                          double river_depth,
                          double ocean_salinity);
    
    /**
     * @brief Set up river head boundary
     */
    void setRiverHeadBoundary(double discharge, double fresh_temperature);
    
    /**
     * @brief Set up ocean mouth boundary with tides
     */
    void setOceanBoundary(const TidalForcing& tides,
                         double ocean_salinity,
                         double ocean_temperature);
    
    /**
     * @brief Compute salt intrusion length
     */
    double computeSaltIntrusionLength(double salinity_threshold = 2.0);
    
    /**
     * @brief Compute estuarine Richardson number
     */
    double computeEstuarineRichardson();
    
    /**
     * @brief Classify estuary type
     */
    enum EstuaryType { SALT_WEDGE, PARTIALLY_MIXED, WELL_MIXED };
    EstuaryType classifyEstuary();
    
    /**
     * @brief Compute tidal excursion
     */
    double computeTidalExcursion(double x_station);
    
    /**
     * @brief Compute flushing time
     */
    double computeFlushingTime();
    
private:
    double river_discharge;
    double ocean_salinity;
    double ocean_temperature;
};

// =============================================================================
// Helper Utilities
// =============================================================================

namespace HydrodynamicUtils {
    /**
     * @brief Compute tidal harmonic analysis
     * @param time Time series
     * @param eta Elevation time series
     * @param constituents List of constituents to extract
     * @return Map of constituent to (amplitude, phase)
     */
    std::map<TidalConstituent, std::pair<double, double>>
    tidalHarmonicAnalysis(const std::vector<double>& time,
                         const std::vector<double>& eta,
                         const std::vector<TidalConstituent>& constituents);
    
    /**
     * @brief Get tidal constituent properties
     */
    void getConstituentProperties(TidalConstituent type,
                                 double& period_hours,
                                 double& frequency);
    
    /**
     * @brief Compute astronomical argument
     */
    double computeAstronomicalArgument(TidalConstituent type, double julian_day);
    
    /**
     * @brief Compute nodal corrections
     */
    void computeNodalCorrections(TidalConstituent type, double year,
                                double& f, double& u);
    
    /**
     * @brief Manning to Cd conversion
     */
    double manningToDragCoefficient(double n, double h);
    
    /**
     * @brief Cd to Manning conversion
     */
    double dragCoefficientToManning(double Cd, double h);
    
    /**
     * @brief Compute wave-current bottom stress
     */
    void computeWaveCurrentStress(double u, double v, double H,
                                 double Hs, double Tp, double wave_dir,
                                 double& tau_x, double& tau_y);
    
    /**
     * @brief Load bathymetry from NetCDF
     */
    bool loadBathymetryNetCDF(const std::string& filename,
                             double lon_min, double lon_max,
                             double lat_min, double lat_max,
                             BathymetryGrid& grid);
    
    /**
     * @brief Generate idealized estuary bathymetry
     */
    BathymetryGrid generateEstuaryBathymetry(double length, double width,
                                            double max_depth,
                                            double river_depth);
}

// =============================================================================
// Configuration Parsing
// =============================================================================

/**
 * @brief Parse hydrodynamic configuration from file
 */
struct HydrodynamicConfigReader {
    static bool parseConfig(const std::string& filename,
                           HydrodynamicConfig& config);
    
    static bool parseWindForcing(const std::map<std::string, std::string>& section,
                                WindForcing& wind);
    
    static bool parseTidalForcing(const std::map<std::string, std::string>& section,
                                 TidalForcing& tidal);
    
    static bool parseRiverInputs(const std::map<std::string, std::string>& section,
                                std::vector<RiverDischarge>& rivers);
    
    static bool parseBoundaryConditions(
        const std::map<std::string, std::string>& section,
        std::vector<OpenBoundaryCondition>& boundaries);
};

} // namespace FSRM

#endif // HYDRODYNAMIC_MODEL_HPP
