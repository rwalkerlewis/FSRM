/**
 * @file TsunamiModel.hpp
 * @brief Tsunami modeling with oceanic coupling for earthquake-generated waves
 *
 * Implements tsunami generation, propagation, and coastal inundation:
 * - Nonlinear Shallow Water Equations (NSWE) solver
 * - Coseismic seafloor displacement using Okada (1985) dislocation model
 * - Coupling with solid earth earthquake rupture models
 * - Bathymetry and topography handling (GEBCO, ETOPO, Shuttle Radar)
 * - Runup and inundation calculations
 * - Tide gauge and DART buoy synthetic recordings
 * - Multi-grid nesting for coastal refinement
 *
 * Supports realistic fault geometries:
 * - Cascadia Subduction Zone (CSZ)
 * - Alaska-Aleutian Subduction Zone
 * - Japan Trench
 * - Chile-Peru Trench
 * - Sumatra-Andaman
 *
 * @author FSRM Development Team
 */

#ifndef TSUNAMI_MODEL_HPP
#define TSUNAMI_MODEL_HPP

#include "ReservoirSim.hpp"
#include "SeismicSource.hpp"
#include "FaultModel.hpp"
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
class DiscontinuousGalerkin;

/**
 * @brief Numerical scheme for shallow water equations
 */
enum class SWESolver {
    GODUNOV,           // Godunov scheme with exact Riemann solver
    ROE,               // Roe approximate Riemann solver
    HLLC,              // Harten-Lax-van Leer-Contact solver
    CENTRAL,           // Central differences (for linear)
    DG_RKDG,           // Discontinuous Galerkin with Runge-Kutta
    FWAVE              // F-wave method (GeoClaw style)
};

/**
 * @brief Wetting/drying treatment for inundation
 */
enum class WettingDryingMethod {
    NONE,              // No wetting/drying (deep water only)
    THIN_LAYER,        // Thin layer cutoff
    NEGATIVE_DEPTH,    // Allow negative depth with correction
    MOMENTUM_CONSERVING // Momentum-conserving approach
};

/**
 * @brief Friction model for bottom stress
 */
enum class BottomFrictionModel {
    NONE,              // No friction (frictionless)
    CHEZY,             // Chézy coefficient
    MANNING,           // Manning's n roughness
    QUADRATIC,         // Quadratic drag law
    VARIABLE_MANNING   // Spatially variable Manning's n
};

/**
 * @brief Bathymetry/topography data source
 */
enum class BathymetrySource {
    CONSTANT,          // Constant depth (for testing)
    ANALYTIC,          // Analytic profile (conical island, etc.)
    GEBCO,             // GEBCO 2023 global bathymetry
    ETOPO,             // ETOPO1/ETOPO2 global relief
    SRTM,              // Shuttle Radar Topography Mission
    LIDAR,             // High-resolution coastal LiDAR
    USER_FILE          // User-provided NetCDF/ASCII file
};

/**
 * @brief Subfault representation for finite fault tsunami source
 */
struct TsunamiSubfault {
    // Location (center of subfault, geographic coordinates)
    double longitude;      // degrees East
    double latitude;       // degrees North
    double depth;          // km below sea level (to top of fault)
    
    // Dimensions
    double length;         // km (along strike)
    double width;          // km (down dip)
    
    // Orientation
    double strike;         // degrees clockwise from North
    double dip;            // degrees from horizontal
    double rake;           // degrees (slip direction)
    
    // Slip
    double slip;           // meters
    double rise_time;      // seconds
    double rupture_time;   // seconds (delay from rupture initiation)
    
    // Optional: variable slip
    double slip_along_strike;  // meters/km gradient
    double slip_down_dip;      // meters/km gradient
    
    TsunamiSubfault() :
        longitude(-125.0), latitude(45.0), depth(10.0),
        length(50.0), width(50.0),
        strike(0.0), dip(15.0), rake(90.0),
        slip(5.0), rise_time(30.0), rupture_time(0.0),
        slip_along_strike(0.0), slip_down_dip(0.0) {}
};

/**
 * @brief Cascadia Subduction Zone fault model parameters
 * 
 * Based on published geodetic and geological studies:
 * - Flück et al. (1997) - Geometry
 * - Wang et al. (2003) - Coupling
 * - Witter et al. (2013) - Paleoseismic constraints
 * - USGS SLAB2.0 model
 */
struct CascadiaFaultModel {
    // Fault extent
    double north_latitude;     // Northern terminus (e.g., 50.0° - Vancouver Island)
    double south_latitude;     // Southern terminus (e.g., 40.0° - Cape Mendocino)
    double trench_longitude;   // Approximate trench position (-126° to -124°)
    
    // Geometry
    double average_strike;     // degrees (approximately N-S, ~350°)
    double shallow_dip;        // degrees (5-10° at trench)
    double deep_dip;           // degrees (15-20° at depth)
    double transition_depth;   // km (where dip steepens)
    
    // Rupture parameters
    double rupture_width;      // km (down-dip extent, typically 80-150 km)
    double average_slip;       // meters (15-20 m for M9.0)
    double peak_slip;          // meters (up to 40 m near trench)
    
    // Segmentation
    int num_along_strike;      // Number of subfaults along strike
    int num_down_dip;          // Number of subfaults down dip
    
    // Rupture dynamics
    double rupture_velocity;   // km/s (typically 2.5-3.0)
    double hypocenter_lat;     // degrees (nucleation point)
    double hypocenter_lon;     // degrees
    double hypocenter_depth;   // km
    
    CascadiaFaultModel() :
        north_latitude(50.0), south_latitude(40.0), trench_longitude(-125.0),
        average_strike(350.0), shallow_dip(8.0), deep_dip(18.0), transition_depth(20.0),
        rupture_width(120.0), average_slip(17.0), peak_slip(35.0),
        num_along_strike(50), num_down_dip(10),
        rupture_velocity(2.8), hypocenter_lat(45.0), hypocenter_lon(-124.5), hypocenter_depth(25.0) {}
        
    // Generate subfault array from model parameters
    std::vector<TsunamiSubfault> generateSubfaults() const;
    
    // Get total seismic moment
    double getTotalMoment() const;
    
    // Get moment magnitude
    double getMagnitude() const;
};

/**
 * @brief Bathymetry/topography grid data
 */
struct BathymetryGrid {
    // Grid parameters
    double lon_min, lon_max;    // Longitude range (degrees East)
    double lat_min, lat_max;    // Latitude range (degrees North)
    double dlon, dlat;          // Grid spacing (degrees)
    int nlon, nlat;             // Grid dimensions
    
    // Depth data (positive downward, negative for land)
    std::vector<double> depth;  // nlon × nlat array
    
    // Optional: Manning roughness at each point
    std::vector<double> manning_n;
    
    // Coordinate reference
    std::string projection;     // e.g., "EPSG:4326" for WGS84
    
    BathymetryGrid() :
        lon_min(-130.0), lon_max(-120.0), lat_min(40.0), lat_max(50.0),
        dlon(0.01), dlat(0.01), nlon(1000), nlat(1000) {}
    
    // Load from file
    bool loadFromNetCDF(const std::string& filename);
    bool loadFromASCII(const std::string& filename);
    bool loadFromGeoTIFF(const std::string& filename);
    
    // Generate analytic bathymetry for testing
    void generateConicalIsland(double lon_center, double lat_center,
                              double radius, double island_height,
                              double surrounding_depth);
    void generateContinentalShelf(double shelf_width, double shelf_depth,
                                 double slope_angle, double deep_depth);
    
    // Interpolate depth at a point
    double getDepth(double lon, double lat) const;
    double getManningN(double lon, double lat) const;
    
    // Gradient for slope
    void getGradient(double lon, double lat, double& dzdlon, double& dzdlat) const;
};

/**
 * @brief Tide gauge / DART buoy recording station
 */
struct TsunamiGauge {
    std::string name;
    std::string id;             // Station ID (e.g., NOAA/NDBC ID)
    
    double longitude;
    double latitude;
    double depth;               // Deployment depth (0 for tide gauge, > 0 for DART)
    
    // Record type
    bool is_dart_buoy;          // DART = true, tide gauge = false
    
    // Recorded time series
    std::vector<double> time;
    std::vector<double> eta;    // Sea surface elevation
    std::vector<double> u;      // East velocity (optional)
    std::vector<double> v;      // North velocity (optional)
    
    TsunamiGauge() :
        name(""), id(""), longitude(0), latitude(0), depth(0),
        is_dart_buoy(false) {}
    
    // Output
    void writeASCII(const std::string& filename) const;
    void writeSAC(const std::string& filename) const;
    
    // Get maximum amplitude
    double getMaxAmplitude() const;
    double getFirstArrivalTime(double threshold = 0.01) const;
};

/**
 * @brief Inundation metrics for a coastal location
 */
struct InundationMetrics {
    double max_runup;           // Maximum runup elevation (m above MSL)
    double max_flow_depth;      // Maximum flow depth (m)
    double max_velocity;        // Maximum current velocity (m/s)
    double max_momentum_flux;   // Maximum momentum flux ρhu² (kg/m/s²)
    double inundation_distance; // Horizontal inundation distance (m)
    double arrival_time;        // First arrival time (s)
    double duration;            // Duration of inundation (s)
    
    // Location
    double longitude;
    double latitude;
};

/**
 * @brief Configuration for tsunami simulation
 */
struct TsunamiConfig {
    // Domain
    double lon_min, lon_max;
    double lat_min, lat_max;
    
    // Grid resolution
    double base_resolution;     // degrees (ocean, typically 0.01-0.1)
    double coastal_resolution;  // degrees (near coast, 0.001-0.01)
    
    // Time stepping
    double dt;                  // Time step (s)
    double end_time;            // Simulation duration (s)
    double cfl_number;          // CFL number for adaptive dt
    bool adaptive_timestep;
    
    // Numerical scheme
    SWESolver solver;
    int spatial_order;          // 1, 2, or higher for DG
    bool use_minmod_limiter;
    bool use_positivity_preserving;
    
    // Wetting/drying
    WettingDryingMethod wetting_drying;
    double dry_tolerance;       // Minimum depth to consider "wet" (m)
    
    // Friction
    BottomFrictionModel friction_model;
    double manning_n_ocean;     // Manning's n for ocean (0.025)
    double manning_n_land;      // Manning's n for land (0.03-0.1)
    
    // Dispersion (Boussinesq corrections)
    bool use_dispersion;
    double dispersion_threshold_depth; // Only apply if h > threshold
    
    // Source
    bool use_kinematic_source;  // Time-dependent seafloor motion
    double source_rise_time;    // Rise time for instantaneous source
    
    // Boundary conditions
    std::string bc_ocean;       // "OPEN", "REFLECTIVE", "SPONGE"
    double sponge_width;        // Width of sponge layer (degrees)
    
    // Nested grids
    bool use_nested_grids;
    int num_nesting_levels;
    double nesting_ratio;       // Refinement ratio (typically 3)
    
    // Output
    int output_interval;        // Time steps between outputs
    bool save_max_values;       // Track maximum eta, velocity
    bool save_inundation_map;
    std::string output_format;  // "NETCDF", "HDF5", "ASCII"
    
    TsunamiConfig() :
        lon_min(-130.0), lon_max(-120.0), lat_min(40.0), lat_max(50.0),
        base_resolution(1.0/60.0), coastal_resolution(1.0/600.0),
        dt(0.5), end_time(14400.0), cfl_number(0.8), adaptive_timestep(true),
        solver(SWESolver::FWAVE), spatial_order(2),
        use_minmod_limiter(true), use_positivity_preserving(true),
        wetting_drying(WettingDryingMethod::MOMENTUM_CONSERVING),
        dry_tolerance(0.001),
        friction_model(BottomFrictionModel::MANNING),
        manning_n_ocean(0.025), manning_n_land(0.03),
        use_dispersion(false), dispersion_threshold_depth(50.0),
        use_kinematic_source(true), source_rise_time(60.0),
        bc_ocean("OPEN"), sponge_width(0.5),
        use_nested_grids(false), num_nesting_levels(0), nesting_ratio(3.0),
        output_interval(60), save_max_values(true), save_inundation_map(true),
        output_format("NETCDF") {}
};

/**
 * @brief Okada (1985) dislocation model for seafloor deformation
 * 
 * Computes surface displacement from rectangular fault slip in an
 * elastic half-space. Used to generate tsunami initial conditions.
 */
class OkadaModel {
public:
    OkadaModel();
    
    /**
     * @brief Compute surface displacement from a single rectangular fault
     * 
     * @param lon Longitude of observation point (degrees)
     * @param lat Latitude of observation point (degrees)
     * @param subfault Fault parameters
     * @param ux Output: East displacement (m)
     * @param uy Output: North displacement (m)
     * @param uz Output: Vertical displacement (m)
     */
    void computeDisplacement(double lon, double lat,
                            const TsunamiSubfault& subfault,
                            double& ux, double& uy, double& uz) const;
    
    /**
     * @brief Compute total displacement from multiple subfaults
     */
    void computeTotalDisplacement(double lon, double lat,
                                 const std::vector<TsunamiSubfault>& subfaults,
                                 double& ux, double& uy, double& uz) const;
    
    /**
     * @brief Compute displacement field on a grid
     * 
     * @param grid Bathymetry grid (provides coordinates)
     * @param subfaults Fault model
     * @param uz_field Output vertical displacement field
     */
    void computeDisplacementField(const BathymetryGrid& grid,
                                 const std::vector<TsunamiSubfault>& subfaults,
                                 std::vector<double>& uz_field) const;
    
    /**
     * @brief Kajiura (1963) filter for water column response
     * 
     * Accounts for filtering of short-wavelength seafloor displacement
     * by the water column. Critical for accurate initial conditions.
     * 
     * @param grid Bathymetry grid
     * @param uz_field Seafloor displacement (modified in place)
     */
    void applyKajiuraFilter(const BathymetryGrid& grid,
                           std::vector<double>& uz_field) const;
    
    // Material properties
    void setShearModulus(double mu) { shear_modulus = mu; }
    void setPoissonRatio(double nu) { poisson_ratio = nu; }
    
private:
    double shear_modulus;      // Pa (default: 30 GPa)
    double poisson_ratio;      // dimensionless (default: 0.25)
    
    // Okada (1985) Green's functions
    void okadaGreenFunction(double x, double y, double d,
                           double dip, double L, double W, double slip,
                           double& ux, double& uy, double& uz) const;
    
    // Coordinate transformation helpers
    void geoToLocal(double lon, double lat, double lon0, double lat0,
                   double strike, double& x, double& y) const;
};

/**
 * @brief Shallow Water Equations solver for tsunami propagation
 * 
 * Solves the nonlinear shallow water equations:
 *   ∂h/∂t + ∇·(hu) = 0
 *   ∂(hu)/∂t + ∇·(huu + ½gh²I) = -gh∇b - τ_b/ρ
 * 
 * where h = water depth, u = velocity, b = bathymetry, τ_b = bottom stress
 */
class ShallowWaterSolver {
public:
    ShallowWaterSolver();
    ~ShallowWaterSolver();
    
    /**
     * @brief Initialize solver with grid and bathymetry
     */
    void initialize(const BathymetryGrid& bathymetry,
                   const TsunamiConfig& config);
    
    /**
     * @brief Set initial conditions (sea surface perturbation)
     */
    void setInitialCondition(const std::vector<double>& eta_initial);
    
    /**
     * @brief Set time-dependent seafloor motion (kinematic source)
     */
    void setSeafloorMotion(std::function<void(double, std::vector<double>&)> motion_func);
    
    /**
     * @brief Advance solution by one time step
     * @return Actual time step taken (may differ if adaptive)
     */
    double step();
    
    /**
     * @brief Run simulation to specified end time
     */
    void run(double end_time);
    
    /**
     * @brief Get current solution
     */
    void getSolution(std::vector<double>& eta,
                    std::vector<double>& hu,
                    std::vector<double>& hv) const;
    
    /**
     * @brief Get sea surface elevation at a point
     */
    double getEta(double lon, double lat) const;
    
    /**
     * @brief Get velocity at a point
     */
    void getVelocity(double lon, double lat, double& u, double& v) const;
    
    /**
     * @brief Record gauge time series
     */
    void recordGauges(double t);
    
    /**
     * @brief Get maximum values tracked during simulation
     */
    void getMaximumValues(std::vector<double>& max_eta,
                         std::vector<double>& max_speed,
                         std::vector<double>& max_momentum_flux) const;
    
    /**
     * @brief Compute inundation metrics along a coastal profile
     */
    std::vector<InundationMetrics> computeInundation() const;
    
    // Access gauges
    void addGauge(const TsunamiGauge& gauge);
    const std::vector<TsunamiGauge>& getGauges() const { return gauges; }
    
    // Get current time
    double getCurrentTime() const { return current_time; }
    
private:
    // Configuration
    TsunamiConfig config;
    
    // Grid
    int nx, ny;
    double dx, dy;
    double lon_min, lat_min;
    
    // Bathymetry
    std::vector<double> depth;  // Still water depth (positive downward)
    std::vector<double> manning_n;
    
    // Conserved variables: q = [h, hu, hv]
    std::vector<double> h;      // Water depth
    std::vector<double> hu;     // x-momentum
    std::vector<double> hv;     // y-momentum
    
    // Auxiliary
    std::vector<double> eta;    // Sea surface elevation η = h - d
    
    // Time stepping
    double current_time;
    double dt;
    
    // Maximum value tracking
    std::vector<double> max_eta;
    std::vector<double> max_speed;
    std::vector<double> max_momentum_flux;
    std::vector<double> arrival_time;
    
    // Gauges
    std::vector<TsunamiGauge> gauges;
    
    // Seafloor motion function
    std::function<void(double, std::vector<double>&)> seafloor_motion;
    
    // Numerical methods
    void computeFluxes(const std::vector<double>& q,
                      std::vector<double>& fx, std::vector<double>& fy);
    
    void riemannSolverHLLC(double hL, double huL, double hvL,
                          double hR, double huR, double hvR,
                          double& f_h, double& f_hu, double& f_hv);
    
    void riemannSolverFwave(double hL, double huL, double bL,
                           double hR, double huR, double bR,
                           double& amdq, double& apdq);
    
    void applyBottomFriction(double dt);
    void applyWettingDrying();
    void applyBoundaryConditions();
    
    double computeCFL() const;
    
    // Grid indexing
    int idx(int i, int j) const { return i + j * nx; }
};

/**
 * @brief Coupled tsunami-earthquake simulation
 * 
 * Combines dynamic earthquake rupture with tsunami generation and propagation.
 * Supports both kinematic (prescribed) and dynamic (spontaneous) rupture.
 */
class CoupledTsunamiEarthquake {
public:
    CoupledTsunamiEarthquake();
    
    /**
     * @brief Initialize from fault model and bathymetry
     */
    void initialize(const std::vector<TsunamiSubfault>& fault_model,
                   const BathymetryGrid& bathymetry,
                   const TsunamiConfig& config);
    
    /**
     * @brief Initialize Cascadia M9.0 scenario
     */
    void initializeCascadia(const CascadiaFaultModel& csz_model,
                           const BathymetryGrid& bathymetry,
                           const TsunamiConfig& config);
    
    /**
     * @brief Set up coupling with solid earth model
     * 
     * Uses the SeismicFaultModel for dynamic rupture simulation.
     * Seafloor displacement is updated based on coseismic slip.
     */
    void coupleWithEarthquakeModel(std::shared_ptr<SeismicFaultModel> earthquake);
    
    /**
     * @brief Run coupled simulation
     * 
     * 1. Compute earthquake rupture (or use kinematic model)
     * 2. Compute seafloor displacement
     * 3. Initialize and propagate tsunami
     */
    void run();
    
    /**
     * @brief Get tsunami arrival times at specified locations
     */
    std::map<std::string, double> getArrivalTimes() const;
    
    /**
     * @brief Get maximum wave heights at gauge locations
     */
    std::map<std::string, double> getMaxAmplitudes() const;
    
    /**
     * @brief Get complete gauge records
     */
    const std::vector<TsunamiGauge>& getGaugeRecords() const;
    
    /**
     * @brief Get inundation map
     */
    void getInundationMap(std::vector<double>& max_eta,
                         std::vector<double>& max_flow_depth,
                         std::vector<double>& max_velocity) const;
    
    /**
     * @brief Write output files
     */
    void writeOutput(const std::string& output_dir) const;
    
private:
    // Fault model
    std::vector<TsunamiSubfault> subfaults;
    
    // Okada model for deformation
    OkadaModel okada;
    
    // Shallow water solver
    std::unique_ptr<ShallowWaterSolver> swe_solver;
    
    // Bathymetry
    BathymetryGrid bathymetry;
    
    // Configuration
    TsunamiConfig config;
    
    // Coupled earthquake model (optional)
    std::shared_ptr<SeismicFaultModel> earthquake_model;
    
    // Results
    std::vector<double> seafloor_displacement;
    
    // Helper functions
    void computeSeafloorDisplacement();
    void generateKinematicSource();
};

/**
 * @brief West Coast of North America gauge network
 * 
 * Pre-defined tide gauge and DART buoy locations for tsunami modeling
 */
class WestCoastGaugeNetwork {
public:
    /**
     * @brief Get all NOAA tide gauge stations
     */
    static std::vector<TsunamiGauge> getTideGauges();
    
    /**
     * @brief Get all DART buoy stations
     */
    static std::vector<TsunamiGauge> getDARTBuoys();
    
    /**
     * @brief Get combined network
     */
    static std::vector<TsunamiGauge> getAllStations();
    
    /**
     * @brief Get stations within a region
     */
    static std::vector<TsunamiGauge> getStationsInRegion(
        double lon_min, double lon_max,
        double lat_min, double lat_max);
    
    /**
     * @brief Get major metropolitan area stations
     */
    static std::vector<TsunamiGauge> getMetropolitanStations();
    
    // Specific locations
    static TsunamiGauge crescent_city();     // Very exposed to Cascadia
    static TsunamiGauge astoria();
    static TsunamiGauge westport();
    static TsunamiGauge seattle();
    static TsunamiGauge victoria();
    static TsunamiGauge tofino();
    static TsunamiGauge san_francisco();
    static TsunamiGauge monterey();
    static TsunamiGauge los_angeles();
    static TsunamiGauge san_diego();
    
    // DART buoys (deep ocean)
    static TsunamiGauge dart_46404();        // West of Oregon
    static TsunamiGauge dart_46407();        // West of Washington
    static TsunamiGauge dart_46419();        // North Pacific
};

/**
 * @brief Cascadia Subduction Zone pre-built scenarios
 */
namespace CascadiaScenarios {
    /**
     * @brief Full-margin M9.0+ rupture
     * Based on 1700 CE event (Wang et al., 2013)
     */
    CascadiaFaultModel fullMarginM9();
    
    /**
     * @brief Southern segment rupture (M8.0-8.5)
     * Cape Mendocino to central Oregon
     */
    CascadiaFaultModel southernSegmentM8();
    
    /**
     * @brief Northern segment rupture (M8.0-8.5)
     * Washington to Vancouver Island
     */
    CascadiaFaultModel northernSegmentM8();
    
    /**
     * @brief Central Oregon rupture (M8.0)
     * Localized high-slip patch
     */
    CascadiaFaultModel centralOregonM8();
    
    /**
     * @brief Deep rupture scenario
     * Episodic Tremor and Slip (ETS) zone to trench
     */
    CascadiaFaultModel deepRuptureM8();
    
    /**
     * @brief Worst-case scenario for planning
     * Maximum slip distribution based on paleotsunami data
     */
    CascadiaFaultModel worstCaseScenario();
}

/**
 * @brief Helper functions for tsunami modeling
 */
namespace TsunamiUtils {
    /**
     * @brief Estimate moment magnitude from subfault parameters
     */
    double computeMomentMagnitude(const std::vector<TsunamiSubfault>& subfaults,
                                 double shear_modulus = 30e9);
    
    /**
     * @brief Estimate tsunami wave speed from depth
     */
    double waveSpeed(double depth);  // sqrt(g * depth)
    
    /**
     * @brief Estimate arrival time (great circle distance / wave speed)
     */
    double estimateArrivalTime(double lon_source, double lat_source,
                              double lon_target, double lat_target,
                              double average_depth);
    
    /**
     * @brief Apply Green's law amplification
     * η₂/η₁ = (h₁/h₂)^(1/4)
     */
    double greensLaw(double eta_deep, double depth_deep, double depth_shallow);
    
    /**
     * @brief Compute inundation distance from wave height
     * Uses empirical Synolakis formula
     */
    double estimateInundation(double wave_height, double beach_slope);
    
    /**
     * @brief Great circle distance between two points
     */
    double haversineDistance(double lon1, double lat1, double lon2, double lat2);
    
    /**
     * @brief Load GEBCO bathymetry for a region
     */
    bool loadGEBCO(const std::string& filename,
                  double lon_min, double lon_max,
                  double lat_min, double lat_max,
                  BathymetryGrid& grid);
    
    /**
     * @brief Generate synthetic bathymetry for testing
     */
    BathymetryGrid generateTestBathymetry(
        double lon_min, double lon_max,
        double lat_min, double lat_max,
        double resolution,
        double shelf_depth = 200.0,
        double deep_depth = 3000.0);
}

/**
 * @brief Configuration reader extensions for tsunami
 */
struct TsunamiConfigReader {
    /**
     * @brief Parse tsunami section from config file
     */
    static bool parseTsunamiConfig(const std::map<std::string, std::string>& section,
                                  TsunamiConfig& config);
    
    /**
     * @brief Parse fault model from config file
     */
    static bool parseFaultModel(const std::map<std::string, std::string>& section,
                               std::vector<TsunamiSubfault>& subfaults);
    
    /**
     * @brief Parse Cascadia model from config file
     */
    static bool parseCascadiaModel(const std::map<std::string, std::string>& section,
                                  CascadiaFaultModel& model);
    
    /**
     * @brief Parse gauge locations from config file
     */
    static bool parseGauges(const std::map<std::string, std::string>& section,
                           std::vector<TsunamiGauge>& gauges);
};

} // namespace FSRM

#endif // TSUNAMI_MODEL_HPP
