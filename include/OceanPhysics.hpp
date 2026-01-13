/**
 * @file OceanPhysics.hpp
 * @brief Comprehensive ocean physics modeling components
 *
 * Provides specialized modules for ocean physics processes:
 * - Surface gravity wave modeling (spectral and phase-resolving)
 * - Internal wave dynamics and breaking
 * - Ocean acoustics and sound propagation
 * - Sediment transport and morphodynamics
 * - Thermohaline circulation
 * - Ocean mixing and turbulence
 * - Air-sea interaction
 *
 * Physical parameterizations include:
 * - Wind-wave generation
 * - Wave-wave interactions (quadruplet)
 * - Wave breaking and dissipation
 * - Wave-induced mixing
 * - Stokes drift
 * - Langmuir turbulence
 *
 * @author FSRM Development Team
 */

#ifndef OCEAN_PHYSICS_HPP
#define OCEAN_PHYSICS_HPP

#include "FSRM.hpp"
#include "HydrodynamicModel.hpp"
#include "TsunamiModel.hpp"
#include <vector>
#include <array>
#include <string>
#include <map>
#include <memory>
#include <functional>
#include <complex>

namespace FSRM {

// Forward declarations
class HydrodynamicModel;

// =============================================================================
// Enumerations
// =============================================================================

/**
 * @brief Wave spectrum type
 */
enum class WaveSpectrumType {
    JONSWAP,               // JONSWAP spectrum (developing seas)
    PIERSON_MOSKOWITZ,     // Pierson-Moskowitz (fully developed)
    TMA,                   // TMA spectrum (shallow water modified)
    BRETSCHNEIDER,         // Bretschneider two-parameter
    OCHI_HUBBLE,           // Ochi-Hubble six-parameter (bimodal)
    MEASURED               // User-provided measured spectrum
};

/**
 * @brief Wave breaking model
 */
enum class WaveBreakingModel {
    NONE,                  // No breaking
    WHITECAPPING,          // Deep water whitecapping
    DEPTH_LIMITED,         // Shallow water depth-limited breaking
    BATTJES_JANSSEN,       // Battjes-Janssen (1978)
    THORNTON_GUZA,         // Thornton-Guza (1983)
    ROELVINK               // Roelvink (1993)
};

/**
 * @brief Internal wave model type
 */
enum class InternalWaveModelType {
    LINEAR_MODES,          // Linear normal modes
    WEAKLY_NONLINEAR,      // KdV-type equations
    FULLY_NONLINEAR,       // Non-hydrostatic equations
    PARAMETRIC             // Parameterized effects only
};

/**
 * @brief Sediment transport formula
 */
enum class SedimentTransportFormula {
    ENGELUND_HANSEN,       // Total load
    VAN_RIJN,              // Separate bed/suspended load
    SOULSBY_VAN_RIJN,      // Combined wave-current
    BIJKER,                // Wave-current interaction
    BAILARD,               // Energetics-based
    MEYER_PETER_MULLER     // Bed load only
};

/**
 * @brief Ocean acoustic propagation model
 */
enum class AcousticPropagationModel {
    RAY_TRACING,           // Geometric ray acoustics
    PARABOLIC_EQUATION,    // Parabolic equation method
    NORMAL_MODES,          // Normal mode solution
    WAVENUMBER_INTEGRATION // Full wave solution
};

// =============================================================================
// Wave Model Data Structures
// =============================================================================

/**
 * @brief 2D wave spectrum representation
 */
struct WaveSpectrum2D {
    // Frequency-direction grid
    std::vector<double> frequencies;      // Hz
    std::vector<double> directions;       // radians (nautical convention)
    
    // Energy density E(f, θ) in m²/Hz/rad
    std::vector<std::vector<double>> energy;
    
    // Integral parameters
    double Hs;             // Significant wave height (m)
    double Tp;             // Peak period (s)
    double Tm01;           // Mean period (s)
    double Tm02;           // Zero-crossing period (s)
    double theta_mean;     // Mean wave direction (rad)
    double theta_peak;     // Peak direction (rad)
    double spread;         // Directional spread (rad)
    
    // Methods
    void computeIntegralParameters();
    double getEnergy(double f, double theta) const;  // Interpolated
    double totalEnergy() const;
    void normalize(double target_Hs);
    
    // Generation
    static WaveSpectrum2D generateJONSWAP(double Hs, double Tp, double gamma,
                                          double theta_mean, double spread_deg,
                                          int nf, int nd);
    static WaveSpectrum2D generatePM(double U10, double theta_mean, 
                                     double spread_deg, int nf, int nd);
};

/**
 * @brief Wave field state at a grid point
 */
struct WaveState {
    double Hs;             // Significant wave height (m)
    double Tp;             // Peak period (s)
    double Tm;             // Mean period (s)
    double dir_mean;       // Mean direction (degrees, nautical)
    double dir_peak;       // Peak direction
    double spread;         // Directional spreading (degrees)
    
    // Energy components
    double E_total;        // Total energy (J/m²)
    double E_wind_sea;     // Wind-sea energy
    double E_swell;        // Swell energy
    
    // Radiation stress tensor
    double Sxx, Sxy, Syy;
    
    // Wave-induced mass flux (Stokes drift)
    double U_stokes;       // Surface Stokes drift x (m/s)
    double V_stokes;       // Surface Stokes drift y (m/s)
    
    // Breaking fraction
    double Qb;             // Fraction of breaking waves
};

/**
 * @brief Wave model configuration
 */
struct WaveModelConfig {
    // Grid
    int nx, ny;
    double dx, dy;
    double lon_min, lon_max;
    double lat_min, lat_max;
    
    // Spectral discretization
    int num_frequencies = 36;
    int num_directions = 36;
    double freq_min = 0.03;    // Hz
    double freq_max = 1.0;     // Hz
    
    // Physics
    bool enable_wind_input = true;
    bool enable_whitecapping = true;
    bool enable_quadruplets = true;
    bool enable_depth_breaking = true;
    bool enable_bottom_friction = true;
    bool enable_refraction = true;
    bool enable_diffraction = false;
    
    // Wave breaking
    WaveBreakingModel breaking_model = WaveBreakingModel::BATTJES_JANSSEN;
    double breaking_gamma = 0.73;      // Breaking parameter H/h
    
    // Whitecapping (Komen et al.)
    double whitecap_coeff = 2.36e-5;
    double whitecap_delta = 1.0;
    
    // Bottom friction
    double bottom_friction_coeff = 0.038;  // JONSWAP value
    
    // Time stepping
    double dt = 300.0;         // seconds
    double end_time;
    
    // Output
    int output_interval = 3600;
    std::string output_dir = "output";
};

// =============================================================================
// Wave Model
// =============================================================================

/**
 * @brief Spectral wave model
 *
 * Solves the wave action balance equation:
 *   ∂N/∂t + ∇·(Cg*N) + ∂(Cθ*N)/∂θ + ∂(Cσ*N)/∂σ = S_in + S_nl + S_ds + S_bf
 *
 * Similar to SWAN, WAVEWATCH III approaches.
 */
class WaveModel {
public:
    WaveModel();
    ~WaveModel();
    
    /**
     * @brief Initialize wave model
     */
    void initialize(const WaveModelConfig& config);
    
    /**
     * @brief Set bathymetry
     */
    void setBathymetry(const BathymetryGrid& bathymetry);
    void setBathymetry(const std::vector<double>& depth);
    
    /**
     * @brief Set wind forcing
     */
    void setWind(const std::vector<double>& u10,
                const std::vector<double>& v10);
    void setWindTimeSeries(const std::vector<double>& time,
                          const std::vector<std::vector<double>>& u10,
                          const std::vector<std::vector<double>>& v10);
    
    /**
     * @brief Set current field for wave-current interaction
     */
    void setCurrents(const std::vector<double>& u,
                    const std::vector<double>& v);
    
    /**
     * @brief Set boundary conditions
     */
    void setBoundarySpectrum(int boundary, const WaveSpectrum2D& spectrum);
    void setBoundaryParameters(int boundary, double Hs, double Tp, double dir);
    
    /**
     * @brief Set initial conditions
     */
    void setInitialSpectrum(const WaveSpectrum2D& spectrum);
    void setInitialParameters(double Hs, double Tp, double dir);
    
    /**
     * @brief Time step the model
     */
    double step();
    
    /**
     * @brief Run to specified time
     */
    void run(double end_time);
    
    /**
     * @brief Get wave state at a point
     */
    WaveState getWaveState(int i, int j) const;
    WaveState getWaveState(double lon, double lat) const;
    
    /**
     * @brief Get significant wave height field
     */
    const std::vector<double>& getHs() const { return Hs_field; }
    
    /**
     * @brief Get peak period field
     */
    const std::vector<double>& getTp() const { return Tp_field; }
    
    /**
     * @brief Get mean direction field
     */
    const std::vector<double>& getMeanDirection() const { return dir_field; }
    
    /**
     * @brief Get radiation stress fields
     */
    void getRadiationStress(std::vector<double>& Sxx,
                           std::vector<double>& Sxy,
                           std::vector<double>& Syy) const;
    
    /**
     * @brief Get wave-induced bottom stress
     */
    void getWaveBottomStress(std::vector<double>& tau_w) const;
    
    /**
     * @brief Get Stokes drift
     */
    void getStokesDrift(std::vector<double>& u_stokes,
                       std::vector<double>& v_stokes) const;
    
    /**
     * @brief Compute wave setup/setdown
     */
    void computeWaveSetup(std::vector<double>& setup) const;
    
    /**
     * @brief Get current time
     */
    double getCurrentTime() const { return current_time; }
    
    /**
     * @brief Write output
     */
    void writeOutput(const std::string& filename) const;
    
private:
    WaveModelConfig config;
    
    // Grid
    int nx, ny, nf, nd;
    double dx, dy;
    std::vector<double> freq;          // Frequencies
    std::vector<double> dir;           // Directions
    std::vector<double> df;            // Frequency bandwidth
    std::vector<double> dd;            // Direction bandwidth
    
    // Bathymetry and derived
    std::vector<double> depth;
    std::vector<int> mask;
    
    // Spectral density N(x,y,f,θ) = E/(σ) [action density]
    std::vector<double> action;        // 4D: (nx, ny, nf, nd)
    
    // Group velocity and direction rate
    std::vector<double> Cg_x, Cg_y;    // Group velocity components
    std::vector<double> C_theta;       // Direction rate (refraction)
    std::vector<double> C_sigma;       // Frequency rate (currents)
    
    // Wind
    std::vector<double> u10_field, v10_field;
    
    // Currents
    std::vector<double> u_current, v_current;
    
    // Integrated wave parameters
    std::vector<double> Hs_field, Tp_field, dir_field;
    std::vector<double> Sxx_field, Sxy_field, Syy_field;
    std::vector<double> tau_wave;
    
    // Time
    double current_time;
    
    // Index helper
    int idx(int i, int j, int f, int d) const {
        return i + j * nx + f * nx * ny + d * nx * ny * nf;
    }
    int idx2d(int i, int j) const { return i + j * nx; }
    
    // Physics source terms
    void computeWindInput(std::vector<double>& S_in);
    void computeWhitecapping(std::vector<double>& S_ds);
    void computeQuadruplets(std::vector<double>& S_nl);
    void computeDepthBreaking(std::vector<double>& S_brk);
    void computeBottomFriction(std::vector<double>& S_bf);
    
    // Propagation
    void computeGroupVelocity();
    void computeRefraction();
    void propagate(double dt);
    
    // Integrated parameters
    void computeIntegratedParameters();
    
    // Dispersion relation
    double dispersion(double f, double h) const;
    double groupVelocity(double f, double h) const;
};

// =============================================================================
// Internal Wave Model
// =============================================================================

/**
 * @brief Internal wave dynamics model
 *
 * Models internal gravity wave generation, propagation, and breaking:
 * - Tidal conversion at topography
 * - Wind-generated near-inertial waves
 * - Propagation and reflection
 * - Breaking and mixing
 */
class InternalWaveModel {
public:
    InternalWaveModel();
    
    /**
     * @brief Initialize with stratification profile
     */
    void initialize(int nx, int ny, int nz,
                   double dx, double dy,
                   const std::vector<double>& z_levels,
                   const std::vector<double>& N2_profile);
    
    /**
     * @brief Set density stratification from T,S profiles
     */
    void setStratification(const std::vector<double>& temperature,
                          const std::vector<double>& salinity);
    
    /**
     * @brief Compute normal modes
     */
    void computeNormalModes(int num_modes);
    
    /**
     * @brief Get mode shapes
     */
    void getModeShape(int mode, std::vector<double>& phi) const;
    
    /**
     * @brief Get mode phase speed
     */
    double getModeSpeed(int mode) const;
    
    /**
     * @brief Estimate internal tide generation
     */
    double estimateTidalConversion(const std::vector<double>& bathymetry,
                                   double U_tide, double omega_tide) const;
    
    /**
     * @brief Compute internal wave energy flux
     */
    void computeEnergyFlux(std::vector<double>& Fx, std::vector<double>& Fy) const;
    
    /**
     * @brief Estimate mixing from internal wave breaking
     */
    void estimateMixing(std::vector<double>& Kv_iw) const;
    
    /**
     * @brief Get buoyancy frequency profile
     */
    const std::vector<double>& getBuoyancyFrequency() const { return N2; }
    
private:
    int nx, ny, nz;
    double dx, dy;
    std::vector<double> z;             // Vertical levels
    std::vector<double> N2;            // Buoyancy frequency squared (rad²/s²)
    
    // Normal modes
    int num_modes;
    std::vector<std::vector<double>> mode_shapes;
    std::vector<double> mode_speeds;   // Phase speeds (m/s)
    
    // Internal wave amplitude
    std::vector<double> eta_iw;        // Isopycnal displacement
    std::vector<double> u_iw, v_iw;    // Horizontal velocity
    
    double coriolis_f;
};

// =============================================================================
// Ocean Acoustics Model
// =============================================================================

/**
 * @brief Ocean acoustic propagation model
 *
 * Computes underwater sound propagation:
 * - Sound speed profile effects
 * - Refraction and ducting
 * - Bottom and surface reflection
 * - Transmission loss
 */
class OceanAcousticsModel {
public:
    OceanAcousticsModel();
    
    /**
     * @brief Initialize acoustic environment
     */
    void initialize(double range_max, double depth_max,
                   int nr, int nz);
    
    /**
     * @brief Set sound speed profile
     */
    void setSoundSpeedProfile(const std::vector<double>& z,
                             const std::vector<double>& c);
    
    /**
     * @brief Compute sound speed from T, S, z
     * Uses Mackenzie (1981) or Chen-Millero (1977)
     */
    static double computeSoundSpeed(double T, double S, double z);
    
    /**
     * @brief Set bottom properties
     */
    void setBottomProperties(double density_ratio,
                            double sound_speed_ratio,
                            double attenuation);
    
    /**
     * @brief Set source parameters
     */
    void setSource(double source_depth, double frequency);
    
    /**
     * @brief Compute transmission loss using ray tracing
     */
    void computeRayTracing(int num_rays);
    
    /**
     * @brief Compute transmission loss using parabolic equation
     */
    void computeParabolicEquation();
    
    /**
     * @brief Get transmission loss field (dB re 1m)
     */
    const std::vector<double>& getTransmissionLoss() const { return TL; }
    
    /**
     * @brief Get transmission loss at a point
     */
    double getTransmissionLoss(double range, double depth) const;
    
    /**
     * @brief Compute ray paths
     */
    void getRayPaths(std::vector<std::vector<double>>& ranges,
                    std::vector<std::vector<double>>& depths) const;
    
    /**
     * @brief Estimate detection range for given source level and threshold
     */
    double estimateDetectionRange(double source_level, double noise_level,
                                 double detection_threshold) const;
    
private:
    int nr, nz;
    double dr, dz;
    double range_max, depth_max;
    
    // Sound speed profile
    std::vector<double> z_levels;
    std::vector<double> c_profile;
    
    // Bottom properties
    double rho_bottom;     // Density ratio
    double c_bottom;       // Sound speed ratio
    double alpha_bottom;   // Attenuation (dB/wavelength)
    
    // Source
    double z_source;
    double frequency;
    
    // Transmission loss field
    std::vector<double> TL;
    
    // Ray paths (for visualization)
    std::vector<std::vector<std::pair<double, double>>> ray_paths;
    
    // Interpolate sound speed
    double getSoundSpeed(double depth) const;
    double getSoundSpeedGradient(double depth) const;
};

// =============================================================================
// Sediment Transport Model
// =============================================================================

/**
 * @brief Sediment transport and morphodynamics model
 *
 * Computes coastal sediment transport and bed evolution:
 * - Bed load and suspended load
 * - Wave-current combined transport
 * - Bed morphology changes
 */
class SedimentTransportModel {
public:
    SedimentTransportModel();
    
    /**
     * @brief Initialize model
     */
    void initialize(int nx, int ny, double dx, double dy);
    
    /**
     * @brief Set sediment properties
     */
    void setSedimentProperties(double d50, double d90,
                              double porosity, double density);
    
    /**
     * @brief Set bed composition (for multiple size classes)
     */
    void setBedComposition(const std::vector<double>& size_classes,
                          const std::vector<std::vector<double>>& fractions);
    
    /**
     * @brief Compute transport for given hydrodynamics
     */
    void computeTransport(const std::vector<double>& u,
                         const std::vector<double>& v,
                         const std::vector<double>& h,
                         const WaveModel* waves = nullptr);
    
    /**
     * @brief Update bed morphology
     */
    void updateBed(double dt);
    
    /**
     * @brief Get bed load transport rates
     */
    void getBedLoadTransport(std::vector<double>& qbx,
                            std::vector<double>& qby) const;
    
    /**
     * @brief Get suspended load transport rates
     */
    void getSuspendedLoadTransport(std::vector<double>& qsx,
                                  std::vector<double>& qsy) const;
    
    /**
     * @brief Get bed level change
     */
    const std::vector<double>& getBedChange() const { return dz_bed; }
    
    /**
     * @brief Get current bed level
     */
    const std::vector<double>& getBedLevel() const { return z_bed; }
    
    /**
     * @brief Compute longshore transport rate (CERC formula)
     */
    double computeLongshoreTransport(double Hs, double Tp, double wave_angle,
                                     double depth) const;
    
    /**
     * @brief Compute equilibrium beach profile
     */
    void computeEquilibriumProfile(double A, std::vector<double>& profile) const;
    
private:
    int nx, ny;
    double dx, dy;
    
    // Sediment properties
    double d50;            // Median grain size (m)
    double d90;            // 90th percentile size (m)
    double porosity;       // Bed porosity
    double rho_s;          // Sediment density (kg/m³)
    
    // Transport rates (m²/s)
    std::vector<double> qb_x, qb_y;    // Bed load
    std::vector<double> qs_x, qs_y;    // Suspended load
    
    // Bed level
    std::vector<double> z_bed;
    std::vector<double> dz_bed;        // Bed change
    
    // Transport formula
    SedimentTransportFormula formula = SedimentTransportFormula::SOULSBY_VAN_RIJN;
    
    // Shields parameter
    double computeShieldsParameter(double tau_b, double d) const;
    double criticalShields(double d) const;
    
    // Transport formulas
    void computeVanRijn(const std::vector<double>& tau_b,
                       const std::vector<double>& h);
    void computeSoulsbyVanRijn(const std::vector<double>& u,
                              const std::vector<double>& v,
                              const std::vector<double>& h,
                              const std::vector<double>& tau_w);
};

// =============================================================================
// Thermohaline Model
// =============================================================================

/**
 * @brief Thermohaline circulation model
 *
 * Models temperature and salinity dynamics including:
 * - Surface heat flux
 * - Freshwater flux (E-P)
 * - Deep water formation
 * - Double diffusion
 */
class ThermohalineModel {
public:
    ThermohalineModel();
    
    /**
     * @brief Initialize model
     */
    void initialize(HydrodynamicModel* hydro_model);
    
    /**
     * @brief Set surface heat flux
     */
    void setSurfaceHeatFlux(const std::vector<double>& Q_net);
    
    /**
     * @brief Compute surface heat flux from bulk formulas
     */
    void computeBulkHeatFlux(const std::vector<double>& T_surface,
                            const std::vector<double>& T_air,
                            const std::vector<double>& u10,
                            const std::vector<double>& humidity,
                            const std::vector<double>& cloud_cover,
                            std::vector<double>& Q_net);
    
    /**
     * @brief Set freshwater flux
     */
    void setFreshwaterFlux(const std::vector<double>& E_minus_P);
    
    /**
     * @brief Compute convective adjustment
     */
    void convectiveAdjustment(std::vector<double>& T,
                             std::vector<double>& S);
    
    /**
     * @brief Check for deep water formation
     */
    bool checkDeepWaterFormation(double T_surface, double S_surface,
                                double lat) const;
    
    /**
     * @brief Compute double-diffusive flux
     */
    void computeDoubleDiffusiveFlux(const std::vector<double>& T,
                                   const std::vector<double>& S,
                                   std::vector<double>& F_T,
                                   std::vector<double>& F_S);
    
    /**
     * @brief Get density ratio (Turner angle)
     */
    void computeDensityRatio(const std::vector<double>& T,
                            const std::vector<double>& S,
                            std::vector<double>& Rrho);
    
private:
    HydrodynamicModel* hydro;
    
    // Surface fluxes
    std::vector<double> Q_sw;          // Shortwave radiation (W/m²)
    std::vector<double> Q_lw;          // Longwave radiation
    std::vector<double> Q_sens;        // Sensible heat flux
    std::vector<double> Q_lat;         // Latent heat flux
    std::vector<double> E_P;           // Evaporation - Precipitation (m/s)
    
    // Penetrative shortwave
    bool use_penetrative_sw = true;
    double sw_attenuation_1 = 0.6;     // Fast decay coefficient
    double sw_attenuation_2 = 20.0;    // Slow decay depth scale
    
    // Double diffusion parameters
    double K_T = 1.4e-7;               // Molecular diffusivity of heat
    double K_S = 1.1e-9;               // Molecular diffusivity of salt
};

// =============================================================================
// Air-Sea Interaction Module
// =============================================================================

/**
 * @brief Air-sea interaction and fluxes
 *
 * Computes momentum, heat, and freshwater fluxes at the ocean surface:
 * - Wind stress (with wave-dependent drag)
 * - Heat fluxes (radiative and turbulent)
 * - Evaporation and precipitation
 */
class AirSeaInteraction {
public:
    /**
     * @brief Compute wind stress using COARE algorithm
     */
    static void computeWindStressCOARE(
        const std::vector<double>& u10,
        const std::vector<double>& v10,
        const std::vector<double>& T_air,
        const std::vector<double>& T_sea,
        const std::vector<double>& humidity,
        std::vector<double>& tau_x,
        std::vector<double>& tau_y);
    
    /**
     * @brief Compute heat fluxes
     */
    static void computeHeatFluxes(
        const std::vector<double>& T_air,
        const std::vector<double>& T_sea,
        const std::vector<double>& u10,
        const std::vector<double>& humidity,
        const std::vector<double>& cloud_cover,
        const std::vector<double>& solar_zenith,
        std::vector<double>& Q_sw,
        std::vector<double>& Q_lw,
        std::vector<double>& Q_sens,
        std::vector<double>& Q_lat);
    
    /**
     * @brief Compute evaporation rate
     */
    static void computeEvaporation(
        const std::vector<double>& T_sea,
        const std::vector<double>& T_air,
        const std::vector<double>& u10,
        const std::vector<double>& humidity,
        std::vector<double>& evap);
    
    /**
     * @brief Compute wave-dependent drag coefficient
     */
    static double waveDependentDrag(double u10, double wave_age);
    
    /**
     * @brief Sea surface temperature from skin temperature
     */
    static double skinToSubskinSST(double T_skin, double u10, 
                                  double Q_net, double Q_sw);
};

// =============================================================================
// Ocean Mixing Parameterizations
// =============================================================================

/**
 * @brief Ocean mixing parameterization module
 *
 * Various parameterizations for vertical mixing:
 * - K-Profile Parameterization (KPP)
 * - Mellor-Yamada level 2.5
 * - Generic Length Scale (GLS)
 * - Langmuir turbulence
 */
class OceanMixingModel {
public:
    OceanMixingModel();
    
    /**
     * @brief Compute KPP mixing coefficients
     */
    void computeKPP(const std::vector<double>& T,
                   const std::vector<double>& S,
                   const std::vector<double>& u,
                   const std::vector<double>& v,
                   const std::vector<double>& tau_x,
                   const std::vector<double>& tau_y,
                   const std::vector<double>& Q_net,
                   double f,
                   std::vector<double>& Kv,
                   std::vector<double>& Kt,
                   std::vector<double>& Ks);
    
    /**
     * @brief Compute boundary layer depth (KPP)
     */
    double computeBoundaryLayerDepth(const std::vector<double>& T,
                                    const std::vector<double>& S,
                                    double u_star, double B_f,
                                    double f) const;
    
    /**
     * @brief Compute non-local flux (KPP)
     */
    void computeNonLocalFlux(const std::vector<double>& T,
                            double h_bl, double B_f,
                            std::vector<double>& gamma_T);
    
    /**
     * @brief Compute Langmuir turbulence enhancement
     */
    double langmuirEnhancement(double u_star, double La_t) const;
    
    /**
     * @brief Compute turbulent Langmuir number
     */
    static double computeLangmuirNumber(double u_star, double U_stokes);
    
    /**
     * @brief Apply Mellor-Yamada 2.5 closure
     */
    void computeMellorYamada25(
        const std::vector<double>& T,
        const std::vector<double>& S,
        const std::vector<double>& u,
        const std::vector<double>& v,
        std::vector<double>& q2,
        std::vector<double>& q2l,
        std::vector<double>& Km,
        std::vector<double>& Kh);
    
private:
    // KPP parameters
    double Ri_c = 0.3;         // Critical Richardson number
    double eps = 0.1;          // Surface layer fraction
    double nu_0 = 1.5e-6;      // Molecular viscosity
    double C_v = 1.8;          // Velocity scale coefficient
    double C_s = 98.96;        // Scalar scale coefficient
    
    // MY25 parameters
    double A1 = 0.92;
    double A2 = 0.74;
    double B1 = 16.6;
    double B2 = 10.1;
    double C1 = 0.08;
    double E1 = 1.8;
    double E2 = 1.33;
    double Sq = 0.2;
};

// =============================================================================
// Utility Functions
// =============================================================================

namespace OceanPhysicsUtils {
    /**
     * @brief Compute Coriolis parameter
     */
    double coriolisParameter(double latitude);
    
    /**
     * @brief Compute Ekman depth
     */
    double ekmanDepth(double f, double Az);
    
    /**
     * @brief Compute Rossby deformation radius
     */
    double rossbyRadius(double N, double H, double f);
    
    /**
     * @brief Compute thermal expansion coefficient
     */
    double thermalExpansion(double T, double S, double p);
    
    /**
     * @brief Compute haline contraction coefficient
     */
    double halineContraction(double T, double S, double p);
    
    /**
     * @brief Compute potential temperature
     */
    double potentialTemperature(double T, double S, double p, double p_ref);
    
    /**
     * @brief Compute potential density
     */
    double potentialDensity(double T, double S, double p, double p_ref);
    
    /**
     * @brief Convert wind from 10m to other heights
     */
    double windHeightAdjustment(double u10, double z_target, double z0);
    
    /**
     * @brief Compute saturation specific humidity
     */
    double saturationHumidity(double T, double p);
    
    /**
     * @brief Deep water wave dispersion relation
     */
    double deepWaterWavelength(double T);
    
    /**
     * @brief Shallow water wave dispersion relation (iterative)
     */
    double wavelength(double T, double h, double tol = 1e-6);
    
    /**
     * @brief Ursell number (nonlinearity parameter)
     */
    double ursellNumber(double H, double L, double h);
    
    /**
     * @brief Compute wave setup
     */
    double waveSetup(double Hs, double T, double beach_slope);
}

} // namespace FSRM

#endif // OCEAN_PHYSICS_HPP
