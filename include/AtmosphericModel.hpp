/**
 * @file AtmosphericModel.hpp
 * @brief Comprehensive atmospheric structure and profile models
 * 
 * This file provides detailed atmospheric models that can be configured
 * via external config files. Supports multiple standard atmospheres and
 * custom profiles.
 * 
 * Supported atmosphere models:
 * - US Standard Atmosphere 1976 (US_STANDARD_1976)
 * - ICAO Standard Atmosphere (ICAO_STANDARD)
 * - International Standard Atmosphere (ISA)
 * - NRLMSISE-00 empirical model (NRLMSISE00)
 * - Tropical atmosphere (TROPICAL)
 * - Mid-latitude summer (MIDLATITUDE_SUMMER)
 * - Mid-latitude winter (MIDLATITUDE_WINTER)
 * - Subarctic summer (SUBARCTIC_SUMMER)
 * - Subarctic winter (SUBARCTIC_WINTER)
 * - Custom profile from file (CUSTOM)
 * 
 * Each model provides:
 * - Temperature profile T(z)
 * - Pressure profile P(z)
 * - Density profile ρ(z)
 * - Composition (molecular weights, species concentrations)
 * - Optional: humidity, wind, turbulence profiles
 */

#ifndef ATMOSPHERIC_MODEL_HPP
#define ATMOSPHERIC_MODEL_HPP

#include <string>
#include <vector>
#include <map>
#include <memory>
#include <functional>
#include <cmath>
#include <array>

namespace FSRM {

// =============================================================================
// Physical Constants for Atmospheric Modeling
// =============================================================================
namespace AtmosphereConstants {
    // Universal constants
    constexpr double R_UNIVERSAL = 8.31446261815324;  // J/(mol·K) universal gas constant
    constexpr double AVOGADRO = 6.02214076e23;        // 1/mol
    constexpr double BOLTZMANN = 1.380649e-23;        // J/K
    constexpr double STEFAN_BOLTZMANN = 5.670374419e-8; // W/(m²·K⁴)
    
    // Earth parameters
    constexpr double G0 = 9.80665;                    // m/s² standard gravity
    constexpr double EARTH_RADIUS = 6356766.0;       // m (polar radius for geopotential)
    constexpr double OMEGA_EARTH = 7.2921150e-5;     // rad/s Earth rotation rate
    
    // Standard atmosphere reference values (sea level)
    constexpr double P0_STD = 101325.0;              // Pa
    constexpr double T0_STD = 288.15;                // K (15°C)
    constexpr double RHO0_STD = 1.225;               // kg/m³
    
    // Air composition (dry air by volume)
    constexpr double N2_FRACTION = 0.78084;
    constexpr double O2_FRACTION = 0.20946;
    constexpr double AR_FRACTION = 0.00934;
    constexpr double CO2_FRACTION = 0.000417;        // ~417 ppm (2023 value)
    
    // Molecular weights (g/mol)
    constexpr double M_N2 = 28.0134;
    constexpr double M_O2 = 31.9988;
    constexpr double M_AR = 39.948;
    constexpr double M_CO2 = 44.0095;
    constexpr double M_H2O = 18.01528;
    constexpr double M_AIR_DRY = 28.9644;            // Dry air mean molecular weight
    
    // Specific gas constants (J/(kg·K))
    constexpr double R_AIR_DRY = 287.05287;          // Dry air
    constexpr double R_WATER_VAPOR = 461.5;          // Water vapor
    
    // Thermodynamic properties
    constexpr double CP_AIR = 1004.68506;            // J/(kg·K) specific heat at constant pressure
    constexpr double CV_AIR = 717.62889;             // J/(kg·K) specific heat at constant volume
    constexpr double GAMMA_AIR = 1.4;                // Ratio of specific heats (cp/cv)
    
    // Water properties
    constexpr double LATENT_HEAT_VAP = 2.501e6;      // J/kg at 0°C
    constexpr double LATENT_HEAT_SUB = 2.834e6;      // J/kg sublimation
    
    // Layer boundaries (geopotential altitude, m)
    constexpr double TROPOPAUSE_HEIGHT = 11000.0;
    constexpr double STRATOPAUSE_HEIGHT = 47000.0;
    constexpr double MESOPAUSE_HEIGHT = 86000.0;
    constexpr double THERMOPAUSE_HEIGHT = 500000.0;
}

// =============================================================================
// Atmospheric Layer Definition
// =============================================================================

/**
 * @brief Definition of an atmospheric layer with linear temperature lapse rate
 */
struct AtmosphericLayer {
    std::string name;                // Layer name (e.g., "Troposphere")
    double base_altitude;            // m (geopotential)
    double top_altitude;             // m (geopotential)
    double base_temperature;         // K at base
    double lapse_rate;               // K/m (negative = cooling with altitude)
    double base_pressure;            // Pa at base
    double base_density;             // kg/m³ at base
    
    // Optional composition changes
    double mean_molecular_weight;    // g/mol (can vary with altitude)
    
    // Computed at runtime
    bool isothermal;                 // True if lapse_rate == 0
};

// =============================================================================
// Species Concentration Profile
// =============================================================================

/**
 * @brief Atmospheric species with altitude-dependent concentration
 */
struct AtmosphericSpecies {
    std::string name;                // Species name (e.g., "O3", "H2O")
    std::string formula;             // Chemical formula
    double molecular_weight;         // g/mol
    
    // Concentration profile options
    enum class ProfileType {
        CONSTANT,                    // Constant mixing ratio
        EXPONENTIAL,                 // Exponential decay with altitude
        TABULATED,                   // Tabulated values
        FUNCTION                     // Custom function
    };
    ProfileType profile_type;
    
    // For CONSTANT
    double constant_mixing_ratio;    // Volume mixing ratio
    
    // For EXPONENTIAL
    double surface_concentration;    // kg/m³ or mixing ratio
    double scale_height;             // m
    
    // For TABULATED
    std::vector<double> altitudes;   // m
    std::vector<double> values;      // Concentration values
    
    // For FUNCTION
    std::function<double(double)> profile_function;
    
    // Get concentration at altitude
    double mixingRatio(double altitude) const;
    double numberDensity(double altitude, double total_number_density) const;
    double massDensity(double altitude, double total_density) const;
};

// =============================================================================
// Wind Profile Models
// =============================================================================

/**
 * @brief Wind profile specification
 */
struct WindProfile {
    enum class Type {
        NONE,                        // No wind
        CONSTANT,                    // Constant wind at all altitudes
        LOGARITHMIC,                 // Log-law boundary layer
        POWER_LAW,                   // Power law profile
        EKMAN_SPIRAL,                // Ekman spiral in boundary layer
        TABULATED,                   // Altitude-dependent table
        FUNCTION                     // Custom function
    };
    Type type = Type::NONE;
    
    // Reference values
    double reference_height = 10.0;  // m (standard is 10m)
    double reference_speed = 0.0;    // m/s at reference height
    double reference_direction = 0.0; // degrees from north (meteorological convention)
    
    // For LOGARITHMIC
    double roughness_length = 0.03;  // m (z0)
    double displacement_height = 0.0; // m (d) for urban/forest canopy
    double friction_velocity = 0.0;  // m/s (u*)
    
    // For POWER_LAW
    double power_exponent = 0.14;    // Typical for neutral stability
    
    // For EKMAN_SPIRAL
    double boundary_layer_height = 1000.0;  // m
    double geostrophic_wind = 10.0;  // m/s
    double surface_angle = 45.0;     // degrees from geostrophic
    
    // For TABULATED
    std::vector<double> altitudes;
    std::vector<double> speeds;
    std::vector<double> directions;
    
    // For FUNCTION
    std::function<void(double z, double& u, double& v, double& w)> wind_function;
    
    // Get wind components at altitude
    void getWind(double altitude, double& u, double& v, double& w) const;
    double getSpeed(double altitude) const;
    double getDirection(double altitude) const;
};

// =============================================================================
// Humidity Profile
// =============================================================================

/**
 * @brief Humidity/moisture profile specification
 */
struct HumidityProfile {
    enum class Type {
        NONE,                        // Dry atmosphere
        CONSTANT_RH,                 // Constant relative humidity
        EXPONENTIAL,                 // Exponential decrease with altitude
        TABULATED,                   // Tabulated profile
        FUNCTION                     // Custom function
    };
    Type type = Type::NONE;
    
    // For CONSTANT_RH
    double constant_relative_humidity = 0.0;  // 0-1
    
    // For EXPONENTIAL
    double surface_relative_humidity = 0.7;
    double humidity_scale_height = 2500.0;    // m
    
    // For TABULATED
    std::vector<double> altitudes;
    std::vector<double> relative_humidities;  // 0-1
    
    // For FUNCTION
    std::function<double(double)> humidity_function;  // Returns RH at altitude
    
    // Get humidity at altitude
    double getRelativeHumidity(double altitude) const;
    double getSpecificHumidity(double altitude, double temperature, double pressure) const;
    double getAbsoluteHumidity(double altitude, double temperature, double pressure) const;
    double getVaporPressure(double altitude, double temperature, double pressure) const;
    double getMixingRatio(double altitude, double temperature, double pressure) const;
    
    // Saturation calculations
    static double saturationVaporPressure(double temperature);  // Pa
    static double saturationMixingRatio(double temperature, double pressure);
    static double dewPoint(double temperature, double relative_humidity);
    static double wetBulbTemperature(double temperature, double relative_humidity, double pressure);
};

// =============================================================================
// Turbulence Profile
// =============================================================================

/**
 * @brief Atmospheric turbulence specification
 */
struct TurbulenceProfile {
    enum class Type {
        NONE,
        CONSTANT,                    // Constant TKE
        BOUNDARY_LAYER,              // Standard BL profile
        TABULATED,
        FUNCTION
    };
    Type type = Type::NONE;
    
    // Boundary layer parameters
    double boundary_layer_height = 1000.0;  // m (zi)
    double surface_heat_flux = 0.0;         // W/m² (H0)
    double convective_velocity = 0.0;       // m/s (w*)
    double friction_velocity = 0.3;         // m/s (u*)
    
    // Stability
    double obukhov_length = 1e10;           // m (L), large = neutral
    char pasquill_class = 'D';              // A-F stability class
    
    // For CONSTANT
    double constant_tke = 0.0;              // m²/s²
    double constant_epsilon = 0.0;          // m²/s³ dissipation rate
    
    // For TABULATED
    std::vector<double> altitudes;
    std::vector<double> tke_values;         // m²/s²
    std::vector<double> epsilon_values;     // m²/s³
    
    // Get turbulence at altitude
    double getTKE(double altitude) const;
    double getDissipation(double altitude) const;
    double getEddyViscosity(double altitude) const;
    double getEddyDiffusivity(double altitude) const;
    double getMixingLength(double altitude) const;
};

// =============================================================================
// Atmospheric State at a Point
// =============================================================================

/**
 * @brief Complete atmospheric state at a single point
 */
struct AtmosphericState {
    // Position
    double altitude;                 // m (geometric)
    double geopotential_altitude;    // m
    double latitude;                 // degrees
    double longitude;                // degrees
    
    // Thermodynamic state
    double temperature;              // K
    double pressure;                 // Pa
    double density;                  // kg/m³
    double potential_temperature;    // K
    double virtual_temperature;      // K (accounts for moisture)
    
    // Composition
    double mean_molecular_weight;    // g/mol
    double specific_gas_constant;    // J/(kg·K)
    double gamma;                    // cp/cv ratio
    
    // Moisture
    double relative_humidity;        // 0-1
    double specific_humidity;        // kg/kg
    double mixing_ratio;             // kg/kg
    double vapor_pressure;           // Pa
    double dew_point;                // K
    
    // Derived quantities
    double sound_speed;              // m/s
    double dynamic_viscosity;        // Pa·s
    double kinematic_viscosity;      // m²/s
    double thermal_conductivity;     // W/(m·K)
    double mean_free_path;           // m
    
    // Wind
    double wind_u;                   // m/s (eastward)
    double wind_v;                   // m/s (northward)
    double wind_w;                   // m/s (upward)
    double wind_speed;               // m/s
    double wind_direction;           // degrees
    
    // Turbulence
    double tke;                      // m²/s² turbulent kinetic energy
    double epsilon;                  // m²/s³ dissipation rate
    double eddy_viscosity;           // m²/s
    
    // Stability
    double brunt_vaisala_frequency;  // rad/s (N)
    double richardson_number;        // Ri
    double potential_vorticity;      // PVU
    
    // Radiation-relevant
    double ozone_density;            // kg/m³
    double aerosol_optical_depth;    // dimensionless
    
    // Layer info
    std::string layer_name;
    int layer_index;
};

// =============================================================================
// Base Atmospheric Model Interface
// =============================================================================

/**
 * @brief Abstract base class for atmospheric models
 */
class AtmosphericModelBase {
public:
    virtual ~AtmosphericModelBase() = default;
    
    // Model identification
    virtual std::string getName() const = 0;
    virtual std::string getDescription() const = 0;
    
    // Core profile functions (must be implemented)
    virtual double temperature(double altitude) const = 0;
    virtual double pressure(double altitude) const = 0;
    virtual double density(double altitude) const = 0;
    
    // Optional extended profiles (default implementations provided)
    virtual double meanMolecularWeight(double altitude) const;
    virtual double specificGasConstant(double altitude) const;
    virtual double gamma(double altitude) const;
    
    // Derived quantities (computed from core profiles)
    virtual double soundSpeed(double altitude) const;
    virtual double scaleHeight(double altitude) const;
    virtual double numberDensity(double altitude) const;
    virtual double potentialTemperature(double altitude) const;
    
    // Transport properties
    virtual double dynamicViscosity(double altitude) const;
    virtual double kinematicViscosity(double altitude) const;
    virtual double thermalConductivity(double altitude) const;
    virtual double meanFreePath(double altitude) const;
    
    // Coordinate conversions
    virtual double geometricToGeopotential(double geometric_alt) const;
    virtual double geopotentialToGeometric(double geopotential_alt) const;
    
    // Gravity variation
    virtual double gravity(double altitude) const;
    
    // Get complete state at altitude
    virtual AtmosphericState getState(double altitude) const;
    
    // Layer information
    virtual int getLayerIndex(double altitude) const;
    virtual std::string getLayerName(double altitude) const;
    virtual const std::vector<AtmosphericLayer>& getLayers() const;
    
    // Valid altitude range
    virtual double getMinAltitude() const { return 0.0; }
    virtual double getMaxAltitude() const { return 86000.0; }
    virtual bool isValidAltitude(double altitude) const;
    
    // Wind (if enabled)
    virtual void setWindProfile(const WindProfile& wind);
    virtual void getWind(double altitude, double& u, double& v, double& w) const;
    virtual bool hasWind() const { return has_wind_; }
    
    // Humidity (if enabled)
    virtual void setHumidityProfile(const HumidityProfile& humidity);
    virtual double getRelativeHumidity(double altitude) const;
    virtual double getSpecificHumidity(double altitude) const;
    virtual bool hasHumidity() const { return has_humidity_; }
    
    // Turbulence (if enabled)
    virtual void setTurbulenceProfile(const TurbulenceProfile& turbulence);
    virtual double getTKE(double altitude) const;
    virtual double getEddyViscosity(double altitude) const;
    virtual bool hasTurbulence() const { return has_turbulence_; }
    
    // Configuration
    virtual bool loadFromConfig(const std::string& filename);
    virtual bool saveToConfig(const std::string& filename) const;
    
protected:
    std::vector<AtmosphericLayer> layers_;
    WindProfile wind_profile_;
    HumidityProfile humidity_profile_;
    TurbulenceProfile turbulence_profile_;
    
    bool has_wind_ = false;
    bool has_humidity_ = false;
    bool has_turbulence_ = false;
    
    // Reference values
    double p0_ = AtmosphereConstants::P0_STD;
    double T0_ = AtmosphereConstants::T0_STD;
    double rho0_ = AtmosphereConstants::RHO0_STD;
};

// =============================================================================
// US Standard Atmosphere 1976
// =============================================================================

/**
 * @brief US Standard Atmosphere 1976 model
 * 
 * Valid from 0 to 86 km geometric altitude (0-84.852 km geopotential).
 * Defines temperature as piecewise linear function of geopotential altitude.
 * Pressure and density derived from hydrostatic equation.
 * 
 * Reference: U.S. Standard Atmosphere, 1976, U.S. Government Printing Office
 */
class USStandardAtmosphere1976 : public AtmosphericModelBase {
public:
    USStandardAtmosphere1976();
    
    std::string getName() const override { return "US_STANDARD_1976"; }
    std::string getDescription() const override { 
        return "U.S. Standard Atmosphere 1976 (0-86 km)"; 
    }
    
    double temperature(double altitude) const override;
    double pressure(double altitude) const override;
    double density(double altitude) const override;
    
    double meanMolecularWeight(double altitude) const override;
    
    // Extended model above 86 km (simplified)
    void enableExtendedModel(bool enable) { extended_model_ = enable; }
    double getMaxAltitude() const override { 
        return extended_model_ ? 1000000.0 : 86000.0; 
    }
    
private:
    void initializeLayers();
    void computeBasePressures();
    
    bool extended_model_ = false;
    
    // Layer boundary geopotential altitudes (m)
    static constexpr std::array<double, 8> H_BOUNDS = {
        0.0, 11000.0, 20000.0, 32000.0, 47000.0, 51000.0, 71000.0, 84852.0
    };
    
    // Temperature lapse rates (K/m)
    static constexpr std::array<double, 7> LAPSE_RATES = {
        -0.0065,   // Troposphere
         0.0,      // Tropopause (isothermal)
         0.001,    // Stratosphere lower
         0.0028,   // Stratosphere upper
         0.0,      // Stratopause (isothermal)
        -0.0028,   // Mesosphere lower
        -0.002     // Mesosphere upper
    };
    
    // Base temperatures (K)
    static constexpr std::array<double, 8> T_BASES = {
        288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.946
    };
    
    // Base pressures (computed)
    std::array<double, 8> p_bases_;
};

// =============================================================================
// ICAO Standard Atmosphere
// =============================================================================

/**
 * @brief ICAO Standard Atmosphere
 * 
 * International Civil Aviation Organization standard atmosphere.
 * Very similar to US Standard Atmosphere but with slightly different
 * conventions for aviation applications.
 */
class ICAOStandardAtmosphere : public AtmosphericModelBase {
public:
    ICAOStandardAtmosphere();
    
    std::string getName() const override { return "ICAO_STANDARD"; }
    std::string getDescription() const override { 
        return "ICAO Standard Atmosphere (0-80 km)"; 
    }
    
    double temperature(double altitude) const override;
    double pressure(double altitude) const override;
    double density(double altitude) const override;
    
    // Aviation-specific functions
    double pressureAltitude(double pressure) const;      // Altitude from pressure
    double densityAltitude(double density) const;        // Altitude from density
    double indicatedAltitude(double pressure, double qnh) const;  // With QNH correction
    
private:
    void initializeLayers();
};

// =============================================================================
// NRLMSISE-00 Empirical Model
// =============================================================================

/**
 * @brief NRLMSISE-00 empirical atmosphere model
 * 
 * Naval Research Laboratory Mass Spectrometer and Incoherent Scatter Radar
 * Exosphere model. Valid from ground to exosphere (0-1000 km).
 * Accounts for solar activity, geomagnetic activity, and seasonal variations.
 * 
 * Reference: Picone et al., J. Geophys. Res., 107(A12), 1468, 2002
 */
class NRLMSISE00Atmosphere : public AtmosphericModelBase {
public:
    NRLMSISE00Atmosphere();
    
    std::string getName() const override { return "NRLMSISE00"; }
    std::string getDescription() const override { 
        return "NRLMSISE-00 Empirical Model (0-1000 km)"; 
    }
    
    double temperature(double altitude) const override;
    double pressure(double altitude) const override;
    double density(double altitude) const override;
    
    double meanMolecularWeight(double altitude) const override;
    double getMaxAltitude() const override { return 1000000.0; }
    
    // Species densities
    double N2_density(double altitude) const;
    double O2_density(double altitude) const;
    double O_density(double altitude) const;   // Atomic oxygen
    double He_density(double altitude) const;
    double H_density(double altitude) const;   // Atomic hydrogen
    double Ar_density(double altitude) const;
    double N_density(double altitude) const;   // Atomic nitrogen
    
    // Space weather inputs
    void setF107(double f107);                 // 10.7 cm solar flux
    void setF107Average(double f107a);         // 81-day average
    void setAp(double ap);                     // Geomagnetic index
    void setApArray(const std::array<double, 7>& ap_array);  // Detailed Ap
    
    // Time and location
    void setDateTime(int year, int day_of_year, double seconds_of_day);
    void setLocation(double latitude, double longitude);
    
private:
    void computeProfile();
    
    // Solar and geomagnetic indices
    double f107_ = 150.0;           // Default moderate solar activity
    double f107a_ = 150.0;
    double ap_ = 4.0;               // Default quiet geomagnetic
    std::array<double, 7> ap_array_ = {4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0};
    
    // Time
    int year_ = 2000;
    int day_of_year_ = 172;         // Mid-year
    double seconds_ = 43200.0;      // Noon
    
    // Location
    double latitude_ = 45.0;
    double longitude_ = 0.0;
    
    // Cached profiles
    mutable bool profile_valid_ = false;
    mutable std::vector<double> alt_grid_;
    mutable std::vector<double> T_profile_;
    mutable std::vector<double> rho_profile_;
};

// =============================================================================
// Reference Atmospheres (Tropical, Mid-latitude, Subarctic)
// =============================================================================

/**
 * @brief Reference atmosphere types for different climate zones
 */
enum class ReferenceAtmosphereType {
    TROPICAL,
    MIDLATITUDE_SUMMER,
    MIDLATITUDE_WINTER,
    SUBARCTIC_SUMMER,
    SUBARCTIC_WINTER,
    US_STANDARD  // Same as US Standard 1976
};

/**
 * @brief Reference atmospheres from AFGL/HITRAN
 * 
 * These are the standard reference atmospheres used in radiative transfer
 * calculations, defined by Anderson et al. (1986).
 */
class ReferenceAtmosphere : public AtmosphericModelBase {
public:
    explicit ReferenceAtmosphere(ReferenceAtmosphereType type);
    
    std::string getName() const override;
    std::string getDescription() const override;
    
    double temperature(double altitude) const override;
    double pressure(double altitude) const override;
    double density(double altitude) const override;
    
    // Species profiles (mixing ratios)
    double H2O_mixingRatio(double altitude) const;
    double CO2_mixingRatio(double altitude) const;
    double O3_mixingRatio(double altitude) const;
    double N2O_mixingRatio(double altitude) const;
    double CO_mixingRatio(double altitude) const;
    double CH4_mixingRatio(double altitude) const;
    
    ReferenceAtmosphereType getType() const { return type_; }
    
private:
    void initializeProfiles();
    double interpolateProfile(double altitude, 
                             const std::vector<double>& alts,
                             const std::vector<double>& vals) const;
    
    ReferenceAtmosphereType type_;
    
    // Tabulated profiles
    std::vector<double> altitudes_;      // km
    std::vector<double> temperatures_;   // K
    std::vector<double> pressures_;      // mb
    std::vector<double> h2o_;            // ppmv
    std::vector<double> co2_;            // ppmv
    std::vector<double> o3_;             // ppmv
};

// =============================================================================
// Custom/Tabulated Atmosphere
// =============================================================================

/**
 * @brief Custom atmosphere from tabulated data
 * 
 * Allows defining arbitrary atmospheric profiles from data files or
 * programmatic input.
 */
class TabulatedAtmosphere : public AtmosphericModelBase {
public:
    TabulatedAtmosphere();
    explicit TabulatedAtmosphere(const std::string& config_file);
    
    std::string getName() const override { return name_; }
    std::string getDescription() const override { return description_; }
    
    double temperature(double altitude) const override;
    double pressure(double altitude) const override;
    double density(double altitude) const override;
    
    double getMinAltitude() const override;
    double getMaxAltitude() const override;
    
    // Data input methods
    void setName(const std::string& name) { name_ = name; }
    void setDescription(const std::string& desc) { description_ = desc; }
    
    void setAltitudes(const std::vector<double>& altitudes);
    void setTemperatures(const std::vector<double>& temperatures);
    void setPressures(const std::vector<double>& pressures);
    void setDensities(const std::vector<double>& densities);  // Optional
    
    void addDataPoint(double altitude, double temperature, double pressure,
                      double density = -1.0);  // -1 = compute from ideal gas
    void clearData();
    
    // Load from different formats
    bool loadFromConfig(const std::string& filename) override;
    bool loadFromCSV(const std::string& filename);
    bool loadFromJSON(const std::string& filename);
    
    // Interpolation settings
    enum class InterpolationType {
        LINEAR,
        LOG_LINEAR,       // Linear in log(p), log(rho)
        CUBIC_SPLINE,
        PCHIP             // Piecewise Cubic Hermite
    };
    void setInterpolationType(InterpolationType type) { interp_type_ = type; }
    
    // Extrapolation settings
    enum class ExtrapolationType {
        CONSTANT,         // Hold last value
        LINEAR,           // Linear extrapolation
        EXPONENTIAL,      // Exponential decay (for pressure/density)
        ERROR             // Throw error
    };
    void setExtrapolationType(ExtrapolationType type) { extrap_type_ = type; }
    
private:
    double interpolate(double altitude, const std::vector<double>& values,
                      bool use_log = false) const;
    
    std::string name_ = "CUSTOM";
    std::string description_ = "Custom tabulated atmosphere";
    
    std::vector<double> altitudes_;
    std::vector<double> temperatures_;
    std::vector<double> pressures_;
    std::vector<double> densities_;
    
    InterpolationType interp_type_ = InterpolationType::LINEAR;
    ExtrapolationType extrap_type_ = ExtrapolationType::EXPONENTIAL;
    
    bool has_density_data_ = false;
};

// =============================================================================
// Composite Atmosphere (combines multiple models)
// =============================================================================

/**
 * @brief Composite atmosphere that blends multiple models
 * 
 * Useful for creating smooth transitions between different atmospheric
 * regimes (e.g., standard atmosphere below 86 km, empirical above).
 */
class CompositeAtmosphere : public AtmosphericModelBase {
public:
    CompositeAtmosphere();
    
    std::string getName() const override { return "COMPOSITE"; }
    std::string getDescription() const override;
    
    double temperature(double altitude) const override;
    double pressure(double altitude) const override;
    double density(double altitude) const override;
    
    double getMinAltitude() const override;
    double getMaxAltitude() const override;
    
    // Add component models
    void addModel(std::shared_ptr<AtmosphericModelBase> model,
                  double min_alt, double max_alt,
                  double blend_width = 0.0);  // Smooth transition width
    
    void clearModels();
    
private:
    struct ModelRegion {
        std::shared_ptr<AtmosphericModelBase> model;
        double min_alt;
        double max_alt;
        double blend_width;
    };
    
    std::vector<ModelRegion> models_;
    
    double blendedValue(double altitude,
                       std::function<double(const AtmosphericModelBase*, double)> getter) const;
};

// =============================================================================
// Perturbed Atmosphere (adds perturbations to base model)
// =============================================================================

/**
 * @brief Atmosphere with perturbations added to a base model
 * 
 * Useful for sensitivity studies, adding waves, or creating
 * atmospheric disturbances.
 */
class PerturbedAtmosphere : public AtmosphericModelBase {
public:
    explicit PerturbedAtmosphere(std::shared_ptr<AtmosphericModelBase> base);
    
    std::string getName() const override;
    std::string getDescription() const override;
    
    double temperature(double altitude) const override;
    double pressure(double altitude) const override;
    double density(double altitude) const override;
    
    // Add perturbations
    void addTemperaturePerturbation(std::function<double(double)> dT);
    void addPressurePerturbation(std::function<double(double)> dP);
    void addDensityPerturbation(std::function<double(double)> dRho);
    
    // Convenience: sinusoidal perturbation
    void addGravityWave(double amplitude, double wavelength, double phase = 0.0);
    
    // Convenience: localized perturbation (Gaussian)
    void addLocalizedPerturbation(double amplitude, double center_alt, 
                                  double width, const std::string& variable);
    
    void clearPerturbations();
    
private:
    std::shared_ptr<AtmosphericModelBase> base_model_;
    std::vector<std::function<double(double)>> T_perturbations_;
    std::vector<std::function<double(double)>> P_perturbations_;
    std::vector<std::function<double(double)>> rho_perturbations_;
};

// =============================================================================
// Atmosphere Factory
// =============================================================================

/**
 * @brief Factory for creating atmospheric models
 */
class AtmosphereFactory {
public:
    // Create by name
    static std::shared_ptr<AtmosphericModelBase> create(const std::string& name);
    
    // Create from config file
    static std::shared_ptr<AtmosphericModelBase> createFromConfig(const std::string& filename);
    
    // Create reference atmosphere
    static std::shared_ptr<AtmosphericModelBase> createReference(ReferenceAtmosphereType type);
    
    // Register custom model type
    using ModelCreator = std::function<std::shared_ptr<AtmosphericModelBase>()>;
    static void registerModel(const std::string& name, ModelCreator creator);
    
    // List available models
    static std::vector<std::string> getAvailableModels();
    
private:
    static std::map<std::string, ModelCreator>& getRegistry();
};

// =============================================================================
// Atmosphere Configuration
// =============================================================================

/**
 * @brief Configuration structure for atmospheric models
 */
struct AtmosphereConfig {
    // Model selection
    std::string model_type = "US_STANDARD_1976";
    std::string custom_file = "";
    
    // Reference values (can override defaults)
    double sea_level_pressure = AtmosphereConstants::P0_STD;
    double sea_level_temperature = AtmosphereConstants::T0_STD;
    double sea_level_density = AtmosphereConstants::RHO0_STD;
    
    // Domain
    double min_altitude = 0.0;
    double max_altitude = 86000.0;
    
    // Wind
    bool enable_wind = false;
    WindProfile wind_profile;
    std::string wind_config_file = "";
    
    // Humidity
    bool enable_humidity = false;
    HumidityProfile humidity_profile;
    std::string humidity_config_file = "";
    
    // Turbulence
    bool enable_turbulence = false;
    TurbulenceProfile turbulence_profile;
    std::string turbulence_config_file = "";
    
    // NRLMSISE-00 specific
    double f107 = 150.0;
    double f107a = 150.0;
    double ap = 4.0;
    int year = 2000;
    int day_of_year = 172;
    double local_time = 12.0;  // hours
    double latitude = 45.0;
    double longitude = 0.0;
    
    // Parse from config file
    bool parse(const std::string& filename);
    bool parse(const std::map<std::string, std::string>& config);
    
    // Create model from this configuration
    std::shared_ptr<AtmosphericModelBase> createModel() const;
};

// =============================================================================
// Utility Functions
// =============================================================================

namespace AtmosphereUtils {
    // Coordinate conversions
    double geometricToGeopotential(double geometric_alt);
    double geopotentialToGeometric(double geopotential_alt);
    
    // Gravity
    double gravity(double altitude, double latitude = 45.0);
    
    // Pressure/altitude conversions
    double barometricFormula(double altitude, double T0, double P0, double lapse_rate);
    double hypsometricEquation(double P1, double P2, double T_mean);
    
    // Transport properties (Sutherland's law and others)
    double dynamicViscosity(double temperature);
    double thermalConductivity(double temperature);
    double diffusivity(double temperature, double pressure);
    
    // Sound speed
    double soundSpeed(double temperature, double gamma = AtmosphereConstants::GAMMA_AIR);
    
    // Moisture calculations
    double saturationVaporPressure(double temperature);  // Tetens formula
    double saturationVaporPressureIce(double temperature);
    double virtualTemperature(double temperature, double mixing_ratio);
    double potentialTemperature(double temperature, double pressure);
    double equivalentPotentialTemperature(double temperature, double pressure, double mixing_ratio);
    
    // Stability
    double bruntVaisalaFrequency(double dT_dz, double T, double g = AtmosphereConstants::G0);
    double richardsonNumber(double N2, double dU_dz);
    char pasquillStabilityClass(double wind_speed, double solar_radiation, bool night);
}

} // namespace FSRM

#endif // ATMOSPHERIC_MODEL_HPP
