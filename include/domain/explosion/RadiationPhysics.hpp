/**
 * @file RadiationPhysics.hpp
 * @brief Comprehensive radiation physics for nuclear explosions and radioactive sources
 * 
 * Implements complete radiation transport and effects including:
 * - Prompt radiation (gamma, neutron) from detonation
 * - Fission product generation and decay chains
 * - Fallout formation, lofting, and deposition
 * - Radioactive dispersal in atmosphere
 * - Dose calculations (external, internal, inhalation)
 * - Shielding and attenuation
 * - Monte Carlo radiation transport
 * 
 * Physical basis and references:
 * - Glasstone & Dolan (1977) - Effects of Nuclear Weapons
 * - ICRP Publication 103 - Radiological Protection
 * - NCRP Report 147 - Structural Shielding
 * - England & Rider (1994) - Fission Product Yields
 * - Way & Wigner (1948) - Fission Product Decay
 * - HOTSPOT/NARAC - Fallout modeling codes
 * - ORIGEN - Nuclear transmutation
 * - RASCAL - Radiological Assessment System
 * 
 * @author FSRM Development Team
 */

#ifndef RADIATION_PHYSICS_HPP
#define RADIATION_PHYSICS_HPP

#include <vector>
#include <array>
#include <string>
#include <memory>
#include <functional>
#include <cmath>
#include <map>
#include <random>
#include <complex>

namespace FSRM {

// Forward declarations
class ExtendedAtmosphericModel;
class MushroomCloudModel;

// =============================================================================
// Physical Constants
// =============================================================================
namespace RadiationConstants {
    // Fundamental constants
    constexpr double SPEED_OF_LIGHT = 2.998e8;            // m/s
    constexpr double ELECTRON_MASS = 9.109e-31;           // kg
    constexpr double ELECTRON_MASS_MEV = 0.511;           // MeV/c²
    constexpr double PROTON_MASS = 1.673e-27;             // kg
    constexpr double NEUTRON_MASS = 1.675e-27;            // kg
    constexpr double AVOGADRO = 6.022e23;                 // 1/mol
    constexpr double PLANCK = 6.626e-34;                  // J·s
    constexpr double BOLTZMANN = 1.381e-23;               // J/K
    constexpr double ELECTRON_CHARGE = 1.602e-19;         // C
    
    // Radiation units
    constexpr double JOULES_PER_MEV = 1.602e-13;          // J/MeV
    constexpr double GRAY_PER_RAD = 0.01;                 // Gy/rad
    constexpr double SIEVERT_PER_REM = 0.01;              // Sv/rem
    constexpr double CURIE_PER_BQ = 2.7e-11;              // Ci/Bq
    constexpr double BQ_PER_CURIE = 3.7e10;               // Bq/Ci
    
    // Nuclear yields
    constexpr double JOULES_PER_KT = 4.184e12;            // J/kt TNT
    constexpr double FISSIONS_PER_KT = 1.45e23;           // fissions/kt
    constexpr double FISSION_ENERGY_MEV = 200.0;          // MeV/fission
    
    // Decay constants
    constexpr double LN2 = 0.693147;                      // ln(2)
    
    // Atmospheric properties
    constexpr double AIR_DENSITY_STP = 1.225;             // kg/m³
    constexpr double AIR_MOLAR_MASS = 28.97e-3;           // kg/mol
}

// =============================================================================
// Enumerations
// =============================================================================

/**
 * @brief Radiation type enumeration
 */
enum class RadiationType {
    GAMMA,                    // Electromagnetic radiation
    NEUTRON,                  // Uncharged nucleons
    ALPHA,                    // He-4 nuclei
    BETA_MINUS,               // Electrons
    BETA_PLUS,                // Positrons
    XRAY,                     // Characteristic/bremsstrahlung X-rays
    PROTON,                   // Charged hydrogen nuclei
    HEAVY_ION                 // Heavier charged particles
};

/**
 * @brief Radiation source type
 */
enum class RadiationSourceType {
    PROMPT_GAMMA,             // Immediate gamma from fission
    PROMPT_NEUTRON,           // Immediate neutrons from fission
    DELAYED_GAMMA,            // Gamma from fission product decay
    DELAYED_NEUTRON,          // Delayed neutrons
    ACTIVATED,                // Neutron activation products
    FALLOUT                   // Deposited fission products
};

/**
 * @brief Fission product category
 */
enum class FissionProductCategory {
    NOBLE_GAS,                // Kr, Xe
    HALOGEN,                  // I, Br
    VOLATILE,                 // Cs, Rb, Te
    REFRACTORY,               // Sr, Ba, Zr, rare earths
    ACTINIDE                  // U, Pu, Am, etc.
};

/**
 * @brief Fallout particle formation mechanism
 */
enum class FalloutFormation {
    CONDENSATION,             // Gas-to-particle condensation
    SURFACE_DEPOSITION,       // Deposition on existing particles
    VOLUME_DISTRIBUTION,      // Incorporated into melted/vaporized material
    FRACTIONATION             // Chemical separation
};

/**
 * @brief Deposition mechanism
 */
enum class DepositionMechanism {
    GRAVITATIONAL_SETTLING,   // Terminal velocity settling
    DRY_DEPOSITION,           // Surface contact without precipitation
    WET_DEPOSITION,           // Washout by rain/snow
    IMPACTION,                // Inertial impaction on surfaces
    DIFFUSION                 // Brownian/turbulent diffusion
};

// =============================================================================
// Radionuclide Data
// =============================================================================

/**
 * @brief Radionuclide properties
 */
struct Radionuclide {
    std::string name;                     // e.g., "Cs-137"
    std::string symbol;                   // e.g., "Cs"
    int atomic_number;                    // Z
    int mass_number;                      // A
    bool is_metastable = false;           // Isomeric state
    
    // Decay properties
    double half_life;                     // seconds
    double decay_constant() const {
        return RadiationConstants::LN2 / half_life;
    }
    
    // Decay modes and branching ratios
    struct DecayMode {
        std::string mode;                 // "beta-", "gamma", "alpha", "EC", "IT"
        double branching_ratio;           // Fraction (0-1)
        std::string daughter;             // Daughter nuclide
        double energy;                    // MeV (mean or specific)
        double intensity;                 // Photons/disintegration
    };
    std::vector<DecayMode> decay_modes;
    
    // Gamma emissions
    struct GammaLine {
        double energy;                    // MeV
        double intensity;                 // Photons per disintegration
        double internal_conversion = 0.0; // IC coefficient
    };
    std::vector<GammaLine> gamma_lines;
    
    // Beta spectrum
    double beta_endpoint_energy = 0.0;    // MeV
    double beta_average_energy = 0.0;     // MeV
    
    // Production
    double fission_yield_thermal = 0.0;   // Fraction per thermal fission (U-235)
    double fission_yield_fast = 0.0;      // Fraction per fast fission
    FissionProductCategory category = FissionProductCategory::REFRACTORY;
    
    // Dose coefficients (ICRP)
    double dose_coefficient_inh = 0.0;    // Sv/Bq (inhalation)
    double dose_coefficient_ing = 0.0;    // Sv/Bq (ingestion)
    double dose_rate_factor = 0.0;        // (Sv/h)/(Bq/m²) ground shine
    
    // Compute specific activity
    double specificActivity() const {
        // A = λN = ln(2)/T × (N_A/M)
        return RadiationConstants::LN2 / half_life * 
               RadiationConstants::AVOGADRO / (mass_number * 1e-3);
    }
};

/**
 * @brief Radionuclide database
 */
class RadionuclideDatabase {
public:
    RadionuclideDatabase();
    
    // Load standard libraries
    void loadDefaultFissionProducts();
    void loadFromFile(const std::string& filename);
    
    // Access
    const Radionuclide* get(const std::string& name) const;
    std::vector<std::string> getFissionProducts() const;
    std::vector<std::string> getByCategory(FissionProductCategory cat) const;
    
    // Decay chain queries
    std::vector<std::string> getDecayChain(const std::string& parent) const;
    double getBranchingRatio(const std::string& parent, const std::string& daughter) const;
    
    // Summary quantities
    double totalFissionYield() const;
    double volatileFraction() const;
    
private:
    std::map<std::string, Radionuclide> database_;
    std::map<std::string, std::vector<std::pair<std::string, double>>> decay_chains_;
    
    void addDefaultNuclide(const std::string& name, int Z, int A, double half_life,
                          double fission_yield, FissionProductCategory cat);
};

// =============================================================================
// Prompt Radiation Model
// =============================================================================

/**
 * @brief Prompt radiation from nuclear detonation
 * 
 * Models immediate gamma and neutron emission including:
 * - Fission gamma spectrum
 * - Fission neutron spectrum (Watt/Maxwellian)
 * - Atmospheric transport and buildup
 * - Dose at distance
 */
class PromptRadiationModel {
public:
    PromptRadiationModel();
    
    // Configuration
    void setYield(double yield_kt);
    void setFissionFraction(double fraction);
    void setBurstHeight(double height_m);
    void setFusionEnhancement(double factor);  // For enhanced radiation weapons
    
    // Spectrum models
    void setGammaSpectrum(std::function<double(double)> spectrum);
    void setNeutronSpectrum(std::function<double(double)> spectrum);
    void useDefaultSpectra();
    
    // Energy output
    double totalGammaEnergy() const;          // MeV
    double totalNeutronEnergy() const;        // MeV
    double gammaYieldFraction() const;        // Fraction of yield
    double neutronYieldFraction() const;
    
    // Spectra (normalized, per MeV)
    double gammaSpectrum(double energy_MeV) const;
    double neutronSpectrum(double energy_MeV) const;
    
    // Watt spectrum for fission neutrons
    double wattSpectrum(double E, double a = 0.988, double b = 2.249) const;
    
    // Atmospheric transport
    double gammaAirAttenuation(double energy_MeV, double distance_m) const;
    double neutronAirAttenuation(double energy_MeV, double distance_m) const;
    double gammaBuildupFactor(double energy_MeV, double mfp) const;
    
    // Dose at distance
    double gammaFluence(double r) const;      // photons/m²
    double neutronFluence(double r) const;    // neutrons/m²
    double gammaDose(double r) const;         // Gy
    double neutronDose(double r) const;       // Gy
    double totalPromptDose(double r) const;   // Gy
    
    // Time-dependent dose rate
    double doseRate(double r, double t) const;  // Gy/s
    
    // Shielded dose
    double shieldedGammaDose(double r, double shield_thickness_m,
                            double shield_density, double shield_Z) const;
    double shieldedNeutronDose(double r, double shield_thickness_m,
                              double hydrogen_density) const;
    
private:
    double yield_kt_;
    double fission_fraction_;
    double burst_height_;
    double fusion_enhancement_ = 1.0;
    
    std::function<double(double)> gamma_spectrum_;
    std::function<double(double)> neutron_spectrum_;
    
    // Attenuation coefficients
    double gammaAttenuationCoeff(double energy_MeV) const;
    double neutronAttenuationCoeff(double energy_MeV) const;
    
    // Dose conversion
    double fluenceToDose(double fluence, RadiationType type, double energy) const;
};

// =============================================================================
// Fission Product Inventory
// =============================================================================

/**
 * @brief Time-dependent fission product inventory
 * 
 * Tracks production and decay of fission products including:
 * - Direct fission yields
 * - Decay chain evolution (Bateman equations)
 * - Activity as function of time
 */
class FissionProductInventory {
public:
    FissionProductInventory();
    
    // Configuration
    void setDatabase(std::shared_ptr<RadionuclideDatabase> db);
    void setFissileType(const std::string& type);  // "U235", "Pu239", etc.
    void setFissionYield(double yield_kt);
    
    // Initialize inventory
    void initialize();
    
    // Time evolution
    void evolve(double dt);
    void evolveTo(double time);
    
    double getCurrentTime() const { return current_time_; }
    
    // Activity queries
    double totalActivity() const;                              // Bq
    double activity(const std::string& nuclide) const;        // Bq
    double activityByCategory(FissionProductCategory cat) const;
    
    // Activity at time
    double totalActivity(double t) const;
    double activity(const std::string& nuclide, double t) const;
    
    // Way-Wigner approximation
    double wayWignerActivity(double t) const;  // Bq (total)
    double wayWignerPower(double t) const;     // W (decay power)
    
    // Isotope mass
    double mass(const std::string& nuclide) const;  // kg
    double totalMass() const;                        // kg
    double massByCategory(FissionProductCategory cat) const;
    
    // Dominant isotopes
    std::vector<std::pair<std::string, double>> topActivityContributors(int n = 10) const;
    std::vector<std::pair<std::string, double>> topDoseContributors(int n = 10) const;
    
    // Specific isotopes of interest
    double iodine131Activity() const;
    double cesium137Activity() const;
    double strontium90Activity() const;
    double xenon133Activity() const;
    
    // Export
    void writeInventory(const std::string& filename) const;
    void writeDecayCurves(const std::string& filename, double t_max, int n_points) const;
    
    // Database access
    std::shared_ptr<RadionuclideDatabase> getDatabase() const { return database_; }
    
private:
    std::shared_ptr<RadionuclideDatabase> database_;
    std::string fissile_type_;
    double fission_yield_kt_;
    double current_time_ = 0.0;
    
    // Nuclide activities (Bq)
    std::map<std::string, double> activities_;
    
    // Initial inventories (atoms)
    std::map<std::string, double> initial_atoms_;
    
    // Solve Bateman equations for decay chain
    void solveBateman(const std::string& parent, double dt);
    
    // Matrix exponential for complex chains
    void solveDecayMatrix(double dt);
};

// =============================================================================
// Fallout Formation Model
// =============================================================================

/**
 * @brief Fallout particle formation physics
 * 
 * Models the creation of radioactive particles including:
 * - Vaporization and condensation
 * - Fractionation (volatile vs refractory separation)
 * - Particle size distribution
 * - Activity incorporation
 */
class FalloutFormationModel {
public:
    /**
     * @brief Fallout particle properties
     */
    struct FalloutParticle {
        double diameter;                  // m
        double density;                   // kg/m³
        double activity;                  // Bq
        double specific_activity;         // Bq/kg
        double mass;                      // kg
        
        std::map<std::string, double> isotope_fractions;
        
        double settlingVelocity(double air_density, double air_viscosity) const;
    };
    
    FalloutFormationModel();
    
    // Configuration
    void setYield(double yield_kt);
    void setFissionYield(double fission_kt);
    void setBurstType(const std::string& type);  // "surface", "low_air", "high_air"
    void setGroundMaterial(double density, double melt_temp, double vaporize_temp);
    
    void setInventory(std::shared_ptr<FissionProductInventory> inventory);
    
    // Compute formation
    void computeParticleFormation();
    
    // Particle size distribution
    void setLogNormalDistribution(double median_diameter, double geometric_std);
    double particleSizePDF(double diameter) const;
    double particleSizeCDF(double diameter) const;
    double medianDiameter() const { return median_diameter_; }
    
    // Fractionation
    double volatileFractionInSmallParticles() const;
    double refractoryFractionInLargeParticles() const;
    double fractionationFactor(const std::string& nuclide, double diameter) const;
    
    // Mass and activity
    double totalFalloutMass() const;            // kg
    double falloutMassInCloud() const;          // kg
    double falloutMassInStem() const;           // kg
    double totalFalloutActivity() const;        // Bq
    
    // Particle generation
    std::vector<FalloutParticle> generateParticles(int n_particles) const;
    
    // Soil activation (neutron-induced)
    double soilActivationActivity() const;       // Bq (Na-24, Mn-56, etc.)
    
private:
    double yield_kt_;
    double fission_yield_kt_;
    std::string burst_type_;
    
    // Ground material
    double soil_density_ = 2500.0;
    double soil_melt_temp_ = 1500.0;
    double soil_vaporize_temp_ = 3000.0;
    
    std::shared_ptr<FissionProductInventory> inventory_;
    
    // Particle distribution parameters
    double median_diameter_ = 100e-6;         // 100 μm
    double geometric_std_ = 2.5;
    
    // Mass partitioning
    double mass_in_cloud_ = 0.0;
    double mass_in_stem_ = 0.0;
    double total_mass_ = 0.0;
    
    // Crater and fireball parameters
    double craterMass() const;
    double fireballMass() const;
    double stemEntrainmentMass() const;
    
    // Random number generator for particle generation
    mutable std::mt19937 rng_;
};

// =============================================================================
// Radioactive Dispersal Model
// =============================================================================

/**
 * @brief Atmospheric dispersal of radioactive material
 * 
 * Combines:
 * - Gaussian plume model for near-field
 * - Lagrangian particle tracking for far-field
 * - Wet and dry deposition
 * - Time-dependent meteorology
 */
class RadioactiveDispersalModel {
public:
    /**
     * @brief Deposition record
     */
    struct DepositionRecord {
        double x, y;                      // Location (m)
        double time;                      // Deposition time (s)
        double mass;                      // kg
        double activity;                  // Bq at deposition time
        double particle_diameter;          // m
        DepositionMechanism mechanism;
    };
    
    RadioactiveDispersalModel();
    
    // Configuration
    void setAtmosphere(std::shared_ptr<ExtendedAtmosphericModel> atm);
    void setCloudModel(std::shared_ptr<MushroomCloudModel> cloud);
    void setFalloutModel(std::shared_ptr<FalloutFormationModel> fallout);
    void setInventory(std::shared_ptr<FissionProductInventory> inventory);
    
    // Wind field
    void setConstantWind(double speed, double direction);
    void setWindProfile(std::function<void(double z, double& u, double& v)> wind);
    void setTimeVaryingWind(std::function<void(double t, double z, double& u, double& v)> wind);
    
    // Precipitation
    void setRainfall(double rate_mm_hr);
    void setRainfallField(std::function<double(double x, double y, double t)> rain);
    void setSnowfall(double rate_mm_hr);
    
    // Dispersion parameters
    void setDispersionScheme(const std::string& scheme);  // "gaussian", "lagrangian", "hybrid"
    void setPasquillStability(char stability_class);      // A-F
    void setTurbulentDiffusivity(double Kh, double Kv);
    
    // Deposition parameters
    void setDryDepositionVelocity(double vd);             // m/s
    void setWashoutCoefficient(double lambda);            // 1/s
    
    // Initialize and run
    void initialize();
    void step(double dt);
    void runTo(double time);
    
    double getCurrentTime() const { return current_time_; }
    
    // Airborne concentration
    double airConcentration(double x, double y, double z, double t) const;  // Bq/m³
    double airConcentrationByNuclide(const std::string& nuclide,
                                     double x, double y, double z, double t) const;
    double integratedAirConcentration(double x, double y, double z,
                                      double t_start, double t_end) const;  // Bq·s/m³
    
    // Ground deposition
    double groundDeposition(double x, double y, double t) const;           // Bq/m²
    double groundDepositionByNuclide(const std::string& nuclide,
                                     double x, double y, double t) const;
    double depositionRate(double x, double y, double t) const;             // Bq/(m²·s)
    
    // Plume characteristics
    void getPlumeContour(double t, double level, 
                        std::vector<std::pair<double,double>>& contour) const;
    double plumeWidth(double downwind_distance) const;
    double plumeHeight(double downwind_distance) const;
    double plumeCenterline(double downwind_distance) const;
    
    // Fallout pattern
    void getFalloutContours(double t, const std::vector<double>& levels,
                           std::vector<std::vector<std::pair<double,double>>>& contours) const;
    double falloutArrivalTime(double x, double y) const;
    double peakFalloutRate(double x, double y) const;
    
    // Total deposited
    double totalDepositedActivity() const;                 // Bq
    double totalDepositedMass() const;                     // kg
    
    // Output
    void writeConcentrationField(const std::string& filename, double t) const;
    void writeDepositionField(const std::string& filename, double t) const;
    void writeDepositionHistory(const std::string& filename) const;
    
private:
    std::shared_ptr<ExtendedAtmosphericModel> atmosphere_;
    std::shared_ptr<MushroomCloudModel> cloud_;
    std::shared_ptr<FalloutFormationModel> fallout_model_;
    std::shared_ptr<FissionProductInventory> inventory_;
    
    // Wind
    std::function<void(double, double, double&, double&)> wind_func_;
    bool time_varying_wind_ = false;
    
    // Precipitation
    std::function<double(double, double, double)> rainfall_func_;
    double constant_rainfall_ = 0.0;
    
    // Dispersion parameters
    std::string dispersion_scheme_ = "hybrid";
    char stability_class_ = 'D';
    double Kh_ = 100.0;               // m²/s horizontal
    double Kv_ = 10.0;                // m²/s vertical
    
    // Deposition parameters
    double dry_deposition_velocity_ = 0.01;  // m/s
    double washout_coefficient_ = 1e-4;       // 1/s per mm/hr rain
    
    // State
    double current_time_ = 0.0;
    
    // Particle tracking (Lagrangian)
    struct TrackedParticle {
        double x, y, z;
        double mass;
        double activity;
        double diameter;
        bool deposited = false;
    };
    std::vector<TrackedParticle> particles_;
    
    // Deposition records
    std::vector<DepositionRecord> depositions_;
    
    // Grid for concentration/deposition (Eulerian)
    std::vector<double> x_grid_, y_grid_, z_grid_;
    std::vector<std::vector<std::vector<double>>> concentration_grid_;
    std::vector<std::vector<double>> deposition_grid_;
    
    // Internal methods
    void advectParticles(double dt);
    void diffuseParticles(double dt);
    void settleParticles(double dt);
    void depositParticles(double dt);
    
    // Pasquill-Gifford dispersion
    double sigmaY(double x, char stability) const;
    double sigmaZ(double x, char stability) const;
    
    // Gaussian plume
    double gaussianPlumeConcentration(double x, double y, double z,
                                      double Q, double H, double u) const;
};

// =============================================================================
// Dose Calculation Model
// =============================================================================

/**
 * @brief Radiation dose calculations
 * 
 * Computes doses from multiple pathways:
 * - External (cloud shine, ground shine)
 * - Internal (inhalation, ingestion)
 * - Prompt radiation
 * 
 * Implements ICRP methodology for dose assessment.
 */
class DoseCalculationModel {
public:
    /**
     * @brief Dose result structure
     */
    struct DoseResult {
        // External doses
        double prompt_gamma_dose = 0.0;       // Sv
        double prompt_neutron_dose = 0.0;     // Sv
        double cloud_shine_dose = 0.0;        // Sv
        double ground_shine_dose = 0.0;       // Sv
        
        // Internal doses
        double inhalation_dose = 0.0;         // Sv (committed effective)
        double ingestion_dose = 0.0;          // Sv (committed effective)
        double wound_dose = 0.0;              // Sv
        
        // Totals
        double total_external() const { 
            return prompt_gamma_dose + prompt_neutron_dose + 
                   cloud_shine_dose + ground_shine_dose; 
        }
        double total_internal() const { 
            return inhalation_dose + ingestion_dose + wound_dose; 
        }
        double total_effective_dose() const { 
            return total_external() + total_internal(); 
        }
        
        // Organ doses (example organs)
        double red_marrow_dose = 0.0;
        double thyroid_dose = 0.0;
        double lung_dose = 0.0;
        double gonads_dose = 0.0;
        double skin_dose = 0.0;
    };
    
    /**
     * @brief Shielding configuration
     */
    struct ShieldingConfig {
        double wall_thickness = 0.0;          // m
        double wall_density = 2400.0;         // kg/m³ (concrete)
        double wall_atomic_number = 14.0;     // Effective Z
        double roof_thickness = 0.0;
        double ground_roughness = 1.0;        // Roughness factor
        double building_factor = 1.0;         // Location factor (inside/outside)
        bool underground = false;
        double underground_depth = 0.0;       // m
    };
    
    DoseCalculationModel();
    
    // Configure models
    void setPromptModel(std::shared_ptr<PromptRadiationModel> prompt);
    void setDispersalModel(std::shared_ptr<RadioactiveDispersalModel> dispersal);
    void setInventory(std::shared_ptr<FissionProductInventory> inventory);
    void setDatabase(std::shared_ptr<RadionuclideDatabase> db);
    
    // Receptor parameters
    void setBreathingRate(double rate);       // m³/s (default: 2.5e-4 light activity)
    void setExposureDuration(double duration); // s
    void setShielding(const ShieldingConfig& shielding);
    
    // Single-point dose calculation
    DoseResult calculateDose(double x, double y, double z, double t_start, double t_end) const;
    
    // External doses
    double promptGammaDose(double r) const;
    double promptNeutronDose(double r) const;
    double cloudShineDose(double x, double y, double z, 
                         double t_start, double t_end) const;
    double groundShineDose(double x, double y, double t_start, double t_end) const;
    
    // Internal doses
    double inhalationDose(double x, double y, double z,
                         double t_start, double t_end) const;
    double thyroidInhalationDose(double x, double y, double z,
                                 double t_start, double t_end) const;  // I-131 specific
    
    // Dose rate
    double totalDoseRate(double x, double y, double z, double t) const;  // Sv/h
    double groundShineDoseRate(double x, double y, double t) const;
    double cloudShineDoseRate(double x, double y, double z, double t) const;
    
    // Dose conversion
    double doseConversionFactor(const std::string& nuclide, 
                                RadiationType type, double energy) const;
    double effectiveDosePerUnitIntake(const std::string& nuclide,
                                      const std::string& pathway) const;
    
    // Shielding
    double gammaShieldingFactor(double thickness, double density, double energy) const;
    double buildingProtectionFactor(const ShieldingConfig& config) const;
    double groundRoughnessFactor(double roughness) const;
    
    // Field calculations
    void computeDoseField(double z_height, double t_exposure,
                         const std::vector<double>& x_grid,
                         const std::vector<double>& y_grid,
                         std::vector<std::vector<double>>& dose_field) const;
    
    void writeDoseMap(const std::string& filename, double z, 
                     double t_start, double t_end) const;
    
private:
    std::shared_ptr<PromptRadiationModel> prompt_model_;
    std::shared_ptr<RadioactiveDispersalModel> dispersal_model_;
    std::shared_ptr<FissionProductInventory> inventory_;
    std::shared_ptr<RadionuclideDatabase> database_;
    
    double breathing_rate_ = 2.5e-4;      // m³/s
    double exposure_duration_ = 3600.0;    // s (1 hour default)
    ShieldingConfig shielding_;
    
    // Internal computation helpers
    double integrateGroundShine(double x, double y, double t_start, double t_end) const;
    double integrateCloudShine(double x, double y, double z, 
                               double t_start, double t_end) const;
    double integrateInhalation(double x, double y, double z,
                               double t_start, double t_end) const;
};

// =============================================================================
// Monte Carlo Radiation Transport
// =============================================================================

/**
 * @brief Monte Carlo particle transport for detailed radiation calculations
 * 
 * Provides accurate transport through complex geometry and materials.
 */
class MonteCarloTransport {
public:
    /**
     * @brief Material definition
     */
    struct Material {
        std::string name;
        double density;                       // kg/m³
        std::vector<std::pair<int, double>> composition;  // (Z, weight_fraction)
        
        // Computed properties
        double macroscopicCrossSection(RadiationType type, double energy) const;
    };
    
    /**
     * @brief Geometry region
     */
    struct Region {
        std::string name;
        std::string material_name;
        std::function<bool(double, double, double)> inside;  // Point-in-region test
    };
    
    /**
     * @brief Particle history
     */
    struct ParticleHistory {
        RadiationType type;
        double energy;
        std::array<double, 3> position;
        std::array<double, 3> direction;
        double weight = 1.0;
        int generation = 0;
        bool alive = true;
    };
    
    MonteCarloTransport();
    
    // Geometry setup
    void addMaterial(const Material& mat);
    void addRegion(const Region& region);
    void setSourceRegion(const std::string& region);
    
    // Source definition
    void setPointSource(double x, double y, double z);
    void setDistributedSource(std::function<void(double&, double&, double&)> sample);
    void setSourceSpectrum(RadiationType type, std::function<double()> sample_energy);
    
    // Physics options
    void enablePhotoelectric(bool enable);
    void enableComptonScattering(bool enable);
    void enablePairProduction(bool enable);
    void enableNeutronScattering(bool enable);
    void enableNeutronCapture(bool enable);
    
    // Variance reduction
    void setImportanceMap(std::function<double(double, double, double)> importance);
    void enableSplitting(bool enable, double threshold);
    void enableRussianRoulette(bool enable, double threshold);
    
    // Run simulation
    void run(int n_histories);
    
    // Results
    double fluenceAt(double x, double y, double z) const;
    double doseAt(double x, double y, double z) const;
    double depositionAt(double x, double y, double z) const;
    
    // Statistics
    double relativeError() const;
    int historyCount() const { return total_histories_; }
    
private:
    std::map<std::string, Material> materials_;
    std::vector<Region> regions_;
    std::string source_region_;
    
    // Source
    std::function<void(double&, double&, double&)> source_position_;
    std::function<double()> source_energy_;
    RadiationType source_type_ = RadiationType::GAMMA;
    
    // Physics enables
    bool photoelectric_ = true;
    bool compton_ = true;
    bool pair_production_ = true;
    bool neutron_scatter_ = true;
    bool neutron_capture_ = true;
    
    // Variance reduction
    std::function<double(double, double, double)> importance_map_;
    bool splitting_ = false;
    bool russian_roulette_ = false;
    double split_threshold_ = 2.0;
    double rr_threshold_ = 0.5;
    
    // Results
    int total_histories_ = 0;
    std::vector<double> tally_fluence_;
    std::vector<double> tally_dose_;
    std::vector<double> tally_squared_;
    
    // Tally mesh
    std::vector<double> x_mesh_, y_mesh_, z_mesh_;
    
    // Transport methods
    void transportParticle(ParticleHistory& p);
    double samplePathLength(const Material& mat, RadiationType type, double energy);
    void sampleInteraction(ParticleHistory& p, const Material& mat);
    std::string findRegion(double x, double y, double z) const;
    
    // Random number generator
    std::mt19937_64 rng_;
    std::uniform_real_distribution<double> uniform_;
};

// =============================================================================
// Configuration Parser
// =============================================================================

/**
 * @brief Parse radiation physics configuration
 */
class RadiationPhysicsConfig {
public:
    RadiationPhysicsConfig() = default;
    
    bool parse(const std::string& filename);
    bool parse(const std::map<std::string, std::string>& config);
    
    // Create models
    std::unique_ptr<PromptRadiationModel> createPromptModel() const;
    std::unique_ptr<FissionProductInventory> createInventory() const;
    std::unique_ptr<FalloutFormationModel> createFalloutModel() const;
    std::unique_ptr<RadioactiveDispersalModel> createDispersalModel() const;
    std::unique_ptr<DoseCalculationModel> createDoseModel() const;
    
    // Access configuration
    double getYield() const { return yield_kt_; }
    double getFissionFraction() const { return fission_fraction_; }
    
private:
    double yield_kt_ = 100.0;
    double fission_fraction_ = 0.5;
    double burst_height_ = 0.0;
    std::string burst_type_ = "surface";
    
    // Dispersion
    std::string dispersion_scheme_ = "hybrid";
    double wind_speed_ = 10.0;
    double wind_direction_ = 0.0;
    char stability_class_ = 'D';
    double rainfall_rate_ = 0.0;
    
    // Fallout
    double median_particle_size_ = 100e-6;
    double geometric_std_ = 2.5;
    
    bool configured_ = false;
};

// =============================================================================
// Convenience Functions
// =============================================================================

namespace RadiationUtils {
    // Activity decay
    double decayActivity(double A0, double half_life, double time);
    
    // Way-Wigner approximation
    double wayWignerActivity(double fissions, double time);  // Bq
    double wayWignerPower(double fissions, double time);     // W
    
    // Settling velocity (Stokes)
    double stokesSettlingVelocity(double diameter, double rho_particle, 
                                  double rho_air, double mu_air);
    
    // Gaussian dispersion coefficients
    double pasquillGiffordSigmaY(double x, char stability);
    double pasquillGiffordSigmaZ(double x, char stability);
    
    // Dose unit conversions
    double grayToRem(double Gy);
    double sievertToRem(double Sv);
    double roentgenToGray(double R);
    
    // Half-value layers
    double gammaHVL(double energy_MeV, const std::string& material);
    double neutronHVL(double energy_MeV, const std::string& material);
    
    // Atmospheric stability
    char stabilityFromMeteo(double wind_speed, double solar_rad, bool night);
}

} // namespace FSRM

#endif // RADIATION_PHYSICS_HPP
