/**
 * @file MaterialLibrary.hpp
 * @brief Comprehensive Material Library for Geomechanical and Reservoir Simulations
 * 
 * This library provides an extensive database of material properties for:
 * - Igneous, sedimentary, and metamorphic rocks
 * - Unconsolidated sediments and soils
 * - Minerals and mineral assemblages
 * - Fault zone materials
 * - Wellbore materials (cements, casings)
 * - Reference/calibration materials
 * - Specific geological formations worldwide
 * 
 * All properties are in SI units:
 * - Density: kg/m³
 * - Elastic moduli: Pa (GPa where noted)
 * - Stress/Strength: Pa (MPa where noted)
 * - Permeability: m² (mD where noted for input convenience)
 * - Thermal conductivity: W/(m·K)
 * - Specific heat: J/(kg·K)
 * - Fracture toughness: Pa·√m (MPa·√m where noted)
 * 
 * @author FSRM Development Team
 * @copyright Copyright (c) 2024
 */

#ifndef MATERIAL_LIBRARY_HPP
#define MATERIAL_LIBRARY_HPP

#include "MaterialModel.hpp"
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <optional>
#include <functional>
#include <array>
#include <cmath>
#include <algorithm>
#include <cctype>

namespace FSRM {

// =============================================================================
// Material Property Structures
// =============================================================================

/**
 * @brief Complete set of rock/material mechanical properties
 */
struct MechanicalProperties {
    // Basic elastic properties
    double density;                     // kg/m³
    double youngs_modulus;              // Pa
    double poisson_ratio;               // dimensionless
    double bulk_modulus;                // Pa (computed if not specified)
    double shear_modulus;               // Pa (computed if not specified)
    
    // Strength properties
    double tensile_strength;            // Pa
    double compressive_strength;        // Pa (UCS - Unconfined Compressive Strength)
    double cohesion;                    // Pa
    double friction_angle;              // degrees
    double dilation_angle;              // degrees
    
    // Fracture properties
    double fracture_toughness_KIc;      // Pa·√m (Mode I)
    double fracture_toughness_KIIc;     // Pa·√m (Mode II)
    double fracture_energy;             // J/m²
    
    // Poroelastic properties
    double biot_coefficient;            // dimensionless (0-1)
    double biot_modulus;                // Pa
    double grain_bulk_modulus;          // Pa (solid grain modulus)
    
    // Viscoelastic properties
    double relaxation_time;             // s
    double viscosity;                   // Pa·s
    double long_term_modulus;           // Pa (for SLS model)
    
    // Property ranges for uncertainty
    double youngs_modulus_min;          // Pa
    double youngs_modulus_max;          // Pa
    double density_min;                 // kg/m³
    double density_max;                 // kg/m³
    
    // Pressure/temperature dependence coefficients
    double dE_dP;                       // Pressure dependence of E
    double dE_dT;                       // Temperature dependence of E
    
    MechanicalProperties();
    void computeDerivedProperties();
    bool validate() const;
};

/**
 * @brief Hydraulic/transport properties
 */
struct HydraulicProperties {
    double porosity;                    // fraction (0-1)
    double permeability;                // m² (matrix permeability)
    double permeability_x;              // m² (directional)
    double permeability_y;              // m²
    double permeability_z;              // m²
    double compressibility;             // 1/Pa (pore compressibility)
    double specific_storage;            // 1/m
    double tortuosity;                  // dimensionless
    
    // Relative permeability parameters (for multiphase)
    double Swc;                         // Connate water saturation
    double Sor;                         // Residual oil saturation
    double Sgc;                         // Critical gas saturation
    double nw, no, ng;                  // Corey exponents
    double krw_max, kro_max, krg_max;   // Endpoint relative perms
    
    // Capillary pressure
    double Pc_entry;                    // Entry pressure (Pa)
    double lambda_pc;                   // Pore size distribution
    
    // Property ranges
    double porosity_min;
    double porosity_max;
    double permeability_min;            // m²
    double permeability_max;            // m²
    
    // Stress/pressure dependence
    double permeability_stress_coeff;   // For k = k0 * exp(-c*σ')
    double porosity_pressure_coeff;     // For φ = φ0 * (1 - c*P)
    
    HydraulicProperties();
    void setIsotropicPermeability(double k);
    void setAnisotropicPermeability(double kx, double ky, double kz);
};

/**
 * @brief Thermal properties
 */
struct ThermalMaterialProperties {
    double thermal_conductivity;        // W/(m·K)
    double specific_heat;               // J/(kg·K)
    double thermal_expansion;           // 1/K
    double thermal_diffusivity;         // m²/s (computed)
    
    // Anisotropic conductivity
    double conductivity_parallel;       // W/(m·K)
    double conductivity_perpendicular;  // W/(m·K)
    
    // Temperature dependence
    double dK_dT;                       // Temperature dependence of conductivity
    double dCp_dT;                      // Temperature dependence of specific heat
    
    // Radiogenic heat production
    double heat_production;             // W/m³
    
    ThermalMaterialProperties();
    double computeDiffusivity(double density) const;
};

/**
 * @brief Anisotropic elastic properties
 */
struct AnisotropicProperties {
    bool is_anisotropic;
    std::string symmetry_type;          // "isotropic", "VTI", "HTI", "orthorhombic", "triclinic"
    
    // For VTI/HTI
    double E_vertical;                  // Pa
    double E_horizontal;                // Pa
    double nu_vh;                       // Poisson's ratio (vertical/horizontal)
    double nu_hh;                       // Poisson's ratio (horizontal/horizontal)
    double G_vertical;                  // Pa
    
    // Thomsen parameters
    double epsilon;                     // P-wave anisotropy
    double delta;                       // Near-vertical anisotropy
    double gamma;                       // S-wave anisotropy
    
    // Full stiffness tensor (21 components for triclinic)
    std::array<double, 21> Cij;
    
    AnisotropicProperties();
    void setVTI(double Ev, double Eh, double nuvh, double nuhh, double Gv);
    void setThomsen(double Vp0, double Vs0, double eps, double del, double gam, double rho);
};

/**
 * @brief Wave propagation properties
 */
struct SeismicProperties {
    double Vp;                          // P-wave velocity (m/s)
    double Vs;                          // S-wave velocity (m/s)
    double Vp_Vs_ratio;                 // Vp/Vs ratio
    double Qp;                          // P-wave quality factor
    double Qs;                          // S-wave quality factor
    double acoustic_impedance;          // kg/(m²·s)
    
    SeismicProperties();
    void computeFromElastic(double E, double nu, double rho);
    void computeElastic(double rho, double& E, double& nu) const;
};

/**
 * @brief Complete material definition
 */
struct MaterialDefinition {
    std::string name;
    std::string category;               // "igneous", "sedimentary", etc.
    std::string subcategory;            // "volcanic", "plutonic", etc.
    std::string description;
    std::string reference;              // Literature reference
    
    MechanicalProperties mechanical;
    HydraulicProperties hydraulic;
    ThermalMaterialProperties thermal;
    AnisotropicProperties anisotropic;
    SeismicProperties seismic;
    
    // Material behavior flags
    bool is_porous;
    bool is_fractured;
    bool is_viscoelastic;
    bool is_anisotropic_flag;
    
    // Recommended constitutive model
    std::string recommended_model;      // "LINEAR_ELASTIC", "POROELASTIC", etc.
    
    // Create RockProperties from this definition
    RockProperties toRockProperties() const;
    
    // Create config map for use with existing system
    std::map<std::string, std::string> toConfigMap() const;
};

// =============================================================================
// Geological Formation Definitions
// =============================================================================

/**
 * @brief Specific geological formation properties
 */
struct FormationDefinition {
    std::string name;
    std::string basin;                  // Basin name
    std::string region;                 // Geographic region
    std::string age;                    // Geological age
    std::string lithology;              // Primary rock type(s)
    
    // Depth range
    double depth_min;                   // m (TVD)
    double depth_max;                   // m (TVD)
    
    // Representative properties
    MaterialDefinition material;
    
    // In-situ stress gradients
    double Sv_gradient;                 // MPa/km (vertical stress)
    double SHmax_gradient;              // MPa/km (max horizontal stress)
    double Shmin_gradient;              // MPa/km (min horizontal stress)
    double Pp_gradient;                 // MPa/km (pore pressure)
    double temperature_gradient;        // °C/km
    
    // Formation-specific parameters
    double net_to_gross;                // Net pay fraction
    double water_saturation;            // Initial water saturation
    double oil_saturation;              // Initial oil saturation
    double gas_saturation;              // Initial gas saturation
    
    std::string notes;
};

// =============================================================================
// Material Library Class
// =============================================================================

/**
 * @brief Material library providing access to material database
 * 
 * Usage:
 * @code
 * auto& lib = MaterialLibrary::getInstance();
 * auto granite = lib.getMaterial("granite");
 * auto bakken = lib.getFormation("bakken_shale");
 * @endcode
 */
class MaterialLibrary {
public:
    // Singleton access
    static MaterialLibrary& getInstance();
    
    // Get material by name (case-insensitive)
    std::optional<MaterialDefinition> getMaterial(const std::string& name) const;
    
    // Get formation by name
    std::optional<FormationDefinition> getFormation(const std::string& name) const;
    
    // Get all materials in a category
    std::vector<MaterialDefinition> getMaterialsByCategory(const std::string& category) const;
    
    // Get all formation names
    std::vector<std::string> getFormationNames() const;
    
    // Get all material names
    std::vector<std::string> getMaterialNames() const;
    
    // Search materials by property ranges
    std::vector<MaterialDefinition> searchByYoungsModulus(double E_min, double E_max) const;
    std::vector<MaterialDefinition> searchByDensity(double rho_min, double rho_max) const;
    std::vector<MaterialDefinition> searchByPorosity(double phi_min, double phi_max) const;
    std::vector<MaterialDefinition> searchByPermeability(double k_min, double k_max) const;
    
    // Category lists
    std::vector<std::string> getCategories() const;
    std::vector<std::string> getSubcategories(const std::string& category) const;
    
    // Create custom material based on existing
    MaterialDefinition createVariant(const std::string& base_name,
                                     const std::map<std::string, double>& modifications) const;
    
    // Create a RockProperties object directly
    RockProperties createRockProperties(const std::string& material_name) const;
    
    // Export material to config format
    std::string exportToConfig(const std::string& name) const;
    std::string exportToConfig(const MaterialDefinition& mat) const;
    
    // Material interpolation (for layered models)
    MaterialDefinition interpolate(const MaterialDefinition& mat1,
                                   const MaterialDefinition& mat2,
                                   double fraction) const;
    
    // Depth-dependent property computation
    MaterialDefinition getDepthAdjusted(const std::string& name, double depth) const;
    
    // Add custom material to library
    void addMaterial(const MaterialDefinition& material);
    void addFormation(const FormationDefinition& formation);
    
private:
    MaterialLibrary();
    void initializeLibrary();
    void initializeIgneousRocks();
    void initializeSedimentaryRocks();
    void initializeMetamorphicRocks();
    void initializeUnconsolidatedMaterials();
    void initializeMinerals();
    void initializeFaultMaterials();
    void initializeWellboreMaterials();
    void initializeReferenceMaterials();
    void initializeFormations();
    
    std::string toLower(const std::string& s) const;
    
    std::map<std::string, MaterialDefinition> materials_;
    std::map<std::string, FormationDefinition> formations_;
};

// =============================================================================
// Predefined Material Categories
// =============================================================================

namespace Materials {

// --- IGNEOUS ROCKS ---
namespace Igneous {
    MaterialDefinition Granite();
    MaterialDefinition GraniteWeathered();
    MaterialDefinition Granodiorite();
    MaterialDefinition Diorite();
    MaterialDefinition Gabbro();
    MaterialDefinition Basalt();
    MaterialDefinition BasaltVesicular();
    MaterialDefinition Andesite();
    MaterialDefinition Rhyolite();
    MaterialDefinition Obsidian();
    MaterialDefinition Pumice();
    MaterialDefinition Tuff();
    MaterialDefinition TuffWelded();
    MaterialDefinition Dacite();
    MaterialDefinition Peridotite();
    MaterialDefinition Dunite();
    MaterialDefinition Syenite();
    MaterialDefinition Diabase();
    MaterialDefinition Pegmatite();
    MaterialDefinition Porphyry();
    MaterialDefinition Aplite();
    MaterialDefinition Komatiite();
    MaterialDefinition Phonolite();
    MaterialDefinition Latite();
    MaterialDefinition Monzonite();
    MaterialDefinition Troctolite();
    MaterialDefinition Norite();
    MaterialDefinition Anorthosite();
    MaterialDefinition Kimberlite();
    MaterialDefinition Lamprophyre();
}

// --- SEDIMENTARY ROCKS ---
namespace Sedimentary {
    // Clastic rocks
    MaterialDefinition Sandstone();
    MaterialDefinition SandstoneQuartz();
    MaterialDefinition SandstoneFeldspathic();
    MaterialDefinition SandstoneLithic();
    MaterialDefinition SandstoneArkose();
    MaterialDefinition SandstoneGreywacke();
    MaterialDefinition SandstoneTight();
    MaterialDefinition Siltstone();
    MaterialDefinition Shale();
    MaterialDefinition ShaleOrganic();
    MaterialDefinition ShaleSiliceous();
    MaterialDefinition ShaleCalcareous();
    MaterialDefinition Mudstone();
    MaterialDefinition Claystone();
    MaterialDefinition Conglomerate();
    MaterialDefinition Breccia();
    MaterialDefinition Tillite();
    MaterialDefinition Loess();
    MaterialDefinition Turbidite();
    
    // Carbonate rocks
    MaterialDefinition Limestone();
    MaterialDefinition LimestoneOolitic();
    MaterialDefinition LimestoneMicritic();
    MaterialDefinition LimestoneBioclastic();
    MaterialDefinition LimestoneChalk();
    MaterialDefinition LimestoneReef();
    MaterialDefinition LimestoneTravertine();
    MaterialDefinition Dolomite();
    MaterialDefinition DolomiteSucrosic();
    MaterialDefinition Marl();
    MaterialDefinition Coquina();
    
    // Chemical/Evaporite rocks
    MaterialDefinition Halite();
    MaterialDefinition Gypsum();
    MaterialDefinition Anhydrite();
    MaterialDefinition Sylvite();
    MaterialDefinition Potash();
    MaterialDefinition Trona();
    MaterialDefinition Chert();
    MaterialDefinition Flint();
    MaterialDefinition Diatomite();
    MaterialDefinition Ironstone();
    MaterialDefinition BandedIronFormation();
    MaterialDefinition Phosphorite();
    
    // Organic rocks
    MaterialDefinition Coal();
    MaterialDefinition CoalAnthracite();
    MaterialDefinition CoalBituminous();
    MaterialDefinition CoalLignite();
    MaterialDefinition Peat();
    MaterialDefinition OilShale();
}

// --- METAMORPHIC ROCKS ---
namespace Metamorphic {
    // Foliated
    MaterialDefinition Slate();
    MaterialDefinition Phyllite();
    MaterialDefinition Schist();
    MaterialDefinition SchistMica();
    MaterialDefinition SchistChlorite();
    MaterialDefinition SchistTalc();
    MaterialDefinition SchistBlueschist();
    MaterialDefinition Gneiss();
    MaterialDefinition GneissGranitic();
    MaterialDefinition GneissBanded();
    MaterialDefinition Migmatite();
    MaterialDefinition Mylonite();
    MaterialDefinition Cataclasite();
    MaterialDefinition Amphibolite();
    
    // Non-foliated
    MaterialDefinition Marble();
    MaterialDefinition MarbleCalcite();
    MaterialDefinition MarbleDolomitic();
    MaterialDefinition Quartzite();
    MaterialDefinition Hornfels();
    MaterialDefinition Serpentinite();
    MaterialDefinition Soapstone();
    MaterialDefinition Eclogite();
    MaterialDefinition Granulite();
    MaterialDefinition Greenstone();
    MaterialDefinition Skarn();
    MaterialDefinition Tactite();
}

// --- UNCONSOLIDATED MATERIALS ---
namespace Unconsolidated {
    MaterialDefinition Sand();
    MaterialDefinition SandFine();
    MaterialDefinition SandMedium();
    MaterialDefinition SandCoarse();
    MaterialDefinition SandSilty();
    MaterialDefinition SandClayey();
    MaterialDefinition Gravel();
    MaterialDefinition GravelSandy();
    MaterialDefinition Clay();
    MaterialDefinition ClayKaolinite();
    MaterialDefinition ClayIllite();
    MaterialDefinition ClayMontmorillonite();
    MaterialDefinition ClayBentonite();
    MaterialDefinition Silt();
    MaterialDefinition Loam();
    MaterialDefinition Till();
    MaterialDefinition Alluvium();
    MaterialDefinition ColluviumTalus();
    MaterialDefinition Laterite();
    MaterialDefinition Saprolite();
    MaterialDefinition Regolith();
}

// --- MINERALS ---
namespace Minerals {
    // Silicates
    MaterialDefinition Quartz();
    MaterialDefinition Feldspar();
    MaterialDefinition FeldsparOrthoclase();
    MaterialDefinition FeldsparPlagioclase();
    MaterialDefinition Mica();
    MaterialDefinition MicaMuscovite();
    MaterialDefinition MicaBiotite();
    MaterialDefinition Olivine();
    MaterialDefinition Pyroxene();
    MaterialDefinition Amphibole();
    MaterialDefinition Hornblende();
    MaterialDefinition Garnet();
    MaterialDefinition Kyanite();
    MaterialDefinition Sillimanite();
    MaterialDefinition Andalusite();
    MaterialDefinition Tourmaline();
    MaterialDefinition Zircon();
    MaterialDefinition Epidote();
    MaterialDefinition Chlorite();
    MaterialDefinition Serpentine();
    MaterialDefinition Talc();
    MaterialDefinition Zeolites();
    
    // Carbonates
    MaterialDefinition Calcite();
    MaterialDefinition Aragonite();
    MaterialDefinition DolomiteMinite();
    MaterialDefinition Siderite();
    MaterialDefinition Magneite();
    MaterialDefinition Rhodochrosite();
    
    // Sulfates
    MaterialDefinition GypsumMineral();
    MaterialDefinition AnhydriteMineral();
    MaterialDefinition Barite();
    MaterialDefinition Celestite();
    
    // Halides
    MaterialDefinition HaliteMineral();
    MaterialDefinition Fluorite();
    MaterialDefinition SylviteMineral();
    
    // Oxides
    MaterialDefinition Magnetite();
    MaterialDefinition Hematite();
    MaterialDefinition Ilmenite();
    MaterialDefinition Rutile();
    MaterialDefinition Corundum();
    MaterialDefinition Spinel();
    MaterialDefinition Chromite();
    MaterialDefinition Limonite();
    MaterialDefinition Goethite();
    
    // Sulfides
    MaterialDefinition Pyrite();
    MaterialDefinition Pyrrhotite();
    MaterialDefinition Galena();
    MaterialDefinition Sphalerite();
    MaterialDefinition Chalcopyrite();
    MaterialDefinition Molybdenite();
    
    // Native elements
    MaterialDefinition GraphiteMineral();
    MaterialDefinition SulfurNative();
    MaterialDefinition GoldNative();
    MaterialDefinition CopperNative();
    
    // Clays (as minerals)
    MaterialDefinition Kaolinite();
    MaterialDefinition Illite();
    MaterialDefinition Montmorillonite();
    MaterialDefinition Smectite();
    MaterialDefinition Vermiculite();
    MaterialDefinition Palygorskite();
}

// --- FAULT ZONE MATERIALS ---
namespace FaultMaterials {
    MaterialDefinition FaultGouge();
    MaterialDefinition FaultGougeClayRich();
    MaterialDefinition FaultBreccia();
    MaterialDefinition FaultCataclasite();
    MaterialDefinition FaultMylonite();
    MaterialDefinition FaultPseudotachylyte();
    MaterialDefinition FaultDamageZone();
    MaterialDefinition FaultCore();
    MaterialDefinition SlipZone();
    MaterialDefinition ShearedZone();
}

// --- WELLBORE MATERIALS ---
namespace Wellbore {
    MaterialDefinition CementClassA();
    MaterialDefinition CementClassC();
    MaterialDefinition CementClassG();
    MaterialDefinition CementClassH();
    MaterialDefinition CementFoamed();
    MaterialDefinition CementLightweight();
    MaterialDefinition CementHighDensity();
    MaterialDefinition SteelCasing();
    MaterialDefinition SteelTubing();
    MaterialDefinition SteelDrillPipe();
    MaterialDefinition Proppant();
    MaterialDefinition ProppantSand();
    MaterialDefinition ProppantCeramic();
    MaterialDefinition ProppantResinCoated();
    MaterialDefinition DrillCuttings();
}

// --- REFERENCE/CALIBRATION MATERIALS ---
namespace Reference {
    MaterialDefinition Steel();
    MaterialDefinition SteelCarbon();
    MaterialDefinition SteelStainless();
    MaterialDefinition Aluminum();
    MaterialDefinition Concrete();
    MaterialDefinition ConcreteHighStrength();
    MaterialDefinition Glass();
    MaterialDefinition Acrylic();
    MaterialDefinition PMMA();
    MaterialDefinition Copper();
    MaterialDefinition Titanium();
    MaterialDefinition InconelAlloy();
    MaterialDefinition Water();
    MaterialDefinition Ice();
    MaterialDefinition Air();
    MaterialDefinition Rubber();
    MaterialDefinition Epoxy();
    MaterialDefinition CarbonFiber();
    MaterialDefinition Ceramic();
}

// --- PLANETARY MATERIALS ---
namespace Planetary {
    MaterialDefinition LunarRegolith();
    MaterialDefinition LunarBasalt();
    MaterialDefinition LunarHighlands();
    MaterialDefinition MartianRegolith();
    MaterialDefinition MartianBasalt();
    MaterialDefinition Meteorite();
    MaterialDefinition MeteoeriteIron();
    MaterialDefinition MeteoeriteStony();
    MaterialDefinition AsteroidMaterial();
}

} // namespace Materials

// =============================================================================
// Formation Database
// =============================================================================

namespace Formations {

// --- NORTH AMERICAN FORMATIONS ---
MaterialDefinition BakkenShale();
MaterialDefinition EagleFordShale();
MaterialDefinition MarcellUsShale();
MaterialDefinition UticaShale();
MaterialDefinition HaynesvilleShale();
MaterialDefinition BarnnettShale();
MaterialDefinition WoodfordShale();
MaterialDefinition NiobraraChalk();
MaterialDefinition PermianBasinWolfcamp();
MaterialDefinition PermianBasinBoneSpring();
MaterialDefinition PermianBasinSpraberry();
MaterialDefinition DelawareBasin();
MaterialDefinition MidlandBasin();
MaterialDefinition DenverBasin();
MaterialDefinition WillistonBasin();
MaterialDefinition AnadarkoBasin();
MaterialDefinition FortWorthBasin();
MaterialDefinition AppalachianBasin();
MaterialDefinition GulfCoastMiocene();
MaterialDefinition DeepwaterGOM();
MaterialDefinition AlaskanNorthSlope();
MaterialDefinition MontereyFormation();
MaterialDefinition GreenRiverFormation();

// --- INTERNATIONAL FORMATIONS ---
MaterialDefinition NorthSeaChalk();
MaterialDefinition NorthSeaBrentGroup();
MaterialDefinition NorthSeaJurassic();
MaterialDefinition VacasMuertasShale();
MaterialDefinition SantosBasin();
MaterialDefinition CamposBasin();
MaterialDefinition OffshoreAngola();
MaterialDefinition OffshorNigeria();
MaterialDefinition GhawarField();
MaterialDefinition BurganField();
MaterialDefinition WestSiberiaBasin();
MaterialDefinition SongliaoBasin();
MaterialDefinition SichuanBasin();
MaterialDefinition CooperBasin();
MaterialDefinition BrowseBasin();
MaterialDefinition CarnarvonBasin();
MaterialDefinition PerthBasin();
MaterialDefinition TaranakiBasin();

// --- GEOTHERMAL FORMATIONS ---
MaterialDefinition GeysersGeothermal();
MaterialDefinition SaltonSeaGeothermal();
MaterialDefinition NewberryGeothermal();
MaterialDefinition LardereelloGeothermal();
MaterialDefinition WairakeiGeothermal();
MaterialDefinition ReykjanesGeothermal();
MaterialDefinition SoultzGeothermal();
MaterialDefinition HellisheidiGeothermal();
MaterialDefinition CerroPrivetoGeothermal();

// --- CARBON STORAGE FORMATIONS ---
MaterialDefinition SleipnerAquifer();
MaterialDefinition InSalahFormation();
MaterialDefinition QuestFormation();
MaterialDefinition DecaturFormation();
MaterialDefinition GorgoFormation();
MaterialDefinition UtsiraFormation();
MaterialDefinition MtSimonSandstone();
MaterialDefinition SaukSequence();

// --- MINING/ORE BODIES ---
MaterialDefinition WitwatersrandGold();
MaterialDefinition BushveldComplex();
MaterialDefinition KirunaIronOre();
MaterialDefinition OlympicDamCopper();
MaterialDefinition GrasbergCopper();
MaterialDefinition ChuquicamataCopper();

// --- AQUIFER FORMATIONS ---
MaterialDefinition OgallalaAquifer();
MaterialDefinition EdwardsAquifer();
MaterialDefinition FlordianAquifer();
MaterialDefinition HighPlainsAquifer();
MaterialDefinition CentralValleyAquifer();
MaterialDefinition GreatArtesianBasin();
MaterialDefinition NubianSandstoneAquifer();
MaterialDefinition GuaraniAquifer();

// --- SEISMICALLY ACTIVE REGIONS ---
MaterialDefinition SanAndreassFault();
MaterialDefinition HaywardFault();
MaterialDefinition NewMadridZone();
MaterialDefinition WasatchFault();
MaterialDefinition AlpineFaultNZ();
MaterialDefinition NAFZone();
MaterialDefinition JapanTrench();
MaterialDefinition CascadiaSubduction();

} // namespace Formations

// =============================================================================
// Utility Functions
// =============================================================================

/**
 * @brief Convert milliDarcy to m²
 */
inline double mDToM2(double permeability_mD) {
    return permeability_mD * 9.869233e-16;  // 1 mD = 9.869233e-16 m²
}

/**
 * @brief Convert m² to milliDarcy
 */
inline double m2ToMD(double permeability_m2) {
    return permeability_m2 / 9.869233e-16;
}

/**
 * @brief Convert GPa to Pa
 */
inline double GPaToPa(double modulus_GPa) {
    return modulus_GPa * 1e9;
}

/**
 * @brief Convert MPa to Pa
 */
inline double MPaToPa(double stress_MPa) {
    return stress_MPa * 1e6;
}

/**
 * @brief Compute bulk modulus from E and ν
 */
inline double computeBulkModulus(double E, double nu) {
    return E / (3.0 * (1.0 - 2.0 * nu));
}

/**
 * @brief Compute shear modulus from E and ν
 */
inline double computeShearModulus(double E, double nu) {
    return E / (2.0 * (1.0 + nu));
}

/**
 * @brief Compute P-wave velocity from elastic moduli and density
 */
inline double computeVp(double E, double nu, double rho) {
    double K = computeBulkModulus(E, nu);
    double G = computeShearModulus(E, nu);
    return std::sqrt((K + 4.0*G/3.0) / rho);
}

/**
 * @brief Compute S-wave velocity from shear modulus and density
 */
inline double computeVs(double E, double nu, double rho) {
    double G = computeShearModulus(E, nu);
    return std::sqrt(G / rho);
}

/**
 * @brief Compute Biot coefficient from bulk moduli
 */
inline double computeBiotCoefficient(double K_drained, double K_solid) {
    return 1.0 - K_drained / K_solid;
}

/**
 * @brief Compute Biot modulus
 */
inline double computeBiotModulus(double K_solid, double K_fluid, double porosity, double alpha) {
    return 1.0 / ((alpha - porosity) / K_solid + porosity / K_fluid);
}

/**
 * @brief Estimate permeability from porosity using Kozeny-Carman
 */
inline double estimatePermeabilityKC(double porosity, double grain_diameter) {
    // k = d² * φ³ / (180 * (1-φ)²)
    double phi3 = porosity * porosity * porosity;
    double one_minus_phi = 1.0 - porosity;
    return grain_diameter * grain_diameter * phi3 / (180.0 * one_minus_phi * one_minus_phi);
}

/**
 * @brief Estimate thermal diffusivity
 */
inline double computeThermalDiffusivity(double k, double rho, double cp) {
    return k / (rho * cp);
}

} // namespace FSRM

#endif // MATERIAL_LIBRARY_HPP
