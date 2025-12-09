/**
 * @file MaterialLibrary.cpp
 * @brief Implementation of the Comprehensive Material Library
 * 
 * Part 1: Struct constructors and methods
 */

#include "MaterialLibrary.hpp"
#include <sstream>
#include <iomanip>
#include <stdexcept>

namespace FSRM {

// =============================================================================
// MechanicalProperties Implementation
// =============================================================================

MechanicalProperties::MechanicalProperties()
    : density(2650.0)
    , youngs_modulus(50.0e9)
    , poisson_ratio(0.25)
    , bulk_modulus(0.0)
    , shear_modulus(0.0)
    , tensile_strength(10.0e6)
    , compressive_strength(100.0e6)
    , cohesion(20.0e6)
    , friction_angle(30.0)
    , dilation_angle(10.0)
    , fracture_toughness_KIc(1.0e6)
    , fracture_toughness_KIIc(1.5e6)
    , fracture_energy(100.0)
    , biot_coefficient(0.7)
    , biot_modulus(10.0e9)
    , grain_bulk_modulus(36.0e9)
    , relaxation_time(1.0e6)
    , viscosity(1.0e18)
    , long_term_modulus(40.0e9)
    , youngs_modulus_min(0.0)
    , youngs_modulus_max(0.0)
    , density_min(0.0)
    , density_max(0.0)
    , dE_dP(0.0)
    , dE_dT(-1.0e7)
{
    computeDerivedProperties();
}

void MechanicalProperties::computeDerivedProperties() {
    if (bulk_modulus <= 0.0 && youngs_modulus > 0.0 && poisson_ratio >= -1.0 && poisson_ratio < 0.5) {
        bulk_modulus = computeBulkModulus(youngs_modulus, poisson_ratio);
    }
    if (shear_modulus <= 0.0 && youngs_modulus > 0.0 && poisson_ratio > -1.0) {
        shear_modulus = computeShearModulus(youngs_modulus, poisson_ratio);
    }
    if (youngs_modulus_min <= 0.0) {
        youngs_modulus_min = youngs_modulus * 0.8;
    }
    if (youngs_modulus_max <= 0.0) {
        youngs_modulus_max = youngs_modulus * 1.2;
    }
    if (density_min <= 0.0) {
        density_min = density * 0.95;
    }
    if (density_max <= 0.0) {
        density_max = density * 1.05;
    }
}

bool MechanicalProperties::validate() const {
    if (density <= 0.0) return false;
    if (youngs_modulus <= 0.0) return false;
    if (poisson_ratio <= -1.0 || poisson_ratio >= 0.5) return false;
    if (friction_angle < 0.0 || friction_angle > 90.0) return false;
    if (biot_coefficient < 0.0 || biot_coefficient > 1.0) return false;
    return true;
}

// =============================================================================
// HydraulicProperties Implementation
// =============================================================================

HydraulicProperties::HydraulicProperties()
    : porosity(0.1)
    , permeability(1.0e-15)
    , permeability_x(1.0e-15)
    , permeability_y(1.0e-15)
    , permeability_z(1.0e-15)
    , compressibility(1.0e-9)
    , specific_storage(1.0e-6)
    , tortuosity(2.0)
    , Swc(0.2)
    , Sor(0.2)
    , Sgc(0.05)
    , nw(2.0)
    , no(2.0)
    , ng(2.0)
    , krw_max(0.3)
    , kro_max(1.0)
    , krg_max(0.8)
    , Pc_entry(10000.0)
    , lambda_pc(2.0)
    , porosity_min(0.0)
    , porosity_max(0.0)
    , permeability_min(0.0)
    , permeability_max(0.0)
    , permeability_stress_coeff(1.0e-8)
    , porosity_pressure_coeff(1.0e-10)
{
    if (porosity_min <= 0.0) porosity_min = porosity * 0.8;
    if (porosity_max <= 0.0) porosity_max = porosity * 1.2;
    if (permeability_min <= 0.0) permeability_min = permeability * 0.1;
    if (permeability_max <= 0.0) permeability_max = permeability * 10.0;
}

void HydraulicProperties::setIsotropicPermeability(double k) {
    permeability = k;
    permeability_x = k;
    permeability_y = k;
    permeability_z = k;
}

void HydraulicProperties::setAnisotropicPermeability(double kx, double ky, double kz) {
    permeability_x = kx;
    permeability_y = ky;
    permeability_z = kz;
    permeability = std::cbrt(kx * ky * kz);  // Geometric mean
}

// =============================================================================
// ThermalMaterialProperties Implementation
// =============================================================================

ThermalMaterialProperties::ThermalMaterialProperties()
    : thermal_conductivity(2.5)
    , specific_heat(900.0)
    , thermal_expansion(1.0e-5)
    , thermal_diffusivity(0.0)
    , conductivity_parallel(2.5)
    , conductivity_perpendicular(2.5)
    , dK_dT(-0.001)
    , dCp_dT(0.1)
    , heat_production(1.0e-6)
{
}

double ThermalMaterialProperties::computeDiffusivity(double density) const {
    if (density > 0.0 && specific_heat > 0.0) {
        return thermal_conductivity / (density * specific_heat);
    }
    return 0.0;
}

// =============================================================================
// AnisotropicProperties Implementation
// =============================================================================

AnisotropicProperties::AnisotropicProperties()
    : is_anisotropic(false)
    , symmetry_type("isotropic")
    , E_vertical(50.0e9)
    , E_horizontal(50.0e9)
    , nu_vh(0.25)
    , nu_hh(0.25)
    , G_vertical(20.0e9)
    , epsilon(0.0)
    , delta(0.0)
    , gamma(0.0)
    , Cij{}
{
    Cij.fill(0.0);
}

void AnisotropicProperties::setVTI(double Ev, double Eh, double nuvh, double nuhh, double Gv) {
    is_anisotropic = true;
    symmetry_type = "VTI";
    E_vertical = Ev;
    E_horizontal = Eh;
    nu_vh = nuvh;
    nu_hh = nuhh;
    G_vertical = Gv;
    
    // Compute stiffness tensor components for VTI
    double Gh = Eh / (2.0 * (1.0 + nuhh));
    
    // C33 from vertical Young's modulus
    double denom = 1.0 - nuhh - 2.0 * nuvh * nuvh * Ev / Eh;
    double C33_val = Ev * (1.0 - nuhh) / denom;
    double C11_val = Eh * (1.0 - nuvh * nuvh * Ev / Eh) / ((1.0 + nuhh) * denom);
    double C13_val = Ev * nuvh / denom;
    double C44_val = Gv;
    double C66_val = Gh;
    
    // Store in Voigt notation (upper triangle)
    Cij[0] = C11_val;   // C11
    Cij[1] = C11_val - 2.0 * C66_val;  // C12
    Cij[2] = C13_val;   // C13
    Cij[6] = C11_val;   // C22 = C11
    Cij[7] = C13_val;   // C23 = C13
    Cij[11] = C33_val;  // C33
    Cij[15] = C44_val;  // C44
    Cij[18] = C44_val;  // C55 = C44
    Cij[20] = C66_val;  // C66
}

void AnisotropicProperties::setThomsen(double Vp0, double Vs0, double eps, double del, double gam, double rho) {
    is_anisotropic = true;
    symmetry_type = "VTI";
    epsilon = eps;
    delta = del;
    gamma = gam;
    
    // Compute stiffness from Thomsen parameters
    double C33_val = rho * Vp0 * Vp0;
    double C44_val = rho * Vs0 * Vs0;
    double C11_val = C33_val * (1.0 + 2.0 * eps);
    double C66_val = C44_val * (1.0 + 2.0 * gam);
    
    // C13 from delta (weak anisotropy approximation)
    double C13_val = std::sqrt((C33_val - C44_val) * (C33_val - C44_val + 2.0 * del * C33_val)) - C44_val;
    
    Cij[0] = C11_val;
    Cij[1] = C11_val - 2.0 * C66_val;
    Cij[2] = C13_val;
    Cij[6] = C11_val;
    Cij[7] = C13_val;
    Cij[11] = C33_val;
    Cij[15] = C44_val;
    Cij[18] = C44_val;
    Cij[20] = C66_val;
    
    // Back-calculate engineering constants (approximate)
    E_vertical = C33_val * (1.0 - (2.0 * C13_val * C13_val) / (C11_val * C33_val));
    E_horizontal = C11_val * (1.0 - (C13_val * C13_val) / (C11_val * C33_val));
    G_vertical = C44_val;
    nu_vh = C13_val / (C11_val + C33_val);
    nu_hh = (C11_val - 2.0 * C66_val) / (2.0 * C11_val);
}

// =============================================================================
// SeismicProperties Implementation
// =============================================================================

SeismicProperties::SeismicProperties()
    : Vp(5000.0)
    , Vs(3000.0)
    , Vp_Vs_ratio(1.67)
    , Qp(100.0)
    , Qs(50.0)
    , acoustic_impedance(13250000.0)
{
}

void SeismicProperties::computeFromElastic(double E, double nu, double rho) {
    double K = computeBulkModulus(E, nu);
    double G = computeShearModulus(E, nu);
    
    Vp = std::sqrt((K + 4.0 * G / 3.0) / rho);
    Vs = std::sqrt(G / rho);
    Vp_Vs_ratio = (Vs > 0.0) ? Vp / Vs : 1.732;
    acoustic_impedance = rho * Vp;
}

void SeismicProperties::computeElastic(double rho, double& E, double& nu) const {
    // From Vp and Vs, compute elastic moduli
    double G = rho * Vs * Vs;
    double K = rho * (Vp * Vp - 4.0 * Vs * Vs / 3.0);
    
    // Convert to E and nu
    E = 9.0 * K * G / (3.0 * K + G);
    nu = (3.0 * K - 2.0 * G) / (2.0 * (3.0 * K + G));
}

// =============================================================================
// MaterialDefinition Implementation
// =============================================================================

RockProperties MaterialDefinition::toRockProperties() const {
    RockProperties props;
    props.name = name;
    
    // Create configuration map
    std::map<std::string, std::string> config;
    
    // Determine mechanical model type
    std::string model_type = "LINEAR_ELASTIC";
    if (recommended_model == "POROELASTIC" || is_porous) {
        model_type = "POROELASTIC";
    } else if (recommended_model == "VISCOELASTIC" || is_viscoelastic) {
        model_type = "VISCOELASTIC_SLS";
    } else if (is_anisotropic_flag) {
        model_type = "ANISOTROPIC_VTI";
    } else if (!recommended_model.empty()) {
        model_type = recommended_model;
    }
    
    config["material_model"] = model_type;
    config["density"] = std::to_string(mechanical.density);
    config["youngs_modulus"] = std::to_string(mechanical.youngs_modulus);
    config["poisson_ratio"] = std::to_string(mechanical.poisson_ratio);
    config["porosity"] = std::to_string(hydraulic.porosity);
    
    // Create mechanical model
    if (model_type == "LINEAR_ELASTIC") {
        auto mat = std::make_unique<LinearElasticMaterial>();
        mat->setElasticModuli(mechanical.youngs_modulus, mechanical.poisson_ratio);
        mat->setDensity(mechanical.density);
        mat->setPorosity(hydraulic.porosity);
        props.mechanical = std::move(mat);
    } else if (model_type == "POROELASTIC") {
        auto mat = std::make_unique<PoroelasticMaterial>();
        mat->setDrainedModuli(mechanical.youngs_modulus, mechanical.poisson_ratio);
        mat->setBiotParameters(mechanical.biot_coefficient, mechanical.biot_modulus);
        mat->setPorosity(hydraulic.porosity);
        mat->setPermeability(hydraulic.permeability);
        props.mechanical = std::move(mat);
    } else if (model_type == "VISCOELASTIC_SLS" || model_type == "VISCOELASTIC") {
        auto mat = std::make_unique<ViscoelasticMaterial>(ViscoelasticMaterial::ViscoType::STANDARD_LINEAR_SOLID);
        mat->setInstantaneousModuli(mechanical.youngs_modulus, mechanical.poisson_ratio);
        mat->setRelaxationTime(mechanical.relaxation_time);
        mat->setLongTermModulus(mechanical.long_term_modulus);
        props.mechanical = std::move(mat);
    } else if (model_type == "ANISOTROPIC_VTI") {
        auto mat = std::make_unique<AnisotropicMaterial>(true);
        mat->setVerticalModuli(anisotropic.E_vertical, anisotropic.G_vertical, anisotropic.nu_vh);
        mat->setHorizontalModuli(anisotropic.E_horizontal, anisotropic.E_horizontal / (2.0 * (1.0 + anisotropic.nu_hh)), anisotropic.nu_hh);
        props.mechanical = std::move(mat);
    } else if (model_type == "ELASTOPLASTIC_MC") {
        auto mat = std::make_unique<ElastoplasticMaterial>(FailureCriterion::MOHR_COULOMB);
        mat->setMohrCoulombParameters(mechanical.cohesion, mechanical.friction_angle * M_PI / 180.0, mechanical.dilation_angle * M_PI / 180.0);
        mat->setTensileStrength(mechanical.tensile_strength);
        props.mechanical = std::move(mat);
    } else {
        // Default to linear elastic
        auto mat = std::make_unique<LinearElasticMaterial>();
        mat->setElasticModuli(mechanical.youngs_modulus, mechanical.poisson_ratio);
        mat->setDensity(mechanical.density);
        mat->setPorosity(hydraulic.porosity);
        props.mechanical = std::move(mat);
    }
    
    // Set up permeability model
    props.permeability = std::make_unique<PermeabilityModel>(PermeabilityModel::PermeabilityType::EXPONENTIAL);
    props.permeability->setInitialPermeability(hydraulic.permeability);
    props.permeability->setStressSensitivity(hydraulic.permeability_stress_coeff);
    props.permeability->setBounds(hydraulic.permeability_min, hydraulic.permeability_max);
    
    // Set thermal properties
    props.thermal.thermal_conductivity = thermal.thermal_conductivity;
    props.thermal.specific_heat = thermal.specific_heat;
    props.thermal.thermal_expansion = thermal.thermal_expansion;
    props.thermal.k_T_coeff = thermal.dK_dT;
    props.thermal.cp_T_coeff = thermal.dCp_dT;
    
    // Set fracture properties
    props.fracture.toughness_KIc = mechanical.fracture_toughness_KIc;
    props.fracture.toughness_KIIc = mechanical.fracture_toughness_KIIc;
    props.fracture.fracture_energy = mechanical.fracture_energy;
    props.fracture.tensile_strength = mechanical.tensile_strength;
    
    return props;
}

std::map<std::string, std::string> MaterialDefinition::toConfigMap() const {
    std::map<std::string, std::string> config;
    
    // Metadata
    config["material.name"] = name;
    config["material.category"] = category;
    config["material.subcategory"] = subcategory;
    
    // Mechanical properties
    config["mechanical.density"] = std::to_string(mechanical.density);
    config["mechanical.youngs_modulus"] = std::to_string(mechanical.youngs_modulus);
    config["mechanical.poisson_ratio"] = std::to_string(mechanical.poisson_ratio);
    config["mechanical.bulk_modulus"] = std::to_string(mechanical.bulk_modulus);
    config["mechanical.shear_modulus"] = std::to_string(mechanical.shear_modulus);
    config["mechanical.tensile_strength"] = std::to_string(mechanical.tensile_strength);
    config["mechanical.compressive_strength"] = std::to_string(mechanical.compressive_strength);
    config["mechanical.cohesion"] = std::to_string(mechanical.cohesion);
    config["mechanical.friction_angle"] = std::to_string(mechanical.friction_angle);
    config["mechanical.biot_coefficient"] = std::to_string(mechanical.biot_coefficient);
    config["mechanical.fracture_toughness"] = std::to_string(mechanical.fracture_toughness_KIc);
    
    // Hydraulic properties
    config["hydraulic.porosity"] = std::to_string(hydraulic.porosity);
    config["hydraulic.permeability"] = std::to_string(hydraulic.permeability);
    config["hydraulic.compressibility"] = std::to_string(hydraulic.compressibility);
    
    // Thermal properties
    config["thermal.conductivity"] = std::to_string(thermal.thermal_conductivity);
    config["thermal.specific_heat"] = std::to_string(thermal.specific_heat);
    config["thermal.expansion"] = std::to_string(thermal.thermal_expansion);
    
    // Seismic properties
    config["seismic.vp"] = std::to_string(seismic.Vp);
    config["seismic.vs"] = std::to_string(seismic.Vs);
    config["seismic.qp"] = std::to_string(seismic.Qp);
    config["seismic.qs"] = std::to_string(seismic.Qs);
    
    return config;
}

// =============================================================================
// MaterialLibrary Implementation
// =============================================================================

MaterialLibrary& MaterialLibrary::getInstance() {
    static MaterialLibrary instance;
    return instance;
}

MaterialLibrary::MaterialLibrary() {
    initializeLibrary();
}

std::string MaterialLibrary::toLower(const std::string& s) const {
    std::string result = s;
    std::transform(result.begin(), result.end(), result.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return result;
}

void MaterialLibrary::initializeLibrary() {
    initializeIgneousRocks();
    initializeSedimentaryRocks();
    initializeMetamorphicRocks();
    initializeUnconsolidatedMaterials();
    initializeMinerals();
    initializeFaultMaterials();
    initializeWellboreMaterials();
    initializeReferenceMaterials();
    initializeFormations();
}

std::optional<MaterialDefinition> MaterialLibrary::getMaterial(const std::string& name) const {
    std::string lower_name = toLower(name);
    
    // Try exact match first
    auto it = materials_.find(lower_name);
    if (it != materials_.end()) {
        return it->second;
    }
    
    // Try partial match
    for (const auto& pair : materials_) {
        if (pair.first.find(lower_name) != std::string::npos ||
            lower_name.find(pair.first) != std::string::npos) {
            return pair.second;
        }
    }
    
    return std::nullopt;
}

std::optional<FormationDefinition> MaterialLibrary::getFormation(const std::string& name) const {
    std::string lower_name = toLower(name);
    
    auto it = formations_.find(lower_name);
    if (it != formations_.end()) {
        return it->second;
    }
    
    // Try partial match
    for (const auto& pair : formations_) {
        if (pair.first.find(lower_name) != std::string::npos ||
            lower_name.find(pair.first) != std::string::npos) {
            return pair.second;
        }
    }
    
    return std::nullopt;
}

std::vector<MaterialDefinition> MaterialLibrary::getMaterialsByCategory(const std::string& category) const {
    std::vector<MaterialDefinition> results;
    std::string lower_cat = toLower(category);
    
    for (const auto& pair : materials_) {
        if (toLower(pair.second.category) == lower_cat) {
            results.push_back(pair.second);
        }
    }
    
    return results;
}

std::vector<std::string> MaterialLibrary::getFormationNames() const {
    std::vector<std::string> names;
    for (const auto& pair : formations_) {
        names.push_back(pair.second.name);
    }
    return names;
}

std::vector<std::string> MaterialLibrary::getMaterialNames() const {
    std::vector<std::string> names;
    for (const auto& pair : materials_) {
        names.push_back(pair.second.name);
    }
    return names;
}

std::vector<MaterialDefinition> MaterialLibrary::searchByYoungsModulus(double E_min, double E_max) const {
    std::vector<MaterialDefinition> results;
    for (const auto& pair : materials_) {
        double E = pair.second.mechanical.youngs_modulus;
        if (E >= E_min && E <= E_max) {
            results.push_back(pair.second);
        }
    }
    return results;
}

std::vector<MaterialDefinition> MaterialLibrary::searchByDensity(double rho_min, double rho_max) const {
    std::vector<MaterialDefinition> results;
    for (const auto& pair : materials_) {
        double rho = pair.second.mechanical.density;
        if (rho >= rho_min && rho <= rho_max) {
            results.push_back(pair.second);
        }
    }
    return results;
}

std::vector<MaterialDefinition> MaterialLibrary::searchByPorosity(double phi_min, double phi_max) const {
    std::vector<MaterialDefinition> results;
    for (const auto& pair : materials_) {
        double phi = pair.second.hydraulic.porosity;
        if (phi >= phi_min && phi <= phi_max) {
            results.push_back(pair.second);
        }
    }
    return results;
}

std::vector<MaterialDefinition> MaterialLibrary::searchByPermeability(double k_min, double k_max) const {
    std::vector<MaterialDefinition> results;
    for (const auto& pair : materials_) {
        double k = pair.second.hydraulic.permeability;
        if (k >= k_min && k <= k_max) {
            results.push_back(pair.second);
        }
    }
    return results;
}

std::vector<std::string> MaterialLibrary::getCategories() const {
    std::vector<std::string> categories;
    std::map<std::string, bool> seen;
    
    for (const auto& pair : materials_) {
        if (!seen[pair.second.category]) {
            categories.push_back(pair.second.category);
            seen[pair.second.category] = true;
        }
    }
    
    return categories;
}

std::vector<std::string> MaterialLibrary::getSubcategories(const std::string& category) const {
    std::vector<std::string> subcategories;
    std::map<std::string, bool> seen;
    std::string lower_cat = toLower(category);
    
    for (const auto& pair : materials_) {
        if (toLower(pair.second.category) == lower_cat && !seen[pair.second.subcategory]) {
            subcategories.push_back(pair.second.subcategory);
            seen[pair.second.subcategory] = true;
        }
    }
    
    return subcategories;
}

MaterialDefinition MaterialLibrary::createVariant(const std::string& base_name,
                                                   const std::map<std::string, double>& modifications) const {
    auto base_opt = getMaterial(base_name);
    if (!base_opt) {
        throw std::runtime_error("Base material not found: " + base_name);
    }
    
    MaterialDefinition variant = *base_opt;
    variant.name = base_name + "_variant";
    
    for (const auto& mod : modifications) {
        const std::string& key = mod.first;
        double value = mod.second;
        
        if (key == "density") variant.mechanical.density = value;
        else if (key == "youngs_modulus") variant.mechanical.youngs_modulus = value;
        else if (key == "poisson_ratio") variant.mechanical.poisson_ratio = value;
        else if (key == "porosity") variant.hydraulic.porosity = value;
        else if (key == "permeability") variant.hydraulic.setIsotropicPermeability(value);
        else if (key == "thermal_conductivity") variant.thermal.thermal_conductivity = value;
        else if (key == "cohesion") variant.mechanical.cohesion = value;
        else if (key == "friction_angle") variant.mechanical.friction_angle = value;
        else if (key == "tensile_strength") variant.mechanical.tensile_strength = value;
        else if (key == "compressive_strength") variant.mechanical.compressive_strength = value;
    }
    
    variant.mechanical.computeDerivedProperties();
    variant.seismic.computeFromElastic(variant.mechanical.youngs_modulus,
                                        variant.mechanical.poisson_ratio,
                                        variant.mechanical.density);
    
    return variant;
}

RockProperties MaterialLibrary::createRockProperties(const std::string& material_name) const {
    auto mat_opt = getMaterial(material_name);
    if (!mat_opt) {
        throw std::runtime_error("Material not found: " + material_name);
    }
    return mat_opt->toRockProperties();
}

std::string MaterialLibrary::exportToConfig(const std::string& name) const {
    auto mat_opt = getMaterial(name);
    if (!mat_opt) {
        return "";
    }
    return exportToConfig(*mat_opt);
}

std::string MaterialLibrary::exportToConfig(const MaterialDefinition& mat) const {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6);
    
    oss << "# Material: " << mat.name << "\n";
    oss << "# Category: " << mat.category << " / " << mat.subcategory << "\n";
    oss << "# " << mat.description << "\n";
    if (!mat.reference.empty()) {
        oss << "# Reference: " << mat.reference << "\n";
    }
    oss << "\n";
    
    oss << "[material]\n";
    oss << "name = " << mat.name << "\n";
    oss << "model = " << mat.recommended_model << "\n";
    oss << "\n";
    
    oss << "[mechanical]\n";
    oss << "density = " << mat.mechanical.density << "\n";
    oss << "youngs_modulus = " << mat.mechanical.youngs_modulus << "\n";
    oss << "poisson_ratio = " << mat.mechanical.poisson_ratio << "\n";
    oss << "tensile_strength = " << mat.mechanical.tensile_strength << "\n";
    oss << "compressive_strength = " << mat.mechanical.compressive_strength << "\n";
    oss << "cohesion = " << mat.mechanical.cohesion << "\n";
    oss << "friction_angle = " << mat.mechanical.friction_angle << "\n";
    oss << "biot_coefficient = " << mat.mechanical.biot_coefficient << "\n";
    oss << "fracture_toughness_KIc = " << mat.mechanical.fracture_toughness_KIc << "\n";
    oss << "\n";
    
    oss << "[hydraulic]\n";
    oss << "porosity = " << mat.hydraulic.porosity << "\n";
    oss << "permeability = " << mat.hydraulic.permeability << "\n";
    oss << "compressibility = " << mat.hydraulic.compressibility << "\n";
    oss << "\n";
    
    oss << "[thermal]\n";
    oss << "conductivity = " << mat.thermal.thermal_conductivity << "\n";
    oss << "specific_heat = " << mat.thermal.specific_heat << "\n";
    oss << "thermal_expansion = " << mat.thermal.thermal_expansion << "\n";
    oss << "\n";
    
    oss << "[seismic]\n";
    oss << "vp = " << mat.seismic.Vp << "\n";
    oss << "vs = " << mat.seismic.Vs << "\n";
    oss << "qp = " << mat.seismic.Qp << "\n";
    oss << "qs = " << mat.seismic.Qs << "\n";
    
    return oss.str();
}

MaterialDefinition MaterialLibrary::interpolate(const MaterialDefinition& mat1,
                                                 const MaterialDefinition& mat2,
                                                 double fraction) const {
    MaterialDefinition result;
    double f1 = 1.0 - fraction;
    double f2 = fraction;
    
    result.name = mat1.name + "_" + mat2.name + "_mix";
    result.category = "interpolated";
    result.subcategory = "mixture";
    result.description = "Interpolation of " + mat1.name + " and " + mat2.name;
    
    // Interpolate mechanical properties
    result.mechanical.density = f1 * mat1.mechanical.density + f2 * mat2.mechanical.density;
    result.mechanical.youngs_modulus = f1 * mat1.mechanical.youngs_modulus + f2 * mat2.mechanical.youngs_modulus;
    result.mechanical.poisson_ratio = f1 * mat1.mechanical.poisson_ratio + f2 * mat2.mechanical.poisson_ratio;
    result.mechanical.tensile_strength = f1 * mat1.mechanical.tensile_strength + f2 * mat2.mechanical.tensile_strength;
    result.mechanical.compressive_strength = f1 * mat1.mechanical.compressive_strength + f2 * mat2.mechanical.compressive_strength;
    result.mechanical.cohesion = f1 * mat1.mechanical.cohesion + f2 * mat2.mechanical.cohesion;
    result.mechanical.friction_angle = f1 * mat1.mechanical.friction_angle + f2 * mat2.mechanical.friction_angle;
    result.mechanical.biot_coefficient = f1 * mat1.mechanical.biot_coefficient + f2 * mat2.mechanical.biot_coefficient;
    result.mechanical.fracture_toughness_KIc = f1 * mat1.mechanical.fracture_toughness_KIc + f2 * mat2.mechanical.fracture_toughness_KIc;
    
    // Interpolate hydraulic properties (geometric mean for permeability)
    result.hydraulic.porosity = f1 * mat1.hydraulic.porosity + f2 * mat2.hydraulic.porosity;
    result.hydraulic.permeability = std::pow(mat1.hydraulic.permeability, f1) * std::pow(mat2.hydraulic.permeability, f2);
    
    // Interpolate thermal properties
    result.thermal.thermal_conductivity = f1 * mat1.thermal.thermal_conductivity + f2 * mat2.thermal.thermal_conductivity;
    result.thermal.specific_heat = f1 * mat1.thermal.specific_heat + f2 * mat2.thermal.specific_heat;
    result.thermal.thermal_expansion = f1 * mat1.thermal.thermal_expansion + f2 * mat2.thermal.thermal_expansion;
    
    result.mechanical.computeDerivedProperties();
    result.seismic.computeFromElastic(result.mechanical.youngs_modulus,
                                       result.mechanical.poisson_ratio,
                                       result.mechanical.density);
    
    result.is_porous = mat1.is_porous || mat2.is_porous;
    result.is_fractured = mat1.is_fractured || mat2.is_fractured;
    result.is_viscoelastic = mat1.is_viscoelastic || mat2.is_viscoelastic;
    result.is_anisotropic_flag = mat1.is_anisotropic_flag || mat2.is_anisotropic_flag;
    result.recommended_model = "LINEAR_ELASTIC";
    
    return result;
}

MaterialDefinition MaterialLibrary::getDepthAdjusted(const std::string& name, double depth) const {
    auto mat_opt = getMaterial(name);
    if (!mat_opt) {
        throw std::runtime_error("Material not found: " + name);
    }
    
    MaterialDefinition mat = *mat_opt;
    mat.name = name + "_depth_" + std::to_string(static_cast<int>(depth)) + "m";
    
    // Apply depth-dependent property modifications
    // Typical gradients based on geomechanics literature
    
    // Pressure at depth (hydrostatic + lithostatic effect)
    double pressure = depth * 2650.0 * 9.81;  // Approximate lithostatic
    
    // Young's modulus increases with confining pressure
    if (mat.mechanical.dE_dP != 0.0) {
        mat.mechanical.youngs_modulus += mat.mechanical.dE_dP * pressure;
    } else {
        // Default increase of ~1% per 100m depth
        mat.mechanical.youngs_modulus *= (1.0 + 0.0001 * depth);
    }
    
    // Porosity decreases with depth (compaction)
    double porosity_reduction = 1.0 - std::exp(-depth / 3000.0);  // e-folding depth of 3km
    mat.hydraulic.porosity *= (1.0 - 0.5 * porosity_reduction);
    mat.hydraulic.porosity = std::max(0.001, mat.hydraulic.porosity);
    
    // Permeability decreases exponentially
    mat.hydraulic.permeability *= std::exp(-mat.hydraulic.permeability_stress_coeff * pressure);
    mat.hydraulic.permeability = std::max(1e-22, mat.hydraulic.permeability);
    
    // Temperature effect (assuming 30°C/km gradient)
    double temperature_increase = depth * 0.030;  // °C
    if (mat.mechanical.dE_dT != 0.0) {
        mat.mechanical.youngs_modulus += mat.mechanical.dE_dT * temperature_increase;
    }
    
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus,
                                    mat.mechanical.poisson_ratio,
                                    mat.mechanical.density);
    
    return mat;
}

void MaterialLibrary::addMaterial(const MaterialDefinition& material) {
    std::string key = toLower(material.name);
    materials_[key] = material;
}

void MaterialLibrary::addFormation(const FormationDefinition& formation) {
    std::string key = toLower(formation.name);
    formations_[key] = formation;
}

// =============================================================================
// Igneous Rock Initialization
// =============================================================================

void MaterialLibrary::initializeIgneousRocks() {
    addMaterial(Materials::Igneous::Granite());
    addMaterial(Materials::Igneous::GraniteWeathered());
    addMaterial(Materials::Igneous::Granodiorite());
    addMaterial(Materials::Igneous::Diorite());
    addMaterial(Materials::Igneous::Gabbro());
    addMaterial(Materials::Igneous::Basalt());
    addMaterial(Materials::Igneous::BasaltVesicular());
    addMaterial(Materials::Igneous::Andesite());
    addMaterial(Materials::Igneous::Rhyolite());
    addMaterial(Materials::Igneous::Obsidian());
    addMaterial(Materials::Igneous::Pumice());
    addMaterial(Materials::Igneous::Tuff());
    addMaterial(Materials::Igneous::TuffWelded());
    addMaterial(Materials::Igneous::Dacite());
    addMaterial(Materials::Igneous::Peridotite());
    addMaterial(Materials::Igneous::Dunite());
    addMaterial(Materials::Igneous::Syenite());
    addMaterial(Materials::Igneous::Diabase());
    addMaterial(Materials::Igneous::Pegmatite());
    addMaterial(Materials::Igneous::Porphyry());
    addMaterial(Materials::Igneous::Aplite());
    addMaterial(Materials::Igneous::Komatiite());
    addMaterial(Materials::Igneous::Phonolite());
    addMaterial(Materials::Igneous::Latite());
    addMaterial(Materials::Igneous::Monzonite());
    addMaterial(Materials::Igneous::Troctolite());
    addMaterial(Materials::Igneous::Norite());
    addMaterial(Materials::Igneous::Anorthosite());
    addMaterial(Materials::Igneous::Kimberlite());
    addMaterial(Materials::Igneous::Lamprophyre());
}

void MaterialLibrary::initializeSedimentaryRocks() {
    // Clastic rocks
    addMaterial(Materials::Sedimentary::Sandstone());
    addMaterial(Materials::Sedimentary::SandstoneQuartz());
    addMaterial(Materials::Sedimentary::SandstoneFeldspathic());
    addMaterial(Materials::Sedimentary::SandstoneLithic());
    addMaterial(Materials::Sedimentary::SandstoneArkose());
    addMaterial(Materials::Sedimentary::SandstoneGreywacke());
    addMaterial(Materials::Sedimentary::SandstoneTight());
    addMaterial(Materials::Sedimentary::Siltstone());
    addMaterial(Materials::Sedimentary::Shale());
    addMaterial(Materials::Sedimentary::ShaleOrganic());
    addMaterial(Materials::Sedimentary::ShaleSiliceous());
    addMaterial(Materials::Sedimentary::ShaleCalcareous());
    addMaterial(Materials::Sedimentary::Mudstone());
    addMaterial(Materials::Sedimentary::Claystone());
    addMaterial(Materials::Sedimentary::Conglomerate());
    addMaterial(Materials::Sedimentary::Breccia());
    addMaterial(Materials::Sedimentary::Tillite());
    addMaterial(Materials::Sedimentary::Loess());
    addMaterial(Materials::Sedimentary::Turbidite());
    
    // Carbonate rocks
    addMaterial(Materials::Sedimentary::Limestone());
    addMaterial(Materials::Sedimentary::LimestoneOolitic());
    addMaterial(Materials::Sedimentary::LimestoneMicritic());
    addMaterial(Materials::Sedimentary::LimestoneBioclastic());
    addMaterial(Materials::Sedimentary::LimestoneChalk());
    addMaterial(Materials::Sedimentary::LimestoneReef());
    addMaterial(Materials::Sedimentary::LimestoneTravertine());
    addMaterial(Materials::Sedimentary::Dolomite());
    addMaterial(Materials::Sedimentary::DolomiteSucrosic());
    addMaterial(Materials::Sedimentary::Marl());
    addMaterial(Materials::Sedimentary::Coquina());
    
    // Evaporites and chemical rocks
    addMaterial(Materials::Sedimentary::Halite());
    addMaterial(Materials::Sedimentary::Gypsum());
    addMaterial(Materials::Sedimentary::Anhydrite());
    addMaterial(Materials::Sedimentary::Sylvite());
    addMaterial(Materials::Sedimentary::Potash());
    addMaterial(Materials::Sedimentary::Trona());
    addMaterial(Materials::Sedimentary::Chert());
    addMaterial(Materials::Sedimentary::Flint());
    addMaterial(Materials::Sedimentary::Diatomite());
    addMaterial(Materials::Sedimentary::Ironstone());
    addMaterial(Materials::Sedimentary::BandedIronFormation());
    addMaterial(Materials::Sedimentary::Phosphorite());
    
    // Organic rocks
    addMaterial(Materials::Sedimentary::Coal());
    addMaterial(Materials::Sedimentary::CoalAnthracite());
    addMaterial(Materials::Sedimentary::CoalBituminous());
    addMaterial(Materials::Sedimentary::CoalLignite());
    addMaterial(Materials::Sedimentary::Peat());
    addMaterial(Materials::Sedimentary::OilShale());
}

void MaterialLibrary::initializeMetamorphicRocks() {
    // Foliated
    addMaterial(Materials::Metamorphic::Slate());
    addMaterial(Materials::Metamorphic::Phyllite());
    addMaterial(Materials::Metamorphic::Schist());
    addMaterial(Materials::Metamorphic::SchistMica());
    addMaterial(Materials::Metamorphic::SchistChlorite());
    addMaterial(Materials::Metamorphic::SchistTalc());
    addMaterial(Materials::Metamorphic::SchistBlueschist());
    addMaterial(Materials::Metamorphic::Gneiss());
    addMaterial(Materials::Metamorphic::GneissGranitic());
    addMaterial(Materials::Metamorphic::GneissBanded());
    addMaterial(Materials::Metamorphic::Migmatite());
    addMaterial(Materials::Metamorphic::Mylonite());
    addMaterial(Materials::Metamorphic::Cataclasite());
    addMaterial(Materials::Metamorphic::Amphibolite());
    
    // Non-foliated
    addMaterial(Materials::Metamorphic::Marble());
    addMaterial(Materials::Metamorphic::MarbleCalcite());
    addMaterial(Materials::Metamorphic::MarbleDolomitic());
    addMaterial(Materials::Metamorphic::Quartzite());
    addMaterial(Materials::Metamorphic::Hornfels());
    addMaterial(Materials::Metamorphic::Serpentinite());
    addMaterial(Materials::Metamorphic::Soapstone());
    addMaterial(Materials::Metamorphic::Eclogite());
    addMaterial(Materials::Metamorphic::Granulite());
    addMaterial(Materials::Metamorphic::Greenstone());
    addMaterial(Materials::Metamorphic::Skarn());
    addMaterial(Materials::Metamorphic::Tactite());
}

void MaterialLibrary::initializeUnconsolidatedMaterials() {
    addMaterial(Materials::Unconsolidated::Sand());
    addMaterial(Materials::Unconsolidated::SandFine());
    addMaterial(Materials::Unconsolidated::SandMedium());
    addMaterial(Materials::Unconsolidated::SandCoarse());
    addMaterial(Materials::Unconsolidated::SandSilty());
    addMaterial(Materials::Unconsolidated::SandClayey());
    addMaterial(Materials::Unconsolidated::Gravel());
    addMaterial(Materials::Unconsolidated::GravelSandy());
    addMaterial(Materials::Unconsolidated::Clay());
    addMaterial(Materials::Unconsolidated::ClayKaolinite());
    addMaterial(Materials::Unconsolidated::ClayIllite());
    addMaterial(Materials::Unconsolidated::ClayMontmorillonite());
    addMaterial(Materials::Unconsolidated::ClayBentonite());
    addMaterial(Materials::Unconsolidated::Silt());
    addMaterial(Materials::Unconsolidated::Loam());
    addMaterial(Materials::Unconsolidated::Till());
    addMaterial(Materials::Unconsolidated::Alluvium());
    addMaterial(Materials::Unconsolidated::ColluviumTalus());
    addMaterial(Materials::Unconsolidated::Laterite());
    addMaterial(Materials::Unconsolidated::Saprolite());
    addMaterial(Materials::Unconsolidated::Regolith());
}

void MaterialLibrary::initializeMinerals() {
    // Silicates
    addMaterial(Materials::Minerals::Quartz());
    addMaterial(Materials::Minerals::Feldspar());
    addMaterial(Materials::Minerals::FeldsparOrthoclase());
    addMaterial(Materials::Minerals::FeldsparPlagioclase());
    addMaterial(Materials::Minerals::Mica());
    addMaterial(Materials::Minerals::MicaMuscovite());
    addMaterial(Materials::Minerals::MicaBiotite());
    addMaterial(Materials::Minerals::Olivine());
    addMaterial(Materials::Minerals::Pyroxene());
    addMaterial(Materials::Minerals::Amphibole());
    addMaterial(Materials::Minerals::Hornblende());
    addMaterial(Materials::Minerals::Garnet());
    addMaterial(Materials::Minerals::Kyanite());
    addMaterial(Materials::Minerals::Sillimanite());
    addMaterial(Materials::Minerals::Andalusite());
    addMaterial(Materials::Minerals::Tourmaline());
    addMaterial(Materials::Minerals::Zircon());
    addMaterial(Materials::Minerals::Epidote());
    addMaterial(Materials::Minerals::Chlorite());
    addMaterial(Materials::Minerals::Serpentine());
    addMaterial(Materials::Minerals::Talc());
    addMaterial(Materials::Minerals::Zeolites());
    
    // Carbonates
    addMaterial(Materials::Minerals::Calcite());
    addMaterial(Materials::Minerals::Aragonite());
    addMaterial(Materials::Minerals::DolomiteMinite());
    addMaterial(Materials::Minerals::Siderite());
    addMaterial(Materials::Minerals::Magneite());
    addMaterial(Materials::Minerals::Rhodochrosite());
    
    // Sulfates
    addMaterial(Materials::Minerals::GypsumMineral());
    addMaterial(Materials::Minerals::AnhydriteMineral());
    addMaterial(Materials::Minerals::Barite());
    addMaterial(Materials::Minerals::Celestite());
    
    // Halides
    addMaterial(Materials::Minerals::HaliteMineral());
    addMaterial(Materials::Minerals::Fluorite());
    addMaterial(Materials::Minerals::SylviteMineral());
    
    // Oxides
    addMaterial(Materials::Minerals::Magnetite());
    addMaterial(Materials::Minerals::Hematite());
    addMaterial(Materials::Minerals::Ilmenite());
    addMaterial(Materials::Minerals::Rutile());
    addMaterial(Materials::Minerals::Corundum());
    addMaterial(Materials::Minerals::Spinel());
    addMaterial(Materials::Minerals::Chromite());
    addMaterial(Materials::Minerals::Limonite());
    addMaterial(Materials::Minerals::Goethite());
    
    // Sulfides
    addMaterial(Materials::Minerals::Pyrite());
    addMaterial(Materials::Minerals::Pyrrhotite());
    addMaterial(Materials::Minerals::Galena());
    addMaterial(Materials::Minerals::Sphalerite());
    addMaterial(Materials::Minerals::Chalcopyrite());
    addMaterial(Materials::Minerals::Molybdenite());
    
    // Native elements
    addMaterial(Materials::Minerals::GraphiteMineral());
    addMaterial(Materials::Minerals::SulfurNative());
    addMaterial(Materials::Minerals::GoldNative());
    addMaterial(Materials::Minerals::CopperNative());
    
    // Clays
    addMaterial(Materials::Minerals::Kaolinite());
    addMaterial(Materials::Minerals::Illite());
    addMaterial(Materials::Minerals::Montmorillonite());
    addMaterial(Materials::Minerals::Smectite());
    addMaterial(Materials::Minerals::Vermiculite());
    addMaterial(Materials::Minerals::Palygorskite());
}

void MaterialLibrary::initializeFaultMaterials() {
    addMaterial(Materials::FaultMaterials::FaultGouge());
    addMaterial(Materials::FaultMaterials::FaultGougeClayRich());
    addMaterial(Materials::FaultMaterials::FaultBreccia());
    addMaterial(Materials::FaultMaterials::FaultCataclasite());
    addMaterial(Materials::FaultMaterials::FaultMylonite());
    addMaterial(Materials::FaultMaterials::FaultPseudotachylyte());
    addMaterial(Materials::FaultMaterials::FaultDamageZone());
    addMaterial(Materials::FaultMaterials::FaultCore());
    addMaterial(Materials::FaultMaterials::SlipZone());
    addMaterial(Materials::FaultMaterials::ShearedZone());
}

void MaterialLibrary::initializeWellboreMaterials() {
    addMaterial(Materials::Wellbore::CementClassA());
    addMaterial(Materials::Wellbore::CementClassC());
    addMaterial(Materials::Wellbore::CementClassG());
    addMaterial(Materials::Wellbore::CementClassH());
    addMaterial(Materials::Wellbore::CementFoamed());
    addMaterial(Materials::Wellbore::CementLightweight());
    addMaterial(Materials::Wellbore::CementHighDensity());
    addMaterial(Materials::Wellbore::SteelCasing());
    addMaterial(Materials::Wellbore::SteelTubing());
    addMaterial(Materials::Wellbore::SteelDrillPipe());
    addMaterial(Materials::Wellbore::Proppant());
    addMaterial(Materials::Wellbore::ProppantSand());
    addMaterial(Materials::Wellbore::ProppantCeramic());
    addMaterial(Materials::Wellbore::ProppantResinCoated());
    addMaterial(Materials::Wellbore::DrillCuttings());
}

void MaterialLibrary::initializeReferenceMaterials() {
    addMaterial(Materials::Reference::Steel());
    addMaterial(Materials::Reference::SteelCarbon());
    addMaterial(Materials::Reference::SteelStainless());
    addMaterial(Materials::Reference::Aluminum());
    addMaterial(Materials::Reference::Concrete());
    addMaterial(Materials::Reference::ConcreteHighStrength());
    addMaterial(Materials::Reference::Glass());
    addMaterial(Materials::Reference::Acrylic());
    addMaterial(Materials::Reference::Copper());
    addMaterial(Materials::Reference::Titanium());
    addMaterial(Materials::Reference::InconelAlloy());
    addMaterial(Materials::Reference::Water());
    addMaterial(Materials::Reference::Ice());
    addMaterial(Materials::Reference::Air());
    addMaterial(Materials::Reference::Rubber());
    addMaterial(Materials::Reference::Epoxy());
    addMaterial(Materials::Reference::CarbonFiber());
    addMaterial(Materials::Reference::Ceramic());
}

void MaterialLibrary::initializeFormations() {
    // Formations will use the base material definitions with formation-specific overrides
    // For now, we add a representative subset of key formations
    
    FormationDefinition bakken;
    bakken.name = "Bakken_Shale";
    bakken.basin = "Williston";
    bakken.region = "North Dakota/Montana";
    bakken.age = "Late Devonian-Early Mississippian";
    bakken.lithology = "Organic shale/dolomite/sandstone";
    bakken.depth_min = 2500.0;
    bakken.depth_max = 3500.0;
    bakken.material = Materials::Sedimentary::ShaleOrganic();
    bakken.material.name = "Bakken_Shale";
    bakken.Sv_gradient = 23.0;
    bakken.SHmax_gradient = 20.0;
    bakken.Shmin_gradient = 15.0;
    bakken.Pp_gradient = 10.5;
    bakken.temperature_gradient = 30.0;
    bakken.net_to_gross = 0.15;
    bakken.water_saturation = 0.30;
    bakken.oil_saturation = 0.50;
    bakken.gas_saturation = 0.20;
    bakken.notes = "Three members: Upper Bakken shale, Middle Bakken dolomite/sandstone, Lower Bakken shale";
    addFormation(bakken);
    
    FormationDefinition eagleford;
    eagleford.name = "Eagle_Ford_Shale";
    eagleford.basin = "Western Gulf";
    eagleford.region = "South Texas";
    eagleford.age = "Late Cretaceous";
    eagleford.lithology = "Organic-rich calcareous mudstone";
    eagleford.depth_min = 1200.0;
    eagleford.depth_max = 4300.0;
    eagleford.material = Materials::Sedimentary::ShaleCalcareous();
    eagleford.material.name = "Eagle_Ford_Shale";
    eagleford.Sv_gradient = 22.6;
    eagleford.SHmax_gradient = 19.0;
    eagleford.Shmin_gradient = 14.0;
    eagleford.Pp_gradient = 10.0;
    eagleford.temperature_gradient = 28.0;
    eagleford.net_to_gross = 0.50;
    eagleford.water_saturation = 0.25;
    eagleford.oil_saturation = 0.45;
    eagleford.gas_saturation = 0.30;
    addFormation(eagleford);
    
    FormationDefinition marcellus;
    marcellus.name = "Marcellus_Shale";
    marcellus.basin = "Appalachian";
    marcellus.region = "Pennsylvania/West Virginia/Ohio";
    marcellus.age = "Middle Devonian";
    marcellus.lithology = "Black organic shale";
    marcellus.depth_min = 1200.0;
    marcellus.depth_max = 2700.0;
    marcellus.material = Materials::Sedimentary::ShaleOrganic();
    marcellus.material.name = "Marcellus_Shale";
    marcellus.Sv_gradient = 24.0;
    marcellus.SHmax_gradient = 28.0;
    marcellus.Shmin_gradient = 18.0;
    marcellus.Pp_gradient = 9.8;
    marcellus.temperature_gradient = 25.0;
    marcellus.net_to_gross = 0.60;
    marcellus.water_saturation = 0.20;
    marcellus.oil_saturation = 0.05;
    marcellus.gas_saturation = 0.75;
    addFormation(marcellus);
    
    FormationDefinition permian;
    permian.name = "Permian_Basin_Wolfcamp";
    permian.basin = "Permian";
    permian.region = "West Texas/New Mexico";
    permian.age = "Permian";
    permian.lithology = "Carbonate/shale/sandstone";
    permian.depth_min = 2000.0;
    permian.depth_max = 4000.0;
    permian.material = Materials::Sedimentary::ShaleSiliceous();
    permian.material.name = "Wolfcamp";
    permian.Sv_gradient = 23.5;
    permian.SHmax_gradient = 19.5;
    permian.Shmin_gradient = 14.5;
    permian.Pp_gradient = 10.0;
    permian.temperature_gradient = 28.0;
    permian.net_to_gross = 0.35;
    permian.water_saturation = 0.30;
    permian.oil_saturation = 0.50;
    permian.gas_saturation = 0.20;
    addFormation(permian);
    
    FormationDefinition northsea;
    northsea.name = "North_Sea_Chalk";
    northsea.basin = "North Sea";
    northsea.region = "Norwegian/Danish Sector";
    northsea.age = "Late Cretaceous";
    northsea.lithology = "Chalk";
    northsea.depth_min = 2000.0;
    northsea.depth_max = 3500.0;
    northsea.material = Materials::Sedimentary::LimestoneChalk();
    northsea.material.name = "North_Sea_Chalk";
    northsea.Sv_gradient = 22.0;
    northsea.SHmax_gradient = 18.0;
    northsea.Shmin_gradient = 15.0;
    northsea.Pp_gradient = 10.0;
    northsea.temperature_gradient = 35.0;
    northsea.net_to_gross = 0.70;
    northsea.water_saturation = 0.30;
    northsea.oil_saturation = 0.60;
    northsea.gas_saturation = 0.10;
    addFormation(northsea);
    
    FormationDefinition ghawar;
    ghawar.name = "Ghawar_Arab_D";
    ghawar.basin = "Arabian";
    ghawar.region = "Saudi Arabia";
    ghawar.age = "Jurassic";
    ghawar.lithology = "Limestone/dolomite";
    ghawar.depth_min = 1800.0;
    ghawar.depth_max = 2200.0;
    ghawar.material = Materials::Sedimentary::LimestoneOolitic();
    ghawar.material.name = "Arab_D";
    ghawar.Sv_gradient = 22.5;
    ghawar.SHmax_gradient = 18.0;
    ghawar.Shmin_gradient = 14.0;
    ghawar.Pp_gradient = 9.5;
    ghawar.temperature_gradient = 32.0;
    ghawar.net_to_gross = 0.85;
    ghawar.water_saturation = 0.15;
    ghawar.oil_saturation = 0.80;
    ghawar.gas_saturation = 0.05;
    addFormation(ghawar);
    
    FormationDefinition geysers;
    geysers.name = "Geysers_Geothermal";
    geysers.basin = "Clear Lake";
    geysers.region = "Northern California";
    geysers.age = "Mesozoic-Cenozoic";
    geysers.lithology = "Graywacke/metavolcanics";
    geysers.depth_min = 1000.0;
    geysers.depth_max = 3000.0;
    geysers.material = Materials::Sedimentary::SandstoneGreywacke();
    geysers.material.name = "Geysers_Reservoir";
    geysers.Sv_gradient = 23.0;
    geysers.SHmax_gradient = 30.0;
    geysers.Shmin_gradient = 20.0;
    geysers.Pp_gradient = 3.0;
    geysers.temperature_gradient = 100.0;
    geysers.net_to_gross = 0.40;
    geysers.water_saturation = 0.20;
    geysers.oil_saturation = 0.0;
    geysers.gas_saturation = 0.80;
    geysers.notes = "Vapor-dominated geothermal system";
    addFormation(geysers);
    
    FormationDefinition sleipner;
    sleipner.name = "Sleipner_Utsira";
    sleipner.basin = "North Sea";
    sleipner.region = "Norwegian Sector";
    sleipner.age = "Miocene-Pliocene";
    sleipner.lithology = "Unconsolidated sand";
    sleipner.depth_min = 800.0;
    sleipner.depth_max = 1100.0;
    sleipner.material = Materials::Sedimentary::Sandstone();
    sleipner.material.name = "Utsira_Sand";
    sleipner.material.hydraulic.porosity = 0.35;
    sleipner.material.hydraulic.setIsotropicPermeability(mDToM2(2000.0));
    sleipner.Sv_gradient = 20.0;
    sleipner.SHmax_gradient = 16.0;
    sleipner.Shmin_gradient = 14.0;
    sleipner.Pp_gradient = 10.0;
    sleipner.temperature_gradient = 35.0;
    sleipner.net_to_gross = 0.95;
    sleipner.water_saturation = 1.0;
    sleipner.oil_saturation = 0.0;
    sleipner.gas_saturation = 0.0;
    sleipner.notes = "Major CO2 storage site since 1996";
    addFormation(sleipner);
    
    FormationDefinition ogallala;
    ogallala.name = "Ogallala_Aquifer";
    ogallala.basin = "High Plains";
    ogallala.region = "Central USA";
    ogallala.age = "Miocene-Pliocene";
    ogallala.lithology = "Sand/gravel/silt";
    ogallala.depth_min = 0.0;
    ogallala.depth_max = 150.0;
    ogallala.material = Materials::Unconsolidated::GravelSandy();
    ogallala.material.name = "Ogallala";
    ogallala.material.hydraulic.porosity = 0.30;
    ogallala.material.hydraulic.setIsotropicPermeability(mDToM2(5000.0));
    ogallala.Sv_gradient = 20.0;
    ogallala.SHmax_gradient = 15.0;
    ogallala.Shmin_gradient = 12.0;
    ogallala.Pp_gradient = 9.8;
    ogallala.temperature_gradient = 25.0;
    ogallala.net_to_gross = 0.80;
    ogallala.water_saturation = 1.0;
    ogallala.oil_saturation = 0.0;
    ogallala.gas_saturation = 0.0;
    ogallala.notes = "Major groundwater source for US agriculture";
    addFormation(ogallala);
    
    FormationDefinition sanandreas;
    sanandreas.name = "San_Andreas_Fault";
    sanandreas.basin = "California";
    sanandreas.region = "California";
    sanandreas.age = "Neogene-present";
    sanandreas.lithology = "Fault gouge/breccia/granite";
    sanandreas.depth_min = 0.0;
    sanandreas.depth_max = 15000.0;
    sanandreas.material = Materials::FaultMaterials::FaultGouge();
    sanandreas.material.name = "SAF_Gouge";
    sanandreas.Sv_gradient = 26.0;
    sanandreas.SHmax_gradient = 35.0;
    sanandreas.Shmin_gradient = 20.0;
    sanandreas.Pp_gradient = 10.0;
    sanandreas.temperature_gradient = 25.0;
    sanandreas.net_to_gross = 1.0;
    sanandreas.water_saturation = 0.50;
    sanandreas.oil_saturation = 0.0;
    sanandreas.gas_saturation = 0.0;
    sanandreas.notes = "Major transform fault, strike-slip motion";
    addFormation(sanandreas);
}

// =============================================================================
// Igneous Rock Definitions - namespace Materials::Igneous
// =============================================================================

namespace Materials {
namespace Igneous {

MaterialDefinition Granite() {
    MaterialDefinition mat;
    mat.name = "Granite";
    mat.category = "Igneous";
    mat.subcategory = "Plutonic";
    mat.description = "Coarse-grained felsic intrusive rock composed of quartz, feldspar, and mica";
    mat.reference = "Jaeger et al., 2007; Carmichael, 1989";
    
    mat.mechanical.density = 2650.0;
    mat.mechanical.youngs_modulus = GPaToPa(60.0);
    mat.mechanical.poisson_ratio = 0.25;
    mat.mechanical.tensile_strength = MPaToPa(10.0);
    mat.mechanical.compressive_strength = MPaToPa(170.0);
    mat.mechanical.cohesion = MPaToPa(30.0);
    mat.mechanical.friction_angle = 55.0;
    mat.mechanical.dilation_angle = 15.0;
    mat.mechanical.fracture_toughness_KIc = 1.5e6;
    mat.mechanical.fracture_toughness_KIIc = 2.0e6;
    mat.mechanical.fracture_energy = 50.0;
    mat.mechanical.biot_coefficient = 0.45;
    mat.mechanical.grain_bulk_modulus = GPaToPa(45.0);
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.01;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.001));
    mat.hydraulic.compressibility = 1.0e-11;
    mat.hydraulic.tortuosity = 3.0;
    
    mat.thermal.thermal_conductivity = 3.0;
    mat.thermal.specific_heat = 790.0;
    mat.thermal.thermal_expansion = 8.0e-6;
    mat.thermal.heat_production = 2.5e-6;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 200.0;
    mat.seismic.Qs = 100.0;
    
    mat.is_porous = false;
    mat.is_fractured = false;
    mat.is_viscoelastic = false;
    mat.is_anisotropic_flag = false;
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition GraniteWeathered() {
    MaterialDefinition mat = Granite();
    mat.name = "Granite_Weathered";
    mat.description = "Weathered granite with reduced strength and increased porosity";
    
    mat.mechanical.youngs_modulus = GPaToPa(25.0);
    mat.mechanical.compressive_strength = MPaToPa(50.0);
    mat.mechanical.tensile_strength = MPaToPa(3.0);
    mat.mechanical.cohesion = MPaToPa(8.0);
    mat.mechanical.friction_angle = 40.0;
    mat.mechanical.biot_coefficient = 0.7;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.08;
    mat.hydraulic.setIsotropicPermeability(mDToM2(1.0));
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.is_porous = true;
    mat.recommended_model = "POROELASTIC";
    
    return mat;
}

MaterialDefinition Granodiorite() {
    MaterialDefinition mat;
    mat.name = "Granodiorite";
    mat.category = "Igneous";
    mat.subcategory = "Plutonic";
    mat.description = "Intermediate plutonic rock between granite and diorite";
    mat.reference = "Jaeger et al., 2007";
    
    mat.mechanical.density = 2700.0;
    mat.mechanical.youngs_modulus = GPaToPa(65.0);
    mat.mechanical.poisson_ratio = 0.24;
    mat.mechanical.tensile_strength = MPaToPa(12.0);
    mat.mechanical.compressive_strength = MPaToPa(190.0);
    mat.mechanical.cohesion = MPaToPa(35.0);
    mat.mechanical.friction_angle = 54.0;
    mat.mechanical.fracture_toughness_KIc = 1.6e6;
    mat.mechanical.biot_coefficient = 0.42;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.008;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.0005));
    
    mat.thermal.thermal_conductivity = 2.8;
    mat.thermal.specific_heat = 800.0;
    mat.thermal.thermal_expansion = 7.5e-6;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 180.0;
    mat.seismic.Qs = 90.0;
    
    mat.is_porous = false;
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition Diorite() {
    MaterialDefinition mat;
    mat.name = "Diorite";
    mat.category = "Igneous";
    mat.subcategory = "Plutonic";
    mat.description = "Intermediate intrusive igneous rock";
    mat.reference = "Carmichael, 1989";
    
    mat.mechanical.density = 2800.0;
    mat.mechanical.youngs_modulus = GPaToPa(70.0);
    mat.mechanical.poisson_ratio = 0.26;
    mat.mechanical.tensile_strength = MPaToPa(14.0);
    mat.mechanical.compressive_strength = MPaToPa(200.0);
    mat.mechanical.cohesion = MPaToPa(38.0);
    mat.mechanical.friction_angle = 52.0;
    mat.mechanical.fracture_toughness_KIc = 1.8e6;
    mat.mechanical.biot_coefficient = 0.40;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.005;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.0001));
    
    mat.thermal.thermal_conductivity = 2.5;
    mat.thermal.specific_heat = 820.0;
    mat.thermal.thermal_expansion = 7.0e-6;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 190.0;
    mat.seismic.Qs = 95.0;
    
    mat.is_porous = false;
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition Gabbro() {
    MaterialDefinition mat;
    mat.name = "Gabbro";
    mat.category = "Igneous";
    mat.subcategory = "Plutonic";
    mat.description = "Coarse-grained mafic intrusive rock, plutonic equivalent of basalt";
    mat.reference = "Jaeger et al., 2007";
    
    mat.mechanical.density = 2950.0;
    mat.mechanical.youngs_modulus = GPaToPa(85.0);
    mat.mechanical.poisson_ratio = 0.27;
    mat.mechanical.tensile_strength = MPaToPa(16.0);
    mat.mechanical.compressive_strength = MPaToPa(250.0);
    mat.mechanical.cohesion = MPaToPa(45.0);
    mat.mechanical.friction_angle = 50.0;
    mat.mechanical.fracture_toughness_KIc = 2.2e6;
    mat.mechanical.biot_coefficient = 0.35;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.003;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.00001));
    
    mat.thermal.thermal_conductivity = 2.2;
    mat.thermal.specific_heat = 850.0;
    mat.thermal.thermal_expansion = 6.5e-6;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 250.0;
    mat.seismic.Qs = 125.0;
    
    mat.is_porous = false;
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition Basalt() {
    MaterialDefinition mat;
    mat.name = "Basalt";
    mat.category = "Igneous";
    mat.subcategory = "Volcanic";
    mat.description = "Fine-grained mafic extrusive rock, most common volcanic rock";
    mat.reference = "Schultz, 1993; Jaeger et al., 2007";
    
    mat.mechanical.density = 2850.0;
    mat.mechanical.youngs_modulus = GPaToPa(75.0);
    mat.mechanical.poisson_ratio = 0.26;
    mat.mechanical.tensile_strength = MPaToPa(12.0);
    mat.mechanical.compressive_strength = MPaToPa(200.0);
    mat.mechanical.cohesion = MPaToPa(35.0);
    mat.mechanical.friction_angle = 48.0;
    mat.mechanical.fracture_toughness_KIc = 1.9e6;
    mat.mechanical.biot_coefficient = 0.38;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.02;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.01));
    
    mat.thermal.thermal_conductivity = 1.8;
    mat.thermal.specific_heat = 880.0;
    mat.thermal.thermal_expansion = 6.0e-6;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 150.0;
    mat.seismic.Qs = 75.0;
    
    mat.is_porous = false;
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition BasaltVesicular() {
    MaterialDefinition mat = Basalt();
    mat.name = "Basalt_Vesicular";
    mat.description = "Vesicular basalt with gas bubbles, higher porosity";
    
    mat.mechanical.density = 2400.0;
    mat.mechanical.youngs_modulus = GPaToPa(35.0);
    mat.mechanical.compressive_strength = MPaToPa(80.0);
    mat.mechanical.tensile_strength = MPaToPa(5.0);
    mat.mechanical.biot_coefficient = 0.65;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.20;
    mat.hydraulic.setIsotropicPermeability(mDToM2(100.0));
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.is_porous = true;
    mat.recommended_model = "POROELASTIC";
    
    return mat;
}

MaterialDefinition Andesite() {
    MaterialDefinition mat;
    mat.name = "Andesite";
    mat.category = "Igneous";
    mat.subcategory = "Volcanic";
    mat.description = "Intermediate volcanic rock, common in continental margins";
    mat.reference = "Carmichael, 1989";
    
    mat.mechanical.density = 2650.0;
    mat.mechanical.youngs_modulus = GPaToPa(55.0);
    mat.mechanical.poisson_ratio = 0.25;
    mat.mechanical.tensile_strength = MPaToPa(10.0);
    mat.mechanical.compressive_strength = MPaToPa(150.0);
    mat.mechanical.cohesion = MPaToPa(28.0);
    mat.mechanical.friction_angle = 50.0;
    mat.mechanical.fracture_toughness_KIc = 1.4e6;
    mat.mechanical.biot_coefficient = 0.50;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.05;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.1));
    
    mat.thermal.thermal_conductivity = 2.0;
    mat.thermal.specific_heat = 850.0;
    mat.thermal.thermal_expansion = 7.0e-6;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 120.0;
    mat.seismic.Qs = 60.0;
    
    mat.is_porous = false;
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition Rhyolite() {
    MaterialDefinition mat;
    mat.name = "Rhyolite";
    mat.category = "Igneous";
    mat.subcategory = "Volcanic";
    mat.description = "Felsic volcanic rock, extrusive equivalent of granite";
    mat.reference = "Jaeger et al., 2007";
    
    mat.mechanical.density = 2500.0;
    mat.mechanical.youngs_modulus = GPaToPa(50.0);
    mat.mechanical.poisson_ratio = 0.24;
    mat.mechanical.tensile_strength = MPaToPa(9.0);
    mat.mechanical.compressive_strength = MPaToPa(140.0);
    mat.mechanical.cohesion = MPaToPa(25.0);
    mat.mechanical.friction_angle = 52.0;
    mat.mechanical.fracture_toughness_KIc = 1.3e6;
    mat.mechanical.biot_coefficient = 0.55;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.06;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.05));
    
    mat.thermal.thermal_conductivity = 2.5;
    mat.thermal.specific_heat = 800.0;
    mat.thermal.thermal_expansion = 8.0e-6;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 100.0;
    mat.seismic.Qs = 50.0;
    
    mat.is_porous = false;
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition Obsidian() {
    MaterialDefinition mat;
    mat.name = "Obsidian";
    mat.category = "Igneous";
    mat.subcategory = "Volcanic";
    mat.description = "Volcanic glass formed by rapid cooling of felsic lava";
    mat.reference = "Carmichael, 1989";
    
    mat.mechanical.density = 2400.0;
    mat.mechanical.youngs_modulus = GPaToPa(70.0);
    mat.mechanical.poisson_ratio = 0.20;
    mat.mechanical.tensile_strength = MPaToPa(6.0);
    mat.mechanical.compressive_strength = MPaToPa(200.0);
    mat.mechanical.cohesion = MPaToPa(40.0);
    mat.mechanical.friction_angle = 45.0;
    mat.mechanical.fracture_toughness_KIc = 0.8e6;
    mat.mechanical.biot_coefficient = 0.10;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.001;
    mat.hydraulic.setIsotropicPermeability(mDToM2(1e-8));
    
    mat.thermal.thermal_conductivity = 1.5;
    mat.thermal.specific_heat = 750.0;
    mat.thermal.thermal_expansion = 5.0e-6;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 300.0;
    mat.seismic.Qs = 150.0;
    
    mat.is_porous = false;
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition Pumice() {
    MaterialDefinition mat;
    mat.name = "Pumice";
    mat.category = "Igneous";
    mat.subcategory = "Volcanic";
    mat.description = "Highly vesicular volcanic glass, very low density";
    mat.reference = "Cas & Wright, 1987";
    
    mat.mechanical.density = 800.0;
    mat.mechanical.youngs_modulus = GPaToPa(2.0);
    mat.mechanical.poisson_ratio = 0.20;
    mat.mechanical.tensile_strength = MPaToPa(0.5);
    mat.mechanical.compressive_strength = MPaToPa(5.0);
    mat.mechanical.cohesion = MPaToPa(1.0);
    mat.mechanical.friction_angle = 35.0;
    mat.mechanical.fracture_toughness_KIc = 0.1e6;
    mat.mechanical.biot_coefficient = 0.90;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.70;
    mat.hydraulic.setIsotropicPermeability(mDToM2(10000.0));
    
    mat.thermal.thermal_conductivity = 0.2;
    mat.thermal.specific_heat = 900.0;
    mat.thermal.thermal_expansion = 10.0e-6;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 20.0;
    mat.seismic.Qs = 10.0;
    
    mat.is_porous = true;
    mat.recommended_model = "POROELASTIC";
    
    return mat;
}

MaterialDefinition Tuff() {
    MaterialDefinition mat;
    mat.name = "Tuff";
    mat.category = "Igneous";
    mat.subcategory = "Pyroclastic";
    mat.description = "Consolidated volcanic ash and pyroclastic material";
    mat.reference = "Fisher & Schmincke, 1984";
    
    mat.mechanical.density = 1800.0;
    mat.mechanical.youngs_modulus = GPaToPa(8.0);
    mat.mechanical.poisson_ratio = 0.22;
    mat.mechanical.tensile_strength = MPaToPa(2.0);
    mat.mechanical.compressive_strength = MPaToPa(20.0);
    mat.mechanical.cohesion = MPaToPa(4.0);
    mat.mechanical.friction_angle = 38.0;
    mat.mechanical.fracture_toughness_KIc = 0.3e6;
    mat.mechanical.biot_coefficient = 0.75;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.35;
    mat.hydraulic.setIsotropicPermeability(mDToM2(10.0));
    
    mat.thermal.thermal_conductivity = 0.8;
    mat.thermal.specific_heat = 850.0;
    mat.thermal.thermal_expansion = 8.0e-6;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 40.0;
    mat.seismic.Qs = 20.0;
    
    mat.is_porous = true;
    mat.recommended_model = "POROELASTIC";
    
    return mat;
}

MaterialDefinition TuffWelded() {
    MaterialDefinition mat = Tuff();
    mat.name = "Tuff_Welded";
    mat.description = "Welded tuff with fused glass shards, higher strength";
    
    mat.mechanical.density = 2200.0;
    mat.mechanical.youngs_modulus = GPaToPa(25.0);
    mat.mechanical.compressive_strength = MPaToPa(60.0);
    mat.mechanical.tensile_strength = MPaToPa(5.0);
    mat.mechanical.cohesion = MPaToPa(12.0);
    mat.mechanical.biot_coefficient = 0.55;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.15;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.5));
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    
    return mat;
}

MaterialDefinition Dacite() {
    MaterialDefinition mat;
    mat.name = "Dacite";
    mat.category = "Igneous";
    mat.subcategory = "Volcanic";
    mat.description = "Intermediate to felsic volcanic rock";
    mat.reference = "Carmichael, 1989";
    
    mat.mechanical.density = 2600.0;
    mat.mechanical.youngs_modulus = GPaToPa(52.0);
    mat.mechanical.poisson_ratio = 0.24;
    mat.mechanical.tensile_strength = MPaToPa(9.5);
    mat.mechanical.compressive_strength = MPaToPa(145.0);
    mat.mechanical.cohesion = MPaToPa(26.0);
    mat.mechanical.friction_angle = 51.0;
    mat.mechanical.fracture_toughness_KIc = 1.35e6;
    mat.mechanical.biot_coefficient = 0.52;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.05;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.08));
    
    mat.thermal.thermal_conductivity = 2.2;
    mat.thermal.specific_heat = 830.0;
    mat.thermal.thermal_expansion = 7.5e-6;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 110.0;
    mat.seismic.Qs = 55.0;
    
    mat.is_porous = false;
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition Peridotite() {
    MaterialDefinition mat;
    mat.name = "Peridotite";
    mat.category = "Igneous";
    mat.subcategory = "Ultramafic";
    mat.description = "Ultramafic rock of the upper mantle, olivine-rich";
    mat.reference = "Christensen, 1974";
    
    mat.mechanical.density = 3300.0;
    mat.mechanical.youngs_modulus = GPaToPa(150.0);
    mat.mechanical.poisson_ratio = 0.28;
    mat.mechanical.tensile_strength = MPaToPa(20.0);
    mat.mechanical.compressive_strength = MPaToPa(350.0);
    mat.mechanical.cohesion = MPaToPa(55.0);
    mat.mechanical.friction_angle = 45.0;
    mat.mechanical.fracture_toughness_KIc = 3.0e6;
    mat.mechanical.biot_coefficient = 0.20;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.001;
    mat.hydraulic.setIsotropicPermeability(mDToM2(1e-7));
    
    mat.thermal.thermal_conductivity = 3.5;
    mat.thermal.specific_heat = 1000.0;
    mat.thermal.thermal_expansion = 3.0e-5;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 500.0;
    mat.seismic.Qs = 250.0;
    
    mat.is_porous = false;
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition Dunite() {
    MaterialDefinition mat;
    mat.name = "Dunite";
    mat.category = "Igneous";
    mat.subcategory = "Ultramafic";
    mat.description = "Ultramafic rock >90% olivine";
    mat.reference = "Christensen, 1974";
    
    mat.mechanical.density = 3300.0;
    mat.mechanical.youngs_modulus = GPaToPa(160.0);
    mat.mechanical.poisson_ratio = 0.27;
    mat.mechanical.tensile_strength = MPaToPa(22.0);
    mat.mechanical.compressive_strength = MPaToPa(380.0);
    mat.mechanical.cohesion = MPaToPa(60.0);
    mat.mechanical.friction_angle = 44.0;
    mat.mechanical.fracture_toughness_KIc = 3.2e6;
    mat.mechanical.biot_coefficient = 0.18;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.0008;
    mat.hydraulic.setIsotropicPermeability(mDToM2(1e-8));
    
    mat.thermal.thermal_conductivity = 4.0;
    mat.thermal.specific_heat = 1050.0;
    mat.thermal.thermal_expansion = 2.8e-5;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 550.0;
    mat.seismic.Qs = 275.0;
    
    mat.is_porous = false;
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition Syenite() {
    MaterialDefinition mat;
    mat.name = "Syenite";
    mat.category = "Igneous";
    mat.subcategory = "Plutonic";
    mat.description = "Coarse-grained intrusive rock with alkali feldspar, quartz-poor";
    mat.reference = "Carmichael, 1989";
    
    mat.mechanical.density = 2750.0;
    mat.mechanical.youngs_modulus = GPaToPa(62.0);
    mat.mechanical.poisson_ratio = 0.25;
    mat.mechanical.tensile_strength = MPaToPa(11.0);
    mat.mechanical.compressive_strength = MPaToPa(175.0);
    mat.mechanical.cohesion = MPaToPa(32.0);
    mat.mechanical.friction_angle = 53.0;
    mat.mechanical.fracture_toughness_KIc = 1.55e6;
    mat.mechanical.biot_coefficient = 0.43;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.009;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.0008));
    
    mat.thermal.thermal_conductivity = 2.6;
    mat.thermal.specific_heat = 810.0;
    mat.thermal.thermal_expansion = 7.8e-6;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 170.0;
    mat.seismic.Qs = 85.0;
    
    mat.is_porous = false;
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition Diabase() {
    MaterialDefinition mat;
    mat.name = "Diabase";
    mat.category = "Igneous";
    mat.subcategory = "Hypabyssal";
    mat.description = "Medium-grained mafic intrusive rock (dolerite)";
    mat.reference = "Jaeger et al., 2007";
    
    mat.mechanical.density = 2900.0;
    mat.mechanical.youngs_modulus = GPaToPa(80.0);
    mat.mechanical.poisson_ratio = 0.27;
    mat.mechanical.tensile_strength = MPaToPa(15.0);
    mat.mechanical.compressive_strength = MPaToPa(230.0);
    mat.mechanical.cohesion = MPaToPa(42.0);
    mat.mechanical.friction_angle = 49.0;
    mat.mechanical.fracture_toughness_KIc = 2.1e6;
    mat.mechanical.biot_coefficient = 0.36;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.004;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.00005));
    
    mat.thermal.thermal_conductivity = 2.3;
    mat.thermal.specific_heat = 860.0;
    mat.thermal.thermal_expansion = 6.2e-6;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 220.0;
    mat.seismic.Qs = 110.0;
    
    mat.is_porous = false;
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition Pegmatite() {
    MaterialDefinition mat;
    mat.name = "Pegmatite";
    mat.category = "Igneous";
    mat.subcategory = "Plutonic";
    mat.description = "Very coarse-grained granitic rock";
    mat.reference = "London, 2008";
    
    mat.mechanical.density = 2600.0;
    mat.mechanical.youngs_modulus = GPaToPa(55.0);
    mat.mechanical.poisson_ratio = 0.24;
    mat.mechanical.tensile_strength = MPaToPa(8.0);
    mat.mechanical.compressive_strength = MPaToPa(150.0);
    mat.mechanical.cohesion = MPaToPa(25.0);
    mat.mechanical.friction_angle = 54.0;
    mat.mechanical.fracture_toughness_KIc = 1.2e6;
    mat.mechanical.biot_coefficient = 0.48;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.015;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.01));
    
    mat.thermal.thermal_conductivity = 3.2;
    mat.thermal.specific_heat = 780.0;
    mat.thermal.thermal_expansion = 8.2e-6;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 150.0;
    mat.seismic.Qs = 75.0;
    
    mat.is_porous = false;
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition Porphyry() {
    MaterialDefinition mat;
    mat.name = "Porphyry";
    mat.category = "Igneous";
    mat.subcategory = "Hypabyssal";
    mat.description = "Rock with large crystals in fine-grained groundmass";
    mat.reference = "Sillitoe, 2010";
    
    mat.mechanical.density = 2700.0;
    mat.mechanical.youngs_modulus = GPaToPa(58.0);
    mat.mechanical.poisson_ratio = 0.25;
    mat.mechanical.tensile_strength = MPaToPa(10.5);
    mat.mechanical.compressive_strength = MPaToPa(160.0);
    mat.mechanical.cohesion = MPaToPa(29.0);
    mat.mechanical.friction_angle = 52.0;
    mat.mechanical.fracture_toughness_KIc = 1.45e6;
    mat.mechanical.biot_coefficient = 0.47;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.03;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.05));
    
    mat.thermal.thermal_conductivity = 2.4;
    mat.thermal.specific_heat = 820.0;
    mat.thermal.thermal_expansion = 7.6e-6;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 130.0;
    mat.seismic.Qs = 65.0;
    
    mat.is_porous = false;
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition Aplite() {
    MaterialDefinition mat;
    mat.name = "Aplite";
    mat.category = "Igneous";
    mat.subcategory = "Hypabyssal";
    mat.description = "Fine-grained granitic rock, typically in dikes";
    mat.reference = "Carmichael, 1989";
    
    mat.mechanical.density = 2620.0;
    mat.mechanical.youngs_modulus = GPaToPa(58.0);
    mat.mechanical.poisson_ratio = 0.23;
    mat.mechanical.tensile_strength = MPaToPa(9.5);
    mat.mechanical.compressive_strength = MPaToPa(165.0);
    mat.mechanical.cohesion = MPaToPa(28.0);
    mat.mechanical.friction_angle = 56.0;
    mat.mechanical.fracture_toughness_KIc = 1.4e6;
    mat.mechanical.biot_coefficient = 0.44;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.008;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.0003));
    
    mat.thermal.thermal_conductivity = 3.1;
    mat.thermal.specific_heat = 785.0;
    mat.thermal.thermal_expansion = 8.0e-6;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 180.0;
    mat.seismic.Qs = 90.0;
    
    mat.is_porous = false;
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition Komatiite() {
    MaterialDefinition mat;
    mat.name = "Komatiite";
    mat.category = "Igneous";
    mat.subcategory = "Ultramafic";
    mat.description = "Archean ultramafic volcanic rock, very high MgO";
    mat.reference = "Arndt & Nisbet, 1982";
    
    mat.mechanical.density = 3100.0;
    mat.mechanical.youngs_modulus = GPaToPa(120.0);
    mat.mechanical.poisson_ratio = 0.27;
    mat.mechanical.tensile_strength = MPaToPa(18.0);
    mat.mechanical.compressive_strength = MPaToPa(280.0);
    mat.mechanical.cohesion = MPaToPa(48.0);
    mat.mechanical.friction_angle = 46.0;
    mat.mechanical.fracture_toughness_KIc = 2.5e6;
    mat.mechanical.biot_coefficient = 0.28;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.008;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.001));
    
    mat.thermal.thermal_conductivity = 3.2;
    mat.thermal.specific_heat = 920.0;
    mat.thermal.thermal_expansion = 4.0e-5;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 300.0;
    mat.seismic.Qs = 150.0;
    
    mat.is_porous = false;
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition Phonolite() {
    MaterialDefinition mat;
    mat.name = "Phonolite";
    mat.category = "Igneous";
    mat.subcategory = "Volcanic";
    mat.description = "Alkaline volcanic rock, extrusive equivalent of nepheline syenite";
    mat.reference = "Carmichael, 1989";
    
    mat.mechanical.density = 2550.0;
    mat.mechanical.youngs_modulus = GPaToPa(48.0);
    mat.mechanical.poisson_ratio = 0.24;
    mat.mechanical.tensile_strength = MPaToPa(8.5);
    mat.mechanical.compressive_strength = MPaToPa(130.0);
    mat.mechanical.cohesion = MPaToPa(23.0);
    mat.mechanical.friction_angle = 50.0;
    mat.mechanical.fracture_toughness_KIc = 1.2e6;
    mat.mechanical.biot_coefficient = 0.52;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.04;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.05));
    
    mat.thermal.thermal_conductivity = 2.1;
    mat.thermal.specific_heat = 810.0;
    mat.thermal.thermal_expansion = 7.2e-6;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 100.0;
    mat.seismic.Qs = 50.0;
    
    mat.is_porous = false;
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition Latite() {
    MaterialDefinition mat;
    mat.name = "Latite";
    mat.category = "Igneous";
    mat.subcategory = "Volcanic";
    mat.description = "Intermediate volcanic rock between andesite and trachyte";
    mat.reference = "Carmichael, 1989";
    
    mat.mechanical.density = 2630.0;
    mat.mechanical.youngs_modulus = GPaToPa(53.0);
    mat.mechanical.poisson_ratio = 0.25;
    mat.mechanical.tensile_strength = MPaToPa(9.8);
    mat.mechanical.compressive_strength = MPaToPa(148.0);
    mat.mechanical.cohesion = MPaToPa(27.0);
    mat.mechanical.friction_angle = 51.0;
    mat.mechanical.fracture_toughness_KIc = 1.38e6;
    mat.mechanical.biot_coefficient = 0.51;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.045;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.07));
    
    mat.thermal.thermal_conductivity = 2.15;
    mat.thermal.specific_heat = 835.0;
    mat.thermal.thermal_expansion = 7.3e-6;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 115.0;
    mat.seismic.Qs = 57.0;
    
    mat.is_porous = false;
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition Monzonite() {
    MaterialDefinition mat;
    mat.name = "Monzonite";
    mat.category = "Igneous";
    mat.subcategory = "Plutonic";
    mat.description = "Intermediate plutonic rock between syenite and diorite";
    mat.reference = "Carmichael, 1989";
    
    mat.mechanical.density = 2770.0;
    mat.mechanical.youngs_modulus = GPaToPa(66.0);
    mat.mechanical.poisson_ratio = 0.25;
    mat.mechanical.tensile_strength = MPaToPa(13.0);
    mat.mechanical.compressive_strength = MPaToPa(185.0);
    mat.mechanical.cohesion = MPaToPa(36.0);
    mat.mechanical.friction_angle = 53.0;
    mat.mechanical.fracture_toughness_KIc = 1.7e6;
    mat.mechanical.biot_coefficient = 0.41;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.007;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.0004));
    
    mat.thermal.thermal_conductivity = 2.55;
    mat.thermal.specific_heat = 815.0;
    mat.thermal.thermal_expansion = 7.4e-6;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 175.0;
    mat.seismic.Qs = 87.0;
    
    mat.is_porous = false;
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition Troctolite() {
    MaterialDefinition mat;
    mat.name = "Troctolite";
    mat.category = "Igneous";
    mat.subcategory = "Plutonic";
    mat.description = "Mafic plutonic rock of olivine and plagioclase";
    mat.reference = "Wager & Brown, 1968";
    
    mat.mechanical.density = 2950.0;
    mat.mechanical.youngs_modulus = GPaToPa(90.0);
    mat.mechanical.poisson_ratio = 0.27;
    mat.mechanical.tensile_strength = MPaToPa(17.0);
    mat.mechanical.compressive_strength = MPaToPa(260.0);
    mat.mechanical.cohesion = MPaToPa(47.0);
    mat.mechanical.friction_angle = 48.0;
    mat.mechanical.fracture_toughness_KIc = 2.3e6;
    mat.mechanical.biot_coefficient = 0.33;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.003;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.00002));
    
    mat.thermal.thermal_conductivity = 2.4;
    mat.thermal.specific_heat = 870.0;
    mat.thermal.thermal_expansion = 6.3e-6;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 260.0;
    mat.seismic.Qs = 130.0;
    
    mat.is_porous = false;
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition Norite() {
    MaterialDefinition mat;
    mat.name = "Norite";
    mat.category = "Igneous";
    mat.subcategory = "Plutonic";
    mat.description = "Mafic plutonic rock with orthopyroxene, related to gabbro";
    mat.reference = "Wager & Brown, 1968";
    
    mat.mechanical.density = 2980.0;
    mat.mechanical.youngs_modulus = GPaToPa(88.0);
    mat.mechanical.poisson_ratio = 0.27;
    mat.mechanical.tensile_strength = MPaToPa(16.5);
    mat.mechanical.compressive_strength = MPaToPa(255.0);
    mat.mechanical.cohesion = MPaToPa(46.0);
    mat.mechanical.friction_angle = 49.0;
    mat.mechanical.fracture_toughness_KIc = 2.25e6;
    mat.mechanical.biot_coefficient = 0.34;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.0025;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.000015));
    
    mat.thermal.thermal_conductivity = 2.25;
    mat.thermal.specific_heat = 855.0;
    mat.thermal.thermal_expansion = 6.4e-6;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 255.0;
    mat.seismic.Qs = 127.0;
    
    mat.is_porous = false;
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition Anorthosite() {
    MaterialDefinition mat;
    mat.name = "Anorthosite";
    mat.category = "Igneous";
    mat.subcategory = "Plutonic";
    mat.description = "Plutonic rock >90% plagioclase feldspar";
    mat.reference = "Ashwal, 1993";
    
    mat.mechanical.density = 2750.0;
    mat.mechanical.youngs_modulus = GPaToPa(75.0);
    mat.mechanical.poisson_ratio = 0.26;
    mat.mechanical.tensile_strength = MPaToPa(13.0);
    mat.mechanical.compressive_strength = MPaToPa(200.0);
    mat.mechanical.cohesion = MPaToPa(38.0);
    mat.mechanical.friction_angle = 50.0;
    mat.mechanical.fracture_toughness_KIc = 1.9e6;
    mat.mechanical.biot_coefficient = 0.38;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.005;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.0001));
    
    mat.thermal.thermal_conductivity = 2.0;
    mat.thermal.specific_heat = 840.0;
    mat.thermal.thermal_expansion = 6.8e-6;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 200.0;
    mat.seismic.Qs = 100.0;
    
    mat.is_porous = false;
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition Kimberlite() {
    MaterialDefinition mat;
    mat.name = "Kimberlite";
    mat.category = "Igneous";
    mat.subcategory = "Ultramafic";
    mat.description = "Diamond-bearing ultramafic volcanic rock";
    mat.reference = "Mitchell, 1986";
    
    mat.mechanical.density = 2850.0;
    mat.mechanical.youngs_modulus = GPaToPa(45.0);
    mat.mechanical.poisson_ratio = 0.28;
    mat.mechanical.tensile_strength = MPaToPa(6.0);
    mat.mechanical.compressive_strength = MPaToPa(80.0);
    mat.mechanical.cohesion = MPaToPa(15.0);
    mat.mechanical.friction_angle = 42.0;
    mat.mechanical.fracture_toughness_KIc = 1.0e6;
    mat.mechanical.biot_coefficient = 0.60;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.15;
    mat.hydraulic.setIsotropicPermeability(mDToM2(1.0));
    
    mat.thermal.thermal_conductivity = 2.0;
    mat.thermal.specific_heat = 900.0;
    mat.thermal.thermal_expansion = 8.0e-6;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 60.0;
    mat.seismic.Qs = 30.0;
    
    mat.is_porous = true;
    mat.recommended_model = "POROELASTIC";
    
    return mat;
}

MaterialDefinition Lamprophyre() {
    MaterialDefinition mat;
    mat.name = "Lamprophyre";
    mat.category = "Igneous";
    mat.subcategory = "Hypabyssal";
    mat.description = "Porphyritic mafic to ultramafic intrusive rock";
    mat.reference = "Rock, 1991";
    
    mat.mechanical.density = 2850.0;
    mat.mechanical.youngs_modulus = GPaToPa(60.0);
    mat.mechanical.poisson_ratio = 0.27;
    mat.mechanical.tensile_strength = MPaToPa(10.0);
    mat.mechanical.compressive_strength = MPaToPa(150.0);
    mat.mechanical.cohesion = MPaToPa(28.0);
    mat.mechanical.friction_angle = 46.0;
    mat.mechanical.fracture_toughness_KIc = 1.5e6;
    mat.mechanical.biot_coefficient = 0.50;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.05;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.1));
    
    mat.thermal.thermal_conductivity = 2.3;
    mat.thermal.specific_heat = 880.0;
    mat.thermal.thermal_expansion = 7.0e-6;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 100.0;
    mat.seismic.Qs = 50.0;
    
    mat.is_porous = false;
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

} // namespace Igneous

// =============================================================================
// Sedimentary Rock Definitions - namespace Materials::Sedimentary
// =============================================================================

namespace Sedimentary {

MaterialDefinition Sandstone() {
    MaterialDefinition mat;
    mat.name = "Sandstone";
    mat.category = "Sedimentary";
    mat.subcategory = "Clastic";
    mat.description = "Medium-grained clastic rock, 0.063-2mm grain size";
    mat.reference = "Jaeger et al., 2007; Mavko et al., 2009";
    
    mat.mechanical.density = 2300.0;
    mat.mechanical.youngs_modulus = GPaToPa(20.0);
    mat.mechanical.poisson_ratio = 0.20;
    mat.mechanical.tensile_strength = MPaToPa(4.0);
    mat.mechanical.compressive_strength = MPaToPa(50.0);
    mat.mechanical.cohesion = MPaToPa(8.0);
    mat.mechanical.friction_angle = 35.0;
    mat.mechanical.dilation_angle = 10.0;
    mat.mechanical.fracture_toughness_KIc = 0.5e6;
    mat.mechanical.biot_coefficient = 0.75;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.20;
    mat.hydraulic.setIsotropicPermeability(mDToM2(100.0));
    mat.hydraulic.compressibility = 5.0e-10;
    mat.hydraulic.Swc = 0.20;
    mat.hydraulic.Sor = 0.20;
    mat.hydraulic.nw = 2.5;
    mat.hydraulic.no = 2.0;
    
    mat.thermal.thermal_conductivity = 2.5;
    mat.thermal.specific_heat = 920.0;
    mat.thermal.thermal_expansion = 1.2e-5;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 50.0;
    mat.seismic.Qs = 25.0;
    
    mat.is_porous = true;
    mat.recommended_model = "POROELASTIC";
    
    return mat;
}

MaterialDefinition SandstoneQuartz() {
    MaterialDefinition mat = Sandstone();
    mat.name = "Sandstone_Quartz";
    mat.description = "Quartz-rich sandstone (>95% quartz), well-cemented";
    
    mat.mechanical.density = 2400.0;
    mat.mechanical.youngs_modulus = GPaToPa(30.0);
    mat.mechanical.compressive_strength = MPaToPa(80.0);
    mat.mechanical.tensile_strength = MPaToPa(6.0);
    mat.mechanical.cohesion = MPaToPa(15.0);
    mat.mechanical.friction_angle = 40.0;
    mat.mechanical.biot_coefficient = 0.65;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.12;
    mat.hydraulic.setIsotropicPermeability(mDToM2(50.0));
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    
    return mat;
}

MaterialDefinition SandstoneFeldspathic() {
    MaterialDefinition mat = Sandstone();
    mat.name = "Sandstone_Feldspathic";
    mat.description = "Feldspar-rich sandstone (arkose), partially weathered";
    
    mat.mechanical.density = 2350.0;
    mat.mechanical.youngs_modulus = GPaToPa(18.0);
    mat.mechanical.compressive_strength = MPaToPa(45.0);
    mat.mechanical.biot_coefficient = 0.78;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.22;
    mat.hydraulic.setIsotropicPermeability(mDToM2(150.0));
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    
    return mat;
}

MaterialDefinition SandstoneLithic() {
    MaterialDefinition mat = Sandstone();
    mat.name = "Sandstone_Lithic";
    mat.description = "Lithic sandstone with rock fragments";
    
    mat.mechanical.density = 2400.0;
    mat.mechanical.youngs_modulus = GPaToPa(22.0);
    mat.mechanical.compressive_strength = MPaToPa(55.0);
    mat.mechanical.biot_coefficient = 0.72;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.18;
    mat.hydraulic.setIsotropicPermeability(mDToM2(80.0));
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    
    return mat;
}

MaterialDefinition SandstoneArkose() {
    MaterialDefinition mat = Sandstone();
    mat.name = "Sandstone_Arkose";
    mat.description = "Arkose sandstone with >25% feldspar";
    
    mat.mechanical.density = 2350.0;
    mat.mechanical.youngs_modulus = GPaToPa(17.0);
    mat.mechanical.compressive_strength = MPaToPa(42.0);
    mat.mechanical.cohesion = MPaToPa(7.0);
    mat.mechanical.biot_coefficient = 0.80;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.24;
    mat.hydraulic.setIsotropicPermeability(mDToM2(200.0));
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    
    return mat;
}

MaterialDefinition SandstoneGreywacke() {
    MaterialDefinition mat = Sandstone();
    mat.name = "Sandstone_Greywacke";
    mat.description = "Greywacke - poorly sorted, clay-rich sandstone";
    
    mat.mechanical.density = 2500.0;
    mat.mechanical.youngs_modulus = GPaToPa(35.0);
    mat.mechanical.compressive_strength = MPaToPa(90.0);
    mat.mechanical.tensile_strength = MPaToPa(7.0);
    mat.mechanical.cohesion = MPaToPa(18.0);
    mat.mechanical.friction_angle = 38.0;
    mat.mechanical.biot_coefficient = 0.60;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.08;
    mat.hydraulic.setIsotropicPermeability(mDToM2(1.0));
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    
    return mat;
}

MaterialDefinition SandstoneTight() {
    MaterialDefinition mat = Sandstone();
    mat.name = "Sandstone_Tight";
    mat.description = "Tight sandstone with low permeability, well-cemented";
    
    mat.mechanical.density = 2500.0;
    mat.mechanical.youngs_modulus = GPaToPa(35.0);
    mat.mechanical.compressive_strength = MPaToPa(100.0);
    mat.mechanical.tensile_strength = MPaToPa(8.0);
    mat.mechanical.biot_coefficient = 0.55;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.06;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.1));
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    
    return mat;
}

MaterialDefinition Siltstone() {
    MaterialDefinition mat;
    mat.name = "Siltstone";
    mat.category = "Sedimentary";
    mat.subcategory = "Clastic";
    mat.description = "Fine-grained clastic rock, 0.004-0.063mm grain size";
    mat.reference = "Jaeger et al., 2007";
    
    mat.mechanical.density = 2450.0;
    mat.mechanical.youngs_modulus = GPaToPa(25.0);
    mat.mechanical.poisson_ratio = 0.22;
    mat.mechanical.tensile_strength = MPaToPa(5.0);
    mat.mechanical.compressive_strength = MPaToPa(60.0);
    mat.mechanical.cohesion = MPaToPa(10.0);
    mat.mechanical.friction_angle = 32.0;
    mat.mechanical.fracture_toughness_KIc = 0.6e6;
    mat.mechanical.biot_coefficient = 0.70;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.15;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.5));
    
    mat.thermal.thermal_conductivity = 2.2;
    mat.thermal.specific_heat = 900.0;
    mat.thermal.thermal_expansion = 1.1e-5;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 40.0;
    mat.seismic.Qs = 20.0;
    
    mat.is_porous = true;
    mat.recommended_model = "POROELASTIC";
    
    return mat;
}

MaterialDefinition Shale() {
    MaterialDefinition mat;
    mat.name = "Shale";
    mat.category = "Sedimentary";
    mat.subcategory = "Clastic";
    mat.description = "Fine-grained fissile rock, clay minerals dominant";
    mat.reference = "Sone & Zoback, 2013; Josh et al., 2012";
    
    mat.mechanical.density = 2500.0;
    mat.mechanical.youngs_modulus = GPaToPa(15.0);
    mat.mechanical.poisson_ratio = 0.25;
    mat.mechanical.tensile_strength = MPaToPa(3.0);
    mat.mechanical.compressive_strength = MPaToPa(40.0);
    mat.mechanical.cohesion = MPaToPa(5.0);
    mat.mechanical.friction_angle = 25.0;
    mat.mechanical.fracture_toughness_KIc = 0.4e6;
    mat.mechanical.biot_coefficient = 0.80;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.10;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.001));
    mat.hydraulic.compressibility = 1.0e-9;
    
    mat.thermal.thermal_conductivity = 1.5;
    mat.thermal.specific_heat = 850.0;
    mat.thermal.thermal_expansion = 1.0e-5;
    
    mat.anisotropic.is_anisotropic = true;
    mat.anisotropic.symmetry_type = "VTI";
    mat.anisotropic.E_vertical = GPaToPa(12.0);
    mat.anisotropic.E_horizontal = GPaToPa(18.0);
    mat.anisotropic.nu_vh = 0.20;
    mat.anisotropic.nu_hh = 0.25;
    mat.anisotropic.epsilon = 0.15;
    mat.anisotropic.delta = 0.10;
    mat.anisotropic.gamma = 0.20;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 30.0;
    mat.seismic.Qs = 15.0;
    
    mat.is_porous = true;
    mat.is_anisotropic_flag = true;
    mat.recommended_model = "ANISOTROPIC_VTI";
    
    return mat;
}

MaterialDefinition ShaleOrganic() {
    MaterialDefinition mat = Shale();
    mat.name = "Shale_Organic";
    mat.description = "Organic-rich shale (source rock), TOC > 2%";
    mat.reference = "Sone & Zoback, 2013";
    
    mat.mechanical.density = 2400.0;
    mat.mechanical.youngs_modulus = GPaToPa(12.0);
    mat.mechanical.compressive_strength = MPaToPa(35.0);
    mat.mechanical.biot_coefficient = 0.85;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.08;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.0001));
    
    mat.thermal.thermal_conductivity = 1.2;
    mat.thermal.heat_production = 3.0e-6;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    
    return mat;
}

MaterialDefinition ShaleSiliceous() {
    MaterialDefinition mat = Shale();
    mat.name = "Shale_Siliceous";
    mat.description = "Silica-rich shale, more brittle";
    
    mat.mechanical.density = 2550.0;
    mat.mechanical.youngs_modulus = GPaToPa(25.0);
    mat.mechanical.compressive_strength = MPaToPa(70.0);
    mat.mechanical.tensile_strength = MPaToPa(6.0);
    mat.mechanical.cohesion = MPaToPa(12.0);
    mat.mechanical.friction_angle = 30.0;
    mat.mechanical.biot_coefficient = 0.65;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.06;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.0005));
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    
    return mat;
}

MaterialDefinition ShaleCalcareous() {
    MaterialDefinition mat = Shale();
    mat.name = "Shale_Calcareous";
    mat.description = "Calcareous shale with carbonate content";
    
    mat.mechanical.density = 2550.0;
    mat.mechanical.youngs_modulus = GPaToPa(22.0);
    mat.mechanical.compressive_strength = MPaToPa(60.0);
    mat.mechanical.cohesion = MPaToPa(10.0);
    mat.mechanical.biot_coefficient = 0.70;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.08;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.0008));
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    
    return mat;
}

MaterialDefinition Mudstone() {
    MaterialDefinition mat;
    mat.name = "Mudstone";
    mat.category = "Sedimentary";
    mat.subcategory = "Clastic";
    mat.description = "Fine-grained, non-fissile claystone";
    mat.reference = "Jaeger et al., 2007";
    
    mat.mechanical.density = 2450.0;
    mat.mechanical.youngs_modulus = GPaToPa(12.0);
    mat.mechanical.poisson_ratio = 0.28;
    mat.mechanical.tensile_strength = MPaToPa(2.5);
    mat.mechanical.compressive_strength = MPaToPa(30.0);
    mat.mechanical.cohesion = MPaToPa(4.0);
    mat.mechanical.friction_angle = 22.0;
    mat.mechanical.fracture_toughness_KIc = 0.3e6;
    mat.mechanical.biot_coefficient = 0.85;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.20;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.01));
    
    mat.thermal.thermal_conductivity = 1.4;
    mat.thermal.specific_heat = 860.0;
    mat.thermal.thermal_expansion = 1.0e-5;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 25.0;
    mat.seismic.Qs = 12.0;
    
    mat.is_porous = true;
    mat.recommended_model = "POROELASTIC";
    
    return mat;
}

MaterialDefinition Claystone() {
    MaterialDefinition mat;
    mat.name = "Claystone";
    mat.category = "Sedimentary";
    mat.subcategory = "Clastic";
    mat.description = "Rock composed primarily of clay-sized particles";
    mat.reference = "Horsrud, 2001";
    
    mat.mechanical.density = 2400.0;
    mat.mechanical.youngs_modulus = GPaToPa(8.0);
    mat.mechanical.poisson_ratio = 0.30;
    mat.mechanical.tensile_strength = MPaToPa(1.5);
    mat.mechanical.compressive_strength = MPaToPa(20.0);
    mat.mechanical.cohesion = MPaToPa(3.0);
    mat.mechanical.friction_angle = 18.0;
    mat.mechanical.fracture_toughness_KIc = 0.2e6;
    mat.mechanical.biot_coefficient = 0.90;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.30;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.001));
    mat.hydraulic.compressibility = 2.0e-9;
    
    mat.thermal.thermal_conductivity = 1.2;
    mat.thermal.specific_heat = 870.0;
    mat.thermal.thermal_expansion = 1.5e-5;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 20.0;
    mat.seismic.Qs = 10.0;
    
    mat.is_porous = true;
    mat.recommended_model = "POROELASTIC";
    
    return mat;
}

MaterialDefinition Conglomerate() {
    MaterialDefinition mat;
    mat.name = "Conglomerate";
    mat.category = "Sedimentary";
    mat.subcategory = "Clastic";
    mat.description = "Coarse-grained clastic rock with rounded clasts >2mm";
    mat.reference = "Jaeger et al., 2007";
    
    mat.mechanical.density = 2400.0;
    mat.mechanical.youngs_modulus = GPaToPa(25.0);
    mat.mechanical.poisson_ratio = 0.22;
    mat.mechanical.tensile_strength = MPaToPa(5.0);
    mat.mechanical.compressive_strength = MPaToPa(65.0);
    mat.mechanical.cohesion = MPaToPa(10.0);
    mat.mechanical.friction_angle = 38.0;
    mat.mechanical.fracture_toughness_KIc = 0.7e6;
    mat.mechanical.biot_coefficient = 0.70;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.18;
    mat.hydraulic.setIsotropicPermeability(mDToM2(500.0));
    
    mat.thermal.thermal_conductivity = 2.8;
    mat.thermal.specific_heat = 900.0;
    mat.thermal.thermal_expansion = 1.2e-5;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 45.0;
    mat.seismic.Qs = 22.0;
    
    mat.is_porous = true;
    mat.recommended_model = "POROELASTIC";
    
    return mat;
}

MaterialDefinition Breccia() {
    MaterialDefinition mat;
    mat.name = "Breccia";
    mat.category = "Sedimentary";
    mat.subcategory = "Clastic";
    mat.description = "Coarse-grained rock with angular clasts";
    mat.reference = "Jaeger et al., 2007";
    
    mat.mechanical.density = 2450.0;
    mat.mechanical.youngs_modulus = GPaToPa(22.0);
    mat.mechanical.poisson_ratio = 0.23;
    mat.mechanical.tensile_strength = MPaToPa(4.0);
    mat.mechanical.compressive_strength = MPaToPa(55.0);
    mat.mechanical.cohesion = MPaToPa(8.0);
    mat.mechanical.friction_angle = 35.0;
    mat.mechanical.fracture_toughness_KIc = 0.6e6;
    mat.mechanical.biot_coefficient = 0.72;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.20;
    mat.hydraulic.setIsotropicPermeability(mDToM2(800.0));
    
    mat.thermal.thermal_conductivity = 2.5;
    mat.thermal.specific_heat = 890.0;
    mat.thermal.thermal_expansion = 1.3e-5;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 40.0;
    mat.seismic.Qs = 20.0;
    
    mat.is_porous = true;
    mat.recommended_model = "POROELASTIC";
    
    return mat;
}

MaterialDefinition Tillite() {
    MaterialDefinition mat;
    mat.name = "Tillite";
    mat.category = "Sedimentary";
    mat.subcategory = "Clastic";
    mat.description = "Lithified glacial till, poorly sorted";
    mat.reference = "Jaeger et al., 2007";
    
    mat.mechanical.density = 2500.0;
    mat.mechanical.youngs_modulus = GPaToPa(28.0);
    mat.mechanical.poisson_ratio = 0.24;
    mat.mechanical.tensile_strength = MPaToPa(5.5);
    mat.mechanical.compressive_strength = MPaToPa(70.0);
    mat.mechanical.cohesion = MPaToPa(12.0);
    mat.mechanical.friction_angle = 36.0;
    mat.mechanical.fracture_toughness_KIc = 0.75e6;
    mat.mechanical.biot_coefficient = 0.65;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.12;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.5));
    
    mat.thermal.thermal_conductivity = 2.4;
    mat.thermal.specific_heat = 880.0;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 55.0;
    mat.seismic.Qs = 27.0;
    
    mat.is_porous = true;
    mat.recommended_model = "POROELASTIC";
    
    return mat;
}

MaterialDefinition Loess() {
    MaterialDefinition mat;
    mat.name = "Loess";
    mat.category = "Sedimentary";
    mat.subcategory = "Clastic";
    mat.description = "Wind-deposited silt, weakly cemented";
    mat.reference = "Derbyshire et al., 1995";
    
    mat.mechanical.density = 1600.0;
    mat.mechanical.youngs_modulus = GPaToPa(1.0);
    mat.mechanical.poisson_ratio = 0.30;
    mat.mechanical.tensile_strength = MPaToPa(0.3);
    mat.mechanical.compressive_strength = MPaToPa(3.0);
    mat.mechanical.cohesion = MPaToPa(0.5);
    mat.mechanical.friction_angle = 28.0;
    mat.mechanical.fracture_toughness_KIc = 0.05e6;
    mat.mechanical.biot_coefficient = 0.95;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.45;
    mat.hydraulic.setIsotropicPermeability(mDToM2(10.0));
    
    mat.thermal.thermal_conductivity = 0.8;
    mat.thermal.specific_heat = 850.0;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 15.0;
    mat.seismic.Qs = 8.0;
    
    mat.is_porous = true;
    mat.recommended_model = "POROELASTIC";
    
    return mat;
}

MaterialDefinition Turbidite() {
    MaterialDefinition mat;
    mat.name = "Turbidite";
    mat.category = "Sedimentary";
    mat.subcategory = "Clastic";
    mat.description = "Deep-marine graded sand/mud sequences";
    mat.reference = "Bouma, 1962";
    
    mat.mechanical.density = 2350.0;
    mat.mechanical.youngs_modulus = GPaToPa(18.0);
    mat.mechanical.poisson_ratio = 0.24;
    mat.mechanical.tensile_strength = MPaToPa(3.5);
    mat.mechanical.compressive_strength = MPaToPa(45.0);
    mat.mechanical.cohesion = MPaToPa(7.0);
    mat.mechanical.friction_angle = 32.0;
    mat.mechanical.biot_coefficient = 0.78;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.22;
    mat.hydraulic.setIsotropicPermeability(mDToM2(20.0));
    
    mat.thermal.thermal_conductivity = 2.0;
    mat.thermal.specific_heat = 900.0;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 35.0;
    mat.seismic.Qs = 17.0;
    
    mat.is_porous = true;
    mat.is_anisotropic_flag = true;
    mat.recommended_model = "POROELASTIC";
    
    return mat;
}

// --- CARBONATE ROCKS ---

MaterialDefinition Limestone() {
    MaterialDefinition mat;
    mat.name = "Limestone";
    mat.category = "Sedimentary";
    mat.subcategory = "Carbonate";
    mat.description = "Calcium carbonate rock, various textures";
    mat.reference = "Jaeger et al., 2007; Mavko et al., 2009";
    
    mat.mechanical.density = 2700.0;
    mat.mechanical.youngs_modulus = GPaToPa(50.0);
    mat.mechanical.poisson_ratio = 0.27;
    mat.mechanical.tensile_strength = MPaToPa(8.0);
    mat.mechanical.compressive_strength = MPaToPa(120.0);
    mat.mechanical.cohesion = MPaToPa(20.0);
    mat.mechanical.friction_angle = 40.0;
    mat.mechanical.fracture_toughness_KIc = 1.0e6;
    mat.mechanical.biot_coefficient = 0.60;
    mat.mechanical.grain_bulk_modulus = GPaToPa(77.0);
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.10;
    mat.hydraulic.setIsotropicPermeability(mDToM2(1.0));
    mat.hydraulic.compressibility = 2.0e-10;
    
    mat.thermal.thermal_conductivity = 2.8;
    mat.thermal.specific_heat = 880.0;
    mat.thermal.thermal_expansion = 8.0e-6;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 80.0;
    mat.seismic.Qs = 40.0;
    
    mat.is_porous = true;
    mat.recommended_model = "POROELASTIC";
    
    return mat;
}

MaterialDefinition LimestoneOolitic() {
    MaterialDefinition mat = Limestone();
    mat.name = "Limestone_Oolitic";
    mat.description = "Oolitic limestone with spherical carbonate grains";
    
    mat.mechanical.density = 2550.0;
    mat.mechanical.youngs_modulus = GPaToPa(35.0);
    mat.mechanical.compressive_strength = MPaToPa(80.0);
    mat.mechanical.biot_coefficient = 0.70;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.18;
    mat.hydraulic.setIsotropicPermeability(mDToM2(50.0));
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    
    return mat;
}

MaterialDefinition LimestoneMicritic() {
    MaterialDefinition mat = Limestone();
    mat.name = "Limestone_Micritic";
    mat.description = "Fine-grained micritic limestone (mudstone)";
    
    mat.mechanical.density = 2710.0;
    mat.mechanical.youngs_modulus = GPaToPa(55.0);
    mat.mechanical.compressive_strength = MPaToPa(130.0);
    mat.mechanical.biot_coefficient = 0.55;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.05;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.01));
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    
    return mat;
}

MaterialDefinition LimestoneBioclastic() {
    MaterialDefinition mat = Limestone();
    mat.name = "Limestone_Bioclastic";
    mat.description = "Bioclastic limestone with fossil fragments";
    
    mat.mechanical.density = 2600.0;
    mat.mechanical.youngs_modulus = GPaToPa(40.0);
    mat.mechanical.compressive_strength = MPaToPa(90.0);
    mat.mechanical.biot_coefficient = 0.68;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.15;
    mat.hydraulic.setIsotropicPermeability(mDToM2(20.0));
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    
    return mat;
}

MaterialDefinition LimestoneChalk() {
    MaterialDefinition mat;
    mat.name = "Chalk";
    mat.category = "Sedimentary";
    mat.subcategory = "Carbonate";
    mat.description = "Soft, fine-grained carbonate from coccoliths";
    mat.reference = "Andersen, 1995; Risnes, 2001";
    
    mat.mechanical.density = 2200.0;
    mat.mechanical.youngs_modulus = GPaToPa(5.0);
    mat.mechanical.poisson_ratio = 0.30;
    mat.mechanical.tensile_strength = MPaToPa(1.0);
    mat.mechanical.compressive_strength = MPaToPa(15.0);
    mat.mechanical.cohesion = MPaToPa(2.0);
    mat.mechanical.friction_angle = 30.0;
    mat.mechanical.fracture_toughness_KIc = 0.15e6;
    mat.mechanical.biot_coefficient = 0.90;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.35;
    mat.hydraulic.setIsotropicPermeability(mDToM2(2.0));
    mat.hydraulic.compressibility = 1.0e-9;
    
    mat.thermal.thermal_conductivity = 1.5;
    mat.thermal.specific_heat = 900.0;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 25.0;
    mat.seismic.Qs = 12.0;
    
    mat.is_porous = true;
    mat.recommended_model = "POROELASTIC";
    
    return mat;
}

MaterialDefinition LimestoneReef() {
    MaterialDefinition mat = Limestone();
    mat.name = "Limestone_Reef";
    mat.description = "Reef limestone with coral and algal framework";
    
    mat.mechanical.density = 2500.0;
    mat.mechanical.youngs_modulus = GPaToPa(30.0);
    mat.mechanical.compressive_strength = MPaToPa(70.0);
    mat.mechanical.biot_coefficient = 0.75;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.25;
    mat.hydraulic.setIsotropicPermeability(mDToM2(500.0));
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    
    return mat;
}

MaterialDefinition LimestoneTravertine() {
    MaterialDefinition mat = Limestone();
    mat.name = "Travertine";
    mat.description = "Chemical/biogenic limestone from hot springs";
    
    mat.mechanical.density = 2400.0;
    mat.mechanical.youngs_modulus = GPaToPa(25.0);
    mat.mechanical.compressive_strength = MPaToPa(60.0);
    mat.mechanical.biot_coefficient = 0.72;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.20;
    mat.hydraulic.setIsotropicPermeability(mDToM2(100.0));
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    
    return mat;
}

MaterialDefinition Dolomite() {
    MaterialDefinition mat;
    mat.name = "Dolomite";
    mat.category = "Sedimentary";
    mat.subcategory = "Carbonate";
    mat.description = "Calcium-magnesium carbonate rock";
    mat.reference = "Jaeger et al., 2007";
    
    mat.mechanical.density = 2850.0;
    mat.mechanical.youngs_modulus = GPaToPa(60.0);
    mat.mechanical.poisson_ratio = 0.25;
    mat.mechanical.tensile_strength = MPaToPa(10.0);
    mat.mechanical.compressive_strength = MPaToPa(150.0);
    mat.mechanical.cohesion = MPaToPa(25.0);
    mat.mechanical.friction_angle = 42.0;
    mat.mechanical.fracture_toughness_KIc = 1.2e6;
    mat.mechanical.biot_coefficient = 0.55;
    mat.mechanical.grain_bulk_modulus = GPaToPa(95.0);
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.08;
    mat.hydraulic.setIsotropicPermeability(mDToM2(10.0));
    
    mat.thermal.thermal_conductivity = 3.2;
    mat.thermal.specific_heat = 870.0;
    mat.thermal.thermal_expansion = 7.0e-6;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 100.0;
    mat.seismic.Qs = 50.0;
    
    mat.is_porous = true;
    mat.recommended_model = "POROELASTIC";
    
    return mat;
}

MaterialDefinition DolomiteSucrosic() {
    MaterialDefinition mat = Dolomite();
    mat.name = "Dolomite_Sucrosic";
    mat.description = "Sucrosic dolomite with good porosity/permeability";
    
    mat.mechanical.density = 2700.0;
    mat.mechanical.youngs_modulus = GPaToPa(45.0);
    mat.mechanical.compressive_strength = MPaToPa(100.0);
    mat.mechanical.biot_coefficient = 0.68;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.18;
    mat.hydraulic.setIsotropicPermeability(mDToM2(100.0));
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    
    return mat;
}

MaterialDefinition Marl() {
    MaterialDefinition mat;
    mat.name = "Marl";
    mat.category = "Sedimentary";
    mat.subcategory = "Carbonate";
    mat.description = "Calcareous mudstone, mix of clay and carbonate";
    mat.reference = "Jaeger et al., 2007";
    
    mat.mechanical.density = 2400.0;
    mat.mechanical.youngs_modulus = GPaToPa(15.0);
    mat.mechanical.poisson_ratio = 0.28;
    mat.mechanical.tensile_strength = MPaToPa(3.0);
    mat.mechanical.compressive_strength = MPaToPa(35.0);
    mat.mechanical.cohesion = MPaToPa(5.0);
    mat.mechanical.friction_angle = 28.0;
    mat.mechanical.fracture_toughness_KIc = 0.4e6;
    mat.mechanical.biot_coefficient = 0.80;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.25;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.1));
    
    mat.thermal.thermal_conductivity = 1.8;
    mat.thermal.specific_heat = 880.0;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 30.0;
    mat.seismic.Qs = 15.0;
    
    mat.is_porous = true;
    mat.recommended_model = "POROELASTIC";
    
    return mat;
}

MaterialDefinition Coquina() {
    MaterialDefinition mat;
    mat.name = "Coquina";
    mat.category = "Sedimentary";
    mat.subcategory = "Carbonate";
    mat.description = "Shell debris limestone, high porosity";
    mat.reference = "Tucker & Wright, 1990";
    
    mat.mechanical.density = 2300.0;
    mat.mechanical.youngs_modulus = GPaToPa(20.0);
    mat.mechanical.poisson_ratio = 0.26;
    mat.mechanical.tensile_strength = MPaToPa(3.5);
    mat.mechanical.compressive_strength = MPaToPa(40.0);
    mat.mechanical.cohesion = MPaToPa(6.0);
    mat.mechanical.friction_angle = 35.0;
    mat.mechanical.biot_coefficient = 0.80;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.30;
    mat.hydraulic.setIsotropicPermeability(mDToM2(1000.0));
    
    mat.thermal.thermal_conductivity = 2.0;
    mat.thermal.specific_heat = 880.0;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 30.0;
    mat.seismic.Qs = 15.0;
    
    mat.is_porous = true;
    mat.recommended_model = "POROELASTIC";
    
    return mat;
}

// --- CHEMICAL/EVAPORITE ROCKS ---

MaterialDefinition Halite() {
    MaterialDefinition mat;
    mat.name = "Halite";
    mat.category = "Sedimentary";
    mat.subcategory = "Evaporite";
    mat.description = "Rock salt, NaCl";
    mat.reference = "Jaeger et al., 2007; Carter & Hansen, 1983";
    
    mat.mechanical.density = 2160.0;
    mat.mechanical.youngs_modulus = GPaToPa(30.0);
    mat.mechanical.poisson_ratio = 0.25;
    mat.mechanical.tensile_strength = MPaToPa(2.0);
    mat.mechanical.compressive_strength = MPaToPa(25.0);
    mat.mechanical.cohesion = MPaToPa(1.5);
    mat.mechanical.friction_angle = 30.0;
    mat.mechanical.fracture_toughness_KIc = 0.2e6;
    mat.mechanical.biot_coefficient = 0.05;
    mat.mechanical.viscosity = 1.0e18;
    mat.mechanical.relaxation_time = 1.0e6;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.005;
    mat.hydraulic.setIsotropicPermeability(mDToM2(1e-6));
    
    mat.thermal.thermal_conductivity = 5.5;
    mat.thermal.specific_heat = 920.0;
    mat.thermal.thermal_expansion = 4.0e-5;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 300.0;
    mat.seismic.Qs = 150.0;
    
    mat.is_porous = false;
    mat.is_viscoelastic = true;
    mat.recommended_model = "VISCOELASTIC_SLS";
    
    return mat;
}

MaterialDefinition Gypsum() {
    MaterialDefinition mat;
    mat.name = "Gypsum";
    mat.category = "Sedimentary";
    mat.subcategory = "Evaporite";
    mat.description = "Hydrated calcium sulfate, CaSO4·2H2O";
    mat.reference = "Jaeger et al., 2007";
    
    mat.mechanical.density = 2300.0;
    mat.mechanical.youngs_modulus = GPaToPa(20.0);
    mat.mechanical.poisson_ratio = 0.28;
    mat.mechanical.tensile_strength = MPaToPa(2.5);
    mat.mechanical.compressive_strength = MPaToPa(30.0);
    mat.mechanical.cohesion = MPaToPa(3.0);
    mat.mechanical.friction_angle = 32.0;
    mat.mechanical.fracture_toughness_KIc = 0.25e6;
    mat.mechanical.biot_coefficient = 0.30;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.02;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.001));
    
    mat.thermal.thermal_conductivity = 1.3;
    mat.thermal.specific_heat = 1090.0;
    mat.thermal.thermal_expansion = 2.5e-5;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 100.0;
    mat.seismic.Qs = 50.0;
    
    mat.is_porous = false;
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition Anhydrite() {
    MaterialDefinition mat;
    mat.name = "Anhydrite";
    mat.category = "Sedimentary";
    mat.subcategory = "Evaporite";
    mat.description = "Anhydrous calcium sulfate, CaSO4";
    mat.reference = "Jaeger et al., 2007";
    
    mat.mechanical.density = 2960.0;
    mat.mechanical.youngs_modulus = GPaToPa(70.0);
    mat.mechanical.poisson_ratio = 0.30;
    mat.mechanical.tensile_strength = MPaToPa(8.0);
    mat.mechanical.compressive_strength = MPaToPa(100.0);
    mat.mechanical.cohesion = MPaToPa(15.0);
    mat.mechanical.friction_angle = 35.0;
    mat.mechanical.fracture_toughness_KIc = 0.8e6;
    mat.mechanical.biot_coefficient = 0.20;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.005;
    mat.hydraulic.setIsotropicPermeability(mDToM2(1e-5));
    
    mat.thermal.thermal_conductivity = 4.8;
    mat.thermal.specific_heat = 750.0;
    mat.thermal.thermal_expansion = 2.0e-5;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 200.0;
    mat.seismic.Qs = 100.0;
    
    mat.is_porous = false;
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition Sylvite() {
    MaterialDefinition mat;
    mat.name = "Sylvite";
    mat.category = "Sedimentary";
    mat.subcategory = "Evaporite";
    mat.description = "Potassium chloride, KCl";
    mat.reference = "Jaeger et al., 2007";
    
    mat.mechanical.density = 1990.0;
    mat.mechanical.youngs_modulus = GPaToPa(25.0);
    mat.mechanical.poisson_ratio = 0.26;
    mat.mechanical.tensile_strength = MPaToPa(1.5);
    mat.mechanical.compressive_strength = MPaToPa(20.0);
    mat.mechanical.cohesion = MPaToPa(1.2);
    mat.mechanical.friction_angle = 28.0;
    mat.mechanical.biot_coefficient = 0.08;
    mat.mechanical.viscosity = 5.0e17;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.003;
    mat.hydraulic.setIsotropicPermeability(mDToM2(1e-7));
    
    mat.thermal.thermal_conductivity = 6.5;
    mat.thermal.specific_heat = 690.0;
    mat.thermal.thermal_expansion = 3.8e-5;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 280.0;
    mat.seismic.Qs = 140.0;
    
    mat.is_porous = false;
    mat.is_viscoelastic = true;
    mat.recommended_model = "VISCOELASTIC_SLS";
    
    return mat;
}

MaterialDefinition Potash() {
    MaterialDefinition mat;
    mat.name = "Potash";
    mat.category = "Sedimentary";
    mat.subcategory = "Evaporite";
    mat.description = "Potassium-bearing evaporite ore";
    mat.reference = "Jaeger et al., 2007";
    
    mat.mechanical.density = 2000.0;
    mat.mechanical.youngs_modulus = GPaToPa(22.0);
    mat.mechanical.poisson_ratio = 0.27;
    mat.mechanical.tensile_strength = MPaToPa(1.2);
    mat.mechanical.compressive_strength = MPaToPa(18.0);
    mat.mechanical.cohesion = MPaToPa(1.0);
    mat.mechanical.friction_angle = 26.0;
    mat.mechanical.biot_coefficient = 0.10;
    mat.mechanical.viscosity = 3.0e17;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.005;
    mat.hydraulic.setIsotropicPermeability(mDToM2(1e-6));
    
    mat.thermal.thermal_conductivity = 5.0;
    mat.thermal.specific_heat = 700.0;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 250.0;
    mat.seismic.Qs = 125.0;
    
    mat.is_porous = false;
    mat.is_viscoelastic = true;
    mat.recommended_model = "VISCOELASTIC_SLS";
    
    return mat;
}

MaterialDefinition Trona() {
    MaterialDefinition mat;
    mat.name = "Trona";
    mat.category = "Sedimentary";
    mat.subcategory = "Evaporite";
    mat.description = "Sodium sesquicarbonate, Na3(CO3)(HCO3)·2H2O";
    mat.reference = "Jaeger et al., 2007";
    
    mat.mechanical.density = 2140.0;
    mat.mechanical.youngs_modulus = GPaToPa(15.0);
    mat.mechanical.poisson_ratio = 0.28;
    mat.mechanical.tensile_strength = MPaToPa(1.8);
    mat.mechanical.compressive_strength = MPaToPa(22.0);
    mat.mechanical.cohesion = MPaToPa(2.0);
    mat.mechanical.friction_angle = 30.0;
    mat.mechanical.biot_coefficient = 0.15;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.008;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.0001));
    
    mat.thermal.thermal_conductivity = 1.8;
    mat.thermal.specific_heat = 1050.0;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 150.0;
    mat.seismic.Qs = 75.0;
    
    mat.is_porous = false;
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition Chert() {
    MaterialDefinition mat;
    mat.name = "Chert";
    mat.category = "Sedimentary";
    mat.subcategory = "Chemical";
    mat.description = "Microcrystalline/cryptocrystalline silica";
    mat.reference = "Jaeger et al., 2007";
    
    mat.mechanical.density = 2600.0;
    mat.mechanical.youngs_modulus = GPaToPa(80.0);
    mat.mechanical.poisson_ratio = 0.15;
    mat.mechanical.tensile_strength = MPaToPa(5.0);
    mat.mechanical.compressive_strength = MPaToPa(200.0);
    mat.mechanical.cohesion = MPaToPa(30.0);
    mat.mechanical.friction_angle = 50.0;
    mat.mechanical.fracture_toughness_KIc = 0.7e6;
    mat.mechanical.biot_coefficient = 0.20;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.02;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.0001));
    
    mat.thermal.thermal_conductivity = 4.0;
    mat.thermal.specific_heat = 740.0;
    mat.thermal.thermal_expansion = 1.0e-5;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 300.0;
    mat.seismic.Qs = 150.0;
    
    mat.is_porous = false;
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition Flint() {
    MaterialDefinition mat = Chert();
    mat.name = "Flint";
    mat.description = "Nodular silica in chalk, very fine-grained";
    
    mat.mechanical.density = 2650.0;
    mat.mechanical.youngs_modulus = GPaToPa(85.0);
    mat.mechanical.compressive_strength = MPaToPa(220.0);
    mat.mechanical.fracture_toughness_KIc = 0.65e6;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.01;
    mat.hydraulic.setIsotropicPermeability(mDToM2(1e-6));
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    
    return mat;
}

MaterialDefinition Diatomite() {
    MaterialDefinition mat;
    mat.name = "Diatomite";
    mat.category = "Sedimentary";
    mat.subcategory = "Chemical";
    mat.description = "Biogenic silica from diatom tests";
    mat.reference = "Jaeger et al., 2007";
    
    mat.mechanical.density = 600.0;
    mat.mechanical.youngs_modulus = GPaToPa(0.5);
    mat.mechanical.poisson_ratio = 0.25;
    mat.mechanical.tensile_strength = MPaToPa(0.2);
    mat.mechanical.compressive_strength = MPaToPa(2.0);
    mat.mechanical.cohesion = MPaToPa(0.3);
    mat.mechanical.friction_angle = 30.0;
    mat.mechanical.biot_coefficient = 0.95;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.70;
    mat.hydraulic.setIsotropicPermeability(mDToM2(100.0));
    
    mat.thermal.thermal_conductivity = 0.15;
    mat.thermal.specific_heat = 900.0;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 10.0;
    mat.seismic.Qs = 5.0;
    
    mat.is_porous = true;
    mat.recommended_model = "POROELASTIC";
    
    return mat;
}

MaterialDefinition Ironstone() {
    MaterialDefinition mat;
    mat.name = "Ironstone";
    mat.category = "Sedimentary";
    mat.subcategory = "Chemical";
    mat.description = "Iron-rich sedimentary rock";
    mat.reference = "Jaeger et al., 2007";
    
    mat.mechanical.density = 3200.0;
    mat.mechanical.youngs_modulus = GPaToPa(50.0);
    mat.mechanical.poisson_ratio = 0.25;
    mat.mechanical.tensile_strength = MPaToPa(8.0);
    mat.mechanical.compressive_strength = MPaToPa(100.0);
    mat.mechanical.cohesion = MPaToPa(18.0);
    mat.mechanical.friction_angle = 38.0;
    mat.mechanical.fracture_toughness_KIc = 1.0e6;
    mat.mechanical.biot_coefficient = 0.50;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.10;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.5));
    
    mat.thermal.thermal_conductivity = 3.5;
    mat.thermal.specific_heat = 650.0;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 80.0;
    mat.seismic.Qs = 40.0;
    
    mat.is_porous = true;
    mat.recommended_model = "POROELASTIC";
    
    return mat;
}

MaterialDefinition BandedIronFormation() {
    MaterialDefinition mat;
    mat.name = "BIF";
    mat.category = "Sedimentary";
    mat.subcategory = "Chemical";
    mat.description = "Precambrian banded iron formation";
    mat.reference = "Jaeger et al., 2007";
    
    mat.mechanical.density = 3500.0;
    mat.mechanical.youngs_modulus = GPaToPa(100.0);
    mat.mechanical.poisson_ratio = 0.22;
    mat.mechanical.tensile_strength = MPaToPa(15.0);
    mat.mechanical.compressive_strength = MPaToPa(200.0);
    mat.mechanical.cohesion = MPaToPa(35.0);
    mat.mechanical.friction_angle = 45.0;
    mat.mechanical.fracture_toughness_KIc = 2.0e6;
    mat.mechanical.biot_coefficient = 0.25;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.02;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.001));
    
    mat.thermal.thermal_conductivity = 5.0;
    mat.thermal.specific_heat = 700.0;
    
    mat.anisotropic.is_anisotropic = true;
    mat.anisotropic.symmetry_type = "VTI";
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 150.0;
    mat.seismic.Qs = 75.0;
    
    mat.is_porous = false;
    mat.is_anisotropic_flag = true;
    mat.recommended_model = "ANISOTROPIC_VTI";
    
    return mat;
}

MaterialDefinition Phosphorite() {
    MaterialDefinition mat;
    mat.name = "Phosphorite";
    mat.category = "Sedimentary";
    mat.subcategory = "Chemical";
    mat.description = "Phosphate-rich sedimentary rock";
    mat.reference = "Jaeger et al., 2007";
    
    mat.mechanical.density = 2700.0;
    mat.mechanical.youngs_modulus = GPaToPa(35.0);
    mat.mechanical.poisson_ratio = 0.26;
    mat.mechanical.tensile_strength = MPaToPa(5.0);
    mat.mechanical.compressive_strength = MPaToPa(70.0);
    mat.mechanical.cohesion = MPaToPa(12.0);
    mat.mechanical.friction_angle = 35.0;
    mat.mechanical.biot_coefficient = 0.60;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.15;
    mat.hydraulic.setIsotropicPermeability(mDToM2(1.0));
    
    mat.thermal.thermal_conductivity = 2.0;
    mat.thermal.specific_heat = 800.0;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 50.0;
    mat.seismic.Qs = 25.0;
    
    mat.is_porous = true;
    mat.recommended_model = "POROELASTIC";
    
    return mat;
}

// --- ORGANIC ROCKS ---

MaterialDefinition Coal() {
    MaterialDefinition mat;
    mat.name = "Coal";
    mat.category = "Sedimentary";
    mat.subcategory = "Organic";
    mat.description = "Combustible sedimentary rock, general";
    mat.reference = "Jaeger et al., 2007";
    
    mat.mechanical.density = 1400.0;
    mat.mechanical.youngs_modulus = GPaToPa(4.0);
    mat.mechanical.poisson_ratio = 0.32;
    mat.mechanical.tensile_strength = MPaToPa(1.0);
    mat.mechanical.compressive_strength = MPaToPa(15.0);
    mat.mechanical.cohesion = MPaToPa(2.0);
    mat.mechanical.friction_angle = 30.0;
    mat.mechanical.fracture_toughness_KIc = 0.2e6;
    mat.mechanical.biot_coefficient = 0.80;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.10;
    mat.hydraulic.setIsotropicPermeability(mDToM2(1.0));
    
    mat.thermal.thermal_conductivity = 0.3;
    mat.thermal.specific_heat = 1300.0;
    mat.thermal.thermal_expansion = 3.0e-5;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 15.0;
    mat.seismic.Qs = 8.0;
    
    mat.is_porous = true;
    mat.recommended_model = "POROELASTIC";
    
    return mat;
}

MaterialDefinition CoalAnthracite() {
    MaterialDefinition mat = Coal();
    mat.name = "Coal_Anthracite";
    mat.description = "Highest rank coal, >86% carbon";
    
    mat.mechanical.density = 1550.0;
    mat.mechanical.youngs_modulus = GPaToPa(8.0);
    mat.mechanical.compressive_strength = MPaToPa(30.0);
    mat.mechanical.tensile_strength = MPaToPa(2.5);
    mat.mechanical.biot_coefficient = 0.60;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.03;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.01));
    
    mat.thermal.thermal_conductivity = 0.5;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    
    return mat;
}

MaterialDefinition CoalBituminous() {
    MaterialDefinition mat = Coal();
    mat.name = "Coal_Bituminous";
    mat.description = "Medium rank coal, 45-86% carbon";
    
    mat.mechanical.density = 1350.0;
    mat.mechanical.youngs_modulus = GPaToPa(3.5);
    mat.mechanical.compressive_strength = MPaToPa(12.0);
    mat.mechanical.biot_coefficient = 0.85;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.08;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.5));
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    
    return mat;
}

MaterialDefinition CoalLignite() {
    MaterialDefinition mat = Coal();
    mat.name = "Coal_Lignite";
    mat.description = "Lowest rank coal, brown coal";
    
    mat.mechanical.density = 1200.0;
    mat.mechanical.youngs_modulus = GPaToPa(1.5);
    mat.mechanical.compressive_strength = MPaToPa(5.0);
    mat.mechanical.tensile_strength = MPaToPa(0.5);
    mat.mechanical.biot_coefficient = 0.92;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.25;
    mat.hydraulic.setIsotropicPermeability(mDToM2(5.0));
    
    mat.thermal.thermal_conductivity = 0.2;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    
    return mat;
}

MaterialDefinition Peat() {
    MaterialDefinition mat;
    mat.name = "Peat";
    mat.category = "Sedimentary";
    mat.subcategory = "Organic";
    mat.description = "Partially decomposed organic matter";
    mat.reference = "Hobbs, 1986";
    
    mat.mechanical.density = 1100.0;
    mat.mechanical.youngs_modulus = GPaToPa(0.1);
    mat.mechanical.poisson_ratio = 0.35;
    mat.mechanical.tensile_strength = MPaToPa(0.05);
    mat.mechanical.compressive_strength = MPaToPa(0.5);
    mat.mechanical.cohesion = MPaToPa(0.1);
    mat.mechanical.friction_angle = 25.0;
    mat.mechanical.biot_coefficient = 0.98;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.85;
    mat.hydraulic.setIsotropicPermeability(mDToM2(1000.0));
    
    mat.thermal.thermal_conductivity = 0.1;
    mat.thermal.specific_heat = 3500.0;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 5.0;
    mat.seismic.Qs = 2.0;
    
    mat.is_porous = true;
    mat.recommended_model = "POROELASTIC";
    
    return mat;
}

MaterialDefinition OilShale() {
    MaterialDefinition mat;
    mat.name = "Oil_Shale";
    mat.category = "Sedimentary";
    mat.subcategory = "Organic";
    mat.description = "Kerogen-rich fine-grained rock";
    mat.reference = "Sone & Zoback, 2013";
    
    mat.mechanical.density = 2200.0;
    mat.mechanical.youngs_modulus = GPaToPa(10.0);
    mat.mechanical.poisson_ratio = 0.28;
    mat.mechanical.tensile_strength = MPaToPa(2.0);
    mat.mechanical.compressive_strength = MPaToPa(25.0);
    mat.mechanical.cohesion = MPaToPa(4.0);
    mat.mechanical.friction_angle = 25.0;
    mat.mechanical.fracture_toughness_KIc = 0.3e6;
    mat.mechanical.biot_coefficient = 0.85;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.12;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.0001));
    
    mat.thermal.thermal_conductivity = 1.0;
    mat.thermal.specific_heat = 1200.0;
    mat.thermal.heat_production = 5.0e-6;
    
    mat.anisotropic.is_anisotropic = true;
    mat.anisotropic.symmetry_type = "VTI";
    mat.anisotropic.epsilon = 0.20;
    mat.anisotropic.delta = 0.12;
    mat.anisotropic.gamma = 0.25;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 20.0;
    mat.seismic.Qs = 10.0;
    
    mat.is_porous = true;
    mat.is_anisotropic_flag = true;
    mat.recommended_model = "ANISOTROPIC_VTI";
    
    return mat;
}

} // namespace Sedimentary

// =============================================================================
// Metamorphic Rock Definitions - namespace Materials::Metamorphic
// =============================================================================

namespace Metamorphic {

MaterialDefinition Slate() {
    MaterialDefinition mat;
    mat.name = "Slate";
    mat.category = "Metamorphic";
    mat.subcategory = "Foliated";
    mat.description = "Very fine-grained foliated rock from shale/mudstone";
    mat.reference = "Jaeger et al., 2007";
    
    mat.mechanical.density = 2750.0;
    mat.mechanical.youngs_modulus = GPaToPa(40.0);
    mat.mechanical.poisson_ratio = 0.20;
    mat.mechanical.tensile_strength = MPaToPa(15.0);
    mat.mechanical.compressive_strength = MPaToPa(150.0);
    mat.mechanical.cohesion = MPaToPa(25.0);
    mat.mechanical.friction_angle = 45.0;
    mat.mechanical.fracture_toughness_KIc = 0.8e6;
    mat.mechanical.biot_coefficient = 0.35;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.005;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.0001));
    
    mat.thermal.thermal_conductivity = 2.0;
    mat.thermal.specific_heat = 850.0;
    mat.thermal.thermal_expansion = 9.0e-6;
    
    mat.anisotropic.is_anisotropic = true;
    mat.anisotropic.symmetry_type = "VTI";
    mat.anisotropic.E_vertical = GPaToPa(35.0);
    mat.anisotropic.E_horizontal = GPaToPa(45.0);
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 100.0;
    mat.seismic.Qs = 50.0;
    
    mat.is_anisotropic_flag = true;
    mat.recommended_model = "ANISOTROPIC_VTI";
    
    return mat;
}

MaterialDefinition Phyllite() {
    MaterialDefinition mat;
    mat.name = "Phyllite";
    mat.category = "Metamorphic";
    mat.subcategory = "Foliated";
    mat.description = "Fine-grained foliated rock between slate and schist";
    mat.reference = "Jaeger et al., 2007";
    
    mat.mechanical.density = 2750.0;
    mat.mechanical.youngs_modulus = GPaToPa(45.0);
    mat.mechanical.poisson_ratio = 0.22;
    mat.mechanical.tensile_strength = MPaToPa(12.0);
    mat.mechanical.compressive_strength = MPaToPa(140.0);
    mat.mechanical.cohesion = MPaToPa(22.0);
    mat.mechanical.friction_angle = 42.0;
    mat.mechanical.fracture_toughness_KIc = 0.9e6;
    mat.mechanical.biot_coefficient = 0.38;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.008;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.001));
    
    mat.thermal.thermal_conductivity = 2.2;
    mat.thermal.specific_heat = 860.0;
    mat.thermal.thermal_expansion = 8.5e-6;
    
    mat.anisotropic.is_anisotropic = true;
    mat.anisotropic.symmetry_type = "VTI";
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 90.0;
    mat.seismic.Qs = 45.0;
    
    mat.is_anisotropic_flag = true;
    mat.recommended_model = "ANISOTROPIC_VTI";
    
    return mat;
}

MaterialDefinition Schist() {
    MaterialDefinition mat;
    mat.name = "Schist";
    mat.category = "Metamorphic";
    mat.subcategory = "Foliated";
    mat.description = "Medium-grade foliated rock with visible minerals";
    mat.reference = "Jaeger et al., 2007";
    
    mat.mechanical.density = 2800.0;
    mat.mechanical.youngs_modulus = GPaToPa(50.0);
    mat.mechanical.poisson_ratio = 0.22;
    mat.mechanical.tensile_strength = MPaToPa(10.0);
    mat.mechanical.compressive_strength = MPaToPa(120.0);
    mat.mechanical.cohesion = MPaToPa(18.0);
    mat.mechanical.friction_angle = 40.0;
    mat.mechanical.fracture_toughness_KIc = 1.0e6;
    mat.mechanical.biot_coefficient = 0.40;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.01;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.01));
    
    mat.thermal.thermal_conductivity = 2.5;
    mat.thermal.specific_heat = 880.0;
    mat.thermal.thermal_expansion = 1.0e-5;
    
    mat.anisotropic.is_anisotropic = true;
    mat.anisotropic.symmetry_type = "VTI";
    mat.anisotropic.E_vertical = GPaToPa(40.0);
    mat.anisotropic.E_horizontal = GPaToPa(60.0);
    mat.anisotropic.epsilon = 0.10;
    mat.anisotropic.delta = 0.05;
    mat.anisotropic.gamma = 0.15;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 80.0;
    mat.seismic.Qs = 40.0;
    
    mat.is_anisotropic_flag = true;
    mat.recommended_model = "ANISOTROPIC_VTI";
    
    return mat;
}

MaterialDefinition SchistMica() {
    MaterialDefinition mat = Schist();
    mat.name = "Schist_Mica";
    mat.description = "Mica-rich schist with high foliation";
    
    mat.mechanical.density = 2750.0;
    mat.mechanical.youngs_modulus = GPaToPa(45.0);
    mat.mechanical.compressive_strength = MPaToPa(100.0);
    mat.mechanical.computeDerivedProperties();
    
    mat.anisotropic.E_vertical = GPaToPa(35.0);
    mat.anisotropic.E_horizontal = GPaToPa(55.0);
    mat.anisotropic.epsilon = 0.15;
    mat.anisotropic.gamma = 0.20;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    
    return mat;
}

MaterialDefinition SchistChlorite() {
    MaterialDefinition mat = Schist();
    mat.name = "Schist_Chlorite";
    mat.description = "Chlorite-rich schist (greenschist)";
    
    mat.mechanical.density = 2850.0;
    mat.mechanical.youngs_modulus = GPaToPa(55.0);
    mat.mechanical.compressive_strength = MPaToPa(130.0);
    mat.mechanical.computeDerivedProperties();
    
    mat.thermal.thermal_conductivity = 2.8;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    
    return mat;
}

MaterialDefinition SchistTalc() {
    MaterialDefinition mat = Schist();
    mat.name = "Schist_Talc";
    mat.description = "Talc-rich schist, very soft";
    
    mat.mechanical.density = 2700.0;
    mat.mechanical.youngs_modulus = GPaToPa(20.0);
    mat.mechanical.compressive_strength = MPaToPa(40.0);
    mat.mechanical.tensile_strength = MPaToPa(3.0);
    mat.mechanical.cohesion = MPaToPa(5.0);
    mat.mechanical.friction_angle = 30.0;
    mat.mechanical.biot_coefficient = 0.55;
    mat.mechanical.computeDerivedProperties();
    
    mat.thermal.thermal_conductivity = 3.0;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    
    return mat;
}

MaterialDefinition SchistBlueschist() {
    MaterialDefinition mat = Schist();
    mat.name = "Blueschist";
    mat.description = "High-P/low-T metamorphic rock with glaucophane";
    mat.reference = "Evans, 1990";
    
    mat.mechanical.density = 3000.0;
    mat.mechanical.youngs_modulus = GPaToPa(70.0);
    mat.mechanical.compressive_strength = MPaToPa(180.0);
    mat.mechanical.tensile_strength = MPaToPa(14.0);
    mat.mechanical.biot_coefficient = 0.30;
    mat.mechanical.computeDerivedProperties();
    
    mat.thermal.thermal_conductivity = 2.8;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 120.0;
    mat.seismic.Qs = 60.0;
    
    return mat;
}

MaterialDefinition Gneiss() {
    MaterialDefinition mat;
    mat.name = "Gneiss";
    mat.category = "Metamorphic";
    mat.subcategory = "Foliated";
    mat.description = "High-grade banded metamorphic rock";
    mat.reference = "Jaeger et al., 2007";
    
    mat.mechanical.density = 2700.0;
    mat.mechanical.youngs_modulus = GPaToPa(60.0);
    mat.mechanical.poisson_ratio = 0.24;
    mat.mechanical.tensile_strength = MPaToPa(12.0);
    mat.mechanical.compressive_strength = MPaToPa(180.0);
    mat.mechanical.cohesion = MPaToPa(30.0);
    mat.mechanical.friction_angle = 50.0;
    mat.mechanical.fracture_toughness_KIc = 1.5e6;
    mat.mechanical.biot_coefficient = 0.40;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.008;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.001));
    
    mat.thermal.thermal_conductivity = 2.8;
    mat.thermal.specific_heat = 820.0;
    mat.thermal.thermal_expansion = 8.0e-6;
    mat.thermal.heat_production = 2.0e-6;
    
    mat.anisotropic.is_anisotropic = true;
    mat.anisotropic.symmetry_type = "VTI";
    mat.anisotropic.E_vertical = GPaToPa(55.0);
    mat.anisotropic.E_horizontal = GPaToPa(65.0);
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 150.0;
    mat.seismic.Qs = 75.0;
    
    mat.is_anisotropic_flag = true;
    mat.recommended_model = "ANISOTROPIC_VTI";
    
    return mat;
}

MaterialDefinition GneissGranitic() {
    MaterialDefinition mat = Gneiss();
    mat.name = "Gneiss_Granitic";
    mat.description = "Granitic composition gneiss (orthogneiss)";
    
    mat.mechanical.density = 2680.0;
    mat.mechanical.youngs_modulus = GPaToPa(65.0);
    mat.mechanical.compressive_strength = MPaToPa(190.0);
    mat.mechanical.computeDerivedProperties();
    
    mat.thermal.heat_production = 2.5e-6;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    
    return mat;
}

MaterialDefinition GneissBanded() {
    MaterialDefinition mat = Gneiss();
    mat.name = "Gneiss_Banded";
    mat.description = "Strongly banded gneiss with high anisotropy";
    
    mat.anisotropic.E_vertical = GPaToPa(50.0);
    mat.anisotropic.E_horizontal = GPaToPa(70.0);
    mat.anisotropic.epsilon = 0.12;
    mat.anisotropic.delta = 0.06;
    mat.anisotropic.gamma = 0.18;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    
    return mat;
}

MaterialDefinition Migmatite() {
    MaterialDefinition mat;
    mat.name = "Migmatite";
    mat.category = "Metamorphic";
    mat.subcategory = "Foliated";
    mat.description = "Partially melted gneiss with leucosome/melanosome";
    mat.reference = "Jaeger et al., 2007";
    
    mat.mechanical.density = 2720.0;
    mat.mechanical.youngs_modulus = GPaToPa(55.0);
    mat.mechanical.poisson_ratio = 0.24;
    mat.mechanical.tensile_strength = MPaToPa(11.0);
    mat.mechanical.compressive_strength = MPaToPa(165.0);
    mat.mechanical.cohesion = MPaToPa(28.0);
    mat.mechanical.friction_angle = 48.0;
    mat.mechanical.fracture_toughness_KIc = 1.4e6;
    mat.mechanical.biot_coefficient = 0.42;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.01;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.005));
    
    mat.thermal.thermal_conductivity = 2.7;
    mat.thermal.specific_heat = 830.0;
    mat.thermal.heat_production = 2.2e-6;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 140.0;
    mat.seismic.Qs = 70.0;
    
    mat.is_anisotropic_flag = true;
    mat.recommended_model = "ANISOTROPIC_VTI";
    
    return mat;
}

MaterialDefinition Mylonite() {
    MaterialDefinition mat;
    mat.name = "Mylonite";
    mat.category = "Metamorphic";
    mat.subcategory = "Foliated";
    mat.description = "Fault rock formed by ductile shearing";
    mat.reference = "Sibson, 1977";
    
    mat.mechanical.density = 2700.0;
    mat.mechanical.youngs_modulus = GPaToPa(40.0);
    mat.mechanical.poisson_ratio = 0.25;
    mat.mechanical.tensile_strength = MPaToPa(6.0);
    mat.mechanical.compressive_strength = MPaToPa(100.0);
    mat.mechanical.cohesion = MPaToPa(12.0);
    mat.mechanical.friction_angle = 35.0;
    mat.mechanical.fracture_toughness_KIc = 0.8e6;
    mat.mechanical.biot_coefficient = 0.55;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.02;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.1));
    
    mat.thermal.thermal_conductivity = 2.4;
    mat.thermal.specific_heat = 850.0;
    
    mat.anisotropic.is_anisotropic = true;
    mat.anisotropic.symmetry_type = "VTI";
    mat.anisotropic.epsilon = 0.20;
    mat.anisotropic.delta = 0.10;
    mat.anisotropic.gamma = 0.25;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 60.0;
    mat.seismic.Qs = 30.0;
    
    mat.is_fractured = true;
    mat.is_anisotropic_flag = true;
    mat.recommended_model = "ANISOTROPIC_VTI";
    
    return mat;
}

MaterialDefinition Cataclasite() {
    MaterialDefinition mat;
    mat.name = "Cataclasite";
    mat.category = "Metamorphic";
    mat.subcategory = "Foliated";
    mat.description = "Fault rock formed by brittle deformation";
    mat.reference = "Sibson, 1977";
    
    mat.mechanical.density = 2600.0;
    mat.mechanical.youngs_modulus = GPaToPa(25.0);
    mat.mechanical.poisson_ratio = 0.26;
    mat.mechanical.tensile_strength = MPaToPa(3.0);
    mat.mechanical.compressive_strength = MPaToPa(50.0);
    mat.mechanical.cohesion = MPaToPa(6.0);
    mat.mechanical.friction_angle = 32.0;
    mat.mechanical.fracture_toughness_KIc = 0.5e6;
    mat.mechanical.biot_coefficient = 0.70;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.08;
    mat.hydraulic.setIsotropicPermeability(mDToM2(10.0));
    
    mat.thermal.thermal_conductivity = 2.0;
    mat.thermal.specific_heat = 870.0;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 40.0;
    mat.seismic.Qs = 20.0;
    
    mat.is_porous = true;
    mat.is_fractured = true;
    mat.recommended_model = "POROELASTIC";
    
    return mat;
}

MaterialDefinition Amphibolite() {
    MaterialDefinition mat;
    mat.name = "Amphibolite";
    mat.category = "Metamorphic";
    mat.subcategory = "Foliated";
    mat.description = "Metamorphosed mafic rock with amphibole + plagioclase";
    mat.reference = "Jaeger et al., 2007";
    
    mat.mechanical.density = 2950.0;
    mat.mechanical.youngs_modulus = GPaToPa(80.0);
    mat.mechanical.poisson_ratio = 0.26;
    mat.mechanical.tensile_strength = MPaToPa(15.0);
    mat.mechanical.compressive_strength = MPaToPa(220.0);
    mat.mechanical.cohesion = MPaToPa(40.0);
    mat.mechanical.friction_angle = 48.0;
    mat.mechanical.fracture_toughness_KIc = 2.0e6;
    mat.mechanical.biot_coefficient = 0.32;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.005;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.0001));
    
    mat.thermal.thermal_conductivity = 2.6;
    mat.thermal.specific_heat = 880.0;
    mat.thermal.thermal_expansion = 7.0e-6;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 180.0;
    mat.seismic.Qs = 90.0;
    
    mat.is_anisotropic_flag = true;
    mat.recommended_model = "ANISOTROPIC_VTI";
    
    return mat;
}

// --- NON-FOLIATED ---

MaterialDefinition Marble() {
    MaterialDefinition mat;
    mat.name = "Marble";
    mat.category = "Metamorphic";
    mat.subcategory = "Non-foliated";
    mat.description = "Metamorphosed limestone/dolomite";
    mat.reference = "Jaeger et al., 2007";
    
    mat.mechanical.density = 2700.0;
    mat.mechanical.youngs_modulus = GPaToPa(55.0);
    mat.mechanical.poisson_ratio = 0.28;
    mat.mechanical.tensile_strength = MPaToPa(10.0);
    mat.mechanical.compressive_strength = MPaToPa(140.0);
    mat.mechanical.cohesion = MPaToPa(25.0);
    mat.mechanical.friction_angle = 42.0;
    mat.mechanical.fracture_toughness_KIc = 1.2e6;
    mat.mechanical.biot_coefficient = 0.50;
    mat.mechanical.grain_bulk_modulus = GPaToPa(77.0);
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.01;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.001));
    
    mat.thermal.thermal_conductivity = 2.7;
    mat.thermal.specific_heat = 880.0;
    mat.thermal.thermal_expansion = 1.1e-5;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 120.0;
    mat.seismic.Qs = 60.0;
    
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition MarbleCalcite() {
    MaterialDefinition mat = Marble();
    mat.name = "Marble_Calcite";
    mat.description = "Calcite marble (metamorphosed limestone)";
    
    mat.mechanical.density = 2700.0;
    mat.mechanical.compressive_strength = MPaToPa(130.0);
    mat.mechanical.computeDerivedProperties();
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    
    return mat;
}

MaterialDefinition MarbleDolomitic() {
    MaterialDefinition mat = Marble();
    mat.name = "Marble_Dolomitic";
    mat.description = "Dolomitic marble (metamorphosed dolomite)";
    
    mat.mechanical.density = 2850.0;
    mat.mechanical.youngs_modulus = GPaToPa(65.0);
    mat.mechanical.compressive_strength = MPaToPa(160.0);
    mat.mechanical.grain_bulk_modulus = GPaToPa(95.0);
    mat.mechanical.computeDerivedProperties();
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    
    return mat;
}

MaterialDefinition Quartzite() {
    MaterialDefinition mat;
    mat.name = "Quartzite";
    mat.category = "Metamorphic";
    mat.subcategory = "Non-foliated";
    mat.description = "Metamorphosed quartz sandstone";
    mat.reference = "Jaeger et al., 2007";
    
    mat.mechanical.density = 2650.0;
    mat.mechanical.youngs_modulus = GPaToPa(80.0);
    mat.mechanical.poisson_ratio = 0.15;
    mat.mechanical.tensile_strength = MPaToPa(15.0);
    mat.mechanical.compressive_strength = MPaToPa(250.0);
    mat.mechanical.cohesion = MPaToPa(40.0);
    mat.mechanical.friction_angle = 55.0;
    mat.mechanical.fracture_toughness_KIc = 1.8e6;
    mat.mechanical.biot_coefficient = 0.25;
    mat.mechanical.grain_bulk_modulus = GPaToPa(37.0);
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.005;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.0001));
    
    mat.thermal.thermal_conductivity = 5.0;
    mat.thermal.specific_heat = 740.0;
    mat.thermal.thermal_expansion = 1.3e-5;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 200.0;
    mat.seismic.Qs = 100.0;
    
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition Hornfels() {
    MaterialDefinition mat;
    mat.name = "Hornfels";
    mat.category = "Metamorphic";
    mat.subcategory = "Non-foliated";
    mat.description = "Contact metamorphic rock, fine-grained";
    mat.reference = "Jaeger et al., 2007";
    
    mat.mechanical.density = 2800.0;
    mat.mechanical.youngs_modulus = GPaToPa(65.0);
    mat.mechanical.poisson_ratio = 0.22;
    mat.mechanical.tensile_strength = MPaToPa(14.0);
    mat.mechanical.compressive_strength = MPaToPa(200.0);
    mat.mechanical.cohesion = MPaToPa(35.0);
    mat.mechanical.friction_angle = 50.0;
    mat.mechanical.fracture_toughness_KIc = 1.6e6;
    mat.mechanical.biot_coefficient = 0.35;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.008;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.001));
    
    mat.thermal.thermal_conductivity = 2.8;
    mat.thermal.specific_heat = 840.0;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 160.0;
    mat.seismic.Qs = 80.0;
    
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition Serpentinite() {
    MaterialDefinition mat;
    mat.name = "Serpentinite";
    mat.category = "Metamorphic";
    mat.subcategory = "Non-foliated";
    mat.description = "Hydrated ultramafic rock, serpentine minerals";
    mat.reference = "Christensen, 1978";
    
    mat.mechanical.density = 2550.0;
    mat.mechanical.youngs_modulus = GPaToPa(40.0);
    mat.mechanical.poisson_ratio = 0.30;
    mat.mechanical.tensile_strength = MPaToPa(5.0);
    mat.mechanical.compressive_strength = MPaToPa(80.0);
    mat.mechanical.cohesion = MPaToPa(10.0);
    mat.mechanical.friction_angle = 32.0;
    mat.mechanical.fracture_toughness_KIc = 0.8e6;
    mat.mechanical.biot_coefficient = 0.60;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.03;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.1));
    
    mat.thermal.thermal_conductivity = 2.5;
    mat.thermal.specific_heat = 1000.0;
    mat.thermal.thermal_expansion = 1.2e-5;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 50.0;
    mat.seismic.Qs = 25.0;
    
    mat.recommended_model = "POROELASTIC";
    
    return mat;
}

MaterialDefinition Soapstone() {
    MaterialDefinition mat;
    mat.name = "Soapstone";
    mat.category = "Metamorphic";
    mat.subcategory = "Non-foliated";
    mat.description = "Talc-rich metamorphic rock, very soft";
    mat.reference = "Jaeger et al., 2007";
    
    mat.mechanical.density = 2750.0;
    mat.mechanical.youngs_modulus = GPaToPa(15.0);
    mat.mechanical.poisson_ratio = 0.28;
    mat.mechanical.tensile_strength = MPaToPa(2.0);
    mat.mechanical.compressive_strength = MPaToPa(30.0);
    mat.mechanical.cohesion = MPaToPa(3.0);
    mat.mechanical.friction_angle = 25.0;
    mat.mechanical.fracture_toughness_KIc = 0.3e6;
    mat.mechanical.biot_coefficient = 0.65;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.01;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.001));
    
    mat.thermal.thermal_conductivity = 6.0;
    mat.thermal.specific_heat = 1020.0;
    mat.thermal.thermal_expansion = 1.0e-5;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 80.0;
    mat.seismic.Qs = 40.0;
    
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition Eclogite() {
    MaterialDefinition mat;
    mat.name = "Eclogite";
    mat.category = "Metamorphic";
    mat.subcategory = "Non-foliated";
    mat.description = "High-P metamorphic rock, garnet + omphacite";
    mat.reference = "Christensen, 1974";
    
    mat.mechanical.density = 3450.0;
    mat.mechanical.youngs_modulus = GPaToPa(140.0);
    mat.mechanical.poisson_ratio = 0.25;
    mat.mechanical.tensile_strength = MPaToPa(25.0);
    mat.mechanical.compressive_strength = MPaToPa(350.0);
    mat.mechanical.cohesion = MPaToPa(60.0);
    mat.mechanical.friction_angle = 50.0;
    mat.mechanical.fracture_toughness_KIc = 3.0e6;
    mat.mechanical.biot_coefficient = 0.20;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.002;
    mat.hydraulic.setIsotropicPermeability(mDToM2(1e-6));
    
    mat.thermal.thermal_conductivity = 3.2;
    mat.thermal.specific_heat = 850.0;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 400.0;
    mat.seismic.Qs = 200.0;
    
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition Granulite() {
    MaterialDefinition mat;
    mat.name = "Granulite";
    mat.category = "Metamorphic";
    mat.subcategory = "Non-foliated";
    mat.description = "High-T metamorphic rock, anhydrous minerals";
    mat.reference = "Jaeger et al., 2007";
    
    mat.mechanical.density = 2800.0;
    mat.mechanical.youngs_modulus = GPaToPa(75.0);
    mat.mechanical.poisson_ratio = 0.24;
    mat.mechanical.tensile_strength = MPaToPa(16.0);
    mat.mechanical.compressive_strength = MPaToPa(220.0);
    mat.mechanical.cohesion = MPaToPa(38.0);
    mat.mechanical.friction_angle = 52.0;
    mat.mechanical.fracture_toughness_KIc = 1.9e6;
    mat.mechanical.biot_coefficient = 0.30;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.003;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.0001));
    
    mat.thermal.thermal_conductivity = 2.8;
    mat.thermal.specific_heat = 830.0;
    mat.thermal.heat_production = 0.5e-6;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 200.0;
    mat.seismic.Qs = 100.0;
    
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition Greenstone() {
    MaterialDefinition mat;
    mat.name = "Greenstone";
    mat.category = "Metamorphic";
    mat.subcategory = "Non-foliated";
    mat.description = "Low-grade metamorphosed mafic volcanic rock";
    mat.reference = "Jaeger et al., 2007";
    
    mat.mechanical.density = 2900.0;
    mat.mechanical.youngs_modulus = GPaToPa(60.0);
    mat.mechanical.poisson_ratio = 0.26;
    mat.mechanical.tensile_strength = MPaToPa(12.0);
    mat.mechanical.compressive_strength = MPaToPa(170.0);
    mat.mechanical.cohesion = MPaToPa(30.0);
    mat.mechanical.friction_angle = 46.0;
    mat.mechanical.fracture_toughness_KIc = 1.5e6;
    mat.mechanical.biot_coefficient = 0.45;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.02;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.01));
    
    mat.thermal.thermal_conductivity = 2.4;
    mat.thermal.specific_heat = 890.0;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 130.0;
    mat.seismic.Qs = 65.0;
    
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition Skarn() {
    MaterialDefinition mat;
    mat.name = "Skarn";
    mat.category = "Metamorphic";
    mat.subcategory = "Non-foliated";
    mat.description = "Contact-metasomatic calc-silicate rock";
    mat.reference = "Meinert et al., 2005";
    
    mat.mechanical.density = 3200.0;
    mat.mechanical.youngs_modulus = GPaToPa(70.0);
    mat.mechanical.poisson_ratio = 0.26;
    mat.mechanical.tensile_strength = MPaToPa(13.0);
    mat.mechanical.compressive_strength = MPaToPa(190.0);
    mat.mechanical.cohesion = MPaToPa(33.0);
    mat.mechanical.friction_angle = 48.0;
    mat.mechanical.fracture_toughness_KIc = 1.7e6;
    mat.mechanical.biot_coefficient = 0.40;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.03;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.1));
    
    mat.thermal.thermal_conductivity = 3.0;
    mat.thermal.specific_heat = 820.0;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 150.0;
    mat.seismic.Qs = 75.0;
    
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition Tactite() {
    MaterialDefinition mat = Skarn();
    mat.name = "Tactite";
    mat.description = "Lime-bearing hornfels, contact metamorphism";
    
    mat.mechanical.density = 3000.0;
    mat.mechanical.youngs_modulus = GPaToPa(65.0);
    mat.mechanical.compressive_strength = MPaToPa(175.0);
    mat.mechanical.computeDerivedProperties();
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    
    return mat;
}

} // namespace Metamorphic

// =============================================================================
// Unconsolidated Materials - namespace Materials::Unconsolidated
// =============================================================================

namespace Unconsolidated {

MaterialDefinition Sand() {
    MaterialDefinition mat;
    mat.name = "Sand";
    mat.category = "Unconsolidated";
    mat.subcategory = "Granular";
    mat.description = "Loose sand, medium grain size";
    mat.reference = "Das, 2010; Jaeger et al., 2007";
    
    mat.mechanical.density = 1700.0;
    mat.mechanical.youngs_modulus = GPaToPa(0.05);
    mat.mechanical.poisson_ratio = 0.30;
    mat.mechanical.tensile_strength = MPaToPa(0.0);
    mat.mechanical.compressive_strength = MPaToPa(0.5);
    mat.mechanical.cohesion = MPaToPa(0.0);
    mat.mechanical.friction_angle = 33.0;
    mat.mechanical.biot_coefficient = 1.0;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.35;
    mat.hydraulic.setIsotropicPermeability(mDToM2(10000.0));
    mat.hydraulic.compressibility = 1.0e-8;
    
    mat.thermal.thermal_conductivity = 0.3;
    mat.thermal.specific_heat = 800.0;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 10.0;
    mat.seismic.Qs = 5.0;
    
    mat.is_porous = true;
    mat.recommended_model = "ELASTOPLASTIC_MC";
    
    return mat;
}

MaterialDefinition SandFine() {
    MaterialDefinition mat = Sand();
    mat.name = "Sand_Fine";
    mat.description = "Fine sand, 0.1-0.25mm";
    
    mat.mechanical.friction_angle = 30.0;
    
    mat.hydraulic.porosity = 0.40;
    mat.hydraulic.setIsotropicPermeability(mDToM2(1000.0));
    
    return mat;
}

MaterialDefinition SandMedium() {
    MaterialDefinition mat = Sand();
    mat.name = "Sand_Medium";
    mat.description = "Medium sand, 0.25-0.5mm";
    
    mat.mechanical.friction_angle = 33.0;
    
    mat.hydraulic.porosity = 0.35;
    mat.hydraulic.setIsotropicPermeability(mDToM2(10000.0));
    
    return mat;
}

MaterialDefinition SandCoarse() {
    MaterialDefinition mat = Sand();
    mat.name = "Sand_Coarse";
    mat.description = "Coarse sand, 0.5-2mm";
    
    mat.mechanical.friction_angle = 36.0;
    
    mat.hydraulic.porosity = 0.32;
    mat.hydraulic.setIsotropicPermeability(mDToM2(50000.0));
    
    return mat;
}

MaterialDefinition SandSilty() {
    MaterialDefinition mat = Sand();
    mat.name = "Sand_Silty";
    mat.description = "Silty sand with fines";
    
    mat.mechanical.density = 1800.0;
    mat.mechanical.friction_angle = 28.0;
    mat.mechanical.cohesion = MPaToPa(0.005);
    
    mat.hydraulic.porosity = 0.38;
    mat.hydraulic.setIsotropicPermeability(mDToM2(100.0));
    
    return mat;
}

MaterialDefinition SandClayey() {
    MaterialDefinition mat = Sand();
    mat.name = "Sand_Clayey";
    mat.description = "Clayey sand with clay matrix";
    
    mat.mechanical.density = 1900.0;
    mat.mechanical.friction_angle = 25.0;
    mat.mechanical.cohesion = MPaToPa(0.010);
    
    mat.hydraulic.porosity = 0.35;
    mat.hydraulic.setIsotropicPermeability(mDToM2(10.0));
    
    return mat;
}

MaterialDefinition Gravel() {
    MaterialDefinition mat;
    mat.name = "Gravel";
    mat.category = "Unconsolidated";
    mat.subcategory = "Granular";
    mat.description = "Loose gravel, 2-64mm particles";
    mat.reference = "Das, 2010";
    
    mat.mechanical.density = 1900.0;
    mat.mechanical.youngs_modulus = GPaToPa(0.1);
    mat.mechanical.poisson_ratio = 0.28;
    mat.mechanical.tensile_strength = MPaToPa(0.0);
    mat.mechanical.compressive_strength = MPaToPa(1.0);
    mat.mechanical.cohesion = MPaToPa(0.0);
    mat.mechanical.friction_angle = 40.0;
    mat.mechanical.biot_coefficient = 1.0;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.30;
    mat.hydraulic.setIsotropicPermeability(mDToM2(100000.0));
    
    mat.thermal.thermal_conductivity = 0.5;
    mat.thermal.specific_heat = 780.0;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 15.0;
    mat.seismic.Qs = 8.0;
    
    mat.is_porous = true;
    mat.recommended_model = "ELASTOPLASTIC_MC";
    
    return mat;
}

MaterialDefinition GravelSandy() {
    MaterialDefinition mat = Gravel();
    mat.name = "Gravel_Sandy";
    mat.description = "Sandy gravel with sand matrix";
    
    mat.mechanical.friction_angle = 38.0;
    
    mat.hydraulic.porosity = 0.28;
    mat.hydraulic.setIsotropicPermeability(mDToM2(50000.0));
    
    return mat;
}

MaterialDefinition Clay() {
    MaterialDefinition mat;
    mat.name = "Clay";
    mat.category = "Unconsolidated";
    mat.subcategory = "Cohesive";
    mat.description = "Soft to stiff clay";
    mat.reference = "Das, 2010; Skempton, 1985";
    
    mat.mechanical.density = 1800.0;
    mat.mechanical.youngs_modulus = GPaToPa(0.02);
    mat.mechanical.poisson_ratio = 0.40;
    mat.mechanical.tensile_strength = MPaToPa(0.01);
    mat.mechanical.compressive_strength = MPaToPa(0.1);
    mat.mechanical.cohesion = MPaToPa(0.025);
    mat.mechanical.friction_angle = 20.0;
    mat.mechanical.biot_coefficient = 1.0;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.50;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.001));
    mat.hydraulic.compressibility = 1.0e-7;
    
    mat.thermal.thermal_conductivity = 1.0;
    mat.thermal.specific_heat = 900.0;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 8.0;
    mat.seismic.Qs = 4.0;
    
    mat.is_porous = true;
    mat.recommended_model = "ELASTOPLASTIC_MC";
    
    return mat;
}

MaterialDefinition ClayKaolinite() {
    MaterialDefinition mat = Clay();
    mat.name = "Clay_Kaolinite";
    mat.description = "Kaolinite-dominant clay";
    
    mat.mechanical.friction_angle = 25.0;
    mat.mechanical.cohesion = MPaToPa(0.020);
    
    mat.hydraulic.porosity = 0.45;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.01));
    
    return mat;
}

MaterialDefinition ClayIllite() {
    MaterialDefinition mat = Clay();
    mat.name = "Clay_Illite";
    mat.description = "Illite-dominant clay";
    
    mat.mechanical.friction_angle = 22.0;
    mat.mechanical.cohesion = MPaToPa(0.022);
    
    mat.hydraulic.porosity = 0.48;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.005));
    
    return mat;
}

MaterialDefinition ClayMontmorillonite() {
    MaterialDefinition mat = Clay();
    mat.name = "Clay_Montmorillonite";
    mat.description = "Montmorillonite-dominant (swelling) clay";
    
    mat.mechanical.friction_angle = 12.0;
    mat.mechanical.cohesion = MPaToPa(0.015);
    mat.mechanical.youngs_modulus = GPaToPa(0.01);
    
    mat.hydraulic.porosity = 0.60;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.0001));
    
    return mat;
}

MaterialDefinition ClayBentonite() {
    MaterialDefinition mat = Clay();
    mat.name = "Clay_Bentonite";
    mat.description = "Bentonite clay (high smectite)";
    
    mat.mechanical.density = 2000.0;
    mat.mechanical.friction_angle = 10.0;
    mat.mechanical.cohesion = MPaToPa(0.012);
    mat.mechanical.youngs_modulus = GPaToPa(0.008);
    
    mat.hydraulic.porosity = 0.65;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.00001));
    
    return mat;
}

MaterialDefinition Silt() {
    MaterialDefinition mat;
    mat.name = "Silt";
    mat.category = "Unconsolidated";
    mat.subcategory = "Cohesive";
    mat.description = "Silt, particle size 0.002-0.063mm";
    mat.reference = "Das, 2010";
    
    mat.mechanical.density = 1700.0;
    mat.mechanical.youngs_modulus = GPaToPa(0.03);
    mat.mechanical.poisson_ratio = 0.35;
    mat.mechanical.tensile_strength = MPaToPa(0.005);
    mat.mechanical.compressive_strength = MPaToPa(0.08);
    mat.mechanical.cohesion = MPaToPa(0.010);
    mat.mechanical.friction_angle = 28.0;
    mat.mechanical.biot_coefficient = 1.0;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.42;
    mat.hydraulic.setIsotropicPermeability(mDToM2(1.0));
    
    mat.thermal.thermal_conductivity = 0.8;
    mat.thermal.specific_heat = 850.0;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 10.0;
    mat.seismic.Qs = 5.0;
    
    mat.is_porous = true;
    mat.recommended_model = "ELASTOPLASTIC_MC";
    
    return mat;
}

MaterialDefinition Loam() {
    MaterialDefinition mat;
    mat.name = "Loam";
    mat.category = "Unconsolidated";
    mat.subcategory = "Mixed";
    mat.description = "Loam soil, mixed sand/silt/clay";
    mat.reference = "Das, 2010";
    
    mat.mechanical.density = 1600.0;
    mat.mechanical.youngs_modulus = GPaToPa(0.02);
    mat.mechanical.poisson_ratio = 0.35;
    mat.mechanical.tensile_strength = MPaToPa(0.003);
    mat.mechanical.compressive_strength = MPaToPa(0.05);
    mat.mechanical.cohesion = MPaToPa(0.008);
    mat.mechanical.friction_angle = 30.0;
    mat.mechanical.biot_coefficient = 1.0;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.45;
    mat.hydraulic.setIsotropicPermeability(mDToM2(50.0));
    
    mat.thermal.thermal_conductivity = 0.6;
    mat.thermal.specific_heat = 900.0;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 8.0;
    mat.seismic.Qs = 4.0;
    
    mat.is_porous = true;
    mat.recommended_model = "ELASTOPLASTIC_MC";
    
    return mat;
}

MaterialDefinition Till() {
    MaterialDefinition mat;
    mat.name = "Till";
    mat.category = "Unconsolidated";
    mat.subcategory = "Mixed";
    mat.description = "Glacial till, unsorted";
    mat.reference = "Boulton & Paul, 1976";
    
    mat.mechanical.density = 2100.0;
    mat.mechanical.youngs_modulus = GPaToPa(0.1);
    mat.mechanical.poisson_ratio = 0.32;
    mat.mechanical.tensile_strength = MPaToPa(0.01);
    mat.mechanical.compressive_strength = MPaToPa(0.5);
    mat.mechanical.cohesion = MPaToPa(0.015);
    mat.mechanical.friction_angle = 32.0;
    mat.mechanical.biot_coefficient = 0.95;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.28;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.1));
    
    mat.thermal.thermal_conductivity = 1.5;
    mat.thermal.specific_heat = 850.0;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 15.0;
    mat.seismic.Qs = 8.0;
    
    mat.is_porous = true;
    mat.recommended_model = "ELASTOPLASTIC_MC";
    
    return mat;
}

MaterialDefinition Alluvium() {
    MaterialDefinition mat;
    mat.name = "Alluvium";
    mat.category = "Unconsolidated";
    mat.subcategory = "Fluvial";
    mat.description = "River-deposited sediments";
    mat.reference = "Das, 2010";
    
    mat.mechanical.density = 1800.0;
    mat.mechanical.youngs_modulus = GPaToPa(0.04);
    mat.mechanical.poisson_ratio = 0.33;
    mat.mechanical.tensile_strength = MPaToPa(0.0);
    mat.mechanical.compressive_strength = MPaToPa(0.2);
    mat.mechanical.cohesion = MPaToPa(0.005);
    mat.mechanical.friction_angle = 32.0;
    mat.mechanical.biot_coefficient = 1.0;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.35;
    mat.hydraulic.setIsotropicPermeability(mDToM2(1000.0));
    
    mat.thermal.thermal_conductivity = 0.8;
    mat.thermal.specific_heat = 850.0;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 12.0;
    mat.seismic.Qs = 6.0;
    
    mat.is_porous = true;
    mat.recommended_model = "ELASTOPLASTIC_MC";
    
    return mat;
}

MaterialDefinition ColluviumTalus() {
    MaterialDefinition mat;
    mat.name = "Colluvium_Talus";
    mat.category = "Unconsolidated";
    mat.subcategory = "Mass-wasting";
    mat.description = "Slope debris and talus";
    mat.reference = "Das, 2010";
    
    mat.mechanical.density = 1900.0;
    mat.mechanical.youngs_modulus = GPaToPa(0.08);
    mat.mechanical.poisson_ratio = 0.30;
    mat.mechanical.tensile_strength = MPaToPa(0.0);
    mat.mechanical.compressive_strength = MPaToPa(0.3);
    mat.mechanical.cohesion = MPaToPa(0.0);
    mat.mechanical.friction_angle = 35.0;
    mat.mechanical.biot_coefficient = 1.0;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.38;
    mat.hydraulic.setIsotropicPermeability(mDToM2(10000.0));
    
    mat.thermal.thermal_conductivity = 0.5;
    mat.thermal.specific_heat = 800.0;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 10.0;
    mat.seismic.Qs = 5.0;
    
    mat.is_porous = true;
    mat.recommended_model = "ELASTOPLASTIC_MC";
    
    return mat;
}

MaterialDefinition Laterite() {
    MaterialDefinition mat;
    mat.name = "Laterite";
    mat.category = "Unconsolidated";
    mat.subcategory = "Residual";
    mat.description = "Iron/aluminum-rich tropical residual soil";
    mat.reference = "Gidigasu, 1976";
    
    mat.mechanical.density = 2000.0;
    mat.mechanical.youngs_modulus = GPaToPa(0.08);
    mat.mechanical.poisson_ratio = 0.32;
    mat.mechanical.tensile_strength = MPaToPa(0.02);
    mat.mechanical.compressive_strength = MPaToPa(0.4);
    mat.mechanical.cohesion = MPaToPa(0.030);
    mat.mechanical.friction_angle = 28.0;
    mat.mechanical.biot_coefficient = 0.90;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.40;
    mat.hydraulic.setIsotropicPermeability(mDToM2(10.0));
    
    mat.thermal.thermal_conductivity = 1.2;
    mat.thermal.specific_heat = 850.0;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 15.0;
    mat.seismic.Qs = 8.0;
    
    mat.is_porous = true;
    mat.recommended_model = "ELASTOPLASTIC_MC";
    
    return mat;
}

MaterialDefinition Saprolite() {
    MaterialDefinition mat;
    mat.name = "Saprolite";
    mat.category = "Unconsolidated";
    mat.subcategory = "Residual";
    mat.description = "Weathered bedrock retaining rock structure";
    mat.reference = "Deere & Patton, 1971";
    
    mat.mechanical.density = 2100.0;
    mat.mechanical.youngs_modulus = GPaToPa(0.15);
    mat.mechanical.poisson_ratio = 0.30;
    mat.mechanical.tensile_strength = MPaToPa(0.05);
    mat.mechanical.compressive_strength = MPaToPa(1.0);
    mat.mechanical.cohesion = MPaToPa(0.040);
    mat.mechanical.friction_angle = 32.0;
    mat.mechanical.biot_coefficient = 0.85;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.35;
    mat.hydraulic.setIsotropicPermeability(mDToM2(5.0));
    
    mat.thermal.thermal_conductivity = 1.5;
    mat.thermal.specific_heat = 850.0;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 20.0;
    mat.seismic.Qs = 10.0;
    
    mat.is_porous = true;
    mat.recommended_model = "ELASTOPLASTIC_MC";
    
    return mat;
}

MaterialDefinition Regolith() {
    MaterialDefinition mat;
    mat.name = "Regolith";
    mat.category = "Unconsolidated";
    mat.subcategory = "Residual";
    mat.description = "General weathered surface layer";
    mat.reference = "Scott & Pain, 2008";
    
    mat.mechanical.density = 1800.0;
    mat.mechanical.youngs_modulus = GPaToPa(0.05);
    mat.mechanical.poisson_ratio = 0.32;
    mat.mechanical.tensile_strength = MPaToPa(0.01);
    mat.mechanical.compressive_strength = MPaToPa(0.3);
    mat.mechanical.cohesion = MPaToPa(0.010);
    mat.mechanical.friction_angle = 30.0;
    mat.mechanical.biot_coefficient = 0.95;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.40;
    mat.hydraulic.setIsotropicPermeability(mDToM2(100.0));
    
    mat.thermal.thermal_conductivity = 0.8;
    mat.thermal.specific_heat = 850.0;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 12.0;
    mat.seismic.Qs = 6.0;
    
    mat.is_porous = true;
    mat.recommended_model = "ELASTOPLASTIC_MC";
    
    return mat;
}

} // namespace Unconsolidated

// =============================================================================
// Minerals - namespace Materials::Minerals
// =============================================================================

namespace Minerals {

MaterialDefinition Quartz() {
    MaterialDefinition mat;
    mat.name = "Quartz";
    mat.category = "Minerals";
    mat.subcategory = "Silicates";
    mat.description = "Crystalline silica, SiO2";
    mat.reference = "Bass, 1995; Mavko et al., 2009";
    
    mat.mechanical.density = 2650.0;
    mat.mechanical.youngs_modulus = GPaToPa(95.0);
    mat.mechanical.poisson_ratio = 0.08;
    mat.mechanical.tensile_strength = MPaToPa(50.0);
    mat.mechanical.compressive_strength = MPaToPa(1100.0);
    mat.mechanical.cohesion = MPaToPa(100.0);
    mat.mechanical.friction_angle = 60.0;
    mat.mechanical.fracture_toughness_KIc = 1.0e6;
    mat.mechanical.biot_coefficient = 0.0;
    mat.mechanical.grain_bulk_modulus = GPaToPa(37.0);
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.0;
    mat.hydraulic.setIsotropicPermeability(0.0);
    
    mat.thermal.thermal_conductivity = 7.7;
    mat.thermal.specific_heat = 740.0;
    mat.thermal.thermal_expansion = 1.2e-5;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 1000.0;
    mat.seismic.Qs = 500.0;
    
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition Feldspar() {
    MaterialDefinition mat;
    mat.name = "Feldspar";
    mat.category = "Minerals";
    mat.subcategory = "Silicates";
    mat.description = "Alkali/plagioclase feldspar group";
    mat.reference = "Bass, 1995";
    
    mat.mechanical.density = 2600.0;
    mat.mechanical.youngs_modulus = GPaToPa(70.0);
    mat.mechanical.poisson_ratio = 0.28;
    mat.mechanical.tensile_strength = MPaToPa(30.0);
    mat.mechanical.compressive_strength = MPaToPa(400.0);
    mat.mechanical.cohesion = MPaToPa(60.0);
    mat.mechanical.friction_angle = 55.0;
    mat.mechanical.grain_bulk_modulus = GPaToPa(55.0);
    mat.mechanical.computeDerivedProperties();
    
    mat.thermal.thermal_conductivity = 2.3;
    mat.thermal.specific_heat = 730.0;
    mat.thermal.thermal_expansion = 5.0e-6;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 800.0;
    mat.seismic.Qs = 400.0;
    
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition FeldsparOrthoclase() {
    MaterialDefinition mat = Feldspar();
    mat.name = "Feldspar_Orthoclase";
    mat.description = "Potassium feldspar, KAlSi3O8";
    mat.mechanical.density = 2560.0;
    mat.mechanical.youngs_modulus = GPaToPa(67.0);
    mat.mechanical.grain_bulk_modulus = GPaToPa(48.0);
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    return mat;
}

MaterialDefinition FeldsparPlagioclase() {
    MaterialDefinition mat = Feldspar();
    mat.name = "Feldspar_Plagioclase";
    mat.description = "Plagioclase feldspar series";
    mat.mechanical.density = 2650.0;
    mat.mechanical.youngs_modulus = GPaToPa(75.0);
    mat.mechanical.grain_bulk_modulus = GPaToPa(62.0);
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    return mat;
}

MaterialDefinition Mica() {
    MaterialDefinition mat;
    mat.name = "Mica";
    mat.category = "Minerals";
    mat.subcategory = "Silicates";
    mat.description = "Phyllosilicate mica group";
    mat.reference = "Bass, 1995";
    
    mat.mechanical.density = 2850.0;
    mat.mechanical.youngs_modulus = GPaToPa(60.0);
    mat.mechanical.poisson_ratio = 0.25;
    mat.mechanical.biot_coefficient = 0.0;
    mat.mechanical.computeDerivedProperties();
    
    mat.thermal.thermal_conductivity = 2.0;
    mat.thermal.specific_heat = 820.0;
    
    mat.anisotropic.is_anisotropic = true;
    mat.anisotropic.symmetry_type = "VTI";
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.is_anisotropic_flag = true;
    mat.recommended_model = "ANISOTROPIC_VTI";
    
    return mat;
}

MaterialDefinition MicaMuscovite() {
    MaterialDefinition mat = Mica();
    mat.name = "Mica_Muscovite";
    mat.description = "White mica, KAl2(AlSi3O10)(OH)2";
    mat.mechanical.density = 2830.0;
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    return mat;
}

MaterialDefinition MicaBiotite() {
    MaterialDefinition mat = Mica();
    mat.name = "Mica_Biotite";
    mat.description = "Black mica, K(Mg,Fe)3(AlSi3O10)(OH)2";
    mat.mechanical.density = 3000.0;
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    return mat;
}

MaterialDefinition Olivine() {
    MaterialDefinition mat;
    mat.name = "Olivine";
    mat.category = "Minerals";
    mat.subcategory = "Silicates";
    mat.description = "(Mg,Fe)2SiO4 - forsterite to fayalite";
    mat.reference = "Bass, 1995; Christensen, 1974";
    
    mat.mechanical.density = 3350.0;
    mat.mechanical.youngs_modulus = GPaToPa(195.0);
    mat.mechanical.poisson_ratio = 0.24;
    mat.mechanical.grain_bulk_modulus = GPaToPa(128.0);
    mat.mechanical.computeDerivedProperties();
    
    mat.thermal.thermal_conductivity = 5.0;
    mat.thermal.specific_heat = 820.0;
    mat.thermal.thermal_expansion = 2.8e-5;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 1200.0;
    mat.seismic.Qs = 600.0;
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition Pyroxene() {
    MaterialDefinition mat;
    mat.name = "Pyroxene";
    mat.category = "Minerals";
    mat.subcategory = "Silicates";
    mat.description = "Single-chain inosilicate";
    mat.reference = "Bass, 1995";
    
    mat.mechanical.density = 3300.0;
    mat.mechanical.youngs_modulus = GPaToPa(160.0);
    mat.mechanical.poisson_ratio = 0.25;
    mat.mechanical.grain_bulk_modulus = GPaToPa(108.0);
    mat.mechanical.computeDerivedProperties();
    
    mat.thermal.thermal_conductivity = 4.5;
    mat.thermal.specific_heat = 800.0;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition Amphibole() {
    MaterialDefinition mat;
    mat.name = "Amphibole";
    mat.category = "Minerals";
    mat.subcategory = "Silicates";
    mat.description = "Double-chain inosilicate";
    mat.reference = "Bass, 1995";
    
    mat.mechanical.density = 3200.0;
    mat.mechanical.youngs_modulus = GPaToPa(120.0);
    mat.mechanical.poisson_ratio = 0.28;
    mat.mechanical.grain_bulk_modulus = GPaToPa(85.0);
    mat.mechanical.computeDerivedProperties();
    
    mat.thermal.thermal_conductivity = 3.5;
    mat.thermal.specific_heat = 850.0;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition Hornblende() {
    MaterialDefinition mat = Amphibole();
    mat.name = "Hornblende";
    mat.description = "Common amphibole mineral";
    mat.mechanical.density = 3100.0;
    mat.mechanical.youngs_modulus = GPaToPa(115.0);
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    return mat;
}

MaterialDefinition Garnet() {
    MaterialDefinition mat;
    mat.name = "Garnet";
    mat.category = "Minerals";
    mat.subcategory = "Silicates";
    mat.description = "Nesosilicate garnet group";
    mat.reference = "Bass, 1995";
    
    mat.mechanical.density = 3800.0;
    mat.mechanical.youngs_modulus = GPaToPa(240.0);
    mat.mechanical.poisson_ratio = 0.27;
    mat.mechanical.grain_bulk_modulus = GPaToPa(170.0);
    mat.mechanical.computeDerivedProperties();
    
    mat.thermal.thermal_conductivity = 3.2;
    mat.thermal.specific_heat = 740.0;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition Kyanite() {
    MaterialDefinition mat;
    mat.name = "Kyanite";
    mat.category = "Minerals";
    mat.subcategory = "Silicates";
    mat.description = "Al2SiO5 polymorph";
    mat.mechanical.density = 3600.0;
    mat.mechanical.youngs_modulus = GPaToPa(200.0);
    mat.mechanical.poisson_ratio = 0.25;
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition Sillimanite() {
    MaterialDefinition mat = Kyanite();
    mat.name = "Sillimanite";
    mat.description = "Al2SiO5 high-T polymorph";
    mat.mechanical.density = 3240.0;
    mat.mechanical.youngs_modulus = GPaToPa(170.0);
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    return mat;
}

MaterialDefinition Andalusite() {
    MaterialDefinition mat = Kyanite();
    mat.name = "Andalusite";
    mat.description = "Al2SiO5 low-P polymorph";
    mat.mechanical.density = 3150.0;
    mat.mechanical.youngs_modulus = GPaToPa(160.0);
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    return mat;
}

MaterialDefinition Tourmaline() {
    MaterialDefinition mat;
    mat.name = "Tourmaline";
    mat.category = "Minerals";
    mat.subcategory = "Silicates";
    mat.description = "Boron silicate mineral group";
    mat.mechanical.density = 3100.0;
    mat.mechanical.youngs_modulus = GPaToPa(150.0);
    mat.mechanical.poisson_ratio = 0.26;
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition Zircon() {
    MaterialDefinition mat;
    mat.name = "Zircon";
    mat.category = "Minerals";
    mat.subcategory = "Silicates";
    mat.description = "ZrSiO4";
    mat.mechanical.density = 4650.0;
    mat.mechanical.youngs_modulus = GPaToPa(280.0);
    mat.mechanical.poisson_ratio = 0.24;
    mat.mechanical.computeDerivedProperties();
    mat.thermal.heat_production = 5.0e-6;
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition Epidote() {
    MaterialDefinition mat;
    mat.name = "Epidote";
    mat.category = "Minerals";
    mat.subcategory = "Silicates";
    mat.description = "Ca2(Al,Fe)3(SiO4)(Si2O7)O(OH)";
    mat.mechanical.density = 3400.0;
    mat.mechanical.youngs_modulus = GPaToPa(160.0);
    mat.mechanical.poisson_ratio = 0.28;
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition Chlorite() {
    MaterialDefinition mat;
    mat.name = "Chlorite";
    mat.category = "Minerals";
    mat.subcategory = "Silicates";
    mat.description = "Phyllosilicate chlorite group";
    mat.mechanical.density = 2800.0;
    mat.mechanical.youngs_modulus = GPaToPa(50.0);
    mat.mechanical.poisson_ratio = 0.30;
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition Serpentine() {
    MaterialDefinition mat;
    mat.name = "Serpentine";
    mat.category = "Minerals";
    mat.subcategory = "Silicates";
    mat.description = "Mg3Si2O5(OH)4 serpentine group";
    mat.mechanical.density = 2550.0;
    mat.mechanical.youngs_modulus = GPaToPa(50.0);
    mat.mechanical.poisson_ratio = 0.32;
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition Talc() {
    MaterialDefinition mat;
    mat.name = "Talc";
    mat.category = "Minerals";
    mat.subcategory = "Silicates";
    mat.description = "Mg3Si4O10(OH)2 - softest mineral";
    mat.mechanical.density = 2750.0;
    mat.mechanical.youngs_modulus = GPaToPa(20.0);
    mat.mechanical.poisson_ratio = 0.30;
    mat.mechanical.tensile_strength = MPaToPa(1.0);
    mat.mechanical.compressive_strength = MPaToPa(20.0);
    mat.mechanical.computeDerivedProperties();
    mat.thermal.thermal_conductivity = 6.0;
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition Zeolites() {
    MaterialDefinition mat;
    mat.name = "Zeolites";
    mat.category = "Minerals";
    mat.subcategory = "Silicates";
    mat.description = "Hydrated aluminosilicate framework";
    mat.mechanical.density = 2200.0;
    mat.mechanical.youngs_modulus = GPaToPa(35.0);
    mat.mechanical.poisson_ratio = 0.25;
    mat.mechanical.computeDerivedProperties();
    mat.hydraulic.porosity = 0.30;
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.is_porous = true;
    mat.recommended_model = "POROELASTIC";
    return mat;
}

// Carbonates
MaterialDefinition Calcite() {
    MaterialDefinition mat;
    mat.name = "Calcite";
    mat.category = "Minerals";
    mat.subcategory = "Carbonates";
    mat.description = "CaCO3 - rhombohedral";
    mat.reference = "Bass, 1995";
    mat.mechanical.density = 2710.0;
    mat.mechanical.youngs_modulus = GPaToPa(84.0);
    mat.mechanical.poisson_ratio = 0.32;
    mat.mechanical.grain_bulk_modulus = GPaToPa(77.0);
    mat.mechanical.computeDerivedProperties();
    mat.thermal.thermal_conductivity = 3.3;
    mat.thermal.thermal_expansion = 2.5e-5;
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition Aragonite() {
    MaterialDefinition mat = Calcite();
    mat.name = "Aragonite";
    mat.description = "CaCO3 - orthorhombic polymorph";
    mat.mechanical.density = 2930.0;
    mat.mechanical.youngs_modulus = GPaToPa(90.0);
    mat.mechanical.grain_bulk_modulus = GPaToPa(65.0);
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    return mat;
}

MaterialDefinition DolomiteMinite() {
    MaterialDefinition mat;
    mat.name = "Dolomite_Mineral";
    mat.category = "Minerals";
    mat.subcategory = "Carbonates";
    mat.description = "CaMg(CO3)2";
    mat.mechanical.density = 2870.0;
    mat.mechanical.youngs_modulus = GPaToPa(117.0);
    mat.mechanical.poisson_ratio = 0.30;
    mat.mechanical.grain_bulk_modulus = GPaToPa(95.0);
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition Siderite() {
    MaterialDefinition mat;
    mat.name = "Siderite";
    mat.category = "Minerals";
    mat.subcategory = "Carbonates";
    mat.description = "FeCO3";
    mat.mechanical.density = 3960.0;
    mat.mechanical.youngs_modulus = GPaToPa(130.0);
    mat.mechanical.poisson_ratio = 0.28;
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition Magneite() {
    MaterialDefinition mat;
    mat.name = "Magnesite";
    mat.category = "Minerals";
    mat.subcategory = "Carbonates";
    mat.description = "MgCO3";
    mat.mechanical.density = 3000.0;
    mat.mechanical.youngs_modulus = GPaToPa(120.0);
    mat.mechanical.poisson_ratio = 0.28;
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition Rhodochrosite() {
    MaterialDefinition mat;
    mat.name = "Rhodochrosite";
    mat.category = "Minerals";
    mat.subcategory = "Carbonates";
    mat.description = "MnCO3";
    mat.mechanical.density = 3700.0;
    mat.mechanical.youngs_modulus = GPaToPa(100.0);
    mat.mechanical.poisson_ratio = 0.30;
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

// Sulfates
MaterialDefinition GypsumMineral() {
    MaterialDefinition mat;
    mat.name = "Gypsum_Mineral";
    mat.category = "Minerals";
    mat.subcategory = "Sulfates";
    mat.description = "CaSO4·2H2O";
    mat.mechanical.density = 2320.0;
    mat.mechanical.youngs_modulus = GPaToPa(42.0);
    mat.mechanical.poisson_ratio = 0.28;
    mat.mechanical.grain_bulk_modulus = GPaToPa(42.0);
    mat.mechanical.computeDerivedProperties();
    mat.thermal.thermal_conductivity = 1.3;
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition AnhydriteMineral() {
    MaterialDefinition mat;
    mat.name = "Anhydrite_Mineral";
    mat.category = "Minerals";
    mat.subcategory = "Sulfates";
    mat.description = "CaSO4";
    mat.mechanical.density = 2980.0;
    mat.mechanical.youngs_modulus = GPaToPa(92.0);
    mat.mechanical.poisson_ratio = 0.30;
    mat.mechanical.grain_bulk_modulus = GPaToPa(55.0);
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition Barite() {
    MaterialDefinition mat;
    mat.name = "Barite";
    mat.category = "Minerals";
    mat.subcategory = "Sulfates";
    mat.description = "BaSO4";
    mat.mechanical.density = 4500.0;
    mat.mechanical.youngs_modulus = GPaToPa(85.0);
    mat.mechanical.poisson_ratio = 0.29;
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition Celestite() {
    MaterialDefinition mat;
    mat.name = "Celestite";
    mat.category = "Minerals";
    mat.subcategory = "Sulfates";
    mat.description = "SrSO4";
    mat.mechanical.density = 3970.0;
    mat.mechanical.youngs_modulus = GPaToPa(78.0);
    mat.mechanical.poisson_ratio = 0.28;
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

// Halides
MaterialDefinition HaliteMineral() {
    MaterialDefinition mat;
    mat.name = "Halite_Mineral";
    mat.category = "Minerals";
    mat.subcategory = "Halides";
    mat.description = "NaCl";
    mat.mechanical.density = 2160.0;
    mat.mechanical.youngs_modulus = GPaToPa(37.0);
    mat.mechanical.poisson_ratio = 0.25;
    mat.mechanical.grain_bulk_modulus = GPaToPa(24.0);
    mat.mechanical.viscosity = 1.0e18;
    mat.mechanical.computeDerivedProperties();
    mat.thermal.thermal_conductivity = 5.5;
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.is_viscoelastic = true;
    mat.recommended_model = "VISCOELASTIC_SLS";
    return mat;
}

MaterialDefinition Fluorite() {
    MaterialDefinition mat;
    mat.name = "Fluorite";
    mat.category = "Minerals";
    mat.subcategory = "Halides";
    mat.description = "CaF2";
    mat.mechanical.density = 3180.0;
    mat.mechanical.youngs_modulus = GPaToPa(110.0);
    mat.mechanical.poisson_ratio = 0.26;
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition SylviteMineral() {
    MaterialDefinition mat;
    mat.name = "Sylvite_Mineral";
    mat.category = "Minerals";
    mat.subcategory = "Halides";
    mat.description = "KCl";
    mat.mechanical.density = 1990.0;
    mat.mechanical.youngs_modulus = GPaToPa(28.0);
    mat.mechanical.poisson_ratio = 0.26;
    mat.mechanical.grain_bulk_modulus = GPaToPa(18.0);
    mat.mechanical.viscosity = 5.0e17;
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.is_viscoelastic = true;
    mat.recommended_model = "VISCOELASTIC_SLS";
    return mat;
}

// Oxides
MaterialDefinition Magnetite() {
    MaterialDefinition mat;
    mat.name = "Magnetite";
    mat.category = "Minerals";
    mat.subcategory = "Oxides";
    mat.description = "Fe3O4";
    mat.mechanical.density = 5200.0;
    mat.mechanical.youngs_modulus = GPaToPa(230.0);
    mat.mechanical.poisson_ratio = 0.28;
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition Hematite() {
    MaterialDefinition mat;
    mat.name = "Hematite";
    mat.category = "Minerals";
    mat.subcategory = "Oxides";
    mat.description = "Fe2O3";
    mat.mechanical.density = 5260.0;
    mat.mechanical.youngs_modulus = GPaToPa(240.0);
    mat.mechanical.poisson_ratio = 0.27;
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition Ilmenite() {
    MaterialDefinition mat;
    mat.name = "Ilmenite";
    mat.category = "Minerals";
    mat.subcategory = "Oxides";
    mat.description = "FeTiO3";
    mat.mechanical.density = 4700.0;
    mat.mechanical.youngs_modulus = GPaToPa(180.0);
    mat.mechanical.poisson_ratio = 0.27;
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition Rutile() {
    MaterialDefinition mat;
    mat.name = "Rutile";
    mat.category = "Minerals";
    mat.subcategory = "Oxides";
    mat.description = "TiO2";
    mat.mechanical.density = 4250.0;
    mat.mechanical.youngs_modulus = GPaToPa(280.0);
    mat.mechanical.poisson_ratio = 0.26;
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition Corundum() {
    MaterialDefinition mat;
    mat.name = "Corundum";
    mat.category = "Minerals";
    mat.subcategory = "Oxides";
    mat.description = "Al2O3 (sapphire/ruby)";
    mat.mechanical.density = 4000.0;
    mat.mechanical.youngs_modulus = GPaToPa(400.0);
    mat.mechanical.poisson_ratio = 0.23;
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition Spinel() {
    MaterialDefinition mat;
    mat.name = "Spinel";
    mat.category = "Minerals";
    mat.subcategory = "Oxides";
    mat.description = "MgAl2O4";
    mat.mechanical.density = 3580.0;
    mat.mechanical.youngs_modulus = GPaToPa(280.0);
    mat.mechanical.poisson_ratio = 0.26;
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition Chromite() {
    MaterialDefinition mat;
    mat.name = "Chromite";
    mat.category = "Minerals";
    mat.subcategory = "Oxides";
    mat.description = "FeCr2O4";
    mat.mechanical.density = 4800.0;
    mat.mechanical.youngs_modulus = GPaToPa(230.0);
    mat.mechanical.poisson_ratio = 0.28;
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition Limonite() {
    MaterialDefinition mat;
    mat.name = "Limonite";
    mat.category = "Minerals";
    mat.subcategory = "Oxides";
    mat.description = "FeO(OH)·nH2O (amorphous)";
    mat.mechanical.density = 4000.0;
    mat.mechanical.youngs_modulus = GPaToPa(100.0);
    mat.mechanical.poisson_ratio = 0.30;
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition Goethite() {
    MaterialDefinition mat;
    mat.name = "Goethite";
    mat.category = "Minerals";
    mat.subcategory = "Oxides";
    mat.description = "FeO(OH)";
    mat.mechanical.density = 4300.0;
    mat.mechanical.youngs_modulus = GPaToPa(120.0);
    mat.mechanical.poisson_ratio = 0.29;
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

// Sulfides
MaterialDefinition Pyrite() {
    MaterialDefinition mat;
    mat.name = "Pyrite";
    mat.category = "Minerals";
    mat.subcategory = "Sulfides";
    mat.description = "FeS2 (fool's gold)";
    mat.mechanical.density = 5020.0;
    mat.mechanical.youngs_modulus = GPaToPa(300.0);
    mat.mechanical.poisson_ratio = 0.16;
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition Pyrrhotite() {
    MaterialDefinition mat;
    mat.name = "Pyrrhotite";
    mat.category = "Minerals";
    mat.subcategory = "Sulfides";
    mat.description = "Fe(1-x)S";
    mat.mechanical.density = 4600.0;
    mat.mechanical.youngs_modulus = GPaToPa(150.0);
    mat.mechanical.poisson_ratio = 0.28;
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition Galena() {
    MaterialDefinition mat;
    mat.name = "Galena";
    mat.category = "Minerals";
    mat.subcategory = "Sulfides";
    mat.description = "PbS";
    mat.mechanical.density = 7600.0;
    mat.mechanical.youngs_modulus = GPaToPa(85.0);
    mat.mechanical.poisson_ratio = 0.30;
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition Sphalerite() {
    MaterialDefinition mat;
    mat.name = "Sphalerite";
    mat.category = "Minerals";
    mat.subcategory = "Sulfides";
    mat.description = "(Zn,Fe)S";
    mat.mechanical.density = 4050.0;
    mat.mechanical.youngs_modulus = GPaToPa(120.0);
    mat.mechanical.poisson_ratio = 0.28;
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition Chalcopyrite() {
    MaterialDefinition mat;
    mat.name = "Chalcopyrite";
    mat.category = "Minerals";
    mat.subcategory = "Sulfides";
    mat.description = "CuFeS2";
    mat.mechanical.density = 4200.0;
    mat.mechanical.youngs_modulus = GPaToPa(140.0);
    mat.mechanical.poisson_ratio = 0.27;
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition Molybdenite() {
    MaterialDefinition mat;
    mat.name = "Molybdenite";
    mat.category = "Minerals";
    mat.subcategory = "Sulfides";
    mat.description = "MoS2";
    mat.mechanical.density = 4700.0;
    mat.mechanical.youngs_modulus = GPaToPa(40.0);
    mat.mechanical.poisson_ratio = 0.25;
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

// Native elements
MaterialDefinition GraphiteMineral() {
    MaterialDefinition mat;
    mat.name = "Graphite";
    mat.category = "Minerals";
    mat.subcategory = "Native";
    mat.description = "Native carbon (layered)";
    mat.mechanical.density = 2200.0;
    mat.mechanical.youngs_modulus = GPaToPa(15.0);
    mat.mechanical.poisson_ratio = 0.30;
    mat.mechanical.computeDerivedProperties();
    mat.thermal.thermal_conductivity = 150.0;
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.is_anisotropic_flag = true;
    mat.recommended_model = "ANISOTROPIC_VTI";
    return mat;
}

MaterialDefinition SulfurNative() {
    MaterialDefinition mat;
    mat.name = "Sulfur";
    mat.category = "Minerals";
    mat.subcategory = "Native";
    mat.description = "Native sulfur";
    mat.mechanical.density = 2070.0;
    mat.mechanical.youngs_modulus = GPaToPa(15.0);
    mat.mechanical.poisson_ratio = 0.28;
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition GoldNative() {
    MaterialDefinition mat;
    mat.name = "Gold";
    mat.category = "Minerals";
    mat.subcategory = "Native";
    mat.description = "Native gold";
    mat.mechanical.density = 19300.0;
    mat.mechanical.youngs_modulus = GPaToPa(78.0);
    mat.mechanical.poisson_ratio = 0.44;
    mat.mechanical.computeDerivedProperties();
    mat.thermal.thermal_conductivity = 318.0;
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition CopperNative() {
    MaterialDefinition mat;
    mat.name = "Copper";
    mat.category = "Minerals";
    mat.subcategory = "Native";
    mat.description = "Native copper";
    mat.mechanical.density = 8960.0;
    mat.mechanical.youngs_modulus = GPaToPa(128.0);
    mat.mechanical.poisson_ratio = 0.34;
    mat.mechanical.computeDerivedProperties();
    mat.thermal.thermal_conductivity = 401.0;
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

// Clay minerals
MaterialDefinition Kaolinite() {
    MaterialDefinition mat;
    mat.name = "Kaolinite";
    mat.category = "Minerals";
    mat.subcategory = "Clays";
    mat.description = "Al2Si2O5(OH)4";
    mat.mechanical.density = 2600.0;
    mat.mechanical.youngs_modulus = GPaToPa(30.0);
    mat.mechanical.poisson_ratio = 0.32;
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition Illite() {
    MaterialDefinition mat;
    mat.name = "Illite";
    mat.category = "Minerals";
    mat.subcategory = "Clays";
    mat.description = "K,H3O(Al,Mg,Fe)2(Si,Al)4O10[(OH)2,(H2O)]";
    mat.mechanical.density = 2750.0;
    mat.mechanical.youngs_modulus = GPaToPa(35.0);
    mat.mechanical.poisson_ratio = 0.30;
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition Montmorillonite() {
    MaterialDefinition mat;
    mat.name = "Montmorillonite";
    mat.category = "Minerals";
    mat.subcategory = "Clays";
    mat.description = "(Na,Ca)0.33(Al,Mg)2(Si4O10)(OH)2·nH2O";
    mat.mechanical.density = 2500.0;
    mat.mechanical.youngs_modulus = GPaToPa(15.0);
    mat.mechanical.poisson_ratio = 0.35;
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition Smectite() {
    MaterialDefinition mat;
    mat.name = "Smectite";
    mat.category = "Minerals";
    mat.subcategory = "Clays";
    mat.description = "Smectite group clays (swelling)";
    mat.mechanical.density = 2400.0;
    mat.mechanical.youngs_modulus = GPaToPa(12.0);
    mat.mechanical.poisson_ratio = 0.36;
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition Vermiculite() {
    MaterialDefinition mat;
    mat.name = "Vermiculite";
    mat.category = "Minerals";
    mat.subcategory = "Clays";
    mat.description = "Hydrous phyllosilicate";
    mat.mechanical.density = 2500.0;
    mat.mechanical.youngs_modulus = GPaToPa(18.0);
    mat.mechanical.poisson_ratio = 0.32;
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition Palygorskite() {
    MaterialDefinition mat;
    mat.name = "Palygorskite";
    mat.category = "Minerals";
    mat.subcategory = "Clays";
    mat.description = "(Mg,Al)2Si4O10(OH)·4H2O (attapulgite)";
    mat.mechanical.density = 2300.0;
    mat.mechanical.youngs_modulus = GPaToPa(15.0);
    mat.mechanical.poisson_ratio = 0.30;
    mat.mechanical.computeDerivedProperties();
    mat.hydraulic.porosity = 0.40;
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.is_porous = true;
    mat.recommended_model = "POROELASTIC";
    return mat;
}

} // namespace Minerals

// =============================================================================
// Fault Zone Materials - namespace Materials::FaultMaterials
// =============================================================================

namespace FaultMaterials {

MaterialDefinition FaultGouge() {
    MaterialDefinition mat;
    mat.name = "Fault_Gouge";
    mat.category = "FaultMaterials";
    mat.subcategory = "Gouge";
    mat.description = "Fine-grained fault gouge";
    mat.reference = "Sibson, 1977";
    
    mat.mechanical.density = 2200.0;
    mat.mechanical.youngs_modulus = GPaToPa(5.0);
    mat.mechanical.poisson_ratio = 0.30;
    mat.mechanical.tensile_strength = MPaToPa(0.5);
    mat.mechanical.compressive_strength = MPaToPa(10.0);
    mat.mechanical.cohesion = MPaToPa(1.0);
    mat.mechanical.friction_angle = 25.0;
    mat.mechanical.biot_coefficient = 0.90;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.20;
    mat.hydraulic.setIsotropicPermeability(mDToM2(10.0));
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 15.0;
    mat.seismic.Qs = 8.0;
    
    mat.is_porous = true;
    mat.is_fractured = true;
    mat.recommended_model = "ELASTOPLASTIC_MC";
    
    return mat;
}

MaterialDefinition FaultGougeClayRich() {
    MaterialDefinition mat = FaultGouge();
    mat.name = "Fault_Gouge_ClayRich";
    mat.description = "Clay-rich fault gouge";
    mat.mechanical.friction_angle = 15.0;
    mat.mechanical.cohesion = MPaToPa(0.5);
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.01));
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    return mat;
}

MaterialDefinition FaultBreccia() {
    MaterialDefinition mat;
    mat.name = "Fault_Breccia";
    mat.category = "FaultMaterials";
    mat.subcategory = "Breccia";
    mat.description = "Fault breccia with angular clasts";
    mat.reference = "Sibson, 1977";
    
    mat.mechanical.density = 2400.0;
    mat.mechanical.youngs_modulus = GPaToPa(15.0);
    mat.mechanical.poisson_ratio = 0.28;
    mat.mechanical.tensile_strength = MPaToPa(2.0);
    mat.mechanical.compressive_strength = MPaToPa(40.0);
    mat.mechanical.cohesion = MPaToPa(5.0);
    mat.mechanical.friction_angle = 35.0;
    mat.mechanical.biot_coefficient = 0.75;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.15;
    mat.hydraulic.setIsotropicPermeability(mDToM2(100.0));
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 25.0;
    mat.seismic.Qs = 12.0;
    
    mat.is_porous = true;
    mat.is_fractured = true;
    mat.recommended_model = "POROELASTIC";
    
    return mat;
}

MaterialDefinition FaultCataclasite() {
    MaterialDefinition mat;
    mat.name = "Fault_Cataclasite";
    mat.category = "FaultMaterials";
    mat.subcategory = "Cataclasite";
    mat.description = "Cataclastic fault rock";
    mat.mechanical.density = 2500.0;
    mat.mechanical.youngs_modulus = GPaToPa(20.0);
    mat.mechanical.poisson_ratio = 0.26;
    mat.mechanical.cohesion = MPaToPa(8.0);
    mat.mechanical.friction_angle = 38.0;
    mat.mechanical.biot_coefficient = 0.65;
    mat.mechanical.computeDerivedProperties();
    mat.hydraulic.porosity = 0.10;
    mat.hydraulic.setIsotropicPermeability(mDToM2(5.0));
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.is_porous = true;
    mat.recommended_model = "POROELASTIC";
    return mat;
}

MaterialDefinition FaultMylonite() {
    MaterialDefinition mat;
    mat.name = "Fault_Mylonite";
    mat.category = "FaultMaterials";
    mat.subcategory = "Mylonite";
    mat.description = "Ductile fault rock";
    mat.mechanical.density = 2700.0;
    mat.mechanical.youngs_modulus = GPaToPa(35.0);
    mat.mechanical.poisson_ratio = 0.25;
    mat.mechanical.cohesion = MPaToPa(15.0);
    mat.mechanical.friction_angle = 40.0;
    mat.mechanical.biot_coefficient = 0.50;
    mat.mechanical.computeDerivedProperties();
    mat.hydraulic.porosity = 0.05;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.1));
    mat.anisotropic.is_anisotropic = true;
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.is_anisotropic_flag = true;
    mat.recommended_model = "ANISOTROPIC_VTI";
    return mat;
}

MaterialDefinition FaultPseudotachylyte() {
    MaterialDefinition mat;
    mat.name = "Fault_Pseudotachylyte";
    mat.category = "FaultMaterials";
    mat.subcategory = "Melt";
    mat.description = "Friction melt rock from seismic slip";
    mat.mechanical.density = 2750.0;
    mat.mechanical.youngs_modulus = GPaToPa(60.0);
    mat.mechanical.poisson_ratio = 0.22;
    mat.mechanical.cohesion = MPaToPa(30.0);
    mat.mechanical.friction_angle = 50.0;
    mat.mechanical.biot_coefficient = 0.30;
    mat.mechanical.computeDerivedProperties();
    mat.hydraulic.porosity = 0.01;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.001));
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition FaultDamageZone() {
    MaterialDefinition mat;
    mat.name = "Fault_DamageZone";
    mat.category = "FaultMaterials";
    mat.subcategory = "DamageZone";
    mat.description = "Fractured rock around fault";
    mat.mechanical.density = 2500.0;
    mat.mechanical.youngs_modulus = GPaToPa(25.0);
    mat.mechanical.poisson_ratio = 0.27;
    mat.mechanical.cohesion = MPaToPa(10.0);
    mat.mechanical.friction_angle = 40.0;
    mat.mechanical.biot_coefficient = 0.70;
    mat.mechanical.computeDerivedProperties();
    mat.hydraulic.porosity = 0.08;
    mat.hydraulic.setIsotropicPermeability(mDToM2(50.0));
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.is_porous = true;
    mat.is_fractured = true;
    mat.recommended_model = "POROELASTIC";
    return mat;
}

MaterialDefinition FaultCore() {
    MaterialDefinition mat;
    mat.name = "Fault_Core";
    mat.category = "FaultMaterials";
    mat.subcategory = "Core";
    mat.description = "Central fault zone";
    mat.mechanical.density = 2300.0;
    mat.mechanical.youngs_modulus = GPaToPa(8.0);
    mat.mechanical.poisson_ratio = 0.30;
    mat.mechanical.cohesion = MPaToPa(2.0);
    mat.mechanical.friction_angle = 28.0;
    mat.mechanical.biot_coefficient = 0.85;
    mat.mechanical.computeDerivedProperties();
    mat.hydraulic.porosity = 0.18;
    mat.hydraulic.setIsotropicPermeability(mDToM2(1.0));
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.is_porous = true;
    mat.is_fractured = true;
    mat.recommended_model = "ELASTOPLASTIC_MC";
    return mat;
}

MaterialDefinition SlipZone() {
    MaterialDefinition mat;
    mat.name = "Slip_Zone";
    mat.category = "FaultMaterials";
    mat.subcategory = "Slip";
    mat.description = "Active slip surface material";
    mat.mechanical.density = 2200.0;
    mat.mechanical.youngs_modulus = GPaToPa(3.0);
    mat.mechanical.poisson_ratio = 0.32;
    mat.mechanical.cohesion = MPaToPa(0.2);
    mat.mechanical.friction_angle = 18.0;
    mat.mechanical.biot_coefficient = 0.95;
    mat.mechanical.computeDerivedProperties();
    mat.hydraulic.porosity = 0.25;
    mat.hydraulic.setIsotropicPermeability(mDToM2(5.0));
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.is_porous = true;
    mat.is_fractured = true;
    mat.recommended_model = "ELASTOPLASTIC_MC";
    return mat;
}

MaterialDefinition ShearedZone() {
    MaterialDefinition mat;
    mat.name = "Sheared_Zone";
    mat.category = "FaultMaterials";
    mat.subcategory = "Shear";
    mat.description = "Sheared rock material";
    mat.mechanical.density = 2400.0;
    mat.mechanical.youngs_modulus = GPaToPa(12.0);
    mat.mechanical.poisson_ratio = 0.28;
    mat.mechanical.cohesion = MPaToPa(4.0);
    mat.mechanical.friction_angle = 32.0;
    mat.mechanical.biot_coefficient = 0.75;
    mat.mechanical.computeDerivedProperties();
    mat.hydraulic.porosity = 0.12;
    mat.hydraulic.setIsotropicPermeability(mDToM2(10.0));
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.is_porous = true;
    mat.is_fractured = true;
    mat.recommended_model = "POROELASTIC";
    return mat;
}

} // namespace FaultMaterials

// =============================================================================
// Wellbore Materials - namespace Materials::Wellbore
// =============================================================================

namespace Wellbore {

MaterialDefinition CementClassA() {
    MaterialDefinition mat;
    mat.name = "Cement_ClassA";
    mat.category = "Wellbore";
    mat.subcategory = "Cement";
    mat.description = "API Class A Portland cement";
    mat.reference = "API RP 10B";
    
    mat.mechanical.density = 1900.0;
    mat.mechanical.youngs_modulus = GPaToPa(10.0);
    mat.mechanical.poisson_ratio = 0.20;
    mat.mechanical.tensile_strength = MPaToPa(3.0);
    mat.mechanical.compressive_strength = MPaToPa(30.0);
    mat.mechanical.cohesion = MPaToPa(5.0);
    mat.mechanical.friction_angle = 35.0;
    mat.mechanical.fracture_toughness_KIc = 0.5e6;
    mat.mechanical.biot_coefficient = 0.60;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.15;
    mat.hydraulic.setIsotropicPermeability(mDToM2(0.01));
    
    mat.thermal.thermal_conductivity = 1.0;
    mat.thermal.specific_heat = 880.0;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.is_porous = true;
    mat.recommended_model = "POROELASTIC";
    
    return mat;
}

MaterialDefinition CementClassC() {
    MaterialDefinition mat = CementClassA();
    mat.name = "Cement_ClassC";
    mat.description = "API Class C high early strength cement";
    mat.mechanical.compressive_strength = MPaToPa(35.0);
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    return mat;
}

MaterialDefinition CementClassG() {
    MaterialDefinition mat = CementClassA();
    mat.name = "Cement_ClassG";
    mat.description = "API Class G basic oil well cement";
    mat.mechanical.density = 1920.0;
    mat.mechanical.youngs_modulus = GPaToPa(12.0);
    mat.mechanical.compressive_strength = MPaToPa(35.0);
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    return mat;
}

MaterialDefinition CementClassH() {
    MaterialDefinition mat = CementClassA();
    mat.name = "Cement_ClassH";
    mat.description = "API Class H basic oil well cement (coarser)";
    mat.mechanical.density = 1980.0;
    mat.mechanical.youngs_modulus = GPaToPa(14.0);
    mat.mechanical.compressive_strength = MPaToPa(40.0);
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    return mat;
}

MaterialDefinition CementFoamed() {
    MaterialDefinition mat = CementClassA();
    mat.name = "Cement_Foamed";
    mat.description = "Lightweight foamed cement";
    mat.mechanical.density = 1200.0;
    mat.mechanical.youngs_modulus = GPaToPa(3.0);
    mat.mechanical.compressive_strength = MPaToPa(8.0);
    mat.mechanical.biot_coefficient = 0.80;
    mat.mechanical.computeDerivedProperties();
    mat.hydraulic.porosity = 0.45;
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    return mat;
}

MaterialDefinition CementLightweight() {
    MaterialDefinition mat = CementClassA();
    mat.name = "Cement_Lightweight";
    mat.description = "Lightweight cement with additives";
    mat.mechanical.density = 1400.0;
    mat.mechanical.youngs_modulus = GPaToPa(5.0);
    mat.mechanical.compressive_strength = MPaToPa(12.0);
    mat.mechanical.computeDerivedProperties();
    mat.hydraulic.porosity = 0.30;
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    return mat;
}

MaterialDefinition CementHighDensity() {
    MaterialDefinition mat = CementClassA();
    mat.name = "Cement_HighDensity";
    mat.description = "High density cement with weighting agents";
    mat.mechanical.density = 2300.0;
    mat.mechanical.youngs_modulus = GPaToPa(18.0);
    mat.mechanical.compressive_strength = MPaToPa(50.0);
    mat.mechanical.biot_coefficient = 0.45;
    mat.mechanical.computeDerivedProperties();
    mat.hydraulic.porosity = 0.08;
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    return mat;
}

MaterialDefinition SteelCasing() {
    MaterialDefinition mat;
    mat.name = "Steel_Casing";
    mat.category = "Wellbore";
    mat.subcategory = "Casing";
    mat.description = "Standard oilfield casing steel";
    mat.reference = "API 5CT";
    
    mat.mechanical.density = 7850.0;
    mat.mechanical.youngs_modulus = GPaToPa(207.0);
    mat.mechanical.poisson_ratio = 0.30;
    mat.mechanical.tensile_strength = MPaToPa(550.0);
    mat.mechanical.compressive_strength = MPaToPa(550.0);
    mat.mechanical.biot_coefficient = 0.0;
    mat.mechanical.computeDerivedProperties();
    
    mat.thermal.thermal_conductivity = 50.0;
    mat.thermal.specific_heat = 500.0;
    mat.thermal.thermal_expansion = 1.2e-5;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    
    return mat;
}

MaterialDefinition SteelTubing() {
    MaterialDefinition mat = SteelCasing();
    mat.name = "Steel_Tubing";
    mat.description = "Production tubing steel";
    mat.mechanical.tensile_strength = MPaToPa(520.0);
    mat.mechanical.compressive_strength = MPaToPa(520.0);
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    return mat;
}

MaterialDefinition SteelDrillPipe() {
    MaterialDefinition mat = SteelCasing();
    mat.name = "Steel_DrillPipe";
    mat.description = "Drill pipe steel";
    mat.mechanical.tensile_strength = MPaToPa(690.0);
    mat.mechanical.compressive_strength = MPaToPa(690.0);
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    return mat;
}

MaterialDefinition Proppant() {
    MaterialDefinition mat;
    mat.name = "Proppant";
    mat.category = "Wellbore";
    mat.subcategory = "Proppant";
    mat.description = "Generic proppant pack";
    mat.mechanical.density = 2650.0;
    mat.mechanical.youngs_modulus = GPaToPa(20.0);
    mat.mechanical.poisson_ratio = 0.25;
    mat.mechanical.biot_coefficient = 0.90;
    mat.mechanical.computeDerivedProperties();
    mat.hydraulic.porosity = 0.35;
    mat.hydraulic.setIsotropicPermeability(mDToM2(100000.0));
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.is_porous = true;
    mat.recommended_model = "POROELASTIC";
    return mat;
}

MaterialDefinition ProppantSand() {
    MaterialDefinition mat = Proppant();
    mat.name = "Proppant_Sand";
    mat.description = "Natural sand proppant";
    mat.mechanical.density = 2650.0;
    mat.mechanical.compressive_strength = MPaToPa(50.0);
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    return mat;
}

MaterialDefinition ProppantCeramic() {
    MaterialDefinition mat = Proppant();
    mat.name = "Proppant_Ceramic";
    mat.description = "Ceramic proppant, high strength";
    mat.mechanical.density = 3200.0;
    mat.mechanical.youngs_modulus = GPaToPa(50.0);
    mat.mechanical.compressive_strength = MPaToPa(100.0);
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    return mat;
}

MaterialDefinition ProppantResinCoated() {
    MaterialDefinition mat = Proppant();
    mat.name = "Proppant_ResinCoated";
    mat.description = "Resin-coated sand proppant";
    mat.mechanical.density = 2550.0;
    mat.mechanical.youngs_modulus = GPaToPa(25.0);
    mat.mechanical.compressive_strength = MPaToPa(65.0);
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    return mat;
}

MaterialDefinition DrillCuttings() {
    MaterialDefinition mat;
    mat.name = "Drill_Cuttings";
    mat.category = "Wellbore";
    mat.subcategory = "Cuttings";
    mat.description = "Drilled rock cuttings";
    mat.mechanical.density = 2000.0;
    mat.mechanical.youngs_modulus = GPaToPa(0.5);
    mat.mechanical.poisson_ratio = 0.30;
    mat.mechanical.biot_coefficient = 1.0;
    mat.mechanical.computeDerivedProperties();
    mat.hydraulic.porosity = 0.40;
    mat.hydraulic.setIsotropicPermeability(mDToM2(10000.0));
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.is_porous = true;
    mat.recommended_model = "ELASTOPLASTIC_MC";
    return mat;
}

} // namespace Wellbore

// =============================================================================
// Reference Materials - namespace Materials::Reference
// =============================================================================

namespace Reference {

MaterialDefinition Steel() {
    MaterialDefinition mat;
    mat.name = "Steel";
    mat.category = "Reference";
    mat.subcategory = "Metals";
    mat.description = "Generic structural steel";
    mat.mechanical.density = 7850.0;
    mat.mechanical.youngs_modulus = GPaToPa(200.0);
    mat.mechanical.poisson_ratio = 0.30;
    mat.mechanical.tensile_strength = MPaToPa(400.0);
    mat.mechanical.biot_coefficient = 0.0;
    mat.mechanical.computeDerivedProperties();
    mat.thermal.thermal_conductivity = 50.0;
    mat.thermal.specific_heat = 500.0;
    mat.thermal.thermal_expansion = 1.2e-5;
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition SteelCarbon() {
    MaterialDefinition mat = Steel();
    mat.name = "Steel_Carbon";
    mat.description = "Carbon steel";
    mat.mechanical.tensile_strength = MPaToPa(450.0);
    return mat;
}

MaterialDefinition SteelStainless() {
    MaterialDefinition mat = Steel();
    mat.name = "Steel_Stainless";
    mat.description = "Stainless steel 316";
    mat.mechanical.density = 8000.0;
    mat.mechanical.youngs_modulus = GPaToPa(193.0);
    mat.mechanical.tensile_strength = MPaToPa(580.0);
    mat.mechanical.computeDerivedProperties();
    mat.thermal.thermal_conductivity = 16.0;
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    return mat;
}

MaterialDefinition Aluminum() {
    MaterialDefinition mat;
    mat.name = "Aluminum";
    mat.category = "Reference";
    mat.subcategory = "Metals";
    mat.description = "Aluminum alloy";
    mat.mechanical.density = 2700.0;
    mat.mechanical.youngs_modulus = GPaToPa(70.0);
    mat.mechanical.poisson_ratio = 0.33;
    mat.mechanical.tensile_strength = MPaToPa(300.0);
    mat.mechanical.biot_coefficient = 0.0;
    mat.mechanical.computeDerivedProperties();
    mat.thermal.thermal_conductivity = 237.0;
    mat.thermal.specific_heat = 900.0;
    mat.thermal.thermal_expansion = 2.3e-5;
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition Concrete() {
    MaterialDefinition mat;
    mat.name = "Concrete";
    mat.category = "Reference";
    mat.subcategory = "Construction";
    mat.description = "Standard concrete";
    mat.mechanical.density = 2400.0;
    mat.mechanical.youngs_modulus = GPaToPa(30.0);
    mat.mechanical.poisson_ratio = 0.20;
    mat.mechanical.tensile_strength = MPaToPa(3.0);
    mat.mechanical.compressive_strength = MPaToPa(30.0);
    mat.mechanical.biot_coefficient = 0.50;
    mat.mechanical.computeDerivedProperties();
    mat.thermal.thermal_conductivity = 1.7;
    mat.thermal.specific_heat = 880.0;
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.is_porous = true;
    mat.recommended_model = "POROELASTIC";
    return mat;
}

MaterialDefinition ConcreteHighStrength() {
    MaterialDefinition mat = Concrete();
    mat.name = "Concrete_HighStrength";
    mat.description = "High strength concrete";
    mat.mechanical.youngs_modulus = GPaToPa(45.0);
    mat.mechanical.compressive_strength = MPaToPa(80.0);
    mat.mechanical.tensile_strength = MPaToPa(5.0);
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    return mat;
}

MaterialDefinition Glass() {
    MaterialDefinition mat;
    mat.name = "Glass";
    mat.category = "Reference";
    mat.subcategory = "Brittle";
    mat.description = "Soda-lime glass";
    mat.mechanical.density = 2500.0;
    mat.mechanical.youngs_modulus = GPaToPa(70.0);
    mat.mechanical.poisson_ratio = 0.22;
    mat.mechanical.tensile_strength = MPaToPa(50.0);
    mat.mechanical.compressive_strength = MPaToPa(1000.0);
    mat.mechanical.biot_coefficient = 0.0;
    mat.mechanical.computeDerivedProperties();
    mat.thermal.thermal_conductivity = 1.0;
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition Acrylic() {
    MaterialDefinition mat;
    mat.name = "Acrylic";
    mat.category = "Reference";
    mat.subcategory = "Polymers";
    mat.description = "Acrylic (PMMA)";
    mat.mechanical.density = 1180.0;
    mat.mechanical.youngs_modulus = GPaToPa(3.2);
    mat.mechanical.poisson_ratio = 0.37;
    mat.mechanical.tensile_strength = MPaToPa(70.0);
    mat.mechanical.biot_coefficient = 0.0;
    mat.mechanical.computeDerivedProperties();
    mat.thermal.thermal_conductivity = 0.19;
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition PMMA() {
    return Acrylic();
}

MaterialDefinition Copper() {
    MaterialDefinition mat;
    mat.name = "Copper_Metal";
    mat.category = "Reference";
    mat.subcategory = "Metals";
    mat.description = "Pure copper";
    mat.mechanical.density = 8960.0;
    mat.mechanical.youngs_modulus = GPaToPa(128.0);
    mat.mechanical.poisson_ratio = 0.34;
    mat.mechanical.tensile_strength = MPaToPa(220.0);
    mat.mechanical.biot_coefficient = 0.0;
    mat.mechanical.computeDerivedProperties();
    mat.thermal.thermal_conductivity = 401.0;
    mat.thermal.specific_heat = 385.0;
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition Titanium() {
    MaterialDefinition mat;
    mat.name = "Titanium";
    mat.category = "Reference";
    mat.subcategory = "Metals";
    mat.description = "Titanium alloy";
    mat.mechanical.density = 4500.0;
    mat.mechanical.youngs_modulus = GPaToPa(116.0);
    mat.mechanical.poisson_ratio = 0.32;
    mat.mechanical.tensile_strength = MPaToPa(900.0);
    mat.mechanical.biot_coefficient = 0.0;
    mat.mechanical.computeDerivedProperties();
    mat.thermal.thermal_conductivity = 22.0;
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition InconelAlloy() {
    MaterialDefinition mat;
    mat.name = "Inconel_625";
    mat.category = "Reference";
    mat.subcategory = "Metals";
    mat.description = "Inconel 625 nickel alloy";
    mat.mechanical.density = 8440.0;
    mat.mechanical.youngs_modulus = GPaToPa(208.0);
    mat.mechanical.poisson_ratio = 0.30;
    mat.mechanical.tensile_strength = MPaToPa(900.0);
    mat.mechanical.biot_coefficient = 0.0;
    mat.mechanical.computeDerivedProperties();
    mat.thermal.thermal_conductivity = 9.8;
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition Water() {
    MaterialDefinition mat;
    mat.name = "Water";
    mat.category = "Reference";
    mat.subcategory = "Fluids";
    mat.description = "Pure water at 20°C";
    mat.mechanical.density = 1000.0;
    mat.mechanical.youngs_modulus = GPaToPa(2.2);
    mat.mechanical.poisson_ratio = 0.5;
    mat.mechanical.biot_coefficient = 1.0;
    mat.mechanical.computeDerivedProperties();
    mat.thermal.thermal_conductivity = 0.6;
    mat.thermal.specific_heat = 4186.0;
    mat.seismic.Vp = 1500.0;
    mat.seismic.Vs = 0.0;
    mat.seismic.Qp = 10000.0;
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition Ice() {
    MaterialDefinition mat;
    mat.name = "Ice";
    mat.category = "Reference";
    mat.subcategory = "Other";
    mat.description = "Ice Ih at 0°C";
    mat.mechanical.density = 917.0;
    mat.mechanical.youngs_modulus = GPaToPa(9.3);
    mat.mechanical.poisson_ratio = 0.33;
    mat.mechanical.tensile_strength = MPaToPa(1.0);
    mat.mechanical.biot_coefficient = 0.0;
    mat.mechanical.computeDerivedProperties();
    mat.thermal.thermal_conductivity = 2.2;
    mat.thermal.specific_heat = 2090.0;
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition Air() {
    MaterialDefinition mat;
    mat.name = "Air";
    mat.category = "Reference";
    mat.subcategory = "Fluids";
    mat.description = "Air at STP";
    mat.mechanical.density = 1.225;
    mat.mechanical.youngs_modulus = 1.42e5;
    mat.mechanical.poisson_ratio = 0.5;
    mat.mechanical.biot_coefficient = 1.0;
    mat.mechanical.computeDerivedProperties();
    mat.thermal.thermal_conductivity = 0.026;
    mat.thermal.specific_heat = 1005.0;
    mat.seismic.Vp = 343.0;
    mat.seismic.Vs = 0.0;
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition Rubber() {
    MaterialDefinition mat;
    mat.name = "Rubber";
    mat.category = "Reference";
    mat.subcategory = "Polymers";
    mat.description = "Natural rubber";
    mat.mechanical.density = 1100.0;
    mat.mechanical.youngs_modulus = GPaToPa(0.01);
    mat.mechanical.poisson_ratio = 0.49;
    mat.mechanical.tensile_strength = MPaToPa(25.0);
    mat.mechanical.biot_coefficient = 0.0;
    mat.mechanical.computeDerivedProperties();
    mat.thermal.thermal_conductivity = 0.15;
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.is_viscoelastic = true;
    mat.recommended_model = "VISCOELASTIC_SLS";
    return mat;
}

MaterialDefinition Epoxy() {
    MaterialDefinition mat;
    mat.name = "Epoxy";
    mat.category = "Reference";
    mat.subcategory = "Polymers";
    mat.description = "Epoxy resin";
    mat.mechanical.density = 1200.0;
    mat.mechanical.youngs_modulus = GPaToPa(3.5);
    mat.mechanical.poisson_ratio = 0.35;
    mat.mechanical.tensile_strength = MPaToPa(60.0);
    mat.mechanical.biot_coefficient = 0.0;
    mat.mechanical.computeDerivedProperties();
    mat.thermal.thermal_conductivity = 0.2;
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition CarbonFiber() {
    MaterialDefinition mat;
    mat.name = "CarbonFiber";
    mat.category = "Reference";
    mat.subcategory = "Composites";
    mat.description = "Carbon fiber composite";
    mat.mechanical.density = 1600.0;
    mat.mechanical.youngs_modulus = GPaToPa(230.0);
    mat.mechanical.poisson_ratio = 0.30;
    mat.mechanical.tensile_strength = MPaToPa(3500.0);
    mat.mechanical.biot_coefficient = 0.0;
    mat.mechanical.computeDerivedProperties();
    mat.thermal.thermal_conductivity = 20.0;
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.is_anisotropic_flag = true;
    mat.recommended_model = "ANISOTROPIC_VTI";
    return mat;
}

MaterialDefinition Ceramic() {
    MaterialDefinition mat;
    mat.name = "Ceramic";
    mat.category = "Reference";
    mat.subcategory = "Brittle";
    mat.description = "Alumina ceramic";
    mat.mechanical.density = 3900.0;
    mat.mechanical.youngs_modulus = GPaToPa(380.0);
    mat.mechanical.poisson_ratio = 0.22;
    mat.mechanical.tensile_strength = MPaToPa(300.0);
    mat.mechanical.compressive_strength = MPaToPa(2500.0);
    mat.mechanical.biot_coefficient = 0.0;
    mat.mechanical.computeDerivedProperties();
    mat.thermal.thermal_conductivity = 30.0;
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

} // namespace Reference

// =============================================================================
// Planetary Materials - namespace Materials::Planetary
// =============================================================================

namespace Planetary {

MaterialDefinition LunarRegolith() {
    MaterialDefinition mat;
    mat.name = "Lunar_Regolith";
    mat.category = "Planetary";
    mat.subcategory = "Lunar";
    mat.description = "Lunar surface regolith";
    mat.reference = "Carrier et al., 1991";
    
    mat.mechanical.density = 1500.0;
    mat.mechanical.youngs_modulus = GPaToPa(0.02);
    mat.mechanical.poisson_ratio = 0.30;
    mat.mechanical.cohesion = MPaToPa(0.001);
    mat.mechanical.friction_angle = 40.0;
    mat.mechanical.biot_coefficient = 1.0;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.45;
    mat.hydraulic.setIsotropicPermeability(mDToM2(100.0));
    
    mat.thermal.thermal_conductivity = 0.01;
    mat.thermal.specific_heat = 700.0;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 3000.0;
    
    mat.is_porous = true;
    mat.recommended_model = "ELASTOPLASTIC_MC";
    
    return mat;
}

MaterialDefinition LunarBasalt() {
    MaterialDefinition mat;
    mat.name = "Lunar_Basalt";
    mat.category = "Planetary";
    mat.subcategory = "Lunar";
    mat.description = "Lunar mare basalt";
    mat.mechanical.density = 3100.0;
    mat.mechanical.youngs_modulus = GPaToPa(80.0);
    mat.mechanical.poisson_ratio = 0.25;
    mat.mechanical.compressive_strength = MPaToPa(200.0);
    mat.mechanical.biot_coefficient = 0.30;
    mat.mechanical.computeDerivedProperties();
    mat.hydraulic.porosity = 0.01;
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 5000.0;
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition LunarHighlands() {
    MaterialDefinition mat;
    mat.name = "Lunar_Highlands";
    mat.category = "Planetary";
    mat.subcategory = "Lunar";
    mat.description = "Lunar highlands anorthosite";
    mat.mechanical.density = 2800.0;
    mat.mechanical.youngs_modulus = GPaToPa(70.0);
    mat.mechanical.poisson_ratio = 0.26;
    mat.mechanical.biot_coefficient = 0.35;
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.seismic.Qp = 4000.0;
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition MartianRegolith() {
    MaterialDefinition mat;
    mat.name = "Martian_Regolith";
    mat.category = "Planetary";
    mat.subcategory = "Mars";
    mat.description = "Martian surface regolith";
    mat.reference = "Moore et al., 1999";
    
    mat.mechanical.density = 1200.0;
    mat.mechanical.youngs_modulus = GPaToPa(0.015);
    mat.mechanical.poisson_ratio = 0.32;
    mat.mechanical.cohesion = MPaToPa(0.005);
    mat.mechanical.friction_angle = 35.0;
    mat.mechanical.biot_coefficient = 1.0;
    mat.mechanical.computeDerivedProperties();
    
    mat.hydraulic.porosity = 0.50;
    mat.hydraulic.setIsotropicPermeability(mDToM2(500.0));
    
    mat.thermal.thermal_conductivity = 0.02;
    mat.thermal.specific_heat = 800.0;
    
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    
    mat.is_porous = true;
    mat.recommended_model = "ELASTOPLASTIC_MC";
    
    return mat;
}

MaterialDefinition MartianBasalt() {
    MaterialDefinition mat;
    mat.name = "Martian_Basalt";
    mat.category = "Planetary";
    mat.subcategory = "Mars";
    mat.description = "Martian volcanic basalt";
    mat.mechanical.density = 2900.0;
    mat.mechanical.youngs_modulus = GPaToPa(70.0);
    mat.mechanical.poisson_ratio = 0.26;
    mat.mechanical.biot_coefficient = 0.40;
    mat.mechanical.computeDerivedProperties();
    mat.hydraulic.porosity = 0.10;
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.is_porous = true;
    mat.recommended_model = "POROELASTIC";
    return mat;
}

MaterialDefinition Meteorite() {
    MaterialDefinition mat;
    mat.name = "Meteorite";
    mat.category = "Planetary";
    mat.subcategory = "Meteorite";
    mat.description = "Generic chondrite meteorite";
    mat.mechanical.density = 3500.0;
    mat.mechanical.youngs_modulus = GPaToPa(60.0);
    mat.mechanical.poisson_ratio = 0.28;
    mat.mechanical.biot_coefficient = 0.30;
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition MeteoeriteIron() {
    MaterialDefinition mat;
    mat.name = "Meteorite_Iron";
    mat.category = "Planetary";
    mat.subcategory = "Meteorite";
    mat.description = "Iron meteorite";
    mat.mechanical.density = 7800.0;
    mat.mechanical.youngs_modulus = GPaToPa(200.0);
    mat.mechanical.poisson_ratio = 0.30;
    mat.mechanical.biot_coefficient = 0.0;
    mat.mechanical.computeDerivedProperties();
    mat.thermal.thermal_conductivity = 40.0;
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition MeteoeriteStony() {
    MaterialDefinition mat;
    mat.name = "Meteorite_Stony";
    mat.category = "Planetary";
    mat.subcategory = "Meteorite";
    mat.description = "Stony meteorite (chondrite)";
    mat.mechanical.density = 3400.0;
    mat.mechanical.youngs_modulus = GPaToPa(55.0);
    mat.mechanical.poisson_ratio = 0.27;
    mat.mechanical.biot_coefficient = 0.35;
    mat.mechanical.computeDerivedProperties();
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.recommended_model = "LINEAR_ELASTIC";
    return mat;
}

MaterialDefinition AsteroidMaterial() {
    MaterialDefinition mat;
    mat.name = "Asteroid_Material";
    mat.category = "Planetary";
    mat.subcategory = "Asteroid";
    mat.description = "Generic asteroid rubble pile";
    mat.mechanical.density = 2000.0;
    mat.mechanical.youngs_modulus = GPaToPa(0.1);
    mat.mechanical.poisson_ratio = 0.30;
    mat.mechanical.cohesion = MPaToPa(0.0001);
    mat.mechanical.friction_angle = 35.0;
    mat.mechanical.biot_coefficient = 1.0;
    mat.mechanical.computeDerivedProperties();
    mat.hydraulic.porosity = 0.40;
    mat.seismic.computeFromElastic(mat.mechanical.youngs_modulus, mat.mechanical.poisson_ratio, mat.mechanical.density);
    mat.is_porous = true;
    mat.recommended_model = "ELASTOPLASTIC_MC";
    return mat;
}

} // namespace Planetary
} // namespace Materials

// =============================================================================
// Formation Material Definitions - namespace Formations
// =============================================================================

namespace Formations {

// North American Formations
MaterialDefinition BakkenShale() {
    auto mat = Materials::Sedimentary::ShaleOrganic();
    mat.name = "Bakken_Shale";
    mat.description = "Bakken Formation organic-rich shale, Williston Basin";
    return mat;
}

MaterialDefinition EagleFordShale() {
    auto mat = Materials::Sedimentary::ShaleCalcareous();
    mat.name = "Eagle_Ford_Shale";
    mat.description = "Eagle Ford Formation, South Texas";
    return mat;
}

MaterialDefinition MarcellUsShale() {
    auto mat = Materials::Sedimentary::ShaleOrganic();
    mat.name = "Marcellus_Shale";
    mat.description = "Marcellus Formation, Appalachian Basin";
    return mat;
}

MaterialDefinition UticaShale() {
    auto mat = Materials::Sedimentary::ShaleOrganic();
    mat.name = "Utica_Shale";
    mat.description = "Utica/Point Pleasant Formation";
    return mat;
}

MaterialDefinition HaynesvilleShale() {
    auto mat = Materials::Sedimentary::ShaleOrganic();
    mat.name = "Haynesville_Shale";
    mat.description = "Haynesville Formation, East Texas/Louisiana";
    return mat;
}

MaterialDefinition BarnnettShale() {
    auto mat = Materials::Sedimentary::ShaleSiliceous();
    mat.name = "Barnett_Shale";
    mat.description = "Barnett Formation, Fort Worth Basin";
    return mat;
}

MaterialDefinition WoodfordShale() {
    auto mat = Materials::Sedimentary::ShaleOrganic();
    mat.name = "Woodford_Shale";
    mat.description = "Woodford Formation, Oklahoma";
    return mat;
}

MaterialDefinition NiobraraChalk() {
    auto mat = Materials::Sedimentary::LimestoneChalk();
    mat.name = "Niobrara_Chalk";
    mat.description = "Niobrara Formation, DJ Basin";
    return mat;
}

MaterialDefinition PermianBasinWolfcamp() {
    auto mat = Materials::Sedimentary::ShaleSiliceous();
    mat.name = "Wolfcamp";
    mat.description = "Wolfcamp Formation, Permian Basin";
    return mat;
}

MaterialDefinition PermianBasinBoneSpring() {
    auto mat = Materials::Sedimentary::Limestone();
    mat.name = "Bone_Spring";
    mat.description = "Bone Spring Formation, Delaware Basin";
    return mat;
}

MaterialDefinition PermianBasinSpraberry() {
    auto mat = Materials::Sedimentary::SandstoneTight();
    mat.name = "Spraberry";
    mat.description = "Spraberry Formation, Midland Basin";
    return mat;
}

MaterialDefinition DelawareBasin() {
    auto mat = Materials::Sedimentary::Limestone();
    mat.name = "Delaware_Basin";
    mat.description = "Delaware Basin carbonates";
    return mat;
}

MaterialDefinition MidlandBasin() {
    auto mat = Materials::Sedimentary::Limestone();
    mat.name = "Midland_Basin";
    mat.description = "Midland Basin carbonates";
    return mat;
}

MaterialDefinition DenverBasin() {
    auto mat = Materials::Sedimentary::Sandstone();
    mat.name = "Denver_Basin";
    mat.description = "Denver Basin clastic rocks";
    return mat;
}

MaterialDefinition WillistonBasin() {
    auto mat = Materials::Sedimentary::Limestone();
    mat.name = "Williston_Basin";
    mat.description = "Williston Basin carbonates";
    return mat;
}

MaterialDefinition AnadarkoBasin() {
    auto mat = Materials::Sedimentary::Sandstone();
    mat.name = "Anadarko_Basin";
    mat.description = "Anadarko Basin mixed clastic/carbonate";
    return mat;
}

MaterialDefinition FortWorthBasin() {
    auto mat = Materials::Sedimentary::Shale();
    mat.name = "Fort_Worth_Basin";
    mat.description = "Fort Worth Basin shales";
    return mat;
}

MaterialDefinition AppalachianBasin() {
    auto mat = Materials::Sedimentary::Shale();
    mat.name = "Appalachian_Basin";
    mat.description = "Appalachian Basin shales";
    return mat;
}

MaterialDefinition GulfCoastMiocene() {
    auto mat = Materials::Sedimentary::Sandstone();
    mat.name = "Gulf_Coast_Miocene";
    mat.description = "Gulf Coast Miocene sands";
    return mat;
}

MaterialDefinition DeepwaterGOM() {
    auto mat = Materials::Sedimentary::Turbidite();
    mat.name = "Deepwater_GOM";
    mat.description = "Deepwater Gulf of Mexico turbidites";
    return mat;
}

MaterialDefinition AlaskanNorthSlope() {
    auto mat = Materials::Sedimentary::Sandstone();
    mat.name = "Alaskan_North_Slope";
    mat.description = "North Slope Alaska reservoirs";
    return mat;
}

MaterialDefinition MontereyFormation() {
    auto mat = Materials::Sedimentary::ShaleSiliceous();
    mat.name = "Monterey_Formation";
    mat.description = "Monterey Formation, California";
    return mat;
}

MaterialDefinition GreenRiverFormation() {
    auto mat = Materials::Sedimentary::OilShale();
    mat.name = "Green_River_Formation";
    mat.description = "Green River oil shale, Colorado/Utah/Wyoming";
    return mat;
}

// International Formations
MaterialDefinition NorthSeaChalk() {
    auto mat = Materials::Sedimentary::LimestoneChalk();
    mat.name = "North_Sea_Chalk";
    mat.description = "North Sea Chalk Group";
    return mat;
}

MaterialDefinition NorthSeaBrentGroup() {
    auto mat = Materials::Sedimentary::Sandstone();
    mat.name = "Brent_Group";
    mat.description = "Brent Group sandstones, North Sea";
    return mat;
}

MaterialDefinition NorthSeaJurassic() {
    auto mat = Materials::Sedimentary::Sandstone();
    mat.name = "North_Sea_Jurassic";
    mat.description = "North Sea Jurassic reservoirs";
    return mat;
}

MaterialDefinition VacasMuertasShale() {
    auto mat = Materials::Sedimentary::ShaleOrganic();
    mat.name = "Vaca_Muerta";
    mat.description = "Vaca Muerta Formation, Argentina";
    return mat;
}

MaterialDefinition SantosBasin() {
    auto mat = Materials::Sedimentary::Limestone();
    mat.name = "Santos_Basin";
    mat.description = "Santos Basin pre-salt carbonates, Brazil";
    return mat;
}

MaterialDefinition CamposBasin() {
    auto mat = Materials::Sedimentary::Limestone();
    mat.name = "Campos_Basin";
    mat.description = "Campos Basin, Brazil";
    return mat;
}

MaterialDefinition OffshoreAngola() {
    auto mat = Materials::Sedimentary::Turbidite();
    mat.name = "Offshore_Angola";
    mat.description = "Offshore Angola turbidites";
    return mat;
}

MaterialDefinition OffshorNigeria() {
    auto mat = Materials::Sedimentary::Turbidite();
    mat.name = "Offshore_Nigeria";
    mat.description = "Niger Delta turbidites";
    return mat;
}

MaterialDefinition GhawarField() {
    auto mat = Materials::Sedimentary::LimestoneOolitic();
    mat.name = "Ghawar_Arab_D";
    mat.description = "Ghawar Field Arab-D reservoir";
    return mat;
}

MaterialDefinition BurganField() {
    auto mat = Materials::Sedimentary::Sandstone();
    mat.name = "Burgan_Field";
    mat.description = "Greater Burgan Field, Kuwait";
    return mat;
}

MaterialDefinition WestSiberiaBasin() {
    auto mat = Materials::Sedimentary::Sandstone();
    mat.name = "West_Siberia";
    mat.description = "West Siberian Basin reservoirs";
    return mat;
}

MaterialDefinition SongliaoBasin() {
    auto mat = Materials::Sedimentary::Sandstone();
    mat.name = "Songliao_Basin";
    mat.description = "Songliao Basin, China";
    return mat;
}

MaterialDefinition SichuanBasin() {
    auto mat = Materials::Sedimentary::Shale();
    mat.name = "Sichuan_Basin";
    mat.description = "Sichuan Basin shales, China";
    return mat;
}

MaterialDefinition CooperBasin() {
    auto mat = Materials::Sedimentary::Sandstone();
    mat.name = "Cooper_Basin";
    mat.description = "Cooper Basin, Australia";
    return mat;
}

MaterialDefinition BrowseBasin() {
    auto mat = Materials::Sedimentary::Sandstone();
    mat.name = "Browse_Basin";
    mat.description = "Browse Basin, offshore Australia";
    return mat;
}

MaterialDefinition CarnarvonBasin() {
    auto mat = Materials::Sedimentary::Sandstone();
    mat.name = "Carnarvon_Basin";
    mat.description = "Carnarvon Basin, offshore Australia";
    return mat;
}

MaterialDefinition PerthBasin() {
    auto mat = Materials::Sedimentary::Sandstone();
    mat.name = "Perth_Basin";
    mat.description = "Perth Basin, Western Australia";
    return mat;
}

MaterialDefinition TaranakiBasin() {
    auto mat = Materials::Sedimentary::Sandstone();
    mat.name = "Taranaki_Basin";
    mat.description = "Taranaki Basin, New Zealand";
    return mat;
}

// Geothermal Formations
MaterialDefinition GeysersGeothermal() {
    auto mat = Materials::Sedimentary::SandstoneGreywacke();
    mat.name = "Geysers_Geothermal";
    mat.description = "The Geysers geothermal field";
    return mat;
}

MaterialDefinition SaltonSeaGeothermal() {
    auto mat = Materials::Sedimentary::Sandstone();
    mat.name = "Salton_Sea_Geothermal";
    mat.description = "Salton Sea geothermal area";
    return mat;
}

MaterialDefinition NewberryGeothermal() {
    auto mat = Materials::Igneous::Basalt();
    mat.name = "Newberry_Geothermal";
    mat.description = "Newberry Volcano EGS site";
    return mat;
}

MaterialDefinition LardereelloGeothermal() {
    auto mat = Materials::Metamorphic::Marble();
    mat.name = "Larderello_Geothermal";
    mat.description = "Larderello geothermal field, Italy";
    return mat;
}

MaterialDefinition WairakeiGeothermal() {
    auto mat = Materials::Igneous::Rhyolite();
    mat.name = "Wairakei_Geothermal";
    mat.description = "Wairakei geothermal field, NZ";
    return mat;
}

MaterialDefinition ReykjanesGeothermal() {
    auto mat = Materials::Igneous::Basalt();
    mat.name = "Reykjanes_Geothermal";
    mat.description = "Reykjanes geothermal field, Iceland";
    return mat;
}

MaterialDefinition SoultzGeothermal() {
    auto mat = Materials::Igneous::Granite();
    mat.name = "Soultz_EGS";
    mat.description = "Soultz-sous-Forêts EGS, France";
    return mat;
}

MaterialDefinition HellisheidiGeothermal() {
    auto mat = Materials::Igneous::Basalt();
    mat.name = "Hellisheidi_Geothermal";
    mat.description = "Hellisheidi geothermal plant, Iceland";
    return mat;
}

MaterialDefinition CerroPrivetoGeothermal() {
    auto mat = Materials::Sedimentary::Sandstone();
    mat.name = "Cerro_Prieto";
    mat.description = "Cerro Prieto geothermal field, Mexico";
    return mat;
}

// Carbon Storage
MaterialDefinition SleipnerAquifer() {
    auto mat = Materials::Sedimentary::Sandstone();
    mat.name = "Sleipner_Utsira";
    mat.description = "Utsira Formation, Sleipner CO2 storage";
    mat.hydraulic.porosity = 0.35;
    mat.hydraulic.setIsotropicPermeability(mDToM2(2000.0));
    return mat;
}

MaterialDefinition InSalahFormation() {
    auto mat = Materials::Sedimentary::Sandstone();
    mat.name = "In_Salah";
    mat.description = "Krechba Formation, In Salah, Algeria";
    return mat;
}

MaterialDefinition QuestFormation() {
    auto mat = Materials::Sedimentary::Sandstone();
    mat.name = "Quest_Basal_Cambrian";
    mat.description = "Basal Cambrian Sands, Quest CCS, Canada";
    return mat;
}

MaterialDefinition DecaturFormation() {
    auto mat = Materials::Sedimentary::Sandstone();
    mat.name = "Decatur_Mt_Simon";
    mat.description = "Mt. Simon Sandstone, Illinois Basin CCS";
    return mat;
}

MaterialDefinition GorgoFormation() {
    auto mat = Materials::Sedimentary::Sandstone();
    mat.name = "Gorgon_Dupuy";
    mat.description = "Dupuy Formation, Gorgon CO2 storage, Australia";
    return mat;
}

MaterialDefinition UtsiraFormation() {
    auto mat = Materials::Sedimentary::Sandstone();
    mat.name = "Utsira_Formation";
    mat.description = "Utsira Sand, North Sea";
    mat.hydraulic.porosity = 0.35;
    mat.hydraulic.setIsotropicPermeability(mDToM2(2000.0));
    return mat;
}

MaterialDefinition MtSimonSandstone() {
    auto mat = Materials::Sedimentary::SandstoneQuartz();
    mat.name = "Mt_Simon_Sandstone";
    mat.description = "Mt. Simon Sandstone, Illinois Basin";
    return mat;
}

MaterialDefinition SaukSequence() {
    auto mat = Materials::Sedimentary::Sandstone();
    mat.name = "Sauk_Sequence";
    mat.description = "Sauk Megasequence basal sandstones";
    return mat;
}

// Mining/Ore Bodies
MaterialDefinition WitwatersrandGold() {
    auto mat = Materials::Sedimentary::Conglomerate();
    mat.name = "Witwatersrand_Gold";
    mat.description = "Witwatersrand gold-bearing conglomerates";
    return mat;
}

MaterialDefinition BushveldComplex() {
    auto mat = Materials::Igneous::Norite();
    mat.name = "Bushveld_Complex";
    mat.description = "Bushveld Igneous Complex, South Africa";
    return mat;
}

MaterialDefinition KirunaIronOre() {
    auto mat = Materials::Minerals::Magnetite();
    mat.name = "Kiruna_Iron_Ore";
    mat.description = "Kiruna magnetite deposits, Sweden";
    return mat;
}

MaterialDefinition OlympicDamCopper() {
    auto mat = Materials::Igneous::Granite();
    mat.name = "Olympic_Dam";
    mat.description = "Olympic Dam IOCG deposit, Australia";
    return mat;
}

MaterialDefinition GrasbergCopper() {
    auto mat = Materials::Igneous::Diorite();
    mat.name = "Grasberg_Copper";
    mat.description = "Grasberg porphyry copper, Indonesia";
    return mat;
}

MaterialDefinition ChuquicamataCopper() {
    auto mat = Materials::Igneous::Porphyry();
    mat.name = "Chuquicamata_Copper";
    mat.description = "Chuquicamata porphyry copper, Chile";
    return mat;
}

// Aquifer Formations
MaterialDefinition OgallalaAquifer() {
    auto mat = Materials::Unconsolidated::GravelSandy();
    mat.name = "Ogallala_Aquifer";
    mat.description = "High Plains Aquifer";
    mat.hydraulic.porosity = 0.30;
    mat.hydraulic.setIsotropicPermeability(mDToM2(5000.0));
    return mat;
}

MaterialDefinition EdwardsAquifer() {
    auto mat = Materials::Sedimentary::Limestone();
    mat.name = "Edwards_Aquifer";
    mat.description = "Edwards Aquifer, Texas";
    mat.hydraulic.porosity = 0.20;
    return mat;
}

MaterialDefinition FlordianAquifer() {
    auto mat = Materials::Sedimentary::Limestone();
    mat.name = "Floridan_Aquifer";
    mat.description = "Floridan Aquifer System";
    mat.hydraulic.porosity = 0.25;
    return mat;
}

MaterialDefinition HighPlainsAquifer() {
    auto mat = Materials::Unconsolidated::Sand();
    mat.name = "High_Plains_Aquifer";
    mat.description = "High Plains Aquifer (Ogallala)";
    return mat;
}

MaterialDefinition CentralValleyAquifer() {
    auto mat = Materials::Unconsolidated::Alluvium();
    mat.name = "Central_Valley_Aquifer";
    mat.description = "California Central Valley aquifer";
    return mat;
}

MaterialDefinition GreatArtesianBasin() {
    auto mat = Materials::Sedimentary::Sandstone();
    mat.name = "Great_Artesian_Basin";
    mat.description = "Great Artesian Basin, Australia";
    return mat;
}

MaterialDefinition NubianSandstoneAquifer() {
    auto mat = Materials::Sedimentary::Sandstone();
    mat.name = "Nubian_Sandstone_Aquifer";
    mat.description = "Nubian Sandstone Aquifer System, North Africa";
    return mat;
}

MaterialDefinition GuaraniAquifer() {
    auto mat = Materials::Sedimentary::Sandstone();
    mat.name = "Guarani_Aquifer";
    mat.description = "Guarani Aquifer, South America";
    return mat;
}

// Seismically Active Regions
MaterialDefinition SanAndreassFault() {
    auto mat = Materials::FaultMaterials::FaultGouge();
    mat.name = "San_Andreas_Fault";
    mat.description = "San Andreas Fault zone material";
    return mat;
}

MaterialDefinition HaywardFault() {
    auto mat = Materials::FaultMaterials::FaultGouge();
    mat.name = "Hayward_Fault";
    mat.description = "Hayward Fault zone material";
    return mat;
}

MaterialDefinition NewMadridZone() {
    auto mat = Materials::Unconsolidated::Alluvium();
    mat.name = "New_Madrid_Zone";
    mat.description = "New Madrid Seismic Zone sediments";
    return mat;
}

MaterialDefinition WasatchFault() {
    auto mat = Materials::FaultMaterials::FaultBreccia();
    mat.name = "Wasatch_Fault";
    mat.description = "Wasatch Fault zone material";
    return mat;
}

MaterialDefinition AlpineFaultNZ() {
    auto mat = Materials::FaultMaterials::FaultMylonite();
    mat.name = "Alpine_Fault_NZ";
    mat.description = "Alpine Fault zone, New Zealand";
    return mat;
}

MaterialDefinition NAFZone() {
    auto mat = Materials::FaultMaterials::FaultGouge();
    mat.name = "North_Anatolian_Fault";
    mat.description = "North Anatolian Fault Zone, Turkey";
    return mat;
}

MaterialDefinition JapanTrench() {
    auto mat = Materials::Sedimentary::Turbidite();
    mat.name = "Japan_Trench";
    mat.description = "Japan Trench subduction zone sediments";
    return mat;
}

MaterialDefinition CascadiaSubduction() {
    auto mat = Materials::Sedimentary::Turbidite();
    mat.name = "Cascadia_Subduction";
    mat.description = "Cascadia Subduction Zone sediments";
    return mat;
}

} // namespace Formations

} // namespace FSRM
