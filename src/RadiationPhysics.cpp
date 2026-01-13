/**
 * @file RadiationPhysics.cpp
 * @brief Implementation of radiation physics for nuclear explosions
 * 
 * See RadiationPhysics.hpp for detailed documentation.
 */

#include "RadiationPhysics.hpp"
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <random>

namespace FSRM {

using namespace RadiationConstants;

// =============================================================================
// RadionuclideDatabase Implementation
// =============================================================================

RadionuclideDatabase::RadionuclideDatabase() {
    loadDefaultFissionProducts();
}

void RadionuclideDatabase::loadDefaultFissionProducts() {
    // Add key fission products with their properties
    // Data based on England & Rider (1994) and ICRP databases
    
    // Format: name, Z, A, half_life (s), fission_yield, category
    
    // Noble gases
    addDefaultNuclide("Kr-85", 36, 85, 3.394e8, 0.00286, FissionProductCategory::NOBLE_GAS);
    addDefaultNuclide("Kr-88", 36, 88, 1.022e4, 0.0357, FissionProductCategory::NOBLE_GAS);
    addDefaultNuclide("Xe-133", 54, 133, 4.528e5, 0.0671, FissionProductCategory::NOBLE_GAS);
    addDefaultNuclide("Xe-135", 54, 135, 3.291e4, 0.0661, FissionProductCategory::NOBLE_GAS);
    
    // Halogens
    addDefaultNuclide("I-131", 53, 131, 6.949e5, 0.0290, FissionProductCategory::HALOGEN);
    addDefaultNuclide("I-132", 53, 132, 8.28e3, 0.0427, FissionProductCategory::HALOGEN);
    addDefaultNuclide("I-133", 53, 133, 7.49e4, 0.0671, FissionProductCategory::HALOGEN);
    addDefaultNuclide("I-135", 53, 135, 2.37e4, 0.0632, FissionProductCategory::HALOGEN);
    addDefaultNuclide("Br-83", 35, 83, 8.64e3, 0.00538, FissionProductCategory::HALOGEN);
    
    // Volatile elements
    addDefaultNuclide("Cs-134", 55, 134, 6.507e7, 0.0, FissionProductCategory::VOLATILE);  // Activation product
    addDefaultNuclide("Cs-137", 55, 137, 9.467e8, 0.0627, FissionProductCategory::VOLATILE);
    addDefaultNuclide("Rb-86", 37, 86, 1.612e6, 0.00194, FissionProductCategory::VOLATILE);
    addDefaultNuclide("Te-132", 52, 132, 2.77e5, 0.0427, FissionProductCategory::VOLATILE);
    
    // Refractory elements
    addDefaultNuclide("Sr-89", 38, 89, 4.363e6, 0.0474, FissionProductCategory::REFRACTORY);
    addDefaultNuclide("Sr-90", 38, 90, 9.08e8, 0.0577, FissionProductCategory::REFRACTORY);
    addDefaultNuclide("Y-90", 39, 90, 2.304e5, 0.0, FissionProductCategory::REFRACTORY);  // From Sr-90 decay
    addDefaultNuclide("Y-91", 39, 91, 5.06e6, 0.0584, FissionProductCategory::REFRACTORY);
    addDefaultNuclide("Zr-95", 40, 95, 5.53e6, 0.0650, FissionProductCategory::REFRACTORY);
    addDefaultNuclide("Nb-95", 41, 95, 3.02e6, 0.0, FissionProductCategory::REFRACTORY);  // From Zr-95
    addDefaultNuclide("Mo-99", 42, 99, 2.375e5, 0.0614, FissionProductCategory::REFRACTORY);
    addDefaultNuclide("Tc-99m", 43, 99, 2.17e4, 0.0, FissionProductCategory::REFRACTORY);
    addDefaultNuclide("Ru-103", 44, 103, 3.39e6, 0.0305, FissionProductCategory::REFRACTORY);
    addDefaultNuclide("Ru-106", 44, 106, 3.22e7, 0.00401, FissionProductCategory::REFRACTORY);
    addDefaultNuclide("Ba-140", 56, 140, 1.101e6, 0.0625, FissionProductCategory::REFRACTORY);
    addDefaultNuclide("La-140", 57, 140, 1.45e5, 0.0, FissionProductCategory::REFRACTORY);
    addDefaultNuclide("Ce-141", 58, 141, 2.81e6, 0.0585, FissionProductCategory::REFRACTORY);
    addDefaultNuclide("Ce-144", 58, 144, 2.46e7, 0.0546, FissionProductCategory::REFRACTORY);
    addDefaultNuclide("Pr-143", 59, 143, 1.17e6, 0.0597, FissionProductCategory::REFRACTORY);
    addDefaultNuclide("Nd-147", 60, 147, 9.50e5, 0.0225, FissionProductCategory::REFRACTORY);
    addDefaultNuclide("Pm-147", 61, 147, 8.27e7, 0.0, FissionProductCategory::REFRACTORY);
    
    // Actinides (for completeness)
    addDefaultNuclide("U-237", 92, 237, 5.83e5, 0.0, FissionProductCategory::ACTINIDE);
    addDefaultNuclide("Np-239", 93, 239, 2.035e5, 0.0, FissionProductCategory::ACTINIDE);
    addDefaultNuclide("Pu-239", 94, 239, 7.59e11, 0.0, FissionProductCategory::ACTINIDE);
    
    // Set up decay chains
    decay_chains_["Kr-88"] = {{"Rb-88", 1.0}};
    decay_chains_["Rb-88"] = {{"Sr-88", 1.0}};
    decay_chains_["Xe-133"] = {{"Cs-133", 1.0}};
    decay_chains_["Xe-135"] = {{"Cs-135", 1.0}};
    decay_chains_["I-131"] = {{"Xe-131", 1.0}};
    decay_chains_["I-132"] = {{"Xe-132", 1.0}};
    decay_chains_["I-133"] = {{"Xe-133", 1.0}};
    decay_chains_["I-135"] = {{"Xe-135", 1.0}};
    decay_chains_["Te-132"] = {{"I-132", 1.0}};
    decay_chains_["Sr-90"] = {{"Y-90", 1.0}};
    decay_chains_["Y-90"] = {{"Zr-90", 1.0}};
    decay_chains_["Zr-95"] = {{"Nb-95", 0.99}};
    decay_chains_["Mo-99"] = {{"Tc-99m", 0.876}, {"Tc-99", 0.124}};
    decay_chains_["Ba-140"] = {{"La-140", 1.0}};
    decay_chains_["Ce-144"] = {{"Pr-144", 1.0}};
}

void RadionuclideDatabase::addDefaultNuclide(const std::string& name, int Z, int A, 
                                             double half_life, double fission_yield,
                                             FissionProductCategory cat) {
    Radionuclide nuclide;
    nuclide.name = name;
    nuclide.atomic_number = Z;
    nuclide.mass_number = A;
    nuclide.half_life = half_life;
    nuclide.fission_yield_thermal = fission_yield;
    nuclide.fission_yield_fast = fission_yield * 0.95;  // Approximate
    nuclide.category = cat;
    
    // Extract symbol
    size_t dash = name.find('-');
    nuclide.symbol = name.substr(0, dash);
    
    // Set dose coefficients (simplified - real values from ICRP)
    switch (cat) {
        case FissionProductCategory::HALOGEN:
            nuclide.dose_coefficient_inh = 2.0e-8;  // Sv/Bq (I-131 like)
            nuclide.dose_coefficient_ing = 2.2e-8;
            nuclide.dose_rate_factor = 6.0e-13;     // (Sv/h)/(Bq/m²)
            break;
        case FissionProductCategory::VOLATILE:
            nuclide.dose_coefficient_inh = 3.9e-8;  // Sv/Bq (Cs-137 like)
            nuclide.dose_coefficient_ing = 1.3e-8;
            nuclide.dose_rate_factor = 1.1e-12;
            break;
        case FissionProductCategory::NOBLE_GAS:
            nuclide.dose_coefficient_inh = 0.0;     // Noble gases not inhaled
            nuclide.dose_coefficient_ing = 0.0;
            nuclide.dose_rate_factor = 4.0e-14;    // Submersion dose
            break;
        default:
            nuclide.dose_coefficient_inh = 1.0e-8;
            nuclide.dose_coefficient_ing = 1.0e-8;
            nuclide.dose_rate_factor = 5.0e-13;
    }
    
    // Add gamma lines for major emitters
    if (name == "Cs-137") {
        nuclide.gamma_lines.push_back({0.662, 0.85, 0.09});  // 662 keV, 85% intensity
        nuclide.beta_endpoint_energy = 0.512;
        nuclide.beta_average_energy = 0.188;
    } else if (name == "I-131") {
        nuclide.gamma_lines.push_back({0.364, 0.81, 0.0});
        nuclide.gamma_lines.push_back({0.637, 0.072, 0.0});
        nuclide.beta_endpoint_energy = 0.606;
        nuclide.beta_average_energy = 0.192;
    } else if (name == "Sr-90") {
        nuclide.beta_endpoint_energy = 0.546;
        nuclide.beta_average_energy = 0.196;
    } else if (name == "Ba-140") {
        nuclide.gamma_lines.push_back({0.537, 0.24, 0.0});
        nuclide.beta_endpoint_energy = 1.02;
        nuclide.beta_average_energy = 0.365;
    }
    
    database_[name] = nuclide;
}

const Radionuclide* RadionuclideDatabase::get(const std::string& name) const {
    auto it = database_.find(name);
    if (it != database_.end()) {
        return &it->second;
    }
    return nullptr;
}

std::vector<std::string> RadionuclideDatabase::getFissionProducts() const {
    std::vector<std::string> products;
    for (const auto& pair : database_) {
        if (pair.second.fission_yield_thermal > 0) {
            products.push_back(pair.first);
        }
    }
    return products;
}

std::vector<std::string> RadionuclideDatabase::getByCategory(FissionProductCategory cat) const {
    std::vector<std::string> result;
    for (const auto& pair : database_) {
        if (pair.second.category == cat) {
            result.push_back(pair.first);
        }
    }
    return result;
}

std::vector<std::string> RadionuclideDatabase::getDecayChain(const std::string& parent) const {
    std::vector<std::string> chain;
    chain.push_back(parent);
    
    auto it = decay_chains_.find(parent);
    while (it != decay_chains_.end() && !it->second.empty()) {
        std::string daughter = it->second[0].first;  // First daughter
        chain.push_back(daughter);
        it = decay_chains_.find(daughter);
    }
    
    return chain;
}

// =============================================================================
// PromptRadiationModel Implementation
// =============================================================================

PromptRadiationModel::PromptRadiationModel()
    : yield_kt_(100.0), fission_fraction_(0.5), burst_height_(0.0) {
    useDefaultSpectra();
}

void PromptRadiationModel::setYield(double yield_kt) {
    yield_kt_ = yield_kt;
}

void PromptRadiationModel::setFissionFraction(double fraction) {
    fission_fraction_ = fraction;
}

void PromptRadiationModel::setBurstHeight(double height_m) {
    burst_height_ = height_m;
}

void PromptRadiationModel::setFusionEnhancement(double factor) {
    fusion_enhancement_ = factor;
}

void PromptRadiationModel::useDefaultSpectra() {
    // Fission gamma spectrum (simplified exponential)
    gamma_spectrum_ = [](double E) {
        if (E < 0.1 || E > 10.0) return 0.0;
        // Peak around 0.7 MeV with exponential tails
        return std::exp(-std::pow(E - 0.7, 2) / 0.5) + 0.3 * std::exp(-E / 2.0);
    };
    
    // Fission neutron spectrum (Watt)
    neutron_spectrum_ = [this](double E) {
        return wattSpectrum(E);
    };
}

double PromptRadiationModel::wattSpectrum(double E, double a, double b) const {
    // N(E) = C * exp(-E/a) * sinh(sqrt(b*E))
    // For U-235 thermal fission: a ≈ 0.988 MeV, b ≈ 2.249 MeV^-1
    if (E <= 0) return 0.0;
    
    double C = 0.4527;  // Normalization constant
    return C * std::exp(-E / a) * std::sinh(std::sqrt(b * E));
}

double PromptRadiationModel::totalGammaEnergy() const {
    // About 3.5% of fission energy appears as prompt gamma
    double fission_energy_J = fission_fraction_ * yield_kt_ * JOULES_PER_KT;
    return 0.035 * fission_energy_J / JOULES_PER_MEV;  // MeV
}

double PromptRadiationModel::totalNeutronEnergy() const {
    // About 3.0% of fission energy appears as prompt neutrons
    double fission_energy_J = fission_fraction_ * yield_kt_ * JOULES_PER_KT;
    return 0.030 * fission_energy_J / JOULES_PER_MEV;  // MeV
}

double PromptRadiationModel::gammaYieldFraction() const {
    return 0.035 * fission_fraction_;
}

double PromptRadiationModel::neutronYieldFraction() const {
    return 0.030 * fission_fraction_ * fusion_enhancement_;
}

double PromptRadiationModel::gammaSpectrum(double energy_MeV) const {
    return gamma_spectrum_(energy_MeV);
}

double PromptRadiationModel::neutronSpectrum(double energy_MeV) const {
    return neutron_spectrum_(energy_MeV);
}

double PromptRadiationModel::gammaAirAttenuation(double energy_MeV, double distance_m) const {
    // Linear attenuation coefficient in air
    double mu = gammaAttenuationCoeff(energy_MeV);
    
    // Include buildup factor for thick shields
    double mfp = 1.0 / mu;
    double x = distance_m;
    double B = gammaBuildupFactor(energy_MeV, x / mfp);
    
    return B * std::exp(-mu * x);
}

double PromptRadiationModel::gammaAttenuationCoeff(double energy_MeV) const {
    // Mass attenuation coefficient for air (cm²/g)
    // Simplified fit to NIST data
    double mu_mass;
    if (energy_MeV < 0.1) {
        mu_mass = 0.15 * std::pow(energy_MeV / 0.1, -2.5);
    } else if (energy_MeV < 3.0) {
        mu_mass = 0.08 * std::pow(energy_MeV, -0.3);
    } else {
        mu_mass = 0.02;
    }
    
    // Convert to linear attenuation coefficient (1/m)
    double rho_air = AIR_DENSITY_STP;  // kg/m³
    return mu_mass * rho_air * 100.0;  // 1/m (factor 100 for cm² to m²)
}

double PromptRadiationModel::neutronAttenuationCoeff(double energy_MeV) const {
    // Approximate neutron removal cross-section in air
    double sigma_r;  // cm²/g
    
    if (energy_MeV < 0.5) {
        sigma_r = 0.05;
    } else if (energy_MeV < 5.0) {
        sigma_r = 0.04;
    } else {
        sigma_r = 0.03;
    }
    
    double rho_air = AIR_DENSITY_STP;
    return sigma_r * rho_air * 100.0;  // 1/m
}

double PromptRadiationModel::gammaBuildupFactor(double energy_MeV, double mfp) const {
    // Taylor buildup factor approximation
    // B(E, x) = 1 + α·μx for small μx
    double alpha = 1.0;  // Simplified
    if (energy_MeV < 0.5) alpha = 0.5;
    if (energy_MeV > 2.0) alpha = 1.5;
    
    if (mfp < 5.0) {
        return 1.0 + alpha * mfp;
    } else {
        // Saturates for thick shields
        return 1.0 + alpha * 5.0 + (mfp - 5.0) * 0.1;
    }
}

double PromptRadiationModel::neutronAirAttenuation(double energy_MeV, double distance_m) const {
    double mu = neutronAttenuationCoeff(energy_MeV);
    return std::exp(-mu * distance_m);
}

double PromptRadiationModel::gammaFluence(double r) const {
    // Total gamma fluence at distance r
    // Φ = N_γ / (4πr²) * exp(-μr) * B
    
    double N_gamma = totalGammaEnergy() / 1.0;  // Photons (assume 1 MeV average)
    
    // Account for height of burst
    double slant_range = std::sqrt(r * r + burst_height_ * burst_height_);
    
    // Average attenuation at 1 MeV
    double atten = gammaAirAttenuation(1.0, slant_range);
    
    return N_gamma / (4.0 * M_PI * slant_range * slant_range) * atten;
}

double PromptRadiationModel::neutronFluence(double r) const {
    // Total neutron fluence at distance r
    double mean_energy = 2.0;  // MeV average for fission neutrons
    double N_neutron = fission_fraction_ * yield_kt_ * FISSIONS_PER_KT * 2.5;  // ~2.5 neutrons/fission
    
    double slant_range = std::sqrt(r * r + burst_height_ * burst_height_);
    double atten = neutronAirAttenuation(mean_energy, slant_range);
    
    return N_neutron / (4.0 * M_PI * slant_range * slant_range) * atten * fusion_enhancement_;
}

double PromptRadiationModel::gammaDose(double r) const {
    // Convert fluence to dose
    // D = Φ * (μ_en/ρ) * E / density_tissue
    
    double fluence = gammaFluence(r);
    double mean_energy_MeV = 1.0;
    double mu_en_rho = 0.03;  // cm²/g for tissue at 1 MeV
    
    // Dose in Gy
    double dose_Gy = fluence * mu_en_rho * mean_energy_MeV * JOULES_PER_MEV * 1e4 / 1.0;  // 1 g/cm³ tissue
    
    return dose_Gy;
}

double PromptRadiationModel::neutronDose(double r) const {
    // Neutron dose with RBE
    double fluence = neutronFluence(r);
    double mean_energy_MeV = 2.0;
    
    // Kerma factor for tissue (approximate)
    double kerma_factor = 3.0e-11;  // Gy·cm²/neutron at 2 MeV
    
    // Dose
    double dose_Gy = fluence * kerma_factor * 1e4;  // m² to cm²
    
    // RBE is already somewhat included in kerma factor, but can add modifier
    return dose_Gy;
}

double PromptRadiationModel::totalPromptDose(double r) const {
    return gammaDose(r) + neutronDose(r);
}

double PromptRadiationModel::doseRate(double r, double t) const {
    // Most prompt radiation arrives within microseconds
    // Model as exponential pulse
    double tau = 1e-6;  // 1 μs characteristic time
    
    if (t < 0) return 0.0;
    if (t > 10.0 * tau) return 0.0;
    
    double total_dose = totalPromptDose(r);
    return total_dose / tau * std::exp(-t / tau);
}

double PromptRadiationModel::shieldedGammaDose(double r, double shield_thickness_m,
                                               double shield_density, double /*shield_Z*/) const {
    // Tenth-value layer for gamma at 1 MeV
    // TVL ≈ 30 cm for concrete
    double mu_shield = shield_density / 2400.0 * 0.077;  // Scale from concrete
    double shield_factor = std::exp(-mu_shield * shield_thickness_m * 100.0);
    
    return gammaDose(r) * shield_factor;
}

double PromptRadiationModel::shieldedNeutronDose(double r, double shield_thickness_m,
                                                 double /*hydrogen_density*/) const {
    // Neutron shielding - hydrogen-rich materials
    // TVL ≈ 10-15 cm for polyethylene
    double mu_shield = 0.1;  // 1/cm approximate
    double shield_factor = std::exp(-mu_shield * shield_thickness_m * 100.0);
    
    return neutronDose(r) * shield_factor;
}

// =============================================================================
// FissionProductInventory Implementation
// =============================================================================

FissionProductInventory::FissionProductInventory() 
    : fission_yield_kt_(100.0), current_time_(0.0) {
    database_ = std::make_shared<RadionuclideDatabase>();
}

void FissionProductInventory::setDatabase(std::shared_ptr<RadionuclideDatabase> db) {
    database_ = db;
}

void FissionProductInventory::setFissileType(const std::string& type) {
    fissile_type_ = type;
}

void FissionProductInventory::setFissionYield(double yield_kt) {
    fission_yield_kt_ = yield_kt;
}

void FissionProductInventory::initialize() {
    activities_.clear();
    initial_atoms_.clear();
    
    // Total number of fissions
    double total_fissions = fission_yield_kt_ * FISSIONS_PER_KT;
    
    // Initialize each fission product
    auto products = database_->getFissionProducts();
    for (const auto& name : products) {
        const Radionuclide* nuc = database_->get(name);
        if (!nuc) continue;
        
        // Initial atoms from fission yield
        double N0 = total_fissions * nuc->fission_yield_thermal;
        initial_atoms_[name] = N0;
        
        // Initial activity (Bq)
        double lambda = nuc->decay_constant();
        activities_[name] = lambda * N0;
    }
    
    current_time_ = 0.0;
}

void FissionProductInventory::evolve(double dt) {
    evolveTo(current_time_ + dt);
}

void FissionProductInventory::evolveTo(double time) {
    if (time <= current_time_) return;
    
    double dt = time - current_time_;
    
    // Simple decay for now (Bateman for chains would be more accurate)
    for (auto& pair : activities_) {
        const Radionuclide* nuc = database_->get(pair.first);
        if (!nuc) continue;
        
        double lambda = nuc->decay_constant();
        pair.second *= std::exp(-lambda * dt);
    }
    
    current_time_ = time;
}

double FissionProductInventory::totalActivity() const {
    double total = 0.0;
    for (const auto& pair : activities_) {
        total += pair.second;
    }
    return total;
}

double FissionProductInventory::activity(const std::string& nuclide) const {
    auto it = activities_.find(nuclide);
    if (it != activities_.end()) {
        return it->second;
    }
    return 0.0;
}

double FissionProductInventory::totalActivity(double t) const {
    // Activity at arbitrary time using Way-Wigner
    return wayWignerActivity(t);
}

double FissionProductInventory::wayWignerActivity(double t) const {
    // Way-Wigner approximation: A(t) = 1.4 × 10²³ × P × t^(-1.2) Bq
    // where P is fission power in MW and t is time in seconds after fission
    // For weapons: use total fissions
    
    if (t < 0.1) t = 0.1;  // Avoid singularity
    
    double total_fissions = fission_yield_kt_ * FISSIONS_PER_KT;
    
    // At 1 second: A ≈ 3.7 × 10²² × kt fissions
    double A_1s = 3.7e22 * fission_yield_kt_;
    
    // t^-1.2 decay
    return A_1s * std::pow(t, -1.2);
}

double FissionProductInventory::wayWignerPower(double t) const {
    // Decay heat power (W)
    // P = 0.066 × P_0 × [(t-t_s)^(-0.2) - (t-t_0)^(-0.2)]
    // Simplified: P ≈ 0.066 × P_0 × t^(-0.2)
    
    if (t < 0.1) t = 0.1;
    
    double P_0 = fission_yield_kt_ * JOULES_PER_KT;  // Total energy
    return 0.066 * P_0 * std::pow(t, -0.2);
}

double FissionProductInventory::iodine131Activity() const {
    return activity("I-131");
}

double FissionProductInventory::cesium137Activity() const {
    return activity("Cs-137");
}

double FissionProductInventory::strontium90Activity() const {
    return activity("Sr-90");
}

double FissionProductInventory::xenon133Activity() const {
    return activity("Xe-133");
}

std::vector<std::pair<std::string, double>> 
FissionProductInventory::topActivityContributors(int n) const {
    std::vector<std::pair<std::string, double>> result;
    for (const auto& pair : activities_) {
        result.push_back(pair);
    }
    
    std::sort(result.begin(), result.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });
    
    if (static_cast<int>(result.size()) > n) {
        result.resize(n);
    }
    return result;
}

void FissionProductInventory::writeInventory(const std::string& filename) const {
    std::ofstream out(filename);
    
    out << "# Fission Product Inventory\n";
    out << "# Time: " << current_time_ << " s\n";
    out << "# Total Activity: " << totalActivity() << " Bq\n";
    out << "#\n";
    out << "# Nuclide  Half-life(s)  Activity(Bq)  Fraction\n";
    
    double total = totalActivity();
    
    for (const auto& pair : activities_) {
        const Radionuclide* nuc = database_->get(pair.first);
        if (!nuc) continue;
        
        out << std::setw(10) << pair.first
            << std::setw(14) << std::scientific << nuc->half_life
            << std::setw(14) << pair.second
            << std::setw(12) << std::fixed << pair.second / total
            << "\n";
    }
}

// =============================================================================
// FalloutFormationModel Implementation
// =============================================================================

FalloutFormationModel::FalloutFormationModel()
    : yield_kt_(100.0), fission_yield_kt_(50.0), burst_type_("surface") {
    rng_.seed(std::random_device{}());
}

void FalloutFormationModel::setYield(double yield_kt) {
    yield_kt_ = yield_kt;
}

void FalloutFormationModel::setFissionYield(double fission_kt) {
    fission_yield_kt_ = fission_kt;
}

void FalloutFormationModel::setBurstType(const std::string& type) {
    burst_type_ = type;
}

void FalloutFormationModel::setGroundMaterial(double density, double melt_temp, double vaporize_temp) {
    soil_density_ = density;
    soil_melt_temp_ = melt_temp;
    soil_vaporize_temp_ = vaporize_temp;
}

void FalloutFormationModel::setInventory(std::shared_ptr<FissionProductInventory> inventory) {
    inventory_ = inventory;
}

void FalloutFormationModel::computeParticleFormation() {
    // Compute masses
    mass_in_cloud_ = fireballMass();
    mass_in_stem_ = stemEntrainmentMass();
    
    if (burst_type_ == "surface") {
        total_mass_ = craterMass() + mass_in_cloud_ + mass_in_stem_;
    } else {
        total_mass_ = mass_in_cloud_;
    }
}

double FalloutFormationModel::particleSizePDF(double diameter) const {
    // Log-normal distribution
    double ln_d = std::log(diameter);
    double ln_median = std::log(median_diameter_);
    double ln_sigma = std::log(geometric_std_);
    
    double exponent = -0.5 * std::pow((ln_d - ln_median) / ln_sigma, 2);
    return 1.0 / (diameter * ln_sigma * std::sqrt(2.0 * M_PI)) * std::exp(exponent);
}

double FalloutFormationModel::particleSizeCDF(double diameter) const {
    double ln_d = std::log(diameter);
    double ln_median = std::log(median_diameter_);
    double ln_sigma = std::log(geometric_std_);
    
    return 0.5 * (1.0 + std::erf((ln_d - ln_median) / (ln_sigma * std::sqrt(2.0))));
}

void FalloutFormationModel::setLogNormalDistribution(double median_diameter, double geometric_std) {
    median_diameter_ = median_diameter;
    geometric_std_ = geometric_std;
}

double FalloutFormationModel::volatileFractionInSmallParticles() const {
    // Volatile species condense preferentially on small particles
    // Fractionation factor increases for smaller particles
    return 0.3;  // 30% of volatiles in smallest size class
}

double FalloutFormationModel::refractoryFractionInLargeParticles() const {
    // Refractory species are volume-distributed in large particles
    return 0.7;  // 70% in large particles
}

double FalloutFormationModel::fractionationFactor(const std::string& nuclide, double diameter) const {
    const Radionuclide* nuc = nullptr;
    if (inventory_ && inventory_->getDatabase()) {
        nuc = inventory_->getDatabase()->get(nuclide);
    }
    if (!nuc) return 1.0;
    
    // Fractionation depends on volatility
    double d_ref = median_diameter_;
    
    switch (nuc->category) {
        case FissionProductCategory::NOBLE_GAS:
            return 0.0;  // Not in particles
            
        case FissionProductCategory::HALOGEN:
        case FissionProductCategory::VOLATILE:
            // More in small particles (surface condensation)
            return std::pow(d_ref / diameter, 0.5);
            
        case FissionProductCategory::REFRACTORY:
            // Volume distributed
            return 1.0;
            
        default:
            return 1.0;
    }
}

double FalloutFormationModel::totalFalloutMass() const {
    return total_mass_;
}

double FalloutFormationModel::falloutMassInCloud() const {
    return mass_in_cloud_;
}

double FalloutFormationModel::falloutMassInStem() const {
    return mass_in_stem_;
}

double FalloutFormationModel::totalFalloutActivity() const {
    if (!inventory_) return 0.0;
    
    // Activity attached to particles (excluding noble gases)
    double total = inventory_->totalActivity();
    double noble_gas_fraction = 0.15;  // ~15% of activity in noble gases
    
    return total * (1.0 - noble_gas_fraction);
}

std::vector<FalloutFormationModel::FalloutParticle> 
FalloutFormationModel::generateParticles(int n_particles) const {
    std::vector<FalloutParticle> particles;
    particles.reserve(n_particles);
    
    std::lognormal_distribution<double> size_dist(std::log(median_diameter_), 
                                                  std::log(geometric_std_));
    
    double total_activity = totalFalloutActivity();
    double mass_per_particle = total_mass_ / n_particles;
    
    for (int i = 0; i < n_particles; ++i) {
        FalloutParticle p;
        
        // Sample diameter from log-normal
        p.diameter = size_dist(rng_);
        p.diameter = std::max(1e-7, std::min(p.diameter, 1e-2));  // Limit range
        
        // Density depends on size (small = condensate, large = soil)
        if (p.diameter < 10e-6) {
            p.density = 2000.0;  // Condensate
        } else {
            p.density = soil_density_;
        }
        
        p.mass = mass_per_particle;
        
        // Activity proportional to surface area for volatiles, volume for refractories
        double d_ref = median_diameter_;
        double activity_factor = 0.5 * std::pow(p.diameter / d_ref, 2) + 
                                0.5 * std::pow(p.diameter / d_ref, 3);
        activity_factor /= (0.5 + 0.5);  // Normalize
        
        p.activity = total_activity / n_particles * activity_factor;
        p.specific_activity = p.activity / p.mass;
        
        particles.push_back(p);
    }
    
    return particles;
}

double FalloutFormationModel::craterMass() const {
    // Mass from crater ejecta (surface burst)
    if (burst_type_ != "surface") return 0.0;
    
    // Crater radius ~ 30 * W^0.3 meters
    double R_crater = 30.0 * std::pow(yield_kt_, 0.3);
    double depth = R_crater * 0.3;  // Depth/radius ratio
    
    // Parabolic crater volume
    double V = 0.5 * M_PI * R_crater * R_crater * depth;
    
    return V * soil_density_ * 0.1;  // ~10% becomes fine fallout
}

double FalloutFormationModel::fireballMass() const {
    // Mass of air in fireball
    double R_fb = 66.0 * std::pow(yield_kt_, 0.4);  // Fireball radius
    double V = 4.0/3.0 * M_PI * R_fb * R_fb * R_fb;
    
    return V * AIR_DENSITY_STP * 0.1;  // ~10% becomes particulate
}

double FalloutFormationModel::stemEntrainmentMass() const {
    // Mass entrained into stem (surface burst only)
    if (burst_type_ != "surface") return 0.0;
    
    double R_stem = 0.3 * 66.0 * std::pow(yield_kt_, 0.4);
    double H_stem = 5000.0;  // 5 km stem height
    
    double V = M_PI * R_stem * R_stem * H_stem;
    return V * AIR_DENSITY_STP * 0.01 + soil_density_ * R_stem * R_stem * 100.0;
}

double FalloutFormationModel::FalloutParticle::settlingVelocity(
    double air_density, double air_viscosity) const {
    // Stokes settling velocity with Cunningham correction for small particles
    double dp = diameter;
    
    // Mean free path of air
    double lambda = 6.6e-8;  // m at STP
    
    // Cunningham correction factor
    double Cc = 1.0 + lambda/dp * (2.514 + 0.8 * std::exp(-0.55 * dp / lambda));
    
    // Stokes velocity
    double v_stokes = (density - air_density) * 9.81 * dp * dp / (18.0 * air_viscosity);
    
    return v_stokes * Cc;
}

// =============================================================================
// RadioactiveDispersalModel Implementation
// =============================================================================

RadioactiveDispersalModel::RadioactiveDispersalModel() {
    // Default constant wind
    wind_func_ = [](double /*t*/, double z, double& u, double& v) {
        u = 10.0;  // 10 m/s eastward
        v = 0.0;
        if (z > 1000.0) u += 0.005 * (z - 1000.0);  // Wind shear
    };
    
    rainfall_func_ = [](double /*x*/, double /*y*/, double /*t*/) { return 0.0; };
}

void RadioactiveDispersalModel::setAtmosphere(std::shared_ptr<ExtendedAtmosphericModel> atm) {
    atmosphere_ = atm;
}

void RadioactiveDispersalModel::setCloudModel(std::shared_ptr<MushroomCloudModel> cloud) {
    cloud_ = cloud;
}

void RadioactiveDispersalModel::setFalloutModel(std::shared_ptr<FalloutFormationModel> fallout) {
    fallout_model_ = fallout;
}

void RadioactiveDispersalModel::setInventory(std::shared_ptr<FissionProductInventory> inventory) {
    inventory_ = inventory;
}

void RadioactiveDispersalModel::setConstantWind(double speed, double direction) {
    double u = speed * std::sin(direction * M_PI / 180.0);
    double v = speed * std::cos(direction * M_PI / 180.0);
    
    wind_func_ = [u, v](double /*t*/, double /*z*/, double& u_out, double& v_out) {
        u_out = u;
        v_out = v;
    };
}

void RadioactiveDispersalModel::setWindProfile(std::function<void(double z, double& u, double& v)> wind) {
    wind_func_ = [wind](double /*t*/, double z, double& u, double& v) {
        wind(z, u, v);
    };
}

void RadioactiveDispersalModel::setTimeVaryingWind(
    std::function<void(double t, double z, double& u, double& v)> wind) {
    wind_func_ = wind;
    time_varying_wind_ = true;
}

void RadioactiveDispersalModel::setRainfall(double rate_mm_hr) {
    constant_rainfall_ = rate_mm_hr;
    rainfall_func_ = [rate_mm_hr](double /*x*/, double /*y*/, double /*t*/) { 
        return rate_mm_hr; 
    };
}

void RadioactiveDispersalModel::setDispersionScheme(const std::string& scheme) {
    dispersion_scheme_ = scheme;
}

void RadioactiveDispersalModel::setPasquillStability(char stability_class) {
    stability_class_ = stability_class;
}

void RadioactiveDispersalModel::setTurbulentDiffusivity(double Kh, double Kv) {
    Kh_ = Kh;
    Kv_ = Kv;
}

void RadioactiveDispersalModel::setDryDepositionVelocity(double vd) {
    dry_deposition_velocity_ = vd;
}

void RadioactiveDispersalModel::setWashoutCoefficient(double lambda) {
    washout_coefficient_ = lambda;
}

void RadioactiveDispersalModel::initialize() {
    particles_.clear();
    depositions_.clear();
    
    // Create grid for Eulerian storage
    int nx = 200, ny = 200, nz = 50;
    double Lx = 200000.0, Ly = 200000.0, Lz = 20000.0;
    
    x_grid_.resize(nx);
    y_grid_.resize(ny);
    z_grid_.resize(nz);
    
    for (int i = 0; i < nx; ++i) x_grid_[i] = -Lx/2 + i * Lx / (nx - 1);
    for (int j = 0; j < ny; ++j) y_grid_[j] = -Ly/2 + j * Ly / (ny - 1);
    for (int k = 0; k < nz; ++k) z_grid_[k] = k * Lz / (nz - 1);
    
    concentration_grid_.resize(nx, std::vector<std::vector<double>>(
        ny, std::vector<double>(nz, 0.0)));
    deposition_grid_.resize(nx, std::vector<double>(ny, 0.0));
    
    // Initialize particles from fallout model
    if (fallout_model_) {
        auto fp = fallout_model_->generateParticles(10000);
        
        for (const auto& f : fp) {
            TrackedParticle p;
            // Initial position in cloud
            if (cloud_) {
                // Note: MushroomCloudModel methods would be called here if the class were defined
                // For now, using default positioning since the class is forward-declared only
                p.x = 0.0;
                p.y = 0.0;
                p.z = 10000.0;  // Default cloud height
            } else {
                p.x = 0.0;
                p.y = 0.0;
                p.z = 10000.0;  // Default cloud height
            }
            
            p.mass = f.mass;
            p.activity = f.activity;
            p.diameter = f.diameter;
            p.deposited = false;
            
            particles_.push_back(p);
        }
    }
    
    current_time_ = 0.0;
}

void RadioactiveDispersalModel::step(double dt) {
    advectParticles(dt);
    diffuseParticles(dt);
    settleParticles(dt);
    depositParticles(dt);
    
    current_time_ += dt;
}

void RadioactiveDispersalModel::runTo(double time) {
    double dt = 60.0;  // 1 minute time step
    while (current_time_ < time) {
        step(dt);
    }
}

void RadioactiveDispersalModel::advectParticles(double dt) {
    for (auto& p : particles_) {
        if (p.deposited) continue;
        
        double u, v;
        wind_func_(current_time_, p.z, u, v);
        
        p.x += u * dt;
        p.y += v * dt;
    }
}

void RadioactiveDispersalModel::diffuseParticles(double dt) {
    std::normal_distribution<double> normal(0.0, 1.0);
    std::mt19937 rng(std::random_device{}());
    
    double sigma_h = std::sqrt(2.0 * Kh_ * dt);
    double sigma_v = std::sqrt(2.0 * Kv_ * dt);
    
    for (auto& p : particles_) {
        if (p.deposited) continue;
        
        p.x += sigma_h * normal(rng);
        p.y += sigma_h * normal(rng);
        p.z += sigma_v * normal(rng);
        
        // Reflect at ground
        if (p.z < 0) p.z = -p.z;
    }
}

void RadioactiveDispersalModel::settleParticles(double dt) {
    for (auto& p : particles_) {
        if (p.deposited) continue;
        
        // Settling velocity (Stokes for small particles)
        // Note: ExtendedAtmosphericModel methods would be called here if the class were defined
        // For now, using default density since the class is forward-declared only
        double rho_air = AIR_DENSITY_STP;  // Default to STP density
        double mu_air = 1.8e-5;
        double rho_p = 2500.0;  // Particle density
        
        double v_settle;
        if (p.diameter < 50e-6) {
            // Stokes
            v_settle = (rho_p - rho_air) * 9.81 * p.diameter * p.diameter / (18.0 * mu_air);
        } else {
            // Larger particles
            v_settle = 0.1 * p.diameter * 1e4;  // Rough approximation
        }
        
        p.z -= v_settle * dt;
    }
}

void RadioactiveDispersalModel::depositParticles(double dt) {
    for (auto& p : particles_) {
        if (p.deposited) continue;
        
        // Check for ground contact
        if (p.z <= 0) {
            p.deposited = true;
            
            // Record deposition
            DepositionRecord rec;
            rec.x = p.x;
            rec.y = p.y;
            rec.time = current_time_;
            rec.mass = p.mass;
            rec.activity = p.activity;
            rec.particle_diameter = p.diameter;
            rec.mechanism = DepositionMechanism::GRAVITATIONAL_SETTLING;
            
            depositions_.push_back(rec);
            
            // Update grid
            int ix = static_cast<int>((p.x - x_grid_.front()) / (x_grid_[1] - x_grid_[0]));
            int iy = static_cast<int>((p.y - y_grid_.front()) / (y_grid_[1] - y_grid_[0]));
            
            if (ix >= 0 && ix < static_cast<int>(x_grid_.size()) &&
                iy >= 0 && iy < static_cast<int>(y_grid_.size())) {
                double cell_area = (x_grid_[1] - x_grid_[0]) * (y_grid_[1] - y_grid_[0]);
                deposition_grid_[ix][iy] += p.activity / cell_area;
            }
            
            continue;
        }
        
        // Dry deposition (probability-based)
        if (p.z < 100.0) {
            double P_deposit = dry_deposition_velocity_ * dt / p.z;
            std::uniform_real_distribution<double> dist(0.0, 1.0);
            std::mt19937 rng(std::random_device{}());
            
            if (dist(rng) < P_deposit) {
                p.deposited = true;
                
                DepositionRecord rec;
                rec.x = p.x;
                rec.y = p.y;
                rec.time = current_time_;
                rec.mass = p.mass;
                rec.activity = p.activity;
                rec.particle_diameter = p.diameter;
                rec.mechanism = DepositionMechanism::DRY_DEPOSITION;
                
                depositions_.push_back(rec);
            }
        }
        
        // Wet deposition (washout by rain)
        double rain_rate = rainfall_func_(p.x, p.y, current_time_);
        if (rain_rate > 0) {
            double lambda = washout_coefficient_ * rain_rate;
            double P_washout = 1.0 - std::exp(-lambda * dt);
            
            std::uniform_real_distribution<double> dist(0.0, 1.0);
            std::mt19937 rng(std::random_device{}());
            
            if (dist(rng) < P_washout) {
                p.deposited = true;
                
                DepositionRecord rec;
                rec.x = p.x;
                rec.y = p.y;
                rec.time = current_time_;
                rec.mass = p.mass;
                rec.activity = p.activity;
                rec.particle_diameter = p.diameter;
                rec.mechanism = DepositionMechanism::WET_DEPOSITION;
                
                depositions_.push_back(rec);
            }
        }
    }
}

double RadioactiveDispersalModel::airConcentration(double x, double y, double z, double /*t*/) const {
    // Find nearby particles and compute concentration
    double concentration = 0.0;
    double sigma = 1000.0;  // Kernel width
    
    for (const auto& p : particles_) {
        if (p.deposited) continue;
        
        double dx = x - p.x;
        double dy = y - p.y;
        double dz = z - p.z;
        double r2 = dx*dx + dy*dy + dz*dz;
        
        // Gaussian kernel
        double kernel = std::exp(-r2 / (2.0 * sigma * sigma));
        concentration += p.activity * kernel / std::pow(sigma * std::sqrt(2.0 * M_PI), 3);
    }
    
    return concentration;
}

double RadioactiveDispersalModel::integratedAirConcentration(double x, double y, double z,
                                                             double t_start, double t_end) const {
    // Integrate air concentration over time using trapezoidal rule
    // This provides the time-integrated concentration (Bq·s/m³) needed for dose calculations
    
    const int num_steps = 100;  // Number of integration steps
    double dt = (t_end - t_start) / num_steps;
    
    if (dt <= 0) {
        return 0.0;
    }
    
    double integrated = 0.0;
    double conc_prev = airConcentration(x, y, z, t_start);
    
    for (int i = 1; i <= num_steps; ++i) {
        double t = t_start + i * dt;
        double conc = airConcentration(x, y, z, t);
        
        // Trapezoidal rule
        integrated += 0.5 * (conc_prev + conc) * dt;
        conc_prev = conc;
    }
    
    return integrated;
}

double RadioactiveDispersalModel::groundDeposition(double x, double y, double /*t*/) const {
    // Find grid cell
    int ix = static_cast<int>((x - x_grid_.front()) / (x_grid_[1] - x_grid_[0]));
    int iy = static_cast<int>((y - y_grid_.front()) / (y_grid_[1] - y_grid_[0]));
    
    if (ix >= 0 && ix < static_cast<int>(x_grid_.size()) &&
        iy >= 0 && iy < static_cast<int>(y_grid_.size())) {
        return deposition_grid_[ix][iy];
    }
    
    return 0.0;
}

double RadioactiveDispersalModel::totalDepositedActivity() const {
    double total = 0.0;
    for (const auto& rec : depositions_) {
        total += rec.activity;
    }
    return total;
}

double RadioactiveDispersalModel::sigmaY(double x, char stability) const {
    // Pasquill-Gifford horizontal dispersion
    double a, b;
    switch (stability) {
        case 'A': a = 0.22; b = 0.0001; break;
        case 'B': a = 0.16; b = 0.0001; break;
        case 'C': a = 0.11; b = 0.0001; break;
        case 'D': a = 0.08; b = 0.0001; break;
        case 'E': a = 0.06; b = 0.0001; break;
        case 'F': a = 0.04; b = 0.0001; break;
        default:  a = 0.08; b = 0.0001; break;
    }
    return a * x * std::pow(1.0 + b * x, -0.5);
}

double RadioactiveDispersalModel::sigmaZ(double x, char stability) const {
    // Pasquill-Gifford vertical dispersion
    double a, b;
    switch (stability) {
        case 'A': a = 0.20; b = 0.0; break;
        case 'B': a = 0.12; b = 0.0; break;
        case 'C': a = 0.08; b = 0.0002; break;
        case 'D': a = 0.06; b = 0.0015; break;
        case 'E': a = 0.03; b = 0.0003; break;
        case 'F': a = 0.016; b = 0.0003; break;
        default:  a = 0.06; b = 0.0015; break;
    }
    return a * x * std::pow(1.0 + b * x, -0.5);
}

void RadioactiveDispersalModel::writeDepositionField(const std::string& filename, double /*t*/) const {
    std::ofstream out(filename);
    
    out << "# Fallout Deposition Field\n";
    out << "# Time: " << current_time_ << " s\n";
    out << "# x(m)  y(m)  activity(Bq/m²)\n";
    
    for (size_t i = 0; i < x_grid_.size(); ++i) {
        for (size_t j = 0; j < y_grid_.size(); ++j) {
            if (deposition_grid_[i][j] > 0) {
                out << x_grid_[i] << "  " << y_grid_[j] << "  " 
                    << deposition_grid_[i][j] << "\n";
            }
        }
    }
}

// =============================================================================
// DoseCalculationModel Implementation
// =============================================================================

DoseCalculationModel::DoseCalculationModel() {}

void DoseCalculationModel::setPromptModel(std::shared_ptr<PromptRadiationModel> prompt) {
    prompt_model_ = prompt;
}

void DoseCalculationModel::setDispersalModel(std::shared_ptr<RadioactiveDispersalModel> dispersal) {
    dispersal_model_ = dispersal;
}

void DoseCalculationModel::setInventory(std::shared_ptr<FissionProductInventory> inventory) {
    inventory_ = inventory;
}

void DoseCalculationModel::setDatabase(std::shared_ptr<RadionuclideDatabase> db) {
    database_ = db;
}

void DoseCalculationModel::setBreathingRate(double rate) {
    breathing_rate_ = rate;
}

void DoseCalculationModel::setExposureDuration(double duration) {
    exposure_duration_ = duration;
}

void DoseCalculationModel::setShielding(const ShieldingConfig& shielding) {
    shielding_ = shielding;
}

DoseCalculationModel::DoseResult 
DoseCalculationModel::calculateDose(double x, double y, double z,
                                    double t_start, double t_end) const {
    DoseResult result;
    
    // Distance from ground zero (assumed at origin)
    double r = std::sqrt(x*x + y*y + z*z);
    
    // Prompt doses (if time includes t=0)
    if (t_start <= 0 && prompt_model_) {
        result.prompt_gamma_dose = promptGammaDose(r);
        result.prompt_neutron_dose = promptNeutronDose(r);
    }
    
    // Cloud shine dose
    result.cloud_shine_dose = cloudShineDose(x, y, z, t_start, t_end);
    
    // Ground shine dose
    result.ground_shine_dose = groundShineDose(x, y, t_start, t_end);
    
    // Inhalation dose
    result.inhalation_dose = inhalationDose(x, y, z, t_start, t_end);
    
    // Thyroid dose (I-131 specific)
    result.thyroid_dose = thyroidInhalationDose(x, y, z, t_start, t_end);
    
    return result;
}

double DoseCalculationModel::promptGammaDose(double r) const {
    if (!prompt_model_) return 0.0;
    
    double dose = prompt_model_->gammaDose(r);
    
    // Apply shielding
    if (shielding_.wall_thickness > 0) {
        dose *= gammaShieldingFactor(shielding_.wall_thickness,
                                    shielding_.wall_density, 1.0);
    }
    
    return dose;
}

double DoseCalculationModel::promptNeutronDose(double r) const {
    if (!prompt_model_) return 0.0;
    return prompt_model_->neutronDose(r);
}

double DoseCalculationModel::cloudShineDose(double x, double y, double z,
                                            double t_start, double t_end) const {
    if (!dispersal_model_) return 0.0;
    return integrateCloudShine(x, y, z, t_start, t_end);
}

double DoseCalculationModel::groundShineDose(double x, double y,
                                             double t_start, double t_end) const {
    if (!dispersal_model_) return 0.0;
    return integrateGroundShine(x, y, t_start, t_end);
}

double DoseCalculationModel::inhalationDose(double x, double y, double z,
                                           double t_start, double t_end) const {
    if (!dispersal_model_) return 0.0;
    return integrateInhalation(x, y, z, t_start, t_end);
}

double DoseCalculationModel::thyroidInhalationDose(double x, double y, double z,
                                                   double t_start, double t_end) const {
    // Thyroid dose from I-131 inhalation
    if (!dispersal_model_ || !database_) return 0.0;
    
    const Radionuclide* I131 = database_->get("I-131");
    if (!I131) return 0.0;
    
    double integrated_conc = dispersal_model_->integratedAirConcentration(x, y, z, t_start, t_end);
    double intake = integrated_conc * breathing_rate_;
    
    // Thyroid dose coefficient for I-131: ~2×10⁻⁷ Sv/Bq
    return intake * 2.0e-7;
}

double DoseCalculationModel::totalDoseRate(double x, double y, double z, double t) const {
    return groundShineDoseRate(x, y, t) + cloudShineDoseRate(x, y, z, t);
}

double DoseCalculationModel::groundShineDoseRate(double x, double y, double t) const {
    if (!dispersal_model_) return 0.0;
    
    double deposition = dispersal_model_->groundDeposition(x, y, t);  // Bq/m²
    
    // Dose rate factor ~5×10⁻¹³ (Sv/h)/(Bq/m²) average
    double dose_rate_factor = 5.0e-13;
    
    // Apply ground roughness reduction
    double roughness_factor = groundRoughnessFactor(shielding_.ground_roughness);
    
    return deposition * dose_rate_factor * roughness_factor;
}

double DoseCalculationModel::cloudShineDoseRate(double x, double y, double z, double t) const {
    if (!dispersal_model_) return 0.0;
    
    double concentration = dispersal_model_->airConcentration(x, y, z, t);  // Bq/m³
    
    // Semi-infinite cloud dose rate factor ~4×10⁻¹⁴ (Sv/h)/(Bq/m³) average
    double dose_rate_factor = 4.0e-14;
    
    return concentration * dose_rate_factor;
}

double DoseCalculationModel::gammaShieldingFactor(double thickness, double density, 
                                                  double energy) const {
    // Tenth-value layer in cm for concrete at various energies
    double TVL;
    if (energy < 0.5) {
        TVL = 15.0;
    } else if (energy < 1.5) {
        TVL = 25.0;
    } else {
        TVL = 35.0;
    }
    
    // Adjust for density
    TVL *= 2400.0 / density;
    
    double n_TVL = thickness * 100.0 / TVL;
    return std::pow(0.1, n_TVL);
}

double DoseCalculationModel::buildingProtectionFactor(const ShieldingConfig& config) const {
    // Protection factor from building shielding
    // Simple model: exponential reduction
    
    double PF = 1.0;
    
    // Wall shielding
    if (config.wall_thickness > 0) {
        PF *= gammaShieldingFactor(config.wall_thickness, config.wall_density, 0.7);
    }
    
    // Building factor (accounts for reduced dose inside)
    PF *= config.building_factor;
    
    // Underground
    if (config.underground) {
        double soil_TVL = 0.20;  // ~20 cm per TVL for soil
        PF *= std::pow(0.1, config.underground_depth / soil_TVL);
    }
    
    return PF;
}

double DoseCalculationModel::groundRoughnessFactor(double roughness) const {
    // Ground roughness reduces exposure by ~0.7 for typical terrain
    return 0.5 + 0.5 / roughness;
}

double DoseCalculationModel::integrateGroundShine(double x, double y, 
                                                  double t_start, double t_end) const {
    // Numerical integration of ground shine dose rate
    double dose = 0.0;
    double dt = (t_end - t_start) / 100.0;
    
    for (double t = t_start; t < t_end; t += dt) {
        dose += groundShineDoseRate(x, y, t) * dt;
    }
    
    // Convert from Sv/h·s to Sv
    return dose / 3600.0;
}

double DoseCalculationModel::integrateCloudShine(double x, double y, double z,
                                                 double t_start, double t_end) const {
    double dose = 0.0;
    double dt = (t_end - t_start) / 100.0;
    
    for (double t = t_start; t < t_end; t += dt) {
        dose += cloudShineDoseRate(x, y, z, t) * dt;
    }
    
    return dose / 3600.0;
}

double DoseCalculationModel::integrateInhalation(double x, double y, double z,
                                                 double t_start, double t_end) const {
    double integrated_conc = dispersal_model_->integratedAirConcentration(x, y, z, t_start, t_end);
    double intake = integrated_conc * breathing_rate_;
    
    // Average dose coefficient ~3×10⁻⁸ Sv/Bq
    return intake * 3.0e-8;
}

// =============================================================================
// RadiationUtils Implementation
// =============================================================================

namespace RadiationUtils {

double decayActivity(double A0, double half_life, double time) {
    return A0 * std::exp(-LN2 * time / half_life);
}

double wayWignerActivity(double fissions, double time) {
    if (time < 0.1) time = 0.1;
    return 1.4e23 * fissions * std::pow(time, -1.2) / FISSIONS_PER_KT;
}

double wayWignerPower(double fissions, double time) {
    if (time < 0.1) time = 0.1;
    return 2.66e6 * fissions * std::pow(time, -1.2) / FISSIONS_PER_KT;
}

double stokesSettlingVelocity(double diameter, double rho_particle,
                              double rho_air, double mu_air) {
    return (rho_particle - rho_air) * 9.81 * diameter * diameter / (18.0 * mu_air);
}

double pasquillGiffordSigmaY(double x, char stability) {
    double a, p;
    switch (stability) {
        case 'A': a = 0.22; p = 0.894; break;
        case 'B': a = 0.16; p = 0.894; break;
        case 'C': a = 0.11; p = 0.894; break;
        case 'D': a = 0.08; p = 0.894; break;
        case 'E': a = 0.06; p = 0.894; break;
        case 'F': a = 0.04; p = 0.894; break;
        default:  a = 0.08; p = 0.894; break;
    }
    return a * std::pow(x, p);
}

double pasquillGiffordSigmaZ(double x, char stability) {
    double a, p;
    switch (stability) {
        case 'A': a = 0.20; p = 0.894; break;
        case 'B': a = 0.12; p = 0.894; break;
        case 'C': a = 0.08; p = 0.894; break;
        case 'D': a = 0.06; p = 0.894; break;
        case 'E': a = 0.03; p = 0.894; break;
        case 'F': a = 0.016; p = 0.894; break;
        default:  a = 0.06; p = 0.894; break;
    }
    return a * std::pow(x, p);
}

double grayToRem(double Gy) {
    return Gy * 100.0;
}

double sievertToRem(double Sv) {
    return Sv * 100.0;
}

double roentgenToGray(double R) {
    return R * 0.00876;  // For air
}

char stabilityFromMeteo(double wind_speed, double solar_rad, bool night) {
    if (night) {
        if (wind_speed < 3.0) return 'F';
        if (wind_speed < 5.0) return 'E';
        return 'D';
    }
    
    if (solar_rad > 700) {
        if (wind_speed < 2.0) return 'A';
        if (wind_speed < 5.0) return 'B';
        return 'C';
    }
    
    if (solar_rad > 350) {
        if (wind_speed < 3.0) return 'B';
        if (wind_speed < 5.0) return 'C';
        return 'D';
    }
    
    return 'D';  // Default neutral
}

}  // namespace RadiationUtils

// =============================================================================
// Configuration Parser Implementation
// =============================================================================

bool RadiationPhysicsConfig::parse(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Cannot open config file: " << filename << std::endl;
        return false;
    }
    
    std::map<std::string, std::string> config;
    std::string line;
    
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#' || line[0] == '[') continue;
        
        size_t eq = line.find('=');
        if (eq == std::string::npos) continue;
        
        std::string key = line.substr(0, eq);
        std::string value = line.substr(eq + 1);
        
        while (!key.empty() && (key.back() == ' ' || key.back() == '\t')) key.pop_back();
        while (!value.empty() && (value.front() == ' ' || value.front() == '\t')) value.erase(0, 1);
        
        config[key] = value;
    }
    
    return parse(config);
}

bool RadiationPhysicsConfig::parse(const std::map<std::string, std::string>& config) {
    auto getDouble = [&](const std::string& key, double def) {
        auto it = config.find(key);
        if (it != config.end()) {
            try { return std::stod(it->second); } catch (...) {}
        }
        return def;
    };
    
    auto getChar = [&](const std::string& key, char def) {
        auto it = config.find(key);
        if (it != config.end() && !it->second.empty()) {
            return it->second[0];
        }
        return def;
    };
    
    yield_kt_ = getDouble("yield_kt", 100.0);
    fission_fraction_ = getDouble("fission_fraction", 0.5);
    burst_height_ = getDouble("burst_height", 0.0);
    
    auto it = config.find("burst_type");
    if (it != config.end()) burst_type_ = it->second;
    
    it = config.find("dispersion_scheme");
    if (it != config.end()) dispersion_scheme_ = it->second;
    
    wind_speed_ = getDouble("wind_speed", 10.0);
    wind_direction_ = getDouble("wind_direction", 0.0);
    stability_class_ = getChar("stability_class", 'D');
    rainfall_rate_ = getDouble("rainfall_rate", 0.0);
    
    median_particle_size_ = getDouble("median_particle_size", 100e-6);
    geometric_std_ = getDouble("particle_size_gsd", 2.5);
    
    configured_ = true;
    return true;
}

std::unique_ptr<PromptRadiationModel> RadiationPhysicsConfig::createPromptModel() const {
    auto model = std::make_unique<PromptRadiationModel>();
    model->setYield(yield_kt_);
    model->setFissionFraction(fission_fraction_);
    model->setBurstHeight(burst_height_);
    return model;
}

std::unique_ptr<FissionProductInventory> RadiationPhysicsConfig::createInventory() const {
    auto inventory = std::make_unique<FissionProductInventory>();
    inventory->setFissionYield(yield_kt_ * fission_fraction_);
    inventory->initialize();
    return inventory;
}

std::unique_ptr<FalloutFormationModel> RadiationPhysicsConfig::createFalloutModel() const {
    auto model = std::make_unique<FalloutFormationModel>();
    model->setYield(yield_kt_);
    model->setFissionYield(yield_kt_ * fission_fraction_);
    model->setBurstType(burst_type_);
    model->setLogNormalDistribution(median_particle_size_, geometric_std_);
    return model;
}

std::unique_ptr<RadioactiveDispersalModel> RadiationPhysicsConfig::createDispersalModel() const {
    auto model = std::make_unique<RadioactiveDispersalModel>();
    model->setConstantWind(wind_speed_, wind_direction_);
    model->setPasquillStability(stability_class_);
    model->setDispersionScheme(dispersion_scheme_);
    if (rainfall_rate_ > 0) {
        model->setRainfall(rainfall_rate_);
    }
    return model;
}

std::unique_ptr<DoseCalculationModel> RadiationPhysicsConfig::createDoseModel() const {
    return std::make_unique<DoseCalculationModel>();
}

} // namespace FSRM
