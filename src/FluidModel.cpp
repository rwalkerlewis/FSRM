#include "FluidModel.hpp"
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>

namespace FSRM {

// =============================================================================
// Utility functions
// =============================================================================

static double parseDouble(const std::map<std::string, std::string>& config,
                         const std::string& key, double default_val) {
    auto it = config.find(key);
    if (it != config.end() && !it->second.empty()) {
        try {
            return std::stod(it->second);
        } catch (...) {
            return default_val;
        }
    }
    return default_val;
}

static std::string parseString(const std::map<std::string, std::string>& config,
                               const std::string& key, const std::string& default_val) {
    auto it = config.find(key);
    if (it != config.end() && !it->second.empty()) {
        return it->second;
    }
    return default_val;
}

static std::vector<double> parseDoubleArray(const std::string& str) {
    std::vector<double> result;
    std::stringstream ss(str);
    std::string item;
    while (std::getline(ss, item, ',')) {
        // Trim whitespace
        item.erase(0, item.find_first_not_of(" \t"));
        item.erase(item.find_last_not_of(" \t") + 1);
        if (!item.empty()) {
            try {
                result.push_back(std::stod(item));
            } catch (...) {}
        }
    }
    return result;
}

// =============================================================================
// SinglePhaseFluid Implementation
// =============================================================================

SinglePhaseFluid::SinglePhaseFluid() 
    : FluidModelBase(FluidType::SINGLE_PHASE),
      density_ref(1000.0),
      P_reference(1e5),
      T_reference(293.15),
      viscosity_ref(0.001),
      compressibility(1e-9),
      thermal_expansion(2e-4),
      viscosity_P_coeff(0.0),
      viscosity_T_coeff(0.0),
      viscosity_model("constant") {}

double SinglePhaseFluid::getDensity(double P, double T) const {
    // Compressible fluid with thermal expansion
    double dP = P - P_reference;
    double dT = T - T_reference;
    return density_ref * (1.0 + compressibility * dP - thermal_expansion * dT);
}

double SinglePhaseFluid::getViscosity(double P, double T) const {
    if (viscosity_model == "constant") {
        return viscosity_ref;
    } else if (viscosity_model == "exponential") {
        // Exponential pressure and temperature dependence
        double dP = P - P_reference;
        double dT = T - T_reference;
        return viscosity_ref * std::exp(viscosity_P_coeff * dP - viscosity_T_coeff * dT);
    } else if (viscosity_model == "arrhenius") {
        // Arrhenius temperature dependence
        double activation_energy = 15000.0;  // J/mol (typical for water)
        double R = 8.314;
        return viscosity_ref * std::exp(activation_energy / R * (1.0/T - 1.0/T_reference));
    }
    return viscosity_ref;
}

double SinglePhaseFluid::getCompressibility(double P, double T) const {
    (void)P; (void)T;  // Part of interface - used in pressure-dependent models
    return compressibility;
}

void SinglePhaseFluid::configure(const std::map<std::string, std::string>& config) {
    density_ref = parseDouble(config, "density", 1000.0);
    P_reference = parseDouble(config, "reference_pressure", 1e5);
    T_reference = parseDouble(config, "reference_temperature", 293.15);
    viscosity_ref = parseDouble(config, "viscosity", 0.001);
    compressibility = parseDouble(config, "compressibility", 1e-9);
    thermal_expansion = parseDouble(config, "thermal_expansion", 2e-4);
    viscosity_model = parseString(config, "viscosity_model", "constant");
    viscosity_P_coeff = parseDouble(config, "viscosity_pressure_coeff", 0.0);
    viscosity_T_coeff = parseDouble(config, "viscosity_temperature_coeff", 0.0);
}

void SinglePhaseFluid::setDensity(double rho_ref, double P_ref, double T_ref) {
    density_ref = rho_ref;
    P_reference = P_ref;
    T_reference = T_ref;
}

void SinglePhaseFluid::setViscosity(double mu_ref, double P_ref, double T_ref) {
    (void)P_ref;  // Reserved for future use
    (void)T_ref;  // Reserved for future use
    viscosity_ref = mu_ref;
    // Allow separate reference for viscosity
}

void SinglePhaseFluid::setCompressibility(double c_total) {
    compressibility = c_total;
}

void SinglePhaseFluid::setThermalExpansion(double beta) {
    thermal_expansion = beta;
}

void SinglePhaseFluid::setViscosityModel(const std::string& model) {
    viscosity_model = model;
}

// =============================================================================
// BlackOilFluid Implementation
// =============================================================================

BlackOilFluid::BlackOilFluid()
    : FluidModelBase(FluidType::BLACK_OIL),
      pvt_correlation(PVTCorrelation::STANDING),
      P_std(101325.0),
      T_std(288.7),
      oil_density_std(850.0),
      oil_viscosity_dead(0.005),
      oil_API(35.0),
      oil_compressibility(1.5e-9),
      gas_density_std(1.0),
      gas_gravity(0.7),
      gas_compressibility(1e-8),
      water_density_std(1020.0),
      water_viscosity(0.0005),
      water_salinity(30000.0),
      water_compressibility(4.5e-10),
      Rs_max(150.0),
      Pb(15e6),
      Swc(0.2),
      Sor(0.2),
      Sgc(0.05),
      nw(3.0),
      no(2.0),
      ng(2.0),
      krw_max(0.5),
      kro_max(1.0),
      krg_max(0.8),
      Pc_entry(10000.0),
      lambda_pc(2.0) {}

double BlackOilFluid::getDensity(double P, double T) const {
    (void)T;  // Temperature handled through PVT correlations
    // Return oil density by default
    return getOilDensity(P, getSolutionGOR(P));
}

double BlackOilFluid::getViscosity(double P, double T) const {
    (void)T;  // Temperature handled through PVT correlations
    return getOilViscosity(P, getSolutionGOR(P));
}

double BlackOilFluid::getCompressibility(double P, double T) const {
    (void)P; (void)T;  // Part of interface - used in pressure-dependent models
    return oil_compressibility;
}

void BlackOilFluid::configure(const std::map<std::string, std::string>& config) {
    // PVT correlation
    std::string corr = parseString(config, "pvt_correlation", "STANDING");
    if (corr == "STANDING") pvt_correlation = PVTCorrelation::STANDING;
    else if (corr == "VASQUEZ_BEGGS") pvt_correlation = PVTCorrelation::VASQUEZ_BEGGS;
    else if (corr == "GLASO") pvt_correlation = PVTCorrelation::GLASO;
    
    // Oil properties
    oil_density_std = parseDouble(config, "oil_density_std", 850.0);
    oil_viscosity_dead = parseDouble(config, "oil_viscosity_dead", 0.005);
    oil_API = parseDouble(config, "oil_api", 35.0);
    oil_compressibility = parseDouble(config, "oil_compressibility", 1.5e-9);
    
    // Gas properties
    gas_density_std = parseDouble(config, "gas_density_std", 1.0);
    gas_gravity = parseDouble(config, "gas_gravity", 0.7);
    gas_compressibility = parseDouble(config, "gas_compressibility", 1e-8);
    
    // Water properties
    water_density_std = parseDouble(config, "water_density_std", 1020.0);
    water_viscosity = parseDouble(config, "water_viscosity", 0.0005);
    water_salinity = parseDouble(config, "water_salinity", 30000.0);
    water_compressibility = parseDouble(config, "water_compressibility", 4.5e-10);
    
    // Solution GOR
    Rs_max = parseDouble(config, "solution_gor", 150.0);
    Pb = parseDouble(config, "bubble_point_pressure", 15e6);
    
    // Corey parameters
    Swc = parseDouble(config, "swc", 0.2);
    Sor = parseDouble(config, "sor", 0.2);
    Sgc = parseDouble(config, "sgc", 0.05);
    nw = parseDouble(config, "corey_nw", 3.0);
    no = parseDouble(config, "corey_no", 2.0);
    ng = parseDouble(config, "corey_ng", 2.0);
    krw_max = parseDouble(config, "krw_max", 0.5);
    kro_max = parseDouble(config, "kro_max", 1.0);
    krg_max = parseDouble(config, "krg_max", 0.8);
    
    // Capillary pressure
    Pc_entry = parseDouble(config, "pc_entry", 10000.0);
    lambda_pc = parseDouble(config, "lambda_pc", 2.0);
}

double BlackOilFluid::getOilDensity(double P, double Rs) const {
    // Oil density accounting for dissolved gas
    double Bo = getOilFVF(P, Rs);
    double rho_o_surface = oil_density_std;
    double rho_g_surface = gas_density_std;
    
    // Density = (rho_o_std + Rs * rho_g_std) / Bo
    return (rho_o_surface + Rs * rho_g_surface * 0.001) / Bo;
}

double BlackOilFluid::getGasDensity(double P) const {
    // Ideal gas law approximation with compressibility
    double z = 1.0;  // Compressibility factor (simplified)
    double T = 350.0;  // Reservoir temperature (K)
    double Mw = 16.04;  // Methane molecular weight
    double R = 8.314;
    
    return P * Mw / (z * R * T);
}

double BlackOilFluid::getWaterDensity(double P) const {
    double dP = P - P_std;
    return water_density_std * (1.0 + water_compressibility * dP);
}

double BlackOilFluid::getOilViscosity(double P, double Rs) const {
    switch (pvt_correlation) {
        case PVTCorrelation::STANDING:
            return standingMuoLive(P, Rs);
        case PVTCorrelation::VASQUEZ_BEGGS:
            return vazquesBeggsViscosity(P, Rs);
        default:
            return standingMuoLive(P, Rs);
    }
}

double BlackOilFluid::getGasViscosity(double P) const {
    // Lee-Gonzalez-Eakin correlation (simplified)
    double T = 350.0;  // K
    double rho_g = getGasDensity(P);
    double Mw = 16.04 * gas_gravity;
    
    double K = (9.4 + 0.02 * Mw) * std::pow(T, 1.5) / (209 + 19 * Mw + T);
    double X = 3.5 + 986 / T + 0.01 * Mw;
    double Y = 2.4 - 0.2 * X;
    
    return K * std::exp(X * std::pow(rho_g / 1000.0, Y)) * 1e-7;  // Pa·s
}

double BlackOilFluid::getWaterViscosity(double P) const {
    // Slight pressure dependence
    return water_viscosity * (1.0 - 1.5e-9 * (P - P_std));
}

double BlackOilFluid::getOilFVF(double P, double Rs) const {
    switch (pvt_correlation) {
        case PVTCorrelation::STANDING:
            return standingBo(P, Rs);
        case PVTCorrelation::VASQUEZ_BEGGS:
            return vazquezBeggsBo(P, Rs);
        default:
            return standingBo(P, Rs);
    }
}

double BlackOilFluid::getGasFVF(double P) const {
    // Bg = (z * T * P_std) / (P * T_std)
    double z = 1.0;  // Compressibility factor
    double T = 350.0;
    return (z * T * P_std) / (P * T_std);
}

double BlackOilFluid::getWaterFVF(double P) const {
    double dP = P - P_std;
    return 1.0 + water_compressibility * dP;  // Very close to 1.0
}

double BlackOilFluid::getSolutionGOR(double P) const {
    switch (pvt_correlation) {
        case PVTCorrelation::STANDING:
            return standingRs(P);
        case PVTCorrelation::VASQUEZ_BEGGS:
            return vazquezBeggsRs(P);
        default:
            return standingRs(P);
    }
}

double BlackOilFluid::getBubblePointPressure() const {
    return Pb;
}

double BlackOilFluid::getKrw(double Sw) const {
    if (Sw <= Swc) return 0.0;
    if (Sw >= 1.0 - Sor) return krw_max;
    
    double Sw_norm = (Sw - Swc) / (1.0 - Swc - Sor);
    return krw_max * std::pow(Sw_norm, nw);
}

double BlackOilFluid::getKro(double So, double Sg) const {
    (void)Sg;  // Reserved for three-phase correlations
    if (So <= Sor) return 0.0;
    
    double So_norm = (So - Sor) / (1.0 - Swc - Sor - Sgc);
    if (So_norm > 1.0) So_norm = 1.0;
    if (So_norm < 0.0) return 0.0;
    
    return kro_max * std::pow(So_norm, no);
}

double BlackOilFluid::getKrg(double Sg) const {
    if (Sg <= Sgc) return 0.0;
    if (Sg >= 1.0 - Swc - Sor) return krg_max;
    
    double Sg_norm = (Sg - Sgc) / (1.0 - Swc - Sor - Sgc);
    return krg_max * std::pow(Sg_norm, ng);
}

double BlackOilFluid::getPcow(double Sw) const {
    // Brooks-Corey model
    if (Sw <= Swc) return Pc_entry * 10.0;  // Cap at high value
    
    double Se = (Sw - Swc) / (1.0 - Swc);
    return Pc_entry * std::pow(Se, -1.0/lambda_pc);
}

double BlackOilFluid::getPcog(double Sg) const {
    if (Sg <= 0.0) return 0.0;
    return 0.5 * Pc_entry * std::pow(Sg, -1.0/lambda_pc);
}

void BlackOilFluid::setPVTCorrelation(PVTCorrelation corr) {
    pvt_correlation = corr;
}

void BlackOilFluid::setOilProperties(double rho_std, double mu_dead, double API) {
    oil_density_std = rho_std;
    oil_viscosity_dead = mu_dead;
    oil_API = API;
}

void BlackOilFluid::setGasProperties(double rho_std, double gamma_g) {
    gas_density_std = rho_std;
    gas_gravity = gamma_g;
}

void BlackOilFluid::setWaterProperties(double rho_std, double mu, double salinity) {
    water_density_std = rho_std;
    water_viscosity = mu;
    water_salinity = salinity;
}

void BlackOilFluid::setSolutionGOR(double Rs, double bubble_pt) {
    Rs_max = Rs;
    Pb = bubble_pt;
}

void BlackOilFluid::setCoreyParameters(double swc, double sor, double sgc,
                                       double nw_in, double no_in, double ng_in,
                                       double krw, double kro, double krg) {
    Swc = swc;
    Sor = sor;
    Sgc = sgc;
    nw = nw_in;
    no = no_in;
    ng = ng_in;
    krw_max = krw;
    kro_max = kro;
    krg_max = krg;
}

// Standing correlations
double BlackOilFluid::standingRs(double P) const {
    if (P >= Pb) return Rs_max;
    
    // Standing correlation for Rs
    double T = 350.0;  // K
    double T_F = (T - 273.15) * 9.0/5.0 + 32.0;  // Convert to °F
    double P_psia = P / 6894.76;  // Convert to psia
    
    double yg = gas_gravity;
    double Rs = yg * std::pow(P_psia / 18.2 + 1.4, 1.2048) * 
                std::pow(10.0, 0.0125 * oil_API - 0.00091 * T_F);
    
    return std::min(Rs, Rs_max);
}

double BlackOilFluid::standingBo(double P, double Rs) const {
    double T = 350.0;  // K
    double T_F = (T - 273.15) * 9.0/5.0 + 32.0;
    
    double yg = gas_gravity;
    double yo = oil_density_std / 1000.0;  // Specific gravity
    
    double Bo = 0.9759 + 0.00012 * std::pow(Rs * std::pow(yg/yo, 0.5) + 
                                           1.25 * T_F, 1.2);
    
    // Undersaturated correction
    if (P > Pb) {
        double co = oil_compressibility;
        Bo = Bo * std::exp(-co * (P - Pb));
    }
    
    return Bo;
}

double BlackOilFluid::standingMuoLive(double P, double Rs) const {
    // Beggs-Robinson correlation for live oil viscosity
    double T = 350.0;  // K
    double T_F = (T - 273.15) * 9.0/5.0 + 32.0;
    
    // Dead oil viscosity (Beggs-Robinson)
    double Z = 3.0324 - 0.02023 * oil_API;
    double mu_od = std::pow(10.0, std::pow(10.0, Z) * std::pow(T_F, -1.163)) - 1.0;
    
    // Live oil viscosity
    double A = 10.715 * std::pow(Rs + 100.0, -0.515);
    double B = 5.44 * std::pow(Rs + 150.0, -0.338);
    
    double mu_ob = A * std::pow(mu_od, B);
    
    // Undersaturated correction
    if (P > Pb) {
        double m = 2.6 * std::pow(P / 1e6, 1.187) * std::exp(-11.513 - 8.98e-5 * (P/1e6));
        mu_ob = mu_ob * std::pow(P / Pb, m);
    }
    
    return std::max(mu_ob, 1e-5);  // Minimum viscosity
}

double BlackOilFluid::vazquezBeggsRs(double P) const {
    if (P >= Pb) return Rs_max;
    
    double T = 350.0;
    double T_F = (T - 273.15) * 9.0/5.0 + 32.0;
    double P_psia = P / 6894.76;
    
    double C1, C2, C3;
    if (oil_API <= 30.0) {
        C1 = 0.0362;
        C2 = 1.0937;
        C3 = 25.724;
    } else {
        C1 = 0.0178;
        C2 = 1.187;
        C3 = 23.931;
    }
    
    double Rs = C1 * gas_gravity * std::pow(P_psia, C2) * 
                std::exp(C3 * oil_API / (T_F + 460.0));
    
    return std::min(Rs, Rs_max);
}

double BlackOilFluid::vazquezBeggsBo(double P, double Rs) const {
    (void)P;  // Pressure dependency captured in Rs correlation
    double T = 350.0;
    double T_F = (T - 273.15) * 9.0/5.0 + 32.0;
    
    double C1, C2, C3;
    if (oil_API <= 30.0) {
        C1 = 4.677e-4;
        C2 = 1.751e-5;
        C3 = -1.811e-8;
    } else {
        C1 = 4.67e-4;
        C2 = 1.1e-5;
        C3 = 1.337e-9;
    }
    
    double Bo = 1.0 + C1 * Rs + (T_F - 60.0) * (oil_API / gas_gravity) * 
                (C2 + C3 * Rs);
    
    return Bo;
}

double BlackOilFluid::vazquesBeggsViscosity(double P, double Rs) const {
    // Similar to Standing, using Vasquez-Beggs Rs
    return standingMuoLive(P, Rs);
}

// =============================================================================
// CompositionalFluid Implementation
// =============================================================================

CompositionalFluid::CompositionalFluid(int num_components)
    : FluidModelBase(FluidType::COMPOSITIONAL),
      nc(num_components),
      eos_type(EOSType::PENG_ROBINSON),
      last_P(-1.0),
      last_T(-1.0) {
    
    // Initialize vectors
    component_names.resize(nc);
    Mw.resize(nc, 16.04);
    Tc.resize(nc, 190.6);
    Pc.resize(nc, 4.6e6);
    omega.resize(nc, 0.011);
    Vc.resize(nc, 0.0);
    z_global.resize(nc, 1.0 / nc);
    
    // Initialize BIP matrix
    kij.resize(nc, std::vector<double>(nc, 0.0));
}

double CompositionalFluid::getDensity(double P, double T) const {
    auto result = flash(P, T);
    return result.L * result.rho_L + result.V * result.rho_V;
}

double CompositionalFluid::getViscosity(double P, double T) const {
    auto result = flash(P, T);
    double mu_L = getLiquidViscosity(P, T, result.x);
    double mu_V = getVaporViscosity(P, T, result.y);
    return result.L * mu_L + result.V * mu_V;
}

double CompositionalFluid::getCompressibility(double P, double T) const {
    // Numerical approximation
    double dP = P * 0.001;
    double rho1 = getDensity(P - dP, T);
    double rho2 = getDensity(P + dP, T);
    return (rho2 - rho1) / (2.0 * dP * getDensity(P, T));
}

void CompositionalFluid::configure(const std::map<std::string, std::string>& config) {
    // EOS type
    std::string eos = parseString(config, "eos", "PENG_ROBINSON");
    if (eos == "PENG_ROBINSON" || eos == "PR") {
        eos_type = EOSType::PENG_ROBINSON;
    } else if (eos == "SRK") {
        eos_type = EOSType::SRK;
    }
    
    // Component properties (comma-separated)
    auto mw_vals = parseDoubleArray(parseString(config, "component_mw", ""));
    auto tc_vals = parseDoubleArray(parseString(config, "component_tc", ""));
    auto pc_vals = parseDoubleArray(parseString(config, "component_pc", ""));
    auto omega_vals = parseDoubleArray(parseString(config, "component_omega", ""));
    auto z_vals = parseDoubleArray(parseString(config, "composition", ""));
    
    if (!mw_vals.empty()) {
        nc = static_cast<int>(mw_vals.size());
        Mw = mw_vals;
    }
    if (!tc_vals.empty() && tc_vals.size() == static_cast<size_t>(nc)) Tc = tc_vals;
    if (!pc_vals.empty() && pc_vals.size() == static_cast<size_t>(nc)) Pc = pc_vals;
    if (!omega_vals.empty() && omega_vals.size() == static_cast<size_t>(nc)) omega = omega_vals;
    if (!z_vals.empty() && z_vals.size() == static_cast<size_t>(nc)) z_global = z_vals;
    
    // Resize other vectors
    component_names.resize(nc);
    Vc.resize(nc, 0.0);
    kij.resize(nc, std::vector<double>(nc, 0.0));
}

void CompositionalFluid::addComponent(const std::string& name, double mw, double tc,
                                      double pc, double om, double vc) {
    int idx = component_names.size();
    if (idx >= nc) {
        nc++;
        component_names.resize(nc);
        Mw.resize(nc);
        Tc.resize(nc);
        Pc.resize(nc);
        omega.resize(nc);
        Vc.resize(nc);
        z_global.resize(nc);
        kij.resize(nc, std::vector<double>(nc, 0.0));
    }
    
    component_names[idx] = name;
    Mw[idx] = mw;
    Tc[idx] = tc;
    Pc[idx] = pc;
    omega[idx] = om;
    Vc[idx] = vc;
}

void CompositionalFluid::setComposition(const std::vector<double>& z) {
    if (z.size() == static_cast<size_t>(nc)) {
        z_global = z;
    }
}

void CompositionalFluid::setEOS(EOSType type) {
    eos_type = type;
}

CompositionalFluid::FlashResult CompositionalFluid::flash(double P, double T) const {
    return flash(P, T, z_global);
}

CompositionalFluid::FlashResult CompositionalFluid::flash(double P, double T,
                                                          const std::vector<double>& z) const {
    FlashResult result;
    result.x.resize(nc);
    result.y.resize(nc);
    result.converged = false;
    
    // Wilson K-value initial guess
    std::vector<double> K(nc);
    for (int i = 0; i < nc; i++) {
        K[i] = (Pc[i] / P) * std::exp(5.373 * (1.0 + omega[i]) * (1.0 - Tc[i] / T));
    }
    
    // Rachford-Rice iteration
    double V = 0.5;  // Initial vapor fraction guess
    
    for (int iter = 0; iter < 100; iter++) {
        // Rachford-Rice function and derivative
        double f = 0.0;
        double df = 0.0;
        
        for (int i = 0; i < nc; i++) {
            double denom = 1.0 + V * (K[i] - 1.0);
            f += z[i] * (K[i] - 1.0) / denom;
            df -= z[i] * (K[i] - 1.0) * (K[i] - 1.0) / (denom * denom);
        }
        
        double dV = -f / df;
        V += dV;
        V = std::max(0.0, std::min(1.0, V));
        
        if (std::abs(dV) < 1e-10) {
            result.converged = true;
            break;
        }
        
        // Update K-values using EOS (simplified - full implementation would iterate)
        auto phi_L = getFugacityCoefficients(P, T, z, false);
        auto phi_V = getFugacityCoefficients(P, T, z, true);
        
        for (int i = 0; i < nc; i++) {
            K[i] = phi_L[i] / phi_V[i];
        }
    }
    
    result.V = V;
    result.L = 1.0 - V;
    
    // Calculate phase compositions
    for (int i = 0; i < nc; i++) {
        double denom = 1.0 + V * (K[i] - 1.0);
        result.x[i] = z[i] / denom;
        result.y[i] = K[i] * z[i] / denom;
    }
    
    // Calculate densities
    result.rho_L = getLiquidDensity(P, T, result.x);
    result.rho_V = getVaporDensity(P, T, result.y);
    
    return result;
}

double CompositionalFluid::getLiquidDensity(double P, double T, 
                                            const std::vector<double>& x) const {
    double Z = getCompressibilityFactor(P, T, x, false);
    
    // Average molecular weight
    double Mw_avg = 0.0;
    for (int i = 0; i < nc; i++) {
        Mw_avg += x[i] * Mw[i];
    }
    
    double R = 8.314;
    return P * Mw_avg / (Z * R * T * 1000.0);  // kg/m³
}

double CompositionalFluid::getVaporDensity(double P, double T,
                                           const std::vector<double>& y) const {
    double Z = getCompressibilityFactor(P, T, y, true);
    
    double Mw_avg = 0.0;
    for (int i = 0; i < nc; i++) {
        Mw_avg += y[i] * Mw[i];
    }
    
    double R = 8.314;
    return P * Mw_avg / (Z * R * T * 1000.0);
}

double CompositionalFluid::getLiquidViscosity(double P, double T,
                                              const std::vector<double>& x) const {
    // Lohrenz-Bray-Clark correlation (simplified)
    double rho = getLiquidDensity(P, T, x);
    return 0.001 * std::exp(0.001 * rho);  // Simplified correlation
}

double CompositionalFluid::getVaporViscosity(double P, double T,
                                             const std::vector<double>& y) const {
    // Lee-Gonzalez-Eakin correlation
    double Mw_avg = 0.0;
    for (int i = 0; i < nc; i++) {
        Mw_avg += y[i] * Mw[i];
    }
    
    double rho = getVaporDensity(P, T, y);
    double K = (9.4 + 0.02 * Mw_avg) * std::pow(T, 1.5) / (209.0 + 19.0 * Mw_avg + T);
    double X = 3.5 + 986.0 / T + 0.01 * Mw_avg;
    double Y = 2.4 - 0.2 * X;
    
    return K * std::exp(X * std::pow(rho / 1000.0, Y)) * 1e-7;
}

double CompositionalFluid::getCompressibilityFactor(double P, double T,
                                                    const std::vector<double>& comp,
                                                    bool vapor) const {
    double a_mix = getMixturea(comp, T);
    double b_mix = getMixtureb(comp);
    
    double R = 8.314;
    double A = a_mix * P / (R * R * T * T);
    double B = b_mix * P / (R * T);
    
    // Solve cubic: Z³ - (1-B)Z² + (A-3B²-2B)Z - (AB-B²-B³) = 0
    auto roots = solveCubic(A, B);
    
    // Select appropriate root
    if (vapor) {
        return *std::max_element(roots.begin(), roots.end());
    } else {
        double Z_min = *std::min_element(roots.begin(), roots.end());
        return Z_min > 0.0 ? Z_min : roots[0];
    }
}

std::vector<double> CompositionalFluid::getFugacityCoefficients(double P, double T,
                                                                 const std::vector<double>& comp,
                                                                 bool vapor) const {
    std::vector<double> phi(nc);
    
    double Z = getCompressibilityFactor(P, T, comp, vapor);
    double a_mix = getMixturea(comp, T);
    double b_mix = getMixtureb(comp);
    
    double R = 8.314;
    double A = a_mix * P / (R * R * T * T);
    double B = b_mix * P / (R * T);
    
    for (int i = 0; i < nc; i++) {
        double ai = getPRa(i, T);
        double bi = getPRb(i);
        
        // Sum for mixing rule derivative
        double sum_aij = 0.0;
        for (int j = 0; j < nc; j++) {
            double aj = getPRa(j, T);
            sum_aij += comp[j] * std::sqrt(ai * aj) * (1.0 - kij[i][j]);
        }
        
        double sqrt2 = std::sqrt(2.0);
        double ln_phi = bi / b_mix * (Z - 1.0) - std::log(Z - B)
                       - A / (2.0 * sqrt2 * B) * (2.0 * sum_aij / a_mix - bi / b_mix)
                       * std::log((Z + (1.0 + sqrt2) * B) / (Z + (1.0 - sqrt2) * B));
        
        phi[i] = std::exp(ln_phi);
    }
    
    return phi;
}

double CompositionalFluid::getPRa(int i, double T) const {
    double R = 8.314;
    double kappa = 0.37464 + 1.54226 * omega[i] - 0.26992 * omega[i] * omega[i];
    double alpha = std::pow(1.0 + kappa * (1.0 - std::sqrt(T / Tc[i])), 2);
    return 0.45724 * R * R * Tc[i] * Tc[i] / Pc[i] * alpha;
}

double CompositionalFluid::getPRb(int i) const {
    double R = 8.314;
    return 0.07780 * R * Tc[i] / Pc[i];
}

double CompositionalFluid::getMixturea(const std::vector<double>& comp, double T) const {
    double a_mix = 0.0;
    for (int i = 0; i < nc; i++) {
        for (int j = 0; j < nc; j++) {
            double ai = getPRa(i, T);
            double aj = getPRa(j, T);
            a_mix += comp[i] * comp[j] * std::sqrt(ai * aj) * (1.0 - kij[i][j]);
        }
    }
    return a_mix;
}

double CompositionalFluid::getMixtureb(const std::vector<double>& comp) const {
    double b_mix = 0.0;
    for (int i = 0; i < nc; i++) {
        b_mix += comp[i] * getPRb(i);
    }
    return b_mix;
}

std::vector<double> CompositionalFluid::solveCubic(double A, double B) const {
    // Solve Z³ + c2*Z² + c1*Z + c0 = 0 for Peng-Robinson
    double c2 = -(1.0 - B);
    double c1 = A - 3.0 * B * B - 2.0 * B;
    double c0 = -(A * B - B * B - B * B * B);
    
    // Cardano's method
    double p = c1 - c2 * c2 / 3.0;
    double q = c0 - c1 * c2 / 3.0 + 2.0 * c2 * c2 * c2 / 27.0;
    
    double D = q * q / 4.0 + p * p * p / 27.0;
    
    std::vector<double> roots;
    
    if (D > 0) {
        // One real root
        double u = std::cbrt(-q / 2.0 + std::sqrt(D));
        double v = std::cbrt(-q / 2.0 - std::sqrt(D));
        roots.push_back(u + v - c2 / 3.0);
    } else {
        // Three real roots
        double theta = std::acos(-q / 2.0 * std::sqrt(-27.0 / (p * p * p)));
        double r = 2.0 * std::sqrt(-p / 3.0);
        
        roots.push_back(r * std::cos(theta / 3.0) - c2 / 3.0);
        roots.push_back(r * std::cos((theta + 2.0 * M_PI) / 3.0) - c2 / 3.0);
        roots.push_back(r * std::cos((theta + 4.0 * M_PI) / 3.0) - c2 / 3.0);
    }
    
    return roots;
}

// =============================================================================
// BrineFluid Implementation
// =============================================================================

BrineFluid::BrineFluid()
    : FluidModelBase(FluidType::BRINE),
      salinity(0.035),  // 35 g/kg seawater
      ionic_strength(0.7),
      pure_water_density(1000.0),
      pure_water_viscosity(0.001) {}

double BrineFluid::getDensity(double P, double T) const {
    // Batzle and Wang (1992) correlation
    double S = salinity * 1000.0;  // ppm
    double T_C = T - 273.15;
    double P_MPa = P / 1e6;
    
    double rho_w = pure_water_density + 
                   ((-80.0 * T_C - 3.3 * T_C * T_C + 0.00175 * T_C * T_C * T_C 
                    + 489.0 * P_MPa - 2.0 * T_C * P_MPa) / 1000.0);
    
    double rho_b = rho_w + S * (0.668 + 0.44 * S + 1e-6 * 
                   (300.0 * P_MPa - 2400.0 * P_MPa * S 
                    + T_C * (80.0 + 3.0 * T_C - 3300.0 * S - 13.0 * P_MPa + 47.0 * P_MPa * S)));
    
    return rho_b;
}

double BrineFluid::getViscosity(double P, double T) const {
    (void)P;  // Pressure effect on viscosity is minimal for typical reservoir conditions
    // Kestin et al. (1981) correlation
    double T_C = T - 273.15;
    double m = salinity / 0.058443;  // Molality (NaCl)
    
    double mu_w = pure_water_viscosity * 
                  std::exp(-1.94 + 274.0 / (T_C + 133.0));
    
    double A = 0.0816 * m - 0.0122 * m * m + 0.000128 * m * m * m;
    double B = 0.0184 * m - 0.00317 * m * m + 0.000142 * m * m * m;
    
    return mu_w * (1.0 + A + B * T_C);
}

double BrineFluid::getCompressibility(double P, double T) const {
    (void)P; (void)T;  // Part of interface - simplified model
    // Slightly higher than pure water
    return 4.5e-10 * (1.0 + 0.5 * salinity);
}

void BrineFluid::configure(const std::map<std::string, std::string>& config) {
    salinity = parseDouble(config, "salinity", 0.035);
    ionic_strength = parseDouble(config, "ionic_strength", 0.7);
}

void BrineFluid::setSalinity(double S) {
    salinity = S;
}

void BrineFluid::setIonicStrength(double I) {
    ionic_strength = I;
}

double BrineFluid::getActivityCoefficient() const {
    // Debye-Hückel approximation
    double A = 0.509;  // at 25°C
    return std::pow(10.0, -A * std::sqrt(ionic_strength) / (1.0 + std::sqrt(ionic_strength)));
}

double BrineFluid::getOsmoticPressure(double T) const {
    // Van 't Hoff equation
    double R = 8.314;
    double n = salinity / 0.058443;  // Moles per kg
    return n * R * T * 1000.0;  // Pa
}

// =============================================================================
// CO2Fluid Implementation
// =============================================================================

CO2Fluid::CO2Fluid() : FluidModelBase(FluidType::CO2) {}

double CO2Fluid::getDensity(double P, double T) const {
    // Span-Wagner EOS (simplified)
    // Reduced temperature Tr = T / Tc_CO2 is used within getReducedDensity
    double rho_r = getReducedDensity(P, T);
    return rho_r * rhoc_CO2;
}

double CO2Fluid::getViscosity(double P, double T) const {
    // Fenghour et al. (1998) correlation (simplified)
    double rho = getDensity(P, T);
    double Tr_vis = T / Tc_CO2;
    
    // Zero-density viscosity
    double mu_0 = 1.0065e-6 * std::sqrt(T) / (1.0 + 0.24 / Tr_vis);
    
    // Excess viscosity
    double rho_r = rho / rhoc_CO2;
    double mu_ex = 2.2e-6 * rho_r * rho_r;
    
    return mu_0 + mu_ex;
}

double CO2Fluid::getCompressibility(double P, double T) const {
    double dP = P * 0.001;
    double rho1 = getDensity(P - dP, T);
    double rho2 = getDensity(P + dP, T);
    return (rho2 - rho1) / (2.0 * dP * getDensity(P, T));
}

void CO2Fluid::configure(const std::map<std::string, std::string>& config) {
    (void)config;  // CO2 properties are fixed, configuration reserved for future extensions
    // CO2 has fixed critical properties
}

bool CO2Fluid::isSupercritical(double P, double T) const {
    return P > Pc_CO2 && T > Tc_CO2;
}

bool CO2Fluid::isLiquid(double P, double T) const {
    if (isSupercritical(P, T)) return false;
    // Simplified vapor pressure curve
    double Psat = Pc_CO2 * std::exp(7.0 * (1.0 - Tc_CO2 / T));
    return P > Psat;
}

bool CO2Fluid::isGas(double P, double T) const {
    return !isLiquid(P, T) && !isSupercritical(P, T);
}

double CO2Fluid::getSolubilityInWater(double P, double T) const {
    // Henry's law with temperature dependence
    double H = 3.4e9 * std::exp(-2400.0 * (1.0/T - 1.0/298.15));  // Pa
    return P / H;  // Mole fraction
}

double CO2Fluid::getSolubilityInOil(double P, double T) const {
    // Higher than in water
    return getSolubilityInWater(P, T) * 5.0;
}

double CO2Fluid::getReducedDensity(double P, double T) const {
    // Iterative solution using PR-like approach
    double Tr = T / Tc_CO2;
    double Pr = P / Pc_CO2;
    
    // Initial guess
    double rho_r = Pr / (Tr * 0.27);  // Ideal gas approximation
    
    // Newton iteration (simplified)
    for (int i = 0; i < 10; i++) {
        double f = Pr - rho_r * Tr * (1.0 + 0.1 * rho_r - 0.01 * rho_r * rho_r);
        double df = -Tr * (1.0 + 0.2 * rho_r - 0.03 * rho_r * rho_r);
        rho_r -= f / df;
    }
    
    return rho_r;
}

// =============================================================================
// Factory Function
// =============================================================================

std::unique_ptr<FluidModelBase> createFluidModel(
    const std::string& type_str,
    const std::map<std::string, std::string>& config) {
    
    FluidType type = parseFluidType(type_str);
    std::unique_ptr<FluidModelBase> model;
    
    switch (type) {
        case FluidType::SINGLE_PHASE:
            model = std::make_unique<SinglePhaseFluid>();
            break;
        case FluidType::BLACK_OIL:
            model = std::make_unique<BlackOilFluid>();
            break;
        case FluidType::COMPOSITIONAL:
            model = std::make_unique<CompositionalFluid>();
            break;
        case FluidType::BRINE:
            model = std::make_unique<BrineFluid>();
            break;
        case FluidType::CO2:
            model = std::make_unique<CO2Fluid>();
            break;
        default:
            model = std::make_unique<SinglePhaseFluid>();
    }
    
    model->configure(config);
    return model;
}

FluidType parseFluidType(const std::string& type_str) {
    std::string s = type_str;
    std::transform(s.begin(), s.end(), s.begin(), ::toupper);
    
    if (s == "SINGLE_PHASE" || s == "SINGLE" || s == "SINGLE_COMPONENT") {
        return FluidType::SINGLE_PHASE;
    } else if (s == "BLACK_OIL" || s == "BLACKOIL") {
        return FluidType::BLACK_OIL;
    } else if (s == "COMPOSITIONAL" || s == "EOS") {
        return FluidType::COMPOSITIONAL;
    } else if (s == "BRINE" || s == "SALINE") {
        return FluidType::BRINE;
    } else if (s == "CO2" || s == "CARBON_DIOXIDE") {
        return FluidType::CO2;
    } else if (s == "DEAD_OIL") {
        return FluidType::DEAD_OIL;
    } else if (s == "DRY_GAS") {
        return FluidType::DRY_GAS;
    }
    
    return FluidType::SINGLE_PHASE;
}

} // namespace FSRM
