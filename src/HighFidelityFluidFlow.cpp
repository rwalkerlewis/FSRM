#include "HighFidelityFluidFlow.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>

namespace FSRM {
namespace HighFidelity {

// =============================================================================
// ForchheimerFlow Implementation
// =============================================================================

ForchheimerFlow::ForchheimerFlow() {}

void ForchheimerFlow::setParameters(const Parameters& params) {
    params_ = params;
}

void ForchheimerFlow::configure(const std::map<std::string, std::string>& config) {
    auto it = config.find("beta");
    if (it != config.end()) params_.beta = std::stod(it->second);
    
    it = config.find("min_permeability");
    if (it != config.end()) params_.min_permeability = std::stod(it->second);
    
    it = config.find("critical_reynolds");
    if (it != config.end()) params_.critical_Re = std::stod(it->second);
    
    it = config.find("beta_model");
    if (it != config.end()) {
        const std::string& model = it->second;
        if (model == "constant") params_.beta_model = BetaCoefficientModel::CONSTANT;
        else if (model == "cooke") params_.beta_model = BetaCoefficientModel::COOKE;
        else if (model == "geertsma") params_.beta_model = BetaCoefficientModel::GEERTSMA;
        else if (model == "frederick_graves") params_.beta_model = BetaCoefficientModel::FREDERICK_GRAVES;
        else if (model == "evans_civan") params_.beta_model = BetaCoefficientModel::EVANS_CIVAN;
        else if (model == "thauvin_mohanty") params_.beta_model = BetaCoefficientModel::THAUVIN_MOHANTY;
    }
}

double ForchheimerFlow::calculateVelocity(double grad_P, double k, double mu, double rho) const {
    if (std::abs(grad_P) < 1e-30) return 0.0;
    return solveForchheimer(std::abs(grad_P), k, mu, rho) * (grad_P > 0 ? 1.0 : -1.0);
}

double ForchheimerFlow::solveForchheimer(double grad_P, double k, double mu, double rho) const {
    // Solve: grad_P = (mu/k)*v + beta*rho*v^2
    // Quadratic: beta*rho*v^2 + (mu/k)*v - grad_P = 0
    
    double a = params_.beta * rho;
    double b = mu / k;
    double c = -grad_P;
    
    if (std::abs(a) < 1e-30) {
        // Linear (Darcy) case
        return grad_P * k / mu;
    }
    
    // Quadratic formula: v = (-b + sqrt(b^2 + 4*a*grad_P)) / (2*a)
    double discriminant = b * b - 4.0 * a * c;
    if (discriminant < 0) return 0.0;
    
    return (-b + std::sqrt(discriminant)) / (2.0 * a);
}

double ForchheimerFlow::getApparentPermeability(double k, double mu, double rho, double velocity) const {
    double Fo = getForchheimerNumber(k, mu, rho, velocity);
    return k / (1.0 + Fo);
}

double ForchheimerFlow::getForchheimerNumber(double k, double mu, double rho, double velocity) const {
    return params_.beta * rho * k * std::abs(velocity) / mu;
}

double ForchheimerFlow::calculateBeta(double k, double phi) const {
    switch (params_.beta_model) {
        case BetaCoefficientModel::CONSTANT:
            return params_.beta;
        case BetaCoefficientModel::COOKE:
            return betaCooke(k);
        case BetaCoefficientModel::GEERTSMA:
            return betaGeertsma(k, phi);
        case BetaCoefficientModel::FREDERICK_GRAVES:
            return betaFrederickGraves(k);
        case BetaCoefficientModel::EVANS_CIVAN:
            return betaEvansCivan(k, phi);
        case BetaCoefficientModel::THAUVIN_MOHANTY:
            return betaThauvinMohanty(k, phi);
        default:
            return params_.beta;
    }
}

double ForchheimerFlow::betaCooke(double k) const {
    // Cooke (1973): β = 4.85e5 * k^(-0.5)
    // k in m², β in 1/m
    return 4.85e5 * std::pow(k, -0.5);
}

double ForchheimerFlow::betaGeertsma(double k, double phi) const {
    // Geertsma (1974): β = 4.5e5 / (k^0.5 * φ^(5.5))
    return 4.5e5 / (std::sqrt(k) * std::pow(phi, 5.5));
}

double ForchheimerFlow::betaFrederickGraves(double k) const {
    // Frederick & Graves: β = 1.55e5 * k^(-0.449)
    return 1.55e5 * std::pow(k, -0.449);
}

double ForchheimerFlow::betaEvansCivan(double k, double phi) const {
    // Evans & Civan: β = 1.485e9 / (k * φ)
    return 1.485e9 / (k * phi);
}

double ForchheimerFlow::betaThauvinMohanty(double k, double phi) const {
    // Thauvin & Mohanty (1998): β = 1.75e12 * τ / (k * φ)
    // Assuming tortuosity τ ≈ 1/φ for simplicity
    double tau = 1.0 / phi;
    return 1.75e12 * tau / (k * phi);
}

double ForchheimerFlow::getEffectiveMobility(double k,
                                             const std::vector<double>& kr,
                                             const std::vector<double>& mu,
                                             const std::vector<double>& rho,
                                             const std::vector<double>& S,
                                             double grad_P) const {
    // Total mobility without non-Darcy
    double lambda_t = 0.0;
    for (size_t i = 0; i < kr.size(); ++i) {
        lambda_t += kr[i] / mu[i];
    }
    
    // Darcy velocity
    double v_darcy = k * lambda_t * std::abs(grad_P);
    
    // Average density
    double rho_avg = 0.0;
    for (size_t i = 0; i < rho.size(); ++i) {
        rho_avg += S[i] * rho[i];
    }
    
    // Non-Darcy correction
    double Fo = getForchheimerNumber(k, 1.0/lambda_t, rho_avg, v_darcy);
    
    return lambda_t / (1.0 + Fo);
}

double ForchheimerFlow::getReynoldsNumber(double k, double phi, double rho, 
                                          double mu, double velocity) const {
    // Pore Reynolds number: Re = ρ·v·d_p / μ
    // Using hydraulic diameter: d_p ≈ sqrt(k/φ)
    double d_p = std::sqrt(k / phi);
    return rho * std::abs(velocity) * d_p / mu;
}

bool ForchheimerFlow::isNonDarcySignificant(double k, double mu, double rho, 
                                            double velocity, double threshold) const {
    double Fo = getForchheimerNumber(k, mu, rho, velocity);
    return Fo > threshold;
}


// =============================================================================
// KlinkenbergCorrection Implementation
// =============================================================================

KlinkenbergCorrection::KlinkenbergCorrection() {}

void KlinkenbergCorrection::setParameters(const Parameters& params) {
    params_ = params;
}

void KlinkenbergCorrection::configure(const std::map<std::string, std::string>& config) {
    auto it = config.find("slip_factor");
    if (it != config.end()) params_.slip_factor = std::stod(it->second);
    
    it = config.find("reference_pressure");
    if (it != config.end()) params_.reference_pressure = std::stod(it->second);
    
    it = config.find("use_florence_model");
    if (it != config.end()) params_.use_florence_model = (it->second == "true");
    
    it = config.find("pore_diameter");
    if (it != config.end()) params_.pore_diameter = std::stod(it->second);
}

double KlinkenbergCorrection::getApparentPermeability(double k_infinity, double P, 
                                                      double T, double Mw) const {
    if (params_.use_florence_model) {
        // Florence (2007) model considering Knudsen diffusion
        double Kn = getKnudsenNumber(P, T, Mw, params_.pore_diameter);
        double b = calculateSlipFactor(params_.pore_diameter/2, P, T, Mw);
        return k_infinity * (1.0 + b/P) * (1.0 + 4.0*Kn/(1.0 + Kn));
    } else {
        // Standard Klinkenberg: k_app = k_∞ * (1 + b/P)
        return k_infinity * (1.0 + params_.slip_factor / P);
    }
}

double KlinkenbergCorrection::calculateSlipFactor(double pore_radius, double P, 
                                                  double T, double Mw) const {
    // b = 4*c*λ*P/r where c ≈ 1 for straight capillaries
    double lambda = meanFreePath(P, T, Mw);
    double c = 1.0;  // Slip flow constant
    return 4.0 * c * lambda * P / pore_radius;
}

double KlinkenbergCorrection::getKnudsenNumber(double P, double T, double Mw, 
                                               double pore_diameter) const {
    double lambda = meanFreePath(P, T, Mw);
    return lambda / pore_diameter;
}

double KlinkenbergCorrection::meanFreePath(double P, double T, double Mw) const {
    // λ = k_B*T / (√2 * π * d² * P)
    // Using approximation: λ ≈ μ*√(π*R*T/(2*Mw)) / P
    // For air at STP: λ ≈ 68 nm
    const double k_B = 1.38e-23;  // Boltzmann constant
    const double d = 3.7e-10;     // Molecular diameter (N2)
    
    return k_B * T / (std::sqrt(2.0) * M_PI * d * d * P);
}


// =============================================================================
// DualPorosityModel Implementation
// =============================================================================

DualPorosityModel::DualPorosityModel() {
    shape_factor_ = 1.0;
}

void DualPorosityModel::setMatrixProperties(const MatrixProperties& props) {
    matrix_ = props;
    updateShapeFactor();
}

void DualPorosityModel::setFractureProperties(const FractureProperties& props) {
    fracture_ = props;
    updateShapeFactor();
}

void DualPorosityModel::setParameters(const Parameters& params) {
    params_ = params;
    updateShapeFactor();
}

void DualPorosityModel::configure(const std::map<std::string, std::string>& config) {
    // Matrix properties
    auto it = config.find("matrix_porosity");
    if (it != config.end()) matrix_.porosity = std::stod(it->second);
    
    it = config.find("matrix_permeability");
    if (it != config.end()) matrix_.permeability = std::stod(it->second);
    
    it = config.find("matrix_block_size");
    if (it != config.end()) {
        double size = std::stod(it->second);
        matrix_.block_size_x = size;
        matrix_.block_size_y = size;
        matrix_.block_size_z = size;
    }
    
    // Fracture properties
    it = config.find("fracture_porosity");
    if (it != config.end()) fracture_.porosity = std::stod(it->second);
    
    it = config.find("fracture_permeability");
    if (it != config.end()) fracture_.permeability = std::stod(it->second);
    
    it = config.find("fracture_spacing");
    if (it != config.end()) fracture_.spacing = std::stod(it->second);
    
    // Parameters
    it = config.find("dual_permeability");
    if (it != config.end()) params_.dual_permeability = (it->second == "true");
    
    it = config.find("shape_model");
    if (it != config.end()) {
        const std::string& model = it->second;
        if (model == "kazemi") params_.shape_model = ShapeFactorModel::KAZEMI;
        else if (model == "coats") params_.shape_model = ShapeFactorModel::COATS;
        else if (model == "gilman") params_.shape_model = ShapeFactorModel::GILMAN;
        else if (model == "lemonnier") params_.shape_model = ShapeFactorModel::LEMONNIER;
    }
    
    updateShapeFactor();
}

void DualPorosityModel::updateShapeFactor() {
    switch (params_.shape_model) {
        case ShapeFactorModel::KAZEMI:
            shape_factor_ = shapeFactor_Kazemi();
            break;
        case ShapeFactorModel::COATS:
            shape_factor_ = shapeFactor_Coats();
            break;
        case ShapeFactorModel::GILMAN:
            shape_factor_ = shapeFactor_Gilman();
            break;
        case ShapeFactorModel::LEMONNIER:
            shape_factor_ = shapeFactor_Lemonnier();
            break;
        default:
            shape_factor_ = shapeFactor_Kazemi();
    }
    shape_factor_ *= params_.shape_factor_multiplier;
}

double DualPorosityModel::calculateShapeFactor() const {
    return shape_factor_;
}

double DualPorosityModel::shapeFactor_Kazemi() const {
    // Kazemi et al. (1976): σ = 4*(1/Lx² + 1/Ly² + 1/Lz²)
    double Lx = matrix_.block_size_x;
    double Ly = matrix_.block_size_y;
    double Lz = matrix_.block_size_z;
    return 4.0 * (1.0/(Lx*Lx) + 1.0/(Ly*Ly) + 1.0/(Lz*Lz));
}

double DualPorosityModel::shapeFactor_Coats() const {
    // Coats (1989): σ = 8/L² for single block size
    double L = (matrix_.block_size_x + matrix_.block_size_y + matrix_.block_size_z) / 3.0;
    return 8.0 / (L * L);
}

double DualPorosityModel::shapeFactor_Gilman() const {
    // Gilman & Kazemi (1983): σ = 4*(n²+1)/L² for n subdivisions
    double L = matrix_.block_size_x;
    int n = 1;  // Single block
    return 4.0 * (n*n + 1) / (L * L);
}

double DualPorosityModel::shapeFactor_Lemonnier() const {
    // Lemonnier & Bourbiaux: includes gravity correction
    // For now, use Kazemi as base
    return shapeFactor_Kazemi();
}

double DualPorosityModel::calculateTransferRate(double P_matrix, double P_fracture, 
                                                 double mu) const {
    // T_mf = σ · (k_m/μ) · (P_m - P_f)
    return shape_factor_ * (matrix_.permeability / mu) * (P_matrix - P_fracture);
}

double DualPorosityModel::calculatePhaseTransfer(double P_m, double P_f,
                                                 double kr_m, double kr_f,
                                                 double mu, double rho, 
                                                 double depth_diff) const {
    // Upstream weighting for relative permeability
    double kr = (P_m > P_f) ? kr_m : kr_f;
    
    // Pressure difference including gravity
    double dP = P_m - P_f;
    if (params_.include_gravity) {
        const double g = 9.81;
        dP -= rho * g * depth_diff;
    }
    
    return shape_factor_ * (matrix_.permeability * kr / mu) * dP;
}

double DualPorosityModel::calculateCapillaryTransfer(double Sw_m, double Sw_f,
                                                     double Pc_m, double Pc_f,
                                                     double mu_w, double mu_o) const {
    if (!params_.include_capillary) return 0.0;
    
    // Counter-current imbibition
    // Simplified model: transfer proportional to capillary pressure difference
    double dPc = Pc_m - Pc_f;
    
    // Harmonic mean mobility
    double lambda = 2.0 / (mu_w + mu_o);
    
    return shape_factor_ * matrix_.permeability * lambda * dPc;
}

void DualPorosityModel::initializeMINC(MINCState& state, double P_init, double S_init) const {
    int n = params_.num_minc_shells;
    state.shell_pressures.resize(n, P_init);
    state.shell_saturations.resize(n, S_init);
    state.shell_volumes.resize(n);
    
    // Shell volumes (nested spheres)
    double L = matrix_.block_size_x;
    double total_volume = L * L * L;
    
    for (int i = 0; i < n; ++i) {
        double r_outer = L/2 * (1.0 - static_cast<double>(i)/n);
        double r_inner = L/2 * (1.0 - static_cast<double>(i+1)/n);
        state.shell_volumes[i] = (4.0/3.0) * M_PI * (std::pow(r_outer,3) - std::pow(r_inner,3));
    }
}

void DualPorosityModel::updateMINC(MINCState& state, double P_fracture, double S_fracture,
                                    double mu, double dt) const {
    int n = params_.num_minc_shells;
    
    // Flow from outer shells to inner
    for (int i = n-1; i >= 0; --i) {
        double P_outer = (i == n-1) ? P_fracture : state.shell_pressures[i+1];
        double dP = P_outer - state.shell_pressures[i];
        
        // Transfer between shells
        double L_shell = matrix_.block_size_x / n;
        double sigma_shell = 3.0 / (L_shell * L_shell);
        
        double transfer = sigma_shell * (matrix_.permeability / mu) * dP * dt;
        
        // Update pressure (simplified)
        state.shell_pressures[i] += transfer / (matrix_.porosity * matrix_.compressibility);
    }
}

double DualPorosityModel::getMINCTransferRate(const MINCState& state, 
                                               double P_fracture, double mu) const {
    // Transfer from outermost shell to fracture
    int n = params_.num_minc_shells;
    double L_shell = matrix_.block_size_x / n;
    double sigma_shell = 3.0 / (L_shell * L_shell);
    
    return sigma_shell * (matrix_.permeability / mu) * 
           (state.shell_pressures[n-1] - P_fracture);
}


// =============================================================================
// TriplePorosityModel Implementation
// =============================================================================

TriplePorosityModel::TriplePorosityModel() {}

void TriplePorosityModel::setVugProperties(const VugProperties& props) {
    vugs_ = props;
}

void TriplePorosityModel::setDualPorosityBase(const DualPorosityModel& dp) {
    dual_porosity_ = dp;
    
    // Calculate additional shape factors for vug interactions
    double L_vug = std::cbrt(vugs_.porosity / vugs_.connectivity);  // Characteristic vug size
    shape_factor_mv_ = 4.0 / (L_vug * L_vug);
    shape_factor_fv_ = 8.0 / (L_vug * L_vug);  // Higher connectivity to fractures
}

std::array<double, 3> TriplePorosityModel::calculateAllTransfers(double P_m, double P_f, 
                                                                  double P_v, double mu) const {
    std::array<double, 3> transfers;
    
    // Matrix-Fracture
    transfers[0] = dual_porosity_.calculateTransferRate(P_m, P_f, mu);
    
    // Fracture-Vug
    transfers[1] = shape_factor_fv_ * (vugs_.permeability / mu) * (P_f - P_v);
    
    // Matrix-Vug (usually smaller)
    transfers[2] = shape_factor_mv_ * (dual_porosity_.getMatrixProperties().permeability / mu) * 
                   (P_m - P_v) * vugs_.connectivity;
    
    return transfers;
}


// =============================================================================
// NonIsothermalMultiphaseFlow Implementation
// =============================================================================

NonIsothermalMultiphaseFlow::NonIsothermalMultiphaseFlow() {}

void NonIsothermalMultiphaseFlow::setThermalProperties(const ThermalProperties& props) {
    thermal_ = props;
}

void NonIsothermalMultiphaseFlow::setParameters(const Parameters& params) {
    params_ = params;
}

void NonIsothermalMultiphaseFlow::configure(const std::map<std::string, std::string>& config) {
    auto it = config.find("rock_thermal_conductivity");
    if (it != config.end()) thermal_.rock_thermal_conductivity = std::stod(it->second);
    
    it = config.find("rock_specific_heat");
    if (it != config.end()) thermal_.rock_specific_heat = std::stod(it->second);
    
    it = config.find("include_conduction");
    if (it != config.end()) params_.include_conduction = (it->second == "true");
    
    it = config.find("include_convection");
    if (it != config.end()) params_.include_convection = (it->second == "true");
    
    it = config.find("include_phase_change");
    if (it != config.end()) params_.include_phase_change = (it->second == "true");
    
    it = config.find("reference_temperature");
    if (it != config.end()) params_.reference_temperature = std::stod(it->second);
}

double NonIsothermalMultiphaseFlow::getEffectiveConductivity(double phi, double Sw, 
                                                             double So, double Sg) const {
    // Volume-weighted harmonic/arithmetic mean
    double k_rock = thermal_.rock_thermal_conductivity;
    double k_fluid = Sw * thermal_.water_thermal_conductivity +
                     So * thermal_.oil_thermal_conductivity +
                     Sg * thermal_.gas_thermal_conductivity;
    
    // Geometric mean (common for porous media)
    return std::pow(k_rock, 1.0 - phi) * std::pow(k_fluid, phi);
}

double NonIsothermalMultiphaseFlow::getEffectiveHeatCapacity(double phi, double Sw, 
                                                             double So, double Sg,
                                                             double rho_w, double rho_o, 
                                                             double rho_g) const {
    double rhoc_rock = thermal_.rock_density * thermal_.rock_specific_heat;
    double rhoc_water = rho_w * thermal_.water_specific_heat;
    double rhoc_oil = rho_o * thermal_.oil_specific_heat;
    double rhoc_gas = rho_g * thermal_.gas_specific_heat;
    
    double rhoc_fluid = Sw * rhoc_water + So * rhoc_oil + Sg * rhoc_gas;
    
    return (1.0 - phi) * rhoc_rock + phi * rhoc_fluid;
}

std::array<double, 3> NonIsothermalMultiphaseFlow::calculateHeatFlux(
    double phi, double Sw, double So, double Sg,
    const std::array<double, 3>& grad_T,
    const std::array<double, 3>& vel_w,
    const std::array<double, 3>& vel_o,
    const std::array<double, 3>& vel_g,
    double rho_w, double rho_o, double rho_g,
    double T) const {
    
    std::array<double, 3> flux = {0.0, 0.0, 0.0};
    
    // Conductive flux: q_cond = -k_eff * ∇T
    if (params_.include_conduction) {
        double k_eff = getEffectiveConductivity(phi, Sw, So, Sg);
        for (int i = 0; i < 3; ++i) {
            flux[i] -= k_eff * grad_T[i];
        }
    }
    
    // Convective flux: q_conv = Σ ρ_i * h_i * v_i
    if (params_.include_convection) {
        double h_w = calculateEnthalpy(T, thermal_.water_specific_heat);
        double h_o = calculateEnthalpy(T, thermal_.oil_specific_heat);
        double h_g = calculateEnthalpy(T, thermal_.gas_specific_heat, 
                                       thermal_.latent_heat_vaporization);
        
        for (int i = 0; i < 3; ++i) {
            flux[i] += rho_w * h_w * vel_w[i];
            flux[i] += rho_o * h_o * vel_o[i];
            flux[i] += rho_g * h_g * vel_g[i];
        }
    }
    
    return flux;
}

double NonIsothermalMultiphaseFlow::calculateEnthalpy(double T, double cp, double latent) const {
    // h = cp * (T - T_ref) + latent
    return cp * (T - params_.reference_temperature) + latent;
}

double NonIsothermalMultiphaseFlow::getViscosity(double mu_ref, double T_ref, 
                                                  double T, double E_a) const {
    // Andrade equation: μ = A * exp(B/T)
    // μ/μ_ref = exp(E_a/R * (1/T - 1/T_ref))
    const double R = 8.314;  // Gas constant
    return mu_ref * std::exp(E_a / R * (1.0/T - 1.0/T_ref));
}

double NonIsothermalMultiphaseFlow::getDensity(double rho_ref, double T_ref, 
                                               double T, double beta) const {
    return rho_ref / (1.0 + beta * (T - T_ref));
}

double NonIsothermalMultiphaseFlow::calculatePhaseChangeRate(double P, double T, 
                                                             double Sw, double Sg) const {
    if (!params_.include_phase_change) return 0.0;
    
    double P_sat = saturationPressure(T);
    
    // Simple model: evaporation if P < P_sat and liquid present
    if (P < P_sat && Sw > 0.01) {
        // Evaporation rate proportional to departure from equilibrium
        return 1e-3 * Sw * (P_sat - P) / P_sat;
    } else if (P > P_sat && Sg > 0.01) {
        // Condensation rate
        return -1e-3 * Sg * (P - P_sat) / P_sat;
    }
    
    return 0.0;
}

double NonIsothermalMultiphaseFlow::saturationPressure(double T) const {
    // Clausius-Clapeyron approximation
    // ln(P/P_ref) = -L/R * (1/T - 1/T_ref)
    const double P_ref = 1e5;  // 1 atm
    const double T_ref = thermal_.boiling_point;
    const double L = thermal_.latent_heat_vaporization;
    const double R = 461.5;  // Water vapor gas constant
    
    return P_ref * std::exp(-L/R * (1.0/T - 1.0/T_ref));
}

double NonIsothermalMultiphaseFlow::getJouleThomsonCoefficient(double P, double T, 
                                                               double rho, double cp,
                                                               double thermal_expansion) const {
    // μ_JT = (1/cp) * (T*β/ρ - 1/ρ) = (T*β - 1)/(ρ*cp)
    return (T * thermal_expansion - 1.0) / (rho * cp);
}


// =============================================================================
// DynamicRelativePermeability Implementation
// =============================================================================

DynamicRelativePermeability::DynamicRelativePermeability() {}

void DynamicRelativePermeability::setParameters(const Parameters& params) {
    params_ = params;
}

void DynamicRelativePermeability::setBaseRelPerm(std::shared_ptr<ThreePhaseRelPerm> base) {
    base_relperm_ = base;
}

void DynamicRelativePermeability::configure(const std::map<std::string, std::string>& config) {
    auto it = config.find("include_capillary_number");
    if (it != config.end()) params_.include_capillary_number = (it->second == "true");
    
    it = config.find("critical_capillary_number");
    if (it != config.end()) params_.critical_capillary_number = std::stod(it->second);
    
    it = config.find("include_fingering");
    if (it != config.end()) params_.include_fingering = (it->second == "true");
    
    it = config.find("fingering_exponent");
    if (it != config.end()) params_.fingering_exponent = std::stod(it->second);
}

double DynamicRelativePermeability::calculateCapillaryNumber(double mu, double v, 
                                                             double sigma) const {
    return mu * std::abs(v) / sigma;
}

void DynamicRelativePermeability::getModifiedEndpoints(double Ca, double& Sor, 
                                                       double& Swc) const {
    // CDC model: Sor(Ca) = Sor_0 * f(Ca/Ca_c)
    // f(x) = (1 + (1/x)^n)^(-1) → 0 at high Ca
    
    double x = Ca / params_.critical_capillary_number;
    double f = 1.0 / (1.0 + std::pow(1.0/std::max(x, 1e-10), params_.cdc_exponent));
    
    double Sor_0 = 0.2;  // Base residual oil
    double Swc_0 = 0.2;  // Base connate water
    
    Sor = Sor_0 * (1.0 - f) + params_.Sor_at_high_Ca * f;
    Swc = Swc_0 * (1.0 - f) + params_.Swc_at_high_Ca * f;
}

std::array<double, 3> DynamicRelativePermeability::getRelPerm(
    double Sw, double So, double Sg,
    double velocity, double mu_w, double mu_o,
    double sigma) const {
    
    // Get base relative permeabilities
    std::array<double, 3> kr;
    if (base_relperm_) {
        kr = base_relperm_->calculate(Sw, So, Sg);
    } else {
        // Simple Corey model fallback
        double Swc = 0.2, Sor = 0.2;
        double Sw_norm = (Sw - Swc) / (1.0 - Swc - Sor);
        Sw_norm = std::max(0.0, std::min(1.0, Sw_norm));
        kr[0] = std::pow(Sw_norm, 2);
        kr[1] = std::pow(1.0 - Sw_norm, 2);
        kr[2] = 0.0;  // No gas
    }
    
    // Apply capillary number effects
    if (params_.include_capillary_number) {
        double Ca = calculateCapillaryNumber(mu_w, velocity, sigma);
        double Sor_mod, Swc_mod;
        getModifiedEndpoints(Ca, Sor_mod, Swc_mod);
        
        // Scale relative permeabilities with modified endpoints
        // This is a simplified approach
        double endpoint_factor = 1.0 + (Ca / params_.critical_capillary_number);
        endpoint_factor = std::min(endpoint_factor, 2.0);
        
        kr[0] *= endpoint_factor;
        kr[1] *= endpoint_factor;
    }
    
    // Apply fingering correction if mobility ratio is high
    if (params_.include_fingering) {
        double M = getMobilityRatio(kr[0], kr[1], mu_w, mu_o);
        if (M > params_.mobility_ratio_threshold) {
            applyFingeringCorrection(M, kr);
        }
    }
    
    // Rate effects
    if (params_.include_rate_effects) {
        double rate_factor = std::pow(std::abs(velocity) / params_.reference_velocity, 
                                      params_.rate_exponent);
        rate_factor = std::max(0.5, std::min(2.0, rate_factor));
        kr[0] *= rate_factor;
        kr[1] /= rate_factor;  // Opposite effect on oil
    }
    
    return kr;
}

double DynamicRelativePermeability::getEffectiveViscosity(double mu_o, double mu_w, 
                                                          double fw, double omega) const {
    // Todd-Longstaff mixing model
    double mu_mix = std::pow(mu_o, 1.0 - fw) * std::pow(mu_w, fw);
    return std::pow(mu_o, 1.0 - omega) * std::pow(mu_mix, omega);
}

double DynamicRelativePermeability::getMobilityRatio(double kr_w, double kr_o, 
                                                     double mu_w, double mu_o) const {
    if (kr_o < 1e-10 || mu_o < 1e-30) return 1e10;
    return (kr_w / mu_w) / (kr_o / mu_o);
}

void DynamicRelativePermeability::applyFingeringCorrection(double M, 
                                                           std::array<double, 3>& kr) const {
    // Reduce water mobility to account for fingering inefficiency
    double omega = params_.fingering_exponent;
    double correction = std::pow(M, -omega);
    kr[0] *= correction;
}


// =============================================================================
// MiscibleFlowModel Implementation
// =============================================================================

MiscibleFlowModel::MiscibleFlowModel() {}

void MiscibleFlowModel::setParameters(const Parameters& params) {
    params_ = params;
}

void MiscibleFlowModel::configure(const std::map<std::string, std::string>& config) {
    auto it = config.find("MMP");
    if (it != config.end()) params_.MMP = std::stod(it->second);
    
    it = config.find("mixing_parameter");
    if (it != config.end()) params_.mixing_parameter = std::stod(it->second);
    
    it = config.find("longitudinal_dispersivity");
    if (it != config.end()) params_.longitudinal_dispersivity = std::stod(it->second);
    
    it = config.find("transverse_dispersivity");
    if (it != config.end()) params_.transverse_dispersivity = std::stod(it->second);
}

MiscibilityState MiscibleFlowModel::getMiscibilityState(double P) const {
    if (P >= params_.MMP) {
        return MiscibilityState::FIRST_CONTACT;
    } else if (P >= params_.MMP - params_.miscibility_transition_width) {
        return MiscibilityState::NEAR_MISCIBLE;
    } else {
        return MiscibilityState::IMMISCIBLE;
    }
}

double MiscibleFlowModel::getMiscibilityFactor(double P) const {
    if (P >= params_.MMP) return 1.0;
    
    double P_lower = params_.MMP - params_.miscibility_transition_width;
    if (P <= P_lower) return 0.0;
    
    // Smooth transition
    double x = (P - P_lower) / params_.miscibility_transition_width;
    return 3.0*x*x - 2.0*x*x*x;  // Smoothstep
}

void MiscibleFlowModel::getEffectiveProperties(double S_s, double S_o, double P, double T,
                                               double rho_s, double rho_o, 
                                               double mu_s, double mu_o,
                                               double& rho_eff, double& mu_eff) const {
    double F = getMiscibilityFactor(P);
    double omega = params_.mixing_parameter;
    
    if (F < 1e-10) {
        // Immiscible
        rho_eff = (S_s * rho_s + S_o * rho_o) / (S_s + S_o + 1e-10);
        mu_eff = std::pow(mu_s, S_s/(S_s+S_o+1e-10)) * 
                 std::pow(mu_o, S_o/(S_s+S_o+1e-10));
    } else {
        // Miscible or near-miscible
        double x_s = S_s / (S_s + S_o + 1e-10);
        double x_o = 1.0 - x_s;
        
        // Todd-Longstaff for density
        double rho_mix = x_s * rho_s + x_o * rho_o;
        rho_eff = F * rho_mix + (1.0 - F) * (S_s*rho_s + S_o*rho_o)/(S_s+S_o+1e-10);
        
        // Todd-Longstaff for viscosity
        double mu_mix = std::pow(mu_s, x_s) * std::pow(mu_o, x_o);
        double mu_s_eff = std::pow(mu_s, 1.0-omega) * std::pow(mu_mix, omega);
        double mu_o_eff = std::pow(mu_o, 1.0-omega) * std::pow(mu_mix, omega);
        
        mu_eff = F * mu_mix + (1.0 - F) * std::pow(mu_s_eff, x_s) * std::pow(mu_o_eff, x_o);
    }
}

double MiscibleFlowModel::getDispersionCoefficient(double velocity, bool longitudinal) const {
    double alpha = longitudinal ? params_.longitudinal_dispersivity : 
                                  params_.transverse_dispersivity;
    return params_.molecular_diffusion + alpha * std::abs(velocity);
}

std::array<double, 3> MiscibleFlowModel::getRelPermWithMiscibility(
    double Sw, double So, double Sg, double P,
    const std::array<double, 3>& kr_immisc) const {
    
    double F = getMiscibilityFactor(P);
    
    // Miscible relative permeabilities (linear in saturation)
    std::array<double, 3> kr_misc;
    kr_misc[0] = Sw;
    kr_misc[1] = So;
    kr_misc[2] = Sg;
    
    // Interpolate
    std::array<double, 3> kr;
    for (int i = 0; i < 3; ++i) {
        kr[i] = F * kr_misc[i] + (1.0 - F) * kr_immisc[i];
    }
    
    return kr;
}


// =============================================================================
// Configuration Parsing
// =============================================================================

void HighFidelityFlowConfig::parseConfig(const std::map<std::string, std::string>& config) {
    auto it = config.find("enable_high_fidelity_flow");
    if (it != config.end()) enable_high_fidelity = (it->second == "true");
    
    it = config.find("enable_non_darcy");
    if (it != config.end()) enable_non_darcy = (it->second == "true");
    
    it = config.find("enable_dual_porosity");
    if (it != config.end()) enable_dual_porosity = (it->second == "true");
    
    it = config.find("enable_dual_permeability");
    if (it != config.end()) enable_dual_permeability = (it->second == "true");
    
    it = config.find("enable_triple_porosity");
    if (it != config.end()) enable_triple_porosity = (it->second == "true");
    
    it = config.find("enable_non_isothermal");
    if (it != config.end()) enable_non_isothermal = (it->second == "true");
    
    it = config.find("enable_dynamic_relperm");
    if (it != config.end()) enable_dynamic_relperm = (it->second == "true");
    
    it = config.find("enable_miscible");
    if (it != config.end()) enable_miscible = (it->second == "true");
    
    it = config.find("enable_klinkenberg");
    if (it != config.end()) enable_klinkenberg = (it->second == "true");
}

bool HighFidelityFlowConfig::validate(std::string& error_msg) const {
    // Check for incompatible options
    if (enable_dual_permeability && !enable_dual_porosity) {
        error_msg = "Dual permeability requires dual porosity to be enabled";
        return false;
    }
    
    if (enable_triple_porosity && !enable_dual_porosity) {
        error_msg = "Triple porosity requires dual porosity to be enabled";
        return false;
    }
    
    return true;
}


// =============================================================================
// Factory Functions
// =============================================================================

std::unique_ptr<ForchheimerFlow> createNonDarcyModel(
    const std::map<std::string, std::string>& config) {
    auto model = std::make_unique<ForchheimerFlow>();
    model->configure(config);
    return model;
}

std::unique_ptr<DualPorosityModel> createDualPorosityModel(
    const std::map<std::string, std::string>& config) {
    auto model = std::make_unique<DualPorosityModel>();
    model->configure(config);
    return model;
}

std::unique_ptr<NonIsothermalMultiphaseFlow> createThermalFlowModel(
    const std::map<std::string, std::string>& config) {
    auto model = std::make_unique<NonIsothermalMultiphaseFlow>();
    model->configure(config);
    return model;
}

} // namespace HighFidelity
} // namespace FSRM
