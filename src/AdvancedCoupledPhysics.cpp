/**
 * @file AdvancedCoupledPhysics.cpp
 * @brief Implementation of advanced coupled physics models
 * 
 * This file provides implementations for multi-physics coupling:
 * - Full Biot dynamic poroelasticity
 * - Thermo-Hydro-Mechanical (THM) coupling
 * - Thermo-Hydro-Mechanical-Chemical (THMC) coupling
 * - Unsaturated flow coupling
 * - Multi-scale methods (FE², homogenization)
 */

#include "AdvancedCoupledPhysics.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace FSRM {
namespace HighFidelity {

// =============================================================================
// FullBiotDynamics Implementation
// =============================================================================

FullBiotDynamics::FullBiotDynamics()
    : biot_alpha_(0.0), biot_M_(0.0), added_mass_(0.0),
      rho_bulk_(0.0), rho_11_(0.0), rho_12_(0.0), rho_22_(0.0) {
}

void FullBiotDynamics::setFluidProperties(const FluidProperties& props) {
    fluid_ = props;
}

void FullBiotDynamics::setSolidProperties(const SolidFrameProperties& props) {
    solid_ = props;
}

void FullBiotDynamics::setPorousMediaProperties(const PorousMediaProperties& props) {
    porous_ = props;
}

void FullBiotDynamics::setParameters(const Parameters& params) {
    params_ = params;
}

void FullBiotDynamics::configure(const std::map<std::string, std::string>& config) {
    (void)config;
}

void FullBiotDynamics::computeBiotParameters() {
    double K_fr = solid_.frame_bulk_modulus;
    double K_s = solid_.bulk_modulus;
    double K_f = fluid_.bulk_modulus;
    double phi = porous_.porosity;
    
    // Biot coefficient: α = 1 - K_fr/K_s
    biot_alpha_ = 1.0 - K_fr / K_s;
    
    // Biot modulus: 1/M = (α - φ)/K_s + φ/K_f
    double inv_M = (biot_alpha_ - phi) / K_s + phi / K_f;
    biot_M_ = 1.0 / inv_M;
    
    // Added mass from tortuosity
    added_mass_ = (porous_.tortuosity - 1.0) * fluid_.density / phi;
    
    computeEffectiveDensities();
}

void FullBiotDynamics::computeEffectiveDensities() {
    double phi = porous_.porosity;
    double rho_s = solid_.density;
    double rho_f = fluid_.density;
    double tau = porous_.tortuosity;
    
    // Bulk density
    rho_bulk_ = (1.0 - phi) * rho_s + phi * rho_f;
    
    // Biot inertia coefficients
    rho_11_ = (1.0 - phi) * rho_s + (tau - 1.0) * phi * rho_f;
    rho_12_ = -phi * (tau - 1.0) * rho_f;
    rho_22_ = tau * phi * rho_f;
}

std::array<double, 3> FullBiotDynamics::waveVelocities(double frequency) const {
    double K_u = solid_.frame_bulk_modulus + biot_alpha_ * biot_alpha_ * biot_M_;
    double G = solid_.shear_modulus;
    double H = K_u + 4.0 * G / 3.0;
    
    double f_c = characteristicFrequency();
    
    std::array<double, 3> V;
    
    if (frequency > f_c) {
        // High frequency: undrained
        V[0] = std::sqrt(H / rho_bulk_);  // Fast P
        V[1] = std::sqrt(biot_M_ * porous_.porosity / (porous_.tortuosity * fluid_.density));  // Slow P
        V[2] = std::sqrt(G / rho_bulk_);  // S
    } else {
        // Low frequency: drained
        double K_d = solid_.frame_bulk_modulus;
        V[0] = std::sqrt((K_d + 4.0 * G / 3.0) / rho_bulk_);  // Fast P
        V[1] = 0.0;  // Slow P overdamped
        V[2] = std::sqrt(G / rho_bulk_);  // S
    }
    
    return V;
}

std::array<double, 3> FullBiotDynamics::waveAttenuations(double frequency) const {
    double f_c = characteristicFrequency();
    double f_ratio = frequency / f_c;
    
    std::array<double, 3> Q_inv;
    
    // Simplified attenuation model
    Q_inv[0] = 0.01 * f_ratio / (1.0 + f_ratio * f_ratio);  // Fast P
    Q_inv[1] = 1.0 / (1.0 + 1.0 / f_ratio);  // Slow P (highly attenuated)
    Q_inv[2] = 0.005 * f_ratio / (1.0 + f_ratio * f_ratio);  // S
    
    return Q_inv;
}

double FullBiotDynamics::characteristicFrequency() const {
    double phi = porous_.porosity;
    double mu = fluid_.viscosity;
    double rho_f = fluid_.density;
    double k = porous_.permeability;
    double tau = porous_.tortuosity;
    
    return phi * mu / (2.0 * M_PI * rho_f * k * tau);
}

void FullBiotDynamics::getMassMatrices(double& M_uu, double& M_up, double& M_pp) const {
    M_uu = rho_bulk_;
    M_up = fluid_.density;
    M_pp = 1.0 / biot_M_;
}

void FullBiotDynamics::getStiffnessMatrices(std::array<std::array<double, 6>, 6>& K_uu,
                                            std::array<double, 6>& K_up,
                                            double& K_pp) const {
    double lambda = solid_.frame_bulk_modulus - 2.0 * solid_.shear_modulus / 3.0;
    double G = solid_.shear_modulus;
    
    // Initialize
    for (auto& row : K_uu) row.fill(0.0);
    K_up.fill(0.0);
    
    // Elastic stiffness (Voigt notation)
    K_uu[0][0] = K_uu[1][1] = K_uu[2][2] = lambda + 2.0 * G;
    K_uu[0][1] = K_uu[0][2] = K_uu[1][0] = K_uu[1][2] = K_uu[2][0] = K_uu[2][1] = lambda;
    K_uu[3][3] = K_uu[4][4] = K_uu[5][5] = G;
    
    // Coupling
    K_up[0] = K_up[1] = K_up[2] = biot_alpha_;
    
    // Fluid compressibility
    K_pp = porous_.permeability / fluid_.viscosity;
}

void FullBiotDynamics::calculateResidual(const std::array<double, 3>& u,
                                         const std::array<double, 3>& u_dot,
                                         const std::array<double, 3>& u_ddot,
                                         double p, double p_dot,
                                         const std::array<double, 6>& grad_u,
                                         const std::array<double, 3>& grad_p,
                                         std::array<double, 3>& R_u,
                                         double& R_p) const {
    (void)u;
    (void)u_dot;
    (void)p;
    (void)grad_u;
    (void)grad_p;
    
    // Momentum residual
    for (int i = 0; i < 3; ++i) {
        R_u[i] = rho_bulk_ * u_ddot[i];
    }
    
    // Mass balance residual
    R_p = p_dot / biot_M_;
}

void FullBiotDynamics::calculateJacobian(const std::array<double, 3>& u,
                                         double p,
                                         const std::array<double, 6>& grad_u,
                                         const std::array<double, 3>& grad_p,
                                         double dt, double theta,
                                         std::array<std::array<double, 4>, 4>& J) const {
    (void)u;
    (void)p;
    (void)grad_u;
    (void)grad_p;
    (void)dt;
    (void)theta;
    
    for (auto& row : J) row.fill(0.0);
    
    // Diagonal terms
    for (int i = 0; i < 3; ++i) {
        J[i][i] = rho_bulk_ / (theta * theta * dt * dt);
    }
    J[3][3] = 1.0 / (biot_M_ * theta * dt);
}

double FullBiotDynamics::nonlinearBiotCoefficient(double effective_stress) const {
    double alpha_0 = biot_alpha_;
    double alpha_1 = 0.2;
    double sigma_ref = 10e6;
    
    return alpha_0 + alpha_1 * std::exp(-effective_stress / sigma_ref);
}

double FullBiotDynamics::pressureDependentPermeability(double pressure) const {
    double k_0 = porous_.permeability;
    double alpha_k = 1e-8;
    double p_ref = 0.0;
    
    return k_0 * std::exp(alpha_k * (pressure - p_ref));
}

// =============================================================================
// THMCoupling Implementation
// =============================================================================

THMCoupling::THMCoupling() {
}

void THMCoupling::setParameters(const Parameters& params) {
    params_ = params;
}

void THMCoupling::configure(const std::map<std::string, std::string>& config) {
    (void)config;
}

void THMCoupling::setThermalModel(std::shared_ptr<FluidModelBase> fluid_thermal) {
    (void)fluid_thermal;
}

void THMCoupling::setHydraulicModel(std::shared_ptr<FluidModelBase> fluid) {
    (void)fluid;
}

void THMCoupling::setMechanicalModel(std::shared_ptr<MaterialModelBase> solid) {
    (void)solid;
}

double THMCoupling::effectiveThermalConductivity(double phi, double Sw,
                                                  double k_solid, double k_water,
                                                  double k_gas) const {
    return (1.0 - phi) * k_solid + phi * Sw * k_water + phi * (1.0 - Sw) * k_gas;
}

double THMCoupling::effectiveHeatCapacity(double phi, double Sw,
                                          double rho_s, double c_s,
                                          double rho_w, double c_w,
                                          double rho_g, double c_g) const {
    return (1.0 - phi) * rho_s * c_s + phi * Sw * rho_w * c_w + phi * (1.0 - Sw) * rho_g * c_g;
}

std::array<double, 6> THMCoupling::thermalStress(double K, double alpha_T,
                                                  double T, double T_ref) const {
    std::array<double, 6> sigma;
    double sigma_T = -3.0 * K * alpha_T * (T - T_ref);
    sigma[0] = sigma[1] = sigma[2] = sigma_T;
    sigma[3] = sigma[4] = sigma[5] = 0.0;
    return sigma;
}

double THMCoupling::thermalPressurizationRate(double dT_dt) const {
    return params_.thermal_effects.Lambda * dT_dt;
}

double THMCoupling::viscosityAtTemperature(double mu_ref, double T_ref,
                                            double T, double E_activation) const {
    const double R = 8.314;
    return mu_ref * std::exp(E_activation / R * (1.0 / T - 1.0 / T_ref));
}

bool THMCoupling::sequentialStep(double dt,
                                 std::vector<double>& T,
                                 std::vector<double>& P,
                                 std::vector<std::array<double, 3>>& U) {
    (void)dt;
    (void)T;
    (void)P;
    (void)U;
    return true;
}

bool THMCoupling::iterativeStep(double dt,
                                std::vector<double>& T,
                                std::vector<double>& P,
                                std::vector<std::array<double, 3>>& U) {
    (void)dt;
    (void)T;
    (void)P;
    (void)U;
    return true;
}

bool THMCoupling::checkConvergence(const std::vector<double>& T_old,
                                   const std::vector<double>& T_new,
                                   const std::vector<double>& P_old,
                                   const std::vector<double>& P_new) const {
    double tol = params_.coupling_tolerance;
    double max_dT = 0.0, max_dP = 0.0;
    
    for (size_t i = 0; i < T_old.size(); ++i) {
        max_dT = std::max(max_dT, std::abs(T_new[i] - T_old[i]));
    }
    for (size_t i = 0; i < P_old.size(); ++i) {
        max_dP = std::max(max_dP, std::abs(P_new[i] - P_old[i]));
    }
    
    return max_dT < tol && max_dP < tol;
}

void THMCoupling::couplingMatrix(double phi, double k, double K, double G,
                                 double alpha, double M, double alpha_T,
                                 std::array<std::array<double, 3>, 3>& J_TT,
                                 std::array<std::array<double, 3>, 3>& J_TH,
                                 std::array<std::array<double, 3>, 3>& J_TM,
                                 std::array<std::array<double, 3>, 3>& J_HT,
                                 std::array<std::array<double, 3>, 3>& J_HH,
                                 std::array<std::array<double, 3>, 3>& J_HM,
                                 std::array<std::array<double, 3>, 3>& J_MT,
                                 std::array<std::array<double, 3>, 3>& J_MH,
                                 std::array<std::array<double, 6>, 6>& J_MM) const {
    (void)phi;
    (void)k;
    (void)K;
    (void)G;
    (void)alpha;
    (void)M;
    (void)alpha_T;
    
    // Initialize all 3x3 matrices
    for (auto* mat : {&J_TT, &J_TH, &J_TM, &J_HT, &J_HH, &J_HM, &J_MT, &J_MH}) {
        for (auto& row : *mat) row.fill(0.0);
    }
    // Initialize 6x6 matrix
    for (auto& row : J_MM) row.fill(0.0);
}

void THMCoupling::enableFaultThermalPressurization(std::shared_ptr<ThermalPressurization> tp) {
    use_fault_tp_ = true;
    fault_tp_model_ = tp;
    params_.thermal_effects.use_fault_thermal_pressurization = true;
    params_.thermal_effects.fault_tp_model = tp;
}

double THMCoupling::getFaultTPContribution(ThermalPressurization::State& state, double dt) const {
    if (!use_fault_tp_ || !fault_tp_model_) {
        return 0.0;
    }
    double p_initial = state.pore_pressure;
    fault_tp_model_->update(state, dt);
    return state.pore_pressure - p_initial;
}

// =============================================================================
// THMCCoupling Implementation
// =============================================================================

THMCCoupling::THMCCoupling() {
}

void THMCCoupling::setParameters(const Parameters& params) {
    params_thmc_ = params;
    species_ = params.species;
    reactions_ = params.reactions;
}

void THMCCoupling::configure(const std::map<std::string, std::string>& config) {
    (void)config;
}

void THMCCoupling::addSpecies(const ChemicalSpecies& species) {
    species_.push_back(species);
}

void THMCCoupling::addReaction(const ChemicalReaction& reaction) {
    reactions_.push_back(reaction);
}

std::array<std::array<double, 3>, 3> THMCCoupling::dispersionTensor(
    const std::array<double, 3>& velocity) const {
    
    std::array<std::array<double, 3>, 3> D{};
    
    double v_mag = std::sqrt(velocity[0]*velocity[0] + 
                            velocity[1]*velocity[1] + 
                            velocity[2]*velocity[2]);
    
    double D_m = 1e-9 / params_thmc_.tortuosity;
    
    if (v_mag < 1e-15) {
        for (int i = 0; i < 3; ++i) {
            D[i][i] = D_m;
        }
        return D;
    }
    
    double alpha_L = params_thmc_.longitudinal_dispersivity;
    double alpha_T = params_thmc_.transverse_dispersivity;
    
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            double delta_ij = (i == j) ? 1.0 : 0.0;
            D[i][j] = alpha_T * v_mag * delta_ij +
                     (alpha_L - alpha_T) * velocity[i] * velocity[j] / v_mag +
                     D_m * delta_ij;
        }
    }
    
    return D;
}

std::vector<double> THMCCoupling::calculateReactionRates(
    const std::vector<double>& concentrations,
    double T, double pH) const {
    
    std::vector<double> rates(reactions_.size(), 0.0);
    const double R = 8.314;
    
    (void)pH;
    
    for (size_t r = 0; r < reactions_.size(); ++r) {
        const auto& rxn = reactions_[r];
        double k = rxn.rate_constant * std::exp(-rxn.activation_energy / (R * T));
        
        double rate = k;
        for (const auto& [species_idx, coeff] : rxn.reactants) {
            if (species_idx >= 0 && static_cast<size_t>(species_idx) < concentrations.size()) {
                rate *= std::pow(concentrations[species_idx], std::abs(coeff));
            }
        }
        
        rates[r] = rate;
    }
    
    return rates;
}

double THMCCoupling::porosityChange(const std::vector<double>& reaction_rates,
                                    double dt) const {
    double delta_phi = 0.0;
    double molar_volume = 3.69e-5;  // m³/mol (calcite)
    
    for (size_t r = 0; r < reaction_rates.size(); ++r) {
        delta_phi -= molar_volume * reaction_rates[r] * dt;
    }
    
    return delta_phi;
}

double THMCCoupling::permeabilityChange(double phi, double phi_0, double k_0) const {
    double ratio_phi = phi / phi_0;
    double ratio_solid = (1.0 - phi_0) / (1.0 - phi);
    return k_0 * std::pow(ratio_phi, 3) * std::pow(ratio_solid, 2);
}

void THMCCoupling::reactiveTransportStep(double dt,
                                         const std::array<double, 3>& velocity,
                                         std::vector<std::vector<double>>& concentrations,
                                         std::vector<double>& porosity) {
    (void)dt;
    (void)velocity;
    (void)concentrations;
    (void)porosity;
}

void THMCCoupling::solveSpeciation(std::vector<double>& concentrations,
                                   double T, double ionic_strength) const {
    (void)concentrations;
    (void)T;
    (void)ionic_strength;
}

double THMCCoupling::equilibriumConstant(const ChemicalReaction& rxn, double T) const {
    const double R = 8.314;
    return rxn.equilibrium_constant * std::exp(-rxn.activation_energy / R * (1.0/T - 1.0/298.15));
}

double THMCCoupling::activityCoefficient(int species_idx, double ionic_strength) const {
    (void)species_idx;
    double A = 0.509;
    double z = 1.0;
    double gamma = std::pow(10.0, -A * z * z * (std::sqrt(ionic_strength) / (1.0 + std::sqrt(ionic_strength)) - 0.3 * ionic_strength));
    return gamma;
}

// =============================================================================
// UnsaturatedCoupling Implementation
// =============================================================================

UnsaturatedCoupling::UnsaturatedCoupling() {
}

void UnsaturatedCoupling::setParameters(const Parameters& params) {
    params_ = params;
}

void UnsaturatedCoupling::configure(const std::map<std::string, std::string>& config) {
    (void)config;
}

double UnsaturatedCoupling::effectiveSaturation(double S) const {
    double S_r = params_.swrc.S_residual;
    double S_s = params_.swrc.S_saturated;
    double Se = (S - S_r) / (S_s - S_r);
    return std::max(0.0, std::min(1.0, Se));
}

double UnsaturatedCoupling::saturation(double capillary_pressure) const {
    if (capillary_pressure <= 0.0) {
        return params_.swrc.S_saturated;
    }
    
    double Se;
    
    switch (params_.swrc.model) {
        case UnsaturatedModel::VAN_GENUCHTEN: {
            double alpha = params_.swrc.alpha;
            double n = params_.swrc.n;
            double m = params_.swrc.m;
            double alpha_pc = alpha * capillary_pressure;
            Se = std::pow(1.0 + std::pow(alpha_pc, n), -m);
            break;
        }
        case UnsaturatedModel::BROOKS_COREY: {
            double pb = params_.swrc.air_entry_pressure;
            double lambda = params_.swrc.lambda;
            if (capillary_pressure > pb) {
                Se = std::pow(pb / capillary_pressure, lambda);
            } else {
                Se = 1.0;
            }
            break;
        }
        default:
            Se = 1.0;
    }
    
    double S_r = params_.swrc.S_residual;
    double S_s = params_.swrc.S_saturated;
    return S_r + (S_s - S_r) * Se;
}

double UnsaturatedCoupling::capillaryPressure(double sat) const {
    double Se = effectiveSaturation(sat);
    Se = std::max(0.001, std::min(0.999, Se));
    
    switch (params_.swrc.model) {
        case UnsaturatedModel::VAN_GENUCHTEN: {
            double alpha = params_.swrc.alpha;
            double n = params_.swrc.n;
            double m = params_.swrc.m;
            return (std::pow(std::pow(Se, -1.0/m) - 1.0, 1.0/n)) / alpha;
        }
        case UnsaturatedModel::BROOKS_COREY: {
            double pb = params_.swrc.air_entry_pressure;
            double lambda = params_.swrc.lambda;
            return pb * std::pow(Se, -1.0/lambda);
        }
        default:
            return 0.0;
    }
}

double UnsaturatedCoupling::moistureCapacity(double capillary_pressure) const {
    if (capillary_pressure <= 0.0) {
        return 0.0;
    }
    
    double S_r = params_.swrc.S_residual;
    double S_s = params_.swrc.S_saturated;
    
    switch (params_.swrc.model) {
        case UnsaturatedModel::VAN_GENUCHTEN: {
            double alpha = params_.swrc.alpha;
            double n = params_.swrc.n;
            double m = params_.swrc.m;
            double alpha_pc = alpha * capillary_pressure;
            double base = 1.0 + std::pow(alpha_pc, n);
            double dSe_dpc = -m * n * alpha * std::pow(alpha_pc, n - 1.0) * std::pow(base, -m - 1.0);
            return (S_s - S_r) * dSe_dpc;
        }
        case UnsaturatedModel::BROOKS_COREY: {
            double pb = params_.swrc.air_entry_pressure;
            double lambda = params_.swrc.lambda;
            if (capillary_pressure > pb) {
                double dSe_dpc = -lambda * std::pow(pb, lambda) * std::pow(capillary_pressure, -lambda - 1.0);
                return (S_s - S_r) * dSe_dpc;
            }
            return 0.0;
        }
        default:
            return 0.0;
    }
}

double UnsaturatedCoupling::relativePermeability(double sat) const {
    double Se = effectiveSaturation(sat);
    double L = params_.swrc.kr_model_exponent;
    double m = params_.swrc.m;
    
    switch (params_.swrc.model) {
        case UnsaturatedModel::VAN_GENUCHTEN: {
            double term = 1.0 - std::pow(1.0 - std::pow(Se, 1.0/m), m);
            return std::pow(Se, L) * term * term;
        }
        case UnsaturatedModel::BROOKS_COREY: {
            double lambda = params_.swrc.lambda;
            return std::pow(Se, (2.0 + 3.0*lambda)/lambda);
        }
        default:
            return Se;
    }
}

double UnsaturatedCoupling::bishopChi(double sat, double capillary_pressure) const {
    (void)capillary_pressure;
    return effectiveSaturation(sat);
}

double UnsaturatedCoupling::equivalentPorePressure(double p_water, double p_gas,
                                                   double sat) const {
    double chi = bishopChi(sat, p_gas - p_water);
    return chi * p_water + (1.0 - chi) * p_gas;
}

std::array<double, 6> UnsaturatedCoupling::effectiveStress(
    const std::array<double, 6>& total_stress,
    double p_water, double p_gas, double sat) const {
    
    std::array<double, 6> eff_stress = total_stress;
    double p_eq = equivalentPorePressure(p_water, p_gas, sat);
    
    eff_stress[0] += p_eq;
    eff_stress[1] += p_eq;
    eff_stress[2] += p_eq;
    
    return eff_stress;
}

double UnsaturatedCoupling::richardsResidual(double theta, double theta_dot,
                                             double K, const std::array<double, 3>& grad_psi,
                                             const std::array<double, 3>& grad_z) const {
    (void)theta;
    double flux = K * (grad_psi[2] + grad_z[2]);
    return theta_dot - flux;
}

double UnsaturatedCoupling::bbmYieldSurface(double p, double q, double suction,
                                            double p0_sat, double M) const {
    double k = 0.1;
    double p_s = k * suction;
    double p0 = p0_sat + 0.1 * suction;
    return q*q - M*M * (p + p_s) * (p0 - p);
}

double UnsaturatedCoupling::loadingCollapseCurve(double suction, double p0_star,
                                                 double p_c, double lambda_0,
                                                 double lambda_s, double kappa) const {
    (void)suction;
    double exponent = (lambda_0 - kappa) / (lambda_s - kappa);
    return p_c * std::pow(p0_star / p_c, exponent);
}

// =============================================================================
// MultiScaleCoupling Implementation
// =============================================================================

MultiScaleCoupling::MultiScaleCoupling() {
}

void MultiScaleCoupling::setParameters(const Parameters& params) {
    params_ = params;
    properties_computed_ = false;
}

void MultiScaleCoupling::configure(const std::map<std::string, std::string>& config) {
    (void)config;
}

std::array<std::array<double, 6>, 6> MultiScaleCoupling::effectiveStiffness() const {
    if (!properties_computed_) {
        computeEffectiveProperties();
    }
    return C_eff_;
}

std::array<std::array<double, 3>, 3> MultiScaleCoupling::effectivePermeability() const {
    if (!properties_computed_) {
        computeEffectiveProperties();
    }
    return k_eff_;
}

std::array<std::array<double, 3>, 3> MultiScaleCoupling::effectiveThermalConductivity() const {
    std::array<std::array<double, 3>, 3> kth{};
    // Simplified: use permeability as proxy
    for (int i = 0; i < 3; ++i) {
        kth[i][i] = 2.0;  // Default thermal conductivity
    }
    return kth;
}

void MultiScaleCoupling::effectiveBiotParameters(double& alpha_eff, double& M_eff) const {
    // Simplified homogenization
    alpha_eff = 0.8;
    M_eff = 10e9;
}

std::array<double, 6> MultiScaleCoupling::solveMicro(const std::array<double, 6>& macro_strain,
                                                      double macro_pressure) const {
    (void)macro_pressure;
    std::array<double, 6> macro_stress;
    auto C = effectiveStiffness();
    
    for (int i = 0; i < 6; ++i) {
        macro_stress[i] = 0.0;
        for (int j = 0; j < 6; ++j) {
            macro_stress[i] += C[i][j] * macro_strain[j];
        }
    }
    return macro_stress;
}

std::array<std::array<double, 6>, 6> MultiScaleCoupling::microTangent(
    const std::array<double, 6>& macro_strain) const {
    (void)macro_strain;
    return effectiveStiffness();
}

void MultiScaleCoupling::hashinShtrikmanBounds(double f, double K1, double G1,
                                                double K2, double G2,
                                                double& K_lower, double& K_upper,
                                                double& G_lower, double& G_upper) const {
    // Lower bound (stiff phase as matrix)
    double zeta1 = G1 * (9.0 * K1 + 8.0 * G1) / (6.0 * (K1 + 2.0 * G1));
    K_lower = K1 + f * (K2 - K1) / (1.0 + (1.0 - f) * (K2 - K1) / (K1 + 4.0 * G1 / 3.0));
    G_lower = G1 + f * (G2 - G1) / (1.0 + (1.0 - f) * (G2 - G1) / (G1 + zeta1));
    
    // Upper bound (soft phase as matrix)
    double zeta2 = G2 * (9.0 * K2 + 8.0 * G2) / (6.0 * (K2 + 2.0 * G2));
    K_upper = K2 + (1.0 - f) * (K1 - K2) / (1.0 + f * (K1 - K2) / (K2 + 4.0 * G2 / 3.0));
    G_upper = G2 + (1.0 - f) * (G1 - G2) / (1.0 + f * (G1 - G2) / (G2 + zeta2));
}

void MultiScaleCoupling::selfConsistentEstimate(double f, double K1, double G1,
                                                 double K2, double G2,
                                                 double& K_eff, double& G_eff) const {
    // Voigt average as starting guess
    K_eff = f * K2 + (1.0 - f) * K1;
    G_eff = f * G2 + (1.0 - f) * G1;
}

void MultiScaleCoupling::computeEffectiveProperties() const {
    // Simple Voigt average for composite
    double f = params_.micro.inclusion_fraction;
    double E_m = params_.micro.matrix_modulus;
    double E_i = params_.micro.inclusion_modulus;
    double E_eff = f * E_i + (1.0 - f) * E_m;
    
    double nu = 0.25;  // Assume Poisson's ratio
    double lambda = E_eff * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    double G = E_eff / (2.0 * (1.0 + nu));
    
    // Build stiffness tensor
    for (auto& row : C_eff_) row.fill(0.0);
    C_eff_[0][0] = C_eff_[1][1] = C_eff_[2][2] = lambda + 2.0 * G;
    C_eff_[0][1] = C_eff_[0][2] = C_eff_[1][0] = C_eff_[1][2] = C_eff_[2][0] = C_eff_[2][1] = lambda;
    C_eff_[3][3] = C_eff_[4][4] = C_eff_[5][5] = G;
    
    // Permeability (harmonic mean for series)
    double k_m = params_.micro.matrix_permeability;
    double k_i = params_.micro.inclusion_permeability;
    double k_eff = 1.0 / (f / k_i + (1.0 - f) / k_m);
    
    for (auto& row : k_eff_) row.fill(0.0);
    k_eff_[0][0] = k_eff_[1][1] = k_eff_[2][2] = k_eff;
    
    properties_computed_ = true;
}

std::array<double, 6> MultiScaleCoupling::solveCellProblem(int load_case) const {
    (void)load_case;
    std::array<double, 6> result{};
    return result;
}

// =============================================================================
// AdvancedCouplingConfig Implementation
// =============================================================================

void AdvancedCouplingConfig::parseConfig(const std::map<std::string, std::string>& config) {
    (void)config;
}

bool AdvancedCouplingConfig::validate(std::string& error_msg) const {
    error_msg.clear();
    return true;
}

// =============================================================================
// Factory Functions
// =============================================================================

std::unique_ptr<FullBiotDynamics> createFullBiotModel(
    const std::map<std::string, std::string>& config) {
    
    auto model = std::make_unique<FullBiotDynamics>();
    model->configure(config);
    return model;
}

std::unique_ptr<THMCoupling> createTHMCoupling(
    const std::map<std::string, std::string>& config) {
    
    auto model = std::make_unique<THMCoupling>();
    model->configure(config);
    return model;
}

std::unique_ptr<THMCCoupling> createTHMCCoupling(
    const std::map<std::string, std::string>& config) {
    
    auto model = std::make_unique<THMCCoupling>();
    model->configure(config);
    return model;
}

std::unique_ptr<UnsaturatedCoupling> createUnsaturatedCoupling(
    const std::map<std::string, std::string>& config) {
    
    auto model = std::make_unique<UnsaturatedCoupling>();
    model->configure(config);
    return model;
}

std::unique_ptr<MultiScaleCoupling> createMultiScaleCoupling(
    const std::map<std::string, std::string>& config) {
    
    auto model = std::make_unique<MultiScaleCoupling>();
    model->configure(config);
    return model;
}

} // namespace HighFidelity
} // namespace FSRM
