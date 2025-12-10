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

namespace fsrm {
namespace coupled {

// =============================================================================
// FullBiotDynamics Implementation
// =============================================================================

FullBiotDynamics::FullBiotDynamics()
    : formulation_(BiotFormulation::U_P), high_frequency_(false),
      rho_solid_(2600.0), rho_fluid_(1000.0), porosity_(0.2),
      permeability_(1e-15), fluid_viscosity_(0.001),
      bulk_modulus_solid_(36e9), shear_modulus_(30e9),
      bulk_modulus_fluid_(2.2e9), bulk_modulus_grain_(40e9),
      biot_coefficient_(0.8), biot_modulus_(10e9), tortuosity_(2.0) {
}

void FullBiotDynamics::computeBiotCoefficients() {
    // Biot coefficient: alpha = 1 - K_drained / K_grain
    // Biot modulus: 1/M = (alpha - n)/K_s + n/K_f
    
    double Kd = bulk_modulus_solid_;
    double Ks = bulk_modulus_grain_;
    double Kf = bulk_modulus_fluid_;
    double n = porosity_;
    
    biot_coefficient_ = 1.0 - Kd / Ks;
    
    double inv_M = (biot_coefficient_ - n) / Ks + n / Kf;
    biot_modulus_ = 1.0 / inv_M;
}

void FullBiotDynamics::computeWaveVelocities(double& Vp_fast, double& Vp_slow,
                                             double& Vs) const {
    // Biot theory wave velocities
    double rho_b = (1.0 - porosity_) * rho_solid_ + porosity_ * rho_fluid_;
    double alpha = biot_coefficient_;
    double M = biot_modulus_;
    double K = bulk_modulus_solid_;
    double G = shear_modulus_;
    double n = porosity_;
    double rho_f = rho_fluid_;
    double tau = tortuosity_;
    
    // Undrained bulk modulus
    double K_u = K + alpha * alpha * M;
    
    // H = K_u + 4G/3
    double H = K_u + 4.0 * G / 3.0;
    
    // Mass coefficients
    double rho_11 = (1.0 - n) * rho_solid_ + (tau - 1.0) * n * rho_f;
    double rho_12 = -n * (tau - 1.0) * rho_f;
    double rho_22 = tau * n * rho_f;
    
    // For high frequency (no fluid flow)
    if (high_frequency_) {
        // Fast P-wave (undrained)
        double lambda_u = K_u - 2.0 * G / 3.0;
        Vp_fast = std::sqrt((lambda_u + 2.0 * G) / rho_b);
        
        // Slow P-wave
        double gamma = n * rho_f / rho_b;
        double M_eff = M * (1.0 - alpha * gamma);
        Vp_slow = std::sqrt(M_eff * n / (tau * rho_f));
        
        // S-wave
        Vs = std::sqrt(G / rho_b);
    }
    else {
        // Low frequency (drained behavior)
        Vp_fast = std::sqrt((K + 4.0 * G / 3.0) / rho_b);
        Vp_slow = 0.0;  // Slow wave is overdamped at low frequency
        Vs = std::sqrt(G / rho_b);
    }
}

void FullBiotDynamics::computeMomentumBalance(const double u[3], const double u_dot[3],
                                              const double u_ddot[3], const double p,
                                              const double p_dot, double residual[3]) const {
    // Full Biot momentum: rho_b * u_ddot + rho_f * w_ddot = div(sigma) + b
    // where w is relative fluid displacement
    
    double rho_b = (1.0 - porosity_) * rho_solid_ + porosity_ * rho_fluid_;
    
    // Simplified residual (without spatial derivatives - those come from FEM)
    for (int i = 0; i < 3; ++i) {
        residual[i] = rho_b * u_ddot[i];  // Inertia term
    }
    
    // Add fluid coupling if using u-w formulation
    if (formulation_ == BiotFormulation::U_W) {
        // Additional terms for relative fluid displacement
    }
}

void FullBiotDynamics::computeMassBalance(const double u[3], const double u_dot[3],
                                          double p, double p_dot,
                                          double& residual) const {
    // Mass balance: (1/M) * p_dot + alpha * div(u_dot) + div(q) = Q
    // where q = -(k/mu) * grad(p)
    
    double alpha = biot_coefficient_;
    double M = biot_modulus_;
    
    // Storage term
    residual = p_dot / M;
    
    // Solid coupling term (div(u_dot) comes from FEM)
    // For point evaluation: residual += alpha * div(u_dot)
}

void FullBiotDynamics::computeCouplingTerms(const double strain[6], double p,
                                            double effective_stress[6],
                                            double& fluid_content) const {
    // Effective stress: sigma'_ij = sigma_ij + alpha * p * delta_ij
    // Fluid content: zeta = alpha * eps_v + p / M
    
    double alpha = biot_coefficient_;
    double M = biot_modulus_;
    
    // Volumetric strain
    double eps_v = strain[0] + strain[1] + strain[2];
    
    // Effective stress (assumes compressive positive for pressure)
    for (int i = 0; i < 6; ++i) {
        effective_stress[i] = -alpha * p * (i < 3 ? 1.0 : 0.0);
    }
    
    // Fluid content change
    fluid_content = alpha * eps_v + p / M;
}

void FullBiotDynamics::assembleSystemMatrix(double M_matrix[], double C_matrix[],
                                            double K_matrix[], int ndof) const {
    // Assemble coupled Biot system matrices
    // [M_uu   M_up] [u_ddot]   [C_uu   C_up] [u_dot]   [K_uu   K_up] [u]   [f]
    // [M_pu   M_pp] [p_ddot] + [C_pu   C_pp] [p_dot] + [K_pu   K_pp] [p] = [q]
    
    // This is a simplified placeholder - actual assembly requires element integration
    double rho_b = (1.0 - porosity_) * rho_solid_ + porosity_ * rho_fluid_;
    double alpha = biot_coefficient_;
    double M_biot = biot_modulus_;
    double k = permeability_;
    double mu = fluid_viscosity_;
    
    // Mass matrix diagonal (simplified)
    for (int i = 0; i < ndof; ++i) {
        M_matrix[i * ndof + i] = rho_b;  // Displacement DOFs
    }
    
    // Damping from Darcy flow
    for (int i = 0; i < ndof; ++i) {
        C_matrix[i * ndof + i] = k / mu;  // Pressure DOFs (hydraulic conductivity)
    }
    
    // Stiffness includes elastic + coupling terms
    // K_uu = stiffness matrix
    // K_up = alpha * B_vol^T (coupling)
    // K_pu = alpha * B_vol (coupling transpose)
    // K_pp = permeability/viscosity
}

// =============================================================================
// THMCoupling Implementation
// =============================================================================

THMCoupling::THMCoupling()
    : coupling_type_(THMCouplingType::ITERATIVE), max_iterations_(20),
      tolerance_(1e-6), relaxation_factor_(0.8),
      thermal_expansion_solid_(1e-5), thermal_expansion_fluid_(2e-4),
      viscosity_temp_coeff_(0.02), storage_compression_coeff_(1e-9),
      stress_permeability_coeff_(5e-8), strain_permeability_exp_(3.0),
      plastic_heating_fraction_(0.9), thermal_pressurization_coeff_(0.1e6) {
}

void THMCoupling::computeThermoMechanicalStress(const double strain[6], double dT,
                                                double thermal_stress[6]) const {
    // Thermal stress: sigma_T = -alpha_T * K * dT * I
    // In tensor form: sigma_T_ij = -3 * alpha * K * dT * delta_ij
    
    double alpha = thermal_expansion_solid_;
    // Assume bulk modulus K (would be passed in practice)
    double K = 36e9;
    
    double sigma_T = -3.0 * alpha * K * dT;
    
    thermal_stress[0] = sigma_T;
    thermal_stress[1] = sigma_T;
    thermal_stress[2] = sigma_T;
    thermal_stress[3] = 0.0;
    thermal_stress[4] = 0.0;
    thermal_stress[5] = 0.0;
}

void THMCoupling::computeThermoHydraulicExpansion(double dT, double& pore_pressure_change,
                                                  double& fluid_density_change) const {
    // Thermal pressurization: dp = Lambda * dT (undrained)
    pore_pressure_change = thermal_pressurization_coeff_ * dT;
    
    // Fluid density change: d(rho)/rho = -beta_f * dT
    fluid_density_change = -thermal_expansion_fluid_ * dT;
}

void THMCoupling::computeStressPermeabilityChange(const double stress[6],
                                                  const double strain[6],
                                                  double k_base,
                                                  double& k_new) const {
    // Stress-dependent permeability models:
    // Option 1: k = k0 * exp(-a * sigma_mean)
    // Option 2: k = k0 * (phi/phi0)^n (Kozeny-Carman)
    
    double sigma_mean = (stress[0] + stress[1] + stress[2]) / 3.0;
    double eps_v = strain[0] + strain[1] + strain[2];
    
    // Model 1: Exponential dependence on mean stress
    double k1 = k_base * std::exp(-stress_permeability_coeff_ * sigma_mean);
    
    // Model 2: Kozeny-Carman (porosity change from strain)
    double phi0 = 0.2;  // Reference porosity
    double phi = phi0 + eps_v * (1.0 - phi0);  // Linear approximation
    phi = std::max(0.01, std::min(0.5, phi));
    double k2 = k_base * std::pow(phi / phi0, strain_permeability_exp_);
    
    // Combine models
    k_new = std::min(k1, k2);
}

void THMCoupling::computePlasticHeating(const double plastic_strain_rate[6],
                                        const double stress[6],
                                        double& heat_source) const {
    // Plastic work: W_p = sigma : d_eps_p
    // Heat source: Q = chi * W_p (Taylor-Quinney coefficient)
    
    double work = 0.0;
    for (int i = 0; i < 3; ++i) {
        work += stress[i] * plastic_strain_rate[i];
    }
    for (int i = 3; i < 6; ++i) {
        work += 2.0 * stress[i] * plastic_strain_rate[i];  // Shear terms
    }
    
    heat_source = plastic_heating_fraction_ * work;
}

void THMCoupling::iterativeCoupling(FieldState& state, const double dt,
                                    int& iterations, bool& converged) const {
    // Iterative staggered THM coupling
    // 1. Solve thermal problem (T^{n+1})
    // 2. Solve hydraulic problem (p^{n+1}) with T^{n+1}
    // 3. Solve mechanical problem (u^{n+1}) with T^{n+1}, p^{n+1}
    // 4. Check convergence, iterate if needed
    
    converged = false;
    
    for (iterations = 0; iterations < max_iterations_; ++iterations) {
        // Store previous iteration values
        double T_prev[3] = {state.temperature[0], state.temperature[1], state.temperature[2]};
        double p_prev[3] = {state.pressure[0], state.pressure[1], state.pressure[2]};
        double u_prev[9];
        for (int i = 0; i < 9; ++i) u_prev[i] = state.displacement[i];
        
        // Thermal solve (placeholder - actual solve would be more complex)
        // Updates state.temperature
        
        // Hydraulic solve with thermal feedback
        // Updates state.pressure
        double dT = state.temperature[0] - T_prev[0];  // Representative
        double dp_thermal, drho;
        computeThermoHydraulicExpansion(dT, dp_thermal, drho);
        for (int i = 0; i < 3; ++i) {
            state.pressure[i] += relaxation_factor_ * dp_thermal;
        }
        
        // Mechanical solve with thermal and hydraulic feedback
        // Updates state.displacement
        
        // Check convergence
        double norm_T = 0.0, norm_p = 0.0, norm_u = 0.0;
        double diff_T = 0.0, diff_p = 0.0, diff_u = 0.0;
        
        for (int i = 0; i < 3; ++i) {
            diff_T += std::pow(state.temperature[i] - T_prev[i], 2);
            norm_T += std::pow(state.temperature[i], 2);
            diff_p += std::pow(state.pressure[i] - p_prev[i], 2);
            norm_p += std::pow(state.pressure[i], 2);
        }
        for (int i = 0; i < 9; ++i) {
            diff_u += std::pow(state.displacement[i] - u_prev[i], 2);
            norm_u += std::pow(state.displacement[i], 2);
        }
        
        double rel_err = std::sqrt(diff_T + diff_p + diff_u) / 
                        (std::sqrt(norm_T + norm_p + norm_u) + 1e-30);
        
        if (rel_err < tolerance_) {
            converged = true;
            break;
        }
    }
}

void THMCoupling::enableFaultThermalPressurization(
    std::shared_ptr<FSRM::ThermalPressurization> tp) {
    // Store the fault TP model in the thermal effects configuration
    // This allows using the detailed 1D diffusion model for fault zones
    // instead of the simpler bulk Lambda coefficient
    fault_tp_model_ = tp;
    use_fault_tp_ = (tp != nullptr);
}

double THMCoupling::getFaultTPContribution(
    FSRM::ThermalPressurization::State& state, double dt) const {
    
    if (!use_fault_tp_ || !fault_tp_model_) {
        return 0.0;
    }
    
    // Store initial pore pressure
    double p_initial = state.pore_pressure;
    
    // Update the ThermalPressurization state using the detailed model
    // This solves the 1D diffusion equations across the fault zone
    fault_tp_model_->update(state, dt);
    
    // Return the pore pressure change from the detailed fault TP model
    return state.pore_pressure - p_initial;
}

// =============================================================================
// THMCCoupling Implementation
// =============================================================================

THMCCoupling::THMCCoupling()
    : reference_temperature_(298.15) {
    reaction_rates_.resize(1, 1e-6);
    activation_energies_.resize(1, 50e3);
}

void THMCCoupling::computeReactionRate(double temperature, double concentration,
                                       int species_idx, double& rate) const {
    // Arrhenius reaction rate: k = A * exp(-E_a / RT) * C
    const double R = 8.314;
    
    double k0 = reaction_rates_[species_idx];
    double Ea = activation_energies_[species_idx];
    
    rate = k0 * std::exp(-Ea / (R * temperature)) * concentration;
}

void THMCCoupling::computePorosityChange(double dissolved_mass, double precipitated_mass,
                                         double mineral_density, double& porosity_change) const {
    // Porosity change from dissolution/precipitation
    // d(phi) = (dissolved - precipitated) / rho_mineral
    
    porosity_change = (dissolved_mass - precipitated_mass) / mineral_density;
}

void THMCCoupling::computePermeabilityFromPorosity(double porosity, double porosity_ref,
                                                   double k_ref, double& k_new) const {
    // Kozeny-Carman relationship
    // k/k0 = (phi/phi0)^3 * ((1-phi0)/(1-phi))^2
    
    double ratio_phi = porosity / porosity_ref;
    double ratio_solid = (1.0 - porosity_ref) / (1.0 - porosity);
    
    k_new = k_ref * std::pow(ratio_phi, 3) * std::pow(ratio_solid, 2);
}

void THMCCoupling::computeChemicalStress(double concentration_change,
                                         double swelling_coefficient,
                                         double stress[6]) const {
    // Chemical swelling stress (similar to thermal expansion)
    // sigma_chem = -K_swell * dC * I
    
    double sigma_c = -swelling_coefficient * concentration_change;
    
    stress[0] = sigma_c;
    stress[1] = sigma_c;
    stress[2] = sigma_c;
    stress[3] = 0.0;
    stress[4] = 0.0;
    stress[5] = 0.0;
}

// =============================================================================
// UnsaturatedFlowCoupling Implementation
// =============================================================================

UnsaturatedFlowCoupling::UnsaturatedFlowCoupling()
    : swrc_model_(SWRCModel::VAN_GENUCHTEN),
      vg_alpha_(0.01), vg_n_(1.5), vg_m_(0.0), residual_saturation_(0.1),
      saturated_permeability_(1e-12), air_entry_pressure_(1000.0) {
    // Compute m from n for VG model
    vg_m_ = 1.0 - 1.0 / vg_n_;
}

double UnsaturatedFlowCoupling::computeSaturation(double suction) const {
    // Suction is positive (psi = p_air - p_water > 0 for unsaturated)
    
    if (suction <= 0.0) {
        return 1.0;  // Saturated
    }
    
    double Se;  // Effective saturation
    
    switch (swrc_model_) {
        case SWRCModel::VAN_GENUCHTEN: {
            // Se = [1 + (alpha * psi)^n]^{-m}
            double alpha_psi = vg_alpha_ * suction;
            Se = std::pow(1.0 + std::pow(alpha_psi, vg_n_), -vg_m_);
            break;
        }
        case SWRCModel::BROOKS_COREY: {
            // Se = (psi_e / psi)^lambda for psi > psi_e, else 1
            if (suction > air_entry_pressure_) {
                double lambda = 2.0;  // Pore size distribution index
                Se = std::pow(air_entry_pressure_ / suction, lambda);
            } else {
                Se = 1.0;
            }
            break;
        }
        default:
            Se = 1.0;
    }
    
    // Convert to actual saturation
    return residual_saturation_ + (1.0 - residual_saturation_) * Se;
}

double UnsaturatedFlowCoupling::computeRelativePermeability(double saturation) const {
    // Effective saturation
    double Se = (saturation - residual_saturation_) / (1.0 - residual_saturation_);
    Se = std::max(0.0, std::min(1.0, Se));
    
    double kr;
    
    switch (swrc_model_) {
        case SWRCModel::VAN_GENUCHTEN: {
            // Mualem model: kr = Se^0.5 * [1 - (1 - Se^{1/m})^m]^2
            double term = 1.0 - std::pow(1.0 - std::pow(Se, 1.0/vg_m_), vg_m_);
            kr = std::sqrt(Se) * term * term;
            break;
        }
        case SWRCModel::BROOKS_COREY: {
            // kr = Se^{(2+3*lambda)/lambda}
            double lambda = 2.0;
            kr = std::pow(Se, (2.0 + 3.0*lambda)/lambda);
            break;
        }
        default:
            kr = Se;
    }
    
    return std::max(1e-10, kr);  // Prevent zero permeability
}

double UnsaturatedFlowCoupling::computeSpecificMoistureCapacity(double suction) const {
    // C = -d(theta)/d(psi) = -n * d(Se)/d(psi)
    // where n is porosity
    
    double n = 0.3;  // Porosity (would be passed in practice)
    double dSe_dpsi;
    
    if (suction <= 0.0) {
        return 0.0;
    }
    
    switch (swrc_model_) {
        case SWRCModel::VAN_GENUCHTEN: {
            double alpha_psi = vg_alpha_ * suction;
            double base = 1.0 + std::pow(alpha_psi, vg_n_);
            double Se = std::pow(base, -vg_m_);
            dSe_dpsi = -vg_m_ * vg_n_ * vg_alpha_ * 
                       std::pow(alpha_psi, vg_n_ - 1.0) * 
                       std::pow(base, -vg_m_ - 1.0);
            break;
        }
        case SWRCModel::BROOKS_COREY: {
            double lambda = 2.0;
            if (suction > air_entry_pressure_) {
                dSe_dpsi = -lambda * std::pow(air_entry_pressure_, lambda) * 
                           std::pow(suction, -lambda - 1.0);
            } else {
                dSe_dpsi = 0.0;
            }
            break;
        }
        default:
            dSe_dpsi = 0.0;
    }
    
    return -n * (1.0 - residual_saturation_) * dSe_dpsi;
}

void UnsaturatedFlowCoupling::computeBishopStress(const double total_stress[6],
                                                  double suction, double saturation,
                                                  double effective_stress[6]) const {
    // Bishop's effective stress for unsaturated soils:
    // sigma' = sigma - u_a + chi*(u_a - u_w)
    // where chi is usually taken as Se
    
    double Se = (saturation - residual_saturation_) / (1.0 - residual_saturation_);
    Se = std::max(0.0, std::min(1.0, Se));
    
    // Assume u_a = 0 (atmospheric), then:
    // sigma'_ij = sigma_ij + chi * psi * delta_ij
    
    for (int i = 0; i < 6; ++i) {
        effective_stress[i] = total_stress[i];
    }
    
    // Add suction contribution to normal stresses
    double suction_stress = Se * suction;
    effective_stress[0] += suction_stress;
    effective_stress[1] += suction_stress;
    effective_stress[2] += suction_stress;
}

void UnsaturatedFlowCoupling::computeCoupledRichardsBiot(const double displacement[3],
                                                        double suction, double saturation,
                                                        double& residual_flow,
                                                        double residual_mech[3]) const {
    // Coupled Richards equation + Biot mechanics
    
    // Flow residual: Storage + Darcy flux
    double C = computeSpecificMoistureCapacity(suction);
    double kr = computeRelativePermeability(saturation);
    
    // Storage term (simplified)
    residual_flow = C;  // Would multiply by d(psi)/dt
    
    // Mechanical coupling
    double alpha_biot = 0.8;
    double Se = (saturation - residual_saturation_) / (1.0 - residual_saturation_);
    
    // Body force from suction (simplified)
    for (int i = 0; i < 3; ++i) {
        residual_mech[i] = alpha_biot * Se * suction;  // Coupling term
    }
}

// =============================================================================
// Factory Functions
// =============================================================================

std::unique_ptr<FullBiotDynamics> createFullBiotSolver(
    const AdvancedCouplingConfig& config) {
    
    auto solver = std::make_unique<FullBiotDynamics>();
    
    if (config.biot_formulation == "u_p") {
        solver->setFormulation(BiotFormulation::U_P);
    } else if (config.biot_formulation == "u_w") {
        solver->setFormulation(BiotFormulation::U_W);
    } else if (config.biot_formulation == "u_p_w") {
        solver->setFormulation(BiotFormulation::U_P_W);
    }
    
    solver->setHighFrequency(config.biot_high_frequency);
    solver->setTortuosity(config.biot_tortuosity);
    
    return solver;
}

std::unique_ptr<THMCoupling> createTHMSolver(const AdvancedCouplingConfig& config) {
    auto solver = std::make_unique<THMCoupling>();
    
    if (config.thm_coupling_type == "sequential") {
        solver->setCouplingType(THMCouplingType::SEQUENTIAL);
    } else if (config.thm_coupling_type == "iterative") {
        solver->setCouplingType(THMCouplingType::ITERATIVE);
    } else if (config.thm_coupling_type == "monolithic") {
        solver->setCouplingType(THMCouplingType::MONOLITHIC);
    }
    
    solver->setMaxIterations(config.max_coupling_iterations);
    solver->setTolerance(config.coupling_tolerance);
    solver->setRelaxationFactor(config.relaxation_factor);
    
    solver->setThermalExpansionSolid(config.thermal_expansion_solid);
    solver->setThermalExpansionFluid(config.thermal_expansion_fluid);
    solver->setThermalPressurizationCoeff(config.thermal_pressurization_coeff);
    
    return solver;
}

std::unique_ptr<UnsaturatedFlowCoupling> createUnsaturatedSolver(
    const AdvancedCouplingConfig& config) {
    
    auto solver = std::make_unique<UnsaturatedFlowCoupling>();
    
    if (config.swrc_model == "van_genuchten") {
        solver->setSWRCModel(SWRCModel::VAN_GENUCHTEN);
    } else if (config.swrc_model == "brooks_corey") {
        solver->setSWRCModel(SWRCModel::BROOKS_COREY);
    }
    
    solver->setVGParameters(config.vg_alpha, config.vg_n);
    solver->setResidualSaturation(config.residual_saturation);
    
    return solver;
}

} // namespace coupled
} // namespace fsrm
