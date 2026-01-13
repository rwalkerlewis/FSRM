/**
 * @file HighFidelityNumerics.cpp
 * @brief Implementation of advanced numerical methods using PETSc
 * 
 * All heavy computational work is delegated to PETSc:
 * - Time integration: PETSc TS
 * - Preconditioning: PETSc PC
 * - Linear algebra: PETSc Vec/Mat
 */

#include "HighFidelityNumerics.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <sstream>

namespace FSRM {
namespace HighFidelity {

// =============================================================================
// Utility Functions (using PETSc types)
// =============================================================================

namespace {

inline PetscReal dot3(const PetscReal a[3], const PetscReal b[3]) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline PetscReal norm3(const PetscReal v[3]) {
    return PetscSqrtReal(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

} // anonymous namespace


// =============================================================================
// SUPGStabilization Implementation
// =============================================================================

SUPGStabilization::SUPGStabilization() {}

void SUPGStabilization::setParameters(const Parameters& params) {
    params_ = params;
}

void SUPGStabilization::configure(const std::map<std::string, std::string>& config) {
    auto it = config.find("supg_method");
    if (it != config.end()) {
        const std::string& val = it->second;
        if (val == "none") params_.method = StabilizationMethod::NONE;
        else if (val == "supg") params_.method = StabilizationMethod::SUPG;
        else if (val == "gls") params_.method = StabilizationMethod::GLS;
        else if (val == "pspg") params_.method = StabilizationMethod::PSPG;
        else if (val == "vms") params_.method = StabilizationMethod::VMS;
        else if (val == "subgrid") params_.method = StabilizationMethod::SUBGRID_VISCOSITY;
    }
    
    it = config.find("supg_tau_formula");
    if (it != config.end()) {
        const std::string& val = it->second;
        if (val == "shakib") params_.tau_formula = TauFormula::SHAKIB;
        else if (val == "tezduyar") params_.tau_formula = TauFormula::TEZDUYAR;
        else if (val == "codina") params_.tau_formula = TauFormula::CODINA;
        else if (val == "franca_valentin") params_.tau_formula = TauFormula::FRANCA_VALENTIN;
        else if (val == "optimal") params_.tau_formula = TauFormula::OPTIMAL;
    }
    
    it = config.find("supg_tau_multiplier");
    if (it != config.end()) params_.tau_multiplier = std::stod(it->second);
    
    it = config.find("supg_crosswind_diffusion");
    if (it != config.end()) params_.crosswind_diffusion = std::stod(it->second);
    
    it = config.find("supg_include_transient");
    if (it != config.end()) params_.include_transient = (it->second == "true" || it->second == "1");
    
    it = config.find("supg_shock_capturing");
    if (it != config.end()) params_.shock_capturing = (it->second == "true" || it->second == "1");
    
    it = config.find("supg_shock_capture_coeff");
    if (it != config.end()) params_.shock_capture_coeff = std::stod(it->second);
}

PetscReal SUPGStabilization::calculateTau(const PetscReal velocity[3],
                                          PetscReal diffusivity, PetscReal reaction,
                                          PetscReal h, PetscReal dt) const {
    PetscReal tau = 0.0;
    
    switch (params_.tau_formula) {
        case TauFormula::SHAKIB:
            tau = tau_Shakib(velocity, diffusivity, reaction, h, dt);
            break;
        case TauFormula::TEZDUYAR:
            tau = tau_Tezduyar(velocity, diffusivity, h, dt);
            break;
        case TauFormula::CODINA:
            tau = tau_Codina(velocity, diffusivity, reaction, h, dt);
            break;
        default:
            tau = tau_Shakib(velocity, diffusivity, reaction, h, dt);
            break;
    }
    
    return params_.tau_multiplier * tau;
}

PetscReal SUPGStabilization::tau_Shakib(const PetscReal velocity[3],
                                        PetscReal diffusivity, PetscReal reaction,
                                        PetscReal h, PetscReal dt) const {
    PetscReal vel_mag = norm3(velocity);
    
    PetscReal term1 = 0.0;
    if (params_.include_transient && dt > 0.0) {
        term1 = 4.0 / (dt * dt);
    }
    
    PetscReal term2 = (h > 0.0) ? PetscSqr(2.0 * vel_mag / h) : 0.0;
    PetscReal term3 = (h > 0.0 && diffusivity > 0.0) ? 
                      PetscSqr(9.0 * 4.0 * diffusivity / (h * h)) : 0.0;
    PetscReal term4 = params_.include_reaction ? reaction * reaction : 0.0;
    
    PetscReal sum = term1 + term2 + term3 + term4;
    return (sum > PETSC_SMALL) ? 1.0 / PetscSqrtReal(sum) : 0.0;
}

PetscReal SUPGStabilization::tau_Tezduyar(const PetscReal velocity[3],
                                          PetscReal diffusivity, PetscReal h, 
                                          PetscReal dt) const {
    PetscReal vel_mag = norm3(velocity);
    
    PetscReal term1 = (params_.include_transient && dt > 0.0) ? PetscSqr(2.0 / dt) : 0.0;
    PetscReal term2 = (h > 0.0) ? PetscSqr(2.0 * vel_mag / h) : 0.0;
    PetscReal term3 = (h > 0.0) ? PetscSqr(4.0 * diffusivity / (h * h)) : 0.0;
    
    PetscReal sum = term1 + term2 + term3;
    return (sum > PETSC_SMALL) ? 1.0 / PetscSqrtReal(sum) : 0.0;
}

PetscReal SUPGStabilization::tau_Codina(const PetscReal velocity[3],
                                        PetscReal diffusivity, PetscReal reaction,
                                        PetscReal h, PetscReal dt) const {
    PetscReal tau = PETSC_MAX_REAL;
    PetscReal vel_mag = norm3(velocity);
    
    if (params_.include_transient && dt > 0.0) {
        tau = PetscMin(tau, dt / 2.0);
    }
    if (vel_mag > PETSC_SMALL && h > 0.0) {
        tau = PetscMin(tau, h / (2.0 * vel_mag));
    }
    if (diffusivity > PETSC_SMALL && h > 0.0) {
        tau = PetscMin(tau, h * h / (4.0 * diffusivity));
    }
    if (params_.include_reaction && reaction > PETSC_SMALL) {
        tau = PetscMin(tau, 1.0 / reaction);
    }
    
    return (tau < PETSC_MAX_REAL / 2.0) ? tau : 0.0;
}

PetscReal SUPGStabilization::pecletNumber(PetscReal velocity_magnitude, PetscReal h,
                                          PetscReal diffusivity) const {
    return (diffusivity > PETSC_SMALL) ? velocity_magnitude * h / (2.0 * diffusivity) : PETSC_MAX_REAL;
}

PetscReal SUPGStabilization::optimalDiffusivity(PetscReal velocity_magnitude, PetscReal h,
                                                 PetscReal diffusivity) const {
    PetscReal Pe = pecletNumber(velocity_magnitude, h, diffusivity);
    if (Pe < 1e-3) return 0.0;
    
    PetscReal coth_Pe = (PetscExpReal(2.0 * Pe) + 1.0) / (PetscExpReal(2.0 * Pe) - 1.0);
    PetscReal alpha = coth_Pe - 1.0 / Pe;
    
    return 0.5 * velocity_magnitude * h * alpha;
}

void SUPGStabilization::addStabilization(const PetscReal velocity[3],
                                         const PetscReal test_gradient[3],
                                         PetscReal residual_strong,
                                         PetscReal tau,
                                         PetscReal* stabilization_term) const {
    if (params_.method == StabilizationMethod::NONE) {
        *stabilization_term = 0.0;
        return;
    }
    
    PetscReal a_grad_v = dot3(velocity, test_gradient);
    *stabilization_term = tau * a_grad_v * residual_strong;
}

void SUPGStabilization::addGLSStabilization(PetscReal L_adjoint_v,
                                            PetscReal residual_strong,
                                            PetscReal tau,
                                            PetscReal* stabilization_term) const {
    *stabilization_term = tau * L_adjoint_v * residual_strong;
}

PetscReal SUPGStabilization::shockCapturingDiffusivity(const PetscReal grad_u[3],
                                                        PetscReal u, PetscReal h,
                                                        PetscReal residual) const {
    if (!params_.shock_capturing) return 0.0;
    
    PetscReal grad_mag = norm3(grad_u);
    if (grad_mag < PETSC_SMALL) return 0.0;
    
    return params_.shock_capture_coeff * h * PetscAbsReal(residual) / grad_mag;
}

PetscReal SUPGStabilization::pspgTerm(const PetscReal grad_p[3],
                                      const PetscReal grad_q[3],
                                      PetscReal tau) const {
    return tau * dot3(grad_p, grad_q);
}


// =============================================================================
// PETScTimeIntegration Implementation
// =============================================================================

PETScTimeIntegration::PETScTimeIntegration(MPI_Comm comm)
    : comm_(comm), ts_(nullptr), ts_created_(PETSC_FALSE) {}

PETScTimeIntegration::~PETScTimeIntegration() {
    if (ts_created_ && ts_) {
        TSDestroy(&ts_);
    }
}

void PETScTimeIntegration::setParameters(const Parameters& params) {
    params_ = params;
}

void PETScTimeIntegration::configure(const std::map<std::string, std::string>& config) {
    auto it = config.find("time_scheme");
    if (it != config.end()) {
        const std::string& val = it->second;
        if (val == "backward_euler" || val == "be") params_.scheme = TimeScheme::BACKWARD_EULER;
        else if (val == "crank_nicolson" || val == "cn") params_.scheme = TimeScheme::CRANK_NICOLSON;
        else if (val == "theta") params_.scheme = TimeScheme::THETA_METHOD;
        else if (val == "rk4") params_.scheme = TimeScheme::RK4;
        else if (val == "ssprk3") params_.scheme = TimeScheme::SSPRK3;
        else if (val == "bdf2") params_.scheme = TimeScheme::BDF2;
        else if (val == "bdf3") params_.scheme = TimeScheme::BDF3;
        else if (val == "bdf4") params_.scheme = TimeScheme::BDF4;
        else if (val == "arkimex2") params_.scheme = TimeScheme::ARKIMEX2;
        else if (val == "arkimex3") params_.scheme = TimeScheme::ARKIMEX3;
        else if (val == "arkimex4") params_.scheme = TimeScheme::ARKIMEX4;
        else if (val == "generalized_alpha") params_.scheme = TimeScheme::GENERALIZED_ALPHA;
        else if (val == "newmark" || val == "newmark_beta") params_.scheme = TimeScheme::NEWMARK_BETA;
    }
    
    it = config.find("time_theta");
    if (it != config.end()) params_.theta = std::stod(it->second);
    
    it = config.find("newmark_beta");
    if (it != config.end()) params_.newmark_beta = std::stod(it->second);
    
    it = config.find("newmark_gamma");
    if (it != config.end()) params_.newmark_gamma = std::stod(it->second);
    
    it = config.find("rho_infinity");
    if (it != config.end()) params_.rho_infinity = std::stod(it->second);
    
    it = config.find("time_adaptive");
    if (it != config.end()) params_.adaptive = (it->second == "true" || it->second == "1") ? PETSC_TRUE : PETSC_FALSE;
    
    it = config.find("time_atol");
    if (it != config.end()) params_.atol = std::stod(it->second);
    
    it = config.find("time_rtol");
    if (it != config.end()) params_.rtol = std::stod(it->second);
}

PetscErrorCode PETScTimeIntegration::createTS(DM dm) {
    PetscFunctionBeginUser;
    
    PetscCall(TSCreate(comm_, &ts_));
    ts_created_ = PETSC_TRUE;
    
    if (dm) {
        PetscCall(TSSetDM(ts_, dm));
    }
    
    PetscCall(configureTSType());
    PetscCall(configureAdaptivity());
    
    // Allow command-line override
    PetscCall(TSSetFromOptions(ts_));
    
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PETScTimeIntegration::configureTSType() {
    PetscFunctionBeginUser;
    
    switch (params_.scheme) {
        case TimeScheme::BACKWARD_EULER:
            PetscCall(TSSetType(ts_, TSTHETA));
            PetscCall(TSThetaSetTheta(ts_, 1.0));
            break;
            
        case TimeScheme::CRANK_NICOLSON:
            PetscCall(TSSetType(ts_, TSTHETA));
            PetscCall(TSThetaSetTheta(ts_, 0.5));
            break;
            
        case TimeScheme::THETA_METHOD:
            PetscCall(TSSetType(ts_, TSTHETA));
            PetscCall(TSThetaSetTheta(ts_, params_.theta));
            break;
            
        case TimeScheme::RK1FE:
            PetscCall(TSSetType(ts_, TSRK));
            PetscCall(TSRKSetType(ts_, TSRK1FE));
            break;
            
        case TimeScheme::RK4:
            PetscCall(TSSetType(ts_, TSRK));
            PetscCall(TSRKSetType(ts_, TSRK4));
            break;
            
        case TimeScheme::RK5DP:
            PetscCall(TSSetType(ts_, TSRK));
            PetscCall(TSRKSetType(ts_, TSRK5DP));
            break;
            
        case TimeScheme::SSPRK3:
            PetscCall(TSSetType(ts_, TSSSP));
            break;
            
        case TimeScheme::ARKIMEX1:
            PetscCall(TSSetType(ts_, TSARKIMEX));
            PetscCall(TSARKIMEXSetType(ts_, TSARKIMEX1BEE));
            break;
            
        case TimeScheme::ARKIMEX2:
            PetscCall(TSSetType(ts_, TSARKIMEX));
            PetscCall(TSARKIMEXSetType(ts_, TSARKIMEXA2));
            break;
            
        case TimeScheme::ARKIMEX3:
            PetscCall(TSSetType(ts_, TSARKIMEX));
            PetscCall(TSARKIMEXSetType(ts_, TSARKIMEX3));
            break;
            
        case TimeScheme::ARKIMEX4:
            PetscCall(TSSetType(ts_, TSARKIMEX));
            PetscCall(TSARKIMEXSetType(ts_, TSARKIMEX4));
            break;
            
        case TimeScheme::ARKIMEX5:
            PetscCall(TSSetType(ts_, TSARKIMEX));
            PetscCall(TSARKIMEXSetType(ts_, TSARKIMEX5));
            break;
            
        case TimeScheme::BDF1:
            PetscCall(TSSetType(ts_, TSBDF));
            PetscCall(TSBDFSetOrder(ts_, 1));
            break;
            
        case TimeScheme::BDF2:
            PetscCall(TSSetType(ts_, TSBDF));
            PetscCall(TSBDFSetOrder(ts_, 2));
            break;
            
        case TimeScheme::BDF3:
            PetscCall(TSSetType(ts_, TSBDF));
            PetscCall(TSBDFSetOrder(ts_, 3));
            break;
            
        case TimeScheme::BDF4:
            PetscCall(TSSetType(ts_, TSBDF));
            PetscCall(TSBDFSetOrder(ts_, 4));
            break;
            
        case TimeScheme::GENERALIZED_ALPHA:
            PetscCall(TSSetType(ts_, TSALPHA));
            PetscCall(TSAlphaSetRadius(ts_, params_.rho_infinity));
            break;
            
        case TimeScheme::GENERALIZED_ALPHA2:
            PetscCall(TSSetType(ts_, TSALPHA2));
            PetscCall(TSAlpha2SetRadius(ts_, params_.rho_infinity));
            break;
            
        case TimeScheme::NEWMARK_BETA:
            // PETSc doesn't have native Newmark, use Alpha2 with specific params
            PetscCall(TSSetType(ts_, TSALPHA2));
            // Set parameters to match Newmark-beta
            // For Newmark: rho_inf determines damping
            PetscCall(TSAlpha2SetRadius(ts_, params_.rho_infinity));
            break;
            
        default:
            PetscCall(TSSetType(ts_, TSBDF));
            PetscCall(TSBDFSetOrder(ts_, 2));
            break;
    }
    
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PETScTimeIntegration::configureAdaptivity() {
    PetscFunctionBeginUser;
    
    TSAdapt adapt;
    PetscCall(TSGetAdapt(ts_, &adapt));
    
    if (params_.adaptive) {
        PetscCall(TSAdaptSetType(adapt, TSADAPTBASIC));
        PetscCall(TSSetTolerances(ts_, params_.atol, NULL, params_.rtol, NULL));
    } else {
        PetscCall(TSAdaptSetType(adapt, TSADAPTNONE));
    }
    
    PetscCall(TSAdaptSetStepLimits(adapt, params_.dt_min, params_.dt_max));
    PetscCall(TSSetMaxSteps(ts_, params_.max_steps));
    
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PETScTimeIntegration::setRHSFunction(TSRHSFunction func, void* ctx) {
    PetscFunctionBeginUser;
    PetscCall(TSSetRHSFunction(ts_, NULL, func, ctx));
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PETScTimeIntegration::setRHSJacobian(TSRHSJacobian func, Mat J, Mat P, void* ctx) {
    PetscFunctionBeginUser;
    PetscCall(TSSetRHSJacobian(ts_, J, P, func, ctx));
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PETScTimeIntegration::setIFunction(TSIFunction func, Vec F, void* ctx) {
    PetscFunctionBeginUser;
    PetscCall(TSSetIFunction(ts_, F, func, ctx));
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PETScTimeIntegration::setIJacobian(TSIJacobian func, Mat J, Mat P, void* ctx) {
    PetscFunctionBeginUser;
    PetscCall(TSSetIJacobian(ts_, J, P, func, ctx));
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PETScTimeIntegration::setI2Function(TSI2Function func, Vec F, void* ctx) {
    PetscFunctionBeginUser;
    PetscCall(TSSetI2Function(ts_, F, func, ctx));
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PETScTimeIntegration::setI2Jacobian(TSI2Jacobian func, Mat J, Mat P, void* ctx) {
    PetscFunctionBeginUser;
    PetscCall(TSSetI2Jacobian(ts_, J, P, func, ctx));
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PETScTimeIntegration::solve(Vec u, PetscReal t0, PetscReal tf) {
    PetscFunctionBeginUser;
    PetscCall(TSSetTime(ts_, t0));
    PetscCall(TSSetMaxTime(ts_, tf));
    PetscCall(TSSetExactFinalTime(ts_, TS_EXACTFINALTIME_MATCHSTEP));
    PetscCall(TSSolve(ts_, u));
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PETScTimeIntegration::step(Vec u) {
    PetscFunctionBeginUser;
    PetscCall(TSStep(ts_));
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PETScTimeIntegration::getTime(PetscReal* t) const {
    PetscFunctionBeginUser;
    PetscCall(TSGetTime(ts_, t));
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PETScTimeIntegration::getTimeStep(PetscReal* dt) const {
    PetscFunctionBeginUser;
    PetscCall(TSGetTimeStep(ts_, dt));
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PETScTimeIntegration::setTimeStep(PetscReal dt) {
    PetscFunctionBeginUser;
    PetscCall(TSSetTimeStep(ts_, dt));
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscInt PETScTimeIntegration::getOrder() const {
    switch (params_.scheme) {
        case TimeScheme::BACKWARD_EULER:
        case TimeScheme::RK1FE:
        case TimeScheme::BDF1:
        case TimeScheme::ARKIMEX1:
            return 1;
        case TimeScheme::CRANK_NICOLSON:
        case TimeScheme::RK2A:
        case TimeScheme::BDF2:
        case TimeScheme::ARKIMEX2:
        case TimeScheme::GENERALIZED_ALPHA:
        case TimeScheme::GENERALIZED_ALPHA2:
        case TimeScheme::NEWMARK_BETA:
            return 2;
        case TimeScheme::RK3:
        case TimeScheme::SSPRK3:
        case TimeScheme::BDF3:
        case TimeScheme::ARKIMEX3:
            return 3;
        case TimeScheme::RK4:
        case TimeScheme::BDF4:
        case TimeScheme::ARKIMEX4:
            return 4;
        case TimeScheme::RK5F:
        case TimeScheme::RK5DP:
        case TimeScheme::BDF5:
        case TimeScheme::ARKIMEX5:
            return 5;
        case TimeScheme::BDF6:
            return 6;
        default:
            return 1;
    }
}

PetscBool PETScTimeIntegration::isAStable() const {
    switch (params_.scheme) {
        case TimeScheme::BACKWARD_EULER:
        case TimeScheme::CRANK_NICOLSON:
        case TimeScheme::THETA_METHOD:
        case TimeScheme::BDF1:
        case TimeScheme::BDF2:
        case TimeScheme::BDF3:
        case TimeScheme::ARKIMEX2:
        case TimeScheme::ARKIMEX3:
        case TimeScheme::ARKIMEX4:
        case TimeScheme::GENERALIZED_ALPHA:
        case TimeScheme::GENERALIZED_ALPHA2:
        case TimeScheme::NEWMARK_BETA:
            return PETSC_TRUE;
        default:
            return PETSC_FALSE;
    }
}

PetscBool PETScTimeIntegration::isLStable() const {
    switch (params_.scheme) {
        case TimeScheme::BACKWARD_EULER:
        case TimeScheme::BDF1:
        case TimeScheme::BDF2:
            return PETSC_TRUE;
        default:
            return PETSC_FALSE;
    }
}


// =============================================================================
// PAdaptivityManager Implementation
// =============================================================================

PAdaptivityManager::PAdaptivityManager() {}

void PAdaptivityManager::setParameters(const Parameters& params) {
    params_ = params;
}

void PAdaptivityManager::configure(const std::map<std::string, std::string>& config) {
    auto it = config.find("p_adapt_min_order");
    if (it != config.end()) params_.min_order = std::stoi(it->second);
    
    it = config.find("p_adapt_max_order");
    if (it != config.end()) params_.max_order = std::stoi(it->second);
    
    it = config.find("p_adapt_default_order");
    if (it != config.end()) params_.default_order = std::stoi(it->second);
    
    it = config.find("p_adapt_smoothness_inc");
    if (it != config.end()) params_.smoothness_threshold_inc = std::stod(it->second);
    
    it = config.find("p_adapt_smoothness_dec");
    if (it != config.end()) params_.smoothness_threshold_dec = std::stod(it->second);
}

PetscErrorCode PAdaptivityManager::smoothnessIndicator(Vec legendre_coeffs, 
                                                        PetscReal* indicator) const {
    PetscFunctionBeginUser;
    
    PetscInt n;
    PetscCall(VecGetSize(legendre_coeffs, &n));
    
    if (n < 2) {
        *indicator = 1.0;
        PetscFunctionReturn(PETSC_SUCCESS);
    }
    
    // Use PETSc Vec operations for efficiency
    const PetscScalar* coeffs;
    PetscCall(VecGetArrayRead(legendre_coeffs, &coeffs));
    
    PetscReal sum_low = 0.0, sum_high = 0.0;
    PetscInt mid = n / 2;
    
    for (PetscInt i = 0; i < mid; ++i) {
        sum_low += PetscRealPart(coeffs[i]) * PetscRealPart(coeffs[i]);
    }
    for (PetscInt i = mid; i < n; ++i) {
        sum_high += PetscRealPart(coeffs[i]) * PetscRealPart(coeffs[i]);
    }
    
    PetscCall(VecRestoreArrayRead(legendre_coeffs, &coeffs));
    
    if (sum_low < PETSC_SMALL) {
        *indicator = 1.0;
    } else {
        PetscReal ratio = sum_high / sum_low;
        *indicator = PetscExpReal(-ratio);
    }
    
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscInt PAdaptivityManager::decideOrderChange(PetscReal smoothness, 
                                                PetscReal error_indicator,
                                                PetscInt current_order) const {
    if (smoothness > params_.smoothness_threshold_inc &&
        error_indicator > params_.p_refinement_threshold &&
        current_order < params_.max_order) {
        return 1;
    }
    
    if ((smoothness < params_.smoothness_threshold_dec ||
         error_indicator < params_.p_refinement_threshold * 0.01) &&
        current_order > params_.min_order) {
        return -1;
    }
    
    return 0;
}

PetscErrorCode PAdaptivityManager::changeOrder(DM dm_old, Vec u_old,
                                                PetscInt new_order, 
                                                DM* dm_new, Vec* u_new) const {
    PetscFunctionBeginUser;
    
    // Clone the DM and set new polynomial order
    PetscCall(DMClone(dm_old, dm_new));
    
    // Get FE space and set new order
    PetscFE fe;
    PetscCall(DMGetField(*dm_new, 0, NULL, (PetscObject*)&fe));
    
    // Create new FE with different order
    PetscInt dim, num_comp;
    PetscCall(PetscFEGetSpatialDimension(fe, &dim));
    PetscCall(PetscFEGetNumComponents(fe, &num_comp));
    
    PetscFE fe_new;
    PetscCall(PetscFECreateDefault(PetscObjectComm((PetscObject)*dm_new), 
                                   dim, num_comp, PETSC_FALSE, NULL, new_order, &fe_new));
    PetscCall(DMSetField(*dm_new, 0, NULL, (PetscObject)fe_new));
    PetscCall(DMCreateDS(*dm_new));
    PetscCall(PetscFEDestroy(&fe_new));
    
    // Create new solution vector
    PetscCall(DMCreateGlobalVector(*dm_new, u_new));
    
    // Project solution from old to new using DMCreateInterpolation
    Mat interp;
    Vec scale;
    PetscCall(DMCreateInterpolation(dm_old, *dm_new, &interp, &scale));
    PetscCall(MatInterpolate(interp, u_old, *u_new));
    PetscCall(MatDestroy(&interp));
    PetscCall(VecDestroy(&scale));
    
    PetscFunctionReturn(PETSC_SUCCESS);
}


// =============================================================================
// PETScPhysicsPreconditioner Implementation
// =============================================================================

PETScPhysicsPreconditioner::PETScPhysicsPreconditioner(MPI_Comm comm)
    : comm_(comm) {}

PETScPhysicsPreconditioner::~PETScPhysicsPreconditioner() {
    for (IS& is : field_is_) {
        if (is) ISDestroy(&is);
    }
}

void PETScPhysicsPreconditioner::setParameters(const Parameters& params) {
    params_ = params;
}

void PETScPhysicsPreconditioner::configure(const std::map<std::string, std::string>& config) {
    auto it = config.find("precond_type");
    if (it != config.end()) {
        const std::string& val = it->second;
        if (val == "bjacobi") params_.type = PreconditionerType::BJACOBI;
        else if (val == "asm") params_.type = PreconditionerType::ASM;
        else if (val == "fieldsplit_additive") params_.type = PreconditionerType::FIELDSPLIT_ADDITIVE;
        else if (val == "fieldsplit_multiplicative") params_.type = PreconditionerType::FIELDSPLIT_MULTIPLICATIVE;
        else if (val == "fieldsplit_schur") params_.type = PreconditionerType::FIELDSPLIT_SCHUR;
        else if (val == "mg") params_.type = PreconditionerType::MG;
        else if (val == "gamg") params_.type = PreconditionerType::GAMG;
        else if (val == "hypre") params_.type = PreconditionerType::HYPRE_BOOMERAMG;
    }
    
    it = config.find("precond_mg_levels");
    if (it != config.end()) params_.mg_levels = std::stoi(it->second);
    
    it = config.find("precond_mg_smoother");
    if (it != config.end()) params_.mg_smoother = it->second;
}

PetscErrorCode PETScPhysicsPreconditioner::setup(KSP ksp, Mat A) {
    PetscFunctionBeginUser;
    
    PC pc;
    PetscCall(KSPGetPC(ksp, &pc));
    
    PetscCall(setPCType(pc));
    
    switch (params_.type) {
        case PreconditionerType::FIELDSPLIT_ADDITIVE:
        case PreconditionerType::FIELDSPLIT_MULTIPLICATIVE:
        case PreconditionerType::FIELDSPLIT_SYMMETRIC:
        case PreconditionerType::FIELDSPLIT_SCHUR:
            PetscCall(setupFieldSplit(pc));
            break;
        case PreconditionerType::MG:
        case PreconditionerType::GAMG:
            PetscCall(setupMultigrid(pc));
            break;
        default:
            break;
    }
    
    PetscCall(PCSetFromOptions(pc));
    
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PETScPhysicsPreconditioner::setPCType(PC pc) {
    PetscFunctionBeginUser;
    
    switch (params_.type) {
        case PreconditionerType::BJACOBI:
            PetscCall(PCSetType(pc, PCBJACOBI));
            break;
        case PreconditionerType::ASM:
            PetscCall(PCSetType(pc, PCASM));
            break;
        case PreconditionerType::GASM:
            PetscCall(PCSetType(pc, PCGASM));
            break;
        case PreconditionerType::FIELDSPLIT_ADDITIVE:
        case PreconditionerType::FIELDSPLIT_MULTIPLICATIVE:
        case PreconditionerType::FIELDSPLIT_SYMMETRIC:
        case PreconditionerType::FIELDSPLIT_SCHUR:
            PetscCall(PCSetType(pc, PCFIELDSPLIT));
            break;
        case PreconditionerType::MG:
            PetscCall(PCSetType(pc, PCMG));
            break;
        case PreconditionerType::GAMG:
            PetscCall(PCSetType(pc, PCGAMG));
            break;
        case PreconditionerType::HYPRE_BOOMERAMG:
            PetscCall(PCSetType(pc, PCHYPRE));
            PetscCall(PCHYPRESetType(pc, "boomeramg"));
            break;
        case PreconditionerType::LU:
            PetscCall(PCSetType(pc, PCLU));
            break;
        case PreconditionerType::CHOLESKY:
            PetscCall(PCSetType(pc, PCCHOLESKY));
            break;
        default:
            PetscCall(PCSetType(pc, PCILU));
            break;
    }
    
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PETScPhysicsPreconditioner::addField(const std::string& name, IS is,
                                                     PetscInt num_components,
                                                     const std::string& sub_pc) {
    PetscFunctionBeginUser;
    
    FieldInfo info;
    info.name = name;
    info.index_set = is;
    info.num_components = num_components;
    info.sub_pc_type = sub_pc;
    
    params_.fields.push_back(info);
    field_is_.push_back(is);
    
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PETScPhysicsPreconditioner::createFieldISFromDM(DM dm, 
                                        const std::vector<std::string>& field_names) {
    PetscFunctionBeginUser;
    
    PetscInt num_fields;
    IS *field_is_array = NULL;
    PetscCall(DMCreateFieldIS(dm, &num_fields, NULL, &field_is_array));
    
    // Copy to our vector
    for (PetscInt i = 0; i < num_fields && static_cast<size_t>(i) < field_names.size(); ++i) {
        field_is_.push_back(field_is_array[i]);
    }
    
    // Free the array (but not the IS objects themselves)
    PetscCall(PetscFree(field_is_array));
    
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PETScPhysicsPreconditioner::setupFieldSplit(PC pc) {
    PetscFunctionBeginUser;
    
    // Set split type
    switch (params_.type) {
        case PreconditionerType::FIELDSPLIT_ADDITIVE:
            PetscCall(PCFieldSplitSetType(pc, PC_COMPOSITE_ADDITIVE));
            break;
        case PreconditionerType::FIELDSPLIT_MULTIPLICATIVE:
            PetscCall(PCFieldSplitSetType(pc, PC_COMPOSITE_MULTIPLICATIVE));
            break;
        case PreconditionerType::FIELDSPLIT_SYMMETRIC:
            PetscCall(PCFieldSplitSetType(pc, PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE));
            break;
        case PreconditionerType::FIELDSPLIT_SCHUR:
            PetscCall(PCFieldSplitSetType(pc, PC_COMPOSITE_SCHUR));
            PetscCall(setupSchurComplement(pc));
            break;
        default:
            break;
    }
    
    // Add fields
    for (size_t i = 0; i < params_.fields.size(); ++i) {
        const FieldInfo& field = params_.fields[i];
        PetscCall(PCFieldSplitSetIS(pc, field.name.c_str(), field.index_set));
    }
    
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PETScPhysicsPreconditioner::setupSchurComplement(PC pc) {
    PetscFunctionBeginUser;
    
    // Set Schur factorization type
    PCFieldSplitSchurFactType fact_type;
    switch (params_.schur_fact) {
        case SchurFactType::DIAG:
            fact_type = PC_FIELDSPLIT_SCHUR_FACT_DIAG;
            break;
        case SchurFactType::LOWER:
            fact_type = PC_FIELDSPLIT_SCHUR_FACT_LOWER;
            break;
        case SchurFactType::UPPER:
            fact_type = PC_FIELDSPLIT_SCHUR_FACT_UPPER;
            break;
        case SchurFactType::FULL:
        default:
            fact_type = PC_FIELDSPLIT_SCHUR_FACT_FULL;
            break;
    }
    PetscCall(PCFieldSplitSetSchurFactType(pc, fact_type));
    
    // Set Schur preconditioner type
    PCFieldSplitSchurPreType pre_type;
    switch (params_.schur_pre) {
        case SchurPreType::SELF:
            pre_type = PC_FIELDSPLIT_SCHUR_PRE_SELF;
            break;
        case SchurPreType::SELFP:
            pre_type = PC_FIELDSPLIT_SCHUR_PRE_SELFP;
            break;
        case SchurPreType::A11:
            pre_type = PC_FIELDSPLIT_SCHUR_PRE_A11;
            break;
        case SchurPreType::USER:
            pre_type = PC_FIELDSPLIT_SCHUR_PRE_USER;
            break;
        case SchurPreType::FULL:
        default:
            pre_type = PC_FIELDSPLIT_SCHUR_PRE_FULL;
            break;
    }
    PetscCall(PCFieldSplitSetSchurPre(pc, pre_type, NULL));
    
    // Set Schur scale
    PetscCall(PCFieldSplitSetSchurScale(pc, params_.schur_scale));
    
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PETScPhysicsPreconditioner::setupMultigrid(PC pc) {
    PetscFunctionBeginUser;
    
    if (params_.type == PreconditionerType::MG) {
        PetscCall(PCMGSetLevels(pc, params_.mg_levels, NULL));
        PetscCall(PCMGSetCycleType(pc, PC_MG_CYCLE_V));
        PetscCall(PCMGSetNumberSmooth(pc, params_.mg_smooth_down));
    } else if (params_.type == PreconditionerType::GAMG) {
        PetscCall(PCGAMGSetNSmooths(pc, params_.gamg_agg_nsmooths));
        PetscCall(PCGAMGSetThreshold(pc, &params_.gamg_threshold, 1));
    }
    
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PETScPhysicsPreconditioner::getPC(KSP ksp, PC* pc) {
    PetscFunctionBeginUser;
    PetscCall(KSPGetPC(ksp, pc));
    PetscFunctionReturn(PETSC_SUCCESS);
}


// =============================================================================
// HighFidelityNumericsConfig Implementation
// =============================================================================

void HighFidelityNumericsConfig::parseConfig(const std::map<std::string, std::string>& config) {
    auto it = config.find("enable_high_fidelity_numerics");
    if (it != config.end()) {
        enable_high_fidelity_numerics = (it->second == "true" || it->second == "1");
    }
    
    it = config.find("enable_supg");
    if (it != config.end()) {
        enable_supg = (it->second == "true" || it->second == "1");
    }
    
    it = config.find("enable_advanced_time");
    if (it != config.end()) {
        enable_advanced_time = (it->second == "true" || it->second == "1");
    }
    
    it = config.find("enable_p_adapt");
    if (it != config.end()) {
        enable_p_adapt = (it->second == "true" || it->second == "1");
    }
    
    it = config.find("enable_physics_precond");
    if (it != config.end()) {
        enable_physics_precond = (it->second == "true" || it->second == "1");
    }
}

PetscErrorCode HighFidelityNumericsConfig::validate(std::string& error_msg) const {
    PetscFunctionBeginUser;
    error_msg.clear();
    
    if (enable_advanced_time) {
        if (time_params.theta < 0.0 || time_params.theta > 1.0) {
            error_msg = "Time integration theta must be in [0, 1]";
            PetscFunctionReturn(PETSC_ERR_ARG_OUTOFRANGE);
        }
        if (time_params.rho_infinity < 0.0 || time_params.rho_infinity > 1.0) {
            error_msg = "Generalized-alpha rho_infinity must be in [0, 1]";
            PetscFunctionReturn(PETSC_ERR_ARG_OUTOFRANGE);
        }
    }
    
    if (enable_p_adapt) {
        if (p_adapt_params.min_order < 1) {
            error_msg = "p-adaptivity min_order must be >= 1";
            PetscFunctionReturn(PETSC_ERR_ARG_OUTOFRANGE);
        }
        if (p_adapt_params.max_order < p_adapt_params.min_order) {
            error_msg = "p-adaptivity max_order must be >= min_order";
            PetscFunctionReturn(PETSC_ERR_ARG_OUTOFRANGE);
        }
    }
    
    PetscFunctionReturn(PETSC_SUCCESS);
}


// =============================================================================
// Factory Functions
// =============================================================================

std::unique_ptr<SUPGStabilization> createSUPGStabilization(
    const std::map<std::string, std::string>& config) {
    
    auto supg = std::make_unique<SUPGStabilization>();
    supg->configure(config);
    return supg;
}

std::unique_ptr<PETScTimeIntegration> createTimeIntegration(
    MPI_Comm comm, const std::map<std::string, std::string>& config) {
    
    auto ts = std::make_unique<PETScTimeIntegration>(comm);
    ts->configure(config);
    return ts;
}

std::unique_ptr<PETScPhysicsPreconditioner> createPhysicsPreconditioner(
    MPI_Comm comm, const std::map<std::string, std::string>& config) {
    
    auto pc = std::make_unique<PETScPhysicsPreconditioner>(comm);
    pc->configure(config);
    return pc;
}

} // namespace HighFidelity
} // namespace FSRM
