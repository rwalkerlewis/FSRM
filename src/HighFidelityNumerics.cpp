/**
 * @file HighFidelityNumerics.cpp
 * @brief Implementation of advanced numerical methods for high-fidelity simulations
 * 
 * Note: DG methods are implemented in DiscontinuousGalerkin.cpp
 *       AMR methods are implemented in AdaptiveMeshRefinement.cpp
 */

#include "HighFidelityNumerics.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <sstream>

namespace FSRM {
namespace HighFidelity {

// =============================================================================
// Utility Functions
// =============================================================================

namespace {

inline double sign(double x) {
    if (x > 0.0) return 1.0;
    if (x < 0.0) return -1.0;
    return 0.0;
}

inline double dot3(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline double norm3(const std::array<double, 3>& v) {
    return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

inline double clamp(double val, double min_val, double max_val) {
    return std::max(min_val, std::min(max_val, val));
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
    if (it != config.end()) {
        params_.tau_multiplier = std::stod(it->second);
    }
    
    it = config.find("supg_crosswind_diffusion");
    if (it != config.end()) {
        params_.crosswind_diffusion = std::stod(it->second);
    }
    
    it = config.find("supg_include_transient");
    if (it != config.end()) {
        params_.include_transient = (it->second == "true" || it->second == "1");
    }
    
    it = config.find("supg_time_scale");
    if (it != config.end()) {
        params_.time_scale = std::stod(it->second);
    }
    
    it = config.find("supg_shock_capturing");
    if (it != config.end()) {
        params_.shock_capturing = (it->second == "true" || it->second == "1");
    }
    
    it = config.find("supg_shock_capture_coeff");
    if (it != config.end()) {
        params_.shock_capture_coeff = std::stod(it->second);
    }
}

double SUPGStabilization::calculateTau(const std::array<double, 3>& velocity,
                                       double diffusivity, double reaction,
                                       double h, double dt) const {
    double tau = 0.0;
    
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
        case TauFormula::FRANCA_VALENTIN:
        case TauFormula::OPTIMAL:
        default:
            tau = tau_Shakib(velocity, diffusivity, reaction, h, dt);
            break;
    }
    
    return params_.tau_multiplier * tau;
}

double SUPGStabilization::tau_Shakib(const std::array<double, 3>& velocity,
                                     double diffusivity, double reaction,
                                     double h, double dt) const {
    // Shakib formula (1991):
    // τ = (4/dt² + (2|a|/h)² + 9(4κ/h²)² + σ²)^(-1/2)
    
    double vel_mag = norm3(velocity);
    
    double term1 = 0.0;
    if (params_.include_transient && dt > 0.0) {
        term1 = 4.0 / (dt * dt);
    }
    
    double term2 = 0.0;
    if (h > 0.0) {
        term2 = std::pow(2.0 * vel_mag / h, 2);
    }
    
    double term3 = 0.0;
    if (h > 0.0 && diffusivity > 0.0) {
        term3 = std::pow(9.0 * 4.0 * diffusivity / (h * h), 2);
    }
    
    double term4 = 0.0;
    if (params_.include_reaction) {
        term4 = reaction * reaction;
    }
    
    double sum = term1 + term2 + term3 + term4;
    if (sum < 1e-30) {
        return 0.0;
    }
    
    return 1.0 / std::sqrt(sum);
}

double SUPGStabilization::tau_Tezduyar(const std::array<double, 3>& velocity,
                                       double diffusivity, double h, double dt) const {
    // Tezduyar formula for incompressible flow:
    // τ = ((2/dt)² + (2|a|/h)² + (4ν/h²)²)^(-1/2)
    
    double vel_mag = norm3(velocity);
    
    double term1 = 0.0;
    if (params_.include_transient && dt > 0.0) {
        term1 = std::pow(2.0 / dt, 2);
    }
    
    double term2 = 0.0;
    if (h > 0.0) {
        term2 = std::pow(2.0 * vel_mag / h, 2);
    }
    
    double term3 = 0.0;
    if (h > 0.0) {
        term3 = std::pow(4.0 * diffusivity / (h * h), 2);
    }
    
    double sum = term1 + term2 + term3;
    if (sum < 1e-30) {
        return 0.0;
    }
    
    return 1.0 / std::sqrt(sum);
}

double SUPGStabilization::tau_Codina(const std::array<double, 3>& velocity,
                                     double diffusivity, double reaction,
                                     double h, double dt) const {
    // Codina (projection-based):
    // τ = min(c1·dt/2, c2·h/(2|a|), c3·h²/(4κ), c4/σ)
    
    const double c1 = 1.0;
    const double c2 = 1.0;
    const double c3 = 1.0;
    const double c4 = 1.0;
    
    double tau = 1e10;
    
    if (params_.include_transient && dt > 0.0) {
        tau = std::min(tau, c1 * dt / 2.0);
    }
    
    double vel_mag = norm3(velocity);
    if (vel_mag > 1e-12 && h > 0.0) {
        tau = std::min(tau, c2 * h / (2.0 * vel_mag));
    }
    
    if (diffusivity > 1e-12 && h > 0.0) {
        tau = std::min(tau, c3 * h * h / (4.0 * diffusivity));
    }
    
    if (params_.include_reaction && reaction > 1e-12) {
        tau = std::min(tau, c4 / reaction);
    }
    
    if (tau > 1e9) {
        tau = 0.0;
    }
    
    return tau;
}

double SUPGStabilization::pecletNumber(double velocity_magnitude, double h,
                                       double diffusivity) const {
    if (diffusivity < 1e-30) {
        return 1e10;  // Pure advection
    }
    return velocity_magnitude * h / (2.0 * diffusivity);
}

double SUPGStabilization::optimalDiffusivity(double velocity_magnitude, double h,
                                              double diffusivity) const {
    double Pe = pecletNumber(velocity_magnitude, h, diffusivity);
    
    if (Pe < 1e-3) {
        return 0.0;
    }
    
    double coth_Pe = (std::exp(2.0 * Pe) + 1.0) / (std::exp(2.0 * Pe) - 1.0);
    double alpha = coth_Pe - 1.0 / Pe;
    
    return 0.5 * velocity_magnitude * h * alpha;
}

void SUPGStabilization::addStabilization(const std::array<double, 3>& velocity,
                                         const std::array<double, 3>& test_gradient,
                                         double residual_strong,
                                         double tau,
                                         double& stabilization_term) const {
    if (params_.method == StabilizationMethod::NONE) {
        stabilization_term = 0.0;
        return;
    }
    
    double a_grad_v = dot3(velocity, test_gradient);
    stabilization_term = tau * a_grad_v * residual_strong;
}

void SUPGStabilization::addGLSStabilization(double L_adjoint_v,
                                            double residual_strong,
                                            double tau,
                                            double& stabilization_term) const {
    stabilization_term = tau * L_adjoint_v * residual_strong;
}

double SUPGStabilization::shockCapturingDiffusivity(const std::array<double, 3>& grad_u,
                                                     double u, double h,
                                                     double residual) const {
    if (!params_.shock_capturing) {
        return 0.0;
    }
    
    double grad_mag = norm3(grad_u);
    if (grad_mag < 1e-30) {
        return 0.0;
    }
    
    double kappa_DC = params_.shock_capture_coeff * h * std::abs(residual) / grad_mag;
    
    return kappa_DC;
}

double SUPGStabilization::pspgTerm(const std::array<double, 3>& grad_p,
                                   const std::array<double, 3>& grad_q,
                                   double tau) const {
    return tau * dot3(grad_p, grad_q);
}


// =============================================================================
// AdvancedTimeIntegration Implementation
// =============================================================================

AdvancedTimeIntegration::AdvancedTimeIntegration() {
    params_.newmark_beta = 0.25;
    params_.newmark_gamma = 0.5;
}

void AdvancedTimeIntegration::setParameters(const Parameters& params) {
    params_ = params;
    
    if (params_.scheme == TimeScheme::GENERALIZED_ALPHA) {
        computeGeneralizedAlphaCoeffs();
    }
}

void AdvancedTimeIntegration::configure(const std::map<std::string, std::string>& config) {
    auto it = config.find("time_scheme");
    if (it != config.end()) {
        const std::string& val = it->second;
        if (val == "forward_euler") params_.scheme = TimeScheme::FORWARD_EULER;
        else if (val == "backward_euler") params_.scheme = TimeScheme::BACKWARD_EULER;
        else if (val == "crank_nicolson") params_.scheme = TimeScheme::CRANK_NICOLSON;
        else if (val == "theta") params_.scheme = TimeScheme::THETA_METHOD;
        else if (val == "rk2") params_.scheme = TimeScheme::RK2;
        else if (val == "rk4") params_.scheme = TimeScheme::RK4;
        else if (val == "ssprk3") params_.scheme = TimeScheme::SSPRK3;
        else if (val == "dirk2") params_.scheme = TimeScheme::DIRK2;
        else if (val == "dirk3") params_.scheme = TimeScheme::DIRK3;
        else if (val == "esdirk4") params_.scheme = TimeScheme::ESDIRK4;
        else if (val == "bdf2") params_.scheme = TimeScheme::BDF2;
        else if (val == "bdf3") params_.scheme = TimeScheme::BDF3;
        else if (val == "bdf4") params_.scheme = TimeScheme::BDF4;
        else if (val == "newmark" || val == "newmark_beta") params_.scheme = TimeScheme::NEWMARK_BETA;
        else if (val == "hht" || val == "hht_alpha") params_.scheme = TimeScheme::HHT_ALPHA;
        else if (val == "generalized_alpha") params_.scheme = TimeScheme::GENERALIZED_ALPHA;
        else if (val == "bathe") params_.scheme = TimeScheme::BATHE;
    }
    
    it = config.find("time_theta");
    if (it != config.end()) {
        params_.theta = std::stod(it->second);
    }
    
    it = config.find("newmark_beta");
    if (it != config.end()) {
        params_.newmark_beta = std::stod(it->second);
    }
    
    it = config.find("newmark_gamma");
    if (it != config.end()) {
        params_.newmark_gamma = std::stod(it->second);
    }
    
    it = config.find("rho_infinity");
    if (it != config.end()) {
        params_.rho_infinity = std::stod(it->second);
    }
    
    it = config.find("time_adaptive");
    if (it != config.end()) {
        params_.adaptive = (it->second == "true" || it->second == "1");
    }
    
    it = config.find("time_error_tolerance");
    if (it != config.end()) {
        params_.error_tolerance = std::stod(it->second);
    }
    
    if (params_.scheme == TimeScheme::GENERALIZED_ALPHA) {
        computeGeneralizedAlphaCoeffs();
    }
}

int AdvancedTimeIntegration::getOrder() const {
    switch (params_.scheme) {
        case TimeScheme::FORWARD_EULER:
        case TimeScheme::BACKWARD_EULER:
            return 1;
        case TimeScheme::CRANK_NICOLSON:
        case TimeScheme::RK2:
        case TimeScheme::DIRK2:
        case TimeScheme::BDF2:
        case TimeScheme::NEWMARK_BETA:
        case TimeScheme::HHT_ALPHA:
        case TimeScheme::GENERALIZED_ALPHA:
            return 2;
        case TimeScheme::SSPRK3:
        case TimeScheme::DIRK3:
        case TimeScheme::BDF3:
            return 3;
        case TimeScheme::RK4:
        case TimeScheme::ESDIRK4:
        case TimeScheme::BDF4:
        case TimeScheme::BATHE:
            return 4;
        case TimeScheme::THETA_METHOD:
            return (std::abs(params_.theta - 0.5) < 1e-10) ? 2 : 1;
        default:
            return 1;
    }
}

bool AdvancedTimeIntegration::isAStable() const {
    switch (params_.scheme) {
        case TimeScheme::FORWARD_EULER:
        case TimeScheme::RK2:
        case TimeScheme::RK4:
        case TimeScheme::SSPRK3:
        case TimeScheme::ADAMS_BASHFORTH_2:
        case TimeScheme::BDF4:
            return false;
        case TimeScheme::BACKWARD_EULER:
        case TimeScheme::CRANK_NICOLSON:
        case TimeScheme::DIRK2:
        case TimeScheme::DIRK3:
        case TimeScheme::ESDIRK4:
        case TimeScheme::BDF2:
        case TimeScheme::BDF3:
        case TimeScheme::ADAMS_MOULTON_2:
        case TimeScheme::NEWMARK_BETA:
        case TimeScheme::HHT_ALPHA:
        case TimeScheme::GENERALIZED_ALPHA:
        case TimeScheme::BATHE:
            return true;
        case TimeScheme::THETA_METHOD:
            return params_.theta >= 0.5;
        default:
            return false;
    }
}

bool AdvancedTimeIntegration::isLStable() const {
    switch (params_.scheme) {
        case TimeScheme::BACKWARD_EULER:
        case TimeScheme::DIRK2:
        case TimeScheme::DIRK3:
        case TimeScheme::ESDIRK4:
        case TimeScheme::BDF2:
        case TimeScheme::BDF3:
            return true;
        case TimeScheme::CRANK_NICOLSON:
        case TimeScheme::ADAMS_MOULTON_2:
            return false;
        case TimeScheme::THETA_METHOD:
            return params_.theta > 0.5;
        case TimeScheme::NEWMARK_BETA:
        case TimeScheme::HHT_ALPHA:
        case TimeScheme::GENERALIZED_ALPHA:
            return params_.newmark_gamma >= 0.5 + params_.newmark_beta;
        default:
            return false;
    }
}

void AdvancedTimeIntegration::getButcherTableau(std::vector<std::vector<double>>& A,
                                                std::vector<double>& b,
                                                std::vector<double>& c) const {
    switch (params_.scheme) {
        case TimeScheme::FORWARD_EULER:
            A = {{0.0}};
            b = {1.0};
            c = {0.0};
            break;
        case TimeScheme::BACKWARD_EULER:
            A = {{1.0}};
            b = {1.0};
            c = {1.0};
            break;
        case TimeScheme::RK2:
            A = {{0.0, 0.0}, {0.5, 0.0}};
            b = {0.0, 1.0};
            c = {0.0, 0.5};
            break;
        case TimeScheme::SSPRK3:
            A = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.25, 0.25, 0.0}};
            b = {1.0/6.0, 1.0/6.0, 2.0/3.0};
            c = {0.0, 1.0, 0.5};
            break;
        case TimeScheme::RK4:
            A = {{0.0, 0.0, 0.0, 0.0},
                 {0.5, 0.0, 0.0, 0.0},
                 {0.0, 0.5, 0.0, 0.0},
                 {0.0, 0.0, 1.0, 0.0}};
            b = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};
            c = {0.0, 0.5, 0.5, 1.0};
            break;
        case TimeScheme::DIRK2: {
            double gamma = 1.0 - 1.0/std::sqrt(2.0);
            A = {{gamma, 0.0}, {1.0 - gamma, gamma}};
            b = {1.0 - gamma, gamma};
            c = {gamma, 1.0};
            break;
        }
        default:
            A = {{0.0}};
            b = {1.0};
            c = {0.0};
    }
}

void AdvancedTimeIntegration::initializeNewmark(const std::vector<double>& u0,
                                                const std::vector<double>& v0,
                                                const std::function<std::vector<double>(
                                                    const std::vector<double>&,
                                                    const std::vector<double>&)>& compute_rhs,
                                                std::vector<double>& a0) const {
    a0 = compute_rhs(u0, v0);
}

void AdvancedTimeIntegration::newmarkPredict(const std::vector<double>& u_n,
                                             const std::vector<double>& v_n,
                                             const std::vector<double>& a_n,
                                             double dt,
                                             std::vector<double>& u_pred,
                                             std::vector<double>& v_pred) const {
    double beta = params_.newmark_beta;
    double gamma = params_.newmark_gamma;
    
    size_t n = u_n.size();
    u_pred.resize(n);
    v_pred.resize(n);
    
    for (size_t i = 0; i < n; ++i) {
        u_pred[i] = u_n[i] + dt * v_n[i] + dt * dt * (0.5 - beta) * a_n[i];
        v_pred[i] = v_n[i] + dt * (1.0 - gamma) * a_n[i];
    }
}

void AdvancedTimeIntegration::newmarkCorrect(const std::vector<double>& u_pred,
                                             const std::vector<double>& v_pred,
                                             const std::vector<double>& a_new,
                                             double dt,
                                             std::vector<double>& u_new,
                                             std::vector<double>& v_new) const {
    double beta = params_.newmark_beta;
    double gamma = params_.newmark_gamma;
    
    size_t n = u_pred.size();
    u_new.resize(n);
    v_new.resize(n);
    
    for (size_t i = 0; i < n; ++i) {
        u_new[i] = u_pred[i] + beta * dt * dt * a_new[i];
        v_new[i] = v_pred[i] + gamma * dt * a_new[i];
    }
}

void AdvancedTimeIntegration::effectiveNewmarkMatrix(double dt, double mass, double damping,
                                                     double stiffness, double& M_eff) const {
    double beta = params_.newmark_beta;
    double gamma = params_.newmark_gamma;
    
    M_eff = mass + gamma * dt * damping + beta * dt * dt * stiffness;
}

void AdvancedTimeIntegration::computeGeneralizedAlphaCoeffs() {
    double rho = params_.rho_infinity;
    
    params_.alpha_m = (2.0 * rho - 1.0) / (rho + 1.0);
    params_.alpha_f = rho / (rho + 1.0);
    params_.newmark_gamma = 0.5 - params_.alpha_m + params_.alpha_f;
    params_.newmark_beta = 0.25 * std::pow(1.0 - params_.alpha_m + params_.alpha_f, 2);
}

double AdvancedTimeIntegration::estimateError(const std::vector<double>& u_high,
                                              const std::vector<double>& u_low) const {
    if (u_high.size() != u_low.size() || u_high.empty()) {
        return 0.0;
    }
    
    double error_sq = 0.0;
    double norm_sq = 0.0;
    
    for (size_t i = 0; i < u_high.size(); ++i) {
        double diff = u_high[i] - u_low[i];
        error_sq += diff * diff;
        norm_sq += u_high[i] * u_high[i];
    }
    
    if (norm_sq < 1e-30) {
        return std::sqrt(error_sq);
    }
    
    return std::sqrt(error_sq / norm_sq);
}

double AdvancedTimeIntegration::computeNewTimestep(double dt_current, double error,
                                                   double tolerance) const {
    if (error < 1e-30) {
        return dt_current * params_.max_dt_factor;
    }
    
    int order = getOrder();
    double factor = params_.safety_factor * std::pow(tolerance / error, 1.0 / (order + 1));
    
    factor = clamp(factor, params_.min_dt_factor, params_.max_dt_factor);
    
    return dt_current * factor;
}

void AdvancedTimeIntegration::bdfCoefficients(int order, std::vector<double>& alpha,
                                              double& beta) const {
    alpha.clear();
    
    switch (order) {
        case 1:
            alpha = {1.0, -1.0};
            beta = 1.0;
            break;
        case 2:
            alpha = {3.0/2.0, -2.0, 1.0/2.0};
            beta = 1.0;
            break;
        case 3:
            alpha = {11.0/6.0, -3.0, 3.0/2.0, -1.0/3.0};
            beta = 1.0;
            break;
        case 4:
            alpha = {25.0/12.0, -4.0, 3.0, -4.0/3.0, 1.0/4.0};
            beta = 1.0;
            break;
        default:
            alpha = {1.0, -1.0};
            beta = 1.0;
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
    if (it != config.end()) {
        params_.min_order = std::stoi(it->second);
    }
    
    it = config.find("p_adapt_max_order");
    if (it != config.end()) {
        params_.max_order = std::stoi(it->second);
    }
    
    it = config.find("p_adapt_default_order");
    if (it != config.end()) {
        params_.default_order = std::stoi(it->second);
    }
    
    it = config.find("p_adapt_smoothness_inc");
    if (it != config.end()) {
        params_.smoothness_threshold_inc = std::stod(it->second);
    }
    
    it = config.find("p_adapt_smoothness_dec");
    if (it != config.end()) {
        params_.smoothness_threshold_dec = std::stod(it->second);
    }
}

double PAdaptivityManager::smoothnessIndicator(const std::vector<double>& legendre_coeffs) const {
    if (legendre_coeffs.size() < 2) {
        return 1.0;
    }
    
    double sum_high = 0.0;
    double sum_low = 0.0;
    
    size_t mid = legendre_coeffs.size() / 2;
    
    for (size_t i = 0; i < mid; ++i) {
        sum_low += legendre_coeffs[i] * legendre_coeffs[i];
    }
    for (size_t i = mid; i < legendre_coeffs.size(); ++i) {
        sum_high += legendre_coeffs[i] * legendre_coeffs[i];
    }
    
    if (sum_low < 1e-30) {
        return 1.0;
    }
    
    double ratio = sum_high / sum_low;
    return std::exp(-ratio);
}

int PAdaptivityManager::decideOrderChange(double smoothness, double error_indicator,
                                          int current_order) const {
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

void PAdaptivityManager::increaseOrder(const std::vector<double>& coeffs_p,
                                       std::vector<double>& coeffs_p1,
                                       int p) const {
    size_t new_size = coeffs_p.size() + (p + 2);
    coeffs_p1.resize(new_size, 0.0);
    
    std::copy(coeffs_p.begin(), coeffs_p.end(), coeffs_p1.begin());
}

void PAdaptivityManager::decreaseOrder(const std::vector<double>& coeffs_p,
                                       std::vector<double>& coeffs_p1,
                                       int p) const {
    if (p <= params_.min_order) {
        coeffs_p1 = coeffs_p;
        return;
    }
    
    size_t new_size = coeffs_p.size() - p;
    if (new_size < 1) new_size = 1;
    
    coeffs_p1.resize(new_size);
    std::copy(coeffs_p.begin(), coeffs_p.begin() + new_size, coeffs_p1.begin());
}


// =============================================================================
// PhysicsBasedPreconditioner Implementation
// =============================================================================

PhysicsBasedPreconditioner::PhysicsBasedPreconditioner() {}

void PhysicsBasedPreconditioner::setParameters(const Parameters& params) {
    params_ = params;
}

void PhysicsBasedPreconditioner::configure(const std::map<std::string, std::string>& config) {
    auto it = config.find("precond_type");
    if (it != config.end()) {
        const std::string& val = it->second;
        if (val == "block_jacobi") params_.type = CoupledPreconditioner::BLOCK_JACOBI;
        else if (val == "block_gauss_seidel") params_.type = CoupledPreconditioner::BLOCK_GAUSS_SEIDEL;
        else if (val == "schur") params_.type = CoupledPreconditioner::SCHUR_COMPLEMENT;
        else if (val == "field_split") params_.type = CoupledPreconditioner::FIELD_SPLIT;
        else if (val == "physics") params_.type = CoupledPreconditioner::PHYSICS_BASED;
        else if (val == "multigrid" || val == "mg") params_.type = CoupledPreconditioner::MULTIGRID;
    }
    
    it = config.find("precond_split_type");
    if (it != config.end()) {
        params_.split_type = it->second;
    }
    
    it = config.find("precond_schur_precond");
    if (it != config.end()) {
        params_.schur_precond = it->second;
    }
    
    it = config.find("precond_mg_levels");
    if (it != config.end()) {
        params_.mg_levels = std::stoi(it->second);
    }
    
    it = config.find("precond_mg_smoother");
    if (it != config.end()) {
        params_.mg_smoother = it->second;
    }
}

void PhysicsBasedPreconditioner::applyBlockJacobi(
    const std::vector<std::vector<double>>& blocks,
    const std::vector<double>& rhs,
    std::vector<double>& solution) const {
    
    if (blocks.empty() || rhs.empty()) {
        return;
    }
    
    size_t n_blocks = blocks.size();
    size_t total_size = rhs.size();
    solution.resize(total_size, 0.0);
    
    size_t block_size = total_size / n_blocks;
    
    for (size_t b = 0; b < n_blocks; ++b) {
        size_t start = b * block_size;
        size_t end = (b == n_blocks - 1) ? total_size : (b + 1) * block_size;
        size_t size = end - start;
        
        const std::vector<double>& block = blocks[b];
        if (block.size() < size * size) continue;
        
        for (size_t i = 0; i < size; ++i) {
            double diag = block[i * size + i];
            if (std::abs(diag) > 1e-30) {
                solution[start + i] = rhs[start + i] / diag;
            }
        }
    }
}

void PhysicsBasedPreconditioner::applyBlockGaussSeidel(
    const std::vector<std::vector<double>>& diag_blocks,
    const std::vector<std::vector<double>>& off_diag,
    const std::vector<double>& rhs,
    std::vector<double>& solution) const {
    
    if (diag_blocks.empty() || rhs.empty()) {
        return;
    }
    
    size_t n_blocks = diag_blocks.size();
    size_t total_size = rhs.size();
    solution.resize(total_size, 0.0);
    
    size_t block_size = total_size / n_blocks;
    
    for (size_t b = 0; b < n_blocks; ++b) {
        size_t start = b * block_size;
        size_t end = (b == n_blocks - 1) ? total_size : (b + 1) * block_size;
        size_t size = end - start;
        
        std::vector<double> r(size);
        for (size_t i = 0; i < size; ++i) {
            r[i] = rhs[start + i];
        }
        
        const std::vector<double>& block = diag_blocks[b];
        if (block.size() >= size * size) {
            for (size_t i = 0; i < size; ++i) {
                double diag = block[i * size + i];
                if (std::abs(diag) > 1e-30) {
                    solution[start + i] = r[i] / diag;
                }
            }
        }
    }
}

void PhysicsBasedPreconditioner::applySchurComplement(
    const std::vector<std::vector<double>>& A,
    const std::vector<std::vector<double>>& B,
    const std::vector<std::vector<double>>& C,
    const std::vector<std::vector<double>>& D,
    const std::vector<double>& f,
    const std::vector<double>& g,
    std::vector<double>& x,
    std::vector<double>& y) const {
    
    size_t nx = f.size();
    size_t ny = g.size();
    
    if (nx == 0 || ny == 0) {
        return;
    }
    
    std::vector<double> z(nx, 0.0);
    if (!A.empty() && A[0].size() >= nx) {
        for (size_t i = 0; i < nx; ++i) {
            double diag = A[0][i * nx + i];
            if (std::abs(diag) > 1e-30) {
                z[i] = f[i] / diag;
            }
        }
    }
    
    std::vector<double> g_mod = g;
    
    y.resize(ny, 0.0);
    if (!D.empty() && D[0].size() >= ny) {
        for (size_t i = 0; i < ny; ++i) {
            double diag = D[0][i * ny + i];
            if (std::abs(diag) > 1e-30) {
                y[i] = g_mod[i] / diag;
            }
        }
    }
    
    std::vector<double> f_mod = f;
    x.resize(nx, 0.0);
    if (!A.empty() && A[0].size() >= nx) {
        for (size_t i = 0; i < nx; ++i) {
            double diag = A[0][i * nx + i];
            if (std::abs(diag) > 1e-30) {
                x[i] = f_mod[i] / diag;
            }
        }
    }
}

void PhysicsBasedPreconditioner::setupFieldSplit(PC pc, const std::vector<IS>& index_sets) const {
    // PETSc configuration would go here
    (void)pc;
    (void)index_sets;
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

bool HighFidelityNumericsConfig::validate(std::string& error_msg) const {
    error_msg.clear();
    
    if (enable_advanced_time) {
        if (time_params.theta < 0.0 || time_params.theta > 1.0) {
            error_msg = "Time integration theta must be in [0, 1]";
            return false;
        }
        if (time_params.newmark_beta < 0.0) {
            error_msg = "Newmark beta must be >= 0";
            return false;
        }
        if (time_params.rho_infinity < 0.0 || time_params.rho_infinity > 1.0) {
            error_msg = "Generalized-alpha rho_infinity must be in [0, 1]";
            return false;
        }
    }
    
    if (enable_p_adapt) {
        if (p_adapt_params.min_order < 1) {
            error_msg = "p-adaptivity min_order must be >= 1";
            return false;
        }
        if (p_adapt_params.max_order < p_adapt_params.min_order) {
            error_msg = "p-adaptivity max_order must be >= min_order";
            return false;
        }
    }
    
    return true;
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

std::unique_ptr<AdvancedTimeIntegration> createTimeIntegration(
    const std::map<std::string, std::string>& config) {
    
    auto time = std::make_unique<AdvancedTimeIntegration>();
    time->configure(config);
    return time;
}

// Note: For DG, use createDGSolver() from DiscontinuousGalerkin.hpp
// Note: For AMR, use AdaptiveMeshRefinement class from AdaptiveMeshRefinement.hpp

} // namespace HighFidelity
} // namespace FSRM
