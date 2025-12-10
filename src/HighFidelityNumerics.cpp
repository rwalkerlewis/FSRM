/**
 * @file HighFidelityNumerics.cpp
 * @brief Implementation of advanced numerical methods for high-fidelity simulations
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
// DGStabilization Implementation
// =============================================================================

DGStabilization::DGStabilization() {}

void DGStabilization::setParameters(const Parameters& params) {
    params_ = params;
}

void DGStabilization::configure(const std::map<std::string, std::string>& config) {
    auto it = config.find("dg_flux_type");
    if (it != config.end()) {
        const std::string& val = it->second;
        if (val == "central") params_.flux_type = DGFluxType::CENTRAL;
        else if (val == "upwind") params_.flux_type = DGFluxType::UPWIND;
        else if (val == "lax_friedrichs") params_.flux_type = DGFluxType::LAX_FRIEDRICHS;
        else if (val == "roe") params_.flux_type = DGFluxType::ROE;
        else if (val == "hllc") params_.flux_type = DGFluxType::HLLC;
        else if (val == "interior_penalty") params_.flux_type = DGFluxType::INTERIOR_PENALTY;
        else if (val == "nipg") params_.flux_type = DGFluxType::NIPG;
        else if (val == "bassi_rebay") params_.flux_type = DGFluxType::BASSI_REBAY;
    }
    
    it = config.find("dg_penalty_type");
    if (it != config.end()) {
        const std::string& val = it->second;
        if (val == "sipg") params_.penalty_type = DGPenaltyType::SIPG;
        else if (val == "nipg") params_.penalty_type = DGPenaltyType::NIPG;
        else if (val == "iipg") params_.penalty_type = DGPenaltyType::IIPG;
        else if (val == "br1") params_.penalty_type = DGPenaltyType::BASSI_REBAY_1;
        else if (val == "br2") params_.penalty_type = DGPenaltyType::BASSI_REBAY_2;
        else if (val == "ldg") params_.penalty_type = DGPenaltyType::LDG;
    }
    
    it = config.find("dg_sigma");
    if (it != config.end()) {
        params_.sigma = std::stod(it->second);
    }
    
    it = config.find("dg_auto_penalty");
    if (it != config.end()) {
        params_.auto_penalty = (it->second == "true" || it->second == "1");
    }
    
    it = config.find("dg_polynomial_order");
    if (it != config.end()) {
        params_.polynomial_order = std::stoi(it->second);
    }
    
    it = config.find("dg_use_limiter");
    if (it != config.end()) {
        params_.use_limiter = (it->second == "true" || it->second == "1");
    }
    
    it = config.find("dg_limiter_type");
    if (it != config.end()) {
        params_.limiter_type = it->second;
    }
    
    it = config.find("dg_tvb_constant");
    if (it != config.end()) {
        params_.tvb_constant = std::stod(it->second);
    }
    
    it = config.find("dg_upwind_fraction");
    if (it != config.end()) {
        params_.upwind_fraction = std::stod(it->second);
    }
}

double DGStabilization::calculatePenalty(double h, int p) const {
    if (params_.auto_penalty) {
        // Standard formula: σ = C * (p+1)^2 / h
        double C = params_.penalty_scaling * params_.sigma;
        return C * (p + 1) * (p + 1) / h;
    }
    return params_.sigma;
}

double DGStabilization::advectionFlux(double u_plus, double u_minus,
                                      const std::array<double, 3>& velocity,
                                      const std::array<double, 3>& normal) const {
    double a_dot_n = dot3(velocity, normal);
    
    switch (params_.flux_type) {
        case DGFluxType::CENTRAL:
            // Central: F = a · {{u}} = a · 0.5·(u⁺ + u⁻)
            return a_dot_n * 0.5 * (u_plus + u_minus);
            
        case DGFluxType::UPWIND:
        default: {
            // Full upwind
            if (a_dot_n >= 0.0) {
                return a_dot_n * u_minus;  // Use upwind value
            } else {
                return a_dot_n * u_plus;
            }
        }
    }
    
    // Blended central/upwind based on upwind_fraction
    double central = a_dot_n * 0.5 * (u_plus + u_minus);
    double upwind;
    if (a_dot_n >= 0.0) {
        upwind = a_dot_n * u_minus;
    } else {
        upwind = a_dot_n * u_plus;
    }
    return params_.upwind_fraction * upwind + (1.0 - params_.upwind_fraction) * central;
}

double DGStabilization::diffusionFlux(double u_plus, double u_minus,
                                      const std::array<double, 3>& grad_u_plus,
                                      const std::array<double, 3>& grad_u_minus,
                                      double diffusivity, double h,
                                      const std::array<double, 3>& normal) const {
    // Average gradient · normal
    std::array<double, 3> avg_grad = {
        0.5 * (grad_u_plus[0] + grad_u_minus[0]),
        0.5 * (grad_u_plus[1] + grad_u_minus[1]),
        0.5 * (grad_u_plus[2] + grad_u_minus[2])
    };
    double avg_grad_n = dot3(avg_grad, normal);
    
    // Jump in solution
    double jump_u = u_plus - u_minus;
    
    // Penalty parameter
    double sigma = calculatePenalty(h, params_.polynomial_order);
    
    double symmetry_factor = 0.0;
    switch (params_.penalty_type) {
        case DGPenaltyType::SIPG:
            symmetry_factor = -1.0;  // Symmetric
            break;
        case DGPenaltyType::NIPG:
            symmetry_factor = 1.0;   // Non-symmetric
            break;
        case DGPenaltyType::IIPG:
            symmetry_factor = 0.0;   // Incomplete
            break;
        default:
            symmetry_factor = -1.0;
    }
    
    // SIPG flux: {{κ∇u}}·n - σ/h·[[u]] + symmetry·{{κ∇v}}·[[u]]
    // Here we return the numerical flux contribution
    return diffusivity * (avg_grad_n - sigma * jump_u / h);
}

double DGStabilization::laxFriedrichsFlux(double u_plus, double u_minus,
                                          double flux_plus, double flux_minus,
                                          double max_eigenvalue) const {
    // F = 0.5·(F⁺ + F⁻) - 0.5·λ_max·(u⁺ - u⁻)
    return 0.5 * (flux_plus + flux_minus) - 0.5 * max_eigenvalue * (u_plus - u_minus);
}

void DGStabilization::applyLimiter(std::vector<double>& solution,
                                   const std::vector<std::array<double, 3>>& gradients,
                                   const std::vector<double>& cell_sizes) const {
    if (!params_.use_limiter) return;
    
    size_t n_cells = solution.size();
    if (n_cells == 0) return;
    
    // For now, implement a simple gradient limiting
    // More sophisticated WENO limiters would go here
    
    for (size_t i = 0; i < n_cells; ++i) {
        double h = cell_sizes[i];
        double grad_mag = norm3(gradients[i]);
        
        if (params_.limiter_type == "minmod") {
            // Apply minmod limiting
            // This is a simplified version; full implementation would need neighbor info
            // Just demonstrate the structure
        } else if (params_.limiter_type == "tvb") {
            // TVB limiting allows larger slopes based on M·h²
            double M = params_.tvb_constant;
            // If |slope| < M·h², no limiting needed
        }
    }
}

double DGStabilization::minmod(double a, double b) const {
    if (a * b <= 0.0) {
        return 0.0;
    }
    if (std::abs(a) < std::abs(b)) {
        return a;
    }
    return b;
}

double DGStabilization::minmod(double a, double b, double c) const {
    if (a > 0.0 && b > 0.0 && c > 0.0) {
        return std::min({a, b, c});
    }
    if (a < 0.0 && b < 0.0 && c < 0.0) {
        return std::max({a, b, c});
    }
    return 0.0;
}

double DGStabilization::tvbMinmod(double a, double b, double c, double h) const {
    double M = params_.tvb_constant;
    if (std::abs(a) <= M * h * h) {
        return a;
    }
    return minmod(a, b, c);
}

void DGStabilization::faceIntegral(int face_id,
                                   const std::vector<double>& u_plus,
                                   const std::vector<double>& u_minus,
                                   const std::array<double, 3>& normal,
                                   double face_area, double h_plus, double h_minus,
                                   std::vector<double>& residual_plus,
                                   std::vector<double>& residual_minus) const {
    // Compute face contributions for DG assembly
    // This is a general template; actual implementation depends on equation type
    
    size_t n_nodes = u_plus.size();
    if (residual_plus.size() != n_nodes) {
        residual_plus.resize(n_nodes, 0.0);
    }
    if (residual_minus.size() != n_nodes) {
        residual_minus.resize(n_nodes, 0.0);
    }
    
    // Compute average element size for penalty
    double h = 0.5 * (h_plus + h_minus);
    double sigma = calculatePenalty(h, params_.polynomial_order);
    
    // For each quadrature point on face (simplified to single point here)
    for (size_t i = 0; i < n_nodes; ++i) {
        // Jump term contribution
        double jump = u_plus[i] - u_minus[i];
        
        // Add to residuals (signs follow convention: normal points from minus to plus)
        residual_plus[i] -= sigma * jump * face_area / n_nodes;
        residual_minus[i] += sigma * jump * face_area / n_nodes;
    }
}


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
        // Low Peclet: no artificial diffusion needed
        return 0.0;
    }
    
    // Optimal artificial diffusivity:
    // κ_art = 0.5·|a|·h·(coth(Pe) - 1/Pe)
    double coth_Pe = (std::exp(2.0 * Pe) + 1.0) / (std::exp(2.0 * Pe) - 1.0);
    double alpha = coth_Pe - 1.0 / Pe;  // ∈ (0, 1)
    
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
    
    // SUPG: τ·(a·∇v)·R
    double a_grad_v = dot3(velocity, test_gradient);
    stabilization_term = tau * a_grad_v * residual_strong;
}

void SUPGStabilization::addGLSStabilization(double L_adjoint_v,
                                            double residual_strong,
                                            double tau,
                                            double& stabilization_term) const {
    // GLS: τ·L*(v)·R
    // Where L* is the adjoint operator
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
    
    // Discontinuity-capturing diffusivity:
    // κ_DC = C · h · |R| / |∇u|
    double kappa_DC = params_.shock_capture_coeff * h * std::abs(residual) / grad_mag;
    
    return kappa_DC;
}

double SUPGStabilization::pspgTerm(const std::array<double, 3>& grad_p,
                                   const std::array<double, 3>& grad_q,
                                   double tau) const {
    // PSPG: τ·∇p·∇q
    return tau * dot3(grad_p, grad_q);
}


// =============================================================================
// AdvancedTimeIntegration Implementation
// =============================================================================

AdvancedTimeIntegration::AdvancedTimeIntegration() {
    // Set default Newmark parameters for trapezoidal rule
    params_.newmark_beta = 0.25;
    params_.newmark_gamma = 0.5;
}

void AdvancedTimeIntegration::setParameters(const Parameters& params) {
    params_ = params;
    
    // Auto-compute generalized-α coefficients if needed
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
    
    it = config.find("time_safety_factor");
    if (it != config.end()) {
        params_.safety_factor = std::stod(it->second);
    }
    
    // Compute derived coefficients
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
        case TimeScheme::ADAMS_BASHFORTH_2:
        case TimeScheme::ADAMS_MOULTON_2:
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
            // θ = 0.5 gives 2nd order, otherwise 1st
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
        case TimeScheme::BDF4:  // BDF4 is not A-stable (outside A(73°))
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
            return false;  // A-stable but not L-stable
            
        case TimeScheme::THETA_METHOD:
            return params_.theta > 0.5;
            
        case TimeScheme::NEWMARK_BETA:
        case TimeScheme::HHT_ALPHA:
        case TimeScheme::GENERALIZED_ALPHA:
            // L-stable with appropriate parameters
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
            // s = 1
            A = {{0.0}};
            b = {1.0};
            c = {0.0};
            break;
            
        case TimeScheme::BACKWARD_EULER:
            // s = 1, implicit
            A = {{1.0}};
            b = {1.0};
            c = {1.0};
            break;
            
        case TimeScheme::RK2:
            // Midpoint method
            A = {{0.0, 0.0},
                 {0.5, 0.0}};
            b = {0.0, 1.0};
            c = {0.0, 0.5};
            break;
            
        case TimeScheme::SSPRK3:
            // Strong stability preserving RK3
            A = {{0.0, 0.0, 0.0},
                 {1.0, 0.0, 0.0},
                 {0.25, 0.25, 0.0}};
            b = {1.0/6.0, 1.0/6.0, 2.0/3.0};
            c = {0.0, 1.0, 0.5};
            break;
            
        case TimeScheme::RK4:
            // Classic RK4
            A = {{0.0, 0.0, 0.0, 0.0},
                 {0.5, 0.0, 0.0, 0.0},
                 {0.0, 0.5, 0.0, 0.0},
                 {0.0, 0.0, 1.0, 0.0}};
            b = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};
            c = {0.0, 0.5, 0.5, 1.0};
            break;
            
        case TimeScheme::DIRK2:
            // 2-stage DIRK (A-L-stable)
            {
                double gamma = 1.0 - 1.0/std::sqrt(2.0);
                A = {{gamma, 0.0},
                     {1.0 - gamma, gamma}};
                b = {1.0 - gamma, gamma};
                c = {gamma, 1.0};
            }
            break;
            
        case TimeScheme::DIRK3:
            // 3-stage DIRK (Alexander's method)
            {
                double alpha = 0.4358665215;  // Root of cubic
                double tau2 = (1.0 + alpha) / 2.0;
                double b1 = -(6.0*alpha*alpha - 16.0*alpha + 1.0) / 4.0;
                double b2 = (6.0*alpha*alpha - 20.0*alpha + 5.0) / 4.0;
                A = {{alpha, 0.0, 0.0},
                     {tau2 - alpha, alpha, 0.0},
                     {b1, b2, alpha}};
                b = {b1, b2, alpha};
                c = {alpha, tau2, 1.0};
            }
            break;
            
        default:
            // Default to forward Euler
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
    // Compute initial acceleration from equilibrium: M·a0 = F - C·v0 - K·u0
    // This is passed through the compute_rhs function
    a0 = compute_rhs(u0, v0);
}

void AdvancedTimeIntegration::newmarkPredict(const std::vector<double>& u_n,
                                             const std::vector<double>& v_n,
                                             const std::vector<double>& a_n,
                                             double dt,
                                             std::vector<double>& u_pred,
                                             std::vector<double>& v_pred) const {
    // Newmark predictor:
    // ũ = u_n + dt·v_n + dt²·(0.5 - β)·a_n
    // ṽ = v_n + dt·(1 - γ)·a_n
    
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
    // Newmark corrector:
    // u_{n+1} = ũ + β·dt²·a_{n+1}
    // v_{n+1} = ṽ + γ·dt·a_{n+1}
    
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
    // Effective mass for Newmark:
    // M_eff = M + γ·dt·C + β·dt²·K
    
    double beta = params_.newmark_beta;
    double gamma = params_.newmark_gamma;
    
    M_eff = mass + gamma * dt * damping + beta * dt * dt * stiffness;
}

void AdvancedTimeIntegration::computeGeneralizedAlphaCoeffs() {
    // Compute coefficients from spectral radius at infinity
    double rho = params_.rho_infinity;
    
    params_.alpha_m = (2.0 * rho - 1.0) / (rho + 1.0);
    params_.alpha_f = rho / (rho + 1.0);
    params_.newmark_gamma = 0.5 - params_.alpha_m + params_.alpha_f;
    params_.newmark_beta = 0.25 * std::pow(1.0 - params_.alpha_m + params_.alpha_f, 2);
}

double AdvancedTimeIntegration::estimateError(const std::vector<double>& u_high,
                                              const std::vector<double>& u_low) const {
    // Error estimate from difference of high and low order solutions
    // ||u_high - u_low||
    
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
    
    return std::sqrt(error_sq / norm_sq);  // Relative error
}

double AdvancedTimeIntegration::computeNewTimestep(double dt_current, double error,
                                                   double tolerance) const {
    if (error < 1e-30) {
        return dt_current * params_.max_dt_factor;
    }
    
    // Standard PI controller for time step adaptation
    int order = getOrder();
    double factor = params_.safety_factor * std::pow(tolerance / error, 1.0 / (order + 1));
    
    // Apply limits
    factor = clamp(factor, params_.min_dt_factor, params_.max_dt_factor);
    
    return dt_current * factor;
}

void AdvancedTimeIntegration::bdfCoefficients(int order, std::vector<double>& alpha,
                                              double& beta) const {
    // BDF coefficients: Σ α_k·u_{n+1-k} = β·dt·f(u_{n+1})
    
    alpha.clear();
    
    switch (order) {
        case 1:
            // BDF1 = Backward Euler
            alpha = {1.0, -1.0};
            beta = 1.0;
            break;
            
        case 2:
            // BDF2
            alpha = {3.0/2.0, -2.0, 1.0/2.0};
            beta = 1.0;
            break;
            
        case 3:
            // BDF3
            alpha = {11.0/6.0, -3.0, 3.0/2.0, -1.0/3.0};
            beta = 1.0;
            break;
            
        case 4:
            // BDF4
            alpha = {25.0/12.0, -4.0, 3.0, -4.0/3.0, 1.0/4.0};
            beta = 1.0;
            break;
            
        default:
            // Default to BDF1
            alpha = {1.0, -1.0};
            beta = 1.0;
    }
}


// =============================================================================
// AdaptiveMeshRefinementManager Implementation
// =============================================================================

AdaptiveMeshRefinementManager::AdaptiveMeshRefinementManager() {}

void AdaptiveMeshRefinementManager::setParameters(const Parameters& params) {
    params_ = params;
}

void AdaptiveMeshRefinementManager::configure(const std::map<std::string, std::string>& config) {
    auto it = config.find("amr_indicator");
    if (it != config.end()) {
        const std::string& val = it->second;
        if (val == "gradient") params_.indicator = RefinementIndicator::GRADIENT;
        else if (val == "residual") params_.indicator = RefinementIndicator::RESIDUAL;
        else if (val == "recovery" || val == "zz") params_.indicator = RefinementIndicator::RECOVERY;
        else if (val == "hessian") params_.indicator = RefinementIndicator::HESSIAN;
        else if (val == "jumps") params_.indicator = RefinementIndicator::JUMPS;
        else if (val == "goal_oriented" || val == "dwr") params_.indicator = RefinementIndicator::GOAL_ORIENTED;
        else if (val == "physics") params_.indicator = RefinementIndicator::PHYSICS_BASED;
    }
    
    it = config.find("amr_strategy");
    if (it != config.end()) {
        const std::string& val = it->second;
        if (val == "fixed_fraction") params_.strategy = RefinementStrategy::FIXED_FRACTION;
        else if (val == "fixed_number") params_.strategy = RefinementStrategy::FIXED_NUMBER;
        else if (val == "threshold") params_.strategy = RefinementStrategy::THRESHOLD;
        else if (val == "equilibration") params_.strategy = RefinementStrategy::EQUILIBRATION;
        else if (val == "kelly") params_.strategy = RefinementStrategy::KELLY;
    }
    
    it = config.find("amr_refine_fraction");
    if (it != config.end()) {
        params_.refine_fraction = std::stod(it->second);
    }
    
    it = config.find("amr_coarsen_fraction");
    if (it != config.end()) {
        params_.coarsen_fraction = std::stod(it->second);
    }
    
    it = config.find("amr_refine_threshold");
    if (it != config.end()) {
        params_.refine_threshold = std::stod(it->second);
    }
    
    it = config.find("amr_coarsen_threshold");
    if (it != config.end()) {
        params_.coarsen_threshold = std::stod(it->second);
    }
    
    it = config.find("amr_max_level");
    if (it != config.end()) {
        params_.max_refinement_level = std::stoi(it->second);
    }
    
    it = config.find("amr_min_element_size");
    if (it != config.end()) {
        params_.min_element_size = std::stod(it->second);
    }
    
    it = config.find("amr_ensure_conforming");
    if (it != config.end()) {
        params_.ensure_conforming = (it->second == "true" || it->second == "1");
    }
}

double AdaptiveMeshRefinementManager::gradientIndicator(const std::array<double, 3>& gradient,
                                                        double h) const {
    // η_K = h_K · ||∇u||_K
    return h * norm3(gradient);
}

double AdaptiveMeshRefinementManager::residualIndicator(double element_residual,
                                                        const std::vector<double>& face_jumps,
                                                        double h) const {
    // η_K² = h²·||R||² + h·Σ||J||²
    
    double eta_sq = h * h * element_residual * element_residual;
    
    for (double jump : face_jumps) {
        eta_sq += h * jump * jump;
    }
    
    return std::sqrt(eta_sq);
}

double AdaptiveMeshRefinementManager::recoveryIndicator(const std::array<double, 3>& raw_gradient,
                                                        const std::array<double, 3>& recovered_gradient) const {
    // Zienkiewicz-Zhu: ||∇u_h - G(∇u_h)||
    std::array<double, 3> diff = {
        raw_gradient[0] - recovered_gradient[0],
        raw_gradient[1] - recovered_gradient[1],
        raw_gradient[2] - recovered_gradient[2]
    };
    return norm3(diff);
}

double AdaptiveMeshRefinementManager::goalOrientedIndicator(double primal_residual,
                                                            double dual_weight) const {
    // DWR: η_K = |R(u_h)·z_h|
    return std::abs(primal_residual * dual_weight);
}

std::vector<int> AdaptiveMeshRefinementManager::markElements(
    const std::vector<double>& indicators) const {
    
    std::vector<int> flags(indicators.size(), 0);
    
    if (indicators.empty()) {
        return flags;
    }
    
    switch (params_.strategy) {
        case RefinementStrategy::FIXED_FRACTION:
            fixedFractionStrategy(indicators, flags);
            break;
            
        case RefinementStrategy::THRESHOLD:
        case RefinementStrategy::KELLY:
            thresholdStrategy(indicators, flags);
            break;
            
        case RefinementStrategy::FIXED_NUMBER:
        case RefinementStrategy::EQUILIBRATION:
        default:
            fixedFractionStrategy(indicators, flags);
            break;
    }
    
    return flags;
}

void AdaptiveMeshRefinementManager::fixedFractionStrategy(
    const std::vector<double>& indicators,
    std::vector<int>& flags) const {
    
    size_t n = indicators.size();
    
    // Sort indicators to find thresholds
    std::vector<double> sorted_indicators = indicators;
    std::sort(sorted_indicators.begin(), sorted_indicators.end());
    
    // Find threshold for refinement (top refine_fraction%)
    size_t refine_idx = static_cast<size_t>((1.0 - params_.refine_fraction) * n);
    double refine_threshold = sorted_indicators[std::min(refine_idx, n - 1)];
    
    // Find threshold for coarsening (bottom coarsen_fraction%)
    size_t coarsen_idx = static_cast<size_t>(params_.coarsen_fraction * n);
    double coarsen_threshold = sorted_indicators[std::min(coarsen_idx, n - 1)];
    
    // Mark elements
    for (size_t i = 0; i < n; ++i) {
        if (indicators[i] >= refine_threshold) {
            flags[i] = 1;   // Refine
        } else if (indicators[i] <= coarsen_threshold) {
            flags[i] = -1;  // Coarsen
        } else {
            flags[i] = 0;   // Keep
        }
    }
}

void AdaptiveMeshRefinementManager::thresholdStrategy(
    const std::vector<double>& indicators,
    std::vector<int>& flags) const {
    
    for (size_t i = 0; i < indicators.size(); ++i) {
        if (indicators[i] > params_.refine_threshold) {
            flags[i] = 1;   // Refine
        } else if (indicators[i] < params_.coarsen_threshold) {
            flags[i] = -1;  // Coarsen
        } else {
            flags[i] = 0;   // Keep
        }
    }
}

void AdaptiveMeshRefinementManager::enforceConformity(
    std::vector<int>& refinement_flags,
    const std::vector<std::vector<int>>& neighbors) const {
    
    if (!params_.ensure_conforming) {
        return;
    }
    
    // Iteratively propagate refinement to avoid hanging nodes
    bool changed = true;
    int max_iterations = 10;
    int iteration = 0;
    
    while (changed && iteration < max_iterations) {
        changed = false;
        ++iteration;
        
        for (size_t i = 0; i < refinement_flags.size(); ++i) {
            if (refinement_flags[i] == 1) {
                // This element will be refined
                // Check neighbors
                for (int neighbor : neighbors[i]) {
                    if (neighbor >= 0 && neighbor < static_cast<int>(refinement_flags.size())) {
                        // If neighbor is marked for coarsening, cancel that
                        if (refinement_flags[neighbor] == -1) {
                            refinement_flags[neighbor] = 0;
                            changed = true;
                        }
                        // Propagate refinement if level difference would be too large
                        // (In a full implementation, we'd track levels per element)
                    }
                }
            }
        }
    }
}

double AdaptiveMeshRefinementManager::globalErrorEstimate(
    const std::vector<double>& indicators) const {
    
    // ||e||² ≈ Σ η_K²
    double sum_sq = 0.0;
    for (double eta : indicators) {
        sum_sq += eta * eta;
    }
    return std::sqrt(sum_sq);
}

double AdaptiveMeshRefinementManager::effectivityIndex(double estimated_error,
                                                       double true_error) const {
    if (true_error < 1e-30) {
        return 1.0;  // Cannot compute
    }
    return estimated_error / true_error;
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
    // Smoothness based on decay of Legendre coefficients
    // Fast decay = smooth solution
    
    if (legendre_coeffs.size() < 2) {
        return 1.0;  // Assume smooth
    }
    
    // Compute decay rate of coefficients
    // For smooth functions: |a_k| ~ C·q^k with q < 1
    // For non-smooth: slower decay
    
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
    
    // Ratio of high to low order coefficients
    // Small ratio = smooth, large ratio = rough
    double ratio = sum_high / sum_low;
    
    // Map to [0, 1] smoothness indicator
    // 1 = very smooth, 0 = rough
    return std::exp(-ratio);
}

int PAdaptivityManager::decideOrderChange(double smoothness, double error_indicator,
                                          int current_order) const {
    // Increase order if:
    // - Solution is smooth (high smoothness indicator)
    // - AND error is above threshold
    
    if (smoothness > params_.smoothness_threshold_inc &&
        error_indicator > params_.p_refinement_threshold &&
        current_order < params_.max_order) {
        return 1;  // Increase p
    }
    
    // Decrease order if:
    // - Solution is rough (low smoothness indicator)
    // - OR very low error (no need for high accuracy)
    
    if ((smoothness < params_.smoothness_threshold_dec ||
         error_indicator < params_.p_refinement_threshold * 0.01) &&
        current_order > params_.min_order) {
        return -1;  // Decrease p
    }
    
    return 0;  // Keep current order
}

void PAdaptivityManager::increaseOrder(const std::vector<double>& coeffs_p,
                                       std::vector<double>& coeffs_p1,
                                       int p) const {
    // Project from order p to p+1
    // Simply pad with zeros (L2 projection would preserve existing coefficients)
    
    size_t new_size = coeffs_p.size() + (p + 2);  // Approximate increase
    coeffs_p1.resize(new_size, 0.0);
    
    std::copy(coeffs_p.begin(), coeffs_p.end(), coeffs_p1.begin());
}

void PAdaptivityManager::decreaseOrder(const std::vector<double>& coeffs_p,
                                       std::vector<double>& coeffs_p1,
                                       int p) const {
    // Project from order p to p-1
    // Truncate high-order coefficients
    
    if (p <= params_.min_order) {
        coeffs_p1 = coeffs_p;
        return;
    }
    
    size_t new_size = coeffs_p.size() - p;  // Approximate decrease
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
    
    it = config.find("precond_mg_smooth_sweeps");
    if (it != config.end()) {
        params_.mg_smooth_sweeps = std::stoi(it->second);
    }
}

void PhysicsBasedPreconditioner::applyBlockJacobi(
    const std::vector<std::vector<double>>& blocks,
    const std::vector<double>& rhs,
    std::vector<double>& solution) const {
    
    // Block Jacobi: Solve each block independently
    // x_i = D_i^(-1) * b_i
    
    if (blocks.empty() || rhs.empty()) {
        return;
    }
    
    size_t n_blocks = blocks.size();
    size_t total_size = rhs.size();
    solution.resize(total_size, 0.0);
    
    // Compute block sizes
    size_t block_size = total_size / n_blocks;
    
    for (size_t b = 0; b < n_blocks; ++b) {
        size_t start = b * block_size;
        size_t end = (b == n_blocks - 1) ? total_size : (b + 1) * block_size;
        size_t size = end - start;
        
        // For each block, solve D_i * x_i = b_i
        // Simplified: assuming diagonal blocks are stored as dense matrices
        // In practice, this would use LU decomposition or iterative solve
        
        const std::vector<double>& block = blocks[b];
        if (block.size() < size * size) continue;
        
        // Simple diagonal solve (for demonstration)
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
    
    // Block Gauss-Seidel: Sequential solve with coupling
    // For block i: D_i * x_i = b_i - Σ_{j<i} L_ij * x_j - Σ_{j>i} U_ij * x_j^{old}
    
    if (diag_blocks.empty() || rhs.empty()) {
        return;
    }
    
    size_t n_blocks = diag_blocks.size();
    size_t total_size = rhs.size();
    solution.resize(total_size, 0.0);
    
    size_t block_size = total_size / n_blocks;
    
    // Forward sweep
    for (size_t b = 0; b < n_blocks; ++b) {
        size_t start = b * block_size;
        size_t end = (b == n_blocks - 1) ? total_size : (b + 1) * block_size;
        size_t size = end - start;
        
        // Compute modified RHS: r_i = b_i - Σ A_ij * x_j
        std::vector<double> r(size);
        for (size_t i = 0; i < size; ++i) {
            r[i] = rhs[start + i];
        }
        
        // Subtract contributions from previously solved blocks
        // (In full implementation, would multiply by off-diagonal blocks)
        
        // Solve D_i * x_i = r_i (simplified diagonal solve)
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
    
    // Schur complement method for saddle-point systems:
    // [A  B] [x]   [f]
    // [C  D] [y] = [g]
    //
    // S = D - C * A^(-1) * B  (Schur complement)
    // S * y = g - C * A^(-1) * f
    // A * x = f - B * y
    
    size_t nx = f.size();
    size_t ny = g.size();
    
    if (nx == 0 || ny == 0) {
        return;
    }
    
    // This is a simplified implementation
    // Full implementation would need actual matrix operations
    
    // Step 1: Solve A * z = f
    std::vector<double> z(nx, 0.0);
    // Simplified: assume A is diagonal
    if (!A.empty() && A[0].size() >= nx) {
        for (size_t i = 0; i < nx; ++i) {
            double diag = A[0][i * nx + i];
            if (std::abs(diag) > 1e-30) {
                z[i] = f[i] / diag;
            }
        }
    }
    
    // Step 2: Compute modified RHS: g_mod = g - C * z
    std::vector<double> g_mod = g;
    // (Would need actual C * z computation)
    
    // Step 3: Solve S * y = g_mod
    // Simplified: assume S ≈ D is diagonal
    y.resize(ny, 0.0);
    if (!D.empty() && D[0].size() >= ny) {
        for (size_t i = 0; i < ny; ++i) {
            double diag = D[0][i * ny + i];
            if (std::abs(diag) > 1e-30) {
                y[i] = g_mod[i] / diag;
            }
        }
    }
    
    // Step 4: Solve A * x = f - B * y
    std::vector<double> f_mod = f;
    // (Would need actual B * y computation)
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
    // Configure PETSc PCFieldSplit
    // This sets up block preconditioning based on physics fields
    
    // Note: This requires PETSc headers and linking
    // The actual implementation would look like:
    /*
    PCSetType(pc, PCFIELDSPLIT);
    
    for (size_t i = 0; i < index_sets.size(); ++i) {
        std::string name = (i < params_.block_names.size()) ? 
                          params_.block_names[i] : std::to_string(i);
        PCFieldSplitSetIS(pc, name.c_str(), index_sets[i]);
    }
    
    if (params_.split_type == "additive") {
        PCFieldSplitSetType(pc, PC_COMPOSITE_ADDITIVE);
    } else if (params_.split_type == "multiplicative") {
        PCFieldSplitSetType(pc, PC_COMPOSITE_MULTIPLICATIVE);
    } else if (params_.split_type == "schur") {
        PCFieldSplitSetType(pc, PC_COMPOSITE_SCHUR);
        PCFieldSplitSetSchurPre(pc, PC_FIELDSPLIT_SCHUR_PRE_SELFP, NULL);
    }
    */
    
    // Placeholder - actual PETSc calls would go here
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
    
    // DG stabilization
    it = config.find("enable_dg");
    if (it != config.end()) {
        enable_dg = (it->second == "true" || it->second == "1");
    }
    
    // SUPG stabilization
    it = config.find("enable_supg");
    if (it != config.end()) {
        enable_supg = (it->second == "true" || it->second == "1");
    }
    
    // Advanced time integration
    it = config.find("enable_advanced_time");
    if (it != config.end()) {
        enable_advanced_time = (it->second == "true" || it->second == "1");
    }
    
    // AMR
    it = config.find("enable_amr");
    if (it != config.end()) {
        enable_amr = (it->second == "true" || it->second == "1");
    }
    
    // p-adaptivity
    it = config.find("enable_p_adapt");
    if (it != config.end()) {
        enable_p_adapt = (it->second == "true" || it->second == "1");
    }
    
    // Physics-based preconditioning
    it = config.find("enable_physics_precond");
    if (it != config.end()) {
        enable_physics_precond = (it->second == "true" || it->second == "1");
    }
    
    // Parse sub-configurations
    // (Individual parameters are parsed when the specific classes configure themselves)
}

bool HighFidelityNumericsConfig::validate(std::string& error_msg) const {
    error_msg.clear();
    
    // Check for incompatible combinations
    if (enable_dg && enable_supg) {
        // DG and SUPG can coexist but need careful treatment
        // For now, allow both
    }
    
    // Validate AMR parameters
    if (enable_amr) {
        if (amr_params.refine_fraction < 0.0 || amr_params.refine_fraction > 1.0) {
            error_msg = "AMR refine_fraction must be in [0, 1]";
            return false;
        }
        if (amr_params.coarsen_fraction < 0.0 || amr_params.coarsen_fraction > 1.0) {
            error_msg = "AMR coarsen_fraction must be in [0, 1]";
            return false;
        }
        if (amr_params.refine_fraction + amr_params.coarsen_fraction > 1.0) {
            error_msg = "AMR refine_fraction + coarsen_fraction must be <= 1";
            return false;
        }
    }
    
    // Validate time integration parameters
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
    
    // Validate p-adaptivity parameters
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

std::unique_ptr<DGStabilization> createDGStabilization(
    const std::map<std::string, std::string>& config) {
    
    auto dg = std::make_unique<DGStabilization>();
    dg->configure(config);
    return dg;
}

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

std::unique_ptr<AdaptiveMeshRefinementManager> createAMRManager(
    const std::map<std::string, std::string>& config) {
    
    auto amr = std::make_unique<AdaptiveMeshRefinementManager>();
    amr->configure(config);
    return amr;
}

} // namespace HighFidelity
} // namespace FSRM
