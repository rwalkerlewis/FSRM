/**
 * @file PlasticityModel.cpp
 * @brief Full implementation of plasticity models with return mapping algorithms
 * 
 * Based on Simo & Hughes (1998) "Computational Inelasticity"
 * and de Souza Neto et al. (2008) "Computational Methods for Plasticity"
 */

#include "PlasticityModel.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>

namespace FSRM {

// =============================================================================
// Utility functions for stress tensor operations
// =============================================================================

namespace {

// Compute first invariant I1 = tr(σ)
double computeI1(const std::array<double, 6>& stress) {
    return stress[0] + stress[1] + stress[2];
}

// Compute second invariant of deviatoric stress J2
double computeJ2(const std::array<double, 6>& stress) {
    double I1 = computeI1(stress);
    double p = I1 / 3.0;
    
    // Deviatoric stress: s = σ - p*I
    double s11 = stress[0] - p;
    double s22 = stress[1] - p;
    double s33 = stress[2] - p;
    double s12 = stress[3];
    double s13 = stress[4];
    double s23 = stress[5];
    
    // J2 = 0.5 * s_ij * s_ij
    return 0.5 * (s11*s11 + s22*s22 + s33*s33) + 
           s12*s12 + s13*s13 + s23*s23;
}

// Compute √(2*J2) - equivalent shear stress
double computeEqShearStress(const std::array<double, 6>& stress) {
    double J2 = computeJ2(stress);
    return std::sqrt(2.0 * J2);
}

// Compute von Mises equivalent stress √(3*J2)
double computeVonMises(const std::array<double, 6>& stress) {
    double J2 = computeJ2(stress);
    return std::sqrt(3.0 * J2);
}

// Compute mean stress p = I1/3
double computeMeanStress(const std::array<double, 6>& stress) {
    return computeI1(stress) / 3.0;
}

// Compute deviatoric stress
void computeDeviatoricStress(const std::array<double, 6>& stress,
                             std::array<double, 6>& deviatoric) {
    double p = computeMeanStress(stress);
    deviatoric[0] = stress[0] - p;
    deviatoric[1] = stress[1] - p;
    deviatoric[2] = stress[2] - p;
    deviatoric[3] = stress[3];
    deviatoric[4] = stress[4];
    deviatoric[5] = stress[5];
}

// Compute principal stresses (eigenvalues)
void computePrincipalStresses(const std::array<double, 6>& stress,
                              double& sigma1, double& sigma2, double& sigma3) {
    // For general 3D stress tensor, use characteristic equation
    // σ³ - I1*σ² + I2*σ - I3 = 0
    
    double I1 = computeI1(stress);
    
    // I2 = 0.5*(I1² - tr(σ²))
    double I2 = stress[0]*stress[1] + stress[1]*stress[2] + stress[2]*stress[0]
              - stress[3]*stress[3] - stress[4]*stress[4] - stress[5]*stress[5];
    
    // I3 = det(σ)
    double I3 = stress[0]*(stress[1]*stress[2] - stress[5]*stress[5])
              - stress[3]*(stress[3]*stress[2] - stress[5]*stress[4])
              + stress[4]*(stress[3]*stress[5] - stress[1]*stress[4]);
    
    // Solve cubic equation using Cardano's formula
    double a = I1 / 3.0;
    double q = I2 / 3.0 - a * a;
    double r = (I2 * a - I3) / 2.0 - a * a * a;
    
    double disc = q * q * q + r * r;
    
    if (disc > 0) {
        // One real root (shouldn't happen for symmetric tensor)
        double s = std::cbrt(r + std::sqrt(disc));
        double t = std::cbrt(r - std::sqrt(disc));
        sigma1 = sigma2 = sigma3 = a + s + t;
    } else {
        // Three real roots
        double theta = std::acos(r / std::sqrt(-q * q * q)) / 3.0;
        double sqrtq = 2.0 * std::sqrt(-q);
        
        sigma1 = sqrtq * std::cos(theta) + a;
        sigma2 = sqrtq * std::cos(theta - 2.0 * M_PI / 3.0) + a;
        sigma3 = sqrtq * std::cos(theta + 2.0 * M_PI / 3.0) + a;
        
        // Sort: σ1 ≥ σ2 ≥ σ3
        if (sigma1 < sigma2) std::swap(sigma1, sigma2);
        if (sigma2 < sigma3) std::swap(sigma2, sigma3);
        if (sigma1 < sigma2) std::swap(sigma1, sigma2);
    }
}

// Compute norm of strain/stress tensor
double tensorNorm(const std::array<double, 6>& tensor) {
    return std::sqrt(tensor[0]*tensor[0] + tensor[1]*tensor[1] + tensor[2]*tensor[2] +
                    2.0*(tensor[3]*tensor[3] + tensor[4]*tensor[4] + tensor[5]*tensor[5]));
}

// Apply elastic modulus: σ = C : ε
void applyElasticModulus(const std::array<std::array<double, 6>, 6>& C,
                        const std::array<double, 6>& strain,
                        std::array<double, 6>& stress) {
    for (int i = 0; i < 6; ++i) {
        stress[i] = 0.0;
        for (int j = 0; j < 6; ++j) {
            stress[i] += C[i][j] * strain[j];
        }
    }
}

// Parse double from config
double parseDouble(const std::map<std::string, std::string>& config,
                  const std::string& key, double default_val) {
    auto it = config.find(key);
    if (it != config.end() && !it->second.empty()) {
        try {
            return std::stod(it->second);
        } catch (...) {}
    }
    return default_val;
}

// Parse string from config
std::string parseString(const std::map<std::string, std::string>& config,
                       const std::string& key, const std::string& default_val) {
    auto it = config.find(key);
    return (it != config.end() && !it->second.empty()) ? it->second : default_val;
}

} // anonymous namespace

// =============================================================================
// DruckerPragerModel Implementation
// =============================================================================

DruckerPragerModel::DruckerPragerModel() {
    updateInternalParameters();
}

void DruckerPragerModel::setParameters(const Parameters& p) {
    params = p;
    updateInternalParameters();
}

void DruckerPragerModel::updateInternalParameters() {
    // Convert Mohr-Coulomb parameters to Drucker-Prager
    // For plane strain matching (intermediate):
    // α = 2sin(φ) / (√3(3 - sin(φ)))
    // κ = 6c·cos(φ) / (√3(3 - sin(φ)))
    
    double sin_phi = std::sin(params.friction_angle);
    double cos_phi = std::cos(params.friction_angle);
    
    // Alternative: inner cone matching (conservative)
    alpha = 6.0 * sin_phi / (std::sqrt(3.0) * (3.0 - sin_phi));
    kappa = 6.0 * params.cohesion * cos_phi / (std::sqrt(3.0) * (3.0 - sin_phi));
}

double DruckerPragerModel::yieldFunction(const std::array<double, 6>& stress,
                                         double hardening_variable) const {
    // F = √J2 + α*I1 - κ ≤ 0
    double I1 = computeI1(stress);
    double J2 = computeJ2(stress);
    double sqrtJ2 = std::sqrt(std::max(0.0, J2));
    
    // Current yield strength with hardening
    double kappa_h = kappa + hardeningFunction(hardening_variable);
    
    return sqrtJ2 + alpha * I1 - kappa_h;
}

void DruckerPragerModel::flowDirection(const std::array<double, 6>& stress,
                                       std::array<double, 6>& flow_dir) const {
    // For non-associative flow, use dilation angle instead of friction angle
    double sin_psi = std::sin(params.dilation_angle);
    double alpha_g = 6.0 * sin_psi / (std::sqrt(3.0) * (3.0 - sin_psi));
    
    // Flow direction: ∂G/∂σ = s/(2√J2) + α_g * I
    double J2 = computeJ2(stress);
    double sqrtJ2 = std::sqrt(std::max(1e-20, J2));
    
    std::array<double, 6> dev;
    computeDeviatoricStress(stress, dev);
    
    // Normal to deviatoric part
    for (int i = 0; i < 3; ++i) {
        flow_dir[i] = dev[i] / (2.0 * sqrtJ2) + alpha_g / 3.0;
    }
    for (int i = 3; i < 6; ++i) {
        flow_dir[i] = dev[i] / (2.0 * sqrtJ2);
    }
    
    // Normalize
    double norm = tensorNorm(flow_dir);
    if (norm > 1e-15) {
        for (int i = 0; i < 6; ++i) {
            flow_dir[i] /= norm;
        }
    }
}

double DruckerPragerModel::hardeningFunction(double accumulated_plastic_strain) const {
    double eps_p = accumulated_plastic_strain;
    
    if (params.hardening_modulus > 0) {
        // Linear isotropic hardening
        double h = params.hardening_modulus * eps_p;
        
        // Saturation
        double max_hardening = params.saturation_yield - kappa;
        return std::min(h, max_hardening);
    }
    else if (params.softening_modulus < 0) {
        // Strain softening
        double softening = -params.softening_modulus * eps_p;
        
        // Minimum strength
        double residual_kappa = params.residual_cohesion * 
                               6.0 * std::cos(params.residual_friction) /
                               (std::sqrt(3.0) * (3.0 - std::sin(params.residual_friction)));
        double max_softening = kappa - residual_kappa;
        
        return -std::min(softening, max_softening);
    }
    
    return 0.0;  // Perfect plasticity
}

void DruckerPragerModel::returnMapping(std::array<double, 6>& stress,
                                       std::array<double, 6>& plastic_strain,
                                       double& accumulated_plastic_strain,
                                       const std::array<double, 6>& strain_increment,
                                       const std::array<std::array<double, 6>, 6>& C,
                                       bool& plastic_loading) {
    // Elastic predictor
    std::array<double, 6> delta_stress;
    applyElasticModulus(C, strain_increment, delta_stress);
    
    std::array<double, 6> trial_stress;
    for (int i = 0; i < 6; ++i) {
        trial_stress[i] = stress[i] + delta_stress[i];
    }
    
    // Check yield condition
    double F_trial = yieldFunction(trial_stress, accumulated_plastic_strain);
    
    if (F_trial <= 0.0) {
        // Elastic step
        stress = trial_stress;
        plastic_loading = false;
        return;
    }
    
    // Plastic corrector - return mapping
    plastic_loading = true;
    
    // Get elastic moduli (assume isotropic)
    double K = (C[0][0] + C[0][1] + C[0][2]) / 3.0;  // Bulk modulus
    double G = (C[0][0] - C[0][1]) / 2.0;           // Shear modulus
    
    // Trial invariants
    double I1_trial = computeI1(trial_stress);
    double J2_trial = computeJ2(trial_stress);
    double sqrtJ2_trial = std::sqrt(std::max(1e-20, J2_trial));
    
    // Flow rule parameters
    double sin_psi = std::sin(params.dilation_angle);
    double alpha_g = 6.0 * sin_psi / (std::sqrt(3.0) * (3.0 - sin_psi));
    
    // Solve for plastic multiplier using Newton-Raphson
    // F = sqrtJ2 + α*I1 - κ = 0
    // sqrtJ2_n+1 = sqrtJ2_trial - G*Δγ
    // I1_n+1 = I1_trial - 9K*α_g*Δγ
    
    double delta_gamma = 0.0;
    double residual = F_trial;
    int max_iter = 20;
    double tol = 1e-10;
    
    for (int iter = 0; iter < max_iter; ++iter) {
        // Current values
        double sqrtJ2 = sqrtJ2_trial - G * delta_gamma;
        double I1 = I1_trial - 9.0 * K * alpha_g * delta_gamma;
        double eps_p_new = accumulated_plastic_strain + delta_gamma;
        
        // Yield function
        double kappa_h = kappa + hardeningFunction(eps_p_new);
        residual = sqrtJ2 + alpha * I1 - kappa_h;
        
        if (std::abs(residual) < tol) break;
        
        // Derivative of F w.r.t. Δγ
        double dF_dgamma = -G - 9.0 * K * alpha * alpha_g - params.hardening_modulus;
        
        // Newton update
        double d_gamma = -residual / dF_dgamma;
        delta_gamma += d_gamma;
        delta_gamma = std::max(0.0, delta_gamma);
    }
    
    // Update stress
    double sqrtJ2 = sqrtJ2_trial - G * delta_gamma;
    double I1 = I1_trial - 9.0 * K * alpha_g * delta_gamma;
    
    // Reconstruct stress from invariants
    double p = I1 / 3.0;
    
    // Deviatoric stress scaled
    std::array<double, 6> s_trial;
    computeDeviatoricStress(trial_stress, s_trial);
    
    double scale = (sqrtJ2_trial > 1e-15) ? sqrtJ2 / sqrtJ2_trial : 0.0;
    
    stress[0] = scale * s_trial[0] + p;
    stress[1] = scale * s_trial[1] + p;
    stress[2] = scale * s_trial[2] + p;
    stress[3] = scale * s_trial[3];
    stress[4] = scale * s_trial[4];
    stress[5] = scale * s_trial[5];
    
    // Update plastic strain
    std::array<double, 6> flow_dir;
    flowDirection(stress, flow_dir);
    
    for (int i = 0; i < 6; ++i) {
        plastic_strain[i] += delta_gamma * flow_dir[i];
    }
    
    accumulated_plastic_strain += delta_gamma;
}

// =============================================================================
// VonMisesModel Implementation
// =============================================================================

VonMisesModel::VonMisesModel() {}

void VonMisesModel::setParameters(const Parameters& p) {
    params = p;
}

double VonMisesModel::yieldFunction(const std::array<double, 6>& stress,
                                    double hardening_variable) const {
    // F = √(3*J2) - σ_y(κ) ≤ 0
    double J2 = computeJ2(stress);
    double sigma_eq = std::sqrt(3.0 * J2);
    
    // Current yield stress with hardening
    double sigma_y = params.yield_stress + params.hardening_modulus * hardening_variable;
    sigma_y = std::min(sigma_y, params.saturation_stress);
    
    return sigma_eq - sigma_y;
}

void VonMisesModel::returnMapping(std::array<double, 6>& stress,
                                  std::array<double, 6>& plastic_strain,
                                  double& accumulated_plastic_strain,
                                  const std::array<double, 6>& strain_increment,
                                  const std::array<std::array<double, 6>, 6>& C,
                                  bool& plastic_loading) {
    // Elastic predictor
    std::array<double, 6> delta_stress;
    applyElasticModulus(C, strain_increment, delta_stress);
    
    std::array<double, 6> trial_stress;
    for (int i = 0; i < 6; ++i) {
        trial_stress[i] = stress[i] + delta_stress[i];
    }
    
    // Check yield
    double F_trial = yieldFunction(trial_stress, accumulated_plastic_strain);
    
    if (F_trial <= 0.0) {
        stress = trial_stress;
        plastic_loading = false;
        return;
    }
    
    plastic_loading = true;
    
    // Von Mises return mapping (radial return)
    double G = (C[0][0] - C[0][1]) / 2.0;  // Shear modulus
    
    // Trial deviatoric stress
    std::array<double, 6> s_trial;
    computeDeviatoricStress(trial_stress, s_trial);
    
    double J2_trial = computeJ2(trial_stress);
    double sigma_eq_trial = std::sqrt(3.0 * J2_trial);
    
    // Current yield stress
    double sigma_y = params.yield_stress + 
                    params.hardening_modulus * accumulated_plastic_strain;
    sigma_y = std::min(sigma_y, params.saturation_stress);
    
    // Plastic multiplier
    double delta_gamma = (sigma_eq_trial - sigma_y) / 
                        (3.0 * G + params.hardening_modulus);
    
    // Update deviatoric stress (radial return)
    double scale = sigma_y / sigma_eq_trial;
    
    double p = computeMeanStress(trial_stress);
    
    stress[0] = scale * s_trial[0] + p;
    stress[1] = scale * s_trial[1] + p;
    stress[2] = scale * s_trial[2] + p;
    stress[3] = scale * s_trial[3];
    stress[4] = scale * s_trial[4];
    stress[5] = scale * s_trial[5];
    
    // Update plastic strain
    double norm_s = std::sqrt(2.0 * J2_trial);
    if (norm_s > 1e-15) {
        for (int i = 0; i < 6; ++i) {
            plastic_strain[i] += delta_gamma * s_trial[i] / norm_s;
        }
    }
    
    accumulated_plastic_strain += delta_gamma;
}

// =============================================================================
// MohrCoulombModel Implementation
// =============================================================================

MohrCoulombModel::MohrCoulombModel() {}

void MohrCoulombModel::setParameters(const Parameters& p) {
    params = p;
}

double MohrCoulombModel::yieldFunction(const std::array<double, 6>& stress,
                                       double hardening_variable) const {
    // Mohr-Coulomb in principal stress form:
    // F = (σ1 - σ3)/2 + (σ1 + σ3)/2 * sin(φ) - c * cos(φ) ≤ 0
    
    double sigma1, sigma2, sigma3;
    computePrincipalStresses(stress, sigma1, sigma2, sigma3);
    
    double sin_phi = std::sin(params.friction_angle);
    double cos_phi = std::cos(params.friction_angle);
    
    return (sigma1 - sigma3) / 2.0 + 
           (sigma1 + sigma3) / 2.0 * sin_phi - 
           params.cohesion * cos_phi;
}

void MohrCoulombModel::returnMapping(std::array<double, 6>& stress,
                                     std::array<double, 6>& plastic_strain,
                                     double& accumulated_plastic_strain,
                                     const std::array<double, 6>& strain_increment,
                                     const std::array<std::array<double, 6>, 6>& C,
                                     bool& plastic_loading) {
    // Elastic predictor
    std::array<double, 6> delta_stress;
    applyElasticModulus(C, strain_increment, delta_stress);
    
    std::array<double, 6> trial_stress;
    for (int i = 0; i < 6; ++i) {
        trial_stress[i] = stress[i] + delta_stress[i];
    }
    
    // Check yield
    double F_trial = yieldFunction(trial_stress, accumulated_plastic_strain);
    
    if (F_trial <= 0.0) {
        stress = trial_stress;
        plastic_loading = false;
        return;
    }
    
    plastic_loading = true;
    
    // Mohr-Coulomb return mapping in principal stress space
    double sigma1_t, sigma2_t, sigma3_t;
    computePrincipalStresses(trial_stress, sigma1_t, sigma2_t, sigma3_t);
    
    double sin_phi = std::sin(params.friction_angle);
    double cos_phi = std::cos(params.friction_angle);
    double sin_psi = std::sin(params.dilation_angle);
    
    // Elastic moduli
    double K = (C[0][0] + C[0][1] + C[0][2]) / 3.0;
    double G = (C[0][0] - C[0][1]) / 2.0;
    
    // Return to main plane (σ1 - σ3 plane)
    double a = 1.0 + sin_phi;
    double b = 1.0 - sin_phi;
    double a_g = 1.0 + sin_psi;
    double b_g = 1.0 - sin_psi;
    
    // Plastic multiplier
    double denom = 4.0 * G + (a*a_g + b*b_g) * K / 3.0;
    double delta_gamma = F_trial / denom;
    
    // Updated principal stresses
    double sigma1 = sigma1_t - 2.0*G*(a_g/2.0)*delta_gamma - K*a_g*delta_gamma/3.0;
    double sigma3 = sigma3_t + 2.0*G*(b_g/2.0)*delta_gamma - K*b_g*delta_gamma/3.0;
    double sigma2 = sigma2_t - K*(a_g - b_g)*delta_gamma/3.0;
    
    // Check for return to edge/apex (simplified)
    double c_cos_phi = params.cohesion * cos_phi;
    
    // Tension cutoff
    if (sigma1 > params.tension_cutoff) {
        sigma1 = params.tension_cutoff;
    }
    
    // Reconstruct stress tensor (simplified - assumes principal directions unchanged)
    // Full implementation would track principal directions
    double p = (sigma1 + sigma2 + sigma3) / 3.0;
    double scale = std::sqrt(2.0 * computeJ2(trial_stress));
    
    std::array<double, 6> s_trial;
    computeDeviatoricStress(trial_stress, s_trial);
    
    // Approximate reconstruction
    double new_J2 = ((sigma1 - p)*(sigma1 - p) + 
                    (sigma2 - p)*(sigma2 - p) + 
                    (sigma3 - p)*(sigma3 - p)) / 2.0;
    
    double ratio = (scale > 1e-15) ? std::sqrt(new_J2 / computeJ2(trial_stress)) : 0.0;
    
    stress[0] = ratio * s_trial[0] + p;
    stress[1] = ratio * s_trial[1] + p;
    stress[2] = ratio * s_trial[2] + p;
    stress[3] = ratio * s_trial[3];
    stress[4] = ratio * s_trial[4];
    stress[5] = ratio * s_trial[5];
    
    accumulated_plastic_strain += delta_gamma;
}

// =============================================================================
// CapModel Implementation
// =============================================================================

CapModel::CapModel() {}

void CapModel::setParameters(const Parameters& p) {
    params = p;
}

double CapModel::shearYield(const std::array<double, 6>& stress, double kappa) const {
    // Drucker-Prager type shear surface
    double sin_phi = std::sin(params.friction_angle);
    double cos_phi = std::cos(params.friction_angle);
    
    double alpha = 6.0 * sin_phi / (std::sqrt(3.0) * (3.0 - sin_phi));
    double k = 6.0 * params.cohesion * cos_phi / (std::sqrt(3.0) * (3.0 - sin_phi));
    
    double I1 = computeI1(stress);
    double J2 = computeJ2(stress);
    
    return std::sqrt(J2) + alpha * I1 - k;
}

double CapModel::capYield(const std::array<double, 6>& stress, double X) const {
    // Elliptical cap surface
    double I1 = computeI1(stress);
    double J2 = computeJ2(stress);
    
    double R = params.cap_eccentricity;
    
    // Cap center
    double L = params.cohesion / std::tan(params.friction_angle);  // Intersection with I1 axis
    
    // Cap function: √J2² + ((I1 - L)/R)² - (X - L) ≤ 0
    double f = std::sqrt(J2 + std::pow((I1 - L) / R, 2)) - (X - L) / R;
    
    return f;
}

void CapModel::returnMapping(std::array<double, 6>& stress,
                             std::array<double, 6>& plastic_strain,
                             double& accumulated_plastic_strain,
                             double& cap_position,
                             const std::array<double, 6>& strain_increment,
                             const std::array<std::array<double, 6>, 6>& C,
                             bool& plastic_loading) {
    // Elastic predictor
    std::array<double, 6> delta_stress;
    applyElasticModulus(C, strain_increment, delta_stress);
    
    std::array<double, 6> trial_stress;
    for (int i = 0; i < 6; ++i) {
        trial_stress[i] = stress[i] + delta_stress[i];
    }
    
    // Check both yield surfaces
    double F_shear = shearYield(trial_stress, 0);
    double F_cap = capYield(trial_stress, cap_position);
    
    if (F_shear <= 0.0 && F_cap <= 0.0) {
        stress = trial_stress;
        plastic_loading = false;
        return;
    }
    
    plastic_loading = true;
    
    // Determine which surface is active
    double G = (C[0][0] - C[0][1]) / 2.0;
    double K = (C[0][0] + C[0][1] + C[0][2]) / 3.0;
    
    if (F_shear > 0.0 && F_cap <= 0.0) {
        // Return to shear surface (similar to Drucker-Prager)
        double sin_phi = std::sin(params.friction_angle);
        double sin_psi = std::sin(params.friction_angle);  // Assume associative for cap
        
        double alpha = 6.0 * sin_phi / (std::sqrt(3.0) * (3.0 - sin_phi));
        double alpha_g = 6.0 * sin_psi / (std::sqrt(3.0) * (3.0 - sin_psi));
        
        double I1_t = computeI1(trial_stress);
        double J2_t = computeJ2(trial_stress);
        double sqrtJ2_t = std::sqrt(std::max(1e-20, J2_t));
        
        double delta_gamma = F_shear / (G + 9.0 * K * alpha * alpha_g);
        
        double sqrtJ2 = sqrtJ2_t - G * delta_gamma;
        double I1 = I1_t - 9.0 * K * alpha_g * delta_gamma;
        
        double p = I1 / 3.0;
        double scale = (sqrtJ2_t > 1e-15) ? sqrtJ2 / sqrtJ2_t : 0.0;
        
        std::array<double, 6> s_t;
        computeDeviatoricStress(trial_stress, s_t);
        
        stress[0] = scale * s_t[0] + p;
        stress[1] = scale * s_t[1] + p;
        stress[2] = scale * s_t[2] + p;
        stress[3] = scale * s_t[3];
        stress[4] = scale * s_t[4];
        stress[5] = scale * s_t[5];
        
        accumulated_plastic_strain += delta_gamma;
    }
    else if (F_cap > 0.0) {
        // Return to cap surface
        double I1_t = computeI1(trial_stress);
        double J2_t = computeJ2(trial_stress);
        
        double R = params.cap_eccentricity;
        double L = params.cohesion / std::tan(params.friction_angle);
        
        // Simplified cap return (radial in I1-√J2 space)
        double dist = std::sqrt(J2_t + std::pow((I1_t - L) / R, 2));
        double cap_radius = (cap_position - L) / R;
        
        if (dist > 1e-15) {
            double scale = cap_radius / dist;
            
            double sqrtJ2 = std::sqrt(J2_t) * scale;
            double I1 = L + (I1_t - L) * scale;
            
            double p = I1 / 3.0;
            double dev_scale = (J2_t > 1e-20) ? sqrtJ2 / std::sqrt(J2_t) : 0.0;
            
            std::array<double, 6> s_t;
            computeDeviatoricStress(trial_stress, s_t);
            
            stress[0] = dev_scale * s_t[0] + p;
            stress[1] = dev_scale * s_t[1] + p;
            stress[2] = dev_scale * s_t[2] + p;
            stress[3] = dev_scale * s_t[3];
            stress[4] = dev_scale * s_t[4];
            stress[5] = dev_scale * s_t[5];
        } else {
            stress = trial_stress;
        }
        
        // Harden cap position
        double eps_v_p = accumulated_plastic_strain;  // Volumetric plastic strain
        double delta_X = params.cap_hardening_modulus * (1.0 - cap_position / params.initial_cap_position);
        cap_position += delta_X * 0.01;  // Increment
    }
}

// =============================================================================
// PlasticityModel (Unified Interface) Implementation
// =============================================================================

PlasticityModel::PlasticityModel(YieldCriterion c)
    : enabled(true), criterion(c), 
      flow_rule(FlowRule::NON_ASSOCIATIVE),
      hardening_law(HardeningLaw::NONE) {
    
    drucker_prager = std::make_unique<DruckerPragerModel>();
    von_mises = std::make_unique<VonMisesModel>();
    mohr_coulomb = std::make_unique<MohrCoulombModel>();
    cap_model = std::make_unique<CapModel>();
}

void PlasticityModel::setYieldCriterion(YieldCriterion c) {
    criterion = c;
}

void PlasticityModel::setFlowRule(FlowRule rule) {
    flow_rule = rule;
}

void PlasticityModel::setHardeningLaw(HardeningLaw law) {
    hardening_law = law;
}

void PlasticityModel::configure(const std::map<std::string, std::string>& config) {
    enabled = parseString(config, "plasticity_enabled", "false") == "true";
    
    std::string crit_str = parseString(config, "yield_criterion", "drucker_prager");
    if (crit_str == "von_mises") criterion = YieldCriterion::VON_MISES;
    else if (crit_str == "mohr_coulomb") criterion = YieldCriterion::MOHR_COULOMB;
    else if (crit_str == "cap") criterion = YieldCriterion::CAP_MODEL;
    else criterion = YieldCriterion::DRUCKER_PRAGER;
    
    // Configure Drucker-Prager
    DruckerPragerModel::Parameters dp_params;
    dp_params.friction_angle = parseDouble(config, "friction_angle", 30.0) * M_PI / 180.0;
    dp_params.cohesion = parseDouble(config, "cohesion", 1e6);
    dp_params.dilation_angle = parseDouble(config, "dilation_angle", 10.0) * M_PI / 180.0;
    dp_params.tensile_strength = parseDouble(config, "tensile_strength", 0.1e6);
    dp_params.hardening_modulus = parseDouble(config, "hardening_modulus", 0.0);
    dp_params.softening_modulus = parseDouble(config, "softening_modulus", 0.0);
    drucker_prager->setParameters(dp_params);
    
    // Configure von Mises
    VonMisesModel::Parameters vm_params;
    vm_params.yield_stress = parseDouble(config, "yield_stress", 10e6);
    vm_params.hardening_modulus = parseDouble(config, "hardening_modulus", 0.0);
    von_mises->setParameters(vm_params);
    
    // Configure Mohr-Coulomb
    MohrCoulombModel::Parameters mc_params;
    mc_params.friction_angle = dp_params.friction_angle;
    mc_params.cohesion = dp_params.cohesion;
    mc_params.dilation_angle = dp_params.dilation_angle;
    mohr_coulomb->setParameters(mc_params);
}

void PlasticityModel::setDruckerPragerParameters(const DruckerPragerModel::Parameters& params) {
    drucker_prager->setParameters(params);
}

void PlasticityModel::setVonMisesParameters(const VonMisesModel::Parameters& params) {
    von_mises->setParameters(params);
}

void PlasticityModel::setMohrCoulombParameters(const MohrCoulombModel::Parameters& params) {
    mohr_coulomb->setParameters(params);
}

void PlasticityModel::setCapModelParameters(const CapModel::Parameters& params) {
    cap_model->setParameters(params);
}

void PlasticityModel::integrateStress(std::array<double, 6>& stress,
                                      const std::array<double, 6>& strain_increment,
                                      const std::array<std::array<double, 6>, 6>& elastic_modulus,
                                      PlasticityState& state) {
    if (!enabled) {
        // Pure elastic
        std::array<double, 6> delta_stress;
        applyElasticModulus(elastic_modulus, strain_increment, delta_stress);
        for (int i = 0; i < 6; ++i) {
            stress[i] += delta_stress[i];
        }
        state.is_plastic = false;
        return;
    }
    
    bool plastic_loading = false;
    
    switch (criterion) {
        case YieldCriterion::VON_MISES:
            von_mises->returnMapping(stress, state.plastic_strain,
                                    state.accumulated_plastic_strain,
                                    strain_increment, elastic_modulus,
                                    plastic_loading);
            break;
            
        case YieldCriterion::MOHR_COULOMB:
            mohr_coulomb->returnMapping(stress, state.plastic_strain,
                                       state.accumulated_plastic_strain,
                                       strain_increment, elastic_modulus,
                                       plastic_loading);
            break;
            
        case YieldCriterion::CAP_MODEL:
            cap_model->returnMapping(stress, state.plastic_strain,
                                    state.accumulated_plastic_strain,
                                    state.cap_position,
                                    strain_increment, elastic_modulus,
                                    plastic_loading);
            break;
            
        case YieldCriterion::DRUCKER_PRAGER:
        default:
            drucker_prager->returnMapping(stress, state.plastic_strain,
                                         state.accumulated_plastic_strain,
                                         strain_increment, elastic_modulus,
                                         plastic_loading);
            break;
    }
    
    state.is_plastic = plastic_loading;
    state.yield_function_value = getYieldFunctionValue(stress, state);
}

bool PlasticityModel::isYielding(const std::array<double, 6>& stress,
                                 const PlasticityState& state) const {
    return getYieldFunctionValue(stress, state) >= 0.0;
}

double PlasticityModel::getYieldFunctionValue(const std::array<double, 6>& stress,
                                              const PlasticityState& state) const {
    switch (criterion) {
        case YieldCriterion::VON_MISES:
            return von_mises->yieldFunction(stress, state.accumulated_plastic_strain);
            
        case YieldCriterion::MOHR_COULOMB:
            return mohr_coulomb->yieldFunction(stress, state.accumulated_plastic_strain);
            
        case YieldCriterion::CAP_MODEL:
            return std::max(cap_model->shearYield(stress, 0),
                          cap_model->capYield(stress, state.cap_position));
            
        case YieldCriterion::DRUCKER_PRAGER:
        default:
            return drucker_prager->yieldFunction(stress, state.accumulated_plastic_strain);
    }
}

void PlasticityModel::getConsistentTangent(const std::array<double, 6>& stress,
                                           const PlasticityState& state,
                                           const std::array<std::array<double, 6>, 6>& elastic_modulus,
                                           std::array<std::array<double, 6>, 6>& tangent_modulus) {
    if (!state.is_plastic) {
        // Elastic tangent
        tangent_modulus = elastic_modulus;
        return;
    }
    
    // Algorithmic consistent tangent for Drucker-Prager (simplified)
    // C^alg = C^e - (C^e : n ⊗ n : C^e) / (n : C^e : n + H)
    
    // Get flow direction
    std::array<double, 6> n;
    drucker_prager->flowDirection(stress, n);
    
    // Compute C^e : n
    std::array<double, 6> Ce_n;
    for (int i = 0; i < 6; ++i) {
        Ce_n[i] = 0.0;
        for (int j = 0; j < 6; ++j) {
            Ce_n[i] += elastic_modulus[i][j] * n[j];
        }
    }
    
    // Compute n : C^e : n
    double n_Ce_n = 0.0;
    for (int i = 0; i < 6; ++i) {
        n_Ce_n += n[i] * Ce_n[i];
    }
    
    // Add hardening modulus
    double H = drucker_prager->getParameters().hardening_modulus;
    double denom = n_Ce_n + H;
    
    if (std::abs(denom) < 1e-15) {
        tangent_modulus = elastic_modulus;
        return;
    }
    
    // C^alg = C^e - (C^e : n) ⊗ (n : C^e) / denom
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            tangent_modulus[i][j] = elastic_modulus[i][j] - Ce_n[i] * Ce_n[j] / denom;
        }
    }
}

// =============================================================================
// DamageModel Implementation
// =============================================================================

DamageModel::DamageModel() {}

void DamageModel::setParameters(const Parameters& p) {
    params = p;
}

void DamageModel::updateDamage(double& damage,
                               const std::array<double, 6>& stress,
                               const std::array<double, 6>& strain,
                               double dt) {
    // Scalar damage model
    // Damage variable D ∈ [0, D_max]
    // Effective stress: σ_eff = σ / (1 - D)
    
    double strain_equiv = tensorNorm(strain);
    
    if (strain_equiv <= params.damage_threshold) {
        // Below threshold - possible healing
        if (params.healing_rate > 0) {
            damage -= params.healing_rate * dt;
            damage = std::max(0.0, damage);
        }
        return;
    }
    
    // Damage evolution
    double damage_rate = params.damage_exponent * 
                        std::pow(strain_equiv - params.damage_threshold, 
                                params.damage_exponent - 1);
    
    damage += damage_rate * dt;
    damage = std::min(damage, params.maximum_damage);
}

void DamageModel::applyDamage(std::array<std::array<double, 6>, 6>& elastic_modulus,
                              double damage) const {
    double reduction = 1.0 - damage;
    
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            elastic_modulus[i][j] *= reduction;
        }
    }
}

// =============================================================================
// PlasticityConfig Implementation
// =============================================================================

void PlasticityConfig::parseConfig(const std::map<std::string, std::string>& config) {
    enabled = parseString(config, "plasticity_enabled", "false") == "true";
    yield_criterion = parseString(config, "yield_criterion", "drucker_prager");
    flow_rule = parseString(config, "flow_rule", "non_associative");
    hardening_law = parseString(config, "hardening_law", "none");
    
    friction_angle = parseDouble(config, "friction_angle", 30.0);
    cohesion = parseDouble(config, "cohesion", 1e6);
    dilation_angle = parseDouble(config, "dilation_angle", 10.0);
    tensile_strength = parseDouble(config, "tensile_strength", 0.1e6);
    
    hardening_modulus = parseDouble(config, "hardening_modulus", 0.0);
    softening_modulus = parseDouble(config, "softening_modulus", 0.0);
    residual_cohesion = parseDouble(config, "residual_cohesion", 0.1e6);
}

PlasticityModel PlasticityConfig::createModel() const {
    YieldCriterion crit = YieldCriterion::DRUCKER_PRAGER;
    if (yield_criterion == "von_mises") crit = YieldCriterion::VON_MISES;
    else if (yield_criterion == "mohr_coulomb") crit = YieldCriterion::MOHR_COULOMB;
    else if (yield_criterion == "cap") crit = YieldCriterion::CAP_MODEL;
    
    PlasticityModel model(crit);
    model.enable(enabled);
    
    DruckerPragerModel::Parameters dp_params;
    dp_params.friction_angle = friction_angle * M_PI / 180.0;
    dp_params.cohesion = cohesion;
    dp_params.dilation_angle = dilation_angle * M_PI / 180.0;
    dp_params.tensile_strength = tensile_strength;
    dp_params.hardening_modulus = hardening_modulus;
    dp_params.softening_modulus = softening_modulus;
    dp_params.residual_cohesion = residual_cohesion;
    model.setDruckerPragerParameters(dp_params);
    
    return model;
}

} // namespace FSRM
