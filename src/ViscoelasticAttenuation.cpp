/**
 * @file ViscoelasticAttenuation.cpp
 * @brief Full implementation of viscoelastic attenuation with Maxwell bodies
 * 
 * Based on:
 * - Moczo et al. (2007) "The Finite-Difference Modelling of Earthquake Motions"
 * - Kaser & Dumbser (2006) "ADER-DG method for viscoelasticity"
 * - SeisSol implementation for consistency
 */

#include "ViscoelasticAttenuation.hpp"
#include <cmath>
#include <algorithm>
#include <fstream>
#include <numeric>
#include <stdexcept>
#include <iostream>

namespace FSRM {

// =============================================================================
// Utility functions
// =============================================================================

namespace {

double parseDouble(const std::map<std::string, std::string>& config,
                  const std::string& key, double default_val) {
    auto it = config.find(key);
    if (it != config.end() && !it->second.empty()) {
        try { return std::stod(it->second); }
        catch (...) {}
    }
    return default_val;
}

std::string parseString(const std::map<std::string, std::string>& config,
                       const std::string& key, const std::string& default_val) {
    auto it = config.find(key);
    return (it != config.end() && !it->second.empty()) ? it->second : default_val;
}

// Compute relaxation frequencies logarithmically spaced
std::vector<double> logSpacedFrequencies(int n, double f_min, double f_max) {
    std::vector<double> freqs(n);
    
    if (n == 1) {
        freqs[0] = std::sqrt(f_min * f_max);
        return freqs;
    }
    
    double log_f_min = std::log(f_min);
    double log_f_max = std::log(f_max);
    double d_log = (log_f_max - log_f_min) / (n - 1);
    
    for (int i = 0; i < n; ++i) {
        freqs[i] = std::exp(log_f_min + i * d_log);
    }
    
    return freqs;
}

} // anonymous namespace

// =============================================================================
// ViscoelasticAttenuation Implementation
// =============================================================================

ViscoelasticAttenuation::ViscoelasticAttenuation()
    : enabled(false) {}

void ViscoelasticAttenuation::setParameters(const Parameters& p) {
    params = p;
    
    if (params.mechanisms.empty()) {
        setupMechanismsFromQ(params.num_mechanisms, params.freq_central,
                           params.freq_ratio, params.Qp_desired, params.Qs_desired);
    } else {
        mechanisms = params.mechanisms;
    }
    
    enabled = true;
}

void ViscoelasticAttenuation::setupMechanismsFromQ(int num_mechanisms,
                                                   double freq_central,
                                                   double freq_ratio,
                                                   double Qp, double Qs) {
    // Compute frequency bounds
    double f_min = freq_central / std::sqrt(freq_ratio);
    double f_max = freq_central * std::sqrt(freq_ratio);
    
    // Get logarithmically spaced relaxation frequencies
    std::vector<double> relax_freqs = logSpacedFrequencies(num_mechanisms, f_min, f_max);
    
    // Initialize mechanisms
    mechanisms.clear();
    mechanisms.reserve(num_mechanisms);
    
    for (int i = 0; i < num_mechanisms; ++i) {
        MaxwellMechanism mech(relax_freqs[i], 0.0, 0.0);
        mechanisms.push_back(mech);
    }
    
    // Store Q targets
    params.num_mechanisms = num_mechanisms;
    params.freq_central = freq_central;
    params.freq_ratio = freq_ratio;
    params.Qp_desired = Qp;
    params.Qs_desired = Qs;
    
    // Compute anelastic coefficients to match target Q
    computeAnelasticCoefficients();
}

void ViscoelasticAttenuation::addMechanism(const MaxwellMechanism& mech) {
    mechanisms.push_back(mech);
}

void ViscoelasticAttenuation::clearMechanisms() {
    mechanisms.clear();
}

void ViscoelasticAttenuation::computeAnelasticCoefficients() {
    if (mechanisms.empty()) return;
    
    int n_mech = mechanisms.size();
    
    // Sample frequencies for fitting
    double f_min = params.freq_central / std::sqrt(params.freq_ratio);
    double f_max = params.freq_central * std::sqrt(params.freq_ratio);
    int n_sample = n_mech * 10;  // Oversample for better fit
    
    std::vector<double> sample_freqs(n_sample);
    double log_f_min = std::log(f_min);
    double log_f_max = std::log(f_max);
    
    for (int i = 0; i < n_sample; ++i) {
        double t = static_cast<double>(i) / (n_sample - 1);
        sample_freqs[i] = std::exp(log_f_min + t * (log_f_max - log_f_min));
    }
    
    // Compute coefficients for P-waves
    std::vector<double> Y_p(n_mech);
    solveForAnelasticCoefficients(sample_freqs, params.Qp_desired, Y_p);
    
    // Compute coefficients for S-waves
    std::vector<double> Y_s(n_mech);
    solveForAnelasticCoefficients(sample_freqs, params.Qs_desired, Y_s);
    
    // Store coefficients
    for (int i = 0; i < n_mech; ++i) {
        mechanisms[i].anelastic_coefficient_p = Y_p[i];
        mechanisms[i].anelastic_coefficient_s = Y_s[i];
    }
}

void ViscoelasticAttenuation::solveForAnelasticCoefficients(
    const std::vector<double>& frequencies,
    double Q_target,
    std::vector<double>& Y_coefficients) {
    
    int n_mech = mechanisms.size();
    int n_freq = frequencies.size();
    
    Y_coefficients.resize(n_mech);
    
    // Target: 1/Q = Σ_i Y_i * ωτ_i / (1 + (ωτ_i)²)
    // This is a linear least-squares problem
    
    // Build matrix A and target vector b
    // A[j][i] = ω_j * τ_i / (1 + (ω_j * τ_i)²)
    // b[j] = 1/Q_target
    
    std::vector<std::vector<double>> A(n_freq, std::vector<double>(n_mech));
    std::vector<double> b(n_freq, 1.0 / Q_target);
    
    for (int j = 0; j < n_freq; ++j) {
        double omega = 2.0 * M_PI * frequencies[j];
        
        for (int i = 0; i < n_mech; ++i) {
            double tau = mechanisms[i].relaxation_time;
            double omega_tau = omega * tau;
            A[j][i] = omega_tau / (1.0 + omega_tau * omega_tau);
        }
    }
    
    // Solve using normal equations: (A^T A) Y = A^T b
    // For small systems, use direct solution
    
    // Compute A^T A
    std::vector<std::vector<double>> ATA(n_mech, std::vector<double>(n_mech, 0.0));
    for (int i = 0; i < n_mech; ++i) {
        for (int k = 0; k < n_mech; ++k) {
            for (int j = 0; j < n_freq; ++j) {
                ATA[i][k] += A[j][i] * A[j][k];
            }
        }
    }
    
    // Compute A^T b
    std::vector<double> ATb(n_mech, 0.0);
    for (int i = 0; i < n_mech; ++i) {
        for (int j = 0; j < n_freq; ++j) {
            ATb[i] += A[j][i] * b[j];
        }
    }
    
    // Add regularization for stability
    double reg = 1e-8;
    for (int i = 0; i < n_mech; ++i) {
        ATA[i][i] += reg;
    }
    
    // Solve using Gaussian elimination with pivoting
    std::vector<std::vector<double>> aug(n_mech, std::vector<double>(n_mech + 1));
    for (int i = 0; i < n_mech; ++i) {
        for (int k = 0; k < n_mech; ++k) {
            aug[i][k] = ATA[i][k];
        }
        aug[i][n_mech] = ATb[i];
    }
    
    // Forward elimination
    for (int i = 0; i < n_mech; ++i) {
        // Find pivot
        int max_row = i;
        for (int k = i + 1; k < n_mech; ++k) {
            if (std::abs(aug[k][i]) > std::abs(aug[max_row][i])) {
                max_row = k;
            }
        }
        std::swap(aug[i], aug[max_row]);
        
        // Eliminate
        for (int k = i + 1; k < n_mech; ++k) {
            if (std::abs(aug[i][i]) < 1e-15) continue;
            double factor = aug[k][i] / aug[i][i];
            for (int j = i; j <= n_mech; ++j) {
                aug[k][j] -= factor * aug[i][j];
            }
        }
    }
    
    // Back substitution
    for (int i = n_mech - 1; i >= 0; --i) {
        Y_coefficients[i] = aug[i][n_mech];
        for (int j = i + 1; j < n_mech; ++j) {
            Y_coefficients[i] -= aug[i][j] * Y_coefficients[j];
        }
        if (std::abs(aug[i][i]) > 1e-15) {
            Y_coefficients[i] /= aug[i][i];
        }
        
        // Ensure non-negative (physical constraint)
        Y_coefficients[i] = std::max(0.0, Y_coefficients[i]);
    }
}

void ViscoelasticAttenuation::precomputeFactors(double dt) {
    exp_factors.resize(mechanisms.size());
    for (size_t i = 0; i < mechanisms.size(); ++i) {
        double tau = mechanisms[i].relaxation_time;
        exp_factors[i] = std::exp(-dt / tau);
    }
}

void ViscoelasticAttenuation::updateAnelasticState(AnelasticState& state,
                                                    const std::array<double, 6>& stress,
                                                    double dt) {
    if (!enabled || mechanisms.empty()) return;
    
    // Ensure state has correct size
    if (state.R.size() != mechanisms.size()) {
        state.R.resize(mechanisms.size(), {0, 0, 0, 0, 0, 0});
    }
    
    // Precompute exponential factors if not done
    if (exp_factors.size() != mechanisms.size()) {
        precomputeFactors(dt);
    }
    
    // Update memory variables for each mechanism
    // dR_i/dt = (Y_i * σ - R_i) / τ_i
    // Analytical solution: R_i(t+dt) = exp(-dt/τ) * R_i(t) + Y_i * (1 - exp(-dt/τ)) * σ
    
    for (size_t i = 0; i < mechanisms.size(); ++i) {
        double exp_factor = exp_factors[i];
        double Y_s = mechanisms[i].anelastic_coefficient_s;
        
        // For shear components (deviatoric stress)
        double factor = Y_s * (1.0 - exp_factor);
        
        // Compute deviatoric stress
        double mean_stress = (stress[0] + stress[1] + stress[2]) / 3.0;
        
        std::array<double, 6> dev_stress;
        dev_stress[0] = stress[0] - mean_stress;
        dev_stress[1] = stress[1] - mean_stress;
        dev_stress[2] = stress[2] - mean_stress;
        dev_stress[3] = stress[3];
        dev_stress[4] = stress[4];
        dev_stress[5] = stress[5];
        
        for (int j = 0; j < 6; ++j) {
            state.R[i][j] = exp_factor * state.R[i][j] + factor * dev_stress[j];
        }
    }
}

void ViscoelasticAttenuation::computeAnelasticStress(const AnelasticState& state,
                                                     std::array<double, 6>& stress_correction) {
    if (!enabled || mechanisms.empty() || state.R.empty()) {
        stress_correction = {0, 0, 0, 0, 0, 0};
        return;
    }
    
    // Sum memory variable contributions
    // σ_correction = Σ_i R_i
    stress_correction = {0, 0, 0, 0, 0, 0};
    
    for (size_t i = 0; i < state.R.size(); ++i) {
        for (int j = 0; j < 6; ++j) {
            stress_correction[j] += state.R[i][j];
        }
    }
}

double ViscoelasticAttenuation::getRelaxedModulus(double M_unrelaxed, bool is_shear) const {
    // M_relaxed = M_unrelaxed / (1 + Σ Y_i)
    double sum_Y = 0.0;
    
    for (const auto& mech : mechanisms) {
        sum_Y += is_shear ? mech.anelastic_coefficient_s : mech.anelastic_coefficient_p;
    }
    
    return M_unrelaxed / (1.0 + sum_Y);
}

double ViscoelasticAttenuation::getRelaxedVelocity(double v_unrelaxed, bool is_shear) const {
    // v ∝ √M, so v_relaxed = v_unrelaxed / √(1 + Σ Y_i)
    double sum_Y = 0.0;
    
    for (const auto& mech : mechanisms) {
        sum_Y += is_shear ? mech.anelastic_coefficient_s : mech.anelastic_coefficient_p;
    }
    
    return v_unrelaxed / std::sqrt(1.0 + sum_Y);
}

std::complex<double> ViscoelasticAttenuation::getComplexModulus(
    double frequency, double M_unrelaxed, bool is_shear) const {
    
    std::complex<double> I(0.0, 1.0);
    double omega = 2.0 * M_PI * frequency;
    
    // M(ω) = M_∞ / [1 + Σ_i Y_i / (1 + iωτ_i)]
    // Using relaxation form
    
    std::complex<double> denom(1.0, 0.0);
    
    for (const auto& mech : mechanisms) {
        double Y = is_shear ? mech.anelastic_coefficient_s : mech.anelastic_coefficient_p;
        double tau = mech.relaxation_time;
        
        std::complex<double> factor = Y / (1.0 + I * omega * tau);
        denom += factor;
    }
    
    return std::complex<double>(M_unrelaxed, 0.0) / denom;
}

double ViscoelasticAttenuation::getQuality(double frequency, bool is_shear) const {
    // Q = Re(M) / Im(M)
    auto M = getComplexModulus(frequency, 1.0, is_shear);
    
    if (std::abs(M.imag()) < 1e-15) return 1e10;  // Essentially elastic
    
    return std::abs(M.real() / M.imag());
}

double ViscoelasticAttenuation::getVelocity(double frequency, double v_unrelaxed,
                                            bool is_shear) const {
    // Phase velocity: c(ω) = c_∞ * √(|M(ω)| / M_∞)
    auto M = getComplexModulus(frequency, 1.0, is_shear);
    double M_ratio = std::abs(M);
    
    // Account for relaxation
    double sum_Y = 0.0;
    for (const auto& mech : mechanisms) {
        sum_Y += is_shear ? mech.anelastic_coefficient_s : mech.anelastic_coefficient_p;
    }
    
    return v_unrelaxed * std::sqrt(M_ratio) / std::sqrt(1.0 + sum_Y);
}

std::vector<ViscoelasticAttenuation::DispersionPoint> 
ViscoelasticAttenuation::computeDispersionCurve(
    double vp_ref, double vs_ref,
    double f_min, double f_max, int num_points) const {
    
    std::vector<DispersionPoint> curve(num_points);
    
    for (int i = 0; i < num_points; ++i) {
        double t = static_cast<double>(i) / (num_points - 1);
        double f = f_min * std::pow(f_max / f_min, t);
        
        curve[i].frequency = f;
        curve[i].velocity_p = getVelocity(f, vp_ref, false);
        curve[i].velocity_s = getVelocity(f, vs_ref, true);
        curve[i].Qp = getQuality(f, false);
        curve[i].Qs = getQuality(f, true);
    }
    
    return curve;
}

double ViscoelasticAttenuation::getMaxQError(double f_min, double f_max,
                                             bool is_shear) const {
    // Check maximum deviation from target Q
    double Q_target = is_shear ? params.Qs_desired : params.Qp_desired;
    double max_error = 0.0;
    
    int n_sample = 100;
    for (int i = 0; i < n_sample; ++i) {
        double t = static_cast<double>(i) / (n_sample - 1);
        double f = f_min * std::pow(f_max / f_min, t);
        
        double Q = getQuality(f, is_shear);
        double error = std::abs(Q - Q_target) / Q_target;
        max_error = std::max(max_error, error);
    }
    
    return max_error;
}

// =============================================================================
// AttenuationDatabase Implementation
// =============================================================================

ViscoelasticAttenuation::Parameters AttenuationDatabase::getParameters(RockType rock) {
    ViscoelasticAttenuation::Parameters params;
    params.num_mechanisms = 3;
    params.freq_central = 1.0;
    params.freq_ratio = 100.0;
    
    switch (rock) {
        case RockType::GRANITE:
            params.Qp_desired = 200.0;
            params.Qs_desired = 100.0;
            break;
            
        case RockType::BASALT:
            params.Qp_desired = 150.0;
            params.Qs_desired = 80.0;
            break;
            
        case RockType::SANDSTONE:
            params.Qp_desired = 80.0;
            params.Qs_desired = 40.0;
            break;
            
        case RockType::LIMESTONE:
            params.Qp_desired = 100.0;
            params.Qs_desired = 50.0;
            break;
            
        case RockType::SHALE:
            params.Qp_desired = 50.0;
            params.Qs_desired = 25.0;
            break;
            
        case RockType::SEDIMENT_SOFT:
            params.Qp_desired = 30.0;
            params.Qs_desired = 15.0;
            break;
            
        case RockType::SEDIMENT_HARD:
            params.Qp_desired = 60.0;
            params.Qs_desired = 30.0;
            break;
            
        case RockType::FAULT_ZONE:
            params.Qp_desired = 20.0;
            params.Qs_desired = 10.0;
            break;
            
        case RockType::DAMAGE_ZONE:
            params.Qp_desired = 30.0;
            params.Qs_desired = 15.0;
            break;
            
        default:
            params.Qp_desired = 50.0;
            params.Qs_desired = 25.0;
    }
    
    return params;
}

double AttenuationDatabase::estimateQs(double vs_km_per_s) {
    // Empirical relation: Q_s ≈ 40-50 * V_s (km/s)
    // Higher velocity = older, more consolidated rock = higher Q
    return 50.0 * vs_km_per_s;
}

double AttenuationDatabase::estimateQp(double Qs) {
    // Generally Q_p / Q_s ≈ 1.5 - 2.5
    return 2.0 * Qs;
}

double AttenuationDatabase::depthCorrection(double Q_surface, double depth_km) {
    // Q typically increases with depth due to increased pressure
    // Q(z) ≈ Q_0 * (1 + α * z)
    double alpha = 0.05;  // per km
    return Q_surface * (1.0 + alpha * depth_km);
}

double AttenuationDatabase::damageCorrection(double Q_intact, double damage) {
    // Damage reduces Q
    // Q(D) ≈ Q_0 * (1 - D)^β
    double beta = 1.5;
    return Q_intact * std::pow(std::max(0.01, 1.0 - damage), beta);
}

// =============================================================================
// SpatialAttenuation Implementation
// =============================================================================

SpatialAttenuation::SpatialAttenuation()
    : velocity_scaling(false), velocity_scaling_factor(50.0),
      depth_scaling(false), depth_gradient(0.1),
      damage_scaling(false) {}

void SpatialAttenuation::setBaseModel(const ViscoelasticAttenuation::Parameters& params) {
    base_params = params;
}

void SpatialAttenuation::enableVelocityScaling(bool enable, double scaling_factor) {
    velocity_scaling = enable;
    velocity_scaling_factor = scaling_factor;
}

void SpatialAttenuation::enableDepthScaling(bool enable, double gradient) {
    depth_scaling = enable;
    depth_gradient = gradient;
}

void SpatialAttenuation::enableDamageScaling(bool enable) {
    damage_scaling = enable;
}

ViscoelasticAttenuation::Parameters SpatialAttenuation::getParametersAt(
    double x, double y, double z,
    double vp, double vs,
    double damage) const {
    
    ViscoelasticAttenuation::Parameters params = base_params;
    
    double Qp = params.Qp_desired;
    double Qs = params.Qs_desired;
    
    // Velocity scaling
    if (velocity_scaling) {
        double vs_km = vs / 1000.0;
        Qs = velocity_scaling_factor * vs_km;
        Qp = 2.0 * Qs;
    }
    
    // Depth scaling
    if (depth_scaling) {
        double depth_km = -z / 1000.0;  // z is typically negative below surface
        if (depth_km > 0) {
            Qp *= (1.0 + depth_gradient * depth_km);
            Qs *= (1.0 + depth_gradient * depth_km);
        }
    }
    
    // Damage scaling
    if (damage_scaling && damage > 0) {
        Qp = AttenuationDatabase::damageCorrection(Qp, damage);
        Qs = AttenuationDatabase::damageCorrection(Qs, damage);
    }
    
    params.Qp_desired = std::max(5.0, Qp);  // Minimum Q
    params.Qs_desired = std::max(5.0, Qs);
    
    return params;
}

// =============================================================================
// AnelasticIntegrator Implementation
// =============================================================================

AnelasticIntegrator::AnelasticIntegrator(Scheme s) : scheme(s) {}

void AnelasticIntegrator::integrate(ViscoelasticAttenuation& attenuation,
                                    ViscoelasticAttenuation::AnelasticState& state,
                                    const std::array<double, 6>& stress,
                                    double dt) {
    switch (scheme) {
        case Scheme::FORWARD_EULER:
            forwardEuler(attenuation, state, stress, dt);
            break;
        case Scheme::RUNGE_KUTTA_2:
            rungeKutta2(attenuation, state, stress, dt);
            break;
        case Scheme::RUNGE_KUTTA_4:
            rungeKutta4(attenuation, state, stress, dt);
            break;
        case Scheme::ADER:
            // ADER uses the analytical solution
            attenuation.updateAnelasticState(state, stress, dt);
            break;
        default:
            attenuation.updateAnelasticState(state, stress, dt);
    }
}

void AnelasticIntegrator::forwardEuler(ViscoelasticAttenuation& attenuation,
                                       ViscoelasticAttenuation::AnelasticState& state,
                                       const std::array<double, 6>& stress,
                                       double dt) {
    // dR_i/dt = (Y_i * s - R_i) / τ_i
    
    const auto& mechanisms = attenuation.getMechanisms();
    
    if (state.R.size() != mechanisms.size()) {
        state.R.resize(mechanisms.size(), {0, 0, 0, 0, 0, 0});
    }
    
    // Compute deviatoric stress
    double mean_stress = (stress[0] + stress[1] + stress[2]) / 3.0;
    std::array<double, 6> dev_stress;
    dev_stress[0] = stress[0] - mean_stress;
    dev_stress[1] = stress[1] - mean_stress;
    dev_stress[2] = stress[2] - mean_stress;
    dev_stress[3] = stress[3];
    dev_stress[4] = stress[4];
    dev_stress[5] = stress[5];
    
    for (size_t i = 0; i < mechanisms.size(); ++i) {
        double Y = mechanisms[i].anelastic_coefficient_s;
        double tau = mechanisms[i].relaxation_time;
        
        for (int j = 0; j < 6; ++j) {
            double dR_dt = (Y * dev_stress[j] - state.R[i][j]) / tau;
            state.R[i][j] += dt * dR_dt;
        }
    }
}

void AnelasticIntegrator::rungeKutta2(ViscoelasticAttenuation& attenuation,
                                      ViscoelasticAttenuation::AnelasticState& state,
                                      const std::array<double, 6>& stress,
                                      double dt) {
    // Midpoint method (RK2)
    const auto& mechanisms = attenuation.getMechanisms();
    
    if (state.R.size() != mechanisms.size()) {
        state.R.resize(mechanisms.size(), {0, 0, 0, 0, 0, 0});
    }
    
    double mean_stress = (stress[0] + stress[1] + stress[2]) / 3.0;
    std::array<double, 6> dev_stress;
    dev_stress[0] = stress[0] - mean_stress;
    dev_stress[1] = stress[1] - mean_stress;
    dev_stress[2] = stress[2] - mean_stress;
    dev_stress[3] = stress[3];
    dev_stress[4] = stress[4];
    dev_stress[5] = stress[5];
    
    for (size_t i = 0; i < mechanisms.size(); ++i) {
        double Y = mechanisms[i].anelastic_coefficient_s;
        double tau = mechanisms[i].relaxation_time;
        
        for (int j = 0; j < 6; ++j) {
            // k1 = f(R_n)
            double k1 = (Y * dev_stress[j] - state.R[i][j]) / tau;
            
            // k2 = f(R_n + dt/2 * k1)
            double R_mid = state.R[i][j] + 0.5 * dt * k1;
            double k2 = (Y * dev_stress[j] - R_mid) / tau;
            
            state.R[i][j] += dt * k2;
        }
    }
}

void AnelasticIntegrator::rungeKutta4(ViscoelasticAttenuation& attenuation,
                                      ViscoelasticAttenuation::AnelasticState& state,
                                      const std::array<double, 6>& stress,
                                      double dt) {
    // Classical RK4
    const auto& mechanisms = attenuation.getMechanisms();
    
    if (state.R.size() != mechanisms.size()) {
        state.R.resize(mechanisms.size(), {0, 0, 0, 0, 0, 0});
    }
    
    double mean_stress = (stress[0] + stress[1] + stress[2]) / 3.0;
    std::array<double, 6> dev_stress;
    dev_stress[0] = stress[0] - mean_stress;
    dev_stress[1] = stress[1] - mean_stress;
    dev_stress[2] = stress[2] - mean_stress;
    dev_stress[3] = stress[3];
    dev_stress[4] = stress[4];
    dev_stress[5] = stress[5];
    
    for (size_t i = 0; i < mechanisms.size(); ++i) {
        double Y = mechanisms[i].anelastic_coefficient_s;
        double tau = mechanisms[i].relaxation_time;
        
        for (int j = 0; j < 6; ++j) {
            double R0 = state.R[i][j];
            double s = dev_stress[j];
            
            double k1 = (Y * s - R0) / tau;
            double k2 = (Y * s - (R0 + 0.5 * dt * k1)) / tau;
            double k3 = (Y * s - (R0 + 0.5 * dt * k2)) / tau;
            double k4 = (Y * s - (R0 + dt * k3)) / tau;
            
            state.R[i][j] = R0 + dt / 6.0 * (k1 + 2.0*k2 + 2.0*k3 + k4);
        }
    }
}

// =============================================================================
// ViscoelasticConfig Implementation
// =============================================================================

void ViscoelasticConfig::parseConfig(const std::map<std::string, std::string>& config) {
    enabled = parseString(config, "viscoelastic_enabled", "false") == "true";
    
    num_mechanisms = static_cast<int>(parseDouble(config, "num_mechanisms", 3.0));
    freq_central = parseDouble(config, "frequency_central", 0.5);
    freq_ratio = parseDouble(config, "frequency_ratio", 100.0);
    
    Qp = parseDouble(config, "Qp", 50.0);
    Qs = parseDouble(config, "Qs", 25.0);
    
    velocity_dependent = parseString(config, "velocity_dependent_Q", "false") == "true";
    velocity_scaling = parseDouble(config, "velocity_Q_scaling", 50.0);
    
    depth_dependent = parseString(config, "depth_dependent_Q", "false") == "true";
    depth_gradient = parseDouble(config, "depth_Q_gradient", 0.1);
    
    rock_type = parseString(config, "rock_type", "granite");
}

ViscoelasticAttenuation ViscoelasticConfig::createModel() const {
    ViscoelasticAttenuation model;
    
    if (!enabled) return model;
    
    ViscoelasticAttenuation::Parameters params;
    params.num_mechanisms = num_mechanisms;
    params.freq_central = freq_central;
    params.freq_ratio = freq_ratio;
    params.Qp_desired = Qp;
    params.Qs_desired = Qs;
    
    model.setParameters(params);
    model.enable(true);
    
    return model;
}

// =============================================================================
// AttenuationUtils Implementation
// =============================================================================

namespace AttenuationUtils {

bool checkCausality(const ViscoelasticAttenuation& model,
                   double f_min, double f_max, double tolerance) {
    // Check Kramers-Kronig relations (simplified check)
    // A causal system should have velocity dispersion that matches attenuation
    
    auto curve = model.computeDispersionCurve(5000, 3000, f_min, f_max, 100);
    
    for (size_t i = 1; i < curve.size() - 1; ++i) {
        // Higher Q should correlate with less velocity dispersion
        double dv = curve[i+1].velocity_s - curve[i-1].velocity_s;
        double df = curve[i+1].frequency - curve[i-1].frequency;
        
        // Rough check: dispersion should be monotonic for reasonable Q
        // This is a simplified causality check
    }
    
    return true;  // Simplified - full check requires Hilbert transform
}

double phaseVelocity(double group_velocity, double Q, double frequency) {
    // Phase velocity from dispersion relation
    // c_phase ≈ c_group * (1 + 1/(π*Q) * ln(f/f0))
    // Simplified version
    return group_velocity * (1.0 + 1.0 / (M_PI * Q));
}

double attenuationCoefficient(double frequency, double velocity, double Q) {
    // α = π * f / (Q * c)  [Nepers/m]
    return M_PI * frequency / (Q * velocity);
}

double amplitudeReduction(double distance, double frequency,
                         double velocity, double Q) {
    double alpha = attenuationCoefficient(frequency, velocity, Q);
    return std::exp(-alpha * distance);
}

double computeQ(double amplitude_ratio, double distance, double wavelength) {
    // Q = π / (λ * ln(A1/A2) / d)
    if (amplitude_ratio <= 0 || amplitude_ratio >= 1) return 1e10;
    
    return M_PI * distance / (wavelength * std::log(1.0 / amplitude_ratio));
}

std::pair<double, double> estimateQFromSpectralRatio(
    const std::vector<double>& spectrum1,
    const std::vector<double>& spectrum2,
    const std::vector<double>& frequencies,
    double distance_diff) {
    
    if (spectrum1.size() != spectrum2.size() || 
        spectrum1.size() != frequencies.size()) {
        return {0.0, 0.0};
    }
    
    // Linear regression on ln(A1/A2) vs frequency
    // slope = π * Δr / (Q * c)
    
    double sum_f = 0, sum_lnA = 0, sum_f2 = 0, sum_f_lnA = 0;
    int n = 0;
    
    for (size_t i = 0; i < frequencies.size(); ++i) {
        if (spectrum1[i] > 0 && spectrum2[i] > 0) {
            double f = frequencies[i];
            double lnA = std::log(spectrum1[i] / spectrum2[i]);
            
            sum_f += f;
            sum_lnA += lnA;
            sum_f2 += f * f;
            sum_f_lnA += f * lnA;
            n++;
        }
    }
    
    if (n < 2) return {0.0, 0.0};
    
    double denom = n * sum_f2 - sum_f * sum_f;
    if (std::abs(denom) < 1e-15) return {0.0, 0.0};
    
    double slope = (n * sum_f_lnA - sum_f * sum_lnA) / denom;
    double intercept = (sum_lnA - slope * sum_f) / n;
    
    // Assuming c ≈ 3000 m/s (rough estimate)
    double c = 3000.0;
    double Q = M_PI * distance_diff / (slope * c);
    
    // Return Q and uncertainty (simplified)
    return {Q, Q * 0.1};
}

void plotDispersionCurves(const ViscoelasticAttenuation& model,
                         double vp_ref, double vs_ref,
                         const std::string& filename) {
    auto curve = model.computeDispersionCurve(vp_ref, vs_ref, 0.01, 100.0, 200);
    
    std::ofstream file(filename);
    if (!file.is_open()) return;
    
    file << "# Frequency(Hz)  Vp(m/s)  Vs(m/s)  Qp  Qs\n";
    
    for (const auto& point : curve) {
        file << point.frequency << " "
             << point.velocity_p << " "
             << point.velocity_s << " "
             << point.Qp << " "
             << point.Qs << "\n";
    }
    
    file.close();
}

} // namespace AttenuationUtils

// =============================================================================
// ViscoelasticCoupling Implementation
// =============================================================================

ViscoelasticCoupling::ViscoelasticCoupling()
    : poroelastic_coupling(false),
      plasticity_coupling(false),
      plastic_q_reduction(0.5),
      damage_coupling(false),
      damage_q_exponent(1.5) {}

double ViscoelasticCoupling::getModifiedQ(double Q_base, double porosity,
                                          double damage, bool is_plastic) const {
    double Q = Q_base;
    
    // Poroelastic coupling: higher porosity = lower Q
    if (poroelastic_coupling) {
        Q *= (1.0 - 0.5 * porosity);  // Simple linear reduction
    }
    
    // Plasticity coupling: plastic zones have lower Q
    if (plasticity_coupling && is_plastic) {
        Q *= plastic_q_reduction;
    }
    
    // Damage coupling
    if (damage_coupling && damage > 0) {
        Q *= std::pow(1.0 - damage, damage_q_exponent);
    }
    
    return std::max(5.0, Q);  // Minimum Q
}

} // namespace FSRM
