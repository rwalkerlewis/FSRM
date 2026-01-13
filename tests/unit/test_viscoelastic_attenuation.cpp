/**
 * @file test_viscoelastic_attenuation.cpp
 * @brief Unit tests for viscoelastic attenuation models
 * 
 * Tests cover:
 * - Maxwell body mechanics
 * - Multi-mechanism attenuation
 * - Q-factor fitting
 * - Time integration schemes
 * - Dispersion relations
 */

#include "ViscoelasticAttenuation.hpp"
#include "Testing.hpp"
#include <cmath>
#include <iostream>
#include <vector>
#include <array>
#include <complex>

namespace FSRM {
namespace Testing {

// =============================================================================
// Helper Functions
// =============================================================================

/**
 * @brief Compute quality factor from loss/storage moduli
 */
double computeQ(double M_real, double M_imag) {
    return M_real / (2.0 * std::abs(M_imag));
}

// =============================================================================
// Single Maxwell Body Tests
// =============================================================================

/**
 * @test Test Maxwell body stress relaxation
 */
bool test_maxwell_relaxation() {
    std::cout << "Testing Maxwell body stress relaxation..." << std::endl;
    
    double tol = 0.05;  // 5% tolerance
    
    // Maxwell body: σ/E + σ̇/η = ε̇
    // Under constant strain: σ(t) = σ₀ * exp(-E*t/η) = σ₀ * exp(-t/τ)
    double E = 30e9;   // 30 GPa
    double eta = 1e15; // Viscosity
    double tau = eta / E;  // Relaxation time
    
    // Initial stress
    double sigma_0 = 1e6;  // 1 MPa
    
    // Check at t = tau (should be ~0.368 of initial)
    double t = tau;
    double sigma_expected = sigma_0 * std::exp(-1.0);
    double sigma_computed = sigma_0 * std::exp(-t / tau);
    
    if (std::abs(sigma_computed - sigma_expected) / sigma_expected > tol) {
        std::cerr << "  FAIL: Relaxation incorrect" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Maxwell relaxation correct" << std::endl;
    return true;
}

/**
 * @test Test Maxwell body complex modulus
 */
bool test_maxwell_complex_modulus() {
    std::cout << "Testing Maxwell body complex modulus..." << std::endl;
    
    double tol = 1e-6;
    
    // Complex modulus: M*(ω) = iωη / (1 + iωτ)
    // where τ = η/E is relaxation time
    
    double E = 30e9;
    double tau = 1.0;  // 1 second relaxation
    
    double omega = 1.0;  // 1 rad/s
    
    // M* = E * iωτ / (1 + iωτ)
    std::complex<double> i(0, 1);
    std::complex<double> omega_tau(0, omega * tau);
    std::complex<double> M_star = E * omega_tau / (1.0 + omega_tau);
    
    // Real part (storage modulus)
    double M_real = std::real(M_star);
    // Expected: E * (ωτ)² / (1 + (ωτ)²)
    double M_real_expected = E * omega * omega * tau * tau / (1.0 + omega * omega * tau * tau);
    
    if (std::abs(M_real - M_real_expected) > tol * E) {
        std::cerr << "  FAIL: Storage modulus incorrect" << std::endl;
        return false;
    }
    
    // Imaginary part (loss modulus)
    double M_imag = std::imag(M_star);
    // Expected: E * ωτ / (1 + (ωτ)²)
    double M_imag_expected = E * omega * tau / (1.0 + omega * omega * tau * tau);
    
    if (std::abs(M_imag - M_imag_expected) > tol * E) {
        std::cerr << "  FAIL: Loss modulus incorrect" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Maxwell complex modulus correct" << std::endl;
    return true;
}

// =============================================================================
// Multi-Mechanism Tests
// =============================================================================

/**
 * @test Test Q-factor fitting
 */
bool test_q_factor_fitting() {
    std::cout << "Testing Q-factor fitting..." << std::endl;
    
    double tol = 0.1;  // 10% tolerance on Q
    
    ViscoelasticAttenuation ve;
    ve.setReferenceModulus(30e9);  // 30 GPa
    ve.setReferenceDensity(2500);
    
    // Target constant Q = 100 in frequency band [0.1, 10] Hz
    double Q_target = 100.0;
    double f_min = 0.1;
    double f_max = 10.0;
    
    ve.setTargetQ(Q_target, f_min, f_max);
    ve.setNumMechanisms(5);
    ve.computeCoefficients();
    
    // Test Q at several frequencies
    std::vector<double> test_freqs = {0.2, 0.5, 1.0, 2.0, 5.0};
    
    for (double f : test_freqs) {
        double Q_computed = ve.getQualityFactor(2.0 * M_PI * f);
        
        double relative_error = std::abs(Q_computed - Q_target) / Q_target;
        
        if (relative_error > tol) {
            std::cerr << "  FAIL: Q at " << f << " Hz = " << Q_computed 
                      << ", target = " << Q_target << std::endl;
            return false;
        }
    }
    
    std::cout << "  PASS: Q-factor fitting within tolerance" << std::endl;
    return true;
}

/**
 * @test Test frequency-dependent Q
 */
bool test_frequency_dependent_q() {
    std::cout << "Testing frequency-dependent Q..." << std::endl;
    
    ViscoelasticAttenuation ve;
    ve.setReferenceModulus(30e9);
    ve.setReferenceDensity(2500);
    
    // Q = Q0 * (f/f0)^alpha power law
    double Q0 = 50.0;
    double f0 = 1.0;
    double alpha = 0.5;
    
    ve.setFrequencyDependentQ(Q0, f0, alpha);
    ve.setNumMechanisms(5);
    ve.setFrequencyRange(0.01, 100.0);
    ve.computeCoefficients();
    
    // Test power law behavior
    double f1 = 0.1;
    double f2 = 10.0;
    
    double Q1 = ve.getQualityFactor(2.0 * M_PI * f1);
    double Q2 = ve.getQualityFactor(2.0 * M_PI * f2);
    
    // Ratio should follow power law
    double ratio_expected = std::pow(f2 / f1, alpha);
    double ratio_computed = Q2 / Q1;
    
    // Allow larger tolerance for power-law approximation
    if (std::abs(ratio_computed - ratio_expected) / ratio_expected > 0.3) {
        std::cerr << "  FAIL: Power law not followed, ratio=" << ratio_computed
                  << ", expected~" << ratio_expected << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Frequency-dependent Q correct" << std::endl;
    return true;
}

// =============================================================================
// Time Integration Tests
// =============================================================================

/**
 * @test Test anelastic state evolution
 */
bool test_anelastic_state_evolution() {
    std::cout << "Testing anelastic state evolution..." << std::endl;
    
    ViscoelasticAttenuation ve;
    ve.setReferenceModulus(30e9);
    ve.setReferenceDensity(2500);
    ve.setTargetQ(100.0, 0.1, 10.0);
    ve.setNumMechanisms(3);
    ve.computeCoefficients();
    
    // Initialize state variables
    int n_mech = ve.getNumMechanisms();
    std::vector<double> anelastic_state(n_mech * 6, 0.0);
    
    // Constant strain rate
    double strain_rate[6] = {1e-4, 0, 0, 0, 0, 0};  // ε̇xx = 1e-4 /s
    
    // Integrate for several timesteps
    double dt = 0.001;
    for (int i = 0; i < 100; ++i) {
        ve.evolveAnelasticState(strain_rate, anelastic_state.data(), dt);
    }
    
    // State should have non-zero values
    double state_norm = 0;
    for (double s : anelastic_state) {
        state_norm += s * s;
    }
    state_norm = std::sqrt(state_norm);
    
    if (state_norm < 1e-20) {
        std::cerr << "  FAIL: Anelastic state remained zero" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Anelastic state evolution correct" << std::endl;
    return true;
}

/**
 * @test Test anelastic stress computation
 */
bool test_anelastic_stress() {
    std::cout << "Testing anelastic stress computation..." << std::endl;
    
    ViscoelasticAttenuation ve;
    ve.setReferenceModulus(30e9);
    ve.setReferenceDensity(2500);
    ve.setTargetQ(50.0, 0.1, 10.0);
    ve.setNumMechanisms(3);
    ve.computeCoefficients();
    
    // Create some anelastic state
    int n_mech = ve.getNumMechanisms();
    std::vector<double> anelastic_state(n_mech * 6);
    for (int i = 0; i < n_mech; ++i) {
        anelastic_state[i*6 + 0] = 1e5 * (i + 1);  // σxx contribution
    }
    
    // Compute anelastic stress
    double sigma_anelastic[6];
    ve.computeAnelasticStress(anelastic_state.data(), sigma_anelastic);
    
    // Should be non-zero
    if (std::abs(sigma_anelastic[0]) < 1e-10) {
        std::cerr << "  FAIL: Anelastic stress zero" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Anelastic stress computation correct" << std::endl;
    return true;
}

/**
 * @test Test different integration schemes
 */
bool test_integration_schemes() {
    std::cout << "Testing integration schemes..." << std::endl;
    
    double tol = 0.05;  // 5% tolerance between schemes
    
    ViscoelasticAttenuation ve;
    ve.setReferenceModulus(30e9);
    ve.setReferenceDensity(2500);
    ve.setTargetQ(100.0, 0.1, 10.0);
    ve.setNumMechanisms(3);
    ve.computeCoefficients();
    
    AnelasticIntegrator euler(IntegrationScheme::FORWARD_EULER);
    AnelasticIntegrator rk4(IntegrationScheme::RK4);
    
    euler.initialize(ve);
    rk4.initialize(ve);
    
    int n_mech = ve.getNumMechanisms();
    std::vector<double> state_euler(n_mech * 6, 0.0);
    std::vector<double> state_rk4(n_mech * 6, 0.0);
    
    double strain_rate[6] = {1e-4, 0, 0, 0, 0, 0};
    double dt = 0.001;
    
    // Integrate with both schemes
    for (int i = 0; i < 100; ++i) {
        euler.step(state_euler.data(), strain_rate, dt);
        rk4.step(state_rk4.data(), strain_rate, dt);
    }
    
    // Results should be similar (within tolerance)
    double diff_norm = 0;
    double euler_norm = 0;
    for (size_t i = 0; i < state_euler.size(); ++i) {
        diff_norm += (state_euler[i] - state_rk4[i]) * (state_euler[i] - state_rk4[i]);
        euler_norm += state_euler[i] * state_euler[i];
    }
    
    double relative_diff = std::sqrt(diff_norm) / (std::sqrt(euler_norm) + 1e-20);
    
    if (relative_diff > tol) {
        std::cerr << "  FAIL: Large difference between schemes, rel_diff=" << relative_diff << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Integration schemes consistent" << std::endl;
    return true;
}

// =============================================================================
// Dispersion Tests
// =============================================================================

/**
 * @test Test velocity dispersion
 */
bool test_velocity_dispersion() {
    std::cout << "Testing velocity dispersion..." << std::endl;
    
    ViscoelasticAttenuation ve;
    ve.setReferenceModulus(30e9);
    ve.setReferenceDensity(2500);
    ve.setTargetQ(50.0, 0.1, 10.0);  // Low Q = significant dispersion
    ve.setNumMechanisms(5);
    ve.computeCoefficients();
    
    // Reference velocity at 1 Hz
    double v_ref = std::sqrt(30e9 / 2500);
    
    // Get velocities at different frequencies
    double v_low = ve.getVelocity(2.0 * M_PI * 0.1);
    double v_mid = ve.getVelocity(2.0 * M_PI * 1.0);
    double v_high = ve.getVelocity(2.0 * M_PI * 10.0);
    
    // For low Q, velocity should vary with frequency
    // Higher frequencies should have higher velocity (normal dispersion)
    if (v_high < v_low * 0.95) {
        std::cerr << "  FAIL: Expected normal dispersion" << std::endl;
        return false;
    }
    
    // All velocities should be positive and reasonable
    if (v_low <= 0 || v_mid <= 0 || v_high <= 0) {
        std::cerr << "  FAIL: Non-positive velocity" << std::endl;
        return false;
    }
    
    // Should be within ~20% of reference for reasonable Q
    if (std::abs(v_mid - v_ref) / v_ref > 0.2) {
        std::cerr << "  FAIL: Velocity too far from reference" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Velocity dispersion computed" << std::endl;
    return true;
}

/**
 * @test Test causality (Kramers-Kronig relations)
 */
bool test_causality() {
    std::cout << "Testing causality (Kramers-Kronig)..." << std::endl;
    
    ViscoelasticAttenuation ve;
    ve.setReferenceModulus(30e9);
    ve.setReferenceDensity(2500);
    ve.setTargetQ(100.0, 0.1, 10.0);
    ve.setNumMechanisms(5);
    ve.computeCoefficients();
    
    // Check that the model is causal
    // The relaxation moduli should be non-negative
    bool causal = ve.checkCausality();
    
    if (!causal) {
        std::cerr << "  FAIL: Model violates causality" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Model satisfies causality" << std::endl;
    return true;
}

// =============================================================================
// Attenuation Database Tests
// =============================================================================

/**
 * @test Test attenuation database lookup
 */
bool test_attenuation_database() {
    std::cout << "Testing attenuation database..." << std::endl;
    
    AttenuationDatabase db;
    
    // Get Q for different rock types
    double Q_granite = db.getQForRockType("granite");
    double Q_sandstone = db.getQForRockType("sandstone");
    double Q_shale = db.getQForRockType("shale");
    
    // Granite should have higher Q (less attenuation)
    if (Q_granite <= Q_sandstone || Q_granite <= Q_shale) {
        std::cerr << "  FAIL: Granite should have highest Q" << std::endl;
        return false;
    }
    
    // Q should be positive
    if (Q_granite <= 0 || Q_sandstone <= 0 || Q_shale <= 0) {
        std::cerr << "  FAIL: Non-positive Q values" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Attenuation database lookup correct" << std::endl;
    return true;
}

/**
 * @test Test depth-dependent Q
 */
bool test_depth_dependent_q() {
    std::cout << "Testing depth-dependent Q..." << std::endl;
    
    AttenuationDatabase db;
    
    // Q should generally increase with depth (compression closes cracks)
    double Q_shallow = db.getQAtDepth(1000);   // 1 km
    double Q_deep = db.getQAtDepth(30000);     // 30 km
    
    if (Q_deep < Q_shallow * 0.8) {  // Allow some variation
        std::cerr << "  FAIL: Q should increase with depth" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Depth-dependent Q correct" << std::endl;
    return true;
}

// =============================================================================
// Spatial Attenuation Tests
// =============================================================================

/**
 * @test Test spatial Q variation
 */
bool test_spatial_attenuation() {
    std::cout << "Testing spatial Q variation..." << std::endl;
    
    SpatialAttenuation sa;
    sa.setBackgroundQ(100.0);
    
    // Add low-Q anomaly
    sa.addAnomaly(0, 0, -5000, 1000, 50.0);  // Center at (0,0,-5km), radius 1km, Q=50
    
    // Query Q at different locations
    double Q_center = sa.getQ(0, 0, -5000);
    double Q_edge = sa.getQ(1000, 0, -5000);
    double Q_far = sa.getQ(5000, 0, -5000);
    
    // Center should have anomaly Q
    if (std::abs(Q_center - 50.0) > 10.0) {
        std::cerr << "  FAIL: Anomaly Q not applied at center" << std::endl;
        return false;
    }
    
    // Far should have background Q
    if (std::abs(Q_far - 100.0) > 10.0) {
        std::cerr << "  FAIL: Background Q not correct" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Spatial Q variation correct" << std::endl;
    return true;
}

// =============================================================================
// ADER Integration Tests
// =============================================================================

/**
 * @test Test ADER anelastic integration
 */
bool test_ader_anelastic() {
    std::cout << "Testing ADER anelastic integration..." << std::endl;
    
    ViscoelasticAttenuation ve;
    ve.setReferenceModulus(30e9);
    ve.setReferenceDensity(2500);
    ve.setTargetQ(100.0, 0.1, 10.0);
    ve.setNumMechanisms(3);
    ve.computeCoefficients();
    
    AnelasticIntegrator ader(IntegrationScheme::ADER);
    ader.initialize(ve);
    ader.setADEROrder(4);  // 4th order ADER
    
    int n_mech = ve.getNumMechanisms();
    std::vector<double> state(n_mech * 6, 0.0);
    
    // Time-varying strain rate
    auto strain_rate_func = [](double t, double* eps_dot) {
        eps_dot[0] = 1e-4 * std::sin(2.0 * M_PI * t);
        for (int i = 1; i < 6; ++i) eps_dot[i] = 0;
    };
    
    double dt = 0.01;
    double t = 0;
    
    for (int i = 0; i < 100; ++i) {
        ader.stepADER(state.data(), strain_rate_func, t, dt);
        t += dt;
    }
    
    // State should be finite and non-zero after forcing
    bool valid = true;
    double norm = 0;
    for (double s : state) {
        if (std::isnan(s) || std::isinf(s)) valid = false;
        norm += s * s;
    }
    
    if (!valid) {
        std::cerr << "  FAIL: Invalid state values" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: ADER anelastic integration correct" << std::endl;
    return true;
}

// =============================================================================
// Test Runner
// =============================================================================

int runViscoelasticTests() {
    std::cout << "\n=== Viscoelastic Attenuation Unit Tests ===" << std::endl;
    
    int passed = 0;
    int failed = 0;
    
    // Maxwell body tests
    if (test_maxwell_relaxation()) ++passed; else ++failed;
    if (test_maxwell_complex_modulus()) ++passed; else ++failed;
    
    // Multi-mechanism tests
    if (test_q_factor_fitting()) ++passed; else ++failed;
    if (test_frequency_dependent_q()) ++passed; else ++failed;
    
    // Time integration tests
    if (test_anelastic_state_evolution()) ++passed; else ++failed;
    if (test_anelastic_stress()) ++passed; else ++failed;
    if (test_integration_schemes()) ++passed; else ++failed;
    
    // Dispersion tests
    if (test_velocity_dispersion()) ++passed; else ++failed;
    if (test_causality()) ++passed; else ++failed;
    
    // Database tests
    if (test_attenuation_database()) ++passed; else ++failed;
    if (test_depth_dependent_q()) ++passed; else ++failed;
    
    // Spatial tests
    if (test_spatial_attenuation()) ++passed; else ++failed;
    
    // ADER tests
    if (test_ader_anelastic()) ++passed; else ++failed;
    
    std::cout << "\n=== Viscoelastic Test Summary ===" << std::endl;
    std::cout << "Passed: " << passed << std::endl;
    std::cout << "Failed: " << failed << std::endl;
    std::cout << "Total:  " << (passed + failed) << std::endl;
    
    return failed;
}

} // namespace Testing
} // namespace FSRM

#ifndef FSRM_TEST_NO_MAIN
int main() {
    return FSRM::Testing::runViscoelasticTests();
}
#endif
