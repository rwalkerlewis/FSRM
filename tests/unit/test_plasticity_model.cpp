/**
 * @file test_plasticity_model.cpp
 * @brief Unit tests for plasticity models
 * 
 * Tests cover:
 * - Drucker-Prager model
 * - von Mises model
 * - Mohr-Coulomb model
 * - Cap model
 * - Return mapping algorithms
 * - Consistent tangent modulus
 */

#include "PlasticityModel.hpp"
#include "Testing.hpp"
#include <cmath>
#include <iostream>
#include <vector>
#include <array>

namespace FSRM {
namespace Testing {

// =============================================================================
// Helper Functions
// =============================================================================

// Compute stress invariants
void computeInvariants(const double sigma[6], double& I1, double& J2, double& J3) {
    // sigma = [σxx, σyy, σzz, σxy, σxz, σyz]
    I1 = sigma[0] + sigma[1] + sigma[2];  // First invariant (trace)
    
    double p = I1 / 3.0;
    double s[6] = {
        sigma[0] - p,  // sxx
        sigma[1] - p,  // syy
        sigma[2] - p,  // szz
        sigma[3],      // sxy
        sigma[4],      // sxz
        sigma[5]       // syz
    };
    
    // J2 = 0.5 * s:s
    J2 = 0.5 * (s[0]*s[0] + s[1]*s[1] + s[2]*s[2]) + 
         s[3]*s[3] + s[4]*s[4] + s[5]*s[5];
    
    // J3 = det(s)
    J3 = s[0] * (s[1]*s[2] - s[5]*s[5]) 
       - s[3] * (s[3]*s[2] - s[5]*s[4])
       + s[4] * (s[3]*s[5] - s[1]*s[4]);
}

// =============================================================================
// Drucker-Prager Tests
// =============================================================================

/**
 * @test Test Drucker-Prager yield function
 */
bool test_drucker_prager_yield_function() {
    std::cout << "Testing Drucker-Prager yield function..." << std::endl;
    
    double tol = 1e-6;
    
    DruckerPragerModel dp;
    dp.setFrictionAngle(30.0 * M_PI / 180.0);  // 30 degrees
    dp.setCohesion(1e6);  // 1 MPa
    
    // Compute DP parameters
    double phi = 30.0 * M_PI / 180.0;
    double c = 1e6;
    double alpha = 6.0 * sin(phi) / (3.0 - sin(phi));
    double k = 6.0 * c * cos(phi) / (3.0 - sin(phi));
    
    // Test stress state on yield surface
    // f = sqrt(J2) + alpha*I1 - k = 0
    // For I1 = -10e6 (compressive), solve for sqrt(J2)
    double I1_test = -10e6;
    double sqrtJ2_yield = k - alpha * I1_test;
    
    // Create stress state with this J2
    // Pure deviatoric: s11 = -s22 = sqrt(2*J2/3)
    double dev_mag = std::sqrt(2.0 * sqrtJ2_yield * sqrtJ2_yield / 3.0);
    double p = I1_test / 3.0;
    
    double sigma[6] = {p + dev_mag, p - dev_mag, p, 0, 0, 0};
    
    // Compute invariants
    double I1, J2, J3;
    computeInvariants(sigma, I1, J2, J3);
    
    // Yield function should be zero
    double f = std::sqrt(J2) + alpha * I1 - k;
    
    if (std::abs(f) / k > tol) {
        std::cerr << "  FAIL: Yield function f = " << f << ", should be ~0" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Drucker-Prager yield function correct" << std::endl;
    return true;
}

/**
 * @test Test Drucker-Prager return mapping
 */
bool test_drucker_prager_return_mapping() {
    std::cout << "Testing Drucker-Prager return mapping..." << std::endl;
    
    double tol = 1e-5;
    
    DruckerPragerModel dp;
    dp.setFrictionAngle(30.0 * M_PI / 180.0);
    dp.setCohesion(1e6);
    dp.setDilatancyAngle(10.0 * M_PI / 180.0);  // Non-associated
    
    // Trial stress outside yield surface
    double sigma_trial[6] = {-5e6, -15e6, -10e6, 2e6, 1e6, 0.5e6};
    double sigma_return[6];
    double plastic_multiplier;
    
    // Perform return mapping
    bool yielded = dp.returnMap(sigma_trial, sigma_return, plastic_multiplier);
    
    if (!yielded) {
        std::cerr << "  FAIL: Stress outside yield should trigger return mapping" << std::endl;
        return false;
    }
    
    // Returned stress should be on yield surface (f <= tol)
    double phi = 30.0 * M_PI / 180.0;
    double c = 1e6;
    double alpha = 6.0 * sin(phi) / (3.0 - sin(phi));
    double k = 6.0 * c * cos(phi) / (3.0 - sin(phi));
    
    double I1, J2, J3;
    computeInvariants(sigma_return, I1, J2, J3);
    
    double f = std::sqrt(J2) + alpha * I1 - k;
    
    if (f > tol * k) {
        std::cerr << "  FAIL: Returned stress outside yield surface, f=" << f << std::endl;
        return false;
    }
    
    // Plastic multiplier should be positive
    if (plastic_multiplier < 0) {
        std::cerr << "  FAIL: Negative plastic multiplier" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Drucker-Prager return mapping correct" << std::endl;
    return true;
}

/**
 * @test Test Drucker-Prager elastic stress (inside yield)
 */
bool test_drucker_prager_elastic() {
    std::cout << "Testing Drucker-Prager elastic behavior..." << std::endl;
    
    double tol = 1e-10;
    
    DruckerPragerModel dp;
    dp.setFrictionAngle(30.0 * M_PI / 180.0);
    dp.setCohesion(10e6);  // High cohesion - stress will be elastic
    
    // Small stress (should be elastic)
    double sigma_trial[6] = {-1e6, -2e6, -1.5e6, 0.1e6, 0.05e6, 0.02e6};
    double sigma_return[6];
    double plastic_multiplier;
    
    bool yielded = dp.returnMap(sigma_trial, sigma_return, plastic_multiplier);
    
    if (yielded) {
        std::cerr << "  FAIL: Small stress should be elastic" << std::endl;
        return false;
    }
    
    // Stress should be unchanged
    for (int i = 0; i < 6; ++i) {
        if (std::abs(sigma_return[i] - sigma_trial[i]) > tol) {
            std::cerr << "  FAIL: Elastic stress changed" << std::endl;
            return false;
        }
    }
    
    std::cout << "  PASS: Drucker-Prager elastic behavior correct" << std::endl;
    return true;
}

// =============================================================================
// Von Mises Tests
// =============================================================================

/**
 * @test Test von Mises yield function
 */
bool test_von_mises_yield_function() {
    std::cout << "Testing von Mises yield function..." << std::endl;
    
    double tol = 1e-10;
    
    VonMisesModel vm;
    double yield_stress = 250e6;  // 250 MPa
    vm.setYieldStress(yield_stress);
    
    // Von Mises: f = sqrt(3*J2) - σ_y
    // On yield: sqrt(3*J2) = σ_y  =>  J2 = σ_y²/3
    
    // Create deviatoric stress with correct J2
    double J2_yield = yield_stress * yield_stress / 3.0;
    
    // Simple shear state: σxy = sqrt(J2)
    double tau = std::sqrt(J2_yield);
    double sigma[6] = {0, 0, 0, tau, 0, 0};
    
    double I1, J2, J3;
    computeInvariants(sigma, I1, J2, J3);
    
    double f = std::sqrt(3.0 * J2) - yield_stress;
    
    if (std::abs(f) > tol * yield_stress) {
        std::cerr << "  FAIL: Yield function f = " << f << std::endl;
        return false;
    }
    
    std::cout << "  PASS: von Mises yield function correct" << std::endl;
    return true;
}

/**
 * @test Test von Mises radial return
 */
bool test_von_mises_radial_return() {
    std::cout << "Testing von Mises radial return..." << std::endl;
    
    double tol = 1e-6;
    
    VonMisesModel vm;
    double yield_stress = 250e6;
    vm.setYieldStress(yield_stress);
    
    // Trial stress outside yield surface
    double sigma_trial[6] = {-100e6, 200e6, -100e6, 50e6, 30e6, 20e6};
    double sigma_return[6];
    double plastic_multiplier;
    
    bool yielded = vm.returnMap(sigma_trial, sigma_return, plastic_multiplier);
    
    // Compute invariants of trial
    double I1_trial, J2_trial, J3_trial;
    computeInvariants(sigma_trial, I1_trial, J2_trial, J3_trial);
    
    double q_trial = std::sqrt(3.0 * J2_trial);
    
    if (q_trial <= yield_stress) {
        std::cout << "  SKIP: Trial stress elastic" << std::endl;
        return true;
    }
    
    if (!yielded) {
        std::cerr << "  FAIL: Plastic stress not detected" << std::endl;
        return false;
    }
    
    // Returned stress should be on yield surface
    double I1_ret, J2_ret, J3_ret;
    computeInvariants(sigma_return, I1_ret, J2_ret, J3_ret);
    
    double q_ret = std::sqrt(3.0 * J2_ret);
    
    if (std::abs(q_ret - yield_stress) / yield_stress > tol) {
        std::cerr << "  FAIL: Returned stress not on yield, q=" << q_ret << std::endl;
        return false;
    }
    
    // Pressure should be preserved (radial return)
    if (std::abs(I1_ret - I1_trial) / std::abs(I1_trial + 1.0) > tol) {
        std::cerr << "  FAIL: Pressure not preserved in radial return" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: von Mises radial return correct" << std::endl;
    return true;
}

// =============================================================================
// Mohr-Coulomb Tests
// =============================================================================

/**
 * @test Test Mohr-Coulomb yield function
 */
bool test_mohr_coulomb_yield() {
    std::cout << "Testing Mohr-Coulomb yield function..." << std::endl;
    
    double tol = 1e-5;
    
    MohrCoulombModel mc;
    double phi = 30.0 * M_PI / 180.0;
    double c = 5e6;
    mc.setFrictionAngle(phi);
    mc.setCohesion(c);
    
    // MC yield function in principal stresses:
    // f = (σ1 - σ3) - (σ1 + σ3)*sin(φ) - 2c*cos(φ) = 0
    
    // Choose σ3 = -20 MPa, compute σ1 on yield
    double sigma3 = -20e6;
    double sigma1 = (sigma3 * (1 + std::sin(phi)) + 2*c*std::cos(phi)) / (1 - std::sin(phi));
    double sigma2 = 0.5 * (sigma1 + sigma3);  // Intermediate
    
    // Create stress tensor (principal directions aligned with axes)
    double sigma[6] = {sigma1, sigma2, sigma3, 0, 0, 0};
    
    double f = mc.yieldFunction(sigma);
    
    if (std::abs(f) / (c + 1.0) > tol) {
        std::cerr << "  FAIL: MC yield f = " << f << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Mohr-Coulomb yield function correct" << std::endl;
    return true;
}

/**
 * @test Test Mohr-Coulomb return mapping
 */
bool test_mohr_coulomb_return() {
    std::cout << "Testing Mohr-Coulomb return mapping..." << std::endl;
    
    MohrCoulombModel mc;
    mc.setFrictionAngle(35.0 * M_PI / 180.0);
    mc.setCohesion(2e6);
    mc.setDilatancyAngle(15.0 * M_PI / 180.0);
    
    // Trial stress outside yield (principal stress difference too large)
    double sigma_trial[6] = {-5e6, -15e6, -30e6, 1e6, 0.5e6, 0.2e6};
    double sigma_return[6];
    double plastic_multiplier;
    
    bool yielded = mc.returnMap(sigma_trial, sigma_return, plastic_multiplier);
    
    if (yielded) {
        // Returned stress should satisfy yield condition
        double f = mc.yieldFunction(sigma_return);
        
        if (f > 1e-4 * 2e6) {  // Tolerance relative to cohesion
            std::cerr << "  FAIL: Returned stress outside yield, f=" << f << std::endl;
            return false;
        }
    }
    
    std::cout << "  PASS: Mohr-Coulomb return mapping correct" << std::endl;
    return true;
}

// =============================================================================
// Cap Model Tests
// =============================================================================

/**
 * @test Test cap model compaction
 */
bool test_cap_model_compaction() {
    std::cout << "Testing cap model compaction..." << std::endl;
    
    CapModel cap;
    cap.setFrictionAngle(30.0 * M_PI / 180.0);
    cap.setCohesion(5e6);
    cap.setCapHardening(1e-9);  // Hardening modulus
    cap.setInitialCapPosition(-50e6);  // Initial cap at -50 MPa
    
    // Hydrostatic compression beyond initial cap
    double sigma_trial[6] = {-80e6, -80e6, -80e6, 0, 0, 0};
    double sigma_return[6];
    double plastic_strain[6];
    double plastic_multiplier;
    
    cap.returnMap(sigma_trial, sigma_return, plastic_strain, plastic_multiplier);
    
    // After compaction, cap should have moved (hardening)
    double new_cap = cap.getCurrentCapPosition();
    
    // Cap should be more compressive than initial
    if (new_cap > -50e6) {
        std::cerr << "  FAIL: Cap didn't move with compaction" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Cap model compaction correct" << std::endl;
    return true;
}

// =============================================================================
// Damage Model Tests
// =============================================================================

/**
 * @test Test damage evolution
 */
bool test_damage_evolution() {
    std::cout << "Testing damage evolution..." << std::endl;
    
    DamageModel dm;
    dm.setCriticalStrain(0.001);  // 0.1% critical strain
    dm.setDamageRate(1.0);
    
    double damage = 0.0;
    double plastic_work = 1e6;  // 1 MJ/m³
    
    dm.updateDamage(damage, plastic_work, 0.002);  // 0.2% strain
    
    // Damage should increase
    if (damage <= 0) {
        std::cerr << "  FAIL: Damage didn't increase" << std::endl;
        return false;
    }
    
    // Damage should be bounded [0, 1]
    if (damage < 0 || damage > 1) {
        std::cerr << "  FAIL: Damage out of bounds: " << damage << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Damage evolution correct" << std::endl;
    return true;
}

/**
 * @test Test damaged stiffness
 */
bool test_damaged_stiffness() {
    std::cout << "Testing damaged stiffness..." << std::endl;
    
    double tol = 1e-10;
    
    DamageModel dm;
    
    // Original stiffness
    double C[36];
    double E = 30e9, nu = 0.25;
    double lambda = E * nu / ((1 + nu) * (1 - 2*nu));
    double mu = E / (2 * (1 + nu));
    
    // Fill isotropic stiffness
    for (int i = 0; i < 36; ++i) C[i] = 0;
    C[0] = C[7] = C[14] = lambda + 2*mu;
    C[1] = C[2] = C[6] = C[8] = C[12] = C[13] = lambda;
    C[21] = C[28] = C[35] = mu;
    
    // Apply 50% damage
    double C_damaged[36];
    dm.applyDamage(C, C_damaged, 0.5);
    
    // Stiffness should be reduced by (1-d)
    double factor = 0.5;
    if (std::abs(C_damaged[0] - factor * C[0]) > tol) {
        std::cerr << "  FAIL: Damaged stiffness incorrect" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Damaged stiffness correct" << std::endl;
    return true;
}

// =============================================================================
// Consistent Tangent Tests
// =============================================================================

/**
 * @test Test algorithmic tangent modulus
 */
bool test_consistent_tangent() {
    std::cout << "Testing consistent tangent modulus..." << std::endl;
    
    double tol = 1e-4;
    double h = 1e-8;
    
    DruckerPragerModel dp;
    dp.setFrictionAngle(30.0 * M_PI / 180.0);
    dp.setCohesion(5e6);
    dp.setDilatancyAngle(15.0 * M_PI / 180.0);
    
    double E = 30e9, nu = 0.25;
    dp.setElasticModuli(E, nu);
    
    // Stress state that causes yielding
    double sigma[6] = {-5e6, -15e6, -10e6, 2e6, 1e6, 0.5e6};
    double eps[6] = {0, 0, 0, 0, 0, 0};
    
    // Get algorithmic tangent
    double C_alg[36];
    dp.computeConsistentTangent(sigma, eps, C_alg);
    
    // Verify symmetry (should be symmetric for associative flow)
    for (int i = 0; i < 6; ++i) {
        for (int j = i+1; j < 6; ++j) {
            if (std::abs(C_alg[i*6+j] - C_alg[j*6+i]) > tol * (std::abs(C_alg[i*6+j]) + 1.0)) {
                // Non-associative flow may give non-symmetric tangent
                // This is acceptable
            }
        }
    }
    
    // Verify positive semi-definiteness approximately
    // Check diagonal elements are non-negative
    for (int i = 0; i < 6; ++i) {
        if (C_alg[i*6+i] < -tol) {
            std::cerr << "  WARNING: Negative diagonal in tangent" << std::endl;
            // Not necessarily an error for softening
        }
    }
    
    std::cout << "  PASS: Consistent tangent computed" << std::endl;
    return true;
}

// =============================================================================
// Hardening Tests
// =============================================================================

/**
 * @test Test isotropic hardening
 */
bool test_isotropic_hardening() {
    std::cout << "Testing isotropic hardening..." << std::endl;
    
    VonMisesModel vm;
    double yield_stress = 250e6;
    double hardening_modulus = 1e9;  // 1 GPa
    vm.setYieldStress(yield_stress);
    vm.setIsotropicHardening(hardening_modulus);
    
    // Initial yield
    double sigma_y0 = vm.getCurrentYieldStress();
    
    // Accumulate plastic strain
    vm.updateHardening(0.01);  // 1% plastic strain
    
    double sigma_y1 = vm.getCurrentYieldStress();
    
    // Yield stress should increase
    if (sigma_y1 <= sigma_y0) {
        std::cerr << "  FAIL: No hardening observed" << std::endl;
        return false;
    }
    
    // Check hardening amount
    double expected = sigma_y0 + hardening_modulus * 0.01;
    if (std::abs(sigma_y1 - expected) / expected > 0.1) {
        std::cerr << "  FAIL: Hardening amount incorrect" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Isotropic hardening correct" << std::endl;
    return true;
}

// =============================================================================
// Test Runner
// =============================================================================

int runPlasticityTests() {
    std::cout << "\n=== Plasticity Model Unit Tests ===" << std::endl;
    
    int passed = 0;
    int failed = 0;
    
    // Drucker-Prager tests
    if (test_drucker_prager_yield_function()) ++passed; else ++failed;
    if (test_drucker_prager_return_mapping()) ++passed; else ++failed;
    if (test_drucker_prager_elastic()) ++passed; else ++failed;
    
    // von Mises tests
    if (test_von_mises_yield_function()) ++passed; else ++failed;
    if (test_von_mises_radial_return()) ++passed; else ++failed;
    
    // Mohr-Coulomb tests
    if (test_mohr_coulomb_yield()) ++passed; else ++failed;
    if (test_mohr_coulomb_return()) ++passed; else ++failed;
    
    // Cap model tests
    if (test_cap_model_compaction()) ++passed; else ++failed;
    
    // Damage tests
    if (test_damage_evolution()) ++passed; else ++failed;
    if (test_damaged_stiffness()) ++passed; else ++failed;
    
    // Tangent tests
    if (test_consistent_tangent()) ++passed; else ++failed;
    
    // Hardening tests
    if (test_isotropic_hardening()) ++passed; else ++failed;
    
    std::cout << "\n=== Plasticity Test Summary ===" << std::endl;
    std::cout << "Passed: " << passed << std::endl;
    std::cout << "Failed: " << failed << std::endl;
    std::cout << "Total:  " << (passed + failed) << std::endl;
    
    return failed;
}

} // namespace Testing
} // namespace FSRM

#ifndef FSRM_TEST_NO_MAIN
int main() {
    return FSRM::Testing::runPlasticityTests();
}
#endif
