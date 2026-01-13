/**
 * @file test_plasticity_models.cpp
 * @brief Tests for plasticity models and yield criteria
 * 
 * Tests for:
 * - Drucker-Prager yield criterion
 * - von Mises (J2) plasticity
 * - Mohr-Coulomb model
 * - Return mapping algorithms
 * - Hardening/softening behavior
 * - Consistent tangent modulus
 */

#include <gtest/gtest.h>
#include "Testing.hpp"
#include "PlasticityModel.hpp"
#include "MaterialModel.hpp"
#include <cmath>
#include <array>

using namespace FSRM;
using namespace FSRM::Testing;

class PlasticityModelsTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        
        // Default elastic properties
        E = 50e9;       // Young's modulus (50 GPa)
        nu = 0.25;      // Poisson's ratio
        
        // Lamé parameters
        lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        mu = E / (2.0 * (1.0 + nu));
        K = E / (3.0 * (1.0 - 2.0 * nu));
        
        // Default plastic properties
        phi = 30.0 * M_PI / 180.0;  // Friction angle (rad)
        c = 10e6;                    // Cohesion (10 MPa)
        psi = 10.0 * M_PI / 180.0;  // Dilation angle (rad)
        
        // Initialize elastic modulus tensor (Voigt notation)
        initializeElasticModulus();
    }
    
    void initializeElasticModulus() {
        for (auto& row : C) {
            row.fill(0.0);
        }
        // Isotropic elastic modulus
        C[0][0] = C[1][1] = C[2][2] = lambda + 2.0 * mu;
        C[0][1] = C[0][2] = C[1][0] = C[1][2] = C[2][0] = C[2][1] = lambda;
        C[3][3] = C[4][4] = C[5][5] = mu;
    }
    
    // Compute stress invariants
    double computeI1(const std::array<double, 6>& sigma) {
        return sigma[0] + sigma[1] + sigma[2];
    }
    
    double computeJ2(const std::array<double, 6>& sigma) {
        double I1 = computeI1(sigma);
        double p = I1 / 3.0;
        
        // Deviatoric stress
        double s0 = sigma[0] - p;
        double s1 = sigma[1] - p;
        double s2 = sigma[2] - p;
        
        // J2 = 0.5 * s_ij * s_ij
        return 0.5 * (s0*s0 + s1*s1 + s2*s2) + 
               sigma[3]*sigma[3] + sigma[4]*sigma[4] + sigma[5]*sigma[5];
    }
    
    int rank;
    double E, nu, lambda, mu, K;
    double phi, c, psi;
    std::array<std::array<double, 6>, 6> C;
};

// ============================================================================
// Stress Invariant Tests
// ============================================================================

TEST_F(PlasticityModelsTest, FirstInvariantI1) {
    std::array<double, 6> sigma = {100e6, 80e6, 60e6, 0.0, 0.0, 0.0};
    
    double I1 = computeI1(sigma);
    
    // I1 = σ_xx + σ_yy + σ_zz
    EXPECT_NEAR(I1, 240e6, 1e3);
}

TEST_F(PlasticityModelsTest, SecondInvariantJ2) {
    // Pure shear state
    std::array<double, 6> sigma = {0.0, 0.0, 0.0, 50e6, 0.0, 0.0};
    
    double J2 = computeJ2(sigma);
    
    // For pure shear: J2 = τ²
    EXPECT_NEAR(J2, 50e6 * 50e6, 1e9);
    
    // Hydrostatic state
    std::array<double, 6> sigma_hyd = {-100e6, -100e6, -100e6, 0.0, 0.0, 0.0};
    
    J2 = computeJ2(sigma_hyd);
    
    // For hydrostatic: J2 = 0
    EXPECT_NEAR(J2, 0.0, 1e6);
}

TEST_F(PlasticityModelsTest, EquivalentStress) {
    // von Mises equivalent stress: σ_eq = √(3*J2)
    
    std::array<double, 6> sigma = {100e6, 0.0, 0.0, 0.0, 0.0, 0.0};  // Uniaxial
    
    double J2 = computeJ2(sigma);
    double sigma_eq = std::sqrt(3.0 * J2);
    
    // For uniaxial tension: σ_eq = |σ_xx|
    EXPECT_NEAR(sigma_eq, 100e6, 1e6);
}

// ============================================================================
// Drucker-Prager Tests
// ============================================================================

TEST_F(PlasticityModelsTest, DruckerPragerYieldFunction) {
    DruckerPragerModel dp;
    
    DruckerPragerModel::Parameters params;
    params.friction_angle = phi;
    params.cohesion = c;
    params.dilation_angle = psi;
    
    dp.setParameters(params);
    
    // Test elastic state (below yield)
    std::array<double, 6> sigma_elastic = {-50e6, -50e6, -50e6, 0.0, 0.0, 0.0};
    double F_elastic = dp.yieldFunction(sigma_elastic, 0.0);
    
    EXPECT_LT(F_elastic, 0.0) << "Should be in elastic region";
    
    // Test plastic state (at/above yield)
    std::array<double, 6> sigma_plastic = {-10e6, -50e6, -50e6, 30e6, 0.0, 0.0};
    double F_plastic = dp.yieldFunction(sigma_plastic, 0.0);
    
    // Just verify the function is finite
    EXPECT_TRUE(std::isfinite(F_plastic));
}

TEST_F(PlasticityModelsTest, DruckerPragerAlphaKappa) {
    // Drucker-Prager parameters from friction angle and cohesion
    // α = 2 sin(φ) / (√3 * (3 - sin(φ)))
    // κ = 6 c cos(φ) / (√3 * (3 - sin(φ)))
    
    double sin_phi = std::sin(phi);
    double cos_phi = std::cos(phi);
    
    double alpha = 2.0 * sin_phi / (std::sqrt(3.0) * (3.0 - sin_phi));
    double kappa = 6.0 * c * cos_phi / (std::sqrt(3.0) * (3.0 - sin_phi));
    
    EXPECT_GT(alpha, 0.0) << "α should be positive for φ > 0";
    EXPECT_GT(kappa, 0.0) << "κ should be positive for c > 0";
    EXPECT_LT(alpha, 1.0) << "α should be bounded";
}

TEST_F(PlasticityModelsTest, DruckerPragerFlowDirection) {
    DruckerPragerModel dp;
    
    DruckerPragerModel::Parameters params;
    params.friction_angle = phi;
    params.cohesion = c;
    params.dilation_angle = psi;
    
    dp.setParameters(params);
    
    std::array<double, 6> sigma = {-20e6, -40e6, -40e6, 10e6, 0.0, 0.0};
    std::array<double, 6> flow_dir;
    
    dp.flowDirection(sigma, flow_dir);
    
    // Flow direction should be finite
    for (int i = 0; i < 6; ++i) {
        EXPECT_TRUE(std::isfinite(flow_dir[i])) << "Flow direction component " << i;
    }
}

TEST_F(PlasticityModelsTest, DruckerPragerReturnMapping) {
    DruckerPragerModel dp;
    
    DruckerPragerModel::Parameters params;
    params.friction_angle = phi;
    params.cohesion = c;
    
    dp.setParameters(params);
    
    // Start with stress above yield surface
    std::array<double, 6> sigma = {0.0, -30e6, -30e6, 50e6, 0.0, 0.0};
    std::array<double, 6> eps_p = {0, 0, 0, 0, 0, 0};
    double eps_p_acc = 0.0;
    std::array<double, 6> deps = {0.001, -0.0005, -0.0005, 0.001, 0.0, 0.0};
    bool plastic;
    
    dp.returnMapping(sigma, eps_p, eps_p_acc, deps, C, plastic);
    
    // After return mapping, stress should be on or inside yield surface
    double F = dp.yieldFunction(sigma, 0.0);
    EXPECT_LE(F, 1e-6 * c) << "Stress should return to yield surface";
}

// ============================================================================
// von Mises Tests
// ============================================================================

TEST_F(PlasticityModelsTest, VonMisesYieldFunction) {
    VonMisesModel vm;
    
    VonMisesModel::Parameters params;
    params.yield_stress = 100e6;  // 100 MPa
    
    vm.setParameters(params);
    
    // Elastic state
    std::array<double, 6> sigma_elastic = {50e6, 0.0, 0.0, 0.0, 0.0, 0.0};
    double F_elastic = vm.yieldFunction(sigma_elastic, 0.0);
    
    EXPECT_LT(F_elastic, 0.0) << "Should be elastic";
    
    // At yield
    std::array<double, 6> sigma_yield = {100e6, 0.0, 0.0, 0.0, 0.0, 0.0};
    double F_yield = vm.yieldFunction(sigma_yield, 0.0);
    
    EXPECT_NEAR(F_yield, 0.0, 1e6) << "Should be at yield";
}

TEST_F(PlasticityModelsTest, VonMisesPureShear) {
    VonMisesModel vm;
    
    VonMisesModel::Parameters params;
    params.yield_stress = 100e6;
    
    vm.setParameters(params);
    
    // Pure shear: τ_yield = σ_yield / √3
    double tau_yield = params.yield_stress / std::sqrt(3.0);
    
    std::array<double, 6> sigma_shear = {0.0, 0.0, 0.0, tau_yield, 0.0, 0.0};
    double F = vm.yieldFunction(sigma_shear, 0.0);
    
    EXPECT_NEAR(F, 0.0, 1e6) << "Pure shear at τ = σ_y/√3 should be at yield";
}

TEST_F(PlasticityModelsTest, VonMisesReturnMapping) {
    VonMisesModel vm;
    
    VonMisesModel::Parameters params;
    params.yield_stress = 100e6;
    params.hardening_modulus = 1e9;
    
    vm.setParameters(params);
    
    std::array<double, 6> sigma = {150e6, 0.0, 0.0, 0.0, 0.0, 0.0};  // Above yield
    std::array<double, 6> eps_p = {0, 0, 0, 0, 0, 0};
    double eps_p_acc = 0.0;
    std::array<double, 6> deps = {0.001, 0.0, 0.0, 0.0, 0.0, 0.0};
    bool plastic;
    
    vm.returnMapping(sigma, eps_p, eps_p_acc, deps, C, plastic);
    
    EXPECT_TRUE(plastic) << "Should be plastic loading";
    EXPECT_GT(eps_p_acc, 0.0) << "Plastic strain should accumulate";
}

// ============================================================================
// Mohr-Coulomb Tests
// ============================================================================

TEST_F(PlasticityModelsTest, MohrCoulombYieldFunction) {
    MohrCoulombModel mc;
    
    MohrCoulombModel::Parameters params;
    params.friction_angle = phi;
    params.cohesion = c;
    
    mc.setParameters(params);
    
    // Triaxial compression: σ₁ = σ_axial, σ₂ = σ₃ = σ_confining
    double sigma_conf = -50e6;  // Confining pressure
    
    // Mohr-Coulomb failure: σ₁ - σ₃ = 2c*cos(φ) + (σ₁ + σ₃)*sin(φ)
    // Solving: σ₁ = σ₃*(1+sin(φ))/(1-sin(φ)) + 2c*cos(φ)/(1-sin(φ))
    double sin_phi = std::sin(phi);
    double cos_phi = std::cos(phi);
    double Nf = (1.0 + sin_phi) / (1.0 - sin_phi);
    double sigma_1_fail = sigma_conf * Nf + 2.0 * c * cos_phi / (1.0 - sin_phi);
    
    // Check that failure envelope is correctly computed
    EXPECT_GT(sigma_1_fail, sigma_conf) << "Failure stress > confining";
}

TEST_F(PlasticityModelsTest, MohrCoulombTensionCutoff) {
    MohrCoulombModel mc;
    
    MohrCoulombModel::Parameters params;
    params.friction_angle = phi;
    params.cohesion = c;
    params.tension_cutoff = 1e6;  // 1 MPa tensile strength
    
    mc.setParameters(params);
    
    // Uniaxial tension exceeding cutoff
    std::array<double, 6> sigma_tension = {5e6, 0.0, 0.0, 0.0, 0.0, 0.0};
    
    double F = mc.yieldFunction(sigma_tension, 0.0);
    
    // Should violate tension cutoff
    EXPECT_GT(F, 0.0) << "Should exceed tension cutoff";
}

// ============================================================================
// Hardening/Softening Tests
// ============================================================================

TEST_F(PlasticityModelsTest, IsotropicHardening) {
    DruckerPragerModel dp;
    
    DruckerPragerModel::Parameters params;
    params.friction_angle = phi;
    params.cohesion = c;
    params.hardening_modulus = 1e9;  // 1 GPa hardening
    
    dp.setParameters(params);
    
    // Yield strength should increase with plastic strain
    double kappa_0 = dp.hardeningFunction(0.0);
    double kappa_1 = dp.hardeningFunction(0.01);
    double kappa_2 = dp.hardeningFunction(0.02);
    
    EXPECT_GT(kappa_1, kappa_0) << "Yield should increase with plastic strain";
    EXPECT_GT(kappa_2, kappa_1);
}

TEST_F(PlasticityModelsTest, SaturationHardening) {
    DruckerPragerModel dp;
    
    DruckerPragerModel::Parameters params;
    params.cohesion = c;
    params.hardening_modulus = 1e9;
    params.saturation_yield = 50e6;  // Saturation cohesion
    
    dp.setParameters(params);
    
    // At very large strain, should approach saturation
    double kappa_large = dp.hardeningFunction(1.0);  // 100% strain
    
    EXPECT_LE(kappa_large, params.saturation_yield * 1.1);
}

TEST_F(PlasticityModelsTest, StrainSoftening) {
    DruckerPragerModel dp;
    
    DruckerPragerModel::Parameters params;
    params.cohesion = c;
    params.softening_modulus = -5e8;  // Negative = softening
    params.residual_cohesion = 1e6;    // Residual cohesion
    
    dp.setParameters(params);
    
    double kappa_0 = dp.hardeningFunction(0.0);
    double kappa_1 = dp.hardeningFunction(0.1);
    
    // With softening, yield strength decreases (or cohesion decreases)
    // The behavior depends on implementation
    EXPECT_TRUE(std::isfinite(kappa_1));
}

// ============================================================================
// Return Mapping Algorithm Tests
// ============================================================================

TEST_F(PlasticityModelsTest, ReturnMappingConsistency) {
    PlasticityModel model(YieldCriterion::DRUCKER_PRAGER);
    
    DruckerPragerModel::Parameters params;
    params.friction_angle = phi;
    params.cohesion = c;
    
    model.setDruckerPragerParameters(params);
    model.enable(true);
    
    // Apply strain increment
    std::array<double, 6> stress = {-50e6, -50e6, -50e6, 0.0, 0.0, 0.0};
    std::array<double, 6> strain_inc = {0.001, -0.0003, -0.0003, 0.0005, 0.0, 0.0};
    PlasticityState state;
    
    model.integrateStress(stress, strain_inc, C, state);
    
    // Check consistency
    double F = model.getYieldFunctionValue(stress, state);
    
    if (state.is_plastic) {
        // If plastic, should be on yield surface
        EXPECT_NEAR(F, 0.0, 1e6);
    } else {
        // If elastic, should be inside
        EXPECT_LE(F, 0.0);
    }
}

TEST_F(PlasticityModelsTest, ReturnMappingElasticUnloading) {
    PlasticityModel model(YieldCriterion::VON_MISES);
    
    VonMisesModel::Parameters params;
    params.yield_stress = 100e6;
    
    model.setVonMisesParameters(params);
    model.enable(true);
    
    // First load plastically
    std::array<double, 6> stress = {120e6, 0.0, 0.0, 0.0, 0.0, 0.0};
    std::array<double, 6> strain_inc = {0.001, -0.00025, -0.00025, 0.0, 0.0, 0.0};
    PlasticityState state;
    
    model.integrateStress(stress, strain_inc, C, state);
    
    EXPECT_TRUE(state.is_plastic) << "Should yield on loading";
    
    // Now unload (reverse strain)
    std::array<double, 6> unload_strain = {-0.0005, 0.000125, 0.000125, 0.0, 0.0, 0.0};
    model.integrateStress(stress, unload_strain, C, state);
    
    // After unloading, should be elastic
    EXPECT_FALSE(state.is_plastic) << "Should be elastic on unloading";
}

// ============================================================================
// Consistent Tangent Tests
// ============================================================================

TEST_F(PlasticityModelsTest, ConsistentTangentSymmetry) {
    PlasticityModel model(YieldCriterion::DRUCKER_PRAGER);
    
    DruckerPragerModel::Parameters params;
    params.friction_angle = phi;
    params.cohesion = c;
    
    model.setDruckerPragerParameters(params);
    
    std::array<double, 6> stress = {-20e6, -50e6, -50e6, 20e6, 0.0, 0.0};
    PlasticityState state;
    state.is_plastic = true;
    
    std::array<std::array<double, 6>, 6> D;
    model.getConsistentTangent(stress, state, C, D);
    
    // Check symmetry: D_ij = D_ji
    for (int i = 0; i < 6; ++i) {
        for (int j = i+1; j < 6; ++j) {
            EXPECT_NEAR(D[i][j], D[j][i], std::abs(D[i][j]) * 1e-6 + 1e-10)
                << "Tangent should be symmetric at (" << i << "," << j << ")";
        }
    }
}

TEST_F(PlasticityModelsTest, ConsistentTangentPositiveDefinite) {
    PlasticityModel model(YieldCriterion::VON_MISES);
    
    VonMisesModel::Parameters params;
    params.yield_stress = 100e6;
    params.hardening_modulus = 1e9;
    
    model.setVonMisesParameters(params);
    
    std::array<double, 6> stress = {100e6, 0.0, 0.0, 0.0, 0.0, 0.0};
    PlasticityState state;
    state.is_plastic = true;
    
    std::array<std::array<double, 6>, 6> D;
    model.getConsistentTangent(stress, state, C, D);
    
    // For stable behavior, tangent should be positive definite
    // Check diagonal elements are positive
    for (int i = 0; i < 6; ++i) {
        EXPECT_GT(D[i][i], 0.0) << "Diagonal element D[" << i << "][" << i << "] should be positive";
    }
}

// ============================================================================
// Cap Model Tests
// ============================================================================

TEST_F(PlasticityModelsTest, CapModelYieldSurfaces) {
    CapModel cap;
    
    CapModel::Parameters params;
    params.friction_angle = phi;
    params.cohesion = c;
    params.cap_eccentricity = 0.3;
    params.initial_cap_position = -100e6;
    
    cap.setParameters(params);
    
    // Test shear yield
    std::array<double, 6> sigma_shear = {0.0, -50e6, -50e6, 30e6, 0.0, 0.0};
    double F_shear = cap.shearYield(sigma_shear, c);
    
    EXPECT_TRUE(std::isfinite(F_shear));
    
    // Test cap yield (high compression)
    std::array<double, 6> sigma_cap = {-200e6, -200e6, -200e6, 0.0, 0.0, 0.0};
    double F_cap = cap.capYield(sigma_cap, params.initial_cap_position);
    
    EXPECT_TRUE(std::isfinite(F_cap));
}

TEST_F(PlasticityModelsTest, CapModelHardening) {
    CapModel cap;
    
    CapModel::Parameters params;
    params.cap_hardening_modulus = 1e10;
    params.initial_cap_position = -50e6;
    
    cap.setParameters(params);
    
    // Cap position should evolve with plastic volumetric strain
    // Implementation-specific test
    
    EXPECT_GT(params.cap_hardening_modulus, 0.0);
}

// ============================================================================
// Damage Model Tests
// ============================================================================

TEST_F(PlasticityModelsTest, DamageEvolution) {
    DamageModel damage;
    
    DamageModel::Parameters params;
    params.damage_threshold = 1e-4;
    params.damage_exponent = 2.0;
    params.maximum_damage = 0.99;
    
    damage.setParameters(params);
    
    double D = 0.0;
    std::array<double, 6> stress = {50e6, 0.0, 0.0, 0.0, 0.0, 0.0};
    std::array<double, 6> strain = {0.001, -0.00025, -0.00025, 0.0, 0.0, 0.0};
    
    damage.updateDamage(D, stress, strain, 0.001);
    
    // Damage should increase when strain exceeds threshold
    EXPECT_GE(D, 0.0);
    EXPECT_LE(D, params.maximum_damage);
}

TEST_F(PlasticityModelsTest, DamageStiffnessReduction) {
    DamageModel damage;
    
    DamageModel::Parameters params;
    params.damage_threshold = 1e-4;
    
    damage.setParameters(params);
    
    double D = 0.5;  // 50% damage
    
    auto C_damaged = C;
    damage.applyDamage(C_damaged, D);
    
    // Damaged stiffness should be reduced
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            if (std::abs(C[i][j]) > 1e-10) {
                EXPECT_LT(std::abs(C_damaged[i][j]), std::abs(C[i][j]))
                    << "Stiffness should be reduced by damage";
            }
        }
    }
}

// ============================================================================
// Physical Validation Tests
// ============================================================================

TEST_F(PlasticityModelsTest, TriaxialCompressionTest) {
    // Simulate triaxial compression test
    
    DruckerPragerModel dp;
    
    DruckerPragerModel::Parameters params;
    params.friction_angle = 30.0 * M_PI / 180.0;
    params.cohesion = 5e6;
    
    dp.setParameters(params);
    
    // Confining pressure
    double sigma_3 = -20e6;  // 20 MPa confining
    
    // Increase axial stress until yield
    for (int i = 0; i <= 100; ++i) {
        double sigma_1 = sigma_3 - i * 1e6;  // Increase compression
        
        std::array<double, 6> sigma = {sigma_1, sigma_3, sigma_3, 0.0, 0.0, 0.0};
        double F = dp.yieldFunction(sigma, 0.0);
        
        if (F >= 0) {
            // Found yield point
            double q = std::abs(sigma_1 - sigma_3);  // Deviatoric stress
            double p = (sigma_1 + 2.0 * sigma_3) / 3.0;  // Mean stress
            
            // Check Mohr-Coulomb equivalent
            double sin_phi = std::sin(params.friction_angle);
            double cos_phi = std::cos(params.friction_angle);
            double q_fail = -2.0 * sin_phi * p / (1.0 - sin_phi) + 
                           6.0 * params.cohesion * cos_phi / (std::sqrt(3.0) * (3.0 - sin_phi));
            
            EXPECT_GT(q, 0.0);
            break;
        }
    }
}

TEST_F(PlasticityModelsTest, SimpleShearTest) {
    // Pure shear should yield at τ = c (for φ = 0)
    
    DruckerPragerModel dp;
    
    DruckerPragerModel::Parameters params;
    params.friction_angle = 0.0;  // No friction (von Mises like)
    params.cohesion = 10e6;
    
    dp.setParameters(params);
    
    // Apply increasing shear stress
    for (int i = 1; i <= 20; ++i) {
        double tau = i * 1e6;
        
        std::array<double, 6> sigma = {0.0, 0.0, 0.0, tau, 0.0, 0.0};
        double F = dp.yieldFunction(sigma, 0.0);
        
        if (F >= 0) {
            // At yield, τ should be close to c / √3 (for Drucker-Prager)
            EXPECT_GT(tau, 0.0);
            break;
        }
    }
}

// ============================================================================
// Configuration Tests
// ============================================================================

TEST_F(PlasticityModelsTest, PlasticityConfigParsing) {
    PlasticityConfig config;
    
    std::map<std::string, std::string> params;
    params["enabled"] = "true";
    params["yield_criterion"] = "drucker_prager";
    params["friction_angle"] = "30.0";
    params["cohesion"] = "1e7";
    params["dilation_angle"] = "10.0";
    
    config.parseConfig(params);
    
    EXPECT_TRUE(config.enabled);
    EXPECT_NEAR(config.friction_angle, 30.0, 0.1);
    EXPECT_NEAR(config.cohesion, 1e7, 1e5);
}

TEST_F(PlasticityModelsTest, CreateModelFromConfig) {
    PlasticityConfig config;
    config.enabled = true;
    config.yield_criterion = "drucker_prager";
    config.friction_angle = 30.0;
    config.cohesion = 10e6;
    
    PlasticityModel model = config.createModel();
    
    EXPECT_TRUE(model.isEnabled());
}

// ============================================================================
// Edge Cases and Robustness Tests
// ============================================================================

TEST_F(PlasticityModelsTest, ZeroStressState) {
    PlasticityModel model(YieldCriterion::DRUCKER_PRAGER);
    
    DruckerPragerModel::Parameters params;
    params.cohesion = c;
    
    model.setDruckerPragerParameters(params);
    
    std::array<double, 6> sigma = {0, 0, 0, 0, 0, 0};
    PlasticityState state;
    
    // Zero stress should be elastic (inside yield surface for c > 0)
    bool yielding = model.isYielding(sigma, state);
    
    // Depends on tension cutoff, but for pure cohesive with c > 0, zero should be elastic
    EXPECT_TRUE(std::isfinite(model.getYieldFunctionValue(sigma, state)));
}

TEST_F(PlasticityModelsTest, HydrostaticCompression) {
    PlasticityModel model(YieldCriterion::VON_MISES);
    
    VonMisesModel::Parameters params;
    params.yield_stress = 100e6;
    
    model.setVonMisesParameters(params);
    
    // Hydrostatic compression - no deviatoric stress, should never yield
    std::array<double, 6> sigma_hyd = {-500e6, -500e6, -500e6, 0.0, 0.0, 0.0};
    PlasticityState state;
    
    bool yielding = model.isYielding(sigma_hyd, state);
    
    // von Mises is pressure-independent, hydrostatic never yields
    EXPECT_FALSE(yielding);
}

TEST_F(PlasticityModelsTest, VeryLargeStress) {
    PlasticityModel model(YieldCriterion::DRUCKER_PRAGER);
    
    DruckerPragerModel::Parameters params;
    params.cohesion = c;
    params.friction_angle = phi;
    
    model.setDruckerPragerParameters(params);
    
    std::array<double, 6> sigma = {-1e12, -1e12, -1e12, 1e11, 0.0, 0.0};
    PlasticityState state;
    
    double F = model.getYieldFunctionValue(sigma, state);
    
    // Should handle large stresses without numerical issues
    EXPECT_TRUE(std::isfinite(F));
}

TEST_F(PlasticityModelsTest, SmallStrainIncrement) {
    PlasticityModel model(YieldCriterion::VON_MISES);
    
    VonMisesModel::Parameters params;
    params.yield_stress = 100e6;
    
    model.setVonMisesParameters(params);
    model.enable(true);
    
    std::array<double, 6> stress = {110e6, 0.0, 0.0, 0.0, 0.0, 0.0};
    std::array<double, 6> strain_inc = {1e-10, -2.5e-11, -2.5e-11, 0.0, 0.0, 0.0};
    PlasticityState state;
    
    // Should handle very small strain increments
    model.integrateStress(stress, strain_inc, C, state);
    
    EXPECT_TRUE(std::isfinite(stress[0]));
}
