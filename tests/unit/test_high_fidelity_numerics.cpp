/**
 * @file test_high_fidelity_numerics.cpp
 * @brief Unit tests for HighFidelityNumerics module
 * 
 * Tests cover:
 * - SUPG/GLS stabilization parameter calculations
 * - PETSc time integration wrapper
 * - p-adaptivity smoothness indicators
 * - Physics-based preconditioning setup
 */

#include <gtest/gtest.h>
#include "HighFidelityNumerics.hpp"
#include <cmath>

using namespace FSRM::HighFidelity;

// =============================================================================
// SUPG Stabilization Tests
// =============================================================================

class SUPGStabilizationTest : public ::testing::Test {
protected:
    SUPGStabilization supg;
    
    void SetUp() override {
        SUPGStabilization::Parameters params;
        params.method = StabilizationMethod::SUPG;
        params.tau_formula = TauFormula::SHAKIB;
        params.tau_multiplier = 1.0;
        params.include_transient = true;
        params.include_reaction = true;
        supg.setParameters(params);
    }
};

TEST_F(SUPGStabilizationTest, TauShakibFormula) {
    // Test case: advection-dominated flow
    PetscReal velocity[3] = {1.0, 0.0, 0.0};
    PetscReal diffusivity = 0.001;  // Small diffusion
    PetscReal reaction = 0.0;
    PetscReal h = 0.1;
    PetscReal dt = 0.01;
    
    PetscReal tau = supg.calculateTau(velocity, diffusivity, reaction, h, dt);
    
    // τ should be positive and finite
    EXPECT_GT(tau, 0.0);
    EXPECT_LT(tau, 1e10);
    
    // For advection-dominated: τ ≈ h / (2|v|)
    PetscReal expected_tau_advection = h / (2.0 * 1.0);
    EXPECT_NEAR(tau, expected_tau_advection, expected_tau_advection * 0.5);
}

TEST_F(SUPGStabilizationTest, TauDiffusionDominated) {
    // Test case: diffusion-dominated flow
    PetscReal velocity[3] = {0.001, 0.0, 0.0};
    PetscReal diffusivity = 1.0;  // Large diffusion
    PetscReal reaction = 0.0;
    PetscReal h = 0.1;
    PetscReal dt = 0.0;  // Steady-state
    
    PetscReal tau = supg.calculateTau(velocity, diffusivity, reaction, h, dt);
    
    // τ should be smaller for diffusion-dominated
    EXPECT_GT(tau, 0.0);
}

TEST_F(SUPGStabilizationTest, PecletNumber) {
    // Test Peclet number calculation
    PetscReal Pe = supg.pecletNumber(1.0, 0.1, 0.01);  // v=1, h=0.1, κ=0.01
    
    // Pe = |v|·h / (2κ) = 1·0.1 / (2·0.01) = 5
    EXPECT_NEAR(Pe, 5.0, 1e-10);
}

TEST_F(SUPGStabilizationTest, PecletNumberZeroDiffusivity) {
    // With zero diffusivity, Pe should be very large
    PetscReal Pe = supg.pecletNumber(1.0, 0.1, 0.0);
    EXPECT_GT(Pe, 1e6);
}

TEST_F(SUPGStabilizationTest, OptimalDiffusivity) {
    // Test optimal artificial diffusivity
    PetscReal kappa_art = supg.optimalDiffusivity(1.0, 0.1, 0.001);
    
    // Should be positive for advection-dominated flow
    EXPECT_GT(kappa_art, 0.0);
    // Should be less than 0.5·|v|·h (maximum)
    EXPECT_LT(kappa_art, 0.5 * 1.0 * 0.1);
}

TEST_F(SUPGStabilizationTest, StabilizationTerm) {
    PetscReal velocity[3] = {1.0, 0.0, 0.0};
    PetscReal test_gradient[3] = {1.0, 0.0, 0.0};
    PetscReal residual_strong = 1.0;
    PetscReal tau = 0.05;
    PetscReal stab_term = 0.0;
    
    supg.addStabilization(velocity, test_gradient, residual_strong, tau, &stab_term);
    
    // stab = τ · (v·∇φ) · R = 0.05 · 1.0 · 1.0 = 0.05
    EXPECT_NEAR(stab_term, 0.05, 1e-10);
}

TEST_F(SUPGStabilizationTest, ShockCapturing) {
    SUPGStabilization::Parameters params;
    params.shock_capturing = true;
    params.shock_capture_coeff = 0.5;
    supg.setParameters(params);
    
    PetscReal grad_u[3] = {10.0, 0.0, 0.0};  // Large gradient
    PetscReal u = 1.0;
    PetscReal h = 0.1;
    PetscReal residual = 1.0;
    
    PetscReal kappa_dc = supg.shockCapturingDiffusivity(grad_u, u, h, residual);
    
    // Should be positive
    EXPECT_GT(kappa_dc, 0.0);
}

TEST_F(SUPGStabilizationTest, PSPGTerm) {
    PetscReal grad_p[3] = {1.0, 2.0, 0.0};
    PetscReal grad_q[3] = {1.0, 0.0, 0.0};
    PetscReal tau = 0.01;
    
    PetscReal pspg = supg.pspgTerm(grad_p, grad_q, tau);
    
    // pspg = τ · (∇p · ∇q) = 0.01 · (1·1 + 2·0) = 0.01
    EXPECT_NEAR(pspg, 0.01, 1e-10);
}

TEST_F(SUPGStabilizationTest, Configure) {
    std::map<std::string, std::string> config;
    config["supg_method"] = "gls";
    config["supg_tau_formula"] = "tezduyar";
    config["supg_tau_multiplier"] = "2.0";
    
    SUPGStabilization supg2;
    supg2.configure(config);
    
    // Just verify it doesn't crash - internal state not directly accessible
    PetscReal velocity[3] = {1.0, 0.0, 0.0};
    PetscReal tau = supg2.calculateTau(velocity, 0.01, 0.0, 0.1, 0.01);
    EXPECT_GT(tau, 0.0);
}


// =============================================================================
// PETSc Time Integration Tests
// =============================================================================

class PETScTimeIntegrationTest : public ::testing::Test {
protected:
    MPI_Comm comm;
    
    void SetUp() override {
        comm = PETSC_COMM_WORLD;
    }
};

TEST_F(PETScTimeIntegrationTest, GetOrderBDF) {
    PETScTimeIntegration ts(comm);
    PETScTimeIntegration::Parameters params;
    
    params.scheme = TimeScheme::BDF1;
    ts.setParameters(params);
    EXPECT_EQ(ts.getOrder(), 1);
    
    params.scheme = TimeScheme::BDF2;
    ts.setParameters(params);
    EXPECT_EQ(ts.getOrder(), 2);
    
    params.scheme = TimeScheme::BDF3;
    ts.setParameters(params);
    EXPECT_EQ(ts.getOrder(), 3);
    
    params.scheme = TimeScheme::BDF4;
    ts.setParameters(params);
    EXPECT_EQ(ts.getOrder(), 4);
}

TEST_F(PETScTimeIntegrationTest, GetOrderRK) {
    PETScTimeIntegration ts(comm);
    PETScTimeIntegration::Parameters params;
    
    params.scheme = TimeScheme::RK4;
    ts.setParameters(params);
    EXPECT_EQ(ts.getOrder(), 4);
    
    params.scheme = TimeScheme::SSPRK3;
    ts.setParameters(params);
    EXPECT_EQ(ts.getOrder(), 3);
}

TEST_F(PETScTimeIntegrationTest, StabilityPropertiesBDF) {
    PETScTimeIntegration ts(comm);
    PETScTimeIntegration::Parameters params;
    
    // BDF2 is A-stable and L-stable
    params.scheme = TimeScheme::BDF2;
    ts.setParameters(params);
    EXPECT_TRUE(ts.isAStable());
    EXPECT_TRUE(ts.isLStable());
    
    // Backward Euler is A-stable and L-stable
    params.scheme = TimeScheme::BACKWARD_EULER;
    ts.setParameters(params);
    EXPECT_TRUE(ts.isAStable());
    EXPECT_TRUE(ts.isLStable());
}

TEST_F(PETScTimeIntegrationTest, StabilityPropertiesRK) {
    PETScTimeIntegration ts(comm);
    PETScTimeIntegration::Parameters params;
    
    // Explicit RK4 is NOT A-stable
    params.scheme = TimeScheme::RK4;
    ts.setParameters(params);
    EXPECT_FALSE(ts.isAStable());
    EXPECT_FALSE(ts.isLStable());
}

TEST_F(PETScTimeIntegrationTest, StabilityPropertiesTheta) {
    PETScTimeIntegration ts(comm);
    PETScTimeIntegration::Parameters params;
    
    // Crank-Nicolson (theta=0.5) is A-stable but NOT L-stable
    params.scheme = TimeScheme::CRANK_NICOLSON;
    ts.setParameters(params);
    EXPECT_TRUE(ts.isAStable());
    // CN is actually not L-stable (but still A-stable)
}

TEST_F(PETScTimeIntegrationTest, Configure) {
    PETScTimeIntegration ts(comm);
    
    std::map<std::string, std::string> config;
    config["time_scheme"] = "bdf2";
    config["time_adaptive"] = "true";
    config["time_atol"] = "1e-8";
    config["time_rtol"] = "1e-6";
    
    ts.configure(config);
    
    // Verify order is 2 for BDF2
    EXPECT_EQ(ts.getOrder(), 2);
}


// =============================================================================
// p-Adaptivity Tests
// =============================================================================

class PAdaptivityTest : public ::testing::Test {
protected:
    PAdaptivityManager p_adapt;
    
    void SetUp() override {
        PAdaptivityManager::Parameters params;
        params.min_order = 1;
        params.max_order = 5;
        params.default_order = 2;
        params.smoothness_threshold_inc = 0.8;
        params.smoothness_threshold_dec = 0.2;
        params.p_refinement_threshold = 0.5;
        p_adapt.setParameters(params);
    }
};

TEST_F(PAdaptivityTest, DecideOrderChangeIncrease) {
    // High smoothness + high error -> increase order
    PetscInt change = p_adapt.decideOrderChange(0.9, 0.6, 2);
    EXPECT_EQ(change, 1);
}

TEST_F(PAdaptivityTest, DecideOrderChangeDecrease) {
    // Low smoothness -> decrease order
    PetscInt change = p_adapt.decideOrderChange(0.1, 0.3, 3);
    EXPECT_EQ(change, -1);
}

TEST_F(PAdaptivityTest, DecideOrderChangeKeep) {
    // Medium smoothness and error -> keep order
    PetscInt change = p_adapt.decideOrderChange(0.5, 0.3, 2);
    EXPECT_EQ(change, 0);
}

TEST_F(PAdaptivityTest, DecideOrderChangeAtMaxOrder) {
    // At max order, can't increase further
    PetscInt change = p_adapt.decideOrderChange(0.9, 0.6, 5);
    EXPECT_EQ(change, 0);
}

TEST_F(PAdaptivityTest, DecideOrderChangeAtMinOrder) {
    // At min order, can't decrease further
    PetscInt change = p_adapt.decideOrderChange(0.1, 0.3, 1);
    EXPECT_EQ(change, 0);
}

TEST_F(PAdaptivityTest, Configure) {
    std::map<std::string, std::string> config;
    config["p_adapt_min_order"] = "2";
    config["p_adapt_max_order"] = "8";
    config["p_adapt_smoothness_inc"] = "0.9";
    
    PAdaptivityManager p_adapt2;
    p_adapt2.configure(config);
    
    // Should not crash and accept configuration
    SUCCEED();
}


// =============================================================================
// Preconditioner Tests
// =============================================================================

class PreconditionerTest : public ::testing::Test {
protected:
    MPI_Comm comm;
    
    void SetUp() override {
        comm = PETSC_COMM_WORLD;
    }
};

TEST_F(PreconditionerTest, Configure) {
    PETScPhysicsPreconditioner pc(comm);
    
    std::map<std::string, std::string> config;
    config["precond_type"] = "fieldsplit_schur";
    config["precond_mg_levels"] = "4";
    config["precond_mg_smoother"] = "jacobi";
    
    pc.configure(config);
    
    // Should not crash
    SUCCEED();
}


// =============================================================================
// Configuration Validation Tests
// =============================================================================

class ConfigValidationTest : public ::testing::Test {
protected:
    HighFidelityNumericsConfig config;
};

TEST_F(ConfigValidationTest, ValidConfiguration) {
    config.enable_high_fidelity_numerics = true;
    config.enable_advanced_time = true;
    config.time_params.scheme = TimeScheme::BDF2;
    config.time_params.theta = 0.5;
    config.time_params.rho_infinity = 0.5;
    
    config.enable_p_adapt = true;
    config.p_adapt_params.min_order = 1;
    config.p_adapt_params.max_order = 5;
    
    std::string error_msg;
    PetscErrorCode ierr = config.validate(error_msg);
    
    EXPECT_EQ(ierr, PETSC_SUCCESS);
    EXPECT_TRUE(error_msg.empty());
}

TEST_F(ConfigValidationTest, InvalidThetaRange) {
    config.enable_advanced_time = true;
    config.time_params.theta = 1.5;  // Invalid: > 1
    
    std::string error_msg;
    PetscErrorCode ierr = config.validate(error_msg);
    
    EXPECT_NE(ierr, PETSC_SUCCESS);
    EXPECT_FALSE(error_msg.empty());
}

TEST_F(ConfigValidationTest, InvalidPAdaptOrder) {
    config.enable_p_adapt = true;
    config.p_adapt_params.min_order = 5;
    config.p_adapt_params.max_order = 3;  // Invalid: max < min
    
    std::string error_msg;
    PetscErrorCode ierr = config.validate(error_msg);
    
    EXPECT_NE(ierr, PETSC_SUCCESS);
    EXPECT_FALSE(error_msg.empty());
}

TEST_F(ConfigValidationTest, ParseConfig) {
    std::map<std::string, std::string> cfg;
    cfg["enable_high_fidelity_numerics"] = "true";
    cfg["enable_supg"] = "true";
    cfg["enable_advanced_time"] = "true";
    cfg["enable_p_adapt"] = "false";
    cfg["enable_physics_precond"] = "true";
    
    config.parseConfig(cfg);
    
    EXPECT_TRUE(config.enable_high_fidelity_numerics);
    EXPECT_TRUE(config.enable_supg);
    EXPECT_TRUE(config.enable_advanced_time);
    EXPECT_FALSE(config.enable_p_adapt);
    EXPECT_TRUE(config.enable_physics_precond);
}


// =============================================================================
// Factory Function Tests
// =============================================================================

TEST(FactoryTest, CreateSUPGStabilization) {
    std::map<std::string, std::string> config;
    config["supg_method"] = "supg";
    config["supg_tau_formula"] = "shakib";
    
    auto supg = createSUPGStabilization(config);
    
    ASSERT_NE(supg, nullptr);
    
    // Verify it works
    PetscReal velocity[3] = {1.0, 0.0, 0.0};
    PetscReal tau = supg->calculateTau(velocity, 0.01, 0.0, 0.1, 0.01);
    EXPECT_GT(tau, 0.0);
}

TEST(FactoryTest, CreateTimeIntegration) {
    std::map<std::string, std::string> config;
    config["time_scheme"] = "bdf2";
    
    auto ts = createTimeIntegration(PETSC_COMM_WORLD, config);
    
    ASSERT_NE(ts, nullptr);
    EXPECT_EQ(ts->getOrder(), 2);
}

TEST(FactoryTest, CreatePhysicsPreconditioner) {
    std::map<std::string, std::string> config;
    config["precond_type"] = "gamg";
    
    auto pc = createPhysicsPreconditioner(PETSC_COMM_WORLD, config);
    
    ASSERT_NE(pc, nullptr);
}
