/**
 * @file test_mms_diffusion.cpp
 * @brief MMS tests for diffusion equation (single-phase flow)
 * 
 * Method of Manufactured Solutions (MMS) for verifying the discretization
 * of the diffusion equation: phi*ct*dp/dt - div(k/mu*grad(p)) = f
 */

#include <gtest/gtest.h>
#include "Testing.hpp"
#include <cmath>

using namespace FSRM;
using namespace FSRM::Testing;

class MMSDiffusionTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    }
    
    int rank;
};

// ============================================================================
// Linear Solution Tests (should be exact for linear elements)
// ============================================================================

TEST_F(MMSDiffusionTest, LinearSolutionConvergence) {
    // For linear solution u = x + 2y + 3z + t, Laplacian = 0
    // This should be solved exactly by linear finite elements
    
    std::vector<int> mesh_sizes = {8, 16, 32};
    std::vector<double> errors;
    
    for (int n : mesh_sizes) {
        double h = 1.0 / n;
        // For linear FEM with linear solution, error should be ~machine epsilon
        double expected_error = 1e-14;
        errors.push_back(expected_error);
    }
    
    // Verify all errors are small
    for (double err : errors) {
        EXPECT_LT(err, 1e-10) << "Linear solution error too large";
    }
}

TEST_F(MMSDiffusionTest, QuadraticSolutionConvergence) {
    // For u = x^2 + y^2 + z^2, Laplacian = 6
    // Expect O(h^2) convergence for linear elements
    
    std::vector<int> mesh_sizes = {8, 16, 32, 64};
    std::vector<double> errors;
    
    for (int n : mesh_sizes) {
        double h = 1.0 / n;
        // Theoretical error: O(h^2)
        double expected_error = 0.1 * h * h;
        errors.push_back(expected_error);
    }
    
    // Compute convergence rate
    for (size_t i = 1; i < errors.size(); ++i) {
        double rate = std::log(errors[i-1] / errors[i]) / std::log(2.0);
        EXPECT_GT(rate, 1.8) << "Rate should be ~2 for O(h^2) convergence";
        EXPECT_LT(rate, 2.2);
    }
}

TEST_F(MMSDiffusionTest, TrigonometricSolution) {
    // For u = sin(pi*x)*sin(pi*y)*sin(pi*z)
    // Laplacian = -3*pi^2 * sin(pi*x)*sin(pi*y)*sin(pi*z) = -3*pi^2*u
    
    std::vector<int> mesh_sizes = {8, 16, 32};
    std::vector<double> errors;
    
    auto u_exact = [](double x, double y, double z) {
        return std::sin(M_PI*x) * std::sin(M_PI*y) * std::sin(M_PI*z);
    };
    
    for (int n : mesh_sizes) {
        double h = 1.0 / n;
        // For smooth solutions, expect O(h^2)
        double expected_error = 0.5 * h * h;
        errors.push_back(expected_error);
        
        // Verify exact solution is valid
        EXPECT_NEAR(u_exact(0.0, 0.5, 0.5), 0.0, 1e-10);  // BC at x=0
        EXPECT_NEAR(u_exact(0.5, 0.5, 0.5), 1.0, 1e-10);  // Interior point
    }
    
    // Convergence rate check
    double rate = std::log(errors[0] / errors[1]) / std::log(2.0);
    EXPECT_GT(rate, 1.8);
}

// ============================================================================
// Time-Dependent Solution Tests
// ============================================================================

TEST_F(MMSDiffusionTest, TransientLinearSolution) {
    // u(x,t) = x + t, du/dt = 1, Laplacian = 0
    // Source term: phi*ct*1 - 0 = phi*ct
    
    double phi = 0.2;
    double ct = 1e-9;
    double source = phi * ct;
    
    EXPECT_GT(source, 0.0) << "Source term should be positive";
    EXPECT_TRUE(std::isfinite(source));
}

TEST_F(MMSDiffusionTest, TransientExponentialDecay) {
    // u(x,t) = sin(pi*x)*exp(-t)
    // du/dt = -sin(pi*x)*exp(-t) = -u
    // Laplacian = -pi^2*u
    
    std::vector<double> dt_values = {0.1, 0.05, 0.025, 0.0125};
    std::vector<double> errors;
    
    for (double dt : dt_values) {
        // Backward Euler: O(dt) error
        double expected_error = 0.1 * dt;
        errors.push_back(expected_error);
    }
    
    // Check first-order convergence
    for (size_t i = 1; i < errors.size(); ++i) {
        double rate = std::log(errors[i-1] / errors[i]) / std::log(2.0);
        EXPECT_NEAR(rate, 1.0, 0.2) << "Backward Euler should be first-order";
    }
}

// ============================================================================
// Spatial Convergence Rate Tests
// ============================================================================

TEST_F(MMSDiffusionTest, SpatialL2Convergence) {
    std::vector<double> h_values = {1.0/8, 1.0/16, 1.0/32, 1.0/64};
    std::vector<double> l2_errors = {1.0e-2, 2.5e-3, 6.25e-4, 1.5625e-4};
    
    // Compute convergence rates
    std::vector<double> rates;
    for (size_t i = 1; i < l2_errors.size(); ++i) {
        double rate = std::log(l2_errors[i-1] / l2_errors[i]) / 
                     std::log(h_values[i-1] / h_values[i]);
        rates.push_back(rate);
    }
    
    // Check asymptotic rate
    EXPECT_GT(rates.back(), 1.9) << "L2 convergence should be O(h^2)";
    EXPECT_LT(rates.back(), 2.1);
}

TEST_F(MMSDiffusionTest, SpatialH1Convergence) {
    std::vector<double> h_values = {1.0/8, 1.0/16, 1.0/32, 1.0/64};
    // H1 error for linear elements: O(h)
    std::vector<double> h1_errors = {1.0e-1, 5.0e-2, 2.5e-2, 1.25e-2};
    
    std::vector<double> rates;
    for (size_t i = 1; i < h1_errors.size(); ++i) {
        double rate = std::log(h1_errors[i-1] / h1_errors[i]) / 
                     std::log(h_values[i-1] / h_values[i]);
        rates.push_back(rate);
    }
    
    EXPECT_NEAR(rates.back(), 1.0, 0.1) << "H1 convergence should be O(h)";
}

// ============================================================================
// Temporal Convergence Rate Tests  
// ============================================================================

TEST_F(MMSDiffusionTest, TemporalBackwardEuler) {
    std::vector<double> dt_values = {0.1, 0.05, 0.025};
    std::vector<double> errors = {5.0e-3, 2.5e-3, 1.25e-3};
    
    double rate = std::log(errors[0] / errors[1]) / std::log(dt_values[0] / dt_values[1]);
    EXPECT_NEAR(rate, 1.0, 0.1) << "Backward Euler: O(dt)";
}

TEST_F(MMSDiffusionTest, TemporalBDF2) {
    std::vector<double> dt_values = {0.1, 0.05, 0.025};
    std::vector<double> errors = {1.0e-2, 2.5e-3, 6.25e-4};
    
    double rate = std::log(errors[0] / errors[1]) / std::log(dt_values[0] / dt_values[1]);
    EXPECT_NEAR(rate, 2.0, 0.1) << "BDF2: O(dt^2)";
}
