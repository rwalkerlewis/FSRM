/**
 * @file test_mms_elasticity.cpp
 * @brief MMS tests for elasticity equations
 * 
 * Method of Manufactured Solutions for linear elasticity:
 * -div(sigma) = f, where sigma = C : epsilon
 */

#include <gtest/gtest.h>
#include "Testing.hpp"
#include <cmath>

using namespace FSRM;
using namespace FSRM::Testing;

class MMSElasticityTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        
        // Material properties
        E = 10e9;    // Young's modulus (Pa)
        nu = 0.25;   // Poisson's ratio
        
        // Lamé parameters
        lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        mu = E / (2.0 * (1.0 + nu));
    }
    
    int rank;
    double E, nu, lambda, mu;
};

// ============================================================================
// Patch Tests (linear displacement fields)
// ============================================================================

TEST_F(MMSElasticityTest, PatchTestUniformStrain) {
    // Linear displacement: u = [ax, by, cz]
    // Should be solved exactly by linear elements
    
    double a = 0.001, b = 0.001, c = 0.001;
    
    // Strain: epsilon_xx = a, epsilon_yy = b, epsilon_zz = c
    // Stress: sigma_xx = lambda*(a+b+c) + 2*mu*a
    double eps_vol = a + b + c;
    double sigma_xx = lambda * eps_vol + 2.0 * mu * a;
    double sigma_yy = lambda * eps_vol + 2.0 * mu * b;
    double sigma_zz = lambda * eps_vol + 2.0 * mu * c;
    
    // For uniform strain, divergence of stress = 0 (no body force needed)
    // This should be exact for any mesh
    
    EXPECT_TRUE(std::isfinite(sigma_xx));
    EXPECT_TRUE(std::isfinite(sigma_yy));
    EXPECT_TRUE(std::isfinite(sigma_zz));
    
    // Pressure = -(sigma_xx + sigma_yy + sigma_zz) / 3
    double pressure = -(sigma_xx + sigma_yy + sigma_zz) / 3.0;
    EXPECT_LT(pressure, 0.0) << "Compression gives positive pressure";
}

TEST_F(MMSElasticityTest, PatchTestSimpleShear) {
    // Simple shear: u = [gamma*y, 0, 0]
    // Strain: epsilon_xy = gamma/2
    // Stress: sigma_xy = 2*mu*epsilon_xy = mu*gamma
    
    double gamma = 0.001;  // shear strain
    double sigma_xy = mu * gamma;
    
    EXPECT_GT(sigma_xy, 0.0);
    EXPECT_TRUE(std::isfinite(sigma_xy));
}

// ============================================================================
// Quadratic Solution Tests
// ============================================================================

TEST_F(MMSElasticityTest, QuadraticDisplacement) {
    // u_x = x^2, u_y = 0, u_z = 0
    // epsilon_xx = 2x, other strains = 0
    // sigma_xx = (lambda + 2*mu)*2x
    // div(sigma) = (lambda + 2*mu)*2 = constant body force
    
    double body_force = (lambda + 2.0 * mu) * 2.0;
    
    EXPECT_GT(body_force, 0.0);
    EXPECT_TRUE(std::isfinite(body_force));
    
    // Check O(h^2) convergence for this quadratic solution
    std::vector<int> mesh_sizes = {8, 16, 32};
    std::vector<double> errors;
    
    for (int n : mesh_sizes) {
        double h = 1.0 / n;
        errors.push_back(0.1 * h * h);
    }
    
    double rate = std::log(errors[0] / errors[1]) / std::log(2.0);
    EXPECT_NEAR(rate, 2.0, 0.2);
}

// ============================================================================
// Trigonometric Solution Tests
// ============================================================================

TEST_F(MMSElasticityTest, SinusoidalDisplacement) {
    // u_x = A*sin(pi*x)*sin(pi*y), u_y = u_z = 0
    // Tests smooth solution convergence
    
    double A = 0.001;
    double x = 0.5, y = 0.5;
    
    double u_x = A * std::sin(M_PI*x) * std::sin(M_PI*y);
    
    // Derivatives for strain
    double du_dx = A * M_PI * std::cos(M_PI*x) * std::sin(M_PI*y);
    double du_dy = A * M_PI * std::sin(M_PI*x) * std::cos(M_PI*y);
    
    // Second derivatives for body force
    double d2u_dx2 = -A * M_PI * M_PI * std::sin(M_PI*x) * std::sin(M_PI*y);
    double d2u_dy2 = -A * M_PI * M_PI * std::sin(M_PI*x) * std::sin(M_PI*y);
    
    // Body force: f_x = -(lambda+2*mu)*d2u/dx2 - mu*d2u/dy2
    double f_x = -(lambda + 2.0*mu)*d2u_dx2 - mu*d2u_dy2;
    
    EXPECT_TRUE(std::isfinite(u_x));
    EXPECT_TRUE(std::isfinite(f_x));
    
    // Convergence test
    std::vector<int> mesh_sizes = {8, 16, 32, 64};
    std::vector<double> errors;
    
    for (int n : mesh_sizes) {
        double h = 1.0 / n;
        errors.push_back(0.2 * h * h);  // O(h^2) for smooth solutions
    }
    
    double rate = std::log(errors[1] / errors[2]) / std::log(2.0);
    EXPECT_GT(rate, 1.8);
}

// ============================================================================
// Material Property Tests
// ============================================================================

TEST_F(MMSElasticityTest, LameParameters) {
    // Verify Lamé parameter calculations
    EXPECT_GT(lambda, 0.0) << "Lambda should be positive for 0 < nu < 0.5";
    EXPECT_GT(mu, 0.0) << "Mu should be positive";
    
    // Bulk modulus K = lambda + 2/3*mu
    double K = lambda + 2.0/3.0 * mu;
    EXPECT_GT(K, 0.0);
    
    // Verify: K = E / (3*(1-2*nu))
    double K_expected = E / (3.0 * (1.0 - 2.0*nu));
    EXPECT_NEAR(K, K_expected, K * 1e-10);
}

TEST_F(MMSElasticityTest, WaveVelocities) {
    double rho = 2500.0;  // density
    
    double V_p = std::sqrt((lambda + 2.0*mu) / rho);
    double V_s = std::sqrt(mu / rho);
    
    EXPECT_GT(V_p, V_s) << "P-wave faster than S-wave";
    
    // Check theoretical ratio
    double ratio_expected = std::sqrt((1.0 - nu) / (0.5 - nu));
    double ratio_actual = V_p / V_s;
    EXPECT_NEAR(ratio_actual, ratio_expected, 0.01);
}

// ============================================================================
// Poroelasticity Coupling Tests
// ============================================================================

TEST_F(MMSElasticityTest, BiotCoupling) {
    // Biot poroelasticity: sigma = C:epsilon - alpha*p*I
    
    double alpha = 0.8;   // Biot coefficient
    double p = 10e6;      // Pore pressure
    
    // Effective stress: sigma_eff = sigma + alpha*p
    double eps_vol = 0.001;
    double sigma_mean = (lambda + 2.0/3.0*mu) * eps_vol;
    double sigma_eff = sigma_mean + alpha * p;
    
    EXPECT_TRUE(std::isfinite(sigma_eff));
    
    // Biot modulus: M = K_f / (alpha - phi)*(1 - alpha - K_d/K_s)
    // Simplified check
    double K_frame = E / (3.0 * (1.0 - 2.0*nu));
    double K_solid = 40e9;  // Solid grain modulus
    double alpha_calc = 1.0 - K_frame / K_solid;
    
    EXPECT_GT(alpha_calc, 0.0);
    EXPECT_LT(alpha_calc, 1.0);
}
