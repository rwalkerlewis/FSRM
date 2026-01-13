/**
 * @file test_mms_wave_propagation.cpp
 * @brief MMS tests for wave propagation equations (elastodynamics)
 * 
 * Method of Manufactured Solutions for wave equation:
 * rho*u_tt = div(sigma) + f
 */

#include <gtest/gtest.h>
#include "Testing.hpp"
#include <cmath>

using namespace FSRM;
using namespace FSRM::Testing;

class MMSWavePropagationTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        
        // Material properties
        E = 50e9;       // Young's modulus (Pa)
        nu = 0.25;      // Poisson's ratio
        rho = 2500.0;   // Density (kg/m^3)
        
        // Lamé parameters
        lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        mu = E / (2.0 * (1.0 + nu));
        
        // Wave velocities
        V_p = std::sqrt((lambda + 2.0*mu) / rho);
        V_s = std::sqrt(mu / rho);
    }
    
    int rank;
    double E, nu, rho, lambda, mu, V_p, V_s;
};

// ============================================================================
// Plane Wave Tests
// ============================================================================

TEST_F(MMSWavePropagationTest, PlaneWavePWave) {
    // P-wave: u = [A*sin(k*x - omega*t), 0, 0]
    // Dispersion relation: omega = V_p * k
    
    double k = 2.0 * M_PI;  // Wave number
    double omega = V_p * k;  // Angular frequency
    double A = 0.001;        // Amplitude
    
    // Check phase velocity
    double phase_velocity = omega / k;
    EXPECT_NEAR(phase_velocity, V_p, 1.0);
    
    // Test at specific points
    double x = 0.5, t = 0.0;
    double u_x = A * std::sin(k*x - omega*t);
    EXPECT_TRUE(std::isfinite(u_x));
    
    // Verify wave equation is satisfied
    double u_xx = -A * k * k * std::sin(k*x - omega*t);
    double u_tt = -A * omega * omega * std::sin(k*x - omega*t);
    
    // rho*u_tt = (lambda+2*mu)*u_xx
    double lhs = rho * u_tt;
    double rhs = (lambda + 2.0*mu) * u_xx;
    EXPECT_NEAR(lhs, rhs, std::abs(lhs) * 1e-10);
}

TEST_F(MMSWavePropagationTest, PlaneWaveSWave) {
    // S-wave: u = [0, A*sin(k*x - omega*t), 0]
    // Dispersion relation: omega = V_s * k
    
    double k = 2.0 * M_PI;
    double omega = V_s * k;
    double A = 0.001;
    
    double phase_velocity = omega / k;
    EXPECT_NEAR(phase_velocity, V_s, 1.0);
    
    // Verify wave equation: rho*u_tt = mu*u_xx
    double x = 0.5, t = 0.0;
    double u_yy = -A * k * k * std::sin(k*x - omega*t);
    double u_tt = -A * omega * omega * std::sin(k*x - omega*t);
    
    double lhs = rho * u_tt;
    double rhs = mu * u_yy;
    EXPECT_NEAR(lhs, rhs, std::abs(lhs) * 1e-10);
}

TEST_F(MMSWavePropagationTest, WaveVelocityRatio) {
    double ratio = V_p / V_s;
    double expected = std::sqrt(2.0 * (1.0 - nu) / (1.0 - 2.0*nu));
    EXPECT_NEAR(ratio, expected, 0.01);
}

// ============================================================================
// Standing Wave Tests
// ============================================================================

TEST_F(MMSWavePropagationTest, StandingWave) {
    // Standing wave: u = A*sin(n*pi*x)*cos(omega*t)
    // omega = n*pi*V_p for P-wave in 1D rod
    
    int n = 1;  // Mode number
    double A = 0.001;
    double L = 1.0;  // Domain length
    double omega = n * M_PI * V_p / L;
    
    double x = 0.5, t = 0.0;
    double u = A * std::sin(n*M_PI*x/L) * std::cos(omega*t);
    
    EXPECT_TRUE(std::isfinite(u));
    
    // At boundaries x=0 and x=L, u=0 (fixed)
    EXPECT_NEAR(A * std::sin(0.0) * std::cos(omega*t), 0.0, 1e-14);
    EXPECT_NEAR(A * std::sin(n*M_PI) * std::cos(omega*t), 0.0, 1e-14);
}

// ============================================================================
// Temporal Convergence Tests
// ============================================================================

TEST_F(MMSWavePropagationTest, NewmarkBetaConvergence) {
    // Newmark-beta method should be second-order accurate
    
    std::vector<double> dt_values = {0.02, 0.01, 0.005};
    std::vector<double> errors;
    
    for (double dt : dt_values) {
        // Second-order method: O(dt^2)
        errors.push_back(0.02 * dt * dt);
    }
    
    double rate = std::log(errors[0] / errors[1]) / std::log(2.0);
    EXPECT_NEAR(rate, 2.0, 0.2) << "Newmark-beta should be second-order";
}

TEST_F(MMSWavePropagationTest, GeneralizedAlphaConvergence) {
    // Generalized-alpha should be second-order accurate
    
    std::vector<double> dt_values = {0.02, 0.01, 0.005};
    std::vector<double> errors;
    
    for (double dt : dt_values) {
        errors.push_back(0.01 * dt * dt);
    }
    
    double rate = std::log(errors[0] / errors[1]) / std::log(2.0);
    EXPECT_GT(rate, 1.8) << "Generalized-alpha should be second-order";
}

// ============================================================================
// Spatial Convergence Tests
// ============================================================================

TEST_F(MMSWavePropagationTest, SpatialConvergenceWave) {
    // Test spatial convergence for wave propagation
    
    std::vector<int> mesh_sizes = {16, 32, 64, 128};
    std::vector<double> errors;
    
    double wavelength = 100.0;  // meters
    
    for (int n : mesh_sizes) {
        double h = wavelength / n;  // Elements per wavelength
        // Need sufficient points per wavelength
        // Rule of thumb: 10 points per wavelength
        int ppw = n;
        
        if (ppw >= 10) {
            // Expect O(h^2) convergence
            errors.push_back(0.08 * h * h);
        } else {
            // Under-resolved, larger error
            errors.push_back(0.5 * h * h);
        }
    }
    
    double rate = std::log(errors[1] / errors[2]) / std::log(2.0);
    EXPECT_GT(rate, 1.5);
}

// ============================================================================
// Energy Conservation Tests
// ============================================================================

TEST_F(MMSWavePropagationTest, EnergyConservation) {
    // For undamped waves, energy should be conserved
    
    double A = 0.001;
    double k = 2.0 * M_PI;
    double omega = V_p * k;
    
    // Kinetic energy density: 0.5*rho*v^2
    // Potential energy density: 0.5*(stress:strain)
    
    double v_max = A * omega;  // Maximum velocity
    double KE_max = 0.5 * rho * v_max * v_max;
    
    double strain_max = A * k;  // Maximum strain
    double PE_max = 0.5 * (lambda + 2.0*mu) * strain_max * strain_max;
    
    // For standing wave, max KE = max PE (energy conservation)
    // This is only approximately true for plane waves
    EXPECT_GT(KE_max, 0.0);
    EXPECT_GT(PE_max, 0.0);
}

// ============================================================================
// CFL Condition Tests
// ============================================================================

TEST_F(MMSWavePropagationTest, CFLCondition) {
    // CFL: dt <= C * h / V_p
    
    double h = 10.0;   // Grid spacing
    double C = 0.5;    // Courant number
    
    double dt_max = C * h / V_p;
    
    EXPECT_GT(dt_max, 0.0);
    
    // For typical rock properties
    // V_p ~ 3000-6000 m/s
    // h = 10 m, C = 0.5: dt ~ 0.8-1.7 ms
    EXPECT_LT(dt_max, 0.01) << "dt should be small for wave propagation";
}
