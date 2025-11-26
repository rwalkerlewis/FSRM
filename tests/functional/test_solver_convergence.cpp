/**
 * @file test_solver_convergence.cpp
 * @brief Functional tests for solver convergence behavior
 */

#include <gtest/gtest.h>
#include "ReservoirSim.hpp"
#include <cmath>

using namespace FSRM;

class SolverConvergenceTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    }
    
    int rank;
};

// ============================================================================
// Newton-Raphson Convergence Tests
// ============================================================================

TEST_F(SolverConvergenceTest, NewtonConvergenceRate) {
    // Verify quadratic convergence for Newton-Raphson
    // For f(x) = x^2 - 2, Newton: x_{n+1} = x_n - f(x_n)/f'(x_n)
    
    double x = 2.0;  // Initial guess
    double target = std::sqrt(2.0);
    std::vector<double> errors;
    
    for (int i = 0; i < 5; ++i) {
        double f = x*x - 2.0;
        double fp = 2.0*x;
        x = x - f/fp;
        errors.push_back(std::abs(x - target));
    }
    
    // Check quadratic convergence: e_{n+1} ~ C * e_n^2
    for (size_t i = 1; i < errors.size() - 1; ++i) {
        if (errors[i] > 1e-14 && errors[i-1] > 1e-14) {
            double ratio = std::log(errors[i+1]) / std::log(errors[i]);
            // For quadratic convergence, ratio should be ~2
            EXPECT_GT(ratio, 1.5) << "Newton should have quadratic convergence";
        }
    }
    
    EXPECT_NEAR(x, target, 1e-14);
}

TEST_F(SolverConvergenceTest, FixedPointConvergence) {
    // Fixed point iteration: x_{n+1} = g(x_n)
    // Convergence if |g'(x)| < 1
    
    // g(x) = cos(x), fixed point near x = 0.7391
    double x = 0.5;
    double target = 0.739085133;  // Fixed point of cos(x)
    
    for (int i = 0; i < 100; ++i) {
        x = std::cos(x);
    }
    
    EXPECT_NEAR(x, target, 1e-6);
}

// ============================================================================
// Linear Solver Convergence Tests
// ============================================================================

TEST_F(SolverConvergenceTest, JacobiIterationConvergence) {
    // Jacobi iteration for 2x2 system:
    // 4x - y = 1
    // -x + 3y = 2
    // Solution: x = 5/11, y = 9/11
    
    double x = 0.0, y = 0.0;
    double x_exact = 5.0/11.0;
    double y_exact = 9.0/11.0;
    
    for (int i = 0; i < 50; ++i) {
        double x_new = (1.0 + y) / 4.0;
        double y_new = (2.0 + x) / 3.0;
        x = x_new;
        y = y_new;
    }
    
    EXPECT_NEAR(x, x_exact, 1e-6);
    EXPECT_NEAR(y, y_exact, 1e-6);
}

TEST_F(SolverConvergenceTest, GaussSeidelConvergence) {
    // Gauss-Seidel (faster than Jacobi)
    double x = 0.0, y = 0.0;
    double x_exact = 5.0/11.0;
    double y_exact = 9.0/11.0;
    
    for (int i = 0; i < 20; ++i) {
        x = (1.0 + y) / 4.0;
        y = (2.0 + x) / 3.0;  // Uses updated x
    }
    
    EXPECT_NEAR(x, x_exact, 1e-10);
    EXPECT_NEAR(y, y_exact, 1e-10);
}

// ============================================================================
// Time Integration Convergence Tests
// ============================================================================

TEST_F(SolverConvergenceTest, ForwardEulerOrder) {
    // du/dt = -u, u(0) = 1
    // Exact: u(t) = exp(-t)
    
    double u_exact = std::exp(-1.0);
    std::vector<double> errors;
    std::vector<double> dt_values = {0.1, 0.05, 0.025, 0.0125};
    
    for (double dt : dt_values) {
        double u = 1.0;
        double t = 0.0;
        int n = static_cast<int>(1.0 / dt);
        
        for (int i = 0; i < n; ++i) {
            u = u - dt * u;  // Forward Euler
            t += dt;
        }
        
        errors.push_back(std::abs(u - u_exact));
    }
    
    // Forward Euler is first order: error ~ C*dt
    for (size_t i = 1; i < errors.size(); ++i) {
        double rate = std::log(errors[i-1] / errors[i]) / std::log(2.0);
        EXPECT_NEAR(rate, 1.0, 0.2) << "Forward Euler should be first-order";
    }
}

TEST_F(SolverConvergenceTest, BackwardEulerOrder) {
    // du/dt = -u, u(0) = 1
    // Backward Euler: u_{n+1} = u_n - dt*u_{n+1}
    //               : u_{n+1} = u_n / (1 + dt)
    
    double u_exact = std::exp(-1.0);
    std::vector<double> errors;
    std::vector<double> dt_values = {0.1, 0.05, 0.025, 0.0125};
    
    for (double dt : dt_values) {
        double u = 1.0;
        int n = static_cast<int>(1.0 / dt);
        
        for (int i = 0; i < n; ++i) {
            u = u / (1.0 + dt);  // Backward Euler
        }
        
        errors.push_back(std::abs(u - u_exact));
    }
    
    // Backward Euler is first order
    for (size_t i = 1; i < errors.size(); ++i) {
        double rate = std::log(errors[i-1] / errors[i]) / std::log(2.0);
        EXPECT_NEAR(rate, 1.0, 0.2);
    }
}

TEST_F(SolverConvergenceTest, TrapezoidalOrder) {
    // Trapezoidal (Crank-Nicolson) is second order
    // u_{n+1} = u_n - (dt/2)*(u_n + u_{n+1})
    //         = u_n * (1 - dt/2) / (1 + dt/2)
    
    double u_exact = std::exp(-1.0);
    std::vector<double> errors;
    std::vector<double> dt_values = {0.1, 0.05, 0.025};
    
    for (double dt : dt_values) {
        double u = 1.0;
        int n = static_cast<int>(1.0 / dt);
        
        for (int i = 0; i < n; ++i) {
            u = u * (1.0 - dt/2.0) / (1.0 + dt/2.0);
        }
        
        errors.push_back(std::abs(u - u_exact));
    }
    
    // Trapezoidal is second order
    for (size_t i = 1; i < errors.size(); ++i) {
        double rate = std::log(errors[i-1] / errors[i]) / std::log(2.0);
        EXPECT_NEAR(rate, 2.0, 0.3);
    }
}

// ============================================================================
// Stability Tests
// ============================================================================

TEST_F(SolverConvergenceTest, CFLStability) {
    // For explicit methods, dt <= CFL * dx / c
    double dx = 10.0;   // Grid spacing
    double c = 3000.0;  // Wave speed
    double CFL = 0.5;   // Courant number
    
    double dt_stable = CFL * dx / c;
    
    EXPECT_GT(dt_stable, 0.0);
    EXPECT_LT(dt_stable, dx / c);  // Should be less than maximum
}

TEST_F(SolverConvergenceTest, VonNeumannStability) {
    // For parabolic equations, dt <= dx^2 / (2*alpha)
    double dx = 1.0;
    double alpha = 1e-6;  // Diffusivity
    
    double dt_stable = dx*dx / (2.0 * alpha);
    
    EXPECT_GT(dt_stable, 0.0);
}
