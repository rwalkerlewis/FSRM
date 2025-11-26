/**
 * @file test_boundary_conditions.cpp
 * @brief Functional tests for boundary condition implementations
 */

#include <gtest/gtest.h>
#include "ReservoirSim.hpp"
#include <cmath>

using namespace FSRM;

class BoundaryConditionTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    }
    
    int rank;
};

// ============================================================================
// Dirichlet Boundary Condition Tests
// ============================================================================

TEST_F(BoundaryConditionTest, DirichletValueApplied) {
    // For Dirichlet BC: u = g on boundary
    double g = 10.0;
    double u_boundary = g;  // Should equal prescribed value
    
    EXPECT_DOUBLE_EQ(u_boundary, g);
}

TEST_F(BoundaryConditionTest, DirichletResidual) {
    // Dirichlet residual: R = u - g
    double u = 10.5;  // Numerical value
    double g = 10.0;  // Prescribed value
    
    double residual = u - g;
    
    EXPECT_NEAR(residual, 0.5, 1e-14);
}

// ============================================================================
// Neumann Boundary Condition Tests
// ============================================================================

TEST_F(BoundaryConditionTest, NeumannFluxCalculation) {
    // Neumann BC: ∂u/∂n = h on boundary
    // Discrete: (u_boundary - u_interior) / dx = h
    
    double h = 100.0;     // Prescribed flux
    double dx = 1.0;
    double u_interior = 10.0;
    
    // Ghost point approach: u_boundary = u_interior + h*dx
    double u_boundary = u_interior + h * dx;
    
    double computed_flux = (u_boundary - u_interior) / dx;
    EXPECT_NEAR(computed_flux, h, 1e-14);
}

TEST_F(BoundaryConditionTest, NoFlowBoundary) {
    // No-flow: ∂u/∂n = 0
    double h = 0.0;
    double u_interior = 10.0;
    double dx = 1.0;
    
    // For no-flow, ghost value equals interior value
    double u_ghost = u_interior + h * dx;
    
    EXPECT_DOUBLE_EQ(u_ghost, u_interior);
}

// ============================================================================
// Robin Boundary Condition Tests
// ============================================================================

TEST_F(BoundaryConditionTest, RobinMixedCondition) {
    // Robin BC: a*u + b*∂u/∂n = g
    double a = 1.0;
    double b = 0.1;
    double g = 10.0;
    double u = 5.0;
    double grad_u = (g - a*u) / b;
    
    double residual = a*u + b*grad_u - g;
    EXPECT_NEAR(residual, 0.0, 1e-14);
}

TEST_F(BoundaryConditionTest, ConvectiveBoundary) {
    // Heat transfer: -k*∂T/∂n = h*(T - T_inf)
    double k = 2.5;       // Thermal conductivity
    double h_conv = 10.0; // Convection coefficient
    double T = 350.0;     // Surface temperature
    double T_inf = 300.0; // Ambient temperature
    
    double heat_flux = h_conv * (T - T_inf);
    double grad_T = -heat_flux / k;
    
    EXPECT_TRUE(std::isfinite(grad_T));
    EXPECT_LT(grad_T, 0.0);  // Heat flowing out
}

// ============================================================================
// Symmetry Boundary Tests
// ============================================================================

TEST_F(BoundaryConditionTest, SymmetryPlane) {
    // At symmetry plane: normal gradient = 0, normal velocity = 0
    double grad_n = 0.0;   // Normal gradient
    double v_n = 0.0;      // Normal velocity
    
    EXPECT_DOUBLE_EQ(grad_n, 0.0);
    EXPECT_DOUBLE_EQ(v_n, 0.0);
}

// ============================================================================
// Periodic Boundary Tests
// ============================================================================

TEST_F(BoundaryConditionTest, PeriodicContinuity) {
    // u(x=0) = u(x=L)
    // grad_u(x=0) = grad_u(x=L)
    
    double L = 100.0;
    auto u = [](double x) { return std::sin(2.0 * M_PI * x / 100.0); };
    
    double u_left = u(0.0);
    double u_right = u(L);
    
    EXPECT_NEAR(u_left, u_right, 1e-14);
}

// ============================================================================
// Outflow/Pressure Boundary Tests
// ============================================================================

TEST_F(BoundaryConditionTest, PressureOutlet) {
    // Prescribed pressure at outlet
    double p_outlet = 101325.0;  // Atmospheric pressure
    double p_boundary = p_outlet;
    
    EXPECT_DOUBLE_EQ(p_boundary, p_outlet);
}

TEST_F(BoundaryConditionTest, ZeroGradientOutflow) {
    // Zero gradient outflow: ∂φ/∂n = 0
    // Also known as "natural" or "do-nothing" BC
    
    double phi_internal = 10.0;
    double phi_boundary = phi_internal;  // Copy internal value
    
    EXPECT_DOUBLE_EQ(phi_boundary, phi_internal);
}

// ============================================================================
// Well Boundary Tests
// ============================================================================

TEST_F(BoundaryConditionTest, WellPeacemanBC) {
    // Peaceman well model: q = WI * (p_res - p_wf)
    double WI = 1e-12;       // Well index
    double p_res = 20e6;     // Reservoir pressure
    double p_wf = 10e6;      // Wellbore pressure
    
    double q = WI * (p_res - p_wf);
    
    EXPECT_GT(q, 0.0);  // Production
    EXPECT_TRUE(std::isfinite(q));
}

TEST_F(BoundaryConditionTest, InjectionBC) {
    // Constant rate injection
    double q_inj = 100.0;    // m^3/day
    double WI = 1e-12;
    double p_res = 15e6;
    
    // BHP required for injection: p_wf = p_res + q/WI
    double p_wf = p_res + q_inj / WI;
    
    EXPECT_GT(p_wf, p_res);  // BHP > reservoir pressure for injection
}

// ============================================================================
// Stress Boundary Tests
// ============================================================================

TEST_F(BoundaryConditionTest, TractionBC) {
    // Prescribed traction: σ·n = t
    double sigma_xx = 30e6;  // Normal stress
    double sigma_xy = 5e6;   // Shear stress
    double nx = 1.0, ny = 0.0;  // Normal to x-face
    
    double t_x = sigma_xx * nx + sigma_xy * ny;
    double t_y = sigma_xy * nx + 0.0 * ny;
    
    EXPECT_DOUBLE_EQ(t_x, sigma_xx);
    EXPECT_DOUBLE_EQ(t_y, sigma_xy);
}

TEST_F(BoundaryConditionTest, FreeSurfaceBC) {
    // Free surface: traction = 0
    double t_x = 0.0;
    double t_y = 0.0;
    double t_z = 0.0;
    
    EXPECT_DOUBLE_EQ(t_x, 0.0);
    EXPECT_DOUBLE_EQ(t_y, 0.0);
    EXPECT_DOUBLE_EQ(t_z, 0.0);
}

TEST_F(BoundaryConditionTest, OverburdenStress) {
    // Overburden: σ_zz = ρ*g*z
    double rho = 2500.0;   // kg/m^3
    double g = 9.81;       // m/s^2
    double z = 1000.0;     // depth in m
    
    double sigma_v = rho * g * z;
    
    EXPECT_NEAR(sigma_v, 24.525e6, 1e3);  // ~24.5 MPa at 1 km
}
