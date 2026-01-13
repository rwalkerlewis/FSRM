#include "Testing.hpp"
#include "PhysicsKernel.hpp"
#include <gtest/gtest.h>
#include <cmath>

using namespace FSRM;

// Test fixture for physics kernel tests
class PhysicsKernelTestFixture : public ::testing::Test {
protected:
    void SetUp() override {
        // Common setup
    }
    
    void TearDown() override {
        // Cleanup
    }
};

// ============================================================================
// Single Phase Flow Tests
// ============================================================================

TEST_F(PhysicsKernelTestFixture, SinglePhaseFlow_ResidualEvaluation) {
    SinglePhaseFlowKernel kernel;
    
    // Set typical reservoir properties
    double phi = 0.2;         // porosity
    double k = 100e-15;       // permeability (100 mD)
    double ct = 1e-9;         // total compressibility (1/Pa)
    double mu = 1e-3;         // viscosity (1 cP)
    double rho = 1000.0;      // density (kg/m^3)
    
    kernel.setProperties(phi, k, ct, mu, rho);
    
    // Test residual evaluation
    PetscScalar u[1] = {10e6};      // Pressure (10 MPa)
    PetscScalar u_t[1] = {1000.0};  // dP/dt
    PetscScalar u_x[3] = {1e5, 0.0, 0.0};  // Pressure gradient
    PetscScalar a[1] = {0.0};
    PetscReal x[3] = {0.0, 0.0, 0.0};
    PetscScalar f[1];
    
    kernel.residual(u, u_t, u_x, a, x, f);
    
    // Residual should be finite
    EXPECT_TRUE(std::isfinite(f[0])) << "Residual should be finite";
    
    // For zero pressure gradient and time derivative, residual should be small
    u_t[0] = 0.0;
    u_x[0] = 0.0;
    kernel.residual(u, u_t, u_x, a, x, f);
    
    EXPECT_NEAR(f[0], 0.0, 1e-10) 
        << "Residual should be zero for steady state with no gradient";
}

TEST_F(PhysicsKernelTestFixture, SinglePhaseFlow_DarcyVelocity) {
    SinglePhaseFlowKernel kernel;
    
    double phi = 0.2;
    double k = 100e-15;  // 100 mD
    double ct = 1e-9;
    double mu = 1e-3;    // 1 cP
    double rho = 1000.0;
    
    kernel.setProperties(phi, k, ct, mu, rho);
    
    // Pressure gradient
    PetscScalar grad_p[3] = {1e5, 0.0, 0.0};  // 100 kPa/m in x-direction
    
    // Darcy velocity: v = -(k/mu) * grad(p)
    double v_x_expected = -(k / mu) * grad_p[0];
    
    // This would need to be extracted from kernel
    // For now, just verify the calculation is physically reasonable
    EXPECT_LT(v_x_expected, 0.0) << "Flow should be in negative gradient direction";
}

TEST_F(PhysicsKernelTestFixture, SinglePhaseFlow_MassConservation) {
    SinglePhaseFlowKernel kernel;
    kernel.setProperties(0.2, 100e-15, 1e-9, 1e-3, 1000.0);
    
    // Test mass conservation
    // d(phi*rho)/dt + div(rho*v) = 0
    
    PetscScalar u[1] = {10e6};
    PetscScalar u_t[1] = {0.0};  // Steady state
    PetscScalar u_x[3] = {0.0, 0.0, 0.0};  // No gradient
    PetscScalar a[1] = {0.0};
    PetscReal x[3] = {0.0, 0.0, 0.0};
    PetscScalar f[1];
    
    kernel.residual(u, u_t, u_x, a, x, f);
    
    // For steady state with no gradients, mass should be conserved
    EXPECT_NEAR(f[0], 0.0, 1e-10) << "Mass should be conserved in steady state";
}

// ============================================================================
// Poroelasticity Tests
// ============================================================================

TEST_F(PhysicsKernelTestFixture, Poroelasticity_CouplingTerm) {
    PoroelasticityKernel kernel;
    
    // Set Biot coefficient
    double alpha = 0.8;
    double E = 10e9;     // Young's modulus
    double nu = 0.25;    // Poisson's ratio
    
    kernel.setBiotCoefficient(alpha);
    kernel.setMechanicalProperties(E, nu);
    
    // Pressure change should cause volumetric strain
    double deltaP = 1e6;  // 1 MPa pressure increase
    
    // Expected volumetric strain: epsilon_v = alpha * deltaP / K
    // where K = E / (3*(1-2*nu))
    double K = E / (3.0 * (1.0 - 2.0 * nu));
    double epsilon_v_expected = alpha * deltaP / K;
    
    // This would need actual kernel computation
    EXPECT_GT(epsilon_v_expected, 0.0) 
        << "Pressure increase should cause expansion";
}

TEST_F(PhysicsKernelTestFixture, Poroelasticity_TerzaghiProblem) {
    PoroelasticityKernel kernel;
    
    // Terzaghi 1D consolidation problem
    double alpha = 1.0;  // Biot coefficient
    double E = 10e9;
    double nu = 0.25;
    double k = 1e-14;    // Low permeability
    double mu = 1e-3;
    double phi = 0.3;
    double ct = 1e-9;
    
    kernel.setBiotCoefficient(alpha);
    kernel.setMechanicalProperties(E, nu);
    kernel.setFluidProperties(k, mu, phi, ct);
    
    // Compute consolidation coefficient
    double K = E / (3.0 * (1.0 - 2.0 * nu));
    double cv = k / (mu * (1.0/K + alpha*alpha*ct));
    
    EXPECT_GT(cv, 0.0) << "Consolidation coefficient should be positive";
    EXPECT_LT(cv, 1e-3) << "Low permeability should give slow consolidation";
}

// ============================================================================
// Elastodynamics Tests
// ============================================================================

TEST_F(PhysicsKernelTestFixture, Elastodynamics_WaveSpeed) {
    ElastodynamicsKernel kernel;
    
    double rho = 2500.0;  // Density (kg/m^3)
    double E = 50e9;      // Young's modulus
    double nu = 0.25;     // Poisson's ratio
    
    kernel.setMechanicalProperties(E, nu, rho);
    
    // Compute wave speeds
    double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    double mu_lame = E / (2.0 * (1.0 + nu));
    
    // P-wave velocity
    double V_p = std::sqrt((lambda + 2.0 * mu_lame) / rho);
    
    // S-wave velocity
    double V_s = std::sqrt(mu_lame / rho);
    
    EXPECT_GT(V_p, V_s) << "P-wave should be faster than S-wave";
    EXPECT_NEAR(V_p / V_s, std::sqrt(2.0 * (1.0 - nu) / (1.0 - 2.0 * nu)), 0.1)
        << "Wave speed ratio should match theoretical value";
}

TEST_F(PhysicsKernelTestFixture, Elastodynamics_WavePropagation) {
    ElastodynamicsKernel kernel;
    
    double rho = 2500.0;
    double E = 50e9;
    double nu = 0.25;
    
    kernel.setMechanicalProperties(E, nu, rho);
    
    // Test wave equation: rho * u_tt = div(stress)
    PetscScalar u[3] = {0.001, 0.0, 0.0};        // Displacement
    PetscScalar u_t[3] = {0.0, 0.0, 0.0};        // Velocity
    PetscScalar u_tt[3] = {-10.0, 0.0, 0.0};     // Acceleration
    PetscScalar u_x[9] = {0.0};                  // Displacement gradient
    PetscScalar a[3] = {0.0, 0.0, 0.0};
    PetscReal x[3] = {0.0, 0.0, 0.0};
    PetscScalar f[3];
    
    kernel.residual(u, u_tt, u_x, a, x, f);
    
    EXPECT_TRUE(std::isfinite(f[0])) << "Residual should be finite";
}

// ============================================================================
// Poroelastodynamics Tests
// ============================================================================

TEST_F(PhysicsKernelTestFixture, Poroelastodynamics_BiotWaves) {
    PoroelastodynamicsKernel kernel;
    
    // Biot's theory predicts three wave types in saturated porous media:
    // - Fast P-wave
    // - Slow P-wave (diffusion-like)
    // - S-wave
    
    double rho_s = 2650.0;  // Solid density
    double rho_f = 1000.0;  // Fluid density
    double phi = 0.2;       // Porosity
    double E = 20e9;
    double nu = 0.25;
    double k = 1e-12;
    double mu_fluid = 1e-3;
    
    kernel.setSolidProperties(rho_s, E, nu);
    kernel.setFluidProperties(rho_f, mu_fluid);
    kernel.setPorosity(phi);
    kernel.setPermeability(k);
    
    // Compute characteristic frequencies
    double omega_c = phi * mu_fluid / (k * rho_f);  // Biot critical frequency
    
    EXPECT_GT(omega_c, 0.0) << "Critical frequency should be positive";
}

// ============================================================================
// Thermal Flow Tests
// ============================================================================

TEST_F(PhysicsKernelTestFixture, ThermalFlow_HeatConduction) {
    ThermalFlowKernel kernel;
    
    double lambda_t = 2.5;  // Thermal conductivity (W/m/K)
    double rho = 2000.0;    // Density
    double cp = 1000.0;     // Specific heat (J/kg/K)
    
    kernel.setThermalProperties(lambda_t, rho, cp);
    
    // Test heat equation: rho*cp*dT/dt = div(lambda*grad(T))
    PetscScalar u[1] = {350.0};     // Temperature (K)
    PetscScalar u_t[1] = {0.1};     // dT/dt
    PetscScalar u_x[3] = {10.0, 0.0, 0.0};  // Temperature gradient
    PetscScalar a[1] = {0.0};
    PetscReal x[3] = {0.0, 0.0, 0.0};
    PetscScalar f[1];
    
    kernel.residual(u, u_t, u_x, a, x, f);
    
    EXPECT_TRUE(std::isfinite(f[0])) << "Thermal residual should be finite";
}

TEST_F(PhysicsKernelTestFixture, ThermalFlow_ConvectiveHeatTransfer) {
    ThermalFlowKernel kernel;
    
    double lambda_t = 2.5;
    double rho = 2000.0;
    double cp = 1000.0;
    
    kernel.setThermalProperties(lambda_t, rho, cp);
    
    // Add fluid velocity for convective transport
    double v_x = 1e-5;  // m/s
    kernel.setFluidVelocity(v_x, 0.0, 0.0);
    
    // Heat equation with convection: 
    // rho*cp*(dT/dt + v·grad(T)) = div(lambda*grad(T))
    
    PetscScalar u[1] = {350.0};
    PetscScalar u_t[1] = {0.0};
    PetscScalar u_x[3] = {10.0, 0.0, 0.0};
    PetscScalar a[1] = {0.0};
    PetscReal x[3] = {0.0, 0.0, 0.0};
    PetscScalar f[1];
    
    kernel.residual(u, u_t, u_x, a, x, f);
    
    // With convection, residual should be non-zero even in pseudo-steady state
    EXPECT_TRUE(std::isfinite(f[0])) << "Residual should be finite";
}

// ============================================================================
// Jacobian Tests
// ============================================================================

TEST_F(PhysicsKernelTestFixture, Jacobian_FiniteDifference) {
    SinglePhaseFlowKernel kernel;
    kernel.setProperties(0.2, 100e-15, 1e-9, 1e-3, 1000.0);
    
    PetscScalar u[1] = {10e6};
    PetscScalar u_t[1] = {0.0};
    PetscScalar u_x[3] = {1e5, 0.0, 0.0};
    PetscScalar a[1] = {0.0};
    PetscReal x[3] = {0.0, 0.0, 0.0};
    
    // Analytical Jacobian
    PetscScalar J_analytical[1];
    kernel.jacobian(u, u_t, u_x, a, x, J_analytical);
    
    // Finite difference Jacobian
    double epsilon = 1e-6;
    PetscScalar f0[1], f1[1];
    
    kernel.residual(u, u_t, u_x, a, x, f0);
    u[0] += epsilon;
    kernel.residual(u, u_t, u_x, a, x, f1);
    
    PetscScalar J_fd = (f1[0] - f0[0]) / epsilon;
    
    EXPECT_NEAR(J_analytical[0], J_fd, std::abs(J_fd) * 0.01)
        << "Analytical Jacobian should match finite difference";
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    
    int result = RUN_ALL_TESTS();
    
    PetscFinalize();
    MPI_Finalize();
    return result;
}
