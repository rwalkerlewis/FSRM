/**
 * @file test_physics_kernels.cpp
 * @brief Unit tests for physics kernel classes
 */

#include <gtest/gtest.h>
#include "PhysicsKernel.hpp"
#include "ReservoirSim.hpp"
#include <cmath>

using namespace FSRM;

class PhysicsKernelBaseTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    }
    
    int rank;
};

// ============================================================================
// SinglePhaseFlowKernel Tests
// ============================================================================

TEST_F(PhysicsKernelBaseTest, SinglePhaseFlowCreate) {
    SinglePhaseFlowKernel kernel;
    EXPECT_EQ(kernel.getType(), PhysicsType::FLUID_FLOW);
    EXPECT_EQ(kernel.getNumFields(), 1);
    EXPECT_EQ(kernel.getNumComponents(0), 1);
}

TEST_F(PhysicsKernelBaseTest, SinglePhaseFlowSetProperties) {
    SinglePhaseFlowKernel kernel;
    
    // Should not throw
    EXPECT_NO_THROW(kernel.setProperties(0.2, 100e-15, 1e-9, 1e-3, 1000.0));
}

TEST_F(PhysicsKernelBaseTest, SinglePhaseFlowResidualFinite) {
    SinglePhaseFlowKernel kernel;
    kernel.setProperties(0.2, 100e-15, 1e-9, 1e-3, 1000.0);
    
    PetscScalar u[1] = {10e6};
    PetscScalar u_t[1] = {1000.0};
    PetscScalar u_x[3] = {1e5, 0.0, 0.0};
    PetscScalar a[1] = {0.0};
    PetscReal x[3] = {0.0, 0.0, 0.0};
    PetscScalar f[1] = {0.0};
    
    kernel.residual(u, u_t, u_x, a, x, f);
    
    EXPECT_TRUE(std::isfinite(f[0])) << "Residual should be finite";
}

TEST_F(PhysicsKernelBaseTest, SinglePhaseFlowSteadyState) {
    SinglePhaseFlowKernel kernel;
    kernel.setProperties(0.2, 100e-15, 1e-9, 1e-3, 1000.0);
    
    PetscScalar u[1] = {10e6};
    PetscScalar u_t[1] = {0.0};  // Steady state
    PetscScalar u_x[3] = {0.0, 0.0, 0.0};  // No gradient
    PetscScalar a[1] = {0.0};
    PetscReal x[3] = {0.0, 0.0, 0.0};
    PetscScalar f[1] = {999.0};
    
    kernel.residual(u, u_t, u_x, a, x, f);
    
    // For steady state with no gradient, residual should be small
    EXPECT_NEAR(f[0], 0.0, 1e-6);
}

// ============================================================================
// GeomechanicsKernel Tests
// ============================================================================

TEST_F(PhysicsKernelBaseTest, GeomechanicsCreate) {
    GeomechanicsKernel kernel(SolidModelType::ELASTIC);
    EXPECT_EQ(kernel.getType(), PhysicsType::GEOMECHANICS);
    EXPECT_EQ(kernel.getNumFields(), 1);
    EXPECT_EQ(kernel.getNumComponents(0), 3);  // 3D displacement
}

TEST_F(PhysicsKernelBaseTest, GeomechanicsSetProperties) {
    GeomechanicsKernel kernel(SolidModelType::ELASTIC);
    EXPECT_NO_THROW(kernel.setMaterialProperties(10e9, 0.25, 2500.0));
}

TEST_F(PhysicsKernelBaseTest, GeomechanicsResidualFinite) {
    GeomechanicsKernel kernel(SolidModelType::ELASTIC);
    kernel.setMaterialProperties(10e9, 0.25, 2500.0);
    
    PetscScalar u[3] = {0.001, 0.0, 0.0};
    PetscScalar u_t[3] = {0.0, 0.0, 0.0};
    PetscScalar u_x[9] = {0.001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    PetscScalar a[3] = {0.0, 0.0, 0.0};
    PetscReal x[3] = {0.0, 0.0, 0.0};
    PetscScalar f[3] = {0.0, 0.0, 0.0};
    
    kernel.residual(u, u_t, u_x, a, x, f);
    
    EXPECT_TRUE(std::isfinite(f[0]));
    EXPECT_TRUE(std::isfinite(f[1]));
    EXPECT_TRUE(std::isfinite(f[2]));
}

// ============================================================================
// ElastodynamicsKernel Tests
// ============================================================================

TEST_F(PhysicsKernelBaseTest, ElastodynamicsCreate) {
    ElastodynamicsKernel kernel;
    EXPECT_EQ(kernel.getType(), PhysicsType::ELASTODYNAMICS);
    EXPECT_EQ(kernel.getNumFields(), 1);
    EXPECT_EQ(kernel.getNumComponents(0), 3);
}

TEST_F(PhysicsKernelBaseTest, ElastodynamicsSetProperties) {
    ElastodynamicsKernel kernel;
    EXPECT_NO_THROW(kernel.setMaterialProperties(50e9, 0.25, 2500.0));
}

// ============================================================================
// ThermalKernel Tests
// ============================================================================

TEST_F(PhysicsKernelBaseTest, ThermalCreate) {
    ThermalKernel kernel;
    EXPECT_EQ(kernel.getType(), PhysicsType::THERMAL);
    EXPECT_EQ(kernel.getNumFields(), 1);
    EXPECT_EQ(kernel.getNumComponents(0), 1);
}

TEST_F(PhysicsKernelBaseTest, ThermalSetProperties) {
    ThermalKernel kernel;
    EXPECT_NO_THROW(kernel.setThermalProperties(2.5, 2500.0, 1000.0));
}

TEST_F(PhysicsKernelBaseTest, ThermalResidualFinite) {
    ThermalKernel kernel;
    kernel.setThermalProperties(2.5, 2500.0, 1000.0);
    
    PetscScalar u[1] = {350.0};
    PetscScalar u_t[1] = {0.1};
    PetscScalar u_x[3] = {10.0, 0.0, 0.0};
    PetscScalar a[1] = {0.0};
    PetscReal x[3] = {0.0, 0.0, 0.0};
    PetscScalar f[1] = {0.0};
    
    kernel.residual(u, u_t, u_x, a, x, f);
    
    EXPECT_TRUE(std::isfinite(f[0]));
}

// ============================================================================
// PoroelastodynamicsKernel Tests
// ============================================================================

TEST_F(PhysicsKernelBaseTest, PoroelastodynamicsCreate) {
    PoroelastodynamicsKernel kernel;
    EXPECT_EQ(kernel.getType(), PhysicsType::POROELASTODYNAMICS);
    EXPECT_EQ(kernel.getNumFields(), 2);  // Displacement + Pressure
    EXPECT_EQ(kernel.getNumComponents(0), 3);  // Displacement
    EXPECT_EQ(kernel.getNumComponents(1), 1);  // Pressure
}

TEST_F(PhysicsKernelBaseTest, PoroelastodynamicsSetProperties) {
    PoroelastodynamicsKernel kernel;
    EXPECT_NO_THROW(kernel.setMaterialProperties(20e9, 0.25, 2650.0, 0.2));
    EXPECT_NO_THROW(kernel.setFluidProperties(1000.0, 1e-3, 2.2e9));
    EXPECT_NO_THROW(kernel.setBiotParameters(0.8, 10e9));
}

// ============================================================================
// Physics Calculations Verification
// ============================================================================

TEST_F(PhysicsKernelBaseTest, DarcyVelocityCalculation) {
    // Verify Darcy's law: v = -(k/mu) * grad(p)
    double k = 100e-15;  // 100 mD
    double mu = 1e-3;    // 1 cP
    double grad_p = 1e5; // 100 kPa/m
    
    double v_darcy = -(k / mu) * grad_p;
    
    EXPECT_LT(v_darcy, 0.0) << "Flow in negative gradient direction";
    EXPECT_TRUE(std::isfinite(v_darcy));
}

TEST_F(PhysicsKernelBaseTest, LameParameterCalculation) {
    double E = 10e9;    // Young's modulus
    double nu = 0.25;   // Poisson's ratio
    
    double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    double mu = E / (2.0 * (1.0 + nu));
    
    EXPECT_GT(lambda, 0.0);
    EXPECT_GT(mu, 0.0);
    
    // Verify: K = lambda + 2/3*mu
    double K = lambda + 2.0/3.0 * mu;
    double K_expected = E / (3.0 * (1.0 - 2.0*nu));
    EXPECT_NEAR(K, K_expected, K * 1e-10);
}

TEST_F(PhysicsKernelBaseTest, WaveVelocityCalculation) {
    double E = 50e9;
    double nu = 0.25;
    double rho = 2500.0;
    
    double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    double mu = E / (2.0 * (1.0 + nu));
    
    double V_p = std::sqrt((lambda + 2.0*mu) / rho);
    double V_s = std::sqrt(mu / rho);
    
    EXPECT_GT(V_p, V_s) << "P-wave faster than S-wave";
    EXPECT_GT(V_p, 3000.0);
    EXPECT_GT(V_s, 1500.0);
}
