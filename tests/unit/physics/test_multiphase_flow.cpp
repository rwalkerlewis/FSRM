/**
 * @file test_multiphase_flow.cpp
 * @brief Unit tests for MultiphaseFlowKernel
 */
#include <gtest/gtest.h>
#include "physics/MultiphaseFlowKernel.hpp"
#include "core/FSRM.hpp"
#include <cmath>
#include <mpi.h>

using namespace FSRM;

class MultiphaseFlowKernelTest : public ::testing::Test {
protected:
    void SetUp() override { MPI_Comm_rank(PETSC_COMM_WORLD, &rank); }
    int rank = 0;
};

TEST_F(MultiphaseFlowKernelTest, TypeAndFieldLayout) {
    MultiphaseFlowKernel kernel;
    EXPECT_EQ(kernel.getType(), PhysicsType::FLUID_FLOW);
    EXPECT_EQ(kernel.getNumFields(), 2);
    EXPECT_EQ(kernel.getNumComponents(0), 1);
    EXPECT_EQ(kernel.getNumComponents(1), 1);
}

TEST_F(MultiphaseFlowKernelTest, ResidualFinite) {
    MultiphaseFlowKernel kernel;
    kernel.setRockProperties(0.15, 200e-15, 200e-15, 100e-15);
    kernel.setWaterProperties(1100.0, 8e-4, 4e-10);
    kernel.setCO2Properties(700.0, 5e-5, 1e-8);
    kernel.setRelPermParams(0.2, 0.05, 4.0, 2.0, 1.0, 1.0);

    PetscScalar u[2] = {20e6, 0.6};
    PetscScalar u_t[2] = {1e3, 0.001};
    PetscScalar u_x[6] = {};
    PetscScalar a_data[1] = {1.0};
    PetscReal x[3] = {0.5, 0.5, 0.5};
    PetscScalar f[2] = {};

    kernel.residual(u, u_t, u_x, a_data, x, f);
    EXPECT_TRUE(std::isfinite(PetscRealPart(f[0])));
    EXPECT_TRUE(std::isfinite(PetscRealPart(f[1])));
}

TEST_F(MultiphaseFlowKernelTest, JacobianDiagonalNonZero) {
    MultiphaseFlowKernel kernel;
    kernel.setRockProperties(0.15, 200e-15, 200e-15, 100e-15);

    PetscScalar u[2] = {20e6, 0.5};
    PetscScalar u_t[2] = {0.0, 0.0};
    PetscScalar u_x[6] = {};
    PetscScalar a_data[1] = {1.0};
    PetscReal x[3] = {};
    PetscScalar J[4] = {};

    kernel.jacobian(u, u_t, u_x, a_data, x, J);
    EXPECT_NE(PetscRealPart(J[0]), 0.0);  // ∂f0/∂P
    EXPECT_NE(PetscRealPart(J[3]), 0.0);  // ∂f1/∂Sw
}

TEST_F(MultiphaseFlowKernelTest, AccumulationOnlyWithZeroGradient) {
    // With zero pressure gradient, flux terms should be zero
    // Residual should only contain accumulation
    MultiphaseFlowKernel kernel;

    PetscScalar u[2] = {20e6, 0.5};
    PetscScalar u_t[2] = {1000.0, 0.01};
    PetscScalar u_x[6] = {0, 0, 0, 0, 0, 0};
    PetscScalar a_data[1] = {1.0};
    PetscReal x[3] = {};
    PetscScalar f[2] = {};

    kernel.residual(u, u_t, u_x, a_data, x, f);
    // f[0] should be phi * ct * P_t (non-zero)
    EXPECT_NE(PetscRealPart(f[0]), 0.0);
    // f[1] should be phi * S_t (non-zero)
    EXPECT_NE(PetscRealPart(f[1]), 0.0);
}
