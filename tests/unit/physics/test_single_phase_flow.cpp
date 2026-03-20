/**
 * @file test_single_phase_flow.cpp
 * @brief Unit tests for SinglePhaseFlowKernel
 */

#include <gtest/gtest.h>
#include "physics/PhysicsKernel.hpp"
#include "core/FSRM.hpp"
#include <cmath>
#include <mpi.h>

using namespace FSRM;

class SinglePhaseFlowKernelTest : public ::testing::Test {
protected:
    void SetUp() override { MPI_Comm_rank(PETSC_COMM_WORLD, &rank); }
    int rank = 0;
};

TEST_F(SinglePhaseFlowKernelTest, TypeAndFieldLayout) {
    SinglePhaseFlowKernel kernel;
    EXPECT_EQ(kernel.getType(), PhysicsType::FLUID_FLOW);
    EXPECT_EQ(kernel.getNumFields(), 1);
    EXPECT_EQ(kernel.getNumComponents(0), 1);
}

TEST_F(SinglePhaseFlowKernelTest, SetPropertiesResidualJacobian) {
    SinglePhaseFlowKernel kernel;
    kernel.setProperties(0.2, 1e-13, 1e-9, 0.001, 1000.0);

    PetscScalar u[1] = {};
    PetscScalar u_t[1] = {};
    PetscScalar u_x[3] = {};
    PetscScalar a[1] = {1.0};
    PetscReal x[3] = {0.5, 0.5, 0.5};
    PetscScalar f[1] = {};
    PetscScalar J[1] = {};

    kernel.residual(u, u_t, u_x, a, x, f);
    EXPECT_TRUE(std::isfinite(PetscRealPart(f[0])));

    u_x[0] = 1e5;
    u_x[1] = 0.0;
    u_x[2] = 0.0;
    f[0] = {};
    kernel.residual(u, u_t, u_x, a, x, f);
    EXPECT_NE(PetscRealPart(f[0]), 0.0);

    u_x[0] = {};
    J[0] = {};
    kernel.jacobian(u, u_t, u_x, a, x, J);
    EXPECT_NE(PetscRealPart(J[0]), 0.0);
}
