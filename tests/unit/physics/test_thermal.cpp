/**
 * @file test_thermal.cpp
 * @brief Unit tests for ThermalKernel
 */

#include <gtest/gtest.h>
#include "physics/PhysicsKernel.hpp"
#include "core/FSRM.hpp"
#include <cmath>
#include <mpi.h>

using namespace FSRM;

class ThermalKernelTest : public ::testing::Test {
protected:
    void SetUp() override { MPI_Comm_rank(PETSC_COMM_WORLD, &rank); }
    int rank = 0;
};

TEST_F(ThermalKernelTest, TypeAndFields) {
    ThermalKernel kernel;
    EXPECT_EQ(kernel.getType(), PhysicsType::THERMAL);
    EXPECT_EQ(kernel.getNumFields(), 1);
    EXPECT_EQ(kernel.getNumComponents(0), 1);
}

TEST_F(ThermalKernelTest, ResidualJacobianSmoke) {
    ThermalKernel kernel;
    kernel.setThermalProperties(2.5, 2500.0, 900.0);

    PetscScalar u[1] = {};
    PetscScalar u_t[1] = {};
    PetscScalar u_x[3] = {};
    PetscScalar a[1] = {1.0};
    PetscReal x[3] = {0.5, 0.5, 0.5};
    PetscScalar f[1] = {};
    PetscScalar J[10] = {};

    kernel.residual(u, u_t, u_x, a, x, f);
    EXPECT_TRUE(std::isfinite(PetscRealPart(f[0])));

    u_x[0] = 100.0;
    f[0] = {};
    kernel.residual(u, u_t, u_x, a, x, f);
    EXPECT_NE(PetscRealPart(f[0]), 0.0);

    u_x[0] = {};
    for (int i = 0; i < 10; ++i) {
        J[i] = {};
    }
    kernel.jacobian(u, u_t, u_x, a, x, J);
    EXPECT_NE(PetscRealPart(J[0]), 0.0);
}
