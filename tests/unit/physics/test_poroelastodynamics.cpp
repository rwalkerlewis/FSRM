/**
 * @file test_poroelastodynamics.cpp
 * @brief Unit tests for PoroelastodynamicsKernel
 */

#include <gtest/gtest.h>
#include "physics/PhysicsKernel.hpp"
#include "core/FSRM.hpp"
#include <cmath>
#include <mpi.h>

using namespace FSRM;

class PoroelastodynamicsKernelTest : public ::testing::Test {
protected:
    void SetUp() override { MPI_Comm_rank(PETSC_COMM_WORLD, &rank); }
    int rank = 0;
};

TEST_F(PoroelastodynamicsKernelTest, TypeAndFieldLayout) {
    PoroelastodynamicsKernel kernel;
    EXPECT_EQ(kernel.getType(), PhysicsType::POROELASTODYNAMICS);
    EXPECT_EQ(kernel.getNumFields(), 2);
    EXPECT_EQ(kernel.getNumComponents(0), 3);
    EXPECT_EQ(kernel.getNumComponents(1), 1);
}

TEST_F(PoroelastodynamicsKernelTest, ResidualJacobianSmoke) {
    PoroelastodynamicsKernel kernel;
    kernel.setMaterialProperties(10e9, 0.25, 2500.0, 0.2);
    kernel.setFluidProperties(1000.0, 0.001, 2.2e9);
    kernel.setBiotParameters(0.9, 1e10);

    PetscScalar u[4] = {};
    PetscScalar u_t[4] = {};
    PetscScalar u_x[12] = {};
    PetscScalar a[1] = {1.0};
    PetscReal x[3] = {0.5, 0.5, 0.5};
    PetscScalar f[4] = {};
    PetscScalar J[196] = {};

    kernel.residual(u, u_t, u_x, a, x, f);
    for (int i = 0; i < 4; ++i) {
        EXPECT_TRUE(std::isfinite(PetscRealPart(f[i])));
    }

    u_x[0] = 0.001;
    u_x[9] = 1e5;
    for (int i = 0; i < 4; ++i) {
        f[i] = {};
    }
    kernel.residual(u, u_t, u_x, a, x, f);
    double s = 0.0;
    for (int i = 0; i < 4; ++i) {
        s += std::abs(PetscRealPart(f[i]));
    }
    EXPECT_GT(s, 0.0);

    for (int i = 0; i < 12; ++i) {
        u_x[i] = {};
    }
    for (int i = 0; i < 196; ++i) {
        J[i] = {};
    }
    kernel.jacobian(u, u_t, u_x, a, x, J);
    double jsum = 0.0;
    for (int i = 0; i < 196; ++i) {
        jsum += std::abs(PetscRealPart(J[i]));
    }
    EXPECT_GT(jsum, 0.0);
}
