/**
 * @file test_elastodynamics.cpp
 * @brief Unit tests for ElastodynamicsKernel
 */

#include <gtest/gtest.h>
#include "physics/PhysicsKernel.hpp"
#include "core/FSRM.hpp"
#include <cmath>
#include <mpi.h>

using namespace FSRM;

class ElastodynamicsKernelTest : public ::testing::Test {
protected:
    void SetUp() override { MPI_Comm_rank(PETSC_COMM_WORLD, &rank); }
    int rank = 0;
};

TEST_F(ElastodynamicsKernelTest, TypeFieldsAndProperties) {
    ElastodynamicsKernel kernel;
    EXPECT_EQ(kernel.getType(), PhysicsType::ELASTODYNAMICS);
    EXPECT_EQ(kernel.getNumFields(), 1);
    EXPECT_EQ(kernel.getNumComponents(0), 3);
    kernel.setMaterialProperties(10e9, 0.25, 2500.0);
    kernel.setWaveProperties(5000.0, 3000.0, 100.0);
    kernel.setDamping(0.01, 0.001);
}

TEST_F(ElastodynamicsKernelTest, ResidualJacobian) {
    ElastodynamicsKernel kernel;
    kernel.setMaterialProperties(10e9, 0.25, 2500.0);
    kernel.setWaveProperties(5000.0, 3000.0, 100.0);
    kernel.setDamping(0.01, 0.001);

    PetscScalar u[3] = {};
    PetscScalar u_t[3] = {};
    PetscScalar u_x[9] = {};
    PetscScalar a[1] = {1.0};
    PetscReal x[3] = {0.5, 0.5, 0.5};
    PetscScalar f[3] = {};
    PetscScalar J[100] = {};

    kernel.residual(u, u_t, u_x, a, x, f);
    EXPECT_TRUE(std::isfinite(PetscRealPart(f[0])));

    u_x[0] = 0.001;
    f[0] = f[1] = f[2] = {};
    kernel.residual(u, u_t, u_x, a, x, f);
    double fn = 0.0;
    for (int i = 0; i < 3; ++i) {
        fn += std::abs(PetscRealPart(f[i]));
    }
    EXPECT_GT(fn, 0.0);

    u_x[0] = {};
    for (int i = 0; i < 100; ++i) {
        J[i] = {};
    }
    kernel.jacobian(u, u_t, u_x, a, x, J);
    double jsum = 0.0;
    for (int i = 0; i < 100; ++i) {
        jsum += std::abs(PetscRealPart(J[i]));
    }
    EXPECT_GT(jsum, 0.0);
}
