/**
 * @file test_geomechanics.cpp
 * @brief Unit tests for GeomechanicsKernel
 */

#include <gtest/gtest.h>
#include "physics/PhysicsKernel.hpp"
#include "core/FSRM.hpp"
#include <cmath>
#include <mpi.h>

using namespace FSRM;

class GeomechanicsKernelTest : public ::testing::Test {
protected:
    void SetUp() override { MPI_Comm_rank(PETSC_COMM_WORLD, &rank); }
    int rank = 0;
};

TEST_F(GeomechanicsKernelTest, ElasticKernelBasics) {
    GeomechanicsKernel kernel(SolidModelType::ELASTIC);
    EXPECT_EQ(kernel.getType(), PhysicsType::GEOMECHANICS);
    EXPECT_EQ(kernel.getNumFields(), 1);
    EXPECT_EQ(kernel.getNumComponents(0), 3);
}

TEST_F(GeomechanicsKernelTest, ResidualAndJacobian) {
    GeomechanicsKernel kernel(SolidModelType::ELASTIC);
    kernel.setMaterialProperties(10e9, 0.25, 2500.0);

    PetscScalar u[3] = {};
    PetscScalar u_t[3] = {};
    PetscScalar u_x[9] = {};
    PetscScalar a[1] = {1.0};
    PetscReal x[3] = {0.5, 0.5, 0.5};
    PetscScalar f[3] = {};
    PetscScalar J[81] = {};

    kernel.residual(u, u_t, u_x, a, x, f);
    EXPECT_TRUE(std::isfinite(PetscRealPart(f[0])));

    u_x[0] = 0.01;
    f[0] = f[1] = f[2] = {};
    kernel.residual(u, u_t, u_x, a, x, f);
    const double fn = std::hypot(PetscRealPart(f[0]),
                                 std::hypot(PetscRealPart(f[1]), PetscRealPart(f[2])));
    EXPECT_GT(fn, 0.0);

    u_x[0] = {};
    for (int i = 0; i < 81; ++i) {
        J[i] = {};
    }
    kernel.jacobian(u, u_t, u_x, a, x, J);
    double jn = 0.0;
    for (int i = 0; i < 81; ++i) {
        jn += std::abs(PetscRealPart(J[i]));
    }
    EXPECT_GT(jn, 0.0);
}
