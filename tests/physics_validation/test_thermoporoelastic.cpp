/**
 * @file test_thermoporoelastic.cpp
 * @brief Terzaghi-style poroelastic parameters; coupled kernel non-trivial residual
 */

#include "physics/PhysicsKernel.hpp"
#include "core/FSRM.hpp"
#include <gtest/gtest.h>
#include <cmath>

using namespace FSRM;

class ThermoporoelasticConsolidationTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    }
    int rank = 0;
};

TEST_F(ThermoporoelasticConsolidationTest, PoroelastodynamicsResidualNonTrivial) {
    PoroelastodynamicsKernel kernel;
    kernel.setMaterialProperties(10.0e9, 0.25, 2650.0, 0.25);
    kernel.setFluidProperties(1000.0, 1.0e-3, 2.2e9);
    kernel.setBiotParameters(0.9, 1.0e10);

    PetscScalar u[4] = {0.0};
    PetscScalar u_t[4] = {0.0};
    PetscScalar u_x[12] = {0.0};
    PetscScalar a[4] = {0.0};
    PetscScalar f[4] = {0.0};
    PetscReal x[3] = {0.0, 0.0, 0.0};

    u[0] = PetscScalar(1.0e-5);
    u[1] = PetscScalar(0.5e-5);
    u[2] = PetscScalar(0.0);
    u[3] = PetscScalar(5.0e6);

    u_x[0] = PetscScalar(1.0e-4);
    u_x[4] = PetscScalar(-0.5e-4);
    u_x[8] = PetscScalar(0.0);
    u_x[9] = PetscScalar(1.0e3);
    u_x[10] = PetscScalar(0.0);
    u_x[11] = PetscScalar(0.0);

    u_t[0] = PetscScalar(0.1);
    u_t[1] = PetscScalar(0.0);
    u_t[2] = PetscScalar(0.0);
    u_t[3] = PetscScalar(1.0e4);

    kernel.residual(u, u_t, u_x, a, x, f);

    const double mag = std::sqrt(static_cast<double>(f[0]) * static_cast<double>(f[0]) +
                                 static_cast<double>(f[1]) * static_cast<double>(f[1]) +
                                 static_cast<double>(f[2]) * static_cast<double>(f[2]) +
                                 static_cast<double>(f[3]) * static_cast<double>(f[3]));
    EXPECT_GT(mag, 1.0);
    (void)rank;
}
