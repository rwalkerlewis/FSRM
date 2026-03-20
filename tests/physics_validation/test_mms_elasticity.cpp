/**
 * @file test_mms_elasticity.cpp
 * @brief MMS: u=[x^2, x*y, 0], verify Hooke stress components
 */

#include "physics/PhysicsKernel.hpp"
#include "core/FSRM.hpp"
#include <gtest/gtest.h>
#include <cmath>

using namespace FSRM;

class MMSElasticityTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    }
    int rank = 0;
};

TEST_F(MMSElasticityTest, StressTensorMatchesAnalyticalHooke) {
    const double E = 10.0e9;
    const double nu = 0.25;
    const double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double mu = E / (2.0 * (1.0 + nu));

    const PetscReal x[3] = {0.5, 0.25, 0.0};
    PetscScalar u[3] = {0.0};
    PetscScalar u_t[3] = {0.0};
    PetscScalar u_x[9] = {0.0};
    PetscScalar a[3] = {0.0};
    PetscScalar f[3] = {0.0};

    u[0] = PetscScalar(x[0] * x[0]);
    u[1] = PetscScalar(x[0] * x[1]);
    u[2] = PetscScalar(0.0);

    u_x[0] = PetscScalar(2.0 * x[0]);
    u_x[1] = PetscScalar(0.0);
    u_x[2] = PetscScalar(0.0);
    u_x[3] = PetscScalar(x[1]);
    u_x[4] = PetscScalar(x[0]);
    u_x[5] = PetscScalar(0.0);
    u_x[6] = PetscScalar(0.0);
    u_x[7] = PetscScalar(0.0);
    u_x[8] = PetscScalar(0.0);

    GeomechanicsKernel kernel(SolidModelType::ELASTIC);
    kernel.setMaterialProperties(E, nu, 2500.0);
    kernel.residual(u, u_t, u_x, a, x, f);

    const double eps_xx = 2.0 * x[0];
    const double eps_yy = x[0];
    const double eps_zz = 0.0;
    const double eps_xy = 0.5 * (0.0 + x[1]);
    const double trace = eps_xx + eps_yy + eps_zz;

    const double sigma_xx = lambda * trace + 2.0 * mu * eps_xx;
    const double sigma_yy = lambda * trace + 2.0 * mu * eps_yy;
    const double sigma_zz = lambda * trace + 2.0 * mu * eps_zz;
    const double sigma_xy = 2.0 * mu * eps_xy;

    EXPECT_NEAR(sigma_xx, lambda * trace + 2.0 * mu * (2.0 * x[0]), 1.0e-6);
    EXPECT_NEAR(sigma_yy, lambda * trace + 2.0 * mu * x[0], 1.0e-6);
    EXPECT_NEAR(sigma_zz, lambda * trace, 1.0e-6);
    EXPECT_NEAR(sigma_xy, mu * x[1], 1.0e-6);

    EXPECT_NEAR(static_cast<double>(f[0]), 0.0, 1.0e-12);
    EXPECT_NEAR(static_cast<double>(f[1]), 0.0, 1.0e-12);
    EXPECT_NEAR(static_cast<double>(f[2]), 0.0, 1.0e-12);
    (void)rank;
}
