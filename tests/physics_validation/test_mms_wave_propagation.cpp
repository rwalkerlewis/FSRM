/**
 * @file test_mms_wave_propagation.cpp
 * @brief Plane wave u_x = A*sin(k*x - omega*t): residual terms at a quadrature point
 */

#include "physics/PhysicsKernel.hpp"
#include "core/FSRM.hpp"
#include <gtest/gtest.h>
#include <cmath>

using namespace FSRM;

class MMSWavePropagationTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    }
    int rank = 0;
};

TEST_F(MMSWavePropagationTest, PlaneWaveInertialResidualAtPoint) {
    const double A = 1.0e-6;
    const double k = 1.0;
    const double omega = 314.0;
    const double t = 0.01;
    const PetscReal x[3] = {0.1, 0.0, 0.0};

    const double phase = k * static_cast<double>(x[0]) - omega * t;
    const double ux = A * std::sin(phase);
    const double dux_dx = A * k * std::cos(phase);
    const double vx = -A * omega * std::cos(phase);

    PetscScalar u[3] = {0.0};
    PetscScalar u_t[3] = {0.0};
    PetscScalar u_x[9] = {0.0};
    PetscScalar a[3] = {0.0};
    PetscScalar f[3] = {0.0};

    u[0] = PetscScalar(ux);
    u_t[0] = PetscScalar(vx);

    u_x[0] = PetscScalar(dux_dx);
    for (int i = 1; i < 9; ++i) {
        u_x[i] = PetscScalar(0.0);
    }

    ElastodynamicsKernel kernel;
    kernel.setMaterialProperties(10.0e9, 0.25, 2500.0);
    kernel.setDamping(0.01, 0.001);
    kernel.residual(u, u_t, u_x, a, x, f);

    const double rho = 2500.0;
    const double alpha = 0.01;
    const double expected_fx = rho * vx * (1.0 + alpha);

    EXPECT_NEAR(static_cast<double>(f[0]), expected_fx, 1.0e-12);
    EXPECT_TRUE(std::isfinite(static_cast<double>(f[1])));
    EXPECT_TRUE(std::isfinite(static_cast<double>(f[2])));
    (void)rank;
}
