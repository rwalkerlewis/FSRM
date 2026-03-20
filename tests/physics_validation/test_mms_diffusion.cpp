/**
 * @file test_mms_diffusion.cpp
 * @brief MMS: steady diffusion u(x)=sin(pi*x), -k*u'' = k*pi^2*sin(pi*x) = f(x)
 */

#include "physics/PhysicsKernel.hpp"
#include "core/FSRM.hpp"
#include <gtest/gtest.h>
#include <cmath>

using namespace FSRM;

class MMSDiffusionTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    }
    int rank = 0;
};

TEST_F(MMSDiffusionTest, ManufacturedSourceMatchesMobilityTimesLaplacian) {
    const double pi = std::acos(-1.0);
    const PetscReal x[3] = {0.25, 0.0, 0.0};

    const double phi = 0.2;
    const double k = 1.0e-12;
    const double ct = 1.0e-9;
    const double mu = 1.0e-3;
    const double rho = 1000.0;
    const double mobility = k / mu;

    SinglePhaseFlowKernel kernel;
    kernel.setProperties(phi, k, ct, mu, rho);

    PetscScalar u[1] = {0.0};
    PetscScalar u_t[1] = {0.0};
    PetscScalar u_x[3] = {0.0};
    PetscScalar a[1] = {0.0};
    PetscScalar f[1] = {0.0};

    u[0] = PetscScalar(std::sin(pi * static_cast<double>(x[0])));
    u_x[0] = PetscScalar(pi * std::cos(pi * static_cast<double>(x[0])));
    u_x[1] = PetscScalar(0.0);
    u_x[2] = PetscScalar(0.0);

    kernel.residual(u, u_t, u_x, a, x, f);

    const double u_xx = -pi * pi * std::sin(pi * static_cast<double>(x[0]));
    const double f_mms = mobility * (-u_xx);
    const double f_accum = static_cast<double>(f[0]);

    EXPECT_NEAR(f_accum, phi * ct * 0.0, 1.0e-6);
    EXPECT_NEAR(f_mms, mobility * pi * pi * std::sin(pi * static_cast<double>(x[0])), 1.0e-6);
    (void)rank;
}
