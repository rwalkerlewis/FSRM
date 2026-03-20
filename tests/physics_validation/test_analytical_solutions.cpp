/**
 * @file test_analytical_solutions.cpp
 * @brief 1D steady-state linear pressure: zero accumulation, zero Laplacian
 */

#include "physics/PhysicsKernel.hpp"
#include "core/FSRM.hpp"
#include <gtest/gtest.h>
#include <cmath>

using namespace FSRM;

class AnalyticalSolutionsTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    }
    int rank = 0;
};

TEST_F(AnalyticalSolutionsTest, SteadyLinearPressureZeroPointResidual) {
    const double p0 = 10.0e6;
    const double p1 = 20.0e6;
    const double L = 100.0;
    const PetscReal x[3] = {35.0, 0.0, 0.0};

    const double p = p0 + (p1 - p0) * (static_cast<double>(x[0]) / L);
    const double dp_dx = (p1 - p0) / L;

    SinglePhaseFlowKernel kernel;
    kernel.setProperties(0.2, 1.0e-13, 1.0e-9, 1.0e-3, 1000.0);

    PetscScalar u[1] = {PetscScalar(p)};
    PetscScalar u_t[1] = {0.0};
    PetscScalar u_x[3] = {PetscScalar(dp_dx), 0.0, 0.0};
    PetscScalar a[1] = {0.0};
    PetscScalar f[1] = {0.0};

    kernel.residual(u, u_t, u_x, a, x, f);

    EXPECT_NEAR(static_cast<double>(f[0]), 0.0, 1.0e-6);
    (void)rank;
}
