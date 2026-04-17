/**
 * @file test_cohesive_fault_kernel.cpp
 * @brief Tests for CohesiveFaultKernel PetscDS callbacks and mode switching
 */

#include <gtest/gtest.h>
#include <cmath>
#include <petscsys.h>

#include "physics/CohesiveFaultKernel.hpp"

using namespace FSRM;

class CohesiveFaultKernelTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    }

    int rank = 0;
};

TEST_F(CohesiveFaultKernelTest, LockedLagrangeConstraintReturnsDisplacementJump) {
    constexpr PetscInt dim = 3;
    PetscInt uOff[] = {0, 6};
    PetscScalar u[9] = {};
    // Negative side at rest; positive side displaced in x
    u[3] = 2e-4;
    u[4] = -1e-4;
    u[5] = 3e-5;

    PetscScalar f[3] = {};
    PetscReal n[] = {0.0, 0.0, 1.0};
    // Provide constants up to COHESIVE_CONST_COUNT; slot 24 = 0.0 (locked)
    PetscScalar constants[CohesiveFaultKernel::COHESIVE_CONST_COUNT] = {};
    constants[CohesiveFaultKernel::COHESIVE_CONST_MODE] = 0.0;  // locked
    constants[CohesiveFaultKernel::COHESIVE_CONST_MU_F] = 0.6;

    CohesiveFaultKernel::f0_lagrange_constraint(dim, 0, 0, uOff, nullptr, u, nullptr,
                                                 nullptr, nullptr, nullptr, nullptr, nullptr,
                                                 nullptr, 0.0, nullptr, n,
                                                 CohesiveFaultKernel::COHESIVE_CONST_COUNT,
                                                 constants, f);

    EXPECT_NEAR(PetscRealPart(f[0]), PetscRealPart(u[3] - u[0]), 1e-15);
    EXPECT_NEAR(PetscRealPart(f[1]), PetscRealPart(u[4] - u[1]), 1e-15);
    EXPECT_NEAR(PetscRealPart(f[2]), PetscRealPart(u[5] - u[2]), 1e-15);
}

TEST_F(CohesiveFaultKernelTest, ModeSwitchingUpdatesIsLockedFlag) {
    CohesiveFaultKernel kernel;
    EXPECT_TRUE(kernel.isLocked());
    kernel.setMode(true);
    EXPECT_TRUE(kernel.isLocked());
    kernel.setMode(false);
    EXPECT_FALSE(kernel.isLocked());
}

TEST_F(CohesiveFaultKernelTest, SlippingModeUsesFrictionStrengthInConstraint) {
    constexpr PetscInt dim = 3;
    PetscInt uOff[] = {0, 6};
    // Tangential slip along x, fault normal +z
    PetscScalar u[9] = {};
    u[3] = 1e-3;  // positive side moves in x
    u[6] = 0.0;
    u[7] = 0.0;
    u[8] = -50e6;  // compressive normal traction carried by lambda_z (sign per implementation)

    PetscScalar f[3] = {};
    PetscReal n[] = {0.0, 0.0, 1.0};
    // Provide constants up to COHESIVE_CONST_COUNT; slot 24 = 1.0 (slipping), slot 25 = mu_f
    PetscScalar constants[CohesiveFaultKernel::COHESIVE_CONST_COUNT] = {};
    constants[CohesiveFaultKernel::COHESIVE_CONST_MODE] = 1.0;  // slipping
    constants[CohesiveFaultKernel::COHESIVE_CONST_MU_F] = 0.6;

    CohesiveFaultKernel::f0_lagrange_constraint(dim, 0, 0, uOff, nullptr, u, nullptr,
                                                 nullptr, nullptr, nullptr, nullptr, nullptr,
                                                 nullptr, 0.0, nullptr, n,
                                                 CohesiveFaultKernel::COHESIVE_CONST_COUNT,
                                                 constants, f);

    // With only x-displacement jump, tangential slip is along x; lambda_t ≈ (lambda_x, lambda_y, 0)
    // Implementation uses lambda for shear balance — residual should be finite and not match locked jump
    const double f0_locked_x = u[3] - u[0];
    EXPECT_GT(std::abs(PetscRealPart(f[0] - f0_locked_x)), 1e-20)
        << "Slipping branch should differ from pure displacement-jump (locked) residual";
}
