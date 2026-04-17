/**
 * @file test_fracture_propagation.cpp
 * @brief Phase 4 propagation criterion validation checks
 */

#include <gtest/gtest.h>

#include "numerics/PetscFEHydrofrac.hpp"

using namespace FSRM;

class FracturePropagationTest : public ::testing::Test
{
};

TEST_F(FracturePropagationTest, CohesiveStrengthScalesWithInverseSqrtH)
{
  const PetscReal k_ic = 1.5e6;
  const PetscReal t_h1 = PetscFEHydrofrac::cohesiveStrengthFromKIc(k_ic, 0.5);
  const PetscReal t_h2 = PetscFEHydrofrac::cohesiveStrengthFromKIc(k_ic, 2.0);

  EXPECT_GT(t_h1, t_h2);
  EXPECT_GT(t_h2, 0.0);

  const PetscReal expected_ratio = PetscSqrtReal(2.0 / 0.5);
  EXPECT_NEAR(t_h1 / t_h2, expected_ratio, 1e-12);
}

TEST_F(FracturePropagationTest, OpeningCriterionTriggersAboveTensileStrength)
{
  const PetscReal tensile_strength = 5.0e6;

  EXPECT_EQ(PetscFEHydrofrac::shouldOpenFractureCell(4.9e6, tensile_strength), PETSC_FALSE);
  EXPECT_EQ(PetscFEHydrofrac::shouldOpenFractureCell(5.0e6, tensile_strength), PETSC_FALSE);
  EXPECT_EQ(PetscFEHydrofrac::shouldOpenFractureCell(5.1e6, tensile_strength), PETSC_TRUE);
}

TEST_F(FracturePropagationTest, InvalidFractureToughnessInputReturnsZero)
{
  EXPECT_EQ(PetscFEHydrofrac::cohesiveStrengthFromKIc(0.0, 1.0), 0.0);
  EXPECT_EQ(PetscFEHydrofrac::cohesiveStrengthFromKIc(1.0e6, 0.0), 0.0);
}
