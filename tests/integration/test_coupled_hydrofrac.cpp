/**
 * @file test_coupled_hydrofrac.cpp
 * @brief Phase 3 coupling validation checks for hydrofracture utilities
 */

#include <gtest/gtest.h>

#include "numerics/PetscFEHydrofrac.hpp"

using namespace FSRM;

class CoupledHydrofracTest : public ::testing::Test
{
};

TEST_F(CoupledHydrofracTest, WidthScalingIncreasesWithInjectionTime)
{
  const PetscReal mu = 1.0e-3;
  const PetscReal q = 0.02;
  const PetscReal eprime = 20.0e9;

  const PetscReal w_t1 = PetscFEHydrofrac::pknWidthScaling(mu, q, 10.0, eprime);
  const PetscReal w_t2 = PetscFEHydrofrac::pknWidthScaling(mu, q, 40.0, eprime);

  EXPECT_GT(w_t1, 0.0);
  EXPECT_GT(w_t2, w_t1);
}

TEST_F(CoupledHydrofracTest, WidthScalingMatchesQuarterPowerRatio)
{
  const PetscReal mu = 2.0e-3;
  const PetscReal q = 0.03;
  const PetscReal eprime = 18.0e9;

  const PetscReal w1 = PetscFEHydrofrac::pknWidthScaling(mu, q, 5.0, eprime);
  const PetscReal w2 = PetscFEHydrofrac::pknWidthScaling(mu, q, 80.0, eprime);

  const PetscReal ratio = w2 / w1;
  const PetscReal expected = PetscPowReal(80.0 / 5.0, 0.25);
  EXPECT_NEAR(ratio, expected, 1e-12);
}

TEST_F(CoupledHydrofracTest, InvalidInputsReturnZeroWidth)
{
  EXPECT_EQ(PetscFEHydrofrac::pknWidthScaling(0.0, 1.0, 1.0, 1.0), 0.0);
  EXPECT_EQ(PetscFEHydrofrac::pknWidthScaling(1.0, -1.0, 1.0, 1.0), 0.0);
  EXPECT_EQ(PetscFEHydrofrac::pknWidthScaling(1.0, 1.0, 0.0, 1.0), 0.0);
  EXPECT_EQ(PetscFEHydrofrac::pknWidthScaling(1.0, 1.0, 1.0, 0.0), 0.0);
}
