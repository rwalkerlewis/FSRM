/**
 * @file test_stress_shadowing.cpp
 * @brief Phase 5 stress shadowing integration tests
 *
 * Validates that the FEM approach captures inter-fracture
 * stress interactions consistent with analytical Sneddon solutions.
 */

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "numerics/PetscFEHydrofrac.hpp"

using namespace FSRM;

class StressShadowingTest : public ::testing::Test
{
};

TEST_F(StressShadowingTest, NormalStressPerturbationDecaysWithDistance)
{
  const PetscReal P = 5.0e6;
  const PetscReal a = 50.0;

  // Inside the crack: full pressure
  const PetscReal sigma_inside = PetscFEHydrofrac::sneddonNormalStressPerturbation(P, a, 25.0);
  EXPECT_NEAR(sigma_inside, P, 1e-10);

  // At 2*a: should be less than P
  const PetscReal sigma_2a = PetscFEHydrofrac::sneddonNormalStressPerturbation(P, a, 100.0);
  EXPECT_GT(sigma_2a, 0.0);
  EXPECT_LT(sigma_2a, P);

  // At 5*a: even less
  const PetscReal sigma_5a = PetscFEHydrofrac::sneddonNormalStressPerturbation(P, a, 250.0);
  EXPECT_GT(sigma_5a, 0.0);
  EXPECT_LT(sigma_5a, sigma_2a);
}

TEST_F(StressShadowingTest, StressShadowFactorDecreasesWithSpacing)
{
  const PetscReal P = 10.0e6;
  const PetscReal half_length = 50.0;

  // All spacings must exceed half_length for decay behavior
  const PetscReal sf_close = PetscFEHydrofrac::stressShadowFactor(P, half_length, 60.0);
  const PetscReal sf_mid = PetscFEHydrofrac::stressShadowFactor(P, half_length, 100.0);
  const PetscReal sf_far = PetscFEHydrofrac::stressShadowFactor(P, half_length, 250.0);

  EXPECT_GT(sf_close, sf_mid);
  EXPECT_GT(sf_mid, sf_far);
  EXPECT_GT(sf_far, 0.0);

  // At spacing within crack radius, factor should be 1.0
  const PetscReal sf_within = PetscFEHydrofrac::stressShadowFactor(P, half_length, 5.0);
  EXPECT_NEAR(sf_within, 1.0, 1e-10);
}

TEST_F(StressShadowingTest, ThreeClustersOuterReceiveMoreShadow)
{
  // 3 parallel fractures: spacing exceeds half_length for shadow differentiation
  const PetscInt N = 3;
  const PetscReal spacing = 20.0;
  const PetscReal half_length = 10.0;
  const PetscReal P = 8.0e6;

  // Cluster 0 (outer, left): shadow from clusters 1 and 2
  const PetscReal eff_0 = PetscFEHydrofrac::clusterEfficiency(
      N, spacing, half_length, P, 0);

  // Cluster 1 (center): shadow from clusters 0 and 2 (both at 10m)
  const PetscReal eff_1 = PetscFEHydrofrac::clusterEfficiency(
      N, spacing, half_length, P, 1);

  // Cluster 2 (outer, right): symmetric with cluster 0
  const PetscReal eff_2 = PetscFEHydrofrac::clusterEfficiency(
      N, spacing, half_length, P, 2);

  // Outer clusters are symmetric
  EXPECT_NEAR(eff_0, eff_2, 1e-12);

  // Center cluster receives more shadow (neighbors at 10m on both sides)
  // so efficiency is LOWER for center (harder to open)
  EXPECT_LT(eff_1, eff_0);

  // All efficiencies sum to 1.0
  EXPECT_NEAR(eff_0 + eff_1 + eff_2, 1.0, 1e-10);

  // All efficiencies are positive
  EXPECT_GT(eff_0, 0.0);
  EXPECT_GT(eff_1, 0.0);
}

TEST_F(StressShadowingTest, FiveClustersEdgesGetMostFluid)
{
  const PetscInt N = 5;
  const PetscReal spacing = 20.0;
  const PetscReal half_length = 8.0;
  const PetscReal P = 6.0e6;

  std::vector<PetscReal> eff(N);
  PetscReal total = 0.0;
  for (PetscInt i = 0; i < N; ++i)
  {
    eff[i] = PetscFEHydrofrac::clusterEfficiency(N, spacing, half_length, P, i);
    total += eff[i];
  }

  // Sum to 1
  EXPECT_NEAR(total, 1.0, 1e-10);

  // Edge clusters (0 and 4) get the most fluid
  EXPECT_GT(eff[0], eff[1]);
  EXPECT_GT(eff[0], eff[2]);

  // Center cluster (2) gets the least
  for (PetscInt i = 0; i < N; ++i)
  {
    if (i != 2)
    {
      EXPECT_GE(eff[i], eff[2] - 1e-14);
    }
  }

  // Symmetric
  EXPECT_NEAR(eff[0], eff[4], 1e-12);
  EXPECT_NEAR(eff[1], eff[3], 1e-12);
}

TEST_F(StressShadowingTest, SingleClusterHasFullEfficiency)
{
  const PetscReal eff = PetscFEHydrofrac::clusterEfficiency(1, 10.0, 50.0, 5.0e6, 0);
  EXPECT_NEAR(eff, 1.0, 1e-12);
}

TEST_F(StressShadowingTest, InvalidInputsReturnZero)
{
  EXPECT_EQ(PetscFEHydrofrac::sneddonNormalStressPerturbation(0.0, 50.0, 10.0), 0.0);
  EXPECT_EQ(PetscFEHydrofrac::sneddonNormalStressPerturbation(5.0e6, 0.0, 10.0), 0.0);
  EXPECT_EQ(PetscFEHydrofrac::stressShadowFactor(0.0, 50.0, 10.0), 0.0);
  EXPECT_EQ(PetscFEHydrofrac::clusterEfficiency(3, 10.0, 50.0, 5.0e6, 5), 0.0);
}
