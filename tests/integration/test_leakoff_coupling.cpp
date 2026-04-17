/**
 * @file test_leakoff_coupling.cpp
 * @brief Phase 8 Carter leak-off coupling integration tests
 *
 * Validates Carter leak-off rate, cumulative volume,
 * and sqrt(t) scaling for fracture-to-formation coupling.
 */

#include <gtest/gtest.h>

#include <cmath>

#include "numerics/PetscFEHydrofrac.hpp"

using namespace FSRM;

class LeakoffCouplingTest : public ::testing::Test
{
};

TEST_F(LeakoffCouplingTest, CarterRateDecaysWithSqrtTime)
{
  const PetscReal c_l = 1.0e-4;   // m/sqrt(s)
  const PetscReal t_open = 0.0;

  const PetscReal q1 = PetscFEHydrofrac::carterLeakoffRate(c_l, 1.0, t_open);
  const PetscReal q4 = PetscFEHydrofrac::carterLeakoffRate(c_l, 4.0, t_open);
  const PetscReal q16 = PetscFEHydrofrac::carterLeakoffRate(c_l, 16.0, t_open);

  // q(t) = C_L / sqrt(t), so q(4)/q(1) = 1/2, q(16)/q(1) = 1/4
  EXPECT_NEAR(q4 / q1, 0.5, 1e-12);
  EXPECT_NEAR(q16 / q1, 0.25, 1e-12);
}

TEST_F(LeakoffCouplingTest, CumulativeLeakoffGrowsWithSqrtTime)
{
  const PetscReal c_l = 2.0e-4;
  const PetscReal t_open = 0.0;

  const PetscReal v1 = PetscFEHydrofrac::carterCumulativeLeakoff(c_l, 1.0, t_open);
  const PetscReal v4 = PetscFEHydrofrac::carterCumulativeLeakoff(c_l, 4.0, t_open);
  const PetscReal v9 = PetscFEHydrofrac::carterCumulativeLeakoff(c_l, 9.0, t_open);

  // V = 2*C_L*sqrt(t), so V(4)/V(1) = 2, V(9)/V(1) = 3
  EXPECT_NEAR(v4 / v1, 2.0, 1e-12);
  EXPECT_NEAR(v9 / v1, 3.0, 1e-12);

  // Absolute value check: V(1) = 2 * 2e-4 * 1 = 4e-4
  EXPECT_NEAR(v1, 4.0e-4, 1e-12);
}

TEST_F(LeakoffCouplingTest, TotalVolumeScalesWithArea)
{
  const PetscReal c_l = 1.0e-4;
  const PetscReal t = 100.0;
  const PetscReal t_open = 0.0;

  const PetscReal vol_1 = PetscFEHydrofrac::totalLeakoffVolume(c_l, t, t_open, 1.0);
  const PetscReal vol_10 = PetscFEHydrofrac::totalLeakoffVolume(c_l, t, t_open, 10.0);

  EXPECT_NEAR(vol_10 / vol_1, 10.0, 1e-10);
}

TEST_F(LeakoffCouplingTest, NoLeakoffBeforeExposure)
{
  const PetscReal c_l = 1.0e-4;

  EXPECT_EQ(PetscFEHydrofrac::carterLeakoffRate(c_l, 5.0, 10.0), 0.0);
  EXPECT_EQ(PetscFEHydrofrac::carterCumulativeLeakoff(c_l, 5.0, 10.0), 0.0);
  EXPECT_EQ(PetscFEHydrofrac::totalLeakoffVolume(c_l, 5.0, 10.0, 100.0), 0.0);
}

TEST_F(LeakoffCouplingTest, DelayedOpeningShiftsTimeline)
{
  const PetscReal c_l = 1.0e-4;
  const PetscReal t_open = 50.0;

  // At t=51, effective exposure = 1s
  const PetscReal q_51 = PetscFEHydrofrac::carterLeakoffRate(c_l, 51.0, t_open);
  // At t=1 with t_open=0, same exposure
  const PetscReal q_1 = PetscFEHydrofrac::carterLeakoffRate(c_l, 1.0, 0.0);

  EXPECT_NEAR(q_51, q_1, 1e-12);
}

TEST_F(LeakoffCouplingTest, LeakoffRateAtExposureTimeIsZero)
{
  // At t = t_open exactly, sqrt(0) -> infinite rate, but our function returns 0
  // because (time <= time_open)
  const PetscReal q = PetscFEHydrofrac::carterLeakoffRate(1.0e-4, 10.0, 10.0);
  EXPECT_EQ(q, 0.0);
}

TEST_F(LeakoffCouplingTest, InvalidCoefficientReturnsZero)
{
  EXPECT_EQ(PetscFEHydrofrac::carterLeakoffRate(0.0, 10.0, 0.0), 0.0);
  EXPECT_EQ(PetscFEHydrofrac::carterLeakoffRate(-1.0, 10.0, 0.0), 0.0);
  EXPECT_EQ(PetscFEHydrofrac::carterCumulativeLeakoff(0.0, 10.0, 0.0), 0.0);
  EXPECT_EQ(PetscFEHydrofrac::totalLeakoffVolume(1.0e-4, 10.0, 0.0, 0.0), 0.0);
}
