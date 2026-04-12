/**
 * @file test_induced_seismicity.cpp
 * @brief Phase 7 induced seismicity integration tests
 *
 * Validates moment tensor computation from fracture events
 * and microseismic magnitude estimation.
 */

#include <gtest/gtest.h>

#include <cmath>

#include "numerics/PetscFEHydrofrac.hpp"

using namespace FSRM;

class InducedSeismicityTest : public ::testing::Test
{
};

TEST_F(InducedSeismicityTest, ScalarMomentFromFractureEvent)
{
  // Typical hydrofrac microseismic event:
  // mu = 30 GPa, slip = 1 mm, area = 10 m^2
  const PetscReal mu = 30.0e9;
  const PetscReal slip = 1.0e-3;
  const PetscReal area = 10.0;

  const PetscReal m0 = PetscFEHydrofrac::scalarMoment(mu, slip, area);
  EXPECT_NEAR(m0, 3.0e8, 1e-6);
}

TEST_F(InducedSeismicityTest, MomentMagnitudeInMicroseismicRange)
{
  // M0 = 3e8 N.m -> Mw = (2/3)*log10(3e8) - 6.07
  const PetscReal m0 = 3.0e8;
  const PetscReal mw = PetscFEHydrofrac::momentMagnitude(m0);

  // Expected: (2/3)*8.477 - 6.07 = 5.651 - 6.07 = -0.42
  EXPECT_NEAR(mw, -0.42, 0.1);
  EXPECT_EQ(PetscFEHydrofrac::isMicroseismicRange(mw), PETSC_TRUE);
}

TEST_F(InducedSeismicityTest, LargerEventsProduceLargerMagnitudes)
{
  const PetscReal mu = 30.0e9;

  // Small event: 0.1 mm slip, 1 m^2
  const PetscReal m0_small = PetscFEHydrofrac::scalarMoment(mu, 1.0e-4, 1.0);
  const PetscReal mw_small = PetscFEHydrofrac::momentMagnitude(m0_small);

  // Large event: 10 mm slip, 100 m^2
  const PetscReal m0_large = PetscFEHydrofrac::scalarMoment(mu, 1.0e-2, 100.0);
  const PetscReal mw_large = PetscFEHydrofrac::momentMagnitude(m0_large);

  EXPECT_GT(m0_large, m0_small);
  EXPECT_GT(mw_large, mw_small);
}

TEST_F(InducedSeismicityTest, MicroseismicRangeBoundary)
{
  // Mw = -3.0 is outside (boundary)
  EXPECT_EQ(PetscFEHydrofrac::isMicroseismicRange(-3.0), PETSC_FALSE);
  // Mw = -2.9 is inside
  EXPECT_EQ(PetscFEHydrofrac::isMicroseismicRange(-2.9), PETSC_TRUE);
  // Mw = 0.9 is inside
  EXPECT_EQ(PetscFEHydrofrac::isMicroseismicRange(0.9), PETSC_TRUE);
  // Mw = 1.0 is outside (boundary)
  EXPECT_EQ(PetscFEHydrofrac::isMicroseismicRange(1.0), PETSC_FALSE);
}

TEST_F(InducedSeismicityTest, TypicalHydrofracEventMagnitudes)
{
  const PetscReal mu = 30.0e9;

  // Sweep of typical hydrofrac events
  struct EventParams
  {
    PetscReal slip;
    PetscReal area;
    PetscReal expected_mw_min;
    PetscReal expected_mw_max;
  };

  // Common values from Davies et al. (2013) and Warpinski (2009)
  EventParams events[] = {
      {1e-4, 1.0, -3.0, -1.0},    // Very small, single cell
      {1e-3, 10.0, -1.5, 0.5},    // Typical fracture tip event
      {5e-3, 50.0, 0.0, 2.0},     // Large induced event
  };

  for (const auto& ev : events)
  {
    const PetscReal m0 = PetscFEHydrofrac::scalarMoment(mu, ev.slip, ev.area);
    const PetscReal mw = PetscFEHydrofrac::momentMagnitude(m0);
    EXPECT_GE(mw, ev.expected_mw_min) << "slip=" << ev.slip << " area=" << ev.area;
    EXPECT_LE(mw, ev.expected_mw_max) << "slip=" << ev.slip << " area=" << ev.area;
  }
}

TEST_F(InducedSeismicityTest, InvalidMomentReturnsNegativeMagnitude)
{
  EXPECT_LT(PetscFEHydrofrac::momentMagnitude(0.0), -90.0);
  EXPECT_LT(PetscFEHydrofrac::momentMagnitude(-1.0), -90.0);
  EXPECT_EQ(PetscFEHydrofrac::scalarMoment(0.0, 1.0, 1.0), 0.0);
  EXPECT_EQ(PetscFEHydrofrac::scalarMoment(1.0, 0.0, 1.0), 0.0);
  EXPECT_EQ(PetscFEHydrofrac::scalarMoment(1.0, 1.0, 0.0), 0.0);
}
