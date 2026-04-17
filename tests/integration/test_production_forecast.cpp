/**
 * @file test_production_forecast.cpp
 * @brief Phase 9 production forecasting integration tests
 *
 * Validates Arps decline curve analysis, cumulative production,
 * and fracture productivity index for post-frac production.
 */

#include <gtest/gtest.h>

#include <cmath>

#include "numerics/PetscFEHydrofrac.hpp"

using namespace FSRM;

class ProductionForecastTest : public ::testing::Test
{
};

TEST_F(ProductionForecastTest, ArpsDeclineRateDecreasesWithTime)
{
  const PetscReal q_i = 100.0;     // 100 bbl/day
  const PetscReal D_i = 0.1;       // 10%/year nominal decline
  const PetscReal b = 0.8;         // Hyperbolic

  const PetscReal q0 = PetscFEHydrofrac::arpsDeclineRate(q_i, D_i, b, 0.0);
  const PetscReal q1 = PetscFEHydrofrac::arpsDeclineRate(q_i, D_i, b, 365.0);
  const PetscReal q5 = PetscFEHydrofrac::arpsDeclineRate(q_i, D_i, b, 1825.0);

  EXPECT_NEAR(q0, q_i, 1e-10);
  EXPECT_LT(q1, q0);
  EXPECT_LT(q5, q1);
  EXPECT_GT(q5, 0.0);
}

TEST_F(ProductionForecastTest, ExponentialDeclineMatchesAnalytical)
{
  const PetscReal q_i = 200.0;
  const PetscReal D_i = 0.01;      // per day
  const PetscReal b = 0.0;         // Exponential
  const PetscReal t = 100.0;       // days

  const PetscReal q = PetscFEHydrofrac::arpsDeclineRate(q_i, D_i, b, t);
  const PetscReal expected = q_i * std::exp(-D_i * t);
  EXPECT_NEAR(q, expected, 1e-8);
}

TEST_F(ProductionForecastTest, CumulativeProductionIncreases)
{
  const PetscReal q_i = 100.0;
  const PetscReal D_i = 0.05;
  const PetscReal b = 0.8;

  const PetscReal np1 = PetscFEHydrofrac::arpsCumulativeProduction(q_i, D_i, b, 30.0);
  const PetscReal np2 = PetscFEHydrofrac::arpsCumulativeProduction(q_i, D_i, b, 365.0);
  const PetscReal np3 = PetscFEHydrofrac::arpsCumulativeProduction(q_i, D_i, b, 3650.0);

  EXPECT_GT(np1, 0.0);
  EXPECT_GT(np2, np1);
  EXPECT_GT(np3, np2);
}

TEST_F(ProductionForecastTest, ExponentialCumulativeMatchesAnalytical)
{
  const PetscReal q_i = 200.0;
  const PetscReal D_i = 0.01;
  const PetscReal b = 0.0;
  const PetscReal t = 100.0;

  const PetscReal np = PetscFEHydrofrac::arpsCumulativeProduction(q_i, D_i, b, t);
  const PetscReal expected = (q_i / D_i) * (1.0 - std::exp(-D_i * t));
  EXPECT_NEAR(np, expected, 1e-6);
}

TEST_F(ProductionForecastTest, HarmonicDeclineSpecialCase)
{
  const PetscReal q_i = 100.0;
  const PetscReal D_i = 0.02;
  const PetscReal b = 1.0;         // Harmonic
  const PetscReal t = 50.0;

  const PetscReal q = PetscFEHydrofrac::arpsDeclineRate(q_i, D_i, b, t);
  const PetscReal expected = q_i / (1.0 + D_i * t);
  EXPECT_NEAR(q, expected, 1e-8);
}

TEST_F(ProductionForecastTest, HarmonicCumulativeMatchesAnalytical)
{
  const PetscReal q_i = 100.0;
  const PetscReal D_i = 0.02;
  const PetscReal b = 1.0;
  const PetscReal t = 50.0;

  const PetscReal np = PetscFEHydrofrac::arpsCumulativeProduction(q_i, D_i, b, t);
  const PetscReal expected = (q_i / D_i) * std::log(1.0 + D_i * t);
  EXPECT_NEAR(np, expected, 1e-6);
}

TEST_F(ProductionForecastTest, ProductivityIndexPositive)
{
  // Propped fracture: k_f = 100 Darcy = 100e-12 m^2
  // w_f = 5mm = 0.005 m, mu = 1e-3 Pa.s
  // r_e = 100 m, r_w = 0.1 m
  const PetscReal J = PetscFEHydrofrac::fractureProductivityIndex(
      100.0e-12, 0.005, 1.0e-3, 100.0, 0.1);

  EXPECT_GT(J, 0.0);
}

TEST_F(ProductionForecastTest, ProductivityIndexScalesWithConductivity)
{
  const PetscReal mu = 1.0e-3;
  const PetscReal re = 100.0;
  const PetscReal rw = 0.1;

  const PetscReal J1 = PetscFEHydrofrac::fractureProductivityIndex(
      10.0e-12, 0.005, mu, re, rw);
  const PetscReal J2 = PetscFEHydrofrac::fractureProductivityIndex(
      20.0e-12, 0.005, mu, re, rw);

  EXPECT_NEAR(J2 / J1, 2.0, 1e-10);
}

TEST_F(ProductionForecastTest, InvalidInputsReturnZero)
{
  EXPECT_EQ(PetscFEHydrofrac::arpsDeclineRate(0.0, 0.1, 0.5, 10.0), 0.0);
  EXPECT_EQ(PetscFEHydrofrac::arpsDeclineRate(100.0, 0.0, 0.5, 10.0), 0.0);
  EXPECT_EQ(PetscFEHydrofrac::arpsCumulativeProduction(100.0, 0.1, 0.5, 0.0), 0.0);
  EXPECT_EQ(PetscFEHydrofrac::fractureProductivityIndex(0.0, 0.005, 1e-3, 100.0, 0.1), 0.0);
  EXPECT_EQ(PetscFEHydrofrac::fractureProductivityIndex(1e-12, 0.005, 1e-3, 100.0, 100.0), 0.0);
}
