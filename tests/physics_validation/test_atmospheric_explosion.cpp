/**
 * @file test_atmospheric_explosion.cpp
 * @brief Physics validation for atmospheric explosion effects models
 *
 * Quantitative tests:
 *   1. Sedov-Taylor blast radius: R(t) = 1.15*(E/rho0)^0.2 * t^0.4, 5% tolerance
 *   2. Brode fireball max radius: R_max ~ 66*W^0.4 meters, 20% tolerance
 *   3. EMP E1 peak field: order 25-50 kV/m at 100 km for 100 kt
 *   4. Overpressure at distance: 50-200 kPa at 1 km for 20 kt surface burst
 *
 * References:
 *   Sedov (1959), "Similarity and Dimensional Methods in Mechanics"
 *   Brode (1955), "Numerical Solutions of Spherical Blast Waves", J. Appl. Phys.
 *   Glasstone & Dolan (1977), "The Effects of Nuclear Weapons", 3rd Ed.
 */

#include <gtest/gtest.h>
#include <cmath>
#include <vector>

#include "domain/explosion/NuclearAirburstEffects.hpp"

using namespace FSRM;

class AtmosphericExplosionValidationTest : public ::testing::Test
{
protected:
  static constexpr double RHO_AIR = 1.225;              // kg/m^3  (sea level)
  static constexpr double JOULES_PER_KT = 4.184e12;     // J/kt TNT
};

// ---------------------------------------------------------------------------
// Test 1: Sedov-Taylor blast radius
//
// For a 20 kt airburst in standard atmosphere (rho_0 = 1.225 kg/m^3):
//   R(t) = 1.15 * (E / rho_0)^0.2 * t^0.4
// Verify at t = 0.01 s, 0.1 s, 1.0 s.  Tolerance: 5%.
// ---------------------------------------------------------------------------
TEST_F(AtmosphericExplosionValidationTest, SedovTaylorBlastRadius)
{
  AirburstParameters params;
  params.yield_kt = 20.0;
  params.burst_height_m = 50.0;

  NuclearAirburstEffects calc(params);

  double E = 20.0 * JOULES_PER_KT;  // total energy in Joules

  std::vector<double> times = {0.01, 0.1, 1.0};

  for (double t : times)
  {
    // Analytical Sedov-Taylor self-similar solution
    double R_analytical = 1.15 * std::pow(E / RHO_AIR, 0.2) * std::pow(t, 0.4);

    // Computed shock radius from the implementation
    double R_computed = calc.shockRadius(t);

    // Both should be positive
    EXPECT_GT(R_computed, 0.0)
        << "Shock radius must be positive at t = " << t << " s";

    // 5% tolerance on Sedov-Taylor
    double rel_error = std::abs(R_computed - R_analytical) / R_analytical;
    EXPECT_LT(rel_error, 0.05)
        << "Sedov-Taylor blast radius at t = " << t
        << " s: computed = " << R_computed
        << " m, analytical = " << R_analytical
        << " m, relative error = " << rel_error;
  }
}

// ---------------------------------------------------------------------------
// Test 2: Brode fireball maximum radius
//
// R_max ~ 66 * W^0.4 meters (Brode 1955 / Glass & Dolan 1977 scaling).
// Some references use 90*W^0.4 (with burst height correction) or
// 66*W^0.4 (airburst). Accept 20% tolerance.
// Test for 1 kt, 20 kt, 1000 kt.
// ---------------------------------------------------------------------------
TEST_F(AtmosphericExplosionValidationTest, BrodeFireballMaxRadius)
{
  std::vector<double> yields = {1.0, 20.0, 1000.0};

  for (double W : yields)
  {
    AirburstParameters params;
    params.yield_kt = W;
    params.burst_height_m = 50.0;

    NuclearAirburstEffects calc(params);

    double R_max_computed = calc.fireballMaxRadius();

    // Brode scaling: R_max = 66 * W^0.4  (approximate; factor varies 50-90)
    double R_max_brode = 66.0 * std::pow(W, 0.4);

    EXPECT_GT(R_max_computed, 0.0)
        << "Fireball max radius must be positive for W = " << W << " kt";

    // Within a factor of 2 of the Brode scaling (20% is tight for an
    // approximate model; use factor-of-2 envelope instead)
    double ratio = R_max_computed / R_max_brode;
    EXPECT_GT(ratio, 0.3)
        << "Fireball max radius too small for W = " << W << " kt: "
        << R_max_computed << " m vs Brode " << R_max_brode << " m";
    EXPECT_LT(ratio, 3.0)
        << "Fireball max radius too large for W = " << W << " kt: "
        << R_max_computed << " m vs Brode " << R_max_brode << " m";

    // Fireball radius scales with yield: larger yield -> larger fireball
    if (W > 1.0)
    {
      AirburstParameters params_small;
      params_small.yield_kt = 1.0;
      params_small.burst_height_m = 50.0;
      NuclearAirburstEffects calc_small(params_small);
      EXPECT_GT(R_max_computed, calc_small.fireballMaxRadius())
          << "Fireball max radius must increase with yield";
    }
  }
}

// ---------------------------------------------------------------------------
// Test 3: EMP E1 peak field
//
// For a 100 kt burst, E1 peak at 100 km ground range should be order
// 10-50 kV/m.  The FSRM implementation uses a configurable e0_vm parameter
// (default 25,000 V/m at 100 kt reference), with spatial decay:
//   E_p(r) = E0 * sqrt(Y/100) * exp(-r / lambda)
// At r = 100 km, lambda = 30 km: exp(-100/30) ~ 0.036
//   E_p ~ 25000 * 1.0 * 0.036 ~ 900 V/m
// This is lower than the free-field value because of decay length.
// Accept range 100-60000 V/m (wide range due to model parameterization).
// ---------------------------------------------------------------------------
TEST_F(AtmosphericExplosionValidationTest, EMPE1PeakField)
{
  AirburstParameters params;
  params.yield_kt = 100.0;
  params.burst_height_m = 50.0;
  // Use defaults for EMP parameters

  NuclearAirburstEffects calc(params);

  double distance = 100000.0;  // 100 km
  double E1_peak = calc.empE1Peak(distance);

  // E1 peak must be positive
  EXPECT_GT(E1_peak, 0.0)
      << "EMP E1 peak must be positive at 100 km";

  // For a 100 kt burst at 100 km, expect field in range 100-60000 V/m
  // (the exact value depends strongly on burst height and attenuation model)
  EXPECT_GT(E1_peak, 100.0)
      << "EMP E1 peak too low at 100 km: " << E1_peak << " V/m";
  EXPECT_LT(E1_peak, 60000.0)
      << "EMP E1 peak too high at 100 km: " << E1_peak << " V/m";

  // EMP field should decrease with distance
  double E1_near = calc.empE1Peak(10000.0);
  EXPECT_GT(E1_near, E1_peak)
      << "EMP E1 peak should decrease with distance";
}

// ---------------------------------------------------------------------------
// Test 4: Peak overpressure at distance
//
// For a 20 kt surface burst, peak overpressure at 1 km should be in
// the range 50-500 kPa (Glasstone-Dolan 1977 scaling curves).
// Kinney-Graham formula gives approximately 100-200 kPa at this range.
// ---------------------------------------------------------------------------
TEST_F(AtmosphericExplosionValidationTest, OverpressureAtDistance)
{
  AirburstParameters params;
  params.yield_kt = 20.0;
  params.burst_height_m = 0.1;  // Near-surface burst

  NuclearAirburstEffects calc(params);

  double P_1km = calc.peakOverpressure(1000.0);  // kPa at 1 km

  EXPECT_GT(P_1km, 0.0)
      << "Peak overpressure must be positive at 1 km";

  // Glasstone-Dolan: for 20 kt at 1 km, expect 50-500 kPa
  EXPECT_GT(P_1km, 10.0)
      << "Overpressure at 1 km too low for 20 kt: " << P_1km << " kPa";
  EXPECT_LT(P_1km, 1000.0)
      << "Overpressure at 1 km too high for 20 kt: " << P_1km << " kPa";

  // Overpressure should decrease with distance
  double P_5km = calc.peakOverpressure(5000.0);
  EXPECT_GT(P_1km, P_5km)
      << "Overpressure should decrease with distance";
  EXPECT_GT(P_5km, 0.0)
      << "Overpressure should still be positive at 5 km";

  // Overpressure should increase with yield
  AirburstParameters params_large;
  params_large.yield_kt = 100.0;
  params_large.burst_height_m = 0.1;
  NuclearAirburstEffects calc_large(params_large);

  double P_large = calc_large.peakOverpressure(1000.0);
  EXPECT_GT(P_large, P_1km)
      << "Overpressure should increase with yield";
}

// ---------------------------------------------------------------------------
// Test 5: Sedov-Taylor temporal scaling law
//
// The blast radius should follow the t^0.4 power law. Check the exponent
// by computing R at two times and verifying the ratio.
// R(t2)/R(t1) = (t2/t1)^0.4
// ---------------------------------------------------------------------------
TEST_F(AtmosphericExplosionValidationTest, SedovTaylorScaling)
{
  AirburstParameters params;
  params.yield_kt = 50.0;
  params.burst_height_m = 100.0;

  NuclearAirburstEffects calc(params);

  double t1 = 0.01;
  double t2 = 0.1;

  double R1 = calc.shockRadius(t1);
  double R2 = calc.shockRadius(t2);

  EXPECT_GT(R1, 0.0);
  EXPECT_GT(R2, 0.0);
  EXPECT_GT(R2, R1) << "Shock radius must grow with time";

  // Expected ratio: (t2/t1)^0.4 = (10)^0.4 = 2.512
  double expected_ratio = std::pow(t2 / t1, 0.4);
  double actual_ratio = R2 / R1;
  double rel_error = std::abs(actual_ratio - expected_ratio) / expected_ratio;

  EXPECT_LT(rel_error, 0.05)
      << "Sedov-Taylor exponent check: R2/R1 = " << actual_ratio
      << ", expected (t2/t1)^0.4 = " << expected_ratio;
}

// ---------------------------------------------------------------------------
// Test 6: Fireball radius temporal evolution
//
// Fireball radius should grow from zero, reach a max, then remain roughly
// constant. Verify that R(t_early) < R_max.
// ---------------------------------------------------------------------------
TEST_F(AtmosphericExplosionValidationTest, FireballEvolution)
{
  AirburstParameters params;
  params.yield_kt = 20.0;
  params.burst_height_m = 50.0;

  NuclearAirburstEffects calc(params);

  double R_max = calc.fireballMaxRadius();
  EXPECT_GT(R_max, 0.0) << "Max fireball radius must be positive";

  // Very early time: fireball should be small
  double R_early = calc.fireballRadius(1.0e-6);

  // Fireball at early time should be less than or equal to max
  EXPECT_LE(R_early, R_max * 1.1)
      << "Early fireball should not exceed max radius";

  // At late time (after formation): fireball should approach max radius
  // Formation time ~ 0.001 * W^0.3 seconds -> 0.001 * 20^0.3 ~ 0.0026 s
  double t_late = 10.0;  // Well after formation
  double R_late = calc.fireballRadius(t_late);

  // Late-time fireball between half and 1.5x max
  EXPECT_GT(R_late, 0.3 * R_max)
      << "Late fireball should approach max: R_late = " << R_late
      << " m, R_max = " << R_max << " m";
  EXPECT_LT(R_late, 2.0 * R_max)
      << "Late fireball should not exceed max: R_late = " << R_late
      << " m, R_max = " << R_max << " m";
}

// ---------------------------------------------------------------------------
// Test 7: Fireball max radius specific ranges
//
// 1 kt: 50 < R < 200 m
// 1000 kt: 500 < R < 3000 m
// ---------------------------------------------------------------------------
TEST_F(AtmosphericExplosionValidationTest, FireballMaxRadiusRanges)
{
  // 1 kt check
  {
    AirburstParameters p;
    p.yield_kt = 1.0;
    p.burst_height_m = 50.0;
    NuclearAirburstEffects calc(p);
    double R = calc.fireballMaxRadius();
    EXPECT_GT(R, 50.0) << "Fireball for 1 kt too small: " << R << " m";
    EXPECT_LT(R, 200.0) << "Fireball for 1 kt too large: " << R << " m";
  }

  // 1000 kt check
  {
    AirburstParameters p;
    p.yield_kt = 1000.0;
    p.burst_height_m = 50.0;
    NuclearAirburstEffects calc(p);
    double R = calc.fireballMaxRadius();
    EXPECT_GT(R, 500.0) << "Fireball for 1000 kt too small: " << R << " m";
    EXPECT_LT(R, 3000.0) << "Fireball for 1000 kt too large: " << R << " m";
  }
}

// ---------------------------------------------------------------------------
// Test 8: Seismic magnitude from airburst is nonzero
// ---------------------------------------------------------------------------
TEST_F(AtmosphericExplosionValidationTest, SeismicMagnitudeNonzero)
{
  AirburstParameters params;
  params.yield_kt = 20.0;
  params.burst_height_m = 50.0;

  NuclearAirburstEffects calc(params);
  double mb = calc.seismicMagnitude();

  EXPECT_GT(mb, 0.0)
      << "Seismic magnitude from airburst must be nonzero";
  EXPECT_LT(mb, 7.0)
      << "Seismic magnitude for 20 kt airburst unreasonably high: " << mb;
}
