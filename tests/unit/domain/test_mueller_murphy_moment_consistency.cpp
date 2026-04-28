/**
 * @file test_mueller_murphy_moment_consistency.cpp
 * @brief Unit test (M1.2 partial): Mueller-Murphy scalar-moment derived diagnostic
 *
 * Pure analytical test, no FEM, no mesh.  Configures a MuellerMurphySource
 * with the Sedan 1962 parameters (104 kt, 194 m depth of burial, alluvium
 * medium properties from models/sedan_1962.vel) and checks two analytical
 * consistency conditions:
 *
 *   1. Numerical integration of momentRate(t) from 0 to ~10/fc reproduces
 *      params.scalar_moment() within 1 percent.
 *   2. params.body_wave_magnitude() falls inside [4.0, 5.5] for the Sedan
 *      yield (USGS reports mb ~ 4.75 for Sedan).
 *
 * Reference: docs/NEM_ROADMAP.md M1.2.  This test does not exercise the
 * source-to-FEM coupling path; that is M1.3.
 */

#include <gtest/gtest.h>

#include <cmath>

#include "domain/explosion/ExplosionImpactPhysics.hpp"

using namespace FSRM;

namespace
{

// Sedan 1962 working point: alluvium medium (top layer of
// models/sedan_1962.vel: vp = 1.5 km/s, vs = 0.7 km/s, rho = 1900 kg/m^3).
constexpr double kSedanYieldKt        = 104.0;
constexpr double kSedanDepthOfBurial  = 194.0;
constexpr double kAlluviumDensity     = 1900.0;
constexpr double kAlluviumVp          = 1500.0;
constexpr double kAlluviumVs          =  700.0;

}  // namespace

class MuellerMurphyMomentConsistencyTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    params_.yield_kt        = kSedanYieldKt;
    params_.depth_of_burial = kSedanDepthOfBurial;
    source_.setParameters(params_);
    source_.setMediumProperties(kAlluviumDensity, kAlluviumVp, kAlluviumVs);
  }

  NuclearSourceParameters params_;
  MuellerMurphySource     source_;
};

// Test 1: numerical integration of momentRate(t) recovers the analytical
// scalar moment.  The Mueller-Murphy moment-rate function decays well
// within ~10/fc from the origin (the Physics.MuellerMurphy validation
// suite already verifies that |momentRate(10 * rise_time)| < 1% of peak).
//
// Tolerance: 1 percent, per docs/NEM_ROADMAP.md M1.2.
TEST_F(MuellerMurphyMomentConsistencyTest, NumericalMomentIntegralMatchesScalarMoment)
{
  const double M0_analytical = params_.scalar_moment();
  ASSERT_GT(M0_analytical, 0.0) << "Scalar moment must be positive";
  ASSERT_TRUE(std::isfinite(M0_analytical));

  const double fc = source_.getCornerFrequency();
  ASSERT_GT(fc, 0.0) << "Corner frequency must be positive";

  // Integrate from 0 to t_end = 10 / fc using composite Simpson's rule.
  // Use enough steps that the time step is much smaller than the rise
  // time (~0.55 / fc): n_steps = 8000 gives dt = (10 / fc) / 8000 which
  // is ~ 7e-4 / fc, well below any feature in the moment rate.
  const double t_end   = 10.0 / fc;
  const int    n_steps = 8000;            // even, required for Simpson's
  const double dt      = t_end / static_cast<double>(n_steps);

  double integral = 0.0;
  for (int i = 0; i <= n_steps; ++i)
  {
    const double t  = static_cast<double>(i) * dt;
    const double mr = source_.momentRate(t);
    double w = 0.0;
    if (i == 0 || i == n_steps)
    {
      w = 1.0;
    }
    else if (i % 2 == 1)
    {
      w = 4.0;
    }
    else
    {
      w = 2.0;
    }
    integral += w * mr;
  }
  integral *= dt / 3.0;

  const double rel_err = std::abs(integral - M0_analytical) / M0_analytical;

  EXPECT_LT(rel_err, 0.01)
      << "Integral of momentRate(t) from 0 to 10/fc = " << integral
      << " N*m, scalar_moment() = " << M0_analytical << " N*m, "
      << "relative error = " << rel_err
      << " (expected < 0.01).  Either the Mueller-Murphy moment-rate "
         "shape and the closed-form M0 are inconsistent, or the rise-"
         "time normalization is wrong.";
}

// Test 2: body-wave magnitude for 104 kt at 194 m in alluvium is in a
// physically defensible range.  USGS catalogs Sedan at mb ~ 4.75; the
// roadmap accepts any value in [4.0, 5.5] (generous but bounded so the
// test catches a sign error or a unit error in scalar_moment()).
TEST_F(MuellerMurphyMomentConsistencyTest, BodyWaveMagnitudeInSedanRange)
{
  const double mb = params_.body_wave_magnitude();
  EXPECT_TRUE(std::isfinite(mb)) << "mb must be finite";
  EXPECT_GE(mb, 4.0) << "mb = " << mb << " is below the Sedan range [4.0, 5.5]";
  EXPECT_LE(mb, 5.5) << "mb = " << mb << " is above the Sedan range [4.0, 5.5]";
}
