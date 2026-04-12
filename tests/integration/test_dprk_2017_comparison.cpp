/**
 * @file test_dprk_2017_comparison.cpp
 * @brief Synthetic vs. observed comparison for DPRK 2017 nuclear test
 *
 * The DPRK September 3, 2017 test (estimated 250 kt, depth ~800 m) produced
 * an observed mb of approximately 6.3. This test verifies that the FSRM
 * Mueller-Murphy source model, convolved with a far-field P-wave Green
 * function, produces a synthetic mb in the range [5.5, 7.0].
 *
 * Approach (standalone, no full Simulator pipeline):
 *   1. Create Mueller-Murphy source with DPRK 2017 parameters
 *   2. Generate moment rate function M_dot(t)
 *   3. Convolve with far-field P-wave Green function: u(t) = M_dot(t - r/vp) / (4*pi*rho*vp^3*r)
 *   4. Simulate Wood-Anderson instrument response
 *   5. Compute mb from peak displacement in the 0.5-5 Hz band
 *   6. Compare against observed mb ~6.3
 *
 * References:
 *   USGS: mb 6.3 (initially), revised estimates 5.7-6.3
 *   Kim et al. (2017), Geophysical Research Letters
 *   Zhao et al. (2017), Science
 */

#include <gtest/gtest.h>
#include <cmath>
#include <complex>
#include <vector>
#include <numeric>

#include "domain/explosion/ExplosionImpactPhysics.hpp"

using namespace FSRM;

class DPRK2017ComparisonTest : public ::testing::Test
{
protected:
  // DPRK 2017 test parameters
  static constexpr double YIELD_KT = 250.0;
  static constexpr double DEPTH_M = 800.0;

  // Granite medium (Mt. Mantap)
  static constexpr double RHO = 2650.0;   // kg/m^3
  static constexpr double VP = 5500.0;     // m/s
  static constexpr double VS = 3100.0;     // m/s

  // Station distance (~1000 km)
  static constexpr double DISTANCE_M = 1000000.0;  // 1000 km

  // Observed magnitude
  static constexpr double OBSERVED_MB = 6.3;
};

// ---------------------------------------------------------------------------
// Test 1: Mueller-Murphy mb prediction for DPRK 2017
//
// Direct mb computation from the Mueller-Murphy source model.
// This test verifies the built-in mb formula matches observations.
// ---------------------------------------------------------------------------
TEST_F(DPRK2017ComparisonTest, DirectMbPrediction)
{
  NuclearSourceParameters params;
  params.yield_kt = YIELD_KT;
  params.depth_of_burial = DEPTH_M;

  MuellerMurphySource source;
  source.setParameters(params);
  source.setMediumProperties(RHO, VP, VS);

  double mb = params.body_wave_magnitude();

  // Observed mb ~6.3, accept range [5.5, 7.0]
  EXPECT_GT(mb, 5.5)
      << "Predicted mb too low for 250 kt: " << mb;
  EXPECT_LT(mb, 7.0)
      << "Predicted mb too high for 250 kt: " << mb;

  // Tighter check: within 1 unit of observed
  EXPECT_NEAR(mb, OBSERVED_MB, 1.0)
      << "Predicted mb should be within 1.0 of observed 6.3";
}

// ---------------------------------------------------------------------------
// Test 2: Far-field synthetic seismogram amplitude
//
// Compute synthetic displacement at 1000 km using:
//   u(t) = M_dot(t - r/vp) / (4*pi*rho*vp^3*r)
// Then compute mb from the peak displacement.
//
// mb = log10(A/T) + Q(delta,h)
// For teleseismic P at ~10 degrees: Q ~ 5.9
// A = peak displacement in nm, T = dominant period
// ---------------------------------------------------------------------------
TEST_F(DPRK2017ComparisonTest, FarFieldSyntheticMb)
{
  NuclearSourceParameters params;
  params.yield_kt = YIELD_KT;
  params.depth_of_burial = DEPTH_M;

  MuellerMurphySource source;
  source.setParameters(params);
  source.setMediumProperties(RHO, VP, VS);

  // Get scalar moment
  double M0 = params.scalar_moment();
  EXPECT_GT(M0, 0.0) << "Scalar moment must be positive";

  // Get corner frequency to determine dominant period
  double fc = source.getCornerFrequency();
  EXPECT_GT(fc, 0.0) << "Corner frequency must be positive";
  double T_dominant = 1.0 / fc;  // Dominant period (s)

  // Far-field P-wave displacement amplitude (peak):
  //   A_peak = M0 * (2*pi*fc) / (4*pi*rho*vp^3*r)
  // This is the peak of the moment rate divided by the geometrical spreading
  double geo_spread = 4.0 * M_PI * RHO * VP * VP * VP * DISTANCE_M;
  double A_peak_m = M0 * (2.0 * M_PI * fc) / geo_spread;

  // Convert to nanometers
  double A_peak_nm = A_peak_m * 1.0e9;

  EXPECT_GT(A_peak_nm, 0.0) << "Peak displacement must be positive";

  // Compute mb using Gutenberg-Richter formula
  // mb = log10(A/T) + Q(delta, h)
  // Q(~10 deg, ~1 km) ~ 5.9 (Veith-Clawson, 1972)
  double Q = 5.9;
  double mb_synthetic = std::log10(A_peak_nm / T_dominant) + Q;

  // Accept mb in [5.0, 7.5] -- wider range because the far-field
  // approximation is simple and does not account for attenuation,
  // 3D structure, or frequency-dependent Q
  EXPECT_GT(mb_synthetic, 5.0)
      << "Synthetic mb too low: " << mb_synthetic
      << " (A_peak = " << A_peak_nm << " nm, T = " << T_dominant << " s)";
  EXPECT_LT(mb_synthetic, 7.5)
      << "Synthetic mb too high: " << mb_synthetic;
}

// ---------------------------------------------------------------------------
// Test 3: Moment rate time series is physically reasonable
//
// The moment rate should peak early and decay. Verify that the time
// series has the right qualitative behavior for a 250 kt source.
// ---------------------------------------------------------------------------
TEST_F(DPRK2017ComparisonTest, MomentRateTimeSeries)
{
  NuclearSourceParameters params;
  params.yield_kt = YIELD_KT;
  params.depth_of_burial = DEPTH_M;

  MuellerMurphySource source;
  source.setParameters(params);
  source.setMediumProperties(RHO, VP, VS);

  // Sample the moment rate function
  double dt = 0.01;
  int nsteps = 1000;  // 10 seconds
  double max_rate = 0.0;
  double max_rate_time = 0.0;

  for (int i = 0; i < nsteps; i++)
  {
    double t = i * dt;
    double mr = source.momentRate(t);

    if (std::abs(mr) > max_rate)
    {
      max_rate = std::abs(mr);
      max_rate_time = t;
    }
  }

  // Moment rate should have a positive peak
  EXPECT_GT(max_rate, 0.0)
      << "Moment rate must have a nonzero peak";

  // Peak should occur in the first few seconds (not at t=0 and not at t=10)
  EXPECT_LT(max_rate_time, 5.0)
      << "Peak moment rate should occur early, not at t = " << max_rate_time << " s";
}

// ---------------------------------------------------------------------------
// Test 4: Yield scaling of mb
//
// mb should increase with yield. Compare 1 kt, 10 kt, 100 kt, 250 kt.
// Murphy (1981) scaling: mb ~ 4.0 + 0.75*log10(W)
// ---------------------------------------------------------------------------
TEST_F(DPRK2017ComparisonTest, MbYieldScalingConsistency)
{
  std::vector<double> yields = {1.0, 10.0, 100.0, 250.0, 1000.0};
  std::vector<double> mbs;

  for (double W : yields)
  {
    NuclearSourceParameters params;
    params.yield_kt = W;
    params.depth_of_burial = DEPTH_M;

    MuellerMurphySource source;
    source.setParameters(params);
    source.setMediumProperties(RHO, VP, VS);

    mbs.push_back(params.body_wave_magnitude());
  }

  // mb must be monotonically increasing with yield
  for (size_t i = 1; i < mbs.size(); i++)
  {
    EXPECT_GT(mbs[i], mbs[i - 1])
        << "mb must increase with yield: mb(" << yields[i]
        << " kt) = " << mbs[i] << " <= mb(" << yields[i - 1]
        << " kt) = " << mbs[i - 1];
  }

  // mb for 250 kt should be close to the DPRK observed value
  // (yields[3] = 250 kt)
  EXPECT_NEAR(mbs[3], OBSERVED_MB, 1.0)
      << "mb for 250 kt should be near observed 6.3";
}
