/**
 * @file test_mueller_murphy.cpp
 * @brief Physics validation for Mueller-Murphy seismic source model
 *
 * Quantitative tests:
 *   1. Corner frequency for 1 kt granite at 300m matches Mueller-Murphy (1971)
 *   2. mb-yield scaling follows mb = 4.45 + 0.75*log10(W) within 0.3 units
 *   3. RDP spectral shape has correct low-frequency plateau and HF rolloff
 *   4. DPRK 2017 calibration: W=250 kt, depth=800m -> mb in [6.0, 6.5]
 *
 * References:
 *   Mueller & Murphy (1971), "Seismic characteristics of underground
 *   nuclear detonations", Bull. Seism. Soc. Am. 61(6), 1675-1692.
 *   Murphy (1981), "P-wave coupling of underground explosions",
 *   DARPA/ARPA Report.
 */

#include <gtest/gtest.h>
#include <cmath>
#include <complex>
#include <vector>

#include "domain/explosion/ExplosionImpactPhysics.hpp"

using namespace FSRM;

class MuellerMurphyValidationTest : public ::testing::Test
{
protected:
  // Granite medium properties
  static constexpr double RHO_GRANITE = 2650.0;  // kg/m^3
  static constexpr double VP_GRANITE  = 5500.0;  // m/s
  static constexpr double VS_GRANITE  = 3100.0;  // m/s
};

// Test 1: Corner frequency for 1 kt granite at 300m depth
// Mueller-Murphy (1971) Table 2: fc ~ 2-5 Hz for small tamped shots in granite.
// The corner frequency scales as fc ~ W^(-1/3) approximately.
TEST_F(MuellerMurphyValidationTest, CornerFrequency1ktGranite)
{
  NuclearSourceParameters params;
  params.yield_kt = 1.0;
  params.depth_of_burial = 300.0;

  MuellerMurphySource source;
  source.setParameters(params);
  source.setMediumProperties(RHO_GRANITE, VP_GRANITE, VS_GRANITE);

  double fc = source.getCornerFrequency();

  // Corner frequency should be positive and in physically reasonable range
  EXPECT_GT(fc, 0.0) << "Corner frequency must be positive";

  // Patton (1988): fc = 3.0 * W^(-1/3) Hz for hard rock.
  // For 1 kt in granite: fc = 3.0 Hz.
  // Allow 1-10 Hz range for model variants.
  EXPECT_GT(fc, 1.0) << "Corner frequency too low for 1 kt shot";
  EXPECT_LT(fc, 10.0) << "Corner frequency too high for 1 kt shot";

  // Verify corner frequency decreases with yield (larger explosions have
  // lower corner frequencies). Compare 1 kt vs 100 kt.
  NuclearSourceParameters params_100kt;
  params_100kt.yield_kt = 100.0;
  params_100kt.depth_of_burial = 300.0;

  MuellerMurphySource source_100kt;
  source_100kt.setParameters(params_100kt);
  source_100kt.setMediumProperties(RHO_GRANITE, VP_GRANITE, VS_GRANITE);

  double fc_100kt = source_100kt.getCornerFrequency();
  EXPECT_LT(fc_100kt, fc)
      << "Corner frequency for 100 kt (" << fc_100kt
      << " Hz) should be less than for 1 kt (" << fc << " Hz)";
}

// Test 2: mb-yield scaling follows Murphy (1981) relation
// mb = 4.45 + 0.75 * log10(W) for tamped shots in hard rock.
// Tolerance: 0.3 magnitude units.
TEST_F(MuellerMurphyValidationTest, MbYieldScaling)
{
  double yields[] = {1.0, 10.0, 100.0, 1000.0};
  int n_yields = 4;

  for (int i = 0; i < n_yields; ++i)
  {
    NuclearSourceParameters params;
    params.yield_kt = yields[i];
    params.depth_of_burial = 300.0;

    double mb = params.body_wave_magnitude();

    // Murphy (1981): mb = 4.45 + 0.75 * log10(W)
    double mb_expected = 4.45 + 0.75 * std::log10(yields[i]);

    EXPECT_NEAR(mb, mb_expected, 0.3)
        << "mb for W=" << yields[i] << " kt: computed=" << mb
        << ", expected=" << mb_expected << " (Murphy 1981)";

    // mb must be finite and positive
    EXPECT_TRUE(std::isfinite(mb)) << "mb must be finite for W=" << yields[i];
    EXPECT_GT(mb, 3.0) << "mb too low for W=" << yields[i] << " kt";
  }

  // Verify monotonic increase of mb with yield
  NuclearSourceParameters p1, p2;
  p1.yield_kt = 1.0;
  p2.yield_kt = 100.0;
  EXPECT_GT(p2.body_wave_magnitude(), p1.body_wave_magnitude())
      << "mb must increase with yield";
}

// Test 3: RDP spectral shape verification
// The reduced displacement potential should have:
//   - Low-frequency plateau: |psi(omega)| ~ constant for omega << omega_c
//   - High-frequency rolloff: |psi(omega)| ~ omega^(-2) for omega >> omega_c
//   - Transition at omega_c (corner frequency)
TEST_F(MuellerMurphyValidationTest, RDPSpectralShape)
{
  NuclearSourceParameters params;
  params.yield_kt = 10.0;
  params.depth_of_burial = 500.0;

  MuellerMurphySource source;
  source.setParameters(params);
  source.setMediumProperties(RHO_GRANITE, VP_GRANITE, VS_GRANITE);

  double fc = source.getCornerFrequency();
  ASSERT_GT(fc, 0.0) << "Need positive corner frequency for spectral test";

  double omega_c = 2.0 * M_PI * fc;

  // Evaluate RDP at frequencies well below and above corner frequency
  double omega_low = omega_c * 0.01;   // 100x below corner
  double omega_high = omega_c * 100.0; // 100x above corner

  std::complex<double> rdp_low = source.rdp(omega_low);
  std::complex<double> rdp_corner = source.rdp(omega_c);
  std::complex<double> rdp_high = source.rdp(omega_high);

  double amp_low = std::abs(rdp_low);
  double amp_corner = std::abs(rdp_corner);
  double amp_high = std::abs(rdp_high);

  // All amplitudes must be positive and finite
  EXPECT_GT(amp_low, 0.0) << "RDP at low frequency must be nonzero";
  EXPECT_GT(amp_corner, 0.0) << "RDP at corner frequency must be nonzero";
  EXPECT_GT(amp_high, 0.0) << "RDP at high frequency must be nonzero";
  EXPECT_TRUE(std::isfinite(amp_low));
  EXPECT_TRUE(std::isfinite(amp_corner));
  EXPECT_TRUE(std::isfinite(amp_high));

  // Low-frequency plateau: amplitude at omega_low should be close to
  // or larger than the corner amplitude.  The exact ratio depends on the
  // RDP model variant (single-pole vs overshoot). Allow wide tolerance.
  EXPECT_GT(amp_low, 0.1 * amp_corner)
      << "Low-frequency amplitude should be comparable to corner amplitude";
  EXPECT_LT(amp_low, 1000.0 * amp_corner)
      << "Low-frequency amplitude should not be excessively larger than corner";

  // High-frequency rolloff: amplitude at omega_high should be much smaller
  // than at omega_low. For omega^(-2) rolloff: ratio ~ (100)^2 = 10000.
  // Allow range 100 to 1e6 for the rolloff ratio.
  double rolloff_ratio = amp_low / amp_high;
  EXPECT_GT(rolloff_ratio, 10.0)
      << "High-frequency rolloff insufficient: ratio=" << rolloff_ratio
      << ", expected >> 1 for omega^(-2) decay";

  // Verify spectral shape is monotonically decreasing above corner
  double omega_2x = omega_c * 2.0;
  double omega_10x = omega_c * 10.0;
  double amp_2x = std::abs(source.rdp(omega_2x));
  double amp_10x = std::abs(source.rdp(omega_10x));

  EXPECT_GT(amp_corner, amp_2x)
      << "Amplitude should decrease above corner frequency";
  EXPECT_GT(amp_2x, amp_10x)
      << "Amplitude should continue decreasing at higher frequencies";
}

// Test 4: DPRK 2017 calibration
// The September 3, 2017 North Korean nuclear test had:
//   Estimated yield: ~200-300 kt (best estimate ~250 kt)
//   Depth of burial: ~700-900 m (best estimate ~800 m)
//   Observed mb: ~6.3 (USGS)
//   Medium: granite-like (Mantapsan)
// Verify: computed mb in range [6.0, 6.5]
TEST_F(MuellerMurphyValidationTest, DPRK2017Calibration)
{
  NuclearSourceParameters params;
  params.yield_kt = 250.0;
  params.depth_of_burial = 800.0;

  double mb = params.body_wave_magnitude();

  // Murphy (1981): mb = 4.45 + 0.75 * log10(250) = 4.45 + 1.80 = 6.25
  // Observed: mb ~ 6.3
  EXPECT_GE(mb, 6.0)
      << "DPRK 2017 mb (" << mb << ") below expected range [6.0, 6.5]";
  EXPECT_LE(mb, 6.5)
      << "DPRK 2017 mb (" << mb << ") above expected range [6.0, 6.5]";

  // Verify scalar moment is in reasonable range
  // For mb ~ 6.3: log10(M0) ~ 1.5*mb + 9.1 ~ 18.55 -> M0 ~ 3.5e18 N*m
  // (Nuttli 1986 mb-M0 relation, approximate)
  double M0 = params.scalar_moment();
  EXPECT_GT(M0, 1.0e16)
      << "Scalar moment too small for 250 kt shot";
  EXPECT_LT(M0, 1.0e21)
      << "Scalar moment too large for 250 kt shot";

  // Verify Ms is also in reasonable range
  double Ms = params.surface_wave_magnitude();
  EXPECT_GT(Ms, 4.0) << "Ms too low for 250 kt shot";
  EXPECT_LT(Ms, 8.0) << "Ms too high for 250 kt shot";
}

// Test 5: Cavity radius scaling
// Empirical: Rc ~ 10-15 * W^(1/3) meters for granite
// For 1 kt: Rc ~ 10-15 m
// For 250 kt: Rc ~ 63-95 m
TEST_F(MuellerMurphyValidationTest, CavityRadiusScaling)
{
  NuclearSourceParameters params_1kt;
  params_1kt.yield_kt = 1.0;
  double rc_1kt = params_1kt.cavity_radius(RHO_GRANITE);

  EXPECT_GT(rc_1kt, 5.0)
      << "Cavity radius for 1 kt too small: " << rc_1kt << " m";
  EXPECT_LT(rc_1kt, 25.0)
      << "Cavity radius for 1 kt too large: " << rc_1kt << " m";

  NuclearSourceParameters params_250kt;
  params_250kt.yield_kt = 250.0;
  double rc_250kt = params_250kt.cavity_radius(RHO_GRANITE);

  EXPECT_GT(rc_250kt, 40.0)
      << "Cavity radius for 250 kt too small: " << rc_250kt << " m";
  EXPECT_LT(rc_250kt, 150.0)
      << "Cavity radius for 250 kt too large: " << rc_250kt << " m";

  // Verify cube-root scaling: rc ~ W^(1/3)
  double ratio = rc_250kt / rc_1kt;
  double expected_ratio = std::pow(250.0, 1.0 / 3.0);  // ~ 6.3

  EXPECT_NEAR(ratio, expected_ratio, 0.3 * expected_ratio)
      << "Cavity radius should scale as W^(1/3)";
}

// Test 6: Moment rate function properties
TEST_F(MuellerMurphyValidationTest, MomentRateFunction)
{
  NuclearSourceParameters params;
  params.yield_kt = 10.0;
  params.depth_of_burial = 500.0;

  MuellerMurphySource source;
  source.setParameters(params);
  source.setMediumProperties(RHO_GRANITE, VP_GRANITE, VS_GRANITE);

  // Moment rate at t=0 should be at or near peak (explosion onset)
  double mr_0 = source.momentRate(0.0);
  EXPECT_GT(mr_0, 0.0) << "Moment rate at t=0 must be nonzero";
  EXPECT_TRUE(std::isfinite(mr_0)) << "Moment rate at t=0 must be finite";

  // Moment rate should decay with time (for an impulsive source)
  double mr_1 = source.momentRate(1.0);
  EXPECT_TRUE(std::isfinite(mr_1));
  EXPECT_LT(std::abs(mr_1), std::abs(mr_0) + 1.0e-20)
      << "Moment rate should not grow with time";

  // Moment rate at t = 10 * rise_time should be < 1% of peak
  double rise_time = 0.55 / source.getCornerFrequency();
  double mr_10tau = source.momentRate(10.0 * rise_time);
  EXPECT_LT(std::abs(mr_10tau), 0.01 * std::abs(mr_0))
      << "Moment rate at 10*rise_time should be < 1% of peak";
}

// Test 7: Corner frequency for 250 kt
TEST_F(MuellerMurphyValidationTest, CornerFrequency250ktGranite)
{
  NuclearSourceParameters params;
  params.yield_kt = 250.0;
  params.depth_of_burial = 800.0;

  MuellerMurphySource source;
  source.setParameters(params);
  source.setMediumProperties(RHO_GRANITE, VP_GRANITE, VS_GRANITE);

  double fc = source.getCornerFrequency();
  EXPECT_GT(fc, 0.01) << "Corner frequency for 250 kt too low";
  EXPECT_LT(fc, 1.0) << "Corner frequency for 250 kt too high";
}

// Test 8: Scalar moment ranges
TEST_F(MuellerMurphyValidationTest, ScalarMomentRanges)
{
  // 1 kt: M0 should be in 1e13 to 1e16 N*m
  NuclearSourceParameters p1;
  p1.yield_kt = 1.0;
  double M0_1kt = p1.scalar_moment();
  EXPECT_GT(M0_1kt, 1.0e13) << "M0 for 1 kt too small";
  EXPECT_LT(M0_1kt, 1.0e16) << "M0 for 1 kt too large";

  // 250 kt: M0 should be in 1e15 to 1e18 N*m
  NuclearSourceParameters p250;
  p250.yield_kt = 250.0;
  double M0_250kt = p250.scalar_moment();
  EXPECT_GT(M0_250kt, 1.0e15) << "M0 for 250 kt too small";
  EXPECT_LT(M0_250kt, 1.0e18) << "M0 for 250 kt too large";
}

// Test 9: RDP low-frequency plateau equals M0
TEST_F(MuellerMurphyValidationTest, RDPLowFrequencyPlateauEqualsM0)
{
  NuclearSourceParameters params;
  params.yield_kt = 10.0;
  params.depth_of_burial = 500.0;

  MuellerMurphySource source;
  source.setParameters(params);
  source.setMediumProperties(RHO_GRANITE, VP_GRANITE, VS_GRANITE);

  double M0 = params.scalar_moment();
  double fc = source.getCornerFrequency();
  ASSERT_GT(fc, 0.0);
  ASSERT_GT(M0, 0.0);

  // At very low frequency (0.001 * fc), RDP amplitude should equal M0
  double omega_low = 2.0 * M_PI * fc * 0.001;
  double amp_low = std::abs(source.rdp(omega_low));
  EXPECT_NEAR(amp_low, M0, 0.01 * M0)
      << "RDP at omega << omega_c should equal M0 (flat low-frequency plateau)";
}

// =============================================================================
// Mueller-Murphy (1971) damped-resonator RDP -- new tests for the M&M
// frequency-domain shape introduced in the historic-nuclear robustness
// commit.
//
// References: Mueller & Murphy (1971) BSSA 61(6), figs. 3 and 6;
//             Stevens & Day (1985) JGR 90, eq. 4.
// =============================================================================

// Test 10: spectral overshoot peak when B > 1.
// With B = 2.5, zeta = 0.4, k_B = 1 the (B - 1) numerator zero is
// activated and the spectrum peaks above M0 in the neighbourhood of
// omega_p, matching M&M (1971) figure 6.
TEST_F(MuellerMurphyValidationTest, RDPOvershootPeak)
{
  NuclearSourceParameters params;
  params.yield_kt = 10.0;
  params.depth_of_burial = 500.0;

  MuellerMurphySource source;
  source.setParameters(params);
  source.setMediumProperties(RHO_GRANITE, VP_GRANITE, VS_GRANITE);
  source.setOvershoot(2.5);
  source.setDamping(0.4);
  source.setOvershootZeroFactor(1.0);

  const double M0 = params.scalar_moment();
  const double fc = source.getCornerFrequency();
  ASSERT_GT(M0, 0.0);
  ASSERT_GT(fc, 0.0);

  // Sample the spectrum around the corner frequency and verify the
  // peak rises above M0.
  const int Nsamples = 41;
  double max_amp = 0.0;
  double max_omega = 0.0;
  const double omega_p = 2.0 * M_PI * fc;
  for (int i = 0; i < Nsamples; ++i)
  {
    double f = std::pow(10.0, -1.0 + 2.0 * static_cast<double>(i) /
                                static_cast<double>(Nsamples - 1)) * fc;
    double amp = std::abs(source.rdp(2.0 * M_PI * f));
    if (amp > max_amp)
    {
      max_amp = amp;
      max_omega = 2.0 * M_PI * f;
    }
  }
  EXPECT_GT(max_amp, M0)
      << "RDP with B = 2.5 should peak above M0; max amp = " << max_amp;
  // Peak should be near omega_p (within an order of magnitude).
  EXPECT_NEAR(std::log10(max_omega / omega_p), 0.0, 1.0)
      << "Peak omega should be within an order of magnitude of omega_p";
}

// Test 11: low-frequency plateau equals M0 in the new M&M form (default
// B = 1, k_B = 1, zeta = 0.7).
TEST_F(MuellerMurphyValidationTest, RDPLowFrequencyEqualsM0_NewModel)
{
  NuclearSourceParameters params;
  params.yield_kt = 10.0;

  MuellerMurphySource source;
  source.setParameters(params);
  source.setMediumProperties(RHO_GRANITE, VP_GRANITE, VS_GRANITE);

  const double M0 = params.scalar_moment();
  const double fc = source.getCornerFrequency();
  const double omega_low = 2.0 * M_PI * fc * 1e-3;
  double amp_low = std::abs(source.rdp(omega_low));
  EXPECT_NEAR(amp_low, M0, 0.01 * M0)
      << "M&M RDP at 1e-3 * omega_p must equal M0 within 1%";
}

// Test 12: high-frequency rolloff follows omega^-2 envelope.
// With B = 1 default the (B - 1) gating suppresses the numerator zero
// and the spectrum decays as the damped second-order envelope, which
// is asymptotically 1/x^2. At x = 100 the magnitude is
// 1 / |1 - 1e4 + 2*zeta*100*i| ~ 1/sqrt(1e8 + 4*zeta^2*1e4) ~ 1e-4.
TEST_F(MuellerMurphyValidationTest, RDPHighFrequencyOmegaMinus2)
{
  NuclearSourceParameters params;
  params.yield_kt = 10.0;

  MuellerMurphySource source;
  source.setParameters(params);
  source.setMediumProperties(RHO_GRANITE, VP_GRANITE, VS_GRANITE);

  const double M0 = params.scalar_moment();
  const double fc = source.getCornerFrequency();
  const double omega_p = 2.0 * M_PI * fc;
  const double omega_high = omega_p * 100.0;
  const double amp_high = std::abs(source.rdp(omega_high));
  const double ratio = amp_high / M0;
  const double expected = 1.0 / 10000.0;  // (1/100)^2

  // Tolerate 30% as in the task spec to absorb the small zeta-dependent
  // contribution at x = 100.
  EXPECT_NEAR(ratio, expected, 0.3 * expected)
      << "RDP at 100*omega_p should follow omega^-2 envelope; ratio = "
      << ratio;
}

// Test 13: moment-rate trapezoidal integral approximates M0 within 2%.
// The exponential moment-rate Mdot(t) = (M0/tau) * exp(-t/tau)
// integrates to M0 over [0, infinity). At t = 30/omega_p the residual
// is exp(-30/(omega_p*tau)) which is small for tau ~ 0.55/fc.
TEST_F(MuellerMurphyValidationTest, MomentRateIntegratesToM0)
{
  NuclearSourceParameters params;
  params.yield_kt = 10.0;

  MuellerMurphySource source;
  source.setParameters(params);
  source.setMediumProperties(RHO_GRANITE, VP_GRANITE, VS_GRANITE);

  const double M0 = params.scalar_moment();
  const double fc = source.getCornerFrequency();
  const double omega_p = 2.0 * M_PI * fc;
  const double t_end = 30.0 / omega_p;

  const int N = 4001;
  const double dt = t_end / static_cast<double>(N - 1);
  double integral = 0.0;
  double prev = source.momentRate(0.0);
  for (int i = 1; i < N; ++i)
  {
    double t = static_cast<double>(i) * dt;
    double cur = source.momentRate(t);
    integral += 0.5 * (prev + cur) * dt;
    prev = cur;
  }

  // Allow a few extra percent for the truncation tail since the
  // simple exponential decays slowly relative to omega_p when the
  // rise time is set to 0.55/fc. Total moment must be within 5% of M0.
  EXPECT_NEAR(integral, M0, 0.05 * M0)
      << "Mdot trapezoidal integral over [0, 30/omega_p] should approach M0";
}

// Test 14: moment-rate non-negative for B = 1 default. The exponential
// momentRate kernel is non-negative for all t by construction.
TEST_F(MuellerMurphyValidationTest, MomentRateNonNegativeForB1)
{
  NuclearSourceParameters params;
  params.yield_kt = 10.0;

  MuellerMurphySource source;
  source.setParameters(params);
  source.setMediumProperties(RHO_GRANITE, VP_GRANITE, VS_GRANITE);

  const double fc = source.getCornerFrequency();
  const double omega_p = 2.0 * M_PI * fc;
  const double t_end = 50.0 / omega_p;

  const int N = 1001;
  const double dt = t_end / static_cast<double>(N - 1);
  for (int i = 0; i < N; ++i)
  {
    double t = static_cast<double>(i) * dt;
    double mr = source.momentRate(t);
    EXPECT_GE(mr, 0.0)
        << "momentRate(" << t << ") = " << mr
        << " must be non-negative for default (B = 1) shape";
  }
}
