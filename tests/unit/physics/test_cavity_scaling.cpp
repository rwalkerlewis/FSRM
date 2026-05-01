/**
 * @file test_cavity_scaling.cpp
 * @brief Medium-aware cavity radius scaling for underground nuclear shots
 *
 * Verifies the literature-backed coefficients in
 * NuclearSourceParameters::cavityCoefficient and the medium-aware
 * cavity_radius overload added in Commit 1 of the historic-nuclear
 * robustness work.
 *
 * Coefficient table (m / kt^(1/3)):
 *   GRANITE   11    Boardman et al. (1964) JGR 69(16); Closmann (1969).
 *   TUFF      18    Closmann (1969); Glenn & Goldstein (1994).
 *   SALT      16    Glenn & Goldstein (1994); Salmon, Sterling tests.
 *   ALLUVIUM  22    Closmann (1969); weakly cemented sediments.
 *   SHALE     14    Carter et al. (1968) Gasbuggy; Lewis Shale estimate.
 *   GENERIC   12    Legacy mid-range value used by single-arg overload.
 *
 * The single-argument legacy overload cavity_radius(rho) must continue
 * to return the GENERIC coefficient result, and the cube-root scaling
 * (Rc ~ W^(1/3)) must hold for every medium.
 */

#include <gtest/gtest.h>
#include <cmath>

#include "domain/explosion/ExplosionImpactPhysics.hpp"

using namespace FSRM;
using MediumType = NuclearSourceParameters::MediumType;

class CavityScalingTest : public ::testing::Test
{
protected:
  static constexpr double GRANITE_RHO = 2650.0;  // kg/m^3
  static constexpr double TUFF_RHO    = 2300.0;
  static constexpr double SALT_RHO    = 2200.0;
  static constexpr double ALLUVIUM_RHO = 1900.0;
  static constexpr double SHALE_RHO    = 2500.0;
};

// 1 kt at granite reference density: cavity radius in [9, 13] m.
TEST_F(CavityScalingTest, GraniteOneKt)
{
  NuclearSourceParameters params;
  params.yield_kt = 1.0;
  double rc = params.cavity_radius(GRANITE_RHO, MediumType::GRANITE);
  EXPECT_GE(rc, 9.0)  << "Granite 1 kt cavity radius too small: " << rc;
  EXPECT_LE(rc, 13.0) << "Granite 1 kt cavity radius too large: " << rc;
}

// 1 kt in tuff: cavity radius in [16, 20] m.
TEST_F(CavityScalingTest, TuffOneKt)
{
  NuclearSourceParameters params;
  params.yield_kt = 1.0;
  // Use reference density so the cube-root density correction is unity
  // and we test only the medium coefficient.
  double rc = params.cavity_radius(GRANITE_RHO, MediumType::TUFF);
  EXPECT_GE(rc, 16.0) << "Tuff 1 kt cavity radius too small: " << rc;
  EXPECT_LE(rc, 20.0) << "Tuff 1 kt cavity radius too large: " << rc;
}

// 1 kt in salt: cavity radius in [14, 18] m.
TEST_F(CavityScalingTest, SaltOneKt)
{
  NuclearSourceParameters params;
  params.yield_kt = 1.0;
  double rc = params.cavity_radius(GRANITE_RHO, MediumType::SALT);
  EXPECT_GE(rc, 14.0) << "Salt 1 kt cavity radius too small: " << rc;
  EXPECT_LE(rc, 18.0) << "Salt 1 kt cavity radius too large: " << rc;
}

// Cube-root scaling holds for every medium within 5%.
TEST_F(CavityScalingTest, CubeRootScalingAllMedia)
{
  const std::vector<MediumType> media = {
    MediumType::GRANITE, MediumType::TUFF, MediumType::SALT,
    MediumType::ALLUVIUM, MediumType::SHALE, MediumType::GENERIC
  };

  const double W_small = 1.0;
  const double W_large = 250.0;
  const double expected_ratio = std::cbrt(W_large / W_small);

  for (auto medium : media)
  {
    NuclearSourceParameters p_small, p_large;
    p_small.yield_kt = W_small;
    p_large.yield_kt = W_large;

    double rc_small = p_small.cavity_radius(GRANITE_RHO, medium);
    double rc_large = p_large.cavity_radius(GRANITE_RHO, medium);

    ASSERT_GT(rc_small, 0.0);
    double ratio = rc_large / rc_small;
    EXPECT_NEAR(ratio, expected_ratio, 0.05 * expected_ratio)
        << "Cube-root scaling violated for medium index "
        << static_cast<int>(medium)
        << ": ratio=" << ratio << " expected=" << expected_ratio;
  }
}

// Backward-compatibility regression: legacy single-argument overload
// returns the GENERIC (C = 12) result and continues to satisfy the
// historical [5, 25] / [40, 150] ranges from
// MuellerMurphyValidationTest.CavityRadiusScaling.
TEST_F(CavityScalingTest, LegacyOverloadMatchesGeneric)
{
  NuclearSourceParameters p_1kt;
  p_1kt.yield_kt = 1.0;

  double rc_legacy = p_1kt.cavity_radius(GRANITE_RHO);
  double rc_generic =
      p_1kt.cavity_radius(GRANITE_RHO, MediumType::GENERIC);
  EXPECT_DOUBLE_EQ(rc_legacy, rc_generic)
      << "Legacy single-arg overload must match GENERIC medium overload";

  // Also verify the historic loose bounds still hold.
  EXPECT_GE(rc_legacy, 5.0);
  EXPECT_LE(rc_legacy, 25.0);

  NuclearSourceParameters p_250kt;
  p_250kt.yield_kt = 250.0;
  double rc_legacy_250 = p_250kt.cavity_radius(GRANITE_RHO);
  EXPECT_GE(rc_legacy_250, 40.0);
  EXPECT_LE(rc_legacy_250, 150.0);
}

// Verify the static cavityCoefficient table directly.
TEST_F(CavityScalingTest, CoefficientTableValues)
{
  EXPECT_DOUBLE_EQ(NuclearSourceParameters::cavityCoefficient(MediumType::GRANITE),  11.0);
  EXPECT_DOUBLE_EQ(NuclearSourceParameters::cavityCoefficient(MediumType::TUFF),     18.0);
  EXPECT_DOUBLE_EQ(NuclearSourceParameters::cavityCoefficient(MediumType::SALT),     16.0);
  EXPECT_DOUBLE_EQ(NuclearSourceParameters::cavityCoefficient(MediumType::ALLUVIUM), 22.0);
  EXPECT_DOUBLE_EQ(NuclearSourceParameters::cavityCoefficient(MediumType::SHALE),    14.0);
  EXPECT_DOUBLE_EQ(NuclearSourceParameters::cavityCoefficient(MediumType::GENERIC),  12.0);
}

// Pass-2: crushed_zone_radius(medium) routes the medium argument through
// cavity_radius. Salt (C = 16) and granite (C = 11) at the same yield
// must produce a crushed-zone ratio matching the coefficient ratio
// (16/11) within 1% (the legacy single-arg crushed_zone_radius hard-
// coded the GENERIC coefficient and ignored medium entirely).
TEST_F(CavityScalingTest, CrushedZoneScalesWithMedium)
{
  NuclearSourceParameters params;
  params.yield_kt = 1.0;

  const double rc_salt = params.crushed_zone_radius(MediumType::SALT);
  const double rc_granite = params.crushed_zone_radius(MediumType::GRANITE);
  ASSERT_GT(rc_salt, 0.0);
  ASSERT_GT(rc_granite, 0.0);

  const double ratio = rc_salt / rc_granite;
  const double expected = 16.0 / 11.0;
  EXPECT_NEAR(ratio, expected, 0.01 * expected)
      << "crushed_zone_radius(SALT)/crushed_zone_radius(GRANITE) "
      << "must equal C_SALT/C_GRANITE = 16/11 within 1%";
}

// Pass-2: same scaling guarantee for the fractured-zone radius.
TEST_F(CavityScalingTest, FracturedZoneScalesWithMedium)
{
  NuclearSourceParameters params;
  params.yield_kt = 1.0;

  const double rf_salt = params.fractured_zone_radius(MediumType::SALT);
  const double rf_granite = params.fractured_zone_radius(MediumType::GRANITE);
  ASSERT_GT(rf_salt, 0.0);
  ASSERT_GT(rf_granite, 0.0);

  const double ratio = rf_salt / rf_granite;
  const double expected = 16.0 / 11.0;
  EXPECT_NEAR(ratio, expected, 0.01 * expected)
      << "fractured_zone_radius(SALT)/fractured_zone_radius(GRANITE) "
      << "must equal C_SALT/C_GRANITE = 16/11 within 1%";
}

// Pass-2 backward-compat: legacy zero-argument crushed_zone_radius and
// fractured_zone_radius continue to use the GENERIC coefficient
// (C = 12) so existing callers that did not pass a medium see
// unchanged behaviour. Regression check guards against accidental
// behavior shift in the legacy overloads.
TEST_F(CavityScalingTest, LegacyCrushedZoneStable)
{
  NuclearSourceParameters params;
  params.yield_kt = 1.0;

  const double rc_legacy = params.crushed_zone_radius();
  const double rc_generic = params.crushed_zone_radius(MediumType::GENERIC);
  EXPECT_DOUBLE_EQ(rc_legacy, rc_generic)
      << "Legacy zero-arg crushed_zone_radius must match GENERIC overload";

  const double rf_legacy = params.fractured_zone_radius();
  const double rf_generic = params.fractured_zone_radius(MediumType::GENERIC);
  EXPECT_DOUBLE_EQ(rf_legacy, rf_generic)
      << "Legacy zero-arg fractured_zone_radius must match GENERIC overload";

  // Numerical sanity: 3 * cavity_radius and 10 * cavity_radius for 1 kt
  // at GENERIC (C = 12, rho_ref = 2650): cavity_radius is 12 m so
  // crushed = 36 m and fractured = 120 m exactly.
  const double rc_cavity = params.cavity_radius(2650.0, MediumType::GENERIC);
  EXPECT_DOUBLE_EQ(rc_legacy, 3.0 * rc_cavity);
  EXPECT_DOUBLE_EQ(rf_legacy, 10.0 * rc_cavity);
}

// MuellerMurphySource::setMedium plumbs the medium into
// computeDerivedQuantities so the scalar moment reflects the medium-aware
// cavity radius.
TEST_F(CavityScalingTest, MuellerMurphyHonorsSetMedium)
{
  NuclearSourceParameters params;
  params.yield_kt = 10.0;

  MuellerMurphySource generic_src;
  generic_src.setParameters(params);
  generic_src.setMediumProperties(GRANITE_RHO, 5500.0, 3100.0);
  // No setMedium call -> defaults to GENERIC

  MuellerMurphySource tuff_src;
  tuff_src.setParameters(params);
  tuff_src.setMediumProperties(GRANITE_RHO, 5500.0, 3100.0);
  tuff_src.setMedium(MediumType::TUFF);

  // Larger cavity coefficient for TUFF (18 vs 12) gives a larger M0 by
  // factor (18/12)^3 = 3.375. Tolerate 1% slack for any internal
  // arithmetic differences.
  double M0_generic = std::abs(generic_src.rdp(0.0));
  double M0_tuff    = std::abs(tuff_src.rdp(0.0));
  double ratio = M0_tuff / M0_generic;
  double expected_ratio = std::pow(18.0 / 12.0, 3.0);
  EXPECT_NEAR(ratio, expected_ratio, 0.01 * expected_ratio)
      << "setMedium must rescale scalar moment via medium-aware cavity radius";
}
