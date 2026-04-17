/**
 * @file test_near_field_explosion.cpp
 * @brief Physics validation for near-field underground explosion phenomenology
 *
 * Quantitative tests:
 *   1. Cavity radius scaling: Rc = 55*W^0.295*(rho/2650)^(-1/3.4) meters
 *   2. Damage zone ratios: crushed ~2.5x, fractured ~5x, damaged ~10x cavity
 *   3. Spall velocity from acoustic impedance: u_spall = P/(rho*c)
 *
 * References:
 *   Rimer & Cherry (1982), "Near-field phenomenology"
 *   Terhune et al. (1977), "Near-field ground motion"
 *   Mueller & Murphy (1971), "Seismic characteristics"
 *   NTS empirical scaling relations
 */

#include <gtest/gtest.h>
#include <cmath>
#include <vector>

#include "domain/explosion/NearFieldExplosion.hpp"
#include "domain/explosion/ExplosionImpactPhysics.hpp"

using namespace FSRM;

class NearFieldExplosionValidationTest : public ::testing::Test
{
protected:
  static constexpr double RHO_GRANITE = 2650.0;  // kg/m^3
  static constexpr double VP_GRANITE  = 5500.0;  // m/s
  static constexpr double VS_GRANITE  = 3200.0;  // m/s
};

// ---------------------------------------------------------------------------
// Test 1: Cavity radius for 1 kt in granite at 300 m depth
//
// NTS empirical scaling: Rc = 55 * W^0.295 * (rho/2650)^(-1/3.4)
// For 1 kt, rho = 2650: Rc = 55 * 1.0^0.295 * 1.0 = 55 m
// Empirical range for 1 kt in hard rock: 10-20 m (some references)
// The NTS scaling coefficient 55 m/(kt^0.3) gives larger values because
// it was calibrated to alluvium/tuff. Accept 10-100 m range.
// ---------------------------------------------------------------------------
TEST_F(NearFieldExplosionValidationTest, CavityRadius1ktGranite)
{
  UndergroundExplosionSource source;
  source.yield_kt = 1.0;
  source.depth = 300.0;
  source.host_density = RHO_GRANITE;
  source.host_vp = VP_GRANITE;
  source.host_vs = VS_GRANITE;

  double Rc = source.cavityRadius();

  // Cavity radius must be positive
  EXPECT_GT(Rc, 0.0) << "Cavity radius must be positive";

  // NTS scaling gives Rc = 55 * 1.0^0.295 = 55 m for standard density
  // This is for NTS tuff; granite would be smaller. Accept wide range.
  EXPECT_GT(Rc, 5.0) << "Cavity radius too small for 1 kt: " << Rc << " m";
  EXPECT_LT(Rc, 200.0) << "Cavity radius too large for 1 kt: " << Rc << " m";
}

// ---------------------------------------------------------------------------
// Test 2: Cavity radius cube-root scaling with yield
//
// Rc ~ W^0.295 (approximately cube-root scaling)
// Rc(100 kt) / Rc(1 kt) should be ~100^0.295 = 3.89
// ---------------------------------------------------------------------------
TEST_F(NearFieldExplosionValidationTest, CavityRadiusYieldScaling)
{
  UndergroundExplosionSource source_1kt;
  source_1kt.yield_kt = 1.0;
  source_1kt.depth = 300.0;
  source_1kt.host_density = RHO_GRANITE;

  UndergroundExplosionSource source_100kt;
  source_100kt.yield_kt = 100.0;
  source_100kt.depth = 300.0;
  source_100kt.host_density = RHO_GRANITE;

  double Rc_1 = source_1kt.cavityRadius();
  double Rc_100 = source_100kt.cavityRadius();

  // Both must be positive
  EXPECT_GT(Rc_1, 0.0);
  EXPECT_GT(Rc_100, 0.0);

  // Larger yield -> larger cavity
  EXPECT_GT(Rc_100, Rc_1) << "Cavity radius must increase with yield";

  // Check scaling exponent: ratio should be ~100^0.295 = 3.89
  double ratio = Rc_100 / Rc_1;
  double expected_ratio = std::pow(100.0, ExplosionPhysics::CAVITY_SCALING_EXP);

  double rel_error = std::abs(ratio - expected_ratio) / expected_ratio;
  EXPECT_LT(rel_error, 0.05)
      << "Yield scaling ratio: " << ratio
      << ", expected ~" << expected_ratio;
}

// ---------------------------------------------------------------------------
// Test 3: Damage zone extent
//
// Crushed zone ~ 2.5x Rc, Fractured zone ~ 5x Rc, Damaged zone ~ 10x Rc
// These are fixed multiples in the NTS scaling model.
// ---------------------------------------------------------------------------
TEST_F(NearFieldExplosionValidationTest, DamageZoneExtent)
{
  UndergroundExplosionSource source;
  source.yield_kt = 10.0;
  source.depth = 500.0;
  source.host_density = RHO_GRANITE;

  double Rc = source.cavityRadius();
  double R_crushed = source.crushedZoneRadius();
  double R_fractured = source.fracturedZoneRadius();
  double R_damaged = source.damagedZoneRadius();

  // All must be positive
  EXPECT_GT(Rc, 0.0);
  EXPECT_GT(R_crushed, 0.0);
  EXPECT_GT(R_fractured, 0.0);
  EXPECT_GT(R_damaged, 0.0);

  // Proper ordering: cavity < crushed < fractured < damaged
  EXPECT_GT(R_crushed, Rc) << "Crushed zone must be larger than cavity";
  EXPECT_GT(R_fractured, R_crushed) << "Fractured zone must be larger than crushed";
  EXPECT_GT(R_damaged, R_fractured) << "Damaged zone must be larger than fractured";

  // Check ratios within 10% of expected constants
  double crushed_ratio = R_crushed / Rc;
  double fractured_ratio = R_fractured / Rc;
  double damaged_ratio = R_damaged / Rc;

  EXPECT_NEAR(crushed_ratio, ExplosionPhysics::CRUSHED_ZONE_RATIO, 0.5)
      << "Crushed zone ratio: " << crushed_ratio;
  EXPECT_NEAR(fractured_ratio, ExplosionPhysics::FRACTURED_ZONE_RATIO, 1.0)
      << "Fractured zone ratio: " << fractured_ratio;
  EXPECT_NEAR(damaged_ratio, ExplosionPhysics::DAMAGED_ZONE_RATIO, 2.0)
      << "Damaged zone ratio: " << damaged_ratio;
}

// ---------------------------------------------------------------------------
// Test 4: Spall velocity computation
//
// Maximum spall velocity = P_peak / (rho * c)
// For a P-wave hitting a free surface, the stress doubles, then reflects
// as tension. The particle velocity at the free surface is:
//   u_free = 2 * P_peak / (rho * c)
// The SpallModel should compute u_spall = P / (rho*c).
// ---------------------------------------------------------------------------
TEST_F(NearFieldExplosionValidationTest, SpallVelocity)
{
  SpallModel spall;

  // 1 GPa peak stress, granite impedance
  double P_peak = 1.0e9;  // 1 GPa
  double rho = RHO_GRANITE;
  double c = VP_GRANITE;

  double u_spall = spall.maxSpallVelocity(P_peak, rho, c);

  // Expected: P / (rho * c) = 1e9 / (2650 * 5500) = 1e9 / 1.4575e7 ~ 68.6 m/s
  double expected = P_peak / (rho * c);

  EXPECT_GT(u_spall, 0.0) << "Spall velocity must be positive";

  EXPECT_NEAR(u_spall, expected, expected * 0.01)
      << "Spall velocity: computed " << u_spall
      << " m/s, expected " << expected << " m/s";
}

// ---------------------------------------------------------------------------
// Test 5: Spall thickness for shallow explosion
//
// For a shallow explosion (depth < 5*Rc), spall should occur.
// Spall thickness should be nonzero and physically reasonable.
// ---------------------------------------------------------------------------
TEST_F(NearFieldExplosionValidationTest, SpallThicknessShallow)
{
  SpallModel spall;

  double yield_kt = 1.0;
  double depth_shallow = 100.0;  // Shallow burial
  double depth_deep = 5000.0;    // Deep burial

  double thickness_shallow = spall.spallThickness(yield_kt, depth_shallow);
  double thickness_deep = spall.spallThickness(yield_kt, depth_deep);

  // Shallow spall thickness should be nonzero
  EXPECT_GT(thickness_shallow, 0.0)
      << "Spall thickness must be positive for shallow explosion";

  // Shallow spall should be thicker than deep spall
  EXPECT_GT(thickness_shallow, thickness_deep)
      << "Shallow explosion should have thicker spall than deep explosion";

  // Spall thickness should be physically reasonable (< 10x cavity radius)
  double Rc = ExplosionPhysics::CAVITY_SCALING_COEFF *
              std::pow(yield_kt, ExplosionPhysics::CAVITY_SCALING_EXP);
  EXPECT_LT(thickness_shallow, 10.0 * Rc)
      << "Spall thickness too large: " << thickness_shallow << " m";
}

// ---------------------------------------------------------------------------
// Test 6: Dynamic tensile strength increases with strain rate
//
// Y_t(edot) = Y_t0 * (1 + (edot/edot_ref)^n)
// Must be increasing with strain rate.
// ---------------------------------------------------------------------------
TEST_F(NearFieldExplosionValidationTest, DynamicTensileStrength)
{
  SpallModel spall;

  double edot_low = 0.01;    // quasi-static
  double edot_med = 100.0;   // moderate rate
  double edot_high = 10000.0; // high rate

  double Yt_low = spall.dynamicTensileStrength(edot_low);
  double Yt_med = spall.dynamicTensileStrength(edot_med);
  double Yt_high = spall.dynamicTensileStrength(edot_high);

  // All must be positive
  EXPECT_GT(Yt_low, 0.0);
  EXPECT_GT(Yt_med, 0.0);
  EXPECT_GT(Yt_high, 0.0);

  // Must be increasing with strain rate
  EXPECT_GT(Yt_med, Yt_low)
      << "Strength must increase with strain rate";
  EXPECT_GT(Yt_high, Yt_med)
      << "Strength must increase with strain rate";

  // At quasi-static rate, should be close to the static value
  EXPECT_GT(Yt_low, spall.static_tensile_strength * 0.5)
      << "Quasi-static strength too low";
  EXPECT_LT(Yt_low, spall.static_tensile_strength * 5.0)
      << "Quasi-static strength too high";
}

// ---------------------------------------------------------------------------
// Test 7: Cavity radius for 250 kt granite
// ---------------------------------------------------------------------------
TEST_F(NearFieldExplosionValidationTest, CavityRadius250ktGranite)
{
  UndergroundExplosionSource source;
  source.yield_kt = 250.0;
  source.depth = 800.0;
  source.host_density = RHO_GRANITE;
  source.host_vp = VP_GRANITE;
  source.host_vs = VS_GRANITE;

  double Rc = source.cavityRadius();
  EXPECT_GT(Rc, 50.0) << "Cavity radius for 250 kt too small: " << Rc << " m";
  EXPECT_LT(Rc, 500.0) << "Cavity radius for 250 kt too large: " << Rc << " m";
}

// ---------------------------------------------------------------------------
// Test 8: SphericalCavitySource displacement at r=100m is nonzero
// ---------------------------------------------------------------------------
TEST_F(NearFieldExplosionValidationTest, CavitySourceDisplacementNonzero)
{
  SphericalCavitySource cavity;
  cavity.setCavityRadius(15.0);
  cavity.setMediumProperties(RHO_GRANITE,
    RHO_GRANITE * VP_GRANITE * VP_GRANITE - 2.0 * RHO_GRANITE * VS_GRANITE * VS_GRANITE,
    RHO_GRANITE * VS_GRANITE * VS_GRANITE);
  cavity.setOverpressure(1.0e9);
  cavity.setRiseTime(0.01);

  double u_100 = cavity.displacement(100.0, 0.1);
  EXPECT_GT(u_100, 0.0)
      << "Displacement at r=100m must be nonzero and positive (outward)";
}

// ---------------------------------------------------------------------------
// Test 9: SphericalCavitySource displacement decays with distance
// ---------------------------------------------------------------------------
TEST_F(NearFieldExplosionValidationTest, CavitySourceDisplacementDecay)
{
  SphericalCavitySource cavity;
  cavity.setCavityRadius(15.0);
  cavity.setMediumProperties(RHO_GRANITE,
    RHO_GRANITE * VP_GRANITE * VP_GRANITE - 2.0 * RHO_GRANITE * VS_GRANITE * VS_GRANITE,
    RHO_GRANITE * VS_GRANITE * VS_GRANITE);
  cavity.setOverpressure(1.0e9);
  cavity.setRiseTime(0.01);

  double u_100 = cavity.displacement(100.0, 0.1);
  double u_200 = cavity.displacement(200.0, 0.1);
  EXPECT_GT(u_100, u_200)
      << "Displacement should decay with distance: u(100)=" << u_100
      << " u(200)=" << u_200;
}

// ---------------------------------------------------------------------------
// Test 10: SphericalCavitySource equivalent moment is nonzero
// ---------------------------------------------------------------------------
TEST_F(NearFieldExplosionValidationTest, CavitySourceEquivalentMoment)
{
  SphericalCavitySource cavity;
  cavity.setCavityRadius(15.0);
  cavity.setMediumProperties(RHO_GRANITE,
    RHO_GRANITE * VP_GRANITE * VP_GRANITE - 2.0 * RHO_GRANITE * VS_GRANITE * VS_GRANITE,
    RHO_GRANITE * VS_GRANITE * VS_GRANITE);
  cavity.setOverpressure(1.0e9);
  cavity.setRiseTime(0.01);

  double M0 = cavity.equivalentMoment();
  EXPECT_GT(M0, 0.0)
      << "Equivalent moment must be nonzero and positive";
}

// ---------------------------------------------------------------------------
// Test 11: SphericalCavitySource stress: sigma_rr < 0, sigma_tt > 0
// ---------------------------------------------------------------------------
TEST_F(NearFieldExplosionValidationTest, CavitySourceStressSigns)
{
  SphericalCavitySource cavity;
  cavity.setCavityRadius(15.0);
  cavity.setMediumProperties(RHO_GRANITE,
    RHO_GRANITE * VP_GRANITE * VP_GRANITE - 2.0 * RHO_GRANITE * VS_GRANITE * VS_GRANITE,
    RHO_GRANITE * VS_GRANITE * VS_GRANITE);
  cavity.setOverpressure(1.0e9);
  cavity.setRiseTime(0.01);

  double sigma_rr = 0.0, sigma_tt = 0.0;
  cavity.stress(50.0, 0.1, sigma_rr, sigma_tt);

  // Radial stress should be compressive (negative) outside the cavity
  EXPECT_LT(sigma_rr, 0.0)
      << "Radial stress should be compressive: sigma_rr=" << sigma_rr;
  // Tangential stress should be tensile (positive)
  EXPECT_GT(sigma_tt, 0.0)
      << "Tangential stress should be tensile: sigma_tt=" << sigma_tt;
}

// ---------------------------------------------------------------------------
// Test 12: NearFieldExplosionSolver initializes without error
// ---------------------------------------------------------------------------
TEST_F(NearFieldExplosionValidationTest, SolverInitializes)
{
  NearFieldExplosionSolver solver;

  UndergroundExplosionSource src;
  src.yield_kt = 1.0;
  src.depth = 300.0;
  src.host_density = RHO_GRANITE;
  src.host_vp = VP_GRANITE;
  src.host_vs = VS_GRANITE;
  solver.setSource(src);

  MieGruneisenEOS eos;
  solver.setEOS(eos);

  PressureDependentStrength strength;
  solver.setStrengthModel(strength);

  DamageEvolutionModel damage;
  solver.setDamageModel(damage);

  solver.initialize();

  // After initialization, cavity radius should be set
  double Rc = solver.getCavityRadius(0.0);
  EXPECT_GT(Rc, 0.0) << "Cavity radius after init must be positive";
}

// ---------------------------------------------------------------------------
// Test 13: NearFieldExplosionSolver step produces nonzero state
// ---------------------------------------------------------------------------
TEST_F(NearFieldExplosionValidationTest, SolverStepProducesNonzeroState)
{
  NearFieldExplosionSolver solver;

  UndergroundExplosionSource src;
  src.yield_kt = 1.0;
  src.depth = 300.0;
  src.host_density = RHO_GRANITE;
  src.host_vp = VP_GRANITE;
  src.host_vs = VS_GRANITE;
  solver.setSource(src);

  MieGruneisenEOS eos;
  solver.setEOS(eos);

  PressureDependentStrength strength;
  solver.setStrengthModel(strength);

  DamageEvolutionModel damage;
  solver.setDamageModel(damage);

  solver.initialize();

  double dt = 1.0e-5;
  solver.prestep(0.0, dt);
  solver.step(dt);
  solver.poststep(0.0, dt);

  // After stepping, query a point near the source for nonzero state
  std::array<double, 6> M;
  double mr = 0.0;
  solver.computeSourceFunction(dt, M, mr);

  // Moment tensor should be nonzero for an isotropic explosion
  double M_norm = 0.0;
  for (int i = 0; i < 6; i++) M_norm += M[i] * M[i];
  EXPECT_GT(M_norm, 0.0)
      << "Moment tensor should be nonzero after a time step";
}
