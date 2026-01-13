/**
 * @file test_near_field_explosion.cpp
 * @brief Unit tests for underground nuclear test near-field modeling
 */
 
#include <gtest/gtest.h>
 
#include "NearFieldExplosion.hpp"
 
#include <cmath>
#include <map>
#include <string>
 
using namespace FSRM;
 
namespace {
double relTol(double x, double rel = 1e-9, double abs_min = 1e-12) {
  return std::max(abs_min, std::abs(x) * rel);
}
}  // namespace
 
TEST(NearFieldExplosionTest, UndergroundSourceScalingAndZones) {
  UndergroundExplosionSource src;
  src.yield_kt = 150.0;
  src.host_density = 2650.0;
 
  const double Rc = src.cavityRadius();
  ASSERT_GT(Rc, 0.0);
 
  const double expected_Rc =
    ExplosionPhysics::CAVITY_SCALING_COEFF *
    std::pow(src.yield_kt, ExplosionPhysics::CAVITY_SCALING_EXP) *
    std::pow(src.host_density / 2650.0, -1.0 / 3.4);
  EXPECT_NEAR(Rc, expected_Rc, relTol(expected_Rc));
 
  EXPECT_GT(src.crushedZoneRadius(), Rc);
  EXPECT_GT(src.fracturedZoneRadius(), src.crushedZoneRadius());
  EXPECT_GT(src.damagedZoneRadius(), src.fracturedZoneRadius());
 
  const double Pc = src.initialCavityPressure();
  EXPECT_GT(Pc, 0.0);
 
  // Seismic scaling invariants
  EXPECT_GT(src.scalarMoment(), 0.0);
  EXPECT_GT(src.bodyWaveMagnitude(), 0.0);
}
 
TEST(ShockAttenuationModelTest, PeakPressureIsMonotoneDecreasing) {
  UndergroundExplosionSource src;
  src.yield_kt = 20.0;
  src.host_density = 2650.0;
  src.host_vp = 5500.0;
 
  ShockAttenuationModel shock;
  shock.setFromSource(src);
 
  ASSERT_GT(shock.r0, 0.0);
  ASSERT_GT(shock.P0, 0.0);
 
  EXPECT_NEAR(shock.peakPressure(shock.r0), shock.P0, relTol(shock.P0));
 
  const double P2 = shock.peakPressure(2.0 * shock.r0);
  const double P4 = shock.peakPressure(4.0 * shock.r0);
  const double P8 = shock.peakPressure(8.0 * shock.r0);
 
  EXPECT_GT(P2, P4);
  EXPECT_GT(P4, P8);
  EXPECT_GT(P2, 0.0);
  EXPECT_GT(shock.arrivalTime(8.0 * shock.r0), 0.0);
}
 
TEST(RDPSeismicSourceTest, CornerFrequencyAndTimeFunctionBehavior) {
  UndergroundExplosionSource src;
  src.yield_kt = 8.0;
  src.host_density = 2650.0;
  src.host_vp = 5500.0;
  src.host_vs = 3200.0;
 
  RDPSeismicSource rdp;
  rdp.initialize(src);
 
  const double expected_fc = 2.5 * std::pow(src.yield_kt, -1.0 / 3.0);
  EXPECT_NEAR(rdp.corner_frequency, expected_fc, relTol(expected_fc, 1e-12));
  EXPECT_GT(rdp.psi_infinity, 0.0);
 
  EXPECT_DOUBLE_EQ(rdp.psi(0.0), 0.0);
  EXPECT_DOUBLE_EQ(rdp.psiDot(0.0), 0.0);
 
  const double tau = 1.0 / (2.0 * M_PI * rdp.corner_frequency);
  const double psi_far = rdp.psi(10.0 * tau);
  EXPECT_NEAR(psi_far, rdp.psi_infinity * rdp.overshoot, relTol(rdp.psi_infinity * rdp.overshoot, 1e-3));
 
  std::array<double, 6> M{};
  rdp.momentTensor(10.0 * tau, M, 0.7, 0.25);
  EXPECT_DOUBLE_EQ(M[3], 0.0);
  EXPECT_DOUBLE_EQ(M[4], 0.0);
  EXPECT_DOUBLE_EQ(M[5], 0.0);
}
 
TEST(NearFieldExplosionConfigTest, ParsesKeyValueConfigAndBuildsSolver) {
  NearFieldExplosionConfig cfg;
  std::map<std::string, std::string> kv{
    {"yield_kt", "12.5"},
    {"depth_of_burial", "800.0"},
    {"source_x", "10.0"},
    {"source_y", "20.0"},
    {"host_density", "2700.0"},
    {"p_wave_velocity", "6000.0"},
    {"s_wave_velocity", "3464.0"},
    {"overshoot_factor", "1.15"},
  };
 
  ASSERT_TRUE(cfg.parse(kv));
  auto solver = cfg.createSolver();
  ASSERT_NE(solver, nullptr);
 
  const auto& src = solver->getSource();
  EXPECT_NEAR(src.yield_kt, 12.5, 1e-12);
  EXPECT_NEAR(src.depth, 800.0, 1e-12);
  EXPECT_NEAR(src.location[0], 10.0, 1e-12);
  EXPECT_NEAR(src.location[1], 20.0, 1e-12);
  EXPECT_NEAR(src.location[2], -800.0, 1e-12);
  EXPECT_NEAR(src.host_density, 2700.0, 1e-12);
  EXPECT_NEAR(src.host_vp, 6000.0, 1e-12);
  EXPECT_NEAR(src.host_vs, 3464.0, 1e-12);
}
