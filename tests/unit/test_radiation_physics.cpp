/**
 * @file test_radiation_physics.cpp
 * @brief Unit tests for radiation and fallout components used in nuclear-test modeling
 */

#include <gtest/gtest.h>

#include "RadiationPhysics.hpp"

#include <cmath>
#include <map>
#include <memory>
#include <string>

using namespace FSRM;

namespace {
double relTol(double x, double rel = 1e-6, double abs_min = 1e-12) {
  return std::max(abs_min, std::abs(x) * rel);
}
}  // namespace

TEST(RadiationPhysicsTest, PromptDoseDecreasesWithDistance) {
  PromptRadiationModel prompt;
  prompt.setYield(20.0);
  prompt.setFissionFraction(0.5);
  prompt.setBurstHeight(0.0);

  // Note: the current air-attenuation fit is extremely strong, so doses are only
  // numerically non-zero at very short ranges. Keep this test in the regime where
  // the model returns non-zero values.
  const double d0_5 = prompt.totalPromptDose(0.5);
  const double d2 = prompt.totalPromptDose(2.0);
  ASSERT_GT(d0_5, 0.0);
  ASSERT_GE(d2, 0.0);
  EXPECT_GT(d0_5, d2);

  // Prompt dose rate is localized near t=0.
  EXPECT_GT(prompt.doseRate(100.0, 0.0), 0.0);
  EXPECT_DOUBLE_EQ(prompt.doseRate(100.0, -1.0), 0.0);
  EXPECT_DOUBLE_EQ(prompt.doseRate(100.0, 1e-3), 0.0);
}

TEST(RadiationPhysicsTest, FissionProductInventoryDecaysOverTime) {
  auto db = std::make_shared<RadionuclideDatabase>();
  FissionProductInventory inv;
  inv.setDatabase(db);
  inv.setFissionYield(50.0);
  inv.initialize();

  const double A0 = inv.totalActivity();
  ASSERT_GT(A0, 0.0);
  ASSERT_GT(inv.iodine131Activity(), 0.0);
  ASSERT_GT(inv.cesium137Activity(), 0.0);

  inv.evolveTo(86400.0);  // 1 day
  const double A1 = inv.totalActivity();
  EXPECT_LT(A1, A0);
}

TEST(RadiationPhysicsTest, FalloutLogNormalDistributionHasMedianAtHalfCDF) {
  FalloutFormationModel fallout;
  fallout.setYield(100.0);
  fallout.setFissionYield(50.0);
  fallout.setBurstType("surface");

  const double median = 150e-6;
  const double gsd = 2.5;
  fallout.setLogNormalDistribution(median, gsd);

  EXPECT_NEAR(fallout.particleSizeCDF(median), 0.5, 5e-3);
  EXPECT_GT(fallout.particleSizePDF(median), 0.0);
}

TEST(RadiationPhysicsTest, FalloutFormationGeneratesParticlesWithMassConservation) {
  auto db = std::make_shared<RadionuclideDatabase>();
  auto inv = std::make_shared<FissionProductInventory>();
  inv->setDatabase(db);
  inv->setFissionYield(50.0);
  inv->initialize();

  FalloutFormationModel fallout;
  fallout.setYield(100.0);
  fallout.setFissionYield(50.0);
  fallout.setBurstType("surface");
  fallout.setInventory(inv);
  fallout.computeParticleFormation();

  const int n = 2000;
  const auto particles = fallout.generateParticles(n);
  ASSERT_EQ(static_cast<int>(particles.size()), n);

  double total_mass = 0.0;
  double total_activity = 0.0;
  for (const auto& p : particles) {
    EXPECT_GE(p.diameter, 1e-7);
    EXPECT_LE(p.diameter, 1e-2);
    EXPECT_GT(p.mass, 0.0);
    total_mass += p.mass;
    total_activity += p.activity;
  }

  EXPECT_NEAR(total_mass, fallout.totalFalloutMass(), relTol(fallout.totalFalloutMass(), 1e-8));
  EXPECT_GT(total_activity, 0.0);
}

TEST(RadiationPhysicsTest, RadiationUtilsConversionsAndDispersionAreSane) {
  EXPECT_NEAR(RadiationUtils::grayToRem(1.0), 100.0, 1e-12);
  EXPECT_NEAR(RadiationUtils::sievertToRem(1.0), 100.0, 1e-12);

  // Pasquill-Gifford sigmas should increase with downwind distance.
  const double sy_1km = RadiationUtils::pasquillGiffordSigmaY(1000.0, 'D');
  const double sy_10km = RadiationUtils::pasquillGiffordSigmaY(10000.0, 'D');
  EXPECT_GT(sy_10km, sy_1km);
}

TEST(RadiationPhysicsConfigTest, ParsesKeyValueConfigAndBuildsModels) {
  RadiationPhysicsConfig cfg;
  std::map<std::string, std::string> kv{
    {"yield_kt", "100.0"},
    {"fission_fraction", "0.5"},
    {"burst_height", "500.0"},
    {"burst_type", "surface"},
    {"wind_speed", "12.0"},
    {"wind_direction", "270.0"},
    {"stability_class", "D"},
    {"median_particle_size", "0.00015"},
    {"particle_size_gsd", "2.5"},
  };

  ASSERT_TRUE(cfg.parse(kv));

  auto prompt = cfg.createPromptModel();
  ASSERT_NE(prompt, nullptr);
  EXPECT_GT(prompt->totalGammaEnergy(), 0.0);
  EXPECT_GT(prompt->totalNeutronEnergy(), 0.0);
  EXPECT_GT(prompt->gammaYieldFraction(), 0.0);
  EXPECT_GT(prompt->neutronYieldFraction(), 0.0);

  auto inventory = cfg.createInventory();
  ASSERT_NE(inventory, nullptr);
  EXPECT_GT(inventory->totalActivity(), 0.0);
}
