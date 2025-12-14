/**
 * @file test_atmospheric_explosion.cpp
 * @brief Unit tests for atmospheric nuclear detonation modeling components
 */

#include <gtest/gtest.h>

#include "AtmosphericExplosion.hpp"

#include <cmath>
#include <map>
#include <memory>
#include <string>

using namespace FSRM;

namespace {
double relTol(double x, double rel = 1e-6, double abs_min = 1e-9) {
  return std::max(abs_min, std::abs(x) * rel);
}
}  // namespace

TEST(AtmosphericExplosionTest, StandardAtmosphereHasExpectedSeaLevelValues) {
  ExtendedAtmosphericModel atm;
  atm.setStandardAtmosphere();

  EXPECT_NEAR(atm.temperature(0.0), AtmosphericExplosionConstants::SEA_LEVEL_TEMPERATURE, 1e-6);
  EXPECT_NEAR(atm.pressure(0.0), AtmosphericExplosionConstants::SEA_LEVEL_PRESSURE, 1e-2);

  // Density is derived from ideal gas and will be close to the constant.
  const double rho0 = atm.density(0.0);
  EXPECT_NEAR(rho0, AtmosphericExplosionConstants::SEA_LEVEL_DENSITY, 1e-3);

  // Monotone decrease with height (basic sanity).
  EXPECT_LT(atm.pressure(1000.0), atm.pressure(0.0));
  EXPECT_LT(atm.density(1000.0), atm.density(0.0));
}

TEST(AtmosphericExplosionTest, MushroomCloudInitialRadiusAndTimesFollowScaling) {
  auto atm = std::make_shared<ExtendedAtmosphericModel>();
  atm->setStandardAtmosphere();

  MushroomCloudModel cloud;
  cloud.setAtmosphere(atm);
  cloud.setYield(500.0);
  cloud.setBurstHeight(500.0);
  cloud.initialize();

  const double expected_Rmax = 66.0 * std::pow(500.0, 0.4);
  EXPECT_NEAR(cloud.cloudRadius(), 0.1 * expected_Rmax, relTol(0.1 * expected_Rmax));

  const double expected_tform = 0.001 * std::pow(500.0, 0.4);
  EXPECT_NEAR(cloud.fireballMaxTime(), expected_tform, relTol(expected_tform, 1e-12));

  const double expected_Hstab = (10.0 + 3.5 * std::pow(500.0, 0.4)) * 1000.0;
  EXPECT_NEAR(cloud.stabilizationHeight(), expected_Hstab, relTol(expected_Hstab, 1e-12));
}

TEST(AtmosphericExplosionConfigTest, ParsesKeyValueConfigAndCreatesSolver) {
  AtmosphericExplosionConfig cfg;
  std::map<std::string, std::string> kv{
    {"yield_kt", "500.0"},
    {"burst_height", "500.0"},
    {"location_x", "50000.0"},
    {"location_y", "50000.0"},
    {"location_z", "500.0"},
    {"fission_fraction", "0.5"},
    {"enable_instabilities", "true"},
    {"enable_mushroom_cloud", "true"},
    {"enable_debris_transport", "true"},
  };

  ASSERT_TRUE(cfg.parse(kv));
  auto solver = cfg.createSolver();
  ASSERT_NE(solver, nullptr);

  solver->initialize();
  auto cloud = solver->getCloudModel();
  ASSERT_NE(cloud, nullptr);

  const double expected_Rmax = 66.0 * std::pow(cfg.getYield(), 0.4);
  EXPECT_NEAR(cloud->cloudRadius(), 0.1 * expected_Rmax, relTol(0.1 * expected_Rmax));
}
