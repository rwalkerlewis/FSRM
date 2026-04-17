/**
 * @file test_proppant_transport.cpp
 * @brief Phase 6 proppant transport integration tests
 *
 * Validates Stokes settling, bridging criterion, pack aperture,
 * and mass conservation for proppant in hydraulic fractures.
 */

#include <gtest/gtest.h>

#include <cmath>

#include "numerics/PetscFEHydrofrac.hpp"

using namespace FSRM;

class ProppantTransportTest : public ::testing::Test
{
};

TEST_F(ProppantTransportTest, StokesSettlingVelocityPositiveForDenseParticles)
{
  // 20/40 mesh sand: d = 0.5 mm, rho_p = 2650 kg/m^3
  // Slickwater: rho_f = 1000 kg/m^3, mu = 1e-3 Pa.s
  const PetscReal d = 5.0e-4;
  const PetscReal rho_p = 2650.0;
  const PetscReal rho_f = 1000.0;
  const PetscReal mu = 1.0e-3;

  const PetscReal v_s = PetscFEHydrofrac::stokesSettlingVelocity(d, rho_p, rho_f, mu);

  // v_s = d^2 * (rho_p - rho_f) * g / (18*mu)
  //     = (5e-4)^2 * 1650 * 9.81 / (18 * 1e-3)
  //     = 2.5e-7 * 1650 * 9.81 / 0.018
  //     = 0.2247 m/s
  EXPECT_GT(v_s, 0.0);
  EXPECT_NEAR(v_s, 0.2247, 0.01);
}

TEST_F(ProppantTransportTest, SettlingVelocityScalesWithDiameterSquared)
{
  const PetscReal rho_p = 2650.0;
  const PetscReal rho_f = 1000.0;
  const PetscReal mu = 1.0e-3;

  const PetscReal v1 = PetscFEHydrofrac::stokesSettlingVelocity(1e-4, rho_p, rho_f, mu);
  const PetscReal v2 = PetscFEHydrofrac::stokesSettlingVelocity(2e-4, rho_p, rho_f, mu);

  const PetscReal ratio = v2 / v1;
  EXPECT_NEAR(ratio, 4.0, 1e-10);
}

TEST_F(ProppantTransportTest, SettlingVelocityNegativeForBuoyantParticles)
{
  // Hollow ceramic: rho_p = 800 kg/m^3 in brine rho_f = 1100 kg/m^3
  const PetscReal v_s = PetscFEHydrofrac::stokesSettlingVelocity(
      5e-4, 800.0, 1100.0, 1.0e-3);
  EXPECT_LT(v_s, 0.0);
}

TEST_F(ProppantTransportTest, BridgingOccursBelowCriticalRatio)
{
  const PetscReal d = 5.0e-4;
  const PetscReal bridging_ratio = 3.0;

  // Narrow aperture: w/d = 2.0 < 3.0 -> bridging
  EXPECT_EQ(PetscFEHydrofrac::proppantBridging(1.0e-3, d, bridging_ratio), PETSC_TRUE);

  // Wide aperture: w/d = 10.0 > 3.0 -> no bridging
  EXPECT_EQ(PetscFEHydrofrac::proppantBridging(5.0e-3, d, bridging_ratio), PETSC_FALSE);

  // At exactly the ratio: w/d = 3.0 -> no bridging (not strictly less)
  EXPECT_EQ(PetscFEHydrofrac::proppantBridging(1.5e-3, d, bridging_ratio), PETSC_FALSE);
}

TEST_F(ProppantTransportTest, PackMinApertureEqualsParticleDiameter)
{
  const PetscReal d = 5.0e-4;
  const PetscReal phi = 0.38;

  const PetscReal w_min = PetscFEHydrofrac::proppantPackMinAperture(d, phi);
  EXPECT_NEAR(w_min, d, 1e-12);
}

TEST_F(ProppantTransportTest, MassBalanceRatioNearUnityForConservation)
{
  // Perfect conservation scenario
  const PetscReal c_in = 100.0;       // 100 kg/m^3 slurry
  const PetscReal Q = 0.01;           // 10 L/s = 0.01 m^3/s
  const PetscReal t_inj = 600.0;      // 10 minutes
  const PetscReal V_f = 6.0;          // 6 m^3 fracture volume
  const PetscReal c_avg = 100.0;      // all mass retained

  const PetscReal balance = PetscFEHydrofrac::proppantMassBalance(
      c_in, Q, t_inj, V_f, c_avg);
  EXPECT_NEAR(balance, 1.0, 1e-10);
}

TEST_F(ProppantTransportTest, MassBalanceDetectsLeakage)
{
  // Some mass leaked off
  const PetscReal c_in = 100.0;
  const PetscReal Q = 0.01;
  const PetscReal t_inj = 600.0;
  const PetscReal V_f = 6.0;
  const PetscReal c_avg = 80.0;       // 80% retained

  const PetscReal balance = PetscFEHydrofrac::proppantMassBalance(
      c_in, Q, t_inj, V_f, c_avg);
  EXPECT_LT(balance, 1.0);
  EXPECT_GT(balance, 0.0);
}

TEST_F(ProppantTransportTest, InvalidInputsReturnZeroOrNegative)
{
  EXPECT_EQ(PetscFEHydrofrac::stokesSettlingVelocity(0.0, 2650.0, 1000.0, 1.0e-3), 0.0);
  EXPECT_EQ(PetscFEHydrofrac::stokesSettlingVelocity(1.0e-4, 2650.0, 1000.0, 0.0), 0.0);
  EXPECT_EQ(PetscFEHydrofrac::proppantPackMinAperture(0.0, 0.38), 0.0);
  EXPECT_LT(PetscFEHydrofrac::proppantMassBalance(100.0, 0.0, 600.0, 6.0, 100.0), 0.0);
}
