/**
 * @file test_scec_tpv5.cpp
 * @brief SCEC TPV5 (vertical strike-slip, rate-state) benchmark validation
 *
 * TPV5 is a 2D strike-slip fault with linear slip-weakening friction.
 *
 * Status:
 *   - Fault mesh splitting: VERIFIED via Unit.FaultMeshManager (32 cohesive cells)
 *   - Friction laws (slip-weakening, rate-state): VERIFIED via Unit.FaultMechanics
 *   - CohesiveFaultKernel registration: VERIFIED via Unit.CohesiveFaultKernel
 *   - Full dynamic rupture end-to-end: NOT YET WIRED (requires CohesiveFaultKernel
 *     integration into the Simulator IFunction/IJacobian pipeline)
 *   - Slip rate comparison against SCEC reference: PENDING (requires end-to-end run)
 *
 * When full dynamic rupture is implemented, the BenchmarkSlipRateWithinFivePercent
 * test should be activated to compare slip rate time histories at on-fault stations
 * against the Day et al. reference solution with 5% tolerance.
 */

#include <gtest/gtest.h>
#include <petscsys.h>
#include <cmath>

#include "numerics/FaultMeshManager.hpp"
#include "physics/CohesiveFaultKernel.hpp"

class SCECTPV5Test : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    }

    int rank = 0;
};

// ---------------------------------------------------------------------------
// Test 1: TPV5 fault setup components
//
// Verify the building blocks needed for SCEC TPV5:
//   - FaultMeshManager can locate fault planes
//   - CohesiveFaultKernel can be constructed with slip-weakening parameters
//   - TPV5 material parameters are physically reasonable
// ---------------------------------------------------------------------------
TEST_F(SCECTPV5Test, FaultInfrastructureReady)
{
  // TPV5 material parameters (SCEC specification)
  const double mu_s = 0.677;     // Static friction coefficient
  const double mu_d = 0.525;     // Dynamic friction coefficient
  const double Dc = 0.40;        // Critical slip distance (m)
  const double rho = 2670.0;     // Density (kg/m^3)
  const double vs = 3464.0;      // S-wave velocity (m/s)

  // Verify the shear modulus matches specification
  double G = rho * vs * vs;      // Should be ~32 GPa
  EXPECT_GT(G, 30.0e9) << "Shear modulus too low for TPV5";
  EXPECT_LT(G, 35.0e9) << "Shear modulus too high for TPV5";

  // Verify friction drop is positive (mu_s > mu_d)
  EXPECT_GT(mu_s, mu_d)
      << "Static friction must exceed dynamic friction for slip weakening";

  // Verify Dc is positive
  EXPECT_GT(Dc, 0.0) << "Critical slip distance must be positive";

  // CohesiveFaultKernel can be constructed
  FSRM::CohesiveFaultKernel kernel;

  // FaultMeshManager can be constructed (requires MPI communicator)
  FSRM::FaultMeshManager fmm(PETSC_COMM_WORLD);

  // Both constructions succeed (no crash) -- this is the minimum viability
  SUCCEED() << "Fault infrastructure components can be instantiated";
}

TEST_F(SCECTPV5Test, BenchmarkSlipRateWithinFivePercent) {
    GTEST_SKIP() << "Full SCEC TPV5 simulation not wired in CI; future test will compare slip "
                    "rate at recording stations to the SCEC reference within 5%. "
                    "Prerequisite: CohesiveFaultKernel integration into Simulator IFunction pipeline.";
}

TEST_F(SCECTPV5Test, DocumentedVerificationTargets) {
    GTEST_SKIP()
        << "TPV5 validation targets: (1) spontaneous rupture on a vertical strike-slip fault "
           "with rate-and-state friction, (2) slip-rate time series at prescribed fault "
           "receivers vs. published benchmark, (3) pass criterion max|DeltaV|/V_ref < 5% on peaks "
           "or equivalent integrated metric. Requires simulation driver + reference data.";
}
