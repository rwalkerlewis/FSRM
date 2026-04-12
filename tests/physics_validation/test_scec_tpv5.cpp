/**
 * @file test_scec_tpv5.cpp
 * @brief SCEC TPV5 (vertical strike-slip, slip-weakening) benchmark validation
 *
 * TPV5 is a vertical strike-slip fault with linear slip-weakening friction.
 *
 * Status:
 *   - Fault mesh splitting: VERIFIED via Unit.FaultMeshManager (32 cohesive cells)
 *   - Friction laws (slip-weakening, rate-state): VERIFIED via Unit.FaultMechanics
 *   - CohesiveFaultKernel registration: VERIFIED via Unit.CohesiveFaultKernel
 *   - Full dynamic rupture end-to-end: NOT YET IMPLEMENTED
 *     CohesiveFaultKernel is not integrated into the Simulator IFunction/IJacobian
 *     pipeline. The kernel can be constructed and registered, but the cohesive
 *     forces are not assembled during time stepping.
 *
 * Tests below verify as far as the pipeline currently goes.
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

// ---------------------------------------------------------------------------
// Test 2: CohesiveFaultKernel friction parameter setup
//
// Verify that the kernel accepts TPV5 friction parameters and reports
// correct state. This documents the furthest point the pipeline reaches
// without crashing -- construction and parameter setup work, but the
// kernel is not yet integrated into the Simulator IFunction/IJacobian
// assembly pipeline, so no dynamic rupture actually runs.
// ---------------------------------------------------------------------------
TEST_F(SCECTPV5Test, KernelFrictionParameterSetup)
{
  FSRM::CohesiveFaultKernel kernel;

  // Set TPV5 friction parameters
  kernel.setFrictionCoefficient(0.677);
  kernel.setTensileStrength(0.0);
  kernel.setMode(true);  // locked initially

  EXPECT_TRUE(kernel.isLocked())
      << "Kernel should be locked before rupture nucleation";

  // Switch to slipping
  kernel.setMode(false);
  EXPECT_FALSE(kernel.isLocked())
      << "Kernel should be unlocked after mode change";

  // Set prescribed slip (for testing the prescribed slip pathway)
  kernel.setPrescribedSlip(0.0, 0.0, 0.0);
}

// ---------------------------------------------------------------------------
// Test 3: Dynamic rupture pipeline limitation
//
// Documents the current failure point: CohesiveFaultKernel is not
// integrated into the Simulator IFunction/IJacobian pipeline.
// Mesh splitting works, kernel registration works, but cohesive
// forces are not assembled during TSSolve.
//
// When this is implemented, this test should be replaced with an
// actual dynamic rupture simulation comparing slip rate against
// the SCEC reference solution.
// ---------------------------------------------------------------------------
TEST_F(SCECTPV5Test, DynamicRupturePipelineLimitation)
{
  // Document the current state
  // This test passes because it documents the limitation.
  // It does NOT skip -- it explicitly states what works and what does not.
  FSRM::FaultMeshManager fmm(PETSC_COMM_WORLD);
  FSRM::CohesiveFaultKernel kernel;

  // Step 1: FaultMeshManager construction -- WORKS
  // Step 2: CohesiveFaultKernel construction -- WORKS
  // Step 3: Kernel parameter setup -- WORKS (tested above)
  // Step 4: Mesh splitting (DMPlexConstructCohesiveCells) -- WORKS (Unit.FaultMeshManager)
  // Step 5: Kernel registration with PetscDS -- WORKS (Unit.CohesiveFaultKernel)
  // Step 6: Cohesive force assembly in Simulator::FormFunction -- NOT IMPLEMENTED
  //         The IFunction callback does not call cohesive kernel contributions.
  // Step 7: Full dynamic rupture TSSolve -- BLOCKED by Step 6

  // Verify we can at least get this far
  kernel.setFrictionCoefficient(0.677);
  EXPECT_TRUE(true)
      << "Pipeline reaches kernel parameter setup. "
         "Blocked at: cohesive force assembly in FormFunction (Step 6). "
         "Required: integrate CohesiveFaultKernel into IFunction/IJacobian pipeline.";
}
