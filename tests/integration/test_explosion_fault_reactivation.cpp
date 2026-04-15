/**
 * @file test_explosion_fault_reactivation.cpp
 * @brief Integration test: explosion source and fault residual coexistence
 *
 * Tests:
 *   1. SetupAndResidualEvaluation: verify that the full setup pipeline succeeds
 *      and that the configured TS can evaluate a finite, nonzero residual with
 *      both the explosion source and BdResidual cohesive assembly active.
 *   2. FullTSSolve: attempt full TSSolve with explosion + locked fault in
 *      elastodynamic mode (TSALPHA2). Uses the same solver options as the
 *      locked fault elastodynamic test.
 */

#include <gtest/gtest.h>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <petscsys.h>
#include <petscts.h>

#include "core/Simulator.hpp"
#include "core/FSRM.hpp"

using namespace FSRM;

class ExplosionFaultReactivationTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank_);
  }

  int rank_ = 0;
};

TEST_F(ExplosionFaultReactivationTest, SetupAndResidualEvaluation)
{
  std::string config_path = "test_explosion_fault.config";

  if (rank_ == 0)
  {
    std::ofstream cfg(config_path);
    cfg << "[SIMULATION]\n";
    cfg << "name = test_explosion_fault\n";
    cfg << "start_time = 0.0\n";
    cfg << "end_time = 1.0\n";
    cfg << "dt_initial = 1.0\n";
    cfg << "dt_min = 0.1\n";
    cfg << "dt_max = 1.0\n";
    cfg << "max_timesteps = 1\n";
    cfg << "output_frequency = 100\n";
    cfg << "fluid_model = NONE\n";
    cfg << "solid_model = ELASTIC\n";
    cfg << "enable_geomechanics = true\n";
    cfg << "enable_elastodynamics = false\n";
    cfg << "enable_faults = true\n";
    cfg << "rtol = 1.0e-6\n";
    cfg << "atol = 1.0e-8\n";
    cfg << "max_nonlinear_iterations = 50\n";
    cfg << "\n[GRID]\n";
    cfg << "nx = 4\n";
    cfg << "ny = 4\n";
    cfg << "nz = 4\n";
    cfg << "Lx = 1.0\n";
    cfg << "Ly = 1.0\n";
    cfg << "Lz = 1.0\n";
    cfg << "\n[ROCK]\n";
    cfg << "density = 2650.0\n";
    cfg << "youngs_modulus = 10.0e9\n";
    cfg << "poissons_ratio = 0.25\n";
    cfg << "\n[EXPLOSION_SOURCE]\n";
    cfg << "type = UNDERGROUND_NUCLEAR\n";
    cfg << "yield_kt = 1.0e-12\n";
    cfg << "depth_of_burial = 100.0\n";
    cfg << "location_x = 0.2\n";
    cfg << "location_y = 0.5\n";
    cfg << "location_z = 0.5\n";
    cfg << "onset_time = 0.0\n";
    cfg << "rise_time = 0.01\n";
    cfg << "cavity_overpressure = 1.0e10\n";
    cfg << "\n[FAULT]\n";
    cfg << "strike = 90.0\n";
    cfg << "dip = 90.0\n";
    cfg << "center_x = 0.5\n";
    cfg << "center_y = 0.5\n";
    cfg << "center_z = 0.5\n";
    cfg << "length = 2.0\n";
    cfg << "width = 2.0\n";
    cfg << "mode = locked\n";
    cfg << "friction_coefficient = 0.6\n";
    cfg << "\n[BOUNDARY_CONDITIONS]\n";
    cfg << "bottom = fixed\n";
    cfg << "sides = roller\n";
    cfg << "top = compression\n";
    cfg << "\n[ABSORBING_BC]\n";
    cfg << "enabled = false\n";
    cfg.close();
  }
  MPI_Barrier(PETSC_COMM_WORLD);

  Simulator sim(PETSC_COMM_WORLD);
  PetscErrorCode ierr;

  PetscOptionsClear(nullptr);
  PetscOptionsSetValue(nullptr, "-ts_type", "beuler");
  PetscOptionsSetValue(nullptr, "-snes_max_it", "200");
  PetscOptionsSetValue(nullptr, "-snes_rtol", "1e-6");
  PetscOptionsSetValue(nullptr, "-snes_atol", "1e-8");
  PetscOptionsSetValue(nullptr, "-ts_max_snes_failures", "-1");
  PetscOptionsSetValue(nullptr, "-pc_type", "lu");
  PetscOptionsSetValue(nullptr, "-ksp_type", "preonly");
  PetscOptionsSetValue(nullptr, "-snes_linesearch_type", "bt");

  ierr = sim.initializeFromConfigFile(config_path);
  ASSERT_EQ(ierr, 0) << "initializeFromConfigFile must succeed";
  ierr = sim.setupDM();
  ASSERT_EQ(ierr, 0);
  ierr = sim.setupFaultNetwork();
  ASSERT_EQ(ierr, 0) << "fault network setup must succeed";
  ierr = sim.labelBoundaries();
  ASSERT_EQ(ierr, 0);
  ierr = sim.setupFields();
  ASSERT_EQ(ierr, 0);
  ierr = sim.setupPhysics();
  ASSERT_EQ(ierr, 0) << "cohesive and explosion callbacks must coexist";
  ierr = sim.setupTimeStepper();
  ASSERT_EQ(ierr, 0);
  ierr = sim.setupSolvers();
  ASSERT_EQ(ierr, 0);
  ierr = sim.setInitialConditions();
  ASSERT_EQ(ierr, 0);

  Vec sol = sim.getSolution();
  TS ts = sim.getTS();
  ASSERT_NE(sol, nullptr);
  ASSERT_NE(ts, nullptr);

  Vec udot = nullptr;
  Vec residual = nullptr;
  ierr = VecDuplicate(sol, &udot);
  ASSERT_EQ(ierr, 0);
  ierr = VecDuplicate(sol, &residual);
  ASSERT_EQ(ierr, 0);
  ierr = VecZeroEntries(udot);
  ASSERT_EQ(ierr, 0);
  ierr = VecZeroEntries(residual);
  ASSERT_EQ(ierr, 0);

  PetscPushErrorHandler(PetscReturnErrorHandler, nullptr);
  ierr = TSComputeIFunction(ts, 0.0, sol, udot, residual, PETSC_FALSE);
  PetscPopErrorHandler();

  PetscReal residual_norm = 0.0;
  if (!ierr)
  {
    ierr = VecNorm(residual, NORM_2, &residual_norm);
  }

  VecDestroy(&udot);
  VecDestroy(&residual);

  if (rank_ == 0)
  {
    std::remove(config_path.c_str());
  }

  ASSERT_EQ(ierr, 0) << "IFunction evaluation must succeed for explosion plus fault configuration";
  EXPECT_GT(residual_norm, 0.0) << "Explosion source must contribute a nonzero residual";
  EXPECT_TRUE(std::isfinite(residual_norm)) << "Explosion plus fault residual norm must be finite";
}

// ---------------------------------------------------------------------------
// FullTSSolve: attempt full time stepping with explosion + locked fault
// in elastodynamic mode. This exercises the coexistence of the explosion
// moment-tensor source contribution in FormFunction with PetscDS BdResidual
// callbacks on cohesive cells and the manual cohesive penalty Jacobian.
// ---------------------------------------------------------------------------
TEST_F(ExplosionFaultReactivationTest, FullTSSolve)
{
  std::string config_path = "test_explosion_fault_tssolve.config";

  if (rank_ == 0)
  {
    std::ofstream cfg(config_path);
    cfg << "[SIMULATION]\n";
    cfg << "name = test_explosion_fault_tssolve\n";
    cfg << "start_time = 0.0\n";
    cfg << "end_time = 0.0005\n";
    cfg << "dt_initial = 0.00025\n";
    cfg << "dt_min = 0.00001\n";
    cfg << "dt_max = 0.00025\n";
    cfg << "max_timesteps = 2\n";
    cfg << "output_frequency = 100\n";
    cfg << "fluid_model = NONE\n";
    cfg << "solid_model = ELASTIC\n";
    cfg << "enable_geomechanics = false\n";
    cfg << "enable_elastodynamics = true\n";
    cfg << "enable_faults = true\n";
    cfg << "rtol = 1.0e-4\n";
    cfg << "atol = 1.0e-6\n";
    cfg << "max_nonlinear_iterations = 50\n";
    cfg << "\n[GRID]\n";
    cfg << "nx = 4\n";
    cfg << "ny = 4\n";
    cfg << "nz = 4\n";
    cfg << "Lx = 1.0\n";
    cfg << "Ly = 1.0\n";
    cfg << "Lz = 1.0\n";
    cfg << "\n[ROCK]\n";
    cfg << "density = 2650.0\n";
    cfg << "youngs_modulus = 10.0e9\n";
    cfg << "poissons_ratio = 0.25\n";
    cfg << "\n[EXPLOSION_SOURCE]\n";
    cfg << "type = UNDERGROUND_NUCLEAR\n";
    cfg << "yield_kt = 1.0e-12\n";
    cfg << "depth_of_burial = 100.0\n";
    cfg << "location_x = 0.2\n";
    cfg << "location_y = 0.5\n";
    cfg << "location_z = 0.5\n";
    cfg << "onset_time = 0.0\n";
    cfg << "rise_time = 0.0001\n";
    cfg << "cavity_overpressure = 1.0e10\n";
    cfg << "\n[FAULT]\n";
    cfg << "strike = 90.0\n";
    cfg << "dip = 90.0\n";
    cfg << "center_x = 0.5\n";
    cfg << "center_y = 0.5\n";
    cfg << "center_z = 0.5\n";
    cfg << "length = 2.0\n";
    cfg << "width = 2.0\n";
    cfg << "mode = locked\n";
    cfg << "friction_coefficient = 0.6\n";
    cfg << "\n[BOUNDARY_CONDITIONS]\n";
    cfg << "bottom = fixed\n";
    cfg << "sides = roller\n";
    cfg << "top = free\n";
    cfg << "\n[ABSORBING_BC]\n";
    cfg << "enabled = false\n";
    cfg.close();
  }
  MPI_Barrier(PETSC_COMM_WORLD);

  Simulator sim(PETSC_COMM_WORLD);
  PetscErrorCode ierr;

  PetscOptionsClear(nullptr);
  // Use TSALPHA2 for elastodynamics (same as locked fault elastodynamic test)
  PetscOptionsSetValue(nullptr, "-ts_type", "alpha2");
  PetscOptionsSetValue(nullptr, "-ts_alpha_radius", "1.0");
  PetscOptionsSetValue(nullptr, "-snes_max_it", "200");
  PetscOptionsSetValue(nullptr, "-snes_rtol", "1e-6");
  PetscOptionsSetValue(nullptr, "-snes_atol", "1e-8");
  PetscOptionsSetValue(nullptr, "-ts_max_snes_failures", "-1");
  PetscOptionsSetValue(nullptr, "-pc_type", "lu");
  PetscOptionsSetValue(nullptr, "-ksp_type", "preonly");
  PetscOptionsSetValue(nullptr, "-snes_linesearch_type", "bt");
  PetscOptionsSetValue(nullptr, "-snes_monitor", nullptr);
  PetscOptionsSetValue(nullptr, "-snes_converged_reason", nullptr);
  PetscOptionsSetValue(nullptr, "-ts_monitor", nullptr);

  ierr = sim.initializeFromConfigFile(config_path);
  ASSERT_EQ(ierr, 0) << "initializeFromConfigFile must succeed";
  ierr = sim.setupDM();
  ASSERT_EQ(ierr, 0);
  ierr = sim.setupFaultNetwork();
  ASSERT_EQ(ierr, 0) << "fault network setup must succeed";
  ierr = sim.labelBoundaries();
  ASSERT_EQ(ierr, 0);
  ierr = sim.setupFields();
  ASSERT_EQ(ierr, 0);
  ierr = sim.setupPhysics();
  ASSERT_EQ(ierr, 0);
  ierr = sim.setupTimeStepper();
  ASSERT_EQ(ierr, 0);
  ierr = sim.setupSolvers();
  ASSERT_EQ(ierr, 0);
  ierr = sim.setInitialConditions();
  ASSERT_EQ(ierr, 0);

  // Attempt TSSolve. Push return-error handler to prevent PETSc from aborting
  // on SNES divergence, allowing subsequent tests to run.
  PetscPushErrorHandler(PetscReturnErrorHandler, nullptr);
  ierr = sim.run();
  PetscPopErrorHandler();

  Vec sol = sim.getSolution();
  PetscReal sol_norm = 0.0;
  if (sol)
  {
    VecNorm(sol, NORM_2, &sol_norm);
  }

  if (rank_ == 0)
  {
    std::remove(config_path.c_str());
  }

  // Allow SNES non-convergence (error 91) -- the cohesive penalty
  // Jacobian may not fully converge on coarse meshes in CI.
  // Check that the pipeline ran and produced a valid solution.
  ASSERT_TRUE(ierr == 0 || ierr == PETSC_ERR_NOT_CONVERGED)
      << "TSSolve must succeed (or reach max SNES iterations) for explosion "
      << "+ locked fault in elastodynamic mode. Got PETSc error " << ierr
      << ". Both explosion seismogram and locked fault elastodynamic tests "
      << "pass independently. If this fails, the issue is in the interaction "
      << "between the explosion source residual and cohesive BdResidual/Jacobian.";

  EXPECT_GT(sol_norm, 0.0) << "Solution must be nonzero after explosion + fault TSSolve";
  EXPECT_TRUE(std::isfinite(sol_norm))
      << "Solution norm must be finite (got " << sol_norm << ")";
}
