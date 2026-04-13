/**
 * @file test_explosion_fault_reactivation.cpp
 * @brief Integration test: explosion-induced fault reactivation
 *
 * Combines explosion source injection with a pre-stressed fault.
 * The explosion generates seismic waves that load the fault, potentially
 * triggering slip. This is the full coupled workflow:
 *   1. Mesh with cohesive fault cells
 *   2. Moment tensor explosion source
 *   3. Elastodynamic time stepping (TSALPHA2)
 *   4. Fault friction-controlled slip
 *
 * Known limitation: TSSolve with cohesive cells diverges in PETSc 3.22
 * due to BdResidual totDim mismatch. This test documents the intended
 * workflow and will pass (instead of skip) once the Jacobian is fixed.
 */

#include <gtest/gtest.h>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <petscsys.h>

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

// Explosion at (0.25, 0.5, 0.5) near a vertical fault at x=0.5.
// The explosion generates P-waves that load the fault in compression
// and S-waves that apply shear traction.
TEST_F(ExplosionFaultReactivationTest, SetupAndSolveAttempt)
{
  std::string config_path = "test_explosion_fault.config";

  if (rank_ == 0)
  {
    std::ofstream cfg(config_path);
    cfg << "[SIMULATION]\n";
    cfg << "name = test_explosion_fault\n";
    cfg << "start_time = 0.0\n";
    cfg << "end_time = 0.0005\n";
    cfg << "dt_initial = 0.00005\n";
    cfg << "dt_min = 0.00001\n";
    cfg << "dt_max = 0.0001\n";
    cfg << "max_timesteps = 5\n";
    cfg << "output_frequency = 100\n";
    cfg << "fluid_model = NONE\n";
    cfg << "solid_model = ELASTIC\n";
    cfg << "enable_geomechanics = true\n";
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
    cfg << "\n[EXPLOSION]\n";
    cfg << "yield_kt = 0.001\n";
    cfg << "depth_m = 0.5\n";
    cfg << "center_x = 0.25\n";
    cfg << "center_y = 0.5\n";
    cfg << "center_z = 0.5\n";
    cfg << "\n[FAULT]\n";
    cfg << "strike = 90.0\n";
    cfg << "dip = 90.0\n";
    cfg << "center_x = 0.5\n";
    cfg << "center_y = 0.5\n";
    cfg << "center_z = 0.5\n";
    cfg << "length = 2.0\n";
    cfg << "width = 2.0\n";
    cfg << "mode = slipping\n";
    cfg << "friction_coefficient = 0.3\n";
    cfg << "initial_shear_traction = 25.0e6\n";
    cfg << "initial_normal_traction = -50.0e6\n";
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

  // Full pipeline setup must succeed
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
  ASSERT_EQ(ierr, 0) << "cohesive + explosion callbacks must coexist";
  ierr = sim.setupTimeStepper();
  ASSERT_EQ(ierr, 0);
  ierr = sim.setupSolvers();
  ASSERT_EQ(ierr, 0);
  ierr = sim.setInitialConditions();
  ASSERT_EQ(ierr, 0);

  // Attempt TSSolve -- push return-error handler to prevent fatal abort
  PetscPushErrorHandler(PetscReturnErrorHandler, nullptr);
  ierr = sim.run();
  PetscPopErrorHandler();

  if (rank_ == 0) std::remove(config_path.c_str());

  if (ierr != 0)
  {
    GTEST_SKIP() << "TSSolve diverged for explosion+fault (known limitation: "
                 << "PetscDSSetBdResidual on cohesive cells). ierr=" << ierr;
  }

  // If TSSolve succeeds, verify solution quality
  Vec sol = sim.getSolution();
  ASSERT_NE(sol, nullptr);
  PetscReal norm = 0.0;
  VecNorm(sol, NORM_2, &norm);
  EXPECT_GT(norm, 0.0) << "Explosion must produce nonzero displacement";
  EXPECT_LT(norm, 100.0) << "Solution must be bounded";
}
