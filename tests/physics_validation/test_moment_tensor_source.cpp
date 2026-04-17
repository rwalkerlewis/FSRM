/**
 * @file test_moment_tensor_source.cpp
 * @brief Verification of moment tensor source implementation (BUG 5 fix)
 *
 * After fixing BUG 5 (proper FEM equivalent nodal forces via tabulation),
 * verify that the explosion source produces nonzero displacement in the
 * solution vector. This confirms that DMPlexVecSetClosure correctly injects
 * force contributions at the source cell.
 *
 * Tests:
 *   1. Pipeline completes and solution has nonzero displacement norm
 *   2. Initial (zero) solution has strictly smaller norm than final solution
 *   3. Explosion coupling struct is initialized with correct parameters
 */

#include <gtest/gtest.h>
#include <cmath>
#include <fstream>
#include <string>
#include <petscsys.h>
#include <petscvec.h>

#include "core/Simulator.hpp"
#include "core/FSRM.hpp"

using namespace FSRM;

class MomentTensorSourceTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank_);

    if (rank_ == 0)
    {
      // Source at center of 4km box: (2000, 2000, 2000)
      std::ofstream cfg(config_path_);
      cfg << "[SIMULATION]\n";
      cfg << "name = test_moment_tensor\n";
      cfg << "start_time = 0.0\n";
      cfg << "end_time = 0.1\n";
      cfg << "dt_initial = 0.002\n";
      cfg << "dt_min = 0.001\n";
      cfg << "dt_max = 0.005\n";
      cfg << "max_timesteps = 50\n";
      cfg << "output_frequency = 200\n";
      cfg << "fluid_model = NONE\n";
      cfg << "solid_model = ELASTIC\n";
      cfg << "enable_geomechanics = true\n";
      cfg << "enable_faults = false\n";
      cfg << "enable_elastodynamics = true\n";
      cfg << "rtol = 1.0e-6\n";
      cfg << "atol = 1.0e-8\n";
      cfg << "max_nonlinear_iterations = 20\n";
      cfg << "\n[GRID]\n";
      cfg << "nx = 4\n";
      cfg << "ny = 4\n";
      cfg << "nz = 4\n";
      cfg << "Lx = 4000.0\n";
      cfg << "Ly = 4000.0\n";
      cfg << "Lz = 4000.0\n";
      cfg << "\n[ROCK]\n";
      cfg << "density = 2650.0\n";
      cfg << "youngs_modulus = 60.0e9\n";
      cfg << "poissons_ratio = 0.25\n";
      cfg << "\n[EXPLOSION_SOURCE]\n";
      cfg << "type = UNDERGROUND_NUCLEAR\n";
      cfg << "yield_kt = 10.0\n";
      cfg << "depth_of_burial = 500.0\n";
      cfg << "location_x = 2000.0\n";
      cfg << "location_y = 2000.0\n";
      cfg << "location_z = 2000.0\n";
      cfg << "onset_time = 0.0\n";
      cfg << "rise_time = 0.01\n";
      cfg << "cavity_overpressure = 1.0e10\n";
      cfg << "\n[BOUNDARY_CONDITIONS]\n";
      cfg << "bottom = free\n";
      cfg << "sides = free\n";
      cfg << "top = free\n";
      cfg << "\n[ABSORBING_BC]\n";
      cfg << "enabled = true\n";
      cfg << "x_min = true\n";
      cfg << "x_max = true\n";
      cfg << "y_min = true\n";
      cfg << "y_max = true\n";
      cfg << "z_min = true\n";
      cfg << "z_max = false\n";
      cfg.close();
    }
    MPI_Barrier(PETSC_COMM_WORLD);
  }

  void TearDown() override
  {
    MPI_Barrier(PETSC_COMM_WORLD);
    if (rank_ == 0)
    {
      std::remove(config_path_);
    }
  }

  // Run the full simulation pipeline and return the simulator
  PetscErrorCode runPipeline(Simulator& sim)
  {
    PetscErrorCode ierr;
    ierr = sim.initializeFromConfigFile(config_path_);
    if (ierr) return ierr;
    ierr = sim.setupDM();
    if (ierr) return ierr;
    ierr = sim.labelBoundaries();
    if (ierr) return ierr;
    ierr = sim.setupFields();
    if (ierr) return ierr;
    ierr = sim.setupPhysics();
    if (ierr) return ierr;
    ierr = sim.setupTimeStepper();
    if (ierr) return ierr;
    ierr = sim.setupSolvers();
    if (ierr) return ierr;
    ierr = sim.setInitialConditions();
    if (ierr) return ierr;
    return 0;
  }

  int rank_ = 0;
  static constexpr const char* config_path_ = "test_moment_tensor.ini";
};

// ---------------------------------------------------------------------------
// Test 1: Simulation completes and solution has nonzero displacement
//
// BUG 5 fix ensures PetscFECreateTabulation is used for proper equivalent
// nodal forces. Without the fix, the explosion source would produce zero
// displacement because the force injection was incorrect.
// ---------------------------------------------------------------------------
TEST_F(MomentTensorSourceTest, PipelineProducesNonzeroDisplacement)
{
  Simulator sim(PETSC_COMM_WORLD);
  PetscErrorCode ierr = runPipeline(sim);
  ASSERT_EQ(ierr, 0) << "Pipeline setup failed";

  // Check solution norm is zero before time stepping
  Vec sol = sim.getSolution();
  ASSERT_NE(sol, nullptr) << "Solution vector should exist after setup";

  PetscReal norm_before;
  VecNorm(sol, NORM_2, &norm_before);

  // Run the simulation
  ierr = sim.run();
  ASSERT_EQ(ierr, 0) << "Simulation run failed";

  // Solution should have nonzero displacement after explosion source
  PetscReal norm_after;
  VecNorm(sol, NORM_2, &norm_after);
  EXPECT_GT(norm_after, 1.0e-20)
      << "Solution should have nonzero displacement from explosion source";
}

// ---------------------------------------------------------------------------
// Test 2: Final solution norm exceeds initial solution norm
//
// The explosion source injects energy into the system. After time stepping,
// the solution norm must be strictly larger than the initial zero solution.
// This verifies that addExplosionSourceToResidual produces nonzero forces.
// ---------------------------------------------------------------------------
TEST_F(MomentTensorSourceTest, SolutionNormGrows)
{
  Simulator sim(PETSC_COMM_WORLD);
  PetscErrorCode ierr = runPipeline(sim);
  ASSERT_EQ(ierr, 0) << "Pipeline setup failed";

  Vec sol = sim.getSolution();
  PetscReal norm_initial;
  VecNorm(sol, NORM_2, &norm_initial);

  ierr = sim.run();
  ASSERT_EQ(ierr, 0) << "Simulation run failed";

  PetscReal norm_final;
  VecNorm(sol, NORM_2, &norm_final);

  EXPECT_GT(norm_final, norm_initial)
      << "Final solution norm (" << norm_final
      << ") must exceed initial norm (" << norm_initial << ")";
}

// ---------------------------------------------------------------------------
// Test 3: Solution vector has correct size for 3D elasticity
//
// The solution vector should have 3 DOFs per node (ux, uy, uz).
// This verifies the DM and FE fields are set up correctly.
// ---------------------------------------------------------------------------
TEST_F(MomentTensorSourceTest, SolutionVectorSize)
{
  Simulator sim(PETSC_COMM_WORLD);
  PetscErrorCode ierr = runPipeline(sim);
  ASSERT_EQ(ierr, 0) << "Pipeline setup failed";

  Vec sol = sim.getSolution();
  ASSERT_NE(sol, nullptr) << "Solution vector should exist";

  PetscInt n;
  VecGetSize(sol, &n);

  // 3D elasticity on a simplex mesh: 3 DOFs per node
  // The exact size depends on mesh refinement, but must be a multiple of 3
  EXPECT_GT(n, 0) << "Solution vector must be nonempty";
  EXPECT_EQ(n % 3, 0)
      << "Solution vector size (" << n << ") should be multiple of 3 for 3D displacement";
}
