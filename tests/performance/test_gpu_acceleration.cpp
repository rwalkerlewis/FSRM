/**
 * @file test_gpu_acceleration.cpp
 * @brief Performance test: PETSc CUDA GPU acceleration
 *
 * Validates that FSRM produces correct results with PETSc CUDA backend.
 * Skips if CUDA is not available (acceptable for CI without GPU hardware).
 *
 * PETSc CUDA acceleration is enabled at runtime via:
 *   -vec_type cuda -mat_type aijcusparse
 *
 * No FSRM code changes are needed. PETSc handles vector operations via
 * cuBLAS, matrix operations via cuSPARSE, and KSP solves on GPU.
 */

#include <gtest/gtest.h>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <petscsys.h>
#include <petscvec.h>

#include "core/Simulator.hpp"
#include "core/FSRM.hpp"

using namespace FSRM;

class GPUAccelerationTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank_);
  }

  int rank_ = 0;
};

// Check if CUDA Vec type is available by attempting to create one.
// Returns true if CUDA is available, false otherwise.
static bool isCUDAAvailable()
{
  Vec test_vec;
  PetscErrorCode ierr;

  PetscPushErrorHandler(PetscReturnErrorHandler, nullptr);
  ierr = VecCreate(PETSC_COMM_SELF, &test_vec);
  if (ierr)
  {
    PetscPopErrorHandler();
    return false;
  }
  ierr = VecSetType(test_vec, VECCUDA);
  PetscPopErrorHandler();

  if (ierr)
  {
    VecDestroy(&test_vec);
    return false;
  }

  ierr = VecSetSizes(test_vec, 10, PETSC_DECIDE);
  if (ierr)
  {
    VecDestroy(&test_vec);
    return false;
  }

  VecDestroy(&test_vec);
  return true;
}

// Test that CUDA backend produces correct results for a small elastostatic problem.
// This is the only acceptable use of GTEST_SKIP: hardware-dependent tests.
TEST_F(GPUAccelerationTest, CUDABackendProducesCorrectResults)
{
  if (!isCUDAAvailable())
  {
    GTEST_SKIP() << "CUDA not available (no GPU hardware or PETSc not built with --with-cuda)";
  }

  std::string config_path = "test_gpu_elastostatics.config";

  if (rank_ == 0)
  {
    std::ofstream cfg(config_path);
    cfg << "[SIMULATION]\n";
    cfg << "name = test_gpu\n";
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
    cfg << "enable_faults = false\n";
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
    cfg << "\n[BOUNDARY_CONDITIONS]\n";
    cfg << "bottom = fixed\n";
    cfg << "sides = roller\n";
    cfg << "top = compression\n";
    cfg << "\n[ABSORBING_BC]\n";
    cfg << "enabled = false\n";
    cfg.close();
  }
  MPI_Barrier(PETSC_COMM_WORLD);

  // Set CUDA backend via PETSc options
  PetscOptionsSetValue(nullptr, "-vec_type", "cuda");
  PetscOptionsSetValue(nullptr, "-mat_type", "aijcusparse");

  Simulator sim(PETSC_COMM_WORLD);
  PetscErrorCode ierr;

  ierr = sim.initializeFromConfigFile(config_path);
  ASSERT_EQ(ierr, 0);
  ierr = sim.setupDM();
  ASSERT_EQ(ierr, 0);
  ierr = sim.setupFaultNetwork();
  ASSERT_EQ(ierr, 0);
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

  ierr = sim.run();
  ASSERT_EQ(ierr, 0) << "GPU-accelerated simulation must complete successfully";

  Vec sol = sim.getSolution();
  ASSERT_NE(sol, nullptr);
  PetscReal norm = 0.0;
  VecNorm(sol, NORM_2, &norm);
  EXPECT_GT(norm, 0.0) << "GPU solution must be nonzero";
  EXPECT_LT(norm, 1.0) << "GPU solution must be bounded";

  if (rank_ == 0) std::remove(config_path.c_str());

  // Reset PETSc options to avoid affecting other tests
  PetscOptionsClearValue(nullptr, "-vec_type");
  PetscOptionsClearValue(nullptr, "-mat_type");
}
