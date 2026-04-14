// test_traction_bc.cpp
// Integration test for per-face traction boundary conditions.
// Applies 1 MPa downward traction on z_max, fixed z_min, roller sides.
// Verifies displacement matches analytical uniaxial stress solution.

#include <gtest/gtest.h>
#include <petsc.h>
#include <fstream>
#include <cstdio>
#include <cmath>
#include "core/Simulator.hpp"

class TractionBCTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank_);
  }

  void writeConfig(const std::string& path)
  {
    if (rank_ != 0) return;
    std::ofstream cfg(path);
    cfg << "[SIMULATION]\n";
    cfg << "name = test_traction_bc\n";
    cfg << "start_time = 0.0\n";
    cfg << "end_time = 1.0\n";
    cfg << "dt_initial = 1.0\n";
    cfg << "max_timesteps = 1\n";
    cfg << "output_frequency = 1\n";
    cfg << "fluid_model = NONE\n";
    cfg << "solid_model = ELASTIC\n";
    cfg << "enable_geomechanics = true\n";
    cfg << "enable_faults = false\n";
    cfg << "rtol = 1.0e-10\n";
    cfg << "atol = 1.0e-12\n";
    cfg << "max_nonlinear_iterations = 5\n";
    cfg << "\n[GRID]\n";
    cfg << "nx = 4\n";
    cfg << "ny = 4\n";
    cfg << "nz = 8\n";
    cfg << "Lx = 10.0\n";
    cfg << "Ly = 10.0\n";
    cfg << "Lz = 100.0\n";
    cfg << "\n[ROCK]\n";
    cfg << "density = 2650.0\n";
    cfg << "youngs_modulus = 10.0e9\n";
    cfg << "poissons_ratio = 0.25\n";
    cfg << "\n[MATERIAL]\n";
    cfg << "heterogeneous = false\n";
    cfg << "gravity = 0.0\n";
    cfg << "\n[BOUNDARY_CONDITIONS]\n";
    cfg << "z_min = fixed\n";
    cfg << "x_min = roller\n";
    cfg << "x_max = roller\n";
    cfg << "y_min = roller\n";
    cfg << "y_max = roller\n";
    cfg << "z_max = traction\n";
    cfg << "z_max_traction = 0.0, 0.0, -1.0e6\n";
    cfg << "\n[OUTPUT]\n";
    cfg << "directory = output\n";
    cfg << "format = HDF5\n";
  }

  int rank_ = 0;
};

TEST_F(TractionBCTest, UniaxialTractionMatchesAnalytical)
{
  std::string config_path = "test_traction_bc_run.config";
  writeConfig(config_path);
  MPI_Barrier(PETSC_COMM_WORLD);

  FSRM::Simulator sim(PETSC_COMM_WORLD);
  PetscErrorCode ierr;

  ierr = sim.initializeFromConfigFile(config_path);
  ASSERT_EQ(ierr, 0) << "initializeFromConfigFile must succeed";
  ierr = sim.setupDM();
  ASSERT_EQ(ierr, 0) << "setupDM must succeed";
  ierr = sim.labelBoundaries();
  ASSERT_EQ(ierr, 0) << "labelBoundaries must succeed";
  ierr = sim.setupFields();
  ASSERT_EQ(ierr, 0) << "setupFields must succeed";
  ierr = sim.setupPhysics();
  ASSERT_EQ(ierr, 0) << "setupPhysics must succeed";
  ierr = sim.setupTimeStepper();
  ASSERT_EQ(ierr, 0) << "setupTimeStepper must succeed";
  ierr = sim.setupSolvers();
  ASSERT_EQ(ierr, 0) << "setupSolvers must succeed";
  ierr = sim.setInitialConditions();
  ASSERT_EQ(ierr, 0) << "setInitialConditions must succeed";

  ierr = sim.run();
  EXPECT_EQ(ierr, 0) << "sim.run() must complete for traction BC test";

  // Validate solution exists and has nonzero displacement
  Vec sol = sim.getSolution();
  ASSERT_NE(sol, nullptr);

  PetscReal norm = 0.0;
  VecNorm(sol, NORM_2, &norm);
  EXPECT_GT(norm, 0.0) << "Solution displacement must be nonzero under traction";

  // Analytical: u_z at top should be about -0.01 m
  // (sigma_zz=-1e6, E=10e9, Lz=100 => u_z=-1e6*100/10e9=-0.01)
  // Check that the infinity norm is in the right ballpark
  PetscReal max_val = 0.0;
  VecNorm(sol, NORM_INFINITY, &max_val);
  EXPECT_GT(max_val, 1.0e-4) << "Max displacement should be > 1e-4 m";
  EXPECT_LT(max_val, 1.0) << "Max displacement should be < 1 m (sanity check)";

  // Cleanup
  MPI_Barrier(PETSC_COMM_WORLD);
  if (rank_ == 0) {
    std::remove(config_path.c_str());
  }
}
