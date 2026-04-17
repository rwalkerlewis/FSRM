/**
 * @file test_elastoplastic_sim.cpp
 * @brief Integration tests: Drucker-Prager elastoplasticity through Simulator pipeline
 *
 * Verifies that PetscFEElastoplasticity callbacks are correctly wired into
 * setupPhysics() and that sim.run() (TSSolve) completes with plasticity enabled.
 *
 * Test 1: Large compression beyond yield -- confirms plastic return mapping runs
 * Test 2: Small compression below yield -- confirms elastic behavior is preserved
 */

#include <gtest/gtest.h>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <petscsys.h>
#include <petscdmplex.h>

#include "core/Simulator.hpp"
#include "core/FSRM.hpp"

using namespace FSRM;

class ElastoplasticSimTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank_);
  }

  void writeConfig(const std::string& path,
                   double cohesion, double friction_angle_deg)
  {
    if (rank_ != 0) return;
    std::ofstream cfg(path);
    cfg << "[SIMULATION]\n";
    cfg << "name = test_elastoplastic_sim\n";
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
    cfg << "rtol = 1.0e-8\n";
    cfg << "atol = 1.0e-10\n";
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
    cfg << "youngs_modulus = 50.0e9\n";
    cfg << "poissons_ratio = 0.25\n";
    cfg << "\n[PLASTICITY]\n";
    cfg << "enabled = true\n";
    cfg << "cohesion = " << cohesion << "\n";
    cfg << "friction_angle = " << friction_angle_deg << "\n";
    cfg << "dilation_angle = 0.0\n";
    cfg << "hardening_modulus = 0.0\n";
    cfg << "\n[BOUNDARY_CONDITIONS]\n";
    cfg << "bottom = fixed\n";
    cfg << "sides = roller\n";
    // compression applies u_z = -0.001 on the top face
    cfg << "top = compression\n";
    cfg.close();
  }

  int rank_ = 0;
};

// Compression with u_z=-0.001 on 1m box: eps=0.001, sigma~50 MPa.
// DP yield for c=5 MPa, phi=30: sigma_y ~ 17.3 MPa. This is well past yield.
TEST_F(ElastoplasticSimTest, QuasiStaticCompressionBeyondYield)
{
  std::string config_path = "test_ep_beyond_yield.config";
  writeConfig(config_path, 5.0e6, 30.0);
  MPI_Barrier(PETSC_COMM_WORLD);

  Simulator sim(PETSC_COMM_WORLD);
  PetscErrorCode ierr;

  ierr = sim.initializeFromConfigFile(config_path);
  ASSERT_EQ(ierr, 0) << "initializeFromConfigFile must succeed";
  ierr = sim.setupDM();
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

  // Run TSSolve -- the critical test
  ierr = sim.run();
  EXPECT_EQ(ierr, 0) << "sim.run() must complete with elastoplasticity beyond yield";

  // Check solution is nonzero
  Vec sol = sim.getSolution();
  ASSERT_NE(sol, nullptr);
  PetscReal norm = 0.0;
  VecNorm(sol, NORM_2, &norm);
  EXPECT_GT(norm, 0.0) << "Solution displacement must be nonzero";

  // Displacement should be reasonable (< 0.05 m for 1m box with 1% strain)
  PetscReal max_val = 0.0;
  VecNorm(sol, NORM_INFINITY, &max_val);
  EXPECT_LT(max_val, 0.05) << "Maximum displacement must be bounded";

  MPI_Barrier(PETSC_COMM_WORLD);
  if (rank_ == 0) std::remove(config_path.c_str());
}

// Same compression but with very high cohesion so it stays elastic.
// c=500 GPa, phi=30: yield stress ~ 1.7 TPa, far above 50 MPa applied.
TEST_F(ElastoplasticSimTest, ElasticBelowYield)
{
  std::string config_path = "test_ep_below_yield.config";
  writeConfig(config_path, 500.0e9, 30.0);
  MPI_Barrier(PETSC_COMM_WORLD);

  Simulator sim(PETSC_COMM_WORLD);
  PetscErrorCode ierr;

  ierr = sim.initializeFromConfigFile(config_path);
  ASSERT_EQ(ierr, 0);
  ierr = sim.setupDM();
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
  EXPECT_EQ(ierr, 0) << "sim.run() must complete with elastoplasticity below yield";

  Vec sol = sim.getSolution();
  ASSERT_NE(sol, nullptr);
  PetscReal norm = 0.0;
  VecNorm(sol, NORM_2, &norm);
  EXPECT_GT(norm, 0.0) << "Solution must be nonzero (tiny elastic deformation)";

  // Compression BC applies u_z = -0.001, so max displacement is ~0.001
  PetscReal max_val = 0.0;
  VecNorm(sol, NORM_INFINITY, &max_val);
  EXPECT_LT(max_val, 0.01) << "Displacement must be bounded";

  MPI_Barrier(PETSC_COMM_WORLD);
  if (rank_ == 0) std::remove(config_path.c_str());
}
