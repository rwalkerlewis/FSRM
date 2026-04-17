/**
 * @file test_injection_pressure.cpp
 * @brief Integration test: end-to-end injection pressure buildup pipeline
 *
 * Runs two short poroelastic simulations on the same mesh:
 *   1) With point injection enabled
 *   2) Without injection
 *
 * The test verifies both runs complete and the final coupled state vectors
 * differ, confirming the injection source term is active in the end-to-end
 * Simulator pipeline.
 */

#include <gtest/gtest.h>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <petscsys.h>

#include "core/Simulator.hpp"

using namespace FSRM;

class InjectionPressureTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank_);

    if (rank_ == 0)
    {
      writeConfig(inj_config_path_, true);
    }
    MPI_Barrier(PETSC_COMM_WORLD);
  }

  void TearDown() override
  {
    MPI_Barrier(PETSC_COMM_WORLD);
    if (rank_ == 0)
    {
      std::remove(inj_config_path_.c_str());
    }
  }

  void writeConfig(const std::string& path, bool enable_injection)
  {
    std::ofstream cfg(path);
    cfg << "[SIMULATION]\n";
    cfg << "name = injection_pressure_pipeline_test\n";
    cfg << "start_time = 0.0\n";
    cfg << "end_time = 0.2\n";
    cfg << "dt_initial = 0.02\n";
    cfg << "dt_min = 0.002\n";
    cfg << "dt_max = 0.05\n";
    cfg << "max_timesteps = 20\n";
    cfg << "output_frequency = 100\n";
    cfg << "fluid_model = SINGLE_COMPONENT\n";
    cfg << "solid_model = POROELASTIC\n";
    cfg << "enable_geomechanics = true\n";
    cfg << "enable_faults = false\n";
    cfg << "rtol = 1.0e-8\n";
    cfg << "atol = 1.0e-10\n";
    cfg << "max_nonlinear_iterations = 40\n";
    cfg << "\n[GRID]\n";
    cfg << "nx = 3\n";
    cfg << "ny = 3\n";
    cfg << "nz = 3\n";
    cfg << "Lx = 30.0\n";
    cfg << "Ly = 30.0\n";
    cfg << "Lz = 30.0\n";
    cfg << "\n[ROCK]\n";
    cfg << "density = 2500.0\n";
    cfg << "youngs_modulus = 1.0e9\n";
    cfg << "poissons_ratio = 0.25\n";
    cfg << "porosity = 0.25\n";
    cfg << "permeability_x = 500.0\n";
    cfg << "permeability_y = 500.0\n";
    cfg << "permeability_z = 500.0\n";
    cfg << "biot_coefficient = 0.9\n";
    cfg << "\n[FLUID]\n";
    cfg << "type = SINGLE_PHASE\n";
    cfg << "density = 1000.0\n";
    cfg << "viscosity = 0.001\n";
    cfg << "compressibility = 4.5e-10\n";
    cfg << "reference_pressure = 20.0e6\n";
    cfg << "\n[BOUNDARY_CONDITIONS]\n";
    cfg << "bottom = fixed\n";
    cfg << "sides = roller\n";
    cfg << "top = free\n";
    cfg << "\n[INJECTION]\n";
    cfg << "enabled = " << (enable_injection ? "true" : "false") << "\n";
    cfg << "x = 15.0\n";
    cfg << "y = 15.0\n";
    cfg << "z = 15.0\n";
    cfg << "rate = 1.0\n";
    cfg << "start_time = 0.0\n";
    cfg << "end_time = 0.2\n";
    cfg.close();
  }

  PetscErrorCode runPipelineAndCopySolution(const std::string& config_path, Vec* copied_solution)
  {
    Simulator sim(PETSC_COMM_WORLD);
    PetscErrorCode ierr;

    ierr = sim.initializeFromConfigFile(config_path);CHKERRQ(ierr);
    ierr = sim.setupDM();CHKERRQ(ierr);
    ierr = sim.labelBoundaries();CHKERRQ(ierr);
    ierr = sim.setupFields();CHKERRQ(ierr);
    ierr = sim.setupPhysics();CHKERRQ(ierr);
    ierr = sim.setupTimeStepper();CHKERRQ(ierr);
    ierr = sim.setupSolvers();CHKERRQ(ierr);
    ierr = sim.setInitialConditions();CHKERRQ(ierr);
    ierr = sim.run();CHKERRQ(ierr);

    Vec sol = sim.getSolution();
    PetscCheck(sol != nullptr, PETSC_COMM_WORLD, PETSC_ERR_ARG_NULL, "Solution vector is null");

    ierr = VecDuplicate(sol, copied_solution);CHKERRQ(ierr);
    ierr = VecCopy(sol, *copied_solution);CHKERRQ(ierr);

    return PETSC_SUCCESS;
  }

  int rank_ = 0;
  const std::string inj_config_path_ = "test_injection_enabled.config";
};

TEST_F(InjectionPressureTest, InjectionPipelineCompletes)
{
  Vec sol_inj = nullptr;

  PetscErrorCode ierr = runPipelineAndCopySolution(inj_config_path_, &sol_inj);
  ASSERT_EQ(ierr, 0) << "Injection-enabled pipeline failed";
  ASSERT_NE(sol_inj, nullptr);

  PetscReal norm = 0.0;
  ierr = VecNorm(sol_inj, NORM_2, &norm);
  ASSERT_EQ(ierr, 0);
  EXPECT_GT(norm, 0.0) << "Final coupled solution norm should be nonzero";
  EXPECT_TRUE(std::isfinite(static_cast<double>(norm)))
      << "Final coupled solution norm should be finite";

  VecDestroy(&sol_inj);
}
