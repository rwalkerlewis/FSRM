/**
 * @file test_fault_absorbing_coexist.cpp
 * @brief Integration test for cohesive + absorbing boundary callback coexistence
 */

#include <gtest/gtest.h>
#include <cstdio>
#include <fstream>
#include <petscsys.h>

#include "core/Simulator.hpp"

using namespace FSRM;

class FaultAbsorbingCoexistTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank_);
  }

  int rank_ = 0;
};

TEST_F(FaultAbsorbingCoexistTest, SetupPhysicsSupportsBothCallbackFamilies)
{
  const std::string config_path = "test_fault_absorbing_coexist.config";

  if (rank_ == 0)
  {
    std::ofstream cfg(config_path);
    cfg << "[SIMULATION]\n";
    cfg << "name = test_fault_absorbing_coexist\n";
    cfg << "start_time = 0.0\n";
    cfg << "end_time = 0.001\n";
    cfg << "dt_initial = 0.0001\n";
    cfg << "dt_min = 0.00001\n";
    cfg << "dt_max = 0.001\n";
    cfg << "max_timesteps = 4\n";
    cfg << "output_frequency = 100\n";
    cfg << "fluid_model = NONE\n";
    cfg << "solid_model = ELASTIC\n";
    cfg << "enable_geomechanics = true\n";
    cfg << "enable_elastodynamics = true\n";
    cfg << "enable_faults = true\n";
    cfg << "rtol = 1.0e-4\n";
    cfg << "atol = 1.0e-6\n";
    cfg << "max_nonlinear_iterations = 40\n";

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

  Simulator sim(PETSC_COMM_WORLD);
  PetscErrorCode ierr = sim.initializeFromConfigFile(config_path);
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
  EXPECT_EQ(ierr, 0) << "cohesive and absorbing callbacks should coexist via region DS";

  if (!ierr)
  {
    ierr = sim.setupTimeStepper();
    EXPECT_EQ(ierr, 0);
  }

  if (rank_ == 0)
  {
    std::remove(config_path.c_str());
  }
}
