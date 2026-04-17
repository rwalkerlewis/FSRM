/**
 * @file test_output_file.cpp
 * @brief Integration test: verify that writeOutput produces an HDF5 file
 *
 * Runs a minimal uniaxial compression simulation (1 timestep) and checks
 * that an HDF5 output file is created and has nonzero size.
 */

#include "core/Simulator.hpp"
#include "core/FSRM.hpp"
#include <gtest/gtest.h>
#include <fstream>
#include <cstdio>
#include <sys/stat.h>

using namespace FSRM;

class OutputFileTest : public ::testing::Test {
protected:
  void SetUp() override
  {
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // Create output directory
    if (rank == 0) {
      mkdir(output_dir_, 0755);
    }
    MPI_Barrier(PETSC_COMM_WORLD);

    // Write config file
    if (rank == 0) {
      std::ofstream cfg(config_path_);
      cfg << "[SIMULATION]\n"
          << "name = output_test\n"
          << "start_time = 0.0\n"
          << "end_time = 1.0\n"
          << "dt_initial = 1.0\n"
          << "max_timesteps = 1\n"
          << "output_frequency = 1\n"
          << "fluid_model = NONE\n"
          << "solid_model = ELASTIC\n"
          << "enable_geomechanics = true\n"
          << "enable_faults = false\n"
          << "output_format = HDF5\n"
          << "rtol = 1.0e-10\n"
          << "atol = 1.0e-12\n"
          << "max_nonlinear_iterations = 5\n"
          << "\n"
          << "[GRID]\n"
          << "nx = 2\n"
          << "ny = 2\n"
          << "nz = 4\n"
          << "Lx = 10.0\n"
          << "Ly = 10.0\n"
          << "Lz = 100.0\n"
          << "\n"
          << "[ROCK]\n"
          << "density = 2650.0\n"
          << "youngs_modulus = 10.0e9\n"
          << "poissons_ratio = 0.25\n"
          << "\n"
          << "[OUTPUT]\n"
          << "output_directory = " << output_dir_ << "\n";
      cfg.close();
    }
    MPI_Barrier(PETSC_COMM_WORLD);
  }

  void TearDown() override
  {
    MPI_Barrier(PETSC_COMM_WORLD);
    if (rank == 0) {
      std::remove(h5_path_);
      std::remove(config_path_);
      rmdir(output_dir_);
    }
  }

  int rank = 0;
  static constexpr const char* config_path_ = "test_output_file.config";
  static constexpr const char* output_dir_ = "test_output_dir";
  static constexpr const char* h5_path_ = "test_output_dir/solution.h5";
};

TEST_F(OutputFileTest, HDF5OutputFileCreated)
{
  Simulator sim(PETSC_COMM_WORLD);
  PetscErrorCode ierr;

  ierr = sim.initializeFromConfigFile(config_path_);
  ASSERT_EQ(ierr, 0) << "initializeFromConfigFile failed";

  ierr = sim.setupDM();
  ASSERT_EQ(ierr, 0) << "setupDM failed";

  ierr = sim.labelBoundaries();
  ASSERT_EQ(ierr, 0) << "labelBoundaries failed";

  ierr = sim.setupFields();
  ASSERT_EQ(ierr, 0) << "setupFields failed";

  ierr = sim.setupPhysics();
  ASSERT_EQ(ierr, 0) << "setupPhysics failed";

  ierr = sim.setupTimeStepper();
  ASSERT_EQ(ierr, 0) << "setupTimeStepper failed";

  ierr = sim.setupSolvers();
  ASSERT_EQ(ierr, 0) << "setupSolvers failed";

  ierr = sim.setInitialConditions();
  ASSERT_EQ(ierr, 0) << "setInitialConditions failed";

  ierr = sim.run();
  ASSERT_EQ(ierr, 0) << "run failed";

  // Verify HDF5 output file exists and has nonzero size
  if (rank == 0) {
    struct stat st;
    int stat_rc = stat(h5_path_, &st);
    ASSERT_EQ(stat_rc, 0) << "HDF5 output file does not exist: " << h5_path_;
    EXPECT_GT(st.st_size, 0) << "HDF5 output file is empty";
  }
}
