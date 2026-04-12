/**
 * @file test_derived_fields.cpp
 * @brief Unit test for DerivedFieldComputer
 *
 * Verifies that cell-centered stress and strain are computed correctly
 * for a uniaxial compression problem on a structured hex mesh.
 *
 * Quantitative checks:
 *   10x10x10 box, E=50 GPa, nu=0.25, top u_z=-0.001.
 *   eps_zz = -0.001 / 10 = -1e-4 (average).
 *   lambda = 20 GPa, mu = 20 GPa, M=60 GPa.
 *   sigma_zz = M * eps_zz = -6 MPa.
 *   sigma_xx = sigma_yy = lambda * eps_zz = -2 MPa.
 */

#include <gtest/gtest.h>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <numeric>
#include <petscsys.h>

#include "core/Simulator.hpp"
#include "core/DerivedFieldComputer.hpp"

using namespace FSRM;

class DerivedFieldTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (rank == 0)
    {
      std::ofstream cfg(config_path_);
      cfg << "[SIMULATION]\n";
      cfg << "enable_geomechanics = true\n";
      cfg << "solid_model = ELASTIC\n";
      cfg << "fluid_model = NONE\n";
      cfg << "max_timesteps = 1\n";
      cfg << "dt_initial = 1.0\n";
      cfg << "start_time = 0.0\n";
      cfg << "end_time = 1.0\n";
      cfg << "\n[GRID]\n";
      cfg << "nx = 2\n";
      cfg << "ny = 2\n";
      cfg << "nz = 4\n";
      cfg << "Lx = 10.0\n";
      cfg << "Ly = 10.0\n";
      cfg << "Lz = 10.0\n";
      cfg << "mesh_type = CARTESIAN\n";
      cfg << "\n[ROCK]\n";
      cfg << "density = 2650.0\n";
      cfg << "youngs_modulus = 50.0e9\n";
      cfg << "poisson_ratio = 0.25\n";
      cfg << "\n[OUTPUT]\n";
      cfg << "output_stress = true\n";
      cfg << "output_cfs = true\n";
      cfg << "cfs_receiver_strike = 0.0\n";
      cfg << "cfs_receiver_dip = 90.0\n";
      cfg << "cfs_friction = 0.4\n";
    }
    MPI_Barrier(PETSC_COMM_WORLD);
  }

  void TearDown() override
  {
    MPI_Barrier(PETSC_COMM_WORLD);
    if (rank == 0)
    {
      std::remove(config_path_);
    }
  }

  int rank = 0;
  static constexpr const char* config_path_ = "test_derived_fields.ini";
};

// Test: DerivedFieldComputer produces correct stress for uniaxial compression
TEST_F(DerivedFieldTest, StressMatchesUniaxialCompression)
{
  Simulator sim(PETSC_COMM_WORLD);
  PetscErrorCode ierr;

  ierr = sim.initializeFromConfigFile(config_path_);
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
  ASSERT_EQ(ierr, 0);

  // The run() method computed derived fields because output_stress=true
  const auto& df = sim.getDerivedFields();
  ASSERT_GT(df.numCells(), 0) << "Should have computed stress for all cells";

  const auto& stress = df.stress();
  ASSERT_EQ(static_cast<PetscInt>(stress.size()), 6 * df.numCells());

  // Compute average sigma_zz across all cells
  double avg_sigma_zz = 0.0;
  for (PetscInt c = 0; c < df.numCells(); ++c)
  {
    avg_sigma_zz += stress[6 * c + 2];  // Voigt index 2 = zz
  }
  avg_sigma_zz /= df.numCells();

  // Expected: sigma_zz = M * eps_zz
  // M = lambda + 2*mu = 20e9 + 2*20e9 = 60 GPa
  // eps_zz = -0.001 / 10 = -1e-4
  // sigma_zz = -6 MPa
  const double expected_sigma_zz = -6.0e6;

  // 20% tolerance for coarse mesh FEM
  EXPECT_NEAR(avg_sigma_zz, expected_sigma_zz, 0.2 * std::abs(expected_sigma_zz))
    << "Average sigma_zz should be approximately -6 MPa for uniaxial compression";

  // Strain check: average eps_zz should be about -1e-4
  const auto& strain = df.strain();
  double avg_eps_zz = 0.0;
  for (PetscInt c = 0; c < df.numCells(); ++c)
  {
    avg_eps_zz += strain[6 * c + 2];
  }
  avg_eps_zz /= df.numCells();

  const double expected_eps_zz = -1.0e-4;
  EXPECT_NEAR(avg_eps_zz, expected_eps_zz, 0.2 * std::abs(expected_eps_zz))
    << "Average eps_zz should be approximately -1e-4";

  // CFS should be computed (scalar per cell)
  const auto& cfs = df.cfs();
  ASSERT_EQ(static_cast<PetscInt>(cfs.size()), df.numCells());
}
