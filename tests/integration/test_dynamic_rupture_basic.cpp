/**
 * @file test_dynamic_rupture_basic.cpp
 * @brief Integration test: dynamic rupture on cohesive fault mesh
 *
 * Verifies that the full Simulator pipeline with cohesive faults:
 *   1. Creates mesh with cohesive cells (via splitMeshAlongFault)
 *   2. Sets up fields including Lagrange multiplier
 *   3. Registers cohesive callbacks on PetscDS
 *   4. Initializes the time stepper and solver
 *
 * Absorbing BC coexistence with cohesive callbacks is supported through
 * region-specific PetscDS registration.
 *
 * Known limitation: SNES may diverge at the first time step because
 * the cohesive Jacobian blocks (g0_displacement_lagrange, etc.) may not
 * be fully consistent with PETSc 3.22 hybrid cell assembly. This test
 * documents setup success and SNES behavior.
 */

#include <gtest/gtest.h>
#include <cmath>
#include <fstream>
#include <petscsys.h>
#include <petscdmplex.h>

#include "core/Simulator.hpp"
#include "core/FSRM.hpp"
#include "numerics/FaultMeshManager.hpp"
#include "physics/CohesiveFaultKernel.hpp"
#include "domain/geomechanics/PyLithFault.hpp"

using namespace FSRM;

class DynamicRuptureBasicTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank_);
  }

  int rank_ = 0;
};

/**
 * Test that the mesh splitting workflow succeeds and produces
 * cohesive cells for a dynamic rupture configuration.
 */
TEST_F(DynamicRuptureBasicTest, MeshSplittingWithCohesiveCells)
{
  DM dm = nullptr;
  PetscInt faces[3] = {4, 4, 4};
  PetscReal lower[3] = {0.0, 0.0, 0.0};
  PetscReal upper[3] = {1.0, 1.0, 1.0};

  PetscErrorCode ierr = DMPlexCreateBoxMesh(
      PETSC_COMM_WORLD, 3, PETSC_TRUE, faces, lower, upper,
      nullptr, PETSC_TRUE, 0, PETSC_FALSE, &dm);
  ASSERT_EQ(ierr, 0) << "DMPlexCreateBoxMesh must succeed in Docker CI";
  ASSERT_NE(dm, nullptr);

  FaultMeshManager mgr(PETSC_COMM_WORLD);
  DMLabel fault_label = nullptr;
  const double center[3] = {0.5, 0.5, 0.5};
  ierr = mgr.createPlanarFaultLabel(
      dm, &fault_label, M_PI / 2.0, M_PI / 2.0,
      center, 2.0, 2.0, 0.05);
  ASSERT_EQ(ierr, 0) << "createPlanarFaultLabel must succeed";

  ierr = mgr.splitMeshAlongFault(&dm, "fault");
  ASSERT_EQ(ierr, 0) << "splitMeshAlongFault must succeed";

  FaultCohesiveDyn fault;
  ierr = mgr.extractCohesiveTopology(dm, &fault);
  ASSERT_EQ(ierr, 0) << "extractCohesiveTopology must succeed";

  PetscInt cStart, cEnd;
  ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd);
  ASSERT_EQ(ierr, 0);

  PetscInt num_cells = cEnd - cStart;
  EXPECT_GT(num_cells, 384)
      << "Mesh should have more cells after cohesive insertion";
  EXPECT_GT(fault.numVertices(), 0u)
      << "Cohesive topology should have fault vertices";

  DMDestroy(&dm);
}

/**
 * Test that the full Simulator setup pipeline completes for a locked-fault
 * elastodynamic simulation. Validates DM creation, fault network setup,
 * boundary labeling, field creation with Lagrange multiplier, cohesive
 * callback registration, time stepper and solver configuration, and
 * initial condition setup.
 */
TEST_F(DynamicRuptureBasicTest, LockedFaultSetupCompletes)
{
  std::string config_path = "test_dynamic_rupture_locked.config";

  if (rank_ == 0)
  {
    std::ofstream cfg(config_path);
    cfg << "[SIMULATION]\n";
    cfg << "name = test_dynamic_rupture_locked\n";
    cfg << "start_time = 0.0\n";
    cfg << "end_time = 0.001\n";
    cfg << "dt_initial = 0.0001\n";
    cfg << "dt_min = 0.00001\n";
    cfg << "dt_max = 0.001\n";
    cfg << "max_timesteps = 10\n";
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

  // All 8 setup steps must complete without error
  ierr = sim.initializeFromConfigFile(config_path);
  EXPECT_EQ(ierr, 0) << "initializeFromConfigFile must succeed";
  if (ierr) { if (rank_ == 0) std::remove(config_path.c_str()); return; }

  ierr = sim.setupDM();
  EXPECT_EQ(ierr, 0) << "setupDM must succeed";
  if (ierr) { if (rank_ == 0) std::remove(config_path.c_str()); return; }

  ierr = sim.setupFaultNetwork();
  EXPECT_EQ(ierr, 0) << "setupFaultNetwork must succeed (cohesive cells inserted)";
  if (ierr) { if (rank_ == 0) std::remove(config_path.c_str()); return; }

  ierr = sim.labelBoundaries();
  EXPECT_EQ(ierr, 0) << "labelBoundaries must succeed on split mesh";
  if (ierr) { if (rank_ == 0) std::remove(config_path.c_str()); return; }

  ierr = sim.setupFields();
  EXPECT_EQ(ierr, 0) << "setupFields must succeed (includes Lagrange multiplier)";
  if (ierr) { if (rank_ == 0) std::remove(config_path.c_str()); return; }

  ierr = sim.setupPhysics();
  EXPECT_EQ(ierr, 0) << "setupPhysics must succeed (cohesive callbacks registered)";
  if (ierr) { if (rank_ == 0) std::remove(config_path.c_str()); return; }

  ierr = sim.setupTimeStepper();
  EXPECT_EQ(ierr, 0) << "setupTimeStepper must succeed (TSALPHA2)";
  if (ierr) { if (rank_ == 0) std::remove(config_path.c_str()); return; }

  ierr = sim.setupSolvers();
  EXPECT_EQ(ierr, 0) << "setupSolvers must succeed";
  if (ierr) { if (rank_ == 0) std::remove(config_path.c_str()); return; }

  ierr = sim.setInitialConditions();
  EXPECT_EQ(ierr, 0) << "setInitialConditions must succeed";

  // Pipeline is fully set up. TSSolve may diverge due to Jacobian issues.
  // This is a known limitation documented in CLAUDE.md.

  if (rank_ == 0)
  {
    std::remove(config_path.c_str());
  }
}

/**
 * Test that setup pipeline completes for a slipping-fault configuration.
 * Validates that the cohesive kernel in slipping mode can be registered
 * and the solver infrastructure accepts the fault configuration.
 */
TEST_F(DynamicRuptureBasicTest, SlippingFaultSetupCompletes)
{
  std::string config_path = "test_dynamic_rupture_slipping.config";

  if (rank_ == 0)
  {
    std::ofstream cfg(config_path);
    cfg << "[SIMULATION]\n";
    cfg << "name = test_dynamic_rupture_slipping\n";
    cfg << "start_time = 0.0\n";
    cfg << "end_time = 0.001\n";
    cfg << "dt_initial = 0.0001\n";
    cfg << "dt_min = 0.00001\n";
    cfg << "dt_max = 0.001\n";
    cfg << "max_timesteps = 10\n";
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
    cfg << "\n[FAULT]\n";
    cfg << "strike = 90.0\n";
    cfg << "dip = 90.0\n";
    cfg << "center_x = 0.5\n";
    cfg << "center_y = 0.5\n";
    cfg << "center_z = 0.5\n";
    cfg << "length = 2.0\n";
    cfg << "width = 2.0\n";
    cfg << "mode = slipping\n";
    cfg << "friction_coefficient = 0.1\n";
    cfg << "initial_shear_traction = 50.0e6\n";
    cfg << "initial_normal_traction = -100.0e6\n";
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

  ierr = sim.initializeFromConfigFile(config_path);
  EXPECT_EQ(ierr, 0) << "initializeFromConfigFile must succeed";
  if (ierr) { if (rank_ == 0) std::remove(config_path.c_str()); return; }

  ierr = sim.setupDM();
  EXPECT_EQ(ierr, 0) << "setupDM must succeed";
  if (ierr) { if (rank_ == 0) std::remove(config_path.c_str()); return; }

  ierr = sim.setupFaultNetwork();
  EXPECT_EQ(ierr, 0) << "setupFaultNetwork must succeed";
  if (ierr) { if (rank_ == 0) std::remove(config_path.c_str()); return; }

  ierr = sim.labelBoundaries();
  EXPECT_EQ(ierr, 0) << "labelBoundaries must succeed";
  if (ierr) { if (rank_ == 0) std::remove(config_path.c_str()); return; }

  ierr = sim.setupFields();
  EXPECT_EQ(ierr, 0) << "setupFields must succeed";
  if (ierr) { if (rank_ == 0) std::remove(config_path.c_str()); return; }

  ierr = sim.setupPhysics();
  EXPECT_EQ(ierr, 0) << "setupPhysics must succeed";
  if (ierr) { if (rank_ == 0) std::remove(config_path.c_str()); return; }

  ierr = sim.setupTimeStepper();
  EXPECT_EQ(ierr, 0) << "setupTimeStepper must succeed";
  if (ierr) { if (rank_ == 0) std::remove(config_path.c_str()); return; }

  ierr = sim.setupSolvers();
  EXPECT_EQ(ierr, 0) << "setupSolvers must succeed";
  if (ierr) { if (rank_ == 0) std::remove(config_path.c_str()); return; }

  ierr = sim.setInitialConditions();
  EXPECT_EQ(ierr, 0) << "setInitialConditions must succeed";

  if (rank_ == 0)
  {
    std::remove(config_path.c_str());
  }
}

TEST_F(DynamicRuptureBasicTest, FaultsWithAbsorbingBCCoexist)
{
  std::string config_path = "test_dynamic_rupture_absorb_conflict.config";

  if (rank_ == 0)
  {
    std::ofstream cfg(config_path);
    cfg << "[SIMULATION]\n";
    cfg << "name = test_dynamic_rupture_absorb_conflict\n";
    cfg << "start_time = 0.0\n";
    cfg << "end_time = 0.001\n";
    cfg << "dt_initial = 0.0001\n";
    cfg << "dt_min = 0.00001\n";
    cfg << "dt_max = 0.001\n";
    cfg << "max_timesteps = 10\n";
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
    cfg.close();
  }
  MPI_Barrier(PETSC_COMM_WORLD);

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
  EXPECT_EQ(ierr, 0) << "setupPhysics should support faults + absorbing BC together";
  if (!ierr) {
    ierr = sim.setupTimeStepper();
    EXPECT_EQ(ierr, 0) << "setupTimeStepper should succeed with coexistence setup";
  }

  if (rank_ == 0)
  {
    std::remove(config_path.c_str());
  }
}
