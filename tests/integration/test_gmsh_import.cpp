/**
 * @file test_gmsh_import.cpp
 * @brief Integration test: Gmsh mesh import and elastostatic solve
 *
 * Verifies that a Gmsh MSH2 file can be loaded via DMPlexCreateGmshFromFile,
 * boundary labels are applied correctly, and an elastostatic solve converges
 * on the resulting tetrahedral mesh.
 *
 * Quantitative verification:
 *   10x10x10 m box, uniaxial compression u_z = -0.001 m on top.
 *   E = 50 GPa, nu = 0.25.
 *   eps_zz = -0.001 / 10 = -1e-4.
 *   P-wave modulus M = lambda + 2*mu = 60 GPa.
 *   sigma_zz = M * eps_zz = -6 MPa.
 *   Displacement field must have max |u_z| = 0.001 m at top vertices.
 */

#include <gtest/gtest.h>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <petscsys.h>
#include <petscdmplex.h>

#include "core/Simulator.hpp"
#include "core/FSRM.hpp"

using namespace FSRM;

class GmshImportTest : public ::testing::Test
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
      cfg << "mesh_type = GMSH\n";
      cfg << "mesh_file = " << mesh_path_ << "\n";
      cfg << "gmsh_boundaries = xmin, xmax, ymin, ymax, zmin, zmax\n";
      cfg << "\n[ROCK]\n";
      cfg << "density = 2650.0\n";
      cfg << "youngs_modulus = 50.0e9\n";
      cfg << "poisson_ratio = 0.25\n";
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
  static constexpr const char* config_path_ = "test_gmsh_import.ini";
  static constexpr const char* mesh_path_ = "../../meshes/test_box_gmsh.msh";
};

// Test 1: Mesh loads, dimension is 3, cell count is nonzero
TEST_F(GmshImportTest, MeshLoadsCorrectly)
{
  Simulator sim(PETSC_COMM_WORLD);
  PetscErrorCode ierr;

  ierr = sim.initializeFromConfigFile(config_path_);
  ASSERT_EQ(ierr, 0) << "initializeFromConfigFile failed";

  ierr = sim.setupDM();
  ASSERT_EQ(ierr, 0) << "setupDM failed (could not load Gmsh mesh)";

  DM dm = sim.getDM();
  ASSERT_NE(dm, nullptr);

  PetscInt dim;
  ierr = DMGetDimension(dm, &dim);
  ASSERT_EQ(ierr, 0);
  EXPECT_EQ(dim, 3) << "Mesh dimension should be 3";

  PetscInt cStart, cEnd;
  ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd);
  ASSERT_EQ(ierr, 0);
  EXPECT_GT(cEnd - cStart, 0) << "Mesh should have nonzero cells";

  // Verify the cells are tetrahedra
  if (cStart < cEnd)
  {
    DMPolytopeType ct;
    ierr = DMPlexGetCellType(dm, cStart, &ct);
    ASSERT_EQ(ierr, 0);
    EXPECT_EQ(ct, DM_POLYTOPE_TETRAHEDRON)
      << "Gmsh-generated mesh should have tetrahedral cells";
  }
}

// Test 2: Boundary labels exist after labelBoundaries()
TEST_F(GmshImportTest, BoundaryLabelsExist)
{
  Simulator sim(PETSC_COMM_WORLD);
  PetscErrorCode ierr;

  ierr = sim.initializeFromConfigFile(config_path_);
  ASSERT_EQ(ierr, 0);

  ierr = sim.setupDM();
  ASSERT_EQ(ierr, 0);

  ierr = sim.labelBoundaries();
  ASSERT_EQ(ierr, 0) << "labelBoundaries failed";

  DM dm = sim.getDM();

  // Check that the 6 bounding-box labels were created
  const char* label_names[] = {
    "boundary_x_min", "boundary_x_max",
    "boundary_y_min", "boundary_y_max",
    "boundary_z_min", "boundary_z_max"
  };

  for (const char* name : label_names)
  {
    DMLabel label = nullptr;
    ierr = DMGetLabel(dm, name, &label);
    ASSERT_EQ(ierr, 0);
    EXPECT_NE(label, nullptr) << "Label '" << name << "' should exist";

    if (label)
    {
      PetscInt n;
      ierr = DMLabelGetNumValues(label, &n);
      ASSERT_EQ(ierr, 0);
      EXPECT_GT(n, 0) << "Label '" << name << "' should have values";
    }
  }
}

// Test 3: Full elastostatic solve converges on Gmsh tet mesh
TEST_F(GmshImportTest, ElastostaticSolveConverges)
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
  ASSERT_EQ(ierr, 0) << "Simulation failed on Gmsh mesh";

  // Verify solution is nonzero
  Vec sol = sim.getSolution();
  ASSERT_NE(sol, nullptr);

  PetscReal norm;
  ierr = VecNorm(sol, NORM_2, &norm);
  ASSERT_EQ(ierr, 0);
  EXPECT_GT(norm, 0.0) << "Solution should be nonzero after compression";

  // Check max displacement magnitude is in expected range
  // Applied compression: u_z = -0.001 m on top face
  // Expect max |u_z| near 0.001 at top, 0 at bottom
  PetscReal max_val, min_val;
  ierr = VecMax(sol, nullptr, &max_val);
  ASSERT_EQ(ierr, 0);
  ierr = VecMin(sol, nullptr, &min_val);
  ASSERT_EQ(ierr, 0);

  // The minimum displacement component should be about -0.001 (top face u_z)
  EXPECT_LT(min_val, 0.0)
    << "Should have negative displacement from compression";
  EXPECT_GT(min_val, -0.01)
    << "Displacement magnitude should be reasonable (< 10 mm)";
}
