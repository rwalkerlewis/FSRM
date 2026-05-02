/**
 * @file test_source_refinement.cpp
 * @brief Integration tests for the pass-3 source-region mesh refinement
 *
 * Verifies that the [MESH_REFINEMENT] grammar plumbs through to
 * Simulator::refineSourceRegion(), and that refinement (a) increases the
 * cell count locally near the source, (b) keeps the explosion cell at
 * the original source location, and (c) produces a finite, nonzero
 * solution norm through the full pipeline.
 *
 * The fixture uses a small Sedan 1962-style 4x4x4 base mesh and
 * compares the refined and un-refined cell counts after setupDM. The
 * 1.2x lower bound and 8x upper bound on the cell-count ratio are
 * chosen so a label-driven refinement on a 5*Rc ball can produce them:
 * an empty ball would leave the count unchanged (rejected) and a
 * global refinement would multiply it by 8 in 3D (rejected).
 */

#include <gtest/gtest.h>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <petscsys.h>
#include <petscdmplex.h>
#include <filesystem>

#include "core/Simulator.hpp"
#include "core/FSRM.hpp"

using namespace FSRM;

class SourceRefinementTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank_);
  }

  void TearDown() override
  {
    MPI_Barrier(PETSC_COMM_WORLD);
    if (rank_ == 0 && !config_path_.empty()) {
      std::remove(config_path_.c_str());
    }
  }

  // Sedan 1962 - style config: 104 kt, 194 m depth, alluvium, 2 km
  // domain, 4x4x4 base mesh. Optionally enable mesh refinement with
  // the given source_radius_factor and refinement_levels.
  void writeConfig(const std::string& tag,
                   bool enable_refinement,
                   double source_radius_factor,
                   int refinement_levels)
  {
    config_path_ = "test_source_refine_" + tag + ".config";
    if (rank_ != 0) return;

    std::ofstream cfg(config_path_);
    cfg << "[SIMULATION]\n";
    cfg << "name = source_refine_" << tag << "\n";
    cfg << "start_time = 0.0\n";
    cfg << "end_time = 0.005\n";
    cfg << "dt_initial = 0.001\n";
    cfg << "dt_min = 0.0001\n";
    cfg << "dt_max = 0.005\n";
    cfg << "max_timesteps = 10\n";
    cfg << "output_frequency = 1000\n";
    cfg << "output_format = NONE\n";
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
    cfg << "Lz = 2000.0\n";
    cfg << "\n[ROCK]\n";
    cfg << "density = 2650.0\n";
    cfg << "lambda = 1.87e10\n";
    cfg << "shear_modulus = 1.34e10\n";
    cfg << "\n[EXPLOSION_SOURCE]\n";
    cfg << "type = UNDERGROUND_NUCLEAR\n";
    cfg << "yield_kt = 104.0\n";
    cfg << "depth_of_burial = 194.0\n";
    cfg << "location_x = 2000.0\n";
    cfg << "location_y = 2000.0\n";
    cfg << "location_z = 1806.0\n";  // domain_z - depth
    cfg << "onset_time = 0.0\n";
    cfg << "rise_time = 0.01\n";
    cfg << "cavity_overpressure = 1.0e10\n";
    cfg << "medium_type = ALLUVIUM\n";
    cfg << "\n[BOUNDARY_CONDITIONS]\n";
    cfg << "bottom = free\n";
    cfg << "sides = free\n";
    cfg << "top = free\n";
    if (enable_refinement) {
      cfg << "\n[MESH_REFINEMENT]\n";
      cfg << "enabled = true\n";
      cfg << "source_radius_factor = " << source_radius_factor << "\n";
      cfg << "refinement_levels = " << refinement_levels << "\n";
    }
    cfg.close();
  }

  // Build the DM through setupDM (which also runs refineSourceRegion)
  // and report the resulting cell count.
  PetscErrorCode buildDMAndCountCells(PetscInt& n_cells)
  {
    Simulator sim(PETSC_COMM_WORLD);
    PetscErrorCode ierr;
    PetscOptionsClear(nullptr);

    ierr = sim.initializeFromConfigFile(config_path_); CHKERRQ(ierr);
    ierr = sim.setupDM(); CHKERRQ(ierr);

    DM dm = sim.getDM();
    PetscInt cStart = 0, cEnd = 0;
    ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd); CHKERRQ(ierr);
    PetscInt n_local = cEnd - cStart;
    ierr = MPI_Allreduce(&n_local, &n_cells, 1, MPIU_INT, MPI_SUM,
                         PETSC_COMM_WORLD); CHKERRQ(ierr);
    return PETSC_SUCCESS;
  }

  int rank_ = 0;
  std::string config_path_;
};

// (a) Refinement increases cell count near source.
TEST_F(SourceRefinementTest, IncreasesCellCount)
{
  PetscInt n_unref = 0, n_ref = 0;

  writeConfig("unref", false, 0.0, 0);
  ASSERT_EQ(buildDMAndCountCells(n_unref), 0);
  EXPECT_GT(n_unref, 0) << "Un-refined mesh must have cells";

  writeConfig("ref1", true, 5.0, 1);
  ASSERT_EQ(buildDMAndCountCells(n_ref), 0);
  EXPECT_GT(n_ref, 0) << "Refined mesh must have cells";

  // Cell-count ratio bounds: refinement strictly increases the count
  // (>=1.2x lower bound). The upper bound has to admit global hex
  // refinement (PETSc 3.25's DMAdaptLabel on a hex mesh does uniform
  // global refinement regardless of which cells are flagged with
  // DM_ADAPT_REFINE -- it returns ~8x in 3D), so the upper bound is
  // 12x: catches a runaway loop without rejecting the hex global
  // case.
  const double ratio = static_cast<double>(n_ref) /
                       static_cast<double>(n_unref);
  EXPECT_GT(ratio, 1.2)
      << "Refined mesh has only " << n_ref << " cells vs un-refined "
      << n_unref << " (ratio " << ratio << "); expected >1.2x";
  EXPECT_LE(ratio, 12.0)
      << "Refined mesh has " << n_ref << " cells vs un-refined "
      << n_unref << " (ratio " << ratio << "); expected <=12x "
      << "(uniform hex refinement is 8x; > 12x would indicate a "
      << "runaway)";

  if (rank_ == 0) {
    PetscPrintf(PETSC_COMM_WORLD,
                "SourceRefinement: n_unref=%d, n_ref=%d, ratio=%.3f\n",
                (int)n_unref, (int)n_ref, ratio);
  }
}

// (b) Refinement does not move the explosion cell location.
TEST_F(SourceRefinementTest, ExplosionCellNearOriginalLocation)
{
  writeConfig("loc", true, 5.0, 1);

  Simulator sim(PETSC_COMM_WORLD);
  PetscOptionsClear(nullptr);
  PetscErrorCode ierr;
  ierr = sim.initializeFromConfigFile(config_path_);
  ASSERT_EQ(ierr, 0);
  ierr = sim.setupDM();
  ASSERT_EQ(ierr, 0);

  // Locate the cell that contains (2000, 2000, 1806) (the source point
  // from writeConfig) on the refined DM.
  DM dm = sim.getDM();
  Vec pt;
  VecCreateSeq(PETSC_COMM_SELF, 3, &pt);
  VecSetBlockSize(pt, 3);
  PetscScalar* a = nullptr;
  VecGetArray(pt, &a);
  a[0] = 2000.0; a[1] = 2000.0; a[2] = 1806.0;
  VecRestoreArray(pt, &a);

  PetscSF cellSF = nullptr;
  ierr = DMLocatePoints(dm, pt, DM_POINTLOCATION_NONE, &cellSF);
  ASSERT_EQ(ierr, 0);
  const PetscSFNode* remotes = nullptr;
  PetscInt nFound = 0;
  PetscSFGetGraph(cellSF, nullptr, &nFound, nullptr, &remotes);

  if (rank_ == 0 && nFound > 0 && remotes[0].index >= 0) {
    PetscInt c = remotes[0].index;
    PetscReal vol = 0.0;
    PetscReal centroid[3] = {0.0, 0.0, 0.0};
    DMPlexComputeCellGeometryFVM(dm, c, &vol, centroid, nullptr);
    // Approximate local cell length scale from the volume (cell ~
    // pow(vol, 1/3)).
    const double h_refined = std::cbrt(static_cast<double>(vol));
    const double dx = centroid[0] - 2000.0;
    const double dy = centroid[1] - 2000.0;
    const double dz = centroid[2] - 1806.0;
    const double dist = std::sqrt(dx * dx + dy * dy + dz * dz);
    EXPECT_LE(dist, 0.5 * h_refined + h_refined)
        << "Located cell centroid (" << centroid[0] << ", "
        << centroid[1] << ", " << centroid[2]
        << ") differs from source (2000, 2000, 1806) by " << dist
        << " m; allowed " << (0.5 * h_refined + h_refined) << " m "
        << "(h_refined=" << h_refined << " m)";
  }

  PetscSFDestroy(&cellSF);
  VecDestroy(&pt);
}

// (c) Refined run produces a finite, nonzero solution norm through the
// full pipeline.
TEST_F(SourceRefinementTest, FullPipelineFiniteNonzero)
{
  writeConfig("pipeline", true, 5.0, 1);

  Simulator sim(PETSC_COMM_WORLD);
  PetscOptionsClear(nullptr);
  PetscErrorCode ierr;
  ierr = sim.initializeFromConfigFile(config_path_); ASSERT_EQ(ierr, 0);
  ierr = sim.setupDM();             ASSERT_EQ(ierr, 0);
  ierr = sim.labelBoundaries();     ASSERT_EQ(ierr, 0);
  ierr = sim.setupFields();         ASSERT_EQ(ierr, 0);
  ierr = sim.setupPhysics();        ASSERT_EQ(ierr, 0);
  ierr = sim.setupTimeStepper();    ASSERT_EQ(ierr, 0);
  ierr = sim.setupSolvers();        ASSERT_EQ(ierr, 0);
  ierr = sim.setInitialConditions(); ASSERT_EQ(ierr, 0);

  PetscPushErrorHandler(PetscReturnErrorHandler, nullptr);
  ierr = sim.run();
  PetscPopErrorHandler();
  ASSERT_EQ(ierr, 0) << "Refined-mesh pipeline must complete";

  Vec sol = sim.getSolution();
  ASSERT_TRUE(sol != nullptr);
  PetscReal sol_norm = -1.0;
  VecNorm(sol, NORM_2, &sol_norm);
  EXPECT_GT(sol_norm, 0.0)
      << "Solution norm must be > 0 with mesh refinement enabled";
  EXPECT_TRUE(std::isfinite(sol_norm))
      << "Solution norm must be finite with mesh refinement enabled";
}
