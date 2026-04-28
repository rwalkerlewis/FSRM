/**
 * @file test_sedan_1962_mesh.cpp
 * @brief Integration test: load the Sedan 1962 near-field Gmsh mesh
 *
 * This is a *mesh* test, not a *physics* test.  It verifies that the
 * pre-built `meshes/historical/sedan_1962_nearfield.msh` (produced by
 * `scripts/meshing/build_sedan_1962_mesh.sh`) loads through the same
 * DMPlexCreateGmshFromFile path used by Simulator::setupDM and
 * satisfies the size / grading sanity bounds documented in
 * `docs/NEM_ROADMAP.md` Milestone M1.1.
 *
 * Bounds checked (see roadmap M1.1):
 *   - tet count in [50000, 5000000]
 *   - minimum tet edge length in [25, 75] m
 *   - maximum tet edge length in [250, 750] m
 *   - at least one tet whose centroid is within 100 m of the shot
 *     point (0, 0, -194) and whose longest edge is below 100 m,
 *     proving that the high-resolution sphere is anchored on the
 *     shot location.
 */

#include <gtest/gtest.h>
#include <petscsys.h>
#include <petscdmplex.h>

#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace
{

// Resolve the mesh path relative to the build directory used by ctest.
// The integration tests are run from `<repo>/build` (see CLAUDE.md and
// other tests/integration/test_*.cpp), so the relative path matches the
// convention already established by test_gmsh_import.cpp.
constexpr const char* kSedanMeshPath =
    "../../meshes/historical/sedan_1962_nearfield.msh";

constexpr PetscReal kSrcX = 0.0;
constexpr PetscReal kSrcY = 0.0;
constexpr PetscReal kSrcZ = -194.0;

}  // namespace

class Sedan1962MeshTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    PetscErrorCode ierr = DMPlexCreateGmshFromFile(
        PETSC_COMM_WORLD, kSedanMeshPath, PETSC_TRUE, &dm_);
    ASSERT_EQ(ierr, PETSC_SUCCESS)
        << "DMPlexCreateGmshFromFile failed for " << kSedanMeshPath;
    ASSERT_NE(dm_, nullptr);
  }

  void TearDown() override
  {
    if (dm_)
    {
      DMDestroy(&dm_);
      dm_ = nullptr;
    }
  }

  // Returns coordinates of the four vertices of cell c (assumed tet) into
  // pts (flat 12-element array).  Uses DMPlexVecGetClosure on the
  // coordinate Vec, the standard PETSc idiom.
  PetscErrorCode getTetVertices(PetscInt c, std::array<PetscReal, 12>& pts) const
  {
    Vec        coord_vec;
    DM         coord_dm;
    PetscSection coord_section;
    PetscScalar* closure = nullptr;
    PetscInt     csize = 0;

    PetscErrorCode ierr;
    ierr = DMGetCoordinateDM(dm_, &coord_dm); CHKERRQ(ierr);
    ierr = DMGetCoordinatesLocal(dm_, &coord_vec); CHKERRQ(ierr);
    ierr = DMGetCoordinateSection(dm_, &coord_section); CHKERRQ(ierr);
    ierr = DMPlexVecGetClosure(coord_dm, coord_section, coord_vec, c,
                               &csize, &closure); CHKERRQ(ierr);
    if (csize != 12)
    {
      DMPlexVecRestoreClosure(coord_dm, coord_section, coord_vec, c,
                              &csize, &closure);
      return PETSC_ERR_PLIB;
    }
    for (PetscInt i = 0; i < 12; ++i)
    {
      pts[static_cast<size_t>(i)] = PetscRealPart(closure[i]);
    }
    ierr = DMPlexVecRestoreClosure(coord_dm, coord_section, coord_vec, c,
                                   &csize, &closure); CHKERRQ(ierr);
    return PETSC_SUCCESS;
  }

  static PetscReal edgeLength(const PetscReal* a, const PetscReal* b)
  {
    const PetscReal dx = a[0] - b[0];
    const PetscReal dy = a[1] - b[1];
    const PetscReal dz = a[2] - b[2];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
  }

  DM dm_ = nullptr;
};

// Test 1: mesh dimension is 3 and cells are tetrahedra.
TEST_F(Sedan1962MeshTest, MeshDimensionAndCellType)
{
  PetscInt dim = 0;
  ASSERT_EQ(DMGetDimension(dm_, &dim), PETSC_SUCCESS);
  EXPECT_EQ(dim, 3);

  PetscInt cStart = 0, cEnd = 0;
  ASSERT_EQ(DMPlexGetHeightStratum(dm_, 0, &cStart, &cEnd), PETSC_SUCCESS);
  ASSERT_GT(cEnd - cStart, 0);

  DMPolytopeType ct;
  ASSERT_EQ(DMPlexGetCellType(dm_, cStart, &ct), PETSC_SUCCESS);
  EXPECT_EQ(ct, DM_POLYTOPE_TETRAHEDRON)
      << "Sedan near-field mesh must be tetrahedral";
}

// Test 2: tet count and per-tet edge-length bounds.
TEST_F(Sedan1962MeshTest, TetCountAndEdgeLengthBounds)
{
  PetscInt cStart = 0, cEnd = 0;
  ASSERT_EQ(DMPlexGetHeightStratum(dm_, 0, &cStart, &cEnd), PETSC_SUCCESS);
  const PetscInt n_tets = cEnd - cStart;

  // Sanity bounds from docs/NEM_ROADMAP.md M1.1.
  EXPECT_GT(n_tets, 50000)
      << "Sedan near-field mesh has " << n_tets
      << " tets; expected > 50000.  Regenerate via "
         "scripts/meshing/build_sedan_1962_mesh.sh.";
  EXPECT_LT(n_tets, 5000000)
      << "Sedan near-field mesh has " << n_tets
      << " tets; expected < 5000000.";

  PetscReal min_edge = std::numeric_limits<PetscReal>::max();
  PetscReal max_edge = 0.0;

  for (PetscInt c = cStart; c < cEnd; ++c)
  {
    std::array<PetscReal, 12> pts{};
    PetscErrorCode ierr = getTetVertices(c, pts);
    ASSERT_EQ(ierr, PETSC_SUCCESS) << "Failed to read vertices of cell " << c;

    for (int i = 0; i < 4; ++i)
    {
      for (int j = i + 1; j < 4; ++j)
      {
        const PetscReal e = edgeLength(&pts[3 * i], &pts[3 * j]);
        if (e < min_edge) min_edge = e;
        if (e > max_edge) max_edge = e;
      }
    }
  }

  // Minimum edge: 50 m target with +-50% slack for grading.
  EXPECT_GE(min_edge, 25.0) << "Minimum tet edge " << min_edge << " m below 25 m";
  EXPECT_LE(min_edge, 75.0) << "Minimum tet edge " << min_edge << " m above 75 m";

  // Maximum edge: graded boundary cells should not exceed 750 m.
  EXPECT_GE(max_edge, 250.0)
      << "Maximum tet edge " << max_edge << " m below 250 m -- mesh appears uniformly fine";
  EXPECT_LE(max_edge, 750.0)
      << "Maximum tet edge " << max_edge << " m above 750 m -- grading is too coarse";
}

// Test 3: high-resolution sphere is anchored on the shot point.
// At least one tet whose centroid is within 100 m of (0, 0, -194) must
// have a longest edge below 100 m, confirming the local refinement
// took effect at the source location.
TEST_F(Sedan1962MeshTest, FineCellNearSourcePoint)
{
  PetscInt cStart = 0, cEnd = 0;
  ASSERT_EQ(DMPlexGetHeightStratum(dm_, 0, &cStart, &cEnd), PETSC_SUCCESS);

  PetscInt count_fine_near_src = 0;

  for (PetscInt c = cStart; c < cEnd; ++c)
  {
    std::array<PetscReal, 12> pts{};
    PetscErrorCode ierr = getTetVertices(c, pts);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    PetscReal cx = 0.0, cy = 0.0, cz = 0.0;
    for (int i = 0; i < 4; ++i)
    {
      cx += pts[3 * i + 0];
      cy += pts[3 * i + 1];
      cz += pts[3 * i + 2];
    }
    cx *= 0.25; cy *= 0.25; cz *= 0.25;

    const PetscReal dx = cx - kSrcX;
    const PetscReal dy = cy - kSrcY;
    const PetscReal dz = cz - kSrcZ;
    const PetscReal d  = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (d > 100.0) continue;

    PetscReal local_max = 0.0;
    for (int i = 0; i < 4; ++i)
    {
      for (int j = i + 1; j < 4; ++j)
      {
        const PetscReal e = edgeLength(&pts[3 * i], &pts[3 * j]);
        if (e > local_max) local_max = e;
      }
    }
    if (local_max < 100.0)
    {
      ++count_fine_near_src;
    }
  }

  EXPECT_GT(count_fine_near_src, 0)
      << "No fine (max-edge < 100 m) tet found within 100 m of the shot "
         "point (0, 0, -194).  The high-resolution sphere is missing or "
         "misaligned.";
}
