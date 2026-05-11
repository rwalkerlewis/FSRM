/**
 * @file test_source3d_ball_mesh.cpp
 * @brief Pass-13a (axis-1b foundation) unit tests for Source3DBallMesh.
 *        Verifies TetGen .node/.ele parsing, DMPlex construction, and
 *        DMPlexDistribute partitioning.
 *
 * Fixture meshes:
 *  - tests/data/source_ball/single_tet/  (1 tet, 4 nodes)
 *  - tests/data/source_ball/cube_5tet/   (5 tets, 8 nodes)
 *
 * Both fixtures are hand-written and do not require TetGen to
 * regenerate.
 */

#include <gtest/gtest.h>

#include <petscdmplex.h>
#include <petscdmlabel.h>
#include <petscsys.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "domain/explosion/Source3DBallMesh.hpp"
#include "domain/explosion/Source3DBallImpl.hpp"
#include "domain/explosion/Source3DBall.hpp"

using FSRM::Source3DBallMesh;
using FSRM::Source3DBallImpl;
using FSRM::Source3DBallConfig;
using FSRM::makeSource3DBall;

namespace
{

/// Resolve the fixture mesh basename relative to the source tree.
/// Tests are typically run from build/, so the fixtures live one
/// directory up. Honour FSRM_TEST_DATA_DIR for out-of-tree builds.
std::string fixtureBasename(const std::string& fixture_name)
{
    const char* env = std::getenv("FSRM_TEST_DATA_DIR");
    std::string root;
    if (env && *env) {
        root = env;
    } else {
        // Search a few standard locations so the test works whether
        // the binary is invoked from build/, build/tests/, or the
        // repo root.
        std::vector<std::string> candidates = {
            "../tests/data/source_ball",
            "../../tests/data/source_ball",
            "tests/data/source_ball",
            "./tests/data/source_ball",
        };
        for (const auto& c : candidates) {
            std::filesystem::path p(c);
            if (std::filesystem::exists(p)) {
                root = c;
                break;
            }
        }
        if (root.empty()) {
            root = "../tests/data/source_ball";  // best-effort default
        }
    }
    return root + "/" + fixture_name + "/source_ball";
}

}  // namespace

class SourceBallMeshTest : public ::testing::Test
{
};

TEST_F(SourceBallMeshTest, MeshLoadFromTetGen_SingleTet)
{
    Source3DBallMesh mesh;
    EXPECT_NO_THROW(mesh.loadFromTetGen(PETSC_COMM_WORLD,
                                        fixtureBasename("single_tet")));

    EXPECT_EQ(mesh.numGlobalCells(), 1);
    EXPECT_NE(mesh.getDM(), nullptr);

    // Local cells: at MPI=1 the lone rank owns the cell. At MPI>1
    // exactly one rank owns it (DMPlexDistribute with overlap=0).
    PetscInt local_total = 0;
    PetscInt my_local = mesh.numLocalCells();
    MPI_Allreduce(&my_local, &local_total, 1, MPIU_INT, MPI_SUM,
                  PETSC_COMM_WORLD);
    EXPECT_EQ(local_total, 1)
        << "Sum of locally-owned cells across ranks must equal 1.";
}

TEST_F(SourceBallMeshTest, MeshLoadFromTetGen_Cube5Tet)
{
    Source3DBallMesh mesh;
    ASSERT_NO_THROW(mesh.loadFromTetGen(PETSC_COMM_WORLD,
                                        fixtureBasename("cube_5tet")));

    EXPECT_EQ(mesh.numGlobalCells(), 5);

    PetscInt my_local = mesh.numLocalCells();
    PetscInt local_total = 0;
    MPI_Allreduce(&my_local, &local_total, 1, MPIU_INT, MPI_SUM,
                  PETSC_COMM_WORLD);
    EXPECT_EQ(local_total, 5)
        << "Sum of locally-owned cells across ranks must equal 5.";
}

TEST_F(SourceBallMeshTest, VertexMarkerLabelSurvivesDistribute)
{
    Source3DBallMesh mesh;
    ASSERT_NO_THROW(mesh.loadFromTetGen(PETSC_COMM_WORLD,
                                        fixtureBasename("cube_5tet")));

    // The cube_5tet fixture marks vertex 1 (origin) as
    // INNER_CAVITY_SURFACE (marker 2) and the other 7 vertices as
    // OUTER_ELASTIC_SURFACE (marker 1). After distribution the global
    // counts must reduce to those targets.
    PetscInt my_cavity = mesh.numCavityMarkedVertices();
    PetscInt my_elastic = mesh.numElasticMarkedVertices();
    PetscInt total_cavity = 0;
    PetscInt total_elastic = 0;
    MPI_Allreduce(&my_cavity, &total_cavity, 1, MPIU_INT, MPI_SUM,
                  PETSC_COMM_WORLD);
    MPI_Allreduce(&my_elastic, &total_elastic, 1, MPIU_INT, MPI_SUM,
                  PETSC_COMM_WORLD);

    PetscMPIInt size = 1;
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    if (size == 1) {
        EXPECT_EQ(total_cavity, 1);
        EXPECT_EQ(total_elastic, 7);
    } else {
        // Under MPI distribution, vertices on the partition boundary
        // are duplicated, so the sum can exceed the global count.
        // The cavity vertex must appear at least once; the elastic
        // count must be at least the global elastic vertex count.
        EXPECT_GE(total_cavity, 1);
        EXPECT_GE(total_elastic, 7);
    }
}

TEST_F(SourceBallMeshTest, GeometryRadiiSane)
{
    Source3DBallMesh mesh;
    ASSERT_NO_THROW(mesh.loadFromTetGen(PETSC_COMM_WORLD,
                                        fixtureBasename("cube_5tet")));

    // Cube_5tet fixture: origin vertex at radius 0, far corner at
    // sqrt(3) ~ 1.7321. Local extrema vary by partition; reduce
    // across ranks for a global check. Empty ranks contribute neutral
    // values (+inf for min, 0 for max).
    double local_min = mesh.localMinVertexRadius();
    double local_max = mesh.localMaxVertexRadius();
    if (mesh.numLocalVertices() == 0) {
        local_min = std::numeric_limits<double>::infinity();
        local_max = 0.0;
    }
    double global_min = local_min;
    double global_max = local_max;
    MPI_Allreduce(&local_min, &global_min, 1, MPI_DOUBLE, MPI_MIN,
                  PETSC_COMM_WORLD);
    MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX,
                  PETSC_COMM_WORLD);

    EXPECT_NEAR(global_min, 0.0, 1.0e-9);
    EXPECT_NEAR(global_max, std::sqrt(3.0), 1.0e-6);
}

TEST_F(SourceBallMeshTest, MissingFileThrows)
{
    Source3DBallMesh mesh;
    EXPECT_THROW(
        mesh.loadFromTetGen(PETSC_COMM_WORLD,
                            "tests/data/source_ball/__no_such_basename__"),
        std::runtime_error);
}

// =============================================================================
// Source3DBallImpl skeleton
// =============================================================================

class SourceBallImplFoundationTest : public ::testing::Test
{
};

TEST_F(SourceBallImplFoundationTest, FactoryReturnsValidImpl)
{
    Source3DBallConfig cfg;
    cfg.outer_radius_m = 1.0;
    cfg.cavity_radius_m = 0.0;  // cube fixture has origin vertex
    cfg.mesh_path = "";  // empty: skip mesh load

    auto ball = makeSource3DBall(cfg);
    ASSERT_NE(ball, nullptr);
    const std::string nm = ball->name();
    // Pass-13c: physics tag bumps to the validation tag. The
    // "Source3DBallImpl" prefix and a "pass13c" substring must appear
    // so callers can regression-check the implementation version.
    EXPECT_NE(nm.find("pass13c"), std::string::npos) << nm;
    EXPECT_NE(nm.find("Source3DBallImpl"), std::string::npos) << nm;
}

TEST_F(SourceBallImplFoundationTest, InitializeWithEmptyMeshPathSucceeds)
{
    Source3DBallConfig cfg;
    cfg.mesh_path = "";  // foundation no-op mode

    Source3DBallImpl ball;
    EXPECT_NO_THROW(ball.initialize(cfg));
    EXPECT_FALSE(ball.isInitialized());
    EXPECT_EQ(ball.mesh(), nullptr);
}

TEST_F(SourceBallImplFoundationTest, InitializeWithFixtureMeshSucceeds)
{
    Source3DBallConfig cfg;
    cfg.mesh_path = fixtureBasename("cube_5tet");

    Source3DBallImpl ball;
    EXPECT_NO_THROW(ball.initialize(cfg));
    EXPECT_TRUE(ball.isInitialized());
    ASSERT_NE(ball.mesh(), nullptr);
    EXPECT_EQ(ball.mesh()->numGlobalCells(), 5);
}

TEST_F(SourceBallImplFoundationTest, StepThrowsInFoundationNoOpMode)
{
    // Pass-13b: step() works on a real loaded mesh. In the foundation
    // no-op mode (empty mesh_path) step still throws because there is
    // no per-cell state allocated. The throw message must reference
    // mesh_path so callers know to populate it.
    Source3DBallConfig cfg;
    cfg.mesh_path = "";
    Source3DBallImpl ball;
    ball.initialize(cfg);

    bool caught = false;
    std::string what;
    try {
        ball.step(1.0e-6);
    } catch (const std::runtime_error& e) {
        caught = true;
        what = e.what();
    }
    EXPECT_TRUE(caught);
    EXPECT_NE(what.find("mesh_path"), std::string::npos) << what;
}

TEST_F(SourceBallImplFoundationTest, StepWithFixtureMeshAdvancesState)
{
    // Pass-13b: with a real loaded mesh, step() returns without
    // throwing and advances the per-cell state. We exercise the
    // smallest fixture so the gate stays fast.
    Source3DBallConfig cfg;
    cfg.mesh_path = fixtureBasename("cube_5tet");
    cfg.asymmetric_overburden = false;  // skip the IC for this gate
    cfg.host_density_kg_per_m3 = 2700.0;
    cfg.bulk_modulus_K_pa = 30.0e9;
    cfg.shear_modulus_G_pa = 18.0e9;
    cfg.dp3d_alpha = 0.30;
    cfg.dp3d_k_pa = 70.0e6;
    Source3DBallImpl ball;
    EXPECT_NO_THROW(ball.initialize(cfg));
    EXPECT_TRUE(ball.isInitialized());

    EXPECT_NO_THROW(ball.step(1.0e-9));
    EXPECT_GT(ball.numLocalCells(), 0);
    // No strain rate set -> stress unchanged from IC (zeros here).
    const auto& sigma = ball.sigmaCells();
    for (const auto& s : sigma) {
        for (int i = 0; i < 6; ++i) {
            EXPECT_EQ(s[i], 0.0);
        }
    }
}

TEST_F(SourceBallImplFoundationTest, MomentTensorIsZeroBeforeStep)
{
    Source3DBallConfig cfg;
    cfg.mesh_path = "";
    Source3DBallImpl ball;
    ball.initialize(cfg);

    std::array<double, 6> M{};
    ball.getMomentTensor(M);
    for (int i = 0; i < 6; ++i) {
        EXPECT_EQ(M[i], 0.0);
    }
    std::array<double, 6> Mdot{};
    ball.getMomentRateTensor(Mdot);
    for (int i = 0; i < 6; ++i) {
        EXPECT_EQ(Mdot[i], 0.0);
    }
}

TEST_F(SourceBallImplFoundationTest, NodalFEMRadiationDisciplineThrows)
{
    Source3DBallConfig cfg;
    cfg.mesh_path = "";
    cfg.radiation_discretization =
        Source3DBallConfig::RadiationDiscretization::FEM_NODAL;

    Source3DBallImpl ball;
    bool caught = false;
    std::string what;
    try {
        ball.initialize(cfg);
    } catch (const std::runtime_error& e) {
        caught = true;
        what = e.what();
    }
    EXPECT_TRUE(caught);
    EXPECT_NE(what.find("pass-15"), std::string::npos) << what;
}

// =============================================================================
// BackwardCompat: pass-11 throw is gone on factory instantiation
// =============================================================================

TEST(BackwardCompatTest, Source3DBallScaffoldThrowGoneOnInstantiation)
{
    Source3DBallConfig cfg;
    EXPECT_NO_THROW({
        auto ball = makeSource3DBall(cfg);
        ASSERT_NE(ball, nullptr);
    }) << "Pass-11 scaffold's throw-on-construct must be replaced; "
          "pass-13a foundation returns a real Source3DBallImpl.";
}
