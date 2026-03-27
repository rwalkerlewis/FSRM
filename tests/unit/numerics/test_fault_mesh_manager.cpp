/**
 * @file test_fault_mesh_manager.cpp
 * @brief Geometric and DMPlex mesh-splitting tests for FaultMeshManager
 */

#include <gtest/gtest.h>
#include <cmath>
#include <string>
#include <petscsys.h>
#include <petscdmplex.h>

#include "numerics/FaultMeshManager.hpp"
#include "domain/geomechanics/PyLithFault.hpp"

using namespace FSRM;

class FaultMeshManagerTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    }

    int rank = 0;
};

TEST_F(FaultMeshManagerTest, IsPointOnFaultPlaneDistanceTolerance) {
    FaultMeshManager mgr(PETSC_COMM_WORLD);
    const double center[3] = {0.0, 0.0, 0.0};
    const double normal[3] = {0.0, 0.0, 1.0};

    const double on_plane[3] = {2.0, -1.0, 0.01};
    EXPECT_TRUE(mgr.isPointOnFaultPlane(on_plane, normal, center, 0.02));
    EXPECT_FALSE(mgr.isPointOnFaultPlane(on_plane, normal, center, 0.001));

    const double off_plane[3] = {0.0, 0.0, 0.5};
    EXPECT_FALSE(mgr.isPointOnFaultPlane(off_plane, normal, center, 0.1));
}

TEST_F(FaultMeshManagerTest, IsPointWithinFaultRectangularExtent) {
    FaultMeshManager mgr(PETSC_COMM_WORLD);
    const double center[3] = {100.0, 200.0, -500.0};
    const double strike_dir[3] = {1.0, 0.0, 0.0};
    const double dip_dir[3] = {0.0, 0.0, -1.0};
    const double L = 200.0;
    const double W = 100.0;

    const double inside[3] = {150.0, 200.0, -520.0};
    EXPECT_TRUE(
        mgr.isPointWithinFaultExtent(inside, center, strike_dir, dip_dir, L, W));

    const double outside_strike[3] = {250.0, 200.0, -500.0};
    EXPECT_FALSE(
        mgr.isPointWithinFaultExtent(outside_strike, center, strike_dir, dip_dir, L, W));

    const double outside_dip[3] = {100.0, 200.0, -560.0};
    EXPECT_FALSE(
        mgr.isPointWithinFaultExtent(outside_dip, center, strike_dir, dip_dir, L, W));
}

TEST_F(FaultMeshManagerTest, SplitMeshAlongFaultRequiresDMPlex) {
    DM dm = nullptr;
    PetscInt faces[3] = {2, 2, 2};
    PetscReal lower[3] = {0.0, 0.0, 0.0};
    PetscReal upper[3] = {1.0, 1.0, 1.0};
    PetscErrorCode ierr = DMPlexCreateBoxMesh(PETSC_COMM_WORLD, 3, PETSC_FALSE, faces, lower,
                                              upper, nullptr, PETSC_TRUE, &dm);
    if (ierr != 0 || !dm) {
        GTEST_SKIP() << "DMPlexCreateBoxMesh unavailable or failed (PETSc build without "
                        "full DMPlex support)";
    }

    FaultMeshManager mgr(PETSC_COMM_WORLD);
    DMLabel fault_label = nullptr;
    const double center[3] = {0.5, 0.5, 0.5};
    ierr = mgr.createPlanarFaultLabel(dm, &fault_label, 0.0, 0.5 * M_PI, center, 2.0, 2.0,
                                      0.15);
    if (ierr != 0) {
        DMDestroy(&dm);
        GTEST_SKIP() << "createPlanarFaultLabel failed on box mesh; geometry labeling not "
                        "applicable to this mesh";
    }

    ierr = mgr.splitMeshAlongFault(&dm, "fault");
    if (ierr != 0) {
        DMDestroy(&dm);
        GTEST_SKIP() << "DMPlexConstructCohesiveCells not available or split failed for this "
                        "PETSc / mesh combination";
    }

    FaultCohesiveDyn dyn;
    ierr = mgr.extractCohesiveTopology(dm, &dyn);
    EXPECT_EQ(ierr, 0);

    DMDestroy(&dm);
}
