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

    // Wrappers for private methods (friend access doesn't extend to TEST_F derived classes)
    static bool callIsPointOnFaultPlane(const FaultMeshManager& mgr,
                                        const double point[3], const double normal[3],
                                        const double center[3], double tolerance) {
        return mgr.isPointOnFaultPlane(point, normal, center, tolerance);
    }

    static bool callIsPointWithinFaultExtent(const FaultMeshManager& mgr,
                                              const double point[3], const double center[3],
                                              const double strike_dir[3], const double dip_dir[3],
                                              double length, double width) {
        return mgr.isPointWithinFaultExtent(point, center, strike_dir, dip_dir, length, width);
    }

    int rank = 0;
};

TEST_F(FaultMeshManagerTest, IsPointOnFaultPlaneDistanceTolerance) {
    FaultMeshManager mgr(PETSC_COMM_WORLD);
    const double center[3] = {0.0, 0.0, 0.0};
    const double normal[3] = {0.0, 0.0, 1.0};

    const double on_plane[3] = {2.0, -1.0, 0.01};
    EXPECT_TRUE(callIsPointOnFaultPlane(mgr, on_plane, normal, center, 0.02));
    EXPECT_FALSE(callIsPointOnFaultPlane(mgr, on_plane, normal, center, 0.001));

    const double off_plane[3] = {0.0, 0.0, 0.5};
    EXPECT_FALSE(callIsPointOnFaultPlane(mgr, off_plane, normal, center, 0.1));
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
        callIsPointWithinFaultExtent(mgr, inside, center, strike_dir, dip_dir, L, W));

    const double outside_strike[3] = {250.0, 200.0, -500.0};
    EXPECT_FALSE(
        callIsPointWithinFaultExtent(mgr, outside_strike, center, strike_dir, dip_dir, L, W));

    const double outside_dip[3] = {100.0, 200.0, -560.0};
    EXPECT_FALSE(
        callIsPointWithinFaultExtent(mgr, outside_dip, center, strike_dir, dip_dir, L, W));
}

TEST_F(FaultMeshManagerTest, SplitMeshAlongFaultRequiresDMPlex) {
    DM dm = nullptr;
    PetscInt faces[3] = {4, 4, 4};
    PetscReal lower[3] = {0.0, 0.0, 0.0};
    PetscReal upper[3] = {1.0, 1.0, 1.0};
    // Use simplex mesh (PETSC_TRUE) - cohesive cells work reliably on simplex in PETSc 3.22
    PetscErrorCode ierr = DMPlexCreateBoxMesh(PETSC_COMM_WORLD, 3, PETSC_TRUE, faces, lower,
                                              upper, nullptr, PETSC_TRUE, 0, PETSC_FALSE, &dm);
    ASSERT_EQ(ierr, 0) << "DMPlexCreateBoxMesh must succeed in Docker CI";
    ASSERT_NE(dm, nullptr);

    FaultMeshManager mgr(PETSC_COMM_WORLD);
    DMLabel fault_label = nullptr;
    const double center[3] = {0.5, 0.5, 0.5};
    // Simplex mesh face centroids are not at exact grid coords, use smaller tolerance
    ierr = mgr.createPlanarFaultLabel(dm, &fault_label, 0.0, 0.5 * M_PI, center, 2.0, 2.0,
                                      0.05);
    ASSERT_EQ(ierr, 0) << "createPlanarFaultLabel must succeed";

    ierr = mgr.splitMeshAlongFault(&dm, "fault");
    ASSERT_EQ(ierr, 0) << "splitMeshAlongFault must succeed";

    FaultCohesiveDyn dyn;
    ierr = mgr.extractCohesiveTopology(dm, &dyn);
    EXPECT_EQ(ierr, 0);

    DMDestroy(&dm);
}
