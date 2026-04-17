/**
 * @file test_cohesive_interface_label.cpp
 * @brief Phase 1 unit test for FSRM::faults::topology::createInterfacesLabel.
 *
 * Builds a minimal 1D DMPlex with two bulk segments, splits it along the
 * interior vertex using DMPlexConstructCohesiveCells, and verifies that
 * the "cohesive interface" label is created with a stratum covering the
 * expected cohesive closure points.
 */

#include <gtest/gtest.h>

#include <petscsys.h>
#include <petscdm.h>
#include <petscdmplex.h>

#include "faults/FaultCohesive.hpp"
#include "faults/TopologyOps.hpp"

namespace {

class CohesiveInterfaceLabelTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    }

    /**
     * Build a 2D two-triangle DM with a fault label on the shared interior
     * edge, then call DMPlexConstructCohesiveCells to produce a DM
     * containing one cohesive (hybrid) cell. Returned DM is owned by the
     * caller. DMPlexConstructCohesiveCells requires dim >= 2.
     */
    PetscErrorCode buildTwoCellDMWithCohesive(DM *dm_out) {
        PetscFunctionBeginUser;
        PetscErrorCode ierr;

        // 2D plex: two triangles sharing edge v0-v2.
        //   v3(0,1) --- v2(1,1)
        //     |   \      |
        //     |    \     |
        //     |     \    |
        //   v0(0,0) --- v1(1,0)
        DM dm_base = nullptr;
        const PetscInt dim = 2;
        const PetscInt numCells = 2;
        const PetscInt numVertices = 4;
        const PetscInt numCorners = 3;
        PetscInt cells[6] = {0, 1, 2,  0, 2, 3};
        PetscReal coords[8] = {0.0, 0.0,  1.0, 0.0,  1.0, 1.0,  0.0, 1.0};

        ierr = DMPlexCreateFromCellListPetsc(PETSC_COMM_SELF,
                                             dim,
                                             numCells,
                                             numVertices,
                                             numCorners,
                                             PETSC_FALSE /*interpolate*/,
                                             cells,
                                             dim /*spaceDim*/,
                                             coords,
                                             &dm_base); CHKERRQ(ierr);

        // Interpolate so edges exist as a stratum
        DM dm_interp = nullptr;
        ierr = DMPlexInterpolate(dm_base, &dm_interp); CHKERRQ(ierr);
        ierr = DMDestroy(&dm_base); CHKERRQ(ierr);

        // Locate the interior edge: the only edge whose support contains
        // both cells. Walk the height-1 stratum (edges in 2D).
        PetscInt eStart = 0, eEnd = 0;
        ierr = DMPlexGetHeightStratum(dm_interp, 1, &eStart, &eEnd); CHKERRQ(ierr);

        PetscInt interior_edge = -1;
        for (PetscInt e = eStart; e < eEnd; ++e) {
            PetscInt supportSize = 0;
            ierr = DMPlexGetSupportSize(dm_interp, e, &supportSize); CHKERRQ(ierr);
            if (supportSize == 2) {
                interior_edge = e;
                break;
            }
        }
        PetscCheck(interior_edge >= 0, PETSC_COMM_SELF, PETSC_ERR_PLIB,
                   "No interior edge found in two-triangle plex");

        // Label the interior edge as fault.
        ierr = DMCreateLabel(dm_interp, "fault"); CHKERRQ(ierr);
        DMLabel faultLabel = nullptr;
        ierr = DMGetLabel(dm_interp, "fault", &faultLabel); CHKERRQ(ierr);
        ierr = DMLabelSetValue(faultLabel, interior_edge, 1); CHKERRQ(ierr);

        // Follow the PyLith cohesive-cell workflow used by
        // FaultMeshManager::splitMeshAlongFault (src/numerics/FaultMeshManager.cpp):
        //   (a) duplicate + complete the fault label
        //   (b) build a fault submesh
        //   (c) orient and extract the subpoint map
        //   (d) derive an oriented cohesive label, clear cell stratum
        //   (e) LabelCohesiveComplete with the fault submesh
        //   (f) ConstructCohesiveCells
        DMLabel completed_label = nullptr;
        ierr = DMLabelDuplicate(faultLabel, &completed_label); CHKERRQ(ierr);
        ierr = DMPlexLabelComplete(dm_interp, completed_label); CHKERRQ(ierr);

        DM fault_dm = nullptr;
        ierr = DMPlexCreateSubmesh(dm_interp, completed_label, 1, PETSC_TRUE,
                                   &fault_dm); CHKERRQ(ierr);
        ierr = DMLabelDestroy(&completed_label); CHKERRQ(ierr);
        ierr = DMPlexOrient(fault_dm); CHKERRQ(ierr);

        DMLabel subpoint_map = nullptr;
        ierr = DMPlexGetSubpointMap(fault_dm, &subpoint_map); CHKERRQ(ierr);

        DMLabel cohesive_label = nullptr;
        ierr = DMLabelDuplicate(subpoint_map, &cohesive_label); CHKERRQ(ierr);
        ierr = DMLabelClearStratum(cohesive_label, dim); CHKERRQ(ierr);
        ierr = DMPlexOrientLabel(dm_interp, cohesive_label); CHKERRQ(ierr);
        ierr = DMPlexLabelCohesiveComplete(dm_interp, cohesive_label, NULL, 0,
                                           PETSC_FALSE, PETSC_FALSE, fault_dm);
        CHKERRQ(ierr);

        DM dm_split = nullptr;
        ierr = DMPlexConstructCohesiveCells(dm_interp, cohesive_label, nullptr,
                                            &dm_split); CHKERRQ(ierr);

        ierr = DMLabelDestroy(&cohesive_label); CHKERRQ(ierr);
        ierr = DMDestroy(&fault_dm); CHKERRQ(ierr);
        ierr = DMDestroy(&dm_interp); CHKERRQ(ierr);

        *dm_out = dm_split;
        PetscFunctionReturn(PETSC_SUCCESS);
    }

    int rank = 0;
};

TEST_F(CohesiveInterfaceLabelTest, CreatesLabelWithNonEmptyStratumOnTwoCellFault) {
    DM dm = nullptr;
    ASSERT_EQ(buildTwoCellDMWithCohesive(&dm), PETSC_SUCCESS);
    ASSERT_NE(dm, nullptr);

    // Label must not exist before the call
    PetscBool hasBefore = PETSC_TRUE;
    ASSERT_EQ(DMHasLabel(dm, "cohesive interface", &hasBefore), PETSC_SUCCESS);
    EXPECT_EQ(hasBefore, PETSC_FALSE);

    DMLabel label = nullptr;
    ASSERT_EQ(FSRM::faults::topology::createInterfacesLabel(dm, nullptr, &label),
              PETSC_SUCCESS);
    EXPECT_NE(label, nullptr);

    // Label must now exist via DMGetLabel as well
    DMLabel fetched = nullptr;
    ASSERT_EQ(DMGetLabel(dm, "cohesive interface", &fetched), PETSC_SUCCESS);
    EXPECT_NE(fetched, nullptr);
    EXPECT_EQ(fetched, label);

    // The stratum with value 1 must contain at least one point. In the 1D
    // two-cell-plus-cohesive geometry, the cohesive cell itself plus its
    // closure points (two bounding vertices in the duplicated layer)
    // produce at least one labelled point at each of height 0 and height 1.
    IS stratumIS = nullptr;
    ASSERT_EQ(DMLabelGetStratumIS(label, 1, &stratumIS), PETSC_SUCCESS);
    ASSERT_NE(stratumIS, nullptr);

    PetscInt nLabelled = 0;
    ASSERT_EQ(ISGetLocalSize(stratumIS, &nLabelled), PETSC_SUCCESS);
    EXPECT_GT(nLabelled, 0);
    ASSERT_EQ(ISDestroy(&stratumIS), PETSC_SUCCESS);

    // There must be exactly one hybrid cell (cohesive prism of the two
    // bulk segments) in the DM at height 0. DMPlexGetSimplexOrBoxCells
    // reports pMax as the index where bulk ends and hybrid begins, so
    // pEnd - pMax = 1 for our geometry.
    PetscInt pStart = 0, pEnd = 0, pMax = 0;
    ASSERT_EQ(DMPlexGetHeightStratum(dm, 0, &pStart, &pEnd), PETSC_SUCCESS);
    ASSERT_EQ(DMPlexGetSimplexOrBoxCells(dm, 0, NULL, &pMax), PETSC_SUCCESS);
    EXPECT_EQ(pEnd - pMax, 1);

    ASSERT_EQ(DMDestroy(&dm), PETSC_SUCCESS);
}

TEST_F(CohesiveInterfaceLabelTest, IsIdempotentAcrossMultipleCalls) {
    DM dm = nullptr;
    ASSERT_EQ(buildTwoCellDMWithCohesive(&dm), PETSC_SUCCESS);

    DMLabel first = nullptr;
    ASSERT_EQ(FSRM::faults::topology::createInterfacesLabel(dm, nullptr, &first),
              PETSC_SUCCESS);
    ASSERT_NE(first, nullptr);

    DMLabel second = nullptr;
    ASSERT_EQ(FSRM::faults::topology::createInterfacesLabel(dm, nullptr, &second),
              PETSC_SUCCESS);
    EXPECT_EQ(first, second);

    ASSERT_EQ(DMDestroy(&dm), PETSC_SUCCESS);
}

// FaultCohesive is abstract; declare a minimal concrete subclass to exercise
// the createCohesiveInterfaceLabel wrapper and the accessor.
class DummyFault : public FSRM::faults::FaultCohesive {
public:
    PetscErrorCode initialize(DM, const FSRM::FaultConfig &) override {
        return PETSC_SUCCESS;
    }
};

TEST_F(CohesiveInterfaceLabelTest, FaultCohesiveWrapperCachesLabelPointer) {
    DM dm = nullptr;
    ASSERT_EQ(buildTwoCellDMWithCohesive(&dm), PETSC_SUCCESS);

    DummyFault fault;
    EXPECT_EQ(fault.getCohesiveInterfaceLabel(), nullptr);

    ASSERT_EQ(fault.createCohesiveInterfaceLabel(dm, nullptr), PETSC_SUCCESS);
    EXPECT_NE(fault.getCohesiveInterfaceLabel(), nullptr);

    // Accessor must match what DMGetLabel returns for "cohesive interface"
    DMLabel fromDM = nullptr;
    ASSERT_EQ(DMGetLabel(dm, "cohesive interface", &fromDM), PETSC_SUCCESS);
    EXPECT_EQ(fault.getCohesiveInterfaceLabel(), fromDM);

    ASSERT_EQ(DMDestroy(&dm), PETSC_SUCCESS);
}

}  // namespace
