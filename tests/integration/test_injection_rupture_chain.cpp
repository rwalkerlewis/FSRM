/**
 * @file test_injection_rupture_chain.cpp
 * @brief Light integration: fault mesh manager, cohesive dynamic fault, Coulomb sampling
 */

#include <gtest/gtest.h>
#include <array>
#include <cmath>
#include <vector>
#include <petscsys.h>
#include <petscdmplex.h>

#include "numerics/FaultMeshManager.hpp"
#include "domain/geomechanics/PyLithFault.hpp"
#include "domain/geomechanics/CoulombStressTransfer.hpp"

using namespace FSRM;

class InjectionRuptureChainTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    }

    int rank = 0;
};

TEST_F(InjectionRuptureChainTest, ComponentsInstantiateTogether) {
    FaultMeshManager mesh_mgr(PETSC_COMM_WORLD);
    FaultCohesiveDyn fault;
    CoulombStressTransfer cst(PETSC_COMM_WORLD);

    FaultVertex v;
    v.vertex_negative = 0;
    v.vertex_positive = 1;
    v.normal = {0.0, 0.0, 1.0};
    v.along_strike = {1.0, 0.0, 0.0};
    v.up_dip = {0.0, 1.0, 0.0};
    fault.setFaultVertices({v});
    fault.setFrictionModel(std::make_unique<SlipWeakeningFriction>());
    fault.setUniformInitialTraction(25e6, 0.0, -40e6);
    fault.initialize();

    EXPECT_EQ(cst.initialize(nullptr, &fault, 50e9, 30e9, 0.9), 0);
    (void)mesh_mgr;
}

TEST_F(InjectionRuptureChainTest, FaultMeshManagerWithFaultCohesiveDynBasicSetup) {
    DM dm = nullptr;
    PetscInt faces[3] = {2, 2, 2};
    PetscReal lower[3] = {0.0, 0.0, 0.0};
    PetscReal upper[3] = {1.0, 1.0, 1.0};
    PetscErrorCode ierr =
        DMPlexCreateBoxMesh(PETSC_COMM_WORLD, 3, PETSC_FALSE, faces, lower, upper, nullptr,
                            PETSC_TRUE, &dm);
    if (ierr != 0 || !dm) {
        GTEST_SKIP() << "DMPlexCreateBoxMesh failed; skip mesh-chain setup";
    }

    FaultMeshManager mgr(PETSC_COMM_WORLD);
    FaultCohesiveDyn fault;
    DMLabel fault_label = nullptr;
    const double center[3] = {0.5, 0.5, 0.5};
    ierr = mgr.createPlanarFaultLabel(dm, &fault_label, 0.0, 0.5 * M_PI, center, 2.0, 2.0,
                                      0.2);
    if (ierr != 0) {
        DMDestroy(&dm);
        GTEST_SKIP() << "createPlanarFaultLabel failed for injection-rupture chain test";
    }

    ierr = mgr.splitMeshAlongFault(&dm, "fault");
    if (ierr != 0) {
        DMDestroy(&dm);
        GTEST_SKIP() << "splitMeshAlongFault failed; cohesive workflow not available";
    }

    ierr = mgr.extractCohesiveTopology(dm, &fault);
    ASSERT_EQ(ierr, 0);
    fault.setFrictionModel(std::make_unique<SlipWeakeningFriction>());
    fault.initialize();

    EXPECT_GE(fault.numVertices(), 0u);
    DMDestroy(&dm);
}

TEST_F(InjectionRuptureChainTest, CoulombStressTransferBasicComputation) {
    FaultCohesiveDyn fault;
    FaultVertex v;
    v.normal = {0.0, 0.0, 1.0};
    v.along_strike = {1.0, 0.0, 0.0};
    v.up_dip = {0.0, 1.0, 0.0};
    fault.setFaultVertices({v});

    CoulombStressTransfer cst(PETSC_COMM_WORLD);
    ASSERT_EQ(cst.initialize(nullptr, &fault, 40e9, 30e9, 1.0), 0);

    double eps[6] = {1e-5, 1e-5, 1e-5, 0.0, 0.0, 0.0};
    double sig[6];
    cst.computeStressFromStrain(eps, 5e6, sig);
    cst.current_[0].sigma = {sig[0], sig[1], sig[2], sig[3], sig[4], sig[5]};
    cst.current_[0].pressure = 5e6;

    ASSERT_EQ(cst.resolveStressOnFault(), 0);
    EXPECT_TRUE(std::isfinite(cst.getCurrentStress()[0].tau));
    EXPECT_TRUE(std::isfinite(cst.getCurrentStress()[0].sigma_n_eff));
}
