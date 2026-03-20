/**
 * @file test_boundary_conditions.cpp
 * @brief Smoke tests for BoundaryConditionManager and Dirichlet/Neumann BC types
 */

#include <gtest/gtest.h>
#include "numerics/BoundaryConditions.hpp"
#include <memory>
#include <mpi.h>

using namespace FSRM;

class BoundaryConditionManagerTest : public ::testing::Test {
protected:
    void SetUp() override { MPI_Comm_rank(PETSC_COMM_WORLD, &rank); }
    int rank = 0;
};

TEST_F(BoundaryConditionManagerTest, InstantiateManager) {
    BoundaryConditionManager mgr;
    (void)mgr;
}

TEST_F(BoundaryConditionManagerTest, AddDirichletAndNeumannViaBoundaryCondition) {
    BoundaryConditionManager mgr;
    mgr.addBoundaryCondition(1, std::make_unique<DirichletBC>());
    mgr.addBoundaryCondition(2, std::make_unique<NeumannBC>());
    SUCCEED();
}
