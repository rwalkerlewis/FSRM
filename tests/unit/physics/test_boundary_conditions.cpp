/**
 * @file test_boundary_conditions.cpp
 * @brief Unit tests for BoundaryConditionManager and BC types
 *
 * Tests:
 *   1. Manager construction
 *   2. Add Dirichlet and Neumann BCs (no crash)
 *   3. DirichletBC and NeumannBC can be constructed
 *   4. Manager PML flag starts false
 *   5. Multiple BC registrations do not crash
 */

#include <gtest/gtest.h>
#include "numerics/BoundaryConditions.hpp"
#include <petsc.h>
#include <memory>
#include <mpi.h>

using namespace FSRM;

class BoundaryConditionManagerTest : public ::testing::Test {
protected:
    void SetUp() override { MPI_Comm_rank(PETSC_COMM_WORLD, &rank); }
    int rank = 0;
};

// Test 1: Manager construction
TEST_F(BoundaryConditionManagerTest, InstantiateManager)
{
  BoundaryConditionManager mgr;
  EXPECT_FALSE(mgr.hasPML()) << "New manager should not have PML by default";
}

// Test 2: Add Dirichlet and Neumann BCs without crash
TEST_F(BoundaryConditionManagerTest, AddBCsNoCrash)
{
  BoundaryConditionManager mgr;
  mgr.addBoundaryCondition(1, std::make_unique<DirichletBC>());
  mgr.addBoundaryCondition(2, std::make_unique<NeumannBC>());
  // No crash means both BC types can be registered
  EXPECT_FALSE(mgr.hasPML()) << "PML should still be false after adding BCs";
}

// Test 3: DirichletBC and NeumannBC construction
TEST_F(BoundaryConditionManagerTest, BCTypeConstruction)
{
  auto dirichlet = std::make_unique<DirichletBC>();
  auto neumann = std::make_unique<NeumannBC>();
  EXPECT_NE(dirichlet, nullptr) << "DirichletBC should construct successfully";
  EXPECT_NE(neumann, nullptr) << "NeumannBC should construct successfully";
}

// Test 4: Multiple BC tags register without error
TEST_F(BoundaryConditionManagerTest, MultipleBCRegistration)
{
  BoundaryConditionManager mgr;

  // Register 6 BCs (one per box face)
  for (int i = 1; i <= 6; i++)
  {
    mgr.addBoundaryCondition(i, std::make_unique<DirichletBC>());
  }

  // Verify no crash and PML still unset
  EXPECT_FALSE(mgr.hasPML());
}

// Test 5: Default BC can be set
TEST_F(BoundaryConditionManagerTest, SetDefaultBC)
{
  BoundaryConditionManager mgr;
  mgr.setDefaultBC(std::make_unique<NeumannBC>());
  // No crash means default BC was accepted
  EXPECT_FALSE(mgr.hasPML());
}
