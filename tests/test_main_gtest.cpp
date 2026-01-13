/**
 * @file test_main_gtest.cpp
 * @brief Main entry point for GTest test suite
 */

#include <gtest/gtest.h>
#include <mpi.h>
#include <petsc.h>

int main(int argc, char **argv) {
    // Initialize MPI
    MPI_Init(&argc, &argv);
    
    // Initialize PETSc (required for simulation tests)
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    
    // Initialize GTest
    ::testing::InitGoogleTest(&argc, argv);
    
    // Get MPI rank for output control
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // Only print from rank 0
    if (rank != 0) {
        ::testing::TestEventListeners& listeners =
            ::testing::UnitTest::GetInstance()->listeners();
        delete listeners.Release(listeners.default_result_printer());
    }
    
    // Run all tests
    int result = RUN_ALL_TESTS();
    
    // Finalize PETSc and MPI
    PetscFinalize();
    MPI_Finalize();
    
    return result;
}
