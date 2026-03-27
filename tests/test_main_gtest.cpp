/**
 * @file test_main_gtest.cpp
 * @brief Main entry point for GTest test suite
 *
 * PetscInitialize handles MPI initialization internally.
 * No separate MPI_Init/MPI_Finalize needed.
 */

#include <gtest/gtest.h>
#include <petscsys.h>

int main(int argc, char **argv) {
    // Initialize PETSc (also initializes MPI)
    PetscInitialize(&argc, &argv, nullptr, nullptr);

    // Initialize GTest
    ::testing::InitGoogleTest(&argc, argv);

    // Get MPI rank for output control
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // Only print from rank 0
    if (rank != 0) {
        ::testing::TestEventListeners& listeners =
            ::testing::UnitTest::GetInstance()->listeners();
        delete listeners.Release(listeners.default_result_printer());
    }

    // Run all tests
    int result = RUN_ALL_TESTS();

    // Finalize PETSc (also finalizes MPI)
    PetscFinalize();

    return result;
}
