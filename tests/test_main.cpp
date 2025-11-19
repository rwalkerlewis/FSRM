#include "Testing.hpp"
#include <mpi.h>
#include <iostream>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    
    if (rank == 0) {
        std::cout << "\n";
        std::cout << "========================================\n";
        std::cout << "  ReservoirSim Test Suite\n";
        std::cout << "========================================\n";
        std::cout << "\n";
    }
    
    // Create test suite
    ResSim::Testing::TestSuite suite("ReservoirSim");
    
    // Add unit tests
    suite.addTest(std::make_shared<ResSim::Testing::EclipseIOTest>());
    suite.addTest(std::make_shared<ResSim::Testing::PhysicsKernelTest>(ResSim::PhysicsType::FLUID_FLOW));
    suite.addTest(std::make_shared<ResSim::Testing::WellModelTest>());
    
    // Run all tests
    auto results = suite.runAll();
    
    if (rank == 0) {
        std::cout << "\n";
        std::cout << "========================================\n";
        std::cout << "  Test Results\n";
        std::cout << "========================================\n";
        std::cout << "Total: " << results.total_tests << "\n";
        std::cout << "Passed: " << results.passed << "\n";
        std::cout << "Failed: " << results.failed << "\n";
        std::cout << "Time: " << results.total_time << " s\n";
        std::cout << "\n";
        
        // Generate report
        suite.generateReport(results, "test_report.txt");
    }
    
    PetscFinalize();
    MPI_Finalize();
    
    return (results.failed == 0) ? 0 : 1;
}
