#include "Testing.hpp"
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    
    // Get timestamp for log file
    std::time_t now = std::time(nullptr);
    char timestamp[100];
    std::strftime(timestamp, sizeof(timestamp), "%Y-%m-%d_%H-%M-%S", std::localtime(&now));
    std::string log_filename = std::string("test_results_") + timestamp + ".log";
    
    std::ofstream logfile;
    if (rank == 0) {
        logfile.open(log_filename);
        
        auto print_header = [&](std::ostream& out) {
            out << "\n";
            out << "========================================\n";
            out << "  FSRM Test Suite\n";
            out << "  " << timestamp << "\n";
            out << "========================================\n";
            out << "\n";
        };
        
        print_header(std::cout);
        print_header(logfile);
    }
    
    // Create test suite
    FSRM::Testing::TestSuite suite("FSRM");
    
    // Add unit tests for all physics types
    if (rank == 0) {
        auto print_section = [&](std::ostream& out, const std::string& title) {
            out << "--- " << title << " ---\n";
        };
        
        print_section(std::cout, "Adding Physics Tests");
        print_section(logfile, "Adding Physics Tests");
    }
    
    // Test all physics types
    std::vector<FSRM::PhysicsType> physics_types = {
        FSRM::PhysicsType::FLUID_FLOW,
        FSRM::PhysicsType::GEOMECHANICS,
        FSRM::PhysicsType::THERMAL,
        FSRM::PhysicsType::ELASTODYNAMICS,
        FSRM::PhysicsType::POROELASTODYNAMICS
    };
    
    std::map<FSRM::PhysicsType, std::string> physics_names = {
        {FSRM::PhysicsType::FLUID_FLOW, "Fluid Flow"},
        {FSRM::PhysicsType::GEOMECHANICS, "Geomechanics"},
        {FSRM::PhysicsType::THERMAL, "Thermal"},
        {FSRM::PhysicsType::ELASTODYNAMICS, "Elastodynamics"},
        {FSRM::PhysicsType::POROELASTODYNAMICS, "Poroelastodynamics"}
    };
    
    for (auto physics : physics_types) {
        if (rank == 0) {
            std::cout << "  - " << physics_names[physics] << " kernel test\n";
            logfile << "  - " << physics_names[physics] << " kernel test\n";
        }
        suite.addTest(std::make_shared<FSRM::Testing::PhysicsKernelTest>(physics));
    }
    
    // Add other unit tests
    if (rank == 0) {
        std::cout << "  - EclipseIO test\n";
        std::cout << "  - WellModel test\n";
        std::cout << "  - MMS Convergence test\n";
        std::cout << "  - Spatial Convergence test\n";
        std::cout << "  - Temporal Convergence test\n";
        std::cout << "  - Detailed Physics test\n\n";
        logfile << "  - EclipseIO test\n";
        logfile << "  - WellModel test\n";
        logfile << "  - MMS Convergence test\n";
        logfile << "  - Spatial Convergence test\n";
        logfile << "  - Temporal Convergence test\n";
        logfile << "  - Detailed Physics test\n\n";
    }
    suite.addTest(std::make_shared<FSRM::Testing::EclipseIOTest>());
    suite.addTest(std::make_shared<FSRM::Testing::WellModelTest>());
    
    // Add MMS convergence tests (declared in test_mms_simple.cpp)
    extern std::shared_ptr<FSRM::Testing::UnitTest> createSimpleMMSTest();
    extern std::shared_ptr<FSRM::Testing::UnitTest> createSpatialConvergenceTest();
    extern std::shared_ptr<FSRM::Testing::UnitTest> createTemporalConvergenceTest();
    extern std::shared_ptr<FSRM::Testing::UnitTest> createDetailedPhysicsTest();
    
    suite.addTest(createSimpleMMSTest());
    suite.addTest(createSpatialConvergenceTest());
    suite.addTest(createTemporalConvergenceTest());
    suite.addTest(createDetailedPhysicsTest());
    
    // Run all tests
    if (rank == 0) {
        std::cout << "Running tests...\n\n";
        logfile << "Running tests...\n\n";
    }
    
    auto results = suite.runAll();
    
    if (rank == 0) {
        auto print_results = [&](std::ostream& out) {
            out << "\n";
            out << "========================================\n";
            out << "  Test Results Summary\n";
            out << "========================================\n";
            out << "Total Tests:  " << results.total_tests << "\n";
            out << "Passed:       " << results.passed << " (" 
                << std::fixed << std::setprecision(1)
                << (100.0 * results.passed / results.total_tests) << "%)\n";
            out << "Failed:       " << results.failed << "\n";
            out << "Total Time:   " << std::fixed << std::setprecision(3) 
                << results.total_time << " s\n";
            out << "\n";
            
            if (results.failed > 0) {
                out << "Failed Tests:\n";
                for (size_t i = 0; i < results.failed_test_names.size(); ++i) {
                    out << "  - " << results.failed_test_names[i] << "\n";
                    out << "    Error: " << results.error_messages[i] << "\n";
                }
                out << "\n";
            }
            
            out << "Detailed report saved to: test_report.txt\n";
            out << "Test log saved to: " << log_filename << "\n";
            out << "\n";
        };
        
        print_results(std::cout);
        print_results(logfile);
        
        // Generate detailed report
        suite.generateReport(results, "test_report.txt");
        
        logfile.close();
    }
    
    PetscFinalize();
    MPI_Finalize();
    
    return (results.failed == 0) ? 0 : 1;
}
