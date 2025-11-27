/*
 * SPE1 Benchmark: First SPE Comparative Solution Project
 * 
 * Configuration-driven simulation of the classic SPE1 benchmark.
 * All parameters specified in config file for reproducibility.
 * 
 * Classic 3-phase black oil problem with gas dissolution:
 * - 10x10x3 Cartesian grid
 * - Single production well
 * - Gas injection well
 * 
 * Reference: Odeh, A.S., "Comparison of Solutions to a Three-Dimensional
 *            Black-Oil Reservoir Simulation Problem", JPT (1981)
 * 
 * Usage:
 *   mpirun -np 4 ./spe1 -c config/spe1_benchmark.config
 */

#include "Simulator.hpp"
#include <iostream>

static char help[] = "SPE1: First SPE Comparative Solution Project\n"
                     "Usage: mpirun -np N ./spe1 -c <config_file>\n\n";

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, help);
    
    MPI_Comm comm = PETSC_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    // Get config file from command line
    char config_file[PETSC_MAX_PATH_LEN] = "config/spe1_benchmark.config";
    PetscOptionsGetString(nullptr, nullptr, "-c", config_file,
                         sizeof(config_file), nullptr);
    
    if (rank == 0) {
        std::cout << "========================================\n";
        std::cout << "  SPE1 Benchmark\n";
        std::cout << "  Configuration-Driven\n";
        std::cout << "========================================\n\n";
        std::cout << "Config file: " << config_file << "\n\n";
        std::cout << "Problem Description:\n";
        std::cout << "  • 3-phase black oil model\n";
        std::cout << "  • 10x10x3 grid cells\n";
        std::cout << "  • Gas dissolution\n";
        std::cout << "  • Producer + Gas injector\n\n";
    }
    
    // Create and run simulation from config
    FSRM::Simulator sim(comm);
    
    PetscErrorCode ierr = sim.initializeFromConfigFile(config_file);
    CHKERRQ(ierr);
    
    ierr = sim.setupDM(); CHKERRQ(ierr);
    ierr = sim.setupFields(); CHKERRQ(ierr);
    ierr = sim.setupPhysics(); CHKERRQ(ierr);
    ierr = sim.setMaterialProperties(); CHKERRQ(ierr);
    ierr = sim.setInitialConditions(); CHKERRQ(ierr);
    ierr = sim.setupTimeStepper(); CHKERRQ(ierr);
    ierr = sim.setupSolvers(); CHKERRQ(ierr);
    
    if (rank == 0) {
        std::cout << "Running SPE1 simulation...\n\n";
    }
    
    ierr = sim.run(); CHKERRQ(ierr);
    
    // Post-processing
    sim.writeSummary();
    sim.computePerformanceMetrics();
    sim.generatePlots();
    
    if (rank == 0) {
        std::cout << "\n========================================\n";
        std::cout << "  SPE1 Completed Successfully\n";
        std::cout << "========================================\n\n";
        std::cout << "Results available in: output/spe1/\n";
        std::cout << "  • Pressure maps: spe1_pressure_*.vtu\n";
        std::cout << "  • Saturation maps: spe1_saturation_*.vtu\n";
        std::cout << "  • Summary data: spe1_SUMMARY.txt\n\n";
    }
    
    PetscFinalize();
    return 0;
}
