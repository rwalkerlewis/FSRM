/*
 * SPE3: Gas Cycling Project
 * 
 * Configuration-driven simulation of the SPE3 benchmark.
 * Tests compositional modeling with gas cycling.
 * 
 * Key Features:
 * - 9x9x4 grid (324 cells)
 * - 4-component compositional model
 * - Gas injection with cycling
 * - Pressure maintenance operations
 * 
 * Reference: Kenyon, D.E. and Behie, G.A., "Third SPE Comparative 
 *            Solution Project: Gas Cycling of Retrograde Condensate 
 *            Reservoirs", JPT (1987)
 * 
 * Usage:
 *   mpirun -np 4 ./spe3 -c config/spe3_benchmark.config
 */

#include "Simulator.hpp"
#include <iostream>

static char help[] = "SPE3: Gas Cycling Benchmark\n"
                     "Usage: mpirun -np N ./spe3 -c <config_file>\n\n";

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, help);
    
    MPI_Comm comm = PETSC_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    // Get config file from command line
    char config_file[PETSC_MAX_PATH_LEN] = "config/spe3_benchmark.config";
    PetscOptionsGetString(nullptr, nullptr, "-c", config_file,
                         sizeof(config_file), nullptr);
    
    if (rank == 0) {
        std::cout << "========================================\n";
        std::cout << "  SPE3 Benchmark - Gas Cycling\n";
        std::cout << "  Configuration-Driven\n";
        std::cout << "========================================\n\n";
        std::cout << "Config file: " << config_file << "\n\n";
        std::cout << "Problem Description:\n";
        std::cout << "  • 4-component compositional model\n";
        std::cout << "  • 9x9x4 grid cells\n";
        std::cout << "  • Gas cycling operations\n";
        std::cout << "  • Pressure maintenance\n\n";
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
        std::cout << "Running SPE3 simulation...\n\n";
    }
    
    ierr = sim.run(); CHKERRQ(ierr);
    
    // Post-processing
    sim.writeSummary();
    sim.computePerformanceMetrics();
    sim.generatePlots();
    
    if (rank == 0) {
        std::cout << "\n========================================\n";
        std::cout << "  SPE3 Completed Successfully\n";
        std::cout << "========================================\n\n";
        std::cout << "Results available in: output/spe3/\n";
    }
    
    PetscFinalize();
    return 0;
}
