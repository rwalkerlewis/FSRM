/*
 * SPE9: Model of a North Sea Reservoir
 * 
 * Configuration-driven simulation of the SPE9 benchmark.
 * Tests black-oil simulation with complex geology.
 * 
 * Key Features:
 * - 24x25x15 grid (9,000 cells)
 * - Black-oil model with solution gas
 * - Multiple producing and injecting wells
 * - Heterogeneous permeability
 * - Gas coning effects
 * 
 * Reference: Killough, J.E., "Ninth SPE Comparative Solution Project:
 *            A Reexamination of Black-Oil Simulation", SPE 29110 (1995)
 * 
 * Usage:
 *   mpirun -np 8 ./spe9 -c config/spe9_benchmark.config
 */

#include "Simulator.hpp"
#include <iostream>

static char help[] = "SPE9: North Sea Reservoir Benchmark\n"
                     "Usage: mpirun -np N ./spe9 -c <config_file>\n\n";

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, help);
    
    MPI_Comm comm = PETSC_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    // Get config file from command line
    char config_file[PETSC_MAX_PATH_LEN] = "config/spe9_benchmark.config";
    PetscOptionsGetString(nullptr, nullptr, "-c", config_file,
                         sizeof(config_file), nullptr);
    
    if (rank == 0) {
        std::cout << "========================================\n";
        std::cout << "  SPE9 Benchmark - North Sea Model\n";
        std::cout << "  Configuration-Driven\n";
        std::cout << "========================================\n\n";
        std::cout << "Config file: " << config_file << "\n\n";
        std::cout << "Problem Description:\n";
        std::cout << "  • Black-oil model\n";
        std::cout << "  • 24x25x15 grid (9,000 cells)\n";
        std::cout << "  • Multiple wells\n";
        std::cout << "  • Heterogeneous reservoir\n\n";
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
        std::cout << "Running SPE9 simulation...\n\n";
    }
    
    ierr = sim.run(); CHKERRQ(ierr);
    
    // Post-processing
    sim.writeSummary();
    sim.computePerformanceMetrics();
    sim.generatePlots();
    
    if (rank == 0) {
        std::cout << "\n========================================\n";
        std::cout << "  SPE9 Completed Successfully\n";
        std::cout << "========================================\n\n";
        std::cout << "Results available in: output/spe9/\n";
    }
    
    PetscFinalize();
    return 0;
}
