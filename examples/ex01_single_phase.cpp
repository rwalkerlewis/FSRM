/*
 * Example 1: Single-Phase Flow
 * 
 * Demonstrates:
 * - Basic single-phase Darcy flow
 * - Structured grid setup
 * - Well injection/production
 * - Pressure transient response
 * 
 * This example is configuration-driven. All parameters are specified
 * in the config file (default: config/single_phase.config).
 * 
 * Usage:
 *   mpirun -np 4 ./ex01_single_phase -c config/single_phase.config
 */

#include "Simulator.hpp"
#include <iostream>

static char help[] = "Example 1: Single-phase flow in homogeneous reservoir\n"
                     "Usage: mpirun -np N ./ex01_single_phase -c <config_file>\n\n";

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, help);
    
    MPI_Comm comm = PETSC_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    // Get config file from command line (default: config/single_phase.config)
    char config_file[PETSC_MAX_PATH_LEN] = "config/single_phase.config";
    PetscOptionsGetString(nullptr, nullptr, "-c", config_file,
                         sizeof(config_file), nullptr);
    
    if (rank == 0) {
        std::cout << "================================================\n";
        std::cout << "  Example 1: Single-Phase Flow\n";
        std::cout << "================================================\n\n";
        std::cout << "Config file: " << config_file << "\n\n";
        std::cout << "This example demonstrates:\n";
        std::cout << "  - Basic single-phase Darcy flow\n";
        std::cout << "  - Structured grid setup\n";
        std::cout << "  - Well injection/production\n";
        std::cout << "  - Pressure transient response\n\n";
    }
    
    // Create simulator
    FSRM::Simulator sim(comm);
    
    // Initialize entirely from config file
    PetscErrorCode ierr = sim.initializeFromConfigFile(config_file);
    CHKERRQ(ierr);
    
    // Setup
    ierr = sim.setupDM(); CHKERRQ(ierr);
    ierr = sim.setupFields(); CHKERRQ(ierr);
    ierr = sim.setupPhysics(); CHKERRQ(ierr);
    ierr = sim.setMaterialProperties(); CHKERRQ(ierr);
    ierr = sim.setInitialConditions(); CHKERRQ(ierr);
    ierr = sim.setupTimeStepper(); CHKERRQ(ierr);
    ierr = sim.setupSolvers(); CHKERRQ(ierr);
    
    // Run simulation
    if (rank == 0) {
        std::cout << "Starting simulation...\n\n";
    }
    
    ierr = sim.run(); CHKERRQ(ierr);
    
    // Generate outputs
    sim.writeSummary();
    sim.computePerformanceMetrics();
    sim.generatePlots();
    
    if (rank == 0) {
        std::cout << "\n================================================\n";
        std::cout << "  Simulation Complete\n";
        std::cout << "================================================\n\n";
        std::cout << "To modify parameters, edit the config file:\n";
        std::cout << "  " << config_file << "\n\n";
    }
    
    PetscFinalize();
    return 0;
}
