/*
 * Example: Buckley-Leverett 2D Waterflooding
 * 
 * Configuration-driven simulation of immiscible two-phase displacement.
 * All physics parameters are specified in the config file.
 * 
 * Demonstrates:
 * - Two-phase immiscible flow (Buckley-Leverett theory)
 * - Fractional flow and shock front
 * - Saturation tracking
 * 
 * Usage:
 *   mpirun -np 4 ./ex_buckley_leverett_2d -c config/buckley_leverett_2d.config
 */

#include "Simulator.hpp"
#include <iostream>

static char help[] = "Example: Buckley-Leverett 2D waterflooding\n"
                     "Usage: mpirun -np N ./ex_buckley_leverett_2d -c <config_file>\n\n";

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, help);
    
    MPI_Comm comm = PETSC_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    // Get config file from command line
    char config_file[PETSC_MAX_PATH_LEN] = "config/buckley_leverett_2d.config";
    PetscOptionsGetString(nullptr, nullptr, "-c", config_file,
                         sizeof(config_file), nullptr);
    
    if (rank == 0) {
        std::cout << "================================================\n";
        std::cout << "  Buckley-Leverett 2D Waterflooding\n";
        std::cout << "  Configuration-Driven Simulation\n";
        std::cout << "================================================\n\n";
        std::cout << "Config file: " << config_file << "\n\n";
        std::cout << "This example demonstrates:\n";
        std::cout << "  • Two-phase immiscible displacement\n";
        std::cout << "  • Fractional flow theory\n";
        std::cout << "  • Shock front propagation\n";
        std::cout << "  • Saturation tracking\n\n";
    }
    
    // Create and run simulation entirely from config
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
        std::cout << "Starting Buckley-Leverett simulation...\n\n";
    }
    
    ierr = sim.run(); CHKERRQ(ierr);
    
    // Post-processing
    sim.writeSummary();
    sim.computePerformanceMetrics();
    sim.generatePlots();
    
    if (rank == 0) {
        std::cout << "\n================================================\n";
        std::cout << "  Simulation Complete\n";
        std::cout << "================================================\n\n";
        std::cout << "Output files in: output/buckley_leverett/\n";
        std::cout << "  • Saturation maps (*.vtu)\n";
        std::cout << "  • Pressure distribution\n\n";
        std::cout << "To modify parameters:\n";
        std::cout << "  Edit " << config_file << "\n";
        std::cout << "  Key sections: [FLUID], [ROCK], [WELL*]\n\n";
    }
    
    PetscFinalize();
    return 0;
}
