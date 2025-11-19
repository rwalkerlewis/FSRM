/*
 * Example: Configuration File Driven Simulation
 * 
 * Demonstrates how to run simulations entirely from config files
 * without any hard-coded parameters in the source code.
 */

#include "Simulator.hpp"
#include <iostream>

static char help[] = "Example: Configuration file driven simulation\n\n";

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, help);
    
    MPI_Comm comm = PETSC_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    // Get config file from command line
    char config_file[PETSC_MAX_PATH_LEN] = "config/default.config";
    PetscOptionsGetString(nullptr, nullptr, "-c", config_file,
                         sizeof(config_file), nullptr);
    
    if (rank == 0) {
        std::cout << "================================================\n";
        std::cout << "  Configuration-Driven Simulation Example\n";
        std::cout << "================================================\n\n";
        std::cout << "Config file: " << config_file << "\n\n";
        std::cout << "This example shows how ALL simulation parameters\n";
        std::cout << "can be specified in an editable config file:\n";
        std::cout << "  - Rock properties\n";
        std::cout << "  - Fluid properties\n";
        std::cout << "  - Well locations and controls\n";
        std::cout << "  - Fracture properties\n";
        std::cout << "  - Fault properties\n";
        std::cout << "  - Particle/proppant properties\n";
        std::cout << "  - Boundary conditions\n";
        std::cout << "  - Initial conditions\n";
        std::cout << "  - All physics settings\n\n";
        std::cout << "No recompilation needed to change parameters!\n\n";
    }
    
    // Create simulator
    ResSim::Simulator sim(comm);
    
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
    
    // Run
    if (rank == 0) {
        std::cout << "Starting simulation...\n\n";
    }
    
    ierr = sim.run(); CHKERRQ(ierr);
    
    // Output
    sim.writeSummary();
    sim.computePerformanceMetrics();
    sim.generatePlots();
    
    if (rank == 0) {
        std::cout << "\n================================================\n";
        std::cout << "  Simulation Complete\n";
        std::cout << "================================================\n\n";
        std::cout << "To run with different parameters:\n";
        std::cout << "1. Copy and edit a config file:\n";
        std::cout << "   cp config/default.config my_custom.config\n";
        std::cout << "   nano my_custom.config\n\n";
        std::cout << "2. Run with your config:\n";
        std::cout << "   mpirun -np 4 ./ex_config_driven -c my_custom.config\n\n";
        std::cout << "Available example configs:\n";
        std::cout << "  - config/default.config\n";
        std::cout << "  - config/shale_reservoir.config\n";
        std::cout << "  - config/geothermal.config\n";
        std::cout << "  - config/induced_seismicity.config\n";
        std::cout << "  - config/co2_storage.config\n";
        std::cout << "  - config/hydraulic_fracturing.config\n";
        std::cout << "  - config/complete_template.config (all options)\n\n";
    }
    
    PetscFinalize();
    return 0;
}
