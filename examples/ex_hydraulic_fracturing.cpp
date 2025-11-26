/*
 * Advanced Example: Hydraulic Fracturing
 * 
 * Demonstrates:
 * - Hydraulic fracture propagation (LEFM-based model)
 * - Proppant transport and placement
 * - Leak-off to reservoir
 * - Fracture closure post-shut-in
 * - Production from fractured well
 * 
 * This example is configuration-driven. All parameters are specified
 * in the config file (default: config/hydraulic_fracturing.config).
 * 
 * The fracture mechanics and material properties are from the generic
 * FSRM libraries (MaterialModel, FractureModel).
 * 
 * Usage:
 *   mpirun -np 4 ./ex_hydraulic_fracturing -c config/hydraulic_fracturing.config
 */

#include "Simulator.hpp"
#include <iostream>

static char help[] = "Example: Hydraulic fracturing simulation\n"
                     "Usage: mpirun -np N ./ex_hydraulic_fracturing -c <config_file>\n\n";

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, help);
    
    MPI_Comm comm = PETSC_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    // Get config file from command line
    char config_file[PETSC_MAX_PATH_LEN] = "config/hydraulic_fracturing.config";
    PetscOptionsGetString(nullptr, nullptr, "-c", config_file,
                         sizeof(config_file), nullptr);
    
    if (rank == 0) {
        std::cout << "================================================\n";
        std::cout << "  Hydraulic Fracturing Example\n";
        std::cout << "  Configuration-Driven Simulation\n";
        std::cout << "================================================\n\n";
        std::cout << "Config file: " << config_file << "\n\n";
        std::cout << "This simulation uses:\n";
        std::cout << "  • LinearElasticMaterial from MaterialModel library\n";
        std::cout << "  • LEFM-based fracture propagation\n";
        std::cout << "  • Proppant transport model\n\n";
    }
    
    // Create simulator
    FSRM::Simulator sim(comm);
    
    // Initialize entirely from config file
    // The ConfigReader automatically:
    // - Creates material model from [ROCK] section
    // - Creates fracture model from [FRACTURE*] sections
    // - Sets up well from [WELL*] sections
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
        std::cout << "Starting hydraulic fracturing simulation...\n\n";
    }
    
    ierr = sim.run(); CHKERRQ(ierr);
    
    // Post-processing
    sim.writeSummary();
    sim.computePerformanceMetrics();
    sim.generatePlots();
    
    if (rank == 0) {
        std::cout << "\n================================================\n";
        std::cout << "  Fracturing Complete\n";
        std::cout << "================================================\n\n";
        std::cout << "Output files:\n";
        std::cout << "  - Fracture geometry: output/hydraulic_fracturing/*.vtu\n";
        std::cout << "  - Pressure field: output/hydraulic_fracturing/*.vtu\n\n";
        std::cout << "To modify parameters:\n";
        std::cout << "  - Edit " << config_file << "\n";
        std::cout << "  - Key sections: [ROCK], [FLUID], [FRACTURE1], [WELL1]\n\n";
    }
    
    PetscFinalize();
    return 0;
}
