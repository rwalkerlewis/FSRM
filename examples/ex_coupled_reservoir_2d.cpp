/*
 * Example: 2D Coupled Reservoir Simulation
 * 
 * Configuration-driven simulation of coupled flow and geomechanics.
 * All physics parameters are specified in the config file.
 * 
 * Demonstrates:
 * - Single-phase compressible flow
 * - Linear poroelasticity (Biot theory)
 * - Subsidence from pressure depletion
 * - Permeability changes from compaction
 * 
 * Usage:
 *   mpirun -np 4 ./ex_coupled_reservoir_2d -c config/coupled_reservoir_2d.config
 */

#include "Simulator.hpp"
#include <iostream>

static char help[] = "Example: 2D coupled flow-geomechanics simulation\n"
                     "Usage: mpirun -np N ./ex_coupled_reservoir_2d -c <config_file>\n\n";

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, help);
    
    MPI_Comm comm = PETSC_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    // Get config file from command line
    char config_file[PETSC_MAX_PATH_LEN] = "config/coupled_reservoir_2d.config";
    PetscOptionsGetString(nullptr, nullptr, "-c", config_file,
                         sizeof(config_file), nullptr);
    
    if (rank == 0) {
        std::cout << "================================================\n";
        std::cout << "  2D Coupled Reservoir Simulation\n";
        std::cout << "  Flow + Geomechanics\n";
        std::cout << "================================================\n\n";
        std::cout << "Config file: " << config_file << "\n\n";
        std::cout << "This simulation uses:\n";
        std::cout << "  • Single-phase compressible flow\n";
        std::cout << "  • Biot poroelasticity\n";
        std::cout << "  • Stress-dependent permeability\n";
        std::cout << "  • Multiple wells (injector + producers)\n\n";
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
        std::cout << "Starting coupled simulation...\n\n";
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
        std::cout << "Output files in: output/coupled_2d/\n";
        std::cout << "  • Pressure distribution (*.vtu)\n";
        std::cout << "  • Displacement field\n";
        std::cout << "  • Permeability changes\n\n";
        std::cout << "To modify parameters:\n";
        std::cout << "  Edit " << config_file << "\n";
        std::cout << "  Key sections: [ROCK], [FLUID], [WELL*]\n\n";
    }
    
    PetscFinalize();
    return 0;
}
