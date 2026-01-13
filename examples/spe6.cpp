/*
 * SPE6: Dual-Porosity / Dual-Permeability Problem
 * 
 * Configuration-driven simulation of the SPE6 benchmark.
 * Tests naturally fractured reservoir simulation.
 * 
 * All physics, materials, and dual-porosity parameters are in the config file.
 * See config/spe6_benchmark.config for the complete parameter set.
 * 
 * Key Features (from config):
 * - Dual-porosity / dual-permeability formulation
 * - Matrix-fracture transfer functions
 * - Shape factor sensitivity
 * 
 * Reference: Firoozabadi, A. and Thomas, L.K., "Sixth SPE Comparative
 *            Solution Project: Dual-Porosity Simulators", JPT (1990)
 * 
 * Usage:
 *   mpirun -np 4 ./spe6 -c config/spe6_benchmark.config
 */

#include "Simulator.hpp"
#include <iostream>

static char help[] = "SPE6: Dual-Porosity Benchmark\n"
                     "Usage: mpirun -np N ./spe6 -c <config_file>\n\n";

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, help);
    
    MPI_Comm comm = PETSC_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    char config_file[PETSC_MAX_PATH_LEN] = "config/spe6_benchmark.config";
    PetscOptionsGetString(nullptr, nullptr, "-c", config_file,
                         sizeof(config_file), nullptr);
    
    if (rank == 0) {
        std::cout << "========================================\n";
        std::cout << "  SPE6 Benchmark - Dual-Porosity\n";
        std::cout << "  Configuration-Driven\n";
        std::cout << "========================================\n\n";
        std::cout << "Config file: " << config_file << "\n\n";
        std::cout << "All parameters loaded from configuration file.\n\n";
    }
    
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
        std::cout << "Running SPE6 simulation...\n\n";
    }
    
    ierr = sim.run(); CHKERRQ(ierr);
    
    sim.writeSummary();
    sim.computePerformanceMetrics();
    sim.generatePlots();
    
    if (rank == 0) {
        std::cout << "\n========================================\n";
        std::cout << "  SPE6 Completed Successfully\n";
        std::cout << "========================================\n\n";
        std::cout << "Results available in output directory.\n";
    }
    
    PetscFinalize();
    return 0;
}
