/*
 * SPE11: CO2 Storage Benchmark
 * 
 * Configuration-driven simulation of the SPE11 benchmark.
 * Tests CO2 geological storage simulation capabilities.
 * 
 * All physics, materials, and CO2 properties are in the config file.
 * See config/spe11_benchmark.config for the complete parameter set.
 * 
 * Key Features (from config):
 * - CO2-brine two-phase flow
 * - Dissolution and trapping mechanisms
 * - Buoyancy-driven migration
 * - Long-term storage security
 * 
 * SPE11 includes three test cases:
 * - Case A: 2D vertical cross-section (small)
 * - Case B: 2D vertical cross-section (large)
 * - Case C: 3D full-field model
 * 
 * Reference: Nordbotten, J.M. et al., "The 11th Society of Petroleum
 *            Engineers Comparative Solution Project: Problem Definition",
 *            SPE Journal (2024)
 *            https://spe11-group.github.io/
 * 
 * Usage:
 *   mpirun -np 8 ./spe11 -c config/spe11_benchmark.config
 *   mpirun -np 8 ./spe11 -c config/spe11a_benchmark.config  # Case A
 *   mpirun -np 8 ./spe11 -c config/spe11b_benchmark.config  # Case B
 *   mpirun -np 32 ./spe11 -c config/spe11c_benchmark.config # Case C
 */

#include "Simulator.hpp"
#include <iostream>

static char help[] = "SPE11: CO2 Storage Benchmark\n"
                     "Usage: mpirun -np N ./spe11 -c <config_file>\n\n";

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, help);
    
    MPI_Comm comm = PETSC_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    char config_file[PETSC_MAX_PATH_LEN] = "config/spe11_benchmark.config";
    PetscOptionsGetString(nullptr, nullptr, "-c", config_file,
                         sizeof(config_file), nullptr);
    
    if (rank == 0) {
        std::cout << "========================================\n";
        std::cout << "  SPE11 Benchmark - CO2 Storage\n";
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
        std::cout << "Running SPE11 simulation...\n\n";
    }
    
    ierr = sim.run(); CHKERRQ(ierr);
    
    sim.writeSummary();
    sim.computePerformanceMetrics();
    sim.generatePlots();
    
    if (rank == 0) {
        std::cout << "\n========================================\n";
        std::cout << "  SPE11 Completed Successfully\n";
        std::cout << "========================================\n\n";
        std::cout << "Results available in output directory.\n";
        std::cout << "Compare with reference solutions at:\n";
        std::cout << "  https://spe11-group.github.io/\n\n";
    }
    
    PetscFinalize();
    return 0;
}
