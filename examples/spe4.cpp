/*
 * SPE4: Equation of State Testing Problem
 * 
 * Configuration-driven simulation of the SPE4 benchmark.
 * Tests cubic equation of state (EOS) for compositional modeling.
 * 
 * All physics, materials, and EOS parameters are in the config file.
 * See config/spe4_benchmark.config for the complete parameter set.
 * 
 * Key Features (from config):
 * - Peng-Robinson EOS
 * - Multi-component flash calculations
 * - Phase behavior testing
 * 
 * Reference: Whitson, C.H. and Brulé, M.R., "Phase Behavior", SPE (2000)
 * 
 * Usage:
 *   mpirun -np 4 ./spe4 -c config/spe4_benchmark.config
 */

#include "Simulator.hpp"
#include <iostream>

static char help[] = "SPE4: Equation of State Testing Benchmark\n"
                     "Usage: mpirun -np N ./spe4 -c <config_file>\n\n";

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, help);
    
    MPI_Comm comm = PETSC_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    char config_file[PETSC_MAX_PATH_LEN] = "config/spe4_benchmark.config";
    PetscOptionsGetString(nullptr, nullptr, "-c", config_file,
                         sizeof(config_file), nullptr);
    
    if (rank == 0) {
        std::cout << "========================================\n";
        std::cout << "  SPE4 Benchmark - EOS Testing\n";
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
        std::cout << "Running SPE4 simulation...\n\n";
    }
    
    ierr = sim.run(); CHKERRQ(ierr);
    
    sim.writeSummary();
    sim.computePerformanceMetrics();
    sim.generatePlots();
    
    if (rank == 0) {
        std::cout << "\n========================================\n";
        std::cout << "  SPE4 Completed Successfully\n";
        std::cout << "========================================\n\n";
        std::cout << "Results available in output directory.\n";
    }
    
    PetscFinalize();
    return 0;
}
