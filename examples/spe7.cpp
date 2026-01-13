/*
 * SPE7: Horizontal Well Problem
 * 
 * Configuration-driven simulation of the SPE7 benchmark.
 * Tests horizontal well modeling in heterogeneous reservoirs.
 * 
 * All physics, materials, and well parameters are in the config file.
 * See config/spe7_benchmark.config for the complete parameter set.
 * 
 * Key Features (from config):
 * - Horizontal well trajectory
 * - Near-wellbore inflow modeling
 * - Heterogeneous permeability
 * 
 * Reference: Nghiem, L.X. et al., "Seventh SPE Comparative Solution Project:
 *            Modeling of Horizontal Wells in Reservoir Simulation", 
 *            SPE 21221 (1991)
 * 
 * Usage:
 *   mpirun -np 8 ./spe7 -c config/spe7_benchmark.config
 */

#include "Simulator.hpp"
#include <iostream>

static char help[] = "SPE7: Horizontal Well Benchmark\n"
                     "Usage: mpirun -np N ./spe7 -c <config_file>\n\n";

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, help);
    
    MPI_Comm comm = PETSC_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    char config_file[PETSC_MAX_PATH_LEN] = "config/spe7_benchmark.config";
    PetscOptionsGetString(nullptr, nullptr, "-c", config_file,
                         sizeof(config_file), nullptr);
    
    if (rank == 0) {
        std::cout << "========================================\n";
        std::cout << "  SPE7 Benchmark - Horizontal Well\n";
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
        std::cout << "Running SPE7 simulation...\n\n";
    }
    
    ierr = sim.run(); CHKERRQ(ierr);
    
    sim.writeSummary();
    sim.computePerformanceMetrics();
    sim.generatePlots();
    
    if (rank == 0) {
        std::cout << "\n========================================\n";
        std::cout << "  SPE7 Completed Successfully\n";
        std::cout << "========================================\n\n";
        std::cout << "Results available in output directory.\n";
    }
    
    PetscFinalize();
    return 0;
}
