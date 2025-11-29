/*
 * SCEC TPV5: Strike-Slip Benchmark with Heterogeneous Initial Stress
 * 
 * Configuration-driven simulation of SCEC TPV5 benchmark problem.
 * Tests spontaneous dynamic rupture propagation on a vertical strike-slip fault.
 * 
 * All physics, materials, and boundary conditions are specified in the config file.
 * See config/scec_tpv5.config for the complete parameter set.
 * 
 * Key Features (from config):
 * - 3D elastic half-space
 * - Vertical planar fault (strike-slip)
 * - Heterogeneous initial stress with asperity
 * - Spontaneous rupture nucleation
 * - Linear slip-weakening friction
 * 
 * Reference: SCEC Dynamic Rupture Code Verification Exercise
 *            Harris et al., Seism. Res. Lett. (2009)
 *            https://strike.scec.org/cvws/
 * 
 * Usage:
 *   mpirun -np 8 ./scec_tpv5 -c config/scec_tpv5.config
 */

#include "Simulator.hpp"
#include <iostream>

static char help[] = "SCEC TPV5: Strike-Slip Dynamic Rupture Benchmark\n"
                     "Usage: mpirun -np N ./scec_tpv5 -c <config_file>\n\n";

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, help);
    
    MPI_Comm comm = PETSC_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    // Get config file from command line
    char config_file[PETSC_MAX_PATH_LEN] = "config/scec_tpv5.config";
    PetscOptionsGetString(nullptr, nullptr, "-c", config_file,
                         sizeof(config_file), nullptr);
    
    if (rank == 0) {
        std::cout << "========================================\n";
        std::cout << "  SCEC TPV5 Benchmark\n";
        std::cout << "  Dynamic Rupture on Strike-Slip Fault\n";
        std::cout << "========================================\n\n";
        std::cout << "Config file: " << config_file << "\n\n";
        std::cout << "All parameters loaded from configuration file.\n";
        std::cout << "Edit the config file to modify simulation settings.\n\n";
    }
    
    // Create and run simulation entirely from config
    FSRM::Simulator sim(comm);
    
    PetscErrorCode ierr = sim.initializeFromConfigFile(config_file);
    CHKERRQ(ierr);
    
    // Setup (fault, initial stress, etc. configured from file)
    ierr = sim.setupDM(); CHKERRQ(ierr);
    ierr = sim.setupFields(); CHKERRQ(ierr);
    ierr = sim.setupPhysics(); CHKERRQ(ierr);
    ierr = sim.setMaterialProperties(); CHKERRQ(ierr);
    ierr = sim.setInitialConditions(); CHKERRQ(ierr);
    ierr = sim.setupTimeStepper(); CHKERRQ(ierr);
    ierr = sim.setupSolvers(); CHKERRQ(ierr);
    
    if (rank == 0) {
        std::cout << "Running SCEC TPV5 simulation...\n\n";
    }
    
    ierr = sim.run(); CHKERRQ(ierr);
    
    // Post-processing
    sim.writeSummary();
    sim.computePerformanceMetrics();
    sim.generatePlots();
    
    if (rank == 0) {
        std::cout << "\n========================================\n";
        std::cout << "  SCEC TPV5 Completed Successfully\n";
        std::cout << "========================================\n\n";
        std::cout << "Results available in output directory.\n";
        std::cout << "Compare with SCEC reference solutions:\n";
        std::cout << "  https://strike.scec.org/cvws/tpv5docs.html\n\n";
    }
    
    PetscFinalize();
    return 0;
}
