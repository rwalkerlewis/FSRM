/*
 * SCEC TPV10: Strike-Slip Benchmark with Branching Fault
 * 
 * Configuration-driven simulation of SCEC TPV10 benchmark problem.
 * Tests rupture propagation onto a branching fault.
 * 
 * All physics, materials, and fault geometry are specified in the config file.
 * See config/scec_tpv10.config for the complete parameter set.
 * 
 * Reference: SCEC Dynamic Rupture Code Verification Exercise
 *            Harris et al., Seism. Res. Lett. (2009)
 *            https://strike.scec.org/cvws/
 * 
 * Usage:
 *   mpirun -np 16 ./scec_tpv10 -c config/scec_tpv10.config
 */

#include "Simulator.hpp"
#include <iostream>

static char help[] = "SCEC TPV10: Dynamic Rupture with Branching Fault\n"
                     "Usage: mpirun -np N ./scec_tpv10 -c <config_file>\n\n";

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, help);
    
    MPI_Comm comm = PETSC_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    char config_file[PETSC_MAX_PATH_LEN] = "config/scec_tpv10.config";
    PetscOptionsGetString(nullptr, nullptr, "-c", config_file,
                         sizeof(config_file), nullptr);
    
    if (rank == 0) {
        std::cout << "========================================\n";
        std::cout << "  SCEC TPV10 Benchmark\n";
        std::cout << "  Branching Fault Dynamic Rupture\n";
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
        std::cout << "Running SCEC TPV10 simulation...\n\n";
    }
    
    ierr = sim.run(); CHKERRQ(ierr);
    
    sim.writeSummary();
    sim.computePerformanceMetrics();
    sim.generatePlots();
    
    if (rank == 0) {
        std::cout << "\n========================================\n";
        std::cout << "  SCEC TPV10 Completed Successfully\n";
        std::cout << "========================================\n\n";
        std::cout << "Results available in output directory.\n";
    }
    
    PetscFinalize();
    return 0;
}
