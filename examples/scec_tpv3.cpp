/*
 * SCEC TPV3: 2D In-Plane Strike-Slip Benchmark
 * 
 * Configuration-driven simulation of SCEC TPV3 benchmark problem.
 * Tests 2D spontaneous rupture on a vertical strike-slip fault.
 * 
 * All physics, materials, and fault parameters are in the config file.
 * See config/scec_tpv3.config for the complete parameter set.
 * 
 * Key Features (from config):
 * - 2D plane strain model
 * - Vertical strike-slip fault
 * - Linear slip-weakening friction
 * - Nucleation zone with forced rupture
 * 
 * Reference: SCEC Dynamic Rupture Code Verification Exercise
 *            Harris et al., Seism. Res. Lett. (2009)
 *            https://strike.scec.org/cvws/
 * 
 * Usage:
 *   mpirun -np 4 ./scec_tpv3 -c config/scec_tpv3.config
 */

#include "Simulator.hpp"
#include <iostream>

static char help[] = "SCEC TPV3: 2D Strike-Slip Dynamic Rupture Benchmark\n"
                     "Usage: mpirun -np N ./scec_tpv3 -c <config_file>\n\n";

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, help);
    
    MPI_Comm comm = PETSC_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    char config_file[PETSC_MAX_PATH_LEN] = "config/scec_tpv3.config";
    PetscOptionsGetString(nullptr, nullptr, "-c", config_file,
                         sizeof(config_file), nullptr);
    
    if (rank == 0) {
        std::cout << "========================================\n";
        std::cout << "  SCEC TPV3 Benchmark\n";
        std::cout << "  2D In-Plane Strike-Slip Rupture\n";
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
        std::cout << "Running SCEC TPV3 simulation...\n\n";
    }
    
    ierr = sim.run(); CHKERRQ(ierr);
    
    sim.writeSummary();
    sim.computePerformanceMetrics();
    sim.generatePlots();
    
    if (rank == 0) {
        std::cout << "\n========================================\n";
        std::cout << "  SCEC TPV3 Completed Successfully\n";
        std::cout << "========================================\n\n";
        std::cout << "Results available in output directory.\n";
        std::cout << "Compare with SCEC reference solutions:\n";
        std::cout << "  https://strike.scec.org/cvws/tpv3docs.html\n\n";
    }
    
    PetscFinalize();
    return 0;
}
