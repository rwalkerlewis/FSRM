/*
 * SCEC TPV24/25: Japan Benchmark with Rate-State Friction
 * 
 * Configuration-driven simulation of SCEC TPV24/25 benchmark problems.
 * Tests Japanese benchmark problem with rate-and-state friction.
 * 
 * TPV24: Homogeneous initial stress
 * TPV25: Heterogeneous initial stress (use scec_tpv25.config)
 * 
 * All physics, materials, and fault parameters are in the config file.
 * See config/scec_tpv24.config for the complete parameter set.
 * 
 * Key Features (from config):
 * - 3D elastic half-space
 * - Vertical strike-slip fault
 * - Rate-and-state friction with flash heating
 * - Strong velocity weakening
 * 
 * Reference: SCEC-Japan Collaborative Benchmark
 *            https://strike.scec.org/cvws/
 * 
 * Usage:
 *   mpirun -np 16 ./scec_tpv24 -c config/scec_tpv24.config
 *   mpirun -np 16 ./scec_tpv24 -c config/scec_tpv25.config
 */

#include "Simulator.hpp"
#include <iostream>

static char help[] = "SCEC TPV24/25: Japan Benchmark with Rate-State\n"
                     "Usage: mpirun -np N ./scec_tpv24 -c <config_file>\n\n";

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, help);
    
    MPI_Comm comm = PETSC_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    char config_file[PETSC_MAX_PATH_LEN] = "config/scec_tpv24.config";
    PetscOptionsGetString(nullptr, nullptr, "-c", config_file,
                         sizeof(config_file), nullptr);
    
    if (rank == 0) {
        std::cout << "========================================\n";
        std::cout << "  SCEC TPV24/25 Benchmark\n";
        std::cout << "  Japan Rate-State Problem\n";
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
        std::cout << "Running SCEC TPV24/25 simulation...\n\n";
    }
    
    ierr = sim.run(); CHKERRQ(ierr);
    
    sim.writeSummary();
    sim.computePerformanceMetrics();
    sim.generatePlots();
    
    if (rank == 0) {
        std::cout << "\n========================================\n";
        std::cout << "  SCEC TPV24/25 Completed Successfully\n";
        std::cout << "========================================\n\n";
        std::cout << "Results available in output directory.\n";
        std::cout << "Compare with SCEC reference solutions:\n";
        std::cout << "  https://strike.scec.org/cvws/tpv24_25docs.html\n\n";
    }
    
    PetscFinalize();
    return 0;
}
