/*
 * SCEC TPV10: Strike-Slip Benchmark with Branching Fault
 * 
 * Configuration-driven simulation of SCEC TPV10 benchmark problem.
 * Tests rupture propagation onto a branching fault.
 * 
 * Key Features:
 * - Main strike-slip fault with 30° branch
 * - Tests whether rupture propagates onto branch
 * - Realistic stress conditions
 * - Rate-and-state friction on branch
 * 
 * Reference: SCEC Dynamic Rupture Code Verification Exercise
 *            Harris et al., Seism. Res. Lett. (2009)
 *            https://strike.scec.org/cvws/
 * 
 * Domain:
 * - 48 km × 48 km × 24 km
 * - Main fault: 30 km × 15 km
 * - Branch fault: 15 km × 15 km at 30° angle
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
        std::cout << "Problem Description:\n";
        std::cout << "  • Main fault + 30° branch fault\n";
        std::cout << "  • Tests rupture branch jumping\n";
        std::cout << "  • Rate-and-state friction\n";
        std::cout << "  • Complex fault geometry\n\n";
    }
    
    FSRM::Simulator sim(comm);
    
    PetscErrorCode ierr = sim.initializeFromConfigFile(config_file);
    CHKERRQ(ierr);
    
    ierr = sim.setupDM(); CHKERRQ(ierr);
    ierr = sim.setupFields(); CHKERRQ(ierr);
    ierr = sim.setupPhysics(); CHKERRQ(ierr);
    ierr = sim.setMaterialProperties(); CHKERRQ(ierr);
    ierr = sim.setInitialStress(); CHKERRQ(ierr);
    
    // Setup main fault and branch
    ierr = sim.setupFault(); CHKERRQ(ierr);
    ierr = sim.setupBranchFault(); CHKERRQ(ierr);
    
    ierr = sim.setInitialConditions(); CHKERRQ(ierr);
    ierr = sim.setupTimeStepper(); CHKERRQ(ierr);
    ierr = sim.setupSolvers(); CHKERRQ(ierr);
    
    if (rank == 0) {
        std::cout << "Running SCEC TPV10 simulation...\n\n";
    }
    
    ierr = sim.run(); CHKERRQ(ierr);
    
    sim.writeSummary();
    sim.computeRuptureMetrics();
    sim.computeBranchingMetrics();  // Branch activation time, branch slip
    sim.generatePlots();
    
    if (rank == 0) {
        std::cout << "\n========================================\n";
        std::cout << "  SCEC TPV10 Completed Successfully\n";
        std::cout << "========================================\n\n";
        std::cout << "Results available in: output/scec_tpv10/\n";
    }
    
    PetscFinalize();
    return 0;
}
