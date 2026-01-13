/*
 * SCEC LOH.2: Layer Over Halfspace - Heterogeneous Wave Test
 * 
 * Configuration-driven simulation of SCEC LOH.2 benchmark problem.
 * Tests wave propagation in a more complex layered medium.
 * 
 * All physics, materials, sources, and receivers are in the config file.
 * See config/scec_loh2.config for the complete parameter set.
 * 
 * Key Features (from config):
 * - Multiple velocity layers
 * - Double-couple point source at depth
 * - Surface wave generation
 * - Comparison with semi-analytical solutions
 * 
 * Reference: SCEC Community Velocity Model (CVM) Tests
 *            https://strike.scec.org/scecpedia/LOH
 * 
 * Usage:
 *   mpirun -np 8 ./scec_loh2 -c config/scec_loh2.config
 */

#include "Simulator.hpp"
#include <iostream>

static char help[] = "SCEC LOH.2: Layer Over Halfspace - Heterogeneous\n"
                     "Usage: mpirun -np N ./scec_loh2 -c <config_file>\n\n";

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, help);
    
    MPI_Comm comm = PETSC_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    char config_file[PETSC_MAX_PATH_LEN] = "config/scec_loh2.config";
    PetscOptionsGetString(nullptr, nullptr, "-c", config_file,
                         sizeof(config_file), nullptr);
    
    if (rank == 0) {
        std::cout << "========================================\n";
        std::cout << "  SCEC LOH.2 Benchmark\n";
        std::cout << "  Layer Over Halfspace - Heterogeneous\n";
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
        std::cout << "Running SCEC LOH.2 simulation...\n\n";
    }
    
    ierr = sim.run(); CHKERRQ(ierr);
    
    sim.writeSummary();
    sim.computePerformanceMetrics();
    sim.generatePlots();
    
    if (rank == 0) {
        std::cout << "\n========================================\n";
        std::cout << "  SCEC LOH.2 Completed Successfully\n";
        std::cout << "========================================\n\n";
        std::cout << "Results available in output directory.\n";
        std::cout << "Compare with SCEC reference solutions:\n";
        std::cout << "  https://strike.scec.org/scecpedia/LOH.2\n\n";
    }
    
    PetscFinalize();
    return 0;
}
