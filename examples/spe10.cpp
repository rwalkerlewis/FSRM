/*
 * SPE10: Comparative Solution Project Model 2
 * 
 * Configuration-driven simulation of the SPE10 benchmark.
 * Tests scalability with large-scale heterogeneous permeability.
 * 
 * Key Features:
 * - 60x220x85 grid (1.1 million cells)
 * - Extreme permeability variations (8 orders of magnitude)
 * - Two-phase flow (water-oil)
 * - Tests upscaling and numerical methods
 * 
 * Reference: Christie, M.A. and Blunt, M.J., "Tenth SPE Comparative
 *            Solution Project: A Comparison of Upscaling Techniques",
 *            SPE 66599 (2001)
 * 
 * Usage:
 *   mpirun -np 32 ./spe10 -c config/spe10_benchmark.config
 */

#include "Simulator.hpp"
#include <iostream>

static char help[] = "SPE10: Large-Scale Heterogeneous Benchmark\n"
                     "Usage: mpirun -np N ./spe10 -c <config_file>\n\n";

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, help);
    
    MPI_Comm comm = PETSC_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    // Get config file from command line
    char config_file[PETSC_MAX_PATH_LEN] = "config/spe10_benchmark.config";
    PetscOptionsGetString(nullptr, nullptr, "-c", config_file,
                         sizeof(config_file), nullptr);
    
    if (rank == 0) {
        std::cout << "========================================\n";
        std::cout << "  SPE10 Benchmark - Large Heterogeneous\n";
        std::cout << "  Configuration-Driven\n";
        std::cout << "========================================\n\n";
        std::cout << "Config file: " << config_file << "\n\n";
        std::cout << "Problem Description:\n";
        std::cout << "  • Two-phase flow (water-oil)\n";
        std::cout << "  • 60x220x85 grid (1.1M cells)\n";
        std::cout << "  • 8 orders of magnitude permeability\n";
        std::cout << "  • Scalability test case\n\n";
    }
    
    // Create and run simulation from config
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
        std::cout << "Running SPE10 simulation...\n\n";
        std::cout << "Note: This is a large-scale benchmark.\n";
        std::cout << "      Recommended: >= 32 MPI processes\n\n";
    }
    
    ierr = sim.run(); CHKERRQ(ierr);
    
    // Post-processing
    sim.writeSummary();
    sim.computePerformanceMetrics();
    sim.generatePlots();
    
    if (rank == 0) {
        std::cout << "\n========================================\n";
        std::cout << "  SPE10 Completed Successfully\n";
        std::cout << "========================================\n\n";
        std::cout << "Results available in: output/spe10/\n";
    }
    
    PetscFinalize();
    return 0;
}
