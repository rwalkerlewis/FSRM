/*
 * Example: Wave Propagation with Dynamic Effects
 * 
 * Configuration-driven simulation of elastic and poroelastic wave propagation.
 * Demonstrates multiple scenarios that can be selected via config file.
 * 
 * Features:
 * - Elastodynamic wave propagation
 * - Poroelastodynamic waves (Biot theory)
 * - Static-to-dynamic triggering
 * - Permeability changes from transient waves
 * 
 * Usage:
 *   mpirun -np 4 ./ex_wave_propagation -c config/elastodynamic_waves.config
 *   mpirun -np 4 ./ex_wave_propagation -c config/poroelastodynamic_waves.config
 *   mpirun -np 4 ./ex_wave_propagation -c config/static_triggered_seismicity.config
 */

#include "Simulator.hpp"
#include <iostream>

static char help[] = "Example: Wave propagation simulations\n"
                     "Usage: mpirun -np N ./ex_wave_propagation -c <config_file>\n\n"
                     "Available configs:\n"
                     "  config/elastodynamic_waves.config\n"
                     "  config/poroelastodynamic_waves.config\n"
                     "  config/static_triggered_seismicity.config\n"
                     "  config/wave_permeability_enhancement.config\n\n";

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, help);
    
    MPI_Comm comm = PETSC_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    // Get config file from command line
    char config_file[PETSC_MAX_PATH_LEN] = "config/elastodynamic_waves.config";
    PetscOptionsGetString(nullptr, nullptr, "-c", config_file,
                         sizeof(config_file), nullptr);
    
    if (rank == 0) {
        std::cout << "======================================\n";
        std::cout << "  Wave Propagation Simulation\n";
        std::cout << "======================================\n\n";
        std::cout << "Config file: " << config_file << "\n\n";
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
        std::cout << "Starting wave propagation simulation...\n\n";
    }
    
    ierr = sim.run(); CHKERRQ(ierr);
    
    // Post-processing
    sim.writeSummary();
    sim.computePerformanceMetrics();
    sim.generatePlots();
    
    if (rank == 0) {
        std::cout << "\n======================================\n";
        std::cout << "  Simulation Complete\n";
        std::cout << "======================================\n\n";
        std::cout << "To run different wave scenarios:\n";
        std::cout << "  • Elastodynamics: -c config/elastodynamic_waves.config\n";
        std::cout << "  • Poroelastodynamics: -c config/poroelastodynamic_waves.config\n";
        std::cout << "  • Triggered seismicity: -c config/static_triggered_seismicity.config\n";
        std::cout << "  • Permeability enhancement: -c config/wave_permeability_enhancement.config\n\n";
    }
    
    PetscFinalize();
    return 0;
}
