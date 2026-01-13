/*
 * Induced Seismicity with Implicit-Explicit Time Integration
 * 
 * Configuration-driven simulation demonstrating automatic transitions between
 * quasi-static (implicit) and dynamic (explicit) time integration for modeling
 * injection-induced seismicity.
 * 
 * All IMEX settings are configured in the config file's [IMEX] section.
 * See config/induced_seismicity_imex.config for the complete parameter set.
 * 
 * Usage:
 *   mpirun -np 4 ./induced_seismicity_imex -c config/induced_seismicity_imex.config
 */

#include "Simulator.hpp"
#include <iostream>

static char help[] = "Induced Seismicity with IMEX Time Integration\n"
                     "Usage: mpirun -np N ./induced_seismicity_imex -c <config_file>\n\n";

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, help);
    
    MPI_Comm comm = PETSC_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    char config_file[PETSC_MAX_PATH_LEN] = "config/induced_seismicity_imex.config";
    PetscOptionsGetString(nullptr, nullptr, "-c", config_file,
                         sizeof(config_file), nullptr);
    
    if (rank == 0) {
        std::cout << "============================================================\n";
        std::cout << "  Induced Seismicity with IMEX Time Integration\n";
        std::cout << "============================================================\n\n";
        std::cout << "Config file: " << config_file << "\n\n";
        std::cout << "All parameters loaded from configuration file.\n";
        std::cout << "IMEX transitions are handled automatically based on:\n";
        std::cout << "  - Coulomb failure criterion\n";
        std::cout << "  - Kinetic energy thresholds\n";
        std::cout << "  - Slip rate monitoring\n\n";
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
        std::cout << "Starting simulation...\n\n";
    }
    
    ierr = sim.run(); CHKERRQ(ierr);
    
    sim.writeSummary();
    sim.computePerformanceMetrics();
    sim.generatePlots();
    
    if (rank == 0) {
        std::cout << "\n============================================================\n";
        std::cout << "  Simulation Completed Successfully\n";
        std::cout << "============================================================\n\n";
        std::cout << "Results available in output directory.\n\n";
    }
    
    PetscFinalize();
    return 0;
}
