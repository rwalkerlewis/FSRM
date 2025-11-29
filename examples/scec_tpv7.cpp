/*
 * SCEC TPV7: 3D Strike-Slip with Velocity-Strengthening Barrier
 * 
 * Configuration-driven benchmark simulation.
 * All parameters specified in config/scec_tpv7.config
 * 
 * Reference: https://strike.scec.org/cvws/tpv7docs.html
 */

#include "Simulator.hpp"
#include <iostream>

static char help[] = "SCEC TPV7: Strike-Slip with VS Barrier\n";

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, help);
    
    MPI_Comm comm = PETSC_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    char config_file[PETSC_MAX_PATH_LEN] = "config/scec_tpv7.config";
    PetscOptionsGetString(nullptr, nullptr, "-c", config_file,
                         sizeof(config_file), nullptr);
    
    if (rank == 0) {
        std::cout << "SCEC TPV7 Benchmark - VS Barrier\n";
        std::cout << "Config: " << config_file << "\n\n";
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
    ierr = sim.run(); CHKERRQ(ierr);
    sim.computePerformanceMetrics();
    
    PetscFinalize();
    return 0;
}
