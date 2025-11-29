/*
 * SCEC TPV21: Shallow Dip-Slip - Rate-State Friction
 * 
 * Configuration-driven benchmark simulation.
 * All parameters specified in config/scec_tpv21.config
 * 
 * Reference: https://strike.scec.org/cvws/tpv20_21docs.html
 */

#include "Simulator.hpp"
#include <iostream>

static char help[] = "SCEC TPV21: Shallow Dip-Slip Rate-State\n";

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, help);
    
    MPI_Comm comm = PETSC_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    char config_file[PETSC_MAX_PATH_LEN] = "config/scec_tpv21.config";
    PetscOptionsGetString(nullptr, nullptr, "-c", config_file,
                         sizeof(config_file), nullptr);
    
    if (rank == 0) {
        std::cout << "SCEC TPV21 Benchmark - Shallow Dip-Slip Rate-State\n";
        std::cout << "Config: " << config_file << "\n\n";
    }
    
    FSRM::Simulator sim(comm);
    
    sim.initializeFromConfigFile(config_file);
    sim.setupDM();
    sim.setupFields();
    sim.setupPhysics();
    sim.setMaterialProperties();
    sim.setInitialConditions();
    sim.setupTimeStepper();
    sim.setupSolvers();
    sim.run();
    sim.computePerformanceMetrics();
    
    PetscFinalize();
    return 0;
}
