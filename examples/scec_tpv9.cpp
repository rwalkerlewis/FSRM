/*
 * SCEC TPV9: 3D Strike-Slip with Depth-Dependent Stress
 * 
 * Configuration-driven benchmark simulation.
 * All parameters specified in config/scec_tpv9.config
 * 
 * Reference: https://strike.scec.org/cvws/tpv9docs.html
 */

#include "Simulator.hpp"
#include <iostream>

static char help[] = "SCEC TPV9: Depth-Dependent Stress\n";

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, help);
    
    MPI_Comm comm = PETSC_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    char config_file[PETSC_MAX_PATH_LEN] = "config/scec_tpv9.config";
    PetscOptionsGetString(nullptr, nullptr, "-c", config_file,
                         sizeof(config_file), nullptr);
    
    if (rank == 0) {
        std::cout << "SCEC TPV9 Benchmark - Depth-Dependent Stress\n";
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
