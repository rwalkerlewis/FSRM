/*
 * Advanced Example: Induced Seismicity with Black Oil Formulation
 * 
 * Demonstrates:
 * - Fully coupled poroelastic simulation
 * - Black oil formulation (oil/water/gas phases)
 * - Pre-existing fault modeling with rate-and-state friction
 * - High-pressure wastewater injection
 * - Stress evolution and fault reactivation
 * - Seismic event detection and cataloging
 * 
 * Scenario: Deep wastewater injection near a critically stressed basement fault
 * This example simulates induced seismicity from Class II disposal wells,
 * a major environmental concern in areas with intensive oil/gas operations.
 * 
 * This example is configuration-driven. All physics parameters are defined
 * in the config file (default: config/induced_seismicity.config).
 * 
 * The fluid models (BlackOilFluid), material models (PoroelasticMaterial),
 * and fault models (RateStateFriction) are from the generic FSRM libraries.
 * 
 * Usage:
 *   mpirun -np 4 ./ex_induced_seismicity -c config/induced_seismicity.config
 */

#include "Simulator.hpp"
#include <iostream>

static char help[] = "Advanced Example: Induced seismicity simulation\n"
                     "Usage: mpirun -np N ./ex_induced_seismicity -c <config_file>\n\n";

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, help);
    
    MPI_Comm comm = PETSC_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    // Get config file from command line
    char config_file[PETSC_MAX_PATH_LEN] = "config/induced_seismicity.config";
    PetscOptionsGetString(nullptr, nullptr, "-c", config_file,
                         sizeof(config_file), nullptr);
    
    if (rank == 0) {
        std::cout << "\n";
        std::cout << "════════════════════════════════════════════════════════════════\n";
        std::cout << "  INDUCED SEISMICITY SIMULATION\n";
        std::cout << "  Configuration-Driven with Generic Physics Libraries\n";
        std::cout << "════════════════════════════════════════════════════════════════\n\n";
        std::cout << "Config file: " << config_file << "\n\n";
        std::cout << "This simulation uses:\n";
        std::cout << "  • BlackOilFluid model from FluidModel library\n";
        std::cout << "  • PoroelasticMaterial from MaterialModel library\n";
        std::cout << "  • RateStateFriction from FaultModel library\n\n";
        std::cout << "Physical Setup (from config):\n";
        std::cout << "  • Deep basement fault with rate-state friction\n";
        std::cout << "  • Class II wastewater injection well\n";
        std::cout << "  • Three-phase flow (oil/water/gas)\n";
        std::cout << "  • Fully coupled poroelasticity\n\n";
    }
    
    // Create simulator
    FSRM::Simulator sim(comm);
    
    // Initialize entirely from config file
    // The ConfigReader automatically:
    // - Creates BlackOilFluid from [FLUID] section
    // - Creates PoroelasticMaterial from [ROCK] section
    // - Creates FaultNetwork with RateStateFriction from [FAULT*] sections
    PetscErrorCode ierr = sim.initializeFromConfigFile(config_file);
    CHKERRQ(ierr);
    
    // Setup
    ierr = sim.setupDM(); CHKERRQ(ierr);
    ierr = sim.setupFields(); CHKERRQ(ierr);
    ierr = sim.setupPhysics(); CHKERRQ(ierr);
    ierr = sim.setMaterialProperties(); CHKERRQ(ierr);
    ierr = sim.setInitialConditions(); CHKERRQ(ierr);
    ierr = sim.setupTimeStepper(); CHKERRQ(ierr);
    ierr = sim.setupSolvers(); CHKERRQ(ierr);
    
    if (rank == 0) {
        std::cout << "Solver initialized. Running simulation...\n\n";
    }
    
    // Run simulation
    ierr = sim.run(); CHKERRQ(ierr);
    
    // Post-processing
    sim.writeSummary();
    sim.computePerformanceMetrics();
    sim.generatePlots();
    
    if (rank == 0) {
        std::cout << "\n════════════════════════════════════════════════════════════════\n";
        std::cout << "  SIMULATION COMPLETE\n";
        std::cout << "════════════════════════════════════════════════════════════════\n\n";
        std::cout << "Output files:\n";
        std::cout << "  - Pressure/stress fields: output/induced_seismicity/*.vtu\n";
        std::cout << "  - Seismic catalog: output/seismic_catalog.csv\n\n";
        std::cout << "To modify simulation parameters:\n";
        std::cout << "  - Edit " << config_file << "\n";
        std::cout << "  - Re-run without recompilation\n\n";
        std::cout << "Key configuration sections:\n";
        std::cout << "  [FLUID] - Black oil PVT properties\n";
        std::cout << "  [ROCK]  - Poroelastic material properties\n";
        std::cout << "  [FAULT1]- Fault geometry and friction law\n";
        std::cout << "  [SEISMICITY] - Event detection settings\n\n";
    }
    
    PetscFinalize();
    return 0;
}
