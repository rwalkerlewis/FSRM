/* 
 * Example 1: Single-Phase Flow
 * 
 * Demonstrates:
 * - Basic single-phase Darcy flow
 * - Structured grid setup
 * - Well injection/production
 * - Pressure transient response
 */

#include "Simulator.hpp"
#include <iostream>

static char help[] = "Example 1: Single-phase flow in homogeneous reservoir\n\n";

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, help);
    
    MPI_Comm comm = PETSC_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    if (rank == 0) {
        std::cout << "Running Example 1: Single-Phase Flow\n";
        std::cout << "====================================\n\n";
    }
    
    // Create simulator
    ResSim::Simulator sim(comm);
    
    // Configure simulation
    ResSim::SimulationConfig config;
    config.start_time = 0.0;
    config.end_time = 3600.0 * 24.0;  // 24 hours
    config.dt_initial = 3600.0;        // 1 hour
    config.output_frequency = 1;
    config.fluid_model = ResSim::FluidModelType::SINGLE_COMPONENT;
    config.enable_geomechanics = false;
    config.enable_thermal = false;
    
    // Grid configuration
    ResSim::GridConfig grid;
    grid.nx = 20;
    grid.ny = 20;
    grid.nz = 5;
    grid.Lx = 1000.0;  // 1 km
    grid.Ly = 1000.0;
    grid.Lz = 50.0;
    
    // Initialize
    sim.initialize(config);
    sim.setupDM();
    sim.setupFields();
    sim.setupPhysics();
    
    // Add production well at center
    sim.addWell("PROD1", ResSim::WellType::PRODUCER);
    sim.setWellControl("PROD1", 100.0);  // 100 m³/day
    
    // Add injection well at corner
    sim.addWell("INJ1", ResSim::WellType::INJECTOR);
    sim.setWellControl("INJ1", 150.0);  // 150 m³/day
    
    // Set initial conditions
    sim.setInitialConditions();
    
    // Setup solvers
    sim.setupTimeStepper();
    sim.setupSolvers();
    
    // Run simulation
    if (rank == 0) {
        std::cout << "Starting simulation...\n";
    }
    
    PetscErrorCode ierr = sim.run(); CHKERRQ(ierr);
    
    // Generate outputs
    sim.writeSummary();
    sim.computePerformanceMetrics();
    sim.generatePlots();
    
    if (rank == 0) {
        std::cout << "\nSimulation completed successfully!\n";
        std::cout << "Output files: output/ex01_*\n";
    }
    
    PetscFinalize();
    return 0;
}
