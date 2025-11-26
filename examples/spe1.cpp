/*
 * SPE1 Benchmark: First SPE Comparative Solution Project
 * 
 * Classic 3-phase black oil problem with gas dissolution
 * 10x10x3 Cartesian grid
 * Single production well, gas injection well
 * 
 * Reference: Odeh, A.S., "Comparison of Solutions to a Three-Dimensional
 *            Black-Oil Reservoir Simulation Problem", JPT (1981)
 */

#include "Simulator.hpp"
#include <iostream>

static char help[] = "SPE1: First SPE Comparative Solution Project\n\n";

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, help);
    
    MPI_Comm comm = PETSC_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    if (rank == 0) {
        std::cout << "========================================\n";
        std::cout << "  SPE1 Benchmark\n";
        std::cout << "========================================\n\n";
        std::cout << "Problem Description:\n";
        std::cout << "- 3-phase black oil model\n";
        std::cout << "- 10x10x3 grid cells\n";
        std::cout << "- Gas dissolution\n";
        std::cout << "- Producer + Gas injector\n\n";
    }
    
    // Create simulator
    FSRM::Simulator sim(comm);
    
    // Load SPE1 data file
    PetscErrorCode ierr = sim.loadEclipseInput("SPE1.DATA");
    if (ierr) {
        if (rank == 0) {
            std::cout << "Note: SPE1.DATA not found. Using programmatic setup.\n\n";
        }
        
        // Manual configuration
        FSRM::SimulationConfig config;
        config.start_time = 0.0;
        config.end_time = 365.25 * 10.0 * 86400.0;  // 10 years
        config.dt_initial = 30.0 * 86400.0;          // 30 days
        config.fluid_model = FSRM::FluidModelType::BLACK_OIL;
        config.output_frequency = 10;
        
        sim.initialize(config);
        
        // Grid: 10x10x3
        FSRM::GridConfig grid;
        grid.nx = 10;
        grid.ny = 10;
        grid.nz = 3;
        grid.Lx = 304.8;  // 1000 ft in meters
        grid.Ly = 304.8;
        grid.Lz = 6.096;  // 20 ft
        
        sim.setupDM();
        sim.setupFields();
        sim.setupPhysics();
        
        // Producer at (10,10,3)
        sim.addWell("PRODUCER", FSRM::WellType::PRODUCER);
        sim.setWellControl("PRODUCER", 6.178e-4);  // 20,000 STB/day to m³/s
        
        // Injector at (1,1,1)
        sim.addWell("INJECTOR", FSRM::WellType::INJECTOR);
        sim.setWellControl("INJECTOR", 3.089e-4);  // 10,000 MSCF/day
        
        sim.setInitialConditions();
        sim.setupTimeStepper();
        sim.setupSolvers();
    }
    
    // Run simulation
    if (rank == 0) {
        std::cout << "Running SPE1 simulation...\n\n";
    }
    
    ierr = sim.run(); CHKERRQ(ierr);
    
    // Post-processing
    if (rank == 0) {
        std::cout << "\nGenerating output and plots...\n";
    }
    
    sim.writeSummary();
    sim.computePerformanceMetrics();
    sim.generatePlots();
    
    if (rank == 0) {
        std::cout << "\n========================================\n";
        std::cout << "  SPE1 Completed Successfully\n";
        std::cout << "========================================\n";
        std::cout << "\nResults available in: output/spe1_*\n";
        std::cout << "\nKey outputs:\n";
        std::cout << "- Pressure maps: output/spe1_pressure_*.vtu\n";
        std::cout << "- Production history: output/spe1_production.png\n";
        std::cout << "- Summary data: output/spe1_SUMMARY.txt\n\n";
    }
    
    PetscFinalize();
    return 0;
}
