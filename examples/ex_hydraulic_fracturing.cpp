/*
 * Advanced Example: Hydraulic Fracturing
 * 
 * Demonstrates:
 * - Hydraulic fracture propagation (P3D model)
 * - Proppant transport and placement
 * - Leak-off to reservoir
 * - Fracture closure post-shut-in
 * - Production from fractured well
 */

#include "Simulator.hpp"
#include "FractureModel.hpp"
#include <iostream>

static char help[] = "Example: Hydraulic fracturing simulation\n\n";

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, help);
    
    MPI_Comm comm = PETSC_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    if (rank == 0) {
        std::cout << "================================================\n";
        std::cout << "  Hydraulic Fracturing Example\n";
        std::cout << "================================================\n\n";
    }
    
    // Create simulator
    ResSim::Simulator sim(comm);
    
    // Configuration
    ResSim::SimulationConfig config;
    config.start_time = 0.0;
    config.end_time = 7200.0;  // 2 hours of pumping + shut-in
    config.dt_initial = 10.0;   // 10 seconds
    config.dt_min = 1.0;
    config.dt_max = 60.0;
    
    config.fluid_model = ResSim::FluidModelType::SINGLE_COMPONENT;
    config.solid_model = ResSim::SolidModelType::ELASTIC;
    config.enable_geomechanics = true;
    config.enable_fractures = true;
    config.enable_particle_transport = true;
    
    sim.initialize(config);
    
    // Grid
    ResSim::GridConfig grid;
    grid.nx = 40;
    grid.ny = 40;
    grid.nz = 20;
    grid.Lx = 400.0;  // 400m
    grid.Ly = 400.0;
    grid.Lz = 100.0;
    
    sim.setupDM();
    sim.setupFields();
    sim.setupPhysics();
    
    // Enable coupling
    sim.enableCoupling(ResSim::PhysicsType::FLUID_FLOW,
                      ResSim::PhysicsType::GEOMECHANICS);
    sim.enableCoupling(ResSim::PhysicsType::FLUID_FLOW,
                      ResSim::PhysicsType::FRACTURE_PROPAGATION);
    sim.enableCoupling(ResSim::PhysicsType::FRACTURE_PROPAGATION,
                      ResSim::PhysicsType::PARTICLE_TRANSPORT);
    
    // Create hydraulic fracture
    if (rank == 0) {
        std::cout << "Initializing hydraulic fracture...\n";
    }
    
    // Initial fracture at wellbore
    std::vector<double> frac_coords = {
        200.0, 200.0, 50.0,  // Center position
        0.0, 1.0, 0.0        // Normal vector (E-W fracture)
    };
    
    sim.addFracture(ResSim::FractureType::INDUCED_HYDRAULIC, frac_coords);
    
    // Add horizontal well (injection point)
    sim.addWell("FRAC_WELL", ResSim::WellType::INJECTOR);
    
    // Fracturing schedule:
    // Stage 1: Pad (0-30 min) - low viscosity, no proppant
    // Stage 2: Proppant ramp (30-90 min) - increasing proppant concentration
    // Stage 3: Flush (90-120 min) - clear fluid
    // Stage 4: Shut-in and flowback
    
    if (rank == 0) {
        std::cout << "\nFracturing Schedule:\n";
        std::cout << "  0-30 min:   Pad stage (clean fluid)\n";
        std::cout << "  30-90 min:  Proppant stages (0-8 ppg)\n";
        std::cout << "  90-120 min: Flush\n";
        std::cout << "  >120 min:   Shut-in and monitoring\n\n";
    }
    
    // Set injection rate: 0.05 m³/s (about 20 bpm)
    sim.setWellControl("FRAC_WELL", 0.05);
    
    // Set initial conditions
    sim.setInitialConditions();
    sim.setupTimeStepper();
    sim.setupSolvers();
    
    // Run simulation
    if (rank == 0) {
        std::cout << "Starting hydraulic fracturing simulation...\n\n";
    }
    
    PetscErrorCode ierr = sim.run(); CHKERRQ(ierr);
    
    // Post-processing
    if (rank == 0) {
        std::cout << "\n================================================\n";
        std::cout << "  Fracturing Complete - Analyzing Results\n";
        std::cout << "================================================\n\n";
    }
    
    sim.writeSummary();
    sim.computePerformanceMetrics();
    sim.generatePlots();
    
    if (rank == 0) {
        std::cout << "Output files:\n";
        std::cout << "- Fracture geometry: output/hf_fracture_*.vtu\n";
        std::cout << "- Proppant distribution: output/hf_proppant_*.vtu\n";
        std::cout << "- Pressure field: output/hf_pressure_*.vtu\n";
        std::cout << "- Stress field: output/hf_stress_*.vtu\n";
        std::cout << "- Treatment plot: output/hf_treatment.png\n";
        std::cout << "- Net pressure: output/hf_net_pressure.png\n\n";
        
        std::cout << "Key Metrics:\n";
        std::cout << "- Final fracture length\n";
        std::cout << "- Final fracture height\n";
        std::cout << "- Average fracture width\n";
        std::cout << "- Proppant-packed length\n";
        std::cout << "- Estimated conductivity\n";
        std::cout << "- Leak-off volume\n\n";
    }
    
    PetscFinalize();
    return 0;
}
