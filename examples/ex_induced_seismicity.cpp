/*
 * Advanced Example: Induced Seismicity
 * 
 * Demonstrates:
 * - Fully coupled fluid flow and geomechanics
 * - Poroelastic effects
 * - Pre-existing fault modeling
 * - Rate-and-state friction
 * - Seismic moment calculation
 * - Stress monitoring
 * 
 * Scenario: High-pressure injection near a critically stressed fault
 */

#include "Simulator.hpp"
#include "FractureModel.hpp"
#include <iostream>

static char help[] = "Example: Induced seismicity from fluid injection\n\n";

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, help);
    
    MPI_Comm comm = PETSC_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    if (rank == 0) {
        std::cout << "================================================\n";
        std::cout << "  Induced Seismicity Example\n";
        std::cout << "================================================\n\n";
        std::cout << "Simulating:\n";
        std::cout << "- High-pressure water injection\n";
        std::cout << "- Poroelastic stress changes\n";
        std::cout << "- Fault reactivation\n";
        std::cout << "- Seismic event detection\n\n";
    }
    
    // Create simulator
    ResSim::Simulator sim(comm);
    
    // Configuration for fully coupled simulation
    ResSim::SimulationConfig config;
    config.start_time = 0.0;
    config.end_time = 365.25 * 86400.0;  // 1 year
    config.dt_initial = 3600.0;           // 1 hour
    config.dt_min = 60.0;                 // 1 minute (for rapid events)
    config.dt_max = 86400.0;              // 1 day
    config.enable_adaptive_timestepping = true;
    
    // Enable coupled physics
    config.fluid_model = ResSim::FluidModelType::SINGLE_COMPONENT;
    config.solid_model = ResSim::SolidModelType::POROELASTIC;
    config.enable_geomechanics = true;
    config.enable_faults = true;
    config.enable_tidal_forces = true;  // Tidal triggering
    
    // Tighter tolerances for coupled problem
    config.rtol = 1e-8;
    config.atol = 1e-10;
    
    sim.initialize(config);
    
    // Grid: Fine mesh near fault
    ResSim::GridConfig grid;
    grid.nx = 50;
    grid.ny = 50;
    grid.nz = 30;
    grid.Lx = 5000.0;  // 5 km
    grid.Ly = 5000.0;
    grid.Lz = 3000.0;  // 3 km depth
    
    sim.setupDM();
    sim.setupFields();
    sim.setupPhysics();
    
    // Enable coupling between fluid flow and geomechanics
    sim.enableCoupling(ResSim::PhysicsType::FLUID_FLOW, 
                      ResSim::PhysicsType::GEOMECHANICS);
    
    // Add fault model
    if (rank == 0) {
        std::cout << "Setting up fault model...\n";
    }
    
    ResSim::FaultModel fault;
    
    // Fault geometry: strike-slip fault
    std::vector<double> strike = {1.0, 0.0, 0.0};  // E-W strike
    std::vector<double> dip = {0.0, 0.0, 1.0};     // Vertical
    fault.setFaultPlane(strike, dip, 2000.0, 1500.0);
    
    // Friction properties
    fault.setFrictionCoefficient(0.6, 0.4);  // static, dynamic
    fault.setCohesion(1.0e6);  // 1 MPa
    
    // Enable rate-and-state friction for realistic slip behavior
    fault.enableRateStateFriction(true);
    fault.setRateStateParameters(0.010, 0.015, 0.001);  // a, b, Dc
    // Note: b > a means velocity weakening (unstable)
    
    // Initial stress state (critically stressed)
    // In-situ stress: sigma_v = rho*g*z
    // sigma_H = 1.2 * sigma_v (tectonic compression)
    // sigma_h = 0.8 * sigma_v
    
    // Add high-rate injection well
    if (rank == 0) {
        std::cout << "Setting up injection well near fault...\n";
    }
    
    sim.addWell("INJ_DEEP", ResSim::WellType::INJECTOR);
    // Position: 500m from fault at 2km depth
    sim.setWellControl("INJ_DEEP", 0.1);  // 100 L/s = 0.1 m³/s
    
    // Observation wells for pressure monitoring
    sim.addWell("OBS1", ResSim::WellType::OBSERVATION);  // Near fault
    sim.addWell("OBS2", ResSim::WellType::OBSERVATION);  // Far field
    
    // Set initial conditions
    sim.setInitialConditions();
    
    // Setup solvers with options for coupled problem
    sim.setupTimeStepper();
    sim.setupSolvers();
    
    if (rank == 0) {
        std::cout << "\nStarting coupled simulation...\n";
        std::cout << "Monitoring for:\n";
        std::cout << "- Pressure buildup\n";
        std::cout << "- Stress changes\n";
        std::cout << "- Fault slip events\n";
        std::cout << "- Seismic moments\n\n";
    }
    
    // Run simulation
    PetscErrorCode ierr = sim.run(); CHKERRQ(ierr);
    
    // Post-processing
    if (rank == 0) {
        std::cout << "\n================================================\n";
        std::cout << "  Simulation Complete - Analyzing Results\n";
        std::cout << "================================================\n\n";
    }
    
    sim.writeSummary();
    sim.computePerformanceMetrics();
    sim.generatePlots();
    
    if (rank == 0) {
        std::cout << "Output files generated:\n";
        std::cout << "- Pressure evolution: output/induced_seismicity_pressure_*.vtu\n";
        std::cout << "- Stress field: output/induced_seismicity_stress_*.vtu\n";
        std::cout << "- Displacement: output/induced_seismicity_displacement_*.vtu\n";
        std::cout << "- Seismic catalog: output/induced_seismicity_events.txt\n";
        std::cout << "- Pressure-time plot: output/induced_seismicity_pressure_history.png\n";
        std::cout << "- Moment magnitude: output/induced_seismicity_magnitude_time.png\n\n";
        
        std::cout << "Visualization:\n";
        std::cout << "  paraview output/induced_seismicity_*.vtu\n\n";
        
        std::cout << "Key Questions Addressed:\n";
        std::cout << "1. How far does pressure propagate from injection point?\n";
        std::cout << "2. What is the stress perturbation on the fault?\n";
        std::cout << "3. When does fault slip initiate?\n";
        std::cout << "4. What is the magnitude of induced events?\n";
        std::cout << "5. What is the role of tidal stresses in triggering?\n\n";
    }
    
    PetscFinalize();
    return 0;
}
