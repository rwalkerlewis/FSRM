/*
 * Example: Wave Propagation with Dynamic Permeability Changes
 * 
 * This example demonstrates:
 * 1. Elastodynamic wave propagation
 * 2. Poroelastodynamic waves (Biot theory)
 * 3. Static-to-dynamic triggering
 * 4. Permeability changes from transient waves
 */

#include "Simulator.hpp"
#include "PhysicsKernel.hpp"
#include "Visualization.hpp"
#include <iostream>
#include <cmath>

using namespace ResSim;

// Example 1: Pure elastodynamic waves
PetscErrorCode runElastodynamicWaves(MPI_Comm comm) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    PetscPrintf(comm, "\n=== Example 1: Elastodynamic Wave Propagation ===\n");
    
    // Create simulator
    Simulator sim(comm);
    
    // Configure simulation
    SimulationConfig config;
    config.start_time = 0.0;
    config.end_time = 2.0;        // 2 seconds
    config.dt_initial = 0.0001;   // 0.1 ms
    config.dt_min = 0.00001;
    config.dt_max = 0.001;
    config.max_timesteps = 20000;
    config.output_frequency = 100;
    config.output_file = "elastodynamic_waves";
    config.output_format = "VTK";
    
    // Enable dynamic mode
    config.use_dynamic_mode = true;
    config.enable_elastodynamics = true;
    config.enable_geomechanics = false;
    config.solid_model = SolidModelType::ELASTODYNAMIC;
    
    // Grid configuration
    GridConfig grid;
    grid.nx = 100;
    grid.ny = 100;
    grid.nz = 50;
    grid.Lx = 1000.0;
    grid.Ly = 1000.0;
    grid.Lz = 500.0;
    
    // Material properties (granite)
    MaterialProperties mat;
    mat.youngs_modulus = 50.0e9;
    mat.poisson_ratio = 0.25;
    mat.density = 2700.0;
    mat.p_wave_velocity = 5500.0;
    mat.s_wave_velocity = 3200.0;
    mat.quality_factor = 200.0;
    mat.damping_alpha = 0.001;
    mat.damping_beta = 0.0001;
    
    // Initialize simulation
    ierr = sim.initialize(config); CHKERRQ(ierr);
    
    PetscPrintf(comm, "P-wave velocity: %.1f m/s\n", mat.p_wave_velocity);
    PetscPrintf(comm, "S-wave velocity: %.1f m/s\n", mat.s_wave_velocity);
    PetscPrintf(comm, "Quality factor: %.1f\n", mat.quality_factor);
    
    // Run simulation
    ierr = sim.run(); CHKERRQ(ierr);
    
    PetscPrintf(comm, "Elastodynamic simulation completed.\n");
    
    PetscFunctionReturn(0);
}

// Example 2: Poroelastodynamic waves with permeability changes
PetscErrorCode runPoroelastodynamicWaves(MPI_Comm comm) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    PetscPrintf(comm, "\n=== Example 2: Poroelastodynamic Waves ===\n");
    
    Simulator sim(comm);
    
    // Configuration
    SimulationConfig config;
    config.start_time = 0.0;
    config.end_time = 5.0;
    config.dt_initial = 0.0001;
    config.dt_min = 0.00001;
    config.dt_max = 0.001;
    config.max_timesteps = 50000;
    config.output_frequency = 100;
    config.output_file = "poroelastodynamic_waves";
    config.output_format = "HDF5";
    
    // Enable poroelastodynamics
    config.use_dynamic_mode = true;
    config.enable_poroelastodynamics = true;
    config.solid_model = SolidModelType::POROELASTODYNAMIC;
    config.fluid_model = FluidModelType::SINGLE_COMPONENT;
    
    // Enable dynamic permeability changes
    config.enable_dynamic_permeability_change = true;
    config.permeability_sensitivity = 2.0;
    config.permeability_recovery_time = 100.0;
    
    // Grid
    GridConfig grid;
    grid.nx = 80;
    grid.ny = 80;
    grid.nz = 40;
    grid.Lx = 800.0;
    grid.Ly = 800.0;
    grid.Lz = 400.0;
    
    // Material properties (sandstone)
    MaterialProperties mat;
    mat.youngs_modulus = 15.0e9;
    mat.poisson_ratio = 0.2;
    mat.density = 2300.0;
    mat.porosity = 0.25;
    mat.permeability_x = 100.0;  // mD
    mat.permeability_y = 100.0;
    mat.permeability_z = 50.0;
    mat.biot_coefficient = 0.9;
    
    // Wave properties
    mat.p_wave_velocity = 3500.0;
    mat.s_wave_velocity = 2000.0;
    mat.quality_factor = 50.0;
    mat.damping_alpha = 0.01;
    mat.damping_beta = 0.001;
    
    // Permeability dynamics
    mat.permeability_strain_coeff = 1.0e-7;
    mat.permeability_stress_coeff = 5.0e-16;
    mat.min_permeability = 10.0;
    mat.max_permeability = 1000.0;
    
    // Fluid properties (water)
    FluidProperties fluid;
    fluid.density = 1000.0;
    fluid.viscosity = 0.001;
    fluid.compressibility = 4.5e-10;
    
    ierr = sim.initialize(config); CHKERRQ(ierr);
    
    PetscPrintf(comm, "Porosity: %.2f\n", mat.porosity);
    PetscPrintf(comm, "Initial permeability: %.1f mD\n", mat.permeability_x);
    PetscPrintf(comm, "Max permeability: %.1f mD\n", mat.max_permeability);
    PetscPrintf(comm, "Recovery time: %.1f s\n", config.permeability_recovery_time);
    
    ierr = sim.run(); CHKERRQ(ierr);
    
    PetscPrintf(comm, "Poroelastodynamic simulation completed.\n");
    
    PetscFunctionReturn(0);
}

// Example 3: Static-to-dynamic triggering (stress-induced seismicity)
PetscErrorCode runStaticTriggeredSeismicity(MPI_Comm comm) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    PetscPrintf(comm, "\n=== Example 3: Static-to-Dynamic Triggering ===\n");
    
    Simulator sim(comm);
    
    // Configuration
    SimulationConfig config;
    config.start_time = 0.0;
    config.end_time = 86400.0;    // 1 day
    config.dt_initial = 60.0;      // 1 minute for quasi-static
    config.dt_min = 0.001;         // 1 ms for dynamic event
    config.dt_max = 3600.0;        // 1 hour
    config.max_timesteps = 10000;
    config.output_frequency = 10;
    config.output_file = "triggered_seismicity";
    config.output_format = "HDF5";
    config.enable_adaptive_timestepping = true;
    
    // Static-to-dynamic triggering
    config.use_dynamic_mode = false;          // Start quasi-static
    config.use_static_triggering = true;
    config.dynamic_trigger_threshold = 5.0e6; // 5 MPa
    config.dynamic_event_duration = 10.0;     // 10 second event
    
    // Enable poroelastodynamics for realistic behavior
    config.enable_poroelastodynamics = true;
    config.solid_model = SolidModelType::POROELASTODYNAMIC;
    config.fluid_model = FluidModelType::SINGLE_COMPONENT;
    config.enable_faults = true;
    
    // Permeability changes during seismic event
    config.enable_dynamic_permeability_change = true;
    config.permeability_sensitivity = 3.0;
    config.permeability_recovery_time = 3600.0;  // 1 hour recovery
    
    // Grid
    GridConfig grid;
    grid.nx = 60;
    grid.ny = 60;
    grid.nz = 30;
    grid.Lx = 3000.0;
    grid.Ly = 3000.0;
    grid.Lz = 1500.0;
    
    // Material properties (basement rock with critically stressed fault)
    MaterialProperties mat;
    mat.youngs_modulus = 35.0e9;
    mat.poisson_ratio = 0.25;
    mat.density = 2650.0;
    mat.porosity = 0.05;
    mat.permeability_x = 1.0;  // Low permeability
    mat.biot_coefficient = 0.7;
    
    mat.p_wave_velocity = 5000.0;
    mat.s_wave_velocity = 2900.0;
    mat.quality_factor = 100.0;
    mat.damping_alpha = 0.02;
    mat.damping_beta = 0.002;
    
    // Strong permeability enhancement during seismic slip
    mat.permeability_strain_coeff = 5.0e-7;
    mat.permeability_stress_coeff = 1.0e-15;
    mat.min_permeability = 0.1;
    mat.max_permeability = 100.0;  // 100x increase possible!
    
    FluidProperties fluid;
    fluid.density = 1020.0;
    fluid.viscosity = 0.0008;
    
    ierr = sim.initialize(config); CHKERRQ(ierr);
    
    PetscPrintf(comm, "Trigger threshold: %.2f MPa\n", 
               config.dynamic_trigger_threshold / 1.0e6);
    PetscPrintf(comm, "Dynamic event duration: %.1f s\n", 
               config.dynamic_event_duration);
    PetscPrintf(comm, "Max permeability enhancement: %.1fx\n",
               mat.max_permeability / mat.permeability_x);
    
    ierr = sim.run(); CHKERRQ(ierr);
    
    PetscPrintf(comm, "Static-triggered seismicity simulation completed.\n");
    
    PetscFunctionReturn(0);
}

// Example 4: Wave-induced permeability enhancement
PetscErrorCode runWavePermeabilityEnhancement(MPI_Comm comm) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    PetscPrintf(comm, "\n=== Example 4: Wave-Induced Permeability Enhancement ===\n");
    
    Simulator sim(comm);
    
    // Configuration
    SimulationConfig config;
    config.start_time = 0.0;
    config.end_time = 1000.0;     // 1000 seconds
    config.dt_initial = 0.01;      // 10 ms
    config.dt_min = 0.001;
    config.dt_max = 1.0;
    config.max_timesteps = 50000;
    config.output_frequency = 100;
    config.output_file = "wave_perm_enhancement";
    config.output_format = "HDF5";
    
    config.use_dynamic_mode = true;
    config.enable_poroelastodynamics = true;
    config.solid_model = SolidModelType::POROELASTODYNAMIC;
    
    // Focus on permeability dynamics
    config.enable_dynamic_permeability_change = true;
    config.permeability_sensitivity = 3.0;     // Very sensitive
    config.permeability_recovery_time = 500.0; // 500 s recovery
    
    // Grid
    GridConfig grid;
    grid.nx = 100;
    grid.ny = 100;
    grid.nz = 50;
    grid.Lx = 500.0;
    grid.Ly = 500.0;
    grid.Lz = 250.0;
    
    // Fractured reservoir rock
    MaterialProperties mat;
    mat.youngs_modulus = 20.0e9;
    mat.poisson_ratio = 0.22;
    mat.density = 2450.0;
    mat.porosity = 0.15;
    mat.permeability_x = 50.0;
    mat.permeability_y = 50.0;
    mat.permeability_z = 25.0;
    mat.biot_coefficient = 0.85;
    
    mat.p_wave_velocity = 4000.0;
    mat.s_wave_velocity = 2300.0;
    mat.quality_factor = 60.0;
    mat.damping_alpha = 0.015;
    mat.damping_beta = 0.0015;
    
    // Strong permeability-strain coupling (key feature!)
    mat.permeability_strain_coeff = 2.0e-6;
    mat.permeability_stress_coeff = 1.0e-14;
    mat.min_permeability = 5.0;
    mat.max_permeability = 500.0;  // 10x enhancement
    
    FluidProperties fluid;
    fluid.density = 850.0;         // Oil
    fluid.viscosity = 0.005;       // 5 cP
    
    ierr = sim.initialize(config); CHKERRQ(ierr);
    
    PetscPrintf(comm, "Application: Enhanced oil recovery via seismic stimulation\n");
    PetscPrintf(comm, "Initial permeability: %.1f mD\n", mat.permeability_x);
    PetscPrintf(comm, "Enhanced permeability: %.1f mD\n", mat.max_permeability);
    PetscPrintf(comm, "Enhancement factor: %.1fx\n", 
               mat.max_permeability / mat.permeability_x);
    PetscPrintf(comm, "Recovery time: %.1f s\n", config.permeability_recovery_time);
    
    ierr = sim.run(); CHKERRQ(ierr);
    
    PetscPrintf(comm, "Wave permeability enhancement simulation completed.\n");
    
    PetscFunctionReturn(0);
}

// Main program
int main(int argc, char **argv) {
    PetscErrorCode ierr;
    
    ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);
    
    MPI_Comm comm = PETSC_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    if (rank == 0) {
        std::cout << "======================================\n";
        std::cout << "Wave Propagation Examples\n";
        std::cout << "======================================\n";
    }
    
    // Parse command line to select example
    PetscBool run_all = PETSC_FALSE;
    PetscInt example = 1;
    
    ierr = PetscOptionsGetInt(nullptr, nullptr, "-example", &example, nullptr); CHKERRQ(ierr);
    ierr = PetscOptionsGetBool(nullptr, nullptr, "-run_all", &run_all, nullptr); CHKERRQ(ierr);
    
    if (run_all) {
        ierr = runElastodynamicWaves(comm); CHKERRQ(ierr);
        ierr = runPoroelastodynamicWaves(comm); CHKERRQ(ierr);
        ierr = runStaticTriggeredSeismicity(comm); CHKERRQ(ierr);
        ierr = runWavePermeabilityEnhancement(comm); CHKERRQ(ierr);
    } else {
        switch(example) {
            case 1:
                ierr = runElastodynamicWaves(comm); CHKERRQ(ierr);
                break;
            case 2:
                ierr = runPoroelastodynamicWaves(comm); CHKERRQ(ierr);
                break;
            case 3:
                ierr = runStaticTriggeredSeismicity(comm); CHKERRQ(ierr);
                break;
            case 4:
                ierr = runWavePermeabilityEnhancement(comm); CHKERRQ(ierr);
                break;
            default:
                if (rank == 0) {
                    std::cerr << "Invalid example number. Choose 1-4.\n";
                }
        }
    }
    
    if (rank == 0) {
        std::cout << "\n======================================\n";
        std::cout << "All simulations completed successfully!\n";
        std::cout << "======================================\n";
        std::cout << "\nUsage:\n";
        std::cout << "  -example <1-4>  : Run specific example\n";
        std::cout << "  -run_all        : Run all examples\n";
        std::cout << "\nExamples:\n";
        std::cout << "  1. Elastodynamic waves\n";
        std::cout << "  2. Poroelastodynamic waves\n";
        std::cout << "  3. Static-triggered seismicity\n";
        std::cout << "  4. Wave permeability enhancement\n";
    }
    
    ierr = PetscFinalize(); CHKERRQ(ierr);
    return 0;
}
