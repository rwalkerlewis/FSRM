/*
 * Induced Seismicity with Implicit-Explicit Time Integration
 * 
 * This example demonstrates automatic transition between quasi-static
 * (implicit) and dynamic (explicit) time integration for modeling
 * injection-induced seismicity.
 * 
 * Physical Scenario:
 * - Wastewater injection increases pore pressure over months
 * - Pore pressure diffuses toward a critically-stressed fault
 * - When fault fails, dynamic rupture occurs (seconds)
 * - Seismic waves propagate through the domain
 * - After waves decay, return to quasi-static modeling
 * 
 * Usage:
 *   mpirun -np 4 ./induced_seismicity_imex config/induced_seismicity_imex.config
 */

#include <petscsys.h>
#include "Simulator.hpp"
#include "ConfigReader.hpp"
#include "ImplicitExplicitTransition.hpp"
#include <iostream>

using namespace FSRM;

int main(int argc, char** argv) {
    PetscErrorCode ierr;
    
    // Initialize PETSc
    ierr = PetscInitialize(&argc, &argv, nullptr, 
        "Induced Seismicity with IMEX Time Integration\n"
        "Usage: mpirun -np N induced_seismicity_imex config_file.config\n"); 
    if (ierr) return ierr;
    
    MPI_Comm comm = PETSC_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    // Get config file from command line
    std::string config_file = "config/induced_seismicity_imex.config";
    if (argc > 1) {
        config_file = argv[1];
    }
    
    if (rank == 0) {
        std::cout << "\n";
        std::cout << "============================================================\n";
        std::cout << "  Induced Seismicity with IMEX Time Integration\n";
        std::cout << "============================================================\n";
        std::cout << "\n";
        std::cout << "This simulation demonstrates automatic transitions between:\n";
        std::cout << "  - IMPLICIT mode: Quasi-static pore pressure diffusion\n";
        std::cout << "  - EXPLICIT mode: Dynamic rupture and wave propagation\n";
        std::cout << "\n";
        std::cout << "Configuration file: " << config_file << "\n";
        std::cout << "\n";
    }
    
    // Create simulator
    Simulator sim(comm);
    
    // Load configuration
    ierr = sim.initializeFromConfigFile(config_file); CHKERRQ(ierr);
    
    // Parse IMEX configuration
    ConfigReader reader;
    if (reader.loadFile(config_file)) {
        ConfigReader::IMEXConfig imex_config;
        if (reader.parseIMEXConfig(imex_config)) {
            if (rank == 0) {
                std::cout << "IMEX Configuration:\n";
                std::cout << "  Enabled: " << (imex_config.enabled ? "yes" : "no") << "\n";
                std::cout << "  Initial mode: " << imex_config.initial_mode << "\n";
                std::cout << "  Trigger type: " << imex_config.trigger_type << "\n";
                std::cout << "  Settling type: " << imex_config.settling_type << "\n";
                std::cout << "\n";
                std::cout << "  Implicit dt: [" << imex_config.implicit_dt_min 
                         << ", " << imex_config.implicit_dt_max << "] s\n";
                std::cout << "  Explicit dt: [" << imex_config.explicit_dt_min 
                         << ", " << imex_config.explicit_dt_max << "] s\n";
                std::cout << "\n";
                std::cout << "  Min dynamic duration: " << imex_config.min_dynamic_duration << " s\n";
                std::cout << "  Max dynamic duration: " << imex_config.max_dynamic_duration << " s\n";
                std::cout << "\n";
            }
            
            // Setup simulation components
            ierr = sim.setupDM(); CHKERRQ(ierr);
            ierr = sim.setupFields(); CHKERRQ(ierr);
            ierr = sim.setupPhysics(); CHKERRQ(ierr);
            ierr = sim.setupTimeStepper(); CHKERRQ(ierr);
            ierr = sim.setupSolvers(); CHKERRQ(ierr);
            
            // Setup fault network for seismicity
            ierr = sim.setupFaultNetwork(); CHKERRQ(ierr);
            
            // Setup IMEX transition manager
            ierr = sim.setupIMEXTransition(imex_config); CHKERRQ(ierr);
            
            // Set initial conditions
            ierr = sim.setInitialConditions(); CHKERRQ(ierr);
            ierr = sim.setMaterialProperties(); CHKERRQ(ierr);
            
            if (rank == 0) {
                std::cout << "============================================================\n";
                std::cout << "  Starting simulation...\n";
                std::cout << "============================================================\n";
                std::cout << "\n";
                std::cout << "Watch for IMEX transitions in the output:\n";
                std::cout << "  'IMEX TRANSITION: IMPLICIT -> EXPLICIT' = fault failure\n";
                std::cout << "  'IMEX TRANSITION: EXPLICIT -> IMPLICIT' = event settled\n";
                std::cout << "\n";
            }
            
            // Run simulation
            ierr = sim.run(); CHKERRQ(ierr);
            
            // Write final outputs
            ierr = sim.writeSummary(); CHKERRQ(ierr);
            
            if (rank == 0) {
                std::cout << "\n";
                std::cout << "============================================================\n";
                std::cout << "  Simulation completed\n";
                std::cout << "============================================================\n";
                std::cout << "\n";
                std::cout << "Output files:\n";
                std::cout << "  - Visualization: output/induced_seismicity_imex_*.vtu\n";
                std::cout << "  - IMEX log: output/imex_transitions.log\n";
                std::cout << "  - Seismic catalog: output/seismic_catalog.csv\n";
                std::cout << "\n";
            }
        } else {
            if (rank == 0) {
                std::cerr << "Warning: No [IMEX] section found in config file.\n";
                std::cerr << "Running without IMEX transition (standard time stepping).\n";
            }
            
            // Run without IMEX
            ierr = sim.setupDM(); CHKERRQ(ierr);
            ierr = sim.setupFields(); CHKERRQ(ierr);
            ierr = sim.setupPhysics(); CHKERRQ(ierr);
            ierr = sim.setupTimeStepper(); CHKERRQ(ierr);
            ierr = sim.setupSolvers(); CHKERRQ(ierr);
            ierr = sim.setInitialConditions(); CHKERRQ(ierr);
            ierr = sim.setMaterialProperties(); CHKERRQ(ierr);
            ierr = sim.run(); CHKERRQ(ierr);
        }
    } else {
        SETERRQ(comm, PETSC_ERR_FILE_OPEN, "Failed to load configuration file");
    }
    
    ierr = PetscFinalize();
    return ierr;
}

/*
 * IMEX Transition Flow:
 * 
 *                    +------------------+
 *                    | Start (IMPLICIT) |
 *                    +--------+---------+
 *                             |
 *                             v
 *               +-------------+-------------+
 *               | Quasi-static time step    |
 *               | dt ~ hours-days           |
 *               | Solve: ∇·σ' - α∇p = 0    |<---------+
 *               +-------------+-------------+          |
 *                             |                        |
 *                             v                        |
 *                    +--------+--------+               |
 *                    | Check trigger   |               |
 *                    | CFF >= 0 ?      |               |
 *                    +--------+--------+               |
 *                             |                        |
 *                     No      |      Yes               |
 *                    +--------+--------+               |
 *                    |                 |               |
 *                    v                 v               |
 *               +----+----+   +--------+--------+      |
 *               | Continue|   | Switch to       |      |
 *               | implicit|   | EXPLICIT mode   |      |
 *               +---------+   +--------+--------+      |
 *                                      |               |
 *                                      v               |
 *               +----------------------+------------+  |
 *               | Dynamic time step                 |  |
 *               | dt ~ ms (CFL limited)             |  |
 *               | Solve: ρ ü = ∇·σ - α∇p           |  |
 *               +----------------------+------------+  |
 *                                      |               |
 *                                      v               |
 *                             +--------+--------+      |
 *                             | Check settling  |      |
 *                             | KE < threshold? |      |
 *                             +--------+--------+      |
 *                                      |               |
 *                              No      |     Yes       |
 *                             +--------+--------+      |
 *                             |                 |      |
 *                             v                 v      |
 *                    +--------+----+   +--------+------+
 *                    | Continue    |   | Switch back   |
 *                    | explicit    |   | to IMPLICIT   +-----+
 *                    +-------------+   +---------------+
 * 
 * Key Physics:
 * 
 * IMPLICIT (Quasi-static):
 *   - Darcy flow: q = -(k/μ)∇p
 *   - Poroelasticity: ∇·σ' = α∇p
 *   - No inertia terms
 *   - Large time steps (hours-days)
 * 
 * EXPLICIT (Dynamic):
 *   - Full momentum equation: ρ ü + damping = ∇·σ - α∇p
 *   - Wave propagation: P-waves, S-waves, Biot waves
 *   - Fault friction: rate-state law
 *   - CFL-limited time steps (milliseconds)
 * 
 * Transition Criteria:
 * 
 *   IMPLICIT -> EXPLICIT:
 *   - Coulomb Failure: CFF = τ - μ(σn - p) - c >= 0
 *   - Slip rate: V > 10⁻³ m/s (seismic threshold)
 *   - Velocity: |v| > 10⁻⁴ m/s
 * 
 *   EXPLICIT -> IMPLICIT:
 *   - Kinetic energy: KE < 1% of peak
 *   - Velocity: |v| < 10⁻⁶ m/s
 *   - Slip rate: V < 10⁻⁹ m/s (interseismic)
 *   - All criteria combined
 */
