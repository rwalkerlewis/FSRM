/*
 * SCEC TPV5: Strike-Slip Benchmark with Heterogeneous Initial Stress
 * 
 * Configuration-driven simulation of SCEC TPV5 benchmark problem.
 * Tests spontaneous dynamic rupture propagation on a vertical strike-slip fault.
 * 
 * Key Features:
 * - 3D elastic half-space
 * - Vertical planar fault (strike-slip)
 * - Heterogeneous initial stress with asperity
 * - Spontaneous rupture nucleation
 * - Linear slip-weakening friction
 * 
 * Reference: SCEC Dynamic Rupture Code Verification Exercise
 *            Harris et al., Seism. Res. Lett. (2009)
 *            https://strike.scec.org/cvws/
 * 
 * Domain:
 * - 48 km × 48 km × 24 km (along-strike × across-fault × depth)
 * - Fault: 30 km × 15 km in center
 * - Nucleation zone: 3 km × 3 km square asperity
 * 
 * Usage:
 *   mpirun -np 8 ./scec_tpv5 -c config/scec_tpv5.config
 */

#include "Simulator.hpp"
#include <iostream>

static char help[] = "SCEC TPV5: Strike-Slip Dynamic Rupture Benchmark\n"
                     "Usage: mpirun -np N ./scec_tpv5 -c <config_file>\n\n";

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, help);
    
    MPI_Comm comm = PETSC_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    // Get config file from command line
    char config_file[PETSC_MAX_PATH_LEN] = "config/scec_tpv5.config";
    PetscOptionsGetString(nullptr, nullptr, "-c", config_file,
                         sizeof(config_file), nullptr);
    
    if (rank == 0) {
        std::cout << "========================================\n";
        std::cout << "  SCEC TPV5 Benchmark\n";
        std::cout << "  Dynamic Rupture on Strike-Slip Fault\n";
        std::cout << "========================================\n\n";
        std::cout << "Config file: " << config_file << "\n\n";
        std::cout << "Problem Description:\n";
        std::cout << "  • 3D elastic half-space\n";
        std::cout << "  • Vertical strike-slip fault (30 km × 15 km)\n";
        std::cout << "  • Heterogeneous initial stress\n";
        std::cout << "  • Linear slip-weakening friction\n";
        std::cout << "  • Spontaneous rupture nucleation\n\n";
        
        std::cout << "Physics:\n";
        std::cout << "  • V_p = 6000 m/s (P-wave velocity)\n";
        std::cout << "  • V_s = 3464 m/s (S-wave velocity)\n";
        std::cout << "  • ρ = 2670 kg/m³ (density)\n";
        std::cout << "  • μ_s = 0.677 (static friction)\n";
        std::cout << "  • μ_d = 0.525 (dynamic friction)\n";
        std::cout << "  • D_c = 0.40 m (slip-weakening distance)\n\n";
    }
    
    // Create and run simulation from config
    FSRM::Simulator sim(comm);
    
    PetscErrorCode ierr = sim.initializeFromConfigFile(config_file);
    CHKERRQ(ierr);
    
    ierr = sim.setupDM(); CHKERRQ(ierr);
    ierr = sim.setupFields(); CHKERRQ(ierr);
    ierr = sim.setupPhysics(); CHKERRQ(ierr);
    ierr = sim.setMaterialProperties(); CHKERRQ(ierr);
    
    // Set heterogeneous initial stress with nucleation asperity
    ierr = sim.setInitialStress(); CHKERRQ(ierr);
    
    // Setup fault with slip-weakening friction
    ierr = sim.setupFault(); CHKERRQ(ierr);
    
    ierr = sim.setInitialConditions(); CHKERRQ(ierr);
    ierr = sim.setupTimeStepper(); CHKERRQ(ierr);
    ierr = sim.setupSolvers(); CHKERRQ(ierr);
    
    if (rank == 0) {
        std::cout << "Running SCEC TPV5 simulation...\n";
        std::cout << "Simulating 12 seconds of rupture propagation\n\n";
    }
    
    ierr = sim.run(); CHKERRQ(ierr);
    
    // Post-processing specific to dynamic rupture
    sim.writeSummary();
    sim.computeRuptureMetrics();  // Rupture speed, slip distribution, etc.
    sim.computeSeismograms();     // Synthetic seismograms at stations
    sim.generatePlots();
    
    if (rank == 0) {
        std::cout << "\n========================================\n";
        std::cout << "  SCEC TPV5 Completed Successfully\n";
        std::cout << "========================================\n\n";
        std::cout << "Results available in: output/scec_tpv5/\n";
        std::cout << "  • Slip distribution: tpv5_slip_*.vtu\n";
        std::cout << "  • Rupture time: tpv5_rupture_time.vtu\n";
        std::cout << "  • Seismograms: tpv5_stations.txt\n";
        std::cout << "  • Summary data: tpv5_SUMMARY.txt\n\n";
        std::cout << "Compare with SCEC reference solutions:\n";
        std::cout << "  https://strike.scec.org/cvws/tpv5docs.html\n\n";
    }
    
    PetscFinalize();
    return 0;
}
