/*
 * SCEC LOH.1: Layer Over Halfspace - Wave Propagation Test
 * 
 * Configuration-driven simulation of SCEC LOH.1 benchmark problem.
 * Tests wave propagation in a simple layered medium.
 * 
 * Key Features:
 * - 1 km layer over elastic halfspace
 * - Point source at depth
 * - Tests P-wave, S-wave, and surface wave propagation
 * - Comparison with analytical/reference solutions
 * 
 * Reference: SCEC Community Velocity Model (CVM) Tests
 *            Olsen et al., BSSA (2006)
 *            https://strike.scec.org/scecpedia/LOH
 * 
 * Domain:
 * - 30 km × 30 km × 17 km
 * - Layer 1: 0-1 km depth (V_p = 4000 m/s, V_s = 2000 m/s)
 * - Halfspace: >1 km depth (V_p = 6000 m/s, V_s = 3464 m/s)
 * - Point source: 2 km depth, explosion
 * 
 * Stations:
 * - 10 receivers at surface along radial lines
 * - Record 3-component velocity
 * 
 * Usage:
 *   mpirun -np 8 ./scec_loh1 -c config/scec_loh1.config
 */

#include "Simulator.hpp"
#include <iostream>

static char help[] = "SCEC LOH.1: Layer Over Halfspace Wave Propagation\n"
                     "Usage: mpirun -np N ./scec_loh1 -c <config_file>\n\n";

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, help);
    
    MPI_Comm comm = PETSC_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    char config_file[PETSC_MAX_PATH_LEN] = "config/scec_loh1.config";
    PetscOptionsGetString(nullptr, nullptr, "-c", config_file,
                         sizeof(config_file), nullptr);
    
    if (rank == 0) {
        std::cout << "========================================\n";
        std::cout << "  SCEC LOH.1 Benchmark\n";
        std::cout << "  Layer Over Halfspace - Wave Propagation\n";
        std::cout << "========================================\n\n";
        std::cout << "Config file: " << config_file << "\n\n";
        std::cout << "Problem Description:\n";
        std::cout << "  • 1 km layer over elastic halfspace\n";
        std::cout << "  • Point explosion source at 2 km depth\n";
        std::cout << "  • 10 surface receivers\n";
        std::cout << "  • Tests wave propagation accuracy\n\n";
        
        std::cout << "Material Properties:\n";
        std::cout << "  Layer (0-1 km):\n";
        std::cout << "    V_p = 4000 m/s, V_s = 2000 m/s, ρ = 2600 kg/m³\n";
        std::cout << "  Halfspace (>1 km):\n";
        std::cout << "    V_p = 6000 m/s, V_s = 3464 m/s, ρ = 2700 kg/m³\n\n";
    }
    
    FSRM::Simulator sim(comm);
    
    PetscErrorCode ierr = sim.initializeFromConfigFile(config_file);
    CHKERRQ(ierr);
    
    ierr = sim.setupDM(); CHKERRQ(ierr);
    ierr = sim.setupFields(); CHKERRQ(ierr);
    ierr = sim.setupPhysics(); CHKERRQ(ierr);
    
    // Setup layered material properties
    ierr = sim.setLayeredMaterialProperties(); CHKERRQ(ierr);
    
    // Setup point source (explosion)
    ierr = sim.setupPointSource(); CHKERRQ(ierr);
    
    // Setup receivers for seismograms
    ierr = sim.setupReceivers(); CHKERRQ(ierr);
    
    ierr = sim.setInitialConditions(); CHKERRQ(ierr);
    ierr = sim.setupTimeStepper(); CHKERRQ(ierr);
    ierr = sim.setupSolvers(); CHKERRQ(ierr);
    
    if (rank == 0) {
        std::cout << "Running SCEC LOH.1 simulation...\n";
        std::cout << "Simulating 9 seconds of wave propagation\n\n";
    }
    
    ierr = sim.run(); CHKERRQ(ierr);
    
    // Post-processing
    sim.writeSummary();
    sim.computeSeismograms();
    sim.compareWithReference("data/scec_loh1_reference.txt");
    sim.computeWaveArrivalTimes();
    sim.generatePlots();
    
    if (rank == 0) {
        std::cout << "\n========================================\n";
        std::cout << "  SCEC LOH.1 Completed Successfully\n";
        std::cout << "========================================\n\n";
        std::cout << "Results available in: output/scec_loh1/\n";
        std::cout << "  • Seismograms: loh1_stations.txt\n";
        std::cout << "  • Wave snapshots: loh1_wave_*.vtu\n";
        std::cout << "  • Comparison: loh1_comparison.png\n\n";
        std::cout << "Compare with SCEC reference solutions:\n";
        std::cout << "  https://strike.scec.org/scecpedia/LOH.1\n\n";
    }
    
    PetscFinalize();
    return 0;
}
