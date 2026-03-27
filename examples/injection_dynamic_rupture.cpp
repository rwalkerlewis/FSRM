/*
 * Injection-to-Seismogram Dynamic Rupture Chain
 *
 * End-to-end simulation of the complete induced seismicity chain:
 *
 * CO2/fluid injection
 *   -> pore pressure diffusion (implicit, days-months)
 *     -> Biot effective stress change
 *       -> Coulomb failure on fault (delta_CFS > 0)
 *         -> IMEX transition to explicit
 *           -> split-mesh dynamic rupture with rate-state friction
 *             -> seismic wave propagation
 *               -> synthetic seismograms at monitoring stations
 *                 -> IMEX transition back to implicit
 *                   -> continued injection
 *
 * Every link is physically connected through the PETSc solver infrastructure.
 * The fault ruptures spontaneously as a cohesive zone in the finite element mesh.
 *
 * Key components:
 * - FaultMeshManager: splits mesh along fault with DMPlexConstructCohesiveCells
 * - CohesiveFaultKernel: PetscDS callbacks for locked/friction constraint
 * - CoulombStressTransfer: FEM stress sampling at fault vertices
 * - ImplicitExplicitTransitionManager: automatic mode switching
 * - SeismometerNetwork: synthetic seismogram output (SAC/miniSEED)
 *
 * Usage:
 *   mpirun -np 4 ./injection_dynamic_rupture -c config/injection_dynamic_rupture_chain.config
 *
 * References:
 *   Segall & Lu (2015): Injection-induced seismicity, Reviews of Geophysics
 *   Aagaard et al. (2013): Domain decomposition for fault slip in FEM models
 *   Harris et al. (2018): Verifying dynamic earthquake rupture codes, SRL
 */

#include "core/Simulator.hpp"
#include <iostream>

static char help[] = "Injection-to-Seismogram Dynamic Rupture Chain\n"
                     "Usage: mpirun -np N ./injection_dynamic_rupture -c <config_file>\n\n";

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, help);

    MPI_Comm comm = PETSC_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);

    char config_file[PETSC_MAX_PATH_LEN] = "config/injection_dynamic_rupture_chain.config";
    PetscOptionsGetString(nullptr, nullptr, "-c", config_file,
                         sizeof(config_file), nullptr);

    if (rank == 0) {
        std::cout << "============================================================\n";
        std::cout << "  Injection-to-Seismogram Dynamic Rupture Chain\n";
        std::cout << "============================================================\n\n";
        std::cout << "Config file: " << config_file << "\n\n";
        std::cout << "Simulation chain:\n";
        std::cout << "  1. Fluid injection (implicit, quasi-static)\n";
        std::cout << "  2. Pore pressure diffusion + Biot coupling\n";
        std::cout << "  3. Coulomb failure detection on fault\n";
        std::cout << "  4. IMEX transition: implicit -> explicit\n";
        std::cout << "  5. Split-mesh dynamic rupture (rate-state friction)\n";
        std::cout << "  6. Seismic wave propagation\n";
        std::cout << "  7. Synthetic seismograms at monitoring stations\n";
        std::cout << "  8. IMEX transition: explicit -> implicit\n";
        std::cout << "  9. Continued injection\n\n";
    }

    FSRM::Simulator sim(comm);

    PetscErrorCode ierr = sim.initializeFromConfigFile(config_file);
    CHKERRQ(ierr);

    // Standard Simulator pipeline
    ierr = sim.setupDM(); CHKERRQ(ierr);
    ierr = sim.setupFields(); CHKERRQ(ierr);
    ierr = sim.setupPhysics(); CHKERRQ(ierr);
    ierr = sim.setMaterialProperties(); CHKERRQ(ierr);
    ierr = sim.setInitialConditions(); CHKERRQ(ierr);
    ierr = sim.setupTimeStepper(); CHKERRQ(ierr);
    ierr = sim.setupSolvers(); CHKERRQ(ierr);

    if (rank == 0) {
        std::cout << "Starting simulation...\n\n";
    }

    ierr = sim.run(); CHKERRQ(ierr);

    // Post-processing
    sim.writeSummary();
    sim.computePerformanceMetrics();

    if (rank == 0) {
        std::cout << "\n============================================================\n";
        std::cout << "  Simulation Completed Successfully\n";
        std::cout << "============================================================\n\n";
        std::cout << "Outputs:\n";
        std::cout << "  - Seismograms: output/*.sac\n";
        std::cout << "  - Field data:  output/*.h5\n";
        std::cout << "  - IMEX log:    imex_transitions.log\n";
        std::cout << "  - Summary:     simulation_summary.txt\n\n";
    }

    PetscFinalize();
    return 0;
}
