#include "ReservoirSim.hpp"
#include "Simulator.hpp"
#include <petsc.h>
#include <iostream>
#include <string>

static char help[] = "ReservoirSim - Comprehensive Petroleum Reservoir Simulator\n"
                    "Usage: reservoirsim -i input.DATA [options]\n\n"
                    "Options:\n"
                    "  -i <file>         Input file (Eclipse format)\n"
                    "  -o <file>         Output file prefix\n"
                    "  -format <type>    Output format (ECLIPSE, VTK, HDF5)\n"
                    "  -ts_type <type>   Time stepping: euler, beuler, bdf, arkimex\n"
                    "  -snes_type <type> Nonlinear solver: newtonls, newtontr, nrichardson\n"
                    "  -ksp_type <type>  Linear solver: gmres, bcgs, cg\n"
                    "  -pc_type <type>   Preconditioner: ilu, asm, gamg, hypre\n\n"
                    "Example:\n"
                    "  mpirun -np 4 reservoirsim -i SPE1.DATA -o output/spe1\n\n";

int main(int argc, char** argv) {
    PetscErrorCode ierr;
    
    // Initialize PETSc
    ierr = PetscInitialize(&argc, &argv, nullptr, help); CHKERRQ(ierr);
    
    {
        MPI_Comm comm = PETSC_COMM_WORLD;
        int rank;
        MPI_Comm_rank(comm, &rank);
        
        // Parse command line arguments
        char input_file[PETSC_MAX_PATH_LEN] = "";
        char output_file[PETSC_MAX_PATH_LEN] = "output";
        char output_format[256] = "ECLIPSE";
        PetscBool input_provided = PETSC_FALSE;
        
        ierr = PetscOptionsGetString(nullptr, nullptr, "-i", input_file, 
                                     sizeof(input_file), &input_provided); CHKERRQ(ierr);
        ierr = PetscOptionsGetString(nullptr, nullptr, "-o", output_file, 
                                     sizeof(output_file), nullptr); CHKERRQ(ierr);
        ierr = PetscOptionsGetString(nullptr, nullptr, "-format", output_format, 
                                     sizeof(output_format), nullptr); CHKERRQ(ierr);
        
        if (!input_provided) {
            if (rank == 0) {
                PetscPrintf(comm, "Error: Input file required. Use -i <filename>\n");
                PetscPrintf(comm, "Run with -help for usage information\n");
            }
            ierr = PetscFinalize();
            return 1;
        }
        
        if (rank == 0) {
            PetscPrintf(comm, "\n");
            PetscPrintf(comm, "============================================================\n");
            PetscPrintf(comm, "  ReservoirSim - Coupled Petroleum Reservoir Simulator\n");
            PetscPrintf(comm, "  Version 1.0.0\n");
            PetscPrintf(comm, "============================================================\n");
            PetscPrintf(comm, "\n");
            PetscPrintf(comm, "Input file:    %s\n", input_file);
            PetscPrintf(comm, "Output prefix: %s\n", output_file);
            PetscPrintf(comm, "Output format: %s\n", output_format);
            PetscPrintf(comm, "\n");
        }
        
        try {
            // Create simulator
            ResSim::Simulator sim(comm);
            
            // Configure simulation
            ResSim::SimulationConfig config;
            config.input_file = input_file;
            config.output_file = output_file;
            config.output_format = output_format;
            
            if (rank == 0) {
                PetscPrintf(comm, "Initializing simulator...\n");
            }
            
            // Initialize
            ierr = sim.initialize(config); CHKERRQ(ierr);
            
            if (rank == 0) {
                PetscPrintf(comm, "Loading input file...\n");
            }
            
            // Load Eclipse input
            ierr = sim.loadEclipseInput(input_file); CHKERRQ(ierr);
            
            if (rank == 0) {
                PetscPrintf(comm, "Setting up problem...\n");
            }
            
            // Setup
            ierr = sim.setupDM(); CHKERRQ(ierr);
            ierr = sim.setupFields(); CHKERRQ(ierr);
            ierr = sim.setupPhysics(); CHKERRQ(ierr);
            ierr = sim.setupTimeStepper(); CHKERRQ(ierr);
            ierr = sim.setupSolvers(); CHKERRQ(ierr);
            
            if (rank == 0) {
                PetscPrintf(comm, "Setting initial conditions...\n");
            }
            
            ierr = sim.setInitialConditions(); CHKERRQ(ierr);
            
            if (rank == 0) {
                PetscPrintf(comm, "\n");
                PetscPrintf(comm, "Starting simulation...\n");
                PetscPrintf(comm, "------------------------------------------------------------\n");
            }
            
            // Run simulation
            double start_time = MPI_Wtime();
            ierr = sim.run(); CHKERRQ(ierr);
            double end_time = MPI_Wtime();
            
            if (rank == 0) {
                PetscPrintf(comm, "------------------------------------------------------------\n");
                PetscPrintf(comm, "Simulation completed successfully!\n");
                PetscPrintf(comm, "Total wall time: %.2f seconds\n", end_time - start_time);
                PetscPrintf(comm, "\n");
                PetscPrintf(comm, "Generating output...\n");
            }
            
            // Write final summary
            ierr = sim.writeSummary(); CHKERRQ(ierr);
            ierr = sim.computePerformanceMetrics(); CHKERRQ(ierr);
            ierr = sim.generatePlots(); CHKERRQ(ierr);
            ierr = sim.writeStatistics(); CHKERRQ(ierr);
            
            if (rank == 0) {
                PetscPrintf(comm, "Output files written to: %s.*\n", output_file);
                PetscPrintf(comm, "\n");
                PetscPrintf(comm, "============================================================\n");
            }
            
        } catch (const std::exception& e) {
            if (rank == 0) {
                PetscPrintf(comm, "\nError: %s\n", e.what());
            }
            ierr = PetscFinalize();
            return 1;
        }
    }
    
    // Finalize PETSc
    ierr = PetscFinalize();
    return 0;
}
