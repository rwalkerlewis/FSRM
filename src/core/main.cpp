#include "core/FSRM.hpp"
#include "core/Simulator.hpp"
#include "core/ConfigReader.hpp"
#include <petsc.h>
#include <iostream>
#include <string>

static char help[] = "FSRM - Full Service Reservoir Model\n"
                    "Coupled multiphysics simulator for nuclear explosion monitoring,\n"
                    "seismic wave propagation, dynamic fault rupture, and THM poroelasticity.\n\n"
                    "Usage: fsrm -c <config_file> [PETSc options]\n\n"
                    "Options:\n"
                    "  -c <file>         Configuration file (.config)\n"
                    "  -o <file>         Output file prefix\n"
                    "  -ts_type <type>   Time stepping: beuler, bdf, alpha2\n"
                    "  -snes_type <type> Nonlinear solver: newtonls, newtontr\n"
                    "  -ksp_type <type>  Linear solver: gmres, bcgs, cg\n"
                    "  -pc_type <type>   Preconditioner: ilu, asm, gamg, hypre\n\n"
                    "Examples:\n"
                    "  ./fsrm -c config/examples/uniaxial_compression.config\n"
                    "  mpirun -np 4 ./fsrm -c config/examples/punggye_ri_layered.config\n\n"
                    "  # Generate template configuration\n"
                    "  ./fsrm -generate_config my_config.config\n\n";

int main(int argc, char** argv) {
    PetscErrorCode ierr;
    
    // Initialize PETSc
    ierr = PetscInitialize(&argc, &argv, nullptr, help); CHKERRQ(ierr);
    
    {
        MPI_Comm comm = PETSC_COMM_WORLD;
        int rank;
        MPI_Comm_rank(comm, &rank);
        
        // Check for config file generation
        char generate_config[PETSC_MAX_PATH_LEN] = "";
        PetscBool gen_config;
        ierr = PetscOptionsGetString(nullptr, nullptr, "-generate_config", generate_config,
                                     sizeof(generate_config), &gen_config); CHKERRQ(ierr);
        
        if (gen_config) {
            if (rank == 0) {
                FSRM::ConfigReader::generateTemplate(generate_config);
                PetscPrintf(comm, "Configuration template written to: %s\n", generate_config);
                PetscPrintf(comm, "Edit this file to customize your simulation.\n");
            }
            ierr = PetscFinalize();
            return 0;
        }
        
        // Parse command line arguments
        char config_file[PETSC_MAX_PATH_LEN] = "";
        char output_file[PETSC_MAX_PATH_LEN] = "output";
        PetscBool config_provided = PETSC_FALSE;
        
        ierr = PetscOptionsGetString(nullptr, nullptr, "-c", config_file,
                                     sizeof(config_file), &config_provided); CHKERRQ(ierr);
        ierr = PetscOptionsGetString(nullptr, nullptr, "-o", output_file, 
                                     sizeof(output_file), nullptr); CHKERRQ(ierr);
        
        if (!config_provided) {
            if (rank == 0) {
                PetscPrintf(comm, "Error: Configuration file (-c) required\n");
                PetscPrintf(comm, "Run with -help for usage information\n");
                PetscPrintf(comm, "Generate template: ./fsrm -generate_config template.config\n");
            }
            ierr = PetscFinalize();
            return 1;
        }
        
        if (rank == 0) {
            PetscPrintf(comm, "\n");
            PetscPrintf(comm, "============================================================\n");
            PetscPrintf(comm, "  FSRM - Full Service Reservoir Model\n");
            PetscPrintf(comm, "  Coupled Multiphysics Simulator\n");
            PetscPrintf(comm, "  Version 1.0.0\n");
            PetscPrintf(comm, "============================================================\n");
            PetscPrintf(comm, "\n");
            PetscPrintf(comm, "Config file:   %s\n", config_file);
            PetscPrintf(comm, "Output prefix: %s\n", output_file);
            PetscPrintf(comm, "\n");
        }
        
        try {
            // Create simulator
            FSRM::Simulator sim(comm);
            
            if (rank == 0) {
                PetscPrintf(comm, "Initializing simulator...\n");
            }
            
            // Initialize from config file
            ierr = sim.initializeFromConfigFile(config_file); CHKERRQ(ierr);
            
            if (rank == 0) {
                PetscPrintf(comm, "Setting up problem...\n");
            }
            
            // Setup
            ierr = sim.setupDM(); CHKERRQ(ierr);
            ierr = sim.setupFaultNetwork(); CHKERRQ(ierr);  // Must be after setupDM, before setupFields
            ierr = sim.labelBoundaries(); CHKERRQ(ierr);     // Label boundary faces after mesh setup
            ierr = sim.setupFields(); CHKERRQ(ierr);         // Calls setupBoundaryConditions internally
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
            }
            
            // Write final summary and seismometer output
            ierr = sim.writeSummary(); CHKERRQ(ierr);
            
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
