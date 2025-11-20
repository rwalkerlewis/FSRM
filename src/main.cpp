#include "ReservoirSim.hpp"
#include "Simulator.hpp"
#include "ConfigReader.hpp"
#include <petsc.h>
#include <iostream>
#include <string>

static char help[] = "ReservoirSim - Comprehensive Petroleum Reservoir Simulator\n"
                    "Usage: reservoirsim [options]\n\n"
                    "Options:\n"
                    "  -i <file>         Input file (Eclipse format .DATA)\n"
                    "  -c <file>         Configuration file (.config)\n"
                    "  -o <file>         Output file prefix\n"
                    "  -format <type>    Output format (ECLIPSE, VTK, HDF5)\n"
                    "  -ts_type <type>   Time stepping: euler, beuler, bdf, arkimex\n"
                    "  -snes_type <type> Nonlinear solver: newtonls, newtontr, nrichardson\n"
                    "  -ksp_type <type>  Linear solver: gmres, bcgs, cg\n"
                    "  -pc_type <type>   Preconditioner: ilu, asm, gamg, hypre\n\n"
                    "Examples:\n"
                    "  # Use configuration file (recommended)\n"
                    "  mpirun -np 4 reservoirsim -c config/shale_reservoir.config\n\n"
                    "  # Use Eclipse input\n"
                    "  mpirun -np 4 reservoirsim -i SPE1.DATA -o output/spe1\n\n"
                    "  # Generate template configuration\n"
                    "  reservoirsim -generate_config my_config.config\n\n";

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
                ResSim::ConfigReader::generateTemplate(generate_config);
                PetscPrintf(comm, "Configuration template written to: %s\n", generate_config);
                PetscPrintf(comm, "Edit this file to customize your simulation.\n");
            }
            ierr = PetscFinalize();
            return 0;
        }
        
        // Parse command line arguments
        char config_file[PETSC_MAX_PATH_LEN] = "";
        char input_file[PETSC_MAX_PATH_LEN] = "";
        char output_file[PETSC_MAX_PATH_LEN] = "output";
        char output_format[256] = "VTK";
        PetscBool config_provided = PETSC_FALSE;
        PetscBool input_provided = PETSC_FALSE;
        
        ierr = PetscOptionsGetString(nullptr, nullptr, "-c", config_file,
                                     sizeof(config_file), &config_provided); CHKERRQ(ierr);
        ierr = PetscOptionsGetString(nullptr, nullptr, "-i", input_file, 
                                     sizeof(input_file), &input_provided); CHKERRQ(ierr);
        ierr = PetscOptionsGetString(nullptr, nullptr, "-o", output_file, 
                                     sizeof(output_file), nullptr); CHKERRQ(ierr);
        ierr = PetscOptionsGetString(nullptr, nullptr, "-format", output_format, 
                                     sizeof(output_format), nullptr); CHKERRQ(ierr);
        
        if (!config_provided && !input_provided) {
            if (rank == 0) {
                PetscPrintf(comm, "Error: Configuration file (-c) or Eclipse input file (-i) required\n");
                PetscPrintf(comm, "Run with -help for usage information\n");
                PetscPrintf(comm, "Generate template: reservoirsim -generate_config template.config\n");
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
            if (config_provided) {
                PetscPrintf(comm, "Config file:   %s\n", config_file);
            }
            if (input_provided) {
                PetscPrintf(comm, "Input file:    %s\n", input_file);
            }
            PetscPrintf(comm, "Output prefix: %s\n", output_file);
            PetscPrintf(comm, "Output format: %s\n", output_format);
            PetscPrintf(comm, "\n");
        }
        
        try {
            // Create simulator
            ResSim::Simulator sim(comm);
            
            if (rank == 0) {
                PetscPrintf(comm, "Initializing simulator...\n");
            }
            
            // Initialize from config file or Eclipse input
            if (config_provided) {
                ierr = sim.initializeFromConfigFile(config_file); CHKERRQ(ierr);
            } else {
                // Configure simulation from Eclipse input
                ResSim::SimulationConfig config;
                config.input_file = input_file;
                config.output_file = output_file;
                config.output_format = output_format;
                
                ierr = sim.initialize(config); CHKERRQ(ierr);
                
                if (rank == 0) {
                    PetscPrintf(comm, "Loading Eclipse input file...\n");
                }
                
                // Load Eclipse input
                ierr = sim.loadEclipseInput(input_file); CHKERRQ(ierr);
            }
            
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
