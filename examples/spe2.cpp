/**
 * @file spe2.cpp
 * @brief SPE2 - Three-Phase Coning Problem
 * 
 * SPE Comparative Solution Project - Problem 2
 * 
 * Description:
 * - Three-phase (water, oil, gas) flow
 * - Radial grid around producing well
 * - Gas coning from gas cap
 * - Water coning from aquifer
 * - Critical rate analysis
 * 
 * Grid: 10 radial × 18 vertical layers (180 cells)
 * Duration: 1825 days (5 years)
 * 
 * Reference:
 * - Aziz et al., "Simulation of Three-Phase Coning", SPE (1985)
 */

#include "Simulator.hpp"
#include "ConfigReader.hpp"
#include <iostream>
#include <mpi.h>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    try {
        if (rank == 0) {
            std::cout << "========================================\n";
            std::cout << "SPE2: Three-Phase Coning Problem\n";
            std::cout << "========================================\n\n";
        }
        
        // Read configuration
        std::string config_file = (argc > 2 && std::string(argv[1]) == "-c") ?
            argv[2] : "config/spe2_benchmark.config";
        
        ConfigReader config(config_file);
        
        // Create simulator
        FSRM::Simulator sim(config);
        
        // Setup problem
        sim.initialize();
        
        if (rank == 0) {
            std::cout << "Problem Setup:\n";
            std::cout << "  Grid: 10 (radial) × 18 (vertical)\n";
            std::cout << "  Phases: Water, Oil, Gas\n";
            std::cout << "  Initial: Gas cap + oil zone + aquifer\n";
            std::cout << "  Well: Center, layers 9-14\n";
            std::cout << "  Focus: Coning behavior\n\n";
        }
        
        // Time stepping
        double t = 0.0;
        double t_end = 1825.0 * 86400.0;  // 5 years
        double dt = 1.0 * 86400.0;         // 1 day initial
        int step = 0;
        
        while (t < t_end) {
            // Solve timestep
            sim.solveTimestep(dt);
            
            // Output
            if (step % 30 == 0 && rank == 0) {
                std::cout << "Time = " << t/86400.0 << " days\n";
                
                // Monitor gas and water production
                double q_oil = sim.getWellOilRate(0);
                double q_gas = sim.getWellGasRate(0);
                double q_water = sim.getWellWaterRate(0);
                double GOR = (q_oil > 0) ? q_gas / q_oil : 0.0;
                double WC = (q_oil + q_water > 0) ? 
                    q_water / (q_oil + q_water) * 100.0 : 0.0;
                
                std::cout << "  Oil rate: " << q_oil * 86400.0 << " m³/day\n";
                std::cout << "  GOR: " << GOR << " (gas coning indicator)\n";
                std::cout << "  Water cut: " << WC << " % (water coning indicator)\n\n";
            }
            
            // Adaptive time stepping
            if (sim.converged()) {
                dt = std::min(dt * 1.2, 10.0 * 86400.0);
            } else {
                dt = std::max(dt * 0.5, 0.1 * 86400.0);
            }
            
            t += dt;
            step++;
        }
        
        if (rank == 0) {
            std::cout << "\n========================================\n";
            std::cout << "SPE2 Simulation Complete\n";
            std::cout << "Total timesteps: " << step << "\n";
            std::cout << "========================================\n";
        }
        
    } catch (const std::exception& e) {
        if (rank == 0) {
            std::cerr << "Error: " << e.what() << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    MPI_Finalize();
    return 0;
}
