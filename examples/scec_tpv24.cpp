/**
 * @file scec_tpv24.cpp
 * @brief SCEC TPV24 - Dynamic Triggering Benchmark
 * 
 * Description:
 * - Two parallel faults
 * - Rupture on first fault triggers second fault
 * - Stress transfer and wave interaction
 * - Nucleation by passing waves
 * 
 * Grid: 192×192×96 cells (48×48×24 km)
 * Duration: 20 seconds
 * 
 * Reference:
 * - Lotto et al., "Should earthquake ruptures be modeled as 
 *   a constrained stochastic process?", JGR (2017)
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
            std::cout << "SCEC TPV24: Dynamic Triggering\n";
            std::cout << "========================================\n\n";
        }
        
        std::string config_file = (argc > 2 && std::string(argv[1]) == "-c") ?
            argv[2] : "config/scec_tpv24.config";
        
        ConfigReader config(config_file);
        FSRM::Simulator sim(config);
        sim.initialize();
        
        if (rank == 0) {
            std::cout << "Problem Setup:\n";
            std::cout << "  Domain: 48×48×24 km\n";
            std::cout << "  Grid: 192×192×96 cells\n";
            std::cout << "  Faults: 2 parallel strike-slip faults\n";
            std::cout << "  Mechanism: Fault 1 triggers Fault 2\n";
            std::cout << "  Duration: 20 seconds\n\n";
        }
        
        double t = 0.0;
        double t_end = 20.0;
        double dt = 0.001;
        int step = 0;
        bool fault2_triggered = false;
        double trigger_time = 0.0;
        
        while (t < t_end) {
            sim.solveTimestep(dt);
            
            // Check if second fault has nucleated
            if (!fault2_triggered && sim.isSlipping(1)) {  // Fault index 1
                fault2_triggered = true;
                trigger_time = t;
                if (rank == 0) {
                    std::cout << "\n*** FAULT 2 TRIGGERED at t = " << trigger_time << " s ***\n\n";
                }
            }
            
            if (step % 200 == 0 && rank == 0) {
                std::cout << "Time = " << t << " s\n";
                std::cout << "  Fault 1 slip: " << sim.getMaxSlip(0) << " m\n";
                std::cout << "  Fault 2 slip: " << sim.getMaxSlip(1) << " m\n";
                std::cout << "  Fault 2 status: " << (fault2_triggered ? "ACTIVE" : "dormant") << "\n\n";
            }
            
            t += dt;
            step++;
        }
        
        if (rank == 0) {
            std::cout << "\nSCEC TPV24 Complete - Total timesteps: " << step << "\n";
            if (fault2_triggered) {
                std::cout << "Fault 2 triggered at t = " << trigger_time << " s\n";
                std::cout << "Delay time: " << trigger_time << " s\n";
            } else {
                std::cout << "WARNING: Fault 2 did not trigger\n";
            }
        }
        
    } catch (const std::exception& e) {
        if (rank == 0) std::cerr << "Error: " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    MPI_Finalize();
    return 0;
}
