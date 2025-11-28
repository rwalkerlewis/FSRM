/**
 * @file spe11.cpp
 * @brief SPE11 - CO2 Storage CSP (Comparative Solution Project)
 * 
 * Description:
 * - CO2 injection into saline aquifer
 * - Upscaling study (fine vs coarse grids)
 * - Long-term storage (50 years)
 * - Grid refinement analysis
 * 
 * Grid: Variable (case A: 840 cells, case B: 8,400 cells, case C: 168,000 cells)
 * Duration: 50 years
 * 
 * Reference:
 * - Nordbotten et al., "The 11th SPE CSP on CO2 Storage", 
 *   Computational Geosciences (2021)
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
            std::cout << "SPE11: CO2 Storage CSP\n";
            std::cout << "========================================\n\n";
        }
        
        std::string config_file = (argc > 2 && std::string(argv[1]) == "-c") ?
            argv[2] : "config/spe11_benchmark.config";
        
        ConfigReader config(config_file);
        FSRM::Simulator sim(config);
        sim.initialize();
        
        if (rank == 0) {
            std::cout << "Problem Setup:\n";
            std::cout << "  Domain: 2.8 km × 1.2 km × 1.2 km\n";
            std::cout << "  Phases: Water + CO2 (supercritical)\n";
            std::cout << "  Wells: 1 injector at bottom\n";
            std::cout << "  Duration: 50 years (25 injection + 25 post-injection)\n";
            std::cout << "  Focus: Long-term plume migration\n\n";
        }
        
        double t = 0.0;
        double t_end = 50.0 * 365.25 * 86400.0;  // 50 years
        double t_injection = 25.0 * 365.25 * 86400.0;  // 25 years injection
        double dt = 10.0 * 86400.0;  // 10 days
        int step = 0;
        
        while (t < t_end) {
            // Stop injection after 25 years
            if (t >= t_injection && t - dt < t_injection) {
                if (rank == 0) {
                    std::cout << "\n*** Injection stopped at t = 25 years ***\n\n";
                }
                sim.setWellRate(0, 0.0);  // Stop injector
            }
            
            sim.solveTimestep(dt);
            
            if (step % 365 == 0 && rank == 0) {
                double years = t / (365.25 * 86400.0);
                std::cout << "Time = " << years << " years\n";
                std::cout << "  CO2 injected: " << sim.getTotalCO2Injected() << " tonnes\n";
                std::cout << "  Mobile CO2: " << sim.getMobileCO2Mass() << " tonnes\n";
                std::cout << "  Dissolved CO2: " << sim.getDissolvedCO2Mass() << " tonnes\n";
                std::cout << "  Plume extent: " << sim.getPlumeLateralExtent() << " m\n\n";
            }
            
            dt = sim.converged() ? std::min(dt * 1.5, 100.0 * 86400.0) : 
                                   std::max(dt * 0.5, 1.0 * 86400.0);
            t += dt;
            step++;
        }
        
        if (rank == 0) {
            std::cout << "\nSPE11 Complete - Total timesteps: " << step << "\n";
            std::cout << "Final plume statistics:\n";
            std::cout << "  Total CO2 stored: " << sim.getTotalCO2Injected() << " tonnes\n";
            std::cout << "  Maximum migration: " << sim.getPlumeLateralExtent() << " m\n";
        }
        
    } catch (const std::exception& e) {
        if (rank == 0) std::cerr << "Error: " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    MPI_Finalize();
    return 0;
}
