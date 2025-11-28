/**
 * @file spe5.cpp
 * @brief SPE5 - Sixth SPE Comparative Solution Project
 * 
 * Description:
 * - Four-component compositional simulation
 * - Volatile oil and gas condensate
 * - Phase behavior modeling
 * - One injector, four producers
 * 
 * Grid: 7×7×3 (147 cells)
 * Duration: 1,500 days
 * 
 * Reference:
 * - Killough & Kossack, "Fifth Comparative Solution Project", SPE 16000 (1987)
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
            std::cout << "SPE5: Compositional Volatile Oil/Gas\n";
            std::cout << "========================================\n\n";
        }
        
        std::string config_file = (argc > 2 && std::string(argv[1]) == "-c") ?
            argv[2] : "config/spe5_benchmark.config";
        
        ConfigReader config(config_file);
        FSRM::Simulator sim(config);
        sim.initialize();
        
        if (rank == 0) {
            std::cout << "Problem Setup:\n";
            std::cout << "  Grid: 7×7×3 (147 cells)\n";
            std::cout << "  Components: C1, C3, C6, C10\n";
            std::cout << "  Wells: 1 injector (center), 4 producers (corners)\n";
            std::cout << "  Duration: 1,500 days\n\n";
        }
        
        double t = 0.0;
        double t_end = 1500.0 * 86400.0;
        double dt = 1.0 * 86400.0;
        int step = 0;
        
        while (t < t_end) {
            sim.solveTimestep(dt);
            
            if (step % 30 == 0 && rank == 0) {
                std::cout << "Time = " << t/86400.0 << " days\n";
                std::cout << "  Oil production: " << sim.getTotalOilProduction() << " m³\n";
                std::cout << "  Gas production: " << sim.getTotalGasProduction() << " m³\n\n";
            }
            
            dt = sim.converged() ? std::min(dt * 1.2, 10.0 * 86400.0) : 
                                   std::max(dt * 0.5, 0.1 * 86400.0);
            t += dt;
            step++;
        }
        
        if (rank == 0) {
            std::cout << "\nSPE5 Complete - Total timesteps: " << step << "\n";
        }
        
    } catch (const std::exception& e) {
        if (rank == 0) std::cerr << "Error: " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    MPI_Finalize();
    return 0;
}
