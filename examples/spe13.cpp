/**
 * @file spe13.cpp
 * @brief SPE13 - Well Control and Constraints
 * 
 * Description:
 * - Three-phase black oil
 * - Complex well control switching
 * - Rate and pressure constraints
 * - Group controls
 * 
 * Grid: 24×25×15 (9,000 cells) - same as SPE9
 * Duration: Variable based on well controls
 * 
 * Reference:
 * - Killough, "Ninth SPE CSP", SPE 29110 (1995)
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
            std::cout << "SPE13: Well Control and Constraints\n";
            std::cout << "========================================\n\n";
        }
        
        std::string config_file = (argc > 2 && std::string(argv[1]) == "-c") ?
            argv[2] : "config/spe13_benchmark.config";
        
        ConfigReader config(config_file);
        FSRM::Simulator sim(config);
        sim.initialize();
        
        if (rank == 0) {
            std::cout << "Problem Setup:\n";
            std::cout << "  Grid: 24×25×15 (9,000 cells)\n";
            std::cout << "  Wells: 26 total (producers and injectors)\n";
            std::cout << "  Controls: Mixed rate/BHP constraints\n";
            std::cout << "  Focus: Well control switching logic\n\n";
        }
        
        double t = 0.0;
        double t_end = 3000.0 * 86400.0;  // ~8 years
        double dt = 1.0 * 86400.0;
        int step = 0;
        int control_switches = 0;
        
        while (t < t_end) {
            // Check and apply well constraints
            int switches = sim.applyWellConstraints();
            control_switches += switches;
            
            sim.solveTimestep(dt);
            
            if (step % 30 == 0 && rank == 0) {
                std::cout << "Time = " << t/86400.0 << " days\n";
                std::cout << "  Field oil rate: " << sim.getFieldOilRate() * 86400.0 << " m³/day\n";
                std::cout << "  Field water rate: " << sim.getFieldWaterRate() * 86400.0 << " m³/day\n";
                std::cout << "  Field gas rate: " << sim.getFieldGasRate() * 86400.0 << " m³/day\n";
                std::cout << "  Control switches: " << control_switches << "\n\n";
            }
            
            dt = sim.converged() ? std::min(dt * 1.2, 10.0 * 86400.0) : 
                                   std::max(dt * 0.5, 0.1 * 86400.0);
            t += dt;
            step++;
        }
        
        if (rank == 0) {
            std::cout << "\nSPE13 Complete\n";
            std::cout << "  Total timesteps: " << step << "\n";
            std::cout << "  Total control switches: " << control_switches << "\n";
            std::cout << "  Average switches/year: " << control_switches / (t / (365.25 * 86400.0)) << "\n";
        }
        
    } catch (const std::exception& e) {
        if (rank == 0) std::cerr << "Error: " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    MPI_Finalize();
    return 0;
}
