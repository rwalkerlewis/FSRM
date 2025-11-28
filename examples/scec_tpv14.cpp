/**
 * @file scec_tpv14.cpp
 * @brief SCEC TPV14 - Bimaterial Fault Benchmark
 * 
 * Description:
 * - Material contrast across fault
 * - Asymmetric rupture propagation
 * - Preferred rupture direction
 * - Stress rotation effects
 * 
 * Grid: 240×192×96 cells (60×48×24 km)
 * Duration: 15 seconds
 * 
 * Reference:
 * - Dalguer & Day, "Staggered-grid split-node method for spontaneous
 *   rupture simulation", JGR (2007)
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
            std::cout << "SCEC TPV14: Bimaterial Fault\n";
            std::cout << "========================================\n\n";
        }
        
        std::string config_file = (argc > 2 && std::string(argv[1]) == "-c") ?
            argv[2] : "config/scec_tpv14.config";
        
        ConfigReader config(config_file);
        FSRM::Simulator sim(config);
        sim.initialize();
        
        if (rank == 0) {
            std::cout << "Problem Setup:\n";
            std::cout << "  Domain: 60×48×24 km\n";
            std::cout << "  Grid: 240×192×96 cells\n";
            std::cout << "  Fault: Vertical strike-slip\n";
            std::cout << "  Material: Bimaterial contrast (ρ, μ different)\n";
            std::cout << "  Duration: 15 seconds\n\n";
            std::cout << "Expected: Asymmetric rupture propagation\n\n";
        }
        
        double t = 0.0;
        double t_end = 15.0;
        double dt = 0.001;
        int step = 0;
        
        while (t < t_end) {
            sim.solveTimestep(dt);
            
            if (step % 200 == 0) {
                double rupture_extent_plus = sim.getRuptureExtent("+");
                double rupture_extent_minus = sim.getRuptureExtent("-");
                double asymmetry = rupture_extent_plus / rupture_extent_minus;
                
                if (rank == 0) {
                    std::cout << "Time = " << t << " s\n";
                    std::cout << "  Rupture extent (+): " << rupture_extent_plus << " km\n";
                    std::cout << "  Rupture extent (-): " << rupture_extent_minus << " km\n";
                    std::cout << "  Asymmetry ratio: " << asymmetry << "\n\n";
                }
            }
            
            t += dt;
            step++;
        }
        
        if (rank == 0) {
            std::cout << "\nSCEC TPV14 Complete - Total timesteps: " << step << "\n";
        }
        
    } catch (const std::exception& e) {
        if (rank == 0) std::cerr << "Error: " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    MPI_Finalize();
    return 0;
}
