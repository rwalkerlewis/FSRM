/**
 * @file scec_loh2.cpp
 * @brief SCEC LOH.2 - Basin Edge Effects
 * 
 * Description:
 * - Wave propagation near sedimentary basin edge
 * - Surface waves and basin amplification
 * - Point explosion source
 * - Edge-generated waves
 * 
 * Grid: 200×200×100 cells (40×40×20 km)
 * Duration: 20 seconds
 * 
 * Reference:
 * - Olsen et al., "Strong shaking in Los Angeles expected from 
 *   southern San Andreas earthquake", GRL (2006)
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
            std::cout << "SCEC LOH.2: Basin Edge Effects\n";
            std::cout << "========================================\n\n";
        }
        
        std::string config_file = (argc > 2 && std::string(argv[1]) == "-c") ?
            argv[2] : "config/scec_loh2.config";
        
        ConfigReader config(config_file);
        FSRM::Simulator sim(config);
        sim.initialize();
        
        if (rank == 0) {
            std::cout << "Problem Setup:\n";
            std::cout << "  Domain: 40×40×20 km\n";
            std::cout << "  Grid: 200×200×100 cells\n";
            std::cout << "  Basin: Sedimentary layer near edge\n";
            std::cout << "  Source: Point explosion at depth\n";
            std::cout << "  Duration: 20 seconds\n\n";
            std::cout << "Focus: Basin edge wave generation\n\n";
        }
        
        double t = 0.0;
        double t_end = 20.0;
        double dt = 0.002;  // 2 ms
        int step = 0;
        
        // Receiver locations
        std::vector<double> receiver_x = {10000, 15000, 20000, 25000, 30000};  // m
        
        while (t < t_end) {
            sim.solveTimestep(dt);
            
            if (step % 100 == 0 && rank == 0) {
                std::cout << "Time = " << t << " s\n";
                
                // Record at receivers
                for (size_t i = 0; i < receiver_x.size(); ++i) {
                    double vel = sim.getSurfaceVelocity(receiver_x[i], 20000.0);
                    std::cout << "  Receiver " << i+1 << " (x=" << receiver_x[i]/1000.0 
                              << " km): v = " << vel << " m/s\n";
                }
                std::cout << "\n";
            }
            
            t += dt;
            step++;
        }
        
        if (rank == 0) {
            std::cout << "\nSCEC LOH.2 Complete - Total timesteps: " << step << "\n";
            std::cout << "Basin amplification effects recorded at surface receivers\n";
        }
        
    } catch (const std::exception& e) {
        if (rank == 0) std::cerr << "Error: " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    MPI_Finalize();
    return 0;
}
