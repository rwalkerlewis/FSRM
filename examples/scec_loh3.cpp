/**
 * @file scec_loh3.cpp
 * @brief SCEC LOH.3 - Layered Medium Wave Propagation
 * 
 * Description:
 * - Multiple horizontal layers with varying properties
 * - Interface reflections and transmissions
 * - Point explosion source
 * - Wave trapping in low-velocity layers
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
            std::cout << "SCEC LOH.3: Layered Medium\n";
            std::cout << "========================================\n\n";
        }
        
        std::string config_file = (argc > 2 && std::string(argv[1]) == "-c") ?
            argv[2] : "config/scec_loh3.config";
        
        ConfigReader config(config_file);
        FSRM::Simulator sim(config);
        sim.initialize();
        
        if (rank == 0) {
            std::cout << "Problem Setup:\n";
            std::cout << "  Domain: 40×40×20 km\n";
            std::cout << "  Grid: 200×200×100 cells\n";
            std::cout << "  Layers: 5 horizontal layers\n";
            std::cout << "    Layer 1 (0-2km): c_p=6000 m/s, c_s=3464 m/s\n";
            std::cout << "    Layer 2 (2-5km): c_p=5500 m/s, c_s=3175 m/s\n";
            std::cout << "    Layer 3 (5-10km): c_p=5000 m/s, c_s=2887 m/s\n";
            std::cout << "    Layer 4 (10-15km): c_p=4500 m/s, c_s=2598 m/s\n";
            std::cout << "    Layer 5 (15-20km): c_p=4000 m/s, c_s=2309 m/s\n";
            std::cout << "  Source: Point explosion at 10 km depth\n";
            std::cout << "  Duration: 20 seconds\n\n";
        }
        
        double t = 0.0;
        double t_end = 20.0;
        double dt = 0.002;
        int step = 0;
        
        // Vertical array of receivers
        std::vector<double> receiver_depths = {0, 2000, 5000, 10000, 15000};  // m
        
        while (t < t_end) {
            sim.solveTimestep(dt);
            
            if (step % 100 == 0 && rank == 0) {
                std::cout << "Time = " << t << " s\n";
                
                // Record at depth receivers
                for (size_t i = 0; i < receiver_depths.size(); ++i) {
                    double vel = sim.getVerticalVelocity(20000.0, 20000.0, receiver_depths[i]);
                    int layer = (i < 1) ? 1 : (i < 2) ? 2 : (i < 3) ? 3 : (i < 4) ? 4 : 5;
                    std::cout << "  Depth " << receiver_depths[i]/1000.0 << " km (Layer " << layer 
                              << "): v_z = " << vel << " m/s\n";
                }
                std::cout << "\n";
            }
            
            t += dt;
            step++;
        }
        
        if (rank == 0) {
            std::cout << "\nSCEC LOH.3 Complete - Total timesteps: " << step << "\n";
            std::cout << "Layer interface effects recorded at depth receivers\n";
        }
        
    } catch (const std::exception& e) {
        if (rank == 0) std::cerr << "Error: " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    MPI_Finalize();
    return 0;
}
