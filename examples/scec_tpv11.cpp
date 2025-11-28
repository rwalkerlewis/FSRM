/**
 * @file scec_tpv11.cpp
 * @brief SCEC TPV11 - Supershear Rupture Benchmark
 * 
 * Description:
 * - Dynamic rupture propagating at supershear velocities
 * - Rupture speed exceeds shear wave speed
 * - Low strength fault
 * - Mach front formation
 * 
 * Grid: 192×192×96 cells (48×48×24 km)
 * Duration: 12 seconds
 * 
 * Reference:
 * - Harris et al., "The SCEC/USGS Dynamic Earthquake Rupture Code 
 *   Verification Exercise", Seismological Research Letters (2011)
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
            std::cout << "SCEC TPV11: Supershear Rupture\n";
            std::cout << "========================================\n\n";
        }
        
        std::string config_file = (argc > 2 && std::string(argv[1]) == "-c") ?
            argv[2] : "config/scec_tpv11.config";
        
        ConfigReader config(config_file);
        FSRM::Simulator sim(config);
        sim.initialize();
        
        if (rank == 0) {
            std::cout << "Problem Setup:\n";
            std::cout << "  Domain: 48×48×24 km\n";
            std::cout << "  Grid: 192×192×96 cells\n";
            std::cout << "  Fault: Vertical strike-slip\n";
            std::cout << "  Friction: Low strength (supershear enabled)\n";
            std::cout << "  Duration: 12 seconds\n\n";
            std::cout << "Expected rupture speed: > c_s (supershear)\n\n";
        }
        
        double t = 0.0;
        double t_end = 12.0;
        double dt = 0.001;  // 1 ms
        int step = 0;
        double max_slip = 0.0;
        double rupture_speed = 0.0;
        
        while (t < t_end) {
            sim.solveTimestep(dt);
            
            if (step % 100 == 0) {
                max_slip = sim.getMaxSlip();
                rupture_speed = sim.getRuptureSpeed();
                
                if (rank == 0) {
                    std::cout << "Time = " << t << " s\n";
                    std::cout << "  Max slip: " << max_slip << " m\n";
                    std::cout << "  Rupture speed: " << rupture_speed << " m/s\n";
                    std::cout << "  Speed/c_s: " << rupture_speed / 3000.0 << "\n";
                    if (rupture_speed > 3000.0) {
                        std::cout << "  >>> SUPERSHEAR regime <<<\n";
                    }
                    std::cout << "\n";
                }
            }
            
            t += dt;
            step++;
        }
        
        if (rank == 0) {
            std::cout << "\nSCEC TPV11 Complete - Total timesteps: " << step << "\n";
            std::cout << "Maximum rupture speed: " << rupture_speed << " m/s\n";
            std::cout << "Ratio to shear wave speed: " << rupture_speed / 3000.0 << "\n";
        }
        
    } catch (const std::exception& e) {
        if (rank == 0) std::cerr << "Error: " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    MPI_Finalize();
    return 0;
}
