/**
 * @file pinatubo_plinian.cpp
 * @brief Pinatubo 1991 VEI-6 Plinian Eruption Simulation
 *
 * This example simulates a large Plinian eruption similar to the 1991 
 * eruption of Mount Pinatubo, Philippines. The eruption produced:
 *   - 35 km high eruption column
 *   - ~5 km³ of dacite tephra
 *   - Pyroclastic flows reaching 16 km from vent
 *   - Major SO₂ release affecting global climate
 *
 * The simulation includes:
 *   1. Magma chamber depressurization
 *   2. Conduit flow with fragmentation
 *   3. Eruption column dynamics
 *   4. Pyroclastic density currents from column collapse
 *   5. Tephra dispersal
 *   6. SO₂ emission
 *   7. Ground deformation
 *
 * Historical Context:
 *   Mount Pinatubo erupted on June 15, 1991, in one of the largest
 *   volcanic eruptions of the 20th century. The eruption ejected
 *   approximately 10 km³ of material, reduced global temperatures
 *   by ~0.5°C for 2 years, and produced lahars that affected
 *   hundreds of thousands of people for years afterward.
 *
 * Usage:
 *   mpirun -np 8 ./pinatubo_plinian [options]
 *
 * Options:
 *   -end_time <s>     Simulation duration (default: 36000 = 10 hours)
 *   -output <dir>     Output directory (default: output/pinatubo)
 *   -vei <int>        Scale to VEI 4-7 (default: 6)
 *
 * @author FSRM Development Team
 */

#include "VolcanoModel.hpp"
#include "ConfigReader.hpp"
#include <iostream>
#include <iomanip>
#include <cstring>
#include <chrono>

using namespace FSRM;

static const char* help = 
    "================================================================\n"
    "  Pinatubo 1991 VEI-6 Plinian Eruption Simulation\n"
    "================================================================\n"
    "\n"
    "Simulates a large Plinian eruption similar to Mount Pinatubo 1991.\n"
    "\n"
    "Usage: mpirun -np N ./pinatubo_plinian [options]\n"
    "\n"
    "Options:\n"
    "  -end_time <s>      Simulation duration in seconds (default: 36000)\n"
    "  -output <dir>      Output directory\n"
    "  -vei <int>         Scale to VEI 4, 5, 6, or 7 (default: 6)\n"
    "  -help              Show this message\n"
    "\n";

/**
 * @brief Print simulation banner
 */
void printBanner(int rank) {
    if (rank != 0) return;
    
    std::cout << "\n";
    std::cout << "╔══════════════════════════════════════════════════════════════╗\n";
    std::cout << "║                                                              ║\n";
    std::cout << "║       MOUNT PINATUBO 1991 PLINIAN ERUPTION SIMULATION        ║\n";
    std::cout << "║                                                              ║\n";
    std::cout << "║     VEI-6 Eruption: Column Height ~35 km, 5 km³ Tephra       ║\n";
    std::cout << "║                                                              ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════╝\n";
    std::cout << "\n";
    std::cout << "Historical Context:\n";
    std::cout << "  • Eruption date: June 15, 1991\n";
    std::cout << "  • VEI: 6 (Colossal)\n";
    std::cout << "  • Column height: ~35 km\n";
    std::cout << "  • Tephra volume: ~10 km³ (5 km³ DRE)\n";
    std::cout << "  • SO₂ released: ~20 Mt\n";
    std::cout << "  • Global cooling: ~0.5°C for 2 years\n";
    std::cout << "\n";
}

/**
 * @brief Scale parameters by VEI
 */
void scaleByVEI(int vei, double& volume_km3, double& MER, double& column_height) {
    // VEI scaling based on Newhall & Self (1982)
    switch(vei) {
        case 4:
            volume_km3 = 0.1;
            MER = 1e7;
            column_height = 15000.0;
            break;
        case 5:
            volume_km3 = 1.0;
            MER = 5e7;
            column_height = 25000.0;
            break;
        case 6:
            volume_km3 = 10.0;
            MER = 2e8;
            column_height = 35000.0;
            break;
        case 7:
            volume_km3 = 100.0;
            MER = 1e9;
            column_height = 45000.0;
            break;
        default:
            volume_km3 = 10.0;
            MER = 2e8;
            column_height = 35000.0;
    }
}

/**
 * @brief Setup Pinatubo-style magma chamber
 */
MagmaChamberGeometry setupChamber() {
    MagmaChamberGeometry chamber;
    chamber.x_center = 0.0;
    chamber.y_center = 0.0;
    chamber.depth = 6000.0;          // 6 km depth
    chamber.volume = 40e9;           // 40 km³ (pre-eruption)
    chamber.semi_axis_a = 3000.0;    // 3 km horizontal
    chamber.semi_axis_b = 3000.0;
    chamber.semi_axis_c = 2000.0;    // 2 km vertical
    chamber.shape = DeformationSourceType::OBLATE_SPHEROID;
    return chamber;
}

/**
 * @brief Setup Pinatubo dacite magma properties
 */
MagmaProperties setupMagma() {
    MagmaProperties magma;
    magma.composition = MagmaComposition::DACITE;
    
    // Pinatubo dacite composition (Pallister et al., 1996)
    magma.SiO2 = 64.5;
    magma.Al2O3 = 16.0;
    magma.FeO = 4.5;
    magma.MgO = 2.2;
    magma.CaO = 5.5;
    magma.Na2O = 4.3;
    magma.K2O = 2.0;
    
    // High volatile content (pre-eruption)
    magma.H2O_total = 6.0;           // 6 wt% water
    magma.CO2_total = 0.05;
    magma.SO2_total = 0.2;           // High sulfur
    
    // Temperature and crystallinity
    magma.temperature = 1053.15;     // 780°C (relatively cool)
    magma.crystal_fraction = 0.35;   // ~35% crystite
    
    return magma;
}

/**
 * @brief Setup volcanic conduit
 */
ConduitGeometry setupConduit() {
    ConduitGeometry conduit;
    conduit.depth = 6000.0;          // 6 km to chamber
    conduit.length = 6000.0;
    conduit.radius_base = 50.0;      // 50 m at chamber
    conduit.radius_vent = 100.0;     // 100 m at vent (flaring)
    return conduit;
}

/**
 * @brief Simple topography for Pinatubo region
 */
double pinatuboTopography(double x, double y) {
    // Simplified conical volcano
    double r = std::sqrt(x*x + y*y);
    double summit_height = 1745.0;   // Pre-eruption height
    double base_radius = 15000.0;    // 15 km radius
    
    if (r < base_radius) {
        // Conical profile
        return summit_height * (1.0 - r/base_radius);
    }
    return 0.0;
}

/**
 * @brief Print eruption progress
 */
void printProgress(const CoupledVolcanoSystem& volcano, double time, int rank) {
    if (rank != 0) return;
    
    static double last_print = -3600.0;
    if (time - last_print < 300.0 && time > 0) return;  // Print every 5 minutes
    last_print = time;
    
    double hours = time / 3600.0;
    double MER = volcano.getMassEruptionRate();
    double column_height = volcano.getColumnHeight();
    
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "t = " << hours << " hr: ";
    std::cout << "MER = " << std::scientific << std::setprecision(2) << MER << " kg/s, ";
    std::cout << "Column = " << std::fixed << std::setprecision(1) << column_height/1000.0 << " km";
    
    if (volcano.isErupting()) {
        std::cout << " [ERUPTING]";
    }
    std::cout << "\n";
}

int main(int argc, char** argv) {
    // Initialize MPI
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    // Parse command line
    std::string output_dir = "output/pinatubo";
    double end_time = 36000.0;  // 10 hours
    int vei = 6;
    
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-help") == 0 || strcmp(argv[i], "--help") == 0) {
            if (rank == 0) std::cout << help;
            MPI_Finalize();
            return 0;
        } else if (strcmp(argv[i], "-end_time") == 0 && i + 1 < argc) {
            end_time = std::stod(argv[++i]);
        } else if (strcmp(argv[i], "-output") == 0 && i + 1 < argc) {
            output_dir = argv[++i];
        } else if (strcmp(argv[i], "-vei") == 0 && i + 1 < argc) {
            vei = std::stoi(argv[++i]);
        }
    }
    
    printBanner(rank);
    
    auto start_wall = std::chrono::high_resolution_clock::now();
    
    try {
        // =====================================================================
        // 1. Setup Volcanic System
        // =====================================================================
        
        // Get scaled parameters
        double volume_km3, MER, column_height;
        scaleByVEI(vei, volume_km3, MER, column_height);
        
        if (rank == 0) {
            std::cout << "Simulation Parameters:\n";
            std::cout << "  VEI: " << vei << "\n";
            std::cout << "  Target volume: " << volume_km3 << " km³\n";
            std::cout << "  Peak MER: " << std::scientific << MER << " kg/s\n";
            std::cout << "  Expected column height: " << std::fixed << column_height/1000.0 << " km\n";
            std::cout << "\n";
        }
        
        // Setup components
        MagmaChamberGeometry chamber = setupChamber();
        MagmaProperties magma = setupMagma();
        ConduitGeometry conduit = setupConduit();
        
        // Configure the system
        VolcanoSystemConfig config;
        config.enable_magma_chamber = true;
        config.enable_conduit_flow = true;
        config.enable_eruption_column = true;
        config.enable_pdc = true;
        config.enable_deformation = true;
        config.enable_seismicity = true;
        config.enable_gas_emission = true;
        config.enable_tephra = true;
        config.couple_chamber_conduit = true;
        config.couple_conduit_column = true;
        config.couple_column_pdc = true;
        config.couple_chamber_deformation = true;
        config.end_time = end_time;
        
        // =====================================================================
        // 2. Initialize and Run
        // =====================================================================
        
        CoupledVolcanoSystem volcano;
        volcano.initialize(config, chamber, magma, conduit);
        volcano.setTopography(pinatuboTopography);
        
        // Trigger Plinian eruption
        if (rank == 0) {
            std::cout << "═══════════════════════════════════════════════════════════════\n";
            std::cout << "                    TRIGGERING ERUPTION                         \n";
            std::cout << "═══════════════════════════════════════════════════════════════\n\n";
        }
        
        volcano.triggerEruption(EruptionType::PLINIAN);
        
        // Run with progress output
        double dt = 10.0;  // 10 second timestep
        double current_time = 0.0;
        
        while (current_time < end_time) {
            volcano.update(dt);
            current_time += dt;
            printProgress(volcano, current_time, rank);
        }
        
        // =====================================================================
        // 3. Output Results
        // =====================================================================
        
        auto end_wall = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(
            end_wall - start_wall).count();
        
        if (rank == 0) {
            std::cout << "\n";
            std::cout << "═══════════════════════════════════════════════════════════════\n";
            std::cout << "                    SIMULATION COMPLETE                          \n";
            std::cout << "═══════════════════════════════════════════════════════════════\n\n";
            
            std::cout << "Wall time: " << duration << " seconds\n";
            std::cout << "Speed-up: " << end_time / duration << "x real-time\n\n";
            
            std::cout << "Writing output to: " << output_dir << "\n";
        }
        
        volcano.writeOutput(output_dir);
        
        // Print summary
        if (rank == 0) {
            auto monitoring = volcano.getMonitoringTimeSeries();
            
            double max_MER = 0.0;
            double max_column = 0.0;
            double total_SO2 = 0.0;
            
            for (const auto& m : monitoring) {
                max_MER = std::max(max_MER, m.mass_eruption_rate);
                max_column = std::max(max_column, m.column_height);
                total_SO2 += m.SO2_emission_rate * 10.0;  // kg (10s interval)
            }
            
            std::cout << "\nEruption Summary:\n";
            std::cout << "  Peak MER: " << std::scientific << max_MER << " kg/s\n";
            std::cout << "  Max column height: " << std::fixed << std::setprecision(1) 
                      << max_column/1000.0 << " km\n";
            std::cout << "  Total SO₂: " << total_SO2/1e12 << " Mt\n\n";
            
            std::cout << "Output files:\n";
            std::cout << "  • " << output_dir << "/monitoring_timeseries.dat\n";
            std::cout << "  • " << output_dir << "/deformation_field.dat\n";
            std::cout << "  • " << output_dir << "/tephra_deposit.dat\n";
            std::cout << "  • " << output_dir << "/pdc_inundation.dat\n";
            std::cout << "  • " << output_dir << "/gas_concentration.dat\n";
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error on rank " << rank << ": " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    MPI_Finalize();
    return 0;
}
