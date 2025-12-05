/**
 * @file kilauea_effusive.cpp
 * @brief Kilauea-style Effusive Basaltic Eruption Simulation
 *
 * This example simulates an effusive basaltic eruption similar to the
 * 2018 Kilauea Lower East Rift Zone eruption in Hawaii. The eruption
 * features:
 *   - Low-viscosity basaltic lava flows
 *   - Fire fountaining (Hawaiian-style)
 *   - Lava flow propagation over complex topography
 *   - SO₂ and laze (lava haze) production
 *   - Ground deformation from magma transport
 *
 * The simulation includes:
 *   1. Magma chamber and rift zone plumbing
 *   2. Lava flow rheology and spreading
 *   3. Thermal evolution and crust formation
 *   4. Volcanic gas emissions
 *   5. Ground deformation monitoring
 *
 * Historical Context:
 *   The 2018 Kilauea eruption began on May 3, 2018, with fissure
 *   eruptions in Leilani Estates. Over 3 months, the eruption
 *   produced ~1 km³ of lava, destroyed 716 structures, and created
 *   new land. The summit experienced dramatic collapse events.
 *
 * Usage:
 *   mpirun -np 4 ./kilauea_effusive [options]
 *
 * Options:
 *   -duration <days>  Simulation duration in days (default: 30)
 *   -output <dir>     Output directory
 *   -effusion <m3/s>  Effusion rate in m³/s DRE (default: 100)
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
    "  Kilauea-style Effusive Basaltic Eruption Simulation\n"
    "================================================================\n"
    "\n"
    "Simulates an effusive basaltic eruption with lava flows.\n"
    "\n"
    "Usage: mpirun -np N ./kilauea_effusive [options]\n"
    "\n"
    "Options:\n"
    "  -duration <days>   Simulation duration in days (default: 30)\n"
    "  -output <dir>      Output directory (default: output/kilauea)\n"
    "  -effusion <m3/s>   Peak effusion rate (default: 100)\n"
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
    std::cout << "║          KILAUEA EFFUSIVE ERUPTION SIMULATION                ║\n";
    std::cout << "║                                                              ║\n";
    std::cout << "║        Hawaiian-Style Basaltic Lava Flow Modeling            ║\n";
    std::cout << "║                                                              ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════╝\n";
    std::cout << "\n";
    std::cout << "Historical Context (2018 Kilauea):\n";
    std::cout << "  • Eruption start: May 3, 2018\n";
    std::cout << "  • Location: Lower East Rift Zone (Leilani Estates)\n";
    std::cout << "  • Duration: ~3 months\n";
    std::cout << "  • Lava volume: ~1 km³\n";
    std::cout << "  • Peak effusion rate: ~200 m³/s\n";
    std::cout << "  • Structures destroyed: 716\n";
    std::cout << "\n";
}

/**
 * @brief Setup Kilauea-style magma reservoir (shallow)
 */
MagmaChamberGeometry setupKilaueaChamber() {
    MagmaChamberGeometry chamber;
    chamber.x_center = 0.0;
    chamber.y_center = 0.0;
    chamber.depth = 2000.0;          // 2 km depth (shallow)
    chamber.volume = 5e9;            // 5 km³
    chamber.semi_axis_a = 2000.0;    // 2 km horizontal
    chamber.semi_axis_b = 2000.0;
    chamber.semi_axis_c = 500.0;     // 500 m vertical (sill-like)
    chamber.shape = DeformationSourceType::OBLATE_SPHEROID;
    return chamber;
}

/**
 * @brief Setup Kilauea basalt properties
 */
MagmaProperties setupBasalt() {
    MagmaProperties magma;
    magma.composition = MagmaComposition::BASALT;
    
    // Kilauea tholeiitic basalt
    magma.SiO2 = 50.0;
    magma.Al2O3 = 13.5;
    magma.FeO = 11.5;
    magma.MgO = 7.5;
    magma.CaO = 11.0;
    magma.Na2O = 2.3;
    magma.K2O = 0.5;
    
    // Moderate volatile content
    magma.H2O_total = 0.5;           // 0.5 wt% water (low)
    magma.CO2_total = 0.03;
    magma.SO2_total = 0.15;          // High sulfur
    
    // Hot and fluid
    magma.temperature = 1423.15;     // 1150°C
    magma.crystal_fraction = 0.05;   // Few crystals (olivine)
    
    return magma;
}

/**
 * @brief Setup volcanic rift zone conduit
 */
ConduitGeometry setupRiftConduit() {
    ConduitGeometry conduit;
    conduit.depth = 2000.0;          // 2 km to reservoir
    conduit.length = 2000.0;
    conduit.radius_base = 5.0;       // 5 m dike
    conduit.radius_vent = 20.0;      // 20 m at fissure
    return conduit;
}

/**
 * @brief Simple topography for Hawaiian volcano flank
 * 
 * Creates a gently sloping shield volcano profile
 */
double hawaiianTopography(double x, double y) {
    // Gentle shield volcano slope (~5 degrees)
    double r = std::sqrt(x*x + y*y);
    double summit_height = 1200.0;   // Summit elevation
    double base_radius = 50000.0;    // 50 km radius
    double slope = 0.05;             // ~3 degrees
    
    // Rift zone runs in x direction, slope in y direction
    double rift_y = 0.0;             // Rift zone location
    double base_elevation = 100.0;   // Coastal elevation
    
    // Distance from rift zone to coast
    double dist_to_coast = 10000.0;  // 10 km to ocean
    
    // Slope from rift zone to ocean
    if (y > rift_y) {
        // Upslope from rift
        return base_elevation + (summit_height - base_elevation) * 
               std::min(1.0, (y - rift_y) / 10000.0);
    } else {
        // Downslope to ocean
        double progress = (rift_y - y) / dist_to_coast;
        if (progress > 1.0) {
            return 0.0;  // Ocean
        }
        return base_elevation * (1.0 - progress);
    }
}

/**
 * @brief Print lava flow progress
 */
void printFlowProgress(const CoupledVolcanoSystem& volcano, 
                       const LavaFlowModel& lava,
                       double time, int rank) {
    if (rank != 0) return;
    
    static double last_print = -86400.0;
    if (time - last_print < 3600.0 && time > 0) return;  // Print every hour
    last_print = time;
    
    double days = time / 86400.0;
    double flow_area = lava.getFlowArea() / 1e6;        // km²
    double flow_length = lava.getFlowLength() / 1000.0; // km
    double flow_volume = lava.getFlowVolume() / 1e6;    // Mm³
    
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "Day " << days << ": ";
    std::cout << "Area = " << flow_area << " km², ";
    std::cout << "Length = " << flow_length << " km, ";
    std::cout << "Volume = " << std::setprecision(3) << flow_volume << " Mm³\n";
}

int main(int argc, char** argv) {
    // Initialize MPI
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    // Parse command line
    std::string output_dir = "output/kilauea";
    double duration_days = 30.0;
    double effusion_rate = 100.0;  // m³/s
    
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-help") == 0 || strcmp(argv[i], "--help") == 0) {
            if (rank == 0) std::cout << help;
            MPI_Finalize();
            return 0;
        } else if (strcmp(argv[i], "-duration") == 0 && i + 1 < argc) {
            duration_days = std::stod(argv[++i]);
        } else if (strcmp(argv[i], "-output") == 0 && i + 1 < argc) {
            output_dir = argv[++i];
        } else if (strcmp(argv[i], "-effusion") == 0 && i + 1 < argc) {
            effusion_rate = std::stod(argv[++i]);
        }
    }
    
    printBanner(rank);
    
    double end_time = duration_days * 86400.0;  // Convert to seconds
    
    auto start_wall = std::chrono::high_resolution_clock::now();
    
    try {
        // =====================================================================
        // 1. Setup Volcanic System
        // =====================================================================
        
        if (rank == 0) {
            std::cout << "Simulation Parameters:\n";
            std::cout << "  Duration: " << duration_days << " days\n";
            std::cout << "  Effusion rate: " << effusion_rate << " m³/s\n";
            std::cout << "  Expected volume: " << effusion_rate * end_time / 1e9 
                      << " km³\n";
            std::cout << "\n";
        }
        
        // Setup components
        MagmaChamberGeometry chamber = setupKilaueaChamber();
        MagmaProperties magma = setupBasalt();
        ConduitGeometry conduit = setupRiftConduit();
        
        // Configure system - effusive eruption
        VolcanoSystemConfig config;
        config.enable_magma_chamber = true;
        config.enable_conduit_flow = true;
        config.enable_eruption_column = false;  // No significant column
        config.enable_lava_flow = true;         // Main focus
        config.enable_deformation = true;
        config.enable_seismicity = true;
        config.enable_gas_emission = true;
        config.enable_tephra = false;           // Minimal tephra
        config.couple_chamber_conduit = true;
        config.couple_chamber_deformation = true;
        config.end_time = end_time;
        
        // =====================================================================
        // 2. Initialize Lava Flow Model
        // =====================================================================
        
        LavaSourceParameters lava_source;
        lava_source.effusion_rate = effusion_rate;
        lava_source.temperature = magma.temperature;
        lava_source.vent_x = 0.0;
        lava_source.vent_y = 0.0;  // On rift zone
        lava_source.vent_z = 100.0;
        lava_source.initial_thickness = 5.0;
        
        LavaFlowModel lava;
        lava.initialize(lava_source, magma);
        lava.setTopography(hawaiianTopography);
        
        // Initialize coupled system
        CoupledVolcanoSystem volcano;
        volcano.initialize(config, chamber, magma, conduit);
        volcano.setTopography(hawaiianTopography);
        
        // =====================================================================
        // 3. Run Simulation
        // =====================================================================
        
        if (rank == 0) {
            std::cout << "═══════════════════════════════════════════════════════════════\n";
            std::cout << "                  STARTING EFFUSIVE ERUPTION                    \n";
            std::cout << "═══════════════════════════════════════════════════════════════\n\n";
        }
        
        volcano.triggerEruption(EruptionType::EFFUSIVE);
        
        // Run with hourly output
        double dt = 60.0;  // 1 minute timestep
        double current_time = 0.0;
        
        while (current_time < end_time) {
            volcano.update(dt);
            lava.update(dt);
            current_time += dt;
            printFlowProgress(volcano, lava, current_time, rank);
        }
        
        // =====================================================================
        // 4. Output Results
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
            
            // Final statistics
            double final_area = lava.getFlowArea() / 1e6;
            double final_length = lava.getFlowLength() / 1000.0;
            double final_volume = lava.getFlowVolume() / 1e9;
            
            std::cout << "Eruption Summary:\n";
            std::cout << "  Duration: " << duration_days << " days\n";
            std::cout << "  Lava volume: " << std::fixed << std::setprecision(3) 
                      << final_volume << " km³\n";
            std::cout << "  Flow area: " << std::setprecision(1) << final_area << " km²\n";
            std::cout << "  Flow length: " << final_length << " km\n";
            std::cout << "\n";
            
            std::cout << "Writing output to: " << output_dir << "\n";
        }
        
        volcano.writeOutput(output_dir);
        
        if (rank == 0) {
            std::cout << "\nOutput files:\n";
            std::cout << "  • " << output_dir << "/lava_flow_extent.dat\n";
            std::cout << "  • " << output_dir << "/lava_thickness.dat\n";
            std::cout << "  • " << output_dir << "/temperature_field.dat\n";
            std::cout << "  • " << output_dir << "/deformation_field.dat\n";
            std::cout << "  • " << output_dir << "/SO2_emission.dat\n";
            std::cout << "  • " << output_dir << "/monitoring_timeseries.dat\n";
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error on rank " << rank << ": " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    MPI_Finalize();
    return 0;
}
