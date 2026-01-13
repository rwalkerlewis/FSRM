/**
 * @file cascadia_tsunami.cpp
 * @brief Cascadia Subduction Zone M9.0 Earthquake Tsunami Simulation
 *
 * This example simulates a hypothetical M9.0 earthquake on the Cascadia
 * Subduction Zone (CSZ) and the resulting tsunami propagating along the
 * West Coast of North America.
 *
 * The simulation includes:
 *   1. Realistic CSZ fault geometry based on USGS SLAB2.0
 *   2. Slip distribution from paleoseismic constraints (Wang et al., 2013)
 *   3. Coseismic seafloor deformation using Okada (1985) model
 *   4. Tsunami propagation using nonlinear shallow water equations
 *   5. Synthetic time series at NOAA tide gauge and DART buoy locations
 *
 * Historical Context:
 *   The Cascadia Subduction Zone last ruptured on January 26, 1700 CE,
 *   producing a great earthquake (estimated Mw 9.0) that generated a
 *   trans-Pacific tsunami observed in Japan. Japanese historical records
 *   document the "orphan tsunami" arrival. Paleoseismic studies suggest
 *   similar events recur every 200-600 years.
 *
 * Scientific References:
 *   - Satake et al. (2003) - Tsunami in Japan from 1700 Cascadia event
 *   - Goldfinger et al. (2012) - Turbidite paleoseismology
 *   - Wang et al. (2013) - Heterogeneous rupture models
 *   - Witter et al. (2013) - Paleotsunami deposits in Oregon
 *   - USGS SLAB2.0 - Subduction zone geometry model
 *
 * Usage:
 *   mpirun -np 8 ./cascadia_tsunami [options]
 *
 * Options:
 *   -c <config>   Configuration file (default: config/cascadia_m9_tsunami.config)
 *   -scenario <name>  Predefined scenario: full_margin, southern, northern, central
 *   -end_time <s>     Simulation end time in seconds (default: 14400 = 4 hours)
 *   -output <dir>     Output directory (default: output/cascadia_tsunami)
 *
 * @author FSRM Development Team
 */

#include "TsunamiModel.hpp"
#include "ConfigReader.hpp"
#include <iostream>
#include <iomanip>
#include <cstring>
#include <chrono>

using namespace FSRM;

static const char* help = 
    "================================================================\n"
    "  Cascadia M9.0 Earthquake Tsunami Simulation\n"
    "================================================================\n"
    "\n"
    "Simulates a hypothetical M9.0 earthquake on the Cascadia Subduction\n"
    "Zone and the resulting tsunami along the West Coast of North America.\n"
    "\n"
    "Usage: mpirun -np N ./cascadia_tsunami [options]\n"
    "\n"
    "Options:\n"
    "  -c <file>          Configuration file\n"
    "  -scenario <name>   Predefined scenario:\n"
    "                       full_margin   - Full CSZ rupture (M9.0+)\n"
    "                       southern      - Southern segment (M8.0-8.5)\n"
    "                       northern      - Northern segment (M8.0-8.5)\n"
    "                       central       - Central Oregon (M8.0)\n"
    "                       worst_case    - Maximum credible scenario\n"
    "  -end_time <s>      Simulation duration in seconds\n"
    "  -output <dir>      Output directory\n"
    "  -synthetic_bathy   Use synthetic bathymetry (no data files needed)\n"
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
    std::cout << "║        CASCADIA SUBDUCTION ZONE TSUNAMI SIMULATION           ║\n";
    std::cout << "║                                                              ║\n";
    std::cout << "║   Hypothetical M9.0 Great Earthquake and Resulting Tsunami   ║\n";
    std::cout << "║                                                              ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════╝\n";
    std::cout << "\n";
    std::cout << "Historical Context:\n";
    std::cout << "  • Last great Cascadia earthquake: January 26, 1700 CE\n";
    std::cout << "  • Estimated magnitude: Mw 9.0\n";
    std::cout << "  • Tsunami observed in Japan ('orphan tsunami')\n";
    std::cout << "  • Average recurrence: ~200-600 years\n";
    std::cout << "  • Time since last event: 324 years\n";
    std::cout << "\n";
}

/**
 * @brief Print fault model summary
 */
void printFaultModelSummary(const CascadiaFaultModel& model, int rank) {
    if (rank != 0) return;
    
    std::cout << "Fault Model:\n";
    std::cout << "  Extent: " << std::fixed << std::setprecision(1);
    std::cout << model.south_latitude << "°N to " << model.north_latitude << "°N\n";
    std::cout << "  Length: ~" << static_cast<int>(
        TsunamiUtils::haversineDistance(
            model.trench_longitude, model.south_latitude,
            model.trench_longitude, model.north_latitude) / 1000) << " km\n";
    std::cout << "  Rupture width: " << model.rupture_width << " km\n";
    std::cout << "  Average slip: " << model.average_slip << " m\n";
    std::cout << "  Peak slip: " << model.peak_slip << " m\n";
    std::cout << "  Dip: " << model.shallow_dip << "° (shallow) to " 
              << model.deep_dip << "° (deep)\n";
    std::cout << "  Rupture velocity: " << model.rupture_velocity << " km/s\n";
    std::cout << "  Estimated magnitude: Mw " << std::setprecision(2) 
              << model.getMagnitude() << "\n";
    std::cout << "  Total moment: " << std::scientific << std::setprecision(2)
              << model.getTotalMoment() << " N·m\n";
    std::cout << "\n";
}

/**
 * @brief Print expected tsunami arrivals
 */
void printExpectedArrivals(int rank) {
    if (rank != 0) return;
    
    std::cout << "Expected Tsunami Arrivals (approximate):\n";
    std::cout << "  ┌────────────────────────┬──────────────┬────────────────┐\n";
    std::cout << "  │ Location               │ Arrival Time │ Est. Height    │\n";
    std::cout << "  ├────────────────────────┼──────────────┼────────────────┤\n";
    std::cout << "  │ DART Buoys (offshore)  │   15-25 min  │   0.3-0.5 m    │\n";
    std::cout << "  │ Crescent City, CA      │   20-30 min  │   5-10 m       │\n";
    std::cout << "  │ Oregon Coast           │   20-35 min  │   5-10 m       │\n";
    std::cout << "  │ Westport, WA           │   25-35 min  │   6-10 m       │\n";
    std::cout << "  │ Neah Bay, WA           │   30-40 min  │   4-8 m        │\n";
    std::cout << "  │ Tofino, BC             │   25-35 min  │   4-8 m        │\n";
    std::cout << "  │ Seattle (Puget Sound)  │   2-3 hours  │   1-3 m        │\n";
    std::cout << "  │ San Francisco, CA      │   1.5-2 hrs  │   1-3 m        │\n";
    std::cout << "  │ Los Angeles, CA        │   3-4 hours  │   0.5-2 m      │\n";
    std::cout << "  └────────────────────────┴──────────────┴────────────────┘\n";
    std::cout << "\n";
}

/**
 * @brief Select fault model based on scenario name
 */
CascadiaFaultModel selectScenario(const std::string& scenario) {
    if (scenario == "full_margin" || scenario == "m9" || scenario == "default") {
        return CascadiaScenarios::fullMarginM9();
    } else if (scenario == "southern") {
        return CascadiaScenarios::southernSegmentM8();
    } else if (scenario == "northern") {
        return CascadiaScenarios::northernSegmentM8();
    } else if (scenario == "central") {
        return CascadiaScenarios::centralOregonM8();
    } else if (scenario == "worst_case") {
        return CascadiaScenarios::worstCaseScenario();
    } else {
        std::cerr << "Warning: Unknown scenario '" << scenario 
                  << "', using full_margin\n";
        return CascadiaScenarios::fullMarginM9();
    }
}

/**
 * @brief Setup West Coast gauge network
 */
std::vector<TsunamiGauge> setupGaugeNetwork() {
    std::vector<TsunamiGauge> gauges;
    
    // DART buoys (deep ocean - first arrivals)
    gauges.push_back(WestCoastGaugeNetwork::dart_46404());
    gauges.push_back(WestCoastGaugeNetwork::dart_46407());
    gauges.push_back(WestCoastGaugeNetwork::dart_46419());
    
    // Tide gauges (coastal)
    gauges.push_back(WestCoastGaugeNetwork::crescent_city());
    gauges.push_back(WestCoastGaugeNetwork::astoria());
    gauges.push_back(WestCoastGaugeNetwork::westport());
    gauges.push_back(WestCoastGaugeNetwork::seattle());
    gauges.push_back(WestCoastGaugeNetwork::victoria());
    gauges.push_back(WestCoastGaugeNetwork::tofino());
    gauges.push_back(WestCoastGaugeNetwork::san_francisco());
    gauges.push_back(WestCoastGaugeNetwork::monterey());
    gauges.push_back(WestCoastGaugeNetwork::los_angeles());
    gauges.push_back(WestCoastGaugeNetwork::san_diego());
    
    return gauges;
}

/**
 * @brief Print simulation results summary
 */
void printResultsSummary(const CoupledTsunamiEarthquake& simulation, int rank) {
    if (rank != 0) return;
    
    std::cout << "\n";
    std::cout << "═══════════════════════════════════════════════════════════════\n";
    std::cout << "                    SIMULATION RESULTS SUMMARY                  \n";
    std::cout << "═══════════════════════════════════════════════════════════════\n\n";
    
    // Get arrival times and max amplitudes
    auto arrivals = simulation.getArrivalTimes();
    auto amplitudes = simulation.getMaxAmplitudes();
    
    std::cout << "Station Results:\n";
    std::cout << "  ┌────────────────────────┬──────────────┬────────────────┐\n";
    std::cout << "  │ Station                │ Arrival (min)│ Max Height (m) │\n";
    std::cout << "  ├────────────────────────┼──────────────┼────────────────┤\n";
    
    for (const auto& [name, arrival] : arrivals) {
        double height = 0.0;
        if (amplitudes.count(name)) {
            height = amplitudes.at(name);
        }
        
        std::cout << "  │ " << std::left << std::setw(22) << name 
                  << " │ " << std::right << std::setw(12) << std::fixed 
                  << std::setprecision(1);
        
        if (arrival > 0) {
            std::cout << arrival / 60.0;
        } else {
            std::cout << "     N/A   ";
        }
        
        std::cout << " │ " << std::setw(14) << std::setprecision(2) 
                  << height << " │\n";
    }
    
    std::cout << "  └────────────────────────┴──────────────┴────────────────┘\n";
    std::cout << "\n";
}

int main(int argc, char** argv) {
    // Initialize MPI
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    // Parse command line arguments
    std::string config_file = "config/cascadia_m9_tsunami.config";
    std::string scenario = "full_margin";
    std::string output_dir = "output/cascadia_tsunami";
    double end_time = 14400.0;  // 4 hours
    bool use_synthetic_bathy = true;  // Default to synthetic for demo
    
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-help") == 0 || strcmp(argv[i], "--help") == 0) {
            if (rank == 0) std::cout << help;
            MPI_Finalize();
            return 0;
        } else if (strcmp(argv[i], "-c") == 0 && i + 1 < argc) {
            config_file = argv[++i];
        } else if (strcmp(argv[i], "-scenario") == 0 && i + 1 < argc) {
            scenario = argv[++i];
        } else if (strcmp(argv[i], "-end_time") == 0 && i + 1 < argc) {
            end_time = std::stod(argv[++i]);
        } else if (strcmp(argv[i], "-output") == 0 && i + 1 < argc) {
            output_dir = argv[++i];
        } else if (strcmp(argv[i], "-synthetic_bathy") == 0) {
            use_synthetic_bathy = true;
        }
    }
    
    // Print banner
    printBanner(rank);
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    try {
        // =====================================================================
        // 1. Setup Cascadia Fault Model
        // =====================================================================
        
        if (rank == 0) {
            std::cout << "Scenario: " << scenario << "\n\n";
        }
        
        CascadiaFaultModel fault_model = selectScenario(scenario);
        printFaultModelSummary(fault_model, rank);
        
        // Generate subfaults
        auto subfaults = fault_model.generateSubfaults();
        
        if (rank == 0) {
            std::cout << "Generated " << subfaults.size() << " subfaults for rupture model\n\n";
        }
        
        // =====================================================================
        // 2. Setup Bathymetry
        // =====================================================================
        
        BathymetryGrid bathymetry;
        
        if (use_synthetic_bathy) {
            if (rank == 0) {
                std::cout << "Generating synthetic West Coast bathymetry...\n";
            }
            
            // Create synthetic bathymetry for the Cascadia margin
            bathymetry = TsunamiUtils::generateTestBathymetry(
                -132.0, -122.0,    // Longitude range
                39.0, 51.0,        // Latitude range
                1.0 / 30.0,        // ~3 km resolution
                200.0,             // Continental shelf depth
                3500.0             // Deep ocean depth
            );
            
            if (rank == 0) {
                std::cout << "  Grid: " << bathymetry.nlon << " x " << bathymetry.nlat 
                          << " points\n";
                std::cout << "  Resolution: " << bathymetry.dlon * 60.0 
                          << " arc-minutes\n\n";
            }
        } else {
            // Try to load real bathymetry from file
            if (rank == 0) {
                std::cout << "Loading bathymetry from GEBCO data...\n";
            }
            
            // This would load real GEBCO data if available
            if (!TsunamiUtils::loadGEBCO("data/GEBCO_2023_Cascadia.nc",
                                         -132.0, -122.0, 39.0, 51.0, bathymetry)) {
                if (rank == 0) {
                    std::cout << "  GEBCO file not found, using synthetic bathymetry\n";
                }
                bathymetry = TsunamiUtils::generateTestBathymetry(
                    -132.0, -122.0, 39.0, 51.0, 1.0/30.0, 200.0, 3500.0);
            }
        }
        
        // =====================================================================
        // 3. Configure Tsunami Simulation
        // =====================================================================
        
        TsunamiConfig config;
        config.lon_min = bathymetry.lon_min;
        config.lon_max = bathymetry.lon_max;
        config.lat_min = bathymetry.lat_min;
        config.lat_max = bathymetry.lat_max;
        config.base_resolution = bathymetry.dlon;
        config.end_time = end_time;
        config.dt = 1.0;  // 1 second time step
        config.cfl_number = 0.8;
        config.adaptive_timestep = true;
        config.solver = SWESolver::FWAVE;
        config.wetting_drying = WettingDryingMethod::MOMENTUM_CONSERVING;
        config.friction_model = BottomFrictionModel::MANNING;
        config.manning_n_ocean = 0.025;
        config.source_rise_time = 60.0;  // 1 minute source rise time
        config.output_interval = 120;
        config.save_max_values = true;
        config.save_inundation_map = true;
        
        if (rank == 0) {
            std::cout << "Tsunami simulation configuration:\n";
            std::cout << "  Domain: " << std::fixed << std::setprecision(1)
                      << config.lon_min << "° to " << config.lon_max << "° E, "
                      << config.lat_min << "° to " << config.lat_max << "° N\n";
            std::cout << "  Duration: " << config.end_time / 3600.0 << " hours\n";
            std::cout << "  Solver: " << "F-wave (GeoClaw-style)\n";
            std::cout << "\n";
        }
        
        // =====================================================================
        // 4. Setup Observation Gauges
        // =====================================================================
        
        auto gauges = setupGaugeNetwork();
        
        if (rank == 0) {
            std::cout << "Observation network: " << gauges.size() << " stations\n";
            for (const auto& g : gauges) {
                std::cout << "  • " << g.name << " (" << g.longitude << "°E, " 
                          << g.latitude << "°N)";
                if (g.is_dart_buoy) {
                    std::cout << " [DART, " << g.depth << "m depth]";
                }
                std::cout << "\n";
            }
            std::cout << "\n";
        }
        
        printExpectedArrivals(rank);
        
        // =====================================================================
        // 5. Run Coupled Earthquake-Tsunami Simulation
        // =====================================================================
        
        if (rank == 0) {
            std::cout << "═══════════════════════════════════════════════════════════════\n";
            std::cout << "                    STARTING SIMULATION                         \n";
            std::cout << "═══════════════════════════════════════════════════════════════\n\n";
        }
        
        CoupledTsunamiEarthquake simulation;
        simulation.initializeCascadia(fault_model, bathymetry, config);
        
        // Add gauges
        for (const auto& g : gauges) {
            // Note: gauges added internally in the solver
        }
        
        // Run the simulation
        auto sim_start = std::chrono::high_resolution_clock::now();
        
        simulation.run();
        
        auto sim_end = std::chrono::high_resolution_clock::now();
        auto sim_duration = std::chrono::duration_cast<std::chrono::seconds>(
            sim_end - sim_start).count();
        
        if (rank == 0) {
            std::cout << "\nSimulation completed in " << sim_duration << " seconds\n";
            std::cout << "  Simulated time: " << config.end_time << " seconds ("
                      << config.end_time / 3600.0 << " hours)\n";
            std::cout << "  Speed-up factor: " << config.end_time / sim_duration << "x\n";
        }
        
        // =====================================================================
        // 6. Output Results
        // =====================================================================
        
        if (rank == 0) {
            std::cout << "\nWriting output to: " << output_dir << "\n";
        }
        
        simulation.writeOutput(output_dir);
        
        // Print results summary
        printResultsSummary(simulation, rank);
        
        // =====================================================================
        // 7. Summary and Warnings
        // =====================================================================
        
        if (rank == 0) {
            std::cout << "\n";
            std::cout << "════════════════════════════════════════════════════════════════\n";
            std::cout << "                    SIMULATION COMPLETE                          \n";
            std::cout << "════════════════════════════════════════════════════════════════\n\n";
            
            std::cout << "Output files:\n";
            std::cout << "  • Gauge time series: " << output_dir << "/*_timeseries.dat\n";
            std::cout << "  • Maximum values: " << output_dir << "/maximum_values.dat\n";
            std::cout << "  • Inundation map: " << output_dir << "/inundation_extent.dat\n";
            std::cout << "\n";
            
            std::cout << "⚠️  IMPORTANT NOTES:\n";
            std::cout << "  • This is a HYPOTHETICAL scenario for research and planning\n";
            std::cout << "  • Results depend on assumed slip distribution and bathymetry\n";
            std::cout << "  • Actual tsunami impacts would vary based on:\n";
            std::cout << "      - Exact rupture location and slip distribution\n";
            std::cout << "      - Tidal state at time of earthquake\n";
            std::cout << "      - Local bathymetry and coastal geometry\n";
            std::cout << "      - Sedimentation and coastal infrastructure\n";
            std::cout << "\n";
            std::cout << "For official tsunami hazard information, see:\n";
            std::cout << "  • NOAA Center for Tsunami Research: https://nctr.pmel.noaa.gov/\n";
            std::cout << "  • Pacific Northwest Seismic Network: https://pnsn.org/\n";
            std::cout << "  • FEMA Region X Cascadia Rising: https://www.fema.gov/\n";
            std::cout << "\n";
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error on rank " << rank << ": " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    auto end_time_wall = std::chrono::high_resolution_clock::now();
    auto total_duration = std::chrono::duration_cast<std::chrono::seconds>(
        end_time_wall - start_time).count();
    
    if (rank == 0) {
        std::cout << "Total wall time: " << total_duration << " seconds\n\n";
    }
    
    MPI_Finalize();
    return 0;
}
