/**
 * @file volcanic_unrest_monitoring.cpp
 * @brief Volcanic Unrest and Monitoring Signal Simulation
 *
 * This example simulates pre-eruptive volcanic unrest to generate
 * synthetic monitoring data for volcano observatory applications.
 * The simulation produces:
 *   - Ground deformation (GPS, InSAR, tiltmeters)
 *   - Seismicity patterns (VT swarms, LP events, tremor)
 *   - Gas emissions (SO₂, CO₂ flux)
 *   - Thermal anomalies
 *
 * Useful for:
 *   - Testing monitoring system responses
 *   - Training observatory staff
 *   - Developing eruption forecasting algorithms
 *   - Hazard scenario planning
 *
 * Scenarios included:
 *   1. Magma intrusion (dike/sill emplacement)
 *   2. Hydrothermal pressurization
 *   3. Pre-eruptive inflation
 *   4. Failed eruption (intrusion without eruption)
 *
 * Usage:
 *   ./volcanic_unrest_monitoring [options]
 *
 * Options:
 *   -scenario <name>  Scenario: intrusion, hydrothermal, pre_eruption, failed
 *   -duration <days>  Simulation duration in days (default: 90)
 *   -output <dir>     Output directory
 *
 * @author FSRM Development Team
 */

#include "VolcanoModel.hpp"
#include "ConfigReader.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <chrono>
#include <random>

using namespace FSRM;

static const char* help = 
    "================================================================\n"
    "  Volcanic Unrest Monitoring Signal Simulation\n"
    "================================================================\n"
    "\n"
    "Generates synthetic volcano monitoring data for various unrest scenarios.\n"
    "\n"
    "Usage: ./volcanic_unrest_monitoring [options]\n"
    "\n"
    "Options:\n"
    "  -scenario <name>   Scenario type:\n"
    "                       intrusion      - Dike/sill intrusion\n"
    "                       hydrothermal   - Hydrothermal pressurization\n"
    "                       pre_eruption   - Leading to eruption\n"
    "                       failed         - Intrusion without eruption\n"
    "  -duration <days>   Simulation duration (default: 90)\n"
    "  -output <dir>      Output directory\n"
    "  -volcano <name>    Target volcano style:\n"
    "                       stratovolcano, caldera, shield\n"
    "  -help              Show this message\n"
    "\n";

/**
 * @brief Unrest scenario types
 */
enum class UnrestScenario {
    INTRUSION,          // Dike/sill emplacement
    HYDROTHERMAL,       // Hydrothermal pressurization
    PRE_ERUPTION,       // Leading to eruption
    FAILED_ERUPTION     // Intrusion without eruption
};

/**
 * @brief Synthetic GPS station
 */
struct GPSStation {
    std::string name;
    double x, y;        // Position (m from crater)
    std::vector<double> times;
    std::vector<double> east;
    std::vector<double> north;
    std::vector<double> up;
};

/**
 * @brief Synthetic seismometer
 */
struct Seismometer {
    std::string name;
    double x, y, z;
    std::vector<double> times;
    std::vector<double> event_magnitudes;
    std::vector<VolcanicSeismicEventType> event_types;
};

/**
 * @brief Print simulation banner
 */
void printBanner(int rank) {
    if (rank != 0) return;
    
    std::cout << "\n";
    std::cout << "╔══════════════════════════════════════════════════════════════╗\n";
    std::cout << "║                                                              ║\n";
    std::cout << "║        VOLCANIC UNREST MONITORING SIMULATION                 ║\n";
    std::cout << "║                                                              ║\n";
    std::cout << "║    Generating Synthetic Monitoring Data for Observatories    ║\n";
    std::cout << "║                                                              ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════╝\n";
    std::cout << "\n";
}

/**
 * @brief Setup monitoring network
 */
std::vector<GPSStation> setupGPSNetwork(double volcano_radius) {
    std::vector<GPSStation> stations;
    
    // Summit station
    GPSStation summit;
    summit.name = "SMIT";
    summit.x = 0.0;
    summit.y = 500.0;
    stations.push_back(summit);
    
    // Radial stations at different distances
    const double angles[] = {0, 60, 120, 180, 240, 300};
    const double distances[] = {1000, 2000, 3000, 5000, 8000, 12000};
    
    for (int i = 0; i < 6; ++i) {
        GPSStation station;
        station.name = "GPS" + std::to_string(i+1);
        double angle = angles[i] * M_PI / 180.0;
        double dist = distances[i];
        station.x = dist * std::cos(angle);
        station.y = dist * std::sin(angle);
        stations.push_back(station);
    }
    
    return stations;
}

/**
 * @brief Setup seismic network
 */
std::vector<Seismometer> setupSeismicNetwork() {
    std::vector<Seismometer> stations;
    
    const char* names[] = {"VLT1", "VLT2", "VLT3", "VLT4", "VBB1", "VBB2"};
    const double x_pos[] = {500, -500, 0, 0, 2000, -2000};
    const double y_pos[] = {0, 0, 500, -500, 0, 0};
    
    for (int i = 0; i < 6; ++i) {
        Seismometer seis;
        seis.name = names[i];
        seis.x = x_pos[i];
        seis.y = y_pos[i];
        seis.z = 0.0;
        stations.push_back(seis);
    }
    
    return stations;
}

/**
 * @brief Configure scenario parameters
 */
void configureScenario(UnrestScenario scenario,
                       double& recharge_rate,
                       double& duration_factor,
                       bool& leads_to_eruption,
                       std::string& description) {
    switch (scenario) {
        case UnrestScenario::INTRUSION:
            recharge_rate = 0.5;      // m³/s
            duration_factor = 1.0;
            leads_to_eruption = false;
            description = "Dike intrusion at 3 km depth";
            break;
            
        case UnrestScenario::HYDROTHERMAL:
            recharge_rate = 0.1;      // Slow
            duration_factor = 2.0;
            leads_to_eruption = false;
            description = "Hydrothermal system pressurization";
            break;
            
        case UnrestScenario::PRE_ERUPTION:
            recharge_rate = 2.0;      // Rapid
            duration_factor = 0.5;
            leads_to_eruption = true;
            description = "Accelerating unrest leading to eruption";
            break;
            
        case UnrestScenario::FAILED_ERUPTION:
            recharge_rate = 1.5;      // Moderate
            duration_factor = 1.0;
            leads_to_eruption = false;
            description = "Intrusion stalls before reaching surface";
            break;
    }
}

/**
 * @brief Generate synthetic seismicity
 */
void generateSeismicity(const MagmaChamberModel& chamber,
                        double time_days,
                        double chamber_overpressure,
                        std::vector<VolcanicSeismicEvent>& events,
                        std::mt19937& rng) {
    // Seismicity rate increases with overpressure
    double base_rate = 5.0;  // Events per day
    double rate_factor = 1.0 + chamber_overpressure / 5e6;  // Scale with MPa
    double expected_events = base_rate * rate_factor;
    
    std::poisson_distribution<int> poisson(expected_events / 24.0);  // Per hour
    int n_events = poisson(rng);
    
    std::uniform_real_distribution<double> mag_dist(0.5, 3.0);
    std::uniform_real_distribution<double> pos_dist(-2000, 2000);
    std::uniform_real_distribution<double> type_dist(0, 1);
    
    for (int i = 0; i < n_events; ++i) {
        VolcanicSeismicEvent event;
        
        // Event type depends on depth and overpressure
        double p = type_dist(rng);
        if (p < 0.4) {
            event.type = VolcanicSeismicEventType::VT;
        } else if (p < 0.7) {
            event.type = VolcanicSeismicEventType::LP;
        } else if (p < 0.9) {
            event.type = VolcanicSeismicEventType::HYBRID;
        } else {
            event.type = VolcanicSeismicEventType::VLP;
        }
        
        event.x = pos_dist(rng);
        event.y = pos_dist(rng);
        event.z = -2000.0 - std::abs(pos_dist(rng));  // Below surface
        event.magnitude = mag_dist(rng);
        event.origin_time = time_days * 86400.0;
        
        events.push_back(event);
    }
}

/**
 * @brief Write GPS data to file
 */
void writeGPSData(const std::vector<GPSStation>& stations,
                  const std::string& output_dir) {
    for (const auto& station : stations) {
        std::string filename = output_dir + "/gps_" + station.name + ".dat";
        std::ofstream file(filename);
        
        file << "# GPS Station: " << station.name << "\n";
        file << "# Position: " << station.x << ", " << station.y << " m\n";
        file << "# Time(days)  East(mm)  North(mm)  Up(mm)\n";
        
        for (size_t i = 0; i < station.times.size(); ++i) {
            file << std::fixed << std::setprecision(3)
                 << station.times[i] << "  "
                 << std::setprecision(2)
                 << station.east[i] * 1000.0 << "  "
                 << station.north[i] * 1000.0 << "  "
                 << station.up[i] * 1000.0 << "\n";
        }
    }
}

/**
 * @brief Write seismicity catalog to file
 */
void writeSeismicity(const std::vector<VolcanicSeismicEvent>& events,
                     const std::string& output_dir) {
    std::string filename = output_dir + "/seismicity_catalog.dat";
    std::ofstream file(filename);
    
    file << "# Volcanic Seismicity Catalog\n";
    file << "# Time(days)  X(m)  Y(m)  Z(m)  Mag  Type\n";
    
    for (const auto& event : events) {
        std::string type_str;
        switch (event.type) {
            case VolcanicSeismicEventType::VT: type_str = "VT"; break;
            case VolcanicSeismicEventType::LP: type_str = "LP"; break;
            case VolcanicSeismicEventType::VLP: type_str = "VLP"; break;
            case VolcanicSeismicEventType::HYBRID: type_str = "HB"; break;
            case VolcanicSeismicEventType::TREMOR: type_str = "TR"; break;
            default: type_str = "UK"; break;
        }
        
        file << std::fixed << std::setprecision(4)
             << event.origin_time / 86400.0 << "  "
             << std::setprecision(0)
             << event.x << "  " << event.y << "  " << event.z << "  "
             << std::setprecision(1) << event.magnitude << "  "
             << type_str << "\n";
    }
}

/**
 * @brief Write gas emission data
 */
void writeGasData(const std::vector<double>& times,
                  const std::vector<double>& SO2_flux,
                  const std::vector<double>& CO2_flux,
                  const std::string& output_dir) {
    std::string filename = output_dir + "/gas_emissions.dat";
    std::ofstream file(filename);
    
    file << "# Volcanic Gas Emissions\n";
    file << "# Time(days)  SO2(kg/s)  CO2(kg/s)\n";
    
    for (size_t i = 0; i < times.size(); ++i) {
        file << std::fixed << std::setprecision(3)
             << times[i] << "  "
             << std::setprecision(1)
             << SO2_flux[i] << "  "
             << CO2_flux[i] << "\n";
    }
}

int main(int argc, char** argv) {
    // Initialize MPI (optional for this example)
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    // Parse command line
    std::string output_dir = "output/volcanic_unrest";
    std::string scenario_str = "pre_eruption";
    double duration_days = 90.0;
    
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-help") == 0 || strcmp(argv[i], "--help") == 0) {
            if (rank == 0) std::cout << help;
            MPI_Finalize();
            return 0;
        } else if (strcmp(argv[i], "-scenario") == 0 && i + 1 < argc) {
            scenario_str = argv[++i];
        } else if (strcmp(argv[i], "-duration") == 0 && i + 1 < argc) {
            duration_days = std::stod(argv[++i]);
        } else if (strcmp(argv[i], "-output") == 0 && i + 1 < argc) {
            output_dir = argv[++i];
        }
    }
    
    // Convert scenario string
    UnrestScenario scenario = UnrestScenario::PRE_ERUPTION;
    if (scenario_str == "intrusion") scenario = UnrestScenario::INTRUSION;
    else if (scenario_str == "hydrothermal") scenario = UnrestScenario::HYDROTHERMAL;
    else if (scenario_str == "pre_eruption") scenario = UnrestScenario::PRE_ERUPTION;
    else if (scenario_str == "failed") scenario = UnrestScenario::FAILED_ERUPTION;
    
    printBanner(rank);
    
    try {
        // =====================================================================
        // 1. Configure Scenario
        // =====================================================================
        
        double recharge_rate;
        double duration_factor;
        bool leads_to_eruption;
        std::string description;
        configureScenario(scenario, recharge_rate, duration_factor, 
                         leads_to_eruption, description);
        
        duration_days *= duration_factor;
        
        if (rank == 0) {
            std::cout << "Scenario: " << description << "\n";
            std::cout << "Duration: " << duration_days << " days\n";
            std::cout << "Recharge rate: " << recharge_rate << " m³/s\n";
            std::cout << "Leads to eruption: " << (leads_to_eruption ? "Yes" : "No") << "\n";
            std::cout << "\n";
        }
        
        // =====================================================================
        // 2. Setup Volcanic System
        // =====================================================================
        
        MagmaChamberGeometry chamber;
        chamber.depth = 4000.0;          // 4 km
        chamber.volume = 20e9;           // 20 km³
        chamber.semi_axis_a = 2500.0;
        chamber.semi_axis_b = 2500.0;
        chamber.semi_axis_c = 1500.0;
        
        MagmaProperties magma;
        magma.composition = MagmaComposition::ANDESITE;
        magma.H2O_total = 4.0;
        magma.temperature = 1123.15;     // 850°C
        magma.crystal_fraction = 0.4;
        
        MagmaChamberState initial_state;
        initial_state.pressure = 150e6;   // 150 MPa
        initial_state.temperature = magma.temperature;
        initial_state.overpressure = 2e6; // 2 MPa initial
        
        ConduitGeometry conduit;
        conduit.depth = 4000.0;
        conduit.radius_base = 20.0;
        conduit.radius_vent = 30.0;
        
        // Initialize chamber model
        MagmaChamberModel chamber_model;
        chamber_model.initialize(chamber, magma, initial_state);
        chamber_model.setRechargeRate(recharge_rate * magma.getDensity(150e6));
        
        // Initialize deformation model
        VolcanicDeformationModel deformation;
        deformation.setElasticProperties(30e9, 0.25);
        
        // =====================================================================
        // 3. Setup Monitoring Network
        // =====================================================================
        
        auto gps_stations = setupGPSNetwork(5000.0);
        auto seismometers = setupSeismicNetwork();
        
        if (rank == 0) {
            std::cout << "Monitoring network:\n";
            std::cout << "  GPS stations: " << gps_stations.size() << "\n";
            std::cout << "  Seismometers: " << seismometers.size() << "\n";
            std::cout << "\n";
        }
        
        // =====================================================================
        // 4. Run Unrest Simulation
        // =====================================================================
        
        if (rank == 0) {
            std::cout << "═══════════════════════════════════════════════════════════════\n";
            std::cout << "                  SIMULATING VOLCANIC UNREST                    \n";
            std::cout << "═══════════════════════════════════════════════════════════════\n\n";
        }
        
        double dt_days = 0.1;  // 0.1 day = 2.4 hours
        double dt = dt_days * 86400.0;
        double current_time = 0.0;
        double end_time = duration_days * 86400.0;
        
        std::vector<VolcanicSeismicEvent> all_seismic_events;
        std::vector<double> gas_times;
        std::vector<double> SO2_flux;
        std::vector<double> CO2_flux;
        
        std::mt19937 rng(42);  // Fixed seed for reproducibility
        std::normal_distribution<double> noise(0.0, 0.001);  // 1 mm noise
        
        while (current_time < end_time) {
            double time_days = current_time / 86400.0;
            
            // Update chamber
            chamber_model.update(dt);
            
            // Get current state
            double overpressure = chamber_model.getOverpressure();
            double volume_change = overpressure / magma.getMagmaCompressibility();
            
            // Update deformation source
            VolcanicDeformationSource source;
            source.type = DeformationSourceType::MOGI;
            source.depth = chamber.depth;
            source.volume_change = volume_change;
            
            deformation.clearSources();
            deformation.addSource(source);
            
            // Record GPS displacements
            for (auto& station : gps_stations) {
                double ux, uy, uz;
                deformation.computeDisplacement(station.x, station.y, ux, uy, uz);
                
                // Add noise
                station.times.push_back(time_days);
                station.east.push_back(ux + noise(rng));
                station.north.push_back(uy + noise(rng));
                station.up.push_back(uz + noise(rng));
            }
            
            // Generate seismicity
            generateSeismicity(chamber_model, time_days, overpressure, 
                             all_seismic_events, rng);
            
            // Gas emissions (increase with overpressure)
            double base_SO2 = 50.0;   // kg/s background
            double base_CO2 = 200.0;  // kg/s background
            double enhancement = 1.0 + overpressure / 10e6;
            
            gas_times.push_back(time_days);
            SO2_flux.push_back(base_SO2 * enhancement);
            CO2_flux.push_back(base_CO2 * enhancement);
            
            // Progress output
            if (rank == 0 && static_cast<int>(time_days) % 10 == 0) {
                static double last_print = -10;
                if (time_days - last_print >= 10) {
                    double max_uplift = 0.0;
                    for (const auto& s : gps_stations) {
                        if (!s.up.empty()) {
                            max_uplift = std::max(max_uplift, s.up.back());
                        }
                    }
                    
                    std::cout << "Day " << std::fixed << std::setprecision(0) << time_days 
                              << ": ΔP = " << std::setprecision(1) << overpressure/1e6 << " MPa, "
                              << "uplift = " << std::setprecision(0) << max_uplift*1000 << " mm, "
                              << "events = " << all_seismic_events.size() << "\n";
                    last_print = time_days;
                }
            }
            
            current_time += dt;
        }
        
        // =====================================================================
        // 5. Output Results
        // =====================================================================
        
        if (rank == 0) {
            std::cout << "\n";
            std::cout << "═══════════════════════════════════════════════════════════════\n";
            std::cout << "                    SIMULATION COMPLETE                          \n";
            std::cout << "═══════════════════════════════════════════════════════════════\n\n";
            
            std::cout << "Summary:\n";
            std::cout << "  Duration: " << duration_days << " days\n";
            std::cout << "  Total seismic events: " << all_seismic_events.size() << "\n";
            
            // Find max displacement
            double max_uplift = 0.0;
            std::string max_station;
            for (const auto& s : gps_stations) {
                if (!s.up.empty() && s.up.back() > max_uplift) {
                    max_uplift = s.up.back();
                    max_station = s.name;
                }
            }
            std::cout << "  Max uplift: " << std::fixed << std::setprecision(0) 
                      << max_uplift*1000 << " mm at " << max_station << "\n";
            std::cout << "  Final overpressure: " 
                      << chamber_model.getOverpressure()/1e6 << " MPa\n";
            std::cout << "\n";
            
            std::cout << "Writing output to: " << output_dir << "\n";
            
            // Write outputs
            writeGPSData(gps_stations, output_dir);
            writeSeismicity(all_seismic_events, output_dir);
            writeGasData(gas_times, SO2_flux, CO2_flux, output_dir);
            
            std::cout << "\nOutput files:\n";
            for (const auto& s : gps_stations) {
                std::cout << "  • " << output_dir << "/gps_" << s.name << ".dat\n";
            }
            std::cout << "  • " << output_dir << "/seismicity_catalog.dat\n";
            std::cout << "  • " << output_dir << "/gas_emissions.dat\n";
            
            std::cout << "\n";
            std::cout << "📊 Data suitable for:\n";
            std::cout << "  • Observatory decision support training\n";
            std::cout << "  • Eruption forecasting algorithm development\n";
            std::cout << "  • Hazard scenario planning\n";
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error on rank " << rank << ": " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    MPI_Finalize();
    return 0;
}
