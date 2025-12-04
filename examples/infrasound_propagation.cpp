/**
 * @file infrasound_propagation.cpp
 * @brief Example: Infrasound propagation from an explosion source
 * 
 * This example demonstrates:
 * - Setting up atmospheric profiles with wind
 * - Configuring infrasound sources (explosions)
 * - Including topographic effects
 * - Computing propagation using ray tracing and parabolic equation
 * - Generating synthetic waveforms at receivers
 * - Unit conversions for input/output
 * 
 * All calculations are performed in SI base units (m, kg, s).
 * User can specify input values in any supported unit.
 * 
 * Usage:
 *   mpirun -np 4 ./infrasound_propagation [config_file]
 * 
 * Default config: config/infrasound_example.config
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <memory>
#include <cmath>

#include "AtmosphericInfrasound.hpp"
#include "UnitSystem.hpp"
#include "ConfigReader.hpp"

using namespace FSRM;

/**
 * @brief Print atmospheric profile summary
 */
void printAtmosphereProfile(const AtmosphericProfile& atm) {
    std::cout << "\n=== Atmospheric Profile ===" << std::endl;
    std::cout << std::setw(12) << "Altitude(km)" 
              << std::setw(12) << "Temp(°C)"
              << std::setw(12) << "c(m/s)"
              << std::setw(12) << "Wind E(m/s)"
              << std::setw(12) << "Wind N(m/s)"
              << std::setw(12) << "ceff(m/s)" << std::endl;
    std::cout << std::string(72, '-') << std::endl;
    
    UnitSystem& units = UnitSystemManager::getInstance();
    
    for (double z = 0; z <= 120000; z += 10000) {
        AtmosphericState state = atm.getState(z);
        double temp_C = state.temperature - 273.15;  // K to °C
        double z_km = units.fromBase(z, "km");
        double ceff = state.effectiveSoundSpeed(0.0);  // Eastward
        
        std::cout << std::fixed << std::setprecision(1)
                  << std::setw(12) << z_km
                  << std::setw(12) << temp_C
                  << std::setw(12) << state.sound_speed
                  << std::setw(12) << state.wind_east
                  << std::setw(12) << state.wind_north
                  << std::setw(12) << ceff << std::endl;
    }
}

/**
 * @brief Compute and print transmission loss
 */
void computeTransmissionLoss(
    ParabolicEquationSolver& pe,
    const std::vector<double>& frequencies,
    double receiver_range,
    double receiver_height
) {
    UnitSystem& units = UnitSystemManager::getInstance();
    
    std::cout << "\n=== Transmission Loss at " 
              << units.fromBase(receiver_range, "km") << " km ===" << std::endl;
    std::cout << std::setw(15) << "Frequency(Hz)" 
              << std::setw(15) << "TL(dB)" << std::endl;
    std::cout << std::string(30, '-') << std::endl;
    
    for (double f : frequencies) {
        pe.setFrequency(f);
        pe.propagate();
        double TL = pe.transmissionLoss(receiver_range, receiver_height);
        
        std::cout << std::fixed << std::setprecision(3)
                  << std::setw(15) << f
                  << std::setw(15) << TL << std::endl;
    }
}

/**
 * @brief Trace eigenrays and print summary
 */
void traceAndPrintRays(
    InfrasoundRayTracer& tracer,
    int max_bounces
) {
    std::cout << "\n=== Eigenray Analysis ===" << std::endl;
    
    auto rays = tracer.findEigenrays(max_bounces, 0.5);
    
    std::cout << "Found " << rays.size() << " eigenrays" << std::endl;
    std::cout << std::setw(8) << "Ray#"
              << std::setw(15) << "Travel Time(s)"
              << std::setw(15) << "Distance(km)"
              << std::setw(10) << "Bounces"
              << std::setw(15) << "Amplitude" << std::endl;
    std::cout << std::string(63, '-') << std::endl;
    
    UnitSystem& units = UnitSystemManager::getInstance();
    
    int i = 0;
    for (const auto& ray : rays) {
        if (ray.reached_receiver) {
            std::cout << std::fixed << std::setprecision(2)
                      << std::setw(8) << ++i
                      << std::setw(15) << ray.travel_time
                      << std::setw(15) << units.fromBase(ray.total_distance, "km")
                      << std::setw(10) << ray.num_bounces
                      << std::setw(15) << ray.final_amplitude << std::endl;
        }
    }
}

/**
 * @brief Main entry point
 */
int main(int argc, char* argv[]) {
    // Initialize PETSc (for future solver use)
    // PetscInitialize(&argc, &argv, nullptr, nullptr);
    
    std::cout << "========================================" << std::endl;
    std::cout << "  FSRM Infrasound Propagation Example" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "\nAll calculations in SI base units (m, kg, s)" << std::endl;
    std::cout << "Unit conversion handled automatically\n" << std::endl;
    
    // Get unit system
    UnitSystem& units = UnitSystemManager::getInstance();
    
    // =========================================================================
    // 1. Configure source parameters
    // =========================================================================
    std::cout << "=== Source Configuration ===" << std::endl;
    
    InfrasoundSourceParams source_params;
    source_params.type = InfrasoundSourceType::EXPLOSION;
    
    // Demonstrate unit conversion: specify yield in kilotons
    source_params.yield_kt = 1.0;  // 1 kt explosion
    
    // Source location (can specify in km, will convert to m internally)
    source_params.x = units.toBase(0.0, "km");
    source_params.y = units.toBase(0.0, "km");
    source_params.z = units.toBase(0.0, "m");  // Surface burst
    
    // Source timing
    source_params.onset_time = 0.0;
    
    // Estimate dominant frequency from yield (empirical relation)
    // T (period) ≈ 0.4 * W^(1/3) seconds, where W is yield in kt
    source_params.dominant_frequency = 1.0 / (0.4 * std::pow(source_params.yield_kt, 1.0/3.0));
    source_params.duration = 0.4 * std::pow(source_params.yield_kt, 1.0/3.0);
    
    std::cout << "  Type: Surface Explosion" << std::endl;
    std::cout << "  Yield: " << source_params.yield_kt << " kt TNT" << std::endl;
    std::cout << "  Location: (" << source_params.x << ", " 
              << source_params.y << ", " << source_params.z << ") m" << std::endl;
    std::cout << "  Dominant frequency: " << std::fixed << std::setprecision(2) 
              << source_params.dominant_frequency << " Hz" << std::endl;
    
    InfrasoundSource source(source_params);
    
    // =========================================================================
    // 2. Configure atmospheric profile
    // =========================================================================
    std::cout << "\n=== Atmospheric Configuration ===" << std::endl;
    
    auto atmosphere = std::make_shared<AtmosphericProfile>(AtmosphereModel::US_STANDARD_1976);
    
    // Add realistic wind profile (stratospheric jet)
    atmosphere->setWindProfile([](double z, double& u, double& v) {
        // Simple parametric wind model
        // Tropospheric westerlies peak around 10 km
        // Stratospheric polar vortex/jet around 50-60 km
        if (z < 12000) {
            // Troposphere: increasing westerly
            u = 20.0 * (z / 12000.0);
            v = 5.0 * std::sin(M_PI * z / 12000.0);
        } else if (z < 50000) {
            // Stratosphere: decreasing then increasing
            double factor = 1.0 - 2.0 * std::abs(z - 30000.0) / 40000.0;
            u = 20.0 + 30.0 * factor;
            v = 2.0;
        } else if (z < 90000) {
            // Upper stratosphere and mesosphere
            u = 20.0 * std::exp(-(z - 50000.0) / 20000.0);
            v = 0.0;
        } else {
            // Thermosphere
            u = 5.0;
            v = 0.0;
        }
    });
    
    // Print profile
    printAtmosphereProfile(*atmosphere);
    
    // =========================================================================
    // 3. Configure topography (flat terrain for this example)
    // =========================================================================
    std::cout << "\n=== Topography ===" << std::endl;
    
    auto topography = std::make_shared<TopographyModel>();
    topography->setFlatTerrain(0.0);  // Sea level
    
    std::cout << "  Model: Flat terrain at sea level" << std::endl;
    
    // For a realistic example with mountains, uncomment:
    // topography->setTerrainFunction([](double x, double y) {
    //     // Gaussian mountain 100 km east of source
    //     double x0 = 100000.0, y0 = 0.0;
    //     double sigma = 20000.0;
    //     double height = 3000.0;
    //     double r2 = (x - x0) * (x - x0) + (y - y0) * (y - y0);
    //     return height * std::exp(-r2 / (2.0 * sigma * sigma));
    // });
    
    // =========================================================================
    // 4. Configure receiver array
    // =========================================================================
    std::cout << "\n=== Receiver Array ===" << std::endl;
    
    std::vector<double> receiver_ranges;
    double receiver_azimuth = 0.0;  // Eastward (radians)
    double receiver_height = 1.0;    // 1 m above ground
    
    // Create array from 50 km to 500 km
    for (double r = 50000; r <= 500000; r += 50000) {
        receiver_ranges.push_back(r);
    }
    
    std::cout << "  " << receiver_ranges.size() << " receivers from " 
              << units.fromBase(receiver_ranges.front(), "km") << " km to "
              << units.fromBase(receiver_ranges.back(), "km") << " km" << std::endl;
    std::cout << "  Azimuth: " << units.fromBase(receiver_azimuth, "deg") << "°" << std::endl;
    std::cout << "  Height: " << receiver_height << " m above ground" << std::endl;
    
    // =========================================================================
    // 5. Ray tracing analysis
    // =========================================================================
    std::cout << "\n=== Ray Tracing ===" << std::endl;
    
    InfrasoundRayTracer tracer;
    tracer.setAtmosphere(atmosphere);
    tracer.setTopography(topography);
    tracer.setSource(source_params.x, source_params.y, source_params.z);
    tracer.setMaxDistance(units.toBase(600, "km"));
    tracer.setStepSize(100.0);  // 100 m step
    
    // Trace rays to receiver at 200 km
    double rx = units.toBase(200, "km");
    tracer.setReceiver(rx, 0.0, receiver_height);
    
    traceAndPrintRays(tracer, 5);
    
    // =========================================================================
    // 6. Parabolic equation propagation
    // =========================================================================
    std::cout << "\n=== Parabolic Equation Propagation ===" << std::endl;
    
    ParabolicEquationSolver pe;
    pe.setAtmosphere(atmosphere);
    pe.setTopography(topography);
    pe.setSource(source);
    pe.setGrid(50.0, 500.0, 120000.0, 600000.0);  // dz, dr, z_max, r_max
    
    // Compute transmission loss at different frequencies
    std::vector<double> frequencies = {0.1, 0.5, 1.0, 2.0, 5.0, 10.0};
    computeTransmissionLoss(pe, frequencies, rx, receiver_height);
    
    // =========================================================================
    // 7. Generate synthetic waveforms
    // =========================================================================
    std::cout << "\n=== Synthetic Waveforms ===" << std::endl;
    
    // Use 1 Hz for waveform generation
    pe.setFrequency(source_params.dominant_frequency);
    pe.propagate();
    
    std::vector<double> times, pressures;
    pe.getTimeSeries(rx, receiver_height, times, pressures);
    
    if (!times.empty()) {
        std::cout << "Waveform at " << units.fromBase(rx, "km") << " km:" << std::endl;
        std::cout << "  Samples: " << times.size() << std::endl;
        std::cout << "  Time range: " << times.front() << " - " 
                  << times.back() << " s" << std::endl;
        
        // Find peak amplitude
        double peak_pressure = 0.0;
        double peak_time = 0.0;
        for (size_t i = 0; i < pressures.size(); ++i) {
            if (std::abs(pressures[i]) > std::abs(peak_pressure)) {
                peak_pressure = pressures[i];
                peak_time = times[i];
            }
        }
        std::cout << "  Peak pressure: " << peak_pressure << " Pa at t = " 
                  << peak_time << " s" << std::endl;
        
        // Save waveform to file
        std::ofstream waveform_file("infrasound_waveform.txt");
        waveform_file << "# Infrasound waveform at " << units.fromBase(rx, "km") << " km\n";
        waveform_file << "# Time (s)    Pressure (Pa)\n";
        for (size_t i = 0; i < times.size(); ++i) {
            waveform_file << std::scientific << std::setprecision(6)
                         << times[i] << "  " << pressures[i] << "\n";
        }
        waveform_file.close();
        std::cout << "  Waveform saved to: infrasound_waveform.txt" << std::endl;
    }
    
    // =========================================================================
    // 8. Summary of propagation characteristics
    // =========================================================================
    std::cout << "\n=== Propagation Summary ===" << std::endl;
    
    // Check for atmospheric ducts
    bool has_strat_duct = atmosphere->hasAcousticDuct(0, 60000, 0);
    bool has_thermo_duct = atmosphere->hasAcousticDuct(80000, 120000, 0);
    
    std::cout << "Atmospheric ducting:" << std::endl;
    std::cout << "  Stratospheric duct (0-60 km): " 
              << (has_strat_duct ? "YES" : "NO") << std::endl;
    std::cout << "  Thermospheric duct (80-120 km): " 
              << (has_thermo_duct ? "YES" : "NO") << std::endl;
    
    // Estimate detection range for IMS-class station
    double noise_level = 0.01;  // Pa (typical microbarom noise)
    std::cout << "\nDetection estimates (SNR > 1):" << std::endl;
    std::cout << "  Noise floor: " << noise_level << " Pa" << std::endl;
    std::cout << "  For 1 kt explosion:" << std::endl;
    std::cout << "    Stratospheric returns: ~1000 km typical" << std::endl;
    std::cout << "    Thermospheric returns: ~500 km typical" << std::endl;
    
    // =========================================================================
    // 9. Unit conversion demonstration
    // =========================================================================
    std::cout << "\n=== Unit Conversion Examples ===" << std::endl;
    
    // Demonstrate the unit system
    double pressure_pa = 100.0;  // 100 Pa
    
    std::cout << "Pressure: " << pressure_pa << " Pa" << std::endl;
    std::cout << "  = " << units.fromBase(pressure_pa, "mbar") << " mbar" << std::endl;
    std::cout << "  = " << units.fromBase(pressure_pa, "psi") << " psi" << std::endl;
    
    double distance_m = 200000.0;  // 200 km
    std::cout << "\nDistance: " << distance_m << " m" << std::endl;
    std::cout << "  = " << units.fromBase(distance_m, "km") << " km" << std::endl;
    std::cout << "  = " << units.fromBase(distance_m, "mi") << " miles" << std::endl;
    
    // Parsing values with units
    double value_si = units.parseAndConvertToBase("500 km");
    std::cout << "\nParsed '500 km' = " << value_si << " m (SI)" << std::endl;
    
    value_si = units.parseAndConvertToBase("1.5 mbar");
    std::cout << "Parsed '1.5 mbar' = " << value_si << " Pa (SI)" << std::endl;
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "  Example completed successfully!" << std::endl;
    std::cout << "========================================" << std::endl;
    
    // Cleanup
    // PetscFinalize();
    
    return 0;
}
