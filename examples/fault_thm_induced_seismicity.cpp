/**
 * @file fault_thm_induced_seismicity.cpp
 * @brief Example: THM-coupled fault for induced seismicity simulation
 * 
 * This example demonstrates how to use the FaultCohesiveTHM class to simulate
 * induced seismicity during fluid injection, with full thermo-hydro-mechanical
 * (THM) coupling including:
 * 
 * - Pore pressure diffusion and its effect on effective stress
 * - Thermal pressurization during slip
 * - Slip-dependent permeability enhancement
 * - Custom property update callbacks
 * 
 * Physical scenario:
 * - Pre-existing fault at seismogenic depth (5 km)
 * - Fluid injection increases pore pressure
 * - Fault destabilizes when effective normal stress decreases
 * - Slip causes permeability enhancement (creating flow conduit)
 * - Frictional heating induces thermal pressurization
 * 
 * Usage:
 *   ./fault_thm_induced_seismicity [config_file]
 */

#include "PyLithFault.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <memory>
#include <vector>
#include <map>
#include <sstream>

using namespace FSRM;

// =============================================================================
// Custom Property Updater for Induced Seismicity
// =============================================================================

/**
 * @brief Custom property updater with injection-specific physics
 * 
 * This updater implements:
 * 1. Slip-dependent permeability enhancement (orders of magnitude increase)
 * 2. Stress-dependent permeability (fault closure under load)
 * 3. Dilatancy-induced porosity changes
 * 4. Thermal weakening of friction
 */
class InducedSeismicityPropertyUpdater : public FaultPropertyUpdateCallback {
public:
    // Parameters
    double slip_enhancement_coeff = 50.0;      // k increases rapidly with slip
    double max_permeability_ratio = 1e4;       // Max k/k0
    double stress_sensitivity = 0.2;           // Stress sensitivity
    double reference_stress = 50e6;            // Reference effective stress
    double dilatancy_coeff = 0.05;             // Porosity change per unit slip
    double thermal_weakening_temp = 500.0;     // Temperature for weakening (K above ambient)
    double ambient_temperature = 423.15;       // 150°C at 5 km depth
    
    FaultPermeabilityTensor updatePermeability(
        const FaultMultiPhysicsState& state,
        const FaultHydraulicProperties& /*props*/,
        double /*dt*/) override {
        
        FaultPermeabilityTensor k = state.current_permeability;
        
        // Slip enhancement: k = k0 * min(exp(α * slip), max_ratio)
        double slip_factor = std::exp(slip_enhancement_coeff * state.slip_total);
        slip_factor = std::min(slip_factor, max_permeability_ratio);
        
        // Stress effect: k = k * exp(-β * σ'_n / σ_ref)
        // Low effective stress (high pore pressure) increases permeability
        double stress_factor = 1.0;
        if (state.effective_normal_stress > 0) {
            stress_factor = std::exp(-stress_sensitivity * 
                          state.effective_normal_stress / reference_stress);
        } else {
            // Tensile effective stress - fault opening
            stress_factor = 10.0;  // Large increase
        }
        
        // Combined effect
        double total_factor = slip_factor * stress_factor;
        
        k.k_strike = k.k_strike_ref * total_factor;
        k.k_dip = k.k_dip_ref * total_factor;
        k.k_normal = k.k_normal_ref * total_factor;
        
        // Bounds
        k.k_strike = std::clamp(k.k_strike, 1e-20, 1e-10);
        k.k_dip = std::clamp(k.k_dip, 1e-20, 1e-10);
        k.k_normal = std::clamp(k.k_normal, 1e-20, 1e-13);
        
        return k;
    }
    
    double updatePorosity(
        const FaultMultiPhysicsState& state,
        const FaultHydraulicProperties& props,
        double /*dt*/) override {
        
        // Dilatancy: porosity increases with slip
        double phi = props.porosity_ref + dilatancy_coeff * state.slip_total;
        
        // Compaction under stress
        if (state.effective_normal_stress > 0) {
            double compaction = 0.01 * state.effective_normal_stress / reference_stress;
            phi -= compaction;
        }
        
        return std::clamp(phi, 0.01, 0.5);
    }
    
    double updateFriction(
        const FaultMultiPhysicsState& state,
        double base_friction) override {
        
        // Thermal weakening
        double T_excess = state.temperature - ambient_temperature;
        if (T_excess > 0) {
            // Exponential weakening
            double weakening = 1.0 - 0.8 * (1.0 - std::exp(-T_excess / thermal_weakening_temp));
            return base_friction * weakening;
        }
        
        return base_friction;
    }
};

// =============================================================================
// Configuration Parser
// =============================================================================

std::map<std::string, std::string> parseConfig(const std::string& filename) {
    std::map<std::string, std::string> config;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Warning: Could not open config file: " << filename << std::endl;
        return config;
    }
    
    std::string line;
    while (std::getline(file, line)) {
        // Skip comments and empty lines
        if (line.empty() || line[0] == '#' || line[0] == ';') continue;
        
        size_t eq_pos = line.find('=');
        if (eq_pos == std::string::npos) continue;
        
        std::string key = line.substr(0, eq_pos);
        std::string value = line.substr(eq_pos + 1);
        
        // Trim whitespace
        auto trim = [](std::string& s) {
            s.erase(0, s.find_first_not_of(" \t"));
            s.erase(s.find_last_not_of(" \t") + 1);
        };
        trim(key);
        trim(value);
        
        config[key] = value;
    }
    
    return config;
}

double getConfigDouble(const std::map<std::string, std::string>& config,
                      const std::string& key, double default_val) {
    auto it = config.find(key);
    if (it != config.end() && !it->second.empty()) {
        try { return std::stod(it->second); } catch (...) {}
    }
    return default_val;
}

std::string getConfigString(const std::map<std::string, std::string>& config,
                           const std::string& key, const std::string& default_val) {
    auto it = config.find(key);
    return (it != config.end() && !it->second.empty()) ? it->second : default_val;
}

// =============================================================================
// Output Functions
// =============================================================================

void writeTimeSeries(const std::string& filename,
                    const std::vector<double>& time,
                    const std::vector<double>& slip,
                    const std::vector<double>& slip_rate,
                    const std::vector<double>& pore_pressure,
                    const std::vector<double>& permeability,
                    const std::vector<double>& temperature) {
    std::ofstream out(filename);
    out << "# Time series output for THM fault\n";
    out << "# time(s)  slip(m)  slip_rate(m/s)  pore_pressure(MPa)  permeability(m^2)  temperature(K)\n";
    
    for (size_t i = 0; i < time.size(); ++i) {
        out << std::scientific << std::setprecision(6)
            << time[i] << "  "
            << slip[i] << "  "
            << slip_rate[i] << "  "
            << pore_pressure[i] / 1e6 << "  "
            << permeability[i] << "  "
            << temperature[i] << "\n";
    }
}

void writeSnapshot(const std::string& filename,
                  const FaultCohesiveTHM& fault,
                  double time) {
    std::ofstream out(filename);
    out << "# Fault snapshot at t = " << time << " s\n";
    out << "# x(m)  slip(m)  slip_rate(m/s)  shear_stress(MPa)  eff_normal(MPa)  "
        << "pore_pressure(MPa)  temperature(K)  permeability(m^2)  porosity\n";
    
    const auto& hyd_field = fault.getHydraulicField();
    (void)fault.getThermalField();  // Field available but not currently used in output
    
    for (size_t i = 0; i < hyd_field.size(); ++i) {
        const auto& state = fault.getMultiPhysicsState(i);
        auto k = fault.getPermeabilityAt(i);
        double phi = fault.getPorosityAt(i);
        
        out << std::scientific << std::setprecision(6)
            << (i * 100.0) << "  "  // x position (assume 100m spacing)
            << state.slip_total << "  "
            << state.slip_rate << "  "
            << state.shear_stress / 1e6 << "  "
            << state.effective_normal_stress / 1e6 << "  "
            << state.pore_pressure / 1e6 << "  "
            << state.temperature << "  "
            << k.k_strike << "  "
            << phi << "\n";
    }
}

void writeGnuplotScript(const std::string& base_name) {
    std::ofstream out(base_name + "_plot.gp");
    out << R"(# Gnuplot script for THM fault induced seismicity visualization
set terminal pngcairo enhanced size 1200,900 font 'Arial,10'

# Time series plots
set output ')" << base_name << R"(_timeseries.png'
set multiplot layout 3,2 title 'THM Fault: Induced Seismicity Time Series' font ',14'

# Slip
set xlabel 'Time (s)'
set ylabel 'Slip (m)'
set grid
plot ')" << base_name << R"(_timeseries.dat' u 1:2 w l lw 2 title 'Slip'

# Slip rate
set ylabel 'Slip Rate (m/s)'
set logscale y
plot ')" << base_name << R"(_timeseries.dat' u 1:($3>1e-15?$3:1e-15) w l lw 2 title 'Slip Rate'
unset logscale y

# Pore pressure
set ylabel 'Pore Pressure (MPa)'
plot ')" << base_name << R"(_timeseries.dat' u 1:4 w l lw 2 title 'Pore Pressure'

# Permeability
set ylabel 'Permeability (m^2)'
set logscale y
plot ')" << base_name << R"(_timeseries.dat' u 1:5 w l lw 2 title 'Permeability'
unset logscale y

# Temperature
set ylabel 'Temperature (K)'
plot ')" << base_name << R"(_timeseries.dat' u 1:6 w l lw 2 title 'Temperature'

unset multiplot

# Snapshot plots
set output ')" << base_name << R"(_snapshot.png'
set multiplot layout 2,2 title 'THM Fault: Spatial Distribution (Final State)' font ',14'

# Slip profile
set xlabel 'Distance along fault (m)'
set ylabel 'Slip (m)'
plot ')" << base_name << R"(_snapshot.dat' u 1:2 w l lw 2 title 'Slip'

# Stress profile
set ylabel 'Stress (MPa)'
plot ')" << base_name << R"(_snapshot.dat' u 1:4 w l lw 2 title 'Shear Stress', \
     ')" << base_name << R"(_snapshot.dat' u 1:5 w l lw 2 title 'Effective Normal Stress'

# Pore pressure and temperature
set ylabel 'Pore Pressure (MPa)'
set y2label 'Temperature (K)'
set y2tics
plot ')" << base_name << R"(_snapshot.dat' u 1:6 w l lw 2 title 'Pore Pressure' axes x1y1, \
     ')" << base_name << R"(_snapshot.dat' u 1:7 w l lw 2 title 'Temperature' axes x1y2

# Permeability and porosity
set ylabel 'Permeability (m^2)'
set y2label 'Porosity'
set logscale y
plot ')" << base_name << R"(_snapshot.dat' u 1:8 w l lw 2 title 'Permeability' axes x1y1, \
     ')" << base_name << R"(_snapshot.dat' u 1:9 w l lw 2 title 'Porosity' axes x1y2
unset logscale y

unset multiplot
print 'Plots generated: )" << base_name << R"(_timeseries.png, )" << base_name << R"(_snapshot.png'
)";
    out.close();
    std::cout << "Gnuplot script written to: " << base_name << "_plot.gp" << std::endl;
}

// =============================================================================
// Main Simulation
// =============================================================================

int main(int argc, char* argv[]) {
    std::cout << "========================================" << std::endl;
    std::cout << "  THM Fault: Induced Seismicity Example" << std::endl;
    std::cout << "========================================\n" << std::endl;
    
    // Load configuration
    std::string config_file = (argc > 1) ? argv[1] : "config/fault_thm_induced_seismicity.config";
    auto config = parseConfig(config_file);
    
    // Simulation parameters
    double fault_length = getConfigDouble(config, "fault_length", 5000.0);      // 5 km
    double fault_depth = getConfigDouble(config, "fault_depth", 5000.0);        // 5 km depth
    int num_vertices = static_cast<int>(getConfigDouble(config, "num_vertices", 51));
    double end_time = getConfigDouble(config, "end_time", 100.0);               // 100 s
    double dt = getConfigDouble(config, "dt", 0.01);                            // 10 ms time step
    int output_interval = static_cast<int>(getConfigDouble(config, "output_interval", 100));
    
    // Material properties
    double shear_modulus = getConfigDouble(config, "shear_modulus", 30e9);      // 30 GPa
    double bulk_density = getConfigDouble(config, "density", 2700.0);           // 2700 kg/m³
    
    // Initial stress state (at 5 km depth)
    double initial_normal_stress = getConfigDouble(config, "initial_normal_stress", 120e6);  // 120 MPa
    double initial_shear_stress = getConfigDouble(config, "initial_shear_stress", 60e6);     // 60 MPa
    double initial_pore_pressure = getConfigDouble(config, "initial_pore_pressure", 50e6);   // 50 MPa hydrostatic
    
    // Injection parameters
    double injection_rate = getConfigDouble(config, "injection_pressure_rate", 1e6);  // 1 MPa/s pressure increase
    (void)getConfigDouble(config, "injection_radius", 500.0);     // 500 m influence radius - available for future use
    double injection_center = fault_length / 2.0;  // Center of fault
    
    std::cout << "Configuration:" << std::endl;
    std::cout << "  Fault length: " << fault_length/1000 << " km" << std::endl;
    std::cout << "  Fault depth: " << fault_depth/1000 << " km" << std::endl;
    std::cout << "  Initial stress: σn = " << initial_normal_stress/1e6 
              << " MPa, τ = " << initial_shear_stress/1e6 << " MPa" << std::endl;
    std::cout << "  Initial pore pressure: " << initial_pore_pressure/1e6 << " MPa" << std::endl;
    std::cout << "  Injection rate: " << injection_rate/1e6 << " MPa/s" << std::endl;
    std::cout << std::endl;
    
    // Create THM fault with configuration
    auto fault = std::make_unique<FaultCohesiveTHM>();
    
    // Merge simulation config with fault config
    config["fault_coupling_mode"] = "thm";
    config["friction_model"] = getConfigString(config, "friction_model", "rate_state_aging");
    config["friction_coefficient"] = getConfigString(config, "friction_coefficient", "0.6");
    config["reference_velocity"] = getConfigString(config, "reference_velocity", "1e-6");
    config["characteristic_slip"] = getConfigString(config, "characteristic_slip", "0.01");
    config["a_parameter"] = getConfigString(config, "a_parameter", "0.01");
    config["b_parameter"] = getConfigString(config, "b_parameter", "0.015");
    config["fault_porosity"] = getConfigString(config, "fault_porosity", "0.1");
    config["fault_permeability_parallel"] = getConfigString(config, "fault_permeability_parallel", "1e-15");
    config["fault_permeability_perpendicular"] = getConfigString(config, "fault_permeability_perpendicular", "1e-17");
    config["fault_thickness"] = getConfigString(config, "fault_thickness", "0.01");
    config["fault_biot_coefficient"] = getConfigString(config, "fault_biot_coefficient", "0.9");
    config["permeability_update_model"] = getConfigString(config, "permeability_update_model", "combined");
    config["fault_reference_temperature"] = std::to_string(423.15);  // 150°C at depth
    
    fault->configure(config);
    
    // Set up custom property updater
    auto custom_updater = std::make_shared<InducedSeismicityPropertyUpdater>();
    fault->setPropertyUpdater(custom_updater);
    
    // Create fault vertices (along-strike direction)
    std::vector<FaultVertex> vertices(num_vertices);
    double dx = fault_length / (num_vertices - 1);
    for (int i = 0; i < num_vertices; ++i) {
        vertices[i].vertex_id = i;
        vertices[i].coords = {i * dx, 0.0, -fault_depth};
        vertices[i].area = dx * 1.0;  // 1 m width
    }
    fault->setFaultVertices(vertices);
    
    // Set initial stress state
    DynamicFaultState initial_state;
    initial_state.effective_normal = initial_normal_stress;
    initial_state.shear_stress = initial_shear_stress;
    
    std::vector<DynamicFaultState> states(num_vertices, initial_state);
    fault->setInitialState(states);
    
    // Initialize the fault
    fault->initialize();
    
    // Set initial pore pressure field
    std::vector<double> pore_pressure(num_vertices, initial_pore_pressure);
    fault->setPorePressureField(pore_pressure);
    
    std::cout << "Fault initialized with " << num_vertices << " vertices" << std::endl;
    std::cout << "Initial effective normal stress: " 
              << (initial_normal_stress - 0.9 * initial_pore_pressure) / 1e6 << " MPa" << std::endl;
    std::cout << "Friction coefficient: " << 0.6 << std::endl;
    std::cout << "Critical shear stress: " 
              << 0.6 * (initial_normal_stress - 0.9 * initial_pore_pressure) / 1e6 << " MPa" << std::endl;
    std::cout << std::endl;
    
    // Storage for time series output (at injection point)
    std::vector<double> ts_time, ts_slip, ts_slip_rate, ts_pore_pressure, ts_permeability, ts_temperature;
    int center_vertex = num_vertices / 2;
    
    // Time-stepping loop
    std::cout << "Starting simulation..." << std::endl;
    double time = 0.0;
    int step = 0;
    
    bool event_occurred = false;
    double event_time = 0.0;
    double max_slip_rate = 0.0;
    
    while (time < end_time) {
        // Update pore pressure from injection
        // Diffusion-like pressure increase centered at injection point
        for (int i = 0; i < num_vertices; ++i) {
            double x = i * dx;
            double dist = std::abs(x - injection_center);
            
            // Gaussian pressure perturbation spreading from injection
            double diffusivity = 0.1;  // m²/s
            double diff_length = std::sqrt(4.0 * diffusivity * (time + 1.0));
            double pressure_increase = injection_rate * time * 
                                      std::exp(-dist * dist / (diff_length * diff_length));
            
            pore_pressure[i] = initial_pore_pressure + pressure_increase;
        }
        fault->setPorePressureField(pore_pressure);
        
        // Pre-step: update effective stress
        fault->prestep(time, dt);
        
        // Compute residual and update (simplified - in real code would use Newton solver)
        // Here we just let the fault update based on current state
        fault->computeResidual(nullptr, nullptr);
        
        // Post-step: update properties
        fault->poststep(time, dt);
        
        // Store time series data
        if (step % output_interval == 0) {
            const auto& state = fault->getMultiPhysicsState(center_vertex);
            auto k = fault->getPermeabilityAt(center_vertex);
            
            ts_time.push_back(time);
            ts_slip.push_back(state.slip_total);
            ts_slip_rate.push_back(state.slip_rate);
            ts_pore_pressure.push_back(state.pore_pressure);
            ts_permeability.push_back(k.k_strike);
            ts_temperature.push_back(state.temperature);
            
            // Track maximum slip rate
            if (state.slip_rate > max_slip_rate) {
                max_slip_rate = state.slip_rate;
            }
            
            // Detect seismic event (slip rate > 0.1 m/s)
            if (!event_occurred && state.slip_rate > 0.1) {
                event_occurred = true;
                event_time = time;
                std::cout << "  Seismic event detected at t = " << time << " s" << std::endl;
                std::cout << "    Pore pressure: " << state.pore_pressure/1e6 << " MPa" << std::endl;
                std::cout << "    Effective normal stress: " 
                          << state.effective_normal_stress/1e6 << " MPa" << std::endl;
            }
        }
        
        // Progress output
        if (step % (output_interval * 10) == 0) {
            std::cout << "  t = " << std::setw(8) << std::fixed << std::setprecision(2) 
                      << time << " s, p_max = " << std::setprecision(1)
                      << (*std::max_element(pore_pressure.begin(), pore_pressure.end()))/1e6 
                      << " MPa" << std::endl;
        }
        
        time += dt;
        ++step;
    }
    
    std::cout << "\nSimulation completed." << std::endl;
    std::cout << "  Total time: " << end_time << " s" << std::endl;
    std::cout << "  Total steps: " << step << std::endl;
    std::cout << "  Max slip rate: " << max_slip_rate << " m/s" << std::endl;
    if (event_occurred) {
        std::cout << "  Event time: " << event_time << " s" << std::endl;
        
        // Compute moment magnitude
        const auto& final_state = fault->getMultiPhysicsState(center_vertex);
        double total_slip = final_state.slip_total;
        double fault_area = fault_length * 1.0;  // 1 m width for 2D
        double Mw = FaultUtils::computeMomentMagnitude(total_slip, fault_area, shear_modulus);
        std::cout << "  Final slip: " << total_slip << " m" << std::endl;
        std::cout << "  Moment magnitude: Mw " << std::fixed << std::setprecision(1) << Mw << std::endl;
    }
    
    // Write output
    std::string base_name = "fault_thm_induced";
    writeTimeSeries(base_name + "_timeseries.dat", ts_time, ts_slip, ts_slip_rate,
                   ts_pore_pressure, ts_permeability, ts_temperature);
    writeSnapshot(base_name + "_snapshot.dat", *fault, time);
    writeGnuplotScript(base_name);
    
    std::cout << "\nOutput files written:" << std::endl;
    std::cout << "  " << base_name << "_timeseries.dat" << std::endl;
    std::cout << "  " << base_name << "_snapshot.dat" << std::endl;
    std::cout << "  " << base_name << "_plot.gp" << std::endl;
    
    std::cout << "\nTo visualize, run: gnuplot " << base_name << "_plot.gp" << std::endl;
    
    return 0;
}
