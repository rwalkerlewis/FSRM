/**
 * @file pylith_strike_slip_example.cpp
 * @brief PyLith-style fault modeling example: Strike-slip earthquake
 * 
 * This example demonstrates how to use the PyLith-style fault representation
 * to model a strike-slip earthquake with:
 * 
 * 1. KINEMATIC FAULT: Prescribed slip history (circular rupture propagation)
 *    - Brune slip time function
 *    - Rupture propagates from hypocenter
 *    - Computes slip, slip rate, and tractions over time
 * 
 * 2. DYNAMIC FAULT: Spontaneous rupture with friction
 *    - Slip-weakening or rate-state friction
 *    - Traction perturbation to nucleate rupture
 *    - Self-consistent rupture propagation
 * 
 * Reference:
 *   - Aagaard et al. (2013): PyLith domain decomposition approach
 *   - Harris et al. (2018): SCEC/USGS dynamic rupture benchmarks
 * 
 * Usage:
 *   ./pylith_strike_slip_example [options]
 * 
 * Options:
 *   -type <kinematic|dynamic>  Fault type (default: kinematic)
 *   -friction <sw|rs>          Friction model for dynamic (default: sw)
 *   -output <dir>              Output directory (default: output)
 *   -time <seconds>            Simulation end time (default: 10.0)
 */

#include "PyLithFault.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <sys/stat.h>

using namespace FSRM;

// =============================================================================
// Configuration
// =============================================================================

struct SimConfig {
    // Fault type
    bool is_kinematic = true;
    std::string friction_type = "slip_weakening";
    
    // Domain
    double fault_length = 20000.0;   // 20 km
    double fault_width = 15000.0;    // 15 km
    double depth_top = 0.0;          // Surface rupturing
    
    // Grid
    int nx = 41;   // Along strike
    int nz = 31;   // Down dip
    
    // Material properties
    double shear_modulus = 30.0e9;   // 30 GPa
    double poisson_ratio = 0.25;
    double density = 2700.0;         // kg/m³
    double shear_velocity = 3464.0;  // √(μ/ρ) m/s
    
    // Kinematic rupture parameters
    double final_slip = 2.0;         // 2 m slip
    double rise_time = 1.0;          // 1 s rise time
    double rupture_velocity = 2800.0; // 0.8 * Vs
    double hypo_x = 0.0;             // Hypocenter along strike
    double hypo_z = -10000.0;        // Hypocenter depth (10 km)
    double rake = 180.0;             // Pure left-lateral (degrees)
    
    // Dynamic rupture parameters
    double initial_shear = 70.0e6;   // 70 MPa initial shear stress
    double initial_normal = -120.0e6; // 120 MPa normal stress (compression)
    double static_friction = 0.677;  // Static friction coefficient
    double dynamic_friction = 0.525; // Dynamic friction coefficient
    double critical_slip = 0.4;      // D_c = 0.4 m
    
    // Perturbation to nucleate
    double pert_amplitude = 11.6e6;  // 11.6 MPa
    double pert_duration = 1.0;      // 1 s
    double pert_radius = 3000.0;     // 3 km nucleation zone
    
    // Time stepping
    double dt = 0.01;                // 10 ms time step
    double end_time = 10.0;          // 10 s simulation
    
    // Output
    std::string output_dir = "output";
    int output_interval = 10;        // Every 10 steps
};

// =============================================================================
// Helper Functions
// =============================================================================

void createOutputDirectory(const std::string& dir) {
    #ifdef _WIN32
    _mkdir(dir.c_str());
    #else
    mkdir(dir.c_str(), 0755);
    #endif
}

void parseCommandLine(int argc, char** argv, SimConfig& config) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        
        if (arg == "-type" && i + 1 < argc) {
            std::string type = argv[++i];
            config.is_kinematic = (type == "kinematic" || type == "kin");
        } else if (arg == "-friction" && i + 1 < argc) {
            config.friction_type = argv[++i];
        } else if (arg == "-output" && i + 1 < argc) {
            config.output_dir = argv[++i];
        } else if (arg == "-time" && i + 1 < argc) {
            config.end_time = std::stod(argv[++i]);
        } else if (arg == "-slip" && i + 1 < argc) {
            config.final_slip = std::stod(argv[++i]);
        } else if (arg == "-h" || arg == "--help") {
            std::cout << "PyLith-Style Strike-Slip Earthquake Example\n\n";
            std::cout << "Usage: " << argv[0] << " [options]\n\n";
            std::cout << "Options:\n";
            std::cout << "  -type <kinematic|dynamic>  Fault type (default: kinematic)\n";
            std::cout << "  -friction <sw|rs>          Friction for dynamic (default: sw)\n";
            std::cout << "  -output <dir>              Output directory (default: output)\n";
            std::cout << "  -time <seconds>            End time (default: 10.0)\n";
            std::cout << "  -slip <meters>             Final slip for kinematic (default: 2.0)\n";
            std::cout << "  -h, --help                 Show this help\n";
            exit(0);
        }
    }
}

// =============================================================================
// Create Fault Geometry
// =============================================================================

std::vector<FaultVertex> createFaultMesh(const SimConfig& config) {
    std::vector<FaultVertex> vertices;
    
    double dx = config.fault_length / (config.nx - 1);
    double dz = config.fault_width / (config.nz - 1);
    
    // Strike direction (N-S fault striking N)
    std::array<double, 3> strike = {0.0, 1.0, 0.0};
    // Down-dip direction (vertical fault)
    std::array<double, 3> dip = {0.0, 0.0, -1.0};
    // Normal (east-facing)
    std::array<double, 3> normal = {1.0, 0.0, 0.0};
    
    int vertex_id = 0;
    for (int j = 0; j < config.nz; ++j) {
        for (int i = 0; i < config.nx; ++i) {
            FaultVertex v;
            v.vertex_id = vertex_id;
            v.vertex_negative = vertex_id;  // Will be assigned by mesh
            v.vertex_positive = vertex_id + config.nx * config.nz;
            v.vertex_lagrange = vertex_id + 2 * config.nx * config.nz;
            
            // Coordinates (centered on hypocenter along-strike)
            v.coords[0] = 0.0;  // On fault plane (x = 0)
            v.coords[1] = -config.fault_length / 2.0 + i * dx;  // Along strike
            v.coords[2] = config.depth_top - j * dz;  // Depth (negative down)
            
            v.normal = normal;
            v.along_strike = strike;
            v.up_dip = dip;
            v.area = dx * dz;  // Cell area
            
            vertices.push_back(v);
            ++vertex_id;
        }
    }
    
    return vertices;
}

// =============================================================================
// Run Kinematic Simulation
// =============================================================================

void runKinematicSimulation(const SimConfig& config) {
    std::cout << "\n=== Kinematic Fault Simulation ===" << std::endl;
    std::cout << "Fault dimensions: " << config.fault_length/1000 << " km × " 
              << config.fault_width/1000 << " km" << std::endl;
    std::cout << "Final slip: " << config.final_slip << " m" << std::endl;
    std::cout << "Rise time: " << config.rise_time << " s" << std::endl;
    std::cout << "Rupture velocity: " << config.rupture_velocity << " m/s" << std::endl;
    
    // Create fault
    FaultCohesiveKin fault;
    fault.setName("strike_slip_kinematic");
    
    // Set vertices
    auto vertices = createFaultMesh(config);
    fault.setFaultVertices(vertices);
    
    // Set circular rupture
    double rake_rad = config.rake * M_PI / 180.0;
    fault.setCircularRupture(
        0.0, config.hypo_x, config.hypo_z,  // Hypocenter (x, y, z)
        config.rupture_velocity,
        config.final_slip,
        config.rise_time,
        rake_rad
    );
    
    // Use Brune slip time function
    fault.setSlipTimeFnType(SlipTimeFnType::BRUNE);
    
    // Initialize
    fault.initialize();
    
    // Compute expected magnitude
    double area = config.fault_length * config.fault_width;
    double Mw = FaultUtils::computeMomentMagnitude(config.final_slip, area, config.shear_modulus);
    std::cout << "Expected magnitude: Mw " << std::fixed << std::setprecision(1) << Mw << std::endl;
    
    // Time stepping
    std::cout << "\nRunning simulation..." << std::endl;
    
    int num_steps = static_cast<int>(config.end_time / config.dt);
    
    // Output file for time series at hypocenter
    std::ofstream ts_file(config.output_dir + "/kinematic_time_series.dat");
    ts_file << "# Time series at hypocenter\n";
    ts_file << "# time(s) slip_ll(m) slip_rate(m/s)\n";
    
    // Find hypocenter vertex
    size_t hypo_idx = 0;
    double min_dist = 1e30;
    for (size_t i = 0; i < vertices.size(); ++i) {
        double dy = vertices[i].coords[1] - config.hypo_x;
        double dz = vertices[i].coords[2] - config.hypo_z;
        double dist = std::sqrt(dy*dy + dz*dz);
        if (dist < min_dist) {
            min_dist = dist;
            hypo_idx = i;
        }
    }
    
    for (int step = 0; step <= num_steps; ++step) {
        double t = step * config.dt;
        
        fault.prestep(t, config.dt);
        
        // Get slip at hypocenter
        double slip_ll, slip_rev, slip_open;
        double rate_ll, rate_rev, rate_open;
        fault.getSlipAtTime(t, hypo_idx, slip_ll, slip_rev, slip_open);
        fault.getSlipRateAtTime(t, hypo_idx, rate_ll, rate_rev, rate_open);
        
        ts_file << std::fixed << std::setprecision(4)
                << t << " " << slip_ll << " " << rate_ll << "\n";
        
        // Write snapshot
        if (step % config.output_interval == 0) {
            std::stringstream ss;
            ss << config.output_dir << "/kinematic_slip_" 
               << std::setfill('0') << std::setw(4) << step << ".dat";
            fault.writeOutput(ss.str(), t);
        }
        
        fault.poststep(t, config.dt);
        
        if (step % 100 == 0) {
            std::cout << "  t = " << std::fixed << std::setprecision(2) << t 
                      << " s, max slip = " << slip_ll << " m" << std::endl;
        }
    }
    
    ts_file.close();
    
    // Write final slip distribution
    std::ofstream final_file(config.output_dir + "/kinematic_final_slip.dat");
    final_file << "# Final slip distribution\n";
    final_file << "# y(m) z(m) slip_ll(m) slip_rev(m) slip_time(s)\n";
    
    for (size_t i = 0; i < vertices.size(); ++i) {
        double slip_ll, slip_rev, slip_open;
        fault.getSlipAtTime(config.end_time, i, slip_ll, slip_rev, slip_open);
        const auto& params = fault.getSlipParams(i);
        
        final_file << vertices[i].coords[1] << " " << vertices[i].coords[2] << " "
                   << slip_ll << " " << slip_rev << " " << params.slip_time << "\n";
    }
    final_file.close();
    
    std::cout << "\nKinematic simulation complete." << std::endl;
    std::cout << "Output written to: " << config.output_dir << "/" << std::endl;
}

// =============================================================================
// Run Dynamic Simulation
// =============================================================================

void runDynamicSimulation(const SimConfig& config) {
    std::cout << "\n=== Dynamic Fault Simulation ===" << std::endl;
    std::cout << "Friction model: " << config.friction_type << std::endl;
    std::cout << "Initial shear stress: " << config.initial_shear/1e6 << " MPa" << std::endl;
    std::cout << "Initial normal stress: " << -config.initial_normal/1e6 << " MPa" << std::endl;
    
    // Create fault
    FaultCohesiveDyn fault;
    fault.setName("strike_slip_dynamic");
    
    // Set vertices
    auto vertices = createFaultMesh(config);
    fault.setFaultVertices(vertices);
    
    // Set initial traction
    fault.setUniformInitialTraction(config.initial_shear, 0.0, config.initial_normal);
    
    // Set friction model
    if (config.friction_type == "rs" || config.friction_type == "rate_state") {
        auto rs = std::make_unique<RateStateFrictionAging>();
        rs->setDirectEffect(0.01);
        rs->setEvolutionEffect(0.015);
        rs->setCriticalSlipDistance(0.01);
        rs->setReferenceFriction(0.6);
        rs->setReferenceSlipRate(1e-6);
        fault.setFrictionModel(std::move(rs));
        std::cout << "Using rate-and-state friction (aging law)" << std::endl;
    } else {
        auto sw = std::make_unique<SlipWeakeningFriction>();
        sw->setStaticCoefficient(config.static_friction);
        sw->setDynamicCoefficient(config.dynamic_friction);
        sw->setCriticalSlipDistance(config.critical_slip);
        fault.setFrictionModel(std::move(sw));
        std::cout << "Using slip-weakening friction" << std::endl;
        std::cout << "  μ_s = " << config.static_friction << ", μ_d = " << config.dynamic_friction << std::endl;
        std::cout << "  D_c = " << config.critical_slip << " m" << std::endl;
    }
    
    // Set nucleation perturbation (Gaussian shear stress increase at hypocenter)
    fault.setTractionPerturbation(
        TractionPerturbationType::GAUSSIAN,
        config.pert_amplitude,
        config.pert_duration,
        0.0, config.hypo_x, config.hypo_z,  // Center at hypocenter
        config.pert_radius
    );
    std::cout << "Nucleation perturbation: " << config.pert_amplitude/1e6 << " MPa over " 
              << config.pert_duration << " s" << std::endl;
    
    // Initialize
    fault.initialize();
    
    // Compute stress ratio
    double tau0 = std::abs(config.initial_shear);
    double sigma_n = std::abs(config.initial_normal);
    double strength = config.static_friction * sigma_n;
    std::cout << "Initial stress ratio τ₀/τ_strength = " << tau0/strength << std::endl;
    
    // Time stepping
    std::cout << "\nRunning simulation..." << std::endl;
    
    int num_steps = static_cast<int>(config.end_time / config.dt);
    
    // Output file for time series
    std::ofstream ts_file(config.output_dir + "/dynamic_time_series.dat");
    ts_file << "# Dynamic rupture time series\n";
    ts_file << "# time(s) max_slip(m) max_slip_rate(m/s) ruptured(0/1)\n";
    
    double first_rupture_time = -1.0;
    
    for (int step = 0; step <= num_steps; ++step) {
        double t = step * config.dt;
        
        fault.prestep(t, config.dt);
        
        // Mock residual/Jacobian (actual implementation would couple to FE solver)
        std::vector<double> solution(vertices.size() * 6, 0.0);
        std::vector<double> residual(vertices.size() * 6, 0.0);
        fault.computeResidual(solution.data(), residual.data());
        
        fault.poststep(t, config.dt);
        
        // Track rupture progress
        double max_slip = fault.getMaxSlip();
        double max_rate = fault.getMaxSlipRate();
        bool has_ruptured = fault.hasRuptured();
        
        if (has_ruptured && first_rupture_time < 0) {
            first_rupture_time = fault.getRuptureTime();
            std::cout << "  Rupture initiated at t = " << first_rupture_time << " s" << std::endl;
        }
        
        ts_file << std::fixed << std::setprecision(4)
                << t << " " << max_slip << " " << max_rate << " " << has_ruptured << "\n";
        
        // Write snapshot
        if (step % config.output_interval == 0) {
            std::stringstream ss;
            ss << config.output_dir << "/dynamic_slip_" 
               << std::setfill('0') << std::setw(4) << step << ".dat";
            fault.writeOutput(ss.str(), t);
        }
        
        if (step % 100 == 0) {
            std::cout << "  t = " << std::fixed << std::setprecision(2) << t 
                      << " s, max slip = " << max_slip << " m, max rate = " << max_rate << " m/s" << std::endl;
        }
    }
    
    ts_file.close();
    
    // Summary
    std::cout << "\nDynamic simulation complete." << std::endl;
    if (first_rupture_time > 0) {
        std::cout << "Rupture initiated at: " << first_rupture_time << " s" << std::endl;
    } else {
        std::cout << "No rupture occurred (stress too low or perturbation insufficient)" << std::endl;
    }
    std::cout << "Final max slip: " << fault.getMaxSlip() << " m" << std::endl;
    std::cout << "Frictional work: " << fault.getFrictionalWork() << " J" << std::endl;
    std::cout << "Output written to: " << config.output_dir << "/" << std::endl;
}

// =============================================================================
// Main
// =============================================================================

int main(int argc, char** argv) {
    std::cout << "============================================================" << std::endl;
    std::cout << "  PyLith-Style Strike-Slip Earthquake Example" << std::endl;
    std::cout << "============================================================" << std::endl;
    
    SimConfig config;
    parseCommandLine(argc, argv, config);
    
    // Create output directory
    createOutputDirectory(config.output_dir);
    
    // Print configuration
    std::cout << "\nConfiguration:" << std::endl;
    std::cout << "  Fault type: " << (config.is_kinematic ? "kinematic" : "dynamic") << std::endl;
    std::cout << "  End time: " << config.end_time << " s" << std::endl;
    std::cout << "  Time step: " << config.dt << " s" << std::endl;
    std::cout << "  Output: " << config.output_dir << std::endl;
    
    // Run appropriate simulation
    if (config.is_kinematic) {
        runKinematicSimulation(config);
    } else {
        runDynamicSimulation(config);
    }
    
    // Write gnuplot script for visualization
    std::ofstream gp(config.output_dir + "/plot_results.gp");
    gp << "# Gnuplot script for visualizing results\n";
    gp << "set terminal pngcairo size 1200,600 enhanced font 'Arial,12'\n";
    gp << "\n";
    gp << "# Time series plot\n";
    gp << "set output 'time_series.png'\n";
    gp << "set xlabel 'Time (s)'\n";
    gp << "set ylabel 'Slip (m)' textcolor 'blue'\n";
    gp << "set y2label 'Slip rate (m/s)' textcolor 'red'\n";
    gp << "set y2tics\n";
    
    if (config.is_kinematic) {
        gp << "plot 'kinematic_time_series.dat' u 1:2 w l lw 2 lc 'blue' title 'Slip', \\\n";
        gp << "     '' u 1:3 w l lw 2 lc 'red' axes x1y2 title 'Slip rate'\n";
    } else {
        gp << "plot 'dynamic_time_series.dat' u 1:2 w l lw 2 lc 'blue' title 'Max Slip', \\\n";
        gp << "     '' u 1:3 w l lw 2 lc 'red' axes x1y2 title 'Max Slip Rate'\n";
    }
    gp.close();
    
    std::cout << "\n============================================================" << std::endl;
    std::cout << "  Example Complete" << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << "\nTo visualize results:" << std::endl;
    std::cout << "  cd " << config.output_dir << std::endl;
    std::cout << "  gnuplot plot_results.gp" << std::endl;
    
    return 0;
}
