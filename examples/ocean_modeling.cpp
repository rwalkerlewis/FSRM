/**
 * @file ocean_modeling.cpp
 * @brief Example: Comprehensive ocean hydrodynamics and physics modeling
 *
 * This example demonstrates the ocean physics and hydrodynamic modeling
 * capabilities of FSRM including:
 *
 * 1. Ocean Circulation Modeling
 *    - 3D primitive equations
 *    - Wind-driven currents
 *    - Thermohaline dynamics
 *    - Tidal forcing
 *
 * 2. Coastal Hydrodynamics
 *    - 2D shallow water equations
 *    - Storm surge prediction
 *    - Wetting/drying (inundation)
 *
 * 3. Wave Modeling
 *    - Spectral wave model
 *    - Wave-current interaction
 *    - Radiation stress
 *
 * 4. Estuarine Dynamics
 *    - Salt intrusion
 *    - River-ocean mixing
 *
 * @author FSRM Development Team
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <memory>

#include "FSRM.hpp"
#include "HydrodynamicModel.hpp"
#include "OceanPhysics.hpp"
#include "TsunamiModel.hpp"
#include "ConfigReader.hpp"

using namespace FSRM;

// =============================================================================
// Example 1: Wind-Driven Ocean Circulation
// =============================================================================

void runOceanCirculationExample() {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "Example 1: Wind-Driven Ocean Circulation\n";
    std::cout << std::string(70, '=') << "\n\n";
    
    // Configure the hydrodynamic model
    HydrodynamicConfig config;
    
    // Domain: Idealized basin (1000 km x 1000 km)
    config.model_type = HydrodynamicModelType::PRIMITIVE_3D;
    config.nx = 100;
    config.ny = 100;
    config.nz = 20;
    config.lon_min = 0.0;
    config.lon_max = 10.0;  // ~1000 km at mid-latitudes
    config.lat_min = 30.0;
    config.lat_max = 40.0;
    
    // Vertical configuration
    config.vertical_type = VerticalCoordinateType::SIGMA;
    config.num_sigma_levels = 20;
    config.theta_s = 5.0;
    config.theta_b = 0.4;
    
    // Time stepping
    config.time_scheme = TimeSteppingScheme::SPLIT_EXPLICIT;
    config.dt = 600.0;           // 10 minute baroclinic step
    config.barotropic_steps = 30;
    config.cfl_number = 0.8;
    config.adaptive_timestep = true;
    config.end_time = 86400.0 * 10;  // 10 days
    
    // Physical parameters
    config.g = 9.81;
    config.rho_0 = 1025.0;
    config.use_coriolis = true;
    config.use_beta_plane = true;
    config.latitude_reference = 35.0;
    config.f0 = 2.0 * 7.29e-5 * std::sin(35.0 * M_PI / 180.0);
    config.beta = 2.0e-11;
    
    // Mixing
    config.turbulence_closure = TurbulenceClosureType::SMAGORINSKY;
    config.horizontal_viscosity = 100.0;
    config.horizontal_diffusivity = 50.0;
    config.vertical_viscosity = 1e-4;
    config.vertical_diffusivity = 1e-5;
    
    // Bottom friction
    config.friction_type = BottomFrictionType::QUADRATIC;
    config.bottom_drag_coefficient = 0.0025;
    
    // Equation of state
    config.use_baroclinic = true;
    config.use_nonlinear_eos = true;
    
    // Output
    config.output_interval = 86400;  // Daily output
    config.output_dir = "output/ocean_circulation";
    
    // Initialize model
    HydrodynamicModel model;
    model.initialize(config);
    
    // Set bathymetry (flat bottom basin with coastal shelves)
    std::vector<double> depth(config.nx * config.ny);
    for (int j = 0; j < config.ny; ++j) {
        for (int i = 0; i < config.nx; ++i) {
            // Flat 4000m depth with sloping shelves
            double x_frac = static_cast<double>(i) / (config.nx - 1);
            double y_frac = static_cast<double>(j) / (config.ny - 1);
            
            double base_depth = 4000.0;
            
            // Western shelf
            if (x_frac < 0.1) {
                base_depth = 100.0 + 3900.0 * (x_frac / 0.1);
            }
            // Eastern shelf
            else if (x_frac > 0.9) {
                base_depth = 100.0 + 3900.0 * ((1.0 - x_frac) / 0.1);
            }
            
            depth[i + j * config.nx] = base_depth;
        }
    }
    model.setBathymetry(depth);
    
    // Set initial temperature and salinity (stratified)
    std::vector<double> T_init(config.nx * config.ny * config.nz);
    std::vector<double> S_init(config.nx * config.ny * config.nz);
    
    for (int k = 0; k < config.nz; ++k) {
        double z_frac = static_cast<double>(k) / (config.nz - 1);  // 0 = bottom, 1 = surface
        
        // Temperature profile (cold at bottom, warm at surface)
        double T = 4.0 + 14.0 * z_frac;  // 4°C at bottom, 18°C at surface
        
        // Add thermocline
        double thermocline_depth = 0.7;  // 70% depth from bottom
        if (z_frac > thermocline_depth) {
            double above_frac = (z_frac - thermocline_depth) / (1.0 - thermocline_depth);
            T = 10.0 + 8.0 * above_frac;  // Sharp increase through thermocline
        }
        
        // Salinity profile (slight increase with depth)
        double S = 35.0 + 0.5 * (1.0 - z_frac);
        
        for (int j = 0; j < config.ny; ++j) {
            for (int i = 0; i < config.nx; ++i) {
                int idx = i + j * config.nx + k * config.nx * config.ny;
                T_init[idx] = T;
                S_init[idx] = S;
            }
        }
    }
    model.setInitialTemperature(T_init);
    model.setInitialSalinity(S_init);
    
    // Set wind forcing (sinusoidal wind stress - classic double-gyre pattern)
    WindForcing wind;
    wind.u10.resize(config.nx * config.ny);
    wind.v10.resize(config.nx * config.ny);
    wind.use_bulk_formula = true;
    
    for (int j = 0; j < config.ny; ++j) {
        double y_frac = static_cast<double>(j) / (config.ny - 1);
        
        // Westerly wind with sinusoidal meridional variation
        // Creates subtropical and subpolar gyres
        double u10 = -0.1 * std::cos(M_PI * y_frac);  // Eastward at edges, westward in middle
        
        for (int i = 0; i < config.nx; ++i) {
            int idx = i + j * config.nx;
            wind.u10[idx] = u10 * 10.0;  // Scale to ~10 m/s
            wind.v10[idx] = 0.0;
        }
    }
    model.setWindForcing(wind);
    
    // Set open boundary conditions (simple radiation)
    OpenBoundaryCondition bc_west, bc_east, bc_south, bc_north;
    bc_west.location = OpenBoundaryCondition::WEST;
    bc_west.type = OpenBoundaryCondition::RADIATION;
    bc_east.location = OpenBoundaryCondition::EAST;
    bc_east.type = OpenBoundaryCondition::RADIATION;
    bc_south.location = OpenBoundaryCondition::SOUTH;
    bc_south.type = OpenBoundaryCondition::RADIATION;
    bc_north.location = OpenBoundaryCondition::NORTH;
    bc_north.type = OpenBoundaryCondition::RADIATION;
    
    model.addOpenBoundary(bc_west);
    model.addOpenBoundary(bc_east);
    model.addOpenBoundary(bc_south);
    model.addOpenBoundary(bc_north);
    
    // Run simulation
    std::cout << "Running ocean circulation simulation...\n";
    std::cout << "Domain: " << config.nx << " x " << config.ny << " x " << config.nz << "\n";
    std::cout << "Duration: " << config.end_time / 86400.0 << " days\n\n";
    
    double output_time = 0.0;
    int output_count = 0;
    
    model.run(config.end_time, [&](double t, const HydrodynamicModel& m) {
        double KE = m.computeKineticEnergy();
        
        const auto& eta = m.getElevation();
        double max_eta = *std::max_element(eta.begin(), eta.end());
        double min_eta = *std::min_element(eta.begin(), eta.end());
        
        std::cout << "Day " << std::fixed << std::setprecision(1) << t / 86400.0 
                  << ": η = [" << std::setprecision(3) << min_eta 
                  << ", " << max_eta << "] m, KE = " 
                  << std::scientific << KE / 1e15 << " PJ\n";
        
        output_count++;
    });
    
    std::cout << "\nOcean circulation simulation completed.\n";
    
    // Compute transport
    double transport = model.computeVolumeTransport(config.nx/2, 0, config.nx/2, config.ny-1);
    std::cout << "Mid-basin meridional transport: " << transport / 1e6 << " Sv\n";
}

// =============================================================================
// Example 2: Storm Surge Simulation
// =============================================================================

void runStormSurgeExample() {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "Example 2: Storm Surge and Coastal Flooding\n";
    std::cout << std::string(70, '=') << "\n\n";
    
    // Configure coastal hydrodynamics model
    HydrodynamicConfig config;
    
    // Domain: Idealized bay with coastal plain
    config.model_type = HydrodynamicModelType::SHALLOW_WATER_2D;
    config.nx = 200;
    config.ny = 100;
    config.nz = 1;
    config.lon_min = -1.0;
    config.lon_max = 1.0;   // ~200 km
    config.lat_min = 29.0;
    config.lat_max = 30.0;  // ~100 km
    
    // Time stepping (short for wetting/drying)
    config.time_scheme = TimeSteppingScheme::RK3;
    config.dt = 2.0;
    config.cfl_number = 0.7;
    config.adaptive_timestep = true;
    config.end_time = 86400.0 * 2;  // 2 days
    
    // Physical parameters
    config.use_coriolis = true;
    config.latitude_reference = 29.5;
    
    // Friction (higher for marshes/land)
    config.friction_type = BottomFrictionType::MANNING;
    config.manning_n = 0.025;
    
    // Wetting/drying
    config.enable_wetting_drying = true;
    config.dry_depth = 0.05;
    config.min_depth = 0.1;
    
    // Output
    config.output_interval = 900;  // Every 15 minutes
    config.output_dir = "output/storm_surge";
    
    // Initialize coastal model
    CoastalHydrodynamicsModel model;
    model.initializeStormSurge(config);
    
    // Set bathymetry/topography
    // Create idealized bay with continental shelf and coastal plain
    std::vector<double> depth(config.nx * config.ny);
    for (int j = 0; j < config.ny; ++j) {
        double y_frac = static_cast<double>(j) / (config.ny - 1);
        
        for (int i = 0; i < config.nx; ++i) {
            double x_frac = static_cast<double>(i) / (config.nx - 1);
            
            // Base profile: deep ocean to south, shallow/land to north
            double base;
            if (y_frac < 0.3) {
                // Continental shelf
                base = 50.0 + 150.0 * (0.3 - y_frac) / 0.3;
            } else if (y_frac < 0.5) {
                // Shallow coastal waters / bay
                base = 5.0 + 45.0 * (0.5 - y_frac) / 0.2;
            } else if (y_frac < 0.7) {
                // Low-lying coastal plain
                base = -2.0 + 7.0 * (0.7 - y_frac) / 0.2;  // -2 to 5 m elevation
            } else {
                // Higher ground
                base = -5.0 - 10.0 * (y_frac - 0.7) / 0.3;  // Negative = above water
            }
            
            // Add bay entrance in center
            if (y_frac > 0.3 && y_frac < 0.6 && std::abs(x_frac - 0.5) < 0.1) {
                base = std::max(base, 10.0);  // Channel through barrier
            }
            
            depth[i + j * config.nx] = base;
        }
    }
    model.setBathymetry(depth);
    
    // Set up hurricane track (idealized landfall)
    // Storm approaches from south, makes landfall, moves inland
    std::vector<double> track_time = {0, 12, 24, 36, 48};  // hours
    std::vector<double> track_lon = {0.0, 0.0, 0.0, 0.0, 0.0};
    std::vector<double> track_lat = {27.0, 28.0, 29.5, 30.5, 31.5};
    std::vector<double> pressure = {96000, 95000, 96000, 98000, 100000};  // Pa
    std::vector<double> max_wind = {45, 50, 45, 35, 20};  // m/s
    std::vector<double> rmw = {40, 35, 40, 50, 60};  // km
    
    // Convert hours to seconds
    for (auto& t : track_time) t *= 3600.0;
    
    model.setHurricaneTrack(track_time, track_lon, track_lat, pressure, max_wind, rmw);
    
    // Set tidal boundary conditions
    TidalForcing tides;
    tides.setupPrimaryConstituents();
    // Scale amplitudes for Gulf of Mexico
    for (auto& c : tides.constituents) {
        c.amplitude *= 0.3;  // ~0.3m tidal range
    }
    model.setTidalForcing(tides);
    
    // Set open boundaries
    OpenBoundaryCondition bc_south;
    bc_south.location = OpenBoundaryCondition::SOUTH;
    bc_south.type = OpenBoundaryCondition::FLATHER;
    bc_south.tidal = tides;
    bc_south.sponge_width = 20000.0;
    model.addOpenBoundary(bc_south);
    
    OpenBoundaryCondition bc_west, bc_east;
    bc_west.location = OpenBoundaryCondition::WEST;
    bc_west.type = OpenBoundaryCondition::SPONGE;
    bc_west.sponge_width = 10000.0;
    bc_east.location = OpenBoundaryCondition::EAST;
    bc_east.type = OpenBoundaryCondition::SPONGE;
    bc_east.sponge_width = 10000.0;
    model.addOpenBoundary(bc_west);
    model.addOpenBoundary(bc_east);
    
    // Run simulation
    std::cout << "Running storm surge simulation...\n";
    std::cout << "Domain: " << config.nx << " x " << config.ny << "\n";
    std::cout << "Duration: " << config.end_time / 3600.0 << " hours\n";
    std::cout << "Hurricane: Cat 2-3, making landfall at t=24h\n\n";
    
    double t = 0.0;
    double dt_output = 3600.0;  // Hourly progress
    double next_output = 0.0;
    
    while (t < config.end_time) {
        // Update wind and pressure from hurricane
        std::vector<double> u10, v10, p_atm;
        model.computeHollandWinds(t, u10, v10);
        model.computeHollandPressure(t, p_atm);
        
        // Set forcing
        WindForcing wind;
        wind.u10 = u10;
        wind.v10 = v10;
        model.setWindForcing(wind);
        
        AtmosphericPressureForcing pressure_forcing;
        pressure_forcing.pressure = p_atm;
        model.setAtmosphericPressure(pressure_forcing);
        
        // Step model
        double dt = model.step();
        t += dt;
        
        // Progress output
        if (t >= next_output) {
            const auto& eta = model.getElevation();
            double max_surge = *std::max_element(eta.begin(), eta.end());
            
            // Count inundated cells
            int inundated = 0;
            for (int j = config.ny / 2; j < config.ny; ++j) {
                for (int i = 0; i < config.nx; ++i) {
                    int idx = i + j * config.nx;
                    if (depth[idx] < 0 && eta[idx] > std::abs(depth[idx]) + 0.1) {
                        inundated++;
                    }
                }
            }
            
            std::cout << "Hour " << std::fixed << std::setprecision(0) << t / 3600.0
                      << ": Max surge = " << std::setprecision(2) << max_surge 
                      << " m, Inundated cells = " << inundated << "\n";
            
            next_output += dt_output;
        }
    }
    
    // Get maximum inundation
    std::vector<double> max_eta, max_depth, max_vel;
    model.computeMaxInundation(max_eta, max_depth, max_vel);
    
    double peak_surge = *std::max_element(max_eta.begin(), max_eta.end());
    double peak_depth = *std::max_element(max_depth.begin(), max_depth.end());
    double peak_vel = *std::max_element(max_vel.begin(), max_vel.end());
    
    std::cout << "\nStorm surge simulation completed.\n";
    std::cout << "Peak surge: " << peak_surge << " m\n";
    std::cout << "Max inundation depth: " << peak_depth << " m\n";
    std::cout << "Max current velocity: " << peak_vel << " m/s\n";
}

// =============================================================================
// Example 3: Spectral Wave Modeling
// =============================================================================

void runWaveModelExample() {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "Example 3: Spectral Wave Model\n";
    std::cout << std::string(70, '=') << "\n\n";
    
    // Configure wave model
    WaveModelConfig config;
    
    // Domain
    config.nx = 100;
    config.ny = 100;
    config.lon_min = -130.0;
    config.lon_max = -120.0;
    config.lat_min = 35.0;
    config.lat_max = 45.0;
    
    // Spectral discretization
    config.num_frequencies = 25;
    config.num_directions = 36;
    config.freq_min = 0.04;
    config.freq_max = 0.5;
    
    // Physics
    config.enable_wind_input = true;
    config.enable_whitecapping = true;
    config.enable_quadruplets = false;  // Simplify for example
    config.enable_depth_breaking = true;
    config.enable_bottom_friction = true;
    config.enable_refraction = true;
    
    // Breaking
    config.breaking_model = WaveBreakingModel::BATTJES_JANSSEN;
    config.breaking_gamma = 0.73;
    
    // Time stepping
    config.dt = 300.0;  // 5 minutes
    config.end_time = 86400.0;  // 1 day
    config.output_interval = 3600;
    
    // Initialize wave model
    WaveModel model;
    model.initialize(config);
    
    // Set bathymetry (continental shelf)
    std::vector<double> depth(config.nx * config.ny);
    for (int j = 0; j < config.ny; ++j) {
        for (int i = 0; i < config.nx; ++i) {
            double x_frac = static_cast<double>(i) / (config.nx - 1);
            
            // Depth profile: deep ocean to shallow coast
            double d;
            if (x_frac < 0.7) {
                d = 3000.0;  // Deep ocean
            } else if (x_frac < 0.9) {
                d = 3000.0 - 2800.0 * (x_frac - 0.7) / 0.2;  // Continental slope
            } else {
                d = 200.0 - 180.0 * (x_frac - 0.9) / 0.1;  // Shelf
            }
            
            depth[i + j * config.nx] = d;
        }
    }
    model.setBathymetry(depth);
    
    // Set wind forcing (uniform westerly wind)
    std::vector<double> u10(config.nx * config.ny, 15.0);  // 15 m/s
    std::vector<double> v10(config.nx * config.ny, 0.0);
    model.setWind(u10, v10);
    
    // Set initial wave conditions (small background swell)
    model.setInitialParameters(0.5, 10.0, 270.0);  // 0.5m, 10s, from west
    
    // Set boundary spectrum (incoming swell from west)
    auto boundary_spectrum = WaveSpectrum2D::generateJONSWAP(
        3.0, 12.0, 3.3, 270.0, 30.0,
        config.num_frequencies, config.num_directions);
    model.setBoundarySpectrum(0, boundary_spectrum);  // West boundary
    
    // Run simulation
    std::cout << "Running spectral wave model...\n";
    std::cout << "Domain: " << config.nx << " x " << config.ny << "\n";
    std::cout << "Spectrum: " << config.num_frequencies << " frequencies x " 
              << config.num_directions << " directions\n";
    std::cout << "Wind: 15 m/s from west\n\n";
    
    double t = 0.0;
    while (t < config.end_time) {
        double dt = model.step();
        t += dt;
        
        if (static_cast<int>(t) % 3600 == 0) {
            const auto& Hs = model.getHs();
            const auto& Tp = model.getTp();
            const auto& dir = model.getMeanDirection();
            
            // Find values at mid-domain
            int mid_idx = config.nx / 2 + (config.ny / 2) * config.nx;
            
            std::cout << "Hour " << t / 3600.0 
                      << ": Hs = " << std::fixed << std::setprecision(2) << Hs[mid_idx]
                      << " m, Tp = " << Tp[mid_idx]
                      << " s, Dir = " << dir[mid_idx] << "°\n";
        }
    }
    
    // Output final state
    const auto& Hs = model.getHs();
    double max_Hs = *std::max_element(Hs.begin(), Hs.end());
    
    std::cout << "\nWave model simulation completed.\n";
    std::cout << "Maximum significant wave height: " << max_Hs << " m\n";
    
    // Get radiation stress for circulation coupling
    std::vector<double> Sxx, Sxy, Syy;
    model.getRadiationStress(Sxx, Sxy, Syy);
    
    double max_Sxx = *std::max_element(Sxx.begin(), Sxx.end());
    std::cout << "Maximum radiation stress Sxx: " << max_Sxx / 1000.0 << " kN/m\n";
}

// =============================================================================
// Example 4: Estuarine Dynamics
// =============================================================================

void runEstuarineExample() {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "Example 4: Estuarine Salt Intrusion\n";
    std::cout << std::string(70, '=') << "\n\n";
    
    // Configure estuarine model
    HydrodynamicConfig config;
    
    // Domain: Idealized estuary
    config.model_type = HydrodynamicModelType::PRIMITIVE_3D;
    config.nx = 200;   // Along-estuary
    config.ny = 20;    // Cross-estuary
    config.nz = 10;    // Vertical
    config.lon_min = 0.0;
    config.lon_max = 0.5;   // ~50 km long
    config.lat_min = 40.0;
    config.lat_max = 40.05; // ~5 km wide
    
    // Vertical
    config.vertical_type = VerticalCoordinateType::SIGMA;
    config.num_sigma_levels = 10;
    
    // Time stepping
    config.time_scheme = TimeSteppingScheme::SPLIT_EXPLICIT;
    config.dt = 60.0;
    config.barotropic_steps = 10;
    config.end_time = 86400.0 * 2;  // 2 tidal days
    
    // Physics
    config.use_coriolis = false;  // Small estuary, ignore Coriolis
    config.use_baroclinic = true;
    config.use_nonlinear_eos = false;  // Linear for simplicity
    
    // Mixing
    config.turbulence_closure = TurbulenceClosureType::MELLOR_YAMADA_25;
    config.horizontal_viscosity = 10.0;
    config.vertical_viscosity = 1e-3;
    
    // Friction
    config.friction_type = BottomFrictionType::QUADRATIC;
    config.bottom_drag_coefficient = 0.003;
    
    // Output
    config.output_interval = 1800;  // Every 30 minutes
    
    // Initialize estuarine model
    EstuarineModel model;
    model.initializeEstuary(config, 5000.0, 10.0, 35.0);  // 5km wide, 10m deep, 35 PSU ocean
    
    // Set bathymetry (channel deepening toward mouth)
    std::vector<double> depth(config.nx * config.ny);
    for (int j = 0; j < config.ny; ++j) {
        double y_frac = static_cast<double>(j) / (config.ny - 1);
        double cross_depth = 1.0 - 4.0 * std::pow(y_frac - 0.5, 2);  // Parabolic cross-section
        
        for (int i = 0; i < config.nx; ++i) {
            double x_frac = static_cast<double>(i) / (config.nx - 1);
            
            // Depth increases toward mouth (east)
            double along_depth = 5.0 + 10.0 * x_frac;
            
            depth[i + j * config.nx] = along_depth * cross_depth;
            
            // Land at edges
            if (cross_depth < 0.3) {
                depth[i + j * config.nx] = -2.0;  // Land
            }
        }
    }
    model.setBathymetry(depth);
    
    // Set river head boundary
    model.setRiverHeadBoundary(500.0, 15.0);  // 500 m³/s, 15°C
    
    // Set ocean boundary with tides
    TidalForcing tides;
    tides.setupM2Only();
    tides.constituents[0].amplitude = 1.0;  // 1m M2 tide
    model.setOceanBoundary(tides, 35.0, 12.0);  // 35 PSU, 12°C
    
    // Run simulation
    std::cout << "Running estuarine dynamics simulation...\n";
    std::cout << "Estuary length: 50 km\n";
    std::cout << "River discharge: 500 m³/s\n";
    std::cout << "Tidal range: 2 m (M2)\n\n";
    
    double t = 0.0;
    double M2_period = 12.42 * 3600.0;
    
    while (t < config.end_time) {
        model.step();
        t = model.getCurrentTime();
        
        // Output every tidal quarter period
        if (std::fmod(t, M2_period / 4.0) < 60.0) {
            double salt_intrusion = model.computeSaltIntrusionLength(2.0);
            auto type = model.classifyEstuary();
            
            std::string type_str;
            switch (type) {
                case EstuarineModel::SALT_WEDGE: type_str = "Salt wedge"; break;
                case EstuarineModel::PARTIALLY_MIXED: type_str = "Partially mixed"; break;
                case EstuarineModel::WELL_MIXED: type_str = "Well mixed"; break;
            }
            
            std::cout << "Hour " << std::fixed << std::setprecision(1) << t / 3600.0
                      << ": Salt intrusion = " << salt_intrusion 
                      << " km, Type = " << type_str << "\n";
        }
    }
    
    // Final statistics
    double intrusion = model.computeSaltIntrusionLength();
    double flushing = model.computeFlushingTime();
    double Ri = model.computeEstuarineRichardson();
    
    std::cout << "\nEstuarine simulation completed.\n";
    std::cout << "Salt intrusion length (2 PSU): " << intrusion << " km\n";
    std::cout << "Flushing time: " << flushing / 86400.0 << " days\n";
    std::cout << "Estuarine Richardson number: " << Ri << "\n";
}

// =============================================================================
// Example 5: Coupled Wave-Current System
// =============================================================================

void runCoupledWaveCurrentExample() {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "Example 5: Coupled Wave-Current Interaction\n";
    std::cout << std::string(70, '=') << "\n\n";
    
    // This example demonstrates coupling between wave and circulation models
    // for nearshore applications where radiation stress is important
    
    std::cout << "Setting up coupled wave-current system...\n\n";
    
    // Simplified demonstration of coupling concepts
    int nx = 50, ny = 50;
    double dx = 100.0, dy = 100.0;  // 100m resolution
    
    // Create simple wave field propagating onshore
    std::vector<double> Hs(nx * ny);
    std::vector<double> dir(nx * ny);
    
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int idx = i + j * nx;
            double x_frac = static_cast<double>(i) / (nx - 1);
            
            // Wave height decreases due to breaking near shore
            double depth = 20.0 * (1.0 - x_frac);  // Sloping beach
            double Hs_deep = 2.0;
            
            if (depth < 0.78 * Hs_deep) {
                // Breaking zone
                Hs[idx] = 0.78 * depth;
            } else {
                // Shoaling (simplified)
                Hs[idx] = Hs_deep * std::pow(20.0 / std::max(depth, 1.0), 0.25);
            }
            
            dir[idx] = 270.0;  // Waves from west
        }
    }
    
    // Compute radiation stress gradient
    std::vector<double> dSxx_dx(nx * ny, 0.0);
    double rho = 1025.0;
    double g = 9.81;
    
    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {
            int idx = i + j * nx;
            int idx_e = (i + 1) + j * nx;
            int idx_w = (i - 1) + j * nx;
            
            // Radiation stress Sxx ≈ (1/8) ρ g H² (1 + cos²θ)
            // For shore-normal waves (θ = 0): Sxx = (3/16) ρ g H²
            double Sxx_e = (3.0/16.0) * rho * g * Hs[idx_e] * Hs[idx_e];
            double Sxx_w = (3.0/16.0) * rho * g * Hs[idx_w] * Hs[idx_w];
            
            dSxx_dx[idx] = (Sxx_e - Sxx_w) / (2.0 * dx);
        }
    }
    
    // Compute wave setup from radiation stress gradient balance
    // ∂η/∂x = -(1/(ρgh)) ∂Sxx/∂x
    std::vector<double> wave_setup(nx * ny, 0.0);
    
    for (int j = 0; j < ny; ++j) {
        // Integrate from offshore
        wave_setup[0 + j * nx] = 0.0;
        for (int i = 1; i < nx; ++i) {
            int idx = i + j * nx;
            double depth = 20.0 * (1.0 - static_cast<double>(i) / (nx - 1));
            depth = std::max(depth, 0.5);
            
            double deta_dx = -dSxx_dx[idx] / (rho * g * depth);
            wave_setup[idx] = wave_setup[(i-1) + j * nx] + deta_dx * dx;
        }
    }
    
    // Compute longshore current from radiation stress Sxy gradient
    // (simplified - full solution requires momentum balance)
    double max_setup = *std::max_element(wave_setup.begin(), wave_setup.end());
    
    std::cout << "Results for shore-normal waves (Hs = 2m offshore):\n";
    std::cout << "Maximum wave setup at shoreline: " << max_setup * 100 << " cm\n";
    std::cout << "Breaking zone width: ~" << int(2.0 / 0.78 / (20.0 / (nx * dx)) * dx) << " m\n";
    
    // Estimate longshore current using Longuet-Higgins formula
    double theta_break = 10.0 * M_PI / 180.0;  // 10 degree wave angle
    double C_f = 0.01;  // Friction coefficient
    double gamma_b = 0.78;
    double h_b = 2.0 / gamma_b;  // Breaking depth
    
    double V_ls = (5.0/16.0) * M_PI * std::sqrt(g * h_b) * gamma_b * std::sin(theta_break) * std::cos(theta_break) / C_f;
    
    std::cout << "Estimated longshore current (10° wave angle): " << V_ls << " m/s\n";
    
    std::cout << "\nCoupled wave-current demonstration completed.\n";
}

// =============================================================================
// Main Function
// =============================================================================

int main(int argc, char* argv[]) {
    std::cout << "================================================\n";
    std::cout << "FSRM Ocean Hydrodynamics and Physics Examples\n";
    std::cout << "================================================\n";
    std::cout << "\nThis program demonstrates various ocean modeling\n";
    std::cout << "capabilities including circulation, waves, storm surge,\n";
    std::cout << "and estuarine dynamics.\n";
    
    // Run examples
    try {
        // Example 1: Ocean circulation
        runOceanCirculationExample();
        
        // Example 2: Storm surge
        runStormSurgeExample();
        
        // Example 3: Wave model
        runWaveModelExample();
        
        // Example 4: Estuary
        runEstuarineExample();
        
        // Example 5: Wave-current coupling
        runCoupledWaveCurrentExample();
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    std::cout << "\n================================================\n";
    std::cout << "All ocean physics examples completed successfully!\n";
    std::cout << "================================================\n";
    
    return 0;
}
