/**
 * @file HydrodynamicModel.cpp
 * @brief Implementation of hydrodynamic modeling for ocean and coastal simulations
 */

#include "HydrodynamicModel.hpp"
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <numeric>

namespace FSRM {

// Physical constants
static const double G = 9.81;                    // Gravitational acceleration (m/s²)
static const double OMEGA = 7.2921e-5;           // Earth rotation rate (rad/s)
static const double EARTH_RADIUS = 6371000.0;    // Earth radius (m)
static const double DEG_TO_RAD = M_PI / 180.0;
static const double RAD_TO_DEG = 180.0 / M_PI;
static const double RHO_AIR = 1.225;             // Air density (kg/m³)
static const double RHO_WATER = 1025.0;          // Seawater density (kg/m³)

// =============================================================================
// WindForcing Implementation
// =============================================================================

void WindForcing::computeWindStress(double t, const std::vector<double>& x,
                                    const std::vector<double>& y,
                                    std::vector<double>& tau_x,
                                    std::vector<double>& tau_y) const {
    (void)y;  // Wind stress assumed spatially uniform (no y-dependence)
    size_t n = x.size();
    tau_x.resize(n);
    tau_y.resize(n);
    
    // Find time index for interpolation
    size_t idx = 0;
    double alpha = 0.0;
    
    if (!time.empty() && t > time[0]) {
        for (size_t i = 0; i < time.size() - 1; ++i) {
            if (t >= time[i] && t <= time[i + 1]) {
                idx = i;
                alpha = (t - time[i]) / (time[i + 1] - time[i]);
                break;
            }
        }
    }
    
    for (size_t i = 0; i < n; ++i) {
        double u10_val, v10_val;
        
        if (time.empty()) {
            // Constant wind
            u10_val = u10.empty() ? 0.0 : u10[std::min(i, u10.size() - 1)];
            v10_val = v10.empty() ? 0.0 : v10[std::min(i, v10.size() - 1)];
        } else {
            // Time-varying wind
            size_t si = std::min(i, u10_series[idx].size() - 1);
            u10_val = (1 - alpha) * u10_series[idx][si] + alpha * u10_series[idx + 1][si];
            v10_val = (1 - alpha) * v10_series[idx][si] + alpha * v10_series[idx + 1][si];
        }
        
        double wind_speed = std::sqrt(u10_val * u10_val + v10_val * v10_val);
        double Cd = use_bulk_formula ? dragCoefficientLargePond(wind_speed) : drag_coefficient;
        
        // Wind stress: τ = ρ_air * Cd * |U10| * U10
        tau_x[i] = air_density * Cd * wind_speed * u10_val;
        tau_y[i] = air_density * Cd * wind_speed * v10_val;
    }
}

double WindForcing::dragCoefficientLargePond(double wind_speed) {
    // Large and Pond (1981) drag coefficient
    if (wind_speed < 11.0) {
        return 1.2e-3;
    } else if (wind_speed < 25.0) {
        return (0.49 + 0.065 * wind_speed) * 1e-3;
    } else {
        // Cap at high wind speeds
        return 2.5e-3;
    }
}

double WindForcing::dragCoefficientCOARE(double wind_speed, double sst,
                                         double air_temp, double humidity) {
    (void)humidity;  // Full COARE uses humidity; simplified version does not
    // Simplified COARE 3.0 algorithm
    double Cd = 1.0e-3;
    
    // Neutral drag coefficient
    if (wind_speed < 10.0) {
        Cd = 1.0e-3 * (0.61 + 0.063 * wind_speed);
    } else {
        Cd = 1.0e-3 * (0.61 + 0.063 * 10.0 + 0.066 * (wind_speed - 10.0));
    }
    
    // Stability correction (simplified)
    double dT = sst - air_temp;
    if (dT > 0) {
        // Unstable - increase Cd
        Cd *= (1.0 + 0.02 * dT);
    } else {
        // Stable - decrease Cd
        Cd *= (1.0 + 0.01 * dT);
    }
    
    return std::max(0.5e-3, std::min(Cd, 3.0e-3));
}

// =============================================================================
// AtmosphericPressureForcing Implementation
// =============================================================================

double AtmosphericPressureForcing::computeInverseBarometer(double p,
                                                           double rho_water) const {
    if (!enable_inverse_barometer) return 0.0;
    
    // Sea level change ≈ -(p - p_ref) / (ρ * g)
    // Approximately 1 cm per hPa
    return -(p - reference_pressure) / (rho_water * G);
}

// =============================================================================
// TidalForcing Implementation
// =============================================================================

double TidalForcing::getTidalElevation(double t, double lon, double lat) const {
    (void)lon; (void)lat;  // Spatial variation requires nodal correction data
    double eta = 0.0;
    
    for (const auto& c : constituents) {
        double phase_rad = c.phase * DEG_TO_RAD;
        double arg = c.frequency * t + c.equilibrium_argument + c.nodal_factor - phase_rad;
        eta += c.amplitude * std::cos(arg);
    }
    
    return eta;
}

void TidalForcing::getTidalCurrent(double t, double lon, double lat,
                                   double& u, double& v) const {
    (void)lon;  // Spatial variation requires full tidal atlas
    // Simplified tidal current estimation
    // In reality, this requires additional constituent data (major/minor axes, inclination)
    u = v = 0.0;
    
    // Coriolis parameter for depth-averaged current estimate
    (void)(2.0 * OMEGA * std::sin(lat * DEG_TO_RAD));  // f, used in full implementation
    double h = 100.0;  // Assume 100m depth for basic estimate
    
    for (const auto& c : constituents) {
        double phase_rad = c.phase * DEG_TO_RAD;
        double arg = c.frequency * t + c.equilibrium_argument + c.nodal_factor - phase_rad;
        
        // Simple estimate: u ≈ (g/h) * η / ω * sin(ωt + φ)
        double u_amp = (G / h) * c.amplitude / c.frequency;
        u += u_amp * std::sin(arg);
    }
}

bool TidalForcing::loadFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) return false;
    
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        std::istringstream iss(line);
        std::string name;
        double amp, phase;
        
        if (iss >> name >> amp >> phase) {
            Constituent c;
            c.amplitude = amp;
            c.phase = phase;
            c.nodal_factor = 1.0;
            c.equilibrium_argument = 0.0;
            
            if (name == "M2") {
                c.type = TidalConstituent::M2;
                c.frequency = 2.0 * M_PI / (12.4206 * 3600.0);
            } else if (name == "S2") {
                c.type = TidalConstituent::S2;
                c.frequency = 2.0 * M_PI / (12.0 * 3600.0);
            } else if (name == "N2") {
                c.type = TidalConstituent::N2;
                c.frequency = 2.0 * M_PI / (12.6583 * 3600.0);
            } else if (name == "K1") {
                c.type = TidalConstituent::K1;
                c.frequency = 2.0 * M_PI / (23.9345 * 3600.0);
            } else if (name == "O1") {
                c.type = TidalConstituent::O1;
                c.frequency = 2.0 * M_PI / (25.8193 * 3600.0);
            }
            
            constituents.push_back(c);
        }
    }
    
    return !constituents.empty();
}

void TidalForcing::setupM2Only() {
    constituents.clear();
    Constituent c;
    c.type = TidalConstituent::M2;
    c.amplitude = 1.0;  // To be set by user
    c.phase = 0.0;
    c.frequency = 2.0 * M_PI / (12.4206 * 3600.0);
    c.nodal_factor = 1.0;
    c.equilibrium_argument = 0.0;
    constituents.push_back(c);
}

void TidalForcing::setupPrimaryConstituents() {
    constituents.clear();
    
    // M2 - Principal lunar semidiurnal
    constituents.push_back({TidalConstituent::M2, 1.0, 0.0,
                           2.0 * M_PI / (12.4206 * 3600.0), 1.0, 0.0});
    
    // S2 - Principal solar semidiurnal
    constituents.push_back({TidalConstituent::S2, 0.46, 0.0,
                           2.0 * M_PI / (12.0 * 3600.0), 1.0, 0.0});
    
    // N2 - Larger lunar elliptic
    constituents.push_back({TidalConstituent::N2, 0.19, 0.0,
                           2.0 * M_PI / (12.6583 * 3600.0), 1.0, 0.0});
    
    // K1 - Lunar diurnal
    constituents.push_back({TidalConstituent::K1, 0.58, 0.0,
                           2.0 * M_PI / (23.9345 * 3600.0), 1.0, 0.0});
    
    // O1 - Principal lunar diurnal
    constituents.push_back({TidalConstituent::O1, 0.41, 0.0,
                           2.0 * M_PI / (25.8193 * 3600.0), 1.0, 0.0});
}

void TidalForcing::setupFullConstituents() {
    setupPrimaryConstituents();
    
    // K2 - Lunisolar semidiurnal
    constituents.push_back({TidalConstituent::K2, 0.13, 0.0,
                           2.0 * M_PI / (11.9672 * 3600.0), 1.0, 0.0});
    
    // P1 - Principal solar diurnal
    constituents.push_back({TidalConstituent::P1, 0.19, 0.0,
                           2.0 * M_PI / (24.0659 * 3600.0), 1.0, 0.0});
    
    // Q1 - Larger lunar elliptic diurnal
    constituents.push_back({TidalConstituent::Q1, 0.08, 0.0,
                           2.0 * M_PI / (26.8684 * 3600.0), 1.0, 0.0});
    
    // M4 - Shallow water overtide
    constituents.push_back({TidalConstituent::M4, 0.02, 0.0,
                           4.0 * M_PI / (12.4206 * 3600.0), 1.0, 0.0});
}

// =============================================================================
// RiverDischarge Implementation
// =============================================================================

double RiverDischarge::getDischarge(double t) const {
    if (time.empty() || discharge_series.empty()) {
        return discharge;
    }
    
    // Linear interpolation
    for (size_t i = 0; i < time.size() - 1; ++i) {
        if (t >= time[i] && t <= time[i + 1]) {
            double alpha = (t - time[i]) / (time[i + 1] - time[i]);
            return (1 - alpha) * discharge_series[i] + alpha * discharge_series[i + 1];
        }
    }
    
    return discharge_series.back();
}

// =============================================================================
// HydrodynamicModel Implementation
// =============================================================================

HydrodynamicModel::HydrodynamicModel()
    : current_time(0.0), dt(0.0), dt_baro(0.0), step_count(0) {}

HydrodynamicModel::~HydrodynamicModel() = default;

void HydrodynamicModel::initialize(const HydrodynamicConfig& cfg) {
    config = cfg;
    
    nx = cfg.nx;
    ny = cfg.ny;
    nz = (cfg.model_type == HydrodynamicModelType::SHALLOW_WATER_2D) ? 1 : cfg.nz;
    
    dx = (cfg.lon_max - cfg.lon_min) * EARTH_RADIUS * DEG_TO_RAD * 
         std::cos(0.5 * (cfg.lat_min + cfg.lat_max) * DEG_TO_RAD) / (nx - 1);
    dy = (cfg.lat_max - cfg.lat_min) * EARTH_RADIUS * DEG_TO_RAD / (ny - 1);
    
    lon_min = cfg.lon_min;
    lat_min = cfg.lat_min;
    
    dt = cfg.dt;
    dt_baro = cfg.dt_barotropic > 0 ? cfg.dt_barotropic : dt / cfg.barotropic_steps;
    
    // Allocate 2D arrays
    int n2d = nx * ny;
    eta.resize(n2d, 0.0);
    u_bar.resize(n2d, 0.0);
    v_bar.resize(n2d, 0.0);
    h.resize(n2d, 100.0);  // Default 100m depth
    H.resize(n2d, 100.0);
    mask.resize(n2d, 1);
    
    eta_old.resize(n2d, 0.0);
    u_bar_old.resize(n2d, 0.0);
    v_bar_old.resize(n2d, 0.0);
    
    deta_dt.resize(n2d, 0.0);
    du_dt.resize(n2d, 0.0);
    dv_dt.resize(n2d, 0.0);
    
    tau_x.resize(n2d, 0.0);
    tau_y.resize(n2d, 0.0);
    p_atm.resize(n2d, 101325.0);
    
    // Allocate 3D arrays if needed
    if (nz > 1) {
        int n3d = nx * ny * nz;
        u.resize(n3d, 0.0);
        v.resize(n3d, 0.0);
        w.resize(n3d, 0.0);
        temperature.resize(n3d, 15.0);  // Default 15°C
        salinity.resize(n3d, 35.0);     // Default 35 PSU
        density.resize(n3d, RHO_WATER);
        
        u_old.resize(n3d, 0.0);
        v_old.resize(n3d, 0.0);
        
        dT_dt.resize(n3d, 0.0);
        dS_dt.resize(n3d, 0.0);
        
        dP_dx.resize(n3d, 0.0);
        dP_dy.resize(n3d, 0.0);
        
        Kv.resize(n3d, cfg.vertical_diffusivity);
        Av.resize(n3d, cfg.vertical_viscosity);
        
        computeSigmaCoordinates();
    }
    
    Kh.resize(n2d, cfg.horizontal_diffusivity);
    Ah.resize(n2d, cfg.horizontal_viscosity);
    
    current_time = 0.0;
    step_count = 0;
    
    std::cout << "HydrodynamicModel initialized:\n"
              << "  Grid: " << nx << " x " << ny << " x " << nz << "\n"
              << "  Resolution: " << dx/1000.0 << " km x " << dy/1000.0 << " km\n"
              << "  Time step: " << dt << " s (baroclinic), " << dt_baro << " s (barotropic)\n";
}

void HydrodynamicModel::setBathymetry(const BathymetryGrid& bathymetry) {
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            double lon = lon_min + i * (config.lon_max - config.lon_min) / (nx - 1);
            double lat = lat_min + j * (config.lat_max - config.lat_min) / (ny - 1);
            
            double depth = bathymetry.getDepth(lon, lat);
            int idx = idx2d(i, j);
            
            h[idx] = std::max(config.min_depth, depth);
            H[idx] = h[idx] + eta[idx];
            mask[idx] = (depth > 0) ? 1 : 0;
        }
    }
}

void HydrodynamicModel::setBathymetry(const std::vector<double>& depth) {
    if (depth.size() != static_cast<size_t>(nx * ny)) {
        std::cerr << "Error: Bathymetry size mismatch\n";
        return;
    }
    
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int idx = idx2d(i, j);
            h[idx] = std::max(config.min_depth, depth[idx]);
            H[idx] = h[idx] + eta[idx];
            mask[idx] = (depth[idx] > 0) ? 1 : 0;
        }
    }
}

void HydrodynamicModel::setInitialElevation(const std::vector<double>& eta_init) {
    eta = eta_init;
    eta_old = eta_init;
    
    // Update total depth
    for (int i = 0; i < nx * ny; ++i) {
        H[i] = h[i] + eta[i];
    }
}

void HydrodynamicModel::setInitialVelocity(const std::vector<double>& u_init,
                                            const std::vector<double>& v_init) {
    u_bar = u_init;
    v_bar = v_init;
    u_bar_old = u_init;
    v_bar_old = v_init;
    
    if (nz > 1) {
        // Initialize 3D velocities from depth-averaged
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                int idx2 = idx2d(i, j);
                for (int k = 0; k < nz; ++k) {
                    int idx3 = idx3d(i, j, k);
                    u[idx3] = u_bar[idx2];
                    v[idx3] = v_bar[idx2];
                }
            }
        }
        u_old = u;
        v_old = v;
    }
}

void HydrodynamicModel::setInitialTemperature(const std::vector<double>& T) {
    temperature = T;
    
    // Update density
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            for (int k = 0; k < nz; ++k) {
                int idx = idx3d(i, j, k);
                double p = -z_at_sigma(i, j, k) / 10.0;  // Approximate pressure in dbar
                density[idx] = computeDensity(temperature[idx], salinity[idx], p);
            }
        }
    }
}

void HydrodynamicModel::setInitialSalinity(const std::vector<double>& S) {
    salinity = S;
    
    // Update density
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            for (int k = 0; k < nz; ++k) {
                int idx = idx3d(i, j, k);
                double p = -z_at_sigma(i, j, k) / 10.0;
                density[idx] = computeDensity(temperature[idx], salinity[idx], p);
            }
        }
    }
}

void HydrodynamicModel::addOpenBoundary(const OpenBoundaryCondition& bc) {
    open_boundaries.push_back(bc);
}

void HydrodynamicModel::setTidalForcing(const TidalForcing& tidal) {
    tidal_forcing = tidal;
}

void HydrodynamicModel::addRiverInput(const RiverDischarge& river) {
    rivers.push_back(river);
}

void HydrodynamicModel::setWindForcing(const WindForcing& wind) {
    wind_forcing = wind;
}

void HydrodynamicModel::setAtmosphericPressure(const AtmosphericPressureForcing& pressure) {
    pressure_forcing = pressure;
}

void HydrodynamicModel::computeSigmaCoordinates() {
    sigma.resize(nz);
    
    // Stretched sigma coordinates (Song and Haidvogel, 1994)
    for (int k = 0; k < nz; ++k) {
        double s = -1.0 + (k + 0.5) / nz;  // -1 < s < 0
        
        // Surface stretching
        double C_s = (1.0 - std::cosh(config.theta_s * s)) / 
                    (std::cosh(config.theta_s) - 1.0);
        
        // Bottom stretching
        double C_b = (std::exp(config.theta_b * C_s) - 1.0) / 
                    (1.0 - std::exp(-config.theta_b));
        
        sigma[k] = C_b;
    }
}

double HydrodynamicModel::z_at_sigma(int i, int j, int k) const {
    if (nz == 1) return -0.5 * H[idx2d(i, j)];
    
    int idx = idx2d(i, j);
    // Note: config.h_c (critical depth) used in stretched sigma coordinates
    
    // z = eta + (eta + h) * sigma  [for simple sigma]
    // For stretched: z = h*s + (h_c*s + h*C)/(h_c + h) * (eta - eta0)
    double depth = H[idx];
    return sigma[k] * depth + eta[idx];
}

double HydrodynamicModel::dz_at_sigma(int i, int j, int k) const {
    if (nz == 1) return H[idx2d(i, j)];
    
    int idx = idx2d(i, j);
    double depth = H[idx];
    
    // Layer thickness
    if (k == 0) {
        return (sigma[0] - (-1.0)) * depth;
    } else if (k == nz - 1) {
        return (0.0 - sigma[nz - 1]) * depth;
    } else {
        return (sigma[k] - sigma[k - 1]) * depth * 0.5 +
               (sigma[k + 1] - sigma[k]) * depth * 0.5;
    }
}

double HydrodynamicModel::f_coriolis(int j) const {
    if (!config.use_coriolis) return 0.0;
    
    double lat = lat_min + j * (config.lat_max - config.lat_min) / (ny - 1);
    double f = 2.0 * OMEGA * std::sin(lat * DEG_TO_RAD);
    
    if (config.use_beta_plane) {
        // Use beta-plane approximation: f = f0 + beta * y
        double y = (lat - config.latitude_reference) * DEG_TO_RAD * EARTH_RADIUS;
        f = config.f0 + config.beta * y;
    }
    
    return f;
}

double HydrodynamicModel::cell_area(int i, int j) const {
    (void)i;  // Cell area depends only on latitude (j), not longitude (i)
    double lat = lat_min + j * (config.lat_max - config.lat_min) / (ny - 1);
    return dx * dy * std::cos(lat * DEG_TO_RAD);
}

double HydrodynamicModel::step() {
    double dt_actual = dt;
    
    // Compute CFL condition for adaptive time stepping
    if (config.adaptive_timestep) {
        double max_speed = 0.0;
        for (int j = 1; j < ny - 1; ++j) {
            for (int i = 1; i < nx - 1; ++i) {
                int idx = idx2d(i, j);
                if (mask[idx] == 0) continue;
                
                double c = std::sqrt(G * H[idx]);  // Gravity wave speed
                double u_speed = std::sqrt(u_bar[idx] * u_bar[idx] + 
                                          v_bar[idx] * v_bar[idx]);
                max_speed = std::max(max_speed, c + u_speed);
            }
        }
        
        if (max_speed > 0) {
            double dt_cfl = config.cfl_number * std::min(dx, dy) / max_speed;
            dt_actual = std::min(dt, dt_cfl);
        }
    }
    
    // Store old values for leapfrog or Adams-Bashforth
    std::swap(eta_old, eta);
    std::swap(u_bar_old, u_bar);
    std::swap(v_bar_old, v_bar);
    if (nz > 1) {
        std::swap(u_old, u);
        std::swap(v_old, v);
    }
    eta = eta_old;
    u_bar = u_bar_old;
    v_bar = v_bar_old;
    if (nz > 1) {
        u = u_old;
        v = v_old;
    }
    
    // Split-explicit time stepping
    if (config.time_scheme == TimeSteppingScheme::SPLIT_EXPLICIT && nz > 1) {
        // Multiple barotropic steps
        int n_baro = static_cast<int>(dt_actual / dt_baro);
        for (int n = 0; n < n_baro; ++n) {
            stepBarotropic();
        }
        // Single baroclinic step
        stepBaroclinic();
    } else {
        // 2D mode or non-split
        stepBarotropic();
    }
    
    // Apply wetting/drying
    if (config.enable_wetting_drying) {
        applyWettingDrying();
    }
    
    current_time += dt_actual;
    step_count++;
    
    return dt_actual;
}

void HydrodynamicModel::stepBarotropic() {
    // Compute surface forcing
    computeSurfaceForcing();
    
    // Compute advection of depth-averaged momentum
    computeAdvection();
    
    // Compute Coriolis
    computeCoriolis();
    
    // Compute horizontal mixing
    computeHorizontalMixing();
    
    // Compute bottom friction
    computeBottomFriction();
    
    // Compute barotropic pressure gradient
    computeBarotropicPressureGradient();
    
    // Update velocities
    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {
            int idx = idx2d(i, j);
            if (mask[idx] == 0) continue;
            
            u_bar[idx] += dt_baro * du_dt[idx];
            v_bar[idx] += dt_baro * dv_dt[idx];
        }
    }
    
    // Update elevation (continuity equation)
    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {
            int idx = idx2d(i, j);
            if (mask[idx] == 0) continue;
            
            // ∂η/∂t = -∇·(H*u_bar)
            int idxE = idx2d(i + 1, j);
            int idxW = idx2d(i - 1, j);
            int idxN = idx2d(i, j + 1);
            int idxS = idx2d(i, j - 1);
            
            double flux_x = (H[idxE] * u_bar[idxE] - H[idxW] * u_bar[idxW]) / (2.0 * dx);
            double flux_y = (H[idxN] * v_bar[idxN] - H[idxS] * v_bar[idxS]) / (2.0 * dy);
            
            eta[idx] -= dt_baro * (flux_x + flux_y);
        }
    }
    
    // Update total depth
    for (int idx = 0; idx < nx * ny; ++idx) {
        H[idx] = h[idx] + eta[idx];
        if (H[idx] < config.dry_depth) {
            H[idx] = config.dry_depth;
        }
    }
    
    // Apply boundary conditions
    applyBarotropicBoundaryConditions();
}

void HydrodynamicModel::stepBaroclinic() {
    if (nz <= 1) return;
    
    // Compute baroclinic pressure gradient
    computeBaroclinicPressureGradient();
    
    // Solve tracer equations (temperature, salinity)
    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {
            for (int k = 0; k < nz; ++k) {
                int idx = idx3d(i, j, k);
                int idx2 = idx2d(i, j);
                if (mask[idx2] == 0) continue;
                
                // Simple upwind advection for tracers
                int idxE = idx3d(i + 1, j, k);
                int idxW = idx3d(i - 1, j, k);
                int idxN = idx3d(i, j + 1, k);
                int idxS = idx3d(i, j - 1, k);
                
                double dT_dx = (u[idx] > 0) ? 
                    (temperature[idx] - temperature[idxW]) / dx :
                    (temperature[idxE] - temperature[idx]) / dx;
                double dT_dy = (v[idx] > 0) ?
                    (temperature[idx] - temperature[idxS]) / dy :
                    (temperature[idxN] - temperature[idx]) / dy;
                
                dT_dt[idx] = -u[idx] * dT_dx - v[idx] * dT_dy;
                
                double dS_dx = (u[idx] > 0) ?
                    (salinity[idx] - salinity[idxW]) / dx :
                    (salinity[idxE] - salinity[idx]) / dx;
                double dS_dy = (v[idx] > 0) ?
                    (salinity[idx] - salinity[idxS]) / dy :
                    (salinity[idxN] - salinity[idx]) / dy;
                
                dS_dt[idx] = -u[idx] * dS_dx - v[idx] * dS_dy;
            }
        }
    }
    
    // Solve vertical diffusion implicitly
    solveVerticalDiffusion();
    
    // Update tracers
    for (int idx = 0; idx < nx * ny * nz; ++idx) {
        temperature[idx] += dt * dT_dt[idx];
        salinity[idx] += dt * dS_dt[idx];
        
        // Bounds checking
        temperature[idx] = std::max(-2.0, std::min(temperature[idx], 40.0));
        salinity[idx] = std::max(0.0, std::min(salinity[idx], 42.0));
    }
    
    // Update density
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            for (int k = 0; k < nz; ++k) {
                int idx = idx3d(i, j, k);
                double p = -z_at_sigma(i, j, k) / 10.0;
                density[idx] = computeDensity(temperature[idx], salinity[idx], p);
            }
        }
    }
    
    // Update 3D velocities
    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {
            int idx2 = idx2d(i, j);
            if (mask[idx2] == 0) continue;
            
            for (int k = 0; k < nz; ++k) {
                int idx = idx3d(i, j, k);
                
                // Add baroclinic pressure gradient
                u[idx] += dt * (-dP_dx[idx] / config.rho_0);
                v[idx] += dt * (-dP_dy[idx] / config.rho_0);
            }
        }
    }
    
    // Compute vertical velocity from continuity
    computeVerticalVelocity();
}

void HydrodynamicModel::computeBarotropicPressureGradient() {
    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {
            int idx = idx2d(i, j);
            if (mask[idx] == 0) {
                du_dt[idx] = dv_dt[idx] = 0.0;
                continue;
            }
            
            int idxE = idx2d(i + 1, j);
            int idxW = idx2d(i - 1, j);
            int idxN = idx2d(i, j + 1);
            int idxS = idx2d(i, j - 1);
            
            // ∂η/∂x and ∂η/∂y
            double deta_dx = (eta[idxE] - eta[idxW]) / (2.0 * dx);
            double deta_dy = (eta[idxN] - eta[idxS]) / (2.0 * dy);
            
            // Pressure gradient force: -g * ∂η/∂x, -g * ∂η/∂y
            du_dt[idx] += -G * deta_dx;
            dv_dt[idx] += -G * deta_dy;
            
            // Add atmospheric pressure forcing (inverse barometer)
            double dp_dx = (p_atm[idxE] - p_atm[idxW]) / (2.0 * dx);
            double dp_dy = (p_atm[idxN] - p_atm[idxS]) / (2.0 * dy);
            du_dt[idx] += -dp_dx / (config.rho_0);
            dv_dt[idx] += -dp_dy / (config.rho_0);
        }
    }
}

void HydrodynamicModel::computeBaroclinicPressureGradient() {
    if (nz <= 1) return;
    
    // Integrate density anomaly from surface
    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {
            int idx2 = idx2d(i, j);
            if (mask[idx2] == 0) continue;
            
            // Horizontal density gradients (3D indices computed in loop)
            double P_int = 0.0;  // Integrated pressure
            
            for (int k = nz - 1; k >= 0; --k) {
                int idx = idx3d(i, j, k);
                int idxE3 = idx3d(i + 1, j, k);
                int idxW3 = idx3d(i - 1, j, k);
                int idxN3 = idx3d(i, j + 1, k);
                int idxS3 = idx3d(i, j - 1, k);
                
                double dz = dz_at_sigma(i, j, k);
                
                // Integrate: P = ∫ ρ g dz
                P_int += density[idx] * G * dz;
                
                // Horizontal gradient at this level
                double drho_dx = (density[idxE3] - density[idxW3]) / (2.0 * dx);
                double drho_dy = (density[idxN3] - density[idxS3]) / (2.0 * dy);
                
                // Baroclinic pressure gradient
                dP_dx[idx] = G * drho_dx * dz;
                dP_dy[idx] = G * drho_dy * dz;
            }
        }
    }
}

void HydrodynamicModel::applyBarotropicBoundaryConditions() {
    // Apply boundary conditions for each open boundary
    for (const auto& bc : open_boundaries) {
        switch (bc.location) {
            case OpenBoundaryCondition::WEST:
                for (int j = 0; j < ny; ++j) {
                    int idx = idx2d(0, j);
                    if (bc.type == OpenBoundaryCondition::CLAMPED) {
                        eta[idx] = tidal_forcing.getTidalElevation(current_time);
                        u_bar[idx] = 0.0;
                    } else if (bc.type == OpenBoundaryCondition::RADIATION) {
                        double c = std::sqrt(G * H[idx]);
                        eta[idx] += -c * dt_baro / dx * (eta[idx] - eta[idx2d(1, j)]);
                    }
                }
                break;
                
            case OpenBoundaryCondition::EAST:
                for (int j = 0; j < ny; ++j) {
                    int idx = idx2d(nx - 1, j);
                    if (bc.type == OpenBoundaryCondition::CLAMPED) {
                        eta[idx] = tidal_forcing.getTidalElevation(current_time);
                        u_bar[idx] = 0.0;
                    } else if (bc.type == OpenBoundaryCondition::RADIATION) {
                        double c = std::sqrt(G * H[idx]);
                        eta[idx] += c * dt_baro / dx * (eta[idx2d(nx - 2, j)] - eta[idx]);
                    }
                }
                break;
                
            case OpenBoundaryCondition::SOUTH:
                for (int i = 0; i < nx; ++i) {
                    int idx = idx2d(i, 0);
                    if (bc.type == OpenBoundaryCondition::CLAMPED) {
                        eta[idx] = tidal_forcing.getTidalElevation(current_time);
                        v_bar[idx] = 0.0;
                    } else if (bc.type == OpenBoundaryCondition::RADIATION) {
                        double c = std::sqrt(G * H[idx]);
                        eta[idx] += -c * dt_baro / dy * (eta[idx] - eta[idx2d(i, 1)]);
                    }
                }
                break;
                
            case OpenBoundaryCondition::NORTH:
                for (int i = 0; i < nx; ++i) {
                    int idx = idx2d(i, ny - 1);
                    if (bc.type == OpenBoundaryCondition::CLAMPED) {
                        eta[idx] = tidal_forcing.getTidalElevation(current_time);
                        v_bar[idx] = 0.0;
                    } else if (bc.type == OpenBoundaryCondition::RADIATION) {
                        double c = std::sqrt(G * H[idx]);
                        eta[idx] += c * dt_baro / dy * (eta[idx2d(i, ny - 2)] - eta[idx]);
                    }
                }
                break;
        }
    }
    
    // Land boundaries (no-slip or free-slip)
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int idx = idx2d(i, j);
            if (mask[idx] == 0) {
                eta[idx] = 0.0;
                u_bar[idx] = 0.0;
                v_bar[idx] = 0.0;
            }
        }
    }
}

void HydrodynamicModel::solveVerticalDiffusion() {
    if (nz <= 1) return;
    
    // Tridiagonal solver for implicit vertical diffusion
    std::vector<double> a(nz), b(nz), c(nz), d(nz), x(nz);
    
    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {
            int idx2 = idx2d(i, j);
            if (mask[idx2] == 0) continue;
            
            // Solve for temperature
            for (int k = 0; k < nz; ++k) {
                double dz = dz_at_sigma(i, j, k);
                double Kv_k = Kv[idx3d(i, j, k)];
                
                if (k == 0) {
                    // Surface boundary: flux = 0 (or specified heat flux)
                    a[k] = 0.0;
                    b[k] = 1.0 + dt * Kv_k / (dz * dz);
                    c[k] = -dt * Kv_k / (dz * dz);
                } else if (k == nz - 1) {
                    // Bottom boundary: flux = 0
                    a[k] = -dt * Kv_k / (dz * dz);
                    b[k] = 1.0 + dt * Kv_k / (dz * dz);
                    c[k] = 0.0;
                } else {
                    a[k] = -dt * Kv_k / (dz * dz);
                    b[k] = 1.0 + 2.0 * dt * Kv_k / (dz * dz);
                    c[k] = -dt * Kv_k / (dz * dz);
                }
                
                d[k] = temperature[idx3d(i, j, k)];
            }
            
            // Thomas algorithm
            for (int k = 1; k < nz; ++k) {
                double w = a[k] / b[k - 1];
                b[k] -= w * c[k - 1];
                d[k] -= w * d[k - 1];
            }
            
            x[nz - 1] = d[nz - 1] / b[nz - 1];
            for (int k = nz - 2; k >= 0; --k) {
                x[k] = (d[k] - c[k] * x[k + 1]) / b[k];
            }
            
            for (int k = 0; k < nz; ++k) {
                temperature[idx3d(i, j, k)] = x[k];
            }
            
            // Similar for salinity...
        }
    }
}

void HydrodynamicModel::computeAdvection() {
    // MUSCL reconstruction with minmod limiter
    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {
            int idx = idx2d(i, j);
            if (mask[idx] == 0) continue;
            
            // x-direction momentum advection
            int idxE = idx2d(i + 1, j);
            int idxW = idx2d(i - 1, j);
            
            double flux_x = 0.0;
            if (u_bar[idx] > 0) {
                double qL = u_bar[idxW];
                double qC = u_bar[idx];
                double qR = u_bar[idxE];
                double q_face = muscl_reconstruct(qL, qC, qR, 'R');
                flux_x = u_bar[idx] * q_face;
            } else {
                double qL = u_bar[idx];
                double qC = u_bar[idxE];
                int idxEE = idx2d(std::min(i + 2, nx - 1), j);
                double qR = u_bar[idxEE];
                double q_face = muscl_reconstruct(qL, qC, qR, 'L');
                flux_x = u_bar[idx] * q_face;
            }
            
            // y-direction momentum advection
            int idxN = idx2d(i, j + 1);
            int idxS = idx2d(i, j - 1);
            
            double flux_y = 0.0;
            if (v_bar[idx] > 0) {
                double qL = u_bar[idxS];
                double qC = u_bar[idx];
                double qR = u_bar[idxN];
                double q_face = muscl_reconstruct(qL, qC, qR, 'R');
                flux_y = v_bar[idx] * q_face;
            } else {
                double qL = u_bar[idx];
                double qC = u_bar[idxN];
                int idxNN = idx2d(i, std::min(j + 2, ny - 1));
                double qR = u_bar[idxNN];
                double q_face = muscl_reconstruct(qL, qC, qR, 'L');
                flux_y = v_bar[idx] * q_face;
            }
            
            // du/dt from advection
            double adv_u = (flux_x - u_bar[idxW] * u_bar[idxW]) / dx +
                          (flux_y - v_bar[idxS] * u_bar[idxS]) / dy;
            
            // Similar for v...
            du_dt[idx] = -adv_u;
        }
    }
}

double HydrodynamicModel::minmod(double a, double b) const {
    if (a * b <= 0) return 0.0;
    return (std::abs(a) < std::abs(b)) ? a : b;
}

double HydrodynamicModel::muscl_reconstruct(double qL, double qC, double qR, char side) const {
    double dqL = qC - qL;
    double dqR = qR - qC;
    double slope = minmod(dqL, dqR);
    
    if (side == 'L') {
        return qC - 0.5 * slope;
    } else {
        return qC + 0.5 * slope;
    }
}

void HydrodynamicModel::computeCoriolis() {
    for (int j = 1; j < ny - 1; ++j) {
        double f = f_coriolis(j);
        
        for (int i = 1; i < nx - 1; ++i) {
            int idx = idx2d(i, j);
            if (mask[idx] == 0) continue;
            
            // Coriolis: du/dt = +f*v, dv/dt = -f*u
            du_dt[idx] += f * v_bar[idx];
            dv_dt[idx] += -f * u_bar[idx];
        }
    }
}

void HydrodynamicModel::computeHorizontalMixing() {
    for (int j = 2; j < ny - 2; ++j) {
        for (int i = 2; i < nx - 2; ++i) {
            int idx = idx2d(i, j);
            if (mask[idx] == 0) continue;
            
            double Ah_local = Ah[idx];
            
            // Laplacian of u
            double lap_u = (u_bar[idx2d(i + 1, j)] - 2.0 * u_bar[idx] + u_bar[idx2d(i - 1, j)]) / (dx * dx) +
                          (u_bar[idx2d(i, j + 1)] - 2.0 * u_bar[idx] + u_bar[idx2d(i, j - 1)]) / (dy * dy);
            
            double lap_v = (v_bar[idx2d(i + 1, j)] - 2.0 * v_bar[idx] + v_bar[idx2d(i - 1, j)]) / (dx * dx) +
                          (v_bar[idx2d(i, j + 1)] - 2.0 * v_bar[idx] + v_bar[idx2d(i, j - 1)]) / (dy * dy);
            
            du_dt[idx] += Ah_local * lap_u;
            dv_dt[idx] += Ah_local * lap_v;
        }
    }
}

void HydrodynamicModel::computeBottomFriction() {
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int idx = idx2d(i, j);
            if (mask[idx] == 0) continue;
            
            double speed = std::sqrt(u_bar[idx] * u_bar[idx] + v_bar[idx] * v_bar[idx]);
            if (speed < 1e-10) continue;
            
            double tau_b_x = 0.0, tau_b_y = 0.0;
            
            switch (config.friction_type) {
                case BottomFrictionType::LINEAR: {
                    double r = 0.001;  // Linear friction coefficient
                    tau_b_x = r * u_bar[idx];
                    tau_b_y = r * v_bar[idx];
                    break;
                }
                
                case BottomFrictionType::QUADRATIC: {
                    double Cd = config.bottom_drag_coefficient;
                    tau_b_x = Cd * speed * u_bar[idx] / H[idx];
                    tau_b_y = Cd * speed * v_bar[idx] / H[idx];
                    break;
                }
                
                case BottomFrictionType::MANNING: {
                    double n = config.manning_n;
                    double Cf = G * n * n / std::pow(H[idx], 1.0/3.0);
                    tau_b_x = Cf * speed * u_bar[idx] / H[idx];
                    tau_b_y = Cf * speed * v_bar[idx] / H[idx];
                    break;
                }
                
                case BottomFrictionType::LOG_LAYER: {
                    double z0 = config.roughness_length;
                    double kappa = 0.41;  // von Karman constant
                    double z_ref = 0.1 * H[idx];  // Reference height
                    double Cd = std::pow(kappa / std::log(z_ref / z0), 2);
                    tau_b_x = Cd * speed * u_bar[idx] / H[idx];
                    tau_b_y = Cd * speed * v_bar[idx] / H[idx];
                    break;
                }
                
                default:
                    break;
            }
            
            // Apply friction (explicit form; implicit factor available for stability)
            // factor = 1.0 / (1.0 + dt * sqrt(tau_x² + tau_y²) / speed)
            du_dt[idx] -= tau_b_x;
            dv_dt[idx] -= tau_b_y;
        }
    }
}

void HydrodynamicModel::computeSurfaceForcing() {
    // Compute wind stress
    std::vector<double> x(nx * ny), y(nx * ny);
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int idx = idx2d(i, j);
            x[idx] = lon_min + i * (config.lon_max - config.lon_min) / (nx - 1);
            y[idx] = lat_min + j * (config.lat_max - config.lat_min) / (ny - 1);
        }
    }
    
    wind_forcing.computeWindStress(current_time, x, y, tau_x, tau_y);
    
    // Apply wind stress as surface forcing
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int idx = idx2d(i, j);
            if (mask[idx] == 0) continue;
            
            // Wind stress contribution to depth-averaged velocity
            // τ_wind / (ρ * H)
            du_dt[idx] += tau_x[idx] / (config.rho_0 * H[idx]);
            dv_dt[idx] += tau_y[idx] / (config.rho_0 * H[idx]);
        }
    }
}

void HydrodynamicModel::applyWettingDrying() {
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int idx = idx2d(i, j);
            
            H[idx] = h[idx] + eta[idx];
            
            if (H[idx] < config.dry_depth) {
                // Cell is dry
                mask[idx] = 0;
                u_bar[idx] = 0.0;
                v_bar[idx] = 0.0;
                H[idx] = config.dry_depth;
            } else {
                mask[idx] = 1;
            }
        }
    }
}

void HydrodynamicModel::computeVerticalVelocity() {
    if (nz <= 1) return;
    
    // Compute omega (vertical velocity in sigma coordinates) from continuity
    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {
            int idx2 = idx2d(i, j);
            if (mask[idx2] == 0) continue;
            
            // Integrate from bottom
            w[idx3d(i, j, 0)] = 0.0;  // No flow through bottom
            
            for (int k = 1; k < nz; ++k) {
                double dz = dz_at_sigma(i, j, k);
                
                // ∂w/∂z = -∂u/∂x - ∂v/∂y
                int idx = idx3d(i, j, k);
                int idxE = idx3d(i + 1, j, k);
                int idxW = idx3d(i - 1, j, k);
                int idxN = idx3d(i, j + 1, k);
                int idxS = idx3d(i, j - 1, k);
                
                double div = (u[idxE] - u[idxW]) / (2.0 * dx) +
                            (v[idxN] - v[idxS]) / (2.0 * dy);
                
                w[idx] = w[idx3d(i, j, k - 1)] - div * dz;
            }
        }
    }
}

// =============================================================================
// Equation of State
// =============================================================================

double HydrodynamicModel::computeDensity(double T, double S, double p) const {
    if (config.use_nonlinear_eos) {
        return densityUNESCO(T, S, p);
    } else {
        return densityLinear(T, S);
    }
}

double HydrodynamicModel::densityLinear(double T, double S) const {
    // Linear equation of state
    double T0 = 10.0;
    double S0 = 35.0;
    double rho0 = 1025.0;
    double alpha = 2.0e-4;  // Thermal expansion (1/K)
    double beta = 7.5e-4;   // Haline contraction (1/PSU)
    
    return rho0 * (1.0 - alpha * (T - T0) + beta * (S - S0));
}

double HydrodynamicModel::densityUNESCO(double T, double S, double p) const {
    // UNESCO 1983 equation of state (simplified)
    // Full implementation would be much longer
    
    double T2 = T * T;
    double T3 = T2 * T;
    double T4 = T3 * T;
    double T5 = T4 * T;
    
    // Pure water density at atmospheric pressure
    double rho_w = 999.842594 + 6.793952e-2 * T - 9.095290e-3 * T2 
                 + 1.001685e-4 * T3 - 1.120083e-6 * T4 + 6.536332e-9 * T5;
    
    // Salinity contribution
    double S15 = std::sqrt(S * S * S);
    double A = 8.24493e-1 - 4.0899e-3 * T + 7.6438e-5 * T2 
             - 8.2467e-7 * T3 + 5.3875e-9 * T4;
    double B = -5.72466e-3 + 1.0227e-4 * T - 1.6546e-6 * T2;
    double C = 4.8314e-4;
    
    double rho_0 = rho_w + A * S + B * S15 + C * S * S;
    
    // Pressure correction (simplified)
    if (p > 0) {
        double K = 19652.21 + 148.4206 * T - 2.327105 * T2 
                 + 1.360477e-2 * T3 - 5.155288e-5 * T4;
        rho_0 = rho_0 / (1.0 - p / K);
    }
    
    return rho_0;
}

double HydrodynamicModel::densityTEOS10(double T, double S, double p) const {
    // Simplified TEOS-10
    return densityUNESCO(T, S, p);
}

// =============================================================================
// Solution Access
// =============================================================================

double HydrodynamicModel::getElevation(double lon, double lat) const {
    double fi = (lon - lon_min) / (config.lon_max - config.lon_min) * (nx - 1);
    double fj = (lat - lat_min) / (config.lat_max - config.lat_min) * (ny - 1);
    
    int i = std::max(0, std::min(static_cast<int>(fi), nx - 2));
    int j = std::max(0, std::min(static_cast<int>(fj), ny - 2));
    
    double u = fi - i;
    double v = fj - j;
    
    return (1 - u) * (1 - v) * eta[idx2d(i, j)] +
           u * (1 - v) * eta[idx2d(i + 1, j)] +
           (1 - u) * v * eta[idx2d(i, j + 1)] +
           u * v * eta[idx2d(i + 1, j + 1)];
}

void HydrodynamicModel::getDepthAveragedVelocity(std::vector<double>& u_out,
                                                  std::vector<double>& v_out) const {
    u_out = u_bar;
    v_out = v_bar;
}

void HydrodynamicModel::getTotalDepth(std::vector<double>& H_out) const {
    H_out = H;
}

void HydrodynamicModel::run(double end_time) {
    while (current_time < end_time) {
        step();
        
        if (step_count % 100 == 0) {
            std::cout << "Time: " << current_time << " s" << std::endl;
        }
    }
}

void HydrodynamicModel::run(double end_time,
                           std::function<void(double, const HydrodynamicModel&)> callback) {
    double next_output = current_time + config.output_interval;
    
    while (current_time < end_time) {
        step();
        
        if (current_time >= next_output) {
            callback(current_time, *this);
            next_output += config.output_interval;
        }
    }
}

// =============================================================================
// Derived Quantities
// =============================================================================

void HydrodynamicModel::computeStreamfunction(std::vector<double>& psi) const {
    psi.resize(nx * ny, 0.0);
    
    // Integrate from southwest corner
    for (int j = 1; j < ny; ++j) {
        psi[idx2d(0, j)] = psi[idx2d(0, j - 1)] + H[idx2d(0, j)] * u_bar[idx2d(0, j)] * dy;
    }
    
    for (int j = 0; j < ny; ++j) {
        for (int i = 1; i < nx; ++i) {
            psi[idx2d(i, j)] = psi[idx2d(i - 1, j)] - H[idx2d(i, j)] * v_bar[idx2d(i, j)] * dx;
        }
    }
}

void HydrodynamicModel::computeVorticity(std::vector<double>& vort) const {
    vort.resize(nx * ny, 0.0);
    
    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {
            int idx = idx2d(i, j);
            if (mask[idx] == 0) continue;
            
            // ζ = ∂v/∂x - ∂u/∂y
            double dv_dx = (v_bar[idx2d(i + 1, j)] - v_bar[idx2d(i - 1, j)]) / (2.0 * dx);
            double du_dy = (u_bar[idx2d(i, j + 1)] - u_bar[idx2d(i, j - 1)]) / (2.0 * dy);
            vort[idx] = dv_dx - du_dy;
        }
    }
}

double HydrodynamicModel::computeKineticEnergy() const {
    double KE = 0.0;
    
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int idx = idx2d(i, j);
            if (mask[idx] == 0) continue;
            
            double speed2 = u_bar[idx] * u_bar[idx] + v_bar[idx] * v_bar[idx];
            KE += 0.5 * config.rho_0 * H[idx] * speed2 * cell_area(i, j);
        }
    }
    
    return KE;
}

double HydrodynamicModel::computeVolumeTransport(int i_start, int j_start,
                                                  int i_end, int j_end) const {
    double transport = 0.0;
    
    // Transport through section (simple line integral)
    if (i_start == i_end) {
        // Meridional section
        for (int j = j_start; j <= j_end; ++j) {
            int idx = idx2d(i_start, j);
            transport += H[idx] * u_bar[idx] * dy;
        }
    } else {
        // Zonal section
        for (int i = i_start; i <= i_end; ++i) {
            int idx = idx2d(i, j_start);
            transport += H[idx] * v_bar[idx] * dx;
        }
    }
    
    return transport;  // m³/s
}

void HydrodynamicModel::computeMixedLayerDepth(std::vector<double>& mld,
                                                double delta_T) const {
    mld.resize(nx * ny, 0.0);
    
    if (nz <= 1) return;
    
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int idx2 = idx2d(i, j);
            if (mask[idx2] == 0) continue;
            
            double T_surf = temperature[idx3d(i, j, nz - 1)];
            
            // Find depth where T differs from surface by delta_T
            for (int k = nz - 2; k >= 0; --k) {
                double T_k = temperature[idx3d(i, j, k)];
                if (std::abs(T_surf - T_k) > delta_T) {
                    mld[idx2] = -z_at_sigma(i, j, k);
                    break;
                }
            }
            
            if (mld[idx2] == 0) {
                mld[idx2] = H[idx2];  // Mixed to bottom
            }
        }
    }
}

// =============================================================================
// I/O
// =============================================================================

void HydrodynamicModel::writeOutput(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open " << filename << " for writing\n";
        return;
    }
    
    file << "# HydrodynamicModel output at t = " << current_time << " s\n";
    file << "# nx=" << nx << " ny=" << ny << " nz=" << nz << "\n";
    file << "# lon lat h eta u_bar v_bar\n";
    
    file << std::fixed << std::setprecision(6);
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int idx = idx2d(i, j);
            double lon = lon_min + i * (config.lon_max - config.lon_min) / (nx - 1);
            double lat = lat_min + j * (config.lat_max - config.lat_min) / (ny - 1);
            
            file << lon << " " << lat << " " << h[idx] << " " 
                 << eta[idx] << " " << u_bar[idx] << " " << v_bar[idx] << "\n";
        }
    }
}

void HydrodynamicModel::printDiagnostics() const {
    // Compute diagnostics
    double KE = computeKineticEnergy();
    double max_eta = *std::max_element(eta.begin(), eta.end());
    double min_eta = *std::min_element(eta.begin(), eta.end());
    double max_speed = 0.0;
    
    for (int idx = 0; idx < nx * ny; ++idx) {
        if (mask[idx] == 0) continue;
        double speed = std::sqrt(u_bar[idx] * u_bar[idx] + v_bar[idx] * v_bar[idx]);
        max_speed = std::max(max_speed, speed);
    }
    
    std::cout << "==============================\n";
    std::cout << "Time: " << current_time << " s (step " << step_count << ")\n";
    std::cout << "Sea level: min=" << min_eta << " m, max=" << max_eta << " m\n";
    std::cout << "Max speed: " << max_speed << " m/s\n";
    std::cout << "Kinetic energy: " << KE / 1e12 << " TJ\n";
    std::cout << "==============================\n";
}

// =============================================================================
// CoastalHydrodynamicsModel Implementation
// =============================================================================

CoastalHydrodynamicsModel::CoastalHydrodynamicsModel() {}

void CoastalHydrodynamicsModel::initializeStormSurge(const HydrodynamicConfig& cfg) {
    initialize(cfg);
    
    max_eta.resize(nx * ny, -1e10);
    max_depth.resize(nx * ny, 0.0);
    max_speed.resize(nx * ny, 0.0);
}

void CoastalHydrodynamicsModel::setHurricaneTrack(
    const std::vector<double>& time,
    const std::vector<double>& lon,
    const std::vector<double>& lat,
    const std::vector<double>& pressure,
    const std::vector<double>& max_wind,
    const std::vector<double>& rmw) {
    
    track_time = time;
    track_lon = lon;
    track_lat = lat;
    central_pressure = pressure;
    max_wind_speed = max_wind;
    radius_max_wind = rmw;
}

void CoastalHydrodynamicsModel::computeHollandWinds(double t,
                                                    std::vector<double>& u10,
                                                    std::vector<double>& v10) {
    // Interpolate hurricane position and intensity
    size_t idx = 0;
    double alpha = 0.0;
    
    for (size_t i = 0; i < track_time.size() - 1; ++i) {
        if (t >= track_time[i] && t <= track_time[i + 1]) {
            idx = i;
            alpha = (t - track_time[i]) / (track_time[i + 1] - track_time[i]);
            break;
        }
    }
    
    double center_lon = (1 - alpha) * track_lon[idx] + alpha * track_lon[idx + 1];
    double center_lat = (1 - alpha) * track_lat[idx] + alpha * track_lat[idx + 1];
    double Pc = (1 - alpha) * central_pressure[idx] + alpha * central_pressure[idx + 1];
    double Vmax = (1 - alpha) * max_wind_speed[idx] + alpha * max_wind_speed[idx + 1];
    double Rmax = (1 - alpha) * radius_max_wind[idx] + alpha * radius_max_wind[idx + 1];
    
    double dP = ambient_pressure - Pc;
    
    // Holland B parameter
    double B = Vmax * Vmax * RHO_AIR * std::exp(1.0) / dP;
    B = std::max(1.0, std::min(B, 2.5));
    
    u10.resize(nx * ny);
    v10.resize(nx * ny);
    
    double f = 2.0 * OMEGA * std::sin(center_lat * DEG_TO_RAD);
    
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int idx2 = idx2d(i, j);
            
            double lon = lon_min + i * (config.lon_max - config.lon_min) / (nx - 1);
            double lat = lat_min + j * (config.lat_max - config.lat_min) / (ny - 1);
            
            // Distance from storm center
            double dx_dist = (lon - center_lon) * DEG_TO_RAD * EARTH_RADIUS * std::cos(lat * DEG_TO_RAD);
            double dy_dist = (lat - center_lat) * DEG_TO_RAD * EARTH_RADIUS;
            double r = std::sqrt(dx_dist * dx_dist + dy_dist * dy_dist);
            
            if (r < 1000) r = 1000;  // Avoid singularity at center
            
            // Holland (1980) wind profile
            double rr = Rmax * 1000.0 / r;
            double V = std::sqrt(B * dP / RHO_AIR * std::pow(rr, B) * std::exp(1.0 - std::pow(rr, B)) +
                                r * r * f * f / 4.0) - r * f / 2.0;
            
            // Inflow angle (cross-isobar angle)
            double theta = std::atan2(dy_dist, dx_dist);
            double inflow_angle = 20.0 * DEG_TO_RAD;  // Typical value
            
            // Wind components (cyclonic rotation in Northern Hemisphere)
            u10[idx2] = -V * std::sin(theta + inflow_angle);
            v10[idx2] = V * std::cos(theta + inflow_angle);
        }
    }
}

void CoastalHydrodynamicsModel::computeHollandPressure(double t, std::vector<double>& p) {
    // Similar interpolation as winds
    size_t idx = 0;
    double alpha = 0.0;
    
    for (size_t i = 0; i < track_time.size() - 1; ++i) {
        if (t >= track_time[i] && t <= track_time[i + 1]) {
            idx = i;
            alpha = (t - track_time[i]) / (track_time[i + 1] - track_time[i]);
            break;
        }
    }
    
    double center_lon = (1 - alpha) * track_lon[idx] + alpha * track_lon[idx + 1];
    double center_lat = (1 - alpha) * track_lat[idx] + alpha * track_lat[idx + 1];
    double Pc = (1 - alpha) * central_pressure[idx] + alpha * central_pressure[idx + 1];
    double Rmax = (1 - alpha) * radius_max_wind[idx] + alpha * radius_max_wind[idx + 1];
    
    double dP = ambient_pressure - Pc;
    double B = holland_B;
    
    p.resize(nx * ny);
    
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int idx2 = idx2d(i, j);
            
            double lon = lon_min + i * (config.lon_max - config.lon_min) / (nx - 1);
            double lat = lat_min + j * (config.lat_max - config.lat_min) / (ny - 1);
            
            double dx_dist = (lon - center_lon) * DEG_TO_RAD * EARTH_RADIUS * std::cos(lat * DEG_TO_RAD);
            double dy_dist = (lat - center_lat) * DEG_TO_RAD * EARTH_RADIUS;
            double r = std::sqrt(dx_dist * dx_dist + dy_dist * dy_dist);
            
            // Holland pressure profile
            double rr = Rmax * 1000.0 / std::max(r, 1000.0);
            p[idx2] = Pc + dP * std::exp(-std::pow(rr, B));
        }
    }
}

void CoastalHydrodynamicsModel::computeMaxInundation(
    std::vector<double>& max_eta_out,
    std::vector<double>& max_depth_out,
    std::vector<double>& max_velocity_out) {
    
    max_eta_out = max_eta;
    max_depth_out = max_depth;
    max_velocity_out = max_speed;
}

// =============================================================================
// EstuarineModel Implementation
// =============================================================================

EstuarineModel::EstuarineModel() 
    : river_discharge(100.0), ocean_salinity(35.0), ocean_temperature(15.0) {}

void EstuarineModel::initializeEstuary(const HydrodynamicConfig& cfg,
                                       double river_width,
                                       double river_depth,
                                       double ocean_sal) {
    (void)river_width; (void)river_depth;  // TODO: Use for mesh refinement near river mouth
    this->ocean_salinity = ocean_sal;
    initialize(cfg);
    this->ocean_salinity = ocean_salinity;
    
    // Initialize salinity field (linear gradient as first guess)
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            for (int k = 0; k < nz; ++k) {
                double x_frac = static_cast<double>(i) / (nx - 1);
                salinity[idx3d(i, j, k)] = ocean_salinity * x_frac;
            }
        }
    }
}

double EstuarineModel::computeSaltIntrusionLength(double salinity_threshold) {
    // Find the upstream extent where bottom salinity > threshold
    for (int i = 0; i < nx; ++i) {
        int j = ny / 2;  // Mid-channel
        double S_bottom = salinity[idx3d(i, j, 0)];
        
        if (S_bottom < salinity_threshold) {
            // Intrusion ends here
            return i * dx / 1000.0;  // km
        }
    }
    
    return nx * dx / 1000.0;  // Full length
}

double EstuarineModel::computeEstuarineRichardson() {
    // Ri_E = g β ΔS H / u²
    // where β ≈ 0.0008 is haline contraction coefficient
    
    double beta_s = 0.0008;
    double delta_S = ocean_salinity;
    double mean_H = 0.0;
    double mean_u2 = 0.0;
    int count = 0;
    
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int idx = idx2d(i, j);
            if (mask[idx] == 0) continue;
            
            mean_H += H[idx];
            mean_u2 += u_bar[idx] * u_bar[idx];
            count++;
        }
    }
    
    if (count == 0 || mean_u2 < 1e-10) return 1e10;
    
    mean_H /= count;
    mean_u2 /= count;
    
    return G * beta_s * delta_S * mean_H / mean_u2;
}

EstuarineModel::EstuaryType EstuarineModel::classifyEstuary() {
    double Ri = computeEstuarineRichardson();
    
    if (Ri > 0.8) {
        return SALT_WEDGE;
    } else if (Ri > 0.08) {
        return PARTIALLY_MIXED;
    } else {
        return WELL_MIXED;
    }
}

double EstuarineModel::computeFlushingTime() {
    // Simple tidal prism method
    double V_total = 0.0;
    double V_tidal = 0.0;
    double mean_eta_range = 2.0;  // Assume 2m tidal range
    
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int idx = idx2d(i, j);
            if (mask[idx] == 0) continue;
            
            V_total += H[idx] * cell_area(i, idx);
            V_tidal += mean_eta_range * cell_area(i, idx);
        }
    }
    
    // Flushing time ≈ V / (P * r) where r is return flow factor
    double r = 0.7;
    double T_tide = 12.42 * 3600.0;  // M2 period in seconds
    
    return V_total / (V_tidal * r) * T_tide;
}

// =============================================================================
// HydrodynamicUtils Implementation
// =============================================================================

namespace HydrodynamicUtils {

std::map<TidalConstituent, std::pair<double, double>>
tidalHarmonicAnalysis(const std::vector<double>& time,
                     const std::vector<double>& eta,
                     const std::vector<TidalConstituent>& constituents) {
    std::map<TidalConstituent, std::pair<double, double>> results;
    
    size_t n = time.size();
    if (n < 2) return results;
    
    for (const auto& type : constituents) {
        double period, omega;
        getConstituentProperties(type, period, omega);
        
        // Least squares harmonic analysis
        double A = 0.0, B = 0.0;
        
        for (size_t i = 0; i < n; ++i) {
            double arg = omega * time[i];
            A += eta[i] * std::cos(arg);
            B += eta[i] * std::sin(arg);
        }
        
        A *= 2.0 / n;
        B *= 2.0 / n;
        
        double amplitude = std::sqrt(A * A + B * B);
        double phase = std::atan2(B, A) * RAD_TO_DEG;
        if (phase < 0) phase += 360.0;
        
        results[type] = {amplitude, phase};
    }
    
    return results;
}

void getConstituentProperties(TidalConstituent type,
                             double& period_hours,
                             double& frequency) {
    switch (type) {
        case TidalConstituent::M2:
            period_hours = 12.4206;
            break;
        case TidalConstituent::S2:
            period_hours = 12.0;
            break;
        case TidalConstituent::N2:
            period_hours = 12.6583;
            break;
        case TidalConstituent::K2:
            period_hours = 11.9672;
            break;
        case TidalConstituent::K1:
            period_hours = 23.9345;
            break;
        case TidalConstituent::O1:
            period_hours = 25.8193;
            break;
        case TidalConstituent::P1:
            period_hours = 24.0659;
            break;
        case TidalConstituent::Q1:
            period_hours = 26.8684;
            break;
        case TidalConstituent::M4:
            period_hours = 6.2103;
            break;
        case TidalConstituent::MS4:
            period_hours = 6.1033;
            break;
        default:
            period_hours = 12.0;
    }
    
    frequency = 2.0 * M_PI / (period_hours * 3600.0);
}

double manningToDragCoefficient(double n, double h) {
    return G * n * n / std::pow(h, 1.0/3.0);
}

double dragCoefficientToManning(double Cd, double h) {
    return std::sqrt(Cd * std::pow(h, 1.0/3.0) / G);
}

BathymetryGrid generateEstuaryBathymetry(double length, double width,
                                         double max_depth,
                                         double river_depth) {
    BathymetryGrid grid;
    
    int nx = 200;
    int ny = 50;
    
    grid.nlon = nx;
    grid.nlat = ny;
    grid.lon_min = 0.0;
    grid.lon_max = length / 111000.0;  // Convert m to degrees
    grid.lat_min = 0.0;
    grid.lat_max = width / 111000.0;
    grid.dlon = (grid.lon_max - grid.lon_min) / (nx - 1);
    grid.dlat = (grid.lat_max - grid.lat_min) / (ny - 1);
    
    grid.depth.resize(nx * ny);
    
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            // Cross-channel profile (parabolic)
            double y_frac = 2.0 * (j - ny / 2.0) / ny;
            double cross_factor = 1.0 - y_frac * y_frac;
            
            // Along-channel profile (increasing toward mouth)
            double x_frac = static_cast<double>(i) / (nx - 1);
            double along_depth = river_depth + (max_depth - river_depth) * x_frac;
            
            grid.depth[i + j * nx] = along_depth * cross_factor;
            
            // Land at edges
            if (cross_factor < 0.1) {
                grid.depth[i + j * nx] = -5.0;  // Land
            }
        }
    }
    
    return grid;
}

} // namespace HydrodynamicUtils

} // namespace FSRM
