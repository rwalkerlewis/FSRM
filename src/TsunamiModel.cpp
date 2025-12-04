/**
 * @file TsunamiModel.cpp
 * @brief Implementation of tsunami modeling with oceanic coupling
 */

#include "TsunamiModel.hpp"
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <numeric>

namespace FSRM {

// Physical constants
static const double G = 9.81;              // Gravitational acceleration (m/s²)
static const double EARTH_RADIUS = 6371000.0; // Earth radius (m)
static const double DEG_TO_RAD = M_PI / 180.0;
static const double RAD_TO_DEG = 180.0 / M_PI;

// =============================================================================
// CascadiaFaultModel Implementation
// =============================================================================

std::vector<TsunamiSubfault> CascadiaFaultModel::generateSubfaults() const {
    std::vector<TsunamiSubfault> subfaults;
    
    // Total along-strike length (km)
    double total_length = TsunamiUtils::haversineDistance(
        trench_longitude, south_latitude,
        trench_longitude, north_latitude) / 1000.0;
    
    // Subfault dimensions
    double subfault_length = total_length / num_along_strike;
    double subfault_width = rupture_width / num_down_dip;
    
    // Latitude step
    double lat_step = (north_latitude - south_latitude) / num_along_strike;
    
    // Generate subfaults
    for (int i = 0; i < num_along_strike; ++i) {
        for (int j = 0; j < num_down_dip; ++j) {
            TsunamiSubfault sf;
            
            // Location (center of subfault)
            sf.latitude = south_latitude + (i + 0.5) * lat_step;
            
            // Compute dip at this depth
            double depth_to_center = (j + 0.5) * subfault_width * std::sin(deep_dip * DEG_TO_RAD);
            double dip = (depth_to_center < transition_depth) ? shallow_dip : deep_dip;
            sf.dip = dip;
            
            // Horizontal distance from trench
            double horiz_dist = 0;
            for (int k = 0; k <= j; ++k) {
                double d = (k < j) ? shallow_dip : dip;
                horiz_dist += subfault_width * std::cos(d * DEG_TO_RAD);
            }
            
            // Longitude (east of trench)
            sf.longitude = trench_longitude + horiz_dist / 111.0; // Approximate km to degrees
            
            // Depth to top of subfault
            sf.depth = 0;
            for (int k = 0; k < j; ++k) {
                double d = (k == 0) ? shallow_dip : deep_dip;
                sf.depth += subfault_width * std::sin(d * DEG_TO_RAD);
            }
            
            sf.strike = average_strike;
            sf.rake = 90.0;  // Pure thrust
            sf.length = subfault_length;
            sf.width = subfault_width;
            
            // Slip distribution - higher near trench
            double depth_factor = 1.0 - 0.5 * (j / static_cast<double>(num_down_dip));
            double along_strike_factor = 1.0 + 0.3 * std::sin(M_PI * i / num_along_strike);
            sf.slip = average_slip * depth_factor * along_strike_factor;
            
            // Limit to peak slip
            sf.slip = std::min(sf.slip, peak_slip);
            
            // Rupture time from hypocenter
            double dist_from_hypo = TsunamiUtils::haversineDistance(
                sf.longitude, sf.latitude,
                hypocenter_lon, hypocenter_lat) / 1000.0;
            sf.rupture_time = dist_from_hypo / rupture_velocity;
            
            // Rise time scales with slip
            sf.rise_time = 5.0 + 20.0 * (sf.slip / peak_slip);
            
            subfaults.push_back(sf);
        }
    }
    
    return subfaults;
}

double CascadiaFaultModel::getTotalMoment() const {
    auto subfaults = generateSubfaults();
    double moment = 0.0;
    double mu = 30.0e9;  // Shear modulus (Pa)
    
    for (const auto& sf : subfaults) {
        double area = sf.length * 1000.0 * sf.width * 1000.0;  // m²
        moment += mu * area * sf.slip;
    }
    
    return moment;
}

double CascadiaFaultModel::getMagnitude() const {
    double M0 = getTotalMoment();
    return (2.0/3.0) * (std::log10(M0) - 9.1);
}

// =============================================================================
// OkadaModel Implementation
// =============================================================================

OkadaModel::OkadaModel() 
    : shear_modulus(30.0e9), poisson_ratio(0.25) {}

void OkadaModel::geoToLocal(double lon, double lat, double lon0, double lat0,
                            double strike, double& x, double& y) const {
    // Convert to local Cartesian coordinates relative to fault center
    double dlon = (lon - lon0) * DEG_TO_RAD;
    double dlat = (lat - lat0) * DEG_TO_RAD;
    
    // Approximate conversion to meters
    double lat_rad = lat0 * DEG_TO_RAD;
    double x_raw = EARTH_RADIUS * dlon * std::cos(lat_rad);
    double y_raw = EARTH_RADIUS * dlat;
    
    // Rotate to fault-aligned coordinate system
    double strike_rad = strike * DEG_TO_RAD;
    double cos_s = std::cos(strike_rad);
    double sin_s = std::sin(strike_rad);
    
    x = x_raw * sin_s + y_raw * cos_s;  // Along-strike
    y = -x_raw * cos_s + y_raw * sin_s; // Perpendicular to strike
}

void OkadaModel::okadaGreenFunction(double x, double y, double d,
                                    double dip, double L, double W, double slip,
                                    double& ux, double& uy, double& uz) const {
    // Okada (1985) analytical solution for surface displacement
    // from a rectangular dislocation in an elastic half-space
    
    double dip_rad = dip * DEG_TO_RAD;
    double sin_d = std::sin(dip_rad);
    double cos_d = std::cos(dip_rad);
    
    // Lamé parameters
    double lambda = 2.0 * shear_modulus * poisson_ratio / (1.0 - 2.0 * poisson_ratio);
    double alpha = (lambda + shear_modulus) / (lambda + 2.0 * shear_modulus);
    
    // Compute at corners of fault
    double p = y * cos_d + d * sin_d;
    double q = y * sin_d - d * cos_d;
    
    auto chinnery = [&](double xi, double eta) -> std::array<double, 3> {
        double R = std::sqrt(xi * xi + eta * eta + q * q);
        double X = std::sqrt(xi * xi + q * q);
        
        // Avoid singularities
        double eps = 1e-6;
        if (R < eps) R = eps;
        if (X < eps) X = eps;
        
        double R_eta = R + eta;
        if (std::abs(R_eta) < eps) R_eta = eps;
        
        // Displacement components for strike-slip (here modified for dip-slip)
        double I1 = 0, I3 = 0, I4 = 0, I5 = 0;
        
        if (std::abs(cos_d) > eps) {
            I5 = alpha * 2.0 / cos_d * std::atan(
                (eta * (X + q * cos_d) + X * (R + X) * sin_d) /
                (xi * (R + X) * cos_d));
            I4 = alpha / cos_d * (std::log(R + eta) - sin_d * std::log(R + q));
            I3 = alpha * (1.0 / cos_d * eta / (R + eta) - I4);
            I1 = alpha * (-1.0 / cos_d * xi / (R + eta)) - tan(dip_rad) * I5;
        } else {
            // Vertical fault
            I1 = -alpha / 2.0 * xi * q / (R + eta) / (R + eta);
            I3 = alpha / 2.0 * (eta / (R + eta) + eta * R / (R + eta) / (R + eta) 
                 - std::log(R + eta));
            I4 = -alpha * q / (R + eta);
            I5 = -alpha * xi * sin_d / (R + eta);
        }
        
        double I2 = alpha * (-std::log(R + eta)) - I3;
        
        // Dip-slip displacement (rake = 90)
        double u_xi = -slip / (2.0 * M_PI) * (
            q / R - I3 * sin_d * cos_d);
        double u_eta = -slip / (2.0 * M_PI) * (
            eta * q / (R * (R + xi)) + cos_d * std::atan(xi * eta / (q * R)) 
            - I1 * sin_d * cos_d);
        double u_q = -slip / (2.0 * M_PI) * (
            q * q / (R * (R + xi)) - I4 * sin_d * cos_d);
        
        return {u_xi, u_eta, u_q};
    };
    
    // Chinnery's notation: evaluate at four corners
    auto c1 = chinnery(x, p);
    auto c2 = chinnery(x, p - W);
    auto c3 = chinnery(x - L, p);
    auto c4 = chinnery(x - L, p - W);
    
    // Apply Chinnery's relation
    std::array<double, 3> u_local;
    for (int i = 0; i < 3; ++i) {
        u_local[i] = c1[i] - c2[i] - c3[i] + c4[i];
    }
    
    // Transform to geographic coordinates
    // (simplified - assumes strike is approximately N-S)
    ux = u_local[0];
    uy = u_local[1] * cos_d + u_local[2] * sin_d;
    uz = -u_local[1] * sin_d + u_local[2] * cos_d;
}

void OkadaModel::computeDisplacement(double lon, double lat,
                                     const TsunamiSubfault& subfault,
                                     double& ux, double& uy, double& uz) const {
    // Convert observation point to local coordinates
    double x, y;
    geoToLocal(lon, lat, subfault.longitude, subfault.latitude,
               subfault.strike, x, y);
    
    // Fault dimensions in meters
    double L = subfault.length * 1000.0;
    double W = subfault.width * 1000.0;
    double d = subfault.depth * 1000.0;
    
    // Compute using Okada Green's function
    okadaGreenFunction(x, y, d, subfault.dip, L, W, subfault.slip, ux, uy, uz);
}

void OkadaModel::computeTotalDisplacement(double lon, double lat,
                                          const std::vector<TsunamiSubfault>& subfaults,
                                          double& ux, double& uy, double& uz) const {
    ux = uy = uz = 0.0;
    
    for (const auto& sf : subfaults) {
        double dux, duy, duz;
        computeDisplacement(lon, lat, sf, dux, duy, duz);
        ux += dux;
        uy += duy;
        uz += duz;
    }
}

void OkadaModel::computeDisplacementField(const BathymetryGrid& grid,
                                          const std::vector<TsunamiSubfault>& subfaults,
                                          std::vector<double>& uz_field) const {
    int n = grid.nlon * grid.nlat;
    uz_field.resize(n, 0.0);
    
    #pragma omp parallel for
    for (int j = 0; j < grid.nlat; ++j) {
        for (int i = 0; i < grid.nlon; ++i) {
            double lon = grid.lon_min + i * grid.dlon;
            double lat = grid.lat_min + j * grid.dlat;
            
            double ux, uy, uz;
            computeTotalDisplacement(lon, lat, subfaults, ux, uy, uz);
            
            uz_field[i + j * grid.nlon] = uz;
        }
    }
}

void OkadaModel::applyKajiuraFilter(const BathymetryGrid& grid,
                                    std::vector<double>& uz_field) const {
    // Kajiura (1963) filter accounts for water column smoothing
    // Transfer function: H(k) = 1/cosh(kh)
    // This is a simplified spatial-domain approximation
    
    int nx = grid.nlon;
    int ny = grid.nlat;
    
    std::vector<double> filtered(uz_field.size());
    
    // Characteristic filter width based on local water depth
    double dx = grid.dlon * EARTH_RADIUS * DEG_TO_RAD * std::cos(0.5 * (grid.lat_min + grid.lat_max) * DEG_TO_RAD);
    double dy = grid.dlat * EARTH_RADIUS * DEG_TO_RAD;
    
    #pragma omp parallel for
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            double h = std::max(10.0, grid.depth[i + j * nx]);  // Water depth
            double filter_width = h / 3.0;  // Approximate filter width
            
            int r = static_cast<int>(filter_width / std::min(dx, dy));
            r = std::max(1, std::min(r, 10));  // Limit filter size
            
            double sum_w = 0;
            double sum_uz = 0;
            
            for (int jj = std::max(0, j - r); jj <= std::min(ny - 1, j + r); ++jj) {
                for (int ii = std::max(0, i - r); ii <= std::min(nx - 1, i + r); ++ii) {
                    double dist = std::sqrt((ii - i) * (ii - i) * dx * dx + 
                                           (jj - j) * (jj - j) * dy * dy);
                    double w = std::exp(-dist * dist / (2.0 * filter_width * filter_width));
                    sum_w += w;
                    sum_uz += w * uz_field[ii + jj * nx];
                }
            }
            
            filtered[i + j * nx] = sum_uz / sum_w;
        }
    }
    
    uz_field = filtered;
}

// =============================================================================
// BathymetryGrid Implementation
// =============================================================================

double BathymetryGrid::getDepth(double lon, double lat) const {
    // Bilinear interpolation
    double fi = (lon - lon_min) / dlon;
    double fj = (lat - lat_min) / dlat;
    
    int i = static_cast<int>(fi);
    int j = static_cast<int>(fj);
    
    // Clamp to valid range
    i = std::max(0, std::min(i, nlon - 2));
    j = std::max(0, std::min(j, nlat - 2));
    
    double u = fi - i;
    double v = fj - j;
    
    double d00 = depth[i + j * nlon];
    double d10 = depth[(i + 1) + j * nlon];
    double d01 = depth[i + (j + 1) * nlon];
    double d11 = depth[(i + 1) + (j + 1) * nlon];
    
    return (1 - u) * (1 - v) * d00 + u * (1 - v) * d10 +
           (1 - u) * v * d01 + u * v * d11;
}

double BathymetryGrid::getManningN(double lon, double lat) const {
    if (manning_n.empty()) return 0.025;  // Default ocean value
    
    double fi = (lon - lon_min) / dlon;
    double fj = (lat - lat_min) / dlat;
    
    int i = static_cast<int>(fi);
    int j = static_cast<int>(fj);
    
    i = std::max(0, std::min(i, nlon - 1));
    j = std::max(0, std::min(j, nlat - 1));
    
    return manning_n[i + j * nlon];
}

void BathymetryGrid::getGradient(double lon, double lat, 
                                  double& dzdlon, double& dzdlat) const {
    double h = 0.001;  // Small step in degrees
    
    double z_plus_lon = getDepth(lon + h, lat);
    double z_minus_lon = getDepth(lon - h, lat);
    double z_plus_lat = getDepth(lon, lat + h);
    double z_minus_lat = getDepth(lon, lat - h);
    
    dzdlon = (z_plus_lon - z_minus_lon) / (2.0 * h);
    dzdlat = (z_plus_lat - z_minus_lat) / (2.0 * h);
}

void BathymetryGrid::generateContinentalShelf(double shelf_width, double shelf_depth,
                                               double slope_angle, double deep_depth) {
    depth.resize(nlon * nlat);
    
    double slope_width = (deep_depth - shelf_depth) / std::tan(slope_angle * DEG_TO_RAD);
    
    for (int j = 0; j < nlat; ++j) {
        for (int i = 0; i < nlon; ++i) {
            // Distance from western edge (coast)
            double dist = (i / static_cast<double>(nlon - 1)) * 
                         (lon_max - lon_min) * 111.0;  // km (approximate)
            
            double d;
            if (dist < 0) {
                d = -10.0;  // Land (10 m elevation)
            } else if (dist < shelf_width) {
                d = shelf_depth * (dist / shelf_width);
            } else if (dist < shelf_width + slope_width) {
                double frac = (dist - shelf_width) / slope_width;
                d = shelf_depth + frac * (deep_depth - shelf_depth);
            } else {
                d = deep_depth;
            }
            
            depth[i + j * nlon] = d;
        }
    }
}

bool BathymetryGrid::loadFromASCII(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) return false;
    
    // Read header
    std::string line;
    std::map<std::string, double> header;
    
    while (std::getline(file, line) && !line.empty()) {
        std::istringstream iss(line);
        std::string key;
        double value;
        if (iss >> key >> value) {
            header[key] = value;
        }
        if (header.size() >= 6) break;
    }
    
    nlon = static_cast<int>(header["ncols"]);
    nlat = static_cast<int>(header["nrows"]);
    lon_min = header["xllcorner"];
    lat_min = header["yllcorner"];
    dlon = header["cellsize"];
    dlat = header["cellsize"];
    lon_max = lon_min + nlon * dlon;
    lat_max = lat_min + nlat * dlat;
    
    // Read data
    depth.resize(nlon * nlat);
    for (int j = nlat - 1; j >= 0; --j) {  // ASCII grids are typically top-to-bottom
        for (int i = 0; i < nlon; ++i) {
            file >> depth[i + j * nlon];
        }
    }
    
    return true;
}

// =============================================================================
// TsunamiGauge Implementation
// =============================================================================

void TsunamiGauge::writeASCII(const std::string& filename) const {
    std::ofstream file(filename);
    file << "# Tsunami gauge: " << name << " (" << id << ")\n";
    file << "# Location: " << longitude << " E, " << latitude << " N\n";
    file << "# Time (s), Eta (m)";
    if (!u.empty()) file << ", U (m/s), V (m/s)";
    file << "\n";
    
    file << std::fixed << std::setprecision(4);
    for (size_t i = 0; i < time.size(); ++i) {
        file << time[i] << ", " << eta[i];
        if (!u.empty()) {
            file << ", " << u[i] << ", " << v[i];
        }
        file << "\n";
    }
}

double TsunamiGauge::getMaxAmplitude() const {
    if (eta.empty()) return 0.0;
    return *std::max_element(eta.begin(), eta.end());
}

double TsunamiGauge::getFirstArrivalTime(double threshold) const {
    for (size_t i = 0; i < eta.size(); ++i) {
        if (std::abs(eta[i]) > threshold) {
            return time[i];
        }
    }
    return -1.0;  // No arrival
}

// =============================================================================
// ShallowWaterSolver Implementation
// =============================================================================

ShallowWaterSolver::ShallowWaterSolver() 
    : current_time(0.0), dt(0.5) {}

ShallowWaterSolver::~ShallowWaterSolver() = default;

void ShallowWaterSolver::initialize(const BathymetryGrid& bathymetry,
                                    const TsunamiConfig& cfg) {
    config = cfg;
    
    // Grid setup
    double lat_avg = 0.5 * (cfg.lat_min + cfg.lat_max);
    dx = cfg.base_resolution * EARTH_RADIUS * DEG_TO_RAD * std::cos(lat_avg * DEG_TO_RAD);
    dy = cfg.base_resolution * EARTH_RADIUS * DEG_TO_RAD;
    
    nx = static_cast<int>((cfg.lon_max - cfg.lon_min) / cfg.base_resolution);
    ny = static_cast<int>((cfg.lat_max - cfg.lat_min) / cfg.base_resolution);
    
    lon_min = cfg.lon_min;
    lat_min = cfg.lat_min;
    
    // Allocate arrays
    int n = nx * ny;
    depth.resize(n);
    manning_n.resize(n);
    h.resize(n, 0.0);
    hu.resize(n, 0.0);
    hv.resize(n, 0.0);
    eta.resize(n, 0.0);
    
    max_eta.resize(n, -1e10);
    max_speed.resize(n, 0.0);
    max_momentum_flux.resize(n, 0.0);
    arrival_time.resize(n, -1.0);
    
    // Copy bathymetry (interpolated to our grid)
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            double lon = lon_min + i * cfg.base_resolution;
            double lat = lat_min + j * cfg.base_resolution;
            depth[idx(i, j)] = bathymetry.getDepth(lon, lat);
            manning_n[idx(i, j)] = bathymetry.getManningN(lon, lat);
        }
    }
    
    current_time = 0.0;
    dt = cfg.dt;
}

void ShallowWaterSolver::setInitialCondition(const std::vector<double>& eta_initial) {
    if (eta_initial.size() != eta.size()) {
        std::cerr << "Error: Initial condition size mismatch\n";
        return;
    }
    
    for (size_t i = 0; i < eta.size(); ++i) {
        eta[i] = eta_initial[i];
        // h = depth + eta for wet cells
        double d = depth[i];
        if (d > 0) {  // Ocean
            h[i] = d + eta[i];
        } else {  // Land
            h[i] = std::max(0.0, -d + eta[i]);
        }
        hu[i] = 0.0;
        hv[i] = 0.0;
    }
}

void ShallowWaterSolver::setSeafloorMotion(
    std::function<void(double, std::vector<double>&)> motion_func) {
    seafloor_motion = motion_func;
}

void ShallowWaterSolver::riemannSolverHLLC(double hL, double huL, double hvL,
                                           double hR, double huR, double hvR,
                                           double& f_h, double& f_hu, double& f_hv) {
    // HLLC approximate Riemann solver
    double eps = 1e-10;
    
    // Handle dry states
    if (hL < eps && hR < eps) {
        f_h = f_hu = f_hv = 0.0;
        return;
    }
    
    double uL = (hL > eps) ? huL / hL : 0.0;
    double uR = (hR > eps) ? huR / hR : 0.0;
    double vL = (hL > eps) ? hvL / hL : 0.0;
    double vR = (hR > eps) ? hvR / hR : 0.0;
    
    double aL = std::sqrt(G * std::max(eps, hL));
    double aR = std::sqrt(G * std::max(eps, hR));
    
    // Wave speed estimates (Einfeldt)
    double hBar = 0.5 * (hL + hR);
    double uBar = (aL * uL + aR * uR) / (aL + aR);
    double aBar = std::sqrt(G * hBar);
    
    double SL = std::min(uL - aL, uBar - aBar);
    double SR = std::max(uR + aR, uBar + aBar);
    
    // Star state
    double SM = (SL * hR * (uR - SR) - SR * hL * (uL - SL)) /
                (hR * (uR - SR) - hL * (uL - SL) + eps);
    
    // Flux calculation
    auto flux = [&](double h_, double hu_, double hv_) -> std::array<double, 3> {
        double u_ = (h_ > eps) ? hu_ / h_ : 0.0;
        return {
            hu_,
            hu_ * u_ + 0.5 * G * h_ * h_,
            hv_ * u_
        };
    };
    
    auto FL = flux(hL, huL, hvL);
    auto FR = flux(hR, huR, hvR);
    
    if (SL >= 0) {
        f_h = FL[0];
        f_hu = FL[1];
        f_hv = FL[2];
    } else if (SR <= 0) {
        f_h = FR[0];
        f_hu = FR[1];
        f_hv = FR[2];
    } else if (SM >= 0) {
        double hSL = hL * (SL - uL) / (SL - SM);
        f_h = FL[0] + SL * (hSL - hL);
        f_hu = FL[1] + SL * (hSL * SM - huL);
        f_hv = FL[2] + SL * (hSL * vL - hvL);
    } else {
        double hSR = hR * (SR - uR) / (SR - SM);
        f_h = FR[0] + SR * (hSR - hR);
        f_hu = FR[1] + SR * (hSR * SM - huR);
        f_hv = FR[2] + SR * (hSR * vR - hvR);
    }
}

double ShallowWaterSolver::computeCFL() const {
    double max_speed_local = 0.0;
    
    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {
            int k = idx(i, j);
            if (h[k] < config.dry_tolerance) continue;
            
            double u_local = hu[k] / h[k];
            double v_local = hv[k] / h[k];
            double c = std::sqrt(G * h[k]);
            
            double speed = std::sqrt(u_local * u_local + v_local * v_local) + c;
            max_speed_local = std::max(max_speed_local, speed);
        }
    }
    
    if (max_speed_local > 0) {
        return config.cfl_number * std::min(dx, dy) / max_speed_local;
    }
    return config.dt;
}

void ShallowWaterSolver::applyBottomFriction(double dt_local) {
    if (config.friction_model == BottomFrictionModel::NONE) return;
    
    #pragma omp parallel for
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int k = idx(i, j);
            if (h[k] < config.dry_tolerance) continue;
            
            double u_local = hu[k] / h[k];
            double v_local = hv[k] / h[k];
            double speed = std::sqrt(u_local * u_local + v_local * v_local);
            
            if (speed < 1e-10) continue;
            
            double n = manning_n[k];
            double Cf = G * n * n / std::pow(h[k], 1.0/3.0);
            
            // Semi-implicit friction
            double factor = 1.0 / (1.0 + dt_local * Cf * speed / h[k]);
            hu[k] *= factor;
            hv[k] *= factor;
        }
    }
}

void ShallowWaterSolver::applyWettingDrying() {
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int k = idx(i, j);
            
            if (h[k] < config.dry_tolerance) {
                h[k] = 0.0;
                hu[k] = 0.0;
                hv[k] = 0.0;
            }
        }
    }
}

void ShallowWaterSolver::applyBoundaryConditions() {
    // Open boundaries (radiation condition)
    for (int j = 0; j < ny; ++j) {
        // Left boundary
        h[idx(0, j)] = h[idx(1, j)];
        hu[idx(0, j)] = hu[idx(1, j)];
        hv[idx(0, j)] = hv[idx(1, j)];
        
        // Right boundary
        h[idx(nx-1, j)] = h[idx(nx-2, j)];
        hu[idx(nx-1, j)] = hu[idx(nx-2, j)];
        hv[idx(nx-1, j)] = hv[idx(nx-2, j)];
    }
    
    for (int i = 0; i < nx; ++i) {
        // Bottom boundary
        h[idx(i, 0)] = h[idx(i, 1)];
        hu[idx(i, 0)] = hu[idx(i, 1)];
        hv[idx(i, 0)] = hv[idx(i, 1)];
        
        // Top boundary
        h[idx(i, ny-1)] = h[idx(i, ny-2)];
        hu[idx(i, ny-1)] = hu[idx(i, ny-2)];
        hv[idx(i, ny-1)] = hv[idx(i, ny-2)];
    }
}

double ShallowWaterSolver::step() {
    // Adaptive time step
    double dt_local = config.adaptive_timestep ? computeCFL() : config.dt;
    dt_local = std::min(dt_local, config.dt);
    
    int n = nx * ny;
    std::vector<double> h_new(n), hu_new(n), hv_new(n);
    
    // Copy current state
    h_new = h;
    hu_new = hu;
    hv_new = hv;
    
    // Seafloor motion source
    if (seafloor_motion) {
        std::vector<double> dz(n);
        seafloor_motion(current_time + dt_local, dz);
        
        // Add to water depth
        for (int i = 0; i < n; ++i) {
            h_new[i] += dz[i];
        }
    }
    
    // X-direction fluxes
    #pragma omp parallel for
    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {
            int k = idx(i, j);
            int kL = idx(i - 1, j);
            int kR = idx(i + 1, j);
            
            // Left interface flux
            double f_h_L, f_hu_L, f_hv_L;
            riemannSolverHLLC(h[kL], hu[kL], hv[kL],
                             h[k], hu[k], hv[k],
                             f_h_L, f_hu_L, f_hv_L);
            
            // Right interface flux
            double f_h_R, f_hu_R, f_hv_R;
            riemannSolverHLLC(h[k], hu[k], hv[k],
                             h[kR], hu[kR], hv[kR],
                             f_h_R, f_hu_R, f_hv_R);
            
            // Update
            h_new[k] -= dt_local / dx * (f_h_R - f_h_L);
            hu_new[k] -= dt_local / dx * (f_hu_R - f_hu_L);
            hv_new[k] -= dt_local / dx * (f_hv_R - f_hv_L);
            
            // Bathymetry source term
            double b_L = depth[kL] > 0 ? -depth[kL] : depth[kL];
            double b_R = depth[kR] > 0 ? -depth[kR] : depth[kR];
            double h_bar = 0.5 * (h[kL] + h[kR]);
            hu_new[k] -= dt_local / dx * G * h_bar * (b_R - b_L) / (2.0 * dx);
        }
    }
    
    // Y-direction fluxes
    #pragma omp parallel for
    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {
            int k = idx(i, j);
            int kB = idx(i, j - 1);
            int kT = idx(i, j + 1);
            
            // Bottom interface flux
            double g_h_B, g_hu_B, g_hv_B;
            riemannSolverHLLC(h[kB], hv[kB], hu[kB],
                             h[k], hv[k], hu[k],
                             g_h_B, g_hv_B, g_hu_B);
            
            // Top interface flux
            double g_h_T, g_hu_T, g_hv_T;
            riemannSolverHLLC(h[k], hv[k], hu[k],
                             h[kT], hv[kT], hu[kT],
                             g_h_T, g_hv_T, g_hu_T);
            
            // Update
            h_new[k] -= dt_local / dy * (g_h_T - g_h_B);
            hu_new[k] -= dt_local / dy * (g_hu_T - g_hu_B);
            hv_new[k] -= dt_local / dy * (g_hv_T - g_hv_B);
            
            // Bathymetry source term
            double b_B = depth[kB] > 0 ? -depth[kB] : depth[kB];
            double b_T = depth[kT] > 0 ? -depth[kT] : depth[kT];
            double h_bar = 0.5 * (h[kB] + h[kT]);
            hv_new[k] -= dt_local / dy * G * h_bar * (b_T - b_B) / (2.0 * dy);
        }
    }
    
    // Update solution
    h = h_new;
    hu = hu_new;
    hv = hv_new;
    
    // Apply friction
    applyBottomFriction(dt_local);
    
    // Wetting/drying
    applyWettingDrying();
    
    // Boundary conditions
    applyBoundaryConditions();
    
    // Update eta
    for (int i = 0; i < n; ++i) {
        double d = depth[i];
        if (d > 0) {  // Ocean
            eta[i] = h[i] - d;
        } else {  // Land
            eta[i] = h[i] + d;  // Flow depth above land surface
        }
    }
    
    // Track maximum values
    for (int i = 0; i < n; ++i) {
        max_eta[i] = std::max(max_eta[i], eta[i]);
        
        if (h[i] > config.dry_tolerance) {
            double u_local = hu[i] / h[i];
            double v_local = hv[i] / h[i];
            double speed = std::sqrt(u_local * u_local + v_local * v_local);
            double mflux = 1000.0 * h[i] * speed * speed;  // ρhu²
            
            max_speed[i] = std::max(max_speed[i], speed);
            max_momentum_flux[i] = std::max(max_momentum_flux[i], mflux);
            
            if (arrival_time[i] < 0 && std::abs(eta[i]) > 0.01) {
                arrival_time[i] = current_time;
            }
        }
    }
    
    current_time += dt_local;
    return dt_local;
}

void ShallowWaterSolver::run(double end_time) {
    while (current_time < end_time) {
        double dt_taken = step();
        recordGauges(current_time);
        
        // Progress output
        if (static_cast<int>(current_time) % 300 == 0) {
            std::cout << "Time: " << current_time << " s" << std::endl;
        }
    }
}

void ShallowWaterSolver::recordGauges(double t) {
    for (auto& gauge : gauges) {
        double eta_val = getEta(gauge.longitude, gauge.latitude);
        double u_val, v_val;
        getVelocity(gauge.longitude, gauge.latitude, u_val, v_val);
        
        gauge.time.push_back(t);
        gauge.eta.push_back(eta_val);
        gauge.u.push_back(u_val);
        gauge.v.push_back(v_val);
    }
}

double ShallowWaterSolver::getEta(double lon, double lat) const {
    double fi = (lon - lon_min) / config.base_resolution;
    double fj = (lat - lat_min) / config.base_resolution;
    
    int i = std::max(0, std::min(static_cast<int>(fi), nx - 1));
    int j = std::max(0, std::min(static_cast<int>(fj), ny - 1));
    
    return eta[idx(i, j)];
}

void ShallowWaterSolver::getVelocity(double lon, double lat, 
                                      double& u_out, double& v_out) const {
    double fi = (lon - lon_min) / config.base_resolution;
    double fj = (lat - lat_min) / config.base_resolution;
    
    int i = std::max(0, std::min(static_cast<int>(fi), nx - 1));
    int j = std::max(0, std::min(static_cast<int>(fj), ny - 1));
    int k = idx(i, j);
    
    if (h[k] > config.dry_tolerance) {
        u_out = hu[k] / h[k];
        v_out = hv[k] / h[k];
    } else {
        u_out = v_out = 0.0;
    }
}

void ShallowWaterSolver::getSolution(std::vector<double>& eta_out,
                                      std::vector<double>& hu_out,
                                      std::vector<double>& hv_out) const {
    eta_out = eta;
    hu_out = hu;
    hv_out = hv;
}

void ShallowWaterSolver::getMaximumValues(std::vector<double>& max_eta_out,
                                           std::vector<double>& max_speed_out,
                                           std::vector<double>& max_mflux_out) const {
    max_eta_out = max_eta;
    max_speed_out = max_speed;
    max_mflux_out = max_momentum_flux;
}

void ShallowWaterSolver::addGauge(const TsunamiGauge& gauge) {
    gauges.push_back(gauge);
}

// =============================================================================
// CoupledTsunamiEarthquake Implementation
// =============================================================================

CoupledTsunamiEarthquake::CoupledTsunamiEarthquake() {}

void CoupledTsunamiEarthquake::initialize(
    const std::vector<TsunamiSubfault>& fault_model,
    const BathymetryGrid& bathy,
    const TsunamiConfig& cfg) {
    
    subfaults = fault_model;
    bathymetry = bathy;
    config = cfg;
    
    swe_solver = std::make_unique<ShallowWaterSolver>();
    swe_solver->initialize(bathymetry, config);
}

void CoupledTsunamiEarthquake::initializeCascadia(
    const CascadiaFaultModel& csz_model,
    const BathymetryGrid& bathy,
    const TsunamiConfig& cfg) {
    
    subfaults = csz_model.generateSubfaults();
    bathymetry = bathy;
    config = cfg;
    
    std::cout << "Cascadia model: " << subfaults.size() << " subfaults\n";
    std::cout << "Magnitude: Mw " << csz_model.getMagnitude() << "\n";
    std::cout << "Total moment: " << csz_model.getTotalMoment() << " N·m\n";
    
    swe_solver = std::make_unique<ShallowWaterSolver>();
    swe_solver->initialize(bathymetry, config);
}

void CoupledTsunamiEarthquake::computeSeafloorDisplacement() {
    int n = bathymetry.nlon * bathymetry.nlat;
    seafloor_displacement.resize(n);
    
    okada.computeDisplacementField(bathymetry, subfaults, seafloor_displacement);
    okada.applyKajiuraFilter(bathymetry, seafloor_displacement);
    
    // Statistics
    double max_up = *std::max_element(seafloor_displacement.begin(), 
                                       seafloor_displacement.end());
    double min_up = *std::min_element(seafloor_displacement.begin(), 
                                       seafloor_displacement.end());
    
    std::cout << "Seafloor displacement: max uplift = " << max_up 
              << " m, max subsidence = " << -min_up << " m\n";
}

void CoupledTsunamiEarthquake::generateKinematicSource() {
    // Create a time-dependent source function
    auto source_func = [this](double t, std::vector<double>& dz) {
        // Get final displacement
        if (seafloor_displacement.empty()) {
            dz.assign(bathymetry.nlon * bathymetry.nlat, 0.0);
            return;
        }
        
        // Ramp up over source rise time
        double factor = 0.0;
        if (t < config.source_rise_time) {
            factor = 0.5 * (1.0 - std::cos(M_PI * t / config.source_rise_time));
        } else {
            factor = 1.0;
        }
        
        // Previous factor for incremental displacement
        double prev_factor = 0.0;
        if (t > 0.5) {  // dt step
            double prev_t = t - 0.5;
            if (prev_t < config.source_rise_time) {
                prev_factor = 0.5 * (1.0 - std::cos(M_PI * prev_t / config.source_rise_time));
            } else {
                prev_factor = 1.0;
            }
        }
        
        double dfactor = factor - prev_factor;
        
        dz.resize(seafloor_displacement.size());
        for (size_t i = 0; i < seafloor_displacement.size(); ++i) {
            dz[i] = dfactor * seafloor_displacement[i];
        }
    };
    
    swe_solver->setSeafloorMotion(source_func);
}

void CoupledTsunamiEarthquake::run() {
    std::cout << "Computing seafloor displacement...\n";
    computeSeafloorDisplacement();
    
    std::cout << "Setting up kinematic source...\n";
    generateKinematicSource();
    
    // Set initial condition to zero
    std::vector<double> eta0(bathymetry.nlon * bathymetry.nlat, 0.0);
    swe_solver->setInitialCondition(eta0);
    
    std::cout << "Running tsunami propagation to " << config.end_time << " s...\n";
    swe_solver->run(config.end_time);
    
    std::cout << "Simulation complete.\n";
}

std::map<std::string, double> CoupledTsunamiEarthquake::getArrivalTimes() const {
    std::map<std::string, double> arrivals;
    for (const auto& gauge : swe_solver->getGauges()) {
        arrivals[gauge.name] = gauge.getFirstArrivalTime();
    }
    return arrivals;
}

std::map<std::string, double> CoupledTsunamiEarthquake::getMaxAmplitudes() const {
    std::map<std::string, double> amplitudes;
    for (const auto& gauge : swe_solver->getGauges()) {
        amplitudes[gauge.name] = gauge.getMaxAmplitude();
    }
    return amplitudes;
}

const std::vector<TsunamiGauge>& CoupledTsunamiEarthquake::getGaugeRecords() const {
    return swe_solver->getGauges();
}

void CoupledTsunamiEarthquake::writeOutput(const std::string& output_dir) const {
    // Write gauge records
    for (const auto& gauge : swe_solver->getGauges()) {
        std::string filename = output_dir + "/" + gauge.id + "_timeseries.dat";
        gauge.writeASCII(filename);
    }
    
    // Write maximum values
    std::vector<double> max_eta, max_speed, max_mflux;
    swe_solver->getMaximumValues(max_eta, max_speed, max_mflux);
    
    std::ofstream max_file(output_dir + "/maximum_values.dat");
    max_file << "# Maximum tsunami values\n";
    max_file << "# lon lat max_eta(m) max_speed(m/s) max_momentum_flux(kg/m/s2)\n";
    
    for (int j = 0; j < bathymetry.nlat; ++j) {
        for (int i = 0; i < bathymetry.nlon; ++i) {
            int k = i + j * bathymetry.nlon;
            double lon = bathymetry.lon_min + i * bathymetry.dlon;
            double lat = bathymetry.lat_min + j * bathymetry.dlat;
            
            max_file << std::fixed << std::setprecision(4)
                     << lon << " " << lat << " "
                     << max_eta[k] << " " << max_speed[k] << " "
                     << max_mflux[k] << "\n";
        }
    }
}

// =============================================================================
// WestCoastGaugeNetwork Implementation
// =============================================================================

TsunamiGauge WestCoastGaugeNetwork::crescent_city() {
    TsunamiGauge g;
    g.name = "Crescent City";
    g.id = "9419750";
    g.longitude = -124.183;
    g.latitude = 41.745;
    g.depth = 0;
    g.is_dart_buoy = false;
    return g;
}

TsunamiGauge WestCoastGaugeNetwork::astoria() {
    TsunamiGauge g;
    g.name = "Astoria";
    g.id = "9439040";
    g.longitude = -123.768;
    g.latitude = 46.208;
    g.depth = 0;
    g.is_dart_buoy = false;
    return g;
}

TsunamiGauge WestCoastGaugeNetwork::westport() {
    TsunamiGauge g;
    g.name = "Westport";
    g.id = "9441102";
    g.longitude = -124.105;
    g.latitude = 46.904;
    g.depth = 0;
    g.is_dart_buoy = false;
    return g;
}

TsunamiGauge WestCoastGaugeNetwork::seattle() {
    TsunamiGauge g;
    g.name = "Seattle";
    g.id = "9447130";
    g.longitude = -122.339;
    g.latitude = 47.602;
    g.depth = 0;
    g.is_dart_buoy = false;
    return g;
}

TsunamiGauge WestCoastGaugeNetwork::victoria() {
    TsunamiGauge g;
    g.name = "Victoria";
    g.id = "7120";
    g.longitude = -123.371;
    g.latitude = 48.425;
    g.depth = 0;
    g.is_dart_buoy = false;
    return g;
}

TsunamiGauge WestCoastGaugeNetwork::tofino() {
    TsunamiGauge g;
    g.name = "Tofino";
    g.id = "8615";
    g.longitude = -125.913;
    g.latitude = 49.154;
    g.depth = 0;
    g.is_dart_buoy = false;
    return g;
}

TsunamiGauge WestCoastGaugeNetwork::san_francisco() {
    TsunamiGauge g;
    g.name = "San Francisco";
    g.id = "9414290";
    g.longitude = -122.466;
    g.latitude = 37.807;
    g.depth = 0;
    g.is_dart_buoy = false;
    return g;
}

TsunamiGauge WestCoastGaugeNetwork::monterey() {
    TsunamiGauge g;
    g.name = "Monterey";
    g.id = "9413450";
    g.longitude = -121.889;
    g.latitude = 36.605;
    g.depth = 0;
    g.is_dart_buoy = false;
    return g;
}

TsunamiGauge WestCoastGaugeNetwork::los_angeles() {
    TsunamiGauge g;
    g.name = "Los Angeles";
    g.id = "9410660";
    g.longitude = -118.272;
    g.latitude = 33.720;
    g.depth = 0;
    g.is_dart_buoy = false;
    return g;
}

TsunamiGauge WestCoastGaugeNetwork::san_diego() {
    TsunamiGauge g;
    g.name = "San Diego";
    g.id = "9410170";
    g.longitude = -117.174;
    g.latitude = 32.714;
    g.depth = 0;
    g.is_dart_buoy = false;
    return g;
}

TsunamiGauge WestCoastGaugeNetwork::dart_46404() {
    TsunamiGauge g;
    g.name = "DART 46404";
    g.id = "46404";
    g.longitude = -128.894;
    g.latitude = 45.857;
    g.depth = 2785;
    g.is_dart_buoy = true;
    return g;
}

TsunamiGauge WestCoastGaugeNetwork::dart_46407() {
    TsunamiGauge g;
    g.name = "DART 46407";
    g.id = "46407";
    g.longitude = -128.742;
    g.latitude = 42.641;
    g.depth = 3284;
    g.is_dart_buoy = true;
    return g;
}

TsunamiGauge WestCoastGaugeNetwork::dart_46419() {
    TsunamiGauge g;
    g.name = "DART 46419";
    g.id = "46419";
    g.longitude = -129.584;
    g.latitude = 48.764;
    g.depth = 2762;
    g.is_dart_buoy = true;
    return g;
}

std::vector<TsunamiGauge> WestCoastGaugeNetwork::getTideGauges() {
    return {
        crescent_city(),
        astoria(),
        westport(),
        seattle(),
        victoria(),
        tofino(),
        san_francisco(),
        monterey(),
        los_angeles(),
        san_diego()
    };
}

std::vector<TsunamiGauge> WestCoastGaugeNetwork::getDARTBuoys() {
    return {
        dart_46404(),
        dart_46407(),
        dart_46419()
    };
}

std::vector<TsunamiGauge> WestCoastGaugeNetwork::getAllStations() {
    auto gauges = getTideGauges();
    auto darts = getDARTBuoys();
    gauges.insert(gauges.end(), darts.begin(), darts.end());
    return gauges;
}

// =============================================================================
// CascadiaScenarios Implementation
// =============================================================================

namespace CascadiaScenarios {

CascadiaFaultModel fullMarginM9() {
    CascadiaFaultModel model;
    model.north_latitude = 50.0;      // Vancouver Island
    model.south_latitude = 40.5;      // Cape Mendocino
    model.trench_longitude = -125.0;
    model.average_strike = 350.0;
    model.shallow_dip = 8.0;
    model.deep_dip = 18.0;
    model.transition_depth = 20.0;
    model.rupture_width = 120.0;
    model.average_slip = 17.0;
    model.peak_slip = 35.0;
    model.num_along_strike = 60;
    model.num_down_dip = 12;
    model.rupture_velocity = 2.8;
    model.hypocenter_lat = 45.0;
    model.hypocenter_lon = -124.5;
    model.hypocenter_depth = 25.0;
    return model;
}

CascadiaFaultModel southernSegmentM8() {
    CascadiaFaultModel model;
    model.north_latitude = 44.0;
    model.south_latitude = 40.5;
    model.trench_longitude = -124.5;
    model.average_strike = 350.0;
    model.shallow_dip = 10.0;
    model.deep_dip = 20.0;
    model.transition_depth = 18.0;
    model.rupture_width = 100.0;
    model.average_slip = 8.0;
    model.peak_slip = 15.0;
    model.num_along_strike = 25;
    model.num_down_dip = 10;
    model.rupture_velocity = 2.5;
    model.hypocenter_lat = 42.0;
    model.hypocenter_lon = -124.3;
    model.hypocenter_depth = 22.0;
    return model;
}

CascadiaFaultModel northernSegmentM8() {
    CascadiaFaultModel model;
    model.north_latitude = 50.0;
    model.south_latitude = 46.0;
    model.trench_longitude = -125.5;
    model.average_strike = 350.0;
    model.shallow_dip = 7.0;
    model.deep_dip = 15.0;
    model.transition_depth = 22.0;
    model.rupture_width = 110.0;
    model.average_slip = 10.0;
    model.peak_slip = 20.0;
    model.num_along_strike = 30;
    model.num_down_dip = 11;
    model.rupture_velocity = 2.6;
    model.hypocenter_lat = 48.0;
    model.hypocenter_lon = -125.0;
    model.hypocenter_depth = 20.0;
    return model;
}

CascadiaFaultModel centralOregonM8() {
    CascadiaFaultModel model;
    model.north_latitude = 46.0;
    model.south_latitude = 43.0;
    model.trench_longitude = -125.0;
    model.average_strike = 355.0;
    model.shallow_dip = 9.0;
    model.deep_dip = 18.0;
    model.transition_depth = 20.0;
    model.rupture_width = 100.0;
    model.average_slip = 12.0;
    model.peak_slip = 25.0;
    model.num_along_strike = 20;
    model.num_down_dip = 10;
    model.rupture_velocity = 2.7;
    model.hypocenter_lat = 44.5;
    model.hypocenter_lon = -124.7;
    model.hypocenter_depth = 23.0;
    return model;
}

CascadiaFaultModel worstCaseScenario() {
    CascadiaFaultModel model = fullMarginM9();
    model.average_slip = 22.0;
    model.peak_slip = 45.0;
    model.num_along_strike = 80;
    model.num_down_dip = 15;
    return model;
}

} // namespace CascadiaScenarios

// =============================================================================
// TsunamiUtils Implementation
// =============================================================================

namespace TsunamiUtils {

double computeMomentMagnitude(const std::vector<TsunamiSubfault>& subfaults,
                              double shear_modulus) {
    double moment = 0.0;
    for (const auto& sf : subfaults) {
        double area = sf.length * 1000.0 * sf.width * 1000.0;  // m²
        moment += shear_modulus * area * sf.slip;
    }
    return (2.0/3.0) * (std::log10(moment) - 9.1);
}

double waveSpeed(double depth) {
    return std::sqrt(G * depth);
}

double estimateArrivalTime(double lon_source, double lat_source,
                           double lon_target, double lat_target,
                           double average_depth) {
    double dist = haversineDistance(lon_source, lat_source, lon_target, lat_target);
    double speed = waveSpeed(average_depth);
    return dist / speed;
}

double greensLaw(double eta_deep, double depth_deep, double depth_shallow) {
    return eta_deep * std::pow(depth_deep / depth_shallow, 0.25);
}

double estimateInundation(double wave_height, double beach_slope) {
    // Synolakis (1987) runup formula for breaking waves
    double tan_beta = beach_slope;
    return 2.831 * std::pow(tan_beta, -0.5) * wave_height;
}

double haversineDistance(double lon1, double lat1, double lon2, double lat2) {
    double dlon = (lon2 - lon1) * DEG_TO_RAD;
    double dlat = (lat2 - lat1) * DEG_TO_RAD;
    
    double a = std::sin(dlat/2) * std::sin(dlat/2) +
               std::cos(lat1 * DEG_TO_RAD) * std::cos(lat2 * DEG_TO_RAD) *
               std::sin(dlon/2) * std::sin(dlon/2);
    double c = 2 * std::atan2(std::sqrt(a), std::sqrt(1-a));
    
    return EARTH_RADIUS * c;
}

BathymetryGrid generateTestBathymetry(double lon_min, double lon_max,
                                       double lat_min, double lat_max,
                                       double resolution,
                                       double shelf_depth,
                                       double deep_depth) {
    BathymetryGrid grid;
    
    grid.lon_min = lon_min;
    grid.lon_max = lon_max;
    grid.lat_min = lat_min;
    grid.lat_max = lat_max;
    grid.dlon = resolution;
    grid.dlat = resolution;
    
    grid.nlon = static_cast<int>((lon_max - lon_min) / resolution);
    grid.nlat = static_cast<int>((lat_max - lat_min) / resolution);
    
    grid.depth.resize(grid.nlon * grid.nlat);
    
    // Generate continental shelf profile
    double coast_lon = lon_max - 0.5;  // Coast near eastern edge
    double shelf_width = 50.0;  // km
    double slope_width = 30.0;  // km
    
    for (int j = 0; j < grid.nlat; ++j) {
        for (int i = 0; i < grid.nlon; ++i) {
            double lon = lon_min + i * resolution;
            
            // Distance from coast in km (approximate)
            double lat_avg = 0.5 * (lat_min + lat_max);
            double dist = (coast_lon - lon) * 111.0 * std::cos(lat_avg * DEG_TO_RAD);
            
            double d;
            if (dist < 0) {
                // Land
                d = -10.0 + dist * 0.1;  // Gentle slope inland
            } else if (dist < shelf_width) {
                // Continental shelf
                d = shelf_depth * dist / shelf_width;
            } else if (dist < shelf_width + slope_width) {
                // Continental slope
                double frac = (dist - shelf_width) / slope_width;
                d = shelf_depth + frac * (deep_depth - shelf_depth);
            } else {
                // Deep ocean
                d = deep_depth;
            }
            
            grid.depth[i + j * grid.nlon] = d;
        }
    }
    
    return grid;
}

} // namespace TsunamiUtils

} // namespace FSRM
