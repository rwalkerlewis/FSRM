/**
 * @file OceanPhysics.cpp
 * @brief Implementation of ocean physics modeling components
 */

#include "OceanPhysics.hpp"
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
static const double CP_WATER = 4000.0;           // Specific heat of water (J/kg/K)
static const double STEFAN_BOLTZMANN = 5.67e-8;  // Stefan-Boltzmann constant
static const double L_EVAP = 2.5e6;              // Latent heat of evaporation (J/kg)

// =============================================================================
// WaveSpectrum2D Implementation
// =============================================================================

void WaveSpectrum2D::computeIntegralParameters() {
    double m0 = 0.0, m1 = 0.0, m2 = 0.0;
    double Ex = 0.0, Ey = 0.0;
    double E_peak = 0.0;
    int peak_f = 0, peak_d = 0;
    
    int nf = frequencies.size();
    int nd = directions.size();
    
    if (nf == 0 || nd == 0) return;
    
    double df = (nf > 1) ? frequencies[1] - frequencies[0] : 0.01;
    double dd = (nd > 1) ? directions[1] - directions[0] : M_PI / 18.0;
    
    for (int f = 0; f < nf; ++f) {
        for (int d = 0; d < nd; ++d) {
            double E = energy[f][d];
            double freq = frequencies[f];
            double dir = directions[d];
            
            m0 += E * df * dd;
            m1 += freq * E * df * dd;
            m2 += freq * freq * E * df * dd;
            
            Ex += E * std::cos(dir) * df * dd;
            Ey += E * std::sin(dir) * df * dd;
            
            if (E > E_peak) {
                E_peak = E;
                peak_f = f;
                peak_d = d;
            }
        }
    }
    
    Hs = 4.0 * std::sqrt(m0);
    Tm01 = (m1 > 0) ? m0 / m1 : 0.0;
    Tm02 = (m2 > 0) ? std::sqrt(m0 / m2) : 0.0;
    Tp = (peak_f < nf) ? 1.0 / frequencies[peak_f] : 0.0;
    
    theta_mean = std::atan2(Ey, Ex);
    theta_peak = (peak_d < nd) ? directions[peak_d] : 0.0;
    
    // Directional spread
    double sum_cos2 = 0.0;
    for (int f = 0; f < nf; ++f) {
        for (int d = 0; d < nd; ++d) {
            double dtheta = directions[d] - theta_mean;
            sum_cos2 += energy[f][d] * std::cos(2.0 * dtheta) * df * dd;
        }
    }
    spread = (m0 > 0) ? std::sqrt(2.0 * (1.0 - sum_cos2 / m0)) : 0.0;
}

double WaveSpectrum2D::getEnergy(double f, double theta) const {
    // Bilinear interpolation
    int nf = frequencies.size();
    int nd = directions.size();
    
    if (nf == 0 || nd == 0) return 0.0;
    
    // Find frequency index
    int fi = 0;
    for (int i = 0; i < nf - 1; ++i) {
        if (f >= frequencies[i] && f < frequencies[i + 1]) {
            fi = i;
            break;
        }
    }
    
    // Find direction index
    int di = 0;
    for (int i = 0; i < nd - 1; ++i) {
        if (theta >= directions[i] && theta < directions[i + 1]) {
            di = i;
            break;
        }
    }
    
    double df = frequencies[fi + 1] - frequencies[fi];
    double dd = directions[di + 1] - directions[di];
    
    double u = (f - frequencies[fi]) / df;
    double v = (theta - directions[di]) / dd;
    
    return (1 - u) * (1 - v) * energy[fi][di] +
           u * (1 - v) * energy[fi + 1][di] +
           (1 - u) * v * energy[fi][di + 1] +
           u * v * energy[fi + 1][di + 1];
}

double WaveSpectrum2D::totalEnergy() const {
    double E_total = 0.0;
    int nf = frequencies.size();
    int nd = directions.size();
    
    if (nf == 0 || nd == 0) return 0.0;
    
    double df = (nf > 1) ? frequencies[1] - frequencies[0] : 0.01;
    double dd = (nd > 1) ? directions[1] - directions[0] : M_PI / 18.0;
    
    for (int f = 0; f < nf; ++f) {
        for (int d = 0; d < nd; ++d) {
            E_total += energy[f][d] * df * dd;
        }
    }
    
    return E_total;
}

void WaveSpectrum2D::normalize(double target_Hs) {
    double current_Hs = 4.0 * std::sqrt(totalEnergy());
    if (current_Hs > 0) {
        double scale = (target_Hs / current_Hs) * (target_Hs / current_Hs);
        for (auto& row : energy) {
            for (auto& val : row) {
                val *= scale;
            }
        }
    }
}

WaveSpectrum2D WaveSpectrum2D::generateJONSWAP(double Hs, double Tp, double gamma,
                                               double theta_mean, double spread_deg,
                                               int nf, int nd) {
    WaveSpectrum2D spec;
    
    spec.frequencies.resize(nf);
    spec.directions.resize(nd);
    spec.energy.resize(nf, std::vector<double>(nd));
    
    double fp = 1.0 / Tp;
    double f_min = 0.5 * fp;
    double f_max = 5.0 * fp;
    
    for (int i = 0; i < nf; ++i) {
        spec.frequencies[i] = f_min + i * (f_max - f_min) / (nf - 1);
    }
    
    for (int i = 0; i < nd; ++i) {
        spec.directions[i] = -M_PI + i * 2.0 * M_PI / nd;
    }
    
    // JONSWAP parameters
    double alpha = 0.0624 / (0.23 + 0.0336 * gamma - 0.185 / (1.9 + gamma));
    alpha *= std::pow(Hs * fp * fp / G, 2);
    
    double spread_rad = spread_deg * DEG_TO_RAD;
    
    for (int f = 0; f < nf; ++f) {
        double freq = spec.frequencies[f];
        double sigma = (freq < fp) ? 0.07 : 0.09;
        
        // PM part
        double S_pm = alpha * G * G * std::pow(2.0 * M_PI, -4) * std::pow(freq, -5) *
                     std::exp(-1.25 * std::pow(fp / freq, 4));
        
        // Peak enhancement
        double gamma_factor = std::pow(gamma, 
            std::exp(-0.5 * std::pow((freq - fp) / (sigma * fp), 2)));
        
        double S_f = S_pm * gamma_factor;
        
        // Directional spreading (cosine-squared)
        for (int d = 0; d < nd; ++d) {
            double dtheta = spec.directions[d] - theta_mean * DEG_TO_RAD;
            while (dtheta > M_PI) dtheta -= 2.0 * M_PI;
            while (dtheta < -M_PI) dtheta += 2.0 * M_PI;
            
            double D = 0.0;
            if (std::abs(dtheta) < M_PI / 2.0) {
                double s = 2.0 / (spread_rad * spread_rad);  // Spreading parameter
                D = std::pow(std::cos(dtheta / 2.0), 2.0 * s);
            }
            
            spec.energy[f][d] = S_f * D;
        }
    }
    
    // Normalize to target Hs
    spec.normalize(Hs);
    spec.computeIntegralParameters();
    
    return spec;
}

WaveSpectrum2D WaveSpectrum2D::generatePM(double U10, double theta_mean,
                                          double spread_deg, int nf, int nd) {
    // Pierson-Moskowitz spectrum
    double alpha_pm = 0.0081;
    double beta_pm = 0.74;
    double fp = 0.877 * G / (2.0 * M_PI * U10);
    double Hs = 0.22 * U10 * U10 / G;
    double Tp = 1.0 / fp;
    
    // Use JONSWAP with gamma=1 (which gives PM)
    return generateJONSWAP(Hs, Tp, 1.0, theta_mean, spread_deg, nf, nd);
}

// =============================================================================
// WaveModel Implementation
// =============================================================================

WaveModel::WaveModel() : current_time(0.0) {}

WaveModel::~WaveModel() = default;

void WaveModel::initialize(const WaveModelConfig& cfg) {
    config = cfg;
    
    nx = cfg.nx;
    ny = cfg.ny;
    nf = cfg.num_frequencies;
    nd = cfg.num_directions;
    
    dx = (cfg.lon_max - cfg.lon_min) / (nx - 1) * EARTH_RADIUS * DEG_TO_RAD *
         std::cos(0.5 * (cfg.lat_min + cfg.lat_max) * DEG_TO_RAD);
    dy = (cfg.lat_max - cfg.lat_min) / (ny - 1) * EARTH_RADIUS * DEG_TO_RAD;
    
    // Set up frequency grid (logarithmic)
    freq.resize(nf);
    df.resize(nf);
    double f_ratio = std::pow(cfg.freq_max / cfg.freq_min, 1.0 / (nf - 1));
    for (int i = 0; i < nf; ++i) {
        freq[i] = cfg.freq_min * std::pow(f_ratio, i);
        if (i > 0) {
            df[i] = freq[i] - freq[i - 1];
        } else {
            df[i] = freq[1] - freq[0];
        }
    }
    
    // Direction grid (equidistant)
    dir.resize(nd);
    dd.resize(nd);
    for (int i = 0; i < nd; ++i) {
        dir[i] = 2.0 * M_PI * i / nd;
        dd[i] = 2.0 * M_PI / nd;
    }
    
    // Allocate arrays
    int n2d = nx * ny;
    int n4d = n2d * nf * nd;
    
    action.resize(n4d, 0.0);
    
    depth.resize(n2d, 100.0);
    mask.resize(n2d, 1);
    
    Cg_x.resize(n4d);
    Cg_y.resize(n4d);
    C_theta.resize(n4d, 0.0);
    C_sigma.resize(n4d, 0.0);
    
    u10_field.resize(n2d, 0.0);
    v10_field.resize(n2d, 0.0);
    u_current.resize(n2d, 0.0);
    v_current.resize(n2d, 0.0);
    
    Hs_field.resize(n2d, 0.0);
    Tp_field.resize(n2d, 0.0);
    dir_field.resize(n2d, 0.0);
    Sxx_field.resize(n2d, 0.0);
    Sxy_field.resize(n2d, 0.0);
    Syy_field.resize(n2d, 0.0);
    tau_wave.resize(n2d, 0.0);
    
    current_time = 0.0;
    
    std::cout << "WaveModel initialized:\n"
              << "  Grid: " << nx << " x " << ny << "\n"
              << "  Frequencies: " << nf << " (" << cfg.freq_min << " - " << cfg.freq_max << " Hz)\n"
              << "  Directions: " << nd << "\n";
}

void WaveModel::setBathymetry(const std::vector<double>& d) {
    depth = d;
    for (size_t i = 0; i < depth.size(); ++i) {
        mask[i] = (depth[i] > 0) ? 1 : 0;
    }
    computeGroupVelocity();
}

void WaveModel::setBathymetry(const BathymetryGrid& bathymetry) {
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            double lon = config.lon_min + i * (config.lon_max - config.lon_min) / (nx - 1);
            double lat = config.lat_min + j * (config.lat_max - config.lat_min) / (ny - 1);
            
            int idx2 = idx2d(i, j);
            depth[idx2] = std::max(0.1, bathymetry.getDepth(lon, lat));
            mask[idx2] = (depth[idx2] > 0) ? 1 : 0;
        }
    }
    computeGroupVelocity();
}

void WaveModel::setWind(const std::vector<double>& u10,
                       const std::vector<double>& v10) {
    u10_field = u10;
    v10_field = v10;
}

void WaveModel::setCurrents(const std::vector<double>& u,
                           const std::vector<double>& v) {
    u_current = u;
    v_current = v;
}

void WaveModel::setInitialSpectrum(const WaveSpectrum2D& spectrum) {
    // Initialize all grid points with the same spectrum
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            if (mask[idx2d(i, j)] == 0) continue;
            
            for (int f = 0; f < nf; ++f) {
                double omega = 2.0 * M_PI * freq[f];
                
                for (int d = 0; d < nd; ++d) {
                    double E = spectrum.getEnergy(freq[f], dir[d]);
                    // Convert energy to action: N = E / ω
                    action[idx(i, j, f, d)] = E / omega;
                }
            }
        }
    }
    
    computeIntegratedParameters();
}

void WaveModel::setInitialParameters(double Hs, double Tp, double dir_deg) {
    auto spectrum = WaveSpectrum2D::generateJONSWAP(Hs, Tp, 3.3, dir_deg, 30.0, nf, nd);
    setInitialSpectrum(spectrum);
}

double WaveModel::dispersion(double f, double h) const {
    // Solve ω² = gk tanh(kh)
    double omega = 2.0 * M_PI * f;
    double omega2 = omega * omega;
    
    // Initial guess (deep water)
    double k = omega2 / G;
    
    // Newton-Raphson iteration
    for (int iter = 0; iter < 20; ++iter) {
        double tanh_kh = std::tanh(k * h);
        double F = omega2 - G * k * tanh_kh;
        double dF = -G * (tanh_kh + k * h * (1.0 - tanh_kh * tanh_kh));
        
        double dk = -F / dF;
        k += dk;
        
        if (std::abs(dk) < 1e-10 * k) break;
    }
    
    return k;
}

double WaveModel::groupVelocity(double f, double h) const {
    double k = dispersion(f, h);
    double omega = 2.0 * M_PI * f;
    double kh = k * h;
    
    // Cg = ω/k * (1/2) * (1 + 2kh/sinh(2kh))
    double n = 0.5 * (1.0 + 2.0 * kh / std::sinh(2.0 * kh));
    return omega / k * n;
}

void WaveModel::computeGroupVelocity() {
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int idx2 = idx2d(i, j);
            double h = depth[idx2];
            
            for (int f = 0; f < nf; ++f) {
                double Cg = groupVelocity(freq[f], h);
                
                for (int d = 0; d < nd; ++d) {
                    int idx4 = idx(i, j, f, d);
                    Cg_x[idx4] = Cg * std::cos(dir[d]);
                    Cg_y[idx4] = Cg * std::sin(dir[d]);
                }
            }
        }
    }
}

void WaveModel::computeRefraction() {
    if (!config.enable_refraction) return;
    
    // Compute direction rate of change due to depth gradient
    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {
            int idx2 = idx2d(i, j);
            
            // Depth gradient
            double dh_dx = (depth[idx2d(i + 1, j)] - depth[idx2d(i - 1, j)]) / (2.0 * dx);
            double dh_dy = (depth[idx2d(i, j + 1)] - depth[idx2d(i, j - 1)]) / (2.0 * dy);
            
            double h = depth[idx2];
            
            for (int f = 0; f < nf; ++f) {
                double k = dispersion(freq[f], h);
                double kh = k * h;
                double n = 0.5 * (1.0 + 2.0 * kh / std::sinh(2.0 * kh));
                double Cg = groupVelocity(freq[f], h);
                double c = 2.0 * M_PI * freq[f] / k;
                
                // dC/dh
                double dC_dh = c / (2.0 * h * n);
                
                for (int d = 0; d < nd; ++d) {
                    int idx4 = idx(i, j, f, d);
                    
                    double cos_d = std::cos(dir[d]);
                    double sin_d = std::sin(dir[d]);
                    
                    // C_theta = (c/c_g) * (sin(θ) * dc/dx - cos(θ) * dc/dy)
                    // where dc/dx = (dC/dh) * (dh/dx)
                    double dc_dx = dC_dh * dh_dx;
                    double dc_dy = dC_dh * dh_dy;
                    
                    C_theta[idx4] = (c / Cg) * (sin_d * dc_dx - cos_d * dc_dy);
                }
            }
        }
    }
}

void WaveModel::computeWindInput(std::vector<double>& S_in) {
    // Janssen (1991) wind input source term
    const double kappa = 0.41;  // von Karman constant
    const double beta_max = 1.2;
    
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int idx2 = idx2d(i, j);
            if (mask[idx2] == 0) continue;
            
            double u10 = u10_field[idx2];
            double v10 = v10_field[idx2];
            double wind_speed = std::sqrt(u10 * u10 + v10 * v10);
            double wind_dir = std::atan2(v10, u10);
            
            if (wind_speed < 0.1) continue;
            
            double h = depth[idx2];
            
            // Friction velocity
            double Cd = 0.001 * (0.8 + 0.065 * wind_speed);
            double u_star = std::sqrt(Cd) * wind_speed;
            
            for (int f = 0; f < nf; ++f) {
                double omega = 2.0 * M_PI * freq[f];
                double k = dispersion(freq[f], h);
                double c = omega / k;
                
                for (int d = 0; d < nd; ++d) {
                    int idx4 = idx(i, j, f, d);
                    
                    // Wind-wave angle
                    double cos_angle = std::cos(dir[d] - wind_dir);
                    
                    if (cos_angle > 0) {
                        // Miles parameter
                        double W = u_star * cos_angle / c;
                        double beta = beta_max * std::max(0.0, W - 1.0) * (u_star / c);
                        
                        S_in[idx4] = beta * omega * action[idx4];
                    } else {
                        S_in[idx4] = 0.0;
                    }
                }
            }
        }
    }
}

void WaveModel::computeWhitecapping(std::vector<double>& S_ds) {
    // Komen et al. (1984) whitecapping dissipation
    const double Cds = config.whitecap_coeff;
    const double delta = config.whitecap_delta;
    
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int idx2 = idx2d(i, j);
            if (mask[idx2] == 0) continue;
            
            double h = depth[idx2];
            
            // Mean steepness
            double m0 = 0.0, m1 = 0.0;
            for (int f = 0; f < nf; ++f) {
                double k = dispersion(freq[f], h);
                for (int d = 0; d < nd; ++d) {
                    double E = action[idx(i, j, f, d)] * 2.0 * M_PI * freq[f];
                    m0 += E * df[f] * dd[d];
                    m1 += k * E * df[f] * dd[d];
                }
            }
            
            double k_mean = (m0 > 1e-10) ? m1 / m0 : 0.01;
            double steepness = k_mean * std::sqrt(m0);
            double s_PM = 0.00302;  // PM steepness
            
            for (int f = 0; f < nf; ++f) {
                double omega = 2.0 * M_PI * freq[f];
                double k = dispersion(freq[f], h);
                
                for (int d = 0; d < nd; ++d) {
                    int idx4 = idx(i, j, f, d);
                    
                    double gamma = Cds * std::pow(steepness / s_PM, 2) *
                                  std::pow(k / k_mean, delta);
                    
                    S_ds[idx4] = -gamma * omega * action[idx4];
                }
            }
        }
    }
}

void WaveModel::computeDepthBreaking(std::vector<double>& S_brk) {
    // Battjes-Janssen (1978) depth-induced breaking
    const double gamma_br = config.breaking_gamma;
    const double alpha = 1.0;
    
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int idx2 = idx2d(i, j);
            if (mask[idx2] == 0) continue;
            
            double h = depth[idx2];
            double Hs = Hs_field[idx2];
            double Hmax = gamma_br * h;
            
            // Breaking probability
            double Qb = 0.0;
            if (Hs > 0.5 * Hmax) {
                // Solve Qb from: 1 - Qb = (Hs/Hmax)² ln(Qb)
                double ratio = Hs / Hmax;
                Qb = std::exp(-std::pow(Hmax / Hs, 2));
            }
            
            // Total energy dissipation rate
            double D_br = 0.25 * alpha * Qb * RHO_WATER * G * Hmax * Hmax / 
                         (Tp_field[idx2] > 0 ? Tp_field[idx2] : 10.0);
            
            // Distribute over spectrum
            double E_total = 0.0;
            for (int f = 0; f < nf; ++f) {
                for (int d = 0; d < nd; ++d) {
                    E_total += action[idx(i, j, f, d)] * 2.0 * M_PI * freq[f] * df[f] * dd[d];
                }
            }
            
            if (E_total > 1e-10) {
                for (int f = 0; f < nf; ++f) {
                    double omega = 2.0 * M_PI * freq[f];
                    for (int d = 0; d < nd; ++d) {
                        int idx4 = idx(i, j, f, d);
                        double E = action[idx4] * omega;
                        S_brk[idx4] = -D_br * (E / E_total) / omega;
                    }
                }
            }
        }
    }
}

void WaveModel::computeBottomFriction(std::vector<double>& S_bf) {
    // JONSWAP bottom friction
    const double C_bf = config.bottom_friction_coeff;
    
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int idx2 = idx2d(i, j);
            if (mask[idx2] == 0) continue;
            
            double h = depth[idx2];
            if (h > 100.0) continue;  // Deep water - no bottom friction
            
            for (int f = 0; f < nf; ++f) {
                double omega = 2.0 * M_PI * freq[f];
                double k = dispersion(freq[f], h);
                double kh = k * h;
                
                // Bottom friction coefficient
                double Cf = C_bf / std::sinh(kh) / std::sinh(kh);
                
                for (int d = 0; d < nd; ++d) {
                    int idx4 = idx(i, j, f, d);
                    S_bf[idx4] = -Cf * omega * action[idx4];
                }
            }
        }
    }
}

void WaveModel::propagate(double dt) {
    // First-order upwind advection
    std::vector<double> action_new = action;
    
    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {
            int idx2 = idx2d(i, j);
            if (mask[idx2] == 0) continue;
            
            for (int f = 0; f < nf; ++f) {
                for (int d = 0; d < nd; ++d) {
                    int idx4 = idx(i, j, f, d);
                    
                    // Spatial advection
                    double Cgx = Cg_x[idx4];
                    double Cgy = Cg_y[idx4];
                    
                    double dN_dx, dN_dy;
                    if (Cgx > 0) {
                        dN_dx = (action[idx4] - action[idx(i - 1, j, f, d)]) / dx;
                    } else {
                        dN_dx = (action[idx(i + 1, j, f, d)] - action[idx4]) / dx;
                    }
                    
                    if (Cgy > 0) {
                        dN_dy = (action[idx4] - action[idx(i, j - 1, f, d)]) / dy;
                    } else {
                        dN_dy = (action[idx(i, j + 1, f, d)] - action[idx4]) / dy;
                    }
                    
                    action_new[idx4] -= dt * (Cgx * dN_dx + Cgy * dN_dy);
                    
                    // Directional advection (refraction)
                    if (config.enable_refraction) {
                        double Ctheta = C_theta[idx4];
                        int d_up = (d + nd - 1) % nd;
                        int d_down = (d + 1) % nd;
                        
                        double dN_dtheta;
                        if (Ctheta > 0) {
                            dN_dtheta = (action[idx4] - action[idx(i, j, f, d_up)]) / dd[d];
                        } else {
                            dN_dtheta = (action[idx(i, j, f, d_down)] - action[idx4]) / dd[d];
                        }
                        
                        action_new[idx4] -= dt * Ctheta * dN_dtheta;
                    }
                }
            }
        }
    }
    
    action = action_new;
}

double WaveModel::step() {
    double dt = config.dt;
    
    // Compute source terms
    std::vector<double> S_in(action.size(), 0.0);
    std::vector<double> S_ds(action.size(), 0.0);
    std::vector<double> S_brk(action.size(), 0.0);
    std::vector<double> S_bf(action.size(), 0.0);
    
    if (config.enable_wind_input) computeWindInput(S_in);
    if (config.enable_whitecapping) computeWhitecapping(S_ds);
    if (config.enable_depth_breaking) computeDepthBreaking(S_brk);
    if (config.enable_bottom_friction) computeBottomFriction(S_bf);
    
    // Apply source terms
    for (size_t i = 0; i < action.size(); ++i) {
        action[i] += dt * (S_in[i] + S_ds[i] + S_brk[i] + S_bf[i]);
        action[i] = std::max(0.0, action[i]);
    }
    
    // Propagation
    computeRefraction();
    propagate(dt);
    
    // Update integrated parameters
    computeIntegratedParameters();
    
    current_time += dt;
    return dt;
}

void WaveModel::run(double end_time) {
    while (current_time < end_time) {
        step();
    }
}

void WaveModel::computeIntegratedParameters() {
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int idx2 = idx2d(i, j);
            if (mask[idx2] == 0) {
                Hs_field[idx2] = 0.0;
                Tp_field[idx2] = 0.0;
                dir_field[idx2] = 0.0;
                continue;
            }
            
            double m0 = 0.0, m1 = 0.0;
            double Ex = 0.0, Ey = 0.0;
            double E_peak = 0.0;
            int peak_f = 0;
            double Sxx = 0.0, Sxy = 0.0, Syy = 0.0;
            
            double h = depth[idx2];
            
            for (int f = 0; f < nf; ++f) {
                double omega = 2.0 * M_PI * freq[f];
                double k = dispersion(freq[f], h);
                double kh = k * h;
                double n = 0.5 * (1.0 + 2.0 * kh / std::sinh(2.0 * kh));
                double Cg = omega / k * n;
                double c = omega / k;
                
                for (int d = 0; d < nd; ++d) {
                    double E = action[idx(i, j, f, d)] * omega;
                    
                    m0 += E * df[f] * dd[d];
                    m1 += freq[f] * E * df[f] * dd[d];
                    
                    Ex += E * std::cos(dir[d]) * df[f] * dd[d];
                    Ey += E * std::sin(dir[d]) * df[f] * dd[d];
                    
                    if (E > E_peak) {
                        E_peak = E;
                        peak_f = f;
                    }
                    
                    // Radiation stress
                    double cos_d = std::cos(dir[d]);
                    double sin_d = std::sin(dir[d]);
                    double E_freq = E * df[f] * dd[d];
                    
                    Sxx += E_freq * (n * (cos_d * cos_d + 1) - 0.5);
                    Sxy += E_freq * n * cos_d * sin_d;
                    Syy += E_freq * (n * (sin_d * sin_d + 1) - 0.5);
                }
            }
            
            Hs_field[idx2] = 4.0 * std::sqrt(m0);
            Tp_field[idx2] = (freq[peak_f] > 0) ? 1.0 / freq[peak_f] : 10.0;
            dir_field[idx2] = std::atan2(Ey, Ex) * RAD_TO_DEG;
            if (dir_field[idx2] < 0) dir_field[idx2] += 360.0;
            
            Sxx_field[idx2] = RHO_WATER * G * Sxx;
            Sxy_field[idx2] = RHO_WATER * G * Sxy;
            Syy_field[idx2] = RHO_WATER * G * Syy;
        }
    }
}

WaveState WaveModel::getWaveState(int i, int j) const {
    WaveState state;
    int idx2 = idx2d(i, j);
    
    state.Hs = Hs_field[idx2];
    state.Tp = Tp_field[idx2];
    state.dir_mean = dir_field[idx2];
    state.Sxx = Sxx_field[idx2];
    state.Sxy = Sxy_field[idx2];
    state.Syy = Syy_field[idx2];
    
    return state;
}

void WaveModel::getRadiationStress(std::vector<double>& Sxx,
                                   std::vector<double>& Sxy,
                                   std::vector<double>& Syy) const {
    Sxx = Sxx_field;
    Sxy = Sxy_field;
    Syy = Syy_field;
}

void WaveModel::getStokesDrift(std::vector<double>& u_stokes,
                               std::vector<double>& v_stokes) const {
    u_stokes.resize(nx * ny);
    v_stokes.resize(nx * ny);
    
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int idx2 = idx2d(i, j);
            if (mask[idx2] == 0) {
                u_stokes[idx2] = v_stokes[idx2] = 0.0;
                continue;
            }
            
            double Us = 0.0, Vs = 0.0;
            double h = depth[idx2];
            
            for (int f = 0; f < nf; ++f) {
                double omega = 2.0 * M_PI * freq[f];
                double k = dispersion(freq[f], h);
                
                for (int d = 0; d < nd; ++d) {
                    double E = action[idx(i, j, f, d)] * omega * df[f] * dd[d];
                    
                    // Stokes drift: U_s = ωk E / (ρg)
                    double U_s_mag = omega * k * E / (RHO_WATER * G);
                    Us += U_s_mag * std::cos(dir[d]);
                    Vs += U_s_mag * std::sin(dir[d]);
                }
            }
            
            u_stokes[idx2] = Us;
            v_stokes[idx2] = Vs;
        }
    }
}

void WaveModel::getWaveBottomStress(std::vector<double>& tau_w) const {
    tau_w.resize(nx * ny);
    
    // Wave friction factor (typical range: 0.01-0.03)
    // Value of 0.015 is representative for sandy/muddy bottoms
    const double WAVE_FRICTION_FACTOR = 0.015;
    
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int idx2 = idx2d(i, j);
            if (mask[idx2] == 0) {
                tau_w[idx2] = 0.0;
                continue;
            }
            
            double h = depth[idx2];
            double Hs = Hs_field[idx2];
            double Tp = Tp_field[idx2];
            
            if (Tp <= 0.0 || Hs <= 0.0) {
                tau_w[idx2] = 0.0;
                continue;
            }
            
            // Compute representative wave orbital velocity at bottom
            double omega = 2.0 * M_PI / Tp;
            double k = dispersion(1.0 / Tp, h);
            
            // RMS orbital velocity amplitude at bottom
            double U_orb = 0.5 * Hs * omega / std::sinh(k * h);
            
            // Bottom stress from wave friction
            // τ_w = (1/2) * ρ * f_w * U_orb²
            tau_w[idx2] = 0.5 * RHO_WATER * WAVE_FRICTION_FACTOR * U_orb * U_orb;
        }
    }
}

// =============================================================================
// InternalWaveModel Implementation
// =============================================================================

InternalWaveModel::InternalWaveModel() 
    : nx(0), ny(0), nz(0), dx(0), dy(0), num_modes(0), coriolis_f(1e-4) {}

void InternalWaveModel::initialize(int nx_, int ny_, int nz_,
                                   double dx_, double dy_,
                                   const std::vector<double>& z_levels,
                                   const std::vector<double>& N2_profile) {
    nx = nx_;
    ny = ny_;
    nz = nz_;
    dx = dx_;
    dy = dy_;
    
    z = z_levels;
    N2 = N2_profile;
    
    eta_iw.resize(nx * ny * nz, 0.0);
    u_iw.resize(nx * ny * nz, 0.0);
    v_iw.resize(nx * ny * nz, 0.0);
}

void InternalWaveModel::setStratification(const std::vector<double>& temperature,
                                          const std::vector<double>& salinity) {
    // Compute N² from T, S profiles
    // N² = -(g/ρ) * dρ/dz
    
    N2.resize(nz);
    double g = 9.81;
    double rho0 = 1025.0;
    double alpha = 2.0e-4;  // Thermal expansion
    double beta = 7.5e-4;   // Haline contraction
    
    for (int k = 1; k < nz - 1; ++k) {
        double dT_dz = (temperature[k + 1] - temperature[k - 1]) / (z[k + 1] - z[k - 1]);
        double dS_dz = (salinity[k + 1] - salinity[k - 1]) / (z[k + 1] - z[k - 1]);
        
        // N² = g * (α * dT/dz - β * dS/dz)
        N2[k] = g * (alpha * dT_dz - beta * dS_dz);
        N2[k] = std::max(1e-8, N2[k]);  // Ensure stability
    }
    
    N2[0] = N2[1];
    N2[nz - 1] = N2[nz - 2];
}

void InternalWaveModel::computeNormalModes(int n_modes) {
    num_modes = n_modes;
    mode_shapes.resize(n_modes);
    mode_speeds.resize(n_modes);
    
    // Solve eigenvalue problem for vertical modes:
    // d²φ/dz² + (N²/c²) φ = 0
    // with φ = 0 at surface and bottom
    
    // Use finite difference approximation
    std::vector<double> dz(nz - 1);
    for (int k = 0; k < nz - 1; ++k) {
        dz[k] = z[k + 1] - z[k];
    }
    
    // Simple estimation based on mean N²
    double N2_mean = 0.0;
    for (int k = 0; k < nz; ++k) {
        N2_mean += N2[k];
    }
    N2_mean /= nz;
    
    double H = std::abs(z[0] - z[nz - 1]);
    
    for (int m = 0; m < n_modes; ++m) {
        // Mode m has (m+1) zero crossings
        mode_speeds[m] = std::sqrt(N2_mean) * H / ((m + 1) * M_PI);
        
        // Simple sinusoidal mode shape
        mode_shapes[m].resize(nz);
        for (int k = 0; k < nz; ++k) {
            double zeta = (z[k] - z[nz - 1]) / H;  // 0 at bottom, 1 at surface
            mode_shapes[m][k] = std::sin((m + 1) * M_PI * zeta);
        }
    }
}

void InternalWaveModel::getModeShape(int mode, std::vector<double>& phi) const {
    if (mode >= 0 && mode < num_modes) {
        phi = mode_shapes[mode];
    }
}

double InternalWaveModel::getModeSpeed(int mode) const {
    if (mode >= 0 && mode < num_modes) {
        return mode_speeds[mode];
    }
    return 0.0;
}

double InternalWaveModel::estimateTidalConversion(
    const std::vector<double>& bathymetry,
    double U_tide, double omega_tide) const {
    // Simplified tidal conversion estimate
    // F ∝ ρ N U² |∇h|² / ω
    
    double F_total = 0.0;
    double N_mean = std::sqrt(std::accumulate(N2.begin(), N2.end(), 0.0) / N2.size());
    
    // Would need bathymetry gradient calculation...
    // Simplified return
    return 0.0;
}

void InternalWaveModel::estimateMixing(std::vector<double>& Kv_iw) const {
    // Mixing from internal wave breaking
    // Osborn (1980): Kv = Γ ε / N²
    // where Γ ≈ 0.2 is mixing efficiency
    
    double Gamma = 0.2;
    double epsilon_bg = 1e-10;  // Background dissipation rate
    
    Kv_iw.resize(nz);
    for (int k = 0; k < nz; ++k) {
        Kv_iw[k] = Gamma * epsilon_bg / std::max(N2[k], 1e-8);
        Kv_iw[k] = std::max(1e-6, Kv_iw[k]);  // Minimum diffusivity
    }
}

// =============================================================================
// OceanAcousticsModel Implementation
// =============================================================================

OceanAcousticsModel::OceanAcousticsModel()
    : nr(0), nz(0), dr(0), dz(0), range_max(0), depth_max(0),
      rho_bottom(1.5), c_bottom(1.1), alpha_bottom(0.5),
      z_source(50.0), frequency(100.0) {}

void OceanAcousticsModel::initialize(double r_max, double d_max, int n_r, int n_z) {
    range_max = r_max;
    depth_max = d_max;
    nr = n_r;
    nz = n_z;
    dr = r_max / (nr - 1);
    dz = d_max / (nz - 1);
    
    TL.resize(nr * nz, 0.0);
}

void OceanAcousticsModel::setSoundSpeedProfile(const std::vector<double>& z_in,
                                               const std::vector<double>& c_in) {
    z_levels = z_in;
    c_profile = c_in;
}

double OceanAcousticsModel::computeSoundSpeed(double T, double S, double z) {
    // Mackenzie (1981) equation
    double c = 1448.96 + 4.591 * T - 5.304e-2 * T * T + 2.374e-4 * T * T * T
             + 1.340 * (S - 35.0) + 1.630e-2 * z + 1.675e-7 * z * z
             - 1.025e-2 * T * (S - 35.0) - 7.139e-13 * T * z * z * z;
    return c;
}

void OceanAcousticsModel::setBottomProperties(double density_ratio,
                                              double sound_speed_ratio,
                                              double attenuation) {
    rho_bottom = density_ratio;
    c_bottom = sound_speed_ratio;
    alpha_bottom = attenuation;
}

void OceanAcousticsModel::setSource(double source_depth, double freq) {
    z_source = source_depth;
    frequency = freq;
}

double OceanAcousticsModel::getSoundSpeed(double depth) const {
    if (z_levels.empty()) return 1500.0;  // Default
    
    // Linear interpolation
    for (size_t i = 0; i < z_levels.size() - 1; ++i) {
        if (depth >= z_levels[i] && depth <= z_levels[i + 1]) {
            double alpha = (depth - z_levels[i]) / (z_levels[i + 1] - z_levels[i]);
            return (1 - alpha) * c_profile[i] + alpha * c_profile[i + 1];
        }
    }
    
    return c_profile.back();
}

void OceanAcousticsModel::computeRayTracing(int num_rays) {
    // Simple ray tracing implementation
    double theta_min = -20.0 * DEG_TO_RAD;
    double theta_max = 20.0 * DEG_TO_RAD;
    double d_theta = (theta_max - theta_min) / (num_rays - 1);
    
    ray_paths.clear();
    ray_paths.resize(num_rays);
    
    for (int ray = 0; ray < num_rays; ++ray) {
        double theta = theta_min + ray * d_theta;
        
        double x = 0.0, z = z_source;
        double c0 = getSoundSpeed(z_source);
        double xi = std::cos(theta) / c0;  // Snell's law parameter
        
        ray_paths[ray].clear();
        ray_paths[ray].push_back({x, z});
        
        double ds = 10.0;  // Step size
        
        while (x < range_max && z > 0 && z < depth_max) {
            double c = getSoundSpeed(z);
            double c_grad = (getSoundSpeed(z + 1) - getSoundSpeed(z - 1)) / 2.0;
            
            double cos_theta = c * xi;
            if (std::abs(cos_theta) > 1.0) break;  // Total reflection
            
            double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
            if (theta < 0) sin_theta = -sin_theta;
            
            x += ds * cos_theta;
            z += ds * sin_theta;
            
            // Update angle using Snell's law
            theta = std::asin(sin_theta);
            
            ray_paths[ray].push_back({x, z});
            
            // Bottom reflection
            if (z >= depth_max) {
                z = 2 * depth_max - z;
                theta = -theta;
            }
            // Surface reflection
            if (z <= 0) {
                z = -z;
                theta = -theta;
            }
        }
    }
    
    // Compute transmission loss from ray density
    // (simplified - would need proper amplitude calculation)
    for (int ir = 0; ir < nr; ++ir) {
        double r = ir * dr;
        for (int iz = 0; iz < nz; ++iz) {
            double z = iz * dz;
            
            // Spreading loss
            double TL_spread = (r > 1.0) ? 10.0 * std::log10(r) + 10.0 * std::log10(z + 1) : 0.0;
            
            // Absorption (frequency-dependent)
            double alpha_freq = 0.003 + 0.1 * frequency / 1000.0 * frequency / 1000.0;
            double TL_absorption = alpha_freq * r / 1000.0;
            
            TL[ir + iz * nr] = TL_spread + TL_absorption;
        }
    }
}

double OceanAcousticsModel::getTransmissionLoss(double range, double depth) const {
    int ir = static_cast<int>(range / dr);
    int iz = static_cast<int>(depth / dz);
    
    ir = std::max(0, std::min(ir, nr - 1));
    iz = std::max(0, std::min(iz, nz - 1));
    
    return TL[ir + iz * nr];
}

double OceanAcousticsModel::estimateDetectionRange(double source_level,
                                                   double noise_level,
                                                   double detection_threshold) const {
    // Find range where SL - TL = NL + DT
    double target_TL = source_level - noise_level - detection_threshold;
    
    for (int ir = nr - 1; ir >= 0; --ir) {
        int iz = static_cast<int>(z_source / dz);
        if (TL[ir + iz * nr] < target_TL) {
            return ir * dr;
        }
    }
    
    return 0.0;
}

// =============================================================================
// SedimentTransportModel Implementation
// =============================================================================

SedimentTransportModel::SedimentTransportModel()
    : nx(0), ny(0), dx(0), dy(0),
      d50(0.0002), d90(0.0005), porosity(0.4), rho_s(2650.0) {}

void SedimentTransportModel::initialize(int nx_, int ny_, double dx_, double dy_) {
    nx = nx_;
    ny = ny_;
    dx = dx_;
    dy = dy_;
    
    int n = nx * ny;
    qb_x.resize(n, 0.0);
    qb_y.resize(n, 0.0);
    qs_x.resize(n, 0.0);
    qs_y.resize(n, 0.0);
    z_bed.resize(n, 0.0);
    dz_bed.resize(n, 0.0);
}

void SedimentTransportModel::setSedimentProperties(double d50_, double d90_,
                                                   double porosity_, double density_) {
    d50 = d50_;
    d90 = d90_;
    porosity = porosity_;
    rho_s = density_;
}

double SedimentTransportModel::computeShieldsParameter(double tau_b, double d) const {
    return tau_b / ((rho_s - RHO_WATER) * G * d);
}

double SedimentTransportModel::criticalShields(double d) const {
    // Soulsby-Whitehouse formula
    double D_star = d * std::pow((rho_s / RHO_WATER - 1.0) * G / (1e-6 * 1e-6), 1.0/3.0);
    return 0.30 / (1.0 + 1.2 * D_star) + 0.055 * (1.0 - std::exp(-0.020 * D_star));
}

void SedimentTransportModel::computeTransport(const std::vector<double>& u,
                                              const std::vector<double>& v,
                                              const std::vector<double>& h,
                                              const WaveModel* waves) {
    std::vector<double> tau_w(nx * ny, 0.0);
    if (waves != nullptr) {
        waves->getWaveBottomStress(tau_w);
    }
    
    computeSoulsbyVanRijn(u, v, h, tau_w);
}

void SedimentTransportModel::computeSoulsbyVanRijn(
    const std::vector<double>& u,
    const std::vector<double>& v,
    const std::vector<double>& h,
    const std::vector<double>& tau_w) {
    
    double Asb = 0.005 * h[0] * std::pow(d50 / h[0] / (1.2 * 0.4), 1.2);
    double Ass = 0.012 * d50 * std::pow(d50 * 0.4, -0.6);
    
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int idx = i + j * nx;
            
            double speed = std::sqrt(u[idx] * u[idx] + v[idx] * v[idx]);
            if (speed < 0.01) {
                qb_x[idx] = qb_y[idx] = 0.0;
                qs_x[idx] = qs_y[idx] = 0.0;
                continue;
            }
            
            // Critical velocity
            double D_star = d50 * std::pow((rho_s / RHO_WATER - 1.0) * G / (1e-6 * 1e-6), 1.0/3.0);
            double U_cr;
            if (D_star <= 4.0) {
                U_cr = 0.19 * std::pow(D_star, 0.1) * std::log10(4.0 * h[idx] / d90);
            } else {
                U_cr = 8.5 * std::pow(D_star, 0.6) * std::log10(4.0 * h[idx] / d90);
            }
            
            // Combined wave-current velocity
            double U_cw = speed;
            if (tau_w[idx] > 0) {
                U_cw = std::sqrt(speed * speed + 0.018 * tau_w[idx] / RHO_WATER);
            }
            
            if (U_cw < U_cr) {
                qb_x[idx] = qb_y[idx] = 0.0;
                qs_x[idx] = qs_y[idx] = 0.0;
                continue;
            }
            
            // Total transport
            double q_total = Asb * std::pow(std::max(0.0, U_cw - U_cr), 2.4) +
                            Ass * std::pow(std::max(0.0, U_cw - U_cr), 2.4);
            
            // Direction
            double dir_factor_x = u[idx] / speed;
            double dir_factor_y = v[idx] / speed;
            
            qb_x[idx] = 0.5 * q_total * dir_factor_x;
            qb_y[idx] = 0.5 * q_total * dir_factor_y;
            qs_x[idx] = 0.5 * q_total * dir_factor_x;
            qs_y[idx] = 0.5 * q_total * dir_factor_y;
        }
    }
}

void SedimentTransportModel::updateBed(double dt) {
    // Bed update from Exner equation: ∂z_b/∂t = -1/(1-p) * ∇·q
    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {
            int idx = i + j * nx;
            
            double dqx_dx = (qb_x[i + 1 + j * nx] + qs_x[i + 1 + j * nx] -
                            qb_x[i - 1 + j * nx] - qs_x[i - 1 + j * nx]) / (2.0 * dx);
            double dqy_dy = (qb_y[i + (j + 1) * nx] + qs_y[i + (j + 1) * nx] -
                            qb_y[i + (j - 1) * nx] - qs_y[i + (j - 1) * nx]) / (2.0 * dy);
            
            dz_bed[idx] = -dt / (1.0 - porosity) * (dqx_dx + dqy_dy);
            z_bed[idx] += dz_bed[idx];
        }
    }
}

double SedimentTransportModel::computeLongshoreTransport(double Hs, double Tp,
                                                         double wave_angle,
                                                         double depth) const {
    // CERC formula
    double K = 0.39;  // Empirical coefficient
    double rho = RHO_WATER;
    double g = G;
    
    double alpha_b = wave_angle * DEG_TO_RAD;  // Wave angle at breaking
    double Hb = 0.78 * depth;  // Breaking wave height
    
    double P_ls = (rho * g * Hb * Hb) / 16.0 * std::sqrt(g * Hb / 0.78) *
                 std::sin(2.0 * alpha_b);
    
    double Q = K * P_ls / ((rho_s - rho) * g * (1.0 - porosity));
    
    return Q;  // m³/s
}

// =============================================================================
// AirSeaInteraction Implementation
// =============================================================================

void AirSeaInteraction::computeWindStressCOARE(
    const std::vector<double>& u10,
    const std::vector<double>& v10,
    const std::vector<double>& T_air,
    const std::vector<double>& T_sea,
    const std::vector<double>& humidity,
    std::vector<double>& tau_x,
    std::vector<double>& tau_y) {
    
    size_t n = u10.size();
    tau_x.resize(n);
    tau_y.resize(n);
    
    for (size_t i = 0; i < n; ++i) {
        double wind_speed = std::sqrt(u10[i] * u10[i] + v10[i] * v10[i]);
        
        // Neutral drag coefficient
        double Cdn;
        if (wind_speed < 10.0) {
            Cdn = 1.0e-3 * (0.61 + 0.063 * wind_speed);
        } else {
            Cdn = 1.0e-3 * (0.61 + 0.063 * 10.0 + 0.066 * (wind_speed - 10.0));
        }
        
        // Stability correction (simplified Monin-Obukhov)
        double dT = T_sea[i] - T_air[i];
        double Cd;
        if (dT > 0) {
            // Unstable
            Cd = Cdn * (1.0 + 0.02 * dT);
        } else {
            // Stable
            Cd = Cdn * (1.0 + 0.01 * dT);
        }
        
        Cd = std::max(0.5e-3, std::min(Cd, 3.0e-3));
        
        tau_x[i] = RHO_AIR * Cd * wind_speed * u10[i];
        tau_y[i] = RHO_AIR * Cd * wind_speed * v10[i];
    }
}

void AirSeaInteraction::computeHeatFluxes(
    const std::vector<double>& T_air,
    const std::vector<double>& T_sea,
    const std::vector<double>& u10,
    const std::vector<double>& humidity,
    const std::vector<double>& cloud_cover,
    const std::vector<double>& solar_zenith,
    std::vector<double>& Q_sw,
    std::vector<double>& Q_lw,
    std::vector<double>& Q_sens,
    std::vector<double>& Q_lat) {
    
    size_t n = T_air.size();
    Q_sw.resize(n);
    Q_lw.resize(n);
    Q_sens.resize(n);
    Q_lat.resize(n);
    
    double S0 = 1361.0;  // Solar constant (W/m²)
    
    for (size_t i = 0; i < n; ++i) {
        double wind_speed = std::sqrt(u10[i] * u10[i]);
        
        // Shortwave radiation
        double cos_zenith = std::cos(solar_zenith[i] * DEG_TO_RAD);
        if (cos_zenith > 0) {
            double atmos_trans = 0.7;
            double cloud_factor = 1.0 - 0.65 * cloud_cover[i];
            Q_sw[i] = S0 * cos_zenith * atmos_trans * cloud_factor * (1.0 - 0.06);  // 6% albedo
        } else {
            Q_sw[i] = 0.0;
        }
        
        // Longwave radiation (Berliand formula)
        double T_sea_K = T_sea[i] + 273.15;
        double e_a = humidity[i] * 6.11 * std::exp(17.27 * T_air[i] / (T_air[i] + 237.3));
        double epsilon = 0.97;
        Q_lw[i] = -epsilon * STEFAN_BOLTZMANN * std::pow(T_sea_K, 4) *
                 (0.39 - 0.05 * std::sqrt(e_a)) * (1.0 - 0.7 * cloud_cover[i]);
        
        // Sensible heat flux (bulk formula)
        double Ch = 1.0e-3 * (1.0 + 0.01 * wind_speed);
        Q_sens[i] = RHO_AIR * CP_WATER * Ch * wind_speed * (T_sea[i] - T_air[i]);
        
        // Latent heat flux
        double q_sat = 0.622 * 6.11 * std::exp(17.27 * T_sea[i] / (T_sea[i] + 237.3)) / 1013.25;
        double q_air = humidity[i] * 0.622 * 6.11 * std::exp(17.27 * T_air[i] / (T_air[i] + 237.3)) / 1013.25;
        double Ce = Ch;  // Dalton number ≈ Stanton number
        Q_lat[i] = RHO_AIR * L_EVAP * Ce * wind_speed * (q_sat - q_air);
    }
}

double AirSeaInteraction::waveDependentDrag(double u10, double wave_age) {
    // Wave age = c_p / u_*
    // Young waves: high drag, old waves: lower drag
    
    double Cd_0 = 1.0e-3 * (0.61 + 0.063 * u10);
    
    if (wave_age < 20) {
        // Young, developing waves
        return Cd_0 * (1.0 + 0.1 * (20.0 - wave_age) / 20.0);
    } else {
        // Mature waves
        return Cd_0;
    }
}

// =============================================================================
// OceanMixingModel Implementation
// =============================================================================

OceanMixingModel::OceanMixingModel() {}

void OceanMixingModel::computeKPP(const std::vector<double>& T,
                                  const std::vector<double>& S,
                                  const std::vector<double>& u,
                                  const std::vector<double>& v,
                                  const std::vector<double>& tau_x,
                                  const std::vector<double>& tau_y,
                                  const std::vector<double>& Q_net,
                                  double f,
                                  std::vector<double>& Kv,
                                  std::vector<double>& Kt,
                                  std::vector<double>& Ks) {
    size_t n = T.size();
    Kv.resize(n);
    Kt.resize(n);
    Ks.resize(n);
    
    // Simplified KPP implementation
    // Full implementation would require boundary layer depth calculation
    // and non-local flux computation
    
    for (size_t i = 0; i < n; ++i) {
        // Background diffusivity
        Kv[i] = 1e-4;
        Kt[i] = 1e-5;
        Ks[i] = 1e-5;
        
        // Would add boundary layer enhancement here...
    }
}

double OceanMixingModel::computeBoundaryLayerDepth(
    const std::vector<double>& T,
    const std::vector<double>& S,
    double u_star, double B_f, double f) const {
    
    // Simplified boundary layer depth estimation
    // Full KPP uses Ri_bulk criterion
    
    double w_s = std::sqrt(u_star * u_star + 0.1 * B_f);
    double h_bl = 0.7 * u_star / std::abs(f);
    
    return std::max(10.0, h_bl);
}

double OceanMixingModel::langmuirEnhancement(double u_star, double La_t) const {
    // Enhancement factor for Langmuir turbulence
    // Based on McWilliams et al. (1997)
    
    if (La_t < 0.3) {
        return 1.0 + 2.0 * (0.3 - La_t) / 0.3;
    }
    return 1.0;
}

double OceanMixingModel::computeLangmuirNumber(double u_star, double U_stokes) {
    // Turbulent Langmuir number La_t = sqrt(u_*/U_s)
    if (U_stokes > 0.001) {
        return std::sqrt(u_star / U_stokes);
    }
    return 10.0;  // Large number = no Langmuir effect
}

// =============================================================================
// OceanPhysicsUtils Implementation
// =============================================================================

namespace OceanPhysicsUtils {

double coriolisParameter(double latitude) {
    return 2.0 * OMEGA * std::sin(latitude * DEG_TO_RAD);
}

double ekmanDepth(double f, double Az) {
    return M_PI * std::sqrt(2.0 * Az / std::abs(f));
}

double rossbyRadius(double N, double H, double f) {
    return N * H / std::abs(f);
}

double thermalExpansion(double T, double S, double p) {
    // Approximate thermal expansion coefficient
    return 2.0e-4 * (1.0 + 0.01 * (T - 10.0));
}

double halineContraction(double T, double S, double p) {
    // Approximate haline contraction coefficient
    return 7.5e-4;
}

double potentialTemperature(double T, double S, double p, double p_ref) {
    // Adiabatic lapse rate ≈ 0.15°C per 1000 dbar
    double gamma = 0.15e-3;
    return T - gamma * (p - p_ref);
}

double saturationHumidity(double T, double p) {
    // Saturation specific humidity
    double e_sat = 6.11 * std::exp(17.27 * T / (T + 237.3));
    return 0.622 * e_sat / (p - 0.378 * e_sat);
}

double deepWaterWavelength(double T) {
    return G * T * T / (2.0 * M_PI);
}

double wavelength(double T, double h, double tol) {
    double L_deep = deepWaterWavelength(T);
    double L = L_deep;
    
    // Newton iteration
    for (int i = 0; i < 50; ++i) {
        double k = 2.0 * M_PI / L;
        double F = L - L_deep * std::tanh(k * h);
        double dF = 1.0 + L_deep * k * h / (std::cosh(k * h) * std::cosh(k * h)) / L;
        
        double dL = F / dF;
        L -= dL;
        
        if (std::abs(dL) < tol * L) break;
    }
    
    return L;
}

double ursellNumber(double H, double L, double h) {
    return H * L * L / (h * h * h);
}

double waveSetup(double Hs, double T, double beach_slope) {
    // Simplified wave setup estimate
    // Full calculation requires radiation stress gradient
    double L = wavelength(T, 10.0);  // Assume 10m depth
    double k = 2.0 * M_PI / L;
    
    // Approximate setup at shoreline
    return 0.15 * Hs;
}

} // namespace OceanPhysicsUtils

} // namespace FSRM
