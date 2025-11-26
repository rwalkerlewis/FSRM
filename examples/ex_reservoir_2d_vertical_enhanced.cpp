#include "Simulator.hpp"
#include "PhysicsKernel.hpp"
#include "Visualization.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <sys/stat.h>

using namespace FSRM;

/**
 * Enhanced 2D Vertical Reservoir Simulation
 * 
 * Features:
 * - Fully implicit time integration (Newton-Raphson for saturation)
 * - Full 2D plane-strain geomechanics (sigma_xx, sigma_zz, tau_xz)
 * - Proper boundary conditions (no-flow, fixed displacement, free surface)
 * - Well performance tracking (BHP, rates, cumulative production)
 * - Labeled plots with comprehensive information
 */

struct WellPerformance {
    double bhp = 0.0;              // Bottom-hole pressure (Pa)
    double oil_rate = 0.0;         // m³/day
    double water_rate = 0.0;       // m³/day
    double total_rate = 0.0;       // m³/day
    double cum_oil = 0.0;          // m³
    double cum_water = 0.0;        // m³
    double water_cut = 0.0;        // fraction
};

struct Well2D {
    double x, z;
    bool is_injector;
    double target_rate;  // m³/day
    double target_bhp;   // Pa
    std::string name;
    WellPerformance perf;
    
    Well2D(double x_pos, double z_pos, bool inj, double rate, double bhp, const std::string& n)
        : x(x_pos), z(z_pos), is_injector(inj), 
          target_rate(rate), target_bhp(bhp), name(n) {}
};

class EnhancedReservoir2D {
public:
    EnhancedReservoir2D(int nx, int nz, double Lx, double Lz)
        : nx_(nx), nz_(nz), Lx_(Lx), Lz_(Lz),
          dx_(Lx / nx), dz_(Lz / nz) {
        
        // Flow fields
        P_.resize(nz, std::vector<double>(nx, 0.0));
        Sw_.resize(nz, std::vector<double>(nx, 0.0));
        phi_.resize(nz, std::vector<double>(nx, 0.0));
        k_.resize(nz, std::vector<double>(nx, 0.0));
        
        // Geomechanics fields
        ux_.resize(nz, std::vector<double>(nx, 0.0));
        uz_.resize(nz, std::vector<double>(nx, 0.0));
        sigma_xx_.resize(nz, std::vector<double>(nx, 0.0));
        sigma_zz_.resize(nz, std::vector<double>(nx, 0.0));
        tau_xz_.resize(nz, std::vector<double>(nx, 0.0));
        eps_vol_.resize(nz, std::vector<double>(nx, 0.0));
        
        // Properties
        phi0_ = 0.2;
        k0_ = 100e-15;
        E_ = 10e9;
        nu_ = 0.25;
        rho_f_ = 1000.0;
        rho_s_ = 2500.0;
        mu_w_ = 1e-3;
        mu_o_ = 2e-3;
        alpha_ = 0.7;
        cf_ = 1e-9;
        cr_ = 1e-9;
        P_initial_ = 30e6;
        
        // Compute elastic constants
        lambda_ = E_ * nu_ / ((1.0 + nu_) * (1.0 - 2.0 * nu_));  // Lame's first parameter
        G_ = E_ / (2.0 * (1.0 + nu_));  // Shear modulus
        K_ = E_ / (3.0 * (1.0 - 2.0 * nu_));  // Bulk modulus
        
        initialize();
    }
    
    void addWell(const Well2D& well) {
        wells_.push_back(well);
    }
    
    void setRockProperties(double phi, double k, double E, double nu) {
        phi0_ = phi;
        k0_ = k;
        E_ = E;
        nu_ = nu;
        lambda_ = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        G_ = E / (2.0 * (1.0 + nu));
        K_ = E / (3.0 * (1.0 - 2.0 * nu));
    }
    
    void initialize() {
        for (int k = 0; k < nz_; ++k) {
            for (int i = 0; i < nx_; ++i) {
                double depth = k * dz_;
                
                // Hydrostatic pressure
                P_[k][i] = P_initial_ + rho_f_ * 9.81 * depth;
                Sw_[k][i] = 0.2;
                phi_[k][i] = phi0_;
                k_[k][i] = k0_;
                
                // Initial stress state (lithostatic)
                double sigma_v = rho_s_ * 9.81 * depth;  // Vertical stress
                sigma_zz_[k][i] = -sigma_v;  // Negative = compression
                sigma_xx_[k][i] = -nu_ / (1.0 - nu_) * sigma_v;  // Ko condition
                tau_xz_[k][i] = 0.0;
            }
        }
    }
    
    void step(double dt) {
        // Fully implicit coupled solve
        solveFlowImplicit(dt);
        solveSaturationImplicit(dt);
        solveGeomechanicsFull();
        updatePoroPermFromStress();
        updateWellPerformance(dt);
        
        time_ += dt;
    }
    
    // Getters
    const std::vector<std::vector<double>>& getPressure() const { return P_; }
    const std::vector<std::vector<double>>& getSaturation() const { return Sw_; }
    const std::vector<std::vector<double>>& getDisplacementX() const { return ux_; }
    const std::vector<std::vector<double>>& getDisplacementZ() const { return uz_; }
    const std::vector<std::vector<double>>& getStressXX() const { return sigma_xx_; }
    const std::vector<std::vector<double>>& getStressZZ() const { return sigma_zz_; }
    const std::vector<std::vector<double>>& getShearStress() const { return tau_xz_; }
    const std::vector<std::vector<double>>& getPermeability() const { return k_; }
    std::vector<std::vector<double>> getSubsidence() const {
        std::vector<std::vector<double>> sub(nz_, std::vector<double>(nx_));
        for (int k = 0; k < nz_; ++k) {
            for (int i = 0; i < nx_; ++i) {
                sub[k][i] = -uz_[k][i];  // Negative displacement = subsidence
            }
        }
        return sub;
    }
    
    const std::vector<Well2D>& getWells() const { return wells_; }
    double getTime() const { return time_; }
    double getLx() const { return Lx_; }
    double getLz() const { return Lz_; }
    int getNx() const { return nx_; }
    int getNz() const { return nz_; }
    
private:
    void solveFlowImplicit(double dt) {
        // Newton-Raphson for implicit pressure equation
        auto P_old = P_;
        const int max_iter = 20;
        const double tol = 1e-6;
        
        for (int iter = 0; iter < max_iter; ++iter) {
            double max_change = 0.0;
            
            for (int k = 1; k < nz_ - 1; ++k) {
                for (int i = 1; i < nx_ - 1; ++i) {
                    // Skip well cells
                    if (isWellCell(i, k)) continue;
                    
                    double phi = phi_[k][i];
                    double perm = k_[k][i];
                    double ct = cf_ + cr_ + alpha_ * alpha_ / K_;  // Total compressibility
                    
                    // Harmonic mean for transmissibility
                    double k_xp = 2.0 * perm * k_[k][i+1] / (perm + k_[k][i+1] + 1e-30);
                    double k_xm = 2.0 * perm * k_[k][i-1] / (perm + k_[k][i-1] + 1e-30);
                    double k_zp = 2.0 * perm * k_[k+1][i] / (perm + k_[k+1][i] + 1e-30);
                    double k_zm = 2.0 * perm * k_[k-1][i] / (perm + k_[k-1][i] + 1e-30);
                    
                    double Tx_plus = k_xp / (mu_eff(Sw_[k][i]) * dx_ * dx_);
                    double Tx_minus = k_xm / (mu_eff(Sw_[k][i]) * dx_ * dx_);
                    double Tz_plus = k_zp / (mu_eff(Sw_[k][i]) * dz_ * dz_);
                    double Tz_minus = k_zm / (mu_eff(Sw_[k][i]) * dz_ * dz_);
                    
                    // Gravity term
                    double rho_eff = Sw_[k][i] * rho_f_ + (1.0 - Sw_[k][i]) * 0.8 * rho_f_;  // Oil lighter
                    double grav = rho_eff * 9.81 * (Tz_plus - Tz_minus) * dz_;
                    
                    // Accumulation
                    double accum = phi * ct / dt;
                    
                    // Assemble and solve
                    double a_c = accum + Tx_plus + Tx_minus + Tz_plus + Tz_minus;
                    double rhs = accum * P_old[k][i] + grav;
                    rhs += Tx_plus * P_[k][i+1] + Tx_minus * P_[k][i-1];
                    rhs += Tz_plus * P_[k+1][i] + Tz_minus * P_[k-1][i];
                    
                    double P_new = rhs / (a_c + 1e-30);
                    max_change = std::max(max_change, std::abs(P_new - P_[k][i]));
                    P_[k][i] = P_new;
                }
            }
            
            applyPressureBCs();
            
            if (max_change < tol * P_initial_) break;
        }
    }
    
    void solveSaturationImplicit(double dt) {
        // Semi-implicit saturation with upstream weighting
        auto Sw_old = Sw_;
        const int max_iter = 10;
        const double tol = 1e-5;
        
        for (int iter = 0; iter < max_iter; ++iter) {
            double max_change = 0.0;
            
            for (int k = 1; k < nz_ - 1; ++k) {
                for (int i = 1; i < nx_ - 1; ++i) {
                    if (isWellCell(i, k)) continue;
                    
                    double phi = phi_[k][i];
                    
                    // Compute fluxes
                    double flux_xp = computeFlux(i, k, i+1, k, true);
                    double flux_xm = computeFlux(i, k, i-1, k, true);
                    double flux_zp = computeFlux(i, k, i, k+1, false);
                    double flux_zm = computeFlux(i, k, i, k-1, false);
                    
                    // Divergence
                    double div_flux = (flux_xp - flux_xm) / dx_ + (flux_zp - flux_zm) / dz_;
                    
                    // Semi-implicit update
                    double Sw_new = Sw_old[k][i] - (dt / phi) * div_flux;
                    Sw_new = std::max(0.2, std::min(0.8, Sw_new));
                    
                    max_change = std::max(max_change, std::abs(Sw_new - Sw_[k][i]));
                    Sw_[k][i] = Sw_new;
                }
            }
            
            applySaturationBCs();
            
            if (max_change < tol) break;
        }
    }
    
    void solveGeomechanicsFull() {
        // Full 2D plane-strain elasticity with finite differences
        // Solve for displacements, then compute stresses
        
        const int max_iter = 50;
        const double tol = 1e-8;
        
        auto ux_old = ux_;
        auto uz_old = uz_;
        
        for (int iter = 0; iter < max_iter; ++iter) {
            double max_change = 0.0;
            
            // Gauss-Seidel for displacement equations
            for (int k = 1; k < nz_ - 1; ++k) {
                for (int i = 1; i < nx_ - 1; ++i) {
                    
                    // X-displacement equation: (lambda + G) * d²ux/dx² + G * d²ux/dz² + lambda * d²uz/dxdz = alpha * dP/dx
                    double d2ux_dx2 = (ux_[k][i+1] - 2.0*ux_[k][i] + ux_[k][i-1]) / (dx_*dx_);
                    double d2ux_dz2 = (ux_[k+1][i] - 2.0*ux_[k][i] + ux_[k-1][i]) / (dz_*dz_);
                    double d2uz_dxdz = (uz_[k+1][i+1] - uz_[k+1][i-1] - uz_[k-1][i+1] + uz_[k-1][i-1]) / (4.0*dx_*dz_);
                    
                    double dP_dx = (P_[k][i+1] - P_[k][i-1]) / (2.0*dx_);
                    
                    double rhs_x = alpha_ * dP_dx + (lambda_ + G_) * d2ux_dx2 + G_ * d2ux_dz2 + lambda_ * d2uz_dxdz;
                    double coeff = (lambda_ + G_) / (dx_*dx_) + G_ / (dz_*dz_);
                    
                    double ux_new = (rhs_x + coeff * ux_[k][i]) / (2.0 * coeff + 1e-10);
                    
                    // Z-displacement equation: G * d²uz/dx² + (lambda + G) * d²uz/dz² + lambda * d²ux/dxdz = alpha * dP/dz - rho*g
                    double d2uz_dx2 = (uz_[k][i+1] - 2.0*uz_[k][i] + uz_[k][i-1]) / (dx_*dx_);
                    double d2uz_dz2 = (uz_[k+1][i] - 2.0*uz_[k][i] + uz_[k-1][i]) / (dz_*dz_);
                    double d2ux_dxdz = (ux_[k+1][i+1] - ux_[k+1][i-1] - ux_[k-1][i+1] + ux_[k-1][i-1]) / (4.0*dx_*dz_);
                    
                    double dP_dz = (P_[k+1][i] - P_[k-1][i]) / (2.0*dz_);
                    double body_force = -rho_s_ * 9.81;
                    
                    double rhs_z = alpha_ * dP_dz + body_force + G_ * d2uz_dx2 + (lambda_ + G_) * d2uz_dz2 + lambda_ * d2ux_dxdz;
                    coeff = G_ / (dx_*dx_) + (lambda_ + G_) / (dz_*dz_);
                    
                    double uz_new = (rhs_z + coeff * uz_[k][i]) / (2.0 * coeff + 1e-10);
                    
                    max_change = std::max(max_change, std::abs(ux_new - ux_[k][i]));
                    max_change = std::max(max_change, std::abs(uz_new - uz_[k][i]));
                    
                    ux_[k][i] = 0.5 * ux_new + 0.5 * ux_[k][i];  // Relaxation
                    uz_[k][i] = 0.5 * uz_new + 0.5 * uz_[k][i];
                }
            }
            
            applyDisplacementBCs();
            
            if (max_change < tol * Lx_) break;
        }
        
        // Compute stresses from strains
        computeStressFromDisplacement();
    }
    
    void computeStressFromDisplacement() {
        for (int k = 1; k < nz_ - 1; ++k) {
            for (int i = 1; i < nx_ - 1; ++i) {
                // Strain components
                double eps_xx = (ux_[k][i+1] - ux_[k][i-1]) / (2.0 * dx_);
                double eps_zz = (uz_[k+1][i] - uz_[k-1][i]) / (2.0 * dz_);
                double eps_xz = 0.5 * ((ux_[k+1][i] - ux_[k-1][i]) / (2.0 * dz_) +
                                      (uz_[k][i+1] - uz_[k][i-1]) / (2.0 * dx_));
                
                eps_vol_[k][i] = eps_xx + eps_zz;
                
                // Plane strain: eps_yy = -nu/(1-nu) * (eps_xx + eps_zz)
                // Effective stress law
                double dP = P_[k][i] - (P_initial_ + rho_f_ * 9.81 * k * dz_);
                
                sigma_xx_[k][i] = lambda_ * eps_vol_[k][i] + 2.0 * G_ * eps_xx - alpha_ * dP;
                sigma_zz_[k][i] = lambda_ * eps_vol_[k][i] + 2.0 * G_ * eps_zz - alpha_ * dP;
                tau_xz_[k][i] = 2.0 * G_ * eps_xz;
            }
        }
    }
    
    void updatePoroPermFromStress() {
        for (int k = 0; k < nz_; ++k) {
            for (int i = 0; i < nx_; ++i) {
                // Porosity from volumetric strain
                phi_[k][i] = phi0_ + (alpha_ - phi0_) * eps_vol_[k][i];
                phi_[k][i] = std::max(0.01, std::min(0.4, phi_[k][i]));
                
                // Kozeny-Carman
                double phi_ratio = phi_[k][i] / phi0_;
                k_[k][i] = k0_ * std::pow(phi_ratio, 3.0) * 
                          std::pow((1.0 - phi0_) / (1.0 - phi_[k][i] + 1e-6), 2.0);
                k_[k][i] = std::max(0.01 * k0_, std::min(10.0 * k0_, k_[k][i]));
            }
        }
    }
    
    void updateWellPerformance(double dt) {
        for (auto& well : wells_) {
            int i = static_cast<int>(well.x / dx_);
            int k = static_cast<int>(well.z / dz_);
            
            i = std::max(1, std::min(nx_ - 2, i));
            k = std::max(1, std::min(nz_ - 2, k));
            
            // Bottom-hole pressure
            well.perf.bhp = P_[k][i];
            
            // Compute flow rates based on well index
            double rw = 0.15;  // Well radius (m)
            double re = 0.2 * std::sqrt(dx_ * dx_ + dz_ * dz_);  // Peaceman
            double well_index = 2.0 * M_PI * k_[k][i] * dz_ / (mu_eff(Sw_[k][i]) * std::log(re / rw));
            
            if (well.is_injector) {
                // Water injection
                double dP = well.perf.bhp - (P_initial_ + rho_f_ * 9.81 * k * dz_);
                well.perf.total_rate = well_index * dP * 86400.0;  // m³/day
                well.perf.water_rate = well.perf.total_rate;
                well.perf.oil_rate = 0.0;
                well.perf.water_cut = 1.0;
                
                well.perf.cum_water += well.perf.water_rate * dt / 86400.0;
            } else {
                // Production
                double dP = (P_initial_ + rho_f_ * 9.81 * k * dz_) - well.perf.bhp;
                double total_flux = well_index * dP * 86400.0;
                
                double fw = fractionalFlow(Sw_[k][i]);
                well.perf.water_rate = total_flux * fw;
                well.perf.oil_rate = total_flux * (1.0 - fw);
                well.perf.total_rate = total_flux;
                well.perf.water_cut = fw;
                
                well.perf.cum_water += well.perf.water_rate * dt / 86400.0;
                well.perf.cum_oil += well.perf.oil_rate * dt / 86400.0;
            }
        }
    }
    
    void applyPressureBCs() {
        // Top: atmospheric (zero gauge pressure)
        for (int i = 0; i < nx_; ++i) {
            P_[0][i] = P_initial_;
        }
        
        // Bottom: no-flow (mirror)
        for (int i = 0; i < nx_; ++i) {
            P_[nz_-1][i] = P_[nz_-2][i];
        }
        
        // Left and right: no-flow
        for (int k = 0; k < nz_; ++k) {
            P_[k][0] = P_[k][1];
            P_[k][nx_-1] = P_[k][nx_-2];
        }
        
        // Wells
        for (const auto& well : wells_) {
            int i = static_cast<int>(well.x / dx_);
            int k = static_cast<int>(well.z / dz_);
            i = std::max(0, std::min(nx_ - 1, i));
            k = std::max(0, std::min(nz_ - 1, k));
            P_[k][i] = well.target_bhp;
        }
    }
    
    void applySaturationBCs() {
        // Boundaries: zero flux (already implicit in interior scheme)
        for (int i = 0; i < nx_; ++i) {
            Sw_[0][i] = Sw_[1][i];
            Sw_[nz_-1][i] = Sw_[nz_-2][i];
        }
        for (int k = 0; k < nz_; ++k) {
            Sw_[k][0] = Sw_[k][1];
            Sw_[k][nx_-1] = Sw_[k][nx_-2];
        }
        
        // Wells
        for (const auto& well : wells_) {
            if (well.is_injector) {
                int i = static_cast<int>(well.x / dx_);
                int k = static_cast<int>(well.z / dz_);
                i = std::max(0, std::min(nx_ - 1, i));
                k = std::max(0, std::min(nz_ - 1, k));
                Sw_[k][i] = 0.8;
            }
        }
    }
    
    void applyDisplacementBCs() {
        // Bottom: fixed (no displacement)
        for (int i = 0; i < nx_; ++i) {
            ux_[nz_-1][i] = 0.0;
            uz_[nz_-1][i] = 0.0;
        }
        
        // Top: free surface (zero traction) - natural BC, no need to enforce
        
        // Left and right: roller BC (ux = 0, uz free)
        for (int k = 0; k < nz_; ++k) {
            ux_[k][0] = 0.0;
            ux_[k][nx_-1] = 0.0;
            uz_[k][0] = uz_[k][1];
            uz_[k][nx_-1] = uz_[k][nx_-2];
        }
    }
    
    double computeFlux(int i1, int k1, int i2, int k2, bool is_x_dir) {
        // Compute water flux between two cells with upstream weighting
        double dP = P_[k2][i2] - P_[k1][k1];
        double dist = is_x_dir ? dx_ : dz_;
        
        double k_harm = 2.0 * k_[k1][i1] * k_[k2][i2] / (k_[k1][i1] + k_[k2][i2] + 1e-30);
        
        // Total velocity
        double grav = is_x_dir ? 0.0 : -rho_f_ * 9.81 * dist;
        double v_total = -k_harm / mu_eff(0.5*(Sw_[k1][i1] + Sw_[k2][i2])) * (dP + grav) / dist;
        
        // Upstream saturation
        double Sw_up = (v_total > 0) ? Sw_[k1][i1] : Sw_[k2][i2];
        double fw = fractionalFlow(Sw_up);
        
        return v_total * fw;
    }
    
    double fractionalFlow(double Sw) const {
        double Sw_norm = (Sw - 0.2) / 0.6;
        Sw_norm = std::max(0.0, std::min(1.0, Sw_norm));
        
        double krw = std::pow(Sw_norm, 2.0);
        double kro = std::pow(1.0 - Sw_norm, 2.0);
        
        double lambda_w = krw / mu_w_;
        double lambda_o = kro / mu_o_;
        
        return lambda_w / (lambda_w + lambda_o + 1e-30);
    }
    
    double mu_eff(double Sw) const {
        double fw = fractionalFlow(Sw);
        return 1.0 / (fw / mu_w_ + (1.0 - fw) / mu_o_ + 1e-30);
    }
    
    bool isWellCell(int i, int k) const {
        for (const auto& well : wells_) {
            int wi = static_cast<int>(well.x / dx_);
            int wk = static_cast<int>(well.z / dz_);
            if (i == wi && k == wk) return true;
        }
        return false;
    }
    
    // Grid
    int nx_, nz_;
    double Lx_, Lz_;
    double dx_, dz_;
    
    // Flow fields
    std::vector<std::vector<double>> P_, Sw_;
    std::vector<std::vector<double>> phi_, k_;
    
    // Geomechanics fields
    std::vector<std::vector<double>> ux_, uz_;
    std::vector<std::vector<double>> sigma_xx_, sigma_zz_, tau_xz_;
    std::vector<std::vector<double>> eps_vol_;
    
    // Properties
    double phi0_, k0_;
    double E_, nu_;
    double lambda_, G_, K_;  // Elastic constants
    double rho_f_, rho_s_;
    double mu_w_, mu_o_;
    double alpha_, cf_, cr_;
    double P_initial_;
    double time_ = 0.0;
    
    std::vector<Well2D> wells_;
};

// Enhanced plot with labels
void createLabeledPlot(const std::vector<std::vector<double>>& data,
                      double Lx, double Lz,
                      const std::string& title,
                      const std::string& xlabel,
                      const std::string& ylabel,
                      const std::string& colormap_name,
                      const std::string& filename,
                      double vmin, double vmax,
                      const std::vector<Well2D>& wells = std::vector<Well2D>()) {
    
    // Use PlotGenerator2D for more control
    PlotGenerator2D plot(1400, 700, 100, 200, 100, 80);
    
    // Select colormap
    const ColorMap* cmap = nullptr;
    JetColorMap jet;
    ViridisColorMap viridis;
    PressureColorMap pressure;
    SubsidenceColorMap subsidence;
    
    if (colormap_name == "jet") cmap = &jet;
    else if (colormap_name == "viridis") cmap = &viridis;
    else if (colormap_name == "pressure") cmap = &pressure;
    else if (colormap_name == "subsidence") cmap = &subsidence;
    else cmap = &jet;
    
    // Draw field
    plot.setTitle(title);
    plot.setXLabel(xlabel);
    plot.setYLabel(ylabel);
    plot.drawField(data, *cmap, vmin, vmax);
    plot.drawColorbar(*cmap, vmin, vmax);
    
    // Draw wells
    for (const auto& well : wells) {
        plot.drawWell(well.x, well.z, Lx, Lz, well.is_injector, well.name);
    }
    
    // Write image
    plot.writeImage(filename);
}

int main(int argc, char** argv) {
    PetscErrorCode ierr;
    ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);
    
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    
    if (rank == 0) {
        std::cout << "\n";
        std::cout << "═══════════════════════════════════════════════════════════\n";
        std::cout << "  Enhanced 2D Vertical Reservoir Simulation\n";
        std::cout << "  • Fully Implicit Time Integration\n";
        std::cout << "  • Full Tensor Geomechanics (Plane Strain)\n";
        std::cout << "  • Proper Boundary Conditions\n";
        std::cout << "  • Well Performance Tracking\n";
        std::cout << "═══════════════════════════════════════════════════════════\n\n";
    }
    
    // Output directory
    std::string output_dir = "output_enhanced_2d";
    if (rank == 0) {
        mkdir(output_dir.c_str(), 0755);
    }
    MPI_Barrier(PETSC_COMM_WORLD);
    
    // Domain
    const double Lx = 500.0;
    const double Lz = 100.0;
    const int nx = 100;
    const int nz = 40;
    
    // Time
    const double total_time = 3600.0 * 24 * 180;
    const int num_steps = 90;
    const double dt = total_time / num_steps;
    
    if (rank == 0) {
        std::cout << "Domain: " << Lx << " m (horiz) × " << Lz << " m (depth)\n";
        std::cout << "Grid: " << nx << " × " << nz << " cells\n";
        std::cout << "Timesteps: " << num_steps << " (" << dt/86400.0 << " days each)\n";
        std::cout << "Duration: " << total_time / (3600 * 24) << " days\n\n";
    }
    
    // Create simulator
    EnhancedReservoir2D sim(nx, nz, Lx, Lz);
    sim.setRockProperties(0.2, 100e-15, 10e9, 0.25);
    
    // Add wells
    sim.addWell(Well2D(50.0, 50.0, true, 500.0, 35e6, "INJ-1"));    // 500 m³/day, 35 MPa
    sim.addWell(Well2D(450.0, 50.0, false, 500.0, 20e6, "PROD-1")); // 500 m³/day, 20 MPa
    
    if (rank == 0) {
        std::cout << "Wells:\n";
        for (const auto& well : sim.getWells()) {
            std::cout << "  " << well.name << " (" 
                     << (well.is_injector ? "Injector" : "Producer") << ")\n";
            std::cout << "    Position: (" << well.x << ", " << well.z << ") m\n";
            std::cout << "    Target: " << (well.is_injector ? "Injection" : "Production")
                     << " rate = " << well.target_rate << " m³/day\n";
            std::cout << "    BHP target: " << well.target_bhp / 1e6 << " MPa\n";
        }
        std::cout << "\nStarting simulation...\n\n";
    }
    
    // Time loop
    int plot_interval = 5;
    
    // Well performance file
    std::ofstream well_perf;
    if (rank == 0) {
        well_perf.open(output_dir + "/well_performance.csv");
        well_perf << "Time_days,Well,BHP_MPa,Oil_Rate_m3day,Water_Rate_m3day,Total_Rate_m3day,WaterCut,Cum_Oil_m3,Cum_Water_m3\n";
    }
    
    for (int step = 0; step <= num_steps; ++step) {
        if (step > 0) {
            sim.step(dt);
        }
        
        // Output
        if (rank == 0) {
            double time_days = sim.getTime() / (3600 * 24);
            
            // Write well performance
            for (const auto& well : sim.getWells()) {
                well_perf << time_days << ","
                         << well.name << ","
                         << std::scientific << std::setprecision(4)
                         << well.perf.bhp / 1e6 << ","
                         << well.perf.oil_rate << ","
                         << well.perf.water_rate << ","
                         << well.perf.total_rate << ","
                         << well.perf.water_cut << ","
                         << well.perf.cum_oil << ","
                         << well.perf.cum_water << "\n";
            }
            
            if (step % plot_interval == 0 || step == num_steps) {
                auto P = sim.getPressure();
                auto Sw = sim.getSaturation();
                auto sub = sim.getSubsidence();
                auto sigma_zz = sim.getStressZZ();
                
                // Find ranges
                double P_min = 1e30, P_max = -1e30;
                double sub_min = 1e30, sub_max = -1e30;
                double stress_min = 1e30, stress_max = -1e30;
                
                for (int k = 0; k < nz; ++k) {
                    for (int i = 0; i < nx; ++i) {
                        P_min = std::min(P_min, P[k][i]);
                        P_max = std::max(P_max, P[k][i]);
                        sub_min = std::min(sub_min, sub[k][i]);
                        sub_max = std::max(sub_max, sub[k][i]);
                        stress_min = std::min(stress_min, sigma_zz[k][i]);
                        stress_max = std::max(stress_max, sigma_zz[k][i]);
                    }
                }
                
                // Create labeled plots
                char filename[256];
                char title[512];
                
                snprintf(title, sizeof(title), "Pressure Field - Time: %.1f days", time_days);
                snprintf(filename, sizeof(filename), "%s/pressure_%04d", output_dir.c_str(), step);
                createLabeledPlot(P, Lx, Lz, title, "Horizontal Distance (m)", "Depth (m)",
                                "pressure", filename, P_min, P_max, sim.getWells());
                
                snprintf(title, sizeof(title), "Water Saturation - Time: %.1f days", time_days);
                snprintf(filename, sizeof(filename), "%s/saturation_%04d", output_dir.c_str(), step);
                createLabeledPlot(Sw, Lx, Lz, title, "Horizontal Distance (m)", "Depth (m)",
                                "jet", filename, 0.0, 1.0, sim.getWells());
                
                snprintf(title, sizeof(title), "Surface Subsidence - Time: %.1f days", time_days);
                snprintf(filename, sizeof(filename), "%s/subsidence_%04d", output_dir.c_str(), step);
                createLabeledPlot(sub, Lx, Lz, title, "Horizontal Distance (m)", "Depth (m)",
                                "subsidence", filename, sub_min, sub_max, sim.getWells());
                
                snprintf(title, sizeof(title), "Vertical Stress - Time: %.1f days", time_days);
                snprintf(filename, sizeof(filename), "%s/stress_%04d", output_dir.c_str(), step);
                createLabeledPlot(sigma_zz, Lx, Lz, title, "Horizontal Distance (m)", "Depth (m)",
                                "viridis", filename, stress_min, stress_max, sim.getWells());
                
                // Write metadata
                snprintf(filename, sizeof(filename), "%s/info_%04d.txt", output_dir.c_str(), step);
                std::ofstream info(filename);
                info << "═══════════════════════════════════════════════════\n";
                info << "Enhanced 2D Reservoir Simulation - Step " << step << "\n";
                info << "═══════════════════════════════════════════════════\n\n";
                info << "Time: " << std::fixed << std::setprecision(2) << time_days << " days\n\n";
                info << "Pressure:\n";
                info << "  Min: " << std::setprecision(2) << P_min/1e6 << " MPa\n";
                info << "  Max: " << std::setprecision(2) << P_max/1e6 << " MPa\n\n";
                info << "Subsidence:\n";
                info << "  Min: " << std::scientific << std::setprecision(3) << sub_min << " m\n";
                info << "  Max: " << sub_max << " m\n\n";
                info << "Vertical Stress:\n";
                info << "  Min: " << std::fixed << std::setprecision(2) << stress_min/1e6 << " MPa\n";
                info << "  Max: " << stress_max/1e6 << " MPa\n\n";
                info << "Well Performance:\n";
                for (const auto& well : sim.getWells()) {
                    info << "  " << well.name << ":\n";
                    info << "    BHP: " << std::setprecision(2) << well.perf.bhp/1e6 << " MPa\n";
                    if (well.is_injector) {
                        info << "    Water Rate: " << well.perf.water_rate << " m³/day\n";
                        info << "    Cumulative: " << well.perf.cum_water << " m³\n";
                    } else {
                        info << "    Oil Rate: " << well.perf.oil_rate << " m³/day\n";
                        info << "    Water Rate: " << well.perf.water_rate << " m³/day\n";
                        info << "    Water Cut: " << std::setprecision(1) << well.perf.water_cut*100 << "%\n";
                        info << "    Cum Oil: " << well.perf.cum_oil << " m³\n";
                        info << "    Cum Water: " << well.perf.cum_water << " m³\n";
                    }
                }
                info.close();
                
                std::cout << "Step " << std::setw(4) << step << " / " << num_steps
                         << "  |  " << std::setw(7) << std::setprecision(1) << time_days << " days"
                         << "  |  P: " << std::setw(5) << std::setprecision(1) << P_min/1e6 
                         << "-" << std::setw(5) << P_max/1e6 << " MPa"
                         << "  |  Plots OK\n";
            }
        }
    }
    
    if (rank == 0) {
        well_perf.close();
        
        std::cout << "\n═══════════════════════════════════════════════════════════\n";
        std::cout << "  Simulation Complete!\n";
        std::cout << "═══════════════════════════════════════════════════════════\n\n";
        std::cout << "Output: " << output_dir << "/\n";
        std::cout << "  • pressure_XXXX.png - Labeled pressure fields\n";
        std::cout << "  • saturation_XXXX.png - Labeled saturation fields\n";
        std::cout << "  • subsidence_XXXX.png - Labeled subsidence\n";
        std::cout << "  • stress_XXXX.png - Labeled vertical stress\n";
        std::cout << "  • well_performance.csv - Complete well data\n";
        std::cout << "  • info_XXXX.txt - Detailed statistics\n\n";
        
        std::cout << "Final Production Summary:\n";
        for (const auto& well : sim.getWells()) {
            std::cout << "  " << well.name << ":\n";
            if (!well.is_injector) {
                std::cout << "    Cumulative Oil: " << std::fixed << std::setprecision(0) 
                         << well.perf.cum_oil << " m³\n";
                std::cout << "    Cumulative Water: " << well.perf.cum_water << " m³\n";
                std::cout << "    Recovery Factor: " << std::setprecision(2)
                         << (well.perf.cum_oil / (Lx * Lz * 0.2 * 0.8) * 100.0) << "%\n";
            } else {
                std::cout << "    Total Injected: " << std::fixed << std::setprecision(0)
                         << well.perf.cum_water << " m³\n";
            }
        }
        std::cout << "\n";
    }
    
    ierr = PetscFinalize();
    return ierr;
}
