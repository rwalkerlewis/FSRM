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

using namespace ResSim;

/**
 * 2D Vertical Section Reservoir Simulation
 * 
 * Demonstrates coupled fluid flow + geomechanics in a vertical slice
 * - Numerically stable implicit solver
 * - 1 injector (left) + 1 producer (right)
 * - Vertical cross-section showing depth effects
 * - Uses generic visualization from Visualization class
 */

struct Well2D {
    double x, z;  // Position (horizontal, vertical)
    bool is_injector;
    double rate;  // m³/s
    std::string name;
    
    Well2D(double x_pos, double z_pos, bool inj, double r, const std::string& n)
        : x(x_pos), z(z_pos), is_injector(inj), rate(r), name(n) {}
};

class StableReservoir2D {
public:
    StableReservoir2D(int nx, int nz, double Lx, double Lz)
        : nx_(nx), nz_(nz), Lx_(Lx), Lz_(Lz),
          dx_(Lx / nx), dz_(Lz / nz) {
        
        // Initialize fields
        P_.resize(nz, std::vector<double>(nx, 0.0));
        Sw_.resize(nz, std::vector<double>(nx, 0.0));
        ux_.resize(nz, std::vector<double>(nx, 0.0));
        uz_.resize(nz, std::vector<double>(nx, 0.0));
        subsidence_.resize(nz, std::vector<double>(nx, 0.0));
        phi_.resize(nz, std::vector<double>(nx, 0.0));
        k_.resize(nz, std::vector<double>(nx, 0.0));
        
        // Default properties
        phi0_ = 0.2;
        k0_ = 100e-15;  // 100 mD
        E_ = 10e9;      // 10 GPa
        nu_ = 0.25;
        rho_f_ = 1000.0;
        mu_ = 1e-3;
        alpha_ = 0.7;
        cf_ = 1e-9;     // Fluid compressibility
        cr_ = 1e-9;     // Rock compressibility
        
        P_initial_ = 30e6;  // 30 MPa
        
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
    }
    
    void initialize() {
        for (int k = 0; k < nz_; ++k) {
            for (int i = 0; i < nx_; ++i) {
                double depth = k * dz_;
                
                // Hydrostatic pressure
                P_[k][i] = P_initial_ + rho_f_ * 9.81 * depth;
                
                // Initial water saturation
                Sw_[k][i] = 0.2;  // Connate
                
                // Rock properties
                phi_[k][i] = phi0_;
                k_[k][i] = k0_;
            }
        }
    }
    
    void step(double dt) {
        // Implicit pressure solve with better stability
        solvePressureImplicit(dt);
        
        // Explicit saturation transport (with CFL limit)
        double cfl_dt = computeCFLTimestep();
        if (dt > cfl_dt) {
            int n_sub = static_cast<int>(std::ceil(dt / cfl_dt));
            double sub_dt = dt / n_sub;
            for (int sub = 0; sub < n_sub; ++sub) {
                solveSaturationExplicit(sub_dt);
            }
        } else {
            solveSaturationExplicit(dt);
        }
        
        // Geomechanics
        solveGeomechanics();
        
        // Update properties
        updatePoroPermFromGeomechanics();
        
        time_ += dt;
    }
    
    const std::vector<std::vector<double>>& getPressure() const { return P_; }
    const std::vector<std::vector<double>>& getSaturation() const { return Sw_; }
    const std::vector<std::vector<double>>& getSubsidence() const { return subsidence_; }
    const std::vector<std::vector<double>>& getPermeability() const { return k_; }
    const std::vector<Well2D>& getWells() const { return wells_; }
    
    double getTime() const { return time_; }
    double getLx() const { return Lx_; }
    double getLz() const { return Lz_; }
    int getNx() const { return nx_; }
    int getNz() const { return nz_; }
    
private:
    void solvePressureImplicit(double dt) {
        // Use Gauss-Seidel iteration for implicit solve
        auto P_old = P_;
        
        const int max_iter = 50;
        const double tol = 1e-6;
        
        for (int iter = 0; iter < max_iter; ++iter) {
            double max_change = 0.0;
            
            for (int k = 1; k < nz_ - 1; ++k) {
                for (int i = 1; i < nx_ - 1; ++i) {
                    double phi = phi_[k][i];
                    double perm = k_[k][i];
                    double ct = cf_ + cr_;
                    
                    // Transmissibilities
                    double Tx_plus = 2.0 * perm * k_[k][i+1] / (perm + k_[k][i+1] + 1e-30) / (mu_ * dx_ * dx_);
                    double Tx_minus = 2.0 * perm * k_[k][i-1] / (perm + k_[k][i-1] + 1e-30) / (mu_ * dx_ * dx_);
                    double Tz_plus = 2.0 * perm * k_[k+1][i] / (perm + k_[k+1][i] + 1e-30) / (mu_ * dz_ * dz_);
                    double Tz_minus = 2.0 * perm * k_[k-1][i] / (perm + k_[k-1][i] + 1e-30) / (mu_ * dz_ * dz_);
                    
                    // Gravity term in z-direction
                    double gravity_flux = (Tz_plus - Tz_minus) * rho_f_ * 9.81 * dz_ / 2.0;
                    
                    // Accumulation
                    double accum_coeff = phi * ct / dt;
                    
                    // Linear system coefficients
                    double a_c = accum_coeff + Tx_plus + Tx_minus + Tz_plus + Tz_minus;
                    double rhs = accum_coeff * P_old[k][i] + gravity_flux;
                    rhs += Tx_plus * P_[k][i+1] + Tx_minus * P_[k][i-1];
                    rhs += Tz_plus * P_[k+1][i] + Tz_minus * P_[k-1][i];
                    
                    double P_new = rhs / (a_c + 1e-30);
                    
                    max_change = std::max(max_change, std::abs(P_new - P_[k][i]));
                    P_[k][i] = P_new;
                }
            }
            
            // Apply well conditions
            applyWellBCs();
            
            // Apply boundary conditions
            applyPressureBCs();
            
            if (max_change < tol * P_initial_) break;
        }
    }
    
    void solveSaturationExplicit(double dt) {
        auto Sw_new = Sw_;
        
        for (int k = 1; k < nz_ - 1; ++k) {
            for (int i = 1; i < nx_ - 1; ++i) {
                double phi = phi_[k][i];
                double perm = k_[k][i];
                
                // Compute phase velocities (Darcy)
                double dP_dx = (P_[k][i+1] - P_[k][i-1]) / (2.0 * dx_);
                double dP_dz = (P_[k+1][i] - P_[k-1][i]) / (2.0 * dz_);
                
                // Water fractional flow
                double fw = fractionalFlow(Sw_[k][i]);
                double fw_left = fractionalFlow(Sw_[k][i-1]);
                double fw_right = fractionalFlow(Sw_[k][i+1]);
                
                // Total velocity (simplified)
                double vx = -(perm / mu_) * dP_dx;
                double vz = -(perm / mu_) * (dP_dz - rho_f_ * 9.81);
                
                // Upwind saturation advection in x
                double flux_x = 0.0;
                if (vx > 0) {
                    flux_x = vx * fw_left;
                } else {
                    flux_x = vx * fw_right;
                }
                
                double dSw_dx = (flux_x) / dx_;
                
                // Update saturation
                Sw_new[k][i] = Sw_[k][i] - (dt / phi) * dSw_dx;
                
                // Clamp to physical bounds
                Sw_new[k][i] = std::max(0.2, std::min(0.8, Sw_new[k][i]));
            }
        }
        
        Sw_ = Sw_new;
        
        // Apply saturation BCs from wells
        applySaturationBCs();
    }
    
    double computeCFLTimestep() const {
        double max_velocity = 0.0;
        
        for (int k = 1; k < nz_ - 1; ++k) {
            for (int i = 1; i < nx_ - 1; ++i) {
                double dP_dx = std::abs(P_[k][i+1] - P_[k][i-1]) / (2.0 * dx_);
                double vx = (k_[k][i] / mu_) * dP_dx;
                max_velocity = std::max(max_velocity, std::abs(vx));
            }
        }
        
        double cfl_dt = 0.5 * dx_ * phi0_ / (max_velocity + 1e-20);
        return std::min(cfl_dt, 1e6);  // Cap at reasonable value
    }
    
    void solveGeomechanics() {
        // Simple uniaxial compaction
        for (int k = 0; k < nz_; ++k) {
            for (int i = 0; i < nx_; ++i) {
                double dP = P_[k][i] - (P_initial_ + rho_f_ * 9.81 * k * dz_);
                
                // Bulk modulus
                double K = E_ / (3.0 * (1.0 - 2.0 * nu_));
                
                // Vertical strain
                double eps_v = alpha_ * dP / K;
                
                // Subsidence (accumulate from bottom)
                subsidence_[k][i] = eps_v * dz_;
                
                // Small horizontal displacement
                double eps_h = -nu_ / (1.0 - nu_) * eps_v;
                ux_[k][i] = eps_h * (i * dx_ - Lx_ / 2.0);
                uz_[k][i] = -subsidence_[k][i];
            }
        }
    }
    
    void updatePoroPermFromGeomechanics() {
        for (int k = 0; k < nz_; ++k) {
            for (int i = 0; i < nx_; ++i) {
                double dP = P_[k][i] - (P_initial_ + rho_f_ * 9.81 * k * dz_);
                
                // Porosity from compressibility
                phi_[k][i] = phi0_ * (1.0 + cr_ * dP);
                phi_[k][i] = std::max(0.01, std::min(0.4, phi_[k][i]));
                
                // Kozeny-Carman
                double phi_ratio = phi_[k][i] / phi0_;
                k_[k][i] = k0_ * std::pow(phi_ratio, 3.0) * 
                          std::pow((1.0 - phi0_) / (1.0 - phi_[k][i] + 1e-6), 2.0);
            }
        }
    }
    
    void applyWellBCs() {
        for (const auto& well : wells_) {
            int i = static_cast<int>(well.x / dx_);
            int k = static_cast<int>(well.z / dz_);
            
            i = std::max(0, std::min(nx_ - 1, i));
            k = std::max(0, std::min(nz_ - 1, k));
            
            if (well.is_injector) {
                P_[k][i] = P_initial_ + 5e6;  // 5 MPa overpressure
            } else {
                P_[k][i] = P_initial_ - 10e6;  // 10 MPa drawdown
            }
        }
    }
    
    void applySaturationBCs() {
        for (const auto& well : wells_) {
            int i = static_cast<int>(well.x / dx_);
            int k = static_cast<int>(well.z / dz_);
            
            i = std::max(0, std::min(nx_ - 1, i));
            k = std::max(0, std::min(nz_ - 1, k));
            
            if (well.is_injector) {
                Sw_[k][i] = 0.8;  // Inject water
            }
        }
    }
    
    void applyPressureBCs() {
        // No-flow on top and bottom
        for (int i = 0; i < nx_; ++i) {
            P_[0][i] = P_[1][i];
            P_[nz_-1][i] = P_[nz_-2][i];
        }
        
        // No-flow on left and right (except wells)
        for (int k = 0; k < nz_; ++k) {
            P_[k][0] = P_[k][1];
            P_[k][nx_-1] = P_[k][nx_-2];
        }
    }
    
    double fractionalFlow(double Sw) const {
        double Sw_norm = (Sw - 0.2) / 0.6;
        Sw_norm = std::max(0.0, std::min(1.0, Sw_norm));
        
        double krw = std::pow(Sw_norm, 2.0);
        double kro = std::pow(1.0 - Sw_norm, 2.0);
        
        double mu_w = mu_;
        double mu_o = 2.0 * mu_;
        
        double lambda_w = krw / mu_w;
        double lambda_o = kro / mu_o;
        
        if (lambda_w + lambda_o < 1e-12) return 0.0;
        return lambda_w / (lambda_w + lambda_o);
    }
    
    int nx_, nz_;
    double Lx_, Lz_;
    double dx_, dz_;
    
    std::vector<std::vector<double>> P_, Sw_;
    std::vector<std::vector<double>> ux_, uz_;
    std::vector<std::vector<double>> subsidence_;
    std::vector<std::vector<double>> phi_, k_;
    
    double phi0_, k0_;
    double E_, nu_;
    double rho_f_, mu_;
    double alpha_, cf_, cr_;
    double P_initial_;
    double time_ = 0.0;
    
    std::vector<Well2D> wells_;
};

int main(int argc, char** argv) {
    PetscErrorCode ierr;
    ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);
    
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    
    if (rank == 0) {
        std::cout << "\n";
        std::cout << "═══════════════════════════════════════════════════\n";
        std::cout << "  2D Vertical Reservoir Section\n";
        std::cout << "  Flow + Wells + Geomechanics\n";
        std::cout << "═══════════════════════════════════════════════════\n\n";
    }
    
    // Output directory
    std::string output_dir = "output_reservoir_2d_vertical";
    if (rank == 0) {
        mkdir(output_dir.c_str(), 0755);
    }
    MPI_Barrier(PETSC_COMM_WORLD);
    
    // Domain (vertical cross-section)
    const double Lx = 500.0;  // Horizontal extent (m)
    const double Lz = 100.0;  // Vertical extent (m)
    const int nx = 150;
    const int nz = 40;
    
    // Time
    const double total_time = 3600.0 * 24 * 180;  // 180 days
    const int num_steps = 100;
    const double dt = total_time / num_steps;
    
    if (rank == 0) {
        std::cout << "Domain: " << Lx << " m (horiz) × " << Lz << " m (vert)\n";
        std::cout << "Grid: " << nx << " × " << nz << " cells\n";
        std::cout << "Timesteps: " << num_steps << "\n";
        std::cout << "Duration: " << total_time / (3600 * 24) << " days\n\n";
    }
    
    // Create simulator
    StableReservoir2D sim(nx, nz, Lx, Lz);
    
    // Set properties
    sim.setRockProperties(0.2, 100e-15, 10e9, 0.25);
    
    // Add wells (injector left, producer right)
    sim.addWell(Well2D(50.0, 50.0, true, 0.01, "INJ-1"));
    sim.addWell(Well2D(450.0, 50.0, false, -0.01, "PROD-1"));
    
    if (rank == 0) {
        std::cout << "Wells:\n";
        for (const auto& well : sim.getWells()) {
            std::cout << "  " << well.name << " (" 
                     << (well.is_injector ? "Injector" : "Producer") << ")\n";
            std::cout << "    Position: (" << well.x << ", " << well.z << ") m\n";
        }
        std::cout << "\nStarting simulation...\n\n";
    }
    
    // Visualization
    Visualization viz;
    viz.setOutputDirectory(output_dir);
    
    // Time loop
    int plot_interval = 5;
    
    for (int step = 0; step <= num_steps; ++step) {
        if (step > 0) {
            sim.step(dt);
        }
        
        // Output plots
        if (rank == 0 && (step % plot_interval == 0 || step == num_steps)) {
            double time_days = sim.getTime() / (3600 * 24);
            
            // Get field data
            auto P = sim.getPressure();
            auto Sw = sim.getSaturation();
            auto sub = sim.getSubsidence();
            auto k = sim.getPermeability();
            
            // Find value ranges
            double P_min = 1e30, P_max = -1e30;
            double sub_min = 1e30, sub_max = -1e30;
            double k_min = 1e30, k_max = -1e30;
            
            for (int kk = 0; kk < nz; ++kk) {
                for (int i = 0; i < nx; ++i) {
                    P_min = std::min(P_min, P[kk][i]);
                    P_max = std::max(P_max, P[kk][i]);
                    sub_min = std::min(sub_min, sub[kk][i]);
                    sub_max = std::max(sub_max, sub[kk][i]);
                    k_min = std::min(k_min, k[kk][i]);
                    k_max = std::max(k_max, k[kk][i]);
                }
            }
            
            // Plot pressure
            char filename[256];
            snprintf(filename, sizeof(filename), "pressure_%04d", step);
            viz.plot2DField(P, Lx, Lz, "Pressure Field", "pressure",
                          filename, P_min, P_max);
            
            // Plot saturation
            snprintf(filename, sizeof(filename), "saturation_%04d", step);
            viz.plot2DField(Sw, Lx, Lz, "Water Saturation", "jet",
                          filename, 0.0, 1.0);
            
            // Plot subsidence
            snprintf(filename, sizeof(filename), "subsidence_%04d", step);
            viz.plot2DField(sub, Lx, Lz, "Subsidence", "subsidence",
                          filename, sub_min, sub_max);
            
            // Write metadata
            snprintf(filename, sizeof(filename), "%s/info_%04d.txt",
                    output_dir.c_str(), step);
            std::ofstream info(filename);
            info << "Step: " << step << "\n";
            info << "Time: " << std::fixed << std::setprecision(2) 
                 << time_days << " days\n";
            info << "Pressure: " << std::setprecision(2) 
                 << P_min/1e6 << " - " << P_max/1e6 << " MPa\n";
            info << "Subsidence: " << std::scientific << std::setprecision(3)
                 << sub_min << " - " << sub_max << " m\n";
            info << "Permeability: " << std::setprecision(1)
                 << k_min*1e15 << " - " << k_max*1e15 << " mD\n";
            info.close();
            
            std::cout << "Step " << std::setw(4) << step << " / " << num_steps
                     << "  |  Time: " << std::setw(7) << std::setprecision(1) 
                     << time_days << " days  |  Plots generated\n";
        }
    }
    
    if (rank == 0) {
        std::cout << "\n═══════════════════════════════════════════════════\n";
        std::cout << "  Simulation Complete!\n";
        std::cout << "═══════════════════════════════════════════════════\n\n";
        std::cout << "Output: " << output_dir << "/\n";
        std::cout << "  • pressure_XXXX.png\n";
        std::cout << "  • saturation_XXXX.png\n";
        std::cout << "  • subsidence_XXXX.png\n";
        std::cout << "  • info_XXXX.txt\n\n";
    }
    
    ierr = PetscFinalize();
    return ierr;
}
