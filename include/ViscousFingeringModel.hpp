/**
 * @file ViscousFingeringModel.hpp
 * @brief Viscous fingering physics model
 * 
 * Implements unstable displacement when low viscosity fluid
 * displaces high viscosity fluid (unfavorable mobility ratio).
 * 
 * Features:
 * - Saffman-Taylor instability
 * - Perturbed interface tracking
 * - Mixing zone modeling
 * - Mobility ratio effects
 * - Finger width prediction
 */

#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include <random>

namespace FSRM {

/**
 * @class ViscousFingeringModel
 * @brief Models viscous instabilities in two-phase displacement
 */
class ViscousFingeringModel {
public:
    /**
     * Constructor
     * @param nx Number of cells in x-direction
     * @param ny Number of cells in y-direction
     * @param Lx Domain length in x (m)
     * @param Ly Domain length in y (m)
     */
    ViscousFingeringModel(int nx, int ny, double Lx, double Ly)
        : nx_(nx), ny_(ny), Lx_(Lx), Ly_(Ly),
          dx_(Lx / nx), dy_(Ly / ny) {
        
        // Initialize fields
        saturation_.resize(nx * ny, 0.0);
        velocity_x_.resize(nx * ny, 0.0);
        velocity_y_.resize(nx * ny, 0.0);
        concentration_.resize(nx * ny, 0.0);
        
        // Initialize random number generator for perturbations
        rng_.seed(12345);
    }
    
    /**
     * Set fluid properties
     * @param mu_inj Injected fluid viscosity (Pa·s)
     * @param mu_disp Displaced fluid viscosity (Pa·s)
     * @param rho_inj Injected fluid density (kg/m³)
     * @param rho_disp Displaced fluid density (kg/m³)
     */
    void setFluidProperties(double mu_inj, double mu_disp,
                           double rho_inj, double rho_disp) {
        mu_injected_ = mu_inj;
        mu_displaced_ = mu_disp;
        rho_injected_ = rho_inj;
        rho_displaced_ = rho_disp;
        
        // Mobility ratio M = λ_inj / λ_disp
        // For equal rel perms: M = μ_disp / μ_inj
        M_ = mu_displaced_ / mu_injected_;
        
        // Stability: M > 1 is stable (viscous displacing less viscous)
        //           M < 1 is unstable (fingers form)
        is_unstable_ = (M_ < 1.0);
    }
    
    /**
     * Set reservoir properties
     * @param k Permeability (m²)
     * @param phi Porosity
     * @param S_wc Connate water saturation
     * @param S_or Residual oil saturation
     */
    void setReservoirProperties(double k, double phi, double S_wc, double S_or) {
        permeability_ = k;
        porosity_ = phi;
        S_wc_ = S_wc;
        S_or_ = S_or;
    }
    
    /**
     * Initialize saturation field with perturbation
     * @param amplitude Perturbation amplitude (0-1)
     * @param wavelength Perturbation wavelength (m)
     */
    void initializeWithPerturbation(double amplitude, double wavelength) {
        std::uniform_real_distribution<double> dist(-amplitude, amplitude);
        
        for (int j = 0; j < ny_; ++j) {
            for (int i = 0; i < nx_; ++i) {
                int idx = j * nx_ + i;
                
                // Base saturation (interface at x = 0.1*Lx)
                double x = i * dx_;
                double base_sat = (x < 0.1 * Lx_) ? 1.0 : 0.0;
                
                // Add random perturbation
                double perturb = dist(rng_);
                
                // Add sinusoidal perturbation in y-direction
                double y = j * dy_;
                double k_wave = 2.0 * M_PI / wavelength;
                double sin_perturb = amplitude * std::sin(k_wave * y);
                
                saturation_[idx] = std::max(0.0, std::min(1.0, 
                    base_sat + perturb + sin_perturb));
            }
        }
    }
    
    /**
     * Compute fractional flow function
     * @param S Saturation
     * @return Fractional flow fw
     */
    double fractionalFlow(double S) const {
        // Normalized saturation
        double S_e = (S - S_wc_) / (1.0 - S_wc_ - S_or_);
        S_e = std::max(0.0, std::min(1.0, S_e));
        
        // Corey rel perms
        double k_rw = S_e * S_e;
        double k_ro = (1.0 - S_e) * (1.0 - S_e);
        
        // Mobility ratio in terms of rel perms
        double lambda_w = k_rw / mu_injected_;
        double lambda_o = k_ro / mu_displaced_;
        double lambda_t = lambda_w + lambda_o;
        
        if (lambda_t < 1.0e-10) return 0.0;
        
        return lambda_w / lambda_t;
    }
    
    /**
     * Compute fractional flow derivative
     * @param S Saturation
     * @return Derivative df/dS
     */
    double fractionalFlowDerivative(double S) const {
        double dS = 1.0e-6;
        double f1 = fractionalFlow(S + dS);
        double f0 = fractionalFlow(S - dS);
        return (f1 - f0) / (2.0 * dS);
    }
    
    /**
     * Time step saturation equation with fingering
     * @param dt Time step (s)
     * @param q_total Total flow rate (m³/s)
     */
    void timeStep(double dt, double q_total) {
        // Darcy velocity
        double A = Ly_ * 1.0;  // Cross-sectional area (assuming unit depth)
        double u_total = q_total / A;
        
        std::vector<double> saturation_new = saturation_;
        
        // Upwind finite volume for saturation transport
        for (int j = 1; j < ny_ - 1; ++j) {
            for (int i = 1; i < nx_ - 1; ++i) {
                int idx = j * nx_ + i;
                
                // Compute fractional flow
                double f_w = fractionalFlow(saturation_[idx]);
                
                // Velocity (including instability amplification)
                double v_x = u_total * f_w / porosity_;
                
                // Add dispersion for unstable case
                if (is_unstable_) {
                    // Add cross-flow due to fingering
                    double finger_velocity = 0.1 * v_x * (1.0 - M_);
                    
                    // Check neighbors for gradients
                    double dS_dx = (saturation_[idx + 1] - saturation_[idx - 1]) / (2.0 * dx_);
                    double dS_dy = (saturation_[idx + nx_] - saturation_[idx - nx_]) / (2.0 * dy_);
                    
                    // Add finger-induced flux
                    v_x += finger_velocity * dS_dx;
                    double v_y = finger_velocity * dS_dy;
                    
                    velocity_y_[idx] = v_y;
                }
                
                velocity_x_[idx] = v_x;
                
                // Upwind flux
                double flux_x_left = (v_x > 0) ? 
                    v_x * saturation_[idx - 1] : v_x * saturation_[idx];
                double flux_x_right = (v_x > 0) ? 
                    v_x * saturation_[idx] : v_x * saturation_[idx + 1];
                
                // Update saturation
                double dS_dt = -(flux_x_right - flux_x_left) / dx_;
                saturation_new[idx] = saturation_[idx] + dt * dS_dt;
                
                // Bound saturation
                saturation_new[idx] = std::max(S_wc_, std::min(1.0 - S_or_, saturation_new[idx]));
            }
        }
        
        saturation_ = saturation_new;
        time_ += dt;
    }
    
    /**
     * Get finger penetration depth
     * @return Penetration depth (m)
     */
    double getFingerPenetration() const {
        double max_penetration = 0.0;
        
        for (int j = 0; j < ny_; ++j) {
            for (int i = 0; i < nx_; ++i) {
                int idx = j * nx_ + i;
                if (saturation_[idx] > 0.5) {
                    double x = i * dx_;
                    max_penetration = std::max(max_penetration, x);
                }
            }
        }
        
        return max_penetration;
    }
    
    /**
     * Get mixing zone width (where saturation varies)
     * @return Mixing width (m)
     */
    double getMixingZoneWidth() const {
        std::vector<double> x_low, x_high;
        
        for (int j = 0; j < ny_; ++j) {
            double x_02 = -1.0, x_08 = -1.0;
            for (int i = 0; i < nx_; ++i) {
                int idx = j * nx_ + i;
                double x = i * dx_;
                
                if (saturation_[idx] >= 0.2 && x_02 < 0) x_02 = x;
                if (saturation_[idx] >= 0.8 && x_08 < 0) x_08 = x;
            }
            
            if (x_02 >= 0) x_low.push_back(x_02);
            if (x_08 >= 0) x_high.push_back(x_08);
        }
        
        if (x_low.empty() || x_high.empty()) return 0.0;
        
        double avg_low = 0.0, avg_high = 0.0;
        for (double x : x_low) avg_low += x;
        for (double x : x_high) avg_high += x;
        avg_low /= x_low.size();
        avg_high /= x_high.size();
        
        return avg_high - avg_low;
    }
    
    /**
     * Compute sweep efficiency
     * @return Sweep efficiency (0-1)
     */
    double getSweepEfficiency() const {
        double total_volume = 0.0;
        double swept_volume = 0.0;
        
        for (int j = 0; j < ny_; ++j) {
            for (int i = 0; i < nx_; ++i) {
                int idx = j * nx_ + i;
                total_volume += 1.0;
                if (saturation_[idx] > 0.5) swept_volume += 1.0;
            }
        }
        
        return swept_volume / total_volume;
    }
    
    // Getters
    const std::vector<double>& getSaturation() const { return saturation_; }
    const std::vector<double>& getVelocityX() const { return velocity_x_; }
    const std::vector<double>& getVelocityY() const { return velocity_y_; }
    double getMobilityRatio() const { return M_; }
    bool isUnstable() const { return is_unstable_; }
    double getCurrentTime() const { return time_; }
    int getNx() const { return nx_; }
    int getNy() const { return ny_; }
    double getDx() const { return dx_; }
    double getDy() const { return dy_; }
    
private:
    // Grid
    int nx_, ny_;
    double Lx_, Ly_;
    double dx_, dy_;
    
    // Fields
    std::vector<double> saturation_;
    std::vector<double> velocity_x_;
    std::vector<double> velocity_y_;
    std::vector<double> concentration_;
    
    // Fluid properties
    double mu_injected_ = 0.001;
    double mu_displaced_ = 0.01;
    double rho_injected_ = 1000.0;
    double rho_displaced_ = 850.0;
    double M_ = 1.0;  // Mobility ratio
    bool is_unstable_ = false;
    
    // Rock properties
    double permeability_ = 100.0e-15;
    double porosity_ = 0.25;
    double S_wc_ = 0.2;
    double S_or_ = 0.2;
    
    // Time
    double time_ = 0.0;
    
    // Random number generator for perturbations
    std::mt19937 rng_;
};

} // namespace FSRM
