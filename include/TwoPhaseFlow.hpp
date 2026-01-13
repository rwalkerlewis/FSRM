#ifndef TWO_PHASE_FLOW_HPP
#define TWO_PHASE_FLOW_HPP

#include "FSRM.hpp"
#include <vector>
#include <memory>
#include <functional>

namespace FSRM {

/**
 * @brief Relative permeability model
 */
struct RelPermModel {
    double Swc = 0.2;       ///< Connate water saturation
    double Sor = 0.2;       ///< Residual oil saturation
    double krw_max = 1.0;   ///< Max water relative perm
    double kro_max = 1.0;   ///< Max oil relative perm
    double nw = 2.0;        ///< Water Corey exponent
    double no = 2.0;        ///< Oil Corey exponent
    
    double normalizedSaturation(double Sw) const {
        return (Sw - Swc) / (1.0 - Swc - Sor);
    }
    
    double krw(double Sw) const {
        double Sw_norm = normalizedSaturation(Sw);
        Sw_norm = std::max(0.0, std::min(1.0, Sw_norm));
        return krw_max * std::pow(Sw_norm, nw);
    }
    
    double kro(double Sw) const {
        double Sw_norm = normalizedSaturation(Sw);
        Sw_norm = std::max(0.0, std::min(1.0, Sw_norm));
        return kro_max * std::pow(1.0 - Sw_norm, no);
    }
};

/**
 * @brief Fractional flow calculation
 */
class FractionalFlow {
public:
    FractionalFlow() = default;
    
    void setRelPermModel(const RelPermModel& model) { relperm_ = model; }
    void setViscosities(double mu_w, double mu_o) { mu_w_ = mu_w; mu_o_ = mu_o; }
    
    /**
     * @brief Calculate water fractional flow
     * @param Sw Water saturation
     * @return Fractional flow of water (0-1)
     */
    double fw(double Sw) const {
        double krw = relperm_.krw(Sw);
        double kro = relperm_.kro(Sw);
        
        double lambda_w = krw / mu_w_;
        double lambda_o = kro / mu_o_;
        
        double total = lambda_w + lambda_o;
        if (total < 1e-12) return 0.0;
        
        return lambda_w / total;
    }
    
    /**
     * @brief Calculate derivative of fractional flow
     * @param Sw Water saturation
     * @return dfw/dSw
     */
    double dfw_dSw(double Sw) const {
        double eps = 1e-6;
        double Sw_plus = std::min(1.0 - relperm_.Sor, Sw + eps);
        double Sw_minus = std::max(relperm_.Swc, Sw - eps);
        return (fw(Sw_plus) - fw(Sw_minus)) / (Sw_plus - Sw_minus);
    }
    
    /**
     * @brief Calculate frontal advance velocity (Buckley-Leverett)
     * @param Sw Water saturation
     * @param u_total Total Darcy velocity
     * @param phi Porosity
     * @return Saturation front velocity
     */
    double frontVelocity(double Sw, double u_total, double phi) const {
        return u_total * dfw_dSw(Sw) / phi;
    }
    
private:
    RelPermModel relperm_;
    double mu_w_ = 0.001;  // Pa·s (water)
    double mu_o_ = 0.005;  // Pa·s (oil)
};

/**
 * @brief Two-phase immiscible flow solver (Buckley-Leverett type)
 * 
 * Solves the saturation equation for immiscible displacement using
 * either IMPES (Implicit Pressure Explicit Saturation) or fully implicit.
 */
class TwoPhaseFlowSolver {
public:
    TwoPhaseFlowSolver() = default;
    
    // Setup
    void setGrid(int nx, int ny, int nz, double dx, double dy, double dz);
    void setRockProperties(double porosity, double permeability);
    void setFluidProperties(double mu_w, double mu_o, double rho_w, double rho_o);
    void setRelPermModel(const RelPermModel& model);
    void setInitialSaturation(double Sw_init);
    void setBoundaryInjection(double Sw_inj, double rate);
    
    // Solve
    void step(double dt);
    
    // Access
    const std::vector<double>& getPressure() const { return pressure_; }
    const std::vector<double>& getSaturation() const { return saturation_; }
    double getTime() const { return time_; }
    
private:
    // Grid
    int nx_ = 0, ny_ = 0, nz_ = 0;
    double dx_ = 1.0, dy_ = 1.0, dz_ = 1.0;
    int ncells_ = 0;
    
    // Rock properties
    double phi_ = 0.2;
    double k_ = 1e-13;  // m²
    
    // Fluid properties
    double mu_w_ = 0.001;
    double mu_o_ = 0.005;
    double rho_w_ = 1000.0;
    double rho_o_ = 850.0;
    
    // Relative permeability
    FractionalFlow frac_flow_;
    
    // State
    std::vector<double> pressure_;
    std::vector<double> saturation_;
    std::vector<double> saturation_old_;
    double time_ = 0.0;
    
    // Boundary conditions
    double Sw_inj_ = 0.8;
    double inj_rate_ = 0.0;
    
    // Helpers
    int idx(int i, int j, int k) const { return i + j * nx_ + k * nx_ * ny_; }
    void solvePressure();
    void updateSaturation(double dt);
};

/**
 * @brief Coupled flow-geomechanics solver for compaction and subsidence
 */
class CoupledFlowGeomechanicsSolver {
public:
    CoupledFlowGeomechanicsSolver() = default;
    
    // Setup
    void setGrid(int nx, int ny, double Lx, double Ly);
    void setRockProperties(double phi, double k, double E, double nu, double biot = 0.7);
    void setFluidProperties(double rho, double mu, double ct);
    void setInitialPressure(double P0);
    
    // Wells
    void addInjector(double x, double y, double rate);
    void addProducer(double x, double y, double rate);
    
    // Solve
    void step(double dt);
    
    // Access
    const std::vector<std::vector<double>>& getPressure() const { return pressure_; }
    const std::vector<std::vector<double>>& getSubsidence() const { return subsidence_; }
    const std::vector<std::vector<double>>& getPorosity() const { return porosity_; }
    const std::vector<std::vector<double>>& getPermeability() const { return permeability_; }
    double getTime() const { return time_; }
    
private:
    // Grid
    int nx_ = 0, ny_ = 0;
    double Lx_ = 0, Ly_ = 0;
    double dx_ = 0, dy_ = 0;
    
    // Rock properties
    double phi0_ = 0.2;
    double k0_ = 1e-13;
    double E_ = 10e9;
    double nu_ = 0.25;
    double alpha_ = 0.7;  // Biot coefficient
    
    // Fluid properties
    double rho_f_ = 1000.0;
    double mu_ = 0.001;
    double ct_ = 1e-9;
    
    // Initial conditions
    double P_initial_ = 30e6;
    
    // Wells
    struct Well {
        double x, y;
        bool is_injector;
        double rate;
    };
    std::vector<Well> wells_;
    
    // State
    std::vector<std::vector<double>> pressure_;
    std::vector<std::vector<double>> porosity_;
    std::vector<std::vector<double>> permeability_;
    std::vector<std::vector<double>> subsidence_;
    std::vector<std::vector<double>> displacement_x_;
    std::vector<std::vector<double>> displacement_y_;
    double time_ = 0.0;
    
    // Solver methods
    void solveFlow(double dt);
    void solveGeomechanics();
    void updatePoroPermFromGeomechanics();
    void applyWellConditions(std::vector<std::vector<double>>& P);
};

} // namespace FSRM

#endif // TWO_PHASE_FLOW_HPP
