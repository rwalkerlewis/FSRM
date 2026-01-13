#include "TwoPhaseFlow.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

namespace FSRM {

// ============================================================================
// TwoPhaseFlowSolver Implementation
// ============================================================================

void TwoPhaseFlowSolver::setGrid(int nx, int ny, int nz, double dx, double dy, double dz) {
    nx_ = nx;
    ny_ = ny;
    nz_ = nz;
    dx_ = dx;
    dy_ = dy;
    dz_ = dz;
    ncells_ = nx * ny * nz;
    
    pressure_.resize(ncells_, 0.0);
    saturation_.resize(ncells_, 0.0);
    saturation_old_.resize(ncells_, 0.0);
}

void TwoPhaseFlowSolver::setRockProperties(double porosity, double permeability) {
    phi_ = porosity;
    k_ = permeability;
}

void TwoPhaseFlowSolver::setFluidProperties(double mu_w, double mu_o, double rho_w, double rho_o) {
    mu_w_ = mu_w;
    mu_o_ = mu_o;
    rho_w_ = rho_w;
    rho_o_ = rho_o;
    frac_flow_.setViscosities(mu_w, mu_o);
}

void TwoPhaseFlowSolver::setRelPermModel(const RelPermModel& model) {
    frac_flow_.setRelPermModel(model);
}

void TwoPhaseFlowSolver::setInitialSaturation(double Sw_init) {
    std::fill(saturation_.begin(), saturation_.end(), Sw_init);
    saturation_old_ = saturation_;
}

void TwoPhaseFlowSolver::setBoundaryInjection(double Sw_inj, double rate) {
    Sw_inj_ = Sw_inj;
    inj_rate_ = rate;
}

void TwoPhaseFlowSolver::step(double dt) {
    saturation_old_ = saturation_;
    
    // IMPES: solve pressure then update saturation
    solvePressure();
    updateSaturation(dt);
    
    time_ += dt;
}

void TwoPhaseFlowSolver::solvePressure() {
    // Simplified single-phase pressure solve for now
    // In production, this would be a full multiphase pressure solve
}

void TwoPhaseFlowSolver::updateSaturation(double dt) {
    // Upwind finite volume for saturation equation
    // dSw/dt + (1/phi) * d(u*fw)/dx = 0
    
    std::vector<double> Sw_new = saturation_;
    
    for (int k = 0; k < nz_; ++k) {
        for (int j = 0; j < ny_; ++j) {
            for (int i = 0; i < nx_; ++i) {
                int cell = idx(i, j, k);
                
                // Injection boundary at i=0
                if (i == 0) {
                    Sw_new[cell] = Sw_inj_;
                    continue;
                }
                
                // Production boundary at i=nx-1 (outflow)
                if (i == nx_ - 1) {
                    Sw_new[cell] = saturation_[idx(i-1, j, k)];
                    continue;
                }
                
                // Interior cell - upwind scheme
                double fw_left = frac_flow_.fw(saturation_[idx(i-1, j, k)]);
                double fw_right = frac_flow_.fw(saturation_[cell]);
                
                // Total velocity (Darcy) - simplified
                double u_total = inj_rate_;  // Constant injection rate
                
                double flux_in = u_total * fw_left;
                double flux_out = u_total * fw_right;
                
                double dSw = -(dt / (phi_ * dx_)) * (flux_out - flux_in);
                Sw_new[cell] = saturation_[cell] + dSw;
                
                // Clamp to physical bounds
                Sw_new[cell] = std::max(0.2, std::min(0.8, Sw_new[cell]));
            }
        }
    }
    
    saturation_ = Sw_new;
}

// ============================================================================
// CoupledFlowGeomechanicsSolver Implementation
// ============================================================================

void CoupledFlowGeomechanicsSolver::setGrid(int nx, int ny, double Lx, double Ly) {
    nx_ = nx;
    ny_ = ny;
    Lx_ = Lx;
    Ly_ = Ly;
    dx_ = Lx / nx;
    dy_ = Ly / ny;
    
    // Initialize 2D arrays
    pressure_.resize(ny, std::vector<double>(nx, 0.0));
    porosity_.resize(ny, std::vector<double>(nx, 0.0));
    permeability_.resize(ny, std::vector<double>(nx, 0.0));
    subsidence_.resize(ny, std::vector<double>(nx, 0.0));
    displacement_x_.resize(ny, std::vector<double>(nx, 0.0));
    displacement_y_.resize(ny, std::vector<double>(nx, 0.0));
}

void CoupledFlowGeomechanicsSolver::setRockProperties(double phi, double k, double E, double nu, double biot) {
    phi0_ = phi;
    k0_ = k;
    E_ = E;
    nu_ = nu;
    alpha_ = biot;
    
    // Initialize porosity and permeability
    for (int j = 0; j < ny_; ++j) {
        for (int i = 0; i < nx_; ++i) {
            porosity_[j][i] = phi0_;
            permeability_[j][i] = k0_;
        }
    }
}

void CoupledFlowGeomechanicsSolver::setFluidProperties(double rho, double mu, double ct) {
    rho_f_ = rho;
    mu_ = mu;
    ct_ = ct;
}

void CoupledFlowGeomechanicsSolver::setInitialPressure(double P0) {
    P_initial_ = P0;
    for (int j = 0; j < ny_; ++j) {
        for (int i = 0; i < nx_; ++i) {
            pressure_[j][i] = P0;
        }
    }
}

void CoupledFlowGeomechanicsSolver::addInjector(double x, double y, double rate) {
    wells_.push_back({x, y, true, rate});
}

void CoupledFlowGeomechanicsSolver::addProducer(double x, double y, double rate) {
    wells_.push_back({x, y, false, std::abs(rate)});
}

void CoupledFlowGeomechanicsSolver::step(double dt) {
    solveFlow(dt);
    solveGeomechanics();
    updatePoroPermFromGeomechanics();
    time_ += dt;
}

void CoupledFlowGeomechanicsSolver::solveFlow(double dt) {
    auto P_new = pressure_;
    
    // Simple explicit finite difference
    for (int j = 1; j < ny_ - 1; ++j) {
        for (int i = 1; i < nx_ - 1; ++i) {
            double P = pressure_[j][i];
            double Px_plus = pressure_[j][i+1];
            double Px_minus = pressure_[j][i-1];
            double Py_plus = pressure_[j+1][i];
            double Py_minus = pressure_[j-1][i];
            
            double phi = porosity_[j][i];
            double k = permeability_[j][i];
            
            // Darcy flux
            double qx = -(k / mu_) * (Px_plus - Px_minus) / (2 * dx_);
            double qy = -(k / mu_) * (Py_plus - Py_minus) / (2 * dy_);
            
            // Divergence
            double div_q = (qx / dx_) + (qy / dy_);
            
            // Pressure update (compressible flow)
            double dP = -(dt / (phi * ct_)) * div_q;
            P_new[j][i] = P + dP;
        }
    }
    
    // Apply well conditions
    applyWellConditions(P_new);
    
    pressure_ = P_new;
}

void CoupledFlowGeomechanicsSolver::solveGeomechanics() {
    // Compute subsidence from pressure depletion
    double K = E_ / (3 * (1 - 2 * nu_));  // Bulk modulus
    
    for (int j = 0; j < ny_; ++j) {
        for (int i = 0; i < nx_; ++i) {
            double dP = pressure_[j][i] - P_initial_;
            
            // Uniaxial compaction model
            double volumetric_strain = alpha_ * dP / K;
            
            // Vertical displacement (subsidence)
            double cell_height = Ly_ / ny_;
            double dz = volumetric_strain * cell_height;
            
            subsidence_[j][i] = -dz;  // Negative is downward
            
            // Horizontal displacement
            double lateral_strain = nu_ / (1 - nu_) * volumetric_strain;
            displacement_x_[j][i] = lateral_strain * (i * dx_ - Lx_ / 2);
            displacement_y_[j][i] = lateral_strain * (j * dy_ - Ly_ / 2);
        }
    }
}

void CoupledFlowGeomechanicsSolver::updatePoroPermFromGeomechanics() {
    double K = E_ / (3 * (1 - 2 * nu_));
    
    for (int j = 0; j < ny_; ++j) {
        for (int i = 0; i < nx_; ++i) {
            double dP = pressure_[j][i] - P_initial_;
            
            // Porosity change from pressure
            double volumetric_strain = alpha_ * dP / K;
            porosity_[j][i] = phi0_ * (1 + volumetric_strain);
            
            // Permeability change (Kozeny-Carman type)
            double phi = porosity_[j][i];
            double phi_ratio = phi / phi0_;
            permeability_[j][i] = k0_ * std::pow(phi_ratio, 3) * 
                                  std::pow((1 - phi0_) / (1 - phi), 2);
        }
    }
}

void CoupledFlowGeomechanicsSolver::applyWellConditions(std::vector<std::vector<double>>& P) {
    for (const auto& well : wells_) {
        int i = static_cast<int>(well.x / dx_);
        int j = static_cast<int>(well.y / dy_);
        
        i = std::max(0, std::min(nx_ - 1, i));
        j = std::max(0, std::min(ny_ - 1, j));
        
        if (well.is_injector) {
            P[j][i] = P_initial_ + 5e6;  // 5 MPa overpressure
        } else {
            P[j][i] = P_initial_ - 10e6;  // 10 MPa drawdown
        }
    }
}

} // namespace FSRM
