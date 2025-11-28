#include "TracerModel.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <iostream>
#include <iomanip>

namespace FSRM {

// ============================================================================
// TracerBreakthroughCurve Implementation
// ============================================================================

void TracerBreakthroughCurve::computeStatistics(double total_injected_mass) {
    if (times.empty() || concentrations.empty()) return;
    
    // Find breakthrough time (first arrival)
    breakthrough_time = times.front();
    for (size_t i = 0; i < concentrations.size(); ++i) {
        if (concentrations[i] > 0.01 * *std::max_element(concentrations.begin(), 
                                                          concentrations.end())) {
            breakthrough_time = times[i];
            break;
        }
    }
    
    // Find peak
    auto max_it = std::max_element(concentrations.begin(), concentrations.end());
    peak_concentration = *max_it;
    peak_time = times[max_it - concentrations.begin()];
    
    // Calculate mean residence time (first moment)
    double m0 = 0.0, m1 = 0.0, m2 = 0.0;
    for (size_t i = 1; i < times.size(); ++i) {
        double dt = times[i] - times[i-1];
        double c_avg = 0.5 * (concentrations[i] + concentrations[i-1]);
        double t_avg = 0.5 * (times[i] + times[i-1]);
        
        m0 += c_avg * dt;
        m1 += c_avg * t_avg * dt;
        m2 += c_avg * t_avg * t_avg * dt;
    }
    
    mean_residence_time = (m0 > 0.0) ? m1 / m0 : 0.0;
    
    // Variance (second central moment)
    variance = (m0 > 0.0) ? (m2 / m0 - mean_residence_time * mean_residence_time) : 0.0;
    
    // Recovery fraction
    double recovered_mass = cumulative_mass.empty() ? m0 : cumulative_mass.back();
    recovery_fraction = (total_injected_mass > 0.0) 
                      ? recovered_mass / total_injected_mass : 0.0;
}

// ============================================================================
// TracerSolver Implementation
// ============================================================================

TracerSolver::TracerSolver()
    : nx_(0), ny_(0), nz_(0),
      dx_(1.0), dy_(1.0), dz_(1.0),
      ncells_(0), rock_density_(2650.0), time_(0.0) {}

void TracerSolver::addTracer(const TracerProperties& props) {
    tracer_props_[props.name] = props;
    concentrations_[props.name] = std::vector<double>(ncells_, 0.0);
}

void TracerSolver::addInjection(const std::string& tracer_name, 
                                 const TracerInjection& inj) {
    injections_[tracer_name].push_back(inj);
}

void TracerSolver::setGrid(int nx, int ny, int nz, 
                            double dx, double dy, double dz) {
    nx_ = nx; ny_ = ny; nz_ = nz;
    dx_ = dx; dy_ = dy; dz_ = dz;
    ncells_ = nx * ny * nz;
    
    // Resize all arrays
    porosity_.resize(ncells_, 0.2);
    permeability_.resize(ncells_, 1e-13);
    vx_.resize(ncells_, 0.0);
    vy_.resize(ncells_, 0.0);
    vz_.resize(ncells_, 0.0);
    Sw_.resize(ncells_, 1.0);
    So_.resize(ncells_, 0.0);
    Sg_.resize(ncells_, 0.0);
    
    // Resize concentration arrays
    for (auto& [name, conc] : concentrations_) {
        conc.resize(ncells_, 0.0);
    }
}

void TracerSolver::setPorosity(const std::vector<double>& phi) {
    if (phi.size() == static_cast<size_t>(ncells_)) {
        porosity_ = phi;
    }
}

void TracerSolver::setPermeability(const std::vector<double>& k) {
    if (k.size() == static_cast<size_t>(ncells_)) {
        permeability_ = k;
    }
}

void TracerSolver::setRockDensity(double rho) {
    rock_density_ = rho;
}

void TracerSolver::setVelocityField(const std::vector<double>& vx,
                                     const std::vector<double>& vy,
                                     const std::vector<double>& vz) {
    if (vx.size() == static_cast<size_t>(ncells_)) vx_ = vx;
    if (vy.size() == static_cast<size_t>(ncells_)) vy_ = vy;
    if (vz.size() == static_cast<size_t>(ncells_)) vz_ = vz;
}

void TracerSolver::setSaturationField(const std::vector<double>& Sw,
                                       const std::vector<double>& So,
                                       const std::vector<double>& Sg) {
    if (Sw.size() == static_cast<size_t>(ncells_)) Sw_ = Sw;
    if (So.size() == static_cast<size_t>(ncells_)) So_ = So;
    if (Sg.size() == static_cast<size_t>(ncells_)) Sg_ = Sg;
}

const std::vector<double>& TracerSolver::getConcentration(
    const std::string& tracer_name) const {
    static std::vector<double> empty;
    auto it = concentrations_.find(tracer_name);
    return (it != concentrations_.end()) ? it->second : empty;
}

double TracerSolver::getConcentration(const std::string& tracer_name,
                                       int i, int j, int k) const {
    auto it = concentrations_.find(tracer_name);
    if (it == concentrations_.end()) return 0.0;
    
    int index = idx(i, j, k);
    if (index < 0 || index >= ncells_) return 0.0;
    
    return it->second[index];
}

void TracerSolver::addObservationWell(const std::string& name, int i, int j, int k) {
    ObsPoint obs;
    obs.name = name;
    obs.i = i; obs.j = j; obs.k = k;
    obs_points_.push_back(obs);
}

TracerBreakthroughCurve TracerSolver::getBreakthroughCurve(
    const std::string& tracer_name,
    const std::string& well_name) const {
    
    for (const auto& obs : obs_points_) {
        if (obs.name == well_name) {
            auto it = obs.curves.find(tracer_name);
            if (it != obs.curves.end()) {
                return it->second;
            }
        }
    }
    return TracerBreakthroughCurve();
}

void TracerSolver::applyInjections(double t, double dt) {
    for (auto& [tracer_name, inj_list] : injections_) {
        for (const auto& inj : inj_list) {
            // Check if injection is active
            if (t < inj.start_time || t > inj.end_time) continue;
            
            // Find well location (simplified - assume index 0)
            // In real implementation, would look up well location
            int well_idx = 0;  // Placeholder
            
            auto& conc = concentrations_[tracer_name];
            if (well_idx < ncells_) {
                if (inj.use_mass_rate) {
                    // Add mass
                    double cell_volume = dx_ * dy_ * dz_ * porosity_[well_idx];
                    double S = 1.0;  // Assume water tracer for simplicity
                    if (tracer_props_[tracer_name].phase == TracerPhase::WATER) {
                        S = Sw_[well_idx];
                    }
                    double added_conc = inj.mass_rate * dt / (cell_volume * S);
                    conc[well_idx] += added_conc;
                } else {
                    // Set concentration
                    conc[well_idx] = inj.concentration;
                }
            }
        }
    }
}

void TracerSolver::recordObservations() {
    for (auto& obs : obs_points_) {
        int index = idx(obs.i, obs.j, obs.k);
        if (index < 0 || index >= ncells_) continue;
        
        for (const auto& [tracer_name, conc] : concentrations_) {
            auto& curve = obs.curves[tracer_name];
            curve.tracer_name = tracer_name;
            curve.location_name = obs.name;
            curve.times.push_back(time_);
            curve.concentrations.push_back(conc[index]);
            
            // Update cumulative
            if (curve.cumulative_mass.empty()) {
                curve.cumulative_mass.push_back(0.0);
            } else {
                double dt = curve.times.size() > 1 
                          ? curve.times.back() - curve.times[curve.times.size()-2] 
                          : 0.0;
                curve.cumulative_mass.push_back(
                    curve.cumulative_mass.back() + conc[index] * dt
                );
            }
        }
    }
}

void TracerSolver::writeVTK(const std::string& filename, int step) const {
    std::ofstream file(filename);
    if (!file.is_open()) return;
    
    file << "# vtk DataFile Version 3.0\n";
    file << "Tracer concentrations step " << step << "\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_POINTS\n";
    file << "DIMENSIONS " << nx_ << " " << ny_ << " " << nz_ << "\n";
    file << "ORIGIN 0 0 0\n";
    file << "SPACING " << dx_ << " " << dy_ << " " << dz_ << "\n";
    file << "POINT_DATA " << ncells_ << "\n";
    
    for (const auto& [name, conc] : concentrations_) {
        file << "SCALARS " << name << " float 1\n";
        file << "LOOKUP_TABLE default\n";
        for (int k = 0; k < nz_; ++k) {
            for (int j = 0; j < ny_; ++j) {
                for (int i = 0; i < nx_; ++i) {
                    file << conc[idx(i, j, k)] << "\n";
                }
            }
        }
    }
    
    file.close();
}

void TracerSolver::writeSummary(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) return;
    
    file << "Tracer Simulation Summary\n";
    file << "=========================\n\n";
    file << "Time: " << time_ << " s\n\n";
    
    for (const auto& [name, props] : tracer_props_) {
        file << "Tracer: " << name << "\n";
        file << "  Phase: ";
        switch (props.phase) {
            case TracerPhase::WATER: file << "Water\n"; break;
            case TracerPhase::OIL: file << "Oil\n"; break;
            case TracerPhase::GAS: file << "Gas\n"; break;
            default: file << "Free\n"; break;
        }
        
        const auto& conc = concentrations_.at(name);
        double min_c = *std::min_element(conc.begin(), conc.end());
        double max_c = *std::max_element(conc.begin(), conc.end());
        double sum_c = std::accumulate(conc.begin(), conc.end(), 0.0);
        
        file << "  Min concentration: " << min_c << "\n";
        file << "  Max concentration: " << max_c << "\n";
        file << "  Total mass in domain: " << sum_c * dx_ * dy_ * dz_ << "\n\n";
    }
    
    file.close();
}

// ============================================================================
// ExplicitTracerSolver Implementation
// ============================================================================

ExplicitTracerSolver::ExplicitTracerSolver()
    : max_courant_(0.5) {}

double ExplicitTracerSolver::getMaxTimeStep() const {
    // CFL condition: dt < dx / v
    double max_v = 0.0;
    for (int i = 0; i < ncells_; ++i) {
        double v = std::sqrt(vx_[i]*vx_[i] + vy_[i]*vy_[i] + vz_[i]*vz_[i]);
        max_v = std::max(max_v, v);
    }
    
    double min_dx = std::min({dx_, dy_, dz_});
    
    if (max_v < 1e-20) return 1e10;
    return max_courant_ * min_dx / max_v;
}

void ExplicitTracerSolver::step(double dt) {
    // Apply injection conditions
    applyInjections(time_, dt);
    
    // Loop over all tracers
    for (const auto& [name, props] : tracer_props_) {
        // Advection
        advectionStep(name, dt);
        
        // Dispersion
        dispersionStep(name, dt);
        
        // Reactions (decay, adsorption)
        reactionStep(name, dt);
    }
    
    // Update time
    time_ += dt;
    
    // Record at observation points
    recordObservations();
}

void ExplicitTracerSolver::advectionStep(const std::string& tracer_name, double dt) {
    auto& C = concentrations_[tracer_name];
    const auto& props = tracer_props_[tracer_name];
    
    std::vector<double> C_new = C;
    
    for (int k = 0; k < nz_; ++k) {
        for (int j = 0; j < ny_; ++j) {
            for (int i = 0; i < nx_; ++i) {
                int index = idx(i, j, k);
                
                // Get phase saturation for this tracer
                double S = 1.0;
                if (props.phase == TracerPhase::WATER) S = Sw_[index];
                else if (props.phase == TracerPhase::OIL) S = So_[index];
                else if (props.phase == TracerPhase::GAS) S = Sg_[index];
                
                if (S < 1e-10) continue;
                
                double phi = porosity_[index];
                double R = props.getRetardationFactor(phi, rock_density_);
                
                // Upwind scheme
                double flux_x = 0.0, flux_y = 0.0, flux_z = 0.0;
                
                // X-direction
                if (vx_[index] > 0 && i > 0) {
                    flux_x = vx_[index] * C[idx(i-1, j, k)];
                } else if (vx_[index] < 0 && i < nx_-1) {
                    flux_x = vx_[index] * C[idx(i+1, j, k)];
                }
                flux_x -= vx_[index] * C[index];
                
                // Y-direction
                if (vy_[index] > 0 && j > 0) {
                    flux_y = vy_[index] * C[idx(i, j-1, k)];
                } else if (vy_[index] < 0 && j < ny_-1) {
                    flux_y = vy_[index] * C[idx(i, j+1, k)];
                }
                flux_y -= vy_[index] * C[index];
                
                // Z-direction
                if (vz_[index] > 0 && k > 0) {
                    flux_z = vz_[index] * C[idx(i, j, k-1)];
                } else if (vz_[index] < 0 && k < nz_-1) {
                    flux_z = vz_[index] * C[idx(i, j, k+1)];
                }
                flux_z -= vz_[index] * C[index];
                
                // Update
                double dC = dt / (phi * S * R) * (
                    flux_x / dx_ + flux_y / dy_ + flux_z / dz_
                );
                
                C_new[index] = std::max(0.0, C[index] + dC);
            }
        }
    }
    
    C = C_new;
}

void ExplicitTracerSolver::dispersionStep(const std::string& tracer_name, double dt) {
    auto& C = concentrations_[tracer_name];
    const auto& props = tracer_props_[tracer_name];
    
    std::vector<double> C_new = C;
    
    for (int k = 0; k < nz_; ++k) {
        for (int j = 0; j < ny_; ++j) {
            for (int i = 0; i < nx_; ++i) {
                int index = idx(i, j, k);
                
                double phi = porosity_[index];
                double v = std::sqrt(vx_[index]*vx_[index] + 
                                    vy_[index]*vy_[index] + 
                                    vz_[index]*vz_[index]);
                double D = props.getDispersion(v, phi);
                
                // Central difference for Laplacian
                double d2C_dx2 = 0.0, d2C_dy2 = 0.0, d2C_dz2 = 0.0;
                
                if (i > 0 && i < nx_-1) {
                    d2C_dx2 = (C[idx(i+1, j, k)] - 2*C[index] + C[idx(i-1, j, k)]) 
                            / (dx_ * dx_);
                }
                if (j > 0 && j < ny_-1) {
                    d2C_dy2 = (C[idx(i, j+1, k)] - 2*C[index] + C[idx(i, j-1, k)]) 
                            / (dy_ * dy_);
                }
                if (k > 0 && k < nz_-1) {
                    d2C_dz2 = (C[idx(i, j, k+1)] - 2*C[index] + C[idx(i, j, k-1)]) 
                            / (dz_ * dz_);
                }
                
                double laplacian = d2C_dx2 + d2C_dy2 + d2C_dz2;
                C_new[index] += D * dt * laplacian;
                C_new[index] = std::max(0.0, C_new[index]);
            }
        }
    }
    
    C = C_new;
}

void ExplicitTracerSolver::reactionStep(const std::string& tracer_name, double dt) {
    auto& C = concentrations_[tracer_name];
    const auto& props = tracer_props_[tracer_name];
    
    // Decay
    if (props.decay_constant > 0.0) {
        double decay_factor = std::exp(-props.decay_constant * dt);
        for (int i = 0; i < ncells_; ++i) {
            C[i] *= decay_factor;
        }
    }
    
    // Reactive tracer
    if (props.reaction_rate > 0.0 && !props.reaction_product.empty()) {
        auto it = concentrations_.find(props.reaction_product);
        if (it != concentrations_.end()) {
            auto& C_product = it->second;
            double reaction_factor = 1.0 - std::exp(-props.reaction_rate * dt);
            
            for (int i = 0; i < ncells_; ++i) {
                double reacted = C[i] * reaction_factor;
                C[i] -= reacted;
                C_product[i] += reacted;
            }
        }
    }
}

// ============================================================================
// ImplicitTracerSolver Implementation
// ============================================================================

ImplicitTracerSolver::ImplicitTracerSolver()
    : max_iterations_(100), tolerance_(1e-8) {}

void ImplicitTracerSolver::step(double dt) {
    applyInjections(time_, dt);
    
    for (const auto& [name, props] : tracer_props_) {
        solveAdvectionDiffusion(name, dt);
    }
    
    time_ += dt;
    recordObservations();
}

void ImplicitTracerSolver::solveAdvectionDiffusion(const std::string& tracer_name, 
                                                    double dt) {
    // Simplified implicit solver using iterative method
    auto& C = concentrations_[tracer_name];
    const auto& props = tracer_props_[tracer_name];
    
    std::vector<double> C_old = C;
    
    for (int iter = 0; iter < max_iterations_; ++iter) {
        std::vector<double> C_new = C_old;
        double max_change = 0.0;
        
        for (int k = 1; k < nz_-1; ++k) {
            for (int j = 1; j < ny_-1; ++j) {
                for (int i = 1; i < nx_-1; ++i) {
                    int index = idx(i, j, k);
                    
                    double phi = porosity_[index];
                    double v = std::sqrt(vx_[index]*vx_[index] + 
                                        vy_[index]*vy_[index] + 
                                        vz_[index]*vz_[index]);
                    double D = props.getDispersion(v, phi);
                    double R = props.getRetardationFactor(phi, rock_density_);
                    
                    // Implicit central difference
                    double ax = D * dt / (dx_ * dx_ * phi * R);
                    double ay = D * dt / (dy_ * dy_ * phi * R);
                    double az = D * dt / (dz_ * dz_ * phi * R);
                    
                    double denom = 1.0 + 2.0 * (ax + ay + az);
                    
                    C_new[index] = (C_old[index] + 
                                   ax * (C[idx(i+1,j,k)] + C[idx(i-1,j,k)]) +
                                   ay * (C[idx(i,j+1,k)] + C[idx(i,j-1,k)]) +
                                   az * (C[idx(i,j,k+1)] + C[idx(i,j,k-1)])) / denom;
                    
                    max_change = std::max(max_change, 
                                          std::abs(C_new[index] - C[index]));
                }
            }
        }
        
        C = C_new;
        
        if (max_change < tolerance_) break;
    }
}

// ============================================================================
// StreamlineTracerSolver Implementation
// ============================================================================

StreamlineTracerSolver::StreamlineTracerSolver()
    : num_streamlines_(100) {}

void StreamlineTracerSolver::step(double dt) {
    applyInjections(time_, dt);
    
    // Trace streamlines if not done
    if (streamlines_.empty()) {
        traceStreamlines();
    }
    
    // Solve 1D problem along each streamline
    for (const auto& [name, props] : tracer_props_) {
        for (auto& sl : streamlines_) {
            solve1D(sl, name, dt);
        }
    }
    
    time_ += dt;
    recordObservations();
}

void StreamlineTracerSolver::traceStreamlines() {
    streamlines_.clear();
    
    // Simple streamline tracing from injection points
    // In a real implementation, would use proper streamline tracing
    
    for (int n = 0; n < num_streamlines_; ++n) {
        Streamline sl;
        
        // Start from random point on inlet face
        double y = (double(n) + 0.5) / num_streamlines_ * ny_ * dy_;
        double z = 0.5 * nz_ * dz_;
        double x = 0.0;
        
        sl.x.push_back(x);
        sl.y.push_back(y);
        sl.z.push_back(z);
        sl.tof.push_back(0.0);
        
        // Trace through domain
        while (x < (nx_ - 1) * dx_) {
            int i = std::min(nx_-1, std::max(0, int(x / dx_)));
            int j = std::min(ny_-1, std::max(0, int(y / dy_)));
            int k = std::min(nz_-1, std::max(0, int(z / dz_)));
            int index = idx(i, j, k);
            
            double vx_local = vx_[index];
            double vy_local = vy_[index];
            double vz_local = vz_[index];
            double v_mag = std::sqrt(vx_local*vx_local + vy_local*vy_local + vz_local*vz_local);
            
            if (v_mag < 1e-20) break;
            
            // Step along streamline
            double ds = 0.1 * dx_;  // Step size
            double dt_step = ds / v_mag;
            
            x += vx_local / v_mag * ds;
            y += vy_local / v_mag * ds;
            z += vz_local / v_mag * ds;
            
            // Keep in bounds
            y = std::max(0.0, std::min(y, (ny_-1) * dy_));
            z = std::max(0.0, std::min(z, (nz_-1) * dz_));
            
            sl.x.push_back(x);
            sl.y.push_back(y);
            sl.z.push_back(z);
            sl.tof.push_back(sl.tof.back() + dt_step);
            sl.cells.push_back(index);
        }
        
        if (sl.x.size() > 1) {
            streamlines_.push_back(sl);
        }
    }
}

void StreamlineTracerSolver::solve1D(Streamline& sl, 
                                      const std::string& tracer_name, 
                                      double dt) {
    if (sl.x.size() < 2) return;
    
    auto& C = concentrations_[tracer_name];
    const auto& props = tracer_props_[tracer_name];
    
    // 1D advection along streamline using time-of-flight
    std::vector<double> C_sl(sl.x.size(), 0.0);
    
    // Map concentrations to streamline
    for (size_t n = 0; n < sl.cells.size(); ++n) {
        C_sl[n] = C[sl.cells[n]];
    }
    
    // Solve 1D advection in time-of-flight coordinates
    // Simple upwind in TOF space
    std::vector<double> C_sl_new = C_sl;
    
    for (size_t n = 1; n < sl.x.size(); ++n) {
        double d_tof = sl.tof[n] - sl.tof[n-1];
        if (d_tof > 1e-20) {
            C_sl_new[n] = C_sl[n-1];  // Perfect advection in TOF coordinates
        }
    }
    
    // Map back to grid
    for (size_t n = 0; n < sl.cells.size(); ++n) {
        // Average with existing value (multiple streamlines per cell)
        C[sl.cells[n]] = 0.5 * (C[sl.cells[n]] + C_sl_new[n]);
    }
}

// ============================================================================
// TracerAnalysis Implementation
// ============================================================================

double TracerAnalysis::estimateSorFromSWCT(const TracerBreakthroughCurve& reactive_curve,
                                            const TracerBreakthroughCurve& passive_curve,
                                            double partition_coeff) {
    // SWCT method: Sor = K * (t_reactive - t_passive) / (t_passive * (1 + K))
    // where K is the oil-water partition coefficient
    
    // Use mean residence times
    double t_r = reactive_curve.mean_residence_time;
    double t_p = passive_curve.mean_residence_time;
    
    if (t_p <= 0.0 || partition_coeff <= 0.0) return 0.0;
    
    double beta = t_r / t_p;
    double Sor = (beta - 1.0) / (partition_coeff * (1.0 + beta));
    
    return std::max(0.0, std::min(Sor, 1.0));
}

std::map<std::string, double> TracerAnalysis::computeConnectivity(
    const std::vector<TracerBreakthroughCurve>& curves) {
    
    std::map<std::string, double> connectivity;
    
    // Calculate total recovered mass
    double total_recovery = 0.0;
    for (const auto& curve : curves) {
        total_recovery += curve.recovery_fraction;
    }
    
    if (total_recovery <= 0.0) return connectivity;
    
    // Fraction of tracer recovered at each location
    for (const auto& curve : curves) {
        connectivity[curve.location_name] = curve.recovery_fraction / total_recovery;
    }
    
    return connectivity;
}

double TracerAnalysis::estimateSweptVolume(const TracerBreakthroughCurve& curve,
                                            double injection_rate) {
    // Swept volume = injection rate * mean residence time
    return injection_rate * curve.mean_residence_time;
}

double TracerAnalysis::computeLorenzCoefficient(const TracerBreakthroughCurve& curve) {
    // Lorenz coefficient from flow capacity - storage capacity plot
    // Higher values indicate more heterogeneity
    
    if (curve.times.empty() || curve.concentrations.empty()) return 0.0;
    
    // Create F-Φ curve (flow capacity vs storage capacity)
    std::vector<double> F, Phi;
    
    double sum_C = 0.0;
    double sum_t = 0.0;
    
    for (size_t i = 0; i < curve.times.size(); ++i) {
        sum_C += curve.concentrations[i];
        sum_t += curve.times[i];
    }
    
    if (sum_C <= 0.0 || sum_t <= 0.0) return 0.0;
    
    double cumul_C = 0.0;
    double cumul_t = 0.0;
    
    for (size_t i = 0; i < curve.times.size(); ++i) {
        cumul_C += curve.concentrations[i];
        cumul_t += curve.times[i];
        
        F.push_back(cumul_C / sum_C);
        Phi.push_back(cumul_t / sum_t);
    }
    
    // Lorenz coefficient = 2 * (area between F-Φ curve and diagonal)
    double area = 0.0;
    for (size_t i = 1; i < F.size(); ++i) {
        double dPhi = Phi[i] - Phi[i-1];
        double avg_F = 0.5 * (F[i] + F[i-1]);
        double avg_diag = 0.5 * (Phi[i] + Phi[i-1]);
        area += (avg_F - avg_diag) * dPhi;
    }
    
    return std::abs(2.0 * area);
}

double TracerAnalysis::fitPecletNumber(const TracerBreakthroughCurve& curve,
                                        double distance) {
    // Fit Pe = v*L/D from breakthrough curve variance
    // For 1D advection-dispersion: σ² = 2*L²/Pe
    
    if (curve.variance <= 0.0 || distance <= 0.0) return 1e10;
    
    double sigma_normalized = std::sqrt(curve.variance) / curve.mean_residence_time;
    double Pe = 2.0 / (sigma_normalized * sigma_normalized);
    
    return Pe;
}

// ============================================================================
// TracerManager Implementation
// ============================================================================

TracerManager::TracerManager()
    : time_(0.0), initialized_(false) {
    solver_ = std::make_unique<ExplicitTracerSolver>();
}

void TracerManager::configure(const std::map<std::string, std::string>& config) {
    // Parse configuration
    if (config.count("solver_type")) {
        const std::string& type = config.at("solver_type");
        if (type == "IMPLICIT") {
            solver_ = std::make_unique<ImplicitTracerSolver>();
        } else if (type == "STREAMLINE") {
            solver_ = std::make_unique<StreamlineTracerSolver>();
        } else {
            solver_ = std::make_unique<ExplicitTracerSolver>();
        }
    }
}

void TracerManager::setGrid(int nx, int ny, int nz, 
                            double dx, double dy, double dz) {
    solver_->setGrid(nx, ny, nz, dx, dy, dz);
    initialized_ = true;
}

void TracerManager::setPorosity(const std::vector<double>& phi) {
    solver_->setPorosity(phi);
}

void TracerManager::setPermeability(const std::vector<double>& k) {
    solver_->setPermeability(k);
}

void TracerManager::setRockDensity(double rho) {
    solver_->setRockDensity(rho);
}

void TracerManager::addWaterTracer(const std::string& name, TracerBehavior behavior) {
    TracerProperties props;
    props.name = name;
    props.phase = TracerPhase::WATER;
    props.behavior = behavior;
    tracers_.push_back(props);
    solver_->addTracer(props);
}

void TracerManager::addOilTracer(const std::string& name, TracerBehavior behavior) {
    TracerProperties props;
    props.name = name;
    props.phase = TracerPhase::OIL;
    props.behavior = behavior;
    tracers_.push_back(props);
    solver_->addTracer(props);
}

void TracerManager::addGasTracer(const std::string& name, TracerBehavior behavior) {
    TracerProperties props;
    props.name = name;
    props.phase = TracerPhase::GAS;
    props.behavior = behavior;
    tracers_.push_back(props);
    solver_->addTracer(props);
}

void TracerManager::addPartitioningTracer(const std::string& name, double Kow) {
    TracerProperties props;
    props.name = name;
    props.phase = TracerPhase::WATER;
    props.behavior = TracerBehavior::PARTITIONING;
    props.partition_coeff_ow = Kow;
    tracers_.push_back(props);
    solver_->addTracer(props);
}

void TracerManager::addRadioactiveTracer(const std::string& name, double half_life_days) {
    TracerProperties props;
    props.name = name;
    props.phase = TracerPhase::WATER;
    props.behavior = TracerBehavior::DECAYING;
    props.half_life = half_life_days * 86400.0;  // Convert to seconds
    props.decay_constant = std::log(2.0) / props.half_life;
    tracers_.push_back(props);
    solver_->addTracer(props);
}

void TracerManager::addInjection(const std::string& tracer_name,
                                  const std::string& well_name,
                                  double concentration,
                                  double start_time, double end_time) {
    TracerInjection inj;
    inj.well_name = well_name;
    inj.concentration = concentration;
    inj.start_time = start_time;
    inj.end_time = end_time;
    inj.use_mass_rate = false;
    solver_->addInjection(tracer_name, inj);
}

void TracerManager::addSlugInjection(const std::string& tracer_name,
                                      const std::string& well_name,
                                      double mass, double duration) {
    TracerInjection inj;
    inj.well_name = well_name;
    inj.mass_rate = mass / duration;
    inj.start_time = time_;
    inj.end_time = time_ + duration;
    inj.use_mass_rate = true;
    solver_->addInjection(tracer_name, inj);
}

void TracerManager::updateFlowField(const std::vector<double>& vx,
                                     const std::vector<double>& vy,
                                     const std::vector<double>& vz,
                                     const std::vector<double>& Sw,
                                     const std::vector<double>& So,
                                     const std::vector<double>& Sg) {
    solver_->setVelocityField(vx, vy, vz);
    solver_->setSaturationField(Sw, So, Sg);
}

void TracerManager::step(double dt) {
    solver_->step(dt);
    time_ += dt;
}

std::vector<double> TracerManager::getConcentration(const std::string& tracer_name) const {
    return solver_->getConcentration(tracer_name);
}

TracerBreakthroughCurve TracerManager::getBreakthrough(const std::string& tracer_name,
                                                        const std::string& location) const {
    return solver_->getBreakthroughCurve(tracer_name, location);
}

void TracerManager::addObservationWell(const std::string& name, int i, int j, int k) {
    solver_->addObservationWell(name, i, j, k);
}

void TracerManager::writeOutput(const std::string& prefix, int step) const {
    solver_->writeVTK(prefix + "_" + std::to_string(step) + ".vtk", step);
}

void TracerManager::writeSummary(const std::string& filename) const {
    solver_->writeSummary(filename);
}

double TracerManager::estimateResidualOil(const std::string& partitioning_tracer,
                                           const std::string& passive_tracer,
                                           const std::string& well) const {
    auto part_curve = getBreakthrough(partitioning_tracer, well);
    auto pass_curve = getBreakthrough(passive_tracer, well);
    
    // Find partition coefficient
    double Kow = 1.0;
    for (const auto& t : tracers_) {
        if (t.name == partitioning_tracer) {
            Kow = t.partition_coeff_ow;
            break;
        }
    }
    
    return TracerAnalysis::estimateSorFromSWCT(part_curve, pass_curve, Kow);
}

} // namespace FSRM
