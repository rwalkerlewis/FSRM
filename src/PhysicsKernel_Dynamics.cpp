/*
 * Dynamic Wave Physics Kernels
 * 
 * Implementation of elastodynamics and poroelastodynamics with:
 * - Full wave propagation with inertia
 * - Static-to-dynamic triggering for stress-induced seismicity
 * - Dynamic permeability changes from transient waves
 */

#include "PhysicsKernel.hpp"
#include <cmath>
#include <algorithm>

namespace ResSim {

// ============================================================================
// Elastodynamics Kernel
// ============================================================================

ElastodynamicsKernel::ElastodynamicsKernel()
    : PhysicsKernel(PhysicsType::ELASTODYNAMICS),
      youngs_modulus(10e9), poisson_ratio(0.25), density(2500.0),
      p_wave_velocity(5000.0), s_wave_velocity(3000.0), quality_factor(100.0),
      damping_alpha(0.01), damping_beta(0.001),
      enable_static_triggering(false), trigger_threshold(1.0e6),
      event_duration(10.0), trigger_time(-1.0), event_active(false) {}

PetscErrorCode ElastodynamicsKernel::setup(DM dm, PetscFE fe) {
    PetscFunctionBeginUser;
    // Setup finite element spaces for dynamic problem
    PetscFunctionReturn(0);
}

void ElastodynamicsKernel::residual(const PetscScalar u[], const PetscScalar u_t[],
                                    const PetscScalar u_x[], const PetscScalar a[],
                                    const PetscReal x[], PetscScalar f[]) {
    // Elastodynamics equation: ρ ∂²u/∂t² = ∇·σ + ρg
    // In first-order form with v = ∂u/∂t:
    // ρ ∂v/∂t = ∇·σ + ρg - damping
    
    // u[0..2] = displacement (ux, uy, uz)
    // u_t[0..2] = velocity (vx, vy, vz) 
    // u_x[0..8] = displacement gradients
    
    // Compute strain tensor from displacement gradient
    double eps_xx = u_x[0];  // du_x/dx
    double eps_yy = u_x[4];  // du_y/dy
    double eps_zz = u_x[8];  // du_z/dz
    double eps_xy = 0.5 * (u_x[1] + u_x[3]);  // 0.5*(du_x/dy + du_y/dx)
    double eps_xz = 0.5 * (u_x[2] + u_x[6]);
    double eps_yz = 0.5 * (u_x[5] + u_x[7]);
    
    // Lame parameters from E and ν
    double lambda = youngs_modulus * poisson_ratio / 
                   ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio));
    double mu = youngs_modulus / (2.0 * (1.0 + poisson_ratio));
    
    // Stress tensor (linear elastic)
    double trace = eps_xx + eps_yy + eps_zz;
    double sigma_xx = lambda * trace + 2.0 * mu * eps_xx;
    double sigma_yy = lambda * trace + 2.0 * mu * eps_yy;
    double sigma_zz = lambda * trace + 2.0 * mu * eps_zz;
    double sigma_xy = 2.0 * mu * eps_xy;
    double sigma_xz = 2.0 * mu * eps_xz;
    double sigma_yz = 2.0 * mu * eps_yz;
    
    // Inertial term: ρ ∂²u/∂t² (using second time derivative)
    // In weak form, this becomes mass matrix times acceleration
    // Here u_t is velocity, so we need acceleration
    // This is handled by the time integrator
    
    // f0 term: accumulation (inertia)
    f[0] = density * u_t[0];  // ρ * vx (mass times velocity derivative)
    f[1] = density * u_t[1];  // ρ * vy
    f[2] = density * u_t[2];  // ρ * vz
    
    // Add Rayleigh damping: C = α*M + β*K
    // Mass proportional damping: α * ρ * v
    f[0] += damping_alpha * density * u_t[0];
    f[1] += damping_alpha * density * u_t[1];
    f[2] += damping_alpha * density * u_t[2];
    
    // Stiffness proportional damping is handled in f1 term
    
    // Note: f1 term (stress divergence) would be handled separately
    // in the actual PETSc weak form assembly
}

void ElastodynamicsKernel::jacobian(const PetscScalar u[], const PetscScalar u_t[],
                                    const PetscScalar u_x[], const PetscScalar a[],
                                    const PetscReal x[], PetscScalar J[]) {
    // Jacobian for elastodynamics
    // ∂F/∂u_t = ρ * I + α * ρ * I (with shift parameter 'a' from time integrator)
    
    double mass_term = density * (1.0 + damping_alpha);
    
    // Diagonal blocks for each component
    J[0] = mass_term * a[0];   // ∂F_x/∂vx
    J[4] = mass_term * a[0];   // ∂F_y/∂vy
    J[8] = mass_term * a[0];   // ∂F_z/∂vz
    
    // Off-diagonal terms from elasticity would be added in full Jacobian
}

void ElastodynamicsKernel::setMaterialProperties(double E, double nu, double rho) {
    youngs_modulus = E;
    poisson_ratio = nu;
    density = rho;
    
    // Compute wave velocities from elastic constants
    double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    double mu = E / (2.0 * (1.0 + nu));
    
    p_wave_velocity = std::sqrt((lambda + 2.0 * mu) / rho);
    s_wave_velocity = std::sqrt(mu / rho);
}

void ElastodynamicsKernel::setWaveProperties(double vp, double vs, double Q) {
    p_wave_velocity = vp;
    s_wave_velocity = vs;
    quality_factor = Q;
    
    // Set damping based on quality factor
    // Q = ω / (2 * ζ) where ζ is damping ratio
    // Simplified: use Q to estimate damping
    double omega_ref = 2.0 * M_PI;  // Reference frequency 1 Hz
    damping_alpha = omega_ref / (2.0 * Q);
}

void ElastodynamicsKernel::setDamping(double alpha, double beta) {
    damping_alpha = alpha;
    damping_beta = beta;
}

void ElastodynamicsKernel::setStaticTriggeringMode(bool enable, double threshold, double duration) {
    enable_static_triggering = enable;
    trigger_threshold = threshold;
    event_duration = duration;
}

bool ElastodynamicsKernel::checkStaticTrigger(const PetscScalar stress[], double current_time) {
    if (!enable_static_triggering || event_active) {
        return false;
    }
    
    // Compute von Mises stress
    double s_xx = stress[0];
    double s_yy = stress[1];
    double s_zz = stress[2];
    double s_xy = stress[3];
    double s_xz = stress[4];
    double s_yz = stress[5];
    
    double s_mean = (s_xx + s_yy + s_zz) / 3.0;
    double dev_xx = s_xx - s_mean;
    double dev_yy = s_yy - s_mean;
    double dev_zz = s_zz - s_mean;
    
    double von_mises = std::sqrt(1.5 * (dev_xx*dev_xx + dev_yy*dev_yy + dev_zz*dev_zz +
                                        2.0 * (s_xy*s_xy + s_xz*s_xz + s_yz*s_yz)));
    
    if (von_mises >= trigger_threshold) {
        trigger_time = current_time;
        event_active = true;
        return true;
    }
    
    return false;
}

bool ElastodynamicsKernel::isInDynamicEvent(double current_time) const {
    if (!enable_static_triggering || !event_active) {
        return false;
    }
    
    return (current_time - trigger_time) <= event_duration;
}

// ============================================================================
// Poroelastodynamics Kernel
// ============================================================================

PoroelastodynamicsKernel::PoroelastodynamicsKernel()
    : PhysicsKernel(PhysicsType::POROELASTODYNAMICS),
      youngs_modulus(10e9), poisson_ratio(0.25), density_solid(2500.0), porosity(0.2),
      density_fluid(1000.0), viscosity(0.001), bulk_modulus_fluid(2.2e9),
      biot_coefficient(1.0), biot_modulus(1e10),
      p_wave_fast(5500.0), s_wave(3000.0), p_wave_slow(1500.0),
      damping_alpha(0.01), damping_beta(0.001),
      permeability(100e-15), permeability_initial(100e-15),
      enable_dynamic_k(false), k_strain_coeff(1e-7), k_stress_coeff(1e-15),
      k_recovery_time(100.0),
      enable_static_triggering(false), trigger_threshold(1.0e6),
      event_duration(10.0), trigger_time(-1.0), event_active(false) {}

PetscErrorCode PoroelastodynamicsKernel::setup(DM dm, PetscFE fe) {
    PetscFunctionBeginUser;
    // Setup finite element spaces for coupled problem
    PetscFunctionReturn(0);
}

void PoroelastodynamicsKernel::residual(const PetscScalar u[], const PetscScalar u_t[],
                                        const PetscScalar u_x[], const PetscScalar a[],
                                        const PetscReal x[], PetscScalar f[]) {
    // Biot's poroelastodynamics equations:
    // 1) Momentum balance: ρ ∂²u/∂t² + ρ_f ∂²w/∂t² = ∇·σ - α∇p
    // 2) Darcy's law (dynamic): ρ_f ∂²u/∂t² + (ρ_f/φ) ∂²w/∂t² = -∇p - (μ/k)∂w/∂t
    // 3) Mass conservation: ∂/∂t(α∇·u + p/M) + ∇·(∂w/∂t) = 0
    // 
    // where:
    //   u = solid displacement
    //   w = relative fluid displacement w = φ(u_f - u_s)
    //   p = pore pressure
    //   ρ = bulk density = (1-φ)ρ_s + φρ_f
    
    // Field ordering: u[0..2] = displacement, u[3] = pressure
    
    double ux = u[0], uy = u[1], uz = u[2];
    double p = u[3];
    
    // Displacement gradients
    double eps_xx = u_x[0];
    double eps_yy = u_x[4];
    double eps_zz = u_x[8];
    double eps_xy = 0.5 * (u_x[1] + u_x[3]);
    double eps_xz = 0.5 * (u_x[2] + u_x[6]);
    double eps_yz = 0.5 * (u_x[5] + u_x[7]);
    
    // Pressure gradient
    double dp_dx = u_x[9];
    double dp_dy = u_x[10];
    double dp_dz = u_x[11];
    
    // Material parameters
    double lambda = youngs_modulus * poisson_ratio / 
                   ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio));
    double mu = youngs_modulus / (2.0 * (1.0 + poisson_ratio));
    double rho_bulk = (1.0 - porosity) * density_solid + porosity * density_fluid;
    
    // Effective stress
    double trace = eps_xx + eps_yy + eps_zz;
    double sigma_xx = lambda * trace + 2.0 * mu * eps_xx - biot_coefficient * p;
    double sigma_yy = lambda * trace + 2.0 * mu * eps_yy - biot_coefficient * p;
    double sigma_zz = lambda * trace + 2.0 * mu * eps_zz - biot_coefficient * p;
    
    // === Momentum equation ===
    // Inertia: ρ ∂²u/∂t²
    f[0] = rho_bulk * u_t[0];  // x-component
    f[1] = rho_bulk * u_t[1];  // y-component
    f[2] = rho_bulk * u_t[2];  // z-component
    
    // Damping
    f[0] += damping_alpha * rho_bulk * u_t[0];
    f[1] += damping_alpha * rho_bulk * u_t[1];
    f[2] += damping_alpha * rho_bulk * u_t[2];
    
    // === Fluid mass conservation (pressure equation) ===
    // ∂/∂t(α∇·u + p/M) + ∇·q = 0
    // where q = -(k/μ)∇p (Darcy flux)
    
    // Storage term
    double storage = biot_coefficient * (u_t[0] + u_t[1] + u_t[2]) + u_t[3] / biot_modulus;
    f[3] = storage;
    
    // Darcy flow with dynamic permeability
    double k_current = permeability;
    if (enable_dynamic_k) {
        // Compute permeability change based on strain/stress
        double strain_amplitude = std::sqrt(eps_xx*eps_xx + eps_yy*eps_yy + eps_zz*eps_zz);
        double stress_amplitude = std::sqrt(sigma_xx*sigma_xx + sigma_yy*sigma_yy + sigma_zz*sigma_zz);
        
        double k_change = k_strain_coeff * strain_amplitude + k_stress_coeff * stress_amplitude;
        k_current = permeability_initial * (1.0 + k_change);
        
        // Clamp to reasonable bounds
        k_current = std::max(1e-18, std::min(k_current, 1e-10));
    }
    
    double mobility = k_current / viscosity;
    
    // Add flux divergence (handled in f1 term in actual implementation)
    // f[3] += mobility * (dp_dx² + dp_dy² + dp_dz²)
}

void PoroelastodynamicsKernel::jacobian(const PetscScalar u[], const PetscScalar u_t[],
                                        const PetscScalar u_x[], const PetscScalar a[],
                                        const PetscReal x[], PetscScalar J[]) {
    // Jacobian for poroelastodynamics (4x4 system in 1D, larger in 3D)
    
    double rho_bulk = (1.0 - porosity) * density_solid + porosity * density_fluid;
    double mass_term = rho_bulk * (1.0 + damping_alpha);
    
    // Displacement-displacement blocks (3x3)
    J[0] = mass_term * a[0];   // ∂F_ux/∂vx
    J[5] = mass_term * a[0];   // ∂F_uy/∂vy
    J[10] = mass_term * a[0];  // ∂F_uz/∂vz
    
    // Pressure-pressure block
    J[15] = (1.0 / biot_modulus) * a[0];  // ∂F_p/∂p_t
    
    // Coupling terms would be added for full Jacobian
}

void PoroelastodynamicsKernel::setMaterialProperties(double E, double nu, double rho, double phi) {
    youngs_modulus = E;
    poisson_ratio = nu;
    density_solid = rho;
    porosity = phi;
}

void PoroelastodynamicsKernel::setFluidProperties(double rho_f, double mu, double K_f) {
    density_fluid = rho_f;
    viscosity = mu;
    bulk_modulus_fluid = K_f;
}

void PoroelastodynamicsKernel::setBiotParameters(double alpha, double M) {
    biot_coefficient = alpha;
    biot_modulus = M;
}

void PoroelastodynamicsKernel::setWaveProperties(double vp_fast, double vs, double vp_slow) {
    p_wave_fast = vp_fast;
    s_wave = vs;
    p_wave_slow = vp_slow;
}

void PoroelastodynamicsKernel::setDamping(double alpha, double beta) {
    damping_alpha = alpha;
    damping_beta = beta;
}

void PoroelastodynamicsKernel::setStaticTriggeringMode(bool enable, double threshold, double duration) {
    enable_static_triggering = enable;
    trigger_threshold = threshold;
    event_duration = duration;
}

void PoroelastodynamicsKernel::enableDynamicPermeabilityChange(bool enable, double strain_coeff,
                                                               double stress_coeff, double recovery_time) {
    enable_dynamic_k = enable;
    k_strain_coeff = strain_coeff;
    k_stress_coeff = stress_coeff;
    k_recovery_time = recovery_time;
}

double PoroelastodynamicsKernel::computePermeabilityChange(const PetscScalar u_x[], 
                                                           const PetscScalar stress[],
                                                           double k_initial, double dt) {
    if (!enable_dynamic_k) {
        return k_initial;
    }
    
    // Compute volumetric strain
    double eps_vol = u_x[0] + u_x[4] + u_x[8];
    
    // Compute deviatoric strain
    double eps_mean = eps_vol / 3.0;
    double eps_dev_xx = u_x[0] - eps_mean;
    double eps_dev_yy = u_x[4] - eps_mean;
    double eps_dev_zz = u_x[8] - eps_mean;
    double eps_dev_mag = std::sqrt(eps_dev_xx*eps_dev_xx + eps_dev_yy*eps_dev_yy + 
                                   eps_dev_zz*eps_dev_zz);
    
    // Compute mean stress
    double stress_mean = (stress[0] + stress[1] + stress[2]) / 3.0;
    
    // Permeability change models:
    // 1) Volumetric strain model (pore opening/closing)
    double k_from_vol_strain = k_initial * std::exp(k_strain_coeff * eps_vol);
    
    // 2) Shear dilation model
    double k_from_shear = k_initial * (1.0 + 0.1 * eps_dev_mag);
    
    // 3) Effective stress model
    double k_from_stress = k_initial * std::exp(-k_stress_coeff * stress_mean);
    
    // Combined model (weighted average or max)
    double k_new = std::max({k_from_vol_strain, k_from_shear, k_from_stress});
    
    // Time-dependent recovery toward initial permeability
    if (k_recovery_time > 0.0) {
        double recovery_factor = std::exp(-dt / k_recovery_time);
        k_new = k_initial + (k_new - k_initial) * recovery_factor;
    }
    
    // Clamp to physical bounds
    k_new = std::max(1e-18, std::min(k_new, 1e-10));  // m²
    
    return k_new;
}

bool PoroelastodynamicsKernel::checkStaticTrigger(const PetscScalar stress[], double current_time) {
    if (!enable_static_triggering || event_active) {
        return false;
    }
    
    // Use Coulomb failure criterion
    double s_xx = stress[0];
    double s_yy = stress[1];
    double s_zz = stress[2];
    double s_xy = stress[3];
    double s_xz = stress[4];
    double s_yz = stress[5];
    
    // Maximum shear stress
    double tau_max = std::sqrt(s_xy*s_xy + s_xz*s_xz + s_yz*s_yz);
    
    if (tau_max >= trigger_threshold) {
        trigger_time = current_time;
        event_active = true;
        return true;
    }
    
    return false;
}

bool PoroelastodynamicsKernel::isInDynamicEvent(double current_time) const {
    if (!enable_static_triggering || !event_active) {
        return false;
    }
    
    return (current_time - trigger_time) <= event_duration;
}

// ============================================================================
// Dynamic Permeability Model
// ============================================================================

DynamicPermeabilityModel::DynamicPermeabilityModel()
    : k_initial(100e-15), strain_coefficient(1e-7), stress_coefficient(1e-15),
      recovery_tau(100.0), k_minimum(1e-18), k_maximum(1e-10) {}

void DynamicPermeabilityModel::setParameters(double k0, double strain_coeff, double stress_coeff,
                                             double recovery_time, double k_min, double k_max) {
    k_initial = k0;
    strain_coefficient = strain_coeff;
    stress_coefficient = stress_coeff;
    recovery_tau = recovery_time;
    k_minimum = k_min;
    k_maximum = k_max;
}

double DynamicPermeabilityModel::computeInstantaneousChange(double strain_amplitude, 
                                                            double stress_amplitude) {
    // Instantaneous permeability change during wave passage
    // Based on empirical observations from earthquake studies
    
    // Strain-induced change (dominates at high frequency)
    double delta_k_strain = strain_coefficient * strain_amplitude;
    
    // Stress-induced change (static component)
    double delta_k_stress = stress_coefficient * stress_amplitude;
    
    // Combined effect (nonlinear)
    double k_new = k_initial * (1.0 + delta_k_strain + delta_k_stress);
    
    return std::max(k_minimum, std::min(k_new, k_maximum));
}

double DynamicPermeabilityModel::computeTimeEvolution(double k_current, double k0, double dt) {
    // Exponential recovery to initial permeability
    if (recovery_tau <= 0.0) {
        return k_current;
    }
    
    double k_new = k0 + (k_current - k0) * std::exp(-dt / recovery_tau);
    return k_new;
}

double DynamicPermeabilityModel::strainModel(double volumetric_strain) {
    // Permeability change from volumetric strain
    // k/k0 = exp(β * ε_vol)
    // where β ~ 100-1000 for rocks
    
    double beta = 100.0;  // Strain sensitivity parameter
    double k_ratio = std::exp(beta * volumetric_strain);
    
    return k_initial * k_ratio;
}

double DynamicPermeabilityModel::stressModel(double effective_stress) {
    // Permeability reduction with effective stress
    // k/k0 = exp(-γ * σ')
    
    double gamma = 1e-8;  // Stress sensitivity (1/Pa)
    double k_ratio = std::exp(-gamma * effective_stress);
    
    return k_initial * k_ratio;
}

double DynamicPermeabilityModel::shearDilationModel(double shear_strain) {
    // Permeability increase from shear dilation
    // Often observed in fault zones during seismic slip
    
    double dilation_coefficient = 10.0;  // Enhancement factor
    double k_enhanced = k_initial * (1.0 + dilation_coefficient * std::abs(shear_strain));
    
    return std::min(k_enhanced, k_maximum);
}

double DynamicPermeabilityModel::crackOpeningModel(double normal_stress) {
    // Permeability change from crack opening/closure
    // Parallel plate model: k ~ w³ where w is aperture
    
    double stress_ref = 10e6;  // Reference stress (10 MPa)
    double aperture_ratio = std::pow(1.0 - normal_stress / stress_ref, 0.33);
    
    if (aperture_ratio < 0.1) aperture_ratio = 0.1;  // Minimum aperture
    
    double k_ratio = std::pow(aperture_ratio, 3.0);  // Cubic law
    
    return k_initial * k_ratio;
}

} // namespace ResSim
