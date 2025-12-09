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

namespace FSRM {

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
    (void)dm;
    (void)fe;
    PetscFunctionBeginUser;
    // Setup finite element spaces for dynamic problem
    PetscFunctionReturn(0);
}

void ElastodynamicsKernel::residual(const PetscScalar u[], const PetscScalar u_t[],
                                    const PetscScalar u_x[], const PetscScalar a[],
                                    const PetscReal x[], PetscScalar f[]) {
    // Suppress unused parameter warnings - these are part of the standard kernel interface
    (void)u;
    (void)a;
    (void)x;
    
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
    
    // Suppress unused variable warnings - stress components used in f1 flux term
    (void)sigma_xx; (void)sigma_yy; (void)sigma_zz;
    (void)sigma_xy; (void)sigma_xz; (void)sigma_yz;
    
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
    // Suppress unused parameter warnings - these are part of the standard kernel interface
    (void)u;
    (void)u_t;
    (void)u_x;
    (void)x;
    
    // Full Jacobian for elastodynamics equation:
    // ρ * ∂²u/∂t² = ∇·σ - C * ∂u/∂t
    //
    // With Rayleigh damping C = α*M + β*K where:
    // M = mass matrix (ρ * I)
    // K = stiffness matrix (∂σ/∂ε)
    //
    // The Jacobian has two main parts:
    // 1. g0: Jacobian w.r.t. solution u and time derivative ∂u/∂t
    //    - Mass: ρ * I * shift (where shift = a[0] from time integrator)
    //    - Mass damping: α * ρ * I * shift
    // 2. g3: Jacobian w.r.t. displacement gradient ∇u
    //    - Stiffness: elasticity tensor C_ijkl
    //    - Stiffness damping: β * C_ijkl * shift
    
    double shift = a[0];  // Time integration shift parameter
    
    // Lamé parameters
    double lambda = youngs_modulus * poisson_ratio / 
                   ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio));
    double mu = youngs_modulus / (2.0 * (1.0 + poisson_ratio));
    
    // Initialize full Jacobian to zero
    // For 3D elastodynamics: 3 components × 3 gradients × 3 directions = 81 entries for g3
    // Plus 9 entries for g0 (3×3 mass matrix)
    for (int i = 0; i < 81; ++i) J[i] = 0.0;
    
    // ==================================
    // g0: Mass matrix Jacobian
    // ==================================
    // ∂(ρ*∂²u/∂t²)/∂(∂u/∂t) = ρ * I * shift
    // With mass-proportional damping: (1 + α) * ρ * I * shift
    
    double mass_term = density * (1.0 + damping_alpha) * shift;
    
    // Diagonal blocks for mass matrix (3x3 identity scaled by mass_term)
    J[0] = mass_term;   // ∂F_x/∂(∂u_x/∂t)
    J[4] = mass_term;   // ∂F_y/∂(∂u_y/∂t)
    J[8] = mass_term;   // ∂F_z/∂(∂u_z/∂t)
    
    // ==================================
    // g3: Stiffness matrix Jacobian
    // ==================================
    // The elasticity tensor for isotropic material:
    // C_ijkl = λ * δ_ij * δ_kl + μ * (δ_ik * δ_jl + δ_il * δ_jk)
    //
    // Stress-strain relation: σ_ij = C_ijkl * ε_kl
    // where ε_kl = 0.5 * (∂u_k/∂x_l + ∂u_l/∂x_k)
    //
    // For gradient-based formulation in PETSc, we need:
    // ∂σ_ij/∂(∂u_m/∂x_n) = C_ijmn + C_ijnm (due to strain symmetry)
    //
    // With stiffness-proportional damping: (1 + β*shift) * C
    
    double stiffness_factor = 1.0 + damping_beta * shift;
    
    // Effective elastic constants
    double lambda_eff = lambda * stiffness_factor;
    double mu_eff = mu * stiffness_factor;
    
    // Layout: J[i*9 + j] for i,j ∈ {0..8} where
    // i = stress component index, j = gradient component index
    // Stress components: σ_xx(0), σ_yy(1), σ_zz(2), σ_xy(3), σ_xz(4), σ_yz(5), σ_yx(6), σ_zx(7), σ_zy(8)
    // Gradient components: ∂u_x/∂x(0), ∂u_x/∂y(1), ∂u_x/∂z(2), ∂u_y/∂x(3), ...
    
    // Actually use 9x9 for stress-gradient Jacobian starting at J[9]
    
    // σ_xx row (stress index 0)
    J[9 + 0*9 + 0] = lambda_eff + 2.0*mu_eff;  // ∂σ_xx/∂(∂u_x/∂x)
    J[9 + 0*9 + 4] = lambda_eff;                // ∂σ_xx/∂(∂u_y/∂y)
    J[9 + 0*9 + 8] = lambda_eff;                // ∂σ_xx/∂(∂u_z/∂z)
    
    // σ_yy row (stress index 1)
    J[9 + 1*9 + 0] = lambda_eff;                // ∂σ_yy/∂(∂u_x/∂x)
    J[9 + 1*9 + 4] = lambda_eff + 2.0*mu_eff;  // ∂σ_yy/∂(∂u_y/∂y)
    J[9 + 1*9 + 8] = lambda_eff;                // ∂σ_yy/∂(∂u_z/∂z)
    
    // σ_zz row (stress index 2)
    J[9 + 2*9 + 0] = lambda_eff;                // ∂σ_zz/∂(∂u_x/∂x)
    J[9 + 2*9 + 4] = lambda_eff;                // ∂σ_zz/∂(∂u_y/∂y)
    J[9 + 2*9 + 8] = lambda_eff + 2.0*mu_eff;  // ∂σ_zz/∂(∂u_z/∂z)
    
    // σ_xy row (stress index 3): σ_xy = 2*μ * ε_xy = μ*(∂u_x/∂y + ∂u_y/∂x)
    J[9 + 3*9 + 1] = mu_eff;  // ∂σ_xy/∂(∂u_x/∂y)
    J[9 + 3*9 + 3] = mu_eff;  // ∂σ_xy/∂(∂u_y/∂x)
    
    // σ_xz row (stress index 4): σ_xz = 2*μ * ε_xz = μ*(∂u_x/∂z + ∂u_z/∂x)
    J[9 + 4*9 + 2] = mu_eff;  // ∂σ_xz/∂(∂u_x/∂z)
    J[9 + 4*9 + 6] = mu_eff;  // ∂σ_xz/∂(∂u_z/∂x)
    
    // σ_yz row (stress index 5): σ_yz = 2*μ * ε_yz = μ*(∂u_y/∂z + ∂u_z/∂y)
    J[9 + 5*9 + 5] = mu_eff;  // ∂σ_yz/∂(∂u_y/∂z)
    J[9 + 5*9 + 7] = mu_eff;  // ∂σ_yz/∂(∂u_z/∂y)
    
    // Symmetric shear components (σ_yx = σ_xy, etc.)
    J[9 + 6*9 + 1] = mu_eff;  // ∂σ_yx/∂(∂u_x/∂y)
    J[9 + 6*9 + 3] = mu_eff;  // ∂σ_yx/∂(∂u_y/∂x)
    
    J[9 + 7*9 + 2] = mu_eff;  // ∂σ_zx/∂(∂u_x/∂z)
    J[9 + 7*9 + 6] = mu_eff;  // ∂σ_zx/∂(∂u_z/∂x)
    
    J[9 + 8*9 + 5] = mu_eff;  // ∂σ_zy/∂(∂u_y/∂z)
    J[9 + 8*9 + 7] = mu_eff;  // ∂σ_zy/∂(∂u_z/∂y)
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
    (void)dm;
    (void)fe;
    PetscFunctionBeginUser;
    // Setup finite element spaces for coupled problem
    PetscFunctionReturn(0);
}

void PoroelastodynamicsKernel::residual(const PetscScalar u[], const PetscScalar u_t[],
                                        const PetscScalar u_x[], const PetscScalar a[],
                                        const PetscReal x[], PetscScalar f[]) {
    // Suppress unused parameter warnings - these are part of the standard kernel interface
    (void)a;
    (void)x;
    
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
    (void)ux; (void)uy; (void)uz;  // Used for boundary conditions in full implementation
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
    (void)dp_dx; (void)dp_dy; (void)dp_dz;  // Used in f1 flux term
    
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
    (void)mobility;  // Used in f1 flux term
    
    // Add flux divergence (handled in f1 term in actual implementation)
    // f[3] += mobility * (dp_dx² + dp_dy² + dp_dz²)
}

void PoroelastodynamicsKernel::jacobian(const PetscScalar u[], const PetscScalar u_t[],
                                        const PetscScalar u_x[], const PetscScalar a[],
                                        const PetscReal x[], PetscScalar J[]) {
    // Suppress unused parameter warnings - these are part of the standard kernel interface
    (void)u;
    (void)u_t;
    (void)u_x;
    (void)x;
    
    // Full Jacobian for Biot poroelastodynamics
    // System of equations (4 fields: ux, uy, uz, p):
    //
    // 1) Momentum balance (3 equations):
    //    ρ_bulk * ∂²u/∂t² = ∇·σ_eff - α*∇p + damping
    //    where σ_eff = λ*tr(ε)*I + 2*μ*ε (drained effective stress)
    //
    // 2) Mass conservation (1 equation):
    //    (1/M)*∂p/∂t + α*∂(∇·u)/∂t + ∇·q = 0
    //    where q = -(k/μ)*∇p (Darcy flux)
    //
    // Full Jacobian is organized as a 4×4 block matrix:
    // [ J_uu  J_up ]   [ 3×3   3×1 ]
    // [ J_pu  J_pp ]   [ 1×3   1×1 ]
    //
    // With gradient terms (g3), total size is much larger
    
    double shift = a[0];  // Time integration shift parameter
    
    // Material parameters
    double lambda = youngs_modulus * poisson_ratio / 
                   ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio));
    double mu = youngs_modulus / (2.0 * (1.0 + poisson_ratio));
    double rho_bulk = (1.0 - porosity) * density_solid + porosity * density_fluid;
    
    // Initialize 16×16 Jacobian to zero (4 fields × 4 gradients for g0/g2)
    // Plus gradient-gradient terms for g3
    // Total: g0(4×4) + g1(4×4×3) + g2(4×4×3) + g3(4×3×4×3) - simplified to essential terms
    // Using 196 entries for 4-field × (4-field × 12-gradients) compact storage
    for (int i = 0; i < 196; ++i) J[i] = 0.0;
    
    // ==================================
    // g0: Jacobian w.r.t. solution fields
    // ==================================
    // Layout: J[0..15] for 4×4 solution Jacobian
    
    // Displacement-displacement blocks (3×3) - mass matrix with damping
    double mass_term = rho_bulk * (1.0 + damping_alpha) * shift;
    
    J[0] = mass_term;    // ∂F_ux/∂(∂u_x/∂t)
    J[5] = mass_term;    // ∂F_uy/∂(∂u_y/∂t)
    J[10] = mass_term;   // ∂F_uz/∂(∂u_z/∂t)
    
    // Displacement-pressure coupling (3×1) - Biot coupling
    // From momentum equation: -α*∇p contributes to ∂F_u/∂p
    // In weak form: ∫ α*p*∇·δu dV → Jacobian has α on u-p coupling
    J[3] = -biot_coefficient;   // ∂F_ux/∂p (from α*∂p/∂x term)
    J[7] = -biot_coefficient;   // ∂F_uy/∂p (from α*∂p/∂y term)
    J[11] = -biot_coefficient;  // ∂F_uz/∂p (from α*∂p/∂z term)
    
    // Pressure-displacement coupling (1×3) - Biot storage coupling
    // From mass conservation: α*∂(∇·u)/∂t
    // ∂F_p/∂(∂u/∂t) = α * ∇·(δu) → contributes through gradient terms
    J[12] = biot_coefficient * shift;   // ∂F_p/∂(∂u_x/∂t) - coupled through ∂u_x/∂x
    J[13] = biot_coefficient * shift;   // ∂F_p/∂(∂u_y/∂t) - coupled through ∂u_y/∂y
    J[14] = biot_coefficient * shift;   // ∂F_p/∂(∂u_z/∂t) - coupled through ∂u_z/∂z
    
    // Pressure-pressure block (1×1) - storage term
    // (1/M)*∂p/∂t → Jacobian = (1/M)*shift
    J[15] = (1.0 / biot_modulus) * shift;
    
    // ==================================
    // g3: Jacobian w.r.t. gradients
    // ==================================
    // Starting at J[16] for gradient Jacobian
    
    // Displacement gradient Jacobian (stiffness matrix)
    // Same structure as elastodynamics, but with undrained modulus for fast waves
    
    double stiffness_factor = 1.0 + damping_beta * shift;
    double lambda_eff = lambda * stiffness_factor;
    double mu_eff = mu * stiffness_factor;
    
    // Undrained bulk modulus adjustment for poroelastic effects
    // K_undrained = K_drained + α²*M
    double K_drained = lambda + 2.0*mu/3.0;
    double K_undrained = K_drained + biot_coefficient * biot_coefficient * biot_modulus;
    double lambda_undrained = K_undrained - 2.0*mu/3.0;
    
    // Use undrained modulus for dynamic response (fast waves)
    lambda_eff = lambda_undrained * stiffness_factor;
    
    // Elasticity tensor - displacement gradient to stress
    // σ_xx, σ_yy, σ_zz rows
    J[16 + 0*12 + 0] = lambda_eff + 2.0*mu_eff;  // ∂σ_xx/∂(∂u_x/∂x)
    J[16 + 0*12 + 4] = lambda_eff;                // ∂σ_xx/∂(∂u_y/∂y)
    J[16 + 0*12 + 8] = lambda_eff;                // ∂σ_xx/∂(∂u_z/∂z)
    
    J[16 + 1*12 + 0] = lambda_eff;                // ∂σ_yy/∂(∂u_x/∂x)
    J[16 + 1*12 + 4] = lambda_eff + 2.0*mu_eff;  // ∂σ_yy/∂(∂u_y/∂y)
    J[16 + 1*12 + 8] = lambda_eff;                // ∂σ_yy/∂(∂u_z/∂z)
    
    J[16 + 2*12 + 0] = lambda_eff;                // ∂σ_zz/∂(∂u_x/∂x)
    J[16 + 2*12 + 4] = lambda_eff;                // ∂σ_zz/∂(∂u_y/∂y)
    J[16 + 2*12 + 8] = lambda_eff + 2.0*mu_eff;  // ∂σ_zz/∂(∂u_z/∂z)
    
    // Shear stress rows (off-diagonal strain components)
    // σ_xy = 2*μ*ε_xy = μ*(∂u_x/∂y + ∂u_y/∂x)
    J[16 + 0*12 + 1] = mu_eff;  // Contribution to F_x from ∂u_x/∂y
    J[16 + 0*12 + 3] = mu_eff;  // Contribution to F_x from ∂u_y/∂x
    J[16 + 1*12 + 1] = mu_eff;  // Contribution to F_y from ∂u_x/∂y
    J[16 + 1*12 + 3] = mu_eff;  // Contribution to F_y from ∂u_y/∂x
    
    // σ_xz = 2*μ*ε_xz = μ*(∂u_x/∂z + ∂u_z/∂x)
    J[16 + 0*12 + 2] = mu_eff;  // Contribution to F_x from ∂u_x/∂z
    J[16 + 0*12 + 6] = mu_eff;  // Contribution to F_x from ∂u_z/∂x
    J[16 + 2*12 + 2] = mu_eff;  // Contribution to F_z from ∂u_x/∂z
    J[16 + 2*12 + 6] = mu_eff;  // Contribution to F_z from ∂u_z/∂x
    
    // σ_yz = 2*μ*ε_yz = μ*(∂u_y/∂z + ∂u_z/∂y)
    J[16 + 1*12 + 5] = mu_eff;  // Contribution to F_y from ∂u_y/∂z
    J[16 + 1*12 + 7] = mu_eff;  // Contribution to F_y from ∂u_z/∂y
    J[16 + 2*12 + 5] = mu_eff;  // Contribution to F_z from ∂u_y/∂z
    J[16 + 2*12 + 7] = mu_eff;  // Contribution to F_z from ∂u_z/∂y
    
    // Pressure equation - Darcy flux Jacobian
    // q = -(k/μ)*∇p, so ∂q/∂(∇p) = -(k/μ)*I
    // In weak form: ∫ q·∇v dV → gradient-gradient Jacobian
    
    double mobility = permeability / viscosity;
    
    // Dynamic permeability correction (if enabled)
    if (enable_dynamic_k) {
        // Simplified: just use base permeability for Jacobian
        // Full implementation would include ∂k/∂ε terms
        mobility = permeability_initial / viscosity;
    }
    
    // Pressure gradient Jacobian (row 3, columns 9-11 for ∂p/∂x, ∂p/∂y, ∂p/∂z)
    J[16 + 3*12 + 9] = mobility;   // ∂q_x/∂(∂p/∂x)
    J[16 + 3*12 + 10] = mobility;  // ∂q_y/∂(∂p/∂y)
    J[16 + 3*12 + 11] = mobility;  // ∂q_z/∂(∂p/∂z)
    
    // Cross-coupling: pressure gradient affects momentum (Terzaghi effective stress)
    // ∂(σ - α*p*I)/∂p = -α*I → contributes to g2 (not g3)
    // This is handled in the g0 block above
    
    // Cross-coupling: displacement divergence affects pressure equation
    // ∂/∂t(α*∇·u) → contributes α*shift to gradient terms
    // These are storage coupling terms already included in g0
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
    // double s_yy = stress[1];  // Not used in 2D failure criterion
    // double s_zz = stress[2];  // Not used in 2D failure criterion
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

} // namespace FSRM
