#ifdef USE_CUDA

#include "GPUKernels.cuh"
#include <cmath>
#include <algorithm>

namespace FSRM {
namespace GPU {

// ============================================================================
// Single-Phase Flow GPU Kernels Implementation
// ============================================================================

__global__ void computeAccumulation(
    const double* pressure_old,
    const double* pressure_new,
    const double* porosity,
    const double* compressibility,
    double dt,
    double* accumulation,
    int n_cells
) {
    int idx = getGlobalThreadId1D();
    if (idx >= n_cells) return;
    
    double dP = (pressure_new[idx] - pressure_old[idx]) / dt;
    accumulation[idx] = porosity[idx] * compressibility[idx] * dP;
}

__global__ void computeDarcyFlux(
    const double* pressure,
    const double* permeability,
    const double* viscosity,
    const int* connectivity,
    double* flux,
    int n_faces,
    int dim
) {
    int face_idx = getGlobalThreadId1D();
    if (face_idx >= n_faces) return;
    
    // Get neighboring cells
    int cell_left = connectivity[face_idx * 2];
    int cell_right = connectivity[face_idx * 2 + 1];
    
    if (cell_left < 0 || cell_right < 0) {
        flux[face_idx] = 0.0;
        return;
    }
    
    // Pressure gradient
    double dP = pressure[cell_right] - pressure[cell_left];
    
    // Harmonic average of permeability
    double k_left = permeability[cell_left];
    double k_right = permeability[cell_right];
    double k_face = 2.0 * k_left * k_right / (k_left + k_right + 1e-30);
    
    // Arithmetic average of viscosity
    double mu = 0.5 * (viscosity[cell_left] + viscosity[cell_right]);
    
    // Darcy flux: q = -k/mu * dP/dx
    flux[face_idx] = -k_face / mu * dP;
}

__global__ void updatePressure(
    double* pressure,
    const double* accumulation,
    const double* flux_divergence,
    double dt,
    int n_cells
) {
    int idx = getGlobalThreadId1D();
    if (idx >= n_cells) return;
    
    // Explicit update: P_new = P_old + dt * (flux_div - accumulation)
    pressure[idx] += dt * (flux_divergence[idx] - accumulation[idx]);
}

// ============================================================================
// Elastodynamics GPU Kernels Implementation
// ============================================================================

__global__ void computeStrain(
    const double* displacement,
    double* strain,
    const int* connectivity,
    int n_cells,
    int dim
) {
    int idx = getGlobalThreadId1D();
    if (idx >= n_cells) return;
    
    // Simple finite difference approximation for strain
    // strain_xx = du_x/dx, strain_yy = du_y/dy, strain_zz = du_z/dz
    // strain_xy = 0.5 * (du_x/dy + du_y/dx), etc.
    
    // For structured grid, compute gradients
    int i = idx % connectivity[0];
    int j = (idx / connectivity[0]) % connectivity[1];
    int k = idx / (connectivity[0] * connectivity[1]);
    
    double dx = 1.0; // Grid spacing (should be passed as parameter)
    
    if (dim >= 1) {
        // x-direction strain
        if (i < connectivity[0] - 1) {
            int idx_next = (k * connectivity[1] + j) * connectivity[0] + (i + 1);
            strain[idx * 6 + 0] = (displacement[idx_next * 3 + 0] - displacement[idx * 3 + 0]) / dx;
        }
    }
    
    if (dim >= 2) {
        // y-direction strain
        if (j < connectivity[1] - 1) {
            int idx_next = (k * connectivity[1] + (j + 1)) * connectivity[0] + i;
            strain[idx * 6 + 1] = (displacement[idx_next * 3 + 1] - displacement[idx * 3 + 1]) / dx;
        }
    }
    
    if (dim >= 3) {
        // z-direction strain
        if (k < connectivity[2] - 1) {
            int idx_next = ((k + 1) * connectivity[1] + j) * connectivity[0] + i;
            strain[idx * 6 + 2] = (displacement[idx_next * 3 + 2] - displacement[idx * 3 + 2]) / dx;
        }
    }
    
    // Shear strains (simplified)
    strain[idx * 6 + 3] = 0.0; // xy
    strain[idx * 6 + 4] = 0.0; // xz
    strain[idx * 6 + 5] = 0.0; // yz
}

__global__ void computeStress(
    const double* strain,
    const double* youngs_modulus,
    const double* poisson_ratio,
    double* stress,
    int n_cells,
    int dim
) {
    int idx = getGlobalThreadId1D();
    if (idx >= n_cells) return;
    
    double E = youngs_modulus[idx];
    double nu = poisson_ratio[idx];
    
    // Lame parameters
    double lambda, mu;
    computeLameParameters(E, nu, lambda, mu);
    
    // Strain components
    double eps_xx = strain[idx * 6 + 0];
    double eps_yy = strain[idx * 6 + 1];
    double eps_zz = strain[idx * 6 + 2];
    double eps_xy = strain[idx * 6 + 3];
    double eps_xz = strain[idx * 6 + 4];
    double eps_yz = strain[idx * 6 + 5];
    
    // Volumetric strain
    double trace = eps_xx + eps_yy + eps_zz;
    
    // Stress tensor (Hooke's law)
    stress[idx * 6 + 0] = lambda * trace + 2.0 * mu * eps_xx; // sigma_xx
    stress[idx * 6 + 1] = lambda * trace + 2.0 * mu * eps_yy; // sigma_yy
    stress[idx * 6 + 2] = lambda * trace + 2.0 * mu * eps_zz; // sigma_zz
    stress[idx * 6 + 3] = 2.0 * mu * eps_xy; // sigma_xy
    stress[idx * 6 + 4] = 2.0 * mu * eps_xz; // sigma_xz
    stress[idx * 6 + 5] = 2.0 * mu * eps_yz; // sigma_yz
}

__global__ void computeStressDivergence(
    const double* stress,
    const int* connectivity,
    double* stress_div,
    int n_cells,
    int dim
) {
    int idx = getGlobalThreadId1D();
    if (idx >= n_cells) return;
    
    // Compute divergence of stress tensor
    // div(sigma) = [d_sigma_xx/dx + d_sigma_xy/dy + d_sigma_xz/dz, ...]
    
    double dx = 1.0; // Grid spacing
    
    int i = idx % connectivity[0];
    int j = (idx / connectivity[0]) % connectivity[1];
    int k = idx / (connectivity[0] * connectivity[1]);
    
    // X-component
    double div_x = 0.0;
    if (i > 0 && i < connectivity[0] - 1) {
        int idx_prev = (k * connectivity[1] + j) * connectivity[0] + (i - 1);
        int idx_next = (k * connectivity[1] + j) * connectivity[0] + (i + 1);
        div_x += (stress[idx_next * 6 + 0] - stress[idx_prev * 6 + 0]) / (2.0 * dx);
    }
    stress_div[idx * 3 + 0] = div_x;
    
    // Y-component
    double div_y = 0.0;
    if (j > 0 && j < connectivity[1] - 1) {
        int idx_prev = (k * connectivity[1] + (j - 1)) * connectivity[0] + i;
        int idx_next = (k * connectivity[1] + (j + 1)) * connectivity[0] + i;
        div_y += (stress[idx_next * 6 + 1] - stress[idx_prev * 6 + 1]) / (2.0 * dx);
    }
    stress_div[idx * 3 + 1] = div_y;
    
    // Z-component
    double div_z = 0.0;
    if (k > 0 && k < connectivity[2] - 1) {
        int idx_prev = ((k - 1) * connectivity[1] + j) * connectivity[0] + i;
        int idx_next = ((k + 1) * connectivity[1] + j) * connectivity[0] + i;
        div_z += (stress[idx_next * 6 + 2] - stress[idx_prev * 6 + 2]) / (2.0 * dx);
    }
    stress_div[idx * 3 + 2] = div_z;
}

__global__ void updateVelocity(
    double* velocity,
    const double* stress_div,
    const double* density,
    double dt,
    int n_cells,
    int dim
) {
    int idx = getGlobalThreadId1D();
    if (idx >= n_cells) return;
    
    double rho = density[idx];
    
    // v_new = v_old + dt * (div(sigma) / rho)
    for (int d = 0; d < dim; ++d) {
        velocity[idx * 3 + d] += dt * stress_div[idx * 3 + d] / rho;
    }
}

__global__ void updateDisplacement(
    double* displacement,
    const double* velocity,
    double dt,
    int n_cells,
    int dim
) {
    int idx = getGlobalThreadId1D();
    if (idx >= n_cells) return;
    
    // u_new = u_old + dt * v
    for (int d = 0; d < dim; ++d) {
        displacement[idx * 3 + d] += dt * velocity[idx * 3 + d];
    }
}

__global__ void applyRayleighDamping(
    double* force,
    const double* velocity,
    const double* displacement,
    const double* density,
    const double* stiffness,
    double alpha,
    double beta,
    int n_cells,
    int dim
) {
    int idx = getGlobalThreadId1D();
    if (idx >= n_cells) return;
    
    // Rayleigh damping: F_damp = -alpha*M*v - beta*K*u
    // M is mass matrix (diagonal: rho*V)
    // K is stiffness matrix
    
    double rho = density[idx];
    double k = stiffness[idx];
    
    for (int d = 0; d < dim; ++d) {
        double damp_force = -alpha * rho * velocity[idx * 3 + d] 
                          - beta * k * displacement[idx * 3 + d];
        force[idx * 3 + d] += damp_force;
    }
}

// ============================================================================
// Poroelastodynamics GPU Kernels Implementation
// ============================================================================

__global__ void computePoroelasticForce(
    const double* displacement,
    const double* pressure,
    const double* velocity_solid,
    const double* youngs_modulus,
    const double* poisson_ratio,
    const double* biot_coefficient,
    const double* porosity,
    const double* density_solid,
    const double* density_fluid,
    double* force_solid,
    const int* connectivity,
    int n_cells,
    int dim
) {
    int idx = getGlobalThreadId1D();
    if (idx >= n_cells) return;
    
    // Compute effective stress and poroelastic coupling
    double E = youngs_modulus[idx];
    double nu = poisson_ratio[idx];
    double alpha = biot_coefficient[idx];
    double phi = porosity[idx];
    double rho_s = density_solid[idx];
    double rho_f = density_fluid[idx];
    
    // Bulk density
    double rho = (1.0 - phi) * rho_s + phi * rho_f;
    
    // This is simplified - full implementation would compute:
    // 1. Effective stress: sigma' = sigma - alpha*p*I
    // 2. Divergence of effective stress
    // 3. Coupling with fluid pressure gradient
    
    // Placeholder: inertial term
    for (int d = 0; d < dim; ++d) {
        force_solid[idx * 3 + d] = rho * velocity_solid[idx * 3 + d];
    }
}

__global__ void computeFluidMassBalance(
    const double* pressure,
    const double* displacement,
    const double* permeability,
    const double* viscosity,
    const double* porosity,
    const double* biot_coefficient,
    const double* biot_modulus,
    double dt,
    double* pressure_rate,
    const int* connectivity,
    int n_cells
) {
    int idx = getGlobalThreadId1D();
    if (idx >= n_cells) return;
    
    // Fluid mass balance:
    // d/dt(alpha*div(u) + p/M) + div(q) = 0
    // where q = -(k/mu)*grad(p)
    
    double phi = porosity[idx];
    double alpha = biot_coefficient[idx];
    double M = biot_modulus[idx];
    double k = permeability[idx];
    double mu = viscosity[idx];
    
    // Storage term
    double storage = alpha * 0.0 + pressure[idx] / M; // div(u) term needs to be computed
    
    // Darcy flux divergence (simplified)
    double flux_div = 0.0;
    
    pressure_rate[idx] = -(flux_div + storage / dt);
}

__global__ void updateDynamicPermeability(
    double* permeability,
    const double* permeability_initial,
    const double* strain,
    const double* stress,
    double strain_coefficient,
    double stress_coefficient,
    double recovery_rate,
    double dt,
    double k_min,
    double k_max,
    int n_cells
) {
    int idx = getGlobalThreadId1D();
    if (idx >= n_cells) return;
    
    double k0 = permeability_initial[idx];
    
    // Volumetric strain
    double eps_vol = strain[idx * 6 + 0] + strain[idx * 6 + 1] + strain[idx * 6 + 2];
    
    // Mean stress
    double s_mean = (stress[idx * 6 + 0] + stress[idx * 6 + 1] + stress[idx * 6 + 2]) / 3.0;
    
    // Permeability change models
    double k_strain = k0 * exp(strain_coefficient * eps_vol);
    double k_stress = k0 * exp(-stress_coefficient * s_mean);
    
    // Combined model
    double k_new = max(k_strain, k_stress);
    
    // Time-dependent recovery
    if (recovery_rate > 0.0) {
        double recovery = exp(-dt * recovery_rate);
        k_new = k0 + (k_new - k0) * recovery;
    }
    
    // Clamp to bounds
    permeability[idx] = max(k_min, min(k_new, k_max));
}

__global__ void checkStaticTrigger(
    const double* stress,
    double threshold,
    bool* triggered,
    int* trigger_location,
    int n_cells
) {
    int idx = getGlobalThreadId1D();
    if (idx >= n_cells) return;
    
    // Compute von Mises stress
    double s_xx = stress[idx * 6 + 0];
    double s_yy = stress[idx * 6 + 1];
    double s_zz = stress[idx * 6 + 2];
    double s_xy = stress[idx * 6 + 3];
    double s_xz = stress[idx * 6 + 4];
    double s_yz = stress[idx * 6 + 5];
    
    double von_mises = computeVonMisesStress(s_xx, s_yy, s_zz, s_xy, s_xz, s_yz);
    
    if (von_mises >= threshold) {
        *triggered = true;
        atomicMin(trigger_location, idx); // Record first trigger location
    }
}

// ============================================================================
// Utility Kernels Implementation
// ============================================================================

__global__ void vectorAdd(const double* a, const double* b, double* c, int n) {
    int idx = getGlobalThreadId1D();
    if (idx >= n) return;
    c[idx] = a[idx] + b[idx];
}

__global__ void vectorScale(const double* a, double scale, double* b, int n) {
    int idx = getGlobalThreadId1D();
    if (idx >= n) return;
    b[idx] = scale * a[idx];
}

__global__ void vectorDot(const double* a, const double* b, double* result, int n) {
    // Shared memory for reduction
    __shared__ double sdata[256];
    
    int idx = getGlobalThreadId1D();
    int tid = threadIdx.x;
    
    // Load data into shared memory
    double val = 0.0;
    if (idx < n) {
        val = a[idx] * b[idx];
    }
    sdata[tid] = val;
    __syncthreads();
    
    // Reduction in shared memory
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            sdata[tid] += sdata[tid + s];
        }
        __syncthreads();
    }
    
    // Write result for this block
    if (tid == 0) {
        atomicAdd(result, sdata[0]);
    }
}

__global__ void computeWaveVelocities(
    const double* youngs_modulus,
    const double* poisson_ratio,
    const double* density,
    double* vp,
    double* vs,
    int n_cells
) {
    int idx = getGlobalThreadId1D();
    if (idx >= n_cells) return;
    
    double E = youngs_modulus[idx];
    double nu = poisson_ratio[idx];
    double rho = density[idx];
    
    // Lame parameters
    double lambda, mu;
    computeLameParameters(E, nu, lambda, mu);
    
    // Wave velocities
    vp[idx] = sqrt((lambda + 2.0 * mu) / rho); // P-wave
    vs[idx] = sqrt(mu / rho);                    // S-wave
}

} // namespace GPU
} // namespace FSRM

#endif // USE_CUDA
