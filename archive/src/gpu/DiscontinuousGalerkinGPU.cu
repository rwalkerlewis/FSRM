/**
 * @file DiscontinuousGalerkinGPU.cu
 * @brief CUDA implementations of GPU-accelerated DG wave equation kernels
 *
 * Implements the kernels declared in DiscontinuousGalerkinGPU.cuh for
 * elastic and poroelastic wave propagation using the ADER-DG method.
 *
 * Each kernel is parallelised over elements (volume/update kernels) or
 * faces (flux kernels). The DG method is inherently element-local, making
 * it an excellent fit for GPU execution.
 */

#ifdef USE_CUDA

#include "numerics/DiscontinuousGalerkinGPU.cuh"
#include <cmath>
#include <cstdio>

namespace FSRM {
namespace GPU {
namespace DG {

// ============================================================================
// Error Checking
// ============================================================================

#ifndef CUDA_CHECK
#define CUDA_CHECK(call) \
    do { \
        cudaError_t err = call; \
        if (err != cudaSuccess) { \
            fprintf(stderr, "CUDA error at %s:%d: %s\n", \
                    __FILE__, __LINE__, cudaGetErrorString(err)); \
        } \
    } while(0)
#endif

// ============================================================================
// Device Helper Functions
// ============================================================================

/**
 * @brief Compute elastic flux F_d(Q) in direction d for a single DOF
 *
 * For the velocity-stress formulation of the elastic wave equation:
 *
 * Q = (σ_xx, σ_yy, σ_zz, σ_xy, σ_yz, σ_xz, v_x, v_y, v_z)
 *
 * F_x(Q):           F_y(Q):           F_z(Q):
 *  -(λ+2μ) v_x       -λ v_y            -λ v_z
 *  -λ v_x            -(λ+2μ) v_y       -λ v_z
 *  -λ v_x            -λ v_y           -(λ+2μ) v_z
 *  -μ v_y            -μ v_x             0
 *   0                -μ v_z            -μ v_y
 *  -μ v_z             0                -μ v_x
 *  -σ_xx/ρ           -σ_xy/ρ           -σ_xz/ρ
 *  -σ_xy/ρ           -σ_yy/ρ           -σ_yz/ρ
 *  -σ_xz/ρ           -σ_yz/ρ           -σ_zz/ρ
 */
__device__ void elasticFlux(
    const double* Q_local,    // 9-component state vector at one DOF
    double lambda, double mu, double rho_inv,
    double* Fx, double* Fy, double* Fz
) {
    // Extract state
    double sxx = Q_local[0], syy = Q_local[1], szz = Q_local[2];
    double sxy = Q_local[3], syz = Q_local[4], sxz = Q_local[5];
    double vx  = Q_local[6], vy  = Q_local[7], vz  = Q_local[8];

    double lam2mu = lambda + 2.0 * mu;

    // F_x
    Fx[0] = -lam2mu * vx;
    Fx[1] = -lambda * vx;
    Fx[2] = -lambda * vx;
    Fx[3] = -mu * vy;
    Fx[4] = 0.0;
    Fx[5] = -mu * vz;
    Fx[6] = -sxx * rho_inv;
    Fx[7] = -sxy * rho_inv;
    Fx[8] = -sxz * rho_inv;

    // F_y
    Fy[0] = -lambda * vy;
    Fy[1] = -lam2mu * vy;
    Fy[2] = -lambda * vy;
    Fy[3] = -mu * vx;
    Fy[4] = -mu * vz;
    Fy[5] = 0.0;
    Fy[6] = -sxy * rho_inv;
    Fy[7] = -syy * rho_inv;
    Fy[8] = -syz * rho_inv;

    // F_z
    Fz[0] = -lambda * vz;
    Fz[1] = -lambda * vz;
    Fz[2] = -lam2mu * vz;
    Fz[3] = 0.0;
    Fz[4] = -mu * vy;
    Fz[5] = -mu * vx;
    Fz[6] = -sxz * rho_inv;
    Fz[7] = -syz * rho_inv;
    Fz[8] = -szz * rho_inv;
}

/**
 * @brief Compute P-wave and S-wave velocities from Lamé parameters
 */
__device__ void waveVelocities(double lambda, double mu, double rho,
                                double& Vp, double& Vs) {
    Vp = sqrt((lambda + 2.0 * mu) / rho);
    Vs = sqrt(mu / rho);
}

// ============================================================================
// Volume Integral — Elastic Wave Equation
// ============================================================================

__global__ void computeVolumeIntegralElastic(
    const double* __restrict__ Q,
    const ElementMaterial* __restrict__ material,
    const double* __restrict__ stiffness_x,
    const double* __restrict__ stiffness_y,
    const double* __restrict__ stiffness_z,
    const double* __restrict__ inv_jacobian,
    const double* __restrict__ det_jacobian,
    double* __restrict__ R_vol,
    int n_elem,
    int n_dof
) {
    int elem = blockIdx.x * blockDim.x + threadIdx.x;
    if (elem >= n_elem) return;

    const int nv = N_VAR_ELASTIC;  // 9

    // Material for this element
    double lambda = material[elem].lambda;
    double mu_val = material[elem].mu;
    double rho_inv = 1.0 / material[elem].rho;

    // Inverse Jacobian for coordinate transformation (reference → physical)
    const double* J_inv = &inv_jacobian[elem * 9];
    double det_J = det_jacobian[elem];

    // Pointer into element's solution and residual
    const double* Q_e = &Q[elem * n_dof * nv];
    double* R_e = &R_vol[elem * n_dof * nv];

    // Zero the residual for this element
    for (int i = 0; i < n_dof * nv; i++) {
        R_e[i] = 0.0;
    }

    // For each DOF, compute the physical flux, transform to reference coords,
    // and accumulate via stiffness matrices
    for (int j = 0; j < n_dof; j++) {
        // Local state at DOF j
        double Q_local[N_VAR_ELASTIC];
        for (int v = 0; v < nv; v++) {
            Q_local[v] = Q_e[j * nv + v];
        }

        // Compute physical fluxes
        double Fx[N_VAR_ELASTIC], Fy[N_VAR_ELASTIC], Fz[N_VAR_ELASTIC];
        elasticFlux(Q_local, lambda, mu_val, rho_inv, Fx, Fy, Fz);

        // Transform to reference coordinates: F̃_ξ = J⁻¹ · F_phys
        // F̃_ξ = J⁻¹[0,0]*Fx + J⁻¹[0,1]*Fy + J⁻¹[0,2]*Fz
        // F̃_η = J⁻¹[1,0]*Fx + J⁻¹[1,1]*Fy + J⁻¹[1,2]*Fz
        // F̃_ζ = J⁻¹[2,0]*Fx + J⁻¹[2,1]*Fy + J⁻¹[2,2]*Fz
        for (int v = 0; v < nv; v++) {
            double F_xi   = J_inv[0] * Fx[v] + J_inv[1] * Fy[v] + J_inv[2] * Fz[v];
            double F_eta  = J_inv[3] * Fx[v] + J_inv[4] * Fy[v] + J_inv[5] * Fz[v];
            double F_zeta = J_inv[6] * Fx[v] + J_inv[7] * Fy[v] + J_inv[8] * Fz[v];

            // Multiply by stiffness matrices and accumulate:
            // R[i,v] += K_ξ[i,j] * F̃_ξ[j,v] + K_η[i,j] * F̃_η[j,v] + K_ζ[i,j] * F̃_ζ[j,v]
            for (int i = 0; i < n_dof; i++) {
                R_e[i * nv + v] += det_J * (
                    stiffness_x[i * n_dof + j] * F_xi +
                    stiffness_y[i * n_dof + j] * F_eta +
                    stiffness_z[i * n_dof + j] * F_zeta
                );
            }
        }
    }
}

// ============================================================================
// Volume Integral — Poroelastic Wave Equation
// ============================================================================

__global__ void computeVolumeIntegralPoroelastic(
    const double* __restrict__ Q,
    const ElementBiotMaterial* __restrict__ material,
    const double* __restrict__ stiffness_x,
    const double* __restrict__ stiffness_y,
    const double* __restrict__ stiffness_z,
    const double* __restrict__ inv_jacobian,
    const double* __restrict__ det_jacobian,
    double* __restrict__ R_vol,
    int n_elem,
    int n_dof
) {
    int elem = blockIdx.x * blockDim.x + threadIdx.x;
    if (elem >= n_elem) return;

    const int nv = N_VAR_POROELASTIC;  // 13

    // Biot parameters
    const ElementBiotMaterial& mat = material[elem];
    double lambda_u = mat.lambda_u;
    double mu_val = mat.mu;
    double alpha = mat.alpha;
    double M = mat.M;
    double rho_bulk = (1.0 - mat.porosity) * mat.rho_solid + mat.porosity * mat.rho_fluid;
    double rho_inv = 1.0 / rho_bulk;
    double kappa = mat.kappa;

    const double* J_inv = &inv_jacobian[elem * 9];
    double det_J = det_jacobian[elem];

    const double* Q_e = &Q[elem * n_dof * nv];
    double* R_e = &R_vol[elem * n_dof * nv];

    // Zero residual
    for (int i = 0; i < n_dof * nv; i++) {
        R_e[i] = 0.0;
    }

    // Poroelastic flux: extends elastic flux with pore pressure coupling
    // Q = (σ'_xx, σ'_yy, σ'_zz, σ'_xy, σ'_yz, σ'_xz, v_x, v_y, v_z, p, q_x, q_y, q_z)
    // where σ' is effective stress, p is pore pressure, q is Darcy flux
    for (int j = 0; j < n_dof; j++) {
        double Q_local[N_VAR_POROELASTIC];
        for (int v = 0; v < nv; v++) {
            Q_local[v] = Q_e[j * nv + v];
        }

        // Effective stress components
        double sxx = Q_local[0], syy = Q_local[1], szz = Q_local[2];
        double sxy = Q_local[3], syz = Q_local[4], sxz = Q_local[5];
        double vx = Q_local[6], vy = Q_local[7], vz = Q_local[8];
        double p = Q_local[9];
        // q_x, q_y, q_z at indices 10, 11, 12

        double lam2mu = lambda_u + 2.0 * mu_val;

        // Solid momentum flux (includes pore pressure: total stress = σ' - α p I)
        double Fx[N_VAR_POROELASTIC] = {0};
        double Fy[N_VAR_POROELASTIC] = {0};
        double Fz[N_VAR_POROELASTIC] = {0};

        // Stress rate equations (undrained)
        Fx[0] = -lam2mu * vx;  Fy[0] = -lambda_u * vy;  Fz[0] = -lambda_u * vz;
        Fx[1] = -lambda_u * vx; Fy[1] = -lam2mu * vy;   Fz[1] = -lambda_u * vz;
        Fx[2] = -lambda_u * vx; Fy[2] = -lambda_u * vy;  Fz[2] = -lam2mu * vz;
        Fx[3] = -mu_val * vy;   Fy[3] = -mu_val * vx;    Fz[3] = 0.0;
        Fx[4] = 0.0;            Fy[4] = -mu_val * vz;     Fz[4] = -mu_val * vy;
        Fx[5] = -mu_val * vz;   Fy[5] = 0.0;              Fz[5] = -mu_val * vx;

        // Momentum equations: ρ ∂v/∂t = ∇·(σ' - αpI)
        Fx[6] = -(sxx - alpha * p) * rho_inv;
        Fy[6] = -sxy * rho_inv;
        Fz[6] = -sxz * rho_inv;

        Fx[7] = -sxy * rho_inv;
        Fy[7] = -(syy - alpha * p) * rho_inv;
        Fz[7] = -syz * rho_inv;

        Fx[8] = -sxz * rho_inv;
        Fy[8] = -syz * rho_inv;
        Fz[8] = -(szz - alpha * p) * rho_inv;

        // Pressure equation: (1/M) ∂p/∂t + α ∇·v + ∇·q = 0
        Fx[9] = -alpha * M * vx;
        Fy[9] = -alpha * M * vy;
        Fz[9] = -alpha * M * vz;

        // Darcy flux: q = -κ ∇p  (no explicit flux, handled via stiffness)
        Fx[10] = -kappa * p;  Fy[10] = 0.0;         Fz[10] = 0.0;
        Fx[11] = 0.0;         Fy[11] = -kappa * p;   Fz[11] = 0.0;
        Fx[12] = 0.0;         Fy[12] = 0.0;          Fz[12] = -kappa * p;

        // Transform and accumulate (same as elastic)
        for (int v = 0; v < nv; v++) {
            double F_xi   = J_inv[0] * Fx[v] + J_inv[1] * Fy[v] + J_inv[2] * Fz[v];
            double F_eta  = J_inv[3] * Fx[v] + J_inv[4] * Fy[v] + J_inv[5] * Fz[v];
            double F_zeta = J_inv[6] * Fx[v] + J_inv[7] * Fy[v] + J_inv[8] * Fz[v];

            for (int i = 0; i < n_dof; i++) {
                R_e[i * nv + v] += det_J * (
                    stiffness_x[i * n_dof + j] * F_xi +
                    stiffness_y[i * n_dof + j] * F_eta +
                    stiffness_z[i * n_dof + j] * F_zeta
                );
            }
        }
    }
}

// ============================================================================
// Numerical Flux — Rusanov / Lax-Friedrichs
// ============================================================================

__global__ void computeNumericalFluxRusanov(
    const double* __restrict__ Q,
    const ElementMaterial* __restrict__ material,
    const int* __restrict__ face_neighbors,
    const int* __restrict__ face_dof_map,
    const double* __restrict__ face_normals,
    const double* __restrict__ face_area,
    const double* __restrict__ mass_inv,
    double* __restrict__ R_surf,
    int n_interior_faces,
    int n_dof,
    int n_var
) {
    int face = blockIdx.x * blockDim.x + threadIdx.x;
    if (face >= n_interior_faces) return;

    // Get left and right element indices
    int elem_L = face_neighbors[face * 2];
    int elem_R = face_neighbors[face * 2 + 1];

    // Face normal (outward from left element)
    double nx = face_normals[face * 3];
    double ny = face_normals[face * 3 + 1];
    double nz = face_normals[face * 3 + 2];
    double area = face_area[face];

    // Material properties
    double Vp_L, Vs_L, Vp_R, Vs_R;
    waveVelocities(material[elem_L].lambda, material[elem_L].mu,
                    material[elem_L].rho, Vp_L, Vs_L);
    waveVelocities(material[elem_R].lambda, material[elem_R].mu,
                    material[elem_R].rho, Vp_R, Vs_R);

    // Maximum wave speed for Rusanov flux
    double max_speed = fmax(Vp_L, Vp_R);

    // For each face DOF (simplified: use first DOF on face)
    // In full implementation, loop over face quadrature points
    // For now, use element-averaged values
    const double* Q_L = &Q[elem_L * n_dof * n_var];
    const double* Q_R = &Q[elem_R * n_dof * n_var];
    double* R_L = &R_surf[elem_L * n_dof * n_var];
    double* R_R = &R_surf[elem_R * n_dof * n_var];

    // Compute average state at face (using first DOF as proxy)
    double rho_inv_L = 1.0 / material[elem_L].rho;
    double rho_inv_R = 1.0 / material[elem_R].rho;

    double Fx_L[N_VAR_ELASTIC], Fy_L[N_VAR_ELASTIC], Fz_L[N_VAR_ELASTIC];
    double Fx_R[N_VAR_ELASTIC], Fy_R[N_VAR_ELASTIC], Fz_R[N_VAR_ELASTIC];

    elasticFlux(Q_L, material[elem_L].lambda, material[elem_L].mu, rho_inv_L,
                Fx_L, Fy_L, Fz_L);
    elasticFlux(Q_R, material[elem_R].lambda, material[elem_R].mu, rho_inv_R,
                Fx_R, Fy_R, Fz_R);

    // Rusanov flux: F* = 0.5*(F_L + F_R)·n - 0.5*|λ_max|*(Q_R - Q_L)
    for (int v = 0; v < n_var; v++) {
        double Fn_L = Fx_L[v] * nx + Fy_L[v] * ny + Fz_L[v] * nz;
        double Fn_R = Fx_R[v] * nx + Fy_R[v] * ny + Fz_R[v] * nz;

        double flux = 0.5 * (Fn_L + Fn_R) - 0.5 * max_speed * (Q_R[v] - Q_L[v]);
        flux *= area;

        // Distribute to left element (subtract) and right element (add)
        // Using atomicAdd for thread safety (multiple faces share elements)
        atomicAdd(&R_L[v], -flux);
        atomicAdd(&R_R[v],  flux);
    }
}

// ============================================================================
// Free Surface Boundary Condition
// ============================================================================

__global__ void applyFreeSurfaceBC(
    double* __restrict__ Q,
    const int* __restrict__ boundary_faces,
    const int* __restrict__ boundary_face_elem,
    const int* __restrict__ face_dof_map,
    const double* __restrict__ face_normals,
    const double* __restrict__ face_area,
    const double* __restrict__ mass_inv,
    double* __restrict__ R_surf,
    int n_boundary_faces,
    int n_dof,
    int n_var
) {
    int bface = blockIdx.x * blockDim.x + threadIdx.x;
    if (bface >= n_boundary_faces) return;

    int elem = boundary_face_elem[bface];
    double nx = face_normals[bface * 3];
    double ny = face_normals[bface * 3 + 1];
    double nz = face_normals[bface * 3 + 2];
    double area = face_area[bface];

    // Mirror state: reflect velocity, negate normal traction
    // For free surface σ·n = 0, which means:
    //   σ⁺_xx = -σ⁻_xx (for component normal to surface)
    //   v⁺ = v⁻ (velocity unchanged)
    const double* Q_e = &Q[elem * n_dof * n_var];
    double* R_e = &R_surf[elem * n_dof * n_var];

    // Use first DOF (simplified)
    double Q_mirror[N_VAR_ELASTIC];
    for (int v = 0; v < n_var && v < N_VAR_ELASTIC; v++) {
        Q_mirror[v] = Q_e[v];
    }

    // Negate stress components (free surface: σ·n = 0 → mirror with sign flip)
    // For a z-normal surface: σ_xz → -σ_xz, σ_yz → -σ_yz, σ_zz → -σ_zz
    // General case: compute traction t = σ·n and subtract 2t
    double traction[3];
    traction[0] = Q_mirror[0]*nx + Q_mirror[3]*ny + Q_mirror[5]*nz;  // σ_xx*nx + σ_xy*ny + σ_xz*nz
    traction[1] = Q_mirror[3]*nx + Q_mirror[1]*ny + Q_mirror[4]*nz;  // σ_xy*nx + σ_yy*ny + σ_yz*nz
    traction[2] = Q_mirror[5]*nx + Q_mirror[4]*ny + Q_mirror[2]*nz;  // σ_xz*nx + σ_yz*ny + σ_zz*nz

    // Mirror stress to enforce zero traction
    Q_mirror[0] -= 2.0 * traction[0] * nx;
    Q_mirror[1] -= 2.0 * traction[1] * ny;
    Q_mirror[2] -= 2.0 * traction[2] * nz;
    Q_mirror[3] -= traction[0] * ny + traction[1] * nx;
    Q_mirror[4] -= traction[1] * nz + traction[2] * ny;
    Q_mirror[5] -= traction[0] * nz + traction[2] * nx;

    // Velocity unchanged at free surface
    // Q_mirror[6,7,8] = Q_e[6,7,8] (already set)

    // Compute boundary flux contribution
    for (int v = 0; v < n_var; v++) {
        double jump = Q_mirror[v] - Q_e[v];
        double boundary_flux = -0.5 * jump * area;  // Simplified
        atomicAdd(&R_e[v], boundary_flux);
    }
}

// ============================================================================
// ADER Predictor
// ============================================================================

__global__ void computeADERPredictor(
    const double* __restrict__ Q,
    const ElementMaterial* __restrict__ material,
    const double* __restrict__ stiffness_x,
    const double* __restrict__ stiffness_y,
    const double* __restrict__ stiffness_z,
    const double* __restrict__ inv_jacobian,
    double dt,
    double* __restrict__ Q_predicted,
    int n_elem,
    int n_dof,
    int n_var,
    int ader_order
) {
    int elem = blockIdx.x * blockDim.x + threadIdx.x;
    if (elem >= n_elem) return;

    double lambda = material[elem].lambda;
    double mu_val = material[elem].mu;
    double rho_inv = 1.0 / material[elem].rho;
    const double* J_inv = &inv_jacobian[elem * 9];

    const double* Q_e = &Q[elem * n_dof * n_var];
    double* Qp_e = &Q_predicted[elem * n_dof * n_var];

    // Initialise predicted solution with current solution
    for (int i = 0; i < n_dof * n_var; i++) {
        Qp_e[i] = Q_e[i];
    }

    // Cauchy-Kovalewski procedure: compute successive time derivatives
    // dQ^(k)/dt^k from the PDE and accumulate the Taylor series
    //
    // Q̃ = Q + Δt · L(Q) + Δt²/2 · L(L(Q)) + ...
    //
    // For efficiency, we use temporary buffers within shared memory if possible.
    // For now, use registers (limited to low orders).

    // Temporary: time derivative at current level
    double dQdt[N_VAR_ELASTIC];  // Only handle elastic for now
    double factorial = 1.0;
    double dt_power = 1.0;

    // We build time derivatives iteratively
    // Level 0: dQ/dt = L(Q)  where L is the spatial operator
    // We store the "current derivative" in a local buffer and keep
    // accumulating into Q_predicted

    // For simplicity at order 1 (forward Euler equivalent), just compute L(Q):
    for (int k = 1; k <= ader_order; k++) {
        dt_power *= dt;
        factorial *= k;
        double coeff = dt_power / factorial;

        // For each DOF, compute the spatial operator contribution
        // This is a simplified version; full implementation would use
        // the complete stiffness matrix multiplication
        for (int j = 0; j < n_dof; j++) {
            double Q_local[N_VAR_ELASTIC];
            for (int v = 0; v < n_var && v < N_VAR_ELASTIC; v++) {
                Q_local[v] = Q_e[j * n_var + v];
            }

            double Fx[N_VAR_ELASTIC], Fy[N_VAR_ELASTIC], Fz[N_VAR_ELASTIC];
            elasticFlux(Q_local, lambda, mu_val, rho_inv, Fx, Fy, Fz);

            // Approximate time derivative from flux divergence
            for (int v = 0; v < n_var && v < N_VAR_ELASTIC; v++) {
                // Simplified: use diagonal stiffness approximation
                double div_F = J_inv[0] * Fx[v] + J_inv[4] * Fy[v] + J_inv[8] * Fz[v];
                Qp_e[j * n_var + v] += coeff * div_F;
            }
        }
    }
}

// ============================================================================
// Solution Update
// ============================================================================

__global__ void updateSolutionDG(
    double* __restrict__ Q,
    const double* __restrict__ R_vol,
    const double* __restrict__ R_surf,
    const double* __restrict__ R_src,
    const double* __restrict__ mass_inv,
    const double* __restrict__ det_jacobian,
    double dt,
    int n_elem,
    int n_dof,
    int n_var
) {
    int elem = blockIdx.x * blockDim.x + threadIdx.x;
    if (elem >= n_elem) return;

    double det_J_inv = 1.0 / det_jacobian[elem];

    double* Q_e = &Q[elem * n_dof * n_var];
    const double* Rv = &R_vol[elem * n_dof * n_var];
    const double* Rs = &R_surf[elem * n_dof * n_var];

    for (int i = 0; i < n_dof; i++) {
        double m_inv = mass_inv[i] * det_J_inv;

        for (int v = 0; v < n_var; v++) {
            int idx = i * n_var + v;
            double rhs = Rv[idx] + Rs[idx];

            // Add source term if provided
            if (R_src != nullptr) {
                rhs += R_src[idx];  // R_src is indexed globally, not per-element here
            }

            Q_e[idx] += dt * m_inv * rhs;
        }
    }
}

// ============================================================================
// Source Term Injection
// ============================================================================

__global__ void applySourceTermDG(
    double* __restrict__ R_src,
    int source_elem_id,
    int source_dof_id,
    const double* __restrict__ moment_tensor,
    double source_amplitude,
    const double* __restrict__ basis_at_source,
    double det_jacobian,
    int n_dof,
    int n_var
) {
    // Single-thread kernel (only one source element)
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid != 0) return;

    double* R_e = &R_src[source_elem_id * n_dof * n_var];
    double det_J_inv = 1.0 / det_jacobian;

    // Moment tensor components: Mxx, Myy, Mzz, Mxy, Myz, Mxz
    double Mxx = moment_tensor[0] * source_amplitude;
    double Myy = moment_tensor[1] * source_amplitude;
    double Mzz = moment_tensor[2] * source_amplitude;
    double Mxy = moment_tensor[3] * source_amplitude;
    double Myz = moment_tensor[4] * source_amplitude;
    double Mxz = moment_tensor[5] * source_amplitude;

    // Distribute over DOFs weighted by basis functions
    for (int i = 0; i < n_dof; i++) {
        double phi = basis_at_source[i] * det_J_inv;

        // Source injects into stress rate equations
        R_e[i * n_var + 0] += Mxx * phi;  // σ_xx rate
        R_e[i * n_var + 1] += Myy * phi;  // σ_yy rate
        R_e[i * n_var + 2] += Mzz * phi;  // σ_zz rate
        R_e[i * n_var + 3] += Mxy * phi;  // σ_xy rate
        R_e[i * n_var + 4] += Myz * phi;  // σ_yz rate
        R_e[i * n_var + 5] += Mxz * phi;  // σ_xz rate
    }
}

__global__ void applyMomentTensorSourceDG(
    double* __restrict__ R_src,
    int source_elem_id,
    const double* __restrict__ moment_rate_tensor,
    const double* __restrict__ basis_at_source,
    double det_jacobian,
    int n_dof,
    int n_var
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid != 0) return;

    double* R_e = &R_src[source_elem_id * n_dof * n_var];
    double det_J_inv = 1.0 / det_jacobian;

    for (int i = 0; i < n_dof; i++) {
        double phi = basis_at_source[i] * det_J_inv;
        for (int v = 0; v < 6 && v < n_var; v++) {
            R_e[i * n_var + v] += moment_rate_tensor[v] * phi;
        }
    }
}

// ============================================================================
// PML Absorbing Boundary
// ============================================================================

__global__ void applyPMLDamping(
    double* __restrict__ Q,
    double* __restrict__ Q_pml,
    const int* __restrict__ pml_elem_ids,
    const double* __restrict__ pml_damping,
    double dt,
    int n_pml_elem,
    int n_dof,
    int n_var,
    int n_pml_var
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_pml_elem) return;

    int elem = pml_elem_ids[idx];

    // Damping coefficients for this element (one per spatial direction)
    double dx = pml_damping[idx * 3 + 0];
    double dy = pml_damping[idx * 3 + 1];
    double dz = pml_damping[idx * 3 + 2];

    double* Q_e = &Q[elem * n_dof * n_var];
    double* P_e = &Q_pml[idx * n_dof * n_pml_var];

    // Apply exponential damping to all DOFs
    for (int i = 0; i < n_dof; i++) {
        // Damp stress components
        for (int v = 0; v < 6 && v < n_var; v++) {
            double damping = dx + dy + dz;  // Total damping (simplified)
            double decay = exp(-damping * dt);
            Q_e[i * n_var + v] *= decay;
        }

        // Damp velocity components directionally
        if (n_var >= 9) {
            Q_e[i * n_var + 6] *= exp(-dx * dt);  // v_x damped by x-PML
            Q_e[i * n_var + 7] *= exp(-dy * dt);  // v_y damped by y-PML
            Q_e[i * n_var + 8] *= exp(-dz * dt);  // v_z damped by z-PML
        }

        // Update PML memory variables (for split-field PML)
        if (n_pml_var > 0) {
            for (int p = 0; p < n_pml_var; p++) {
                double total_d = dx + dy + dz;
                P_e[i * n_pml_var + p] *= exp(-total_d * dt);
            }
        }
    }
}

// ============================================================================
// Local Time Stepping
// ============================================================================

__global__ void computeLocalTimestep(
    const ElementMaterial* __restrict__ material,
    const double* __restrict__ elem_size,
    double cfl,
    double* __restrict__ dt_local,
    int n_elem
) {
    int elem = blockIdx.x * blockDim.x + threadIdx.x;
    if (elem >= n_elem) return;

    double Vp, Vs;
    waveVelocities(material[elem].lambda, material[elem].mu,
                    material[elem].rho, Vp, Vs);

    double h = elem_size[elem];
    dt_local[elem] = cfl * h / Vp;
}

__global__ void synchroniseLTSClusters(
    double* __restrict__ Q,
    const int* __restrict__ cluster_ids,
    const int* __restrict__ cluster_boundary_faces,
    const double* __restrict__ flux_buffer,
    int n_boundary_faces,
    int n_dof,
    int n_var
) {
    int face = blockIdx.x * blockDim.x + threadIdx.x;
    if (face >= n_boundary_faces) return;

    // Apply buffered flux from finer cluster to coarser cluster boundary
    // This ensures conservation at LTS cluster interfaces
    int face_id = cluster_boundary_faces[face];
    // Simplified: direct flux application would go here
    // Full implementation requires cluster-level bookkeeping
    (void)Q; (void)cluster_ids; (void)flux_buffer;
    (void)face_id; (void)n_dof; (void)n_var;
}

// ============================================================================
// Receiver Data Extraction
// ============================================================================

__global__ void extractReceiverData(
    const double* __restrict__ Q,
    const int* __restrict__ receiver_elem,
    const double* __restrict__ receiver_basis,
    double* __restrict__ receiver_data,
    int n_receivers,
    int n_dof,
    int n_var
) {
    int rcv = blockIdx.x * blockDim.x + threadIdx.x;
    if (rcv >= n_receivers) return;

    int elem = receiver_elem[rcv];
    const double* Q_e = &Q[elem * n_dof * n_var];
    const double* phi = &receiver_basis[rcv * n_dof];
    double* out = &receiver_data[rcv * n_var];

    // Interpolate: Q(x_rcv) = ∑_i φ_i(x_rcv) Q_i
    for (int v = 0; v < n_var; v++) {
        double val = 0.0;
        for (int i = 0; i < n_dof; i++) {
            val += phi[i] * Q_e[i * n_var + v];
        }
        out[v] = val;
    }
}

// ============================================================================
// Utility Kernels
// ============================================================================

__global__ void zeroArray(double* data, size_t n) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) data[idx] = 0.0;
}

__global__ void extractComponent(
    const double* __restrict__ Q,
    int component,
    double* __restrict__ output,
    int n_elem,
    int n_dof,
    int n_var
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int total = n_elem * n_dof;
    if (idx >= total) return;

    int elem = idx / n_dof;
    int dof = idx % n_dof;
    output[idx] = Q[elem * n_dof * n_var + dof * n_var + component];
}

__global__ void computeFieldNorm(
    const double* __restrict__ Q,
    int component,
    double* __restrict__ partial_norms,
    int n_elem,
    int n_dof,
    int n_var
) {
    extern __shared__ double sdata[];

    int tid = threadIdx.x;
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int total = n_elem * n_dof;

    double val = 0.0;
    if (idx < total) {
        int elem = idx / n_dof;
        int dof = idx % n_dof;
        double q = Q[elem * n_dof * n_var + dof * n_var + component];
        val = q * q;
    }
    sdata[tid] = val;
    __syncthreads();

    // Parallel reduction
    for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            sdata[tid] += sdata[tid + s];
        }
        __syncthreads();
    }

    if (tid == 0) {
        partial_norms[blockIdx.x] = sdata[0];
    }
}

} // namespace DG
} // namespace GPU
} // namespace FSRM

#endif // USE_CUDA
