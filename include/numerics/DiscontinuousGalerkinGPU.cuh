/**
 * @file DiscontinuousGalerkinGPU.cuh
 * @brief CUDA kernels for GPU-accelerated Discontinuous Galerkin wave equation solver
 *
 * Provides GPU-parallel implementations of the core ADER-DG operations for
 * elastic and poroelastic wave propagation:
 *
 *   - Volume integral computation (stiffness × flux)
 *   - Surface integral / numerical flux (Riemann solver at element interfaces)
 *   - ADER predictor step (Cauchy-Kovalewski local time integration)
 *   - Solution update (inverse mass matrix application)
 *   - Source term injection (explosion / moment tensor)
 *   - PML absorbing boundary damping
 *   - Local time stepping synchronisation
 *
 * The DG formulation is embarrassingly parallel over elements — each CUDA thread
 * processes one element, making this an ideal GPU workload.
 *
 * ## Elastic Wave Equation (velocity-stress formulation)
 *
 * The 3-D elastic wave equation in first-order form:
 *
 *   ∂σ/∂t = C : ∇v            (Hooke's law rate form)
 *   ρ ∂v/∂t = ∇·σ + f         (momentum balance)
 *
 * where σ is the stress tensor (6 independent components in Voigt notation),
 * v is the velocity vector (3 components), C is the elastic stiffness tensor,
 * ρ is density, and f is the body force (source term).
 *
 * Conservative variable vector Q = (σ_xx, σ_yy, σ_zz, σ_xy, σ_yz, σ_xz, v_x, v_y, v_z)
 * has 9 components per element DOF.
 *
 * ## Memory Layout
 *
 * Element-major ordering:  data[elem * n_dof * n_var + dof * n_var + var]
 * This layout ensures coalesced memory access when threads are mapped to elements.
 */

#ifndef DISCONTINUOUS_GALERKIN_GPU_CUH
#define DISCONTINUOUS_GALERKIN_GPU_CUH

#ifdef USE_CUDA

#include <cuda_runtime.h>
#include <cuda.h>

namespace FSRM {
namespace GPU {
namespace DG {

// ============================================================================
// Constants
// ============================================================================

/// Maximum polynomial order supported on GPU
constexpr int MAX_DG_ORDER = 6;

/// Maximum DOFs per element for 3D hexahedra at MAX_DG_ORDER: (p+1)^3
constexpr int MAX_DOFS_PER_ELEM = (MAX_DG_ORDER + 1) * (MAX_DG_ORDER + 1) * (MAX_DG_ORDER + 1);

/// Number of conservative variables for elastic wave equation (σ_6 + v_3)
constexpr int N_VAR_ELASTIC = 9;

/// Number of conservative variables for poroelastic wave equation (σ_6 + v_3 + p + q_3)
constexpr int N_VAR_POROELASTIC = 13;

/// Number of faces per hexahedral element
constexpr int FACES_PER_HEX = 6;

/// Default CUDA block size for element-parallel kernels
constexpr int DG_BLOCK_SIZE = 128;

// ============================================================================
// Element Material Properties (device-side structure)
// ============================================================================

/**
 * @brief Per-element material data for elastic wave propagation
 */
struct ElementMaterial {
    double rho;             ///< Density (kg/m³)
    double lambda;          ///< First Lamé parameter (Pa)
    double mu;              ///< Shear modulus / second Lamé parameter (Pa)
    double Qp;              ///< P-wave quality factor (attenuation)
    double Qs;              ///< S-wave quality factor (attenuation)
};

/**
 * @brief Per-element Biot poroelastic parameters
 */
struct ElementBiotMaterial {
    double rho_solid;       ///< Solid grain density
    double rho_fluid;       ///< Fluid density
    double lambda_u;        ///< Undrained first Lamé parameter
    double mu;              ///< Shear modulus (same drained/undrained)
    double alpha;           ///< Biot coefficient
    double M;               ///< Biot modulus
    double kappa;           ///< Permeability / viscosity (m²/Pa·s)
    double porosity;        ///< Porosity
};

// ============================================================================
// Volume Integral Kernels
// ============================================================================

/**
 * @brief Compute the DG volume integral for elastic wave propagation
 *
 * For each element (one thread per element), evaluates:
 *   R_vol[e] = ∑_q w_q |J_e| ∇φ(x_q) · F(Q(x_q))
 *
 * where F is the flux function for the elastic wave equation and
 * ∇φ are basis function gradients at quadrature points.
 *
 * This is equivalent to applying the stiffness matrices K_ξ, K_η, K_ζ
 * to the flux evaluated at DOF locations (nodal DG with collocated
 * quadrature).
 *
 * @param Q            Solution array [n_elem × n_dof × N_VAR_ELASTIC]
 * @param material     Per-element material properties [n_elem]
 * @param stiffness_x  Reference stiffness matrix K_ξ [n_dof × n_dof]
 * @param stiffness_y  Reference stiffness matrix K_η [n_dof × n_dof]
 * @param stiffness_z  Reference stiffness matrix K_ζ [n_dof × n_dof]
 * @param inv_jacobian Inverse Jacobian per element [n_elem × 9] (3×3 row-major)
 * @param det_jacobian Jacobian determinant per element [n_elem]
 * @param R_vol        Output volume residual [n_elem × n_dof × N_VAR_ELASTIC]
 * @param n_elem       Number of elements
 * @param n_dof        DOFs per element
 */
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
);

/**
 * @brief Compute the DG volume integral for poroelastic wave propagation
 *
 * Same structure as elastic but with Biot-coupled equations:
 *   σ' evolution, fluid mass balance, Darcy flux, etc.
 */
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
);

// ============================================================================
// Numerical Flux / Surface Integral Kernels
// ============================================================================

/**
 * @brief Compute the Godunov / Rusanov numerical flux at element interfaces
 *
 * For each interior face, computes the upwind flux:
 *   F* = 0.5 * (F(Q⁻) + F(Q⁺)) - 0.5 * |A| * (Q⁺ - Q⁻)
 *
 * where |A| is the maximum wave speed (Rusanov/Lax-Friedrichs flux).
 * For the elastic wave equation, |A| = max(Vp_left, Vp_right).
 *
 * The surface integral contribution is then:
 *   R_surf[e] = ∑_f ∑_q w_q φ(x_q) F*(Q⁻, Q⁺, n) dS
 *
 * @param Q              Solution [n_elem × n_dof × n_var]
 * @param material       Material properties [n_elem]
 * @param face_neighbors Face connectivity [n_interior_faces × 2] (left/right elem)
 * @param face_dof_map   DOF mapping for face quadrature [n_face_dof]
 * @param face_normals   Outward face normals [n_interior_faces × 3]
 * @param face_area      Face areas [n_interior_faces]
 * @param mass_inv       Inverse mass matrix diagonal [n_dof] (for diagonal mass)
 * @param R_surf         Output surface residual [n_elem × n_dof × n_var]
 * @param n_interior_faces Number of interior faces
 * @param n_dof          DOFs per element
 * @param n_var          Variables per DOF
 */
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
);

/**
 * @brief Apply free-surface boundary condition
 *
 * At the free surface, traction vanishes: σ·n = 0.
 * This is enforced weakly through the numerical flux by using a
 * mirror state Q⁺ that negates the stress components and mirrors velocity.
 */
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
);

// ============================================================================
// ADER Predictor Step
// ============================================================================

/**
 * @brief Compute the ADER (Cauchy-Kovalewski) local time predictor
 *
 * For each element (independently, no communication needed), compute
 * a space-time Taylor expansion of the solution:
 *
 *   Q̃(x,t) = ∑_{k=0}^{p} (Δt)^k / k! · ∂^k Q / ∂t^k (x, tⁿ)
 *
 * Time derivatives are computed recursively using the PDE:
 *   ∂Q/∂t = L(Q)    ← spatial operator from DG discretisation
 *   ∂²Q/∂t² = L(∂Q/∂t)
 *   etc.
 *
 * This yields a high-order-in-time prediction without global communication,
 * which is then used in the numerical flux computation.
 *
 * @param Q            Current solution [n_elem × n_dof × n_var]
 * @param material     Material properties [n_elem]
 * @param stiffness_x  Stiffness matrices [n_dof × n_dof]
 * @param stiffness_y
 * @param stiffness_z
 * @param inv_jacobian  Inverse Jacobian [n_elem × 9]
 * @param dt            Time step
 * @param Q_predicted   Output predicted solution [n_elem × n_dof × n_var]
 * @param n_elem        Number of elements
 * @param n_dof         DOFs per element
 * @param n_var         Variables per DOF
 * @param ader_order    Number of ADER derivatives (= DG order)
 */
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
);

// ============================================================================
// Solution Update
// ============================================================================

/**
 * @brief Update element DOFs: Q^{n+1} = Q^n + Δt M⁻¹ (R_vol + R_surf + R_src)
 *
 * Applies the inverse mass matrix and accumulates all residual contributions.
 * For nodal DG with diagonal mass matrix, this is a simple element-wise operation.
 *
 * @param Q          Solution [n_elem × n_dof × n_var] (updated in-place)
 * @param R_vol      Volume integral residual
 * @param R_surf     Surface integral residual
 * @param R_src      Source term residual (can be NULL if no source)
 * @param mass_inv   Inverse mass matrix diagonal [n_dof]
 * @param det_jacobian Jacobian determinant [n_elem]
 * @param dt         Time step
 * @param n_elem     Number of elements
 * @param n_dof      DOFs per element
 * @param n_var      Variables per DOF
 */
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
);

// ============================================================================
// Source Term Injection
// ============================================================================

/**
 * @brief Inject an explosion source via equivalent body forces
 *
 * Implements the Mueller-Murphy RDP source as an equivalent moment tensor
 * applied at the source element. The source time function S(t) is evaluated
 * at the current time and the resulting moment rate tensor Ṁ is distributed
 * over the source element's DOFs.
 *
 * For an isotropic explosion: M_ij = M_0 · S(t) · δ_ij
 *   → injects pressure (σ_xx = σ_yy = σ_zz) at the source point.
 *
 * @param R_src           Source residual [n_elem × n_dof × n_var] (accumulated)
 * @param source_elem_id  Index of element containing the source
 * @param source_dof_id   Nearest DOF within source element
 * @param moment_tensor   6-component moment tensor rate [Mxx, Myy, Mzz, Mxy, Myz, Mxz]
 * @param source_amplitude Current source amplitude (S(t) × M_0)
 * @param basis_at_source  Basis function values at source point [n_dof]
 * @param det_jacobian     Jacobian determinant at source element
 * @param n_dof            DOFs per element
 * @param n_var            Variables per DOF
 */
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
);

/**
 * @brief Inject a moment tensor source (general double-couple or CLVD)
 *
 * Supports full 6-component moment tensor for modelling:
 * - Isotropic (explosion)
 * - CLVD (compensated linear vector dipole)
 * - Double-couple (earthquake/tectonic release)
 * - Any combination
 */
__global__ void applyMomentTensorSourceDG(
    double* __restrict__ R_src,
    int source_elem_id,
    const double* __restrict__ moment_rate_tensor,
    const double* __restrict__ basis_at_source,
    double det_jacobian,
    int n_dof,
    int n_var
);

// ============================================================================
// PML Absorbing Boundary
// ============================================================================

/**
 * @brief Apply Perfectly Matched Layer (PML) damping to boundary elements
 *
 * For elements within the PML region, applies coordinate-stretching damping:
 *   σ̃_i = σ_i · (1 + d_i(x)/iω)
 *
 * In the time domain this becomes auxiliary differential equations that
 * damp outgoing waves exponentially within the PML layer.
 *
 * @param Q              Solution [n_elem × n_dof × n_var]
 * @param Q_pml          PML auxiliary variables [n_pml_elem × n_dof × n_pml_var]
 * @param pml_elem_ids   Indices of elements in PML region [n_pml_elem]
 * @param pml_damping    Damping profile σ(x) for each PML element [n_pml_elem × 3]
 * @param dt             Time step
 * @param n_pml_elem     Number of PML elements
 * @param n_dof          DOFs per element
 * @param n_var          Solution variables per DOF
 * @param n_pml_var      PML auxiliary variables per DOF
 */
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
);

// ============================================================================
// Local Time Stepping
// ============================================================================

/**
 * @brief Compute element-local CFL-stable timestep
 *
 * For each element, computes:
 *   dt_e = CFL × h_e / max(Vp_e)
 *
 * where h_e is the element size (minimum edge length / inscribed sphere radius)
 * and Vp_e is the P-wave velocity.
 *
 * @param material     Material properties [n_elem]
 * @param elem_size    Element characteristic size [n_elem]
 * @param cfl          CFL number
 * @param dt_local     Output per-element timestep [n_elem]
 * @param n_elem       Number of elements
 */
__global__ void computeLocalTimestep(
    const ElementMaterial* __restrict__ material,
    const double* __restrict__ elem_size,
    double cfl,
    double* __restrict__ dt_local,
    int n_elem
);

/**
 * @brief Synchronise elements at LTS cluster boundaries
 *
 * For local time stepping, elements are grouped into clusters with
 * timestep ratios of 2:1. At cluster boundaries, elements must exchange
 * flux data at the finest timestep resolution.
 */
__global__ void synchroniseLTSClusters(
    double* __restrict__ Q,
    const int* __restrict__ cluster_ids,
    const int* __restrict__ cluster_boundary_faces,
    const double* __restrict__ flux_buffer,
    int n_boundary_faces,
    int n_dof,
    int n_var
);

// ============================================================================
// Receiver / Seismogram Extraction
// ============================================================================

/**
 * @brief Extract wavefield values at receiver locations
 *
 * For each receiver, interpolates the solution Q at the receiver position
 * using the element's basis functions.
 *
 * @param Q              Solution [n_elem × n_dof × n_var]
 * @param receiver_elem  Element containing each receiver [n_receivers]
 * @param receiver_basis Basis function values at each receiver [n_receivers × n_dof]
 * @param receiver_data  Output: interpolated solution [n_receivers × n_var]
 * @param n_receivers    Number of receivers
 * @param n_dof          DOFs per element
 * @param n_var          Variables per DOF
 */
__global__ void extractReceiverData(
    const double* __restrict__ Q,
    const int* __restrict__ receiver_elem,
    const double* __restrict__ receiver_basis,
    double* __restrict__ receiver_data,
    int n_receivers,
    int n_dof,
    int n_var
);

// ============================================================================
// Utility / Launch Helpers
// ============================================================================

/**
 * @brief Zero-initialise a device array
 */
__global__ void zeroArray(double* data, size_t n);

/**
 * @brief Copy solution for output (device to device subset copy)
 */
__global__ void extractComponent(
    const double* __restrict__ Q,
    int component,
    double* __restrict__ output,
    int n_elem,
    int n_dof,
    int n_var
);

/**
 * @brief Compute L2 norm of a field component (reduction kernel)
 */
__global__ void computeFieldNorm(
    const double* __restrict__ Q,
    int component,
    double* __restrict__ partial_norms,
    int n_elem,
    int n_dof,
    int n_var
);

// ============================================================================
// Launch Configuration Helpers
// ============================================================================

/**
 * @brief Get optimal grid dimensions for element-parallel kernels
 */
inline dim3 getDGGrid(int n_elem, int block_size = DG_BLOCK_SIZE) {
    return dim3((n_elem + block_size - 1) / block_size);
}

/**
 * @brief Get optimal grid dimensions for face-parallel kernels
 */
inline dim3 getFaceGrid(int n_faces, int block_size = DG_BLOCK_SIZE) {
    return dim3((n_faces + block_size - 1) / block_size);
}

} // namespace DG
} // namespace GPU
} // namespace FSRM

#endif // USE_CUDA

#endif // DISCONTINUOUS_GALERKIN_GPU_CUH
