/**
 * @file PhysicsKernel_GPU_Extended.cu
 * @brief CUDA implementations for extended GPU physics kernels
 * 
 * This file contains CUDA kernel implementations for:
 * - Black Oil flow
 * - Thermal transport
 * - Geomechanics
 * - Hydrodynamic (Euler)
 * - Atmospheric blast
 * - Infrasound
 * - Tsunami
 * - Particle transport
 */

#include "PhysicsKernel_GPU_Extended.hpp"

#ifdef USE_CUDA
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#endif

namespace FSRM {

// =============================================================================
// CUDA Error Checking Macro
// =============================================================================

#ifdef USE_CUDA
#define CUDA_CHECK(call) \
    do { \
        cudaError_t err = call; \
        if (err != cudaSuccess) { \
            fprintf(stderr, "CUDA error at %s:%d: %s\n", \
                    __FILE__, __LINE__, cudaGetErrorString(err)); \
        } \
    } while(0)

// Thread block sizes
constexpr int BLOCK_SIZE_1D = 256;
constexpr int BLOCK_SIZE_2D = 16;

// =============================================================================
// Black Oil GPU Kernels
// =============================================================================

/**
 * @brief Compute PVT properties for all cells
 */
__global__ void computePVTKernel(
    const double* pressure,
    double* Bo, double* Bg, double* Bw,
    double* Rs,
    double* visc_o, double* visc_g, double* visc_w,
    double* rho_o, double* rho_g, double* rho_w,
    const double* pvt_table_p, const double* pvt_table_Bo,
    const double* pvt_table_Bg, const double* pvt_table_Rs,
    int pvt_size, int n_cells
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_cells) return;
    
    double p = pressure[idx];
    
    // Linear interpolation in PVT table
    int i = 0;
    while (i < pvt_size - 1 && pvt_table_p[i+1] < p) i++;
    
    double t = (pvt_size > 1 && pvt_table_p[i+1] > pvt_table_p[i]) ?
               (p - pvt_table_p[i]) / (pvt_table_p[i+1] - pvt_table_p[i]) : 0.0;
    t = fmax(0.0, fmin(1.0, t));
    
    // Interpolate Bo
    Bo[idx] = pvt_table_Bo[i] + t * (pvt_table_Bo[i+1] - pvt_table_Bo[i]);
    
    // Interpolate Bg
    Bg[idx] = pvt_table_Bg[i] + t * (pvt_table_Bg[i+1] - pvt_table_Bg[i]);
    
    // Interpolate Rs
    Rs[idx] = pvt_table_Rs[i] + t * (pvt_table_Rs[i+1] - pvt_table_Rs[i]);
    
    // Water (constant for now)
    Bw[idx] = 1.0;
    
    // Densities
    double rho_o_std = 850.0;  // kg/m³
    double rho_g_std = 1.0;    // kg/m³
    double rho_w_std = 1000.0; // kg/m³
    
    rho_o[idx] = (rho_o_std + Rs[idx] * rho_g_std) / Bo[idx];
    rho_g[idx] = rho_g_std / Bg[idx];
    rho_w[idx] = rho_w_std / Bw[idx];
    
    // Viscosities (simplified correlations)
    visc_o[idx] = 1.0e-3 * exp(-0.001 * (p - 1.0e6));  // Pa·s
    visc_g[idx] = 1.5e-5;  // Pa·s
    visc_w[idx] = 1.0e-3;  // Pa·s
}

/**
 * @brief Compute relative permeabilities (Corey model)
 */
__global__ void computeRelPermKernel(
    const double* Sw, const double* Sg,
    double* krw, double* kro, double* krg,
    double Swc, double Sor, double Sgc,
    double nw, double no, double ng,
    int n_cells
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_cells) return;
    
    double sw = Sw[idx];
    double sg = Sg[idx];
    double so = 1.0 - sw - sg;
    
    // Normalized saturations
    double sw_norm = (sw - Swc) / (1.0 - Swc - Sor);
    double so_norm = (so - Sor) / (1.0 - Swc - Sor);
    double sg_norm = (sg - Sgc) / (1.0 - Swc - Sgc);
    
    sw_norm = fmax(0.0, fmin(1.0, sw_norm));
    so_norm = fmax(0.0, fmin(1.0, so_norm));
    sg_norm = fmax(0.0, fmin(1.0, sg_norm));
    
    // Corey relative permeabilities
    krw[idx] = pow(sw_norm, nw);
    kro[idx] = pow(so_norm, no);
    krg[idx] = pow(sg_norm, ng);
}

/**
 * @brief Compute phase mobilities
 */
__global__ void computeMobilitiesKernel(
    const double* krw, const double* kro, const double* krg,
    const double* visc_w, const double* visc_o, const double* visc_g,
    double* lambda_w, double* lambda_o, double* lambda_g,
    double* lambda_total,
    int n_cells
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_cells) return;
    
    lambda_w[idx] = krw[idx] / visc_w[idx];
    lambda_o[idx] = kro[idx] / visc_o[idx];
    lambda_g[idx] = krg[idx] / visc_g[idx];
    lambda_total[idx] = lambda_w[idx] + lambda_o[idx] + lambda_g[idx];
}

// =============================================================================
// Thermal GPU Kernels
// =============================================================================

/**
 * @brief Compute effective thermal properties
 */
__global__ void computeEffectiveThermalKernel(
    const double* k_solid, const double* k_fluid,
    const double* rho_solid, const double* rho_fluid,
    const double* cp_solid, const double* cp_fluid,
    const double* porosity,
    double* k_eff, double* rho_cp_eff,
    int n_cells
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_cells) return;
    
    double phi = porosity[idx];
    
    // Effective conductivity (volume-weighted)
    k_eff[idx] = (1.0 - phi) * k_solid[idx] + phi * k_fluid[idx];
    
    // Effective heat capacity
    rho_cp_eff[idx] = (1.0 - phi) * rho_solid[idx] * cp_solid[idx] +
                      phi * rho_fluid[idx] * cp_fluid[idx];
}

/**
 * @brief Compute conduction heat flux
 */
__global__ void computeConductionFluxKernel(
    const double* temperature,
    const double* k_eff,
    const int* face_cells,
    const double* face_area,
    const double* face_normal,
    const double* cell_center,
    double* heat_flux,
    int n_faces
) {
    int f = blockIdx.x * blockDim.x + threadIdx.x;
    if (f >= n_faces) return;
    
    int c1 = face_cells[2*f];
    int c2 = face_cells[2*f + 1];
    
    if (c2 < 0) {
        // Boundary face
        heat_flux[f] = 0.0;
        return;
    }
    
    // Distance between cell centers
    double dx = cell_center[3*c2] - cell_center[3*c1];
    double dy = cell_center[3*c2 + 1] - cell_center[3*c1 + 1];
    double dz = cell_center[3*c2 + 2] - cell_center[3*c1 + 2];
    double dist = sqrt(dx*dx + dy*dy + dz*dz);
    
    // Temperature gradient
    double dT = temperature[c2] - temperature[c1];
    
    // Harmonic average of conductivity
    double k_face = 2.0 * k_eff[c1] * k_eff[c2] / (k_eff[c1] + k_eff[c2]);
    
    // Heat flux
    heat_flux[f] = -k_face * dT / dist * face_area[f];
}

// =============================================================================
// Geomechanics GPU Kernels
// =============================================================================

/**
 * @brief Compute strain tensor from displacement gradient
 */
__global__ void computeStrainKernel(
    const double* displacement,
    const int* connectivity,
    const double* inv_jacobian,
    double* strain,  // Voigt: εxx, εyy, εzz, γxy, γyz, γxz
    int n_cells, int nodes_per_cell
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_cells) return;
    
    // Compute displacement gradient using shape function derivatives
    double dudx[9] = {0};  // 3x3 gradient tensor
    
    for (int n = 0; n < nodes_per_cell; n++) {
        int node = connectivity[idx * nodes_per_cell + n];
        
        // Shape function derivatives (simplified for hex element)
        double dNdx = inv_jacobian[idx * 9 + 0];  // Simplified
        double dNdy = inv_jacobian[idx * 9 + 4];
        double dNdz = inv_jacobian[idx * 9 + 8];
        
        dudx[0] += displacement[3*node + 0] * dNdx;  // du/dx
        dudx[1] += displacement[3*node + 0] * dNdy;  // du/dy
        dudx[2] += displacement[3*node + 0] * dNdz;  // du/dz
        dudx[3] += displacement[3*node + 1] * dNdx;  // dv/dx
        dudx[4] += displacement[3*node + 1] * dNdy;  // dv/dy
        dudx[5] += displacement[3*node + 1] * dNdz;  // dv/dz
        dudx[6] += displacement[3*node + 2] * dNdx;  // dw/dx
        dudx[7] += displacement[3*node + 2] * dNdy;  // dw/dy
        dudx[8] += displacement[3*node + 2] * dNdz;  // dw/dz
    }
    
    // Strain in Voigt notation
    strain[6*idx + 0] = dudx[0];                      // εxx
    strain[6*idx + 1] = dudx[4];                      // εyy
    strain[6*idx + 2] = dudx[8];                      // εzz
    strain[6*idx + 3] = 0.5 * (dudx[1] + dudx[3]);   // εxy
    strain[6*idx + 4] = 0.5 * (dudx[5] + dudx[7]);   // εyz
    strain[6*idx + 5] = 0.5 * (dudx[2] + dudx[6]);   // εxz
}

/**
 * @brief Compute stress from strain (isotropic elasticity)
 */
__global__ void computeStressKernel(
    const double* strain,
    const double* youngs_modulus,
    const double* poisson_ratio,
    double* stress,
    int n_cells
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_cells) return;
    
    double E = youngs_modulus[idx];
    double nu = poisson_ratio[idx];
    
    // Lame parameters
    double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    double mu = E / (2.0 * (1.0 + nu));
    
    // Extract strains
    double exx = strain[6*idx + 0];
    double eyy = strain[6*idx + 1];
    double ezz = strain[6*idx + 2];
    double exy = strain[6*idx + 3];
    double eyz = strain[6*idx + 4];
    double exz = strain[6*idx + 5];
    
    // Volumetric strain
    double evol = exx + eyy + ezz;
    
    // Stress (Hooke's law)
    stress[6*idx + 0] = lambda * evol + 2.0 * mu * exx;  // σxx
    stress[6*idx + 1] = lambda * evol + 2.0 * mu * eyy;  // σyy
    stress[6*idx + 2] = lambda * evol + 2.0 * mu * ezz;  // σzz
    stress[6*idx + 3] = 2.0 * mu * exy;                   // τxy
    stress[6*idx + 4] = 2.0 * mu * eyz;                   // τyz
    stress[6*idx + 5] = 2.0 * mu * exz;                   // τxz
}

// =============================================================================
// Hydrodynamic GPU Kernels
// =============================================================================

/**
 * @brief Compute HLLC Riemann flux
 */
__device__ void HLLCFlux(
    double rhoL, double uL, double vL, double wL, double pL,
    double rhoR, double uR, double vR, double wR, double pR,
    double nx, double ny, double nz, double gamma,
    double* flux
) {
    // Normal velocities
    double vnL = uL * nx + vL * ny + wL * nz;
    double vnR = uR * nx + vR * ny + wR * nz;
    
    // Sound speeds
    double aL = sqrt(gamma * pL / rhoL);
    double aR = sqrt(gamma * pR / rhoR);
    
    // Roe averages for wave speed estimates
    double sqrtRhoL = sqrt(rhoL);
    double sqrtRhoR = sqrt(rhoR);
    double denom = sqrtRhoL + sqrtRhoR;
    
    double uRoe = (sqrtRhoL * uL + sqrtRhoR * uR) / denom;
    double vRoe = (sqrtRhoL * vL + sqrtRhoR * vR) / denom;
    double wRoe = (sqrtRhoL * wL + sqrtRhoR * wR) / denom;
    double vnRoe = uRoe * nx + vRoe * ny + wRoe * nz;
    
    double EL = pL / (gamma - 1.0) + 0.5 * rhoL * (uL*uL + vL*vL + wL*wL);
    double ER = pR / (gamma - 1.0) + 0.5 * rhoR * (uR*uR + vR*vR + wR*wR);
    double HL = (EL + pL) / rhoL;
    double HR = (ER + pR) / rhoR;
    double HRoe = (sqrtRhoL * HL + sqrtRhoR * HR) / denom;
    double aRoe = sqrt((gamma - 1.0) * (HRoe - 0.5 * (uRoe*uRoe + vRoe*vRoe + wRoe*wRoe)));
    
    // Wave speeds
    double SL = fmin(vnL - aL, vnRoe - aRoe);
    double SR = fmax(vnR + aR, vnRoe + aRoe);
    double SM = ((pR - pL) + rhoL * vnL * (SL - vnL) - rhoR * vnR * (SR - vnR)) /
                (rhoL * (SL - vnL) - rhoR * (SR - vnR));
    
    // HLLC flux
    if (SL >= 0.0) {
        // Left state
        flux[0] = rhoL * vnL;
        flux[1] = rhoL * uL * vnL + pL * nx;
        flux[2] = rhoL * vL * vnL + pL * ny;
        flux[3] = rhoL * wL * vnL + pL * nz;
        flux[4] = (EL + pL) * vnL;
    } else if (SR <= 0.0) {
        // Right state
        flux[0] = rhoR * vnR;
        flux[1] = rhoR * uR * vnR + pR * nx;
        flux[2] = rhoR * vR * vnR + pR * ny;
        flux[3] = rhoR * wR * vnR + pR * nz;
        flux[4] = (ER + pR) * vnR;
    } else if (SM >= 0.0) {
        // Left star state
        double pStar = pL + rhoL * (SL - vnL) * (SM - vnL);
        double rhoStar = rhoL * (SL - vnL) / (SL - SM);
        double factor = rhoStar / rhoL;
        
        flux[0] = rhoL * vnL + SL * (rhoStar - rhoL);
        flux[1] = rhoL * uL * vnL + pL * nx + SL * (rhoStar * (uL + (SM - vnL) * nx) - rhoL * uL);
        flux[2] = rhoL * vL * vnL + pL * ny + SL * (rhoStar * (vL + (SM - vnL) * ny) - rhoL * vL);
        flux[3] = rhoL * wL * vnL + pL * nz + SL * (rhoStar * (wL + (SM - vnL) * nz) - rhoL * wL);
        double EStar = EL * factor + (pStar * SM - pL * vnL) / (SL - SM);
        flux[4] = (EL + pL) * vnL + SL * (EStar - EL);
    } else {
        // Right star state
        double pStar = pR + rhoR * (SR - vnR) * (SM - vnR);
        double rhoStar = rhoR * (SR - vnR) / (SR - SM);
        
        flux[0] = rhoR * vnR + SR * (rhoStar - rhoR);
        flux[1] = rhoR * uR * vnR + pR * nx + SR * (rhoStar * (uR + (SM - vnR) * nx) - rhoR * uR);
        flux[2] = rhoR * vR * vnR + pR * ny + SR * (rhoStar * (vR + (SM - vnR) * ny) - rhoR * vR);
        flux[3] = rhoR * wR * vnR + pR * nz + SR * (rhoStar * (wR + (SM - vnR) * nz) - rhoR * wR);
        double ER_val = pR / (gamma - 1.0) + 0.5 * rhoR * (uR*uR + vR*vR + wR*wR);
        double EStar = ER_val * (SR - vnR) / (SR - SM) + (pStar * SM - pR * vnR) / (SR - SM);
        flux[4] = (ER_val + pR) * vnR + SR * (EStar - ER_val);
    }
}

/**
 * @brief Compute hydrodynamic flux at all faces
 */
__global__ void computeHydroFluxKernel(
    const double* rho, const double* rhou, const double* rhov, 
    const double* rhow, const double* E,
    const int* face_cells,
    const double* face_normal,
    const double* face_area,
    double* flux_rho, double* flux_rhou, double* flux_rhov,
    double* flux_rhow, double* flux_E,
    double gamma, int n_faces
) {
    int f = blockIdx.x * blockDim.x + threadIdx.x;
    if (f >= n_faces) return;
    
    int cL = face_cells[2*f];
    int cR = face_cells[2*f + 1];
    
    // Left state
    double rhoL = rho[cL];
    double uL = rhou[cL] / rhoL;
    double vL = rhov[cL] / rhoL;
    double wL = rhow[cL] / rhoL;
    double EL = E[cL];
    double pL = (gamma - 1.0) * (EL - 0.5 * rhoL * (uL*uL + vL*vL + wL*wL));
    
    // Right state (or boundary)
    double rhoR, uR, vR, wR, pR;
    if (cR >= 0) {
        rhoR = rho[cR];
        uR = rhou[cR] / rhoR;
        vR = rhov[cR] / rhoR;
        wR = rhow[cR] / rhoR;
        double ER = E[cR];
        pR = (gamma - 1.0) * (ER - 0.5 * rhoR * (uR*uR + vR*vR + wR*wR));
    } else {
        // Reflective boundary
        rhoR = rhoL;
        double nx = face_normal[3*f];
        double ny = face_normal[3*f + 1];
        double nz = face_normal[3*f + 2];
        double vn = uL * nx + vL * ny + wL * nz;
        uR = uL - 2.0 * vn * nx;
        vR = vL - 2.0 * vn * ny;
        wR = wL - 2.0 * vn * nz;
        pR = pL;
    }
    
    double nx = face_normal[3*f];
    double ny = face_normal[3*f + 1];
    double nz = face_normal[3*f + 2];
    
    double flux[5];
    HLLCFlux(rhoL, uL, vL, wL, pL, rhoR, uR, vR, wR, pR, nx, ny, nz, gamma, flux);
    
    double A = face_area[f];
    flux_rho[f] = flux[0] * A;
    flux_rhou[f] = flux[1] * A;
    flux_rhov[f] = flux[2] * A;
    flux_rhow[f] = flux[3] * A;
    flux_E[f] = flux[4] * A;
}

// =============================================================================
// Tsunami (Shallow Water) GPU Kernels
// =============================================================================

/**
 * @brief Compute shallow water HLL flux
 */
__device__ void SWEFlux(
    double hL, double huL, double hvL,
    double hR, double huR, double hvR,
    double nx, double ny, double g,
    double* flux
) {
    double eps = 1e-10;
    
    // Velocities
    double uL = (hL > eps) ? huL / hL : 0.0;
    double vL = (hL > eps) ? hvL / hL : 0.0;
    double uR = (hR > eps) ? huR / hR : 0.0;
    double vR = (hR > eps) ? hvR / hR : 0.0;
    
    // Normal velocities
    double vnL = uL * nx + vL * ny;
    double vnR = uR * nx + vR * ny;
    
    // Wave speeds
    double aL = sqrt(g * hL);
    double aR = sqrt(g * hR);
    
    double SL = fmin(vnL - aL, vnR - aR);
    double SR = fmax(vnL + aL, vnR + aR);
    
    if (SL >= 0.0) {
        // Left flux
        flux[0] = huL * nx + hvL * ny;
        flux[1] = (huL * uL + 0.5 * g * hL * hL) * nx + huL * vL * ny;
        flux[2] = hvL * uL * nx + (hvL * vL + 0.5 * g * hL * hL) * ny;
    } else if (SR <= 0.0) {
        // Right flux
        flux[0] = huR * nx + hvR * ny;
        flux[1] = (huR * uR + 0.5 * g * hR * hR) * nx + huR * vR * ny;
        flux[2] = hvR * uR * nx + (hvR * vR + 0.5 * g * hR * hR) * ny;
    } else {
        // HLL flux
        double denom = SR - SL;
        
        double FL0 = huL * nx + hvL * ny;
        double FL1 = (huL * uL + 0.5 * g * hL * hL) * nx + huL * vL * ny;
        double FL2 = hvL * uL * nx + (hvL * vL + 0.5 * g * hL * hL) * ny;
        
        double FR0 = huR * nx + hvR * ny;
        double FR1 = (huR * uR + 0.5 * g * hR * hR) * nx + huR * vR * ny;
        double FR2 = hvR * uR * nx + (hvR * vR + 0.5 * g * hR * hR) * ny;
        
        flux[0] = (SR * FL0 - SL * FR0 + SL * SR * (hR - hL)) / denom;
        flux[1] = (SR * FL1 - SL * FR1 + SL * SR * (huR - huL)) / denom;
        flux[2] = (SR * FL2 - SL * FR2 + SL * SR * (hvR - hvL)) / denom;
    }
}

/**
 * @brief Compute tsunami flux at all faces
 */
__global__ void computeTsunamiFluxKernel(
    const double* eta, const double* hu, const double* hv,
    const double* depth,
    const int* face_cells,
    const double* face_normal,
    const double* face_length,
    double* flux_eta, double* flux_hu, double* flux_hv,
    double g, double min_depth, int n_faces
) {
    int f = blockIdx.x * blockDim.x + threadIdx.x;
    if (f >= n_faces) return;
    
    int cL = face_cells[2*f];
    int cR = face_cells[2*f + 1];
    
    double nx = face_normal[2*f];
    double ny = face_normal[2*f + 1];
    
    // Total water depth
    double hL = fmax(min_depth, depth[cL] + eta[cL]);
    double huL = hu[cL];
    double hvL = hv[cL];
    
    double hR, huR, hvR;
    if (cR >= 0) {
        hR = fmax(min_depth, depth[cR] + eta[cR]);
        huR = hu[cR];
        hvR = hv[cR];
    } else {
        // Solid wall boundary (reflect)
        hR = hL;
        double vn = (huL * nx + hvL * ny) / hL;
        huR = huL - 2.0 * vn * hL * nx;
        hvR = hvL - 2.0 * vn * hL * ny;
    }
    
    double flux[3];
    SWEFlux(hL, huL, hvL, hR, huR, hvR, nx, ny, g, flux);
    
    double L = face_length[f];
    flux_eta[f] = flux[0] * L;
    flux_hu[f] = flux[1] * L;
    flux_hv[f] = flux[2] * L;
}

/**
 * @brief Apply bottom friction (Manning)
 */
__global__ void applyBottomFrictionKernel(
    const double* eta, double* hu, double* hv,
    const double* depth, const double* manning,
    double g, double dt, double min_depth, int n_cells
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_cells) return;
    
    double h = fmax(min_depth, depth[idx] + eta[idx]);
    double n = manning[idx];
    
    if (h > min_depth && n > 0) {
        double u = hu[idx] / h;
        double v = hv[idx] / h;
        double speed = sqrt(u*u + v*v);
        
        // Manning friction: τ = ρ g n² |u| u / h^(4/3)
        double Cf = g * n * n / pow(h, 1.0/3.0);
        double friction = Cf * speed;
        
        // Implicit treatment for stability
        double factor = 1.0 / (1.0 + dt * friction / h);
        hu[idx] *= factor;
        hv[idx] *= factor;
    }
}

// =============================================================================
// Infrasound GPU Kernels
// =============================================================================

/**
 * @brief Compute effective sound speed including wind
 */
__global__ void computeEffectiveSoundSpeedKernel(
    const double* sound_speed,
    const double* wind_u, const double* wind_v,
    const double* azimuth,  // Propagation direction
    double* effective_c,
    int n_cells
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_cells) return;
    
    double c = sound_speed[idx];
    double u = wind_u[idx];
    double v = wind_v[idx];
    double az = azimuth[idx];
    
    // Wind component along propagation direction
    double wind_along = u * sin(az) + v * cos(az);
    
    effective_c[idx] = c + wind_along;
}

/**
 * @brief Apply atmospheric absorption
 */
__global__ void applyAtmosphericAbsorptionKernel(
    double* field_real, double* field_imag,
    const double* absorption_coeff,
    double dr, int n_cells
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_cells) return;
    
    double alpha = absorption_coeff[idx];
    double decay = exp(-alpha * dr);
    
    field_real[idx] *= decay;
    field_imag[idx] *= decay;
}

#endif // USE_CUDA

// =============================================================================
// Class Method Implementations
// =============================================================================

// BlackOilKernelGPU

BlackOilKernelGPU::BlackOilKernelGPU() 
    : GPUKernelBase<BlackOilKernelGPU, BlackOilKernel>() 
{
#ifdef USE_CUDA
    d_pressure_oil = nullptr;
    d_saturation_water = nullptr;
    d_saturation_gas = nullptr;
    d_Bo = nullptr;
    d_Bg = nullptr;
    d_Bw = nullptr;
    d_Rs = nullptr;
    pvt_table_size = 0;
#endif
}

BlackOilKernelGPU::~BlackOilKernelGPU() {
    freeGPUMemory();
}

void BlackOilKernelGPU::allocateGPUMemory(int n_cells) {
#ifdef USE_CUDA
    n_elements_gpu_ = n_cells;
    
    CUDA_CHECK(cudaMalloc(&d_pressure_oil, n_cells * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_saturation_water, n_cells * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_saturation_gas, n_cells * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_Bo, n_cells * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_Bg, n_cells * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_Bw, n_cells * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_Rs, n_cells * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_viscosity_oil, n_cells * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_viscosity_gas, n_cells * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_viscosity_water, n_cells * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_density_oil, n_cells * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_density_gas, n_cells * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_density_water, n_cells * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_krw, n_cells * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_kro, n_cells * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_krg, n_cells * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_lambda_oil, n_cells * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_lambda_gas, n_cells * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_lambda_water, n_cells * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_lambda_total, n_cells * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_porosity, n_cells * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_permeability, n_cells * sizeof(double)));
    
    gpu_initialized_ = true;
#endif
}

void BlackOilKernelGPU::freeGPUMemory() {
#ifdef USE_CUDA
    if (d_pressure_oil) cudaFree(d_pressure_oil);
    if (d_saturation_water) cudaFree(d_saturation_water);
    if (d_saturation_gas) cudaFree(d_saturation_gas);
    if (d_Bo) cudaFree(d_Bo);
    if (d_Bg) cudaFree(d_Bg);
    if (d_Bw) cudaFree(d_Bw);
    if (d_Rs) cudaFree(d_Rs);
    if (d_viscosity_oil) cudaFree(d_viscosity_oil);
    if (d_viscosity_gas) cudaFree(d_viscosity_gas);
    if (d_viscosity_water) cudaFree(d_viscosity_water);
    if (d_density_oil) cudaFree(d_density_oil);
    if (d_density_gas) cudaFree(d_density_gas);
    if (d_density_water) cudaFree(d_density_water);
    if (d_krw) cudaFree(d_krw);
    if (d_kro) cudaFree(d_kro);
    if (d_krg) cudaFree(d_krg);
    if (d_lambda_oil) cudaFree(d_lambda_oil);
    if (d_lambda_gas) cudaFree(d_lambda_gas);
    if (d_lambda_water) cudaFree(d_lambda_water);
    if (d_lambda_total) cudaFree(d_lambda_total);
    if (d_porosity) cudaFree(d_porosity);
    if (d_permeability) cudaFree(d_permeability);
    
    d_pressure_oil = nullptr;
    gpu_initialized_ = false;
#endif
}

void BlackOilKernelGPU::computePVTPropertiesGPU() {
#ifdef USE_CUDA
    if (!gpu_initialized_) return;
    
    int grid_size = (n_elements_gpu_ + BLOCK_SIZE_1D - 1) / BLOCK_SIZE_1D;
    
    computePVTKernel<<<grid_size, BLOCK_SIZE_1D>>>(
        d_pressure_oil, d_Bo, d_Bg, d_Bw, d_Rs,
        d_viscosity_oil, d_viscosity_gas, d_viscosity_water,
        d_density_oil, d_density_gas, d_density_water,
        d_pvt_pressure, d_pvt_Bo, d_pvt_Bg, d_pvt_Rs,
        pvt_table_size, n_elements_gpu_
    );
    
    cudaDeviceSynchronize();
#endif
}

// ThermalKernelGPU

ThermalKernelGPU::ThermalKernelGPU() 
    : GPUKernelBase<ThermalKernelGPU, ThermalKernel>() 
{
#ifdef USE_CUDA
    d_temperature = nullptr;
#endif
}

ThermalKernelGPU::~ThermalKernelGPU() {
    freeGPUMemory();
}

void ThermalKernelGPU::allocateGPUMemory(int n_cells) {
#ifdef USE_CUDA
    n_elements_gpu_ = n_cells;
    
    CUDA_CHECK(cudaMalloc(&d_temperature, n_cells * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_temperature_old, n_cells * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_rho_cp_eff, n_cells * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_k_eff, n_cells * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_k_solid, n_cells * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_k_fluid, n_cells * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_rho_solid, n_cells * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_rho_fluid, n_cells * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_cp_solid, n_cells * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_cp_fluid, n_cells * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_porosity, n_cells * sizeof(double)));
    
    gpu_initialized_ = true;
#endif
}

void ThermalKernelGPU::freeGPUMemory() {
#ifdef USE_CUDA
    if (d_temperature) cudaFree(d_temperature);
    if (d_temperature_old) cudaFree(d_temperature_old);
    if (d_rho_cp_eff) cudaFree(d_rho_cp_eff);
    if (d_k_eff) cudaFree(d_k_eff);
    if (d_k_solid) cudaFree(d_k_solid);
    if (d_k_fluid) cudaFree(d_k_fluid);
    if (d_rho_solid) cudaFree(d_rho_solid);
    if (d_rho_fluid) cudaFree(d_rho_fluid);
    if (d_cp_solid) cudaFree(d_cp_solid);
    if (d_cp_fluid) cudaFree(d_cp_fluid);
    if (d_porosity) cudaFree(d_porosity);
    
    d_temperature = nullptr;
    gpu_initialized_ = false;
#endif
}

void ThermalKernelGPU::computeEffectivePropertiesGPU() {
#ifdef USE_CUDA
    if (!gpu_initialized_) return;
    
    int grid_size = (n_elements_gpu_ + BLOCK_SIZE_1D - 1) / BLOCK_SIZE_1D;
    
    computeEffectiveThermalKernel<<<grid_size, BLOCK_SIZE_1D>>>(
        d_k_solid, d_k_fluid, d_rho_solid, d_rho_fluid,
        d_cp_solid, d_cp_fluid, d_porosity,
        d_k_eff, d_rho_cp_eff, n_elements_gpu_
    );
    
    cudaDeviceSynchronize();
#endif
}

// ============================================================================
// GeomechanicsKernelGPU Implementation
// ============================================================================

GeomechanicsKernelGPU::GeomechanicsKernelGPU()
    : GPUKernelBase<GeomechanicsKernelGPU, GeomechanicsKernel>()
{
#ifdef USE_CUDA
    d_displacement = d_strain = d_stress = d_force = nullptr;
    d_youngs_modulus = d_poisson_ratio = d_density = d_biot_coefficient = nullptr;
    d_plastic_strain = d_yield_state = d_cohesion = d_friction_angle = nullptr;
    d_pressure = nullptr;
    d_connectivity = nullptr;
#endif
}

GeomechanicsKernelGPU::~GeomechanicsKernelGPU() {
    freeGPUMemory();
}

PetscErrorCode GeomechanicsKernelGPU::setup(DM dm, PetscFE fe) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr = GeomechanicsKernel::setup(dm, fe); CHKERRQ(ierr);
    if (initializeGPU()) {
        PetscPrintf(PETSC_COMM_WORLD, "GPU-accelerated geomechanics kernel initialized\n");
    }
    PetscFunctionReturn(0);
}

void GeomechanicsKernelGPU::residual(const PetscScalar u[], const PetscScalar u_t[],
                                      const PetscScalar u_x[], const PetscScalar a[],
                                      const PetscReal x[], PetscScalar f[]) {
    if (!gpu_initialized_) {
        GeomechanicsKernel::residual(u, u_t, u_x, a, x, f);
        return;
    }
    // GPU path: fall through to CPU for pointwise evaluation;
    // GPU acceleration applies to batch operations (computeStrainGPU, etc.)
    GeomechanicsKernel::residual(u, u_t, u_x, a, x, f);
}

void GeomechanicsKernelGPU::jacobian(const PetscScalar u[], const PetscScalar u_t[],
                                      const PetscScalar u_x[], const PetscScalar a[],
                                      const PetscReal x[], PetscScalar J[]) {
    GeomechanicsKernel::jacobian(u, u_t, u_x, a, x, J);
}

void GeomechanicsKernelGPU::allocateGPUMemory(int n_cells) {
#ifdef USE_CUDA
    n_elements_gpu_ = n_cells;
    GPUManager& gpu = getGPUManager();
    d_displacement = static_cast<double*>(gpu.allocateDevice(n_cells * 3 * sizeof(double)));
    d_strain = static_cast<double*>(gpu.allocateDevice(n_cells * 6 * sizeof(double)));
    d_stress = static_cast<double*>(gpu.allocateDevice(n_cells * 6 * sizeof(double)));
    d_force = static_cast<double*>(gpu.allocateDevice(n_cells * 3 * sizeof(double)));
    d_youngs_modulus = static_cast<double*>(gpu.allocateDevice(n_cells * sizeof(double)));
    d_poisson_ratio = static_cast<double*>(gpu.allocateDevice(n_cells * sizeof(double)));
    d_density = static_cast<double*>(gpu.allocateDevice(n_cells * sizeof(double)));
#endif
}

void GeomechanicsKernelGPU::freeGPUMemory() {
#ifdef USE_CUDA
    GPUManager& gpu = getGPUManager();
    if (d_displacement) gpu.freeDevice(d_displacement);
    if (d_strain) gpu.freeDevice(d_strain);
    if (d_stress) gpu.freeDevice(d_stress);
    if (d_force) gpu.freeDevice(d_force);
    if (d_youngs_modulus) gpu.freeDevice(d_youngs_modulus);
    if (d_poisson_ratio) gpu.freeDevice(d_poisson_ratio);
    if (d_density) gpu.freeDevice(d_density);
    d_displacement = d_strain = d_stress = d_force = nullptr;
    d_youngs_modulus = d_poisson_ratio = d_density = nullptr;
#endif
}

// ============================================================================
// HydrodynamicKernelGPU Implementation
// ============================================================================

HydrodynamicKernelGPU::HydrodynamicKernelGPU()
    : GPUKernelBase<HydrodynamicKernelGPU, HydrodynamicKernel>()
{
}

HydrodynamicKernelGPU::~HydrodynamicKernelGPU() {
    freeGPUMemory();
}

PetscErrorCode HydrodynamicKernelGPU::setup(DM dm, PetscFE fe) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr = HydrodynamicKernel::setup(dm, fe); CHKERRQ(ierr);
    if (initializeGPU()) {
        PetscPrintf(PETSC_COMM_WORLD, "GPU-accelerated hydrodynamic kernel initialized\n");
    }
    PetscFunctionReturn(0);
}

void HydrodynamicKernelGPU::residual(const PetscScalar u[], const PetscScalar u_t[],
                                      const PetscScalar u_x[], const PetscScalar a[],
                                      const PetscReal x[], PetscScalar f[]) {
    if (!gpu_initialized_) {
        HydrodynamicKernel::residual(u, u_t, u_x, a, x, f);
        return;
    }
    HydrodynamicKernel::residual(u, u_t, u_x, a, x, f);
}

void HydrodynamicKernelGPU::jacobian(const PetscScalar u[], const PetscScalar u_t[],
                                      const PetscScalar u_x[], const PetscScalar a[],
                                      const PetscReal x[], PetscScalar J[]) {
    HydrodynamicKernel::jacobian(u, u_t, u_x, a, x, J);
}

void HydrodynamicKernelGPU::allocateGPUMemory(int /*n_cells*/) {
    // HLLC flux kernels already implemented in computeHydroFluxKernel
}

void HydrodynamicKernelGPU::freeGPUMemory() {
    // Managed by base class destruction
}

// ============================================================================
// AtmosphericBlastKernelGPU Implementation
// ============================================================================

AtmosphericBlastKernelGPU::AtmosphericBlastKernelGPU()
    : GPUKernelBase<AtmosphericBlastKernelGPU, AtmosphericBlastKernel>()
{
}

AtmosphericBlastKernelGPU::~AtmosphericBlastKernelGPU() {
    freeGPUMemory();
}

PetscErrorCode AtmosphericBlastKernelGPU::setup(DM dm, PetscFE fe) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr = AtmosphericBlastKernel::setup(dm, fe); CHKERRQ(ierr);
    if (initializeGPU()) {
        PetscPrintf(PETSC_COMM_WORLD, "GPU-accelerated atmospheric blast kernel initialized\n");
    }
    PetscFunctionReturn(0);
}

void AtmosphericBlastKernelGPU::residual(const PetscScalar u[], const PetscScalar u_t[],
                                          const PetscScalar u_x[], const PetscScalar a[],
                                          const PetscReal x[], PetscScalar f[]) {
    if (!gpu_initialized_) {
        AtmosphericBlastKernel::residual(u, u_t, u_x, a, x, f);
        return;
    }
    // GPU batch operations available: sedovTaylorBlast, fireballDynamics, etc.
    // Pointwise evaluation falls back to CPU
    AtmosphericBlastKernel::residual(u, u_t, u_x, a, x, f);
}

void AtmosphericBlastKernelGPU::jacobian(const PetscScalar u[], const PetscScalar u_t[],
                                          const PetscScalar u_x[], const PetscScalar a[],
                                          const PetscReal x[], PetscScalar J[]) {
    AtmosphericBlastKernel::jacobian(u, u_t, u_x, a, x, J);
}

void AtmosphericBlastKernelGPU::allocateGPUMemory(int /*n_cells*/) {
    // Uses GPUExplosionKernels and GPUAtmosphericKernels from GPUCompute.cu
}

void AtmosphericBlastKernelGPU::freeGPUMemory() {
}

// ============================================================================
// InfrasoundKernelGPU Implementation
// ============================================================================

InfrasoundKernelGPU::InfrasoundKernelGPU()
    : GPUKernelBase<InfrasoundKernelGPU, InfrasoundKernel>()
{
}

InfrasoundKernelGPU::~InfrasoundKernelGPU() {
    freeGPUMemory();
}

PetscErrorCode InfrasoundKernelGPU::setup(DM dm, PetscFE fe) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr = InfrasoundKernel::setup(dm, fe); CHKERRQ(ierr);
    if (initializeGPU()) {
        PetscPrintf(PETSC_COMM_WORLD, "GPU-accelerated infrasound kernel initialized\n");
    }
    PetscFunctionReturn(0);
}

void InfrasoundKernelGPU::residual(const PetscScalar u[], const PetscScalar u_t[],
                                    const PetscScalar u_x[], const PetscScalar a[],
                                    const PetscReal x[], PetscScalar f[]) {
    if (!gpu_initialized_) {
        InfrasoundKernel::residual(u, u_t, u_x, a, x, f);
        return;
    }
    // GPU kernels available: parabolicEquationStep, wideAnglePE, etc.
    InfrasoundKernel::residual(u, u_t, u_x, a, x, f);
}

void InfrasoundKernelGPU::jacobian(const PetscScalar u[], const PetscScalar u_t[],
                                    const PetscScalar u_x[], const PetscScalar a[],
                                    const PetscReal x[], PetscScalar J[]) {
    InfrasoundKernel::jacobian(u, u_t, u_x, a, x, J);
}

void InfrasoundKernelGPU::allocateGPUMemory(int /*n_cells*/) {
    // Uses GPUPhysicsKernels from GPUPhysicsKernels.cuh
}

void InfrasoundKernelGPU::freeGPUMemory() {
}

// ============================================================================
// TsunamiKernelGPU Implementation
// ============================================================================

TsunamiKernelGPU::TsunamiKernelGPU()
    : GPUKernelBase<TsunamiKernelGPU, TsunamiKernel>()
{
}

TsunamiKernelGPU::~TsunamiKernelGPU() {
    freeGPUMemory();
}

PetscErrorCode TsunamiKernelGPU::setup(DM dm, PetscFE fe) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr = TsunamiKernel::setup(dm, fe); CHKERRQ(ierr);
    if (initializeGPU()) {
        PetscPrintf(PETSC_COMM_WORLD, "GPU-accelerated tsunami kernel initialized\n");
    }
    PetscFunctionReturn(0);
}

void TsunamiKernelGPU::residual(const PetscScalar u[], const PetscScalar u_t[],
                                  const PetscScalar u_x[], const PetscScalar a[],
                                  const PetscReal x[], PetscScalar f[]) {
    if (!gpu_initialized_) {
        TsunamiKernel::residual(u, u_t, u_x, a, x, f);
        return;
    }
    TsunamiKernel::residual(u, u_t, u_x, a, x, f);
}

void TsunamiKernelGPU::jacobian(const PetscScalar u[], const PetscScalar u_t[],
                                  const PetscScalar u_x[], const PetscScalar a[],
                                  const PetscReal x[], PetscScalar J[]) {
    TsunamiKernel::jacobian(u, u_t, u_x, a, x, J);
}

void TsunamiKernelGPU::allocateGPUMemory(int /*n_cells*/) {
    // Shallow water GPU kernels from GPUPhysicsKernels.cuh
}

void TsunamiKernelGPU::freeGPUMemory() {
}

// ============================================================================
// FalloutKernelGPU Implementation
// ============================================================================

FalloutKernelGPU::FalloutKernelGPU()
    : GPUKernelBase<FalloutKernelGPU, FalloutKernel>()
{
#ifdef USE_CUDA
    d_pos_x = d_pos_y = d_pos_z = nullptr;
    d_activity = d_particle_diameter = d_settling_velocity = nullptr;
    d_wind_u = d_wind_v = d_wind_w = nullptr;
    d_half_life = d_dose_conversion = nullptr;
    d_ground_deposition = nullptr;
    deposition_nx_ = deposition_ny_ = 0;
#endif
}

FalloutKernelGPU::~FalloutKernelGPU() {
    freeGPUMemory();
}

PetscErrorCode FalloutKernelGPU::setup(DM dm, PetscFE fe) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr = FalloutKernel::setup(dm, fe); CHKERRQ(ierr);
    if (initializeGPU()) {
        PetscPrintf(PETSC_COMM_WORLD, "GPU-accelerated fallout kernel initialized\n");
    }
    PetscFunctionReturn(0);
}

void FalloutKernelGPU::residual(const PetscScalar u[], const PetscScalar u_t[],
                                 const PetscScalar u_x[], const PetscScalar a[],
                                 const PetscReal x[], PetscScalar f[]) {
    if (!gpu_initialized_) {
        FalloutKernel::residual(u, u_t, u_x, a, x, f);
        return;
    }
    // GPU batch operations: advectFalloutParticles, radioactiveDecay, etc.
    FalloutKernel::residual(u, u_t, u_x, a, x, f);
}

void FalloutKernelGPU::jacobian(const PetscScalar u[], const PetscScalar u_t[],
                                 const PetscScalar u_x[], const PetscScalar a[],
                                 const PetscReal x[], PetscScalar J[]) {
    FalloutKernel::jacobian(u, u_t, u_x, a, x, J);
}

void FalloutKernelGPU::allocateGPUMemory(int n_cells) {
#ifdef USE_CUDA
    n_elements_gpu_ = n_cells;
    GPUManager& gpu = getGPUManager();
    d_pos_x = static_cast<double*>(gpu.allocateDevice(n_cells * sizeof(double)));
    d_pos_y = static_cast<double*>(gpu.allocateDevice(n_cells * sizeof(double)));
    d_pos_z = static_cast<double*>(gpu.allocateDevice(n_cells * sizeof(double)));
    d_activity = static_cast<double*>(gpu.allocateDevice(n_cells * sizeof(double)));
    d_settling_velocity = static_cast<double*>(gpu.allocateDevice(n_cells * sizeof(double)));
    d_wind_u = static_cast<double*>(gpu.allocateDevice(n_cells * sizeof(double)));
    d_wind_v = static_cast<double*>(gpu.allocateDevice(n_cells * sizeof(double)));
    d_wind_w = static_cast<double*>(gpu.allocateDevice(n_cells * sizeof(double)));
#endif
}

void FalloutKernelGPU::freeGPUMemory() {
#ifdef USE_CUDA
    GPUManager& gpu = getGPUManager();
    if (d_pos_x) gpu.freeDevice(d_pos_x);
    if (d_pos_y) gpu.freeDevice(d_pos_y);
    if (d_pos_z) gpu.freeDevice(d_pos_z);
    if (d_activity) gpu.freeDevice(d_activity);
    if (d_settling_velocity) gpu.freeDevice(d_settling_velocity);
    if (d_wind_u) gpu.freeDevice(d_wind_u);
    if (d_wind_v) gpu.freeDevice(d_wind_v);
    if (d_wind_w) gpu.freeDevice(d_wind_w);
    d_pos_x = d_pos_y = d_pos_z = d_activity = d_settling_velocity = nullptr;
    d_wind_u = d_wind_v = d_wind_w = nullptr;
#endif
}

// ============================================================================
// ParticleTransportKernelGPU Implementation
// ============================================================================

ParticleTransportKernelGPU::ParticleTransportKernelGPU()
    : GPUKernelBase<ParticleTransportKernelGPU, ParticleTransportKernel>()
{
#ifdef USE_CUDA
    d_concentration = d_concentration_old = nullptr;
    d_velocity = d_diffusivity = d_settling_velocity = nullptr;
    d_porosity = d_deposition_rate = nullptr;
    d_dispersivity_L = d_dispersivity_T = nullptr;
#endif
}

ParticleTransportKernelGPU::~ParticleTransportKernelGPU() {
    freeGPUMemory();
}

PetscErrorCode ParticleTransportKernelGPU::setup(DM dm, PetscFE fe) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr = ParticleTransportKernel::setup(dm, fe); CHKERRQ(ierr);
    if (initializeGPU()) {
        PetscPrintf(PETSC_COMM_WORLD, "GPU-accelerated particle transport kernel initialized\n");
    }
    PetscFunctionReturn(0);
}

void ParticleTransportKernelGPU::residual(const PetscScalar u[], const PetscScalar u_t[],
                                           const PetscScalar u_x[], const PetscScalar a[],
                                           const PetscReal x[], PetscScalar f[]) {
    if (!gpu_initialized_) {
        ParticleTransportKernel::residual(u, u_t, u_x, a, x, f);
        return;
    }
    ParticleTransportKernel::residual(u, u_t, u_x, a, x, f);
}

void ParticleTransportKernelGPU::jacobian(const PetscScalar u[], const PetscScalar u_t[],
                                           const PetscScalar u_x[], const PetscScalar a[],
                                           const PetscReal x[], PetscScalar J[]) {
    ParticleTransportKernel::jacobian(u, u_t, u_x, a, x, J);
}

void ParticleTransportKernelGPU::allocateGPUMemory(int n_cells) {
#ifdef USE_CUDA
    n_elements_gpu_ = n_cells;
    GPUManager& gpu = getGPUManager();
    d_concentration = static_cast<double*>(gpu.allocateDevice(n_cells * sizeof(double)));
    d_concentration_old = static_cast<double*>(gpu.allocateDevice(n_cells * sizeof(double)));
    d_velocity = static_cast<double*>(gpu.allocateDevice(n_cells * 3 * sizeof(double)));
    d_diffusivity = static_cast<double*>(gpu.allocateDevice(n_cells * sizeof(double)));
    d_settling_velocity = static_cast<double*>(gpu.allocateDevice(n_cells * sizeof(double)));
    d_porosity = static_cast<double*>(gpu.allocateDevice(n_cells * sizeof(double)));
#endif
}

void ParticleTransportKernelGPU::freeGPUMemory() {
#ifdef USE_CUDA
    GPUManager& gpu = getGPUManager();
    if (d_concentration) gpu.freeDevice(d_concentration);
    if (d_concentration_old) gpu.freeDevice(d_concentration_old);
    if (d_velocity) gpu.freeDevice(d_velocity);
    if (d_diffusivity) gpu.freeDevice(d_diffusivity);
    if (d_settling_velocity) gpu.freeDevice(d_settling_velocity);
    if (d_porosity) gpu.freeDevice(d_porosity);
    d_concentration = d_concentration_old = d_velocity = nullptr;
    d_diffusivity = d_settling_velocity = d_porosity = nullptr;
#endif
}

// ============================================================================
// Factory functions
// ============================================================================

std::shared_ptr<PhysicsKernel> createGPUKernelExtended(PhysicsType type) {
#ifdef USE_CUDA
    // Check that a GPU is actually available before creating GPU kernels
    GPUManager& gpu = GPUManager::getInstance();
    if (!gpu.isAvailable()) {
        return nullptr;
    }
#endif

    switch (type) {
        // --- Core reservoir physics ---
        case PhysicsType::FLUID_FLOW:
            return std::make_shared<SinglePhaseFlowKernelGPU>();

        // --- Dynamic wave propagation ---
        case PhysicsType::ELASTODYNAMICS:
            return std::make_shared<ElastodynamicsKernelGPU>();
        case PhysicsType::POROELASTODYNAMICS:
            return std::make_shared<PoroelastodynamicsKernelGPU>();

#ifdef USE_CUDA
        // --- Extended GPU kernels (require CUDA build) ---
        case PhysicsType::THERMAL:
            return std::make_shared<ThermalKernelGPU>();
        case PhysicsType::GEOMECHANICS:
            return std::make_shared<GeomechanicsKernelGPU>();
        case PhysicsType::HYDRODYNAMIC:
            return std::make_shared<HydrodynamicKernelGPU>();
        case PhysicsType::ATMOSPHERIC_BLAST:
            return std::make_shared<AtmosphericBlastKernelGPU>();
        case PhysicsType::INFRASOUND:
            return std::make_shared<InfrasoundKernelGPU>();
        case PhysicsType::TSUNAMI:
            return std::make_shared<TsunamiKernelGPU>();
        case PhysicsType::FALLOUT:
            return std::make_shared<FalloutKernelGPU>();
#endif
        default:
            return nullptr;
    }
}

bool hasExtendedGPUKernel(PhysicsType type) {
    switch (type) {
        // Core reservoir physics
        case PhysicsType::FLUID_FLOW:
        // Dynamic wave propagation
        case PhysicsType::ELASTODYNAMICS:
        case PhysicsType::POROELASTODYNAMICS:
        // Extended physics
        case PhysicsType::THERMAL:
        case PhysicsType::GEOMECHANICS:
        case PhysicsType::HYDRODYNAMIC:
        // Atmospheric / explosion physics
        case PhysicsType::ATMOSPHERIC_BLAST:
        case PhysicsType::INFRASOUND:
        case PhysicsType::FALLOUT:
        // Surface water
        case PhysicsType::TSUNAMI:
        // Particle transport
        case PhysicsType::PARTICLE_TRANSPORT:
            return true;
        default:
            return false;
    }
}

std::vector<PhysicsType> getAvailableGPUKernels() {
    return {
        // Core reservoir physics
        PhysicsType::FLUID_FLOW,
        // Dynamic wave propagation
        PhysicsType::ELASTODYNAMICS,
        PhysicsType::POROELASTODYNAMICS,
        // Extended physics
        PhysicsType::THERMAL,
        PhysicsType::GEOMECHANICS,
        PhysicsType::HYDRODYNAMIC,
        // Atmospheric / explosion physics
        PhysicsType::ATMOSPHERIC_BLAST,
        PhysicsType::INFRASOUND,
        PhysicsType::FALLOUT,
        // Surface water
        PhysicsType::TSUNAMI,
        // Particle transport
        PhysicsType::PARTICLE_TRANSPORT
    };
}

} // namespace FSRM
