/**
 * @file BoundaryConditions.cpp
 * @brief Implementation of boundary conditions for seismic wave propagation
 */

#include "BoundaryConditions.hpp"
#include "DiscontinuousGalerkin.hpp"
#include <cmath>
#include <algorithm>
#include <complex>
#include <stdexcept>
#include <iostream>

namespace FSRM {

// =============================================================================
// FreeSurfaceBC Implementation
// =============================================================================

FreeSurfaceBC::FreeSurfaceBC()
    : use_characteristic_decomposition(true) {}

void FreeSurfaceBC::applyToResidual(const BoundaryFace& face,
                                    const double* u, double* residual,
                                    double t) const {
    // Free surface: enforce σ·n = 0
    // For velocity-stress formulation, the stress components are set
    // to cancel the incoming wave (reflection)
    
    const double* n = face.normal.data();
    
    // Stress components: σxx, σyy, σzz, σxy, σxz, σyz (indices 3-8 in u)
    // Velocity components: vx, vy, vz (indices 0-2 in u)
    
    // Compute traction t = σ·n
    double traction[3];
    computeSurfaceTraction(&u[3], n, traction);
    
    // Apply zero traction constraint
    // This modifies the residual to enforce t = 0
    double penalty = 1e6;  // Large penalty for enforcement
    
    residual[0] -= penalty * traction[0] * n[0];
    residual[1] -= penalty * traction[1] * n[1];
    residual[2] -= penalty * traction[2] * n[2];
}

void FreeSurfaceBC::computeBoundaryFlux(const BoundaryFace& face,
                                        const double* u_int, double* flux_n,
                                        double t) const {
    // For free surface, compute flux using ghost state
    double u_ext[9];
    getBoundaryState(face, u_int, u_ext, t);
    
    const double* n = face.normal.data();
    
    // Average flux (central flux for free surface)
    // Using velocity-stress formulation:
    // F·n = [σ·n; v⊗n] (simplified)
    
    double rho = face.rho;
    double lambda = rho * (face.vp * face.vp - 2.0 * face.vs * face.vs);
    double mu = rho * face.vs * face.vs;
    
    // Interior velocity and stress
    double vx_int = u_int[0], vy_int = u_int[1], vz_int = u_int[2];
    double sxx_int = u_int[3], syy_int = u_int[4], szz_int = u_int[5];
    double sxy_int = u_int[6], sxz_int = u_int[7], syz_int = u_int[8];
    
    // Exterior (ghost) velocity and stress
    double vx_ext = u_ext[0], vy_ext = u_ext[1], vz_ext = u_ext[2];
    double sxx_ext = u_ext[3], syy_ext = u_ext[4], szz_ext = u_ext[5];
    double sxy_ext = u_ext[6], sxz_ext = u_ext[7], syz_ext = u_ext[8];
    
    // Average states
    double vx = 0.5 * (vx_int + vx_ext);
    double vy = 0.5 * (vy_int + vy_ext);
    double vz = 0.5 * (vz_int + vz_ext);
    
    double sxx = 0.5 * (sxx_int + sxx_ext);
    double syy = 0.5 * (syy_int + syy_ext);
    double szz = 0.5 * (szz_int + szz_ext);
    double sxy = 0.5 * (sxy_int + sxy_ext);
    double sxz = 0.5 * (sxz_int + sxz_ext);
    double syz = 0.5 * (syz_int + syz_ext);
    
    // Flux for velocity equations: F_v = -σ·n / ρ
    double tn_x = sxx * n[0] + sxy * n[1] + sxz * n[2];
    double tn_y = sxy * n[0] + syy * n[1] + syz * n[2];
    double tn_z = sxz * n[0] + syz * n[1] + szz * n[2];
    
    flux_n[0] = -tn_x / rho;
    flux_n[1] = -tn_y / rho;
    flux_n[2] = -tn_z / rho;
    
    // Flux for stress equations: F_σ = -(λ(∇·v)I + μ(∇v + ∇v^T))
    double vn = vx * n[0] + vy * n[1] + vz * n[2];
    
    flux_n[3] = -(lambda + 2.0 * mu) * vn * n[0] * n[0] - lambda * vn;
    flux_n[4] = -(lambda + 2.0 * mu) * vn * n[1] * n[1] - lambda * vn;
    flux_n[5] = -(lambda + 2.0 * mu) * vn * n[2] * n[2] - lambda * vn;
    flux_n[6] = -mu * vn * n[0] * n[1];
    flux_n[7] = -mu * vn * n[0] * n[2];
    flux_n[8] = -mu * vn * n[1] * n[2];
}

void FreeSurfaceBC::getBoundaryState(const BoundaryFace& face,
                                     const double* u_int, double* u_ext,
                                     double t) const {
    // Mirror state for free surface
    // Velocity: v_ext = v_int (symmetric)
    // Stress: σ_ext = mirror(σ_int) such that σ·n = 0
    
    const double* n = face.normal.data();
    
    // Copy velocity (symmetric reflection)
    u_ext[0] = u_int[0];
    u_ext[1] = u_int[1];
    u_ext[2] = u_int[2];
    
    // Mirror stress
    applyMirrorStress(&u_int[3], n, &u_ext[3]);
}

void FreeSurfaceBC::applyMirrorStress(const double* sigma_int, const double* n,
                                      double* sigma_ext) const {
    // Mirror stress tensor across free surface
    // σ_ext = σ_int - 2(σ_int·n)⊗n
    
    // Interior stress tensor
    double sxx = sigma_int[0], syy = sigma_int[1], szz = sigma_int[2];
    double sxy = sigma_int[3], sxz = sigma_int[4], syz = sigma_int[5];
    
    // Compute σ·n
    double tn_x = sxx * n[0] + sxy * n[1] + sxz * n[2];
    double tn_y = sxy * n[0] + syy * n[1] + syz * n[2];
    double tn_z = sxz * n[0] + syz * n[1] + szz * n[2];
    
    // Mirror: σ_ext = σ_int - 2(t⊗n)
    sigma_ext[0] = sxx - 2.0 * tn_x * n[0];  // σxx
    sigma_ext[1] = syy - 2.0 * tn_y * n[1];  // σyy
    sigma_ext[2] = szz - 2.0 * tn_z * n[2];  // σzz
    sigma_ext[3] = sxy - (tn_x * n[1] + tn_y * n[0]);  // σxy
    sigma_ext[4] = sxz - (tn_x * n[2] + tn_z * n[0]);  // σxz
    sigma_ext[5] = syz - (tn_y * n[2] + tn_z * n[1]);  // σyz
}

void FreeSurfaceBC::computeSurfaceTraction(const double* sigma, const double* n,
                                           double* traction) const {
    // t = σ·n
    double sxx = sigma[0], syy = sigma[1], szz = sigma[2];
    double sxy = sigma[3], sxz = sigma[4], syz = sigma[5];
    
    traction[0] = sxx * n[0] + sxy * n[1] + sxz * n[2];
    traction[1] = sxy * n[0] + syy * n[1] + syz * n[2];
    traction[2] = sxz * n[0] + syz * n[1] + szz * n[2];
}

// =============================================================================
// PMLBoundary Implementation
// =============================================================================

PMLBoundary::PMLBoundary()
    : thickness(1.0), damping_order(3), reflection_coef(1e-6),
      alpha_max(M_PI), kappa_max(1.0),
      x_min(0), x_max(1), y_min(0), y_max(1), z_min(0), z_max(1),
      rho_ref(2500.0), vp_ref(5000.0), vs_ref(3000.0), d_max(0) {}

void PMLBoundary::setThickness(double d) {
    thickness = d;
}

void PMLBoundary::setOrder(int m) {
    damping_order = m;
}

void PMLBoundary::setReflectionCoefficient(double R) {
    reflection_coef = R;
    // Compute optimal d_max based on reflection coefficient
    // d_max = -(n+1) * c * ln(R) / (2 * L)
    d_max = -(damping_order + 1) * vp_ref * std::log(reflection_coef) / (2.0 * thickness);
}

void PMLBoundary::setCFSParameters(double alpha, double kappa) {
    alpha_max = alpha;
    kappa_max = kappa;
}

void PMLBoundary::setDomainBounds(double xmin, double xmax,
                                  double ymin, double ymax,
                                  double zmin, double zmax) {
    x_min = xmin; x_max = xmax;
    y_min = ymin; y_max = ymax;
    z_min = zmin; z_max = zmax;
}

void PMLBoundary::setMaterialProperties(double rho, double vp, double vs) {
    rho_ref = rho;
    vp_ref = vp;
    vs_ref = vs;
    
    // Recompute d_max
    d_max = -(damping_order + 1) * vp_ref * std::log(reflection_coef) / (2.0 * thickness);
}

bool PMLBoundary::isInPML(double x, double y, double z) const {
    return (x < x_min + thickness) || (x > x_max - thickness) ||
           (y < y_min + thickness) || (y > y_max - thickness) ||
           (z < z_min + thickness) || (z > z_max - thickness);
}

int PMLBoundary::getPMLDirection(double x, double y, double z) const {
    int dir = 0;
    if (x < x_min + thickness || x > x_max - thickness) dir |= 1;
    if (y < y_min + thickness || y > y_max - thickness) dir |= 2;
    if (z < z_min + thickness || z > z_max - thickness) dir |= 4;
    return dir;
}

void PMLBoundary::getDampingCoefficients(double x, double y, double z,
                                         double& d_x, double& d_y, double& d_z) const {
    d_x = d_y = d_z = 0.0;
    
    // X direction
    if (x < x_min + thickness) {
        double dist = (x_min + thickness - x) / thickness;
        d_x = dampingProfile(dist);
    } else if (x > x_max - thickness) {
        double dist = (x - (x_max - thickness)) / thickness;
        d_x = dampingProfile(dist);
    }
    
    // Y direction
    if (y < y_min + thickness) {
        double dist = (y_min + thickness - y) / thickness;
        d_y = dampingProfile(dist);
    } else if (y > y_max - thickness) {
        double dist = (y - (y_max - thickness)) / thickness;
        d_y = dampingProfile(dist);
    }
    
    // Z direction
    if (z < z_min + thickness) {
        double dist = (z_min + thickness - z) / thickness;
        d_z = dampingProfile(dist);
    } else if (z > z_max - thickness) {
        double dist = (z - (z_max - thickness)) / thickness;
        d_z = dampingProfile(dist);
    }
}

void PMLBoundary::getStretchingFunctions(double x, double y, double z,
                                         std::complex<double>& s_x,
                                         std::complex<double>& s_y,
                                         std::complex<double>& s_z) const {
    double d_x, d_y, d_z;
    getDampingCoefficients(x, y, z, d_x, d_y, d_z);
    
    // Complex frequency-shifted stretching function
    // s(ξ) = κ(ξ) + d(ξ) / (α(ξ) + iω)
    // For time-domain: s(ξ) ≈ κ(ξ) + d(ξ)/α(ξ) when α >> ω
    
    double dist_x = 0, dist_y = 0, dist_z = 0;
    
    if (x < x_min + thickness) dist_x = (x_min + thickness - x) / thickness;
    else if (x > x_max - thickness) dist_x = (x - (x_max - thickness)) / thickness;
    
    if (y < y_min + thickness) dist_y = (y_min + thickness - y) / thickness;
    else if (y > y_max - thickness) dist_y = (y - (y_max - thickness)) / thickness;
    
    if (z < z_min + thickness) dist_z = (z_min + thickness - z) / thickness;
    else if (z > z_max - thickness) dist_z = (z - (z_max - thickness)) / thickness;
    
    double kappa_x = kappaProfile(dist_x);
    double kappa_y = kappaProfile(dist_y);
    double kappa_z = kappaProfile(dist_z);
    
    double alpha_x = alphaProfile(dist_x);
    double alpha_y = alphaProfile(dist_y);
    double alpha_z = alphaProfile(dist_z);
    
    // Using ω = 1 for simplicity (would be actual angular frequency)
    std::complex<double> i(0, 1);
    double omega = 2.0 * M_PI;  // Nominal frequency
    
    s_x = kappa_x + d_x / (alpha_x + i * omega);
    s_y = kappa_y + d_y / (alpha_y + i * omega);
    s_z = kappa_z + d_z / (alpha_z + i * omega);
}

double PMLBoundary::dampingProfile(double dist) const {
    // Polynomial damping profile
    // d(ξ) = d_max * (ξ/L)^n
    dist = std::max(0.0, std::min(1.0, dist));
    return d_max * std::pow(dist, damping_order);
}

double PMLBoundary::alphaProfile(double dist) const {
    // Linear decrease of α from α_max to 0
    dist = std::max(0.0, std::min(1.0, dist));
    return alpha_max * (1.0 - dist);
}

double PMLBoundary::kappaProfile(double dist) const {
    // κ increases from 1 to κ_max
    dist = std::max(0.0, std::min(1.0, dist));
    return 1.0 + (kappa_max - 1.0) * std::pow(dist, damping_order);
}

void PMLBoundary::applyToResidual(const BoundaryFace& face,
                                  const double* u, double* residual,
                                  double t) const {
    // PML is applied to volume, not just boundary
    // This function handles the outer boundary of the PML
    
    // At the outer PML boundary, use simple absorbing BC
    // (Lysmer-Kuhlemeyer style dashpot)
    
    const double* n = face.normal.data();
    double rho = face.rho > 0 ? face.rho : rho_ref;
    double vp = face.vp > 0 ? face.vp : vp_ref;
    double vs = face.vs > 0 ? face.vs : vs_ref;
    
    // Velocity at boundary
    double vx = u[0], vy = u[1], vz = u[2];
    
    // Normal and tangential velocities
    double vn = vx * n[0] + vy * n[1] + vz * n[2];
    double vtx = vx - vn * n[0];
    double vty = vy - vn * n[1];
    double vtz = vz - vn * n[2];
    
    // Absorbing traction: t = -ρ*vp*vn - ρ*vs*vt
    double tn = -rho * vp * vn;
    double ttx = -rho * vs * vtx;
    double tty = -rho * vs * vty;
    double ttz = -rho * vs * vtz;
    
    // Add to residual
    residual[0] += (tn * n[0] + ttx) / rho;
    residual[1] += (tn * n[1] + tty) / rho;
    residual[2] += (tn * n[2] + ttz) / rho;
}

void PMLBoundary::computeBoundaryFlux(const BoundaryFace& face,
                                      const double* u_int, double* flux_n,
                                      double t) const {
    double u_ext[9];
    getBoundaryState(face, u_int, u_ext, t);
    
    // Use Riemann solver for flux
    // For simplicity, use central flux
    for (int i = 0; i < 9; ++i) {
        flux_n[i] = 0.5 * (u_int[i] + u_ext[i]);
    }
}

void PMLBoundary::getBoundaryState(const BoundaryFace& face,
                                   const double* u_int, double* u_ext,
                                   double t) const {
    // Characteristic-based outgoing wave condition
    const double* n = face.normal.data();
    double rho = face.rho > 0 ? face.rho : rho_ref;
    double vp = face.vp > 0 ? face.vp : vp_ref;
    double vs = face.vs > 0 ? face.vs : vs_ref;
    
    // Interior state
    double vx_int = u_int[0], vy_int = u_int[1], vz_int = u_int[2];
    
    // Normal velocity
    double vn = vx_int * n[0] + vy_int * n[1] + vz_int * n[2];
    
    // For outgoing wave: set exterior state to absorb
    u_ext[0] = 0;
    u_ext[1] = 0;
    u_ext[2] = 0;
    
    // Zero stress in ghost region
    for (int i = 3; i < 9; ++i) {
        u_ext[i] = 0;
    }
}

int PMLBoundary::getNumStateVariables() const {
    // For split-field PML with 9 field variables (3 velocity + 6 stress)
    // Need 2 auxiliary fields per direction per variable
    // Total: 9 * 3 * 2 = 54 (for 3D)
    return 54;
}

void PMLBoundary::initializeAuxiliaryVariables(int num_pml_elements) {
    int n_aux = getNumStateVariables();
    aux_vars.resize(num_pml_elements);
    for (auto& elem_aux : aux_vars) {
        elem_aux.resize(n_aux, 0.0);
    }
}

void PMLBoundary::computeAuxiliaryRHS(int elem_id, const double* u,
                                      const double* aux, double* aux_rhs) {
    // Split-field PML auxiliary variable evolution
    // For each field variable q_i:
    // ∂ψ_x/∂t = -d_x * ψ_x + ∂q/∂x * (s_x - 1)
    
    // Simplified: just apply damping to auxiliary variables
    int n_aux = getNumStateVariables();
    for (int i = 0; i < n_aux; ++i) {
        // Simple exponential decay
        aux_rhs[i] = -d_max * aux[i];
    }
}

void PMLBoundary::updateAuxiliaryVariables(int elem_id, const double* aux_new) {
    if (elem_id < 0 || elem_id >= static_cast<int>(aux_vars.size())) return;
    
    int n_aux = getNumStateVariables();
    for (int i = 0; i < n_aux; ++i) {
        aux_vars[elem_id][i] = aux_new[i];
    }
}

void PMLBoundary::getPMLJacobian(int elem_id, double* jacobian) const {
    // Modified Jacobian for PML
    // The flux Jacobian is stretched by the PML functions
    // A_pml = s_y * s_z / s_x * A_x + s_x * s_z / s_y * A_y + s_x * s_y / s_z * A_z
    
    // Return identity for now (full implementation would compute stretched Jacobians)
    int n = 9;
    for (int i = 0; i < n * n; ++i) {
        jacobian[i] = 0.0;
    }
    for (int i = 0; i < n; ++i) {
        jacobian[i * n + i] = 1.0;
    }
}

// =============================================================================
// ClaytonEngquistBC Implementation
// =============================================================================

ClaytonEngquistBC::ClaytonEngquistBC()
    : order(1), rho(2500.0), vp(5000.0), vs(3000.0) {
    lambda = rho * (vp * vp - 2.0 * vs * vs);
    mu = rho * vs * vs;
}

void ClaytonEngquistBC::setOrder(int ord) {
    order = std::max(1, std::min(2, ord));
}

void ClaytonEngquistBC::setMaterialProperties(double r, double p, double s) {
    rho = r;
    vp = p;
    vs = s;
    lambda = rho * (vp * vp - 2.0 * vs * vs);
    mu = rho * vs * vs;
}

void ClaytonEngquistBC::applyToResidual(const BoundaryFace& face,
                                        const double* u, double* residual,
                                        double t) const {
    // First-order Clayton-Engquist: t = -Z · v
    double traction[3];
    applyFirstOrder(face, u, traction);
    
    double local_rho = face.rho > 0 ? face.rho : rho;
    
    // Add traction contribution to momentum residual
    residual[0] -= traction[0] / local_rho;
    residual[1] -= traction[1] / local_rho;
    residual[2] -= traction[2] / local_rho;
}

void ClaytonEngquistBC::computeBoundaryFlux(const BoundaryFace& face,
                                            const double* u_int, double* flux_n,
                                            double t) const {
    double u_ext[9];
    getBoundaryState(face, u_int, u_ext, t);
    
    const double* n = face.normal.data();
    double local_rho = face.rho > 0 ? face.rho : rho;
    double local_vp = face.vp > 0 ? face.vp : vp;
    double local_vs = face.vs > 0 ? face.vs : vs;
    
    // Upwind flux based on characteristic decomposition
    // P-wave characteristic: λ = vp
    // S-wave characteristic: λ = vs
    
    double vx_int = u_int[0], vy_int = u_int[1], vz_int = u_int[2];
    double vx_ext = u_ext[0], vy_ext = u_ext[1], vz_ext = u_ext[2];
    
    // Normal velocity
    double vn_int = vx_int * n[0] + vy_int * n[1] + vz_int * n[2];
    double vn_ext = vx_ext * n[0] + vy_ext * n[1] + vz_ext * n[2];
    
    // Use upwind (interior) state for outgoing waves
    double vn = (vn_int > 0) ? vn_int : vn_ext;
    double vx = (vn_int > 0) ? vx_int : vx_ext;
    double vy = (vn_int > 0) ? vy_int : vy_ext;
    double vz = (vn_int > 0) ? vz_int : vz_ext;
    
    // Compute flux
    flux_n[0] = vx;
    flux_n[1] = vy;
    flux_n[2] = vz;
    
    // Stress flux (simplified)
    for (int i = 3; i < 9; ++i) {
        flux_n[i] = 0.5 * (u_int[i] + u_ext[i]);
    }
}

void ClaytonEngquistBC::getBoundaryState(const BoundaryFace& face,
                                         const double* u_int, double* u_ext,
                                         double t) const {
    // Characteristic-based non-reflecting BC
    const double* n = face.normal.data();
    
    double local_rho = face.rho > 0 ? face.rho : rho;
    double local_vp = face.vp > 0 ? face.vp : vp;
    double local_vs = face.vs > 0 ? face.vs : vs;
    
    double Zp = local_rho * local_vp;
    double Zs = local_rho * local_vs;
    
    // Interior velocity
    double vx = u_int[0], vy = u_int[1], vz = u_int[2];
    double vn = vx * n[0] + vy * n[1] + vz * n[2];
    
    // Interior stress (if available)
    double sxx = u_int[3], syy = u_int[4], szz = u_int[5];
    double sxy = u_int[6], sxz = u_int[7], syz = u_int[8];
    
    // Normal traction
    double tn_x = sxx * n[0] + sxy * n[1] + sxz * n[2];
    double tn_y = sxy * n[0] + syy * n[1] + syz * n[2];
    double tn_z = sxz * n[0] + syz * n[1] + szz * n[2];
    double tn = tn_x * n[0] + tn_y * n[1] + tn_z * n[2];
    
    // Characteristic variables for P-wave: W_p = t_n + Z_p * v_n
    // For non-reflecting: W_p^- = 0 -> t_n = -Z_p * v_n
    
    // Exterior state that gives zero incoming characteristic
    double tn_ext = -Zp * vn;
    
    // Set exterior velocity (same as interior for outgoing)
    u_ext[0] = vx;
    u_ext[1] = vy;
    u_ext[2] = vz;
    
    // Compute exterior stress from exterior traction
    // This is a simplified version
    u_ext[3] = tn_ext * n[0] * n[0];
    u_ext[4] = tn_ext * n[1] * n[1];
    u_ext[5] = tn_ext * n[2] * n[2];
    u_ext[6] = tn_ext * n[0] * n[1];
    u_ext[7] = tn_ext * n[0] * n[2];
    u_ext[8] = tn_ext * n[1] * n[2];
}

void ClaytonEngquistBC::computeImpedanceMatrix(const double* n, double* Z) const {
    // Impedance matrix Z such that t = Z · v
    // For isotropic material:
    // Z = ρ * vp * (n ⊗ n) + ρ * vs * (I - n ⊗ n)
    
    double Zp = rho * vp;
    double Zs = rho * vs;
    
    // Z_ij = Zs * δ_ij + (Zp - Zs) * n_i * n_j
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            Z[i * 3 + j] = (i == j ? Zs : 0.0) + (Zp - Zs) * n[i] * n[j];
        }
    }
}

void ClaytonEngquistBC::applyFirstOrder(const BoundaryFace& face,
                                        const double* v, double* traction) const {
    // First-order paraxial: t = -Z · v
    double Z[9];
    computeImpedanceMatrix(face.normal.data(), Z);
    
    traction[0] = -(Z[0] * v[0] + Z[1] * v[1] + Z[2] * v[2]);
    traction[1] = -(Z[3] * v[0] + Z[4] * v[1] + Z[5] * v[2]);
    traction[2] = -(Z[6] * v[0] + Z[7] * v[1] + Z[8] * v[2]);
}

void ClaytonEngquistBC::applySecondOrder(const BoundaryFace& face,
                                         const double* v, const double* dv_dt,
                                         double* traction) const {
    // Second-order includes time derivatives
    // Not fully implemented - use first order
    applyFirstOrder(face, v, traction);
}

// =============================================================================
// LysmerKuhlemeyer Implementation
// =============================================================================

LysmerKuhlemeyer::LysmerKuhlemeyer()
    : rho(2500.0), vp(5000.0), vs(3000.0) {
    impedance_p = rho * vp;
    impedance_s = rho * vs;
}

void LysmerKuhlemeyer::setMaterialProperties(double r, double p, double s) {
    rho = r;
    vp = p;
    vs = s;
    impedance_p = rho * vp;
    impedance_s = rho * vs;
}

void LysmerKuhlemeyer::applyToResidual(const BoundaryFace& face,
                                       const double* u, double* residual,
                                       double t) const {
    double traction[3];
    computeAbsorbingTraction(u, face.normal.data(), traction);
    
    double local_rho = face.rho > 0 ? face.rho : rho;
    
    residual[0] -= traction[0] / local_rho;
    residual[1] -= traction[1] / local_rho;
    residual[2] -= traction[2] / local_rho;
}

void LysmerKuhlemeyer::computeBoundaryFlux(const BoundaryFace& face,
                                           const double* u_int, double* flux_n,
                                           double t) const {
    double u_ext[9];
    getBoundaryState(face, u_int, u_ext, t);
    
    for (int i = 0; i < 9; ++i) {
        flux_n[i] = 0.5 * (u_int[i] + u_ext[i]);
    }
}

void LysmerKuhlemeyer::getBoundaryState(const BoundaryFace& face,
                                        const double* u_int, double* u_ext,
                                        double t) const {
    // Set exterior state for absorbing boundary
    const double* n = face.normal.data();
    
    // Exterior velocity: zero (incoming wave = 0)
    u_ext[0] = 0;
    u_ext[1] = 0;
    u_ext[2] = 0;
    
    // Exterior stress: absorbing condition
    for (int i = 3; i < 9; ++i) {
        u_ext[i] = -u_int[i];  // Anti-symmetric for absorption
    }
}

void LysmerKuhlemeyer::computeAbsorbingTraction(const double* v, const double* n,
                                                double* traction) const {
    // t = -ρ * vp * vn * n - ρ * vs * vt
    double vn = v[0] * n[0] + v[1] * n[1] + v[2] * n[2];
    
    // Normal traction
    double tn = -impedance_p * vn;
    
    // Tangential velocity
    double vtx = v[0] - vn * n[0];
    double vty = v[1] - vn * n[1];
    double vtz = v[2] - vn * n[2];
    
    // Tangential traction
    double ttx = -impedance_s * vtx;
    double tty = -impedance_s * vty;
    double ttz = -impedance_s * vtz;
    
    traction[0] = tn * n[0] + ttx;
    traction[1] = tn * n[1] + tty;
    traction[2] = tn * n[2] + ttz;
}

// =============================================================================
// DirichletBC Implementation
// =============================================================================

DirichletBC::DirichletBC() {
    setZero();
}

void DirichletBC::setPrescribedValue(
    std::function<void(double, double, double, double, double*)> func) {
    prescribed_func = func;
}

void DirichletBC::setZero() {
    prescribed_func = [](double, double, double, double, double* u) {
        for (int i = 0; i < 9; ++i) u[i] = 0.0;
    };
}

void DirichletBC::applyToResidual(const BoundaryFace& face,
                                  const double* u, double* residual,
                                  double t) const {
    // Strong enforcement: u = u_prescribed
    double u_prescribed[9];
    prescribed_func(face.centroid[0], face.centroid[1], face.centroid[2], 
                   t, u_prescribed);
    
    double penalty = 1e10;
    for (int i = 0; i < 9; ++i) {
        residual[i] += penalty * (u_prescribed[i] - u[i]);
    }
}

void DirichletBC::computeBoundaryFlux(const BoundaryFace& face,
                                      const double* u_int, double* flux_n,
                                      double t) const {
    double u_ext[9];
    getBoundaryState(face, u_int, u_ext, t);
    
    for (int i = 0; i < 9; ++i) {
        flux_n[i] = 0.5 * (u_int[i] + u_ext[i]);
    }
}

void DirichletBC::getBoundaryState(const BoundaryFace& face,
                                   const double* u_int, double* u_ext,
                                   double t) const {
    // Ghost state: 2 * u_bc - u_int (for strong enforcement)
    double u_bc[9];
    prescribed_func(face.centroid[0], face.centroid[1], face.centroid[2],
                   t, u_bc);
    
    for (int i = 0; i < 9; ++i) {
        u_ext[i] = 2.0 * u_bc[i] - u_int[i];
    }
}

// =============================================================================
// NeumannBC Implementation
// =============================================================================

NeumannBC::NeumannBC()
    : is_constant(true), constant_traction({0, 0, 0}) {}

void NeumannBC::setPrescribedTraction(
    std::function<void(double, double, double, double, double*)> func) {
    traction_func = func;
    is_constant = false;
}

void NeumannBC::setConstantTraction(double tx, double ty, double tz) {
    constant_traction = {tx, ty, tz};
    is_constant = true;
    traction_func = [tx, ty, tz](double, double, double, double, double* t) {
        t[0] = tx; t[1] = ty; t[2] = tz;
    };
}

void NeumannBC::applyToResidual(const BoundaryFace& face,
                                const double* u, double* residual,
                                double t) const {
    double traction[3];
    if (is_constant) {
        traction[0] = constant_traction[0];
        traction[1] = constant_traction[1];
        traction[2] = constant_traction[2];
    } else {
        traction_func(face.centroid[0], face.centroid[1], face.centroid[2],
                     t, traction);
    }
    
    double local_rho = face.rho > 0 ? face.rho : 2500.0;
    
    // Add traction to momentum equation
    residual[0] += traction[0] / local_rho;
    residual[1] += traction[1] / local_rho;
    residual[2] += traction[2] / local_rho;
}

void NeumannBC::computeBoundaryFlux(const BoundaryFace& face,
                                    const double* u_int, double* flux_n,
                                    double t) const {
    // For Neumann BC, flux is specified directly
    double traction[3];
    if (is_constant) {
        traction[0] = constant_traction[0];
        traction[1] = constant_traction[1];
        traction[2] = constant_traction[2];
    } else {
        traction_func(face.centroid[0], face.centroid[1], face.centroid[2],
                     t, traction);
    }
    
    flux_n[0] = traction[0];
    flux_n[1] = traction[1];
    flux_n[2] = traction[2];
    
    for (int i = 3; i < 9; ++i) {
        flux_n[i] = u_int[i];  // Use interior stress
    }
}

void NeumannBC::getBoundaryState(const BoundaryFace& face,
                                 const double* u_int, double* u_ext,
                                 double t) const {
    // Copy interior state
    for (int i = 0; i < 9; ++i) {
        u_ext[i] = u_int[i];
    }
}

// =============================================================================
// PeriodicBC Implementation
// =============================================================================

PeriodicBC::PeriodicBC()
    : periodic_dir(0), period(1.0) {}

void PeriodicBC::setPeriodicMapping(
    const std::vector<std::pair<int, int>>& face_pairs) {
    periodic_pairs.clear();
    for (const auto& pair : face_pairs) {
        periodic_pairs[pair.first] = pair.second;
        periodic_pairs[pair.second] = pair.first;
    }
}

void PeriodicBC::setPeriodicDirection(int dir, double p) {
    periodic_dir = dir;
    period = p;
}

void PeriodicBC::applyToResidual(const BoundaryFace& face,
                                 const double* u, double* residual,
                                 double t) const {
    // Periodic BC enforces u(x) = u(x + L)
    // Handled through face matching, not residual modification
}

void PeriodicBC::computeBoundaryFlux(const BoundaryFace& face,
                                     const double* u_int, double* flux_n,
                                     double t) const {
    // Standard interior flux (periodic face is treated as internal)
    for (int i = 0; i < 9; ++i) {
        flux_n[i] = u_int[i];
    }
}

void PeriodicBC::getBoundaryState(const BoundaryFace& face,
                                  const double* u_int, double* u_ext,
                                  double t) const {
    // State from matched face
    // In practice, this would look up the matched face's solution
    for (int i = 0; i < 9; ++i) {
        u_ext[i] = u_int[i];  // Placeholder - actual impl gets from matched face
    }
}

// =============================================================================
// SymmetryBC Implementation
// =============================================================================

SymmetryBC::SymmetryBC(bool is_anti)
    : antisymmetric(is_anti) {}

void SymmetryBC::applyToResidual(const BoundaryFace& face,
                                 const double* u, double* residual,
                                 double t) const {
    // Symmetry: v_n = 0 (normal velocity zero)
    // Anti-symmetry: v_t = 0 (tangential velocity zero)
    
    const double* n = face.normal.data();
    double vn = u[0] * n[0] + u[1] * n[1] + u[2] * n[2];
    
    double penalty = 1e8;
    
    if (!antisymmetric) {
        // Enforce v_n = 0
        residual[0] += penalty * vn * n[0];
        residual[1] += penalty * vn * n[1];
        residual[2] += penalty * vn * n[2];
    } else {
        // Enforce v_t = 0
        double vtx = u[0] - vn * n[0];
        double vty = u[1] - vn * n[1];
        double vtz = u[2] - vn * n[2];
        
        residual[0] += penalty * vtx;
        residual[1] += penalty * vty;
        residual[2] += penalty * vtz;
    }
}

void SymmetryBC::computeBoundaryFlux(const BoundaryFace& face,
                                     const double* u_int, double* flux_n,
                                     double t) const {
    double u_ext[9];
    getBoundaryState(face, u_int, u_ext, t);
    
    for (int i = 0; i < 9; ++i) {
        flux_n[i] = 0.5 * (u_int[i] + u_ext[i]);
    }
}

void SymmetryBC::getBoundaryState(const BoundaryFace& face,
                                  const double* u_int, double* u_ext,
                                  double t) const {
    const double* n = face.normal.data();
    
    // Reflect velocity
    reflectVelocity(u_int, n, u_ext);
    
    // Reflect stress
    reflectStress(&u_int[3], n, &u_ext[3]);
}

void SymmetryBC::reflectVelocity(const double* v_int, const double* n,
                                 double* v_ext) const {
    double vn = v_int[0] * n[0] + v_int[1] * n[1] + v_int[2] * n[2];
    
    if (!antisymmetric) {
        // Symmetry: v_ext = v_int - 2*vn*n (reflect normal component)
        v_ext[0] = v_int[0] - 2.0 * vn * n[0];
        v_ext[1] = v_int[1] - 2.0 * vn * n[1];
        v_ext[2] = v_int[2] - 2.0 * vn * n[2];
    } else {
        // Anti-symmetry: v_ext = -v_int + 2*vn*n (reflect tangential component)
        v_ext[0] = -v_int[0] + 2.0 * vn * n[0];
        v_ext[1] = -v_int[1] + 2.0 * vn * n[1];
        v_ext[2] = -v_int[2] + 2.0 * vn * n[2];
    }
}

void SymmetryBC::reflectStress(const double* sigma_int, const double* n,
                               double* sigma_ext) const {
    // Stress reflection similar to velocity
    // For symmetry plane, normal-tangent shear components flip sign
    
    double sxx = sigma_int[0], syy = sigma_int[1], szz = sigma_int[2];
    double sxy = sigma_int[3], sxz = sigma_int[4], syz = sigma_int[5];
    
    // Simplified: just copy for now
    sigma_ext[0] = sxx;
    sigma_ext[1] = syy;
    sigma_ext[2] = szz;
    sigma_ext[3] = antisymmetric ? -sxy : sxy;
    sigma_ext[4] = antisymmetric ? -sxz : sxz;
    sigma_ext[5] = antisymmetric ? -syz : syz;
}

// =============================================================================
// BoundaryConditionManager Implementation
// =============================================================================

BoundaryConditionManager::BoundaryConditionManager()
    : dg(nullptr) {}

BoundaryConditionManager::~BoundaryConditionManager() = default;

void BoundaryConditionManager::initialize(DiscontinuousGalerkin* dg_ptr) {
    dg = dg_ptr;
}

void BoundaryConditionManager::addBoundaryCondition(
    int tag, std::unique_ptr<BoundaryConditionBase> bc) {
    bc_by_tag[tag] = std::move(bc);
}

void BoundaryConditionManager::setDefaultBC(
    std::unique_ptr<BoundaryConditionBase> bc) {
    default_bc = std::move(bc);
}

void BoundaryConditionManager::classifyBoundaryFaces(
    const std::vector<int>& face_tags) {
    // This would be called after mesh is loaded
    // Classify each boundary face based on its tag
    
    all_boundary_faces.clear();
    faces_by_type.clear();
    face_bc.resize(face_tags.size());
    
    for (size_t i = 0; i < face_tags.size(); ++i) {
        BoundaryFace face;
        face.element_id = i;  // Simplified - actual impl gets from mesh
        
        int tag = face_tags[i];
        
        auto it = bc_by_tag.find(tag);
        if (it != bc_by_tag.end()) {
            face_bc[i] = it->second.get();
        } else {
            face_bc[i] = default_bc.get();
        }
        
        all_boundary_faces.push_back(face);
    }
}

void BoundaryConditionManager::setupPML(double pml_thickness, int directions) {
    pml = std::make_unique<PMLBoundary>();
    pml->setThickness(pml_thickness);
    
    // Would need domain bounds from mesh
    // pml->setDomainBounds(...)
}

void BoundaryConditionManager::applyBoundaryConditions(
    double* residual, const double* solution, double t) {
    
    for (size_t i = 0; i < all_boundary_faces.size(); ++i) {
        if (face_bc[i]) {
            const BoundaryFace& face = all_boundary_faces[i];
            // Get solution at face
            const double* u_face = &solution[face.element_id * 9];  // Simplified
            double* res_face = &residual[face.element_id * 9];
            
            face_bc[i]->applyToResidual(face, u_face, res_face, t);
        }
    }
}

void BoundaryConditionManager::computeBoundaryFluxes(
    const double* solution, double* flux, double t) {
    
    for (size_t i = 0; i < all_boundary_faces.size(); ++i) {
        if (face_bc[i]) {
            const BoundaryFace& face = all_boundary_faces[i];
            const double* u_face = &solution[face.element_id * 9];
            double* flux_face = &flux[i * 9];
            
            face_bc[i]->computeBoundaryFlux(face, u_face, flux_face, t);
        }
    }
}

const std::vector<BoundaryFace>& BoundaryConditionManager::getBoundaryFaces(
    BoundaryType type) const {
    static std::vector<BoundaryFace> empty;
    auto it = faces_by_type.find(type);
    return (it != faces_by_type.end()) ? it->second : empty;
}

const std::vector<BoundaryFace>& BoundaryConditionManager::getAllBoundaryFaces() const {
    return all_boundary_faces;
}

BoundaryConditionBase* BoundaryConditionManager::getBCForFace(int face_id) const {
    if (face_id < 0 || face_id >= static_cast<int>(face_bc.size())) {
        return nullptr;
    }
    return face_bc[face_id];
}

void BoundaryConditionManager::updateMaterialProperties() {
    // Update material properties on boundary faces from mesh/material data
}

int BoundaryConditionManager::getNumPMLAuxVariables() const {
    if (!pml) return 0;
    return pml->getNumStateVariables() * pml_elements.size();
}

void BoundaryConditionManager::evolvePMLAuxVariables(
    const double* solution, double dt) {
    if (!pml) return;
    
    // Evolve PML auxiliary variables using simple Euler
    for (size_t i = 0; i < pml_elements.size(); ++i) {
        int elem = pml_elements[i];
        const double* u = &solution[elem * 9];  // Simplified indexing
        
        std::vector<double> aux(pml->getNumStateVariables());
        std::vector<double> aux_rhs(pml->getNumStateVariables());
        
        pml->computeAuxiliaryRHS(i, u, aux.data(), aux_rhs.data());
        
        for (int j = 0; j < pml->getNumStateVariables(); ++j) {
            aux[j] += dt * aux_rhs[j];
        }
        
        pml->updateAuxiliaryVariables(i, aux.data());
    }
}

// =============================================================================
// BoundaryConditionFactory Implementation
// =============================================================================

std::unique_ptr<BoundaryConditionBase> BoundaryConditionFactory::create(
    const std::string& type) {
    
    if (type == "free_surface") {
        return std::make_unique<FreeSurfaceBC>();
    } else if (type == "absorbing_pml" || type == "pml") {
        return std::make_unique<PMLBoundary>();
    } else if (type == "absorbing_ce" || type == "clayton_engquist") {
        return std::make_unique<ClaytonEngquistBC>();
    } else if (type == "absorbing_lysmer" || type == "lysmer") {
        return std::make_unique<LysmerKuhlemeyer>();
    } else if (type == "dirichlet" || type == "fixed") {
        return std::make_unique<DirichletBC>();
    } else if (type == "neumann" || type == "traction") {
        return std::make_unique<NeumannBC>();
    } else if (type == "periodic") {
        return std::make_unique<PeriodicBC>();
    } else if (type == "symmetry") {
        return std::make_unique<SymmetryBC>(false);
    } else if (type == "antisymmetry") {
        return std::make_unique<SymmetryBC>(true);
    }
    
    // Default to free surface
    return std::make_unique<FreeSurfaceBC>();
}

std::unique_ptr<BoundaryConditionBase> BoundaryConditionFactory::createFromConfig(
    const std::map<std::string, std::string>& config,
    const std::string& prefix) {
    
    auto get_string = [&](const std::string& key, const std::string& def) {
        auto it = config.find(prefix + key);
        return (it != config.end()) ? it->second : def;
    };
    
    auto get_double = [&](const std::string& key, double def) {
        auto it = config.find(prefix + key);
        if (it != config.end()) {
            try { return std::stod(it->second); }
            catch (...) { return def; }
        }
        return def;
    };
    
    std::string type = get_string("type", "free_surface");
    auto bc = create(type);
    
    // Configure material properties
    double rho = get_double("density", 2500.0);
    double vp = get_double("vp", 5000.0);
    double vs = get_double("vs", 3000.0);
    
    // Type-specific configuration
    if (type == "absorbing_ce" || type == "clayton_engquist") {
        auto* ce = dynamic_cast<ClaytonEngquistBC*>(bc.get());
        if (ce) {
            ce->setOrder(static_cast<int>(get_double("order", 1.0)));
            ce->setMaterialProperties(rho, vp, vs);
        }
    } else if (type == "absorbing_lysmer" || type == "lysmer") {
        auto* lk = dynamic_cast<LysmerKuhlemeyer*>(bc.get());
        if (lk) {
            lk->setMaterialProperties(rho, vp, vs);
        }
    }
    
    return bc;
}

std::unique_ptr<FreeSurfaceBC> BoundaryConditionFactory::createFreeSurface() {
    return std::make_unique<FreeSurfaceBC>();
}

std::unique_ptr<BoundaryConditionBase> BoundaryConditionFactory::createAbsorbing(
    const std::string& method,
    double rho, double vp, double vs) {
    
    if (method == "pml") {
        auto bc = std::make_unique<PMLBoundary>();
        bc->setMaterialProperties(rho, vp, vs);
        return bc;
    } else if (method == "clayton_engquist" || method == "ce") {
        auto bc = std::make_unique<ClaytonEngquistBC>();
        bc->setMaterialProperties(rho, vp, vs);
        return bc;
    } else {
        auto bc = std::make_unique<LysmerKuhlemeyer>();
        bc->setMaterialProperties(rho, vp, vs);
        return bc;
    }
}

std::unique_ptr<PMLBoundary> BoundaryConditionFactory::createPML(
    double thickness, double reflection_coef,
    double x_min, double x_max,
    double y_min, double y_max,
    double z_min, double z_max) {
    
    auto pml = std::make_unique<PMLBoundary>();
    pml->setThickness(thickness);
    pml->setReflectionCoefficient(reflection_coef);
    pml->setDomainBounds(x_min, x_max, y_min, y_max, z_min, z_max);
    return pml;
}

// =============================================================================
// BoundaryConditionConfig Implementation
// =============================================================================

void BoundaryConditionConfig::parseConfig(
    const std::map<std::string, std::string>& config) {
    
    auto get_string = [&](const std::string& key, const std::string& def) {
        auto it = config.find(key);
        return (it != config.end()) ? it->second : def;
    };
    
    auto get_double = [&](const std::string& key, double def) {
        auto it = config.find(key);
        if (it != config.end()) {
            try { return std::stod(it->second); }
            catch (...) { return def; }
        }
        return def;
    };
    
    auto get_int = [&](const std::string& key, int def) {
        auto it = config.find(key);
        if (it != config.end()) {
            try { return std::stoi(it->second); }
            catch (...) { return def; }
        }
        return def;
    };
    
    default_bc_type = get_string("default_boundary_condition", "free_surface");
    absorbing_method = get_string("absorbing_method", "pml");
    absorbing_order = get_int("absorbing_order", 1);
    
    pml_thickness = get_double("pml_thickness", 10.0);
    pml_reflection_coef = get_double("pml_reflection_coefficient", 1e-6);
    pml_damping_order = get_int("pml_damping_order", 3);
    pml_alpha_max = get_double("pml_alpha_max", M_PI);
    pml_kappa_max = get_double("pml_kappa_max", 1.0);
    
    boundary_density = get_double("boundary_density", 2500.0);
    boundary_vp = get_double("boundary_vp", 5000.0);
    boundary_vs = get_double("boundary_vs", 3000.0);
    
    use_characteristic_bc = get_string("use_characteristic_bc", "true") == "true";
    
    // Parse BC types for specific tags (format: "bc_type_1=free_surface")
    for (const auto& kv : config) {
        if (kv.first.substr(0, 8) == "bc_type_") {
            try {
                int tag = std::stoi(kv.first.substr(8));
                bc_types[tag] = kv.second;
            } catch (...) {}
        }
    }
}

void BoundaryConditionConfig::setupBoundaryConditions(
    BoundaryConditionManager& bcm) const {
    
    // Set default BC
    bcm.setDefaultBC(BoundaryConditionFactory::create(default_bc_type));
    
    // Add BCs for specific tags
    for (const auto& kv : bc_types) {
        bcm.addBoundaryCondition(kv.first, 
            BoundaryConditionFactory::create(kv.second));
    }
    
    // Setup PML if requested
    if (absorbing_method == "pml") {
        bcm.setupPML(pml_thickness);
        if (bcm.getPML()) {
            bcm.getPML()->setReflectionCoefficient(pml_reflection_coef);
            bcm.getPML()->setOrder(pml_damping_order);
            bcm.getPML()->setCFSParameters(pml_alpha_max, pml_kappa_max);
            bcm.getPML()->setMaterialProperties(boundary_density, 
                                               boundary_vp, boundary_vs);
        }
    }
}

} // namespace FSRM
