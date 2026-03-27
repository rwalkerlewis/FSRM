/**
 * @file MultiphaseFlowKernel.cpp
 * @brief PetscDS-compatible multiphase flow residual and Jacobian callbacks
 *
 * Implements the coupled pressure-saturation weak form for two-phase
 * (water/brine + CO2) flow in porous media. All static callback functions
 * follow the exact PetscPointFunc / PetscPointJac signatures.
 *
 * Pressure equation (field 0):
 *   f0: φ * ct_eff * dP/dt   (accumulation with effective compressibility)
 *   f1: Σ_α k * kr_α/μ_α * (∇P - ρ_α*g*ê_z)  (total Darcy flux)
 *
 * Saturation equation (field 1):
 *   f0: φ * dSw/dt   (water phase accumulation)
 *   f1: k * kr_w/μ_w * (∇P - ρ_w*g*ê_z)  (water phase Darcy flux)
 */

#include "physics/MultiphaseFlowKernel.hpp"
#include <cmath>
#include <algorithm>

namespace FSRM {

// ============================================================================
// Constructor and setup
// ============================================================================

MultiphaseFlowKernel::MultiphaseFlowKernel()
    : PhysicsKernel(PhysicsType::FLUID_FLOW),
      porosity_(0.15), perm_x_(200e-15), perm_y_(200e-15), perm_z_(100e-15),
      water_density_(1100.0), water_viscosity_(8e-4), water_compressibility_(4e-10),
      co2_density_(700.0), co2_viscosity_(5e-5), co2_compressibility_(1e-8),
      Swr_(0.2), Sgr_(0.05), nw_(4.0), ng_(2.0), krw_max_(1.0), krg_max_(1.0) {}

PetscErrorCode MultiphaseFlowKernel::setup(DM dm, PetscFE fe) {
    (void)dm;
    (void)fe;
    PetscFunctionBeginUser;
    PetscFunctionReturn(0);
}

void MultiphaseFlowKernel::residual(const PetscScalar u[], const PetscScalar u_t[],
                                    const PetscScalar u_x[], const PetscScalar a[],
                                    const PetscReal x[], PetscScalar f[]) {
    (void)a; (void)x;
    // Simplified pointwise evaluation for unit testing
    // u[0] = P, u[1] = Sw
    double Sw = PetscRealPart(u[1]);
    double ct = Sw * water_compressibility_ + (1.0 - Sw) * co2_compressibility_;

    // Pressure residual: accumulation
    f[0] = porosity_ * ct * u_t[0];

    // Saturation residual: accumulation
    f[1] = porosity_ * u_t[1];

    // Pointwise test path: flux uses PetscDS callbacks with a full constants[]
    // array; avoid computeKrw/Sw with nullptr here.
    (void)u_x;
}

void MultiphaseFlowKernel::jacobian(const PetscScalar u[], const PetscScalar u_t[],
                                    const PetscScalar u_x[], const PetscScalar a[],
                                    const PetscReal x[], PetscScalar J[]) {
    (void)u_t; (void)u_x; (void)x;
    double shift = PetscRealPart(a[0]);
    double Sw = PetscRealPart(u[1]);
    double ct = Sw * water_compressibility_ + (1.0 - Sw) * co2_compressibility_;

    // 2x2 Jacobian
    J[0] = porosity_ * ct * shift;         // ∂f0/∂P
    J[1] = porosity_ * (water_compressibility_ - co2_compressibility_)
           * PetscRealPart(u_t[0]);         // ∂f0/∂Sw
    J[2] = 0.0;                             // ∂f1/∂P
    J[3] = porosity_ * shift;              // ∂f1/∂Sw
}

// ============================================================================
// Property setters
// ============================================================================

void MultiphaseFlowKernel::setRockProperties(double phi, double kx, double ky, double kz) {
    porosity_ = phi; perm_x_ = kx; perm_y_ = ky; perm_z_ = kz;
}

void MultiphaseFlowKernel::setWaterProperties(double rho, double mu, double cw) {
    water_density_ = rho; water_viscosity_ = mu; water_compressibility_ = cw;
}

void MultiphaseFlowKernel::setCO2Properties(double rho, double mu, double cg) {
    co2_density_ = rho; co2_viscosity_ = mu; co2_compressibility_ = cg;
}

void MultiphaseFlowKernel::setRelPermParams(double Swr, double Sgr, double nw,
                                             double ng, double krw_max, double krg_max) {
    Swr_ = Swr; Sgr_ = Sgr; nw_ = nw; ng_ = ng;
    krw_max_ = krw_max; krg_max_ = krg_max;
}

// ============================================================================
// PetscDS Residual: Pressure equation
// ============================================================================

void MultiphaseFlowKernel::f0_pressure(
    PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[],
    const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, const PetscReal x[],
    PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
    (void)dim; (void)Nf; (void)NfAux;
    (void)uOff_x; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)x;

    const PetscScalar Sw  = u[uOff[1]];
    const PetscScalar P_t = u_t[uOff[0]];
    const PetscScalar S_t = u_t[uOff[1]];

    const PetscScalar phi = (numConstants > 0) ? constants[0] : 0.15;
    const PetscScalar cw  = (numConstants > 8) ? constants[8] : 4e-10;
    const PetscScalar cg  = (numConstants > 9) ? constants[9] : 1e-8;
    const PetscScalar rho_w = (numConstants > 6) ? constants[6] : 1100.0;
    const PetscScalar rho_g = (numConstants > 7) ? constants[7] : 700.0;

    // Effective compressibility: ct = Sw*cw + Sg*cg
    PetscScalar Sg = 1.0 - Sw;
    PetscScalar ct_eff = Sw * cw + Sg * cg;

    // f0 = φ*ct*dP/dt + φ*(ρw - ρg)*dSw/dt  (density-difference coupling)
    f0[0] = phi * ct_eff * P_t + phi * (rho_w - rho_g) * S_t * 1e-9;
    // Note: the density-difference term is scaled by 1e-9 for dimensional consistency
    // In a full implementation, this would use specific volume differences
}

void MultiphaseFlowKernel::f1_pressure(
    PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[],
    const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, const PetscReal x[],
    PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{
    (void)Nf; (void)NfAux;
    (void)uOff; (void)u; (void)u_t;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)x;

    const PetscScalar Sw = u[uOff[1]];

    // Permeabilities
    const PetscScalar kx = (numConstants > 1) ? constants[1] : 200e-15;
    const PetscScalar ky = (numConstants > 2) ? constants[2] : 200e-15;
    const PetscScalar kz = (numConstants > 3) ? constants[3] : 100e-15;

    // Viscosities
    const PetscScalar mu_w = (numConstants > 4) ? constants[4] : 8e-4;
    const PetscScalar mu_g = (numConstants > 5) ? constants[5] : 5e-5;

    // Densities
    const PetscScalar rho_w = (numConstants > 6) ? constants[6] : 1100.0;
    const PetscScalar rho_g = (numConstants > 7) ? constants[7] : 700.0;

    // Gravity
    const PetscScalar grav = (numConstants > 16) ? constants[16] : 9.81;

    // Relative permeabilities
    double krw = computeKrw(Sw, constants);
    double krg = computeKrg(Sw, constants);

    // Total mobility per direction
    double mob_w_x = kx * krw / mu_w;
    double mob_g_x = kx * krg / mu_g;
    double mob_w_y = ky * krw / mu_w;
    double mob_g_y = ky * krg / mu_g;
    double mob_w_z = kz * krw / mu_w;
    double mob_g_z = kz * krg / mu_g;

    // Pressure gradient components
    for (PetscInt d = 0; d < dim; ++d) {
        PetscScalar dP_dd = u_x[uOff_x[0] + d];

        double mob_w_d, mob_g_d;
        if (d == 0) { mob_w_d = mob_w_x; mob_g_d = mob_g_x; }
        else if (d == 1) { mob_w_d = mob_w_y; mob_g_d = mob_g_y; }
        else { mob_w_d = mob_w_z; mob_g_d = mob_g_z; }

        // Gravity acts in z-direction (last dimension)
        PetscScalar grav_w = (d == dim - 1) ? rho_w * grav : 0.0;
        PetscScalar grav_g = (d == dim - 1) ? rho_g * grav : 0.0;

        // Total flux = water flux + CO2 flux
        f1[d] = mob_w_d * (dP_dd - grav_w) + mob_g_d * (dP_dd - grav_g);
    }
}

// ============================================================================
// PetscDS Residual: Saturation equation
// ============================================================================

void MultiphaseFlowKernel::f0_saturation(
    PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[],
    const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, const PetscReal x[],
    PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
    (void)dim; (void)Nf; (void)NfAux;
    (void)uOff_x; (void)u; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)x;

    const PetscScalar S_t = u_t[uOff[1]];
    const PetscScalar phi = (numConstants > 0) ? constants[0] : 0.15;

    f0[0] = phi * S_t;
}

void MultiphaseFlowKernel::f1_saturation(
    PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[],
    const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, const PetscReal x[],
    PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{
    (void)Nf; (void)NfAux;
    (void)u; (void)u_t;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)x;

    const PetscScalar Sw = u[uOff[1]];

    const PetscScalar kx = (numConstants > 1) ? constants[1] : 200e-15;
    const PetscScalar ky = (numConstants > 2) ? constants[2] : 200e-15;
    const PetscScalar kz = (numConstants > 3) ? constants[3] : 100e-15;
    const PetscScalar mu_w = (numConstants > 4) ? constants[4] : 8e-4;
    const PetscScalar rho_w = (numConstants > 6) ? constants[6] : 1100.0;
    const PetscScalar grav = (numConstants > 16) ? constants[16] : 9.81;

    double krw = computeKrw(Sw, constants);

    for (PetscInt d = 0; d < dim; ++d) {
        PetscScalar dP_dd = u_x[uOff_x[0] + d];

        double k_d;
        if (d == 0) k_d = kx;
        else if (d == 1) k_d = ky;
        else k_d = kz;

        double mob_w_d = k_d * krw / mu_w;

        PetscScalar grav_w = (d == dim - 1) ? rho_w * grav : 0.0;

        f1[d] = mob_w_d * (dP_dd - grav_w);
    }
}

// ============================================================================
// PetscDS Jacobian: (pressure, pressure) block
// ============================================================================

void MultiphaseFlowKernel::g0_pp(
    PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[],
    const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, PetscReal u_tShift, const PetscReal x[],
    PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
{
    (void)dim; (void)Nf; (void)NfAux;
    (void)uOff_x; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)x;

    const PetscScalar Sw  = u[uOff[1]];
    const PetscScalar phi = (numConstants > 0) ? constants[0] : 0.15;
    const PetscScalar cw  = (numConstants > 8) ? constants[8] : 4e-10;
    const PetscScalar cg  = (numConstants > 9) ? constants[9] : 1e-8;

    PetscScalar Sg = 1.0 - Sw;
    PetscScalar ct_eff = Sw * cw + Sg * cg;

    // ∂f0_p/∂P = φ * ct_eff * shift
    g0[0] = phi * ct_eff * u_tShift;
}

void MultiphaseFlowKernel::g3_pp(
    PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[],
    const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, PetscReal u_tShift, const PetscReal x[],
    PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
{
    (void)Nf; (void)NfAux;
    (void)uOff_x; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)u_tShift; (void)x;

    const PetscScalar Sw = u[uOff[1]];

    const PetscScalar kx = (numConstants > 1) ? constants[1] : 200e-15;
    const PetscScalar ky = (numConstants > 2) ? constants[2] : 200e-15;
    const PetscScalar kz = (numConstants > 3) ? constants[3] : 100e-15;
    const PetscScalar mu_w = (numConstants > 4) ? constants[4] : 8e-4;
    const PetscScalar mu_g = (numConstants > 5) ? constants[5] : 5e-5;

    double krw = computeKrw(Sw, constants);
    double krg = computeKrg(Sw, constants);

    // Total mobility tensor (diagonal for axis-aligned permeability)
    double k_vals[3] = {PetscRealPart(kx), PetscRealPart(ky), PetscRealPart(kz)};
    double mob_total[3];
    for (int d = 0; d < dim; ++d) {
        double k_d = (d < 3) ? k_vals[d] : k_vals[0];
        mob_total[d] = k_d * (krw / PetscRealPart(mu_w) + krg / PetscRealPart(mu_g));
    }

    // g3 is dim×dim tensor: ∂f1_i/∂(∂P/∂x_j)
    for (PetscInt i = 0; i < dim; ++i) {
        for (PetscInt j = 0; j < dim; ++j) {
            g3[i * dim + j] = (i == j) ? mob_total[i] : 0.0;
        }
    }
}

// ============================================================================
// PetscDS Jacobian: (pressure, saturation) cross-coupling
// ============================================================================

void MultiphaseFlowKernel::g0_ps(
    PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[],
    const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, PetscReal u_tShift, const PetscReal x[],
    PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
{
    (void)dim; (void)Nf; (void)NfAux;
    (void)uOff_x; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)x;

    const PetscScalar P_t = u_t[uOff[0]];
    const PetscScalar phi = (numConstants > 0) ? constants[0] : 0.15;
    const PetscScalar cw  = (numConstants > 8) ? constants[8] : 4e-10;
    const PetscScalar cg  = (numConstants > 9) ? constants[9] : 1e-8;

    // ∂f0_p/∂Sw = φ * (cw - cg) * dP/dt * shift_factor
    // + φ * (ρw - ρg) * 1e-9 * shift  (from density-difference coupling)
    g0[0] = phi * (cw - cg) * P_t + phi * u_tShift * 0.0;
    // The second term is zero because f0_p's dSw/dt term differentiates to
    // phi * (rho_w - rho_g) * 1e-9 * u_tShift, which is negligible
    // compared to the compressibility coupling.
    // For completeness:
    const PetscScalar rho_w = (numConstants > 6) ? constants[6] : 1100.0;
    const PetscScalar rho_g = (numConstants > 7) ? constants[7] : 700.0;
    g0[0] += phi * (rho_w - rho_g) * 1e-9 * u_tShift;
}

// ============================================================================
// PetscDS Jacobian: (saturation, saturation)
// ============================================================================

void MultiphaseFlowKernel::g0_ss(
    PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[],
    const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, PetscReal u_tShift, const PetscReal x[],
    PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
{
    (void)dim; (void)Nf; (void)NfAux;
    (void)uOff; (void)uOff_x; (void)u; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)x;

    const PetscScalar phi = (numConstants > 0) ? constants[0] : 0.15;

    // ∂f0_s/∂Sw = φ * shift
    g0[0] = phi * u_tShift;
}

// ============================================================================
// PetscDS Jacobian: (saturation, pressure) — water mobility tensor
// ============================================================================

void MultiphaseFlowKernel::g3_sp(
    PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[],
    const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, PetscReal u_tShift, const PetscReal x[],
    PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
{
    (void)Nf; (void)NfAux;
    (void)uOff_x; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)u_tShift; (void)x;

    const PetscScalar Sw = u[uOff[1]];

    const PetscScalar kx = (numConstants > 1) ? constants[1] : 200e-15;
    const PetscScalar ky = (numConstants > 2) ? constants[2] : 200e-15;
    const PetscScalar kz = (numConstants > 3) ? constants[3] : 100e-15;
    const PetscScalar mu_w = (numConstants > 4) ? constants[4] : 8e-4;

    double krw = computeKrw(Sw, constants);

    double k_vals[3] = {PetscRealPart(kx), PetscRealPart(ky), PetscRealPart(kz)};

    // ∂f1_s/∂(∇P) = water mobility tensor
    for (PetscInt i = 0; i < dim; ++i) {
        for (PetscInt j = 0; j < dim; ++j) {
            if (i == j) {
                double k_d = (i < 3) ? k_vals[i] : k_vals[0];
                g3[i * dim + j] = k_d * krw / PetscRealPart(mu_w);
            } else {
                g3[i * dim + j] = 0.0;
            }
        }
    }
}

} // namespace FSRM
