/**
 * @file PetscFEThermal.cpp
 * @brief PetscDS pointwise callbacks for thermal diffusion and thermoelastic coupling
 */

#include "numerics/PetscFEThermal.hpp"
#include "numerics/PetscFEElasticityAux.hpp"

namespace FSRM {
namespace PetscFEThermal {

// ---------------------------------------------------------------------------
// Helper: safe constant read with bounds check
// ---------------------------------------------------------------------------
static inline PetscReal getConst(PetscInt numConstants, const PetscScalar c[],
                                 PetscInt idx, PetscReal fallback) {
    return (numConstants > idx) ? PetscRealPart(c[idx]) : fallback;
}

// =========================================================================
// f0_thermal: rho * Cp * dT/dt
// =========================================================================
void f0_thermal(PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[],
    const PetscScalar u_x[], const PetscInt aOff[],
    const PetscInt aOff_x[], const PetscScalar a[],
    const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, const PetscReal x[], PetscInt numConstants,
    const PetscScalar constants[], PetscScalar f0[])
{
    (void)dim; (void)Nf; (void)NfAux;
    (void)uOff_x; (void)u; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)x;

    // Thermal field index in solution
    PetscInt thermal_idx = static_cast<PetscInt>(getConst(numConstants, constants,
                                                          THERMAL_CONST_FIELD_IDX, static_cast<PetscReal>(Nf - 1)));
    PetscReal rho = getConst(numConstants, constants, 2, 2650.0);
    PetscReal Cp  = getConst(numConstants, constants, THERMAL_CONST_CP, 800.0);

    // dT/dt
    f0[0] = rho * Cp * u_t[uOff[thermal_idx]];
}

// =========================================================================
// f1_thermal: kappa * grad(T)
// =========================================================================
void f1_thermal(PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[],
    const PetscScalar u_x[], const PetscInt aOff[],
    const PetscInt aOff_x[], const PetscScalar a[],
    const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, const PetscReal x[], PetscInt numConstants,
    const PetscScalar constants[], PetscScalar f1[])
{
    (void)Nf; (void)NfAux;
    (void)uOff; (void)u; (void)u_t;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)x;

    PetscInt thermal_idx = static_cast<PetscInt>(getConst(numConstants, constants,
                                                          THERMAL_CONST_FIELD_IDX, static_cast<PetscReal>(Nf - 1)));
    PetscReal kappa = getConst(numConstants, constants, THERMAL_CONST_KAPPA, 3.0);

    // grad(T) -- temperature gradient
    for (PetscInt d = 0; d < dim; ++d) {
        f1[d] = kappa * u_x[uOff_x[thermal_idx] + d];
    }
}

// =========================================================================
// g0_thermal_thermal: rho * Cp * shift
// =========================================================================
void g0_thermal_thermal(PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[],
    const PetscScalar u_x[], const PetscInt aOff[],
    const PetscInt aOff_x[], const PetscScalar a[],
    const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, PetscReal u_tShift, const PetscReal x[],
    PetscInt numConstants, const PetscScalar constants[],
    PetscScalar g0[])
{
    (void)dim; (void)Nf; (void)NfAux;
    (void)uOff; (void)uOff_x;
    (void)u; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)x;

    PetscReal rho = getConst(numConstants, constants, 2, 2650.0);
    PetscReal Cp  = getConst(numConstants, constants, THERMAL_CONST_CP, 800.0);

    g0[0] = rho * Cp * u_tShift;
}

// =========================================================================
// g3_thermal_thermal: kappa * delta_ij
// =========================================================================
void g3_thermal_thermal(PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[],
    const PetscScalar u_x[], const PetscInt aOff[],
    const PetscInt aOff_x[], const PetscScalar a[],
    const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, PetscReal u_tShift, const PetscReal x[],
    PetscInt numConstants, const PetscScalar constants[],
    PetscScalar g3[])
{
    (void)Nf; (void)NfAux;
    (void)uOff; (void)uOff_x;
    (void)u; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)u_tShift; (void)x;

    PetscReal kappa = getConst(numConstants, constants, THERMAL_CONST_KAPPA, 3.0);

    // For a scalar field, g3 is dim x dim
    for (PetscInt d = 0; d < dim; ++d) {
        for (PetscInt e = 0; e < dim; ++e) {
            g3[d * dim + e] = (d == e) ? kappa : 0.0;
        }
    }
}

// =========================================================================
// f1_thermoelastic: elastic stress with thermal strain
// sigma_ij = lambda*(eps_kk - (dim)*alpha_T*dT)*delta_ij
//          + 2*mu*(eps_ij - alpha_T*dT*delta_ij)
// where dT = T - T_ref
// =========================================================================
void f1_thermoelastic(PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[],
    const PetscScalar u_x[], const PetscInt aOff[],
    const PetscInt aOff_x[], const PetscScalar a[],
    const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, const PetscReal x[], PetscInt numConstants,
    const PetscScalar constants[], PetscScalar f1[])
{
    (void)Nf; (void)NfAux;
    (void)u; (void)u_t;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)x;

    // Material parameters from constants
    const PetscScalar lambda = (numConstants > 0) ? constants[0] : 1e10;
    const PetscScalar mu     = (numConstants > 1) ? constants[1] : 1e10;
    PetscReal alpha_T = getConst(numConstants, constants, THERMAL_CONST_ALPHA_T, 0.0);
    PetscReal T_ref   = getConst(numConstants, constants, THERMAL_CONST_T_REF, 293.0);
    PetscInt thermal_idx = static_cast<PetscInt>(getConst(numConstants, constants,
                                                          THERMAL_CONST_FIELD_IDX, static_cast<PetscReal>(Nf - 1)));

    // Read temperature from solution
    PetscScalar T_val = u[uOff[thermal_idx]];
    PetscScalar dT = T_val - T_ref;

    // Displacement gradient (field 0 for NONE fluid model)
    // The displacement field index is assumed to be 0 when there are no fluid fields.
    PetscInt disp_off = uOff_x[0];

    // Compute volumetric strain: eps_kk = du_x/dx + du_y/dy + du_z/dz
    PetscScalar eps_kk = 0.0;
    for (PetscInt d = 0; d < dim; ++d) {
        eps_kk += u_x[disp_off + d * dim + d];
    }

    // Thermal volumetric strain: eps_thermal_kk = dim * alpha_T * dT (3D: 3*alpha_T*dT)
    PetscScalar eps_thermal_kk = static_cast<PetscScalar>(dim) * alpha_T * dT;

    // sigma_ij = lambda * (eps_kk - eps_thermal_kk) * delta_ij
    //          + 2 * mu * (eps_ij - alpha_T * dT * delta_ij)
    for (PetscInt c = 0; c < dim; ++c) {
        for (PetscInt d = 0; d < dim; ++d) {
            PetscScalar eps_cd = 0.5 * (u_x[disp_off + c * dim + d] +
                                        u_x[disp_off + d * dim + c]);

            PetscScalar sigma_cd = 2.0 * mu * (eps_cd - ((c == d) ? alpha_T * dT : 0.0));
            if (c == d) {
                sigma_cd += lambda * (eps_kk - eps_thermal_kk);
            }
            f1[c * dim + d] = sigma_cd;
        }
    }
}

// =========================================================================
// f1_thermoelastic_aux: same but reads lambda, mu from aux fields
// =========================================================================
void f1_thermoelastic_aux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[],
    const PetscScalar u_x[], const PetscInt aOff[],
    const PetscInt aOff_x[], const PetscScalar a[],
    const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, const PetscReal x[], PetscInt numConstants,
    const PetscScalar constants[], PetscScalar f1[])
{
    (void)Nf; (void)NfAux;
    (void)u; (void)u_t;
    (void)aOff_x; (void)a_t; (void)a_x;
    (void)t; (void)x;

    // Material parameters from auxiliary fields
    const PetscScalar lambda = a[aOff[AUX_LAMBDA]];
    const PetscScalar mu     = a[aOff[AUX_MU]];
    PetscReal alpha_T = getConst(numConstants, constants, THERMAL_CONST_ALPHA_T, 0.0);
    PetscReal T_ref   = getConst(numConstants, constants, THERMAL_CONST_T_REF, 293.0);
    PetscInt thermal_idx = static_cast<PetscInt>(getConst(numConstants, constants,
                                                          THERMAL_CONST_FIELD_IDX, static_cast<PetscReal>(Nf - 1)));

    // Read temperature from solution
    PetscScalar T_val = u[uOff[thermal_idx]];
    PetscScalar dT = T_val - T_ref;

    // Displacement gradient
    PetscInt disp_off = uOff_x[0];

    // Compute volumetric strain
    PetscScalar eps_kk = 0.0;
    for (PetscInt d = 0; d < dim; ++d) {
        eps_kk += u_x[disp_off + d * dim + d];
    }

    PetscScalar eps_thermal_kk = static_cast<PetscScalar>(dim) * alpha_T * dT;

    for (PetscInt c = 0; c < dim; ++c) {
        for (PetscInt d = 0; d < dim; ++d) {
            PetscScalar eps_cd = 0.5 * (u_x[disp_off + c * dim + d] +
                                        u_x[disp_off + d * dim + c]);

            PetscScalar sigma_cd = 2.0 * mu * (eps_cd - ((c == d) ? alpha_T * dT : 0.0));
            if (c == d) {
                sigma_cd += lambda * (eps_kk - eps_thermal_kk);
            }
            f1[c * dim + d] = sigma_cd;
        }
    }
}

} // namespace PetscFEThermal
} // namespace FSRM
