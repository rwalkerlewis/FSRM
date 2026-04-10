#include "numerics/PetscFEPoroelasticity.hpp"

namespace FSRM {
namespace PetscFEPoroelasticity {

// =============================================================================
// Pressure equation (Field 0)
// Biot's diffusion equation:
//   (1/M)*dp/dt + alpha*div(du/dt) + div((k/mu)*grad(p)) = 0
// =============================================================================

void f0_pressure(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                 const PetscInt uOff[], const PetscInt uOff_x[],
                 const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscInt aOff[],
                 const PetscInt aOff_x[], const PetscScalar a[],
                 const PetscScalar a_x[], const PetscScalar a_t[],
                 PetscReal t,
                 const PetscReal x[], PetscInt numConstants,
                 const PetscScalar constants[], PetscScalar f0[]) {
    PetscFunctionBeginUser;

    // Suppress unused parameter warnings
    (void)Nf; (void)NfAux;
    (void)u; (void)u_x;
    (void)aOff; (void)aOff_x;
    (void)a; (void)a_x; (void)a_t;
    (void)t; (void)x;

    // Constants: [22]=biot_coefficient, [23]=1/M (inverse Biot modulus)
    const PetscScalar alpha = (numConstants > 22) ? constants[22] : 1.0;
    const PetscScalar M_inv = (numConstants > 23) ? constants[23] : 1e-10;

    // Field 0: pressure (scalar)
    const PetscScalar p_t = u_t[uOff[0]];

    // Field 1: displacement (vector, 3 components)
    // Compute div(u_t) = sum_c du_c/dx_c over time derivative
    PetscScalar div_u_t = 0.0;
    for (PetscInt c = 0; c < dim; ++c) {
        div_u_t += u_t[uOff[1] + c];  // This represents du_c/dt at a point
    }

    // Note: In the full implementation, div(u_t) should actually be computed from
    // the gradient of u_t. However, PETSc's weak form integrates against test functions,
    // and the coupling term alpha*div(du/dt) appears in the strong form.
    // In the weak form, after integration by parts, this becomes:
    //   int( (1/M)*p_t * v + alpha*(du/dt):grad(v) ) dx
    // The alpha*div(du/dt) coupling goes into f1 via integration by parts:
    //   int(alpha*div(u_t)*v) = -int(alpha*u_t . grad(v))
    // So f1[d] -= alpha*u_t[d], with Jacobian in g2_pu.

    // For the f0 term (no integration by parts):
    // f0 = (1/M)*dp/dt
    f0[0] = M_inv * p_t;

    PetscFunctionReturnVoid();
}

void f1_pressure(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                 const PetscInt uOff[], const PetscInt uOff_x[],
                 const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscInt aOff[],
                 const PetscInt aOff_x[], const PetscScalar a[],
                 const PetscScalar a_x[], const PetscScalar a_t[],
                 PetscReal t,
                 const PetscReal x[], PetscInt numConstants,
                 const PetscScalar constants[], PetscScalar f1[]) {
    PetscFunctionBeginUser;

    // Suppress unused parameter warnings
    (void)Nf; (void)NfAux;
    (void)u;
    (void)aOff; (void)aOff_x;
    (void)a; (void)a_x; (void)a_t;
    (void)t; (void)x;

    // Constants: [4]=permeability, [10]=fluid_viscosity, [22]=biot_coefficient
    const PetscScalar k     = (numConstants > 4) ? constants[4] : 1e-13;  // m^2
    const PetscScalar mu_f  = (numConstants > 10) ? constants[10] : 1e-3;   // Pa*s
    const PetscScalar alpha = (numConstants > 22) ? constants[22] : 1.0;

    // Mobility: k/mu
    const PetscScalar mobility = k / mu_f;

    // Weak form of pressure equation after integration by parts on Biot coupling:
    //   int[(1/M)*p_t * v] + int[(k/mu)*grad(p) - alpha*u_t] . grad(v) = 0
    // f1[d] = (k/mu)*dp/dx_d - alpha*du_d/dt
    for (PetscInt d = 0; d < dim; ++d) {
        f1[d] = mobility * u_x[uOff_x[0] + d];
        if (u_t) {
            f1[d] -= alpha * u_t[uOff[1] + d];
        }
    }

    PetscFunctionReturnVoid();
}

// =============================================================================
// Displacement equation (Field 1)
// Quasi-static momentum balance:
//   -div(sigma - alpha*p*I) = 0
// where sigma = lambda*eps_kk*I + 2*mu*eps is the effective stress
// =============================================================================

void f0_displacement(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                     const PetscInt uOff[], const PetscInt uOff_x[],
                     const PetscScalar u[], const PetscScalar u_t[],
                     const PetscScalar u_x[], const PetscInt aOff[],
                     const PetscInt aOff_x[], const PetscScalar a[],
                     const PetscScalar a_x[], const PetscScalar a_t[],
                     PetscReal t,
                     const PetscReal x[], PetscInt numConstants,
                     const PetscScalar constants[], PetscScalar f0[]) {
    PetscFunctionBeginUser;

    // Suppress unused parameter warnings
    (void)Nf; (void)NfAux;
    (void)uOff; (void)uOff_x;
    (void)u; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x;
    (void)a; (void)a_x; (void)a_t;
    (void)t; (void)x;
    (void)numConstants; (void)constants;

    // Quasi-static: no body force, no inertia
    for (PetscInt c = 0; c < dim; ++c) {
        f0[c] = 0.0;
    }

    PetscFunctionReturnVoid();
}

void f1_displacement(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                     const PetscInt uOff[], const PetscInt uOff_x[],
                     const PetscScalar u[], const PetscScalar u_t[],
                     const PetscScalar u_x[], const PetscInt aOff[],
                     const PetscInt aOff_x[], const PetscScalar a[],
                     const PetscScalar a_x[], const PetscScalar a_t[],
                     PetscReal t,
                     const PetscReal x[], PetscInt numConstants,
                     const PetscScalar constants[], PetscScalar f1[]) {
    PetscFunctionBeginUser;

    // Suppress unused parameter warnings
    (void)Nf; (void)NfAux;
    (void)u_t;
    (void)aOff; (void)aOff_x;
    (void)a; (void)a_x; (void)a_t;
    (void)t; (void)x;

    // Constants: [0]=lambda, [1]=mu, [22]=biot_coefficient
    const PetscScalar lambda = (numConstants > 0) ? constants[0] : 1e10;
    const PetscScalar mu     = (numConstants > 1) ? constants[1] : 1e10;
    const PetscScalar alpha  = (numConstants > 22) ? constants[22] : 1.0;

    // Field 0: pressure (scalar)
    const PetscScalar p = u[uOff[0]];

    // Field 1: displacement (vector with 3 components)
    // For vector field: u_x[uOff_x[1] + c*dim + d] = d(u_c)/dx_d

    // Compute volumetric strain: eps_kk = du_x/dx + du_y/dy + du_z/dz
    PetscScalar eps_kk = 0.0;
    for (PetscInt d = 0; d < dim; ++d) {
        eps_kk += u_x[uOff_x[1] + d*dim + d];  // d(u_d)/dx_d
    }

    // Compute stress tensor and populate f1
    // f1[c*dim + d] = sigma_{c,d} - alpha*p*delta_{c,d}
    // where sigma_{c,d} = lambda*eps_kk*delta_{c,d} + 2*mu*eps_{c,d}
    // and eps_{c,d} = 0.5*(du_c/dx_d + du_d/dx_c)

    for (PetscInt c = 0; c < dim; ++c) {
        for (PetscInt d = 0; d < dim; ++d) {
            // Strain component: eps_{c,d} = 0.5*(du_c/dx_d + du_d/dx_c)
            PetscScalar eps_cd = 0.5 * (u_x[uOff_x[1] + c*dim + d] +
                                        u_x[uOff_x[1] + d*dim + c]);

            // Effective stress: sigma_{c,d} = lambda*eps_kk*delta_{c,d} + 2*mu*eps_{c,d}
            PetscScalar sigma_cd = 2.0 * mu * eps_cd;
            if (c == d) {
                sigma_cd += lambda * eps_kk;
            }

            // Total stress includes Biot coupling: sigma_total = sigma_eff - alpha*p*I
            // f1[c*dim + d] represents the flux term in weak form
            if (c == d) {
                f1[c*dim + d] = sigma_cd - alpha * p;
            } else {
                f1[c*dim + d] = sigma_cd;
            }
        }
    }

    PetscFunctionReturnVoid();
}

// =============================================================================
// Jacobian callbacks
// =============================================================================

// g0_pp: d(f0_pressure)/d(p) = (1/M)*shift
void g0_pp(PetscInt dim, PetscInt Nf, PetscInt NfAux,
           const PetscInt uOff[], const PetscInt uOff_x[],
           const PetscScalar u[], const PetscScalar u_t[],
           const PetscScalar u_x[], const PetscInt aOff[],
           const PetscInt aOff_x[], const PetscScalar a[],
           const PetscScalar a_x[], const PetscScalar a_t[],
           PetscReal t, PetscReal u_tShift,
           const PetscReal x[], PetscInt numConstants,
           const PetscScalar constants[], PetscScalar g0[]) {
    PetscFunctionBeginUser;

    // Suppress unused parameter warnings
    (void)dim; (void)Nf; (void)NfAux;
    (void)uOff; (void)uOff_x;
    (void)u; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x;
    (void)a; (void)a_x; (void)a_t;
    (void)t; (void)x;

    // Constants: [23]=1/M (inverse Biot modulus)
    const PetscScalar M_inv = (numConstants > 23) ? constants[23] : 1e-10;

    // g0 = d(f0)/d(u) * shift
    // Since f0 = (1/M)*p_t, we have d(f0)/d(p) = (1/M)*d(p_t)/d(p) = (1/M)*shift
    g0[0] = M_inv * u_tShift;

    PetscFunctionReturnVoid();
}

// g3_pp: d(f1_pressure)/d(grad_p) = (k/mu)*delta_{dd}
void g3_pp(PetscInt dim, PetscInt Nf, PetscInt NfAux,
           const PetscInt uOff[], const PetscInt uOff_x[],
           const PetscScalar u[], const PetscScalar u_t[],
           const PetscScalar u_x[], const PetscInt aOff[],
           const PetscInt aOff_x[], const PetscScalar a[],
           const PetscScalar a_x[], const PetscScalar a_t[],
           PetscReal t, PetscReal u_tShift,
           const PetscReal x[], PetscInt numConstants,
           const PetscScalar constants[], PetscScalar g3[]) {
    PetscFunctionBeginUser;

    // Suppress unused parameter warnings
    (void)Nf; (void)NfAux;
    (void)uOff; (void)uOff_x;
    (void)u; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x;
    (void)a; (void)a_x; (void)a_t;
    (void)t; (void)u_tShift; (void)x;

    // Constants: [4]=permeability, [10]=fluid_viscosity (mu_w)
    const PetscScalar k    = (numConstants > 4) ? constants[4] : 1e-13;
    const PetscScalar mu_f = (numConstants > 10) ? constants[10] : 1e-3;

    const PetscScalar mobility = k / mu_f;

    // g3[d1*dim + d2] = d(f1[d1])/d(grad_p[d2])
    // Since f1[d] = mobility * grad_p[d], we have:
    // d(f1[d1])/d(grad_p[d2]) = mobility * delta_{d1,d2}

    for (PetscInt d1 = 0; d1 < dim; ++d1) {
        for (PetscInt d2 = 0; d2 < dim; ++d2) {
            g3[d1*dim + d2] = (d1 == d2) ? mobility : 0.0;
        }
    }

    PetscFunctionReturnVoid();
}

// g2_pu: d(f1_pressure)/d(u_t) = -alpha*shift*delta_{dc}
// Biot coupling in pressure residual: f1[d] -= alpha*u_t[d]
// Jacobian: d(f1[d])/d(u_t_c) * shift = -alpha*shift*delta_{dc}
void g2_pu(PetscInt dim, PetscInt Nf, PetscInt NfAux,
           const PetscInt uOff[], const PetscInt uOff_x[],
           const PetscScalar u[], const PetscScalar u_t[],
           const PetscScalar u_x[], const PetscInt aOff[],
           const PetscInt aOff_x[], const PetscScalar a[],
           const PetscScalar a_x[], const PetscScalar a_t[],
           PetscReal t, PetscReal u_tShift,
           const PetscReal x[], PetscInt numConstants,
           const PetscScalar constants[], PetscScalar g2[]) {
    PetscFunctionBeginUser;

    // Suppress unused parameter warnings
    (void)Nf; (void)NfAux;
    (void)uOff; (void)uOff_x;
    (void)u; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x;
    (void)a; (void)a_x; (void)a_t;
    (void)t; (void)x;

    // Constants: [22]=biot_coefficient
    const PetscScalar alpha = (numConstants > 22) ? constants[22] : 1.0;

    // g2[d*dim + c] = shift * d(f1_pressure[d])/d(u_t_c)
    // f1[d] = (k/mu)*dp/dx_d - alpha*u_t[d]
    // d(f1[d])/d(u_t_c) = -alpha*delta_{dc}
    // With shift: g2[d*dim+c] = -alpha*u_tShift*delta_{dc}

    for (PetscInt d = 0; d < dim; ++d) {
        for (PetscInt c = 0; c < dim; ++c) {
            g2[d*dim + c] = (d == c) ? -alpha * u_tShift : 0.0;
        }
    }

    PetscFunctionReturnVoid();
}

// g2_up: d(f1_displacement)/d(p) = -alpha*delta_{cd}
// Coupling from pressure to displacement momentum balance
void g2_up(PetscInt dim, PetscInt Nf, PetscInt NfAux,
           const PetscInt uOff[], const PetscInt uOff_x[],
           const PetscScalar u[], const PetscScalar u_t[],
           const PetscScalar u_x[], const PetscInt aOff[],
           const PetscInt aOff_x[], const PetscScalar a[],
           const PetscScalar a_x[], const PetscScalar a_t[],
           PetscReal t, PetscReal u_tShift,
           const PetscReal x[], PetscInt numConstants,
           const PetscScalar constants[], PetscScalar g2[]) {
    PetscFunctionBeginUser;

    // Suppress unused parameter warnings
    (void)Nf; (void)NfAux;
    (void)uOff; (void)uOff_x;
    (void)u; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x;
    (void)a; (void)a_x; (void)a_t;
    (void)t; (void)u_tShift; (void)x;

    // Constants: [22]=biot_coefficient
    const PetscScalar alpha = (numConstants > 22) ? constants[22] : 1.0;

    // g2[c*dim + d] = d(f1_displacement[c*dim + d])/d(p)
    // Since f1[c*dim + d] includes -alpha*p*delta_{c,d}, we have:
    // d(f1[c*dim + d])/d(p) = -alpha*delta_{c,d}

    for (PetscInt c = 0; c < dim; ++c) {
        for (PetscInt d = 0; d < dim; ++d) {
            g2[c*dim + d] = (c == d) ? -alpha : 0.0;
        }
    }

    PetscFunctionReturnVoid();
}

// g3_uu: d(f1_displacement)/d(grad_u)
// Elasticity stiffness tensor
void g3_uu(PetscInt dim, PetscInt Nf, PetscInt NfAux,
           const PetscInt uOff[], const PetscInt uOff_x[],
           const PetscScalar u[], const PetscScalar u_t[],
           const PetscScalar u_x[], const PetscInt aOff[],
           const PetscInt aOff_x[], const PetscScalar a[],
           const PetscScalar a_x[], const PetscScalar a_t[],
           PetscReal t, PetscReal u_tShift,
           const PetscReal x[], PetscInt numConstants,
           const PetscScalar constants[], PetscScalar g3[]) {
    PetscFunctionBeginUser;

    // Suppress unused parameter warnings
    (void)Nf; (void)NfAux;
    (void)uOff; (void)uOff_x;
    (void)u; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x;
    (void)a; (void)a_x; (void)a_t;
    (void)t; (void)u_tShift; (void)x;

    // Constants: [0]=lambda, [1]=mu
    const PetscScalar lambda = (numConstants > 0) ? constants[0] : 1e10;
    const PetscScalar mu     = (numConstants > 1) ? constants[1] : 1e10;

    // Jacobian of elasticity: d(sigma_{i,d1})/d(du_j/dx_{d2})
    // = C_{i,d1,j,d2}
    // = lambda*delta_{i,d1}*delta_{j,d2} + mu*(delta_{ij}*delta_{d1,d2} + delta_{i,d2}*delta_{d1,j})
    //
    // g3[i*dim*dim*dim + j*dim*dim + d1*dim + d2] = C_{i,d1,j,d2}
    //
    // This is the standard isotropic elasticity stiffness tensor.

    for (PetscInt i = 0; i < dim; ++i) {
        for (PetscInt j = 0; j < dim; ++j) {
            for (PetscInt d1 = 0; d1 < dim; ++d1) {
                for (PetscInt d2 = 0; d2 < dim; ++d2) {
                    PetscScalar val = 0.0;

                    // lambda * delta_{i,d1} * delta_{j,d2}
                    if (i == d1 && j == d2) {
                        val += lambda;
                    }

                    // mu * delta_{i,j} * delta_{d1,d2}
                    if (i == j && d1 == d2) {
                        val += mu;
                    }

                    // mu * delta_{i,d2} * delta_{d1,j}
                    if (i == d2 && d1 == j) {
                        val += mu;
                    }

                    g3[i*dim*dim*dim + j*dim*dim + d1*dim + d2] = val;
                }
            }
        }
    }

    PetscFunctionReturnVoid();
}

} // namespace PetscFEPoroelasticity
} // namespace FSRM
