#include "numerics/PetscFEElasticity.hpp"

namespace FSRM {
namespace PetscFEElasticity {

// =============================================================================
// Elastostatics: -div(sigma) = f_body
// where sigma_{ij} = lambda * eps_kk * delta_ij + 2*mu * eps_ij
// and eps_ij = 0.5*(du_i/dx_j + du_j/dx_i)
// =============================================================================

void f0_elastostatics(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                      const PetscInt uOff[], const PetscInt uOff_x[],
                      const PetscScalar u[], const PetscScalar u_t[],
                      const PetscScalar u_x[], const PetscInt aOff[],
                      const PetscInt aOff_x[], const PetscScalar a[],
                      const PetscScalar a_x[], const PetscScalar a_t[],
                      PetscReal t,
                      const PetscReal x[], PetscInt numConstants,
                      const PetscScalar constants[], PetscScalar f0[]) {
    // Suppress unused parameter warnings - part of PETSc pointwise function interface
    (void)dim; (void)Nf; (void)NfAux;
    (void)uOff; (void)uOff_x;
    (void)u; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x;
    (void)a; (void)a_x; (void)a_t;
    (void)t; (void)x;
    (void)numConstants; (void)constants;

    // Quasi-static elasticity: no body force, no inertia
    // f0[c] = 0 for all components c
    for (PetscInt c = 0; c < dim; ++c) {
        f0[c] = 0.0;
    }
}

void f1_elastostatics(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                      const PetscInt uOff[], const PetscInt uOff_x[],
                      const PetscScalar u[], const PetscScalar u_t[],
                      const PetscScalar u_x[], const PetscInt aOff[],
                      const PetscInt aOff_x[], const PetscScalar a[],
                      const PetscScalar a_x[], const PetscScalar a_t[],
                      PetscReal t,
                      const PetscReal x[], PetscInt numConstants,
                      const PetscScalar constants[], PetscScalar f1[]) {
    // Suppress unused parameter warnings
    (void)Nf; (void)NfAux;
    (void)uOff; (void)u; (void)u_t;
    (void)aOff; (void)aOff_x;
    (void)a; (void)a_x; (void)a_t;
    (void)t; (void)x;

    // Constants: [0]=lambda, [1]=mu
    const PetscScalar lambda = (numConstants > 0) ? constants[0] : 1e10;
    const PetscScalar mu     = (numConstants > 1) ? constants[1] : 1e10;

    // For a vector field with nc components, gradient layout is:
    // u_x[uOff_x[f] + c*dim + d] = d(u_c)/dx_d
    // where c is component index, d is spatial direction

    // Compute volumetric strain: eps_kk = du_x/dx + du_y/dy + du_z/dz
    PetscScalar eps_kk = 0.0;
    for (PetscInt d = 0; d < dim; ++d) {
        eps_kk += u_x[uOff_x[0] + d*dim + d];  // d(u_d)/dx_d
    }

    // Compute stress tensor: sigma_{ij} = lambda*eps_kk*delta_ij + 2*mu*eps_ij
    // where eps_ij = 0.5*(du_i/dx_j + du_j/dx_i)
    //
    // f1[c*dim + d] = sigma_{c,d} for component c, direction d

    for (PetscInt c = 0; c < dim; ++c) {
        for (PetscInt d = 0; d < dim; ++d) {
            // Strain component: eps_{c,d} = 0.5*(du_c/dx_d + du_d/dx_c)
            PetscScalar eps_cd = 0.5 * (u_x[uOff_x[0] + c*dim + d] +
                                        u_x[uOff_x[0] + d*dim + c]);

            // Stress: sigma_{c,d} = lambda*eps_kk*delta_{c,d} + 2*mu*eps_{c,d}
            PetscScalar sigma_cd = 2.0 * mu * eps_cd;
            if (c == d) {
                sigma_cd += lambda * eps_kk;
            }

            f1[c*dim + d] = sigma_cd;
        }
    }
}

void g3_elastostatics(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                      const PetscInt uOff[], const PetscInt uOff_x[],
                      const PetscScalar u[], const PetscScalar u_t[],
                      const PetscScalar u_x[], const PetscInt aOff[],
                      const PetscInt aOff_x[], const PetscScalar a[],
                      const PetscScalar a_x[], const PetscScalar a_t[],
                      PetscReal t, PetscReal u_tShift,
                      const PetscReal x[], PetscInt numConstants,
                      const PetscScalar constants[], PetscScalar g3[]) {
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
}

// =============================================================================
// Elastodynamics: rho*u_tt - div(sigma) = f_body
// Adds inertial term
// =============================================================================

void f0_elastodynamics(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[],
                       const PetscScalar u[], const PetscScalar u_t[],
                       const PetscScalar u_x[], const PetscInt aOff[],
                       const PetscInt aOff_x[], const PetscScalar a[],
                       const PetscScalar a_x[], const PetscScalar a_t[],
                       PetscReal t,
                       const PetscReal x[], PetscInt numConstants,
                       const PetscScalar constants[], PetscScalar f0[]) {
    // Suppress unused parameter warnings
    (void)Nf; (void)NfAux;
    (void)uOff_x; (void)u; (void)u_x;
    (void)aOff; (void)aOff_x;
    (void)a; (void)a_x; (void)a_t;
    (void)t; (void)x;

    // Constants: [2]=rho (density)
    const PetscScalar rho = (numConstants > 2) ? constants[2] : 2700.0;

    // Inertial term: rho * u_tt
    // In PETSc implicit time stepping (TSBEULER, etc), u_t is the time derivative
    // For second-order systems, TS expects the residual F = M*u_tt + K*u - f
    // where u_tt is handled internally by the time stepper
    //
    // For now, we use first-order formulation with velocity as auxiliary variable,
    // or rely on TS's TSALPHA2 for proper second-order dynamics.
    // In the simple case with implicit Euler, f0[c] = rho * u_t[c] gives mass*velocity term.

    for (PetscInt c = 0; c < dim; ++c) {
        f0[c] = rho * u_t[uOff[0] + c];
    }
}

void f1_elastodynamics(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[],
                       const PetscScalar u[], const PetscScalar u_t[],
                       const PetscScalar u_x[], const PetscInt aOff[],
                       const PetscInt aOff_x[], const PetscScalar a[],
                       const PetscScalar a_x[], const PetscScalar a_t[],
                       PetscReal t,
                       const PetscReal x[], PetscInt numConstants,
                       const PetscScalar constants[], PetscScalar f1[]) {
    // Same as elastostatics - stress tensor
    f1_elastostatics(dim, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
                     aOff, aOff_x, a, a_x, a_t, t, x,
                     numConstants, constants, f1);
}

void g0_elastodynamics(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[],
                       const PetscScalar u[], const PetscScalar u_t[],
                       const PetscScalar u_x[], const PetscInt aOff[],
                       const PetscInt aOff_x[], const PetscScalar a[],
                       const PetscScalar a_x[], const PetscScalar a_t[],
                       PetscReal t, PetscReal u_tShift,
                       const PetscReal x[], PetscInt numConstants,
                       const PetscScalar constants[], PetscScalar g0[]) {
    // Suppress unused parameter warnings
    (void)Nf; (void)NfAux;
    (void)uOff; (void)uOff_x;
    (void)u; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x;
    (void)a; (void)a_x; (void)a_t;
    (void)t; (void)x;

    // Constants: [2]=rho (density)
    const PetscScalar rho = (numConstants > 2) ? constants[2] : 2700.0;

    // Mass matrix: d(f0)/d(u_t) = rho * u_tShift * I
    // where u_tShift is the time integration parameter (e.g., 1/dt for backward Euler)
    //
    // g0[i*dim + j] = d(f0[i])/d(u_t[j]) * u_tShift
    //               = rho * delta_{i,j} * u_tShift

    for (PetscInt i = 0; i < dim; ++i) {
        for (PetscInt j = 0; j < dim; ++j) {
            g0[i*dim + j] = (i == j) ? (rho * u_tShift) : 0.0;
        }
    }
}

void g3_elastodynamics(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[],
                       const PetscScalar u[], const PetscScalar u_t[],
                       const PetscScalar u_x[], const PetscInt aOff[],
                       const PetscInt aOff_x[], const PetscScalar a[],
                       const PetscScalar a_x[], const PetscScalar a_t[],
                       PetscReal t, PetscReal u_tShift,
                       const PetscReal x[], PetscInt numConstants,
                       const PetscScalar constants[], PetscScalar g3[]) {
    // Same as elastostatics - stiffness tensor
    g3_elastostatics(dim, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
                     aOff, aOff_x, a, a_x, a_t, t, u_tShift, x,
                     numConstants, constants, g3);
}

// =============================================================================
// I2-form callbacks for second-order-in-time dynamics (TSALPHA2)
// For use with TSSetI2Function / TSSetI2Jacobian
// =============================================================================

void f0_elastodynamics_I2(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                          const PetscInt uOff[], const PetscInt uOff_x[],
                          const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_tt[],
                          const PetscScalar u_x[], const PetscInt aOff[],
                          const PetscInt aOff_x[], const PetscScalar a[],
                          const PetscScalar a_x[], const PetscScalar a_t[],
                          PetscReal t,
                          const PetscReal x[], PetscInt numConstants,
                          const PetscScalar constants[], PetscScalar f0[]) {
    // Suppress unused parameter warnings
    (void)Nf; (void)NfAux;
    (void)uOff_x; (void)u; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x;
    (void)a; (void)a_x; (void)a_t;
    (void)t; (void)x;

    // Constants: [2]=rho (density)
    const PetscScalar rho = (numConstants > 2) ? constants[2] : 2700.0;

    // Mass term: rho * u_tt
    // F(t, U, U_t, U_tt) = M*U_tt + K*U - f
    // f0[c] = rho * u_tt[c] for each component c
    for (PetscInt c = 0; c < dim; ++c) {
        f0[c] = rho * u_tt[uOff[0] + c];
    }
}

void f1_elastodynamics_I2(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                          const PetscInt uOff[], const PetscInt uOff_x[],
                          const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_tt[],
                          const PetscScalar u_x[], const PetscInt aOff[],
                          const PetscInt aOff_x[], const PetscScalar a[],
                          const PetscScalar a_x[], const PetscScalar a_t[],
                          PetscReal t,
                          const PetscReal x[], PetscInt numConstants,
                          const PetscScalar constants[], PetscScalar f1[]) {
    // Suppress unused parameter warnings
    (void)Nf; (void)NfAux;
    (void)uOff; (void)u_t; (void)u_tt;
    (void)aOff; (void)aOff_x;
    (void)a; (void)a_x; (void)a_t;
    (void)t; (void)x;

    // Constants: [0]=lambda, [1]=mu
    const PetscScalar lambda = (numConstants > 0) ? constants[0] : 1e10;
    const PetscScalar mu     = (numConstants > 1) ? constants[1] : 1e10;

    // Stiffness term: -div(sigma)
    // This is identical to elastostatics f1 - compute stress tensor
    // Compute volumetric strain: eps_kk = du_x/dx + du_y/dy + du_z/dz
    PetscScalar eps_kk = 0.0;
    for (PetscInt d = 0; d < dim; ++d) {
        eps_kk += u_x[uOff_x[0] + d*dim + d];  // d(u_d)/dx_d
    }

    // Compute stress tensor: sigma_{ij} = lambda*eps_kk*delta_ij + 2*mu*eps_ij
    for (PetscInt c = 0; c < dim; ++c) {
        for (PetscInt d = 0; d < dim; ++d) {
            // Strain component: eps_{c,d} = 0.5*(du_c/dx_d + du_d/dx_c)
            PetscScalar eps_cd = 0.5 * (u_x[uOff_x[0] + c*dim + d] +
                                        u_x[uOff_x[0] + d*dim + c]);

            // Stress: sigma_{c,d} = lambda*eps_kk*delta_{c,d} + 2*mu*eps_{c,d}
            PetscScalar sigma_cd = 2.0 * mu * eps_cd;
            if (c == d) {
                sigma_cd += lambda * eps_kk;
            }

            f1[c*dim + d] = sigma_cd;
        }
    }
}

void g0_elastodynamics_I2(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                          const PetscInt uOff[], const PetscInt uOff_x[],
                          const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_tt[],
                          const PetscScalar u_x[], const PetscInt aOff[],
                          const PetscInt aOff_x[], const PetscScalar a[],
                          const PetscScalar a_x[], const PetscScalar a_t[],
                          PetscReal t, PetscReal u_tShift, PetscReal u_ttShift,
                          const PetscReal x[], PetscInt numConstants,
                          const PetscScalar constants[], PetscScalar g0[]) {
    // Suppress unused parameter warnings
    (void)Nf; (void)NfAux;
    (void)uOff; (void)uOff_x;
    (void)u; (void)u_t; (void)u_tt; (void)u_x;
    (void)aOff; (void)aOff_x;
    (void)a; (void)a_x; (void)a_t;
    (void)t; (void)u_tShift; (void)x;

    // Constants: [2]=rho (density)
    const PetscScalar rho = (numConstants > 2) ? constants[2] : 2700.0;

    // Mass matrix Jacobian: d(f0)/d(u_tt) = rho * u_ttShift * I
    // g0[i*dim + j] = d(f0[i])/d(u_tt[j]) * u_ttShift
    //               = rho * delta_{i,j} * u_ttShift
    for (PetscInt i = 0; i < dim; ++i) {
        for (PetscInt j = 0; j < dim; ++j) {
            g0[i*dim + j] = (i == j) ? (rho * u_ttShift) : 0.0;
        }
    }
}

void g3_elastodynamics_I2(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                          const PetscInt uOff[], const PetscInt uOff_x[],
                          const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_tt[],
                          const PetscScalar u_x[], const PetscInt aOff[],
                          const PetscInt aOff_x[], const PetscScalar a[],
                          const PetscScalar a_x[], const PetscScalar a_t[],
                          PetscReal t, PetscReal u_tShift, PetscReal u_ttShift,
                          const PetscReal x[], PetscInt numConstants,
                          const PetscScalar constants[], PetscScalar g3[]) {
    // Suppress unused parameter warnings
    (void)Nf; (void)NfAux;
    (void)uOff; (void)uOff_x;
    (void)u; (void)u_t; (void)u_tt; (void)u_x;
    (void)aOff; (void)aOff_x;
    (void)a; (void)a_x; (void)a_t;
    (void)t; (void)u_tShift; (void)u_ttShift; (void)x;

    // Constants: [0]=lambda, [1]=mu
    const PetscScalar lambda = (numConstants > 0) ? constants[0] : 1e10;
    const PetscScalar mu     = (numConstants > 1) ? constants[1] : 1e10;

    // Stiffness tensor Jacobian (same as elastostatics)
    // g3[i*dim*dim*dim + j*dim*dim + d1*dim + d2] = C_{i,d1,j,d2}
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
}

} // namespace PetscFEElasticity
} // namespace FSRM
