/**
 * @file CohesiveFaultKernel.cpp
 * @brief PetscDS-compatible pointwise callbacks for cohesive fault cells
 *
 * Implements the fault constraint equations:
 * - LOCKED: displacement jump = 0 (Lagrange multiplier enforces zero slip)
 * - SLIPPING: traction = friction_strength * slip_direction
 *
 * The callbacks follow PETSc's hybrid cell integration convention where
 * the solution values from both sides of the fault are available through
 * the u[] array, and the fault normal is provided via n[].
 */

#include "physics/CohesiveFaultKernel.hpp"
#include "domain/geomechanics/PyLithFault.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace FSRM {

CohesiveFaultKernel::CohesiveFaultKernel()
    : locked_(true), mu_friction_(0.6), friction_model_(nullptr) {}

void CohesiveFaultKernel::setMode(bool is_locked) {
    locked_ = is_locked;
}

void CohesiveFaultKernel::setFrictionCoefficient(double mu_f) {
    mu_friction_ = mu_f;
}

void CohesiveFaultKernel::setTensileStrength(double T_s) {
    tensile_strength_ = T_s;
}

void CohesiveFaultKernel::setPrescribedSlip(double sx, double sy, double sz) {
    prescribed_slip_[0] = sx;
    prescribed_slip_[1] = sy;
    prescribed_slip_[2] = sz;
    prescribed_slip_mode_ = true;
}

void CohesiveFaultKernel::setSlipWeakeningParams(double mu_s, double mu_d, double D_c) {
    mu_static_ = mu_s;
    mu_dynamic_ = mu_d;
    critical_slip_distance_ = D_c;
    use_slip_weakening_ = true;
}

void CohesiveFaultKernel::setFrictionModel(FaultFrictionModel* model) {
    friction_model_ = model;
}

PetscErrorCode CohesiveFaultKernel::registerWithDS(PetscDS prob,
                                                     int displacement_field,
                                                     int lagrange_field) {
    PetscFunctionBeginUser;
    (void)prob;
    (void)displacement_field;
    (void)lagrange_field;

    // Session 4: the PetscDSSetBdResidual registration path is deprecated.
    // The cohesive residual / Jacobian are driven via
    // DMPlexComputeResidualHybridByKey / DMPlexComputeJacobianHybridByKey
    // after explicit PetscWeakFormSetIndexBdResidual / ...BdJacobian
    // registrations (see registerHybridWithDS). Calling this legacy entry
    // point is a no-op so no stale PetscDS callbacks are attached and so
    // double-integration cannot occur alongside the hybrid driver.
    std::fprintf(stderr,
                 "[CohesiveFaultKernel] WARNING: registerWithDS is deprecated; "
                 "use registerHybridWithDS (Session 4 hybrid-by-key path).\n");
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode CohesiveFaultKernel::registerHybridWithDS(PetscDS prob,
                                                          PetscWeakForm wf,
                                                          DMLabel material_label,
                                                          DMLabel cohesive_label,
                                                          PetscInt displacement_field,
                                                          PetscInt lagrange_field) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

    // -------------------------------------------------------------------------
    // Propagate current mode and friction coefficient into the PetscDS constant
    // array, matching the pre-Session 4 registerWithDS behaviour.
    // -------------------------------------------------------------------------
    PetscInt nconst_old = 0;
    const PetscScalar *old_c = nullptr;
    ierr = PetscDSGetConstants(prob, &nconst_old, &old_c); CHKERRQ(ierr);

    PetscInt nconst_new = (nconst_old >= COHESIVE_CONST_COUNT)
                          ? nconst_old : COHESIVE_CONST_COUNT;
    std::vector<PetscScalar> new_c(static_cast<std::size_t>(nconst_new), 0.0);
    for (PetscInt i = 0; i < nconst_old; ++i) new_c[i] = old_c[i];
    // Mode slot: 0 locked, 1 slipping, 2 prescribed_slip.
    double mode_val = locked_ ? 0.0 : 1.0;
    if (prescribed_slip_mode_) mode_val = 2.0;
    new_c[COHESIVE_CONST_MODE] = mode_val;
    new_c[COHESIVE_CONST_MU_F] = static_cast<PetscScalar>(mu_friction_);
    new_c[COHESIVE_CONST_TENSILE_STRENGTH] = static_cast<PetscScalar>(tensile_strength_);
    if (prescribed_slip_mode_) {
        new_c[COHESIVE_CONST_PRESCRIBED_SLIP_X] =
            static_cast<PetscScalar>(prescribed_slip_[0]);
        new_c[COHESIVE_CONST_PRESCRIBED_SLIP_Y] =
            static_cast<PetscScalar>(prescribed_slip_[1]);
        new_c[COHESIVE_CONST_PRESCRIBED_SLIP_Z] =
            static_cast<PetscScalar>(prescribed_slip_[2]);
    }
    ierr = PetscDSSetConstants(prob, nconst_new, new_c.data()); CHKERRQ(ierr);

    // -------------------------------------------------------------------------
    // Register the three weak forms against explicit (label, value, field)
    // keys. Matches PETSc ex5.c:TestAssembly lines 1194-1202.
    // -------------------------------------------------------------------------
    ierr = PetscWeakFormSetIndexBdResidual(
        wf, material_label, 1, displacement_field, 0,
        0, f0_hybrid_u_neg, 0, nullptr); CHKERRQ(ierr);
    ierr = PetscWeakFormSetIndexBdResidual(
        wf, material_label, 2, displacement_field, 0,
        0, f0_hybrid_u_pos, 0, nullptr); CHKERRQ(ierr);
    ierr = PetscWeakFormSetIndexBdResidual(
        wf, cohesive_label, 1, lagrange_field, 0,
        0, f0_hybrid_lambda, 0, nullptr); CHKERRQ(ierr);

    // Jacobians:
    //   (material, 1, disp, lambda) -> g0 = -I  (neg side)
    //   (material, 2, disp, lambda) -> g0 = +I  (pos side)
    //   (cohesive, 1, lambda, disp) -> block layout [-I | +I] in the combined
    //                                   u_neg / u_pos trial space.
    //   (cohesive, 1, lambda, lambda) -> zero for locked / prescribed, non-zero
    //                                     for slipping (ported from the old
    //                                     g0_lagrange_lagrange).
    ierr = PetscWeakFormSetIndexBdJacobian(
        wf, material_label, 1, displacement_field, lagrange_field, 0,
        0, g0_hybrid_u_lambda_neg, 0, nullptr, 0, nullptr, 0, nullptr); CHKERRQ(ierr);
    ierr = PetscWeakFormSetIndexBdJacobian(
        wf, material_label, 2, displacement_field, lagrange_field, 0,
        0, g0_hybrid_u_lambda_pos, 0, nullptr, 0, nullptr, 0, nullptr); CHKERRQ(ierr);
    ierr = PetscWeakFormSetIndexBdJacobian(
        wf, cohesive_label, 1, lagrange_field, displacement_field, 0,
        0, g0_hybrid_lambda_u, 0, nullptr, 0, nullptr, 0, nullptr); CHKERRQ(ierr);
    ierr = PetscWeakFormSetIndexBdJacobian(
        wf, cohesive_label, 1, lagrange_field, lagrange_field, 0,
        0, g0_hybrid_lambda_lambda, 0, nullptr, 0, nullptr, 0, nullptr); CHKERRQ(ierr);

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode CohesiveFaultKernel::updateAuxiliaryState(DM dm,
                                                          FaultCohesiveDyn* fault) {
    PetscFunctionBeginUser;
    (void)dm;
    (void)fault;
    // In a full implementation, this would copy per-vertex friction state
    // (theta, slip_rate, cumulative_slip) into a PetscDS auxiliary field
    // so the pointwise callbacks can access them via the a[] parameter.
    //
    // For the initial implementation, the friction state is accessed
    // through the FaultCohesiveDyn object directly in the standalone
    // computeResidual/computeJacobian methods.
    PetscFunctionReturn(0);
}

// ============================================================================
// Static PetscDS callback implementations
// ============================================================================

void CohesiveFaultKernel::f0_displacement_cohesive(
    PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[],
    const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[],
    const PetscScalar a[], const PetscScalar a_t[],
    const PetscScalar a_x[],
    PetscReal t, const PetscReal x[],
    const PetscReal n[],
    PetscInt numConstants,
    const PetscScalar constants[],
    PetscScalar f[]) {

    (void)Nf; (void)NfAux; (void)uOff_x; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)x; (void)numConstants; (void)constants;

    // The Lagrange multiplier (traction) is applied as a force to
    // the + and - sides of the fault.
    // For the negative side: f_u = -lambda (traction acts opposite)
    // For the positive side: f_u = +lambda
    //
    // In PETSc's hybrid cell convention, this callback is called once
    // for the cohesive cell, and the result is distributed to both sides.
    // The Lagrange multiplier field values are in u[uOff[lagrange_field]...]
    //
    // For a single-component displacement field per direction:
    // uOff[0] = displacement, uOff[1] = lagrange
    // The n[] provides the fault normal direction.

    // Apply traction from Lagrange multiplier along the normal direction
    for (PetscInt d = 0; d < dim; ++d) {
        // Lagrange multiplier acts as traction: ±lambda on each side
        // The sign convention depends on which side this is being evaluated for.
        // For the integrated constraint, we contribute lambda_d to the displacement.
        f[d] = u[uOff[1] + d];
    }
}

void CohesiveFaultKernel::f0_lagrange_constraint(
    PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[],
    const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[],
    const PetscScalar a[], const PetscScalar a_t[],
    const PetscScalar a_x[],
    PetscReal t, const PetscReal x[],
    const PetscReal n[],
    PetscInt numConstants,
    const PetscScalar constants[],
    PetscScalar f[]) {

    (void)Nf; (void)NfAux; (void)uOff_x; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)x;

    // Determine mode from dedicated constant slot (0 = locked, 1 = slipping).
    // COHESIVE_CONST_MODE = 25 to avoid aliasing PoroelasticSolver constants.
    PetscScalar mode = (numConstants > COHESIVE_CONST_MODE)
                       ? constants[COHESIVE_CONST_MODE] : 0.0;

    if (PetscRealPart(mode) < 0.5) {
        // LOCKED: constraint = displacement jump = u+ - u- = 0
        // With optional tensile failure: if the Lagrange multiplier
        // (normal traction) exceeds the tensile strength, release
        // the normal constraint to allow opening.
        PetscScalar T_s = (numConstants > COHESIVE_CONST_TENSILE_STRENGTH)
                          ? constants[COHESIVE_CONST_TENSILE_STRENGTH] : 0.0;

        if (PetscRealPart(T_s) > 0.0) {
            // Check Lagrange multiplier for tensile failure
            PetscScalar lambda_n = 0.0;
            for (PetscInt d = 0; d < dim; ++d) {
                lambda_n += u[uOff[1] + d] * n[d];
            }

            if (PetscRealPart(lambda_n) > PetscRealPart(T_s)) {
                // Tensile failure: release constraint, enforce zero traction
                // f = lambda (so solver drives lambda -> 0)
                for (PetscInt d = 0; d < dim; ++d) {
                    f[d] = u[uOff[1] + d];
                }
            } else {
                // Still locked: zero displacement jump
                for (PetscInt d = 0; d < dim; ++d) {
                    f[d] = u[uOff[0] + dim + d] - u[uOff[0] + d];
                }
            }
        } else {
            // Standard locked: zero displacement jump (no tensile check)
            for (PetscInt d = 0; d < dim; ++d) {
                f[d] = u[uOff[0] + dim + d] - u[uOff[0] + d];
            }
        }
    } else {
        // SLIPPING: friction-governed constraint
        // f_lambda = lambda - mu_f * sigma_n * slip_direction
        // where lambda is the Lagrange multiplier (traction),
        // mu_f is the friction coefficient, sigma_n is the effective
        // normal stress, and slip_direction = slip_rate / |slip_rate|.
        //
        // For the PetscDS callback, we use a simplified form:
        // f = lambda_tangential - tau_strength * slip_hat
        // This requires friction state from auxiliary fields.

        // Compute displacement jump
        PetscScalar slip[3] = {0.0, 0.0, 0.0};
        for (PetscInt d = 0; d < dim; ++d) {
            slip[d] = u[uOff[0] + dim + d] - u[uOff[0] + d];
        }

        // Normal component of slip
        PetscScalar slip_n = 0.0;
        for (PetscInt d = 0; d < dim; ++d) {
            slip_n += slip[d] * n[d];
        }

        // Tangential slip
        PetscScalar slip_t[3];
        for (PetscInt d = 0; d < dim; ++d) {
            slip_t[d] = slip[d] - slip_n * n[d];
        }
        PetscScalar slip_mag = 0.0;
        for (PetscInt d = 0; d < dim; ++d) {
            slip_mag += PetscRealPart(slip_t[d]) * PetscRealPart(slip_t[d]);
        }
        slip_mag = std::sqrt(PetscRealPart(slip_mag));

        // Lagrange multiplier (traction)
        PetscScalar lambda_n = 0.0;
        for (PetscInt d = 0; d < dim; ++d) {
            lambda_n += u[uOff[1] + d] * n[d];
        }

        // Friction coefficient: constant or slip-weakening
        PetscReal friction_model = (numConstants > COHESIVE_CONST_FRICTION_MODEL)
                                   ? PetscRealPart(constants[COHESIVE_CONST_FRICTION_MODEL]) : 0.0;
        PetscReal mu_f_val;
        if (friction_model > 0.5) {
            // Linear slip-weakening: mu = mu_s - (mu_s - mu_d) * min(|slip_t|, Dc) / Dc
            PetscReal mu_s = (numConstants > COHESIVE_CONST_MU_S)
                             ? PetscRealPart(constants[COHESIVE_CONST_MU_S]) : 0.677;
            PetscReal mu_d = (numConstants > COHESIVE_CONST_MU_D)
                             ? PetscRealPart(constants[COHESIVE_CONST_MU_D]) : 0.525;
            PetscReal D_c = (numConstants > COHESIVE_CONST_DC)
                            ? PetscRealPart(constants[COHESIVE_CONST_DC]) : 0.40;
            PetscReal slip_mag_real = PetscRealPart(slip_mag);
            mu_f_val = mu_s - (mu_s - mu_d) * PetscMin(slip_mag_real, D_c) / D_c;
        } else {
            mu_f_val = (numConstants > COHESIVE_CONST_MU_F)
                       ? PetscRealPart(constants[COHESIVE_CONST_MU_F]) : 0.6;
        }

        // Friction strength: tau_f = mu_f * |sigma_n_eff|
        // sigma_n_eff is the normal traction (compressive positive)
        PetscScalar tau_f = mu_f_val * std::abs(PetscRealPart(-lambda_n));

        // Tangential traction from Lagrange multiplier
        PetscScalar lambda_t[3];
        for (PetscInt d = 0; d < dim; ++d) {
            lambda_t[d] = u[uOff[1] + d] - lambda_n * n[d];
        }
        PetscScalar lambda_t_mag = 0.0;
        for (PetscInt d = 0; d < dim; ++d) {
            lambda_t_mag += PetscRealPart(lambda_t[d]) * PetscRealPart(lambda_t[d]);
        }
        lambda_t_mag = std::sqrt(PetscRealPart(lambda_t_mag));

        // Constraint: tangential traction magnitude = friction strength
        // Normal: non-penetration (slip_n >= 0)
        for (PetscInt d = 0; d < dim; ++d) {
            if (slip_mag > 1e-15) {
                // Friction constraint in tangential direction
                PetscScalar slip_hat = slip_t[d] / slip_mag;
                f[d] = lambda_t[d] - tau_f * slip_hat;
            } else {
                // No slip: traction must be within friction cone
                f[d] = lambda_t[d];
                if (lambda_t_mag > PetscRealPart(tau_f) + 1e-10) {
                    f[d] = lambda_t[d] * (1.0 - tau_f / lambda_t_mag);
                } else {
                    f[d] = 0.0;  // Within friction cone
                }
            }
            // Normal direction: enforce non-penetration
            f[d] += (slip_n < 0.0 ? -slip_n : 0.0) * n[d];
        }
    }
}

void CohesiveFaultKernel::g0_displacement_lagrange(
    PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[],
    const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[],
    const PetscScalar a[], const PetscScalar a_t[],
    const PetscScalar a_x[],
    PetscReal t, PetscReal u_tShift,
    const PetscReal x[], const PetscReal n[],
    PetscInt numConstants,
    const PetscScalar constants[],
    PetscScalar g0[]) {

    (void)Nf; (void)NfAux; (void)uOff; (void)uOff_x;
    (void)u; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)u_tShift; (void)x; (void)n;
    (void)numConstants; (void)constants;

    // d(f_u)/d(lambda) = I (identity)
    // The traction from the Lagrange multiplier directly enters the
    // displacement residual, so the Jacobian is the identity matrix.
    for (PetscInt i = 0; i < dim; ++i) {
        for (PetscInt j = 0; j < dim; ++j) {
            g0[i * dim + j] = (i == j) ? 1.0 : 0.0;
        }
    }
}

void CohesiveFaultKernel::g0_lagrange_displacement(
    PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[],
    const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[],
    const PetscScalar a[], const PetscScalar a_t[],
    const PetscScalar a_x[],
    PetscReal t, PetscReal u_tShift,
    const PetscReal x[], const PetscReal n[],
    PetscInt numConstants,
    const PetscScalar constants[],
    PetscScalar g0[]) {

    (void)Nf; (void)NfAux; (void)uOff_x;
    (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)u_tShift; (void)x;

    PetscScalar mode = (numConstants > COHESIVE_CONST_MODE)
                       ? constants[COHESIVE_CONST_MODE] : 0.0;

    if (PetscRealPart(mode) < 0.5) {
        // LOCKED: d(f_lambda)/d(u) = I
        // The locked constraint is f_lambda = u+ - u-
        // PETSc hybrid cell convention handles the +/- sign.
        for (PetscInt i = 0; i < dim; ++i) {
            for (PetscInt j = 0; j < dim; ++j) {
                g0[i * dim + j] = (i == j) ? 1.0 : 0.0;
            }
        }
    } else {
        // SLIPPING: d(C_d)/d(u_e) where the constraint is:
        //   C_d = lambda_t_d - tau_f * s_hat_d + max(-slip_n, 0) * n_d
        //
        // d(C_d)/d(u_e) = -tau_f * d(s_hat_d)/d(u_e) - H(-slip_n) * n_d * n_e
        // where d(s_hat_d)/d(u_e) = (P_de - s_hat_d * s_hat_e) / |slip_t|
        // and P_de = delta_de - n_d * n_e (tangential projector)

        // Compute displacement jump (slip)
        PetscScalar slip[3] = {0.0, 0.0, 0.0};
        for (PetscInt d = 0; d < dim; ++d) {
            slip[d] = u[uOff[0] + dim + d] - u[uOff[0] + d];
        }

        // Normal component of slip
        PetscReal slip_n = 0.0;
        for (PetscInt d = 0; d < dim; ++d) {
            slip_n += PetscRealPart(slip[d]) * n[d];
        }

        // Tangential slip
        PetscReal slip_t[3] = {0.0, 0.0, 0.0};
        PetscReal slip_t_mag2 = 0.0;
        for (PetscInt d = 0; d < dim; ++d) {
            slip_t[d] = PetscRealPart(slip[d]) - slip_n * n[d];
            slip_t_mag2 += slip_t[d] * slip_t[d];
        }
        PetscReal slip_t_mag = std::sqrt(slip_t_mag2);

        // Slip direction
        PetscReal s_hat[3] = {0.0, 0.0, 0.0};
        if (slip_t_mag > 1e-15) {
            for (PetscInt d = 0; d < dim; ++d) {
                s_hat[d] = slip_t[d] / slip_t_mag;
            }
        }

        // Friction coefficient: constant or slip-weakening
        PetscReal friction_model = (numConstants > COHESIVE_CONST_FRICTION_MODEL)
                                   ? PetscRealPart(constants[COHESIVE_CONST_FRICTION_MODEL]) : 0.0;
        PetscReal mu_s = 0.0, mu_d = 0.0, D_c = 1.0;
        PetscReal mu_f;
        if (friction_model > 0.5) {
            mu_s = (numConstants > COHESIVE_CONST_MU_S)
                   ? PetscRealPart(constants[COHESIVE_CONST_MU_S]) : 0.677;
            mu_d = (numConstants > COHESIVE_CONST_MU_D)
                   ? PetscRealPart(constants[COHESIVE_CONST_MU_D]) : 0.525;
            D_c = (numConstants > COHESIVE_CONST_DC)
                  ? PetscRealPart(constants[COHESIVE_CONST_DC]) : 0.40;
            mu_f = mu_s - (mu_s - mu_d) * std::min(slip_t_mag, D_c) / D_c;
        } else {
            mu_f = (numConstants > COHESIVE_CONST_MU_F)
                   ? PetscRealPart(constants[COHESIVE_CONST_MU_F]) : 0.6;
        }

        PetscReal lambda_n = 0.0;
        for (PetscInt d = 0; d < dim; ++d) {
            lambda_n += PetscRealPart(u[uOff[1] + d]) * n[d];
        }
        PetscReal sigma_n_comp = std::max(0.0, -lambda_n);
        PetscReal tau_f = mu_f * sigma_n_comp;

        // Slip-weakening derivative: d(mu_f)/d(|slip_t|) = -(mu_s - mu_d)/Dc
        // when |slip_t| < Dc, 0 otherwise
        PetscReal dmu_dslip = 0.0;
        if (friction_model > 0.5 && slip_t_mag < D_c) {
            dmu_dslip = -(mu_s - mu_d) / D_c;
        }

        for (PetscInt i = 0; i < dim; ++i) {
            for (PetscInt j = 0; j < dim; ++j) {
                PetscReal P_ij = ((i == j) ? 1.0 : 0.0) - n[i] * n[j];
                PetscReal dS_ij = 0.0;
                if (slip_t_mag > 1e-15) {
                    dS_ij = (P_ij - s_hat[i] * s_hat[j]) / slip_t_mag;
                }
                g0[i * dim + j] = -tau_f * dS_ij;

                // Slip-weakening contribution: d(tau_f * s_hat_i)/d(u_j)
                // includes d(mu_f)/d(|slip_t|) * d(|slip_t|)/d(u_j) * sigma_n * s_hat_i
                // d(|slip_t|)/d(u_j) = s_hat_j, so:
                // contribution = -dmu_dslip * sigma_n_comp * s_hat_i * s_hat_j
                if (dmu_dslip != 0.0 && slip_t_mag > 1e-15) {
                    g0[i * dim + j] += -dmu_dslip * sigma_n_comp * s_hat[i] * s_hat[j];
                }

                // Normal non-penetration: d(max(-slip_n,0)*n_d)/d(u_e)
                // = -H(-slip_n) * n_d * n_e (derivative of max(-slip_n,0) w.r.t. slip_e)
                if (slip_n < 0.0) {
                    g0[i * dim + j] += -n[i] * n[j];
                }
            }
        }
    }
}

void CohesiveFaultKernel::g0_lagrange_lagrange(
    PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[],
    const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[],
    const PetscScalar a[], const PetscScalar a_t[],
    const PetscScalar a_x[],
    PetscReal t, PetscReal u_tShift,
    const PetscReal x[], const PetscReal n[],
    PetscInt numConstants,
    const PetscScalar constants[],
    PetscScalar g0[]) {

    (void)Nf; (void)NfAux; (void)uOff_x;
    (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)u_tShift; (void)x;

    PetscScalar mode = (numConstants > COHESIVE_CONST_MODE)
                       ? constants[COHESIVE_CONST_MODE] : 0.0;

    if (PetscRealPart(mode) < 0.5) {
        // LOCKED: d(f_lambda)/d(lambda) = 0
        // The constraint (u+ - u-) does not depend on lambda.
        for (PetscInt i = 0; i < dim; ++i) {
            for (PetscInt j = 0; j < dim; ++j) {
                g0[i * dim + j] = 0.0;
            }
        }
    } else {
        // SLIPPING: d(C_d)/d(lambda_e)
        // The constraint is:
        //   C_d = lambda_t_d - tau_f * s_hat_d + max(-slip_n, 0) * n_d
        //
        // d(C_d)/d(lambda_e):
        //   d(lambda_t_d)/d(lambda_e) = P_de = delta_de - n_d * n_e
        //   d(-tau_f * s_hat_d)/d(lambda_e) = -mu_f * sign(-lambda_n) * s_hat_d * (-n_e)
        //     since tau_f = mu_f * |compressive sigma_n| = mu_f * max(0, -lambda_n)
        //     d(tau_f)/d(lambda_e) = mu_f * sign(-lambda_n) * (-n_e)
        //                          = -mu_f * H(-lambda_n) * n_e  (where H is Heaviside, lambda_n = lam . n)
        //   d(max(-slip_n,0)*n_d)/d(lambda_e) = 0 (slip does not depend on lambda)
        //
        // So: d(C_d)/d(lambda_e) = P_de + mu_f * H(sigma_n_comp > 0) * s_hat_d * n_e

        // Compute slip for the slip direction
        PetscScalar slip[3] = {0.0, 0.0, 0.0};
        for (PetscInt d = 0; d < dim; ++d) {
            slip[d] = u[uOff[0] + dim + d] - u[uOff[0] + d];
        }

        // Normal component of slip
        PetscReal slip_n_val = 0.0;
        for (PetscInt d = 0; d < dim; ++d) {
            slip_n_val += PetscRealPart(slip[d]) * n[d];
        }

        // Tangential slip and direction
        PetscReal slip_t[3] = {0.0, 0.0, 0.0};
        PetscReal slip_t_mag2 = 0.0;
        for (PetscInt d = 0; d < dim; ++d) {
            slip_t[d] = PetscRealPart(slip[d]) - slip_n_val * n[d];
            slip_t_mag2 += slip_t[d] * slip_t[d];
        }
        PetscReal slip_t_mag = std::sqrt(slip_t_mag2);

        PetscReal s_hat[3] = {0.0, 0.0, 0.0};
        if (slip_t_mag > 1e-15) {
            for (PetscInt d = 0; d < dim; ++d) {
                s_hat[d] = slip_t[d] / slip_t_mag;
            }
        }

        // Friction coefficient: constant or slip-weakening
        PetscReal friction_model = (numConstants > COHESIVE_CONST_FRICTION_MODEL)
                                   ? PetscRealPart(constants[COHESIVE_CONST_FRICTION_MODEL]) : 0.0;
        PetscReal mu_f;
        if (friction_model > 0.5) {
            PetscReal mu_s = (numConstants > COHESIVE_CONST_MU_S)
                             ? PetscRealPart(constants[COHESIVE_CONST_MU_S]) : 0.677;
            PetscReal mu_d = (numConstants > COHESIVE_CONST_MU_D)
                             ? PetscRealPart(constants[COHESIVE_CONST_MU_D]) : 0.525;
            PetscReal D_c = (numConstants > COHESIVE_CONST_DC)
                            ? PetscRealPart(constants[COHESIVE_CONST_DC]) : 0.40;
            mu_f = mu_s - (mu_s - mu_d) * std::min(slip_t_mag, D_c) / D_c;
        } else {
            mu_f = (numConstants > COHESIVE_CONST_MU_F)
                   ? PetscRealPart(constants[COHESIVE_CONST_MU_F]) : 0.6;
        }

        PetscReal lambda_n = 0.0;
        for (PetscInt d = 0; d < dim; ++d) {
            lambda_n += PetscRealPart(u[uOff[1] + d]) * n[d];
        }
        PetscReal sigma_n_comp = std::max(0.0, -lambda_n);

        for (PetscInt i = 0; i < dim; ++i) {
            for (PetscInt j = 0; j < dim; ++j) {
                // Tangential projector: d(lambda_t_d)/d(lambda_e)
                PetscReal P_ij = ((i == j) ? 1.0 : 0.0) - n[i] * n[j];
                g0[i * dim + j] = P_ij;

                // Friction strength derivative: d(-tau_f * s_hat_d)/d(lambda_e)
                // = mu_f * s_hat_d * n_e (when under compression and slipping)
                if (sigma_n_comp > 0.0 && slip_t_mag > 1e-15) {
                    g0[i * dim + j] += mu_f * s_hat[i] * n[j];
                }
            }
        }
    }
}

PetscErrorCode CohesiveFaultKernel::registerPrescribedSlipWithDS(
    PetscDS prob, int displacement_field, int lagrange_field) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

    // Read existing constants, expand to COHESIVE_CONST_COUNT, patch prescribed slip slots
    PetscInt nconst_old = 0;
    const PetscScalar *old_c = nullptr;
    ierr = PetscDSGetConstants(prob, &nconst_old, &old_c); CHKERRQ(ierr);

    PetscInt nconst_new = (nconst_old >= COHESIVE_CONST_COUNT)
                          ? nconst_old : COHESIVE_CONST_COUNT;
    std::vector<PetscScalar> new_c(static_cast<std::size_t>(nconst_new), 0.0);
    for (PetscInt i = 0; i < nconst_old; ++i) new_c[i] = old_c[i];
    new_c[COHESIVE_CONST_MODE] = 2.0;  // prescribed_slip mode
    new_c[COHESIVE_CONST_MU_F] = static_cast<PetscScalar>(mu_friction_);
    new_c[COHESIVE_CONST_TENSILE_STRENGTH] = static_cast<PetscScalar>(tensile_strength_);
    new_c[COHESIVE_CONST_PRESCRIBED_SLIP_X] = static_cast<PetscScalar>(prescribed_slip_[0]);
    new_c[COHESIVE_CONST_PRESCRIBED_SLIP_Y] = static_cast<PetscScalar>(prescribed_slip_[1]);
    new_c[COHESIVE_CONST_PRESCRIBED_SLIP_Z] = static_cast<PetscScalar>(prescribed_slip_[2]);
    ierr = PetscDSSetConstants(prob, nconst_new, new_c.data()); CHKERRQ(ierr);

    // Register residual: displacement still gets traction from lambda
    ierr = PetscDSSetBdResidual(prob, displacement_field,
                                 f0_displacement_cohesive, nullptr); CHKERRQ(ierr);
    // Lagrange multiplier gets prescribed slip constraint
    ierr = PetscDSSetBdResidual(prob, lagrange_field,
                                 f0_prescribed_slip, nullptr); CHKERRQ(ierr);

    // Jacobian blocks are identical to the locked case:
    // d(f_u)/d(lambda) = I, d(f_lambda)/d(u) = I, d(f_lambda)/d(lambda) = 0
    ierr = PetscDSSetBdJacobian(prob, displacement_field, lagrange_field,
                                 g0_displacement_lagrange, nullptr,
                                 nullptr, nullptr); CHKERRQ(ierr);
    ierr = PetscDSSetBdJacobian(prob, lagrange_field, displacement_field,
                                 g0_lagrange_displacement, nullptr,
                                 nullptr, nullptr); CHKERRQ(ierr);
    ierr = PetscDSSetBdJacobian(prob, lagrange_field, lagrange_field,
                                 g0_lagrange_lagrange, nullptr,
                                 nullptr, nullptr); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

// ============================================================================
// Prescribed slip callback
// ============================================================================

void CohesiveFaultKernel::f0_prescribed_slip(
    PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[],
    const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[],
    const PetscScalar a[], const PetscScalar a_t[],
    const PetscScalar a_x[],
    PetscReal t, const PetscReal x[],
    const PetscReal n[],
    PetscInt numConstants,
    const PetscScalar constants[],
    PetscScalar f[]) {

    (void)Nf; (void)NfAux; (void)uOff_x; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)x; (void)n;

    // Prescribed slip constraint: f_lambda = (u+ - u-) - delta_prescribed
    // u[uOff[0]..uOff[0]+dim-1] = u_minus (negative side)
    // u[uOff[0]+dim..uOff[0]+2*dim-1] = u_plus (positive side)
    // constants[COHESIVE_CONST_PRESCRIBED_SLIP_X..Z] = delta
    for (PetscInt d = 0; d < dim; ++d) {
        PetscScalar jump = u[uOff[0] + dim + d] - u[uOff[0] + d];
        PetscScalar delta = (numConstants > COHESIVE_CONST_PRESCRIBED_SLIP_X + d)
                            ? constants[COHESIVE_CONST_PRESCRIBED_SLIP_X + d] : 0.0;
        f[d] = jump - delta;
    }
}

// ============================================================================
// PETSc hybrid-driver callbacks (Session 4)
//
// Layout conventions (see PETSc src/dm/impls/plex/tests/ex5.c TestAssembly):
//   - For the displacement weak form (key material/1 or 2):
//       u[uOff[0] + c] = displacement component c on that side
//       u[uOff[1] + c] = Lagrange multiplier component c
//   - For the Lagrange weak form (key cohesive/1):
//       Nc          = uOff[2] - uOff[1] = dim
//       u[c]        = negative-side displacement, c in [0, Nc)
//       u[Nc + c]   = positive-side displacement, c in [0, Nc)
//       u[uOff[1] + c] = Lagrange multiplier
//
// Residual sign convention: we enforce (u_pos - u_neg) = delta_prescribed,
// so f0_hybrid_lambda[c] = u[Nc+c] - u[c] - delta[c] for the prescribed /
// locked (delta=0) branches. The Jacobian dR_lambda/du is then [-I | +I].
// The displacement callbacks match PETSc ex5.c exactly: -lambda on the
// negative side, +lambda on the positive side.
// ============================================================================

void CohesiveFaultKernel::f0_hybrid_u_neg(
    PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[],
    const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[],
    const PetscScalar a[], const PetscScalar a_t[],
    const PetscScalar a_x[],
    PetscReal t, const PetscReal x[],
    const PetscReal n[],
    PetscInt numConstants,
    const PetscScalar constants[],
    PetscScalar f0[]) {

    (void)Nf; (void)NfAux; (void)uOff_x; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)x; (void)n; (void)numConstants; (void)constants;

    for (PetscInt c = 0; c < dim; ++c) {
        f0[c] = -u[uOff[1] + c];
    }
}

void CohesiveFaultKernel::f0_hybrid_u_pos(
    PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[],
    const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[],
    const PetscScalar a[], const PetscScalar a_t[],
    const PetscScalar a_x[],
    PetscReal t, const PetscReal x[],
    const PetscReal n[],
    PetscInt numConstants,
    const PetscScalar constants[],
    PetscScalar f0[]) {

    (void)Nf; (void)NfAux; (void)uOff_x; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)x; (void)n; (void)numConstants; (void)constants;

    for (PetscInt c = 0; c < dim; ++c) {
        f0[c] = u[uOff[1] + c];
    }
}

void CohesiveFaultKernel::f0_hybrid_lambda(
    PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[],
    const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[],
    const PetscScalar a[], const PetscScalar a_t[],
    const PetscScalar a_x[],
    PetscReal t, const PetscReal x[],
    const PetscReal n[],
    PetscInt numConstants,
    const PetscScalar constants[],
    PetscScalar f0[]) {

    (void)Nf; (void)uOff_x; (void)u_t; (void)u_x;
    (void)aOff_x; (void)a_t; (void)a_x;
    (void)t; (void)x;

    const PetscInt Nc = uOff[2] - uOff[1];

    // Decode mode. Slot layout matches the PetscDS constants array populated
    // in Simulator::setupPhysics.
    PetscReal mode = (numConstants > COHESIVE_CONST_MODE)
                     ? PetscRealPart(constants[COHESIVE_CONST_MODE]) : 0.0;

    if (mode < 0.5) {
        // LOCKED: jump = 0.
        for (PetscInt c = 0; c < Nc; ++c) {
            f0[c] = u[Nc + c] - u[c];
        }
    } else if (mode > 1.5) {
        // PRESCRIBED SLIP. Session 16 attempted the PyLith f0l_slip
        // pattern (read slip from auxiliary field a[aOff[0]..aOff[0]+2]
        // in fault-local coordinates and rotate via refDir1/refDir2 to
        // global XYZ) but PETSc 3.25's hybrid driver (plexfem.c:5644 ff.)
        // enforces closure-vs-DS totDim / tabulation point consistency
        // between the main Lagrange field and any cohesive-key aux field
        // that this session could not satisfy with a single shared aux
        // DM. Cell-constant degree-0 slip avoids the totDim check but
        // triggers the tabulation-point-count mismatch in
        // PetscFEIntegrateHybridResidual_Basic; label-restricted surface
        // FE avoids the tabulation mismatch but the section closure on
        // cohesive prism cells does not match the FE totDim.
        //
        // The prescribed-slip branch therefore retains the constants
        // array path for Session 16. The aux-field plumbing
        // (faultAuxDM_/faultAuxVec_/updateFaultAuxiliary) stays in tree
        // as groundwork for a future session that reconstructs the aux
        // as a PyLith-style shared Field on the main DM instead of a
        // clone. refDir1/refDir2 slots (80..85) are populated by
        // setupPhysics so the rotation is ready when the aux path lands.
        for (PetscInt c = 0; c < Nc; ++c) {
            PetscScalar delta =
                (numConstants > COHESIVE_CONST_PRESCRIBED_SLIP_X + c)
                ? constants[COHESIVE_CONST_PRESCRIBED_SLIP_X + c] : 0.0;
            f0[c] = u[Nc + c] - u[c] - delta;
        }
    } else {
        // SLIPPING (Coulomb friction, semi-smooth Newton). Ported from the
        // pre-Session-4 f0_lagrange_constraint, with slip = u_pos - u_neg
        // resolved via the hybrid (u[c], u[Nc+c]) layout.
        PetscScalar slip[3] = {0.0, 0.0, 0.0};
        for (PetscInt d = 0; d < Nc; ++d) {
            slip[d] = u[Nc + d] - u[d];
        }

        PetscReal slip_n = 0.0;
        for (PetscInt d = 0; d < Nc; ++d) {
            slip_n += PetscRealPart(slip[d]) * n[d];
        }

        PetscReal slip_t[3] = {0.0, 0.0, 0.0};
        PetscReal slip_t_mag2 = 0.0;
        for (PetscInt d = 0; d < Nc; ++d) {
            slip_t[d] = PetscRealPart(slip[d]) - slip_n * n[d];
            slip_t_mag2 += slip_t[d] * slip_t[d];
        }
        PetscReal slip_mag = std::sqrt(slip_t_mag2);

        PetscReal lambda_n = 0.0;
        for (PetscInt d = 0; d < Nc; ++d) {
            lambda_n += PetscRealPart(u[uOff[1] + d]) * n[d];
        }

        PetscReal friction_model =
            (numConstants > COHESIVE_CONST_FRICTION_MODEL)
            ? PetscRealPart(constants[COHESIVE_CONST_FRICTION_MODEL]) : 0.0;
        PetscReal mu_f_val;
        if (friction_model > 0.5) {
            PetscReal mu_s = (numConstants > COHESIVE_CONST_MU_S)
                             ? PetscRealPart(constants[COHESIVE_CONST_MU_S]) : 0.677;
            PetscReal mu_d = (numConstants > COHESIVE_CONST_MU_D)
                             ? PetscRealPart(constants[COHESIVE_CONST_MU_D]) : 0.525;
            PetscReal D_c = (numConstants > COHESIVE_CONST_DC)
                            ? PetscRealPart(constants[COHESIVE_CONST_DC]) : 0.40;
            mu_f_val = mu_s - (mu_s - mu_d) * PetscMin(slip_mag, D_c) / D_c;
        } else {
            mu_f_val = (numConstants > COHESIVE_CONST_MU_F)
                       ? PetscRealPart(constants[COHESIVE_CONST_MU_F]) : 0.6;
        }

        PetscReal tau_f = mu_f_val * std::abs(-lambda_n);

        PetscReal lambda_t[3] = {0.0, 0.0, 0.0};
        PetscReal lambda_t_mag2 = 0.0;
        for (PetscInt d = 0; d < Nc; ++d) {
            lambda_t[d] = PetscRealPart(u[uOff[1] + d]) - lambda_n * n[d];
            lambda_t_mag2 += lambda_t[d] * lambda_t[d];
        }
        PetscReal lambda_t_mag = std::sqrt(lambda_t_mag2);

        for (PetscInt d = 0; d < Nc; ++d) {
            PetscScalar component;
            if (slip_mag > 1e-15) {
                PetscReal slip_hat = slip_t[d] / slip_mag;
                component = lambda_t[d] - tau_f * slip_hat;
            } else {
                if (lambda_t_mag > tau_f + 1e-10) {
                    component = lambda_t[d] * (1.0 - tau_f / lambda_t_mag);
                } else {
                    component = 0.0;
                }
            }
            component += (slip_n < 0.0 ? -slip_n : 0.0) * n[d];
            f0[d] = component;
        }
    }
}

void CohesiveFaultKernel::g0_hybrid_u_lambda_neg(
    PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[],
    const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[],
    const PetscScalar a[], const PetscScalar a_t[],
    const PetscScalar a_x[],
    PetscReal t, PetscReal u_tShift,
    const PetscReal x[], const PetscReal n[],
    PetscInt numConstants,
    const PetscScalar constants[],
    PetscScalar g0[]) {

    (void)Nf; (void)NfAux; (void)uOff; (void)uOff_x;
    (void)u; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)u_tShift; (void)x; (void)n;
    (void)numConstants; (void)constants;

    for (PetscInt c = 0; c < dim; ++c) {
        g0[c * dim + c] = -1.0;
    }
}

void CohesiveFaultKernel::g0_hybrid_u_lambda_pos(
    PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[],
    const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[],
    const PetscScalar a[], const PetscScalar a_t[],
    const PetscScalar a_x[],
    PetscReal t, PetscReal u_tShift,
    const PetscReal x[], const PetscReal n[],
    PetscInt numConstants,
    const PetscScalar constants[],
    PetscScalar g0[]) {

    (void)Nf; (void)NfAux; (void)uOff; (void)uOff_x;
    (void)u; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)u_tShift; (void)x; (void)n;
    (void)numConstants; (void)constants;

    for (PetscInt c = 0; c < dim; ++c) {
        g0[c * dim + c] = 1.0;
    }
}

void CohesiveFaultKernel::g0_hybrid_lambda_u(
    PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[],
    const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[],
    const PetscScalar a[], const PetscScalar a_t[],
    const PetscScalar a_x[],
    PetscReal t, PetscReal u_tShift,
    const PetscReal x[], const PetscReal n[],
    PetscInt numConstants,
    const PetscScalar constants[],
    PetscScalar g0[]) {

    (void)Nf; (void)NfAux; (void)uOff_x;
    (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)u_tShift; (void)x;

    const PetscInt Nc = uOff[2] - uOff[1];

    PetscReal mode = (numConstants > COHESIVE_CONST_MODE)
                     ? PetscRealPart(constants[COHESIVE_CONST_MODE]) : 0.0;

    if (mode < 0.5 || mode > 1.5) {
        // LOCKED or PRESCRIBED: d(u_pos - u_neg - delta)/d(u_neg) = -I,
        //                        d(...)/d(u_pos) = +I.
        for (PetscInt c = 0; c < Nc; ++c) {
            g0[c * Nc + c]            = -1.0;
            g0[Nc * Nc + c * Nc + c]  = 1.0;
        }
    } else {
        // SLIPPING: port the semi-smooth Newton tangent from the pre-Session-4
        // g0_lagrange_displacement. d(slip)/d(u_neg) = -I, d(slip)/d(u_pos) =
        // +I, so the neg-side block is -G and the pos-side block is +G, where
        // G = d(constraint)/d(slip).
        PetscScalar slip[3] = {0.0, 0.0, 0.0};
        for (PetscInt d = 0; d < Nc; ++d) {
            slip[d] = u[Nc + d] - u[d];
        }
        PetscReal slip_n = 0.0;
        for (PetscInt d = 0; d < Nc; ++d) {
            slip_n += PetscRealPart(slip[d]) * n[d];
        }
        PetscReal slip_t[3] = {0.0, 0.0, 0.0};
        PetscReal slip_t_mag2 = 0.0;
        for (PetscInt d = 0; d < Nc; ++d) {
            slip_t[d] = PetscRealPart(slip[d]) - slip_n * n[d];
            slip_t_mag2 += slip_t[d] * slip_t[d];
        }
        PetscReal slip_t_mag = std::sqrt(slip_t_mag2);

        PetscReal s_hat[3] = {0.0, 0.0, 0.0};
        if (slip_t_mag > 1e-15) {
            for (PetscInt d = 0; d < Nc; ++d) {
                s_hat[d] = slip_t[d] / slip_t_mag;
            }
        }

        PetscReal friction_model =
            (numConstants > COHESIVE_CONST_FRICTION_MODEL)
            ? PetscRealPart(constants[COHESIVE_CONST_FRICTION_MODEL]) : 0.0;
        PetscReal mu_s = 0.0, mu_d = 0.0, D_c = 1.0;
        PetscReal mu_f;
        if (friction_model > 0.5) {
            mu_s = (numConstants > COHESIVE_CONST_MU_S)
                   ? PetscRealPart(constants[COHESIVE_CONST_MU_S]) : 0.677;
            mu_d = (numConstants > COHESIVE_CONST_MU_D)
                   ? PetscRealPart(constants[COHESIVE_CONST_MU_D]) : 0.525;
            D_c = (numConstants > COHESIVE_CONST_DC)
                  ? PetscRealPart(constants[COHESIVE_CONST_DC]) : 0.40;
            mu_f = mu_s - (mu_s - mu_d) * std::min(slip_t_mag, D_c) / D_c;
        } else {
            mu_f = (numConstants > COHESIVE_CONST_MU_F)
                   ? PetscRealPart(constants[COHESIVE_CONST_MU_F]) : 0.6;
        }

        PetscReal lambda_n = 0.0;
        for (PetscInt d = 0; d < Nc; ++d) {
            lambda_n += PetscRealPart(u[uOff[1] + d]) * n[d];
        }
        PetscReal sigma_n_comp = std::max(0.0, -lambda_n);
        PetscReal tau_f = mu_f * sigma_n_comp;

        PetscReal dmu_dslip = 0.0;
        if (friction_model > 0.5 && slip_t_mag < D_c) {
            dmu_dslip = -(mu_s - mu_d) / D_c;
        }

        for (PetscInt i = 0; i < Nc; ++i) {
            for (PetscInt j = 0; j < Nc; ++j) {
                PetscReal P_ij = ((i == j) ? 1.0 : 0.0) - n[i] * n[j];
                PetscReal dS_ij = 0.0;
                if (slip_t_mag > 1e-15) {
                    dS_ij = (P_ij - s_hat[i] * s_hat[j]) / slip_t_mag;
                }
                PetscReal G_ij = -tau_f * dS_ij;
                if (dmu_dslip != 0.0 && slip_t_mag > 1e-15) {
                    G_ij += -dmu_dslip * sigma_n_comp * s_hat[i] * s_hat[j];
                }
                if (slip_n < 0.0) {
                    G_ij += -n[i] * n[j];
                }
                g0[i * Nc + j]            = -G_ij;
                g0[Nc * Nc + i * Nc + j]  =  G_ij;
            }
        }
    }
}

void CohesiveFaultKernel::g0_hybrid_lambda_lambda(
    PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[],
    const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[],
    const PetscScalar a[], const PetscScalar a_t[],
    const PetscScalar a_x[],
    PetscReal t, PetscReal u_tShift,
    const PetscReal x[], const PetscReal n[],
    PetscInt numConstants,
    const PetscScalar constants[],
    PetscScalar g0[]) {

    (void)Nf; (void)NfAux; (void)uOff_x;
    (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)u_tShift; (void)x;

    const PetscInt Nc = uOff[2] - uOff[1];
    (void)dim;

    PetscReal mode = (numConstants > COHESIVE_CONST_MODE)
                     ? PetscRealPart(constants[COHESIVE_CONST_MODE]) : 0.0;

    for (PetscInt i = 0; i < Nc; ++i) {
        for (PetscInt j = 0; j < Nc; ++j) {
            g0[i * Nc + j] = 0.0;
        }
    }

    if (mode < 0.5 || mode > 1.5) {
        // LOCKED or PRESCRIBED: no dependency on lambda.
        return;
    }

    // SLIPPING: ported from the pre-Session-4 g0_lagrange_lagrange.
    PetscScalar slip[3] = {0.0, 0.0, 0.0};
    for (PetscInt d = 0; d < Nc; ++d) {
        slip[d] = u[Nc + d] - u[d];
    }
    PetscReal slip_n_val = 0.0;
    for (PetscInt d = 0; d < Nc; ++d) {
        slip_n_val += PetscRealPart(slip[d]) * n[d];
    }
    PetscReal slip_t[3] = {0.0, 0.0, 0.0};
    PetscReal slip_t_mag2 = 0.0;
    for (PetscInt d = 0; d < Nc; ++d) {
        slip_t[d] = PetscRealPart(slip[d]) - slip_n_val * n[d];
        slip_t_mag2 += slip_t[d] * slip_t[d];
    }
    PetscReal slip_t_mag = std::sqrt(slip_t_mag2);

    PetscReal s_hat[3] = {0.0, 0.0, 0.0};
    if (slip_t_mag > 1e-15) {
        for (PetscInt d = 0; d < Nc; ++d) {
            s_hat[d] = slip_t[d] / slip_t_mag;
        }
    }

    PetscReal friction_model =
        (numConstants > COHESIVE_CONST_FRICTION_MODEL)
        ? PetscRealPart(constants[COHESIVE_CONST_FRICTION_MODEL]) : 0.0;
    PetscReal mu_f;
    if (friction_model > 0.5) {
        PetscReal mu_s = (numConstants > COHESIVE_CONST_MU_S)
                         ? PetscRealPart(constants[COHESIVE_CONST_MU_S]) : 0.677;
        PetscReal mu_d = (numConstants > COHESIVE_CONST_MU_D)
                         ? PetscRealPart(constants[COHESIVE_CONST_MU_D]) : 0.525;
        PetscReal D_c = (numConstants > COHESIVE_CONST_DC)
                        ? PetscRealPart(constants[COHESIVE_CONST_DC]) : 0.40;
        mu_f = mu_s - (mu_s - mu_d) * std::min(slip_t_mag, D_c) / D_c;
    } else {
        mu_f = (numConstants > COHESIVE_CONST_MU_F)
               ? PetscRealPart(constants[COHESIVE_CONST_MU_F]) : 0.6;
    }

    PetscReal lambda_n = 0.0;
    for (PetscInt d = 0; d < Nc; ++d) {
        lambda_n += PetscRealPart(u[uOff[1] + d]) * n[d];
    }
    PetscReal sigma_n_comp = std::max(0.0, -lambda_n);

    for (PetscInt i = 0; i < Nc; ++i) {
        for (PetscInt j = 0; j < Nc; ++j) {
            PetscReal P_ij = ((i == j) ? 1.0 : 0.0) - n[i] * n[j];
            g0[i * Nc + j] = P_ij;
            if (sigma_n_comp > 0.0 && slip_t_mag > 1e-15) {
                g0[i * Nc + j] += mu_f * s_hat[i] * n[j];
            }
        }
    }
}

} // namespace FSRM
