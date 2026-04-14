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

void CohesiveFaultKernel::setFrictionModel(FaultFrictionModel* model) {
    friction_model_ = model;
}

PetscErrorCode CohesiveFaultKernel::registerWithDS(PetscDS prob,
                                                     int displacement_field,
                                                     int lagrange_field) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

    // -------------------------------------------------------------------------
    // Propagate current mode and friction coefficient into the PetscDS constant
    // array.  We read whatever constants are already there, expand to at least
    // COHESIVE_CONST_COUNT slots, patch the two dedicated cohesive slots, and
    // write back.  This avoids overwriting PoroelasticSolver constants (0-23)
    // while ensuring the static callbacks always see the correct mode.
    // -------------------------------------------------------------------------
    PetscInt        nconst_old = 0;
    const PetscScalar *old_c   = nullptr;
    ierr = PetscDSGetConstants(prob, &nconst_old, &old_c); CHKERRQ(ierr);

    PetscInt nconst_new = (nconst_old >= COHESIVE_CONST_COUNT)
                          ? nconst_old : COHESIVE_CONST_COUNT;
    std::vector<PetscScalar> new_c(static_cast<std::size_t>(nconst_new), 0.0);
    for (PetscInt i = 0; i < nconst_old; ++i) new_c[i] = old_c[i];
    new_c[COHESIVE_CONST_MODE] = locked_ ? 0.0 : 1.0;
    new_c[COHESIVE_CONST_MU_F] = static_cast<PetscScalar>(mu_friction_);
    new_c[COHESIVE_CONST_TENSILE_STRENGTH] = static_cast<PetscScalar>(tensile_strength_);
    ierr = PetscDSSetConstants(prob, nconst_new, new_c.data()); CHKERRQ(ierr);

    // Register residual callbacks for cohesive cells
    // These will be called by PETSc's hybrid cell integration
    ierr = PetscDSSetBdResidual(prob, displacement_field,
                                 f0_displacement_cohesive, nullptr); CHKERRQ(ierr);
    ierr = PetscDSSetBdResidual(prob, lagrange_field,
                                 f0_lagrange_constraint, nullptr); CHKERRQ(ierr);

    // Register Jacobian callbacks for cohesive cells
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

        // Friction coefficient from dedicated constant slot (COHESIVE_CONST_MU_F = 25)
        PetscScalar mu_f = (numConstants > COHESIVE_CONST_MU_F)
                           ? constants[COHESIVE_CONST_MU_F] : 0.6;

        // Friction strength: tau_f = mu_f * |sigma_n_eff|
        // sigma_n_eff is the normal traction (compressive positive)
        PetscScalar tau_f = mu_f * std::abs(PetscRealPart(-lambda_n));

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

        // Friction coefficient and strength
        PetscReal mu_f = (numConstants > COHESIVE_CONST_MU_F)
                         ? PetscRealPart(constants[COHESIVE_CONST_MU_F]) : 0.6;
        PetscReal lambda_n = 0.0;
        for (PetscInt d = 0; d < dim; ++d) {
            lambda_n += PetscRealPart(u[uOff[1] + d]) * n[d];
        }
        PetscReal sigma_n_comp = std::max(0.0, -lambda_n);
        PetscReal tau_f = mu_f * sigma_n_comp;

        for (PetscInt i = 0; i < dim; ++i) {
            for (PetscInt j = 0; j < dim; ++j) {
                PetscReal P_ij = ((i == j) ? 1.0 : 0.0) - n[i] * n[j];
                PetscReal dS_ij = 0.0;
                if (slip_t_mag > 1e-15) {
                    dS_ij = (P_ij - s_hat[i] * s_hat[j]) / slip_t_mag;
                }
                g0[i * dim + j] = -tau_f * dS_ij;

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

        // Normal traction and friction
        PetscReal mu_f = (numConstants > COHESIVE_CONST_MU_F)
                         ? PetscRealPart(constants[COHESIVE_CONST_MU_F]) : 0.6;
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

} // namespace FSRM
