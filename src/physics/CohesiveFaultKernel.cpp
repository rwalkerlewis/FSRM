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

namespace FSRM {

CohesiveFaultKernel::CohesiveFaultKernel()
    : locked_(true), friction_model_(nullptr) {}

void CohesiveFaultKernel::setMode(bool is_locked) {
    locked_ = is_locked;
}

void CohesiveFaultKernel::setFrictionModel(FaultFrictionModel* model) {
    friction_model_ = model;
}

PetscErrorCode CohesiveFaultKernel::registerWithDS(PetscDS prob,
                                                     int displacement_field,
                                                     int lagrange_field) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

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

    // Determine mode from constants[0]: 0 = locked, 1 = slipping
    PetscScalar mode = (numConstants > 0) ? constants[0] : 0.0;

    if (PetscRealPart(mode) < 0.5) {
        // LOCKED: constraint = displacement jump = u+ - u- = 0
        // In PETSc's hybrid convention, the displacement values from both
        // sides are encoded in the u[] array. The exact layout depends on
        // the PETSc version. For the standard convention:
        // u[uOff[0]..] = displacement on negative side
        // u[uOff[0] + dim..] = displacement on positive side (shifted)
        // The normal n[] points from - to + side.
        //
        // For the locked constraint, we enforce: slip = u+ - u- = 0
        for (PetscInt d = 0; d < dim; ++d) {
            // Displacement jump projected onto fault coordinates
            // For simplicity, enforce zero jump in all components
            // A more sophisticated version would project onto
            // (normal, strike, dip) components
            f[d] = u[uOff[0] + dim + d] - u[uOff[0] + d];
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

        // Friction coefficient from constants
        PetscScalar mu_f = (numConstants > 1) ? constants[1] : 0.6;

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

    (void)Nf; (void)NfAux; (void)uOff; (void)uOff_x;
    (void)u; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)u_tShift; (void)x; (void)n;
    (void)numConstants; (void)constants;

    // For LOCKED: d(f_lambda)/d(u) = +I for positive side, -I for negative
    // The locked constraint is f_lambda = u+ - u-
    // d(f_lambda)/d(u+) = +I, d(f_lambda)/d(u-) = -I
    // In the combined form for the hybrid cell:
    for (PetscInt i = 0; i < dim; ++i) {
        for (PetscInt j = 0; j < dim; ++j) {
            g0[i * dim + j] = (i == j) ? 1.0 : 0.0;
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

    (void)Nf; (void)NfAux; (void)uOff; (void)uOff_x;
    (void)u; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)u_tShift; (void)x; (void)n;

    PetscScalar mode = (numConstants > 0) ? constants[0] : 0.0;

    // For LOCKED: d(f_lambda)/d(lambda) = 0 (constraint doesn't depend on lambda)
    // For SLIPPING: d(f_lambda)/d(lambda) = d(tau)/d(lambda) (traction direction)
    for (PetscInt i = 0; i < dim; ++i) {
        for (PetscInt j = 0; j < dim; ++j) {
            if (PetscRealPart(mode) < 0.5) {
                g0[i * dim + j] = 0.0;
            } else {
                // For slipping mode, the Jacobian depends on the friction
                // law and current state. Use identity as a preconditioner.
                g0[i * dim + j] = (i == j) ? 1.0 : 0.0;
            }
        }
    }
}

} // namespace FSRM
