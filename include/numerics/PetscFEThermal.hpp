#ifndef PETSCFE_THERMAL_HPP
#define PETSCFE_THERMAL_HPP

/**
 * @file PetscFEThermal.hpp
 * @brief PetscDS pointwise callbacks for thermal diffusion and thermoelastic coupling
 *
 * Temperature PDE (heat equation):
 *   rho * C_p * dT/dt - div(kappa * grad(T)) = Q_heat
 *
 * Thermoelastic coupling: thermal stress
 *   sigma = C : (eps - alpha_T * (T - T_ref) * I)
 *
 * Constants layout (unified array, slots 74-78):
 *   [74] = thermal conductivity kappa (W/(m*K))
 *   [75] = specific heat C_p (J/(kg*K))
 *   [76] = thermal expansion coefficient alpha_T (1/K)
 *   [77] = reference temperature T_ref (K)
 *   [78] = thermal field index in solution (used by thermoelastic callback)
 */

#include <petsc.h>

namespace FSRM {
namespace PetscFEThermal {

// Constants indices for thermal parameters
static constexpr PetscInt THERMAL_CONST_KAPPA   = 74;
static constexpr PetscInt THERMAL_CONST_CP      = 75;
static constexpr PetscInt THERMAL_CONST_ALPHA_T = 76;
static constexpr PetscInt THERMAL_CONST_T_REF   = 77;
static constexpr PetscInt THERMAL_CONST_FIELD_IDX = 78;

// =========================================================================
// Pure thermal PDE callbacks (heat equation on the temperature field)
// =========================================================================

/**
 * @brief f0 residual: rho * Cp * dT/dt (accumulation term)
 *
 * The temperature field is assumed to be the field at uOff index matching
 * the thermal field index stored in constants[THERMAL_CONST_FIELD_IDX].
 * rho is read from constants[2], Cp from constants[THERMAL_CONST_CP].
 */
void f0_thermal(PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[],
    const PetscScalar u_x[], const PetscInt aOff[],
    const PetscInt aOff_x[], const PetscScalar a[],
    const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, const PetscReal x[], PetscInt numConstants,
    const PetscScalar constants[], PetscScalar f0[]);

/**
 * @brief f1 residual: kappa * grad(T) (diffusion flux)
 */
void f1_thermal(PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[],
    const PetscScalar u_x[], const PetscInt aOff[],
    const PetscInt aOff_x[], const PetscScalar a[],
    const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, const PetscReal x[], PetscInt numConstants,
    const PetscScalar constants[], PetscScalar f1[]);

/**
 * @brief g0 Jacobian: rho * Cp * shift (time derivative coupling)
 */
void g0_thermal_thermal(PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[],
    const PetscScalar u_x[], const PetscInt aOff[],
    const PetscInt aOff_x[], const PetscScalar a[],
    const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, PetscReal u_tShift, const PetscReal x[],
    PetscInt numConstants, const PetscScalar constants[],
    PetscScalar g0[]);

/**
 * @brief g3 Jacobian: kappa * delta_ij (diffusion stiffness)
 */
void g3_thermal_thermal(PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[],
    const PetscScalar u_x[], const PetscInt aOff[],
    const PetscInt aOff_x[], const PetscScalar a[],
    const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, PetscReal u_tShift, const PetscReal x[],
    PetscInt numConstants, const PetscScalar constants[],
    PetscScalar g3[]);

// =========================================================================
// Thermoelastic coupling: stress with thermal strain
// =========================================================================

/**
 * @brief f1 for displacement with thermal strain subtraction
 *
 * sigma_ij = lambda * (eps_kk - 3*alpha_T*dT) * delta_ij
 *          + 2 * mu * (eps_ij - alpha_T * dT * delta_ij)
 *
 * where dT = T - T_ref. Temperature is read from the solution at the
 * thermal field index. Material parameters from constants[0]=lambda,
 * constants[1]=mu.
 */
void f1_thermoelastic(PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[],
    const PetscScalar u_x[], const PetscInt aOff[],
    const PetscInt aOff_x[], const PetscScalar a[],
    const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, const PetscReal x[], PetscInt numConstants,
    const PetscScalar constants[], PetscScalar f1[]);

/**
 * @brief f1 for displacement with thermal strain using auxiliary fields
 *
 * Same as f1_thermoelastic but reads lambda, mu from aux fields.
 * Temperature still from solution field.
 */
void f1_thermoelastic_aux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[],
    const PetscScalar u_x[], const PetscInt aOff[],
    const PetscInt aOff_x[], const PetscScalar a[],
    const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, const PetscReal x[], PetscInt numConstants,
    const PetscScalar constants[], PetscScalar f1[]);

} // namespace PetscFEThermal
} // namespace FSRM

#endif // PETSCFE_THERMAL_HPP
