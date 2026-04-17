#ifndef PETSCFE_ELASTICITY_GRAVITY_HPP
#define PETSCFE_ELASTICITY_GRAVITY_HPP

/**
 * @file PetscFEElasticityGravity.hpp
 * @brief Gravity body force callback for elastostatics
 *
 * Provides a gravity body force callback that reads density (rho)
 * from the auxiliary field DM and gravity magnitude from constants[0].
 * This is a public interface to the gravity functionality already
 * present in PetscFEElasticityAux.
 *
 * The auxiliary field layout follows AuxFieldIndex:
 *   a[aOff[AUX_LAMBDA]] = lambda
 *   a[aOff[AUX_MU]]     = mu
 *   a[aOff[AUX_RHO]]    = rho
 *
 * Constants:
 *   constants[0] = gravity magnitude (m/s^2), positive downward
 *
 * Reference: For a gravity-loaded column with height H:
 *   sigma_zz(z) = -rho * g * (H - z)
 *   sigma_xx = sigma_yy = K0 * sigma_zz
 *   where K0 = nu / (1 - nu) for elastic uniaxial strain
 */

#include <petsc.h>

namespace FSRM
{

class PetscFEElasticityGravity
{
public:
  /**
   * Gravity body force residual f0.
   * Reads rho from aux field, gravity from constants[0].
   * f0[dim-1] = -rho * g  (gravity in -z direction)
   */
  static void f0_gravity(PetscInt dim, PetscInt Nf, PetscInt NfAux,
      const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
      const PetscScalar u_t[], const PetscScalar u_x[],
      const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[],
      const PetscScalar a_t[], const PetscScalar a_x[],
      PetscReal t, const PetscReal x[], PetscInt numConstants,
      const PetscScalar constants[], PetscScalar f0[]);

  /**
   * Stress divergence residual f1 using aux lambda, mu.
   * Delegates to PetscFEElasticityAux::f1_elastostatics_aux.
   */
  static void f1_stress(PetscInt dim, PetscInt Nf, PetscInt NfAux,
      const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
      const PetscScalar u_t[], const PetscScalar u_x[],
      const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[],
      const PetscScalar a_t[], const PetscScalar a_x[],
      PetscReal t, const PetscReal x[], PetscInt numConstants,
      const PetscScalar constants[], PetscScalar f1[]);

  /**
   * Elasticity Jacobian g3 using aux lambda, mu.
   * Delegates to PetscFEElasticityAux::g3_elastostatics_aux.
   */
  static void g3_elasticity(PetscInt dim, PetscInt Nf, PetscInt NfAux,
      const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
      const PetscScalar u_t[], const PetscScalar u_x[],
      const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[],
      const PetscScalar a_t[], const PetscScalar a_x[],
      PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants,
      const PetscScalar constants[], PetscScalar g3[]);
};

} // namespace FSRM

#endif // PETSCFE_ELASTICITY_GRAVITY_HPP
