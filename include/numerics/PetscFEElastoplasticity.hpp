#ifndef PETSCFE_ELASTOPLASTICITY_HPP
#define PETSCFE_ELASTOPLASTICITY_HPP

#include <petsc.h>

namespace FSRM
{

/**
 * @brief Elastoplasticity callbacks for Drucker-Prager return mapping
 *        within PetscDS FEM assembly.
 *
 * Reads material properties from auxiliary fields (same as PetscFEElasticityAux):
 *   a[aOff[0]] = lambda (Pa)
 *   a[aOff[1]] = mu (Pa)
 *   a[aOff[2]] = rho (kg/m^3)
 *   a[aOff[3..8]] = plastic strain (Voigt: xx,yy,zz,yz,xz,xy)
 *   a[aOff[4]] = accumulated plastic strain (scalar)
 *
 * Yield parameters from constants array:
 *   constants[EP_CONST_COHESION] = cohesion (Pa)
 *   constants[EP_CONST_FRICTION_ANGLE] = friction angle (rad)
 *   constants[EP_CONST_DILATION_ANGLE] = dilation angle (rad)
 *   constants[EP_CONST_HARDENING_MODULUS] = hardening modulus (Pa)
 */
class PetscFEElastoplasticity
{
public:
  // Constants array indices for elastoplasticity parameters
  // Placed after cohesive fault (25-30) and hydrofrac (31) constants
  enum
  {
    EP_CONST_COHESION          = 32,
    EP_CONST_FRICTION_ANGLE    = 33,
    EP_CONST_DILATION_ANGLE    = 34,
    EP_CONST_HARDENING_MODULUS = 35,
    EP_CONST_COUNT             = 36
  };

  // Aux field indices (extends PetscFEElasticityAux)
  enum
  {
    EP_AUX_LAMBDA              = 0,
    EP_AUX_MU                  = 1,
    EP_AUX_RHO                 = 2,
    EP_AUX_PLASTIC_STRAIN      = 3,  // 6-component vector
    EP_AUX_ACC_PLASTIC_STRAIN  = 4,  // scalar
    EP_NUM_AUX_FIELDS          = 5
  };

  /**
   * @brief f1 residual: elastoplastic stress via Drucker-Prager return mapping
   *
   * Computes trial stress sigma_trial = lambda*tr(eps)*I + 2*mu*eps,
   * checks yield, applies return mapping if needed, returns corrected stress.
   */
  static void f1_elastoplastic_aux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
      const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
      const PetscScalar u_t[], const PetscScalar u_x[],
      const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[],
      const PetscScalar a_t[], const PetscScalar a_x[],
      PetscReal t, const PetscReal x[], PetscInt numConstants,
      const PetscScalar constants[], PetscScalar f1[]);

  /**
   * @brief g3 Jacobian: elastic tangent stiffness (approximation)
   *
   * Uses the elastic stiffness tensor as an approximation to the
   * consistent tangent. This is first-order accurate and sufficient
   * for initial convergence testing.
   */
  static void g3_elastoplastic_aux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
      const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
      const PetscScalar u_t[], const PetscScalar u_x[],
      const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[],
      const PetscScalar a_t[], const PetscScalar a_x[],
      PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants,
      const PetscScalar constants[], PetscScalar g3[]);
};

} // namespace FSRM

#endif // PETSCFE_ELASTOPLASTICITY_HPP
