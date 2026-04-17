#ifndef PETSCFE_ELASTICITY_AUX_HPP
#define PETSCFE_ELASTICITY_AUX_HPP

#include <petsc.h>

namespace FSRM
{

// Auxiliary field indices for heterogeneous material properties
enum AuxFieldIndex
{
  AUX_LAMBDA = 0,
  AUX_MU     = 1,
  AUX_RHO    = 2,
  NUM_AUX_FIELDS = 3
};

/**
 * @brief Elasticity callbacks that read material properties from auxiliary fields
 *        instead of the constants array.
 *
 * For homogeneous fallback, the constants array is used:
 *   constants[0] = gravity magnitude (0 if disabled)
 *
 * Auxiliary fields (per-cell):
 *   a[aOff[AUX_LAMBDA]] = lambda (first Lame parameter, Pa)
 *   a[aOff[AUX_MU]]     = mu (shear modulus, Pa)
 *   a[aOff[AUX_RHO]]    = rho (density, kg/m^3)
 */
class PetscFEElasticityAux
{
public:
  // Residual: body force (gravity) using aux rho
  static void f0_elastostatics_aux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
      const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
      const PetscScalar u_t[], const PetscScalar u_x[],
      const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[],
      const PetscScalar a_t[], const PetscScalar a_x[],
      PetscReal t, const PetscReal x[], PetscInt numConstants,
      const PetscScalar constants[], PetscScalar f0[]);

  // Residual: stress divergence using aux lambda, mu
  static void f1_elastostatics_aux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
      const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
      const PetscScalar u_t[], const PetscScalar u_x[],
      const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[],
      const PetscScalar a_t[], const PetscScalar a_x[],
      PetscReal t, const PetscReal x[], PetscInt numConstants,
      const PetscScalar constants[], PetscScalar f1[]);

  // Jacobian: g3 block using aux lambda, mu
  static void g3_elastostatics_aux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
      const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
      const PetscScalar u_t[], const PetscScalar u_x[],
      const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[],
      const PetscScalar a_t[], const PetscScalar a_x[],
      PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants,
      const PetscScalar constants[], PetscScalar g3[]);

  // Elastodynamics residual: inertia (rho * u_t) using aux rho
  static void f0_elastodynamics_aux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
      const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
      const PetscScalar u_t[], const PetscScalar u_x[],
      const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[],
      const PetscScalar a_t[], const PetscScalar a_x[],
      PetscReal t, const PetscReal x[], PetscInt numConstants,
      const PetscScalar constants[], PetscScalar f0[]);

  // Elastodynamics Jacobian: mass matrix (rho * shift * I) using aux rho
  static void g0_elastodynamics_aux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
      const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
      const PetscScalar u_t[], const PetscScalar u_x[],
      const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[],
      const PetscScalar a_t[], const PetscScalar a_x[],
      PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants,
      const PetscScalar constants[], PetscScalar g0[]);
};

} // namespace FSRM

#endif // PETSCFE_ELASTICITY_AUX_HPP
