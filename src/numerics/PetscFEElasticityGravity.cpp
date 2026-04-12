#include "numerics/PetscFEElasticityGravity.hpp"
#include "numerics/PetscFEElasticityAux.hpp"

namespace FSRM
{

// Gravity body force: delegates to PetscFEElasticityAux::f0_elastostatics_aux
void PetscFEElasticityGravity::f0_gravity(PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
    const PetscScalar u_t[], const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[],
    const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, const PetscReal x[], PetscInt numConstants,
    const PetscScalar constants[], PetscScalar f0[])
{
  PetscFEElasticityAux::f0_elastostatics_aux(dim, Nf, NfAux,
      uOff, uOff_x, u, u_t, u_x,
      aOff, aOff_x, a, a_t, a_x,
      t, x, numConstants, constants, f0);
}

// Stress divergence: delegates to PetscFEElasticityAux::f1_elastostatics_aux
void PetscFEElasticityGravity::f1_stress(PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
    const PetscScalar u_t[], const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[],
    const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, const PetscReal x[], PetscInt numConstants,
    const PetscScalar constants[], PetscScalar f1[])
{
  PetscFEElasticityAux::f1_elastostatics_aux(dim, Nf, NfAux,
      uOff, uOff_x, u, u_t, u_x,
      aOff, aOff_x, a, a_t, a_x,
      t, x, numConstants, constants, f1);
}

// Elasticity Jacobian: delegates to PetscFEElasticityAux::g3_elastostatics_aux
void PetscFEElasticityGravity::g3_elasticity(PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
    const PetscScalar u_t[], const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[],
    const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants,
    const PetscScalar constants[], PetscScalar g3[])
{
  PetscFEElasticityAux::g3_elastostatics_aux(dim, Nf, NfAux,
      uOff, uOff_x, u, u_t, u_x,
      aOff, aOff_x, a, a_t, a_x,
      t, u_tShift, x, numConstants, constants, g3);
}

} // namespace FSRM
