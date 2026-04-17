#include "numerics/AbsorbingBC.hpp"
#include "numerics/PetscFEElasticityAux.hpp"

namespace FSRM
{

void AbsorbingBC::f0_absorbing(PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
    const PetscScalar u_t[], const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[],
    const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, const PetscReal x[], const PetscReal n[],
    PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  // Read material properties from aux fields or constants
  PetscScalar lambda, mu, rho;
  if (NfAux >= 3 && a != NULL) {
    lambda = a[aOff[AUX_LAMBDA]];
    mu     = a[aOff[AUX_MU]];
    rho    = a[aOff[AUX_RHO]];
  } else {
    lambda = constants[0];
    mu     = constants[1];
    rho    = constants[2];
  }

  // Zero output
  for (PetscInt d = 0; d < dim; ++d) f0[d] = 0.0;

  // No velocity means no absorbing traction
  if (u_t == NULL) return;

  // Compute wave speeds
  const PetscReal cp = PetscSqrtReal(PetscRealPart((lambda + 2.0 * mu) / rho));
  const PetscReal cs = PetscSqrtReal(PetscRealPart(mu / rho));

  // Compute n dot v
  PetscScalar ndotv = 0.0;
  for (PetscInt d = 0; d < dim; ++d)
  {
    ndotv += n[d] * u_t[uOff[0] + d];
  }

  // Absorbing traction:
  //   f0_i = rho * cp * (n.v) * n_i + rho * cs * (v_i - (n.v) * n_i)
  for (PetscInt d = 0; d < dim; ++d)
  {
    f0[d] = rho * cp * ndotv * n[d]
          + rho * cs * (u_t[uOff[0] + d] - ndotv * n[d]);
  }
}

void AbsorbingBC::g0_absorbing(PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
    const PetscScalar u_t[], const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[],
    const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, PetscReal u_tShift, const PetscReal x[], const PetscReal n[],
    PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
{
  // Read material properties from aux fields or constants
  PetscScalar lambda, mu, rho;
  if (NfAux >= 3 && a != NULL) {
    lambda = a[aOff[AUX_LAMBDA]];
    mu     = a[aOff[AUX_MU]];
    rho    = a[aOff[AUX_RHO]];
  } else {
    lambda = constants[0];
    mu     = constants[1];
    rho    = constants[2];
  }

  // Compute wave speeds
  const PetscReal cp = PetscSqrtReal(PetscRealPart((lambda + 2.0 * mu) / rho));
  const PetscReal cs = PetscSqrtReal(PetscRealPart(mu / rho));

  // Jacobian of absorbing traction wrt u_t:
  //   g0_{ij} = u_tShift * (rho * cp * n_i * n_j + rho * cs * (delta_ij - n_i * n_j))
  for (PetscInt i = 0; i < dim; ++i)
  {
    for (PetscInt j = 0; j < dim; ++j)
    {
      g0[i * dim + j] = u_tShift * (
        rho * cp * n[i] * n[j]
        + rho * cs * ((i == j ? 1.0 : 0.0) - n[i] * n[j])
      );
    }
  }
}

} // namespace FSRM
