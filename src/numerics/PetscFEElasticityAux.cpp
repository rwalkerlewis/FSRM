#include "numerics/PetscFEElasticityAux.hpp"

namespace FSRM
{

// Residual f0: body force from gravity using auxiliary rho
// constants[0] = gravity magnitude (positive value, e.g. 9.81; 0 if disabled)
//
// PETSc FEM convention: the weak form residual is
//   F = integral(f0 . v) + integral(f1 : grad(v)) = 0
// Integration by parts gives: div(sigma) = f0 (strong form).
// The equilibrium PDE is: div(sigma) + b = 0, where b = -rho*g*ez (gravity).
// So f0 = -b = rho*g*ez, i.e. f0[dim-1] = +rho*g (positive upward).
void PetscFEElasticityAux::f0_elastostatics_aux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
    const PetscScalar u_t[], const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[],
    const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, const PetscReal x[], PetscInt numConstants,
    const PetscScalar constants[], PetscScalar f0[])
{
  const PetscScalar rho = a[aOff[AUX_RHO]];
  for (PetscInt d = 0; d < dim; ++d) f0[d] = 0.0;
  // Gravity body force: f0 = -b = rho*g (positive, since b = -rho*g*ez)
  if (numConstants > 0 && PetscRealPart(constants[0]) != 0.0)
  {
    f0[dim - 1] = rho * constants[0];
  }
}

// Residual f1: stress divergence sigma = lambda*tr(eps)*I + 2*mu*eps
// Displacement is field 0 for elastostatics (fluid_model = NONE)
void PetscFEElasticityAux::f1_elastostatics_aux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
    const PetscScalar u_t[], const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[],
    const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, const PetscReal x[], PetscInt numConstants,
    const PetscScalar constants[], PetscScalar f1[])
{
  const PetscScalar lambda = a[aOff[AUX_LAMBDA]];
  const PetscScalar mu     = a[aOff[AUX_MU]];
  const PetscInt    off    = uOff_x[0];

  // Compute volumetric strain
  PetscScalar trace = 0.0;
  for (PetscInt d = 0; d < dim; ++d)
  {
    trace += u_x[off + d * dim + d];
  }

  // Assemble stress tensor: sigma_ij = lambda * eps_kk * delta_ij + 2*mu*eps_ij
  for (PetscInt i = 0; i < dim; ++i)
  {
    for (PetscInt j = 0; j < dim; ++j)
    {
      const PetscScalar eps_ij = 0.5 * (u_x[off + i * dim + j] + u_x[off + j * dim + i]);
      f1[i * dim + j] = 2.0 * mu * eps_ij;
      if (i == j) f1[i * dim + j] += lambda * trace;
    }
  }
}

// Jacobian g3: elasticity tensor C_ijkl = lambda*delta_ij*delta_kl + mu*(delta_ik*delta_jl + delta_il*delta_jk)
void PetscFEElasticityAux::g3_elastostatics_aux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
    const PetscScalar u_t[], const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[],
    const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants,
    const PetscScalar constants[], PetscScalar g3[])
{
  const PetscScalar lambda = a[aOff[AUX_LAMBDA]];
  const PetscScalar mu     = a[aOff[AUX_MU]];

  for (PetscInt i = 0; i < dim; ++i)
  {
    for (PetscInt j = 0; j < dim; ++j)
    {
      for (PetscInt k = 0; k < dim; ++k)
      {
        for (PetscInt l = 0; l < dim; ++l)
        {
          g3[((i * dim + j) * dim + k) * dim + l] =
              lambda * (i == j ? 1.0 : 0.0) * (k == l ? 1.0 : 0.0)
            + mu * ((i == k ? 1.0 : 0.0) * (j == l ? 1.0 : 0.0)
                  + (i == l ? 1.0 : 0.0) * (j == k ? 1.0 : 0.0));
        }
      }
    }
  }
}

// Elastodynamics residual f0: inertia term rho * u_t using auxiliary rho
void PetscFEElasticityAux::f0_elastodynamics_aux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
    const PetscScalar u_t[], const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[],
    const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, const PetscReal x[], PetscInt numConstants,
    const PetscScalar constants[], PetscScalar f0[])
{
  (void)Nf; (void)NfAux; (void)uOff_x; (void)u; (void)u_x;
  (void)aOff_x; (void)a_t; (void)a_x; (void)t; (void)x; (void)numConstants;

  const PetscScalar rho = a[aOff[AUX_RHO]];

  for (PetscInt c = 0; c < dim; ++c)
  {
    f0[c] = rho * u_t[uOff[0] + c];
  }
}

// Elastodynamics Jacobian g0: mass matrix rho * shift * I using auxiliary rho
void PetscFEElasticityAux::g0_elastodynamics_aux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
    const PetscScalar u_t[], const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[],
    const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants,
    const PetscScalar constants[], PetscScalar g0[])
{
  (void)Nf; (void)NfAux; (void)uOff; (void)uOff_x;
  (void)u; (void)u_t; (void)u_x;
  (void)aOff_x; (void)a_t; (void)a_x; (void)t; (void)x; (void)numConstants;

  const PetscScalar rho = a[aOff[AUX_RHO]];

  for (PetscInt i = 0; i < dim; ++i)
  {
    for (PetscInt j = 0; j < dim; ++j)
    {
      g0[i * dim + j] = (i == j) ? (rho * u_tShift) : 0.0;
    }
  }
}

} // namespace FSRM
