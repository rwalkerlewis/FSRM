#ifndef FSRM_ABSORBING_BC_HPP
#define FSRM_ABSORBING_BC_HPP

#include <petsc.h>

namespace FSRM
{

/**
 * @brief First-order Clayton-Engquist absorbing boundary conditions
 *
 * On each boundary face, the absorbing condition adds a damping traction:
 *   t_i = rho * cp * (n_i * n_j) * v_j + rho * cs * (delta_ij - n_i * n_j) * v_j
 *
 * where v_j is particle velocity (u_t), n_i is outward normal,
 * cp = sqrt((lambda + 2*mu) / rho) is P-wave speed,
 * cs = sqrt(mu / rho) is S-wave speed.
 *
 * Material properties come from auxiliary fields if NfAux >= 3,
 * otherwise from constants[0..2] (lambda, mu, rho).
 */
class AbsorbingBC
{
public:
  /**
   * @brief Boundary residual f0: absorbing traction
   *
   * PetscBdPointFunc signature. Applied to DM_BC_NATURAL boundaries.
   */
  static void f0_absorbing(PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
    const PetscScalar u_t[], const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[],
    const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, const PetscReal x[], const PetscReal n[],
    PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[]);

  /**
   * @brief Boundary Jacobian g0: linearization of absorbing traction wrt u_t
   *
   * PetscBdPointJac signature. For implicit time stepping.
   * g0_{ij} = u_tShift * (rho * cp * n_i * n_j + rho * cs * (delta_ij - n_i * n_j))
   */
  static void g0_absorbing(PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
    const PetscScalar u_t[], const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[],
    const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, PetscReal u_tShift, const PetscReal x[], const PetscReal n[],
    PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[]);
};

} // namespace FSRM

#endif // FSRM_ABSORBING_BC_HPP
