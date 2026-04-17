/**
 * @file PetscFEElastoplasticity.cpp
 * @brief Drucker-Prager elastoplasticity callbacks for PetscDS FEM assembly
 *
 * Implements the return mapping algorithm inline within the f1 residual
 * callback. The trial stress is computed from the total strain (u_x),
 * the Drucker-Prager yield function is evaluated, and if the stress
 * exceeds yield, the return mapping projects it back to the yield surface.
 *
 * The plastic strain history is read from auxiliary fields. Note that
 * updating the aux fields after a converged time step requires a
 * separate post-solve callback (not implemented here).
 */

#include "numerics/PetscFEElastoplasticity.hpp"
#include <cmath>
#include <algorithm>

namespace FSRM
{

void PetscFEElastoplasticity::f1_elastoplastic_aux(PetscInt dim, PetscInt Nf,
    PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
    const PetscScalar u_t[], const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[],
    const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, const PetscReal x[], PetscInt numConstants,
    const PetscScalar constants[], PetscScalar f1[])
{
  // Read material properties from aux fields
  const PetscReal lambda = PetscRealPart(a[aOff[EP_AUX_LAMBDA]]);
  const PetscReal mu     = PetscRealPart(a[aOff[EP_AUX_MU]]);
  const PetscInt  off    = uOff_x[0];

  // Read yield parameters from constants
  const PetscReal cohesion        = (numConstants > EP_CONST_COHESION)
      ? PetscRealPart(constants[EP_CONST_COHESION]) : 1.0e6;
  const PetscReal friction_angle  = (numConstants > EP_CONST_FRICTION_ANGLE)
      ? PetscRealPart(constants[EP_CONST_FRICTION_ANGLE]) : 30.0 * M_PI / 180.0;
  const PetscReal dilation_angle  = (numConstants > EP_CONST_DILATION_ANGLE)
      ? PetscRealPart(constants[EP_CONST_DILATION_ANGLE]) : 0.0;
  const PetscReal hardening_mod   = (numConstants > EP_CONST_HARDENING_MODULUS)
      ? PetscRealPart(constants[EP_CONST_HARDENING_MODULUS]) : 0.0;

  // Drucker-Prager parameters (inner cone matching)
  const PetscReal sin_phi = PetscSinReal(friction_angle);
  const PetscReal cos_phi = PetscCosReal(friction_angle);
  const PetscReal alpha_dp = 6.0 * sin_phi / (PetscSqrtReal(3.0) * (3.0 - sin_phi));
  const PetscReal kappa_dp = 6.0 * cohesion * cos_phi
      / (PetscSqrtReal(3.0) * (3.0 - sin_phi));

  const PetscReal sin_psi = PetscSinReal(dilation_angle);
  const PetscReal alpha_g = (PetscAbsReal(3.0 - sin_psi) > 1e-15)
      ? 6.0 * sin_psi / (PetscSqrtReal(3.0) * (3.0 - sin_psi)) : 0.0;

  // Read accumulated plastic strain from aux field (if available)
  PetscReal acc_eps_p = 0.0;
  if (NfAux > EP_AUX_ACC_PLASTIC_STRAIN)
  {
    acc_eps_p = PetscRealPart(a[aOff[EP_AUX_ACC_PLASTIC_STRAIN]]);
  }

  // Compute strain tensor from displacement gradient
  // eps_ij = 0.5 * (u_{i,j} + u_{j,i})
  // Compute trial stress: sigma_ij = lambda * eps_kk * delta_ij + 2*mu*eps_ij
  PetscReal trace = 0.0;
  for (PetscInt d = 0; d < dim; ++d)
  {
    trace += PetscRealPart(u_x[off + d * dim + d]);
  }

  // Build trial stress in Voigt notation for 3D: [xx, yy, zz, yz, xz, xy]
  PetscReal trial_stress[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  if (dim == 3)
  {
    PetscReal eps[6];
    eps[0] = PetscRealPart(u_x[off + 0]);              // eps_xx
    eps[1] = PetscRealPart(u_x[off + 4]);              // eps_yy
    eps[2] = PetscRealPart(u_x[off + 8]);              // eps_zz
    eps[3] = 0.5 * (PetscRealPart(u_x[off + 7]) + PetscRealPart(u_x[off + 5])); // eps_yz
    eps[4] = 0.5 * (PetscRealPart(u_x[off + 6]) + PetscRealPart(u_x[off + 2])); // eps_xz
    eps[5] = 0.5 * (PetscRealPart(u_x[off + 3]) + PetscRealPart(u_x[off + 1])); // eps_xy

    trial_stress[0] = lambda * trace + 2.0 * mu * eps[0];
    trial_stress[1] = lambda * trace + 2.0 * mu * eps[1];
    trial_stress[2] = lambda * trace + 2.0 * mu * eps[2];
    trial_stress[3] = 2.0 * mu * eps[3];
    trial_stress[4] = 2.0 * mu * eps[4];
    trial_stress[5] = 2.0 * mu * eps[5];
  }
  else
  {
    // 2D: just compute sigma and write to f1 (elastic only for now)
    for (PetscInt i = 0; i < dim; ++i)
    {
      for (PetscInt j = 0; j < dim; ++j)
      {
        const PetscScalar eps_ij = 0.5
            * (u_x[off + i * dim + j] + u_x[off + j * dim + i]);
        f1[i * dim + j] = 2.0 * mu * eps_ij;
        if (i == j) f1[i * dim + j] += lambda * trace;
      }
    }
    return;
  }

  // Evaluate Drucker-Prager yield function
  // F = sqrt(J2) + alpha * I1 - kappa_h
  PetscReal I1_trial = trial_stress[0] + trial_stress[1] + trial_stress[2];
  PetscReal p_trial = I1_trial / 3.0;

  // Deviatoric stress
  PetscReal s[6];
  s[0] = trial_stress[0] - p_trial;
  s[1] = trial_stress[1] - p_trial;
  s[2] = trial_stress[2] - p_trial;
  s[3] = trial_stress[3];
  s[4] = trial_stress[4];
  s[5] = trial_stress[5];

  // J2 = 0.5 * s_ij * s_ij
  PetscReal J2_trial = 0.5 * (s[0]*s[0] + s[1]*s[1] + s[2]*s[2])
      + s[3]*s[3] + s[4]*s[4] + s[5]*s[5];
  PetscReal sqrtJ2_trial = PetscSqrtReal(PetscMax(0.0, J2_trial));

  PetscReal kappa_h = kappa_dp + hardening_mod * acc_eps_p;
  PetscReal F_trial = sqrtJ2_trial + alpha_dp * I1_trial - kappa_h;

  if (F_trial <= 0.0)
  {
    // Elastic: write trial stress to f1
    // f1 is in tensor form [dim x dim]
    f1[0]          = trial_stress[0];  // sigma_xx
    f1[dim + 1]    = trial_stress[1];  // sigma_yy
    f1[2*dim + 2]  = trial_stress[2];  // sigma_zz
    f1[dim + 2]    = trial_stress[3];  // sigma_yz
    f1[2*dim + 1]  = trial_stress[3];  // sigma_zy
    f1[2]          = trial_stress[4];  // sigma_xz
    f1[2*dim]      = trial_stress[4];  // sigma_zx
    f1[1]          = trial_stress[5];  // sigma_xy
    f1[dim]        = trial_stress[5];  // sigma_yx
    return;
  }

  // Plastic return mapping (Drucker-Prager spectral return)
  PetscReal K = lambda + 2.0 * mu / 3.0;  // Bulk modulus
  PetscReal G = mu;                         // Shear modulus

  PetscReal delta_gamma = 0.0;
  PetscReal residual = F_trial;
  const int max_iter = 20;
  const PetscReal tol = 1e-10;

  for (int iter = 0; iter < max_iter; ++iter)
  {
    PetscReal sqrtJ2 = sqrtJ2_trial - G * delta_gamma;
    PetscReal I1 = I1_trial - 9.0 * K * alpha_g * delta_gamma;
    PetscReal eps_p_new = acc_eps_p + delta_gamma;
    PetscReal kappa_iter = kappa_dp + hardening_mod * eps_p_new;
    residual = sqrtJ2 + alpha_dp * I1 - kappa_iter;

    if (PetscAbsReal(residual) < tol) break;

    PetscReal dF_dgamma = -G - 9.0 * K * alpha_dp * alpha_g - hardening_mod;
    PetscReal d_gamma = -residual / dF_dgamma;
    delta_gamma += d_gamma;
    delta_gamma = PetscMax(0.0, delta_gamma);
  }

  // Reconstruct corrected stress
  PetscReal sqrtJ2_corr = sqrtJ2_trial - G * delta_gamma;
  PetscReal I1_corr = I1_trial - 9.0 * K * alpha_g * delta_gamma;
  PetscReal p_corr = I1_corr / 3.0;

  PetscReal scale = (sqrtJ2_trial > 1e-15) ? sqrtJ2_corr / sqrtJ2_trial : 0.0;

  PetscReal sigma[6];
  sigma[0] = scale * s[0] + p_corr;
  sigma[1] = scale * s[1] + p_corr;
  sigma[2] = scale * s[2] + p_corr;
  sigma[3] = scale * s[3];
  sigma[4] = scale * s[4];
  sigma[5] = scale * s[5];

  // Write corrected stress to f1 [dim x dim] tensor
  f1[0]          = sigma[0];  // sigma_xx
  f1[dim + 1]    = sigma[1];  // sigma_yy
  f1[2*dim + 2]  = sigma[2];  // sigma_zz
  f1[dim + 2]    = sigma[3];  // sigma_yz
  f1[2*dim + 1]  = sigma[3];  // sigma_zy
  f1[2]          = sigma[4];  // sigma_xz
  f1[2*dim]      = sigma[4];  // sigma_zx
  f1[1]          = sigma[5];  // sigma_xy
  f1[dim]        = sigma[5];  // sigma_yx
}

void PetscFEElastoplasticity::g3_elastoplastic_aux(PetscInt dim, PetscInt Nf,
    PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
    const PetscScalar u_t[], const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[],
    const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants,
    const PetscScalar constants[], PetscScalar g3[])
{
  // Use elastic tangent as approximation to consistent tangent.
  // This gives linear convergence but is simple and correct.
  const PetscScalar lambda = a[aOff[EP_AUX_LAMBDA]];
  const PetscScalar mu     = a[aOff[EP_AUX_MU]];

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

} // namespace FSRM
