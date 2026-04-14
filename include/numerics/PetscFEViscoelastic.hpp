#ifndef PETSC_FE_VISCOELASTIC_HPP
#define PETSC_FE_VISCOELASTIC_HPP

/**
 * @file PetscFEViscoelastic.hpp
 * @brief PetscDS pointwise callbacks for viscoelastic attenuation
 *
 * Implements the generalized Maxwell body (GMB) for seismic attenuation.
 * The stress-strain relation is:
 *
 *   sigma(t) = C_relaxed : eps(t) + sum_i R_i(t)
 *
 * where R_i are memory variables satisfying:
 *
 *   dR_i/dt = -R_i / tau_i + delta_mu_i * deps/dt
 *
 * The memory variables R_i are stored as auxiliary fields (6 components each,
 * symmetric tensor in Voigt notation: xx, yy, zz, yz, xz, xy).
 *
 * Auxiliary field layout (all scalar degree-0 Lagrange):
 *   aOff[0] = lambda (relaxed first Lame parameter)
 *   aOff[1] = mu (relaxed shear modulus)
 *   aOff[2] = rho (density)
 *   aOff[3]  = R_0_xx,  aOff[4]  = R_0_yy,  ...  aOff[8]  = R_0_xy
 *   aOff[9]  = R_1_xx,  aOff[10] = R_1_yy,  ...  aOff[14] = R_1_xy
 *   ...
 *   aOff[3 + m*6 + v] = R_m component v (Voigt: 0=xx,1=yy,2=zz,3=yz,4=xz,5=xy)
 *
 * Constants layout (starting at VISCO_CONST_BASE = 54):
 *   [54] = num_mechanisms N
 *   [55..59] = tau_i (relaxation times, up to 5 mechanisms)
 *   [60..64] = delta_mu_i (shear modulus weights)
 *   [65..69] = delta_kappa_i (bulk modulus weights)
 */

#include <petsc.h>

namespace FSRM {
namespace PetscFEViscoelastic {

// Constants array layout for viscoelastic parameters
static constexpr PetscInt VISCO_CONST_BASE     = 54;
static constexpr PetscInt VISCO_CONST_N        = 54;  // Number of mechanisms
static constexpr PetscInt VISCO_CONST_TAU_BASE = 55;  // tau_0..tau_4
static constexpr PetscInt VISCO_CONST_DMU_BASE = 60;  // delta_mu_0..delta_mu_4
static constexpr PetscInt VISCO_CONST_DK_BASE  = 65;  // delta_kappa_0..delta_kappa_4
static constexpr PetscInt VISCO_CONST_COUNT    = 70;
static constexpr PetscInt VISCO_MAX_MECHANISMS = 5;

// Auxiliary field indices (memory variables start after lambda, mu, rho)
static constexpr PetscInt VISCO_AUX_MEMORY_BASE = 3;  // R_0 starts at aOff[3]

/**
 * @brief Compute relaxation times and anelastic weights from Q values.
 *
 * Log-spaces N relaxation frequencies between f_min and f_max, then
 * solves a least-squares problem to find weights that produce the
 * desired constant Q over the bandwidth.
 *
 * @param N Number of mechanisms (1..5)
 * @param f_min Low frequency bound (Hz)
 * @param f_max High frequency bound (Hz)
 * @param Q_target Desired quality factor
 * @param tau Output: relaxation times (size N)
 * @param delta_mu Output: shear modulus weights (size N)
 */
void computeMechanismWeights(int N, double f_min, double f_max, double Q_target,
                             double tau[], double delta_mu[]);

/**
 * @brief f1 callback: viscoelastic stress = C_relaxed : eps + sum_i R_i
 *
 * Reads relaxed moduli (lambda, mu) from auxiliary fields and memory
 * variables R_i from auxiliary fields 3..2+N. Computes the elastic
 * trial stress plus the viscoelastic correction from stored memory
 * variables.
 */
void f1_viscoelastic_aux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                         const PetscInt uOff[], const PetscInt uOff_x[],
                         const PetscScalar u[], const PetscScalar u_t[],
                         const PetscScalar u_x[], const PetscInt aOff[],
                         const PetscInt aOff_x[], const PetscScalar a[],
                         const PetscScalar a_x[], const PetscScalar a_t[],
                         PetscReal t,
                         const PetscReal x[], PetscInt numConstants,
                         const PetscScalar constants[], PetscScalar f1[]);

/**
 * @brief g3 callback: viscoelastic tangent modulus (unrelaxed approximation)
 *
 * Uses the unrelaxed modulus C_unrelaxed = C_relaxed + sum_i delta_C_i
 * as the tangent for implicit time integration. This is exact for
 * infinitesimal time steps and a good approximation for moderate steps.
 */
void g3_viscoelastic_aux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                         const PetscInt uOff[], const PetscInt uOff_x[],
                         const PetscScalar u[], const PetscScalar u_t[],
                         const PetscScalar u_x[], const PetscInt aOff[],
                         const PetscInt aOff_x[], const PetscScalar a[],
                         const PetscScalar a_x[], const PetscScalar a_t[],
                         PetscReal t, PetscReal u_tShift,
                         const PetscReal x[], PetscInt numConstants,
                         const PetscScalar constants[], PetscScalar g3[]);

} // namespace PetscFEViscoelastic
} // namespace FSRM

#endif // PETSC_FE_VISCOELASTIC_HPP
