#ifndef PETSC_FE_ELASTICITY_HPP
#define PETSC_FE_ELASTICITY_HPP

#include <petsc.h>

namespace FSRM {
namespace PetscFEElasticity {

// -----------------------------------------------------------------------------
// PETScFE (PetscDS) pointwise callbacks for linear elasticity
// -----------------------------------------------------------------------------

// Elastostatics: -div(sigma) = f_body
// where sigma_{ij} = lambda * eps_kk * delta_ij + 2*mu * eps_ij
// and eps_ij = 0.5*(du_i/dx_j + du_j/dx_i)
//
// Constants layout:
//   constants[0] = lambda (first Lame parameter)
//   constants[1] = mu (shear modulus, second Lame parameter)
//   constants[2] = rho (density, used by elastodynamics)
//
// Field layout: single vector field with 3 components (ux, uy, uz)
//   u[uOff[f] + c] = component c of field f
//   u_x[uOff_x[f] + c*dim + d] = d(u_c)/dx_d

// Elastostatics residual
void f0_elastostatics(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                      const PetscInt uOff[], const PetscInt uOff_x[],
                      const PetscScalar u[], const PetscScalar u_t[],
                      const PetscScalar u_x[], const PetscInt aOff[],
                      const PetscInt aOff_x[], const PetscScalar a[],
                      const PetscScalar a_x[], const PetscScalar a_t[],
                      PetscReal t,
                      const PetscReal x[], PetscInt numConstants,
                      const PetscScalar constants[], PetscScalar f0[]);

void f1_elastostatics(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                      const PetscInt uOff[], const PetscInt uOff_x[],
                      const PetscScalar u[], const PetscScalar u_t[],
                      const PetscScalar u_x[], const PetscInt aOff[],
                      const PetscInt aOff_x[], const PetscScalar a[],
                      const PetscScalar a_x[], const PetscScalar a_t[],
                      PetscReal t,
                      const PetscReal x[], PetscInt numConstants,
                      const PetscScalar constants[], PetscScalar f1[]);

// Elastostatics Jacobian
void g3_elastostatics(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                      const PetscInt uOff[], const PetscInt uOff_x[],
                      const PetscScalar u[], const PetscScalar u_t[],
                      const PetscScalar u_x[], const PetscInt aOff[],
                      const PetscInt aOff_x[], const PetscScalar a[],
                      const PetscScalar a_x[], const PetscScalar a_t[],
                      PetscReal t, PetscReal u_tShift,
                      const PetscReal x[], PetscInt numConstants,
                      const PetscScalar constants[], PetscScalar g3[]);

// Elastodynamics: rho*u_tt - div(sigma) = f_body
// Adds inertial term to f0
void f0_elastodynamics(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[],
                       const PetscScalar u[], const PetscScalar u_t[],
                       const PetscScalar u_x[], const PetscInt aOff[],
                       const PetscInt aOff_x[], const PetscScalar a[],
                       const PetscScalar a_x[], const PetscScalar a_t[],
                       PetscReal t,
                       const PetscReal x[], PetscInt numConstants,
                       const PetscScalar constants[], PetscScalar f0[]);

void f1_elastodynamics(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[],
                       const PetscScalar u[], const PetscScalar u_t[],
                       const PetscScalar u_x[], const PetscInt aOff[],
                       const PetscInt aOff_x[], const PetscScalar a[],
                       const PetscScalar a_x[], const PetscScalar a_t[],
                       PetscReal t,
                       const PetscReal x[], PetscInt numConstants,
                       const PetscScalar constants[], PetscScalar f1[]);

// Elastodynamics Jacobian (adds mass matrix g0 term)
void g0_elastodynamics(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[],
                       const PetscScalar u[], const PetscScalar u_t[],
                       const PetscScalar u_x[], const PetscInt aOff[],
                       const PetscInt aOff_x[], const PetscScalar a[],
                       const PetscScalar a_x[], const PetscScalar a_t[],
                       PetscReal t, PetscReal u_tShift,
                       const PetscReal x[], PetscInt numConstants,
                       const PetscScalar constants[], PetscScalar g0[]);

void g3_elastodynamics(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[],
                       const PetscScalar u[], const PetscScalar u_t[],
                       const PetscScalar u_x[], const PetscInt aOff[],
                       const PetscInt aOff_x[], const PetscScalar a[],
                       const PetscScalar a_x[], const PetscScalar a_t[],
                       PetscReal t, PetscReal u_tShift,
                       const PetscReal x[], PetscInt numConstants,
                       const PetscScalar constants[], PetscScalar g3[]);

// -----------------------------------------------------------------------------
// I2-form callbacks for second-order-in-time dynamics (TSALPHA2)
// For use with TSSetI2Function / TSSetI2Jacobian
// F(t, U, U_t, U_tt) = M*U_tt + K*U - f = 0
// -----------------------------------------------------------------------------

// I2Function residual: rho*u_tt - div(sigma) = 0
void f0_elastodynamics_I2(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                          const PetscInt uOff[], const PetscInt uOff_x[],
                          const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_tt[],
                          const PetscScalar u_x[], const PetscInt aOff[],
                          const PetscInt aOff_x[], const PetscScalar a[],
                          const PetscScalar a_x[], const PetscScalar a_t[],
                          PetscReal t,
                          const PetscReal x[], PetscInt numConstants,
                          const PetscScalar constants[], PetscScalar f0[]);

void f1_elastodynamics_I2(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                          const PetscInt uOff[], const PetscInt uOff_x[],
                          const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_tt[],
                          const PetscScalar u_x[], const PetscInt aOff[],
                          const PetscInt aOff_x[], const PetscScalar a[],
                          const PetscScalar a_x[], const PetscScalar a_t[],
                          PetscReal t,
                          const PetscReal x[], PetscInt numConstants,
                          const PetscScalar constants[], PetscScalar f1[]);

// I2Jacobian: mass matrix and stiffness matrix
void g0_elastodynamics_I2(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                          const PetscInt uOff[], const PetscInt uOff_x[],
                          const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_tt[],
                          const PetscScalar u_x[], const PetscInt aOff[],
                          const PetscInt aOff_x[], const PetscScalar a[],
                          const PetscScalar a_x[], const PetscScalar a_t[],
                          PetscReal t, PetscReal u_tShift, PetscReal u_ttShift,
                          const PetscReal x[], PetscInt numConstants,
                          const PetscScalar constants[], PetscScalar g0[]);

void g3_elastodynamics_I2(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                          const PetscInt uOff[], const PetscInt uOff_x[],
                          const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_tt[],
                          const PetscScalar u_x[], const PetscInt aOff[],
                          const PetscInt aOff_x[], const PetscScalar a[],
                          const PetscScalar a_x[], const PetscScalar a_t[],
                          PetscReal t, PetscReal u_tShift, PetscReal u_ttShift,
                          const PetscReal x[], PetscInt numConstants,
                          const PetscScalar constants[], PetscScalar g3[]);

} // namespace PetscFEElasticity
} // namespace FSRM

#endif // PETSC_FE_ELASTICITY_HPP
