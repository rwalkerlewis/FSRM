#ifndef PETSC_FE_FLUID_FLOW_HPP
#define PETSC_FE_FLUID_FLOW_HPP

#include <petsc.h>

namespace FSRM {
namespace PetscFEFluidFlow {

// -----------------------------------------------------------------------------
// PETScFE (PetscDS) pointwise callbacks for fluid flow
// -----------------------------------------------------------------------------

// Single-phase pressure: phi*ct*dP/dt - div( (k/mu) * grad(P) ) = 0
void f0_SinglePhase(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[],
                    const PetscScalar u[], const PetscScalar u_t[],
                    const PetscScalar u_x[], const PetscInt aOff[],
                    const PetscInt aOff_x[], const PetscScalar a[],
                    const PetscScalar a_x[], const PetscScalar a_t[],
                    PetscReal t,
                    const PetscReal x[], PetscInt numConstants,
                    const PetscScalar constants[], PetscScalar f0[]);

void f1_SinglePhase(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[],
                    const PetscScalar u[], const PetscScalar u_t[],
                    const PetscScalar u_x[], const PetscInt aOff[],
                    const PetscInt aOff_x[], const PetscScalar a[],
                    const PetscScalar a_x[], const PetscScalar a_t[],
                    PetscReal t,
                    const PetscReal x[], PetscInt numConstants,
                    const PetscScalar constants[], PetscScalar f1[]);

// Multiphase (black-oil) pressure:
//   phi * ct_eff(Sw,Sg) * dP/dt - div( K * lambda_t(Sw,Sg) * grad(P) ) = 0
void f0_BlackOilPressure(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                         const PetscInt uOff[], const PetscInt uOff_x[],
                         const PetscScalar u[], const PetscScalar u_t[],
                         const PetscScalar u_x[], const PetscInt aOff[],
                         const PetscInt aOff_x[], const PetscScalar a[],
                         const PetscScalar a_x[], const PetscScalar a_t[],
                         PetscReal t,
                         const PetscReal x[], PetscInt numConstants,
                         const PetscScalar constants[], PetscScalar f0[]);

void f1_BlackOilPressure(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                         const PetscInt uOff[], const PetscInt uOff_x[],
                         const PetscScalar u[], const PetscScalar u_t[],
                         const PetscScalar u_x[], const PetscInt aOff[],
                         const PetscInt aOff_x[], const PetscScalar a[],
                         const PetscScalar a_x[], const PetscScalar a_t[],
                         PetscReal t,
                         const PetscReal x[], PetscInt numConstants,
                         const PetscScalar constants[], PetscScalar f1[]);

// Jacobians for pressure equation (block (P,P))
void g0_BlackOilPressurePressure(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                                 const PetscInt uOff[], const PetscInt uOff_x[],
                                 const PetscScalar u[], const PetscScalar u_t[],
                                 const PetscScalar u_x[], const PetscInt aOff[],
                                 const PetscInt aOff_x[], const PetscScalar a[],
                                 const PetscScalar a_x[], const PetscScalar a_t[],
                                 PetscReal t, PetscReal u_tShift,
                                 const PetscReal x[], PetscInt numConstants,
                                 const PetscScalar constants[], PetscScalar g0[]);

void g3_BlackOilPressurePressure(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                                 const PetscInt uOff[], const PetscInt uOff_x[],
                                 const PetscScalar u[], const PetscScalar u_t[],
                                 const PetscScalar u_x[], const PetscInt aOff[],
                                 const PetscInt aOff_x[], const PetscScalar a[],
                                 const PetscScalar a_x[], const PetscScalar a_t[],
                                 PetscReal t, PetscReal u_tShift,
                                 const PetscReal x[], PetscInt numConstants,
                                 const PetscScalar constants[], PetscScalar g3[]);

// Safe placeholders for saturation fields (keep constant in time)
void f0_BlackOilSw(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                   const PetscInt uOff[], const PetscInt uOff_x[],
                   const PetscScalar u[], const PetscScalar u_t[],
                   const PetscScalar u_x[], const PetscInt aOff[],
                   const PetscInt aOff_x[], const PetscScalar a[],
                   const PetscScalar a_x[], const PetscScalar a_t[],
                   PetscReal t,
                   const PetscReal x[], PetscInt numConstants,
                   const PetscScalar constants[], PetscScalar f0[]);

void f0_BlackOilSg(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                   const PetscInt uOff[], const PetscInt uOff_x[],
                   const PetscScalar u[], const PetscScalar u_t[],
                   const PetscScalar u_x[], const PetscInt aOff[],
                   const PetscInt aOff_x[], const PetscScalar a[],
                   const PetscScalar a_x[], const PetscScalar a_t[],
                   PetscReal t,
                   const PetscReal x[], PetscInt numConstants,
                   const PetscScalar constants[], PetscScalar f0[]);

void f1_Zero(PetscInt dim, PetscInt Nf, PetscInt NfAux,
             const PetscInt uOff[], const PetscInt uOff_x[],
             const PetscScalar u[], const PetscScalar u_t[],
             const PetscScalar u_x[], const PetscInt aOff[],
             const PetscInt aOff_x[], const PetscScalar a[],
             const PetscScalar a_x[], const PetscScalar a_t[],
             PetscReal t,
             const PetscReal x[], PetscInt numConstants,
             const PetscScalar constants[], PetscScalar f1[]);

void g0_Identity(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                 const PetscInt uOff[], const PetscInt uOff_x[],
                 const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscInt aOff[],
                 const PetscInt aOff_x[], const PetscScalar a[],
                 const PetscScalar a_x[], const PetscScalar a_t[],
                 PetscReal t, PetscReal u_tShift,
                 const PetscReal x[], PetscInt numConstants,
                 const PetscScalar constants[], PetscScalar g0[]);

} // namespace PetscFEFluidFlow
} // namespace FSRM

#endif // PETSC_FE_FLUID_FLOW_HPP

