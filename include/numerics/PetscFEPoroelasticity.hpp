#ifndef PETSC_FE_POROELASTICITY_HPP
#define PETSC_FE_POROELASTICITY_HPP

#include <petsc.h>

namespace FSRM {
namespace PetscFEPoroelasticity {

// -----------------------------------------------------------------------------
// PETScFE (PetscDS) pointwise callbacks for Biot poroelasticity
// -----------------------------------------------------------------------------

// Biot poroelasticity: coupled fluid flow and solid deformation
//
// Field layout:
//   Field 0: pressure (scalar, 1 component)
//   Field 1: displacement (vector, 3 components: ux, uy, uz)
//
// UNIFIED PetscDS CONSTANTS LAYOUT (set once in Simulator::setupPhysics):
//   [0]=lambda [1]=mu [2]=rho_s [3]=phi [4]=kx [5]=ky [6]=kz
//   [7]=cw [8]=co [9]=cg [10]=mu_w [11]=mu_o [12]=mu_g
//   [13]=Swr [14]=Sor [15]=Sgr [16]=nw [17]=no [18]=ng
//   [19]=krw0 [20]=kro0 [21]=krg0 [22]=biot_alpha [23]=1/M [24]=rho_f
//   Total: 25 constants
//
// Poroelasticity-specific indices used in callbacks:
//   [0]  = lambda (first Lame parameter)
//   [1]  = mu (shear modulus)
//   [2]  = rho_solid
//   [3]  = porosity
//   [4]  = permeability kx (isotropic assumption, m^2)
//   [10] = mu_w (water/fluid viscosity, Pa*s)
//   [22] = biot_alpha (Biot coefficient)
//   [23] = 1/M (inverse Biot modulus)
//   [24] = rho_f (fluid density)
//
// Governing equations:
//   Pressure:     (1/M)*dp/dt + alpha*div(du/dt) + div((k/mu)*grad(p)) = 0
//   Displacement: -div(sigma - alpha*p*I) = 0  (quasi-static)
//
// where:
//   sigma_{ij} = lambda*eps_kk*delta_ij + 2*mu*eps_ij  (effective stress)
//   eps_ij = 0.5*(du_i/dx_j + du_j/dx_i)               (strain)
//
// PETSc gradient layout for vector field 1 (displacement):
//   u_x[uOff_x[1] + c*dim + d] = d(u_c)/dx_d
//
// PETSc gradient layout for scalar field 0 (pressure):
//   u_x[uOff_x[0] + d] = dp/dx_d

// -----------------------------------------------------------------------------
// Pressure equation residuals (Field 0)
// -----------------------------------------------------------------------------

// f0: (1/M)*dp/dt + alpha*div(du/dt)
void f0_pressure(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                 const PetscInt uOff[], const PetscInt uOff_x[],
                 const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscInt aOff[],
                 const PetscInt aOff_x[], const PetscScalar a[],
                 const PetscScalar a_x[], const PetscScalar a_t[],
                 PetscReal t,
                 const PetscReal x[], PetscInt numConstants,
                 const PetscScalar constants[], PetscScalar f0[]);

// f1: (k/mu)*grad(p)
void f1_pressure(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                 const PetscInt uOff[], const PetscInt uOff_x[],
                 const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscInt aOff[],
                 const PetscInt aOff_x[], const PetscScalar a[],
                 const PetscScalar a_x[], const PetscScalar a_t[],
                 PetscReal t,
                 const PetscReal x[], PetscInt numConstants,
                 const PetscScalar constants[], PetscScalar f1[]);

// -----------------------------------------------------------------------------
// Displacement equation residuals (Field 1)
// -----------------------------------------------------------------------------

// f0: 0 (quasi-static, no body force)
void f0_displacement(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                     const PetscInt uOff[], const PetscInt uOff_x[],
                     const PetscScalar u[], const PetscScalar u_t[],
                     const PetscScalar u_x[], const PetscInt aOff[],
                     const PetscInt aOff_x[], const PetscScalar a[],
                     const PetscScalar a_x[], const PetscScalar a_t[],
                     PetscReal t,
                     const PetscReal x[], PetscInt numConstants,
                     const PetscScalar constants[], PetscScalar f0[]);

// f1: sigma_{cd} - alpha*p*delta_{cd}
void f1_displacement(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                     const PetscInt uOff[], const PetscInt uOff_x[],
                     const PetscScalar u[], const PetscScalar u_t[],
                     const PetscScalar u_x[], const PetscInt aOff[],
                     const PetscInt aOff_x[], const PetscScalar a[],
                     const PetscScalar a_x[], const PetscScalar a_t[],
                     PetscReal t,
                     const PetscReal x[], PetscInt numConstants,
                     const PetscScalar constants[], PetscScalar f1[]);

// -----------------------------------------------------------------------------
// Jacobian callbacks
// -----------------------------------------------------------------------------

// g0_pp: d(f0_pressure)/d(p) = (1/M)*shift
void g0_pp(PetscInt dim, PetscInt Nf, PetscInt NfAux,
           const PetscInt uOff[], const PetscInt uOff_x[],
           const PetscScalar u[], const PetscScalar u_t[],
           const PetscScalar u_x[], const PetscInt aOff[],
           const PetscInt aOff_x[], const PetscScalar a[],
           const PetscScalar a_x[], const PetscScalar a_t[],
           PetscReal t, PetscReal u_tShift,
           const PetscReal x[], PetscInt numConstants,
           const PetscScalar constants[], PetscScalar g0[]);

// g3_pp: d(f1_pressure)/d(grad_p) = (k/mu)*delta_{dd}
void g3_pp(PetscInt dim, PetscInt Nf, PetscInt NfAux,
           const PetscInt uOff[], const PetscInt uOff_x[],
           const PetscScalar u[], const PetscScalar u_t[],
           const PetscScalar u_x[], const PetscInt aOff[],
           const PetscInt aOff_x[], const PetscScalar a[],
           const PetscScalar a_x[], const PetscScalar a_t[],
           PetscReal t, PetscReal u_tShift,
           const PetscReal x[], PetscInt numConstants,
           const PetscScalar constants[], PetscScalar g3[]);

// g2_pu: d(f1_pressure)/d(u_t) = -alpha*shift*delta_{dc}
// Biot coupling: f1[d] -= alpha*u_t[d], registered in g2 slot (test_grad, trial_basis)
void g2_pu(PetscInt dim, PetscInt Nf, PetscInt NfAux,
           const PetscInt uOff[], const PetscInt uOff_x[],
           const PetscScalar u[], const PetscScalar u_t[],
           const PetscScalar u_x[], const PetscInt aOff[],
           const PetscInt aOff_x[], const PetscScalar a[],
           const PetscScalar a_x[], const PetscScalar a_t[],
           PetscReal t, PetscReal u_tShift,
           const PetscReal x[], PetscInt numConstants,
           const PetscScalar constants[], PetscScalar g2[]);

// g2_up: d(f1_displacement)/d(p) = -alpha*delta_{cd}
// Coupling from pressure to displacement equation
void g2_up(PetscInt dim, PetscInt Nf, PetscInt NfAux,
           const PetscInt uOff[], const PetscInt uOff_x[],
           const PetscScalar u[], const PetscScalar u_t[],
           const PetscScalar u_x[], const PetscInt aOff[],
           const PetscInt aOff_x[], const PetscScalar a[],
           const PetscScalar a_x[], const PetscScalar a_t[],
           PetscReal t, PetscReal u_tShift,
           const PetscReal x[], PetscInt numConstants,
           const PetscScalar constants[], PetscScalar g2[]);

// g3_uu: d(f1_displacement)/d(grad_u) = C_{cdef}
// Elasticity stiffness tensor
void g3_uu(PetscInt dim, PetscInt Nf, PetscInt NfAux,
           const PetscInt uOff[], const PetscInt uOff_x[],
           const PetscScalar u[], const PetscScalar u_t[],
           const PetscScalar u_x[], const PetscInt aOff[],
           const PetscInt aOff_x[], const PetscScalar a[],
           const PetscScalar a_x[], const PetscScalar a_t[],
           PetscReal t, PetscReal u_tShift,
           const PetscReal x[], PetscInt numConstants,
           const PetscScalar constants[], PetscScalar g3[]);

} // namespace PetscFEPoroelasticity
} // namespace FSRM

#endif // PETSC_FE_POROELASTICITY_HPP
