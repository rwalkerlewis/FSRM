/**
 * @file FaultCohesiveKinKernels.cpp
 * @brief Locked-fault pointwise kernels for Phase 2 commit B.
 *
 * Scope: locked fault only, d(x,t) = 0. Prescribed slip (KinSrc) and
 * friction-governed faulting land in later phases. The file contains
 * only the kernel math; registration happens in FaultCohesiveKin.
 */

#include "faults/fekernels/FaultCohesiveKinKernels.hpp"

namespace FSRM {
namespace faults {
namespace fekernels {

// The disp field is index 0 in the uOff array passed to these kernels;
// the Lagrange field is index 1. This ordering matches how
// FaultCohesiveKin::_setKernels registers the field-key pair (see the
// PetscWeakFormAddBd* calls: fieldI/fieldJ are the whole-DM indices,
// but uOff[] is indexed by the DS field enumeration which PETSc
// renumbers to [i_disp, i_lagrange] = [0, 1] on the cohesive region DS).
//
// We do not hardcode field indices below; each kernel reads the correct
// offset through uOff. The region DS constructed by DMSetField with a
// label presents exactly two active fields on cohesive cells: disp and
// lagrange, in that order.
static constexpr PetscInt DISP_FIELD = 0;
static constexpr PetscInt LAGRANGE_FIELD = 1;

// ---------------------------------------------------------------------
// Residual kernels
// ---------------------------------------------------------------------

void FaultCohesiveKinKernels::f0u_neg(
    PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[],
    const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[],
    const PetscScalar a[], const PetscScalar a_t[],
    const PetscScalar a_x[],
    PetscReal t, const PetscReal x[], const PetscReal n[],
    PetscInt numConstants, const PetscScalar constants[],
    PetscScalar f[])
{
    (void)Nf; (void)NfAux; (void)uOff_x; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)x; (void)n; (void)numConstants; (void)constants;

    // f_0^{u^-}_i = + lambda_i
    for (PetscInt d = 0; d < dim; ++d) {
        f[d] = +u[uOff[LAGRANGE_FIELD] + d];
    }
}

void FaultCohesiveKinKernels::f0u_pos(
    PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[],
    const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[],
    const PetscScalar a[], const PetscScalar a_t[],
    const PetscScalar a_x[],
    PetscReal t, const PetscReal x[], const PetscReal n[],
    PetscInt numConstants, const PetscScalar constants[],
    PetscScalar f[])
{
    (void)Nf; (void)NfAux; (void)uOff_x; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)x; (void)n; (void)numConstants; (void)constants;

    // f_0^{u^+}_i = - lambda_i
    for (PetscInt d = 0; d < dim; ++d) {
        f[d] = -u[uOff[LAGRANGE_FIELD] + d];
    }
}

void FaultCohesiveKinKernels::f0l(
    PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[],
    const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[],
    const PetscScalar a[], const PetscScalar a_t[],
    const PetscScalar a_x[],
    PetscReal t, const PetscReal x[], const PetscReal n[],
    PetscInt numConstants, const PetscScalar constants[],
    PetscScalar f[])
{
    (void)Nf; (void)NfAux; (void)uOff_x; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)x; (void)n; (void)numConstants; (void)constants;

    // Locked-fault scope: d(x,t) = 0.
    // f_0^{lambda}_i = (u^-_i - u^+_i) + 0
    const PetscInt uOff_disp = uOff[DISP_FIELD];
    for (PetscInt d = 0; d < dim; ++d) {
        const PetscScalar u_neg = u[uOff_disp + d];
        const PetscScalar u_pos = u[uOff_disp + dim + d];
        f[d] = u_neg - u_pos;
    }
}

// ---------------------------------------------------------------------
// Jacobian g0 kernels
// ---------------------------------------------------------------------

void FaultCohesiveKinKernels::Jf0ul_neg(
    PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[],
    const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[],
    const PetscScalar a[], const PetscScalar a_t[],
    const PetscScalar a_x[],
    PetscReal t, PetscReal u_tShift,
    const PetscReal x[], const PetscReal n[],
    PetscInt numConstants, const PetscScalar constants[],
    PetscScalar g0[])
{
    (void)Nf; (void)NfAux; (void)uOff; (void)uOff_x;
    (void)u; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)u_tShift; (void)x; (void)n;
    (void)numConstants; (void)constants;

    // d(+lambda_i)/d(lambda_j) = +delta_{ij}
    // Layout: g0[Nc_test * Nc_trial] = g0[dim * dim]
    for (PetscInt i = 0; i < dim; ++i) {
        for (PetscInt j = 0; j < dim; ++j) {
            g0[i * dim + j] = (i == j) ? +1.0 : 0.0;
        }
    }
}

void FaultCohesiveKinKernels::Jf0ul_pos(
    PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[],
    const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[],
    const PetscScalar a[], const PetscScalar a_t[],
    const PetscScalar a_x[],
    PetscReal t, PetscReal u_tShift,
    const PetscReal x[], const PetscReal n[],
    PetscInt numConstants, const PetscScalar constants[],
    PetscScalar g0[])
{
    (void)Nf; (void)NfAux; (void)uOff; (void)uOff_x;
    (void)u; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)u_tShift; (void)x; (void)n;
    (void)numConstants; (void)constants;

    // d(-lambda_i)/d(lambda_j) = -delta_{ij}
    for (PetscInt i = 0; i < dim; ++i) {
        for (PetscInt j = 0; j < dim; ++j) {
            g0[i * dim + j] = (i == j) ? -1.0 : 0.0;
        }
    }
}

void FaultCohesiveKinKernels::Jf0lu(
    PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[],
    const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[],
    const PetscScalar a[], const PetscScalar a_t[],
    const PetscScalar a_x[],
    PetscReal t, PetscReal u_tShift,
    const PetscReal x[], const PetscReal n[],
    PetscInt numConstants, const PetscScalar constants[],
    PetscScalar g0[])
{
    (void)Nf; (void)NfAux; (void)uOff; (void)uOff_x;
    (void)u; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)u_tShift; (void)x; (void)n;
    (void)numConstants; (void)constants;

    // f_0^{lambda}_i = u^-_i - u^+_i
    // Trial space is the hybrid displacement closure with 2*dim components:
    // [neg_0..neg_{dim-1}, pos_0..pos_{dim-1}].
    // g0 layout: g0[i * (2*dim) + j]
    //   j in [0, dim)         -> d(f_l_i)/d(u_neg_j) = +delta_{ij}
    //   j in [dim, 2*dim)     -> d(f_l_i)/d(u_pos_{j-dim}) = -delta_{ij-dim}
    const PetscInt trial = 2 * dim;
    for (PetscInt i = 0; i < dim; ++i) {
        for (PetscInt j = 0; j < dim; ++j) {
            g0[i * trial + j]       = (i == j) ? +1.0 : 0.0;  // neg block
            g0[i * trial + dim + j] = (i == j) ? -1.0 : 0.0;  // pos block
        }
    }
}

}  // namespace fekernels
}  // namespace faults
}  // namespace FSRM
