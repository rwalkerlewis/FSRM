#ifndef FSRM_FAULTS_FEKERNELS_FAULT_COHESIVE_KIN_KERNELS_HPP
#define FSRM_FAULTS_FEKERNELS_FAULT_COHESIVE_KIN_KERNELS_HPP

#include <petsc.h>
#include <petscds.h>

namespace FSRM {
namespace faults {
namespace fekernels {

/**
 * @brief PyLith-style pointwise PetscFE kernels for a kinematic cohesive fault.
 *
 * Registered via PetscWeakFormAddBd{Residual,Jacobian} with three distinct
 * (label, value, part) keys so the hybrid-cell integrator dispatches each
 * kernel to the correct side of the cohesive cell:
 *
 *   part = 0  -> NEG side displacement residual / Jacobian
 *   part = 1  -> POS side displacement residual / Jacobian
 *   part = 2  -> FAULT (Lagrange multiplier) residual / Jacobian
 *
 * Scope (Phase 2 commit B): locked fault only. The slip field d(x,t) in
 * f0l is held at zero. KinSrc / prescribed slip lands in Phase 3.
 *
 * Sign convention (matches PyLith
 * libsrc/pylith/fekernels/FaultCohesiveKin.hh @ 6accf7e7):
 *
 *   f_0^{u^-}_i = +lambda_i
 *   f_0^{u^+}_i = -lambda_i
 *   f_0^{lambda}_i = (u^-_i - u^+_i) + d_i(x,t)
 *
 * so the converged solution satisfies u^+ - u^- = d. The hybrid closure
 * packs neg-side displacement at uOff[i_disp] + [0, dim), positive-side
 * at uOff[i_disp] + [dim, 2*dim), and Lagrange at uOff[i_lagrange].
 *
 * Jacobian blocks (all scalar, identity in spaceDim x spaceDim):
 *
 *   J_{u^-, lambda} = +I        (part = 0)
 *   J_{u^+, lambda} = -I        (part = 1)
 *   J_{lambda, u}   = [+I, -I]  (part = 2; neg block +I, pos block -I)
 *   J_{lambda, lambda} = 0      (not registered)
 */
class FaultCohesiveKinKernels
{
public:
    // -----------------------------------------------------------------
    // Residual kernels
    // -----------------------------------------------------------------

    /** f_0^{u^-} = +lambda; registered on (label, value, i_disp, part=0). */
    static void f0u_neg(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                        const PetscInt uOff[], const PetscInt uOff_x[],
                        const PetscScalar u[], const PetscScalar u_t[],
                        const PetscScalar u_x[],
                        const PetscInt aOff[], const PetscInt aOff_x[],
                        const PetscScalar a[], const PetscScalar a_t[],
                        const PetscScalar a_x[],
                        PetscReal t, const PetscReal x[], const PetscReal n[],
                        PetscInt numConstants, const PetscScalar constants[],
                        PetscScalar f[]);

    /** f_0^{u^+} = -lambda; registered on (label, value, i_disp, part=1). */
    static void f0u_pos(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                        const PetscInt uOff[], const PetscInt uOff_x[],
                        const PetscScalar u[], const PetscScalar u_t[],
                        const PetscScalar u_x[],
                        const PetscInt aOff[], const PetscInt aOff_x[],
                        const PetscScalar a[], const PetscScalar a_t[],
                        const PetscScalar a_x[],
                        PetscReal t, const PetscReal x[], const PetscReal n[],
                        PetscInt numConstants, const PetscScalar constants[],
                        PetscScalar f[]);

    /**
     * f_0^{lambda} = (u^- - u^+) + d(x,t); for the locked case d = 0.
     * Registered on (label, value, i_lagrange, part=2).
     */
    static void f0l(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[],
                    const PetscScalar u[], const PetscScalar u_t[],
                    const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[],
                    const PetscScalar a[], const PetscScalar a_t[],
                    const PetscScalar a_x[],
                    PetscReal t, const PetscReal x[], const PetscReal n[],
                    PetscInt numConstants, const PetscScalar constants[],
                    PetscScalar f[]);

    // -----------------------------------------------------------------
    // Jacobian g0 kernels. All other g*-slots are zero.
    // -----------------------------------------------------------------

    /** d(f_0^{u^-})/d(lambda) = +I; registered on part=0 with (i_disp, i_lagrange). */
    static void Jf0ul_neg(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                          const PetscInt uOff[], const PetscInt uOff_x[],
                          const PetscScalar u[], const PetscScalar u_t[],
                          const PetscScalar u_x[],
                          const PetscInt aOff[], const PetscInt aOff_x[],
                          const PetscScalar a[], const PetscScalar a_t[],
                          const PetscScalar a_x[],
                          PetscReal t, PetscReal u_tShift,
                          const PetscReal x[], const PetscReal n[],
                          PetscInt numConstants, const PetscScalar constants[],
                          PetscScalar g0[]);

    /** d(f_0^{u^+})/d(lambda) = -I; registered on part=1 with (i_disp, i_lagrange). */
    static void Jf0ul_pos(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                          const PetscInt uOff[], const PetscInt uOff_x[],
                          const PetscScalar u[], const PetscScalar u_t[],
                          const PetscScalar u_x[],
                          const PetscInt aOff[], const PetscInt aOff_x[],
                          const PetscScalar a[], const PetscScalar a_t[],
                          const PetscScalar a_x[],
                          PetscReal t, PetscReal u_tShift,
                          const PetscReal x[], const PetscReal n[],
                          PetscInt numConstants, const PetscScalar constants[],
                          PetscScalar g0[]);

    /**
     * d(f_0^{lambda})/d(u) writes both neg (+I) and pos (-I) blocks into
     * the 2*dim-wide trial space of the hybrid closure.
     * Registered on part=2 with (i_lagrange, i_disp).
     */
    static void Jf0lu(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                      const PetscInt uOff[], const PetscInt uOff_x[],
                      const PetscScalar u[], const PetscScalar u_t[],
                      const PetscScalar u_x[],
                      const PetscInt aOff[], const PetscInt aOff_x[],
                      const PetscScalar a[], const PetscScalar a_t[],
                      const PetscScalar a_x[],
                      PetscReal t, PetscReal u_tShift,
                      const PetscReal x[], const PetscReal n[],
                      PetscInt numConstants, const PetscScalar constants[],
                      PetscScalar g0[]);
};

}  // namespace fekernels
}  // namespace faults
}  // namespace FSRM

#endif  // FSRM_FAULTS_FEKERNELS_FAULT_COHESIVE_KIN_KERNELS_HPP
