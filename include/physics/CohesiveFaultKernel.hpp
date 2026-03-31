#ifndef COHESIVE_FAULT_KERNEL_HPP
#define COHESIVE_FAULT_KERNEL_HPP

/**
 * @file CohesiveFaultKernel.hpp
 * @brief PetscDS-compatible pointwise callbacks for cohesive fault cells
 *
 * Provides residual and Jacobian callbacks for the Lagrange multiplier
 * constraint on cohesive cells, supporting both:
 *   - LOCKED mode: zero displacement jump (u+ - u- = 0)
 *   - SLIPPING mode: friction-governed traction constraint
 *
 * These callbacks are registered with PetscDS for hybrid/cohesive cells
 * and are evaluated by PETSc's DMPlexSNESComputeResidualHybrid.
 *
 * References:
 *   Aagaard et al. (2013): Domain decomposition for fault slip in FEM
 *   PETSc source: dm/impls/plex/plexsnes.c (hybrid cell integration)
 */

#include <petscdmplex.h>
#include <petscds.h>

namespace FSRM {

// Forward declarations
class FaultFrictionModel;
class FaultCohesiveDyn;

/**
 * @brief Manages PetscDS callback registration for cohesive fault cells
 *
 * Constants layout (shared PetscDS constant array):
 *   Slots 0-24 are reserved for fluid/elasticity/poroelasticity (unified layout).
 *   Slot COHESIVE_CONST_MODE (25): 0.0 = locked, 1.0 = slipping.
 *   Slot COHESIVE_CONST_MU_F  (26): Coulomb friction coefficient μ_f.
 *
 * registerWithDS() reads the existing constants array, patches these two
 * dedicated slots, and writes the updated array back with PetscDSSetConstants.
 * The static callbacks then read from these fixed indices, so the mode is
 * always consistent with the last call to setMode() / setFrictionCoefficient().
 */
class CohesiveFaultKernel {
public:
    // Dedicated constant indices for the cohesive kernel.
    // Unified constants: [0-24] = fluid/elasticity/poroelasticity, [25-26] = cohesive fault
    static constexpr PetscInt COHESIVE_CONST_MODE  = 25;  // 0=locked, 1=slipping (was 24)
    static constexpr PetscInt COHESIVE_CONST_MU_F  = 26;  // friction coefficient (was 25)
    static constexpr PetscInt COHESIVE_CONST_COUNT = 27;  // total constant slots needed (was 26)

    CohesiveFaultKernel();

    /**
     * @brief Set the fault constraint mode
     * @param is_locked If true, enforce zero displacement jump (locked fault)
     *                  If false, enforce friction-governed traction
     *
     * NOTE: Call registerWithDS() after setMode() to propagate the change to
     * the PetscDS constants read by the static callbacks.
     */
    void setMode(bool is_locked);

    /**
     * @brief Set the Coulomb friction coefficient used in slipping mode
     * @param mu_f Friction coefficient (dimensionless, typically 0.4–0.85)
     *
     * NOTE: Call registerWithDS() after setFrictionCoefficient() to propagate
     * the change to the PetscDS constants read by the static callbacks.
     */
    void setFrictionCoefficient(double mu_f);

    /**
     * @brief Set the friction model for slipping constraint
     */
    void setFrictionModel(FaultFrictionModel* model);

    /**
     * @brief Register residual and Jacobian callbacks with PetscDS and
     *        update the PetscDS constant slots reserved for the cohesive kernel.
     *
     * Reads the existing constants from @p prob, patches slots
     * COHESIVE_CONST_MODE and COHESIVE_CONST_MU_F with the current mode and
     * friction coefficient, then calls PetscDSSetConstants to write them back.
     * This ensures the static callbacks observe the correct mode on the next
     * residual/Jacobian evaluation.
     *
     * @param prob The PetscDS problem
     * @param displacement_field Field index for displacement
     * @param lagrange_field Field index for Lagrange multiplier
     * @return PETSc error code
     */
    PetscErrorCode registerWithDS(PetscDS prob,
                                   int displacement_field,
                                   int lagrange_field);

    /**
     * @brief Update auxiliary state data (friction state) before residual evaluation
     */
    PetscErrorCode updateAuxiliaryState(DM dm, FaultCohesiveDyn* fault);

    bool isLocked() const { return locked_; }

    // =========================================================================
    // Static PetscDS pointwise callback functions
    // =========================================================================

    /**
     * @brief Residual for displacement equation on cohesive cells
     *
     * Applies the Lagrange multiplier traction to the + and - sides:
     *   f_u(+side) = +lambda * area
     *   f_u(-side) = -lambda * area
     */
    static void f0_displacement_cohesive(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                                          const PetscInt uOff[], const PetscInt uOff_x[],
                                          const PetscScalar u[], const PetscScalar u_t[],
                                          const PetscScalar u_x[],
                                          const PetscInt aOff[], const PetscInt aOff_x[],
                                          const PetscScalar a[], const PetscScalar a_t[],
                                          const PetscScalar a_x[],
                                          PetscReal t, const PetscReal x[],
                                          const PetscReal n[],
                                          PetscInt numConstants,
                                          const PetscScalar constants[],
                                          PetscScalar f[]);

    /**
     * @brief Residual for Lagrange multiplier (constraint equation)
     *
     * LOCKED: f_lambda = u+ - u- (displacement jump = 0)
     * SLIPPING: f_lambda = |tau| - f(V,theta)*sigma_n_eff (friction law)
     *
     * The mode is selected via constants[COHESIVE_CONST_MODE] (0 = locked, 1 = slipping).
     * The friction coefficient is read from constants[COHESIVE_CONST_MU_F].
     */
    static void f0_lagrange_constraint(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                                        const PetscInt uOff[], const PetscInt uOff_x[],
                                        const PetscScalar u[], const PetscScalar u_t[],
                                        const PetscScalar u_x[],
                                        const PetscInt aOff[], const PetscInt aOff_x[],
                                        const PetscScalar a[], const PetscScalar a_t[],
                                        const PetscScalar a_x[],
                                        PetscReal t, const PetscReal x[],
                                        const PetscReal n[],
                                        PetscInt numConstants,
                                        const PetscScalar constants[],
                                        PetscScalar f[]);

    /**
     * @brief Jacobian: d(f_u)/d(lambda) = +/- I
     */
    static void g0_displacement_lagrange(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                                          const PetscInt uOff[], const PetscInt uOff_x[],
                                          const PetscScalar u[], const PetscScalar u_t[],
                                          const PetscScalar u_x[],
                                          const PetscInt aOff[], const PetscInt aOff_x[],
                                          const PetscScalar a[], const PetscScalar a_t[],
                                          const PetscScalar a_x[],
                                          PetscReal t, PetscReal u_tShift,
                                          const PetscReal x[], const PetscReal n[],
                                          PetscInt numConstants,
                                          const PetscScalar constants[],
                                          PetscScalar g0[]);

    /**
     * @brief Jacobian: d(f_lambda)/d(u) = I (locked) or d(friction)/d(u) (slipping)
     */
    static void g0_lagrange_displacement(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                                          const PetscInt uOff[], const PetscInt uOff_x[],
                                          const PetscScalar u[], const PetscScalar u_t[],
                                          const PetscScalar u_x[],
                                          const PetscInt aOff[], const PetscInt aOff_x[],
                                          const PetscScalar a[], const PetscScalar a_t[],
                                          const PetscScalar a_x[],
                                          PetscReal t, PetscReal u_tShift,
                                          const PetscReal x[], const PetscReal n[],
                                          PetscInt numConstants,
                                          const PetscScalar constants[],
                                          PetscScalar g0[]);

    /**
     * @brief Jacobian: d(f_lambda)/d(lambda) = 0 (locked) or d(friction)/d(lambda) (slipping)
     */
    static void g0_lagrange_lagrange(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                                      const PetscInt uOff[], const PetscInt uOff_x[],
                                      const PetscScalar u[], const PetscScalar u_t[],
                                      const PetscScalar u_x[],
                                      const PetscInt aOff[], const PetscInt aOff_x[],
                                      const PetscScalar a[], const PetscScalar a_t[],
                                      const PetscScalar a_x[],
                                      PetscReal t, PetscReal u_tShift,
                                      const PetscReal x[], const PetscReal n[],
                                      PetscInt numConstants,
                                      const PetscScalar constants[],
                                      PetscScalar g0[]);

private:
    bool locked_ = true;
    double mu_friction_ = 0.6;  // Coulomb friction coefficient for slipping mode
    FaultFrictionModel* friction_model_ = nullptr;
};

} // namespace FSRM

#endif // COHESIVE_FAULT_KERNEL_HPP
