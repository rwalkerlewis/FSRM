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
    static constexpr PetscInt COHESIVE_CONST_MODE  = 25;  // 0=locked, 1=slipping, 2=prescribed_slip
    static constexpr PetscInt COHESIVE_CONST_MU_F  = 26;  // friction coefficient (was 25)
    static constexpr PetscInt COHESIVE_CONST_TENSILE_STRENGTH = 27;  // tensile strength (Pa)
    static constexpr PetscInt COHESIVE_CONST_PRESCRIBED_SLIP_X = 28;  // prescribed slip x
    static constexpr PetscInt COHESIVE_CONST_PRESCRIBED_SLIP_Y = 29;  // prescribed slip y
    static constexpr PetscInt COHESIVE_CONST_PRESCRIBED_SLIP_Z = 30;  // prescribed slip z
    static constexpr PetscInt COHESIVE_CONST_COUNT = 31;  // minimum constant slots for base cohesive

    // Slip-weakening friction constants (slots 70-73, after viscoelastic 54-69)
    static constexpr PetscInt COHESIVE_CONST_FRICTION_MODEL = 70;  // 0=constant, 1=slip_weakening
    static constexpr PetscInt COHESIVE_CONST_MU_S = 71;            // static friction coefficient
    static constexpr PetscInt COHESIVE_CONST_MU_D = 72;            // dynamic friction coefficient
    static constexpr PetscInt COHESIVE_CONST_DC = 73;              // critical slip distance (m)

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
     * @brief Set the tensile strength for cohesive fracture opening
     * @param T_s Tensile strength in Pa (0 = no tensile failure check)
     */
    void setTensileStrength(double T_s);

    /**
     * @brief Set the prescribed slip vector in Cartesian coordinates
     * @param sx Prescribed displacement jump in x
     * @param sy Prescribed displacement jump in y
     * @param sz Prescribed displacement jump in z
     */
    void setPrescribedSlip(double sx, double sy, double sz);

    /**
     * @brief Read the configured prescribed slip vector.
     * @param[out] sx X component (Cartesian).
     * @param[out] sy Y component (Cartesian).
     * @param[out] sz Z component (Cartesian).
     */
    void getPrescribedSlip(double& sx, double& sy, double& sz) const {
        sx = prescribed_slip_[0];
        sy = prescribed_slip_[1];
        sz = prescribed_slip_[2];
    }

    /**
     * @brief Set slip-weakening friction parameters
     * @param mu_s Static friction coefficient
     * @param mu_d Dynamic friction coefficient
     * @param D_c Critical slip distance (m)
     */
    void setSlipWeakeningParams(double mu_s, double mu_d, double D_c);

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
     * @brief Register prescribed-slip residual and Jacobian callbacks with PetscDS
     *
     * Similar to registerWithDS() but registers f0_prescribed_slip instead of
     * f0_lagrange_constraint, enforcing u+ - u- = delta_prescribed.
     * The prescribed slip vector is stored in constants slots 28-30.
     *
     * @param prob The PetscDS problem
     * @param displacement_field Field index for displacement
     * @param lagrange_field Field index for Lagrange multiplier
     * @return PETSc error code
     */
    PetscErrorCode registerPrescribedSlipWithDS(PetscDS prob,
                                                 int displacement_field,
                                                 int lagrange_field);

    /**
     * @brief Update auxiliary state data (friction state) before residual evaluation
     */
    PetscErrorCode updateAuxiliaryState(DM dm, FaultCohesiveDyn* fault);

    bool isLocked() const { return locked_; }
    bool isPrescribedSlip() const { return prescribed_slip_mode_; }

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
     * @brief Residual for Lagrange multiplier: prescribed slip constraint
     *
     * Enforces u+ - u- = delta_prescribed:
     *   f_lambda = (u+ - u-) - delta
     * where delta is read from constants[COHESIVE_CONST_PRESCRIBED_SLIP_X..Z].
     */
    static void f0_prescribed_slip(PetscInt dim, PetscInt Nf, PetscInt NfAux,
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
    bool prescribed_slip_mode_ = false;
    double mu_friction_ = 0.6;  // Coulomb friction coefficient for slipping mode
    double tensile_strength_ = 0.0;  // Tensile strength in Pa (0 = no check)
    double prescribed_slip_[3] = {0.0, 0.0, 0.0};  // Cartesian slip vector
    double mu_static_ = 0.677;           // Static friction coefficient (slip-weakening)
    double mu_dynamic_ = 0.525;          // Dynamic friction coefficient (slip-weakening)
    double critical_slip_distance_ = 0.40;  // Critical slip distance in m
    bool use_slip_weakening_ = false;    // True if slip-weakening friction model is active
    FaultFrictionModel* friction_model_ = nullptr;
};

} // namespace FSRM

#endif // COHESIVE_FAULT_KERNEL_HPP
