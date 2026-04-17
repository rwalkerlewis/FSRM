#ifndef FSRM_FAULTS_FAULT_COHESIVE_KIN_HPP
#define FSRM_FAULTS_FAULT_COHESIVE_KIN_HPP

#include "faults/FaultCohesive.hpp"

#include <petsc.h>
#include <petscdm.h>
#include <petscds.h>

namespace FSRM {
namespace faults {

/**
 * @brief Kinematic cohesive fault (PyLith FaultCohesiveKin analog).
 *
 * Phase 2 commit B2 scope: locked fault only.
 *
 * Wiring sequence:
 *   1. Simulator creates FaultCohesiveKin.
 *   2. Simulator creates disp + Lagrange PetscFE and records field indices.
 *   3. Before DMCreateDS: restrictLagrangeToInterface() creates the
 *      "cohesive interface" DMLabel and calls DMSetField to confine the
 *      Lagrange field to cohesive cells.
 *   4. After DMCreateDS: setKernels() gets the cohesive-cell DS via
 *      DMGetCellDS(cMax) and registers kernels using
 *      PetscWeakFormSetIndexBdResidual / PetscWeakFormSetIndexBdJacobian
 *      (ex5.c pattern). It also builds the PetscFormKey[3] array and an
 *      IS over cohesive cells, which computeHybridResidual /
 *      computeHybridJacobian use to call DMPlexComputeResidualHybridByKey /
 *      DMPlexComputeJacobianHybridByKey explicitly from FormFunction /
 *      FormJacobian.
 *
 * Mirrors PyLith libsrc/pylith/faults/FaultCohesiveKin.{hh,cc}
 * @ geodynamics/pylith 6accf7e7395aa792444ad1b489d455e41ece34bc.
 */
class FaultCohesiveKin : public FaultCohesive
{
public:
    FaultCohesiveKin() = default;
    ~FaultCohesiveKin() override;

    PetscErrorCode initialize(DM dm, const FaultConfig& config) override;

    PetscErrorCode restrictLagrangeToInterface(DM dm, DMLabel fault_label,
                                               PetscInt i_lagrange,
                                               PetscObject fe_lagrange);

    /**
     * @brief Register locked-fault kernels on the cohesive-cell DS weak form.
     *
     * Follows the ex5.c (PETSc) pattern: DMGetCellDS(dm, cMax) to get the
     * region DS, PetscWeakFormSetIndexBdResidual/BdJacobian to register
     * kernels keyed by (cohesive_interface_label_, value={1,2,3}), and
     * ISCreateStride for the cohesive cell range.
     *
     * Keys:
     *   keys_[0] = {label, 1, i_disp,     0}  → f0u_neg (neg side)
     *   keys_[1] = {label, 2, i_disp,     0}  → f0u_pos (pos side)
     *   keys_[2] = {label, 3, i_lagrange, 0}  → f0l     (constraint)
     */
    PetscErrorCode setKernels(DM dm, PetscInt i_disp, PetscInt i_lagrange);

    /** @return true after setKernels() succeeds. */
    bool kernelsSet() const { return kernels_set_; }

    /**
     * @brief Compute cohesive-cell residual contributions.
     *
     * Calls DMPlexComputeResidualHybridByKey with the stored keys_ and
     * cohesive_cells_is_. Must be called on local (scatter-scatter) vectors.
     */
    PetscErrorCode computeHybridResidual(DM dm, PetscReal t,
                                         Vec locU, Vec locU_t,
                                         PetscReal s, Vec locF,
                                         void* ctx);

    /**
     * @brief Compute cohesive-cell Jacobian contributions.
     *
     * Calls DMPlexComputeJacobianHybridByKey with the stored keys_ and
     * cohesive_cells_is_.
     */
    PetscErrorCode computeHybridJacobian(DM dm, PetscReal t, PetscReal a,
                                          Vec locU, Vec locU_t,
                                          Mat J, Mat P,
                                          void* ctx);

private:
    PetscFormKey keys_[3]       = {};
    IS           cohesive_cells_is_ = nullptr;
    bool         kernels_set_   = false;
};

}  // namespace faults
}  // namespace FSRM

#endif  // FSRM_FAULTS_FAULT_COHESIVE_KIN_HPP
