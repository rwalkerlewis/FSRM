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
 * Phase 2 commit B scope: locked fault only. The initialize() override
 * and wiring into Simulator land in commit B2; this header exists in
 * commit B1 so the class and its static kernel registrations compile
 * standalone without affecting the active assembly path.
 *
 * Intended wiring (activated in B2):
 *   1. Simulator creates FaultCohesiveKin.
 *   2. Simulator creates the disp and lagrange PetscFE objects and
 *      records their field indices i_disp, i_lagrange.
 *   3. Before DMCreateDS, Simulator calls
 *        restrictLagrangeToInterface(dm, i_lagrange, fe_lagrange)
 *      which (a) calls createCohesiveInterfaceLabel from the base and
 *      (b) registers the Lagrange field on that label via DMSetField.
 *   4. After DMCreateDS, Simulator calls setKernels(dm, i_disp,
 *      i_lagrange) which registers the PyLith three-sided pointwise
 *      kernels via PetscWeakFormAddBdResidual / PetscWeakFormAddBdJacobian.
 *
 * Mirrors PyLith libsrc/pylith/faults/FaultCohesiveKin.{hh,cc}
 * @ geodynamics/pylith 6accf7e7395aa792444ad1b489d455e41ece34bc.
 */
class FaultCohesiveKin : public FaultCohesive
{
public:
    FaultCohesiveKin() = default;
    ~FaultCohesiveKin() override = default;

    /**
     * @brief Create the cohesive interface DMLabel and restrict the
     *        Lagrange field to live only on cohesive cells.
     *
     * Must be called BEFORE DMCreateDS. Wraps
     *   FaultCohesive::createCohesiveInterfaceLabel
     *   DMSetField(dm, i_lagrange, cohesive_interface_label, fe_lagrange)
     *
     * @param dm            DMPlex after DMPlexConstructCohesiveCells.
     * @param fault_label   DMLabel used to build the cohesive interface
     *                      label. Pass NULL to accept the default.
     * @param i_lagrange    Field index of the Lagrange multiplier.
     * @param fe_lagrange   PetscFE object for the Lagrange field. Passed
     *                      as PetscObject to match the DMSetField API.
     */
    PetscErrorCode restrictLagrangeToInterface(DM dm, DMLabel fault_label,
                                               PetscInt i_lagrange,
                                               PetscObject fe_lagrange);

    /**
     * @brief Register the locked-fault pointwise kernels on the DS with
     *        the three-sided (part=0, part=1, part=2) key pattern.
     *
     * Must be called AFTER DMCreateDS, so the PetscWeakForm owned by
     * the DS is available. The cohesive interface label must already
     * be created (via @ref restrictLagrangeToInterface or the base
     * @ref FaultCohesive::createCohesiveInterfaceLabel).
     *
     * @param dm            DMPlex with the DS already created.
     * @param i_disp        Field index of displacement.
     * @param i_lagrange    Field index of Lagrange multiplier.
     */
    PetscErrorCode setKernels(DM dm, PetscInt i_disp, PetscInt i_lagrange);
};

}  // namespace faults
}  // namespace FSRM

#endif  // FSRM_FAULTS_FAULT_COHESIVE_KIN_HPP
