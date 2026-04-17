/**
 * @file FaultCohesiveKin.cpp
 * @brief Kinematic cohesive fault class (PyLith FaultCohesiveKin analog).
 *
 * Phase 2 commit B1 scope: standalone class + kernel registration helpers.
 * Nothing in Simulator calls this yet; wiring lands in B2.
 */

#include "faults/FaultCohesiveKin.hpp"
#include "faults/fekernels/FaultCohesiveKinKernels.hpp"

namespace FSRM {
namespace faults {

PetscErrorCode FaultCohesiveKin::restrictLagrangeToInterface(
    DM dm, DMLabel fault_label, PetscInt i_lagrange, PetscObject fe_lagrange)
{
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

    // Create (or fetch) the cohesive interface label on the DM via the
    // base-class helper; this caches the DMLabel in
    // cohesive_interface_label_.
    ierr = createCohesiveInterfaceLabel(dm, fault_label); CHKERRQ(ierr);

    // Restrict the Lagrange multiplier field to live only on cells
    // tagged by the cohesive interface label. PyLith applies this right
    // before DMCreateDS so the region-DS for the label presents the
    // Lagrange field as active, and the bulk-DS does not.
    ierr = DMSetField(dm, i_lagrange, cohesive_interface_label_, fe_lagrange);
    CHKERRQ(ierr);

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode FaultCohesiveKin::setKernels(DM dm, PetscInt i_disp,
                                            PetscInt i_lagrange)
{
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

    PetscDS ds = nullptr;
    ierr = DMGetDS(dm, &ds); CHKERRQ(ierr);

    PetscWeakForm wf = nullptr;
    ierr = PetscDSGetWeakForm(ds, &wf); CHKERRQ(ierr);

    DMLabel label = cohesive_interface_label_;
    if (label == nullptr) {
        ierr = DMGetLabel(dm, "cohesive interface", &label); CHKERRQ(ierr);
    }
    PetscCheck(label != nullptr, PetscObjectComm((PetscObject)dm),
               PETSC_ERR_ARG_WRONGSTATE,
               "FaultCohesiveKin::setKernels requires the 'cohesive interface' "
               "DMLabel to exist before kernel registration");

    const PetscInt value      = 1;
    const PetscInt part_neg   = 0;
    const PetscInt part_pos   = 1;
    const PetscInt part_fault = 2;

    using K = fekernels::FaultCohesiveKinKernels;

    // Residual: f0 only (no f1 on cohesive face in the quasistatic case).
    ierr = PetscWeakFormAddBdResidual(wf, label, value, i_disp, part_neg,
                                      K::f0u_neg, nullptr); CHKERRQ(ierr);
    ierr = PetscWeakFormAddBdResidual(wf, label, value, i_disp, part_pos,
                                      K::f0u_pos, nullptr); CHKERRQ(ierr);
    ierr = PetscWeakFormAddBdResidual(wf, label, value, i_lagrange, part_fault,
                                      K::f0l, nullptr); CHKERRQ(ierr);

    // Jacobian: g0 only. g1, g2, g3 are all null.
    ierr = PetscWeakFormAddBdJacobian(wf, label, value,
                                      i_disp, i_lagrange, part_neg,
                                      K::Jf0ul_neg, nullptr, nullptr, nullptr);
    CHKERRQ(ierr);
    ierr = PetscWeakFormAddBdJacobian(wf, label, value,
                                      i_disp, i_lagrange, part_pos,
                                      K::Jf0ul_pos, nullptr, nullptr, nullptr);
    CHKERRQ(ierr);
    ierr = PetscWeakFormAddBdJacobian(wf, label, value,
                                      i_lagrange, i_disp, part_fault,
                                      K::Jf0lu, nullptr, nullptr, nullptr);
    CHKERRQ(ierr);

    // J_{lambda,lambda} is exactly zero by design; do not register it.

    PetscFunctionReturn(PETSC_SUCCESS);
}

}  // namespace faults
}  // namespace FSRM
