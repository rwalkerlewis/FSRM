/**
 * @file FaultCohesiveKin.cpp
 * @brief Kinematic cohesive fault class (PyLith FaultCohesiveKin analog).
 *
 * Phase 2 commit B2: wired into Simulator. Kernels are registered on the
 * cohesive-cell DS weak form using the ex5.c PETSc pattern
 * (DMGetCellDS + PetscWeakFormSetIndexBdResidual/BdJacobian). Residual and
 * Jacobian are integrated explicitly in FormFunction/FormJacobian via
 * DMPlexComputeResidualHybridByKey / DMPlexComputeJacobianHybridByKey;
 * DMPlexTSComputeIFunctionFEM does NOT automatically dispatch to the hybrid
 * path.
 */

#include "faults/FaultCohesiveKin.hpp"
#include "faults/fekernels/FaultCohesiveKinKernels.hpp"

#include <petscdmplex.h>

namespace FSRM {
namespace faults {

FaultCohesiveKin::~FaultCohesiveKin()
{
    if (cohesive_cells_is_) {
        ISDestroy(&cohesive_cells_is_);
    }
}

PetscErrorCode FaultCohesiveKin::initialize(DM dm, const FaultConfig& /*config*/)
{
    PetscFunctionBeginUser;
    PetscErrorCode ierr = createCohesiveInterfaceLabel(dm, nullptr); CHKERRQ(ierr);
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode FaultCohesiveKin::restrictLagrangeToInterface(
    DM dm, DMLabel fault_label, PetscInt i_lagrange, PetscObject fe_lagrange)
{
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

    ierr = createCohesiveInterfaceLabel(dm, fault_label); CHKERRQ(ierr);

    ierr = DMSetField(dm, i_lagrange, cohesive_interface_label_, fe_lagrange);
    CHKERRQ(ierr);

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode FaultCohesiveKin::setKernels(DM dm, PetscInt i_disp,
                                            PetscInt i_lagrange)
{
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

    // Locate cohesive cell range [cMax, cEnd).
    PetscInt cMax = 0, cEnd = 0;
    ierr = DMPlexGetSimplexOrBoxCells(dm, 0, nullptr, &cMax); CHKERRQ(ierr);
    ierr = DMPlexGetHeightStratum(dm, 0, nullptr, &cEnd); CHKERRQ(ierr);
    PetscCheck(cMax < cEnd,
               PetscObjectComm((PetscObject)dm), PETSC_ERR_ARG_WRONGSTATE,
               "FaultCohesiveKin::setKernels: no cohesive cells (cMax=%d cEnd=%d)",
               (int)cMax, (int)cEnd);

    // Get the DS for the first cohesive cell (region DS for the cohesive
    // interface label). This is the ex5.c pattern: DMGetCellDS(dm, cMax, &ds).
    PetscDS ds = nullptr;
    ierr = DMGetCellDS(dm, cMax, &ds, nullptr); CHKERRQ(ierr);
    PetscCheck(ds != nullptr,
               PetscObjectComm((PetscObject)dm), PETSC_ERR_ARG_WRONGSTATE,
               "FaultCohesiveKin::setKernels: DMGetCellDS returned null DS for "
               "first cohesive cell %d", (int)cMax);

    PetscWeakForm wf = nullptr;
    ierr = PetscDSGetWeakForm(ds, &wf); CHKERRQ(ierr);

    // The cohesive interface label (value=1 on all cohesive points) is reused
    // as the weak-form dictionary key, with three DISTINCT values to
    // differentiate the three face roles:
    //   value=1: neg-side face  (f0u_neg, Jf0ul_neg)
    //   value=2: pos-side face  (f0u_pos, Jf0ul_pos)
    //   value=3: fault interface (f0l,   Jf0lu)
    // These values do NOT need to be present on any mesh points — they are
    // purely a lookup key in the PetscWeakForm dictionary. (ex5.c uses the
    // "material" label, which has values only on bulk cells, for the same
    // purpose.)
    DMLabel label = cohesive_interface_label_;
    if (!label) {
        ierr = DMGetLabel(dm, "cohesive interface", &label); CHKERRQ(ierr);
    }
    PetscCheck(label != nullptr,
               PetscObjectComm((PetscObject)dm), PETSC_ERR_ARG_WRONGSTATE,
               "FaultCohesiveKin::setKernels: 'cohesive interface' DMLabel not found");

    using K = fekernels::FaultCohesiveKinKernels;

    // Residual kernels — PetscWeakFormSetIndexBdResidual signature:
    //   (wf, label, value, field, part, i0, f0, n1, f1)
    ierr = PetscWeakFormSetIndexBdResidual(wf, label, 1, i_disp,     0, 0,
                                           K::f0u_neg, 0, nullptr); CHKERRQ(ierr);
    ierr = PetscWeakFormSetIndexBdResidual(wf, label, 2, i_disp,     0, 0,
                                           K::f0u_pos, 0, nullptr); CHKERRQ(ierr);
    ierr = PetscWeakFormSetIndexBdResidual(wf, label, 3, i_lagrange, 0, 0,
                                           K::f0l,     0, nullptr); CHKERRQ(ierr);

    // Jacobian kernels — PetscWeakFormSetIndexBdJacobian signature:
    //   (wf, label, value, fieldI, fieldJ, part, i0, g0, n1, g1, n2, g2, n3, g3)
    ierr = PetscWeakFormSetIndexBdJacobian(wf, label, 1,
                                           i_disp, i_lagrange, 0, 0,
                                           K::Jf0ul_neg,
                                           0, nullptr, 0, nullptr, 0, nullptr);
    CHKERRQ(ierr);
    ierr = PetscWeakFormSetIndexBdJacobian(wf, label, 2,
                                           i_disp, i_lagrange, 0, 0,
                                           K::Jf0ul_pos,
                                           0, nullptr, 0, nullptr, 0, nullptr);
    CHKERRQ(ierr);
    ierr = PetscWeakFormSetIndexBdJacobian(wf, label, 3,
                                           i_lagrange, i_disp, 0, 0,
                                           K::Jf0lu,
                                           0, nullptr, 0, nullptr, 0, nullptr);
    CHKERRQ(ierr);

    // Build PetscFormKey array. keys_[i] is applied to cone[i] of each
    // cohesive cell by DMPlexComputeResidualHybridByKey:
    //   cone[0] = neg-side face  → keys_[0]
    //   cone[1] = pos-side face  → keys_[1]
    //   cone[2] = fault interface → keys_[2]
    keys_[0] = {label, 1, i_disp,     0};
    keys_[1] = {label, 2, i_disp,     0};
    keys_[2] = {label, 3, i_lagrange, 0};

    // Build cohesive cells IS [cMax, cEnd).
    if (cohesive_cells_is_) {
        ierr = ISDestroy(&cohesive_cells_is_); CHKERRQ(ierr);
    }
    ierr = ISCreateStride(PETSC_COMM_SELF, cEnd - cMax, cMax, 1,
                          &cohesive_cells_is_); CHKERRQ(ierr);

    // Only enable hybrid integration if the DS has the Lagrange field marked
    // cohesive. In FSRM the Lagrange field is currently added via
    // DMAddField(dm, nullptr, fe) + DMSetField(dm, f, label, fe), which does
    // not trigger PetscDS auto-cohesive marking in DMCreateDS. Without the
    // cohesive flag, DMPlexComputeResidualHybridByKey segfaults because the
    // field jets code paths assume cohesive layout.
    PetscBool is_cohesive = PETSC_FALSE;
    ierr = PetscDSGetCohesive(ds, i_lagrange, &is_cohesive); CHKERRQ(ierr);
    if (!is_cohesive) {
        PetscPrintf(PetscObjectComm((PetscObject)dm),
                    "FaultCohesiveKin: Lagrange field %d NOT marked cohesive in "
                    "region DS (DMAddField was called with nullptr label). "
                    "Hybrid integration disabled; relying on manual cohesive "
                    "assembly path.\n", (int)i_lagrange);
        kernels_set_ = false;
        PetscFunctionReturn(PETSC_SUCCESS);
    }

    kernels_set_ = true;

    if (PetscObjectComm((PetscObject)dm) != MPI_COMM_NULL) {
        PetscPrintf(PetscObjectComm((PetscObject)dm),
                    "FaultCohesiveKin: registered three-sided kernels on fields "
                    "%d (disp) and %d (lagrange), %d cohesive cells\n",
                    (int)i_disp, (int)i_lagrange, (int)(cEnd - cMax));
    }

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode FaultCohesiveKin::computeHybridResidual(DM dm, PetscReal t,
                                                        Vec locU, Vec locU_t,
                                                        PetscReal s, Vec locF,
                                                        void* ctx)
{
    PetscFunctionBeginUser;
    if (!kernels_set_) PetscFunctionReturn(PETSC_SUCCESS);
    PetscErrorCode ierr;
    ierr = DMPlexComputeResidualHybridByKey(dm, keys_, cohesive_cells_is_,
                                            t, locU, locU_t, s, locF, ctx);
    CHKERRQ(ierr);
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode FaultCohesiveKin::computeHybridJacobian(DM dm, PetscReal t,
                                                         PetscReal a,
                                                         Vec locU, Vec locU_t,
                                                         Mat J, Mat P,
                                                         void* ctx)
{
    PetscFunctionBeginUser;
    if (!kernels_set_) PetscFunctionReturn(PETSC_SUCCESS);
    PetscErrorCode ierr;
    ierr = DMPlexComputeJacobianHybridByKey(dm, keys_, cohesive_cells_is_,
                                             t, a, locU, locU_t, J, P, ctx);
    CHKERRQ(ierr);
    PetscFunctionReturn(PETSC_SUCCESS);
}

}  // namespace faults
}  // namespace FSRM
