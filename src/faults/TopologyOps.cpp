#include "faults/TopologyOps.hpp"

namespace FSRM {
namespace faults {
namespace topology {

PetscErrorCode createInterfacesLabel(DM dm,
                                     DMLabel fault_label,
                                     DMLabel *cohesive_interface_label_out)
{
    PetscFunctionBeginUser;
    (void)fault_label;  // reserved for Phase 2 per-fault restriction

    PetscErrorCode ierr;
    const char *labelName = "cohesive interface";

    PetscCheck(cohesive_interface_label_out != nullptr,
               PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "cohesive_interface_label_out must not be NULL");

    PetscBool hasLabel = PETSC_FALSE;
    ierr = DMHasLabel(dm, labelName, &hasLabel); CHKERRQ(ierr);
    if (hasLabel) {
        ierr = DMGetLabel(dm, labelName, cohesive_interface_label_out); CHKERRQ(ierr);
        PetscFunctionReturn(PETSC_SUCCESS);
    }

    PetscInt dim = 0;
    ierr = DMGetDimension(dm, &dim); CHKERRQ(ierr);
    ierr = DMCreateLabel(dm, labelName); CHKERRQ(ierr);
    ierr = DMGetLabel(dm, labelName, cohesive_interface_label_out); CHKERRQ(ierr);

    // Walk each height stratum. DMPlexGetSimplexOrBoxCells returns the split
    // index pMax below which points are bulk (simplex/box) and at-or-above
    // which points are hybrid (cohesive prism). Tag every hybrid point.
    for (PetscInt iDim = 0; iDim <= dim; ++iDim) {
        PetscInt pStart = 0, pEnd = 0, pMax = 0;
        ierr = DMPlexGetHeightStratum(dm, iDim, &pStart, &pEnd); CHKERRQ(ierr);
        ierr = DMPlexGetSimplexOrBoxCells(dm, iDim, NULL, &pMax); CHKERRQ(ierr);
        for (PetscInt p = pMax; p < pEnd; ++p) {
            ierr = DMLabelSetValue(*cohesive_interface_label_out, p, 1); CHKERRQ(ierr);
        }
    }

    PetscFunctionReturn(PETSC_SUCCESS);
}

}  // namespace topology
}  // namespace faults
}  // namespace FSRM
