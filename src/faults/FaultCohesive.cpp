#include "faults/FaultCohesive.hpp"
#include "faults/TopologyOps.hpp"

namespace FSRM {
namespace faults {

PetscErrorCode FaultCohesive::createCohesiveInterfaceLabel(DM dm, DMLabel fault_label)
{
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    ierr = topology::createInterfacesLabel(dm, fault_label, &cohesive_interface_label_);
    CHKERRQ(ierr);
    PetscFunctionReturn(PETSC_SUCCESS);
}

}  // namespace faults
}  // namespace FSRM
