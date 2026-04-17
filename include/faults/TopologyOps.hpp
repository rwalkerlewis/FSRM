#ifndef FSRM_FAULTS_TOPOLOGY_OPS_HPP
#define FSRM_FAULTS_TOPOLOGY_OPS_HPP

#include <petsc.h>
#include <petscdm.h>
#include <petscdmplex.h>

namespace FSRM {
namespace faults {
namespace topology {

/**
 * @brief Create the "cohesive interface" DMLabel covering cohesive cell closure points.
 *
 * Mirrors PyLith `TopologyOps::getInterfacesLabel`
 * (`libsrc/pylith/topology/TopologyOps.cc:438-467` @ geodynamics/pylith
 * 6accf7e7395aa792444ad1b489d455e41ece34bc). Walks every height stratum on the
 * DMPlex, queries `DMPlexGetSimplexOrBoxCells` for the split between
 * bulk and hybrid points at that height, and tags every point at index
 * >= pMax with value 1 in a label named "cohesive interface".
 *
 * Idempotent: if the label already exists on the DM, returns it without
 * modification.
 *
 * @param[in]  dm                              DMPlex that already contains
 *                                             cohesive cells (post
 *                                             DMPlexConstructCohesiveCells).
 * @param[in]  fault_label                     Reserved for Phase 2 per-fault
 *                                             restriction. Currently unused;
 *                                             pass nullptr. Accepted for API
 *                                             forward compatibility with
 *                                             PyLith's per-fault workflow.
 * @param[out] cohesive_interface_label_out    Filled with the created or
 *                                             fetched DMLabel*. Must not be
 *                                             NULL.
 *
 * @return PETSC_SUCCESS on success.
 */
PetscErrorCode createInterfacesLabel(DM dm,
                                     DMLabel fault_label,
                                     DMLabel *cohesive_interface_label_out);

}  // namespace topology
}  // namespace faults
}  // namespace FSRM

#endif  // FSRM_FAULTS_TOPOLOGY_OPS_HPP
