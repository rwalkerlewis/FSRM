#ifndef FAULT_MESH_MANAGER_HPP
#define FAULT_MESH_MANAGER_HPP

/**
 * @file FaultMeshManager.hpp
 * @brief Mesh surgery for cohesive fault insertion using PETSc DMPlex
 *
 * Wraps DMPlexConstructCohesiveCells() to insert zero-volume cohesive cells
 * along fault surfaces in a DMPlex mesh. This is the foundation for
 * split-node dynamic rupture modeling (Aagaard et al., 2013).
 *
 * Usage sequence:
 *   1. Create/import DMPlex mesh
 *   2. Identify or create fault surface label
 *   3. Split mesh along fault (inserts cohesive cells)
 *   4. Extract cohesive topology for FaultCohesiveDyn
 *   5. Compute fault geometry (normals, areas)
 *   6. THEN distribute mesh (DMPlexDistribute)
 *   7. THEN add fields (DMAddField / DMCreateDS)
 *
 * CRITICAL: splitMeshAlongFault() must be called BEFORE DMPlexDistribute()
 *           and BEFORE any field setup.
 */

#include <petscdmplex.h>
#include <string>
#include <vector>

namespace FSRM {

// Forward declaration
class FaultCohesiveDyn;

/**
 * @brief Manages mesh splitting for cohesive fault insertion
 */
class FaultMeshManager {
    friend class ::FaultMeshManagerTest;

public:
    FaultMeshManager(MPI_Comm comm);

    /**
     * @brief Identify fault surface from an existing DMLabel
     *
     * For Gmsh meshes with named physical groups (e.g., "fault_surface"),
     * looks up the label by name.
     *
     * @param dm The DMPlex mesh
     * @param label_name Name of the label marking fault faces
     * @param[out] fault_label Pointer to the retrieved DMLabel
     * @return PETSc error code
     */
    PetscErrorCode identifyFaultSurface(DM dm, const std::string& label_name,
                                        DMLabel* fault_label);

    /**
     * @brief Create a planar fault label from geometric definition
     *
     * Marks all mesh faces whose centroids lie within tolerance of the
     * fault plane defined by strike, dip, center, and extent.
     *
     * @param dm The DMPlex mesh
     * @param[out] fault_label Created label
     * @param strike Fault strike angle (radians, from north CW)
     * @param dip Fault dip angle (radians, from horizontal)
     * @param center Center point of the fault [x, y, z]
     * @param length Fault length along strike (m)
     * @param width Fault width down-dip (m)
     * @param tolerance Distance tolerance for face selection (m)
     * @return PETSc error code
     */
    PetscErrorCode createPlanarFaultLabel(DM dm, DMLabel* fault_label,
                                          double strike, double dip,
                                          const double center[3],
                                          double length, double width,
                                          double tolerance);

    /**
     * @brief Split the mesh along the labeled fault surface
     *
     * Calls DMPlexConstructCohesiveCells() to insert zero-volume cohesive
     * cells. The DM is modified IN PLACE (old DM destroyed, replaced by new).
     *
     * MUST be called BEFORE DMPlexDistribute() and BEFORE setupFields().
     *
     * @param[in,out] dm Pointer to DM; on output, points to the split DM
     * @param fault_label_name Name of the label marking fault faces
     * @return PETSc error code
     */
    PetscErrorCode splitMeshAlongFault(DM* dm, const std::string& fault_label_name);

    /**
     * @brief Extract cohesive cell topology after mesh splitting
     *
     * Maps DMPlex cohesive cells to FaultVertex structs with
     * vertex_negative, vertex_positive, vertex_lagrange assignments.
     *
     * @param dm The split DMPlex mesh
     * @param fault The FaultCohesiveDyn to populate with topology
     * @return PETSc error code
     */
    PetscErrorCode extractCohesiveTopology(DM dm, FaultCohesiveDyn* fault);

    /**
     * @brief Compute fault vertex areas and normals from split mesh geometry
     *
     * @param dm The split DMPlex mesh
     * @param fault The FaultCohesiveDyn to update with geometry
     * @return PETSc error code
     */
    PetscErrorCode computeFaultGeometry(DM dm, FaultCohesiveDyn* fault);

    /**
     * @brief Get the number of cohesive cells found after splitting
     */
    PetscInt getNumCohesiveCells() const { return num_cohesive_cells_; }

private:
    MPI_Comm comm_;
    int rank_;
    PetscInt num_cohesive_cells_ = 0;

    /**
     * @brief Check if a point lies on the fault plane within tolerance
     */
    bool isPointOnFaultPlane(const double point[3], const double normal[3],
                            const double center[3], double tolerance) const;

    /**
     * @brief Check if a point is within the rectangular fault extent
     */
    bool isPointWithinFaultExtent(const double point[3], const double center[3],
                                  const double strike_dir[3], const double dip_dir[3],
                                  double length, double width) const;
};

} // namespace FSRM

#endif // FAULT_MESH_MANAGER_HPP
