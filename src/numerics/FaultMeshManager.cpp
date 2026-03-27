/**
 * @file FaultMeshManager.cpp
 * @brief Mesh surgery for cohesive fault insertion using PETSc DMPlex
 *
 * Implements the mesh splitting workflow:
 *   1. Label fault faces in the DMPlex mesh
 *   2. Call DMPlexConstructCohesiveCells to insert zero-volume cells
 *   3. Extract topology for FaultCohesiveDyn
 *
 * References:
 *   Aagaard et al. (2013): Domain decomposition for fault slip in FEM models
 *   PETSc DMPlex documentation: DMPlexConstructCohesiveCells
 */

#include "numerics/FaultMeshManager.hpp"
#include "domain/geomechanics/PyLithFault.hpp"
#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>

namespace FSRM {

FaultMeshManager::FaultMeshManager(MPI_Comm comm)
    : comm_(comm), rank_(0), num_cohesive_cells_(0) {
    MPI_Comm_rank(comm_, &rank_);
}

PetscErrorCode FaultMeshManager::identifyFaultSurface(DM dm, const std::string& label_name,
                                                       DMLabel* fault_label) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

    ierr = DMGetLabel(dm, label_name.c_str(), fault_label); CHKERRQ(ierr);
    if (!*fault_label) {
        if (rank_ == 0) {
            PetscPrintf(comm_, "FaultMeshManager: Label '%s' not found in DM\n",
                       label_name.c_str());
        }
        SETERRQ(comm_, PETSC_ERR_ARG_WRONG,
                "Fault label not found in DMPlex mesh");
    }

    if (rank_ == 0) {
        PetscPrintf(comm_, "FaultMeshManager: Found fault label '%s'\n",
                   label_name.c_str());
    }

    PetscFunctionReturn(0);
}

PetscErrorCode FaultMeshManager::createPlanarFaultLabel(DM dm, DMLabel* fault_label,
                                                         double strike, double dip,
                                                         const double center[3],
                                                         double length, double width,
                                                         double tolerance) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    PetscInt dim;

    ierr = DMGetDimension(dm, &dim); CHKERRQ(ierr);

    // Compute fault plane normal and orientation vectors
    // Strike: angle from north (y-axis) clockwise
    // Dip: angle from horizontal
    double cos_s = std::cos(strike);
    double sin_s = std::sin(strike);
    double cos_d = std::cos(dip);
    double sin_d = std::sin(dip);

    // Strike direction (along-strike, horizontal)
    double strike_dir[3] = {sin_s, cos_s, 0.0};

    // Dip direction (down-dip, into the fault)
    double dip_dir[3] = {-cos_s * cos_d, sin_s * cos_d, -sin_d};

    // Fault normal (pointing from - side to + side)
    double normal[3] = {cos_s * sin_d, -sin_s * sin_d, cos_d};

    // Create the label
    ierr = DMCreateLabel(dm, "fault"); CHKERRQ(ierr);
    ierr = DMGetLabel(dm, "fault", fault_label); CHKERRQ(ierr);

    // Iterate over all faces (codimension 1 entities = height 1)
    PetscInt fStart, fEnd;
    ierr = DMPlexGetHeightStratum(dm, 1, &fStart, &fEnd); CHKERRQ(ierr);

    PetscInt num_marked = 0;
    for (PetscInt f = fStart; f < fEnd; ++f) {
        // Compute face centroid
        PetscReal centroid[3] = {0.0, 0.0, 0.0};
        PetscReal vol;
        ierr = DMPlexComputeCellGeometryFVM(dm, f, &vol, centroid, nullptr); CHKERRQ(ierr);

        double point[3] = {centroid[0], centroid[1], centroid[2]};

        // Check if face centroid lies on the fault plane
        if (isPointOnFaultPlane(point, normal, center, tolerance) &&
            isPointWithinFaultExtent(point, center, strike_dir, dip_dir, length, width)) {
            ierr = DMLabelSetValue(*fault_label, f, 1); CHKERRQ(ierr);
            num_marked++;
        }
    }

    if (rank_ == 0) {
        PetscPrintf(comm_, "FaultMeshManager: Marked %d faces on fault plane "
                   "(strike=%.1f°, dip=%.1f°, center=[%.0f,%.0f,%.0f], "
                   "L=%.0f, W=%.0f, tol=%.1f)\n",
                   (int)num_marked, strike * 180.0 / M_PI, dip * 180.0 / M_PI,
                   center[0], center[1], center[2], length, width, tolerance);
    }

    PetscFunctionReturn(0);
}

PetscErrorCode FaultMeshManager::splitMeshAlongFault(DM* dm,
                                                      const std::string& fault_label_name) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

    DMLabel fault_label = nullptr;
    ierr = DMGetLabel(*dm, fault_label_name.c_str(), &fault_label); CHKERRQ(ierr);
    if (!fault_label) {
        SETERRQ(comm_, PETSC_ERR_ARG_WRONG,
                "Fault label not found; call createPlanarFaultLabel or identifyFaultSurface first");
    }

    // Count cells before splitting
    PetscInt cStart_before, cEnd_before;
    ierr = DMPlexGetHeightStratum(*dm, 0, &cStart_before, &cEnd_before); CHKERRQ(ierr);
    PetscInt num_cells_before = cEnd_before - cStart_before;

    if (rank_ == 0) {
        PetscPrintf(comm_, "FaultMeshManager: Splitting mesh along '%s' "
                   "(%d cells before split)\n",
                   fault_label_name.c_str(), (int)num_cells_before);
    }

    // Insert cohesive cells along the fault
    DM sdm = nullptr;
    ierr = DMPlexConstructCohesiveCells(*dm, fault_label, nullptr, &sdm); CHKERRQ(ierr);

    // Replace the old DM with the split DM
    ierr = DMDestroy(dm); CHKERRQ(ierr);
    *dm = sdm;

    // Count cells after splitting
    PetscInt cStart_after, cEnd_after;
    ierr = DMPlexGetHeightStratum(*dm, 0, &cStart_after, &cEnd_after); CHKERRQ(ierr);
    PetscInt num_cells_after = cEnd_after - cStart_after;
    num_cohesive_cells_ = num_cells_after - num_cells_before;

    if (rank_ == 0) {
        PetscPrintf(comm_, "FaultMeshManager: Split complete: %d cells after split "
                   "(%d cohesive cells inserted)\n",
                   (int)num_cells_after, (int)num_cohesive_cells_);
    }

    PetscFunctionReturn(0);
}

PetscErrorCode FaultMeshManager::extractCohesiveTopology(DM dm, FaultCohesiveDyn* fault) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

    if (!fault) {
        SETERRQ(comm_, PETSC_ERR_ARG_NULL, "FaultCohesiveDyn pointer is null");
    }

    PetscInt dim;
    ierr = DMGetDimension(dm, &dim); CHKERRQ(ierr);

    // Identify cohesive cells by their cell type
    PetscInt cStart, cEnd;
    ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd); CHKERRQ(ierr);

    std::vector<FaultVertex> vertices;
    std::vector<CohesiveCell> cells;

    PetscInt cohesive_count = 0;
    for (PetscInt c = cStart; c < cEnd; ++c) {
        DMPolytopeType ct;
        ierr = DMPlexGetCellType(dm, c, &ct); CHKERRQ(ierr);

        // Check for tensor product prism cell types (cohesive cells)
        bool is_cohesive = false;
        if (ct == DM_POLYTOPE_SEG_PRISM_TENSOR ||      // 2D cohesive (segment prism)
            ct == DM_POLYTOPE_TRI_PRISM_TENSOR ||       // 3D cohesive (triangle prism)
            ct == DM_POLYTOPE_QUAD_PRISM_TENSOR) {      // 3D cohesive (quad prism)
            is_cohesive = true;
        }

        if (!is_cohesive) continue;
        cohesive_count++;

        // Get the cone (faces) of the cohesive cell
        const PetscInt* cone = nullptr;
        PetscInt cone_size;
        ierr = DMPlexGetConeSize(dm, c, &cone_size); CHKERRQ(ierr);
        ierr = DMPlexGetCone(dm, c, &cone); CHKERRQ(ierr);

        // For a cohesive cell (tensor product), the cone has:
        //   - First face: negative side face
        //   - Second face: positive side face
        //   - Remaining: lateral faces (connecting - to +)
        if (cone_size < 2) continue;

        PetscInt neg_face = cone[0];
        PetscInt pos_face = cone[1];

        // Get vertices of each face
        const PetscInt* neg_verts = nullptr;
        const PetscInt* pos_verts = nullptr;
        PetscInt neg_nverts, pos_nverts;
        ierr = DMPlexGetConeSize(dm, neg_face, &neg_nverts); CHKERRQ(ierr);
        ierr = DMPlexGetCone(dm, neg_face, &neg_verts); CHKERRQ(ierr);
        ierr = DMPlexGetConeSize(dm, pos_face, &pos_nverts); CHKERRQ(ierr);
        ierr = DMPlexGetCone(dm, pos_face, &pos_verts); CHKERRQ(ierr);

        CohesiveCell cell;
        cell.cell_id = static_cast<int>(c);
        cell.vertices_negative.resize(neg_nverts);
        cell.vertices_positive.resize(pos_nverts);
        for (PetscInt i = 0; i < neg_nverts; ++i) {
            cell.vertices_negative[i] = static_cast<int>(neg_verts[i]);
        }
        for (PetscInt i = 0; i < pos_nverts; ++i) {
            cell.vertices_positive[i] = static_cast<int>(pos_verts[i]);
        }
        cells.push_back(cell);

        // Create FaultVertex entries for the vertex pairs
        // Pair neg_verts[i] with pos_verts[i] (PETSc maintains correspondence)
        PetscInt n_pairs = std::min(neg_nverts, pos_nverts);
        for (PetscInt i = 0; i < n_pairs; ++i) {
            FaultVertex fv;
            fv.vertex_id = static_cast<int>(vertices.size());
            fv.vertex_negative = static_cast<int>(neg_verts[i]);
            fv.vertex_positive = static_cast<int>(pos_verts[i]);
            // Leave Lagrange multiplier DOF unset here; -1 indicates "no DOF assigned".
            // Any Lagrange multiplier field must assign vertex_lagrange during field/constraint setup.
            fv.vertex_lagrange = -1;

            // Get coordinates of the negative vertex
            PetscSection coord_sec;
            Vec coord_vec;
            ierr = DMGetCoordinateSection(dm, &coord_sec); CHKERRQ(ierr);
            ierr = DMGetCoordinatesLocal(dm, &coord_vec); CHKERRQ(ierr);

            const PetscScalar* coords;
            PetscInt cdof, coff;
            ierr = VecGetArrayRead(coord_vec, &coords); CHKERRQ(ierr);
            ierr = PetscSectionGetDof(coord_sec, neg_verts[i], &cdof); CHKERRQ(ierr);
            ierr = PetscSectionGetOffset(coord_sec, neg_verts[i], &coff); CHKERRQ(ierr);

            fv.coords[0] = (cdof > 0) ? PetscRealPart(coords[coff + 0]) : 0.0;
            fv.coords[1] = (cdof > 1) ? PetscRealPart(coords[coff + 1]) : 0.0;
            fv.coords[2] = (cdof > 2) ? PetscRealPart(coords[coff + 2]) : 0.0;

            ierr = VecRestoreArrayRead(coord_vec, &coords); CHKERRQ(ierr);

            // Default orientation (will be computed in computeFaultGeometry)
            fv.normal = {0.0, 0.0, 1.0};
            fv.along_strike = {1.0, 0.0, 0.0};
            fv.up_dip = {0.0, 1.0, 0.0};
            fv.area = 0.0;

            vertices.push_back(fv);
        }
    }

    // Set the topology on the fault object
    fault->setFaultVertices(vertices);
    fault->setCohesiveCells(cells);

    if (rank_ == 0) {
        PetscPrintf(comm_, "FaultMeshManager: Extracted %d cohesive cells, "
                   "%d fault vertices\n",
                   (int)cells.size(), (int)vertices.size());
    }

    PetscFunctionReturn(0);
}

PetscErrorCode FaultMeshManager::computeFaultGeometry(DM dm, FaultCohesiveDyn* fault) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

    if (!fault) {
        SETERRQ(comm_, PETSC_ERR_ARG_NULL, "FaultCohesiveDyn pointer is null");
    }

    PetscInt dim;
    ierr = DMGetDimension(dm, &dim); CHKERRQ(ierr);

    // Get coordinate access
    PetscSection coord_sec;
    Vec coord_vec;
    ierr = DMGetCoordinateSection(dm, &coord_sec); CHKERRQ(ierr);
    ierr = DMGetCoordinatesLocal(dm, &coord_vec); CHKERRQ(ierr);

    const PetscScalar* coords;
    ierr = VecGetArrayRead(coord_vec, &coords); CHKERRQ(ierr);

    auto getCoords = [&](PetscInt vertex, double out[3]) -> PetscErrorCode {
        PetscInt cdof, coff;
        PetscErrorCode err;
        err = PetscSectionGetDof(coord_sec, vertex, &cdof); CHKERRQ(err);
        err = PetscSectionGetOffset(coord_sec, vertex, &coff); CHKERRQ(err);
        out[0] = (cdof > 0) ? PetscRealPart(coords[coff + 0]) : 0.0;
        out[1] = (cdof > 1) ? PetscRealPart(coords[coff + 1]) : 0.0;
        out[2] = (cdof > 2) ? PetscRealPart(coords[coff + 2]) : 0.0;
        return 0;
    };

    // For each cohesive cell, compute the normal from the face geometry
    // and assign to the corresponding fault vertices
    size_t nCells = fault->numCells();
    for (size_t ci = 0; ci < nCells; ++ci) {
        const CohesiveCell& cell = fault->getCell(ci);

        if (cell.vertices_negative.size() < 2) continue;

        // Compute face centroid and normal from negative face vertices
        double centroid[3] = {0.0, 0.0, 0.0};
        size_t nv = cell.vertices_negative.size();
        std::vector<double> vcoords(nv * 3);
        for (size_t i = 0; i < nv; ++i) {
            double vc[3];
            ierr = getCoords(cell.vertices_negative[i], vc); CHKERRQ(ierr);
            vcoords[3 * i + 0] = vc[0];
            vcoords[3 * i + 1] = vc[1];
            vcoords[3 * i + 2] = vc[2];
            centroid[0] += vc[0];
            centroid[1] += vc[1];
            centroid[2] += vc[2];
        }
        centroid[0] /= nv;
        centroid[1] /= nv;
        centroid[2] /= nv;

        // Compute normal using first two edges (cross product)
        double normal[3] = {0.0, 0.0, 1.0};
        double area = 0.0;

        if (nv >= 3) {
            double e1[3] = {vcoords[3] - vcoords[0],
                           vcoords[4] - vcoords[1],
                           vcoords[5] - vcoords[2]};
            double e2[3] = {vcoords[6] - vcoords[0],
                           vcoords[7] - vcoords[1],
                           vcoords[8] - vcoords[2]};

            // Cross product
            normal[0] = e1[1] * e2[2] - e1[2] * e2[1];
            normal[1] = e1[2] * e2[0] - e1[0] * e2[2];
            normal[2] = e1[0] * e2[1] - e1[1] * e2[0];

            double norm = std::sqrt(normal[0] * normal[0] +
                                    normal[1] * normal[1] +
                                    normal[2] * normal[2]);
            if (norm > 1e-15) {
                normal[0] /= norm;
                normal[1] /= norm;
                normal[2] /= norm;
                // Area is half the cross product magnitude for a triangle,
                // full for a parallelogram
                area = (nv == 3) ? 0.5 * norm : norm;
            }
        } else if (nv == 2 && dim == 2) {
            // 2D: "face" is an edge, normal is perpendicular to edge
            double e[3] = {vcoords[3] - vcoords[0],
                          vcoords[4] - vcoords[1],
                          0.0};
            double edge_len = std::sqrt(e[0] * e[0] + e[1] * e[1]);
            if (edge_len > 1e-15) {
                normal[0] = -e[1] / edge_len;
                normal[1] =  e[0] / edge_len;
                normal[2] = 0.0;
                area = edge_len;
            }
        }

        // Compute strike and dip directions from normal
        // Strike direction: horizontal component perpendicular to normal
        double strike_dir[3], dip_dir[3];
        if (std::abs(normal[2]) > 0.999) {
            // Nearly horizontal fault: use x as strike
            strike_dir[0] = 1.0; strike_dir[1] = 0.0; strike_dir[2] = 0.0;
        } else {
            // Strike = cross(normal, z-hat) normalized
            strike_dir[0] = -normal[1];
            strike_dir[1] =  normal[0];
            strike_dir[2] = 0.0;
            double sn = std::sqrt(strike_dir[0] * strike_dir[0] +
                                   strike_dir[1] * strike_dir[1]);
            if (sn > 1e-15) {
                strike_dir[0] /= sn;
                strike_dir[1] /= sn;
            }
        }
        // Dip direction = cross(strike, normal)
        dip_dir[0] = strike_dir[1] * normal[2] - strike_dir[2] * normal[1];
        dip_dir[1] = strike_dir[2] * normal[0] - strike_dir[0] * normal[2];
        dip_dir[2] = strike_dir[0] * normal[1] - strike_dir[1] * normal[0];

        // This is a const-correct workaround: we need mutable access to vertices.
        // FaultCohesive stores vertices internally and provides const access.
        // The proper way is to build a new vertex list and call setFaultVertices.
        // For now, we store results and update at the end.
        (void)centroid;
        (void)area;
        (void)strike_dir;
        (void)dip_dir;
        (void)normal;
    }

    // Re-build vertex list with updated geometry
    // For simplicity, iterate over all vertices and update using cell adjacency.
    // Since multiple cells may share vertices, accumulate and average.
    std::vector<FaultVertex> updated_vertices;
    size_t nv_total = fault->numVertices();
    updated_vertices.reserve(nv_total);
    for (size_t i = 0; i < nv_total; ++i) {
        FaultVertex v = fault->getVertex(i);
        // Keep existing coordinates, set default geometry
        // In a full implementation, compute from adjacent cell normals
        updated_vertices.push_back(v);
    }

    // Update each cell's vertices with the computed geometry
    for (size_t ci = 0; ci < nCells; ++ci) {
        const CohesiveCell& cell = fault->getCell(ci);
        if (cell.vertices_negative.size() < 2) continue;

        // Recompute normal for this cell
        size_t nv_cell = cell.vertices_negative.size();
        std::vector<double> vc(nv_cell * 3);
        for (size_t i = 0; i < nv_cell; ++i) {
            double vcc[3];
            ierr = getCoords(cell.vertices_negative[i], vcc); CHKERRQ(ierr);
            vc[3 * i + 0] = vcc[0];
            vc[3 * i + 1] = vcc[1];
            vc[3 * i + 2] = vcc[2];
        }

        double face_normal[3] = {0, 0, 1};
        double face_area = 0.0;

        if (nv_cell >= 3) {
            double e1[3] = {vc[3] - vc[0], vc[4] - vc[1], vc[5] - vc[2]};
            double e2[3] = {vc[6] - vc[0], vc[7] - vc[1], vc[8] - vc[2]};
            face_normal[0] = e1[1] * e2[2] - e1[2] * e2[1];
            face_normal[1] = e1[2] * e2[0] - e1[0] * e2[2];
            face_normal[2] = e1[0] * e2[1] - e1[1] * e2[0];
            double norm = std::sqrt(face_normal[0] * face_normal[0] +
                                    face_normal[1] * face_normal[1] +
                                    face_normal[2] * face_normal[2]);
            if (norm > 1e-15) {
                face_normal[0] /= norm; face_normal[1] /= norm; face_normal[2] /= norm;
                face_area = (nv_cell == 3) ? 0.5 * norm : norm;
            }
        } else if (nv_cell == 2) {
            double e[3] = {vc[3] - vc[0], vc[4] - vc[1], 0.0};
            double elen = std::sqrt(e[0] * e[0] + e[1] * e[1]);
            if (elen > 1e-15) {
                face_normal[0] = -e[1] / elen;
                face_normal[1] = e[0] / elen;
                face_normal[2] = 0.0;
                face_area = elen;
            }
        }

        // Strike direction
        double sd[3], dd[3];
        if (std::abs(face_normal[2]) > 0.999) {
            sd[0] = 1; sd[1] = 0; sd[2] = 0;
        } else {
            sd[0] = -face_normal[1]; sd[1] = face_normal[0]; sd[2] = 0;
            double sn = std::sqrt(sd[0] * sd[0] + sd[1] * sd[1]);
            if (sn > 1e-15) { sd[0] /= sn; sd[1] /= sn; }
        }
        dd[0] = sd[1] * face_normal[2] - sd[2] * face_normal[1];
        dd[1] = sd[2] * face_normal[0] - sd[0] * face_normal[2];
        dd[2] = sd[0] * face_normal[1] - sd[1] * face_normal[0];

        // Distribute area equally among cell vertices
        double area_per_vert = face_area / std::max((size_t)1, nv_cell);

        // Find which updated_vertices correspond to this cell's negative vertices
        for (size_t i = 0; i < std::min(nv_cell, (size_t)updated_vertices.size()); ++i) {
            int neg_vert = cell.vertices_negative[i];
            for (auto& uv : updated_vertices) {
                if (uv.vertex_negative == neg_vert) {
                    uv.normal = {face_normal[0], face_normal[1], face_normal[2]};
                    uv.along_strike = {sd[0], sd[1], sd[2]};
                    uv.up_dip = {dd[0], dd[1], dd[2]};
                    uv.area += area_per_vert;
                    break;
                }
            }
        }
    }

    // Set updated vertices
    fault->setFaultVertices(updated_vertices);

    ierr = VecRestoreArrayRead(coord_vec, &coords); CHKERRQ(ierr);

    if (rank_ == 0) {
        PetscPrintf(comm_, "FaultMeshManager: Computed geometry for %d fault vertices\n",
                   (int)updated_vertices.size());
    }

    PetscFunctionReturn(0);
}

bool FaultMeshManager::isPointOnFaultPlane(const double point[3], const double normal[3],
                                            const double center[3], double tolerance) const {
    // Signed distance from point to plane
    double dist = (point[0] - center[0]) * normal[0] +
                  (point[1] - center[1]) * normal[1] +
                  (point[2] - center[2]) * normal[2];
    return std::abs(dist) <= tolerance;
}

bool FaultMeshManager::isPointWithinFaultExtent(const double point[3], const double center[3],
                                                 const double strike_dir[3],
                                                 const double dip_dir[3],
                                                 double length, double width) const {
    // Project the vector (point - center) onto strike and dip directions
    double dp[3] = {point[0] - center[0],
                    point[1] - center[1],
                    point[2] - center[2]};

    double along_strike = dp[0] * strike_dir[0] + dp[1] * strike_dir[1] + dp[2] * strike_dir[2];
    double along_dip    = dp[0] * dip_dir[0]    + dp[1] * dip_dir[1]    + dp[2] * dip_dir[2];

    return (std::abs(along_strike) <= 0.5 * length) &&
           (std::abs(along_dip)    <= 0.5 * width);
}

} // namespace FSRM
