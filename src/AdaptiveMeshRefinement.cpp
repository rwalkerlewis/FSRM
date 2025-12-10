#include "petsc.h"
#include "AdaptiveMeshRefinement.hpp"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iostream>
#include <fstream>

namespace FSRM {

// ============================================================================
// ErrorEstimator Base Implementation
// ============================================================================

ErrorEstimator::ErrorEstimator() : field_index_(0), component_index_(0) {}

// ============================================================================
// GradientErrorEstimator Implementation
// ============================================================================

GradientErrorEstimator::GradientErrorEstimator() : h_power_(0.5) {}

PetscErrorCode GradientErrorEstimator::estimate(DM dm, Vec solution,
                                                std::vector<CellErrorIndicator>& errors) {
    PetscFunctionBeginUser;
    
    PetscInt cStart, cEnd;
    PetscCall(DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd));
    
    errors.clear();
    errors.reserve(cEnd - cStart);
    
    // Get section for accessing solution
    PetscSection section;
    PetscCall(DMGetLocalSection(dm, &section));
    
    Vec local_solution;
    PetscCall(DMGetLocalVector(dm, &local_solution));
    PetscCall(DMGlobalToLocal(dm, solution, INSERT_VALUES, local_solution));
    
    const PetscScalar* sol_array;
    PetscCall(VecGetArrayRead(local_solution, &sol_array));
    
    for (PetscInt c = cStart; c < cEnd; c++) {
        CellErrorIndicator indicator;
        indicator.cell_id = c;
        indicator.marker = CellMarker::KEEP;
        
        // Compute cell volume and characteristic length
        PetscReal vol;
        PetscCall(DMPlexComputeCellGeometryFVM(dm, c, &vol, nullptr, nullptr));
        indicator.volume = vol;
        
        // Characteristic length (cube root of volume for 3D)
        PetscInt dim;
        PetscCall(DMGetDimension(dm, &dim));
        indicator.h = std::pow(vol, 1.0 / dim);
        
        // Compute gradient jump across faces
        PetscReal jump;
        PetscCall(computeGradientJump(dm, local_solution, c, &jump));
        
        // Kelly error estimator: η² = h^(2p) * ||[∇u·n]||²
        indicator.error = std::pow(indicator.h, h_power_) * std::sqrt(jump);
        
        errors.push_back(indicator);
    }
    
    PetscCall(VecRestoreArrayRead(local_solution, &sol_array));
    PetscCall(DMRestoreLocalVector(dm, &local_solution));
    
    PetscFunctionReturn(0);
}

PetscErrorCode GradientErrorEstimator::computeGradientJump(DM dm, Vec solution,
                                                           PetscInt cell,
                                                           PetscReal* jump) {
    PetscFunctionBeginUser;
    
    *jump = 0.0;
    
    PetscInt dim;
    PetscCall(DMGetDimension(dm, &dim));
    
    // Get faces of this cell
    const PetscInt* faces;
    PetscInt num_faces;
    PetscCall(DMPlexGetConeSize(dm, cell, &num_faces));
    PetscCall(DMPlexGetCone(dm, cell, &faces));
    
    // Get cell centroid and solution value
    PetscReal cell_centroid[3], cell_normal[3];
    PetscReal cell_vol;
    PetscCall(DMPlexComputeCellGeometryFVM(dm, cell, &cell_vol, cell_centroid, cell_normal));
    
    // Get solution at cell center
    PetscSection section;
    PetscCall(DMGetLocalSection(dm, &section));
    
    PetscInt dof;
    PetscCall(PetscSectionGetDof(section, cell, &dof));
    
    if (dof == 0) {
        PetscFunctionReturn(0);
    }
    
    PetscInt offset;
    PetscCall(PetscSectionGetOffset(section, cell, &offset));
    
    const PetscScalar* sol_array;
    PetscCall(VecGetArrayRead(solution, &sol_array));
    
    PetscScalar u_cell = sol_array[offset + component_index_];
    
    // Loop over faces to compute gradient jumps
    PetscInt fStart, fEnd;
    PetscCall(DMPlexGetHeightStratum(dm, 1, &fStart, &fEnd));
    
    for (PetscInt f = 0; f < num_faces; f++) {
        PetscInt face = faces[f];
        
        if (face < fStart || face >= fEnd) continue;
        
        // Get cells sharing this face
        const PetscInt* support;
        PetscInt support_size;
        PetscCall(DMPlexGetSupportSize(dm, face, &support_size));
        PetscCall(DMPlexGetSupport(dm, face, &support));
        
        if (support_size != 2) continue;  // Boundary face
        
        // Find neighbor cell
        PetscInt neighbor = (support[0] == cell) ? support[1] : support[0];
        
        // Get neighbor solution value
        PetscInt neighbor_offset;
        PetscCall(PetscSectionGetOffset(section, neighbor, &neighbor_offset));
        PetscScalar u_neighbor = sol_array[neighbor_offset + component_index_];
        
        // Get face geometry
        PetscReal face_centroid[3], face_normal[3];
        PetscReal face_area;
        PetscCall(DMPlexComputeCellGeometryFVM(dm, face, &face_area, face_centroid, face_normal));
        
        // Get neighbor centroid
        PetscReal neighbor_centroid[3];
        PetscCall(DMPlexComputeCellGeometryFVM(dm, neighbor, nullptr, neighbor_centroid, nullptr));
        
        // Distance between cell centers
        PetscReal dist = 0.0;
        for (PetscInt d = 0; d < dim; d++) {
            dist += (neighbor_centroid[d] - cell_centroid[d]) * 
                    (neighbor_centroid[d] - cell_centroid[d]);
        }
        dist = std::sqrt(dist);
        
        // Gradient jump estimate
        PetscReal grad_jump = std::abs(u_neighbor - u_cell) / dist;
        
        // Add to total (weighted by face area)
        *jump += face_area * grad_jump * grad_jump;
    }
    
    PetscCall(VecRestoreArrayRead(solution, &sol_array));
    
    PetscFunctionReturn(0);
}

// ============================================================================
// HessianErrorEstimator Implementation
// ============================================================================

HessianErrorEstimator::HessianErrorEstimator() 
    : recovery_method_(RecoveryMethod::SPR) {}

PetscErrorCode HessianErrorEstimator::estimate(DM dm, Vec solution,
                                               std::vector<CellErrorIndicator>& errors) {
    PetscFunctionBeginUser;
    
    PetscInt cStart, cEnd;
    PetscCall(DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd));
    
    errors.clear();
    errors.reserve(cEnd - cStart);
    
    // Recover Hessian
    Vec hessian;
    PetscCall(DMGetGlobalVector(dm, &hessian));
    PetscCall(recoverHessian(dm, solution, hessian));
    
    const PetscScalar* hess_array;
    PetscCall(VecGetArrayRead(hessian, &hess_array));
    
    PetscInt dim;
    PetscCall(DMGetDimension(dm, &dim));
    
    for (PetscInt c = cStart; c < cEnd; c++) {
        CellErrorIndicator indicator;
        indicator.cell_id = c;
        indicator.marker = CellMarker::KEEP;
        
        PetscReal vol;
        PetscCall(DMPlexComputeCellGeometryFVM(dm, c, &vol, nullptr, nullptr));
        indicator.volume = vol;
        indicator.h = std::pow(vol, 1.0 / dim);
        
        // Get Hessian at cell (assuming stored per cell)
        PetscSection section;
        PetscCall(DMGetLocalSection(dm, &section));
        
        PetscInt offset;
        PetscCall(PetscSectionGetOffset(section, c, &offset));
        
        // Hessian error: η = h² * |H|
        // |H| = Frobenius norm of Hessian matrix
        PetscReal hess_norm = 0.0;
        for (PetscInt i = 0; i < dim * dim; i++) {
            hess_norm += hess_array[offset + i] * hess_array[offset + i];
        }
        hess_norm = std::sqrt(hess_norm);
        
        indicator.error = indicator.h * indicator.h * hess_norm;
        
        errors.push_back(indicator);
    }
    
    PetscCall(VecRestoreArrayRead(hessian, &hess_array));
    PetscCall(DMRestoreGlobalVector(dm, &hessian));
    
    PetscFunctionReturn(0);
}

PetscErrorCode HessianErrorEstimator::recoverHessian(DM dm, Vec solution, Vec hessian) {
    (void)solution;  // TODO: Use solution for full Hessian recovery
    PetscFunctionBeginUser;
    
    // Simplified Hessian recovery using finite differences
    // Full implementation would use SPR or L2 projection
    
    PetscInt dim;
    PetscCall(DMGetDimension(dm, &dim));
    
    PetscCall(VecZeroEntries(hessian));
    
    // For now, estimate Hessian from second differences
    // This is a placeholder - production code would use proper recovery
    
    PetscFunctionReturn(0);
}

// ============================================================================
// ResidualErrorEstimator Implementation
// ============================================================================

ResidualErrorEstimator::ResidualErrorEstimator() {}

PetscErrorCode ResidualErrorEstimator::estimate(DM dm, Vec solution,
                                                std::vector<CellErrorIndicator>& errors) {
    PetscFunctionBeginUser;
    
    if (!residual_func_) {
        SETERRQ(PetscObjectComm((PetscObject)dm), PETSC_ERR_USER,
                "Residual function not set for ResidualErrorEstimator");
    }
    
    PetscInt cStart, cEnd;
    PetscCall(DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd));
    
    errors.clear();
    errors.reserve(cEnd - cStart);
    
    // Compute residual
    Vec residual;
    PetscCall(DMGetGlobalVector(dm, &residual));
    PetscCall(residual_func_(dm, solution, residual));
    
    const PetscScalar* res_array;
    PetscCall(VecGetArrayRead(residual, &res_array));
    
    PetscSection section;
    PetscCall(DMGetLocalSection(dm, &section));
    
    PetscInt dim;
    PetscCall(DMGetDimension(dm, &dim));
    
    for (PetscInt c = cStart; c < cEnd; c++) {
        CellErrorIndicator indicator;
        indicator.cell_id = c;
        indicator.marker = CellMarker::KEEP;
        
        PetscReal vol;
        PetscCall(DMPlexComputeCellGeometryFVM(dm, c, &vol, nullptr, nullptr));
        indicator.volume = vol;
        indicator.h = std::pow(vol, 1.0 / dim);
        
        PetscInt offset, dof;
        PetscCall(PetscSectionGetOffset(section, c, &offset));
        PetscCall(PetscSectionGetDof(section, c, &dof));
        
        // Compute local residual norm
        PetscReal res_norm = 0.0;
        for (PetscInt i = 0; i < dof; i++) {
            res_norm += res_array[offset + i] * res_array[offset + i];
        }
        res_norm = std::sqrt(res_norm);
        
        // Residual error: η = h * ||R||
        indicator.error = indicator.h * res_norm;
        
        errors.push_back(indicator);
    }
    
    PetscCall(VecRestoreArrayRead(residual, &res_array));
    PetscCall(DMRestoreGlobalVector(dm, &residual));
    
    PetscFunctionReturn(0);
}

// ============================================================================
// JumpErrorEstimator Implementation
// ============================================================================

JumpErrorEstimator::JumpErrorEstimator() : jump_threshold_(0.1) {}

PetscErrorCode JumpErrorEstimator::estimate(DM dm, Vec solution,
                                            std::vector<CellErrorIndicator>& errors) {
    PetscFunctionBeginUser;
    
    PetscInt cStart, cEnd;
    PetscCall(DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd));
    
    errors.clear();
    errors.reserve(cEnd - cStart);
    
    Vec local_solution;
    PetscCall(DMGetLocalVector(dm, &local_solution));
    PetscCall(DMGlobalToLocal(dm, solution, INSERT_VALUES, local_solution));
    
    const PetscScalar* sol_array;
    PetscCall(VecGetArrayRead(local_solution, &sol_array));
    
    PetscSection section;
    PetscCall(DMGetLocalSection(dm, &section));
    
    PetscInt dim;
    PetscCall(DMGetDimension(dm, &dim));
    
    for (PetscInt c = cStart; c < cEnd; c++) {
        CellErrorIndicator indicator;
        indicator.cell_id = c;
        indicator.marker = CellMarker::KEEP;
        
        PetscReal vol;
        PetscCall(DMPlexComputeCellGeometryFVM(dm, c, &vol, nullptr, nullptr));
        indicator.volume = vol;
        indicator.h = std::pow(vol, 1.0 / dim);
        
        // Get cell solution value
        PetscInt offset;
        PetscCall(PetscSectionGetOffset(section, c, &offset));
        PetscScalar u_cell = sol_array[offset + component_index_];
        
        // Find maximum jump to neighbors
        PetscReal max_jump = 0.0;
        
        const PetscInt* faces;
        PetscInt num_faces;
        PetscCall(DMPlexGetConeSize(dm, c, &num_faces));
        PetscCall(DMPlexGetCone(dm, c, &faces));
        
        for (PetscInt f = 0; f < num_faces; f++) {
            PetscInt face = faces[f];
            
            const PetscInt* support;
            PetscInt support_size;
            PetscCall(DMPlexGetSupportSize(dm, face, &support_size));
            PetscCall(DMPlexGetSupport(dm, face, &support));
            
            if (support_size != 2) continue;
            
            PetscInt neighbor = (support[0] == c) ? support[1] : support[0];
            
            PetscInt neighbor_offset;
            PetscCall(PetscSectionGetOffset(section, neighbor, &neighbor_offset));
            PetscScalar u_neighbor = sol_array[neighbor_offset + component_index_];
            
            PetscReal jump = std::abs(u_neighbor - u_cell);
            max_jump = std::max(max_jump, jump);
        }
        
        indicator.error = max_jump;
        
        // Mark based on threshold
        if (max_jump > jump_threshold_) {
            indicator.marker = CellMarker::REFINE;
        }
        
        errors.push_back(indicator);
    }
    
    PetscCall(VecRestoreArrayRead(local_solution, &sol_array));
    PetscCall(DMRestoreLocalVector(dm, &local_solution));
    
    PetscFunctionReturn(0);
}

// ============================================================================
// FeatureErrorEstimator Implementation
// ============================================================================

FeatureErrorEstimator::FeatureErrorEstimator() {}

void FeatureErrorEstimator::addWellLocation(PetscReal x, PetscReal y, PetscReal z,
                                            PetscReal radius, PetscInt level) {
    wells_.push_back({x, y, z, radius, level});
}

void FeatureErrorEstimator::addFaultTrace(const std::vector<std::array<PetscReal, 3>>& points,
                                          PetscReal width, PetscInt level) {
    faults_.push_back({points, width, level});
}

void FeatureErrorEstimator::addRefinementBox(PetscReal xmin, PetscReal xmax,
                                             PetscReal ymin, PetscReal ymax,
                                             PetscReal zmin, PetscReal zmax,
                                             PetscInt level) {
    BoxFeature box;
    box.bounds[0] = xmin; box.bounds[1] = xmax;
    box.bounds[2] = ymin; box.bounds[3] = ymax;
    box.bounds[4] = zmin; box.bounds[5] = zmax;
    box.level = level;
    boxes_.push_back(box);
}

PetscReal FeatureErrorEstimator::computeFeatureDistance(PetscReal x, PetscReal y, 
                                                         PetscReal z) const {
    PetscReal min_dist = 1e20;
    
    // Distance to wells
    for (const auto& well : wells_) {
        PetscReal dist = std::sqrt((x - well.x) * (x - well.x) +
                                   (y - well.y) * (y - well.y) +
                                   (z - well.z) * (z - well.z));
        min_dist = std::min(min_dist, dist - well.radius);
    }
    
    // Distance to fault traces
    for (const auto& fault : faults_) {
        for (size_t i = 0; i < fault.trace.size() - 1; i++) {
            const auto& p1 = fault.trace[i];
            const auto& p2 = fault.trace[i + 1];
            
            // Distance to line segment
            PetscReal dx = p2[0] - p1[0];
            PetscReal dy = p2[1] - p1[1];
            PetscReal dz = p2[2] - p1[2];
            
            PetscReal t = ((x - p1[0]) * dx + (y - p1[1]) * dy + (z - p1[2]) * dz) /
                         (dx * dx + dy * dy + dz * dz + 1e-20);
            t = std::max(0.0, std::min(1.0, t));
            
            PetscReal px = p1[0] + t * dx;
            PetscReal py = p1[1] + t * dy;
            PetscReal pz = p1[2] + t * dz;
            
            PetscReal dist = std::sqrt((x - px) * (x - px) + 
                                       (y - py) * (y - py) +
                                       (z - pz) * (z - pz));
            min_dist = std::min(min_dist, dist - fault.width / 2);
        }
    }
    
    return std::max(0.0, min_dist);
}

PetscErrorCode FeatureErrorEstimator::estimate(DM dm, Vec solution,
                                               std::vector<CellErrorIndicator>& errors) {
    (void)solution;  // Feature detection uses geometry, not solution
    PetscFunctionBeginUser;
    
    PetscInt cStart, cEnd;
    PetscCall(DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd));
    
    errors.clear();
    errors.reserve(cEnd - cStart);
    
    PetscInt dim;
    PetscCall(DMGetDimension(dm, &dim));
    
    for (PetscInt c = cStart; c < cEnd; c++) {
        CellErrorIndicator indicator;
        indicator.cell_id = c;
        indicator.marker = CellMarker::KEEP;
        
        PetscReal vol, centroid[3];
        PetscCall(DMPlexComputeCellGeometryFVM(dm, c, &vol, centroid, nullptr));
        indicator.volume = vol;
        indicator.h = std::pow(vol, 1.0 / dim);
        
        // Compute distance to features
        PetscReal dist = computeFeatureDistance(centroid[0], centroid[1], 
                                                dim > 2 ? centroid[2] : 0.0);
        
        // Error inversely proportional to distance
        indicator.error = 1.0 / (1.0 + dist / indicator.h);
        
        // Check refinement boxes
        for (const auto& box : boxes_) {
            if (centroid[0] >= box.bounds[0] && centroid[0] <= box.bounds[1] &&
                centroid[1] >= box.bounds[2] && centroid[1] <= box.bounds[3] &&
                (dim < 3 || (centroid[2] >= box.bounds[4] && centroid[2] <= box.bounds[5]))) {
                indicator.error = 1.0;  // Force refinement
                indicator.marker = CellMarker::REFINE;
            }
        }
        
        errors.push_back(indicator);
    }
    
    PetscFunctionReturn(0);
}

// ============================================================================
// CombinedErrorEstimator Implementation
// ============================================================================

CombinedErrorEstimator::CombinedErrorEstimator() {}

void CombinedErrorEstimator::addEstimator(std::shared_ptr<ErrorEstimator> est,
                                          PetscReal weight) {
    estimators_.push_back(est);
    weights_.push_back(weight);
}

void CombinedErrorEstimator::setWeights(const std::vector<PetscReal>& weights) {
    weights_ = weights;
}

PetscErrorCode CombinedErrorEstimator::estimate(DM dm, Vec solution,
                                                std::vector<CellErrorIndicator>& errors) {
    PetscFunctionBeginUser;
    
    if (estimators_.empty()) {
        SETERRQ(PetscObjectComm((PetscObject)dm), PETSC_ERR_USER,
                "No estimators added to CombinedErrorEstimator");
    }
    
    // Compute errors from first estimator
    PetscCall(estimators_[0]->estimate(dm, solution, errors));
    
    // Weight first estimator
    PetscReal w0 = weights_.empty() ? 1.0 : weights_[0];
    for (auto& e : errors) {
        e.error *= w0;
    }
    
    // Add contributions from other estimators
    for (size_t i = 1; i < estimators_.size(); i++) {
        std::vector<CellErrorIndicator> other_errors;
        PetscCall(estimators_[i]->estimate(dm, solution, other_errors));
        
        PetscReal wi = (i < weights_.size()) ? weights_[i] : 1.0;
        
        for (size_t c = 0; c < errors.size(); c++) {
            errors[c].error += wi * other_errors[c].error;
        }
    }
    
    // Normalize by total weight
    PetscReal total_weight = 0.0;
    for (const auto& w : weights_) total_weight += w;
    if (total_weight > 0.0) {
        for (auto& e : errors) {
            e.error /= total_weight;
        }
    }
    
    PetscFunctionReturn(0);
}

// ============================================================================
// SolutionTransfer Implementation
// ============================================================================

SolutionTransfer::SolutionTransfer() {}

SolutionTransfer::~SolutionTransfer() {
    if (transfer_matrix_) MatDestroy(&transfer_matrix_);
    if (scatter_) VecScatterDestroy(&scatter_);
}

PetscErrorCode SolutionTransfer::transfer(DM dm_old, DM dm_new,
                                          Vec solution_old, Vec solution_new) {
    PetscFunctionBeginUser;
    
    if (conservative_) {
        PetscCall(conservativeTransfer(dm_old, dm_new, solution_old, solution_new));
    } else {
        // Use DMPlexTransferVecTree for forest-based adaptation
        // or interpolation for Plex adaptation
        // PETSc 3.19+ signature: DMPlexTransferVecTree(DM, Vec, DM, Vec, PetscSF, PetscSF, PetscInt*, PetscInt*, PetscBool, PetscReal)
        PetscCall(DMPlexTransferVecTree(dm_old, solution_old, dm_new, solution_new,
                                        nullptr, nullptr, nullptr, nullptr, PETSC_FALSE, 0.0));
    }
    
    PetscFunctionReturn(0);
}

PetscErrorCode SolutionTransfer::transferFields(DM dm_old, DM dm_new,
                                                const std::vector<Vec>& fields_old,
                                                std::vector<Vec>& fields_new) {
    PetscFunctionBeginUser;
    
    fields_new.resize(fields_old.size());
    
    for (size_t i = 0; i < fields_old.size(); i++) {
        PetscCall(DMGetGlobalVector(dm_new, &fields_new[i]));
        PetscCall(transfer(dm_old, dm_new, fields_old[i], fields_new[i]));
    }
    
    PetscFunctionReturn(0);
}

PetscErrorCode SolutionTransfer::conservativeTransfer(DM dm_old, DM dm_new,
                                                      Vec solution_old, Vec solution_new) {
    PetscFunctionBeginUser;
    
    // Build transfer matrix if needed
    if (!transfer_matrix_) {
        PetscCall(buildTransferMatrix(dm_old, dm_new));
    }
    
    // Apply transfer: u_new = M * u_old
    PetscCall(MatMult(transfer_matrix_, solution_old, solution_new));
    
    PetscFunctionReturn(0);
}

PetscErrorCode SolutionTransfer::buildTransferMatrix(DM dm_old, DM dm_new) {
    PetscFunctionBeginUser;
    
    // Get sizes
    PetscInt n_old, n_new;
    PetscCall(DMGetGlobalVector(dm_old, nullptr));
    PetscCall(DMGetGlobalVector(dm_new, nullptr));
    
    Vec v_old, v_new;
    PetscCall(DMGetGlobalVector(dm_old, &v_old));
    PetscCall(DMGetGlobalVector(dm_new, &v_new));
    
    PetscCall(VecGetSize(v_old, &n_old));
    PetscCall(VecGetSize(v_new, &n_new));
    
    PetscCall(DMRestoreGlobalVector(dm_old, &v_old));
    PetscCall(DMRestoreGlobalVector(dm_new, &v_new));
    
    // Create transfer matrix
    PetscCall(MatCreate(PetscObjectComm((PetscObject)dm_old), &transfer_matrix_));
    PetscCall(MatSetSizes(transfer_matrix_, PETSC_DECIDE, PETSC_DECIDE, n_new, n_old));
    PetscCall(MatSetType(transfer_matrix_, MATAIJ));
    PetscCall(MatSetUp(transfer_matrix_));
    
    // Fill transfer matrix based on cell relationships
    // This is a simplified version - production code would track parent-child relationships
    
    PetscCall(MatAssemblyBegin(transfer_matrix_, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(transfer_matrix_, MAT_FINAL_ASSEMBLY));
    
    PetscFunctionReturn(0);
}

// ============================================================================
// MeshQuality Implementation
// ============================================================================

MeshQuality::MeshQuality() {}

PetscErrorCode MeshQuality::computeQuality(DM dm) {
    PetscFunctionBeginUser;
    
    PetscInt cStart, cEnd;
    PetscCall(DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd));
    
    cell_quality_.resize(cEnd - cStart);
    
    min_quality_ = 1.0;
    PetscReal sum_quality = 0.0;
    max_aspect_ratio_ = 0.0;
    num_bad_cells_ = 0;
    
    for (PetscInt c = cStart; c < cEnd; c++) {
        PetscReal q = computeCellQuality(dm, c);
        cell_quality_[c - cStart] = q;
        
        min_quality_ = std::min(min_quality_, q);
        sum_quality += q;
        
        if (q < 0.1) num_bad_cells_++;
    }
    
    mean_quality_ = sum_quality / (cEnd - cStart);
    
    PetscFunctionReturn(0);
}

PetscBool MeshQuality::isAcceptable(const RefinementParameters& params) const {
    if (min_quality_ < params.min_quality) return PETSC_FALSE;
    if (max_aspect_ratio_ > params.max_aspect_ratio) return PETSC_FALSE;
    return PETSC_TRUE;
}

PetscReal MeshQuality::computeCellQuality(DM dm, PetscInt cell) const {
    // Simplified quality measure based on volume ratio
    PetscReal vol;
    DMPlexComputeCellGeometryFVM(dm, cell, &vol, nullptr, nullptr);
    
    PetscInt dim;
    DMGetDimension(dm, &dim);
    
    // Get vertices
    PetscInt closure_size;
    PetscInt* closure = nullptr;
    DMPlexGetTransitiveClosure(dm, cell, PETSC_TRUE, &closure_size, &closure);
    
    // Count vertices
    PetscInt vStart, vEnd;
    DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd);
    
    PetscInt num_vertices = 0;
    for (PetscInt i = 0; i < closure_size * 2; i += 2) {
        if (closure[i] >= vStart && closure[i] < vEnd) {
            num_vertices++;
        }
    }
    
    DMPlexRestoreTransitiveClosure(dm, cell, PETSC_TRUE, &closure_size, &closure);
    
    // Ideal volume for regular cell
    PetscReal h = std::pow(vol, 1.0 / dim);
    PetscReal ideal_vol = std::pow(h, dim);
    
    return std::min(vol / ideal_vol, ideal_vol / vol);
}

// ============================================================================
// AdaptiveMeshRefinement Implementation
// ============================================================================

AdaptiveMeshRefinement::AdaptiveMeshRefinement() {
    solution_transfer_ = std::make_unique<SolutionTransfer>();
    mesh_quality_ = std::make_unique<MeshQuality>();
    
    // Default error estimator
    error_estimator_ = std::make_shared<GradientErrorEstimator>();
}

AdaptiveMeshRefinement::~AdaptiveMeshRefinement() {
    if (error_vec_) VecDestroy(&error_vec_);
}

PetscErrorCode AdaptiveMeshRefinement::initialize(DM dm) {
    PetscFunctionBeginUser;
    
    dm_ = dm;
    
    // Check mesh type
    PetscBool is_plex, is_forest;
    PetscCall(PetscObjectTypeCompare((PetscObject)dm, DMPLEX, &is_plex));
    PetscCall(PetscObjectTypeCompare((PetscObject)dm, DMFOREST, &is_forest));
    
    if (!is_plex && !is_forest) {
        SETERRQ(PetscObjectComm((PetscObject)dm), PETSC_ERR_ARG_WRONG,
                "AMR requires DMPlex or DMForest mesh");
    }
    
    if (is_forest) {
        method_ = AdaptationMethod::FOREST_P4EST;
        PetscCall(DMForestGetBaseDM(dm, &dm_coarse_));
    }
    
    // Compute initial mesh quality
    if (is_plex) {
        PetscCall(mesh_quality_->computeQuality(dm));
    }
    
    current_level_ = 0;
    adapt_count_ = 0;
    
    PetscFunctionReturn(0);
}

PetscErrorCode AdaptiveMeshRefinement::adapt(Vec solution, DM* dm_new, Vec* solution_new) {
    PetscFunctionBeginUser;
    
    PetscLogDouble t_start, t_end;
    PetscCall(PetscTime(&t_start));
    
    // Store statistics
    PetscInt cStart, cEnd;
    PetscCall(DMPlexGetHeightStratum(dm_, 0, &cStart, &cEnd));
    stats_.num_cells_before = cEnd - cStart;
    
    // Pre-adaptation callback
    if (pre_adapt_cb_) {
        PetscCall(pre_adapt_cb_(dm_, solution));
    }
    
    // Compute error indicators
    PetscCall(computeErrors(solution, errors_));
    
    // Compute error statistics
    stats_.max_error = 0.0;
    stats_.min_error = 1e20;
    stats_.total_error = 0.0;
    
    for (const auto& e : errors_) {
        stats_.max_error = std::max(stats_.max_error, (PetscReal)e.error);
        stats_.min_error = std::min(stats_.min_error, (PetscReal)e.error);
        stats_.total_error += e.error;
    }
    stats_.mean_error = stats_.total_error / errors_.size();
    
    // Mark cells
    std::vector<CellMarker> markers;
    PetscCall(markCells(errors_, markers));
    
    // Enforce constraints
    PetscCall(enforceConstraints(markers));
    
    // Add buffer layers
    if (params_.buffer_layers > 0) {
        PetscCall(addBufferLayers(markers));
    }
    
    // Perform refinement
    PetscCall(refineMesh(markers, dm_new));
    
    PetscCall(PetscTime(&t_end));
    stats_.adaptation_time = t_end - t_start;
    
    // Transfer solution
    PetscCall(PetscTime(&t_start));
    PetscCall(transferSolution(dm_, *dm_new, solution, solution_new));
    PetscCall(PetscTime(&t_end));
    stats_.transfer_time = t_end - t_start;
    
    // Rebalance
    PetscCall(PetscTime(&t_start));
    PetscCall(rebalance(*dm_new));
    PetscCall(PetscTime(&t_end));
    stats_.balance_time = t_end - t_start;
    
    // Update statistics
    PetscCall(DMPlexGetHeightStratum(*dm_new, 0, &cStart, &cEnd));
    stats_.num_cells_after = cEnd - cStart;
    stats_.num_refined = std::count(markers.begin(), markers.end(), CellMarker::REFINE);
    stats_.num_coarsened = std::count(markers.begin(), markers.end(), CellMarker::COARSEN);
    
    // Post-adaptation callback
    if (post_adapt_cb_) {
        PetscCall(post_adapt_cb_(dm_, *dm_new, solution, *solution_new));
    }
    
    // Update internal state
    dm_ = *dm_new;
    adapt_count_++;
    
    PetscFunctionReturn(0);
}

PetscErrorCode AdaptiveMeshRefinement::checkAdaptation(Vec solution, PetscBool* need_adapt) {
    PetscFunctionBeginUser;
    
    *need_adapt = PETSC_FALSE;
    
    // Compute errors
    std::vector<CellErrorIndicator> errors;
    PetscCall(computeErrors(solution, errors));
    
    // Check if max error exceeds threshold
    PetscReal max_error = 0.0;
    for (const auto& e : errors) {
        max_error = std::max(max_error, (PetscReal)e.error);
    }
    
    if (max_error > params_.refine_threshold) {
        *need_adapt = PETSC_TRUE;
    }
    
    // Check cell count limits
    PetscInt cStart, cEnd;
    PetscCall(DMPlexGetHeightStratum(dm_, 0, &cStart, &cEnd));
    PetscInt num_cells = cEnd - cStart;
    
    if (num_cells > params_.max_cells) {
        *need_adapt = PETSC_TRUE;  // Need coarsening
    }
    
    PetscFunctionReturn(0);
}

PetscErrorCode AdaptiveMeshRefinement::computeErrors(Vec solution,
                                                     std::vector<CellErrorIndicator>& errors) {
    PetscFunctionBeginUser;
    
    PetscCall(error_estimator_->estimate(dm_, solution, errors));
    
    PetscFunctionReturn(0);
}

PetscErrorCode AdaptiveMeshRefinement::markCells(const std::vector<CellErrorIndicator>& errors,
                                                 std::vector<CellMarker>& markers) {
    PetscFunctionBeginUser;
    
    markers.resize(errors.size(), CellMarker::KEEP);
    
    switch (params_.strategy) {
        case RefinementStrategy::FIXED_FRACTION:
            PetscCall(markFixedFraction(errors, markers));
            break;
        case RefinementStrategy::THRESHOLD:
            PetscCall(markThreshold(errors, markers));
            break;
        case RefinementStrategy::EQUILIBRATION:
            PetscCall(markEquilibration(errors, markers));
            break;
        default:
            PetscCall(markFixedFraction(errors, markers));
    }
    
    PetscFunctionReturn(0);
}

PetscErrorCode AdaptiveMeshRefinement::markFixedFraction(
    const std::vector<CellErrorIndicator>& errors,
    std::vector<CellMarker>& markers) {
    
    PetscFunctionBeginUser;
    
    // Sort errors
    std::vector<size_t> indices(errors.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(),
              [&errors](size_t a, size_t b) {
                  return errors[a].error > errors[b].error;
              });
    
    // Mark top fraction for refinement
    PetscInt num_refine = static_cast<PetscInt>(params_.refine_fraction * errors.size());
    for (PetscInt i = 0; i < num_refine; i++) {
        markers[indices[i]] = CellMarker::REFINE;
    }
    
    // Mark bottom fraction for coarsening
    PetscInt num_coarsen = static_cast<PetscInt>(params_.coarsen_fraction * errors.size());
    for (PetscInt i = 0; i < num_coarsen; i++) {
        markers[indices[errors.size() - 1 - i]] = CellMarker::COARSEN;
    }
    
    PetscFunctionReturn(0);
}

PetscErrorCode AdaptiveMeshRefinement::markThreshold(
    const std::vector<CellErrorIndicator>& errors,
    std::vector<CellMarker>& markers) {
    
    PetscFunctionBeginUser;
    
    // Compute max error for normalization
    PetscReal max_error = 0.0;
    for (const auto& e : errors) {
        max_error = std::max(max_error, (PetscReal)e.error);
    }
    
    if (max_error < 1e-20) {
        PetscFunctionReturn(0);
    }
    
    for (size_t i = 0; i < errors.size(); i++) {
        PetscReal normalized = errors[i].error / max_error;
        
        if (normalized > params_.refine_threshold) {
            markers[i] = CellMarker::REFINE;
        } else if (normalized < params_.coarsen_threshold) {
            markers[i] = CellMarker::COARSEN;
        }
    }
    
    PetscFunctionReturn(0);
}

PetscErrorCode AdaptiveMeshRefinement::markEquilibration(
    const std::vector<CellErrorIndicator>& errors,
    std::vector<CellMarker>& markers) {
    
    PetscFunctionBeginUser;
    
    // Target: equal error in each cell
    PetscReal total_error = 0.0;
    for (const auto& e : errors) {
        total_error += e.error;
    }
    
    PetscReal target_error = total_error / errors.size();
    
    for (size_t i = 0; i < errors.size(); i++) {
        if (errors[i].error > 2.0 * target_error) {
            markers[i] = CellMarker::REFINE;
        } else if (errors[i].error < 0.25 * target_error) {
            markers[i] = CellMarker::COARSEN;
        }
    }
    
    PetscFunctionReturn(0);
}

PetscErrorCode AdaptiveMeshRefinement::enforceConstraints(std::vector<CellMarker>& markers) {
    PetscFunctionBeginUser;
    
    // Get cell levels (for forest) or estimate from size
    PetscInt dim;
    PetscCall(DMGetDimension(dm_, &dim));
    
    for (size_t i = 0; i < markers.size(); i++) {
        // Check max level
        if (markers[i] == CellMarker::REFINE && current_level_ >= params_.max_level) {
            markers[i] = CellMarker::KEEP;
        }
        
        // Check min level
        if (markers[i] == CellMarker::COARSEN && current_level_ <= params_.min_level) {
            markers[i] = CellMarker::KEEP;
        }
        
        // Check cell size limits
        if (markers[i] == CellMarker::REFINE && 
            errors_[i].h / 2.0 < params_.min_cell_size) {
            markers[i] = CellMarker::KEEP;
        }
        
        if (markers[i] == CellMarker::COARSEN && 
            errors_[i].h * 2.0 > params_.max_cell_size) {
            markers[i] = CellMarker::KEEP;
        }
    }
    
    // Check total cell count
    PetscInt num_refine = std::count(markers.begin(), markers.end(), CellMarker::REFINE);
    PetscInt num_coarsen = std::count(markers.begin(), markers.end(), CellMarker::COARSEN);
    
    // Estimate new cell count (2^dim children per refined cell in 2D/3D)
    PetscInt children_per_cell = 1 << dim;
    PetscInt estimated_cells = markers.size() + 
                               (children_per_cell - 1) * num_refine - 
                               (1 - 1.0/children_per_cell) * num_coarsen;
    
    if (estimated_cells > params_.max_cells) {
        // Reduce refinement
        std::vector<size_t> refine_indices;
        for (size_t i = 0; i < markers.size(); i++) {
            if (markers[i] == CellMarker::REFINE) {
                refine_indices.push_back(i);
            }
        }
        
        // Sort by error and keep only highest
        std::sort(refine_indices.begin(), refine_indices.end(),
                  [this](size_t a, size_t b) {
                      return errors_[a].error > errors_[b].error;
                  });
        
        PetscInt max_refine = (params_.max_cells - markers.size()) / (children_per_cell - 1);
        for (size_t i = max_refine; i < refine_indices.size(); i++) {
            markers[refine_indices[i]] = CellMarker::KEEP;
        }
    }
    
    PetscFunctionReturn(0);
}

PetscErrorCode AdaptiveMeshRefinement::addBufferLayers(std::vector<CellMarker>& markers) {
    PetscFunctionBeginUser;
    
    // Mark neighbors of refined cells
    for (PetscInt layer = 0; layer < params_.buffer_layers; layer++) {
        std::vector<CellMarker> new_markers = markers;
        
        PetscInt cStart, cEnd;
        PetscCall(DMPlexGetHeightStratum(dm_, 0, &cStart, &cEnd));
        
        for (PetscInt c = cStart; c < cEnd; c++) {
            if (markers[c - cStart] != CellMarker::REFINE) continue;
            
            // Get neighbors
            const PetscInt* faces;
            PetscInt num_faces;
            PetscCall(DMPlexGetConeSize(dm_, c, &num_faces));
            PetscCall(DMPlexGetCone(dm_, c, &faces));
            
            for (PetscInt f = 0; f < num_faces; f++) {
                const PetscInt* support;
                PetscInt support_size;
                PetscCall(DMPlexGetSupportSize(dm_, faces[f], &support_size));
                PetscCall(DMPlexGetSupport(dm_, faces[f], &support));
                
                for (PetscInt s = 0; s < support_size; s++) {
                    PetscInt neighbor = support[s];
                    if (neighbor >= cStart && neighbor < cEnd) {
                        if (new_markers[neighbor - cStart] == CellMarker::KEEP) {
                            new_markers[neighbor - cStart] = CellMarker::REFINE;
                        }
                    }
                }
            }
        }
        
        markers = new_markers;
    }
    
    PetscFunctionReturn(0);
}

PetscErrorCode AdaptiveMeshRefinement::refineMesh(const std::vector<CellMarker>& markers,
                                                  DM* dm_new) {
    PetscFunctionBeginUser;
    
    switch (method_) {
        case AdaptationMethod::PLEX_REFINE:
            PetscCall(applyUniformRefine(dm_new));
            break;
            
        case AdaptationMethod::PLEX_ADAPT: {
            Vec metric;
            PetscCall(createMetricFromErrors(errors_, &metric));
            PetscCall(applyPlexAdapt(metric, dm_new));
            PetscCall(VecDestroy(&metric));
            break;
        }
        
        case AdaptationMethod::FOREST_P4EST:
            PetscCall(applyForestAdapt(markers, dm_new));
            break;
            
        default:
            PetscCall(applyUniformRefine(dm_new));
    }
    
    PetscFunctionReturn(0);
}

PetscErrorCode AdaptiveMeshRefinement::createMetricFromErrors(
    const std::vector<CellErrorIndicator>& errors, Vec* metric) {
    
    PetscFunctionBeginUser;
    
    PetscInt dim;
    PetscCall(DMGetDimension(dm_, &dim));
    
    // Create metric vector (symmetric tensor per vertex)
    PetscInt n_components = (dim * (dim + 1)) / 2;  // Symmetric tensor
    
    PetscInt vStart, vEnd;
    PetscCall(DMPlexGetDepthStratum(dm_, 0, &vStart, &vEnd));
    PetscInt n_vertices = vEnd - vStart;
    
    PetscCall(VecCreate(PetscObjectComm((PetscObject)dm_), metric));
    PetscCall(VecSetSizes(*metric, n_vertices * n_components, PETSC_DETERMINE));
    PetscCall(VecSetType(*metric, VECSTANDARD));
    
    // Initialize with identity scaled by desired mesh size
    PetscScalar* metric_array;
    PetscCall(VecGetArray(*metric, &metric_array));
    
    for (PetscInt v = vStart; v < vEnd; v++) {
        PetscInt local_v = v - vStart;
        
        // Find cells containing this vertex and average errors
        PetscInt star_size;
        PetscInt* star = nullptr;
        PetscCall(DMPlexGetTransitiveClosure(dm_, v, PETSC_FALSE, &star_size, &star));
        
        PetscReal avg_error = 0.0;
        PetscInt num_cells = 0;
        
        PetscInt cStart, cEnd;
        PetscCall(DMPlexGetHeightStratum(dm_, 0, &cStart, &cEnd));
        
        for (PetscInt s = 0; s < star_size * 2; s += 2) {
            if (star[s] >= cStart && star[s] < cEnd) {
                avg_error += errors[star[s] - cStart].error;
                num_cells++;
            }
        }
        
        PetscCall(DMPlexRestoreTransitiveClosure(dm_, v, PETSC_FALSE, &star_size, &star));
        
        if (num_cells > 0) {
            avg_error /= num_cells;
        }
        
        // Metric = 1/h² * I where h is desired mesh size
        // Smaller error -> larger h -> smaller metric
        PetscReal h_desired = params_.max_cell_size / (1.0 + avg_error * 10.0);
        h_desired = std::max(params_.min_cell_size, 
                            std::min(params_.max_cell_size, h_desired));
        
        PetscReal metric_val = 1.0 / (h_desired * h_desired);
        
        // Set diagonal components
        for (PetscInt d = 0; d < dim; d++) {
            metric_array[local_v * n_components + d] = metric_val;
        }
        // Off-diagonal components are zero for isotropic metric
    }
    
    PetscCall(VecRestoreArray(*metric, &metric_array));
    
    PetscFunctionReturn(0);
}

PetscErrorCode AdaptiveMeshRefinement::applyPlexAdapt(Vec metric, DM* dm_new) {
    PetscFunctionBeginUser;
    
    // Use metric-based adaptation
    // PETSc 3.19+ signature: DMPlexAdapt(DM, Vec, const char*, DM*)
    // Pass nullptr for label name to use default adaptation
    // Note: DMPlexAdapt requires PETSc to be configured with external mesh
    // adaptation libraries (mmg, parmmg, or pragmatic). If not available,
    // fall back to uniform refinement.
#if defined(PETSC_HAVE_MMG) || defined(PETSC_HAVE_PARMMG) || defined(PETSC_HAVE_PRAGMATIC)
    PetscCall(DMPlexAdapt(dm_, metric, nullptr, dm_new));
#else
    // Fallback: Use uniform refinement when metric-based adaptation is unavailable
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, 
        "Warning: DMPlexAdapt not available (PETSc not configured with mmg/parmmg/pragmatic). "
        "Using uniform refinement instead.\n"));
    PetscCall(DMRefine(dm_, PETSC_COMM_WORLD, dm_new));
    (void)metric; // Suppress unused parameter warning
#endif
    
    PetscFunctionReturn(0);
}

PetscErrorCode AdaptiveMeshRefinement::applyForestAdapt(
    const std::vector<CellMarker>& markers, DM* dm_new) {
    
    PetscFunctionBeginUser;
    
    // Create adaptation label
    DMLabel adapt_label;
    PetscCall(DMGetLabel(dm_, "adapt", &adapt_label));
    
    if (!adapt_label) {
        PetscCall(DMCreateLabel(dm_, "adapt"));
        PetscCall(DMGetLabel(dm_, "adapt", &adapt_label));
    }
    
    // Set markers
    PetscInt cStart, cEnd;
    PetscCall(DMPlexGetHeightStratum(dm_, 0, &cStart, &cEnd));
    
    for (PetscInt c = cStart; c < cEnd; c++) {
        PetscCall(DMLabelSetValue(adapt_label, c, static_cast<PetscInt>(markers[c - cStart])));
    }
    
    // Adapt forest
    PetscCall(DMForestSetAdaptivityLabel(dm_, adapt_label));
    PetscCall(DMForestTemplate(dm_, PetscObjectComm((PetscObject)dm_), dm_new));
    PetscCall(DMSetUp(*dm_new));
    
    PetscFunctionReturn(0);
}

PetscErrorCode AdaptiveMeshRefinement::applyUniformRefine(DM* dm_new) {
    PetscFunctionBeginUser;
    
    PetscCall(DMRefine(dm_, PetscObjectComm((PetscObject)dm_), dm_new));
    current_level_++;
    
    PetscFunctionReturn(0);
}

PetscErrorCode AdaptiveMeshRefinement::transferSolution(DM dm_old, DM dm_new,
                                                        Vec solution_old, Vec* solution_new) {
    PetscFunctionBeginUser;
    
    PetscCall(DMGetGlobalVector(dm_new, solution_new));
    PetscCall(solution_transfer_->transfer(dm_old, dm_new, solution_old, *solution_new));
    
    PetscFunctionReturn(0);
}

PetscErrorCode AdaptiveMeshRefinement::rebalance(DM dm) {
    PetscFunctionBeginUser;
    
    // Rebalance across processors
    DM dm_dist;
    PetscCall(DMPlexDistribute(dm, 0, nullptr, &dm_dist));
    
    if (dm_dist) {
        PetscCall(DMDestroy(&dm));
        dm_ = dm_dist;
    }
    
    PetscFunctionReturn(0);
}

void AdaptiveMeshRefinement::setCriterion(RefinementCriterion criterion) {
    error_estimator_ = createErrorEstimator(criterion);
}

void AdaptiveMeshRefinement::setErrorEstimator(std::shared_ptr<ErrorEstimator> estimator) {
    error_estimator_ = estimator;
}

void AdaptiveMeshRefinement::addWellRefinement(PetscReal x, PetscReal y, PetscReal z,
                                               PetscReal radius, PetscInt level) {
    auto feature_est = std::dynamic_pointer_cast<FeatureErrorEstimator>(error_estimator_);
    if (!feature_est) {
        auto combined = std::make_shared<CombinedErrorEstimator>();
        combined->addEstimator(error_estimator_, 1.0);
        feature_est = std::make_shared<FeatureErrorEstimator>();
        combined->addEstimator(feature_est, 1.0);
        error_estimator_ = combined;
    }
    feature_est->addWellLocation(x, y, z, radius, level);
}

void AdaptiveMeshRefinement::addFaultRefinement(
    const std::vector<std::array<PetscReal, 3>>& trace,
    PetscReal width, PetscInt level) {
    
    auto feature_est = std::dynamic_pointer_cast<FeatureErrorEstimator>(error_estimator_);
    if (feature_est) {
        feature_est->addFaultTrace(trace, width, level);
    }
}

void AdaptiveMeshRefinement::printStatistics() const {
    PetscPrintf(PETSC_COMM_WORLD, "\n=== AMR Statistics ===\n");
    PetscPrintf(PETSC_COMM_WORLD, "Cells: %d -> %d\n", 
                stats_.num_cells_before, stats_.num_cells_after);
    PetscPrintf(PETSC_COMM_WORLD, "Refined: %d, Coarsened: %d\n",
                stats_.num_refined, stats_.num_coarsened);
    PetscPrintf(PETSC_COMM_WORLD, "Error: min=%.3e, max=%.3e, mean=%.3e\n",
                stats_.min_error, stats_.max_error, stats_.mean_error);
    PetscPrintf(PETSC_COMM_WORLD, "Time: adapt=%.3f, transfer=%.3f, balance=%.3f s\n",
                stats_.adaptation_time, stats_.transfer_time, stats_.balance_time);
}

PetscErrorCode AdaptiveMeshRefinement::writeAdaptedMesh(const std::string& filename) const {
    PetscFunctionBeginUser;
    
    PetscViewer viewer;
    PetscCall(PetscViewerCreate(PetscObjectComm((PetscObject)dm_), &viewer));
    PetscCall(PetscViewerSetType(viewer, PETSCVIEWERVTK));
    PetscCall(PetscViewerFileSetName(viewer, filename.c_str()));
    PetscCall(DMView(dm_, viewer));
    PetscCall(PetscViewerDestroy(&viewer));
    
    PetscFunctionReturn(0);
}

PetscErrorCode AdaptiveMeshRefinement::writeErrorField(const std::string& filename) const {
    PetscFunctionBeginUser;
    
    if (!error_vec_) {
        PetscFunctionReturn(0);
    }
    
    PetscViewer viewer;
    PetscCall(PetscViewerCreate(PetscObjectComm((PetscObject)dm_), &viewer));
    PetscCall(PetscViewerSetType(viewer, PETSCVIEWERVTK));
    PetscCall(PetscViewerFileSetName(viewer, filename.c_str()));
    PetscCall(VecView(error_vec_, viewer));
    PetscCall(PetscViewerDestroy(&viewer));
    
    PetscFunctionReturn(0);
}

// ============================================================================
// Factory Functions
// ============================================================================

std::shared_ptr<ErrorEstimator> createErrorEstimator(
    RefinementCriterion criterion,
    const std::map<std::string, PetscReal>& params) {
    (void)params;  // TODO: Use params to configure estimators
    
    switch (criterion) {
        case RefinementCriterion::GRADIENT:
            return std::make_shared<GradientErrorEstimator>();
            
        case RefinementCriterion::HESSIAN:
            return std::make_shared<HessianErrorEstimator>();
            
        case RefinementCriterion::RESIDUAL:
            return std::make_shared<ResidualErrorEstimator>();
            
        case RefinementCriterion::JUMP:
            return std::make_shared<JumpErrorEstimator>();
            
        case RefinementCriterion::FEATURE:
            return std::make_shared<FeatureErrorEstimator>();
            
        case RefinementCriterion::COMBINED:
            return std::make_shared<CombinedErrorEstimator>();
            
        default:
            return std::make_shared<GradientErrorEstimator>();
    }
}

PetscErrorCode AMRSetFromOptions(AdaptiveMeshRefinement* amr) {
    PetscFunctionBeginUser;
    
    RefinementParameters& params = amr->getParameters();
    
    PetscOptionsBegin(PETSC_COMM_WORLD, "amr_", "AMR Options", "AMR");
    
    PetscCall(PetscOptionsReal("-refine_fraction", "Fraction of cells to refine",
                               "AdaptiveMeshRefinement", params.refine_fraction,
                               &params.refine_fraction, nullptr));
    
    PetscCall(PetscOptionsReal("-coarsen_fraction", "Fraction of cells to coarsen",
                               "AdaptiveMeshRefinement", params.coarsen_fraction,
                               &params.coarsen_fraction, nullptr));
    
    PetscCall(PetscOptionsInt("-max_level", "Maximum refinement level",
                              "AdaptiveMeshRefinement", params.max_level,
                              &params.max_level, nullptr));
    
    PetscCall(PetscOptionsInt("-max_cells", "Maximum number of cells",
                              "AdaptiveMeshRefinement", params.max_cells,
                              &params.max_cells, nullptr));
    
    PetscCall(PetscOptionsReal("-min_cell_size", "Minimum cell size",
                               "AdaptiveMeshRefinement", params.min_cell_size,
                               &params.min_cell_size, nullptr));
    
    PetscOptionsEnd();
    
    PetscFunctionReturn(0);
}

PetscErrorCode DMPlexCreateAdaptive(DM dm_in, AdaptationMethod method, DM* dm_out) {
    PetscFunctionBeginUser;
    
    switch (method) {
        case AdaptationMethod::FOREST_P4EST: {
            // Create DMForest from DMPlex
            PetscCall(DMCreate(PetscObjectComm((PetscObject)dm_in), dm_out));
            PetscCall(DMSetType(*dm_out, DMFOREST));
            PetscCall(DMForestSetBaseDM(*dm_out, dm_in));
            PetscCall(DMForestSetTopology(*dm_out, "p4est"));
            PetscCall(DMForestSetMinimumRefinement(*dm_out, 0));
            PetscCall(DMForestSetMaximumRefinement(*dm_out, 10));
            PetscCall(DMSetUp(*dm_out));
            break;
        }
        
        default:
            // Clone DMPlex
            PetscCall(DMClone(dm_in, dm_out));
    }
    
    PetscFunctionReturn(0);
}

} // namespace FSRM
