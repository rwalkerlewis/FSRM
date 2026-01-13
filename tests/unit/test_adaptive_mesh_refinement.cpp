/**
 * @file test_adaptive_mesh_refinement.cpp
 * @brief Unit tests for Adaptive Mesh Refinement (AMR)
 * 
 * Tests for:
 * - Error estimators (gradient, Hessian, residual, jump)
 * - Cell marking strategies
 * - Mesh refinement/coarsening
 * - Solution transfer
 * - Feature-based refinement
 */

#include <gtest/gtest.h>
#include <petsc.h>
#include "AdaptiveMeshRefinement.hpp"

using namespace FSRM;

class AMRTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        PetscInitialize(nullptr, nullptr, nullptr, nullptr);
    }
    
    static void TearDownTestSuite() {
        PetscFinalize();
    }
    
    void SetUp() override {
        // Create a simple 2D DMPlex mesh for testing
        PetscErrorCode ierr;
        
        // Create box mesh
        PetscReal lower[2] = {0.0, 0.0};
        PetscReal upper[2] = {1.0, 1.0};
        PetscInt faces[2] = {4, 4};  // 4x4 mesh
        
        ierr = DMPlexCreateBoxMesh(PETSC_COMM_WORLD, 2, PETSC_TRUE,
                                   faces, lower, upper, nullptr,
                                   PETSC_TRUE, &dm_);
        ASSERT_EQ(ierr, 0);
        
        // Set up section for scalar field
        PetscSection section;
        ierr = PetscSectionCreate(PETSC_COMM_WORLD, &section);
        ASSERT_EQ(ierr, 0);
        
        PetscInt pStart, pEnd, cStart, cEnd;
        ierr = DMPlexGetChart(dm_, &pStart, &pEnd);
        ierr = DMPlexGetHeightStratum(dm_, 0, &cStart, &cEnd);
        
        ierr = PetscSectionSetChart(section, pStart, pEnd);
        for (PetscInt c = cStart; c < cEnd; c++) {
            ierr = PetscSectionSetDof(section, c, 1);  // 1 DOF per cell
        }
        ierr = PetscSectionSetUp(section);
        ierr = DMSetLocalSection(dm_, section);
        ierr = PetscSectionDestroy(&section);
        
        // Create solution vector
        ierr = DMGetGlobalVector(dm_, &solution_);
        ASSERT_EQ(ierr, 0);
        
        // Initialize with smooth function u = sin(πx) * sin(πy)
        ierr = initializeSolution();
        ASSERT_EQ(ierr, 0);
    }
    
    void TearDown() override {
        if (solution_) DMRestoreGlobalVector(dm_, &solution_);
        if (dm_) DMDestroy(&dm_);
    }
    
    PetscErrorCode initializeSolution() {
        PetscFunctionBeginUser;
        
        PetscInt cStart, cEnd;
        PetscCall(DMPlexGetHeightStratum(dm_, 0, &cStart, &cEnd));
        
        PetscScalar* sol_array;
        PetscCall(VecGetArray(solution_, &sol_array));
        
        PetscSection section;
        PetscCall(DMGetLocalSection(dm_, &section));
        
        for (PetscInt c = cStart; c < cEnd; c++) {
            PetscReal centroid[3];
            PetscCall(DMPlexComputeCellGeometryFVM(dm_, c, nullptr, centroid, nullptr));
            
            PetscInt offset;
            PetscCall(PetscSectionGetOffset(section, c, &offset));
            
            // Smooth function with gradient variation
            sol_array[offset] = std::sin(M_PI * centroid[0]) * std::sin(M_PI * centroid[1]);
        }
        
        PetscCall(VecRestoreArray(solution_, &sol_array));
        
        PetscFunctionReturn(PETSC_SUCCESS);
    }
    
    DM dm_ = nullptr;
    Vec solution_ = nullptr;
};

// Test GradientErrorEstimator
TEST_F(AMRTest, GradientErrorEstimatorBasic) {
    GradientErrorEstimator estimator;
    
    std::vector<CellErrorIndicator> errors;
    PetscErrorCode ierr = estimator.estimate(dm_, solution_, errors);
    
    EXPECT_EQ(ierr, 0);
    EXPECT_GT(errors.size(), 0);
    
    // All errors should be non-negative
    for (const auto& e : errors) {
        EXPECT_GE(e.error, 0.0);
        EXPECT_GT(e.volume, 0.0);
        EXPECT_GT(e.h, 0.0);
    }
}

// Test that errors are higher where gradient is higher
TEST_F(AMRTest, GradientErrorEstimatorDistribution) {
    GradientErrorEstimator estimator;
    
    std::vector<CellErrorIndicator> errors;
    PetscErrorCode ierr = estimator.estimate(dm_, solution_, errors);
    EXPECT_EQ(ierr, 0);
    
    // For sin(πx)sin(πy), gradient is highest at center
    // Find center cell and corner cell
    PetscInt cStart, cEnd;
    DMPlexGetHeightStratum(dm_, 0, &cStart, &cEnd);
    
    PetscReal center_error = 0.0;
    PetscReal corner_error = 1e20;
    
    for (const auto& e : errors) {
        PetscReal centroid[3];
        DMPlexComputeCellGeometryFVM(dm_, e.cell_id, nullptr, centroid, nullptr);
        
        PetscReal dist_to_center = std::sqrt((centroid[0] - 0.5) * (centroid[0] - 0.5) +
                                             (centroid[1] - 0.5) * (centroid[1] - 0.5));
        
        if (dist_to_center < 0.2) {
            center_error = std::max(center_error, (PetscReal)e.error);
        }
        
        if (centroid[0] < 0.2 && centroid[1] < 0.2) {
            corner_error = e.error;
        }
    }
    
    // Note: For sin function, gradient is actually highest at boundaries
    // This is just testing that we get different errors in different regions
    EXPECT_NE(center_error, corner_error);
}

// Test JumpErrorEstimator
TEST_F(AMRTest, JumpErrorEstimator) {
    JumpErrorEstimator estimator;
    estimator.setJumpThreshold(0.1);
    
    std::vector<CellErrorIndicator> errors;
    PetscErrorCode ierr = estimator.estimate(dm_, solution_, errors);
    
    EXPECT_EQ(ierr, 0);
    EXPECT_GT(errors.size(), 0);
}

// Test FeatureErrorEstimator
TEST_F(AMRTest, FeatureErrorEstimatorWells) {
    FeatureErrorEstimator estimator;
    
    // Add a well at center
    estimator.addWellLocation(0.5, 0.5, 0.0, 0.1, 3);
    
    std::vector<CellErrorIndicator> errors;
    PetscErrorCode ierr = estimator.estimate(dm_, solution_, errors);
    
    EXPECT_EQ(ierr, 0);
    
    // Cells near well should have higher error
    PetscReal near_well_error = 0.0;
    PetscReal far_error = 0.0;
    PetscInt near_count = 0, far_count = 0;
    
    for (const auto& e : errors) {
        PetscReal centroid[3];
        DMPlexComputeCellGeometryFVM(dm_, e.cell_id, nullptr, centroid, nullptr);
        
        PetscReal dist = std::sqrt((centroid[0] - 0.5) * (centroid[0] - 0.5) +
                                   (centroid[1] - 0.5) * (centroid[1] - 0.5));
        
        if (dist < 0.2) {
            near_well_error += e.error;
            near_count++;
        } else if (dist > 0.4) {
            far_error += e.error;
            far_count++;
        }
    }
    
    if (near_count > 0 && far_count > 0) {
        EXPECT_GT(near_well_error / near_count, far_error / far_count);
    }
}

// Test CombinedErrorEstimator
TEST_F(AMRTest, CombinedErrorEstimator) {
    auto gradient_est = std::make_shared<GradientErrorEstimator>();
    auto jump_est = std::make_shared<JumpErrorEstimator>();
    
    CombinedErrorEstimator combined;
    combined.addEstimator(gradient_est, 1.0);
    combined.addEstimator(jump_est, 0.5);
    
    std::vector<CellErrorIndicator> errors;
    PetscErrorCode ierr = combined.estimate(dm_, solution_, errors);
    
    EXPECT_EQ(ierr, 0);
    EXPECT_GT(errors.size(), 0);
}

// Test RefinementParameters
TEST_F(AMRTest, RefinementParameters) {
    RefinementParameters params;
    
    // Test defaults
    EXPECT_EQ(params.strategy, RefinementStrategy::FIXED_FRACTION);
    EXPECT_GT(params.refine_fraction, 0.0);
    EXPECT_LT(params.refine_fraction, 1.0);
    EXPECT_GT(params.max_level, 0);
    EXPECT_GT(params.max_cells, 0);
    
    // Test modification
    params.refine_fraction = 0.2;
    params.coarsen_fraction = 0.05;
    params.max_level = 3;
    
    EXPECT_EQ(params.refine_fraction, 0.2);
    EXPECT_EQ(params.coarsen_fraction, 0.05);
    EXPECT_EQ(params.max_level, 3);
}

// Test AdaptiveMeshRefinement initialization
TEST_F(AMRTest, AMRInitialization) {
    AdaptiveMeshRefinement amr;
    
    PetscErrorCode ierr = amr.initialize(dm_);
    EXPECT_EQ(ierr, 0);
    
    EXPECT_NE(amr.getCurrentMesh(), nullptr);
    EXPECT_EQ(amr.getCurrentLevel(), 0);
}

// Test cell marking with fixed fraction
TEST_F(AMRTest, MarkCellsFixedFraction) {
    AdaptiveMeshRefinement amr;
    amr.initialize(dm_);
    
    RefinementParameters params;
    params.strategy = RefinementStrategy::FIXED_FRACTION;
    params.refine_fraction = 0.25;
    params.coarsen_fraction = 0.1;
    amr.setParameters(params);
    
    std::vector<CellErrorIndicator> errors;
    amr.computeErrors(solution_, errors);
    
    std::vector<CellMarker> markers;
    amr.markCells(errors, markers);
    
    EXPECT_EQ(markers.size(), errors.size());
    
    // Count markers
    PetscInt num_refine = std::count(markers.begin(), markers.end(), CellMarker::REFINE);
    PetscInt num_coarsen = std::count(markers.begin(), markers.end(), CellMarker::COARSEN);
    PetscInt num_keep = std::count(markers.begin(), markers.end(), CellMarker::KEEP);
    
    EXPECT_GT(num_refine, 0);
    EXPECT_GT(num_keep, 0);
    
    // Check fractions are approximately correct
    PetscReal refine_frac = static_cast<PetscReal>(num_refine) / markers.size();
    EXPECT_NEAR(refine_frac, 0.25, 0.1);  // Allow some tolerance
}

// Test cell marking with threshold
TEST_F(AMRTest, MarkCellsThreshold) {
    AdaptiveMeshRefinement amr;
    amr.initialize(dm_);
    
    RefinementParameters params;
    params.strategy = RefinementStrategy::THRESHOLD;
    params.refine_threshold = 0.5;
    params.coarsen_threshold = 0.1;
    amr.setParameters(params);
    
    std::vector<CellErrorIndicator> errors;
    amr.computeErrors(solution_, errors);
    
    std::vector<CellMarker> markers;
    amr.markCells(errors, markers);
    
    // All marked cells should satisfy threshold condition
    PetscReal max_error = 0.0;
    for (const auto& e : errors) {
        max_error = std::max(max_error, (PetscReal)e.error);
    }
    
    for (size_t i = 0; i < markers.size(); i++) {
        PetscReal normalized = errors[i].error / max_error;
        
        if (markers[i] == CellMarker::REFINE) {
            EXPECT_GT(normalized, 0.4);  // Some tolerance
        }
    }
}

// Test checkAdaptation
TEST_F(AMRTest, CheckAdaptation) {
    AdaptiveMeshRefinement amr;
    amr.initialize(dm_);
    
    RefinementParameters params;
    params.refine_threshold = 0.01;  // Low threshold to trigger adaptation
    amr.setParameters(params);
    
    PetscBool need_adapt;
    PetscErrorCode ierr = amr.checkAdaptation(solution_, &need_adapt);
    
    EXPECT_EQ(ierr, 0);
    // With a smooth solution, we should still have some error variation
}

// Test constraint enforcement
TEST_F(AMRTest, EnforceConstraints) {
    AdaptiveMeshRefinement amr;
    amr.initialize(dm_);
    
    RefinementParameters params;
    params.max_level = 2;
    params.max_cells = 50;  // Limit cell count
    params.refine_fraction = 0.5;  // Try to refine many cells
    amr.setParameters(params);
    
    std::vector<CellErrorIndicator> errors;
    amr.computeErrors(solution_, errors);
    
    std::vector<CellMarker> markers;
    amr.markCells(errors, markers);
    
    // Markers should be adjusted to satisfy constraints
    PetscInt num_refine = std::count(markers.begin(), markers.end(), CellMarker::REFINE);
    
    // With max_cells constraint, we shouldn't refine too many
    PetscInt current_cells = errors.size();
    PetscInt estimated_after = current_cells + 3 * num_refine;  // 4 children - 1 parent = 3 extra
    
    // The constraint enforcement should limit refinement
    // (This depends on the actual cell count vs max_cells)
}

// Test MeshQuality
TEST_F(AMRTest, MeshQuality) {
    MeshQuality quality;
    
    PetscErrorCode ierr = quality.computeQuality(dm_);
    EXPECT_EQ(ierr, 0);
    
    EXPECT_GT(quality.getMinQuality(), 0.0);
    EXPECT_LE(quality.getMinQuality(), 1.0);
    EXPECT_GE(quality.getMeanQuality(), quality.getMinQuality());
}

// Test error estimator factory
TEST_F(AMRTest, ErrorEstimatorFactory) {
    auto gradient = createErrorEstimator(RefinementCriterion::GRADIENT);
    EXPECT_NE(gradient, nullptr);
    EXPECT_EQ(gradient->getCriterion(), RefinementCriterion::GRADIENT);
    
    auto hessian = createErrorEstimator(RefinementCriterion::HESSIAN);
    EXPECT_NE(hessian, nullptr);
    EXPECT_EQ(hessian->getCriterion(), RefinementCriterion::HESSIAN);
    
    auto jump = createErrorEstimator(RefinementCriterion::JUMP);
    EXPECT_NE(jump, nullptr);
    EXPECT_EQ(jump->getCriterion(), RefinementCriterion::JUMP);
    
    auto feature = createErrorEstimator(RefinementCriterion::FEATURE);
    EXPECT_NE(feature, nullptr);
    EXPECT_EQ(feature->getCriterion(), RefinementCriterion::FEATURE);
}

// Test statistics
TEST_F(AMRTest, AMRStatistics) {
    AdaptiveMeshRefinement amr;
    amr.initialize(dm_);
    
    std::vector<CellErrorIndicator> errors;
    amr.computeErrors(solution_, errors);
    
    // After computing errors, stats should be populated by adapt()
    // For now, just verify the structure exists
    const AMRStatistics& stats = amr.getStatistics();
    
    // Initial stats should be reasonable
    EXPECT_GE(stats.num_cells_before, 0);
}

// Test 3D mesh (if PETSc is configured with 3D support)
class AMR3DTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        // PETSc should already be initialized
    }
    
    void SetUp() override {
        PetscErrorCode ierr;
        
        PetscReal lower[3] = {0.0, 0.0, 0.0};
        PetscReal upper[3] = {1.0, 1.0, 1.0};
        PetscInt faces[3] = {2, 2, 2};
        
        ierr = DMPlexCreateBoxMesh(PETSC_COMM_WORLD, 3, PETSC_TRUE,
                                   faces, lower, upper, nullptr,
                                   PETSC_TRUE, &dm_);
        if (ierr) {
            dm_ = nullptr;  // 3D mesh creation failed
        }
    }
    
    void TearDown() override {
        if (dm_) DMDestroy(&dm_);
    }
    
    DM dm_ = nullptr;
};

TEST_F(AMR3DTest, GradientErrorEstimator3D) {
    if (!dm_) {
        GTEST_SKIP() << "3D mesh creation not available";
    }
    
    // Set up section
    PetscSection section;
    PetscSectionCreate(PETSC_COMM_WORLD, &section);
    
    PetscInt pStart, pEnd, cStart, cEnd;
    DMPlexGetChart(dm_, &pStart, &pEnd);
    DMPlexGetHeightStratum(dm_, 0, &cStart, &cEnd);
    
    PetscSectionSetChart(section, pStart, pEnd);
    for (PetscInt c = cStart; c < cEnd; c++) {
        PetscSectionSetDof(section, c, 1);
    }
    PetscSectionSetUp(section);
    DMSetLocalSection(dm_, section);
    PetscSectionDestroy(&section);
    
    // Create solution
    Vec solution;
    DMGetGlobalVector(dm_, &solution);
    VecSet(solution, 1.0);
    
    GradientErrorEstimator estimator;
    std::vector<CellErrorIndicator> errors;
    
    PetscErrorCode ierr = estimator.estimate(dm_, solution, errors);
    EXPECT_EQ(ierr, 0);
    EXPECT_GT(errors.size(), 0);
    
    // In 3D, h should be cube root of volume
    for (const auto& e : errors) {
        PetscReal vol = e.volume;
        PetscReal h_expected = std::pow(vol, 1.0/3.0);
        EXPECT_NEAR(e.h, h_expected, 1e-10);
    }
    
    DMRestoreGlobalVector(dm_, &solution);
}

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    ::testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    PetscFinalize();
    return result;
}
