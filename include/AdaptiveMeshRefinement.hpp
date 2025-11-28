#ifndef ADAPTIVE_MESH_REFINEMENT_HPP
#define ADAPTIVE_MESH_REFINEMENT_HPP

/**
 * @file AdaptiveMeshRefinement.hpp
 * @brief Adaptive Mesh Refinement (AMR) for unstructured grids using PETSc
 * 
 * Implements h-adaptivity using PETSc's DMPlex and DMForest:
 * - Local refinement based on error estimators
 * - Coarsening for efficiency
 * - Solution transfer between meshes
 * - Dynamic load balancing
 * 
 * Supports multiple refinement criteria:
 * - Gradient-based (pressure, saturation, velocity)
 * - A posteriori error estimators
 * - Feature-based (wells, faults, fronts)
 * - Physics-based (phase boundaries, shocks)
 * 
 * PETSc Components Used:
 * - DMPlex: Unstructured mesh management
 * - DMForest: Octree/quadtree-based AMR
 * - DMPlexAdapt: Metric-based mesh adaptation
 * - VecScatter: Solution transfer
 */

#include <petsc.h>
#include <petscdm.h>
#include <petscdmplex.h>
#include <petscdmforest.h>
#include <petscds.h>
#include <petscsnes.h>

#include <vector>
#include <memory>
#include <string>
#include <map>
#include <functional>
#include <array>

namespace FSRM {

/**
 * @brief Refinement criterion type
 */
enum class RefinementCriterion {
    GRADIENT,           ///< Gradient-based (Kelly error estimator)
    HESSIAN,            ///< Hessian-based for smooth solutions
    RESIDUAL,           ///< Residual-based error estimator
    RECOVERY,           ///< ZZ recovery-based estimator
    JUMP,               ///< Jump-based (interface tracking)
    FEATURE,            ///< Feature-based (wells, faults)
    PHYSICS,            ///< Physics-based (phase fronts)
    COMBINED,           ///< Weighted combination
    CUSTOM              ///< User-defined criterion
};

/**
 * @brief Refinement strategy
 */
enum class RefinementStrategy {
    FIXED_FRACTION,     ///< Refine/coarsen fixed fraction of cells
    FIXED_NUMBER,       ///< Refine/coarsen fixed number of cells
    THRESHOLD,          ///< Refine if error > threshold
    EQUILIBRATION,      ///< Equilibrate error across cells
    OPTIMIZED           ///< Minimize error for given DOF count
};

/**
 * @brief Mesh adaptation method
 */
enum class AdaptationMethod {
    PLEX_REFINE,        ///< DMPlex uniform refinement
    PLEX_ADAPT,         ///< DMPlexAdapt with metric
    FOREST_P4EST,       ///< p4est forest (2D/3D octree)
    FOREST_PUMI,        ///< PUMI/SCOREC mesh adaptation
    PRAGMATIC,          ///< PRAgMaTIc anisotropic adaptation
    MMG,                ///< MMG remeshing library
    CUSTOM              ///< User-provided adaptation
};

/**
 * @brief Cell marking for refinement/coarsening
 */
enum class CellMarker {
    COARSEN = -1,       ///< Mark for coarsening
    KEEP = 0,           ///< No change
    REFINE = 1          ///< Mark for refinement
};

/**
 * @brief Error indicator for a single cell
 */
struct CellErrorIndicator {
    PetscInt cell_id;           ///< Cell index
    PetscReal error;            ///< Error estimate
    PetscReal volume;           ///< Cell volume
    PetscReal h;                ///< Characteristic length
    CellMarker marker;          ///< Refinement decision
    
    // Component-wise errors (for multi-physics)
    PetscReal pressure_error;
    PetscReal saturation_error;
    PetscReal velocity_error;
    PetscReal temperature_error;
};

/**
 * @brief Refinement parameters
 */
struct RefinementParameters {
    // Strategy parameters
    RefinementStrategy strategy = RefinementStrategy::FIXED_FRACTION;
    PetscReal refine_fraction = 0.3;    ///< Fraction to refine
    PetscReal coarsen_fraction = 0.1;   ///< Fraction to coarsen
    PetscReal refine_threshold = 0.5;   ///< Threshold for THRESHOLD strategy
    PetscReal coarsen_threshold = 0.1;  ///< Threshold for coarsening
    
    // Mesh limits
    PetscInt max_level = 5;             ///< Maximum refinement level
    PetscInt min_level = 0;             ///< Minimum refinement level
    PetscInt max_cells = 1000000;       ///< Maximum number of cells
    PetscInt min_cells = 100;           ///< Minimum number of cells
    PetscReal min_cell_size = 1e-6;     ///< Minimum cell size
    PetscReal max_cell_size = 1e6;      ///< Maximum cell size
    
    // Quality constraints
    PetscReal min_quality = 0.1;        ///< Minimum cell quality
    PetscReal max_aspect_ratio = 10.0;  ///< Maximum aspect ratio
    
    // Feature preservation
    PetscBool preserve_boundaries = PETSC_TRUE;
    PetscBool preserve_wells = PETSC_TRUE;
    PetscBool preserve_faults = PETSC_TRUE;
    
    // Buffer zones
    PetscInt buffer_layers = 1;         ///< Layers around refined cells
    
    // Weights for combined criterion
    PetscReal weight_pressure = 1.0;
    PetscReal weight_saturation = 1.0;
    PetscReal weight_velocity = 0.5;
    PetscReal weight_temperature = 0.5;
};

/**
 * @brief AMR statistics
 */
struct AMRStatistics {
    PetscInt num_cells_before;
    PetscInt num_cells_after;
    PetscInt num_refined;
    PetscInt num_coarsened;
    PetscInt current_max_level;
    PetscReal max_error;
    PetscReal min_error;
    PetscReal mean_error;
    PetscReal total_error;
    PetscReal adaptation_time;
    PetscReal transfer_time;
    PetscReal balance_time;
};

/**
 * @brief Forward declarations
 */
class ErrorEstimator;
class SolutionTransfer;
class MeshQuality;

/**
 * @brief Base class for error estimators
 */
class ErrorEstimator {
public:
    ErrorEstimator();
    virtual ~ErrorEstimator() = default;
    
    /**
     * @brief Estimate error for all cells
     * 
     * @param dm Distributed mesh
     * @param solution Solution vector
     * @param errors Output error indicators
     * @return PETSc error code
     */
    virtual PetscErrorCode estimate(DM dm, Vec solution,
                                    std::vector<CellErrorIndicator>& errors) = 0;
    
    /**
     * @brief Get the criterion type
     */
    virtual RefinementCriterion getCriterion() const = 0;
    
    // Configuration
    void setFieldIndex(PetscInt idx) { field_index_ = idx; }
    void setComponentIndex(PetscInt idx) { component_index_ = idx; }
    
protected:
    PetscInt field_index_ = 0;
    PetscInt component_index_ = 0;
};

/**
 * @brief Gradient-based error estimator (Kelly estimator)
 * 
 * η_K² = h_K * Σ_{e∈∂K} ||[∇u · n]||²_{L²(e)}
 * 
 * Measures jump in normal gradient across cell faces.
 */
class GradientErrorEstimator : public ErrorEstimator {
public:
    GradientErrorEstimator();
    
    PetscErrorCode estimate(DM dm, Vec solution,
                           std::vector<CellErrorIndicator>& errors) override;
    
    RefinementCriterion getCriterion() const override { 
        return RefinementCriterion::GRADIENT; 
    }
    
    // Power for h scaling (typically 0.5 for linear elements)
    void setHPower(PetscReal p) { h_power_ = p; }
    
private:
    PetscReal h_power_ = 0.5;
    
    PetscErrorCode computeGradientJump(DM dm, Vec solution, PetscInt cell,
                                       PetscReal* jump);
};

/**
 * @brief Hessian-based error estimator
 * 
 * Uses recovered Hessian for smooth solutions.
 * Better for elliptic problems without shocks.
 */
class HessianErrorEstimator : public ErrorEstimator {
public:
    HessianErrorEstimator();
    
    PetscErrorCode estimate(DM dm, Vec solution,
                           std::vector<CellErrorIndicator>& errors) override;
    
    RefinementCriterion getCriterion() const override { 
        return RefinementCriterion::HESSIAN; 
    }
    
    // Recovery method
    enum class RecoveryMethod { SPR, L2_PROJECTION, PATCH_RECOVERY };
    void setRecoveryMethod(RecoveryMethod method) { recovery_method_ = method; }
    
private:
    RecoveryMethod recovery_method_ = RecoveryMethod::SPR;
    
    PetscErrorCode recoverHessian(DM dm, Vec solution, Vec hessian);
};

/**
 * @brief Residual-based error estimator
 * 
 * η_K² = h_K² ||R(u_h)||²_{L²(K)} + h_K Σ_{e∈∂K} ||[flux]||²_{L²(e)}
 */
class ResidualErrorEstimator : public ErrorEstimator {
public:
    ResidualErrorEstimator();
    
    PetscErrorCode estimate(DM dm, Vec solution,
                           std::vector<CellErrorIndicator>& errors) override;
    
    RefinementCriterion getCriterion() const override { 
        return RefinementCriterion::RESIDUAL; 
    }
    
    // Set residual function
    using ResidualFunction = std::function<PetscErrorCode(DM, Vec, Vec)>;
    void setResidualFunction(ResidualFunction func) { residual_func_ = func; }
    
private:
    ResidualFunction residual_func_;
};

/**
 * @brief Jump-based error estimator for interfaces
 * 
 * Tracks saturation fronts, phase boundaries.
 */
class JumpErrorEstimator : public ErrorEstimator {
public:
    JumpErrorEstimator();
    
    PetscErrorCode estimate(DM dm, Vec solution,
                           std::vector<CellErrorIndicator>& errors) override;
    
    RefinementCriterion getCriterion() const override { 
        return RefinementCriterion::JUMP; 
    }
    
    void setJumpThreshold(PetscReal thresh) { jump_threshold_ = thresh; }
    
private:
    PetscReal jump_threshold_ = 0.1;
};

/**
 * @brief Feature-based refinement around wells and faults
 */
class FeatureErrorEstimator : public ErrorEstimator {
public:
    FeatureErrorEstimator();
    
    PetscErrorCode estimate(DM dm, Vec solution,
                           std::vector<CellErrorIndicator>& errors) override;
    
    RefinementCriterion getCriterion() const override { 
        return RefinementCriterion::FEATURE; 
    }
    
    // Add features requiring refinement
    void addWellLocation(PetscReal x, PetscReal y, PetscReal z, 
                         PetscReal radius, PetscInt level);
    void addFaultTrace(const std::vector<std::array<PetscReal, 3>>& points,
                       PetscReal width, PetscInt level);
    void addRefinementBox(PetscReal xmin, PetscReal xmax,
                          PetscReal ymin, PetscReal ymax,
                          PetscReal zmin, PetscReal zmax, PetscInt level);
    
private:
    struct WellFeature {
        PetscReal x, y, z, radius;
        PetscInt level;
    };
    
    struct FaultFeature {
        std::vector<std::array<PetscReal, 3>> trace;
        PetscReal width;
        PetscInt level;
    };
    
    struct BoxFeature {
        PetscReal bounds[6];  // xmin, xmax, ymin, ymax, zmin, zmax
        PetscInt level;
    };
    
    std::vector<WellFeature> wells_;
    std::vector<FaultFeature> faults_;
    std::vector<BoxFeature> boxes_;
    
    PetscReal computeFeatureDistance(PetscReal x, PetscReal y, PetscReal z) const;
};

/**
 * @brief Combined multi-physics error estimator
 */
class CombinedErrorEstimator : public ErrorEstimator {
public:
    CombinedErrorEstimator();
    
    PetscErrorCode estimate(DM dm, Vec solution,
                           std::vector<CellErrorIndicator>& errors) override;
    
    RefinementCriterion getCriterion() const override { 
        return RefinementCriterion::COMBINED; 
    }
    
    void addEstimator(std::shared_ptr<ErrorEstimator> est, PetscReal weight);
    void setWeights(const std::vector<PetscReal>& weights);
    
private:
    std::vector<std::shared_ptr<ErrorEstimator>> estimators_;
    std::vector<PetscReal> weights_;
};

/**
 * @brief Solution transfer between meshes
 */
class SolutionTransfer {
public:
    SolutionTransfer();
    ~SolutionTransfer();
    
    /**
     * @brief Transfer solution from old mesh to new mesh
     * 
     * @param dm_old Old distributed mesh
     * @param dm_new New adapted mesh
     * @param solution_old Solution on old mesh
     * @param solution_new Solution on new mesh (output)
     * @return PETSc error code
     */
    PetscErrorCode transfer(DM dm_old, DM dm_new,
                            Vec solution_old, Vec solution_new);
    
    /**
     * @brief Transfer multiple fields
     */
    PetscErrorCode transferFields(DM dm_old, DM dm_new,
                                  const std::vector<Vec>& fields_old,
                                  std::vector<Vec>& fields_new);
    
    // Transfer method
    enum class TransferMethod {
        INTERPOLATION,      ///< Finite element interpolation
        PROJECTION,         ///< L2 projection
        INJECTION,          ///< Direct injection (coarsening)
        CONSERVATIVE        ///< Mass-conserving transfer
    };
    
    void setMethod(TransferMethod method) { method_ = method; }
    void setConservative(PetscBool cons) { conservative_ = cons; }
    
private:
    TransferMethod method_ = TransferMethod::INTERPOLATION;
    PetscBool conservative_ = PETSC_TRUE;
    
    Mat transfer_matrix_ = nullptr;
    VecScatter scatter_ = nullptr;
    
    PetscErrorCode buildTransferMatrix(DM dm_old, DM dm_new);
    PetscErrorCode conservativeTransfer(DM dm_old, DM dm_new,
                                        Vec solution_old, Vec solution_new);
};

/**
 * @brief Mesh quality metrics
 */
class MeshQuality {
public:
    MeshQuality();
    
    /**
     * @brief Compute quality metrics for all cells
     */
    PetscErrorCode computeQuality(DM dm);
    
    /**
     * @brief Check if mesh satisfies quality constraints
     */
    PetscBool isAcceptable(const RefinementParameters& params) const;
    
    // Get statistics
    PetscReal getMinQuality() const { return min_quality_; }
    PetscReal getMeanQuality() const { return mean_quality_; }
    PetscReal getMaxAspectRatio() const { return max_aspect_ratio_; }
    PetscInt getNumBadCells() const { return num_bad_cells_; }
    
    // Quality measure type
    enum class QualityMeasure {
        SCALED_JACOBIAN,    ///< Scaled Jacobian (0 to 1)
        SHAPE,              ///< Shape quality
        ASPECT_RATIO,       ///< Aspect ratio
        SKEWNESS,           ///< Skewness
        CONDITION_NUMBER    ///< Condition number of Jacobian
    };
    
    void setQualityMeasure(QualityMeasure measure) { measure_ = measure; }
    
private:
    QualityMeasure measure_ = QualityMeasure::SCALED_JACOBIAN;
    
    PetscReal min_quality_ = 1.0;
    PetscReal mean_quality_ = 1.0;
    PetscReal max_aspect_ratio_ = 1.0;
    PetscInt num_bad_cells_ = 0;
    
    std::vector<PetscReal> cell_quality_;
    
    PetscReal computeCellQuality(DM dm, PetscInt cell) const;
};

/**
 * @brief Main AMR controller class
 * 
 * Coordinates error estimation, mesh adaptation, and solution transfer.
 */
class AdaptiveMeshRefinement {
public:
    AdaptiveMeshRefinement();
    ~AdaptiveMeshRefinement();
    
    /**
     * @brief Initialize AMR with mesh
     * 
     * @param dm Initial distributed mesh (DMPlex or DMForest)
     * @return PETSc error code
     */
    PetscErrorCode initialize(DM dm);
    
    /**
     * @brief Perform mesh adaptation
     * 
     * @param solution Current solution vector
     * @param dm_new Output: new adapted mesh
     * @param solution_new Output: solution on new mesh
     * @return PETSc error code
     */
    PetscErrorCode adapt(Vec solution, DM* dm_new, Vec* solution_new);
    
    /**
     * @brief Check if adaptation is needed
     * 
     * @param solution Current solution
     * @param need_adapt Output: whether adaptation is needed
     * @return PETSc error code
     */
    PetscErrorCode checkAdaptation(Vec solution, PetscBool* need_adapt);
    
    /**
     * @brief Compute error indicators without adapting
     */
    PetscErrorCode computeErrors(Vec solution,
                                 std::vector<CellErrorIndicator>& errors);
    
    /**
     * @brief Mark cells for refinement/coarsening
     */
    PetscErrorCode markCells(const std::vector<CellErrorIndicator>& errors,
                             std::vector<CellMarker>& markers);
    
    /**
     * @brief Refine marked cells
     */
    PetscErrorCode refineMesh(const std::vector<CellMarker>& markers,
                              DM* dm_new);
    
    /**
     * @brief Transfer solution to new mesh
     */
    PetscErrorCode transferSolution(DM dm_old, DM dm_new,
                                    Vec solution_old, Vec* solution_new);
    
    /**
     * @brief Rebalance mesh across processors
     */
    PetscErrorCode rebalance(DM dm);
    
    // Configuration
    void setParameters(const RefinementParameters& params) { params_ = params; }
    RefinementParameters& getParameters() { return params_; }
    
    void setMethod(AdaptationMethod method) { method_ = method; }
    void setCriterion(RefinementCriterion criterion);
    void setErrorEstimator(std::shared_ptr<ErrorEstimator> estimator);
    
    // Add feature-based refinement
    void addWellRefinement(PetscReal x, PetscReal y, PetscReal z,
                           PetscReal radius, PetscInt level);
    void addFaultRefinement(const std::vector<std::array<PetscReal, 3>>& trace,
                            PetscReal width, PetscInt level);
    
    // Callbacks
    using PreAdaptCallback = std::function<PetscErrorCode(DM, Vec)>;
    using PostAdaptCallback = std::function<PetscErrorCode(DM, DM, Vec, Vec)>;
    
    void setPreAdaptCallback(PreAdaptCallback cb) { pre_adapt_cb_ = cb; }
    void setPostAdaptCallback(PostAdaptCallback cb) { post_adapt_cb_ = cb; }
    
    // Statistics
    const AMRStatistics& getStatistics() const { return stats_; }
    void printStatistics() const;
    
    // Access current mesh
    DM getCurrentMesh() const { return dm_; }
    PetscInt getCurrentLevel() const { return current_level_; }
    
    // Output
    PetscErrorCode writeAdaptedMesh(const std::string& filename) const;
    PetscErrorCode writeErrorField(const std::string& filename) const;
    
private:
    DM dm_ = nullptr;                   ///< Current mesh
    DM dm_coarse_ = nullptr;            ///< Coarsest mesh (for forest)
    
    RefinementParameters params_;
    AdaptationMethod method_ = AdaptationMethod::PLEX_ADAPT;
    
    std::shared_ptr<ErrorEstimator> error_estimator_;
    std::unique_ptr<SolutionTransfer> solution_transfer_;
    std::unique_ptr<MeshQuality> mesh_quality_;
    
    std::vector<CellErrorIndicator> errors_;
    Vec error_vec_ = nullptr;           ///< Error field for output
    
    AMRStatistics stats_;
    PetscInt current_level_ = 0;
    PetscInt adapt_count_ = 0;
    
    PreAdaptCallback pre_adapt_cb_;
    PostAdaptCallback post_adapt_cb_;
    
    // Private methods
    PetscErrorCode createMetricFromErrors(const std::vector<CellErrorIndicator>& errors,
                                          Vec* metric);
    PetscErrorCode applyPlexAdapt(Vec metric, DM* dm_new);
    PetscErrorCode applyForestAdapt(const std::vector<CellMarker>& markers,
                                    DM* dm_new);
    PetscErrorCode applyUniformRefine(DM* dm_new);
    PetscErrorCode enforceConstraints(std::vector<CellMarker>& markers);
    PetscErrorCode addBufferLayers(std::vector<CellMarker>& markers);
    
    // Marking strategies
    PetscErrorCode markFixedFraction(const std::vector<CellErrorIndicator>& errors,
                                     std::vector<CellMarker>& markers);
    PetscErrorCode markThreshold(const std::vector<CellErrorIndicator>& errors,
                                 std::vector<CellMarker>& markers);
    PetscErrorCode markEquilibration(const std::vector<CellErrorIndicator>& errors,
                                     std::vector<CellMarker>& markers);
};

/**
 * @brief Factory function to create error estimator
 */
std::shared_ptr<ErrorEstimator> createErrorEstimator(
    RefinementCriterion criterion,
    const std::map<std::string, PetscReal>& params = {});

/**
 * @brief Configure AMR from options
 */
PetscErrorCode AMRSetFromOptions(AdaptiveMeshRefinement* amr);

/**
 * @brief Create AMR-enabled DM from existing mesh
 */
PetscErrorCode DMPlexCreateAdaptive(DM dm_in, AdaptationMethod method, DM* dm_out);

} // namespace FSRM

#endif // ADAPTIVE_MESH_REFINEMENT_HPP
