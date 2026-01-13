/**
 * @file NeuralAMR.hpp
 * @brief Neural Network-Based Adaptive Mesh Refinement
 * 
 * Learns optimal mesh refinement strategies:
 * 
 * - Neural error estimator for refinement indicators
 * - Learned refinement criteria
 * - Predictive mesh adaptation
 * - Multi-physics refinement
 * 
 * Benefits:
 * - Anticipatory refinement before error develops
 * - Physics-aware adaptation
 * - Reduced over-refinement
 * - Feature-tracking capabilities
 */

#ifndef FSRM_NEURAL_AMR_HPP
#define FSRM_NEURAL_AMR_HPP

#include "FourierNeuralOperator.hpp"
#include "NeuralSurrogates.hpp"
#include "AdaptiveMeshRefinement.hpp"
#include <petsc.h>
#include <vector>
#include <memory>

namespace FSRM {
namespace ML {

// =============================================================================
// Neural Error Indicator
// =============================================================================

/**
 * @brief Neural network-based error indicator
 */
class NeuralErrorIndicator {
public:
    struct ErrorIndicatorConfig {
        // Input features
        bool use_solution_values = true;
        bool use_gradient = true;
        bool use_hessian = true;
        bool use_residual = true;
        bool use_geometry = true;
        int neighbor_layers = 2;  // Local stencil size
        
        // Network
        std::vector<int> hidden_layers = {128, 256, 128};
        std::string activation = "relu";
        
        // Output
        enum class OutputType {
            BINARY,          // Refine/don't refine
            REFINEMENT_LEVEL, // 0, 1, 2, ... levels
            CONTINUOUS        // Continuous error estimate
        };
        OutputType output_type = OutputType::CONTINUOUS;
        
        // Training
        double learning_rate = 1e-3;
        int epochs = 200;
    };
    
    NeuralErrorIndicator(const ErrorIndicatorConfig& config);
    
    /**
     * @brief Compute error indicator for each cell
     */
    Tensor computeIndicator(DM dm, Vec solution);
    
    /**
     * @brief Get cells to refine based on indicator
     */
    std::vector<PetscInt> getCellsToRefine(DM dm, Vec solution,
                                           double threshold);
    
    /**
     * @brief Get cells to coarsen
     */
    std::vector<PetscInt> getCellsToCoarsen(DM dm, Vec solution,
                                            double threshold);
    
    /**
     * @brief Batch prediction for multiple states
     */
    std::vector<Tensor> computeIndicatorBatch(
        const std::vector<std::pair<DM, Vec>>& states);
    
    /**
     * @brief Train from reference error data
     */
    void train(const std::vector<std::pair<DM, Vec>>& states,
              const std::vector<Tensor>& reference_errors);
    
    /**
     * @brief Train from refinement decisions
     */
    void trainFromDecisions(
        const std::vector<std::pair<DM, Vec>>& before_states,
        const std::vector<std::pair<DM, Vec>>& after_states,
        const std::vector<Tensor>& errors);
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    ErrorIndicatorConfig config_;
    std::unique_ptr<MLPModel> indicator_net_;
    
    // Feature extraction
    Tensor extractCellFeatures(DM dm, Vec solution, PetscInt cell);
    Tensor extractNeighborFeatures(DM dm, Vec solution, PetscInt cell);
};

// =============================================================================
// Predictive Refinement
// =============================================================================

/**
 * @brief Predict future refinement needs
 */
class PredictiveRefinement {
public:
    struct PredictiveConfig {
        // Prediction horizon
        int lookahead_steps = 10;
        double lookahead_time = 1.0;
        
        // Network
        std::vector<int> hidden_layers = {256, 512, 256};
        bool use_temporal_encoding = true;
        
        // Features
        bool use_velocity_field = true;  // For advection
        bool use_physics_predictors = true;
        
        // Training
        double learning_rate = 1e-3;
        int epochs = 200;
    };
    
    PredictiveRefinement(const PredictiveConfig& config);
    
    /**
     * @brief Predict future error locations
     */
    Tensor predictFutureError(DM dm, Vec current_solution,
                             double current_time, double future_time);
    
    /**
     * @brief Get anticipatory refinement
     */
    std::vector<PetscInt> getAnticipatoryRefinement(
        DM dm, Vec current_solution,
        double current_time, double predict_time,
        double threshold);
    
    /**
     * @brief Predict feature trajectory (e.g., shock position)
     */
    std::vector<Tensor> predictFeatureTrajectory(
        DM dm, Vec solution, int num_steps);
    
    /**
     * @brief Train from time evolution data
     */
    void train(const std::vector<std::vector<std::pair<DM, Vec>>>& trajectories,
              const std::vector<std::vector<double>>& times);
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    PredictiveConfig config_;
    
    // Temporal prediction network
    std::unique_ptr<MLPModel> predictor_;
    
    // Feature tracker
    std::unique_ptr<MLPModel> feature_tracker_;
    
    Tensor temporalEncoding(double current_time, double future_time);
};

// =============================================================================
// Physics-Aware Refinement
// =============================================================================

/**
 * @brief Multi-physics aware refinement
 */
class PhysicsAwareRefinement {
public:
    struct PhysicsConfig {
        // Physics types to consider
        bool flow_physics = true;
        bool geomechanics = true;
        bool fracture_physics = true;
        bool thermal = false;
        
        // Per-physics weights
        std::map<std::string, double> physics_weights;
        
        // Coupling awareness
        bool consider_coupling = true;
        
        // Network
        std::vector<int> hidden_layers = {256, 256};
    };
    
    PhysicsAwareRefinement(const PhysicsConfig& config);
    
    /**
     * @brief Compute physics-aware indicator
     */
    Tensor computePhysicsIndicator(DM dm, Vec solution,
                                   const std::string& physics_type);
    
    /**
     * @brief Combine multiple physics indicators
     */
    Tensor combineIndicators(const std::vector<Tensor>& indicators,
                            const std::vector<std::string>& physics_types);
    
    /**
     * @brief Learn optimal combination weights
     */
    void learnCombinationWeights(
        const std::vector<std::vector<Tensor>>& physics_indicators,
        const std::vector<Tensor>& optimal_refinements);
    
    /**
     * @brief Identify coupled regions needing special treatment
     */
    std::vector<PetscInt> identifyCoupledRegions(DM dm, Vec solution);
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    PhysicsConfig config_;
    
    // Per-physics indicator networks
    std::map<std::string, std::unique_ptr<NeuralErrorIndicator>> physics_indicators_;
    
    // Combination network
    std::unique_ptr<MLPModel> combiner_;
    
    // Coupling detector
    std::unique_ptr<MLPModel> coupling_detector_;
};

// =============================================================================
// Refinement Level Predictor
// =============================================================================

/**
 * @brief Predict optimal refinement level per cell
 */
class RefinementLevelPredictor {
public:
    struct LevelConfig {
        int max_refinement_level = 5;
        
        // Network
        std::vector<int> hidden_layers = {128, 128};
        
        // Cost awareness
        bool consider_cost = true;
        double cost_per_level = 8.0;  // Typical 2D: 4x, 3D: 8x
    };
    
    RefinementLevelPredictor(const LevelConfig& config);
    
    /**
     * @brief Predict optimal refinement level
     */
    std::vector<int> predictLevels(DM dm, Vec solution,
                                   double target_error);
    
    /**
     * @brief Predict with cost constraint
     */
    std::vector<int> predictLevelsWithBudget(DM dm, Vec solution,
                                             double target_error,
                                             double max_cells);
    
    /**
     * @brief Train from optimal refinement sequences
     */
    void train(const std::vector<std::pair<DM, Vec>>& states,
              const std::vector<std::vector<int>>& optimal_levels);
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    LevelConfig config_;
    std::unique_ptr<MLPModel> level_predictor_;
    
    // Cost-aware optimization
    std::vector<int> optimizeForBudget(const std::vector<double>& scores,
                                       double budget);
};

// =============================================================================
// Feature Tracker
// =============================================================================

/**
 * @brief Track and follow features for refinement
 */
class FeatureTracker {
public:
    struct TrackerConfig {
        // Feature types
        enum class FeatureType {
            SHOCK,
            INTERFACE,
            FRACTURE_TIP,
            FAULT_SLIP,
            HIGH_GRADIENT,
            CUSTOM
        };
        std::vector<FeatureType> feature_types;
        
        // Tracking
        double detection_threshold = 0.5;
        double tracking_radius = 0.1;
        
        // Network
        std::vector<int> hidden_layers = {256, 256};
    };
    
    FeatureTracker(const TrackerConfig& config);
    
    /**
     * @brief Detect features in current state
     */
    std::vector<Tensor> detectFeatures(DM dm, Vec solution);
    
    /**
     * @brief Track features between time steps
     */
    void track(DM dm_old, Vec sol_old,
              DM dm_new, Vec sol_new, double dt);
    
    /**
     * @brief Get cells near tracked features
     */
    std::vector<PetscInt> getCellsNearFeatures(DM dm, double radius);
    
    /**
     * @brief Predict feature motion for next step
     */
    std::vector<Tensor> predictFeatureMotion(double dt);
    
    /**
     * @brief Train feature detector
     */
    void trainDetector(const std::vector<std::pair<DM, Vec>>& states,
                      const std::vector<std::vector<Tensor>>& feature_labels);
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    TrackerConfig config_;
    
    // Feature detector
    std::unique_ptr<MLPModel> detector_;
    
    // Motion predictor
    std::unique_ptr<MLPModel> motion_predictor_;
    
    // Current tracked features
    std::vector<Tensor> tracked_features_;
    std::vector<Tensor> feature_velocities_;
};

// =============================================================================
// Unified Neural AMR Controller
// =============================================================================

/**
 * @brief Complete neural AMR system
 */
class NeuralAMRController {
public:
    struct AMRConfig {
        // Component selection
        bool use_neural_error = true;
        bool use_predictive = true;
        bool use_physics_aware = true;
        bool use_feature_tracking = true;
        
        // Refinement parameters
        double refinement_threshold = 0.5;
        double coarsening_threshold = 0.1;
        int max_level = 5;
        int min_level = 0;
        
        // Adaptivity frequency
        int adapt_frequency = 10;  // Time steps between adaptation
        
        // Component configs
        NeuralErrorIndicator::ErrorIndicatorConfig error_config;
        PredictiveRefinement::PredictiveConfig predictive_config;
        PhysicsAwareRefinement::PhysicsConfig physics_config;
        FeatureTracker::TrackerConfig tracker_config;
        
        // Budget constraints
        double max_cell_growth = 2.0;
        size_t max_cells = 1000000;
    };
    
    NeuralAMRController(const AMRConfig& config);
    
    /**
     * @brief Check if adaptation should occur
     */
    bool shouldAdapt(int step, double time);
    
    /**
     * @brief Compute combined refinement indicator
     */
    Tensor computeRefinementIndicator(DM dm, Vec solution, double time);
    
    /**
     * @brief Get refinement/coarsening flags
     */
    void getRefinementFlags(DM dm, Vec solution, double time,
                           std::vector<PetscInt>& refine_cells,
                           std::vector<PetscInt>& coarsen_cells);
    
    /**
     * @brief Perform adaptation
     */
    PetscErrorCode adapt(DM& dm, Vec& solution, double time);
    
    /**
     * @brief Train all components
     */
    void train(const std::vector<std::vector<std::pair<DM, Vec>>>& trajectories,
              const std::vector<std::vector<double>>& times,
              const std::vector<std::vector<double>>& errors);
    
    // Statistics
    struct AMRStats {
        int num_adaptations = 0;
        int total_refinements = 0;
        int total_coarsenings = 0;
        double avg_cells = 0.0;
        double max_cells = 0.0;
        double neural_time = 0.0;
    };
    
    const AMRStats& getStats() const { return stats_; }
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    AMRConfig config_;
    AMRStats stats_;
    
    std::unique_ptr<NeuralErrorIndicator> error_indicator_;
    std::unique_ptr<PredictiveRefinement> predictive_;
    std::unique_ptr<PhysicsAwareRefinement> physics_aware_;
    std::unique_ptr<FeatureTracker> feature_tracker_;
    std::unique_ptr<RefinementLevelPredictor> level_predictor_;
    
    int step_counter_ = 0;
    
    // Combine all indicators
    Tensor combineAllIndicators(DM dm, Vec solution, double time);
    
    // Respect budget constraints
    void enforceConstraints(std::vector<PetscInt>& refine_cells,
                           std::vector<PetscInt>& coarsen_cells,
                           size_t current_cells);
};

// =============================================================================
// Configuration Parsing
// =============================================================================

NeuralErrorIndicator::ErrorIndicatorConfig parseNeuralErrorConfig(
    const std::map<std::string, std::string>& config);
NeuralAMRController::AMRConfig parseNeuralAMRConfig(
    const std::map<std::string, std::string>& config);

} // namespace ML
} // namespace FSRM

#endif // FSRM_NEURAL_AMR_HPP
