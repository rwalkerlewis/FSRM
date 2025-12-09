/**
 * @file MultiFidelityLearning.hpp
 * @brief Multi-Fidelity Machine Learning for Scientific Computing
 * 
 * Combines models of different fidelity levels:
 * 
 * - Low-fidelity: Fast but less accurate (coarse mesh, neural surrogates)
 * - High-fidelity: Accurate but expensive (fine mesh, full physics)
 * 
 * Key approaches:
 * - Additive correction: y_HF = y_LF + Δ_NN
 * - Multiplicative correction: y_HF = y_LF * (1 + Δ_NN)
 * - Residual learning: Learn the discrepancy
 * - Transfer learning: Pre-train on LF, fine-tune on HF
 * - Co-Kriging: Gaussian process multi-fidelity
 * 
 * Benefits: 10-100x cost reduction with minimal accuracy loss
 */

#ifndef FSRM_MULTI_FIDELITY_LEARNING_HPP
#define FSRM_MULTI_FIDELITY_LEARNING_HPP

#include "FourierNeuralOperator.hpp"
#include "NeuralSurrogates.hpp"
#include "BayesianNeuralOperator.hpp"
#include <petsc.h>
#include <vector>
#include <memory>
#include <functional>

namespace FSRM {
namespace ML {

// =============================================================================
// Multi-Fidelity Configuration
// =============================================================================

/**
 * @brief Fidelity level specification
 */
struct FidelityLevel {
    int level;                    // 0 = lowest, higher = more accurate
    std::string name;             // e.g., "coarse_mesh", "fine_mesh"
    double relative_cost;         // Relative computational cost
    double expected_accuracy;     // Expected accuracy (0-1)
    
    // Model specification
    enum class ModelType {
        NEURAL_SURROGATE,
        FNO,
        GNO,
        COARSE_SIMULATION,
        FINE_SIMULATION
    };
    ModelType model_type;
    
    // For simulation-based fidelity
    int mesh_refinement = 1;
    double time_step_factor = 1.0;
    bool use_simplified_physics = false;
};

// =============================================================================
// Additive Correction Model
// =============================================================================

/**
 * @brief Learn additive correction: y_HF ≈ y_LF + Δ(x)
 */
class AdditiveCorrection {
public:
    struct AdditiveConfig {
        // Correction network
        std::vector<int> hidden_layers = {256, 256};
        std::string activation = "gelu";
        
        // Input features
        bool include_lf_output = true;    // Include y_LF in correction input
        bool include_x = true;            // Include original input x
        bool include_gradient = false;    // Include LF gradient
        
        // Training
        double learning_rate = 1e-3;
        int epochs = 200;
        int batch_size = 32;
        
        // Regularization
        double l2_reg = 1e-4;
        bool penalize_large_correction = true;
        double correction_penalty = 0.1;
    };
    
    AdditiveCorrection(const AdditiveConfig& config);
    
    /**
     * @brief Apply correction to low-fidelity output
     */
    Tensor correct(const Tensor& lf_output, const Tensor& input);
    
    /**
     * @brief Train correction model
     * 
     * @param inputs Original inputs
     * @param lf_outputs Low-fidelity predictions
     * @param hf_outputs High-fidelity targets
     */
    void train(const std::vector<Tensor>& inputs,
              const std::vector<Tensor>& lf_outputs,
              const std::vector<Tensor>& hf_outputs);
    
    /**
     * @brief Get correction magnitude
     */
    Tensor getCorrectionMagnitude(const Tensor& lf_output, const Tensor& input);
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    AdditiveConfig config_;
    std::unique_ptr<MLPModel> correction_net_;
};

// =============================================================================
// Multiplicative Correction Model
// =============================================================================

/**
 * @brief Learn multiplicative correction: y_HF ≈ y_LF * (1 + Δ(x))
 */
class MultiplicativeCorrection {
public:
    struct MultiplicativeConfig {
        std::vector<int> hidden_layers = {256, 256};
        std::string activation = "tanh";  // Bounded output
        
        // Bounds on correction factor
        double min_factor = 0.1;
        double max_factor = 10.0;
        
        // Training
        double learning_rate = 1e-3;
        int epochs = 200;
    };
    
    MultiplicativeCorrection(const MultiplicativeConfig& config);
    
    /**
     * @brief Apply multiplicative correction
     */
    Tensor correct(const Tensor& lf_output, const Tensor& input);
    
    /**
     * @brief Get correction factor
     */
    Tensor getCorrectionFactor(const Tensor& lf_output, const Tensor& input);
    
    void train(const std::vector<Tensor>& inputs,
              const std::vector<Tensor>& lf_outputs,
              const std::vector<Tensor>& hf_outputs);
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    MultiplicativeConfig config_;
    std::unique_ptr<MLPModel> factor_net_;
};

// =============================================================================
// Residual Learning Model
// =============================================================================

/**
 * @brief Learn residual between fidelity levels
 */
class ResidualFidelityModel {
public:
    struct ResidualConfig {
        // Multiple fidelity levels
        int num_levels = 3;  // e.g., coarse, medium, fine
        
        // Per-level correction networks
        std::vector<int> hidden_layers = {256, 256};
        
        // Cascaded vs parallel
        bool cascaded = true;  // Each level corrects previous
        
        // Training
        double learning_rate = 1e-3;
        int epochs = 200;
    };
    
    ResidualFidelityModel(const ResidualConfig& config);
    
    /**
     * @brief Predict with cascaded corrections
     */
    Tensor predict(const std::vector<Tensor>& level_outputs,
                  const Tensor& input);
    
    /**
     * @brief Predict up to specified fidelity level
     */
    Tensor predictToLevel(const std::vector<Tensor>& level_outputs,
                         const Tensor& input,
                         int target_level);
    
    /**
     * @brief Train all correction models
     */
    void train(const std::vector<std::vector<Tensor>>& level_outputs,
              const std::vector<Tensor>& inputs,
              const std::vector<Tensor>& hf_targets);
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    ResidualConfig config_;
    std::vector<std::unique_ptr<AdditiveCorrection>> corrections_;
};

// =============================================================================
// Transfer Learning Framework
// =============================================================================

/**
 * @brief Transfer learning between fidelity levels
 */
class TransferLearningMF {
public:
    struct TransferConfig {
        // Base model architecture
        std::vector<int> shared_layers = {256, 512};  // Shared across fidelities
        std::vector<int> task_layers = {128};          // Task-specific
        
        // Pre-training on LF
        int lf_pretrain_epochs = 500;
        double lf_learning_rate = 1e-3;
        
        // Fine-tuning on HF
        int hf_finetune_epochs = 100;
        double hf_learning_rate = 1e-4;
        bool freeze_shared = false;  // Freeze shared layers during HF training
        
        // Progressive training
        bool use_progressive = true;  // Train LF -> medium -> HF progressively
    };
    
    TransferLearningMF(const TransferConfig& config);
    
    /**
     * @brief Pre-train on low-fidelity data
     */
    void pretrain(const std::vector<Tensor>& lf_inputs,
                 const std::vector<Tensor>& lf_outputs);
    
    /**
     * @brief Fine-tune on high-fidelity data
     */
    void finetune(const std::vector<Tensor>& hf_inputs,
                 const std::vector<Tensor>& hf_outputs);
    
    /**
     * @brief Progressive training across all levels
     */
    void trainProgressive(
        const std::vector<std::vector<Tensor>>& all_level_inputs,
        const std::vector<std::vector<Tensor>>& all_level_outputs);
    
    /**
     * @brief Predict at highest fidelity
     */
    Tensor predict(const Tensor& input);
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    TransferConfig config_;
    
    // Shared encoder
    std::unique_ptr<MLPModel> shared_encoder_;
    
    // Task-specific heads for each fidelity level
    std::vector<std::unique_ptr<MLPModel>> task_heads_;
};

// =============================================================================
// Multi-Fidelity Neural Operator
// =============================================================================

/**
 * @brief FNO with multi-fidelity training
 */
class MultiFidelityFNO {
public:
    struct MFFNOConfig : FNOConfig {
        // Multi-fidelity settings
        int num_fidelity_levels = 2;
        
        // Per-level mode counts
        std::vector<int> modes_per_level = {4, 8, 12};
        
        // Correction type
        enum class CorrectionType {
            ADDITIVE,
            MULTIPLICATIVE,
            SPECTRAL  // Correct in Fourier space
        };
        CorrectionType correction_type = CorrectionType::ADDITIVE;
        
        // Data budget allocation
        std::vector<int> samples_per_level = {1000, 100, 10};
    };
    
    MultiFidelityFNO(const MFFNOConfig& config);
    
    /**
     * @brief Train with data from multiple fidelity levels
     */
    void train(const std::vector<std::vector<Tensor>>& level_inputs,
              const std::vector<std::vector<Tensor>>& level_outputs);
    
    /**
     * @brief Predict at highest fidelity
     */
    Tensor predict(const Tensor& input);
    
    /**
     * @brief Predict at specific fidelity level
     */
    Tensor predictAtLevel(const Tensor& input, int level);
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    MFFNOConfig config_;
    
    // FNO at each fidelity level
    std::vector<std::unique_ptr<FNOModel>> level_models_;
    
    // Correction networks between levels
    std::vector<std::unique_ptr<AdditiveCorrection>> corrections_;
};

// =============================================================================
// Multi-Fidelity Gaussian Process (Co-Kriging)
// =============================================================================

/**
 * @brief Gaussian process for multi-fidelity modeling
 */
class CoKrigingModel {
public:
    struct CoKrigingConfig {
        int input_dim = 10;
        int output_dim = 1;
        int num_levels = 2;
        
        // Kernel settings
        std::string kernel = "rbf";  // rbf, matern32, matern52
        std::vector<double> length_scales;
        std::vector<double> variances;
        
        // Correlation between levels
        std::vector<double> rho;  // Autoregressive parameter
        
        // Optimization
        int max_optimization_iter = 100;
    };
    
    CoKrigingModel(const CoKrigingConfig& config);
    
    /**
     * @brief Predict with uncertainty
     */
    void predict(const Tensor& input, Tensor& mean, Tensor& variance);
    
    /**
     * @brief Add data at fidelity level
     */
    void addData(int level, const Tensor& inputs, const Tensor& outputs);
    
    /**
     * @brief Fit model (optimize hyperparameters)
     */
    void fit();
    
    /**
     * @brief Compute expected improvement for acquisition
     */
    double expectedImprovement(const Tensor& input, double best_value);
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    CoKrigingConfig config_;
    
    // Training data per level
    std::vector<std::vector<Tensor>> level_inputs_;
    std::vector<std::vector<Tensor>> level_outputs_;
    
    // Kernel matrices (per level)
    std::vector<Tensor> K_levels_;
    Tensor K_cross_;  // Cross-level covariance
    
    // Hyperparameters
    Tensor hyperparameters_;
    
    Tensor computeKernel(const Tensor& X1, const Tensor& X2, int level);
    void optimizeHyperparameters();
};

// =============================================================================
// Adaptive Fidelity Selection
// =============================================================================

/**
 * @brief Dynamically select appropriate fidelity level
 */
class AdaptiveFidelitySelector {
public:
    struct SelectorConfig {
        // Cost-accuracy trade-off
        double accuracy_weight = 0.5;  // vs cost
        
        // Confidence threshold for using lower fidelity
        double confidence_threshold = 0.9;
        
        // Fallback to high fidelity
        bool always_verify = false;
        double verification_fraction = 0.1;
    };
    
    AdaptiveFidelitySelector(const SelectorConfig& config);
    
    /**
     * @brief Select fidelity level for given input
     */
    int selectFidelity(const Tensor& input,
                      const std::vector<FidelityLevel>& levels,
                      std::function<double(const Tensor&, int)> confidence_estimator);
    
    /**
     * @brief Batch selection with budget constraint
     */
    std::vector<int> selectFidelitiesBatch(
        const std::vector<Tensor>& inputs,
        const std::vector<FidelityLevel>& levels,
        double total_budget);
    
    /**
     * @brief Update selection model based on observed errors
     */
    void updateFromObservations(const Tensor& input,
                               int selected_level,
                               double actual_error);
    
private:
    SelectorConfig config_;
    
    // Learned selection model
    std::unique_ptr<MLPModel> selection_model_;
    
    // Error history for adaptation
    std::vector<std::tuple<Tensor, int, double>> history_;
};

// =============================================================================
// Multi-Fidelity Active Learning
// =============================================================================

/**
 * @brief Active learning for multi-fidelity data acquisition
 */
class MultiFidelityActiveLearning {
public:
    struct ActiveConfig {
        // Budget allocation
        double total_budget = 1000.0;  // In computational units
        
        // Acquisition function
        enum class Acquisition {
            EXPECTED_IMPROVEMENT,
            KNOWLEDGE_GRADIENT,
            INFORMATION_GAIN,
            COST_AWARE_EI
        };
        Acquisition acquisition = Acquisition::COST_AWARE_EI;
        
        // Exploration vs exploitation
        double exploration_weight = 0.1;
    };
    
    MultiFidelityActiveLearning(const ActiveConfig& config);
    
    /**
     * @brief Select next query point and fidelity level
     */
    std::pair<Tensor, int> selectNext(
        const std::vector<FidelityLevel>& levels,
        const Tensor& candidate_points,
        std::function<void(const Tensor&, int, Tensor&, Tensor&)> predict_func);
    
    /**
     * @brief Compute acquisition function value
     */
    double computeAcquisition(const Tensor& point, int level,
                             const Tensor& prediction,
                             const Tensor& uncertainty,
                             double level_cost);
    
    /**
     * @brief Get optimal budget allocation across levels
     */
    std::vector<double> allocateBudget(
        const std::vector<FidelityLevel>& levels,
        double total_budget);
    
private:
    ActiveConfig config_;
    
    double remaining_budget_;
    
    double expectedImprovement(double mean, double std, double best);
    double costAwareEI(double ei, double cost);
};

// =============================================================================
// Unified Multi-Fidelity Controller
// =============================================================================

/**
 * @brief Main interface for multi-fidelity modeling
 */
class MultiFidelityController {
public:
    struct MFConfig {
        // Fidelity levels
        std::vector<FidelityLevel> levels;
        
        // Correction method
        enum class Method {
            ADDITIVE,
            MULTIPLICATIVE,
            RESIDUAL,
            TRANSFER,
            COKRIGING,
            NEURAL_OPERATOR
        };
        Method method = Method::ADDITIVE;
        
        // Adaptive selection
        bool use_adaptive_selection = true;
        
        // Active learning
        bool use_active_learning = false;
        double active_learning_budget = 100.0;
    };
    
    MultiFidelityController(const MFConfig& config);
    
    /**
     * @brief Register forward models for each level
     */
    void registerModel(int level, 
                      std::function<Tensor(const Tensor&)> model);
    
    /**
     * @brief Train multi-fidelity model
     */
    void train(const std::vector<std::vector<Tensor>>& level_inputs,
              const std::vector<std::vector<Tensor>>& level_outputs);
    
    /**
     * @brief Predict with automatic fidelity selection
     */
    Tensor predict(const Tensor& input);
    
    /**
     * @brief Predict at specified fidelity level
     */
    Tensor predictAtLevel(const Tensor& input, int level);
    
    /**
     * @brief Get prediction uncertainty
     */
    Tensor getUncertainty(const Tensor& input);
    
    /**
     * @brief Active learning: get next point to sample
     */
    std::pair<Tensor, int> getNextSamplePoint(
        const Tensor& candidate_points);
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    MFConfig config_;
    
    // Registered models
    std::map<int, std::function<Tensor(const Tensor&)>> models_;
    
    // Correction models
    std::unique_ptr<AdditiveCorrection> additive_;
    std::unique_ptr<MultiplicativeCorrection> multiplicative_;
    std::unique_ptr<ResidualFidelityModel> residual_;
    std::unique_ptr<TransferLearningMF> transfer_;
    std::unique_ptr<CoKrigingModel> cokriging_;
    std::unique_ptr<MultiFidelityFNO> mf_fno_;
    
    // Adaptive selection
    std::unique_ptr<AdaptiveFidelitySelector> selector_;
    
    // Active learning
    std::unique_ptr<MultiFidelityActiveLearning> active_learner_;
};

// =============================================================================
// Configuration Parsing
// =============================================================================

AdditiveCorrection::AdditiveConfig parseAdditiveConfig(
    const std::map<std::string, std::string>& config);
MultiFidelityController::MFConfig parseMFConfig(
    const std::map<std::string, std::string>& config);

} // namespace ML
} // namespace FSRM

#endif // FSRM_MULTI_FIDELITY_LEARNING_HPP
