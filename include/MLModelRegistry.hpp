/**
 * @file MLModelRegistry.hpp
 * @brief Unified Machine Learning Model Registry and Selection System
 * 
 * Provides a single interface for:
 * - Registering and selecting ML models
 * - Model lifecycle management (train, save, load)
 * - Automatic model selection based on problem type
 * - Ensemble and hybrid model combinations
 * - Runtime model switching
 * 
 * Enables seamless integration of ML with classical solvers
 * where ML can be used "in place of or in addition to" existing methods.
 */

#ifndef FSRM_ML_MODEL_REGISTRY_HPP
#define FSRM_ML_MODEL_REGISTRY_HPP

// Include all ML components
#include "FourierNeuralOperator.hpp"
#include "NeuralSurrogates.hpp"
#include "GraphNeuralOperator.hpp"
#include "NeuralReducedOrderModel.hpp"
#include "NeuralLinearAlgebra.hpp"
#include "BayesianNeuralOperator.hpp"
#include "NeuralInversion.hpp"
#include "MultiFidelityLearning.hpp"
#include "NeuralTimestepping.hpp"
#include "NeuralAMR.hpp"

#include <petsc.h>
#include <vector>
#include <memory>
#include <map>
#include <functional>
#include <any>

namespace FSRM {
namespace ML {

// =============================================================================
// Model Categories
// =============================================================================

/**
 * @brief Categories of ML models
 */
enum class ModelCategory {
    // Forward Operators
    FNO,                    // Fourier Neural Operator
    GNO,                    // Graph Neural Operator
    DEEPONET,               // Deep Operator Network
    ROM,                    // Reduced Order Model
    
    // Surrogates
    PVT_SURROGATE,          // PVT property surrogate
    FLASH_SURROGATE,        // Flash calculation surrogate
    RELPERM_SURROGATE,      // Relative permeability surrogate
    EOS_SURROGATE,          // Equation of state surrogate
    
    // Linear Algebra
    NEURAL_PRECONDITIONER,  // Neural preconditioner
    NEURAL_COARSE_GRID,     // Coarse grid correction
    NEURAL_AMG,             // Neural algebraic multigrid
    
    // Time Integration
    NEURAL_TIMESTEP,        // Time step selection
    NEURAL_IMEX,            // IMEX transition
    
    // Mesh Adaptation
    NEURAL_AMR,             // Adaptive mesh refinement
    
    // Uncertainty
    BAYESIAN_NN,            // Bayesian neural network
    DEEP_ENSEMBLE,          // Deep ensemble
    
    // Inversion
    NEURAL_FORWARD_MODEL,   // Forward model surrogate
    NEURAL_ENKF,            // Neural EnKF
    
    // Multi-Fidelity
    MF_CORRECTION,          // Multi-fidelity correction
    MF_COKRIGING            // Co-Kriging
};

/**
 * @brief Model status
 */
enum class ModelStatus {
    UNINITIALIZED,
    INITIALIZED,
    TRAINED,
    LOADED,
    ACTIVE,
    DISABLED
};

// =============================================================================
// Model Registration Entry
// =============================================================================

/**
 * @brief Entry for registered model
 */
struct ModelEntry {
    std::string name;
    ModelCategory category;
    ModelStatus status = ModelStatus::UNINITIALIZED;
    
    // Model pointer (type-erased)
    std::any model;
    
    // Configuration
    std::map<std::string, std::string> config;
    
    // Metadata
    std::string description;
    size_t num_parameters = 0;
    double training_error = 0.0;
    double validation_error = 0.0;
    double speedup_factor = 1.0;
    
    // Usage statistics
    int num_calls = 0;
    double total_time = 0.0;
    double avg_time = 0.0;
    
    // Model path for persistence
    std::string model_path;
};

// =============================================================================
// ML Model Registry
// =============================================================================

/**
 * @brief Central registry for all ML models
 */
class MLModelRegistry {
public:
    static MLModelRegistry& getInstance();
    
    // Prevent copying
    MLModelRegistry(const MLModelRegistry&) = delete;
    MLModelRegistry& operator=(const MLModelRegistry&) = delete;
    
    // ==========================================================================
    // Model Registration
    // ==========================================================================
    
    /**
     * @brief Register FNO model
     */
    void registerFNO(const std::string& name, const FNOConfig& config);
    
    /**
     * @brief Register GNO model
     */
    void registerGNO(const std::string& name, const GNOConfig& config);
    
    /**
     * @brief Register neural surrogate
     */
    void registerPVTSurrogate(const std::string& name, 
                             const NeuralPVTConfig& config);
    void registerFlashSurrogate(const std::string& name,
                               const NeuralFlashConfig& config);
    void registerRelPermSurrogate(const std::string& name,
                                 const NeuralRelPermConfig& config);
    
    /**
     * @brief Register ROM
     */
    void registerROM(const std::string& name, const PODNNConfig& config);
    void registerAutoencoderROM(const std::string& name,
                               const AutoencoderROMConfig& config);
    
    /**
     * @brief Register linear algebra models
     */
    void registerNeuralPreconditioner(const std::string& name,
                                     const NeuralPreconditioner::PrecondConfig& config);
    void registerNeuralCoarseGrid(const std::string& name,
                                 const NeuralCoarseGridConfig& config);
    
    /**
     * @brief Register time stepping models
     */
    void registerNeuralTimeStep(const std::string& name,
                               const NeuralTimeStepPredictor::TimeStepConfig& config);
    void registerNeuralIMEX(const std::string& name,
                           const NeuralIMEXPredictor::IMEXConfig& config);
    
    /**
     * @brief Register AMR model
     */
    void registerNeuralAMR(const std::string& name,
                          const NeuralAMRController::AMRConfig& config);
    
    /**
     * @brief Register uncertainty model
     */
    void registerBayesianNN(const std::string& name, const BNNConfig& config);
    void registerDeepEnsemble(const std::string& name, 
                             const DeepEnsembleConfig& config);
    
    /**
     * @brief Register inversion model
     */
    void registerNeuralEnKF(const std::string& name,
                           const NeuralEnKF::EnKFConfig& config);
    
    /**
     * @brief Register multi-fidelity model
     */
    void registerMultiFidelity(const std::string& name,
                              const MultiFidelityController::MFConfig& config);
    
    /**
     * @brief Generic registration
     */
    template<typename ModelType>
    void registerModel(const std::string& name, ModelCategory category,
                      std::shared_ptr<ModelType> model);
    
    // ==========================================================================
    // Model Access
    // ==========================================================================
    
    /**
     * @brief Get model by name
     */
    template<typename ModelType>
    ModelType* getModel(const std::string& name);
    
    /**
     * @brief Get models by category
     */
    std::vector<std::string> getModelsByCategory(ModelCategory category) const;
    
    /**
     * @brief Check if model exists
     */
    bool hasModel(const std::string& name) const;
    
    /**
     * @brief Get model status
     */
    ModelStatus getModelStatus(const std::string& name) const;
    
    /**
     * @brief Get all registered model names
     */
    std::vector<std::string> getAllModelNames() const;
    
    // ==========================================================================
    // Model Lifecycle
    // ==========================================================================
    
    /**
     * @brief Train model
     */
    void train(const std::string& name,
              const std::vector<Tensor>& inputs,
              const std::vector<Tensor>& outputs);
    
    /**
     * @brief Save model to file
     */
    void save(const std::string& name, const std::string& path);
    
    /**
     * @brief Load model from file
     */
    void load(const std::string& name, const std::string& path);
    
    /**
     * @brief Save all models
     */
    void saveAll(const std::string& directory);
    
    /**
     * @brief Load all models from directory
     */
    void loadAll(const std::string& directory);
    
    /**
     * @brief Enable/disable model
     */
    void setModelStatus(const std::string& name, ModelStatus status);
    
    // ==========================================================================
    // Model Selection
    // ==========================================================================
    
    /**
     * @brief Select best model for task
     */
    std::string selectBestModel(ModelCategory category,
                               const std::map<std::string, double>& criteria = {});
    
    /**
     * @brief Get recommended model for problem type
     */
    std::string recommendModel(const std::string& problem_type);
    
    /**
     * @brief Benchmark models and select fastest
     */
    std::string benchmarkAndSelect(ModelCategory category,
                                  const Tensor& sample_input);
    
    // ==========================================================================
    // Hybrid/Ensemble Management
    // ==========================================================================
    
    /**
     * @brief Create ensemble from multiple models
     */
    void createEnsemble(const std::string& ensemble_name,
                       const std::vector<std::string>& model_names,
                       const std::vector<double>& weights = {});
    
    /**
     * @brief Create hybrid model (ML + classical)
     */
    void createHybrid(const std::string& hybrid_name,
                     const std::string& ml_model,
                     const std::string& classical_method);
    
    // ==========================================================================
    // Statistics and Monitoring
    // ==========================================================================
    
    /**
     * @brief Get model statistics
     */
    const ModelEntry& getModelInfo(const std::string& name) const;
    
    /**
     * @brief Record model call
     */
    void recordCall(const std::string& name, double elapsed_time);
    
    /**
     * @brief Print summary of all models
     */
    void printSummary() const;
    
    /**
     * @brief Reset statistics
     */
    void resetStatistics();
    
    // ==========================================================================
    // Configuration
    // ==========================================================================
    
    /**
     * @brief Load configuration from file
     */
    void loadConfig(const std::string& config_file);
    
    /**
     * @brief Set global ML options
     */
    void setGlobalOption(const std::string& key, const std::string& value);
    
    /**
     * @brief Get global option
     */
    std::string getGlobalOption(const std::string& key) const;
    
private:
    MLModelRegistry() = default;
    
    std::map<std::string, ModelEntry> models_;
    std::map<std::string, std::string> global_options_;
    
    // Default model paths
    std::string model_directory_ = "ml_models";
    
    // Helper for category to string
    static std::string categoryToString(ModelCategory cat);
    static ModelCategory stringToCategory(const std::string& str);
};

// =============================================================================
// ML Solver Wrapper
// =============================================================================

/**
 * @brief Wrapper to use ML models as drop-in replacements
 */
class MLSolverWrapper {
public:
    /**
     * @brief Create wrapper for category of models
     */
    MLSolverWrapper(ModelCategory category);
    
    /**
     * @brief Set active model
     */
    void setActiveModel(const std::string& model_name);
    
    /**
     * @brief Enable/disable ML (falls back to classical)
     */
    void enableML(bool enable = true);
    
    /**
     * @brief Set classical fallback
     */
    void setFallback(std::function<Tensor(const Tensor&)> fallback);
    
    /**
     * @brief Solve using ML or fallback
     */
    Tensor solve(const Tensor& input);
    
    /**
     * @brief Solve with confidence threshold
     */
    Tensor solveWithConfidence(const Tensor& input, double threshold,
                              bool& used_ml);
    
    /**
     * @brief Get whether ML is currently being used
     */
    bool isMLActive() const { return ml_enabled_ && active_model_name_.size() > 0; }
    
private:
    ModelCategory category_;
    std::string active_model_name_;
    bool ml_enabled_ = true;
    
    std::function<Tensor(const Tensor&)> fallback_;
    
    double confidence_threshold_ = 0.8;
};

// =============================================================================
// Helper Functions
// =============================================================================

/**
 * @brief Parse ML model configuration from config file
 */
void parseMLConfig(const std::map<std::string, std::string>& config);

/**
 * @brief Initialize all ML models from configuration
 */
void initializeML(const std::string& config_file);

/**
 * @brief Shutdown ML system (cleanup)
 */
void shutdownML();

/**
 * @brief Get summary string of ML system status
 */
std::string getMLStatus();

// =============================================================================
// Macros for Easy Integration
// =============================================================================

// Register model macro
#define REGISTER_ML_MODEL(name, category, config_type, config) \
    FSRM::ML::MLModelRegistry::getInstance().registerModel<config_type>( \
        name, category, std::make_shared<config_type>(config))

// Get model macro
#define GET_ML_MODEL(name, type) \
    FSRM::ML::MLModelRegistry::getInstance().getModel<type>(name)

// Check if ML available
#define ML_AVAILABLE(name) \
    (FSRM::ML::MLModelRegistry::getInstance().hasModel(name) && \
     FSRM::ML::MLModelRegistry::getInstance().getModelStatus(name) == \
         FSRM::ML::ModelStatus::TRAINED)

// =============================================================================
// Configuration File Format
// =============================================================================

/**
 * Example configuration file (ml_models.config):
 * 
 * [global]
 * model_directory = ml_models
 * use_gpu = true
 * default_precision = float32
 * 
 * [fno_pressure]
 * type = FNO
 * modes = 12
 * width = 64
 * layers = 4
 * learning_rate = 1e-3
 * epochs = 100
 * model_path = ml_models/fno_pressure.model
 * 
 * [pvt_surrogate]
 * type = PVT_SURROGATE
 * hidden_layers = 128,256,256,128
 * activation = gelu
 * p_min = 1e5
 * p_max = 1e8
 * 
 * [neural_amr]
 * type = NEURAL_AMR
 * use_predictive = true
 * use_physics_aware = true
 * max_level = 5
 */

} // namespace ML
} // namespace FSRM

// =============================================================================
// Template Implementations
// =============================================================================

namespace FSRM {
namespace ML {

template<typename ModelType>
void MLModelRegistry::registerModel(const std::string& name, 
                                   ModelCategory category,
                                   std::shared_ptr<ModelType> model) {
    ModelEntry entry;
    entry.name = name;
    entry.category = category;
    entry.status = ModelStatus::INITIALIZED;
    entry.model = model;
    entry.description = categoryToString(category) + " model: " + name;
    
    models_[name] = std::move(entry);
}

template<typename ModelType>
ModelType* MLModelRegistry::getModel(const std::string& name) {
    auto it = models_.find(name);
    if (it == models_.end()) {
        return nullptr;
    }
    
    try {
        auto model_ptr = std::any_cast<std::shared_ptr<ModelType>>(it->second.model);
        return model_ptr.get();
    } catch (const std::bad_any_cast&) {
        return nullptr;
    }
}

} // namespace ML
} // namespace FSRM

#endif // FSRM_ML_MODEL_REGISTRY_HPP
