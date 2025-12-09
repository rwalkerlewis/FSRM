/**
 * @file NeuralInversion.hpp
 * @brief Neural Network Acceleration for Inverse Problems
 * 
 * Accelerates inverse problems and history matching:
 * 
 * - Neural Surrogate Forward Models: Fast forward model evaluation
 * - Gradient Approximation: Learn sensitivities for optimization
 * - Learned EnKF: Neural-accelerated Ensemble Kalman Filter
 * - Physics-Guided Inversion: PINN-based parameter estimation
 * 
 * Applications:
 * - History matching for reservoir simulation
 * - Parameter estimation from observations
 * - Seismic inversion
 * - Material property identification
 */

#ifndef FSRM_NEURAL_INVERSION_HPP
#define FSRM_NEURAL_INVERSION_HPP

#include "FourierNeuralOperator.hpp"
#include "NeuralSurrogates.hpp"
#include "BayesianNeuralOperator.hpp"
#include "NeuralReducedOrderModel.hpp"
#include <petsc.h>
#include <vector>
#include <memory>
#include <functional>

namespace FSRM {
namespace ML {

// =============================================================================
// Neural Forward Model Surrogate
// =============================================================================

/**
 * @brief Neural network surrogate for forward model
 * 
 * Maps: parameters -> observations
 * Much faster than running full simulation
 */
class NeuralForwardModel {
public:
    struct ForwardModelConfig {
        // Input/output dimensions
        int num_parameters = 100;  // Parameter dimension
        int num_observations = 50; // Observation dimension
        
        // Network architecture
        std::vector<int> hidden_layers = {256, 512, 512, 256};
        std::string activation = "gelu";
        
        // Time-dependent forward model
        bool time_dependent = false;
        int num_time_steps = 10;
        
        // Training
        double learning_rate = 1e-3;
        int epochs = 500;
        int batch_size = 64;
        
        // Data augmentation
        double noise_level = 0.01;
        bool use_latin_hypercube = true;
        int num_training_samples = 1000;
    };
    
    NeuralForwardModel(const ForwardModelConfig& config);
    
    /**
     * @brief Evaluate forward model
     */
    Tensor forward(const Tensor& parameters);
    
    /**
     * @brief Evaluate with time
     */
    Tensor forward(const Tensor& parameters, double time);
    
    /**
     * @brief Batch evaluation (more efficient)
     */
    Tensor forwardBatch(const Tensor& parameter_batch);
    
    /**
     * @brief Compute Jacobian (sensitivities)
     */
    Tensor jacobian(const Tensor& parameters);
    
    /**
     * @brief Train from simulation runs
     */
    void train(const std::vector<Tensor>& parameters,
              const std::vector<Tensor>& observations);
    
    /**
     * @brief Generate training data from full simulator
     */
    void generateTrainingData(
        std::function<Tensor(const Tensor&)> full_simulator,
        const Tensor& parameter_lower_bounds,
        const Tensor& parameter_upper_bounds,
        int num_samples);
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
    double getValidationError() const { return validation_error_; }
    
private:
    ForwardModelConfig config_;
    std::unique_ptr<MLPModel> model_;
    
    // Normalization
    Tensor param_mean_, param_std_;
    Tensor obs_mean_, obs_std_;
    
    double validation_error_ = 0.0;
    
    std::vector<Tensor> latinHypercubeSampling(
        const Tensor& lower, const Tensor& upper, int n);
};

// =============================================================================
// Neural Gradient Approximator
// =============================================================================

/**
 * @brief Learn sensitivities for gradient-based optimization
 */
class NeuralGradientApproximator {
public:
    struct GradientConfig {
        int parameter_dim = 100;
        int observation_dim = 50;
        
        // Architecture
        std::vector<int> hidden_layers = {256, 256};
        
        // Training method
        enum class TrainingMethod {
            FINITE_DIFFERENCE,  // Learn from FD approximations
            ADJOINT,            // Learn from adjoint gradients
            DIRECT              // Learn Jacobian directly
        };
        TrainingMethod method = TrainingMethod::FINITE_DIFFERENCE;
        
        double finite_diff_eps = 1e-6;
    };
    
    NeuralGradientApproximator(const GradientConfig& config);
    
    /**
     * @brief Compute gradient of objective w.r.t. parameters
     */
    Tensor computeGradient(const Tensor& parameters,
                          const Tensor& observations,
                          const Tensor& observed_data);
    
    /**
     * @brief Compute full Jacobian matrix
     */
    Tensor computeJacobian(const Tensor& parameters);
    
    /**
     * @brief Train from parameter-gradient pairs
     */
    void train(const std::vector<Tensor>& parameters,
              const std::vector<Tensor>& gradients);
    
    /**
     * @brief Train from forward model using finite differences
     */
    void trainFromForwardModel(
        std::function<Tensor(const Tensor&)> forward_model,
        const std::vector<Tensor>& sample_parameters);
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    GradientConfig config_;
    
    // Option 1: Direct gradient prediction
    std::unique_ptr<MLPModel> gradient_net_;
    
    // Option 2: Jacobian prediction
    std::unique_ptr<MLPModel> jacobian_net_;
};

// =============================================================================
// Neural EnKF
// =============================================================================

/**
 * @brief Neural-accelerated Ensemble Kalman Filter
 */
class NeuralEnKF {
public:
    struct EnKFConfig {
        int ensemble_size = 100;
        int parameter_dim = 100;
        int observation_dim = 50;
        
        // Neural acceleration
        bool use_neural_forward = true;
        bool use_neural_localization = true;
        bool use_neural_inflation = true;
        
        // Localization
        double localization_radius = 0.5;
        
        // Adaptive inflation
        double inflation_factor = 1.05;
        bool adaptive_inflation = true;
        
        // Training
        int epochs = 100;
    };
    
    NeuralEnKF(const EnKFConfig& config);
    
    /**
     * @brief Initialize ensemble
     */
    void initialize(const Tensor& prior_mean, const Tensor& prior_covariance);
    
    /**
     * @brief Forecast step (propagate ensemble)
     */
    void forecast(std::function<Tensor(const Tensor&)> forward_model);
    
    /**
     * @brief Analysis step (assimilate observations)
     */
    void analysis(const Tensor& observations, 
                 const Tensor& observation_covariance);
    
    /**
     * @brief Full update cycle
     */
    void update(const Tensor& observations,
               const Tensor& observation_covariance,
               std::function<Tensor(const Tensor&)> forward_model);
    
    /**
     * @brief Get current ensemble mean
     */
    Tensor getMean() const;
    
    /**
     * @brief Get current ensemble covariance
     */
    Tensor getCovariance() const;
    
    /**
     * @brief Get individual ensemble members
     */
    const std::vector<Tensor>& getEnsemble() const { return ensemble_; }
    
    /**
     * @brief Train neural components
     */
    void trainNeuralComponents(
        const std::vector<Tensor>& historical_parameters,
        const std::vector<Tensor>& historical_observations);
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    EnKFConfig config_;
    
    std::vector<Tensor> ensemble_;
    Tensor ensemble_predictions_;
    
    // Neural components
    std::unique_ptr<NeuralForwardModel> neural_forward_;
    std::unique_ptr<MLPModel> localization_net_;
    std::unique_ptr<MLPModel> inflation_net_;
    
    Tensor computeKalmanGain(const Tensor& innovations,
                            const Tensor& obs_covariance);
    Tensor neuralLocalization(const Tensor& kalman_gain);
    double computeAdaptiveInflation(const Tensor& innovations,
                                    const Tensor& obs_covariance);
};

// =============================================================================
// Physics-Guided Inversion
// =============================================================================

/**
 * @brief PINN-based parameter estimation
 */
class PhysicsGuidedInversion {
public:
    struct InversionConfig {
        // Dimensions
        int parameter_dim = 100;
        int state_dim = 1000;
        int observation_dim = 50;
        
        // Network architecture
        std::vector<int> encoder_layers = {256, 256};
        std::vector<int> physics_layers = {256, 512, 256};
        
        // Physics constraints
        std::string pde_type = "darcy";  // darcy, wave, heat, elastic
        double physics_weight = 0.1;
        
        // Regularization
        double l2_regularization = 1e-4;
        double total_variation = 0.0;
        
        // Optimization
        double learning_rate = 1e-3;
        int max_iterations = 1000;
        double convergence_tol = 1e-6;
    };
    
    PhysicsGuidedInversion(const InversionConfig& config);
    
    /**
     * @brief Invert for parameters from observations
     */
    Tensor invert(const Tensor& observations,
                 const Tensor& observation_locations,
                 const Tensor& initial_guess);
    
    /**
     * @brief Invert with uncertainty
     */
    void invertWithUncertainty(const Tensor& observations,
                              const Tensor& observation_locations,
                              Tensor& parameters,
                              Tensor& uncertainty);
    
    /**
     * @brief Evaluate physics residual
     */
    Tensor computePhysicsResidual(const Tensor& parameters,
                                  const Tensor& state);
    
    /**
     * @brief Compute data misfit
     */
    double computeDataMisfit(const Tensor& parameters,
                            const Tensor& observations);
    
    /**
     * @brief Train physics-informed network
     */
    void train(const std::vector<Tensor>& parameter_samples,
              const std::vector<Tensor>& corresponding_states);
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    InversionConfig config_;
    
    // Encoder: observations -> parameter estimate
    std::unique_ptr<MLPModel> encoder_;
    
    // Physics network: parameters -> state
    std::unique_ptr<MLPModel> physics_net_;
    
    // For uncertainty estimation
    std::unique_ptr<BayesianNeuralNetwork> uncertainty_net_;
    
    // PDE residual computation
    Tensor computeDarcyResidual(const Tensor& k, const Tensor& p);
    Tensor computeWaveResidual(const Tensor& c, const Tensor& u);
};

// =============================================================================
// Iterative Neural Inversion
// =============================================================================

/**
 * @brief Learned iterative inversion scheme
 */
class IterativeNeuralInversion {
public:
    struct IterativeConfig {
        int max_iterations = 50;
        double tolerance = 1e-6;
        
        // Network for update step
        std::vector<int> update_layers = {256, 256};
        
        // Learned step size
        bool learn_step_size = true;
        double initial_step_size = 1.0;
        
        // Momentum
        bool use_momentum = true;
        double momentum = 0.9;
    };
    
    IterativeNeuralInversion(const IterativeConfig& config);
    
    /**
     * @brief Learned iterative inversion
     * 
     * Uses neural network to predict optimal update direction
     * at each iteration, similar to learned optimization.
     */
    Tensor solve(const Tensor& observations,
                std::function<Tensor(const Tensor&)> forward_model,
                const Tensor& initial_parameters);
    
    /**
     * @brief Single iteration step
     */
    Tensor step(const Tensor& current_params,
               const Tensor& current_residual,
               const Tensor& gradient);
    
    /**
     * @brief Train on inversion trajectories
     */
    void train(const std::vector<std::vector<Tensor>>& param_trajectories,
              const std::vector<std::vector<Tensor>>& residual_trajectories);
    
private:
    IterativeConfig config_;
    
    // Update direction predictor
    std::unique_ptr<MLPModel> update_net_;
    
    // Step size predictor
    std::unique_ptr<MLPModel> step_size_net_;
    
    // LSTM for trajectory-aware updates
    std::unique_ptr<LSTMCell> lstm_;
};

// =============================================================================
// History Matching Controller
// =============================================================================

/**
 * @brief Unified interface for history matching
 */
class HistoryMatchingController {
public:
    enum class Method {
        NEURAL_ENKF,
        PHYSICS_GUIDED,
        ITERATIVE_NEURAL,
        GRADIENT_DESCENT,
        HYBRID
    };
    
    struct HistoryMatchConfig {
        Method method = Method::HYBRID;
        
        // Problem definition
        int num_parameters = 100;
        int num_observations = 50;
        int num_time_steps = 10;
        
        // Prior information
        bool has_prior = false;
        
        // Convergence
        double tolerance = 1e-4;
        int max_iterations = 100;
        
        // Component configs
        NeuralEnKF::EnKFConfig enkf_config;
        PhysicsGuidedInversion::InversionConfig physics_config;
        IterativeNeuralInversion::IterativeConfig iterative_config;
    };
    
    HistoryMatchingController(const HistoryMatchConfig& config);
    
    /**
     * @brief Set prior distribution
     */
    void setPrior(const Tensor& prior_mean, const Tensor& prior_covariance);
    
    /**
     * @brief Add observation at time
     */
    void addObservation(double time, const Tensor& observation,
                       const Tensor& observation_error);
    
    /**
     * @brief Set forward model
     */
    void setForwardModel(std::function<Tensor(const Tensor&, double)> model);
    
    /**
     * @brief Run history matching
     */
    Tensor solve(const Tensor& initial_guess);
    
    /**
     * @brief Get uncertainty estimate
     */
    Tensor getUncertainty() const { return uncertainty_; }
    
    /**
     * @brief Get matched trajectory
     */
    std::vector<Tensor> getMatchedTrajectory() const { return matched_trajectory_; }
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    HistoryMatchConfig config_;
    
    // Components
    std::unique_ptr<NeuralEnKF> enkf_;
    std::unique_ptr<PhysicsGuidedInversion> physics_inversion_;
    std::unique_ptr<IterativeNeuralInversion> iterative_;
    std::unique_ptr<NeuralForwardModel> neural_forward_;
    std::unique_ptr<NeuralGradientApproximator> neural_gradient_;
    
    // Problem data
    std::vector<double> observation_times_;
    std::vector<Tensor> observations_;
    std::vector<Tensor> observation_errors_;
    std::function<Tensor(const Tensor&, double)> forward_model_;
    
    Tensor prior_mean_, prior_covariance_;
    Tensor uncertainty_;
    std::vector<Tensor> matched_trajectory_;
    
    // Solution methods
    Tensor solveEnKF(const Tensor& initial);
    Tensor solvePhysicsGuided(const Tensor& initial);
    Tensor solveIterative(const Tensor& initial);
    Tensor solveHybrid(const Tensor& initial);
};

// =============================================================================
// Configuration Parsing
// =============================================================================

NeuralForwardModel::ForwardModelConfig parseForwardModelConfig(
    const std::map<std::string, std::string>& config);
NeuralEnKF::EnKFConfig parseEnKFConfig(
    const std::map<std::string, std::string>& config);

} // namespace ML
} // namespace FSRM

#endif // FSRM_NEURAL_INVERSION_HPP
