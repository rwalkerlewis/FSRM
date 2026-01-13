/**
 * @file BayesianNeuralOperator.hpp
 * @brief Bayesian Neural Networks for Uncertainty Quantification
 * 
 * Provides uncertainty-aware neural operators:
 * 
 * - Bayesian Neural Networks: Weight distributions for uncertainty
 * - Deep Ensembles: Multiple models for prediction variance
 * - Monte Carlo Dropout: Dropout-based uncertainty estimation
 * - Probabilistic FNO: Uncertainty-aware Fourier Neural Operator
 * 
 * Critical for:
 * - Quantifying prediction confidence
 * - Detecting out-of-distribution inputs
 * - Risk-aware decision making
 * - Active learning for data acquisition
 */

#ifndef FSRM_BAYESIAN_NEURAL_OPERATOR_HPP
#define FSRM_BAYESIAN_NEURAL_OPERATOR_HPP

#include "FourierNeuralOperator.hpp"
#include "NeuralSurrogates.hpp"
#include <petsc.h>
#include <vector>
#include <memory>
#include <random>

namespace FSRM {
namespace ML {

// =============================================================================
// Variational Layers for BNN
// =============================================================================

/**
 * @brief Variational linear layer with weight distributions
 */
class VariationalLinear {
public:
    VariationalLinear(int in_features, int out_features, 
                     double prior_std = 1.0);
    
    /**
     * @brief Forward pass with weight sampling
     */
    Tensor forward(const Tensor& x, bool sample = true);
    
    /**
     * @brief KL divergence from prior
     */
    double klDivergence() const;
    
    // Parameters (mean and log variance)
    Tensor weight_mu;
    Tensor weight_logvar;
    Tensor bias_mu;
    Tensor bias_logvar;
    
    // Gradients
    Tensor grad_weight_mu;
    Tensor grad_weight_logvar;
    Tensor grad_bias_mu;
    Tensor grad_bias_logvar;
    
    std::vector<Tensor*> parameters();
    std::vector<Tensor*> gradients();
    void zeroGrad();
    
private:
    int in_features_, out_features_;
    double prior_std_;
    Tensor sampled_weight_, sampled_bias_;
    Tensor cached_input_;
    
    std::mt19937 rng_;
};

// =============================================================================
// Bayesian Neural Network
// =============================================================================

/**
 * @brief Configuration for Bayesian Neural Network
 */
struct BNNConfig {
    std::vector<int> layer_sizes;
    std::string activation = "relu";
    
    // Prior
    double prior_std = 1.0;
    
    // Training
    double learning_rate = 1e-3;
    int epochs = 200;
    int batch_size = 32;
    double kl_weight = 1e-3;  // Weight for KL divergence term
    
    // Inference
    int num_samples = 100;  // MC samples for prediction
};

/**
 * @brief Full Bayesian Neural Network
 */
class BayesianNeuralNetwork {
public:
    BayesianNeuralNetwork(const BNNConfig& config);
    
    /**
     * @brief Forward with uncertainty
     * 
     * @param x Input
     * @param mean Output: Mean prediction
     * @param variance Output: Predictive variance
     */
    void predict(const Tensor& x, Tensor& mean, Tensor& variance);
    
    /**
     * @brief Single sample forward
     */
    Tensor forward(const Tensor& x, bool sample = true);
    
    /**
     * @brief Compute ELBO loss
     */
    double computeELBO(const Tensor& prediction, const Tensor& target);
    
    /**
     * @brief Train with variational inference
     */
    void train(const std::vector<Tensor>& inputs,
              const std::vector<Tensor>& targets);
    
    /**
     * @brief Get epistemic uncertainty (model uncertainty)
     */
    Tensor epistemicUncertainty(const Tensor& x);
    
    /**
     * @brief Get aleatoric uncertainty (data uncertainty)
     */
    Tensor aleatoricUncertainty(const Tensor& x);
    
    /**
     * @brief Check if input is out-of-distribution
     */
    bool isOutOfDistribution(const Tensor& x, double threshold = 2.0);
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    BNNConfig config_;
    std::vector<std::unique_ptr<VariationalLinear>> layers_;
    
    double total_kl_divergence() const;
};

// =============================================================================
// Deep Ensemble
// =============================================================================

/**
 * @brief Configuration for Deep Ensemble
 */
struct DeepEnsembleConfig {
    int num_models = 5;
    
    // Individual model config
    std::vector<int> layer_sizes;
    std::string activation = "relu";
    
    // Training
    double learning_rate = 1e-3;
    int epochs = 100;
    
    // Diversity
    bool use_different_seeds = true;
    bool use_bootstrap = false;  // Bootstrap training data
    double data_fraction = 1.0;  // Fraction of data per model
};

/**
 * @brief Ensemble of neural networks for uncertainty
 */
class DeepEnsemble {
public:
    DeepEnsemble(const DeepEnsembleConfig& config);
    
    /**
     * @brief Predict with mean and variance
     */
    void predict(const Tensor& x, Tensor& mean, Tensor& variance);
    
    /**
     * @brief Get individual model predictions
     */
    std::vector<Tensor> getModelPredictions(const Tensor& x);
    
    /**
     * @brief Train all ensemble members
     */
    void train(const std::vector<Tensor>& inputs,
              const std::vector<Tensor>& targets);
    
    /**
     * @brief Compute disagreement (uncertainty metric)
     */
    double computeDisagreement(const Tensor& x);
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    DeepEnsembleConfig config_;
    std::vector<std::unique_ptr<MLPModel>> models_;
    
    std::vector<int> getBootstrapIndices(int n);
};

// =============================================================================
// Monte Carlo Dropout
// =============================================================================

/**
 * @brief MC Dropout for uncertainty estimation
 */
class MCDropout {
public:
    struct MCDropoutConfig {
        std::vector<int> layer_sizes;
        std::string activation = "relu";
        double dropout_rate = 0.1;
        
        // Training
        double learning_rate = 1e-3;
        int epochs = 100;
        
        // Inference
        int num_mc_samples = 50;
    };
    
    MCDropout(const MCDropoutConfig& config);
    
    /**
     * @brief Predict with uncertainty
     */
    void predict(const Tensor& x, Tensor& mean, Tensor& variance);
    
    /**
     * @brief Single forward pass (dropout always on)
     */
    Tensor forward(const Tensor& x);
    
    /**
     * @brief Train (dropout during training as usual)
     */
    void train(const std::vector<Tensor>& inputs,
              const std::vector<Tensor>& targets);
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    MCDropoutConfig config_;
    std::unique_ptr<MLPModel> model_;
};

// =============================================================================
// Probabilistic Fourier Neural Operator
// =============================================================================

/**
 * @brief FNO with uncertainty quantification
 */
class ProbabilisticFNO {
public:
    struct ProbFNOConfig : FNOConfig {
        // Uncertainty method
        enum class UQMethod {
            BAYESIAN_WEIGHTS,    // Variational weights
            ENSEMBLE,            // Multiple FNOs
            MC_DROPOUT,          // Dropout at inference
            DEEP_EVIDENTIAL      // Evidential deep learning
        };
        UQMethod uq_method = UQMethod::ENSEMBLE;
        
        // For ensemble
        int num_ensemble_members = 5;
        
        // For Bayesian
        double prior_std = 1.0;
        double kl_weight = 1e-4;
        
        // For MC Dropout
        double dropout_rate = 0.1;
        int num_mc_samples = 50;
    };
    
    ProbabilisticFNO(const ProbFNOConfig& config);
    
    /**
     * @brief Predict field with uncertainty
     */
    void predict(const Tensor& input_field, 
                Tensor& mean_prediction,
                Tensor& variance_prediction);
    
    /**
     * @brief Get per-point uncertainty
     */
    Tensor getPointwiseUncertainty(const Tensor& input_field);
    
    /**
     * @brief Get overall prediction confidence
     */
    double getPredictionConfidence(const Tensor& input_field);
    
    /**
     * @brief Train with uncertainty-aware loss
     */
    void train(const std::vector<Tensor>& inputs,
              const std::vector<Tensor>& targets);
    
    /**
     * @brief Active learning: find most informative sample
     */
    int selectNextSample(const std::vector<Tensor>& candidates);
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    ProbFNOConfig config_;
    
    // Ensemble members or single model
    std::vector<std::unique_ptr<FNOModel>> ensemble_;
    std::unique_ptr<BayesianNeuralNetwork> bayesian_fno_;
    
    // For evidential learning
    std::unique_ptr<MLPModel> evidence_head_;
    
    Tensor computeEvidentialUncertainty(const Tensor& input);
};

// =============================================================================
// Uncertainty Propagation
// =============================================================================

/**
 * @brief Propagate uncertainty through physical simulation
 */
class UncertaintyPropagator {
public:
    struct PropagatorConfig {
        enum class Method {
            MONTE_CARLO,            // Full MC sampling
            UNSCENTED_TRANSFORM,    // Sigma point method
            POLYNOMIAL_CHAOS,       // PCE
            NEURAL_SURROGATE        // Learn uncertainty mapping
        };
        Method method = Method::NEURAL_SURROGATE;
        
        // For MC
        int num_mc_samples = 1000;
        
        // For Unscented Transform
        double alpha = 1e-3;
        double beta = 2.0;
        double kappa = 0.0;
        
        // For Neural surrogate
        std::vector<int> surrogate_layers = {128, 256, 128};
    };
    
    UncertaintyPropagator(const PropagatorConfig& config);
    
    /**
     * @brief Propagate input uncertainty through model
     * 
     * @param input_mean Mean of input distribution
     * @param input_covariance Covariance of input
     * @param output_mean Output: Mean of output distribution
     * @param output_covariance Output: Covariance of output
     */
    void propagate(const Tensor& input_mean,
                  const Tensor& input_covariance,
                  std::function<Tensor(const Tensor&)> forward_model,
                  Tensor& output_mean,
                  Tensor& output_covariance);
    
    /**
     * @brief Train neural surrogate for uncertainty propagation
     */
    void trainSurrogate(std::function<Tensor(const Tensor&)> forward_model,
                       const std::vector<Tensor>& input_samples);
    
private:
    PropagatorConfig config_;
    
    // Neural surrogate: maps (mean, variance) -> (mean, variance)
    std::unique_ptr<MLPModel> surrogate_;
    
    void unscentedTransform(const Tensor& mean, const Tensor& cov,
                           std::function<Tensor(const Tensor&)> f,
                           Tensor& out_mean, Tensor& out_cov);
};

// =============================================================================
// Uncertainty-Aware Solver Selection
// =============================================================================

/**
 * @brief Select solver based on uncertainty estimates
 */
class UncertaintyAwareSolverSelector {
public:
    struct SelectorConfig {
        // Uncertainty thresholds
        double low_uncertainty_threshold = 0.1;   // Use fast neural solver
        double high_uncertainty_threshold = 0.5;  // Use accurate classical solver
        
        // Solvers
        std::string fast_solver = "neural_fno";
        std::string accurate_solver = "classical_newton";
        std::string hybrid_solver = "neural_preconditioned";
    };
    
    UncertaintyAwareSolverSelector(const SelectorConfig& config);
    
    /**
     * @brief Select appropriate solver for problem
     */
    std::string selectSolver(const Tensor& problem_state,
                            ProbabilisticFNO& uncertainty_estimator);
    
    /**
     * @brief Adaptive solve with automatic method selection
     */
    Tensor adaptiveSolve(const Tensor& problem_state,
                        std::map<std::string, 
                            std::function<Tensor(const Tensor&)>>& solvers,
                        ProbabilisticFNO& uncertainty_estimator);
    
private:
    SelectorConfig config_;
};

// =============================================================================
// Calibration
// =============================================================================

/**
 * @brief Calibrate uncertainty estimates
 */
class UncertaintyCalibrator {
public:
    /**
     * @brief Temperature scaling for calibration
     */
    static double calibrateTemperature(
        const std::vector<Tensor>& predictions,
        const std::vector<Tensor>& uncertainties,
        const std::vector<Tensor>& ground_truths);
    
    /**
     * @brief Compute Expected Calibration Error
     */
    static double computeECE(
        const std::vector<Tensor>& predictions,
        const std::vector<Tensor>& uncertainties,
        const std::vector<Tensor>& ground_truths,
        int num_bins = 10);
    
    /**
     * @brief Reliability diagram data
     */
    static void reliabilityDiagram(
        const std::vector<Tensor>& predictions,
        const std::vector<Tensor>& uncertainties,
        const std::vector<Tensor>& ground_truths,
        std::vector<double>& confidence_bins,
        std::vector<double>& accuracy_bins,
        int num_bins = 10);
};

// =============================================================================
// Configuration Parsing
// =============================================================================

BNNConfig parseBNNConfig(const std::map<std::string, std::string>& config);
DeepEnsembleConfig parseDeepEnsembleConfig(const std::map<std::string, std::string>& config);

} // namespace ML
} // namespace FSRM

#endif // FSRM_BAYESIAN_NEURAL_OPERATOR_HPP
