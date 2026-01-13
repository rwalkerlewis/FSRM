/**
 * @file NeuralTimestepping.hpp
 * @brief Neural Network-Based Time Step Control
 * 
 * Learns optimal time stepping strategies:
 * 
 * - Adaptive time step selection
 * - IMEX transition prediction
 * - Stability monitoring
 * - Error estimation
 * 
 * Benefits:
 * - Faster simulations through optimal time steps
 * - Improved stability through learned predictions
 * - Reduced failed time steps
 */

#ifndef FSRM_NEURAL_TIMESTEPPING_HPP
#define FSRM_NEURAL_TIMESTEPPING_HPP

#include "FourierNeuralOperator.hpp"
#include "NeuralSurrogates.hpp"
#include <petsc.h>
#include <vector>
#include <memory>
#include <deque>

namespace FSRM {
namespace ML {

// =============================================================================
// Neural Time Step Predictor
// =============================================================================

/**
 * @brief Predict optimal time step size
 */
class NeuralTimeStepPredictor {
public:
    struct TimeStepConfig {
        // Feature extraction
        int state_history_length = 5;       // Past states to consider
        bool use_gradient_features = true;
        bool use_curvature_features = true;
        bool use_physics_features = true;
        
        // Constraints
        double dt_min = 1e-10;
        double dt_max = 1e3;
        double max_change_ratio = 2.0;      // Max ratio dt_new/dt_old
        
        // Network
        std::vector<int> hidden_layers = {128, 128};
        std::string activation = "relu";
        
        // Training
        double learning_rate = 1e-3;
        int epochs = 200;
        
        // Safety factor
        double safety_factor = 0.9;
    };
    
    NeuralTimeStepPredictor(const TimeStepConfig& config);
    
    /**
     * @brief Predict optimal time step
     * 
     * @param state Current solution state
     * @param time Current time
     * @param dt_prev Previous time step
     * @return Predicted optimal time step
     */
    double predictTimeStep(const Tensor& state, double time, double dt_prev);
    
    /**
     * @brief Predict with confidence
     */
    double predictTimeStepWithConfidence(const Tensor& state, double time, 
                                         double dt_prev, double& confidence);
    
    /**
     * @brief Update history with new state
     */
    void updateHistory(const Tensor& state, double time, double dt);
    
    /**
     * @brief Train from simulation data
     * 
     * @param states Solution trajectories
     * @param times Time points
     * @param successful_dts Time steps that succeeded
     * @param failed_dts Time steps that failed
     */
    void train(const std::vector<std::vector<Tensor>>& states,
              const std::vector<std::vector<double>>& times,
              const std::vector<std::vector<double>>& successful_dts,
              const std::vector<std::vector<double>>& failed_dts);
    
    /**
     * @brief Online learning from current simulation
     */
    void onlineUpdate(const Tensor& state, double dt_used, bool success);
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    TimeStepConfig config_;
    std::unique_ptr<MLPModel> predictor_;
    
    // State history
    std::deque<Tensor> state_history_;
    std::deque<double> dt_history_;
    std::deque<bool> success_history_;
    
    // Feature extraction
    Tensor extractFeatures(const Tensor& state, double time, double dt_prev);
    Tensor computeGradientFeatures(const Tensor& state);
    Tensor computeCurvatureFeatures(const Tensor& state);
    Tensor computePhysicsFeatures(const Tensor& state);
};

// =============================================================================
// Neural IMEX Transition Predictor
// =============================================================================

/**
 * @brief Predict optimal IMEX transition timing
 */
class NeuralIMEXPredictor {
public:
    struct IMEXConfig {
        // Transition triggers (classical)
        double stress_threshold = 0.8;
        double slip_rate_threshold = 1e-3;
        double velocity_threshold = 1e-4;
        
        // Neural prediction
        std::vector<int> hidden_layers = {128, 256, 128};
        bool use_lookahead = true;       // Predict future state
        int lookahead_steps = 5;
        
        // Training
        double learning_rate = 1e-3;
        int epochs = 100;
    };
    
    NeuralIMEXPredictor(const IMEXConfig& config);
    
    /**
     * @brief Predict if transition should occur
     * 
     * @param state Current state
     * @param time Current time
     * @param current_mode Current mode (0=implicit, 1=explicit)
     * @return Probability of needing transition
     */
    double predictTransitionProbability(const Tensor& state, double time, 
                                        int current_mode);
    
    /**
     * @brief Predict time to next transition
     */
    double predictTimeToTransition(const Tensor& state, double time, 
                                   int current_mode);
    
    /**
     * @brief Get recommended mode
     */
    int recommendMode(const Tensor& state, double time);
    
    /**
     * @brief Train from transition data
     */
    void train(const std::vector<Tensor>& states,
              const std::vector<double>& times,
              const std::vector<int>& modes,
              const std::vector<bool>& transitions);
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    IMEXConfig config_;
    
    // Transition classifier
    std::unique_ptr<MLPModel> transition_classifier_;
    
    // Time-to-transition regressor
    std::unique_ptr<MLPModel> time_predictor_;
    
    // State forecaster for lookahead
    std::unique_ptr<MLPModel> forecaster_;
    
    Tensor extractIMEXFeatures(const Tensor& state);
};

// =============================================================================
// Neural Stability Monitor
// =============================================================================

/**
 * @brief Monitor and predict numerical stability
 */
class NeuralStabilityMonitor {
public:
    struct StabilityConfig {
        // Stability criteria
        double cfl_target = 0.5;
        double energy_growth_threshold = 1.1;
        
        // Network
        std::vector<int> hidden_layers = {128, 128};
        
        // Detection
        int stability_check_frequency = 10;
        double instability_threshold = 0.9;
    };
    
    NeuralStabilityMonitor(const StabilityConfig& config);
    
    /**
     * @brief Check if current state is stable
     */
    bool isStable(const Tensor& state, double dt);
    
    /**
     * @brief Predict stability probability
     */
    double stabilityProbability(const Tensor& state, double dt);
    
    /**
     * @brief Predict maximum stable time step
     */
    double predictMaxStableTimeStep(const Tensor& state);
    
    /**
     * @brief Get stability features (for diagnostics)
     */
    Tensor getStabilityFeatures(const Tensor& state);
    
    /**
     * @brief Train from stability data
     */
    void train(const std::vector<Tensor>& states,
              const std::vector<double>& dts,
              const std::vector<bool>& stable_flags);
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    StabilityConfig config_;
    
    std::unique_ptr<MLPModel> stability_classifier_;
    std::unique_ptr<MLPModel> dt_regressor_;
    
    Tensor computeSpectralFeatures(const Tensor& state);
    Tensor computeGradientFeatures(const Tensor& state);
};

// =============================================================================
// Neural Error Estimator
// =============================================================================

/**
 * @brief Estimate time stepping error
 */
class NeuralErrorEstimator {
public:
    struct ErrorConfig {
        // Error types
        bool estimate_local_error = true;
        bool estimate_global_error = true;
        
        // Network
        std::vector<int> hidden_layers = {128, 256, 128};
        
        // Error targets
        double local_error_target = 1e-5;
        double global_error_target = 1e-3;
    };
    
    NeuralErrorEstimator(const ErrorConfig& config);
    
    /**
     * @brief Estimate local truncation error
     */
    double estimateLocalError(const Tensor& state, const Tensor& state_prev,
                             double dt);
    
    /**
     * @brief Estimate accumulated global error
     */
    double estimateGlobalError(const Tensor& state, double time);
    
    /**
     * @brief Suggest time step based on error target
     */
    double suggestTimeStep(const Tensor& state, double dt_current,
                          double target_error);
    
    /**
     * @brief Train from high-accuracy reference solutions
     */
    void train(const std::vector<Tensor>& states,
              const std::vector<Tensor>& reference_states,
              const std::vector<double>& dts);
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    ErrorConfig config_;
    
    std::unique_ptr<MLPModel> local_error_net_;
    std::unique_ptr<MLPModel> global_error_net_;
};

// =============================================================================
// Neural Time Integration Scheme
// =============================================================================

/**
 * @brief Learned time integration method
 */
class NeuralTimeIntegrator {
public:
    struct IntegratorConfig {
        // Base method
        std::string base_method = "implicit_euler";  // implicit_euler, bdf2, rk4
        
        // Neural correction
        bool learn_correction = true;
        std::vector<int> correction_layers = {256, 256};
        
        // Multi-step
        int num_steps = 1;  // For multi-step methods
        
        // Training
        double learning_rate = 1e-3;
        int epochs = 200;
    };
    
    NeuralTimeIntegrator(const IntegratorConfig& config);
    
    /**
     * @brief Single time step
     */
    Tensor step(const Tensor& state, 
               std::function<Tensor(const Tensor&)> residual,
               double dt);
    
    /**
     * @brief Step with correction
     */
    Tensor stepWithCorrection(const Tensor& state,
                             std::function<Tensor(const Tensor&)> residual,
                             double dt);
    
    /**
     * @brief Train correction network
     */
    void train(std::function<Tensor(const Tensor&)> residual,
              const std::vector<Tensor>& reference_trajectories,
              const std::vector<double>& dts);
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    IntegratorConfig config_;
    
    std::unique_ptr<MLPModel> correction_net_;
    
    Tensor baseStep(const Tensor& state, 
                   std::function<Tensor(const Tensor&)> residual,
                   double dt);
};

// =============================================================================
// Unified Neural Time Stepping Controller
// =============================================================================

/**
 * @brief Complete neural time stepping system
 */
class NeuralTimeSteppingController {
public:
    struct ControllerConfig {
        // Components to enable
        bool use_neural_dt_selection = true;
        bool use_neural_imex = true;
        bool use_neural_stability = true;
        bool use_neural_error = true;
        bool use_neural_integrator = false;
        
        // Component configs
        NeuralTimeStepPredictor::TimeStepConfig dt_config;
        NeuralIMEXPredictor::IMEXConfig imex_config;
        NeuralStabilityMonitor::StabilityConfig stability_config;
        NeuralErrorEstimator::ErrorConfig error_config;
        NeuralTimeIntegrator::IntegratorConfig integrator_config;
        
        // Fallback to classical
        bool allow_fallback = true;
        int max_retries = 3;
    };
    
    NeuralTimeSteppingController(const ControllerConfig& config);
    
    /**
     * @brief Get recommended time step
     */
    double getRecommendedTimeStep(const Tensor& state, double time,
                                  double dt_prev);
    
    /**
     * @brief Get recommended integration mode
     */
    int getRecommendedMode(const Tensor& state, double time);
    
    /**
     * @brief Full time step with all neural components
     */
    Tensor step(const Tensor& state,
               std::function<Tensor(const Tensor&)> residual,
               double& dt, int& mode);
    
    /**
     * @brief Report step outcome for learning
     */
    void reportOutcome(const Tensor& state, double dt, int mode, bool success);
    
    /**
     * @brief Train all components
     */
    void train(const std::vector<std::vector<Tensor>>& trajectories,
              const std::vector<std::vector<double>>& times,
              const std::vector<std::vector<double>>& dts,
              const std::vector<std::vector<int>>& modes);
    
    // Statistics
    struct TimeSteppingStats {
        int total_steps = 0;
        int successful_steps = 0;
        int failed_steps = 0;
        int mode_transitions = 0;
        double avg_dt = 0.0;
        double neural_time = 0.0;
    };
    
    const TimeSteppingStats& getStats() const { return stats_; }
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    ControllerConfig config_;
    TimeSteppingStats stats_;
    
    std::unique_ptr<NeuralTimeStepPredictor> dt_predictor_;
    std::unique_ptr<NeuralIMEXPredictor> imex_predictor_;
    std::unique_ptr<NeuralStabilityMonitor> stability_monitor_;
    std::unique_ptr<NeuralErrorEstimator> error_estimator_;
    std::unique_ptr<NeuralTimeIntegrator> integrator_;
    
    int current_mode_ = 0;
};

// =============================================================================
// Configuration Parsing
// =============================================================================

NeuralTimeStepPredictor::TimeStepConfig parseTimeStepConfig(
    const std::map<std::string, std::string>& config);
NeuralTimeSteppingController::ControllerConfig parseTimeSteppingControllerConfig(
    const std::map<std::string, std::string>& config);

} // namespace ML
} // namespace FSRM

#endif // FSRM_NEURAL_TIMESTEPPING_HPP
