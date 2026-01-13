/**
 * @file NeuralReducedOrderModel.hpp
 * @brief Neural Network-Enhanced Reduced Order Models
 * 
 * Combines traditional Proper Orthogonal Decomposition (POD) with neural
 * networks to create fast, accurate surrogate models:
 * 
 * - POD-NN: POD basis with neural network for reduced dynamics
 * - POD-LSTM: POD with LSTM for temporal evolution
 * - Autoencoder-ROM: Learned nonlinear reduction
 * - DeepONet-ROM: Deep operator network for parametric systems
 * 
 * Achieves 100-1000x speedup for parametric studies while
 * maintaining 1-5% error relative to full-order simulations.
 */

#ifndef FSRM_NEURAL_REDUCED_ORDER_MODEL_HPP
#define FSRM_NEURAL_REDUCED_ORDER_MODEL_HPP

#include "FourierNeuralOperator.hpp"
#include "NeuralSurrogates.hpp"
#include <petsc.h>
#include <vector>
#include <memory>
#include <complex>

namespace FSRM {
namespace ML {

// =============================================================================
// Forward Declarations
// =============================================================================

class PODBasis;
class PODNeuralROM;
class AutoencoderROM;
class LSTMCell;
class PODLSTM;
class DeepONetROM;
class ROMSolver;

// =============================================================================
// POD Basis Computation
// =============================================================================

/**
 * @brief Proper Orthogonal Decomposition basis
 */
class PODBasis {
public:
    PODBasis() = default;
    
    /**
     * @brief Compute POD basis from snapshots
     * 
     * @param snapshots Matrix of solution snapshots [n_dof x n_snapshots]
     * @param num_modes Number of POD modes to retain
     * @param energy_threshold Alternative: retain modes capturing this fraction of energy
     */
    void compute(const std::vector<Tensor>& snapshots, 
                int num_modes = -1,
                double energy_threshold = 0.99);
    
    /**
     * @brief Compute from PETSc vectors
     */
    void compute(const std::vector<Vec>& snapshots,
                int num_modes = -1,
                double energy_threshold = 0.99);
    
    /**
     * @brief Project high-dimensional state to reduced coordinates
     */
    Tensor project(const Tensor& full_state) const;
    void project(Vec full_state, Tensor& reduced) const;
    
    /**
     * @brief Reconstruct full state from reduced coordinates
     */
    Tensor reconstruct(const Tensor& reduced_state) const;
    void reconstruct(const Tensor& reduced_state, Vec full_state) const;
    
    /**
     * @brief Get projection and reconstruction errors
     */
    double getProjectionError(const Tensor& full_state) const;
    
    // Accessors
    int numModes() const { return num_modes_; }
    int fullDimension() const { return n_dof_; }
    double energyCaptured() const { return energy_captured_; }
    
    const Tensor& getBasis() const { return basis_; }
    const Tensor& getSingularValues() const { return singular_values_; }
    const Tensor& getMean() const { return mean_; }
    
    // Save/Load
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    Tensor basis_;           // [n_dof x num_modes]
    Tensor singular_values_; // [num_modes]
    Tensor mean_;            // [n_dof] - mean state for centering
    
    int n_dof_ = 0;
    int num_modes_ = 0;
    int n_snapshots_ = 0;
    double energy_captured_ = 0.0;
    
    void svd(const Tensor& data, Tensor& U, Tensor& S, Tensor& Vt);
};

// =============================================================================
// POD-Neural Network ROM
// =============================================================================

/**
 * @brief Configuration for POD-NN ROM
 */
struct PODNNConfig {
    // POD settings
    int num_modes = 50;
    double energy_threshold = 0.999;
    
    // Neural network for reduced dynamics
    std::vector<int> hidden_layers = {256, 256, 256};
    std::string activation = "gelu";
    bool use_batch_norm = true;
    
    // Time integration in reduced space
    std::string time_integrator = "implicit_euler";  // implicit_euler, rk4
    
    // Parameter dependence
    int num_parameters = 0;     // Physical parameters affecting solution
    bool parameter_aware = false;
    
    // Training
    double learning_rate = 1e-3;
    int batch_size = 64;
    int epochs = 200;
    
    // Regularization
    double weight_decay = 1e-4;
    bool use_physics_loss = false;
    double physics_weight = 0.1;
};

/**
 * @brief POD-NN Reduced Order Model
 * 
 * Projects high-dimensional dynamics to POD subspace, then uses
 * neural network to learn the reduced dynamics.
 * 
 * Full model: dU/dt = F(U, μ)
 * Reduced: da/dt = f_NN(a, μ)
 * 
 * where U = Φa + U_mean, and Φ is the POD basis.
 */
class PODNeuralROM {
public:
    PODNeuralROM(const PODNNConfig& config);
    
    /**
     * @brief Build ROM from simulation snapshots
     */
    void buildFromSnapshots(const std::vector<Tensor>& snapshots,
                           const std::vector<double>& times,
                           const std::vector<Tensor>& parameters = {});
    
    /**
     * @brief Solve in reduced space
     * 
     * @param initial_condition Full-dimensional initial condition
     * @param t_final Final time
     * @param dt Time step
     * @param parameters Physical parameters
     * @return Tensor Full-dimensional solution at t_final
     */
    Tensor solve(const Tensor& initial_condition, 
                double t_final, double dt,
                const Tensor& parameters = Tensor());
    
    /**
     * @brief Get full trajectory
     */
    std::vector<Tensor> trajectory(const Tensor& initial_condition,
                                   double t_final, double dt,
                                   const Tensor& parameters = Tensor());
    
    /**
     * @brief Evaluate reduced dynamics (for debugging/analysis)
     */
    Tensor evaluateReducedDynamics(const Tensor& reduced_state,
                                   double time,
                                   const Tensor& parameters = Tensor());
    
    /**
     * @brief Train the neural network for reduced dynamics
     */
    void train(const std::vector<Tensor>& reduced_states,
              const std::vector<Tensor>& reduced_velocities,
              const std::vector<Tensor>& parameters = {});
    
    // Access to POD basis
    const PODBasis& getBasis() const { return pod_basis_; }
    
    // Validation
    double validate(const std::vector<Tensor>& test_snapshots,
                   const std::vector<double>& test_times);
    
    // Save/Load
    void save(const std::string& path) const;
    void load(const std::string& path);
    
    // Info
    size_t numParameters() const;
    double getSpeedup() const { return speedup_factor_; }
    
private:
    PODNNConfig config_;
    PODBasis pod_basis_;
    
    // Neural network for reduced dynamics
    std::unique_ptr<MLPModel> dynamics_network_;
    
    // For parameter-dependent problems
    std::unique_ptr<MLPModel> parameter_encoder_;
    
    // Performance metrics
    double speedup_factor_ = 1.0;
    double training_error_ = 0.0;
    double validation_error_ = 0.0;
    
    // Time stepping in reduced space
    Tensor implicitEulerStep(const Tensor& a, double dt, const Tensor& params);
    Tensor rk4Step(const Tensor& a, double t, double dt, const Tensor& params);
};

// =============================================================================
// LSTM Cell for Temporal Learning
// =============================================================================

/**
 * @brief Long Short-Term Memory cell
 */
class LSTMCell {
public:
    LSTMCell(int input_size, int hidden_size);
    
    /**
     * @brief Single LSTM step
     * 
     * @param input Current input
     * @param h_prev Previous hidden state
     * @param c_prev Previous cell state
     * @return pair<h_new, c_new>
     */
    std::pair<Tensor, Tensor> forward(const Tensor& input,
                                      const Tensor& h_prev,
                                      const Tensor& c_prev);
    
    std::pair<Tensor, Tensor> backward(const Tensor& grad_h,
                                       const Tensor& grad_c);
    
    std::vector<Tensor*> parameters();
    std::vector<Tensor*> gradients();
    void zeroGrad();
    
private:
    int input_size_, hidden_size_;
    
    // Gates weights: [input; hidden] -> gate
    Tensor W_f_, U_f_, b_f_;  // Forget gate
    Tensor W_i_, U_i_, b_i_;  // Input gate
    Tensor W_o_, U_o_, b_o_;  // Output gate
    Tensor W_c_, U_c_, b_c_;  // Cell candidate
    
    // Gradients
    Tensor grad_W_f_, grad_U_f_, grad_b_f_;
    Tensor grad_W_i_, grad_U_i_, grad_b_i_;
    Tensor grad_W_o_, grad_U_o_, grad_b_o_;
    Tensor grad_W_c_, grad_U_c_, grad_b_c_;
    
    // Cached for backward
    Tensor cached_input_, cached_h_, cached_c_;
    Tensor cached_f_, cached_i_, cached_o_, cached_c_hat_;
    
    Tensor sigmoid(const Tensor& x);
    Tensor tanh(const Tensor& x);
};

/**
 * @brief POD-LSTM Reduced Order Model
 * 
 * Uses LSTM to learn temporal evolution in reduced space,
 * suitable for capturing long-term dependencies.
 */
class PODLSTM {
public:
    struct PODLSTMConfig {
        int num_modes = 50;
        int lstm_hidden_size = 128;
        int num_lstm_layers = 2;
        int sequence_length = 10;  // History for prediction
        
        double learning_rate = 1e-3;
        int epochs = 100;
        int batch_size = 32;
    };
    
    PODLSTM(const PODLSTMConfig& config);
    
    /**
     * @brief Train on trajectories
     */
    void train(const std::vector<std::vector<Tensor>>& trajectories);
    
    /**
     * @brief Predict next state given history
     */
    Tensor predictNext(const std::vector<Tensor>& history);
    
    /**
     * @brief Generate trajectory
     */
    std::vector<Tensor> generate(const std::vector<Tensor>& initial_sequence,
                                 int num_steps);
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    PODLSTMConfig config_;
    PODBasis pod_basis_;
    
    std::vector<std::unique_ptr<LSTMCell>> lstm_layers_;
    std::unique_ptr<DenseLayer> output_layer_;
    
    std::vector<Tensor> h_states_;
    std::vector<Tensor> c_states_;
};

// =============================================================================
// Autoencoder-Based ROM
// =============================================================================

/**
 * @brief Configuration for Autoencoder ROM
 */
struct AutoencoderROMConfig {
    // Encoder architecture
    std::vector<int> encoder_layers = {1000, 500, 100, 20};
    
    // Latent dimension
    int latent_dim = 20;
    
    // Decoder (mirrored encoder by default)
    std::vector<int> decoder_layers;  // Empty = mirror encoder
    
    // Activation
    std::string activation = "gelu";
    
    // Variational autoencoder
    bool variational = false;
    double kl_weight = 1e-3;
    
    // Dynamics network
    std::vector<int> dynamics_layers = {64, 64};
    
    // Training
    double learning_rate = 1e-3;
    int epochs = 500;
    int batch_size = 64;
};

/**
 * @brief Autoencoder-based nonlinear reduction
 * 
 * Learns nonlinear low-dimensional manifold via autoencoder,
 * then learns dynamics on this manifold.
 */
class AutoencoderROM {
public:
    AutoencoderROM(const AutoencoderROMConfig& config);
    
    /**
     * @brief Train autoencoder on snapshots
     */
    void trainAutoencoder(const std::vector<Tensor>& snapshots);
    
    /**
     * @brief Train dynamics in latent space
     */
    void trainDynamics(const std::vector<Tensor>& snapshots,
                      const std::vector<double>& times);
    
    /**
     * @brief End-to-end training
     */
    void train(const std::vector<Tensor>& snapshots,
              const std::vector<double>& times);
    
    /**
     * @brief Encode to latent space
     */
    Tensor encode(const Tensor& full_state);
    
    /**
     * @brief Decode from latent space
     */
    Tensor decode(const Tensor& latent_state);
    
    /**
     * @brief Solve in latent space
     */
    Tensor solve(const Tensor& initial_condition, 
                double t_final, double dt);
    
    /**
     * @brief Get latent trajectory
     */
    std::vector<Tensor> latentTrajectory(const Tensor& initial_latent,
                                         double t_final, double dt);
    
    // Reconstruction error
    double reconstructionError(const std::vector<Tensor>& test_data);
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    AutoencoderROMConfig config_;
    
    std::unique_ptr<MLPModel> encoder_;
    std::unique_ptr<MLPModel> decoder_;
    std::unique_ptr<MLPModel> dynamics_;
    
    // For VAE
    std::unique_ptr<MLPModel> mu_layer_;
    std::unique_ptr<MLPModel> logvar_layer_;
    
    Tensor reparameterize(const Tensor& mu, const Tensor& logvar);
};

// =============================================================================
// Deep Operator Network ROM
// =============================================================================

/**
 * @brief DeepONet for parametric PDEs
 * 
 * Learns operator mapping: (initial_condition, parameters, location) -> solution
 * 
 * Decomposition: G(u)(y) = sum_k b_k(u) * t_k(y)
 * where b is the branch net (encodes input function)
 * and t is the trunk net (encodes evaluation location)
 */
class DeepONetROM {
public:
    struct DeepONetConfig {
        // Branch network (input function encoder)
        int num_sensors = 100;      // Points to sample input function
        std::vector<int> branch_layers = {100, 128, 128, 64};
        
        // Trunk network (location encoder)
        int coord_dim = 3;          // Spatial dimensions
        std::vector<int> trunk_layers = {3, 128, 128, 64};
        
        // Number of basis functions
        int num_basis = 64;
        
        // Training
        double learning_rate = 1e-3;
        int epochs = 200;
        int batch_size = 256;
    };
    
    DeepONetROM(const DeepONetConfig& config);
    
    /**
     * @brief Evaluate operator at locations
     * 
     * @param input_function Values at sensor locations
     * @param query_locations Where to evaluate output
     * @return Output values at query locations
     */
    Tensor evaluate(const Tensor& input_function,
                   const Tensor& query_locations);
    
    /**
     * @brief Solve PDE given initial/boundary conditions
     */
    Tensor solve(const Tensor& boundary_values,
                const Tensor& domain_points);
    
    /**
     * @brief Train from input-output pairs
     */
    void train(const std::vector<Tensor>& input_functions,
              const std::vector<Tensor>& output_functions,
              const std::vector<Tensor>& locations);
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    DeepONetConfig config_;
    
    std::unique_ptr<MLPModel> branch_net_;  // Encodes input function
    std::unique_ptr<MLPModel> trunk_net_;   // Encodes evaluation location
    
    Tensor branch_bias_;  // Bias term
};

// =============================================================================
// Unified ROM Solver
// =============================================================================

/**
 * @brief Unified interface for all ROM types
 */
class ROMSolver {
public:
    enum class ROMType {
        POD_NN,
        POD_LSTM,
        AUTOENCODER,
        DEEPONET
    };
    
    struct ROMConfig {
        ROMType type = ROMType::POD_NN;
        
        // Common settings
        int latent_dim = 50;
        double energy_threshold = 0.999;
        
        // Type-specific configs (use appropriate one)
        PODNNConfig pod_nn_config;
        PODLSTM::PODLSTMConfig pod_lstm_config;
        AutoencoderROMConfig autoencoder_config;
        DeepONetROM::DeepONetConfig deeponet_config;
    };
    
    ROMSolver(const ROMConfig& config);
    
    /**
     * @brief Build ROM from snapshots
     */
    void build(const std::vector<Vec>& snapshots,
              const std::vector<double>& times,
              const std::vector<double>& parameters = {});
    
    /**
     * @brief Solve with ROM
     */
    PetscErrorCode solve(Vec initial, double t_final, Vec result);
    
    /**
     * @brief Get full trajectory
     */
    PetscErrorCode trajectory(Vec initial, double t_final, double dt,
                             std::vector<Vec>& results);
    
    /**
     * @brief Online error indicator
     */
    double errorIndicator(Vec rom_solution);
    
    // Model I/O
    void saveModel(const std::string& path) const;
    void loadModel(const std::string& path);
    
    // Statistics
    struct ROMStats {
        double compression_ratio = 1.0;
        double speedup = 1.0;
        double training_error = 0.0;
        double validation_error = 0.0;
        int latent_dimension = 0;
        int full_dimension = 0;
        double energy_captured = 0.0;
    };
    
    const ROMStats& getStats() const { return stats_; }
    
private:
    ROMConfig config_;
    ROMStats stats_;
    
    // ROM implementations
    std::unique_ptr<PODNeuralROM> pod_nn_;
    std::unique_ptr<PODLSTM> pod_lstm_;
    std::unique_ptr<AutoencoderROM> autoencoder_;
    std::unique_ptr<DeepONetROM> deeponet_;
    
    // Conversion helpers
    Tensor vecToTensor(Vec v) const;
    void tensorToVec(const Tensor& t, Vec v) const;
};

// =============================================================================
// Hyper-Reduction for Nonlinear Terms
// =============================================================================

/**
 * @brief Discrete Empirical Interpolation Method (DEIM)
 * 
 * Reduces cost of evaluating nonlinear terms in ROM.
 */
class DEIMHyperReduction {
public:
    DEIMHyperReduction() = default;
    
    /**
     * @brief Compute DEIM interpolation points and basis
     * 
     * @param nonlinear_snapshots Snapshots of nonlinear term evaluations
     * @param num_points Number of interpolation points
     */
    void compute(const std::vector<Tensor>& nonlinear_snapshots,
                int num_points);
    
    /**
     * @brief Get interpolation indices
     */
    const std::vector<int>& getIndices() const { return indices_; }
    
    /**
     * @brief Approximate full nonlinear term from samples
     */
    Tensor approximate(const Tensor& sampled_values);
    
    /**
     * @brief Get the DEIM basis
     */
    const Tensor& getBasis() const { return basis_; }
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    std::vector<int> indices_;     // Interpolation point indices
    Tensor basis_;                  // DEIM basis vectors
    Tensor interpolation_matrix_;   // For reconstruction
};

// =============================================================================
// Configuration Parsing
// =============================================================================

PODNNConfig parsePODNNConfig(const std::map<std::string, std::string>& config);
AutoencoderROMConfig parseAutoencoderConfig(const std::map<std::string, std::string>& config);

} // namespace ML
} // namespace FSRM

#endif // FSRM_NEURAL_REDUCED_ORDER_MODEL_HPP
