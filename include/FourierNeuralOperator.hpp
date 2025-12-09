/**
 * @file FourierNeuralOperator.hpp
 * @brief Fourier Neural Operator (FNO) for learning PDE solutions
 * 
 * FNO learns operators between function spaces by parameterizing
 * the integral kernel in Fourier space. This enables efficient
 * learning of complex PDE solutions with mesh-independent inference.
 * 
 * References:
 * - Li et al. (2020). "Fourier Neural Operator for Parametric Partial Differential Equations"
 * - Li et al. (2021). "Neural Operator: Graph Kernel Network for Partial Differential Equations"
 * 
 * Features:
 * - 1D, 2D, and 3D FNO architectures
 * - Multiple Fourier modes for multi-scale learning
 * - GPU acceleration (CUDA/ROCm)
 * - Pre-trained model loading
 * - Online training with simulation data
 * - Hybrid numerical-FNO solving
 */

#ifndef FSRM_FOURIER_NEURAL_OPERATOR_HPP
#define FSRM_FOURIER_NEURAL_OPERATOR_HPP

#include <vector>
#include <array>
#include <complex>
#include <memory>
#include <string>
#include <functional>
#include <random>
#include <fstream>
#include <map>

namespace FSRM {
namespace ML {

// =============================================================================
// Forward Declarations
// =============================================================================

class FNOLayer;
class FNOModel;
class FNOTrainer;
class FNOSolver;

// =============================================================================
// Configuration Structures
// =============================================================================

/**
 * @brief Solver method selection
 */
enum class SolverMethod {
    NUMERICAL,          // Traditional numerical solver (DG, FEM, FD)
    FNO,               // Pure Fourier Neural Operator
    HYBRID_FNO_REFINE, // FNO prediction + numerical refinement
    HYBRID_COARSE_FNO, // Coarse numerical + FNO upscaling
    ENSEMBLE           // Ensemble of numerical and FNO predictions
};

/**
 * @brief FNO architecture configuration
 */
struct FNOConfig {
    // Architecture
    int input_channels = 3;       // Number of input fields
    int output_channels = 3;      // Number of output fields
    int hidden_channels = 64;     // Width of hidden layers
    int num_layers = 4;           // Number of FNO layers
    int num_modes = 16;           // Number of Fourier modes per dimension
    
    // Spatial dimensions
    int spatial_dim = 2;          // 1, 2, or 3
    std::vector<int> resolution;  // Grid resolution [nx, ny, nz]
    
    // Training
    double learning_rate = 1e-3;
    double weight_decay = 1e-4;
    int batch_size = 16;
    int epochs = 100;
    double validation_split = 0.1;
    
    // Activation and normalization
    std::string activation = "gelu";      // gelu, relu, tanh, silu
    bool use_layer_norm = true;
    bool use_residual = true;
    
    // GPU settings
    bool use_gpu = true;
    int gpu_device_id = 0;
    
    // Model I/O
    std::string model_path = "";          // Path to pre-trained model
    bool load_pretrained = false;
    
    // Hybrid settings
    double fno_weight = 0.8;              // Weight for FNO in ensemble
    int refinement_iterations = 3;        // For hybrid refinement
};

/**
 * @brief Training data sample
 */
struct FNOSample {
    std::vector<double> input;    // Input fields (flattened)
    std::vector<double> output;   // Target output fields (flattened)
    std::vector<double> params;   // Optional physics parameters
};

/**
 * @brief Training dataset
 */
struct FNODataset {
    std::vector<FNOSample> samples;
    std::vector<int> input_shape;   // [channels, nx, ny, nz]
    std::vector<int> output_shape;  // [channels, nx, ny, nz]
    
    void addSample(const FNOSample& sample) { samples.push_back(sample); }
    void shuffle(std::mt19937& rng);
    void split(double ratio, FNODataset& train, FNODataset& val) const;
    size_t size() const { return samples.size(); }
};

// =============================================================================
// Tensor Class (Simple implementation for CPU, extended for GPU)
// =============================================================================

/**
 * @brief Simple tensor class for FNO operations
 */
class Tensor {
public:
    std::vector<double> data;
    std::vector<int> shape;
    bool on_gpu = false;
    
    Tensor() = default;
    Tensor(const std::vector<int>& shape);
    Tensor(const std::vector<int>& shape, double value);
    Tensor(const std::vector<int>& shape, const std::vector<double>& data);
    
    // Basic operations
    size_t size() const;
    size_t numel() const;
    int dim() const { return shape.size(); }
    
    // Element access
    double& operator()(int i);
    double& operator()(int i, int j);
    double& operator()(int i, int j, int k);
    double& operator()(int i, int j, int k, int l);
    
    const double& operator()(int i) const;
    const double& operator()(int i, int j) const;
    const double& operator()(int i, int j, int k) const;
    const double& operator()(int i, int j, int k, int l) const;
    
    // Linear index
    double& at(size_t idx) { return data[idx]; }
    const double& at(size_t idx) const { return data[idx]; }
    
    // Reshape
    Tensor reshape(const std::vector<int>& new_shape) const;
    Tensor transpose(int dim1, int dim2) const;
    
    // Math operations
    Tensor operator+(const Tensor& other) const;
    Tensor operator-(const Tensor& other) const;
    Tensor operator*(const Tensor& other) const;  // Element-wise
    Tensor operator*(double scalar) const;
    Tensor& operator+=(const Tensor& other);
    
    // Reductions
    double sum() const;
    double mean() const;
    double max() const;
    double min() const;
    double norm() const;
    
    // GPU transfer
    void toGPU();
    void toCPU();
    
    // Initialization
    void zeros();
    void ones();
    void randn(double mean = 0.0, double std = 1.0);
    void xavier_init();
    void kaiming_init();
    
    // Serialization
    void save(const std::string& path) const;
    void load(const std::string& path);
};

// Complex tensor for Fourier operations
class ComplexTensor {
public:
    std::vector<std::complex<double>> data;
    std::vector<int> shape;
    
    ComplexTensor() = default;
    ComplexTensor(const std::vector<int>& shape);
    
    size_t numel() const;
    
    std::complex<double>& operator()(int i, int j, int k, int l);
    const std::complex<double>& operator()(int i, int j, int k, int l) const;
    
    ComplexTensor operator*(const ComplexTensor& other) const;  // Element-wise
    ComplexTensor& operator+=(const ComplexTensor& other);
};

// =============================================================================
// Activation Functions
// =============================================================================

namespace Activation {
    Tensor relu(const Tensor& x);
    Tensor gelu(const Tensor& x);
    Tensor silu(const Tensor& x);
    Tensor tanh(const Tensor& x);
    
    Tensor relu_backward(const Tensor& x, const Tensor& grad);
    Tensor gelu_backward(const Tensor& x, const Tensor& grad);
    Tensor silu_backward(const Tensor& x, const Tensor& grad);
    Tensor tanh_backward(const Tensor& x, const Tensor& grad);
    
    using ActivationFn = std::function<Tensor(const Tensor&)>;
    using ActivationBackwardFn = std::function<Tensor(const Tensor&, const Tensor&)>;
    
    ActivationFn getActivation(const std::string& name);
    ActivationBackwardFn getActivationBackward(const std::string& name);
}

// =============================================================================
// FFT Operations
// =============================================================================

/**
 * @brief FFT utilities for FNO
 */
class FFT {
public:
    // 1D FFT
    static void fft1d(const double* input, std::complex<double>* output, int n);
    static void ifft1d(const std::complex<double>* input, double* output, int n);
    
    // 2D FFT
    static void fft2d(const Tensor& input, ComplexTensor& output);
    static void ifft2d(const ComplexTensor& input, Tensor& output);
    
    // 3D FFT
    static void fft3d(const Tensor& input, ComplexTensor& output);
    static void ifft3d(const ComplexTensor& input, Tensor& output);
    
    // Real FFT (more efficient for real inputs)
    static void rfft2d(const Tensor& input, ComplexTensor& output);
    static void irfft2d(const ComplexTensor& input, Tensor& output, const std::vector<int>& shape);
    
    // Batched FFT (for batch x channels x spatial)
    static void batchedRfft2d(const Tensor& input, ComplexTensor& output);
    static void batchedIrfft2d(const ComplexTensor& input, Tensor& output, const std::vector<int>& shape);
};

// =============================================================================
// FNO Layers
// =============================================================================

/**
 * @brief Spectral convolution layer
 * 
 * Performs convolution in Fourier space by element-wise multiplication
 * with learnable complex weights.
 */
class SpectralConv2d {
public:
    SpectralConv2d(int in_channels, int out_channels, int modes1, int modes2);
    
    // Forward pass
    Tensor forward(const Tensor& x);
    
    // Backward pass
    Tensor backward(const Tensor& grad_output);
    
    // Parameters
    ComplexTensor weights;  // [in_ch, out_ch, modes1, modes2]
    ComplexTensor grad_weights;
    
    int in_channels, out_channels;
    int modes1, modes2;
    
    // For backward pass
    ComplexTensor last_x_ft;
    std::vector<int> last_input_shape;
    
private:
    ComplexTensor complex_mul(const ComplexTensor& a, const ComplexTensor& b);
};

/**
 * @brief 3D Spectral convolution for volumetric data
 */
class SpectralConv3d {
public:
    SpectralConv3d(int in_channels, int out_channels, int modes1, int modes2, int modes3);
    
    Tensor forward(const Tensor& x);
    Tensor backward(const Tensor& grad_output);
    
    ComplexTensor weights;
    ComplexTensor grad_weights;
    
    int in_channels, out_channels;
    int modes1, modes2, modes3;
};

/**
 * @brief Linear layer (1x1 convolution)
 */
class Linear {
public:
    Linear(int in_features, int out_features, bool bias = true);
    
    Tensor forward(const Tensor& x);
    Tensor backward(const Tensor& grad_output);
    
    Tensor weight;      // [out, in]
    Tensor bias_vec;    // [out]
    Tensor grad_weight;
    Tensor grad_bias;
    
    bool has_bias;
    Tensor last_input;
};

/**
 * @brief Layer normalization
 */
class LayerNorm {
public:
    LayerNorm(int normalized_shape, double eps = 1e-5);
    
    Tensor forward(const Tensor& x);
    Tensor backward(const Tensor& grad_output);
    
    Tensor gamma;  // Scale
    Tensor beta;   // Shift
    Tensor grad_gamma;
    Tensor grad_beta;
    
    int normalized_shape;
    double eps;
    
    // Cached for backward
    Tensor x_normalized;
    Tensor std_inv;
};

/**
 * @brief Single FNO layer (spectral conv + linear + activation)
 */
class FNOLayer {
public:
    FNOLayer(int channels, int modes1, int modes2, 
             const std::string& activation = "gelu",
             bool use_layer_norm = true,
             bool use_residual = true);
    
    Tensor forward(const Tensor& x);
    Tensor backward(const Tensor& grad_output);
    
    // Components
    std::unique_ptr<SpectralConv2d> spectral_conv;
    std::unique_ptr<Linear> linear;
    std::unique_ptr<LayerNorm> norm;
    
    Activation::ActivationFn activation_fn;
    Activation::ActivationBackwardFn activation_backward;
    
    bool use_layer_norm;
    bool use_residual;
    
    // Cached for backward
    Tensor pre_activation;
    Tensor input_cache;
};

// =============================================================================
// FNO Model
// =============================================================================

/**
 * @brief Complete FNO model for learning PDE operators
 */
class FNOModel {
public:
    FNOModel(const FNOConfig& config);
    
    // Forward pass
    Tensor forward(const Tensor& x);
    
    // Backward pass (returns loss gradient w.r.t input)
    Tensor backward(const Tensor& grad_output);
    
    // Prediction (no grad computation)
    Tensor predict(const Tensor& x);
    
    // Loss computation
    double computeLoss(const Tensor& prediction, const Tensor& target);
    Tensor computeLossGrad(const Tensor& prediction, const Tensor& target);
    
    // Parameter access
    std::vector<Tensor*> parameters();
    std::vector<Tensor*> gradients();
    void zeroGrad();
    
    // Model I/O
    void save(const std::string& path) const;
    void load(const std::string& path);
    void exportONNX(const std::string& path) const;
    
    // GPU
    void toGPU();
    void toCPU();
    bool isOnGPU() const { return on_gpu; }
    
    // Info
    size_t numParameters() const;
    void summary() const;
    
private:
    FNOConfig config;
    
    // Layers
    std::unique_ptr<Linear> lifting;      // Input projection
    std::vector<std::unique_ptr<FNOLayer>> fno_layers;
    std::unique_ptr<Linear> projection1;  // Output MLP
    std::unique_ptr<Linear> projection2;
    
    bool on_gpu = false;
    
    // Cached tensors for backward
    Tensor lifted;
    std::vector<Tensor> layer_outputs;
    Tensor proj1_out;
};

// =============================================================================
// Optimizer
// =============================================================================

/**
 * @brief Adam optimizer
 */
class AdamOptimizer {
public:
    AdamOptimizer(std::vector<Tensor*> params, double lr = 1e-3,
                  double beta1 = 0.9, double beta2 = 0.999,
                  double eps = 1e-8, double weight_decay = 0.0);
    
    void step(std::vector<Tensor*> grads);
    void zeroGrad();
    
    void setLearningRate(double lr) { this->lr = lr; }
    double getLearningRate() const { return lr; }
    
private:
    std::vector<Tensor*> params;
    std::vector<Tensor> m;  // First moment
    std::vector<Tensor> v;  // Second moment
    
    double lr, beta1, beta2, eps, weight_decay;
    int t = 0;  // Time step
};

/**
 * @brief Learning rate scheduler
 */
class LRScheduler {
public:
    enum class Type {
        CONSTANT,
        STEP,
        EXPONENTIAL,
        COSINE_ANNEALING,
        REDUCE_ON_PLATEAU
    };
    
    LRScheduler(AdamOptimizer& optimizer, Type type, double initial_lr);
    
    void step(double metric = 0.0);
    double getCurrentLR() const;
    
    // Type-specific settings
    void setStepSchedule(int step_size, double gamma = 0.1);
    void setExponentialSchedule(double gamma);
    void setCosineSchedule(int T_max, double eta_min = 0.0);
    void setReduceOnPlateau(double factor = 0.1, int patience = 10);
    
private:
    AdamOptimizer& optimizer;
    Type type;
    double initial_lr, current_lr;
    int epoch = 0;
    
    // Schedule parameters
    int step_size = 30;
    double gamma = 0.1;
    int T_max = 100;
    double eta_min = 0.0;
    int patience = 10;
    double best_metric = 1e10;
    int patience_counter = 0;
};

// =============================================================================
// FNO Trainer
// =============================================================================

/**
 * @brief Training loop for FNO models
 */
class FNOTrainer {
public:
    FNOTrainer(FNOModel& model, const FNOConfig& config);
    
    // Training
    void train(const FNODataset& train_data, const FNODataset& val_data);
    double trainEpoch(const FNODataset& data);
    double validate(const FNODataset& data);
    
    // Callbacks
    using EpochCallback = std::function<void(int epoch, double train_loss, double val_loss)>;
    void setCallback(EpochCallback cb) { callback = cb; }
    
    // Metrics
    struct TrainingHistory {
        std::vector<double> train_losses;
        std::vector<double> val_losses;
        std::vector<double> learning_rates;
    };
    
    const TrainingHistory& getHistory() const { return history; }
    
    // Early stopping
    void setEarlyStopping(int patience, double min_delta = 1e-4);
    
private:
    FNOModel& model;
    FNOConfig config;
    
    std::unique_ptr<AdamOptimizer> optimizer;
    std::unique_ptr<LRScheduler> scheduler;
    
    TrainingHistory history;
    EpochCallback callback;
    
    // Early stopping
    bool use_early_stopping = false;
    int es_patience = 10;
    double es_min_delta = 1e-4;
    int es_counter = 0;
    double best_val_loss = 1e10;
    
    std::mt19937 rng;
};

// =============================================================================
// Data Generation
// =============================================================================

/**
 * @brief Generate training data from numerical simulations
 */
class FNODataGenerator {
public:
    FNODataGenerator(const FNOConfig& config);
    
    // Generate samples from numerical solver
    void generateFromSimulation(
        std::function<Tensor(const Tensor&)> numerical_solver,
        std::function<Tensor()> initial_condition_sampler,
        int num_samples);
    
    // Add existing simulation data
    void addSimulationResult(const Tensor& input, const Tensor& output);
    void addSimulationResult(const Tensor& input, const Tensor& output,
                            const std::vector<double>& params);
    
    // Access dataset
    FNODataset& getDataset() { return dataset; }
    const FNODataset& getDataset() const { return dataset; }
    
    // Save/load dataset
    void saveDataset(const std::string& path) const;
    void loadDataset(const std::string& path);
    
    // Data augmentation
    void enableAugmentation(bool flip = true, bool rotate = true, bool noise = true);
    
private:
    FNOConfig config;
    FNODataset dataset;
    
    bool augment_flip = false;
    bool augment_rotate = false;
    bool augment_noise = false;
    double noise_std = 0.01;
    
    std::mt19937 rng;
    
    void augmentSample(FNOSample& sample);
};

// =============================================================================
// FNO Solver (Main Interface)
// =============================================================================

/**
 * @brief Main solver interface combining numerical and FNO approaches
 */
class FNOSolver {
public:
    FNOSolver(const FNOConfig& config);
    
    // Set solver method
    void setMethod(SolverMethod method);
    SolverMethod getMethod() const { return method; }
    
    // Set numerical solver for hybrid modes
    using NumericalSolverFn = std::function<Tensor(const Tensor&, double)>;
    void setNumericalSolver(NumericalSolverFn solver);
    
    // Main solve interface
    Tensor solve(const Tensor& initial_condition, double t_final);
    
    // Step-by-step solving
    Tensor step(const Tensor& current_state, double dt);
    
    // Training interface
    void trainFromData(const FNODataset& dataset);
    void trainOnline(const Tensor& input, const Tensor& target);
    
    // Model management
    void loadModel(const std::string& path);
    void saveModel(const std::string& path) const;
    bool isModelTrained() const { return model_trained; }
    
    // Get FNO model for direct access
    FNOModel& getModel() { return *model; }
    const FNOModel& getModel() const { return *model; }
    
    // Get trainer for advanced training control
    FNOTrainer& getTrainer() { return *trainer; }
    
    // Configuration
    void setHybridWeight(double weight) { config.fno_weight = weight; }
    void setRefinementIterations(int iters) { config.refinement_iterations = iters; }
    
    // Performance metrics
    struct SolverStats {
        double fno_time = 0.0;
        double numerical_time = 0.0;
        double total_time = 0.0;
        int fno_calls = 0;
        int numerical_calls = 0;
        double avg_error_estimate = 0.0;
    };
    
    const SolverStats& getStats() const { return stats; }
    void resetStats() { stats = SolverStats(); }
    
private:
    FNOConfig config;
    SolverMethod method = SolverMethod::NUMERICAL;
    
    std::unique_ptr<FNOModel> model;
    std::unique_ptr<FNOTrainer> trainer;
    std::unique_ptr<FNODataGenerator> data_gen;
    
    NumericalSolverFn numerical_solver;
    
    bool model_trained = false;
    SolverStats stats;
    
    // Hybrid solving methods
    Tensor solveHybridRefine(const Tensor& initial_condition, double t_final);
    Tensor solveHybridCoarse(const Tensor& initial_condition, double t_final);
    Tensor solveEnsemble(const Tensor& initial_condition, double t_final);
    
    // Error estimation
    double estimateError(const Tensor& fno_result, const Tensor& numerical_result);
};

// =============================================================================
// Physics-Informed FNO Extensions
// =============================================================================

/**
 * @brief Physics-informed loss functions
 */
class PhysicsInformedLoss {
public:
    // Compute physics residual for wave equation
    static Tensor waveEquationResidual(const Tensor& u, const Tensor& u_prev,
                                       double c, double dx, double dt);
    
    // Compute physics residual for diffusion equation
    static Tensor diffusionResidual(const Tensor& u, const Tensor& u_prev,
                                    double D, double dx, double dt);
    
    // Compute physics residual for elastic wave
    static Tensor elasticWaveResidual(const Tensor& u, const Tensor& stress,
                                      double rho, double lambda, double mu,
                                      double dx, double dt);
    
    // Combined data + physics loss
    static double computeLoss(const Tensor& prediction, const Tensor& target,
                             const Tensor& physics_residual,
                             double data_weight = 1.0, double physics_weight = 0.1);
};

/**
 * @brief Super-resolution FNO for upscaling coarse solutions
 */
class SuperResolutionFNO {
public:
    SuperResolutionFNO(int scale_factor, const FNOConfig& base_config);
    
    Tensor upscale(const Tensor& coarse_solution);
    void train(const FNODataset& paired_data);  // Pairs of (coarse, fine)
    
private:
    int scale_factor;
    std::unique_ptr<FNOModel> model;
    FNOConfig config;
};

// =============================================================================
// Extended Physics-Informed Neural Operator (PINN Extensions)
// =============================================================================

/**
 * @brief Configuration for Physics-Informed FNO
 */
struct PINNConfig : FNOConfig {
    // Physics equation type
    enum class PhysicsType {
        NONE,               // No physics constraints
        POISSON,            // -∇²u = f
        DIFFUSION,          // ∂u/∂t = D∇²u
        ADVECTION,          // ∂u/∂t + v·∇u = 0
        WAVE,               // ∂²u/∂t² = c²∇²u
        NAVIER_STOKES,      // Incompressible NS
        ELASTODYNAMICS,     // Elastic wave equation
        DARCY,              // Darcy flow in porous media
        RESERVOIR,          // Reservoir flow (multiphase)
        COUPLED             // Coupled multi-physics
    };
    PhysicsType physics_type = PhysicsType::NONE;
    
    // Loss weights
    double data_loss_weight = 1.0;
    double physics_loss_weight = 0.1;
    double boundary_loss_weight = 0.1;
    double conservation_loss_weight = 0.05;
    double smoothness_loss_weight = 0.01;
    
    // Collocation points for physics loss
    int num_collocation_points = 1000;
    bool use_random_collocation = true;
    
    // Boundary conditions
    enum class BCType { DIRICHLET, NEUMANN, PERIODIC, MIXED };
    BCType bc_type = BCType::DIRICHLET;
    
    // Conservation laws to enforce
    bool enforce_mass_conservation = false;
    bool enforce_energy_conservation = false;
    bool enforce_momentum_conservation = false;
    
    // Gradient penalty (for well-posedness)
    bool use_gradient_penalty = false;
    double gradient_penalty_weight = 0.01;
    
    // Automatic differentiation
    bool use_autodiff = true;  // Use AD for physics derivatives
    
    // Hard vs soft constraints
    bool use_hard_constraints = false;  // Project onto constraint manifold
    
    // Curriculum learning
    bool use_curriculum = true;
    double physics_weight_start = 0.01;
    double physics_weight_end = 0.5;
    int curriculum_epochs = 50;
};

/**
 * @brief Extended Physics-Informed Loss Functions
 */
class ExtendedPhysicsLoss {
public:
    ExtendedPhysicsLoss(const PINNConfig& config);
    
    /**
     * @brief Compute full PINN loss
     */
    double computeLoss(const Tensor& prediction, const Tensor& target,
                      const Tensor& collocation_points);
    
    /**
     * @brief Individual loss components
     */
    double dataLoss(const Tensor& prediction, const Tensor& target);
    double physicsLoss(const Tensor& prediction, const Tensor& collocation_points);
    double boundaryLoss(const Tensor& prediction, const Tensor& boundary_points,
                       const Tensor& boundary_values);
    double conservationLoss(const Tensor& prediction);
    double smoothnessLoss(const Tensor& prediction);
    
    /**
     * @brief Compute PDE residuals
     */
    Tensor computePoissonResidual(const Tensor& u, const Tensor& f,
                                 double dx, double dy);
    Tensor computeDiffusionResidual(const Tensor& u, const Tensor& u_prev,
                                    double D, double dx, double dy, double dt);
    Tensor computeAdvectionResidual(const Tensor& u, const Tensor& u_prev,
                                    const Tensor& v, double dx, double dy, double dt);
    Tensor computeWaveResidual(const Tensor& u, const Tensor& u_prev,
                               const Tensor& u_pprev, double c,
                               double dx, double dy, double dt);
    Tensor computeNavierStokesResidual(const Tensor& u, const Tensor& v,
                                       const Tensor& p, double nu,
                                       double dx, double dy);
    Tensor computeElastodynamicsResidual(const Tensor& displacement,
                                         const Tensor& displacement_prev,
                                         double rho, double lambda, double mu,
                                         double dx, double dy, double dt);
    Tensor computeDarcyResidual(const Tensor& pressure, const Tensor& permeability,
                               double mu, double dx, double dy);
    Tensor computeReservoirResidual(const Tensor& pressure, const Tensor& saturation,
                                    const Tensor& permeability, const Tensor& porosity,
                                    double dx, double dy, double dt);
    
    /**
     * @brief Compute spatial derivatives
     */
    Tensor gradient(const Tensor& u, double dx, double dy);
    Tensor divergence(const Tensor& u, const Tensor& v, double dx, double dy);
    Tensor laplacian(const Tensor& u, double dx, double dy);
    Tensor curl(const Tensor& u, const Tensor& v, double dx, double dy);
    
    /**
     * @brief Conservation quantities
     */
    double totalMass(const Tensor& rho);
    double totalEnergy(const Tensor& u, const Tensor& v, double rho);
    Tensor momentum(const Tensor& u, const Tensor& v, double rho);
    
private:
    PINNConfig config_;
    
    // Random generator for collocation points
    std::mt19937 rng_;
    
    // Current training state for curriculum
    int current_epoch_ = 0;
    double current_physics_weight_ = 0.01;
    
    void updateCurriculum(int epoch);
};

/**
 * @brief Physics-Informed Fourier Neural Operator
 */
class PhysicsInformedFNO {
public:
    PhysicsInformedFNO(const PINNConfig& config);
    
    /**
     * @brief Forward pass with physics constraints
     */
    Tensor forward(const Tensor& input);
    
    /**
     * @brief Forward pass that returns physics residual
     */
    std::pair<Tensor, Tensor> forwardWithResidual(const Tensor& input);
    
    /**
     * @brief Predict with uncertainty estimate
     */
    void predict(const Tensor& input, Tensor& prediction, Tensor& uncertainty);
    
    /**
     * @brief Train with physics-informed loss
     */
    void train(const FNODataset& data);
    
    /**
     * @brief Train with additional collocation points
     */
    void trainPhysicsInformed(const FNODataset& labeled_data,
                              const Tensor& collocation_points);
    
    /**
     * @brief Enforce hard constraints on output
     */
    Tensor enforceConstraints(const Tensor& prediction);
    
    /**
     * @brief Get physics residual for trained model
     */
    Tensor getPhysicsResidual(const Tensor& prediction);
    
    /**
     * @brief Check if physics constraints are satisfied
     */
    bool checkPhysicsConstraints(const Tensor& prediction, double tolerance = 1e-3);
    
    // Model I/O
    void save(const std::string& path) const;
    void load(const std::string& path);
    
    // Access base FNO
    FNOModel& getModel() { return *base_model_; }
    const FNOModel& getModel() const { return *base_model_; }
    
private:
    PINNConfig config_;
    std::unique_ptr<FNOModel> base_model_;
    std::unique_ptr<ExtendedPhysicsLoss> physics_loss_;
    
    // For hard constraint enforcement
    std::unique_ptr<MLPModel> constraint_projector_;
    
    // Domain information
    double dx_, dy_, dz_;
    
    // Collocation point generator
    Tensor generateCollocationPoints(int n);
};

/**
 * @brief PINN-based Solver that wraps FNO with physics guarantees
 */
class PINNSolver {
public:
    struct SolverConfig {
        PINNConfig pinn_config;
        
        // Solver behavior
        bool verify_physics = true;       // Check physics after solve
        double physics_tolerance = 1e-3;  // Acceptable physics residual
        int max_correction_steps = 5;     // Newton iterations for correction
        
        // Fallback
        bool use_numerical_fallback = true;
        double fallback_threshold = 1e-2;
    };
    
    PINNSolver(const SolverConfig& config);
    
    /**
     * @brief Solve PDE with physics guarantees
     */
    Tensor solve(const Tensor& initial_condition,
                const Tensor& boundary_conditions,
                double t_final, double dt);
    
    /**
     * @brief Single time step
     */
    Tensor step(const Tensor& current_state, double dt);
    
    /**
     * @brief Correct physics violations
     */
    Tensor correctPhysics(const Tensor& prediction, double dt);
    
    /**
     * @brief Get last physics residual
     */
    double getLastPhysicsResidual() const { return last_residual_; }
    
    /**
     * @brief Train solver
     */
    void train(const std::vector<Tensor>& trajectories,
              const std::vector<double>& times);
    
private:
    SolverConfig config_;
    std::unique_ptr<PhysicsInformedFNO> pinn_;
    
    double last_residual_ = 0.0;
    bool last_used_fallback_ = false;
};

// =============================================================================
// 3D FNO Extension
// =============================================================================

/**
 * @brief 3D Fourier Neural Operator
 */
class FNO3D {
public:
    struct FNO3DConfig : FNOConfig {
        int modes_x = 8;
        int modes_y = 8;
        int modes_z = 8;
        int nx = 64;
        int ny = 64;
        int nz = 64;
    };
    
    FNO3D(const FNO3DConfig& config);
    
    /**
     * @brief 3D forward pass
     */
    Tensor forward(const Tensor& input);
    
    /**
     * @brief 3D FFT operations
     */
    ComplexTensor fft3d(const Tensor& x);
    Tensor ifft3d(const ComplexTensor& x);
    
    void train(const std::vector<Tensor>& inputs,
              const std::vector<Tensor>& outputs);
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    FNO3DConfig config_;
    
    // 3D spectral convolution layers
    std::vector<ComplexTensor> spectral_weights_;
    std::vector<Tensor> linear_weights_;
};

// =============================================================================
// Utility Functions
// =============================================================================

namespace FNOUtils {
    // Grid generation
    Tensor createGrid2D(int nx, int ny, double Lx, double Ly);
    Tensor createGrid3D(int nx, int ny, int nz, double Lx, double Ly, double Lz);
    
    // Normalization
    struct Normalizer {
        Tensor mean, std;
        
        void fit(const FNODataset& data);
        Tensor normalize(const Tensor& x) const;
        Tensor denormalize(const Tensor& x) const;
    };
    
    // Interpolation (for variable resolution)
    Tensor interpolate2D(const Tensor& x, int new_nx, int new_ny);
    Tensor interpolate3D(const Tensor& x, int new_nx, int new_ny, int new_nz);
    
    // Model analysis
    void visualizeFilters(const FNOModel& model, const std::string& output_dir);
    void analyzeSpectrum(const Tensor& x, const std::string& output_path);
}

// =============================================================================
// Configuration Parsing
// =============================================================================

/**
 * @brief Parse FNO configuration from config file
 */
FNOConfig parseFNOConfig(const std::map<std::string, std::string>& config_map);

/**
 * @brief Parse solver method from string
 */
SolverMethod parseSolverMethod(const std::string& method_str);

/**
 * @brief Get solver method as string
 */
std::string solverMethodToString(SolverMethod method);

} // namespace ML
} // namespace FSRM

#endif // FSRM_FOURIER_NEURAL_OPERATOR_HPP
