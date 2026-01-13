/**
 * @file NeuralSurrogates.hpp
 * @brief Neural network surrogates for expensive physics computations
 * 
 * Provides ML-accelerated replacements for:
 * - PVT (Pressure-Volume-Temperature) calculations
 * - Flash calculations for compositional simulation
 * - Relative permeability and capillary pressure
 * - Equation of state evaluations
 * - Linear solver preconditioners
 * 
 * Each surrogate can be used standalone or as drop-in replacement
 * for traditional implementations.
 */

#ifndef FSRM_NEURAL_SURROGATES_HPP
#define FSRM_NEURAL_SURROGATES_HPP

#include "FourierNeuralOperator.hpp"
#include <petsc.h>
#include <vector>
#include <array>
#include <memory>
#include <string>
#include <functional>
#include <unordered_map>

namespace FSRM {
namespace ML {

// =============================================================================
// Forward Declarations
// =============================================================================

class NeuralPVTModel;
class NeuralFlashCalculator;
class NeuralRelPermModel;
class NeuralEOSModel;
class NeuralPreconditioner;
class SurrogateTrainer;

// =============================================================================
// Base Neural Surrogate Class
// =============================================================================

/**
 * @brief Base class for all neural surrogates
 */
class NeuralSurrogate {
public:
    virtual ~NeuralSurrogate() = default;
    
    // Model I/O
    virtual void save(const std::string& path) const = 0;
    virtual void load(const std::string& path) = 0;
    virtual void exportONNX(const std::string& path) const = 0;
    
    // Training interface
    virtual void train(const std::vector<std::vector<double>>& inputs,
                      const std::vector<std::vector<double>>& outputs,
                      int epochs = 100, double learning_rate = 1e-3) = 0;
    
    // Validation
    virtual double validate(const std::vector<std::vector<double>>& inputs,
                           const std::vector<std::vector<double>>& outputs) = 0;
    
    // Status
    virtual bool isTrained() const = 0;
    virtual size_t numParameters() const = 0;
    virtual std::string getModelInfo() const = 0;
    
    // GPU control
    virtual void toGPU() = 0;
    virtual void toCPU() = 0;
    virtual bool isOnGPU() const = 0;
    
    // Uncertainty estimation (if available)
    virtual bool hasUncertainty() const { return false; }
    virtual double getLastUncertainty() const { return 0.0; }
    
protected:
    bool trained_ = false;
    bool on_gpu_ = false;
    double last_uncertainty_ = 0.0;
};

// =============================================================================
// Multi-Layer Perceptron Base
// =============================================================================

/**
 * @brief Configuration for MLP architecture
 */
struct MLPConfig {
    std::vector<int> layer_sizes;      // Including input and output
    std::string activation = "gelu";    // relu, gelu, tanh, silu, swish
    bool use_batch_norm = true;
    bool use_dropout = false;
    double dropout_rate = 0.1;
    bool use_residual = true;
    bool use_layer_norm = false;
    
    // Training
    double learning_rate = 1e-3;
    double weight_decay = 1e-5;
    int batch_size = 64;
    int epochs = 100;
    
    // Initialization
    std::string weight_init = "kaiming";  // xavier, kaiming, orthogonal
};

/**
 * @brief Dense layer for MLP
 */
class DenseLayer {
public:
    DenseLayer(int in_features, int out_features, bool bias = true);
    
    Tensor forward(const Tensor& x);
    Tensor backward(const Tensor& grad_output);
    
    Tensor weight;
    Tensor bias_vec;
    Tensor grad_weight;
    Tensor grad_bias;
    
    bool has_bias;
    Tensor last_input;
    
    void initWeights(const std::string& method);
};

/**
 * @brief Batch normalization layer
 */
class BatchNorm1d {
public:
    BatchNorm1d(int num_features, double eps = 1e-5, double momentum = 0.1);
    
    Tensor forward(const Tensor& x, bool training = true);
    Tensor backward(const Tensor& grad_output);
    
    Tensor gamma;
    Tensor beta;
    Tensor running_mean;
    Tensor running_var;
    
    Tensor grad_gamma;
    Tensor grad_beta;
    
private:
    int num_features_;
    double eps_;
    double momentum_;
    Tensor cached_std_;
    Tensor cached_normalized_;
};

/**
 * @brief Dropout layer
 */
class Dropout {
public:
    Dropout(double p = 0.5);
    
    Tensor forward(const Tensor& x, bool training = true);
    Tensor backward(const Tensor& grad_output);
    
private:
    double p_;
    Tensor mask_;
};

/**
 * @brief Multi-layer perceptron model
 */
class MLPModel {
public:
    MLPModel(const MLPConfig& config);
    
    Tensor forward(const Tensor& x);
    Tensor backward(const Tensor& grad_output);
    Tensor predict(const Tensor& x);
    
    // Loss
    double computeLoss(const Tensor& prediction, const Tensor& target);
    Tensor computeLossGrad(const Tensor& prediction, const Tensor& target);
    
    // Parameters
    std::vector<Tensor*> parameters();
    std::vector<Tensor*> gradients();
    void zeroGrad();
    
    // Model I/O
    void save(const std::string& path) const;
    void load(const std::string& path);
    
    // Info
    size_t numParameters() const;
    void summary() const;
    
    // Training mode
    void train() { training_ = true; }
    void eval() { training_ = false; }
    
private:
    MLPConfig config_;
    std::vector<std::unique_ptr<DenseLayer>> dense_layers_;
    std::vector<std::unique_ptr<BatchNorm1d>> bn_layers_;
    std::vector<std::unique_ptr<Dropout>> dropout_layers_;
    
    Activation::ActivationFn activation_fn_;
    Activation::ActivationBackwardFn activation_backward_;
    
    bool training_ = true;
    std::vector<Tensor> layer_outputs_;  // For backward pass
};

// =============================================================================
// Neural PVT Model
// =============================================================================

/**
 * @brief PVT property types
 */
enum class PVTProperty {
    OIL_FVF,           // Formation volume factor Bo
    OIL_VISCOSITY,     // Oil viscosity μo
    OIL_DENSITY,       // Oil density ρo
    OIL_COMPRESSIBILITY, // Oil compressibility co
    GAS_FVF,           // Formation volume factor Bg
    GAS_VISCOSITY,     // Gas viscosity μg
    GAS_DENSITY,       // Gas density ρg
    GAS_Z_FACTOR,      // Gas compressibility factor Z
    WATER_FVF,         // Formation volume factor Bw
    WATER_VISCOSITY,   // Water viscosity μw
    WATER_DENSITY,     // Water density ρw
    SOLUTION_GOR,      // Solution gas-oil ratio Rs
    VAPORIZED_OGR,     // Vaporized oil-gas ratio Rv
    BUBBLE_POINT,      // Bubble point pressure Pb
    DEW_POINT          // Dew point pressure Pd
};

/**
 * @brief Configuration for neural PVT model
 */
struct NeuralPVTConfig {
    // Architecture
    std::vector<int> hidden_layers = {128, 256, 256, 128};
    std::string activation = "gelu";
    bool use_batch_norm = true;
    
    // Input normalization bounds
    double p_min = 1e5;      // Pa
    double p_max = 1e8;      // Pa
    double T_min = 273.15;   // K
    double T_max = 473.15;   // K
    double Rs_min = 0.0;     // sm3/sm3
    double Rs_max = 500.0;   // sm3/sm3
    
    // Output properties to predict
    std::vector<PVTProperty> properties = {
        PVTProperty::OIL_FVF,
        PVTProperty::OIL_VISCOSITY,
        PVTProperty::GAS_FVF,
        PVTProperty::GAS_VISCOSITY,
        PVTProperty::SOLUTION_GOR
    };
    
    // Training
    double learning_rate = 1e-3;
    int batch_size = 256;
    int epochs = 200;
    
    // Constraints
    bool enforce_monotonicity = true;
    bool enforce_positivity = true;
};

/**
 * @brief Neural network model for PVT property prediction
 * 
 * Replaces expensive PVT correlations (Standing, Vazquez-Beggs, etc.)
 * with fast neural network inference.
 * 
 * Input: (pressure, temperature, [composition])
 * Output: All requested PVT properties
 */
class NeuralPVTModel : public NeuralSurrogate {
public:
    NeuralPVTModel(const NeuralPVTConfig& config);
    ~NeuralPVTModel() override = default;
    
    // Main prediction interface
    void predict(double pressure, double temperature, double Rs,
                std::vector<double>& properties);
    
    // Batch prediction (more efficient)
    void predictBatch(const std::vector<double>& pressures,
                     const std::vector<double>& temperatures,
                     const std::vector<double>& Rs_values,
                     std::vector<std::vector<double>>& properties);
    
    // Individual property access
    double getOilFVF(double p, double T, double Rs);
    double getOilViscosity(double p, double T, double Rs);
    double getGasFVF(double p, double T);
    double getGasViscosity(double p, double T);
    double getSolutionGOR(double p, double T);
    double getBubblePoint(double T, double Rs);
    
    // With uncertainty
    void predictWithUncertainty(double pressure, double temperature, double Rs,
                               std::vector<double>& properties,
                               std::vector<double>& uncertainties);
    
    // Training interface
    void train(const std::vector<std::vector<double>>& inputs,
              const std::vector<std::vector<double>>& outputs,
              int epochs, double learning_rate) override;
    
    // Generate training data from correlations
    void generateTrainingData(
        std::function<void(double, double, double, std::vector<double>&)> pvt_func,
        int num_samples);
    
    // Model I/O
    void save(const std::string& path) const override;
    void load(const std::string& path) override;
    void exportONNX(const std::string& path) const override;
    
    // Validation
    double validate(const std::vector<std::vector<double>>& inputs,
                   const std::vector<std::vector<double>>& outputs) override;
    
    // Status
    bool isTrained() const override { return trained_; }
    size_t numParameters() const override;
    std::string getModelInfo() const override;
    
    // GPU
    void toGPU() override;
    void toCPU() override;
    bool isOnGPU() const override { return on_gpu_; }
    
private:
    NeuralPVTConfig config_;
    std::unique_ptr<MLPModel> model_;
    
    // Normalization parameters
    Tensor input_mean_, input_std_;
    Tensor output_mean_, output_std_;
    
    // Helper methods
    Tensor normalizeInput(double p, double T, double Rs);
    std::vector<double> denormalizeOutput(const Tensor& output);
    void enforceConstraints(std::vector<double>& properties);
};

// =============================================================================
// Neural Flash Calculator
// =============================================================================

/**
 * @brief Configuration for neural flash calculator
 */
struct NeuralFlashConfig {
    int num_components = 5;
    std::vector<int> hidden_layers = {256, 512, 512, 256};
    std::string activation = "gelu";
    
    // Component properties for normalization
    std::vector<double> Tc;  // Critical temperatures
    std::vector<double> Pc;  // Critical pressures
    std::vector<double> omega;  // Acentric factors
    
    // Training
    double learning_rate = 1e-3;
    int batch_size = 128;
    int epochs = 500;
    
    // Physics constraints
    bool enforce_mass_balance = true;
    bool enforce_fugacity_equality = false;  // Can be expensive
};

/**
 * @brief Neural network flash calculator
 * 
 * Replaces iterative Rachford-Rice flash calculations with neural prediction.
 * 
 * Input: pressure, temperature, overall composition z[]
 * Output: vapor fraction L, liquid composition x[], vapor composition y[]
 */
class NeuralFlashCalculator : public NeuralSurrogate {
public:
    NeuralFlashCalculator(const NeuralFlashConfig& config);
    ~NeuralFlashCalculator() override = default;
    
    /**
     * @brief Perform flash calculation
     * 
     * @param z Overall composition (mole fractions, sum = 1)
     * @param P Pressure [Pa]
     * @param T Temperature [K]
     * @param L Output: liquid mole fraction
     * @param x Output: liquid phase composition
     * @param y Output: vapor phase composition
     */
    void flash(const std::vector<double>& z, double P, double T,
              double& L, std::vector<double>& x, std::vector<double>& y);
    
    /**
     * @brief Flash with stability check
     */
    void flashWithStability(const std::vector<double>& z, double P, double T,
                           double& L, std::vector<double>& x, std::vector<double>& y,
                           bool& is_two_phase);
    
    /**
     * @brief Batch flash for multiple cells (GPU-accelerated)
     */
    void flashBatch(const std::vector<std::vector<double>>& z_batch,
                   const std::vector<double>& P_batch,
                   const std::vector<double>& T_batch,
                   std::vector<double>& L_batch,
                   std::vector<std::vector<double>>& x_batch,
                   std::vector<std::vector<double>>& y_batch);
    
    /**
     * @brief Get fugacity coefficients (for phase equilibrium check)
     */
    void getFugacityCoefficients(const std::vector<double>& composition,
                                 double P, double T,
                                 std::vector<double>& phi);
    
    // Training
    void train(const std::vector<std::vector<double>>& inputs,
              const std::vector<std::vector<double>>& outputs,
              int epochs, double learning_rate) override;
    
    /**
     * @brief Generate training data from EOS
     */
    void generateTrainingData(
        std::function<void(const std::vector<double>&, double, double,
                          double&, std::vector<double>&, std::vector<double>&)> eos_flash,
        int num_samples);
    
    // Model I/O
    void save(const std::string& path) const override;
    void load(const std::string& path) override;
    void exportONNX(const std::string& path) const override;
    
    double validate(const std::vector<std::vector<double>>& inputs,
                   const std::vector<std::vector<double>>& outputs) override;
    
    bool isTrained() const override { return trained_; }
    size_t numParameters() const override;
    std::string getModelInfo() const override;
    
    void toGPU() override;
    void toCPU() override;
    bool isOnGPU() const override { return on_gpu_; }
    
private:
    NeuralFlashConfig config_;
    std::unique_ptr<MLPModel> phase_split_model_;   // Predicts L
    std::unique_ptr<MLPModel> composition_model_;   // Predicts x, y
    std::unique_ptr<MLPModel> stability_model_;     // Two-phase check
    
    // Normalization
    Tensor z_mean_, z_std_;
    double P_mean_, P_std_;
    double T_mean_, T_std_;
    
    void enforceConstraints(double& L, std::vector<double>& x, 
                           std::vector<double>& y, const std::vector<double>& z);
};

// =============================================================================
// Neural Relative Permeability Model
// =============================================================================

/**
 * @brief Configuration for neural relative permeability
 */
struct NeuralRelPermConfig {
    int num_phases = 3;  // water, oil, gas
    std::vector<int> hidden_layers = {64, 128, 128, 64};
    std::string activation = "softplus";  // Ensures positivity
    
    // Saturation bounds
    double Swc = 0.2;   // Connate water
    double Sor = 0.2;   // Residual oil
    double Sgc = 0.05;  // Critical gas
    
    // Training
    double learning_rate = 1e-3;
    int epochs = 200;
    
    // Constraints
    bool enforce_endpoint_values = true;
    bool enforce_monotonicity = true;
};

/**
 * @brief Neural relative permeability and capillary pressure
 */
class NeuralRelPermModel : public NeuralSurrogate {
public:
    NeuralRelPermModel(const NeuralRelPermConfig& config);
    ~NeuralRelPermModel() override = default;
    
    /**
     * @brief Get relative permeabilities
     * 
     * @param Sw Water saturation
     * @param Sg Gas saturation (So = 1 - Sw - Sg)
     * @param krw Output: water relative permeability
     * @param kro Output: oil relative permeability
     * @param krg Output: gas relative permeability
     */
    void getRelPerm(double Sw, double Sg,
                   double& krw, double& kro, double& krg);
    
    /**
     * @brief Get relative permeabilities with derivatives
     */
    void getRelPermWithDerivatives(double Sw, double Sg,
                                   double& krw, double& kro, double& krg,
                                   double& dkrw_dSw, double& dkro_dSw, double& dkrg_dSw,
                                   double& dkrw_dSg, double& dkro_dSg, double& dkrg_dSg);
    
    /**
     * @brief Batch evaluation for all cells
     */
    void getRelPermBatch(const std::vector<double>& Sw,
                        const std::vector<double>& Sg,
                        std::vector<double>& krw,
                        std::vector<double>& kro,
                        std::vector<double>& krg);
    
    /**
     * @brief Get capillary pressures
     */
    void getCapillaryPressure(double Sw, double Sg,
                             double& Pcow, double& Pcgo);
    
    // Training
    void train(const std::vector<std::vector<double>>& inputs,
              const std::vector<std::vector<double>>& outputs,
              int epochs, double learning_rate) override;
    
    // Generate from Corey/Brooks-Corey
    void generateFromCorey(double nw, double no, double ng,
                          double krw_max, double kro_max, double krg_max,
                          int num_samples);
    
    // Model I/O
    void save(const std::string& path) const override;
    void load(const std::string& path) override;
    void exportONNX(const std::string& path) const override;
    
    double validate(const std::vector<std::vector<double>>& inputs,
                   const std::vector<std::vector<double>>& outputs) override;
    
    bool isTrained() const override { return trained_; }
    size_t numParameters() const override;
    std::string getModelInfo() const override;
    
    void toGPU() override;
    void toCPU() override;
    bool isOnGPU() const override { return on_gpu_; }
    
private:
    NeuralRelPermConfig config_;
    std::unique_ptr<MLPModel> krw_model_;
    std::unique_ptr<MLPModel> kro_model_;
    std::unique_ptr<MLPModel> krg_model_;
    std::unique_ptr<MLPModel> pcow_model_;
    std::unique_ptr<MLPModel> pcgo_model_;
};

// =============================================================================
// Neural Equation of State
// =============================================================================

/**
 * @brief Neural equation of state model
 * 
 * Learns the P-V-T relationship and derivative properties.
 */
class NeuralEOSModel : public NeuralSurrogate {
public:
    struct EOSConfig {
        int num_components = 5;
        std::vector<int> hidden_layers = {256, 512, 512, 256};
        std::string eos_type = "PR";  // Base EOS for training data: PR, SRK, vdW
        
        // Component critical properties
        std::vector<double> Tc;
        std::vector<double> Pc;
        std::vector<double> omega;
        std::vector<double> Vc;
        
        // Binary interaction parameters (flattened upper triangle)
        std::vector<double> kij;
    };
    
    NeuralEOSModel(const EOSConfig& config);
    ~NeuralEOSModel() override = default;
    
    /**
     * @brief Compute compressibility factor Z
     */
    double computeZ(const std::vector<double>& composition, double P, double T,
                   bool vapor_phase = true);
    
    /**
     * @brief Compute molar volume
     */
    double computeMolarVolume(const std::vector<double>& composition, 
                             double P, double T, bool vapor_phase = true);
    
    /**
     * @brief Compute fugacity coefficients
     */
    void computeFugacityCoefficients(const std::vector<double>& composition,
                                     double P, double T, bool vapor_phase,
                                     std::vector<double>& phi);
    
    /**
     * @brief Compute all thermodynamic properties
     */
    struct ThermoProperties {
        double Z;           // Compressibility factor
        double V;           // Molar volume
        double H;           // Enthalpy
        double S;           // Entropy
        double Cp;          // Heat capacity at constant P
        double Cv;          // Heat capacity at constant V
        std::vector<double> phi;  // Fugacity coefficients
    };
    
    ThermoProperties computeProperties(const std::vector<double>& composition,
                                       double P, double T, bool vapor_phase);
    
    // Training
    void train(const std::vector<std::vector<double>>& inputs,
              const std::vector<std::vector<double>>& outputs,
              int epochs, double learning_rate) override;
    
    void save(const std::string& path) const override;
    void load(const std::string& path) override;
    void exportONNX(const std::string& path) const override;
    
    double validate(const std::vector<std::vector<double>>& inputs,
                   const std::vector<std::vector<double>>& outputs) override;
    
    bool isTrained() const override { return trained_; }
    size_t numParameters() const override;
    std::string getModelInfo() const override;
    
    void toGPU() override;
    void toCPU() override;
    bool isOnGPU() const override { return on_gpu_; }
    
private:
    EOSConfig config_;
    std::unique_ptr<MLPModel> z_model_;
    std::unique_ptr<MLPModel> phi_model_;
    std::unique_ptr<MLPModel> thermo_model_;
};

// =============================================================================
// Neural Preconditioner
// =============================================================================

/**
 * @brief Neural network preconditioner for iterative solvers
 * 
 * Learns to approximate the inverse of the Jacobian for faster
 * Krylov solver convergence.
 */
class NeuralPreconditioner : public NeuralSurrogate {
public:
    struct PrecondConfig {
        // Architecture
        std::vector<int> hidden_layers = {512, 1024, 1024, 512};
        std::string activation = "gelu";
        
        // Input/output dimension (set from matrix)
        int n_dof = 0;
        
        // Training
        double learning_rate = 1e-4;
        int batch_size = 32;
        int epochs = 100;
        
        // Preconditioning strategy
        enum class Strategy {
            DIRECT,          // Learn M ≈ A^{-1}
            RESIDUAL,        // Learn correction to ILU
            DEFLATION,       // Learn low-rank deflation
            COARSE_GRID      // Learn coarse grid correction
        };
        Strategy strategy = Strategy::DIRECT;
        
        // Reference preconditioner for residual strategy
        std::string reference_precond = "ILU0";
        
        // Spectral information
        bool use_eigenspace = false;
        int num_eigenvectors = 10;
    };
    
    NeuralPreconditioner(const PrecondConfig& config);
    ~NeuralPreconditioner() override;
    
    /**
     * @brief Apply preconditioner: y = M^{-1} * x
     */
    PetscErrorCode apply(Vec x, Vec y);
    
    /**
     * @brief Apply preconditioner (Tensor interface)
     */
    Tensor apply(const Tensor& residual);
    
    /**
     * @brief Setup with matrix (extracts structure)
     */
    PetscErrorCode setup(Mat A);
    
    /**
     * @brief Update for new matrix (same sparsity pattern)
     */
    PetscErrorCode update(Mat A);
    
    /**
     * @brief Train on matrix-vector pairs
     * 
     * @param A_samples Vector of matrices (different parameter values)
     * @param b_samples Corresponding right-hand sides
     * @param x_samples Corresponding solutions
     */
    void trainFromSamples(const std::vector<Mat>& A_samples,
                         const std::vector<Vec>& b_samples,
                         const std::vector<Vec>& x_samples);
    
    /**
     * @brief Get iteration reduction factor
     */
    double getIterationReduction() const { return iteration_reduction_; }
    
    // Training
    void train(const std::vector<std::vector<double>>& inputs,
              const std::vector<std::vector<double>>& outputs,
              int epochs, double learning_rate) override;
    
    void save(const std::string& path) const override;
    void load(const std::string& path) override;
    void exportONNX(const std::string& path) const override;
    
    double validate(const std::vector<std::vector<double>>& inputs,
                   const std::vector<std::vector<double>>& outputs) override;
    
    bool isTrained() const override { return trained_; }
    size_t numParameters() const override;
    std::string getModelInfo() const override;
    
    void toGPU() override;
    void toCPU() override;
    bool isOnGPU() const override { return on_gpu_; }
    
private:
    PrecondConfig config_;
    std::unique_ptr<MLPModel> precond_model_;
    
    // Matrix information
    int n_rows_, n_cols_, nnz_;
    std::vector<PetscInt> row_ptr_, col_idx_;
    
    // Performance tracking
    double iteration_reduction_ = 1.0;
    int total_applies_ = 0;
    double total_time_ = 0.0;
    
    // Reference preconditioner for residual strategy
    PC reference_pc_;
    
    // Low-rank deflation vectors
    std::vector<Tensor> deflation_vectors_;
};

// =============================================================================
// Surrogate Trainer Utilities
// =============================================================================

/**
 * @brief Training data generator and utilities
 */
class SurrogateTrainer {
public:
    SurrogateTrainer();
    
    /**
     * @brief Generate Latin Hypercube samples
     */
    static std::vector<std::vector<double>> latinHypercube(
        const std::vector<double>& lower_bounds,
        const std::vector<double>& upper_bounds,
        int num_samples);
    
    /**
     * @brief Generate Sobol sequence samples
     */
    static std::vector<std::vector<double>> sobolSequence(
        const std::vector<double>& lower_bounds,
        const std::vector<double>& upper_bounds,
        int num_samples);
    
    /**
     * @brief Cross-validation
     */
    template<typename Model>
    static std::vector<double> crossValidate(
        Model& model,
        const std::vector<std::vector<double>>& inputs,
        const std::vector<std::vector<double>>& outputs,
        int k_folds = 5);
    
    /**
     * @brief Learning rate finder
     */
    template<typename Model>
    static double findOptimalLR(
        Model& model,
        const std::vector<std::vector<double>>& inputs,
        const std::vector<std::vector<double>>& outputs,
        double lr_min = 1e-7, double lr_max = 1.0);
    
    /**
     * @brief Data augmentation with noise
     */
    static void augmentWithNoise(
        std::vector<std::vector<double>>& inputs,
        std::vector<std::vector<double>>& outputs,
        double noise_level = 0.01,
        int augmentation_factor = 2);
};

// =============================================================================
// Configuration Parsing
// =============================================================================

NeuralPVTConfig parseNeuralPVTConfig(const std::map<std::string, std::string>& config);
NeuralFlashConfig parseNeuralFlashConfig(const std::map<std::string, std::string>& config);
NeuralRelPermConfig parseNeuralRelPermConfig(const std::map<std::string, std::string>& config);

} // namespace ML
} // namespace FSRM

#endif // FSRM_NEURAL_SURROGATES_HPP
