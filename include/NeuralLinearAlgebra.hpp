/**
 * @file NeuralLinearAlgebra.hpp
 * @brief AI-Accelerated Linear Algebra Operations
 * 
 * Provides neural network enhanced linear algebra components:
 * 
 * - Neural Preconditioners: Learning optimal preconditioning operators
 * - Neural Coarse-Grid Correction: ML-enhanced multigrid
 * - Learned Matrix Ordering: Optimal ordering for sparse solvers
 * - Neural AMG Coarsening: Learning algebraic multigrid hierarchies
 * - Neural Jacobian Approximation: Fast Jacobian estimation
 * 
 * These components integrate with PETSc's solver infrastructure
 * for seamless use in existing simulations.
 */

#ifndef FSRM_NEURAL_LINEAR_ALGEBRA_HPP
#define FSRM_NEURAL_LINEAR_ALGEBRA_HPP

#include "FourierNeuralOperator.hpp"
#include "NeuralSurrogates.hpp"
#include <petsc.h>
#include <petscksp.h>
#include <petscmat.h>
#include <vector>
#include <memory>
#include <map>

// Forward declarations for SLEPc types (optional dependency)
typedef struct _p_EPS* EPS;

namespace FSRM {
namespace ML {

// =============================================================================
// Neural Coarse-Grid Correction
// =============================================================================

/**
 * @brief Configuration for neural coarse-grid operator
 */
struct NeuralCoarseGridConfig {
    // Coarse grid dimensions (learned mapping)
    int fine_dim = 1000;
    int coarse_dim = 100;
    
    // Network architecture
    std::vector<int> restriction_layers = {256, 256};
    std::vector<int> coarse_solve_layers = {256, 512, 256};
    std::vector<int> prolongation_layers = {256, 256};
    std::string activation = "gelu";
    
    // Training
    double learning_rate = 1e-4;
    int epochs = 200;
    int batch_size = 32;
    
    // Hybrid mode
    bool use_geometric_restriction = false;  // Use geometric + neural
    bool use_exact_coarse_solve = false;     // Use direct solve at coarse
};

/**
 * @brief Neural network-based coarse-grid correction for multigrid
 * 
 * Learns optimal restriction, coarse-grid solve, and prolongation
 * operators to accelerate multigrid convergence.
 */
class NeuralCoarseGridCorrection {
public:
    NeuralCoarseGridCorrection(const NeuralCoarseGridConfig& config);
    
    /**
     * @brief Apply coarse-grid correction
     * 
     * @param residual Fine-grid residual
     * @return Coarse-grid correction prolongated to fine grid
     */
    Tensor apply(const Tensor& residual);
    
    /**
     * @brief PETSc interface
     */
    PetscErrorCode apply(Vec residual, Vec correction);
    
    /**
     * @brief Setup from fine-level operator
     */
    PetscErrorCode setup(Mat A_fine);
    
    /**
     * @brief Train from fine-grid solutions
     */
    void train(const std::vector<Tensor>& residuals,
              const std::vector<Tensor>& corrections);
    
    /**
     * @brief Train with geometric information
     */
    void trainWithGeometry(const std::vector<Tensor>& residuals,
                          const std::vector<Tensor>& corrections,
                          DM dm_fine, DM dm_coarse);
    
    // Access learned operators
    Tensor restrict(const Tensor& fine);
    Tensor prolongate(const Tensor& coarse);
    Tensor coarseSolve(const Tensor& coarse_residual);
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    NeuralCoarseGridConfig config_;
    
    std::unique_ptr<MLPModel> restriction_net_;
    std::unique_ptr<MLPModel> coarse_solve_net_;
    std::unique_ptr<MLPModel> prolongation_net_;
    
    // Geometric operators (optional hybrid mode)
    Tensor geometric_restriction_;
    Tensor geometric_prolongation_;
    
    // Statistics
    double avg_iteration_reduction_ = 1.0;
};

// =============================================================================
// Learned Matrix Ordering
// =============================================================================

/**
 * @brief Neural network for optimal matrix ordering
 * 
 * Learns to predict optimal row/column ordering to minimize
 * fill-in during factorization or improve convergence.
 */
class NeuralMatrixOrdering {
public:
    struct OrderingConfig {
        enum class OrderingType {
            BANDWIDTH_REDUCTION,    // Minimize bandwidth
            FILL_REDUCTION,         // Minimize factorization fill
            CONVERGENCE_OPTIMAL     // Optimal for iterative solver
        };
        OrderingType type = OrderingType::BANDWIDTH_REDUCTION;
        
        // Graph features
        int max_nodes = 10000;
        int num_gnn_layers = 4;
        int hidden_dim = 64;
        
        // Training
        double learning_rate = 1e-3;
        int epochs = 100;
    };
    
    NeuralMatrixOrdering(const OrderingConfig& config);
    
    /**
     * @brief Compute optimal ordering for matrix
     */
    PetscErrorCode computeOrdering(Mat A, IS* ordering);
    
    /**
     * @brief Train on matrices with known optimal orderings
     */
    void train(const std::vector<Mat>& matrices,
              const std::vector<IS>& optimal_orderings);
    
    /**
     * @brief Evaluate ordering quality
     */
    double evaluateOrdering(Mat A, IS ordering);
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    OrderingConfig config_;
    
    // GNN for graph representation
    std::unique_ptr<MLPModel> node_encoder_;
    std::unique_ptr<MLPModel> score_predictor_;
    
    // Convert matrix to graph
    void matrixToGraph(Mat A, Tensor& node_features, Tensor& edge_index);
};

// =============================================================================
// Neural AMG Coarsening
// =============================================================================

/**
 * @brief Learn optimal algebraic multigrid coarsening
 */
class NeuralAMGCoarsening {
public:
    struct AMGConfig {
        // Coarsening strategy
        enum class Strategy {
            STRENGTH_BASED,     // Traditional strength-based
            NEURAL_STRENGTH,    // Learned strength threshold
            NEURAL_CLUSTERING,  // Learned node clustering
            GRAPH_PARTITIONING  // GNN-based partitioning
        };
        Strategy strategy = Strategy::NEURAL_STRENGTH;
        
        // Network
        int hidden_dim = 64;
        int num_layers = 3;
        
        // Training
        double learning_rate = 1e-3;
        int epochs = 100;
    };
    
    NeuralAMGCoarsening(const AMGConfig& config);
    
    /**
     * @brief Compute coarse/fine splitting
     */
    PetscErrorCode computeCFSplitting(Mat A, IS* fine_nodes, IS* coarse_nodes);
    
    /**
     * @brief Compute interpolation operator
     */
    PetscErrorCode computeInterpolation(Mat A, IS coarse_nodes, Mat* P);
    
    /**
     * @brief Train on matrix hierarchy
     */
    void train(const std::vector<Mat>& matrices,
              const std::vector<std::pair<IS, IS>>& cf_splittings);
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    AMGConfig config_;
    
    // Strength prediction network
    std::unique_ptr<MLPModel> strength_predictor_;
    
    // C/F node classifier
    std::unique_ptr<MLPModel> cf_classifier_;
    
    // Interpolation weight predictor
    std::unique_ptr<MLPModel> weight_predictor_;
};

// =============================================================================
// Neural Jacobian Approximation
// =============================================================================

/**
 * @brief Learn to approximate Jacobian matrices efficiently
 */
class NeuralJacobianApproximator {
public:
    struct JacobianConfig {
        // Approximation method
        enum class Method {
            FULL_DENSE,         // Learn dense Jacobian
            LOW_RANK,           // Learn low-rank factors
            SPARSE_PATTERN,     // Learn sparse pattern + values
            BLOCK_DIAGONAL      // Learn block structure
        };
        Method method = Method::LOW_RANK;
        
        // Low-rank settings
        int rank = 50;
        
        // Network
        std::vector<int> hidden_layers = {256, 256};
        std::string activation = "gelu";
        
        // Training
        double learning_rate = 1e-3;
        int epochs = 200;
    };
    
    NeuralJacobianApproximator(const JacobianConfig& config);
    
    /**
     * @brief Approximate Jacobian at given state
     */
    Tensor approximate(const Tensor& state);
    
    /**
     * @brief Apply Jacobian to vector (J * v)
     */
    Tensor apply(const Tensor& state, const Tensor& vector);
    
    /**
     * @brief Apply transpose Jacobian (J^T * v)
     */
    Tensor applyTranspose(const Tensor& state, const Tensor& vector);
    
    /**
     * @brief Train from state-Jacobian pairs
     */
    void train(const std::vector<Tensor>& states,
              const std::vector<Tensor>& jacobians);
    
    /**
     * @brief Train from finite differences
     */
    void trainFromFunction(
        std::function<Tensor(const Tensor&)> residual_func,
        const std::vector<Tensor>& sample_states,
        double epsilon = 1e-6);
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    JacobianConfig config_;
    
    // For full/low-rank method
    std::unique_ptr<MLPModel> jacobian_net_;
    
    // For low-rank: J ≈ U * V^T
    std::unique_ptr<MLPModel> U_net_;
    std::unique_ptr<MLPModel> V_net_;
    
    // For sparse method
    std::unique_ptr<MLPModel> sparsity_net_;
    std::unique_ptr<MLPModel> value_net_;
};

// =============================================================================
// Neural Spectral Solver
// =============================================================================

/**
 * @brief Learn eigenvalue/eigenvector computation
 */
class NeuralSpectralSolver {
public:
    struct SpectralConfig {
        int num_eigenvalues = 10;  // Number of eigenvalues to compute
        
        enum class Target {
            SMALLEST,           // Smallest magnitude
            LARGEST,            // Largest magnitude
            SMALLEST_REAL,      // Smallest real part
            LARGEST_REAL        // Largest real part
        };
        Target target = Target::SMALLEST;
        
        // Network
        std::vector<int> hidden_layers = {256, 512, 256};
        
        // Training
        double learning_rate = 1e-3;
        int epochs = 100;
    };
    
    NeuralSpectralSolver(const SpectralConfig& config);
    
    /**
     * @brief Predict eigenvalues
     */
    Tensor predictEigenvalues(Mat A);
    
    /**
     * @brief Predict eigenvectors
     */
    Tensor predictEigenvectors(Mat A);
    
    /**
     * @brief Provide initial guess for iterative eigensolver
     */
    PetscErrorCode initializeEPS(EPS eps, Mat A);
    
    void train(const std::vector<Mat>& matrices,
              const std::vector<Tensor>& eigenvalues,
              const std::vector<Tensor>& eigenvectors);
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    SpectralConfig config_;
    
    std::unique_ptr<MLPModel> eigenvalue_net_;
    std::unique_ptr<MLPModel> eigenvector_net_;
    
    // Matrix features
    Tensor extractFeatures(Mat A);
};

// =============================================================================
// Neural Deflation Spaces
// =============================================================================

/**
 * @brief Learn optimal deflation vectors for Krylov methods
 */
class NeuralDeflation {
public:
    struct DeflationConfig {
        int num_deflation_vectors = 10;
        
        enum class VectorType {
            EIGENVECTORS,       // Approximate eigenvectors
            RECYCLED,           // Recycled Krylov vectors
            PROBLEM_ADAPTED     // Learned problem-specific
        };
        VectorType type = VectorType::PROBLEM_ADAPTED;
        
        // Network
        std::vector<int> hidden_layers = {256, 256};
        
        // Training
        int epochs = 100;
    };
    
    NeuralDeflation(const DeflationConfig& config);
    
    /**
     * @brief Compute deflation subspace for matrix
     */
    PetscErrorCode computeDeflationSpace(Mat A, Mat* W);
    
    /**
     * @brief Apply deflated preconditioning
     */
    PetscErrorCode applyDeflatedPC(PC pc_base, Mat W, Vec x, Vec y);
    
    /**
     * @brief Train from solve sequences
     */
    void train(const std::vector<Mat>& matrices,
              const std::vector<std::vector<Vec>>& krylov_spaces);
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    DeflationConfig config_;
    
    std::unique_ptr<MLPModel> deflation_net_;
    
    // Cached deflation space
    Mat deflation_space_ = nullptr;
};

// =============================================================================
// Unified Neural Linear Solver
// =============================================================================

/**
 * @brief Complete neural-enhanced linear solver
 */
class NeuralLinearSolver {
public:
    struct SolverConfig {
        // Solver selection
        enum class SolverType {
            NEURAL_PRECOND_ONLY,    // Use neural preconditioner
            NEURAL_COARSE_GRID,     // Neural coarse correction
            NEURAL_AMG,             // Full neural AMG
            HYBRID,                 // Neural + classical
            PURE_NEURAL             // Direct neural solve
        };
        SolverType solver_type = SolverType::HYBRID;
        
        // Component configs
        NeuralCoarseGridConfig coarse_config;
        NeuralAMGCoarsening::AMGConfig amg_config;
        NeuralDeflation::DeflationConfig deflation_config;
        
        // Fallback
        std::string fallback_solver = "gmres";
        std::string fallback_precond = "ilu";
        
        // Convergence
        double rtol = 1e-8;
        int max_iterations = 1000;
        
        // Training
        bool online_learning = false;  // Adapt during simulation
    };
    
    NeuralLinearSolver(const SolverConfig& config);
    
    /**
     * @brief Solve linear system Ax = b
     */
    PetscErrorCode solve(Mat A, Vec b, Vec x);
    
    /**
     * @brief Setup for matrix (can be called for pattern only)
     */
    PetscErrorCode setup(Mat A);
    
    /**
     * @brief Update for new matrix values (same pattern)
     */
    PetscErrorCode update(Mat A);
    
    /**
     * @brief Train all components
     */
    void train(const std::vector<Mat>& matrices,
              const std::vector<Vec>& rhs_vectors,
              const std::vector<Vec>& solutions);
    
    /**
     * @brief Online adaptation from recent solves
     */
    void adaptOnline();
    
    // Statistics
    struct SolverStats {
        int total_iterations = 0;
        double total_time = 0.0;
        double neural_time = 0.0;
        double classical_time = 0.0;
        int num_solves = 0;
        double avg_iteration_reduction = 1.0;
    };
    
    const SolverStats& getStats() const { return stats_; }
    
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    SolverConfig config_;
    SolverStats stats_;
    
    // Components
    std::unique_ptr<NeuralPreconditioner> neural_precond_;
    std::unique_ptr<NeuralCoarseGridCorrection> coarse_grid_;
    std::unique_ptr<NeuralAMGCoarsening> neural_amg_;
    std::unique_ptr<NeuralDeflation> deflation_;
    
    // PETSc objects
    KSP ksp_ = nullptr;
    PC pc_ = nullptr;
    
    // Online learning buffer
    std::vector<std::pair<Tensor, Tensor>> recent_solves_;
    
    void setupPETSc(Mat A);
};

// =============================================================================
// Configuration Parsing
// =============================================================================

NeuralCoarseGridConfig parseNeuralCoarseGridConfig(
    const std::map<std::string, std::string>& config);
    
NeuralLinearSolver::SolverConfig parseNeuralSolverConfig(
    const std::map<std::string, std::string>& config);

} // namespace ML
} // namespace FSRM

#endif // FSRM_NEURAL_LINEAR_ALGEBRA_HPP
