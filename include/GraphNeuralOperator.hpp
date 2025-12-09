/**
 * @file GraphNeuralOperator.hpp
 * @brief Graph Neural Networks for learning PDE solutions on unstructured meshes
 * 
 * Implements message-passing neural networks that operate directly on mesh
 * connectivity, enabling:
 * - Learning on arbitrary unstructured grids
 * - Natural handling of faults and fractures
 * - Adaptive mesh compatibility
 * - Multi-resolution learning
 * 
 * Key architectures:
 * - GraphNeuralOperator (GNO): Graph-based neural operator
 * - MeshGraphNet: For physics simulations on meshes
 * - Multipole GNN: Multi-scale graph networks
 * 
 * References:
 * - Li et al. "Neural Operator: Graph Kernel Network for PDEs"
 * - Pfaff et al. "Learning Mesh-Based Simulation with Graph Networks"
 */

#ifndef FSRM_GRAPH_NEURAL_OPERATOR_HPP
#define FSRM_GRAPH_NEURAL_OPERATOR_HPP

#include "FourierNeuralOperator.hpp"
#include "NeuralSurrogates.hpp"
#include <petsc.h>
#include <vector>
#include <array>
#include <memory>
#include <unordered_map>
#include <unordered_set>

namespace FSRM {
namespace ML {

// =============================================================================
// Forward Declarations
// =============================================================================

class MessagePassingLayer;
class EdgeConvLayer;
class GraphAttentionLayer;
class GraphNeuralOperator;
class MeshGraphNet;
class MultipoleGNN;

// =============================================================================
// Graph Data Structures
// =============================================================================

/**
 * @brief Edge representation for graph
 */
struct Edge {
    int src;            // Source node index
    int dst;            // Destination node index
    std::vector<double> features;  // Edge features (e.g., distance, relative position)
};

/**
 * @brief Graph structure for mesh representation
 */
class MeshGraph {
public:
    MeshGraph() = default;
    
    /**
     * @brief Create from mesh connectivity
     */
    void fromMesh(DM dm);
    
    /**
     * @brief Create from explicit edge list
     */
    void fromEdges(int num_nodes, const std::vector<Edge>& edges);
    
    /**
     * @brief Create k-nearest neighbor graph
     */
    void fromKNN(const std::vector<std::array<double, 3>>& positions, int k);
    
    /**
     * @brief Create radius graph (connect nodes within radius)
     */
    void fromRadius(const std::vector<std::array<double, 3>>& positions, double radius);
    
    // Graph properties
    int numNodes() const { return num_nodes_; }
    int numEdges() const { return edges_.size(); }
    
    // Access
    const std::vector<Edge>& edges() const { return edges_; }
    const std::vector<std::vector<int>>& neighbors() const { return neighbors_; }
    const std::vector<std::array<double, 3>>& nodePositions() const { return positions_; }
    
    // Node features
    void setNodeFeatures(const Tensor& features);
    const Tensor& nodeFeatures() const { return node_features_; }
    Tensor& nodeFeatures() { return node_features_; }
    
    // Edge features
    void setEdgeFeatures(const Tensor& features);
    const Tensor& edgeFeatures() const { return edge_features_; }
    Tensor& edgeFeatures() { return edge_features_; }
    
    // Graph construction helpers
    void computeEdgeFeatures(bool include_distance = true,
                            bool include_relative_pos = true,
                            bool include_normals = false);
    
    // Multi-level graph (for hierarchical GNN)
    void buildCoarseGraph(int coarsening_factor, MeshGraph& coarse);
    
    // Batching
    static MeshGraph batch(const std::vector<MeshGraph>& graphs);
    std::vector<MeshGraph> unbatch(const std::vector<int>& batch_sizes);
    
private:
    int num_nodes_ = 0;
    std::vector<Edge> edges_;
    std::vector<std::vector<int>> neighbors_;
    std::vector<std::array<double, 3>> positions_;
    
    Tensor node_features_;   // [num_nodes, node_feature_dim]
    Tensor edge_features_;   // [num_edges, edge_feature_dim]
    
    std::vector<int> batch_indices_;  // For batched graphs
};

// =============================================================================
// Message Passing Layer
// =============================================================================

/**
 * @brief Configuration for message passing
 */
struct MessagePassingConfig {
    int node_dim = 128;         // Node feature dimension
    int edge_dim = 64;          // Edge feature dimension
    int hidden_dim = 256;       // Hidden dimension
    int output_dim = 128;       // Output node dimension
    
    std::string aggregation = "mean";  // mean, sum, max
    std::string activation = "gelu";
    bool use_layer_norm = true;
    bool use_residual = true;
};

/**
 * @brief Basic message passing layer
 * 
 * For each node i:
 *   m_ij = message(x_i, x_j, e_ij)
 *   x_i' = update(x_i, aggregate({m_ij : j in N(i)}))
 */
class MessagePassingLayer {
public:
    MessagePassingLayer(const MessagePassingConfig& config);
    
    /**
     * @brief Forward pass
     * 
     * @param graph Input graph with node and edge features
     * @return Updated node features
     */
    virtual Tensor forward(const MeshGraph& graph);
    
    /**
     * @brief Backward pass
     */
    virtual Tensor backward(const MeshGraph& graph, const Tensor& grad_output);
    
    // Parameters
    std::vector<Tensor*> parameters();
    std::vector<Tensor*> gradients();
    void zeroGrad();
    
protected:
    MessagePassingConfig config_;
    
    // Message function: MLP(concat(x_i, x_j, e_ij)) -> message
    std::unique_ptr<MLPModel> message_mlp_;
    
    // Update function: MLP(concat(x_i, aggregated_message)) -> x_i'
    std::unique_ptr<MLPModel> update_mlp_;
    
    // Optional layer norm
    std::unique_ptr<LayerNorm> layer_norm_;
    
    // Cached for backward
    Tensor cached_messages_;
    Tensor cached_aggregated_;
    MeshGraph cached_graph_;
    
    // Aggregation functions
    Tensor aggregate(const Tensor& messages, const MeshGraph& graph);
};

/**
 * @brief Edge convolution layer (PointNet++ style)
 */
class EdgeConvLayer : public MessagePassingLayer {
public:
    EdgeConvLayer(const MessagePassingConfig& config);
    
    Tensor forward(const MeshGraph& graph) override;
    
private:
    // EdgeConv uses dynamic graph construction
    bool dynamic_graph_ = false;
    int k_neighbors_ = 20;
};

/**
 * @brief Graph attention layer
 */
class GraphAttentionLayer : public MessagePassingLayer {
public:
    struct GATConfig : MessagePassingConfig {
        int num_heads = 4;
        double dropout = 0.0;
        bool concat_heads = true;
    };
    
    GraphAttentionLayer(const GATConfig& config);
    
    Tensor forward(const MeshGraph& graph) override;
    Tensor backward(const MeshGraph& graph, const Tensor& grad_output) override;
    
private:
    GATConfig gat_config_;
    
    // Attention parameters per head
    std::vector<Tensor> W_query_;
    std::vector<Tensor> W_key_;
    std::vector<Tensor> W_value_;
    Tensor attention_weights_;
    
    Tensor computeAttention(const Tensor& query, const Tensor& key, 
                           const MeshGraph& graph);
};

// =============================================================================
// Graph Neural Operator (GNO)
// =============================================================================

/**
 * @brief Configuration for Graph Neural Operator
 */
struct GNOConfig {
    // Architecture
    int input_dim = 3;          // Input node features
    int output_dim = 3;         // Output node features
    int hidden_dim = 128;       // Hidden dimension
    int edge_dim = 64;          // Edge feature dimension
    int num_layers = 6;         // Number of message passing layers
    
    // Layer type
    enum class LayerType {
        BASIC,          // Basic message passing
        EDGE_CONV,      // Edge convolution
        GAT,            // Graph attention
        TRANSFORMER     // Graph transformer
    };
    LayerType layer_type = LayerType::BASIC;
    
    // Attention (for GAT/Transformer)
    int num_heads = 4;
    
    // Training
    double learning_rate = 1e-3;
    double weight_decay = 1e-5;
    int batch_size = 4;
    int epochs = 100;
    
    // Normalization
    bool use_layer_norm = true;
    bool use_batch_norm = false;
    
    // Residual connections
    bool use_residual = true;
    bool use_skip_connections = true;
    
    // Position encoding
    bool use_positional_encoding = true;
    int pos_encoding_dim = 32;
    
    // GPU
    bool use_gpu = true;
    int gpu_device_id = 0;
};

/**
 * @brief Graph Neural Operator for learning operators on unstructured meshes
 * 
 * Unlike FNO which requires regular grids, GNO operates on arbitrary
 * graph structures, making it suitable for:
 * - Complex geometries with faults and fractures
 * - Adaptive meshes
 * - Unstructured finite element meshes
 */
class GraphNeuralOperator {
public:
    GraphNeuralOperator(const GNOConfig& config);
    
    /**
     * @brief Forward pass
     * 
     * @param graph Mesh graph with input node features
     * @return Tensor of output node features
     */
    Tensor forward(const MeshGraph& graph);
    
    /**
     * @brief Backward pass
     */
    Tensor backward(const Tensor& grad_output);
    
    /**
     * @brief Predict (no gradient computation)
     */
    Tensor predict(const MeshGraph& graph);
    
    /**
     * @brief Predict with mesh data
     */
    Tensor predict(DM dm, Vec solution);
    
    // Loss computation
    double computeLoss(const Tensor& prediction, const Tensor& target);
    Tensor computeLossGrad(const Tensor& prediction, const Tensor& target);
    
    // Training
    void train(const std::vector<MeshGraph>& train_graphs,
              const std::vector<Tensor>& train_targets,
              const std::vector<MeshGraph>& val_graphs,
              const std::vector<Tensor>& val_targets);
    
    double trainEpoch(const std::vector<MeshGraph>& graphs,
                     const std::vector<Tensor>& targets);
    
    double validate(const std::vector<MeshGraph>& graphs,
                   const std::vector<Tensor>& targets);
    
    // Parameters
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
    bool isOnGPU() const { return on_gpu_; }
    
    // Info
    size_t numParameters() const;
    void summary() const;
    
private:
    GNOConfig config_;
    bool on_gpu_ = false;
    
    // Encoder: map input features to hidden dimension
    std::unique_ptr<MLPModel> encoder_;
    
    // Message passing layers
    std::vector<std::unique_ptr<MessagePassingLayer>> mp_layers_;
    
    // Decoder: map hidden features to output
    std::unique_ptr<MLPModel> decoder_;
    
    // Positional encoding
    Tensor positional_encoding_;
    
    // Skip connections
    std::vector<Tensor> skip_outputs_;
    
    // Cached for backward
    MeshGraph cached_graph_;
    Tensor encoded_;
    
    void buildPositionalEncoding(int max_nodes);
};

// =============================================================================
// MeshGraphNet
// =============================================================================

/**
 * @brief MeshGraphNet for physics simulation (Pfaff et al.)
 * 
 * Specialized for mesh-based simulations with:
 * - Node types (interior, boundary, source, etc.)
 * - World-space and mesh-space edges
 * - Multi-step rollout training
 */
class MeshGraphNet {
public:
    struct MGNConfig {
        int node_feature_dim = 16;
        int edge_feature_dim = 8;
        int hidden_dim = 128;
        int num_mp_steps = 15;
        
        // Node types
        int num_node_types = 4;  // e.g., NORMAL, OBSTACLE, INLET, OUTLET
        
        // Edge types
        enum class EdgeType { MESH, WORLD };
        double world_edge_radius = 0.1;
        
        // Rollout
        int rollout_steps = 1;    // Training: usually 1
        int eval_rollout = 100;   // Evaluation: longer rollout
        
        // Noise for training
        double noise_std = 0.003;
    };
    
    MeshGraphNet(const MGNConfig& config);
    
    /**
     * @brief Single step prediction
     * 
     * @param current_state Current node features (positions, velocities, etc.)
     * @param mesh Mesh graph structure
     * @return Next state prediction
     */
    Tensor step(const Tensor& current_state, const MeshGraph& mesh);
    
    /**
     * @brief Multi-step rollout
     */
    std::vector<Tensor> rollout(const Tensor& initial_state, 
                                const MeshGraph& mesh,
                                int num_steps);
    
    /**
     * @brief Training with pushforward trick
     */
    void trainWithRollout(const std::vector<std::vector<Tensor>>& trajectories,
                         const MeshGraph& mesh,
                         int epochs);
    
    // Model I/O
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    MGNConfig config_;
    
    // Encoder networks
    std::unique_ptr<MLPModel> node_encoder_;
    std::unique_ptr<MLPModel> edge_encoder_;
    
    // Processor (message passing)
    std::unique_ptr<GraphNeuralOperator> processor_;
    
    // Decoder
    std::unique_ptr<MLPModel> decoder_;
    
    // Node type embeddings
    Tensor node_type_embeddings_;
    
    // Add noise for training
    void addNoise(Tensor& state, double std);
};

// =============================================================================
// Multipole Graph Neural Network
// =============================================================================

/**
 * @brief Multipole GNN for multi-scale learning
 * 
 * Uses hierarchical graph structure to capture both local and global
 * interactions efficiently (O(N log N) complexity).
 */
class MultipoleGNN {
public:
    struct MultipoleConfig {
        int input_dim = 3;
        int output_dim = 3;
        int hidden_dim = 128;
        
        // Hierarchy levels
        int num_levels = 4;
        int coarsening_ratio = 4;  // Nodes reduce by this factor each level
        
        // Per-level message passing
        int mp_steps_per_level = 2;
        
        // Upsampling/downsampling
        std::string upsample_method = "interpolate";  // interpolate, mlp
        std::string downsample_method = "pool";  // pool, mlp
    };
    
    MultipoleGNN(const MultipoleConfig& config);
    
    /**
     * @brief Forward pass through hierarchy
     */
    Tensor forward(const MeshGraph& graph);
    
    /**
     * @brief Build hierarchical graph structure
     */
    void buildHierarchy(const MeshGraph& fine_graph);
    
    // Training
    void train(const std::vector<MeshGraph>& graphs,
              const std::vector<Tensor>& targets,
              int epochs);
    
    // Model I/O
    void save(const std::string& path) const;
    void load(const std::string& path);
    
private:
    MultipoleConfig config_;
    
    // Hierarchy of graphs
    std::vector<MeshGraph> graph_hierarchy_;
    
    // Prolongation/restriction operators between levels
    std::vector<Tensor> prolongation_weights_;
    std::vector<Tensor> restriction_weights_;
    
    // GNO at each level
    std::vector<std::unique_ptr<GraphNeuralOperator>> level_processors_;
    
    // Cross-level message passing
    std::unique_ptr<MLPModel> upsample_mlp_;
    std::unique_ptr<MLPModel> downsample_mlp_;
    
    Tensor restrict(const Tensor& fine_features, int level);
    Tensor prolongate(const Tensor& coarse_features, int level);
};

// =============================================================================
// Physics-Informed Graph Networks
// =============================================================================

/**
 * @brief Physics-Informed GNN with conservation law enforcement
 */
class PhysicsInformedGNN : public GraphNeuralOperator {
public:
    struct PIGNNConfig : GNOConfig {
        // Physics constraints
        bool enforce_mass_conservation = true;
        bool enforce_momentum_conservation = false;
        bool enforce_energy_conservation = false;
        
        // Constraint weights
        double conservation_weight = 1.0;
        double pde_residual_weight = 0.1;
        
        // PDE type
        std::string pde_type = "poisson";  // poisson, wave, navier_stokes
    };
    
    PhysicsInformedGNN(const PIGNNConfig& config);
    
    /**
     * @brief Forward with physics constraints
     */
    Tensor forward(const MeshGraph& graph);
    
    /**
     * @brief Compute physics loss terms
     */
    double computePhysicsLoss(const MeshGraph& graph, const Tensor& prediction);
    
    /**
     * @brief Compute PDE residual on graph
     */
    Tensor computePDEResidual(const MeshGraph& graph, const Tensor& prediction);
    
private:
    PIGNNConfig pi_config_;
    
    // Conservation enforcement
    void enforceMassConservation(Tensor& output, const MeshGraph& graph);
    
    // Gradient computation on graph (for PDE residual)
    Tensor graphGradient(const Tensor& field, const MeshGraph& graph);
    Tensor graphLaplacian(const Tensor& field, const MeshGraph& graph);
    Tensor graphDivergence(const Tensor& vector_field, const MeshGraph& graph);
};

// =============================================================================
// GNO Solver Interface
// =============================================================================

/**
 * @brief Unified solver interface for GNO models
 */
class GNOSolver {
public:
    enum class GNOType {
        BASIC_GNO,
        MESH_GRAPH_NET,
        MULTIPOLE_GNN,
        PHYSICS_INFORMED_GNN
    };
    
    GNOSolver(GNOType type, const GNOConfig& config);
    
    /**
     * @brief Solve PDE on mesh
     */
    PetscErrorCode solve(DM dm, Vec input, Vec output);
    
    /**
     * @brief Time stepping with GNO
     */
    PetscErrorCode step(DM dm, Vec current, double dt, Vec next);
    
    /**
     * @brief Train from simulation data
     */
    void trainFromSimulations(const std::vector<DM>& meshes,
                             const std::vector<Vec>& inputs,
                             const std::vector<Vec>& outputs);
    
    // Model management
    void loadModel(const std::string& path);
    void saveModel(const std::string& path) const;
    
    // Statistics
    struct SolverStats {
        double inference_time = 0.0;
        int num_calls = 0;
        double avg_error = 0.0;
    };
    
    const SolverStats& getStats() const { return stats_; }
    
private:
    GNOType type_;
    GNOConfig config_;
    
    std::unique_ptr<GraphNeuralOperator> gno_;
    std::unique_ptr<MeshGraphNet> mgn_;
    std::unique_ptr<MultipoleGNN> mpgnn_;
    std::unique_ptr<PhysicsInformedGNN> pignn_;
    
    SolverStats stats_;
    
    // Convert PETSc objects to graph
    MeshGraph meshToGraph(DM dm, Vec data);
    void graphToVec(const Tensor& node_data, Vec output);
};

// =============================================================================
// Utility Functions
// =============================================================================

namespace GraphUtils {

/**
 * @brief Compute edge features from node positions
 */
void computeEdgeFeaturesFromPositions(
    const std::vector<std::array<double, 3>>& positions,
    std::vector<Edge>& edges,
    bool normalize = true);

/**
 * @brief Build graph from DMPlex
 */
MeshGraph buildGraphFromDMPlex(DM dm, bool include_boundary_edges = true);

/**
 * @brief Coarsen graph using clustering
 */
MeshGraph coarsenGraph(const MeshGraph& graph, int target_nodes,
                       const std::string& method = "graclus");

/**
 * @brief Compute graph Laplacian
 */
Tensor computeLaplacian(const MeshGraph& graph, bool normalized = true);

/**
 * @brief Random walk positional encoding
 */
Tensor randomWalkPE(const MeshGraph& graph, int walk_length);

/**
 * @brief Spectral positional encoding (eigenvectors of Laplacian)
 */
Tensor spectralPE(const MeshGraph& graph, int num_eigenvectors);

} // namespace GraphUtils

// =============================================================================
// Configuration Parsing
// =============================================================================

GNOConfig parseGNOConfig(const std::map<std::string, std::string>& config);

} // namespace ML
} // namespace FSRM

#endif // FSRM_GRAPH_NEURAL_OPERATOR_HPP
