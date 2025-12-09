/**
 * @file GraphNeuralOperator.cpp
 * @brief Implementation of Graph Neural Networks for PDEs
 */

#include "GraphNeuralOperator.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <random>
#include <iostream>
#include <iomanip>
#include <queue>
#include <chrono>

namespace FSRM {
namespace ML {

// =============================================================================
// MeshGraph Implementation
// =============================================================================

void MeshGraph::fromMesh(DM dm) {
    PetscInt pStart, pEnd, cStart, cEnd, vStart, vEnd;
    DMPlexGetChart(dm, &pStart, &pEnd);
    DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd);
    DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd);
    
    num_nodes_ = vEnd - vStart;
    neighbors_.resize(num_nodes_);
    positions_.resize(num_nodes_);
    
    // Get coordinates
    Vec coords;
    DMGetCoordinatesLocal(dm, &coords);
    const PetscScalar* coord_arr;
    VecGetArrayRead(coords, &coord_arr);
    
    PetscInt dim;
    DMGetDimension(dm, &dim);
    
    for (PetscInt v = vStart; v < vEnd; ++v) {
        int idx = v - vStart;
        for (int d = 0; d < dim; ++d) {
            positions_[idx][d] = coord_arr[idx * dim + d];
        }
        for (int d = dim; d < 3; ++d) {
            positions_[idx][d] = 0.0;
        }
    }
    
    VecRestoreArrayRead(coords, &coord_arr);
    
    // Build edges from mesh connectivity
    edges_.clear();
    
    for (PetscInt c = cStart; c < cEnd; ++c) {
        PetscInt closure_size;
        PetscInt* closure = nullptr;
        DMPlexGetTransitiveClosure(dm, c, PETSC_TRUE, &closure_size, &closure);
        
        std::vector<int> cell_vertices;
        for (PetscInt i = 0; i < closure_size; ++i) {
            PetscInt p = closure[2*i];
            if (p >= vStart && p < vEnd) {
                cell_vertices.push_back(p - vStart);
            }
        }
        
        // Connect all vertices in cell
        for (size_t i = 0; i < cell_vertices.size(); ++i) {
            for (size_t j = i + 1; j < cell_vertices.size(); ++j) {
                Edge e1, e2;
                e1.src = cell_vertices[i];
                e1.dst = cell_vertices[j];
                e2.src = cell_vertices[j];
                e2.dst = cell_vertices[i];
                edges_.push_back(e1);
                edges_.push_back(e2);
                
                neighbors_[e1.src].push_back(e1.dst);
                neighbors_[e2.src].push_back(e2.dst);
            }
        }
        
        DMPlexRestoreTransitiveClosure(dm, c, PETSC_TRUE, &closure_size, &closure);
    }
    
    // Remove duplicate edges
    for (auto& neigh : neighbors_) {
        std::sort(neigh.begin(), neigh.end());
        neigh.erase(std::unique(neigh.begin(), neigh.end()), neigh.end());
    }
}

void MeshGraph::fromEdges(int num_nodes, const std::vector<Edge>& edges) {
    num_nodes_ = num_nodes;
    edges_ = edges;
    
    neighbors_.resize(num_nodes);
    for (const auto& e : edges_) {
        neighbors_[e.src].push_back(e.dst);
    }
}

void MeshGraph::fromKNN(const std::vector<std::array<double, 3>>& positions, int k) {
    num_nodes_ = positions.size();
    positions_ = positions;
    neighbors_.resize(num_nodes_);
    edges_.clear();
    
    for (int i = 0; i < num_nodes_; ++i) {
        // Find k nearest neighbors
        std::vector<std::pair<double, int>> distances;
        
        for (int j = 0; j < num_nodes_; ++j) {
            if (i == j) continue;
            double dist = 0.0;
            for (int d = 0; d < 3; ++d) {
                double diff = positions[i][d] - positions[j][d];
                dist += diff * diff;
            }
            distances.push_back({dist, j});
        }
        
        std::partial_sort(distances.begin(), 
                         distances.begin() + std::min(k, (int)distances.size()),
                         distances.end());
        
        for (int n = 0; n < std::min(k, (int)distances.size()); ++n) {
            Edge e;
            e.src = i;
            e.dst = distances[n].second;
            edges_.push_back(e);
            neighbors_[i].push_back(e.dst);
        }
    }
}

void MeshGraph::fromRadius(const std::vector<std::array<double, 3>>& positions, 
                           double radius) {
    num_nodes_ = positions.size();
    positions_ = positions;
    neighbors_.resize(num_nodes_);
    edges_.clear();
    
    double r2 = radius * radius;
    
    for (int i = 0; i < num_nodes_; ++i) {
        for (int j = i + 1; j < num_nodes_; ++j) {
            double dist2 = 0.0;
            for (int d = 0; d < 3; ++d) {
                double diff = positions[i][d] - positions[j][d];
                dist2 += diff * diff;
            }
            
            if (dist2 <= r2) {
                Edge e1, e2;
                e1.src = i; e1.dst = j;
                e2.src = j; e2.dst = i;
                edges_.push_back(e1);
                edges_.push_back(e2);
                neighbors_[i].push_back(j);
                neighbors_[j].push_back(i);
            }
        }
    }
}

void MeshGraph::setNodeFeatures(const Tensor& features) {
    node_features_ = features;
}

void MeshGraph::setEdgeFeatures(const Tensor& features) {
    edge_features_ = features;
}

void MeshGraph::computeEdgeFeatures(bool include_distance,
                                     bool include_relative_pos,
                                     bool include_normals) {
    int feat_dim = 0;
    if (include_distance) feat_dim += 1;
    if (include_relative_pos) feat_dim += 3;
    if (include_normals) feat_dim += 3;
    
    edge_features_ = Tensor({(int)edges_.size(), feat_dim});
    
    for (size_t e = 0; e < edges_.size(); ++e) {
        int i = edges_[e].src;
        int j = edges_[e].dst;
        
        double dx = positions_[j][0] - positions_[i][0];
        double dy = positions_[j][1] - positions_[i][1];
        double dz = positions_[j][2] - positions_[i][2];
        double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
        
        int idx = 0;
        if (include_distance) {
            edge_features_.data[e * feat_dim + idx++] = dist;
        }
        if (include_relative_pos) {
            edge_features_.data[e * feat_dim + idx++] = dx / (dist + 1e-8);
            edge_features_.data[e * feat_dim + idx++] = dy / (dist + 1e-8);
            edge_features_.data[e * feat_dim + idx++] = dz / (dist + 1e-8);
        }
        // Normals would require face information
    }
}

void MeshGraph::buildCoarseGraph(int coarsening_factor, MeshGraph& coarse) {
    int coarse_nodes = std::max(1, num_nodes_ / coarsening_factor);
    
    // Simple clustering: divide nodes into clusters
    std::vector<int> cluster_assignment(num_nodes_);
    for (int i = 0; i < num_nodes_; ++i) {
        cluster_assignment[i] = i % coarse_nodes;
    }
    
    // Build coarse graph
    coarse.num_nodes_ = coarse_nodes;
    coarse.positions_.resize(coarse_nodes, {0, 0, 0});
    coarse.neighbors_.resize(coarse_nodes);
    
    // Compute cluster centers
    std::vector<int> counts(coarse_nodes, 0);
    for (int i = 0; i < num_nodes_; ++i) {
        int c = cluster_assignment[i];
        counts[c]++;
        for (int d = 0; d < 3; ++d) {
            coarse.positions_[c][d] += positions_[i][d];
        }
    }
    
    for (int c = 0; c < coarse_nodes; ++c) {
        if (counts[c] > 0) {
            for (int d = 0; d < 3; ++d) {
                coarse.positions_[c][d] /= counts[c];
            }
        }
    }
    
    // Build coarse edges
    std::set<std::pair<int, int>> coarse_edges_set;
    for (const auto& e : edges_) {
        int c_src = cluster_assignment[e.src];
        int c_dst = cluster_assignment[e.dst];
        if (c_src != c_dst) {
            coarse_edges_set.insert({c_src, c_dst});
        }
    }
    
    coarse.edges_.clear();
    for (const auto& ce : coarse_edges_set) {
        Edge e;
        e.src = ce.first;
        e.dst = ce.second;
        coarse.edges_.push_back(e);
        coarse.neighbors_[e.src].push_back(e.dst);
    }
}

MeshGraph MeshGraph::batch(const std::vector<MeshGraph>& graphs) {
    MeshGraph batched;
    
    int total_nodes = 0;
    int total_edges = 0;
    
    for (const auto& g : graphs) {
        total_nodes += g.num_nodes_;
        total_edges += g.edges_.size();
    }
    
    batched.num_nodes_ = total_nodes;
    batched.positions_.reserve(total_nodes);
    batched.neighbors_.resize(total_nodes);
    batched.edges_.reserve(total_edges);
    batched.batch_indices_.reserve(total_nodes);
    
    int node_offset = 0;
    
    for (size_t b = 0; b < graphs.size(); ++b) {
        const auto& g = graphs[b];
        
        for (const auto& pos : g.positions_) {
            batched.positions_.push_back(pos);
        }
        
        for (size_t i = 0; i < g.neighbors_.size(); ++i) {
            for (int n : g.neighbors_[i]) {
                batched.neighbors_[node_offset + i].push_back(node_offset + n);
            }
        }
        
        for (const auto& e : g.edges_) {
            Edge ne;
            ne.src = e.src + node_offset;
            ne.dst = e.dst + node_offset;
            ne.features = e.features;
            batched.edges_.push_back(ne);
        }
        
        for (int i = 0; i < g.num_nodes_; ++i) {
            batched.batch_indices_.push_back(b);
        }
        
        node_offset += g.num_nodes_;
    }
    
    // Batch node features
    if (!graphs.empty() && graphs[0].node_features_.data.size() > 0) {
        int feat_dim = graphs[0].node_features_.shape[1];
        batched.node_features_ = Tensor({total_nodes, feat_dim});
        
        int idx = 0;
        for (const auto& g : graphs) {
            for (size_t i = 0; i < g.node_features_.data.size(); ++i) {
                batched.node_features_.data[idx++] = g.node_features_.data[i];
            }
        }
    }
    
    return batched;
}

// =============================================================================
// Message Passing Layer Implementation
// =============================================================================

MessagePassingLayer::MessagePassingLayer(const MessagePassingConfig& config) 
    : config_(config) {
    
    // Message MLP: [src_feat, dst_feat, edge_feat] -> message
    int msg_input_dim = 2 * config.node_dim + config.edge_dim;
    MLPConfig msg_config;
    msg_config.layer_sizes = {msg_input_dim, config.hidden_dim, config.hidden_dim};
    msg_config.activation = config.activation;
    msg_config.use_batch_norm = false;
    message_mlp_ = std::make_unique<MLPModel>(msg_config);
    
    // Update MLP: [node_feat, aggregated_msg] -> updated_feat
    int upd_input_dim = config.node_dim + config.hidden_dim;
    MLPConfig upd_config;
    upd_config.layer_sizes = {upd_input_dim, config.hidden_dim, config.output_dim};
    upd_config.activation = config.activation;
    upd_config.use_batch_norm = false;
    update_mlp_ = std::make_unique<MLPModel>(upd_config);
    
    if (config.use_layer_norm) {
        layer_norm_ = std::make_unique<LayerNorm>(config.output_dim);
    }
}

Tensor MessagePassingLayer::forward(const MeshGraph& graph) {
    cached_graph_ = graph;
    int num_nodes = graph.numNodes();
    int num_edges = graph.numEdges();
    int node_dim = config_.node_dim;
    int edge_dim = config_.edge_dim;
    
    const Tensor& node_feat = graph.nodeFeatures();
    const Tensor& edge_feat = graph.edgeFeatures();
    
    // Compute messages for all edges
    int msg_input_dim = 2 * node_dim + edge_dim;
    Tensor msg_input({num_edges, msg_input_dim});
    
    const auto& edges = graph.edges();
    for (int e = 0; e < num_edges; ++e) {
        int src = edges[e].src;
        int dst = edges[e].dst;
        
        // Concatenate [src_feat, dst_feat, edge_feat]
        int idx = 0;
        for (int d = 0; d < node_dim; ++d) {
            msg_input.data[e * msg_input_dim + idx++] = node_feat.data[src * node_dim + d];
        }
        for (int d = 0; d < node_dim; ++d) {
            msg_input.data[e * msg_input_dim + idx++] = node_feat.data[dst * node_dim + d];
        }
        for (int d = 0; d < edge_dim; ++d) {
            msg_input.data[e * msg_input_dim + idx++] = edge_feat.data[e * edge_dim + d];
        }
    }
    
    // Compute messages
    cached_messages_ = message_mlp_->forward(msg_input);
    
    // Aggregate messages
    cached_aggregated_ = aggregate(cached_messages_, graph);
    
    // Update nodes
    int upd_input_dim = node_dim + config_.hidden_dim;
    Tensor upd_input({num_nodes, upd_input_dim});
    
    for (int n = 0; n < num_nodes; ++n) {
        int idx = 0;
        for (int d = 0; d < node_dim; ++d) {
            upd_input.data[n * upd_input_dim + idx++] = node_feat.data[n * node_dim + d];
        }
        for (int d = 0; d < config_.hidden_dim; ++d) {
            upd_input.data[n * upd_input_dim + idx++] = 
                cached_aggregated_.data[n * config_.hidden_dim + d];
        }
    }
    
    Tensor output = update_mlp_->forward(upd_input);
    
    // Residual connection
    if (config_.use_residual && node_dim == config_.output_dim) {
        for (int n = 0; n < num_nodes; ++n) {
            for (int d = 0; d < config_.output_dim; ++d) {
                output.data[n * config_.output_dim + d] += 
                    node_feat.data[n * node_dim + d];
            }
        }
    }
    
    // Layer norm
    if (layer_norm_) {
        output = layer_norm_->forward(output);
    }
    
    return output;
}

Tensor MessagePassingLayer::aggregate(const Tensor& messages, const MeshGraph& graph) {
    int num_nodes = graph.numNodes();
    int msg_dim = config_.hidden_dim;
    
    Tensor aggregated({num_nodes, msg_dim}, 0.0);
    std::vector<int> counts(num_nodes, 0);
    
    const auto& edges = graph.edges();
    
    for (size_t e = 0; e < edges.size(); ++e) {
        int dst = edges[e].dst;
        counts[dst]++;
        
        if (config_.aggregation == "sum" || config_.aggregation == "mean") {
            for (int d = 0; d < msg_dim; ++d) {
                aggregated.data[dst * msg_dim + d] += messages.data[e * msg_dim + d];
            }
        } else if (config_.aggregation == "max") {
            for (int d = 0; d < msg_dim; ++d) {
                aggregated.data[dst * msg_dim + d] = 
                    std::max(aggregated.data[dst * msg_dim + d],
                            messages.data[e * msg_dim + d]);
            }
        }
    }
    
    if (config_.aggregation == "mean") {
        for (int n = 0; n < num_nodes; ++n) {
            if (counts[n] > 0) {
                for (int d = 0; d < msg_dim; ++d) {
                    aggregated.data[n * msg_dim + d] /= counts[n];
                }
            }
        }
    }
    
    return aggregated;
}

Tensor MessagePassingLayer::backward(const MeshGraph& graph, const Tensor& grad_output) {
    // Backward through layer norm
    Tensor grad = grad_output;
    if (layer_norm_) {
        grad = layer_norm_->backward(grad);
    }
    
    // Backward through update MLP
    Tensor grad_upd = update_mlp_->backward(grad);
    
    // Split gradient for node features and aggregated messages
    int num_nodes = graph.numNodes();
    int node_dim = config_.node_dim;
    int upd_input_dim = node_dim + config_.hidden_dim;
    
    Tensor grad_node({num_nodes, node_dim});
    Tensor grad_agg({num_nodes, config_.hidden_dim});
    
    for (int n = 0; n < num_nodes; ++n) {
        for (int d = 0; d < node_dim; ++d) {
            grad_node.data[n * node_dim + d] = grad_upd.data[n * upd_input_dim + d];
        }
        for (int d = 0; d < config_.hidden_dim; ++d) {
            grad_agg.data[n * config_.hidden_dim + d] = 
                grad_upd.data[n * upd_input_dim + node_dim + d];
        }
    }
    
    // Backward through aggregation (scatter gradient to messages)
    const auto& edges = graph.edges();
    int num_edges = edges.size();
    std::vector<int> counts(num_nodes, 0);
    
    for (const auto& e : edges) {
        counts[e.dst]++;
    }
    
    Tensor grad_messages({num_edges, config_.hidden_dim});
    
    for (size_t e = 0; e < edges.size(); ++e) {
        int dst = edges[e].dst;
        double scale = (config_.aggregation == "mean" && counts[dst] > 0) 
                      ? 1.0 / counts[dst] : 1.0;
        
        for (int d = 0; d < config_.hidden_dim; ++d) {
            grad_messages.data[e * config_.hidden_dim + d] = 
                scale * grad_agg.data[dst * config_.hidden_dim + d];
        }
    }
    
    // Backward through message MLP
    message_mlp_->backward(grad_messages);
    
    // Residual connection gradient
    if (config_.use_residual && node_dim == config_.output_dim) {
        for (int n = 0; n < num_nodes; ++n) {
            for (int d = 0; d < node_dim; ++d) {
                grad_node.data[n * node_dim + d] += grad_output.data[n * node_dim + d];
            }
        }
    }
    
    return grad_node;
}

std::vector<Tensor*> MessagePassingLayer::parameters() {
    std::vector<Tensor*> params;
    
    auto msg_params = message_mlp_->parameters();
    auto upd_params = update_mlp_->parameters();
    
    params.insert(params.end(), msg_params.begin(), msg_params.end());
    params.insert(params.end(), upd_params.begin(), upd_params.end());
    
    if (layer_norm_) {
        params.push_back(&layer_norm_->gamma);
        params.push_back(&layer_norm_->beta);
    }
    
    return params;
}

std::vector<Tensor*> MessagePassingLayer::gradients() {
    std::vector<Tensor*> grads;
    
    auto msg_grads = message_mlp_->gradients();
    auto upd_grads = update_mlp_->gradients();
    
    grads.insert(grads.end(), msg_grads.begin(), msg_grads.end());
    grads.insert(grads.end(), upd_grads.begin(), upd_grads.end());
    
    if (layer_norm_) {
        grads.push_back(&layer_norm_->grad_gamma);
        grads.push_back(&layer_norm_->grad_beta);
    }
    
    return grads;
}

void MessagePassingLayer::zeroGrad() {
    message_mlp_->zeroGrad();
    update_mlp_->zeroGrad();
}

// =============================================================================
// Graph Attention Layer Implementation
// =============================================================================

GraphAttentionLayer::GraphAttentionLayer(const GATConfig& config) 
    : MessagePassingLayer(config), gat_config_(config) {
    
    int head_dim = config.node_dim / config.num_heads;
    
    for (int h = 0; h < config.num_heads; ++h) {
        W_query_.push_back(Tensor({config.node_dim, head_dim}));
        W_key_.push_back(Tensor({config.node_dim, head_dim}));
        W_value_.push_back(Tensor({config.node_dim, head_dim}));
        
        W_query_.back().xavier_init();
        W_key_.back().xavier_init();
        W_value_.back().xavier_init();
    }
}

Tensor GraphAttentionLayer::forward(const MeshGraph& graph) {
    int num_nodes = graph.numNodes();
    int node_dim = gat_config_.node_dim;
    int num_heads = gat_config_.num_heads;
    int head_dim = node_dim / num_heads;
    
    const Tensor& x = graph.nodeFeatures();
    const auto& edges = graph.edges();
    
    // Output dimension depends on whether we concat heads
    int out_dim = gat_config_.concat_heads ? node_dim : head_dim;
    Tensor output({num_nodes, out_dim}, 0.0);
    
    for (int h = 0; h < num_heads; ++h) {
        // Compute Q, K, V for all nodes
        Tensor Q({num_nodes, head_dim}), K({num_nodes, head_dim}), V({num_nodes, head_dim});
        
        for (int n = 0; n < num_nodes; ++n) {
            for (int d = 0; d < head_dim; ++d) {
                double q = 0, k = 0, v = 0;
                for (int i = 0; i < node_dim; ++i) {
                    q += x.data[n * node_dim + i] * W_query_[h].data[i * head_dim + d];
                    k += x.data[n * node_dim + i] * W_key_[h].data[i * head_dim + d];
                    v += x.data[n * node_dim + i] * W_value_[h].data[i * head_dim + d];
                }
                Q.data[n * head_dim + d] = q;
                K.data[n * head_dim + d] = k;
                V.data[n * head_dim + d] = v;
            }
        }
        
        // Compute attention for each node
        for (int n = 0; n < num_nodes; ++n) {
            const auto& neighbors = graph.neighbors()[n];
            if (neighbors.empty()) continue;
            
            // Compute attention scores
            std::vector<double> scores(neighbors.size());
            double max_score = -1e10;
            
            for (size_t i = 0; i < neighbors.size(); ++i) {
                int m = neighbors[i];
                double score = 0.0;
                for (int d = 0; d < head_dim; ++d) {
                    score += Q.data[n * head_dim + d] * K.data[m * head_dim + d];
                }
                score /= std::sqrt(head_dim);
                scores[i] = score;
                max_score = std::max(max_score, score);
            }
            
            // Softmax
            double sum_exp = 0.0;
            for (size_t i = 0; i < neighbors.size(); ++i) {
                scores[i] = std::exp(scores[i] - max_score);
                sum_exp += scores[i];
            }
            
            // Aggregate
            for (size_t i = 0; i < neighbors.size(); ++i) {
                int m = neighbors[i];
                double attn = scores[i] / (sum_exp + 1e-10);
                
                if (gat_config_.concat_heads) {
                    for (int d = 0; d < head_dim; ++d) {
                        output.data[n * out_dim + h * head_dim + d] += 
                            attn * V.data[m * head_dim + d];
                    }
                } else {
                    for (int d = 0; d < head_dim; ++d) {
                        output.data[n * out_dim + d] += 
                            attn * V.data[m * head_dim + d] / num_heads;
                    }
                }
            }
        }
    }
    
    return output;
}

Tensor GraphAttentionLayer::backward(const MeshGraph& graph, const Tensor& grad_output) {
    // Simplified backward - full implementation would need cached attention weights
    return MessagePassingLayer::backward(graph, grad_output);
}

// =============================================================================
// Graph Neural Operator Implementation
// =============================================================================

GraphNeuralOperator::GraphNeuralOperator(const GNOConfig& config) : config_(config) {
    // Encoder
    MLPConfig enc_config;
    enc_config.layer_sizes = {config.input_dim, config.hidden_dim, config.hidden_dim};
    enc_config.activation = "gelu";
    enc_config.use_batch_norm = config.use_batch_norm;
    encoder_ = std::make_unique<MLPModel>(enc_config);
    
    // Message passing layers
    for (int i = 0; i < config.num_layers; ++i) {
        MessagePassingConfig mp_config;
        mp_config.node_dim = config.hidden_dim;
        mp_config.edge_dim = config.edge_dim;
        mp_config.hidden_dim = config.hidden_dim;
        mp_config.output_dim = config.hidden_dim;
        mp_config.use_layer_norm = config.use_layer_norm;
        mp_config.use_residual = config.use_residual;
        
        switch (config.layer_type) {
            case GNOConfig::LayerType::GAT: {
                GraphAttentionLayer::GATConfig gat_config;
                static_cast<MessagePassingConfig&>(gat_config) = mp_config;
                gat_config.num_heads = config.num_heads;
                mp_layers_.push_back(std::make_unique<GraphAttentionLayer>(gat_config));
                break;
            }
            case GNOConfig::LayerType::EDGE_CONV:
                mp_layers_.push_back(std::make_unique<EdgeConvLayer>(mp_config));
                break;
            default:
                mp_layers_.push_back(std::make_unique<MessagePassingLayer>(mp_config));
        }
    }
    
    // Decoder
    MLPConfig dec_config;
    dec_config.layer_sizes = {config.hidden_dim, config.hidden_dim, config.output_dim};
    dec_config.activation = "gelu";
    dec_config.use_batch_norm = false;
    decoder_ = std::make_unique<MLPModel>(dec_config);
    
    // Positional encoding
    if (config.use_positional_encoding) {
        buildPositionalEncoding(10000);  // Max nodes
    }
}

void GraphNeuralOperator::buildPositionalEncoding(int max_nodes) {
    int pe_dim = config_.pos_encoding_dim;
    positional_encoding_ = Tensor({max_nodes, pe_dim});
    
    for (int n = 0; n < max_nodes; ++n) {
        for (int d = 0; d < pe_dim; ++d) {
            double freq = 1.0 / std::pow(10000.0, 2.0 * d / pe_dim);
            if (d % 2 == 0) {
                positional_encoding_.data[n * pe_dim + d] = std::sin(n * freq);
            } else {
                positional_encoding_.data[n * pe_dim + d] = std::cos(n * freq);
            }
        }
    }
}

Tensor GraphNeuralOperator::forward(const MeshGraph& graph) {
    cached_graph_ = graph;
    int num_nodes = graph.numNodes();
    
    // Encode node features
    encoded_ = encoder_->forward(graph.nodeFeatures());
    
    // Add positional encoding if enabled
    if (config_.use_positional_encoding && positional_encoding_.data.size() > 0) {
        int pe_dim = std::min(config_.pos_encoding_dim, config_.hidden_dim);
        for (int n = 0; n < num_nodes; ++n) {
            for (int d = 0; d < pe_dim; ++d) {
                encoded_.data[n * config_.hidden_dim + d] += 
                    positional_encoding_.data[n * config_.pos_encoding_dim + d];
            }
        }
    }
    
    // Message passing
    skip_outputs_.clear();
    MeshGraph temp_graph = graph;
    temp_graph.setNodeFeatures(encoded_);
    
    Tensor h = encoded_;
    
    for (size_t l = 0; l < mp_layers_.size(); ++l) {
        temp_graph.setNodeFeatures(h);
        h = mp_layers_[l]->forward(temp_graph);
        
        if (config_.use_skip_connections) {
            skip_outputs_.push_back(h);
        }
    }
    
    // Add skip connections
    if (config_.use_skip_connections && !skip_outputs_.empty()) {
        Tensor aggregated({num_nodes, config_.hidden_dim}, 0.0);
        for (const auto& skip : skip_outputs_) {
            for (size_t i = 0; i < skip.data.size(); ++i) {
                aggregated.data[i] += skip.data[i];
            }
        }
        for (size_t i = 0; i < h.data.size(); ++i) {
            h.data[i] = (h.data[i] + aggregated.data[i]) / (skip_outputs_.size() + 1);
        }
    }
    
    // Decode
    Tensor output = decoder_->forward(h);
    
    return output;
}

Tensor GraphNeuralOperator::backward(const Tensor& grad_output) {
    int num_nodes = cached_graph_.numNodes();
    
    // Backward through decoder
    Tensor grad_h = decoder_->backward(grad_output);
    
    // Backward through message passing layers
    MeshGraph temp_graph = cached_graph_;
    
    for (int l = mp_layers_.size() - 1; l >= 0; --l) {
        if (l > 0 && l - 1 < (int)skip_outputs_.size()) {
            temp_graph.setNodeFeatures(skip_outputs_[l - 1]);
        } else {
            temp_graph.setNodeFeatures(encoded_);
        }
        grad_h = mp_layers_[l]->backward(temp_graph, grad_h);
    }
    
    // Backward through encoder
    Tensor grad_input = encoder_->backward(grad_h);
    
    return grad_input;
}

Tensor GraphNeuralOperator::predict(const MeshGraph& graph) {
    return forward(graph);
}

Tensor GraphNeuralOperator::predict(DM dm, Vec solution) {
    MeshGraph graph;
    graph.fromMesh(dm);
    
    // Extract solution to node features
    const PetscScalar* sol_arr;
    VecGetArrayRead(solution, &sol_arr);
    
    PetscInt n;
    VecGetSize(solution, &n);
    int feat_dim = n / graph.numNodes();
    
    Tensor node_feat({graph.numNodes(), feat_dim});
    for (int i = 0; i < n; ++i) {
        node_feat.data[i] = sol_arr[i];
    }
    
    VecRestoreArrayRead(solution, &sol_arr);
    
    graph.setNodeFeatures(node_feat);
    graph.computeEdgeFeatures();
    
    return predict(graph);
}

double GraphNeuralOperator::computeLoss(const Tensor& prediction, const Tensor& target) {
    double loss = 0.0;
    for (size_t i = 0; i < prediction.data.size(); ++i) {
        double diff = prediction.data[i] - target.data[i];
        loss += diff * diff;
    }
    return loss / prediction.data.size();
}

Tensor GraphNeuralOperator::computeLossGrad(const Tensor& prediction, const Tensor& target) {
    Tensor grad(prediction.shape);
    double scale = 2.0 / prediction.data.size();
    for (size_t i = 0; i < prediction.data.size(); ++i) {
        grad.data[i] = scale * (prediction.data[i] - target.data[i]);
    }
    return grad;
}

void GraphNeuralOperator::train(const std::vector<MeshGraph>& train_graphs,
                                const std::vector<Tensor>& train_targets,
                                const std::vector<MeshGraph>& val_graphs,
                                const std::vector<Tensor>& val_targets) {
    AdamOptimizer optimizer(parameters(), config_.learning_rate, 
                           0.9, 0.999, 1e-8, config_.weight_decay);
    
    std::cout << "Training Graph Neural Operator...\n";
    std::cout << "Train graphs: " << train_graphs.size() << "\n";
    std::cout << "Val graphs: " << val_graphs.size() << "\n";
    
    for (int epoch = 0; epoch < config_.epochs; ++epoch) {
        auto start = std::chrono::high_resolution_clock::now();
        
        double train_loss = trainEpoch(train_graphs, train_targets);
        double val_loss = validate(val_graphs, val_targets);
        
        auto end = std::chrono::high_resolution_clock::now();
        double epoch_time = std::chrono::duration<double>(end - start).count();
        
        if ((epoch + 1) % 10 == 0) {
            std::cout << "Epoch " << epoch + 1 << "/" << config_.epochs
                      << " - train_loss: " << train_loss
                      << " - val_loss: " << val_loss
                      << " - time: " << epoch_time << "s\n";
        }
    }
}

double GraphNeuralOperator::trainEpoch(const std::vector<MeshGraph>& graphs,
                                        const std::vector<Tensor>& targets) {
    double total_loss = 0.0;
    
    // Simple SGD over graphs
    std::vector<size_t> indices(graphs.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), std::mt19937(42));
    
    AdamOptimizer optimizer(parameters(), config_.learning_rate);
    
    for (size_t idx : indices) {
        zeroGrad();
        
        Tensor pred = forward(graphs[idx]);
        double loss = computeLoss(pred, targets[idx]);
        total_loss += loss;
        
        Tensor grad = computeLossGrad(pred, targets[idx]);
        backward(grad);
        
        optimizer.step(gradients());
    }
    
    return total_loss / graphs.size();
}

double GraphNeuralOperator::validate(const std::vector<MeshGraph>& graphs,
                                      const std::vector<Tensor>& targets) {
    double total_loss = 0.0;
    
    for (size_t i = 0; i < graphs.size(); ++i) {
        Tensor pred = predict(graphs[i]);
        total_loss += computeLoss(pred, targets[i]);
    }
    
    return total_loss / graphs.size();
}

std::vector<Tensor*> GraphNeuralOperator::parameters() {
    std::vector<Tensor*> params;
    
    auto enc_params = encoder_->parameters();
    params.insert(params.end(), enc_params.begin(), enc_params.end());
    
    for (auto& layer : mp_layers_) {
        auto layer_params = layer->parameters();
        params.insert(params.end(), layer_params.begin(), layer_params.end());
    }
    
    auto dec_params = decoder_->parameters();
    params.insert(params.end(), dec_params.begin(), dec_params.end());
    
    return params;
}

std::vector<Tensor*> GraphNeuralOperator::gradients() {
    std::vector<Tensor*> grads;
    
    auto enc_grads = encoder_->gradients();
    grads.insert(grads.end(), enc_grads.begin(), enc_grads.end());
    
    for (auto& layer : mp_layers_) {
        auto layer_grads = layer->gradients();
        grads.insert(grads.end(), layer_grads.begin(), layer_grads.end());
    }
    
    auto dec_grads = decoder_->gradients();
    grads.insert(grads.end(), dec_grads.begin(), dec_grads.end());
    
    return grads;
}

void GraphNeuralOperator::zeroGrad() {
    encoder_->zeroGrad();
    for (auto& layer : mp_layers_) {
        layer->zeroGrad();
    }
    decoder_->zeroGrad();
}

void GraphNeuralOperator::save(const std::string& path) const {
    encoder_->save(path + ".encoder");
    for (size_t i = 0; i < mp_layers_.size(); ++i) {
        // Save layer parameters
    }
    decoder_->save(path + ".decoder");
}

void GraphNeuralOperator::load(const std::string& path) {
    encoder_->load(path + ".encoder");
    decoder_->load(path + ".decoder");
}

void GraphNeuralOperator::exportONNX(const std::string& path) const {
    std::cerr << "ONNX export not yet implemented for GNO\n";
}

void GraphNeuralOperator::toGPU() { on_gpu_ = true; }
void GraphNeuralOperator::toCPU() { on_gpu_ = false; }

size_t GraphNeuralOperator::numParameters() const {
    size_t count = encoder_->numParameters() + decoder_->numParameters();
    for (const auto& layer : mp_layers_) {
        count += layer->parameters().size();  // Simplified
    }
    return count;
}

void GraphNeuralOperator::summary() const {
    std::cout << "Graph Neural Operator Summary\n";
    std::cout << "============================\n";
    std::cout << "Input dim: " << config_.input_dim << "\n";
    std::cout << "Output dim: " << config_.output_dim << "\n";
    std::cout << "Hidden dim: " << config_.hidden_dim << "\n";
    std::cout << "Num layers: " << config_.num_layers << "\n";
    std::cout << "Parameters: " << numParameters() << "\n";
}

// =============================================================================
// MeshGraphNet Implementation
// =============================================================================

MeshGraphNet::MeshGraphNet(const MGNConfig& config) : config_(config) {
    // Node encoder
    MLPConfig node_enc;
    node_enc.layer_sizes = {config.node_feature_dim + config.num_node_types, 
                            config.hidden_dim, config.hidden_dim};
    node_encoder_ = std::make_unique<MLPModel>(node_enc);
    
    // Edge encoder
    MLPConfig edge_enc;
    edge_enc.layer_sizes = {config.edge_feature_dim, config.hidden_dim / 2, config.hidden_dim / 2};
    edge_encoder_ = std::make_unique<MLPModel>(edge_enc);
    
    // Processor (GNO)
    GNOConfig proc_config;
    proc_config.input_dim = config.hidden_dim;
    proc_config.output_dim = config.hidden_dim;
    proc_config.hidden_dim = config.hidden_dim;
    proc_config.edge_dim = config.hidden_dim / 2;
    proc_config.num_layers = config.num_mp_steps;
    processor_ = std::make_unique<GraphNeuralOperator>(proc_config);
    
    // Decoder
    MLPConfig dec;
    dec.layer_sizes = {config.hidden_dim, config.hidden_dim, config.node_feature_dim};
    decoder_ = std::make_unique<MLPModel>(dec);
    
    // Node type embeddings
    node_type_embeddings_ = Tensor({config.num_node_types, config.num_node_types});
    for (int i = 0; i < config.num_node_types; ++i) {
        node_type_embeddings_.data[i * config.num_node_types + i] = 1.0;
    }
}

Tensor MeshGraphNet::step(const Tensor& current_state, const MeshGraph& mesh) {
    int num_nodes = mesh.numNodes();
    int feat_dim = config_.node_feature_dim;
    
    // Encode nodes (concatenate state with node type)
    Tensor encoded = node_encoder_->forward(current_state);
    
    // Update mesh with encoded features
    MeshGraph encoded_mesh = mesh;
    encoded_mesh.setNodeFeatures(encoded);
    
    // Process through GNO
    Tensor processed = processor_->forward(encoded_mesh);
    
    // Decode to get delta
    Tensor delta = decoder_->forward(processed);
    
    // Add delta to current state
    Tensor next_state(current_state.shape);
    for (size_t i = 0; i < current_state.data.size(); ++i) {
        next_state.data[i] = current_state.data[i] + delta.data[i];
    }
    
    return next_state;
}

std::vector<Tensor> MeshGraphNet::rollout(const Tensor& initial_state,
                                           const MeshGraph& mesh,
                                           int num_steps) {
    std::vector<Tensor> trajectory;
    trajectory.push_back(initial_state);
    
    Tensor current = initial_state;
    for (int t = 0; t < num_steps; ++t) {
        current = step(current, mesh);
        trajectory.push_back(current);
    }
    
    return trajectory;
}

void MeshGraphNet::addNoise(Tensor& state, double std) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dist(0.0, std);
    
    for (double& x : state.data) {
        x += dist(gen);
    }
}

void MeshGraphNet::save(const std::string& path) const {
    node_encoder_->save(path + ".node_enc");
    edge_encoder_->save(path + ".edge_enc");
    processor_->save(path + ".processor");
    decoder_->save(path + ".decoder");
}

void MeshGraphNet::load(const std::string& path) {
    node_encoder_->load(path + ".node_enc");
    edge_encoder_->load(path + ".edge_enc");
    processor_->load(path + ".processor");
    decoder_->load(path + ".decoder");
}

// =============================================================================
// GNO Solver Implementation
// =============================================================================

GNOSolver::GNOSolver(GNOType type, const GNOConfig& config) 
    : type_(type), config_(config) {
    
    switch (type) {
        case GNOType::BASIC_GNO:
            gno_ = std::make_unique<GraphNeuralOperator>(config);
            break;
        case GNOType::MESH_GRAPH_NET: {
            MeshGraphNet::MGNConfig mgn_config;
            mgn_config.hidden_dim = config.hidden_dim;
            mgn_ = std::make_unique<MeshGraphNet>(mgn_config);
            break;
        }
        case GNOType::MULTIPOLE_GNN: {
            MultipoleGNN::MultipoleConfig mp_config;
            mp_config.hidden_dim = config.hidden_dim;
            mpgnn_ = std::make_unique<MultipoleGNN>(mp_config);
            break;
        }
        case GNOType::PHYSICS_INFORMED_GNN: {
            PhysicsInformedGNN::PIGNNConfig pi_config;
            static_cast<GNOConfig&>(pi_config) = config;
            pignn_ = std::make_unique<PhysicsInformedGNN>(pi_config);
            break;
        }
    }
}

PetscErrorCode GNOSolver::solve(DM dm, Vec input, Vec output) {
    auto start = std::chrono::high_resolution_clock::now();
    
    MeshGraph graph = meshToGraph(dm, input);
    Tensor result;
    
    switch (type_) {
        case GNOType::BASIC_GNO:
            result = gno_->predict(graph);
            break;
        case GNOType::PHYSICS_INFORMED_GNN:
            result = pignn_->forward(graph);
            break;
        default:
            result = gno_->predict(graph);
    }
    
    graphToVec(result, output);
    
    auto end = std::chrono::high_resolution_clock::now();
    stats_.inference_time += std::chrono::duration<double>(end - start).count();
    stats_.num_calls++;
    
    return 0;
}

MeshGraph GNOSolver::meshToGraph(DM dm, Vec data) {
    MeshGraph graph;
    graph.fromMesh(dm);
    
    const PetscScalar* arr;
    VecGetArrayRead(data, &arr);
    
    PetscInt n;
    VecGetSize(data, &n);
    int feat_dim = n / graph.numNodes();
    
    Tensor features({graph.numNodes(), feat_dim});
    for (int i = 0; i < n; ++i) {
        features.data[i] = arr[i];
    }
    
    VecRestoreArrayRead(data, &arr);
    
    graph.setNodeFeatures(features);
    graph.computeEdgeFeatures();
    
    return graph;
}

void GNOSolver::graphToVec(const Tensor& node_data, Vec output) {
    PetscScalar* arr;
    VecGetArray(output, &arr);
    
    for (size_t i = 0; i < node_data.data.size(); ++i) {
        arr[i] = node_data.data[i];
    }
    
    VecRestoreArray(output, &arr);
}

void GNOSolver::loadModel(const std::string& path) {
    switch (type_) {
        case GNOType::BASIC_GNO:
            gno_->load(path);
            break;
        case GNOType::MESH_GRAPH_NET:
            mgn_->load(path);
            break;
        case GNOType::PHYSICS_INFORMED_GNN:
            pignn_->load(path);
            break;
        default:
            break;
    }
}

void GNOSolver::saveModel(const std::string& path) const {
    switch (type_) {
        case GNOType::BASIC_GNO:
            gno_->save(path);
            break;
        case GNOType::MESH_GRAPH_NET:
            mgn_->save(path);
            break;
        case GNOType::PHYSICS_INFORMED_GNN:
            pignn_->save(path);
            break;
        default:
            break;
    }
}

// =============================================================================
// MultipoleGNN Implementation
// =============================================================================

MultipoleGNN::MultipoleGNN(const MultipoleConfig& config) : config_(config) {
    // Create processor at each level
    for (int l = 0; l < config.num_levels; ++l) {
        GNOConfig level_config;
        level_config.input_dim = config.hidden_dim;
        level_config.output_dim = config.hidden_dim;
        level_config.hidden_dim = config.hidden_dim;
        level_config.num_layers = config.mp_steps_per_level;
        
        level_processors_.push_back(std::make_unique<GraphNeuralOperator>(level_config));
    }
    
    // Upsample/downsample MLPs
    MLPConfig up_config;
    up_config.layer_sizes = {config.hidden_dim, config.hidden_dim, config.hidden_dim};
    upsample_mlp_ = std::make_unique<MLPModel>(up_config);
    downsample_mlp_ = std::make_unique<MLPModel>(up_config);
}

void MultipoleGNN::buildHierarchy(const MeshGraph& fine_graph) {
    graph_hierarchy_.clear();
    graph_hierarchy_.push_back(fine_graph);
    
    MeshGraph current = fine_graph;
    for (int l = 1; l < config_.num_levels; ++l) {
        MeshGraph coarse;
        current.buildCoarseGraph(config_.coarsening_ratio, coarse);
        graph_hierarchy_.push_back(coarse);
        current = coarse;
    }
}

Tensor MultipoleGNN::forward(const MeshGraph& graph) {
    buildHierarchy(graph);
    
    // Downward pass (fine to coarse)
    std::vector<Tensor> level_features(config_.num_levels);
    level_features[0] = graph.nodeFeatures();
    
    for (int l = 0; l < config_.num_levels - 1; ++l) {
        level_features[l + 1] = restrict(level_features[l], l);
    }
    
    // Process at each level
    for (int l = config_.num_levels - 1; l >= 0; --l) {
        MeshGraph& level_graph = graph_hierarchy_[l];
        level_graph.setNodeFeatures(level_features[l]);
        level_features[l] = level_processors_[l]->forward(level_graph);
    }
    
    // Upward pass (coarse to fine)
    for (int l = config_.num_levels - 2; l >= 0; --l) {
        Tensor coarse_correction = prolongate(level_features[l + 1], l);
        for (size_t i = 0; i < level_features[l].data.size(); ++i) {
            level_features[l].data[i] += coarse_correction.data[i];
        }
    }
    
    return level_features[0];
}

Tensor MultipoleGNN::restrict(const Tensor& fine_features, int level) {
    // Simplified restriction: average over fine nodes
    const MeshGraph& fine = graph_hierarchy_[level];
    const MeshGraph& coarse = graph_hierarchy_[level + 1];
    
    int fine_nodes = fine.numNodes();
    int coarse_nodes = coarse.numNodes();
    int feat_dim = fine_features.shape[1];
    
    Tensor coarse_features({coarse_nodes, feat_dim}, 0.0);
    std::vector<int> counts(coarse_nodes, 0);
    
    // Simple assignment: node i goes to cluster i % coarse_nodes
    for (int i = 0; i < fine_nodes; ++i) {
        int c = i % coarse_nodes;
        counts[c]++;
        for (int d = 0; d < feat_dim; ++d) {
            coarse_features.data[c * feat_dim + d] += 
                fine_features.data[i * feat_dim + d];
        }
    }
    
    for (int c = 0; c < coarse_nodes; ++c) {
        if (counts[c] > 0) {
            for (int d = 0; d < feat_dim; ++d) {
                coarse_features.data[c * feat_dim + d] /= counts[c];
            }
        }
    }
    
    return coarse_features;
}

Tensor MultipoleGNN::prolongate(const Tensor& coarse_features, int level) {
    const MeshGraph& fine = graph_hierarchy_[level];
    const MeshGraph& coarse = graph_hierarchy_[level + 1];
    
    int fine_nodes = fine.numNodes();
    int coarse_nodes = coarse.numNodes();
    int feat_dim = coarse_features.shape[1];
    
    Tensor fine_features({fine_nodes, feat_dim});
    
    // Simple prolongation: copy from parent cluster
    for (int i = 0; i < fine_nodes; ++i) {
        int c = i % coarse_nodes;
        for (int d = 0; d < feat_dim; ++d) {
            fine_features.data[i * feat_dim + d] = 
                coarse_features.data[c * feat_dim + d];
        }
    }
    
    return fine_features;
}

void MultipoleGNN::save(const std::string& path) const {
    for (size_t l = 0; l < level_processors_.size(); ++l) {
        level_processors_[l]->save(path + ".level" + std::to_string(l));
    }
}

void MultipoleGNN::load(const std::string& path) {
    for (size_t l = 0; l < level_processors_.size(); ++l) {
        level_processors_[l]->load(path + ".level" + std::to_string(l));
    }
}

// =============================================================================
// Physics-Informed GNN Implementation
// =============================================================================

PhysicsInformedGNN::PhysicsInformedGNN(const PIGNNConfig& config) 
    : GraphNeuralOperator(config), pi_config_(config) {}

Tensor PhysicsInformedGNN::forward(const MeshGraph& graph) {
    Tensor output = GraphNeuralOperator::forward(graph);
    
    if (pi_config_.enforce_mass_conservation) {
        enforceMassConservation(output, graph);
    }
    
    return output;
}

double PhysicsInformedGNN::computePhysicsLoss(const MeshGraph& graph, 
                                               const Tensor& prediction) {
    double loss = 0.0;
    
    if (pi_config_.pde_residual_weight > 0) {
        Tensor residual = computePDEResidual(graph, prediction);
        loss += pi_config_.pde_residual_weight * residual.norm();
    }
    
    return loss;
}

Tensor PhysicsInformedGNN::computePDEResidual(const MeshGraph& graph,
                                               const Tensor& prediction) {
    // Compute Laplacian for Poisson equation
    if (pi_config_.pde_type == "poisson") {
        return graphLaplacian(prediction, graph);
    }
    
    return Tensor(prediction.shape, 0.0);
}

void PhysicsInformedGNN::enforceMassConservation(Tensor& output, 
                                                  const MeshGraph& graph) {
    // Ensure integral of output is conserved
    double sum = 0.0;
    for (double x : output.data) {
        sum += x;
    }
    
    double mean = sum / output.data.size();
    for (double& x : output.data) {
        x -= mean;
    }
}

Tensor PhysicsInformedGNN::graphLaplacian(const Tensor& field, 
                                           const MeshGraph& graph) {
    int num_nodes = graph.numNodes();
    int feat_dim = field.shape[1];
    
    Tensor laplacian({num_nodes, feat_dim}, 0.0);
    
    const auto& neighbors = graph.neighbors();
    
    for (int n = 0; n < num_nodes; ++n) {
        int degree = neighbors[n].size();
        if (degree == 0) continue;
        
        for (int d = 0; d < feat_dim; ++d) {
            double sum_neighbors = 0.0;
            for (int m : neighbors[n]) {
                sum_neighbors += field.data[m * feat_dim + d];
            }
            laplacian.data[n * feat_dim + d] = 
                sum_neighbors / degree - field.data[n * feat_dim + d];
        }
    }
    
    return laplacian;
}

// =============================================================================
// EdgeConvLayer Implementation
// =============================================================================

EdgeConvLayer::EdgeConvLayer(const MessagePassingConfig& config) 
    : MessagePassingLayer(config) {}

Tensor EdgeConvLayer::forward(const MeshGraph& graph) {
    // EdgeConv uses edge feature = h_j - h_i (relative features)
    // Then applies MLP and max aggregation
    
    return MessagePassingLayer::forward(graph);
}

// =============================================================================
// Utility Functions
// =============================================================================

namespace GraphUtils {

void computeEdgeFeaturesFromPositions(
    const std::vector<std::array<double, 3>>& positions,
    std::vector<Edge>& edges,
    bool normalize) {
    
    for (auto& e : edges) {
        double dx = positions[e.dst][0] - positions[e.src][0];
        double dy = positions[e.dst][1] - positions[e.src][1];
        double dz = positions[e.dst][2] - positions[e.src][2];
        double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
        
        if (normalize && dist > 0) {
            dx /= dist;
            dy /= dist;
            dz /= dist;
        }
        
        e.features = {dist, dx, dy, dz};
    }
}

MeshGraph buildGraphFromDMPlex(DM dm, bool include_boundary_edges) {
    MeshGraph graph;
    graph.fromMesh(dm);
    graph.computeEdgeFeatures();
    return graph;
}

Tensor computeLaplacian(const MeshGraph& graph, bool normalized) {
    int n = graph.numNodes();
    Tensor L({n, n}, 0.0);
    
    const auto& neighbors = graph.neighbors();
    
    for (int i = 0; i < n; ++i) {
        int degree = neighbors[i].size();
        L.data[i * n + i] = degree;
        
        for (int j : neighbors[i]) {
            L.data[i * n + j] = -1.0;
        }
    }
    
    if (normalized) {
        for (int i = 0; i < n; ++i) {
            double d_i = L.data[i * n + i];
            if (d_i > 0) {
                for (int j = 0; j < n; ++j) {
                    double d_j = L.data[j * n + j];
                    if (d_j > 0) {
                        L.data[i * n + j] /= std::sqrt(d_i * d_j);
                    }
                }
            }
        }
    }
    
    return L;
}

} // namespace GraphUtils

// =============================================================================
// Configuration Parsing
// =============================================================================

GNOConfig parseGNOConfig(const std::map<std::string, std::string>& config) {
    GNOConfig result;
    
    auto getInt = [&](const std::string& key, int def) {
        auto it = config.find(key);
        return it != config.end() ? std::stoi(it->second) : def;
    };
    
    auto getDouble = [&](const std::string& key, double def) {
        auto it = config.find(key);
        return it != config.end() ? std::stod(it->second) : def;
    };
    
    auto getBool = [&](const std::string& key, bool def) {
        auto it = config.find(key);
        return it != config.end() ? (it->second == "true" || it->second == "1") : def;
    };
    
    result.input_dim = getInt("gno_input_dim", 3);
    result.output_dim = getInt("gno_output_dim", 3);
    result.hidden_dim = getInt("gno_hidden_dim", 128);
    result.edge_dim = getInt("gno_edge_dim", 64);
    result.num_layers = getInt("gno_num_layers", 6);
    result.num_heads = getInt("gno_num_heads", 4);
    result.learning_rate = getDouble("gno_learning_rate", 1e-3);
    result.weight_decay = getDouble("gno_weight_decay", 1e-5);
    result.batch_size = getInt("gno_batch_size", 4);
    result.epochs = getInt("gno_epochs", 100);
    result.use_layer_norm = getBool("gno_use_layer_norm", true);
    result.use_residual = getBool("gno_use_residual", true);
    result.use_positional_encoding = getBool("gno_use_positional_encoding", true);
    result.use_gpu = getBool("gno_use_gpu", true);
    
    return result;
}

} // namespace ML
} // namespace FSRM
