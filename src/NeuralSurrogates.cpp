/**
 * @file NeuralSurrogates.cpp
 * @brief Implementation of neural network surrogates
 */

#include "NeuralSurrogates.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <random>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>

namespace FSRM {
namespace ML {

// =============================================================================
// Dense Layer Implementation
// =============================================================================

DenseLayer::DenseLayer(int in_features, int out_features, bool bias)
    : has_bias(bias) {
    weight = Tensor({out_features, in_features});
    grad_weight = Tensor({out_features, in_features});
    
    if (has_bias) {
        bias_vec = Tensor({out_features}, 0.0);
        grad_bias = Tensor({out_features});
    }
    
    initWeights("kaiming");
}

void DenseLayer::initWeights(const std::string& method) {
    if (method == "xavier") {
        weight.xavier_init();
    } else if (method == "kaiming") {
        weight.kaiming_init();
    } else if (method == "orthogonal") {
        // Simplified orthogonal init
        weight.randn(0.0, 1.0);
    }
}

Tensor DenseLayer::forward(const Tensor& x) {
    last_input = x;
    
    int in_f = weight.shape[1];
    int out_f = weight.shape[0];
    size_t batch_size = x.numel() / in_f;
    
    Tensor output({(int)batch_size, out_f});
    
    for (size_t b = 0; b < batch_size; ++b) {
        for (int o = 0; o < out_f; ++o) {
            double sum = has_bias ? bias_vec(o) : 0.0;
            for (int i = 0; i < in_f; ++i) {
                sum += weight(o, i) * x.data[b * in_f + i];
            }
            output.data[b * out_f + o] = sum;
        }
    }
    
    return output;
}

Tensor DenseLayer::backward(const Tensor& grad_output) {
    int in_f = weight.shape[1];
    int out_f = weight.shape[0];
    size_t batch_size = grad_output.numel() / out_f;
    
    // Gradient w.r.t. weight
    grad_weight.zeros();
    for (size_t b = 0; b < batch_size; ++b) {
        for (int o = 0; o < out_f; ++o) {
            for (int i = 0; i < in_f; ++i) {
                grad_weight(o, i) += grad_output.data[b * out_f + o] * 
                                    last_input.data[b * in_f + i];
            }
        }
    }
    
    // Gradient w.r.t. bias
    if (has_bias) {
        grad_bias.zeros();
        for (size_t b = 0; b < batch_size; ++b) {
            for (int o = 0; o < out_f; ++o) {
                grad_bias(o) += grad_output.data[b * out_f + o];
            }
        }
    }
    
    // Gradient w.r.t. input
    Tensor grad_input({(int)batch_size, in_f});
    for (size_t b = 0; b < batch_size; ++b) {
        for (int i = 0; i < in_f; ++i) {
            double sum = 0.0;
            for (int o = 0; o < out_f; ++o) {
                sum += grad_output.data[b * out_f + o] * weight(o, i);
            }
            grad_input.data[b * in_f + i] = sum;
        }
    }
    
    return grad_input;
}

// =============================================================================
// Batch Normalization Implementation
// =============================================================================

BatchNorm1d::BatchNorm1d(int num_features, double eps, double momentum)
    : num_features_(num_features), eps_(eps), momentum_(momentum) {
    
    gamma = Tensor({num_features}, 1.0);
    beta = Tensor({num_features}, 0.0);
    running_mean = Tensor({num_features}, 0.0);
    running_var = Tensor({num_features}, 1.0);
    
    grad_gamma = Tensor({num_features});
    grad_beta = Tensor({num_features});
}

Tensor BatchNorm1d::forward(const Tensor& x, bool training) {
    size_t batch_size = x.numel() / num_features_;
    Tensor output(x.shape);
    
    cached_std_ = Tensor({num_features_});
    cached_normalized_ = Tensor(x.shape);
    
    if (training && batch_size > 1) {
        // Compute batch mean and variance
        Tensor batch_mean({num_features_}, 0.0);
        Tensor batch_var({num_features_}, 0.0);
        
        for (int f = 0; f < num_features_; ++f) {
            double sum = 0.0;
            for (size_t b = 0; b < batch_size; ++b) {
                sum += x.data[b * num_features_ + f];
            }
            batch_mean(f) = sum / batch_size;
            
            double var_sum = 0.0;
            for (size_t b = 0; b < batch_size; ++b) {
                double diff = x.data[b * num_features_ + f] - batch_mean(f);
                var_sum += diff * diff;
            }
            batch_var(f) = var_sum / batch_size;
        }
        
        // Update running stats
        for (int f = 0; f < num_features_; ++f) {
            running_mean(f) = (1 - momentum_) * running_mean(f) + momentum_ * batch_mean(f);
            running_var(f) = (1 - momentum_) * running_var(f) + momentum_ * batch_var(f);
        }
        
        // Normalize
        for (int f = 0; f < num_features_; ++f) {
            cached_std_(f) = std::sqrt(batch_var(f) + eps_);
            for (size_t b = 0; b < batch_size; ++b) {
                size_t idx = b * num_features_ + f;
                cached_normalized_.data[idx] = (x.data[idx] - batch_mean(f)) / cached_std_(f);
                output.data[idx] = gamma(f) * cached_normalized_.data[idx] + beta(f);
            }
        }
    } else {
        // Use running stats
        for (int f = 0; f < num_features_; ++f) {
            cached_std_(f) = std::sqrt(running_var(f) + eps_);
            for (size_t b = 0; b < batch_size; ++b) {
                size_t idx = b * num_features_ + f;
                cached_normalized_.data[idx] = (x.data[idx] - running_mean(f)) / cached_std_(f);
                output.data[idx] = gamma(f) * cached_normalized_.data[idx] + beta(f);
            }
        }
    }
    
    return output;
}

Tensor BatchNorm1d::backward(const Tensor& grad_output) {
    size_t batch_size = grad_output.numel() / num_features_;
    Tensor grad_input(grad_output.shape);
    
    // Gradient w.r.t. gamma and beta
    grad_gamma.zeros();
    grad_beta.zeros();
    
    for (int f = 0; f < num_features_; ++f) {
        for (size_t b = 0; b < batch_size; ++b) {
            size_t idx = b * num_features_ + f;
            grad_gamma(f) += grad_output.data[idx] * cached_normalized_.data[idx];
            grad_beta(f) += grad_output.data[idx];
        }
    }
    
    // Simplified gradient w.r.t. input
    for (int f = 0; f < num_features_; ++f) {
        for (size_t b = 0; b < batch_size; ++b) {
            size_t idx = b * num_features_ + f;
            grad_input.data[idx] = grad_output.data[idx] * gamma(f) / cached_std_(f);
        }
    }
    
    return grad_input;
}

// =============================================================================
// Dropout Implementation
// =============================================================================

Dropout::Dropout(double p) : p_(p) {}

Tensor Dropout::forward(const Tensor& x, bool training) {
    if (!training || p_ == 0.0) {
        return x;
    }
    
    Tensor output(x.shape);
    mask_ = Tensor(x.shape);
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::bernoulli_distribution dist(1.0 - p_);
    
    double scale = 1.0 / (1.0 - p_);
    
    for (size_t i = 0; i < x.data.size(); ++i) {
        mask_.data[i] = dist(gen) ? scale : 0.0;
        output.data[i] = x.data[i] * mask_.data[i];
    }
    
    return output;
}

Tensor Dropout::backward(const Tensor& grad_output) {
    Tensor grad_input(grad_output.shape);
    for (size_t i = 0; i < grad_output.data.size(); ++i) {
        grad_input.data[i] = grad_output.data[i] * mask_.data[i];
    }
    return grad_input;
}

// =============================================================================
// MLP Model Implementation
// =============================================================================

MLPModel::MLPModel(const MLPConfig& config) : config_(config) {
    if (config_.layer_sizes.size() < 2) {
        throw std::runtime_error("MLP needs at least input and output layers");
    }
    
    // Create dense layers
    for (size_t i = 0; i < config_.layer_sizes.size() - 1; ++i) {
        dense_layers_.push_back(std::make_unique<DenseLayer>(
            config_.layer_sizes[i], config_.layer_sizes[i + 1], true));
        
        // Batch norm for hidden layers
        if (config_.use_batch_norm && i < config_.layer_sizes.size() - 2) {
            bn_layers_.push_back(std::make_unique<BatchNorm1d>(
                config_.layer_sizes[i + 1]));
        }
        
        // Dropout for hidden layers
        if (config_.use_dropout && i < config_.layer_sizes.size() - 2) {
            dropout_layers_.push_back(std::make_unique<Dropout>(
                config_.dropout_rate));
        }
    }
    
    activation_fn_ = Activation::getActivation(config_.activation);
    activation_backward_ = Activation::getActivationBackward(config_.activation);
}

Tensor MLPModel::forward(const Tensor& x) {
    layer_outputs_.clear();
    Tensor h = x;
    
    for (size_t i = 0; i < dense_layers_.size(); ++i) {
        h = dense_layers_[i]->forward(h);
        layer_outputs_.push_back(h);
        
        // Not last layer
        if (i < dense_layers_.size() - 1) {
            // Batch norm
            if (config_.use_batch_norm && i < bn_layers_.size()) {
                h = bn_layers_[i]->forward(h, training_);
            }
            
            // Activation
            h = activation_fn_(h);
            layer_outputs_.push_back(h);
            
            // Dropout
            if (config_.use_dropout && i < dropout_layers_.size()) {
                h = dropout_layers_[i]->forward(h, training_);
            }
        }
    }
    
    return h;
}

Tensor MLPModel::backward(const Tensor& grad_output) {
    Tensor grad = grad_output;
    
    for (int i = dense_layers_.size() - 1; i >= 0; --i) {
        if (i < (int)dense_layers_.size() - 1) {
            // Dropout backward
            if (config_.use_dropout && i < (int)dropout_layers_.size()) {
                grad = dropout_layers_[i]->backward(grad);
            }
            
            // Activation backward
            size_t layer_idx = i * 2 + 1;  // After dense, after BN
            if (layer_idx < layer_outputs_.size()) {
                grad = activation_backward_(layer_outputs_[layer_idx - 1], grad);
            }
            
            // Batch norm backward
            if (config_.use_batch_norm && i < (int)bn_layers_.size()) {
                grad = bn_layers_[i]->backward(grad);
            }
        }
        
        grad = dense_layers_[i]->backward(grad);
    }
    
    return grad;
}

Tensor MLPModel::predict(const Tensor& x) {
    bool was_training = training_;
    training_ = false;
    Tensor result = forward(x);
    training_ = was_training;
    return result;
}

double MLPModel::computeLoss(const Tensor& prediction, const Tensor& target) {
    double loss = 0.0;
    for (size_t i = 0; i < prediction.data.size(); ++i) {
        double diff = prediction.data[i] - target.data[i];
        loss += diff * diff;
    }
    return loss / prediction.data.size();
}

Tensor MLPModel::computeLossGrad(const Tensor& prediction, const Tensor& target) {
    Tensor grad(prediction.shape);
    double scale = 2.0 / prediction.data.size();
    for (size_t i = 0; i < prediction.data.size(); ++i) {
        grad.data[i] = scale * (prediction.data[i] - target.data[i]);
    }
    return grad;
}

std::vector<Tensor*> MLPModel::parameters() {
    std::vector<Tensor*> params;
    
    for (auto& layer : dense_layers_) {
        params.push_back(&layer->weight);
        if (layer->has_bias) params.push_back(&layer->bias_vec);
    }
    
    for (auto& bn : bn_layers_) {
        params.push_back(&bn->gamma);
        params.push_back(&bn->beta);
    }
    
    return params;
}

std::vector<Tensor*> MLPModel::gradients() {
    std::vector<Tensor*> grads;
    
    for (auto& layer : dense_layers_) {
        grads.push_back(&layer->grad_weight);
        if (layer->has_bias) grads.push_back(&layer->grad_bias);
    }
    
    for (auto& bn : bn_layers_) {
        grads.push_back(&bn->grad_gamma);
        grads.push_back(&bn->grad_beta);
    }
    
    return grads;
}

void MLPModel::zeroGrad() {
    for (auto* g : gradients()) {
        g->zeros();
    }
}

size_t MLPModel::numParameters() const {
    size_t count = 0;
    for (const auto& layer : dense_layers_) {
        count += layer->weight.numel();
        if (layer->has_bias) count += layer->bias_vec.numel();
    }
    for (const auto& bn : bn_layers_) {
        count += bn->gamma.numel() + bn->beta.numel();
    }
    return count;
}

void MLPModel::save(const std::string& path) const {
    std::ofstream file(path, std::ios::binary);
    
    // Save config
    int num_layers = config_.layer_sizes.size();
    file.write(reinterpret_cast<const char*>(&num_layers), sizeof(int));
    file.write(reinterpret_cast<const char*>(config_.layer_sizes.data()), 
               num_layers * sizeof(int));
    
    // Save weights
    for (const auto& layer : dense_layers_) {
        layer->weight.save(path + ".dense" + std::to_string(&layer - &dense_layers_[0]));
    }
}

void MLPModel::load(const std::string& path) {
    for (size_t i = 0; i < dense_layers_.size(); ++i) {
        dense_layers_[i]->weight.load(path + ".dense" + std::to_string(i));
    }
}

void MLPModel::summary() const {
    std::cout << "MLP Model Summary\n";
    std::cout << "================\n";
    std::cout << "Layers: ";
    for (int s : config_.layer_sizes) std::cout << s << " ";
    std::cout << "\nActivation: " << config_.activation << "\n";
    std::cout << "Batch Norm: " << (config_.use_batch_norm ? "Yes" : "No") << "\n";
    std::cout << "Dropout: " << (config_.use_dropout ? std::to_string(config_.dropout_rate) : "No") << "\n";
    std::cout << "Parameters: " << numParameters() << "\n";
}

// =============================================================================
// Neural PVT Model Implementation
// =============================================================================

NeuralPVTModel::NeuralPVTModel(const NeuralPVTConfig& config) : config_(config) {
    // Build MLP config
    MLPConfig mlp_config;
    mlp_config.layer_sizes.push_back(3);  // p, T, Rs
    for (int h : config.hidden_layers) {
        mlp_config.layer_sizes.push_back(h);
    }
    mlp_config.layer_sizes.push_back(config.properties.size());
    mlp_config.activation = config.activation;
    mlp_config.use_batch_norm = config.use_batch_norm;
    
    model_ = std::make_unique<MLPModel>(mlp_config);
    
    // Initialize normalization
    input_mean_ = Tensor({3}, 0.0);
    input_std_ = Tensor({3}, 1.0);
    output_mean_ = Tensor({(int)config.properties.size()}, 0.0);
    output_std_ = Tensor({(int)config.properties.size()}, 1.0);
}

void NeuralPVTModel::predict(double pressure, double temperature, double Rs,
                             std::vector<double>& properties) {
    Tensor input = normalizeInput(pressure, temperature, Rs);
    Tensor output = model_->predict(input);
    properties = denormalizeOutput(output);
    
    if (config_.enforce_positivity) {
        enforceConstraints(properties);
    }
}

Tensor NeuralPVTModel::normalizeInput(double p, double T, double Rs) {
    Tensor input({1, 3});
    input(0, 0) = (p - config_.p_min) / (config_.p_max - config_.p_min);
    input(0, 1) = (T - config_.T_min) / (config_.T_max - config_.T_min);
    input(0, 2) = (Rs - config_.Rs_min) / (config_.Rs_max - config_.Rs_min);
    return input;
}

std::vector<double> NeuralPVTModel::denormalizeOutput(const Tensor& output) {
    std::vector<double> result(output.data.size());
    for (size_t i = 0; i < result.size(); ++i) {
        result[i] = output.data[i] * output_std_.data[i % output_std_.numel()] + 
                   output_mean_.data[i % output_mean_.numel()];
    }
    return result;
}

void NeuralPVTModel::enforceConstraints(std::vector<double>& properties) {
    // Ensure positivity
    for (double& p : properties) {
        p = std::max(p, 1e-10);
    }
    
    // FVF >= 1.0
    if (!config_.properties.empty()) {
        for (size_t i = 0; i < config_.properties.size(); ++i) {
            if (config_.properties[i] == PVTProperty::OIL_FVF ||
                config_.properties[i] == PVTProperty::GAS_FVF ||
                config_.properties[i] == PVTProperty::WATER_FVF) {
                properties[i] = std::max(properties[i], 1.0);
            }
        }
    }
}

void NeuralPVTModel::predictBatch(const std::vector<double>& pressures,
                                  const std::vector<double>& temperatures,
                                  const std::vector<double>& Rs_values,
                                  std::vector<std::vector<double>>& properties) {
    int n = pressures.size();
    Tensor input({n, 3});
    
    for (int i = 0; i < n; ++i) {
        input(i, 0) = (pressures[i] - config_.p_min) / (config_.p_max - config_.p_min);
        input(i, 1) = (temperatures[i] - config_.T_min) / (config_.T_max - config_.T_min);
        input(i, 2) = (Rs_values[i] - config_.Rs_min) / (config_.Rs_max - config_.Rs_min);
    }
    
    Tensor output = model_->predict(input);
    
    properties.resize(n);
    int n_props = config_.properties.size();
    for (int i = 0; i < n; ++i) {
        properties[i].resize(n_props);
        for (int j = 0; j < n_props; ++j) {
            properties[i][j] = output.data[i * n_props + j] * output_std_.data[j] + 
                              output_mean_.data[j];
        }
        if (config_.enforce_positivity) {
            enforceConstraints(properties[i]);
        }
    }
}

double NeuralPVTModel::getOilFVF(double p, double T, double Rs) {
    std::vector<double> props;
    predict(p, T, Rs, props);
    for (size_t i = 0; i < config_.properties.size(); ++i) {
        if (config_.properties[i] == PVTProperty::OIL_FVF) {
            return props[i];
        }
    }
    return 1.0;
}

double NeuralPVTModel::getOilViscosity(double p, double T, double Rs) {
    std::vector<double> props;
    predict(p, T, Rs, props);
    for (size_t i = 0; i < config_.properties.size(); ++i) {
        if (config_.properties[i] == PVTProperty::OIL_VISCOSITY) {
            return props[i];
        }
    }
    return 1e-3;
}

double NeuralPVTModel::getGasFVF(double p, double T) {
    std::vector<double> props;
    predict(p, T, 0.0, props);
    for (size_t i = 0; i < config_.properties.size(); ++i) {
        if (config_.properties[i] == PVTProperty::GAS_FVF) {
            return props[i];
        }
    }
    return 1.0;
}

double NeuralPVTModel::getGasViscosity(double p, double T) {
    std::vector<double> props;
    predict(p, T, 0.0, props);
    for (size_t i = 0; i < config_.properties.size(); ++i) {
        if (config_.properties[i] == PVTProperty::GAS_VISCOSITY) {
            return props[i];
        }
    }
    return 1e-5;
}

double NeuralPVTModel::getSolutionGOR(double p, double T) {
    std::vector<double> props;
    predict(p, T, 0.0, props);
    for (size_t i = 0; i < config_.properties.size(); ++i) {
        if (config_.properties[i] == PVTProperty::SOLUTION_GOR) {
            return props[i];
        }
    }
    return 0.0;
}

double NeuralPVTModel::getBubblePoint(double T, double Rs) {
    std::vector<double> props;
    predict(1e7, T, Rs, props);
    for (size_t i = 0; i < config_.properties.size(); ++i) {
        if (config_.properties[i] == PVTProperty::BUBBLE_POINT) {
            return props[i];
        }
    }
    return 1e7;
}

void NeuralPVTModel::train(const std::vector<std::vector<double>>& inputs,
                           const std::vector<std::vector<double>>& outputs,
                           int epochs, double learning_rate) {
    // Compute normalization statistics
    int n = inputs.size();
    int in_dim = inputs[0].size();
    int out_dim = outputs[0].size();
    
    // Input normalization
    for (int d = 0; d < in_dim; ++d) {
        double sum = 0.0, sum_sq = 0.0;
        for (int i = 0; i < n; ++i) {
            sum += inputs[i][d];
            sum_sq += inputs[i][d] * inputs[i][d];
        }
        input_mean_.data[d] = sum / n;
        input_std_.data[d] = std::sqrt(sum_sq / n - input_mean_.data[d] * input_mean_.data[d]) + 1e-8;
    }
    
    // Output normalization
    output_mean_ = Tensor({out_dim}, 0.0);
    output_std_ = Tensor({out_dim}, 1.0);
    for (int d = 0; d < out_dim; ++d) {
        double sum = 0.0, sum_sq = 0.0;
        for (int i = 0; i < n; ++i) {
            sum += outputs[i][d];
            sum_sq += outputs[i][d] * outputs[i][d];
        }
        output_mean_.data[d] = sum / n;
        output_std_.data[d] = std::sqrt(sum_sq / n - output_mean_.data[d] * output_mean_.data[d]) + 1e-8;
    }
    
    // Create optimizer
    AdamOptimizer optimizer(model_->parameters(), learning_rate);
    
    // Training loop
    std::cout << "Training Neural PVT Model...\n";
    
    for (int epoch = 0; epoch < epochs; ++epoch) {
        double total_loss = 0.0;
        
        // Shuffle indices
        std::vector<int> indices(n);
        std::iota(indices.begin(), indices.end(), 0);
        std::shuffle(indices.begin(), indices.end(), std::mt19937(epoch));
        
        // Mini-batch training
        int batch_size = config_.batch_size;
        for (int b = 0; b < n; b += batch_size) {
            int end = std::min(b + batch_size, n);
            int bs = end - b;
            
            // Build batch tensors
            Tensor batch_input({bs, in_dim});
            Tensor batch_target({bs, out_dim});
            
            for (int i = 0; i < bs; ++i) {
                int idx = indices[b + i];
                for (int d = 0; d < in_dim; ++d) {
                    batch_input.data[i * in_dim + d] = 
                        (inputs[idx][d] - input_mean_.data[d]) / input_std_.data[d];
                }
                for (int d = 0; d < out_dim; ++d) {
                    batch_target.data[i * out_dim + d] = 
                        (outputs[idx][d] - output_mean_.data[d]) / output_std_.data[d];
                }
            }
            
            // Forward
            model_->zeroGrad();
            Tensor pred = model_->forward(batch_input);
            
            // Loss
            double loss = model_->computeLoss(pred, batch_target);
            total_loss += loss * bs;
            
            // Backward
            Tensor grad = model_->computeLossGrad(pred, batch_target);
            model_->backward(grad);
            
            // Update
            optimizer.step(model_->gradients());
        }
        
        if ((epoch + 1) % 10 == 0) {
            std::cout << "Epoch " << epoch + 1 << "/" << epochs 
                      << " - Loss: " << total_loss / n << std::endl;
        }
    }
    
    trained_ = true;
}

void NeuralPVTModel::generateTrainingData(
    std::function<void(double, double, double, std::vector<double>&)> pvt_func,
    int num_samples) {
    
    std::vector<std::vector<double>> inputs, outputs;
    
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> p_dist(config_.p_min, config_.p_max);
    std::uniform_real_distribution<double> T_dist(config_.T_min, config_.T_max);
    std::uniform_real_distribution<double> Rs_dist(config_.Rs_min, config_.Rs_max);
    
    for (int i = 0; i < num_samples; ++i) {
        double p = p_dist(rng);
        double T = T_dist(rng);
        double Rs = Rs_dist(rng);
        
        std::vector<double> props;
        pvt_func(p, T, Rs, props);
        
        inputs.push_back({p, T, Rs});
        outputs.push_back(props);
    }
    
    train(inputs, outputs, config_.epochs, config_.learning_rate);
}

void NeuralPVTModel::save(const std::string& path) const {
    model_->save(path);
    input_mean_.save(path + ".input_mean");
    input_std_.save(path + ".input_std");
    output_mean_.save(path + ".output_mean");
    output_std_.save(path + ".output_std");
}

void NeuralPVTModel::load(const std::string& path) {
    model_->load(path);
    input_mean_.load(path + ".input_mean");
    input_std_.load(path + ".input_std");
    output_mean_.load(path + ".output_mean");
    output_std_.load(path + ".output_std");
    trained_ = true;
}

void NeuralPVTModel::exportONNX(const std::string& path) const {
    // ONNX export would require additional library
    std::cerr << "ONNX export not yet implemented\n";
}

double NeuralPVTModel::validate(const std::vector<std::vector<double>>& inputs,
                                const std::vector<std::vector<double>>& outputs) {
    double total_error = 0.0;
    int n = inputs.size();
    
    for (int i = 0; i < n; ++i) {
        std::vector<double> pred;
        predict(inputs[i][0], inputs[i][1], inputs[i][2], pred);
        
        for (size_t j = 0; j < pred.size() && j < outputs[i].size(); ++j) {
            double rel_error = std::abs(pred[j] - outputs[i][j]) / 
                              (std::abs(outputs[i][j]) + 1e-10);
            total_error += rel_error;
        }
    }
    
    return total_error / (n * outputs[0].size());
}

size_t NeuralPVTModel::numParameters() const {
    return model_->numParameters();
}

std::string NeuralPVTModel::getModelInfo() const {
    return "Neural PVT Model - " + std::to_string(numParameters()) + " parameters";
}

void NeuralPVTModel::toGPU() { on_gpu_ = true; }
void NeuralPVTModel::toCPU() { on_gpu_ = false; }

// =============================================================================
// Neural Flash Calculator Implementation
// =============================================================================

NeuralFlashCalculator::NeuralFlashCalculator(const NeuralFlashConfig& config) 
    : config_(config) {
    
    int nc = config.num_components;
    
    // Phase split model: z, P, T -> L
    MLPConfig split_config;
    split_config.layer_sizes = {nc + 2};
    for (int h : config.hidden_layers) split_config.layer_sizes.push_back(h);
    split_config.layer_sizes.push_back(1);
    split_config.activation = config.activation;
    split_config.use_batch_norm = true;
    phase_split_model_ = std::make_unique<MLPModel>(split_config);
    
    // Composition model: z, P, T, L -> x, y
    MLPConfig comp_config;
    comp_config.layer_sizes = {nc + 3};
    for (int h : config.hidden_layers) comp_config.layer_sizes.push_back(h);
    comp_config.layer_sizes.push_back(2 * nc);
    comp_config.activation = config.activation;
    comp_config.use_batch_norm = true;
    composition_model_ = std::make_unique<MLPModel>(comp_config);
    
    // Stability model: z, P, T -> is_two_phase
    MLPConfig stab_config;
    stab_config.layer_sizes = {nc + 2, 128, 64, 1};
    stab_config.activation = "gelu";
    stability_model_ = std::make_unique<MLPModel>(stab_config);
    
    // Initialize normalization
    z_mean_ = Tensor({nc}, 1.0 / nc);
    z_std_ = Tensor({nc}, 0.2);
    P_mean_ = 1e7; P_std_ = 1e7;
    T_mean_ = 350.0; T_std_ = 50.0;
}

void NeuralFlashCalculator::flash(const std::vector<double>& z, double P, double T,
                                   double& L, std::vector<double>& x, 
                                   std::vector<double>& y) {
    int nc = config_.num_components;
    
    // Normalize input
    Tensor split_input({1, nc + 2});
    for (int i = 0; i < nc; ++i) {
        split_input.data[i] = (z[i] - z_mean_.data[i]) / z_std_.data[i];
    }
    split_input.data[nc] = (P - P_mean_) / P_std_;
    split_input.data[nc + 1] = (T - T_mean_) / T_std_;
    
    // Predict L
    Tensor L_pred = phase_split_model_->predict(split_input);
    L = 1.0 / (1.0 + std::exp(-L_pred.data[0]));  // Sigmoid for [0,1]
    
    // Prepare composition input
    Tensor comp_input({1, nc + 3});
    for (int i = 0; i < nc + 2; ++i) {
        comp_input.data[i] = split_input.data[i];
    }
    comp_input.data[nc + 2] = L;
    
    // Predict x, y
    Tensor comp_pred = composition_model_->predict(comp_input);
    
    x.resize(nc);
    y.resize(nc);
    
    // Softmax normalization to ensure sum = 1
    double sum_x = 0.0, sum_y = 0.0;
    for (int i = 0; i < nc; ++i) {
        x[i] = std::exp(comp_pred.data[i]);
        y[i] = std::exp(comp_pred.data[nc + i]);
        sum_x += x[i];
        sum_y += y[i];
    }
    
    for (int i = 0; i < nc; ++i) {
        x[i] /= sum_x;
        y[i] /= sum_y;
    }
    
    if (config_.enforce_mass_balance) {
        enforceConstraints(L, x, y, z);
    }
}

void NeuralFlashCalculator::flashWithStability(const std::vector<double>& z, 
                                                double P, double T,
                                                double& L, std::vector<double>& x, 
                                                std::vector<double>& y,
                                                bool& is_two_phase) {
    int nc = config_.num_components;
    
    // Check stability first
    Tensor stab_input({1, nc + 2});
    for (int i = 0; i < nc; ++i) {
        stab_input.data[i] = (z[i] - z_mean_.data[i]) / z_std_.data[i];
    }
    stab_input.data[nc] = (P - P_mean_) / P_std_;
    stab_input.data[nc + 1] = (T - T_mean_) / T_std_;
    
    Tensor stab_pred = stability_model_->predict(stab_input);
    is_two_phase = stab_pred.data[0] > 0.0;  // Sigmoid > 0.5
    
    if (is_two_phase) {
        flash(z, P, T, L, x, y);
    } else {
        // Single phase - copy z to appropriate phase
        L = 1.0;
        x = z;
        y = z;
    }
}

void NeuralFlashCalculator::flashBatch(const std::vector<std::vector<double>>& z_batch,
                                        const std::vector<double>& P_batch,
                                        const std::vector<double>& T_batch,
                                        std::vector<double>& L_batch,
                                        std::vector<std::vector<double>>& x_batch,
                                        std::vector<std::vector<double>>& y_batch) {
    int n = z_batch.size();
    int nc = config_.num_components;
    
    // Build batch tensor
    Tensor batch_input({n, nc + 2});
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < nc; ++j) {
            batch_input.data[i * (nc + 2) + j] = 
                (z_batch[i][j] - z_mean_.data[j]) / z_std_.data[j];
        }
        batch_input.data[i * (nc + 2) + nc] = (P_batch[i] - P_mean_) / P_std_;
        batch_input.data[i * (nc + 2) + nc + 1] = (T_batch[i] - T_mean_) / T_std_;
    }
    
    // Predict L
    Tensor L_pred = phase_split_model_->predict(batch_input);
    L_batch.resize(n);
    for (int i = 0; i < n; ++i) {
        L_batch[i] = 1.0 / (1.0 + std::exp(-L_pred.data[i]));
    }
    
    // Build composition input
    Tensor comp_input({n, nc + 3});
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < nc + 2; ++j) {
            comp_input.data[i * (nc + 3) + j] = batch_input.data[i * (nc + 2) + j];
        }
        comp_input.data[i * (nc + 3) + nc + 2] = L_batch[i];
    }
    
    // Predict x, y
    Tensor comp_pred = composition_model_->predict(comp_input);
    
    x_batch.resize(n);
    y_batch.resize(n);
    
    for (int i = 0; i < n; ++i) {
        x_batch[i].resize(nc);
        y_batch[i].resize(nc);
        
        double sum_x = 0.0, sum_y = 0.0;
        for (int j = 0; j < nc; ++j) {
            x_batch[i][j] = std::exp(comp_pred.data[i * 2 * nc + j]);
            y_batch[i][j] = std::exp(comp_pred.data[i * 2 * nc + nc + j]);
            sum_x += x_batch[i][j];
            sum_y += y_batch[i][j];
        }
        
        for (int j = 0; j < nc; ++j) {
            x_batch[i][j] /= sum_x;
            y_batch[i][j] /= sum_y;
        }
        
        if (config_.enforce_mass_balance) {
            enforceConstraints(L_batch[i], x_batch[i], y_batch[i], z_batch[i]);
        }
    }
}

void NeuralFlashCalculator::enforceConstraints(double& L, std::vector<double>& x,
                                                std::vector<double>& y,
                                                const std::vector<double>& z) {
    int nc = z.size();
    
    // Clamp L
    L = std::max(0.0, std::min(1.0, L));
    
    // Enforce mass balance: z_i = L * x_i + (1-L) * y_i
    // Simple correction
    for (int i = 0; i < nc; ++i) {
        double calc_z = L * x[i] + (1.0 - L) * y[i];
        double error = z[i] - calc_z;
        // Distribute error
        x[i] += error * 0.5;
        y[i] += error * 0.5;
        x[i] = std::max(0.0, std::min(1.0, x[i]));
        y[i] = std::max(0.0, std::min(1.0, y[i]));
    }
    
    // Renormalize
    double sum_x = 0.0, sum_y = 0.0;
    for (int i = 0; i < nc; ++i) {
        sum_x += x[i];
        sum_y += y[i];
    }
    for (int i = 0; i < nc; ++i) {
        x[i] /= sum_x;
        y[i] /= sum_y;
    }
}

void NeuralFlashCalculator::train(const std::vector<std::vector<double>>& inputs,
                                   const std::vector<std::vector<double>>& outputs,
                                   int epochs, double learning_rate) {
    // Split inputs: z, P, T
    // Outputs: L, x, y
    
    int n = inputs.size();
    int nc = config_.num_components;
    
    // Compute normalization stats
    z_mean_ = Tensor({nc}, 0.0);
    z_std_ = Tensor({nc}, 1.0);
    P_mean_ = 0.0; P_std_ = 1.0;
    T_mean_ = 0.0; T_std_ = 1.0;
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < nc; ++j) {
            z_mean_.data[j] += inputs[i][j];
        }
        P_mean_ += inputs[i][nc];
        T_mean_ += inputs[i][nc + 1];
    }
    
    for (int j = 0; j < nc; ++j) z_mean_.data[j] /= n;
    P_mean_ /= n;
    T_mean_ /= n;
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < nc; ++j) {
            double d = inputs[i][j] - z_mean_.data[j];
            z_std_.data[j] += d * d;
        }
        P_std_ += (inputs[i][nc] - P_mean_) * (inputs[i][nc] - P_mean_);
        T_std_ += (inputs[i][nc + 1] - T_mean_) * (inputs[i][nc + 1] - T_mean_);
    }
    
    for (int j = 0; j < nc; ++j) {
        z_std_.data[j] = std::sqrt(z_std_.data[j] / n) + 1e-8;
    }
    P_std_ = std::sqrt(P_std_ / n) + 1e-8;
    T_std_ = std::sqrt(T_std_ / n) + 1e-8;
    
    // Train phase split model
    std::cout << "Training phase split model...\n";
    AdamOptimizer split_opt(phase_split_model_->parameters(), learning_rate);
    
    for (int epoch = 0; epoch < epochs; ++epoch) {
        double loss = 0.0;
        
        for (int i = 0; i < n; i += config_.batch_size) {
            int bs = std::min(config_.batch_size, n - i);
            
            Tensor batch_in({bs, nc + 2});
            Tensor batch_out({bs, 1});
            
            for (int b = 0; b < bs; ++b) {
                for (int j = 0; j < nc; ++j) {
                    batch_in.data[b * (nc + 2) + j] = 
                        (inputs[i + b][j] - z_mean_.data[j]) / z_std_.data[j];
                }
                batch_in.data[b * (nc + 2) + nc] = (inputs[i + b][nc] - P_mean_) / P_std_;
                batch_in.data[b * (nc + 2) + nc + 1] = (inputs[i + b][nc + 1] - T_mean_) / T_std_;
                
                // L is first output
                double L = outputs[i + b][0];
                batch_out.data[b] = std::log(L / (1.0 - L + 1e-8));  // Inverse sigmoid
            }
            
            phase_split_model_->zeroGrad();
            Tensor pred = phase_split_model_->forward(batch_in);
            loss += phase_split_model_->computeLoss(pred, batch_out) * bs;
            Tensor grad = phase_split_model_->computeLossGrad(pred, batch_out);
            phase_split_model_->backward(grad);
            split_opt.step(phase_split_model_->gradients());
        }
        
        if ((epoch + 1) % 50 == 0) {
            std::cout << "Phase split epoch " << epoch + 1 << " loss: " << loss / n << "\n";
        }
    }
    
    // Train composition model
    std::cout << "Training composition model...\n";
    AdamOptimizer comp_opt(composition_model_->parameters(), learning_rate);
    
    for (int epoch = 0; epoch < epochs; ++epoch) {
        double loss = 0.0;
        
        for (int i = 0; i < n; i += config_.batch_size) {
            int bs = std::min(config_.batch_size, n - i);
            
            Tensor batch_in({bs, nc + 3});
            Tensor batch_out({bs, 2 * nc});
            
            for (int b = 0; b < bs; ++b) {
                for (int j = 0; j < nc; ++j) {
                    batch_in.data[b * (nc + 3) + j] = 
                        (inputs[i + b][j] - z_mean_.data[j]) / z_std_.data[j];
                }
                batch_in.data[b * (nc + 3) + nc] = (inputs[i + b][nc] - P_mean_) / P_std_;
                batch_in.data[b * (nc + 3) + nc + 1] = (inputs[i + b][nc + 1] - T_mean_) / T_std_;
                batch_in.data[b * (nc + 3) + nc + 2] = outputs[i + b][0];  // L
                
                // x and y as log-softmax target
                for (int j = 0; j < nc; ++j) {
                    batch_out.data[b * 2 * nc + j] = std::log(outputs[i + b][1 + j] + 1e-10);
                    batch_out.data[b * 2 * nc + nc + j] = std::log(outputs[i + b][1 + nc + j] + 1e-10);
                }
            }
            
            composition_model_->zeroGrad();
            Tensor pred = composition_model_->forward(batch_in);
            loss += composition_model_->computeLoss(pred, batch_out) * bs;
            Tensor grad = composition_model_->computeLossGrad(pred, batch_out);
            composition_model_->backward(grad);
            comp_opt.step(composition_model_->gradients());
        }
        
        if ((epoch + 1) % 50 == 0) {
            std::cout << "Composition epoch " << epoch + 1 << " loss: " << loss / n << "\n";
        }
    }
    
    trained_ = true;
}

void NeuralFlashCalculator::save(const std::string& path) const {
    phase_split_model_->save(path + ".split");
    composition_model_->save(path + ".comp");
    stability_model_->save(path + ".stab");
    
    std::ofstream f(path + ".norm", std::ios::binary);
    f.write(reinterpret_cast<const char*>(&P_mean_), sizeof(double));
    f.write(reinterpret_cast<const char*>(&P_std_), sizeof(double));
    f.write(reinterpret_cast<const char*>(&T_mean_), sizeof(double));
    f.write(reinterpret_cast<const char*>(&T_std_), sizeof(double));
    z_mean_.save(path + ".z_mean");
    z_std_.save(path + ".z_std");
}

void NeuralFlashCalculator::load(const std::string& path) {
    phase_split_model_->load(path + ".split");
    composition_model_->load(path + ".comp");
    stability_model_->load(path + ".stab");
    
    std::ifstream f(path + ".norm", std::ios::binary);
    f.read(reinterpret_cast<char*>(&P_mean_), sizeof(double));
    f.read(reinterpret_cast<char*>(&P_std_), sizeof(double));
    f.read(reinterpret_cast<char*>(&T_mean_), sizeof(double));
    f.read(reinterpret_cast<char*>(&T_std_), sizeof(double));
    z_mean_.load(path + ".z_mean");
    z_std_.load(path + ".z_std");
    
    trained_ = true;
}

void NeuralFlashCalculator::exportONNX(const std::string& path) const {
    std::cerr << "ONNX export not yet implemented\n";
}

double NeuralFlashCalculator::validate(const std::vector<std::vector<double>>& inputs,
                                        const std::vector<std::vector<double>>& outputs) {
    double total_error = 0.0;
    int n = inputs.size();
    int nc = config_.num_components;
    
    for (int i = 0; i < n; ++i) {
        std::vector<double> z(inputs[i].begin(), inputs[i].begin() + nc);
        double P = inputs[i][nc];
        double T = inputs[i][nc + 1];
        
        double L;
        std::vector<double> x, y;
        flash(z, P, T, L, x, y);
        
        // L error
        total_error += std::abs(L - outputs[i][0]);
        
        // Composition error
        for (int j = 0; j < nc; ++j) {
            total_error += std::abs(x[j] - outputs[i][1 + j]);
            total_error += std::abs(y[j] - outputs[i][1 + nc + j]);
        }
    }
    
    return total_error / (n * (1 + 2 * nc));
}

size_t NeuralFlashCalculator::numParameters() const {
    return phase_split_model_->numParameters() + 
           composition_model_->numParameters() +
           stability_model_->numParameters();
}

std::string NeuralFlashCalculator::getModelInfo() const {
    return "Neural Flash Calculator - " + std::to_string(numParameters()) + " parameters";
}

void NeuralFlashCalculator::toGPU() { on_gpu_ = true; }
void NeuralFlashCalculator::toCPU() { on_gpu_ = false; }

// =============================================================================
// Neural Relative Permeability Implementation
// =============================================================================

NeuralRelPermModel::NeuralRelPermModel(const NeuralRelPermConfig& config) 
    : config_(config) {
    
    MLPConfig mlp_config;
    mlp_config.layer_sizes = {2};  // Sw, Sg
    for (int h : config.hidden_layers) mlp_config.layer_sizes.push_back(h);
    mlp_config.layer_sizes.push_back(1);
    mlp_config.activation = config.activation;
    
    krw_model_ = std::make_unique<MLPModel>(mlp_config);
    kro_model_ = std::make_unique<MLPModel>(mlp_config);
    krg_model_ = std::make_unique<MLPModel>(mlp_config);
    pcow_model_ = std::make_unique<MLPModel>(mlp_config);
    pcgo_model_ = std::make_unique<MLPModel>(mlp_config);
}

void NeuralRelPermModel::getRelPerm(double Sw, double Sg,
                                     double& krw, double& kro, double& krg) {
    Tensor input({1, 2});
    input(0, 0) = Sw;
    input(0, 1) = Sg;
    
    krw = std::abs(krw_model_->predict(input).data[0]);
    kro = std::abs(kro_model_->predict(input).data[0]);
    krg = std::abs(krg_model_->predict(input).data[0]);
    
    // Clamp to [0, 1]
    krw = std::max(0.0, std::min(1.0, krw));
    kro = std::max(0.0, std::min(1.0, kro));
    krg = std::max(0.0, std::min(1.0, krg));
}

void NeuralRelPermModel::getRelPermWithDerivatives(
    double Sw, double Sg,
    double& krw, double& kro, double& krg,
    double& dkrw_dSw, double& dkro_dSw, double& dkrg_dSw,
    double& dkrw_dSg, double& dkro_dSg, double& dkrg_dSg) {
    
    // Central differences
    double eps = 1e-6;
    
    double krw_p, kro_p, krg_p, krw_m, kro_m, krg_m;
    
    // Get base values
    getRelPerm(Sw, Sg, krw, kro, krg);
    
    // Derivatives w.r.t. Sw
    getRelPerm(Sw + eps, Sg, krw_p, kro_p, krg_p);
    getRelPerm(Sw - eps, Sg, krw_m, kro_m, krg_m);
    dkrw_dSw = (krw_p - krw_m) / (2 * eps);
    dkro_dSw = (kro_p - kro_m) / (2 * eps);
    dkrg_dSw = (krg_p - krg_m) / (2 * eps);
    
    // Derivatives w.r.t. Sg
    getRelPerm(Sw, Sg + eps, krw_p, kro_p, krg_p);
    getRelPerm(Sw, Sg - eps, krw_m, kro_m, krg_m);
    dkrw_dSg = (krw_p - krw_m) / (2 * eps);
    dkro_dSg = (kro_p - kro_m) / (2 * eps);
    dkrg_dSg = (krg_p - krg_m) / (2 * eps);
}

void NeuralRelPermModel::getRelPermBatch(const std::vector<double>& Sw,
                                          const std::vector<double>& Sg,
                                          std::vector<double>& krw,
                                          std::vector<double>& kro,
                                          std::vector<double>& krg) {
    int n = Sw.size();
    Tensor input({n, 2});
    
    for (int i = 0; i < n; ++i) {
        input.data[i * 2] = Sw[i];
        input.data[i * 2 + 1] = Sg[i];
    }
    
    Tensor krw_pred = krw_model_->predict(input);
    Tensor kro_pred = kro_model_->predict(input);
    Tensor krg_pred = krg_model_->predict(input);
    
    krw.resize(n);
    kro.resize(n);
    krg.resize(n);
    
    for (int i = 0; i < n; ++i) {
        krw[i] = std::max(0.0, std::min(1.0, std::abs(krw_pred.data[i])));
        kro[i] = std::max(0.0, std::min(1.0, std::abs(kro_pred.data[i])));
        krg[i] = std::max(0.0, std::min(1.0, std::abs(krg_pred.data[i])));
    }
}

void NeuralRelPermModel::getCapillaryPressure(double Sw, double Sg,
                                               double& Pcow, double& Pcgo) {
    Tensor input({1, 2});
    input(0, 0) = Sw;
    input(0, 1) = Sg;
    
    Pcow = pcow_model_->predict(input).data[0];
    Pcgo = pcgo_model_->predict(input).data[0];
}

void NeuralRelPermModel::generateFromCorey(double nw, double no, double ng,
                                            double krw_max, double kro_max, 
                                            double krg_max, int num_samples) {
    std::vector<std::vector<double>> inputs;
    std::vector<std::vector<double>> outputs_w, outputs_o, outputs_g;
    
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    
    for (int i = 0; i < num_samples; ++i) {
        double Sw = dist(rng);
        double Sg = dist(rng) * (1.0 - Sw);
        double So = 1.0 - Sw - Sg;
        
        inputs.push_back({Sw, Sg});
        
        // Corey formulas
        double Sw_norm = (Sw - config_.Swc) / (1.0 - config_.Swc - config_.Sor);
        double Sg_norm = (Sg - config_.Sgc) / (1.0 - config_.Swc - config_.Sor - config_.Sgc);
        double So_norm = (So - config_.Sor) / (1.0 - config_.Swc - config_.Sor);
        
        Sw_norm = std::max(0.0, std::min(1.0, Sw_norm));
        Sg_norm = std::max(0.0, std::min(1.0, Sg_norm));
        So_norm = std::max(0.0, std::min(1.0, So_norm));
        
        double krw = krw_max * std::pow(Sw_norm, nw);
        double kro = kro_max * std::pow(So_norm, no);
        double krg = krg_max * std::pow(Sg_norm, ng);
        
        outputs_w.push_back({krw});
        outputs_o.push_back({kro});
        outputs_g.push_back({krg});
    }
    
    // Train each model
    int epochs = config_.epochs;
    double lr = config_.learning_rate;
    
    AdamOptimizer opt_w(krw_model_->parameters(), lr);
    AdamOptimizer opt_o(kro_model_->parameters(), lr);
    AdamOptimizer opt_g(krg_model_->parameters(), lr);
    
    for (int epoch = 0; epoch < epochs; ++epoch) {
        for (size_t i = 0; i < inputs.size(); ++i) {
            Tensor in({1, 2}, inputs[i]);
            
            krw_model_->zeroGrad();
            Tensor pred_w = krw_model_->forward(in);
            Tensor tgt_w({1, 1}, outputs_w[i]);
            Tensor grad_w = krw_model_->computeLossGrad(pred_w, tgt_w);
            krw_model_->backward(grad_w);
            opt_w.step(krw_model_->gradients());
            
            kro_model_->zeroGrad();
            Tensor pred_o = kro_model_->forward(in);
            Tensor tgt_o({1, 1}, outputs_o[i]);
            Tensor grad_o = kro_model_->computeLossGrad(pred_o, tgt_o);
            kro_model_->backward(grad_o);
            opt_o.step(kro_model_->gradients());
            
            krg_model_->zeroGrad();
            Tensor pred_g = krg_model_->forward(in);
            Tensor tgt_g({1, 1}, outputs_g[i]);
            Tensor grad_g = krg_model_->computeLossGrad(pred_g, tgt_g);
            krg_model_->backward(grad_g);
            opt_g.step(krg_model_->gradients());
        }
    }
    
    trained_ = true;
}

void NeuralRelPermModel::train(const std::vector<std::vector<double>>& inputs,
                                const std::vector<std::vector<double>>& outputs,
                                int epochs, double learning_rate) {
    // outputs format: [krw, kro, krg]
    std::vector<std::vector<double>> out_w, out_o, out_g;
    for (const auto& o : outputs) {
        out_w.push_back({o[0]});
        out_o.push_back({o[1]});
        out_g.push_back({o[2]});
    }
    
    // Train each model independently
    AdamOptimizer opt_w(krw_model_->parameters(), learning_rate);
    AdamOptimizer opt_o(kro_model_->parameters(), learning_rate);
    AdamOptimizer opt_g(krg_model_->parameters(), learning_rate);
    
    for (int epoch = 0; epoch < epochs; ++epoch) {
        for (size_t i = 0; i < inputs.size(); ++i) {
            Tensor in({1, 2}, inputs[i]);
            
            krw_model_->zeroGrad();
            Tensor pred_w = krw_model_->forward(in);
            Tensor tgt_w({1, 1}, out_w[i]);
            krw_model_->backward(krw_model_->computeLossGrad(pred_w, tgt_w));
            opt_w.step(krw_model_->gradients());
            
            kro_model_->zeroGrad();
            Tensor pred_o = kro_model_->forward(in);
            Tensor tgt_o({1, 1}, out_o[i]);
            kro_model_->backward(kro_model_->computeLossGrad(pred_o, tgt_o));
            opt_o.step(kro_model_->gradients());
            
            krg_model_->zeroGrad();
            Tensor pred_g = krg_model_->forward(in);
            Tensor tgt_g({1, 1}, out_g[i]);
            krg_model_->backward(krg_model_->computeLossGrad(pred_g, tgt_g));
            opt_g.step(krg_model_->gradients());
        }
    }
    
    trained_ = true;
}

void NeuralRelPermModel::save(const std::string& path) const {
    krw_model_->save(path + ".krw");
    kro_model_->save(path + ".kro");
    krg_model_->save(path + ".krg");
    pcow_model_->save(path + ".pcow");
    pcgo_model_->save(path + ".pcgo");
}

void NeuralRelPermModel::load(const std::string& path) {
    krw_model_->load(path + ".krw");
    kro_model_->load(path + ".kro");
    krg_model_->load(path + ".krg");
    pcow_model_->load(path + ".pcow");
    pcgo_model_->load(path + ".pcgo");
    trained_ = true;
}

void NeuralRelPermModel::exportONNX(const std::string& path) const {
    std::cerr << "ONNX export not yet implemented\n";
}

double NeuralRelPermModel::validate(const std::vector<std::vector<double>>& inputs,
                                     const std::vector<std::vector<double>>& outputs) {
    double error = 0.0;
    for (size_t i = 0; i < inputs.size(); ++i) {
        double krw, kro, krg;
        getRelPerm(inputs[i][0], inputs[i][1], krw, kro, krg);
        error += std::abs(krw - outputs[i][0]) + 
                 std::abs(kro - outputs[i][1]) + 
                 std::abs(krg - outputs[i][2]);
    }
    return error / (inputs.size() * 3);
}

size_t NeuralRelPermModel::numParameters() const {
    return krw_model_->numParameters() + kro_model_->numParameters() + 
           krg_model_->numParameters() + pcow_model_->numParameters() + 
           pcgo_model_->numParameters();
}

std::string NeuralRelPermModel::getModelInfo() const {
    return "Neural RelPerm Model - " + std::to_string(numParameters()) + " parameters";
}

void NeuralRelPermModel::toGPU() { on_gpu_ = true; }
void NeuralRelPermModel::toCPU() { on_gpu_ = false; }

// =============================================================================
// Neural Preconditioner Implementation
// =============================================================================

NeuralPreconditioner::NeuralPreconditioner(const PrecondConfig& config) 
    : config_(config), reference_pc_(nullptr) {
    
    if (config.n_dof > 0) {
        MLPConfig mlp_config;
        mlp_config.layer_sizes = {config.n_dof};
        for (int h : config.hidden_layers) mlp_config.layer_sizes.push_back(h);
        mlp_config.layer_sizes.push_back(config.n_dof);
        mlp_config.activation = config.activation;
        mlp_config.use_batch_norm = false;  // Important for preconditioners
        
        precond_model_ = std::make_unique<MLPModel>(mlp_config);
    }
}

NeuralPreconditioner::~NeuralPreconditioner() {
    if (reference_pc_) {
        PCDestroy(&reference_pc_);
    }
}

PetscErrorCode NeuralPreconditioner::setup(Mat A) {
    PetscErrorCode ierr;
    
    ierr = MatGetSize(A, &n_rows_, &n_cols_);CHKERRQ(ierr);
    
    // Recreate model if dimensions changed
    if (config_.n_dof != n_rows_) {
        config_.n_dof = n_rows_;
        
        MLPConfig mlp_config;
        mlp_config.layer_sizes = {n_rows_};
        for (int h : config_.hidden_layers) mlp_config.layer_sizes.push_back(h);
        mlp_config.layer_sizes.push_back(n_rows_);
        mlp_config.activation = config_.activation;
        
        precond_model_ = std::make_unique<MLPModel>(mlp_config);
    }
    
    // Setup reference preconditioner if using residual strategy
    if (config_.strategy == PrecondConfig::Strategy::RESIDUAL) {
        if (!reference_pc_) {
            ierr = PCCreate(PETSC_COMM_SELF, &reference_pc_);CHKERRQ(ierr);
        }
        ierr = PCSetType(reference_pc_, PCILU);CHKERRQ(ierr);
        ierr = PCSetOperators(reference_pc_, A, A);CHKERRQ(ierr);
        ierr = PCSetUp(reference_pc_);CHKERRQ(ierr);
    }
    
    return 0;
}

PetscErrorCode NeuralPreconditioner::update(Mat A) {
    if (config_.strategy == PrecondConfig::Strategy::RESIDUAL && reference_pc_) {
        PetscErrorCode ierr;
        ierr = PCSetOperators(reference_pc_, A, A);CHKERRQ(ierr);
        ierr = PCSetUp(reference_pc_);CHKERRQ(ierr);
    }
    return 0;
}

PetscErrorCode NeuralPreconditioner::apply(Vec x, Vec y) {
    auto start = std::chrono::high_resolution_clock::now();
    
    PetscErrorCode ierr;
    const PetscScalar* x_arr;
    PetscScalar* y_arr;
    
    ierr = VecGetArrayRead(x, &x_arr);CHKERRQ(ierr);
    ierr = VecGetArray(y, &y_arr);CHKERRQ(ierr);
    
    // Convert to tensor
    Tensor input({1, n_rows_});
    for (int i = 0; i < n_rows_; ++i) {
        input.data[i] = x_arr[i];
    }
    
    Tensor output;
    
    switch (config_.strategy) {
        case PrecondConfig::Strategy::DIRECT:
            output = precond_model_->predict(input);
            break;
            
        case PrecondConfig::Strategy::RESIDUAL: {
            // Apply reference preconditioner first
            Vec temp;
            ierr = VecDuplicate(y, &temp);CHKERRQ(ierr);
            ierr = PCApply(reference_pc_, x, temp);CHKERRQ(ierr);
            
            const PetscScalar* temp_arr;
            ierr = VecGetArrayRead(temp, &temp_arr);CHKERRQ(ierr);
            
            // Neural correction
            Tensor temp_tensor({1, n_rows_});
            for (int i = 0; i < n_rows_; ++i) {
                temp_tensor.data[i] = temp_arr[i];
            }
            
            Tensor correction = precond_model_->predict(input);
            output = Tensor({1, n_rows_});
            for (int i = 0; i < n_rows_; ++i) {
                output.data[i] = temp_tensor.data[i] + correction.data[i];
            }
            
            ierr = VecRestoreArrayRead(temp, &temp_arr);CHKERRQ(ierr);
            ierr = VecDestroy(&temp);CHKERRQ(ierr);
            break;
        }
        
        case PrecondConfig::Strategy::DEFLATION:
            // Project onto deflation subspace, then apply MLP
            output = precond_model_->predict(input);
            break;
            
        case PrecondConfig::Strategy::COARSE_GRID:
            output = precond_model_->predict(input);
            break;
    }
    
    // Copy to output
    for (int i = 0; i < n_rows_; ++i) {
        y_arr[i] = output.data[i];
    }
    
    ierr = VecRestoreArrayRead(x, &x_arr);CHKERRQ(ierr);
    ierr = VecRestoreArray(y, &y_arr);CHKERRQ(ierr);
    
    auto end = std::chrono::high_resolution_clock::now();
    total_time_ += std::chrono::duration<double>(end - start).count();
    total_applies_++;
    
    return 0;
}

Tensor NeuralPreconditioner::apply(const Tensor& residual) {
    return precond_model_->predict(residual);
}

void NeuralPreconditioner::trainFromSamples(const std::vector<Mat>& A_samples,
                                             const std::vector<Vec>& b_samples,
                                             const std::vector<Vec>& x_samples) {
    std::vector<std::vector<double>> inputs, outputs;
    
    for (size_t s = 0; s < A_samples.size(); ++s) {
        const PetscScalar* b_arr;
        const PetscScalar* x_arr;
        
        VecGetArrayRead(b_samples[s], &b_arr);
        VecGetArrayRead(x_samples[s], &x_arr);
        
        PetscInt n;
        VecGetSize(b_samples[s], &n);
        
        std::vector<double> in(n), out(n);
        for (int i = 0; i < n; ++i) {
            in[i] = b_arr[i];
            out[i] = x_arr[i];
        }
        
        inputs.push_back(in);
        outputs.push_back(out);
        
        VecRestoreArrayRead(b_samples[s], &b_arr);
        VecRestoreArrayRead(x_samples[s], &x_arr);
    }
    
    train(inputs, outputs, config_.epochs, config_.learning_rate);
}

void NeuralPreconditioner::train(const std::vector<std::vector<double>>& inputs,
                                  const std::vector<std::vector<double>>& outputs,
                                  int epochs, double learning_rate) {
    if (inputs.empty()) return;
    
    int n = inputs.size();
    int dim = inputs[0].size();
    
    AdamOptimizer optimizer(precond_model_->parameters(), learning_rate);
    
    std::cout << "Training neural preconditioner...\n";
    
    for (int epoch = 0; epoch < epochs; ++epoch) {
        double total_loss = 0.0;
        
        for (int i = 0; i < n; i += config_.batch_size) {
            int bs = std::min(config_.batch_size, n - i);
            
            Tensor batch_in({bs, dim});
            Tensor batch_out({bs, dim});
            
            for (int b = 0; b < bs; ++b) {
                for (int d = 0; d < dim; ++d) {
                    batch_in.data[b * dim + d] = inputs[i + b][d];
                    batch_out.data[b * dim + d] = outputs[i + b][d];
                }
            }
            
            precond_model_->zeroGrad();
            Tensor pred = precond_model_->forward(batch_in);
            double loss = precond_model_->computeLoss(pred, batch_out);
            total_loss += loss * bs;
            
            Tensor grad = precond_model_->computeLossGrad(pred, batch_out);
            precond_model_->backward(grad);
            optimizer.step(precond_model_->gradients());
        }
        
        if ((epoch + 1) % 10 == 0) {
            std::cout << "Epoch " << epoch + 1 << " loss: " << total_loss / n << "\n";
        }
    }
    
    trained_ = true;
}

void NeuralPreconditioner::save(const std::string& path) const {
    precond_model_->save(path);
}

void NeuralPreconditioner::load(const std::string& path) {
    precond_model_->load(path);
    trained_ = true;
}

void NeuralPreconditioner::exportONNX(const std::string& path) const {
    std::cerr << "ONNX export not yet implemented\n";
}

double NeuralPreconditioner::validate(const std::vector<std::vector<double>>& inputs,
                                       const std::vector<std::vector<double>>& outputs) {
    double total_error = 0.0;
    int dim = inputs[0].size();
    
    for (size_t i = 0; i < inputs.size(); ++i) {
        Tensor in({1, dim}, inputs[i]);
        Tensor pred = precond_model_->predict(in);
        
        for (int d = 0; d < dim; ++d) {
            total_error += std::abs(pred.data[d] - outputs[i][d]);
        }
    }
    
    return total_error / (inputs.size() * dim);
}

size_t NeuralPreconditioner::numParameters() const {
    return precond_model_ ? precond_model_->numParameters() : 0;
}

std::string NeuralPreconditioner::getModelInfo() const {
    return "Neural Preconditioner - " + std::to_string(numParameters()) + " parameters";
}

void NeuralPreconditioner::toGPU() { on_gpu_ = true; }
void NeuralPreconditioner::toCPU() { on_gpu_ = false; }

// =============================================================================
// Surrogate Trainer Utilities
// =============================================================================

SurrogateTrainer::SurrogateTrainer() {}

std::vector<std::vector<double>> SurrogateTrainer::latinHypercube(
    const std::vector<double>& lower_bounds,
    const std::vector<double>& upper_bounds,
    int num_samples) {
    
    int dim = lower_bounds.size();
    std::vector<std::vector<double>> samples(num_samples, std::vector<double>(dim));
    
    std::mt19937 rng(42);
    
    for (int d = 0; d < dim; ++d) {
        // Create permutation
        std::vector<int> perm(num_samples);
        std::iota(perm.begin(), perm.end(), 0);
        std::shuffle(perm.begin(), perm.end(), rng);
        
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        
        for (int i = 0; i < num_samples; ++i) {
            double u = (perm[i] + dist(rng)) / num_samples;
            samples[i][d] = lower_bounds[d] + u * (upper_bounds[d] - lower_bounds[d]);
        }
    }
    
    return samples;
}

std::vector<std::vector<double>> SurrogateTrainer::sobolSequence(
    const std::vector<double>& lower_bounds,
    const std::vector<double>& upper_bounds,
    int num_samples) {
    
    int dim = lower_bounds.size();
    std::vector<std::vector<double>> samples(num_samples, std::vector<double>(dim));
    
    // Simplified Sobol - use Van der Corput sequence for each dimension
    for (int i = 0; i < num_samples; ++i) {
        for (int d = 0; d < dim; ++d) {
            // Van der Corput in base 2 + d
            int base = 2 + d;
            double u = 0.0;
            double f = 1.0 / base;
            int n = i + 1;
            
            while (n > 0) {
                u += f * (n % base);
                n /= base;
                f /= base;
            }
            
            samples[i][d] = lower_bounds[d] + u * (upper_bounds[d] - lower_bounds[d]);
        }
    }
    
    return samples;
}

void SurrogateTrainer::augmentWithNoise(
    std::vector<std::vector<double>>& inputs,
    std::vector<std::vector<double>>& outputs,
    double noise_level,
    int augmentation_factor) {
    
    size_t n = inputs.size();
    std::mt19937 rng(42);
    std::normal_distribution<double> noise(0.0, noise_level);
    
    for (int aug = 0; aug < augmentation_factor - 1; ++aug) {
        for (size_t i = 0; i < n; ++i) {
            std::vector<double> noisy_input = inputs[i];
            for (double& x : noisy_input) {
                x *= (1.0 + noise(rng));
            }
            inputs.push_back(noisy_input);
            outputs.push_back(outputs[i]);
        }
    }
}

// =============================================================================
// Configuration Parsing
// =============================================================================

NeuralPVTConfig parseNeuralPVTConfig(const std::map<std::string, std::string>& config) {
    NeuralPVTConfig result;
    
    auto getInt = [&](const std::string& key, int def) {
        auto it = config.find(key);
        return it != config.end() ? std::stoi(it->second) : def;
    };
    
    auto getDouble = [&](const std::string& key, double def) {
        auto it = config.find(key);
        return it != config.end() ? std::stod(it->second) : def;
    };
    
    auto getString = [&](const std::string& key, const std::string& def) {
        auto it = config.find(key);
        return it != config.end() ? it->second : def;
    };
    
    result.activation = getString("neural_pvt_activation", "gelu");
    result.learning_rate = getDouble("neural_pvt_lr", 1e-3);
    result.batch_size = getInt("neural_pvt_batch_size", 256);
    result.epochs = getInt("neural_pvt_epochs", 200);
    result.p_min = getDouble("neural_pvt_p_min", 1e5);
    result.p_max = getDouble("neural_pvt_p_max", 1e8);
    result.T_min = getDouble("neural_pvt_T_min", 273.15);
    result.T_max = getDouble("neural_pvt_T_max", 473.15);
    
    return result;
}

NeuralFlashConfig parseNeuralFlashConfig(const std::map<std::string, std::string>& config) {
    NeuralFlashConfig result;
    
    auto getInt = [&](const std::string& key, int def) {
        auto it = config.find(key);
        return it != config.end() ? std::stoi(it->second) : def;
    };
    
    result.num_components = getInt("neural_flash_num_components", 5);
    result.batch_size = getInt("neural_flash_batch_size", 128);
    result.epochs = getInt("neural_flash_epochs", 500);
    
    return result;
}

NeuralRelPermConfig parseNeuralRelPermConfig(const std::map<std::string, std::string>& config) {
    NeuralRelPermConfig result;
    
    auto getDouble = [&](const std::string& key, double def) {
        auto it = config.find(key);
        return it != config.end() ? std::stod(it->second) : def;
    };
    
    result.Swc = getDouble("neural_relperm_Swc", 0.2);
    result.Sor = getDouble("neural_relperm_Sor", 0.2);
    result.Sgc = getDouble("neural_relperm_Sgc", 0.05);
    
    return result;
}

} // namespace ML
} // namespace FSRM
