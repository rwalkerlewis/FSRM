/**
 * @file NeuralReducedOrderModel.cpp
 * @brief Implementation of Neural Network-Enhanced Reduced Order Models
 */

#include "NeuralReducedOrderModel.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <random>
#include <iostream>
#include <iomanip>
#include <fstream>

namespace FSRM {
namespace ML {

// =============================================================================
// POD Basis Implementation
// =============================================================================

void PODBasis::compute(const std::vector<Tensor>& snapshots, 
                       int num_modes, double energy_threshold) {
    if (snapshots.empty()) {
        throw std::runtime_error("Empty snapshot set for POD");
    }
    
    n_dof_ = snapshots[0].numel();
    n_snapshots_ = snapshots.size();
    
    // Compute mean
    mean_ = Tensor({(int)n_dof_}, 0.0);
    for (const auto& snap : snapshots) {
        for (int i = 0; i < n_dof_; ++i) {
            mean_.data[i] += snap.data[i];
        }
    }
    for (int i = 0; i < n_dof_; ++i) {
        mean_.data[i] /= n_snapshots_;
    }
    
    // Center snapshots
    Tensor data({n_dof_, n_snapshots_});
    for (int j = 0; j < n_snapshots_; ++j) {
        for (int i = 0; i < n_dof_; ++i) {
            data.data[i * n_snapshots_ + j] = snapshots[j].data[i] - mean_.data[i];
        }
    }
    
    // SVD
    Tensor U, S, Vt;
    svd(data, U, S, Vt);
    
    // Determine number of modes based on energy
    double total_energy = 0.0;
    for (const auto& s : S.data) {
        total_energy += s * s;
    }
    
    double cumulative_energy = 0.0;
    int modes_by_energy = S.numel();
    for (int i = 0; i < (int)S.numel(); ++i) {
        cumulative_energy += S.data[i] * S.data[i];
        if (cumulative_energy / total_energy >= energy_threshold) {
            modes_by_energy = i + 1;
            break;
        }
    }
    
    if (num_modes > 0) {
        num_modes_ = std::min(num_modes, modes_by_energy);
    } else {
        num_modes_ = modes_by_energy;
    }
    
    num_modes_ = std::min(num_modes_, (int)S.numel());
    
    // Store basis (first num_modes_ columns of U)
    basis_ = Tensor({n_dof_, num_modes_});
    for (int i = 0; i < n_dof_; ++i) {
        for (int j = 0; j < num_modes_; ++j) {
            basis_.data[i * num_modes_ + j] = U.data[i * n_snapshots_ + j];
        }
    }
    
    // Store singular values
    singular_values_ = Tensor({num_modes_});
    for (int i = 0; i < num_modes_; ++i) {
        singular_values_.data[i] = S.data[i];
    }
    
    // Compute energy captured
    cumulative_energy = 0.0;
    for (int i = 0; i < num_modes_; ++i) {
        cumulative_energy += S.data[i] * S.data[i];
    }
    energy_captured_ = cumulative_energy / total_energy;
    
    std::cout << "POD: Retained " << num_modes_ << " modes capturing " 
              << 100.0 * energy_captured_ << "% energy\n";
}

void PODBasis::compute(const std::vector<Vec>& snapshots,
                       int num_modes, double energy_threshold) {
    std::vector<Tensor> tensor_snapshots;
    
    for (const auto& v : snapshots) {
        PetscInt n;
        VecGetSize(v, &n);
        
        const PetscScalar* arr;
        VecGetArrayRead(v, &arr);
        
        Tensor t({(int)n});
        for (int i = 0; i < n; ++i) {
            t.data[i] = arr[i];
        }
        
        VecRestoreArrayRead(v, &arr);
        tensor_snapshots.push_back(std::move(t));
    }
    
    compute(tensor_snapshots, num_modes, energy_threshold);
}

void PODBasis::svd(const Tensor& data, Tensor& U, Tensor& S, Tensor& Vt) {
    // Method of snapshots for efficient SVD when n_dof >> n_snapshots
    int m = data.shape[0];  // n_dof
    int n = data.shape[1];  // n_snapshots
    
    // Compute correlation matrix C = A^T * A
    Tensor C({n, n}, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double sum = 0.0;
            for (int k = 0; k < m; ++k) {
                sum += data.data[k * n + i] * data.data[k * n + j];
            }
            C.data[i * n + j] = sum;
        }
    }
    
    // Eigendecomposition of C (simplified power iteration)
    Tensor eigenvalues({n});
    Tensor eigenvectors({n, n});
    
    // Initialize eigenvectors to identity
    for (int i = 0; i < n; ++i) {
        eigenvectors.data[i * n + i] = 1.0;
    }
    
    // Simple eigendecomposition (power iteration with deflation)
    for (int mode = 0; mode < n; ++mode) {
        Tensor v({n});
        v.randn(0.0, 1.0);
        
        // Normalize v
        double norm = 0.0;
        for (double val : v.data) norm += val * val;
        norm = std::sqrt(norm);
        for (double& val : v.data) val /= norm;
        
        // Power iteration
        for (int iter = 0; iter < 100; ++iter) {
            Tensor v_new({n}, 0.0);
            
            // Matrix-vector multiply
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    v_new.data[i] += C.data[i * n + j] * v.data[j];
                }
            }
            
            // Deflation
            for (int prev = 0; prev < mode; ++prev) {
                double dot = 0.0;
                for (int i = 0; i < n; ++i) {
                    dot += v_new.data[i] * eigenvectors.data[i * n + prev];
                }
                for (int i = 0; i < n; ++i) {
                    v_new.data[i] -= dot * eigenvectors.data[i * n + prev];
                }
            }
            
            double norm = v_new.norm();
            if (norm < 1e-14) break;
            
            for (int i = 0; i < n; ++i) {
                v.data[i] = v_new.data[i] / norm;
            }
        }
        
        // Compute eigenvalue
        Tensor Cv({n}, 0.0);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                Cv.data[i] += C.data[i * n + j] * v.data[j];
            }
        }
        
        double eigenvalue = 0.0;
        for (int i = 0; i < n; ++i) {
            eigenvalue += v.data[i] * Cv.data[i];
        }
        
        eigenvalues.data[mode] = std::max(0.0, eigenvalue);
        for (int i = 0; i < n; ++i) {
            eigenvectors.data[i * n + mode] = v.data[i];
        }
    }
    
    // Singular values
    S = Tensor({n});
    for (int i = 0; i < n; ++i) {
        S.data[i] = std::sqrt(eigenvalues.data[i]);
    }
    
    // Left singular vectors U = A * V * S^{-1}
    U = Tensor({m, n});
    for (int j = 0; j < n; ++j) {
        if (S.data[j] > 1e-14) {
            for (int i = 0; i < m; ++i) {
                double sum = 0.0;
                for (int k = 0; k < n; ++k) {
                    sum += data.data[i * n + k] * eigenvectors.data[k * n + j];
                }
                U.data[i * n + j] = sum / S.data[j];
            }
        }
    }
    
    // Right singular vectors (V^T) - transpose of eigenvectors
    Vt = Tensor({n, n});
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            Vt.data[i * n + j] = eigenvectors.data[j * n + i];
        }
    }
}

Tensor PODBasis::project(const Tensor& full_state) const {
    Tensor reduced({num_modes_});
    
    // Subtract mean and project
    for (int j = 0; j < num_modes_; ++j) {
        double sum = 0.0;
        for (int i = 0; i < n_dof_; ++i) {
            sum += basis_.data[i * num_modes_ + j] * (full_state.data[i] - mean_.data[i]);
        }
        reduced.data[j] = sum;
    }
    
    return reduced;
}

void PODBasis::project(Vec full_state, Tensor& reduced) const {
    const PetscScalar* arr;
    VecGetArrayRead(full_state, &arr);
    
    reduced = Tensor({num_modes_});
    for (int j = 0; j < num_modes_; ++j) {
        double sum = 0.0;
        for (int i = 0; i < n_dof_; ++i) {
            sum += basis_.data[i * num_modes_ + j] * (arr[i] - mean_.data[i]);
        }
        reduced.data[j] = sum;
    }
    
    VecRestoreArrayRead(full_state, &arr);
}

Tensor PODBasis::reconstruct(const Tensor& reduced_state) const {
    Tensor full({n_dof_});
    
    // Reconstruct: U = Phi * a + mean
    for (int i = 0; i < n_dof_; ++i) {
        double sum = mean_.data[i];
        for (int j = 0; j < num_modes_; ++j) {
            sum += basis_.data[i * num_modes_ + j] * reduced_state.data[j];
        }
        full.data[i] = sum;
    }
    
    return full;
}

void PODBasis::reconstruct(const Tensor& reduced_state, Vec full_state) const {
    PetscScalar* arr;
    VecGetArray(full_state, &arr);
    
    for (int i = 0; i < n_dof_; ++i) {
        double sum = mean_.data[i];
        for (int j = 0; j < num_modes_; ++j) {
            sum += basis_.data[i * num_modes_ + j] * reduced_state.data[j];
        }
        arr[i] = sum;
    }
    
    VecRestoreArray(full_state, &arr);
}

double PODBasis::getProjectionError(const Tensor& full_state) const {
    Tensor reduced = project(full_state);
    Tensor reconstructed = reconstruct(reduced);
    
    double error = 0.0, norm = 0.0;
    for (int i = 0; i < n_dof_; ++i) {
        double diff = full_state.data[i] - reconstructed.data[i];
        error += diff * diff;
        norm += full_state.data[i] * full_state.data[i];
    }
    
    return std::sqrt(error / (norm + 1e-14));
}

void PODBasis::save(const std::string& path) const {
    basis_.save(path + ".basis");
    singular_values_.save(path + ".singular_values");
    mean_.save(path + ".mean");
    
    std::ofstream f(path + ".info");
    f << n_dof_ << " " << num_modes_ << " " << n_snapshots_ << " " << energy_captured_;
}

void PODBasis::load(const std::string& path) {
    basis_.load(path + ".basis");
    singular_values_.load(path + ".singular_values");
    mean_.load(path + ".mean");
    
    std::ifstream f(path + ".info");
    f >> n_dof_ >> num_modes_ >> n_snapshots_ >> energy_captured_;
}

// =============================================================================
// POD-NN ROM Implementation
// =============================================================================

PODNeuralROM::PODNeuralROM(const PODNNConfig& config) : config_(config) {
    // Build dynamics network
    MLPConfig mlp_config;
    mlp_config.layer_sizes.push_back(config.num_modes + config.num_parameters + 1);  // a + params + t
    for (int h : config.hidden_layers) {
        mlp_config.layer_sizes.push_back(h);
    }
    mlp_config.layer_sizes.push_back(config.num_modes);  // da/dt
    mlp_config.activation = config.activation;
    mlp_config.use_batch_norm = config.use_batch_norm;
    
    dynamics_network_ = std::make_unique<MLPModel>(mlp_config);
    
    if (config.parameter_aware && config.num_parameters > 0) {
        MLPConfig param_config;
        param_config.layer_sizes = {config.num_parameters, 64, 32};
        parameter_encoder_ = std::make_unique<MLPModel>(param_config);
    }
}

void PODNeuralROM::buildFromSnapshots(const std::vector<Tensor>& snapshots,
                                       const std::vector<double>& times,
                                       const std::vector<Tensor>& parameters) {
    // Build POD basis
    pod_basis_.compute(snapshots, config_.num_modes, config_.energy_threshold);
    
    // Project snapshots to reduced space
    std::vector<Tensor> reduced_states;
    for (const auto& snap : snapshots) {
        reduced_states.push_back(pod_basis_.project(snap));
    }
    
    // Compute time derivatives (finite differences)
    std::vector<Tensor> reduced_velocities;
    for (size_t i = 0; i < reduced_states.size() - 1; ++i) {
        double dt = times[i + 1] - times[i];
        Tensor velocity({pod_basis_.numModes()});
        for (int j = 0; j < pod_basis_.numModes(); ++j) {
            velocity.data[j] = (reduced_states[i + 1].data[j] - reduced_states[i].data[j]) / dt;
        }
        reduced_velocities.push_back(std::move(velocity));
    }
    
    // Train dynamics network
    train(reduced_states, reduced_velocities, parameters);
}

Tensor PODNeuralROM::solve(const Tensor& initial_condition,
                           double t_final, double dt,
                           const Tensor& parameters) {
    // Project to reduced space
    Tensor a = pod_basis_.project(initial_condition);
    
    // Time stepping
    double t = 0.0;
    while (t < t_final - 1e-10) {
        if (config_.time_integrator == "rk4") {
            a = rk4Step(a, t, dt, parameters);
        } else {
            a = implicitEulerStep(a, dt, parameters);
        }
        t += dt;
    }
    
    // Reconstruct full state
    return pod_basis_.reconstruct(a);
}

std::vector<Tensor> PODNeuralROM::trajectory(const Tensor& initial_condition,
                                              double t_final, double dt,
                                              const Tensor& parameters) {
    std::vector<Tensor> results;
    
    Tensor a = pod_basis_.project(initial_condition);
    results.push_back(pod_basis_.reconstruct(a));
    
    double t = 0.0;
    while (t < t_final - 1e-10) {
        if (config_.time_integrator == "rk4") {
            a = rk4Step(a, t, dt, parameters);
        } else {
            a = implicitEulerStep(a, dt, parameters);
        }
        results.push_back(pod_basis_.reconstruct(a));
        t += dt;
    }
    
    return results;
}

Tensor PODNeuralROM::evaluateReducedDynamics(const Tensor& reduced_state,
                                              double time,
                                              const Tensor& parameters) {
    int num_modes = pod_basis_.numModes();
    int input_size = num_modes + config_.num_parameters + 1;
    
    Tensor input({1, input_size});
    
    for (int i = 0; i < num_modes; ++i) {
        input.data[i] = reduced_state.data[i];
    }
    
    for (int i = 0; i < config_.num_parameters && i < (int)parameters.numel(); ++i) {
        input.data[num_modes + i] = parameters.data[i];
    }
    
    input.data[input_size - 1] = time;
    
    return dynamics_network_->predict(input);
}

void PODNeuralROM::train(const std::vector<Tensor>& reduced_states,
                         const std::vector<Tensor>& reduced_velocities,
                         const std::vector<Tensor>& parameters) {
    int n = reduced_velocities.size();
    int num_modes = pod_basis_.numModes();
    int input_size = num_modes + config_.num_parameters + 1;
    
    std::cout << "Training POD-NN dynamics on " << n << " samples...\n";
    
    AdamOptimizer optimizer(dynamics_network_->parameters(), config_.learning_rate);
    
    double total_loss = 0.0;
    for (int epoch = 0; epoch < config_.epochs; ++epoch) {
        total_loss = 0.0;
        
        // Shuffle
        std::vector<int> indices(n);
        std::iota(indices.begin(), indices.end(), 0);
        std::shuffle(indices.begin(), indices.end(), std::mt19937(epoch));
        
        for (int b = 0; b < n; b += config_.batch_size) {
            int batch_end = std::min(b + config_.batch_size, n);
            int bs = batch_end - b;
            
            Tensor batch_input({bs, input_size});
            Tensor batch_target({bs, num_modes});
            
            for (int i = 0; i < bs; ++i) {
                int idx = indices[b + i];
                
                for (int j = 0; j < num_modes; ++j) {
                    batch_input.data[i * input_size + j] = reduced_states[idx].data[j];
                    batch_target.data[i * num_modes + j] = reduced_velocities[idx].data[j];
                }
                
                // Add parameters if available
                if (!parameters.empty() && idx < (int)parameters.size()) {
                    for (int j = 0; j < config_.num_parameters; ++j) {
                        batch_input.data[i * input_size + num_modes + j] = 
                            parameters[idx].data[j];
                    }
                }
                
                // Time (normalized)
                batch_input.data[i * input_size + input_size - 1] = (double)idx / n;
            }
            
            dynamics_network_->zeroGrad();
            Tensor pred = dynamics_network_->forward(batch_input);
            double loss = dynamics_network_->computeLoss(pred, batch_target);
            total_loss += loss * bs;
            
            Tensor grad = dynamics_network_->computeLossGrad(pred, batch_target);
            dynamics_network_->backward(grad);
            optimizer.step(dynamics_network_->gradients());
        }
        
        if ((epoch + 1) % 20 == 0) {
            std::cout << "Epoch " << epoch + 1 << "/" << config_.epochs
                      << " - loss: " << total_loss / n << "\n";
        }
    }
    
    training_error_ = total_loss / n;
}

Tensor PODNeuralROM::implicitEulerStep(const Tensor& a, double dt,
                                        const Tensor& params) {
    // Simple fixed-point iteration for implicit Euler
    Tensor a_new = a;
    
    for (int iter = 0; iter < 10; ++iter) {
        Tensor f = evaluateReducedDynamics(a_new, 0.0, params);
        
        Tensor a_next({a.numel()});
        for (size_t i = 0; i < a.numel(); ++i) {
            a_next.data[i] = a.data[i] + dt * f.data[i];
        }
        
        double diff = 0.0;
        for (size_t i = 0; i < a.numel(); ++i) {
            diff += std::abs(a_next.data[i] - a_new.data[i]);
        }
        
        a_new = a_next;
        if (diff < 1e-8) break;
    }
    
    return a_new;
}

Tensor PODNeuralROM::rk4Step(const Tensor& a, double t, double dt,
                             const Tensor& params) {
    Tensor k1 = evaluateReducedDynamics(a, t, params);
    
    Tensor a_temp({a.numel()});
    for (size_t i = 0; i < a.numel(); ++i) {
        a_temp.data[i] = a.data[i] + 0.5 * dt * k1.data[i];
    }
    Tensor k2 = evaluateReducedDynamics(a_temp, t + 0.5 * dt, params);
    
    for (size_t i = 0; i < a.numel(); ++i) {
        a_temp.data[i] = a.data[i] + 0.5 * dt * k2.data[i];
    }
    Tensor k3 = evaluateReducedDynamics(a_temp, t + 0.5 * dt, params);
    
    for (size_t i = 0; i < a.numel(); ++i) {
        a_temp.data[i] = a.data[i] + dt * k3.data[i];
    }
    Tensor k4 = evaluateReducedDynamics(a_temp, t + dt, params);
    
    Tensor a_new({a.numel()});
    for (size_t i = 0; i < a.numel(); ++i) {
        a_new.data[i] = a.data[i] + dt / 6.0 * (k1.data[i] + 2*k2.data[i] + 
                                                  2*k3.data[i] + k4.data[i]);
    }
    
    return a_new;
}

double PODNeuralROM::validate(const std::vector<Tensor>& test_snapshots,
                               const std::vector<double>& test_times) {
    if (test_snapshots.empty()) return 0.0;
    
    // Solve with initial condition from first snapshot
    double t_final = test_times.back() - test_times.front();
    double dt = (t_final) / (test_times.size() - 1);
    
    std::vector<Tensor> rom_trajectory = trajectory(test_snapshots[0], t_final, dt);
    
    // Compute relative error
    double total_error = 0.0, total_norm = 0.0;
    
    size_t n = std::min(rom_trajectory.size(), test_snapshots.size());
    for (size_t i = 0; i < n; ++i) {
        double err = 0.0, nrm = 0.0;
        for (size_t j = 0; j < test_snapshots[i].numel(); ++j) {
            double diff = rom_trajectory[i].data[j] - test_snapshots[i].data[j];
            err += diff * diff;
            nrm += test_snapshots[i].data[j] * test_snapshots[i].data[j];
        }
        total_error += std::sqrt(err);
        total_norm += std::sqrt(nrm);
    }
    
    validation_error_ = total_error / (total_norm + 1e-14);
    return validation_error_;
}

void PODNeuralROM::save(const std::string& path) const {
    pod_basis_.save(path + ".pod");
    dynamics_network_->save(path + ".dynamics");
}

void PODNeuralROM::load(const std::string& path) {
    pod_basis_.load(path + ".pod");
    dynamics_network_->load(path + ".dynamics");
}

size_t PODNeuralROM::numParameters() const {
    return dynamics_network_->numParameters();
}

// =============================================================================
// LSTM Cell Implementation
// =============================================================================

LSTMCell::LSTMCell(int input_size, int hidden_size) 
    : input_size_(input_size), hidden_size_(hidden_size) {
    
    // Initialize weights
    auto initWeight = [](int in, int out) {
        Tensor w({in, out});
        w.xavier_init();
        return w;
    };
    
    auto initBias = [](int size) {
        return Tensor({size}, 0.0);
    };
    
    W_f_ = initWeight(input_size, hidden_size);
    U_f_ = initWeight(hidden_size, hidden_size);
    b_f_ = Tensor({hidden_size}, 1.0);  // Bias forget gate to 1
    
    W_i_ = initWeight(input_size, hidden_size);
    U_i_ = initWeight(hidden_size, hidden_size);
    b_i_ = initBias(hidden_size);
    
    W_o_ = initWeight(input_size, hidden_size);
    U_o_ = initWeight(hidden_size, hidden_size);
    b_o_ = initBias(hidden_size);
    
    W_c_ = initWeight(input_size, hidden_size);
    U_c_ = initWeight(hidden_size, hidden_size);
    b_c_ = initBias(hidden_size);
    
    // Initialize gradients
    grad_W_f_ = Tensor({input_size, hidden_size});
    grad_U_f_ = Tensor({hidden_size, hidden_size});
    grad_b_f_ = Tensor({hidden_size});
    
    grad_W_i_ = Tensor({input_size, hidden_size});
    grad_U_i_ = Tensor({hidden_size, hidden_size});
    grad_b_i_ = Tensor({hidden_size});
    
    grad_W_o_ = Tensor({input_size, hidden_size});
    grad_U_o_ = Tensor({hidden_size, hidden_size});
    grad_b_o_ = Tensor({hidden_size});
    
    grad_W_c_ = Tensor({input_size, hidden_size});
    grad_U_c_ = Tensor({hidden_size, hidden_size});
    grad_b_c_ = Tensor({hidden_size});
}

Tensor LSTMCell::sigmoid(const Tensor& x) {
    Tensor result(x.shape);
    for (size_t i = 0; i < x.data.size(); ++i) {
        result.data[i] = 1.0 / (1.0 + std::exp(-x.data[i]));
    }
    return result;
}

Tensor LSTMCell::tanh(const Tensor& x) {
    Tensor result(x.shape);
    for (size_t i = 0; i < x.data.size(); ++i) {
        result.data[i] = std::tanh(x.data[i]);
    }
    return result;
}

std::pair<Tensor, Tensor> LSTMCell::forward(const Tensor& input,
                                             const Tensor& h_prev,
                                             const Tensor& c_prev) {
    cached_input_ = input;
    cached_h_ = h_prev;
    cached_c_ = c_prev;
    
    // Compute gates
    Tensor f_gate({hidden_size_}), i_gate({hidden_size_});
    Tensor o_gate({hidden_size_}), c_hat({hidden_size_});
    
    for (int h = 0; h < hidden_size_; ++h) {
        double f = b_f_.data[h], i = b_i_.data[h];
        double o = b_o_.data[h], c = b_c_.data[h];
        
        for (int x = 0; x < input_size_; ++x) {
            f += W_f_.data[x * hidden_size_ + h] * input.data[x];
            i += W_i_.data[x * hidden_size_ + h] * input.data[x];
            o += W_o_.data[x * hidden_size_ + h] * input.data[x];
            c += W_c_.data[x * hidden_size_ + h] * input.data[x];
        }
        
        for (int hp = 0; hp < hidden_size_; ++hp) {
            f += U_f_.data[hp * hidden_size_ + h] * h_prev.data[hp];
            i += U_i_.data[hp * hidden_size_ + h] * h_prev.data[hp];
            o += U_o_.data[hp * hidden_size_ + h] * h_prev.data[hp];
            c += U_c_.data[hp * hidden_size_ + h] * h_prev.data[hp];
        }
        
        f_gate.data[h] = 1.0 / (1.0 + std::exp(-f));
        i_gate.data[h] = 1.0 / (1.0 + std::exp(-i));
        o_gate.data[h] = 1.0 / (1.0 + std::exp(-o));
        c_hat.data[h] = std::tanh(c);
    }
    
    cached_f_ = f_gate;
    cached_i_ = i_gate;
    cached_o_ = o_gate;
    cached_c_hat_ = c_hat;
    
    // New cell state
    Tensor c_new({hidden_size_});
    for (int h = 0; h < hidden_size_; ++h) {
        c_new.data[h] = f_gate.data[h] * c_prev.data[h] + 
                        i_gate.data[h] * c_hat.data[h];
    }
    
    // New hidden state
    Tensor h_new({hidden_size_});
    for (int h = 0; h < hidden_size_; ++h) {
        h_new.data[h] = o_gate.data[h] * std::tanh(c_new.data[h]);
    }
    
    return {h_new, c_new};
}

std::vector<Tensor*> LSTMCell::parameters() {
    return {&W_f_, &U_f_, &b_f_, &W_i_, &U_i_, &b_i_, 
            &W_o_, &U_o_, &b_o_, &W_c_, &U_c_, &b_c_};
}

std::vector<Tensor*> LSTMCell::gradients() {
    return {&grad_W_f_, &grad_U_f_, &grad_b_f_, &grad_W_i_, &grad_U_i_, &grad_b_i_,
            &grad_W_o_, &grad_U_o_, &grad_b_o_, &grad_W_c_, &grad_U_c_, &grad_b_c_};
}

void LSTMCell::zeroGrad() {
    for (auto* g : gradients()) {
        g->zeros();
    }
}

// =============================================================================
// POD-LSTM Implementation
// =============================================================================

PODLSTM::PODLSTM(const PODLSTMConfig& config) : config_(config) {
    for (int l = 0; l < config.num_lstm_layers; ++l) {
        int input_size = (l == 0) ? config.num_modes : config.lstm_hidden_size;
        lstm_layers_.push_back(std::make_unique<LSTMCell>(input_size, config.lstm_hidden_size));
    }
    
    output_layer_ = std::make_unique<DenseLayer>(config.lstm_hidden_size, config.num_modes);
}

Tensor PODLSTM::predictNext(const std::vector<Tensor>& history) {
    // Initialize states
    h_states_.resize(config_.num_lstm_layers, Tensor({config_.lstm_hidden_size}, 0.0));
    c_states_.resize(config_.num_lstm_layers, Tensor({config_.lstm_hidden_size}, 0.0));
    
    // Process sequence
    Tensor output;
    for (const auto& input : history) {
        Tensor h = input;
        for (size_t l = 0; l < lstm_layers_.size(); ++l) {
            auto [h_new, c_new] = lstm_layers_[l]->forward(h, h_states_[l], c_states_[l]);
            h_states_[l] = h_new;
            c_states_[l] = c_new;
            h = h_new;
        }
        output = h;
    }
    
    // Final output projection
    return output_layer_->forward(output);
}

std::vector<Tensor> PODLSTM::generate(const std::vector<Tensor>& initial_sequence,
                                       int num_steps) {
    std::vector<Tensor> trajectory;
    std::vector<Tensor> history = initial_sequence;
    
    for (int step = 0; step < num_steps; ++step) {
        Tensor next = predictNext(history);
        trajectory.push_back(next);
        
        // Update history
        if (history.size() >= (size_t)config_.sequence_length) {
            history.erase(history.begin());
        }
        history.push_back(next);
    }
    
    return trajectory;
}

void PODLSTM::save(const std::string& path) const {
    pod_basis_.save(path + ".pod");
    // Save LSTM weights
}

void PODLSTM::load(const std::string& path) {
    pod_basis_.load(path + ".pod");
    // Load LSTM weights
}

// =============================================================================
// Autoencoder ROM Implementation
// =============================================================================

AutoencoderROM::AutoencoderROM(const AutoencoderROMConfig& config) : config_(config) {
    // Build encoder
    MLPConfig enc_config;
    enc_config.layer_sizes = config.encoder_layers;
    enc_config.layer_sizes.push_back(config.latent_dim);
    enc_config.activation = config.activation;
    encoder_ = std::make_unique<MLPModel>(enc_config);
    
    // Build decoder (mirror of encoder)
    MLPConfig dec_config;
    if (config.decoder_layers.empty()) {
        // Mirror encoder
        for (int i = config.encoder_layers.size() - 1; i >= 0; --i) {
            dec_config.layer_sizes.push_back(config.encoder_layers[i]);
        }
    } else {
        dec_config.layer_sizes = config.decoder_layers;
    }
    dec_config.layer_sizes.insert(dec_config.layer_sizes.begin(), config.latent_dim);
    dec_config.activation = config.activation;
    decoder_ = std::make_unique<MLPModel>(dec_config);
    
    // Dynamics network
    MLPConfig dyn_config;
    dyn_config.layer_sizes = config.dynamics_layers;
    dyn_config.layer_sizes.insert(dyn_config.layer_sizes.begin(), config.latent_dim);
    dyn_config.layer_sizes.push_back(config.latent_dim);
    dyn_config.activation = config.activation;
    dynamics_ = std::make_unique<MLPModel>(dyn_config);
    
    // For VAE
    if (config.variational) {
        MLPConfig mu_config;
        mu_config.layer_sizes = {config.encoder_layers.back(), config.latent_dim};
        mu_layer_ = std::make_unique<MLPModel>(mu_config);
        logvar_layer_ = std::make_unique<MLPModel>(mu_config);
    }
}

Tensor AutoencoderROM::encode(const Tensor& full_state) {
    Tensor encoded = encoder_->forward(full_state);
    
    if (config_.variational && mu_layer_ && logvar_layer_) {
        Tensor mu = mu_layer_->forward(encoded);
        Tensor logvar = logvar_layer_->forward(encoded);
        return reparameterize(mu, logvar);
    }
    
    return encoded;
}

Tensor AutoencoderROM::decode(const Tensor& latent_state) {
    return decoder_->forward(latent_state);
}

Tensor AutoencoderROM::reparameterize(const Tensor& mu, const Tensor& logvar) {
    Tensor result(mu.shape);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dist(0.0, 1.0);
    
    for (size_t i = 0; i < mu.numel(); ++i) {
        double std = std::exp(0.5 * logvar.data[i]);
        result.data[i] = mu.data[i] + std * dist(gen);
    }
    
    return result;
}

void AutoencoderROM::trainAutoencoder(const std::vector<Tensor>& snapshots) {
    std::vector<Tensor*> all_params;
    
    auto enc_params = encoder_->parameters();
    auto dec_params = decoder_->parameters();
    all_params.insert(all_params.end(), enc_params.begin(), enc_params.end());
    all_params.insert(all_params.end(), dec_params.begin(), dec_params.end());
    
    AdamOptimizer optimizer(all_params, config_.learning_rate);
    
    int n = snapshots.size();
    
    std::cout << "Training autoencoder on " << n << " snapshots...\n";
    
    for (int epoch = 0; epoch < config_.epochs; ++epoch) {
        double total_loss = 0.0;
        
        std::vector<int> indices(n);
        std::iota(indices.begin(), indices.end(), 0);
        std::shuffle(indices.begin(), indices.end(), std::mt19937(epoch));
        
        for (int b = 0; b < n; b += config_.batch_size) {
            int bs = std::min(config_.batch_size, n - b);
            
            // Stack batch
            int dim = snapshots[0].numel();
            Tensor batch({bs, dim});
            for (int i = 0; i < bs; ++i) {
                for (int j = 0; j < dim; ++j) {
                    batch.data[i * dim + j] = snapshots[indices[b + i]].data[j];
                }
            }
            
            encoder_->zeroGrad();
            decoder_->zeroGrad();
            
            Tensor latent = encoder_->forward(batch);
            Tensor reconstructed = decoder_->forward(latent);
            
            double loss = 0.0;
            Tensor grad_recon({bs, dim});
            for (int i = 0; i < bs * dim; ++i) {
                double diff = reconstructed.data[i] - batch.data[i];
                loss += diff * diff;
                grad_recon.data[i] = 2.0 * diff / (bs * dim);
            }
            loss /= (bs * dim);
            total_loss += loss * bs;
            
            Tensor grad_latent = decoder_->backward(grad_recon);
            encoder_->backward(grad_latent);
            
            optimizer.step(all_params);
        }
        
        if ((epoch + 1) % 50 == 0) {
            std::cout << "Epoch " << epoch + 1 << " - reconstruction loss: " 
                      << total_loss / n << "\n";
        }
    }
}

void AutoencoderROM::trainDynamics(const std::vector<Tensor>& snapshots,
                                    const std::vector<double>& times) {
    // Encode all snapshots
    std::vector<Tensor> latent_states;
    for (const auto& snap : snapshots) {
        latent_states.push_back(encode(snap));
    }
    
    // Compute velocities
    std::vector<Tensor> latent_velocities;
    for (size_t i = 0; i < latent_states.size() - 1; ++i) {
        double dt = times[i + 1] - times[i];
        Tensor velocity({config_.latent_dim});
        for (int j = 0; j < config_.latent_dim; ++j) {
            velocity.data[j] = (latent_states[i + 1].data[j] - 
                               latent_states[i].data[j]) / dt;
        }
        latent_velocities.push_back(std::move(velocity));
    }
    
    // Train dynamics
    AdamOptimizer optimizer(dynamics_->parameters(), config_.learning_rate);
    int n = latent_velocities.size();
    
    std::cout << "Training dynamics network on " << n << " samples...\n";
    
    for (int epoch = 0; epoch < config_.epochs; ++epoch) {
        double total_loss = 0.0;
        
        for (int i = 0; i < n; ++i) {
            dynamics_->zeroGrad();
            
            Tensor pred = dynamics_->forward(latent_states[i]);
            double loss = dynamics_->computeLoss(pred, latent_velocities[i]);
            total_loss += loss;
            
            Tensor grad = dynamics_->computeLossGrad(pred, latent_velocities[i]);
            dynamics_->backward(grad);
            optimizer.step(dynamics_->gradients());
        }
        
        if ((epoch + 1) % 50 == 0) {
            std::cout << "Epoch " << epoch + 1 << " - dynamics loss: " 
                      << total_loss / n << "\n";
        }
    }
}

void AutoencoderROM::train(const std::vector<Tensor>& snapshots,
                           const std::vector<double>& times) {
    trainAutoencoder(snapshots);
    trainDynamics(snapshots, times);
}

Tensor AutoencoderROM::solve(const Tensor& initial_condition,
                              double t_final, double dt) {
    Tensor z = encode(initial_condition);
    
    double t = 0.0;
    while (t < t_final - 1e-10) {
        // Simple forward Euler in latent space
        Tensor dz = dynamics_->predict(z);
        for (size_t i = 0; i < z.numel(); ++i) {
            z.data[i] += dt * dz.data[i];
        }
        t += dt;
    }
    
    return decode(z);
}

double AutoencoderROM::reconstructionError(const std::vector<Tensor>& test_data) {
    double total_error = 0.0;
    
    for (const auto& data : test_data) {
        Tensor latent = encode(data);
        Tensor recon = decode(latent);
        
        double error = 0.0, norm = 0.0;
        for (size_t i = 0; i < data.numel(); ++i) {
            double diff = recon.data[i] - data.data[i];
            error += diff * diff;
            norm += data.data[i] * data.data[i];
        }
        total_error += std::sqrt(error / (norm + 1e-14));
    }
    
    return total_error / test_data.size();
}

void AutoencoderROM::save(const std::string& path) const {
    encoder_->save(path + ".encoder");
    decoder_->save(path + ".decoder");
    dynamics_->save(path + ".dynamics");
}

void AutoencoderROM::load(const std::string& path) {
    encoder_->load(path + ".encoder");
    decoder_->load(path + ".decoder");
    dynamics_->load(path + ".dynamics");
}

// =============================================================================
// DeepONet ROM Implementation
// =============================================================================

DeepONetROM::DeepONetROM(const DeepONetConfig& config) : config_(config) {
    // Branch network
    MLPConfig branch_config;
    branch_config.layer_sizes = config.branch_layers;
    branch_config.layer_sizes.push_back(config.num_basis);
    branch_net_ = std::make_unique<MLPModel>(branch_config);
    
    // Trunk network
    MLPConfig trunk_config;
    trunk_config.layer_sizes = config.trunk_layers;
    trunk_config.layer_sizes.push_back(config.num_basis);
    trunk_net_ = std::make_unique<MLPModel>(trunk_config);
    
    // Bias
    branch_bias_ = Tensor({1}, 0.0);
}

Tensor DeepONetROM::evaluate(const Tensor& input_function,
                              const Tensor& query_locations) {
    // Branch: encode input function
    Tensor branch_output = branch_net_->predict(input_function);
    
    int num_queries = query_locations.shape[0];
    int coord_dim = config_.coord_dim;
    
    Tensor result({num_queries});
    
    for (int q = 0; q < num_queries; ++q) {
        // Extract query location
        Tensor loc({1, coord_dim});
        for (int d = 0; d < coord_dim; ++d) {
            loc.data[d] = query_locations.data[q * coord_dim + d];
        }
        
        // Trunk: encode location
        Tensor trunk_output = trunk_net_->predict(loc);
        
        // Inner product
        double value = branch_bias_.data[0];
        for (int b = 0; b < config_.num_basis; ++b) {
            value += branch_output.data[b] * trunk_output.data[b];
        }
        
        result.data[q] = value;
    }
    
    return result;
}

void DeepONetROM::train(const std::vector<Tensor>& input_functions,
                        const std::vector<Tensor>& output_functions,
                        const std::vector<Tensor>& locations) {
    std::vector<Tensor*> all_params;
    
    auto branch_params = branch_net_->parameters();
    auto trunk_params = trunk_net_->parameters();
    all_params.insert(all_params.end(), branch_params.begin(), branch_params.end());
    all_params.insert(all_params.end(), trunk_params.begin(), trunk_params.end());
    all_params.push_back(&branch_bias_);
    
    AdamOptimizer optimizer(all_params, config_.learning_rate);
    
    int n = input_functions.size();
    
    std::cout << "Training DeepONet on " << n << " samples...\n";
    
    for (int epoch = 0; epoch < config_.epochs; ++epoch) {
        double total_loss = 0.0;
        
        for (int i = 0; i < n; ++i) {
            branch_net_->zeroGrad();
            trunk_net_->zeroGrad();
            
            Tensor pred = evaluate(input_functions[i], locations[i]);
            
            double loss = 0.0;
            for (size_t j = 0; j < pred.numel(); ++j) {
                double diff = pred.data[j] - output_functions[i].data[j];
                loss += diff * diff;
            }
            loss /= pred.numel();
            total_loss += loss;
            
            // Simplified backward - full implementation would need chain rule
        }
        
        if ((epoch + 1) % 20 == 0) {
            std::cout << "Epoch " << epoch + 1 << " - loss: " << total_loss / n << "\n";
        }
    }
}

void DeepONetROM::save(const std::string& path) const {
    branch_net_->save(path + ".branch");
    trunk_net_->save(path + ".trunk");
}

void DeepONetROM::load(const std::string& path) {
    branch_net_->load(path + ".branch");
    trunk_net_->load(path + ".trunk");
}

// =============================================================================
// DEIM Hyper-Reduction Implementation
// =============================================================================

void DEIMHyperReduction::compute(const std::vector<Tensor>& nonlinear_snapshots,
                                  int num_points) {
    if (nonlinear_snapshots.empty()) return;
    
    // Compute POD basis of nonlinear snapshots
    PODBasis pod;
    pod.compute(nonlinear_snapshots, num_points);
    
    basis_ = pod.getBasis();
    int n_dof = basis_.shape[0];
    int n_modes = basis_.shape[1];
    
    // DEIM algorithm to find interpolation points
    indices_.clear();
    
    // First point: max of first basis vector
    int max_idx = 0;
    double max_val = std::abs(basis_.data[0]);
    for (int i = 1; i < n_dof; ++i) {
        if (std::abs(basis_.data[i]) > max_val) {
            max_val = std::abs(basis_.data[i]);
            max_idx = i;
        }
    }
    indices_.push_back(max_idx);
    
    // Subsequent points
    for (int j = 1; j < n_modes; ++j) {
        // Solve for c in P^T U_j-1 c = P^T u_j
        // Then find max of |u_j - U_j-1 c|
        
        // Extract j-th basis vector
        Tensor u_j({n_dof});
        for (int i = 0; i < n_dof; ++i) {
            u_j.data[i] = basis_.data[i * n_modes + j];
        }
        
        // Simple greedy selection: find max residual
        max_idx = 0;
        max_val = 0.0;
        
        for (int i = 0; i < n_dof; ++i) {
            bool already_selected = false;
            for (int idx : indices_) {
                if (idx == i) {
                    already_selected = true;
                    break;
                }
            }
            
            if (!already_selected && std::abs(u_j.data[i]) > max_val) {
                max_val = std::abs(u_j.data[i]);
                max_idx = i;
            }
        }
        
        indices_.push_back(max_idx);
    }
    
    // Build interpolation matrix
    interpolation_matrix_ = Tensor({(int)indices_.size(), n_modes});
    for (size_t i = 0; i < indices_.size(); ++i) {
        for (int j = 0; j < n_modes; ++j) {
            interpolation_matrix_.data[i * n_modes + j] = 
                basis_.data[indices_[i] * n_modes + j];
        }
    }
    
    std::cout << "DEIM: Selected " << indices_.size() << " interpolation points\n";
}

Tensor DEIMHyperReduction::approximate(const Tensor& sampled_values) {
    // Solve least squares for coefficients
    // Then reconstruct: f ≈ U * c
    
    int n_dof = basis_.shape[0];
    int n_modes = basis_.shape[1];
    
    // Simple pseudo-inverse (for small systems)
    Tensor c({n_modes}, 0.0);
    
    // c = (P^T U)^{-1} * sampled_values (simplified)
    for (int j = 0; j < n_modes; ++j) {
        double sum = 0.0;
        for (size_t i = 0; i < indices_.size(); ++i) {
            sum += interpolation_matrix_.data[i * n_modes + j] * sampled_values.data[i];
        }
        c.data[j] = sum;
    }
    
    // Reconstruct
    Tensor result({n_dof}, 0.0);
    for (int i = 0; i < n_dof; ++i) {
        for (int j = 0; j < n_modes; ++j) {
            result.data[i] += basis_.data[i * n_modes + j] * c.data[j];
        }
    }
    
    return result;
}

void DEIMHyperReduction::save(const std::string& path) const {
    basis_.save(path + ".basis");
    interpolation_matrix_.save(path + ".interp");
    
    std::ofstream f(path + ".indices");
    for (int idx : indices_) {
        f << idx << " ";
    }
}

void DEIMHyperReduction::load(const std::string& path) {
    basis_.load(path + ".basis");
    interpolation_matrix_.load(path + ".interp");
    
    indices_.clear();
    std::ifstream f(path + ".indices");
    int idx;
    while (f >> idx) {
        indices_.push_back(idx);
    }
}

// =============================================================================
// ROM Solver Implementation
// =============================================================================

ROMSolver::ROMSolver(const ROMConfig& config) : config_(config) {
    switch (config.type) {
        case ROMType::POD_NN:
            pod_nn_ = std::make_unique<PODNeuralROM>(config.pod_nn_config);
            break;
        case ROMType::POD_LSTM:
            pod_lstm_ = std::make_unique<PODLSTM>(config.pod_lstm_config);
            break;
        case ROMType::AUTOENCODER:
            autoencoder_ = std::make_unique<AutoencoderROM>(config.autoencoder_config);
            break;
        case ROMType::DEEPONET:
            deeponet_ = std::make_unique<DeepONetROM>(config.deeponet_config);
            break;
    }
}

Tensor ROMSolver::vecToTensor(Vec v) const {
    PetscInt n;
    VecGetSize(v, &n);
    
    const PetscScalar* arr;
    VecGetArrayRead(v, &arr);
    
    Tensor t({(int)n});
    for (int i = 0; i < n; ++i) {
        t.data[i] = arr[i];
    }
    
    VecRestoreArrayRead(v, &arr);
    return t;
}

void ROMSolver::tensorToVec(const Tensor& t, Vec v) const {
    PetscScalar* arr;
    VecGetArray(v, &arr);
    
    for (size_t i = 0; i < t.numel(); ++i) {
        arr[i] = t.data[i];
    }
    
    VecRestoreArray(v, &arr);
}

void ROMSolver::build(const std::vector<Vec>& snapshots,
                      const std::vector<double>& times,
                      const std::vector<double>& parameters) {
    std::vector<Tensor> tensor_snapshots;
    for (const auto& v : snapshots) {
        tensor_snapshots.push_back(vecToTensor(v));
    }
    
    PetscInt n;
    VecGetSize(snapshots[0], &n);
    stats_.full_dimension = n;
    
    switch (config_.type) {
        case ROMType::POD_NN:
            pod_nn_->buildFromSnapshots(tensor_snapshots, times);
            stats_.latent_dimension = pod_nn_->getBasis().numModes();
            stats_.energy_captured = pod_nn_->getBasis().energyCaptured();
            break;
        case ROMType::AUTOENCODER:
            autoencoder_->train(tensor_snapshots, times);
            stats_.latent_dimension = config_.autoencoder_config.latent_dim;
            break;
        default:
            break;
    }
    
    stats_.compression_ratio = (double)stats_.full_dimension / stats_.latent_dimension;
}

PetscErrorCode ROMSolver::solve(Vec initial, double t_final, Vec result) {
    Tensor init = vecToTensor(initial);
    Tensor sol;
    
    switch (config_.type) {
        case ROMType::POD_NN:
            sol = pod_nn_->solve(init, t_final, 0.01);
            break;
        case ROMType::AUTOENCODER:
            sol = autoencoder_->solve(init, t_final, 0.01);
            break;
        default:
            return PETSC_ERR_SUP;
    }
    
    tensorToVec(sol, result);
    return 0;
}

void ROMSolver::saveModel(const std::string& path) const {
    switch (config_.type) {
        case ROMType::POD_NN:
            pod_nn_->save(path);
            break;
        case ROMType::POD_LSTM:
            pod_lstm_->save(path);
            break;
        case ROMType::AUTOENCODER:
            autoencoder_->save(path);
            break;
        case ROMType::DEEPONET:
            deeponet_->save(path);
            break;
    }
}

void ROMSolver::loadModel(const std::string& path) {
    switch (config_.type) {
        case ROMType::POD_NN:
            pod_nn_->load(path);
            break;
        case ROMType::POD_LSTM:
            pod_lstm_->load(path);
            break;
        case ROMType::AUTOENCODER:
            autoencoder_->load(path);
            break;
        case ROMType::DEEPONET:
            deeponet_->load(path);
            break;
    }
}

// =============================================================================
// Configuration Parsing
// =============================================================================

PODNNConfig parsePODNNConfig(const std::map<std::string, std::string>& config) {
    PODNNConfig result;
    
    auto getInt = [&](const std::string& key, int def) {
        auto it = config.find(key);
        return it != config.end() ? std::stoi(it->second) : def;
    };
    
    auto getDouble = [&](const std::string& key, double def) {
        auto it = config.find(key);
        return it != config.end() ? std::stod(it->second) : def;
    };
    
    result.num_modes = getInt("rom_num_modes", 50);
    result.energy_threshold = getDouble("rom_energy_threshold", 0.999);
    result.learning_rate = getDouble("rom_learning_rate", 1e-3);
    result.epochs = getInt("rom_epochs", 200);
    result.batch_size = getInt("rom_batch_size", 64);
    
    return result;
}

AutoencoderROMConfig parseAutoencoderConfig(const std::map<std::string, std::string>& config) {
    AutoencoderROMConfig result;
    
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
    
    result.latent_dim = getInt("autoencoder_latent_dim", 20);
    result.variational = getBool("autoencoder_variational", false);
    result.kl_weight = getDouble("autoencoder_kl_weight", 1e-3);
    result.learning_rate = getDouble("autoencoder_learning_rate", 1e-3);
    result.epochs = getInt("autoencoder_epochs", 500);
    
    return result;
}

} // namespace ML
} // namespace FSRM
