/**
 * @file FourierNeuralOperator.cpp
 * @brief Implementation of Fourier Neural Operator
 */

#include "FourierNeuralOperator.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cassert>

namespace FSRM {
namespace ML {

// =============================================================================
// Tensor Implementation
// =============================================================================

Tensor::Tensor(const std::vector<int>& shape) : shape(shape) {
    size_t n = 1;
    for (int s : shape) n *= s;
    data.resize(n, 0.0);
}

Tensor::Tensor(const std::vector<int>& shape, double value) : shape(shape) {
    size_t n = 1;
    for (int s : shape) n *= s;
    data.resize(n, value);
}

Tensor::Tensor(const std::vector<int>& shape, const std::vector<double>& data)
    : shape(shape), data(data) {}

size_t Tensor::size() const {
    return data.size();
}

size_t Tensor::numel() const {
    size_t n = 1;
    for (int s : shape) n *= s;
    return n;
}

double& Tensor::operator()(int i) {
    return data[i];
}

double& Tensor::operator()(int i, int j) {
    return data[i * shape[1] + j];
}

double& Tensor::operator()(int i, int j, int k) {
    return data[(i * shape[1] + j) * shape[2] + k];
}

double& Tensor::operator()(int i, int j, int k, int l) {
    return data[((i * shape[1] + j) * shape[2] + k) * shape[3] + l];
}

const double& Tensor::operator()(int i) const {
    return data[i];
}

const double& Tensor::operator()(int i, int j) const {
    return data[i * shape[1] + j];
}

const double& Tensor::operator()(int i, int j, int k) const {
    return data[(i * shape[1] + j) * shape[2] + k];
}

const double& Tensor::operator()(int i, int j, int k, int l) const {
    return data[((i * shape[1] + j) * shape[2] + k) * shape[3] + l];
}

Tensor Tensor::reshape(const std::vector<int>& new_shape) const {
    Tensor result;
    result.data = data;
    result.shape = new_shape;
    return result;
}

Tensor Tensor::transpose(int dim1, int dim2) const {
    // Simple transpose for 2D
    if (shape.size() == 2 && dim1 == 0 && dim2 == 1) {
        Tensor result({shape[1], shape[0]});
        for (int i = 0; i < shape[0]; ++i) {
            for (int j = 0; j < shape[1]; ++j) {
                result(j, i) = (*this)(i, j);
            }
        }
        return result;
    }
    // General transpose - swap dimensions
    std::vector<int> new_shape = shape;
    std::swap(new_shape[dim1], new_shape[dim2]);
    Tensor result(new_shape);
    // TODO: General implementation
    return result;
}

Tensor Tensor::operator+(const Tensor& other) const {
    assert(data.size() == other.data.size());
    Tensor result(shape);
    for (size_t i = 0; i < data.size(); ++i) {
        result.data[i] = data[i] + other.data[i];
    }
    return result;
}

Tensor Tensor::operator-(const Tensor& other) const {
    assert(data.size() == other.data.size());
    Tensor result(shape);
    for (size_t i = 0; i < data.size(); ++i) {
        result.data[i] = data[i] - other.data[i];
    }
    return result;
}

Tensor Tensor::operator*(const Tensor& other) const {
    assert(data.size() == other.data.size());
    Tensor result(shape);
    for (size_t i = 0; i < data.size(); ++i) {
        result.data[i] = data[i] * other.data[i];
    }
    return result;
}

Tensor Tensor::operator*(double scalar) const {
    Tensor result(shape);
    for (size_t i = 0; i < data.size(); ++i) {
        result.data[i] = data[i] * scalar;
    }
    return result;
}

Tensor& Tensor::operator+=(const Tensor& other) {
    for (size_t i = 0; i < data.size(); ++i) {
        data[i] += other.data[i];
    }
    return *this;
}

double Tensor::sum() const {
    return std::accumulate(data.begin(), data.end(), 0.0);
}

double Tensor::mean() const {
    return sum() / data.size();
}

double Tensor::max() const {
    return *std::max_element(data.begin(), data.end());
}

double Tensor::min() const {
    return *std::min_element(data.begin(), data.end());
}

double Tensor::norm() const {
    double sum_sq = 0.0;
    for (double x : data) sum_sq += x * x;
    return std::sqrt(sum_sq);
}

void Tensor::zeros() {
    std::fill(data.begin(), data.end(), 0.0);
}

void Tensor::ones() {
    std::fill(data.begin(), data.end(), 1.0);
}

void Tensor::randn(double mean, double std) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dist(mean, std);
    for (double& x : data) x = dist(gen);
}

void Tensor::xavier_init() {
    double scale = std::sqrt(6.0 / (shape[0] + shape[1]));
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(-scale, scale);
    for (double& x : data) x = dist(gen);
}

void Tensor::kaiming_init() {
    double std = std::sqrt(2.0 / shape[0]);
    randn(0.0, std);
}

void Tensor::save(const std::string& path) const {
    std::ofstream file(path, std::ios::binary);
    int ndim = shape.size();
    file.write(reinterpret_cast<const char*>(&ndim), sizeof(int));
    file.write(reinterpret_cast<const char*>(shape.data()), ndim * sizeof(int));
    size_t n = data.size();
    file.write(reinterpret_cast<const char*>(&n), sizeof(size_t));
    file.write(reinterpret_cast<const char*>(data.data()), n * sizeof(double));
}

void Tensor::load(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    int ndim;
    file.read(reinterpret_cast<char*>(&ndim), sizeof(int));
    shape.resize(ndim);
    file.read(reinterpret_cast<char*>(shape.data()), ndim * sizeof(int));
    size_t n;
    file.read(reinterpret_cast<char*>(&n), sizeof(size_t));
    data.resize(n);
    file.read(reinterpret_cast<char*>(data.data()), n * sizeof(double));
}

void Tensor::toGPU() {
    // TODO: CUDA implementation
    on_gpu = true;
}

void Tensor::toCPU() {
    on_gpu = false;
}

// =============================================================================
// ComplexTensor Implementation
// =============================================================================

ComplexTensor::ComplexTensor(const std::vector<int>& shape) : shape(shape) {
    size_t n = 1;
    for (int s : shape) n *= s;
    data.resize(n, {0.0, 0.0});
}

size_t ComplexTensor::numel() const {
    size_t n = 1;
    for (int s : shape) n *= s;
    return n;
}

std::complex<double>& ComplexTensor::operator()(int i, int j, int k, int l) {
    return data[((i * shape[1] + j) * shape[2] + k) * shape[3] + l];
}

const std::complex<double>& ComplexTensor::operator()(int i, int j, int k, int l) const {
    return data[((i * shape[1] + j) * shape[2] + k) * shape[3] + l];
}

ComplexTensor ComplexTensor::operator*(const ComplexTensor& other) const {
    ComplexTensor result(shape);
    for (size_t i = 0; i < data.size(); ++i) {
        result.data[i] = data[i] * other.data[i];
    }
    return result;
}

ComplexTensor& ComplexTensor::operator+=(const ComplexTensor& other) {
    for (size_t i = 0; i < data.size(); ++i) {
        data[i] += other.data[i];
    }
    return *this;
}

// =============================================================================
// Activation Functions
// =============================================================================

namespace Activation {

Tensor relu(const Tensor& x) {
    Tensor result(x.shape);
    for (size_t i = 0; i < x.data.size(); ++i) {
        result.data[i] = std::max(0.0, x.data[i]);
    }
    return result;
}

Tensor gelu(const Tensor& x) {
    // GELU(x) = 0.5 * x * (1 + tanh(sqrt(2/π) * (x + 0.044715 * x³)))
    static const double sqrt_2_pi = std::sqrt(2.0 / M_PI);
    Tensor result(x.shape);
    for (size_t i = 0; i < x.data.size(); ++i) {
        double xi = x.data[i];
        double inner = sqrt_2_pi * (xi + 0.044715 * xi * xi * xi);
        result.data[i] = 0.5 * xi * (1.0 + std::tanh(inner));
    }
    return result;
}

Tensor silu(const Tensor& x) {
    // SiLU(x) = x * sigmoid(x)
    Tensor result(x.shape);
    for (size_t i = 0; i < x.data.size(); ++i) {
        double xi = x.data[i];
        result.data[i] = xi / (1.0 + std::exp(-xi));
    }
    return result;
}

Tensor tanh(const Tensor& x) {
    Tensor result(x.shape);
    for (size_t i = 0; i < x.data.size(); ++i) {
        result.data[i] = std::tanh(x.data[i]);
    }
    return result;
}

Tensor relu_backward(const Tensor& x, const Tensor& grad) {
    Tensor result(x.shape);
    for (size_t i = 0; i < x.data.size(); ++i) {
        result.data[i] = x.data[i] > 0 ? grad.data[i] : 0.0;
    }
    return result;
}

Tensor gelu_backward(const Tensor& x, const Tensor& grad) {
    static const double sqrt_2_pi = std::sqrt(2.0 / M_PI);
    Tensor result(x.shape);
    for (size_t i = 0; i < x.data.size(); ++i) {
        double xi = x.data[i];
        double inner = sqrt_2_pi * (xi + 0.044715 * xi * xi * xi);
        double tanh_inner = std::tanh(inner);
        double sech2 = 1.0 - tanh_inner * tanh_inner;
        double d_inner = sqrt_2_pi * (1.0 + 0.134145 * xi * xi);
        double d_gelu = 0.5 * (1.0 + tanh_inner) + 0.5 * xi * sech2 * d_inner;
        result.data[i] = d_gelu * grad.data[i];
    }
    return result;
}

Tensor silu_backward(const Tensor& x, const Tensor& grad) {
    Tensor result(x.shape);
    for (size_t i = 0; i < x.data.size(); ++i) {
        double xi = x.data[i];
        double sig = 1.0 / (1.0 + std::exp(-xi));
        double d_silu = sig * (1.0 + xi * (1.0 - sig));
        result.data[i] = d_silu * grad.data[i];
    }
    return result;
}

Tensor tanh_backward(const Tensor& x, const Tensor& grad) {
    Tensor result(x.shape);
    for (size_t i = 0; i < x.data.size(); ++i) {
        double t = std::tanh(x.data[i]);
        result.data[i] = (1.0 - t * t) * grad.data[i];
    }
    return result;
}

ActivationFn getActivation(const std::string& name) {
    if (name == "relu") return relu;
    if (name == "gelu") return gelu;
    if (name == "silu") return silu;
    if (name == "tanh") return tanh;
    return gelu;  // Default
}

ActivationBackwardFn getActivationBackward(const std::string& name) {
    if (name == "relu") return relu_backward;
    if (name == "gelu") return gelu_backward;
    if (name == "silu") return silu_backward;
    if (name == "tanh") return tanh_backward;
    return gelu_backward;  // Default
}

} // namespace Activation

// =============================================================================
// FFT Implementation (Cooley-Tukey)
// =============================================================================

void FFT::fft1d(const double* input, std::complex<double>* output, int n) {
    // Base case
    if (n == 1) {
        output[0] = std::complex<double>(input[0], 0.0);
        return;
    }
    
    // Separate even and odd
    std::vector<double> even(n / 2), odd(n / 2);
    for (int i = 0; i < n / 2; ++i) {
        even[i] = input[2 * i];
        odd[i] = input[2 * i + 1];
    }
    
    // Recursive FFT
    std::vector<std::complex<double>> even_ft(n / 2), odd_ft(n / 2);
    fft1d(even.data(), even_ft.data(), n / 2);
    fft1d(odd.data(), odd_ft.data(), n / 2);
    
    // Combine
    for (int k = 0; k < n / 2; ++k) {
        double angle = -2.0 * M_PI * k / n;
        std::complex<double> w(std::cos(angle), std::sin(angle));
        output[k] = even_ft[k] + w * odd_ft[k];
        output[k + n / 2] = even_ft[k] - w * odd_ft[k];
    }
}

void FFT::ifft1d(const std::complex<double>* input, double* output, int n) {
    // IFFT = conjugate -> FFT -> conjugate / n
    std::vector<std::complex<double>> conj_input(n);
    for (int i = 0; i < n; ++i) {
        conj_input[i] = std::conj(input[i]);
    }
    
    std::vector<std::complex<double>> result(n);
    std::vector<double> real_input(n);
    for (int i = 0; i < n; ++i) {
        real_input[i] = conj_input[i].real();
    }
    fft1d(real_input.data(), result.data(), n);
    
    for (int i = 0; i < n; ++i) {
        output[i] = result[i].real() / n;
    }
}

void FFT::fft2d(const Tensor& input, ComplexTensor& output) {
    int batch = input.shape[0];
    int channels = input.shape[1];
    int nx = input.shape[2];
    int ny = input.shape[3];
    
    output = ComplexTensor({batch, channels, nx, ny});
    
    // FFT along each dimension
    for (int b = 0; b < batch; ++b) {
        for (int c = 0; c < channels; ++c) {
            // Row-wise FFT
            std::vector<std::complex<double>> row_ft(ny);
            for (int i = 0; i < nx; ++i) {
                std::vector<double> row(ny);
                for (int j = 0; j < ny; ++j) {
                    row[j] = input(b, c, i, j);
                }
                fft1d(row.data(), row_ft.data(), ny);
                for (int j = 0; j < ny; ++j) {
                    output(b, c, i, j) = row_ft[j];
                }
            }
            
            // Column-wise FFT
            std::vector<std::complex<double>> col(nx), col_ft(nx);
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    col[i] = output(b, c, i, j);
                }
                // FFT on complex values
                for (int i = 0; i < nx; ++i) {
                    col_ft[i] = {0, 0};
                    for (int k = 0; k < nx; ++k) {
                        double angle = -2.0 * M_PI * i * k / nx;
                        std::complex<double> w(std::cos(angle), std::sin(angle));
                        col_ft[i] += col[k] * w;
                    }
                }
                for (int i = 0; i < nx; ++i) {
                    output(b, c, i, j) = col_ft[i];
                }
            }
        }
    }
}

void FFT::ifft2d(const ComplexTensor& input, Tensor& output) {
    int batch = input.shape[0];
    int channels = input.shape[1];
    int nx = input.shape[2];
    int ny = input.shape[3];
    
    output = Tensor({batch, channels, nx, ny});
    
    // IFFT along each dimension
    for (int b = 0; b < batch; ++b) {
        for (int c = 0; c < channels; ++c) {
            // Column-wise IFFT first
            std::vector<std::complex<double>> col(nx), col_ft(nx);
            ComplexTensor temp({1, 1, nx, ny});
            
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    col[i] = input(b, c, i, j);
                }
                for (int i = 0; i < nx; ++i) {
                    col_ft[i] = {0, 0};
                    for (int k = 0; k < nx; ++k) {
                        double angle = 2.0 * M_PI * i * k / nx;
                        std::complex<double> w(std::cos(angle), std::sin(angle));
                        col_ft[i] += col[k] * w;
                    }
                    col_ft[i] /= nx;
                }
                for (int i = 0; i < nx; ++i) {
                    temp(0, 0, i, j) = col_ft[i];
                }
            }
            
            // Row-wise IFFT
            for (int i = 0; i < nx; ++i) {
                std::vector<std::complex<double>> row(ny), row_ft(ny);
                for (int j = 0; j < ny; ++j) {
                    row[j] = temp(0, 0, i, j);
                }
                for (int j = 0; j < ny; ++j) {
                    row_ft[j] = {0, 0};
                    for (int k = 0; k < ny; ++k) {
                        double angle = 2.0 * M_PI * j * k / ny;
                        std::complex<double> w(std::cos(angle), std::sin(angle));
                        row_ft[j] += row[k] * w;
                    }
                    row_ft[j] /= ny;
                }
                for (int j = 0; j < ny; ++j) {
                    output(b, c, i, j) = row_ft[j].real();
                }
            }
        }
    }
}

void FFT::rfft2d(const Tensor& input, ComplexTensor& output) {
    // Real-input FFT (half output due to symmetry)
    fft2d(input, output);
}

void FFT::irfft2d(const ComplexTensor& input, Tensor& output, const std::vector<int>& shape) {
    ifft2d(input, output);
}

// =============================================================================
// SpectralConv2d Implementation
// =============================================================================

SpectralConv2d::SpectralConv2d(int in_channels, int out_channels, int modes1, int modes2)
    : in_channels(in_channels), out_channels(out_channels), modes1(modes1), modes2(modes2) {
    
    // Initialize complex weights
    weights = ComplexTensor({in_channels, out_channels, modes1, modes2});
    grad_weights = ComplexTensor({in_channels, out_channels, modes1, modes2});
    
    // Xavier initialization for complex weights
    double scale = std::sqrt(2.0 / (in_channels + out_channels));
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dist(0.0, scale);
    
    for (size_t i = 0; i < weights.data.size(); ++i) {
        weights.data[i] = std::complex<double>(dist(gen), dist(gen));
    }
}

Tensor SpectralConv2d::forward(const Tensor& x) {
    // x: [batch, in_channels, nx, ny]
    int batch = x.shape[0];
    int nx = x.shape[2];
    int ny = x.shape[3];
    
    last_input_shape = x.shape;
    
    // FFT
    ComplexTensor x_ft({batch, in_channels, nx, ny});
    FFT::fft2d(x, x_ft);
    last_x_ft = x_ft;
    
    // Multiply with weights (only for low-frequency modes)
    ComplexTensor out_ft({batch, out_channels, nx, ny});
    
    for (int b = 0; b < batch; ++b) {
        for (int oc = 0; oc < out_channels; ++oc) {
            for (int i = 0; i < std::min(modes1, nx); ++i) {
                for (int j = 0; j < std::min(modes2, ny); ++j) {
                    std::complex<double> sum(0, 0);
                    for (int ic = 0; ic < in_channels; ++ic) {
                        sum += x_ft(b, ic, i, j) * weights(ic, oc, i, j);
                    }
                    out_ft(b, oc, i, j) = sum;
                }
            }
            // Copy high-frequency modes (zeros for now)
            for (int i = modes1; i < nx; ++i) {
                for (int j = 0; j < ny; ++j) {
                    out_ft(b, oc, i, j) = {0, 0};
                }
            }
            for (int i = 0; i < nx; ++i) {
                for (int j = modes2; j < ny; ++j) {
                    out_ft(b, oc, i, j) = {0, 0};
                }
            }
        }
    }
    
    // IFFT
    Tensor output({batch, out_channels, nx, ny});
    FFT::ifft2d(out_ft, output);
    
    return output;
}

Tensor SpectralConv2d::backward(const Tensor& grad_output) {
    int batch = grad_output.shape[0];
    int nx = grad_output.shape[2];
    int ny = grad_output.shape[3];
    
    // FFT of gradient
    ComplexTensor grad_ft({batch, out_channels, nx, ny});
    FFT::fft2d(grad_output, grad_ft);
    
    // Gradient w.r.t. weights
    for (int ic = 0; ic < in_channels; ++ic) {
        for (int oc = 0; oc < out_channels; ++oc) {
            for (int i = 0; i < modes1; ++i) {
                for (int j = 0; j < modes2; ++j) {
                    std::complex<double> sum(0, 0);
                    for (int b = 0; b < batch; ++b) {
                        sum += std::conj(last_x_ft(b, ic, i, j)) * grad_ft(b, oc, i, j);
                    }
                    grad_weights(ic, oc, i, j) = sum;
                }
            }
        }
    }
    
    // Gradient w.r.t. input
    ComplexTensor grad_x_ft({batch, in_channels, nx, ny});
    
    for (int b = 0; b < batch; ++b) {
        for (int ic = 0; ic < in_channels; ++ic) {
            for (int i = 0; i < nx; ++i) {
                for (int j = 0; j < ny; ++j) {
                    grad_x_ft(b, ic, i, j) = {0, 0};
                }
            }
            for (int i = 0; i < modes1; ++i) {
                for (int j = 0; j < modes2; ++j) {
                    std::complex<double> sum(0, 0);
                    for (int oc = 0; oc < out_channels; ++oc) {
                        sum += grad_ft(b, oc, i, j) * std::conj(weights(ic, oc, i, j));
                    }
                    grad_x_ft(b, ic, i, j) = sum;
                }
            }
        }
    }
    
    Tensor grad_input({batch, in_channels, nx, ny});
    FFT::ifft2d(grad_x_ft, grad_input);
    
    return grad_input;
}

// =============================================================================
// Linear Layer Implementation
// =============================================================================

Linear::Linear(int in_features, int out_features, bool bias)
    : has_bias(bias) {
    
    weight = Tensor({out_features, in_features});
    grad_weight = Tensor({out_features, in_features});
    weight.xavier_init();
    
    if (has_bias) {
        bias_vec = Tensor({out_features}, 0.0);
        grad_bias = Tensor({out_features});
    }
}

Tensor Linear::forward(const Tensor& x) {
    last_input = x;
    
    // x: [batch, ..., in_features] or [batch, in_features]
    // Output: [batch, ..., out_features]
    
    int in_f = weight.shape[1];
    int out_f = weight.shape[0];
    
    // Flatten batch dimensions
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

Tensor Linear::backward(const Tensor& grad_output) {
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
// LayerNorm Implementation
// =============================================================================

LayerNorm::LayerNorm(int normalized_shape, double eps)
    : normalized_shape(normalized_shape), eps(eps) {
    
    gamma = Tensor({normalized_shape}, 1.0);
    beta = Tensor({normalized_shape}, 0.0);
    grad_gamma = Tensor({normalized_shape});
    grad_beta = Tensor({normalized_shape});
}

Tensor LayerNorm::forward(const Tensor& x) {
    // Normalize along last dimension
    size_t batch_size = x.numel() / normalized_shape;
    
    x_normalized = Tensor(x.shape);
    std_inv = Tensor({(int)batch_size});
    
    for (size_t b = 0; b < batch_size; ++b) {
        // Compute mean
        double mean = 0.0;
        for (int i = 0; i < normalized_shape; ++i) {
            mean += x.data[b * normalized_shape + i];
        }
        mean /= normalized_shape;
        
        // Compute variance
        double var = 0.0;
        for (int i = 0; i < normalized_shape; ++i) {
            double diff = x.data[b * normalized_shape + i] - mean;
            var += diff * diff;
        }
        var /= normalized_shape;
        
        // Normalize
        double std_inv_val = 1.0 / std::sqrt(var + eps);
        std_inv(b) = std_inv_val;
        
        for (int i = 0; i < normalized_shape; ++i) {
            x_normalized.data[b * normalized_shape + i] = 
                (x.data[b * normalized_shape + i] - mean) * std_inv_val;
        }
    }
    
    // Apply gamma and beta
    Tensor output(x.shape);
    for (size_t b = 0; b < batch_size; ++b) {
        for (int i = 0; i < normalized_shape; ++i) {
            output.data[b * normalized_shape + i] = 
                gamma(i) * x_normalized.data[b * normalized_shape + i] + beta(i);
        }
    }
    
    return output;
}

Tensor LayerNorm::backward(const Tensor& grad_output) {
    size_t batch_size = grad_output.numel() / normalized_shape;
    
    // Gradient w.r.t. gamma and beta
    grad_gamma.zeros();
    grad_beta.zeros();
    
    for (size_t b = 0; b < batch_size; ++b) {
        for (int i = 0; i < normalized_shape; ++i) {
            grad_gamma(i) += grad_output.data[b * normalized_shape + i] * 
                            x_normalized.data[b * normalized_shape + i];
            grad_beta(i) += grad_output.data[b * normalized_shape + i];
        }
    }
    
    // Gradient w.r.t. input (simplified)
    Tensor grad_input(grad_output.shape);
    for (size_t b = 0; b < batch_size; ++b) {
        for (int i = 0; i < normalized_shape; ++i) {
            grad_input.data[b * normalized_shape + i] = 
                grad_output.data[b * normalized_shape + i] * gamma(i) * std_inv(b);
        }
    }
    
    return grad_input;
}

// =============================================================================
// FNOLayer Implementation
// =============================================================================

FNOLayer::FNOLayer(int channels, int modes1, int modes2,
                   const std::string& activation,
                   bool use_layer_norm, bool use_residual)
    : use_layer_norm(use_layer_norm), use_residual(use_residual) {
    
    spectral_conv = std::make_unique<SpectralConv2d>(channels, channels, modes1, modes2);
    linear = std::make_unique<Linear>(channels, channels);
    
    if (use_layer_norm) {
        norm = std::make_unique<LayerNorm>(channels);
    }
    
    activation_fn = Activation::getActivation(activation);
    activation_backward = Activation::getActivationBackward(activation);
}

Tensor FNOLayer::forward(const Tensor& x) {
    input_cache = x;
    
    // Spectral convolution
    Tensor x1 = spectral_conv->forward(x);
    
    // Linear (1x1 conv)
    // Reshape for linear: [batch, channels, nx, ny] -> [batch*nx*ny, channels]
    int batch = x.shape[0];
    int channels = x.shape[1];
    int nx = x.shape[2];
    int ny = x.shape[3];
    
    Tensor x_flat({batch * nx * ny, channels});
    for (int b = 0; b < batch; ++b) {
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                for (int c = 0; c < channels; ++c) {
                    x_flat((b * nx + i) * ny + j, c) = x(b, c, i, j);
                }
            }
        }
    }
    
    Tensor x2_flat = linear->forward(x_flat);
    
    // Reshape back
    Tensor x2({batch, channels, nx, ny});
    for (int b = 0; b < batch; ++b) {
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                for (int c = 0; c < channels; ++c) {
                    x2(b, c, i, j) = x2_flat((b * nx + i) * ny + j, c);
                }
            }
        }
    }
    
    // Combine
    Tensor combined = x1 + x2;
    
    // Residual connection
    if (use_residual) {
        combined = combined + x;
    }
    
    // Layer norm (if enabled)
    if (use_layer_norm) {
        // Apply per-channel normalization
        // Simplified: normalize entire tensor
        combined = norm->forward(combined.reshape({batch * channels, nx * ny}))
                       .reshape({batch, channels, nx, ny});
    }
    
    // Activation
    pre_activation = combined;
    return activation_fn(combined);
}

Tensor FNOLayer::backward(const Tensor& grad_output) {
    // Backward through activation
    Tensor grad = activation_backward(pre_activation, grad_output);
    
    // Backward through layer norm
    if (use_layer_norm) {
        int batch = grad.shape[0];
        int channels = grad.shape[1];
        int nx = grad.shape[2];
        int ny = grad.shape[3];
        grad = norm->backward(grad.reshape({batch * channels, nx * ny}))
                   .reshape({batch, channels, nx, ny});
    }
    
    // Residual gradient
    Tensor grad_residual = use_residual ? grad : Tensor(grad.shape, 0.0);
    
    // Backward through linear
    int batch = grad.shape[0];
    int channels = grad.shape[1];
    int nx = grad.shape[2];
    int ny = grad.shape[3];
    
    Tensor grad_flat({batch * nx * ny, channels});
    for (int b = 0; b < batch; ++b) {
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                for (int c = 0; c < channels; ++c) {
                    grad_flat((b * nx + i) * ny + j, c) = grad(b, c, i, j);
                }
            }
        }
    }
    
    Tensor grad_linear = linear->backward(grad_flat);
    
    // Reshape
    Tensor grad_linear_4d({batch, channels, nx, ny});
    for (int b = 0; b < batch; ++b) {
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                for (int c = 0; c < channels; ++c) {
                    grad_linear_4d(b, c, i, j) = grad_linear((b * nx + i) * ny + j, c);
                }
            }
        }
    }
    
    // Backward through spectral conv
    Tensor grad_spectral = spectral_conv->backward(grad);
    
    // Combine gradients
    return grad_spectral + grad_linear_4d + grad_residual;
}

// =============================================================================
// FNOModel Implementation
// =============================================================================

FNOModel::FNOModel(const FNOConfig& config) : config(config) {
    // Lifting layer: input_channels -> hidden_channels
    lifting = std::make_unique<Linear>(config.input_channels, config.hidden_channels);
    
    // FNO layers
    for (int i = 0; i < config.num_layers; ++i) {
        fno_layers.push_back(std::make_unique<FNOLayer>(
            config.hidden_channels,
            config.num_modes,
            config.num_modes,
            config.activation,
            config.use_layer_norm,
            config.use_residual
        ));
    }
    
    // Projection layers: hidden_channels -> 128 -> output_channels
    projection1 = std::make_unique<Linear>(config.hidden_channels, 128);
    projection2 = std::make_unique<Linear>(128, config.output_channels);
}

Tensor FNOModel::forward(const Tensor& x) {
    // x: [batch, input_channels, nx, ny]
    int batch = x.shape[0];
    int nx = x.shape[2];
    int ny = x.shape[3];
    
    // Lifting: reshape, apply linear, reshape back
    Tensor x_flat({batch * nx * ny, config.input_channels});
    for (int b = 0; b < batch; ++b) {
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                for (int c = 0; c < config.input_channels; ++c) {
                    x_flat((b * nx + i) * ny + j, c) = x(b, c, i, j);
                }
            }
        }
    }
    
    Tensor lifted_flat = lifting->forward(x_flat);
    
    lifted = Tensor({batch, config.hidden_channels, nx, ny});
    for (int b = 0; b < batch; ++b) {
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                for (int c = 0; c < config.hidden_channels; ++c) {
                    lifted(b, c, i, j) = lifted_flat((b * nx + i) * ny + j, c);
                }
            }
        }
    }
    
    // FNO layers
    layer_outputs.clear();
    Tensor h = lifted;
    for (auto& layer : fno_layers) {
        h = layer->forward(h);
        layer_outputs.push_back(h);
    }
    
    // Projection
    Tensor h_flat({batch * nx * ny, config.hidden_channels});
    for (int b = 0; b < batch; ++b) {
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                for (int c = 0; c < config.hidden_channels; ++c) {
                    h_flat((b * nx + i) * ny + j, c) = h(b, c, i, j);
                }
            }
        }
    }
    
    Tensor p1 = Activation::gelu(projection1->forward(h_flat));
    proj1_out = p1;
    Tensor p2 = projection2->forward(p1);
    
    // Reshape output
    Tensor output({batch, config.output_channels, nx, ny});
    for (int b = 0; b < batch; ++b) {
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                for (int c = 0; c < config.output_channels; ++c) {
                    output(b, c, i, j) = p2((b * nx + i) * ny + j, c);
                }
            }
        }
    }
    
    return output;
}

Tensor FNOModel::backward(const Tensor& grad_output) {
    int batch = grad_output.shape[0];
    int nx = grad_output.shape[2];
    int ny = grad_output.shape[3];
    
    // Flatten gradient
    Tensor grad_flat({batch * nx * ny, config.output_channels});
    for (int b = 0; b < batch; ++b) {
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                for (int c = 0; c < config.output_channels; ++c) {
                    grad_flat((b * nx + i) * ny + j, c) = grad_output(b, c, i, j);
                }
            }
        }
    }
    
    // Backward through projection
    Tensor grad_p2 = projection2->backward(grad_flat);
    Tensor grad_p1 = Activation::gelu_backward(proj1_out, grad_p2);
    Tensor grad_proj = projection1->backward(grad_p1);
    
    // Reshape
    Tensor grad_h({batch, config.hidden_channels, nx, ny});
    for (int b = 0; b < batch; ++b) {
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                for (int c = 0; c < config.hidden_channels; ++c) {
                    grad_h(b, c, i, j) = grad_proj((b * nx + i) * ny + j, c);
                }
            }
        }
    }
    
    // Backward through FNO layers
    for (int i = fno_layers.size() - 1; i >= 0; --i) {
        grad_h = fno_layers[i]->backward(grad_h);
    }
    
    // Backward through lifting
    Tensor grad_lifted_flat({batch * nx * ny, config.hidden_channels});
    for (int b = 0; b < batch; ++b) {
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                for (int c = 0; c < config.hidden_channels; ++c) {
                    grad_lifted_flat((b * nx + i) * ny + j, c) = grad_h(b, c, i, j);
                }
            }
        }
    }
    
    Tensor grad_input_flat = lifting->backward(grad_lifted_flat);
    
    // Reshape to input
    Tensor grad_input({batch, config.input_channels, nx, ny});
    for (int b = 0; b < batch; ++b) {
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                for (int c = 0; c < config.input_channels; ++c) {
                    grad_input(b, c, i, j) = grad_input_flat((b * nx + i) * ny + j, c);
                }
            }
        }
    }
    
    return grad_input;
}

Tensor FNOModel::predict(const Tensor& x) {
    return forward(x);
}

double FNOModel::computeLoss(const Tensor& prediction, const Tensor& target) {
    // MSE loss
    double loss = 0.0;
    for (size_t i = 0; i < prediction.data.size(); ++i) {
        double diff = prediction.data[i] - target.data[i];
        loss += diff * diff;
    }
    return loss / prediction.data.size();
}

Tensor FNOModel::computeLossGrad(const Tensor& prediction, const Tensor& target) {
    // d(MSE)/d(prediction) = 2 * (prediction - target) / n
    Tensor grad(prediction.shape);
    double scale = 2.0 / prediction.data.size();
    for (size_t i = 0; i < prediction.data.size(); ++i) {
        grad.data[i] = scale * (prediction.data[i] - target.data[i]);
    }
    return grad;
}

std::vector<Tensor*> FNOModel::parameters() {
    std::vector<Tensor*> params;
    
    params.push_back(&lifting->weight);
    if (lifting->has_bias) params.push_back(&lifting->bias_vec);
    
    for (auto& layer : fno_layers) {
        // Spectral conv weights (complex - handled separately)
        params.push_back(&layer->linear->weight);
        if (layer->linear->has_bias) params.push_back(&layer->linear->bias_vec);
        if (layer->use_layer_norm) {
            params.push_back(&layer->norm->gamma);
            params.push_back(&layer->norm->beta);
        }
    }
    
    params.push_back(&projection1->weight);
    if (projection1->has_bias) params.push_back(&projection1->bias_vec);
    params.push_back(&projection2->weight);
    if (projection2->has_bias) params.push_back(&projection2->bias_vec);
    
    return params;
}

std::vector<Tensor*> FNOModel::gradients() {
    std::vector<Tensor*> grads;
    
    grads.push_back(&lifting->grad_weight);
    if (lifting->has_bias) grads.push_back(&lifting->grad_bias);
    
    for (auto& layer : fno_layers) {
        grads.push_back(&layer->linear->grad_weight);
        if (layer->linear->has_bias) grads.push_back(&layer->linear->grad_bias);
        if (layer->use_layer_norm) {
            grads.push_back(&layer->norm->grad_gamma);
            grads.push_back(&layer->norm->grad_beta);
        }
    }
    
    grads.push_back(&projection1->grad_weight);
    if (projection1->has_bias) grads.push_back(&projection1->grad_bias);
    grads.push_back(&projection2->grad_weight);
    if (projection2->has_bias) grads.push_back(&projection2->grad_bias);
    
    return grads;
}

void FNOModel::zeroGrad() {
    for (auto* grad : gradients()) {
        grad->zeros();
    }
}

void FNOModel::save(const std::string& path) const {
    std::ofstream file(path, std::ios::binary);
    
    // Save config
    file.write(reinterpret_cast<const char*>(&config.input_channels), sizeof(int));
    file.write(reinterpret_cast<const char*>(&config.output_channels), sizeof(int));
    file.write(reinterpret_cast<const char*>(&config.hidden_channels), sizeof(int));
    file.write(reinterpret_cast<const char*>(&config.num_layers), sizeof(int));
    file.write(reinterpret_cast<const char*>(&config.num_modes), sizeof(int));
    
    // Save parameters
    lifting->weight.save(path + ".lifting");
    for (int i = 0; i < config.num_layers; ++i) {
        fno_layers[i]->linear->weight.save(path + ".layer" + std::to_string(i));
    }
    projection1->weight.save(path + ".proj1");
    projection2->weight.save(path + ".proj2");
}

void FNOModel::load(const std::string& path) {
    // Load parameters
    lifting->weight.load(path + ".lifting");
    for (int i = 0; i < config.num_layers; ++i) {
        fno_layers[i]->linear->weight.load(path + ".layer" + std::to_string(i));
    }
    projection1->weight.load(path + ".proj1");
    projection2->weight.load(path + ".proj2");
}

size_t FNOModel::numParameters() const {
    size_t count = 0;
    count += lifting->weight.numel();
    for (const auto& layer : fno_layers) {
        count += layer->spectral_conv->weights.numel() * 2;  // Complex
        count += layer->linear->weight.numel();
    }
    count += projection1->weight.numel();
    count += projection2->weight.numel();
    return count;
}

void FNOModel::summary() const {
    std::cout << "FNO Model Summary" << std::endl;
    std::cout << "=================" << std::endl;
    std::cout << "Input channels: " << config.input_channels << std::endl;
    std::cout << "Output channels: " << config.output_channels << std::endl;
    std::cout << "Hidden channels: " << config.hidden_channels << std::endl;
    std::cout << "Number of layers: " << config.num_layers << std::endl;
    std::cout << "Fourier modes: " << config.num_modes << std::endl;
    std::cout << "Total parameters: " << numParameters() << std::endl;
}

// =============================================================================
// AdamOptimizer Implementation
// =============================================================================

AdamOptimizer::AdamOptimizer(std::vector<Tensor*> params, double lr,
                             double beta1, double beta2, double eps, double weight_decay)
    : params(params), lr(lr), beta1(beta1), beta2(beta2), eps(eps), weight_decay(weight_decay) {
    
    // Initialize moments
    for (auto* p : params) {
        m.push_back(Tensor(p->shape, 0.0));
        v.push_back(Tensor(p->shape, 0.0));
    }
}

void AdamOptimizer::step(std::vector<Tensor*> grads) {
    t++;
    
    double bias_correction1 = 1.0 - std::pow(beta1, t);
    double bias_correction2 = 1.0 - std::pow(beta2, t);
    
    for (size_t i = 0; i < params.size(); ++i) {
        Tensor* p = params[i];
        Tensor* g = grads[i];
        
        // Weight decay
        if (weight_decay > 0) {
            for (size_t j = 0; j < p->data.size(); ++j) {
                g->data[j] += weight_decay * p->data[j];
            }
        }
        
        // Update moments
        for (size_t j = 0; j < p->data.size(); ++j) {
            m[i].data[j] = beta1 * m[i].data[j] + (1 - beta1) * g->data[j];
            v[i].data[j] = beta2 * v[i].data[j] + (1 - beta2) * g->data[j] * g->data[j];
        }
        
        // Bias-corrected step
        for (size_t j = 0; j < p->data.size(); ++j) {
            double m_hat = m[i].data[j] / bias_correction1;
            double v_hat = v[i].data[j] / bias_correction2;
            p->data[j] -= lr * m_hat / (std::sqrt(v_hat) + eps);
        }
    }
}

void AdamOptimizer::zeroGrad() {
    // Gradients are zeroed in model
}

// =============================================================================
// LRScheduler Implementation
// =============================================================================

LRScheduler::LRScheduler(AdamOptimizer& optimizer, Type type, double initial_lr)
    : optimizer(optimizer), type(type), initial_lr(initial_lr), current_lr(initial_lr) {}

void LRScheduler::step(double metric) {
    epoch++;
    
    switch (type) {
        case Type::CONSTANT:
            break;
            
        case Type::STEP:
            if (epoch % step_size == 0) {
                current_lr *= gamma;
            }
            break;
            
        case Type::EXPONENTIAL:
            current_lr = initial_lr * std::pow(gamma, epoch);
            break;
            
        case Type::COSINE_ANNEALING:
            current_lr = eta_min + 0.5 * (initial_lr - eta_min) * 
                        (1 + std::cos(M_PI * epoch / T_max));
            break;
            
        case Type::REDUCE_ON_PLATEAU:
            if (metric < best_metric) {
                best_metric = metric;
                patience_counter = 0;
            } else {
                patience_counter++;
                if (patience_counter >= patience) {
                    current_lr *= gamma;
                    patience_counter = 0;
                }
            }
            break;
    }
    
    optimizer.setLearningRate(current_lr);
}

double LRScheduler::getCurrentLR() const {
    return current_lr;
}

void LRScheduler::setStepSchedule(int step_size, double gamma) {
    this->step_size = step_size;
    this->gamma = gamma;
}

void LRScheduler::setExponentialSchedule(double gamma) {
    this->gamma = gamma;
}

void LRScheduler::setCosineSchedule(int T_max, double eta_min) {
    this->T_max = T_max;
    this->eta_min = eta_min;
}

void LRScheduler::setReduceOnPlateau(double factor, int patience) {
    this->gamma = factor;
    this->patience = patience;
}

// =============================================================================
// FNOTrainer Implementation
// =============================================================================

FNOTrainer::FNOTrainer(FNOModel& model, const FNOConfig& config)
    : model(model), config(config), rng(std::random_device{}()) {
    
    optimizer = std::make_unique<AdamOptimizer>(
        model.parameters(),
        config.learning_rate,
        0.9, 0.999, 1e-8,
        config.weight_decay
    );
    
    scheduler = std::make_unique<LRScheduler>(
        *optimizer,
        LRScheduler::Type::COSINE_ANNEALING,
        config.learning_rate
    );
    scheduler->setCosineSchedule(config.epochs);
}

void FNOTrainer::train(const FNODataset& train_data, const FNODataset& val_data) {
    std::cout << "Training FNO model..." << std::endl;
    std::cout << "Training samples: " << train_data.size() << std::endl;
    std::cout << "Validation samples: " << val_data.size() << std::endl;
    
    for (int epoch = 0; epoch < config.epochs; ++epoch) {
        auto start = std::chrono::high_resolution_clock::now();
        
        double train_loss = trainEpoch(train_data);
        double val_loss = validate(val_data);
        
        auto end = std::chrono::high_resolution_clock::now();
        double epoch_time = std::chrono::duration<double>(end - start).count();
        
        history.train_losses.push_back(train_loss);
        history.val_losses.push_back(val_loss);
        history.learning_rates.push_back(scheduler->getCurrentLR());
        
        std::cout << "Epoch " << std::setw(3) << epoch + 1 << "/" << config.epochs
                  << " - train_loss: " << std::scientific << std::setprecision(4) << train_loss
                  << " - val_loss: " << val_loss
                  << " - lr: " << scheduler->getCurrentLR()
                  << " - time: " << std::fixed << std::setprecision(1) << epoch_time << "s"
                  << std::endl;
        
        if (callback) {
            callback(epoch, train_loss, val_loss);
        }
        
        // Early stopping
        if (use_early_stopping) {
            if (val_loss < best_val_loss - es_min_delta) {
                best_val_loss = val_loss;
                es_counter = 0;
            } else {
                es_counter++;
                if (es_counter >= es_patience) {
                    std::cout << "Early stopping at epoch " << epoch + 1 << std::endl;
                    break;
                }
            }
        }
        
        scheduler->step(val_loss);
    }
}

double FNOTrainer::trainEpoch(const FNODataset& data) {
    double total_loss = 0.0;
    int num_batches = 0;
    
    // Shuffle indices
    std::vector<size_t> indices(data.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), rng);
    
    // Process in batches
    for (size_t batch_start = 0; batch_start < data.size(); batch_start += config.batch_size) {
        size_t batch_end = std::min(batch_start + config.batch_size, data.size());
        int batch_size = batch_end - batch_start;
        
        // Build batch tensors
        Tensor input(data.input_shape);
        input.shape[0] = batch_size;
        input.data.resize(batch_size * input.numel() / data.input_shape[0]);
        
        Tensor target(data.output_shape);
        target.shape[0] = batch_size;
        target.data.resize(batch_size * target.numel() / data.output_shape[0]);
        
        for (int b = 0; b < batch_size; ++b) {
            const auto& sample = data.samples[indices[batch_start + b]];
            size_t input_size = sample.input.size();
            size_t output_size = sample.output.size();
            
            std::copy(sample.input.begin(), sample.input.end(),
                     input.data.begin() + b * input_size);
            std::copy(sample.output.begin(), sample.output.end(),
                     target.data.begin() + b * output_size);
        }
        
        // Forward
        model.zeroGrad();
        Tensor prediction = model.forward(input);
        
        // Loss
        double loss = model.computeLoss(prediction, target);
        total_loss += loss;
        num_batches++;
        
        // Backward
        Tensor grad = model.computeLossGrad(prediction, target);
        model.backward(grad);
        
        // Update
        optimizer->step(model.gradients());
    }
    
    return total_loss / num_batches;
}

double FNOTrainer::validate(const FNODataset& data) {
    double total_loss = 0.0;
    int num_batches = 0;
    
    for (size_t batch_start = 0; batch_start < data.size(); batch_start += config.batch_size) {
        size_t batch_end = std::min(batch_start + config.batch_size, data.size());
        int batch_size = batch_end - batch_start;
        
        Tensor input(data.input_shape);
        input.shape[0] = batch_size;
        input.data.resize(batch_size * input.numel() / data.input_shape[0]);
        
        Tensor target(data.output_shape);
        target.shape[0] = batch_size;
        target.data.resize(batch_size * target.numel() / data.output_shape[0]);
        
        for (int b = 0; b < batch_size; ++b) {
            const auto& sample = data.samples[batch_start + b];
            size_t input_size = sample.input.size();
            size_t output_size = sample.output.size();
            
            std::copy(sample.input.begin(), sample.input.end(),
                     input.data.begin() + b * input_size);
            std::copy(sample.output.begin(), sample.output.end(),
                     target.data.begin() + b * output_size);
        }
        
        Tensor prediction = model.predict(input);
        double loss = model.computeLoss(prediction, target);
        total_loss += loss;
        num_batches++;
    }
    
    return total_loss / num_batches;
}

void FNOTrainer::setEarlyStopping(int patience, double min_delta) {
    use_early_stopping = true;
    es_patience = patience;
    es_min_delta = min_delta;
}

// =============================================================================
// FNOSolver Implementation
// =============================================================================

FNOSolver::FNOSolver(const FNOConfig& config) : config(config) {
    model = std::make_unique<FNOModel>(config);
    trainer = std::make_unique<FNOTrainer>(*model, config);
    data_gen = std::make_unique<FNODataGenerator>(config);
    
    if (config.load_pretrained && !config.model_path.empty()) {
        loadModel(config.model_path);
        model_trained = true;
    }
}

void FNOSolver::setMethod(SolverMethod method) {
    this->method = method;
}

void FNOSolver::setNumericalSolver(NumericalSolverFn solver) {
    numerical_solver = solver;
}

Tensor FNOSolver::solve(const Tensor& initial_condition, double t_final) {
    auto start = std::chrono::high_resolution_clock::now();
    
    Tensor result;
    
    switch (method) {
        case SolverMethod::NUMERICAL:
            if (!numerical_solver) {
                throw std::runtime_error("Numerical solver not set");
            }
            result = numerical_solver(initial_condition, t_final);
            stats.numerical_calls++;
            break;
            
        case SolverMethod::FNO:
            if (!model_trained) {
                throw std::runtime_error("FNO model not trained");
            }
            result = model->predict(initial_condition);
            stats.fno_calls++;
            break;
            
        case SolverMethod::HYBRID_FNO_REFINE:
            result = solveHybridRefine(initial_condition, t_final);
            break;
            
        case SolverMethod::HYBRID_COARSE_FNO:
            result = solveHybridCoarse(initial_condition, t_final);
            break;
            
        case SolverMethod::ENSEMBLE:
            result = solveEnsemble(initial_condition, t_final);
            break;
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    stats.total_time += std::chrono::duration<double>(end - start).count();
    
    return result;
}

Tensor FNOSolver::step(const Tensor& current_state, double dt) {
    // Single time step
    if (method == SolverMethod::FNO || method == SolverMethod::HYBRID_FNO_REFINE) {
        return model->predict(current_state);
    } else {
        return numerical_solver(current_state, dt);
    }
}

void FNOSolver::trainFromData(const FNODataset& dataset) {
    FNODataset train_data, val_data;
    dataset.split(config.validation_split, val_data, train_data);
    
    trainer->train(train_data, val_data);
    model_trained = true;
}

void FNOSolver::trainOnline(const Tensor& input, const Tensor& target) {
    // Single sample online training
    model->zeroGrad();
    Tensor prediction = model->forward(input);
    Tensor grad = model->computeLossGrad(prediction, target);
    model->backward(grad);
    
    // Need to manually call optimizer
    // This is simplified - real implementation would use proper optimizer
}

void FNOSolver::loadModel(const std::string& path) {
    model->load(path);
    model_trained = true;
}

void FNOSolver::saveModel(const std::string& path) const {
    model->save(path);
}

Tensor FNOSolver::solveHybridRefine(const Tensor& initial_condition, double t_final) {
    // FNO prediction
    auto fno_start = std::chrono::high_resolution_clock::now();
    Tensor fno_result = model->predict(initial_condition);
    auto fno_end = std::chrono::high_resolution_clock::now();
    stats.fno_time += std::chrono::duration<double>(fno_end - fno_start).count();
    stats.fno_calls++;
    
    // Numerical refinement starting from FNO prediction
    auto num_start = std::chrono::high_resolution_clock::now();
    Tensor refined = fno_result;
    double dt = t_final / config.refinement_iterations;
    for (int i = 0; i < config.refinement_iterations; ++i) {
        refined = numerical_solver(refined, dt);
    }
    auto num_end = std::chrono::high_resolution_clock::now();
    stats.numerical_time += std::chrono::duration<double>(num_end - num_start).count();
    stats.numerical_calls++;
    
    return refined;
}

Tensor FNOSolver::solveHybridCoarse(const Tensor& initial_condition, double t_final) {
    // Coarse numerical solution
    auto num_start = std::chrono::high_resolution_clock::now();
    Tensor coarse = numerical_solver(initial_condition, t_final);
    auto num_end = std::chrono::high_resolution_clock::now();
    stats.numerical_time += std::chrono::duration<double>(num_end - num_start).count();
    stats.numerical_calls++;
    
    // FNO upscaling/correction
    auto fno_start = std::chrono::high_resolution_clock::now();
    Tensor correction = model->predict(coarse);
    auto fno_end = std::chrono::high_resolution_clock::now();
    stats.fno_time += std::chrono::duration<double>(fno_end - fno_start).count();
    stats.fno_calls++;
    
    // Add correction
    return coarse + correction * 0.5;
}

Tensor FNOSolver::solveEnsemble(const Tensor& initial_condition, double t_final) {
    // FNO prediction
    auto fno_start = std::chrono::high_resolution_clock::now();
    Tensor fno_result = model->predict(initial_condition);
    auto fno_end = std::chrono::high_resolution_clock::now();
    stats.fno_time += std::chrono::duration<double>(fno_end - fno_start).count();
    stats.fno_calls++;
    
    // Numerical solution
    auto num_start = std::chrono::high_resolution_clock::now();
    Tensor num_result = numerical_solver(initial_condition, t_final);
    auto num_end = std::chrono::high_resolution_clock::now();
    stats.numerical_time += std::chrono::duration<double>(num_end - num_start).count();
    stats.numerical_calls++;
    
    // Weighted average
    Tensor result(fno_result.shape);
    for (size_t i = 0; i < result.data.size(); ++i) {
        result.data[i] = config.fno_weight * fno_result.data[i] + 
                        (1 - config.fno_weight) * num_result.data[i];
    }
    
    // Estimate error
    stats.avg_error_estimate = estimateError(fno_result, num_result);
    
    return result;
}

double FNOSolver::estimateError(const Tensor& fno_result, const Tensor& numerical_result) {
    double error = 0.0;
    double norm = 0.0;
    for (size_t i = 0; i < fno_result.data.size(); ++i) {
        double diff = fno_result.data[i] - numerical_result.data[i];
        error += diff * diff;
        norm += numerical_result.data[i] * numerical_result.data[i];
    }
    return std::sqrt(error / (norm + 1e-10));
}

// =============================================================================
// Data Generation
// =============================================================================

FNODataGenerator::FNODataGenerator(const FNOConfig& config) 
    : config(config), rng(std::random_device{}()) {}

void FNODataGenerator::generateFromSimulation(
    std::function<Tensor(const Tensor&)> numerical_solver,
    std::function<Tensor()> initial_condition_sampler,
    int num_samples) {
    
    std::cout << "Generating " << num_samples << " training samples..." << std::endl;
    
    for (int i = 0; i < num_samples; ++i) {
        Tensor ic = initial_condition_sampler();
        Tensor output = numerical_solver(ic);
        
        addSimulationResult(ic, output);
        
        if ((i + 1) % 100 == 0) {
            std::cout << "  Generated " << (i + 1) << "/" << num_samples << " samples" << std::endl;
        }
    }
    
    // Set dataset shapes
    if (!dataset.samples.empty()) {
        auto& sample = dataset.samples[0];
        // Infer shapes from data
        int total_input = sample.input.size();
        int total_output = sample.output.size();
        
        // Assume 2D for now
        int channels_in = config.input_channels;
        int channels_out = config.output_channels;
        int nx = config.resolution.empty() ? 64 : config.resolution[0];
        int ny = config.resolution.size() > 1 ? config.resolution[1] : nx;
        
        dataset.input_shape = {1, channels_in, nx, ny};
        dataset.output_shape = {1, channels_out, nx, ny};
    }
}

void FNODataGenerator::addSimulationResult(const Tensor& input, const Tensor& output) {
    FNOSample sample;
    sample.input = input.data;
    sample.output = output.data;
    
    if (augment_flip || augment_rotate || augment_noise) {
        augmentSample(sample);
    }
    
    dataset.addSample(sample);
}

void FNODataGenerator::addSimulationResult(const Tensor& input, const Tensor& output,
                                           const std::vector<double>& params) {
    FNOSample sample;
    sample.input = input.data;
    sample.output = output.data;
    sample.params = params;
    
    if (augment_flip || augment_rotate || augment_noise) {
        augmentSample(sample);
    }
    
    dataset.addSample(sample);
}

void FNODataGenerator::saveDataset(const std::string& path) const {
    std::ofstream file(path, std::ios::binary);
    
    size_t n_samples = dataset.samples.size();
    file.write(reinterpret_cast<const char*>(&n_samples), sizeof(size_t));
    
    for (const auto& sample : dataset.samples) {
        size_t input_size = sample.input.size();
        size_t output_size = sample.output.size();
        
        file.write(reinterpret_cast<const char*>(&input_size), sizeof(size_t));
        file.write(reinterpret_cast<const char*>(sample.input.data()), 
                  input_size * sizeof(double));
        
        file.write(reinterpret_cast<const char*>(&output_size), sizeof(size_t));
        file.write(reinterpret_cast<const char*>(sample.output.data()), 
                  output_size * sizeof(double));
    }
}

void FNODataGenerator::loadDataset(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    
    size_t n_samples;
    file.read(reinterpret_cast<char*>(&n_samples), sizeof(size_t));
    
    dataset.samples.clear();
    dataset.samples.reserve(n_samples);
    
    for (size_t i = 0; i < n_samples; ++i) {
        FNOSample sample;
        
        size_t input_size, output_size;
        file.read(reinterpret_cast<char*>(&input_size), sizeof(size_t));
        sample.input.resize(input_size);
        file.read(reinterpret_cast<char*>(sample.input.data()), 
                 input_size * sizeof(double));
        
        file.read(reinterpret_cast<char*>(&output_size), sizeof(size_t));
        sample.output.resize(output_size);
        file.read(reinterpret_cast<char*>(sample.output.data()), 
                 output_size * sizeof(double));
        
        dataset.addSample(sample);
    }
}

void FNODataGenerator::enableAugmentation(bool flip, bool rotate, bool noise) {
    augment_flip = flip;
    augment_rotate = rotate;
    augment_noise = noise;
}

void FNODataGenerator::augmentSample(FNOSample& sample) {
    std::uniform_real_distribution<double> uniform(0, 1);
    
    // Add noise
    if (augment_noise) {
        std::normal_distribution<double> noise(0, noise_std);
        for (double& x : sample.input) {
            x += noise(rng);
        }
    }
    
    // Flip and rotate would require knowing the spatial dimensions
    // Simplified for now
}

// =============================================================================
// Dataset Methods
// =============================================================================

void FNODataset::shuffle(std::mt19937& rng) {
    std::shuffle(samples.begin(), samples.end(), rng);
}

void FNODataset::split(double ratio, FNODataset& set1, FNODataset& set2) const {
    size_t split_idx = static_cast<size_t>(samples.size() * ratio);
    
    set1.samples.assign(samples.begin(), samples.begin() + split_idx);
    set2.samples.assign(samples.begin() + split_idx, samples.end());
    
    set1.input_shape = input_shape;
    set1.output_shape = output_shape;
    set2.input_shape = input_shape;
    set2.output_shape = output_shape;
}

// =============================================================================
// Configuration Parsing
// =============================================================================

FNOConfig parseFNOConfig(const std::map<std::string, std::string>& config_map) {
    FNOConfig config;
    
    auto getInt = [&](const std::string& key, int default_val) {
        auto it = config_map.find(key);
        return it != config_map.end() ? std::stoi(it->second) : default_val;
    };
    
    auto getDouble = [&](const std::string& key, double default_val) {
        auto it = config_map.find(key);
        return it != config_map.end() ? std::stod(it->second) : default_val;
    };
    
    auto getString = [&](const std::string& key, const std::string& default_val) {
        auto it = config_map.find(key);
        return it != config_map.end() ? it->second : default_val;
    };
    
    auto getBool = [&](const std::string& key, bool default_val) {
        auto it = config_map.find(key);
        if (it == config_map.end()) return default_val;
        return it->second == "true" || it->second == "1";
    };
    
    config.input_channels = getInt("fno_input_channels", 3);
    config.output_channels = getInt("fno_output_channels", 3);
    config.hidden_channels = getInt("fno_hidden_channels", 64);
    config.num_layers = getInt("fno_num_layers", 4);
    config.num_modes = getInt("fno_num_modes", 16);
    config.spatial_dim = getInt("fno_spatial_dim", 2);
    
    config.learning_rate = getDouble("fno_learning_rate", 1e-3);
    config.weight_decay = getDouble("fno_weight_decay", 1e-4);
    config.batch_size = getInt("fno_batch_size", 16);
    config.epochs = getInt("fno_epochs", 100);
    config.validation_split = getDouble("fno_validation_split", 0.1);
    
    config.activation = getString("fno_activation", "gelu");
    config.use_layer_norm = getBool("fno_use_layer_norm", true);
    config.use_residual = getBool("fno_use_residual", true);
    
    config.use_gpu = getBool("fno_use_gpu", true);
    config.gpu_device_id = getInt("fno_gpu_device", 0);
    
    config.model_path = getString("fno_model_path", "");
    config.load_pretrained = getBool("fno_load_pretrained", false);
    
    config.fno_weight = getDouble("fno_ensemble_weight", 0.8);
    config.refinement_iterations = getInt("fno_refinement_iterations", 3);
    
    return config;
}

SolverMethod parseSolverMethod(const std::string& method_str) {
    if (method_str == "numerical" || method_str == "NUMERICAL") {
        return SolverMethod::NUMERICAL;
    } else if (method_str == "fno" || method_str == "FNO") {
        return SolverMethod::FNO;
    } else if (method_str == "hybrid_refine" || method_str == "HYBRID_FNO_REFINE") {
        return SolverMethod::HYBRID_FNO_REFINE;
    } else if (method_str == "hybrid_coarse" || method_str == "HYBRID_COARSE_FNO") {
        return SolverMethod::HYBRID_COARSE_FNO;
    } else if (method_str == "ensemble" || method_str == "ENSEMBLE") {
        return SolverMethod::ENSEMBLE;
    }
    return SolverMethod::NUMERICAL;  // Default
}

std::string solverMethodToString(SolverMethod method) {
    switch (method) {
        case SolverMethod::NUMERICAL: return "NUMERICAL";
        case SolverMethod::FNO: return "FNO";
        case SolverMethod::HYBRID_FNO_REFINE: return "HYBRID_FNO_REFINE";
        case SolverMethod::HYBRID_COARSE_FNO: return "HYBRID_COARSE_FNO";
        case SolverMethod::ENSEMBLE: return "ENSEMBLE";
    }
    return "NUMERICAL";
}

// =============================================================================
// Utility Functions
// =============================================================================

namespace FNOUtils {

Tensor createGrid2D(int nx, int ny, double Lx, double Ly) {
    Tensor grid({2, nx, ny});
    
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            grid(0, i, j) = (double)i / (nx - 1) * Lx;
            grid(1, i, j) = (double)j / (ny - 1) * Ly;
        }
    }
    
    return grid;
}

Tensor createGrid3D(int nx, int ny, int nz, double Lx, double Ly, double Lz) {
    Tensor grid({3, nx, ny, nz});
    
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                grid.data[((0 * nx + i) * ny + j) * nz + k] = (double)i / (nx - 1) * Lx;
                grid.data[((1 * nx + i) * ny + j) * nz + k] = (double)j / (ny - 1) * Ly;
                grid.data[((2 * nx + i) * ny + j) * nz + k] = (double)k / (nz - 1) * Lz;
            }
        }
    }
    
    return grid;
}

void Normalizer::fit(const FNODataset& data) {
    if (data.samples.empty()) return;
    
    size_t n = data.samples[0].input.size();
    mean = Tensor({(int)n}, 0.0);
    std = Tensor({(int)n}, 0.0);
    
    // Compute mean
    for (const auto& sample : data.samples) {
        for (size_t i = 0; i < n; ++i) {
            mean.data[i] += sample.input[i];
        }
    }
    for (size_t i = 0; i < n; ++i) {
        mean.data[i] /= data.samples.size();
    }
    
    // Compute std
    for (const auto& sample : data.samples) {
        for (size_t i = 0; i < n; ++i) {
            double diff = sample.input[i] - mean.data[i];
            std.data[i] += diff * diff;
        }
    }
    for (size_t i = 0; i < n; ++i) {
        std.data[i] = std::sqrt(std.data[i] / data.samples.size()) + 1e-8;
    }
}

Tensor Normalizer::normalize(const Tensor& x) const {
    Tensor result(x.shape);
    for (size_t i = 0; i < x.data.size(); ++i) {
        result.data[i] = (x.data[i] - mean.data[i % mean.data.size()]) / 
                         std.data[i % std.data.size()];
    }
    return result;
}

Tensor Normalizer::denormalize(const Tensor& x) const {
    Tensor result(x.shape);
    for (size_t i = 0; i < x.data.size(); ++i) {
        result.data[i] = x.data[i] * std.data[i % std.data.size()] + 
                         mean.data[i % mean.data.size()];
    }
    return result;
}

} // namespace FNOUtils

} // namespace ML
} // namespace FSRM
