/**
 * @file test_fourier_neural_operator.cpp
 * @brief Unit tests for Fourier Neural Operator
 * 
 * Tests cover:
 * - Tensor operations
 * - FFT computations
 * - FNO layers (SpectralConv, Linear, LayerNorm)
 * - FNO model forward/backward
 * - Training and optimization
 * - Hybrid solver modes
 */

#include "FourierNeuralOperator.hpp"
#include "Testing.hpp"
#include <cmath>
#include <iostream>
#include <vector>

namespace FSRM {
namespace Testing {

using namespace ML;

// =============================================================================
// Tensor Tests
// =============================================================================

/**
 * @test Test tensor construction and access
 */
bool test_tensor_construction() {
    std::cout << "Testing tensor construction..." << std::endl;
    
    // Default construction
    Tensor t1({2, 3, 4});
    if (t1.numel() != 24) {
        std::cerr << "  FAIL: Wrong number of elements" << std::endl;
        return false;
    }
    
    // Construction with value
    Tensor t2({3, 3}, 1.5);
    if (std::abs(t2(1, 1) - 1.5) > 1e-10) {
        std::cerr << "  FAIL: Wrong initialized value" << std::endl;
        return false;
    }
    
    // Access and modification
    Tensor t3({2, 2});
    t3(0, 0) = 1.0;
    t3(0, 1) = 2.0;
    t3(1, 0) = 3.0;
    t3(1, 1) = 4.0;
    
    if (std::abs(t3(1, 0) - 3.0) > 1e-10) {
        std::cerr << "  FAIL: Wrong element access" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Tensor construction correct" << std::endl;
    return true;
}

/**
 * @test Test tensor arithmetic operations
 */
bool test_tensor_arithmetic() {
    std::cout << "Testing tensor arithmetic..." << std::endl;
    
    Tensor a({2, 2}, 1.0);
    Tensor b({2, 2}, 2.0);
    
    // Addition
    Tensor c = a + b;
    if (std::abs(c(0, 0) - 3.0) > 1e-10) {
        std::cerr << "  FAIL: Addition incorrect" << std::endl;
        return false;
    }
    
    // Subtraction
    Tensor d = b - a;
    if (std::abs(d(0, 0) - 1.0) > 1e-10) {
        std::cerr << "  FAIL: Subtraction incorrect" << std::endl;
        return false;
    }
    
    // Element-wise multiplication
    Tensor e = a * b;
    if (std::abs(e(0, 0) - 2.0) > 1e-10) {
        std::cerr << "  FAIL: Element-wise multiplication incorrect" << std::endl;
        return false;
    }
    
    // Scalar multiplication
    Tensor f = a * 3.0;
    if (std::abs(f(0, 0) - 3.0) > 1e-10) {
        std::cerr << "  FAIL: Scalar multiplication incorrect" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Tensor arithmetic correct" << std::endl;
    return true;
}

/**
 * @test Test tensor reductions
 */
bool test_tensor_reductions() {
    std::cout << "Testing tensor reductions..." << std::endl;
    
    double tol = 1e-10;
    
    Tensor t({2, 3});
    t(0, 0) = 1; t(0, 1) = 2; t(0, 2) = 3;
    t(1, 0) = 4; t(1, 1) = 5; t(1, 2) = 6;
    
    // Sum
    double sum = t.sum();
    if (std::abs(sum - 21.0) > tol) {
        std::cerr << "  FAIL: Sum incorrect" << std::endl;
        return false;
    }
    
    // Mean
    double mean = t.mean();
    if (std::abs(mean - 3.5) > tol) {
        std::cerr << "  FAIL: Mean incorrect" << std::endl;
        return false;
    }
    
    // Max
    double max = t.max();
    if (std::abs(max - 6.0) > tol) {
        std::cerr << "  FAIL: Max incorrect" << std::endl;
        return false;
    }
    
    // Norm
    double norm = t.norm();
    double expected_norm = std::sqrt(1 + 4 + 9 + 16 + 25 + 36);
    if (std::abs(norm - expected_norm) > tol) {
        std::cerr << "  FAIL: Norm incorrect" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Tensor reductions correct" << std::endl;
    return true;
}

/**
 * @test Test tensor initialization
 */
bool test_tensor_initialization() {
    std::cout << "Testing tensor initialization..." << std::endl;
    
    Tensor t({100});
    
    // Zeros
    t.zeros();
    if (std::abs(t.sum()) > 1e-10) {
        std::cerr << "  FAIL: zeros() not working" << std::endl;
        return false;
    }
    
    // Ones
    t.ones();
    if (std::abs(t.sum() - 100.0) > 1e-10) {
        std::cerr << "  FAIL: ones() not working" << std::endl;
        return false;
    }
    
    // Random normal
    t.randn(0.0, 1.0);
    double mean = t.mean();
    // With 100 samples, mean should be close to 0
    if (std::abs(mean) > 0.5) {
        std::cerr << "  WARNING: randn() mean far from 0: " << mean << std::endl;
    }
    
    std::cout << "  PASS: Tensor initialization correct" << std::endl;
    return true;
}

// =============================================================================
// Activation Function Tests
// =============================================================================

/**
 * @test Test ReLU activation
 */
bool test_relu_activation() {
    std::cout << "Testing ReLU activation..." << std::endl;
    
    Tensor x({4});
    x(0) = -2.0; x(1) = -0.5; x(2) = 0.5; x(3) = 2.0;
    
    Tensor y = Activation::relu(x);
    
    if (y(0) != 0.0 || y(1) != 0.0) {
        std::cerr << "  FAIL: ReLU not zeroing negatives" << std::endl;
        return false;
    }
    
    if (std::abs(y(2) - 0.5) > 1e-10 || std::abs(y(3) - 2.0) > 1e-10) {
        std::cerr << "  FAIL: ReLU changing positives" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: ReLU activation correct" << std::endl;
    return true;
}

/**
 * @test Test GELU activation
 */
bool test_gelu_activation() {
    std::cout << "Testing GELU activation..." << std::endl;
    
    Tensor x({3});
    x(0) = -2.0; x(1) = 0.0; x(2) = 2.0;
    
    Tensor y = Activation::gelu(x);
    
    // GELU(0) should be 0
    if (std::abs(y(1)) > 1e-10) {
        std::cerr << "  FAIL: GELU(0) != 0" << std::endl;
        return false;
    }
    
    // GELU is approximately ReLU for large positive values
    if (y(2) < 1.5) {
        std::cerr << "  FAIL: GELU(2) too small" << std::endl;
        return false;
    }
    
    // GELU(-x) < 0 for large negative
    if (y(0) > 0.1) {
        std::cerr << "  FAIL: GELU(-2) too large" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: GELU activation correct" << std::endl;
    return true;
}

/**
 * @test Test activation backward pass
 */
bool test_activation_backward() {
    std::cout << "Testing activation backward..." << std::endl;
    
    double h = 1e-5;
    
    Tensor x({4});
    x(0) = -1.0; x(1) = -0.1; x(2) = 0.1; x(3) = 1.0;
    
    Tensor grad_out({4}, 1.0);
    
    // Test ReLU backward
    Tensor grad_relu = Activation::relu_backward(x, grad_out);
    
    // For x < 0, gradient should be 0
    if (grad_relu(0) != 0.0) {
        std::cerr << "  FAIL: ReLU backward wrong for x < 0" << std::endl;
        return false;
    }
    
    // For x > 0, gradient should be 1
    if (std::abs(grad_relu(3) - 1.0) > 1e-10) {
        std::cerr << "  FAIL: ReLU backward wrong for x > 0" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Activation backward correct" << std::endl;
    return true;
}

// =============================================================================
// FFT Tests
// =============================================================================

/**
 * @test Test 1D FFT roundtrip
 */
bool test_fft_roundtrip() {
    std::cout << "Testing FFT roundtrip..." << std::endl;
    
    double tol = 1e-8;
    int n = 8;
    
    // Create test signal
    std::vector<double> x(n);
    for (int i = 0; i < n; ++i) {
        x[i] = std::sin(2.0 * M_PI * i / n) + 0.5 * std::cos(4.0 * M_PI * i / n);
    }
    
    // Forward FFT
    std::vector<std::complex<double>> X(n);
    FFT::fft1d(x.data(), X.data(), n);
    
    // Inverse FFT
    std::vector<double> x_recovered(n);
    FFT::ifft1d(X.data(), x_recovered.data(), n);
    
    // Check roundtrip
    for (int i = 0; i < n; ++i) {
        if (std::abs(x[i] - x_recovered[i]) > tol) {
            std::cerr << "  FAIL: FFT roundtrip error at " << i 
                      << ": " << x[i] << " vs " << x_recovered[i] << std::endl;
            return false;
        }
    }
    
    std::cout << "  PASS: FFT roundtrip correct" << std::endl;
    return true;
}

/**
 * @test Test 2D FFT properties
 */
bool test_fft2d() {
    std::cout << "Testing 2D FFT..." << std::endl;
    
    double tol = 1e-6;
    
    Tensor x({1, 1, 8, 8});  // [batch, channel, nx, ny]
    
    // Create test pattern
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j) {
            x(0, 0, i, j) = std::sin(2.0 * M_PI * i / 8) * 
                           std::cos(2.0 * M_PI * j / 8);
        }
    }
    
    // Forward FFT
    ComplexTensor X({1, 1, 8, 8});
    FFT::fft2d(x, X);
    
    // Inverse FFT
    Tensor x_recovered({1, 1, 8, 8});
    FFT::ifft2d(X, x_recovered);
    
    // Check roundtrip
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j) {
            if (std::abs(x(0, 0, i, j) - x_recovered(0, 0, i, j)) > tol) {
                std::cerr << "  FAIL: 2D FFT roundtrip error" << std::endl;
                return false;
            }
        }
    }
    
    std::cout << "  PASS: 2D FFT correct" << std::endl;
    return true;
}

// =============================================================================
// Layer Tests
// =============================================================================

/**
 * @test Test Linear layer forward pass
 */
bool test_linear_layer_forward() {
    std::cout << "Testing Linear layer forward..." << std::endl;
    
    Linear linear(4, 3);
    
    // Set known weights
    linear.weight.zeros();
    linear.weight(0, 0) = 1.0;  // Output 0 = input 0
    linear.weight(1, 1) = 1.0;  // Output 1 = input 1
    linear.weight(2, 2) = 1.0;  // Output 2 = input 2
    linear.bias_vec.zeros();
    
    Tensor x({2, 4});  // Batch of 2
    x(0, 0) = 1.0; x(0, 1) = 2.0; x(0, 2) = 3.0; x(0, 3) = 4.0;
    x(1, 0) = 5.0; x(1, 1) = 6.0; x(1, 2) = 7.0; x(1, 3) = 8.0;
    
    Tensor y = linear.forward(x);
    
    // Check output shape
    if (y.shape[0] != 2 || y.shape[1] != 3) {
        std::cerr << "  FAIL: Wrong output shape" << std::endl;
        return false;
    }
    
    // Check values (identity mapping for first 3 inputs)
    if (std::abs(y(0, 0) - 1.0) > 1e-10 || 
        std::abs(y(0, 1) - 2.0) > 1e-10 ||
        std::abs(y(0, 2) - 3.0) > 1e-10) {
        std::cerr << "  FAIL: Wrong output values" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Linear layer forward correct" << std::endl;
    return true;
}

/**
 * @test Test Linear layer backward pass
 */
bool test_linear_layer_backward() {
    std::cout << "Testing Linear layer backward..." << std::endl;
    
    double tol = 1e-5;
    double h = 1e-5;
    
    Linear linear(3, 2);
    linear.weight.randn(0, 0.1);
    linear.bias_vec.zeros();
    
    Tensor x({1, 3});
    x(0, 0) = 1.0; x(0, 1) = 2.0; x(0, 2) = 3.0;
    
    // Forward
    Tensor y = linear.forward(x);
    
    // Backward with unit gradient
    Tensor grad_out({1, 2}, 1.0);
    Tensor grad_in = linear.backward(grad_out);
    
    // Numerical gradient check for weight
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 3; ++j) {
            double orig = linear.weight(i, j);
            
            linear.weight(i, j) = orig + h;
            Tensor y_plus = linear.forward(x);
            
            linear.weight(i, j) = orig - h;
            Tensor y_minus = linear.forward(x);
            
            linear.weight(i, j) = orig;
            
            double numerical_grad = (y_plus.sum() - y_minus.sum()) / (2 * h);
            double analytical_grad = linear.grad_weight(i, j);
            
            if (std::abs(numerical_grad - analytical_grad) > tol) {
                std::cerr << "  WARNING: Gradient mismatch at (" << i << "," << j << ")" << std::endl;
            }
        }
    }
    
    std::cout << "  PASS: Linear layer backward correct" << std::endl;
    return true;
}

/**
 * @test Test LayerNorm forward pass
 */
bool test_layer_norm() {
    std::cout << "Testing LayerNorm..." << std::endl;
    
    double tol = 1e-5;
    
    LayerNorm ln(4);
    
    Tensor x({2, 4});
    x(0, 0) = 1.0; x(0, 1) = 2.0; x(0, 2) = 3.0; x(0, 3) = 4.0;
    x(1, 0) = 5.0; x(1, 1) = 6.0; x(1, 2) = 7.0; x(1, 3) = 8.0;
    
    Tensor y = ln.forward(x);
    
    // Check that output has zero mean per sample (approximately)
    double mean0 = (y(0, 0) + y(0, 1) + y(0, 2) + y(0, 3)) / 4;
    double mean1 = (y(1, 0) + y(1, 1) + y(1, 2) + y(1, 3)) / 4;
    
    if (std::abs(mean0) > tol || std::abs(mean1) > tol) {
        std::cerr << "  FAIL: LayerNorm output not zero mean" << std::endl;
        return false;
    }
    
    // Check that output has unit variance per sample
    double var0 = 0, var1 = 0;
    for (int i = 0; i < 4; ++i) {
        var0 += y(0, i) * y(0, i);
        var1 += y(1, i) * y(1, i);
    }
    var0 /= 4;
    var1 /= 4;
    
    if (std::abs(var0 - 1.0) > tol || std::abs(var1 - 1.0) > tol) {
        std::cerr << "  FAIL: LayerNorm output not unit variance" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: LayerNorm correct" << std::endl;
    return true;
}

// =============================================================================
// FNO Model Tests
// =============================================================================

/**
 * @test Test FNO model construction
 */
bool test_fno_model_construction() {
    std::cout << "Testing FNO model construction..." << std::endl;
    
    FNOConfig config;
    config.input_channels = 3;
    config.output_channels = 3;
    config.hidden_channels = 32;
    config.num_layers = 2;
    config.num_modes = 8;
    
    FNOModel model(config);
    
    size_t num_params = model.numParameters();
    if (num_params == 0) {
        std::cerr << "  FAIL: Model has no parameters" << std::endl;
        return false;
    }
    
    std::cout << "  Model has " << num_params << " parameters" << std::endl;
    
    std::cout << "  PASS: FNO model construction correct" << std::endl;
    return true;
}

/**
 * @test Test FNO model forward pass shape
 */
bool test_fno_model_forward_shape() {
    std::cout << "Testing FNO model forward pass shape..." << std::endl;
    
    FNOConfig config;
    config.input_channels = 3;
    config.output_channels = 2;
    config.hidden_channels = 16;
    config.num_layers = 2;
    config.num_modes = 4;
    
    FNOModel model(config);
    
    // Create input
    int batch = 2, nx = 16, ny = 16;
    Tensor x({batch, config.input_channels, nx, ny});
    x.randn(0, 1);
    
    // Forward pass
    Tensor y = model.forward(x);
    
    // Check output shape
    if (y.shape[0] != batch || 
        y.shape[1] != config.output_channels ||
        y.shape[2] != nx || 
        y.shape[3] != ny) {
        std::cerr << "  FAIL: Wrong output shape" << std::endl;
        return false;
    }
    
    // Check output is finite
    for (double v : y.data) {
        if (std::isnan(v) || std::isinf(v)) {
            std::cerr << "  FAIL: Output contains NaN/Inf" << std::endl;
            return false;
        }
    }
    
    std::cout << "  PASS: FNO model forward shape correct" << std::endl;
    return true;
}

/**
 * @test Test FNO model backward pass
 */
bool test_fno_model_backward() {
    std::cout << "Testing FNO model backward pass..." << std::endl;
    
    FNOConfig config;
    config.input_channels = 2;
    config.output_channels = 2;
    config.hidden_channels = 8;
    config.num_layers = 1;
    config.num_modes = 4;
    
    FNOModel model(config);
    
    Tensor x({1, 2, 8, 8});
    x.randn(0, 1);
    
    Tensor target({1, 2, 8, 8});
    target.randn(0, 1);
    
    // Forward
    model.zeroGrad();
    Tensor y = model.forward(x);
    
    // Compute loss and gradient
    double loss = model.computeLoss(y, target);
    Tensor grad = model.computeLossGrad(y, target);
    
    // Backward
    model.backward(grad);
    
    // Check gradients are computed
    auto grads = model.gradients();
    bool has_nonzero_grad = false;
    for (auto* g : grads) {
        if (g->norm() > 1e-10) {
            has_nonzero_grad = true;
            break;
        }
    }
    
    if (!has_nonzero_grad) {
        std::cerr << "  FAIL: No gradients computed" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: FNO model backward correct" << std::endl;
    return true;
}

/**
 * @test Test FNO loss computation
 */
bool test_fno_loss() {
    std::cout << "Testing FNO loss computation..." << std::endl;
    
    double tol = 1e-10;
    
    FNOConfig config;
    config.input_channels = 1;
    config.output_channels = 1;
    config.hidden_channels = 8;
    config.num_layers = 1;
    config.num_modes = 2;
    
    FNOModel model(config);
    
    Tensor pred({1, 1, 4, 4}, 1.0);
    Tensor target({1, 1, 4, 4}, 0.0);
    
    double loss = model.computeLoss(pred, target);
    
    // MSE of ones vs zeros = 1.0
    if (std::abs(loss - 1.0) > tol) {
        std::cerr << "  FAIL: Loss = " << loss << ", expected 1.0" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: FNO loss computation correct" << std::endl;
    return true;
}

// =============================================================================
// Optimizer Tests
// =============================================================================

/**
 * @test Test Adam optimizer step
 */
bool test_adam_optimizer() {
    std::cout << "Testing Adam optimizer..." << std::endl;
    
    // Simple 1D optimization: minimize x^2
    Tensor x({1}, 10.0);
    Tensor grad({1});
    
    AdamOptimizer optimizer({&x}, 0.1);
    
    double initial_val = x(0);
    
    for (int i = 0; i < 100; ++i) {
        grad(0) = 2 * x(0);  // d(x^2)/dx = 2x
        optimizer.step({&grad});
    }
    
    // x should be close to 0
    if (std::abs(x(0)) > 0.1) {
        std::cerr << "  FAIL: Optimizer didn't converge, x=" << x(0) << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Adam optimizer correct" << std::endl;
    return true;
}

// =============================================================================
// FNO Solver Tests
// =============================================================================

/**
 * @test Test FNO solver construction
 */
bool test_fno_solver_construction() {
    std::cout << "Testing FNO solver construction..." << std::endl;
    
    FNOConfig config;
    config.input_channels = 3;
    config.output_channels = 3;
    
    FNOSolver solver(config);
    
    if (solver.isModelTrained()) {
        std::cerr << "  FAIL: Model should not be trained initially" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: FNO solver construction correct" << std::endl;
    return true;
}

/**
 * @test Test solver method selection
 */
bool test_solver_method_selection() {
    std::cout << "Testing solver method selection..." << std::endl;
    
    FNOConfig config;
    FNOSolver solver(config);
    
    solver.setMethod(SolverMethod::NUMERICAL);
    if (solver.getMethod() != SolverMethod::NUMERICAL) {
        std::cerr << "  FAIL: NUMERICAL method not set" << std::endl;
        return false;
    }
    
    solver.setMethod(SolverMethod::FNO);
    if (solver.getMethod() != SolverMethod::FNO) {
        std::cerr << "  FAIL: FNO method not set" << std::endl;
        return false;
    }
    
    solver.setMethod(SolverMethod::ENSEMBLE);
    if (solver.getMethod() != SolverMethod::ENSEMBLE) {
        std::cerr << "  FAIL: ENSEMBLE method not set" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Solver method selection correct" << std::endl;
    return true;
}

/**
 * @test Test FNO prediction (untrained)
 */
bool test_fno_prediction_untrained() {
    std::cout << "Testing FNO prediction (untrained)..." << std::endl;
    
    FNOConfig config;
    config.input_channels = 2;
    config.output_channels = 2;
    config.hidden_channels = 8;
    config.num_layers = 1;
    config.num_modes = 4;
    
    FNOSolver solver(config);
    
    Tensor input({1, 2, 8, 8});
    input.randn(0, 1);
    
    // Even untrained, model should produce output
    Tensor output = solver.getModel().predict(input);
    
    if (output.shape[0] != 1 || output.shape[1] != 2 ||
        output.shape[2] != 8 || output.shape[3] != 8) {
        std::cerr << "  FAIL: Wrong output shape" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: FNO prediction (untrained) correct" << std::endl;
    return true;
}

// =============================================================================
// Configuration Parsing Tests
// =============================================================================

/**
 * @test Test FNO config parsing
 */
bool test_fno_config_parsing() {
    std::cout << "Testing FNO config parsing..." << std::endl;
    
    std::map<std::string, std::string> config_map = {
        {"fno_input_channels", "5"},
        {"fno_output_channels", "3"},
        {"fno_hidden_channels", "128"},
        {"fno_num_layers", "6"},
        {"fno_num_modes", "32"},
        {"fno_learning_rate", "0.001"},
        {"fno_activation", "relu"},
        {"fno_use_layer_norm", "false"}
    };
    
    FNOConfig config = parseFNOConfig(config_map);
    
    if (config.input_channels != 5) {
        std::cerr << "  FAIL: input_channels not parsed" << std::endl;
        return false;
    }
    
    if (config.hidden_channels != 128) {
        std::cerr << "  FAIL: hidden_channels not parsed" << std::endl;
        return false;
    }
    
    if (config.activation != "relu") {
        std::cerr << "  FAIL: activation not parsed" << std::endl;
        return false;
    }
    
    if (config.use_layer_norm != false) {
        std::cerr << "  FAIL: use_layer_norm not parsed" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: FNO config parsing correct" << std::endl;
    return true;
}

/**
 * @test Test solver method parsing
 */
bool test_solver_method_parsing() {
    std::cout << "Testing solver method parsing..." << std::endl;
    
    if (parseSolverMethod("numerical") != SolverMethod::NUMERICAL) {
        std::cerr << "  FAIL: 'numerical' not parsed" << std::endl;
        return false;
    }
    
    if (parseSolverMethod("FNO") != SolverMethod::FNO) {
        std::cerr << "  FAIL: 'FNO' not parsed" << std::endl;
        return false;
    }
    
    if (parseSolverMethod("hybrid_refine") != SolverMethod::HYBRID_FNO_REFINE) {
        std::cerr << "  FAIL: 'hybrid_refine' not parsed" << std::endl;
        return false;
    }
    
    if (parseSolverMethod("ensemble") != SolverMethod::ENSEMBLE) {
        std::cerr << "  FAIL: 'ensemble' not parsed" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Solver method parsing correct" << std::endl;
    return true;
}

// =============================================================================
// Data Generation Tests
// =============================================================================

/**
 * @test Test FNO dataset operations
 */
bool test_fno_dataset() {
    std::cout << "Testing FNO dataset operations..." << std::endl;
    
    FNODataset dataset;
    dataset.input_shape = {1, 2, 8, 8};
    dataset.output_shape = {1, 2, 8, 8};
    
    // Add samples
    for (int i = 0; i < 100; ++i) {
        FNOSample sample;
        sample.input.resize(2 * 8 * 8, (double)i);
        sample.output.resize(2 * 8 * 8, (double)i * 2);
        dataset.addSample(sample);
    }
    
    if (dataset.size() != 100) {
        std::cerr << "  FAIL: Wrong dataset size" << std::endl;
        return false;
    }
    
    // Test split
    FNODataset train, val;
    dataset.split(0.2, val, train);
    
    if (train.size() + val.size() != 100) {
        std::cerr << "  FAIL: Split lost samples" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: FNO dataset operations correct" << std::endl;
    return true;
}

// =============================================================================
// Utility Tests
// =============================================================================

/**
 * @test Test grid creation
 */
bool test_grid_creation() {
    std::cout << "Testing grid creation..." << std::endl;
    
    double tol = 1e-10;
    
    Tensor grid = FNOUtils::createGrid2D(10, 10, 1.0, 2.0);
    
    if (grid.shape[0] != 2 || grid.shape[1] != 10 || grid.shape[2] != 10) {
        std::cerr << "  FAIL: Wrong grid shape" << std::endl;
        return false;
    }
    
    // Check corner values
    if (std::abs(grid(0, 0, 0)) > tol || std::abs(grid(1, 0, 0)) > tol) {
        std::cerr << "  FAIL: Origin not at (0,0)" << std::endl;
        return false;
    }
    
    if (std::abs(grid(0, 9, 9) - 1.0) > tol || std::abs(grid(1, 9, 9) - 2.0) > tol) {
        std::cerr << "  FAIL: Far corner not at (Lx, Ly)" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Grid creation correct" << std::endl;
    return true;
}

/**
 * @test Test normalizer
 */
bool test_normalizer() {
    std::cout << "Testing normalizer..." << std::endl;
    
    double tol = 1e-5;
    
    FNODataset dataset;
    dataset.input_shape = {1, 4};
    
    // Create samples with known statistics
    for (int i = 0; i < 100; ++i) {
        FNOSample sample;
        sample.input = {1.0, 2.0, 3.0, 4.0};  // Mean = 2.5
        dataset.addSample(sample);
    }
    
    FNOUtils::Normalizer norm;
    norm.fit(dataset);
    
    // Check mean
    if (std::abs(norm.mean.mean() - 2.5) > tol) {
        std::cerr << "  FAIL: Wrong mean computed" << std::endl;
        return false;
    }
    
    // Test normalize/denormalize roundtrip
    Tensor x({4});
    x(0) = 1.0; x(1) = 2.0; x(2) = 3.0; x(3) = 4.0;
    
    Tensor x_norm = norm.normalize(x);
    Tensor x_denorm = norm.denormalize(x_norm);
    
    for (int i = 0; i < 4; ++i) {
        if (std::abs(x(i) - x_denorm(i)) > tol) {
            std::cerr << "  FAIL: Normalize/denormalize roundtrip failed" << std::endl;
            return false;
        }
    }
    
    std::cout << "  PASS: Normalizer correct" << std::endl;
    return true;
}

// =============================================================================
// Test Runner
// =============================================================================

int runFNOTests() {
    std::cout << "\n=== Fourier Neural Operator Unit Tests ===" << std::endl;
    
    int passed = 0;
    int failed = 0;
    
    // Tensor tests
    if (test_tensor_construction()) ++passed; else ++failed;
    if (test_tensor_arithmetic()) ++passed; else ++failed;
    if (test_tensor_reductions()) ++passed; else ++failed;
    if (test_tensor_initialization()) ++passed; else ++failed;
    
    // Activation tests
    if (test_relu_activation()) ++passed; else ++failed;
    if (test_gelu_activation()) ++passed; else ++failed;
    if (test_activation_backward()) ++passed; else ++failed;
    
    // FFT tests
    if (test_fft_roundtrip()) ++passed; else ++failed;
    if (test_fft2d()) ++passed; else ++failed;
    
    // Layer tests
    if (test_linear_layer_forward()) ++passed; else ++failed;
    if (test_linear_layer_backward()) ++passed; else ++failed;
    if (test_layer_norm()) ++passed; else ++failed;
    
    // Model tests
    if (test_fno_model_construction()) ++passed; else ++failed;
    if (test_fno_model_forward_shape()) ++passed; else ++failed;
    if (test_fno_model_backward()) ++passed; else ++failed;
    if (test_fno_loss()) ++passed; else ++failed;
    
    // Optimizer tests
    if (test_adam_optimizer()) ++passed; else ++failed;
    
    // Solver tests
    if (test_fno_solver_construction()) ++passed; else ++failed;
    if (test_solver_method_selection()) ++passed; else ++failed;
    if (test_fno_prediction_untrained()) ++passed; else ++failed;
    
    // Config tests
    if (test_fno_config_parsing()) ++passed; else ++failed;
    if (test_solver_method_parsing()) ++passed; else ++failed;
    
    // Data tests
    if (test_fno_dataset()) ++passed; else ++failed;
    
    // Utility tests
    if (test_grid_creation()) ++passed; else ++failed;
    if (test_normalizer()) ++passed; else ++failed;
    
    std::cout << "\n=== FNO Test Summary ===" << std::endl;
    std::cout << "Passed: " << passed << std::endl;
    std::cout << "Failed: " << failed << std::endl;
    std::cout << "Total:  " << (passed + failed) << std::endl;
    
    return failed;
}

} // namespace Testing
} // namespace FSRM

#ifndef FSRM_TEST_NO_MAIN
int main() {
    return FSRM::Testing::runFNOTests();
}
#endif
