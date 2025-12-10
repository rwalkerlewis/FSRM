/**
 * @file GPUCompute.hpp
 * @brief Unified GPU computing interface for physics and ML operations
 * 
 * This file provides a comprehensive GPU compute layer that supports:
 * - Tensor operations for neural networks
 * - FFT operations for spectral methods
 * - Physics kernels for simulation
 * - Automatic CPU fallback when GPU unavailable
 * 
 * Supports both CUDA and HIP backends with unified API.
 */

#ifndef GPU_COMPUTE_HPP
#define GPU_COMPUTE_HPP

#include "GPUManager.hpp"
#include <vector>
#include <complex>
#include <memory>
#include <functional>
#include <string>

#ifdef USE_CUDA
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cufft.h>
#include <curand.h>
#endif

#ifdef USE_HIP
#include <hip/hip_runtime.h>
#include <hipblas/hipblas.h>
#include <hipfft/hipfft.h>
#include <hiprand/hiprand.h>
#endif

namespace FSRM {
namespace GPU {

// =============================================================================
// GPU Tensor Class for ML Operations
// =============================================================================

/**
 * @brief GPU-accelerated tensor for neural network operations
 * 
 * Supports automatic memory management, CPU/GPU transfers, and
 * common tensor operations with GPU acceleration.
 */
class GPUTensor {
public:
    // Constructors
    GPUTensor();
    GPUTensor(const std::vector<int>& shape);
    GPUTensor(const std::vector<int>& shape, float value);
    GPUTensor(const std::vector<int>& shape, const float* data, bool on_gpu = false);
    ~GPUTensor();
    
    // Move semantics
    GPUTensor(GPUTensor&& other) noexcept;
    GPUTensor& operator=(GPUTensor&& other) noexcept;
    
    // Copy (creates new allocation)
    GPUTensor(const GPUTensor& other);
    GPUTensor& operator=(const GPUTensor& other);
    
    // Properties
    size_t numel() const;
    int ndim() const { return static_cast<int>(shape_.size()); }
    const std::vector<int>& shape() const { return shape_; }
    int shape(int dim) const { return shape_[dim]; }
    bool isOnGPU() const { return on_gpu_; }
    bool isContiguous() const { return true; }  // Always contiguous
    
    // Data access
    float* data() { return on_gpu_ ? d_data_ : h_data_; }
    const float* data() const { return on_gpu_ ? d_data_ : h_data_; }
    float* gpuData() { return d_data_; }
    const float* gpuData() const { return d_data_; }
    float* cpuData() { return h_data_; }
    const float* cpuData() const { return h_data_; }
    
    // Device transfer
    void toGPU();
    void toCPU();
    void synchronize();  // Wait for pending operations
    
    // Reshape (no data copy if compatible)
    GPUTensor reshape(const std::vector<int>& new_shape) const;
    GPUTensor view(const std::vector<int>& new_shape);
    GPUTensor transpose(int dim1, int dim2) const;
    GPUTensor permute(const std::vector<int>& dims) const;
    
    // Slicing
    GPUTensor slice(int dim, int start, int end) const;
    
    // Basic math operations (return new tensor)
    GPUTensor operator+(const GPUTensor& other) const;
    GPUTensor operator-(const GPUTensor& other) const;
    GPUTensor operator*(const GPUTensor& other) const;  // Element-wise
    GPUTensor operator/(const GPUTensor& other) const;  // Element-wise
    GPUTensor operator*(float scalar) const;
    GPUTensor operator+(float scalar) const;
    
    // In-place operations
    GPUTensor& operator+=(const GPUTensor& other);
    GPUTensor& operator-=(const GPUTensor& other);
    GPUTensor& operator*=(const GPUTensor& other);
    GPUTensor& operator*=(float scalar);
    GPUTensor& add_(float scalar);
    GPUTensor& mul_(float scalar);
    GPUTensor& clamp_(float min_val, float max_val);
    
    // Matrix operations
    GPUTensor matmul(const GPUTensor& other) const;      // Matrix multiplication
    GPUTensor mm(const GPUTensor& other) const;          // Alias for matmul
    GPUTensor bmm(const GPUTensor& other) const;         // Batched matmul
    GPUTensor transpose() const;                          // 2D transpose
    GPUTensor t() const { return transpose(); }
    
    // Reductions
    float sum() const;
    float mean() const;
    float max() const;
    float min() const;
    float norm(int p = 2) const;
    GPUTensor sum(int dim, bool keepdim = false) const;
    GPUTensor mean(int dim, bool keepdim = false) const;
    GPUTensor max(int dim, bool keepdim = false) const;
    GPUTensor min(int dim, bool keepdim = false) const;
    
    // Activations (in-place versions available)
    GPUTensor relu() const;
    GPUTensor sigmoid() const;
    GPUTensor tanh() const;
    GPUTensor gelu() const;
    GPUTensor silu() const;  // Swish
    GPUTensor softmax(int dim = -1) const;
    GPUTensor log_softmax(int dim = -1) const;
    
    // In-place activations
    GPUTensor& relu_();
    GPUTensor& sigmoid_();
    GPUTensor& tanh_();
    GPUTensor& gelu_();
    
    // Element-wise functions
    GPUTensor exp() const;
    GPUTensor log() const;
    GPUTensor sqrt() const;
    GPUTensor pow(float exponent) const;
    GPUTensor abs() const;
    GPUTensor neg() const;
    GPUTensor sin() const;
    GPUTensor cos() const;
    
    // Comparison
    GPUTensor eq(const GPUTensor& other) const;  // Element-wise equal
    GPUTensor ne(const GPUTensor& other) const;
    GPUTensor gt(const GPUTensor& other) const;
    GPUTensor lt(const GPUTensor& other) const;
    GPUTensor ge(const GPUTensor& other) const;
    GPUTensor le(const GPUTensor& other) const;
    
    // Initialization
    void zeros();
    void ones();
    void fill(float value);
    void uniform(float low = 0.0f, float high = 1.0f);
    void normal(float mean = 0.0f, float std = 1.0f);
    void xavier_uniform();
    void xavier_normal();
    void kaiming_uniform();
    void kaiming_normal();
    
    // Factory methods
    static GPUTensor zeros(const std::vector<int>& shape);
    static GPUTensor ones(const std::vector<int>& shape);
    static GPUTensor full(const std::vector<int>& shape, float value);
    static GPUTensor randn(const std::vector<int>& shape);
    static GPUTensor rand(const std::vector<int>& shape);
    static GPUTensor eye(int n);
    static GPUTensor arange(float start, float end, float step = 1.0f);
    static GPUTensor linspace(float start, float end, int steps);
    
    // I/O
    void save(const std::string& path) const;
    static GPUTensor load(const std::string& path);
    void print(const std::string& name = "") const;
    
private:
    std::vector<int> shape_;
    std::vector<int> strides_;
    size_t numel_;
    
    float* h_data_ = nullptr;   // Host (CPU) data
    float* d_data_ = nullptr;   // Device (GPU) data
    bool on_gpu_ = false;
    bool owns_data_ = true;
    
    void allocate();
    void deallocate();
    void computeStrides();
};

// =============================================================================
// Complex Tensor for FFT Operations
// =============================================================================

/**
 * @brief Complex-valued tensor for FFT operations
 */
class GPUComplexTensor {
public:
    GPUComplexTensor();
    GPUComplexTensor(const std::vector<int>& shape);
    ~GPUComplexTensor();
    
    GPUComplexTensor(GPUComplexTensor&& other) noexcept;
    GPUComplexTensor& operator=(GPUComplexTensor&& other) noexcept;
    
    size_t numel() const { return numel_; }
    const std::vector<int>& shape() const { return shape_; }
    bool isOnGPU() const { return on_gpu_; }
    
    // Complex data as interleaved real/imag (float2 or cuComplex)
    void* data() { return on_gpu_ ? d_data_ : h_data_; }
    const void* data() const { return on_gpu_ ? d_data_ : h_data_; }
    
    void toGPU();
    void toCPU();
    
    // Complex operations
    GPUComplexTensor conj() const;
    GPUTensor abs() const;
    GPUTensor angle() const;
    GPUTensor real() const;
    GPUTensor imag() const;
    
    // Element-wise complex multiply
    GPUComplexTensor operator*(const GPUComplexTensor& other) const;
    GPUComplexTensor& operator*=(const GPUComplexTensor& other);
    GPUComplexTensor& operator+=(const GPUComplexTensor& other);
    
private:
    std::vector<int> shape_;
    size_t numel_;
    void* h_data_ = nullptr;
    void* d_data_ = nullptr;
    bool on_gpu_ = false;
};

// =============================================================================
// GPU FFT Operations
// =============================================================================

/**
 * @brief GPU-accelerated FFT operations for spectral neural operators
 */
class GPUFFT {
public:
    GPUFFT();
    ~GPUFFT();
    
    // 1D FFT
    void fft(const GPUTensor& input, GPUComplexTensor& output);
    void ifft(const GPUComplexTensor& input, GPUTensor& output);
    void rfft(const GPUTensor& input, GPUComplexTensor& output);
    void irfft(const GPUComplexTensor& input, GPUTensor& output, int output_size);
    
    // 2D FFT
    void fft2(const GPUTensor& input, GPUComplexTensor& output);
    void ifft2(const GPUComplexTensor& input, GPUTensor& output);
    void rfft2(const GPUTensor& input, GPUComplexTensor& output);
    void irfft2(const GPUComplexTensor& input, GPUTensor& output, const std::vector<int>& output_shape);
    
    // 3D FFT
    void fft3(const GPUTensor& input, GPUComplexTensor& output);
    void ifft3(const GPUComplexTensor& input, GPUTensor& output);
    void rfft3(const GPUTensor& input, GPUComplexTensor& output);
    void irfft3(const GPUComplexTensor& input, GPUTensor& output, const std::vector<int>& output_shape);
    
    // Batched FFT (first dim is batch)
    void batchedRfft2(const GPUTensor& input, GPUComplexTensor& output);
    void batchedIrfft2(const GPUComplexTensor& input, GPUTensor& output, const std::vector<int>& output_shape);
    void batchedRfft3(const GPUTensor& input, GPUComplexTensor& output);
    void batchedIrfft3(const GPUComplexTensor& input, GPUTensor& output, const std::vector<int>& output_shape);
    
    // Spectral convolution (for FNO)
    void spectralConv2d(const GPUComplexTensor& input, const GPUComplexTensor& weights,
                        GPUComplexTensor& output, int modes1, int modes2);
    void spectralConv3d(const GPUComplexTensor& input, const GPUComplexTensor& weights,
                        GPUComplexTensor& output, int modes1, int modes2, int modes3);
    
private:
#ifdef USE_CUDA
    cufftHandle plan_1d_;
    cufftHandle plan_2d_;
    cufftHandle plan_3d_;
    bool plans_created_ = false;
#endif
    
    void createPlan(int n, int batch, int dim);
};

// =============================================================================
// GPU BLAS Operations
// =============================================================================

/**
 * @brief GPU-accelerated BLAS operations
 */
class GPUBLAS {
public:
    static GPUBLAS& getInstance();
    
    // Level 1 BLAS
    void axpy(int n, float alpha, const float* x, float* y);  // y = alpha*x + y
    float dot(int n, const float* x, const float* y);
    float nrm2(int n, const float* x);
    void scal(int n, float alpha, float* x);
    
    // Level 2 BLAS
    void gemv(bool transA, int m, int n, float alpha, const float* A,
              const float* x, float beta, float* y);
    
    // Level 3 BLAS
    void gemm(bool transA, bool transB, int m, int n, int k,
              float alpha, const float* A, const float* B,
              float beta, float* C);
    
    // Batched GEMM
    void gemmBatched(bool transA, bool transB, int m, int n, int k,
                     float alpha, const float** A, const float** B,
                     float beta, float** C, int batch_count);
    
    void gemmStridedBatched(bool transA, bool transB, int m, int n, int k,
                            float alpha, const float* A, int strideA,
                            const float* B, int strideB,
                            float beta, float* C, int strideC, int batch_count);
    
private:
    GPUBLAS();
    ~GPUBLAS();
    
#ifdef USE_CUDA
    cublasHandle_t handle_;
#endif
};

// =============================================================================
// GPU Neural Network Layers
// =============================================================================

/**
 * @brief GPU-accelerated linear (fully connected) layer
 */
class GPULinear {
public:
    GPULinear(int in_features, int out_features, bool bias = true);
    
    GPUTensor forward(const GPUTensor& input);
    GPUTensor backward(const GPUTensor& grad_output);
    
    // Parameters
    GPUTensor& weight() { return weight_; }
    GPUTensor& bias() { return bias_; }
    GPUTensor& gradWeight() { return grad_weight_; }
    GPUTensor& gradBias() { return grad_bias_; }
    
    void zeroGrad();
    
private:
    int in_features_, out_features_;
    bool has_bias_;
    GPUTensor weight_;
    GPUTensor bias_;
    GPUTensor grad_weight_;
    GPUTensor grad_bias_;
    GPUTensor input_cache_;
};

/**
 * @brief GPU-accelerated 2D convolution
 */
class GPUConv2d {
public:
    GPUConv2d(int in_channels, int out_channels, int kernel_size,
              int stride = 1, int padding = 0, bool bias = true);
    
    GPUTensor forward(const GPUTensor& input);
    GPUTensor backward(const GPUTensor& grad_output);
    
    GPUTensor& weight() { return weight_; }
    GPUTensor& bias() { return bias_; }
    
    void zeroGrad();
    
private:
    int in_channels_, out_channels_;
    int kernel_size_, stride_, padding_;
    bool has_bias_;
    GPUTensor weight_;
    GPUTensor bias_;
    GPUTensor grad_weight_;
    GPUTensor grad_bias_;
};

/**
 * @brief GPU-accelerated layer normalization
 */
class GPULayerNorm {
public:
    GPULayerNorm(int normalized_shape, float eps = 1e-5f);
    
    GPUTensor forward(const GPUTensor& input);
    GPUTensor backward(const GPUTensor& grad_output);
    
    GPUTensor& gamma() { return gamma_; }
    GPUTensor& beta() { return beta_; }
    
    void zeroGrad();
    
private:
    int normalized_shape_;
    float eps_;
    GPUTensor gamma_;
    GPUTensor beta_;
    GPUTensor grad_gamma_;
    GPUTensor grad_beta_;
    GPUTensor mean_cache_;
    GPUTensor var_cache_;
    GPUTensor input_cache_;
};

/**
 * @brief GPU-accelerated batch normalization
 */
class GPUBatchNorm {
public:
    GPUBatchNorm(int num_features, float eps = 1e-5f, float momentum = 0.1f);
    
    GPUTensor forward(const GPUTensor& input, bool training = true);
    GPUTensor backward(const GPUTensor& grad_output);
    
    void setTraining(bool training) { training_ = training; }
    
private:
    int num_features_;
    float eps_, momentum_;
    bool training_ = true;
    GPUTensor gamma_;
    GPUTensor beta_;
    GPUTensor running_mean_;
    GPUTensor running_var_;
};

/**
 * @brief GPU-accelerated spectral convolution for FNO
 */
class GPUSpectralConv2d {
public:
    GPUSpectralConv2d(int in_channels, int out_channels, int modes1, int modes2);
    
    GPUTensor forward(const GPUTensor& input);
    GPUTensor backward(const GPUTensor& grad_output);
    
    GPUComplexTensor& weights() { return weights_; }
    
    void zeroGrad();
    
private:
    int in_channels_, out_channels_;
    int modes1_, modes2_;
    GPUComplexTensor weights_;
    GPUComplexTensor grad_weights_;
    GPUFFT fft_;
};

/**
 * @brief GPU-accelerated spectral convolution for 3D FNO
 */
class GPUSpectralConv3d {
public:
    GPUSpectralConv3d(int in_channels, int out_channels, int modes1, int modes2, int modes3);
    
    GPUTensor forward(const GPUTensor& input);
    GPUTensor backward(const GPUTensor& grad_output);
    
private:
    int in_channels_, out_channels_;
    int modes1_, modes2_, modes3_;
    GPUComplexTensor weights_;
    GPUFFT fft_;
};

// =============================================================================
// GPU Adam Optimizer
// =============================================================================

/**
 * @brief GPU-accelerated Adam optimizer
 */
class GPUAdam {
public:
    GPUAdam(std::vector<GPUTensor*> params, float lr = 1e-3f,
            float beta1 = 0.9f, float beta2 = 0.999f,
            float eps = 1e-8f, float weight_decay = 0.0f);
    
    void step();
    void zeroGrad();
    void setLR(float lr) { lr_ = lr; }
    float getLR() const { return lr_; }
    
private:
    std::vector<GPUTensor*> params_;
    std::vector<GPUTensor> m_;  // First moment
    std::vector<GPUTensor> v_;  // Second moment
    float lr_, beta1_, beta2_, eps_, weight_decay_;
    int t_ = 0;
};

// =============================================================================
// Physics GPU Kernels
// =============================================================================

/**
 * @brief GPU kernels for atmospheric physics
 */
class GPUAtmosphericKernels {
public:
    // Temperature/pressure/density profiles
    static void computeAtmosphericState(
        const float* altitudes, int n,
        float* temperature, float* pressure, float* density,
        int model_type);  // 0=US76, 1=ICAO, 2=NRLMSISE
    
    // Wind profile interpolation
    static void interpolateWind(
        const float* altitudes, int n,
        const float* wind_profile_alt, const float* wind_u, const float* wind_v,
        int profile_size, float* wind_out_u, float* wind_out_v);
    
    // Sound speed computation
    static void computeSoundSpeed(
        const float* temperature, const float* gamma, int n,
        float* sound_speed);
    
    // Atmospheric refraction (for infrasound)
    static void computeRefraction(
        const float* sound_speed, const float* wind_u, const float* wind_v,
        const float* ray_dx, const float* ray_dy, const float* ray_dz,
        int n, float* effective_c);
};

/**
 * @brief GPU kernels for radiation physics
 */
class GPURadiationKernels {
public:
    // Blackbody radiation
    static void computePlanckFunction(
        const float* temperature, const float* wavelength,
        int n_temp, int n_wave, float* radiance);
    
    // Absorption/emission
    static void computeAbsorption(
        const float* density, const float* cross_section,
        const float* path_length, int n, float* absorption);
    
    // Radiative transfer (discrete ordinates)
    static void discreteOrdinates(
        const float* source, const float* absorption,
        const float* scattering, const float* weights,
        int nx, int ny, int n_angles, float* intensity);
    
    // Thermal radiation exchange
    static void computeViewFactors(
        const float* positions, const float* normals, const float* areas,
        int n_surfaces, float* view_factors);
    
    // Nuclear fireball radiation
    static void fireballRadiation(
        float yield_kt, float time, float altitude,
        const float* observer_positions, int n_observers,
        float* thermal_flux, float* peak_temperature);
    
    // Fallout particle transport
    static void advectParticles(
        float* positions, float* velocities, const float* settling_vel,
        const float* wind_u, const float* wind_v, const float* wind_w,
        int n_particles, float dt);
    
    // Radioactive decay
    static void computeDecay(
        const float* activity, const float* half_life,
        float dt, int n_species, float* new_activity);
    
    // Dose computation
    static void computeDoseRate(
        const float* activity, const float* dose_conversion,
        const float* distances, int n_sources, int n_points,
        float* dose_rate);
};

/**
 * @brief GPU kernels for explosion/blast physics
 */
class GPUExplosionKernels {
public:
    // Blast wave (Sedov-Taylor)
    static void sedovTaylorBlast(
        float energy, float rho0, float gamma, float time,
        const float* radii, int n, float* pressure, float* density, float* velocity);
    
    // Nuclear fireball evolution
    static void fireballDynamics(
        float yield_kt, float altitude, float time,
        float* radius, float* temperature, float* rise_velocity);
    
    // Rayleigh-Taylor instability growth
    static void rayleighTaylorGrowth(
        float* interface_perturbation, float atwood_number,
        float g, float dt, int nx, int ny);
    
    // Kelvin-Helmholtz instability
    static void kelvinHelmholtzGrowth(
        float* vorticity, const float* velocity_u, const float* velocity_v,
        float density_ratio, float dt, int nx, int ny);
    
    // Mushroom cloud dynamics
    static void mushroomCloudStep(
        float* density, float* temperature, float* velocity_u,
        float* velocity_v, float* velocity_w,
        const float* ambient_density, const float* ambient_temp,
        float g, float dt, int nx, int ny, int nz);
    
    // Ground shock coupling
    static void groundShockCoupling(
        const float* air_pressure, const float* ground_impedance,
        int n_interface, float* ground_velocity, float* ground_stress);
};

/**
 * @brief GPU kernels for fluid instabilities
 */
class GPUInstabilityKernels {
public:
    // General interface tracking (VOF style)
    static void advectVOF(
        float* volume_fraction, const float* velocity_u, const float* velocity_v,
        float dt, int nx, int ny);
    
    // Density-driven instability (buoyancy)
    static void computeBuoyancy(
        const float* density, float reference_density, float g,
        int nx, int ny, int nz, float* buoyancy_force);
    
    // Viscous fingering
    static void viscousFingeringStep(
        float* saturation, const float* pressure,
        float mobility_ratio, float dt, int nx, int ny);
    
    // Richtmyer-Meshkov instability
    static void richtmyerMeshkovGrowth(
        float* interface_amplitude, float atwood_number,
        float shock_mach, float time, int n_modes);
};

// =============================================================================
// GPU Compute Context
// =============================================================================

/**
 * @brief Global GPU compute context and utilities
 */
class GPUComputeContext {
public:
    static GPUComputeContext& getInstance();
    
    // Initialization
    bool initialize();
    bool isInitialized() const { return initialized_; }
    bool isGPUAvailable() const;
    
    // Device info
    int getDeviceCount() const;
    std::string getDeviceName(int device = 0) const;
    size_t getTotalMemory(int device = 0) const;
    size_t getFreeMemory(int device = 0) const;
    
    // Stream management
    void* createStream();
    void destroyStream(void* stream);
    void synchronizeStream(void* stream);
    void synchronizeDevice();
    
    // Memory pool (for efficient allocation)
    void* allocate(size_t bytes);
    void deallocate(void* ptr);
    void enableMemoryPool(bool enable);
    
    // Random number generation
    void setSeed(unsigned long seed);
    void generateUniform(float* data, size_t n, float low = 0.0f, float high = 1.0f);
    void generateNormal(float* data, size_t n, float mean = 0.0f, float std = 1.0f);
    
    // Profiling
    void startTimer();
    float stopTimer();  // Returns ms
    void enableProfiling(bool enable);
    
    // Error handling
    std::string getLastError() const;
    void checkError(const char* msg);
    
private:
    GPUComputeContext();
    ~GPUComputeContext();
    
    bool initialized_ = false;
    
#ifdef USE_CUDA
    curandGenerator_t curand_gen_;
    cudaEvent_t start_event_, stop_event_;
#endif
};

// =============================================================================
// Utility Functions
// =============================================================================

namespace GPUUtils {
    // Element-wise operations (CUDA kernels called internally)
    void add(const float* a, const float* b, float* c, size_t n);
    void sub(const float* a, const float* b, float* c, size_t n);
    void mul(const float* a, const float* b, float* c, size_t n);
    void div(const float* a, const float* b, float* c, size_t n);
    void scale(const float* a, float s, float* b, size_t n);
    void addScalar(const float* a, float s, float* b, size_t n);
    
    // Activations
    void relu(const float* in, float* out, size_t n);
    void sigmoid(const float* in, float* out, size_t n);
    void tanh(const float* in, float* out, size_t n);
    void gelu(const float* in, float* out, size_t n);
    void silu(const float* in, float* out, size_t n);
    void softmax(const float* in, float* out, int batch, int dim);
    
    // Reductions
    float sum(const float* data, size_t n);
    float mean(const float* data, size_t n);
    float max(const float* data, size_t n);
    float min(const float* data, size_t n);
    float norm(const float* data, size_t n, int p = 2);
    
    // Copy operations
    void copy(const float* src, float* dst, size_t n);
    void fill(float* data, float value, size_t n);
    void zeros(float* data, size_t n);
    void ones(float* data, size_t n);
    
    // Complex operations
    void complexMul(const void* a, const void* b, void* c, size_t n);
    void complexConj(const void* a, void* b, size_t n);
    void complexAbs(const void* a, float* b, size_t n);
}

// =============================================================================
// High-level API for Physics Simulations
// =============================================================================

/**
 * @brief Run physics simulation on GPU
 */
class GPUPhysicsSimulator {
public:
    GPUPhysicsSimulator();
    ~GPUPhysicsSimulator();
    
    // Set domain
    void setDomain(int nx, int ny, int nz, float dx, float dy, float dz);
    
    // Allocate fields
    void allocateField(const std::string& name, int components = 1);
    GPUTensor& getField(const std::string& name);
    
    // Physics modules
    void enableAtmospheric(bool enable);
    void enableRadiation(bool enable);
    void enableExplosion(bool enable);
    
    // Time stepping
    void step(float dt);
    void run(float t_final, float dt);
    
    // Output
    void writeOutput(const std::string& filename, int step);
    
private:
    int nx_, ny_, nz_;
    float dx_, dy_, dz_;
    std::map<std::string, GPUTensor> fields_;
    
    bool use_atmospheric_ = false;
    bool use_radiation_ = false;
    bool use_explosion_ = false;
};

} // namespace GPU
} // namespace FSRM

#endif // GPU_COMPUTE_HPP
