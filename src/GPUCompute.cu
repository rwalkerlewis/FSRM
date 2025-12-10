/**
 * @file GPUCompute.cu
 * @brief CUDA implementation of GPU compute operations for physics and ML
 */

#include "GPUCompute.hpp"
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <numeric>

#ifdef USE_CUDA
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cufft.h>
#include <curand.h>
#include <curand_kernel.h>

// Error checking macro
#define CUDA_CHECK(call) \
    do { \
        cudaError_t err = call; \
        if (err != cudaSuccess) { \
            throw std::runtime_error(std::string("CUDA error: ") + cudaGetErrorString(err)); \
        } \
    } while(0)

#define CUBLAS_CHECK(call) \
    do { \
        cublasStatus_t status = call; \
        if (status != CUBLAS_STATUS_SUCCESS) { \
            throw std::runtime_error("cuBLAS error"); \
        } \
    } while(0)

#define CUFFT_CHECK(call) \
    do { \
        cufftResult status = call; \
        if (status != CUFFT_SUCCESS) { \
            throw std::runtime_error("cuFFT error"); \
        } \
    } while(0)

#endif // USE_CUDA

namespace FSRM {
namespace GPU {

// =============================================================================
// CUDA Kernels
// =============================================================================

#ifdef USE_CUDA

// Block size constants
constexpr int BLOCK_SIZE = 256;
constexpr int BLOCK_SIZE_2D = 16;

// -----------------------------------------------------------------------------
// Element-wise Kernels
// -----------------------------------------------------------------------------

__global__ void kernelAdd(const float* a, const float* b, float* c, size_t n) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) c[idx] = a[idx] + b[idx];
}

__global__ void kernelSub(const float* a, const float* b, float* c, size_t n) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) c[idx] = a[idx] - b[idx];
}

__global__ void kernelMul(const float* a, const float* b, float* c, size_t n) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) c[idx] = a[idx] * b[idx];
}

__global__ void kernelDiv(const float* a, const float* b, float* c, size_t n) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) c[idx] = a[idx] / b[idx];
}

__global__ void kernelScale(const float* a, float s, float* b, size_t n) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) b[idx] = a[idx] * s;
}

__global__ void kernelAddScalar(const float* a, float s, float* b, size_t n) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) b[idx] = a[idx] + s;
}

__global__ void kernelFill(float* data, float value, size_t n) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) data[idx] = value;
}

// -----------------------------------------------------------------------------
// Activation Kernels
// -----------------------------------------------------------------------------

__global__ void kernelRelu(const float* in, float* out, size_t n) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) out[idx] = fmaxf(0.0f, in[idx]);
}

__global__ void kernelSigmoid(const float* in, float* out, size_t n) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) out[idx] = 1.0f / (1.0f + expf(-in[idx]));
}

__global__ void kernelTanh(const float* in, float* out, size_t n) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) out[idx] = tanhf(in[idx]);
}

__global__ void kernelGelu(const float* in, float* out, size_t n) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        float x = in[idx];
        // GELU approximation: 0.5 * x * (1 + tanh(sqrt(2/pi) * (x + 0.044715 * x^3)))
        float cdf = 0.5f * (1.0f + tanhf(0.7978845608f * (x + 0.044715f * x * x * x)));
        out[idx] = x * cdf;
    }
}

__global__ void kernelSilu(const float* in, float* out, size_t n) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        float x = in[idx];
        out[idx] = x / (1.0f + expf(-x));
    }
}

__global__ void kernelSoftmax(const float* in, float* out, int batch, int dim) {
    int b = blockIdx.x;
    if (b >= batch) return;
    
    const float* in_row = in + b * dim;
    float* out_row = out + b * dim;
    
    // Find max for numerical stability
    float max_val = in_row[0];
    for (int i = 1; i < dim; i++) {
        max_val = fmaxf(max_val, in_row[i]);
    }
    
    // Compute exp and sum
    float sum = 0.0f;
    for (int i = 0; i < dim; i++) {
        out_row[i] = expf(in_row[i] - max_val);
        sum += out_row[i];
    }
    
    // Normalize
    for (int i = 0; i < dim; i++) {
        out_row[i] /= sum;
    }
}

// -----------------------------------------------------------------------------
// Reduction Kernels
// -----------------------------------------------------------------------------

__global__ void kernelReduceSum(const float* data, float* result, size_t n) {
    __shared__ float sdata[BLOCK_SIZE];
    
    size_t tid = threadIdx.x;
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    sdata[tid] = (idx < n) ? data[idx] : 0.0f;
    __syncthreads();
    
    // Reduction in shared memory
    for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            sdata[tid] += sdata[tid + s];
        }
        __syncthreads();
    }
    
    if (tid == 0) {
        atomicAdd(result, sdata[0]);
    }
}

__global__ void kernelReduceMax(const float* data, float* result, size_t n) {
    __shared__ float sdata[BLOCK_SIZE];
    
    size_t tid = threadIdx.x;
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    sdata[tid] = (idx < n) ? data[idx] : -1e30f;
    __syncthreads();
    
    for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            sdata[tid] = fmaxf(sdata[tid], sdata[tid + s]);
        }
        __syncthreads();
    }
    
    if (tid == 0) {
        // Atomic max for floats (using int atomicMax with reinterpret)
        int* address_as_int = (int*)result;
        int old = *address_as_int, assumed;
        do {
            assumed = old;
            old = atomicCAS(address_as_int, assumed, 
                           __float_as_int(fmaxf(__int_as_float(assumed), sdata[0])));
        } while (assumed != old);
    }
}

__global__ void kernelReduceMin(const float* data, float* result, size_t n) {
    __shared__ float sdata[BLOCK_SIZE];
    
    size_t tid = threadIdx.x;
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    sdata[tid] = (idx < n) ? data[idx] : 1e30f;
    __syncthreads();
    
    for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            sdata[tid] = fminf(sdata[tid], sdata[tid + s]);
        }
        __syncthreads();
    }
    
    if (tid == 0) {
        int* address_as_int = (int*)result;
        int old = *address_as_int, assumed;
        do {
            assumed = old;
            old = atomicCAS(address_as_int, assumed,
                           __float_as_int(fminf(__int_as_float(assumed), sdata[0])));
        } while (assumed != old);
    }
}

// -----------------------------------------------------------------------------
// Complex Kernels (for FFT)
// -----------------------------------------------------------------------------

__global__ void kernelComplexMul(const float2* a, const float2* b, float2* c, size_t n) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        float2 av = a[idx];
        float2 bv = b[idx];
        c[idx].x = av.x * bv.x - av.y * bv.y;
        c[idx].y = av.x * bv.y + av.y * bv.x;
    }
}

__global__ void kernelComplexConj(const float2* a, float2* b, size_t n) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        b[idx].x = a[idx].x;
        b[idx].y = -a[idx].y;
    }
}

__global__ void kernelComplexAbs(const float2* a, float* b, size_t n) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        float2 av = a[idx];
        b[idx] = sqrtf(av.x * av.x + av.y * av.y);
    }
}

// -----------------------------------------------------------------------------
// Layer Normalization Kernel
// -----------------------------------------------------------------------------

__global__ void kernelLayerNorm(const float* input, float* output,
                                const float* gamma, const float* beta,
                                int batch, int dim, float eps) {
    int b = blockIdx.x;
    if (b >= batch) return;
    
    const float* in_row = input + b * dim;
    float* out_row = output + b * dim;
    
    // Compute mean
    float mean = 0.0f;
    for (int i = 0; i < dim; i++) {
        mean += in_row[i];
    }
    mean /= dim;
    
    // Compute variance
    float var = 0.0f;
    for (int i = 0; i < dim; i++) {
        float diff = in_row[i] - mean;
        var += diff * diff;
    }
    var /= dim;
    
    // Normalize
    float inv_std = rsqrtf(var + eps);
    for (int i = 0; i < dim; i++) {
        out_row[i] = gamma[i] * (in_row[i] - mean) * inv_std + beta[i];
    }
}

// -----------------------------------------------------------------------------
// Atmospheric Physics Kernels
// -----------------------------------------------------------------------------

__global__ void kernelUSStandardAtmosphere(const float* altitudes, int n,
                                           float* temperature, float* pressure, float* density) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    
    float h = altitudes[idx];  // meters
    float T, P;
    
    // US Standard Atmosphere 1976
    const float T0 = 288.15f;
    const float P0 = 101325.0f;
    const float g = 9.80665f;
    const float R = 287.05f;
    
    if (h < 11000.0f) {
        // Troposphere
        float L = -0.0065f;
        T = T0 + L * h;
        P = P0 * powf(T / T0, -g / (L * R));
    } else if (h < 20000.0f) {
        // Tropopause (isothermal)
        T = 216.65f;
        float P_11 = 22632.0f;
        P = P_11 * expf(-g * (h - 11000.0f) / (R * T));
    } else if (h < 32000.0f) {
        // Stratosphere
        float L = 0.001f;
        T = 216.65f + L * (h - 20000.0f);
        float P_20 = 5474.9f;
        P = P_20 * powf(T / 216.65f, -g / (L * R));
    } else {
        // Simplified upper atmosphere
        T = 228.65f + 0.0028f * (h - 32000.0f);
        T = fminf(T, 270.65f);
        float scale_h = R * T / g;
        P = 868.0f * expf(-(h - 32000.0f) / scale_h);
    }
    
    temperature[idx] = T;
    pressure[idx] = P;
    density[idx] = P / (R * T);
}

__global__ void kernelComputeSoundSpeed(const float* temperature, const float* gamma_arr,
                                        int n, float* sound_speed) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    
    const float R = 287.05f;
    float gamma = (gamma_arr != nullptr) ? gamma_arr[idx] : 1.4f;
    sound_speed[idx] = sqrtf(gamma * R * temperature[idx]);
}

// -----------------------------------------------------------------------------
// Radiation Physics Kernels
// -----------------------------------------------------------------------------

__global__ void kernelPlanckFunction(const float* temperature, const float* wavelength,
                                     int n_temp, int n_wave, float* radiance) {
    int i_temp = blockIdx.x * blockDim.x + threadIdx.x;
    int i_wave = blockIdx.y * blockDim.y + threadIdx.y;
    
    if (i_temp >= n_temp || i_wave >= n_wave) return;
    
    const float h = 6.62607e-34f;  // Planck constant
    const float c = 2.998e8f;       // Speed of light
    const float k = 1.38065e-23f;   // Boltzmann constant
    
    float T = temperature[i_temp];
    float lambda = wavelength[i_wave];
    
    float c1 = 2.0f * h * c * c;
    float c2 = h * c / k;
    
    float exponent = c2 / (lambda * T);
    float B = c1 / (lambda * lambda * lambda * lambda * lambda * (expf(exponent) - 1.0f));
    
    radiance[i_temp * n_wave + i_wave] = B;
}

__global__ void kernelRadioactiveDecay(const float* activity, const float* half_life,
                                       float dt, int n_species, float* new_activity) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_species) return;
    
    float decay_const = 0.693147f / half_life[idx];  // ln(2) / T_half
    new_activity[idx] = activity[idx] * expf(-decay_const * dt);
}

__global__ void kernelComputeDoseRate(const float* activity, const float* dose_conversion,
                                      const float* source_pos, const float* observer_pos,
                                      int n_sources, int n_observers, float* dose_rate) {
    int obs = blockIdx.x * blockDim.x + threadIdx.x;
    if (obs >= n_observers) return;
    
    float ox = observer_pos[obs * 3];
    float oy = observer_pos[obs * 3 + 1];
    float oz = observer_pos[obs * 3 + 2];
    
    float total_dose = 0.0f;
    for (int src = 0; src < n_sources; src++) {
        float sx = source_pos[src * 3];
        float sy = source_pos[src * 3 + 1];
        float sz = source_pos[src * 3 + 2];
        
        float dx = ox - sx;
        float dy = oy - sy;
        float dz = oz - sz;
        float r2 = dx*dx + dy*dy + dz*dz;
        
        if (r2 > 1e-6f) {
            total_dose += activity[src] * dose_conversion[src] / (4.0f * 3.14159f * r2);
        }
    }
    
    dose_rate[obs] = total_dose;
}

// -----------------------------------------------------------------------------
// Explosion Physics Kernels
// -----------------------------------------------------------------------------

__global__ void kernelSedovTaylor(float energy, float rho0, float gamma, float time,
                                  const float* radii, int n,
                                  float* pressure, float* density, float* velocity) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    
    float r = radii[idx];
    
    // Sedov-Taylor similarity solution
    float alpha = 1.15167f;  // For gamma = 1.4
    float R_shock = alpha * powf(energy * time * time / rho0, 0.2f);
    
    if (r > R_shock) {
        // Ambient conditions
        density[idx] = rho0;
        pressure[idx] = 1e5f;  // Ambient pressure
        velocity[idx] = 0.0f;
    } else {
        float xi = r / R_shock;
        
        // Self-similar profiles (simplified)
        float rho_ratio = ((gamma + 1.0f) / (gamma - 1.0f)) * 
                         powf(1.0f - 0.5f * (gamma - 1.0f) / (gamma + 1.0f) * xi * xi, 
                              1.0f / (gamma - 1.0f));
        density[idx] = rho0 * rho_ratio;
        
        float v_shock = 0.4f * R_shock / time;
        velocity[idx] = v_shock * xi;
        
        float P_shock = 2.0f * rho0 * v_shock * v_shock / (gamma + 1.0f);
        pressure[idx] = P_shock * (1.0f - 0.5f * xi * xi);
    }
}

__global__ void kernelFireballDynamics(float yield_kt, float altitude, float time,
                                       float* radius, float* temperature, float* rise_velocity) {
    // Nuclear fireball scaling laws
    float W = yield_kt;
    
    // Fireball radius (Glasstone & Dolan)
    float R_max = 145.0f * powf(W, 0.4f);  // meters
    float t_max = 0.038f * powf(W, 0.45f); // seconds for max radius
    
    if (time < t_max) {
        *radius = R_max * powf(time / t_max, 0.4f);
        *temperature = 8000.0f * powf(t_max / fmaxf(time, 0.001f), 0.3f);
    } else {
        *radius = R_max;
        *temperature = 3000.0f * powf(t_max / time, 0.5f);
    }
    
    // Buoyant rise velocity
    float T_amb = 288.15f - 0.0065f * altitude;
    float delta_T = *temperature - T_amb;
    float g = 9.81f;
    
    if (delta_T > 0 && time > t_max) {
        *rise_velocity = sqrtf(2.0f * g * (*radius) * delta_T / T_amb);
        *rise_velocity = fminf(*rise_velocity, 100.0f);  // Cap at 100 m/s
    } else {
        *rise_velocity = 0.0f;
    }
}

__global__ void kernelRayleighTaylor(float* perturbation, float atwood, float g, float dt,
                                     int nx, int ny) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    
    if (i >= nx || j >= ny) return;
    int idx = i * ny + j;
    
    // Simplified RT growth: η(t) = η_0 * exp(sqrt(A*g*k) * t)
    // Here we use growth rate directly
    float eta = perturbation[idx];
    float k = 2.0f * 3.14159f / (float)ny;  // Characteristic wavenumber
    float growth_rate = sqrtf(fabsf(atwood) * g * k);
    
    perturbation[idx] = eta * expf(growth_rate * dt);
}

__global__ void kernelMushroomCloud(float* rho, float* T, float* u, float* v, float* w,
                                    const float* rho_amb, const float* T_amb,
                                    float g, float dt, int nx, int ny, int nz) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;
    
    if (i >= nx || j >= ny || k >= nz) return;
    int idx = i + j * nx + k * nx * ny;
    
    // Buoyancy force
    float drho = rho[idx] - rho_amb[k];
    float buoyancy = -g * drho / rho[idx];
    
    // Update vertical velocity
    w[idx] += buoyancy * dt;
    
    // Simple advection (upwind)
    if (w[idx] > 0 && k > 0) {
        int idx_below = i + j * nx + (k-1) * nx * ny;
        rho[idx] -= dt * w[idx] * (rho[idx] - rho[idx_below]) / 100.0f;  // dz = 100m
        T[idx] -= dt * w[idx] * (T[idx] - T[idx_below]) / 100.0f;
    }
    
    // Entrainment (simplified)
    float entrainment_rate = 0.1f;
    rho[idx] += entrainment_rate * dt * (rho_amb[k] - rho[idx]);
    T[idx] += entrainment_rate * dt * (T_amb[k] - T[idx]);
}

// -----------------------------------------------------------------------------
// Particle Transport Kernels
// -----------------------------------------------------------------------------

__global__ void kernelAdvectParticles(float* pos_x, float* pos_y, float* pos_z,
                                      float* vel_x, float* vel_y, float* vel_z,
                                      const float* settling_vel,
                                      const float* wind_u, const float* wind_v, const float* wind_w,
                                      int n_particles, float dt) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_particles) return;
    
    // Update velocity with wind
    vel_x[idx] = wind_u[idx];
    vel_y[idx] = wind_v[idx];
    vel_z[idx] = wind_w[idx] - settling_vel[idx];
    
    // Update position
    pos_x[idx] += vel_x[idx] * dt;
    pos_y[idx] += vel_y[idx] * dt;
    pos_z[idx] += vel_z[idx] * dt;
    
    // Ground check
    if (pos_z[idx] < 0.0f) {
        pos_z[idx] = 0.0f;
        vel_z[idx] = 0.0f;
    }
}

#endif // USE_CUDA

// =============================================================================
// GPUTensor Implementation
// =============================================================================

GPUTensor::GPUTensor() : numel_(0), on_gpu_(false), owns_data_(true) {}

GPUTensor::GPUTensor(const std::vector<int>& shape) : shape_(shape), on_gpu_(false), owns_data_(true) {
    computeStrides();
    numel_ = 1;
    for (int s : shape_) numel_ *= s;
    allocate();
}

GPUTensor::GPUTensor(const std::vector<int>& shape, float value) : GPUTensor(shape) {
    fill(value);
}

GPUTensor::GPUTensor(const std::vector<int>& shape, const float* data, bool on_gpu)
    : shape_(shape), on_gpu_(on_gpu), owns_data_(true) {
    computeStrides();
    numel_ = 1;
    for (int s : shape_) numel_ *= s;
    allocate();
    
#ifdef USE_CUDA
    if (on_gpu && on_gpu_) {
        CUDA_CHECK(cudaMemcpy(d_data_, data, numel_ * sizeof(float), cudaMemcpyDeviceToDevice));
    } else if (!on_gpu && on_gpu_) {
        CUDA_CHECK(cudaMemcpy(d_data_, data, numel_ * sizeof(float), cudaMemcpyHostToDevice));
    } else if (on_gpu && !on_gpu_) {
        CUDA_CHECK(cudaMemcpy(h_data_, data, numel_ * sizeof(float), cudaMemcpyDeviceToHost));
    } else {
        std::copy(data, data + numel_, h_data_);
    }
#else
    std::copy(data, data + numel_, h_data_);
#endif
}

GPUTensor::~GPUTensor() {
    if (owns_data_) deallocate();
}

GPUTensor::GPUTensor(GPUTensor&& other) noexcept
    : shape_(std::move(other.shape_)), strides_(std::move(other.strides_)),
      numel_(other.numel_), h_data_(other.h_data_), d_data_(other.d_data_),
      on_gpu_(other.on_gpu_), owns_data_(other.owns_data_) {
    other.h_data_ = nullptr;
    other.d_data_ = nullptr;
    other.owns_data_ = false;
}

GPUTensor& GPUTensor::operator=(GPUTensor&& other) noexcept {
    if (this != &other) {
        if (owns_data_) deallocate();
        shape_ = std::move(other.shape_);
        strides_ = std::move(other.strides_);
        numel_ = other.numel_;
        h_data_ = other.h_data_;
        d_data_ = other.d_data_;
        on_gpu_ = other.on_gpu_;
        owns_data_ = other.owns_data_;
        other.h_data_ = nullptr;
        other.d_data_ = nullptr;
        other.owns_data_ = false;
    }
    return *this;
}

GPUTensor::GPUTensor(const GPUTensor& other)
    : shape_(other.shape_), strides_(other.strides_), numel_(other.numel_),
      on_gpu_(other.on_gpu_), owns_data_(true) {
    allocate();
#ifdef USE_CUDA
    if (on_gpu_) {
        CUDA_CHECK(cudaMemcpy(d_data_, other.d_data_, numel_ * sizeof(float), cudaMemcpyDeviceToDevice));
    } else {
        std::copy(other.h_data_, other.h_data_ + numel_, h_data_);
    }
#else
    std::copy(other.h_data_, other.h_data_ + numel_, h_data_);
#endif
}

GPUTensor& GPUTensor::operator=(const GPUTensor& other) {
    if (this != &other) {
        if (owns_data_) deallocate();
        shape_ = other.shape_;
        strides_ = other.strides_;
        numel_ = other.numel_;
        on_gpu_ = other.on_gpu_;
        owns_data_ = true;
        allocate();
#ifdef USE_CUDA
        if (on_gpu_) {
            CUDA_CHECK(cudaMemcpy(d_data_, other.d_data_, numel_ * sizeof(float), cudaMemcpyDeviceToDevice));
        } else {
            std::copy(other.h_data_, other.h_data_ + numel_, h_data_);
        }
#else
        std::copy(other.h_data_, other.h_data_ + numel_, h_data_);
#endif
    }
    return *this;
}

void GPUTensor::allocate() {
#ifdef USE_CUDA
    GPUManager& gpu = GPUManager::getInstance();
    if (gpu.isAvailable()) {
        d_data_ = static_cast<float*>(gpu.allocateDevice(numel_ * sizeof(float)));
        on_gpu_ = true;
    } else {
        h_data_ = new float[numel_];
        on_gpu_ = false;
    }
#else
    h_data_ = new float[numel_];
    on_gpu_ = false;
#endif
}

void GPUTensor::deallocate() {
#ifdef USE_CUDA
    if (d_data_) {
        GPUManager::getInstance().freeDevice(d_data_);
        d_data_ = nullptr;
    }
#endif
    if (h_data_) {
        delete[] h_data_;
        h_data_ = nullptr;
    }
}

void GPUTensor::computeStrides() {
    strides_.resize(shape_.size());
    if (!shape_.empty()) {
        strides_.back() = 1;
        for (int i = static_cast<int>(shape_.size()) - 2; i >= 0; i--) {
            strides_[i] = strides_[i + 1] * shape_[i + 1];
        }
    }
}

size_t GPUTensor::numel() const { return numel_; }

void GPUTensor::toGPU() {
#ifdef USE_CUDA
    if (on_gpu_) return;
    GPUManager& gpu = GPUManager::getInstance();
    if (!gpu.isAvailable()) return;
    
    d_data_ = static_cast<float*>(gpu.allocateDevice(numel_ * sizeof(float)));
    CUDA_CHECK(cudaMemcpy(d_data_, h_data_, numel_ * sizeof(float), cudaMemcpyHostToDevice));
    on_gpu_ = true;
#endif
}

void GPUTensor::toCPU() {
#ifdef USE_CUDA
    if (!on_gpu_) return;
    if (!h_data_) h_data_ = new float[numel_];
    CUDA_CHECK(cudaMemcpy(h_data_, d_data_, numel_ * sizeof(float), cudaMemcpyDeviceToHost));
    on_gpu_ = false;
#endif
}

void GPUTensor::synchronize() {
#ifdef USE_CUDA
    CUDA_CHECK(cudaDeviceSynchronize());
#endif
}

void GPUTensor::zeros() { fill(0.0f); }
void GPUTensor::ones() { fill(1.0f); }

void GPUTensor::fill(float value) {
#ifdef USE_CUDA
    if (on_gpu_) {
        int blocks = (numel_ + BLOCK_SIZE - 1) / BLOCK_SIZE;
        kernelFill<<<blocks, BLOCK_SIZE>>>(d_data_, value, numel_);
        CUDA_CHECK(cudaGetLastError());
    } else {
        std::fill(h_data_, h_data_ + numel_, value);
    }
#else
    std::fill(h_data_, h_data_ + numel_, value);
#endif
}

GPUTensor GPUTensor::operator+(const GPUTensor& other) const {
    GPUTensor result(shape_);
#ifdef USE_CUDA
    if (on_gpu_) {
        int blocks = (numel_ + BLOCK_SIZE - 1) / BLOCK_SIZE;
        kernelAdd<<<blocks, BLOCK_SIZE>>>(d_data_, other.d_data_, result.d_data_, numel_);
        CUDA_CHECK(cudaGetLastError());
    } else {
        for (size_t i = 0; i < numel_; i++) {
            result.h_data_[i] = h_data_[i] + other.h_data_[i];
        }
    }
#else
    for (size_t i = 0; i < numel_; i++) {
        result.h_data_[i] = h_data_[i] + other.h_data_[i];
    }
#endif
    return result;
}

GPUTensor GPUTensor::operator*(float scalar) const {
    GPUTensor result(shape_);
#ifdef USE_CUDA
    if (on_gpu_) {
        int blocks = (numel_ + BLOCK_SIZE - 1) / BLOCK_SIZE;
        kernelScale<<<blocks, BLOCK_SIZE>>>(d_data_, scalar, result.d_data_, numel_);
        CUDA_CHECK(cudaGetLastError());
    } else {
        for (size_t i = 0; i < numel_; i++) {
            result.h_data_[i] = h_data_[i] * scalar;
        }
    }
#else
    for (size_t i = 0; i < numel_; i++) {
        result.h_data_[i] = h_data_[i] * scalar;
    }
#endif
    return result;
}

GPUTensor GPUTensor::relu() const {
    GPUTensor result(shape_);
#ifdef USE_CUDA
    if (on_gpu_) {
        int blocks = (numel_ + BLOCK_SIZE - 1) / BLOCK_SIZE;
        kernelRelu<<<blocks, BLOCK_SIZE>>>(d_data_, result.d_data_, numel_);
        CUDA_CHECK(cudaGetLastError());
    } else {
        for (size_t i = 0; i < numel_; i++) {
            result.h_data_[i] = std::max(0.0f, h_data_[i]);
        }
    }
#else
    for (size_t i = 0; i < numel_; i++) {
        result.h_data_[i] = std::max(0.0f, h_data_[i]);
    }
#endif
    return result;
}

GPUTensor GPUTensor::gelu() const {
    GPUTensor result(shape_);
#ifdef USE_CUDA
    if (on_gpu_) {
        int blocks = (numel_ + BLOCK_SIZE - 1) / BLOCK_SIZE;
        kernelGelu<<<blocks, BLOCK_SIZE>>>(d_data_, result.d_data_, numel_);
        CUDA_CHECK(cudaGetLastError());
    } else {
        for (size_t i = 0; i < numel_; i++) {
            float x = h_data_[i];
            float cdf = 0.5f * (1.0f + std::tanh(0.7978845608f * (x + 0.044715f * x * x * x)));
            result.h_data_[i] = x * cdf;
        }
    }
#else
    for (size_t i = 0; i < numel_; i++) {
        float x = h_data_[i];
        float cdf = 0.5f * (1.0f + std::tanh(0.7978845608f * (x + 0.044715f * x * x * x)));
        result.h_data_[i] = x * cdf;
    }
#endif
    return result;
}

float GPUTensor::sum() const {
#ifdef USE_CUDA
    if (on_gpu_) {
        float* d_result;
        CUDA_CHECK(cudaMalloc(&d_result, sizeof(float)));
        CUDA_CHECK(cudaMemset(d_result, 0, sizeof(float)));
        
        int blocks = (numel_ + BLOCK_SIZE - 1) / BLOCK_SIZE;
        kernelReduceSum<<<blocks, BLOCK_SIZE>>>(d_data_, d_result, numel_);
        CUDA_CHECK(cudaGetLastError());
        
        float result;
        CUDA_CHECK(cudaMemcpy(&result, d_result, sizeof(float), cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaFree(d_result));
        return result;
    }
#endif
    return std::accumulate(h_data_, h_data_ + numel_, 0.0f);
}

float GPUTensor::mean() const {
    return sum() / static_cast<float>(numel_);
}

GPUTensor GPUTensor::zeros(const std::vector<int>& shape) {
    return GPUTensor(shape, 0.0f);
}

GPUTensor GPUTensor::ones(const std::vector<int>& shape) {
    return GPUTensor(shape, 1.0f);
}

// =============================================================================
// GPUBLAS Implementation
// =============================================================================

GPUBLAS& GPUBLAS::getInstance() {
    static GPUBLAS instance;
    return instance;
}

GPUBLAS::GPUBLAS() {
#ifdef USE_CUDA
    CUBLAS_CHECK(cublasCreate(&handle_));
#endif
}

GPUBLAS::~GPUBLAS() {
#ifdef USE_CUDA
    cublasDestroy(handle_);
#endif
}

void GPUBLAS::gemm(bool transA, bool transB, int m, int n, int k,
                   float alpha, const float* A, const float* B,
                   float beta, float* C) {
#ifdef USE_CUDA
    cublasOperation_t opA = transA ? CUBLAS_OP_T : CUBLAS_OP_N;
    cublasOperation_t opB = transB ? CUBLAS_OP_T : CUBLAS_OP_N;
    
    int lda = transA ? k : m;
    int ldb = transB ? n : k;
    
    CUBLAS_CHECK(cublasSgemm(handle_, opA, opB, m, n, k,
                             &alpha, A, lda, B, ldb, &beta, C, m));
#endif
}

// =============================================================================
// GPUAtmosphericKernels Implementation
// =============================================================================

void GPUAtmosphericKernels::computeAtmosphericState(
    const float* altitudes, int n,
    float* temperature, float* pressure, float* density,
    int model_type) {
#ifdef USE_CUDA
    int blocks = (n + BLOCK_SIZE - 1) / BLOCK_SIZE;
    
    switch (model_type) {
        case 0:  // US Standard 1976
        default:
            kernelUSStandardAtmosphere<<<blocks, BLOCK_SIZE>>>(
                altitudes, n, temperature, pressure, density);
            break;
    }
    CUDA_CHECK(cudaGetLastError());
#endif
}

void GPUAtmosphericKernels::computeSoundSpeed(
    const float* temperature, const float* gamma, int n,
    float* sound_speed) {
#ifdef USE_CUDA
    int blocks = (n + BLOCK_SIZE - 1) / BLOCK_SIZE;
    kernelComputeSoundSpeed<<<blocks, BLOCK_SIZE>>>(temperature, gamma, n, sound_speed);
    CUDA_CHECK(cudaGetLastError());
#endif
}

// =============================================================================
// GPURadiationKernels Implementation
// =============================================================================

void GPURadiationKernels::computePlanckFunction(
    const float* temperature, const float* wavelength,
    int n_temp, int n_wave, float* radiance) {
#ifdef USE_CUDA
    dim3 blocks((n_temp + BLOCK_SIZE_2D - 1) / BLOCK_SIZE_2D,
                (n_wave + BLOCK_SIZE_2D - 1) / BLOCK_SIZE_2D);
    dim3 threads(BLOCK_SIZE_2D, BLOCK_SIZE_2D);
    
    kernelPlanckFunction<<<blocks, threads>>>(temperature, wavelength, n_temp, n_wave, radiance);
    CUDA_CHECK(cudaGetLastError());
#endif
}

void GPURadiationKernels::computeDecay(
    const float* activity, const float* half_life,
    float dt, int n_species, float* new_activity) {
#ifdef USE_CUDA
    int blocks = (n_species + BLOCK_SIZE - 1) / BLOCK_SIZE;
    kernelRadioactiveDecay<<<blocks, BLOCK_SIZE>>>(activity, half_life, dt, n_species, new_activity);
    CUDA_CHECK(cudaGetLastError());
#endif
}

void GPURadiationKernels::computeDoseRate(
    const float* activity, const float* dose_conversion,
    const float* distances, int n_sources, int n_points,
    float* dose_rate) {
#ifdef USE_CUDA
    int blocks = (n_points + BLOCK_SIZE - 1) / BLOCK_SIZE;
    kernelComputeDoseRate<<<blocks, BLOCK_SIZE>>>(
        activity, dose_conversion, nullptr, distances, n_sources, n_points, dose_rate);
    CUDA_CHECK(cudaGetLastError());
#endif
}

// =============================================================================
// GPUExplosionKernels Implementation
// =============================================================================

void GPUExplosionKernels::sedovTaylorBlast(
    float energy, float rho0, float gamma, float time,
    const float* radii, int n, float* pressure, float* density, float* velocity) {
#ifdef USE_CUDA
    int blocks = (n + BLOCK_SIZE - 1) / BLOCK_SIZE;
    kernelSedovTaylor<<<blocks, BLOCK_SIZE>>>(energy, rho0, gamma, time,
                                               radii, n, pressure, density, velocity);
    CUDA_CHECK(cudaGetLastError());
#endif
}

void GPUExplosionKernels::fireballDynamics(
    float yield_kt, float altitude, float time,
    float* radius, float* temperature, float* rise_velocity) {
#ifdef USE_CUDA
    kernelFireballDynamics<<<1, 1>>>(yield_kt, altitude, time,
                                     radius, temperature, rise_velocity);
    CUDA_CHECK(cudaGetLastError());
#endif
}

void GPUExplosionKernels::rayleighTaylorGrowth(
    float* interface_perturbation, float atwood_number,
    float g, float dt, int nx, int ny) {
#ifdef USE_CUDA
    dim3 blocks((nx + BLOCK_SIZE_2D - 1) / BLOCK_SIZE_2D,
                (ny + BLOCK_SIZE_2D - 1) / BLOCK_SIZE_2D);
    dim3 threads(BLOCK_SIZE_2D, BLOCK_SIZE_2D);
    
    kernelRayleighTaylor<<<blocks, threads>>>(interface_perturbation, atwood_number, g, dt, nx, ny);
    CUDA_CHECK(cudaGetLastError());
#endif
}

// =============================================================================
// GPUComputeContext Implementation
// =============================================================================

GPUComputeContext& GPUComputeContext::getInstance() {
    static GPUComputeContext instance;
    return instance;
}

GPUComputeContext::GPUComputeContext() {
    initialize();
}

GPUComputeContext::~GPUComputeContext() {
#ifdef USE_CUDA
    if (initialized_) {
        curandDestroyGenerator(curand_gen_);
        cudaEventDestroy(start_event_);
        cudaEventDestroy(stop_event_);
    }
#endif
}

bool GPUComputeContext::initialize() {
    if (initialized_) return true;
    
#ifdef USE_CUDA
    int device_count;
    cudaError_t err = cudaGetDeviceCount(&device_count);
    if (err != cudaSuccess || device_count == 0) {
        initialized_ = false;
        return false;
    }
    
    CUDA_CHECK(cudaSetDevice(0));
    CUDA_CHECK(curandCreateGenerator(&curand_gen_, CURAND_RNG_PSEUDO_DEFAULT));
    CUDA_CHECK(curandSetPseudoRandomGeneratorSeed(curand_gen_, 42));
    CUDA_CHECK(cudaEventCreate(&start_event_));
    CUDA_CHECK(cudaEventCreate(&stop_event_));
    
    initialized_ = true;
    return true;
#else
    return false;
#endif
}

bool GPUComputeContext::isGPUAvailable() const {
#ifdef USE_CUDA
    return initialized_;
#else
    return false;
#endif
}

void GPUComputeContext::generateUniform(float* data, size_t n, float low, float high) {
#ifdef USE_CUDA
    if (!initialized_) return;
    curandGenerateUniform(curand_gen_, data, n);
    
    // Scale to [low, high]
    float scale = high - low;
    int blocks = (n + BLOCK_SIZE - 1) / BLOCK_SIZE;
    kernelScale<<<blocks, BLOCK_SIZE>>>(data, scale, data, n);
    kernelAddScalar<<<blocks, BLOCK_SIZE>>>(data, low, data, n);
#endif
}

void GPUComputeContext::generateNormal(float* data, size_t n, float mean, float std) {
#ifdef USE_CUDA
    if (!initialized_) return;
    curandGenerateNormal(curand_gen_, data, n, mean, std);
#endif
}

void GPUComputeContext::startTimer() {
#ifdef USE_CUDA
    cudaEventRecord(start_event_);
#endif
}

float GPUComputeContext::stopTimer() {
#ifdef USE_CUDA
    cudaEventRecord(stop_event_);
    cudaEventSynchronize(stop_event_);
    float ms;
    cudaEventElapsedTime(&ms, start_event_, stop_event_);
    return ms;
#else
    return 0.0f;
#endif
}

// =============================================================================
// GPUUtils Implementation
// =============================================================================

namespace GPUUtils {

void add(const float* a, const float* b, float* c, size_t n) {
#ifdef USE_CUDA
    int blocks = (n + BLOCK_SIZE - 1) / BLOCK_SIZE;
    kernelAdd<<<blocks, BLOCK_SIZE>>>(a, b, c, n);
#endif
}

void relu(const float* in, float* out, size_t n) {
#ifdef USE_CUDA
    int blocks = (n + BLOCK_SIZE - 1) / BLOCK_SIZE;
    kernelRelu<<<blocks, BLOCK_SIZE>>>(in, out, n);
#endif
}

void gelu(const float* in, float* out, size_t n) {
#ifdef USE_CUDA
    int blocks = (n + BLOCK_SIZE - 1) / BLOCK_SIZE;
    kernelGelu<<<blocks, BLOCK_SIZE>>>(in, out, n);
#endif
}

void softmax(const float* in, float* out, int batch, int dim) {
#ifdef USE_CUDA
    kernelSoftmax<<<batch, 1>>>(in, out, batch, dim);
#endif
}

float sum(const float* data, size_t n) {
#ifdef USE_CUDA
    float* d_result;
    cudaMalloc(&d_result, sizeof(float));
    cudaMemset(d_result, 0, sizeof(float));
    
    int blocks = (n + BLOCK_SIZE - 1) / BLOCK_SIZE;
    kernelReduceSum<<<blocks, BLOCK_SIZE>>>(data, d_result, n);
    
    float result;
    cudaMemcpy(&result, d_result, sizeof(float), cudaMemcpyDeviceToHost);
    cudaFree(d_result);
    return result;
#else
    return 0.0f;
#endif
}

void complexMul(const void* a, const void* b, void* c, size_t n) {
#ifdef USE_CUDA
    int blocks = (n + BLOCK_SIZE - 1) / BLOCK_SIZE;
    kernelComplexMul<<<blocks, BLOCK_SIZE>>>(
        static_cast<const float2*>(a),
        static_cast<const float2*>(b),
        static_cast<float2*>(c), n);
#endif
}

} // namespace GPUUtils

} // namespace GPU
} // namespace FSRM
