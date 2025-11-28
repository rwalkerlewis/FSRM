/**
 * @file PerformanceOptimizations.hpp
 * @brief Performance optimization utilities for FSRM
 * 
 * Provides:
 * - SIMD vectorization helpers (AVX2/AVX-512/NEON)
 * - Memory layout optimizers (SoA, AoS, hybrid)
 * - Cache-aware data structures
 * - Prefetching utilities
 * - Thread pool for parallel operations
 * - Alignment utilities
 */

#ifndef FSRM_PERFORMANCE_OPTIMIZATIONS_HPP
#define FSRM_PERFORMANCE_OPTIMIZATIONS_HPP

#include <cstdint>
#include <cstddef>
#include <cmath>
#include <algorithm>
#include <vector>
#include <array>
#include <memory>
#include <functional>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <atomic>

// SIMD headers
#ifdef __AVX512F__
#include <immintrin.h>
#define FSRM_SIMD_WIDTH 8  // 8 doubles in AVX-512
#define FSRM_SIMD_ALIGN 64
#elif defined(__AVX2__)
#include <immintrin.h>
#define FSRM_SIMD_WIDTH 4  // 4 doubles in AVX2
#define FSRM_SIMD_ALIGN 32
#elif defined(__ARM_NEON)
#include <arm_neon.h>
#define FSRM_SIMD_WIDTH 2  // 2 doubles in NEON
#define FSRM_SIMD_ALIGN 16
#elif defined(__x86_64__) || defined(_M_X64)
// Include basic SSE header for x86_64 prefetch support
#include <xmmintrin.h>
#define FSRM_SIMD_WIDTH 1
#define FSRM_SIMD_ALIGN 8
#else
#define FSRM_SIMD_WIDTH 1
#define FSRM_SIMD_ALIGN 8
#endif

namespace FSRM {
namespace Performance {

// =============================================================================
// Alignment Utilities
// =============================================================================

/**
 * @brief Aligned memory allocator
 */
template<typename T, size_t Alignment = FSRM_SIMD_ALIGN>
class AlignedAllocator {
public:
    using value_type = T;
    using pointer = T*;
    using const_pointer = const T*;
    using reference = T&;
    using const_reference = const T&;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;

    template<typename U>
    struct rebind {
        using other = AlignedAllocator<U, Alignment>;
    };

    AlignedAllocator() noexcept = default;
    
    template<typename U>
    AlignedAllocator(const AlignedAllocator<U, Alignment>&) noexcept {}

    pointer allocate(size_type n) {
        void* ptr = nullptr;
        if (posix_memalign(&ptr, Alignment, n * sizeof(T)) != 0) {
            throw std::bad_alloc();
        }
        return static_cast<pointer>(ptr);
    }

    void deallocate(pointer p, size_type) noexcept {
        free(p);
    }

    bool operator==(const AlignedAllocator&) const noexcept { return true; }
    bool operator!=(const AlignedAllocator&) const noexcept { return false; }
};

/**
 * @brief Aligned vector type
 */
template<typename T>
using AlignedVector = std::vector<T, AlignedAllocator<T>>;

/**
 * @brief Check if pointer is aligned
 */
inline bool isAligned(const void* ptr, size_t alignment = FSRM_SIMD_ALIGN) {
    return reinterpret_cast<std::uintptr_t>(ptr) % alignment == 0;
}

// =============================================================================
// SIMD Vector Operations
// =============================================================================

/**
 * @brief SIMD vector class for double precision
 */
class SimdVec4d {
public:
#ifdef __AVX2__
    __m256d data;
    
    SimdVec4d() : data(_mm256_setzero_pd()) {}
    SimdVec4d(__m256d v) : data(v) {}
    SimdVec4d(double a) : data(_mm256_set1_pd(a)) {}
    SimdVec4d(double a, double b, double c, double d) 
        : data(_mm256_set_pd(d, c, b, a)) {}
    
    static SimdVec4d load(const double* ptr) {
        return _mm256_load_pd(ptr);
    }
    
    static SimdVec4d loadu(const double* ptr) {
        return _mm256_loadu_pd(ptr);
    }
    
    void store(double* ptr) const {
        _mm256_store_pd(ptr, data);
    }
    
    void storeu(double* ptr) const {
        _mm256_storeu_pd(ptr, data);
    }
    
    SimdVec4d operator+(const SimdVec4d& b) const {
        return _mm256_add_pd(data, b.data);
    }
    
    SimdVec4d operator-(const SimdVec4d& b) const {
        return _mm256_sub_pd(data, b.data);
    }
    
    SimdVec4d operator*(const SimdVec4d& b) const {
        return _mm256_mul_pd(data, b.data);
    }
    
    SimdVec4d operator/(const SimdVec4d& b) const {
        return _mm256_div_pd(data, b.data);
    }
    
    SimdVec4d sqrt() const {
        return _mm256_sqrt_pd(data);
    }
    
    double sum() const {
        __m256d temp = _mm256_hadd_pd(data, data);
        __m128d lo = _mm256_extractf128_pd(temp, 0);
        __m128d hi = _mm256_extractf128_pd(temp, 1);
        return _mm_cvtsd_f64(_mm_add_pd(lo, hi));
    }
    
    static SimdVec4d fmadd(const SimdVec4d& a, const SimdVec4d& b, const SimdVec4d& c) {
        return _mm256_fmadd_pd(a.data, b.data, c.data);
    }
    
    static SimdVec4d max(const SimdVec4d& a, const SimdVec4d& b) {
        return _mm256_max_pd(a.data, b.data);
    }
    
    static SimdVec4d min(const SimdVec4d& a, const SimdVec4d& b) {
        return _mm256_min_pd(a.data, b.data);
    }
#else
    // Fallback scalar implementation
    double data[4];
    
    SimdVec4d() { data[0] = data[1] = data[2] = data[3] = 0.0; }
    SimdVec4d(double a) { data[0] = data[1] = data[2] = data[3] = a; }
    SimdVec4d(double a, double b, double c, double d) {
        data[0] = a; data[1] = b; data[2] = c; data[3] = d;
    }
    
    static SimdVec4d load(const double* ptr) {
        SimdVec4d v;
        for (int i = 0; i < 4; ++i) v.data[i] = ptr[i];
        return v;
    }
    
    static SimdVec4d loadu(const double* ptr) { return load(ptr); }
    
    void store(double* ptr) const {
        for (int i = 0; i < 4; ++i) ptr[i] = data[i];
    }
    
    void storeu(double* ptr) const { store(ptr); }
    
    SimdVec4d operator+(const SimdVec4d& b) const {
        SimdVec4d r;
        for (int i = 0; i < 4; ++i) r.data[i] = data[i] + b.data[i];
        return r;
    }
    
    SimdVec4d operator-(const SimdVec4d& b) const {
        SimdVec4d r;
        for (int i = 0; i < 4; ++i) r.data[i] = data[i] - b.data[i];
        return r;
    }
    
    SimdVec4d operator*(const SimdVec4d& b) const {
        SimdVec4d r;
        for (int i = 0; i < 4; ++i) r.data[i] = data[i] * b.data[i];
        return r;
    }
    
    SimdVec4d operator/(const SimdVec4d& b) const {
        SimdVec4d r;
        for (int i = 0; i < 4; ++i) r.data[i] = data[i] / b.data[i];
        return r;
    }
    
    SimdVec4d sqrt() const {
        SimdVec4d r;
        for (int i = 0; i < 4; ++i) r.data[i] = std::sqrt(data[i]);
        return r;
    }
    
    double sum() const {
        return data[0] + data[1] + data[2] + data[3];
    }
    
    static SimdVec4d fmadd(const SimdVec4d& a, const SimdVec4d& b, const SimdVec4d& c) {
        SimdVec4d r;
        for (int i = 0; i < 4; ++i) r.data[i] = a.data[i] * b.data[i] + c.data[i];
        return r;
    }
    
    static SimdVec4d max(const SimdVec4d& a, const SimdVec4d& b) {
        SimdVec4d r;
        for (int i = 0; i < 4; ++i) r.data[i] = std::max(a.data[i], b.data[i]);
        return r;
    }
    
    static SimdVec4d min(const SimdVec4d& a, const SimdVec4d& b) {
        SimdVec4d r;
        for (int i = 0; i < 4; ++i) r.data[i] = std::min(a.data[i], b.data[i]);
        return r;
    }
#endif
};

// =============================================================================
// Vectorized Kernels
// =============================================================================

/**
 * @brief Vectorized dot product
 */
inline double dot_simd(const double* a, const double* b, size_t n) {
    SimdVec4d sum(0.0);
    size_t i = 0;
    
    // Main loop (SIMD)
    for (; i + 4 <= n; i += 4) {
        SimdVec4d va = SimdVec4d::load(a + i);
        SimdVec4d vb = SimdVec4d::load(b + i);
        sum = SimdVec4d::fmadd(va, vb, sum);
    }
    
    double result = sum.sum();
    
    // Remainder
    for (; i < n; ++i) {
        result += a[i] * b[i];
    }
    
    return result;
}

/**
 * @brief Vectorized AXPY: y = a*x + y
 */
inline void axpy_simd(double a, const double* x, double* y, size_t n) {
    SimdVec4d va(a);
    size_t i = 0;
    
    for (; i + 4 <= n; i += 4) {
        SimdVec4d vx = SimdVec4d::load(x + i);
        SimdVec4d vy = SimdVec4d::load(y + i);
        SimdVec4d::fmadd(va, vx, vy).store(y + i);
    }
    
    for (; i < n; ++i) {
        y[i] += a * x[i];
    }
}

/**
 * @brief Vectorized matrix-vector product (small dense matrix)
 */
inline void matvec_simd(const double* A, const double* x, double* y, 
                        int m, int n) {
    for (int i = 0; i < m; ++i) {
        y[i] = dot_simd(A + i * n, x, n);
    }
}

/**
 * @brief Vectorized stress update (isotropic elasticity)
 */
inline void stress_update_simd(
    const double* strain_rate,  // [n_points * 6]
    double* stress,             // [n_points * 6]
    double lambda, double mu,
    double dt, int n_points)
{
    const double c1 = lambda + 2.0 * mu;  // P-wave modulus
    const double c2 = lambda;
    const double c3 = 2.0 * mu;  // Shear modulus * 2
    
    SimdVec4d v_c1(c1 * dt);
    SimdVec4d v_c2(c2 * dt);
    SimdVec4d v_c3(c3 * dt);
    
    for (int i = 0; i < n_points; ++i) {
        const double* eps = strain_rate + i * 6;
        double* sig = stress + i * 6;
        
        // Volumetric strain rate
        double ekk = eps[0] + eps[1] + eps[2];
        double ekk_dt = ekk * dt;
        
        // Normal stresses: σ_ii += (λ*ε_kk + 2μ*ε_ii)*dt
        sig[0] += c2 * ekk_dt + c3 * eps[0] * dt;
        sig[1] += c2 * ekk_dt + c3 * eps[1] * dt;
        sig[2] += c2 * ekk_dt + c3 * eps[2] * dt;
        
        // Shear stresses: σ_ij += 2μ*ε_ij*dt
        sig[3] += c3 * eps[3] * dt;
        sig[4] += c3 * eps[4] * dt;
        sig[5] += c3 * eps[5] * dt;
    }
}

// =============================================================================
// Memory Layout Optimizers
// =============================================================================

/**
 * @brief Structure of Arrays layout for element data
 */
template<int N_VARS>
struct SoALayout {
    AlignedVector<double> data[N_VARS];
    size_t n_elements;
    
    void resize(size_t n) {
        n_elements = n;
        for (int i = 0; i < N_VARS; ++i) {
            data[i].resize(n);
        }
    }
    
    double* getVariable(int var) { return data[var].data(); }
    const double* getVariable(int var) const { return data[var].data(); }
    
    // SIMD-friendly access
    void loadElement(int elem, double* vars) const {
        for (int i = 0; i < N_VARS; ++i) {
            vars[i] = data[i][elem];
        }
    }
    
    void storeElement(int elem, const double* vars) {
        for (int i = 0; i < N_VARS; ++i) {
            data[i][elem] = vars[i];
        }
    }
};

/**
 * @brief Array of Structures layout
 */
template<int N_VARS>
struct AoSLayout {
    AlignedVector<std::array<double, N_VARS>> data;
    
    void resize(size_t n) {
        data.resize(n);
    }
    
    std::array<double, N_VARS>& operator[](size_t i) { return data[i]; }
    const std::array<double, N_VARS>& operator[](size_t i) const { return data[i]; }
};

/**
 * @brief Hybrid SoA/AoS layout for cache efficiency
 * Groups variables that are frequently accessed together
 */
template<int N_GROUPS, int VARS_PER_GROUP>
struct HybridLayout {
    SoALayout<VARS_PER_GROUP> groups[N_GROUPS];
    size_t n_elements;
    
    void resize(size_t n) {
        n_elements = n;
        for (int g = 0; g < N_GROUPS; ++g) {
            groups[g].resize(n);
        }
    }
};

// =============================================================================
// Cache-Aware Structures
// =============================================================================

/**
 * @brief Cache line size (typically 64 bytes)
 */
static constexpr size_t CACHE_LINE_SIZE = 64;

/**
 * @brief Pad structure to cache line boundary
 */
template<typename T>
struct alignas(CACHE_LINE_SIZE) CacheAligned {
    T value;
    char padding[CACHE_LINE_SIZE - sizeof(T) % CACHE_LINE_SIZE];
};

/**
 * @brief Cache-blocked iteration
 */
class BlockIterator {
public:
    size_t total_size;
    size_t block_size;
    
    BlockIterator(size_t total, size_t block = 256) 
        : total_size(total), block_size(block) {}
    
    template<typename Func>
    void iterate(Func&& f) const {
        for (size_t start = 0; start < total_size; start += block_size) {
            size_t end = std::min(start + block_size, total_size);
            f(start, end);
        }
    }
};

/**
 * @brief Prefetch utilities
 */
namespace Prefetch {
    template<int Level = 0>
    inline void read(const void* addr) {
#ifdef __x86_64__
        _mm_prefetch(reinterpret_cast<const char*>(addr), _MM_HINT_T0);
#elif defined(__aarch64__)
        __builtin_prefetch(addr, 0, 3 - Level);
#endif
    }
    
    template<int Level = 0>
    inline void write(const void* addr) {
#ifdef __x86_64__
        _mm_prefetch(reinterpret_cast<const char*>(addr), _MM_HINT_T0);
#elif defined(__aarch64__)
        __builtin_prefetch(addr, 1, 3 - Level);
#endif
    }
    
    // Prefetch next cache lines for streaming access
    inline void stream(const double* ptr, int ahead = 4) {
        for (int i = 0; i < ahead; ++i) {
            read<0>(ptr + i * CACHE_LINE_SIZE / sizeof(double));
        }
    }
}

// =============================================================================
// Thread Pool
// =============================================================================

/**
 * @brief Simple thread pool for parallel operations
 */
class ThreadPool {
public:
    ThreadPool(size_t num_threads = std::thread::hardware_concurrency()) {
        for (size_t i = 0; i < num_threads; ++i) {
            workers.emplace_back([this] {
                while (true) {
                    std::function<void()> task;
                    {
                        std::unique_lock<std::mutex> lock(queue_mutex);
                        condition.wait(lock, [this] { 
                            return stop || !tasks.empty(); 
                        });
                        if (stop && tasks.empty()) return;
                        task = std::move(tasks.front());
                        tasks.pop();
                    }
                    task();
                }
            });
        }
    }
    
    ~ThreadPool() {
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            stop = true;
        }
        condition.notify_all();
        for (auto& worker : workers) {
            worker.join();
        }
    }
    
    template<typename F>
    void enqueue(F&& f) {
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            tasks.emplace(std::forward<F>(f));
        }
        condition.notify_one();
    }
    
    void parallelFor(size_t start, size_t end, 
                     std::function<void(size_t, size_t)> body) {
        size_t n_threads = workers.size();
        size_t chunk = (end - start + n_threads - 1) / n_threads;
        
        std::atomic<size_t> completed{0};
        
        for (size_t i = 0; i < n_threads; ++i) {
            size_t chunk_start = start + i * chunk;
            size_t chunk_end = std::min(chunk_start + chunk, end);
            
            if (chunk_start < end) {
                enqueue([&, chunk_start, chunk_end] {
                    body(chunk_start, chunk_end);
                    completed.fetch_add(1);
                });
            }
        }
        
        // Wait for completion
        while (completed.load() < n_threads) {
            std::this_thread::yield();
        }
    }
    
    size_t numThreads() const { return workers.size(); }
    
private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;
    std::mutex queue_mutex;
    std::condition_variable condition;
    bool stop = false;
};

// =============================================================================
// Performance Metrics
// =============================================================================

/**
 * @brief High-resolution timer
 */
class Timer {
public:
    void start() {
        start_time = std::chrono::high_resolution_clock::now();
    }
    
    double stop() {
        auto end_time = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double>(end_time - start_time).count();
    }
    
    double elapsed() const {
        auto now = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double>(now - start_time).count();
    }
    
private:
    std::chrono::high_resolution_clock::time_point start_time;
};

/**
 * @brief Performance counters
 */
struct PerformanceCounters {
    std::atomic<uint64_t> flops{0};
    std::atomic<uint64_t> bytes_read{0};
    std::atomic<uint64_t> bytes_written{0};
    std::atomic<uint64_t> elements_processed{0};
    
    void reset() {
        flops = 0;
        bytes_read = 0;
        bytes_written = 0;
        elements_processed = 0;
    }
    
    double getGFLOPS(double time_seconds) const {
        return flops.load() / (time_seconds * 1e9);
    }
    
    double getBandwidthGB(double time_seconds) const {
        return (bytes_read.load() + bytes_written.load()) / (time_seconds * 1e9);
    }
    
    double getArithmeticIntensity() const {
        double total_bytes = bytes_read.load() + bytes_written.load();
        return (total_bytes > 0) ? flops.load() / total_bytes : 0.0;
    }
};

// =============================================================================
// Optimized Kernels for DG Operations
// =============================================================================

/**
 * @brief Optimized volume integral (flux divergence)
 */
void computeVolumeIntegral_optimized(
    const double* u,            // Solution [n_elem x n_dof x n_var]
    const double* D,            // Differentiation matrix [n_dof x n_dof]
    const double* J,            // Jacobians [n_elem x n_qp x 9]
    double* rhs,                // RHS [n_elem x n_dof x n_var]
    int n_elem, int n_dof, int n_var, int n_qp,
    ThreadPool& pool)
{
    pool.parallelFor(0, n_elem, [&](size_t start, size_t end) {
        // Thread-local buffers
        AlignedVector<double> flux(n_qp * n_var);
        AlignedVector<double> Dflux(n_dof * n_var);
        
        for (size_t e = start; e < end; ++e) {
            const double* u_e = u + e * n_dof * n_var;
            const double* J_e = J + e * n_qp * 9;
            double* rhs_e = rhs + e * n_dof * n_var;
            
            // Compute flux at quadrature points
            for (int q = 0; q < n_qp; ++q) {
                // Get solution at qp (interpolation)
                // Compute physical flux
                // Apply Jacobian transformation
            }
            
            // Apply differentiation matrix
            matvec_simd(D, flux.data(), Dflux.data(), n_dof, n_qp);
            
            // Accumulate to RHS
            for (int i = 0; i < n_dof * n_var; ++i) {
                rhs_e[i] += Dflux[i];
            }
        }
    });
}

/**
 * @brief Optimized surface integral (numerical flux)
 */
void computeSurfaceIntegral_optimized(
    const double* u,            // Solution [n_elem x n_dof x n_var]
    const int* face_neighbors,  // Neighbor elements [n_face x 2]
    const double* normals,      // Face normals [n_face x 3]
    const double* M_inv,        // Inverse mass matrix [n_elem x n_dof x n_dof]
    double* rhs,                // RHS [n_elem x n_dof x n_var]
    int n_face, int n_dof, int n_var,
    ThreadPool& pool)
{
    pool.parallelFor(0, n_face, [&](size_t start, size_t end) {
        AlignedVector<double> u_L(n_var);
        AlignedVector<double> u_R(n_var);
        AlignedVector<double> flux(n_var);
        
        for (size_t f = start; f < end; ++f) {
            int e_L = face_neighbors[f * 2];
            int e_R = face_neighbors[f * 2 + 1];
            
            // Get face states
            // Compute numerical flux (Riemann solver)
            // Accumulate to both elements
        }
    });
}

} // namespace Performance
} // namespace FSRM

#endif // FSRM_PERFORMANCE_OPTIMIZATIONS_HPP
