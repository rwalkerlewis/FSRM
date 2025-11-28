/**
 * @file PerformanceOptimizations.cpp
 * @brief Implementation of performance optimization utilities
 */

#include "PerformanceOptimizations.hpp"
#include <algorithm>
#include <numeric>
#include <cstring>
#include <iostream>

namespace FSRM {
namespace Performance {

// =============================================================================
// Optimized Basis Function Evaluation
// =============================================================================

/**
 * @brief Batch evaluate Lagrange basis at multiple points
 */
void evaluateLagrangeBasis_batch(
    const double* nodes,     // [n_nodes]
    const double* points,    // [n_points]
    double* values,          // [n_points x n_nodes]
    int n_nodes, int n_points)
{
    // Precompute denominators
    AlignedVector<double> denom(n_nodes);
    for (int i = 0; i < n_nodes; ++i) {
        double d = 1.0;
        for (int j = 0; j < n_nodes; ++j) {
            if (i != j) {
                d *= (nodes[i] - nodes[j]);
            }
        }
        denom[i] = 1.0 / d;
    }
    
    // Evaluate at all points
    for (int p = 0; p < n_points; ++p) {
        double x = points[p];
        
        // Compute numerator product
        double prod = 1.0;
        for (int j = 0; j < n_nodes; ++j) {
            prod *= (x - nodes[j]);
        }
        
        // Compute each basis function
        for (int i = 0; i < n_nodes; ++i) {
            double diff = x - nodes[i];
            if (std::abs(diff) < 1e-14) {
                // x is at node i
                for (int k = 0; k < n_nodes; ++k) {
                    values[p * n_nodes + k] = (k == i) ? 1.0 : 0.0;
                }
                break;
            } else {
                values[p * n_nodes + i] = prod / diff * denom[i];
            }
        }
    }
}

/**
 * @brief SIMD-optimized Legendre polynomial evaluation
 */
void evaluateLegendre_simd(
    const double* x,     // [n_points]
    double* P,           // [n_points x (order+1)]
    int n_points, int order)
{
    // P_0 = 1, P_1 = x
    for (int i = 0; i < n_points; i += 4) {
        SimdVec4d one(1.0);
        SimdVec4d xi = SimdVec4d::load(x + i);
        
        one.store(P + i);
        if (order >= 1) {
            xi.store(P + n_points + i);
        }
        
        // Bonnet's recursion: (n+1)P_{n+1} = (2n+1)xP_n - nP_{n-1}
        for (int n = 1; n < order; ++n) {
            SimdVec4d Pn = SimdVec4d::load(P + n * n_points + i);
            SimdVec4d Pnm1 = SimdVec4d::load(P + (n-1) * n_points + i);
            
            double a = (2.0 * n + 1.0) / (n + 1.0);
            double b = (double)n / (n + 1.0);
            
            SimdVec4d va(a);
            SimdVec4d vb(b);
            
            SimdVec4d Pnp1 = va * xi * Pn - vb * Pnm1;
            Pnp1.store(P + (n+1) * n_points + i);
        }
    }
}

// =============================================================================
// Optimized Matrix Operations
// =============================================================================

/**
 * @brief Small matrix-matrix multiply (optimized for DG)
 */
void matmul_small(
    const double* A,  // [m x k]
    const double* B,  // [k x n]
    double* C,        // [m x n]
    int m, int k, int n)
{
    // Clear output
    std::memset(C, 0, m * n * sizeof(double));
    
    // Cache-blocked multiplication
    constexpr int BLOCK = 8;
    
    for (int i0 = 0; i0 < m; i0 += BLOCK) {
        for (int j0 = 0; j0 < n; j0 += BLOCK) {
            for (int p0 = 0; p0 < k; p0 += BLOCK) {
                int i_end = std::min(i0 + BLOCK, m);
                int j_end = std::min(j0 + BLOCK, n);
                int p_end = std::min(p0 + BLOCK, k);
                
                for (int i = i0; i < i_end; ++i) {
                    for (int p = p0; p < p_end; ++p) {
                        double a_ip = A[i * k + p];
                        for (int j = j0; j < j_end; ++j) {
                            C[i * n + j] += a_ip * B[p * n + j];
                        }
                    }
                }
            }
        }
    }
}

/**
 * @brief Batched matrix-vector multiply
 */
void batched_matvec(
    const double* A,     // [batch x m x n]
    const double* x,     // [batch x n]
    double* y,           // [batch x m]
    int batch, int m, int n,
    ThreadPool& pool)
{
    pool.parallelFor(0, batch, [&](size_t start, size_t end) {
        for (size_t b = start; b < end; ++b) {
            const double* A_b = A + b * m * n;
            const double* x_b = x + b * n;
            double* y_b = y + b * m;
            
            for (int i = 0; i < m; ++i) {
                y_b[i] = dot_simd(A_b + i * n, x_b, n);
            }
        }
    });
}

/**
 * @brief In-place matrix transpose
 */
void transpose_inplace(double* A, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            std::swap(A[i * n + j], A[j * n + i]);
        }
    }
}

// =============================================================================
// Optimized Riemann Solver
// =============================================================================

/**
 * @brief Batch Rusanov flux computation
 */
void rusanov_flux_batch(
    const double* u_L,      // [n_face x n_var]
    const double* u_R,      // [n_face x n_var]
    const double* normals,  // [n_face x 3]
    const double* vp,       // [n_face] P-wave speeds
    const double* vs,       // [n_face] S-wave speeds
    double* flux,           // [n_face x n_var]
    int n_face, int n_var)
{
    // Process 4 faces at a time
    for (int f = 0; f < n_face; f += 4) {
        int n_batch = std::min(4, n_face - f);
        
        for (int b = 0; b < n_batch; ++b) {
            int idx = f + b;
            const double* uL = u_L + idx * n_var;
            const double* uR = u_R + idx * n_var;
            double* F = flux + idx * n_var;
            
            // Maximum wave speed
            double lambda_max = vp[idx];
            
            // Rusanov: F = 0.5*(F_L + F_R) - 0.5*lambda_max*(u_R - u_L)
            for (int v = 0; v < n_var; ++v) {
                // Central flux (simplified)
                double F_avg = 0.5 * (uL[v] + uR[v]);
                // Dissipation
                double diss = 0.5 * lambda_max * (uR[v] - uL[v]);
                F[v] = F_avg - diss;
            }
        }
    }
}

// =============================================================================
// Memory Pool
// =============================================================================

/**
 * @brief Simple memory pool for temporary allocations
 */
class MemoryPool {
public:
    explicit MemoryPool(size_t size = 64 * 1024 * 1024)  // 64 MB default
        : pool_size(size), offset(0)
    {
        pool = static_cast<char*>(aligned_alloc(CACHE_LINE_SIZE, size));
    }
    
    ~MemoryPool() {
        free(pool);
    }
    
    void* allocate(size_t size, size_t alignment = CACHE_LINE_SIZE) {
        // Align offset
        size_t aligned_offset = (offset + alignment - 1) & ~(alignment - 1);
        
        if (aligned_offset + size > pool_size) {
            throw std::bad_alloc();
        }
        
        void* ptr = pool + aligned_offset;
        offset = aligned_offset + size;
        return ptr;
    }
    
    template<typename T>
    T* allocate(size_t n) {
        return static_cast<T*>(allocate(n * sizeof(T), alignof(T)));
    }
    
    void reset() {
        offset = 0;
    }
    
    size_t used() const { return offset; }
    size_t available() const { return pool_size - offset; }
    
private:
    char* pool;
    size_t pool_size;
    size_t offset;
};

// Global memory pool (thread-local for safety)
thread_local MemoryPool g_temp_pool;

// =============================================================================
// Optimized Stress-Strain Operations
// =============================================================================

/**
 * @brief Compute strain from velocity gradient (batch)
 */
void compute_strain_batch(
    const double* grad_v,   // [n_points x 9] velocity gradients
    double* strain_rate,    // [n_points x 6] symmetric strain rate
    int n_points)
{
    for (int i = 0; i < n_points; i += 4) {
        // Load velocity gradients for 4 points
        // grad_v layout: [dvx/dx, dvx/dy, dvx/dz, dvy/dx, dvy/dy, dvy/dz, dvz/dx, dvz/dy, dvz/dz]
        
        for (int b = 0; b < std::min(4, n_points - i); ++b) {
            const double* g = grad_v + (i + b) * 9;
            double* e = strain_rate + (i + b) * 6;
            
            // Symmetric strain rate
            e[0] = g[0];                          // e_xx = dv_x/dx
            e[1] = g[4];                          // e_yy = dv_y/dy
            e[2] = g[8];                          // e_zz = dv_z/dz
            e[3] = 0.5 * (g[1] + g[3]);           // e_xy = 0.5*(dv_x/dy + dv_y/dx)
            e[4] = 0.5 * (g[2] + g[6]);           // e_xz = 0.5*(dv_x/dz + dv_z/dx)
            e[5] = 0.5 * (g[5] + g[7]);           // e_yz = 0.5*(dv_y/dz + dv_z/dy)
        }
    }
}

/**
 * @brief Isotropic stress update (batch)
 */
void update_stress_isotropic_batch(
    const double* strain_rate,  // [n_points x 6]
    double* stress,             // [n_points x 6]
    double lambda, double mu, double dt,
    int n_points)
{
    const double c1 = (lambda + 2.0 * mu) * dt;
    const double c2 = lambda * dt;
    const double c3 = 2.0 * mu * dt;
    
    SimdVec4d v_c1(c1);
    SimdVec4d v_c2(c2);
    SimdVec4d v_c3(c3);
    
    for (int i = 0; i < n_points; ++i) {
        const double* e = strain_rate + i * 6;
        double* s = stress + i * 6;
        
        double trace = e[0] + e[1] + e[2];
        
        s[0] += c2 * trace + c3 * e[0];
        s[1] += c2 * trace + c3 * e[1];
        s[2] += c2 * trace + c3 * e[2];
        s[3] += c3 * e[3];
        s[4] += c3 * e[4];
        s[5] += c3 * e[5];
    }
}

/**
 * @brief Anisotropic stress update (batch)
 */
void update_stress_anisotropic_batch(
    const double* strain_rate,  // [n_points x 6]
    const double* C,            // [n_points x 36] or [36] if uniform
    double* stress,             // [n_points x 6]
    double dt, int n_points, bool uniform_C)
{
    if (uniform_C) {
        // Same stiffness for all points
        for (int i = 0; i < n_points; ++i) {
            const double* e = strain_rate + i * 6;
            double* s = stress + i * 6;
            
            for (int j = 0; j < 6; ++j) {
                for (int k = 0; k < 6; ++k) {
                    s[j] += C[j * 6 + k] * e[k] * dt;
                }
            }
        }
    } else {
        // Point-varying stiffness
        for (int i = 0; i < n_points; ++i) {
            const double* e = strain_rate + i * 6;
            const double* Ci = C + i * 36;
            double* s = stress + i * 6;
            
            for (int j = 0; j < 6; ++j) {
                s[j] += dot_simd(Ci + j * 6, e, 6) * dt;
            }
        }
    }
}

// =============================================================================
// Data Layout Conversion
// =============================================================================

/**
 * @brief Convert AoS to SoA
 */
void aos_to_soa(
    const double* aos,  // [n_elem x n_var]
    double** soa,       // [n_var][n_elem]
    int n_elem, int n_var)
{
    for (int v = 0; v < n_var; ++v) {
        for (int e = 0; e < n_elem; ++e) {
            soa[v][e] = aos[e * n_var + v];
        }
    }
}

/**
 * @brief Convert SoA to AoS
 */
void soa_to_aos(
    double** soa,       // [n_var][n_elem]
    double* aos,        // [n_elem x n_var]
    int n_elem, int n_var)
{
    for (int e = 0; e < n_elem; ++e) {
        for (int v = 0; v < n_var; ++v) {
            aos[e * n_var + v] = soa[v][e];
        }
    }
}

// =============================================================================
// Profiling Utilities
// =============================================================================

/**
 * @brief Profiling results
 */
struct ProfileResult {
    std::string name;
    double total_time;
    double min_time;
    double max_time;
    int call_count;
    
    double avgTime() const { return total_time / call_count; }
};

/**
 * @brief Simple profiler
 */
class Profiler {
public:
    void start(const std::string& name) {
        auto& entry = entries[name];
        entry.name = name;
        entry.timer.start();
    }
    
    void stop(const std::string& name) {
        auto& entry = entries[name];
        double t = entry.timer.stop();
        entry.result.total_time += t;
        entry.result.min_time = std::min(entry.result.min_time, t);
        entry.result.max_time = std::max(entry.result.max_time, t);
        entry.result.call_count++;
    }
    
    void report(std::ostream& os = std::cout) const {
        os << "\n=== Profiling Results ===" << std::endl;
        os << std::setw(30) << "Name" 
           << std::setw(12) << "Total (s)"
           << std::setw(12) << "Avg (ms)"
           << std::setw(12) << "Calls"
           << std::endl;
        os << std::string(66, '-') << std::endl;
        
        for (const auto& [name, entry] : entries) {
            os << std::setw(30) << name
               << std::setw(12) << std::fixed << std::setprecision(3) 
               << entry.result.total_time
               << std::setw(12) << entry.result.avgTime() * 1000
               << std::setw(12) << entry.result.call_count
               << std::endl;
        }
    }
    
    void reset() {
        entries.clear();
    }
    
private:
    struct Entry {
        Timer timer;
        ProfileResult result{};
    };
    std::map<std::string, Entry> entries;
};

// Global profiler instance
Profiler g_profiler;

// =============================================================================
// Optimized Sorting for LTS
// =============================================================================

/**
 * @brief Stable bucket sort for LTS cluster assignment
 */
void bucket_sort_elements(
    const double* dt_local,    // [n_elem] local timestep per element
    int* cluster_id,           // [n_elem] output cluster assignment
    int n_elem, int n_clusters, double dt_min)
{
    // Bucket boundaries (powers of 2)
    std::vector<double> boundaries(n_clusters + 1);
    boundaries[0] = dt_min;
    for (int c = 1; c <= n_clusters; ++c) {
        boundaries[c] = dt_min * (1 << c);
    }
    
    // Assign elements to clusters
    for (int e = 0; e < n_elem; ++e) {
        double dt = dt_local[e];
        
        // Binary search for cluster
        int c = 0;
        for (c = 0; c < n_clusters; ++c) {
            if (dt < boundaries[c + 1]) break;
        }
        cluster_id[e] = std::min(c, n_clusters - 1);
    }
}

/**
 * @brief Reorder elements by cluster for cache efficiency
 */
void reorder_by_cluster(
    const int* cluster_id,     // [n_elem]
    int* permutation,          // [n_elem] output permutation
    int* inverse_perm,         // [n_elem] output inverse permutation
    int n_elem, int n_clusters)
{
    // Count elements per cluster
    std::vector<int> counts(n_clusters, 0);
    for (int e = 0; e < n_elem; ++e) {
        counts[cluster_id[e]]++;
    }
    
    // Compute offsets
    std::vector<int> offsets(n_clusters + 1);
    offsets[0] = 0;
    for (int c = 0; c < n_clusters; ++c) {
        offsets[c + 1] = offsets[c] + counts[c];
    }
    
    // Build permutation
    std::vector<int> current_offset = offsets;
    for (int e = 0; e < n_elem; ++e) {
        int c = cluster_id[e];
        int new_pos = current_offset[c]++;
        permutation[new_pos] = e;
        inverse_perm[e] = new_pos;
    }
}

} // namespace Performance
} // namespace FSRM
