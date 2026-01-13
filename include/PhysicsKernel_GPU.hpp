#ifndef PHYSICS_KERNEL_GPU_HPP
#define PHYSICS_KERNEL_GPU_HPP

#include "PhysicsKernel.hpp"
#include "GPUManager.hpp"

#ifdef USE_CUDA
#include "GPUKernels.cuh"
#endif

namespace FSRM {

// =============================================================================
// GPU Kernel Base Infrastructure
// =============================================================================

/**
 * @brief CRTP base class for GPU-accelerated physics kernels
 * 
 * This template provides common GPU functionality through the Curiously 
 * Recurring Template Pattern (CRTP). Derived classes inherit GPU memory
 * management, execution control, and fallback behavior.
 * 
 * @tparam Derived The derived GPU kernel class
 * @tparam CPUBase The CPU base kernel class
 * 
 * ## Usage Example:
 * ```cpp
 * class MyKernelGPU : public GPUKernelBase<MyKernelGPU, MyKernel> {
 *     // Implement GPU-specific methods
 * };
 * ```
 */
template<typename Derived, typename CPUBase>
class GPUKernelBase : public CPUBase {
public:
    using Base = CPUBase;
    
    GPUKernelBase() : CPUBase(), gpu_initialized_(false), n_elements_gpu_(0) {
        // Set execution traits for GPU support
        this->setExecutionTraits(ExecutionTraits::GPUAccelerated(1000));
        this->setCapabilities(
            KernelCapability::ALL_CPU | 
            KernelCapability::GPU_RESIDUAL | 
            KernelCapability::GPU_JACOBIAN
        );
    }
    
    virtual ~GPUKernelBase() {
        if (gpu_initialized_) {
            static_cast<Derived*>(this)->freeGPUMemory();
        }
    }
    
    /// Check if GPU is initialized
    bool isGPUInitialized() const { return gpu_initialized_; }
    
    /// Get number of elements allocated on GPU
    int getGPUElementCount() const { return n_elements_gpu_; }
    
    /**
     * @brief Select execution backend based on problem size and GPU availability
     * 
     * @param n_elements Number of elements to process
     * @return Recommended execution backend
     */
    ExecutionBackend selectBackend(size_t n_elements) const override {
        if (!gpu_initialized_) {
            return ExecutionBackend::CPU;
        }
        if (n_elements >= this->getMinGPUElements()) {
            return ExecutionBackend::CUDA;
        }
        return ExecutionBackend::CPU;
    }
    
protected:
    bool gpu_initialized_;
    int n_elements_gpu_;
    
    /// Initialize GPU resources
    bool initializeGPU() {
        GPUManager& gpu = GPUManager::getInstance();
        if (!gpu.isAvailable()) {
            return false;
        }
        gpu_initialized_ = true;
        return true;
    }
    
    /// Get GPU manager reference
    GPUManager& getGPUManager() {
        return GPUManager::getInstance();
    }
};

/**
 * @brief GPU memory allocation helper with automatic cleanup
 * 
 * Manages a set of related GPU arrays that are allocated and freed together.
 */
class GPUMemoryBlock {
public:
    GPUMemoryBlock() = default;
    ~GPUMemoryBlock() { freeAll(); }
    
    // Disable copy
    GPUMemoryBlock(const GPUMemoryBlock&) = delete;
    GPUMemoryBlock& operator=(const GPUMemoryBlock&) = delete;
    
    /// Allocate a new GPU array and track it
    template<typename T>
    T* allocate(size_t count) {
        GPUManager& gpu = GPUManager::getInstance();
        T* ptr = static_cast<T*>(gpu.allocateDevice(count * sizeof(T)));
        if (ptr) {
            allocations_.push_back(ptr);
        }
        return ptr;
    }
    
    /// Free all tracked allocations
    void freeAll() {
        GPUManager& gpu = GPUManager::getInstance();
        for (void* ptr : allocations_) {
            if (ptr) gpu.freeDevice(ptr);
        }
        allocations_.clear();
    }
    
private:
    std::vector<void*> allocations_;
};

// =============================================================================
// GPU-Accelerated Kernel Implementations
// =============================================================================

/**
 * @brief GPU-accelerated single phase flow kernel
 * 
 * Implements Darcy flow on GPU with accumulation and flux computation.
 * Falls back to CPU implementation when GPU is unavailable or for small
 * problem sizes.
 * 
 * ## GPU Operations:
 * - Accumulation: φ·c_t·∂P/∂t (parallel over cells)
 * - Darcy flux: -k/μ·∇P (parallel over faces)
 * - Pressure update: (parallel over cells)
 */
class SinglePhaseFlowKernelGPU : public GPUKernelBase<SinglePhaseFlowKernelGPU, SinglePhaseFlowKernel> {
public:
    SinglePhaseFlowKernelGPU();
    ~SinglePhaseFlowKernelGPU() override;
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    /**
     * @brief Compute residual with automatic CPU/GPU dispatch
     * 
     * Uses GPU for large problems, CPU for small problems or when
     * GPU is unavailable.
     */
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    /**
     * @brief Compute Jacobian (always uses CPU base class)
     * 
     * Pointwise Jacobian evaluation is CPU-based. GPU acceleration
     * applies to matrix assembly, not individual point evaluation.
     */
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    // =========================================================================
    // GPU Memory Management
    // =========================================================================
    
    /// Allocate GPU memory for given problem size
    void allocateGPUMemory(int n_cells);
    
    /// Free all GPU memory
    void freeGPUMemory();
    
    /// Copy data from host to GPU
    void copyDataToGPU(const double* host_data, size_t size);
    
    /// Copy data from GPU to host
    void copyDataFromGPU(double* host_data, size_t size);
    
    // =========================================================================
    // Batch GPU Operations
    // =========================================================================
    
    /// Compute accumulation term on GPU for all cells
    void computeAccumulationGPU(double dt);
    
    /// Compute Darcy fluxes on GPU for all faces
    void computeFluxGPU();
    
private:
#ifdef USE_CUDA
    // GPU arrays (device pointers)
    double* d_pressure;
    double* d_pressure_old;
    double* d_porosity;
    double* d_permeability;
    double* d_compressibility;
    double* d_viscosity;
    double* d_accumulation;
    double* d_flux;
#endif
};

/**
 * @brief GPU-accelerated elastodynamics kernel
 * 
 * Implements elastic wave propagation on GPU with full inertial dynamics.
 * Suitable for seismic wave modeling, earthquake simulation, and 
 * dynamic fracture problems.
 * 
 * ## Governing Equation:
 * ρ·∂²u/∂t² = ∇·σ + f - C·∂u/∂t
 * 
 * ## GPU Operations:
 * - Strain computation from displacement gradients
 * - Stress computation via Hooke's law
 * - Stress divergence for force balance
 * - Velocity/displacement updates
 * - Rayleigh damping application
 * 
 * ## Performance:
 * - 5-20× speedup over CPU for >10K cells
 * - Optimal for structured grids with regular connectivity
 */
class ElastodynamicsKernelGPU : public GPUKernelBase<ElastodynamicsKernelGPU, ElastodynamicsKernel> {
public:
    ElastodynamicsKernelGPU();
    ~ElastodynamicsKernelGPU() override;
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    /**
     * @brief Compute residual with GPU acceleration for large problems
     */
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    /**
     * @brief Compute Jacobian (full stiffness + mass matrix with damping)
     */
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    // =========================================================================
    // GPU Memory Management
    // =========================================================================
    
    void allocateGPUMemory(int n_cells);
    void freeGPUMemory();
    
    // =========================================================================
    // GPU Batch Operations
    // =========================================================================
    
    /**
     * @brief Execute full elastodynamics timestep on GPU
     * 
     * @param displacement Current displacement field
     * @param velocity Current velocity field
     * @param force Output force vector
     * @param dt Time step size
     * 
     * Computes: strain → stress → ∇·σ → force (with damping)
     */
    void computeOnGPU(const double* displacement, const double* velocity,
                     double* force, double dt);
    
    /// Compute strain tensor on GPU
    void computeStrainGPU();
    
    /// Compute stress tensor on GPU
    void computeStressGPU();
    
    /// Compute stress divergence on GPU
    void computeStressDivergenceGPU();
    
    /// Apply Rayleigh damping on GPU
    void applyDampingGPU(double dt);
    
private:
#ifdef USE_CUDA
    // Solution arrays
    double* d_displacement;     ///< Displacement field (n_cells × 3)
    double* d_velocity;         ///< Velocity field (n_cells × 3)
    
    // Derived quantities
    double* d_strain;           ///< Strain tensor (n_cells × 6, Voigt notation)
    double* d_stress;           ///< Stress tensor (n_cells × 6, Voigt notation)
    double* d_stress_div;       ///< Stress divergence (n_cells × 3)
    double* d_force;            ///< Force vector (n_cells × 3)
    
    // Material properties
    double* d_youngs_modulus;   ///< Young's modulus per cell
    double* d_poisson_ratio;    ///< Poisson's ratio per cell
    double* d_density;          ///< Density per cell
    
    // Mesh connectivity
    int* d_connectivity;        ///< Cell connectivity for gradient computation
#endif
};

/**
 * @brief GPU-accelerated poroelastodynamics kernel
 * 
 * Implements Biot's coupled fluid-solid wave propagation on GPU.
 * Captures three wave types: fast P-wave, slow P-wave (Biot wave II), 
 * and S-wave. Includes dynamic permeability enhancement from wave passage.
 * 
 * ## Governing Equations (Biot Theory):
 * 
 * Momentum balance:
 *   ρ_bulk·∂²u/∂t² = ∇·σ' - α·∇p + damping
 * 
 * Mass conservation:
 *   (1/M)·∂p/∂t + α·∇·(∂u/∂t) + ∇·q = 0
 * 
 * Where:
 * - σ' = effective stress (drained)
 * - α = Biot coefficient
 * - M = Biot modulus
 * - q = Darcy flux = -(k/μ)·∇p
 * 
 * ## GPU Operations:
 * - Coupled solid-fluid momentum
 * - Fluid mass balance with storage
 * - Dynamic permeability updates
 * - Undrained modulus for fast waves
 * 
 * ## Performance:
 * - 5-50× speedup over CPU for >10K cells
 * - Essential for induced seismicity modeling
 */
class PoroelastodynamicsKernelGPU : public GPUKernelBase<PoroelastodynamicsKernelGPU, PoroelastodynamicsKernel> {
public:
    PoroelastodynamicsKernelGPU();
    ~PoroelastodynamicsKernelGPU() override;
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    /**
     * @brief Compute coupled residual with GPU acceleration
     */
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    /**
     * @brief Compute coupled Jacobian (includes Biot coupling terms)
     */
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    // =========================================================================
    // GPU Memory Management
    // =========================================================================
    
    void allocateGPUMemory(int n_cells);
    void freeGPUMemory();
    
    // =========================================================================
    // GPU Batch Operations
    // =========================================================================
    
    /**
     * @brief Execute full poroelastodynamics timestep on GPU
     * 
     * @param displacement Current displacement field
     * @param velocity Current velocity field
     * @param pressure Current pressure field
     * @param force Output solid force vector
     * @param pressure_rate Output pressure time derivative
     * @param dt Time step size
     * 
     * Computes coupled solid-fluid response including:
     * - Effective stress with pore pressure coupling
     * - Fluid storage and Darcy flow
     * - Biot coupling terms
     */
    void computeOnGPU(const double* displacement, const double* velocity,
                     const double* pressure, double* force, double* pressure_rate,
                     double dt);
    
    /**
     * @brief Update permeability based on strain/stress on GPU
     * 
     * @param dt Time step for recovery calculation
     * 
     * Implements dynamic permeability enhancement from wave passage:
     * - Volumetric strain model: k = k₀·exp(β·ε_vol)
     * - Shear dilation model
     * - Time-dependent recovery to initial value
     */
    void updatePermeabilityGPU(double dt);
    
    /// Compute Biot coupling terms on GPU
    void computeBiotCouplingGPU();
    
    /// Compute Darcy flux on GPU
    void computeDarcyFluxGPU();
    
private:
#ifdef USE_CUDA
    // =========================================================================
    // Solid Mechanics Arrays
    // =========================================================================
    double* d_displacement;         ///< Displacement (n_cells × 3)
    double* d_velocity;             ///< Velocity (n_cells × 3)
    double* d_strain;               ///< Strain tensor (n_cells × 6)
    double* d_stress;               ///< Effective stress (n_cells × 6)
    double* d_force_solid;          ///< Solid force (n_cells × 3)
    
    // =========================================================================
    // Fluid Flow Arrays
    // =========================================================================
    double* d_pressure;             ///< Pore pressure (n_cells)
    double* d_pressure_rate;        ///< Pressure rate (n_cells)
    double* d_permeability;         ///< Current permeability (n_cells)
    double* d_permeability_initial; ///< Initial permeability (n_cells)
    double* d_darcy_flux;           ///< Darcy flux (n_faces × 1)
    
    // =========================================================================
    // Material Properties
    // =========================================================================
    double* d_youngs_modulus;       ///< Young's modulus (n_cells)
    double* d_poisson_ratio;        ///< Poisson's ratio (n_cells)
    double* d_biot_coefficient;     ///< Biot coefficient α (n_cells)
    double* d_biot_modulus;         ///< Biot modulus M (n_cells)
    double* d_porosity;             ///< Porosity φ (n_cells)
    double* d_density_solid;        ///< Solid density (n_cells)
    double* d_density_fluid;        ///< Fluid density (n_cells)
    double* d_viscosity;            ///< Fluid viscosity (n_cells)
    
    // =========================================================================
    // Mesh Data
    // =========================================================================
    int* d_connectivity;            ///< Cell connectivity
    int* d_face_cells;              ///< Face-to-cell mapping
#endif
};

// =============================================================================
// GPU Kernel Factory
// =============================================================================

/**
 * @brief Create GPU-accelerated kernel for specified physics type
 * 
 * Factory function that creates the appropriate GPU kernel implementation.
 * Returns nullptr if GPU is not available or physics type doesn't have
 * a GPU implementation.
 * 
 * @param type Physics type to create kernel for
 * @return Shared pointer to GPU kernel, or nullptr if unavailable
 * 
 * ## Supported GPU Kernels:
 * - FLUID_FLOW → SinglePhaseFlowKernelGPU
 * - ELASTODYNAMICS → ElastodynamicsKernelGPU
 * - POROELASTODYNAMICS → PoroelastodynamicsKernelGPU
 */
std::shared_ptr<PhysicsKernel> createGPUKernel(PhysicsType type);

/**
 * @brief Create kernel with automatic CPU/GPU selection
 * 
 * @param type Physics type
 * @param mode Execution mode preference
 * @return Appropriate kernel (GPU if available and requested, else CPU)
 */
std::shared_ptr<PhysicsKernel> createKernel(PhysicsType type, 
                                            GPUExecutionMode mode = GPUExecutionMode::CPU_FALLBACK);

/**
 * @brief Check if GPU kernel is available for physics type
 * @param type Physics type to check
 * @return true if GPU implementation exists and GPU is available
 */
bool hasGPUKernel(PhysicsType type);

// =============================================================================
// GPU Data Transfer Utilities
// =============================================================================

/**
 * @brief RAII wrapper for GPU-CPU data transfers
 * 
 * Manages device memory allocation and provides convenient
 * copy operations. Memory is automatically freed on destruction.
 */
class GPUDataTransfer {
public:
    /**
     * @brief Allocate device memory
     * @param size Size in bytes
     */
    explicit GPUDataTransfer(size_t size);
    ~GPUDataTransfer();
    
    // Disable copy
    GPUDataTransfer(const GPUDataTransfer&) = delete;
    GPUDataTransfer& operator=(const GPUDataTransfer&) = delete;
    
    // Enable move
    GPUDataTransfer(GPUDataTransfer&& other) noexcept;
    GPUDataTransfer& operator=(GPUDataTransfer&& other) noexcept;
    
    /// Copy data from host to device
    void copyToGPU(const void* host_data);
    
    /// Copy data from device to host
    void copyFromGPU(void* host_data);
    
    /// Async copy to GPU (returns immediately)
    void copyToGPUAsync(const void* host_data, void* stream = nullptr);
    
    /// Async copy from GPU (returns immediately)
    void copyFromGPUAsync(void* host_data, void* stream = nullptr);
    
    /// Get device pointer
    void* getDevicePointer() { return device_ptr_; }
    const void* getDevicePointer() const { return device_ptr_; }
    
    /// Get allocation size
    size_t getSize() const { return size_bytes_; }
    
    /// Check if allocated
    bool isValid() const { return device_ptr_ != nullptr; }
    
private:
    void* device_ptr_;
    size_t size_bytes_;
};

// =============================================================================
// GPU Performance Benchmarking
// =============================================================================

/**
 * @brief GPU vs CPU performance benchmarking utility
 * 
 * Provides standardized benchmarks for comparing GPU and CPU
 * kernel performance. Useful for determining optimal execution
 * backend for different problem sizes.
 */
class GPUBenchmark {
public:
    /// Benchmark single-phase flow kernel
    static void runSinglePhaseFlowBenchmark(int n_cells, int n_iterations);
    
    /// Benchmark elastodynamics kernel
    static void runElastodynamicsBenchmark(int n_cells, int n_iterations);
    
    /// Benchmark poroelastodynamics kernel
    static void runPoroelastodynamicsBenchmark(int n_cells, int n_iterations);
    
    /// Print formatted benchmark results
    static void printBenchmarkResults();
    
    /// Get last CPU time in milliseconds
    static double getCPUTime() { return cpu_time_ms_; }
    
    /// Get last GPU time in milliseconds
    static double getGPUTime() { return gpu_time_ms_; }
    
    /// Get speedup factor (CPU_time / GPU_time)
    static double getSpeedup() { return speedup_factor_; }
    
    /**
     * @brief Estimate optimal problem size threshold for GPU
     * @param type Physics kernel type
     * @return Minimum elements for GPU to be beneficial
     */
    static size_t estimateGPUThreshold(PhysicsType type);
    
private:
    static double cpu_time_ms_;
    static double gpu_time_ms_;
    static double speedup_factor_;
};

// =============================================================================
// GPU Execution Statistics
// =============================================================================

/**
 * @brief Runtime statistics for GPU kernel execution
 * 
 * Tracks execution time, memory usage, and throughput for
 * performance monitoring and optimization.
 */
struct GPUExecutionStats {
    double total_gpu_time_ms = 0.0;
    double total_cpu_time_ms = 0.0;
    size_t gpu_calls = 0;
    size_t cpu_fallback_calls = 0;
    size_t total_elements_processed = 0;
    size_t peak_gpu_memory_bytes = 0;
    
    /// Reset all statistics
    void reset() {
        total_gpu_time_ms = total_cpu_time_ms = 0.0;
        gpu_calls = cpu_fallback_calls = 0;
        total_elements_processed = peak_gpu_memory_bytes = 0;
    }
    
    /// Get average time per GPU call
    double avgGPUTime() const { 
        return gpu_calls > 0 ? total_gpu_time_ms / gpu_calls : 0.0; 
    }
    
    /// Get GPU utilization fraction
    double gpuUtilization() const {
        size_t total = gpu_calls + cpu_fallback_calls;
        return total > 0 ? static_cast<double>(gpu_calls) / total : 0.0;
    }
    
    /// Get throughput (elements per millisecond)
    double throughput() const {
        double total_time = total_gpu_time_ms + total_cpu_time_ms;
        return total_time > 0 ? total_elements_processed / total_time : 0.0;
    }
};

/// Global execution statistics (thread-local for safety)
GPUExecutionStats& getGPUStats();

} // namespace FSRM

#endif // PHYSICS_KERNEL_GPU_HPP
