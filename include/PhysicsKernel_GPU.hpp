#ifndef PHYSICS_KERNEL_GPU_HPP
#define PHYSICS_KERNEL_GPU_HPP

#include "PhysicsKernel.hpp"
#include "GPUManager.hpp"

#ifdef USE_CUDA
#include "GPUKernels.cuh"
#endif

namespace FSRM {

// GPU-accelerated single phase flow kernel
class SinglePhaseFlowKernelGPU : public SinglePhaseFlowKernel {
public:
    SinglePhaseFlowKernelGPU();
    ~SinglePhaseFlowKernelGPU() override;
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    // GPU-specific methods
    void allocateGPUMemory(int n_cells);
    void freeGPUMemory();
    void copyDataToGPU(const double* host_data, size_t size);
    void copyDataFromGPU(double* host_data, size_t size);
    
private:
    bool gpu_initialized;
    int n_cells_gpu;
    
#ifdef USE_CUDA
    // GPU arrays
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

// GPU-accelerated elastodynamics kernel
class ElastodynamicsKernelGPU : public ElastodynamicsKernel {
public:
    ElastodynamicsKernelGPU();
    ~ElastodynamicsKernelGPU() override;
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    // GPU-specific methods
    void allocateGPUMemory(int n_cells);
    void freeGPUMemory();
    void computeOnGPU(const double* displacement, const double* velocity,
                     double* force, double dt);
    
private:
    bool gpu_initialized;
    int n_cells_gpu;
    
#ifdef USE_CUDA
    // GPU arrays
    double* d_displacement;
    double* d_velocity;
    double* d_strain;
    double* d_stress;
    double* d_stress_div;
    double* d_force;
    double* d_youngs_modulus;
    double* d_poisson_ratio;
    double* d_density;
    int* d_connectivity;
#endif
};

// GPU-accelerated poroelastodynamics kernel
class PoroelastodynamicsKernelGPU : public PoroelastodynamicsKernel {
public:
    PoroelastodynamicsKernelGPU();
    ~PoroelastodynamicsKernelGPU() override;
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    // GPU-specific methods
    void allocateGPUMemory(int n_cells);
    void freeGPUMemory();
    void computeOnGPU(const double* displacement, const double* velocity,
                     const double* pressure, double* force, double* pressure_rate,
                     double dt);
    
    // Dynamic permeability update on GPU
    void updatePermeabilityGPU(double dt);
    
private:
    bool gpu_initialized;
    int n_cells_gpu;
    
#ifdef USE_CUDA
    // GPU arrays for solid mechanics
    double* d_displacement;
    double* d_velocity;
    double* d_strain;
    double* d_stress;
    double* d_force_solid;
    
    // GPU arrays for fluid flow
    double* d_pressure;
    double* d_pressure_rate;
    double* d_permeability;
    double* d_permeability_initial;
    
    // Material properties on GPU
    double* d_youngs_modulus;
    double* d_poisson_ratio;
    double* d_biot_coefficient;
    double* d_biot_modulus;
    double* d_porosity;
    double* d_density_solid;
    double* d_density_fluid;
    double* d_viscosity;
    
    // Connectivity
    int* d_connectivity;
#endif
};

// Factory function to create GPU-accelerated kernels
std::shared_ptr<PhysicsKernel> createGPUKernel(PhysicsType type);

// Helper class for managing GPU-CPU data transfers
class GPUDataTransfer {
public:
    GPUDataTransfer(size_t size);
    ~GPUDataTransfer();
    
    void copyToGPU(const void* host_data);
    void copyFromGPU(void* host_data);
    
    void* getDevicePointer() { return device_ptr; }
    size_t getSize() const { return size_bytes; }
    
private:
    void* device_ptr;
    size_t size_bytes;
};

// Benchmark GPU vs CPU performance
class GPUBenchmark {
public:
    static void runSinglePhaseFlowBenchmark(int n_cells, int n_iterations);
    static void runElastodynamicsBenchmark(int n_cells, int n_iterations);
    static void runPoroelastodynamicsBenchmark(int n_cells, int n_iterations);
    static void printBenchmarkResults();
    
private:
    static double cpu_time_ms;
    static double gpu_time_ms;
    static double speedup_factor;
};

} // namespace FSRM

#endif // PHYSICS_KERNEL_GPU_HPP
