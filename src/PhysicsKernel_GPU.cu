#ifdef USE_CUDA

#include "PhysicsKernel_GPU.hpp"
#include "GPUKernels.cuh"
#include <iostream>

namespace ResSim {

// ============================================================================
// SinglePhaseFlowKernelGPU Implementation
// ============================================================================

SinglePhaseFlowKernelGPU::SinglePhaseFlowKernelGPU() 
    : SinglePhaseFlowKernel(), gpu_initialized(false), n_cells_gpu(0),
      d_pressure(nullptr), d_pressure_old(nullptr), d_porosity(nullptr),
      d_permeability(nullptr), d_compressibility(nullptr), d_viscosity(nullptr),
      d_accumulation(nullptr), d_flux(nullptr) {
}

SinglePhaseFlowKernelGPU::~SinglePhaseFlowKernelGPU() {
    freeGPUMemory();
}

PetscErrorCode SinglePhaseFlowKernelGPU::setup(DM dm, PetscFE fe) {
    PetscFunctionBeginUser;
    
    // Call base class setup
    PetscErrorCode ierr = SinglePhaseFlowKernel::setup(dm, fe);
    CHKERRQ(ierr);
    
    // Initialize GPU
    GPUManager& gpu = GPUManager::getInstance();
    if (!gpu.isAvailable()) {
        PetscPrintf(PETSC_COMM_WORLD, "Warning: GPU not available, falling back to CPU\n");
        PetscFunctionReturn(0);
    }
    
    gpu_initialized = true;
    PetscPrintf(PETSC_COMM_WORLD, "GPU-accelerated single-phase flow kernel initialized\n");
    
    PetscFunctionReturn(0);
}

void SinglePhaseFlowKernelGPU::allocateGPUMemory(int n_cells) {
    if (!gpu_initialized) return;
    
    GPUManager& gpu = GPUManager::getInstance();
    n_cells_gpu = n_cells;
    
    d_pressure = static_cast<double*>(gpu.allocateDevice(n_cells * sizeof(double)));
    d_pressure_old = static_cast<double*>(gpu.allocateDevice(n_cells * sizeof(double)));
    d_porosity = static_cast<double*>(gpu.allocateDevice(n_cells * sizeof(double)));
    d_permeability = static_cast<double*>(gpu.allocateDevice(n_cells * sizeof(double)));
    d_compressibility = static_cast<double*>(gpu.allocateDevice(n_cells * sizeof(double)));
    d_viscosity = static_cast<double*>(gpu.allocateDevice(n_cells * sizeof(double)));
    d_accumulation = static_cast<double*>(gpu.allocateDevice(n_cells * sizeof(double)));
}

void SinglePhaseFlowKernelGPU::freeGPUMemory() {
    if (!gpu_initialized) return;
    
    GPUManager& gpu = GPUManager::getInstance();
    
    if (d_pressure) gpu.freeDevice(d_pressure);
    if (d_pressure_old) gpu.freeDevice(d_pressure_old);
    if (d_porosity) gpu.freeDevice(d_porosity);
    if (d_permeability) gpu.freeDevice(d_permeability);
    if (d_compressibility) gpu.freeDevice(d_compressibility);
    if (d_viscosity) gpu.freeDevice(d_viscosity);
    if (d_accumulation) gpu.freeDevice(d_accumulation);
    
    d_pressure = d_pressure_old = d_porosity = nullptr;
    d_permeability = d_compressibility = d_viscosity = d_accumulation = nullptr;
}

void SinglePhaseFlowKernelGPU::residual(const PetscScalar u[], const PetscScalar u_t[],
                                        const PetscScalar u_x[], const PetscScalar a[],
                                        const PetscReal x[], PetscScalar f[]) {
    if (!gpu_initialized) {
        // Fallback to CPU
        SinglePhaseFlowKernel::residual(u, u_t, u_x, a, x, f);
        return;
    }
    
    // GPU execution
    // This is a simplified example - full implementation would handle
    // data transfers and kernel launches
    f[0] = porosity * compressibility * u_t[0];
}

void SinglePhaseFlowKernelGPU::jacobian(const PetscScalar u[], const PetscScalar u_t[],
                                        const PetscScalar u_x[], const PetscScalar a[],
                                        const PetscReal x[], PetscScalar J[]) {
    // Always use full Jacobian from base class - GPU acceleration applies to
    // matrix assembly, not to the pointwise Jacobian evaluation
    SinglePhaseFlowKernel::jacobian(u, u_t, u_x, a, x, J);
}

// ============================================================================
// ElastodynamicsKernelGPU Implementation
// ============================================================================

ElastodynamicsKernelGPU::ElastodynamicsKernelGPU()
    : ElastodynamicsKernel(), gpu_initialized(false), n_cells_gpu(0),
      d_displacement(nullptr), d_velocity(nullptr), d_strain(nullptr),
      d_stress(nullptr), d_stress_div(nullptr), d_force(nullptr),
      d_youngs_modulus(nullptr), d_poisson_ratio(nullptr), d_density(nullptr),
      d_connectivity(nullptr) {
}

ElastodynamicsKernelGPU::~ElastodynamicsKernelGPU() {
    freeGPUMemory();
}

PetscErrorCode ElastodynamicsKernelGPU::setup(DM dm, PetscFE fe) {
    PetscFunctionBeginUser;
    
    PetscErrorCode ierr = ElastodynamicsKernel::setup(dm, fe);
    CHKERRQ(ierr);
    
    GPUManager& gpu = GPUManager::getInstance();
    if (!gpu.isAvailable()) {
        PetscPrintf(PETSC_COMM_WORLD, "Warning: GPU not available for elastodynamics\n");
        PetscFunctionReturn(0);
    }
    
    gpu_initialized = true;
    PetscPrintf(PETSC_COMM_WORLD, "GPU-accelerated elastodynamics kernel initialized\n");
    
    PetscFunctionReturn(0);
}

void ElastodynamicsKernelGPU::allocateGPUMemory(int n_cells) {
    if (!gpu_initialized) return;
    
    GPUManager& gpu = GPUManager::getInstance();
    n_cells_gpu = n_cells;
    
    int dim = 3;
    d_displacement = static_cast<double*>(gpu.allocateDevice(n_cells * dim * sizeof(double)));
    d_velocity = static_cast<double*>(gpu.allocateDevice(n_cells * dim * sizeof(double)));
    d_strain = static_cast<double*>(gpu.allocateDevice(n_cells * 6 * sizeof(double)));
    d_stress = static_cast<double*>(gpu.allocateDevice(n_cells * 6 * sizeof(double)));
    d_stress_div = static_cast<double*>(gpu.allocateDevice(n_cells * dim * sizeof(double)));
    d_force = static_cast<double*>(gpu.allocateDevice(n_cells * dim * sizeof(double)));
    d_youngs_modulus = static_cast<double*>(gpu.allocateDevice(n_cells * sizeof(double)));
    d_poisson_ratio = static_cast<double*>(gpu.allocateDevice(n_cells * sizeof(double)));
    d_density = static_cast<double*>(gpu.allocateDevice(n_cells * sizeof(double)));
}

void ElastodynamicsKernelGPU::freeGPUMemory() {
    if (!gpu_initialized) return;
    
    GPUManager& gpu = GPUManager::getInstance();
    
    if (d_displacement) gpu.freeDevice(d_displacement);
    if (d_velocity) gpu.freeDevice(d_velocity);
    if (d_strain) gpu.freeDevice(d_strain);
    if (d_stress) gpu.freeDevice(d_stress);
    if (d_stress_div) gpu.freeDevice(d_stress_div);
    if (d_force) gpu.freeDevice(d_force);
    if (d_youngs_modulus) gpu.freeDevice(d_youngs_modulus);
    if (d_poisson_ratio) gpu.freeDevice(d_poisson_ratio);
    if (d_density) gpu.freeDevice(d_density);
    if (d_connectivity) gpu.freeDevice(d_connectivity);
}

void ElastodynamicsKernelGPU::computeOnGPU(const double* displacement, 
                                           const double* velocity,
                                           double* force, double dt) {
    if (!gpu_initialized) return;
    
    GPUManager& gpu = GPUManager::getInstance();
    int dim = 3;
    
    // Copy data to GPU
    gpu.copyHostToDevice(d_displacement, displacement, n_cells_gpu * dim * sizeof(double));
    gpu.copyHostToDevice(d_velocity, velocity, n_cells_gpu * dim * sizeof(double));
    
    // Launch GPU kernels
    int block_size = GPU::BLOCK_SIZE_1D;
    int grid_size = (n_cells_gpu + block_size - 1) / block_size;
    
    // Compute strain
    GPU::computeStrain<<<grid_size, block_size>>>(
        d_displacement, d_strain, d_connectivity, n_cells_gpu, dim
    );
    
    // Compute stress
    GPU::computeStress<<<grid_size, block_size>>>(
        d_strain, d_youngs_modulus, d_poisson_ratio, d_stress, n_cells_gpu, dim
    );
    
    // Compute stress divergence
    GPU::computeStressDivergence<<<grid_size, block_size>>>(
        d_stress, d_connectivity, d_stress_div, n_cells_gpu, dim
    );
    
    // Update velocity
    GPU::updateVelocity<<<grid_size, block_size>>>(
        d_velocity, d_stress_div, d_density, dt, n_cells_gpu, dim
    );
    
    // Apply damping
    GPU::applyRayleighDamping<<<grid_size, block_size>>>(
        d_force, d_velocity, d_displacement, d_density, d_youngs_modulus,
        damping_alpha, damping_beta, n_cells_gpu, dim
    );
    
    gpu.synchronize();
    
    // Copy results back
    gpu.copyDeviceToHost(force, d_force, n_cells_gpu * dim * sizeof(double));
}

void ElastodynamicsKernelGPU::residual(const PetscScalar u[], const PetscScalar u_t[],
                                       const PetscScalar u_x[], const PetscScalar a[],
                                       const PetscReal x[], PetscScalar f[]) {
    if (!gpu_initialized) {
        ElastodynamicsKernel::residual(u, u_t, u_x, a, x, f);
        return;
    }
    
    // GPU execution - simplified
    f[0] = density * u_t[0];
    f[1] = density * u_t[1];
    f[2] = density * u_t[2];
    
    f[0] += damping_alpha * density * u_t[0];
    f[1] += damping_alpha * density * u_t[1];
    f[2] += damping_alpha * density * u_t[2];
}

void ElastodynamicsKernelGPU::jacobian(const PetscScalar u[], const PetscScalar u_t[],
                                       const PetscScalar u_x[], const PetscScalar a[],
                                       const PetscReal x[], PetscScalar J[]) {
    // Always use full Jacobian from base class - includes mass matrix + full
    // stiffness tensor with Rayleigh damping (α*M + β*K)
    ElastodynamicsKernel::jacobian(u, u_t, u_x, a, x, J);
}

// ============================================================================
// PoroelastodynamicsKernelGPU Implementation
// ============================================================================

PoroelastodynamicsKernelGPU::PoroelastodynamicsKernelGPU()
    : PoroelastodynamicsKernel(), gpu_initialized(false), n_cells_gpu(0) {
    // Initialize all GPU pointers to nullptr
    d_displacement = d_velocity = d_strain = d_stress = d_force_solid = nullptr;
    d_pressure = d_pressure_rate = d_permeability = d_permeability_initial = nullptr;
    d_youngs_modulus = d_poisson_ratio = d_biot_coefficient = d_biot_modulus = nullptr;
    d_porosity = d_density_solid = d_density_fluid = d_viscosity = nullptr;
    d_connectivity = nullptr;
}

PoroelastodynamicsKernelGPU::~PoroelastodynamicsKernelGPU() {
    freeGPUMemory();
}

PetscErrorCode PoroelastodynamicsKernelGPU::setup(DM dm, PetscFE fe) {
    PetscFunctionBeginUser;
    
    PetscErrorCode ierr = PoroelastodynamicsKernel::setup(dm, fe);
    CHKERRQ(ierr);
    
    GPUManager& gpu = GPUManager::getInstance();
    if (!gpu.isAvailable()) {
        PetscPrintf(PETSC_COMM_WORLD, "Warning: GPU not available for poroelastodynamics\n");
        PetscFunctionReturn(0);
    }
    
    gpu_initialized = true;
    PetscPrintf(PETSC_COMM_WORLD, "GPU-accelerated poroelastodynamics kernel initialized\n");
    
    PetscFunctionReturn(0);
}

void PoroelastodynamicsKernelGPU::allocateGPUMemory(int n_cells) {
    if (!gpu_initialized) return;
    
    GPUManager& gpu = GPUManager::getInstance();
    n_cells_gpu = n_cells;
    
    int dim = 3;
    d_displacement = static_cast<double*>(gpu.allocateDevice(n_cells * dim * sizeof(double)));
    d_velocity = static_cast<double*>(gpu.allocateDevice(n_cells * dim * sizeof(double)));
    d_strain = static_cast<double*>(gpu.allocateDevice(n_cells * 6 * sizeof(double)));
    d_stress = static_cast<double*>(gpu.allocateDevice(n_cells * 6 * sizeof(double)));
    d_force_solid = static_cast<double*>(gpu.allocateDevice(n_cells * dim * sizeof(double)));
    d_pressure = static_cast<double*>(gpu.allocateDevice(n_cells * sizeof(double)));
    d_pressure_rate = static_cast<double*>(gpu.allocateDevice(n_cells * sizeof(double)));
    d_permeability = static_cast<double*>(gpu.allocateDevice(n_cells * sizeof(double)));
    d_permeability_initial = static_cast<double*>(gpu.allocateDevice(n_cells * sizeof(double)));
    d_youngs_modulus = static_cast<double*>(gpu.allocateDevice(n_cells * sizeof(double)));
    d_poisson_ratio = static_cast<double*>(gpu.allocateDevice(n_cells * sizeof(double)));
    d_biot_coefficient = static_cast<double*>(gpu.allocateDevice(n_cells * sizeof(double)));
    d_biot_modulus = static_cast<double*>(gpu.allocateDevice(n_cells * sizeof(double)));
    d_porosity = static_cast<double*>(gpu.allocateDevice(n_cells * sizeof(double)));
    d_density_solid = static_cast<double*>(gpu.allocateDevice(n_cells * sizeof(double)));
    d_density_fluid = static_cast<double*>(gpu.allocateDevice(n_cells * sizeof(double)));
    d_viscosity = static_cast<double*>(gpu.allocateDevice(n_cells * sizeof(double)));
}

void PoroelastodynamicsKernelGPU::freeGPUMemory() {
    if (!gpu_initialized) return;
    
    GPUManager& gpu = GPUManager::getInstance();
    
    // Free all allocated GPU memory
    if (d_displacement) gpu.freeDevice(d_displacement);
    if (d_velocity) gpu.freeDevice(d_velocity);
    if (d_strain) gpu.freeDevice(d_strain);
    if (d_stress) gpu.freeDevice(d_stress);
    if (d_force_solid) gpu.freeDevice(d_force_solid);
    if (d_pressure) gpu.freeDevice(d_pressure);
    if (d_pressure_rate) gpu.freeDevice(d_pressure_rate);
    if (d_permeability) gpu.freeDevice(d_permeability);
    if (d_permeability_initial) gpu.freeDevice(d_permeability_initial);
    if (d_youngs_modulus) gpu.freeDevice(d_youngs_modulus);
    if (d_poisson_ratio) gpu.freeDevice(d_poisson_ratio);
    if (d_biot_coefficient) gpu.freeDevice(d_biot_coefficient);
    if (d_biot_modulus) gpu.freeDevice(d_biot_modulus);
    if (d_porosity) gpu.freeDevice(d_porosity);
    if (d_density_solid) gpu.freeDevice(d_density_solid);
    if (d_density_fluid) gpu.freeDevice(d_density_fluid);
    if (d_viscosity) gpu.freeDevice(d_viscosity);
    if (d_connectivity) gpu.freeDevice(d_connectivity);
}

void PoroelastodynamicsKernelGPU::updatePermeabilityGPU(double dt) {
    if (!gpu_initialized || !enable_dynamic_k) return;
    
    GPUManager& gpu = GPUManager::getInstance();
    
    int block_size = GPU::BLOCK_SIZE_1D;
    int grid_size = (n_cells_gpu + block_size - 1) / block_size;
    
    GPU::updateDynamicPermeability<<<grid_size, block_size>>>(
        d_permeability, d_permeability_initial, d_strain, d_stress,
        k_strain_coeff, k_stress_coeff, 1.0 / k_recovery_time, dt,
        1e-18, 1e-10, n_cells_gpu
    );
    
    gpu.synchronize();
}

void PoroelastodynamicsKernelGPU::residual(const PetscScalar u[], const PetscScalar u_t[],
                                           const PetscScalar u_x[], const PetscScalar a[],
                                           const PetscReal x[], PetscScalar f[]) {
    if (!gpu_initialized) {
        PoroelastodynamicsKernel::residual(u, u_t, u_x, a, x, f);
        return;
    }
    
    // Simplified GPU residual evaluation
    double rho_bulk = (1.0 - porosity) * density_solid + porosity * density_fluid;
    
    f[0] = rho_bulk * u_t[0];
    f[1] = rho_bulk * u_t[1];
    f[2] = rho_bulk * u_t[2];
    f[3] = biot_coefficient * (u_t[0] + u_t[1] + u_t[2]) + u_t[3] / biot_modulus;
    
    f[0] += damping_alpha * rho_bulk * u_t[0];
    f[1] += damping_alpha * rho_bulk * u_t[1];
    f[2] += damping_alpha * rho_bulk * u_t[2];
}

void PoroelastodynamicsKernelGPU::jacobian(const PetscScalar u[], const PetscScalar u_t[],
                                           const PetscScalar u_x[], const PetscScalar a[],
                                           const PetscReal x[], PetscScalar J[]) {
    // Always use full Jacobian from base class - includes:
    // - Mass matrix with damping
    // - Full elasticity tensor (undrained moduli for fast waves)
    // - Biot coupling terms (α for u-p and p-u coupling)
    // - Darcy flux Jacobian for pressure gradient
    PoroelastodynamicsKernel::jacobian(u, u_t, u_x, a, x, J);
}

// ============================================================================
// Factory function
// ============================================================================

std::shared_ptr<PhysicsKernel> createGPUKernel(PhysicsType type) {
    GPUManager& gpu = GPUManager::getInstance();
    
    if (!gpu.isAvailable()) {
        return nullptr;
    }
    
    switch (type) {
        case PhysicsType::FLUID_FLOW:
            return std::make_shared<SinglePhaseFlowKernelGPU>();
        case PhysicsType::ELASTODYNAMICS:
            return std::make_shared<ElastodynamicsKernelGPU>();
        case PhysicsType::POROELASTODYNAMICS:
            return std::make_shared<PoroelastodynamicsKernelGPU>();
        default:
            return nullptr;
    }
}

} // namespace ResSim

#endif // USE_CUDA
