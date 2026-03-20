#ifndef GPU_KERNELS_CUH
#define GPU_KERNELS_CUH

#ifdef USE_CUDA

#include <cuda_runtime.h>
#include <cuda.h>

namespace FSRM {
namespace GPU {

// Block and grid dimensions
constexpr int BLOCK_SIZE_1D = 256;
constexpr int BLOCK_SIZE_2D = 16;
constexpr int BLOCK_SIZE_3D = 8;

// ============================================================================
// Single-Phase Flow GPU Kernels
// ============================================================================

// Compute accumulation term: phi * ct * dP/dt
__global__ void computeAccumulation(
    const double* pressure_old,
    const double* pressure_new,
    const double* porosity,
    const double* compressibility,
    double dt,
    double* accumulation,
    int n_cells
);

// Compute Darcy flux: -k/mu * grad(P)
__global__ void computeDarcyFlux(
    const double* pressure,
    const double* permeability,
    const double* viscosity,
    const int* connectivity,
    double* flux,
    int n_faces,
    int dim
);

// Update pressure field
__global__ void updatePressure(
    double* pressure,
    const double* accumulation,
    const double* flux_divergence,
    double dt,
    int n_cells
);

// ============================================================================
// Elastodynamics GPU Kernels
// ============================================================================

// Compute strain tensor from displacement gradients
__global__ void computeStrain(
    const double* displacement,
    double* strain,
    const int* connectivity,
    int n_cells,
    int dim
);

// Compute stress tensor from strain (linear elasticity)
__global__ void computeStress(
    const double* strain,
    const double* youngs_modulus,
    const double* poisson_ratio,
    double* stress,
    int n_cells,
    int dim
);

// Compute stress divergence for elastodynamics
__global__ void computeStressDivergence(
    const double* stress,
    const int* connectivity,
    double* stress_div,
    int n_cells,
    int dim
);

// Update velocity: v_new = v_old + dt * (div(sigma) / rho)
__global__ void updateVelocity(
    double* velocity,
    const double* stress_div,
    const double* density,
    double dt,
    int n_cells,
    int dim
);

// Update displacement: u_new = u_old + dt * v
__global__ void updateDisplacement(
    double* displacement,
    const double* velocity,
    double dt,
    int n_cells,
    int dim
);

// Add Rayleigh damping: F_damp = alpha*M*v + beta*K*u
__global__ void applyRayleighDamping(
    double* force,
    const double* velocity,
    const double* displacement,
    const double* density,
    const double* stiffness,
    double alpha,
    double beta,
    int n_cells,
    int dim
);

// ============================================================================
// Poroelastodynamics GPU Kernels
// ============================================================================

// Coupled solid-fluid momentum equation
__global__ void computePoroelasticForce(
    const double* displacement,
    const double* pressure,
    const double* velocity_solid,
    const double* youngs_modulus,
    const double* poisson_ratio,
    const double* biot_coefficient,
    const double* porosity,
    const double* density_solid,
    const double* density_fluid,
    double* force_solid,
    const int* connectivity,
    int n_cells,
    int dim
);

// Fluid mass conservation with Darcy flow
__global__ void computeFluidMassBalance(
    const double* pressure,
    const double* displacement,
    const double* permeability,
    const double* viscosity,
    const double* porosity,
    const double* biot_coefficient,
    const double* biot_modulus,
    double dt,
    double* pressure_rate,
    const int* connectivity,
    int n_cells
);

// Dynamic permeability update based on strain/stress
__global__ void updateDynamicPermeability(
    double* permeability,
    const double* permeability_initial,
    const double* strain,
    const double* stress,
    double strain_coefficient,
    double stress_coefficient,
    double recovery_rate,
    double dt,
    double k_min,
    double k_max,
    int n_cells
);

// Check for static triggering threshold
__global__ void checkStaticTrigger(
    const double* stress,
    double threshold,
    bool* triggered,
    int* trigger_location,
    int n_cells
);

// ============================================================================
// Black Oil GPU Kernels
// ============================================================================

// Compute oil, water, gas phase saturations
__global__ void updateSaturations(
    double* saturation_oil,
    double* saturation_water,
    double* saturation_gas,
    const double* pressure,
    const double* composition,
    int n_cells
);

// Compute phase densities and viscosities (PVT)
__global__ void computePVT(
    const double* pressure,
    const double* temperature,
    double* density_oil,
    double* density_water,
    double* density_gas,
    double* viscosity_oil,
    double* viscosity_water,
    double* viscosity_gas,
    int n_cells
);

// Compute relative permeabilities
__global__ void computeRelativePermeability(
    const double* saturation_oil,
    const double* saturation_water,
    const double* saturation_gas,
    double* kr_oil,
    double* kr_water,
    double* kr_gas,
    int n_cells
);

// ============================================================================
// Thermal Transport GPU Kernels
// ============================================================================

// Heat conduction: div(k * grad(T))
__global__ void computeHeatFlux(
    const double* temperature,
    const double* thermal_conductivity,
    double* heat_flux,
    const int* connectivity,
    int n_cells,
    int dim
);

// Temperature update: rho*cp*dT/dt = div(k*grad(T)) + Q
__global__ void updateTemperature(
    double* temperature,
    const double* heat_flux_div,
    const double* heat_capacity,
    const double* density,
    const double* heat_source,
    double dt,
    int n_cells
);

// ============================================================================
// Utility Kernels
// ============================================================================

// Vector operations
__global__ void vectorAdd(const double* a, const double* b, double* c, int n);
__global__ void vectorScale(const double* a, double scale, double* b, int n);
__global__ void vectorDot(const double* a, const double* b, double* result, int n);
__global__ void vectorNorm(const double* a, double* result, int n);

// Matrix-vector multiply for sparse matrices
__global__ void sparseMatVec(
    const double* values,
    const int* row_ptr,
    const int* col_idx,
    const double* x,
    double* y,
    int n_rows
);

// Reduction operations
template<typename T>
__global__ void reduceSum(const T* input, T* output, int n);

template<typename T>
__global__ void reduceMax(const T* input, T* output, int n);

template<typename T>
__global__ void reduceMin(const T* input, T* output, int n);

// ============================================================================
// Wave Propagation Kernels
// ============================================================================

// Compute P-wave and S-wave velocities from elastic constants
__global__ void computeWaveVelocities(
    const double* youngs_modulus,
    const double* poisson_ratio,
    const double* density,
    double* vp,
    double* vs,
    int n_cells
);

// Apply absorbing boundary conditions (Perfectly Matched Layer)
__global__ void applyPML(
    double* displacement,
    double* velocity,
    const double* damping_profile,
    int boundary_width,
    int n_cells,
    int dim
);

// Compute seismic moment and energy release
__global__ void computeSeismicMoment(
    const double* displacement,
    const double* stress,
    const double* shear_modulus,
    const double* fault_area,
    double* moment,
    double* energy,
    int n_fault_cells
);

// ============================================================================
// Helper Functions (Device)
// ============================================================================

__device__ inline int getGlobalThreadId1D() {
    return blockIdx.x * blockDim.x + threadIdx.x;
}

__device__ inline int getGlobalThreadId2D(int& i, int& j) {
    i = blockIdx.x * blockDim.x + threadIdx.x;
    j = blockIdx.y * blockDim.y + threadIdx.y;
    return i + j * gridDim.x * blockDim.x;
}

__device__ inline int getGlobalThreadId3D(int& i, int& j, int& k) {
    i = blockIdx.x * blockDim.x + threadIdx.x;
    j = blockIdx.y * blockDim.y + threadIdx.y;
    k = blockIdx.z * blockDim.z + threadIdx.z;
    return i + j * gridDim.x * blockDim.x + k * gridDim.x * blockDim.x * gridDim.y * blockDim.y;
}

// Von Mises stress
__device__ inline double computeVonMisesStress(
    double s_xx, double s_yy, double s_zz,
    double s_xy, double s_xz, double s_yz
) {
    double s_mean = (s_xx + s_yy + s_zz) / 3.0;
    double dev_xx = s_xx - s_mean;
    double dev_yy = s_yy - s_mean;
    double dev_zz = s_zz - s_mean;
    
    return sqrt(1.5 * (dev_xx*dev_xx + dev_yy*dev_yy + dev_zz*dev_zz +
                       2.0 * (s_xy*s_xy + s_xz*s_xz + s_yz*s_yz)));
}

// Lame parameters
__device__ inline void computeLameParameters(
    double E, double nu,
    double& lambda, double& mu
) {
    lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    mu = E / (2.0 * (1.0 + nu));
}

} // namespace GPU
} // namespace FSRM

#endif // USE_CUDA

#endif // GPU_KERNELS_CUH
