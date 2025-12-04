# FSRM Physics and GPU Architecture

This document provides a unified reference for FSRM's physics kernel system and GPU acceleration architecture. All physics kernels are designed to run on both CPU and GPU through a unified execution policy system.

---

## Table of Contents

1. [Architecture Overview](#architecture-overview)
2. [Physics Kernel Design](#physics-kernel-design)
3. [Execution Policies](#execution-policies)
4. [GPU Acceleration](#gpu-acceleration)
5. [Kernel Reference](#kernel-reference)
6. [Configuration](#configuration)
7. [Performance Guidelines](#performance-guidelines)

---

## Architecture Overview

FSRM uses a unified physics kernel architecture where:

1. **All physics equations** are implemented as kernel classes derived from `PhysicsKernel`
2. **Kernels are execution-agnostic** - the same physics logic runs on CPU or GPU
3. **GPU versions inherit from CPU base classes** - ensures consistent behavior
4. **Automatic fallback** - GPU kernels fall back to CPU if GPU unavailable

### Class Hierarchy

```
PhysicsKernel (abstract base)
├── SinglePhaseFlowKernel
│   └── SinglePhaseFlowKernelGPU
├── BlackOilKernel
│   └── BlackOilKernelGPU (planned)
├── CompositionalKernel
│   └── CompositionalKernelGPU (planned)
├── GeomechanicsKernel
│   └── GeomechanicsKernelGPU (planned)
├── ThermalKernel
│   └── ThermalKernelGPU (planned)
├── ElastodynamicsKernel
│   └── ElastodynamicsKernelGPU
├── PoroelastodynamicsKernel
│   └── PoroelastodynamicsKernelGPU
├── ParticleTransportKernel
├── FracturePropagationKernel
└── TidalForcesKernel
```

---

## Physics Kernel Design

### Base Kernel Interface

All physics kernels implement the same interface for PETSc finite element integration:

```cpp
class PhysicsKernel {
public:
    // Core methods required for PETSc FE integration
    virtual PetscErrorCode setup(DM dm, PetscFE fe) = 0;
    
    virtual void residual(const PetscScalar u[], const PetscScalar u_t[],
                         const PetscScalar u_x[], const PetscScalar a[],
                         const PetscReal x[], PetscScalar f[]) = 0;
    
    virtual void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                         const PetscScalar u_x[], const PetscScalar a[],
                         const PetscReal x[], PetscScalar J[]) = 0;
    
    virtual int getNumFields() const = 0;
    virtual int getNumComponents(int field) const = 0;
};
```

### Parameter Descriptions

| Parameter | Description |
|-----------|-------------|
| `u[]` | Solution values at quadrature point |
| `u_t[]` | Time derivatives of solution |
| `u_x[]` | Spatial gradients of solution |
| `a[]` | Auxiliary data (shift parameters from time integrator) |
| `x[]` | Physical coordinates of quadrature point |
| `f[]` | Output: residual vector contribution |
| `J[]` | Output: Jacobian matrix contribution |

### Jacobian Structure

The Jacobian combines accumulation (g0) and flux (g3) terms:

```
J = g0 + g3

g0: ∂f/∂u + shift * ∂f/∂(∂u/∂t)  (solution and time derivative)
g3: ∂(flux)/∂(∇u)                 (gradient terms, e.g., diffusion)
```

For coupled multi-field problems, the Jacobian is a block matrix:

```
J = [ J_11  J_12  ...  J_1n ]
    [ J_21  J_22  ...  J_2n ]
    [ ...   ...   ...  ...  ]
    [ J_n1  J_n2  ...  J_nn ]
```

---

## Execution Policies

### GPU Execution Modes

FSRM supports four GPU execution modes controlled by `gpu_mode` configuration:

| Mode | Description | Use Case |
|------|-------------|----------|
| `CPU_ONLY` | Run entirely on CPU | Testing, small problems |
| `GPU_ONLY` | Run entirely on GPU | Large problems with GPU |
| `HYBRID` | Automatic CPU/GPU load balancing | Multi-GPU + CPU systems |
| `CPU_FALLBACK` | Try GPU, fallback to CPU | Production (default) |

### Execution Selection Logic

```cpp
// Factory function for creating execution-appropriate kernels
std::shared_ptr<PhysicsKernel> createKernel(PhysicsType type, 
                                            GPUExecutionMode mode) {
    if (mode != GPUExecutionMode::CPU_ONLY) {
        auto gpu_kernel = createGPUKernel(type);
        if (gpu_kernel) return gpu_kernel;
    }
    // Fall back to CPU kernel
    return createCPUKernel(type);
}
```

### Automatic GPU Detection

```cpp
GPUManager& gpu = GPUManager::getInstance();
if (gpu.isAvailable()) {
    // GPU kernels will be used
} else {
    // CPU fallback automatically engaged
}
```

---

## GPU Acceleration

### GPU Kernel Structure

GPU kernels inherit from CPU base classes and override execution methods:

```cpp
class SinglePhaseFlowKernelGPU : public SinglePhaseFlowKernel {
public:
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    // GPU-specific methods
    void allocateGPUMemory(int n_cells);
    void freeGPUMemory();
    void copyDataToGPU(const double* host_data, size_t size);
    void copyDataFromGPU(double* host_data, size_t size);
    
private:
    bool gpu_initialized;
    // Device pointers (CUDA/HIP)
#ifdef USE_CUDA
    double* d_pressure;
    double* d_porosity;
    // ... additional device arrays
#endif
};
```

### GPU Memory Management

FSRM uses RAII for GPU memory via `GPUArray<T>`:

```cpp
// Automatic memory management
GPUArray<double> pressure(n_cells);
pressure.copyFromHost(host_pressure);
// Use d_pressure = pressure.data() in kernels
pressure.copyToHost(host_pressure);
// Memory freed automatically on destruction
```

### Available GPU Kernels

| Kernel | CUDA | HIP | Description |
|--------|------|-----|-------------|
| Single-Phase Flow | ✓ | ✓ | Darcy flow accumulation, flux |
| Elastodynamics | ✓ | ✓ | Wave propagation, stress/strain |
| Poroelastodynamics | ✓ | ✓ | Coupled fluid-solid waves |
| Black Oil | Planned | Planned | Three-phase flow |
| Thermal | Planned | Planned | Heat conduction/convection |

### CUDA Kernel Execution

```cpp
// Block and grid configuration
constexpr int BLOCK_SIZE_1D = 256;
int grid_size = (n_cells + BLOCK_SIZE_1D - 1) / BLOCK_SIZE_1D;

// Launch kernel
computeStrain<<<grid_size, BLOCK_SIZE_1D>>>(
    d_displacement, d_strain, d_connectivity, n_cells, dim
);

// Synchronize
gpu.synchronize();
```

---

## Kernel Reference

### Flow Kernels

#### SinglePhaseFlowKernel

**Governing Equation:**
$$\phi c_t \frac{\partial p}{\partial t} = \nabla \cdot \left(\frac{k}{\mu} \nabla p\right) + q$$

**Fields:** 1 (pressure)

**Properties:**
- `porosity` (φ): Volume fraction of pore space
- `permeability` (k): Rock permeability [m²]
- `compressibility` (c_t): Total compressibility [1/Pa]
- `viscosity` (μ): Fluid dynamic viscosity [Pa·s]

**GPU Acceleration:** Full support via `SinglePhaseFlowKernelGPU`

#### BlackOilKernel

**Governing Equations:** Mass conservation for oil, water, gas phases with dissolution

**Fields:** 3 (pressure, water saturation, gas saturation)

**Key Features:**
- PVT correlations (Standing, Vazquez-Beggs)
- Relative permeability (Corey model)
- Capillary pressure
- Solution gas/oil ratio

**GPU Acceleration:** Planned

#### CompositionalKernel

**Governing Equations:** Component mass balance with flash calculations

**Fields:** N+1 (pressure + N-1 compositions)

**Key Features:**
- Peng-Robinson equation of state
- Rachford-Rice flash calculation
- Fugacity coefficient computation

**GPU Acceleration:** Planned

---

### Mechanics Kernels

#### GeomechanicsKernel

**Governing Equation:**
$$\nabla \cdot \boldsymbol{\sigma} + \rho \mathbf{g} = \rho \frac{\partial^2 \mathbf{u}}{\partial t^2}$$

**Fields:** 1 (displacement vector: 3 components)

**Properties:**
- `youngs_modulus` (E): Elastic modulus [Pa]
- `poisson_ratio` (ν): Poisson's ratio [-]
- `density` (ρ): Bulk density [kg/m³]

**Stress-Strain Relation (Isotropic):**
$$\sigma_{ij} = \lambda \epsilon_{kk} \delta_{ij} + 2G \epsilon_{ij}$$

where λ and G are Lamé parameters.

**GPU Acceleration:** Planned

#### ElastodynamicsKernel

**Governing Equation:**
$$\rho \frac{\partial^2 \mathbf{u}}{\partial t^2} = \nabla \cdot \boldsymbol{\sigma} + \mathbf{f} - \mathbf{C} \cdot \dot{\mathbf{u}}$$

**Fields:** 1 (displacement vector: 3 components)

**Additional Properties:**
- `p_wave_velocity` (V_p): P-wave velocity [m/s]
- `s_wave_velocity` (V_s): S-wave velocity [m/s]
- `damping_alpha`: Rayleigh mass damping
- `damping_beta`: Rayleigh stiffness damping
- `quality_factor` (Q): Seismic attenuation factor

**Features:**
- Full inertial dynamics
- Rayleigh damping: C = αM + βK
- Static-to-dynamic triggering
- PML absorbing boundaries

**GPU Acceleration:** Full support via `ElastodynamicsKernelGPU`

#### PoroelastodynamicsKernel

**Governing Equations (Biot Theory):**

Momentum balance:
$$\rho \frac{\partial^2 \mathbf{u}}{\partial t^2} = \nabla \cdot \boldsymbol{\sigma}' - \alpha \nabla p$$

Mass conservation:
$$\frac{1}{M} \frac{\partial p}{\partial t} + \alpha \nabla \cdot \dot{\mathbf{u}} + \nabla \cdot \mathbf{q} = 0$$

**Fields:** 2 (displacement: 3 components, pressure: 1 component)

**Additional Properties:**
- `biot_coefficient` (α): Biot coefficient
- `biot_modulus` (M): Biot modulus [Pa]
- `permeability`: Initial permeability [m²]
- `p_wave_fast`: Fast P-wave velocity
- `p_wave_slow`: Slow P-wave (Biot wave II)

**Features:**
- Coupled solid-fluid dynamics
- Three wave types: fast P, slow P, S waves
- Dynamic permeability enhancement
- Undrained/drained behavior transition

**GPU Acceleration:** Full support via `PoroelastodynamicsKernelGPU`

---

### Thermal Kernel

#### ThermalKernel

**Governing Equation:**
$$(\rho c_p)_{eff} \frac{\partial T}{\partial t} = \nabla \cdot (k_T \nabla T) + Q$$

**Fields:** 1 (temperature)

**Properties:**
- Effective heat capacity: $(\rho c_p)_{eff} = (1-\phi)\rho_s c_{p,s} + \phi \rho_f c_{p,f}$
- Effective conductivity: $k_{T,eff} = (1-\phi) k_{T,s} + \phi k_{T,f}$

**GPU Acceleration:** Planned

---

### Transport Kernels

#### ParticleTransportKernel

**Governing Equation:**
$$\phi \frac{\partial C}{\partial t} + \nabla \cdot (\mathbf{v} C) - \nabla \cdot (D \nabla C) + v_s \frac{\partial C}{\partial z} = 0$$

**Fields:** 1 (concentration)

**Properties:**
- `particle_diameter`: Particle size [m]
- `particle_density`: Particle density [kg/m³]
- `diffusivity`: Molecular diffusivity [m²/s]
- `longitudinal_dispersivity`: Mechanical dispersion [m]
- `transverse_dispersivity`: Transverse dispersion [m]

**Features:**
- Gravitational settling (Stokes velocity)
- Mechanical dispersion
- Bridging/filtration effects

**GPU Acceleration:** Planned

---

### Dynamic Permeability Model

For poroelastodynamics, permeability can change dynamically:

```cpp
class DynamicPermeabilityModel {
    // Permeability enhancement from wave passage
    double k_new = k_initial * (1 + strain_coeff * |ε| + stress_coeff * |σ|);
    
    // Time-dependent recovery
    k_new = k_initial + (k_new - k_initial) * exp(-dt / τ_recovery);
    
    // Physical bounds
    k_new = clamp(k_new, k_min, k_max);
};
```

**Physical Mechanisms:**
1. **Volumetric strain model**: k = k₀ exp(β ε_vol)
2. **Shear dilation model**: k = k₀ (1 + γ |ε_dev|)
3. **Stress model**: k = k₀ exp(-α σ')
4. **Crack opening model**: k ∝ w³ (cubic law)

---

## Configuration

### Enabling GPU Execution

```ini
[SIMULATION]
use_gpu = true
gpu_mode = CPU_FALLBACK       # AUTO, GPU_ONLY, CPU_FALLBACK, HYBRID
gpu_device_id = 0             # GPU device to use
gpu_memory_fraction = 0.8     # Fraction of GPU memory to use
gpu_verbose = false           # Print GPU information

# GPU solver options
use_gpu_preconditioner = true
use_gpu_matrix_assembly = true
pin_host_memory = true        # Faster CPU-GPU transfers
```

### Enabling Dynamic Physics

```ini
[SIMULATION]
# Wave propagation
enable_elastodynamics = true
enable_poroelastodynamics = false

# Static-to-dynamic triggering
use_static_triggering = true
dynamic_trigger_threshold = 1.0e6    # Pa (stress threshold)
dynamic_event_duration = 10.0        # seconds

# Dynamic permeability
enable_dynamic_permeability_change = true
permeability_sensitivity = 1.0
permeability_recovery_time = 100.0   # seconds
```

### Physics Kernel Selection

```ini
[SIMULATION]
fluid_model = BLACK_OIL       # SINGLE_COMPONENT, BLACK_OIL, COMPOSITIONAL
solid_model = POROELASTIC     # ELASTIC, VISCOELASTIC, POROELASTIC, ELASTODYNAMIC
```

---

## Performance Guidelines

### GPU Suitability by Problem Type

| Physics | Grid Size | GPU Benefit |
|---------|-----------|-------------|
| Single-phase flow | >100K cells | 2-5× |
| Black oil | >50K cells | 3-8× |
| Geomechanics | >100K cells | 3-10× |
| Elastodynamics | >10K cells | 5-20× |
| Poroelastodynamics | >10K cells | 5-50× |

### GPU Memory Estimation

```
Memory (bytes) ≈ n_cells × n_fields × n_components × 8 × factor

Factors:
- Single-phase: factor ≈ 10-20
- Black oil: factor ≈ 30-50
- Elastodynamics: factor ≈ 40-80
- Poroelastodynamics: factor ≈ 60-100
```

### Optimal Kernel Launch Configuration

| Problem Dimension | Block Size | Notes |
|-------------------|------------|-------|
| 1D arrays | 256 | Standard for vector ops |
| 2D grids | 16×16 | Good for spatial locality |
| 3D grids | 8×8×8 | Limited by shared memory |

### Data Transfer Optimization

1. **Minimize transfers**: Keep data on GPU across timesteps
2. **Async transfers**: Overlap computation with data movement
3. **Pinned memory**: Use for faster host-device copies
4. **Batched operations**: Combine multiple small kernels

---

## Building with GPU Support

### CUDA Build

```bash
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DENABLE_CUDA=ON \
  -DCMAKE_CUDA_ARCHITECTURES="70;80;86;90"

make -j$(nproc)
```

### HIP/ROCm Build

```bash
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DENABLE_HIP=ON \
  -DCMAKE_HIP_ARCHITECTURES="gfx908;gfx90a"

make -j$(nproc)
```

### Checking GPU Support

```bash
# Check if built with GPU support
./fsrm --version

# Check available GPUs
./fsrm --gpu-info
```

---

## Extending with New Kernels

### Creating a New CPU Kernel

```cpp
class MyPhysicsKernel : public PhysicsKernel {
public:
    MyPhysicsKernel() : PhysicsKernel(PhysicsType::MY_PHYSICS) {}
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    int getNumFields() const override { return 1; }
    int getNumComponents(int field) const override { return 3; }
};
```

### Creating a GPU Version

```cpp
class MyPhysicsKernelGPU : public MyPhysicsKernel {
public:
    MyPhysicsKernelGPU();
    ~MyPhysicsKernelGPU() override;
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override {
        if (!gpu_initialized) {
            MyPhysicsKernel::residual(u, u_t, u_x, a, x, f);
            return;
        }
        // GPU-accelerated implementation
    }
    
private:
    bool gpu_initialized = false;
#ifdef USE_CUDA
    double* d_solution;
    // ... device arrays
#endif
};
```

---

## References

1. Biot, M.A. (1941). General theory of three-dimensional consolidation.
2. Hughes, T.J.R. (2000). The Finite Element Method.
3. Balay et al. PETSc Users Manual. https://petsc.org/
4. NVIDIA CUDA Programming Guide.
5. AMD ROCm Documentation.
