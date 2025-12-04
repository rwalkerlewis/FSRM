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
│
├── Flow Kernels
│   ├── SinglePhaseFlowKernel
│   │   └── SinglePhaseFlowKernelGPU
│   ├── BlackOilKernel
│   │   └── BlackOilKernelGPU (planned)
│   └── CompositionalKernel
│       └── CompositionalKernelGPU (planned)
│
├── Mechanics Kernels
│   ├── GeomechanicsKernel
│   │   └── GeomechanicsKernelGPU (planned)
│   ├── ElastodynamicsKernel
│   │   └── ElastodynamicsKernelGPU
│   └── PoroelastodynamicsKernel
│       └── PoroelastodynamicsKernelGPU
│
├── Thermal Kernels
│   └── ThermalKernel
│       └── ThermalKernelGPU (planned)
│
├── Transport Kernels
│   └── ParticleTransportKernel
│       └── ParticleTransportKernelGPU (planned)
│
├── Explosion/Impact Kernels
│   ├── ExplosionSourceKernel
│   ├── NearFieldDamageKernel
│   ├── HydrodynamicKernel
│   │   └── HydrodynamicKernelGPU (planned)
│   └── CraterFormationKernel
│
├── Atmospheric Kernels
│   ├── AtmosphericBlastKernel
│   │   └── AtmosphericBlastKernelGPU
│   ├── InfrasoundKernel
│   │   └── InfrasoundKernelGPU
│   ├── ThermalRadiationKernel
│   ├── EMPKernel
│   └── FalloutKernel
│
├── Surface/Coupling Kernels
│   ├── TsunamiKernel
│   │   └── TsunamiKernelGPU
│   ├── SurfaceDeformationKernel
│   ├── FracturePropagationKernel
│   └── TidalForcesKernel
│
└── Chemistry Kernels
    └── ChemicalReactionKernel
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
| Black Oil | Planned | Planned | Three-phase flow |
| Geomechanics | Planned | Planned | Static stress/strain |
| Elastodynamics | ✓ | ✓ | Elastic wave propagation |
| Poroelastodynamics | ✓ | ✓ | Coupled fluid-solid waves |
| Thermal | Planned | Planned | Heat conduction/convection |
| Hydrodynamic | Planned | Planned | High-pressure shock flow |
| Atmospheric Blast | ✓ | ✓ | Compressible flow, blast |
| Infrasound | ✓ | ✓ | Low-frequency acoustic |
| Tsunami | ✓ | ✓ | Shallow water waves |

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

### Explosion and Impact Kernels

FSRM includes comprehensive physics models for explosion sources and impact events.

#### ExplosionSourceKernel

**Applications:**
- Nuclear test monitoring and forensics
- Chemical explosion modeling
- Industrial accident analysis

**Physics Models:**
- Mueller-Murphy source model for underground explosions
- Spherical cavity source equivalents
- Near-field damage zone evolution (cavity, crushed, fractured zones)
- Seismic moment and magnitude estimation

**Governing Equations (Cavity Pressure):**
$$\frac{\partial^2 r_c}{\partial t^2} = \frac{1}{\rho} \left( P_c(t) - P_{conf} - \frac{4G}{3} \frac{r_c - r_0}{r_c} \right)$$

**Fields:** 1 (cavity radius or equivalent moment)

**Properties:**
- `yield_kt`: Explosive yield [kt TNT]
- `depth_of_burial`: Burial depth [m]
- `fission_fraction`: Fission yield fraction
- `host_rock_properties`: Density, velocities, strength

**GPU Acceleration:** Planned

#### NearFieldDamageKernel

**Governing Equations:** Damage evolution based on strain and pressure:
$$D = 1 - \exp\left( -\alpha \langle \epsilon - \epsilon_c \rangle^+ \right)$$

**Fields:** 1 (damage parameter, 0-1)

**Features:**
- Cavity zone (D = 1): Complete vaporization/melting
- Crushed zone: Pervasive fracturing
- Fractured zone: Induced fractures, permeability enhancement
- Coupling to permeability via damage-dependent model

**GPU Acceleration:** Planned

#### HydrodynamicKernel

**Governing Equations (Euler):**
$$\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{v}) = 0$$
$$\frac{\partial (\rho \mathbf{v})}{\partial t} + \nabla \cdot (\rho \mathbf{v} \otimes \mathbf{v}) + \nabla p = \rho \mathbf{g}$$
$$\frac{\partial E}{\partial t} + \nabla \cdot ((E+p)\mathbf{v}) = \rho \mathbf{v} \cdot \mathbf{g}$$

**Fields:** 5 (density, 3 momentum components, energy)

**Properties:**
- Equation of state (ideal gas, Mie-Grüneisen, Tillotson)
- Shock capturing via artificial viscosity or limiters
- Multi-material interface tracking

**GPU Acceleration:** Planned (ideal for GPU due to local computations)

#### CraterFormationKernel

**Physics Models:**
- Pi-group crater scaling (Holsapple-Schmidt)
- Z-model excavation flow
- Shock wave attenuation in target
- Ejecta dynamics

**Key Outputs:**
- Transient and final crater dimensions
- Excavation depth and volume
- Ejecta distribution

**GPU Acceleration:** Planned

---

### Atmospheric Kernels

FSRM includes atmospheric physics for blast wave and infrasound propagation.

#### AtmosphericBlastKernel

**Governing Equations (Compressible Euler with stratification):**
$$\frac{\partial \mathbf{U}}{\partial t} + \nabla \cdot \mathbf{F}(\mathbf{U}) = \mathbf{S}$$

where $\mathbf{U} = (\rho, \rho \mathbf{v}, E)$ and $\mathbf{S}$ includes gravity and source terms.

**Fields:** 5 (density, 3 momentum, energy)

**Features:**
- Stratified atmosphere with temperature/density profiles
- Wind effects on propagation
- Ground reflection and Mach stem formation
- Thermal precursor effects

**Properties:**
- `yield_kt`: Explosion yield
- `burst_height`: Height of burst
- `atmosphere_model`: US_STANDARD_1976, ICAO, etc.
- `wind_profile`: Wind velocity vs altitude

**GPU Acceleration:** Full support via `AtmosphericBlastKernelGPU`

#### InfrasoundKernel

**Governing Equations (Linearized acoustics in moving medium):**
$$\frac{\partial p'}{\partial t} + \rho c^2 \nabla \cdot \mathbf{v}' + \mathbf{v}_0 \cdot \nabla p_0' = S$$
$$\rho \frac{\partial \mathbf{v}'}{\partial t} + \nabla p' = 0$$

**Fields:** 2 (pressure perturbation: 1, velocity perturbation: 3)

**Frequency Range:** 0.01 - 20 Hz (infrasound)

**Features:**
- Atmospheric refraction (temperature and wind effects)
- Stratospheric and thermospheric ducting
- Topographic scattering and reflection
- Molecular absorption (O₂, N₂ relaxation)

**Applications:**
- Nuclear test monitoring (CTBT/IMS)
- Volcanic eruption detection
- Bolide/meteorite tracking
- Severe weather monitoring

**Properties:**
- `atmosphere_profile`: Temperature, pressure, wind vs altitude
- `topography_model`: Digital elevation model
- `source_type`: EXPLOSION, VOLCANIC, BOLIDE, EARTHQUAKE
- `attenuation_model`: Classical + molecular relaxation

**Propagation Methods:**
1. **Ray tracing**: Fast, geometric optics approximation
2. **Parabolic equation (PE)**: Accurate for stratified media
3. **Full wave**: Finite difference in 2D/3D

**GPU Acceleration:** Full support via `InfrasoundKernelGPU`

#### ThermalRadiationKernel

**Governing Equation (Radiative transfer):**
$$\frac{1}{c}\frac{\partial I}{\partial t} + \hat{\mathbf{n}} \cdot \nabla I + \kappa I = j$$

**Fields:** 1 (intensity or thermal fluence)

**Features:**
- Fireball evolution model
- Atmospheric transmission
- Time-resolved thermal pulse
- Scaling with yield and burst height

**GPU Acceleration:** Planned

---

### Surface and Coupling Kernels

#### TsunamiKernel

**Governing Equations (Nonlinear shallow water):**
$$\frac{\partial \eta}{\partial t} + \nabla \cdot ((h + \eta) \mathbf{v}) = 0$$
$$\frac{\partial \mathbf{v}}{\partial t} + \mathbf{v} \cdot \nabla \mathbf{v} + g \nabla \eta = -C_d \frac{|\mathbf{v}| \mathbf{v}}{h + \eta}$$

**Fields:** 3 (water elevation, 2 velocity components)

**Features:**
- Seafloor deformation source
- Coriolis effects
- Variable bathymetry
- Inundation modeling

**GPU Acceleration:** Full support

#### SurfaceDeformationKernel

**Purpose:** Tracks ground surface displacement from subsurface processes.

**Coupling:**
- Receives displacement from elastodynamics/geomechanics
- Provides boundary condition for atmosphere/hydrodynamics
- Computes InSAR-observable quantities

---

## Unit System

### Design Philosophy

**All calculations are performed in SI base units:**
- Length: meters (m)
- Mass: kilograms (kg)
- Time: seconds (s)

**User interface supports any units:**
- Input values can be specified with units (e.g., "5000 psi", "100 mD")
- Output can be configured to display in preferred units
- Automatic conversion happens transparently

### Configuration

```ini
[UNITS]
# Input/output system preference
input_system = SI        # SI, FIELD, METRIC
output_system = FIELD    # SI, FIELD, METRIC

# Per-quantity output units (override system default)
pressure_unit = psi
length_unit = ft
time_unit = day
permeability_unit = mD
temperature_unit = degC
```

### Usage in Config Files

```ini
[ROCK]
permeability_x = 150 mD          # → 1.48×10⁻¹³ m² internally
youngs_modulus = 15 GPa          # → 1.5×10¹⁰ Pa internally
density = 2.55 g/cm3             # → 2550 kg/m³ internally

[SOURCE]
yield = 1.0 kt                   # → 4.184×10¹² J internally
depth = 500 ft                   # → 152.4 m internally
```

### Programmatic Usage

```cpp
#include "UnitSystem.hpp"

using namespace FSRM;

// Get unit system instance
UnitSystem& units = UnitSystemManager::getInstance();

// Convert between units
double p_pa = units.convert(5000.0, "psi", "Pa");  // psi → Pa

// Convert to/from SI base
double k_si = units.toBase(100.0, "mD");           // mD → m²
double k_field = units.fromBase(1e-15, "mD");      // m² → mD

// Parse value with unit string
double value = units.parseAndConvertToBase("500 km");  // → 500000 m
```

### Supported Unit Categories

| Category | Common Units | SI Base |
|----------|--------------|---------|
| Length | m, km, ft, mi | m |
| Mass | kg, g, lbm | kg |
| Time | s, min, hr, day | s |
| Pressure | Pa, MPa, psi, bar | Pa |
| Permeability | m², D, mD | m² |
| Viscosity | Pa·s, cP, P | Pa·s |
| Temperature | K, °C, °F | K |
| Energy | J, kJ, BTU, kt TNT | J |
| Velocity | m/s, km/s, ft/s | m/s |
| Density | kg/m³, g/cm³, lbm/ft³ | kg/m³ |

See `docs/UNIT_SYSTEM.md` for complete unit reference.

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
