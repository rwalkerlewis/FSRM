# FSRM User Guide

Complete documentation for running FSRM simulations.

## Table of Contents

1. [Overview](#overview)
2. [Running Simulations](#running-simulations)
3. [Configuration Files](#configuration-files)
4. [Physics Models](#physics-models)
5. [Boundary and Initial Conditions](#boundary-and-initial-conditions)
6. [Wells](#wells)
7. [Fractures](#fractures)
8. [Faults and Seismicity](#faults-and-seismicity)
9. [Output and Visualization](#output-and-visualization)
10. [Performance Optimization](#performance-optimization)
11. [Troubleshooting](#troubleshooting)

---

## Overview

FSRM is a fully-coupled reservoir simulator that supports:

- **Fluid Flow**: Single-phase, black oil, compositional with EOS
- **Geomechanics**: Linear elastic, viscoelastic, poroelastic, elastoplastic
- **Thermal**: Heat conduction and convection
- **Fractures**: Natural DFN, hydraulic fracturing with proppant
- **Faults**: Coulomb and rate-state friction, induced seismicity
- **Dynamics**: Elastic wave propagation, static-to-dynamic triggering

All parameters are configurable via text files - no code changes needed.

---

## Running Simulations

### Basic Execution

```bash
# Using configuration file (recommended)
mpirun -np <N> fsrm -c <config_file>

# Using Eclipse input
mpirun -np <N> fsrm -i <input.DATA> -o <output_prefix>

# Generate template configuration
fsrm -generate_config <template.config>
```

### Command-Line Options

| Option | Description |
|--------|-------------|
| `-c <file>` | Configuration file path |
| `-i <file>` | Eclipse .DATA input file |
| `-o <prefix>` | Output file prefix |
| `-format <type>` | Output format: VTK, HDF5, ECLIPSE |
| `-generate_config <file>` | Generate configuration template |
| `-help` | Show help message |

### PETSc Options

Fine-tune solvers via command line:

```bash
mpirun -np 4 fsrm -c config.config \
  -ts_type beuler \
  -snes_type newtonls \
  -snes_rtol 1e-8 \
  -ksp_type gmres \
  -pc_type ilu \
  -log_view
```

---

## Configuration Files

Configuration files use INI format with sections and key-value pairs.

### Structure

```ini
[SECTION]
key = value                  # Comment
another_key = 1.5e6          # Scientific notation supported
array_key = 1.0, 2.0, 3.0    # Comma-separated arrays
```

### Main Sections

| Section | Description |
|---------|-------------|
| `[SIMULATION]` | Time stepping, physics flags, solver options |
| `[GRID]` | Domain dimensions and mesh |
| `[ROCK]` / `[ROCK1]` | Material properties |
| `[FLUID]` | Fluid properties |
| `[WELL1]`, `[WELL2]`... | Well definitions |
| `[FRACTURE1]`... | Fracture definitions |
| `[FAULT1]`... | Fault definitions |
| `[BC1]`, `[BC2]`... | Boundary conditions |
| `[IC1]`, `[IC2]`... | Initial conditions |
| `[DYNAMICS]` | Wave propagation settings |
| `[SEISMICITY]` | Fault mechanics settings |
| `[OUTPUT]` | Output configuration |

### Example Configuration

```ini
[SIMULATION]
start_time = 0.0
end_time = 31536000.0        # 1 year
dt_initial = 86400.0         # 1 day
fluid_model = BLACK_OIL
solid_model = POROELASTIC
enable_geomechanics = true
enable_faults = true

[GRID]
nx = 50
ny = 50
nz = 20
Lx = 5000.0
Ly = 5000.0
Lz = 1000.0

[ROCK]
constitutive_model = POROELASTIC
porosity = 0.15
permeability_x = 50.0
permeability_z = 5.0
youngs_modulus = 20.0e9
poisson_ratio = 0.25
biot_coefficient = 0.8

[FLUID]
type = BLACK_OIL
oil_density_std = 850.0
solution_gor = 100.0
pvt_correlation = STANDING

[WELL1]
name = INJ1
type = INJECTOR
i = 25
j = 25
k = 10
control_mode = RATE
target_value = 0.02

[FAULT1]
name = MAIN_FAULT
x = 2500.0
y = 2500.0
z = 500.0
strike = 45.0
dip = 60.0
length = 3000.0
width = 500.0
friction_law = RATE_STATE_AGING
```

---

## Physics Models

### Fluid Models

| Model | Config Value | Description |
|-------|--------------|-------------|
| Single-phase | `SINGLE_COMPONENT` | Compressible single-phase flow |
| Black Oil | `BLACK_OIL` | Three-phase (oil/water/gas) with PVT |
| Compositional | `COMPOSITIONAL` | Multi-component with EOS flash |
| Brine | `BRINE` | Saline water with salinity effects |
| CO2 | `CO2` | Supercritical CO2 behavior |

### Solid Models

| Model | Config Value | Description |
|-------|--------------|-------------|
| Linear Elastic | `ELASTIC` | Standard Hookean elasticity |
| Viscoelastic | `VISCOELASTIC` | Maxwell, Kelvin-Voigt, or SLS |
| Poroelastic | `POROELASTIC` | Biot's coupled theory |
| Elastoplastic | `ELASTOPLASTIC` | Mohr-Coulomb or Drucker-Prager |
| Anisotropic | `VTI` / `HTI` | Transverse isotropy |

### Coupling Options

Enable via flags in `[SIMULATION]`:

```ini
enable_geomechanics = true
enable_thermal = true
enable_fractures = true
enable_faults = true
enable_elastodynamics = true
```

---

## Boundary and Initial Conditions

### Boundary Conditions

```ini
[BC1]
type = DIRICHLET              # DIRICHLET, NEUMANN, ROBIN
field = PRESSURE              # PRESSURE, TEMPERATURE, DISPLACEMENT
location = XMIN               # XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX
value = 20.0e6                # Value (Pa for pressure)

[BC2]
type = NEUMANN
field = DISPLACEMENT
location = ZMAX
value = 0.0                   # Zero traction (free surface)
```

### Initial Conditions

```ini
[IC1]
field = PRESSURE
distribution = GRADIENT       # UNIFORM, GRADIENT, FILE
value = 20.0e6                # Base value
gradient = 0, 0, 10000        # Gradient vector (Pa/m)

[IC2]
field = TEMPERATURE
distribution = GRADIENT
value = 300.0                 # Surface temperature (K)
gradient = 0, 0, 0.03         # Geothermal gradient (K/m)
```

---

## Wells

### Well Configuration

```ini
[WELL1]
name = PROD1                  # Well identifier
type = PRODUCER               # PRODUCER, INJECTOR, OBSERVATION
i = 10                        # Grid index (0-based)
j = 10
k = 5
control_mode = RATE           # RATE, BHP, THP
target_value = 0.01           # m³/s (rate) or Pa (pressure)
max_rate = 0.1                # Maximum rate constraint
min_bhp = 5.0e6               # Minimum BHP constraint
diameter = 0.2                # Wellbore diameter (m)
skin = 0.0                    # Skin factor
```

### Well Types

| Type | Description |
|------|-------------|
| `PRODUCER` | Production well (fluid removal) |
| `INJECTOR` | Injection well (fluid addition) |
| `OBSERVATION` | Monitoring well (no flow) |

### Control Modes

| Mode | Description |
|------|-------------|
| `RATE` | Fixed rate (m³/s), limited by BHP |
| `BHP` | Fixed bottomhole pressure |
| `THP` | Fixed tubing head pressure |

---

## Fractures

### Natural Fractures

```ini
[FRACTURE1]
type = NATURAL
location = 500, 500, 50, 0.707, 0.707, 0  # x,y,z, normal_x,y,z
aperture = 0.0001             # meters
permeability = 1.0e-12        # m²
enable_propagation = false
```

### Hydraulic Fractures

```ini
[FRACTURE2]
type = HYDRAULIC
location = 250, 500, 50, 0, 1, 0
aperture = 0.001              # meters
permeability = 1.0e-10        # m²
toughness = 1.0e6             # Pa·m^0.5
energy = 100.0                # J/m²
enable_propagation = true
enable_proppant = true
```

---

## Faults and Seismicity

### Fault Definition

```ini
[FAULT1]
name = MAIN_FAULT
x = 5000.0                    # Center X (m)
y = 0.0                       # Center Y
z = 4000.0                    # Center depth
strike = 0.0                  # Degrees from North
dip = 60.0                    # Degrees from horizontal
length = 2000.0               # Along-strike (m)
width = 1500.0                # Down-dip (m)
```

### Friction Laws

```ini
# Coulomb friction
friction_law = COULOMB
static_friction = 0.6
dynamic_friction = 0.4
cohesion = 1.0e6              # Pa

# Rate-and-state friction
friction_law = RATE_STATE_AGING
rate_state_a = 0.010          # Direct effect
rate_state_b = 0.015          # Evolution effect (b>a = unstable)
rate_state_dc = 1.0e-4        # Critical slip distance (m)
rate_state_v0 = 1.0e-6        # Reference velocity (m/s)
rate_state_f0 = 0.6           # Reference friction
```

### Seismicity Configuration

```ini
[SEISMICITY]
enable = true
friction_law = RATE_STATE_AGING
nucleation_size = 1.0         # Critical patch size (m)
seismic_slip_rate = 1.0e-3    # Seismic threshold (m/s)
b_value = 1.0                 # Gutenberg-Richter b-value
aftershocks = true
stress_transfer = true
```

---

## Output and Visualization

### Output Configuration

```ini
[OUTPUT]
format = VTK                  # VTK, HDF5, ECLIPSE
path = output                 # Directory
frequency = 10                # Every N timesteps

# Fields to output
pressure = true
displacement = true
stress = false
strain = false
velocity = false
permeability = false
temperature = false
saturation = false
fault_slip = true
seismic_catalog = true
```

### Visualization

```bash
# View with ParaView
paraview output/*.vtu

# View seismic catalog
cat output/seismic_catalog.csv
```

---

## Performance Optimization

### Parallel Scaling

```bash
# Scale processes with problem size
# ~10,000 cells per process is optimal
mpirun -np 8 fsrm -c large_simulation.config
```

### Solver Options

```ini
[SIMULATION]
rtol = 1.0e-6
atol = 1.0e-8
max_nonlinear_iterations = 50
max_linear_iterations = 1000
```

### GPU Acceleration

```ini
[SIMULATION]
use_gpu = true
gpu_mode = CPU_FALLBACK
gpu_device_id = 0
gpu_memory_fraction = 0.8
```

```bash
# Build with CUDA
cmake .. -DENABLE_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES="80;86"
make
```

---

## Troubleshooting

### Common Issues

**PETSc not found:**
```bash
export PETSC_DIR=/path/to/petsc
export PETSC_ARCH=arch-linux-c-debug
```

**Convergence failure:**
- Reduce timestep: `dt_initial = 100.0`
- Increase iterations: `max_nonlinear_iterations = 100`
- Check boundary conditions for consistency

**Memory issues:**
- Reduce grid resolution
- Use more MPI processes
- Enable GPU if available

**Slow performance:**
- Increase grid coarseness
- Use adaptive timestepping
- Enable AMG preconditioner: `-pc_type hypre -pc_hypre_type boomeramg`

### Debug Mode

```bash
# Build with debug symbols
cmake .. -DCMAKE_BUILD_TYPE=Debug
make

# Run with verbose output
mpirun -np 4 fsrm -c config.config \
  -snes_monitor \
  -ksp_monitor \
  -ts_monitor \
  -log_view
```

---

## References

- PETSc: https://petsc.org/
- Biot poroelasticity: Wang (2000), Theory of Linear Poroelasticity
- Rate-state friction: Dieterich (1979), Ruina (1983)
- Black oil PVT: Standing (1947), Vazquez-Beggs (1980)
