# Tutorial 02: Configuration System

Master the FSRM configuration file format to set up any simulation without code changes.

## Table of Contents

1. [Configuration File Format](#configuration-file-format)
2. [Section Reference](#section-reference)
3. [Simulation Section](#simulation-section)
4. [Grid Section](#grid-section)
5. [Rock Properties](#rock-properties)
6. [Fluid Properties](#fluid-properties)
7. [Boundary Conditions](#boundary-conditions)
8. [Initial Conditions](#initial-conditions)
9. [Output Configuration](#output-configuration)
10. [Unit System](#unit-system)
11. [Tips and Best Practices](#tips-and-best-practices)

---

## Configuration File Format

FSRM uses an INI-style configuration format that's easy to read and edit.

### Basic Syntax

```ini
# This is a comment
; This is also a comment

[SECTION_NAME]
key = value                    # Inline comments allowed
another_key = 1.5e6            # Scientific notation
array_key = 1.0, 2.0, 3.0      # Comma-separated arrays
string_key = SOME_STRING       # Strings (no quotes needed)
```

### Key Rules

1. **Sections** are enclosed in brackets: `[SECTION]`
2. **Keys** are case-insensitive: `Porosity = 0.2` equals `porosity = 0.2`
3. **Values** are parsed based on context (integer, float, string, array)
4. **Comments** start with `#` or `;`
5. **Whitespace** around `=` is ignored

### Data Types

| Type | Examples |
|------|----------|
| Integer | `10`, `1000`, `-5` |
| Float | `0.25`, `1.5e-6`, `3.14159` |
| Boolean | `true`, `false`, `yes`, `no`, `1`, `0` |
| String | `PRODUCER`, `VTK`, `GMRES` |
| Array | `1.0, 2.0, 3.0` or `0.1 0.2 0.3` |

---

## Section Reference

### Complete Section List

| Section | Purpose | Multiple Allowed |
|---------|---------|------------------|
| `[SIMULATION]` | Core simulation parameters | No |
| `[GRID]` | Domain and mesh settings | No |
| `[ROCK]` | Default rock properties | No |
| `[ROCK1]`, `[ROCK2]`... | Regional rock properties | Yes |
| `[FLUID]` | Fluid model and properties | No |
| `[WELL1]`, `[WELL2]`... | Well definitions | Yes |
| `[FRACTURE1]`, `[FRACTURE2]`... | Fracture definitions | Yes |
| `[FAULT1]`, `[FAULT2]`... | Fault definitions | Yes |
| `[BC1]`, `[BC2]`... | Boundary conditions | Yes |
| `[IC1]`, `[IC2]`... | Initial conditions | Yes |
| `[DYNAMICS]` | Wave propagation settings | No |
| `[SEISMICITY]` | Fault mechanics settings | No |
| `[AMR]` | Adaptive mesh refinement | No |
| `[OUTPUT]` | Output configuration | No |

---

## Simulation Section

The `[SIMULATION]` section controls core simulation behavior.

### Time Stepping

```ini
[SIMULATION]
# Simulation identification
name = my_simulation           # Identifier for output files

# Time parameters (in seconds unless units specified)
start_time = 0.0               # Simulation start
end_time = 86400.0             # End time (1 day = 86400 s)
dt_initial = 3600.0            # Initial timestep (1 hour)
dt_min = 1.0                   # Minimum allowed timestep
dt_max = 86400.0               # Maximum allowed timestep
adaptive_dt = true             # Enable adaptive timestepping
```

### Physics Selection

```ini
[SIMULATION]
# Fluid model options:
# SINGLE_COMPONENT, BLACK_OIL, COMPOSITIONAL, BRINE, CO2
fluid_model = SINGLE_COMPONENT

# Solid model options:
# ELASTIC, VISCOELASTIC, POROELASTIC, ELASTOPLASTIC, VTI, HTI
solid_model = ELASTIC

# Physics coupling flags
enable_geomechanics = false    # Enable mechanical deformation
enable_thermal = false         # Enable heat transfer
enable_fractures = false       # Enable fracture modeling
enable_faults = false          # Enable fault mechanics
enable_elastodynamics = false  # Enable elastic wave propagation
enable_poroelastodynamics = false  # Enable coupled wave propagation
```

### Solver Settings

```ini
[SIMULATION]
# Nonlinear solver tolerances
rtol = 1.0e-6                  # Relative tolerance
atol = 1.0e-8                  # Absolute tolerance
max_nonlinear_iterations = 50  # Max SNES iterations

# Linear solver settings
max_linear_iterations = 1000   # Max KSP iterations
solver_type = gmres            # Options: gmres, cg, bicgstab, preonly
preconditioner = ilu           # Options: ilu, jacobi, asm, hypre, lu
```

### GPU Acceleration

```ini
[SIMULATION]
# GPU settings
use_gpu = false                # Enable GPU computation
gpu_mode = CPU_FALLBACK        # GPU_ONLY, CPU_FALLBACK, AUTO
gpu_device_id = 0              # Which GPU to use
gpu_memory_fraction = 0.8      # Max fraction of GPU memory
```

---

## Grid Section

Define the computational domain and mesh.

### Cartesian Grid

```ini
[GRID]
# All grids use DMPlex for unstructured mesh representation
# This enables consistent mesh operations and AMR capabilities

# Grid type (defines how the mesh is created/loaded)
grid_type = CARTESIAN          # CARTESIAN, GMSH, CORNER_POINT, EXODUS

# Number of cells in each direction (for box meshes)
nx = 50                        # X cells
ny = 50                        # Y cells
nz = 20                        # Z cells

# Domain dimensions (meters)
Lx = 5000.0                    # Domain length X
Ly = 5000.0                    # Domain length Y
Lz = 1000.0                    # Domain length Z

# Domain origin (default 0,0,0)
origin_x = 0.0
origin_y = 0.0
origin_z = 0.0

# All grids use DMPlex internally (required for all simulations)
use_unstructured = true

# Uniform refinement (0 = no refinement)
refinement_level = 0
```

> **Important:** All problems use DMPlex for the spatial domain. This ensures consistent unstructured mesh handling across all simulation types.

### External Mesh File

```ini
[GRID]
grid_type = GMSH               # Use external Gmsh mesh
mesh_file = mesh/reservoir.msh # Path to mesh file
use_unstructured = true        # Required for all grids

# Coordinate transformations
coordinate_system = CARTESIAN  # CARTESIAN, CYLINDRICAL, SPHERICAL
```

### Grid Resolution Guidelines

| Application | Typical Resolution |
|-------------|-------------------|
| Regional flow | 100-500m cells |
| Well vicinity | 1-10m cells |
| Fracture zone | 0.1-1m cells |
| Wave propagation | λ/10 minimum |

---

## Rock Properties

### Default Properties

The `[ROCK]` section sets default properties for the entire domain:

```ini
[ROCK]
# Constitutive model
# Options: LINEAR_ELASTIC, VISCOELASTIC_MAXWELL, VISCOELASTIC_KELVIN,
#          VISCOELASTIC_SLS, POROELASTIC, ELASTOPLASTIC_MC, 
#          ELASTOPLASTIC_DP, VTI, HTI
constitutive_model = LINEAR_ELASTIC

# Basic properties
porosity = 0.20                # Fraction (0-1)
permeability_x = 100.0         # milliDarcy
permeability_y = 100.0         # milliDarcy (default = perm_x)
permeability_z = 10.0          # milliDarcy (often lower)
density = 2650.0               # Grain density kg/m³
compressibility = 1.0e-9       # 1/Pa
```

### Elastic Properties

```ini
[ROCK]
# Method 1: Young's modulus and Poisson's ratio
youngs_modulus = 20.0e9        # Pa (20 GPa)
poisson_ratio = 0.25

# Method 2: Bulk and shear modulus (alternative)
# bulk_modulus = 13.3e9        # Pa
# shear_modulus = 8.0e9        # Pa
```

### Poroelastic Properties

For coupled flow-mechanics:

```ini
[ROCK]
constitutive_model = POROELASTIC

# Biot parameters
biot_coefficient = 0.8         # α (0-1, typically 0.7-1.0)
biot_modulus = 1.0e10          # M (Pa)

# Alternative specification
# undrained_poisson = 0.35     # νu (undrained Poisson ratio)
# skempton = 0.9               # B (Skempton coefficient)
```

### Viscoelastic Properties

For time-dependent deformation:

```ini
[ROCK]
constitutive_model = VISCOELASTIC_MAXWELL

viscoelastic_type = MAXWELL    # MAXWELL, KELVIN_VOIGT, SLS
viscosity = 1.0e18             # Pa·s
relaxation_time = 1.0e12       # s
```

### Plasticity Properties

For rock failure:

```ini
[ROCK]
constitutive_model = ELASTOPLASTIC_MC

failure_criterion = MOHR_COULOMB  # MOHR_COULOMB, DRUCKER_PRAGER
cohesion = 1.0e6               # Pa
friction_angle = 30.0          # degrees
dilation_angle = 0.0           # degrees
tensile_strength = 0.0         # Pa
```

### Anisotropic Properties (VTI/HTI)

For layered or fractured media:

```ini
[ROCK]
constitutive_model = VTI       # Vertical Transverse Isotropy

# Method 1: Full stiffness tensor components
c11 = 50.0e9                   # Pa
c12 = 15.0e9
c13 = 12.0e9
c33 = 40.0e9
c44 = 12.0e9

# Method 2: Thomsen parameters
# thomsen_epsilon = 0.1
# thomsen_delta = 0.05
# thomsen_gamma = 0.08
```

### Regional Properties

Define different rock types for different regions:

```ini
[ROCK1]
name = sandstone
region = 1                     # Region identifier
porosity = 0.25
permeability_x = 200.0
youngs_modulus = 15.0e9

[ROCK2]
name = shale
region = 2
porosity = 0.05
permeability_x = 0.01
youngs_modulus = 30.0e9
```

### Thermal Properties

```ini
[ROCK]
thermal_conductivity = 2.5     # W/(m·K)
specific_heat = 800.0          # J/(kg·K)
thermal_expansion = 1.0e-5     # 1/K
```

---

## Fluid Properties

### Single-Phase

```ini
[FLUID]
type = SINGLE_PHASE

density = 1000.0               # kg/m³ at reference conditions
viscosity = 0.001              # Pa·s (1 cP)
compressibility = 4.5e-10      # 1/Pa
reference_pressure = 1.0e5     # Pa (1 atm)
```

### Black Oil Model

For oil/gas/water systems:

```ini
[FLUID]
type = BLACK_OIL

# Standard condition densities
oil_density_std = 850.0        # kg/m³
gas_density_std = 0.9          # kg/m³
water_density_std = 1000.0     # kg/m³

# Solution gas-oil ratio
solution_gor = 100.0           # scf/stb or sm³/sm³
bubble_point = 15.0e6          # Pa

# PVT correlation
pvt_correlation = STANDING     # STANDING, VASQUEZ_BEGGS, GLASO
```

### Compositional Model

For multi-component systems:

```ini
[FLUID]
type = COMPOSITIONAL

num_components = 3
component_names = C1, C3, C7+

# Critical properties
critical_temperatures = 190.6, 369.8, 560.0    # K
critical_pressures = 4.6e6, 4.25e6, 2.4e6      # Pa
acentric_factors = 0.011, 0.152, 0.35
molar_weights = 16.04, 44.1, 100.0             # kg/kmol

# Equation of state
eos_type = PENG_ROBINSON       # PENG_ROBINSON, SRK, VAN_DER_WAALS
```

### Brine

For saline water systems:

```ini
[FLUID]
type = BRINE

salinity = 50000               # ppm
salt_type = NaCl               # Salt composition
```

### CO2

For CO2 storage simulations:

```ini
[FLUID]
type = CO2

co2_model = SPAN_WAGNER        # Property correlation
```

---

## Boundary Conditions

Define boundary conditions using numbered sections:

### Pressure (Dirichlet)

```ini
[BC1]
type = DIRICHLET               # Fixed value
field = PRESSURE               # Which field
location = XMIN                # XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX
value = 20.0e6                 # Value (Pa)
```

### No-Flow (Neumann)

```ini
[BC2]
type = NEUMANN                 # Fixed flux
field = PRESSURE
location = XMAX
value = 0.0                    # Zero flux = no flow
```

### Temperature

```ini
[BC3]
type = DIRICHLET
field = TEMPERATURE
location = ZMIN
value = 400.0                  # Bottom temperature (K)

[BC4]
type = DIRICHLET
field = TEMPERATURE
location = ZMAX
value = 300.0                  # Top temperature (K)
```

### Displacement

```ini
[BC5]
type = DIRICHLET
field = DISPLACEMENT
location = ZMIN
value = 0.0                    # Fixed base
component = z                  # Which component (x, y, z, all)

[BC6]
type = NEUMANN
field = DISPLACEMENT
location = ZMAX
value = 0.0                    # Zero traction (free surface)
```

### Time-Dependent Boundaries

```ini
[BC7]
type = DIRICHLET
field = PRESSURE
location = XMIN
value = 25.0e6

# Time function options
time_function = ramp           # constant, ramp, step, sinusoidal, table
ramp_start = 0.0               # Start time
ramp_end = 3600.0              # End time
ramp_start_value = 20.0e6      # Initial value
ramp_end_value = 25.0e6        # Final value
```

---

## Initial Conditions

### Uniform Distribution

```ini
[IC1]
field = PRESSURE
distribution = UNIFORM
value = 20.0e6                 # 20 MPa everywhere
```

### Gradient Distribution

```ini
[IC2]
field = PRESSURE
distribution = GRADIENT
value = 20.0e6                 # Base value at origin
gradient = 0.0, 0.0, 10000.0   # Pa/m in (x, y, z)
```

### Temperature with Geothermal Gradient

```ini
[IC3]
field = TEMPERATURE
distribution = GRADIENT
value = 300.0                  # Surface temperature (K)
gradient = 0.0, 0.0, 0.03      # 30°C/km geothermal gradient
```

### From File

```ini
[IC4]
field = PRESSURE
distribution = FILE
file = initial_pressure.dat    # External data file
```

### Saturation (Multi-phase)

```ini
[IC5]
field = SATURATION
phase = WATER
distribution = UNIFORM
value = 0.3                    # 30% water saturation
```

---

## Output Configuration

```ini
[OUTPUT]
# Format: VTK, HDF5, ECLIPSE, ASCII
format = VTK

# Output directory
path = output

# Output frequency
frequency = 10                 # Every 10 timesteps
# Alternative: time-based output
# time_interval = 3600.0       # Every hour

# Fields to output (true/false)
pressure = true
displacement = true
velocity = false
stress = false
strain = false
permeability = false
temperature = false
saturation = false
fault_slip = false
seismic_catalog = false

# Output precision
precision = double             # float or double
compression = true             # Enable compression
```

---

## Unit System

FSRM supports flexible unit specification:

### Default Units (SI)

| Quantity | Default Unit |
|----------|-------------|
| Length | meters (m) |
| Time | seconds (s) |
| Mass | kilograms (kg) |
| Pressure | Pascals (Pa) |
| Temperature | Kelvin (K) |
| Permeability | milliDarcy (mD) |

### Specifying Units

```ini
[SIMULATION]
end_time = 365 days            # Alternative to seconds

[GRID]
Lx = 5 km                      # Instead of 5000.0
Lz = 500 ft                    # Feet also work

[ROCK]
youngs_modulus = 20 GPa        # Instead of 20.0e9
permeability_x = 100 mD        # Explicit units

[FLUID]
viscosity = 1 cP               # Centipoise
```

### Supported Unit Conversions

```ini
# Length
m, km, ft, mile, in

# Time  
s, min, hour, day, year

# Pressure
Pa, kPa, MPa, GPa, bar, psi, atm

# Temperature
K, C, F

# Permeability
mD, D, m2

# Viscosity
Pa.s, cP, mPa.s
```

### Using Unit-Aware Config

Create a config with explicit units:

```ini
# config/with_units_example.config

[SIMULATION]
end_time = 1 year
dt_initial = 1 day

[GRID]
Lx = 10 km
Ly = 10 km
Lz = 2 km

[ROCK]
youngs_modulus = 30 GPa
permeability_x = 50 mD

[FLUID]
viscosity = 0.5 cP
reference_pressure = 1 atm
```

---

## Tips and Best Practices

### 1. Start from Templates

```bash
# Generate a complete template
./fsrm -generate_config my_template.config

# Or copy an existing config
cp config/complete_template.config my_simulation.config
```

### 2. Use Comments Liberally

```ini
# ============================================
# Reservoir Flow Simulation - Field A
# Date: 2024-01-15
# Author: Your Name
# ============================================

[SIMULATION]
name = field_a_production      # Project identifier

# Time settings for 5-year forecast
end_time = 157680000.0         # 5 years in seconds
```

### 3. Validate Before Running

```bash
# Check config syntax (dry run)
./fsrm -c config.config -validate

# Or just run briefly
./fsrm -c config.config -ts_max_steps 1
```

### 4. Use Include Files (Advanced)

```ini
# main.config
[SIMULATION]
name = main_simulation

# Include rock properties from another file
include = rock_properties.config

# Include well definitions
include = wells.config
```

### 5. Common Mistakes to Avoid

```ini
# WRONG: Spaces in section names
[ROCK 1]                       # Use [ROCK1]

# WRONG: Missing sections
[WELL1]
type = PRODUCER
# Missing i, j, k location!

# WRONG: Inconsistent units
permeability_x = 100           # Is this mD or m²?

# RIGHT: Be explicit
permeability_x = 100.0 mD      # Clear units
```

### 6. Parameter Studies

Create a base config and modify programmatically:

```python
# Python script for parameter studies
import subprocess

base_config = """
[SIMULATION]
end_time = 86400.0

[ROCK]
permeability_x = {perm}

[OUTPUT]
path = output/perm_{perm}
"""

for perm in [10, 50, 100, 200, 500]:
    config = base_config.format(perm=perm)
    with open(f'temp_{perm}.config', 'w') as f:
        f.write(config)
    subprocess.run(['mpirun', '-np', '4', './fsrm', 
                    '-c', f'temp_{perm}.config'])
```

---

## Complete Example Configuration

Here's a complete, well-documented configuration file:

```ini
# ==============================================================================
# FSRM Configuration: Reservoir Pressure Depletion Study
# ==============================================================================
# Description: 5-year production simulation with single producer
# Author: Tutorial Example
# ==============================================================================

# ------------------------------------------------------------------------------
# Simulation Control
# ------------------------------------------------------------------------------
[SIMULATION]
name = pressure_depletion_study

# 5-year simulation
start_time = 0.0
end_time = 157680000.0         # 5 years = 5*365*24*3600 seconds

# Time stepping
dt_initial = 86400.0           # Start with 1-day steps
dt_min = 3600.0                # Minimum 1 hour
dt_max = 2592000.0             # Maximum 30 days
adaptive_dt = true

# Physics
fluid_model = SINGLE_COMPONENT
solid_model = ELASTIC
enable_geomechanics = true     # Track compaction

# Solver
rtol = 1.0e-6
atol = 1.0e-8

# ------------------------------------------------------------------------------
# Domain Definition (all grids use DMPlex internally)
# ------------------------------------------------------------------------------
[GRID]
grid_type = CARTESIAN
nx = 50
ny = 50
nz = 10
Lx = 2000.0                    # 2 km × 2 km × 100 m domain
Ly = 2000.0
Lz = 100.0
use_unstructured = true        # All grids use DMPlex

# ------------------------------------------------------------------------------
# Rock Properties
# ------------------------------------------------------------------------------
[ROCK]
constitutive_model = POROELASTIC

# Flow properties
porosity = 0.18
permeability_x = 150.0         # mD
permeability_y = 150.0
permeability_z = 15.0          # 10:1 anisotropy

# Mechanical properties
youngs_modulus = 20.0e9        # 20 GPa
poisson_ratio = 0.25
density = 2600.0

# Poroelastic coupling
biot_coefficient = 0.85
biot_modulus = 8.0e9

# ------------------------------------------------------------------------------
# Fluid Properties
# ------------------------------------------------------------------------------
[FLUID]
type = SINGLE_PHASE
density = 800.0                # Light oil, kg/m³
viscosity = 0.002              # 2 cP
compressibility = 1.5e-9       # 1/Pa
reference_pressure = 25.0e6    # 25 MPa

# ------------------------------------------------------------------------------
# Wells
# ------------------------------------------------------------------------------
[WELL1]
name = PROD-1
type = PRODUCER
i = 25                         # Center of domain
j = 25
k = 5
control_mode = RATE
target_value = 0.01            # 0.01 m³/s ≈ 860 bbl/day
min_bhp = 10.0e6               # Minimum BHP constraint
diameter = 0.2
skin = 0.0

# ------------------------------------------------------------------------------
# Boundary Conditions
# ------------------------------------------------------------------------------
[BC1]
# No-flow on all boundaries (default)
type = NEUMANN
field = PRESSURE
location = ALL
value = 0.0

[BC2]
# Fixed displacement at base
type = DIRICHLET
field = DISPLACEMENT
location = ZMIN
value = 0.0
component = all

# ------------------------------------------------------------------------------
# Initial Conditions
# ------------------------------------------------------------------------------
[IC1]
field = PRESSURE
distribution = GRADIENT
value = 25.0e6                 # 25 MPa at surface
gradient = 0.0, 0.0, 10000.0   # Hydrostatic gradient

# ------------------------------------------------------------------------------
# Output
# ------------------------------------------------------------------------------
[OUTPUT]
format = VTK
path = output/pressure_depletion
frequency = 30                 # Monthly output (approx)
pressure = true
displacement = true
stress = true
permeability = true
```

---

**Previous**: [← Getting Started](01_GETTING_STARTED.md) | **Next**: [Physics Models →](03_PHYSICS_MODELS.md)
