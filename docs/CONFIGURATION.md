# Configuration Reference

FSRM simulations are driven entirely by configuration files, eliminating the need to write custom C++ code for most use cases. This document provides a complete reference for all configuration options.

## File Format

Configuration files use INI-style syntax with sections and key-value pairs:

```ini
# Comment
[SECTION_NAME]
key = value
string_key = "quoted string"
list_key = value1, value2, value3
```

## Core Sections

### [SIMULATION]

Controls overall simulation behavior.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `name` | string | "Simulation" | Simulation identifier |
| `type` | enum | RESERVOIR | RESERVOIR, GEOMECHANICS, WAVE_PROPAGATION, COUPLED |
| `start_time` | double | 0.0 | Start time (seconds or days) |
| `end_time` | double | 1.0 | End time |
| `dt` | double | 0.001 | Initial time step |
| `max_timesteps` | int | 1000 | Maximum time steps |
| `output_frequency` | int | 10 | Steps between outputs |

### [TIME]

Advanced time stepping options.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `method` | enum | BACKWARD_EULER | FORWARD_EULER, BACKWARD_EULER, CRANK_NICOLSON, BDF2 |
| `adaptive` | bool | true | Enable adaptive time stepping |
| `min_dt` | double | 1e-10 | Minimum time step |
| `max_dt` | double | 1e6 | Maximum time step |
| `cfl` | double | 0.5 | CFL number (explicit methods) |
| `dt_growth_factor` | double | 1.2 | Max dt increase per step |
| `dt_reduction_factor` | double | 0.5 | dt reduction on failure |

### [GRID]

Mesh and domain definition.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `nx`, `ny`, `nz` | int | 10 | Grid cells in each direction |
| `Lx`, `Ly`, `Lz` | double | 1000.0 | Domain size (meters) |
| `origin_x`, `origin_y`, `origin_z` | double | 0.0 | Domain origin |
| `mesh_type` | enum | CARTESIAN | CARTESIAN, GMSH, CORNER_POINT, EXODUS |
| `mesh_file` | string | "" | External mesh file path |
| `use_unstructured` | bool | false | Use unstructured grid |
| `input_crs` | string | "" | Input coordinate system (EPSG) |
| `model_crs` | string | "" | Model coordinate system (EPSG) |
| `use_local_coordinates` | bool | true | Apply local origin offset |
| `local_origin_x`, `local_origin_y`, `local_origin_z` | double | 0.0 | Local origin |
| `auto_detect_utm` | bool | false | Auto-detect UTM zone |

### [PHYSICS]

Enable/disable physics modules.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `enable_flow` | bool | true | Enable fluid flow |
| `enable_transport` | bool | false | Enable species transport |
| `enable_geomechanics` | bool | false | Enable geomechanics |
| `enable_thermal` | bool | false | Enable heat transfer |
| `enable_fractures` | bool | false | Enable discrete fractures |
| `enable_faults` | bool | false | Enable fault slip |

## Material Sections

### [FLUID]

Fluid properties.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `type` | enum | SINGLE_COMPONENT | SINGLE_COMPONENT, BLACK_OIL, COMPOSITIONAL |
| `density` | double | 1000.0 | Reference density (kg/m³) |
| `viscosity` | double | 0.001 | Viscosity (Pa·s) |
| `compressibility` | double | 4.5e-10 | Compressibility (1/Pa) |
| `Bo` | double | 1.0 | Oil formation volume factor |
| `Bw` | double | 1.0 | Water formation volume factor |
| `Bg` | double | 1.0 | Gas formation volume factor |
| `Rs` | double | 0.0 | Solution gas-oil ratio |
| `mu_o` | double | 0.001 | Oil viscosity |
| `mu_w` | double | 0.001 | Water viscosity |
| `mu_g` | double | 0.00001 | Gas viscosity |
| `relperm_model` | enum | COREY | LINEAR, COREY, BROOKS_COREY |
| `n_w` | double | 2.0 | Water Corey exponent |
| `n_o` | double | 2.0 | Oil Corey exponent |
| `S_wc` | double | 0.2 | Connate water saturation |
| `S_or` | double | 0.2 | Residual oil saturation |
| `k_rw_max` | double | 0.3 | Max water relative perm |
| `k_ro_max` | double | 1.0 | Max oil relative perm |

### [ROCK]

Rock/matrix properties.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `porosity` | double | 0.2 | Porosity |
| `permeability` | double | 1e-13 | Isotropic permeability (m²) |
| `permeability_x`, `permeability_y`, `permeability_z` | double | 1e-13 | Anisotropic permeability |
| `compressibility` | double | 1e-10 | Rock compressibility (1/Pa) |
| `density` | double | 2650 | Grain density (kg/m³) |
| `thermal_conductivity` | double | 2.5 | Thermal conductivity (W/m·K) |
| `specific_heat` | double | 800 | Specific heat (J/kg·K) |

### [SOLID]

Geomechanical properties.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `model` | enum | ELASTIC | ELASTIC, POROELASTIC, VISCOELASTIC |
| `E` | double | 30e9 | Young's modulus (Pa) |
| `nu` | double | 0.25 | Poisson's ratio |
| `density` | double | 2650 | Solid density (kg/m³) |
| `biot_coefficient` | double | 0.8 | Biot coefficient |
| `lambda` | double | - | Lamé's first parameter |
| `mu` | double | - | Shear modulus |
| `Vp` | double | - | P-wave velocity (m/s) |
| `Vs` | double | - | S-wave velocity (m/s) |

## Wells

Define wells using numbered sections `[WELL1]`, `[WELL2]`, etc.

| Key | Type | Description |
|-----|------|-------------|
| `name` | string | Well identifier |
| `type` | enum | INJECTOR, PRODUCER |
| `x`, `y`, `z` | double | Well location |
| `i`, `j`, `k` | int | Well cell indices (alternative) |
| `radius` | double | Wellbore radius |
| `skin` | double | Skin factor |
| `control_mode` | enum | RATE, BHP |
| `target_value` | double | Rate (m³/s) or BHP (Pa) |
| `fluid` | enum | OIL, WATER, GAS |
| `perforation_start` | double | Top of perforation |
| `perforation_end` | double | Bottom of perforation |

Example:
```ini
[WELL1]
name = INJECTOR-1
type = INJECTOR
x = 100.0
y = 250.0
z = -1500.0
radius = 0.1
control_mode = RATE
target_value = 0.001
fluid = WATER

[WELL2]
name = PRODUCER-1
type = PRODUCER
x = 900.0
y = 250.0
z = -1500.0
control_mode = BHP
target_value = 10.0e6
fluid = OIL
```

## Boundary Conditions

Define BCs using numbered sections `[BC1]`, `[BC2]`, etc.

| Key | Type | Description |
|-----|------|-------------|
| `type` | enum | DIRICHLET, NEUMANN, ROBIN |
| `field` | enum | PRESSURE, TEMPERATURE, DISPLACEMENT, SATURATION |
| `location` | string | XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, or physical group name |
| `value` | double | BC value |
| `component` | int | Vector component (0=x, 1=y, 2=z) |

Example:
```ini
[BC1]
type = DIRICHLET
field = PRESSURE
location = XMIN
value = 20.0e6

[BC2]
type = NEUMANN
field = PRESSURE
location = XMAX, YMIN, YMAX
value = 0.0

[BC3]
type = DIRICHLET
field = DISPLACEMENT
location = ZMIN
value = 0.0
component = 2
```

## Initial Conditions

Define ICs using numbered sections `[IC1]`, `[IC2]`, etc.

| Key | Type | Description |
|-----|------|-------------|
| `field` | enum | PRESSURE, SATURATION_WATER, SATURATION_OIL, TEMPERATURE, etc. |
| `type` | enum | CONSTANT, LINEAR, FUNCTION |
| `value` | double | Constant value |
| `gradient_x`, `gradient_y`, `gradient_z` | double | Gradient components |

Example:
```ini
[IC1]
field = PRESSURE
type = CONSTANT
value = 15.0e6

[IC2]
field = SATURATION_WATER
type = CONSTANT
value = 0.2

[IC3]
field = TEMPERATURE
type = LINEAR
value = 60.0           # At reference point
gradient_z = 0.03      # 30°C/km geothermal gradient
```

## Solver Settings

### [SOLVER]

Linear and nonlinear solver configuration.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `linear_solver` | enum | GMRES | CG, GMRES, BICGSTAB, DIRECT |
| `preconditioner` | enum | ILU | JACOBI, ILU, ASM, GAMG, HYPRE |
| `max_linear_iterations` | int | 1000 | Max linear iterations |
| `linear_tolerance` | double | 1e-8 | Linear solver tolerance |
| `nonlinear_solver` | enum | NEWTON | NEWTON, PICARD |
| `max_nonlinear_iterations` | int | 20 | Max Newton iterations |
| `nonlinear_tolerance` | double | 1e-6 | Nonlinear tolerance |
| `line_search` | bool | true | Enable line search |

### [GPU]

GPU acceleration settings.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `enable` | bool | false | Enable GPU |
| `device_id` | int | 0 | GPU device index |
| `streams` | int | 4 | CUDA streams |

## Output Settings

### [OUTPUT]

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `format` | enum | VTK | VTK, HDF5, ECLIPSE |
| `directory` | string | "output" | Output directory |
| `prefix` | string | "solution" | File prefix |
| `fields` | list | all | Fields to output |
| `binary` | bool | true | Binary output |
| `compression` | bool | false | Enable compression |

Example:
```ini
[OUTPUT]
format = VTK
directory = results
prefix = reservoir_sim
fields = PRESSURE, SATURATION, VELOCITY
binary = true
```

## Example Configurations

### Single-Phase Flow

```ini
[SIMULATION]
name = SinglePhaseReservoir
type = RESERVOIR
end_time = 86400
output_frequency = 100

[GRID]
nx = 50
ny = 50
nz = 10
Lx = 1000
Ly = 500
Lz = 50

[FLUID]
type = SINGLE_COMPONENT
density = 800
viscosity = 0.001
compressibility = 5e-10

[ROCK]
porosity = 0.2
permeability = 1e-13

[WELL1]
name = INJECTOR
type = INJECTOR
x = 100
y = 250
z = -25
control_mode = RATE
target_value = 0.001

[WELL2]
name = PRODUCER
type = PRODUCER
x = 900
y = 250
z = -25
control_mode = BHP
target_value = 10e6

[BC1]
type = NEUMANN
field = PRESSURE
location = XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX
value = 0.0

[IC1]
field = PRESSURE
type = CONSTANT
value = 20e6

[OUTPUT]
format = VTK
directory = output/single_phase
```

### Two-Phase Waterflooding

See `/workspace/config/buckley_leverett_2d.config`

### Coupled Flow-Geomechanics

See `/workspace/config/coupled_reservoir_2d.config`

### Unstructured Mesh with Coordinates

See `/workspace/config/unstructured_gmsh_example.config`

## Running Simulations

All examples are driven by configuration files:

```bash
# Run with default config location
./ex_config_driven config/my_simulation.config

# Run specific examples
./ex_buckley_leverett_2d config/buckley_leverett_2d.config
./ex_coupled_reservoir_2d config/coupled_reservoir_2d.config
./spe1 config/spe1_benchmark.config

# Run in parallel
mpirun -np 4 ./ex_config_driven config/large_simulation.config
```

## See Also

- [Unstructured Meshes](UNSTRUCTURED_MESHES.md)
- [Coordinate Systems](COORDINATE_SYSTEMS.md)
- [Physics Models](PHYSICS.md)
