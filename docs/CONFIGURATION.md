# FSRM Configuration Reference

Complete reference for all configuration file options.

---

## File Format

FSRM uses INI-style configuration files:

```ini
# Comments start with # or ;
[SECTION_NAME]
key = value
numeric = 1.5e6           # Scientific notation
array = 1.0, 2.0, 3.0     # Comma-separated
string = myfile.txt       # Plain strings
bool = true               # true/false
```

---

## [SIMULATION] Section

Core simulation parameters.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `name` | string | "simulation" | Simulation identifier |
| `start_time` | float | 0.0 | Start time (s) |
| `end_time` | float | 1.0 | End time (s) |
| `dt_initial` | float | 0.001 | Initial timestep (s) |
| `dt_min` | float | 1e-10 | Minimum timestep (s) |
| `dt_max` | float | 1e6 | Maximum timestep (s) |
| `adaptive_dt` | bool | true | Enable adaptive timestepping |

### Physics Flags

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `fluid_model` | enum | SINGLE_COMPONENT | Fluid model type |
| `solid_model` | enum | ELASTIC | Solid model type |
| `enable_geomechanics` | bool | false | Enable geomechanics coupling |
| `enable_thermal` | bool | false | Enable thermal effects |
| `enable_fractures` | bool | false | Enable fracture modeling |
| `enable_faults` | bool | false | Enable fault mechanics |
| `enable_elastodynamics` | bool | false | Enable wave propagation |
| `enable_poroelastodynamics` | bool | false | Enable coupled poro-wave |

### Solver Settings

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `rtol` | float | 1e-6 | Relative tolerance |
| `atol` | float | 1e-8 | Absolute tolerance |
| `max_nonlinear_iterations` | int | 50 | Max SNES iterations |
| `max_linear_iterations` | int | 1000 | Max KSP iterations |
| `solver_type` | string | "gmres" | Linear solver (gmres, cg) |
| `preconditioner` | string | "ilu" | Preconditioner (ilu, jacobi, asm, hypre) |

### GPU Settings

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `use_gpu` | bool | false | Enable GPU acceleration |
| `gpu_mode` | enum | CPU_FALLBACK | GPU_ONLY, CPU_FALLBACK, AUTO |
| `gpu_device_id` | int | 0 | GPU device index |
| `gpu_memory_fraction` | float | 0.8 | Max GPU memory fraction |

### Fluid Model Types

| Value | Description |
|-------|-------------|
| `SINGLE_COMPONENT` | Single-phase compressible |
| `BLACK_OIL` | Three-phase black oil |
| `COMPOSITIONAL` | Multi-component with EOS |
| `BRINE` | Saline water |
| `CO2` | Supercritical CO2 |

### Solid Model Types

| Value | Description |
|-------|-------------|
| `ELASTIC` | Linear elastic |
| `VISCOELASTIC` | Viscoelastic (Maxwell/Kelvin-Voigt/SLS) |
| `POROELASTIC` | Biot poroelastic |
| `ELASTOPLASTIC` | With plasticity (Mohr-Coulomb/Drucker-Prager) |
| `VTI` | Vertical transverse isotropy |
| `HTI` | Horizontal transverse isotropy |

---

## [GRID] Section

Domain and mesh configuration.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `nx` | int | 10 | Grid cells in X |
| `ny` | int | 10 | Grid cells in Y |
| `nz` | int | 1 | Grid cells in Z |
| `Lx` | float | 1000.0 | Domain length X (m) |
| `Ly` | float | 1000.0 | Domain length Y (m) |
| `Lz` | float | 100.0 | Domain length Z (m) |
| `origin_x` | float | 0.0 | Domain origin X (m) |
| `origin_y` | float | 0.0 | Domain origin Y (m) |
| `origin_z` | float | 0.0 | Domain origin Z (m) |
| `grid_type` | enum | CARTESIAN | CARTESIAN, CORNER_POINT |
| `refinement_level` | int | 0 | Mesh refinement (0=none) |

---

## [ROCK] Section

Default rock/material properties. Use `[ROCK1]`, `[ROCK2]`, etc. for heterogeneous materials.

### Basic Properties

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `porosity` | float | 0.2 | Porosity (-) |
| `permeability_x` | float | 100.0 | Permeability X (mD) |
| `permeability_y` | float | 100.0 | Permeability Y (mD) |
| `permeability_z` | float | 10.0 | Permeability Z (mD) |
| `density` | float | 2650.0 | Grain density (kg/m³) |
| `compressibility` | float | 1e-9 | Rock compressibility (1/Pa) |

### Elastic Properties

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `constitutive_model` | enum | LINEAR_ELASTIC | Constitutive model |
| `youngs_modulus` | float | 20e9 | Young's modulus (Pa) |
| `poisson_ratio` | float | 0.25 | Poisson's ratio (-) |
| `bulk_modulus` | float | - | Bulk modulus (Pa) |
| `shear_modulus` | float | - | Shear modulus (Pa) |

### Poroelastic Properties

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `biot_coefficient` | float | 1.0 | Biot coefficient (-) |
| `biot_modulus` | float | 1e10 | Biot modulus (Pa) |
| `undrained_poisson` | float | - | Undrained Poisson ratio |
| `skempton` | float | - | Skempton coefficient |

### Viscoelastic Properties

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `viscoelastic_type` | enum | MAXWELL | MAXWELL, KELVIN_VOIGT, SLS |
| `viscosity` | float | 1e18 | Viscosity (Pa·s) |
| `relaxation_time` | float | - | Relaxation time (s) |

### Plasticity Properties

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `failure_criterion` | enum | MOHR_COULOMB | MOHR_COULOMB, DRUCKER_PRAGER |
| `cohesion` | float | 1e6 | Cohesion (Pa) |
| `friction_angle` | float | 30.0 | Friction angle (degrees) |
| `dilation_angle` | float | 0.0 | Dilation angle (degrees) |
| `tensile_strength` | float | 0.0 | Tensile strength (Pa) |

### Anisotropic Properties (VTI/HTI)

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `c11` | float | - | Stiffness C11 (Pa) |
| `c12` | float | - | Stiffness C12 (Pa) |
| `c13` | float | - | Stiffness C13 (Pa) |
| `c33` | float | - | Stiffness C33 (Pa) |
| `c44` | float | - | Stiffness C44 (Pa) |
| `thomsen_epsilon` | float | - | Thomsen epsilon (-) |
| `thomsen_delta` | float | - | Thomsen delta (-) |
| `thomsen_gamma` | float | - | Thomsen gamma (-) |

### Thermal Properties

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `thermal_conductivity` | float | 2.5 | Conductivity (W/m·K) |
| `specific_heat` | float | 800.0 | Specific heat (J/kg·K) |
| `thermal_expansion` | float | 1e-5 | Expansion coefficient (1/K) |

### Permeability Model

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `permeability_model` | enum | CONSTANT | CONSTANT, KOZENY_CARMAN, EXPONENTIAL, CUBIC_LAW |
| `permeability_reference` | float | 100.0 | Reference permeability (mD) |
| `permeability_exponent` | float | 3.0 | Exponential factor |

---

## [FLUID] Section

Fluid model configuration.

### Common Properties

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `type` | enum | SINGLE_PHASE | Fluid type |
| `density` | float | 1000.0 | Reference density (kg/m³) |
| `viscosity` | float | 0.001 | Viscosity (Pa·s) |
| `compressibility` | float | 4.5e-10 | Compressibility (1/Pa) |
| `reference_pressure` | float | 1e5 | Reference pressure (Pa) |

### Black Oil Properties

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `oil_density_std` | float | 850.0 | Oil density at std (kg/m³) |
| `gas_density_std` | float | 0.9 | Gas density at std (kg/m³) |
| `water_density_std` | float | 1000.0 | Water density at std (kg/m³) |
| `solution_gor` | float | 100.0 | Solution GOR (scf/stb) |
| `bubble_point` | float | 15e6 | Bubble point (Pa) |
| `pvt_correlation` | enum | STANDING | PVT correlation |

### PVT Correlations

| Value | Description |
|-------|-------------|
| `STANDING` | Standing correlation |
| `VASQUEZ_BEGGS` | Vazquez-Beggs correlation |
| `GLASO` | Glasø correlation |
| `AL_MARHOUN` | Al-Marhoun correlation |

### Compositional Properties

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `num_components` | int | 3 | Number of components |
| `component_names` | array | - | Component names |
| `critical_temperatures` | array | - | Tc values (K) |
| `critical_pressures` | array | - | Pc values (Pa) |
| `acentric_factors` | array | - | Acentric factors |
| `molar_weights` | array | - | Molar weights (kg/kmol) |
| `eos_type` | enum | PENG_ROBINSON | Equation of state |

### EOS Types

| Value | Description |
|-------|-------------|
| `PENG_ROBINSON` | Peng-Robinson cubic EOS |
| `SRK` | Soave-Redlich-Kwong EOS |
| `VAN_DER_WAALS` | van der Waals EOS |
| `IDEAL` | Ideal gas law |

### Brine Properties

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `salinity` | float | 0.0 | Salinity (ppm or kg/kg) |
| `salt_type` | string | "NaCl" | Salt composition |

### CO2 Properties

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `co2_model` | enum | SPAN_WAGNER | CO2 property model |

---

## [WELLn] Sections

Well definitions. Use `[WELL1]`, `[WELL2]`, etc.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `name` | string | - | Well identifier (required) |
| `type` | enum | PRODUCER | Well type |
| `i` | int | - | Grid index I (required) |
| `j` | int | - | Grid index J (required) |
| `k` | int | - | Grid index K (required) |
| `control_mode` | enum | BHP | Control type |
| `target_value` | float | - | Target (rate or pressure) |
| `max_rate` | float | - | Max rate constraint (m³/s) |
| `min_bhp` | float | - | Min BHP constraint (Pa) |
| `max_bhp` | float | - | Max BHP constraint (Pa) |
| `diameter` | float | 0.1 | Wellbore diameter (m) |
| `skin` | float | 0.0 | Skin factor (-) |
| `wellbore_storage` | float | 0.0 | Wellbore storage (m³/Pa) |
| `fluid` | enum | OIL | Injection fluid type |
| `temperature` | float | - | Injection temperature (K) |

### Well Types

| Value | Description |
|-------|-------------|
| `PRODUCER` | Production well |
| `INJECTOR` | Injection well |
| `OBSERVATION` | Monitoring only |

### Control Modes

| Value | Description |
|-------|-------------|
| `RATE` | Constant surface rate |
| `BHP` | Constant bottomhole pressure |
| `THP` | Constant tubing head pressure |

---

## [FRACTUREn] Sections

Fracture definitions.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `type` | enum | NATURAL | NATURAL, HYDRAULIC |
| `location` | array | - | x, y, z, nx, ny, nz |
| `aperture` | float | 1e-4 | Aperture (m) |
| `permeability` | float | 1e-12 | Permeability (m²) |
| `half_length` | float | - | Half-length (m) |
| `height` | float | - | Height (m) |
| `conductivity` | float | - | Conductivity (m³) |
| `enable_propagation` | bool | false | Enable growth |
| `toughness` | float | 1e6 | Fracture toughness (Pa·√m) |
| `energy` | float | 100 | Fracture energy (J/m²) |
| `enable_proppant` | bool | false | Track proppant |
| `proppant_density` | float | 2650 | Proppant density (kg/m³) |
| `proppant_diameter` | float | 0.001 | Proppant diameter (m) |

---

## [FAULTn] Sections

Fault definitions.

### Geometry

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `name` | string | - | Fault identifier |
| `x` | float | - | Center X (m) |
| `y` | float | - | Center Y (m) |
| `z` | float | - | Center depth (m) |
| `strike` | float | 0.0 | Strike (degrees N) |
| `dip` | float | 90.0 | Dip (degrees) |
| `length` | float | 1000.0 | Along-strike length (m) |
| `width` | float | 500.0 | Down-dip width (m) |
| `rake` | float | 0.0 | Slip direction (degrees) |

### Friction

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `friction_law` | enum | COULOMB | Friction model |
| `static_friction` | float | 0.6 | Static coefficient |
| `dynamic_friction` | float | 0.4 | Dynamic coefficient |
| `cohesion` | float | 0.0 | Cohesion (Pa) |

### Rate-State Friction

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `rate_state_a` | float | 0.01 | Direct effect a |
| `rate_state_b` | float | 0.015 | Evolution effect b |
| `rate_state_dc` | float | 1e-5 | Critical slip Dc (m) |
| `rate_state_v0` | float | 1e-6 | Reference velocity (m/s) |
| `rate_state_f0` | float | 0.6 | Reference friction |

### Friction Laws

| Value | Description |
|-------|-------------|
| `COULOMB` | Static/dynamic Coulomb |
| `RATE_STATE_AGING` | RS with aging law |
| `RATE_STATE_SLIP` | RS with slip law |
| `SLIP_WEAKENING` | Linear slip weakening |

---

## [BCn] Sections

Boundary conditions.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `type` | enum | DIRICHLET | BC type |
| `field` | enum | - | Field to apply (required) |
| `location` | enum | - | Face/region (required) |
| `value` | float | - | BC value |
| `time_function` | string | - | Time-dependent function |

### BC Types

| Value | Description |
|-------|-------------|
| `DIRICHLET` | Fixed value |
| `NEUMANN` | Fixed flux |
| `ROBIN` | Mixed type |
| `PERIODIC` | Periodic |

### Fields

| Value | Description |
|-------|-------------|
| `PRESSURE` | Pore pressure |
| `TEMPERATURE` | Temperature |
| `DISPLACEMENT` | Displacement |
| `SATURATION` | Saturation |
| `CONCENTRATION` | Species concentration |

### Locations

| Value | Description |
|-------|-------------|
| `XMIN` | Left face (X=0) |
| `XMAX` | Right face (X=Lx) |
| `YMIN` | Front face (Y=0) |
| `YMAX` | Back face (Y=Ly) |
| `ZMIN` | Bottom (Z=0) |
| `ZMAX` | Top (Z=Lz) |
| `ALL` | All boundaries |

---

## [ICn] Sections

Initial conditions.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `field` | enum | - | Field (required) |
| `distribution` | enum | UNIFORM | Distribution type |
| `value` | float | - | Base/uniform value |
| `gradient` | array | - | Gradient (dx, dy, dz) |
| `file` | string | - | Data file path |

### Distribution Types

| Value | Description |
|-------|-------------|
| `UNIFORM` | Constant value |
| `GRADIENT` | Linear gradient |
| `FUNCTION` | Analytic function |
| `FILE` | From external file |

---

## [DYNAMICS] Section

Wave propagation and dynamics settings.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `enable` | bool | false | Enable dynamics |
| `wave_type` | enum | ELASTIC | Wave physics type |
| `cfl_number` | float | 0.5 | CFL stability factor |
| `absorbing_layers` | int | 10 | PML thickness (cells) |
| `absorbing_strength` | float | 100.0 | PML damping coefficient |
| `source_type` | enum | MOMENT_TENSOR | Source mechanism |
| `source_location` | array | - | Source x, y, z (m) |
| `source_time_function` | enum | RICKER | Source wavelet |
| `source_frequency` | float | 10.0 | Peak frequency (Hz) |
| `source_amplitude` | float | 1.0 | Source amplitude |
| `dynamic_permeability` | bool | false | Perm changes with wave |
| `permeability_model` | enum | - | Dynamic perm model |

### Wave Types

| Value | Description |
|-------|-------------|
| `ELASTIC` | Elastic wave propagation |
| `POROELASTIC` | Biot's poroelastic waves |
| `ACOUSTIC` | Acoustic (pressure only) |

### Source Types

| Value | Description |
|-------|-------------|
| `MOMENT_TENSOR` | Full moment tensor |
| `POINT_FORCE` | Point force |
| `EXPLOSIVE` | Isotropic expansion |
| `DOUBLE_COUPLE` | Shear dislocation |

---

## [SEISMICITY] Section

Fault mechanics and induced seismicity.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `enable` | bool | false | Enable seismicity |
| `friction_law` | enum | COULOMB | Default friction |
| `nucleation_size` | float | 1.0 | Critical patch (m) |
| `seismic_slip_rate` | float | 1e-3 | Seismic threshold (m/s) |
| `b_value` | float | 1.0 | Gutenberg-Richter b |
| `aftershocks` | bool | true | Enable aftershocks |
| `stress_transfer` | bool | true | Coulomb stress changes |
| `max_magnitude` | float | 6.0 | Maximum allowed M |
| `catalog_file` | string | - | Output catalog path |

---

## [OUTPUT] Section

Output configuration.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `format` | enum | VTK | Output format |
| `path` | string | "output" | Output directory |
| `frequency` | int | 1 | Output every N steps |
| `time_interval` | float | - | Output time interval (s) |

### Field Output Flags

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `pressure` | bool | true | Output pressure |
| `displacement` | bool | false | Output displacement |
| `velocity` | bool | false | Output velocity |
| `stress` | bool | false | Output stress tensor |
| `strain` | bool | false | Output strain tensor |
| `permeability` | bool | false | Output permeability |
| `temperature` | bool | false | Output temperature |
| `saturation` | bool | false | Output saturation |
| `fault_slip` | bool | false | Output fault slip |
| `seismic_catalog` | bool | false | Write event catalog |

### Output Formats

| Value | Description |
|-------|-------------|
| `VTK` | VTK XML format (.vtu) |
| `HDF5` | HDF5 format |
| `ECLIPSE` | Eclipse restart format |
| `ASCII` | Plain text |

---

## Including Other Files

Configuration files can include other files:

```ini
# Include another config
@include "base_settings.config"

# Override specific values
[SIMULATION]
end_time = 86400.0
```

---

## Environment Variables

Reference environment variables:

```ini
[OUTPUT]
path = ${HOME}/simulations/output
```

---

## Generating Templates

Generate a complete template with all options:

```bash
fsrm -generate_config complete_template.config
```

This creates a documented template with default values.
