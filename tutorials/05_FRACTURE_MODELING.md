# Tutorial 05: Fracture and Hydraulic Fracturing Modeling

Model natural fracture networks and hydraulic fracturing operations in FSRM.

## Table of Contents

1. [Fracture Modeling Overview](#fracture-modeling-overview)
2. [Natural Fracture Networks](#natural-fracture-networks)
3. [Hydraulic Fracture Models](#hydraulic-fracture-models)
4. [Fracture Propagation](#fracture-propagation)
5. [Proppant Transport](#proppant-transport)
6. [Leak-off Modeling](#leak-off-modeling)
7. [Multi-Stage Fracturing](#multi-stage-fracturing)
8. [Practical Examples](#practical-examples)

---

## Fracture Modeling Overview

FSRM provides multiple approaches to fracture modeling:

```
┌──────────────────────────────────────────────────────────────┐
│               Fracture Modeling Approaches                   │
├──────────────────────────────────────────────────────────────┤
│                                                               │
│  1. EXPLICIT FRACTURES (DFN)                                 │
│     ├─ Individual fractures in mesh                          │
│     ├─ Aperture and permeability per fracture               │
│     └─ Best for small number of dominant fractures          │
│                                                               │
│  2. DUAL POROSITY/PERMEABILITY                               │
│     ├─ Continuum approach                                    │
│     ├─ Matrix-fracture transfer                              │
│     └─ Best for densely fractured media                     │
│                                                               │
│  3. HYDRAULIC FRACTURING                                     │
│     ├─ PKN, KGD, P3D models                                 │
│     ├─ Propagation mechanics                                 │
│     └─ Proppant transport                                   │
│                                                               │
│  4. COHESIVE ZONE MODEL                                      │
│     ├─ LEFM-based propagation                               │
│     └─ Stress-dependent behavior                             │
│                                                               │
└──────────────────────────────────────────────────────────────┘
```

### Enabling Fractures

```ini
[SIMULATION]
enable_fractures = true
fracture_model = DFN          # DFN, DUAL_POROSITY, HYDRAULIC
```

---

## Natural Fracture Networks

### Single Fracture Definition

```ini
[FRACTURE1]
type = NATURAL                 # NATURAL or HYDRAULIC

# Location: center point (x, y, z) and normal vector (nx, ny, nz)
location = 500.0, 500.0, 50.0, 0.707, 0.707, 0.0

# Hydraulic properties
aperture = 0.0001              # 0.1 mm aperture
permeability = 1.0e-12         # m² (~ 1000 D)

# Mechanical properties
stiffness_normal = 1.0e10      # Pa/m
stiffness_shear = 1.0e9        # Pa/m

# Geometry
half_length = 100.0            # meters
height = 50.0                  # meters

# Propagation
enable_propagation = false     # Static fracture
```

### Fracture Orientation Options

```ini
[FRACTURE1]
# Method 1: Normal vector
location = 500, 500, 50, 0, 0, 1     # Horizontal fracture

# Method 2: Strike and dip
orientation_method = STRIKE_DIP
strike = 45.0                  # degrees from N
dip = 90.0                     # degrees from horizontal

# Method 3: Dip direction and dip
orientation_method = DIP_DIRECTION
dip_direction = 135.0          # degrees from N
dip = 60.0
```

### Stochastic Fracture Network

Generate a random fracture network:

```ini
[FRACTURES]
type = DFN_STOCHASTIC

# Number of fractures
num_fractures = 100

# Size distribution (power law)
size_distribution = POWER_LAW
min_length = 10.0              # meters
max_length = 500.0
exponent = 2.5                 # -D in N(>L) ∝ L^-D

# Orientation distribution
orientation_distribution = FISHER
mean_strike = 45.0             # degrees
mean_dip = 80.0
fisher_kappa = 20.0            # Concentration parameter

# Aperture model
aperture_model = LENGTH_CORRELATED
aperture_coefficient = 1.0e-5  # a = c × L

# Spatial distribution
spatial_distribution = UNIFORM # or CLUSTERED
region = 0, 1000, 0, 1000, 0, 100  # xmin, xmax, ymin, ymax, zmin, zmax

# Random seed (for reproducibility)
seed = 12345
```

### Dual Porosity Model

Continuum approach for densely fractured rock:

```ini
[SIMULATION]
fracture_model = DUAL_POROSITY

[ROCK]
# Matrix properties
porosity = 0.05                # Matrix porosity
permeability_x = 0.01          # Matrix perm (mD)

# Fracture properties  
fracture_porosity = 0.01       # Fracture porosity
fracture_permeability = 100.0  # Fracture perm (mD)

# Shape factor model
shape_factor_model = WARREN_ROOT  # WARREN_ROOT, KAZEMI, COATS
matrix_block_size = 1.0        # meters

# Transfer function
transfer_function = PSEUDO_STEADY_STATE  # or TRANSIENT
```

### Shape Factor Options

| Model | Formula | Best For |
|-------|---------|----------|
| Warren-Root | σ = 12/L² | Simple, fast |
| Kazemi | σ = 4(1/Lx² + 1/Ly² + 1/Lz²) | Anisotropic blocks |
| Coats | σ = 8/L² | History matching |

---

## Hydraulic Fracture Models

### PKN Model (Perkins-Kern-Nordgren)

Long, contained fractures:

```ini
[FRACTURE1]
type = HYDRAULIC
model = PKN                    # PKN, KGD, P3D, RADIAL

# Initial geometry
center = 500.0, 500.0, 50.0
orientation = 0.0, 1.0, 0.0    # Perpendicular to min stress
initial_length = 10.0          # meters
height = 50.0                  # Contained height

# Rock properties
youngs_modulus = 30.0e9
poisson_ratio = 0.25
fracture_toughness = 1.5e6     # Pa√m

# Fluid properties
injection_rate = 0.1           # m³/s
fluid_viscosity = 0.1          # Pa·s (slickwater = 0.001)

# Stress state
min_horizontal_stress = 30.0e6 # σhmin (Pa)
max_horizontal_stress = 35.0e6 # σHmax (Pa)
closure_pressure = 30.0e6      # Equal to σhmin
```

### KGD Model (Kristianovic-Geertsma-de Klerk)

Penny-shaped or short fractures:

```ini
[FRACTURE1]
type = HYDRAULIC
model = KGD

center = 500.0, 500.0, 50.0
orientation = 0.0, 1.0, 0.0
initial_length = 10.0

youngs_modulus = 30.0e9
poisson_ratio = 0.25
fracture_toughness = 1.5e6

injection_rate = 0.05
fluid_viscosity = 0.1
```

### P3D Model (Pseudo-3D)

Height growth with barriers:

```ini
[FRACTURE1]
type = HYDRAULIC
model = P3D

center = 500.0, 500.0, 50.0
orientation = 0.0, 1.0, 0.0
initial_length = 10.0
initial_height = 20.0

# Layered stress profile (enables height growth)
num_layers = 3
layer_depths = 30.0, 55.0, 80.0    # Bottom of each layer
layer_stresses = 32.0e6, 30.0e6, 33.0e6  # Closure stress per layer

# Barrier effect
barrier_stress_contrast = 3.0e6    # Stress jump at barriers

youngs_modulus = 30.0e9
poisson_ratio = 0.25
fracture_toughness = 1.5e6

injection_rate = 0.08
fluid_viscosity = 0.05
```

### Planar 3D Model

Full 3D fracture geometry:

```ini
[FRACTURE1]
type = HYDRAULIC
model = PLANAR_3D

center = 500.0, 500.0, 50.0
orientation = 0.0, 1.0, 0.0

# Mesh resolution for fracture
fracture_nx = 50               # Elements along length
fracture_nz = 30               # Elements along height

youngs_modulus = 30.0e9
poisson_ratio = 0.25
fracture_toughness = 1.5e6

injection_rate = 0.1
fluid_viscosity = 0.03

# Stress field (or use [IC] for stress)
stress_gradient = 0.0, 0.0, 25000.0  # Pa/m in z
```

---

## Fracture Propagation

### LEFM-Based Propagation

Linear Elastic Fracture Mechanics:

```ini
[FRACTURE1]
type = HYDRAULIC
enable_propagation = true

# LEFM parameters
propagation_criterion = KI_CRITICAL  # KI_CRITICAL, ENERGY, COHESIVE

# Mode I (opening) toughness
fracture_toughness = 1.5e6     # KIc (Pa√m)

# Optional: Mode II (shear) toughness
# shear_toughness = 2.0e6      # KIIc (Pa√m)

# Propagation direction
propagation_direction = MAX_HOOP_STRESS  # or TANGENT, ENERGY_RELEASE
```

### Cohesive Zone Model

For more accurate tip behavior:

```ini
[FRACTURE1]
enable_propagation = true
propagation_criterion = COHESIVE

# Cohesive parameters
cohesive_strength = 5.0e6      # Peak traction (Pa)
fracture_energy = 100.0        # Gc (J/m²)

# Traction-separation law
cohesive_law = LINEAR          # LINEAR, EXPONENTIAL, TRAPEZOIDAL

# Characteristic length
critical_opening = 4.0e-5      # δc (m), where Gc = σc × δc / 2
```

### Stress-Dependent Propagation

```ini
[FRACTURE1]
# Re-orientation based on stress field
enable_reorientation = true
reorientation_criterion = PRINCIPAL_STRESS

# Interaction with existing fractures
enable_fracture_interaction = true
interaction_distance = 10.0    # meters
```

---

## Proppant Transport

### Basic Proppant Model

```ini
[FRACTURE1]
type = HYDRAULIC
enable_proppant = true

# Proppant properties
proppant_density = 2650.0      # kg/m³
proppant_diameter = 0.0005     # 0.5 mm (20/40 mesh)
proppant_concentration = 0.2   # Initial volume fraction

# Transport physics
enable_settling = true
settling_model = STOKES        # STOKES, HINDERED
hindering_exponent = 4.65      # Richardson-Zaki

enable_bridging = true
bridging_aperture_ratio = 3.0  # Bridge when w < 3×d

enable_convection = true       # Advection with fluid
enable_diffusion = true
proppant_diffusivity = 1.0e-8  # m²/s
```

### Proppant Schedule

```ini
[FRACTURE1]
# Slickwater pad (no proppant)
stage_1_duration = 600.0       # 10 minutes
stage_1_rate = 0.15            # m³/s
stage_1_concentration = 0.0

# Ramp up proppant
stage_2_duration = 1200.0      # 20 minutes
stage_2_rate = 0.12
stage_2_concentration = 0.05   # 5% proppant

# Main proppant
stage_3_duration = 2400.0      # 40 minutes
stage_3_rate = 0.10
stage_3_concentration = 0.15   # 15% proppant

# Flush
stage_4_duration = 300.0       # 5 minutes
stage_4_rate = 0.15
stage_4_concentration = 0.0
```

### Propped Conductivity

```ini
[FRACTURE1]
# Post-closure conductivity
conductivity_model = PROPPANT_PACK

# Pack properties
proppant_pack_porosity = 0.38
proppant_pack_permeability = 1.0e-10  # m² (reference)

# Stress-dependent conductivity
enable_stress_dependency = true
# k/k0 = exp(-α × σn')
conductivity_decline_coefficient = 5.0e-8  # 1/Pa

# Embedment
enable_embedment = true
embedment_coefficient = 1.0e-5 # m/Pa
```

---

## Leak-off Modeling

### Carter Leak-off

Classical square-root-of-time model:

```ini
[FRACTURE1]
enable_leakoff = true
leakoff_model = CARTER

# Carter leak-off coefficient
leakoff_coefficient = 5.0e-5   # m/√s

# Spurt loss
spurt_loss = 0.001             # m³/m² instantaneous
```

### Pressure-Dependent Leak-off

```ini
[FRACTURE1]
leakoff_model = PRESSURE_DEPENDENT

# k × Δp / (μ × √(π × k × φ × ct × t))
matrix_permeability = 0.001    # mD (for calculation)
matrix_porosity = 0.05
matrix_compressibility = 5.0e-10

# Wall-building coefficient (filter cake)
filter_cake_permeability = 1.0e-18  # m²
filter_cake_thickness = 0.001  # m
```

### Dynamic Leak-off

```ini
[FRACTURE1]
leakoff_model = DYNAMIC

# Coupled with reservoir simulator
# Leak-off computed from pressure difference and local permeability
enable_reservoir_coupling = true
```

---

## Multi-Stage Fracturing

### Sequential Stages

```ini
[SIMULATION]
enable_fractures = true
multistage_fracturing = true

# Overall schedule
num_stages = 5
stage_spacing = 100.0          # meters between stages

[STAGE1]
time_start = 0.0
time_end = 3600.0              # 1 hour per stage
center = 200.0, 500.0, 50.0
injection_rate = 0.1
enable_proppant = true

[STAGE2]
time_start = 7200.0            # Start after Stage 1 + wait
time_end = 10800.0
center = 300.0, 500.0, 50.0
injection_rate = 0.1
enable_proppant = true

# ... repeat for more stages
```

### Stress Shadowing

```ini
[FRACTURES]
enable_stress_shadowing = true

# Method
shadowing_model = SNEDDON      # SNEDDON, DDM, FULL_COUPLING

# Parameters
shadow_distance = 200.0        # meters (influence range)

# Can affect:
# - Propagation direction (turning)
# - Net pressure increase
# - Asymmetric growth
```

### Simultaneous/Zipper Fracs

```ini
[MULTISTAGE]
fracturing_pattern = ZIPPER    # SEQUENTIAL, SIMULTANEOUS, ZIPPER

# For zipper frac (alternating wells)
well_1_stages = 1, 3, 5
well_2_stages = 2, 4, 6

# Timing
inter_stage_delay = 3600.0     # 1 hour between stages
```

---

## Practical Examples

### Shale Horizontal Well with Fractures

```ini
[SIMULATION]
name = shale_msfrac
fluid_model = SINGLE_COMPONENT
enable_geomechanics = true
enable_fractures = true
end_time = 157680000.0         # 5 years

[GRID]
nx = 100
ny = 50
nz = 10
Lx = 2000.0
Ly = 1000.0
Lz = 100.0

[ROCK]
constitutive_model = POROELASTIC
porosity = 0.08
permeability_x = 0.0005        # Ultra-low perm shale
permeability_y = 0.0005
permeability_z = 0.00005
youngs_modulus = 35.0e9
poisson_ratio = 0.25
biot_coefficient = 0.8

[FLUID]
type = SINGLE_PHASE
density = 200.0                # Gas
viscosity = 2.0e-5
compressibility = 5.0e-8

[WELL1]
name = HORZ-1
type = PRODUCER
well_geometry = HORIZONTAL
heel_i = 10
heel_j = 25
heel_k = 5
toe_i = 90
toe_j = 25
toe_k = 5
control_mode = BHP
target_value = 5.0e6

# Create multiple fractures along well
[FRACTURE1]
type = HYDRAULIC
model = P3D
center = 200.0, 500.0, 50.0
enable_propagation = false     # Pre-created
half_length = 150.0
height = 80.0
conductivity = 5.0e-14         # m³

[FRACTURE2]
type = HYDRAULIC
model = P3D
center = 400.0, 500.0, 50.0
half_length = 150.0
height = 80.0
conductivity = 5.0e-14

[FRACTURE3]
type = HYDRAULIC
model = P3D
center = 600.0, 500.0, 50.0
half_length = 150.0
height = 80.0
conductivity = 5.0e-14

[FRACTURE4]
type = HYDRAULIC
model = P3D
center = 800.0, 500.0, 50.0
half_length = 150.0
height = 80.0
conductivity = 5.0e-14

[FRACTURE5]
type = HYDRAULIC
model = P3D
center = 1000.0, 500.0, 50.0
half_length = 150.0
height = 80.0
conductivity = 5.0e-14
```

### Naturally Fractured Carbonate

```ini
[SIMULATION]
name = nfc_waterflood
fluid_model = BLACK_OIL
enable_fractures = true
fracture_model = DUAL_POROSITY

[ROCK]
# Matrix
porosity = 0.08
permeability_x = 1.0           # mD

# Fracture system
fracture_porosity = 0.01
fracture_permeability = 500.0  # mD

# Shape factor
shape_factor_model = KAZEMI
matrix_block_size = 2.0        # 2m blocks

[WELL1]
name = PROD-1
type = PRODUCER
control_mode = RATE
rate_type = LIQUID
target_value = 0.01
# ... location

[WELL2]
name = INJ-1
type = INJECTOR
fluid = WATER
control_mode = RATE
target_value = 0.01
```

### Real-Time Frac Job Simulation

```ini
[SIMULATION]
name = realtime_frac
enable_fractures = true
mode = PUMPING                 # Real-time pumping simulation
end_time = 7200.0              # 2-hour job

[FRACTURE1]
type = HYDRAULIC
model = PLANAR_3D
enable_propagation = true

# Real-time input from pumping data
pumping_schedule_file = pump_data.csv
# Format: time, rate, concentration

# Output for real-time monitoring
output_frequency = 10.0        # Every 10 seconds
output_net_pressure = true
output_fracture_geometry = true
output_proppant_coverage = true
```

---

## Fracture Output and Visualization

### Output Options

```ini
[OUTPUT]
# Fracture-specific outputs
fracture_geometry = true       # VTK files for fractures
fracture_aperture = true
fracture_pressure = true
fracture_proppant = true
fracture_conductivity = true

# DFN outputs
dfn_network = true             # Full network visualization
dfn_connectivity = true        # Connection matrix
```

### Post-Processing

```python
import numpy as np
import matplotlib.pyplot as plt

# Load fracture data
data = np.loadtxt('output/fracture_1_geometry.csv', delimiter=',', skiprows=1)
time = data[:, 0]
length = data[:, 1]
height = data[:, 2]
width = data[:, 3] * 1000  # Convert to mm

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Length vs time
axes[0, 0].plot(time / 60, length)
axes[0, 0].set_xlabel('Time (min)')
axes[0, 0].set_ylabel('Half-Length (m)')
axes[0, 0].set_title('Fracture Length Growth')

# Width profile
axes[0, 1].plot(time / 60, width)
axes[0, 1].set_xlabel('Time (min)')
axes[0, 1].set_ylabel('Max Width (mm)')
axes[0, 1].set_title('Fracture Width')

# Height vs time
axes[1, 0].plot(time / 60, height)
axes[1, 0].set_xlabel('Time (min)')
axes[1, 0].set_ylabel('Height (m)')
axes[1, 0].set_title('Fracture Height')

# Net pressure (if available)
if data.shape[1] > 4:
    net_p = data[:, 4] / 1e6
    axes[1, 1].plot(time / 60, net_p)
    axes[1, 1].set_xlabel('Time (min)')
    axes[1, 1].set_ylabel('Net Pressure (MPa)')
    axes[1, 1].set_title('Treatment Pressure')

plt.tight_layout()
plt.savefig('fracture_analysis.png', dpi=150)
```

---

**Previous**: [← Well Modeling](04_WELL_MODELING.md) | **Next**: [Fault Mechanics →](06_FAULT_MECHANICS.md)
