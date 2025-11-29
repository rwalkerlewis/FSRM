# Tutorial 06: Fault Mechanics and Induced Seismicity

Model fault behavior, friction laws, slip dynamics, and induced seismicity in FSRM.

## Table of Contents

1. [Fault Mechanics Overview](#fault-mechanics-overview)
2. [Fault Geometry Definition](#fault-geometry-definition)
3. [Friction Laws](#friction-laws)
4. [Stress Analysis on Faults](#stress-analysis-on-faults)
5. [Induced Seismicity Modeling](#induced-seismicity-modeling)
6. [Coulomb Stress Transfer](#coulomb-stress-transfer)
7. [Split Node Method](#split-node-method)
8. [IMEX Time Integration](#imex-time-integration)
9. [Practical Examples](#practical-examples)

---

## Fault Mechanics Overview

FSRM provides comprehensive fault mechanics capabilities:

```
┌──────────────────────────────────────────────────────────────┐
│                 Fault Mechanics System                       │
├──────────────────────────────────────────────────────────────┤
│  FRICTION MODELS                                             │
│  ├─ Coulomb (static/dynamic)                                │
│  ├─ Rate-and-state (aging/slip laws)                        │
│  ├─ Slip weakening                                          │
│  └─ Thermal pressurization                                  │
│                                                               │
│  FAULT REPRESENTATION                                        │
│  ├─ Embedded interfaces                                     │
│  ├─ Split nodes (discontinuous displacement)                │
│  └─ Zero-thickness elements                                 │
│                                                               │
│  SEISMICITY                                                   │
│  ├─ Event detection and cataloging                          │
│  ├─ Magnitude estimation                                    │
│  ├─ Coulomb stress transfer                                 │
│  └─ Aftershock triggering                                   │
│                                                               │
│  TIME INTEGRATION                                            │
│  ├─ Implicit (quasi-static)                                 │
│  ├─ Explicit (dynamic)                                      │
│  └─ IMEX (automatic switching)                              │
└──────────────────────────────────────────────────────────────┘
```

### Enabling Fault Mechanics

```ini
[SIMULATION]
enable_faults = true
enable_geomechanics = true     # Required for stress coupling
```

---

## Fault Geometry Definition

### Basic Fault

```ini
[FAULT1]
name = MAIN_FAULT

# Center location (meters)
x = 5000.0                     # X coordinate
y = 0.0                        # Y coordinate
z = 3000.0                     # Depth (positive down)

# Orientation
strike = 0.0                   # Degrees from North (azimuth)
dip = 60.0                     # Degrees from horizontal

# Dimensions
length = 2000.0                # Along-strike length (m)
width = 1000.0                 # Down-dip width (m)

# Optional: rake for slip direction
rake = 0.0                     # 0 = pure strike-slip, 90 = thrust, -90 = normal
```

### Orientation Conventions

```
      N (0°)
       ↑
       │
W (270°)──┼──→ E (90°)
       │
       ↓
      S (180°)

Strike: measured clockwise from N
Dip: measured from horizontal (right-hand rule from strike)
Rake: measured on fault plane (0° = left-lateral, 180° = right-lateral)
```

### Fault Mesh Resolution

```ini
[FAULT1]
# Discretization for split nodes
num_strike_elements = 40       # Along strike
num_dip_elements = 20          # Along dip

# Or specify element size
element_size = 50.0            # meters
```

### Multiple Faults

```ini
[FAULT1]
name = FAULT_A
x = 5000.0
y = 0.0
z = 3000.0
strike = 0.0
dip = 60.0
length = 2000.0
width = 1000.0

[FAULT2]
name = FAULT_B
x = 3000.0
y = 1000.0
z = 2500.0
strike = 45.0
dip = 70.0
length = 1500.0
width = 800.0

# Fault intersection handling
enable_fault_intersection = true
intersection_method = SPLIT    # SPLIT, MERGE, IGNORE
```

---

## Friction Laws

### Coulomb Friction

Classic static/dynamic friction:

```ini
[FAULT1]
friction_law = COULOMB

# Coefficients
static_friction = 0.6          # μs
dynamic_friction = 0.4         # μd
cohesion = 0.0                 # Pa (often 0 for pre-existing faults)

# Slip weakening (optional)
slip_weakening_dc = 0.01       # Critical slip distance (m)
```

**Failure Criterion:**
```
τ ≤ c + μ × σn'
```

### Rate-and-State Friction

For realistic earthquake nucleation and propagation:

```ini
[FAULT1]
friction_law = RATE_STATE_AGING  # or RATE_STATE_SLIP

# Parameters (Dieterich-Ruina)
rate_state_a = 0.010           # Direct effect (velocity strengthening)
rate_state_b = 0.015           # Evolution effect
rate_state_dc = 1.0e-4         # Critical slip distance Dc (m)
rate_state_v0 = 1.0e-6         # Reference velocity V0 (m/s)
rate_state_f0 = 0.6            # Reference friction f0

# Initial state
initial_state = 1.0e6          # θ0 (seconds)
```

**Friction Law:**
```
f = f0 + a × ln(V/V0) + b × ln(V0×θ/Dc)
```

**State Evolution:**
- Aging law: `dθ/dt = 1 - V×θ/Dc`
- Slip law: `dθ/dt = -V×θ/Dc × ln(V×θ/Dc)`

**Stability:**
```
a - b < 0  →  Velocity weakening (unstable, can nucleate earthquakes)
a - b > 0  →  Velocity strengthening (stable, aseismic creep)
```

### Slip Weakening

Linear weakening with slip:

```ini
[FAULT1]
friction_law = SLIP_WEAKENING

static_friction = 0.6
dynamic_friction = 0.4
critical_slip = 0.5            # Dc (m)

# Optional: re-strengthening
enable_healing = true
healing_rate = 1.0e-8          # 1/s
```

### Thermal Pressurization

High slip rate effects:

```ini
[FAULT1]
friction_law = THERMAL_PRESSURIZATION

# Base friction
base_friction = 0.6

# Thermal parameters
shear_zone_width = 0.001       # m (1 mm)
thermal_diffusivity = 1.0e-6   # m²/s
hydraulic_diffusivity = 1.0e-5 # m²/s
specific_heat = 1000.0         # J/(kg·K)
thermal_expansion = 1.0e-4     # 1/K (fluid)

# Weakening parameters
flash_heating_threshold = 0.1  # m/s slip rate
```

---

## Stress Analysis on Faults

### Initial Stress State

Define the regional stress field:

```ini
[IC_STRESS]
# Method 1: Principal stresses
stress_type = PRINCIPAL
sigma_1 = -50.0e6              # Maximum principal (compression negative)
sigma_2 = -40.0e6              # Intermediate
sigma_3 = -30.0e6              # Minimum principal
sigma_1_azimuth = 90.0         # Azimuth of σ1 from N
sigma_1_plunge = 0.0           # Plunge of σ1 from horizontal

# Method 2: Stress tensor components
# stress_type = TENSOR
# sxx = -40.0e6
# syy = -50.0e6
# szz = -30.0e6
# sxy = 5.0e6
# sxz = 0.0
# syz = 0.0

# Method 3: Gradient with depth
# stress_type = ANDERSON
# tectonic_regime = STRIKE_SLIP  # NORMAL, STRIKE_SLIP, THRUST
# vertical_gradient = 25000.0    # Pa/m (≈ρg)
# stress_ratio = 0.8             # K0 = σh/σv
```

### Fault Stress Resolution

Stresses are automatically resolved onto the fault plane:

```
σn = n·σ·n        (normal stress)
τ = |t - (t·n)n|  (shear stress)
```

### Coulomb Failure Function

```ini
[SEISMICITY]
# CFF = τ - μ(σn - P) - c
# Positive CFF = failure imminent/occurred

cff_output = true
cff_threshold = 0.0            # Alert when CFF > threshold
```

---

## Induced Seismicity Modeling

### Seismicity Configuration

```ini
[SIMULATION]
enable_faults = true
enable_geomechanics = true
solid_model = POROELASTIC      # Pore pressure coupling

[SEISMICITY]
enable = true

# Event detection
nucleation_size = 10.0         # Critical patch size (m)
seismic_slip_rate = 1.0e-3     # Threshold slip rate (m/s)

# Magnitude estimation
shear_modulus = 30.0e9         # For moment calculation

# Statistics
b_value = 1.0                  # Gutenberg-Richter b-value

# Additional physics
enable_aftershocks = true
enable_stress_transfer = true

# Event catalog output
catalog_file = output/seismic_catalog.csv

# Maximum allowed magnitude (safety)
max_magnitude = 5.0            # Stop if exceeded
```

### Magnitude Calculation

```
M0 = μ × A × D        (Seismic moment)
Mw = (2/3) × log10(M0) - 6.07  (Moment magnitude)
```

Where:
- μ = shear modulus
- A = rupture area
- D = average slip

### Event Catalog Format

```csv
time,x,y,z,magnitude,moment,slip,area,stress_drop,type
1000.5,5000,0,3000,2.1,1.2e12,0.02,5000,1.5e6,earthquake
1500.2,5050,100,2900,1.8,5.0e11,0.01,2000,1.0e6,earthquake
...
```

---

## Coulomb Stress Transfer

Model how earthquakes affect stress on nearby faults:

### Configuration

```ini
[SEISMICITY]
enable_stress_transfer = true

# Method
stress_transfer_method = OKADA  # OKADA, DDM, FEM

# Parameters
poisson_ratio = 0.25           # For elastic half-space
coefficient_of_friction = 0.4  # For receiver fault CFF

# Cascade triggering
enable_cascade = true
cascade_time_window = 86400.0  # 1 day for aftershocks
```

### Stress Transfer Calculation

```ini
# ΔCFF = Δτ + μ × Δσn
# Positive ΔCFF promotes failure
```

### Aftershock Triggering

```ini
[SEISMICITY]
# Dieterich (1994) seismicity rate model
enable_rate_change = true
Aσ = 1.0e6                     # Effective stress parameter (Pa)
ta = 1.0e8                     # Aftershock duration (s)

# Omori-Utsu decay
omori_c = 100.0                # c-value (s)
omori_p = 1.1                  # p-value
```

---

## Split Node Method

For explicit displacement discontinuities across faults:

### Configuration

```ini
[FAULT1]
# Enable split nodes
use_split_nodes = true

# Method
split_node_method = PENALTY    # PENALTY, LAGRANGE, NITSCHE, MORTAR

# Penalty parameters (for PENALTY method)
penalty_normal = 1.0e12        # Pa/m
penalty_tangent = 1.0e10       # Pa/m

# Contact behavior
allow_separation = false       # Allow tensile opening?
symmetric_contact = true       # Symmetric formulation
```

### Traction Types

```ini
[FAULT1]
use_split_nodes = true

# Traction type
traction_type = FRICTION_DEPENDENT  # PRESCRIBED, FRICTION_DEPENDENT, COHESIVE_ZONE

# For PRESCRIBED:
# prescribed_traction_normal = -30.0e6  # Compression
# prescribed_traction_strike = 5.0e6    # Shear
# prescribed_traction_dip = 0.0

# For COHESIVE_ZONE:
# cohesive_strength = 5.0e6    # Peak traction
# critical_opening = 1.0e-4    # m
# critical_slip = 1.0e-3       # m
```

### Augmented Lagrangian

For robust contact enforcement:

```ini
[FAULT1]
split_node_method = AUGMENTED_LAGRANGIAN

aug_lag_r = 1.0e10             # Penalty-like parameter
aug_lag_max_iter = 10          # Inner iterations
aug_lag_tol = 1.0e-8           # Convergence tolerance
```

---

## IMEX Time Integration

Automatic switching between implicit (slow) and explicit (fast) time integration:

### Why IMEX?

- **Implicit**: Efficient for quasi-static loading (large timesteps)
- **Explicit**: Necessary for dynamic rupture (small timesteps for stability)
- **IMEX**: Automatically switches when earthquake nucleates

### Configuration

```ini
[IMEX]
enabled = true
initial_mode = IMPLICIT        # Start in quasi-static mode

# Implicit settings
implicit_dt_initial = 3600.0   # 1 hour
implicit_dt_min = 1.0          # 1 second
implicit_dt_max = 86400.0      # 1 day
implicit_method = BEULER       # Backward Euler

# Explicit settings
explicit_dt_initial = 1.0e-4   # 0.1 ms
explicit_dt_min = 1.0e-8
explicit_dt_max = 1.0e-2
cfl_factor = 0.5               # CFL stability factor
explicit_method = NEWMARK      # Newmark-beta

# Triggering (implicit → explicit)
trigger_type = COULOMB_FAILURE # STRESS, COULOMB, SLIP_RATE, VELOCITY
coulomb_threshold = 0.0        # Switch when CFF ≥ 0
slip_rate_threshold = 1.0e-3   # m/s

# Settling (explicit → implicit)
settling_type = COMBINED       # TIME, VELOCITY, ENERGY, COMBINED
min_dynamic_duration = 1.0     # Minimum time in explicit (s)
max_dynamic_duration = 60.0    # Maximum time (force switch back)
settling_velocity = 1.0e-6     # m/s threshold
settling_slip_rate = 1.0e-9    # m/s

# Smooth transition
smooth_transition = true
transition_ramp_time = 0.01    # s
```

### Transition Logging

```ini
[IMEX]
log_transitions = true
transition_log_file = imex_transitions.log
```

Log format:
```
Time: 1000.000 s | Transition: IMPLICIT -> EXPLICIT | Trigger: CFF=1.2e6 Pa
Time: 1002.345 s | Transition: EXPLICIT -> IMPLICIT | Settled: Vmax=1.2e-7 m/s
```

---

## Practical Examples

### Injection-Induced Seismicity

```ini
[SIMULATION]
name = induced_seismicity
fluid_model = SINGLE_COMPONENT
solid_model = POROELASTIC
enable_geomechanics = true
enable_faults = true
end_time = 31536000.0          # 1 year

[GRID]
nx = 100
ny = 100
nz = 50
Lx = 10000.0
Ly = 10000.0
Lz = 5000.0

[ROCK]
constitutive_model = POROELASTIC
porosity = 0.15
permeability_x = 100.0
youngs_modulus = 30.0e9
poisson_ratio = 0.25
biot_coefficient = 0.8

[FLUID]
type = SINGLE_PHASE
density = 1000.0
viscosity = 0.001
compressibility = 4.5e-10

# Injection well near fault
[WELL1]
name = INJ-1
type = INJECTOR
fluid = WATER
i = 50
j = 50
k = 25
control_mode = RATE
target_value = 0.05            # 50 L/s = 4300 m³/day
max_bhp = 40.0e6

# Pre-existing fault
[FAULT1]
name = FAULT_A
x = 5500.0                     # 500m from well
y = 5000.0
z = 2500.0
strike = 45.0
dip = 70.0
length = 3000.0
width = 2000.0

friction_law = RATE_STATE_AGING
rate_state_a = 0.005
rate_state_b = 0.008           # b > a: unstable
rate_state_dc = 5.0e-5
rate_state_v0 = 1.0e-9
rate_state_f0 = 0.6

use_split_nodes = true
split_node_method = PENALTY

# Initial stress (critically stressed)
[IC_STRESS]
stress_type = ANDERSON
tectonic_regime = STRIKE_SLIP
vertical_gradient = 25000.0
stress_ratio = 0.8

# Seismicity monitoring
[SEISMICITY]
enable = true
nucleation_size = 5.0
seismic_slip_rate = 1.0e-3
b_value = 1.0
enable_stress_transfer = true
catalog_file = output/induced_catalog.csv
max_magnitude = 4.0

# IMEX for dynamic events
[IMEX]
enabled = true
initial_mode = IMPLICIT
trigger_type = COULOMB_FAILURE
coulomb_threshold = 0.0

[OUTPUT]
pressure = true
displacement = true
stress = true
fault_slip = true
seismic_catalog = true
```

### SCEC Benchmark (TPV3)

Standard earthquake rupture benchmark:

```ini
[SIMULATION]
name = scec_tpv3
solid_model = ELASTIC
enable_faults = true
enable_elastodynamics = true
end_time = 15.0                # 15 seconds of rupture

[GRID]
nx = 300
ny = 300
nz = 150
Lx = 30000.0
Ly = 30000.0
Lz = 15000.0

[ROCK]
youngs_modulus = 80.0e9        # 32 GPa μ, ν=0.25
poisson_ratio = 0.25
density = 2670.0

[FAULT1]
name = VERTICAL_STRIKE_SLIP
x = 15000.0
y = 15000.0
z = 7500.0
strike = 90.0                  # E-W striking
dip = 90.0                     # Vertical
length = 30000.0
width = 15000.0

friction_law = SLIP_WEAKENING
static_friction = 0.677        # τ0/σn0 = 70 MPa / 120 MPa
dynamic_friction = 0.525
critical_slip = 0.4            # Dc = 0.4 m

use_split_nodes = true

# Nucleation zone (forced rupture start)
nucleation_enabled = true
nucleation_x = 0.0             # Along strike from center
nucleation_z = 7500.0          # At depth
nucleation_radius = 3000.0     # 3 km patch
nucleation_traction = 81.6e6   # Force to start rupture

# Initial stress
[IC_STRESS]
stress_type = TENSOR
sxx = -120.0e6                 # Normal stress on fault
syy = 0.0
szz = -120.0e6
sxy = 70.0e6                   # Shear stress (right-lateral)

[DYNAMICS]
enable = true
cfl_number = 0.5
absorbing_layers = 20

[OUTPUT]
frequency = 100
displacement = true
velocity = true
fault_slip = true
```

### Traffic Light Protocol Implementation

```ini
# Monitoring with automatic response
[SEISMICITY]
enable = true

# Traffic light thresholds
tlp_enabled = true
tlp_green_max = 1.0            # M < 1.0: Green, continue
tlp_yellow_max = 2.0           # 1.0 ≤ M < 2.0: Yellow, reduce rate
tlp_red = 3.0                  # M ≥ 3.0: Red, stop injection

# Actions
tlp_yellow_rate_factor = 0.5   # Reduce to 50%
tlp_red_action = SHUTIN        # SHUTIN, REDUCE, ALERT

# Monitoring period
tlp_evaluation_window = 86400.0  # 1 day rolling window
```

---

## Output and Analysis

### Fault Output

```ini
[OUTPUT]
fault_slip = true              # Cumulative slip
fault_slip_rate = true         # Instantaneous slip rate
fault_traction = true          # Normal and shear tractions
fault_state = true             # Rate-state θ
fault_cff = true               # Coulomb failure function
```

### Analysis Scripts

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load seismic catalog
catalog = pd.read_csv('output/seismic_catalog.csv')

# Magnitude-time plot
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Time series
axes[0, 0].stem(catalog['time'] / 86400, catalog['magnitude'], 
                basefmt=' ', linefmt='b-', markerfmt='bo')
axes[0, 0].set_xlabel('Time (days)')
axes[0, 0].set_ylabel('Magnitude')
axes[0, 0].set_title('Seismicity Time Series')

# Gutenberg-Richter
mags = catalog['magnitude']
bins = np.arange(0, mags.max() + 0.5, 0.5)
counts, edges = np.histogram(mags, bins=bins)
cum_counts = np.cumsum(counts[::-1])[::-1]
mag_centers = (edges[:-1] + edges[1:]) / 2

axes[0, 1].semilogy(mag_centers[counts > 0], cum_counts[counts > 0], 'ko-')
axes[0, 1].set_xlabel('Magnitude')
axes[0, 1].set_ylabel('N(≥M)')
axes[0, 1].set_title(f'Gutenberg-Richter (b ≈ {1.0:.2f})')

# Map view
axes[1, 0].scatter(catalog['x'] / 1000, catalog['y'] / 1000, 
                   c=catalog['time'] / 86400, s=catalog['magnitude']**2 * 10,
                   alpha=0.7, cmap='viridis')
axes[1, 0].set_xlabel('X (km)')
axes[1, 0].set_ylabel('Y (km)')
axes[1, 0].set_title('Event Locations')
plt.colorbar(axes[1, 0].collections[0], ax=axes[1, 0], label='Time (days)')

# Depth distribution
axes[1, 1].hist(catalog['z'] / 1000, bins=20, orientation='horizontal')
axes[1, 1].set_xlabel('Count')
axes[1, 1].set_ylabel('Depth (km)')
axes[1, 1].invert_yaxis()
axes[1, 1].set_title('Depth Distribution')

plt.tight_layout()
plt.savefig('seismicity_analysis.png', dpi=150)
```

---

**Previous**: [← Fracture Modeling](05_FRACTURE_MODELING.md) | **Next**: [GPU Acceleration →](07_GPU_ACCELERATION.md)
