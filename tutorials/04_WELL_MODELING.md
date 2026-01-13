# Tutorial 04: Well Modeling

Learn how to configure production wells, injection wells, and advanced well types in FSRM.

## Table of Contents

1. [Well Model Overview](#well-model-overview)
2. [Basic Well Configuration](#basic-well-configuration)
3. [Well Control Modes](#well-control-modes)
4. [Multi-Phase Wells](#multi-phase-wells)
5. [Advanced Well Types](#advanced-well-types)
6. [Well Performance Analysis](#well-performance-analysis)
7. [Practical Examples](#practical-examples)

---

## Well Model Overview

FSRM supports various well types and configurations:

```
┌──────────────────────────────────────────────────────────────┐
│                    Well Modeling System                       │
├──────────────────────────────────────────────────────────────┤
│  WELL TYPES                                                   │
│  ├─ Vertical wells                                           │
│  ├─ Deviated/directional wells                               │
│  ├─ Horizontal wells                                         │
│  └─ Multilateral wells                                       │
│                                                               │
│  CONTROL MODES                                                │
│  ├─ Rate control (liquid, gas, water)                        │
│  ├─ Bottomhole pressure (BHP)                                │
│  └─ Tubing head pressure (THP)                               │
│                                                               │
│  PHYSICS                                                      │
│  ├─ Peaceman well index                                      │
│  ├─ Wellbore hydraulics                                      │
│  ├─ Inflow Performance Relationships (IPR)                   │
│  └─ Multi-phase flow in wellbore                             │
└──────────────────────────────────────────────────────────────┘
```

---

## Basic Well Configuration

### Vertical Producer

```ini
[WELL1]
name = PROD-1                  # Unique identifier
type = PRODUCER                # PRODUCER or INJECTOR

# Location (grid indices, 0-based)
i = 10                         # X index
j = 10                         # Y index
k = 5                          # Z index (completion layer)

# Control
control_mode = RATE            # RATE, BHP, or THP
target_value = 0.01            # m³/s for RATE, Pa for BHP

# Well properties
diameter = 0.2                 # Wellbore diameter (m)
skin = 0.0                     # Skin factor (-2 to +20 typical)
```

### Vertical Injector

```ini
[WELL2]
name = INJ-1
type = INJECTOR

i = 2
j = 2
k = 5

control_mode = RATE
target_value = 0.005           # 0.005 m³/s injection rate

# Injection fluid
fluid = WATER                  # WATER, GAS, STEAM, CO2

# Injection temperature (for thermal simulations)
temperature = 300.0            # K

# Pressure constraint
max_bhp = 40.0e6               # Don't exceed 40 MPa
```

### Observation Well

```ini
[WELL3]
name = OBS-1
type = OBSERVATION             # No flow, just monitoring

i = 15
j = 15
k = 5
```

---

## Well Control Modes

### Rate Control

Control the volumetric flow rate:

```ini
[WELL1]
control_mode = RATE
target_value = 0.01            # m³/s (surface conditions)

# For multi-phase, specify which rate
rate_type = LIQUID             # LIQUID, OIL, GAS, WATER
```

**Rate Type Options:**
| Type | Description |
|------|-------------|
| `LIQUID` | Total liquid rate (oil + water) |
| `OIL` | Oil rate at surface conditions |
| `GAS` | Gas rate at surface conditions |
| `WATER` | Water rate |
| `RESERVOIR_VOIDAGE` | Reservoir volume rate |

### BHP Control

Control bottomhole pressure:

```ini
[WELL1]
control_mode = BHP
target_value = 15.0e6          # 15 MPa bottomhole

# Rate constraints
max_rate = 0.05                # Max rate if BHP allows
```

### THP Control

Control tubing head pressure (requires wellbore model):

```ini
[WELL1]
control_mode = THP
target_value = 5.0e6           # 5 MPa at wellhead

# Enable wellbore hydraulics
enable_wellbore_model = true
wellbore_roughness = 1.5e-5    # Pipe roughness (m)
tubing_id = 0.1                # Inner diameter (m)
```

### Automatic Switching

Wells can switch control modes when limits are hit:

```ini
[WELL1]
control_mode = RATE
target_value = 0.02            # Target 0.02 m³/s

# Switch to BHP if rate can't be achieved
min_bhp = 5.0e6                # Switch when BHP drops to 5 MPa

# Or for injectors
max_bhp = 50.0e6               # Switch when BHP reaches 50 MPa
```

---

## Multi-Phase Wells

### Black Oil Well

```ini
[SIMULATION]
fluid_model = BLACK_OIL

[WELL1]
name = PROD-1
type = PRODUCER
i = 25
j = 25
k = 5

control_mode = RATE
rate_type = OIL                # Control oil rate
target_value = 0.001           # m³/s oil

# Constraints
max_water_cut = 0.95           # Shut if water cut > 95%
max_gor = 5000.0               # Shut if GOR too high
min_bhp = 10.0e6
```

### Gas Injection

```ini
[WELL2]
name = GAS-INJ-1
type = INJECTOR
fluid = GAS

i = 5
j = 5
k = 5

control_mode = RATE
rate_type = GAS
target_value = 1.0             # 1 m³/s gas at surface

# Composition for compositional model
injection_composition = 0.95, 0.03, 0.02  # C1, C2, C3
```

### Water Alternating Gas (WAG)

```ini
[WELL2]
name = WAG-1
type = INJECTOR

i = 5
j = 5
k = 5

# WAG schedule
wag_enabled = true
wag_water_slug = 30 days       # Water injection period
wag_gas_slug = 30 days         # Gas injection period
wag_ratio = 1.0                # Volume ratio

water_rate = 0.01              # m³/s during water slug
gas_rate = 1.0                 # m³/s during gas slug
```

---

## Advanced Well Types

### Multi-Completion Well

Wells with multiple perforated intervals:

```ini
[WELL1]
name = PROD-1
type = PRODUCER

# Multiple completions
num_completions = 3

# Completion 1 (shallow zone)
completion_1_i = 10
completion_1_j = 10
completion_1_k = 2
completion_1_diameter = 0.2
completion_1_skin = 0.0
completion_1_open = true

# Completion 2 (middle zone)
completion_2_i = 10
completion_2_j = 10
completion_2_k = 5
completion_2_diameter = 0.2
completion_2_skin = 1.5
completion_2_open = true

# Completion 3 (deep zone)
completion_3_i = 10
completion_3_j = 10
completion_3_k = 8
completion_3_diameter = 0.2
completion_3_skin = 0.0
completion_3_open = false      # Shut off initially
```

### Horizontal Well

```ini
[WELL1]
name = HORZ-1
type = PRODUCER
well_geometry = HORIZONTAL

# Heel location
heel_i = 5
heel_j = 25
heel_k = 5

# Toe location (horizontal extent)
toe_i = 45
toe_j = 25
toe_k = 5

# Or specify by length and azimuth
# horizontal_length = 1000.0   # meters
# azimuth = 90.0               # degrees from N

control_mode = RATE
target_value = 0.02

# All cells between heel and toe are perforated
```

### Deviated/Directional Well

```ini
[WELL1]
name = DEVIATED-1
type = PRODUCER
well_geometry = DIRECTIONAL

# Surface location
surface_x = 1000.0
surface_y = 1000.0

# Kickoff depth (where deviation starts)
kickoff_depth = 500.0          # meters TVD

# Build rate
build_rate = 3.0               # degrees per 30m

# Target inclination
target_inclination = 60.0      # degrees from vertical

# Target azimuth
target_azimuth = 45.0          # degrees from N

# Target depth
target_tvd = 2000.0            # True vertical depth

# Survey data (alternative to calculated trajectory)
# survey_file = well_survey.csv
```

### Multilateral Well

```ini
[WELL1]
name = MULTI-1
type = PRODUCER
well_geometry = MULTILATERAL

# Main bore
main_heel_i = 5
main_heel_j = 25
main_heel_k = 5
main_toe_i = 25
main_toe_j = 25
main_toe_k = 5

# Lateral 1
lateral_1_junction_md = 500.0  # Junction measured depth
lateral_1_azimuth = 135.0      # Branch direction
lateral_1_length = 500.0       # Branch length
lateral_1_open = true

# Lateral 2
lateral_2_junction_md = 800.0
lateral_2_azimuth = 225.0
lateral_2_length = 400.0
lateral_2_open = true

control_mode = RATE
target_value = 0.03
```

---

## Well Performance Analysis

### Well Index Calculation

FSRM calculates well indices using the Peaceman formula:

```
WI = (2π × kh × h) / (μ × ln(rₑ/rw) + S)
```

Override if needed:

```ini
[WELL1]
# Automatic calculation (default)
well_index_model = PEACEMAN

# Or specify directly
# well_index = 1.0e-11         # m³/(Pa·s)

# Or provide equivalent radius
# equivalent_radius = 100.0    # m
```

### Inflow Performance Relationship (IPR)

```ini
[WELL1]
# IPR model for productivity
ipr_model = VOGEL             # LINEAR, VOGEL, FETKOVITCH

# Vogel parameters (for solution gas drive)
# q/qmax = 1 - 0.2(Pwf/Pr) - 0.8(Pwf/Pr)²
productivity_index = 5.0e-11   # m³/(Pa·s)
```

### Wellbore Hydraulics

For pressure drop in the wellbore:

```ini
[WELL1]
enable_wellbore_model = true

# Geometry
wellbore_depth = 2000.0        # Total depth (m)
tubing_id = 0.1                # Inner diameter (m)
tubing_roughness = 1.5e-5      # Roughness (m)

# Fluid in wellbore
wellbore_fluid_density = 800.0 # kg/m³

# Friction model
friction_model = HAZEN_WILLIAMS  # or DARCY_WEISBACH, MOODY

# For gas wells
enable_gas_lift = false
```

### VFP Tables (Vertical Flow Performance)

```ini
[WELL1]
# Use pre-calculated VFP tables
use_vfp = true
vfp_table_file = vfp/well1.vfp

# Or specify VFP table number (Eclipse format)
# vfp_table_number = 1
```

---

## Well Groups and Scheduling

### Well Groups

Control multiple wells as a group:

```ini
[GROUP1]
name = FIELD
wells = PROD-1, PROD-2, PROD-3

control_mode = RATE
rate_type = OIL
target_value = 0.05            # Group oil rate

# Distribute among wells
allocation = EQUAL             # EQUAL, PRORATED, PRIORITY
```

### Well Scheduling

Time-dependent well operations:

```ini
[WELL1]
name = PROD-1

# Initial state
initial_state = OPEN

# Schedule changes
schedule_enabled = true

# At t = 0 days
event_1_time = 0.0
event_1_action = OPEN
event_1_rate = 0.01

# At t = 100 days, increase rate
event_2_time = 8640000.0       # 100 days in seconds
event_2_action = RATE
event_2_rate = 0.02

# At t = 365 days, convert to injector
event_3_time = 31536000.0
event_3_action = CONVERT
event_3_type = INJECTOR
event_3_fluid = WATER
event_3_rate = 0.015
```

---

## Practical Examples

### Five-Spot Pattern

Classic water injection pattern:

```ini
# Central producer
[WELL1]
name = PROD-1
type = PRODUCER
i = 25
j = 25
k = 5
control_mode = RATE
target_value = 0.02
min_bhp = 10.0e6

# Corner injectors
[WELL2]
name = INJ-1
type = INJECTOR
fluid = WATER
i = 5
j = 5
k = 5
control_mode = RATE
target_value = 0.005
max_bhp = 35.0e6

[WELL3]
name = INJ-2
type = INJECTOR
fluid = WATER
i = 45
j = 5
k = 5
control_mode = RATE
target_value = 0.005
max_bhp = 35.0e6

[WELL4]
name = INJ-3
type = INJECTOR
fluid = WATER
i = 5
j = 45
k = 5
control_mode = RATE
target_value = 0.005
max_bhp = 35.0e6

[WELL5]
name = INJ-4
type = INJECTOR
fluid = WATER
i = 45
j = 45
k = 5
control_mode = RATE
target_value = 0.005
max_bhp = 35.0e6
```

### Geothermal Doublet

Injection-production pair:

```ini
[WELL1]
name = INJECTION
type = INJECTOR
fluid = WATER

i = 10
j = 25
k = 10

control_mode = RATE
target_value = 0.03            # 30 L/s
temperature = 313.0            # 40°C reinjection
max_bhp = 35.0e6

[WELL2]
name = PRODUCTION
type = PRODUCER

i = 40
j = 25
k = 10

control_mode = RATE
target_value = 0.03            # Match injection
min_bhp = 5.0e6
```

### Shale Gas Multi-Stage

Horizontal well with staged fractures:

```ini
[SIMULATION]
enable_fractures = true

[WELL1]
name = SHALE-H1
type = PRODUCER
well_geometry = HORIZONTAL

heel_i = 5
heel_j = 25
heel_k = 5
toe_i = 95
toe_j = 25
toe_k = 5

control_mode = BHP
target_value = 5.0e6           # Low BHP for gas

# Associated fractures (see Fracture Tutorial)
associated_fractures = FRAC-1, FRAC-2, FRAC-3, FRAC-4, FRAC-5
```

---

## Well Output and Analysis

### Output Configuration

```ini
[OUTPUT]
# Well data output
well_data = true
well_output_frequency = 1      # Every timestep

# What to output
well_pressure = true           # BHP
well_rate = true               # Flow rates
well_cumulative = true         # Cumulative production
well_water_cut = true          # Water cut
well_gor = true                # Gas-oil ratio
```

### Post-Processing

Well data is written to CSV files:

```csv
# output/wells/PROD-1.csv
time,bhp,rate_oil,rate_water,rate_gas,cum_oil,cum_water,cum_gas,wcut,gor
0.0,25000000.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
3600.0,24500000.0,0.0095,0.0005,0.5,34.2,1.8,1800.0,0.05,52.6
...
```

Python analysis:

```python
import pandas as pd
import matplotlib.pyplot as plt

# Load well data
well = pd.read_csv('output/wells/PROD-1.csv')

# Plot production rate
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Oil rate
axes[0,0].plot(well['time']/86400, well['rate_oil']*86400)
axes[0,0].set_xlabel('Time (days)')
axes[0,0].set_ylabel('Oil Rate (m³/day)')

# BHP
axes[0,1].plot(well['time']/86400, well['bhp']/1e6)
axes[0,1].set_xlabel('Time (days)')
axes[0,1].set_ylabel('BHP (MPa)')

# Cumulative production
axes[1,0].plot(well['time']/86400, well['cum_oil'])
axes[1,0].set_xlabel('Time (days)')
axes[1,0].set_ylabel('Cumulative Oil (m³)')

# Water cut
axes[1,1].plot(well['time']/86400, well['wcut']*100)
axes[1,1].set_xlabel('Time (days)')
axes[1,1].set_ylabel('Water Cut (%)')

plt.tight_layout()
plt.savefig('well_performance.png')
```

---

## Troubleshooting

### Common Issues

| Issue | Cause | Solution |
|-------|-------|----------|
| Zero rate | BHP at limit | Increase target BHP or reduce constraints |
| Negative rate | Wrong well type | Check PRODUCER vs INJECTOR |
| Very high BHP | Skin too high | Reduce skin factor |
| Convergence issues | Rate too high | Reduce rate or use BHP control |
| No output | Well indices | Check well location is in domain |

### Debugging Tips

```bash
# Verbose well output
mpirun -np 4 ./fsrm -c config.config \
  -well_debug true \
  -well_trace PROD-1
```

Check well indices:

```ini
[OUTPUT]
well_index = true              # Output calculated WI
```

---

**Previous**: [← Physics Models](03_PHYSICS_MODELS.md) | **Next**: [Fracture Modeling →](05_FRACTURE_MODELING.md)
