# Tutorial 10: Advanced Simulations

Combine multiple physics for complex, real-world simulation scenarios.

## Table of Contents

1. [Multi-Physics Coupling](#multi-physics-coupling)
2. [CO2 Storage Simulation](#co2-storage-simulation)
3. [Geothermal Systems](#geothermal-systems)
4. [Shale Gas Production](#shale-gas-production)
5. [Enhanced Oil Recovery](#enhanced-oil-recovery)
6. [Nuclear Waste Disposal](#nuclear-waste-disposal)
7. [Mining and Excavation](#mining-and-excavation)
8. [Workflow Automation](#workflow-automation)

---

## Multi-Physics Coupling

### Coupling Framework

FSRM supports various coupling strategies:

```
┌──────────────────────────────────────────────────────────────┐
│                  Multi-Physics Coupling                       │
├──────────────────────────────────────────────────────────────┤
│                                                               │
│  COUPLING TYPES                                               │
│  ├─ Sequential: Physics solved one after another             │
│  ├─ Iterative: Multiple passes until convergence             │
│  └─ Monolithic: All physics in single system                 │
│                                                               │
│  COMMON COMBINATIONS                                          │
│  ├─ HM: Flow ↔ Mechanics (reservoir compaction)             │
│  ├─ TH: Flow ↔ Thermal (geothermal)                         │
│  ├─ TM: Thermal ↔ Mechanics (thermal stress)                │
│  ├─ THM: Full coupling (geothermal, EGS)                    │
│  └─ THMC: + Chemistry (nuclear waste, CO2)                  │
│                                                               │
│  COUPLING MECHANISMS                                          │
│  ├─ Pore pressure → Effective stress                        │
│  ├─ Deformation → Porosity/permeability change              │
│  ├─ Temperature → Fluid properties                          │
│  ├─ Temperature → Thermal stress                            │
│  └─ Chemistry → Mineral precipitation/dissolution           │
│                                                               │
└──────────────────────────────────────────────────────────────┘
```

### Full THM Configuration

```ini
[SIMULATION]
name = thm_coupled
fluid_model = SINGLE_COMPONENT
solid_model = POROELASTIC
enable_geomechanics = true
enable_thermal = true

# Coupling settings
coupling_type = ITERATIVE      # SEQUENTIAL, ITERATIVE, MONOLITHIC
coupling_tolerance = 1.0e-6
max_coupling_iterations = 10

# Physics flags
flow_mechanics_coupling = true
flow_thermal_coupling = true
mechanics_thermal_coupling = true
```

### Coupling Iteration

```ini
[COUPLING]
# Iterative coupling settings
staggered_scheme = FIXED_STRESS  # FIXED_STRESS, FIXED_STRAIN, UNDRAINED

# Under-relaxation for stability
relaxation_pressure = 0.8
relaxation_displacement = 1.0
relaxation_temperature = 0.9

# Convergence monitoring
convergence_field = PRESSURE     # Primary field to check
```

---

## CO2 Storage Simulation

### Saline Aquifer Storage

```ini
[SIMULATION]
name = co2_storage
fluid_model = CO2
solid_model = POROELASTIC
enable_geomechanics = true
enable_faults = true
end_time = 315360000.0         # 10 years

[GRID]
nx = 100
ny = 100
nz = 20
Lx = 10000.0                   # 10 km
Ly = 10000.0
Lz = 200.0                     # 200m reservoir

[ROCK]
constitutive_model = POROELASTIC
porosity = 0.20
permeability_x = 200.0         # mD
permeability_z = 20.0
youngs_modulus = 15.0e9
poisson_ratio = 0.25
biot_coefficient = 0.85

# Caprock properties
[ROCK1]
name = caprock
region = 1
porosity = 0.02
permeability_x = 0.001         # Very low perm
youngs_modulus = 25.0e9

[FLUID]
type = CO2
co2_model = SPAN_WAGNER

# Brine properties
brine_density = 1100.0         # kg/m³
salinity = 100000              # ppm

# Mutual solubility
enable_dissolution = true
henry_coefficient = 1.5e-5

[WELL1]
name = CO2_INJ_1
type = INJECTOR
fluid = CO2
i = 50
j = 50
k = 10
control_mode = RATE
target_value = 0.5             # 0.5 kg/s CO2
max_bhp = 25.0e6               # Fracture pressure limit

# Nearby fault
[FAULT1]
name = MONITORING_FAULT
x = 7000.0
y = 5000.0
z = 100.0
strike = 30.0
dip = 70.0
length = 5000.0
width = 200.0
friction_law = COULOMB
static_friction = 0.6

# Initial hydrostatic pressure
[IC1]
field = PRESSURE
distribution = GRADIENT
value = 10.0e6                 # Surface pressure
gradient = 0.0, 0.0, 10000.0   # Hydrostatic

# Geothermal gradient
[IC2]
field = TEMPERATURE
distribution = GRADIENT
value = 288.0                  # Surface T (15°C)
gradient = 0.0, 0.0, 0.03      # 30°C/km

[SEISMICITY]
enable = true
max_magnitude = 3.0            # Safety threshold

[OUTPUT]
pressure = true
saturation = true              # CO2 saturation
displacement = true
fault_slip = true
```

### Monitoring Points

```ini
[MONITORING]
# Pressure monitoring wells
monitor_1 = 6000.0, 5000.0, 100.0
monitor_2 = 8000.0, 5000.0, 100.0
monitor_3 = 5000.0, 7000.0, 100.0

# Surface deformation (for InSAR comparison)
surface_displacement = true
surface_nx = 50
surface_ny = 50
```

---

## Geothermal Systems

### Enhanced Geothermal System (EGS)

```ini
[SIMULATION]
name = egs_simulation
fluid_model = SINGLE_COMPONENT
solid_model = POROELASTIC
enable_geomechanics = true
enable_thermal = true
enable_fractures = true
end_time = 946080000.0         # 30 years

[GRID]
nx = 80
ny = 80
nz = 40
Lx = 2000.0
Ly = 2000.0
Lz = 1000.0

# Start at reservoir depth
origin_z = 3000.0              # 3 km depth

[ROCK]
constitutive_model = POROELASTIC
porosity = 0.01                # Crystalline rock
permeability_x = 0.1           # Very low matrix perm
permeability_y = 0.1
permeability_z = 0.01
youngs_modulus = 50.0e9        # Granite
poisson_ratio = 0.25
biot_coefficient = 0.6

# Thermal properties
thermal_conductivity = 3.0     # W/(m·K)
specific_heat = 900.0          # J/(kg·K)
thermal_expansion = 8.0e-6     # 1/K

[FLUID]
type = SINGLE_PHASE
density = 900.0                # Hot water at depth
viscosity = 0.0002             # Low viscosity at high T
compressibility = 5.0e-10
specific_heat = 4200.0
thermal_conductivity = 0.65

# Injection well (cold water)
[WELL1]
name = INJECTION
type = INJECTOR
fluid = WATER
i = 20
j = 40
k = 20
control_mode = RATE
target_value = 0.05            # 50 kg/s
temperature = 323.0            # 50°C injection
max_bhp = 50.0e6

# Production well (hot water)
[WELL2]
name = PRODUCTION
type = PRODUCER
i = 60
j = 40
k = 20
control_mode = RATE
target_value = 0.05
min_bhp = 20.0e6

# Stimulated fracture zone connecting wells
[FRACTURE1]
type = DFN_STOCHASTIC
num_fractures = 200
region = 10, 70, 30, 50, 10, 30  # Between wells

size_distribution = POWER_LAW
min_length = 20.0
max_length = 300.0
exponent = 2.5

aperture_model = LENGTH_CORRELATED
aperture_coefficient = 5.0e-6

# Initial temperature at depth
[IC1]
field = TEMPERATURE
distribution = UNIFORM
value = 473.0                  # 200°C at 3km depth

[IC2]
field = PRESSURE
distribution = UNIFORM
value = 30.0e6                 # Hydrostatic at 3km

[OUTPUT]
frequency = 365                # Yearly output
pressure = true
temperature = true
displacement = true

# Thermal output
heat_extraction_rate = true
thermal_breakthrough = true
```

### Hydrothermal Convection

```ini
[SIMULATION]
name = hydrothermal
enable_thermal = true
fluid_model = SINGLE_COMPONENT

[THERMAL]
enable_convection = true
enable_conduction = true
rayleigh_number_based = true   # Auto-check stability

# Buoyancy
enable_buoyancy = true
thermal_expansion_fluid = 4.0e-4  # 1/K

[BC1]
type = DIRICHLET
field = TEMPERATURE
location = ZMIN
value = 523.0                  # Hot base (250°C)

[BC2]
type = DIRICHLET
field = TEMPERATURE
location = ZMAX
value = 293.0                  # Cool top (20°C)
```

---

## Shale Gas Production

### Multi-Stage Fractured Horizontal Well

```ini
[SIMULATION]
name = shale_gas
fluid_model = SINGLE_COMPONENT
solid_model = POROELASTIC
enable_geomechanics = true
enable_fractures = true
end_time = 315360000.0         # 10 years

[GRID]
nx = 200
ny = 100
nz = 20
Lx = 3000.0                    # 3 km along well
Ly = 1500.0                    # 1.5 km lateral
Lz = 100.0                     # 100m pay zone

[ROCK]
constitutive_model = POROELASTIC
porosity = 0.06
permeability_x = 0.0001        # 100 nD
permeability_y = 0.0001
permeability_z = 0.00001
youngs_modulus = 35.0e9        # Stiff shale
poisson_ratio = 0.25
biot_coefficient = 0.9

# Adsorption for gas
langmuir_volume = 0.003        # m³/kg
langmuir_pressure = 5.0e6      # Pa

[FLUID]
type = SINGLE_PHASE
density = 100.0                # Gas at reservoir conditions
viscosity = 2.0e-5             # Gas viscosity
compressibility = 5.0e-8       # High compressibility

# Horizontal well
[WELL1]
name = HORZ_PROD
type = PRODUCER
well_geometry = HORIZONTAL
heel_i = 20
heel_j = 50
heel_k = 10
toe_i = 180
toe_j = 50
toe_k = 10
control_mode = BHP
target_value = 5.0e6           # Low BHP for gas

# Multi-stage fractures (10 stages)
[FRACTURE1]
type = HYDRAULIC
center = 300.0, 750.0, 50.0
half_length = 200.0
height = 80.0
conductivity = 1.0e-13

[FRACTURE2]
type = HYDRAULIC
center = 600.0, 750.0, 50.0
half_length = 200.0
height = 80.0
conductivity = 1.0e-13

# ... (repeat for all stages at 300m spacing)

[FRACTURE10]
type = HYDRAULIC
center = 2700.0, 750.0, 50.0
half_length = 200.0
height = 80.0
conductivity = 1.0e-13

[IC1]
field = PRESSURE
distribution = UNIFORM
value = 35.0e6                 # Overpressured

[OUTPUT]
pressure = true
rate_gas = true
cumulative_gas = true
```

---

## Enhanced Oil Recovery

### Water-Alternating-Gas (WAG)

```ini
[SIMULATION]
name = wag_eor
fluid_model = BLACK_OIL
end_time = 157680000.0         # 5 years

[ROCK]
porosity = 0.22
permeability_x = 150.0
permeability_z = 15.0

# Relative permeability (three-phase)
[RELPERM]
model = STONE_II               # Three-phase model
Swc = 0.2
Sor_wg = 0.15                  # Residual to WAG
krw_max = 0.35
kro_max = 0.9
krg_max = 0.8

[FLUID]
type = BLACK_OIL
oil_density_std = 850.0
gas_density_std = 0.9
water_density_std = 1000.0
solution_gor = 150.0

# Producer
[WELL1]
name = PROD-1
type = PRODUCER
i = 45
j = 45
k = 5
control_mode = RATE
rate_type = LIQUID
target_value = 0.01

# WAG Injector
[WELL2]
name = WAG-1
type = INJECTOR
i = 5
j = 5
k = 5

wag_enabled = true
wag_water_slug = 2592000.0     # 30 days water
wag_gas_slug = 2592000.0       # 30 days gas
wag_ratio = 1.0

water_rate = 0.01              # m³/s
gas_rate = 1.0                 # m³/s
max_bhp = 30.0e6

[IC1]
field = SATURATION
phase = OIL
distribution = UNIFORM
value = 0.75                   # 75% oil initially

[IC2]
field = SATURATION
phase = WATER
distribution = UNIFORM
value = 0.25                   # 25% connate water
```

### Polymer Flooding

```ini
[SIMULATION]
name = polymer_flood
fluid_model = BLACK_OIL
enable_tracer = true

[TRACER]
name = POLYMER
type = POLYMER

# Polymer properties
viscosity_multiplier = 10.0    # Viscosity increase
adsorption = 0.0001            # kg/kg rock
degradation_rate = 1.0e-8      # 1/s
inaccessible_pv = 0.1          # Fraction inaccessible

[WELL2]
name = POLY_INJ
type = INJECTOR
fluid = WATER
i = 5
j = 5
k = 5
control_mode = RATE
target_value = 0.008

# Polymer schedule
polymer_injection = true
polymer_concentration = 0.001  # kg/m³
polymer_start_time = 31536000  # Start after 1 year
polymer_end_time = 94608000    # End after 3 years
```

---

## Nuclear Waste Disposal

### Deep Geological Repository

```ini
[SIMULATION]
name = nuclear_repository
fluid_model = SINGLE_COMPONENT
solid_model = POROELASTIC
enable_geomechanics = true
enable_thermal = true
end_time = 3153600000000.0     # 100,000 years!

# Special time stepping for long durations
[TIMESTEPPING]
dt_initial = 86400.0           # 1 day
dt_max = 31536000000.0         # 1000 years max step
adaptive_dt = true
growth_factor = 1.5

[GRID]
nx = 50
ny = 50
nz = 50
Lx = 1000.0
Ly = 1000.0
Lz = 500.0

# Host rock (clay/granite)
[ROCK]
constitutive_model = VISCOELASTIC_MAXWELL
porosity = 0.15
permeability_x = 1.0e-6        # Ultra-low perm
youngs_modulus = 10.0e9
poisson_ratio = 0.3
viscosity = 1.0e20             # Creep behavior

thermal_conductivity = 2.0
specific_heat = 800.0
thermal_expansion = 1.0e-5

# Heat source (waste canister)
[HEAT_SOURCE]
type = DECAYING
location = 500.0, 500.0, 250.0
initial_power = 1000.0         # W
decay_constant = 7.3e-10       # 1/s (~30 year half-life)
radius = 2.0                   # Source region

[BC1]
type = DIRICHLET
field = TEMPERATURE
location = ALL
value = 288.0                  # 15°C

[OUTPUT]
# Long-term monitoring
frequency = 100
temperature = true
displacement = true
stress = true

# Safety metrics
max_temperature = true
effective_stress = true
```

---

## Mining and Excavation

### Underground Excavation

```ini
[SIMULATION]
name = tunnel_excavation
solid_model = ELASTOPLASTIC
enable_geomechanics = true
end_time = 31536000.0          # 1 year

[ROCK]
constitutive_model = ELASTOPLASTIC_MC
youngs_modulus = 25.0e9
poisson_ratio = 0.25
density = 2700.0

failure_criterion = MOHR_COULOMB
cohesion = 5.0e6
friction_angle = 35.0
dilation_angle = 10.0

# In-situ stress
[IC_STRESS]
stress_type = ANDERSON
tectonic_regime = STRIKE_SLIP
vertical_gradient = 27000.0    # ρg
stress_ratio = 1.5             # σH/σv

# Excavation sequence
[EXCAVATION]
enabled = true
method = ELEMENT_REMOVAL

# Stage 1: Top heading
stage_1_time = 0.0
stage_1_region = 490, 510, 490, 510, 245, 255

# Stage 2: Bench
stage_2_time = 2592000.0       # 30 days later
stage_2_region = 490, 510, 490, 510, 235, 245

# Support installation
[SUPPORT]
type = SHOTCRETE
thickness = 0.15               # 15 cm
youngs_modulus = 20.0e9
application_time = 86400.0     # 1 day after excavation

[OUTPUT]
displacement = true
stress = true
plastic_strain = true
safety_factor = true           # F = strength/stress
```

---

## Workflow Automation

### Python Workflow Script

```python
#!/usr/bin/env python3
"""
Advanced FSRM Workflow Automation
"""

import os
import subprocess
import numpy as np
import pandas as pd
from jinja2 import Template

# Configuration template
CONFIG_TEMPLATE = """
[SIMULATION]
name = {{ name }}
fluid_model = SINGLE_COMPONENT
solid_model = POROELASTIC
enable_geomechanics = true
end_time = {{ end_time }}

[GRID]
nx = {{ nx }}
ny = {{ ny }}
nz = {{ nz }}
Lx = {{ Lx }}
Ly = {{ Ly }}
Lz = {{ Lz }}

[ROCK]
porosity = {{ porosity }}
permeability_x = {{ perm }}
youngs_modulus = {{ E }}
poisson_ratio = 0.25
biot_coefficient = {{ biot }}

[WELL1]
name = INJ-1
type = INJECTOR
fluid = WATER
i = {{ well_i }}
j = {{ well_j }}
k = {{ well_k }}
control_mode = RATE
target_value = {{ rate }}

[OUTPUT]
path = output/{{ name }}
pressure = true
displacement = true
"""

def run_parameter_study():
    """Run parameter sensitivity study"""
    
    # Parameter ranges
    permeabilities = [10, 50, 100, 200, 500]
    biots = [0.5, 0.7, 0.9]
    
    results = []
    
    for perm in permeabilities:
        for biot in biots:
            name = f"perm_{perm}_biot_{int(biot*10)}"
            
            # Generate config
            template = Template(CONFIG_TEMPLATE)
            config = template.render(
                name=name,
                end_time=86400.0,
                nx=50, ny=50, nz=10,
                Lx=1000.0, Ly=1000.0, Lz=100.0,
                porosity=0.2,
                perm=perm,
                E=20.0e9,
                biot=biot,
                well_i=25, well_j=25, well_k=5,
                rate=0.01
            )
            
            # Write config
            config_file = f"configs/{name}.config"
            os.makedirs("configs", exist_ok=True)
            with open(config_file, 'w') as f:
                f.write(config)
            
            # Run simulation
            print(f"Running: {name}")
            subprocess.run(
                ["mpirun", "-np", "4", "./fsrm", "-c", config_file],
                check=True
            )
            
            # Collect results
            summary = pd.read_csv(f"output/{name}/summary.csv")
            results.append({
                'name': name,
                'perm': perm,
                'biot': biot,
                'max_pressure': summary['pressure_max'].max(),
                'final_displacement': summary['displacement_max'].iloc[-1]
            })
    
    # Save results
    df = pd.DataFrame(results)
    df.to_csv('parameter_study_results.csv', index=False)
    
    return df

def run_optimization():
    """Optimize injection rate for target pressure"""
    from scipy.optimize import minimize_scalar
    
    def objective(rate):
        """Run simulation and return pressure deviation from target"""
        name = f"opt_rate_{rate:.4f}"
        
        # Generate and run simulation
        # ... (similar to above)
        
        # Return objective (deviation from 25 MPa target)
        summary = pd.read_csv(f"output/{name}/summary.csv")
        max_p = summary['pressure_max'].max()
        return (max_p - 25e6)**2
    
    result = minimize_scalar(objective, bounds=(0.001, 0.1), method='bounded')
    print(f"Optimal rate: {result.x:.4f} m³/s")
    return result.x

if __name__ == "__main__":
    # Run parameter study
    results = run_parameter_study()
    print("\nParameter Study Results:")
    print(results)
```

### Batch Processing

```bash
#!/bin/bash
# batch_run.sh - Run multiple simulations

CONFIGS=(
    "config/co2_storage.config"
    "config/geothermal.config"
    "config/shale_reservoir.config"
)

NP=8  # Processors

for config in "${CONFIGS[@]}"; do
    name=$(basename "$config" .config)
    echo "Running: $name"
    
    mpirun -np $NP ./fsrm -c "$config" \
        -o "output/$name" \
        2>&1 | tee "logs/$name.log"
    
    # Post-process
    python scripts/analyze_results.py "output/$name"
done

echo "All simulations complete"
```

---

## Best Practices for Complex Simulations

### 1. Start Simple
```ini
# First: single physics, coarse mesh
# Then: add physics, refine mesh
```

### 2. Use Checkpointing
```ini
[OUTPUT]
checkpoint_frequency = 100
checkpoint_path = checkpoints/
```

### 3. Monitor Convergence
```bash
mpirun -np 8 ./fsrm -c config.config \
    -snes_monitor \
    -ksp_monitor_true_residual \
    -log_view
```

### 4. Validate Against Benchmarks
- SPE benchmarks (reservoirs)
- SCEC benchmarks (earthquakes)
- DECOVALEX benchmarks (THM)

### 5. Document Everything
```ini
# Header in config file:
# ====================
# Project: CO2 Storage Site A
# Author: Your Name
# Date: 2024-01-15
# Description: 10-year injection simulation
# ====================
```

---

**Previous**: [← Wave Propagation](09_WAVE_PROPAGATION.md) | **Back to Index**: [Tutorials Index →](README.md)
