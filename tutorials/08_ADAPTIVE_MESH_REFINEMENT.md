# Tutorial 08: Adaptive Mesh Refinement

Dynamically optimize mesh resolution for accuracy and efficiency.

## Table of Contents

1. [AMR Overview](#amr-overview)
2. [Refinement Criteria](#refinement-criteria)
3. [Refinement Strategies](#refinement-strategies)
4. [Configuration](#configuration)
5. [Feature-Based Refinement](#feature-based-refinement)
6. [Solution Transfer](#solution-transfer)
7. [Load Balancing](#load-balancing)
8. [Practical Examples](#practical-examples)

---

## AMR Overview

Adaptive Mesh Refinement (AMR) automatically adjusts grid resolution where needed:

```
┌─────────────────────────────────────────────────────────────────┐
│                    Adaptive Mesh Refinement                      │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  BEFORE AMR                    AFTER AMR                         │
│  ┌─────────────────┐          ┌─────────────────┐               │
│  │ │ │ │ │ │ │ │ │ │          │   │   │ │ │ │   │               │
│  │─┼─┼─┼─┼─┼─┼─┼─┼─│          │───┼───┼─┼─┼─┼───│               │
│  │ │ │ │ │ │ │ │ │ │          │   │ │ │ │ │ │   │               │
│  │─┼─┼─┼─┼─┼─┼─┼─┼─│          │───┼─┼─┼─┼─┼─┼───│               │
│  │ │ │ │ │ │ │ │ │ │    →     │ │ │ │▓│▓│ │ │ │ │               │
│  │─┼─┼─┼─┼─┼─┼─┼─┼─│          │─┼─┼─┼─┼─┼─┼─┼─┼─│               │
│  │ │ │ │ │ │ │ │ │ │          │   │ │ │ │ │ │   │               │
│  │─┼─┼─┼─┼─┼─┼─┼─┼─│          │───┼─┼─┼─┼─┼─┼───│               │
│  │ │ │ │ │ │ │ │ │ │          │   │   │ │ │ │   │               │
│  └─────────────────┘          └─────────────────┘               │
│   Uniform mesh                 Refined near feature ▓           │
│                                                                  │
│  Benefits:                                                       │
│  • Higher accuracy where needed                                  │
│  • Reduced computational cost                                    │
│  • Dynamic adaptation to solution                                │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
```

### AMR Methods in FSRM

| Method | Description | Best For |
|--------|-------------|----------|
| `PLEX_REFINE` | DMPlex uniform refinement | Structured meshes |
| `PLEX_ADAPT` | Metric-based adaptation | Unstructured meshes |
| `FOREST_P4EST` | Octree/quadtree (p4est) | Large-scale AMR |
| `PRAGMATIC` | Anisotropic adaptation | Complex geometries |
| `MMG` | Remeshing library | Quality-focused |

---

## Refinement Criteria

### Gradient-Based (Kelly Estimator)

Refine where solution gradients are large:

```ini
[AMR]
enabled = true
criterion_type = GRADIENT

# Field to estimate error
primary_field = PRESSURE       # or SATURATION, TEMPERATURE, etc.

# Gradient threshold
gradient_threshold = 1000.0    # Pa/m for pressure
```

**Error Estimate:**
```
η_K² = h_K × Σ ||[∇u · n]||²  (jump in gradient across faces)
```

### Hessian-Based

For smooth solutions, uses second derivatives:

```ini
[AMR]
criterion_type = HESSIAN

# Recovery method for Hessian
hessian_recovery = SPR         # SPR, L2_PROJECTION, PATCH_RECOVERY
```

### Residual-Based

Refine where equation residuals are large:

```ini
[AMR]
criterion_type = RESIDUAL

# Residual threshold
residual_threshold = 1.0e-6
```

### Jump-Based (Interface Tracking)

Track sharp interfaces (saturation fronts, phase boundaries):

```ini
[AMR]
criterion_type = JUMP

# Field with discontinuities
jump_field = SATURATION

# Jump threshold
jump_threshold = 0.1           # 10% saturation change
```

### Physics-Based

Refine based on physical indicators:

```ini
[AMR]
criterion_type = PHYSICS

# Physics indicators
refine_near_wells = true
refine_near_faults = true
refine_at_fronts = true
refine_high_velocity = true
velocity_threshold = 1.0e-5    # m/s
```

### Combined Criteria

Weight multiple indicators:

```ini
[AMR]
criterion_type = COMBINED

# Field weights (sum to 1)
weight_pressure = 0.3
weight_saturation = 0.4
weight_velocity = 0.2
weight_temperature = 0.1
```

---

## Refinement Strategies

### Fixed Fraction

Refine/coarsen a fixed fraction of cells:

```ini
[AMR]
strategy = FIXED_FRACTION

refine_fraction = 0.2          # Refine top 20% error cells
coarsen_fraction = 0.05        # Coarsen bottom 5%
```

### Threshold-Based

Refine/coarsen based on error thresholds:

```ini
[AMR]
strategy = THRESHOLD

refine_threshold = 0.5         # Refine if error > 50% of max
coarsen_threshold = 0.05       # Coarsen if error < 5% of max
```

### Fixed Number

Maintain approximately constant cell count:

```ini
[AMR]
strategy = FIXED_NUMBER

target_cells = 100000          # Aim for 100K cells
tolerance = 0.1                # ±10% variation allowed
```

### Error Equilibration

Distribute error evenly across cells:

```ini
[AMR]
strategy = EQUILIBRATION

target_error = 1.0e-4          # Total error budget
```

---

## Configuration

### Basic AMR Setup

```ini
[SIMULATION]
enable_amr = true

[AMR]
enabled = true
method = PLEX_ADAPT            # Adaptation method

# When to adapt
adapt_every = 10               # Every 10 timesteps
# Or adapt on significant change:
adapt_on_change = true
change_threshold = 0.1         # 10% change triggers adapt

# Criterion
criterion_type = GRADIENT
weight_pressure = 1.0

# Strategy
strategy = FIXED_FRACTION
refine_fraction = 0.2
coarsen_fraction = 0.05
```

### Mesh Limits

```ini
[AMR]
# Refinement levels
max_level = 5                  # Max refinement depth
min_level = 0                  # Min refinement (0 = coarse)

# Cell counts
max_cells = 1000000            # Hard limit
min_cells = 1000               # Minimum cells

# Cell sizes
min_cell_size = 1.0            # meters (finest resolution)
max_cell_size = 1000.0         # meters (coarsest resolution)
```

### Quality Constraints

```ini
[AMR]
# Cell quality
min_quality = 0.2              # Minimum scaled Jacobian
max_aspect_ratio = 10.0        # Maximum aspect ratio

# Quality measure
quality_measure = SCALED_JACOBIAN  # SHAPE, ASPECT_RATIO, SKEWNESS
```

### Boundary and Feature Preservation

```ini
[AMR]
# Preserve important features
preserve_boundaries = true     # Keep boundary resolution
preserve_wells = true          # Maintain well refinement
preserve_faults = true         # Keep fault mesh quality

# Buffer layers around features
buffer_layers = 2              # Transition cells
```

---

## Feature-Based Refinement

### Well Refinement

```ini
[AMR]
# Automatic well refinement
refine_near_wells = true
well_refinement_radius = 100.0 # meters
well_refinement_level = 3      # Refinement depth

# Or specify per well
[AMR_WELL1]
well_name = PROD-1
x = 500.0
y = 500.0
z = 50.0
radius = 150.0
level = 4
```

### Fault Refinement

```ini
[AMR]
refine_near_faults = true
fault_refinement_width = 50.0  # meters from fault
fault_refinement_level = 3

[AMR_FAULT1]
fault_name = MAIN_FAULT
trace = 0.0, 0.0, 2000.0, 2000.0  # x1, y1, x2, y2
width = 100.0
level = 4
```

### Box Refinement

```ini
[AMR_BOX1]
name = injection_zone
xmin = 400.0
xmax = 600.0
ymin = 400.0
ymax = 600.0
zmin = 40.0
zmax = 60.0
level = 3
```

### Front Tracking

```ini
[AMR]
# Track moving fronts
enable_front_tracking = true
front_field = SATURATION
front_value = 0.5              # Track S=0.5 contour
front_refinement_width = 20.0  # meters
front_refinement_level = 3
```

---

## Solution Transfer

When mesh changes, solution must be transferred:

### Transfer Methods

```ini
[AMR]
# Transfer method
transfer_method = INTERPOLATION  # INTERPOLATION, PROJECTION, CONSERVATIVE

# For conservative transfer (mass-preserving)
conservative = true
```

| Method | Description | Use Case |
|--------|-------------|----------|
| `INTERPOLATION` | FE interpolation | General purpose |
| `PROJECTION` | L2 projection | Smooth fields |
| `INJECTION` | Direct copy | Coarsening |
| `CONSERVATIVE` | Mass-preserving | Flow quantities |

### Multi-Field Transfer

```ini
[AMR]
# Fields to transfer
transfer_pressure = true
transfer_saturation = true
transfer_displacement = true
transfer_temperature = true

# Per-field methods
pressure_transfer = CONSERVATIVE
saturation_transfer = CONSERVATIVE
displacement_transfer = INTERPOLATION
```

---

## Load Balancing

Redistribute mesh across processors after refinement:

```ini
[AMR]
# Enable load balancing
balance_enabled = true

# Balancing method
balance_method = PARMETIS      # PARMETIS, PTSCOTCH, SIMPLE

# When to rebalance
rebalance_threshold = 0.2      # Rebalance if imbalance > 20%

# Balancing options
partition_overlap = 1          # Overlap cells for communication
weight_field = PRESSURE        # Use field for weighting (optional)
```

### Partitioning Options

| Method | Description | Best For |
|--------|-------------|----------|
| `PARMETIS` | Graph partitioning | General meshes |
| `PTSCOTCH` | Alternative graph | Large meshes |
| `SIMPLE` | Geometric | Structured meshes |
| `SFC` | Space-filling curve | Octree meshes |

---

## Practical Examples

### Water Flooding with Front Tracking

```ini
[SIMULATION]
name = waterflood_amr
fluid_model = BLACK_OIL
enable_amr = true
end_time = 31536000.0          # 1 year

[GRID]
nx = 50
ny = 50
nz = 5
Lx = 1000.0
Ly = 1000.0
Lz = 50.0

[AMR]
enabled = true
method = PLEX_ADAPT

# Criterion: track saturation front
criterion_type = COMBINED
weight_pressure = 0.2
weight_saturation = 0.8        # Focus on saturation

# Strategy
strategy = FIXED_FRACTION
refine_fraction = 0.15
coarsen_fraction = 0.05

# Limits
max_level = 4
min_cell_size = 5.0            # 5m minimum
max_cells = 200000

# Front tracking
enable_front_tracking = true
front_field = SATURATION
front_value = 0.5              # Water front at Sw = 0.5

# When to adapt
adapt_every = 20               # Every 20 steps

# Transfer
transfer_method = CONSERVATIVE
conservative = true

[WELL1]
name = INJ-1
type = INJECTOR
fluid = WATER
i = 5
j = 5
k = 2
control_mode = RATE
target_value = 0.01

[WELL2]
name = PROD-1
type = PRODUCER
i = 45
j = 45
k = 2
control_mode = RATE
target_value = 0.01
```

### Induced Seismicity with Fault Refinement

```ini
[SIMULATION]
name = seismicity_amr
solid_model = POROELASTIC
enable_geomechanics = true
enable_faults = true
enable_amr = true

[AMR]
enabled = true
method = PLEX_ADAPT

# Combined criteria
criterion_type = COMBINED
weight_pressure = 0.3
weight_velocity = 0.3
weight_stress = 0.4

# Strategy
strategy = THRESHOLD
refine_threshold = 0.3
coarsen_threshold = 0.05

# Limits
max_level = 4
max_cells = 500000

# Feature preservation
preserve_faults = true
preserve_wells = true
buffer_layers = 2

# Fault refinement
[AMR_FAULT1]
fault_name = MAIN_FAULT
trace = 5000.0, 0.0, 5000.0, 10000.0
width = 200.0
level = 4

# Well refinement
[AMR_WELL1]
well_name = INJ-1
x = 4500.0
y = 5000.0
z = 2500.0
radius = 500.0
level = 3

# Load balancing
balance_enabled = true
balance_method = PARMETIS
rebalance_threshold = 0.15
```

### Geothermal with Thermal Fronts

```ini
[SIMULATION]
name = geothermal_amr
enable_thermal = true
enable_amr = true

[AMR]
enabled = true
method = PLEX_ADAPT

# Criteria
criterion_type = COMBINED
weight_pressure = 0.3
weight_temperature = 0.5
weight_velocity = 0.2

# Track thermal front
enable_front_tracking = true
front_field = TEMPERATURE
front_value = 350.0            # Track T = 350K contour

# Strategy
strategy = EQUILIBRATION
target_error = 1.0e-3

max_level = 4
min_cell_size = 10.0

# Feature refinement
[AMR_WELL1]
well_name = INJECTION
x = 1000.0
y = 5000.0
z = 2000.0
radius = 300.0
level = 4

[AMR_WELL2]
well_name = PRODUCTION
x = 9000.0
y = 5000.0
z = 2000.0
radius = 300.0
level = 4
```

---

## AMR Output and Visualization

### Output Configuration

```ini
[OUTPUT]
# AMR-specific outputs
write_mesh = true              # Output adapted mesh
write_errors = true            # Output error field
write_levels = true            # Output refinement levels

mesh_output_format = VTK       # VTK, EXODUS, HDF5
```

### Visualization in ParaView

1. Load mesh files: `output/mesh_*.vtu`
2. Color by `RefinementLevel` to see adaptation
3. Color by `ErrorIndicator` to see refinement drivers
4. Use `Threshold` filter to show only refined regions

### Post-Processing Script

```python
import numpy as np
import matplotlib.pyplot as plt

# Load AMR statistics
stats = np.loadtxt('output/amr_statistics.csv', delimiter=',', skiprows=1)

time = stats[:, 0]
num_cells = stats[:, 1]
max_level = stats[:, 2]
max_error = stats[:, 3]
mean_error = stats[:, 4]

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Cell count
axes[0, 0].plot(time / 86400, num_cells / 1000)
axes[0, 0].set_xlabel('Time (days)')
axes[0, 0].set_ylabel('Cells (thousands)')
axes[0, 0].set_title('Mesh Size Evolution')

# Max refinement level
axes[0, 1].plot(time / 86400, max_level)
axes[0, 1].set_xlabel('Time (days)')
axes[0, 1].set_ylabel('Max Level')
axes[0, 1].set_title('Refinement Depth')

# Error evolution
axes[1, 0].semilogy(time / 86400, max_error, label='Max')
axes[1, 0].semilogy(time / 86400, mean_error, label='Mean')
axes[1, 0].set_xlabel('Time (days)')
axes[1, 0].set_ylabel('Error')
axes[1, 0].legend()
axes[1, 0].set_title('Error Indicators')

# Efficiency (error per cell)
efficiency = mean_error * num_cells
axes[1, 1].plot(time / 86400, efficiency / efficiency[0])
axes[1, 1].set_xlabel('Time (days)')
axes[1, 1].set_ylabel('Relative Total Error')
axes[1, 1].set_title('Efficiency')

plt.tight_layout()
plt.savefig('amr_analysis.png', dpi=150)
```

---

## Best Practices

1. **Start coarse**: Begin with a coarse mesh and let AMR refine
2. **Set conservative limits**: Don't allow too many cells
3. **Use appropriate criteria**: Match refinement to physics
4. **Test transfer**: Verify mass/energy conservation
5. **Monitor balance**: Check load distribution
6. **Profile performance**: AMR has overhead; ensure benefit outweighs cost

### When AMR Helps Most

✅ **Good for:**
- Moving fronts (flooding, injection)
- Localized features (wells, faults)
- Multi-scale problems
- Long-time simulations

❌ **Less beneficial for:**
- Uniform problems
- Very small domains
- Short simulations
- Highly dynamic problems (constant remeshing)

---

**Previous**: [← GPU Acceleration](07_GPU_ACCELERATION.md) | **Next**: [Wave Propagation →](09_WAVE_PROPAGATION.md)
