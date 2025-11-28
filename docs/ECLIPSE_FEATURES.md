# ECLIPSE Reservoir Simulator Features

This document describes the ECLIPSE-compatible features implemented in FSRM (Fuck Stanford Reservoir Model). These features enable industry-standard reservoir simulation workflows.

## Table of Contents

1. [Aquifer Models](#aquifer-models)
2. [Tracer Simulation](#tracer-simulation)
3. [VFP Tables](#vfp-tables)
4. [Group Control](#group-control)
5. [Relative Permeability](#relative-permeability)
6. [Summary Output](#summary-output)
7. [Configuration Examples](#configuration-examples)

---

## Aquifer Models

FSRM supports multiple aquifer modeling approaches for simulating water influx from connected aquifers.

### Supported Aquifer Types

| Type | ECLIPSE Keyword | Description |
|------|-----------------|-------------|
| Carter-Tracy | `AQUCT` | Analytical infinite-acting aquifer using superposition |
| Fetkovich | `AQUFETP` | Pseudo-steady state model for depleting aquifers |
| Numerical | `AQUNUM` | Gridded aquifer cells for complex geometries |
| Constant Pressure | - | Simple pressure boundary condition |
| Constant Flux | - | Fixed influx rate boundary |

### Carter-Tracy Aquifer

The Carter-Tracy model is the most accurate for unsteady-state aquifer response:

```ini
[AQUIFER1]
type = CARTER_TRACY
id = 1
geometry = RADIAL                     # RADIAL, LINEAR, BOTTOM, EDGE
initial_pressure = 25.5e6             # Pa
permeability = 150.0                  # mD
porosity = 0.22
thickness = 30.0                      # m
compressibility = 4.5e-10             # 1/Pa
viscosity = 0.0005                    # Pa·s
inner_radius = 1000.0                 # Reservoir radius (m)
outer_radius = 50000.0                # Outer radius (m)
encroachment_angle = 0.75             # Fraction of circle
```

### Fetkovich Aquifer

Simpler model for pseudo-steady state conditions:

```ini
[AQUIFER2]
type = FETKOVICH
id = 2
initial_pressure = 26e6               # Pa
productivity_index = 5e-7             # m³/s/Pa
initial_volume = 2e9                  # Encroachable water volume (m³)
```

### Aquifer Connections

Connect aquifers to reservoir cells:

```ini
[AQUIFER_CONNECTION1]
aquifer_id = 1
i1 = 1
i2 = 1
j1 = 1
j2 = 20
k1 = 1
k2 = 5
face = I_MINUS                        # I_MINUS, I_PLUS, J_MINUS, J_PLUS, K_MINUS, K_PLUS
trans_mult = 1.0
allow_crossflow = true
```

### C++ API

```cpp
#include "AquiferModel.hpp"

// Create Carter-Tracy aquifer
CarterTracyAquifer aquifer(1);
aquifer.setAquiferProperties(100.0, 0.2, 1e-9, 0.001, 50.0, 0.5);
aquifer.setAquiferRadius(1000.0, 10000.0);
aquifer.setInitialPressure(30e6);

// Calculate influx
double influx = aquifer.calculateInflux(reservoir_pressure, dt);
aquifer.updateState(influx, dt);

// Use manager for multiple aquifers
AquiferManager manager;
manager.addCarterTracyAquifer(1, config);
manager.addFetkovichAquifer(2, config);
double total_influx = manager.calculateTotalInflux(reservoir_pressure, dt);
```

---

## Tracer Simulation

FSRM supports passive and reactive tracer simulation for inter-well testing and flow characterization.

### Tracer Types

| Type | Description |
|------|-------------|
| Passive | Inert tracer that follows fluid flow |
| Adsorbing | Tracer that adsorbs to rock (retardation) |
| Decaying | Radioactive tracer with first-order decay |
| Partitioning | Tracer that partitions between phases (SWCT) |
| Reactive | Tracer with chemical reactions |

### Configuration Example

```ini
[TRACER1]
name = TRACER_W1
phase = WATER                         # WATER, OIL, GAS
behavior = PASSIVE                    # PASSIVE, ADSORBING, DECAYING, PARTITIONING
molecular_weight = 100.0
diffusion_coefficient = 1e-9
dispersivity_L = 1.0
dispersivity_T = 0.1

[TRACER2]
name = TRACER_PART
phase = WATER
behavior = PARTITIONING
partition_coeff_ow = 2.0              # Oil-water partition coefficient

[TRACER_INJECTION1]
tracer = TRACER_W1
well = INJ1
start_time = 0.0
end_time = 86400.0                    # 1 day slug
concentration = 1.0                   # kg/m³
```

### C++ API

```cpp
#include "TracerModel.hpp"

TracerManager manager;
manager.setGrid(nx, ny, nz, dx, dy, dz);
manager.setPorosity(porosity);

// Add tracers
manager.addWaterTracer("TRACER1");
manager.addPartitioningTracer("PART_TRACER", 2.0);  // Kow = 2

// Add injection
manager.addInjection("TRACER1", "INJ1", 1.0, 0.0, 86400.0);

// Simulate
manager.updateFlowField(vx, vy, vz, Sw, So, Sg);
manager.step(dt);

// Analyze breakthrough
auto curve = manager.getBreakthrough("TRACER1", "PROD1");
double Sor = manager.estimateResidualOil("PART_TRACER", "TRACER1", "PROD1");
```

### SWCT Analysis

Single-Well Chemical Tracer tests can estimate residual oil saturation:

```cpp
double Sor = TracerAnalysis::estimateSorFromSWCT(
    partitioning_curve, passive_curve, Kow);
```

---

## VFP Tables

Vertical Flow Performance tables relate wellhead conditions to bottomhole pressure.

### VFPPROD - Production Wells

5-dimensional tables: BHP = f(THP, rate, WCT, GOR, ALQ)

```ini
[VFP_PROD1]
table_number = 1
datum_depth = 2000.0                  # m
flow_type = OIL                       # OIL, LIQ, GAS, WATER
wct_type = WCT                        # WCT, WGR
glr_type = GOR                        # GOR, GLR, OGR
alq_type = GASLIFT                    # NONE, GASLIFT, ESP, CHOKE

# Axis values
thp_values = 1e6, 2e6, 3e6, 4e6
rate_values = 0.001, 0.005, 0.01, 0.05, 0.1
wct_values = 0.0, 0.3, 0.6, 0.9
gor_values = 50, 100, 200, 400
alq_values = 0, 10000, 20000
```

### VFPINJ - Injection Wells

2-dimensional tables: BHP = f(THP, rate)

```ini
[VFP_INJ1]
table_number = 1
datum_depth = 2000.0
flow_type = WATER
thp_values = 1e6, 2e6, 3e6
rate_values = 0.001, 0.01, 0.1
```

### C++ API

```cpp
#include "VFPTables.hpp"

// Generate table from wellbore model
VFPTableGenerator generator;
generator.setTotalDepth(2000.0);
generator.setWellboreRadius(0.1);
generator.setOilProperties(800.0, 0.005);

auto table = generator.generateProductionTable(
    thp_range, rate_range, wct_range, gor_range, alq_range);

// Use for BHP calculation
double bhp = table.calculateBHP(thp, rate, wct, gor, alq);

// Inverse calculation
double rate = table.calculateFlowRate(bhp, thp, wct, gor, alq);

// Gas lift optimization
double optimal_alq = table.optimizeGasLift(bhp, thp, wct, gor, max_alq);
```

---

## Group Control

Hierarchical well group management with production/injection constraints.

### Group Hierarchy (GRUPTREE)

```ini
[GROUP1]
name = PLATFORM_A
parent = FIELD

[GROUP2]
name = REGION_NORTH
parent = PLATFORM_A
wells = PROD1, PROD2

[GROUP3]
name = REGION_SOUTH
parent = PLATFORM_A
wells = PROD3, PROD4
```

### Production Constraints (GCONPROD)

```ini
[GROUP_PROD_CONSTRAINTS]
group = PLATFORM_A
control_mode = ORAT                   # NONE, ORAT, WRAT, GRAT, LRAT, RESV
oil_rate_target = 0.1                 # m³/s
oil_rate_max = 0.12
water_cut_max = 0.95
gas_oil_ratio_max = 500
guide_rate_type = OIL                 # OIL, WATER, GAS, LIQ, COMB, FORM
```

### Injection Constraints (GCONINJE)

```ini
[GROUP_INJ_CONSTRAINTS]
group = INJECTORS
control_mode = VREP                   # RATE, RESV, REIN, VREP
phase = WATER
vrep_target = 1.0                     # 100% voidage replacement
surface_rate_max = 0.05
```

### C++ API

```cpp
#include "GroupControl.hpp"

GroupControlManager manager;

// Build hierarchy
manager.createGroup("PLATFORM", "FIELD");
manager.assignWellToGroup("PROD1", "PLATFORM");

// Set constraints
GroupProdConstraints constraints;
constraints.control_mode = GroupProdControlMode::ORAT;
constraints.oil_rate_target = 0.1;
manager.setProductionConstraints("PLATFORM", constraints);

// Apply controls
auto allocated_rates = manager.applyGroupControls(well_potentials);

// Check violations
auto violations = manager.checkConstraintViolations();
```

---

## Relative Permeability

Comprehensive relative permeability models with hysteresis and three-phase support.

### Two-Phase Models

| Model | Description |
|-------|-------------|
| Corey | Power-law kr = kr_max * S^n |
| Brooks-Corey | Based on pore-size distribution |
| Van Genuchten | Soil physics model |
| LET | Flexible 3-parameter model |
| Tabular | From SWOF/SGOF tables |

### Three-Phase Models

| Model | ECLIPSE | Description |
|-------|---------|-------------|
| Stone I | STONE1 | Probabilistic model |
| Stone II | STONE2 | Modified Stone method |
| Baker | BAKER | Linear interpolation |

### Hysteresis Models

| Model | Description |
|-------|-------------|
| Killough | Scanning curve with curvature parameter |
| Carlson | Parallel curve model |

### Configuration Example

```ini
[RELATIVE_PERMEABILITY]
model = TABULAR                       # COREY, BROOKS_COREY, VAN_GENUCHTEN, TABULAR
three_phase_model = STONE_I           # STONE_I, STONE_II, BAKER

# End-points
swc = 0.20                            # Connate water
sor = 0.25                            # Residual oil
sgc = 0.05                            # Critical gas

# Corey exponents (if using COREY model)
corey_nw = 2.5
corey_no = 2.0
corey_ng = 2.5

# Hysteresis
enable_hysteresis = true
hysteresis_model = KILLOUGH
killough_curvature = 0.1
land_parameter = 2.0
```

### End-Point Scaling

```ini
[ENDPOINT_SCALING]
enable = true
three_point = true                    # Use 3-point vs 2-point scaling

# Cell-specific scaling (from SWATINIT)
use_swatinit = true
swatinit_file = swatinit.dat
```

### C++ API

```cpp
#include "RelativePermeability.hpp"

// Create Corey model
auto corey = std::make_shared<CoreyRelPerm>();
corey->setWettingExponent(2.0);
corey->setNonWettingExponent(2.0);

// Create three-phase model
ThreePhaseRelPerm three_phase;
three_phase.setWaterOilTable(water_oil_table);
three_phase.setGasOilTable(gas_oil_table);
three_phase.setThreePhaseModel(ThreePhaseModel::STONE_I);

auto [kr_w, kr_o, kr_g] = three_phase.calculate(Sw, So, Sg);

// With hysteresis
HysteresisRelPerm hyst(drainage_curve, imbibition_curve, HysteresisModel::KILLOUGH);
hyst.updateState(current_Sw);
double kr = hyst.kr_wetting(current_Sw);
```

---

## Summary Output

ECLIPSE-compatible summary output for production data and analysis.

### Supported Keywords

#### Field Vectors (F*)
- FOPR, FWPR, FGPR, FLPR - Production rates
- FOPT, FWPT, FGPT - Cumulative production
- FWIR, FGIR - Injection rates
- FWIT, FGIT - Cumulative injection
- FWCT, FGOR - Ratios
- FPR - Average field pressure

#### Well Vectors (W*)
- WOPR, WWPR, WGPR - Well production rates
- WOPT, WWPT, WGPT - Well cumulative production
- WBHP, WTHP - Well pressures
- WWCT, WGOR - Well ratios

#### Group Vectors (G*)
- GOPR, GWPR, GGPR - Group production rates
- GOPT, GWPT, GGPT - Group cumulative

#### Block Vectors (B*)
- BPR - Block pressure
- BSWAT, BSGAS, BOSAT - Block saturations

#### Aquifer Vectors (A*)
- AAQP - Aquifer pressure
- AAQR - Aquifer influx rate
- AAQT - Aquifer cumulative influx

### Configuration

```ini
[SUMMARY]
# Field vectors
FOPR FWPR FGPR FOPT FWPT FGPT FPR FWCT FGOR

# Well vectors (for all wells)
WOPR WWPR WGPR WBHP WWCT

# Specific wells
WOPR PROD1
WBHP PROD1

# Block data
BPR 10 10 5

# Aquifer data
AAQP 1
AAQR 1
AAQT 1
```

### C++ API

```cpp
#include "SummaryOutput.hpp"

SummaryOutput output("CASE_NAME");
output.setWells({"PROD1", "PROD2", "INJ1"});
output.addDefaultOutput();
output.initialize();

// Record data each timestep
output.beginStep(time, dt, report_step);
output.recordField(oil_rate, water_rate, gas_rate, 
                   oil_total, water_total, gas_total, avg_pressure);
output.recordWell("PROD1", oil_rate, water_rate, gas_rate, bhp, thp,
                  oil_total, water_total, gas_total);
output.recordAquifer(1, pressure, influx_rate, influx_total);
output.endStep();

// Finalize
output.finalize();
output.writeCSV("summary.csv");
```

---

## Configuration Examples

### Complete Waterflood Example

See `config/waterdrive_aquifer.config` for a complete example demonstrating:
- Carter-Tracy and Fetkovich aquifers
- Multiple producers with group control
- Black oil PVT
- Summary output configuration

### Running Simulations

```bash
# Run with configuration file
mpirun -np 4 fsrm -c config/waterdrive_aquifer.config

# Override parameters
mpirun -np 4 fsrm -c config/waterdrive_aquifer.config \
    --dt-max 86400 \
    --output-frequency 5
```

---

## References

1. Schlumberger, "ECLIPSE Reference Manual"
2. Carter, R.D. and Tracy, G.W., "An Improved Method for Calculating Water Influx", Trans. AIME (1960)
3. Fetkovich, M.J., "A Simplified Approach to Water Influx Calculations", JPT (1971)
4. Stone, H.L., "Probability Model for Estimating Three-Phase Relative Permeability", JPT (1970)
5. Killough, J.E., "Reservoir Simulation with History-Dependent Saturation Functions", SPE J. (1976)
