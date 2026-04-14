# Volcano Benchmark Suite for FSRM

This directory contains verification benchmarks for FSRM's volcanic modeling capabilities.

## Overview

The volcano benchmark suite validates multi-physics volcanic simulations against:
- Analytical solutions (Mogi deformation, Woods column model)
- Laboratory experiments (PDC flume tests, lava rheology)
- Field observations (Pinatubo, Kilauea, Montserrat)
- Community model intercomparisons (IAVCEI benchmarks)

## Benchmark Categories

### 1. Eruption Column Dynamics (volcano_plinian_column.config)
- **Purpose**: Verify buoyant eruption column model
- **Model**: Woods (1988) integral column model
- **Verification**:
  - Column height vs. mass eruption rate
  - Air entrainment coefficient
  - Neutral buoyancy level
  - Column collapse criteria
- **Runtime**: ~5 minutes

### 2. Volcanic Deformation (volcano_mogi_deformation.config)
- **Purpose**: Verify deformation source models
- **Model**: Mogi (1958) point source
- **Verification**:
  - Radial and vertical displacement
  - Tilt components
  - Distance decay (1/r²)
- **Runtime**: ~1 minute (static)

### 3. Lava Flow Propagation (volcano_lava_flow.config)
- **Purpose**: Verify lava flow spreading and cooling
- **Model**: Harris & Rowland (2001) FLOWGO-style
- **Verification**:
  - Flow length and width
  - Temperature decay
  - Crust formation
  - Rheological evolution
- **Runtime**: ~15 minutes

### 4. Pyroclastic Density Currents (volcano_pdc_runout.config)
- **Purpose**: Verify PDC dynamics and runout
- **Model**: Shallow water + suspension equations
- **Verification**:
  - Runout distance
  - Dynamic pressure
  - Deposit thickness
  - Arrival times
- **Runtime**: ~10 minutes

## Usage

### Running Individual Benchmarks

```bash
# Run specific benchmark
./fsrm -c benchmarks/volcano_plinian_column.config

# Run with verification
./fsrm -c benchmarks/volcano_mogi_deformation.config --verify

# Run with visualization
./fsrm -c benchmarks/volcano_lava_flow.config --output output/lava_test
```

### Running Full Suite

```bash
# Run all volcano benchmarks
./run_volcano_suite.sh

# With verification report
./run_volcano_suite.sh --verify --report
```

### Verification Script

```bash
# Compare with reference solutions
python compare_volcano_benchmarks.py volcano_plinian_column
python compare_volcano_benchmarks.py volcano_mogi_deformation
```

## Success Criteria

| Benchmark | Metric | Tolerance |
|-----------|--------|-----------|
| Column Height | Woods (1988) model | ±10% |
| Mogi Displacement | Analytical solution | ±1% |
| Lava Flow Length | Field calibration | ±15% |
| PDC Runout | Lab experiments | ±15% |
| Dynamic Pressure | Field observations | ±30% |

## Benchmark Details

### Plinian Column Benchmark

Tests eruption column dynamics for a VEI-5 scenario:
- Mass eruption rate: 5×10⁷ kg/s
- Exit velocity: 200 m/s
- Exit temperature: 850°C
- Expected column height: ~25 km

Key physics verified:
- Jet thrust phase
- Buoyancy-driven rise
- Air entrainment (α ~ 0.09)
- Thermal equilibration

### Mogi Deformation Benchmark

Tests elastic half-space deformation from spherical source:
- Source depth: 5 km
- Volume change: 10⁶ m³
- Shear modulus: 30 GPa
- Poisson ratio: 0.25

Analytical solution (Mogi, 1958):
```
u_z = (1-ν)/G × ΔV × d / (r² + d²)^(3/2)
u_r = (1-ν)/G × ΔV × r / (r² + d²)^(3/2)
```

### Lava Flow Benchmark

Tests basaltic lava flow on uniform slope:
- Effusion rate: 20 m³/s
- Temperature: 1150°C
- Slope: 5°
- Duration: 24 hours

Key physics verified:
- Gravity-driven spreading
- Temperature-dependent viscosity
- Surface crust formation
- Cooling by radiation and convection

### PDC Runout Benchmark

Tests column-collapse PDC on conical volcano:
- Collapse height: 5 km
- Mass flux: 10⁸ kg/s
- Initial temperature: 500°C
- Particle concentration: 1%

Key physics verified:
- Density stratification
- Air entrainment
- Particle settling
- Topographic channeling

## Expected Results

### Plinian Column

```
Column parameters:
  Maximum height: 24-26 km
  NBL height: 18-22 km
  Entrainment coefficient: 0.08-0.10
  Column type: Sustained (no collapse)
```

### Mogi Deformation

```
Surface displacements at key points:
  r=0:    u_z = 50 mm,   u_r = 0 mm
  r=5km:  u_z = 17.7 mm, u_r = 17.7 mm  
  r=10km: u_z = 4.5 mm,  u_r = 8.9 mm
```

### Lava Flow

```
Flow dimensions after 24 hours:
  Length: 7-9 km
  Width: 150-250 m
  Area: 1.2-1.8 km²
  Volume: ~1.7 Mm³
```

### PDC Runout

```
Runout distances:
  Downslope: 10-14 km
  Cross-slope: 6-10 km
  
Dynamic pressure at 5 km: 10-20 kPa
Arrival time at 5 km: 100-150 s
```

## Reference Solutions

Reference solutions are available in:
```
benchmarks/reference/
├── volcano_column/
│   ├── woods1988_analytical.dat
│   └── iavcei_intercomparison.dat
├── volcano_deformation/
│   ├── mogi1958_analytical.dat
│   └── mctigue1987_finite.dat
├── volcano_lava/
│   ├── flowgo_calibration.dat
│   └── kilauea2018_observed.dat
└── volcano_pdc/
    ├── titan2d_comparison.dat
    └── montserrat_observed.dat
```

## Continuous Integration

Volcano benchmarks are run:
- On every commit to volcano-related files
- Weekly full suite verification
- Before release builds

CI checks:
- All benchmarks complete without errors
- Results within specified tolerances
- No performance regression (>20% slower)

## Adding New Benchmarks

To add a new volcano benchmark:

1. Create config file: `volcano_<name>.config`
2. Add reference solution if available
3. Update verification script
4. Document expected results
5. Run and verify

## Known Limitations

### Column Model
- 1D integral model (no 3D effects)
- Simplified particle settling
- No wind effects currently

### Deformation Model
- Half-space assumption
- Elastic only (no viscoelastic)
- Single source at a time verified

### Lava Flow Model
- No lava tube formation
- Simplified breakout behavior
- 2D only (no vertical structure)

### PDC Model
- Depth-averaged equations
- Simplified particle physics
- No erosion currently

## Scientific References

### Eruption Columns
1. Woods, A.W. (1988). The fluid dynamics and thermodynamics of eruption columns. *Bull. Volcanol.*, 50, 169-193.
2. Bursik, M. (2001). Effect of wind on the rise height of volcanic plumes. *Geophys. Res. Lett.*, 28, 3621-3624.
3. Carazzo, G., et al. (2008). On the rise of turbulent plumes: Quantitative effects of variable entrainment. *J. Geophys. Res.*, 113, B09201.

### Volcanic Deformation
1. Mogi, K. (1958). Relations between the eruptions of various volcanoes and the deformation of the ground surfaces around them. *Bull. Earthq. Res. Inst.*, 36, 99-134.
2. McTigue, D.F. (1987). Elastic stress and deformation near a finite spherical magma body. *J. Geophys. Res.*, 92, 12931-12940.
3. Segall, P. (2010). *Earthquake and Volcano Deformation*. Princeton University Press.

### Lava Flows
1. Harris, A.J.L. & Rowland, S.K. (2001). FLOWGO: A kinematic thermo-rheological model for lava flowing in a channel. *Bull. Volcanol.*, 63, 20-44.
2. Crisci, G.M., et al. (2004). Predicting the Etna 2001 lava flow paths by a cellular automata model. *Environ. Model. Softw.*, 19, 879-890.

### Pyroclastic Density Currents
1. Denlinger, R.P. & Iverson, R.M. (2001). Flow of variably fluidized granular masses across three-dimensional terrain. *J. Geophys. Res.*, 106, 553-566.
2. Patra, A.K., et al. (2005). Parallel adaptive numerical simulation of dry avalanches over natural terrain. *J. Volcanol. Geotherm. Res.*, 139, 1-21.

## Contact

For questions about volcano benchmarks in FSRM:
- Implementation: See `include/VolcanoModel.hpp`
- Documentation: See `docs/VOLCANO_MODELING.md`
- Issues: GitHub issue tracker
