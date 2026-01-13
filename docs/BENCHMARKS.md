# FSRM Benchmarks

This document describes all available benchmarks for validating FSRM's accuracy and performance.

## Overview

FSRM includes comprehensive benchmark suites covering:
- SPE (Society of Petroleum Engineers) comparative solution projects
- SCEC (Southern California Earthquake Center) dynamic rupture benchmarks
- Analytical solutions for verification
- Performance and scalability tests

## SPE Benchmarks

### SPE1 - First SPE Comparative Solution Project

**Description**: Three-phase black oil simulation with gas injection.

**Features**:
- 10x10x3 Cartesian grid
- Black oil model (oil, gas, water)
- Vertical injector and producer
- Initial oil with dissolved gas

**Configuration**: `config/spe1_benchmark.config`

**Executable**: `examples/spe1`

**Run**:
```bash
cd build/examples
./spe1 ../config/spe1_benchmark.config
```

**Reference**: Odeh, A.S. (1981). "Comparison of Solutions to a Three-Dimensional Black-Oil Reservoir Simulation Problem." JPT.

### SPE3 - Third SPE Comparative Solution Project

**Description**: Gas cycling in a gas condensate reservoir.

**Configuration**: `config/spe3_benchmark.config`

**Run**:
```bash
mpirun -np 4 fsrm -c config/spe3_benchmark.config
```

### SPE9 - Ninth SPE Comparative Solution Project

**Description**: Black oil simulation with complex well controls.

**Configuration**: `config/spe9_benchmark.config`

### SPE10 - Tenth SPE Comparative Solution Project

**Description**: Large-scale heterogeneous reservoir (60x220x85 cells).

**Features**:
- Highly heterogeneous permeability field (Tarbert and Upper Ness formations)
- Two-phase flow (oil and water)
- Waterflooding scenario
- 1.1 million cells (full model)

**Configuration**: `config/spe10_benchmark.config`

**Note**: SPE10 data files must be downloaded separately from SPE website.

## SCEC Dynamic Rupture Benchmarks

FSRM implements the SCEC/USGS dynamic rupture benchmark suite for validating earthquake physics.

### TPV5 - Homogeneous Halfspace with Slip-Weakening Friction

**Description**: Vertical strike-slip fault in homogeneous elastic halfspace.

**Features**:
- Linear slip-weakening friction
- Spontaneous rupture nucleation
- Supershear transition

**Configuration**: `config/scec_tpv5.config`

**Executable**: `examples/scec_tpv5`

**Run**:
```bash
cd build/examples
mpirun -np 4 ./scec_tpv5
```

**Reference**: [SCEC TPV5 Description](http://scecdata.usc.edu/cvws/tpv5docs.html)

### TPV10 - Dipping Fault (60°)

**Description**: Strike-slip fault dipping at 60°.

**Features**:
- Complex geometry
- Dip-parallel and dip-perpendicular motion
- Surface breaking fault

**Configuration**: `config/scec_tpv10.config`

**Executable**: `examples/scec_tpv10`

### TPV16 - Heterogeneous Initial Stress

**Description**: Strike-slip fault with spatially variable initial stress.

**Features**:
- Heterogeneous stress field
- Multiple rupture nucleation sites
- Bilateral rupture propagation

**Configuration**: `benchmarks/scec_tpv16.config`

**Executable**: `examples/scec_tpv16`

### Complete SCEC Suite

All SCEC benchmarks are organized in the `benchmarks/` directory:

```bash
# Run complete suite
cd benchmarks
./run_scec_suite.sh --all

# Run specific benchmark
./run_scec_suite.sh --tpv 5

# Run with verification against reference
./run_scec_suite.sh --all --verify
```

**Available Benchmarks**:
- TPV5: Basic slip-weakening
- TPV10: Dipping fault
- TPV13: Branched fault with plasticity
- TPV16: Heterogeneous stress
- TPV34: Thermal pressurization
- TPV101: Rate-and-state friction (aging law)
- TPV104: Rate-and-state friction (strong velocity weakening)

See `benchmarks/SCEC_BENCHMARKS_README.md` for complete documentation.

## Analytical Verification Tests

### Single-Phase Flow

#### Theis Solution
**Description**: Radial flow from a point source in infinite domain.

**Configuration**: Available in test suite

**Verification**: Analytical solution for pressure transient.

#### Buckley-Leverett Solution
**Description**: 1D two-phase waterflooding with sharp front.

**Configuration**: `config/buckley_leverett_2d.config`

**Verification**: Analytical shock front position and saturation profile.

### Geomechanics

#### Terzaghi Consolidation
**Description**: 1D consolidation under constant load.

**Verification**: Analytical solution for pressure dissipation and settlement.

#### Mandel-Cryer Effect
**Description**: 2D poroelastic problem with initial pressure increase.

**Verification**: Analytical solution for pressure and displacement fields.

### Wave Propagation

#### LOH.1 - Layer Over Halfspace
**Description**: Seismic wave propagation through layered medium.

**Features**:
- Vertical velocity contrast
- Reflected and refracted waves
- Surface waves

**Configuration**: `config/scec_loh1.config`

**Executable**: `examples/scec_loh1`

**Reference**: SCEC-USGS verification exercises.

## Performance Benchmarks

### Scalability Tests

Test parallel efficiency with increasing processor counts.

**Configuration**: `config/default.config` with varying grid sizes

**Run**:
```bash
# Weak scaling (constant cells per processor)
for np in 1 2 4 8 16 32; do
  mpirun -np $np fsrm -c config/scaling_weak.config
done

# Strong scaling (constant total cells)
for np in 1 2 4 8 16 32; do
  mpirun -np $np fsrm -c config/scaling_strong.config
done
```

### GPU Performance

Compare CPU vs GPU execution.

**Run**:
```bash
# CPU baseline
mpirun -np 8 fsrm -c config/gpu_test.config

# Single GPU
fsrm -c config/gpu_test.config --use-gpu

# Multiple GPUs
mpirun -np 4 fsrm -c config/gpu_test.config --use-gpu

# Benchmark mode
fsrm -c config/gpu_test.config --benchmark-gpu
```

**Expected Speedups**:
- Single-phase flow: 10-20x
- Elastodynamics: 30-50x
- Poroelastodynamics: 20-40x
- Black oil: 15-25x

See `docs/DEPLOYMENT.md` for GPU setup and `README.md` for performance details.

## Application Examples

### Hydraulic Fracturing

**Configuration**: `config/hydraulic_fracturing.config`

**Features**:
- PKN/KGD fracture models
- Proppant transport
- Leak-off and pressure decline

### Induced Seismicity

**Configuration**: `config/induced_seismicity.config`

**Features**:
- Fluid injection into faulted reservoir
- Rate-and-state friction
- Dynamic rupture triggering

### Enhanced Geothermal System

**Configuration**: `config/geothermal.config`

**Features**:
- High temperature (300°C)
- Low permeability granite
- Thermal-hydraulic-mechanical coupling

### CO2 Storage

**Configuration**: `config/co2_storage.config`

**Features**:
- Deep saline aquifer
- CO2 phase behavior
- Capillary trapping
- Long-term migration

### Shale Reservoir

**Configuration**: `config/shale_reservoir.config`

**Features**:
- Ultra-low permeability (~100 nD)
- Multi-stage hydraulic fracturing
- Gas desorption
- Long horizontal wells

## Running Benchmarks

### Quick Test
```bash
# Single benchmark
mpirun -np 4 fsrm -c config/spe1_benchmark.config

# With verification
mpirun -np 4 fsrm -c config/spe1_benchmark.config --verify
```

### Complete Test Suite
```bash
# All unit tests
cd build
make test

# All SCEC benchmarks
cd benchmarks
./run_scec_suite.sh --all --verify

# Method of Manufactured Solutions tests
cd tests
./run_mms_tests.sh
```

### Visualization
```bash
# VTK output (for ParaView)
paraview output/*.vtu

# HDF5 output (more efficient)
python scripts/hdf5_to_xdmf.py output/
paraview output/solution.xdmf

# Plots (automatically generated)
ls output/*.png
```

## Benchmark Results

### Typical Accuracy

| Benchmark | L2 Error | Convergence Rate |
|-----------|----------|------------------|
| SPE1 | < 1% | N/A (discrete) |
| SCEC TPV5 | < 2% | N/A (reference) |
| Theis | < 0.01% | O(h²) |
| Terzaghi | < 0.1% | O(h²) |
| LOH.1 | < 5% | O(h²) |

### Typical Performance (CPU, 32 cores)

| Problem Size | Cells | Time/Step | Memory |
|--------------|-------|-----------|--------|
| Small | 1,000 | 0.01 s | 10 MB |
| Medium | 100,000 | 1 s | 1 GB |
| Large | 1,000,000 | 20 s | 10 GB |
| Very Large | 10,000,000 | 5 min | 100 GB |

### GPU Speedup (vs 32 CPU cores)

| Problem Size | Single GPU | 4 GPUs |
|--------------|------------|--------|
| Small (1k) | 1-2x | - |
| Medium (100k) | 10-15x | - |
| Large (1M) | 20-30x | 18-25x |
| Very Large (10M) | 30-40x | 60-100x |

## References

### SPE Benchmarks
1. Odeh, A.S. (1981). "Comparison of Solutions to a Three-Dimensional Black-Oil Reservoir Simulation Problem." JPT, January 1981.
2. Christie, M.A. and Blunt, M.J. (2001). "Tenth SPE Comparative Solution Project: A Comparison of Upscaling Techniques." SPE Reservoir Eval. & Eng.

### SCEC Benchmarks
1. Harris, R.A., et al. (2009). "The SCEC/USGS Dynamic Earthquake Rupture Code Verification Exercise." Seismological Research Letters.
2. Harris, R.A., et al. (2018). "A Suite of Exercises for Verifying Dynamic Earthquake Rupture Codes." Seismological Research Letters.

### Analytical Solutions
1. Theis, C.V. (1935). "The relation between the lowering of the Piezometric surface and the rate and duration of discharge of a well using ground-water storage." Trans. AGU.
2. Buckley, S.E. and Leverett, M.C. (1942). "Mechanism of Fluid Displacement in Sands." Trans. AIME.
3. Mandel, J. (1953). "Consolidation Des Sols (Étude Mathématique)." Géotechnique.

## See Also

- [User Guide](USER_GUIDE.md) - Running simulations
- [Configuration Reference](CONFIGURATION.md) - All configuration options
- [Physics Models](PHYSICS_MODELS.md) - Mathematical formulations
- [Development Guide](DEVELOPMENT.md) - Building and testing
