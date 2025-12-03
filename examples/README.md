# FSRM Examples - Config-Driven Approach

## Overview

All FSRM examples now use a **unified config-driven approach**. Instead of having separate executables for each example, there is a single `simulator` executable that reads all simulation parameters from configuration files.

## Benefits

- **No recompilation needed**: Change any parameter by editing the config file
- **Easy parameter exploration**: Copy and modify config files for different scenarios
- **Reproducibility**: Complete simulation setup in a single, version-controllable text file
- **Standardization**: All examples follow the same structure and conventions
- **HDF5 default output**: More efficient storage and I/O compared to VTK

## Quick Start

### 1. Build the simulator

```bash
mkdir build && cd build
cmake ..
make -j
```

This creates the `simulator` executable in `build/examples/`.

### 2. Run an example

```bash
# Run with a specific config file
./simulator -c config/lefm_fracture_growth.config

# Run with MPI (parallel)
mpirun -np 4 ./simulator -c config/hydraulic_fracturing.config

# Run with custom config
./simulator -c /path/to/my_custom.config
```

### 3. Visualize results

**HDF5 outputs (default):**
```bash
# Convert to XDMF for ParaView
python scripts/hdf5_to_xdmf.py output/

# Open in ParaView
paraview output/solution.xdmf
```

**VTK outputs (legacy):**
```bash
# Set output_format = VTK in config file
paraview output/*.vtu
```

## Available Examples

All example configs are in the `config/` directory:

### Development Optimization Examples

#### `well_spacing_optimization.config` ⭐ **FEATURED**
- **Physics**: Well spacing + fracture placement + economic optimization
- **Features**: Multi-well interference, stress shadowing, NPV maximization
- **Runtime**: ~5-30 minutes
- **Use case**: Unconventional field development, completion design optimization

### Basic Flow Examples

#### `single_phase.config`
- **Physics**: Single-phase Darcy flow
- **Features**: Injection/production wells, basic pressure diffusion
- **Runtime**: ~1 minute
- **Use case**: Learning the basics, testing setup

#### `shale_reservoir.config`
- **Physics**: Single-phase flow in tight shale
- **Features**: Ultra-low permeability, long-term production
- **Runtime**: ~5 minutes
- **Use case**: Unconventional reservoir simulation

### Hydraulic Fracturing Examples

#### `hydraulic_fracturing.config`
- **Physics**: Coupled flow + geomechanics + fracture propagation
- **Features**: Proppant transport, fracture network
- **Runtime**: ~10 minutes
- **Use case**: Frac treatment design

#### `lefm_fracture_growth.config` ⭐ **IMPROVED**
- **Physics**: Linear Elastic Fracture Mechanics (LEFM)
- **Features**: Stress intensity factor K_I, toughness-dominated propagation
- **Runtime**: ~15 minutes
- **Use case**: Fracture mechanics validation, academic research
- **Documentation**: Extensive inline comments explaining LEFM theory

**Key features of LEFM example:**
- Stress intensity factor calculation: K_I = Y × σ × √(πa)
- Initiation criterion: K_I = K_Ic (1.0 MPa·m^0.5)
- Propagation: v = C × (K_I - K_Ic) / K_Ic
- Width profile from Sneddon solution
- Carter leak-off model
- Time-dependent fracture geometry evolution

### Seismicity Examples

#### `induced_seismicity.config`
- **Physics**: Poroelasticity + rate-state friction faults
- **Features**: Injection-induced earthquakes, Coulomb stress transfer
- **Runtime**: ~30 minutes
- **Use case**: Seismic hazard assessment

#### `static_triggered_seismicity.config`
- **Physics**: Static stress triggering on fault network
- **Features**: Multiple faults, aftershock sequences
- **Runtime**: ~20 minutes

#### `coulomb_stress_parallel_faults.config` ⭐ **NEW**
- **Physics**: Coulomb stress transfer between parallel faults
- **Features**: ΔCFF calculation, stress shadows, triggering lobes
- **Runtime**: ~1 minute (standalone executable)
- **Use case**: Fault interaction analysis, earthquake triggering studies
- **Output formats**: Snapshot images (PNG/PPM) and ASCII data files

**Key features:**
- Calculate Coulomb stress change (ΔCFF) from slip on source fault
- Compute stress distribution on parallel receiver fault  
- Visualize stress transfer patterns (positive lobes, stress shadows)
- Output as snapshot images for quick visualization
- Output as ASCII files for external analysis (Gnuplot, MATLAB, Python)
- Command-line configuration for parameter sweeps

**Usage:**
```bash
# Run with default parameters
./coulomb_stress_parallel_faults

# Custom slip and separation
./coulomb_stress_parallel_faults -slip 2.0 -sep 3000

# ASCII output only
./coulomb_stress_parallel_faults -output_format ascii

# Snapshot images only
./coulomb_stress_parallel_faults -output_format snapshot

# Both formats (default)
./coulomb_stress_parallel_faults -output_format both
```

### Wave Propagation Examples

#### `wave_propagation.config`
- **Physics**: Elastodynamic wave equation
- **Features**: Seismic source, wave speed visualization
- **Runtime**: ~5 minutes
- **Use case**: Microseismic monitoring simulation

#### `elastodynamic_waves.config`
- **Physics**: Full elastodynamics with advanced source models
- **Features**: P-waves, S-waves, surface waves
- **Runtime**: ~10 minutes

#### `poroelastodynamic_waves.config`
- **Physics**: Coupled poroelastodynamics (Biot waves)
- **Features**: Slow P-wave, fast P-wave, S-wave
- **Runtime**: ~15 minutes

### Energy Applications

#### `geothermal.config`
- **Physics**: Coupled flow + heat transfer
- **Features**: Temperature-dependent properties, thermal recovery
- **Runtime**: ~20 minutes
- **Use case**: Geothermal reservoir management

#### `co2_storage.config`
- **Physics**: Two-phase CO2-brine flow + geomechanics
- **Features**: Capillary pressure, buoyancy-driven migration
- **Runtime**: ~25 minutes
- **Use case**: Carbon sequestration

### SPE Benchmarks (Complete Suite)

#### `spe1_benchmark.config`
- **Physics**: Black oil model (3-phase)
- **Features**: 10×10×3 grid, gas injection, industry standard
- **Runtime**: ~3 minutes
- **Use case**: Basic black oil validation

#### `spe2_benchmark.config`
- **Physics**: Three-phase coning
- **Features**: Radial grid, gas cap, aquifer, critical rate analysis
- **Runtime**: ~10 minutes
- **Use case**: Coning validation

#### `spe3_benchmark.config`
- **Physics**: Compositional (gas cycling)
- **Features**: 9×9×4 grid, 4-component model, pressure maintenance
- **Runtime**: ~15 minutes
- **Use case**: Compositional validation

#### `spe4_benchmark.config`
- **Physics**: Equation of State testing
- **Features**: Peng-Robinson EOS, multi-component flash, phase envelope
- **Runtime**: ~5 minutes
- **Use case**: EOS implementation validation

#### `spe5_benchmark.config`
- **Physics**: Miscible CO2 flooding
- **Features**: 7×7×3 grid, 6-component, first-contact miscibility
- **Runtime**: ~20 minutes
- **Use case**: Miscible flood validation

#### `spe6_benchmark.config`
- **Physics**: Dual-porosity / dual-permeability
- **Features**: Matrix-fracture transfer, Kazemi shape factor
- **Runtime**: ~10 minutes
- **Use case**: Naturally fractured reservoir validation

#### `spe7_benchmark.config`
- **Physics**: Horizontal well modeling
- **Features**: 9×9×6 grid, horizontal trajectory, Peaceman PI
- **Runtime**: ~15 minutes
- **Use case**: Horizontal well validation

#### `spe8_benchmark.config`
- **Physics**: Gridding and upscaling
- **Features**: 40×40×3 fine grid, channel heterogeneity, 5-spot
- **Runtime**: ~20 minutes
- **Use case**: Grid sensitivity analysis

#### `spe9_benchmark.config`
- **Physics**: North Sea reservoir model
- **Features**: 24×25×15 grid (9,000 cells), multiple wells
- **Runtime**: ~30 minutes
- **Use case**: Large-scale black oil validation

#### `spe10_benchmark.config`
- **Physics**: Large heterogeneous problem
- **Features**: 60×220×85 grid (1.1M cells), extreme permeability contrast
- **Runtime**: ~2-4 hours (parallel recommended)
- **Use case**: Scalability and upscaling validation

#### `spe11_benchmark.config` / `spe11a_benchmark.config` / `spe11c_benchmark.config`
- **Physics**: CO2 geological storage
- **Features**: CO2-brine two-phase flow, dissolution, density-driven convection, trapping mechanisms
- **Case A** (`spe11a`): Small 2D test case (2.8 km × 1.2 km) - quick validation
- **Case B** (`spe11`): Standard 2D case (8.4 km × 1.2 km) - default benchmark
- **Case C** (`spe11c`): Full 3D model (8.4 km × 5 km × 1.2 km) - comprehensive testing
- **Duration**: 1000-year simulation with 50-year injection period
- **Runtime**: Case A: ~1 hour, Case B: ~4-8 hours, Case C: ~24+ hours (parallel required)
- **Reference**: Nordbotten et al. SPE Journal (2024), https://spe11-group.github.io/
- **Use case**: Carbon capture and storage (CCS) simulation validation

### SCEC Benchmarks (Complete CVWS Suite)

The complete SCEC Code Verification Workshop Series (CVWS) benchmark suite is available.
Reference: https://strike.scec.org/cvws/benchmark_descriptions.html

#### TPV3-9: Basic Dynamic Rupture

| Config | Physics | Key Features | Runtime |
|--------|---------|--------------|---------|
| `scec_tpv3.config` | 2D in-plane strike-slip | Linear slip-weakening, forced nucleation | ~10 min |
| `scec_tpv4.config` | 3D strike-slip rate-state | VW/VS zones, aging law | ~30 min |
| `scec_tpv5.config` | 3D vertical strike-slip | Heterogeneous stress, slip-weakening | ~20 min |
| `scec_tpv6.config` | 3D bilateral rupture | Symmetric propagation | ~25 min |
| `scec_tpv7.config` | VS barrier | Velocity-strengthening arrest | ~25 min |
| `scec_tpv8.config` | Rate-dependent friction | VW/VS rate dependence | ~30 min |
| `scec_tpv9.config` | Depth-dependent stress | Lithostatic gradient | ~25 min |

#### TPV10-15: Fault Geometry

| Config | Physics | Key Features | Runtime |
|--------|---------|--------------|---------|
| `scec_tpv10.config` | Branching fault | Main + 30° branch | ~45 min |
| `scec_tpv11.config` | Normal faulting | Dip-slip mechanism | ~30 min |
| `scec_tpv12.config` | 60° dipping (SW) | Slip-weakening friction | ~30 min |
| `scec_tpv13.config` | 60° dipping (RS) | Rate-state friction | ~35 min |
| `scec_tpv14.config` | Branch main-dominant | Stress favors main fault | ~45 min |
| `scec_tpv15.config` | Branch branch-dominant | Stress favors branch | ~45 min |

#### TPV16-23: Advanced Physics

| Config | Physics | Key Features | Runtime |
|--------|---------|--------------|---------|
| `scec_tpv16.config` | Rough fault | Self-similar roughness | ~45 min |
| `scec_tpv17.config` | Supershear transition | Subshear → supershear | ~30 min |
| `scec_tpv18.config` | Layered medium | Multiple velocity layers | ~35 min |
| `scec_tpv19.config` | Self-healing pulse | Strong rate-weakening | ~40 min |
| `scec_tpv20.config` | Shallow dip-slip (SW) | Surface rupture, SW friction | ~35 min |
| `scec_tpv21.config` | Shallow dip-slip (RS) | Surface rupture, RS friction | ~40 min |
| `scec_tpv22.config` | Bilateral (SW) | Centered nucleation, SW | ~30 min |
| `scec_tpv23.config` | Bilateral (RS) | Centered nucleation, RS | ~35 min |

#### TPV24-31: Japan & Complex Scenarios

| Config | Physics | Key Features | Runtime |
|--------|---------|--------------|---------|
| `scec_tpv24.config` | Japan homogeneous | Strong velocity weakening | ~1 hr |
| `scec_tpv25.config` | Japan heterogeneous | Flash heating, stress variability | ~1 hr |
| `scec_tpv26.config` | Thermal pressurization | Pore fluid heating | ~45 min |
| `scec_tpv27.config` | TP heterogeneous | Thermal press. + hetero. stress | ~50 min |
| `scec_tpv28.config` | T-junction | Perpendicular fault intersection | ~40 min |
| `scec_tpv29.config` | Topography | Non-planar free surface | ~50 min |
| `scec_tpv30.config` | Heterogeneous stress | Stochastic stress field | ~35 min |
| `scec_tpv31.config` | Friction heterogeneity | Variable friction parameters | ~35 min |

#### TPV32-37: Material Heterogeneity

| Config | Physics | Key Features | Runtime |
|--------|---------|--------------|---------|
| `scec_tpv32.config` | Normal fault (SW) | Free surface, slip-weakening | ~35 min |
| `scec_tpv33.config` | Normal fault (RS) | Free surface, rate-state | ~40 min |
| `scec_tpv34.config` | Off-fault plasticity | Drucker-Prager, homogeneous | ~45 min |
| `scec_tpv35.config` | Plasticity hetero. | Drucker-Prager, heterogeneous | ~50 min |
| `scec_tpv36.config` | Bimaterial (SW) | Material contrast, SW friction | ~40 min |
| `scec_tpv37.config` | Bimaterial (RS) | Material contrast, RS friction | ~45 min |

#### TPV100+ Series: Inelastic Benchmarks

| Config | Physics | Key Features | Runtime |
|--------|---------|--------------|---------|
| `scec_tpv101.config` | Viscoelastic plasticity | Q-factor attenuation, Drucker-Prager | ~1 hr |
| `scec_tpv102.config` | Enhanced plasticity | Depth-dependent cohesion, softening | ~1 hr |
| `scec_tpv103.config` | RS + plasticity | Rate-state friction + off-fault plastic | ~1 hr |
| `scec_tpv104.config` | Viscoplasticity | Perzyna model, rate-dependent plastic | ~1.5 hr |
| `scec_tpv105.config` | Damage mechanics | Lyakhovsky damage, stiffness degradation | ~1 hr |

#### LOH Series: Wave Propagation

| Config | Physics | Key Features | Runtime |
|--------|---------|--------------|---------|
| `scec_loh1.config` | Layer over halfspace | 1 km layer, explosion source | ~15 min |
| `scec_loh2.config` | Heterogeneous layer | Double-couple, Love waves | ~20 min |
| `scec_loh3.config` | Attenuation | Viscoelastic Q, anelastic effects | ~25 min |
| `scec_loh4.config` | Multiple layers | Crustal velocity model | ~30 min |

### Optimization Examples

#### `well_spacing_optimization.config` ⭐ **NEW**
- **Physics**: Well spacing and hydraulic fracture placement optimization
- **Features**: Multi-well interference, stress shadowing, economic NPV analysis
- **Runtime**: ~5-30 minutes (depending on grid search resolution)
- **Use case**: Unconventional reservoir development planning, completion optimization

**Key capabilities:**
- Grid search optimization over well spacing, stage spacing, and fracture geometry
- Well-to-well interference modeling (pressure depletion, frac hits)
- Inter-fracture stress shadowing effects
- Production forecasting with decline curve analysis
- Full economic analysis: NPV, IRR, payout, cost per BOE
- Response surface generation for visualization
- Sensitivity analysis and Monte Carlo uncertainty

**Optimization variables:**
- Inter-well spacing: 150-450 m (500-1500 ft)
- Stage spacing: 30-90 m (100-300 ft)
- Cluster spacing: 6-18 m (20-60 ft)
- Fracture half-length: 75-225 m (250-750 ft)
- Number of wells per drilling unit

**Usage:**
```bash
# Full optimization
./well_spacing_optimization -c config/well_spacing_optimization.config

# Quick mode (coarser grid search)
./well_spacing_optimization -c config/well_spacing_optimization.config -quick

# Evaluate specific scenario
./well_spacing_optimization -c config/well_spacing_optimization.config \
    -well_spacing 250.0 -stage_spacing 50.0 -frac_length 150.0 -eval_only

# Parallel execution
mpirun -np 4 ./well_spacing_optimization -c config/well_spacing_optimization.config
```

**Output files:**
- `optimization_summary.txt` - Human-readable results with recommendations
- `optimization_results.csv` - All scenarios for spreadsheet analysis
- `response_surface.dat` - Data for contour plots
- `plot_response_surface.gp` - Gnuplot visualization script

**Typical optimal ranges (Permian Basin-like):**
- Well spacing: 200-300 m (660-1000 ft)
- Stage spacing: 45-60 m (150-200 ft)
- Cluster spacing: 9-15 m (30-50 ft)
- Fracture half-length: 120-180 m (400-600 ft)

### GPU-Accelerated Examples

#### `gpu_elastodynamics.config`
- **Physics**: GPU-accelerated wave propagation
- **Features**: CUDA kernels for large-scale problems
- **Requirements**: CUDA-capable GPU

#### `gpu_hybrid_large_scale.config`
- **Physics**: Hybrid CPU/GPU computation
- **Features**: Dynamic load balancing
- **Requirements**: Multi-GPU system

## Creating Custom Simulations

### Option 1: Copy and modify existing config

```bash
# Copy an example config
cp config/lefm_fracture_growth.config my_simulation.config

# Edit parameters
nano my_simulation.config

# Run your custom simulation
./simulator -c my_simulation.config
```

### Option 2: Generate template

```bash
# Generate complete template with all options
./simulator --generate-template my_template.config

# Edit and run
nano my_template.config
./simulator -c my_template.config
```

## Configuration File Structure

A typical config file has these sections:

```ini
[SIMULATION]
# Time control, physics models, solver settings
start_time = 0.0
end_time = 3600.0
output_format = HDF5  # HDF5 (default) or VTK (secondary)

[GRID]
# Domain dimensions and mesh resolution
nx = 80
ny = 80
nz = 40

[ROCK]
# Rock properties (porosity, permeability, elasticity)
youngs_modulus = 25.0e9
fracture_toughness = 1.0e6  # For LEFM

[FLUID]
# Fluid properties (density, viscosity)
viscosity = 0.001

[WELL1], [WELL2], ...
# Well locations and controls
type = INJECTOR
target_value = 0.05  # 50 L/s

[FRACTURE1], [FRACTURE2], ...
# Fracture properties and propagation criteria
toughness = 1.0e6
enable_propagation = true

[BC1], [BC2], ...
# Boundary conditions
type = DIRICHLET
field = PRESSURE
value = 25.0e6

[IC1], [IC2], ...
# Initial conditions
field = PRESSURE
value = 25.0e6
```

See `config/complete_template.config` for all available options.

## Output Formats

### HDF5 (Default) ✅

**Advantages:**
- Efficient binary format
- Parallel I/O support
- Large file handling (>2GB)
- Self-describing metadata
- Industry standard

**Visualization:**
```bash
# Generate XDMF wrapper
python scripts/hdf5_to_xdmf.py output/

# Open in ParaView
paraview output/solution.xdmf

# Or use Visit, yt, h5py
```

### VTK (Secondary)

**When to use:**
- Legacy workflows
- Small simulations
- Direct ParaView compatibility

**Usage:**
Set in config file:
```ini
output_format = VTK
```

Then open directly:
```bash
paraview output/*.vtu
```

### ECLIPSE (Specialized)

For reservoir engineering workflows:
```ini
output_format = ECLIPSE
```

## Common Workflows

### Parameter Study

```bash
# Create parameter variations
for k_ic in 0.5 1.0 1.5 2.0; do
    sed "s/fracture_toughness = 1.0e6/fracture_toughness = ${k_ic}e6/" \
        config/lefm_fracture_growth.config > config/lefm_kic_${k_ic}.config
    ./simulator -c config/lefm_kic_${k_ic}.config
done
```

### Restart from Checkpoint

```bash
# Enable checkpointing in config
enable_checkpointing = true
checkpoint_frequency = 100

# Run initial simulation
./simulator -c my_sim.config

# Restart from checkpoint
./simulator -c my_sim.config -restart output/checkpoint_1000.h5
```

### Parallel Runs

```bash
# Strong scaling study
for np in 1 2 4 8 16; do
    mpirun -np $np ./simulator -c config/large_scale.config
done
```

## Troubleshooting

### "Cannot open config file"
```bash
# Use absolute path or correct relative path
./simulator -c $(pwd)/config/lefm_fracture_growth.config
```

### "HDF5 library not found"
```bash
# Rebuild with HDF5 support
cmake -DENABLE_HDF5=ON ..
make
```

### "Simulation diverged"
- Reduce timestep: `dt_initial = 0.1` (smaller value)
- Tighten tolerances: `rtol = 1e-8`
- Check boundary conditions consistency
- Enable adaptive timestepping: `adaptive_timestepping = true`

### "VTK output requested but not available"
- VTK is still supported, just secondary
- Check that VTK libraries are installed
- Use HDF5 instead for better performance

## Example Output

After running `./simulator -c config/lefm_fracture_growth.config`:

```
output/
├── lefm/
│   ├── solution_0000.h5      # Initial state
│   ├── solution_0060.h5      # After 60 timesteps
│   ├── solution_0120.h5
│   ├── ...
│   ├── solution.xdmf         # ParaView wrapper
│   ├── fracture_geometry.csv # Fracture length/height vs time
│   ├── stress_intensity.csv  # K_I evolution
│   └── well_data.csv         # BHP, rates
```

## LEFM Example Deep Dive

The improved `lefm_fracture_growth.config` demonstrates fracture mechanics in detail:

**Phase 1: Pressure Build-up (0-60s)**
- Injection at 50 L/s builds wellbore pressure
- Stress intensity factor K_I increases
- Monitoring for initiation criterion

**Phase 2: Initiation**
- Fracture initiates when K_I = K_Ic
- K_Ic = 1.0 MPa·m^0.5 (typical for shale)
- Small initial notch provides stress concentration

**Phase 3: Propagation (60-3600s)**
- Velocity: v = C × (K_I - K_Ic) / K_Ic
- Width from elasticity: w(x) = (4(1-ν²)/E) × √(a²-x²) × (P - σ_min)
- Leak-off reduces effective injection

**Validation:**
- Compare with PKN model: L(t) ∝ t^(4/5)
- Compare with KGD model: L(t) ∝ t^(2/3)
- Check scaling regimes (toughness vs viscosity dominated)

**Expected Results:**
- Final fracture length: ~100-150 m
- Final fracture height: ~30-50 m
- Max fracture width: ~5-10 mm (at wellbore)
- Net pressure: ~2-4 MPa above σ_min

## References

**LEFM Theory:**
- Sneddon (1946): Crack opening displacement
- Perkins & Kern (1961): PKN model
- Geertsma & de Klerk (1969): KGD model
- Adachi & Detournay (2002): Scaling regimes

**Numerical Methods:**
- PETSc documentation: https://petsc.org/
- HDF5 documentation: https://www.hdfgroup.org/

**Code Documentation:**
- User Guide: `docs/USER_GUIDE.md`
- API Reference: `docs/API_REFERENCE.md`
- Configuration: `docs/CONFIGURATION.md`

## Contributing

To add a new example:

1. Create config file: `config/my_example.config`
2. Document physics and parameters inline
3. Set output_format = HDF5 (default)
4. Test with: `./simulator -c config/my_example.config`
5. Update this README with example description

## Support

- Documentation: `docs/`
- Issues: GitHub issue tracker
- Discussions: GitHub discussions

---

**Note**: All examples now default to HDF5 output for efficiency. VTK remains available as a secondary option for legacy workflows. The config-driven approach eliminates the need for separate example executables and makes parameter exploration much easier.
