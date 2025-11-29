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

### SCEC Benchmarks (Complete Suite)

#### `scec_tpv3.config`
- **Physics**: 2D in-plane strike-slip dynamic rupture
- **Features**: Linear slip-weakening, forced nucleation
- **Runtime**: ~10 minutes
- **Reference**: https://strike.scec.org/cvws/tpv3docs.html

#### `scec_tpv4.config`
- **Physics**: 3D strike-slip with rate-state friction
- **Features**: VW/VS zones, aging law, velocity perturbation
- **Runtime**: ~30 minutes
- **Reference**: https://strike.scec.org/cvws/tpv4docs.html

#### `scec_tpv5.config`
- **Physics**: 3D vertical strike-slip dynamic rupture
- **Features**: Heterogeneous stress, slip-weakening friction
- **Runtime**: ~20 minutes
- **Reference**: https://strike.scec.org/cvws/tpv5docs.html

#### `scec_tpv10.config`
- **Physics**: Branching fault rupture
- **Features**: Main fault + 30° branch, rate-state friction
- **Runtime**: ~45 minutes
- **Reference**: https://strike.scec.org/cvws/

#### `scec_tpv12.config` / `scec_tpv13.config`
- **Physics**: 60° dipping normal fault
- **Features**: Hanging wall/footwall asymmetry, free surface
- **TPV12**: Slip-weakening friction
- **TPV13**: Rate-state friction
- **Runtime**: ~30 minutes
- **Reference**: https://strike.scec.org/cvws/tpv12_13docs.html

#### `scec_tpv16.config`
- **Physics**: Rough fault surface rupture
- **Features**: Self-similar roughness, geometric complexity
- **Runtime**: ~45 minutes
- **Reference**: https://strike.scec.org/cvws/

#### `scec_tpv24.config` / `scec_tpv25.config`
- **Physics**: Japan benchmark with rate-state friction
- **Features**: Strong velocity weakening, flash heating
- **TPV24**: Homogeneous stress
- **TPV25**: Heterogeneous stress
- **Runtime**: ~1 hour
- **Reference**: https://strike.scec.org/cvws/tpv24_25docs.html

#### `scec_tpv101.config`
- **Physics**: Viscoelastic off-fault plasticity
- **Features**: Drucker-Prager plasticity, Q-factor attenuation
- **Runtime**: ~1 hour
- **Reference**: https://strike.scec.org/cvws/tpv101docs.html

#### `scec_loh1.config`
- **Physics**: Layer over halfspace wave propagation
- **Features**: 1 km layer, point explosion source
- **Runtime**: ~15 minutes
- **Reference**: https://strike.scec.org/scecpedia/LOH.1

#### `scec_loh2.config`
- **Physics**: Heterogeneous layer over halfspace
- **Features**: Double-couple source, surface waves, Love waves
- **Runtime**: ~20 minutes
- **Reference**: https://strike.scec.org/scecpedia/LOH.2

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
