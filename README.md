# ReservoirSim - Comprehensive Petroleum Reservoir Simulator

A fully coupled, physics-based simulator for petroleum reservoirs and earth systems, built on PETSc for scalable parallel computing.

## Features

### Physics Models
- **Fluid Flow**: Single-phase, black oil, compositional (EOS)
- **Geomechanics**: Elastic, viscoelastic, poroelastic
- **Thermal**: Heat conduction, convection, phase change
- **Fractures**: Natural fractures (DFN), hydraulic fracturing (PKN/KGD/P3D)
- **Particle Transport**: Proppant settling and bridging
- **Faults**: Rate-and-state friction, induced seismicity
- **Tidal Forces**: Earth tides effect on reservoirs

### Key Capabilities
- Fully coupled multi-physics (THM - Thermo-Hydro-Mechanical)
- Eclipse format I/O (DATA files)
- User-editable configuration files
- Structured and unstructured grids
- Advanced well models (vertical, horizontal, multilateral)
- Stochastic reservoir modeling
- Comprehensive testing framework (MMS, analytical solutions)
- Parallel scalability (MPI)

## Building

### Prerequisites
```bash
# Required
- CMake >= 3.15
- C++17 compiler
- PETSc >= 3.15
- MPI
- HDF5

# Optional
- GTest (for testing)
- gnuplot (for visualization)
- ParaView (for 3D visualization)
```

### Compilation
```bash
# Configure
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release

# Build
make -j4

# Run tests
make test

# Install
sudo make install
```

## Usage

### Using Configuration Files (Recommended)

Configuration files allow you to specify **ALL** simulation parameters without recompiling:
- Rock/formation properties
- Fluid properties  
- Well locations, types, and controls
- Fracture properties and locations
- Fault properties and friction laws
- Particle/proppant properties
- Boundary conditions
- Initial conditions
- All physics settings

#### 1. Generate a template configuration
```bash
# Generate a basic template
reservoirsim -generate_config my_simulation.config

# Or use the complete template with all options
cp config/complete_template.config my_simulation.config
```

#### 2. Edit the configuration file
```ini
[SIMULATION]
start_time = 0.0
end_time = 86400.0        # 1 day
fluid_model = BLACK_OIL
enable_geomechanics = true

[ROCK]
porosity = 0.20           # 20%
permeability_x = 100.0    # milliDarcy
youngs_modulus = 10.0e9   # 10 GPa

[FLUID]
oil_viscosity = 0.005     # Pa·s (5 cP)
solution_GOR = 100.0      # scf/stb

[WELL1]
name = PROD1
type = PRODUCER
i = 10
j = 10
k = 5
control_mode = RATE
target_value = 0.01       # m³/s

[FRACTURE1]
type = HYDRAULIC
location = 500, 500, 50, 0, 1, 0
aperture = 0.001          # m
enable_propagation = true

[BC1]
type = DIRICHLET
field = PRESSURE
location = XMIN
value = 20.0e6            # Pa

[IC1]
field = PRESSURE
distribution = GRADIENT
value = 20.0e6
gradient = 0, 0, 10000    # Hydrostatic
```

#### 3. Run simulation
```bash
mpirun -np 4 reservoirsim -c my_simulation.config
```

**No recompilation needed!** Edit the config file and re-run.

### Using Eclipse Format
```bash
mpirun -np 4 reservoirsim -i SPE1.DATA -o output/spe1
```

### Pre-configured Examples

We provide several ready-to-use configurations:

```bash
# Shale reservoir with hydraulic fracturing
mpirun -np 8 reservoirsim -c config/shale_reservoir.config

# Enhanced geothermal system (EGS)
mpirun -np 8 reservoirsim -c config/geothermal.config

# Induced seismicity from injection
mpirun -np 8 reservoirsim -c config/induced_seismicity.config

# CO2 geological storage
mpirun -np 8 reservoirsim -c config/co2_storage.config
```

## Configuration File Format

Configuration files use a simple INI format with sections:

### [SIMULATION] - Simulation Control
- **Time parameters**: start_time, end_time, dt_initial, dt_min, dt_max
- **Output control**: output_frequency, output_format, checkpointing
- **Solver settings**: rtol, atol, max_iterations
- **Physics selection**: fluid_model, solid_model, enable_* flags

### [GRID] - Mesh Definition
- **Dimensions**: nx, ny, nz (number of cells)
- **Domain size**: Lx, Ly, Lz (meters)
- **Mesh type**: use_unstructured, mesh_file

### [ROCK] - Material Properties
- **Porosity and permeability** (milliDarcy)
- **Mechanical properties**: Young's modulus, Poisson's ratio, density
- **Viscoelastic properties**: relaxation_time, viscosity
- **Thermal properties**: thermal_conductivity, heat_capacity
- **Fracture properties**: fracture_toughness, fracture_energy

### [FLUID] - Fluid Properties
- **Single-phase**: density, viscosity, compressibility
- **Black oil**: oil/gas/water properties, solution GOR
- **Compositional**: component arrays (MW, Tc, Pc, omega)

### [WELL1], [WELL2], ... - Well Definitions
- **Location**: i, j, k (grid indices)
- **Type**: PRODUCER, INJECTOR, OBSERVATION
- **Control**: control_mode (RATE/BHP/THP), target_value
- **Constraints**: max_rate, min_bhp
- **Geometry**: diameter, skin

### [FRACTURE1], [FRACTURE2], ... - Fracture Definitions
- **Type**: NATURAL, HYDRAULIC, INDUCED_THERMAL
- **Location**: x, y, z, nx, ny, nz (center + normal)
- **Properties**: aperture, permeability, toughness, energy
- **Options**: enable_propagation, enable_proppant

### [FAULT1], [FAULT2], ... - Fault Definitions
- **Geometry**: strike, dip, length, width
- **Friction**: static_friction, dynamic_friction, cohesion
- **Rate-state**: use_rate_state, a_parameter, b_parameter, Dc_parameter

### [PARTICLE] - Proppant/Particle Properties
- **Size**: diameter, density
- **Transport**: concentration, diffusivity
- **Physics**: enable_settling, enable_bridging

### [BC1], [BC2], ... - Boundary Conditions
- **Type**: DIRICHLET, NEUMANN, ROBIN
- **Field**: PRESSURE, TEMPERATURE, DISPLACEMENT
- **Location**: XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX
- **Value**: value, gradient

### [IC1], [IC2], ... - Initial Conditions
- **Field**: PRESSURE, TEMPERATURE, SATURATION_OIL, etc.
- **Distribution**: UNIFORM, GRADIENT, FILE
- **Value**: value, gradient, file

### Multiple Rock Types
Define heterogeneous properties using multiple sections:
```ini
[ROCK1]  # Reservoir sandstone
porosity = 0.25
permeability_x = 500.0

[ROCK2]  # Shale barrier
porosity = 0.05
permeability_x = 0.001
```

## Examples

### Basic Examples
- `ex01_single_phase.cpp` - Single-phase Darcy flow
- `ex02_black_oil.cpp` - Three-phase black oil
- `ex03_compositional.cpp` - Multi-component compositional
- `ex04_geomechanics.cpp` - Coupled geomechanics

### SPE Benchmarks
- `spe1.cpp` - SPE1 comparative solution project
- `spe3.cpp` - SPE3 gas cycling
- `spe9.cpp` - SPE9 waterflood
- `spe10.cpp` - SPE10 upscaling

### Advanced Examples
- `ex_hydraulic_fracturing.cpp` - Frac design with proppant
- `ex_induced_seismicity.cpp` - Fault reactivation
- `ex_thermal_recovery.cpp` - Steam injection
- `ex_stochastic_reservoir.cpp` - Monte Carlo uncertainty

### Running Examples
```bash
cd build/examples

# Basic single-phase flow
./ex01_single_phase

# SPE1 benchmark
./spe1

# Hydraulic fracturing
mpirun -np 4 ./ex_hydraulic_fracturing

# Induced seismicity
mpirun -np 8 ./ex_induced_seismicity
```

## Testing

### Unit Tests
```bash
cd build
./tests/run_tests
```

### Method of Manufactured Solutions
```bash
# Convergence tests
cd tests
./run_mms_tests.sh
```

### Analytical Benchmarks
- Theis solution (radial flow)
- Buckley-Leverett (two-phase)
- Mandel-Cryer (poroelasticity)
- Terzaghi (consolidation)

## Visualization

### VTK Output
```bash
# View with ParaView
paraview output/*.vtu
```

### Plots
Automatically generated plots:
- Pressure history
- Production curves
- Convergence behavior
- Scalability studies
- Error analysis

## Material Properties Database

Common reservoir types are pre-configured:

| Reservoir Type | Config File | Key Features |
|----------------|-------------|--------------|
| Conventional | default.config | Standard sandstone |
| Shale | shale_reservoir.config | Ultra-tight, viscoelastic |
| Geothermal | geothermal.config | High temperature, granite |
| Fault Zone | induced_seismicity.config | Rate-state friction |
| Storage | co2_storage.config | Deep saline, multi-component |

### Typical Property Ranges

**Permeability** (in milliDarcy):
- Conventional: 10 - 1000 mD
- Tight gas: 0.001 - 1 mD
- Shale: 0.0001 - 0.001 mD
- Fractured: 0.1 - 100 mD

**Porosity**:
- Sandstone: 0.15 - 0.30
- Carbonate: 0.05 - 0.25
- Shale: 0.02 - 0.08
- Granite: 0.001 - 0.01

**Young's Modulus**:
- Soft rock: 5 - 15 GPa
- Sandstone: 10 - 30 GPa
- Shale: 20 - 40 GPa
- Granite: 40 - 60 GPa

## PETSc Options

Fine-tune solvers via command line:
```bash
# Time stepping
-ts_type beuler          # Backward Euler
-ts_adapt_type basic     # Adaptive time stepping

# Nonlinear solver
-snes_type newtonls      # Newton line search
-snes_linesearch_type bt # Backtracking
-snes_monitor            # Monitor convergence

# Linear solver
-ksp_type gmres          # GMRES
-ksp_rtol 1e-6           # Relative tolerance
-pc_type hypre           # Hypre preconditioner
-pc_hypre_type boomeramg # Algebraic multigrid

# Monitoring
-log_view                # Performance summary
-snes_converged_reason   # Why converged/failed
```

## Performance

Typical performance on modern cluster:
- 1M cells: 4-8 cores, ~10 min/timestep
- 10M cells: 64-128 cores, ~20 min/timestep
- Parallel efficiency: >80% up to 256 cores

## Contributing

Contributions welcome! Areas of interest:
- Additional PVT correlations
- Relative permeability models
- Fracture network algorithms
- Machine learning integration
- GPU acceleration

## Citation

If you use ReservoirSim in your research, please cite:
```bibtex
@software{reservoirsim2024,
  title = {ReservoirSim: Coupled Multi-Physics Reservoir Simulator},
  author = {Your Name},
  year = {2024},
  version = {1.0.0}
}
```

## License

MIT License - see LICENSE file

## Support

- Documentation: See `docs/` directory
- Issues: GitHub issue tracker
- Examples: `examples/` directory
- Configs: `config/` directory

## Acknowledgments

Built with:
- PETSc (Portable, Extensible Toolkit for Scientific Computation)
- MPI (Message Passing Interface)
- HDF5 (Hierarchical Data Format)

Inspired by:
- SPE comparative solution projects
- Stanford SUPRI research
- Various open-source reservoir simulators
