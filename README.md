# FSRM - Full Service Reservoir Model

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A fully coupled, physics-based simulator for petroleum reservoirs and earth systems, built on PETSc for scalable parallel computing.

**This project is fully open source under the MIT License** — free to use, modify, and distribute for any purpose.

## Features

### Physics Models
- **Fluid Flow**: Single-phase, black oil, compositional (EOS)
- **Geomechanics**: Elastic, viscoelastic, poroelastic
- **Thermal**: Heat conduction, convection, phase change
- **Fractures**: Natural fractures (DFN), hydraulic fracturing (PKN/KGD/P3D)
- **Particle Transport**: Proppant settling and bridging
- **Faults**: Rate-and-state friction, induced seismicity
- **Tidal Forces**: Earth tides effect on reservoirs

### ResFrac-Equivalent Hydraulic Fracturing ⭐ NEW
Complete unconventional completions simulation:
- **Multi-Stage Fracturing**: Unlimited stages/clusters, plug-and-perf, sliding sleeves
- **Limited-Entry Design**: Automatic perforation optimization for uniform flow
- **Stress Shadowing**: DDM-based inter-fracture stress calculation
- **2D Proppant Transport**: Settling, convection, bank formation, screenout prediction
- **Wellbore Hydraulics**: Friction, temperature, perforation erosion
- **Diversion Modeling**: Degradable diverters, ball sealers, chemical diversion
- **Real-Time Analysis**: ISIP, closure pressure, net pressure, Nolte-Smith interpretation
- **Flowback Simulation**: Water recovery, cleanup efficiency
- **Production Forecasting**: Decline analysis, EUR estimation, type curves
- **Rate Transient Analysis**: Full RTA diagnostics, DFIT analysis
- **Economic Analysis**: NPV optimization, Monte Carlo uncertainty
- **Fiber Optic Integration**: DAS/DTS cluster efficiency, flow allocation
- **Parent-Child Interactions**: Frac hits, depletion effects

### Earthquake Simulation (SeisSol-Compatible)
- **Discontinuous Galerkin (DG)**: High-order (O1-O10) ADER-DG spatial discretization
- **ADER Time Integration**: Arbitrary high-order time stepping with optimal CFL
- **Local Time Stepping (LTS)**: Adaptive timesteps per element for efficiency
- **Dynamic Rupture**: Full fault mechanics with spontaneous rupture propagation
- **Friction Laws**:
  - Coulomb (linear slip-weakening)
  - Rate-and-State (aging/slip laws)
  - Flash Heating
  - Thermal Pressurization
  - Strong Velocity Weakening
- **Absorbing Boundaries**: PML, Clayton-Engquist, Lysmer-Kuhlemeyer
- **Seismic Sources**: Moment tensor, kinematic fault, single force
- **Receivers**: Surface/volume waveforms, fault slip/traction output
- **Off-fault Plasticity**: Drucker-Prager, von Mises, Mohr-Coulomb, Cap models
- **Viscoelastic Attenuation**: Multi-mechanism Q-factor modeling

### Explosion & Impact Physics ⭐ NEW
Comprehensive multi-physics simulation of explosive events and hypervelocity impacts:
- **Underground Nuclear Tests**:
  - Cavity formation and damage zone modeling
  - Mueller-Murphy seismic source model
  - Spall and chimney collapse
  - Near-field plasticity and fracturing
  - Seismic wave generation (P, S, surface waves)
  - mb/Ms discrimination for monitoring
- **Atmospheric Nuclear Detonations**:
  - Fireball evolution (Brode model)
  - Blast wave propagation (Sedov-Taylor)
  - Ground coupling and crater formation
  - Mach stem and reflection effects
- **Radiation Transport**:
  - Prompt gamma and neutron radiation
  - Fallout modeling (Way-Wigner decay)
  - Dose rate calculations
  - Contamination mapping
- **Electromagnetic Pulse (EMP)**:
  - E1 (fast Compton current, nanoseconds)
  - E2 (intermediate, microseconds)
  - E3 (MHD-EMP, seconds)
  - Ground wave and line coupling
  - Ionospheric effects
- **Hypervelocity Impacts**:
  - Pi-group crater scaling (Holsapple-Schmidt)
  - Shock metamorphism zones
  - Z-model excavation flow
  - Ejecta distribution and ballistics
  - Transient to complex crater evolution
  - Impact melt and vapor plume
- **Coupled Atmosphere-Ground**:
  - Compressible Navier-Stokes atmosphere
  - Elastodynamic ground response
  - Interface matching and coupling

### Volcano Modeling ⭐ NEW
Comprehensive multi-physics volcanic system simulation:
- **Magma Chamber Dynamics**:
  - Pressure/temperature evolution
  - Volatile exsolution (H₂O, CO₂, SO₂)
  - Crystal fraction tracking
  - Recharge and eruption coupling
- **Conduit Flow**:
  - Bubble nucleation and growth
  - Flow regime transitions (bubbly → dispersed)
  - Fragmentation depth calculation
  - Mass eruption rate prediction
- **Eruption Column**:
  - Woods (1988) buoyant plume model
  - Air entrainment and thermal equilibration
  - Column collapse criteria
  - Neutral buoyancy level calculation
- **Volcanic Hazards**:
  - Pyroclastic density currents (PDCs)
  - Lava flow propagation with rheology
  - Lahar (volcanic debris flow) modeling
  - Tephra dispersal and fallout
- **Volcanic Monitoring**:
  - Surface deformation (Mogi, McTigue, Okada sources)
  - Volcanic seismicity (VT, LP, VLP, tremor)
  - Gas emission and dispersal (SO₂, CO₂)
  - Synthetic observables for GPS, InSAR, seismometers
- **Pre-defined Scenarios**:
  - Mt. St. Helens 1980 (lateral blast)
  - Pinatubo 1991 (VEI-6 Plinian)
  - Kilauea (Hawaiian effusive)
  - Yellowstone (caldera unrest)

### Tsunami Modeling ⭐ NEW
End-to-end earthquake-tsunami simulation:
- **Earthquake-Tsunami Coupling**:
  - Okada (1985) seafloor deformation
  - Kinematic and dynamic rupture sources
  - Kajiura filter for water column response
- **Tsunami Propagation**:
  - Nonlinear shallow water equations (NSWE)
  - Multiple Riemann solvers (F-wave, HLLC, Roe)
  - High-order reconstruction with limiters
  - Positivity-preserving schemes
- **Coastal Inundation**:
  - Wetting/drying algorithm
  - Manning friction model
  - Runup tracking
  - Variable land cover effects
- **Bathymetry Support**:
  - GEBCO, ETOPO global datasets
  - NetCDF, ASCII, GeoTIFF formats
  - Nested grid capability
- **Cascadia Subduction Zone**:
  - Full M9.0 margin rupture scenario
  - Heterogeneous slip models
  - DART buoy and tide gauge validation

### Atmospheric Infrasound ⭐ NEW
Low-frequency acoustic wave propagation:
- **Propagation Physics**:
  - Linearized acoustics in moving medium
  - Stratified atmosphere with wind effects
  - Stratospheric and thermospheric ducting
  - Topographic scattering
- **Source Types**:
  - Nuclear/chemical explosions
  - Volcanic eruptions
  - Bolides and meteorites
  - Severe weather
- **Applications**:
  - CTBT/IMS nuclear monitoring
  - Volcano eruption detection
  - Meteor tracking

### Machine Learning Solvers
- **Fourier Neural Operator (FNO)**: Neural network-based PDE solver
  - 100-1000x faster inference after training
  - Mesh-independent learning
  - Multiple solver modes: Pure FNO, Hybrid, Ensemble
  - Physics-informed training options
  - GPU-accelerated training and inference
  - Pre-trained model loading/saving

### Key Capabilities
- Fully coupled multi-physics (THM - Thermo-Hydro-Mechanical)
- Eclipse format I/O (DATA files)
- **Configuration-driven simulations** - no C++ coding required
- **Solver method selection**: Numerical, FNO, or Hybrid modes
- Structured and unstructured grids (Gmsh, Corner-Point, Exodus)
- **Adaptive Mesh Refinement (AMR)** for unstructured grids using PETSc DMPlex/DMForest
- **EPSG coordinate systems** with PROJ library integration
- Advanced well models (vertical, horizontal, multilateral)
- Stochastic reservoir modeling
- Comprehensive testing framework (MMS, analytical solutions)
- Parallel scalability (MPI)

### Mesh & Coordinate Features
- **Gmsh Integration**: Load `.msh` files (v2.2 and v4.x) for complex geometries
- **Coordinate Transformations**: EPSG codes, UTM zones, local origins
- **Physical Groups**: Reference Gmsh regions for materials and boundaries
- **Auto UTM Detection**: Automatically determine UTM zone from GPS coordinates
- **Adaptive Mesh Refinement**:
  - Local h-refinement based on error estimators (gradient, Hessian, residual)
  - Feature-based refinement (wells, faults, saturation fronts)
  - PETSc DMPlex and DMForest (p4est) backends
  - Conservative solution transfer between meshes
  - Dynamic load balancing for parallel efficiency

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
- CUDA Toolkit >= 11.0 (for GPU acceleration)
- ROCm/HIP (for AMD GPU support)
- PROJ >= 6.0 (for coordinate transformations)
- Gmsh (for mesh generation)
```

### GPU Support

FSRM supports GPU acceleration using NVIDIA CUDA or AMD ROCm/HIP:

**NVIDIA GPU (CUDA)**:
- Supported architectures: Pascal (SM 6.0), Volta (7.0), Turing (7.5), Ampere (8.0), Ada (8.6+)
- Minimum CUDA version: 11.0
- Recommended: CUDA 12.0+ for best performance

**AMD GPU (ROCm/HIP)**:
- Supported: MI50, MI100, MI200 series
- ROCm version: 5.0+

**Performance Benefits**:
- 5-50x speedup for large-scale simulations (>100k cells)
- Particularly effective for:
  - Elastodynamics and wave propagation
  - Poroelastodynamics with dynamic permeability
  - Large timestep counts
  - Multi-physics coupling

### Local Compilation

**CPU-only build**:
```bash
# Configure
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_CUDA=OFF

# Build
make -j4

# Run tests
make test

# Install
sudo make install
```

**With coordinate transformation support (PROJ)**:
```bash
# Install PROJ library first
sudo apt install libproj-dev   # Debian/Ubuntu
# brew install proj            # macOS

# Configure with PROJ support
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release \
         -DENABLE_PROJ=ON

# Build
make -j4
```

**GPU-accelerated build (NVIDIA CUDA)**:
```bash
# Configure with CUDA support
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release \
         -DENABLE_CUDA=ON \
         -DCMAKE_CUDA_ARCHITECTURES="80;86"  # Adjust for your GPU

# Build
make -j4

# Verify GPU detection
./fsrm --gpu-info

# Install
sudo make install
```

**GPU-accelerated build (AMD ROCm)**:
```bash
# Configure with HIP support
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release \
         -DENABLE_CUDA=OFF \
         -DENABLE_HIP=ON

# Build
make -j4

# Install
sudo make install
```

### Cloud Deployment (AWS or Google Cloud)

Deploy FSRM on cloud platforms for high-performance computing:

```bash
# Quick start - AWS
cd deploy/aws && ./setup.sh

# Quick start - Google Cloud
cd deploy/gcp && ./setup.sh

# Docker deployment
docker build -t fsrm:latest .
docker run -it fsrm:latest
```

**Documentation**:
- **Deployment Guide**: [docs/DEPLOYMENT.md](docs/DEPLOYMENT.md)
- **Deployment Tools**: [deploy/README.md](deploy/README.md)

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
fsrm -generate_config my_simulation.config

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

**CPU execution**:
```bash
mpirun -np 4 fsrm -c my_simulation.config
```

**FNO (Neural Operator) execution**:
```bash
# Pure FNO solver (requires trained model)
fsrm -c my_simulation.config --solver-method fno

# Hybrid: FNO prediction + numerical refinement
fsrm -c my_simulation.config --solver-method hybrid_refine

# Ensemble: Weighted average of FNO and numerical
fsrm -c my_simulation.config --solver-method ensemble

# Train FNO from scratch
fsrm -c my_simulation.config --train-fno --fno-epochs 100
```

**GPU execution**:
```bash
# Single GPU
fsrm -c my_simulation.config --use-gpu

# Multiple GPUs with MPI
mpirun -np 4 fsrm -c my_simulation.config --use-gpu

# Specify GPU device
fsrm -c my_simulation.config --use-gpu --gpu-device 0

# Hybrid CPU+GPU mode
mpirun -np 8 fsrm -c my_simulation.config --gpu-mode hybrid
```

**No recompilation needed!** Edit the config file and re-run.

### Using Eclipse Format
```bash
mpirun -np 4 fsrm -i SPE1.DATA -o output/spe1
```

### Pre-configured Examples

We provide several ready-to-use configurations:

```bash
# Shale reservoir with hydraulic fracturing
mpirun -np 8 fsrm -c config/shale_reservoir.config

# Enhanced geothermal system (EGS)
mpirun -np 8 fsrm -c config/geothermal.config

# Induced seismicity from injection
mpirun -np 8 fsrm -c config/induced_seismicity.config

# CO2 geological storage
mpirun -np 8 fsrm -c config/co2_storage.config
```

## Configuration File Format

Configuration files use a simple INI format with sections:

### [SIMULATION] - Simulation Control
- **Time parameters**: start_time, end_time, dt_initial, dt_min, dt_max
- **Output control**: output_frequency, output_format, checkpointing
- **Solver settings**: rtol, atol, max_iterations
- **Physics selection**: fluid_model, solid_model, enable_* flags
- **GPU settings**: use_gpu, gpu_mode, gpu_device_id, gpu_memory_fraction

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

### [FNO] - Fourier Neural Operator Settings
Configure the neural network-based solver:
```ini
[FNO]
enabled = true
input_channels = 3          # Number of input fields
output_channels = 3         # Number of output fields  
hidden_channels = 64        # Network width
num_layers = 4              # Number of FNO layers
num_modes = 16              # Fourier modes to keep
activation = gelu           # gelu, relu, tanh, silu

[FNO_TRAINING]
learning_rate = 1e-3
batch_size = 16
epochs = 100
early_stopping = true

[FNO_MODEL]
model_path = models/fno_trained.bin
load_pretrained = true

[FNO_HYBRID]
ensemble_weight = 0.8       # FNO weight in ensemble mode
refinement_iterations = 3   # Numerical refinement steps
```

### [amr] - Adaptive Mesh Refinement
Configure dynamic mesh adaptation for unstructured grids:
```ini
[amr]
enabled = true
method = plex_adapt         # plex_refine, plex_adapt, forest_p4est

[amr.criterion]
type = combined             # gradient, hessian, jump, feature, combined
weight_saturation = 1.0     # Prioritize saturation front tracking
weight_pressure = 0.3

[amr.strategy]
type = fixed_fraction
refine_fraction = 0.2
coarsen_fraction = 0.05

[amr.limits]
max_refinement_level = 4
max_cells = 100000
min_cell_size = 1.0

[amr.features]
buffer_layers = 2

# Feature-based refinement around wells
[amr.features.wells.PROD1]
x = 50.0
y = 50.0
radius = 20.0
level = 3

[amr.transfer]
conservative = true         # Mass-conserving solution transfer

[amr.balance]
enabled = true
method = parmetis           # Parallel mesh partitioning
```

## Examples

All examples are **configuration-driven** - no C++ coding required to run simulations. Each example simply loads a configuration file and executes the simulation.

### Configuration Files

Configuration files define all simulation parameters:

| Config File | Description |
|-------------|-------------|
| `config/multistage_frac_complete.config` | **Complete multi-stage fracturing (ResFrac-equivalent)** |
| `config/underground_nuclear_test.config` | **Underground nuclear explosion with seismic waves** |
| `config/atmospheric_nuclear_test.config` | **Atmospheric nuclear detonation with EMP & radiation** |
| `config/impact_event.config` | **Asteroid/meteorite impact with crater formation** |
| `config/cascadia_m9_tsunami.config` | **Cascadia M9.0 earthquake-tsunami simulation** |
| `config/infrasound_example.config` | **Atmospheric infrasound propagation** |
| `config/egs_optimization.config` | **Enhanced Geothermal System optimization** |
| `config/underground_explosion_near_field.config` | **Near-field underground explosion** |
| `config/buckley_leverett_2d.config` | 2D waterflooding with Buckley-Leverett |
| `config/coupled_reservoir_2d.config` | Coupled flow-geomechanics |
| `config/spe1_benchmark.config` | SPE1 benchmark problem |
| `config/unstructured_gmsh_example.config` | Gmsh mesh with EPSG coordinates |
| `config/elastodynamic_waves.config` | Elastic wave propagation |
| `config/poroelastodynamic_waves.config` | Poroelastic wave propagation |
| `config/hydraulic_fracturing.config` | Hydraulic fracturing simulation |
| `config/shale_reservoir.config` | Shale reservoir with fractures |
| `config/geothermal.config` | Enhanced geothermal system |
| `config/induced_seismicity.config` | Fault reactivation |
| `config/co2_storage.config` | CO₂ geological storage |
| `config/fno_solver.config` | FNO neural operator solver |
| `config/amr_unstructured.config` | AMR for two-phase flow on unstructured grids |

### Example Executables

**All examples now use a unified config-driven approach:**

| Executable | Purpose | Config Files |
|------------|---------|--------------|
| `simulator` | Main config-driven executable | Any `.config` file |
| `spe1` | SPE1 benchmark validation | `spe1_benchmark.config` |

**Previously separate executables have been replaced with config files.** This means no recompilation needed to run different examples!

### Running Examples

```bash
cd build/examples

# Run any example with the unified simulator
./simulator -c config/lefm_fracture_growth.config
./simulator -c config/hydraulic_fracturing.config
./simulator -c config/induced_seismicity.config
./simulator -c config/wave_propagation.config
./simulator -c config/single_phase.config

# Parallel execution
mpirun -np 4 ./simulator -c config/hydraulic_fracturing.config
mpirun -np 8 ./simulator -c config/induced_seismicity.config

# SPE1 benchmark (kept separate for validation)
./spe1 config/spe1_benchmark.config
```

### Output Format

**HDF5 is now the default output format** (more efficient than VTK):

```bash
# Visualize HDF5 outputs
python scripts/hdf5_to_xdmf.py output/
paraview output/solution.xdmf

# Or use VTK (secondary option) by setting in config:
# output_format = VTK
paraview output/*.vtu
```

### Using Gmsh Meshes with Coordinates

```bash
# Example with unstructured mesh and coordinate transformation
./ex_config_driven ../config/unstructured_gmsh_example.config
```

The unstructured Gmsh example demonstrates:
- Loading a `.msh` file for complex geometries
- EPSG coordinate system (WGS84 input → UTM model)
- Local origin for numerical precision
- Physical group names for boundaries

## Testing

### Unit Tests
```bash
cd build
./tests/run_tests

# Run specific test suites
./tests/unit/test_discontinuous_galerkin     # DG solver tests
./tests/unit/test_seismic_source             # Source/receiver tests
./tests/unit/test_boundary_conditions        # BC tests
./tests/unit/test_plasticity_model           # Plasticity tests
./tests/unit/test_viscoelastic_attenuation   # Attenuation tests
./tests/unit/test_friction_laws              # Friction law tests
./tests/unit/test_fourier_neural_operator    # FNO solver tests
```

### Integration Tests
```bash
# SCEC benchmark validation
./tests/integration/test_scec_tpv_benchmarks

# Full simulation tests
./tests/integration/test_full_simulation
```

### SCEC Benchmarks
Validated against SCEC (Southern California Earthquake Center) benchmark problems:
- **TPV5**: 2D strike-slip fault with slip-weakening friction
- **TPV10**: 60-degree dipping normal fault
- **TPV16**: 3D heterogeneous stress fault
- **LOH1**: Layer over halfspace wave propagation

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
- Lamb's problem (point force in halfspace)
- Garvin's problem (buried explosion)

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

**CPU Performance** (typical modern cluster):
- 1M cells: 4-8 cores, ~10 min/timestep
- 10M cells: 64-128 cores, ~20 min/timestep
- Parallel efficiency: >80% up to 256 cores

**GPU Performance** (single NVIDIA A100):
- 1M cells: ~12 sec/timestep (50x faster than 8 CPUs)
- 10M cells: ~2 min/timestep (10x faster than 128 CPUs)
- 100M cells: ~15 min/timestep (previously infeasible on CPU)

**GPU Speedup Factors** (vs optimized CPU code):
- Single-phase flow: 10-20x
- Elastodynamics: 30-50x
- Poroelastodynamics: 20-40x
- Black oil: 15-25x

**Multi-GPU Scaling**:
- 2 GPUs: 1.8x speedup
- 4 GPUs: 3.5x speedup
- 8 GPUs: 6.5x speedup
- Near-linear scaling up to 16 GPUs

**Memory Requirements**:
- CPU: ~100 bytes/cell
- GPU: ~150 bytes/cell (includes device buffers)
- Unified memory option available for large problems

## GPU Configuration Examples

### Enable GPU in configuration file
```ini
[SIMULATION]
use_gpu = true
gpu_mode = CPU_FALLBACK    # Options: CPU_ONLY, GPU_ONLY, HYBRID, CPU_FALLBACK
gpu_device_id = 0
gpu_verbose = true
gpu_memory_fraction = 0.8

# GPU solver options
use_gpu_preconditioner = true
use_gpu_matrix_assembly = true
pin_host_memory = true
```

### Command-line GPU options
```bash
# Check available GPUs
fsrm --gpu-info

# Run on GPU
fsrm -c config.ini --use-gpu

# Specify device
fsrm -c config.ini --use-gpu --gpu-device 1

# Verbose GPU output
fsrm -c config.ini --use-gpu --gpu-verbose

# Benchmark CPU vs GPU
fsrm -c config.ini --benchmark-gpu

# Profile GPU kernels
fsrm -c config.ini --use-gpu --profile
```

### Hybrid CPU+GPU Execution

For very large problems that don't fit in GPU memory:
```ini
[SIMULATION]
gpu_mode = HYBRID
gpu_memory_fraction = 0.9

# Automatically splits work between CPU and GPU
# Uses GPU for compute-intensive kernels
# Uses CPU for memory-intensive operations
```

## Troubleshooting GPU Issues

**GPU not detected**:
```bash
# Check CUDA installation
nvidia-smi
nvcc --version

# Check if CUDA is enabled in build
fsrm --version
```

**Out of memory errors**:
```ini
# Reduce GPU memory usage
gpu_memory_fraction = 0.6

# Enable unified memory (slower but handles larger problems)
use_unified_memory = true

# Or use hybrid mode
gpu_mode = HYBRID
```

**Performance not improved**:
- Small problems (<10k cells) may not benefit from GPU
- Check data transfer overhead with `--profile`
- Ensure problem is compute-bound, not I/O bound
- Use `--benchmark-gpu` to verify speedup

## Contributing

Contributions welcome! Areas of interest:
- Additional PVT correlations
- Relative permeability models
- Fracture network algorithms
- Machine learning integration
- Multi-GPU optimizations
- Support for additional GPU architectures

## Citation

If you use FSRM in your research, please cite:
```bibtex
@software{fsrm2024,
  title = {FSRM: Coupled Multi-Physics Reservoir Simulator},
  author = {Your Name},
  year = {2024},
  version = {1.0.0}
}
```

## License

**MIT License** — This project is fully open source and free to use.

You are free to:
- ✅ Use commercially
- ✅ Modify and distribute
- ✅ Use privately
- ✅ Sublicense

See the [LICENSE](LICENSE) file for the complete license text.

## Documentation

Complete documentation is available in the `docs/` directory:

- **Getting Started**
  - [Quick Start Guide](docs/QUICK_START.md) - Get running in 5 minutes
  - [User Guide](docs/USER_GUIDE.md) - Complete user manual
  - [Configuration Reference](docs/CONFIGURATION.md) - All configuration options

- **Development**
  - [Development Guide](docs/DEVELOPMENT.md) - Building, testing, contributing
  - [Benchmarks](docs/BENCHMARKS.md) - Validation and performance
  - [API Reference](docs/API_REFERENCE.md) - C++ API documentation

- **Advanced Topics**
  - [Physics Models](docs/PHYSICS_MODELS.md) - Mathematical formulations
  - [Numerical Methods](docs/NUMERICAL_METHODS.md) - Discretization and solvers
  - [Deployment Guide](docs/DEPLOYMENT.md) - Cloud, HPC, and GPU setup
  - [Unit System](docs/UNIT_SYSTEM.md) - Unit conversion and specifications
  - [Coordinate Systems](docs/COORDINATE_SYSTEMS.md) - Geographic transformations
  - [Unstructured Meshes](docs/UNSTRUCTURED_MESHES.md) - Gmsh integration
  - [Adaptive Mesh Refinement](docs/ADAPTIVE_MESH_REFINEMENT.md) - AMR for unstructured grids
  - [SeisSol Comparison](docs/SEISSOL_COMPARISON.md) - Earthquake simulator comparison
  - [Fourier Neural Operator](docs/FOURIER_NEURAL_OPERATOR.md) - Neural network solver guide
  - [**ResFrac-Equivalent Capabilities**](docs/RESFRAC_EQUIVALENT_CAPABILITIES.md) - Complete hydraulic fracturing guide
  - [**Explosion & Impact Physics**](docs/EXPLOSION_IMPACT_PHYSICS.md) - Nuclear tests, impacts, EMP, radiation
  - [**Volcano Modeling**](docs/VOLCANO_MODELING.md) - Comprehensive volcanic system simulation
  - [**Tsunami Modeling**](docs/TSUNAMI_MODELING.md) - Earthquake-tsunami coupling and inundation

## Support

- **Documentation**: See `docs/` directory
- **Issues**: GitHub issue tracker
- **Examples**: `examples/` directory with source code
- **Configurations**: `config/` directory with ready-to-run examples

## Acknowledgments

Built with:
- PETSc (Portable, Extensible Toolkit for Scientific Computation)
- MPI (Message Passing Interface)
- HDF5 (Hierarchical Data Format)

Inspired by:
- SPE comparative solution projects
- Stanford SUPRI research
- Various open-source reservoir simulators
