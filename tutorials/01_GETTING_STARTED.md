# Tutorial 01: Getting Started with FSRM

Welcome to the Full-Scale Reservoir Modeling (FSRM) tutorial series! This first tutorial will guide you through installation, basic concepts, and running your first simulation.

## Table of Contents

1. [What is FSRM?](#what-is-fsrm)
2. [Prerequisites](#prerequisites)
3. [Installation](#installation)
4. [Project Structure](#project-structure)
5. [Your First Simulation](#your-first-simulation)
6. [Understanding the Output](#understanding-the-output)
7. [Common Workflows](#common-workflows)
8. [Next Steps](#next-steps)

---

## What is FSRM?

FSRM (Full-Scale Reservoir Modeling) is a high-performance, fully-coupled simulator for:

- **Reservoir Simulation**: Single-phase, black oil, compositional, CO2 storage
- **Geomechanics**: Linear elastic, viscoelastic, poroelastic, elastoplastic
- **Thermal Modeling**: Heat transfer with thermal-hydraulic-mechanical coupling
- **Fracture Mechanics**: Natural fractures (DFN) and hydraulic fracturing
- **Fault Mechanics**: Coulomb and rate-state friction with induced seismicity
- **Wave Propagation**: Elastodynamics and poroelastodynamics

### Key Features

- **Configuration-driven**: Change parameters without recompilation
- **GPU Accelerated**: CUDA/HIP support for large-scale problems
- **Parallel Scalable**: MPI-based for distributed computing
- **Multi-physics Coupling**: Seamless integration of flow, mechanics, thermal
- **Adaptive Mesh Refinement**: Dynamic mesh optimization
- **Industry Standard I/O**: Eclipse-compatible input/output

---

## Prerequisites

### Required Software

| Software | Version | Purpose |
|----------|---------|---------|
| CMake | ≥ 3.15 | Build system |
| C++17 Compiler | GCC 7+, Clang 6+ | Compilation |
| PETSc | ≥ 3.15 | Linear algebra, mesh management |
| MPI | OpenMPI or MPICH | Parallel computing |

### Optional Software

| Software | Version | Purpose |
|----------|---------|---------|
| HDF5 | ≥ 1.10 | Efficient I/O |
| CUDA | ≥ 11.0 | GPU acceleration |
| ParaView | ≥ 5.0 | Visualization |
| Python | ≥ 3.8 | Post-processing scripts |

---

## Installation

### Step 1: Install Dependencies (Ubuntu/Debian)

```bash
# Update package manager
sudo apt update

# Install build essentials
sudo apt install -y cmake g++ gfortran

# Install PETSc and dependencies
sudo apt install -y libpetsc-dev libhdf5-mpi-dev

# Set environment variables
export PETSC_DIR=/usr/lib/petsc
export PETSC_ARCH=""

# Add to ~/.bashrc for persistence
echo 'export PETSC_DIR=/usr/lib/petsc' >> ~/.bashrc
echo 'export PETSC_ARCH=""' >> ~/.bashrc
```

### Step 2: Build FSRM

```bash
# Clone repository (or navigate to your copy)
cd /path/to/FSRM

# Create build directory
mkdir build && cd build

# Configure with CMake
cmake .. -DCMAKE_BUILD_TYPE=Release

# Build (use all available cores)
make -j$(nproc)

# Verify installation
./fsrm -help
```

### Step 3: Verify Installation

```bash
# Run a quick test
./fsrm -c ../config/default.config

# Check for expected output
# Should see simulation progress and completion message
```

---

## Project Structure

Understanding the project layout helps you navigate and customize FSRM:

```
FSRM/
├── config/                 # Configuration files
│   ├── default.config         # Basic single-phase flow
│   ├── shale_reservoir.config # Shale with fractures
│   ├── geothermal.config      # Geothermal simulation
│   ├── induced_seismicity.config # Fault mechanics
│   ├── co2_storage.config     # CO2 sequestration
│   └── complete_template.config # ALL available options
│
├── include/                # Header files (API)
│   ├── Simulator.hpp          # Main simulator class
│   ├── PhysicsKernel.hpp      # Physics implementations
│   ├── WellModel.hpp          # Well modeling
│   ├── FractureModel.hpp      # Fracture mechanics
│   ├── FaultModel.hpp         # Fault mechanics
│   └── ...                    # Other components
│
├── src/                    # Source implementations
│   ├── main.cpp               # Entry point
│   ├── Simulator.cpp          # Simulator implementation
│   └── ...                    # Other sources
│
├── examples/               # Example programs
│   ├── ex_config_driven.cpp   # Config-driven example
│   ├── spe1.cpp              # SPE benchmark 1
│   └── scec_tpv*.cpp         # SCEC earthquake benchmarks
│
├── tests/                  # Unit and integration tests
│
├── docs/                   # Documentation
│   ├── USER_GUIDE.md
│   ├── CONFIGURATION.md
│   └── ...
│
└── tutorials/              # This tutorial series
```

---

## Your First Simulation

Let's run a simple single-phase flow simulation step by step.

### Step 1: Examine the Configuration

First, look at the default configuration:

```bash
# View the default config
cat config/default.config
```

Key sections in the config file:

```ini
# [SIMULATION] - Core parameters
[SIMULATION]
name = single_phase_example
start_time = 0.0
end_time = 86400.0          # 1 day in seconds
dt_initial = 3600.0         # 1 hour initial timestep
fluid_model = SINGLE_COMPONENT
solid_model = ELASTIC

# [GRID] - Domain definition
[GRID]
nx = 20
ny = 20
nz = 5
Lx = 1000.0                 # 1 km domain
Ly = 1000.0
Lz = 100.0

# [ROCK] - Material properties
[ROCK]
porosity = 0.20
permeability_x = 100.0      # milliDarcy
youngs_modulus = 10.0e9     # 10 GPa
poisson_ratio = 0.25

# [FLUID] - Fluid properties
[FLUID]
density = 1000.0            # kg/m³
viscosity = 0.001           # Pa·s (1 cP)
```

### Step 2: Run the Simulation

```bash
# Single processor run
./fsrm -c config/default.config

# Or parallel run with 4 processes
mpirun -np 4 ./fsrm -c config/default.config
```

### Step 3: Monitor Progress

During execution, you'll see output like:

```
================================================
  FSRM - Full-Scale Reservoir Modeling
================================================

Configuration: config/default.config
Processors: 4
Grid: 20 x 20 x 5 = 2000 cells

Initializing...
  Creating mesh...          [OK]
  Setting up fields...      [OK]
  Setting up physics...     [OK]
  Setting initial conditions... [OK]

Running simulation...
  Step     1 | t =    3600.0 s | dt =  3600.0 s | Iterations:  3
  Step     2 | t =    7200.0 s | dt =  3600.0 s | Iterations:  2
  Step     3 | t =   10800.0 s | dt =  3600.0 s | Iterations:  2
  ...
  Step    24 | t =   86400.0 s | dt =  3600.0 s | Iterations:  2

================================================
  Simulation Complete
================================================
  Total time:     86400.0 s
  Wall clock:        12.3 s
  Timesteps:           24
  Output files:    output/
```

---

## Understanding the Output

### Output Directory Structure

```
output/
├── pressure_0000.vtu      # Pressure field at timestep 0
├── pressure_0001.vtu      # Pressure field at timestep 1
├── ...
├── solution.pvd           # ParaView collection file
├── summary.csv            # Time series data
└── log.txt                # Simulation log
```

### Visualizing Results

Using ParaView:

```bash
# Open the collection file
paraview output/solution.pvd

# Or open individual timesteps
paraview output/pressure_*.vtu
```

In ParaView:
1. Click "Apply" to load data
2. Select "Pressure" in the dropdown
3. Click play to animate through timesteps

### Analyzing Summary Data

```python
# Python script to plot results
import pandas as pd
import matplotlib.pyplot as plt

# Load summary data
data = pd.read_csv('output/summary.csv')

# Plot pressure vs time
plt.figure(figsize=(10, 6))
plt.plot(data['time'] / 3600, data['pressure_avg'] / 1e6)
plt.xlabel('Time (hours)')
plt.ylabel('Average Pressure (MPa)')
plt.title('Reservoir Pressure Evolution')
plt.grid(True)
plt.savefig('pressure_plot.png')
plt.show()
```

---

## Common Workflows

### Workflow 1: Modify an Existing Config

```bash
# 1. Copy a template
cp config/default.config my_simulation.config

# 2. Edit parameters
nano my_simulation.config
# Change grid size, properties, etc.

# 3. Run
mpirun -np 4 ./fsrm -c my_simulation.config
```

### Workflow 2: Generate a Template

```bash
# Generate a complete template with all options
./fsrm -generate_config my_complete_config.config

# Edit the generated file
nano my_complete_config.config
```

### Workflow 3: Use Eclipse Input

If you have existing Eclipse input files:

```bash
# Run with Eclipse .DATA file
mpirun -np 4 ./fsrm -i SPE1.DATA -o output/spe1
```

### Workflow 4: Run Parameter Studies

```bash
# Script to run multiple simulations
#!/bin/bash
for perm in 10 50 100 200 500; do
    # Create config with different permeability
    sed "s/permeability_x = .*/permeability_x = $perm/" \
        config/default.config > temp_${perm}.config
    
    # Run simulation
    mpirun -np 4 ./fsrm -c temp_${perm}.config \
        -o output/perm_${perm}
done
```

---

## Command-Line Options

| Option | Description | Example |
|--------|-------------|---------|
| `-c <file>` | Configuration file | `-c config/default.config` |
| `-i <file>` | Eclipse input file | `-i SPE1.DATA` |
| `-o <prefix>` | Output prefix | `-o output/run1` |
| `-format <type>` | Output format | `-format VTK` |
| `-generate_config` | Generate template | `-generate_config template.config` |
| `-help` | Show help | `-help` |

### PETSc Options

Fine-tune solver behavior:

```bash
mpirun -np 4 ./fsrm -c config.config \
  -ts_type beuler \           # Backward Euler time stepping
  -snes_type newtonls \       # Newton line search
  -snes_rtol 1e-8 \           # Relative tolerance
  -ksp_type gmres \           # GMRES linear solver
  -pc_type ilu \              # ILU preconditioner
  -log_view                   # Performance summary
```

---

## Troubleshooting

### PETSc Not Found

```bash
# Check if PETSC_DIR is set
echo $PETSC_DIR

# If not, find PETSc installation
locate petscconf.h
# or
find /usr -name "petscconf.h" 2>/dev/null

# Set environment
export PETSC_DIR=/path/to/petsc
```

### MPI Errors

```bash
# Check MPI installation
which mpirun
mpirun --version

# Try single-processor first
./fsrm -c config/default.config

# Check for SSH issues (multi-node)
ssh localhost hostname
```

### Convergence Issues

If the simulation doesn't converge:

1. **Reduce timestep**:
   ```ini
   dt_initial = 100.0    # Smaller initial step
   dt_max = 3600.0       # Limit maximum step
   ```

2. **Increase iterations**:
   ```ini
   max_nonlinear_iterations = 100
   max_linear_iterations = 2000
   ```

3. **Relax tolerances** (temporarily):
   ```ini
   rtol = 1.0e-4
   atol = 1.0e-6
   ```

4. **Check boundary conditions** for physical consistency

### Memory Issues

```bash
# Check available memory
free -h

# Reduce grid size or use more processes
mpirun -np 8 ./fsrm -c config.config

# Enable GPU if available
[SIMULATION]
use_gpu = true
```

---

## Next Steps

Now that you've run your first simulation, explore these tutorials:

1. **[Tutorial 02: Configuration System](02_CONFIGURATION.md)** - Master the config file format
2. **[Tutorial 03: Physics Models](03_PHYSICS_MODELS.md)** - Understand available physics
3. **[Tutorial 04: Well Modeling](04_WELL_MODELING.md)** - Add production/injection wells
4. **[Tutorial 05: Fracture Modeling](05_FRACTURE_MODELING.md)** - Model fractures and hydraulic fracturing
5. **[Tutorial 06: Fault Mechanics](06_FAULT_MECHANICS.md)** - Induced seismicity analysis

### Practice Exercises

1. **Exercise 1**: Modify `default.config` to simulate a 2-day period with hourly output
2. **Exercise 2**: Double the grid resolution and compare runtime
3. **Exercise 3**: Change rock permeability and observe pressure evolution
4. **Exercise 4**: Run a parallel simulation on 2, 4, and 8 cores and compare performance

---

## Quick Reference

### Minimum Config File

```ini
[SIMULATION]
end_time = 86400.0
fluid_model = SINGLE_COMPONENT

[GRID]
nx = 20
ny = 20
nz = 5
Lx = 1000.0
Ly = 1000.0
Lz = 100.0

[ROCK]
porosity = 0.20
permeability_x = 100.0

[FLUID]
density = 1000.0
viscosity = 0.001
```

### Running Commands

```bash
# Basic run
./fsrm -c config.config

# Parallel run
mpirun -np N ./fsrm -c config.config

# With output directory
mpirun -np N ./fsrm -c config.config -o my_output/

# Verbose solver output
mpirun -np N ./fsrm -c config.config -snes_monitor -ksp_monitor
```

---

**Next Tutorial**: [Configuration System →](02_CONFIGURATION.md)
