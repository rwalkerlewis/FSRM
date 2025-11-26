# FSRM Quick Start Guide

Get running in 5 minutes.

## Prerequisites

- CMake >= 3.15
- C++17 compiler (GCC 7+, Clang 6+)
- PETSc >= 3.15 with MPI
- HDF5 (optional)

## Installation

### 1. Install Dependencies (Ubuntu/Debian)

```bash
# Install PETSc and other dependencies
sudo apt update
sudo apt install -y cmake g++ libpetsc-dev libhdf5-mpi-dev

# Set environment
export PETSC_DIR=/usr/lib/petsc
export PETSC_ARCH=""
```

### 2. Build FSRM

```bash
# Clone and build
git clone https://github.com/your-repo/FSRM.git
cd FSRM
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j4

# Verify
./fsrm -help
```

## Your First Simulation

### Option 1: Use Example Configuration

```bash
# Run with pre-configured example
mpirun -np 4 ./fsrm -c ../config/default.config
```

### Option 2: Generate Your Own Configuration

```bash
# Generate a template
./fsrm -generate_config my_simulation.config

# Edit the configuration file
nano my_simulation.config

# Run your simulation
mpirun -np 4 ./fsrm -c my_simulation.config
```

### Option 3: Use Eclipse Input

```bash
# Run with Eclipse .DATA file
mpirun -np 4 ./fsrm -i SPE1.DATA -o output/spe1
```

## Configuration File Basics

Configuration files use a simple INI format:

```ini
[SIMULATION]
start_time = 0.0
end_time = 86400.0           # 1 day (in seconds)
fluid_model = SINGLE_COMPONENT

[GRID]
nx = 20
ny = 20
nz = 5
Lx = 1000.0                  # meters
Ly = 1000.0
Lz = 100.0

[ROCK]
porosity = 0.20
permeability_x = 100.0       # milliDarcy
youngs_modulus = 10.0e9      # Pa

[FLUID]
density = 1000.0             # kg/m³
viscosity = 0.001            # Pa·s

[WELL1]
name = PROD1
type = PRODUCER
i = 10
j = 10
k = 5
control_mode = RATE
target_value = 0.01          # m³/s
```

## Common Use Cases

### Single-Phase Flow
```bash
mpirun -np 4 ./fsrm -c ../config/default.config
```

### Shale Reservoir with Hydraulic Fracturing
```bash
mpirun -np 8 ./fsrm -c ../config/shale_reservoir.config
```

### Geothermal System
```bash
mpirun -np 8 ./fsrm -c ../config/geothermal.config
```

### Induced Seismicity Analysis
```bash
mpirun -np 8 ./fsrm -c ../config/induced_seismicity.config
```

### CO2 Storage
```bash
mpirun -np 8 ./fsrm -c ../config/co2_storage.config
```

## View Results

```bash
# Results are in VTK format by default
paraview output/*.vtu
```

## Next Steps

- Read the [User Guide](USER_GUIDE.md) for detailed documentation
- See [Configuration Reference](CONFIGURATION.md) for all options
- Check out [Physics Models](PHYSICS_MODELS.md) for theory
- For cloud deployment, see [Deployment Guide](DEPLOYMENT.md)

## Troubleshooting

### PETSc not found
```bash
# Check environment
echo $PETSC_DIR

# Try setting it manually
export PETSC_DIR=/path/to/petsc
```

### MPI errors
```bash
# Check MPI installation
which mpirun
mpirun --version

# Run on single process first
./fsrm -c config.config
```

### Convergence issues
```bash
# Reduce timestep
dt_initial = 100.0    # instead of 3600

# Increase iterations
max_nonlinear_iterations = 100
```
