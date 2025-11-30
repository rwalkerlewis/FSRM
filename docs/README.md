# FSRM Documentation

This directory contains the complete documentation for FSRM (Fully-coupled Subsurface Reservoir Model).

## Documentation Structure

| File | Description |
|------|-------------|
| [QUICK_START.md](QUICK_START.md) | Get started in 5 minutes |
| [USER_GUIDE.md](USER_GUIDE.md) | Complete user manual |
| [CONFIGURATION.md](CONFIGURATION.md) | Configuration file reference |
| [PHYSICS_MODELS.md](PHYSICS_MODELS.md) | Physics and mathematical models |
| [NUMERICAL_METHODS.md](NUMERICAL_METHODS.md) | Numerical approaches and algorithms |
| [API_REFERENCE.md](API_REFERENCE.md) | C++ API documentation |
| [DEPLOYMENT.md](DEPLOYMENT.md) | Cloud, HPC, and GPU deployment |
| [DEVELOPMENT.md](DEVELOPMENT.md) | Building, testing, and contributing |
| [BENCHMARKS.md](BENCHMARKS.md) | Validation and benchmark results |
| [UNIT_SYSTEM.md](UNIT_SYSTEM.md) | Unit conversion and specifications |
| [COORDINATE_SYSTEMS.md](COORDINATE_SYSTEMS.md) | Geographic coordinate transformations |
| [UNSTRUCTURED_MESHES.md](UNSTRUCTURED_MESHES.md) | Gmsh and unstructured mesh support |
| [GMSH_MESH_GUIDE.md](GMSH_MESH_GUIDE.md) | Complete guide to creating Gmsh meshes |
| [ADAPTIVE_MESH_REFINEMENT.md](ADAPTIVE_MESH_REFINEMENT.md) | AMR for unstructured grids |
| [ECLIPSE_FEATURES.md](ECLIPSE_FEATURES.md) | Eclipse reservoir simulator compatibility |
| [FOURIER_NEURAL_OPERATOR.md](FOURIER_NEURAL_OPERATOR.md) | Neural network PDE solver |
| [IMEX_TRANSITION.md](IMEX_TRANSITION.md) | Implicit-Explicit time integration |
| [RESFRAC_EQUIVALENT_CAPABILITIES.md](RESFRAC_EQUIVALENT_CAPABILITIES.md) | Hydraulic fracturing capabilities |
| [SEISSOL_COMPARISON.md](SEISSOL_COMPARISON.md) | Comparison with SeisSol earthquake simulator |

## Quick Links

### For New Users
- **Start Here**: [QUICK_START.md](QUICK_START.md)
- **Running Simulations**: [USER_GUIDE.md](USER_GUIDE.md)
- **Configuration Options**: [CONFIGURATION.md](CONFIGURATION.md)

### For Developers
- **Building**: [DEVELOPMENT.md](DEVELOPMENT.md)
- **Testing**: [BENCHMARKS.md](BENCHMARKS.md)
- **API Reference**: [API_REFERENCE.md](API_REFERENCE.md)

### By Topic
- **Physics**: [PHYSICS_MODELS.md](PHYSICS_MODELS.md)
- **Numerical Methods**: [NUMERICAL_METHODS.md](NUMERICAL_METHODS.md)
- **GPU Acceleration**: [DEPLOYMENT.md](DEPLOYMENT.md#gpu-deployment)
- **Neural Networks**: [FOURIER_NEURAL_OPERATOR.md](FOURIER_NEURAL_OPERATOR.md)
- **Time Integration**: [IMEX_TRANSITION.md](IMEX_TRANSITION.md)
- **Units**: [UNIT_SYSTEM.md](UNIT_SYSTEM.md)
- **Coordinates**: [COORDINATE_SYSTEMS.md](COORDINATE_SYSTEMS.md)
- **Meshes**: [UNSTRUCTURED_MESHES.md](UNSTRUCTURED_MESHES.md), [GMSH_MESH_GUIDE.md](GMSH_MESH_GUIDE.md)
- **AMR**: [ADAPTIVE_MESH_REFINEMENT.md](ADAPTIVE_MESH_REFINEMENT.md)
- **Hydraulic Fracturing**: [RESFRAC_EQUIVALENT_CAPABILITIES.md](RESFRAC_EQUIVALENT_CAPABILITIES.md)
- **Eclipse Compatibility**: [ECLIPSE_FEATURES.md](ECLIPSE_FEATURES.md)

## Key Concepts

### Configuration-Driven Simulations

FSRM is designed so that **all simulation parameters can be specified in text configuration files**. No C++ code changes are needed to:

- Define rock and fluid properties
- Set up wells, fractures, and faults
- Configure physics models and solver options
- Specify boundary and initial conditions
- Select input/output units

Example:
```bash
# Run with a configuration file
mpirun -np 4 fsrm -c config/my_simulation.config
```

### Supported Physics

1. **Fluid Flow**: Single-phase, black oil, compositional
2. **Geomechanics**: Elastic, viscoelastic, poroelastic, plasticity
3. **Thermal**: Heat conduction and convection
4. **Fractures**: Natural fractures (DFN), hydraulic fracturing (PKN/KGD/P3D)
5. **Faults**: Coulomb and rate-state friction, induced seismicity
6. **Dynamics**: Wave propagation, static triggering, dynamic rupture
7. **Coupling**: Full THM (Thermo-Hydro-Mechanical) coupling

### Output Formats

- HDF5 (efficient large-scale data) - **Default**
- VTK (ParaView visualization)
- Eclipse (industry-compatible)

### GPU Acceleration

FSRM supports GPU acceleration for 5-50x speedup on large problems:
- NVIDIA CUDA (11.0+)
- AMD ROCm/HIP (5.0+)
- Multi-GPU with MPI
- Automatic CPU fallback

See [DEPLOYMENT.md](DEPLOYMENT.md) for details.
