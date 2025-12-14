# FSRM Documentation

This directory contains the complete documentation for FSRM (Fully-coupled Subsurface Reservoir Model).

## Documentation Structure

| File | Description |
|------|-------------|
| [QUICK_START.md](QUICK_START.md) | Get started in 5 minutes |
| [USER_GUIDE.md](USER_GUIDE.md) | Complete user manual |
| [CONFIGURATION.md](CONFIGURATION.md) | Configuration file reference |
| [CODE_INTERACTION_DIAGRAMS.md](CODE_INTERACTION_DIAGRAMS.md) | **How the code modules interact (diagrams)** |
| [PHYSICS_MODELS.md](PHYSICS_MODELS.md) | Physics and mathematical models |
| [PHYSICS_AND_GPU_ARCHITECTURE.md](PHYSICS_AND_GPU_ARCHITECTURE.md) | **Unified physics kernel and GPU architecture** |
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
| [EXPLOSION_IMPACT_PHYSICS.md](EXPLOSION_IMPACT_PHYSICS.md) | **Nuclear tests, impacts, EMP, radiation** |
| [VOLCANO_MODELING.md](VOLCANO_MODELING.md) | **Comprehensive volcanic system simulation** |
| [TSUNAMI_MODELING.md](TSUNAMI_MODELING.md) | **Earthquake-tsunami coupling and inundation** |

## Quick Links

### For New Users
- **Start Here**: [QUICK_START.md](QUICK_START.md)
- **Running Simulations**: [USER_GUIDE.md](USER_GUIDE.md)
- **Configuration Options**: [CONFIGURATION.md](CONFIGURATION.md)

### For Developers
- **Building**: [DEVELOPMENT.md](DEVELOPMENT.md)
- **Testing**: [BENCHMARKS.md](BENCHMARKS.md)
- **API Reference**: [API_REFERENCE.md](API_REFERENCE.md)
- **Code Interactions**: [CODE_INTERACTION_DIAGRAMS.md](CODE_INTERACTION_DIAGRAMS.md)

### By Topic
- **Physics**: [PHYSICS_MODELS.md](PHYSICS_MODELS.md)
- **System Architecture**: [CODE_INTERACTION_DIAGRAMS.md](CODE_INTERACTION_DIAGRAMS.md)
- **Physics Kernels & GPU**: [PHYSICS_AND_GPU_ARCHITECTURE.md](PHYSICS_AND_GPU_ARCHITECTURE.md) ⭐
- **Numerical Methods**: [NUMERICAL_METHODS.md](NUMERICAL_METHODS.md)
- **GPU Deployment**: [DEPLOYMENT.md](DEPLOYMENT.md#gpu-configuration)
- **Neural Networks**: [FOURIER_NEURAL_OPERATOR.md](FOURIER_NEURAL_OPERATOR.md)
- **Time Integration**: [IMEX_TRANSITION.md](IMEX_TRANSITION.md)
- **Units**: [UNIT_SYSTEM.md](UNIT_SYSTEM.md)
- **Coordinates**: [COORDINATE_SYSTEMS.md](COORDINATE_SYSTEMS.md)
- **Meshes**: [UNSTRUCTURED_MESHES.md](UNSTRUCTURED_MESHES.md), [GMSH_MESH_GUIDE.md](GMSH_MESH_GUIDE.md)
- **AMR**: [ADAPTIVE_MESH_REFINEMENT.md](ADAPTIVE_MESH_REFINEMENT.md)
- **Hydraulic Fracturing**: [RESFRAC_EQUIVALENT_CAPABILITIES.md](RESFRAC_EQUIVALENT_CAPABILITIES.md)
- **Eclipse Compatibility**: [ECLIPSE_FEATURES.md](ECLIPSE_FEATURES.md)
- **Explosions & Impacts**: [EXPLOSION_IMPACT_PHYSICS.md](EXPLOSION_IMPACT_PHYSICS.md) ⭐
- **Volcanoes**: [VOLCANO_MODELING.md](VOLCANO_MODELING.md) ⭐
- **Tsunamis**: [TSUNAMI_MODELING.md](TSUNAMI_MODELING.md) ⭐

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
8. **Explosions & Impacts**: Underground/atmospheric explosions, crater formation, EMP ⭐
9. **Volcanoes**: Magma chambers, conduit flow, eruption columns, PDCs, lava flows ⭐
10. **Tsunamis**: Shallow water equations, seafloor deformation, inundation ⭐
11. **Infrasound**: Atmospheric acoustic propagation ⭐

### Output Formats

- HDF5 (efficient large-scale data) - **Default**
- VTK (ParaView visualization)
- Eclipse (industry-compatible)

### GPU Acceleration

FSRM supports GPU acceleration for 5-50× speedup on large problems:

| Physics Kernel | GPU Support | Typical Speedup |
|----------------|-------------|-----------------|
| Single-Phase Flow | ✓ | 2-5× |
| Elastodynamics | ✓ | 5-20× |
| Poroelastodynamics | ✓ | 5-50× |
| Black Oil | ✓ | 3-8× |
| Thermal | ✓ | 2-5× |
| Geomechanics | ✓ | 3-10× |
| Hydrodynamic | ✓ | 10-50× |
| Tsunami | ✓ | 10-50× |
| Infrasound | ✓ | 5-30× |
| Atmospheric Blast | ✓ | 5-20× |
| Crater Formation | ✓ | 10-30× |
| Volcano (multi-physics) | ✓ | 5-20× |

**Supported Platforms:**
- NVIDIA CUDA (11.0+)
- AMD ROCm/HIP (5.0+)
- Multi-GPU with MPI
- Automatic CPU fallback

See:
- [PHYSICS_AND_GPU_ARCHITECTURE.md](PHYSICS_AND_GPU_ARCHITECTURE.md) - Kernel architecture
- [DEPLOYMENT.md](DEPLOYMENT.md#gpu-configuration) - GPU setup and configuration

## License

FSRM is open source software licensed under the [MIT License](../LICENSE).
