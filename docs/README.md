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

## Quick Links

- **New Users**: Start with [QUICK_START.md](QUICK_START.md)
- **Running Simulations**: See [USER_GUIDE.md](USER_GUIDE.md)
- **Configuration Options**: See [CONFIGURATION.md](CONFIGURATION.md)
- **Physics Background**: See [PHYSICS_MODELS.md](PHYSICS_MODELS.md)
- **Numerical Methods**: See [NUMERICAL_METHODS.md](NUMERICAL_METHODS.md)
- **Cloud Deployment**: See [DEPLOYMENT.md](DEPLOYMENT.md)

## Key Concepts

### Configuration-Driven Simulations

FSRM is designed so that **all simulation parameters can be specified in text configuration files**. No C++ code changes are needed to:

- Define rock and fluid properties
- Set up wells, fractures, and faults
- Configure physics models and solver options
- Specify boundary and initial conditions

Example:
```bash
# Run with a configuration file
mpirun -np 4 fsrm -c config/my_simulation.config
```

### Supported Physics

1. **Fluid Flow**: Single-phase, black oil, compositional
2. **Geomechanics**: Elastic, viscoelastic, poroelastic, plasticity
3. **Thermal**: Heat conduction and convection
4. **Fractures**: Natural and hydraulic, with propagation
5. **Faults**: Coulomb and rate-state friction, induced seismicity
6. **Dynamics**: Wave propagation, static triggering
7. **Coupling**: Full THM (Thermo-Hydro-Mechanical) coupling

### Output Formats

- VTK (ParaView visualization)
- HDF5 (efficient large-scale data)
- Eclipse (industry-compatible)
