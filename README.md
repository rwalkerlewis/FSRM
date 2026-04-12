# FSRM - Fully-coupled Subsurface Reservoir Model

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![FSRM CI](https://github.com/rwalkerlewis/FSRM/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/rwalkerlewis/FSRM/actions/workflows/ci.yml)
[![Static Analysis](https://github.com/rwalkerlewis/FSRM/actions/workflows/static-analysis.yml/badge.svg)](https://github.com/rwalkerlewis/FSRM/actions/workflows/static-analysis.yml)

A coupled multi-physics simulator for nuclear explosion monitoring, seismic wave propagation, dynamic fault rupture, and poroelasticity, built on PETSc/MPI for parallel unstructured FEM.

**MIT License** -- free to use, modify, and distribute for any purpose.

## Verified Capabilities (All Tests Pass)

The following features have automated tests with quantitative pass/fail criteria.

### Solid Mechanics
- **Elastostatics**: PetscFE pointwise callbacks (f0, f1, g3). Uniaxial compression converges in 1 SNES iteration. Unit tested.
- **Gravity body force**: Density-scaled body force callback. Unit tested and setup-verified.
- **Lithostatic stress**: Full FEM solve of 1D gravitational column. Analytical comparison of vertical and lateral stress profiles (K0 ratio). Passes within 5% tolerance.
- **Boundary conditions**: 6-face bounding box labeling, Dirichlet via DMAddBoundary, section rebuild. Verified with elastostatics.
- **Derived fields**: Cell-centered stress, strain, and Coulomb failure stress from FEM solution (DerivedFieldComputer).

### Wave Propagation
- **Elastodynamics**: Implicit time stepping with TS. Lamb's problem (point force on halfspace) and Garvin's problem (buried explosion Green's function) pass with quantitative error norms.
- **Absorbing boundary conditions**: Clayton-Engquist first-order. Energy absorption efficiency >99% at normal incidence. Unit tested.

### Poroelasticity
- **Biot coupling**: Pressure-displacement f0/f1 callbacks and all 4 Jacobian blocks (g0, g1, g2, g3). Unit tested.
- **Terzaghi consolidation**: 1D consolidation benchmark with analytical solution comparison. Passes with quantitative tolerance.

### Nuclear Explosion Monitoring
- **Mueller-Murphy seismic source**: Corner frequency, mb-yield scaling, reduced displacement potential (RDP) spectrum, cavity radius scaling, moment rate function. 9 tests pass.
- **Explosion seismogram pipeline**: Source injection via moment tensor, elastodynamic wave propagation, seismometer sampling (DMInterpolation), SAC output format. Integration-tested.
- **Moment tensor source injection**: Proper FEM equivalent nodal forces via PetscFECreateTabulation. Verified with nonzero displacement norm and correct solution vector size.
- **DPRK 2017 comparison**: Synthetic body-wave magnitude (mb) vs observed for 250 kt test at Punggye-ri. 4 tests pass.

### Atmospheric Explosion Effects
- **Sedov-Taylor blast**: Blast radius scaling (5% tolerance), t^0.4 temporal exponent. Tested.
- **Brode fireball**: Maximum radius estimation. Tested.
- **EMP E1 peak field**: 100-60000 V/m range validation. Tested.
- **Overpressure**: Distance-dependent overpressure (10-1000 kPa range). Tested.

### Near-field Explosion Phenomenology
- **Cavity radius**: Yield-dependent (W^0.295 scaling). Tested.
- **Damage zones**: Crushed/fractured/damaged ordering and extent ratios. Tested.
- **Spall**: Velocity (P/rho*c formula), thickness depth dependence, dynamic tensile strength rate effects. Tested.

### Fault Mechanics
- **Cohesive cell mesh splitting**: PyLith workflow (DMPlexCreateSubmesh, subpoint map, DMPlexLabelCohesiveComplete, DMPlexConstructCohesiveCells). 32 cohesive cells on 4x4x4 simplex mesh. All tests pass.
- **Friction laws**: Slip-weakening and rate-state (aging law). Unit tested.
- **Coulomb stress transfer**: Hooke stress computation, fault projection, delta-CFS calculation. Unit tested.
- **SCEC TPV5**: Infrastructure verified (parameters, CohesiveFaultKernel construction, FaultMeshManager). Full benchmark solve is a work in progress.

### Plasticity (Partial)
- **Yield function evaluation**: Drucker-Prager, von Mises, Mohr-Coulomb yield surface detection works correctly. Hydrostatic stress confirmed inside all surfaces.
- **Known limitation**: Return mapping algorithms do not produce nonzero plastic strain. Not wired into PETSc FEM pipeline.

### Fluid Flow Callbacks
- **Single-phase and black oil**: PetscFE pointwise callbacks. Unit tested. NOT verified end-to-end in a simulation.

## Planned / In Development

The following features exist as code stubs, partial implementations, or configuration templates but do **not** have passing end-to-end tests. Status labels indicate maturity.

| Feature | Status | Notes |
|---------|--------|-------|
| DG/ADER-DG spatial discretization | Stub | Header/source files exist. Not functional. |
| Local time stepping (LTS) | Stub | Not implemented. |
| PML absorbing boundaries | Stub | Only Clayton-Engquist is tested. |
| Viscoelastic attenuation | Stub | Not wired into FEM pipeline. |
| Off-fault plasticity (return mapping) | Broken | Yield detection works; return mapping does not produce plastic strain. |
| GPU acceleration (CUDA/HIP) | Stub | Build flags exist. No working GPU kernels. |
| Fourier Neural Operator (FNO) solver | Stub | Headers exist. Not functional. |
| Hydraulic fracturing (prototype) | Coded, unverified | Config exists. No end-to-end test. |
| Multi-phase fluid flow (end-to-end) | Coded, unverified | Callbacks unit-tested. No simulation test. |
| Injection source term | Coded, unverified | Config exists (injection_pressure_buildup.config). |
| Gmsh mesh import | Coded, unverified | Config stubs exist. Material region assignment untested. |
| Adaptive Mesh Refinement (AMR) | Stub | Config exists in aspirational. Not functional. |
| Volcano modeling | Stub | Documentation only. No physics implementation. |
| Tsunami modeling | Stub | Documentation only. No physics implementation. |
| Ocean circulation | Stub | Documentation only. No physics implementation. |
| Atmospheric infrasound | Stub | Documentation only. No physics implementation. |
| Radiation transport / fallout | Stub | Documentation only. No physics implementation. |
| Hypervelocity impacts | Stub | Documentation only. No physics implementation. |
| Advanced fracturing (multi-stage, proppant) | Stub | Not implemented. |
| Eclipse format I/O | Stub | Not tested. |
| Compositional EOS fluid | Stub | Not implemented. |
| Thermal coupling | Stub | Not implemented end-to-end. |
| Per-cell material properties | Not implemented | Single constants array for entire mesh. |

## Technology Stack

- **Language**: C++17
- **Build System**: CMake >= 3.15
- **Parallel Computing**: PETSc 3.22.2, MPI
- **FEM**: PETSc DMPlex unstructured finite elements with PetscDS pointwise callbacks
- **I/O**: HDF5
- **Testing**: Google Test, CTest
- **Containers**: Docker (Dockerfile.ci for reproducible builds)

## Building

### Docker (Recommended)

```bash
# Build the CI image
docker build -f Dockerfile.ci -t fsrm-ci:local .

# Build FSRM
docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:local bash -c \
  'mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release \
  -DENABLE_TESTING=ON -DENABLE_CUDA=OFF -DBUILD_EXAMPLES=ON && make -j$(nproc)'

# Run tests
docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:local \
  ctest --output-on-failure

# Interactive shell
docker run --rm -it -v $(pwd):/workspace -w /workspace fsrm-ci:local bash
```

### Native Build

Prerequisites: CMake >= 3.15, PETSc >= 3.15 (tested with 3.22.2), MPI, HDF5, C++17 compiler.

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON
make -j$(nproc)
```

## Running

The executable is `fsrm`. It reads an INI-style configuration file:

```bash
./fsrm config/examples/uniaxial_compression.config
```

Parallel execution:

```bash
mpirun -np 4 ./fsrm config/examples/explosion_seismogram.config
```

### Tested Configuration Files

These configs are in `config/examples/` and correspond to features with passing tests:

| Config File | Description |
|-------------|-------------|
| `uniaxial_compression.config` | Elastostatics benchmark |
| `terzaghi_consolidation.config` | 1D poroelastic consolidation |
| `explosion_seismogram.config` | Explosion source with seismometer output |
| `dprk_2017_quick.config` | DPRK 2017 synthetic test (quick run) |
| `lambs_problem.config` | Point force on halfspace (analytical benchmark) |
| `prescribed_slip_test.config` | Prescribed slip on cohesive fault |
| `fault_compression.config` | Fault under compression |
| `lithostatic_column.config` | Lithostatic stress column |
| `underground_explosion_template.config` | Underground explosion template |
| `minimal_explosion.config` | Minimal explosion test case |

### Aspirational Configuration Files

Configs in `config/aspirational/` reference features that are not yet functional. **Do not use them for production runs.** They serve as design targets for future development.

## Configuration File Format

FSRM uses INI-style configuration files with `[SECTION]` headers:

```ini
[SIMULATION]
start_time = 0.0
end_time = 86400.0        # 1 day in seconds
fluid_model = BLACK_OIL
enable_geomechanics = true

[ROCK]
porosity = 0.20           # 20%
permeability_x = 100.0    # milliDarcy
youngs_modulus = 10.0e9   # 10 GPa (Pa)

[FLUID]
oil_viscosity = 0.005     # Pa*s (5 cP)

[BC1]
type = DIRICHLET
field = PRESSURE
location = XMIN
value = 20.0e6            # Pa
```

See [docs/CONFIGURATION.md](docs/CONFIGURATION.md) for the full reference.

## Testing

```bash
cd build

# Run all tests
ctest --output-on-failure

# Run by label
ctest -L "unit"           # Unit tests
ctest -L "physics"        # Physics validation tests
ctest -L "integration"    # Integration tests

# Run specific test
ctest -R test_name -V
```

### Test Labels
- `unit`: Fast unit tests for individual components
- `functional`: Functional tests
- `physics`: Physics validation with analytical solutions or known benchmarks
- `integration`: Multi-component integration tests
- `performance`: Performance benchmarks

See [docs/TEST_RESULTS.md](docs/TEST_RESULTS.md) for the complete test matrix.

## Project Structure

```
FSRM/
  include/          # Header files (.hpp)
  src/              # Source files (.cpp)
  tests/            # Test files
    unit/           # Unit tests
    functional/     # Functional tests
    integration/    # Integration tests
    physics_validation/  # Physics/MMS tests
    performance/    # Performance tests
  config/
    examples/       # Tested configuration files
    aspirational/   # Untested design-target configs (DO NOT USE)
  docs/             # Documentation
  examples/         # Example executables
  external/         # External dependencies
```

## PETSc Options

Fine-tune solvers via command line:

```bash
-ts_type beuler          # Backward Euler
-snes_type newtonls      # Newton line search
-ksp_type gmres          # GMRES
-pc_type hypre           # Hypre preconditioner
-pc_hypre_type boomeramg # Algebraic multigrid
-snes_monitor            # Monitor convergence
-log_view                # Performance summary
```

## Contributing

Contributions welcome. Areas of interest:
- End-to-end verification of coded-but-unverified features
- Plasticity return mapping fix
- Gmsh mesh import testing
- Additional analytical benchmarks
- Per-cell material properties (auxiliary field pattern)
- Absorbing boundary improvements

## Citation

```bibtex
@software{fsrm2025,
  title = {FSRM: Coupled Multi-Physics Reservoir Simulator},
  author = {Robert Walker},
  year = {2025},
  version = {1.0.0}
}
```

## License

**MIT License** -- free to use, modify, and distribute. See [LICENSE](LICENSE).

## Documentation

- [Quick Start Guide](docs/QUICK_START.md)
- [User Guide](docs/USER_GUIDE.md)
- [Configuration Reference](docs/CONFIGURATION.md)
- [Development Guide](docs/DEVELOPMENT.md)
- [Physics Models](docs/PHYSICS_MODELS.md)
- [Numerical Methods](docs/NUMERICAL_METHODS.md)
- [API Reference](docs/API_REFERENCE.md)
- [Benchmarks](docs/BENCHMARKS.md)
- [Explosion & Impact Physics](docs/EXPLOSION_IMPACT_PHYSICS.md)
- [Unit System](docs/UNIT_SYSTEM.md)
- [Unstructured Meshes](docs/UNSTRUCTURED_MESHES.md)

## Acknowledgments

Built with:
- [PETSc](https://petsc.org/) (Portable, Extensible Toolkit for Scientific Computation)
- MPI (Message Passing Interface)
- HDF5 (Hierarchical Data Format)
