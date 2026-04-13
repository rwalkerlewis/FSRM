# FSRM - Fully-coupled Subsurface Reservoir Model

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![FSRM CI](https://github.com/rwalkerlewis/FSRM/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/rwalkerlewis/FSRM/actions/workflows/ci.yml)
[![Static Analysis](https://github.com/rwalkerlewis/FSRM/actions/workflows/static-analysis.yml/badge.svg)](https://github.com/rwalkerlewis/FSRM/actions/workflows/static-analysis.yml)

A coupled multi-physics simulator for nuclear explosion monitoring, seismic wave propagation, dynamic fault rupture, and poroelasticity, built on PETSc/MPI for parallel unstructured FEM.

**MIT License** -- free to use, modify, and distribute for any purpose.

## Verified Capabilities (91 Docker Tests Pass)

The following features have automated tests with quantitative pass/fail criteria. Every claim below is backed by a specific test name. Run `ctest --output-on-failure` to verify.

### Solid Mechanics
- **Elastostatics**: PetscFE pointwise callbacks (f0, f1, g3). Uniaxial compression converges in 1 SNES iteration. Tests: `Physics.ElastostaticsPatch`, `Physics.LithostaticStress`, `Integration.FullSimulation`.
- **Gravity body force**: Density-scaled body force callback. Lithostatic stress column with K0 ratio within 5% tolerance. Tests: `Physics.GravityLithostatic`.
- **Lithostatic stress**: Full FEM solve of 1D gravitational column. Analytical comparison of vertical and lateral stress profiles. Tests: `Physics.LithostaticStress`.
- **Boundary conditions**: 6-face bounding box labeling, Dirichlet via DMAddBoundary, section rebuild. Tests: `Unit.BoundaryConditions`.
- **Derived fields**: Cell-centered stress, strain, and Coulomb failure stress from FEM solution. Tests: `Integration.DerivedFields`.
- **Layered elastostatics**: Depth-based material layering via auxiliary fields. Tests: `Integration.LayeredElastostatics`.

### Wave Propagation
- **Elastodynamics**: Implicit time stepping with TSALPHA2. Tests: `Physics.LambsProblem` (point force on halfspace), `Physics.GarvinsProblem` (buried explosion Green function).
- **Absorbing boundary conditions**: Clayton-Engquist first-order Lysmer-Kuhlmeyer. Energy absorption >99% at normal incidence. Tests: `Physics.AbsorbingBC`.

### Poroelasticity
- **Biot coupling**: Pressure-displacement f0/f1 callbacks and all 4 Jacobian blocks (g0, g1, g2, g3). Tests: `Physics.TerzaghiConsolidation`.
- **Injection pressure buildup**: Poroelastic Simulator with point injection source. Two-simulation comparison. Tests: `Integration.InjectionPressure`.

### Nuclear Explosion Monitoring
- **Mueller-Murphy seismic source**: Corner frequency, mb-yield scaling, reduced displacement potential (RDP) spectrum, cavity radius scaling, moment rate function. Tests: `Physics.MuellerMurphy` (9 assertions).
- **Explosion seismogram pipeline**: Source injection via moment tensor, elastodynamic wave propagation, seismometer sampling (DMInterpolation), SAC output format. Tests: `Integration.ExplosionSeismogram`.
- **Punggye-ri layered workflow**: Three-layer elastic model with absorbing boundaries, damage-zone degradation, HDF5 output, and SAC seismograms. Tests: `Integration.PunggyeRiLayered`.
- **Moment tensor source injection**: Proper FEM equivalent nodal forces via PetscFECreateTabulation. Tests: `Physics.MomentTensorSource`.
- **DPRK 2017 comparison**: Synthetic body-wave magnitude (mb) vs observed for 250 kt test at Punggye-ri. Tests: `Integration.DPRK2017Comparison` (4 assertions).

### Atmospheric Explosion Effects
- **Sedov-Taylor blast**: Blast radius scaling (5% tolerance), t^0.4 temporal exponent. Tests: `Physics.AtmosphericExplosion`.
- **Brode fireball**: Maximum radius estimation. Tests: `Physics.AtmosphericExplosion`.
- **EMP E1 peak field**: 100-60000 V/m range validation. Tests: `Physics.AtmosphericExplosion`.
- **Overpressure**: Distance-dependent (10-1000 kPa range). Tests: `Physics.AtmosphericExplosion`.

### Near-field Explosion Phenomenology
- **1D Lagrangian solver**: Mie-Gruneisen EOS, damage, spall. Tests: `Physics.NearFieldExplosion`.
- **FEM-coupled damage-zone degradation**: Auxiliary material properties reduced around explosion source, dynamic response changes with damage. Tests: `Physics.ExplosionDamageZone`.

### Gmsh Multi-Material Support
- **Gmsh physical-name import**: PETSc DMPlex loads MSH2 meshes and preserves physical-group labels. Tests: `Integration.GmshImport`.
- **Per-cell material assignment by Gmsh label**: Auxiliary fields populated from `[MATERIAL_REGION_N]` sections keyed by `gmsh_label`. Tests: `Integration.GmshMultiMaterial`.
- **Historical mesh validation**: Gasbuggy layered mesh with three distinct material regions. Tests: `Integration.GasbuggyMesh`.
- **Gmsh nuclear twin workflow**: Compact Gmsh mesh with mapped material regions, underground source, and HDF5 output. Tests: `Integration.NuclearTwinGmsh`.

### Fault Mechanics
- **Cohesive cell mesh splitting**: PyLith workflow (DMPlexCreateSubmesh, subpoint map, DMPlexLabelCohesiveComplete, DMPlexConstructCohesiveCells). Tests: `Functional.DynamicRuptureSetup.MeshSplitting`.
- **Friction laws**: Slip-weakening and rate-state (aging law). Tests: `Unit.FaultMechanics`.
- **Coulomb stress transfer**: Hooke stress computation, fault projection, delta-CFS calculation. Tests: `Unit.CoulombStressTransfer`.
- **Fault plus absorbing BC coexistence**: Region-specific PetscDS assignment via `DMSetRegionDS`. Tests: `Functional.FaultAbsorbingCoexist`, `Functional.DynamicRuptureSetup.AbsorbingCoexist`.
- **SCEC TPV5 infrastructure**: Parameters, CohesiveFaultKernel construction, FaultMeshManager verified. Full benchmark solve is a work in progress. Tests: `Physics.SCEC.TPV5`.
- **Prescribed slip**: Cohesive fault with imposed slip. Tests: `Integration.PrescribedSlip`.
- **Pressurized fracture FEM**: Traction balance with Sneddon aperture validation through TSSolve. Tests: `Integration.PressurizedFractureFEM`.

### Hydraulic Fracturing (Standalone Utilities)
- **PKN width scaling**: Quarter-power width validated for coupled flow-deformation. Tests: `Unit.HydrofracFormulas`.
- **Fracture propagation criterion**: LEFM cohesive strength mapping from K_Ic and opening threshold. Tests: `Unit.FracturePropagation`.
- **Lubrication flow callbacks**: Fracture pressure source and Poiseuille flux with aperture-cubed scaling. Tests: `Unit.FractureFlowCallbacks`.
- **Pressurized fracture callbacks**: Traction balance and Sneddon validation. Tests: `Unit.PressurizedFractureCallbacks`.
- **Multi-cluster stress shadowing**: Sneddon normal-stress perturbation, shadow factor decay, cluster efficiency allocation. Tests: `Physics.StressShadowing` (6 assertions).
- **Induced seismicity**: Scalar moment (M0 = mu * slip * area), Hanks-Kanamori magnitude, microseismic range validation. Tests: `Physics.InducedSeismicity` (6 assertions).
- **Proppant transport**: Stokes settling, bridging criterion, pack minimum aperture, mass conservation. Tests: `Physics.ProppantTransport` (8 assertions).
- **Carter leak-off coupling**: Carter leak-off rate, cumulative volume, delayed opening, area scaling. Tests: `Physics.LeakoffCoupling` (7 assertions).
- **Production forecasting**: Arps hyperbolic/exponential/harmonic decline curves, cumulative production, fracture productivity index. Tests: `Physics.ProductionForecast` (9 assertions).

### Plasticity (Material-Point Level)
- **Drucker-Prager return mapping**: PlasticityModel::integrateStress produces nonzero plastic strain for deviatoric loading. Tests: `Unit.DruckerPragerStandalone`.
- **Yield function evaluation**: Drucker-Prager, von Mises, Mohr-Coulomb yield surface detection. Tests: `Unit.Elastoplasticity`.
- **Known limitation**: Return mapping works at the material-point level only. Not wired into PETSc FEM pipeline (no config flag to enable).

### Fluid Flow Callbacks
- **Single-phase and black oil**: PetscFE pointwise callbacks. Tests: `Unit.SinglePhaseFlow`, `Unit.MultiphaseFlow`. NOT verified end-to-end in a simulation.

### Output
- **HDF5/VTK output**: Solution and derived fields written to HDF5. Tests: `Integration.OutputFile`.
- **Simulation restart**: Checkpoint and restart lifecycle. Tests: `Integration.Restart`.

## What Does Not Work (Known Gaps)

These features exist as code but are NOT functional end-to-end:

| Feature | Status | Notes |
|---------|--------|-------|
| Elastoplasticity FEM coupling | Not wired | Return mapping works standalone; no config flag to enable in Simulator |
| Hydraulic fracturing (full coupled) | Partial | PressurizedFractureFEM passes; lubrication+deformation coupled solve is callback-tested only |
| Dynamic rupture TSSolve | Setup only | Pipeline setup completes; TSSolve never called. SNES may diverge |
| Fault + absorbing coexistence solve | Setup only | Setup succeeds; TSSolve never called |
| End-to-end multiphase flow | Not tested | Callbacks unit-tested; no simulation test |
| DG/ADER-DG spatial discretization | Stub | Source files exist. Not functional |
| GPU acceleration (CUDA/HIP) | Stub | Build flags exist. No working GPU kernels |
| Fourier Neural Operator (FNO) | Stub | Headers exist. Not functional |
| Adaptive Mesh Refinement (AMR) | Stub | Not functional |
| Volcano modeling | Dead code | ~4,200 lines. Never referenced. Never tested |
| Tsunami modeling | Dead code | ~5,000 lines. Never referenced. Never tested |
| Ocean circulation | Dead code | Part of tsunami/ocean module. Never tested |
| Radiation transport | Dead code | ~1,650 lines. Never referenced. Never tested |
| ML/neural operators | Dead code | ~9,200 lines. Never referenced. Never tested |
| Material library | Dead code | ~7,700 lines. Never referenced. Never tested |

### Dead Code Disclosure

Approximately 44,000 lines (~53% of source .cpp files) compile into the library but are never called by the Simulator or any test. This includes stubs for volcano, tsunami, ocean, ML, DG/ADER, radiation, AMR, and an unused material library. See CLAUDE.md for the complete inventory.

## Technology Stack

- **Language**: C++17
- **Build System**: CMake >= 3.15
- **Parallel Computing**: PETSc 3.22.2, MPI
- **FEM**: PETSc DMPlex unstructured finite elements with PetscDS pointwise callbacks
- **I/O**: HDF5
- **Testing**: Google Test, CTest (91 tests across 6 executables)
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

# Run all tests
docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:local \
  ctest -j$(nproc) --output-on-failure

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

# Run all tests in parallel
ctest -j$(nproc) --output-on-failure

# Run by label
ctest -L "unit"                 # 36 unit tests
ctest -L "functional"           # 10 functional tests
ctest -L "physics_validation"   # 23 physics validation tests
ctest -L "integration"          # 18 integration tests
ctest -L "performance"          # 3 performance tests

# Run specific test
ctest -R Physics.TerzaghiConsolidation -V
```

### Test Labels
- `unit`: Standalone formula, callback, and component tests (36 tests)
- `functional`: Setup pipeline verification, no TSSolve (10 tests)
- `physics_validation`: Analytical solutions, FEM-coupled benchmarks (23 tests)
- `integration`: Full Simulator pipeline through TSSolve (18 tests)
- `performance`: Performance benchmarks (3 tests)
- `experimental`: Neural operator stubs (1 test)

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
- End-to-end verification of fluid flow callbacks (simulation test for single-phase/multiphase)
- Wiring elastoplasticity into Simulator (PetscDS callback registration + config flag)
- Dynamic rupture TSSolve convergence (cohesive Jacobian debugging)
- Additional analytical benchmarks
- Absorbing boundary improvements (PML)

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
