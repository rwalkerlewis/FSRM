# FSRM - Full Service Reservoir Model

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![FSRM CI](https://github.com/rwalkerlewis/FSRM/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/rwalkerlewis/FSRM/actions/workflows/ci.yml)

FSRM is a coupled multiphysics simulator for nuclear explosion monitoring, seismic wave propagation, dynamic fault rupture, and THM poroelasticity. Built on PETSc 3.25.0 and MPI, it uses unstructured finite elements (DMPlex) with pointwise PetscDS callbacks for all PDE assembly. The simulator is config-driven: a single executable reads `.config` files that specify geometry, materials, physics, sources, and output.

**MIT License** -- free to use, modify, and distribute for any purpose.

---

## Quick Start

```bash
# Clone
git clone https://github.com/rwalkerlewis/FSRM.git && cd FSRM

# Build (Docker)
docker build -f Dockerfile.ci -t fsrm-ci:local .
docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:local bash -c \
  'mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON -DENABLE_CUDA=OFF && make -j$(nproc)'

# Run an example
docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:local \
  ./fsrm -c ../config/examples/uniaxial_compression.config

# Run tests
docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:local \
  ctest -j$(nproc) --output-on-failure

# Visualize seismogram output (host machine, requires Python)
pip install matplotlib obspy
python3 scripts/plot_seismograms.py build/output/seismograms/
```

---

## Examples

Eighteen runnable examples demonstrate verified capabilities. Each has a `README.md`, `run.sh`, and a config file.

| # | Example | Physics | Run Time |
|---|---------|---------|----------|
| 01 | [Uniaxial Compression](examples/01_uniaxial_compression/) | Elastostatics, Dirichlet BCs | < 1s |
| 02 | [Explosion Seismogram](examples/02_explosion_seismogram/) | Elastodynamics, Mueller-Murphy source, SAC output | ~2s |
| 03 | [Elastoplastic Compression](examples/03_elastoplastic_compression/) | Drucker-Prager plasticity | < 1s |
| 04 | [Locked Fault](examples/04_locked_fault/) | Cohesive cell insertion, locked constraint | < 1s |
| 05 | [Punggye-ri Nuclear Test](examples/05_punggye_ri_nuclear_test/) | Layered geology, explosion, seismograms | ~15s |
| 06 | [Gmsh Multi-Material](examples/06_gmsh_multimaterial/) | Gmsh mesh import, per-region materials | < 1s |
| 07 | [Traction BC](examples/07_traction_bc/) | Per-face Neumann traction BC, analytical verification | < 1s |
| 08 | [Time-Dependent Slip](examples/08_time_dependent_slip/) | Prescribed fault slip with linear time ramp | < 1s |
| 09 | [Gasbuggy 1967](examples/09_gasbuggy_1967/) | 29 kt, 4-layer Lewis Shale, SAC output | ~15s |
| 10 | [Gnome 1961](examples/10_gnome_1961/) | 3.1 kt, 4-layer Salado Salt, SAC output | ~15s |
| 11 | [Sedan 1962](examples/11_sedan_1962/) | 104 kt, 3-layer alluvium, SAC output | ~15s |
| 12 | [Degelen Mountain](examples/12_degelen_mountain/) | 50 kt, 3-layer granite, SAC output | ~15s |
| 13 | [NTS Pahute Mesa](examples/13_nts_pahute_mesa/) | 150 kt, 4-layer tuff, SAC output | ~15s |
| 14 | [Single-Phase Flow](examples/14_single_phase_flow/) | Darcy pressure diffusion, Dirichlet pressure BCs | < 1s |
| 15 | [Viscoelastic Attenuation](examples/15_viscoelastic_attenuation/) | Generalized Maxwell body, Q-factor, seismograms | ~15s |
| 16 | [SCEC TPV5](examples/16_scec_tpv5/) | Dynamic rupture, slip-weakening friction, nucleation | ~300s |
| 17 | [Velocity Model](examples/17_velocity_model/) | Per-cell material from binary velocity file (Vp/Vs/rho) | varies |
| 18 | [Thermal Expansion](examples/18_thermal_expansion/) | THM coupling, thermoelastic stress, uniform heating | < 1s |

---

## Verified Capabilities

Every feature below has automated tests with quantitative pass/fail criteria. Run `ctest --output-on-failure` to verify. 116 tests total.

### Integration-Tested Through TSSolve

| Feature | Test(s) | What It Proves |
|---------|---------|----------------|
| Elastostatics | `Physics.ElastostaticsPatch`, `Physics.LithostaticStress` | Hooke stress, patch test, K0 ratio |
| Elastodynamics | `Physics.LambsProblem`, `Physics.GarvinsProblem` | Wave propagation, analytical error norms |
| Poroelasticity | `Physics.TerzaghiConsolidation` | Biot coupling, analytical consolidation |
| Absorbing BCs | `Physics.AbsorbingBC` | Clayton-Engquist, >99% energy absorption |
| Gravity body force | `Physics.GravityLithostatic` | Lithostatic column, K0 within 5% |
| Moment tensor source | `Physics.MomentTensorSource` | FEM equivalent nodal forces |
| Explosion seismograms | `Integration.ExplosionSeismogram` | Source -> waves -> SAC output |
| Injection pressure | `Integration.InjectionPressure` | Poroelastic injection end-to-end |
| Depth-layered material | `Integration.LayeredElastostatics` | Aux field material assignment |
| Gmsh mesh import | `Integration.GmshImport` | MSH2 physical names, tet cells |
| Gmsh multi-material | `Integration.GmshMultiMaterial`, `Integration.GasbuggyMesh` | Per-label lambda/mu/rho |
| Gmsh nuclear twin | `Integration.NuclearTwinGmsh` | Mapped materials + explosion + HDF5 |
| Explosion damage zones | `Physics.ExplosionDamageZone` | Degraded aux near cavity |
| Punggye-ri layered | `Integration.PunggyeRiLayered` | 3-layer + absorbing + SAC + HDF5 |
| Pressurized fracture | `Integration.PressurizedFractureFEM` | Cohesive traction through TSSolve |
| Elastoplasticity | `Integration.ElastoplasticSim` | Drucker-Prager through TSSolve |
| Locked fault (quasi-static) | `Integration.DynamicRuptureSolve.LockedQuasiStatic` | Manual cohesive assembly |
| Locked fault (elastodynamic) | `Integration.DynamicRuptureSolve.LockedElastodynamic` | Cohesive + TSALPHA2 |
| Prescribed slip | `Integration.DynamicRuptureSolve.PrescribedSlip` | Imposed displacement jump |
| Locked fault transparency | `Physics.LockedFaultTransparency` | Fault slip < 5e-4 |
| Derived fields | `Integration.DerivedFields` | Stress/strain/CFS from solution |
| HDF5/VTK output | `Integration.OutputFile` | PetscViewerHDF5, VTK |
| Restart | `Integration.Restart` | Checkpoint/restore lifecycle |
| DPRK 2017 synthetic mb | `Integration.DPRK2017Comparison` | Synthetic vs observed body-wave magnitude |
| Explosion+fault residual | `Integration.ExplosionFaultReactivation` | Coexistence of moment-tensor and cohesive residual |
| NearField-FEM coupling | `Integration.NearFieldCoupled` | COUPLED_ANALYTIC 1D solver to 3D FEM moment rate |
| Per-face traction BC | `Integration.TractionBC` | Manual assembly, uniaxial analytical |
| Time-dependent slip ramp | `Integration.TimeDependentSlip` | Linear slip ramp with onset/rise time |
| Historic: Gasbuggy 1967 | `Integration.HistoricNuclear.Gasbuggy1967` | 29 kt, 4-layer Lewis Shale, SAC output |
| Historic: Gnome 1961 | `Integration.HistoricNuclear.Gnome1961` | 3.1 kt, 4-layer Salado Salt, SAC output |
| Historic: Sedan 1962 | `Integration.HistoricNuclear.Sedan1962` | 104 kt, 3-layer alluvium, SAC output |
| Historic: Degelen Mountain | `Integration.HistoricNuclear.DegelenMountain` | 50 kt, 3-layer granite, SAC output |
| Historic: NTS Pahute Mesa | `Integration.HistoricNuclear.NtsPahuteMesa` | 150 kt, 4-layer tuff, SAC output |
| Per-cell material from velocity model | `Integration.VelocityModelMaterial` | Binary Vp/Vs/rho grid mapped to mesh via interpolation |
| Thermal diffusion (heat equation) | `Integration.ThermalDiffusion` | Steady-state thermal field through TSSolve |
| Thermoelastic stress (THM) | `Integration.ThermalExpansion` | Thermal expansion coupling, uniform heating, displacement checks |

### Standalone Verified (Correct, Tested, Not FEM-Coupled)

| Feature | Test(s) |
|---------|---------|
| Mueller-Murphy source | `Physics.MuellerMurphy` |
| Near-field explosion (1D Lagrangian) | `Physics.NearFieldExplosion` |
| Atmospheric explosion (Sedov-Taylor, Brode, EMP) | `Physics.AtmosphericExplosion` |
| Drucker-Prager return mapping | `Unit.DruckerPragerStandalone` |
| Friction laws (slip-weakening, rate-state) | `Unit.FaultMechanics` |
| Coulomb stress transfer | `Unit.CoulombStressTransfer` |
| Hydrofrac formulas (Sneddon, PKN, Carter, Arps) | `Unit.HydrofracFormulas`, `Physics.StressShadowing`, etc. |
| Fluid flow callbacks | `Unit.SinglePhaseFlow`, `Unit.MultiphaseFlow` |

---

## Known Gaps

These entries track remaining limitations, work in progress, and recently closed gaps.

| Feature | Status |
|---------|--------|
| Multiphase flow end-to-end | Callbacks unit-tested; no simulation test |
| Hydraulic fracture coupled solve | PressurizedFractureFEM passes; lubrication+deformation not coupled |
| Explosion + fault full TSSolve | Test added (`ExplosionFaultReactivationTest.FullTSSolve` in `Integration.ExplosionFaultReactivation`); pending CI verification |
| Per-cell material from velocity model | DONE (`Integration.VelocityModelMaterial`, Example 17) |
| Thermal coupling | DONE (`Integration.ThermalExpansion`, `Integration.ThermalDiffusion`, Example 18) |
| Radiation transport / fallout | Source code archived; not integrated |

---

## Roadmap

Each item requires PetscDS callbacks integrated into `setupPhysics()`, integration tests through TSSolve, an example config, and visualization. Source code in `archive/src/` may provide a starting point but must be rewritten to use the PetscDS callback pattern.

1. ~~Slipping fault convergence~~ DONE
2. Multiphase flow end-to-end (Buckley-Leverett waterflood)
3. Full coupled hydraulic fracturing (lubrication + deformation)
4. ~~Viscoelastic attenuation~~ DONE
5. ~~Thermal coupling (heat equation + THM Biot)~~ DONE
6. Radiation transport (advection-diffusion for fallout)
7. ~~Per-cell material from velocity model files~~ DONE
8. ~~SCEC TPV5 dynamic rupture benchmark~~ DONE
9. Multi-stage hydraulic fracturing with stress shadowing
10. Production forecasting through propped fracture

---

## Technology Stack

| Component | Version | Role |
|-----------|---------|------|
| C++17 | GCC 11+ | Language standard |
| PETSc | 3.25.0 | FEM assembly, solvers, mesh (DMPlex) |
| MPI | OpenMPI 4+ | Parallelism |
| HDF5 | 1.10+ | Solution output |
| GTest | 1.14+ | Test framework |
| Docker | - | Build environment |

---

## Build Instructions

### Docker (Recommended)

```bash
docker build -f Dockerfile.ci -t fsrm-ci:local .
docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:local bash -c \
  'mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON -DENABLE_CUDA=OFF && make -j$(nproc)'
```

### Native (Requires PETSc 3.25.0)

```bash
export PETSC_DIR=/path/to/petsc-3.25.0
export PETSC_ARCH=arch-linux-c-opt
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON -DENABLE_CUDA=OFF
make -j$(nproc)
ctest -j$(nproc) --output-on-failure
```

### GPU Acceleration (PETSc CUDA)

FSRM supports GPU acceleration via PETSc native CUDA backend. No FSRM source changes needed. Build PETSc with `--with-cuda` and add runtime flags:

```bash
./fsrm -c config.config -vec_type cuda -mat_type aijcusparse -log_view
```

PETSc handles vector operations via cuBLAS, matrix operations via cuSPARSE, and KSP solves on GPU automatically.

---

## Repository Structure

```
src/                    Live source code (~45 files, ~35k lines)
include/                Headers
tests/                  116 automated tests (unit, functional, physics, integration)
config/examples/        Working example configurations
examples/               18 runnable examples with README and run.sh
scripts/                Visualization scripts (Python, reads C++ output)
meshes/                 Gmsh mesh files for examples
archive/                Removed dead code, fake examples, aspirational configs
  archive/src/          ~45k lines of dead/untested source code
  archive/examples/     80+ fake C++ example stubs
  archive/config/       141 aspirational configs for non-existent features
  archive/docs/         11 docs for non-existent features
```

---

## License

MIT License. See [LICENSE](LICENSE).
