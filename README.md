# FSRM - Full Service Reservoir Model

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![FSRM CI](https://github.com/rwalkerlewis/FSRM/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/rwalkerlewis/FSRM/actions/workflows/ci.yml)

FSRM is a coupled multiphysics simulator for nuclear-explosion monitoring,
seismic wave propagation, dynamic fault rupture, and THM poroelasticity.
Built on PETSc 3.25.0 with MPI; uses unstructured finite elements (DMPlex)
with pointwise PetscDS callbacks for all PDE assembly. The simulator is
config-driven: a single `fsrm` executable reads a `.config` text file that
specifies geometry, materials, physics, sources, and output.

The headline track is the historic-nuclear validation program: a fidelity
ladder of source physics, EOS, opacity, and time integration that is
verified gate-by-gate against published seismic observations from twenty-
seven historic underground nuclear tests. See
[docs/AXIS_1A_FIDELITY_REPORT.md](docs/AXIS_1A_FIDELITY_REPORT.md) for the
canonical cross-pass verification result and
[docs/HISTORIC_NUCLEAR_ROADMAP.md](docs/HISTORIC_NUCLEAR_ROADMAP.md) for the
forward roadmap.

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

# Run an example end-to-end (config lives next to the runner)
docker run --rm -v $(pwd):/workspace -w /workspace/examples/20_salmon_1964 fsrm-ci:local \
  ./run.sh

# Run all tests
docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:local \
  ctest -j$(nproc) --output-on-failure
```

For a guided walk through a first simulation, see
[docs/QUICK_START.md](docs/QUICK_START.md).

---

## Fidelity Ladders

The source-physics modules in `src/domain/explosion/` and the time-
integration strategies in `src/numerics/` ship as user-selectable LOW / MED
/ HIGH / HIGHEST tiers. Selection is per config block. Defaults are
backwards-compatible: a config that does not set a tier reproduces the
previous pass byte-for-byte.

| Knob | LOW | MED | HIGH | HIGHEST |
|---|---|---|---|---|
| `radiation_phase` | `ZELDOVICH_RAIZER` | `MARSHAK_GREY` | `MARSHAK_GREY` + `TABULATED_PATCHED` | `MARSHAK_MULTIGROUP` |
| `cavity_eos` | `IDEAL_GAS` | `TILLOTSON` | `TILLOTSON_TABULATED_PATCH` | `TABULATED_FULL` |
| `opacity_model` | `CONSTANT` | `POWER_LAW_ZR` | `TABULATED_PATCHED` | `TABULATED_FULL` |
| `operator_splitting` | `LIE` | `LIE` | `STRANG` | `STRANG_MULTIGROUP` |
| `time_integrator` (hydro) | `EXPLICIT_EULER` | `EXPLICIT_EULER` | `EXPLICIT_EULER` | `RK3_SSP` |
| `time_integrator_diffusion` | `BACKWARD_EULER` | `BACKWARD_EULER` | `CRANK_NICOLSON` | `BDF2` |

See [docs/FIDELITY_LADDER_GUIDE.md](docs/FIDELITY_LADDER_GUIDE.md) for the
single-page guide to picking a tier; see
[docs/EXPLOSION_IMPACT_PHYSICS.md](docs/EXPLOSION_IMPACT_PHYSICS.md) for
the per-pass physics description; see
[docs/AXIS_1A_FIDELITY_REPORT.md](docs/AXIS_1A_FIDELITY_REPORT.md) for the
gate-by-gate verification result.

---

## Examples

The `examples/` directory contains 38 runnable demonstrations of verified
capabilities. Each example is self-contained: `config.config`, `README.md`,
and `run.sh` live in the same directory; running `./run.sh` from the example
directory builds output under `output/`.

The historic-nuclear examples cover every era of underground testing
across all five weapon-state programs:

| US (NTS, Pacific, Alaska) | USSR (Semipalatinsk, Novaya Zemlya) | Other |
|---|---|---|
| 09 Gasbuggy 1967 (29 kt) | 30 Chagan 1965 (140 kt) | 12 Degelen Mountain (50 kt) |
| 10 Gnome 1961 (3.1 kt) | 31 Azgir A1 1966 (1.1 kt) | 32 Pokhran I 1974 (8 kt) |
| 11 Sedan 1962 (104 kt) | | 33-37 DPRK 2006-2016 (kt range) |
| 13 NTS Pahute Mesa (150 kt) | | 38 Lop Nor 1976 (Chinese, 4 Mt) |
| 19 Rainier 1957 (1.7 kt) | | |
| 20 Salmon 1964 (5.3 kt) | | |
| 21 Sterling 1966 (0.38 kt decoupled) | | |
| 22 Long Shot 1965 (80 kt, Amchitka) | | |
| 23 Milrow 1969 (1 Mt, Amchitka) | | |
| 24 Cannikin 1971 (~5 Mt, Amchitka) | | |
| 25 Faultless 1968 (1 Mt, Nevada) | | |
| 26 Baneberry 1970 (10 kt, NTS vent) | | |
| 27 Schooner 1968 (30 kt, Plowshare) | | |
| 28 Rulison 1969 (40 kt, gas stim) | | |
| 29 Rio Blanco 1973 (3 x 33 kt) | | |

Non-historic examples cover bring-up and feature verification:

| # | Example | Physics |
|---|---|---|
| 01 | Uniaxial Compression | Elastostatics, Dirichlet BCs |
| 02 | Explosion Seismogram | Elastodynamics, Mueller-Murphy source, SAC output |
| 03 | Elastoplastic Compression | Drucker-Prager plasticity |
| 04 | Locked Fault | Cohesive cell insertion, locked constraint |
| 05 | Punggye-ri Nuclear Test | Layered geology, explosion, seismograms |
| 06 | Gmsh Multi-Material | Gmsh mesh import, per-region materials |
| 07 | Traction BC | Per-face Neumann traction BC, analytical verification |
| 08 | Time-Dependent Slip | Prescribed fault slip with linear time ramp |
| 14 | Single-Phase Flow | Darcy pressure diffusion, Dirichlet pressure BCs |
| 15 | Viscoelastic Attenuation | Generalized Maxwell body, Q-factor, seismograms |
| 16 | SCEC TPV5 | Dynamic rupture, slip-weakening friction, nucleation |
| 17 | Velocity Model | Per-cell material from binary velocity file (Vp/Vs/rho) |
| 18 | Thermal Expansion | THM coupling, thermoelastic stress, uniform heating |

---

## Verified Capabilities

Every feature has automated tests with quantitative pass/fail criteria.
Run `ctest --output-on-failure` to verify. 116 tests are registered;
110 pass and 6 are documented honest failures (six fault-solver tests
blocked by a PETSc 3.25 BdResidual limitation; see
[docs/SOLVER_STATE.md](docs/SOLVER_STATE.md)). For the per-test detail
see [docs/TEST_RESULTS.md](docs/TEST_RESULTS.md) and CLAUDE.md
"Test Suite". Verified families:

- **Elasticity / elastodynamics**: patch test, Lithostatic stress,
  Lamb's problem, Garvin's problem.
- **Poroelasticity**: Terzaghi consolidation against analytical solution.
- **Absorbing boundaries**: Clayton-Engquist with > 99% energy absorption.
- **Cohesive faults**: locked fault transparency, prescribed slip,
  TSALPHA2 elastodynamic locked fault.
- **Dynamic rupture**: SCEC TPV5 with slip-weakening friction.
- **Source physics**: Mueller-Murphy RDP, near-field Lagrangian solver,
  multi-cell moment-tensor distribution with M0 conservation.
- **Historic-nuclear V&V**: 27 events from Gnome 1961 through DPRK 2017
  with pinned configs and quantitative far-field amplitude / onset /
  polarity / mb gates (see [HISTORIC_NUCLEAR_FIDELITY](docs/HISTORIC_NUCLEAR_FIDELITY.md)).
- **IRIS waveform V&V**: synthetic-vs-observed cross-correlation gates
  on Salmon 1964 with cached station data
  (see [docs/WAVEFORM_VV.md](docs/WAVEFORM_VV.md)).
- **Coupled THM**: heat equation, thermoelastic stress, full Biot.
- **Mesh I/O**: Gmsh per-label material assignment, binary velocity-model
  ingestion, HDF5 / VTK / SAC output.

---

## Pass-11 status

Pass-11 closed axis-1c (implicit time integration for the radiation
diffusion solve, Israeli-Orszag sponge layer for the 549 m free-field BC,
extended ANEOS Hugoniot match) and shipped the axis-1b 3D source ball
scaffolding ([docs/AXIS_1B_DESIGN.md](docs/AXIS_1B_DESIGN.md)). Pass-12
is housekeeping: documentation hygiene, config relocation into example
directories, the DPRK 2017 example, and presentation figures for the
showcase events. The next physics pass (pass-13) implements the 3D source
ball on the pass-11 scaffold.

---

## Technology Stack

| Component | Version | Role |
|---|---|---|
| C++17 | GCC 11+ | Language standard |
| PETSc | 3.25.0 | FEM assembly, solvers, mesh (DMPlex) |
| MPI | OpenMPI 4+ | Parallelism |
| HDF5 | 1.10+ | Solution output |
| GTest | 1.14+ | Test framework |
| Docker | - | Build environment |

---

## Build Instructions

### Docker (recommended)

```bash
docker build -f Dockerfile.ci -t fsrm-ci:local .
docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:local bash -c \
  'mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON -DENABLE_CUDA=OFF && make -j$(nproc)'
```

### Native (requires PETSc 3.25.0)

```bash
export PETSC_DIR=/path/to/petsc-3.25.0
export PETSC_ARCH=arch-linux-c-opt
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON -DENABLE_CUDA=OFF
make -j$(nproc)
ctest -j$(nproc) --output-on-failure
```

### GPU acceleration (PETSc CUDA)

FSRM does not require source changes for GPU. Build PETSc with `--with-cuda`
(see `Dockerfile.cuda`) and add runtime flags:

```bash
./fsrm -c config.config -vec_type cuda -mat_type aijcusparse -log_view
```

PETSc handles vector operations via cuBLAS, matrix operations via cuSPARSE,
and KSP solves on GPU. The PetscDS pointwise callbacks remain on CPU; data
transfer is automatic.

---

## Repository Structure

```
src/                  Live source code (~54 kloc across 63 .cpp files)
include/              Headers (~34 kloc across 73 .hpp files)
tests/                116 automated tests (unit, functional, physics, integration)
config/               Schema-anchor templates (default, complete_template,
                      test_*); per-event configs live next to their examples
examples/             38 runnable examples with config.config, README.md,
                      and run.sh
scripts/              Visualization scripts and waveform fetchers
                      (Python; read C++ output)
tools/                Build-time tooling (figure styles, tabulated EOS/opacity
                      data refresh, waveform V&V infrastructure)
meshes/               Gmsh mesh files for examples
docs/                 Standing reference docs; sessions/ archive of per-
                      session reports; archive/ for retired physics docs
archive/              Removed dead code (~45 kloc) preserved for context
```

---

## License

MIT License. See [LICENSE](LICENSE).
