# FSRM Documentation

Reference documentation for FSRM (Full Service Reservoir Model). Navigate
by topic below; the standing references for current state are
[SOLVER_STATE.md](SOLVER_STATE.md), [HISTORIC_NUCLEAR_FIDELITY.md](HISTORIC_NUCLEAR_FIDELITY.md),
and [AXIS_1A_FIDELITY_REPORT.md](AXIS_1A_FIDELITY_REPORT.md).

## Getting started

| Document | What it covers |
|---|---|
| [QUICK_START.md](QUICK_START.md) | Five-minute build-and-run guide |
| [USER_GUIDE.md](USER_GUIDE.md) | End-to-end user manual; configs, runs, output |
| [DEVELOPMENT.md](DEVELOPMENT.md) | Building, testing, contributing |

## Configuration and API

| Document | What it covers |
|---|---|
| [CONFIGURATION.md](CONFIGURATION.md) | Config-file key reference |
| [FIDELITY_LADDER_GUIDE.md](FIDELITY_LADDER_GUIDE.md) | LOW / MED / HIGH / HIGHEST tier selection across all six ladders |
| [API_REFERENCE.md](API_REFERENCE.md) | Public C++ API |
| [CODE_INTERACTION_DIAGRAMS.md](CODE_INTERACTION_DIAGRAMS.md) | Module interaction diagrams |

## Physics and numerics

| Document | What it covers |
|---|---|
| [EXPLOSION_IMPACT_PHYSICS.md](EXPLOSION_IMPACT_PHYSICS.md) | Underground explosion source physics, per-pass detail |
| [PHYSICS_MODELS.md](PHYSICS_MODELS.md) | Constitutive models (Drucker-Prager, Tillotson, ANEOS, Mihalas-Mihalas) |
| [NUMERICAL_METHODS.md](NUMERICAL_METHODS.md) | Discretization choices (Lagrangian FV, Wilkins AV, multigroup diffusion, Strang, BDF2, RK3-SSP, sponge BC) |
| [WAVEFORM_VV.md](WAVEFORM_VV.md) | IRIS waveform comparison infrastructure |

## Verification

| Document | What it covers |
|---|---|
| [BENCHMARKS.md](BENCHMARKS.md) | Per-gate envelope status across all axes |
| [AXIS_1A_FIDELITY_REPORT.md](AXIS_1A_FIDELITY_REPORT.md) | Canonical cross-pass axis-1a result |
| [HISTORIC_NUCLEAR_FIDELITY.md](HISTORIC_NUCLEAR_FIDELITY.md) | Per-pass detail for the historic-nuclear program |
| [TEST_RESULTS.md](TEST_RESULTS.md) | `ctest` output snapshot |

## Roadmap

| Document | What it covers |
|---|---|
| [HISTORIC_NUCLEAR_ROADMAP.md](HISTORIC_NUCLEAR_ROADMAP.md) | Forward-looking six-axis roadmap |
| [AXIS_1B_DESIGN.md](AXIS_1B_DESIGN.md) | 3D source ball design (pass-13 implementation track) |
| [NUCLEAR_TEST_DIGITAL_TWIN.md](NUCLEAR_TEST_DIGITAL_TWIN.md) | Long-range vision document |

## Solver state

| Document | What it covers |
|---|---|
| [SOLVER_STATE.md](SOLVER_STATE.md) | Current fault-solver state, pass/fail, bottlenecks |
| [PYLITH_REFERENCE.md](PYLITH_REFERENCE.md) | Verified PyLith architecture pins |

## Mesh and coordinates

| Document | What it covers |
|---|---|
| [GMSH_MESH_GUIDE.md](GMSH_MESH_GUIDE.md) | Gmsh mesh authoring guide |
| [UNSTRUCTURED_MESHES.md](UNSTRUCTURED_MESHES.md) | DMPlex unstructured-mesh support |
| [COORDINATE_SYSTEMS.md](COORDINATE_SYSTEMS.md) | Geographic and local coordinates |
| [UNIT_SYSTEM.md](UNIT_SYSTEM.md) | Units and conversions |

## Sessions archive

[sessions/](sessions/) holds per-session work reports from the bring-up
and fault-solver development phases (sessions 2-32). They are point-in-
time and are not refreshed; prefer the standing references above for
current state.

## Archive

[archive/](archive/) holds retired physics docs that no longer describe
active work (NEM roadmap, fault-solver decision trail, PyLith feature-
parity wishlist).
