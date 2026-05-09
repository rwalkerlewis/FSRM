# Benchmarks

The FSRM verification gates fall into three categories: analytical
solutions, historic-nuclear V&V against published observations, and
SCEC dynamic rupture. This document is the gate-by-gate detail behind
the cross-pass summary in
[AXIS_1A_FIDELITY_REPORT.md](AXIS_1A_FIDELITY_REPORT.md). For per-pass
narrative, see [HISTORIC_NUCLEAR_FIDELITY.md](HISTORIC_NUCLEAR_FIDELITY.md).

## Analytical solutions

Each gate has a quantitative `EXPECT_*` assertion against an analytical
or semi-analytical reference. Tests live under
`tests/physics_validation/`.

| Test | Reference | Tolerance |
|---|---|---|
| `Physics.ElastostaticsPatch` | Patch test, exact Hooke stress | 1e-12 nodal |
| `Physics.LithostaticStress` | Closed-form K0 ratio | 5 % |
| `Physics.GravityLithostatic` | Lithostatic column | 5 % K0 |
| `Physics.LambsProblem` | Lamb 1904 surface wave | L2 norm 5 % |
| `Physics.GarvinsProblem` | Garvin 1956 buried explosion | L2 norm 5 % |
| `Physics.TerzaghiConsolidation` | Terzaghi 1925 1D consolidation | 2 % time-evolution |
| `Physics.AbsorbingBC` | Clayton-Engquist energy flux | > 99 % absorbed |
| `Physics.MomentTensorSource` | Aki-Richards 1980 equivalent body force | 5 % node-load |
| `Unit.DruckerPragerStandalone` | Return-mapping fixture | 1e-10 stress |
| `Unit.HydrofracFormulas` | Sneddon 1951, PKN, Carter, Arps | 1e-6 |
| `Physics.LockedFaultTransparency` | Kinematic transparency | slip < 5e-4 |
| `Physics.ViscoelasticRelaxation` | GMB closed-form decay | 5 % |
| `Physics.CohesiveBdResidual` | PetscDS BdResidual on cohesive | 1e-12 |
| `Physics.SCEC.TPV5` | SCEC TPV5 rupture-front benchmark | matches PyLith ref |

## Historic-nuclear V&V (axis 1a)

Per-event quantitative gates against published seismic observations and
near-field measurements. Pinned tier configs ensure pass-N reproducibility.

| Gate | Pass-7 | Pass-8 | Pass-9 | Pass-10 | Pass-11 | Spec target | Residual / next |
|---|---|---|---|---|---|---|---|
| Salmon CavityRadius | f5 | f4 | f3 | f3 | f3 | 5 % | axis-1b 3D source ball |
| Salmon FreeFieldPeakVelocity_166m | f5 | f4 | **f2** | f2 | f2 | f2 | closed |
| Salmon FreeFieldPeakVelocity_322m | f6 | f4 | **f2** | f2 | f2 | f2 | closed |
| Salmon FreeFieldPeakVelocity_549m | skip | skip | skip | f4 | f4 | f2 | impedance BC contributes; sponge BC opt-in. Closure requires spec sweep or axis-1d |
| Salmon FarFieldBodyWaveMagnitude | +/-0.5 | +/-0.3 | +/-0.3 | +/-0.3 | +/-0.3 | +/-0.2 | propagation-path drift, axis-3 |
| Marshak SelfSimilarPureRadiation | n/a | f10 | f8 | f2.5 | **f2** (BDF2) | f2 | closed (HIGHEST tier under axis-1c) |
| Marshak RadiationEnergyConservation | n/a | 25 % | 10 % | 10 % | **2 %** (BDF2) | 2 % | closed (HIGHEST tier under axis-1c) |
| Marshak GreyVsZRComparison | n/a | f5 | f3 | f3 | f3 | f3 | closed |
| Chagan CavityRadius | f6 | f5 | f5 | f5 | f5 | f3 | axis-1b 3D source ball |
| Chagan FarFieldBodyWaveMagnitude | +/-0.5 | +/-0.4 | +/-0.4 | +/-0.4 | +/-0.4 | +/-0.3 | propagation-path drift, axis-3 |
| PokhranI FarFieldBodyWaveMagnitude | +/-0.5 | +/-0.4 | +/-0.4 | +/-0.4 | +/-0.4 | +/-0.3 | regional Murphy 1981 reference drift, axis-4 |
| Granite Hugoniot match (Marsh 1980) | n/a | n/a | n/a | n/a | **30 %, half within** | 5 % | Tillotson refit, axis-4 |
| Salt Hugoniot match (McQueen 1970) | n/a | n/a | n/a | n/a | **30 %, half within** | 5 % | Tillotson refit, axis-4 |

Notation: `f<N>` denotes the achieved peak-amplitude / cavity-radius
envelope factor. `+/-x` denotes the body-wave magnitude (mb) tolerance.
Bold cells denote the pass that closed the gate at spec target.

The pass-by-pass narrative including which gates were retightened, what
work closed them, and what residuals remain is in
[HISTORIC_NUCLEAR_FIDELITY.md](HISTORIC_NUCLEAR_FIDELITY.md).

## IRIS waveform V&V

Synthetic-vs-observed waveform comparison anchored on Salmon 1964 with
Chagan 1965 and Pokhran I 1974 cross-validation. ObsPy refresh tool at
`tools/waveform_vv/refresh.py`. Cached under `tools/waveform_vv/cache/`.

| Gate | Reference | Threshold |
|---|---|---|
| Cross-correlation (0.5-5 Hz) | Cached IRIS waveforms | > 0.7 |
| Peak amplitude envelope | Cached IRIS waveforms | factor 2 |
| Arrival-time delta | Cached IRIS waveforms | < 1 sample at closest station |

CTest label `iris_validation` runs five integration gates that
GTEST_SKIP cleanly when the cache is empty.

## Multi-physics integration

| Test | What it proves |
|---|---|
| `Integration.NearFieldCoupled` | 1D Lagrangian solver to 3D FEM moment-rate handoff |
| `Integration.HistoricNuclear.{Gasbuggy,Gnome,Sedan,Degelen,Pahute}` | End-to-end pinned-config historic events |
| `Integration.HistoricNuclear.FarFieldAmplitudeRegression` | SINGLE_CELL vs UNIFORM_SPHERE distribution regression CSV |
| `Integration.SourceDistribution.{SingleCellLegacyByteIdentical, GaussianM0Conserved, ...}` | M0 conservation, fallback-to-SINGLE_CELL safety |
| `Integration.MPI.SalmonSerialVsParallelEquivalence` | Pass-12 parallel correctness gate |

## SCEC TPV5 dynamic rupture

`Physics.SCEC.TPV5` runs the SCEC TPV5 benchmark with slip-weakening
friction, initial fault stress, and a nucleation patch. Quantitative
matches against the published reference solution at the rupture-front
arrival time and slip-rate magnitude.

## Performance benchmarks

`Performance.{Benchmark, Scaling, Memory, GPU}` are smoke-level
performance tests; they verify that solvers complete in expected
wall-time orders, memory does not grow unbounded, and GPU paths
respond to PETSc CUDA flags. Strong/weak scaling studies are not in
the standard CI label set; see CLAUDE.md "Test Suite" for the
breakdown.

## How to run

```bash
# All benchmarks
ctest -j$(nproc) --output-on-failure

# By label
ctest -L unit
ctest -L physics_validation
ctest -L integration
ctest -L iris_validation
ctest -L performance

# Single test
ctest -R Physics.LambsProblem --output-on-failure
```

## Updating benchmarks

When a gate envelope is tightened, edit
[AXIS_1A_FIDELITY_REPORT.md](AXIS_1A_FIDELITY_REPORT.md) (canonical
cross-pass record), [HISTORIC_NUCLEAR_FIDELITY.md](HISTORIC_NUCLEAR_FIDELITY.md)
(per-pass narrative), and this document (gate-by-gate detail) in lockstep.
The integration test source and the doc must move together.
