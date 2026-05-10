# FSRM User Guide

End-to-end manual for running FSRM simulations and producing seismograms.
For a five-minute build-and-run, see [QUICK_START.md](QUICK_START.md). For
the full config-key reference, see [CONFIGURATION.md](CONFIGURATION.md).

## Contents

1. [What FSRM does](#what-fsrm-does)
2. [Running a simulation](#running-a-simulation)
3. [Running in parallel](#running-in-parallel)
4. [Config-file structure](#config-file-structure)
5. [The source-physics fidelity ladder](#the-source-physics-fidelity-ladder)
6. [Boundary conditions and meshes](#boundary-conditions-and-meshes)
7. [Output and visualization](#output-and-visualization)
8. [IRIS waveform V&V](#iris-waveform-vv)
9. [Troubleshooting](#troubleshooting)

## What FSRM does

FSRM solves coupled multiphysics PDEs on unstructured DMPlex meshes via
PETSc PetscDS pointwise callbacks. Verified physics (each with quantitative
integration tests through TSSolve):

- Linear elasticity (quasi-static and dynamic, TSALPHA2 generalized-α)
- Biot poroelasticity
- Drucker-Prager elastoplasticity
- Generalized-Maxwell viscoelastic attenuation
- Heat equation and full thermo-hydro-mechanical coupling
- Cohesive-cell faults (locked, prescribed-slip, slip-weakening)
- Underground-explosion source physics with a six-rung fidelity ladder
- Clayton-Engquist absorbing boundaries
- Single-phase Darcy flow

Twenty-seven historic underground nuclear tests ship as integration
fixtures with pinned configs (see [HISTORIC_NUCLEAR_FIDELITY.md](HISTORIC_NUCLEAR_FIDELITY.md)).
Pass-12 housekeeping adds DPRK 2017 and other modern events.

## Running a simulation

Each `examples/N_<event>/` is self-contained. Run any one with:

```bash
cd examples/20_salmon_1964
./run.sh
```

The runner builds output in `examples/20_salmon_1964/output/`. It calls
`scripts/run_with_mpi.sh` which wraps `mpirun` with the right
ABI-specific flags (OpenMPI vs MPICH detected at runtime). Override
the rank count with the `MPI_RANKS` env var (default 4 for examples,
8 for showcase events).

You can also call the executable directly:

```bash
mpirun -n 4 ./build/fsrm -c examples/20_salmon_1964/config.config
```

## Running in parallel

The FEM far-field is fully MPI-parallel via PETSc DMPlex. The 1-D radial
source-ball solver (`src/domain/explosion/RadialLagrangian.cpp`) is
deliberately serial: at production resolution (~800 cells x 16 frequency
groups for multigroup) the radial solve is small tridiagonal kernels per
timestep and MPI overhead would exceed any saved time. Pass-13's axis-1b
3-D source ball will introduce source-side MPI parallelism via PETSc.

Recommended rank counts:

| Mesh resolution | Recommended `MPI_RANKS` |
|---|---|
| ~10k cells (tutorial / unit) | 1-2 |
| ~100k cells (CI examples) | 4 |
| ~1M cells (showcase / production) | 8-16 |

`Integration.MPI.SalmonSerialVsParallelEquivalence` verifies that
serial and parallel produce numerically equivalent SAC output (peak
amplitude relative difference under 1 %, cross-correlation above
0.99 in the 0.5-5 Hz band). If a new event you author drifts these
numbers between rank counts, the FEM partitioning has lost
correctness; check rank-aware logic added since.

### Parallel KSP / PC

The PETSc 3.25 build that ships in the FSRM Docker image does not
include MUMPS or SuperLU_DIST, so PETSc's default `-pc_type lu` (KLU)
fails on `MPIAIJ` matrices when running with more than one rank.
`scripts/run_with_mpi.sh` automatically appends
`-pc_type bjacobi -sub_pc_type lu` when `MPI_RANKS > 1` so each rank
does its own block LU. Override or extend by setting
`FSRM_MPI_PETSC_OPTS` in the environment (set to empty string to
disable injection):

```bash
FSRM_MPI_PETSC_OPTS="-pc_type gamg -ksp_type cg" MPI_RANKS=8 ./run.sh
```

## Config-file structure

Configs use INI format. Each section name is uppercase; keys are
lowercase. Comments start with `#`.

```ini
[SECTION]
key = value                  # comment
scientific = 1.5e6           # scientific notation
```

The parser's recognized sections are `SIMULATION`, `GRID`, `MATERIAL`,
`MESH_REFINEMENT`, `EXPLOSION_SOURCE`, `NEAR_FIELD_SOURCE`,
`SOURCE_DISTRIBUTION`, `BOUNDARY_CONDITIONS`, `INITIAL_CONDITIONS`,
`ABSORBING_BC`, `THERMAL`, `VISCOELASTIC`, `PLASTICITY`, `FAULT`,
`FRACTURE_PLANE`, `HYDRAULIC_FRACTURE`, `INJECTION`, `SEISMOMETERS`,
`SEISMICITY`, `OUTPUT`, `NUCLEAR_TRIGGER`, `WAVEFORM_VV`. See
[CONFIGURATION.md](CONFIGURATION.md) for the full key reference.

A minimal explosion-seismogram config has the structure below; the
historic-nuclear configs add layered-material, source-distribution,
and per-event timing detail.

```ini
[SIMULATION]
end_time = 8.0                  # seconds
enable_elastodynamics = true
enable_geomechanics = true

[GRID]
mesh_type = STRUCTURED
nx = 24
ny = 24
nz = 24
Lx = 4000.0
Ly = 4000.0
Lz = 4000.0

[MATERIAL]
density = 2700.0
youngs_modulus = 60.0e9
poisson_ratio = 0.25

[ABSORBING_BC]
enabled = true                  # Clayton-Engquist on all six faces

[EXPLOSION_SOURCE]
yield_kt = 5.3
location = 2000.0, 2000.0, 1300.0
medium_label = SALT
mode = COUPLED_ANALYTIC

[NEAR_FIELD_SOURCE]
mode = DYNAMIC_PLASTIC
solver_kind = RADIAL_LAGRANGIAN
cavity_eos = TILLOTSON

[SEISMOMETERS]
station_count = 3
station_1 = 1000.0, 2000.0, 0.0
station_2 = 2000.0, 1000.0, 0.0
station_3 = 3000.0, 2000.0, 0.0

[OUTPUT]
sac_enabled = true
hdf5_enabled = true
output_dir = output
```

## The source-physics fidelity ladder

The `[NEAR_FIELD_SOURCE]` block exposes six fidelity-ladder sub-keys.
Each maps to a LOW / MED / HIGH / HIGHEST tier; see
[FIDELITY_LADDER_GUIDE.md](FIDELITY_LADDER_GUIDE.md) for the at-a-glance
selection table.

```ini
[NEAR_FIELD_SOURCE]
mode = DYNAMIC_PLASTIC                # KINEMATIC_RDP for closed-form path
solver_kind = RADIAL_LAGRANGIAN       # CLOSED_FORM for pass-5 byte-identical

# Fidelity ladders
radiation_phase = MARSHAK_GREY        # ZELDOVICH_RAIZER (LOW) | MARSHAK_GREY (MED) | MARSHAK_MULTIGROUP (HIGHEST)
cavity_eos = TILLOTSON                # IDEAL_GAS | TILLOTSON | TILLOTSON_TABULATED_PATCH | TABULATED_FULL
opacity_model = POWER_LAW_ZR          # CONSTANT | POWER_LAW_ZR | TABULATED_PATCHED | TABULATED_FULL
operator_splitting = LIE              # LIE | STRANG | STRANG_MULTIGROUP
time_integrator = EXPLICIT_EULER      # EXPLICIT_EULER | TVD_RK2 | RK3_SSP
time_integrator_diffusion = BDF2      # BACKWARD_EULER | CRANK_NICOLSON | BDF2

# Multigroup parameters (only used when radiation_phase = MARSHAK_MULTIGROUP)
radiation_n_groups = 16
radiation_freq_min_hz = 1.0e14
radiation_freq_max_hz = 1.0e18

# Tabulated data paths (only used by TABULATED_* variants)
tabulated_eos_table_path = tools/tabulated_data/tables/eos/granite_aneos.h5
tabulated_opacity_rosseland_path = tools/tabulated_data/tables/opacity/granite_rosseland.h5
tabulated_opacity_planck_path = tools/tabulated_data/tables/opacity/granite_planck.h5
tabulated_eos_blend_lower_pa = 5.0e10
tabulated_eos_blend_upper_pa = 6.0e10
tabulated_opacity_blend_lower_k = 1.0e5
tabulated_opacity_blend_upper_k = 1.26e5

# Pass-11 sponge layer (Israeli-Orszag 1981) and outer-radius override
sponge_layer_enabled = false
radial_outer_radius_m = -1.0          # -1 lets the solver choose

# Cavity geometry (axis-1b scaffold; pass-13 implementation)
cavity_geometry = SPHERICAL           # THREE_DIMENSIONAL throws until pass-13
```

Defaults reproduce the pre-pass-12 behaviour. Selecting an unknown
value yields a warn-and-fall-back to the safe LOW tier on rank 0.

## Boundary conditions and meshes

Three orthogonal mechanisms ship:

- **Clayton-Engquist absorbing**: `[ABSORBING_BC] enabled = true` adds
  first-order absorbing tractions on all six bounding-box faces.
  `Physics.AbsorbingBC` verifies > 99 % energy absorption.
- **Per-face Dirichlet / traction**: `[BOUNDARY_CONDITIONS]` block
  with `face = X_MIN | X_MAX | Y_MIN | ...` and explicit per-component
  vectors. See `examples/07_traction_bc/`.
- **Gmsh per-region material assignment**: `[GRID] mesh_type = GMSH`
  with `mesh_file = path/to/mesh.msh` (MSH2 format). Physical names in
  the mesh map to material labels in the config. See
  `examples/06_gmsh_multimaterial/`.

For sub-cell-resolution sources, `[MESH_REFINEMENT]` adaptively
refines around the explosion location. For per-cell velocity-model
material assignment from binary `Vp/Vs/rho` grids, set
`[MATERIAL] velocity_model_path = path/to/velocity_model.bin`.

## Output and visualization

`[OUTPUT]` controls what the simulator writes. Defaults are fine for
the historic-nuclear examples.

```ini
[OUTPUT]
output_dir = output
hdf5_enabled = true                   # full wavefield, large
hdf5_write_interval = 0.1             # seconds
sac_enabled = true                    # SAC seismograms at SEISMOMETERS
write_near_field_history = true       # 6-component M(t), Mdot(t), R_cav(t)
write_near_field_profile_h5 = true    # per-snapshot radial state
```

After a run, the `output/` directory typically contains:

- `seismograms/<station>.{r,t,z}.sac`: rotated to radial / transverse /
  vertical at each `[SEISMOMETERS]` station.
- `near_field_history.csv`: 14 columns -- time, six components of M(t),
  six of Mdot(t), and the cavity radius.
- `near_field_profile.h5` + `.xdmf`: per-snapshot radial state of the
  source ball (vp, vs, rho, T_matter, T_radiation, plastic strain). The
  XDMF wrapper opens directly in ParaView.
- `solution.h5` + `.xmf`: 3-D wavefield at every output cadence. Open
  the XMF in ParaView to visualize displacement / velocity / stress.

For the showcase events (Sedan 1962, Salmon 1964, Punggye-ri 2017,
Cannikin 1971, Sterling 1966), `figures/regenerate.sh` produces a
six-figure presentation pack from the output directory using shared
matplotlib styling in `tools/figures/figure_style.py`.

## IRIS waveform V&V

The `[WAVEFORM_VV]` block opts the run into IRIS waveform-comparison
gates. Cached station data lives under `tools/waveform_vv/cache/`.
Refresh with `tools/waveform_vv/refresh.py` (requires ObsPy).

```ini
[WAVEFORM_VV]
enabled = true
event_id = SALMON_1964
cache_dir = tools/waveform_vv/cache/salmon_1964
band_lo_hz = 0.5
band_hi_hz = 5.0
correlation_threshold = 0.7
peak_amplitude_envelope_factor = 2.0
```

The integration test family `Physics.WaveformVV.*` runs the comparison
gates. The `iris_validation` ctest label gates the suite as a whole. See
[WAVEFORM_VV.md](WAVEFORM_VV.md) for the cached-data format and the
waveform-comparison metric definitions.

## Troubleshooting

### Run dies with "cavity_geometry=THREE_DIMENSIONAL not implemented"

This is the pass-13 axis-1b scaffold guard. Set `cavity_geometry =
SPHERICAL` (the default).

### `radiation_phase = SN_TRANSPORT` raises a runtime error

`SN_TRANSPORT` is a HIGHEST-tier scaffold; it is not implemented and is
not on the pass-13 roadmap. Use `MARSHAK_MULTIGROUP`.

### "TABULATED_FULL requires a table" warning, results look like power-law fallback

The path you set for `tabulated_eos_table_path` (or the opacity paths)
either does not exist or failed HDF5 parsing. Check
`tools/tabulated_data/README.md` for the expected layout. The solver
emits a one-time stderr warning when it falls back to the analytic.

### Convergence fails / NaN in radial solver

The CFL constant is set automatically. NaNs in the source ball usually
indicate a Tillotson extrapolation past the threshold; raise
`tillotson_extrapolation_warning_threshold_pa` if you genuinely need
to operate above it (default 5e10 Pa).

### `mpirun` errors when running as root in Docker

The wrapper `scripts/run_with_mpi.sh` (sourced by every example
`run.sh`) detects this and applies `--allow-run-as-root`
automatically. Bypass the wrapper only if you know the MPI ABI.

### Six fault tests fail in `ctest`

These are documented in [SOLVER_STATE.md](SOLVER_STATE.md). They sit
behind a PETSc 3.25 BdResidual limitation that FSRM cannot fix from
the application side. They do not block any historic-nuclear run.

## References

- PETSc 3.25.0: https://petsc.org/
- PyLith verified architecture pins: [PYLITH_REFERENCE.md](PYLITH_REFERENCE.md)
- Source physics primary literature: see [EXPLOSION_IMPACT_PHYSICS.md](EXPLOSION_IMPACT_PHYSICS.md) "References"
- Numerical methods: see [NUMERICAL_METHODS.md](NUMERICAL_METHODS.md)
