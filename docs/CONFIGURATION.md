# Configuration Reference

FSRM simulations are driven by a single `.config` text file consumed by
the `fsrm` executable. This document is the per-section, per-key
reference. For at-a-glance fidelity-tier selection, see
[FIDELITY_LADDER_GUIDE.md](FIDELITY_LADDER_GUIDE.md). For an end-to-end
worked example, see [USER_GUIDE.md](USER_GUIDE.md).

## File format

```ini
# Comment
[SECTION_NAME]
key = value
string_key = string with no quotes
list_key = value1, value2, value3
scientific = 1.5e6
```

The parser is in `src/core/ConfigReader.cpp`. Section names are
uppercase, keys are lowercase. Unknown keys are silently ignored;
unknown values for enum keys warn-and-fall-back to the default with
a one-time stderr message on rank 0.

## Sections

The recognized sections are listed below. Unset sections take their
documented defaults; the source of truth for every default is
`src/core/Simulator.cpp::initializeFromConfigFile()`.

| Section | Purpose |
|---|---|
| `[SIMULATION]` | Top-level run settings, physics enable flags |
| `[GRID]` | Mesh definition (structured or Gmsh) |
| `[MATERIAL]` | Per-region elastic / poroelastic / thermal properties |
| `[MESH_REFINEMENT]` | Adaptive refinement around the source |
| `[EXPLOSION_SOURCE]` | Underground-explosion top-level parameters |
| `[NEAR_FIELD_SOURCE]` | Six-rung fidelity-ladder source physics |
| `[SOURCE_DISTRIBUTION]` | How the moment tensor is distributed over cells |
| `[BOUNDARY_CONDITIONS]` | Per-face Dirichlet / traction |
| `[INITIAL_CONDITIONS]` | Per-field initial values |
| `[ABSORBING_BC]` | Clayton-Engquist absorbing on bounding-box faces |
| `[THERMAL]` | Thermal field configuration |
| `[VISCOELASTIC]` | Generalized-Maxwell viscoelastic memory variables |
| `[PLASTICITY]` | Drucker-Prager plasticity parameters |
| `[FAULT]` | Cohesive fault network |
| `[FRACTURE_PLANE]` | Pressurized-fracture formulation |
| `[HYDRAULIC_FRACTURE]` | Hydraulic-fracture coupling |
| `[INJECTION]` | Pressure-injection point source |
| `[SEISMOMETERS]` | Receiver locations and SAC station metadata |
| `[OUTPUT]` | What to write and where |
| `[WAVEFORM_VV]` | IRIS waveform-comparison V&V |
| `[NUCLEAR_TRIGGER]` | Coupling between explosion source and fault |

## `[SIMULATION]`

| Key | Type | Default | Description |
|---|---|---|---|
| `name` | string | "FSRM" | Run identifier |
| `start_time` | double | 0.0 | Simulation start, seconds |
| `end_time` | double | 1.0 | Simulation end, seconds |
| `dt` | double | 0.001 | Initial timestep, seconds |
| `enable_geomechanics` | bool | false | Solid mechanics on |
| `enable_elastodynamics` | bool | false | TSALPHA2 dynamic ODE |
| `enable_poroelasticity` | bool | false | Biot poroelasticity |
| `enable_thermal` | bool | false | Heat equation |
| `enable_faults` | bool | false | Cohesive fault network |
| `enable_viscoelastic` | bool | false | GMB attenuation |

## `[GRID]`

| Key | Type | Default | Description |
|---|---|---|---|
| `mesh_type` | enum | STRUCTURED | `STRUCTURED` or `GMSH` |
| `nx`, `ny`, `nz` | int | 10 | Cells per axis (structured) |
| `Lx`, `Ly`, `Lz` | double | 1.0 | Domain extents (structured), meters |
| `mesh_file` | string | -- | Path to MSH2 file (`mesh_type = GMSH`) |
| `cell_type` | enum | HEX | `HEX` or `TET` (faults require simplices) |

## `[MATERIAL]`

Three input modes are supported:

1. Single homogeneous material (omit per-region keys).
2. Depth-layered: `layer_count = N` plus `layer_<n>_z_min`, `layer_<n>_z_max`,
   and per-layer properties.
3. Per-cell binary velocity model: `velocity_model_path =
   path/to/file.bin` reads (Vp, Vs, rho) trilinearly interpolated to mesh
   centroids.

Common keys:

| Key | Type | Description |
|---|---|---|
| `density` | double | kg/m^3 |
| `youngs_modulus` | double | Pa |
| `poisson_ratio` | double | dimensionless |
| `permeability` | double | m^2 |
| `porosity` | double | 0..1 |
| `biot_coefficient` | double | dimensionless |
| `biot_modulus` | double | Pa |
| `thermal_conductivity` | double | W/(m K) |
| `specific_heat` | double | J/(kg K) |
| `thermal_expansion_coefficient` | double | 1/K |

## `[EXPLOSION_SOURCE]`

| Key | Type | Default | Description |
|---|---|---|---|
| `yield_kt` | double | -- | Equivalent TNT yield in kilotons |
| `location` | 3 doubles | -- | Source coordinates, meters |
| `medium_label` | enum | GRANITE | `GRANITE`, `TUFF`, `SALT`, `ALLUVIUM`, `SHALE` |
| `mode` | enum | COUPLED_ANALYTIC | `COUPLED_ANALYTIC` (Mueller-Murphy + dynamic) or `PROXY` (legacy scalar) |
| `explosion_solve_mode` | enum | COUPLED_ANALYTIC | Same enum as `mode`; legacy alias |

## `[NEAR_FIELD_SOURCE]`

The six-rung fidelity-ladder block. Defaults are pre-pass-12 backwards-
compatible: a config that does not set a key gets the previous-pass
behaviour byte-for-byte.

| Key | Type | Default | Values |
|---|---|---|---|
| `mode` | enum | KINEMATIC_RDP | `KINEMATIC_RDP` (closed-form) or `DYNAMIC_PLASTIC` (1D Lagrangian) |
| `solver_kind` | enum | RADIAL_LAGRANGIAN (in pass-7+ DYNAMIC_PLASTIC) | `CLOSED_FORM` or `RADIAL_LAGRANGIAN` |
| `radiation_phase` | enum | ZELDOVICH_RAIZER | `ZELDOVICH_RAIZER` (LOW), `MARSHAK_GREY` (MED), `MARSHAK_MULTIGROUP` (HIGHEST), `SN_TRANSPORT` (named-only scaffold) |
| `cavity_eos` | enum | TILLOTSON | `IDEAL_GAS` (LOW), `TILLOTSON` (MED), `TILLOTSON_TABULATED_PATCH` (HIGH), `TABULATED_FULL` (HIGHEST) |
| `opacity_model` | enum | POWER_LAW_ZR | `CONSTANT` (LOW), `POWER_LAW_ZR` (MED), `TABULATED_PATCHED` (HIGH), `TABULATED_FULL` (HIGHEST) |
| `operator_splitting` | enum | LIE | `LIE` (LOW/MED), `STRANG` (HIGH), `STRANG_MULTIGROUP` (HIGHEST), `LIE_MULTIGROUP` |
| `time_integrator` | enum | EXPLICIT_EULER | `EXPLICIT_EULER`, `TVD_RK2`, `RK3_SSP` |
| `time_integrator_diffusion` | enum | BACKWARD_EULER | `BACKWARD_EULER`, `CRANK_NICOLSON`, `BDF2` |
| `cavity_geometry` | enum | SPHERICAL | `SPHERICAL` (only working value); `THREE_DIMENSIONAL` is the pass-13 axis-1b scaffold and throws |
| `cavity_initialization` | enum | PHYSICS_BASED | `PHYSICS_BASED` (Newton energy-partition solve) or `MANUAL` |
| `initial_cavity_radius_m` | double | -1 | If `MANUAL`, set this; -1 lets the solver choose |
| `radiation_transition_time_s` | double | -1 | Override Z-R transition time; -1 lets solver choose |
| `tillotson_parameter_set` | string | medium-derived | Override Tillotson parameter set |
| `tillotson_extrapolation_warning_threshold_pa` | double | 5e10 | Warn when Tillotson extrapolated above this |

### Multigroup (only used when `radiation_phase = MARSHAK_MULTIGROUP`)

| Key | Type | Default | Description |
|---|---|---|---|
| `radiation_n_groups` | int | 16 | Number of frequency groups (1..256) |
| `radiation_freq_min_hz` | double | 1e14 | Lowest group frequency |
| `radiation_freq_max_hz` | double | 1e18 | Highest group frequency |
| `radiation_simpson_points` | int | 17 | Per-group Planck-integral Simpson points |
| `output_per_group_radiation` | bool | false | Emit per-group radiation flux to HDF5 |

### Tabulated EOS / opacity (HIGH and HIGHEST tiers)

| Key | Type | Default | Description |
|---|---|---|---|
| `tabulated_eos_table_path` | string | -- | HDF5 path for cavity_eos = TABULATED_* |
| `tabulated_opacity_rosseland_path` | string | -- | HDF5 path for opacity Rosseland |
| `tabulated_opacity_planck_path` | string | -- | HDF5 path for opacity Planck |
| `tabulated_eos_blend_lower_pa` | double | 5e10 | Pressure patch lower bound |
| `tabulated_eos_blend_upper_pa` | double | 6e10 | Pressure patch upper bound |
| `tabulated_opacity_blend_lower_k` | double | 1.0e5 | Temperature patch lower bound |
| `tabulated_opacity_blend_upper_k` | double | 1.26e5 | Temperature patch upper bound |
| `kappa_constant_m2_per_kg` | double | 0 | Used when `opacity_model = CONSTANT` |

### Pass-11 sponge layer and outer-radius override

| Key | Type | Default | Description |
|---|---|---|---|
| `sponge_layer_enabled` | bool | false | Israeli-Orszag 1981 graded-damping sponge |
| `radial_outer_radius_m` | double | -1 | Override outer-domain radius; -1 lets solver choose |

### Newton coupling

| Key | Type | Default | Description |
|---|---|---|---|
| `radiation_max_newton_iter` | int | 10 | Max Newton iterations per timestep |
| `radiation_newton_tolerance` | double | 1e-6 | Newton residual tolerance |
| `radiation_handoff_debounce_steps` | int | 3 | Hydro-to-radiation handoff debounce |

### Diagnostic flag

| Key | Type | Default | Description |
|---|---|---|---|
| `operator_splitting_convergence_diagnostic` | bool | false | Per-substep order-check instrumentation |

## `[SOURCE_DISTRIBUTION]`

Pass-4 grammar that controls how the moment tensor for an
`[EXPLOSION_SOURCE]` is injected into the FEM residual. The default
preserves the pre-pass-4 single-cell injection.

| Key | Type | Default | Description |
|---|---|---|---|
| `mode` | enum | SINGLE_CELL | `SINGLE_CELL`, `GAUSSIAN`, `UNIFORM_SPHERE` |
| `support_radius_factor` | double | 1.0 | Multiplier on cavity radius |
| `gaussian_sigma_factor` | double | 0.5 | Sigma factor for `GAUSSIAN` mode |
| `min_cells` | int | 1 | Fall back to `SINGLE_CELL` if support has fewer cells |

The integrated moment density equals `M0` in all modes; per-cell weights
are normalized via `MPI_Allreduce` to preserve the low-frequency RDP
plateau across distribution choices.

## `[BOUNDARY_CONDITIONS]`

Per-face Dirichlet / traction. Faces are named by axis half-space:
`X_MIN`, `X_MAX`, `Y_MIN`, `Y_MAX`, `Z_MIN`, `Z_MAX`.

```ini
[BOUNDARY_CONDITIONS]
face = X_MIN
type = DIRICHLET
displacement = 0.0, 0.0, 0.0

face = Z_MAX
type = TRACTION
traction = 0.0, 0.0, -1.0e6   # downward 1 MPa
```

## `[ABSORBING_BC]`

| Key | Type | Default | Description |
|---|---|---|---|
| `enabled` | bool | false | Apply Clayton-Engquist on all 6 faces |

> 99 % energy absorption verified by `Physics.AbsorbingBC`.

## `[THERMAL]`

| Key | Type | Default | Description |
|---|---|---|---|
| `enabled` | bool | false | Heat equation on |
| `initial_temperature` | double | 293.0 | K |
| `reference_temperature` | double | 293.0 | K (THM coupling reference) |

## `[VISCOELASTIC]`

| Key | Type | Default | Description |
|---|---|---|---|
| `enabled` | bool | false | GMB attenuation on |
| `n_mechanisms` | int | 3 | Number of relaxation mechanisms |
| `tau_<n>` | double | -- | Relaxation time, seconds |
| `delta_mu_<n>` | double | -- | Shear modulus contribution |
| `delta_kappa_<n>` | double | -- | Bulk modulus contribution |

## `[PLASTICITY]`

| Key | Type | Default | Description |
|---|---|---|---|
| `enabled` | bool | false | Drucker-Prager on |
| `cohesion` | double | -- | Pa |
| `friction_angle` | double | -- | radians |
| `dilation_angle` | double | -- | radians |
| `hardening_modulus` | double | 0.0 | Pa |

## `[FAULT]`

| Key | Type | Default | Description |
|---|---|---|---|
| `enabled` | bool | false | Insert cohesive cells along fault |
| `mode` | enum | LOCKED | `LOCKED`, `PRESCRIBED_SLIP`, `SLIP_WEAKENING`, `RATE_STATE` |
| `friction_model` | enum | CONSTANT | `CONSTANT`, `SLIP_WEAKENING` |
| `mu_s` | double | 0.6 | Static friction (slip-weakening) |
| `mu_d` | double | 0.5 | Dynamic friction |
| `D_c` | double | 0.4 | Critical slip distance, meters |

## `[SEISMOMETERS]`

| Key | Type | Default | Description |
|---|---|---|---|
| `station_count` | int | 0 | Number of stations |
| `station_<n>` | 3 doubles | -- | Station coordinates, meters |
| `station_<n>_name` | string | "S<n>" | Station identifier in SAC headers |
| `sac_sampling_rate_hz` | double | 100.0 | Output sampling rate |

## `[OUTPUT]`

| Key | Type | Default | Description |
|---|---|---|---|
| `output_dir` | string | output | Output directory |
| `hdf5_enabled` | bool | true | Full wavefield HDF5 |
| `hdf5_write_interval` | double | 0.1 | Seconds between snapshots |
| `vtk_enabled` | bool | false | VTK output |
| `sac_enabled` | bool | true | Per-station SAC seismograms |
| `write_near_field_history` | bool | true | M(t), Mdot(t), R_cav(t) CSV |
| `write_near_field_profile_h5` | bool | true | Per-snapshot radial state HDF5 |
| `write_derived_fields` | bool | false | Stress, strain, CFS at each snapshot |

## `[WAVEFORM_VV]`

| Key | Type | Default | Description |
|---|---|---|---|
| `enabled` | bool | false | Run IRIS waveform-comparison gate |
| `event_id` | string | -- | Cache subdirectory name |
| `cache_dir` | string | tools/waveform_vv/cache/<event_id> | SAC cache location |
| `band_lo_hz` | double | 0.5 | Source-physics-dominated low-frequency cutoff |
| `band_hi_hz` | double | 5.0 | High-frequency cutoff |
| `correlation_threshold` | double | 0.7 | Min cross-correlation gate |
| `peak_amplitude_envelope_factor` | double | 2.0 | Envelope factor for peak-amplitude gate |

See [WAVEFORM_VV.md](WAVEFORM_VV.md) for the cached-data format.

## See also

- [USER_GUIDE.md](USER_GUIDE.md) for an end-to-end worked example.
- [FIDELITY_LADDER_GUIDE.md](FIDELITY_LADDER_GUIDE.md) for tier selection.
- [EXPLOSION_IMPACT_PHYSICS.md](EXPLOSION_IMPACT_PHYSICS.md) for the
  per-pass physics description.
- `config/complete_template.config` and `config/default.config` are
  schema-anchor templates checked into the repository.
