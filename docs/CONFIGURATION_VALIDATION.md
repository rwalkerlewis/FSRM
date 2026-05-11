# Strict Configuration Validation

Pass-12 followup 2 (V&V hardening) added a strict configuration
validator to FSRM. Every configuration loaded by the production
`Simulator::initializeFromConfigFile` entry point is now checked
against an in-tree schema before any physics module sees it.
This document describes what the validator catches, how to
extend the schema, and the opt-out mechanisms.

## Why this exists

Pass-12 (PR #127, #128) found two pass-12 bug classes that the
existing test suite did not catch:

1. Configs that parsed cleanly but did not couple. An
   `[EXPLOSION_SOURCE]` block without a `[BOUNDARY_CONDITIONS]`
   block silently fell back to a Dirichlet "compression top,
   fixed bottom" default that suppressed the explosion source.
   Inline `[SEISMOMETERS] station_<n> = ...` keys were silently
   accepted by the parser but never consumed; the simulation
   ran with zero seismometers. An `[OUTPUT]` block was
   accepted at the top level but had no effect because the
   simulator reads output settings from `[SIMULATION]`.
2. Examples that did not run end-to-end under MPI. The
   default `-pc_type lu` (KLU) preconditioner fails on MPIAIJ
   matrices in PETSc 3.25 builds without MUMPS / SuperLU_DIST.
   Until PR #128 patched `scripts/run_with_mpi.sh`, every
   example exited non-zero at `MPI_RANKS > 1`.

The strict validator addresses class 1. The new
`examples_runtime` CTest gate addresses class 2; see
`docs/USER_GUIDE.md` "Parallel KSP / PC".

## What strict mode catches

The validator runs after `ConfigReader::loadFile` returns
successfully, before the parsed data is handed to physics
modules. It rejects four classes of issue.

### 1. Unknown top-level sections

Every `[SECTION]` header must appear in the registry kept in
`src/io/ConfigValidator.cpp`. Sections that follow the
`PREFIX_<n>` convention (such as `[LAYER_3]` or
`[SEISMOMETER_2]`) are matched against the dynamic-prefix
table.

Failing example:

```
[NOT_A_SECTION]
foo = bar
```

Error:

```
Unknown section [NOT_A_SECTION]. If this is a new section,
add it to the registry in src/io/ConfigValidator.cpp
(knownSections or dynamicSections). [...]
```

### 2. Unknown keys within a known section

Each section in the registry enumerates its allowed keys.
Anything outside that set is rejected with a Levenshtein-based
"did you mean" hint when a close match exists.

Failing example:

```
[SIMULATION]
end_tim = 5.0      # typo of end_time
```

Error:

```
Unknown key at config.config: [SIMULATION].end_tim. Did you
mean `end_time`?
```

### 3. Deprecated key forms

A small curated list of forms that previously existed in the
schema but now have replacements. Each rejection points at the
canonical form.

| Section | Key | Replacement |
| --- | --- | --- |
| `[SEISMOMETERS]` | `station_<n>` | Add a `[SEISMOMETER_<n>]` block per station with `sta = <name>` and `location_xyz = x,y,z`. |
| `[SEISMOMETERS]` | `station_count` | Inferred from the number of `[SEISMOMETER_<n>]` blocks. Remove. |
| `[SEISMOMETERS]` | `sac_sampling_rate_hz` | Renamed to `default_sample_rate_hz`. |
| `[SEISMOMETERS]` | `hdf5_enabled` / `sac_enabled` | Use `formats = SAC,HDF5` (comma-separated). |
| `[EXPLOSION]` | (any) | Rename the block to `[EXPLOSION_SOURCE]`; the `[EXPLOSION]` block is silently ignored by the parser. |

### 4. Required-section omissions

Some sections require a partner section so the config does not
fall back to a silently-incorrect default. Today there is one
such rule.

`[EXPLOSION_SOURCE]` requires `[BOUNDARY_CONDITIONS]`. Without
the explicit BC block, the simulator falls back to a Dirichlet
compression-top / fixed-bottom default that suppresses the
explosion source coupling. The validator emits an error
naming both sections and gives an example BC block in the
message.

## Opt-out mechanisms

Strict validation is the default. Two opt-out mechanisms exist
for backward compatibility; neither runs silently.

### Per-file opt-out

Add a `[META]` section to the config file:

```
[META]
strict_validation = false
```

The validator then skips the strict checks and emits a stderr
warning so the opt-out is visible in the run output:

```
[ConfigValidator][warn] config.config: [META] strict_validation =
false. Strict config validation has been disabled for this file.
This is intended only for short-lived debugging; remove the
opt-out before committing.
```

This is the right knob for short-lived debugging or for
`config/complete_template.config`-style schema-anchor files
that intentionally enumerate every option for documentation.

### Global env-var bypass

Set the environment variable `FSRM_DISABLE_STRICT_VALIDATION=1`
before invoking `fsrm`. The validator returns success without
running any checks. Use this only as an emergency override for
legacy callers that produce dirty configs and cannot be
updated immediately. There is no warning emitted in this mode.

## How to extend the schema

The schema lives in one place: `src/io/ConfigValidator.cpp`,
inside the `knownSections()` and `dynamicSections()` static
maps. Adding a new section or key is two edits and a rebuild.

### Adding a new section

1. Add an entry to `knownSections()` keyed by the new section
   name. The value is the `KeySet` of every allowed key in
   that section.
2. If the section is conceptually one of a family (e.g. a new
   `LAYER_<n>` style), add a prefix entry to
   `dynamicSections()` instead. Set `numeric_suffix_only` to
   `true` if the suffix must be an integer.

### Adding a new key

Add the key string to the appropriate section's `KeySet` in
`knownSections()` (or to the prefix's `KeySet` in
`dynamicSections()`). The validator picks it up automatically.

### Adding a new deprecation

Add an entry to `deprecatedForms()` with the offending
section, the offending key (supports `*` suffix wildcard), and
a clear message that names the replacement form.

## How to extend the required-section rules

`sectionsRequiringBCs()` lists sections that force
`[BOUNDARY_CONDITIONS]` to be present. The pattern generalises
trivially: add a new `requires<X>` predicate alongside it for
any other "if A is present then B must be too" rule.

## Coverage

`tests/unit/io/test_config_strict_validation.cpp` ships ten
unit tests, one per rejection class plus two opt-out tests
plus the `AllShippedExampleConfigsValidate` gate that walks
every `*.config` under `examples/` and asserts it validates.
The gate runs in under one second.

## Known broken examples (parallel KSP convergence)

Pass-12 followup 2 added the `examples_runtime` CTest gate
that runs every example end-to-end at `MPI_RANKS=4`. Fifteen
shipped examples failed that gate at the time. Pass-12
followup 3 (the parallel KSP/PC fix, see
[PARALLEL_KSP.md](PARALLEL_KSP.md)) unblocked twelve of them:
ten "SNES diverges step 0" cases shared one root cause -- the
parallel option string `scripts/run_with_mpi.sh` injected
overrode `-pc_type` to `bjacobi` but left the KSP type at the
`KSPPREONLY` that `Simulator::setupSolvers` sets, so each
Newton step was a single block-Jacobi sweep -- and the fix
was to add `-ksp_type gmres` to that string. The remaining
two (`06_gmsh_multimaterial`, `17_velocity_model`) were
mislabelled: 06 had a stale relative `mesh_file` path after
pass-12 relocated its config into the example directory, and
17 never shipped the `velocity_model.bin` it references; both
got example-local fixes.

The three cohesive-cells fault examples (`04_locked_fault`,
`08_time_dependent_slip`, `16_scec_tpv5`) stay in
`EXAMPLE_SMOKE_MPI4_KNOWN_BROKEN`. They are blocked on the
PETSc 3.25 saddle-point solver tracked in
`docs/SOLVER_STATE.md`, which is a separate axis from the
launcher fix.

The list lives in `tests/CMakeLists.txt` as
`EXAMPLE_SMOKE_MPI4_KNOWN_BROKEN`. Removing an entry from
that list flips the gate for that example from inverted-pass
to genuine-pass, so a fix that lands separately is detected
the moment the test starts succeeding.

To inspect the diagnosis:

```bash
# The committed before/after PETSc-verbose run logs:
ls tests/diagnostics/parallel_ksp_baseline/

# Regenerate one (post-fix it shows convergence):
ctest -L diagnostic -R 09_gasbuggy_1967 --output-on-failure
```

The smoke wrapper supports a `config_ci.config` override per
example: if `examples/<N>_<event>/config_ci.config` exists
the wrapper sets `FSRM_CONFIG_OVERRIDE` so a future
`run.sh` enhancement can pick it up. Today only
`examples/16_scec_tpv5/config_ci.config` ships such a
variant.

## Cross-references

* `tests/unit/io/test_config_strict_validation.cpp`: the unit
  test suite for the validator.
* `tests/integration/run_example_smoke.sh`: the smoke wrapper
  used by the `examples_runtime` CTest gate; sets
  `FSRM_FINAL_TIME_OVERRIDE` so each example completes in
  under 30 seconds.
* `docs/CONFIGURATION.md`: full reference for every section
  and key that ships with FSRM, including the
  `FSRM_FINAL_TIME_OVERRIDE` env var added in this pass.
* `docs/HISTORIC_NUCLEAR_FIDELITY.md` 4l: the closeout entry
  for pass-12 followup 2.
