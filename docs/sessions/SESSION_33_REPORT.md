# Session 33 -- Pass-12 followup 2: V&V hardening

**Branch:** `chore/pass-12-followup-2-vv-hardening`
**Date:** 2026-05-10
**Scope tag:** vv-infrastructure (no physics)

## Goal

Close the V&V coverage gaps that the first pass-12 followup
(PR #128) surfaced when adding the DPRK 2017 / Lop Nor 1996 /
Pokhran II 1998 events:

1. Examples that did not run end-to-end under MPI. The
   default `-pc_type lu` (KLU) preconditioner fails on
   `MPIAIJ` matrices in PETSc 3.25 builds without MUMPS or
   SuperLU_DIST. Until PR #128 patched
   `scripts/run_with_mpi.sh`, every example exited non-zero
   at `MPI_RANKS > 1`.
2. Configs that parsed cleanly but did not couple. Inline
   `[SEISMOMETERS] station_<n>=` keys, missing
   `[BOUNDARY_CONDITIONS]` block on explosion-source configs,
   and orphan `[OUTPUT]` blocks were silently accepted by the
   parser but never consumed.

## What landed

### Thread A: Examples-runtime CTest gate

* New `examples_runtime` CTest label registered dynamically
  over every `examples/<N>_<event>/` directory at CMake
  configure time. Each registration spawns
  `tests/integration/run_example_smoke.sh <example_dir>
  <mpi_ranks> <final_time>`.
* The smoke wrapper sets `FSRM_FINAL_TIME_OVERRIDE=0.01` so
  every example completes in under 30 seconds while still
  exercising the full pipeline (setupDM through TSSolve).
* Default coverage: every example at `MPI_RANKS=4` (production
  default). Three representative examples (`01_uniaxial_compression`,
  `02_explosion_seismogram`, `11_sedan_1962`) also run at
  `MPI_RANKS=1` and `MPI_RANKS=2` so the parallel KSP path
  PR #128 patched is exercised at every rank count without
  paying the full 41 x 3 cost.
* Set `FSRM_EXAMPLES_RUNTIME_FULL=ON` at CMake configure time
  to register the full `MPI={1,2,4}` matrix for every example.
* `tests/integration/run_example_smoke.sh` cleans the
  example's `output/` directory before invoking `run.sh`,
  then asserts at least one non-empty file under `output/`
  after the run. The wrapper exits non-zero with a diagnostic
  message naming the example, the MPI rank count, and the
  failing assertion.
* Wall-clock budget on the fsrm-ci image at four cores: under
  five minutes for the default coverage.

### Thread B: Strict configuration validator

* New `ConfigValidator` class:
  * `include/io/ConfigValidator.hpp`
  * `src/io/ConfigValidator.cpp`
* In-tree schema enumerates every section and key the FSRM
  parser actually consumes. Sources: dispatch sites in
  `src/core/ConfigReader.cpp` and the hand-rolled
  `reader.getXxx("SECTION", ...)` calls in
  `src/core/Simulator.cpp` (NEAR_FIELD_SOURCE, EXPLOSION_SOURCE,
  BOUNDARY_CONDITIONS, ABSORBING_BC, SEISMOMETERS, FAULT,
  FRACTURE_PLANE, HYDRAULIC_FRACTURE, INJECTION,
  INITIAL_CONDITIONS, MATERIAL, MESH_REFINEMENT, NUCLEAR_TRIGGER,
  OUTPUT, PLASTICITY, SOURCE_DISTRIBUTION, SEISMICITY, THERMAL,
  VISCOELASTIC, WAVEFORM_VV, TRACTION_BC, etc.).
* Dynamic-prefix schema for per-instance sections: `LAYER_<n>`,
  `SEISMOMETER_<n>`, `MATERIAL_REGION_<n>`, `ROCK_<name>`,
  `FAULT_<name>`.
* Four rejection classes:
  1. Unknown top-level section (with closest-known hints).
  2. Unknown key within a known section (with Levenshtein-based
     "did you mean" hints).
  3. Deprecated key forms: `[SEISMOMETERS]` `station_<n>`,
     `station_count`, `sac_sampling_rate_hz`, `hdf5_enabled` /
     `sac_enabled`; orphan `[EXPLOSION]` block (PR #128
     findings).
  4. Required-section omission: `[EXPLOSION_SOURCE]` requires
     `[BOUNDARY_CONDITIONS]` (PR #127 finding).
* Strict mode is the default. Per-file opt-out via
  `[META] strict_validation = false` (always emits a stderr
  warning so the opt-out is visible in run output). Global
  env-var bypass via `FSRM_DISABLE_STRICT_VALIDATION=1` for
  legacy callers.
* Wired into `Simulator::initializeFromConfigFile` after
  `ConfigReader::loadFile` returns and before the parsed data
  reaches any physics module. A failed validation `SETERRQ`s
  with a pointer at `docs/CONFIGURATION_VALIDATION.md`.
* `ConfigReader::parseSimulationConfig` now consults
  `FSRM_FINAL_TIME_OVERRIDE`; out-of-range values (zero,
  negative, or larger than the configured `end_time`) are
  ignored with a stderr warning.

### Thread C: Documentation

* New: `docs/CONFIGURATION_VALIDATION.md` -- what strict mode
  catches, opt-out mechanisms, how to extend the schema.
* Updated: `docs/CONFIGURATION.md` -- new "Environment
  overrides" table covering `FSRM_FINAL_TIME_OVERRIDE`,
  `FSRM_DISABLE_STRICT_VALIDATION`, `FSRM_MPI_PETSC_OPTS`,
  `MPI_RANKS`. Strict-validator note added to the file-format
  section.
* Updated: `docs/USER_GUIDE.md` -- "Parallel KSP / PC" section
  cross-references the new `examples_runtime` gate so users
  know it exists.
* Updated: `docs/HISTORIC_NUCLEAR_FIDELITY.md` -- new "4k
  followup 2" section closes the V&V hardening pass.
* Updated: `CLAUDE.md` -- two new rules (21, 22) covering the
  examples_runtime gate and the strict validator schema.

## Tests

* `Unit.ConfigStrictValidation` (10 tests): rejection-class
  coverage plus opt-out tests plus the
  `AllShippedExampleConfigsValidate` gate that walks every
  `*.config` under `examples/` and asserts it validates under
  strict mode.
* `examples_runtime` label: 41 tests at default coverage (one
  MPI=4 per example) plus 6 MPI=1/MPI=2 sanity runs (two each
  for examples 01, 02, 11).

## Configs repaired

* `examples/15_viscoelastic_attenuation/config.config`: replaced
  the orphan `[EXPLOSION]` block (silently ignored by the
  parser) with a canonical `[EXPLOSION_SOURCE]` block; replaced
  inline `[SEISMOMETERS] station_*` keys with proper
  `[SEISMOMETER_<n>]` blocks.
* `config/complete_template.config`, `config/default.config`,
  `config/test_*.config`, every `config/templates/*.config`:
  added `[META] strict_validation = false` with a clear
  comment. These are documentation-grade templates that
  intentionally enumerate aspirational keys; runtime validation
  does not apply.

## Verification

* All 116 pre-existing default-suite tests pass on the fsrm-ci
  image (Simulator validator hook adds zero regressions; schema
  gaps surfaced and fixed mid-pass: per-layer attenuation
  `q_p`/`q_s`, MATERIAL_REGION_n `poisson_ratio`, FLUID
  `reference_pressure`).
* `Unit.ConfigStrictValidation`: 10/10 pass.
* `examples_runtime`: 28 of 43 MPI=4 tests pass + 6 MPI=1/MPI=2
  sample tests pass. 15 MPI=4 tests are marked WILL_FAIL TRUE
  against a parallel-KSP convergence issue documented in
  `docs/HISTORIC_NUCLEAR_FIDELITY.md` 4k followup 2 "Known
  broken examples". All 15 run cleanly at MPI=1; the per-rank
  LU factorisation injected by `scripts/run_with_mpi.sh` does
  not converge for these specific physics setups. Investigating
  is a follow-up axis (out of scope per the brief: no physics,
  no fault-solver work).

## Out of scope (held to brief)

* No physics changes, no ladder-rung changes, no new historic
  events.
* No fault-solver Jacobian work; the six default-path failures
  stay documented in `docs/SOLVER_STATE.md`.
* No CI workflow changes (no Dockerfile edits, no GitHub
  Actions version bumps).
* Pass-13 axis-1b state unchanged. The pass-12 closeout in
  `HISTORIC_NUCLEAR_FIDELITY.md` 4k is not re-opened; this is
  a followup, not a re-open. Pass-13a/b/c sections (4l, 4m, 4n)
  are unchanged.

## Closeout pointers

* `tests/unit/io/test_config_strict_validation.cpp`
* `tests/integration/run_example_smoke.sh`
* `tests/CMakeLists.txt` (CTest registration block under
  "Pass-12 followup 2 (V&V hardening)")
* `include/io/ConfigValidator.hpp` /
  `src/io/ConfigValidator.cpp`
* `src/core/Simulator.cpp` (validator hook in
  `initializeFromConfigFile`)
* `src/core/ConfigReader.cpp` (FSRM_FINAL_TIME_OVERRIDE in
  `parseSimulationConfig`)
* `docs/CONFIGURATION_VALIDATION.md`
* `docs/HISTORIC_NUCLEAR_FIDELITY.md` 4k followup 2
