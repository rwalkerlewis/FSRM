# FSRM Nuclear Explosion Monitoring — Honest Baseline (Milestone 0)

This document is the inventory output of Milestone 0 of `docs/NEM_ROADMAP.md`. It is
inventory only. No files are deleted or modified by Milestone 0. The deletions happen
in the Milestone 3 cleanup pass, after at least the Sedan and Hardhat / Piledriver
events have been independently validated against observations.

This document is written from the `nem-roadmap` branch at the head commit recorded
in the session log. It will not be edited again unless that head commit changes.

## 1. Integration test status

Per the roadmap: "Run all integration tests on `main`. Record each test's pass/fail
status, runtime, and any error output."

### 1.1 Build attempt result

Docker is the supported build environment per `CLAUDE.md` rule 1 ("Build and test
in Docker. Always."). The build was attempted with:

```
docker build -f Dockerfile.ci -t fsrm-ci:local .
```

The build failed at the PETSc compile step (Dockerfile.ci stage 7/8) with:

```
#11 [7/8] RUN ./configure ...
#11 0.166 /bin/sh: 1: ./configure: not found
#11 ERROR: process "/bin/sh -c ./configure ... && make ..." did not complete
successfully: exit code: 127
```

Root cause: the preceding stage downloads `petsc-3.25.0.tar.gz` from
`web.cels.anl.gov/projects/petsc/download/release-snapshots/`, untars it, and then
runs `ls -d petsc-3.25* | head -1 | xargs -I{} mv {} petsc-3.25.0 || true`. The
`|| true` masks the failure if the tarball does not produce a `petsc-3.25*` directory
(for example because the URL now serves an HTML redirect or because the snapshot
filename has changed). The next stage `WORKDIR /opt/petsc-3.25.0` then succeeds
silently by creating an empty directory, and `./configure` is not found inside it.

Per roadmap ground rule 5 and the "Final notes for the agent" section: "If a tool,
library, or build step does not work, document the failure in the session log and
stop. Do not paper over it with a stub." Fixing the Dockerfile is not Milestone 0
work; Milestone 0 is inventory. The fix is the next concrete step recorded in the
session log of `docs/NEM_ROADMAP.md`.

### 1.2 Test inventory (static, by file)

Because no working build environment is available on this branch as of Milestone 0,
the integration tests could not be executed and per-test runtime / error output
could not be recorded. The static enumeration of the test suite is:

| Category               | Files | Source directory             |
| ---------------------- | ----: | ---------------------------- |
| Unit                   |    30 | `tests/unit/`                |
| Functional             |     4 | `tests/functional/`          |
| Physics validation     |    25 | `tests/physics_validation/`  |
| Integration            |    42 | `tests/integration/`         |
| Performance            |     4 | `tests/performance/`         |
| Experimental           |     1 | `tests/experimental/`        |
| **Total source files** | **106** | -                          |

The 42 integration test source files are:

```
test_coupled_hydrofrac           test_layered_elastostatics
test_coupled_physics             test_leakoff_coupling
test_derived_fields              test_locked_fault
test_dprk_2017_comparison        test_nearfield_coupled
test_dynamic_rupture_basic       test_nuclear_twin_gmsh
test_dynamic_rupture_solve       test_output_file
test_elastoplastic_sim           test_prescribed_slip
test_explosion_fault_reactivation test_pressurized_fracture_fem
test_explosion_seismogram        test_production_forecast
test_fault_absorbing_coexist     test_proppant_transport
test_fracture_flow               test_punggye_ri_layered
test_fracture_propagation        test_restart
test_full_simulation             test_single_phase_flow
test_gasbuggy_mesh               test_slip_weakening
test_gmsh_import                 test_slipping_fault_solve
test_gmsh_multimaterial          test_stress_shadowing
test_historic_nuclear            test_thermal_diffusion
test_induced_seismicity          test_thermal_expansion
test_injection_pressure          test_time_dependent_slip
test_injection_rupture_chain     test_traction_bc
                                 test_velocity_model_material
                                 test_viscoelastic_wave
```

A pass/fail/runtime table will be filled in once a working build environment is
restored. Until then, no claim can be made about which tests pass.

## 2. Unimplemented headers

Headers under `include/` with no matching `.cpp` file under `src/` (matched by base
filename, recursively). These are candidates for deletion in the Milestone 3 cleanup
pass. They are not deleted by Milestone 0.

| Header path                                                | Notes |
| ---------------------------------------------------------- | ----- |
| `include/numerics/AdaptiveMeshRefinement.hpp`              | AMR header; no matching `src/numerics/AdaptiveMeshRefinement.cpp`. `CLAUDE.md` notes the header is kept only because a test file depends on it. |
| `include/numerics/DiscontinuousGalerkin.hpp`               | DG header; no `.cpp`. Listed explicitly as "dead" in `docs/NEM_ROADMAP.md` and in `CLAUDE.md`. |
| `include/domain/atmospheric/AtmosphericInfrasound.hpp`     | No matching `.cpp`. Atmospheric infrasound is out of scope per the roadmap. |
| `include/domain/geomechanics/FaultModel.hpp`               | No matching `.cpp`. Fault physics is in `CohesiveFaultKernel.cpp`. |
| `include/core/FSRM.hpp`                                    | No matching `.cpp`. Likely an umbrella include header; verify before deletion. |
| `include/core/PhysicsModuleInterface.hpp`                  | No matching `.cpp`. May be a pure-virtual interface; verify before deletion. |
| `include/ml/FourierNeuralOperatorGPU.hpp`                  | GPU FNO header; no `.cpp`. CPU FNO exists in `src/ml/FourierNeuralOperator.cpp`. |
| `include/gpu/GPUCompute.hpp`                               | No matching `.cpp` under `src/gpu/`. The `src/gpu/` directory does not exist; only `include/gpu/`. |
| `include/gpu/GPUManager.hpp`                               | No matching `.cpp`. Same as above. |
| `include/physics/PhysicsKernel_GPU.hpp`                    | No matching `.cpp`. |
| `include/physics/PhysicsKernel_GPU_Extended.hpp`           | No matching `.cpp`. |

Total: 11 headers without matching implementation files.

Note that some of these may be pure interface headers (for example
`PhysicsModuleInterface.hpp`) where the lack of a `.cpp` is correct. Each entry
must be re-evaluated case-by-case at Milestone 3, not deleted blindly.

## 3. Stub implementations

`.cpp` files under `src/` whose function bodies are stubs (return zero or
`PETSC_SUCCESS`, empty body, or every parameter cast to `void` with no logic).

Heuristics like "many `(void)` casts" produce false positives because the PetscFE
pointwise callback files (`PetscFEElasticity.cpp`, `PetscFEPoroelasticity.cpp`,
`PetscFEFluidFlow.cpp`, etc.) legitimately use `(void)` to silence unused-parameter
warnings on the standard PETSc `f0_*`, `f1_*`, `g0_*`, `g3_*` signatures. Those
files have real callback math and are not stubs. To avoid false flagging of
production code, the list below is restricted to files manually inspected and
confirmed to consist almost entirely of empty bodies or trivial returns.

### 3.1 Confirmed stubs (`src/experimental/`)

All five files in `src/experimental/` are stubs. Each method body either returns
`false`, returns `PETSC_SUCCESS`, returns a hard-coded `Tensor`, or is empty.
Example, `src/experimental/NeuralAMR.cpp`, lines 12-38:

```
bool NeuralAMRController::shouldAdapt(int /* step */, double /* time */) {
    return false;
}

PetscErrorCode NeuralAMRController::adapt(DM&, Vec&, double) {
    return PETSC_SUCCESS;
}
```

| File                                          | LOC | Description |
| --------------------------------------------- | --: | ----------- |
| `src/experimental/NeuralAMR.cpp`              |  41 | Returns `false` / `PETSC_SUCCESS`; no neural model. |
| `src/experimental/NeuralInversion.cpp`        |  42 | Stub. |
| `src/experimental/NeuralLinearAlgebra.cpp`    |  53 | Stub; matches `docs/NEM_ROADMAP.md` "What is not real" list. |
| `src/experimental/NeuralTimestepping.cpp`     |  80 | Stub. |
| `src/experimental/MultiFidelityLearning.cpp`  |  44 | Stub. |

Total: 5 files, ~260 LOC.

### 3.2 Other suspects (need manual verification at Milestone 3)

The following files were flagged by the heuristic (high `(void)` cast count or
many trivial returns) but contain at least some real logic when spot-checked. They
are listed here so they are inspected (not deleted) at Milestone 3:

- `src/ml/BayesianNeuralOperator.cpp` (192 LOC)
- `src/ml/GraphNeuralOperator.cpp` (1583 LOC)
- `src/ml/NeuralReducedOrderModel.cpp` (1469 LOC)
- `src/ml/NeuralSurrogates.cpp` (1873 LOC)
- `src/ml/MLModelRegistry.cpp` (767 LOC)
- `src/ml/FourierNeuralOperator.cpp` (2040 LOC) — the roadmap says explicitly
  "The existing FNO is fine; do not extend it on this branch." Inspection only.

A definitive stub-vs-real classification of `src/ml/` requires linking and
running the test suite, which is gated on the Dockerfile fix above.

## 4. Non-runnable Python scripts

Per the roadmap "What is not real" section: "The `fsrm` Python package — imported
throughout `scripts/` but does not exist in the repo."

Verification:

- No directory named `fsrm/`, `python/fsrm/`, or `src/fsrm/` exists in the repo.
- No `pyproject.toml` exists at repo root.
- No `setup.py` or `setup.cfg` exists at repo root.

The following 23 scripts under `scripts/` import from `fsrm.*` and therefore cannot
be executed in any checkout of this repository as-is:

```
scripts/compare_cannikin_1971.py
scripts/fetch_dprk_2017_waveforms.py
scripts/fetch_lop_nor_1996_waveforms.py
scripts/fetch_lop_nor_2020_waveforms.py
scripts/fetch_ojai_2023_waveforms.py
scripts/fetch_us_divider_1992_waveforms.py
scripts/invert_dprk_2017_moment_tensor.py
scripts/invert_lop_nor_2020_moment_tensor.py
scripts/invert_ojai_2023_moment_tensor.py
scripts/invert_us_divider_1992_moment_tensor.py
scripts/model_degelen_1988_nuclear_test.py
scripts/model_dprk_2017_nuclear_test.py
scripts/model_lop_nor_1996_nuclear_test.py
scripts/model_lop_nor_2020_nuclear_test.py
scripts/model_ojai_2023_earthquake.py
scripts/model_us_divider_1992_nuclear_test.py
scripts/plot_dprk_2017_comparison.py
scripts/plot_lop_nor_1996_comparison.py
scripts/run_explosion_example.py
scripts/run_fsrm_divider_1992_far_field.py
scripts/run_fsrm_dprk_2017_far_field.py
scripts/run_fsrm_dprk_2017_gpu.py
scripts/run_lop_nor_1996_numerical.py
```

Submodules these scripts try to import (none of which exist in the repository):

- `fsrm.source_physics`
- `fsrm.velocity_models` (with classes `LopNorVelocityModel`, `PunggyeRiVelocityModel`,
  `NTSVelocityModel`, `GenericGraniteVelocityModel`, `DegalenVelocityModel`,
  `SoCalVelocityModel`)
- `fsrm.propagation` (with helpers `regional_travel_times`, `generate_synthetic`,
  `mt_mantap_collapse_signal`)
- `fsrm.signal_processing` (`time_relative`, `spectral_amplitude`, `envelope`)
- `fsrm.data_fetching` (`fetch_raw_waveforms`, `process_trace`,
  `fetch_waveforms_for_inversion`)
- `fsrm.fetch_plotting`
- `fsrm.inversion`
- `fsrm.run_utils` (`read_config`)

Per roadmap rule 3 ("Cut, don't add, when in doubt") and the explicit out-of-scope
clause, these scripts will be either deleted, marked `BROKEN` at the top of the
file, or replaced by direct standalone scripts (importing only `obspy`, `numpy`,
`matplotlib` per the M1.4 directive). That work happens in Milestone 1 and
Milestone 3, not in Milestone 0.

## 5. Summary of cleanup work owed to Milestone 3

Driven by this baseline, Milestone 3 must:

1. Repair `Dockerfile.ci` so the integration test suite can be executed and a
   real pass/fail table can replace section 1.2 of this document.
2. Re-run all 106 test source files and update section 1.
3. Resolve each entry in section 2 (Unimplemented headers): delete or move to
   `archive/aspirational/`, except for legitimate pure-virtual interfaces.
4. Delete or move to `archive/aspirational/` all five files in
   `src/experimental/`.
5. Either create a real `fsrm` Python package with `pyproject.toml` (in its
   own session, per the roadmap's "What is explicitly out of scope" clause)
   or delete or `BROKEN`-mark the 23 scripts in section 4.
6. Update `README.md` to remove all capability claims that depend on items
   removed by steps 3-5.

No claim in `README.md` should be added or strengthened until Milestone 1
(Sedan 1962) has met its exit criterion.
