# Session Runbook

Standard operating procedures for FSRM sessions: env vars, test subsets,
command recipes, decision-branch template, report structure.

## Standard env-var matrix

| Config | Env vars | Purpose |
|---|---|---|
| Default | (none) | Baseline fault subset 10/16; full suite 110/116 |
| Experimental fieldsplit | `FSRM_ENABLE_SLIP_AUX=1 FSRM_SADDLE_SOLVER=fieldsplit` | Under investigation for PrescribedSlip |
| Experimental GAMG | `FSRM_ENABLE_SLIP_AUX=1 FSRM_SADDLE_SOLVER=gamg` | PyLith-style monolithic GAMG path (Session 22) |
| Debug KSP trace | + `FSRM_KSP_VIEW=1` | Enables `-ksp_view -ksp_monitor -ksp_converged_reason -ksp_error_if_not_converged` |
| Jacobian FD check | + `FSRM_SNES_TEST_JAC=1` | Enables `-snes_test_jacobian -snes_test_jacobian_view`, caps SNES at 1 iter |
| SVD spectrum | + `FSRM_SVD_DEBUG=1` | Swaps LU for SVD, prints singular value spectrum. Caps SNES at 1 iter |
| Session-N diagnostic | `FSRM_S<N>_<TAG>=1` | Per-session diagnostic prints (e.g. `FSRM_S29_BYHAND=1`); remove before commit |

Notes:

- `FSRM_SVD_DEBUG=1` sets `-pc_type svd` at test-harness level BEFORE
  `setupSolvers` runs. The fieldsplit branch at
  `src/core/Simulator.cpp:3932` unconditionally resets `-pc_type fieldsplit`,
  so `FSRM_SADDLE_SOLVER=fieldsplit` + `FSRM_SVD_DEBUG=1` gives the
  fieldsplit path, not SVD. For spectrum measurement, run with
  `FSRM_SVD_DEBUG=1` alone (no `FSRM_SADDLE_SOLVER` override).
- `FSRM_KSP_VIEW=1` also sets `-ksp_error_if_not_converged`, which can
  cause runs to abort early. Expected; inspect the saved log.

## Standard test subsets

| Subset | Regex / flag | Contents | Count |
|---|---|---|---|
| Fault subset | see below | Fault-adjacent tests | 16 ctest entries |
| Full suite | no `-R` | All tests | 116 |
| Single test | e.g. `DynamicRuptureSolve.PrescribedSlip` | One test | 1 |

Fault subset regex (keep on a single line for copy-paste):

```
Physics.LockedFaultTransparency|Physics.CohesiveBdResidual|Integration.DynamicRuptureSolve|Integration.TimeDependentSlip|Integration.SlippingFaultSolve|Integration.SlipWeakeningFault|Physics.SCEC.TPV5|Integration.FaultAbsorbingCoexist|Integration.ExplosionFaultReactivation|Integration.PressurizedFractureFEM|Functional.DynamicRuptureSetup
```

## Standard commands

Replace `NN` with the current session number. Always save the full log to
`/tmp` and never pipe through `head|tail|grep` (Rule 16 -- truncating
loses diagnostic output). Use the `Read` tool on the saved log.

### Build

```bash
docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:main \
  bash -c 'cd build && make -j$(nproc)' 2>&1 | tee /tmp/sNN_build.log
```

### Fault subset, default path

```bash
docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:main \
  ctest --output-on-failure -R \
  'Physics.LockedFaultTransparency|Physics.CohesiveBdResidual|Integration.DynamicRuptureSolve|Integration.TimeDependentSlip|Integration.SlippingFaultSolve|Integration.SlipWeakeningFault|Physics.SCEC.TPV5|Integration.FaultAbsorbingCoexist|Integration.ExplosionFaultReactivation|Integration.PressurizedFractureFEM|Functional.DynamicRuptureSetup' \
  2>&1 | tee /tmp/sNN_fault_default.log
```

### Fault subset, experimental fieldsplit

```bash
docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:main \
  bash -c 'FSRM_ENABLE_SLIP_AUX=1 FSRM_SADDLE_SOLVER=fieldsplit \
    ctest --output-on-failure -R \
    "Physics.LockedFaultTransparency|Physics.CohesiveBdResidual|Integration.DynamicRuptureSolve|Integration.TimeDependentSlip|Integration.SlippingFaultSolve|Integration.SlipWeakeningFault|Physics.SCEC.TPV5|Integration.FaultAbsorbingCoexist|Integration.ExplosionFaultReactivation|Integration.PressurizedFractureFEM|Functional.DynamicRuptureSetup"' \
  2>&1 | tee /tmp/sNN_fault_exp.log
```

### Full suite

```bash
docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:main \
  ctest --output-on-failure 2>&1 | tee /tmp/sNN_full.log
```

### Single test with KSP trace

```bash
docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:main \
  bash -c 'FSRM_ENABLE_SLIP_AUX=1 FSRM_SADDLE_SOLVER=fieldsplit FSRM_KSP_VIEW=1 \
    ctest --output-on-failure -R DynamicRuptureSolve.PrescribedSlip -V' \
  2>&1 | tee /tmp/sNN_prescribed_exp.log
```

### SVD spectrum (single test, no fieldsplit override)

```bash
docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:main \
  bash -c 'FSRM_ENABLE_SLIP_AUX=1 FSRM_SVD_DEBUG=1 \
    ctest --output-on-failure -R DynamicRuptureSolve.PrescribedSlip -V' \
  2>&1 | tee /tmp/sNN_svd.log
```

### Build + full-suite single docker invocation

```bash
docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:main bash -c \
  'cd build && make -j$(nproc) && ctest -j$(nproc) --output-on-failure' \
  2>&1 | tee /tmp/sNN_build_and_test.log
```

## Decision-branch template

Every session prompt that proposes a fix or change should include these
four branches in its decision section. The implementing session must
label its result with exactly one branch.

- **Branch A**: Fix works + no regression -> flip gate (if applicable),
  update `docs/SOLVER_STATE.md` "Current test state" and change log, land.
- **Branch B**: Fix works for target but regresses something else -> keep
  opt-in via env var. Update `docs/SOLVER_STATE.md` change log. Document
  what regressed and why in the session report.
- **Branch C**: Fix does not work -> document what was measured and how
  it falls short. Update `docs/SOLVER_STATE.md` "known bottlenecks" if
  new bottleneck is identified, or "known non-issues" if this confirms
  the attempted approach is ruled out.
- **Branch D**: Build break, crash, or unrelated regression -> revert,
  debug the regression first, rerun. Do not commit a broken state.

## Standard report structure

Every `docs/SESSION_NN_REPORT.md` must contain these sections:

1. **Summary** (<= 20 lines). One paragraph plus headline numbers.
2. **Code changes** with `file:line` references. Absent if doc-only.
3. **Measurement results**. Numbers only where possible; tables preferred.
4. **Fault subset table** (experimental and default, if changed).
5. **Full suite count** (if changed).
6. **Branch verdict** (A/B/C/D from decision template).
7. **What this changes in `docs/SOLVER_STATE.md`**. Explicit diff if
   any (test counts, bottleneck list, etc.).
8. **Session N+1 recommendation**. One paragraph.

Keep reports under ~350 lines. Longer than that usually means the session
did too much; break into two sessions.

## Reading order when resuming work

Before starting a new solver session:

1. `docs/SOLVER_STATE.md` -- current truth.
2. `docs/PYLITH_REFERENCE.md` -- if anything about PyLith comes up.
3. `docs/SESSION_RUNBOOK.md` -- this document, for commands.

Avoid reading every prior `SESSION_NN_REPORT.md` in sequence; the three
standing documents should suffice. Reach for a specific report only when
the standing docs refer you to one.

## Smoke test for this runbook

To verify the ctest command syntax from this document is still valid,
run a single known-passing test:

```bash
docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:main \
  ctest --output-on-failure -R 'Physics.CohesiveBdResidual' -V
```

Expected: `1 test passed, 0 failed.` If this fails, the docker image or
build state is broken and no other commands will produce meaningful
results.
