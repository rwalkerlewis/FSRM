# Session 13 Report

## Goal

Remove the geometric auto-detection of buried edges and match PyLith semantics:
a `buried_edges` label is a user-supplied or mesh-generator-supplied input, not
auto-detected. When the DM has no `buried_edges` label, skip the Lagrange
rim-pin BC entirely and rely on the prescribed-slip constraint to regularize
the Lagrange block.

## Code changes

### `src/core/Simulator.cpp`

`getOrCreateBuriedCohesiveLabel` was rewritten. The previous body scanned
`interfaces_label_`, filtered for `SEG_PRISM_TENSOR` / `POINT_PRISM_TENSOR`
cells, computed the mean coordinate of each cone vertex, and labeled the
prism if any cone midpoint sat on the domain bounding box. That geometric
detection flagged 56 of 56 cohesive edge-prisms for the 4x4x4 prescribed-slip
test, which caused the over-pinning that Session 12 documented.

The new implementation:

1. Calls `DMHasLabel(dm, "buried_edges", ...)`. If the label is absent, logs
   `"Session 13: no 'buried_edges' label on DM; skipping Lagrange rim-pin BC
   (PyLith default for faults that cut the full domain)"` and returns
   `*buriedCohesiveLabel = nullptr`. The existing gate at
   `setupBoundaryConditions` (`if (buried_cohesive_label_ && lagrange_field_idx_ >= 0)`
   around `src/core/Simulator.cpp:7592`) already short-circuits the BC
   registration when the label is null, so no other code needed to change.
2. If the label is present, creates the output `buried_cohesive` label,
   fetches the stratum-1 index set of `buried_edges`, walks the support of
   each labeled point, filters for `SEG_PRISM_TENSOR` / `POINT_PRISM_TENSOR`
   support cells, and for each such cell sets `buried_cohesive = 1` if any
   cone entry is in `buried_edges`. This is the PyLith pattern verbatim
   (`libsrc/pylith/faults/FaultCohesive.cc:444-480`).

No other source files were touched. No config-schema field was added; the
label name `buried_edges` matches PyLith's default and is a mesh-generator
convention. If a user needs a different name in the future, that would be a
follow-on config change, but the current tests exercise the default.

## Diagnostics

All numbers are from `fsrm-ci:main` (PETSc 3.25), commit `16fc46e` + Session
13 edits.

### `buried_edges` label presence on test meshes

| Test | Has `buried_edges` label? | Points in `buried_cohesive_label_` |
|------|:-----------------:|----------:|
| `LockedQuasiStatic`            | NO | 0 |
| `LockedElastodynamic`          | NO | 0 |
| `TimeDependentSlip`            | NO | 0 |
| `PrescribedSlipQuasiStatic`    | NO | 0 |
| `SCEC.TPV5`                    | NO | 0 |
| `LockedFaultTransparency`      | NO | 0 |
| `SlippingFaultSolve`           | NO | 0 |
| `SlipWeakeningFault`           | NO | 0 |
| `PressurizedFractureFEM`       | NO | 0 |

All FSRM test meshes are built at runtime by `DMPlexCreateBoxMesh` and split
by `FaultMeshManager::splitMeshAlongFault`. Neither path writes a
`buried_edges` label, so with Session 13's PyLith-semantic gating, every test
mesh takes the "no input label -> skip rim pin" branch. The expected log line
`Session 13: no 'buried_edges' label on DM; ...` appears on every fault test.

### Per-field constraint counts for PrescribedSlip (no rim pin)

```
Local section: total dof=525, constrained dof=228 (free=297)
  field 0: dof=450, constrained=228, free=222
  field 1: dof=75, constrained=0, free=75
```

Matches the Session 11 baseline exactly (as expected once the rim-pin BC is
dropped).

### PCSVD condition number and smallest singular values for PrescribedSlip

```
SVD: condition number 2.533601461364e+17,
     0 of 297 singular values are (nearly) zero
SVD: smallest singular values: 7.880764282420e-08 9.547879833136e-08
                               1.006747692764e-07 1.345527200289e-07
                               2.021749505977e-07
SVD: largest singular values:  1.699263701422e+10 1.756995922183e+10
                               1.759717590135e+10 1.995735090500e+10
                               1.996671590260e+10
```

Identical to Session 11. The Lagrange block is uniformly rank-deficient again
because no DOFs are pinned at the rim.

### PrescribedSlipQuasiStatic SNES trace

```
0 TS dt 1. time 0.
  0 SNES Function norm 2.118962305111e+01
  1 SNES Function norm 1.400732819250e+01
  DIVERGED_MAX_IT iterations 1        (repeats with TS backoff, never converges)
max_fault_slip = 0
```

The prescribed-slip constraint alone does not regularize the Lagrange block
on PETSc 3.25 with the current solver config (cohesive section reorder +
superlu_dist LU). The three outcomes the prompt enumerated for Step D: this
corresponds to outcome 2 ("KSP diverges with condition number similar to
Session 11"; PETSc returns `DIVERGED_MAX_IT` rather than
`DIVERGED_LINEAR_SOLVE` here, but the mechanism is the same).

### Fault-subset table vs Session 12

Full subset via `ctest -R '...fault-pattern...'` sequential:

| Test | Session 12 | Session 13 |
|------|:----------:|:---------:|
| Physics.SCEC.TPV5                                      | FAIL | FAIL |
| Physics.LockedFaultTransparency                        | FAIL | FAIL |
| Physics.CohesiveBdResidual                             | pass | pass |
| Functional.DynamicRuptureSetup.MeshSplitting           | pass | pass |
| Functional.DynamicRuptureSetup.LockedFault             | pass | pass |
| Functional.DynamicRuptureSetup.FaultGeometryFromConfig | pass | pass |
| Functional.DynamicRuptureSetup.SlippingFault           | pass | pass |
| Functional.DynamicRuptureSetup.AbsorbingCoexist        | pass | pass |
| Integration.PressurizedFractureFEM                     | FAIL | FAIL |
| **Integration.DynamicRuptureSolve.LockedQuasiStatic**  | **PASS** | **FAIL** |
| **Integration.DynamicRuptureSolve.LockedElastodynamic**| **PASS** | **FAIL** |
| Integration.DynamicRuptureSolve.PrescribedSlip         | FAIL | FAIL |
| Integration.ExplosionFaultReactivation                 | pass | pass |
| **Integration.TimeDependentSlip**                      | **PASS** | **FAIL** |
| Integration.SlippingFaultSolve                         | FAIL | FAIL |
| Integration.SlipWeakeningFault                         | FAIL | FAIL |

**Three fault regressions** (Locked*, TimeDependentSlip). **Zero fault tests
moved from FAIL to PASS.**

All three regressing tests show SNES `DIVERGED_LINEAR_SOLVE iterations 0` at
the first Newton step, i.e. the first KSP solve on the rank-deficient
Lagrange block fails before a direction can be computed. Session 12's rim
pin was what made them converge; removing it drops them back to the
Session 11 failure mode.

### Full-suite count vs Session 12

| Run | Failures (of 116) |
|-----|:----------:|
| Session 12 parallel | 14 |
| Session 13 parallel | 14 |

Same numeric count but the composition shifted: three previously-passing
fault tests regressed (`LockedQuasiStatic`, `LockedElastodynamic`,
`TimeDependentSlip`), and three HDF5-parallel flakes that failed in Session 12
passed this run (NearFieldCoupled, Gasbuggy1967, NtsPahuteMesa). Per the
CLAUDE.md caveat, the HDF5 set shuffles between parallel runs due to output
file conflicts, so that swing is noise.

Session 13 parallel failures:

```
40  Physics.SCEC.TPV5                              (fault-related, baseline)
42  Physics.TerzaghiConsolidation                  (HDF5 parallel noise)
43  Physics.AbsorbingBC                            (HDF5 parallel noise)
45  Physics.LithostaticStress                      (HDF5 parallel noise)
57  Physics.LockedFaultTransparency                (fault-related, baseline)
67  Integration.ExplosionSeismogram                (HDF5 parallel noise)
91  Integration.PressurizedFractureFEM             (fault-related, baseline)
93  Integration.DynamicRuptureSolve.LockedQuasiStatic      (REGRESSION)
94  Integration.DynamicRuptureSolve.LockedElastodynamic    (REGRESSION)
95  Integration.DynamicRuptureSolve.PrescribedSlip  (fault-related, baseline)
97  Integration.TractionBC                          (HDF5 parallel noise)
98  Integration.TimeDependentSlip                   (REGRESSION)
100 Integration.SlippingFaultSolve                  (fault-related, baseline)
101 Integration.SlipWeakeningFault                  (fault-related, baseline)
```

## Decision-gate outcome

The prompt defines four gates:

1. PrescribedSlip passes AND no regressions: faults work; propose Session 14
   as friction-Jacobian port for slipping tests.
2. PrescribedSlip passes but regressions in Locked*/TimeDependent*: revert
   Step B or add `buried_edges` to test mesh setup.
3. PrescribedSlip still fails, no regressions: constraint alone is not
   sufficient regularization; Session 14 adds `MatNullSpace` on the Lagrange
   subspace or fieldsplit.
4. Both regress: something in the gating is wrong. Report and stop.

Session 13 lands between gates 3 and 4. PrescribedSlip still fails (gate 3
symptom), but Locked*/TimeDependent* regress (gate 4 symptom). Three tests
that Session 12 moved from FAIL to PASS are back to FAIL. PrescribedSlip's
condition number and SNES behavior are identical to pre-Session-12 Session 11.

The gating itself is correct for the PyLith pattern, which is what the prompt
asked for. The problem is that the test meshes do not carry a `buried_edges`
label, so every fault test takes the "no rim pin" branch. For faults that
cut the full cross-section of the unit cube (no real buried tips), PyLith
would indeed register no rim pin — but PyLith's default solver configuration
would then handle the saddle point differently than what FSRM currently does
on PETSc 3.25 direct LU.

The empirical finding is that the rim-pin that Session 12 introduced, while
over-pinning by the prompt's own diagnosis, was in practice what made
Locked* and TimeDependentSlip converge. Removing it regresses them without
fixing PrescribedSlip.

Per the gate 4 rule, reporting and stopping without adding further
workarounds.

## Session 14 scope proposal

Two choices; the second is closer to PyLith's own solver-level strategy:

1. Revert Step B gating (keep Session 12's geometric auto-detection) and
   treat Session 13 as a design note on the semantic mismatch. This restores
   the three passing fault tests at the cost of continuing to over-pin
   PrescribedSlip. Not recommended.

2. Keep Session 13's PyLith-semantic gating (so no rim pin when no
   `buried_edges` label is present) and handle the resulting Lagrange-block
   rank deficiency at the solver level. PyLith uses
   `PetscObjectCompose(..., "nearnullspace", ...)` on the displacement field
   (`libsrc/pylith/problems/Problem.cc:701-706`), but relies on the
   prescribed/locked constraint plus the solver setup to regularize the
   Lagrange block. Concrete steps:

   - Attach a `MatNullSpace` or `MatSetNearNullSpace` to the Lagrange
     subblock describing the rigid-body modes that the constraint does not
     eliminate.
   - Switch to the `-pc_type fieldsplit -pc_fieldsplit_type schur` path
     (Session 10 flagged it as the final fallback). On the Schur complement,
     the Lagrange block is solved explicitly as a saddle point and direct
     LU on a rank-deficient block is avoided.
   - Alternatively, write a `buried_edges` label at mesh-generation time in
     `FaultMeshManager::splitMeshAlongFault`, derived from the pre-split
     fault submesh's topological boundary (PyLith mesh generators do this
     upstream). For faults that cut the full domain the label would be
     empty; for partial faults it would carry the true buried tips.

Recommendation: pursue option 2, likely the Schur-complement solver switch
first (smallest code delta, matches PyLith for heterogeneous/partial faults
without requiring mesh-generator changes).

## Artefacts

- `/tmp/s13_build.log` — clean rebuild against `fsrm-ci:main`.
- `/tmp/s13_passing.log` — sequential ctest run of the three
  Session-12-passing fault tests (all three regress).
- `/tmp/s13_prescribed.log` — sequential PrescribedSlip test.
- `/tmp/s13_svd.log` — PrescribedSlip with `-pc_svd_monitor`.
- `/tmp/s13_fault.log` — sequential fault-subset ctest run
  (saved via per-tool persisted-output: `bsi7rv33z.txt`).
- `/tmp/s13_test.log` — parallel full-suite ctest run
  (saved via per-tool persisted-output: `bzbnbp1cf.txt`).

## Commit

```bash
git add -A
git commit -m "session 13: gate Lagrange rim-pin on user-supplied buried_edges label (PyLith semantics)"
git push origin local_fix
```
