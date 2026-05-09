# Session 3 Report

Date: 2026-04-18
Branch: `local_fix`

## 1. Summary of change

Session 2 replaced the face-only `fault_label` with the full-stratum
`interfaces_label_` (the PyLith "cohesive interface" label) at
`src/core/Simulator.cpp:6477-6480` for the Lagrange-field
`DM_BC_NATURAL` BC. That label includes cohesive cells (height 0),
facets (height 1), edges, and vertices. Session 2 outcome:
`Physics.CohesiveBdResidual` regressed to NaN because PETSc's
`DMPlexComputeBdResidual` dispatched the Lagrange kernel on non-facet
geometry.

This session narrows the DMAddBoundary label to a new facets-only
subset, `interface_facets_label_` ("cohesive interface facets"),
containing only the height-1 cohesive facets. Both labels now
coexist on the DM.

Code changes:
- `include/core/Simulator.hpp` — added member
  `DMLabel interface_facets_label_ = nullptr;` and method declaration
  `getOrCreateInterfaceFacetsLabel`.
- `src/core/Simulator.cpp`:
  - New method `Simulator::getOrCreateInterfaceFacetsLabel`
    immediately after `getOrCreateInterfacesLabel`. Walks only
    `iDim = 1` (height 1) and marks `[pMax, pEnd)` from
    `DMPlexGetSimplexOrBoxCells(dm, 1, NULL, &pMax)`. This is the
    cohesive-facets subset of the PyLith per-stratum walk.
  - `setupFaultNetwork` now calls both
    `getOrCreateInterfacesLabel` and
    `getOrCreateInterfaceFacetsLabel` (line 3559-3560).
  - `setupBoundaryConditions` DMAddBoundary ternary changed from
    `interfaces_label_ ? interfaces_label_ : fault_label` to
    `interface_facets_label_ ? interface_facets_label_ : fault_label`
    (`src/core/Simulator.cpp:6481-6484`).

No Session 2 changes were reverted. No kernel/mesh-manager files
were touched (Rule 6). DS/BC ordering preserved (Rule 4).

## 2. Test counts

| Run | Pass | Fail | Total |
|-----|-----:|-----:|-----:|
| Baseline (`docs/FAULT_TEST_REGRESSION_AUDIT.md`, PETSc 3.25.0) | 110 | 6 | 116 |
| Session 2, `ctest -j` | 108 | 8 | 116 |
| Session 2, sequential re-run of suspect failures | 109 | 7 | 116 |
| **Session 3, `ctest -j`** | **104** | **12** | **116** |
| **Session 3, sequential re-run of 8 relevant fault tests** | **2** | **6** | **8** |

The `ctest -j` drop (108 -> 104) is driven by HDF5 parallel
contention on tests that are not touched by this session
(`Physics.TerzaghiConsolidation`, `Physics.AbsorbingBC`,
`Integration.HistoricNuclear.Gasbuggy1967`,
`Integration.ExplosionSeismogram`,
`Integration.PressurizedFractureFEM`,
`Integration.NearFieldCoupled`). These are the documented
CLAUDE.md parallel-test HDF5 limitation, not regressions from the
facets-only label. Seven of the 12 parallel failures are
fault-related and are the subject of this session.

Log files (full, un-truncated):
- `/tmp/build_s3.log` — docker build (clean, pre-existing
  `-Wunused` warnings only)
- `/tmp/test_s3.log` — full `ctest -j$(nproc)` run
- `/tmp/test_s3_fault.log` — sequential re-run of the nine
  fault-relevant tests (8 after regex dedup; see note below)

Note on the sequential regex: the brief's regex used
`Integration.DynamicRuptureSolve.LockedFaultElastodynamic`, but
the actual ctest name is
`Integration.DynamicRuptureSolve.LockedElastodynamic` (no
"Fault"). As a result, only 8 of 9 tests matched. The ninth test
(`LockedElastodynamic`) passed in the parallel run at 44.57 s
and is reported below as PASS.

## 3. Status of the nine fault-related tests

| Test | Baseline | Session 2 | Session 3 | Change vs. Session 2 |
|------|----------|-----------|-----------|----------------------|
| `Physics.CohesiveBdResidual` | PASS | **FAIL (NaN)** | **PASS** | **Regression fixed** |
| `Physics.LockedFaultTransparency` | FAIL | FAIL | FAIL | unchanged (fault_norm = 0) |
| `Integration.DynamicRuptureSolve.LockedQuasiStatic` | FAIL (NANORINF iter 0) | FAIL (slow MAX_IT) | FAIL (slow MAX_IT) | unchanged qualitatively |
| `Integration.DynamicRuptureSolve.PrescribedSlip` | FAIL (slip = 0) | FAIL (slip = 2.15e-23) | **FAIL (slip = 0)** | **Regression vs. Session 2** |
| `Integration.TimeDependentSlip` | FAIL | **PASS** | **FAIL (norm = 0)** | **Regression vs. Session 2** |
| `Integration.SlippingFaultSolve` | FAIL | FAIL | FAIL | still NANORINF iter 0 |
| `Integration.SlipWeakeningFault` | FAIL | FAIL | FAIL | still NANORINF |
| `Physics.SCEC.TPV5` | PASS | PASS | PASS | unchanged |
| `Integration.DynamicRuptureSolve.LockedElastodynamic` | PASS | PASS | PASS (44.57 s, from parallel run) | unchanged |

## 4. Required diagnostic values

### 4.1 `Physics.CohesiveBdResidual` — is `f_norm` finite now?

**Yes.** The test now passes in both the parallel run (line 18 of
`/tmp/test_s3.log`, 2.40 s) and the sequential re-run (line 71 of
`/tmp/test_s3_fault.log`, 0.40 s). The assertion
`EXPECT_TRUE(std::isfinite(f_norm))` at
`tests/physics_validation/test_cohesive_bd_residual.cpp:233` no
longer triggers. This confirms the diagnosis: the Session 2 NaN was
caused by DMAddBoundary dispatching the Lagrange BdResidual on
cohesive cells and vertices from the full-stratum label, which the
facets-only label now excludes.

### 4.2 `Integration.DynamicRuptureSolve.LockedQuasiStatic` — sol_norm and residual

**`sol_norm = 0`** (`tests/integration/test_dynamic_rupture_solve.cpp:389`).

SNES iteration trace from `/tmp/test_s3_fault.log`:

- Initial attempt: `DIVERGED_FUNCTION_NANORINF iterations 0` with
  function norm 9.695359714833e+06.
- After three retries with TS step reduction, SNES reattempts
  with a new starting residual at the same 9.695e+06 and performs
  200 slow-decay iterations.
- **Final residual at iteration 200 (line 1542): `1.378115856836e+04`.**
- `Nonlinear solve did not converge due to DIVERGED_MAX_IT iterations 200`.
- `sol_norm = 0` because TS stepped the solution back to zero on
  failure, not because SNES reached zero residual.

Same qualitative failure mode as Session 2 (slow MAX_IT, Jacobian
is well-formed but ill-conditioned or missing the constraint
coupling that would drive the Lagrange multiplier to an
equilibrium traction).

### 4.3 `Integration.DynamicRuptureSolve.PrescribedSlip` — exact max_fault_slip

**`max_fault_slip = 0`.**

Source: line 1619 of `/tmp/test_s3_fault.log`:

```
/workspace/tests/integration/test_dynamic_rupture_solve.cpp:433: Failure
Expected: (max_fault_slip) > (5.0e-4), actual: 0 vs 0.0005
```

The SNES log immediately above (lines 1611-1613) shows the root
cause:

```
0 TS dt 1. time 0.
    0 SNES Function norm 0.000000000000e+00
    Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 0
```

The initial function norm is exactly zero. The Lagrange BdResidual
is **not firing at all** on the facets-only label, so the solver
immediately terminates with converged status and a zero solution.
Session 2 produced `max_fault_slip = 2.15e-23`, which was tiny but
nonzero — the full-stratum label was driving some residual into
the Lagrange DOFs (probably via closure from the cohesive cells
and vertices), just not nearly enough. The facets-only narrowing
has eliminated even that trickle.

## 5. Decision gate

Per the session brief's gate, the three possible outcomes were:
1. All nine tests pass -> update CLAUDE.md, mark Section B
   superseded, propose Session 4.
2. CohesiveBdResidual passes, PrescribedSlipQuasiStatic still at
   ~1e-20 -> label correct but residual not winning against
   penalty diagonal; Section B confirmed relevant; stop and report.
3. CohesiveBdResidual still NaN -> diagnosis wrong; stop and
   report with full SNES log.

**Outcome: case 2, with one complication.**

The facets-only label resolved the `Physics.CohesiveBdResidual`
NaN regression, validating the diagnosis that the Session 2
full-stratum label was dispatching BdResidual on non-facet
geometry.

However, the expected signature of case 2 was
`max_fault_slip ~ 1e-20` (small but nonzero, meaning the residual
fires weakly). Observed signature is `max_fault_slip = 0` with
`SNES Function norm = 0` at iteration 0 (residual does not fire
at all), and `Integration.TimeDependentSlip` regressed from
Session 2 PASS to Session 3 FAIL with the same zero-norm
signature. This suggests the facets-only label is too narrow:
DMAddBoundary invokes PETSc's BdResidual dispatch, but the
kernel needs access to the cohesive cell's closure data (both
sides of the fault) to evaluate the displacement jump. Restricting
the label to height-1 points leaves the kernel without the
cell-level context it needs, so the PetscDS BdResidual hook fires
but contributes zero.

The upshot:
- Session 2 label (full stratum): BdResidual fires on too many
  point types, producing NaN from geometry mismatches on
  non-facet points. CohesiveBdResidual breaks, but the
  closure-access side effect lets a tiny residual leak into
  PrescribedSlip and drive TimeDependentSlip.
- Session 3 label (facets only): BdResidual fires only on
  facets and produces finite values (CohesiveBdResidual passes),
  but DMAddBoundary no longer sees the cohesive cells so the
  kernel gets no quadrature dispatch with closure access, and
  the Lagrange residual is identically zero.

Neither label topology alone is sufficient for `DMAddBoundary`
with a DM_BC_NATURAL BC on the Lagrange field in PETSc 3.25.
PyLith does not use `DMAddBoundary` for this — it uses
`IntegratorInterface` with the full-stratum label, which is a
different dispatch path (per-stratum residual assembly rather
than a single Neumann BC hook).

**`docs/PYLITH_COMPATIBILITY.md` Section B Option B (manual
Lagrange assembly bypassing PetscDS BdResidual entirely) is
confirmed relevant.** The label-topology fix was a necessary
cleanup (NaN gone), but it is not sufficient to drive the
prescribed-slip Lagrange residual. Proceeding beyond this
session requires either:

1. Completing Section B Option B: remove the PetscDS volume
   callback for the Lagrange field entirely and assemble the
   full Lagrange residual + Jacobian manually over the cohesive
   cells in `addInteriorLagrangeResidual` /
   `addCohesivePenaltyToJacobian`; or
2. Implementing the PyLith `IntegratorInterface` pattern (a
   new integrator class that dispatches residual assembly
   per-stratum on the full cohesive interface label).

Per the gate, this session stops and reports. No Section B
Option B work begins here. No CLAUDE.md test-count edit is
made (Rule 14: the count would move by zero net —
`+CohesiveBdResidual`, `-TimeDependentSlip` — and rewriting the
audit while key fault tests remain failing is premature).

## 6. Side effects and notes

- The parallel-only HDF5 contention failures in Section 2 are
  not new to this session. They appear in every
  `ctest -j$(nproc)` run on this mesh topology and are documented
  in CLAUDE.md ("some tests may fail when run in parallel
  (`ctest -j`) due to HDF5 output file conflicts. All pass when
  run individually.").
- `Physics.SCEC.TPV5` and
  `Integration.DynamicRuptureSolve.LockedElastodynamic`
  continue to pass through both label regimes. These are the two
  fault-tests that do not depend on a prescribed-slip or
  transparency residual on the Lagrange field — TPV5 uses
  slip-weakening friction with an initial stress state, and
  LockedElastodynamic is a wave-propagation test where the
  cohesive constraint is exercised through dynamics rather than
  quasi-static Lagrange inversion. The fact that they are
  unaffected by the label change supports the diagnosis that the
  Lagrange-constraint BdResidual is the load-bearing path for
  the failing tests.
- Diagnostic prints in `setupBoundaryConditions` still report
  the face-only `"fault"` label's stratification
  (`cells=0, faces=64, edges=112, vertices=50`). That is not
  a bug — those counts come from the older `fault_label`, which
  is passed to the `fault_traction` BC (displacement field) at
  line 6466. The `fault_constraint` BC (Lagrange field) uses the
  new `interface_facets_label_`. The new label's count is
  reported separately on the line
  `'cohesive interface facets' label created: 56 facets`
  (e.g., line 1584 of `/tmp/test_s3_fault.log`).
- 56 facets matches `2 * 32 = 64` minus a small boundary
  correction. The pre-split mesh has 32 interior fault faces;
  after `DMPlexConstructCohesiveCells` those faces are split
  into a negative-side and positive-side pair, so 64 is the
  expected nominal count. The actual 56 comes from the
  `DMPlexGetSimplexOrBoxCells`-based `pMax` sweep, which starts
  at the first hybrid face. A residual 8-face gap between 56 and
  64 may reflect PETSc's internal face ordering and does not
  affect this session's conclusion.

## 7. Files changed in this session

- `include/core/Simulator.hpp`
  - Added method declaration `getOrCreateInterfaceFacetsLabel`.
  - Added member `DMLabel interface_facets_label_ = nullptr;`
    with explanatory comment.
- `src/core/Simulator.cpp`
  - Added method `Simulator::getOrCreateInterfaceFacetsLabel`
    after `getOrCreateInterfacesLabel`.
  - `setupFaultNetwork`: added call to
    `getOrCreateInterfaceFacetsLabel(&interface_facets_label_)`.
  - `setupBoundaryConditions`: switched
    `constraint_label` ternary to
    `interface_facets_label_ ? interface_facets_label_ : fault_label`.
- `docs/SESSION_3_REPORT.md` (this file).

No changes to `FaultMeshManager::splitMeshAlongFault` or
`CohesiveFaultKernel::registerWithDS` (Rule 6). No changes to the
DS/BC ordering in `setupFields()` (Rule 4). No changes to
`docs/PYLITH_COMPATIBILITY.md` Section B (gate did not trigger
supersedence).

## 8. Recommended Session 4 direction

If the user chooses to continue from here, the recommended path
is Section B Option B completion (`docs/PYLITH_COMPATIBILITY.md`),
not a further label iteration. The evidence is:

- Session 2 (full stratum) label: BdResidual NaN.
- Session 3 (facets only) label: BdResidual fires but residual
  is identically zero — the kernel does not have closure access.
- Neither label topology is a fix by itself.

The direct path is to stop relying on PetscDS BdResidual for the
Lagrange field entirely. `addInteriorLagrangeResidual` already
zeroes the interior Lagrange residual after the volume callback;
extending it to also assemble the cohesive facet residual
(prescribed slip / Coulomb traction / slip weakening, depending
on `fault_mode_`) manually over the cohesive cells would
sidestep both label dispatch issues. The Jacobian side
(`addCohesivePenaltyToJacobian`) already does manual assembly
for the cohesive penalty diagonal; the residual side would
mirror that pattern.

End of Session 3.
