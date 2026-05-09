# Session 31 Report -- Documentation reorganization

## Summary

Session 31 is a documentation-only session. Three standing reference
documents were created so future sessions do not re-read the 30 prior
session reports to resume work. No source, header, or test files touched.
`CLAUDE.md` gained Rules 17-20 (formalizing the "no tuning without
measurement" guidance that Sessions 28/29/30 prompts introduced ad hoc)
and a Documentation Layout section. The known-failing test list in
`CLAUDE.md` was reconciled from the stale 8-test version to the verified
6-test Session 30 measurement.

Branch verdict: N/A (no code changes). Full suite not re-run (no code
change possible); the runbook smoke test (`Physics.CohesiveBdResidual`)
passes.

## Documents created

| File | Lines | Purpose |
|---|---:|---|
| `docs/SOLVER_STATE.md` | 277 | Standing truth for fault-solver state: current test counts, matrix characterization, kernel verification, preconditioner state, KSP stall magnitude, bottlenecks, known non-issues, change log S09-S31, quantitative timeline, discrepancy reconciliation. |
| `docs/PYLITH_REFERENCE.md` | 408 | Verified PyLith architecture pins with file:line refs. PyLith commit SHA `52b3012a0f87992544059917a3cfedee12e7d02e`. Covers weak-form key encoding, aux attachment, interface patches, kernel registration, residual hybrid-compute, rim-pin, null-space, solver defaults, fieldsplit options, tolerances. |
| `docs/SESSION_RUNBOOK.md` | 167 | Env-var matrix, test subsets, command recipes, decision-branch template, standard report structure, smoke test. |

## Documents modified

- `CLAUDE.md`: added Documentation Layout section (between Project Summary
  and Codebase Size); reconciled Test Suite known-failures list to the
  verified Session 30 six-failure set; added Rules 17-20.
- `docs/LAGRANGE_FIX_STATUS.md`: prepended a "HISTORICAL / SUPERSEDED"
  header pointing to `SOLVER_STATE.md` for current state. Kept for
  decision-trail value (retained per Session 31 prompt's conditional
  delete: SOLVER_STATE.md now covers every outstanding item, but the
  detailed Phase 1-4 trail and the Options A/B/C/D analysis remain a
  useful historical index into Sessions 10-14).
- `docs/PYLITH_COMPATIBILITY.md`: prepended a one-line note pointing
  readers to `PYLITH_REFERENCE.md` for verified architecture; explicit
  marker that this document is a feature wishlist, not a reference.

## Discrepancies found while harvesting

Catalogued in `docs/SOLVER_STATE.md` "Known report discrepancies" section
in full. Summary:

1. **CLAUDE.md baseline listed 8 known failures**; Sessions 23/27/30 all
   measure 6 on the default path. `Integration.DynamicRuptureSolve.
   LockedQuasiStatic` and `Integration.TimeDependentSlip` started passing
   at Session 23 (from the displacement FE rename, not a kernel fix).
   Reconciled by rewriting the CLAUDE.md known-failures list.
2. **Session 12 claimed 3 test flips to pass**; Session 15 reports 11/16;
   Session 23 remeasures 10/16. 10/16 is the stable default-path count;
   Session 15's 11/16 is treated as an HDF5 run-order flake (CLAUDE.md
   already warns that parallel ctest can produce HDF5-conflict false
   positives).
3. **Session 27 prompt cites PyLith tolerances `ksp_rtol=1e-14,
   snes_rtol=1e-14, ksp_atol=1e-7, snes_atol=5e-7`**. The PyLith source at
   commit `52b3012a0` has tighter *relative* but looser *absolute*
   tolerances: `ksp_rtol=1e-12, ksp_atol=1e-12, snes_rtol=1e-12,
   snes_atol=1e-9` (`libsrc/pylith/utils/PetscOptions.cc:317-330`). Either
   FSRM's Session 27 values came from a different PyLith commit or from a
   benchmark cfg rather than the programmatic defaults. Not a blocker:
   KSP stall is the bottleneck (Session 30 confirmed), not tolerance,
   so this is documented in `PYLITH_REFERENCE.md` but not reconciled
   in code.
4. **Session 28's row-29 displacement / row-163 single-entry FD diff
   (5.37e+08 cross-fault coupling miss)**: remains unexplained but not
   blocking. Session 29 confirms kernels are correct at u=0 and u!=0.
   The single-entry residual FD diff is treated as a context-dependent
   PETSc FD artifact (same reasoning as the 1.78e-2 Frobenius ratio at
   u=0).

## CLAUDE.md edits

- Added "Documentation Layout" section after "## Project Summary"
  listing SOLVER_STATE.md, PYLITH_REFERENCE.md, SESSION_RUNBOOK.md as
  primary reading order, with LAGRANGE_FIX_STATUS.md and
  PYLITH_COMPATIBILITY.md marked secondary.
- Rewrote "### Test Suite" known-failures from the stale 8-test list
  (with Session 18 details that are superseded) to the verified
  6-test Session 30 set + 3-test experimental regression note.
- Added Rules 17 through 20:
  - Rule 17: No option tuning without diagnostic data.
  - Rule 18: No Jacobian-adjacent modifications in a diagnostic session.
  - Rule 19: Every session that changes solver behavior must append to
    `docs/SOLVER_STATE.md` change log.
  - Rule 20: Read `docs/SOLVER_STATE.md` before solver sessions; consult
    `docs/PYLITH_REFERENCE.md` before re-deriving PyLith facts.

## Verification gate results

- `docs/SOLVER_STATE.md` exists, reads cleanly, contains every item from
  the Session 31 prompt template. No "TBD" placeholders remain (every
  number has a session-report citation).
- `docs/PYLITH_REFERENCE.md` exists. Every claim cites a specific
  `file:line_range`. PyLith commit SHA recorded at top
  (`52b3012a0f87992544059917a3cfedee12e7d02e`).
- `docs/SESSION_RUNBOOK.md` exists. Smoke test (ctest
  `Physics.CohesiveBdResidual`) executed at
  `/tmp/s31_smoke.log` -- 1/1 passed.
- `CLAUDE.md` now has 20 numbered rules (1-20) and a Documentation
  Layout section.
- `git diff --name-only origin/local_fix` returns only `CLAUDE.md`,
  `docs/LAGRANGE_FIX_STATUS.md`, `docs/PYLITH_COMPATIBILITY.md`; no
  `src/`, `include/`, or `tests/` files modified. Untracked additions
  are the three new `docs/` files and this report.

## Proposed Session 32 scope

Resume solver work with the Session 30 verdict's three Schur
investigations as candidates. Every candidate must include a
measurement before a tuning decision (Rule 17):

1. `pc_fieldsplit_schur_precondition`: try `a11` (use the Lagrange
   block itself as the Schur approximation; currently `selfp` computes
   `B diag(A)^-1 B^T` which may miss the right scale). Measure
   `sigma_min(Schur)` before and after.
2. `pc_fieldsplit_schur_factorization_type`: try `full` instead of
   `lower`. Measure KSP iteration count and residual trajectory.
3. Lagrange block scaling at assembly: scale the PetscDS g0_lambda_lambda
   callback by `E/h` so the diagonal matches displacement block scale.
   This is a kernel-math change and will need Rule 5 / Rule 6 /
   Rule 18 review before proceeding. Measure cond(J) before and after;
   verify Locked*/TimeDependent* default-path tests do not regress.

Session 32 should pick ONE candidate, measure first, then change. Do
not batch all three.

## Rule compliance

- Rule 31 (Session 31-specific, no code changes): honored.
  `git diff --name-only origin/local_fix` returns only CLAUDE.md and
  two modified docs; no src/include/tests edits.
- Rule 14: CLAUDE.md updated with current known-failures list and
  Documentation Layout section.
- Rule 16: smoke test log saved to `/tmp/s31_smoke.log`.
- Rule 19 (new): N/A, no solver-state change this session. Change-log
  entry S31 added to `docs/SOLVER_STATE.md` recording the doc
  reorganization.
