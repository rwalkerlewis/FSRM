# Archived Documentation

Documents in this directory describe abandoned, superseded, or dormant work.
They are kept for historical context and to preserve the decision trail.
None of them describe the current state of FSRM.

For current state, see the standing reference docs in `docs/`:

- `docs/SOLVER_STATE.md` for current fault-solver state.
- `docs/HISTORIC_NUCLEAR_FIDELITY.md` for what each pass shipped.
- `docs/AXIS_1A_FIDELITY_REPORT.md` for axis-1a verification.
- `docs/HISTORIC_NUCLEAR_ROADMAP.md` for the forward roadmap.

## Inventory

### Nuclear Explosion Monitoring early roadmap

- `NEM_BASELINE.md` (was `docs/NEM_BASELINE.md`): Milestone 0 honest-baseline
  inventory from the `nem-roadmap` branch. The roadmap was retired when the
  same scope was reorganized as the historic-nuclear axis-1 program and
  formalized in `docs/HISTORIC_NUCLEAR_ROADMAP.md`. The baseline inventory
  was a useful one-time exercise but was not maintained against the active
  branch and is now stale.
- `NEM_ROADMAP.md` (was `docs/NEM_ROADMAP.md`): The branch roadmap that
  produced `NEM_BASELINE`. Superseded by `docs/HISTORIC_NUCLEAR_ROADMAP.md`,
  which expresses the same goal (validated near-field-to-far-field pipeline
  for historic underground nuclear tests anchored to published observations)
  in the current six-axis fidelity-ladder formulation.

### Fault solver decision trail

- `LAGRANGE_FIX_STATUS.md` (was `docs/LAGRANGE_FIX_STATUS.md`): Decision
  trail from Sessions 10-14's exploration of PETSc 3.25 region-DS
  limitations and the rim-pin + manual-interior-assembly strategy. Marked
  HISTORICAL / SUPERSEDED in Session 31. For current solver state see
  `docs/SOLVER_STATE.md`.
- `FAULT_TEST_REGRESSION_AUDIT.md` (was `docs/FAULT_TEST_REGRESSION_AUDIT.md`):
  Session 1.5 baseline audit walking 10 commits back along `local_fix`.
  Establishes that no last-known-good commit exists within the audited
  window. The six fault tests it tracks were subsequently disabled in
  PR #119 because the underlying PETSc 3.25 cohesive-cell BdResidual
  blocker is not addressable from FSRM. Current state in
  `docs/SOLVER_STATE.md`.

### PyLith feature-parity

- `PYLITH_COMPATIBILITY.md` (was `docs/PYLITH_COMPATIBILITY.md`):
  Feature-parity wishlist describing alignments and gaps between FSRM and
  PyLith. With the fault-solver work paused after PR #119, this wishlist no
  longer describes an active development track. The verified architecture
  reference `docs/PYLITH_REFERENCE.md` (still in `docs/`) remains the
  canonical source for PyLith file:line references.
