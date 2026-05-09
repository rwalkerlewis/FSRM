# Session Reports Archive

Numbered per-session work logs from the FSRM development history. These
documents are point-in-time snapshots; they were not refreshed after their
session closed and they will not be brought up to date.

For current state, prefer the standing reference documents in `docs/`:

- `docs/SOLVER_STATE.md` for the current fault-solver state and pass/fail.
- `docs/HISTORIC_NUCLEAR_FIDELITY.md` for what each historic-nuclear pass shipped.
- `docs/AXIS_1A_FIDELITY_REPORT.md` for the cross-pass axis-1a verification result.
- `docs/HISTORIC_NUCLEAR_ROADMAP.md` for the forward roadmap.

The session-process runbook (`SESSION_RUNBOOK.md`) lives here because it is a
contributor process doc, not a user-facing runbook.

## Index

| Session | Topic | Track |
|---:|---|---|
| 2 | Initial baseline (2026-04-18) | bring-up |
| 3 | Bring-up continued (2026-04-18) | bring-up |
| 4 | Material label initialization | bring-up |
| 5 | FE construction change | bring-up |
| 6 | Code changes landed | bring-up |
| 7.5 | Volume dispatch port from PyLith | fault solver |
| 8 | PCSVD configuration for saddle-point Jacobian | fault solver |
| 9 | Jacobian-residual consistency bug | fault solver |
| 10 | PyLith fault-rim Dirichlet BC pattern | fault solver |
| 11 | Lagrange FE degree fix | fault solver |
| 12 | Pin Lagrange DOFs at fault rim | fault solver |
| 13 | Remove buried-edge auto-detection | fault solver |
| 14 | Saddle-point solver for fault-enabled problems | fault solver |
| 15 | PyLith aux-field slip pattern port | fault solver |
| 16 | aux-field plumbing continuation | fault solver |
| 17 | aux-field rollback investigation | fault solver |
| 18 | Unify auxiliary DMs for fault-slip aux field | fault solver |
| 19 | Volume-FE aux slip plumbing | fault solver |
| 20 | Material-aux face quadrature + aux-slip gating probe | fault solver |
| 21 | GAMG + rigid-body near-null-space infrastructure | fault solver |
| 22 | GAMG smoother correction | fault solver |
| 23 | PyLith production fieldsplit Schur replication | fault solver |
| 24 | Fieldsplit Schur selfp + rim-pin ON | fault solver |
| 25 | Geometric rim-pin criterion refinement | fault solver |
| 26 | SNES line-search damping for fieldsplit | fault solver |
| 27 | Tighten KSP/SNES tolerances to PyLith values | fault solver |
| 28 | Jacobian audit: FD check | fault solver |
| 29 | Classify u=0 Jacobian diff; KSP stall scope | fault solver |
| 30 | Sub-DM rigid-body near-null-space for fieldsplit | fault solver |
| 31 | Documentation reorganization | docs |
| 32 | Multi-cell moment-tensor source distribution (pass-4) | historic-nuclear axis-1 |

Pass-5 onward switched from numbered session reports to "passes" tracked in
`docs/HISTORIC_NUCLEAR_FIDELITY.md` and `docs/HISTORIC_NUCLEAR_ROADMAP.md`.
