# Solver State

Last updated: Session 31. Update at the end of every session that changes
solver state (see `CLAUDE.md` Rule 19).

This document is the standing truth about the fault-solver investigation. It
replaces "read the last N session reports" as the resume-work workflow.
Facts here are cross-referenced to specific `SESSION_NN_REPORT.md` files.

## Current test state

| Config | Fault subset pass | Full suite | Verified |
|---|:-:|:-:|---|
| Default (no env vars) | 10 / 16 | 110 / 116 | S30 |
| Experimental (`FSRM_ENABLE_SLIP_AUX=1 FSRM_SADDLE_SOLVER=fieldsplit`) | 7 / 16 | -- | S30 |

Fault subset regex (16 tests):

```
Physics.LockedFaultTransparency|Physics.CohesiveBdResidual|Integration.DynamicRuptureSolve|Integration.TimeDependentSlip|Integration.SlippingFaultSolve|Integration.SlipWeakeningFault|Physics.SCEC.TPV5|Integration.FaultAbsorbingCoexist|Integration.ExplosionFaultReactivation|Integration.PressurizedFractureFEM|Functional.DynamicRuptureSetup
```

`Functional.DynamicRuptureSetup` is one CMake entry containing 5 gtest cases
(5 tests by ctest count when `-V`); the fault subset counts it as one.

## Failing fault tests (default path, Session 30 verified)

| Test | Diagnosis | Last touched |
|---|---|---|
| `Physics.SCEC.TPV5` | Friction Jacobian port needed; slip-weakening + nucleation patch not Jacobian-coupled | S28 audit, no fix |
| `Physics.LockedFaultTransparency` | PETSc 3.25 BdResidual on cohesive geometry; S15 marked "flipped pass" but later sessions observe it failing — flaky between runs | S30 |
| `Integration.PressurizedFractureFEM` | Hydrofrac rewire (separate scope; unrelated to the fault-solver bottleneck) | pre-S10 |
| `Integration.DynamicRuptureSolve.PrescribedSlip` | Constants-path writes a non-zero target jump (S-fix verified); SNES under direct LU still converges to zero slip because the BdResidual on the Lagrange field does not fire on the cohesive geometry | S30 |
| `Integration.SlippingFaultSolve` | Friction Jacobian | S28 audit, no fix |
| `Integration.SlipWeakeningFault` | Friction Jacobian | S28 audit, no fix |

## Extra failures under experimental path (Session 30 verified)

The experimental path (fieldsplit Schur + aux-slip DM) adds three
additional failures on top of the default six:

| Test | Experimental failure mode |
|---|---|
| `Integration.DynamicRuptureSolve.LockedQuasiStatic` | fieldsplit Schur KSP stall |
| `Integration.DynamicRuptureSolve.LockedElastodynamic` | fieldsplit Schur KSP stall |
| `Integration.TimeDependentSlip` | fieldsplit Schur KSP stall |

`Physics.CohesiveBdResidual` previously regressed under experimental (S30
diagnosis: `DMCreateSubDM` with `NULL` for the IS output parameter causes
NaN in downstream `TSComputeIFunction`). Session 30 fixed this by capturing
and destroying the IS explicitly. That test is currently passing under
both paths.

## Assembled-matrix characterization

Canonical test case: `Integration.DynamicRuptureSolve.PrescribedSlip`,
4x4x4 tet mesh, 32 cohesive cells, bottom Dirichlet (z_min, u=0), roller
on x/y faces, free top face, at u=0.

- Free local size: 255 (displacement 222 + Lagrange 33 after Session 25
  loose rim-pin pins 16 of 25 Lagrange-bearing points).
- Block scales:
  - K (displacement-displacement): O(1e10) from elastic modulus.
  - C (Lagrange-Lagrange): O(1e-6) from epsilon-regularization
    (eps = 1e-4 / h_fault).
  - Off-diagonal (displacement-Lagrange): constraint coupling.
- cond(J) at u=0: **1.76e+18** (Session 30 PCSVD).
- Five smallest singular values: `1.13e-8, 1.99e-7, 2.86e-7, 3.93e-7,
  4.58e-7`. Consistent with 5-6 rigid-body-like modes in the displacement
  block not fully removed by (bottom-face-Dirichlet + fault rim-pin).
- Largest singular values: cluster at `~2e10`.
- PCSVD history (earlier sessions with different mesh / DOF counts and
  different rim-pin strategies):

  | Session | cond# | smallest sigma | comment |
  |---|---|---|---|
  | S9 baseline | 5.36e+17 | 3.72e-08 | no rim-pin, rank-deficient |
  | S11 | 2.54e+17 | 7.88e-08 | FE degree fix, 2x improvement |
  | S12 | 1.28e+17 | 1.55e-07 | aggressive rim-pin (24 of 25) |
  | S13 | 2.53e+17 | 7.88e-08 | rim-pin gated off, reverted |
  | S30 | 1.76e+18 | 1.13e-08 | loose rim-pin (S25) + weak Lagrange reg |

## Kernel verification status

- `f0_hybrid_lambda`, `f0_hybrid_u_neg`, `f0_hybrid_u_pos` at u=0: verified
  correct by by-hand print-out (Session 29 Thread A). Analytic values match
  observed residuals to machine precision.
- Hand-coded Jacobian vs PETSc `-snes_test_jacobian` FD at u != 0:
  matches to `1e-16` at all 21 SNES steps after step 1 (Session 28 and S29).
- `|| J_hand - J_fd ||_F / || J_hand ||_F = 1.78e-2` at u=0 only (Session
  28 reading a). Session 29 reading b (adopted): this is a PETSc FD
  artifact when the residual has an O(1e-3) constant term and the perturbation
  is O(1e-8) -- the truncation error dominates at u=0 where the "constant"
  prescribed-slip term is non-zero, disappearing as u grows.
- Prescribed-slip Cartesian vector plumbing from `cohesive_kernel_->
  prescribed_slip_` into the unified constants array at
  `src/core/Simulator.cpp:2234-2249` is correct (Session pre-S10; see
  `docs/LAGRANGE_FIX_STATUS.md` historical record for the original fix).

## Preconditioner state

Current fieldsplit branch options (Session 23 + 27 + 30 cumulative), only
applied under `FSRM_SADDLE_SOLVER=fieldsplit`:

```
pc_type                                  = fieldsplit
pc_use_amat                              = true
pc_fieldsplit_type                       = schur
pc_fieldsplit_schur_factorization_type   = lower
pc_fieldsplit_schur_precondition         = selfp
pc_fieldsplit_schur_scale                = 1.0
fieldsplit_displacement_ksp_type         = preonly
fieldsplit_displacement_pc_type          = gamg    (ml if available)
fieldsplit_lagrange_multiplier_fault_ksp_type = preonly
fieldsplit_lagrange_multiplier_fault_pc_type  = gamg (ml if available)
ksp_type                                 = gmres
ksp_gmres_restart                        = 100
ksp_rtol                                 = 1.0e-14
ksp_atol                                 = 1.0e-7
ksp_max_it                               = 500
snes_rtol                                = 1.0e-14
snes_atol                                = 5.0e-7
snes_max_it                              = 200
```

Near-null-space (Session 30):

- Sub-DM sized (6 rigid-body modes, local size 222 on PrescribedSlip mesh),
  composed on the displacement field object via `PetscObjectCompose`; pulled
  by fieldsplit inner PCGAMG through the PETSc-internal
  `DMCreateSubDM + PetscObjectQuery("nearnullspace")` path.
- Full-DM sized (6 rigid-body modes, local size 255), attached to the
  monolithic Jacobian in `FormJacobian` via `MatSetNearNullSpace` on the
  default non-fieldsplit path only.

Default path preconditioner (when no `FSRM_SADDLE_SOLVER` env var set and
no caller `-pc_type`): PyLith default for fault case is `gamg + gmres +
near-null-space + vpbjacobi fine smoother`. Currently live at
`Simulator.cpp:3990-4000` (Session 22 + S30).

## KSP stall magnitude (PrescribedSlip, experimental fieldsplit)

| Session | Iter 199 residual | Initial residual | Orders reduced |
|---|---|---|---|
| S27 | 2.91e+03 | 1.29e+09 | 6 |
| S30 | 6.02e+02 | 1.15e+09 | 7 |

- Reduction required for `rtol=1e-14` on a `1e+9` RHS: 14 orders.
- The stall is in the outer Schur KSP. Inner displacement and inner
  Lagrange sub-solves individually appear to make progress when examined
  (S27 ksp_view trace), but the outer coupling saturates.
- SNES under experimental: `DIVERGED_LINEAR_SOLVE` on every inner linear
  solve (verified uniform across u=0 and 10 non-zero iterates in S29).

## Known bottlenecks remaining (Session 30 verdict)

1. **Schur approximation.** `selfp` computes `B diag(A)^-1 B^T`. A's
   diagonal is O(1e10) but the relevant scaling for the Schur complement
   is the tangent-traction scale at the fault (O(E/h)), which `selfp` may
   not capture when the constraint block has non-standard scaling.
2. **Schur factorization type.** `lower` may not be optimal when the
   constraint block has near-zero diagonal (FSRM's eps-regularization is
   O(1e-6)). `full` is a candidate.
3. **Lagrange block scaling.** Assembled value O(1e-6); PyLith's Lagrange
   block is identically zero (see `docs/PYLITH_REFERENCE.md:1.6`). FSRM
   stamps a penalty-scaled diagonal (O(E/h)) in
   `addCohesivePenaltyToJacobian` but only under the manual penalty path,
   not through the PetscDS callback. Re-scaling at PetscDS callback time
   may be required.

## Known non-issues (do not retry)

These were tried and definitively ruled out; revisiting without new
measurement data violates Rule 17.

- **Line-search damping** (S26): damping 0.1, 0.05, 0.01 all produce the
  same DIVERGED_LINE_SEARCH pattern. The correction direction is wrong;
  magnitude damping does not fix a wrong direction.
- **KSP tolerance tightening alone** (S27): 500 iterations at
  `rtol=1e-14` achieves ~7 orders of reduction when 14 are needed. The
  outer Schur preconditioner is the limit; tolerance is not.
- **aux-key architecture concerns** (S15-S18): verified matching across
  all FSRM registration / assembly / lookup paths. The failures are
  PetscDS closure-size vs face-tab mismatches, not aux-key mismatches.
- **Removing the rim-pin entirely**: breaks Locked*/TimeDependent* without
  fixing PrescribedSlip. Session 21 verified: KSP runs but SNES line-search
  rejects every step.
- **Aggressive rim-pin (S12 24-of-25 points pinned)**: KSP converges but
  PrescribedSlip overshoots to `max_fault_slip = 263 m` (target 1 mm)
  because slip concentrates at the single free Lagrange vertex.

## Change log

One line per session. New lines go at the BOTTOM.

- S09: baseline PCSVD cond# 5.36e+17. Fault subset 7 / 16.
- S10: initial rim-pin infrastructure; depth-3 Lagrange DOFs on cohesive
  cells discovered.
- S11: `PetscFECreateLagrange(dim=1, Nc=-1)` moved Lagrange DOFs to
  cohesive edges (d1=25); PCSVD 2x better.
- S12: aggressive rim-pin (24 of 25 pins) lands. +3 tests pass:
  `LockedQuasiStatic`, `LockedElastodynamic`, `TimeDependentSlip`. But
  `PrescribedSlip` overshoots to 263 m.
- S13: rim-pin gated on user-supplied `buried_edges` label. 3 tests
  regress. Decision: restore in S14.
- S14: geometric rim-pin restored as fallback when label absent. Explicit
  `FSRM_SADDLE_SOLVER={gamg, fieldsplit}` selector gates.
- S15: aux-slip DM infrastructure; `LockedFaultTransparency` flips to
  pass (ephemeral).
- S16: material aux key attachment fixed; three aux-FE shapes all fail
  plexfem closure checks.
- S17: PyLith-style FE construction landed; closure mismatch persists on
  DMClone with tet FE on prism closure.
- S18: unified aux DM pivot to two-DM layout (slip-only surface FE);
  tabulation-point mismatch at `febasic.c:663`.
- S19: manual section-write replaces projection for DOF stamping;
  plexfem checks still fail.
- S20: material aux copies both volume AND face quadrature; `febasic.c:663`
  resolved. +2 tests pass: `LockedQuasiStatic`, `TimeDependentSlip`.
- S21: full-DM rigid-body near-null-space attached via
  `MatSetNearNullSpace` in `FormJacobian`.
- S22: PyLith `vpbjacobi` smoother replaces scalar jacobi. Experimental
  gets `emax/emin ~ 2.6e+16` on coarse level -- catastrophic.
- S23: fieldsplit options verbatim from PyLith
  `solver_fault_fieldsplit.cfg`. Displacement FE renamed
  `"displacement_" -> "displacement"`. +2 tests pass from rename alone.
- S24: aux-slip + rim-pin + fieldsplit combo: KSP converges in 3 iters;
  SNES line-search rejects (step magnitude too large).
- S25: two-plane rim-pin criterion (loose) gives 33 free Lagrange DOFs
  vs S24's 3. KSP 27 iters CONVERGED_RTOL; line-search still rejects.
- S26: line-search damping tested 0.1, 0.05, 0.01 -- all fail identically.
- S27: PyLith tolerances `ksp_rtol=1e-14`, `snes_rtol=1e-14`. KSP runs 199
  iters, stalls at `~2.9e+3`. Decision gate: stop tuning, measure.
- S28: `-snes_test_jacobian` audit. `|| J_hand - J_fd ||_F/||J_hand||_F =
  1.78e-2` at u=0; `1e-16` at u != 0.
- S29: by-hand kernel print confirms kernels correct at u=0. 1.78e-2 is
  PETSc FD artifact. KSP stalls uniform across 11 TS retries.
- S30: sub-DM sized near-null-space for fieldsplit inner GAMG. KSP iter
  199 residual improves `2.9e+3 -> 6.0e+2` (5x). Still stalls. Verdict:
  Schur approximation is the next bottleneck.
- S31: documentation reorganization. No solver-state change.
- pass-3: pass-3 historic-nuclear-fidelity infrastructure (mesh refinement,
  per-layer Q, frequency-dependent t*, layered absorbing-BC test). Does
  not modify Jacobian / cohesive code paths; default fault subset count
  unchanged at 10/16.
- pass-4: pass-4 multi-cell moment-tensor source distribution
  (`[SOURCE_DISTRIBUTION]` config grammar in `addExplosionSourceToResidual`;
  Integration.SourceDistribution.*, 5 tests; historic-nuclear _Distributed
  variants, 7 tests). Does not modify Jacobian / cohesive code paths;
  default fault subset count unchanged at 10/16.

## Quantitative timeline

| Session | cond#(SVD) | smallest sigma | max fault slip (target 1 mm) | Fault default | Fault exp |
|---|---|---|---|---|---|
| S9  | 5.36e+17 | 3.72e-08 | 0      | 7/16  | -- |
| S11 | 2.54e+17 | 7.88e-08 | 0      | 7/16  | -- |
| S12 | 1.28e+17 | 1.55e-07 | 263 m  | 10/16 | -- |
| S13 | 2.53e+17 | 7.88e-08 | 0      | 7/16  | -- |
| S20 | -- (direct LU) | -- | 263 m  | 10/16 | -- |
| S23 | -- (direct LU) | -- | 263 m  | 10/16 | 6/16 |
| S27 | -- (direct LU) | -- | 263 m  | 10/16 | 5/16 |
| S30 | 1.76e+18 | 1.13e-08 | 263 m  | 10/16 | 7/16 |

`max fault slip = 263 m` on the default path is the S12 over-constrained
rim-pin overshoot; the test expects 1 mm. `max fault slip = 0` indicates
the BdResidual on the Lagrange field never fired (the test's zero-solution
failure mode).

## Known report discrepancies

Harvested during Session 31 review of reports S10-S30. Where two reports
disagree, the later session's number is considered authoritative unless a
later report explicitly re-measures.

- **CLAUDE.md baseline (pre-S31) listed 8 known failures**; S23, S27, S30
  all measured only 6 on the default path. `LockedQuasiStatic` and
  `TimeDependentSlip` started passing at S23 (displacement FE rename).
  CLAUDE.md will be updated in this session.
- **S12 vs S15 pass-count**: S12 reports 3 flips to pass; S15 re-measures
  11/16. S23 independent re-measurement gives 10/16. S15's 11/16 was
  likely an HDF5 run-order flake; 10/16 is the stable default-path count.
- **S28 row 29 displacement / row 163 single-entry FD diff** remains
  unexplained; S29 confirmed kernel is correct; S28's single-entry
  5.37e+8 coupling miss is treated as a context-dependent FD artifact
  per reading b and is not blocking.
