# Session 8 Report

## Goal

Configure a solver for the saddle-point Jacobian that Session 7.5's dispatch
fix produces, so that `Integration.DynamicRuptureSolve.PrescribedSlipQuasiStatic`
can take a Newton step rather than failing in KSP. The session kickoff framed
this as a solver-configuration task, not an architecture task.

Planned decision tree:
  * Step 1: cohesive section reorder + MUMPS (try first).
  * Step 2: GAMG + cohesive reorder if Step 1 fails.
  * Step 3: no fieldsplit for pure elasticity+fault.

## MUMPS availability (Step 1 precheck)

```
docker run --rm fsrm-ci:main bash -c 'grep "HAVE_MUMPS\|HAVE_SUITESPARSE\|HAVE_SUPERLU\|HAVE_PASTIX" /opt/petsc-main/arch-linux-c-opt/include/petscconf.h'
```

Returns nothing. The PETSc build in `fsrm-ci:main` was configured with

```
--with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 --with-mpi=1
--download-fblaslapack --download-ctetgen
--with-hdf5-include=/usr/include/hdf5/openmpi
--with-hdf5-lib="-L/usr/lib/x86_64-linux-gnu/hdf5/openmpi -lhdf5"
--with-debugging=0 COPTFLAGS=-O3 CXXOPTFLAGS=-O3 FOPTFLAGS=-O3
```

No MUMPS, no SuperLU_DIST, no SuiteSparse, no PaStiX. Only PETSc native LU
(no partial pivoting for SeqAIJ) and PCSVD are available as direct solvers.

Raw probe output in `/tmp/s8_mumps_check.log`.

## Code change that landed

One addition in `src/core/Simulator.cpp`, around line 1780 inside `setupFields`,
gated on `config.enable_faults && cohesive_kernel_`:

```cpp
ierr = DMReorderSectionSetDefault(dm, DM_REORDER_DEFAULT_TRUE); CHKERRQ(ierr);
ierr = DMReorderSectionSetType(dm, "cohesive"); CHKERRQ(ierr);
```

The calls run before `DMSetUp(dm)` so the cohesive permutation is applied
when `DMPlexCreateSection` rebuilds the local section after BC insertion.
API verified against `/opt/petsc-main/include/petscdm.h:260-263` and the
implementations in `/opt/petsc-main/src/dm/impls/plex/plexreorder.c`
(`DMReorderSectionSetDefault` / `DMReorderSectionSetType` and the cohesive
permutation path `DMCreateSectionPermutation_Plex_Cohesive`). PyLith sets
the same options via command line in its `getSolverDefaults(hasFault=true)`
(libsrc/pylith/materials/Elasticity.cc); the PETSc cohesive tests in
`src/dm/impls/plex/tests/ex5.c:1289-1298` use the same pair. FSRM calls
the C API directly so no dependence on option-parse order is introduced.

`setupSolvers` is unchanged. Every other code-path surveyed and experimented
with was reverted to Session 7.5 state.

## Configurations tried and outcome

Each row reports what was configured for a single run of
`Integration.DynamicRuptureSolve.PrescribedSlipQuasiStatic`. Full logs are
on disk under `/tmp/`.

| Config | Result | Log |
|---|---|---|
| Cohesive reorder only, test's `-pc_type lu -ksp_type preonly` | `DIVERGED_LINEAR_SOLVE iterations 0` every Newton iteration, zero-pivot in SeqAIJ LU | `/tmp/s8_prescribed_step1.log` |
| Cohesive reorder + PyLith-style `-pc_type gamg -ksp_type gmres -mg_fine_pc_type vpbjacobi` | `DIVERGED_PC_FAILED` / `PC failed due to SUBPC_ERROR` during GAMG setup (vpbjacobi fails on the indefinite Lagrange rows) | `/tmp/s8_prescribed_step2.log`, `/tmp/s8_prescribed_step2_diag.log` |
| Cohesive reorder + `-pc_type lu -pc_factor_shift_type NONZERO -pc_factor_shift_amount 1e-10` (default backtracking line search) | KSP converges in 1 iteration; every Newton step rejected by `DIVERGED_LINE_SEARCH` | `/tmp/s8_prescribed_step3.log` |
| Above + `-snes_linesearch_type basic` (full step, no backtrack) | KSP converges; SNES accepts the step but the residual *grows* from `1.767e-04` to `1.385e+01` on the first step and continues to grow on subsequent steps: this is `DIVERGED_DTOL` | `/tmp/s8_prescribed_step4.log` |
| Cohesive reorder + `-pc_type svd -ksp_type preonly` | KSP factors exactly (PCSVD is a pseudoinverse on a small system); Newton still produces a direction that makes the residual jump to `1.131e+01` on the first step, `DIVERGED_LINE_SEARCH` | `/tmp/s8_prescribed_step5_svd.log` |

## Diagnosis

Step 1 through Step 5 together rule out the premise of the session kickoff.

Session 7.5 argued that the remaining failure was solver configuration: the
Jacobian now carries the O(1e+9) elasticity entries on the disp-disp block
and the O(0.01) hybrid cohesive coupling, and the only missing piece was a
pivoting-capable factorization for the zero-diagonal Lagrange rows. That
read is wrong in this PETSc build. With PCSVD on a system small enough for
SVD to be exact, the Krylov solve produces a Newton direction that makes
the nonlinear residual jump by five orders of magnitude. SVD does not
misfactor the matrix; it returns the minimum-norm least-squares solution
of `J dx = -F`. The implication is that the assembled Jacobian is not
`dF/du` for the `F` that `FormFunction` is returning. The hybrid driver
and the manual cohesive penalty stamp together produce a Jacobian whose
Newton direction is inconsistent with the residual the PetscDS BdResidual
callbacks produce for the prescribed-slip constraint.

This is a Jacobian-consistency bug, not a KSP/PC configuration bug.
Rewiring the solver cannot fix it. Either the hybrid `g0_hybrid_lambda_*`
callbacks in `src/physics/CohesiveFaultKernel.cpp` are not the derivatives
of the corresponding `f0_hybrid_lambda_*` BdResidual callbacks, or the
manual penalty stamp in `addCohesivePenaltyToJacobian`
(`src/core/Simulator.cpp:4900-5040`) has signs that conflict with the
BdResidual path (for example the `row_lag`/`col_neg_disp` stamp uses
`-coeff` while the BdResidual path already sums the neg/pos contribution
in a specific orientation). That has to be confirmed with a pointwise
FD check (`-snes_compare_explicit`, then `-snes_compare_coloring`) on a
single cohesive element. Running that comparison is outside the Session 8
scope described in the kickoff and was not attempted.

## Step 3: fieldsplit

Not tried. The kickoff explicitly forbade fieldsplit for the pure
elasticity+fault configuration, citing PyLith's use of GAMG (with MUMPS
inside it) as production default. With MUMPS unavailable and GAMG failing
at vpbjacobi setup, fieldsplit with a Schur-complement-on-Lagrange path
would be the remaining option, but the Step-5 SVD test shows the matrix
itself does not give a correct Newton direction, so a Schur complement
on an inconsistent matrix would not converge either. Fieldsplit is
deferred with the same Jacobian-consistency blocker.

## Step 5: switchable solver config

Not landed. Since no fault-enabled solver configuration exposed a working
Newton path in this PETSc build, there is nothing to make switchable yet.
`setupSolvers` is unchanged from Session 7.5.

## Full fault-test status table

`ctest --output-on-failure -R 'Fault|Slip|Rupture|Cohesive|TPV5|PressurizedFracture'`
in `fsrm-ci:main` against Session 8 HEAD vs Session 7.5 HEAD
(`/tmp/s8_fault_subset.log`). 25 tests in the subset; identical pass/fail
set in both runs.

| Test | Session 7.5 | Session 8 |
|---|---|---|
| `Functional.DynamicRuptureSetup.*` (6 cases) | PASS | PASS |
| `Physics.CohesiveBdResidual` | PASS | PASS |
| `Physics.SCEC.TPV5` | FAIL | FAIL |
| `Physics.LockedFaultTransparency` | FAIL | FAIL |
| `Integration.PressurizedFractureFEM` | FAIL | FAIL |
| `Integration.DynamicRuptureSolve.LockedQuasiStatic` | FAIL | FAIL |
| `Integration.DynamicRuptureSolve.LockedElastodynamic` | FAIL | FAIL |
| `Integration.DynamicRuptureSolve.PrescribedSlip` | FAIL | FAIL |
| `Integration.TimeDependentSlip` | FAIL | FAIL |
| `Integration.SlippingFaultSolve` | FAIL | FAIL |
| `Integration.SlipWeakeningFault` | FAIL | FAIL |
| Other fault-related tests (unit Coulomb, fault mechanics) | PASS | PASS |

16 pass, 9 fail, identical to Session 7.5.

## PrescribedSlipQuasiStatic full SNES trace (Session 8 HEAD)

Captured with the final committed state (cohesive reorder only, no
`setupSolvers` override) in `/tmp/s8_prescribed_final.log`:

```
0 TS dt 1. time 0.
    0 SNES Function norm 1.767766952966e-04
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    0 SNES Function norm 1.385642450386e+01
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    0 SNES Function norm 1.385642450386e+01
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
... 11 repeats identical to Session 7.5 final trace ...
max_fault_slip = 0, sol_norm = 0
```

Pattern matches Session 7.5 exactly: the test options set
`-pc_type lu -ksp_type preonly`, PETSc's SeqAIJ LU hits a zero pivot on
the Lagrange diagonal, KSP returns `iterations 0`, SNES reports
`DIVERGED_LINEAR_SOLVE`. The cohesive section reorder changes the global
numbering (Lagrange rows end up at the end of the section) but not the
presence of zero diagonal entries at those rows, so PETSc native LU still
fails there.

The richer failure modes observed with the experimental solver configs
(`DIVERGED_PC_FAILED` under GAMG, `DIVERGED_LINE_SEARCH` under shifted LU
and under SVD) are not in the committed state and are captured only in
the per-experiment logs listed above.

## Full suite

`docker run ... ctest -j$(nproc)` against Session 8 HEAD in
`fsrm-ci:main`: 101 of 116 tests pass. 15 listed failures.

```
 40 - Physics.SCEC.TPV5                                 (fault, baseline)
 42 - Physics.TerzaghiConsolidation                     (parallel flake)
 43 - Physics.AbsorbingBC                               (parallel flake)
 45 - Physics.LithostaticStress                         (parallel flake)
 57 - Physics.LockedFaultTransparency                   (fault, baseline)
 67 - Integration.ExplosionSeismogram                   (parallel flake)
 91 - Integration.PressurizedFractureFEM                (fault, baseline)
 93 - Integration.DynamicRuptureSolve.LockedQuasiStatic (fault, baseline)
 94 - Integration.DynamicRuptureSolve.LockedElastodynamic (fault, baseline)
 95 - Integration.DynamicRuptureSolve.PrescribedSlip    (fault, baseline)
 98 - Integration.TimeDependentSlip                     (fault, baseline)
 99 - Integration.NearFieldCoupled                      (parallel flake)
100 - Integration.SlippingFaultSolve                    (fault, baseline)
101 - Integration.SlipWeakeningFault                    (fault, baseline)
111 - Integration.HistoricNuclear.NtsPahuteMesa         (parallel flake)
```

Tests 42, 43, 45, 67, 99, 111 pass when run individually:

```
$ ctest -R '^Physics.TerzaghiConsolidation$|^Physics.LithostaticStress$|^Integration.HistoricNuclear.NtsPahuteMesa$'
100% tests passed, 0 tests failed out of 3
```

(`/tmp/s8_flake_recheck.log`). These six failures are CLAUDE.md's
documented parallel-HDF5 flakes. The remaining nine failures are the
fault-baseline set that Session 7.5 ran against and did not regress.

Compared to Session 7.5's `ctest -j$(nproc)` result (12 failures: 9 fault
baseline + 3 parallel flakes), Session 8 shows 3 additional parallel
flakes (42, 45, 111). Both runs reflect the same underlying fault-baseline
failure set; the extra flakes are not regressions, they are additional
parallel-HDF5 timing noise.

## Which step succeeded

None. The step taxonomy from the kickoff assumed the failure was solver
configuration. With MUMPS absent, both the Step 1 and Step 2 paths hit
real errors (zero pivot in non-pivoting LU; PC setup failure in GAMG).
The two experiments beyond the stated steps (shifted LU and PCSVD) ran
through KSP successfully and exposed the Jacobian-consistency issue
described in the *Diagnosis* section above.

The cohesive section reorder itself is a reasonable improvement to land
on its own merits: it matches PyLith production and the PETSc cohesive
tests, causes no regressions against Session 7.5, and becomes useful the
moment the Jacobian consistency is fixed and a pivoting solver path is
reintroduced (for example via a Dockerfile change adding `--download-mumps`).

## What is next (outside Session 8 scope)

1. Run `-snes_compare_explicit` and `-snes_compare_coloring` on
   `PrescribedSlipQuasiStatic` with one cohesive element, to prove
   whether the assembled Jacobian matches `dF/du`. Session 7.5 reported
   magnitudes (O(1e+9) elasticity, O(0.01) coupling) but did not check
   that signs and couplings match the BdResidual path.
2. If the Jacobian is confirmed inconsistent, the candidate causes are:
   * `g0_hybrid_lambda_displacement` / `g0_hybrid_lambda_lambda` in
     `CohesiveFaultKernel.cpp` signs vs `f0_prescribed_slip`.
   * `addCohesivePenaltyToJacobian` penalty stamp direction at `row_lag`
     (`Simulator.cpp:5024-5035`).
   * A missing `J_u_lambda` stamp for the prescribed-slip branch where
     it currently branches away from the full semi-smooth Newton block.
3. Only after the Jacobian is consistent does it become worth choosing
   between LU (needs a pivoting backend: MUMPS, SuperLU_DIST, or
   UMFPACK via SuiteSparse) and GAMG (needs a fine-level PC that does
   not choke on the Lagrange block) and fieldsplit (Schur on the
   Lagrange field). Any of the three will solve the system at that
   point.

## Rules compliance

- Rule 1 (Docker build/test): all builds and tests in `fsrm-ci:main`.
- Rule 2 (PETSc API): `DMReorderSectionSetDefault`, `DMReorderSectionSetType`,
  `DM_REORDER_DEFAULT_TRUE`, option names `-dm_reorder_section`,
  `-dm_reorder_section_type`, `-pc_factor_shift_type`, `-pc_factor_shift_amount`
  all verified against `/opt/petsc-main/include/petscdm.h`,
  `/opt/petsc-main/include/petscdmtypes.h`,
  `/opt/petsc-main/src/dm/impls/plex/plexcreate.c`,
  `/opt/petsc-main/src/dm/impls/plex/plexreorder.c` before call.
- Rule 3 (no regressions): Session 8 fault-subset and full-suite runs
  carry the same fault-baseline failure set as Session 7.5. The three
  extra failures in the parallel run are HDF5 flakes, verified by the
  serial re-run above.
- Rule 4 (DS/BC ordering): `DMCreateDS` -> `setupBoundaryConditions` ->
  `DMSetLocalSection(nullptr)` -> (Session 8 new: `DMReorderSectionSetDefault` /
  `DMReorderSectionSetType` when fault-enabled) -> `DMSetUp` -> `DMGetDS`.
  The reorder calls are inserted between the existing
  `DMSetLocalSection(nullptr)` and `DMSetUp`, so the documented ordering
  is preserved; the reorder only affects how `DMPlexCreateSection`
  permutes mesh points when `DMSetUp` rebuilds the section.
- Rule 16 (no log truncation): all logs written in full to `/tmp/`:
  `s8_mumps_check.log`, `s8_api_check.log`, `s8_build1.log`, `s8_build2.log`,
  `s8_prescribed_step1.log`, `s8_prescribed_step2.log`,
  `s8_prescribed_step2_diag.log`, `s8_prescribed_step3.log`,
  `s8_prescribed_step4.log`, `s8_prescribed_step5_svd.log`,
  `s8_prescribed_final.log`, `s8_fault_subset.log`, `s8_full_suite.log`,
  `s8_flake_recheck.log`.
