# Parallel KSP / PC diagnostic baseline logs

This directory holds one PETSc-verbose run log per example that was
listed as `SNES diverges step 0` under the `EXAMPLE_SMOKE_MPI4_KNOWN_BROKEN`
gate in `tests/CMakeLists.txt` at the start of pass-12 followup 3.

## Convention

- One file per example: `<N>_<event>.mpi4.log`.
- The file is produced by `tests/diagnostics/capture_parallel_ksp_baseline.sh`,
  which runs the example's `run.sh` at `MPI_RANKS=4` with
  `FSRM_FINAL_TIME_OVERRIDE` set to the example's smoke budget and
  `PETSC_OPTIONS` set to a verbose monitor string (see the script header
  for the exact options).
- The CTest target `Diagnostic.ParallelKSP.<N>_<event>.BaselineMPI4Capture`
  (label `diagnostic`) regenerates the file and asserts it is non-empty.
  The test succeeds whether or not the example itself converges: the point
  is to capture the failure signature, or the post-fix convergence, exactly.

## What the log must contain

`-snes_view` dumps the full SNES / KSP / PC / Mat tree once after the
first `SNESSolve`: the KSP type (`preonly` in the baseline, `gmres`
post-fix), the PC composition (`bjacobi` with a per-rank `lu`), and the
Jacobian sparsity (`rows`/`cols`/`nonzeros` of the `mpiaij` matrix).
`-snes_monitor` / `-ksp_monitor_true_residual` show the residual
trajectory; in the baseline the KSP norm type is `none`
(`CONVERGED_ITS iterations 1` per Newton step, i.e. one block-Jacobi
sweep), post-fix it is a real Krylov solve. `-snes_converged_reason` /
`-ksp_converged_reason` print the terminal reason code: for the baseline
runs it is `SNES ... DIVERGED_MAX_IT` followed by a
`TSStep ... DIVERGED_NONLINEAR_SOLVE` abort; for the post-fix runs it is
`SNES ... CONVERGED_*` and the run reaches `-ts_max_steps 1` and prints
`Simulation completed successfully`. The two example-local data fixes
(06 mesh path, 17 missing `velocity_model.bin`) show a
`Cannot open ... No such file or directory` abort in the baseline log
and normal startup post-fix.

## Files

| File | Baseline reason code | Post-fix outcome |
|---|---|---|
| `06_gmsh_multimaterial.mpi4.log` | see log | converges |
| `09_gasbuggy_1967.mpi4.log` | see log | converges |
| `17_velocity_model.mpi4.log` | see log | converges |
| `23_milrow_1969.mpi4.log` | see log | converges |
| `24_cannikin_1971.mpi4.log` | see log | converges |
| `25_faultless_1968.mpi4.log` | see log | converges |
| `29_rio_blanco_1973.mpi4.log` | see log | converges |
| `33_dprk_2006.mpi4.log` | see log | converges |
| `34_dprk_2009.mpi4.log` | see log | converges |
| `35_dprk_2013.mpi4.log` | see log | converges |
| `36_dprk_2016a.mpi4.log` | see log | converges |
| `37_dprk_2016b.mpi4.log` | see log | converges |

See `docs/PARALLEL_KSP.md` for the diagnosis, the PC sweep, and the fix.
