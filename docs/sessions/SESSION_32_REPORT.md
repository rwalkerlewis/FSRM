# Session 32 Report -- Multi-cell moment-tensor source distribution (pass 4)

## Summary

Pass-4 closes the pass-3 single-cell moment-tensor limitation that
`docs/HISTORIC_NUCLEAR_FIDELITY.md` section 2 identified as the
structural reason the historic-nuclear amplitude envelope could not
shrink below factor 100. The new `[SOURCE_DISTRIBUTION]` config
grammar exposes three injection modes (`SINGLE_CELL` legacy default,
`GAUSSIAN`, `UNIFORM_SPHERE`); the distributed modes spread the
moment tensor over a spherically symmetric support ball with weights
that sum to `M0` globally. The five anchor historic tests gain
`_Distributed` variants asserting at the tightened factor-30
envelope. Baneberry 1970 (previously blocked at factor 100) is
unblocked at factor 30; Sterling 1966 (previously blocked at factor
100) is unblocked at the legacy factor-100 envelope; DPRK 2006
remains blocked because Rc ~ 10 m vs cell-corner-to-centroid scale
~707 m on the 4x4x4 CI mesh.

Headline numbers: full-suite 147 / 154 passing (up from 135 / 142
pre-pass-4 baseline computed by subtracting the 12 new tests
registered in pass 4); the 7 failures are the same 6 pre-existing
fault-solver tests plus the pre-existing
`Unit.MuellerMurphyMomentConsistency.BodyWaveMagnitudeInSedanRange`
that pass-3 already documented as a stale-bound test on `main`. No
regressions.

Branch verdict: **A** (multi-cell distribution works for 6 / 8
documented blocked-or-loose tests, tightens the envelope by 3.3x
for the anchors, and the byte-identical SINGLE_CELL guarantee plus
M0 conservation tests verify backward compatibility).

## Code changes

### Config grammar and runtime dispatch

| File:line | Change |
|---|---|
| `src/core/Simulator.cpp:412-419` | New `DistributedSourceCell` struct in the anonymous namespace; cell index plus pre-computed `weight_volume` (w_c * V_c, normalized) |
| `src/core/Simulator.cpp:451-475` | New fields on `Simulator::ExplosionCoupling`: `source_distribution_mode`, `source_support_radius_factor`, `source_gaussian_sigma_factor`, `source_min_cells`, plus the `dist_cells` cache and DM-pointer guard |
| `src/core/Simulator.cpp:1247-1295` | `[SOURCE_DISTRIBUTION]` parsing in `Simulator::initializeFromConfigFile`; case-insensitive mode validation with single rank-0 warning on unrecognised values |
| `src/core/Simulator.cpp:6709-6794` | `injectMomentTensorIntoCell` static helper -- factored from the original single-cell injection so SINGLE_CELL and distributed modes share the per-cell nodal-force formula |
| `src/core/Simulator.cpp:6796-6905` | `buildDistributedSourceCellsImpl` static helper -- enumerates cells inside the support ball, computes per-cell weight (Gaussian or uniform), MPI_Allreduces the normalization sum so `sum_c w_c V_c = 1` globally; falls back to SINGLE_CELL if `< min_cells` enumerate |
| `src/core/Simulator.cpp:6907-7010` | `addExplosionSourceToResidual` rewritten as a dispatch on `source_distribution_mode`; SINGLE_CELL path is byte-identical to the pre-pass-4 code |

### Tests

| File:line | Change |
|---|---|
| `tests/integration/test_source_distribution.cpp` | New file (387 lines); five tests: `SingleCellLegacyByteIdentical`, `GaussianM0Conserved`, `UniformSphereM0Conserved`, `GaussianBeatsSingleCellOnFineMesh`, `FallbackToSingleCellWhenBallEmpty` |
| `tests/integration/test_historic_nuclear.cpp:88-99` | `writeConfig` gains a `dist_section` parameter (default empty preserves byte-identical config); appended verbatim if non-empty |
| `tests/integration/test_historic_nuclear.cpp:184-187` | `dist_section` insertion in the EXPLOSION_SOURCE block |
| `tests/integration/test_historic_nuclear.cpp:443-457` | `assertFarFieldAndPolarity` gains an `amplitude_envelope_factor` parameter (default 100, distributed variants pass 30 or 100) |
| `tests/integration/test_historic_nuclear.cpp:1239-1320` | `FarFieldAmplitudeRegression` rewritten to record both SINGLE_CELL and UNIFORM_SPHERE_50x rows in `historic_nuclear_regression.csv`; CSV schema gains a `source_distribution_mode` and `ratio` column |
| `tests/integration/test_historic_nuclear.cpp:1325-1545` | New `_Distributed` variants for the five anchor tests, Baneberry 1970, Sterling 1966, and DPRK 2006 |
| `tests/CMakeLists.txt:138` | `integration/test_source_distribution.cpp` added to `INTEGRATION_TEST_SOURCES` |
| `tests/CMakeLists.txt:351-364` | New `add_test` registrations for the seven distributed historic variants (DPRK 2006 left unregistered with documented note) |
| `tests/CMakeLists.txt:357-361, 477-481, 503-507` | New `add_test` registrations for the five `Integration.SourceDistribution.*` tests; all entries added to the integration timeout/label group |

### Documentation

| File:line | Change |
|---|---|
| `docs/HISTORIC_NUCLEAR_FIDELITY.md` | Pass-4 entry added to the introduction; section 1 ("What is verified") gains 9 new rows; section 2 ("What is NOT verified") "Single-cell moment-tensor approximation" replaced by the narrower DPRK 2006 limitation; section 3 ("Deviations") records the per-test envelope choice; new section 4c "Closed in pass 4"; pass-3 "remain open" item moved into 4c |
| `docs/CONFIGURATION.md` | New `[SOURCE_DISTRIBUTION]` section between Initial Conditions and Solver Settings, describing all four keys with defaults and a worked example |
| `README.md:67` | Test count updated from 116 to 154 |
| `README.md:101-110` | New rows for the seven distributed historic variants and the five `Integration.SourceDistribution.*` tests |

## Measurement results

### Anchor-test peak/u_far ratio: SINGLE_CELL vs UNIFORM_SPHERE_50x

The numbers below are the worst-case ratio across each test's BHZ
peak amplitude divided by the Aki & Richards far-field analytic
estimate `u_far ~ 2 * Mdot_peak / (4 pi rho vp^3 R)`. Single-cell
ratios are upper-only (peak too high); distributed ratios swing the
peak below the analytic estimate (peak too low) because the moment
is smeared over a region wider than Rc on the coarse CI mesh. Both
inflation directions are captured by the `[1/N, N]` envelope.

| Test | Single-cell ratio | Distributed (50x) ratio | Distributed envelope |
|---|---:|---:|:---:|
| Gasbuggy 1967 | 9 | < 30 | factor 30 |
| Gnome 1961 | 78 | < 30 | factor 30 |
| Sedan 1962 | 9 | 0.31 | factor 30 |
| Degelen Mountain | < 30 | < 30 | factor 30 |
| NTS Pahute Mesa | 58 | < 30 | factor 30 |
| Baneberry 1970 | 249 (blocked) | < 30 | factor 30 |
| Sterling 1966 | 141 (blocked) | 96 (UNIFORM_SPHERE_100x) | factor 100 |
| DPRK 2006 | 290 (blocked) | 172 (UNIFORM_SPHERE_100x) | still blocked |

The `historic_nuclear_regression.csv` records both rows of the
Sedan 1962 measurement automatically; the table above for the other
tests was reconstructed from the per-test CTest stdout.

### Test count

| Subset | Pre-pass-4 baseline | Post-pass-4 | Delta |
|---|---:|---:|---:|
| Total registered tests | 142 | 154 | +12 |
| Passing | 135 | 147 | +12 |
| Failing (pre-existing) | 7 | 7 | 0 |
| Failing (regressions) | 0 | 0 | 0 |

The 7 pre-existing failures: `Physics.SCEC.TPV5`,
`Physics.LockedFaultTransparency`,
`Integration.PressurizedFractureFEM`,
`Integration.DynamicRuptureSolve.PrescribedSlip`,
`Integration.SlippingFaultSolve`,
`Integration.SlipWeakeningFault`, and
`Unit.MuellerMurphyMomentConsistency.BodyWaveMagnitudeInSedanRange`.
The first six are documented in `docs/SOLVER_STATE.md` "Failing
fault tests"; the last is documented in
`docs/HISTORIC_NUCLEAR_FIDELITY.md` section 4b.

## Fault subset table

Unchanged. Pass-4 does not touch the fault-solver code path
(`src/physics/`, `src/numerics/`, the cohesive Jacobian assembly in
`src/core/Simulator.cpp`, or the rim-pin / aux-slip infrastructure).
The default fault subset remains 10 / 16 (Session 30 verified) and
the experimental fieldsplit path remains 7 / 16. No change to
`docs/SOLVER_STATE.md` "Current test state" table.

## Branch verdict

**Branch A** (fix works, no regression).

The multi-cell moment-tensor distribution closes the pass-3 envelope
limitation for 6 of 8 documented blocked-or-loose historic tests
(the five anchors at factor 30 plus Baneberry 1970 at factor 30 plus
Sterling 1966 at factor 100). DPRK 2006 remains blocked due to a
mesh-resolution constraint (4x4x4 CI cells cannot resolve Rc ~ 10
m); this is documented honestly as the residual case in section 4c
of `HISTORIC_NUCLEAR_FIDELITY.md`. The byte-identical SINGLE_CELL
test (`Integration.SourceDistribution.SingleCellLegacyByteIdentical`)
guarantees backward compatibility for every config that omits the
new section.

## What this changes in `docs/SOLVER_STATE.md`

No change to the fault-solver state table. Adding a single change-log
entry below the existing pass-3 line:

```
- pass-4: pass-4 multi-cell moment-tensor source distribution
  (Integration.SourceDistribution.*, 5 tests; historic-nuclear
  _Distributed variants, 7 tests). Does not modify Jacobian /
  cohesive code paths; default fault subset count unchanged at 10/16.
```

## Session 33 recommendation

Two next-cycle candidates surfaced during pass 4:

1. **Per-test mesh refinement for DPRK 2006**, or equivalent
   1-D-source embedding. The 4x4x4 CI mesh cannot resolve Rc ~ 10 m
   on a 0.7 kt granite shot, but extending the existing
   `[MESH_REFINEMENT]` plumbing to refine to `refinement_levels = 2`
   with `[SOURCE_DISTRIBUTION] mode = GAUSSIAN, support_radius =
   100 * Rc` may close the residual ratio. CI time-budget impact
   needs measurement first.
2. **AK135 1-D Earth ray tracing for the DPRK 2017 synthetic-mb
   test.** Pass-3 documented `t_star_ref = 5.5 s` overstates the
   literature t* by ~7x as a fudge factor compensating for un-modeled
   radiation pattern, free-surface doubling, and station correction.
   Closing this requires routing the Mueller-Murphy moment rate
   through an AK135 source-receiver path with proper P/SV
   decomposition and station correction.

Pick one. Candidate 1 has a smaller scope and reuses pass-3 and
pass-4 infrastructure; candidate 2 is a bigger lift but closes a
deeper documented gap.
