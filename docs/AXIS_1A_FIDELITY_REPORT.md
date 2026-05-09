# Axis-1a Fidelity Report (passes 5-10)

This report is the canonical "where does axis-1a stand" reference for
the historic-nuclear track. It summarises the full axis-1a fidelity
ladder, every gate result across passes 5-10, every named residual,
and the priority queue of follow-up axes.

After pass-10, axis-1a is closed for now; pass-11 rotates to topography
(axis 2). Remaining axis-1 work (axis-1b 3D source ball, S_n transport,
axis-1c implicit-time-stepping for the diffusion solve) is named here
explicitly as the lower-priority backlog.

## Fidelity ladder (closed at every tier after pass-10)

| Ladder              | LOW                | MED            | HIGH                              | HIGHEST              |
| ------------------- | ------------------ | -------------- | --------------------------------- | -------------------- |
| `radiation_phase`   | `ZELDOVICH_RAIZER` | `MARSHAK_GREY` | `MARSHAK_GREY` + `TABULATED_PATCHED` | `MARSHAK_MULTIGROUP` |
| `cavity_eos`        | `IDEAL_GAS`        | `TILLOTSON`    | `TILLOTSON_TABULATED_PATCH`       | `TABULATED_FULL`     |
| `opacity_model`     | `CONSTANT`         | `POWER_LAW_ZR` | `TABULATED_PATCHED`               | `TABULATED_FULL`     |
| `operator_splitting`| `LIE`              | `LIE`          | `STRANG`                          | `STRANG_MULTIGROUP`  |
| `time_integrator`   | `EXPLICIT_EULER`   | `EXPLICIT_EULER` | `EXPLICIT_EULER`                | `RK3_SSP`            |

Every cell has a working implementation as of pass-10. LOW / MED / HIGH
defaults reproduce the corresponding pass-N behaviour byte-for-byte
under the pinned configs of historic-nuclear tests. HIGHEST tier is
opt-in via the new config sub-keys.

## Gate results across passes (pass-9 / pass-10 deltas highlighted)

Each row lists the achieved envelope for the named gate at that pass.
"Spec target" is the strict-spec envelope from the original
HISTORIC_NUCLEAR_ROADMAP.md fidelity table. "Pass-10 final" is the
envelope shipped on this branch. "Residual / next" names the
follow-up axis when the strict spec is not yet reached.

| Gate                              | Pass-7 | Pass-8 | Pass-9 | Pass-10 final | Spec target | Residual / next        |
| --------------------------------- | ------ | ------ | ------ | ------------- | ----------- | ---------------------- |
| Salmon CavityRadius               | factor 5 | factor 4 | factor 3 | factor 3 retained | 5%       | axis-1b 3D source ball |
| Salmon FreeFieldPeakVelocity_166m | factor 5 | factor 4 | factor 2 (active) | factor 2 retained | factor 2 | -- (closed)        |
| Salmon FreeFieldPeakVelocity_322m | factor 6 | factor 4 | factor 2 (active) | factor 2 retained | factor 2 | -- (closed)        |
| Salmon FreeFieldPeakVelocity_549m | skipped | skipped | skipped | factor 4 (named pass-10 work) | factor 2 | RK3_SSP available; gate enabled at factor 4 envelope |
| Salmon FarFieldBodyWaveMagnitude  | +/- 0.5 | +/- 0.3 | +/- 0.3 | +/- 0.3 retained | +/- 0.2 | propagation-path drift, axis-3 layered-medium fidelity |
| Marshak SelfSimilarPureRadiation  | n/a | factor 10 | factor 8 | factor 2.5 (multigroup) | factor 2 | very-early-time BE smearing residual; axis-1c CN/BDF2 |
| Marshak RadiationEnergyConservation | n/a | 25% | 10% | 10% retained | 2% | implicit-Euler matter coupling residual; axis-1c |
| Marshak GreyVsZRComparison         | n/a | factor 5 | factor 3 | factor 3 retained | factor 3 | -- (closed)             |
| Chagan CavityRadius               | factor 6 | factor 5 | factor 5 | factor 5 retained | factor 3 | axis-1b 3D source ball |
| Chagan FarFieldBodyWaveMagnitude  | +/- 0.5 | +/- 0.4 | +/- 0.4 | +/- 0.4 retained | +/- 0.3 | propagation-path drift, axis-3 |
| PokhranI FarFieldBodyWaveMagnitude | +/- 0.5 | +/- 0.4 | +/- 0.4 | +/- 0.4 retained | +/- 0.3 | regional Murphy 1981 reference-value drift, axis-4 regional refit |

## What pass-10 closed

- **MARSHAK_MULTIGROUP transport** (HIGHEST tier on the radiation ladder).
  Per-group analytic opacity model (path A: Mihalas-Mihalas 1984
  sec 82.2 smoothed-continuum bound-bound + Kramers free-free + Thomson
  scattering). Six new physics-validation gates: PlanckIntegralSumsToSigmaT4,
  GroupOpacityAnalyticPathSanity, SelfSimilarPureRadiation,
  RadiationFrontPositionConvergesWithG, GroupSumMatchesGrey,
  GreyVsMultigroupComparison.
- **TABULATED_FULL EOS** (HIGHEST tier on the EOS ladder). Pure-tabulated
  evaluation everywhere with a Tillotson safety net for out-of-table
  queries. The pass-9 sin^2 patch window is removed under this mode.
- **TABULATED_FULL opacity** (HIGHEST tier on the opacity ladder). Pure
  tabulated grey opacity (no Z-R fallback). The dispatch matrix
  `(opacity_model, radiation_phase)` is documented and enforced:
  TABULATED_PATCHED + MARSHAK_MULTIGROUP throws a clear error directing
  the user at TABULATED_FULL or POWER_LAW_ZR.
- **TVD_RK2 (Heun's method) and RK3_SSP (Shu-Osher 1988)** time
  integrators. EXPLICIT_EULER preserved as the pass-9 byte-identical
  default. RK3_SSP unblocks the 549 m free-field gauge by relaxing the
  CFL safety cap that pass-9 documented.
- **Strang inner-substep convergence diagnostic.** New
  operator_splitting_convergence_diagnostic flag exposes the
  inner-CFL substep dt for the OperatorSplittingConvergence test.

## What pass-10 did NOT close (named follow-up)

The pass-10 spec set strict-spec targets for every gate. The targets
not reached under HIGHEST tier on this pass are tracked here with
honest naming of the follow-up axis.

- **Salmon CavityRadius (5% spec target).** Pass-10 retains the factor-3
  envelope. The axis-1a path A (1D radial Lagrangian + multigroup +
  Tillotson tabulated patch) is fundamentally limited by the
  spherical-symmetry assumption: the actual NTS cavity is asymmetric
  under overburden. Closure named **axis-1b** (3D source ball).
- **Marshak SelfSimilarPureRadiation (factor 2 spec target).**
  Pass-10 ships factor 2.5 at late times. The very-early-time
  backward-Euler smearing of the diffusion front is the residual.
  Closure named **axis-1c** (Crank-Nicolson or BDF2 time stepping for
  the per-group diffusion solve; structure-preserving alternative to
  the backward-Euler solve in MarshakRadiationDiffusion.cpp).
- **Marshak RadiationEnergyConservation (2% spec target).** Pass-10
  retains 10%. Same axis-1c follow-up as above.
- **Salmon / Chagan / PokhranI FarFieldBodyWaveMagnitude.** These are
  propagation-path issues, not source-physics. Closure named
  **axis-3** (layered-medium velocity model fidelity; pass-10 uses
  the configured 1D layered medium; pass-9 documented this as
  needing regional refit of the Murphy 1981 reference values for
  PokhranI specifically).
- **549 m FreeFieldPeakVelocity (factor 2 spec target).** Pass-10
  enables this gate at the factor-4 envelope under HIGHEST tier
  (RK3_SSP + radial_outer_radius_m = 700). Closure to factor 2
  named **axis-1d** (3D wavefield extraction from a coupled radial
  source -> 3D far-field FEM mesh; pass-10's surface-integral
  moment-tensor extraction at the elastic radius is the source-side
  bottleneck).

## Next-pass priorities (axis 2 and beyond)

- **Pass-11: axis-2 (topography).** The historic-nuclear roadmap calls
  out free-surface topography as the next dominant fidelity axis after
  axis-1a. The mountain mesa, basin, and valley shots in the historic
  catalogue (Pahute Mesa, Chagan, Faultless) should be re-run with a
  topography-resolved upper boundary.
- **Axis-1b: 3D source ball.** Asymmetric overburden, layered medium
  in the source-region 3D solve. Removes the spherical-symmetry
  assumption in axis-1a. Required to close Salmon CavityRadius and
  Chagan CavityRadius to spec.
- **Axis-1c: implicit time stepping (CN / BDF2).** Crank-Nicolson or
  BDF2 second-order time integration for the (multigroup) diffusion
  solve. Closes the very-early-time BE smearing residual in the
  Marshak gates.
- **Axis-1d: 3D far-field wavefield extraction.** Replace the
  surface-integral moment-tensor extraction at the elastic radius
  with a coupled radial-source -> 3D FEM coupling that lets the
  far-field mesh resolve the wavefield on its own scale. Closes the
  549 m free-field gauge to factor 2.
- **Axis-3: layered-medium fidelity.** Higher-resolution velocity model,
  Q-attenuation as a function of depth, regional refit of magnitude
  reference values. Closes the body-wave magnitude gates.
- **Axis-4 (named only): full ANEOS table coverage** for tuff /
  alluvium media at coverage extending to 1e7 Pa. Pass-10 ships
  pass-9-level coverage; the spec called for extension to 1e8 Pa
  but the underlying ANEOS rerun is a separate body of work.

## References

- Mihalas, D. and Mihalas, B. W. (1984), "Foundations of Radiation
  Hydrodynamics", Oxford University Press. Sec 80 (multigroup
  diffusion); sec 82.2 (line opacity smoothed continuum); sec 96-97
  (operator splitting between matter and radiation).
- Pomraning, G. C. (1973), "The Equations of Radiation Hydrodynamics",
  Pergamon Press. Ch IV (multigroup formulation).
- Zel'dovich, Y. B. and Raizer, Y. P. (1967), "Physics of Shock Waves
  and High-Temperature Hydrodynamic Phenomena", vol I, Academic Press.
  Ch V (free-free Kramers' opacity, frequency dependence).
- Shu, C.-W. and Osher, S. (1988), "Efficient implementation of
  essentially non-oscillatory shock-capturing schemes", J. Comp. Phys
  77, pp 439-471. Strong-stability-preserving Runge-Kutta family.
- Strang, G. (1968), "On the construction and comparison of difference
  schemes", SIAM J. Num. Anal. 5(3), pp 506-517. Operator splitting.
- Larsen, E. W. (1988), "A grey transport acceleration method for
  time-dependent radiative transfer", J. Comp. Phys 78, pp 459-480.
  Linearised Newton on T^4 closure (carried through to pass-10
  multigroup).
- Wilkins, M. L. (1980), "Computer Simulation of Dynamic Phenomena",
  Springer. Artificial viscosity, hourglassing, radial return
  formulation (carried over from pass-6).
- Murphy, J. R. (1981), "Magnitude-yield relations for underground
  nuclear explosions". Reference values for the body-wave magnitude
  gates (axis-3 follow-up named).
