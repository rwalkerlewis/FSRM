# Axis-1a Fidelity Report (passes 5-11)

This report is the canonical "where does axis-1a stand" reference for
the historic-nuclear track. It summarises the full axis-1a fidelity
ladder, every gate result across passes 5-11, every named residual,
and the priority queue of follow-up axes.

After pass-11 (axis-1c implicit time stepping for the diffusion solve
+ outer-BC sponge layer + extended ANEOS validation), axis-1a has
three more gates closed at spec target. Remaining axis-1 work
(axis-1b 3D source ball, axis-1d 3D far-field coupling, axis-3
layered-medium fidelity, axis-4 Tillotson refit) is named here
explicitly as the lower-priority backlog. Pass-12 rotates to
axis-1b (3D source ball) per the pass-11 spec.

## Fidelity ladder (closed at every tier after pass-11)

| Ladder                       | LOW                | MED            | HIGH                              | HIGHEST              |
| ---------------------------- | ------------------ | -------------- | --------------------------------- | -------------------- |
| `radiation_phase`            | `ZELDOVICH_RAIZER` | `MARSHAK_GREY` | `MARSHAK_GREY` + `TABULATED_PATCHED` | `MARSHAK_MULTIGROUP` |
| `cavity_eos`                 | `IDEAL_GAS`        | `TILLOTSON`    | `TILLOTSON_TABULATED_PATCH`       | `TABULATED_FULL`     |
| `opacity_model`              | `CONSTANT`         | `POWER_LAW_ZR` | `TABULATED_PATCHED`               | `TABULATED_FULL`     |
| `operator_splitting`         | `LIE`              | `LIE`          | `STRANG`                          | `STRANG_MULTIGROUP`  |
| `time_integrator`            | `EXPLICIT_EULER`   | `EXPLICIT_EULER` | `EXPLICIT_EULER`                | `RK3_SSP`            |
| `time_integrator_diffusion`  | `BACKWARD_EULER`   | `BACKWARD_EULER` | `CRANK_NICOLSON`                | `BDF2`               |

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

| Gate                              | Pass-7 | Pass-8 | Pass-9 | Pass-10 final | Pass-11 final | Spec target | Residual / next        |
| --------------------------------- | ------ | ------ | ------ | ------------- | ------------- | ----------- | ---------------------- |
| Salmon CavityRadius               | factor 5 | factor 4 | factor 3 | factor 3 retained | factor 3 retained | 5%       | axis-1b 3D source ball |
| Salmon FreeFieldPeakVelocity_166m | factor 5 | factor 4 | factor 2 (active) | factor 2 retained | factor 2 retained | factor 2 | -- (closed)        |
| Salmon FreeFieldPeakVelocity_322m | factor 6 | factor 4 | factor 2 (active) | factor 2 retained | factor 2 retained | factor 2 | -- (closed)        |
| Salmon FreeFieldPeakVelocity_549m | skipped | skipped | skipped | factor 4 (named pass-10 work) | sponge layer available; gate retained at factor 4 with corrected diagnosis | factor 2 | impedance BC contributes; sponge BC ships as opt-in. Closure requires either tighter spec sweep or axis-1d 3D far-field coupling |
| Salmon FarFieldBodyWaveMagnitude  | +/- 0.5 | +/- 0.3 | +/- 0.3 | +/- 0.3 retained | +/- 0.3 retained | +/- 0.2 | propagation-path drift, axis-3 layered-medium fidelity |
| Marshak SelfSimilarPureRadiation  | n/a | factor 10 | factor 8 | factor 2.5 (multigroup) | **factor 2** (BDF2) | factor 2 | -- (closed; HIGHEST tier under axis-1c) |
| Marshak RadiationEnergyConservation | n/a | 25% | 10% | 10% retained | **2%** (BDF2) | 2% | -- (closed; HIGHEST tier under axis-1c) |
| Marshak GreyVsZRComparison         | n/a | factor 5 | factor 3 | factor 3 retained | factor 3 retained | factor 3 | -- (closed)             |
| Chagan CavityRadius               | factor 6 | factor 5 | factor 5 | factor 5 retained | factor 5 retained | factor 3 | axis-1b 3D source ball |
| Chagan FarFieldBodyWaveMagnitude  | +/- 0.5 | +/- 0.4 | +/- 0.4 | +/- 0.4 retained | +/- 0.4 retained | +/- 0.3 | propagation-path drift, axis-3 |
| PokhranI FarFieldBodyWaveMagnitude | +/- 0.5 | +/- 0.4 | +/- 0.4 | +/- 0.4 retained | +/- 0.4 retained | +/- 0.3 | regional Murphy 1981 reference-value drift, axis-4 regional refit |
| Granite Hugoniot match            | n/a | n/a | n/a | n/a | **30%, half points within** | 5% | Tillotson parameter refit, axis-4 |
| Salt Hugoniot match               | n/a | n/a | n/a | n/a | **30%, half points within** | 5% | Tillotson parameter refit, axis-4 |

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

## What pass-11 closed

- **Axis-1c (implicit diffusion time stepping).** New
  DiffusionTimeIntegrator strategy with three concrete subclasses
  (BackwardEuler, CrankNicolson, BDF2). Selected via the new
  `time_integrator_diffusion` config knob; BDF2 is the HIGHEST tier
  default and BACKWARD_EULER preserves pass-10 byte-identical
  behaviour. Closes the Marshak self-similar gate from factor 2.5
  to factor 2 (HIGHEST tier under BDF2) and the radiation-energy
  conservation gate from 10% to 2%.
- **549 m free-field BC triage.** The pass-10 finding ("549 m sits
  at factor 4 vs spec target factor 2") was triaged via a
  radial_outer_radius_m sweep at [700, 1000, 1500, 2000] m. Sweep
  evidence: peak velocity at 549 m drops by factor ~3.4 from
  r_outer = 700 m to 1000 m, indicating the impedance BC at 700 m
  contaminates the gauge. Pass-11 ships the Israeli & Orszag 1981
  graded-damping sponge layer (gated by sponge_layer_enabled) as
  the BC fix.
- **Extended ANEOS Hugoniot validation.** New gates verify the
  granite (Marsh 1980) and salt (McQueen 1970) tabulated EOS
  reproduces published shock-Hugoniot data within 30% at half or
  more sample points. The 5% spec target requires a Tillotson
  parameter refit (named axis-4 follow-up).
- **Axis-1b 3D source ball scaffolding.** Source3DBall.hpp interface,
  cavity_geometry config knob, makeSource3DBall factory, dispatch
  wiring (THREE_DIMENSIONAL throws on selection with a clear
  pass-12 / axis-1b diagnostic), and docs/AXIS_1B_DESIGN.md design
  stub. Pass-12 lands the implementation on this contract.

## What pass-11 did NOT close (named follow-up)

- **Salmon / Chagan CavityRadius (factor 3 / 5 envelopes).** Closure
  named **axis-1b** (3D source ball, scaffolded in pass-11; pass-12
  implementation).
- **Salmon FreeFieldPeakVelocity_549m (factor 4 envelope).** Pass-11
  ships the sponge BC as the fix. Closing the gate to factor 2
  requires either a stricter spec sweep run with the sponge enabled
  in HIGHEST tier or axis-1d 3D far-field coupling.
- **Salmon / Chagan / PokhranI FarFieldBodyWaveMagnitude.** Same
  axis-3 (layered-medium fidelity) follow-up as pass-10.
- **Granite / Salt Hugoniot 5% target.** The pass-9 / pass-10
  Tillotson parameter sets fit individual Hugoniot points to ~30%.
  Closing to 5% requires a Tillotson refit named **axis-4**.

## Next-pass priorities (pass-12 and beyond)

- **Pass-12: axis-1b (3D source ball).** Implement the concrete
  Source3DBall subclass on the pass-11 scaffold. Unstructured-tet
  mesh, 3D Drucker-Prager, asymmetric overburden, surface-integral
  moment-tensor extraction. Closes Salmon and Chagan CavityRadius
  gates. Design stub: docs/AXIS_1B_DESIGN.md.
- **Axis-1d: 3D far-field wavefield extraction.** Replace the
  surface-integral moment-tensor extraction at the elastic radius
  with a coupled radial-source -> 3D FEM coupling. Pass-12 axis-1b
  feeds the asymmetric moment tensor that axis-1d consumes.
- **Axis-2: topography.** Free-surface topography for the mountain
  mesa, basin, and valley shots in the historic catalogue (Pahute
  Mesa, Chagan, Faultless).
- **Axis-3: layered-medium fidelity.** Higher-resolution velocity
  model, Q-attenuation as a function of depth, regional refit of
  magnitude reference values. Closes the body-wave magnitude gates.
- **Axis-4: Tillotson parameter refit.** Per-medium least-squares
  refit of the Tillotson EOS constants (a, b, A, B, alpha, beta,
  E_0, E_iv, E_cv) against published Hugoniot data. Closes the 5%
  Hugoniot match target.

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
