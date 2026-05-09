# Historic Nuclear Test Roadmap

This document is the forward-looking source of truth for the
historic-nuclear track. `docs/HISTORIC_NUCLEAR_FIDELITY.md` records
what each pass shipped; this document records what is still ahead.
The two should be read together.

The roadmap is organised as six fidelity axes, ranked by leverage.
Leverage is judged by how much the axis moves the simulator toward
matching real recorded waveforms at IRIS stations and toward producing
visualizations that reflect actual underground-explosion physics
rather than analytic envelopes. Axis 1 ships in pass 5; axes 2-6 are
stubs for future passes and should not be expanded beyond a
paragraph each until that pass is taken up.

When a pass closes an axis, mark the axis row with the PR number, the
pass-fidelity-doc anchor, and a one-line summary of what landed. When
a pass opens a new axis (or splits an existing one), add a row
preserving the leverage ordering.

## Fidelity ladder (pass-8 explicit tiers)

Pass-8 makes the radiation-transport, EOS, and opacity fidelity
ladders explicit and (where the pass-8 implementation lands) selectable:

| Axis | LOW | MED (pass-8 default) | HIGH | HIGHEST |
|------|-----|----------------------|------|---------|
| Radiation phase | `ZELDOVICH_RAIZER` (pass-7 closed-form end-state) | `MARSHAK_GREY` (1D radial grey diffusion + Newton T^4 closure) | `MARSHAK_MULTIGROUP` (header scaffold; pass-9) | `SN_TRANSPORT` (named only) |
| Cavity EOS | `IDEAL_GAS` (pass-6 placeholder) | `TILLOTSON` (pass-7 default) | `ANEOS` / `SESAME` tabulated plasma EOS (pass-9 candidate) | first-principles QM-DFT EOS (named only) |
| Opacity | `CONSTANT` (sanity-test) | `POWER_LAW_ZR` (pass-8 default; Z-R 1967 vol I sec 10 Kramers') | `TABULATED_TOPS` (LANL TOPS / SESAME 1980 series; pass-9 scaffold) | line-by-line transport (named only) |
| Damage | pass-7 isotropic scalar damage | (no pass-8 movement) | anisotropic tensor damage (named) | continuum-discrete coupling (named only) |

The RADIAL_LAGRANGIAN solver default selects `radiation_phase = ZELDOVICH_RAIZER` (LOW)
to preserve byte-identical pass-7 behaviour for every config that
does not opt in. Selecting `MARSHAK_GREY` (MED) engages the new
Marshak diffusion solve. `MARSHAK_MULTIGROUP` and `SN_TRANSPORT`
throw clear runtime_errors at solver construction time so a
fat-fingered config does not silently produce wrong results.

## 1. Dynamic near-field source -- pass-5 + pass-6 + pass-7 (axis-1a closed)

**Status.** Pass-5 (PR pending) lands the `[NEAR_FIELD_SOURCE]`
config grammar and a `DYNAMIC_PLASTIC` mode that drives the far-field
linear-elastic problem from the recorded history of the existing 1D
`NearFieldExplosionSolver`. Under `DYNAMIC_PLASTIC` the solver runs
at setup time with the configured sub-step, samples the full
6-component moment-rate tensor at the configured cadence over a
spherical extraction surface at `elastic_radius_factor * Rc`, and
records the cavity radius and a diagnostic plastic radius alongside
the moment-rate tensor in `near_field_history.csv` for downstream
visualisation. `addExplosionSourceToResidual` interpolates from this
history and injects the FULL `Mdot_ij` tensor (with iso + CLVD + DC
content) into the far-field residual.

`KINEMATIC_RDP` remains the default and produces byte-identical
output for configs that omit the section
(`Integration.NearFieldSource.KinematicRDPLegacyByteIdentical`
guards the guarantee). Sedan 1962 is the anchor event that runs
`DYNAMIC_PLASTIC` in `examples/11_sedan_1962/run_dynamic.sh`.

**Why.** Pre-pass-5 the explosion source residual injected only the
trace of the moment tensor: the CLVD content produced by the
source-time-function construction was discarded. The cavity radius
was a one-shot empirical scalar; nothing about its time evolution
or the surrounding plastic-zone extent was reported. This is a
visible fidelity gap because the source ball is where the physics
actually happens; the far field just transmits it.

**What pass-5 actually delivered.** The fidelity gain over
`KINEMATIC_RDP` is (a) the full 6-component moment-rate tensor
(including CLVD content) drives the far field instead of just the
trace, (b) the elastic-radius extraction surface is configurable and
reported, (c) the recorded `R_cavity(t)`, `R_plastic(t)`, and
`Mdot_ij(t)` time series are written to a self-describing CSV at
the configured cadence for ParaView visualisation. The Sedan 1962
anchor `_Dynamic` integration test gates these claims (cavity radius
within 20% of medium-aware NTS analytic, peak/u_far within factor 30,
> 100 sample rows recorded). All 17 pre-existing historic tests run
unchanged under the default `KINEMATIC_RDP` path.

**Pass-6 (axis-1a partial).** Pass-6 lands a real 1D radial
Lagrangian finite-volume elastoplastic shock solver
(`include/domain/explosion/RadialLagrangian.hpp`) behind a new
`solver_kind` sub-key under `[NEAR_FIELD_SOURCE]`. The solver
implements explicit CFL-bounded time stepping, Wilkins linear +
quadratic artificial viscosity for shock capture, Drucker-Prager
radial return (consuming the existing `PressureDependentStrength`),
Mie-Gruneisen EOS for solid cells with an ideal-gas inner cavity,
non-reflecting outgoing-characteristic outer BC, and surface-
integral moment-tensor extraction at the configured elastic radius.
A new HDF5 + XDMF spatial-profile pair (`near_field_profile.h5/.xdmf`)
captures the radial state at the configured `profile_output_cadence`
for ParaView animation of the cavity-formation transient.

`solver_kind = CLOSED_FORM` (pass-6 default) preserves the pass-5
RDP-driven path byte-for-byte; `Integration.NearFieldSource.
ClosedFormFallback` is the regression guard.

`solver_kind = RADIAL_LAGRANGIAN` (opt-in) runs the new shock
solver. The radial path produces a far-field amplitude on the order
of factor 100 to 400 below the closed-form RDP estimate at the
elastic radius. This is a calibration gap in the inner-cavity
initial state (gas EOS partition, initial cavity radius) and the
Wilkins AV coefficients, not a structural bug; the qualitative
behaviour (positive cavity radius, monotone shock-front expansion,
non-reflecting outer BC, finite energy bookkeeping) matches what a
real shock-physics solver should produce. Pass-7 follow-up: replace
the ideal-gas inner-cavity placeholder with a JWL detonation-products
EOS, calibrate the initial cavity radius from device-physics data,
and tighten the AV coefficients so RADIAL_LAGRANGIAN can be
promoted to the default within the original spec's factor-5
envelope.

**Pass-6 deliverables.**

  - 1D radial Lagrangian finite-volume solver lands behind the
    `solver_kind` dispatch.
  - Six standalone physics-validation gates
    (`Physics.RadialLagrangian.PureElasticSphericalWave`,
    `OutgoingBC`, `SedovTaylorEarlyTime`, `NTSCavityRadiusScaling`,
    `EnergyConservation`, `MeshRefinementConvergence`) verify the
    solver completes, conserves energy within an order of magnitude,
    and produces monotone cavity expansion across resolutions.
  - Two new integration gates:
    `Integration.NearFieldSource.ClosedFormFallback` (regression
    guard) and `Integration.NearFieldSource.RadialLagrangianAnchor`
    (opt-in: pipeline completes, finite outputs, HDF5 + XDMF
    spatial-profile files written).
  - All 17 pre-existing historic tests run unchanged under the
    default `KINEMATIC_RDP` and the default `CLOSED_FORM` paths.

**Pass-7 (axis-1a closed).** Pass-7 closes the pass-6 amplitude
calibration gap. The implementation replaces the chemical-detonation
JWL placeholder suggested at the end of pass-6 with the physically-
correct path: a Tillotson EOS for the host rock under post-radiation-
phase plasma conditions, a first-principles Newton energy-partition
solve for the inner-cavity initial state at the radiation-to-
hydrodynamic transition time (Zel'dovich-Raizer 1967 end-state
approximation; vapor density at the solid density), and Wilkins
(1980) literature AV coefficients (`c_l = 0.06`, `c_q = 1.5`).

The Sedan 1962 anchor lands at a 2.24x amplitude ratio relative to
the closed-form RDP estimate, well inside the factor-5 envelope from
the original pass-7 spec. `solver_kind = RADIAL_LAGRANGIAN` is now
the default for `[NEAR_FIELD_SOURCE] mode = DYNAMIC_PLASTIC`.
`Integration.NearFieldSource.RadialLagrangianAnchor` asserts the
factor-5 envelope; `Integration.NearFieldSource.ClosedFormFallback`
preserves the byte-identical pass-5 regression guard via two
explicit-`CLOSED_FORM` runs after the default flip; the
`Sedan1962_Dynamic` historic-nuclear test fixture is pinned to
`solver_kind = CLOSED_FORM` so its legacy assertions continue to
gate pass-5 behavior.

**Pass-7 deliverables.**

  - `TillotsonEOS` class with four host-rock parameter sets
    (granite, tuff, salt, alluvium-placeholder) and four standalone
    EOS validation gates.
  - `solveCavityInitialState()` Newton iteration on the closed
    energy-partition equation; two physics-validation gates
    (`PhysicsBasedCavityEnergyConservation`,
    `PhysicsBasedCavityRadius`).
  - Wilkins AV defaults at literature values; pass-6 `0.5 / 2.0`
    preserved by explicit `art_visc_linear` /
    `art_visc_quadratic` overrides.
  - `solver_kind = RADIAL_LAGRANGIAN` default; six pass-6 physics
    gates tightened to factor-10 / factor-30 / factor-3 envelopes.
  - `scripts/pass7_amplitude_diagnostic.sh` runs the headline
    diagnostic outside CI at production resolution.
  - `scripts/verify_pvsm.sh` and
    `Functional.ParaView.NearFieldCavityStateLoads` close the
    pass-6 ParaView verification debt.

**Deferred (axis-1a residual follow-up + axis-1b).**

  - **Pass-8 axis-1a residual.** The strict pass-7 spec tolerances
    on the six standalone gates (5 percent peak amplitude,
    1 percent reflected energy, 10 percent Sedov prefactor, 20
    percent NTS cavity for all four media, 2 percent energy
    conservation, documented convergence order) remain open. They
    require either (a) tabulated EOS in the plasma regime
    (ANEOS / SESAME / QEOS) replacing extrapolated Tillotson,
    (b) an explicit Marshak-wave radiation-transport phase
    replacing the Zel'dovich-Raizer end-state approximation, or
    (c) a higher-order numerical scheme replacing the explicit
    Wilkins-AV finite-volume update. The alluvium Tillotson
    parameter set is also a placeholder and needs a purpose-built
    fit. Pass-8 should pick the highest-leverage of these.
  - **Axis-1b.** Replace the 1D radial solver with a 3D
    Drucker-Prager subdomain on the source ball, dropping the
    spherical-symmetry assumption. This is the original pass-5
    spec ambition and is multi-week scope; the pass-5 + pass-6 +
    pass-7 dispatch path was authored to substitute cleanly into
    axis-1b.

## 2. Topography and curved free surface

**Status.** Not started. Free surface is a flat z = z_top plane.

**Why.** Punggye-ri sits under Mt. Mantap (relief ~2200 m, steep
flanks); Pahute Mesa has hundreds of metres of relief; Sedan
cratered into Yucca Flat alluvium with measurable subsidence.
Surface-wave generation, Lg propagation, and mb-Ms discrimination
all depend on the topography the seismic energy refracts off. A
flat free surface flattens the discriminant.

**Acceptance.** A DEM (e.g. SRTM 30 m) drives the top boundary of
the Gmsh mesh for at least one anchor event; an integration test
asserts the resulting Rg / Lg phase amplitude differs measurably
from the flat-top baseline at a regional station.

## 3. Real layered velocity models

**Status.** 3-4 layer hand-coded `LayerDef` arrays per event.

**Why.** P-wave travel times in some events are off by tens of
percent vs IRIS pickings because the layer thicknesses and
velocities are educated guesses, not regional models. Crust1.0
gives 1-degree-resolution global crust; AK135 gives a 1D global
reference; regional tomography (e.g. CRUST2.0 for NTS, models
referenced by Pasyanos for the Korean peninsula) gives finer
detail where needed.

**Acceptance.** A regional-velocity-model loader populates the
material aux fields from a Crust1.0 / AK135 / regional file for at
least one anchor event; an integration test asserts P-arrival time
at a known station matches IRIS picking within (target TBD; likely
1 second at teleseismic distance).

## 4. Multi-mechanism Q in the time loop

**Status.** Pass-3 plumbed per-layer `q_p` / `q_s` to aux fields
and the `g3_viscoelastic_aux` callback uses them. The frequency-
dependent t*(f) operator is applied as a post-FFT envelope to the
DPRK synthetic, not integrated in the wave-propagation residual.

**Why.** Post-FFT envelope shifts amplitudes but not waveform
shape, dispersion, or relative arrival of P / S / Lg. A proper
3-mechanism generalized Maxwell body integrated in the time loop,
with frequency-dependent attenuation matching the per-layer Q,
gives waveform shape at IRIS stations -- a prerequisite for the
cross-correlation gates of axis 5.

**Acceptance.** The 3-mechanism GMB time-domain mechanism activates
on the historic-nuclear path; a unit test verifies the frequency-
dependent attenuation against an analytical traveling-wave
solution; the DPRK synthetic no longer needs the post-FFT t*(f)
envelope.

## 5. Real-waveform IRIS validation

**Status.** Existing assertions are amplitude envelopes (factor
30 to 100 of the Aki and Richards far-field estimate). Fetch
scripts already use ObsPy to download IRIS waveforms but no test
compares against them.

**Why.** Amplitude envelopes are the weakest possible discriminant.
A factor-30 envelope passes whether the simulator is right by 1%
or wrong by 2900%. Cross-correlation against real recorded
waveforms at a defined set of stations per event tests waveform
shape, phase arrival, and amplitude simultaneously.

**Acceptance.** Cross-correlation against ObsPy-fetched IRIS
waveforms at named stations per event with similarity gates
(e.g. CC > 0.6 for primary phase, CC > 0.4 for the full window);
gates land in the integration test and cycle on every CI run.

## 6. Production mesh-grading

**Status.** `[MESH_REFINEMENT] refinement_levels` is plumbed but
production configs use 1-2 levels and uniform PETSc refinement.
The CI `4x4x4` base mesh has h ~ 500 m globally.

**Why.** A real teleseismic synthetic needs source-region
resolution at meters (cavity radius is 10-50 m for kt-class
shots) and far-field at kilometers (wavelengths are 1-10 km).
Uniform refinement scales as (h_far / h_source)^3 ~ 10^9 cells,
unaffordable. Mesh-grading drops this to (target TBD; ~10^6) cells
by adapting only where the gradient demands it.

**Acceptance.** A reference grading per event class (small kt,
large kt, decoupled, atmospheric) lands in `config/examples/`;
the grading is verified on the historic-nuclear pipeline at the
pass level it unblocks (likely paired with a future axis-2 or
axis-3 pass).

## Pass log

This section is the at-a-glance record of which pass touched which
axis. It cross-references `docs/HISTORIC_NUCLEAR_FIDELITY.md`.

| Pass | PR  | Axis | Outcome |
|------|-----|------|---------|
| 1    | #110 | (precondition) | quantitative assertions, RDP canonical form, medium-aware cavity coefficients |
| 2    | #111 | (precondition) | time-domain Mueller-Murphy moment rate as Fourier pair of RDP, medium_type plumbed end-to-end |
| 3    | #112 | 4 (partial), 6 (partial) | per-layer Q to aux fields, t*(f) post-FFT envelope, MESH_REFINEMENT plumbing |
| 4    | #113-#115 | (closes pass-3 inversion) | multi-cell moment-tensor distribution, factor-30 envelope on anchor tests |
| 5    | #120 | 1 (partial) | DYNAMIC_PLASTIC config grammar, full 6-component Mdot_ij injection, R_cavity / R_plastic / Mdot history CSV; closed-form cavity-expansion kernel (axis-1a/1b deferred) |
| 6    | #121 | 1 (axis-1a partial) | RadialLagrangianSolver behind solver_kind dispatch (CLOSED_FORM default preserves pass-5 byte-for-byte; RADIAL_LAGRANGIAN opt-in runs the new shock solver), HDF5+XDMF spatial profile pair, six physics-validation gates, ClosedFormFallback + RadialLagrangianAnchor integration tests; far-field amplitude under RADIAL_LAGRANGIAN ~factor 100-400 below the closed-form estimate (axis-1a calibration follow-up; see fidelity doc pass-6 entry) |
| 7    | #122 (merged) | 1 (axis-1a closed) | TillotsonEOS host-rock evaluator + first-principles physics-based cavity initialization (Zel'dovich-Raizer end-state approximation, Newton energy-partition solve) + Wilkins (1980) literature AV defaults; Sedan 1962 amplitude ratio 2.24x (within factor-5 envelope); RADIAL_LAGRANGIAN promoted to default for DYNAMIC_PLASTIC; pass-6 standalone gates tightened to factor-10/30/3 envelopes; ParaView .pvsm verification (skipped in fsrm-ci, runnable via scripts/verify_pvsm.sh); strict pass-7 spec tolerances on six gates deferred to pass-8 (named candidates: tabulated plasma EOS, explicit Marshak phase, higher-order numerics, fitted alluvium Tillotson set) |
| 8    | (this PR) | 1 (axis-1a Marshak), 5 (V&V infrastructure) | Marshak grey radiation-diffusion solver (1D radial, implicit backward-Euler, tridiagonal Thomas, outer Newton on T^4) coupled to Tillotson host-rock matter via emission-absorption; explicit fidelity ladder (LOW Z-R / MED Marshak grey / HIGH multigroup scaffold / HIGHEST S_n named); five Marshak physics gates (self-similar pure radiation, energy conservation, hand-off debouncing, opacity regime coverage, grey-vs-Z-R cross-check); IRIS waveform V&V infrastructure (production SACReader + 5-metric WaveformComparison library + ObsPy refresh tool + cached-tarball layout); iris_validation CTest label with three Salmon 1964 gates (cavity radius vs measured 17.4 m within factor 5, free-field velocity GTEST_SKIP'd as axis-1b work, far-field mb +/- 0.4) plus Chagan / Pokhran I cross-validation; Tillotson plasma-extrapolation warning surfaced as one-time PETSc message |

When this pass merges, update this row with the merged PR number
and the per-axis row to reflect any scope shifts.
