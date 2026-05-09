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

## 1. Dynamic near-field source -- pass-5 (partial)

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

**Deferred to pass-N+1.** The pass-5 `DYNAMIC_PLASTIC` path uses
the existing 1D solver whose internal `step()` integrates a
closed-form exponential cavity-expansion kernel rather than a
finite-difference radial momentum + constitutive update. The
strength model, damage model, and EOS data structures plumbed
through the solver are referenced in the recorded diagnostics but
do not drive the cavity expansion itself. Two follow-ups:

  - **Pass-N+1a.** Replace the closed-form kernel with a true 1D
    radial Lagrangian finite-volume solver: explicit shock-friendly
    time-stepping, Drucker-Prager radial-return per cell,
    Mie-Gruneisen EOS evaluation per cell, surface-integral moment
    tensor extraction. The pass-5 dispatch path and history CSV
    do not need to change; only the substitution at the solver
    interface.
  - **Pass-N+1b.** Replace the 1D radial solver with a 3D
    Drucker-Prager subdomain on the source ball, dropping the
    spherical-symmetry assumption. This is the spec-as-written
    ambition of the pass-5 task and is multi-week scope; the pass-5
    config grammar and dispatch path were authored to substitute
    cleanly into either an axis-1a or axis-1b implementation.

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
| 5    | (this PR) | 1 (partial) | DYNAMIC_PLASTIC config grammar, full 6-component Mdot_ij injection, R_cavity / R_plastic / Mdot history CSV; closed-form cavity-expansion kernel (axis-1a/1b deferred) |

When pass-5 merges, update this row with the merged PR number and
the per-axis row to reflect any scope shifts.
