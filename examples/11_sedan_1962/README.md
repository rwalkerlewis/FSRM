# Example 11: Sedan Crater (1962)

## Physics
Elastodynamic wave propagation through a 3-layer geology model representing
NTS Yucca Flat. Models the 104 kt Sedan underground nuclear test from
Operation Storax/Plowshare, detonated at only 194m depth in alluvium.
This shallow burial produced the largest US nuclear crater (390m diameter,
100m deep).

Features used:
- Layered heterogeneous material (aux field assignment)
- Underground nuclear explosion source (Mueller-Murphy)
- Absorbing boundary conditions (Clayton-Engquist)
- Seismometer network with SAC output

## Geology
| Layer | Depth (m) | Vp (m/s) | Vs (m/s) | Density (kg/m3) |
|-------|-----------|----------|----------|-----------------|
| Alluvium | 0-300 | 1800 | 800 | 1700 |
| Welded Tuff | 300-800 | 4200 | 2400 | 2300 |
| Paleozoic Carbonate | 800-5000 | 5500 | 3100 | 2700 |

## Configs

Two configs ship with this example:

- `examples/11_sedan_1962/config.config` -- baseline KINEMATIC_RDP path.
- `examples/11_sedan_1962/config.config` -- pass-5 DYNAMIC_PLASTIC
  path. Pass-7 promoted RADIAL_LAGRANGIAN to the new default for
  DYNAMIC_PLASTIC, but this fixture pins `solver_kind = CLOSED_FORM`
  explicitly so the legacy RDP-driven cavity radius and far-field
  amplitude assertions in `Integration.HistoricNuclear.Sedan1962_Dynamic`
  continue to gate pass-5 byte-identical behavior.

Sedan is the pass-5/6/7 anchor event because the cratering shot
exercises near-surface damage and the alluvium medium puts the
cavity radius at a tractable scale (Rc ~ 87 m via the medium-aware
NTS coefficient) relative to the far-field cell scale (Lz / nz =
500 m). See `docs/HISTORIC_NUCLEAR_FIDELITY.md` "Closed in pass 7"
and `docs/HISTORIC_NUCLEAR_ROADMAP.md` axis 1 for the rationale.

## Expected Output

Baseline (`run.sh`):
- `output/sedan_1962/*.SAC` -- synthetic seismograms at 3 stations

Dynamic-plastic (`run_dynamic.sh`, fixture pinned to
`solver_kind = CLOSED_FORM`):
- `output/sedan_1962_dynamic/*.SAC` -- synthetic seismograms
- `output/sedan_1962_dynamic/near_field_history.csv` -- recorded 1D
  solver moment-rate tensor and cavity / plastic radius time series

Pass-7 RADIAL_LAGRANGIAN opt-in (delete or override the explicit
`solver_kind = CLOSED_FORM` line in the dynamic config):
- All of the above, with the moment-tensor history sourced from the
  1D radial Lagrangian shock solver under the pass-7 defaults
  (Tillotson host-rock EOS, physics-based energy-partition cavity
  initialization, Wilkins 1980 literature AV coefficients).
- `output/sedan_1962_dynamic/near_field_profile.h5` -- pass-6 HDF5
  spatial-profile time series (radial state at the configured
  cadence: r, v_r, rho, p, sigma_rr, sigma_tt, eps_p, damage,
  yield_indicator). Schema in
  `include/domain/explosion/RadialLagrangianOutput.hpp`.
- `output/sedan_1962_dynamic/near_field_profile.xdmf` -- ParaView
  wrapper around the HDF5.

## Running

```bash
# Baseline KINEMATIC_RDP path:
./run.sh

# Pass-5/6/7 DYNAMIC_PLASTIC path (the fixture pins
# solver_kind = CLOSED_FORM):
./run_dynamic.sh

# Pass-7 RADIAL_LAGRANGIAN: open
# examples/11_sedan_1962/config.config and either remove the
# `solver_kind = CLOSED_FORM` line (the new default is
# RADIAL_LAGRANGIAN) or replace it with `solver_kind =
# RADIAL_LAGRANGIAN`. Add `radial_cells = 200` and
# `profile_output_cadence_microseconds = 5000` to capture the
# spatial-profile output.

# Pass-7 amplitude diagnostic outside CI (Sedan 1962 fixture, both
# solver_kinds, prints the peak |M0_iso_dot| ratio):
scripts/pass7_amplitude_diagnostic.sh 200
```

## Visualization

Three ParaView state files in `paraview/` configure rendering of the
dynamic-plastic outputs:

- `paraview/near_field_cavity.pvsm` -- XY chart of cavity / plastic
  radius and isotropic moment rate from `near_field_history.csv`,
  plus the pass-6 spatial-profile XDMF as a sibling source. Series
  visibility, colours, and right-axis assignment for the moment
  rate are pre-configured.
- `paraview/far_field_propagation.pvsm` -- 3D wavefield render of
  `output/solution.h5`.
- `paraview/combined.pvsm` -- multi-view layout of both.

These are hand-authored XML referencing documented ParaView 5.10+
proxy types (CSVReader, XdmfReader, XYChartView,
XYChartRepresentation, RenderView). The XML is structurally valid
but the rendered output is not verified in CI (the FSRM build
environment does not run ParaView headless); see `paraview/README.md`
for the verification status and what to confirm after first open.

## Verified By

- `Integration.HistoricNuclear.Sedan1962` -- baseline layered + explosion + SAC output
- `Integration.HistoricNuclear.Sedan1962_Distributed` -- pass-4 multi-cell
  moment-tensor distribution at the factor-30 envelope
- `Integration.HistoricNuclear.Sedan1962_Dynamic` -- pass-5
  DYNAMIC_PLASTIC source: cavity radius within 20% of medium-aware
  NTS analytic, peak/u_far within factor 30, near_field_history.csv
  emitted with > 100 sample rows. Continues to pass under the pass-6
  default `solver_kind = CLOSED_FORM`.
- `Integration.NearFieldSource.ClosedFormFallback` (pass-6 / pass-7)
  -- byte-identical regression guard. Pass-7 promoted
  RADIAL_LAGRANGIAN to the default for DYNAMIC_PLASTIC, so the test
  now compares two explicit-CLOSED_FORM runs to keep the pass-5
  byte-identical guarantee meaningful.
- `Integration.NearFieldSource.RadialLagrangianAnchor` (pass-7) --
  factor-5 envelope: peak |M0_iso_dot| ratio between
  RADIAL_LAGRANGIAN and CLOSED_FORM is asserted to lie inside
  [0.2, 5.0] on the same Sedan 1962 fixture. Lands at ~ 2.24x at
  pass-7 defaults.
- `Functional.ParaView.NearFieldCavityStateLoads` (pass-7) -- pvpython
  loads `paraview/near_field_cavity.pvsm` cleanly. Recorded as
  GTEST_SKIP under fsrm-ci (no ParaView in the image); developers
  run `scripts/verify_pvsm.sh` standalone for the verification.
