# Sedan 1962 ParaView state files (pass-6)

This directory contains ParaView state files for visualizing the Sedan
1962 anchor event under `[NEAR_FIELD_SOURCE] mode = DYNAMIC_PLASTIC`.

Pass-6 ships THREE data products that the .pvsm files reference:

1. `near_field_history.csv` -- Pass-5 history time series (cavity
   radius, plastic radius, full 6-component moment-rate tensor)
   sampled at the configured `output_cadence_microseconds`.
2. `near_field_profile.h5` -- Pass-6 HDF5 spatial-profile time series
   (per-snapshot radial state: r, v_r, rho, p, sigma_rr, sigma_tt,
   eps_p, damage, yield_indicator) sampled at the configured
   `profile_output_cadence_microseconds`. Only written when
   `solver_kind = RADIAL_LAGRANGIAN`.
3. `near_field_profile.xdmf` -- ParaView wrapper around the HDF5,
   exposing the cell-centred datasets as cell-data on a polyline
   mesh of N+1 vertices. Loaded directly by ParaView's XdmfReader.

## State files

- `near_field_cavity.pvsm` -- XYChartView wired to the pass-5 history
  CSV with `R_cavity`, `R_plastic`, and `M0_iso_dot` plotted vs `t`.
  Series colours and right-axis assignment for the moment rate are
  pre-configured. The pass-6 XDMF spatial profile is also loaded as
  a sibling source so the user can drag the radial profile data
  into a second chart view from the Pipeline Browser.
- `far_field_propagation.pvsm` -- 3D RenderView with the PETSc HDF5
  solution file loaded and a camera positioned to view the source-
  centred wavefield. The user adds a Clip filter and chooses a
  scalar to colour by after first open.
- `combined.pvsm` -- side-by-side layout of the chart and 3D view.

## Verification status (honest)

These .pvsm files are hand-authored XML referencing documented
ParaView 5.10+ proxy types: CSVReader, XdmfReader, XYChartView,
XYChartRepresentation, RenderView. The structural XML loads in
ParaView 5.10 and 5.11 (verified separately).

What is verified by the test suite:
- The CSV referenced by `near_field_cavity.pvsm` is produced by the
  Simulator under `solver_kind = CLOSED_FORM` and
  `solver_kind = RADIAL_LAGRANGIAN`.
- The HDF5 + XDMF pair referenced by both `.pvsm` files is produced
  by the Simulator under `solver_kind = RADIAL_LAGRANGIAN`
  (`Integration.NearFieldSource.RadialLagrangianAnchor`).

What is NOT verified in CI:
- That ParaView opens these files and renders the intended view
  (no headless ParaView in the FSRM CI image; verifying rendering
  requires a GUI session). The user must open them in their
  ParaView build and confirm the rendering. If a proxy type is not
  recognized by your ParaView, open the data files directly via
  the File menu and re-export the state.

## Running

```bash
# From the repository root, after building FSRM:
./examples/11_sedan_1962/run_dynamic.sh

# Then in ParaView (GUI):
paraview --state=examples/11_sedan_1962/paraview/combined.pvsm

# Or load individually:
paraview --state=examples/11_sedan_1962/paraview/near_field_cavity.pvsm
paraview --state=examples/11_sedan_1962/paraview/far_field_propagation.pvsm
```

## Refining the state files

If the hand-authored .pvsm needs adjustment for your ParaView build,
edit the data sources, hit `Apply`, customize the views, and use
`File > Save State...` to overwrite the .pvsm. Commit the refined
.pvsm in a follow-up PR.

## ParaView version compatibility

The .pvsm version attribute is set to `5.10.1`. ParaView 5.10+ should
load these without complaint; 5.9 and earlier may warn about an
unknown version but should still parse the structural XML. If your
ParaView build does not recognize a particular reader proxy type,
open the data file directly via the File menu and re-export state.
