# Sedan 1962 ParaView state files (pass-5)

This directory contains ParaView state files for visualizing the Sedan
1962 anchor event under `[NEAR_FIELD_SOURCE] mode = DYNAMIC_PLASTIC`.

Three .pvsm files are provided, each minimal-but-valid ParaView 5.10+
state-file XML:

- `near_field_cavity.pvsm` -- XY chart view loading
  `../output/sedan_1962_dynamic/near_field_history.csv` (the pass-5
  near-field history CSV: cavity radius, plastic radius, and
  6-component moment-rate tensor as a function of time).
- `far_field_propagation.pvsm` -- 3D render view loading
  `../output/solution.h5` (the PETSc HDF5 mesh + solution file with
  displacement / velocity).
- `combined.pvsm` -- multi-view layout combining the two.

## Honest scope

These .pvsm files are minimal stubs. They configure the data sources,
camera, and view-type defaults, but they do not pre-configure
representations, colormaps, isosurfaces, or station-point glyphs.
They are committed as documentation of how the saved data is meant
to be viewed; the user should open them in ParaView, add the
visual elements that match the comments inside each file, and
re-export the refined state if a more complete .pvsm is desired.

The .pvsm format requires GUI-level proxy / representation /
ColorTransferFunction state to render anything; hand-authoring a
fully-rendered state without ParaView available is fragile. Pass-5
ships the data + the view-type-and-camera scaffolding; pass-N+1 may
ship full curated state files once the dynamic-plastic pipeline has
production runs to render against.

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

When you open a .pvsm in ParaView and want to keep your refinements,
use `File > Save State...` to overwrite the .pvsm with the GUI-edited
version. Commit the refined .pvsm in a follow-up PR.

## What to add after first open

For `near_field_cavity.pvsm`:
- In the Properties panel for the CSVReader, click "Apply".
- Open `Plot Data` (Filters > Plotting > Plot Data) on the CSV.
- Set Series: enable `R_cavity`, `R_plastic`, optionally `M0_iso_dot`.
- Use the right-axis toggle to plot `M0_iso_dot` on a secondary
  axis if the moment-rate magnitudes dwarf the radii in metres.

For `far_field_propagation.pvsm`:
- The HDF5 reader may need to be switched to "Generic HDF5 Reader"
  or "VisIt HDF5 Reader" depending on your ParaView build.
- Apply a `Clip` filter (Filters > Common > Clip) through the source
  centre (15000, 15000, 4806).
- Colour by `displacement` magnitude or `velocity` magnitude.
- Add seismometer station glyphs (sphere sources at):
  - SPALL: (15000, 15000, 5000)
  - STA05: (20000, 15000, 5000)
  - STA10: (25000, 15000, 5000)

For `combined.pvsm`:
- Apply both refinements above; the two views are pre-laid-out
  side-by-side.

## Known limitations

The .pvsm version attribute is set to `5.10.1`. ParaView 5.10+ should
load these without complaint; 5.9 and earlier may warn about an
unknown version but should still parse the structural XML. If your
ParaView build does not recognize a particular reader proxy type
(e.g. `VisItPDBReader`), open the data file directly via the File
menu and re-export state.
