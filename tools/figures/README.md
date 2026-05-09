# tools/figures

Shared infrastructure for the showcase-event figure scripts.

## Contents

- `figure_style.py`: shared matplotlib style (Wong 2011 colorblind-safe
  palette, fixed sizing, plot helpers `plot_seismogram`,
  `plot_radial_profile`, `savefig`). Imported by every per-event
  figure script.
- `requirements.txt`: Python dependencies for figure generation.

## Constraint: visualization, not analysis

The figure scripts read already-validated simulation output (HDF5,
CSV, SAC) and produce PNGs. They do not modify simulation output, do
not reimplement physics, and do not curve-fit results to look better.
They are visualization, not analysis.

If you find a script computing physics quantities, normalizing data
in ways that would change a V&V comparison, or fitting parameters to
match an observed reference, that is a bug. Move that work into the
simulator and persist the result through HDF5 / CSV.

## Generating figures

Each showcase event has its own `figures/` subdirectory under its
example folder, containing the per-event scripts. Run from the example
directory:

```bash
cd examples/20_salmon_1964
./run_showcase.sh         # runs the simulation + regenerates figures
# or just the figures, if output/ already exists:
./figures/regenerate.sh
```

## Reproducibility

The scripts write PNGs at 300 dpi with fixed `rcParams` and the SVG
font-type set to "none" so committed PNGs are deterministic on a
given matplotlib + DejaVu Sans + numpy + h5py + obspy stack. Slight
metadata differences between Linux and macOS PDFs are immaterial.

## Adding a new showcase event

1. Create `examples/N_<event>/figures/scripts/` and copy the six
   plotting scripts from a similar showcase as the starting point.
2. Update the readers / station names / titles to match the new event.
3. Add `examples/N_<event>/figures/regenerate.sh` that calls each
   script with the appropriate args.
4. Add `examples/N_<event>/run_showcase.sh` that builds output/, runs
   the simulation, and calls `figures/regenerate.sh`.
5. Add a "Figure catalog" section to the example's README.md.
