# Example 20: Salmon 1964 (Project Dribble) -- showcase event

The canonical Marshak + IRIS V&V anchor of the historic-nuclear
program. Salmon 1964 is the kt-class radiation-coupled-hydro
calibration event in the open literature: salt host rock, no
significant topography, decades of post-shot drilling and free-field
measurements, and a published moment tensor. The companion event
Sterling 1966 (`examples/21_sterling_1966/`) was detonated in the
Salmon-generated cavity and recovers the decoupling factor.

## Physics summary

| Property | Value |
|---|---|
| Date | 1964-10-22 |
| Site | Tatum Salt Dome, Lamar County, Mississippi |
| Yield | 5.3 kt (coupled / tamped leg of Project Dribble) |
| Depth | 828 m |
| Host rock | Halite (Tatum Salt Dome) |
| Body-wave magnitude | mb 4.9 (Murphy 1981; Stump 1994) |

## What the simulator does for Salmon

Salmon is the showcase event for the source-physics fidelity ladder:

- **LOW tier (`config.config`)**: pass-7 Z-R closed-form end-state +
  Tillotson cavity EOS. Default invocation produces SAC seismograms
  at the configured stations and matches the pinned pass-7 envelope
  byte-for-byte.
- **MED tier (`config_marshak.config`)**: opts into `radiation_phase =
  MARSHAK_GREY` (the explicit grey radiation-diffusion solver) plus
  the `[WAVEFORM_VV]` block pointing at the IRIS cache. This is the
  pass-8 default for Marshak gates.
- **HIGHEST tier (`config_highest.config`)**: opts into multigroup
  radiation transport, BDF2 implicit diffusion, RK3-SSP hydro, and
  STRANG_MULTIGROUP splitting. Reproduces the pass-11 axis-1c result
  on the Marshak self-similar gate (factor 2) and the radiation-
  energy conservation gate (2 %).

## Geology / velocity model

| Layer | Depth (m) | Vp (m/s) | Vs (m/s) | rho (kg/m^3) |
|---|---|---:|---:|---:|
| Surface alluvium / cap soils | 0-200 | 1700 | 900 | 2000 |
| Tertiary sand-shale | 200-500 | 2700 | 1500 | 2400 |
| Anhydrite cap rock | 500-700 | 3500 | 2000 | 2500 |
| Tatum dome halite | 700-2000 | 4500 | 2500 | 2200 |

References: Springer et al. 1968; Patton 1991; Stump et al. 1994.

## V&V gates

| Gate | Pass-11 envelope | Spec target | Status |
|---|---|---|---|
| Salmon CavityRadius | factor 3 | 5 % | open (axis-1b 3D source ball) |
| FreeFieldPeakVelocity_166m | factor 2 | factor 2 | closed (pass-9) |
| FreeFieldPeakVelocity_322m | factor 2 | factor 2 | closed (pass-9) |
| FreeFieldPeakVelocity_549m | factor 4 | factor 2 | open (impedance BC) |
| FarFieldBodyWaveMagnitude | +/-0.3 | +/-0.2 | open (axis-3 path) |
| Marshak SelfSimilar | factor 2 (BDF2) | factor 2 | closed pass-11 |
| Marshak EnergyConservation | 2 % (BDF2) | 2 % | closed pass-11 |
| Marshak GreyVsZR | factor 3 | factor 3 | closed |

Source of truth: [docs/AXIS_1A_FIDELITY_REPORT.md](../../docs/AXIS_1A_FIDELITY_REPORT.md).

## Running the showcase

End-to-end, from a clean tree:

```bash
cd examples/20_salmon_1964
MPI_RANKS=8 ./run_showcase.sh   # default 8 ranks for showcase runs
```

This produces three output trees (`output/`, `output_marshak/`,
`output_highest/`) and runs `figures/regenerate.sh` to produce six
PNGs in `figures/`.

For a single-tier MED run:

```bash
./run_marshak.sh
```

For a single-tier LOW run (default config):

```bash
./run.sh
```

## Output catalog

Each tier produces:

- `seismograms/*.sac`: SAC traces at every configured station.
- `near_field_history.csv`: 6-component M(t), Mdot(t), R_cav(t).
- `near_field_profile.h5` + `.xdmf`: per-snapshot radial state of the
  source ball (open in ParaView via XDMF wrapper).
- `solution.h5` + `.xmf`: full 3-D wavefield.

## Figure catalog

Six PNGs land in `figures/` after `figures/regenerate.sh`:

| File | What it shows | Source data |
|---|---|---|
| `01_geology_cross_section.png` | Layered velocity model + source location | `config.config` |
| `02_cavity_formation_radial_profile.png` | 4-panel radial profile at four times | `output/near_field_profile.h5` |
| `03_moment_tensor_history.png` | M(t) and Mdot(t) for six independent components | `output/near_field_history.csv` |
| `04_synthetic_seismograms_grid.png` | All-station seismogram small multiples | `output/seismograms/*.sac` |
| `05_synthetic_vs_observed.png` | Synthetic vs observed at three named stations with cross-correlation | `output/seismograms/*.sac` + `tools/waveform_vv/cache/salmon_1964/*.sac` |
| `06_fidelity_ladder_comparison.png` | LOW vs MED vs HIGHEST peak velocity overlay | `output*/seismograms/*.sac` |

The figure scripts under `figures/scripts/` import the shared style
infrastructure in `tools/figures/figure_style.py`. They are
visualization, not analysis.

## IRIS waveform cache

```bash
# Outside the FSRM Docker image (image does not ship ObsPy):
pip install -r tools/figures/requirements.txt
python tools/waveform_vv/refresh.py --event Salmon1964
```

Populates `tools/waveform_vv/cache/salmon_1964/` with one SAC per
(station, channel) plus inventory metadata. After refresh,
`ctest -L iris_validation` runs the Salmon waveform-comparison gates.

## Verified by

- `Integration.HistoricNuclear.Salmon1964` -- layered salt + explosion + SAC output
- `Physics.WaveformVV.Salmon1964.*` (under `iris_validation` ctest label) when cache is populated
- `Integration.MPI.SalmonSerialVsParallelEquivalence` -- pass-12 parallel correctness gate

## References

- Springer, D. L., et al. (1968), "Seismic source summary for U.S.
  underground nuclear explosions", BSSA 58.
- Patton, H. J. (1991), "Reassessment of seismic moments and source
  spectra for the Salmon and Sterling events", BSSA 81.
- Stump, B. W., et al. (1994), "The Salmon experiment: review of
  free-field measurements", U.S. Geological Survey.
- Healy, J. H. (1971), "Seismic source mechanism studies of the
  Salmon and Sterling events", USGS Professional Paper 750-D.
- Glenn, L. A. and Goldstein, P. (1994), "Seismic decoupling with
  chemical and nuclear explosions in salt", JGR 99(B6).
- Murphy, J. R. (1981), "P-wave coupling of underground explosions
  in various geologic media", in Identification of Seismic Sources.
- Pomraning, G. C. (1973), "The Equations of Radiation
  Hydrodynamics", Pergamon (Marshak grey diffusion).
- Zel'dovich, Y. B. and Raizer, Y. P. (1967), "Physics of Shock
  Waves and High-Temperature Hydrodynamic Phenomena", vol I.
