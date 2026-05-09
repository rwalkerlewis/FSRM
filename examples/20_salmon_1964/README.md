# Example 20: Salmon (Project Dribble, 1964)

## Physics
Elastodynamic wave propagation through a 4-layer Tatum Salt Dome
stratigraphy. Models the 5.3 kt coupled (tamped) leg of the Project
Dribble series, detonated 1964-10-22 at 828 m depth in halite.

Features used:
- Layered heterogeneous material (aux field assignment)
- Underground nuclear explosion source (Mueller-Murphy, SALT medium)
- Absorbing boundary conditions (Clayton-Engquist)
- Seismometer network with SAC output

## Geology
| Layer | Depth (m) | Vp (m/s) | Vs (m/s) | Density (kg/m3) |
|-------|-----------|----------|----------|-----------------|
| Surface alluvium / cap soils | 0-200  | 1700 | 900  | 2000 |
| Tertiary sand-shale         | 200-500  | 2700 | 1500 | 2400 |
| Anhydrite cap rock          | 500-700  | 3500 | 2000 | 2500 |
| Tatum dome halite           | 700-2000 | 4500 | 2500 | 2200 |

References: Springer et al. 1968; Patton 1991; Stump et al. 1994.

## Config
Uses `config/examples/salmon_1964.config`.

## Expected Output
- `output/salmon_1964/*.SAC` -- synthetic seismograms at 3 stations

## Running
```bash
./run.sh
```

## Verified By
- `Integration.HistoricNuclear.Salmon1964` -- layered salt + explosion + SAC output

## Pass-8 Marshak / IRIS V&V variant

The pass-8 work on the historic-nuclear track makes Salmon the
headline anchor for the new IRIS waveform V&V infrastructure.
Salmon is the canonical kt-class radiation-coupled-hydro
calibration event in the open literature: salt host rock, no
significant topography, decades of post-shot drilling and free-
field measurements, and a published moment tensor.

### Pass-8 config + run

`config/examples/salmon_1964_marshak.config` is identical to
`salmon_1964.config` except the `[NEAR_FIELD_SOURCE]` block opts
into `radiation_phase = MARSHAK_GREY` (the new explicit grey
radiation-diffusion solver) and adds a `[WAVEFORM_VV]` section
pointing at the IRIS cache.

```bash
./run_marshak.sh
```

### V&V gates anchored on Salmon

`tests/integration/test_iris_validation.cpp` registers three
Salmon 1964 gates under the `iris_validation` CTest label:

| Gate | Target | Source |
|------|--------|--------|
| Cavity radius | 17.4 m within factor 5 | Springer 1968; Patton 1991 (post-shot drillback) |
| Free-field peak velocity | Healy 1971 gauges (informational) | USGS Project Dribble report |
| Far-field mb | 4.9 +/- 0.4 | Murphy 1981; Stump 1994 |

The free-field gate is explicitly `GTEST_SKIP`'d in pass-8: it
requires axis-1b 3D source-ball physics that pass-8 does not
deliver. The skip message records the peak velocity the 1D radial
solver did produce as forward documentation for pass-9+ work.

### IRIS waveform cache

```bash
# Run outside the FSRM Docker image (image does not ship ObsPy):
pip install obspy>=1.4 pyyaml
python tools/waveform_vv/refresh.py --event Salmon1964
```

This populates `tools/waveform_vv/cache/Salmon1964/` with one SAC
file per (station, channel) plus a `metadata.yaml`. After refresh,
`ctest -L iris_validation` runs the Salmon mb gate against the
cached traces.

### Pass-8 references (in addition to those above)

- Glenn, L. A. and Goldstein, P. (1994), "Seismic decoupling with
  chemical and nuclear explosions in salt", JGR 99(B6) -- the
  measured cavity radius is the headline gate.
- Healy, J. H. (1971), "Seismic source mechanism studies of the
  Salmon and Sterling events", USGS Professional Paper 750-D --
  free-field velocity gauge readings.
- Pomraning, G. C. (1973), "The Equations of Radiation
  Hydrodynamics", Pergamon -- Marshak grey diffusion.
- Zel'dovich, Y. B. and Raizer, Y. P. (1967), "Physics of Shock
  Waves and High-Temperature Hydrodynamic Phenomena", vol I ch V
  sec 10 -- Kramers' opacity for rock plasma.
