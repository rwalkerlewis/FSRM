# Example 33: DPRK 2006 Nuclear Test

## Physics
Elastodynamic wave propagation through a 3-layer Mt. Mantap granite
column. Models the ~0.7 kt (USGS mb 4.1) DPRK event of 2006-10-09 at
~470 m depth in Punggye-ri, the first DPRK underground test.

Same Mt. Mantap velocity model as the rest of the DPRK series and
matches `config/examples/punggye_ri_layered.config`.

Features used:
- Layered heterogeneous material (aux field assignment)
- Underground nuclear explosion source (Mueller-Murphy, GRANITE medium)
- Absorbing boundary conditions (Clayton-Engquist)
- Seismometer network with SAC output

## Geology
| Layer | Depth (m) | Vp (m/s) | Vs (m/s) | Density (kg/m3) |
|-------|-----------|----------|----------|-----------------|
| Rhyolite / tuff cap   | 0-200    | 3500 | 2020 | 2200 |
| Fractured granite     | 200-1000 | 4500 | 2598 | 2500 |
| Competent granite     | 1000-2000 | 5800 | 3349 | 2650 |

References: Mt. Mantap velocity model from existing
`config/examples/punggye_ri_layered.config` (DPRK 2017 reference). USGS
event listing for the 2006-10-09 detection.

## Config
Uses `config/examples/dprk_2006.config`.

## Expected Output
- `output/dprk_2006/*.SAC` -- synthetic seismograms at 3 stations

## Running
```bash
./run.sh
```

## Verified By
- `Integration.HistoricNuclear.DPRK2006` -- Mt. Mantap granite + SAC output
