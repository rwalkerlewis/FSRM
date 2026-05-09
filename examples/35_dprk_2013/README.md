# Example 35: DPRK 2013 Nuclear Test

## Physics
Elastodynamic wave propagation through the Mt. Mantap granite column.
Models the ~10 kt (USGS mb 5.1) DPRK event of 2013-02-12 at ~600 m
depth in Punggye-ri.

Same Mt. Mantap velocity model as the rest of the DPRK series.

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
`config/examples/punggye_ri_layered.config`. USGS event listing for the
2013-02-12 detection.

## Config
Uses `config/examples/dprk_2013.config`.

## Expected Output
- `output/dprk_2013/*.SAC` -- synthetic seismograms at 3 stations

## Running
```bash
./run.sh
```

## Verified By
- `Integration.HistoricNuclear.DPRK2013` -- Mt. Mantap granite + SAC output
