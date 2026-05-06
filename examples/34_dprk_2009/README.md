# Example 34: DPRK 2009 Nuclear Test

## Physics
Elastodynamic wave propagation through the Mt. Mantap granite column.
Models the ~4 kt (USGS mb 4.5) DPRK event of 2009-05-25 at ~550 m depth
in Punggye-ri.

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
2009-05-25 detection.

## Config
Uses `config/examples/dprk_2009.config`.

## Expected Output
- `output/dprk_2009/*.SAC` -- synthetic seismograms at 3 stations

## Running
```bash
./run.sh
```

## Verified By
- `Integration.HistoricNuclear.DPRK2009` -- Mt. Mantap granite + SAC output
