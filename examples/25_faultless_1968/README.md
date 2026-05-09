# Example 25: Faultless (1968)

## Physics
Elastodynamic wave propagation through a 3-layer
alluvium/Tertiary-volcanics/basement column at Hot Creek Valley NV.
Models the ~1 Mt Faultless event of 1968-01-19 at 975 m depth in the
Central Nevada Test Area. Triggered visible surface fault offset along
pre-existing scarps; ended the CNTA campaign as a single shot.

Features used:
- Layered heterogeneous material (aux field assignment)
- Underground nuclear explosion source (Mueller-Murphy, ALLUVIUM medium)
- Absorbing boundary conditions (Clayton-Engquist)
- Seismometer network with SAC output

## Geology
| Layer | Depth (m) | Vp (m/s) | Vs (m/s) | Density (kg/m3) |
|-------|-----------|----------|----------|-----------------|
| Quaternary alluvium         | 0-300    | 1800 | 900  | 1900 |
| Tertiary volcanics          | 300-1500 | 4000 | 2200 | 2400 |
| Pre-Tertiary basement       | 1500-3000 | 5500 | 3100 | 2700 |

References: Hamilton & Healy 1969; McKeown & Dickey 1969.

## Config
Uses `config/examples/faultless_1968.config`.

## Expected Output
- `output/faultless_1968/*.SAC` -- synthetic seismograms at 3 stations

## Running
```bash
./run.sh
```

## Verified By
- `Integration.HistoricNuclear.Faultless1968` -- megaton CNTA + explosion + SAC output
