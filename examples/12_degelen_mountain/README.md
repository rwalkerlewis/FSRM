# Example 12: Degelen Mountain, Semipalatinsk (Kazakhstan)

## Physics
Elastodynamic wave propagation through a 3-layer granitic geology model
representing the Degelen Mountain testing complex at the Semipalatinsk
(now Semey) Test Site in Kazakhstan. Models a typical 50 kt tunnel test
at 300m depth.

Features used:
- Layered heterogeneous material (aux field assignment)
- Underground nuclear explosion source (Mueller-Murphy)
- Absorbing boundary conditions (Clayton-Engquist)
- Seismometer network with SAC output

## Geology
| Layer | Depth (m) | Vp (m/s) | Vs (m/s) | Density (kg/m3) |
|-------|-----------|----------|----------|-----------------|
| Weathered Granite | 0-100 | 3000 | 1732 | 2200 |
| Fractured Granite | 100-500 | 4800 | 2771 | 2550 |
| Competent Granite/Gneiss | 500-5000 | 6000 | 3464 | 2700 |

## Config
Uses `config/examples/degelen_mountain.config`.

## Expected Output
- `output/degelen_mountain/*.SAC` -- synthetic seismograms at 3 stations

## Running
```bash
./run.sh
```

## Verified By
- `Integration.HistoricNuclear.DegelenMountain` -- layered + explosion + SAC output
