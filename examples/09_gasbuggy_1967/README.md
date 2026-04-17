# Example 09: Project Gasbuggy (1967)

## Physics
Elastodynamic wave propagation through a 4-layer geology model representing
the San Juan Basin, New Mexico. Models the 29 kt Gasbuggy underground
nuclear test from Project Plowshare, detonated at 1280m depth in the Lewis
Shale formation.

Features used:
- Layered heterogeneous material (aux field assignment)
- Underground nuclear explosion source (Mueller-Murphy)
- Absorbing boundary conditions (Clayton-Engquist)
- Seismometer network with SAC output

## Geology
| Layer | Depth (m) | Vp (m/s) | Vs (m/s) | Density (kg/m3) |
|-------|-----------|----------|----------|-----------------|
| Alluvium/Sandstone | 0-200 | 2500 | 1200 | 2100 |
| Mesaverde Sandstone | 200-800 | 4000 | 2300 | 2450 |
| Lewis Shale | 800-1600 | 3500 | 2000 | 2550 |
| Pictured Cliffs Ss | 1600-5000 | 4500 | 2600 | 2500 |

## Config
Uses `config/examples/gasbuggy_1967.config`.

## Expected Output
- `output/gasbuggy_1967/*.SAC` -- synthetic seismograms at 3 stations

## Running
```bash
./run.sh
```

## Verified By
- `Integration.HistoricNuclear.Gasbuggy1967` -- layered + explosion + SAC output
