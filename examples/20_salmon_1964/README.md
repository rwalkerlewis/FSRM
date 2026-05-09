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
