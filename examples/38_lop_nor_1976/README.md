# Example 38: Lop Nor 1976-10-17 (PRC First Underground Test)

## Physics
Elastodynamic wave propagation through a 4-layer Beishan-equivalent
granite column. Models the ~4 kt PRC event of 1976-10-17 at ~250 m
depth -- the first Chinese underground nuclear test.

Features used:
- Layered heterogeneous material (aux field assignment)
- Underground nuclear explosion source (Mueller-Murphy, GRANITE medium)
- Absorbing boundary conditions (Clayton-Engquist)
- Seismometer network with SAC output

## Geology
| Layer | Depth (m) | Vp (m/s) | Vs (m/s) | Density (kg/m3) |
|-------|-----------|----------|----------|-----------------|
| Surface alluvium / talus | 0-50    | 1500 | 750  | 1900 |
| Weathered granite        | 50-300  | 3500 | 1900 | 2400 |
| Competent granite        | 300-1500 | 5500 | 3100 | 2700 |
| Deep granite             | 1500-2000 | 6000 | 3400 | 2750 |

References: Wen et al. 2006; Yang et al. 2003 (Lop Nor basement
velocity structure).

## Config
Uses `config/examples/lop_nor_1976.config`.

## Expected Output
- `output/lop_nor_1976/*.SAC` -- synthetic seismograms at 3 stations

## Running
```bash
./run.sh
```

## Verified By
- `Integration.HistoricNuclear.LopNor1976` -- Lop Nor granite + SAC output
