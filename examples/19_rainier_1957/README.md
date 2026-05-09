# Example 19: Rainier (1957)

## Physics
Elastodynamic wave propagation through a 3-layer tuff stratigraphy at
NTS Area 12 (Rainier Mesa). Models the 1.7 kt Rainier event of
1957-09-19, the first fully contained US underground nuclear test.

Features used:
- Layered heterogeneous material (aux field assignment)
- Underground nuclear explosion source (Mueller-Murphy, TUFF medium)
- Absorbing boundary conditions (Clayton-Engquist)
- Seismometer network with SAC output

## Geology
| Layer | Depth (m) | Vp (m/s) | Vs (m/s) | Density (kg/m3) |
|-------|-----------|----------|----------|-----------------|
| Surface alluvium      | 0-50    | 1500 | 750  | 1900 |
| Welded Rainier Mesa Tuff | 50-200 | 2700 | 1500 | 2200 |
| Bedded tunnel-bed tuff | 200-2000 | 3300 | 1800 | 2300 |

References: Springer & Kinnaman 1971; Carlos 1995.

## Config
Uses `config/examples/rainier_1957.config`.

## Expected Output
- `output/rainier_1957/*.SAC` -- synthetic seismograms at 3 stations

## Running
```bash
./run.sh
```

## Verified By
- `Integration.HistoricNuclear.Rainier1957` -- layered tuff + explosion + SAC output
