# Example 30: Chagan ("Atomic Lake", 1965)

## Physics
Elastodynamic wave propagation through a 4-layer Semipalatinsk Test
Site sandstone-siltstone column with a permafrost surface zone. Models
the 140 kt Chagan event of 1965-01-15 at 178 m depth -- the largest
Soviet peaceful-use cratering experiment, analogous in scale and
purpose to Sedan 1962.

Features used:
- Layered heterogeneous material (aux field assignment)
- Underground nuclear explosion source (Mueller-Murphy, ALLUVIUM medium)
- Absorbing boundary conditions (Clayton-Engquist)
- Seismometer network with SAC output

## Geology
| Layer | Depth (m) | Vp (m/s) | Vs (m/s) | Density (kg/m3) |
|-------|-----------|----------|----------|-----------------|
| Permafrost / alluvium    | 0-50    | 1700 | 850  | 1900 |
| Sandstone                | 50-300  | 2800 | 1500 | 2300 |
| Siltstone                | 300-1000 | 3500 | 2000 | 2400 |
| Pre-Mesozoic basement    | 1000-2000 | 4500 | 2600 | 2600 |

Reference: Adushkin & Spivak 2015 (Semipalatinsk well-log summary;
Soviet PNE cratering shots).

## Config
Uses `config/examples/chagan_1965.config`.

## Expected Output
- `output/chagan_1965/*.SAC` -- synthetic seismograms at 3 stations

## Running
```bash
./run.sh
```

## Verified By
- `Integration.HistoricNuclear.Chagan1965` -- Soviet sandstone cratering + SAC output
