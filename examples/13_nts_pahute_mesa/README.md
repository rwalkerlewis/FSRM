# Example 13: NTS Pahute Mesa, Nevada

## Physics
Elastodynamic wave propagation through a 4-layer volcanic geology model
representing the Pahute Mesa area of the Nevada Test Site. Models a
typical 150 kt shaft test at 600m depth in zeolitized tuff. Pahute Mesa
hosted many of the largest US underground tests.

Features used:
- Layered heterogeneous material (aux field assignment)
- Underground nuclear explosion source (Mueller-Murphy)
- Absorbing boundary conditions (Clayton-Engquist)
- Seismometer network with SAC output

## Geology
| Layer | Depth (m) | Vp (m/s) | Vs (m/s) | Density (kg/m3) |
|-------|-----------|----------|----------|-----------------|
| Alluvium | 0-200 | 2200 | 1100 | 1800 |
| Welded Tuff | 200-400 | 4000 | 2300 | 2300 |
| Zeolitized Tuff | 400-1000 | 3200 | 1800 | 2100 |
| Rhyolite | 1000-5000 | 5200 | 3000 | 2600 |

## Config
Uses `config/examples/nts_pahute_mesa.config`.

## Expected Output
- `output/nts_pahute_mesa/*.SAC` -- synthetic seismograms at 3 stations

## Running
```bash
./run.sh
```

## Verified By
- `Integration.HistoricNuclear.NtsPahuteMesa` -- layered + explosion + SAC output
