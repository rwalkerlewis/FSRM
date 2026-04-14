# Example 11: Sedan Crater (1962)

## Physics
Elastodynamic wave propagation through a 3-layer geology model representing
NTS Yucca Flat. Models the 104 kt Sedan underground nuclear test from
Operation Storax/Plowshare, detonated at only 194m depth in alluvium.
This shallow burial produced the largest US nuclear crater (390m diameter,
100m deep).

Features used:
- Layered heterogeneous material (aux field assignment)
- Underground nuclear explosion source (Mueller-Murphy)
- Absorbing boundary conditions (Clayton-Engquist)
- Seismometer network with SAC output

## Geology
| Layer | Depth (m) | Vp (m/s) | Vs (m/s) | Density (kg/m3) |
|-------|-----------|----------|----------|-----------------|
| Alluvium | 0-300 | 1800 | 800 | 1700 |
| Welded Tuff | 300-800 | 4200 | 2400 | 2300 |
| Paleozoic Carbonate | 800-5000 | 5500 | 3100 | 2700 |

## Config
Uses `config/examples/sedan_1962.config`.

## Expected Output
- `output/sedan_1962/*.SAC` -- synthetic seismograms at 3 stations

## Running
```bash
./run.sh
```

## Verified By
- `Integration.HistoricNuclear.Sedan1962` -- layered + explosion + SAC output
