# Example 28: Rulison (1969)

## Physics
Elastodynamic wave propagation through a 5-layer Piceance Basin
sandstone-shale column. Models the 40 kt Rulison event of 1969-09-10
at 2568 m depth in the tight-gas Williams Fork (Mesaverde) Sandstone
-- the second of three Plowshare gas-stimulation tests.

Features used:
- Layered heterogeneous material (aux field assignment)
- Underground nuclear explosion source (Mueller-Murphy, SHALE medium)
- Absorbing boundary conditions (Clayton-Engquist)
- Seismometer network with SAC output

## Geology
| Layer | Depth (m) | Vp (m/s) | Vs (m/s) | Density (kg/m3) |
|-------|-----------|----------|----------|-----------------|
| Surface clays / colluvium  | 0-100    | 1700 | 850  | 2000 |
| Tertiary fluvial sediments | 100-500  | 2500 | 1400 | 2300 |
| Fort Union Formation       | 500-1500 | 3500 | 2000 | 2400 |
| Williams Fork (Mesaverde)  | 1500-3000 | 4500 | 2500 | 2500 |
| Mancos Shale               | 3000-5500 | 4200 | 2300 | 2500 |

References: Lombard et al. 1971; Reynolds et al. 1971 (Rulison
Piceance Basin geology).

## Config
Uses `config/examples/rulison_1969.config`.

## Expected Output
- `output/rulison_1969/*.SAC` -- synthetic seismograms at 3 stations

## Running
```bash
./run.sh
```

## Verified By
- `Integration.HistoricNuclear.Rulison1969` -- deep Mesaverde + explosion + SAC output
