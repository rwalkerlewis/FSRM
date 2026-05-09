# Example 32: Pokhran-I (Smiling Buddha, 1974)

## Physics
Elastodynamic wave propagation through a 3-layer Trans-Aravalli
sandstone/granite column. Models the 8 kt declared-yield Pokhran-I
event of 1974-05-18 at 107 m depth in Vindhyan sandstone -- the first
Indian nuclear test.

Features used:
- Layered heterogeneous material (aux field assignment)
- Underground nuclear explosion source (Mueller-Murphy, SHALE medium)
- Absorbing boundary conditions (Clayton-Engquist)
- Seismometer network with SAC output

## Geology
| Layer | Depth (m) | Vp (m/s) | Vs (m/s) | Density (kg/m3) |
|-------|-----------|----------|----------|-----------------|
| Aeolian sand / surface dunes  | 0-30   | 1500 | 700  | 1800 |
| Vindhyan sandstone            | 30-200 | 3500 | 1900 | 2400 |
| Trans-Aravalli granite/gneiss | 200-2000 | 5500 | 3100 | 2700 |

Reference: Sikka et al. 2000 (Pokhran Test Range geology).

## Config
Uses `config/examples/pokhran_i_1974.config`.

## Expected Output
- `output/pokhran_i_1974/*.SAC` -- synthetic seismograms at 3 stations

## Running
```bash
./run.sh
```

## Verified By
- `Integration.HistoricNuclear.PokhranI_1974` -- Vindhyan sandstone + SAC output
