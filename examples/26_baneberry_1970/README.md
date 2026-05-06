# Example 26: Baneberry (1970)

## Physics
Elastodynamic wave propagation through a 4-layer Yucca Flat column
including a saturated tuff zone. Models the 10 kt Baneberry event of
1970-12-18 at 278 m depth in NTS Area 8 (U-8d), notable for venting
radioactive gas to the surface within minutes of detonation.

Saturation note: FSRM does not represent pore-water saturation
explicitly; the lowered Vs/Vp values in LAYER_2 approximate the
saturated mechanical state inferred by Carothers et al. 1995.

Features used:
- Layered heterogeneous material (aux field assignment)
- Underground nuclear explosion source (Mueller-Murphy, TUFF medium)
- Absorbing boundary conditions (Clayton-Engquist)
- Seismometer network with SAC output

## Geology
| Layer | Depth (m) | Vp (m/s) | Vs (m/s) | Density (kg/m3) |
|-------|-----------|----------|----------|-----------------|
| Alluvium                | 0-100   | 1500 | 700  | 1800 |
| Saturated tuff          | 100-400 | 2200 | 1100 | 2100 |
| Welded tuff             | 400-1000 | 3500 | 1900 | 2300 |
| Paleozoic basement      | 1000-2000 | 5000 | 2900 | 2700 |

Reference: Carothers et al. 1995 (Baneberry venting analysis).

## Config
Uses `config/examples/baneberry_1970.config`.

## Expected Output
- `output/baneberry_1970/*.SAC` -- synthetic seismograms at 3 stations

## Running
```bash
./run.sh
```

## Verified By
- `Integration.HistoricNuclear.Baneberry1970` -- saturated tuff + explosion + SAC output
