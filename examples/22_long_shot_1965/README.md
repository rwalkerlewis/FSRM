# Example 22: Long Shot (1965)

## Physics
Elastodynamic wave propagation through a 4-layer Aleutian-arc
andesite/breccia stratigraphy at Amchitka Island, AK. Models the 80 kt
Long Shot event of 1965-10-29 at 716 m depth, the first of three
Amchitka shots (Long Shot 1965, Milrow 1969, Cannikin 1971).

Medium choice: GENERIC. The cavity-coefficient table in
`ExplosionImpactPhysics.cpp` does not include andesite; GENERIC
(12 m/kt^(1/3)) is used as a small upward bias acknowledging the
volcanic-arc rheology.

Features used:
- Layered heterogeneous material (aux field assignment)
- Underground nuclear explosion source (Mueller-Murphy, GENERIC medium)
- Absorbing boundary conditions (Clayton-Engquist)
- Seismometer network with SAC output

## Geology
| Layer | Depth (m) | Vp (m/s) | Vs (m/s) | Density (kg/m3) |
|-------|-----------|----------|----------|-----------------|
| Tundra / weathered overburden | 0-50    | 1500 | 750  | 1800 |
| Weathered andesite           | 50-300  | 2800 | 1500 | 2300 |
| Andesite breccia             | 300-800 | 4000 | 2200 | 2500 |
| Competent andesite           | 800-2000 | 4800 | 2700 | 2700 |

Reference: King, Foster, & Bingham 1972 (Amchitka).

## Config
Uses `config/examples/long_shot_1965.config`.

## Expected Output
- `output/long_shot_1965/*.SAC` -- synthetic seismograms at 3 stations

## Running
```bash
./run.sh
```

## Verified By
- `Integration.HistoricNuclear.LongShot1965` -- layered andesite + explosion + SAC output
