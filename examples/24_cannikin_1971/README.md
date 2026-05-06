# Example 24: Cannikin (1971)

## Physics
Elastodynamic wave propagation through a 5-layer Aleutian-arc andesite
column at Amchitka Island, AK. Models the up-to-5 Mt Cannikin event of
1971-11-06 at 1860 m depth -- the largest US underground nuclear test
and the deepest in the Amchitka series. Triggered an mb 6.97 seismic
event recorded teleseismically.

Medium choice: GENERIC (no andesite coefficient in the cavity table).

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
| Andesite breccia             | 300-1200 | 4000 | 2200 | 2500 |
| Competent andesite           | 1200-3000 | 4800 | 2700 | 2700 |
| Deep andesite                | 3000-5000 | 5500 | 3100 | 2800 |

References: King, Foster, & Bingham 1972 (Amchitka); Lay, Wallace, &
Helmberger 1984 (Cannikin teleseismic mb).

## Config
Uses `config/examples/cannikin_1971.config`.

## Expected Output
- `output/cannikin_1971/*.SAC` -- synthetic seismograms at 3 stations

## Running
```bash
./run.sh
```

## Verified By
- `Integration.HistoricNuclear.Cannikin1971` -- multi-megaton andesite + explosion + SAC output
