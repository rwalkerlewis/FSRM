# Example 29: Rio Blanco (1973)

## Physics
Elastodynamic wave propagation through a 5-layer Piceance Basin
sandstone-shale column. Models the Rio Blanco event of 1973-05-17 in
Rio Blanco County, CO -- the third and final Plowshare gas-stimulation
test, originally three 33 kt devices detonated simultaneously at 1780,
1990, and 2040 m depth.

Multi-source simplification: FSRM's [EXPLOSION_SOURCE] section accepts
only one source (the parser uses a single `hasSection("EXPLOSION_SOURCE")`
check at `Simulator.cpp:1162`). The three Rio Blanco devices are
modeled as a single equivalent source with summed yield (99 kt) at the
geometric centroid (1903 m). This loses the vertical-stack moment-tensor
structure but preserves total moment release; the SAC peak amplitude
bound and mb estimate remain valid.

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
| Williams Fork (Mesaverde)  | 1500-2500 | 4500 | 2500 | 2500 |
| Mesaverde Group base       | 2500-5000 | 4200 | 2300 | 2500 |

References: Lombard et al. 1971; Reynolds et al. 1971 (Piceance Basin
Plowshare geology).

## Config
Uses `config/examples/rio_blanco_1973.config`.

## Expected Output
- `output/rio_blanco_1973/*.SAC` -- synthetic seismograms at 3 stations

## Running
```bash
./run.sh
```

## Verified By
- `Integration.HistoricNuclear.RioBlanco1973` -- multi-device single-equivalent + Mesaverde
