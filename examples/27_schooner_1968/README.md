# Example 27: Schooner (1968)

## Physics
Elastodynamic wave propagation through a 4-layer Pahute Mesa column
(welded tuff over basalt over bedded tuff). Models the 30 kt Schooner
event of 1968-12-08 at 111 m depth -- the final Plowshare cratering
shot, well into the cratering regime (scaled depth-of-burial
~36 m/kt^(1/3.4)).

Features used:
- Layered heterogeneous material (aux field assignment)
- Underground nuclear explosion source (Mueller-Murphy, TUFF medium)
- Absorbing boundary conditions (Clayton-Engquist)
- Seismometer network with SAC output

## Geology
| Layer | Depth (m) | Vp (m/s) | Vs (m/s) | Density (kg/m3) |
|-------|-----------|----------|----------|-----------------|
| Surface alluvium / weathered tuff | 0-50    | 1800 | 900  | 2000 |
| Welded tuff                       | 50-300  | 3500 | 2000 | 2300 |
| Basalt                            | 300-1000 | 4500 | 2500 | 2600 |
| Bedded tuff                       | 1000-2000 | 3300 | 1800 | 2300 |

References: Knox & Terhune 1965 (cratering scaling); Glasstone & Dolan
1977 ch. 6 (Plowshare cratering); Closmann 1969 (NTS tuff/basalt
cavity geology).

## Config
Uses `config/examples/schooner_1968.config`.

## Expected Output
- `output/schooner_1968/*.SAC` -- synthetic seismograms at 3 stations

## Running
```bash
./run.sh
```

## Verified By
- `Integration.HistoricNuclear.Schooner1968` -- shallow tuff/basalt + explosion + SAC output
