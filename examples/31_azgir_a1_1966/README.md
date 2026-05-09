# Example 31: Azgir A-1 (1966)

## Physics
Elastodynamic wave propagation through a 3-layer Caspian salt-dome
column. Models the 1.1 kt Azgir A-1 event of 1966-04-22 at 165 m depth
in halite -- the first Soviet salt-cavity test and the basis of the
subsequent Soviet decoupling-experiment series.

Features used:
- Layered heterogeneous material (aux field assignment)
- Underground nuclear explosion source (Mueller-Murphy, SALT medium)
- Absorbing boundary conditions (Clayton-Engquist)
- Seismometer network with SAC output

## Geology
| Layer | Depth (m) | Vp (m/s) | Vs (m/s) | Density (kg/m3) |
|-------|-----------|----------|----------|-----------------|
| Surface alluvium / clays | 0-50    | 1900 | 900  | 2000 |
| Anhydrite cap rock       | 50-150  | 5000 | 2900 | 2700 |
| Halite (salt dome)       | 150-2000 | 4500 | 2500 | 2200 |

References: Sultanov et al. 1999; Adushkin & Spivak 2015 (Caspian
salt-dome geology and Soviet PNE catalog).

## Config
Uses `config/examples/azgir_a1_1966.config`.

## Expected Output
- `output/azgir_a1_1966/*.SAC` -- synthetic seismograms at 3 stations

## Running
```bash
./run.sh
```

## Verified By
- `Integration.HistoricNuclear.AzgirA1_1966` -- Soviet salt-dome + SAC output
