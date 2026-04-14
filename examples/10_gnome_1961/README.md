# Example 10: Project Gnome (1961)

## Physics
Elastodynamic wave propagation through a 4-layer evaporite geology model
representing the Delaware Basin near Carlsbad, New Mexico. Models the 3.1 kt
Gnome underground nuclear test, the first US Plowshare experiment, detonated
at 361m depth in the Salado Salt Formation.

Features used:
- Layered heterogeneous material (aux field assignment)
- Underground nuclear explosion source (Mueller-Murphy)
- Absorbing boundary conditions (Clayton-Engquist)
- Seismometer network with SAC output

## Geology
| Layer | Depth (m) | Vp (m/s) | Vs (m/s) | Density (kg/m3) |
|-------|-----------|----------|----------|-----------------|
| Alluvium/Gatuna | 0-50 | 2000 | 1000 | 1800 |
| Rustler Anhydrite | 50-200 | 5000 | 2900 | 2800 |
| Salado Salt | 200-600 | 4500 | 2500 | 2200 |
| Castile Anhydrite | 600-5000 | 5500 | 3200 | 2900 |

## Config
Uses `config/examples/gnome_1961.config`.

## Expected Output
- `output/gnome_1961/*.SAC` -- synthetic seismograms at 3 stations

## Running
```bash
./run.sh
```

## Verified By
- `Integration.HistoricNuclear.Gnome1961` -- layered + explosion + SAC output
