# Example 21: Sterling (Project Dribble, 1966)

## Physics
Elastodynamic wave propagation through the Tatum Salt Dome
stratigraphy. Models the 0.38 kt decoupled leg of the Project Dribble
series, detonated 1966-12-03 inside the gas-filled cavity left by
Salmon 1964 at 828 m depth.

Decoupling note: FSRM does not represent the pre-existing Salmon
cavity or its gas filling, so the decoupling factor (~70x relative to
a tamped shot of the same yield, per Patton 1991) is implicit in the
reduced 0.38 kt yield rather than in cavity geometry. A realistic
decoupled-shot model would scale the moment by 1/70 below a corner
frequency set by the cavity radius and roll back to coupled above it;
this example does not do that.

Features used:
- Layered heterogeneous material (aux field assignment)
- Underground nuclear explosion source (Mueller-Murphy, SALT medium)
- Absorbing boundary conditions (Clayton-Engquist)
- Seismometer network with SAC output

## Geology
Identical to Salmon 1964 (same emplacement formation):

| Layer | Depth (m) | Vp (m/s) | Vs (m/s) | Density (kg/m3) |
|-------|-----------|----------|----------|-----------------|
| Surface alluvium / cap soils | 0-200  | 1700 | 900  | 2000 |
| Tertiary sand-shale         | 200-500  | 2700 | 1500 | 2400 |
| Anhydrite cap rock          | 500-700  | 3500 | 2000 | 2500 |
| Tatum dome halite           | 700-2000 | 4500 | 2500 | 2200 |

References: Springer et al. 1968; Patton 1991.

## Config
Uses `config/examples/sterling_1966.config`.

## Expected Output
- `output/sterling_1966/*.SAC` -- synthetic seismograms at 3 stations

## Running
```bash
./run.sh
```

## Verified By
- `Integration.HistoricNuclear.Sterling1966` -- layered salt + small-yield decoupled shot
