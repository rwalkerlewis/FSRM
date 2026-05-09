# Example 11: Sedan Crater (1962)

## Physics
Elastodynamic wave propagation through a 3-layer geology model representing
NTS Yucca Flat. Models the 104 kt Sedan underground nuclear test from
Operation Storax/Plowshare, detonated at only 194m depth in alluvium.
This shallow burial produced the largest US nuclear crater (390m diameter,
100m deep).

Features used:
- Layered heterogeneous material (aux field assignment)
- Underground nuclear explosion source (Mueller-Murphy)
- Absorbing boundary conditions (Clayton-Engquist)
- Seismometer network with SAC output

## Geology
| Layer | Depth (m) | Vp (m/s) | Vs (m/s) | Density (kg/m3) |
|-------|-----------|----------|----------|-----------------|
| Alluvium | 0-300 | 1800 | 800 | 1700 |
| Welded Tuff | 300-800 | 4200 | 2400 | 2300 |
| Paleozoic Carbonate | 800-5000 | 5500 | 3100 | 2700 |

## Configs

Two configs ship with this example:

- `config/examples/sedan_1962.config` -- baseline KINEMATIC_RDP path.
- `config/examples/sedan_1962_dynamic.config` -- pass-5 DYNAMIC_PLASTIC
  path; drives the far-field FEM problem from the recorded 1D
  NearFieldExplosionSolver moment-tensor history (full 6 components,
  including CLVD content).

Sedan is the pass-5 anchor event because the cratering shot exercises
near-surface damage and the alluvium medium puts the cavity radius at
a tractable scale (Rc ~ 87 m via the medium-aware NTS coefficient)
relative to the far-field cell scale (Lz / nz = 500 m). See
`docs/HISTORIC_NUCLEAR_FIDELITY.md` "Closed in pass 5" and
`docs/HISTORIC_NUCLEAR_ROADMAP.md` axis 1 for the rationale.

## Expected Output

Baseline (`run.sh`):
- `output/sedan_1962/*.SAC` -- synthetic seismograms at 3 stations

Dynamic-plastic (`run_dynamic.sh`):
- `output/sedan_1962_dynamic/*.SAC` -- synthetic seismograms
- `output/sedan_1962_dynamic/near_field_history.csv` -- recorded 1D
  solver moment-rate tensor and cavity / plastic radius time series

## Running

```bash
# Baseline KINEMATIC_RDP path:
./run.sh

# Pass-5 DYNAMIC_PLASTIC path:
./run_dynamic.sh
```

## Visualization (pass-5)

Three ParaView state files in `paraview/` configure the rendering of
the dynamic-plastic outputs:

- `paraview/near_field_cavity.pvsm` -- XY chart of cavity / plastic
  radius and moment-rate trace over time.
- `paraview/far_field_propagation.pvsm` -- 3D wavefield render of
  `output/solution.h5`.
- `paraview/combined.pvsm` -- multi-view layout of both.

These are minimal hand-authored stubs (the FSRM build environment does
not run ParaView in CI so they cannot be visually verified end-to-end);
see `paraview/README.md` for what to refine after first opening.

## Verified By

- `Integration.HistoricNuclear.Sedan1962` -- baseline layered + explosion + SAC output
- `Integration.HistoricNuclear.Sedan1962_Distributed` -- pass-4 multi-cell
  moment-tensor distribution at the factor-30 envelope
- `Integration.HistoricNuclear.Sedan1962_Dynamic` -- pass-5
  DYNAMIC_PLASTIC source: cavity radius within 20% of medium-aware
  NTS analytic, peak/u_far within factor 30, near_field_history.csv
  emitted with > 100 sample rows
