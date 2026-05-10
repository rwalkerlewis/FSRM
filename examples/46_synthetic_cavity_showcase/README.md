# Example 46: synthetic cavity-formation showcase

Pass-13c showcase example designed for cinematic ParaView animation.
NOT a verification anchor — see `examples/44_lambs_problem/` and
`examples/45_layered_halfspace_explosion/` for the closed-form
verification anchors that pass-13c ships.

## Physics

Synthetic 1 kt explosion in a homogeneous granite halfspace at 1500 m
depth. The depth is chosen to balance two competing requirements:

- Deep enough that the cavity-formation phase finishes before
  significant interaction with the free surface (so the source-region
  animation isn't contaminated by reflected waves).
- Shallow enough that the radiated wavefield reaches the surface
  within the simulation window so the animation shows the surface-wave
  development.

## Geometry

- 8 km x 8 km x 6 km domain.
- Free surface on top; absorbing BCs on sides and bottom.
- Source at (4000, 4000, 4500): centred horizontally, 1500 m below
  the free surface.
- Three surface receivers along +x at 500, 1500, 3000 m offset.

## Resolution

The shipped config uses 32 x 32 x 24 = 24576 elements for CI
tractability. For the full cinematic render, override the GRID block:

```ini
[GRID]
nx = 64
ny = 64
nz = 48
```

That gives ~196608 elements and ~10 minutes wall-clock at MPI=4 on a
modern workstation.

For "production" cinematic quality (~500 k elements in the source
region, ~2 M total per the pass-13c spec), use the staged refinement
script:

```
scripts/refine_mesh.py --base 64 --refine 3
```

(Production-grade output is intentionally outside CI scope.)

## Wavefield output

`[OUTPUT] wavefield_format = HDF5_XDMF` with cadence 25 steps gives
~320 animation frames at the default time horizon (4 s with
~ 0.0005 s `dt_initial`). The XDMF wrapper is refreshed at every
write so partial runs are still ParaView-readable.

The shipped `paraview/wavefield_3view.pvsm` ParaView 5.10+ state
configures three render views:

1. **Source-region** (top-left): cavity-formation animation
   centred at (4000, 4000, 4500) with displacement-magnitude pseudo-
   colour, Viridis colormap.
2. **Far-field** (right): full domain volume render from a side angle
   showing the radiated wavefield.
3. **Surface** (bottom-left): surface displacement contours at z = +Lz
   (top) showing the Rayleigh-wave pattern.

## Validation

Smoke-only: pipeline completes; non-zero solution norm; SAC files
written. No closed-form gate (see Lamb / layered-halfspace examples).

## How to run

```bash
./run.sh
```

Override rank count with `MPI_RANKS=8 ./run.sh`. For the production-
resolution cinematic render, edit the GRID block first.

## References

- Stevens, J. L. and Day, S. M. (1985), "The physical basis of mb:Ms
  and variable frequency magnitude methods", J. Geophys. Res. 90.
- Patton, H. J. (1991), "Decoupling and topological factors at NTS",
  LANL technical report.
