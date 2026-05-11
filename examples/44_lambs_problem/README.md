# Example 44: Lamb's problem (verification anchor)

Vertical step force on the surface of an elastic halfspace
(Lamb 1904). One of two academic verification anchors that pass-13c
adds with closed-form reference solutions; the other is
`examples/45_layered_halfspace_explosion/`.

## Physics

Free-surface response to a localized vertical force at the surface of
an elastic halfspace. The radiated wavefield contains:

1. P-wave: longitudinal body wave, geometric spreading 1/r.
2. S-wave: transverse body wave, geometric spreading 1/r.
3. Rayleigh wave: free-surface wave, geometric spreading 1/sqrt(r),
   dominant at far range.

The classical solution is given in Lamb 1904 (Phil. Trans. Roy. Soc.
London A 203). Eringen & Suhubi 1975 vol II ch 7 gives the closed-form
expressions in modern notation; Mooney 1974 (Bull. Seism. Soc. Am.
64) provides numerical tables.

## Geometry

- 5 km cubic domain.
- Free surface on top (z = +Lz). Absorbing BCs on the sides and
  bottom emulate an infinite halfspace.
- Source at the centre of the top surface, depth-of-burial 50 m to
  approximate a surface-localised force.
- Three surface receivers at ranges 500, 1000, 1500 m along the +x
  axis.

## Validation

The relevant pass-13c gate is the Rayleigh-wave 1/sqrt(r) geometric-
spreading scaling: vertical displacement peak amplitudes at the three
receivers should obey
    |u_z(r1)| / |u_z(r2)|  ~  sqrt(r2 / r1)
within factor 2 (CI mesh resolution).

The full closed-form Eringen-Suhubi expressions are not implemented in
CI; users can compare against Mooney 1974 tables in any post-processing
analysis.

## Pass-13c wavefield output

`[OUTPUT] wavefield_format = HDF5_XDMF` emits a ParaView-readable
time series at `wavefield_cadence_steps = 50`. Open
`output/wavefield.xdmf` in ParaView 5.10+ to animate the wave fronts.
The `paraview/wavefield.pvsm` state file ships a pre-configured view
with displacement-magnitude colormap (Viridis), surface-wave time-keeper,
and a side camera angle to show the P-wave / S-wave / Rayleigh-wave
arrival sequence.

## How to run

```bash
./run.sh
```

Override rank count with `MPI_RANKS=8 ./run.sh`.

## References

- Lamb, H. (1904), "On the propagation of tremors over the surface of
  an elastic solid", Phil. Trans. Roy. Soc. London A 203, pp 1-42.
- Eringen, A. C. and Suhubi, E. S. (1975), "Elastodynamics" vol II,
  Academic Press, ch 7.
- Mooney, H. M. (1974), "Some numerical solutions for Lamb's problem",
  Bull. Seism. Soc. Am. 64, pp 473-491.
