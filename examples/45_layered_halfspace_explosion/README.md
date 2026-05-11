# Example 45: Layered halfspace explosion (verification anchor)

Synthetic isotropic explosion in a 4-layer Earth model. Pass-13c
academic verification anchor; the other is
`examples/44_lambs_problem/`.

## Physics

Wave propagation through a vertically stratified medium with a free
surface and a buried isotropic explosion source. The wavefield
contains:

1. Direct P / S body waves (refracted at layer interfaces).
2. Layer-bounce reflections (PmP, SmS, sP, etc.).
3. Surface waves: Rayleigh phase (Rg) dominated by the sediment
   waveguide; Love phase if the source is non-isotropic (not in this
   example).

The Haskell-Thomson propagator-matrix method (Haskell 1953, Thomson
1950) gives the dispersion relation, group velocities, and arrival
times for each phase. Modern exposition: Aki & Richards 2002 ch 7.

## Geology

| Layer | z_top (m) | z_bot (m) | Vp (m/s) | Vs (m/s) | rho (kg/m^3) |
|------:|----------:|----------:|---------:|---------:|-------------:|
| 1     |    60000  |    58000  |   3000   |   1700   |     2300     |
| 2     |    58000  |    45000  |   5500   |   3200   |     2700     |
| 3     |    45000  |    25000  |   6600   |   3700   |     2900     |
| 4     |    25000  |        0  |   8000   |   4500   |     3300     |

(Mesh z = 60000 is the free surface; z = 0 is the absorbing
bottom.)

## Validation

Headline gate: Rg arrival time at 25 km surface range matches the
Haskell-Thomson group velocity v_Rg ~ 3.0 km/s (sediment-dominant Rg)
within 5 percent.

Expected Rg arrival at 25 km: t_Rg ~ 8.3 s.
Gate envelope: t_observed in [7.9, 8.7] s.

## Pass-13c wavefield output

`[OUTPUT] wavefield_format = HDF5_XDMF` emits the full ParaView
animation. The 100-step cadence at the 30-second time horizon yields
~ 40 frames; opening `output/wavefield.xdmf` in ParaView shows the
P-wave / S-wave / Rg arrival sequence with the layer interfaces
visible as reflection / refraction discontinuities.

## How to run

```bash
./run.sh
```

Override rank count with `MPI_RANKS=8 ./run.sh`.

## References

- Haskell, N. A. (1953), "The dispersion of surface waves on
  multilayered media", Bull. Seism. Soc. Am. 43, pp 17-34.
- Thomson, W. T. (1950), "Transmission of elastic waves through a
  stratified solid medium", J. Appl. Phys. 21, pp 89-93.
- Aki, K. and Richards, P. G. (2002), "Quantitative Seismology"
  2nd ed, ch 7.
