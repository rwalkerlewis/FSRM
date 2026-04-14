# Example 16: SCEC TPV5 Dynamic Rupture Benchmark

## Overview

SCEC TPV5 is a community benchmark for dynamic rupture simulations, defined by the
Southern California Earthquake Center (SCEC) Code Verification for Wave-Propagation
Simulations (CVWS) working group.

Reference: https://strike.scec.org/cvws/cgi-bin/cvws.cgi

## Setup

- **Geometry**: 32 km x 100 m x 17 km domain with vertical strike-slip fault at center
- **Material**: Homogeneous elastic (rho=2670 kg/m^3, Vs=3464 m/s, Vp=6000 m/s)
- **Friction**: Linear slip-weakening (mu_s=0.677, mu_d=0.525, Dc=0.40 m)
- **Initial stress**: sigma_n = -120 MPa, tau = 70 MPa (background)
- **Nucleation**: Overstressed patch (tau=81.6 MPa) at center, radius 1500 m
- **Time stepping**: TSALPHA2, dt=0.002 s, 12 s total

## Running

```bash
cd /workspace/build
./fsrm -c ../examples/16_scec_tpv5/config.config \
  -ts_type alpha2 \
  -snes_max_it 50 \
  -pc_type lu \
  -ksp_type preonly
```

## Expected Behavior

The nucleation patch has shear stress (81.6 MPa) exceeding the static friction
strength (mu_s * |sigma_n| = 0.677 * 120 MPa = 81.24 MPa), so rupture nucleates
immediately and propagates bilaterally along the fault. Slip weakening reduces
friction to mu_d * |sigma_n| = 0.525 * 120 MPa = 63 MPa, allowing dynamic
rupture propagation.
