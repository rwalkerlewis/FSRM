# Example 18: Thermal Expansion

## Overview

This example demonstrates thermo-mechanical coupling (THM) through
free thermal expansion of a 3D box subjected to uniform heating.

## Physics

- Quasi-static geomechanics coupled with thermal diffusion
- Thermoelastic stress: sigma = C : (eps - alpha_T * dT * I)
- Uniform temperature increase: dT = 100 K above T_ref = 293 K

## Material Properties

| Property | Value | Unit |
|----------|-------|------|
| Young's modulus | 10 GPa | Pa |
| Poisson's ratio | 0.25 | - |
| Density | 2650 | kg/m^3 |
| Thermal conductivity | 3.0 | W/(m*K) |
| Specific heat | 800 | J/(kg*K) |
| Thermal expansion | 1e-5 | 1/K |

## Expected Results

With free thermal expansion (fixed bottom, free top and sides):
- Vertical strain: eps = alpha_T * dT = 1e-5 * 100 = 1e-3
- Top displacement: uz = eps * Lz = 1e-3 m
- Stress: approximately zero (no mechanical constraint)

## Running

```bash
cd ../../build
./fsrm -c ../examples/18_thermal_expansion/config.config \
    -ts_type beuler -pc_type lu -ksp_type preonly
```
