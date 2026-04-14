# Example 14: Single-Phase Pressure Diffusion

## Physics

This example demonstrates single-phase Darcy flow in a porous medium. The governing equation is the pressure diffusion equation:

    phi * ct * dP/dt = div( (k/mu) * grad(P) )

where:
- phi = 0.20 (porosity)
- ct = 1e-9 1/Pa (total compressibility)
- k = 100 mD = 9.87e-14 m^2 (permeability)
- mu = 0.001 Pa*s (water viscosity)

## Setup

- Domain: 100m x 1m x 1m (effectively 1D)
- Grid: 40 x 1 x 1 hexahedral cells
- Boundary conditions:
  - x_min: Dirichlet pressure = 20 MPa
  - x_max: Dirichlet pressure = 10 MPa
  - All other faces: no-flow (natural Neumann)
- Initial condition: zero pressure (cold start)
- End time: 200 s

## Expected Results

At steady state, the pressure profile is linear:

    P(x) = 20e6 - (10e6 / 100) * x

The center pressure should approach 15 MPa. The characteristic diffusion time is:

    t_diff = phi * mu * ct * L^2 / k
           = 0.2 * 0.001 * 1e-9 * 100^2 / 9.87e-14
           = 20.3 s

So after 200 s (about 10 diffusion times), the solution should be very close to steady state.

## Running

```bash
docker run --rm -v $(pwd)/../..:/workspace -w /workspace/build fsrm-ci:local \
  ./fsrm -c ../examples/14_single_phase_flow/config.config
```
