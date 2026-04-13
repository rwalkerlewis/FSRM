# Example 03: Elastoplastic Compression

## Physics

Quasi-static compression of a granite column with Drucker-Prager plasticity.

A 1x1x2 m granite column (E=50 GPa, nu=0.25) with Drucker-Prager yield
parameters (cohesion=5 MPa, friction angle=30 deg) is compressed axially.
Under oedometric conditions (roller sides), the effective modulus M = lambda+2mu
produces stress sigma = M*eps that exceeds the yield stress (~17.3 MPa),
demonstrating plastic deformation.

Uses:
- PetscFEElastoplasticity f1/g3 callbacks
- PlasticityModel return mapping (Drucker-Prager)
- Unified constants indices 32-35 for plasticity parameters

## Config

Uses `config/examples/elastoplastic_compression.config`.

## Expected Output

- `output/solution.h5` -- HDF5 solution file with displacement field

The stress state should exceed the Drucker-Prager yield surface, producing
permanent plastic deformation. Below yield, the response matches linear elasticity.

## Running

```bash
./run.sh
```

Or manually:
```bash
cd /path/to/build
./fsrm -c ../config/examples/elastoplastic_compression.config
```

## Verified By

- `Integration.ElastoplasticSim` -- quasi-static beyond/below yield
- `Physics.ElastoplasticBearingCapacity` -- below-yield matches elastic
- `Unit.DruckerPragerStandalone` -- return mapping correctness
