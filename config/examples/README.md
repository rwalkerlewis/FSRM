# FSRM Working Examples

These examples exercise verified code paths and produce physical results.

## uniaxial_compression.config
Quasi-static uniaxial compression of an elastic cube.
Fixed bottom, applied 1mm displacement on top.
Validates elastostatic PetscFE callbacks and Dirichlet BCs.

## fault_compression.config
Same as uniaxial compression but with a locked cohesive fault
at the domain center. Validates fault mesh splitting and
Lagrange multiplier constraints.

## Running
```bash
cd build
./fsrm -c ../config/examples/uniaxial_compression.config
./fsrm -c ../config/examples/fault_compression.config
```
