# FSRM Working Examples

These examples exercise verified code paths and produce physical results.
All examples build and run in the Docker CI environment.

## Elastostatics

### uniaxial_compression.config
Quasi-static uniaxial compression of an elastic cube.
Fixed bottom, applied 1mm displacement on top.
Validates elastostatic PetscFE callbacks and Dirichlet BCs.

### fault_compression.config
Same as uniaxial compression but with a locked cohesive fault
at the domain center. Validates fault mesh splitting and
Lagrange multiplier constraints.

## Poroelasticity

### terzaghi_consolidation.config
1D Terzaghi consolidation: drained top, impermeable sides, fixed bottom.
Couples single-phase fluid flow with linear elasticity via Biot poroelasticity.
Field 0 = pressure, Field 1 = displacement. SNES converges in 1-2 iterations.

### injection_pressure_buildup.config
Poroelastic simulation with a point injection source at the domain center.
Demonstrates pressure buildup from fluid injection (Q = 0.001 m³/s).
Validates injection cell location via DMLocatePoints and residual modification.

## Hydraulic Fracture

### hydraulic_fracture_pkn.config
Standalone PKN (Perkins-Kern-Nordgren) hydraulic fracture propagation.
Demonstrates analytical fracture width, length, and pressure evolution.
Uses the HydraulicFractureModel with PKN geometry.

## Cohesive Fracture

### cohesive_hydraulic_fracture.config
Poroelastic simulation with a pre-defined cohesive fracture plane.
Uses Lagrange multiplier constraints with tensile failure criterion.
**Note**: SNES diverges due to known PETSc hybrid cell assembly limitations.
Mesh surgery, field setup, and boundary conditions all work correctly.

## Explosion Source

### explosion_seismogram.config
Underground nuclear explosion (10 kt, 300 m depth) with 3 virtual seismometers.
Elastodynamics with TSCN timestepper and isotropic moment tensor source.
Uses Mueller-Murphy (1971) source model with cavity mechanics.

### dprk_2017_quick.config
DPRK 2017 Punggye-ri nuclear test (250 kt, 800 m depth in granite).
Quick near-field simulation (5 seconds). 20x20x20 mesh, 20 km domain.
Expected mb ~ 6.25 (observed USGS mb 6.3).

### underground_explosion_template.config
Comprehensive template with detailed comments for customizing underground
explosion simulations. Includes material property tables and scaling rules.

## Running

```bash
cd build

# Elastostatics
./fsrm -c ../config/examples/uniaxial_compression.config -snes_monitor -pc_type lu

# Poroelasticity
./fsrm -c ../config/examples/terzaghi_consolidation.config -ts_monitor -snes_monitor -pc_type lu
./fsrm -c ../config/examples/injection_pressure_buildup.config -ts_monitor -snes_monitor -pc_type lu

# Explosion
./fsrm -c ../config/examples/explosion_seismogram.config -ts_monitor -snes_monitor -pc_type lu
./fsrm -c ../config/examples/dprk_2017_quick.config -ts_monitor -snes_monitor -pc_type lu
```
