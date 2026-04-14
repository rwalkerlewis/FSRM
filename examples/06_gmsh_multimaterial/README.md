# Example 06: Gmsh Multi-Material

## Physics

Elastodynamic wave propagation through a Gmsh-imported mesh with per-region
material properties and an underground explosion source.

The mesh (`meshes/test_two_material.msh`) contains two labeled regions:
- **soft_rock**: E=12 GPa, nu=0.25, rho=2200 kg/m^3
- **hard_rock**: E=75 GPa, nu=0.25, rho=2800 kg/m^3

Material properties are assigned per cell based on Gmsh physical labels using
the auxiliary field pattern. A 5 kt explosion source generates waves that
interact with the material interface.

Uses:
- GmshIO for MSH2 mesh import with physical-name labels
- Per-cell material assignment via `[MATERIAL_REGION_N]` config sections
- Auxiliary field population from gmsh_label mapping
- Mueller-Murphy source model
- Damage-zone material degradation near the explosion cavity

## Config

Uses `config/examples/nuclear_twin_gmsh.config`.
Requires mesh: `meshes/test_two_material.msh`.

## Expected Output

- `output/nuclear_twin_gmsh/solution.h5` -- HDF5 solution snapshots

The wavefield should show different wave speeds in the soft and hard rock
regions, with reflections and refractions at the material interface.

## Running

```bash
./run.sh
```

Or manually:
```bash
cd /path/to/build
./fsrm -c ../config/examples/nuclear_twin_gmsh.config
```

## Verified By

- `Integration.NuclearTwinGmsh` -- mapped materials + explosion + HDF5
- `Integration.GmshMultiMaterial` -- per-label lambda/mu/rho assignment
- `Integration.GmshImport` -- MSH2 physical names, tet cells
