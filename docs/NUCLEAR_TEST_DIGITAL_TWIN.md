# Nuclear Test Digital Twin Workflow

This document describes a verified workflow for running an underground nuclear test digital twin in FSRM using layered materials, absorbing boundaries, and seismometer outputs.

## Scope

This workflow is validated by automated tests in Docker, including:

- `Integration.PunggyeRiLayered`
- `Integration.ExplosionSeismogram`
- `Integration.GmshImport`
- `Integration.GmshMultiMaterial`
- `Integration.GasbuggyMesh`
- `Integration.NuclearTwinGmsh`
- `Physics.ExplosionDamageZone`

## Prerequisites

- Docker image built from `Dockerfile.ci`
- Existing build directory at `build/`
- Example config files in `config/examples/`

Build and test baseline:

```bash
docker build -f Dockerfile.ci -t fsrm-ci:local .
docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:local bash -c \
  'mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON -DENABLE_CUDA=OFF -DBUILD_EXAMPLES=ON && make -j$(nproc)'

docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:local \
  ctest --output-on-failure
```

## Path A: Layered Punggye-ri style model

Use the compact layered example:

- `config/examples/punggye_ri_layered.config`

Run:

```bash
docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:local \
  ./fsrm ../config/examples/punggye_ri_layered.config
```

Expected behavior:

- HDF5 solution snapshots are written to the configured output directory.
- SAC seismograms are produced for listed stations.
- Absorbing boundaries are active on enabled faces.

Important boundary condition requirement:

- If a face uses absorbing boundary conditions, that face must not use a Dirichlet constraint at the same time.
- In practice for these explosion examples, set `bottom = free` and `sides = free` when `z_min`, `x_min/x_max`, `y_min/y_max` are enabled under `[ABSORBING_BC]`.

## Path B: Gmsh multi-material model

FSRM supports per-cell auxiliary material properties from Gmsh physical labels via:

- `[MATERIAL] assignment = gmsh_label`
- `[MATERIAL_REGION_N] gmsh_label = <physical-name>`

Ready-to-run example config:

- `config/examples/nuclear_twin_gmsh.config`

Validated examples and tests:

- `tests/integration/test_gmsh_multimaterial.cpp` with `meshes/test_two_material.msh`
- `tests/integration/test_gasbuggy_mesh.cpp` with `meshes/historical/gasbuggy_layered_3d.msh`

Run targeted checks:

```bash
docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:local \
  ctest --output-on-failure -R 'Integration.GmshMultiMaterial|Integration.GasbuggyMesh'
```

## Path C: Combined Gmsh nuclear twin workflow

This path combines mapped Gmsh material regions with an underground source and HDF5 output in one integration scenario.

Use:

- `config/examples/nuclear_twin_gmsh.config`

Run example:

```bash
docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:local \
  ./fsrm ../config/examples/nuclear_twin_gmsh.config
```

Run validation test:

```bash
docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:local \
  ctest --output-on-failure -R Integration.NuclearTwinGmsh
```

## Damage-zone coupling in FEM run

Explosion near-field degradation of auxiliary material properties is enabled with:

- `[EXPLOSION_SOURCE] apply_damage_zone = true`

Validated by:

- `Physics.ExplosionDamageZone`

Run targeted check:

```bash
docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:local \
  ctest --output-on-failure -R Physics.ExplosionDamageZone
```

## Output validation checklist

After a successful run, verify:

1. HDF5 output file exists and is non-empty.
2. SAC files exist for all configured stations and components.
3. Seismogram amplitudes are non-zero.
4. Station amplitude decay with distance is physically plausible for selected components.

## Notes on current verified status

- The current Docker baseline is 91 passing tests.
- The Punggye-ri layered pipeline is integration-tested with HDF5 and SAC outputs.
- Gmsh label based per-cell material assignment is integration-tested.
- Historical Gasbuggy mesh mapping is integration-tested.
- Combined Gmsh nuclear twin workflow is integration-tested.
- Explosion damage-zone degradation is physics-validated in the FEM pipeline.
