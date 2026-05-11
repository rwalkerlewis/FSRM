# salmon_default 3D source-ball mesh fixture

Two-icosphere shell mesh for the axis-1b 3D source ball, sized for the
Salmon 1964 shot (5.3 kt, Tatum Salt Dome):

- cavity surface radius `R_cavity = 17.0 m`
- elastic surface radius `R_elastic = 70.0 m`
- 935 vertices, ~2972 tetrahedra
- vertex marker convention: `1 = outer elastic surface`,
  `2 = inner cavity surface`, `0 = interior`

Regenerate (TetGen required on the build host, not in CI) via:

```
python3 tools/mesh_generation/build_source_ball_mesh.py \
    --cavity-radius 17.0 --elastic-radius 70.0 \
    --edge-near 5.0 --edge-outer 12.0 \
    --output cache/source_ball_meshes/salmon_default
```

then copy `salmon_default.{node,ele}` to `source_ball.{node,ele}` here.

Shipped under `tests/data/` (the convention for committed mesh
fixtures, alongside `cube_5tet/` and `single_tet/`) so the pass-14a
`HistoricNuclearTest.Salmon3D_WithOverburdenAtMPI4` end-to-end gate
runs in CI without depending on the gitignored
`cache/source_ball_meshes/` directory. The mesh is deterministic
(TetGen on a fixed icosphere tessellation), so re-running the tool
produces the same file.
