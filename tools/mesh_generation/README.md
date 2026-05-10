# Source-ball mesh generation (axis-1b)

Pass-13a (foundation slice) ships the `build_source_ball_mesh.py`
TetGen-driven mesh generator for the axis-1b 3D source ball. The C++
runtime consumes the generated `.node` / `.ele` files via
`Source3DBallMesh`; this directory is for the build-host pre-process
step that produces them.

## Dependency

`build_source_ball_mesh.py` shells out to **TetGen** (Si 2015,
BSD-licensed), which must be on `PATH`. TetGen is not a
build-time CMake dependency and is not in the `fsrm-ci:local`
Docker image. Install it on the build host (the developer workstation
or a CI step that owns mesh regeneration):

```
apt-get install tetgen        # Debian / Ubuntu
brew install tetgen           # macOS
```

or build from source: <https://www.wias-berlin.de/software/tetgen>.

The runtime never invokes TetGen. The Python tool runs once per
mesh-config combination and the resulting `.node` / `.ele` pair is
cached under `cache/source_ball_meshes/<event>_<resolution>/`.

## Usage

```
python3 tools/mesh_generation/build_source_ball_mesh.py \
    --cavity-radius 17.0 \
    --elastic-radius 70.0 \
    --edge-near 5.0 \
    --edge-outer 12.0 \
    --output cache/source_ball_meshes/salmon_default
```

This produces `cache/source_ball_meshes/salmon_default.node` and
`cache/source_ball_meshes/salmon_default.ele`. The `.poly` PSLG that
TetGen consumes is removed unless `--keep-poly` is passed.

Arguments:

| Flag | Meaning |
|------|---------|
| `--cavity-radius` | Inner sphere radius `R_c` in metres. |
| `--elastic-radius` | Outer sphere radius `r_elastic` in metres. |
| `--edge-near` | Target tet edge length on the cavity surface. |
| `--edge-outer` | Target tet edge length on the outer surface. |
| `--quality` | TetGen radius-edge ratio (default `1.4`). |
| `--output` | Output basename (without extension). |
| `--keep-poly` | Preserve the intermediate `.poly` file (debugging). |

## Cache layout

```
cache/source_ball_meshes/
    salmon_default.node           # ~50000 tets at edge_outer=5 m, edge_near=2 m
    salmon_default.ele
    salmon_smoke.node             # ~5000 tets, faster smoke-build
    salmon_smoke.ele
    chagan_default.node
    chagan_default.ele
    ...
```

Each meshfile pair is keyed by event + resolution. Pass-13a does not
yet read this cache from a runtime config field; that wires up in
pass-13b once the hydro substep needs the mesh at production
resolution. Foundation tests use the smaller fixtures committed under
`tests/data/source_ball/`.

## Geometry

Two concentric icospheres tessellate the spherical shell. The inner
icosphere (radius `R_c`) bounds the cavity; the outer icosphere
(radius `r_elastic`) is the elastic-radius outer boundary. A
hole marker at the origin removes the cavity interior. TetGen meshes
the resulting PSLG with quality constraints.

Surface markers in the emitted `.node` file:

| Marker | Surface |
|--------|---------|
| `1` | Outer (elastic-radius) surface |
| `2` | Inner (cavity) surface |
| `0` | Interior vertex |

`Source3DBallMesh` uses these markers to label the corresponding
DMPlex facet sets (`CAVITY_SURFACE`, `ELASTIC_SURFACE`) for the
surface-integral moment-tensor extraction wired in pass-13c.

## Reference fixture meshes (tests/data)

Foundation unit tests do not require TetGen. They consume small
hand-written `.node` / `.ele` fixtures committed under
`tests/data/source_ball/`. Those fixtures are the reference for the
TetGen output format; the runtime parser is exercised against them
without depending on a TetGen install.

## Reference

- Si, H. (2015), "TetGen, a Delaunay-based quality tetrahedral mesh
  generator", ACM Trans. Math. Software 41(2), Article 11.
