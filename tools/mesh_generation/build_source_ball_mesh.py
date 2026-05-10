#!/usr/bin/env python3
"""
build_source_ball_mesh.py
=========================

Pass-13a (axis-1b foundation) TetGen-driven mesh generator for the 3D
source ball used by the axis-1b solver. Produces TetGen `.node` and
`.ele` files that the C++ runtime consumes via Source3DBallMesh.

Architecture decision (pass-13a, see docs/AXIS_1B_DESIGN.md):
    TetGen runs as a build-host pre-process tool, not as a build- or
    runtime CMake dependency. The C++ runtime parses TetGen's plain-
    text output. This keeps fsrm-ci's Docker image free of TetGen and
    keeps the meshes deterministic across CI runs (the mesh files are
    cached under cache/source_ball_meshes/ and reused).

Geometry:
    Two concentric icospheres tessellate the source-ball domain.
        inner: r = R_cavity              (approximate cavity surface)
        outer: r = R_elastic             (elastic-radius outer boundary)
    The shell between them is the meshed volume. A hole marker placed
    at the origin removes the cavity interior.

References:
    Si, H. (2015), "TetGen, a Delaunay-based quality tetrahedral mesh
    generator", ACM Trans. Math. Software 41(2), Article 11.

Usage example:
    python3 tools/mesh_generation/build_source_ball_mesh.py \\
        --cavity-radius 17.0 --elastic-radius 70.0 \\
        --edge-near 5.0 --edge-outer 12.0 \\
        --output cache/source_ball_meshes/salmon_default

This emits cache/source_ball_meshes/salmon_default.{node,ele}.
"""

import argparse
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path


def golden_icosahedron_vertices(radius):
    """Twelve vertices of a regular icosahedron on the sphere of given radius.

    Standard construction from three orthogonal golden rectangles (phi =
    (1 + sqrt(5)) / 2). Returns a list of (x, y, z) tuples.
    """
    phi = (1.0 + math.sqrt(5.0)) / 2.0
    raw = [
        (-1, phi, 0), (1, phi, 0), (-1, -phi, 0), (1, -phi, 0),
        (0, -1, phi), (0, 1, phi), (0, -1, -phi), (0, 1, -phi),
        (phi, 0, -1), (phi, 0, 1), (-phi, 0, -1), (-phi, 0, 1),
    ]
    norm = math.sqrt(1.0 + phi * phi)
    return [(x * radius / norm, y * radius / norm, z * radius / norm)
            for (x, y, z) in raw]


ICOSAHEDRON_FACES = [
    # Twenty triangular faces (CCW outward normal). Standard tabulation.
    (0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11),
    (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6), (7, 1, 8),
    (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9),
    (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1),
]


def subdivide_triangle_grid(vertices, faces, radius, levels):
    """Subdivide each triangular face into 4^levels sub-triangles.

    Each midpoint is normalized to the sphere of given radius. Returns
    (verts_out, faces_out) with deduplicated vertices.
    """
    if levels == 0:
        return list(vertices), list(faces)

    verts_out = list(vertices)
    midpoint_cache = {}

    def midpoint_index(a, b):
        key = (min(a, b), max(a, b))
        if key in midpoint_cache:
            return midpoint_cache[key]
        ax, ay, az = verts_out[a]
        bx, by, bz = verts_out[b]
        mx, my, mz = 0.5 * (ax + bx), 0.5 * (ay + by), 0.5 * (az + bz)
        mlen = math.sqrt(mx * mx + my * my + mz * mz)
        scale = radius / mlen
        verts_out.append((mx * scale, my * scale, mz * scale))
        idx = len(verts_out) - 1
        midpoint_cache[key] = idx
        return idx

    faces_out = []
    for (a, b, c) in faces:
        ab = midpoint_index(a, b)
        bc = midpoint_index(b, c)
        ca = midpoint_index(c, a)
        faces_out.append((a, ab, ca))
        faces_out.append((b, bc, ab))
        faces_out.append((c, ca, bc))
        faces_out.append((ab, bc, ca))

    if levels > 1:
        return subdivide_triangle_grid(verts_out, faces_out, radius, levels - 1)
    return verts_out, faces_out


def icosphere(radius, edge_target):
    """Icosphere whose triangle edge length is approximately edge_target.

    Raw icosahedron edge length is radius * 1.0515. Each subdivision
    halves the edge, so levels = ceil(log2(L_raw / edge_target)).
    """
    raw_edge = radius * 1.0515
    if edge_target <= 0.0:
        levels = 0
    else:
        ratio = raw_edge / edge_target
        if ratio <= 1.0:
            levels = 0
        else:
            levels = max(0, int(math.ceil(math.log2(ratio))))
    base = golden_icosahedron_vertices(radius)
    return subdivide_triangle_grid(base, ICOSAHEDRON_FACES, radius, levels)


def write_poly(path, inner_verts, inner_faces, outer_verts, outer_faces,
               region_target_volume):
    """Write a TetGen .poly PSLG file describing the spherical shell.

    File layout (TetGen 1.6 manual sec 2.1):
        # Node section
        N_nodes 3 0 1
        i x y z marker_i      (marker = 1 for outer surface, 2 for inner)
        ...
        # Facet section
        N_facets 1
        1
        3 i j k
        ...
        # Hole section
        1
        1 0.0 0.0 0.0
        # Region section
        1
        1 region_x region_y region_z 0 region_target_volume
    """
    n_inner = len(inner_verts)
    n_outer = len(outer_verts)
    n_total = n_inner + n_outer
    with open(path, "w") as f:
        f.write(f"# pass-13a source-ball PSLG\n")
        f.write(f"{n_total} 3 0 1\n")
        for i, v in enumerate(inner_verts):
            f.write(f"{i + 1} {v[0]:.10g} {v[1]:.10g} {v[2]:.10g} 2\n")
        for i, v in enumerate(outer_verts):
            f.write(f"{i + 1 + n_inner} {v[0]:.10g} "
                    f"{v[1]:.10g} {v[2]:.10g} 1\n")

        n_facets = len(inner_faces) + len(outer_faces)
        f.write(f"{n_facets} 1\n")
        for face in inner_faces:
            a, b, c = face
            f.write("1 0 2\n")
            f.write(f"3 {a + 1} {b + 1} {c + 1}\n")
        for face in outer_faces:
            a, b, c = face
            f.write("1 0 1\n")
            f.write(f"3 {a + 1 + n_inner} {b + 1 + n_inner} "
                    f"{c + 1 + n_inner}\n")

        f.write("1\n")
        f.write("1 0.0 0.0 0.0\n")

        avg_radius = 0.5 * (
            math.sqrt(sum(c * c for c in inner_verts[0]))
            + math.sqrt(sum(c * c for c in outer_verts[0])))
        rx = avg_radius
        f.write("1\n")
        f.write(f"1 {rx:.10g} 0.0 0.0 1 {region_target_volume:.10g}\n")


def run_tetgen(poly_path, quality):
    """Invoke `tetgen -pq{quality}aA` on the .poly file."""
    tetgen_bin = shutil.which("tetgen")
    if tetgen_bin is None:
        raise FileNotFoundError(
            "TetGen not found on PATH. Install from "
            "https://www.wias-berlin.de/software/tetgen (BSD licensed) "
            "or from your distribution's package manager.")
    cmd = [tetgen_bin, "-pq", str(quality), "-A", "-Y", str(poly_path)]
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        sys.stderr.write(result.stdout)
        sys.stderr.write(result.stderr)
        raise RuntimeError(f"tetgen failed (exit {result.returncode}): "
                           f"{' '.join(cmd)}")
    return result.stdout


def main():
    parser = argparse.ArgumentParser(
        description="Generate a TetGen .node/.ele pair for the axis-1b "
                    "3D source ball.")
    parser.add_argument("--cavity-radius", type=float, required=True,
                        help="Inner sphere radius [m].")
    parser.add_argument("--elastic-radius", type=float, required=True,
                        help="Outer sphere radius [m].")
    parser.add_argument("--edge-near", type=float, required=True,
                        help="Target edge length on the inner surface [m].")
    parser.add_argument("--edge-outer", type=float, required=True,
                        help="Target edge length on the outer surface [m].")
    parser.add_argument("--quality", type=float, default=1.4,
                        help="TetGen radius-edge ratio (default 1.4).")
    parser.add_argument("--output", type=str, required=True,
                        help="Output basename. Emits <output>.node, "
                             "<output>.ele.")
    parser.add_argument("--keep-poly", action="store_true",
                        help="Preserve the intermediate .poly file.")
    args = parser.parse_args()

    if args.cavity_radius <= 0.0 or args.elastic_radius <= args.cavity_radius:
        raise SystemExit("require 0 < cavity_radius < elastic_radius")
    if args.edge_near <= 0.0 or args.edge_outer <= 0.0:
        raise SystemExit("edge sizes must be positive")

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    inner_verts, inner_faces = icosphere(args.cavity_radius, args.edge_near)
    outer_verts, outer_faces = icosphere(args.elastic_radius, args.edge_outer)

    poly_path = out_path.with_suffix(".poly")
    region_target_volume = (args.edge_outer ** 3) / 6.0
    write_poly(poly_path, inner_verts, inner_faces,
               outer_verts, outer_faces, region_target_volume)

    print(f"[build_source_ball_mesh] icosphere counts: "
          f"inner={len(inner_verts)} verts / {len(inner_faces)} tris, "
          f"outer={len(outer_verts)} verts / {len(outer_faces)} tris",
          file=sys.stderr)
    print(f"[build_source_ball_mesh] region volume target: "
          f"{region_target_volume:.6g} m^3", file=sys.stderr)

    log = run_tetgen(poly_path, args.quality)
    sys.stderr.write(log)

    base_no_ext = poly_path.with_suffix("")
    src_node = base_no_ext.with_suffix(".1.node")
    src_ele = base_no_ext.with_suffix(".1.ele")
    dst_node = out_path.with_suffix(".node")
    dst_ele = out_path.with_suffix(".ele")
    if not src_node.exists() or not src_ele.exists():
        raise RuntimeError(
            f"TetGen did not emit expected output files: "
            f"{src_node}, {src_ele}")
    shutil.copy(src_node, dst_node)
    shutil.copy(src_ele, dst_ele)

    extras = list(base_no_ext.parent.glob(f"{base_no_ext.name}.1.*"))
    for f in extras:
        try:
            os.remove(f)
        except OSError:
            pass
    if not args.keep_poly:
        try:
            os.remove(poly_path)
        except OSError:
            pass

    print(f"[build_source_ball_mesh] wrote {dst_node}", file=sys.stderr)
    print(f"[build_source_ball_mesh] wrote {dst_ele}", file=sys.stderr)


if __name__ == "__main__":
    main()
