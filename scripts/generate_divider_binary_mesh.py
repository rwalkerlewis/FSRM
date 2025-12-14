#!/usr/bin/env python3
"""
Generate a *binary* Gmsh .msh for the Divider example.

This script is intentionally simple: it calls the `gmsh` CLI on the provided
`.geo` file and writes a binary MSH v4 mesh.

Usage:
  python3 scripts/generate_divider_binary_mesh.py
  python3 scripts/generate_divider_binary_mesh.py --out meshes/us_divider/divider_complex_geology_3d_bin.msh
  python3 scripts/generate_divider_binary_mesh.py --threads 8

Notes:
  - The Divider geometry is complex; meshing can take a while.
  - If you want ASCII instead, use `--bin 0`.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


def main() -> int:
    repo_root = Path(__file__).resolve().parents[1]
    default_geo = repo_root / "meshes" / "us_divider" / "divider_complex_geology.geo"
    default_out = repo_root / "meshes" / "us_divider" / "divider_complex_geology_3d_bin.msh"

    parser = argparse.ArgumentParser(description="Generate Divider binary Gmsh mesh (.msh).")
    parser.add_argument(
        "--geo",
        type=Path,
        default=default_geo,
        help=f"Path to .geo file (default: {default_geo})",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=default_out,
        help=f"Output .msh path (default: {default_out})",
    )
    parser.add_argument(
        "--format",
        default="msh4",
        choices=["msh2", "msh4"],
        help="Gmsh output format (default: msh4)",
    )
    parser.add_argument(
        "--bin",
        type=int,
        default=1,
        choices=[0, 1],
        help="Write binary (1) or ASCII (0) (default: 1)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="Number of meshing threads (-nt). If omitted, gmsh default is used.",
    )

    args = parser.parse_args()

    gmsh = shutil.which("gmsh")
    if not gmsh:
        print("error: gmsh not found on PATH (install gmsh, or add it to PATH)", file=sys.stderr)
        return 2

    geo = args.geo if args.geo.is_absolute() else (repo_root / args.geo)
    out = args.out if args.out.is_absolute() else (repo_root / args.out)

    if not geo.exists():
        print(f"error: .geo file not found: {geo}", file=sys.stderr)
        return 2

    out.parent.mkdir(parents=True, exist_ok=True)

    cmd: list[str] = [
        gmsh,
        "-3",
        str(geo),
        "-format",
        args.format,
        "-bin",
        str(args.bin),
        "-o",
        str(out),
    ]
    if args.threads is not None:
        cmd.extend(["-nt", str(args.threads)])

    print("running:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    print("wrote:", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

