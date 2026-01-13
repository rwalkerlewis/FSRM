#!/usr/bin/env python3
"""
Generate a Divider-style 3D Gmsh mesh *from scratch* using the Gmsh Python API.

This does NOT consume the existing .geo file. It recreates the Divider example's
main ingredients:
  - stratigraphic layers (boxes)
  - intrusive dome (sphere)
  - channel/lens body (box)
  - chimney/breccia pipe (box)
  - divider + secondary "fault core" corridors (boxes)
  - nested damage halos around an explosion point (classified by radius)
  - physical groups with names matching config/us_divider_1992_complex_geology_gmsh.config
  - boundary physical surfaces: surface, bottom, upstream, downstream, north, south

Output defaults to a *binary* MSH v4 file.

Usage:
  python3 scripts/generate_divider_binary_mesh_from_scratch.py
  python3 scripts/generate_divider_binary_mesh_from_scratch.py --out meshes/us_divider/divider_complex_geology_3d_bin_from_scratch.msh
  python3 scripts/generate_divider_binary_mesh_from_scratch.py --threads 8

Notes:
  - This is a "Divider-style" mesh generator, not a byte-for-byte identical
    reproduction of meshes/us_divider/divider_complex_geology_3d.msh.
  - Physical group NAMES are what FSRM uses for mappings (via $PhysicalNames).
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

import gmsh


def _unique_tags(dimtags: Iterable[Tuple[int, int]], dim: int) -> List[int]:
    out: List[int] = []
    seen: Set[int] = set()
    for d, t in dimtags:
        if d != dim:
            continue
        if t in seen:
            continue
        seen.add(t)
        out.append(t)
    return out


def _bbox_entities(dim: int, xmin: float, ymin: float, zmin: float, xmax: float, ymax: float, zmax: float) -> List[int]:
    dimtags = gmsh.model.getEntitiesInBoundingBox(xmin, ymin, zmin, xmax, ymax, zmax, dim)
    return _unique_tags(dimtags, dim)


def _center_of_mass(dim: int, tag: int) -> Tuple[float, float, float]:
    return gmsh.model.occ.getCenterOfMass(dim, tag)


def _dist(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> float:
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def _add_phys(dim: int, name: str, tags: List[int], phys_id: int) -> None:
    if not tags:
        # Keep generation robust even if a group ends up empty.
        return
    gmsh.model.addPhysicalGroup(dim, tags, phys_id)
    gmsh.model.setPhysicalName(dim, phys_id, name)


def main() -> int:
    repo_root = Path(__file__).resolve().parents[1]
    default_out = repo_root / "meshes" / "us_divider" / "divider_complex_geology_3d_bin_from_scratch.msh"

    p = argparse.ArgumentParser(description="Generate Divider-style binary Gmsh mesh from scratch (Python API).")
    p.add_argument("--out", type=Path, default=default_out, help=f"Output .msh path (default: {default_out})")
    p.add_argument("--threads", type=int, default=None, help="Number of meshing threads.")
    p.add_argument("--bin", type=int, choices=[0, 1], default=1, help="Binary (1) or ASCII (0) output (default: 1)")

    # Mesh sizing (defaults are deliberately coarser than the .geo example to keep runtime reasonable)
    p.add_argument("--lc-far", type=float, default=700.0)
    p.add_argument("--lc-mid", type=float, default=350.0)
    p.add_argument("--lc-near", type=float, default=200.0)
    p.add_argument("--lc-ultra", type=float, default=140.0)

    args = p.parse_args()

    out = args.out if args.out.is_absolute() else (repo_root / args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    if args.threads is not None:
        gmsh.option.setNumber("General.NumThreads", float(args.threads))
        gmsh.option.setNumber("Mesh.MaxNumThreads1D", float(args.threads))
        gmsh.option.setNumber("Mesh.MaxNumThreads2D", float(args.threads))
        gmsh.option.setNumber("Mesh.MaxNumThreads3D", float(args.threads))

    gmsh.model.add("divider_from_scratch")

    # -------------------------------------------------------------------------
    # Geometry parameters (mirrors meshes/us_divider/divider_complex_geology.geo)
    # -------------------------------------------------------------------------
    Lx = 12000.0
    Ly = 8000.0
    zTop = 0.0
    zSoil = -50.0
    zA = -200.0
    zT1 = -550.0
    zT2 = -900.0
    zWT = -1500.0
    zV = -2050.0
    zC = -2700.0
    zM = -3300.0
    zB = -4200.0

    xDiv = 6000.0
    core_t = 80.0
    xDivL = xDiv - core_t / 2.0
    xDivR = xDiv + core_t / 2.0

    xF2_top = 8900.0
    xF2_bot = 9650.0
    core2_t = 120.0
    xF2_min = min(xF2_top, xF2_bot)
    xF2_max = max(xF2_top, xF2_bot)

    dome_x, dome_y, dome_z, dome_r = 8500.0, 4200.0, -1800.0, 700.0

    lens_x, lens_y, lens_z = 7600.0, 2400.0, -700.0
    lens_rx, lens_ry, lens_rz = 1200.0, 450.0, 220.0

    pipe_x, pipe_y = 5650.0, 3950.0
    chimney_r = 180.0
    pipe_z0 = zWT + 5.0
    pipe_z1 = zSoil - 5.0

    GZx, GZy, GZz = 5600.0, 4000.0, -800.0
    R_crushed, R_fractured, R_damaged = 260.0, 650.0, 1300.0

    # -------------------------------------------------------------------------
    # Create base stratigraphic boxes (these overlap; we'll fragment everything)
    # -------------------------------------------------------------------------
    occ = gmsh.model.occ

    # Box(x, y, z, dx, dy, dz)
    v_soil = occ.addBox(0, 0, zSoil, Lx, Ly, zTop - zSoil)
    v_alluv = occ.addBox(0, 0, zA, Lx, Ly, zSoil - zA)
    v_tuff1 = occ.addBox(0, 0, zT1, Lx, Ly, zA - zT1)
    v_tuff2 = occ.addBox(0, 0, zT2, Lx, Ly, zT1 - zT2)
    v_weld = occ.addBox(0, 0, zWT, Lx, Ly, zT2 - zWT)
    v_vit = occ.addBox(0, 0, zV, Lx, Ly, zWT - zV)
    v_carb = occ.addBox(0, 0, zC, Lx, Ly, zV - zC)
    v_meta = occ.addBox(0, 0, zM, Lx, Ly, zC - zM)
    v_base = occ.addBox(0, 0, zB, Lx, Ly, zM - zB)

    # Heterogeneities / special volumes
    v_dome = occ.addSphere(dome_x, dome_y, dome_z, dome_r)
    v_lens = occ.addBox(lens_x - lens_rx, lens_y - lens_ry, lens_z - lens_rz, 2 * lens_rx, 2 * lens_ry, 2 * lens_rz)
    v_pipe = occ.addBox(pipe_x - chimney_r, pipe_y - chimney_r, pipe_z0, 2 * chimney_r, 2 * chimney_r, pipe_z1 - pipe_z0)

    v_divcore = occ.addBox(xDiv - core_t / 2.0, 0, zB, core_t, Ly, zTop - zB)
    v_f2core = occ.addBox(xF2_min - core2_t / 2.0, 0, zB, (xF2_max - xF2_min) + core2_t, Ly, zTop - zB)

    v_damage = occ.addSphere(GZx, GZy, GZz, R_damaged)

    occ.synchronize()

    # Fragment all volumes to make a conforming, disjoint partition
    base = [(3, t) for t in [v_soil, v_alluv, v_tuff1, v_tuff2, v_weld, v_vit, v_carb, v_meta, v_base]]
    tools = [(3, t) for t in [v_dome, v_lens, v_pipe, v_divcore, v_f2core, v_damage]]
    occ.fragment(base, tools, removeObject=True, removeTool=True)
    occ.synchronize()

    # -------------------------------------------------------------------------
    # Physical groups: classify resulting volumes by center-of-mass heuristics
    # -------------------------------------------------------------------------
    vols = _unique_tags(gmsh.model.getEntities(3), 3)

    # Helper flags
    special: Set[int] = set()

    # Precompute centers
    centers: Dict[int, Tuple[float, float, float]] = {t: _center_of_mass(3, t) for t in vols}

    # Identify special volumes first (dome/lens/pipe/cores/damage shells)
    intrusive_dome: List[int] = []
    tuff_channel_lens: List[int] = []
    chimney_pipe: List[int] = []
    divider_fault_core: List[int] = []
    secondary_fault_core: List[int] = []
    damage_crushed: List[int] = []
    damage_fractured: List[int] = []
    damage_damaged: List[int] = []

    dome_center = (dome_x, dome_y, dome_z)
    src_center = (GZx, GZy, GZz)

    for t in vols:
        c = centers[t]

        # Corridor cores (simple x tests)
        if abs(c[0] - xDiv) <= (core_t * 0.55):
            divider_fault_core.append(t)
            special.add(t)
            continue
        if (xF2_min - core2_t * 0.55) <= c[0] <= (xF2_max + core2_t * 0.55):
            # Narrow it slightly by also requiring it's not too close to xDiv corridor.
            if abs(c[0] - xDiv) > core_t:
                secondary_fault_core.append(t)
                special.add(t)
                continue

        # Chimney
        if (pipe_x - chimney_r) <= c[0] <= (pipe_x + chimney_r) and (pipe_y - chimney_r) <= c[1] <= (pipe_y + chimney_r) and pipe_z0 <= c[2] <= pipe_z1:
            chimney_pipe.append(t)
            special.add(t)
            continue

        # Lens
        if (lens_x - lens_rx) <= c[0] <= (lens_x + lens_rx) and (lens_y - lens_ry) <= c[1] <= (lens_y + lens_ry) and (lens_z - lens_rz) <= c[2] <= (lens_z + lens_rz):
            tuff_channel_lens.append(t)
            special.add(t)
            continue

        # Dome
        if _dist(c, dome_center) <= dome_r * 1.05:
            intrusive_dome.append(t)
            special.add(t)
            continue

        # Damage halos by radius band
        r = _dist(c, src_center)
        if r <= R_crushed:
            damage_crushed.append(t)
            special.add(t)
            continue
        if R_crushed < r <= R_fractured:
            damage_fractured.append(t)
            special.add(t)
            continue
        if R_fractured < r <= R_damaged:
            damage_damaged.append(t)
            special.add(t)
            continue

    # Background volumes: classify by stratigraphic z band and upstream/downstream, excluding specials
    soil_upstream: List[int] = []
    soil_downstream: List[int] = []
    alluvium_upstream: List[int] = []
    alluvium_downstream: List[int] = []
    tuff1_upstream: List[int] = []
    tuff1_downstream: List[int] = []
    tuff2_upstream: List[int] = []
    tuff2_downstream: List[int] = []
    welded_tuff_upstream: List[int] = []
    welded_tuff_downstream: List[int] = []
    vitric_tuff_upstream: List[int] = []
    vitric_tuff_downstream: List[int] = []
    carbonate_upstream: List[int] = []
    carbonate_downstream: List[int] = []
    metasediment_upstream: List[int] = []
    metasediment_downstream: List[int] = []
    basement_upstream: List[int] = []
    basement_downstream: List[int] = []

    def side(x: float) -> Optional[str]:
        if x < xDivL:
            return "up"
        if x > xDivR:
            return "down"
        return None  # inside divider corridor; should have been classified as special

    for t in vols:
        if t in special:
            continue
        x, y, z = centers[t]
        s = side(x)
        if s is None:
            continue

        if zSoil <= z <= zTop:
            (soil_upstream if s == "up" else soil_downstream).append(t)
        elif zA <= z < zSoil:
            (alluvium_upstream if s == "up" else alluvium_downstream).append(t)
        elif zT1 <= z < zA:
            (tuff1_upstream if s == "up" else tuff1_downstream).append(t)
        elif zT2 <= z < zT1:
            (tuff2_upstream if s == "up" else tuff2_downstream).append(t)
        elif zWT <= z < zT2:
            (welded_tuff_upstream if s == "up" else welded_tuff_downstream).append(t)
        elif zV <= z < zWT:
            (vitric_tuff_upstream if s == "up" else vitric_tuff_downstream).append(t)
        elif zC <= z < zV:
            (carbonate_upstream if s == "up" else carbonate_downstream).append(t)
        elif zM <= z < zC:
            (metasediment_upstream if s == "up" else metasediment_downstream).append(t)
        elif zB <= z < zM:
            (basement_upstream if s == "up" else basement_downstream).append(t)

    # Add physical volumes with stable IDs (any IDs are fine; names are what matter for FSRM mappings)
    phys = 1
    for name, tags in [
        ("soil_upstream", soil_upstream),
        ("alluvium_upstream", alluvium_upstream),
        ("tuff1_upstream", tuff1_upstream),
        ("tuff2_upstream", tuff2_upstream),
        ("welded_tuff_upstream", welded_tuff_upstream),
        ("vitric_tuff_upstream", vitric_tuff_upstream),
        ("carbonate_upstream", carbonate_upstream),
        ("metasediment_upstream", metasediment_upstream),
        ("basement_upstream", basement_upstream),
        ("soil_downstream", soil_downstream),
        ("alluvium_downstream", alluvium_downstream),
        ("tuff1_downstream", tuff1_downstream),
        ("tuff2_downstream", tuff2_downstream),
        ("welded_tuff_downstream", welded_tuff_downstream),
        ("vitric_tuff_downstream", vitric_tuff_downstream),
        ("carbonate_downstream", carbonate_downstream),
        ("metasediment_downstream", metasediment_downstream),
        ("basement_downstream", basement_downstream),
        ("intrusive_dome", intrusive_dome),
        ("tuff_channel_lens", tuff_channel_lens),
        ("chimney_pipe", chimney_pipe),
        ("divider_fault_core", divider_fault_core),
        ("secondary_fault_core", secondary_fault_core),
        ("damage_crushed", damage_crushed),
        ("damage_fractured", damage_fractured),
        ("damage_damaged", damage_damaged),
    ]:
        _add_phys(3, name, tags, phys)
        phys += 1

    # -------------------------------------------------------------------------
    # Boundary physical surfaces
    # -------------------------------------------------------------------------
    eps = 1.0  # meters; keep robust against meshing/fragment tolerances
    surface = _bbox_entities(2, -eps, -eps, zTop - eps, Lx + eps, Ly + eps, zTop + eps)
    bottom = _bbox_entities(2, -eps, -eps, zB - eps, Lx + eps, Ly + eps, zB + eps)
    upstream = _bbox_entities(2, -eps, -eps, zB - eps, 0.0 + eps, Ly + eps, zTop + eps)
    downstream = _bbox_entities(2, Lx - eps, -eps, zB - eps, Lx + eps, Ly + eps, zTop + eps)
    south = _bbox_entities(2, -eps, -eps, zB - eps, Lx + eps, 0.0 + eps, zTop + eps)
    north = _bbox_entities(2, -eps, Ly - eps, zB - eps, Lx + eps, Ly + eps, zTop + eps)

    # Use higher IDs for surfaces to avoid overlap with volume IDs.
    sphys = 10_000
    for name, tags in [
        ("surface", surface),
        ("bottom", bottom),
        ("upstream", upstream),
        ("downstream", downstream),
        ("south", south),
        ("north", north),
    ]:
        _add_phys(2, name, tags, sphys)
        sphys += 1

    # -------------------------------------------------------------------------
    # Mesh size fields (Ball + Min), similar to divider_complex_geology.geo
    # -------------------------------------------------------------------------
    gmsh.model.mesh.field.add("Ball", 1)
    gmsh.model.mesh.field.setNumber(1, "XCenter", GZx)
    gmsh.model.mesh.field.setNumber(1, "YCenter", GZy)
    gmsh.model.mesh.field.setNumber(1, "ZCenter", GZz)
    gmsh.model.mesh.field.setNumber(1, "Radius", R_damaged + 500.0)
    gmsh.model.mesh.field.setNumber(1, "VIn", args.lc_near)
    gmsh.model.mesh.field.setNumber(1, "VOut", args.lc_far)

    gmsh.model.mesh.field.add("Ball", 2)
    gmsh.model.mesh.field.setNumber(2, "XCenter", GZx)
    gmsh.model.mesh.field.setNumber(2, "YCenter", GZy)
    gmsh.model.mesh.field.setNumber(2, "ZCenter", GZz)
    gmsh.model.mesh.field.setNumber(2, "Radius", R_fractured + 250.0)
    gmsh.model.mesh.field.setNumber(2, "VIn", args.lc_ultra)
    gmsh.model.mesh.field.setNumber(2, "VOut", args.lc_near)

    gmsh.model.mesh.field.add("Ball", 3)
    gmsh.model.mesh.field.setNumber(3, "XCenter", xDiv)
    gmsh.model.mesh.field.setNumber(3, "YCenter", Ly / 2.0)
    gmsh.model.mesh.field.setNumber(3, "ZCenter", (zTop + zB) / 2.0)
    gmsh.model.mesh.field.setNumber(3, "Radius", 1200.0)
    gmsh.model.mesh.field.setNumber(3, "VIn", args.lc_ultra)
    gmsh.model.mesh.field.setNumber(3, "VOut", args.lc_mid)

    gmsh.model.mesh.field.add("Ball", 4)
    gmsh.model.mesh.field.setNumber(4, "XCenter", (xF2_top + xF2_bot) / 2.0)
    gmsh.model.mesh.field.setNumber(4, "YCenter", Ly / 2.0)
    gmsh.model.mesh.field.setNumber(4, "ZCenter", (zTop + zB) / 2.0)
    gmsh.model.mesh.field.setNumber(4, "Radius", 1400.0)
    gmsh.model.mesh.field.setNumber(4, "VIn", args.lc_ultra)
    gmsh.model.mesh.field.setNumber(4, "VOut", args.lc_mid)

    gmsh.model.mesh.field.add("Ball", 7)
    gmsh.model.mesh.field.setNumber(7, "XCenter", dome_x)
    gmsh.model.mesh.field.setNumber(7, "YCenter", dome_y)
    gmsh.model.mesh.field.setNumber(7, "ZCenter", dome_z)
    gmsh.model.mesh.field.setNumber(7, "Radius", 1700.0)
    gmsh.model.mesh.field.setNumber(7, "VIn", args.lc_mid)
    gmsh.model.mesh.field.setNumber(7, "VOut", args.lc_far)

    gmsh.model.mesh.field.add("Ball", 8)
    gmsh.model.mesh.field.setNumber(8, "XCenter", lens_x)
    gmsh.model.mesh.field.setNumber(8, "YCenter", lens_y)
    gmsh.model.mesh.field.setNumber(8, "ZCenter", lens_z)
    gmsh.model.mesh.field.setNumber(8, "Radius", 1600.0)
    gmsh.model.mesh.field.setNumber(8, "VIn", args.lc_near)
    gmsh.model.mesh.field.setNumber(8, "VOut", args.lc_far)

    gmsh.model.mesh.field.add("Ball", 9)
    gmsh.model.mesh.field.setNumber(9, "XCenter", pipe_x)
    gmsh.model.mesh.field.setNumber(9, "YCenter", pipe_y)
    gmsh.model.mesh.field.setNumber(9, "ZCenter", (pipe_z0 + pipe_z1) / 2.0)
    gmsh.model.mesh.field.setNumber(9, "Radius", 750.0)
    gmsh.model.mesh.field.setNumber(9, "VIn", args.lc_ultra)
    gmsh.model.mesh.field.setNumber(9, "VOut", args.lc_mid)

    gmsh.model.mesh.field.add("Min", 20)
    gmsh.model.mesh.field.setNumbers(20, "FieldsList", [1, 2, 3, 4, 7, 8, 9])
    gmsh.model.mesh.field.setAsBackgroundMesh(20)

    # Global meshing options
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", args.lc_ultra)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", args.lc_far)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)  # Delaunay
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.option.setNumber("Mesh.Binary", float(args.bin))

    # Generate + write
    gmsh.model.mesh.generate(3)
    gmsh.write(str(out))
    gmsh.finalize()

    print(f"wrote: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

