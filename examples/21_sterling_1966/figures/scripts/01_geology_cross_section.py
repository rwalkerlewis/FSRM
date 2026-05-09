#!/usr/bin/env python3
"""01_geology_cross_section.py -- Sterling 1966

Read the layered velocity model from config.config and render a labelled
cross-section with the source location marked. Visualization only.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[4]
EXAMPLE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tools" / "figures"))
from figure_style import FIG_SINGLE, PALETTE, apply_style, savefig  # noqa: E402


def parse_config(config_path: Path):
    """Return ([layers], (loc_x, loc_y, loc_z), (Lx, Ly, Lz))."""
    text = config_path.read_text()
    layers = []
    cur = None
    for line in text.splitlines():
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        m = re.match(r"\[([A-Z_0-9]+)\]", line)
        if m:
            sec = m.group(1)
            if sec.startswith("LAYER_"):
                cur = {"name": sec}
                layers.append(cur)
            else:
                cur = {"name": sec}
            continue
        if "=" in line and cur is not None:
            k, v = (s.strip() for s in line.split("=", 1))
            cur[k] = v
    by_name = {l["name"]: l for l in [{"name": x["name"], **x} for x in
                                       (l for l in layers
                                        if l["name"].startswith("LAYER_"))]}
    layers_clean = []
    for n in sorted(by_name.keys()):
        l = by_name[n]
        layers_clean.append({
            "z_top": float(l["z_top"]),
            "z_bottom": float(l["z_bottom"]),
            "rho": float(l.get("rho", 0)),
        })
    grid = next((x for x in [{"name": s["name"], **s} for s in [l for l in
                                                                 [s for s in
                                                                  text.splitlines()] if False]]
                if False), None)
    Lx = Ly = Lz = 0.0
    grid_block = re.search(r"\[GRID\](.*?)\[", text, re.DOTALL)
    if grid_block:
        for k, var in (("Lx", "Lx"), ("Ly", "Ly"), ("Lz", "Lz")):
            m = re.search(rf"^{k}\s*=\s*([0-9.eE+-]+)", grid_block.group(1),
                          re.MULTILINE)
            if m:
                if var == "Lx":
                    Lx = float(m.group(1))
                elif var == "Ly":
                    Ly = float(m.group(1))
                else:
                    Lz = float(m.group(1))
    src = {"location_x": 0.0, "location_y": 0.0, "location_z": 0.0}
    src_block = re.search(r"\[EXPLOSION_SOURCE\](.*?)(\[|$)", text, re.DOTALL)
    if src_block:
        for k in ("location_x", "location_y", "location_z"):
            m = re.search(rf"^{k}\s*=\s*([0-9.eE+-]+)", src_block.group(1),
                          re.MULTILINE)
            if m:
                src[k] = float(m.group(1))
    return layers_clean, (src["location_x"], src["location_y"],
                          src["location_z"]), (Lx, Ly, Lz)


def main():
    apply_style()
    config_path = EXAMPLE / "config.config"
    layers, src, (Lx, Ly, Lz) = parse_config(config_path)

    fig, ax = plt.subplots(figsize=FIG_SINGLE)

    layer_colors = [PALETTE["low_tier"], PALETTE["med_tier"],
                    PALETTE["high_tier"], PALETTE["highest_tier"]]
    cy = 0
    for i, l in enumerate(layers):
        z0 = l["z_bottom"]
        z1 = l["z_top"]
        ax.add_patch(patches.Rectangle((0, z0), Lx, z1 - z0,
                                       facecolor=layer_colors[i % 4],
                                       alpha=0.5,
                                       edgecolor=PALETTE["neutral"]))
        ax.text(Lx / 2, (z0 + z1) / 2,
                f"{l['name'] if 'name' in l else f'L{i+1}'}: rho={l['rho']:.0f}",
                ha="center", va="center", fontsize=9)

    ax.plot(src[0], src[2], "*", color=PALETTE["annotation"],
            markersize=20, label="Source")
    ax.set_xlim(0, Lx)
    ax.set_ylim(0, Lz)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("z (m)")
    ax.set_title("Sterling 1966: layered velocity model + source")
    ax.legend(loc="lower right")

    out = Path(__file__).resolve().parents[1] / "01_geology_cross_section.png"
    savefig(fig, out)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
