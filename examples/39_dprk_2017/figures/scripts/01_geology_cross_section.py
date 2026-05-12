#!/usr/bin/env python3
"""01_geology_cross_section.py -- DPRK 2017 (Punggye-ri, Sixth Test)

Read the layered velocity model from config.config and render a labelled
cross-section with the source location marked. Visualization only.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import matplotlib.patches as patches
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[4]
EXAMPLE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tools" / "figures"))
from figure_style import FIG_SINGLE, PALETTE, apply_style, savefig  # noqa: E402

LAYER_LABELS = {
    "LAYER_1": "Volcanic tuff / overburden\nvp=2400 m/s, rho=2200 kg/m3",
    "LAYER_2": "Competent granite (host rock)\nvp=5800 m/s, rho=2700 kg/m3",
    "LAYER_3": "Pre-Cambrian metamorphic basement\nvp=6300 m/s, rho=2800 kg/m3",
}


def parse_config(config_path: Path):
    """Return ([layers], (loc_x, loc_y, loc_z), (Lx, Ly, Lz))."""
    text = config_path.read_text()
    layers = []
    cur = None
    current_section = None
    for line in text.splitlines():
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        m = re.match(r"\[([A-Z_0-9]+)\]", line)
        if m:
            current_section = m.group(1)
            if current_section.startswith("LAYER_"):
                cur = {"name": current_section}
                layers.append(cur)
            else:
                cur = {"name": current_section}
            continue
        if "=" in line and cur is not None:
            k, v = (s.strip() for s in line.split("=", 1))
            cur[k] = v

    by_name = {}
    for l in layers:
        if l["name"].startswith("LAYER_"):
            by_name[l["name"]] = l

    layers_clean = []
    for n in sorted(by_name.keys()):
        l = by_name[n]
        layers_clean.append({
            "name": n,
            "z_top": float(l["z_top"]),
            "z_bottom": float(l["z_bottom"]),
            "rho": float(l.get("rho", 0)),
        })

    Lx = Ly = Lz = 0.0
    grid_block = re.search(r"\[GRID\](.*?)\[", text, re.DOTALL)
    if grid_block:
        for k in ("Lx", "Ly", "Lz"):
            m = re.search(rf"^{k}\s*=\s*([0-9.eE+-]+)",
                          grid_block.group(1), re.MULTILINE)
            if m:
                if k == "Lx":
                    Lx = float(m.group(1))
                elif k == "Ly":
                    Ly = float(m.group(1))
                else:
                    Lz = float(m.group(1))

    src = {"location_x": 0.0, "location_z": 0.0}
    src_block = re.search(r"\[EXPLOSION_SOURCE\](.*?)(\[|$)", text, re.DOTALL)
    if src_block:
        for k in ("location_x", "location_z", "depth_of_burial"):
            m = re.search(rf"^{k}\s*=\s*([0-9.eE+-]+)",
                          src_block.group(1), re.MULTILINE)
            if m:
                src[k] = float(m.group(1))

    return layers_clean, src, (Lx, Ly, Lz)


def main():
    apply_style()
    config_path = EXAMPLE / "config.config"
    layers, src, (Lx, Ly, Lz) = parse_config(config_path)

    fig, ax = plt.subplots(figsize=FIG_SINGLE)

    layer_colors = [PALETTE["low_tier"], PALETTE["med_tier"],
                    PALETTE["high_tier"], PALETTE["highest_tier"]]
    for i, l in enumerate(layers):
        z0 = l["z_bottom"]
        z1 = l["z_top"]
        ax.add_patch(patches.Rectangle(
            (0, z0), Lx, z1 - z0,
            facecolor=layer_colors[i % 4],
            alpha=0.45,
            edgecolor=PALETTE["neutral"],
            linewidth=0.8,
        ))
        label = LAYER_LABELS.get(l["name"],
                                 f"{l['name']}: rho={l['rho']:.0f} kg/m3")
        ax.text(Lx / 2, (z0 + z1) / 2, label,
                ha="center", va="center", fontsize=8)

    ax.plot(src["location_x"], src["location_z"],
            "*", color=PALETTE["annotation"],
            markersize=18, zorder=5, label="Source (~250 kt)")
    ax.annotate(
        f"  600 m depth\n  250 kt granite",
        xy=(src["location_x"], src["location_z"]),
        xytext=(src["location_x"] + Lx * 0.08, src["location_z"] - Lz * 0.08),
        fontsize=8, color=PALETTE["annotation"],
        arrowprops=dict(arrowstyle="->", color=PALETTE["annotation"],
                        lw=0.8),
    )

    # Mark surface seismometer positions
    seis = {
        "DPRK_EPI": 4000.0,
        "DPRK_NEAR_E": 5500.0,
        "DPRK_FAR_E": 6500.0,
    }
    for name, x in seis.items():
        ax.plot(x, Lz, "v", color=PALETTE["neutral"],
                markersize=7, zorder=4)
        ax.text(x, Lz + Lz * 0.01, name, ha="center", va="bottom",
                fontsize=7, rotation=45)

    ax.set_xlim(0, Lx)
    ax.set_ylim(-Lz * 0.02, Lz * 1.12)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("z (m, 0 = domain bottom)")
    ax.set_title(
        "DPRK 2017 (Punggye-ri, 6th test): layered geology cross-section\n"
        "Mt. Mantap, Kilju County -- 250 kt @ 600 m depth in granite"
    )
    ax.legend(loc="lower right")

    out = Path(__file__).resolve().parents[1] / "01_geology_cross_section.png"
    savefig(fig, out)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
