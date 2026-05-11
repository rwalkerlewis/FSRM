#!/usr/bin/env python3
"""03_moment_tensor_history.py -- Cannikin 1971

M(t) and Mdot(t) for the six independent moment-tensor components.
Reads near_field_history.csv written by the radial Lagrangian solver.
Visualization only.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[4]
EXAMPLE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tools" / "figures"))
from figure_style import FIG_WIDE, PALETTE, apply_style, savefig  # noqa: E402


def main():
    apply_style()
    csv_path = EXAMPLE / "output" / "near_field_history.csv"
    if not csv_path.exists():
        print(f"WARN: {csv_path} not found; run the simulation first.",
              file=sys.stderr)
        return 1

    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    if not rows:
        print(f"WARN: {csv_path} empty", file=sys.stderr)
        return 1

    t = np.array([float(r["time_s"]) for r in rows])
    components = ["xx", "yy", "zz", "xy", "xz", "yz"]
    M = {c: np.array([float(r[f"M_{c}_Nm"]) for r in rows])
         for c in components}
    Mdot = {c: np.array([float(r[f"Mdot_{c}_Nms"]) for r in rows])
            for c in components}

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=FIG_WIDE)
    color_keys = ["synthetic", "observed", "analytic",
                  "low_tier", "med_tier", "high_tier"]
    for c, key in zip(components, color_keys):
        ax1.plot(t, M[c], color=PALETTE[key], label=f"M_{c}")
        ax2.plot(t, Mdot[c], color=PALETTE[key], label=f"Mdot_{c}")
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel("M (N\u00b7m)")
    ax1.set_title("Moment tensor M(t)")
    ax1.legend(fontsize=7, ncol=2)
    ax2.set_xlabel("Time (s)")
    ax2.set_ylabel(r"$\dot{M}$ (N\u00b7m/s)")
    ax2.set_title("Moment-rate $\\dot{M}$(t)")
    ax2.legend(fontsize=7, ncol=2)

    fig.suptitle("Cannikin 1971: 6-component moment tensor history",
                 y=1.02)
    out = Path(__file__).resolve().parents[1] / \
        "03_moment_tensor_history.png"
    savefig(fig, out)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
