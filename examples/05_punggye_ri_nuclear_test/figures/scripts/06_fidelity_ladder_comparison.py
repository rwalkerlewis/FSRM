#!/usr/bin/env python3
"""06_fidelity_ladder_comparison.py -- Punggye-ri 2017

Same event run at LOW, MED, HIGHEST tiers; key output (peak velocity at
the closest station) overlaid. Anchors the "fidelity costs and
trade-offs" presentation slide. Visualization only.

Expected layout: per-tier output directories from run_showcase.sh
landing at output/, output_marshak/, output_highest/. The script does
not run the simulation; it consumes the existing directories.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from obspy import read

REPO = Path(__file__).resolve().parents[4]
EXAMPLE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tools" / "figures"))
from figure_style import FIG_WIDE, PALETTE, apply_style, savefig  # noqa: E402

TIER_DIRS = [
    ("LOW (Z-R)", EXAMPLE / "output", "low_tier"),
    ("MED (Marshak grey)", EXAMPLE / "output_marshak", "med_tier"),
    ("HIGHEST (Marshak multigroup)",
     EXAMPLE / "output_highest", "highest_tier"),
]


def closest_sac(sac_dir):
    cand = list(sac_dir.glob("*.sac"))
    if not cand:
        return None
    # Pick the one with the smallest peak time as a heuristic for "closest"
    cand.sort()
    return cand[0]


def main():
    apply_style()
    fig, ax = plt.subplots(figsize=FIG_WIDE)

    have_any = False
    for label, base, color_key in TIER_DIRS:
        sac_dir = base / "seismograms"
        if not sac_dir.exists():
            print(f"WARN: {sac_dir} missing; skip {label}",
                  file=sys.stderr)
            continue
        sac_path = closest_sac(sac_dir)
        if sac_path is None:
            print(f"WARN: no SAC files in {sac_dir}", file=sys.stderr)
            continue
        st = read(str(sac_path))
        tr = st[0]
        ax.plot(tr.times(), tr.data, color=PALETTE[color_key],
                label=label)
        peak = np.max(np.abs(tr.data))
        ax.axhline(peak, color=PALETTE[color_key],
                   linestyle=":", alpha=0.5)
        have_any = True

    if not have_any:
        print("WARN: no tier outputs available; run run_showcase.sh first",
              file=sys.stderr)
        return 1

    ax.set_xlabel("t (s)")
    ax.set_ylabel("velocity (m/s)")
    ax.set_title("Punggye-ri 2017: peak velocity at closest station, "
                 "by fidelity tier")
    ax.legend(loc="upper right")

    out = Path(__file__).resolve().parents[1] / \
        "06_fidelity_ladder_comparison.png"
    savefig(fig, out)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
