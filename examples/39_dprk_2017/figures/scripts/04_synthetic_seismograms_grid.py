#!/usr/bin/env python3
"""04_synthetic_seismograms_grid.py -- DPRK 2017 (Punggye-ri, Sixth Test)

Time-domain seismograms at all configured stations as a small-multiples
grid. Reads SAC outputs from output/*.sac (DPRK 2017 example writes SAC
files directly into output/, not into output/seismograms/).
Visualization only.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[4]
EXAMPLE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tools" / "figures"))
from figure_style import PALETTE, apply_style, plot_seismogram, savefig  # noqa: E402

# Station display order and colours
STATION_ORDER = ["DPRK_EPI", "DPRK_NEAR_E", "DPRK_NEAR_N",
                 "DPRK_FAR_E", "DPRK_FAR_N"]
COMPONENT_COLORS = {"BHZ": "synthetic", "BHE": "analytic", "BHN": "observed"}

# Station offsets from source (m) -- for annotation
STATION_OFFSETS = {
    "DPRK_EPI":    0,
    "DPRK_NEAR_E": 1500,
    "DPRK_NEAR_N": 1500,
    "DPRK_FAR_E":  2500,
    "DPRK_FAR_N":  2500,
}


def main():
    apply_style()
    sac_dir = EXAMPLE / "output"

    # SAC files land directly in output/: XX.DPRK_EPI.00.BHZ.sac etc.
    sac_files = sorted(sac_dir.glob("*.sac"))
    if not sac_files:
        print(f"WARN: no .sac files in {sac_dir}; run the simulation first.",
              file=sys.stderr)
        return 1

    # Group by station, sorted by STATION_ORDER
    by_station: dict[str, list[Path]] = {}
    for s in STATION_ORDER:
        matches = [f for f in sac_files if f"_{s}." in f.name]
        if matches:
            by_station[s] = sorted(matches)

    if not by_station:
        # Fallback: use all SAC files ungrouped
        by_station = {"all": sac_files}

    n_sta = len(by_station)
    n_comp = max(len(v) for v in by_station.values())
    cols = n_comp
    rows = n_sta

    fig, axes = plt.subplots(rows, cols,
                             figsize=(4.5 * cols, 2.2 * rows),
                             squeeze=False)

    for row_i, (sta, files) in enumerate(by_station.items()):
        for col_i in range(cols):
            ax = axes[row_i][col_i]
            if col_i < len(files):
                f = files[col_i]
                # Determine component from filename
                comp = "BHZ"
                for c in ("BHZ", "BHE", "BHN"):
                    if c in f.name:
                        comp = c
                        break
                color_key = COMPONENT_COLORS.get(comp, "synthetic")
                ret = plot_seismogram(ax, f, label=comp,
                                      color=color_key)
                offset = STATION_OFFSETS.get(sta, "?")
                ax.set_title(
                    f"{sta} | {comp} | Δ={offset} m",
                    fontsize=8,
                )
                ax.set_xlabel("t (s)")
                ax.set_ylabel("displ (m)")
            else:
                ax.axis("off")

    fig.suptitle(
        "DPRK 2017 (Punggye-ri, 6th test): synthetic seismograms\n"
        "250 kt @ 600 m depth in granite, 3-layer geology",
        y=1.01,
    )
    fig.tight_layout()

    out = Path(__file__).resolve().parents[1] / \
        "04_synthetic_seismograms_grid.png"
    savefig(fig, out)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
