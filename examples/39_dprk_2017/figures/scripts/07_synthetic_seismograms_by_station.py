#!/usr/bin/env python3
"""07_synthetic_seismograms_by_station.py -- DPRK 2017 (Punggye-ri, Sixth Test)

Produces one figure per seismic station, each with three stacked panels
(BHE, BHN, BHZ).  Output files are named
  07_seismogram_<STATION>.png
and written next to this script's parent figures/ directory.

SAC files are read from output/ (DPRK 2017 writes them directly there,
not into output/seismograms/).
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[4]
EXAMPLE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tools" / "figures"))
from figure_style import PALETTE, apply_style, plot_seismogram, savefig  # noqa: E402

STATION_ORDER = ["DPRK_EPI", "DPRK_NEAR_E", "DPRK_NEAR_N",
                 "DPRK_FAR_E", "DPRK_FAR_N"]

STATION_OFFSETS = {
    "DPRK_EPI":    0,
    "DPRK_NEAR_E": 1500,
    "DPRK_NEAR_N": 1500,
    "DPRK_FAR_E":  2500,
    "DPRK_FAR_N":  2500,
}

COMPONENT_ORDER = ["BHE", "BHN", "BHZ"]
COMPONENT_COLORS = {"BHE": "analytic", "BHN": "observed", "BHZ": "synthetic"}
COMPONENT_LABELS = {"BHE": "East (BHE)", "BHN": "North (BHN)", "BHZ": "Vertical (BHZ)"}

OUT_DIR = Path(__file__).resolve().parents[1]


def make_station_figure(station: str, files_by_comp: dict[str, Path]) -> int:
    """Render one figure for *station* with one panel per component."""
    fig, axes = plt.subplots(
        len(COMPONENT_ORDER), 1,
        figsize=(9, 7),
        sharex=True,
    )

    for ax, comp in zip(axes, COMPONENT_ORDER):
        f = files_by_comp.get(comp)
        if f is not None and f.exists():
            plot_seismogram(ax, f, label=comp,
                            color=COMPONENT_COLORS.get(comp, "synthetic"))
        else:
            ax.text(0.5, 0.5, f"{comp} not found",
                    ha="center", va="center", transform=ax.transAxes,
                    color="gray", fontsize=10)
            ax.set_xlim(0, 1)

        ax.set_ylabel(f"{COMPONENT_LABELS[comp]}\ndispl (m)", fontsize=9)
        ax.tick_params(labelsize=8)

    axes[-1].set_xlabel("Time (s)", fontsize=10)

    offset_m = STATION_OFFSETS.get(station, "?")
    fig.suptitle(
        f"DPRK 2017 (Punggye-ri, 6th test) -- station {station}\n"
        f"Epicentral distance {offset_m} m  |  250 kt @ 600 m depth",
        fontsize=11,
    )
    fig.tight_layout()

    tag = station.lower().replace(" ", "_")
    out = OUT_DIR / f"07_seismogram_{tag}.png"
    savefig(fig, out)
    print(f"wrote {out}")
    plt.close(fig)
    return 0


def main() -> int:
    apply_style()
    sac_dir = EXAMPLE / "output"

    sac_files = sorted(sac_dir.glob("*.sac"))
    if not sac_files:
        print(f"WARN: no .sac files in {sac_dir}; run the simulation first.",
              file=sys.stderr)
        return 1

    # Group: station -> comp -> Path
    grouped: dict[str, dict[str, Path]] = {s: {} for s in STATION_ORDER}
    unmatched: list[Path] = []
    for f in sac_files:
        matched = False
        for sta in STATION_ORDER:
            if f".{sta}." in f.name:
                for comp in COMPONENT_ORDER:
                    if f".{comp}." in f.name:
                        grouped[sta][comp] = f
                        matched = True
                        break
                if matched:
                    break
        if not matched:
            unmatched.append(f)

    if unmatched:
        print(f"INFO: {len(unmatched)} SAC file(s) did not match any "
              f"known station; skipping.", file=sys.stderr)

    rc = 0
    for sta in STATION_ORDER:
        comps = grouped.get(sta, {})
        if not comps:
            print(f"INFO: no SAC files found for {sta}; skipping.",
                  file=sys.stderr)
            continue
        rc |= make_station_figure(sta, comps)

    return rc


if __name__ == "__main__":
    sys.exit(main() or 0)
