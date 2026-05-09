#!/usr/bin/env python3
"""04_synthetic_seismograms_grid.py -- Sedan 1962

Time-domain seismograms at all configured stations as a small-multiples
grid. Reads SAC outputs. Visualization only.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[4]
EXAMPLE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tools" / "figures"))
from figure_style import PALETTE, apply_style, plot_seismogram, savefig  # noqa: E402


def main():
    apply_style()
    sac_dir = EXAMPLE / "output" / "seismograms"
    if not sac_dir.exists():
        print(f"WARN: {sac_dir} not found; run the simulation first.",
              file=sys.stderr)
        return 1

    sac_files = sorted(sac_dir.glob("*.sac"))
    if not sac_files:
        print(f"WARN: no .sac files in {sac_dir}", file=sys.stderr)
        return 1

    n = len(sac_files)
    cols = min(3, n)
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 2.5 * rows),
                             squeeze=False)

    for i, sac in enumerate(sac_files):
        ax = axes[i // cols][i % cols]
        plot_seismogram(ax, sac, label=sac.stem, color="synthetic")
        ax.set_title(sac.stem, fontsize=9)
        ax.set_xlabel("t (s)")
        ax.set_ylabel("velocity (m/s)")

    # Hide unused axes
    for j in range(n, rows * cols):
        axes[j // cols][j % cols].axis("off")

    fig.suptitle("Sedan 1962: synthetic seismograms (all stations)",
                 y=1.01)
    out = Path(__file__).resolve().parents[1] / \
        "04_synthetic_seismograms_grid.png"
    savefig(fig, out)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
