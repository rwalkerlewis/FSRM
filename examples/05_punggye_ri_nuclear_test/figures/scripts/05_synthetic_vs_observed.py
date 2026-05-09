#!/usr/bin/env python3
"""05_synthetic_vs_observed.py -- Punggye-ri 2017

Synthetic and observed traces overlaid for three named stations, with
cross-correlation values annotated. Reads SAC outputs from output/ for
the synthetic and from tools/waveform_vv/cache/punggye_ri/ for the
observed.

Visualization only: cross-correlation is computed for display, but the
underlying gate is enforced by the C++ waveform_vv test.
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
from figure_style import PALETTE, apply_style, savefig  # noqa: E402

CACHE_DIR = REPO / "tools" / "waveform_vv" / "cache" / "punggye_ri"

# Station identifiers in the cache; pick three as the comparison set.
DEFAULT_STATIONS = ["PUN_1", "PUN_2", "PUN_3"]


def normalized_xcorr(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    n = min(len(a), len(b))
    a = a[:n]
    b = b[:n]
    a = a - a.mean()
    b = b - b.mean()
    den = np.sqrt((a * a).sum() * (b * b).sum())
    if den <= 0:
        return float("nan")
    return float((a * b).sum() / den)


def main():
    apply_style()
    if not CACHE_DIR.exists() or not list(CACHE_DIR.glob("*.sac")):
        print(f"WARN: {CACHE_DIR} empty; populate with "
              f"tools/waveform_vv/refresh.py --event Salmon1964",
              file=sys.stderr)
        return 1

    out_sac_dir = EXAMPLE / "output" / "seismograms"
    if not out_sac_dir.exists():
        print(f"WARN: {out_sac_dir} not found", file=sys.stderr)
        return 1

    fig, axes = plt.subplots(3, 1, figsize=(8, 6), sharex=True)

    for i, name in enumerate(DEFAULT_STATIONS):
        ax = axes[i]
        # Try to find a synthetic SAC matching this name
        cands = list(out_sac_dir.glob(f"*{name}*.sac"))
        cand_obs = list(CACHE_DIR.glob(f"*{name}*.sac"))
        if not cands or not cand_obs:
            ax.text(0.5, 0.5, f"missing data for {name}",
                    ha="center", va="center", transform=ax.transAxes)
            continue
        st_syn = read(str(cands[0]))
        st_obs = read(str(cand_obs[0]))
        # Bandpass match the source-physics-dominated band
        for st in (st_syn, st_obs):
            st[0].filter("bandpass", freqmin=0.5, freqmax=5.0,
                         corners=4, zerophase=True)
        t_syn = st_syn[0].times()
        t_obs = st_obs[0].times()
        ax.plot(t_syn, st_syn[0].data,
                color=PALETTE["synthetic"], label="synthetic")
        ax.plot(t_obs, st_obs[0].data,
                color=PALETTE["observed"], label="observed",
                alpha=0.8)
        cc = normalized_xcorr(st_syn[0].data, st_obs[0].data)
        ax.set_title(f"{name}: cross-correlation = {cc:.3f}",
                     fontsize=10)
        if i == 0:
            ax.legend(loc="upper right", fontsize=8)
        ax.set_ylabel("velocity (m/s)")

    axes[-1].set_xlabel("t (s)")
    fig.suptitle("Punggye-ri 2017: synthetic vs observed (0.5-5 Hz band)",
                 y=1.00)
    out = Path(__file__).resolve().parents[1] / \
        "05_synthetic_vs_observed.png"
    savefig(fig, out)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
