#!/usr/bin/env python3
"""06_fidelity_ladder_comparison.py -- DPRK 2017 (Punggye-ri, Sixth Test)

Peak displacement at the epicentral station (DPRK_EPI.BHZ) overlaid
across fidelity tiers. The DPRK 2017 example ships only the LOW tier
(COUPLED_ANALYTIC). MED and HIGH tiers require re-running with
explosion_solve_mode = RADIAL_LAGRANGIAN and the Marshak radiation
options respectively (see config templates in config/templates/).

If only one tier directory is present, a bar chart of peak amplitude
by component is rendered instead.
Visualization only.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[4]
EXAMPLE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tools" / "figures"))
from figure_style import FIG_WIDE, PALETTE, apply_style, savefig  # noqa: E402

TIER_DIRS = [
    ("LOW (COUPLED_ANALYTIC)", EXAMPLE / "output", "low_tier"),
    ("MED (RADIAL_LAGRANGIAN)", EXAMPLE / "output_radlag", "med_tier"),
    ("HIGH (Marshak grey)", EXAMPLE / "output_marshak", "high_tier"),
]

TARGET_STATION_KEY = "DPRK_EPI"
TARGET_COMPONENT = "BHZ"


def read_sac_simple(path: Path):
    """Return (times, data) using obspy if available, else plain numpy."""
    try:
        from obspy import read as obspy_read  # noqa: PLC0415
        st = obspy_read(str(path))
        return st[0].times(), st[0].data.astype(float)
    except Exception:
        pass
    import struct  # noqa: PLC0415
    raw = path.read_bytes()
    for endian in ("<", ">"):
        delta = struct.unpack_from(f"{endian}f", raw, 0)[0]
        if 0.0001 < delta < 10.0:
            break
    npts = struct.unpack_from(f"{endian}i", raw, 316)[0]
    data = np.frombuffer(raw[632:632 + npts * 4], dtype=f"{endian}f4")
    return np.arange(npts) * delta, data.astype(float)


def main():
    apply_style()
    out = Path(__file__).resolve().parents[1] / \
        "06_fidelity_ladder_comparison.png"

    available = [(label, base, key)
                 for label, base, key in TIER_DIRS
                 if (base).exists() and list(base.glob("*.sac"))]

    if not available:
        print(
            "WARN: no tier output directories found. "
            "Run the simulation to populate output/.",
            file=sys.stderr,
        )
        return 1

    if len(available) == 1:
        # Single tier: render peak-amplitude bar chart by component
        label, base, color_key = available[0]
        sac_files = sorted(base.glob(f"*{TARGET_STATION_KEY}*.sac"))
        if not sac_files:
            sac_files = sorted(base.glob("*.sac"))

        components = []
        peaks = []
        for f in sac_files:
            _, data = read_sac_simple(f)
            components.append(f.stem.split(".")[-1])
            peaks.append(float(np.max(np.abs(data))))

        fig, (ax_trace, ax_bar) = plt.subplots(1, 2, figsize=FIG_WIDE)

        # Trace panel: plot the BHZ trace
        bhz = [f for f in sac_files if TARGET_COMPONENT in f.name]
        if bhz:
            t, d = read_sac_simple(bhz[0])
            ax_trace.plot(t, d, color=PALETTE[color_key], lw=1.2)
            ax_trace.set_xlabel("t (s)")
            ax_trace.set_ylabel("displ (m)")
            ax_trace.set_title(
                f"{TARGET_STATION_KEY}.{TARGET_COMPONENT} -- {label}"
            )

        # Bar chart panel
        x = np.arange(len(components))
        ax_bar.bar(x, peaks, color=PALETTE[color_key], alpha=0.8)
        ax_bar.set_xticks(x)
        ax_bar.set_xticklabels(components, rotation=30, ha="right")
        ax_bar.set_ylabel("peak |displ| (m)")
        ax_bar.set_title("Peak displacement by component")
        for xi, p in zip(x, peaks):
            ax_bar.text(xi, p * 1.02, f"{p:.2e}", ha="center",
                        fontsize=8, color=PALETTE["annotation"])

        fig.suptitle(
            f"DPRK 2017 (Punggye-ri): fidelity tier = {label}\n"
            "MED / HIGH tiers: re-run with RADIAL_LAGRANGIAN / Marshak configs",
            y=1.01,
        )
    else:
        # Multi-tier overlay
        fig, ax = plt.subplots(figsize=FIG_WIDE)
        for label, base, color_key in available:
            cands = sorted(base.glob(
                f"*{TARGET_STATION_KEY}*{TARGET_COMPONENT}*.sac"
            ))
            if not cands:
                continue
            t, d = read_sac_simple(cands[0])
            ax.plot(t, d, color=PALETTE[color_key], lw=1.3, label=label)
            peak = float(np.max(np.abs(d)))
            ax.axhline(peak, color=PALETTE[color_key],
                       ls=":", lw=0.7, alpha=0.6)

        ax.set_xlabel("t (s)")
        ax.set_ylabel("displ (m)")
        ax.set_title(
            f"DPRK 2017: peak displacement at {TARGET_STATION_KEY}.{TARGET_COMPONENT}"
        )
        ax.legend(loc="upper right")
        fig.suptitle(
            "DPRK 2017 (Punggye-ri, 6th test): fidelity tier comparison",
            y=1.01,
        )

    fig.tight_layout()
    savefig(fig, out)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
