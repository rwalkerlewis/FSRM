#!/usr/bin/env python3
"""05_synthetic_vs_observed.py -- DPRK 2017 (Punggye-ri, Sixth Test)

Overlay synthetic displacement seismograms against publicly available
IRIS waveforms for the 2017-09-03 DPRK test. Reads SAC outputs from
output/ for the synthetic and from tools/waveform_vv/cache/dprk_2017/
for the observed.

The DPRK 2017 event is the most-instrumented modern underground nuclear
test, with hundreds of IRIS broadband stations recording mb 6.3.
Reference: Voytan et al. (2019) GRL 46, Tian et al. (2018) GJI 213.

Visualization only: the cross-correlation is computed for display; the
gate is enforced by Integration.DPRK2017Comparison in C++.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[4]
EXAMPLE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tools" / "figures"))
from figure_style import PALETTE, apply_style, savefig  # noqa: E402

CACHE_DIR = REPO / "tools" / "waveform_vv" / "cache" / "dprk_2017"

# Station match pattern: synthetic SAC file substring -> observed SAC substring
STATION_PAIRS = [
    ("DPRK_EPI", "DPRK_EPI",
     "Epicentral (0 m offset, surface)"),
    ("DPRK_NEAR_E", "DPRK_NEAR",
     "Near-field E (1500 m offset)"),
    ("DPRK_FAR_E", "DPRK_FAR",
     "Far-field E (2500 m offset)"),
]


def normalized_xcorr(a: np.ndarray, b: np.ndarray) -> float:
    n = min(len(a), len(b))
    a = np.asarray(a[:n], dtype=float)
    b = np.asarray(b[:n], dtype=float)
    a -= a.mean()
    b -= b.mean()
    den = np.sqrt((a * a).sum() * (b * b).sum())
    if den <= 0:
        return float("nan")
    return float((a * b).sum() / den)


def read_sac(path: Path):
    """Return (times, data) using obspy if available, else plain numpy."""
    try:
        from obspy import read as obspy_read  # noqa: PLC0415
        st = obspy_read(str(path))
        tr = st[0]
        tr.filter("bandpass", freqmin=0.5, freqmax=5.0,
                  corners=4, zerophase=True)
        return tr.times(), tr.data
    except Exception:
        pass
    # Minimal SAC binary reader (big- or little-endian, 70 header floats)
    import struct  # noqa: PLC0415
    raw = path.read_bytes()
    # Check endian via sampling rate sanity
    for endian in ("<", ">"):
        delta = struct.unpack_from(f"{endian}f", raw, 0)[0]
        if 0.0001 < delta < 10.0:
            break
    npts = struct.unpack_from(f"{endian}i", raw, 316)[0]
    data = np.frombuffer(raw[632:632 + npts * 4], dtype=f"{endian}f4")
    times = np.arange(npts) * delta
    return times, data.astype(float)


def main():
    apply_style()
    out_dir = EXAMPLE / "output"
    sac_files = list(out_dir.glob("*.sac"))

    cache_available = CACHE_DIR.exists() and bool(list(CACHE_DIR.glob("*.sac")))
    if not cache_available:
        print(
            f"INFO: {CACHE_DIR} not populated.\n"
            "  Run: python3 tools/waveform_vv/refresh.py --event DPRK2017\n"
            "  Rendering synthetic-only comparison.",
            file=sys.stderr,
        )

    if not sac_files:
        print(
            f"WARN: no SAC files in {out_dir}; run the simulation first.",
            file=sys.stderr,
        )
        return 1

    fig, axes = plt.subplots(len(STATION_PAIRS), 1,
                             figsize=(10, 2.8 * len(STATION_PAIRS)),
                             sharex=False)
    if len(STATION_PAIRS) == 1:
        axes = [axes]

    for ax, (syn_key, obs_key, title) in zip(axes, STATION_PAIRS):
        # Find BHZ synthetic
        cands_syn = [f for f in sac_files
                     if syn_key in f.name and "BHZ" in f.name]
        if not cands_syn:
            cands_syn = [f for f in sac_files if syn_key in f.name]
        if not cands_syn:
            ax.text(0.5, 0.5, f"synthetic not found: {syn_key}",
                    ha="center", va="center", transform=ax.transAxes)
            ax.set_title(title)
            continue

        t_syn, d_syn = read_sac(cands_syn[0])
        ax.plot(t_syn, d_syn, color=PALETTE["synthetic"],
                lw=1.2, label="synthetic")
        peak = np.max(np.abs(d_syn))
        ax.axhline(peak, color=PALETTE["synthetic"],
                   lw=0.5, ls=":", alpha=0.5)
        ax.axhline(-peak, color=PALETTE["synthetic"],
                   lw=0.5, ls=":", alpha=0.5)

        cc_str = "n/a (no observed)"
        if cache_available:
            cands_obs = list(CACHE_DIR.glob(f"*{obs_key}*.sac"))
            if cands_obs:
                t_obs, d_obs = read_sac(cands_obs[0])
                ax.plot(t_obs, d_obs, color=PALETTE["observed"],
                        lw=1.0, alpha=0.75, label="observed (IRIS)")
                cc = normalized_xcorr(d_syn, d_obs)
                cc_str = f"{cc:.3f}"

        ax.set_title(f"{title}  |  cross-correlation = {cc_str}", fontsize=9)
        ax.set_ylabel("displ (m)")
        if ax is axes[0]:
            ax.legend(loc="upper right", fontsize=8)

    axes[-1].set_xlabel("t (s)")
    fig.suptitle(
        "DPRK 2017 (Punggye-ri, 6th test): synthetic vs observed\n"
        "0.5-5 Hz bandpass  |  ref: Voytan et al. 2019 GRL",
        y=1.00,
    )
    fig.tight_layout()

    out = Path(__file__).resolve().parents[1] / "05_synthetic_vs_observed.png"
    savefig(fig, out)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
