"""Shared matplotlib style for FSRM showcase figures.

Provides a fixed colorblind-safe palette (Wong 2011), consistent fonts and
sizing, and a small set of plot helpers used across the per-event showcase
figure scripts. Scripts that read simulation output should import from here
so that slides assembled from any combination of figures look like one
piece of work.

Usage:
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools" / "figures"))
    from figure_style import apply_style, PALETTE, FIG_WIDE, FIG_SINGLE
    apply_style()

The figure scripts read already-validated simulation output (HDF5, CSV,
SAC) and produce PNGs. They do not modify simulation output, do not
reimplement physics, and do not curve-fit results.
"""

from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt

# Wong 2011 colorblind-safe palette (DOI 10.1038/nmeth.1618).
# Names are semantic: pick the role, not the color.
PALETTE = {
    "synthetic": "#0072B2",   # blue
    "observed": "#D55E00",    # vermillion
    "analytic": "#009E73",    # bluish green
    "low_tier": "#CC79A7",    # reddish purple
    "med_tier": "#F0E442",    # yellow
    "high_tier": "#56B4E9",   # sky blue
    "highest_tier": "#0072B2", # blue (same as synthetic, used in pairs)
    "neutral": "#999999",     # gray (background, grids)
    "annotation": "#000000",  # black (text overlays)
}

FIG_SINGLE = (6.0, 4.0)   # single-panel default, inches
FIG_WIDE = (12.0, 4.0)    # multi-panel default
FIG_SQUARE = (6.0, 6.0)   # square (e.g. polar plots, station maps)

DPI = 300


def apply_style() -> None:
    """Apply the FSRM rcParams. Idempotent."""
    mpl.rcParams.update({
        # Fonts
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans"],
        "font.size": 11,
        "axes.titlesize": 12,
        "axes.labelsize": 11,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 9,
        # Grid
        "axes.grid": True,
        "grid.color": PALETTE["neutral"],
        "grid.linewidth": 0.4,
        "grid.alpha": 0.5,
        # Lines
        "lines.linewidth": 1.4,
        "lines.markersize": 4,
        # Spines
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
        # Saving
        "savefig.dpi": DPI,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,
        # Determinism
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    })


def plot_seismogram(ax, sac_path, label=None, color="synthetic",
                    band=None, normalize=False):
    """Plot a single SAC trace onto an axis.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    sac_path : path-like
    label : str | None
    color : str
        Key into PALETTE.
    band : tuple(float, float) | None
        (lo_hz, hi_hz). Optional bandpass filter.
    normalize : bool
        If True, divide by max abs before plotting.
    """
    from obspy import read

    st = read(str(sac_path))
    tr = st[0]
    if band is not None:
        tr.filter("bandpass", freqmin=band[0], freqmax=band[1],
                  corners=4, zerophase=True)
    data = tr.data
    if normalize:
        peak = abs(data).max()
        if peak > 0:
            data = data / peak
    t = tr.times()
    ax.plot(t, data, color=PALETTE[color], label=label or tr.id)


def plot_radial_profile(ax, h5_path, snapshot_index, field, color="med_tier"):
    """Plot one radial-state field at one snapshot index.

    The HDF5 file is the per-snapshot near-field profile written by the
    radial Lagrangian solver (`output/near_field_profile.h5`). The
    expected layout has a `radius_m` 1D array and per-snapshot
    field datasets under `snapshots/<index>/<field>`.
    """
    import h5py

    with h5py.File(str(h5_path), "r") as f:
        radius = f["radius_m"][:]
        snap_group = f[f"snapshots/{snapshot_index:04d}"]
        data = snap_group[field][:]
        t = snap_group.attrs.get("time_s", float("nan"))

    ax.plot(radius, data, color=PALETTE[color],
            label=f"t = {t:.2e} s" if not _isnan(t) else None)


def _isnan(x):
    try:
        return x != x
    except Exception:
        return False


def savefig(fig, path):
    """Save with deterministic settings."""
    fig.savefig(str(path), dpi=DPI, bbox_inches="tight",
                pad_inches=0.05, metadata={"Software": "FSRM"})


__all__ = [
    "PALETTE",
    "FIG_SINGLE",
    "FIG_WIDE",
    "FIG_SQUARE",
    "DPI",
    "apply_style",
    "plot_seismogram",
    "plot_radial_profile",
    "savefig",
]
