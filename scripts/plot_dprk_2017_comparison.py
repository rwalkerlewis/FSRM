#!/usr/bin/env python3
"""
DPRK 2017 Nuclear Test: Observed vs. Mueller-Murphy Synthetic Comparison

Produces a single publication-quality figure showing:
  - Left panel:  Real waveforms (record section, distance-sorted)
  - Right panel: Synthetic waveforms at matching station distances
  - Shared y-axis (distance in km), shared x-axis scale (time in seconds)

This script either:
  1. Loads real DPRK 2017 waveforms already downloaded to figures/dprk_2017/, OR
  2. Fetches them fresh from IRIS FDSN if not available.

Synthetics are generated using the MuellerMurphySource class from
mueller_murphy_demo.py with the Punggye-ri 1-D velocity model.

Usage:
    python scripts/plot_dprk_2017_comparison.py
"""

import os
import sys
import warnings

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import AutoMinorLocator

from obspy import UTCDateTime

warnings.filterwarnings("ignore")

# Add scripts directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fsrm.data_fetching import fetch_raw_waveforms, process_trace
from fsrm.signal_processing import time_relative
from fsrm.propagation import generate_synthetic, regional_travel_times, mt_mantap_collapse_signal
from fsrm.velocity_models import PunggyeRiVelocityModel
from mueller_murphy_demo import MuellerMurphySource


# ============================================================================
# Event parameters
# ============================================================================
EVENT_TIME = UTCDateTime("2017-09-03T03:30:01")
EVENT_LAT = 41.300
EVENT_LON = 129.076
EVENT_DEPTH_KM = 0.76
YIELD_KT = 250.0

FREQMIN = 0.5
FREQMAX = 8.0
PRE_ORIGIN = 60
POST_ORIGIN = 600

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.dirname(SCRIPT_DIR)
OUTDIR = os.path.join(REPO_DIR, "figures", "dprk_2017")

PHASE_VELOCITIES = {
    "Pn": 8.1,
    "Pg": 6.1,
    "Sn": 4.6,
    "Lg": 3.5,
}

STATIONS = [
    ("IC", "MDJ",  "*", "BH", "Mudanjiang, China"),
    ("KG", "INCN", "*", "BH", "Incheon, South Korea"),
    ("IC", "BJT",  "*", "BH", "Baijiatuan, China"),
    ("IC", "HIA",  "*", "BH", "Hailar, China"),
    ("IC", "SSE",  "*", "BH", "Shanghai, China"),
    ("IU", "MAJO", "*", "BH", "Matsushiro, Japan"),
    ("II", "ERM",  "*", "BH", "Erimo, Japan"),
    ("IU", "YSS",  "*", "BH", "Yuzhno, Russia"),
    ("IC", "ENH",  "*", "BH", "Enshi, China"),
    ("IU", "ULN",  "*", "BH", "Ulaanbaatar, Mongolia"),
    ("II", "AAK",  "*", "BH", "Ala Archa, Kyrgyzstan"),
]


# ============================================================================
# Phase colors and styles
# ============================================================================
PHASE_COLORS = {
    "Pn": "#1f77b4", "Pg": "#2ca02c",
    "Sn": "#d62728", "Lg": "#ff7f0e",
}
PHASE_STYLES = {
    "Pn": "--", "Pg": "-.",
    "Sn": "--", "Lg": "-.",
}


def fetch_observed_data():
    """Fetch and process real waveforms from IRIS FDSN."""
    print("  Downloading waveforms from IRIS FDSN ...")
    raw_results = fetch_raw_waveforms(
        EVENT_TIME, EVENT_LAT, EVENT_LON, STATIONS,
        pre_origin=PRE_ORIGIN, post_origin=POST_ORIGIN)

    if not raw_results:
        print("  ERROR: No waveforms retrieved.")
        return []

    processed = []
    for net, sta, dist_km, az, tr, inv, desc in raw_results:
        try:
            tr_proc = process_trace(tr, inv, FREQMIN, FREQMAX)
            processed.append((net, sta, dist_km, az, tr_proc, desc))
        except Exception as exc:
            print(f"  Failed to process {net}.{sta}: {exc}")

    return processed


def generate_synthetics_at_distances(distances_km, model):
    """Generate Mueller-Murphy synthetics at each observed station distance."""
    synthetics = []
    for dist_km in distances_km:
        t, v = generate_synthetic(
            dist_km, YIELD_KT, model,
            decoupling_factor=1.0, depth_km=EVENT_DEPTH_KM,
            dt=0.05, duration=POST_ORIGIN + PRE_ORIGIN,
            fmin=FREQMIN, fmax=FREQMAX,
            collapse_signal=mt_mantap_collapse_signal)
        synthetics.append((t - PRE_ORIGIN, v))
    return synthetics


def plot_comparison(observed, synthetics, distances_km, station_labels):
    """
    Produce the two-panel observed vs. synthetic comparison figure.

    Parameters:
        observed:       list of (times_array, data_array) per station
        synthetics:     list of (times_array, data_array) per station
        distances_km:   list of epicentral distances
        station_labels: list of "NET.STA" strings
    """
    n = len(distances_km)
    fig = plt.figure(figsize=(22, max(12, 1.5 * n)))
    gs = GridSpec(1, 2, figure=fig, wspace=0.08)
    ax_obs = fig.add_subplot(gs[0, 0])
    ax_syn = fig.add_subplot(gs[0, 1])

    trace_height = 0.38

    # Phase moveout lines in both panels
    d_line = np.linspace(0, max(distances_km) * 1.15, 300)
    for ax in [ax_obs, ax_syn]:
        for phase, vel in PHASE_VELOCITIES.items():
            col = PHASE_COLORS.get(phase, "gray")
            ls = PHASE_STYLES.get(phase, "--")
            t_line = d_line / vel
            mask = (t_line >= -PRE_ORIGIN) & (t_line <= POST_ORIGIN)
            if np.any(mask):
                ax.plot(t_line[mask], d_line[mask], color=col, ls=ls, lw=1.0,
                        alpha=0.6, zorder=10)

    # Plot observed waveforms
    for i, ((t_obs, d_obs), dist_km, label) in enumerate(
            zip(observed, distances_km, station_labels)):
        peak = np.max(np.abs(d_obs))
        if peak > 0:
            d_norm = d_obs / peak
        else:
            d_norm = d_obs

        # Distance-scaled spacing
        waveform_height = min(100, max(30, 50))
        d_plot = d_norm * waveform_height + dist_km

        ax_obs.fill_between(t_obs, d_plot, dist_km, where=(d_norm > 0),
                            color="black", alpha=0.3, zorder=4)
        ax_obs.fill_between(t_obs, d_plot, dist_km, where=(d_norm < 0),
                            color="gray", alpha=0.15, zorder=4)
        ax_obs.plot(t_obs, d_plot, "k-", lw=0.5, alpha=0.85, zorder=5)
        ax_obs.text(-PRE_ORIGIN - 5, dist_km, label,
                    fontsize=7.5, fontweight="bold", ha="right", va="center",
                    color="#333", clip_on=False)

    # Plot synthetic waveforms
    for i, ((t_syn, d_syn), dist_km, label) in enumerate(
            zip(synthetics, distances_km, station_labels)):
        peak = np.max(np.abs(d_syn))
        if peak > 0:
            d_norm = d_syn / peak
        else:
            d_norm = d_syn

        waveform_height = min(100, max(30, 50))
        d_plot = d_norm * waveform_height + dist_km

        ax_syn.fill_between(t_syn, d_plot, dist_km, where=(d_norm > 0),
                            color="#1f77b4", alpha=0.3, zorder=4)
        ax_syn.fill_between(t_syn, d_plot, dist_km, where=(d_norm < 0),
                            color="#aec7e8", alpha=0.2, zorder=4)
        ax_syn.plot(t_syn, d_plot, color="#1f77b4", lw=0.5, alpha=0.85, zorder=5)
        ax_syn.text(-PRE_ORIGIN - 5, dist_km, label,
                    fontsize=7.5, fontweight="bold", ha="right", va="center",
                    color="#333", clip_on=False)

    # Common axis formatting
    dist_min = min(distances_km)
    dist_max = max(distances_km)
    dist_pad = (dist_max - dist_min) * 0.08
    y_min = max(0, dist_min - dist_pad - 50)
    y_max = dist_max + dist_pad + 50

    for ax, title_text in [(ax_obs, "Observed (IRIS FDSN)"),
                           (ax_syn, "Mueller-Murphy Synthetic")]:
        ax.set_xlim(-PRE_ORIGIN, POST_ORIGIN)
        ax.set_ylim(y_min, y_max)
        ax.set_xlabel("Time after origin (s)", fontsize=11, fontweight="bold")
        ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        ax.grid(True, which="major", ls="-", alpha=0.15)
        ax.grid(True, which="minor", ls=":", alpha=0.08)
        ax.set_title(title_text, fontsize=13, fontweight="bold")

        # Phase legend
        for phase, vel in PHASE_VELOCITIES.items():
            col = PHASE_COLORS.get(phase, "gray")
            ls = PHASE_STYLES.get(phase, "--")
            ax.plot([], [], color=col, ls=ls, lw=1.5,
                    label=f"{phase} ({vel} km/s)")
        ax.legend(loc="upper right", fontsize=8, framealpha=0.9)

    ax_obs.set_ylabel("Epicentral Distance (km)", fontsize=11, fontweight="bold")
    ax_syn.set_yticklabels([])

    # Source parameters annotation
    src = MuellerMurphySource(YIELD_KT, rho=2700.0, vp=5600.0, depth_m=EVENT_DEPTH_KM * 1000)
    info_text = (
        f"Event: DPRK 2017-09-03 03:30:01 UTC\n"
        f"Yield: ~{YIELD_KT:.0f} kt (fully coupled)\n"
        f"DOB: {EVENT_DEPTH_KM * 1000:.0f} m in granite\n"
        f"mb(pred): {src.mb_from_yield(YIELD_KT):.2f}  |  mb(obs): 6.3\n"
        f"fc: {src.corner_frequency:.3f} Hz\n"
        f"Rc: {src.cavity_radius:.0f} m\n"
        f"Filter: {FREQMIN}-{FREQMAX} Hz BP"
    )
    fig.text(0.5, 0.02, info_text, ha="center", va="bottom",
             fontsize=9, fontfamily="monospace",
             bbox=dict(boxstyle="round,pad=0.5", fc="#fffff0", ec="#cca",
                       alpha=0.9))

    fig.suptitle(
        "DPRK 2017 Nuclear Test: Observed vs. Mueller-Murphy Synthetic\n"
        f"BHZ vertical component, {FREQMIN}-{FREQMAX} Hz bandpass, "
        f"Punggye-ri 1-D velocity model",
        fontsize=14, fontweight="bold", y=0.98)

    plt.tight_layout(rect=[0, 0.08, 1, 0.95])
    outpath = os.path.join(OUTDIR, "observed_vs_synthetic.png")
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {outpath}")


def main():
    os.makedirs(OUTDIR, exist_ok=True)

    print("=" * 72)
    print("  DPRK 2017: Observed vs. Mueller-Murphy Synthetic Comparison")
    print("=" * 72)
    print()

    # Step 1: Get observed data
    print("[1/3] Fetching observed waveforms ...")
    observed_raw = fetch_observed_data()
    if not observed_raw:
        print("ERROR: Could not retrieve any observed data. Exiting.")
        sys.exit(1)
    print(f"       Got {len(observed_raw)} station(s)")
    print()

    # Extract time series and metadata
    distances_km = []
    station_labels = []
    observed_traces = []
    for net, sta, dist_km, az, tr, desc in observed_raw:
        times = time_relative(tr, EVENT_TIME)
        distances_km.append(dist_km)
        station_labels.append(f"{net}.{sta}")
        observed_traces.append((times, tr.data.copy()))

    # Step 2: Generate synthetics at matching distances
    print("[2/3] Generating Mueller-Murphy synthetics ...")
    model = PunggyeRiVelocityModel()
    synthetics = generate_synthetics_at_distances(distances_km, model)
    print(f"       Generated {len(synthetics)} synthetic trace(s)")
    print()

    # Step 3: Plot comparison
    print("[3/3] Plotting observed vs. synthetic comparison ...")
    plot_comparison(observed_traces, synthetics, distances_km, station_labels)
    print()

    print("=" * 72)
    print("  Done. Figure saved to:", OUTDIR)
    print("=" * 72)


if __name__ == "__main__":
    main()
