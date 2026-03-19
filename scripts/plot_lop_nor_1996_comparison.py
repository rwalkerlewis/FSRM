#!/usr/bin/env python3
"""
Lop Nor 1996 Final Chinese Nuclear Test: Observed vs. Synthetic Comparison

Produces a publication-quality multi-panel figure showing:
  - Left panel:  Real waveforms from KNDC archive (record section, distance-sorted)
  - Right panel: Mueller-Murphy synthetic waveforms at matching station distances
  - Additional panels: overlay comparison for key stations, spectral analysis

The KNDC archive provides 26 short-period (SKM/SKD) Z-component recordings
from Kazakhstan and Kyrgyzstan stations at distances of 761–1531 km.

Synthetics are generated using the Mueller-Murphy RDP model with the
Lop Nor 1-D velocity model (Kuruktag / Tarim Basin, Moho 48 km).

Event:  1996-07-29  01:48:57.170 UTC, Lop Nor (~41.7163°N, 88.3748°E)
        45th and final confirmed Chinese underground nuclear test
        mb 4.90 (ISC), estimated yield ~5 kt (coupled in granite)

Usage:
    python scripts/plot_lop_nor_1996_comparison.py
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

from obspy import UTCDateTime, read
from obspy.geodetics import gps2dist_azimuth

warnings.filterwarnings("ignore")

# Add scripts directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fsrm.signal_processing import time_relative, spectral_amplitude, envelope
from fsrm.propagation import generate_synthetic, regional_travel_times
from fsrm.velocity_models import LopNorVelocityModel
from fsrm.source_physics import (
    cavity_radius_empirical, corner_frequency_patton, mb_from_yield,
)


# ============================================================================
# Event parameters
# ============================================================================
EVENT_TIME = UTCDateTime("1996-07-29T01:48:57.170")
EVENT_LAT = 41.7163
EVENT_LON = 88.3748
EVENT_DEPTH_KM = 0.4
YIELD_KT = 5.0

FREQMIN = 0.5
FREQMAX = 8.0
PRE_ORIGIN = 60
POST_ORIGIN = 500

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.dirname(SCRIPT_DIR)
OUTDIR = os.path.join(REPO_DIR, "figures", "lop_nor_1996")

PHASE_VELOCITIES = {
    "Pn": 7.9,
    "Pg": 6.1,
    "Sn": 4.5,
    "Lg": 3.5,
}

PHASE_COLORS = {
    "Pn": "#1f77b4", "Pg": "#2ca02c",
    "Sn": "#d62728", "Lg": "#ff7f0e",
}
PHASE_STYLES = {
    "Pn": "--", "Pg": "-.",
    "Sn": "--", "Lg": "-.",
}


def _parse_site_coords():
    """Parse KNDC .site file to build station -> (lat, lon) lookup."""
    site_file = os.path.join(REPO_DIR, "data", "historic_nuclear",
                             "03.Lopnor", "19960729-0148.site")
    coords = {}
    if not os.path.isfile(site_file):
        return coords
    event_jdate = 1996211
    with open(site_file, "r") as f:
        for line in f:
            parts = line.split()
            if len(parts) < 6:
                continue
            sta = parts[0].strip()
            try:
                ondate = int(parts[1])
                offdate = int(parts[2]) if parts[2] != "-1" else 9999999
                lat = float(parts[3])
                lon = float(parts[4])
            except (ValueError, IndexError):
                continue
            if ondate <= event_jdate <= offdate:
                if sta not in coords:
                    coords[sta] = (lat, lon)
    return coords


def load_kndc_observed():
    """
    Load Z-component SAC files from the KNDC archive.

    Returns list of (net, sta, dist_km, az, trace, description) tuples.
    """
    data_dir = os.path.join(REPO_DIR, "data", "historic_nuclear",
                            "03.Lopnor", "19960729.0148", "wf")
    if not os.path.isdir(data_dir):
        print(f"  WARNING: KNDC directory not found: {data_dir}")
        return []

    site_coords = _parse_site_coords()

    results = []
    sac_files = sorted(f for f in os.listdir(data_dir) if f.endswith(".sac"))
    z_files = [f for f in sac_files if "z1" in f.lower() or "Z1" in f]

    for fname in z_files:
        fpath = os.path.join(data_dir, fname)
        try:
            st = read(fpath, format="SAC")
            tr = st[0]

            sta = tr.stats.sac.get("kstnm", "").strip()
            if not sta:
                base = fname.replace(".sac", "")
                if base[:5].isdigit():
                    sta = base[5:].replace("z1", "").replace("z2", "").upper()
                else:
                    idx = base.find("96211")
                    if idx > 0:
                        sta = base[:idx].upper()
                    else:
                        sta = base[:4].upper()

            stla = tr.stats.sac.get("stla", None)
            stlo = tr.stats.sac.get("stlo", None)
            if stla is not None and stlo is not None and stla != 0.0:
                dist_m, az, _ = gps2dist_azimuth(EVENT_LAT, EVENT_LON, stla, stlo)
                dist_km = dist_m / 1000.0
            elif sta.upper() in site_coords:
                stla, stlo = site_coords[sta.upper()]
                dist_m, az, _ = gps2dist_azimuth(EVENT_LAT, EVENT_LON, stla, stlo)
                dist_km = dist_m / 1000.0
            else:
                dist_km = tr.stats.sac.get("dist", 0.0)
                az = tr.stats.sac.get("az", 0.0)
                if dist_km == 0.0:
                    continue

            # Basic processing
            tr.detrend("demean")
            tr.detrend("linear")
            tr.taper(max_percentage=0.05, type="cosine")
            tr.filter("bandpass", freqmin=FREQMIN, freqmax=FREQMAX,
                       corners=4, zerophase=True)

            results.append(("KZ", sta, dist_km, az, tr, f"{sta} ({dist_km:.0f} km)"))
        except Exception as exc:
            print(f"  WARNING: Could not read {fname}: {exc}")

    results.sort(key=lambda x: x[2])
    return results


def generate_synthetics_at_distances(distances_km, model):
    """Generate Mueller-Murphy synthetics at each observed station distance."""
    synthetics = []
    for dist_km in distances_km:
        t, v = generate_synthetic(
            dist_km, YIELD_KT, model,
            decoupling_factor=1.0, depth_km=EVENT_DEPTH_KM,
            dt=0.05, duration=POST_ORIGIN + PRE_ORIGIN,
            fmin=FREQMIN, fmax=FREQMAX)
        synthetics.append((t - PRE_ORIGIN, v))
    return synthetics


def plot_record_section_comparison(observed, synthetics, distances_km, station_labels):
    """
    Two-panel observed vs. synthetic record section.
    """
    n = len(distances_km)
    fig = plt.figure(figsize=(22, max(12, 1.2 * n)))
    gs = GridSpec(1, 2, figure=fig, wspace=0.08)
    ax_obs = fig.add_subplot(gs[0, 0])
    ax_syn = fig.add_subplot(gs[0, 1])

    # Phase moveout lines
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
        d_norm = d_obs / peak if peak > 0 else d_obs

        waveform_height = min(80, max(25, 40))
        d_plot = d_norm * waveform_height + dist_km

        ax_obs.fill_between(t_obs, d_plot, dist_km, where=(d_norm > 0),
                            color="black", alpha=0.3, zorder=4)
        ax_obs.fill_between(t_obs, d_plot, dist_km, where=(d_norm < 0),
                            color="gray", alpha=0.15, zorder=4)
        ax_obs.plot(t_obs, d_plot, "k-", lw=0.5, alpha=0.85, zorder=5)
        ax_obs.text(-PRE_ORIGIN - 5, dist_km, label,
                    fontsize=6.5, fontweight="bold", ha="right", va="center",
                    color="#333", clip_on=False)

    # Plot synthetic waveforms
    for i, ((t_syn, d_syn), dist_km, label) in enumerate(
            zip(synthetics, distances_km, station_labels)):
        peak = np.max(np.abs(d_syn))
        d_norm = d_syn / peak if peak > 0 else d_syn

        waveform_height = min(80, max(25, 40))
        d_plot = d_norm * waveform_height + dist_km

        ax_syn.fill_between(t_syn, d_plot, dist_km, where=(d_norm > 0),
                            color="#1f77b4", alpha=0.3, zorder=4)
        ax_syn.fill_between(t_syn, d_plot, dist_km, where=(d_norm < 0),
                            color="#aec7e8", alpha=0.2, zorder=4)
        ax_syn.plot(t_syn, d_plot, color="#1f77b4", lw=0.5, alpha=0.85, zorder=5)
        ax_syn.text(-PRE_ORIGIN - 5, dist_km, label,
                    fontsize=6.5, fontweight="bold", ha="right", va="center",
                    color="#333", clip_on=False)

    # Axis formatting
    dist_min = min(distances_km)
    dist_max = max(distances_km)
    dist_pad = (dist_max - dist_min) * 0.08
    y_min = max(0, dist_min - dist_pad - 40)
    y_max = dist_max + dist_pad + 40

    for ax, title_text in [(ax_obs, "Observed (KNDC Archive)"),
                           (ax_syn, "Mueller-Murphy Synthetic")]:
        ax.set_xlim(-PRE_ORIGIN, POST_ORIGIN)
        ax.set_ylim(y_min, y_max)
        ax.set_xlabel("Time after origin (s)", fontsize=11, fontweight="bold")
        ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        ax.grid(True, which="major", ls="-", alpha=0.15)
        ax.grid(True, which="minor", ls=":", alpha=0.08)
        ax.set_title(title_text, fontsize=13, fontweight="bold")

        for phase, vel in PHASE_VELOCITIES.items():
            col = PHASE_COLORS.get(phase, "gray")
            ls = PHASE_STYLES.get(phase, "--")
            ax.plot([], [], color=col, ls=ls, lw=1.5,
                    label=f"{phase} ({vel} km/s)")
        ax.legend(loc="upper right", fontsize=8, framealpha=0.9)

    ax_obs.set_ylabel("Epicentral Distance (km)", fontsize=11, fontweight="bold")
    ax_syn.set_yticklabels([])

    # Info annotation
    fc = corner_frequency_patton(YIELD_KT)
    Rc = cavity_radius_empirical(YIELD_KT, rho=2650.0)
    mb_pred = mb_from_yield(YIELD_KT, 1.0)
    info_text = (
        f"Event: Lop Nor 1996-07-29 01:48:57 UTC\n"
        f"45th & Final Chinese Underground Nuclear Test\n"
        f"Yield: ~{YIELD_KT:.0f} kt (coupled, from mb-yield)\n"
        f"DOB: {EVENT_DEPTH_KM * 1000:.0f} m in granite tunnel\n"
        f"mb(pred): {mb_pred:.2f}  |  mb(obs): 4.90\n"
        f"fc: {fc:.3f} Hz  |  Rc: {Rc:.0f} m\n"
        f"Filter: {FREQMIN}\u2013{FREQMAX} Hz BP, Lop Nor 1-D model"
    )
    fig.text(0.5, 0.02, info_text, ha="center", va="bottom",
             fontsize=9, fontfamily="monospace",
             bbox=dict(boxstyle="round,pad=0.5", fc="#fffff0", ec="#cca",
                       alpha=0.9))

    fig.suptitle(
        "Lop Nor 1996 Final Chinese Nuclear Test: Observed vs. Synthetic\n"
        f"BHZ vertical component, {FREQMIN}\u2013{FREQMAX} Hz bandpass, "
        f"26 KNDC regional stations (761\u20131531 km)",
        fontsize=14, fontweight="bold", y=0.98)

    plt.tight_layout(rect=[0, 0.08, 1, 0.95])
    outpath = os.path.join(OUTDIR, "observed_vs_synthetic.png")
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {outpath}")


def plot_overlay_key_stations(observed_data, model):
    """
    Overlay comparison for key stations — observed (black) vs synthetic (blue).
    """
    key_stas = ["DZHR", "UZB", "PRZ", "SATY", "KURM", "TRG"]
    selected = [(net, sta, dist, az, tr, desc)
                for net, sta, dist, az, tr, desc in observed_data
                if sta in key_stas]

    if not selected:
        print("  No key stations found for overlay plot")
        return

    n = len(selected)
    fig, axes = plt.subplots(n, 2, figsize=(18, 3.0 * n), squeeze=False)

    for i, (net, sta, dist_km, az, tr, desc) in enumerate(selected):
        times_obs = time_relative(tr, EVENT_TIME)
        data_obs = tr.data.copy()

        # Generate synthetic at same distance
        t_syn, v_syn = generate_synthetic(
            dist_km, YIELD_KT, model,
            decoupling_factor=1.0, depth_km=EVENT_DEPTH_KM,
            dt=0.05, duration=600,
            fmin=FREQMIN, fmax=FREQMAX)

        # Left panel: waveform overlay
        ax = axes[i, 0]
        peak_obs = np.max(np.abs(data_obs))
        peak_syn = np.max(np.abs(v_syn))
        if peak_obs > 0:
            obs_norm = data_obs / peak_obs
        else:
            obs_norm = data_obs
        if peak_syn > 0:
            syn_norm = v_syn / peak_syn
        else:
            syn_norm = v_syn

        ax.plot(times_obs, obs_norm, "k-", lw=0.6, alpha=0.8, label="Observed (KNDC)")
        ax.plot(t_syn, syn_norm, "b-", lw=0.6, alpha=0.7, label="Synthetic (M-M)")

        # Phase arrival markers
        tt = regional_travel_times(dist_km, model)
        for ph, (tarr, _, _) in tt.items():
            if 0 < tarr < 500:
                ax.axvline(tarr, color=PHASE_COLORS.get(ph, "gray"),
                           ls="--", lw=0.8, alpha=0.5)

        ax.set_xlim(-30, 400)
        ax.set_ylim(-1.4, 1.4)
        ax.set_ylabel("Normalised velocity", fontsize=8)
        ax.set_title(f"{sta} ({dist_km:.0f} km, az {az:.0f}\u00b0)",
                     fontsize=10, fontweight="bold", loc="left")
        ax.legend(fontsize=7, loc="upper right")
        ax.grid(True, alpha=0.15)

        # Right panel: spectral comparison
        ax = axes[i, 1]
        freq_obs, spec_obs = spectral_amplitude(data_obs, tr.stats.delta)
        freq_syn, spec_syn = spectral_amplitude(v_syn, 0.05)

        m_obs = (freq_obs > 0.2) & (freq_obs < FREQMAX * 1.5)
        m_syn = (freq_syn > 0.2) & (freq_syn < FREQMAX * 1.5)

        smax_obs = np.max(spec_obs[m_obs]) if np.any(m_obs) and np.max(spec_obs[m_obs]) > 0 else 1.0
        smax_syn = np.max(spec_syn[m_syn]) if np.any(m_syn) and np.max(spec_syn[m_syn]) > 0 else 1.0

        ax.semilogy(freq_obs[m_obs], spec_obs[m_obs] / smax_obs,
                     "k-", lw=1, alpha=0.7, label="Observed")
        ax.semilogy(freq_syn[m_syn], spec_syn[m_syn] / smax_syn,
                     "b-", lw=1, alpha=0.7, label="Synthetic")

        fc = corner_frequency_patton(YIELD_KT)
        ax.axvline(fc, color="red", ls=":", lw=1, alpha=0.6,
                   label=f"fc = {fc:.2f} Hz")

        ax.set_xlabel("Frequency (Hz)", fontsize=8)
        ax.set_ylabel("Normalised amplitude", fontsize=8)
        ax.set_title(f"Spectrum: {sta}", fontsize=10, fontweight="bold", loc="left")
        ax.legend(fontsize=7, loc="upper right")
        ax.grid(True, alpha=0.15, which="both")
        ax.set_xlim(0.2, FREQMAX * 1.5)

    axes[-1, 0].set_xlabel("Time after origin (s)", fontsize=9)
    axes[-1, 1].set_xlabel("Frequency (Hz)", fontsize=9)

    fig.suptitle(
        "Lop Nor 1996: Observed vs Synthetic — Key Station Overlays\n"
        f"Mueller-Murphy {YIELD_KT:.0f} kt, {FREQMIN}\u2013{FREQMAX} Hz",
        fontsize=13, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    outpath = os.path.join(OUTDIR, "overlay_key_stations.png")
    fig.savefig(outpath, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {outpath}")


def plot_spectral_summary(observed_data, model):
    """
    Multi-station spectral comparison with Mueller-Murphy source spectrum overlay.
    """
    fig, axes = plt.subplots(2, 1, figsize=(16, 12))

    # (a) All-station spectral comparison
    ax = axes[0]
    cmap = plt.cm.viridis
    n = len(observed_data)
    for i, (net, sta, dist_km, az, tr, desc) in enumerate(observed_data):
        freq, spec = spectral_amplitude(tr.data, tr.stats.delta)
        m = (freq > 0.2) & (freq < FREQMAX * 1.5)
        if np.max(spec[m]) > 0:
            spec_n = spec[m] / np.max(spec[m])
            col = cmap(i / max(n - 1, 1))
            ax.semilogy(freq[m], spec_n, "-", color=col, lw=0.8, alpha=0.6,
                        label=f"{sta} ({dist_km:.0f} km)" if i < 8 else "")

    # Mueller-Murphy source spectrum envelope
    from fsrm.source_physics import mueller_murphy_spectrum
    ff = np.logspace(-1, 1, 200)
    mm_spec = mueller_murphy_spectrum(ff, YIELD_KT)
    if np.max(mm_spec) > 0:
        mm_n = mm_spec / np.max(mm_spec)
        ax.semilogy(ff, mm_n, "r-", lw=2.5, alpha=0.8,
                    label=f"Mueller-Murphy source ({YIELD_KT:.0f} kt)")

    ax.set_xlabel("Frequency (Hz)", fontsize=11)
    ax.set_ylabel("Normalised spectral amplitude", fontsize=11)
    ax.set_title("(a) Observed spectra (all stations) vs Mueller-Murphy source",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=7, ncol=3, loc="upper right")
    ax.grid(True, alpha=0.2, which="both")
    ax.set_xlim(0.2, 15)

    # (b) Distance vs spectral slope (discrimination)
    ax = axes[1]
    distances = []
    spectral_slopes = []
    for net, sta, dist_km, az, tr, desc in observed_data:
        freq, spec = spectral_amplitude(tr.data, tr.stats.delta)
        m = (freq > 1.0) & (freq < 6.0)
        if np.any(m) and np.max(spec[m]) > 0:
            log_f = np.log10(freq[m])
            log_s = np.log10(spec[m] / np.max(spec[m]) + 1e-10)
            if len(log_f) > 2:
                slope = np.polyfit(log_f, log_s, 1)[0]
                distances.append(dist_km)
                spectral_slopes.append(slope)

    if distances:
        ax.scatter(distances, spectral_slopes, c="r", s=50, alpha=0.7,
                   label="KNDC stations", zorder=5)
        ax.axhline(-2.0, color="blue", ls="--", lw=1.5, alpha=0.5,
                   label="Mueller-Murphy f\u207b\u00b2 rolloff")
        ax.axhline(-3.0, color="green", ls="--", lw=1.5, alpha=0.5,
                   label="Brune earthquake f\u207b\u00b3 rolloff")

    ax.set_xlabel("Epicentral Distance (km)", fontsize=11)
    ax.set_ylabel("Spectral slope (log-log)", fontsize=11)
    ax.set_title("(b) High-frequency spectral slope vs distance — discrimination diagnostic",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.2)

    fig.suptitle(
        "Lop Nor 1996 Final Test: Spectral Analysis\n"
        f"26 KNDC stations, {FREQMIN}\u2013{FREQMAX} Hz",
        fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    outpath = os.path.join(OUTDIR, "spectral_analysis.png")
    fig.savefig(outpath, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {outpath}")


def main():
    os.makedirs(OUTDIR, exist_ok=True)

    print("=" * 72)
    print("  Lop Nor 1996: Observed vs. Mueller-Murphy Synthetic Comparison")
    print("  45th & Final Chinese Underground Nuclear Test")
    print("=" * 72)
    print()

    # Step 1: Load observed data from KNDC archive
    print("[1/4] Loading KNDC archive waveforms ...")
    observed_data = load_kndc_observed()
    if not observed_data:
        print("ERROR: Could not load any KNDC data. Exiting.")
        sys.exit(1)
    print(f"       Loaded {len(observed_data)} station(s)")
    print()

    # Extract arrays for record section
    distances_km = []
    station_labels = []
    observed_traces = []
    for net, sta, dist_km, az, tr, desc in observed_data:
        times = time_relative(tr, EVENT_TIME)
        distances_km.append(dist_km)
        station_labels.append(sta)
        observed_traces.append((times, tr.data.copy()))

    # Step 2: Generate synthetics
    print("[2/4] Generating Mueller-Murphy synthetics ...")
    model = LopNorVelocityModel()
    synthetics = generate_synthetics_at_distances(distances_km, model)
    print(f"       Generated {len(synthetics)} synthetic trace(s)")
    print()

    # Step 3: Record section comparison
    print("[3/4] Plotting record section comparison ...")
    plot_record_section_comparison(observed_traces, synthetics,
                                   distances_km, station_labels)

    # Step 4: Key station overlays
    print("[4/4] Plotting key station overlays and spectral analysis ...")
    plot_overlay_key_stations(observed_data, model)
    plot_spectral_summary(observed_data, model)

    print()
    print("=" * 72)
    print(f"  All figures saved to: {OUTDIR}")
    print("=" * 72)


if __name__ == "__main__":
    main()
