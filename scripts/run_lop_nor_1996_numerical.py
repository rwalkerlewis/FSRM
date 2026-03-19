#!/usr/bin/env python3
"""
Lop Nor 1996 — Numerical Wave Propagation: Setup, Validation & Comparison

This script manages the full pipeline for the FSRM numerical elastic wave
propagation simulation of the 1996-07-29 Lop Nor nuclear test (~5 kt coupled):

  Phase 1: Config validation (CFL, memory, resolution, wall-clock estimates)
  Phase 2: Print the mpirun invocation command for the FSRM C++ solver
  Phase 3: Post-processing of FSRM numerical output (if simulation has run)
           → Record section of numerical seismograms
           → Numerical vs analytical (Mueller-Murphy) waveform comparison
           → Spectral comparison (6-panel)
           → Amplitude decay vs geometric spreading
  Phase 4: Far-field extrapolation from domain edge (200 km) to real KNDC
           station distances (761–1531 km) using Q-corrected spreading
  Phase 5: Comparison against observed KNDC archive waveforms
           → Extrapolated numerical vs observed record section
           → Per-station waveform overlay (key stations)
           → Spectral analysis comparison

The FSRM C++ simulation itself is NOT executed by this script; it requires
HPC resources.  Run the solver first, then re-run this script for Phase 3–5.

Usage:
    python scripts/run_lop_nor_1996_numerical.py

Event:  1996-07-29 01:48:57.170 UTC — Lop Nor, China
        45th and final confirmed Chinese underground nuclear test
        ~5 kt coupled in granite tunnel, mb 4.90 (ISC)
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

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, SCRIPT_DIR)

from fsrm.run_utils import (
    validate_config, print_invocation, postprocess_far_field,
    read_sac_seismograms, read_mseed_seismograms,
    extrapolate_to_far_field, get_trace_distance,
)
from fsrm.source_physics import (
    mueller_murphy_spectrum, corner_frequency_patton, cavity_radius_empirical,
    scalar_moment_coupled, mb_from_yield, Mw_from_M0,
)
from fsrm.signal_processing import spectral_amplitude, time_relative, envelope
from fsrm.propagation import generate_synthetic, regional_travel_times
from fsrm.velocity_models import LopNorVelocityModel


# ═══════════════════════════════════════════════════════════════════════════════
# Event Parameters
# ═══════════════════════════════════════════════════════════════════════════════
EVENT_TIME = UTCDateTime("1996-07-29T01:48:57.170")
EVENT_LAT = 41.7163
EVENT_LON = 88.3748
EVENT_DEPTH_KM = 0.4
YIELD_KT = 5.0

FMIN, FMAX = 0.5, 4.0    # Analysis band for numerical (conservative)

CONFIG_PATH = os.path.join(PROJECT_ROOT, "config",
                           "lop_nor_1996_numerical_wave_propagation.config")
SIM_OUTPUT_DIR = os.path.join(PROJECT_ROOT, "output",
                              "lop_nor_1996_numerical", "seismograms")
OUTDIR = os.path.join(PROJECT_ROOT, "figures", "lop_nor_1996_numerical")

EVENT_LABEL = "Lop Nor 1996-07-29, ~5 kt coupled, Kuruktag granite"

GRID = dict(
    nx=800, ny=800, nz=120,
    dx=500.0,
    vp_max=7900.0,
    domain_km=400.0,
    cfl_target=0.45,
    dt_suggested=0.02,
    nprocs_min=32, nprocs_max=64,
    duration=90.0,
)

PHASE_VELOCITIES = {
    "Pn": 7.9, "Pg": 5.8, "Sn": 4.5, "Lg": 3.5,
}
PHASE_COLORS = {
    "Pn": "#1f77b4", "Pg": "#2ca02c", "Sn": "#d62728", "Lg": "#ff7f0e",
}

# KNDC station real distances and domain-edge proxy distances (all at 150 km
# in the simulation domain, matching azimuth)
KNDC_STATION_MAP = {
    # sta: (real_dist_km, domain_dist_km, azimuth_deg)
    "DZHR": (761.0,  150.0, 315.0),
    "UZB":  (785.0,  150.0, 308.0),
    "PRZ":  (829.0,  150.0, 285.0),
    "KURM": (860.0,  150.0, 316.0),
    "TRG":  (899.0,  150.0, 307.0),
    "KUU":  (1013.0, 150.0, 323.0),
    "BRD":  (1531.0, 150.0, 291.0),
    "YUZH": (1519.0, 150.0, 280.0),
}

# Average crustal Q and velocity for far-field extrapolation
Q_AVG = 500.0         # Average intrinsic Q along central Asian paths
VP_AVG = 6.5          # Average crustal P velocity (km/s)


# ═══════════════════════════════════════════════════════════════════════════════
# KNDC Archive Loading (observed data)
# ═══════════════════════════════════════════════════════════════════════════════

def _parse_site_coords():
    """Parse KNDC .site file to build station -> (lat, lon) lookup."""
    site_file = os.path.join(PROJECT_ROOT, "data", "historic_nuclear",
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
    """Load Z-component SAC files from the KNDC archive."""
    data_dir = os.path.join(PROJECT_ROOT, "data", "historic_nuclear",
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
                continue

            tr.detrend("demean")
            tr.detrend("linear")
            tr.taper(max_percentage=0.05, type="cosine")
            tr.filter("bandpass", freqmin=FMIN, freqmax=FMAX,
                       corners=4, zerophase=True)

            results.append(("KZ", sta, dist_km, az, tr, f"{sta} (KNDC, {dist_km:.0f} km)"))
        except Exception as exc:
            print(f"  WARNING: Could not read {fname}: {exc}")

    results.sort(key=lambda x: x[2])
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# Comparison Figures (Numerical vs Observed)
# ═══════════════════════════════════════════════════════════════════════════════

def fig_extrapolated_record_section(extrap_traces, observed_data, model):
    """Two-panel record section: extrapolated numerical vs observed."""
    n_num = len(extrap_traces)
    n_obs = len(observed_data)
    if n_num == 0:
        print("  No extrapolated traces — skipping record section")
        return

    fig = plt.figure(figsize=(24, max(14, 1.5 * max(n_num, n_obs))))
    gs = GridSpec(1, 2, figure=fig, wspace=0.08)

    for panel_idx, (dataset, title_str) in enumerate([
        (extrap_traces, "FSRM Numerical (extrapolated)"),
        (observed_data, "Observed (KNDC archive)"),
    ]):
        ax = fig.add_subplot(gs[0, panel_idx])
        n = len(dataset)
        if n == 0:
            ax.text(0.5, 0.5, "No data", ha="center", va="center",
                    fontsize=14, transform=ax.transAxes)
            ax.set_title(title_str, fontsize=12, fontweight="bold")
            continue

        if panel_idx == 0:
            # Numerical: distance is stored in stats
            dists = [get_trace_distance(tr) for tr in dataset]
            order = np.argsort(dists)
            for yi, idx in enumerate(order):
                tr = dataset[idx]
                dist_km = dists[idx]
                t = np.arange(tr.stats.npts) * tr.stats.delta
                d = tr.data.copy().astype(float)
                peak = np.max(np.abs(d)) if len(d) > 0 else 1.0
                if peak > 0:
                    d /= peak
                d *= 0.4
                ax.plot(t, d + yi, "k-", lw=0.35, alpha=0.85)
                ax.fill_between(t, d + yi, yi, where=(d > 0),
                                color="#444", alpha=0.10)
                sta = getattr(tr.stats, "station", "?")
                ax.text(-5, yi, f"{sta} ({dist_km:.0f} km)",
                        fontsize=7, ha="right", va="center", fontweight="bold",
                        clip_on=False)
                # Phase markers
                for phase, vel in PHASE_VELOCITIES.items():
                    t_arr = dist_km / vel
                    if 0 < t_arr < t[-1]:
                        ax.axvline(t_arr, ymin=(yi - 0.3) / n, ymax=(yi + 0.3) / n,
                                   color=PHASE_COLORS.get(phase, "gray"),
                                   ls="--", lw=0.6, alpha=0.4)
            ax.set_ylim(-0.8, len(dataset) - 0.2)
        else:
            # Observed
            for yi, (net, sta, dist_km, az, tr, desc) in enumerate(dataset):
                times = time_relative(tr, EVENT_TIME)
                d = tr.data.copy()
                peak = np.max(np.abs(d)) if len(d) > 0 else 1.0
                if peak > 0:
                    d /= peak
                d *= 0.4
                ax.plot(times, d + yi, "k-", lw=0.35, alpha=0.85)
                ax.fill_between(times, d + yi, yi, where=(d > 0),
                                color="#444", alpha=0.10)
                ax.text(-5, yi, f"{sta} ({dist_km:.0f} km)",
                        fontsize=7, ha="right", va="center", fontweight="bold",
                        clip_on=False)
            ax.set_ylim(-0.8, n - 0.2)

        ax.set_xlabel("Time (s)", fontsize=11, fontweight="bold")
        ax.set_title(title_str, fontsize=12, fontweight="bold")
        ax.set_yticks([])
        ax.grid(True, axis="x", alpha=0.2)
        ax.xaxis.set_minor_locator(AutoMinorLocator(5))

    fig.suptitle(f"FSRM Numerical vs Observed — {EVENT_LABEL}\n"
                 f"Numerical extrapolated from 150 km to real station distances, "
                 f"Q={Q_AVG:.0f}",
                 fontsize=13, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.94])

    path = os.path.join(OUTDIR, "fig_numerical_vs_observed_record_section.png")
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {path}")


def fig_station_overlay(extrap_traces, observed_data, model):
    """Per-station waveform + spectral overlay for matching KNDC stations."""
    # Build lookup of extrapolated numerical traces by station name
    num_by_sta = {}
    for tr in extrap_traces:
        sta = getattr(tr.stats, "station", "")
        if sta:
            num_by_sta[sta] = tr

    # Build lookup of observed traces
    obs_by_sta = {}
    for net, sta, dist_km, az, tr, desc in observed_data:
        if sta in KNDC_STATION_MAP:
            obs_by_sta[sta] = (dist_km, tr)

    # Find matching stations
    common = sorted(set(num_by_sta.keys()) & set(obs_by_sta.keys()))
    n = len(common)
    if n == 0:
        print("  No matching stations for overlay — skipping")
        return

    fig, axes = plt.subplots(n, 2, figsize=(22, 3.5 * n), squeeze=False)
    TMIN, TMAX = -20, 400

    for i, sta in enumerate(common):
        tr_num = num_by_sta[sta]
        dist_real, tr_obs = obs_by_sta[sta]

        # Numerical waveform
        t_num = np.arange(tr_num.stats.npts) * tr_num.stats.delta
        d_num = tr_num.data.copy().astype(float)
        peak_num = np.max(np.abs(d_num)) if len(d_num) > 0 else 1.0
        if peak_num > 0:
            d_num /= peak_num

        # Observed waveform
        times_obs = time_relative(tr_obs, EVENT_TIME)
        d_obs = tr_obs.data.copy()
        mask = (times_obs >= TMIN) & (times_obs <= TMAX)
        t_obs = times_obs[mask]
        d_obs = d_obs[mask]
        peak_obs = np.max(np.abs(d_obs)) if len(d_obs) > 0 else 1.0
        if peak_obs > 0:
            d_obs /= peak_obs

        # Waveform overlay
        ax = axes[i][0]
        ax.plot(t_obs, d_obs, "k-", lw=0.5, alpha=0.7, label="Observed (KNDC)")
        ax.plot(t_num, d_num, "r-", lw=0.6, alpha=0.6,
                label="FSRM numerical (extrap.)")

        # Phase markers
        tt = regional_travel_times(dist_real, model)
        for ph, (tarr, _, _) in tt.items():
            if TMIN < tarr < TMAX:
                ax.axvline(tarr, color=PHASE_COLORS.get(ph, "gray"),
                           ls="--", lw=0.8, alpha=0.5)
                ax.text(tarr + 2, 0.85, ph, fontsize=7,
                        color=PHASE_COLORS.get(ph, "gray"), fontweight="bold")

        ax.set_xlim(TMIN, TMAX)
        ax.set_ylim(-1.2, 1.2)
        ax.set_title(f"{sta} ({dist_real:.0f} km) — Waveforms",
                     fontsize=10, fontweight="bold")
        ax.legend(fontsize=7, loc="upper right")
        ax.grid(True, alpha=0.2)
        if i == n - 1:
            ax.set_xlabel("Time after origin (s)")

        # Spectral comparison
        ax2 = axes[i][1]
        if len(d_obs) > 0:
            freq_o, sp_o = spectral_amplitude(d_obs, tr_obs.stats.delta)
            m_o = (freq_o > 0.1) & (freq_o < 10)
            if np.any(m_o) and np.max(sp_o[m_o]) > 0:
                ax2.loglog(freq_o[m_o], sp_o[m_o] / np.max(sp_o[m_o]),
                           "k-", lw=1, alpha=0.7, label="Observed")
        if len(d_num) > 0:
            freq_n, sp_n = spectral_amplitude(d_num, tr_num.stats.delta)
            m_n = (freq_n > 0.1) & (freq_n < 10)
            if np.any(m_n) and np.max(sp_n[m_n]) > 0:
                ax2.loglog(freq_n[m_n], sp_n[m_n] / np.max(sp_n[m_n]),
                           "r-", lw=1, alpha=0.7, label="Numerical")

        ax2.set_title(f"{sta} — Spectra", fontsize=10, fontweight="bold")
        ax2.legend(fontsize=7)
        ax2.grid(True, alpha=0.2, which="both")
        ax2.set_xlim(0.1, 10)
        if i == n - 1:
            ax2.set_xlabel("Frequency (Hz)")

    fig.suptitle(f"Per-Station Comparison: FSRM Numerical vs Observed\n"
                 f"{EVENT_LABEL}  |  Q-extrapolated, {FMIN}–{FMAX} Hz",
                 fontsize=13, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    path = os.path.join(OUTDIR, "fig_station_overlay_num_vs_obs.png")
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {path}")


def fig_travel_time_validation(extrap_traces, observed_data, model):
    """Compare P and S arrival times: numerical vs observed vs predicted."""
    if not extrap_traces and not observed_data:
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8))

    # Predicted travel times
    dists_line = np.linspace(100, 1600, 300)
    for phase, vel in PHASE_VELOCITIES.items():
        t_line = dists_line / vel
        col = PHASE_COLORS.get(phase, "gray")
        ax1.plot(dists_line, t_line, color=col, ls="--", lw=2, alpha=0.7,
                 label=f"{phase} ({vel:.1f} km/s)")
        ax2.plot(dists_line, t_line, color=col, ls="--", lw=2, alpha=0.7,
                 label=f"{phase} ({vel:.1f} km/s)")

    # Observed: mark first-break times (estimate from max envelope)
    for net, sta, dist_km, az, tr, desc in observed_data:
        times = time_relative(tr, EVENT_TIME)
        d = tr.data.copy()
        if len(d) == 0:
            continue
        env_d = np.abs(d)
        # Rough P pick: time of first significant energy
        threshold = 0.15 * np.max(env_d)
        above = np.where((env_d > threshold) & (times > 0))[0]
        if len(above) > 0:
            t_pick = times[above[0]]
            ax1.plot(dist_km, t_pick, "ko", ms=4, alpha=0.5, zorder=5)

    ax1.set_xlabel("Distance (km)", fontsize=11)
    ax1.set_ylabel("Travel time (s)", fontsize=11)
    ax1.set_title("Observed P picks vs predicted phases", fontsize=12,
                  fontweight="bold")
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.2)
    ax1.set_xlim(600, 1600)

    # Numerical traces: same analysis
    for tr in extrap_traces:
        dist_km = get_trace_distance(tr)
        t = np.arange(tr.stats.npts) * tr.stats.delta
        d = tr.data.copy().astype(float)
        if len(d) == 0:
            continue
        env_d = np.abs(d)
        threshold = 0.15 * np.max(env_d)
        above = np.where((env_d > threshold) & (t > 0))[0]
        if len(above) > 0:
            t_pick = t[above[0]]
            sta = getattr(tr.stats, "station", "?")
            ax2.plot(dist_km, t_pick, "rs", ms=7, alpha=0.7, zorder=5)

    ax2.set_xlabel("Distance (km)", fontsize=11)
    ax2.set_ylabel("Travel time (s)", fontsize=11)
    ax2.set_title("Numerical P picks (extrapolated) vs predicted phases",
                  fontsize=12, fontweight="bold")
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.2)
    ax2.set_xlim(600, 1600)

    fig.suptitle(f"Travel-Time Validation — {EVENT_LABEL}",
                 fontsize=13, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    path = os.path.join(OUTDIR, "fig_travel_time_validation.png")
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {path}")


# ═══════════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    os.makedirs(OUTDIR, exist_ok=True)
    model = LopNorVelocityModel()
    fc = corner_frequency_patton(YIELD_KT)

    print("=" * 72)
    print("  FSRM Numerical Wave Propagation — Lop Nor 1996-07-29 (~5 kt)")
    print("  Setup, Validation & Comparison with KNDC Observed Data")
    print("=" * 72)
    print()

    # ── Phase 1: Validate config ──
    print("[Phase 1] Validating FSRM configuration ...")
    print()
    validate_config(CONFIG_PATH, **GRID)

    # ── Phase 2: Print invocation ──
    print("[Phase 2] FSRM invocation command:")
    print_invocation(CONFIG_PATH, nprocs=32,
                     nprocs_min=GRID["nprocs_min"],
                     nprocs_max=GRID["nprocs_max"])

    # ── Phase 3: Post-process numerical output ──
    print("[Phase 3] Post-processing FSRM numerical output ...")
    print()
    postprocess_far_field(
        SIM_OUTPUT_DIR, OUTDIR,
        model=model,
        generate_synth_fn=generate_synthetic,
        spectral_fn=spectral_amplitude,
        mm_spectrum_fn=mueller_murphy_spectrum,
        yield_kt=YIELD_KT,
        fc=fc,
        event_label=EVENT_LABEL,
    )

    # ── Phase 4: Load numerical traces & extrapolate to far field ──
    print()
    print("[Phase 4] Far-field extrapolation ...")
    print()

    numerical_traces = read_sac_seismograms(SIM_OUTPUT_DIR)
    if not numerical_traces:
        numerical_traces = read_mseed_seismograms(SIM_OUTPUT_DIR)

    extrap_traces = []
    if numerical_traces:
        # Find azimuthal array traces and extrapolate to real distances
        az_traces = []
        true_dists = []
        domain_dists = []
        for tr in numerical_traces:
            sta = getattr(tr.stats, "station", "")
            if sta in KNDC_STATION_MAP:
                real_d, dom_d, az_deg = KNDC_STATION_MAP[sta]
                az_traces.append(tr)
                true_dists.append(real_d)
                domain_dists.append(dom_d)

        if az_traces:
            print(f"  Found {len(az_traces)} azimuthal station(s) to extrapolate")
            extrap_traces = extrapolate_to_far_field(
                az_traces, true_dists, domain_dists,
                Q_avg=Q_AVG, Vp_avg=VP_AVG, f_ref=1.0)
            print(f"  Extrapolated {len(extrap_traces)} trace(s) to real distances")
        else:
            print("  No azimuthal array traces found in numerical output")
    else:
        print("  No numerical traces available — skipping extrapolation")
        print("  Run the FSRM simulation first, then re-run this script.")

    # ── Phase 5: Compare with KNDC observed data ──
    print()
    print("[Phase 5] Loading observed KNDC archive data ...")
    kndc_data = load_kndc_observed()
    print(f"  Loaded {len(kndc_data)} observed trace(s)")
    print()

    if extrap_traces and kndc_data:
        print("  Generating numerical vs observed comparison figures ...")
        fig_extrapolated_record_section(extrap_traces, kndc_data, model)
        fig_station_overlay(extrap_traces, kndc_data, model)
        fig_travel_time_validation(extrap_traces, kndc_data, model)
    elif kndc_data:
        print("  Generating observed-only figures (no numerical output yet) ...")
        fig_travel_time_validation([], kndc_data, model)
    else:
        print("  No data available for comparison figures")

    print()
    print("=" * 72)
    print(f"  All figures saved to: {OUTDIR}")
    print("=" * 72)


if __name__ == "__main__":
    main()
