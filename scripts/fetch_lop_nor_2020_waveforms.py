#!/usr/bin/env python3
"""
Fetch and display real seismic waveforms from the June 22, 2020 Lop Nor event.

Downloads broadband data from IRIS FDSN for all confirmed stations, removes
instrument response, bandpass-filters, and produces a distance-sorted record
section (waveform gather).

Event:  2020-06-22 ~09:18 UTC, Lop Nor test site (~41.735°N, 88.730°E)
        Detected at PS23 Makanchi (Kazakhstan), mb ≈ 2.75

Usage:
    python scripts/fetch_lop_nor_2020_waveforms.py
"""

import os
import sys
import warnings
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import AutoMinorLocator
from matplotlib.gridspec import GridSpec

from obspy import UTCDateTime, Stream, Trace
from obspy.clients.fdsn import Client
from obspy.geodetics import gps2dist_azimuth
from obspy.signal.trigger import classic_sta_lta

warnings.filterwarnings("ignore", category=UserWarning)

# ──────────────────────────────────────────────────────────────────────────────
# Event parameters
# ──────────────────────────────────────────────────────────────────────────────
EVENT_TIME = UTCDateTime("2020-06-22T09:18:00")
EVENT_LAT = 41.735
EVENT_LON = 88.730
EVENT_DEPTH_KM = 0.3          # shallow, tunnel-depth assumption

# Time windows (seconds relative to origin time)
PRE_ORIGIN = 120              # 2 min before
POST_ORIGIN = 900             # 15 min after  (captures Lg to ~2000 km)

# Processing
FREQMIN = 0.5                 # Hz – lower bound of monitoring band
FREQMAX = 8.0                 # Hz – upper bound of monitoring band
WATER_LEVEL = 60              # dB for instrument-response removal

# Output
OUTDIR = os.path.join(os.path.dirname(os.path.dirname(__file__)),
                      "figures", "lop_nor_2020")

# ──────────────────────────────────────────────────────────────────────────────
# Stations to fetch – all confirmed available on IRIS for this time window
# (net, sta, loc, chan_pref, description)
# ──────────────────────────────────────────────────────────────────────────────
STATIONS = [
    # Closest
    ("IC", "WMQ",  "*", "BH",  "Urumqi, China"),
    # PS23 cluster (main detection)
    ("IU", "MAKZ", "*", "BH",  "Makanchi, Kazakhstan (PS23)"),
    ("KZ", "MKAR", "*", "BH",  "Makanchi Array, Kazakhstan"),
    # Regional
    ("G",  "WUS",  "*", "BH",  "Wushi, China"),
    ("KZ", "PDGK", "*", "BH",  "Podgonoye, Kazakhstan"),
    ("KR", "PRZ",  "*", "BH",  "Karakol, Kyrgyzstan"),
    ("KZ", "KNDC", "*", "BH",  "KNDC Almaty, Kazakhstan"),
    ("II", "AAK",  "*", "BH",  "Ala Archa, Kyrgyzstan"),
    # More distant
    ("II", "KURK", "*", "BH",  "Kurchatov, Kazakhstan"),
    ("IC", "LSA",  "*", "BH",  "Lhasa, Tibet"),
    ("II", "NIL",  "*", "BH",  "Nilore, Pakistan"),
    ("IC", "XAN",  "*", "BH",  "Xi'an, China"),
]


def fetch_waveforms():
    """Download waveforms and station metadata from IRIS."""
    client = Client("IRIS", timeout=60)
    t1 = EVENT_TIME - PRE_ORIGIN
    t2 = EVENT_TIME + POST_ORIGIN

    results = []   # list of (net, sta, distance_km, azimuth, trace_z, description)

    for net, sta, loc, chan_pref, desc in STATIONS:
        chan = f"{chan_pref}Z"
        label = f"{net}.{sta}"
        try:
            st = client.get_waveforms(net, sta, loc, chan, t1, t2)
            if len(st) == 0:
                print(f"  {label:12s}  – empty response, skipping")
                continue

            # Merge if multiple traces
            st.merge(fill_value="interpolate")
            tr = st[0]

            # Get station coordinates
            try:
                inv = client.get_stations(network=net, station=sta,
                                          starttime=EVENT_TIME,
                                          endtime=EVENT_TIME,
                                          level="response",
                                          channel=chan)
            except Exception:
                inv = client.get_stations(network=net, station=sta,
                                          starttime=EVENT_TIME,
                                          endtime=EVENT_TIME,
                                          level="station")
            sta_lat = inv[0][0].latitude
            sta_lon = inv[0][0].longitude
            dist_m, az, baz = gps2dist_azimuth(EVENT_LAT, EVENT_LON,
                                                sta_lat, sta_lon)
            dist_km = dist_m / 1000.0

            # Attach response for removal
            tr.stats.response = inv[0][0][0].response if hasattr(inv[0][0][0], 'response') else None

            results.append((net, sta, dist_km, az, tr, inv, desc))
            print(f"  {label:12s}  dist={dist_km:7.1f} km  az={az:5.1f}°  "
                  f"sr={tr.stats.sampling_rate} Hz  npts={tr.stats.npts}")

        except Exception as exc:
            print(f"  {label:12s}  – {str(exc)[:70]}")

    # Sort by distance
    results.sort(key=lambda x: x[2])
    return results


def process_trace(tr, inv, freqmin, freqmax):
    """Remove instrument response, detrend, taper, bandpass filter."""
    tr = tr.copy()

    # Detrend and taper first
    tr.detrend("demean")
    tr.detrend("linear")
    tr.taper(max_percentage=0.05, type="cosine")

    # Remove instrument response → velocity (m/s)
    try:
        tr.remove_response(inventory=inv, output="VEL",
                           water_level=WATER_LEVEL,
                           pre_filt=[freqmin * 0.5, freqmin,
                                     freqmax, freqmax * 1.25])
    except Exception as exc:
        print(f"    response removal failed for {tr.id}: {exc}")
        # Fall back to simple filter
        tr.detrend("demean")

    # Bandpass
    tr.filter("bandpass", freqmin=freqmin, freqmax=freqmax,
              corners=4, zerophase=True)

    return tr


# ──────────────────────────────────────────────────────────────────────────────
# Plotting helpers
# ──────────────────────────────────────────────────────────────────────────────

def _time_axis(tr, origin):
    """Return time axis in seconds relative to origin."""
    npts = tr.stats.npts
    dt = tr.stats.delta
    start_offset = tr.stats.starttime - origin
    return np.arange(npts) * dt + start_offset


def plot_record_section(results_processed, title_extra="", fname="waveform_gather"):
    """
    Produce a distance-sorted record-section with waveforms normalised per-trace.
    Uses evenly-spaced rows (not true distance) so clustered stations don't overlap.
    Overlays theoretical regional-phase travel-time curves.
    """
    n = len(results_processed)
    if n == 0:
        print("No traces to plot.")
        return

    fig = plt.figure(figsize=(20, max(12, 1.6 * n)))
    gs = GridSpec(1, 24, figure=fig)
    ax_main = fig.add_subplot(gs[0, :20])
    ax_info = fig.add_subplot(gs[0, 20:])
    ax_info.axis("off")

    distances = [r[2] for r in results_processed]

    # Evenly-spaced Y positions so clustered stations don't overlap
    y_positions = np.arange(n)
    trace_height = 0.38   # half-height of each normalised trace in row-units

    # Theoretical travel-time curves (approximate regional velocities)
    phase_velocities = {
        "Pn":  8.1,    # km/s – upper-mantle refraction
        "Pg":  6.1,    # km/s – upper-crustal P
        "Sn":  4.6,    # km/s – upper-mantle S refraction
        "Lg":  3.5,    # km/s – crustal guided wave
    }
    phase_colors = {"Pn": "#1f77b4", "Pg": "#2ca02c",
                    "Sn": "#d62728", "Lg": "#ff7f0e"}
    phase_styles = {"Pn": "--", "Pg": "-.", "Sn": "--", "Lg": "-."}

    # Draw phase move-out lines (interpolate from distance to y-position)
    for phase, vel in phase_velocities.items():
        for i, dist in enumerate(distances):
            t_arr = dist / vel
            ax_main.plot(t_arr, y_positions[i], "d",
                         color=phase_colors[phase], ms=5,
                         alpha=0.5, zorder=3)

    # Plot waveforms
    for i, (net, sta, dist_km, az, tr, inv, desc) in enumerate(results_processed):
        times = _time_axis(tr, EVENT_TIME)
        data = tr.data.copy()

        # Normalise per-trace (peak amplitude = trace_height row-units)
        peak = np.max(np.abs(data))
        if peak > 0:
            data = data / peak * trace_height

        y0 = y_positions[i]

        # Light fill for waveform
        ax_main.fill_between(times, data + y0, y0,
                             where=(data > 0), color="#444", alpha=0.12,
                             zorder=1)
        ax_main.fill_between(times, data + y0, y0,
                             where=(data < 0), color="#444", alpha=0.06,
                             zorder=1)
        ax_main.plot(times, data + y0, "-", lw=0.35,
                     color="k", alpha=0.85, zorder=2)

        # Station label at left – include distance
        ax_main.text(-PRE_ORIGIN - 8, y0,
                     f"{net}.{sta}  ({dist_km:.0f} km)",
                     fontsize=8.5, fontweight="bold", ha="right", va="center",
                     color="#222", clip_on=False)

    # Phase legend via invisible proxy artists
    for phase, vel in phase_velocities.items():
        ax_main.plot([], [], "d", color=phase_colors[phase], ms=6,
                     label=f"{phase}  ({vel} km/s)")

    # Axis formatting
    ax_main.set_xlabel("Time after origin (s)", fontsize=12, fontweight="bold")
    ax_main.set_ylabel("Station  (sorted by distance →)", fontsize=12,
                       fontweight="bold")
    ax_main.set_ylim(-0.8, n - 0.2)
    ax_main.set_xlim(-PRE_ORIGIN, POST_ORIGIN)
    ax_main.set_yticks(y_positions)
    ax_main.set_yticklabels([])          # labels drawn manually above
    ax_main.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax_main.grid(True, which="major", axis="x", ls="-", alpha=0.22)
    ax_main.grid(True, which="minor", axis="x", ls=":", alpha=0.10)
    # Horizontal separator lines
    for y in y_positions:
        ax_main.axhline(y, color="#ccc", lw=0.3, zorder=0)

    ax_main.legend(loc="lower right", fontsize=9, ncol=2,
                   title="Theoretical phase arrivals",
                   title_fontsize=9, framealpha=0.9, edgecolor="#aaa")

    # Info panel
    info_lines = [
        "Event Parameters",
        "─" * 28,
        f"Date:  2020-06-22",
        f"Time:  ~09:18 UTC",
        f"Lat:   {EVENT_LAT:.3f}°N",
        f"Lon:   {EVENT_LON:.3f}°E",
        f"Depth: ~{EVENT_DEPTH_KM} km (est.)",
        f"mb:    ~2.75 (NORSAR)",
        "",
        "Processing",
        "─" * 28,
        f"Filter: {FREQMIN}–{FREQMAX} Hz BP",
        f"Resp:  removed → vel (m/s)",
        f"Norm:  per-trace peak",
        "",
        "Stations (near→far)",
        "─" * 28,
    ]
    for net, sta, dist_km, az, tr, inv, desc in results_processed:
        info_lines.append(f"{net}.{sta:5s} {dist_km:6.0f} km  az {az:5.1f}°")

    ax_info.text(0.05, 0.98, "\n".join(info_lines),
                 transform=ax_info.transAxes, fontfamily="monospace",
                 fontsize=7.5, va="top", ha="left",
                 bbox=dict(boxstyle="round,pad=0.4", fc="#f8f8f8", ec="#bbb"))

    title = (f"Waveform Gather — Lop Nor Event 2020-06-22 09:18 UTC\n"
             f"Vertical component (BHZ), {FREQMIN}–{FREQMAX} Hz bandpass"
             + (f" — {title_extra}" if title_extra else ""))
    fig.suptitle(title, fontsize=14, fontweight="bold", y=0.995)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    outpath = os.path.join(OUTDIR, f"{fname}.png")
    fig.savefig(outpath, dpi=180, bbox_inches="tight")
    print(f"\nSaved: {outpath}")
    plt.close(fig)
    return outpath


def plot_record_section_zoomed(results_processed, fname="waveform_gather_zoomed"):
    """
    Same as record section but zoomed to the expected signal window
    (from 20 s before Pn to 60 s after Lg at each station distance).
    Cuts out later unrelated events.  Time window: -30 s to +350 s.
    """
    n = len(results_processed)
    if n == 0:
        return

    TMIN, TMAX = -30, 400   # focused window

    fig = plt.figure(figsize=(20, max(12, 1.6 * n)))
    gs = GridSpec(1, 24, figure=fig)
    ax_main = fig.add_subplot(gs[0, :20])
    ax_info = fig.add_subplot(gs[0, 20:])
    ax_info.axis("off")

    distances = [r[2] for r in results_processed]
    y_positions = np.arange(n)
    trace_height = 0.38

    phase_velocities = {
        "Pn":  8.1, "Pg": 6.1, "Sn": 4.6, "Lg": 3.5,
    }
    phase_colors = {"Pn": "#1f77b4", "Pg": "#2ca02c",
                    "Sn": "#d62728", "Lg": "#ff7f0e"}

    for phase, vel in phase_velocities.items():
        for i, dist in enumerate(distances):
            t_arr = dist / vel
            if TMIN <= t_arr <= TMAX:
                ax_main.plot(t_arr, y_positions[i], "d",
                             color=phase_colors[phase], ms=6,
                             alpha=0.6, zorder=3)

    for i, (net, sta, dist_km, az, tr, inv, desc) in enumerate(results_processed):
        times = _time_axis(tr, EVENT_TIME)
        data = tr.data.copy()

        # Window
        mask = (times >= TMIN) & (times <= TMAX)
        t_w = times[mask]
        d_w = data[mask]

        peak = np.max(np.abs(d_w)) if len(d_w) > 0 else 1.0
        if peak > 0:
            d_w = d_w / peak * trace_height

        y0 = y_positions[i]
        ax_main.fill_between(t_w, d_w + y0, y0,
                             where=(d_w > 0), color="#444", alpha=0.12, zorder=1)
        ax_main.fill_between(t_w, d_w + y0, y0,
                             where=(d_w < 0), color="#444", alpha=0.06, zorder=1)
        ax_main.plot(t_w, d_w + y0, "-", lw=0.45, color="k", alpha=0.85, zorder=2)

        ax_main.text(TMIN - 4, y0,
                     f"{net}.{sta}  ({dist_km:.0f} km)",
                     fontsize=8.5, fontweight="bold", ha="right", va="center",
                     color="#222", clip_on=False)

    for phase, vel in phase_velocities.items():
        ax_main.plot([], [], "d", color=phase_colors[phase], ms=6,
                     label=f"{phase}  ({vel} km/s)")

    ax_main.set_xlabel("Time after origin (s)", fontsize=12, fontweight="bold")
    ax_main.set_ylabel("Station  (sorted by distance →)", fontsize=12,
                       fontweight="bold")
    ax_main.set_ylim(-0.8, n - 0.2)
    ax_main.set_xlim(TMIN, TMAX)
    ax_main.set_yticks(y_positions)
    ax_main.set_yticklabels([])
    ax_main.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax_main.grid(True, which="major", axis="x", ls="-", alpha=0.22)
    ax_main.grid(True, which="minor", axis="x", ls=":", alpha=0.10)
    for y in y_positions:
        ax_main.axhline(y, color="#ccc", lw=0.3, zorder=0)

    ax_main.legend(loc="lower right", fontsize=9, ncol=2,
                   title="Theoretical phase arrivals",
                   title_fontsize=9, framealpha=0.9, edgecolor="#aaa")

    info_lines = [
        "Zoomed View",
        "─" * 28,
        f"Window:  {TMIN} – {TMAX} s",
        f"Isolates signal from",
        f"later unrelated events",
        "",
        "Event: 2020-06-22 09:18 UTC",
        f"mb ~2.75 (NORSAR est.)",
        "",
        "Phase velocities:",
        f"  Pn  = 8.1 km/s",
        f"  Pg  = 6.1 km/s",
        f"  Sn  = 4.6 km/s",
        f"  Lg  = 3.5 km/s",
    ]
    ax_info.text(0.05, 0.98, "\n".join(info_lines),
                 transform=ax_info.transAxes, fontfamily="monospace",
                 fontsize=8, va="top", ha="left",
                 bbox=dict(boxstyle="round,pad=0.4", fc="#f8f8f8", ec="#bbb"))

    title = (f"Waveform Gather (zoomed) — Lop Nor Event 2020-06-22 09:18 UTC\n"
             f"BHZ, {FREQMIN}–{FREQMAX} Hz — expected signal window only")
    fig.suptitle(title, fontsize=14, fontweight="bold", y=0.995)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    outpath = os.path.join(OUTDIR, f"{fname}.png")
    fig.savefig(outpath, dpi=180, bbox_inches="tight")
    print(f"\nSaved: {outpath}")
    plt.close(fig)
    return outpath


def plot_spectrograms(results_processed, fname="spectrograms"):
    """Plot spectrograms for each station, 2 columns."""
    n = len(results_processed)
    if n == 0:
        return
    ncols = 2
    nrows = (n + 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(18, 3.0 * nrows),
                             squeeze=False)

    for i, (net, sta, dist_km, az, tr, inv, desc) in enumerate(results_processed):
        row, col = divmod(i, ncols)
        ax = axes[row][col]

        times = _time_axis(tr, EVENT_TIME)
        data = tr.data.astype(np.float64)

        # Spectrogram via matplotlib
        sr = tr.stats.sampling_rate
        nfft = int(2 ** np.ceil(np.log2(sr * 4)))   # ~4-s windows
        noverlap = nfft // 2

        # Shift time so t=0 is origin
        ax.specgram(data, NFFT=nfft, Fs=sr, noverlap=noverlap,
                    xextent=(times[0], times[-1]),
                    cmap="inferno", vmin=-200, vmax=-120)

        # Theoretical arrivals
        for phase, vel in {"Pn": 8.1, "Pg": 6.1, "Sn": 4.6, "Lg": 3.5}.items():
            t_arr = dist_km / vel
            ax.axvline(t_arr, color="w", ls="--", lw=0.9, alpha=0.8)
            ax.text(t_arr + 1, FREQMAX * 0.92, phase, color="w",
                    fontsize=7, fontweight="bold", va="top")

        ax.set_ylim(FREQMIN, FREQMAX)
        ax.set_xlim(-PRE_ORIGIN, POST_ORIGIN)
        ax.set_title(f"{net}.{sta}  ({dist_km:.0f} km, az {az:.0f}°)",
                     fontsize=9, fontweight="bold")
        ax.set_ylabel("Freq (Hz)", fontsize=8)
        ax.set_xlabel("Time after origin (s)", fontsize=8)

    # Hide unused axes
    for j in range(n, nrows * ncols):
        row, col = divmod(j, ncols)
        axes[row][col].axis("off")

    fig.suptitle("Spectrograms — Lop Nor Event 2020-06-22 09:18 UTC\n"
                 f"Vertical component, {FREQMIN}–{FREQMAX} Hz",
                 fontsize=12, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    outpath = os.path.join(OUTDIR, f"{fname}.png")
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    print(f"Saved: {outpath}")
    plt.close(fig)
    return outpath


def plot_closeup_ps23(results_processed, fname="closeup_ps23"):
    """
    Detailed multi-panel view of PS23 / Makanchi — the primary detection station.
    Top row:  full 15-min window with envelope
    Bottom row: zoomed on Pn arrival window (60–160 s)
    Also includes WMQ (closest) if available.
    """
    from scipy.signal import hilbert

    # Stations of primary interest for close-up
    target_stas = ["MAKZ", "MKAR", "WMQ"]
    entries = [r for r in results_processed if r[1] in target_stas]
    if not entries:
        print("No PS23/WMQ data available for close-up.")
        return

    n = len(entries)
    fig, axes = plt.subplots(n, 2, figsize=(20, 4.2 * n), squeeze=False)

    for idx, (net, sta, dist_km, az, tr, inv, desc) in enumerate(entries):
        times = _time_axis(tr, EVENT_TIME)
        data = tr.data.copy()

        peak = np.max(np.abs(data))
        if peak == 0:
            peak = 1.0
        norm_data = data / peak

        analytic = hilbert(data)
        envelope = np.abs(analytic) / peak

        # Expected phase arrivals
        phases = [("Pn", 8.1, "#1f77b4"), ("Pg", 6.1, "#2ca02c"),
                  ("Sn", 4.6, "#d62728"), ("Lg", 3.5, "#ff7f0e")]

        # ── Left panel: full window ──
        ax = axes[idx][0]
        ax.plot(times, norm_data, "k-", lw=0.35, alpha=0.7, label="Velocity")
        ax.plot(times, envelope, "r-", lw=0.8, alpha=0.5, label="Envelope")
        ax.plot(times, -envelope, "r-", lw=0.8, alpha=0.5)

        for phase, vel, col in phases:
            t_arr = dist_km / vel
            ax.axvline(t_arr, color=col, ls="--", lw=1.2, alpha=0.7)
            ax.text(t_arr + 3, 0.92, phase, color=col, fontsize=9,
                    fontweight="bold", va="top")

        ax.set_xlim(-30, POST_ORIGIN)
        ax.set_ylim(-1.15, 1.15)
        ax.set_xlabel("Time after origin (s)", fontsize=10)
        ax.set_ylabel("Normalised velocity", fontsize=10)
        ax.set_title(f"{net}.{sta} — {desc} ({dist_km:.0f} km)  |  Full window",
                     fontsize=10, fontweight="bold")
        ax.legend(loc="upper right", fontsize=8)
        ax.grid(True, alpha=0.2)
        ax.xaxis.set_minor_locator(AutoMinorLocator(5))

        # ── Right panel: zoom on P-wave arrival window ──
        ax2 = axes[idx][1]
        # Zoom to Pn-30 .. Lg+60
        t_pn = dist_km / 8.1
        t_lg = dist_km / 3.5
        t_zoom_start = max(0, t_pn - 30)
        t_zoom_end = min(POST_ORIGIN, t_lg + 80)

        mask = (times >= t_zoom_start) & (times <= t_zoom_end)
        t_w = times[mask]
        d_w = norm_data[mask] if np.any(mask) else norm_data[:10]
        e_w = envelope[mask] if np.any(mask) else envelope[:10]

        # Re-normalise in zoom window for better visibility
        zoom_peak = np.max(np.abs(d_w)) if len(d_w) > 0 else 1.0
        if zoom_peak > 0:
            d_w = d_w / zoom_peak
            e_w = e_w / zoom_peak

        ax2.plot(t_w, d_w, "k-", lw=0.5, alpha=0.8, label="Velocity")
        ax2.plot(t_w, e_w, "r-", lw=1.0, alpha=0.55, label="Envelope")
        ax2.plot(t_w, -e_w, "r-", lw=1.0, alpha=0.55)

        for phase, vel, col in phases:
            t_arr = dist_km / vel
            if t_zoom_start <= t_arr <= t_zoom_end:
                ax2.axvline(t_arr, color=col, ls="--", lw=1.3, alpha=0.8)
                ax2.text(t_arr + 1, 0.92, phase, color=col, fontsize=10,
                         fontweight="bold", va="top")

        ax2.set_xlim(t_zoom_start, t_zoom_end)
        ax2.set_ylim(-1.15, 1.15)
        ax2.set_xlabel("Time after origin (s)", fontsize=10)
        ax2.set_ylabel("Normalised (zoom window)", fontsize=10)
        ax2.set_title(f"{net}.{sta}  |  Zoomed: Pn–Lg arrival window",
                      fontsize=10, fontweight="bold")
        ax2.legend(loc="upper right", fontsize=8)
        ax2.grid(True, alpha=0.2)
        ax2.xaxis.set_minor_locator(AutoMinorLocator(5))

    fig.suptitle(
        "Close-up: Key Station Detections — Lop Nor 2020-06-22 09:18 UTC\n"
        f"Bandpass {FREQMIN}–{FREQMAX} Hz, instrument response removed → velocity",
        fontsize=13, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    outpath = os.path.join(OUTDIR, f"{fname}.png")
    fig.savefig(outpath, dpi=180, bbox_inches="tight")
    print(f"Saved: {outpath}")
    plt.close(fig)
    return outpath


# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────

def main():
    os.makedirs(OUTDIR, exist_ok=True)

    print("=" * 72)
    print("  Lop Nor 2020-06-22 Event — Waveform Gather from FDSN/IRIS")
    print("=" * 72)
    print(f"  Origin time:  {EVENT_TIME}")
    print(f"  Location:     {EVENT_LAT:.3f}°N  {EVENT_LON:.3f}°E")
    print(f"  Window:       {PRE_ORIGIN} s before → {POST_ORIGIN} s after origin")
    print(f"  Filter:       {FREQMIN}–{FREQMAX} Hz bandpass")
    print()

    # 1. Download
    print("Downloading waveforms ...")
    raw_results = fetch_waveforms()
    print(f"\n  Retrieved {len(raw_results)} station(s)\n")

    if len(raw_results) == 0:
        print("ERROR: No waveforms retrieved.  Exiting.")
        sys.exit(1)

    # 2. Process
    print("Processing (instrument response removal + bandpass filter) ...")
    processed = []
    for net, sta, dist_km, az, tr, inv, desc in raw_results:
        try:
            tr_proc = process_trace(tr, inv, FREQMIN, FREQMAX)
            processed.append((net, sta, dist_km, az, tr_proc, inv, desc))
            print(f"  ✓ {net}.{sta}")
        except Exception as exc:
            print(f"  ✗ {net}.{sta}: {exc}")

    print(f"\n  Processed {len(processed)} trace(s)\n")

    # 3. Plot record section (full window)
    print("Generating waveform gather (record section – full window) ...")
    plot_record_section(processed,
                        title_extra="instrument-response-corrected velocity",
                        fname="observed_waveform_gather")

    # 4. Plot record section (zoomed on expected signal)
    print("Generating waveform gather (record section – zoomed) ...")
    plot_record_section_zoomed(processed,
                               fname="observed_waveform_gather_zoomed")

    # 5. Plot spectrograms
    print("Generating spectrograms ...")
    plot_spectrograms(processed, fname="observed_spectrograms")

    # 6. Close-up on PS23 / WMQ
    print("Generating PS23/WMQ close-up ...")
    plot_closeup_ps23(processed, fname="observed_closeup_ps23")

    print()
    print("=" * 72)
    print("  Done.  Figures saved to:", OUTDIR)
    print("=" * 72)


if __name__ == "__main__":
    main()
