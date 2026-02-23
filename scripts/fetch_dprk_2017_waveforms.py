#!/usr/bin/env python3
"""
Fetch and display real seismic waveforms from the 2017 DPRK nuclear test.

Downloads broadband data from IRIS FDSN for regional stations, removes
instrument response, bandpass-filters, and produces a distance-sorted record
section (waveform gather).

Event:  2017-09-03 03:30:01 UTC, DPRK nuclear test site (~41.300°N, 129.076°E)
        mb ~6.3, depth ~760 m, fully coupled ~250 kt

Usage:
    python scripts/fetch_dprk_2017_waveforms.py
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
EVENT_TIME = UTCDateTime("2017-09-03T03:30:01")
EVENT_LAT = 41.300
EVENT_LON = 129.076
EVENT_DEPTH_KM = 0.76            # shallow, ~760 m depth

# Time windows (seconds relative to origin time)
PRE_ORIGIN = 60                  # 1 min before (well-detected event)
POST_ORIGIN = 600                # 10 min after

# Processing
FREQMIN = 0.5                    # Hz – lower bound of monitoring band
FREQMAX = 8.0                    # Hz – upper bound of monitoring band
WATER_LEVEL = 60                 # dB for instrument-response removal

# Output
OUTDIR = os.path.join(os.path.dirname(os.path.dirname(__file__)),
                      "figures", "dprk_2017")

# ──────────────────────────────────────────────────────────────────────────────
# Stations to fetch – regional broadband stations for this event
# (net, sta, loc, chan_pref, description)
# ──────────────────────────────────────────────────────────────────────────────
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
        f"Date:  2017-09-03",
        f"Time:  03:30:01 UTC",
        f"Lat:   {EVENT_LAT:.3f}°N",
        f"Lon:   {EVENT_LON:.3f}°E",
        f"Depth: ~{EVENT_DEPTH_KM} km (~760 m)",
        f"mb:    ~6.3",
        f"Yield: ~250 kt (fully coupled)",
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

    title = (f"Waveform Gather — DPRK 2017-09-03 Nuclear Test\n"
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
    Same as record section but zoomed to the expected signal window.
    Time window: -30 s to +300 s.
    """
    n = len(results_processed)
    if n == 0:
        return

    TMIN, TMAX = -30, 300   # focused window

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
        "Event: 2017-09-03 03:30:01 UTC",
        f"mb ~6.3, ~250 kt",
        f"Depth ~760 m",
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

    title = (f"Waveform Gather (zoomed) — DPRK 2017-09-03 Nuclear Test\n"
             f"BHZ, {FREQMIN}–{FREQMAX} Hz — expected signal window only")
    fig.suptitle(title, fontsize=14, fontweight="bold", y=0.995)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    outpath = os.path.join(OUTDIR, f"{fname}.png")
    fig.savefig(outpath, dpi=180, bbox_inches="tight")
    print(f"\nSaved: {outpath}")
    plt.close(fig)
    return outpath


def plot_record_section_distance_scaled(results_processed,
                                        fname="waveform_gather_distance_scaled"):
    """
    Record section with Y-axis scaled to true epicentral distance.

    This allows theoretical phase arrivals to be shown as smooth straight lines
    through the origin (t = d/v), making it easy to identify phases and measure
    apparent velocities.

    Waveform amplitudes are scaled to fit the vertical spacing between stations.
    """
    from scipy.interpolate import interp1d

    n = len(results_processed)
    if n == 0:
        print("No traces to plot.")
        return

    fig = plt.figure(figsize=(18, 14))
    gs = GridSpec(1, 20, figure=fig, wspace=0.05)
    ax_main = fig.add_subplot(gs[0, :17])
    ax_info = fig.add_subplot(gs[0, 17:])
    ax_info.axis("off")

    # Get distance range
    distances = np.array([r[2] for r in results_processed])
    dist_min = np.min(distances)
    dist_max = np.max(distances)

    # Add padding to distance range
    dist_pad = (dist_max - dist_min) * 0.08
    y_min = max(0, dist_min - dist_pad)
    y_max = dist_max + dist_pad

    # Time window settings
    TMIN, TMAX = -20, 500  # seconds relative to origin

    # Theoretical phase velocities (km/s) and colors
    phases = {
        "Pn": {"vel": 8.1, "color": "#1f77b4", "ls": "-", "lw": 1.2},
        "Pg": {"vel": 6.1, "color": "#2ca02c", "ls": "--", "lw": 1.0},
        "Sn": {"vel": 4.6, "color": "#d62728", "ls": "-", "lw": 1.2},
        "Lg": {"vel": 3.5, "color": "#ff7f0e", "ls": "--", "lw": 1.0},
    }

    # Calculate waveform amplitude scaling
    distances_sorted = np.sort(distances)
    if len(distances_sorted) > 1:
        spacings = np.diff(distances_sorted)
        significant_spacings = spacings[spacings > 20]
        if len(significant_spacings) > 0:
            min_spacing = np.min(significant_spacings)
        else:
            min_spacing = (dist_max - dist_min) / n
    else:
        min_spacing = 100

    waveform_height = max(40, min_spacing * 0.40)

    # ── Draw smooth theoretical phase arrival lines ──
    d_line = np.linspace(0, y_max * 1.2, 500)

    for phase_name, props in phases.items():
        vel = props["vel"]
        t_line = d_line / vel

        mask = (t_line >= TMIN) & (t_line <= TMAX)
        if np.any(mask):
            ax_main.plot(t_line[mask], d_line[mask],
                        color=props["color"],
                        ls=props["ls"],
                        lw=props["lw"],
                        alpha=0.85,
                        label=f"{phase_name} ({vel:.1f} km/s)",
                        zorder=10)

    # ── Plot waveforms at their true distances ──
    for i, (net, sta, dist_km, az, tr, inv, desc) in enumerate(results_processed):
        times = _time_axis(tr, EVENT_TIME)
        data = tr.data.copy()

        # Window to time range
        mask = (times >= TMIN) & (times <= TMAX)
        t_plot = times[mask]
        d_plot = data[mask]

        # Normalize and scale
        peak = np.max(np.abs(d_plot))
        if peak > 0:
            d_plot = d_plot / peak * waveform_height

        # Plot waveform centered at this station's distance
        y_waveform = d_plot + dist_km

        # Fill positive and negative wiggles
        ax_main.fill_between(t_plot, y_waveform, dist_km,
                            where=(d_plot > 0),
                            color="black", alpha=0.35, zorder=4)
        ax_main.fill_between(t_plot, y_waveform, dist_km,
                            where=(d_plot < 0),
                            color="gray", alpha=0.20, zorder=4)

        # Plot trace line
        ax_main.plot(t_plot, y_waveform, 'k-', lw=0.6, alpha=0.9, zorder=5)

        # Baseline at station distance
        ax_main.axhline(dist_km, color='#bbbbbb', lw=0.3, zorder=1)

        # Station label on left edge
        ax_main.text(TMIN - 5, dist_km, f"{net}.{sta}",
                    fontsize=8, fontweight="bold",
                    ha="right", va="center", color="#333",
                    clip_on=False)

    # ── Mark phase arrivals on each waveform ──
    for i, (net, sta, dist_km, az, tr, inv, desc) in enumerate(results_processed):
        for phase_name, props in phases.items():
            t_arr = dist_km / props["vel"]
            if TMIN <= t_arr <= TMAX:
                ax_main.plot(t_arr, dist_km, 'o',
                            color=props["color"],
                            markersize=4,
                            markeredgecolor='white',
                            markeredgewidth=0.5,
                            zorder=11)

    # ── Axis formatting ──
    ax_main.set_xlabel("Time after origin (s)", fontsize=12, fontweight="bold")
    ax_main.set_ylabel("Epicentral Distance (km)", fontsize=12, fontweight="bold")
    ax_main.set_xlim(TMIN, TMAX)
    ax_main.set_ylim(y_min, y_max)

    ax_main.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax_main.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax_main.grid(True, which="major", ls="-", alpha=0.15)
    ax_main.grid(True, which="minor", ls=":", alpha=0.08)

    # Legend for phases
    ax_main.legend(loc="upper left", fontsize=10,
                  title="Theoretical Phase Arrivals", title_fontsize=10,
                  framealpha=0.95, edgecolor="#aaa", fancybox=True)

    # ── Info panel ──
    info_lines = [
        "Event Information",
        "═" * 26,
        "",
        f"Date:   2017-09-03",
        f"Time:   03:30:01 UTC",
        f"Lat:    {EVENT_LAT:.3f}°N",
        f"Lon:    {EVENT_LON:.3f}°E",
        f"Depth:  ~{EVENT_DEPTH_KM} km (~760 m)",
        f"mb:     ~6.3",
        f"Yield:  ~250 kt (fully coupled)",
        "",
        "Processing",
        "─" * 26,
        f"Filter: {FREQMIN}–{FREQMAX} Hz",
        f"Response: removed",
        f"Output: velocity",
        "",
        "Phase Velocities",
        "─" * 26,
        f"Pn (mantle P): 8.1 km/s",
        f"Pg (crustal P): 6.1 km/s",
        f"Sn (mantle S): 4.6 km/s",
        f"Lg (crustal S): 3.5 km/s",
        "",
        "Stations",
        "─" * 26,
    ]
    for net, sta, dist_km, az, tr, inv, desc in results_processed:
        info_lines.append(f"{net}.{sta:5s} {dist_km:6.0f} km")

    ax_info.text(0.08, 0.98, "\n".join(info_lines),
                transform=ax_info.transAxes, fontfamily="monospace",
                fontsize=8, va="top", ha="left",
                bbox=dict(boxstyle="round,pad=0.5", fc="#f9f9f4", ec="#bbb"))

    # ── Title ──
    title = ("Waveform Gather — DPRK 2017-09-03 Nuclear Test\n"
             f"Distance-scaled record section with theoretical phase arrivals (BHZ, {FREQMIN}–{FREQMAX} Hz)")
    fig.suptitle(title, fontsize=13, fontweight="bold", y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    outpath = os.path.join(OUTDIR, f"{fname}.png")
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    print(f"\nSaved: {outpath}")
    plt.close(fig)
    return outpath


def plot_record_section_long_period(results_raw,
                                    fname="waveform_gather_long_period"):
    """
    Record section with long-period filtering for moment tensor inversion.

    Uses 0.02-0.10 Hz bandpass (10-50 s period) - the standard band for
    regional moment tensor inversion of surface waves.

    Y-axis is true epicentral distance with smooth phase arrival lines.
    """
    from scipy.signal import butter, filtfilt

    # Long-period filter band for MT inversion
    LP_FREQMIN = 0.02  # Hz (50 s period)
    LP_FREQMAX = 0.10  # Hz (10 s period)

    n = len(results_raw)
    if n == 0:
        print("No traces to plot.")
        return

    # Process traces with long-period filter
    results_processed = []
    print(f"  Applying long-period filter ({LP_FREQMIN}–{LP_FREQMAX} Hz)...")
    for net, sta, dist_km, az, tr_orig, inv, desc in results_raw:
        tr = tr_orig.copy()
        try:
            # Remove instrument response first
            tr.detrend('demean')
            tr.detrend('linear')
            tr.taper(0.05)
            try:
                tr.remove_response(inventory=inv, output='VEL',
                                  water_level=60,
                                  pre_filt=[LP_FREQMIN * 0.5, LP_FREQMIN,
                                           LP_FREQMAX, LP_FREQMAX * 1.5])
            except Exception:
                pass  # Continue without response removal if it fails

            # Apply long-period bandpass filter
            tr.filter('bandpass', freqmin=LP_FREQMIN, freqmax=LP_FREQMAX,
                     corners=3, zerophase=True)
            results_processed.append((net, sta, dist_km, az, tr, inv, desc))
        except Exception as e:
            print(f"    {net}.{sta}: filter failed - {e}")

    if len(results_processed) == 0:
        print("No traces after filtering.")
        return

    fig = plt.figure(figsize=(18, 14))
    gs = GridSpec(1, 20, figure=fig, wspace=0.05)
    ax_main = fig.add_subplot(gs[0, :17])
    ax_info = fig.add_subplot(gs[0, 17:])
    ax_info.axis("off")

    # Get distance range
    distances = np.array([r[2] for r in results_processed])
    dist_min = np.min(distances)
    dist_max = np.max(distances)

    # Add padding to distance range
    dist_pad = (dist_max - dist_min) * 0.08
    y_min = max(0, dist_min - dist_pad)
    y_max = dist_max + dist_pad

    # Time window - longer for surface waves
    TMIN, TMAX = -50, 800  # seconds relative to origin

    # Theoretical phase velocities for long-period waves
    phases = {
        "Pn": {"vel": 8.1, "color": "#1f77b4", "ls": "-", "lw": 1.2},
        "Sn": {"vel": 4.6, "color": "#d62728", "ls": "-", "lw": 1.2},
        "Rayleigh": {"vel": 3.2, "color": "#9467bd", "ls": "-", "lw": 1.2},
        "Love": {"vel": 3.8, "color": "#8c564b", "ls": "--", "lw": 1.0},
    }

    # Calculate waveform amplitude scaling
    distances_sorted = np.sort(distances)
    if len(distances_sorted) > 1:
        spacings = np.diff(distances_sorted)
        significant_spacings = spacings[spacings > 20]
        if len(significant_spacings) > 0:
            min_spacing = np.min(significant_spacings)
        else:
            min_spacing = (dist_max - dist_min) / n
    else:
        min_spacing = 100

    waveform_height = max(50, min_spacing * 0.45)

    # ── Plot waveforms at their true distances ──
    for i, (net, sta, dist_km, az, tr, inv, desc) in enumerate(results_processed):
        times = _time_axis(tr, EVENT_TIME)
        data = tr.data.copy()

        # Window to time range
        mask = (times >= TMIN) & (times <= TMAX)
        t_plot = times[mask]
        d_plot = data[mask]

        # Normalize and scale
        peak = np.max(np.abs(d_plot))
        if peak > 0:
            d_plot = d_plot / peak * waveform_height

        # Plot waveform centered at this station's distance
        y_waveform = d_plot + dist_km

        # Fill positive and negative wiggles
        ax_main.fill_between(t_plot, y_waveform, dist_km,
                            where=(d_plot > 0),
                            color="black", alpha=0.35, zorder=4)
        ax_main.fill_between(t_plot, y_waveform, dist_km,
                            where=(d_plot < 0),
                            color="gray", alpha=0.20, zorder=4)

        # Plot trace line
        ax_main.plot(t_plot, y_waveform, 'k-', lw=0.6, alpha=0.9, zorder=5)

        # Baseline at station distance
        ax_main.axhline(dist_km, color='#bbbbbb', lw=0.3, zorder=1)

        # Station label on left edge
        ax_main.text(TMIN - 10, dist_km, f"{net}.{sta}",
                    fontsize=8, fontweight="bold",
                    ha="right", va="center", color="#333",
                    clip_on=False)

    # ── Draw smooth theoretical phase arrival lines ──
    d_line = np.linspace(0, y_max * 1.2, 500)

    for phase_name, props in phases.items():
        vel = props["vel"]
        t_line = d_line / vel

        mask = (t_line >= TMIN) & (t_line <= TMAX)
        if np.any(mask):
            ax_main.plot(t_line[mask], d_line[mask],
                        color=props["color"],
                        ls=props["ls"],
                        lw=props["lw"],
                        alpha=0.85,
                        label=f"{phase_name} ({vel:.1f} km/s)",
                        zorder=10)

    # ── Mark phase arrivals on each waveform ──
    for i, (net, sta, dist_km, az, tr, inv, desc) in enumerate(results_processed):
        for phase_name, props in phases.items():
            t_arr = dist_km / props["vel"]
            if TMIN <= t_arr <= TMAX:
                ax_main.plot(t_arr, dist_km, 'o',
                            color=props["color"],
                            markersize=4,
                            markeredgecolor='white',
                            markeredgewidth=0.5,
                            zorder=11)

    # ── Axis formatting ──
    ax_main.set_xlabel("Time after origin (s)", fontsize=12, fontweight="bold")
    ax_main.set_ylabel("Epicentral Distance (km)", fontsize=12, fontweight="bold")
    ax_main.set_xlim(TMIN, TMAX)
    ax_main.set_ylim(y_min, y_max)

    ax_main.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax_main.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax_main.grid(True, which="major", ls="-", alpha=0.15)
    ax_main.grid(True, which="minor", ls=":", alpha=0.08)

    # Legend for phases
    ax_main.legend(loc="upper left", fontsize=10,
                  title="Theoretical Phase Arrivals", title_fontsize=10,
                  framealpha=0.95, edgecolor="#aaa", fancybox=True)

    # ── Info panel ──
    info_lines = [
        "Long-Period Filter",
        "═" * 26,
        "",
        f"Bandpass: {LP_FREQMIN}–{LP_FREQMAX} Hz",
        f"Period:   {1/LP_FREQMAX:.0f}–{1/LP_FREQMIN:.0f} s",
        "",
        "Used for moment tensor",
        "inversion (surface waves)",
        "",
        "Event Information",
        "─" * 26,
        f"Date:   2017-09-03",
        f"Time:   03:30:01 UTC",
        f"Lat:    {EVENT_LAT:.3f}°N",
        f"Lon:    {EVENT_LON:.3f}°E",
        f"mb:     ~6.3",
        f"Yield:  ~250 kt (fully coupled)",
        "",
        "Surface Wave Velocities",
        "─" * 26,
        f"Rayleigh: 3.2 km/s",
        f"Love:     3.8 km/s",
        "",
        "Stations",
        "─" * 26,
    ]
    for net, sta, dist_km, az, tr, inv, desc in results_processed:
        info_lines.append(f"{net}.{sta:5s} {dist_km:6.0f} km")

    ax_info.text(0.08, 0.98, "\n".join(info_lines),
                transform=ax_info.transAxes, fontfamily="monospace",
                fontsize=8, va="top", ha="left",
                bbox=dict(boxstyle="round,pad=0.5", fc="#f9f9f4", ec="#bbb"))

    # ── Title ──
    title = ("Waveform Gather — DPRK 2017-09-03 Nuclear Test\n"
             f"Long-period filter ({LP_FREQMIN}–{LP_FREQMAX} Hz) for moment tensor inversion")
    fig.suptitle(title, fontsize=13, fontweight="bold", y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    outpath = os.path.join(OUTDIR, f"{fname}.png")
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
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

    fig.suptitle("Spectrograms — DPRK 2017-09-03 Nuclear Test\n"
                 f"Vertical component, {FREQMIN}–{FREQMAX} Hz",
                 fontsize=12, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    outpath = os.path.join(OUTDIR, f"{fname}.png")
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    print(f"Saved: {outpath}")
    plt.close(fig)
    return outpath


def plot_closeup_ps23(results_processed, fname="closeup_key_stations"):
    """
    Detailed multi-panel view of closest/key stations — MDJ, INCN, MAJO.
    Top row:  full window with envelope
    Bottom row: zoomed on Pn arrival window
    """
    from scipy.signal import hilbert

    target_stas = ["MDJ", "INCN", "MAJO"]
    entries = [r for r in results_processed if r[1] in target_stas]
    if not entries:
        print("No key station data available for close-up.")
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
        "Close-up: Key Station Detections — DPRK 2017-09-03 Nuclear Test\n"
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
    print("  DPRK 2017-09-03 Nuclear Test — Waveform Gather from FDSN/IRIS")
    print("=" * 72)
    print(f"  Origin time:  {EVENT_TIME}")
    print(f"  Location:     {EVENT_LAT:.3f}°N  {EVENT_LON:.3f}°E")
    print(f"  Depth:        ~{EVENT_DEPTH_KM} km (~760 m)")
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

    # 5. Plot distance-scaled record section with smooth phase lines
    print("Generating distance-scaled record section with phase arrivals ...")
    plot_record_section_distance_scaled(processed,
                                        fname="observed_waveform_gather_distance_scaled")

    # 6. Plot long-period filtered record section (for MT inversion band)
    print("Generating long-period filtered record section (0.02–0.10 Hz) ...")
    plot_record_section_long_period(raw_results,
                                    fname="observed_waveform_gather_long_period")

    # 7. Plot spectrograms
    print("Generating spectrograms ...")
    plot_spectrograms(processed, fname="observed_spectrograms")

    # 8. Close-up on key stations (MDJ, INCN, MAJO)
    print("Generating key station close-up ...")
    plot_closeup_ps23(processed, fname="observed_closeup_key_stations")

    print()
    print("=" * 72)
    print("  Done.  Figures saved to:", OUTDIR)
    print("=" * 72)


if __name__ == "__main__":
    main()
