#!/usr/bin/env python3
"""
Plot all seismic data for the alleged 2020-06-22 Lop Nor nuclear test.

Fetches broadband data from IRIS FDSN for all confirmed regional stations
plus the PS23 Makanchi short-period array, and produces:

  1. PyGMT station map (event + all stations including PS23 array elements)
  2. Multi-filter record sections (1-10 Hz, 0.2-8 Hz, 5 Hz HP)
  3. Two-pulse timing display (cut to show ~12 s separation at close stations)
  4. PS23 array waveforms in three filter bands
  5. Array stack at minimum lag (cross-correlation optimised delay-and-sum)

Event:
    2020-06-22 ~09:18 UTC
    Location: 41.735 N, 88.730 E (Lop Nor test site, Kuruktag mountains)
    Detection: PS23 Makanchi array, Kazakhstan (~780 km), mb ~ 2.75
    CTBTO note: Two very small events separated by ~12 seconds

Usage:
    python scripts/plot_lop_nor_2020_all_data.py
"""

import os
import sys
import tempfile
import warnings

import numpy as np
from scipy.signal import hilbert

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import pygmt
import pandas as pd

from obspy import UTCDateTime, Stream
from obspy.clients.fdsn import Client as FDSNClient
from obspy.geodetics import gps2dist_azimuth

warnings.filterwarnings("ignore")

# ═══════════════════════════════════════════════════════════════════════════════
# Event parameters
# ═══════════════════════════════════════════════════════════════════════════════
EVENT_TIME = UTCDateTime("2020-06-22T09:18:00")
EVENT_LAT = 41.735
EVENT_LON = 88.730
EVENT_DEPTH_KM = 0.3
MB_OBSERVED = 2.75

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)
OUTDIR = os.path.join(PROJECT_DIR, "figures", "lop_nor_2020_all_data")

# Filter definitions: (label, kwargs for obspy filter)
FILTERS = [
    ("1\u201310 Hz BP",  dict(type="bandpass", freqmin=1.0, freqmax=10.0,
                              corners=4, zerophase=True)),
    ("0.2\u20138 Hz BP", dict(type="bandpass", freqmin=0.2, freqmax=8.0,
                              corners=4, zerophase=True)),
    ("5 Hz HP",          dict(type="highpass", freq=5.0,
                              corners=4, zerophase=True)),
]

# Phase velocities for move-out lines (km/s)
PHASE_VELOCITIES = {"Pn": 8.1, "Pg": 6.1, "Sn": 4.6, "Lg": 3.5}
PHASE_COLORS = {"Pn": "#1f77b4", "Pg": "#2ca02c",
                "Sn": "#d62728", "Lg": "#ff7f0e"}

# ═══════════════════════════════════════════════════════════════════════════════
# Broadband stations
# ═══════════════════════════════════════════════════════════════════════════════
BB_STATIONS = [
    # (net, sta, loc, chan_pref, description)
    ("IC",  "WMQ",  "*", "BH", "Urumqi, China"),
    ("IU",  "MAKZ", "*", "BH", "Makanchi, Kazakhstan (PS23)"),
    ("KZ",  "MKAR", "*", "BH", "Makanchi Array, Kazakhstan"),
    ("G",   "WUS",  "*", "BH", "Wushi, China"),
    ("KZ",  "PDGK", "*", "BH", "Podgonoye, Kazakhstan"),
    ("KR",  "PRZ",  "*", "BH", "Karakol, Kyrgyzstan"),
    ("KZ",  "KNDC", "*", "BH", "KNDC Almaty, Kazakhstan"),
    ("II",  "AAK",  "*", "BH", "Ala Archa, Kyrgyzstan"),
    ("II",  "KURK", "*", "BH", "Kurchatov, Kazakhstan"),
    ("IC",  "LSA",  "*", "BH", "Lhasa, Tibet"),
    ("II",  "NIL",  "*", "BH", "Nilore, Pakistan"),
    ("IC",  "XAN",  "*", "BH", "Xi'an, China"),
]

# PS23 Makanchi array elements (short-period)
PS23_ELEMENTS = ["MK01", "MK02", "MK03", "MK04", "MK05",
                 "MK06", "MK07", "MK08", "MK09"]
PS23_NETWORK = "KZ"

# Station coordinates for PyGMT (from IRIS metadata)
STATION_COORDS = pd.DataFrame({
    "lon":     [87.6951, 81.9770, 82.2904, 79.2389, 76.9669, 78.4373,
                76.9661, 74.4942, 78.6189, 91.1500, 73.2686, 108.9230],
    "lat":     [43.8144, 46.7928, 46.7939, 41.5508, 43.2220, 42.4665,
                43.2194, 42.6389, 50.7154, 29.7025, 33.6506, 34.0310],
    "net_sta": ["IC.WMQ", "IU.MAKZ", "KZ.MKAR", "G.WUS", "KZ.PDGK",
                "KR.PRZ", "KZ.KNDC", "II.AAK", "II.KURK", "IC.LSA",
                "II.NIL", "IC.XAN"],
    "dist_km": [246, 780, 778, 839, 1020, 1000, 1020, 1179, 1264,
                1350, 1370, 1965],
    "desc":    ["Urumqi", "PS23 Makanchi", "Makanchi Array", "Wushi",
                "Podgonoye", "Karakol", "Almaty", "Ala Archa",
                "Kurchatov", "Lhasa", "Nilore", "Xi'an"],
})


# ═══════════════════════════════════════════════════════════════════════════════
# Data fetching
# ═══════════════════════════════════════════════════════════════════════════════

def fetch_broadband_data():
    """
    Fetch broadband waveforms from IRIS FDSN for all regional stations.

    Returns list of (net, sta, dist_km, az, trace, inv, desc) sorted by
    distance.  Instrument response is removed and output is velocity.
    """
    client = FDSNClient("IRIS", timeout=90)
    pre, post = 120, 900
    t1 = EVENT_TIME - pre
    t2 = EVENT_TIME + post
    results = []

    for net, sta, loc, chan_pref, desc in BB_STATIONS:
        label = f"{net}.{sta}"
        for cpref in [chan_pref, "LH", "SH"]:
            chan = f"{cpref}Z"
            try:
                st = client.get_waveforms(net, sta, loc, chan, t1, t2)
                if len(st) == 0:
                    continue
                st.merge(fill_value="interpolate")
                tr = st[0]
                inv = client.get_stations(
                    network=net, station=sta,
                    starttime=EVENT_TIME, endtime=EVENT_TIME,
                    level="response", channel=chan)
                sta_lat = inv[0][0].latitude
                sta_lon = inv[0][0].longitude
                d_m, az, _ = gps2dist_azimuth(EVENT_LAT, EVENT_LON,
                                               sta_lat, sta_lon)

                # Remove response -> velocity
                tr = tr.copy()
                tr.detrend("demean")
                tr.detrend("linear")
                tr.taper(0.05, type="cosine")
                try:
                    tr.remove_response(inventory=inv, output="VEL",
                                       water_level=60,
                                       pre_filt=[0.01, 0.02, 40, 45])
                except Exception:
                    pass

                results.append((net, sta, d_m / 1000.0, az, tr, inv, desc))
                print(f"  {label:12s}  dist={d_m/1000:.1f} km  "
                      f"sr={tr.stats.sampling_rate} Hz")
                break
            except Exception as exc:
                if cpref == "SH":
                    print(f"  {label:12s}  FAILED: {str(exc)[:60]}")

    results.sort(key=lambda x: x[2])
    return results


def fetch_ps23_array():
    """
    Fetch PS23 Makanchi array data (9 SHZ elements + MKAR BHZ).

    Returns (stream_shz, stream_all, inventory).
    """
    client = FDSNClient("IRIS", timeout=90)
    t1 = EVENT_TIME - 60
    t2 = EVENT_TIME + 300

    sta_list = ",".join(PS23_ELEMENTS + ["MKAR"])
    st = client.get_waveforms(
        network=PS23_NETWORK, station=sta_list,
        location="*", channel="*HZ",
        starttime=t1, endtime=t2)
    print(f"  PS23 array: {len(st)} traces downloaded")

    inv = client.get_stations(
        network=PS23_NETWORK, station=sta_list,
        starttime=t1, endtime=t2, level="channel")

    # Attach coordinates
    for tr in st:
        for net in inv:
            for sta_obj in net:
                if sta_obj.code == tr.stats.station:
                    tr.stats.coordinates = {
                        "latitude": sta_obj.latitude,
                        "longitude": sta_obj.longitude,
                        "elevation": sta_obj.elevation / 1000.0,
                    }
                    break

    st_shz = st.select(channel="SHZ").copy()

    # Also add IU.MAKZ broadband
    st_all = st.copy()
    try:
        st_makz = client.get_waveforms(
            network="IU", station="MAKZ",
            location="00", channel="BHZ",
            starttime=t1, endtime=t2)
        st_all += st_makz
    except Exception:
        pass

    return st_shz, st_all, inv


# ═══════════════════════════════════════════════════════════════════════════════
# Array utilities
# ═══════════════════════════════════════════════════════════════════════════════

def get_array_coords(inventory):
    """Extract station coordinates from inventory."""
    coords = {}
    for net in inventory:
        for sta in net:
            coords[sta.code] = (sta.latitude, sta.longitude, sta.elevation)
    lats = [c[0] for c in coords.values()]
    lons = [c[1] for c in coords.values()]
    return coords, np.mean(lats), np.mean(lons)


def compute_relative_positions(coords, ref_lat, ref_lon):
    """Geographic coords -> km offsets from reference point."""
    positions = {}
    for code, (lat, lon, elev) in coords.items():
        dy = (lat - ref_lat) * 111.19
        dx = (lon - ref_lon) * 111.19 * np.cos(np.radians(ref_lat))
        positions[code] = (dx, dy)
    return positions


def compute_backazimuth(array_lat, array_lon, ev_lat, ev_lon):
    """Backazimuth from array to event (degrees)."""
    dlon = np.radians(ev_lon - array_lon)
    y = np.sin(dlon) * np.cos(np.radians(ev_lat))
    x = (np.cos(np.radians(array_lat)) * np.sin(np.radians(ev_lat))
         - np.sin(np.radians(array_lat)) * np.cos(np.radians(ev_lat))
         * np.cos(dlon))
    return (np.degrees(np.arctan2(y, x)) + 360) % 360


def delay_and_sum_beam(stream, coords, ref_lat, ref_lon, baz_deg, slow_skm):
    """
    Delay-and-sum beamforming.

    Parameters
    ----------
    stream : obspy.Stream
        Array element traces.
    coords : dict
        Station -> (lat, lon, elev_m).
    baz_deg, slow_skm : float
        Backazimuth (degrees) and slowness (s/km).

    Returns
    -------
    beam : ndarray
        Beamformed trace.
    times : ndarray
        Time axis in seconds.
    """
    baz_rad = np.radians(baz_deg)
    positions = compute_relative_positions(coords, ref_lat, ref_lon)
    sx = slow_skm * np.sin(baz_rad)
    sy = slow_skm * np.cos(baz_rad)

    traces = []
    delays = []
    for tr in stream:
        sta = tr.stats.station
        if sta not in positions:
            continue
        dx, dy = positions[sta]
        delay_s = dx * sx + dy * sy
        traces.append(tr.data.astype(float))
        delays.append(delay_s)

    if not traces:
        return np.array([]), np.array([])

    sr = stream[0].stats.sampling_rate
    n = min(len(d) for d in traces)
    beam = np.zeros(n)
    for data, delay in zip(traces, delays):
        shift = int(round(delay * sr))
        beam += np.roll(data[:n], -shift)
    beam /= len(traces)
    return beam, np.arange(n) / sr


def find_minimum_lag_stack(stream, coords, ref_lat, ref_lon, baz_deg,
                           slow_range=(0.05, 0.25), n_slow=80):
    """
    Beamform at a grid of slowness values and return the stack that
    maximises the peak beam amplitude (minimum temporal lag = best alignment).

    Returns (best_beam, best_times, best_slow, all_beams, slownesses).
    """
    slownesses = np.linspace(slow_range[0], slow_range[1], n_slow)
    best_peak = 0
    best_beam = None
    best_times = None
    best_slow = slownesses[0]
    all_beams = []

    for slow in slownesses:
        beam, times = delay_and_sum_beam(
            stream, coords, ref_lat, ref_lon, baz_deg, slow)
        all_beams.append(beam)
        peak = np.max(np.abs(beam)) if len(beam) > 0 else 0
        if peak > best_peak:
            best_peak = peak
            best_beam = beam.copy()
            best_times = times.copy()
            best_slow = slow

    return best_beam, best_times, best_slow, all_beams, slownesses


# ═══════════════════════════════════════════════════════════════════════════════
# Figure 1: PyGMT Station Map
# ═══════════════════════════════════════════════════════════════════════════════

def fig01_pygmt_station_map(outdir):
    """PyGMT map of event and all seismic stations."""
    print("[1/5] PyGMT station map ...")
    region = [65, 115, 25, 55]
    projection = "M18c"

    fig = pygmt.Figure()

    fig.grdimage(
        grid="@earth_relief_30s",
        region=region, projection=projection,
        shading=True, cmap="geo")

    fig.basemap(
        region=region, projection=projection,
        frame=["WSne+tLop Nor 2020-06-22 Event and Seismic Stations",
               "xa10f5g10", "ya10f5g10"])

    fig.coast(
        region=region, projection=projection,
        resolution="l", water="lightblue",
        borders="1/0.5p,gray40", shorelines="thin,gray30",
        area_thresh=5000)

    # Great-circle paths
    for _, sta in STATION_COORDS.iterrows():
        fig.plot(x=[EVENT_LON, sta["lon"]], y=[EVENT_LAT, sta["lat"]],
                 pen="0.5p,gray70,--")

    # Stations
    fig.plot(x=STATION_COORDS["lon"], y=STATION_COORDS["lat"],
             style="t0.35c", fill="blue", pen="0.5p,black")

    # Event star
    fig.plot(x=[EVENT_LON], y=[EVENT_LAT],
             style="a0.5c", fill="red", pen="0.5p,black")

    # Labels
    for _, sta in STATION_COORDS.iterrows():
        fig.text(x=sta["lon"], y=sta["lat"] + 0.6,
                 text=sta["net_sta"], font="6p,Helvetica,black",
                 justify="BC")

    fig.text(x=EVENT_LON, y=EVENT_LAT - 0.8,
             text="Lop Nor", font="8p,Helvetica-Bold,red", justify="TC")

    # Info box
    info_lines = [
        "Event Information:",
        f"Date: 2020-06-22  ~09:18 UTC",
        f"Location: {EVENT_LAT:.3f}@.N, {EVENT_LON:.3f}@.E",
        "Depth: ~300 m (tunnel)",
        f"mb @~\\176@~ {MB_OBSERVED}",
        "Alleged decoupled nuclear test",
    ]
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt',
                                     delete=False) as f:
        for i, line in enumerate(info_lines):
            font = "8p,Helvetica-Bold,black" if i == 0 \
                else "7p,Helvetica,black"
            f.write(f"L {font} L {line}\n")
        info_file = f.name
    fig.legend(spec=info_file, position="JTL+w5.5c+o0.5c/0.5c",
               box="+p0.5p+gwhite@30+s")
    os.unlink(info_file)

    # Symbol legend
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt',
                                     delete=False) as f:
        f.write("S 0.2c a 0.4c red 0.5p,black 0.8c Event Origin\n")
        f.write("S 0.2c t 0.3c blue 0.5p,black 0.8c Seismic Station\n")
        legend_file = f.name
    fig.legend(spec=legend_file, position="JBL+w4c+o0.5c/0.5c",
               box="+p0.5p+gwhite+s")
    os.unlink(legend_file)

    # Scale bar
    fig.basemap(region=region, projection=projection,
                map_scale="g72/27+c41+w500k+f+l")

    # Inset
    with fig.inset(position="jTR+w5c+o0.3c", box="+p0.5p+gwhite"):
        fig.coast(region="g", projection=f"G{EVENT_LON}/{EVENT_LAT}/?",
                  resolution="c", land="gray80", water="lightblue",
                  shorelines="faint", area_thresh=10000)
        fig.plot(x=[EVENT_LON], y=[EVENT_LAT],
                 style="a0.3c", fill="red", pen="0.3p,black")

    p = os.path.join(outdir, "fig01_station_map.png")
    fig.savefig(p, dpi=300)
    print(f"  Saved {p}")


# ═══════════════════════════════════════════════════════════════════════════════
# Figure 2: Multi-filter record sections
# ═══════════════════════════════════════════════════════════════════════════════

def fig02_multifilter_record_section(outdir, bb_results):
    """
    Three-panel record section: each panel uses a different filter.
    All broadband stations, sorted by distance.
    """
    print("[2/5] Multi-filter record sections ...")
    n_sta = len(bb_results)
    if n_sta == 0:
        print("  No broadband data to plot.")
        return

    fig, axes = plt.subplots(1, 3, figsize=(24, max(10, 1.2 * n_sta)),
                             sharey=True)

    for col, (filt_label, filt_kw) in enumerate(FILTERS):
        ax = axes[col]

        for i, (net, sta, dist_km, az, tr_orig, inv, desc) in \
                enumerate(bb_results):
            tr = tr_orig.copy()
            tr.detrend("demean")
            tr.taper(0.05)
            tr.filter(**filt_kw)

            t_rel = np.arange(tr.stats.npts) / tr.stats.sampling_rate
            t_origin = tr.stats.starttime - EVENT_TIME
            t_plot = t_rel + t_origin
            data = tr.data.astype(float)
            peak = np.max(np.abs(data))
            if peak > 0:
                data = data / peak

            ax.plot(t_plot, data * 0.38 + i, "k-", lw=0.4)
            ax.fill_between(t_plot, i, data * 0.38 + i,
                            where=(data > 0), color="red", alpha=0.15)
            ax.fill_between(t_plot, i, data * 0.38 + i,
                            where=(data < 0), color="blue", alpha=0.15)

            if col == 0:
                ax.text(-115, i, f"{net}.{sta}\n({dist_km:.0f} km)",
                        fontsize=7, ha="right", va="center",
                        fontweight="bold")

        # Phase move-out lines
        distances = [r[2] for r in bb_results]
        for phase, vel in PHASE_VELOCITIES.items():
            t_phase = np.array(distances) / vel
            ax.plot(t_phase, np.arange(n_sta), "--",
                    color=PHASE_COLORS.get(phase, "gray"), lw=1.2,
                    alpha=0.7, label=phase)

        ax.set_xlabel("Time after origin (s)", fontsize=11)
        ax.set_xlim(-60, 600)
        ax.set_ylim(-0.5, n_sta - 0.5)
        ax.set_title(filt_label, fontsize=13, fontweight="bold")
        ax.grid(True, alpha=0.15)
        if col == 0:
            ax.legend(fontsize=8, loc="lower right")

    fig.suptitle("Lop Nor 2020-06-22 -- Multi-filter record sections "
                 "(all broadband stations)",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    p = os.path.join(outdir, "fig02_multifilter_record_section.png")
    fig.savefig(p, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {p}")


# ═══════════════════════════════════════════════════════════════════════════════
# Figure 3: Two-pulse timing display
# ═══════════════════════════════════════════════════════════════════════════════

def fig03_two_pulse_timing(outdir, bb_results):
    """
    Zoomed view of closest stations showing the ~12 s two-pulse structure.
    Three rows (one per filter), with the 3-4 closest stations overlaid.
    Time cut: P arrival - 5 s  to  P arrival + 40 s.
    """
    print("[3/5] Two-pulse timing display ...")
    n_sta = len(bb_results)
    if n_sta == 0:
        print("  No broadband data to plot.")
        return

    # Use closest 4 stations (or fewer if not enough)
    n_show = min(4, n_sta)
    close_results = bb_results[:n_show]

    fig, axes = plt.subplots(len(FILTERS), 1, figsize=(18, 5 * len(FILTERS)),
                             sharex=True)
    if len(FILTERS) == 1:
        axes = [axes]

    for row, (filt_label, filt_kw) in enumerate(FILTERS):
        ax = axes[row]

        for i, (net, sta, dist_km, az, tr_orig, inv, desc) in \
                enumerate(close_results):
            tr = tr_orig.copy()
            tr.detrend("demean")
            tr.taper(0.05)
            tr.filter(**filt_kw)

            t_rel = np.arange(tr.stats.npts) / tr.stats.sampling_rate
            t_origin = tr.stats.starttime - EVENT_TIME
            t_plot = t_rel + t_origin
            data = tr.data.astype(float)
            peak = np.max(np.abs(data))
            if peak > 0:
                data = data / peak

            ax.plot(t_plot, data * 0.35 + i, lw=0.8,
                    label=f"{net}.{sta} ({dist_km:.0f} km)")

        # Expected P arrival for the closest station
        dist_close = close_results[0][2]
        p_arrival = dist_close / PHASE_VELOCITIES["Pn"]

        # Mark the two pulses: P and P+12s
        ax.axvline(p_arrival, color="red", ls="--", lw=1.5, alpha=0.7,
                   label="Expected Pn")
        ax.axvline(p_arrival + 12, color="orange", ls=":", lw=1.5,
                   alpha=0.7, label="Pn + 12 s (second pulse)")

        # Shade the two-pulse window
        ax.axvspan(p_arrival - 1, p_arrival + 5, color="red",
                   alpha=0.05)
        ax.axvspan(p_arrival + 10, p_arrival + 18, color="orange",
                   alpha=0.05)

        ax.set_xlim(p_arrival - 10, p_arrival + 50)
        ax.set_ylabel(filt_label, fontsize=11, fontweight="bold")
        ax.legend(fontsize=8, loc="upper right", ncol=2)
        ax.grid(True, alpha=0.2)
        ax.axhline(0, color="gray", lw=0.3)

        if row == 0:
            ax.set_title("Two-pulse timing -- closest stations, "
                         "cut to show ~12 s separation",
                         fontsize=13, fontweight="bold")

    axes[-1].set_xlabel("Time after origin (s)", fontsize=12)

    fig.suptitle("Lop Nor 2020-06-22 -- Two-pulse signature display",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    p = os.path.join(outdir, "fig03_two_pulse_timing.png")
    fig.savefig(p, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {p}")


# ═══════════════════════════════════════════════════════════════════════════════
# Figure 4: PS23 array waveforms in three filter bands
# ═══════════════════════════════════════════════════════════════════════════════

def fig04_array_multifilter(outdir, stream_all):
    """
    PS23 array element waveforms in the three requested filter bands.
    Zoomed to the P-arrival window.
    """
    print("[4/5] PS23 array multi-filter display ...")
    # Expected P at ~780 km via Pn
    p_arrival = 780.0 / PHASE_VELOCITIES["Pn"]

    sorted_codes = sorted(set(tr.stats.station for tr in stream_all))
    n_sta = len(sorted_codes)

    fig, axes = plt.subplots(1, 3, figsize=(24, max(10, 1.0 * n_sta)),
                             sharey=True)

    for col, (filt_label, filt_kw) in enumerate(FILTERS):
        ax = axes[col]
        st_filt = stream_all.copy()
        st_filt.detrend("demean")
        st_filt.taper(0.05)
        st_filt.filter(**filt_kw)

        for i, code in enumerate(sorted_codes):
            trs = st_filt.select(station=code)
            if not trs:
                continue
            # Pick BHZ or SHZ -- prefer SHZ for array elements
            tr = trs[0]
            t_rel = np.arange(tr.stats.npts) / tr.stats.sampling_rate
            t_origin = tr.stats.starttime - EVENT_TIME
            t_plot = t_rel + t_origin
            data = tr.data.astype(float)
            peak = np.max(np.abs(data))
            if peak > 0:
                data = data / peak
            ax.plot(t_plot, data * 0.4 + i, "k-", lw=0.5)

            if col == 0:
                ax.text(p_arrival - 14, i, code, fontsize=7, ha="right",
                        va="center", fontweight="bold")

        ax.axvline(p_arrival, color="red", ls="--", lw=1.5, alpha=0.7,
                   label=f"Expected P ({p_arrival:.0f} s)")
        ax.axvline(p_arrival + 12, color="orange", ls=":", lw=1.5,
                   alpha=0.7, label="P + 12 s")
        ax.set_xlim(p_arrival - 10, p_arrival + 50)
        ax.set_xlabel("Time after origin (s)", fontsize=11)
        ax.set_title(filt_label, fontsize=13, fontweight="bold")
        ax.set_yticks([])
        ax.grid(True, alpha=0.2)
        if col == 0:
            ax.legend(fontsize=8)

    fig.suptitle("PS23 Makanchi array -- Multi-filter display "
                 "(zoomed to P-arrival + two-pulse window)",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    p = os.path.join(outdir, "fig04_array_multifilter.png")
    fig.savefig(p, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {p}")


# ═══════════════════════════════════════════════════════════════════════════════
# Figure 5: Array stack at minimum lag
# ═══════════════════════════════════════════════════════════════════════════════

def fig05_array_stack_minimum_lag(outdir, stream_shz, inventory):
    """
    Beamform PS23 array across a grid of slowness values.
    Show: (a) slowness-power curve, (b) best beam vs single element for each
    filter band, (c) beam envelope.
    """
    print("[5/5] Array stack at minimum lag ...")
    coords, ref_lat, ref_lon = get_array_coords(inventory)
    baz = compute_backazimuth(ref_lat, ref_lon, EVENT_LAT, EVENT_LON)
    p_arrival = 780.0 / PHASE_VELOCITIES["Pn"]

    fig = plt.figure(figsize=(20, 16))
    gs = GridSpec(3, 2, figure=fig, hspace=0.4, wspace=0.3)

    # ── (a) Slowness scan: peak power vs slowness ──
    ax_scan = fig.add_subplot(gs[0, 0])

    # Use 1-10 Hz filtered data for slowness scan
    st_scan = stream_shz.copy()
    st_scan.detrend("demean")
    st_scan.taper(0.05)
    st_scan.filter("bandpass", freqmin=1.0, freqmax=10.0,
                   corners=4, zerophase=True)

    best_beam, best_times, best_slow, all_beams, slownesses = \
        find_minimum_lag_stack(st_scan, coords, ref_lat, ref_lon, baz)

    # Compute peak power in P window for each slowness
    t0_off = st_scan[0].stats.starttime - EVENT_TIME
    peaks = []
    for beam in all_beams:
        if len(beam) == 0:
            peaks.append(0)
            continue
        btimes = np.arange(len(beam)) / st_scan[0].stats.sampling_rate + t0_off
        p_mask = (btimes > p_arrival - 2) & (btimes < p_arrival + 20)
        peaks.append(np.max(np.abs(beam[p_mask])) if np.any(p_mask) else 0)

    ax_scan.plot(slownesses, peaks / max(max(peaks), 1e-30), "k-", lw=2)
    ax_scan.axvline(best_slow, color="red", ls="--", lw=1.5,
                    label=f"Best: {best_slow:.4f} s/km "
                          f"({1/best_slow:.1f} km/s)")
    ax_scan.axvline(1 / 7.9, color="blue", ls=":", lw=1.5, alpha=0.7,
                    label="Pn (7.9 km/s)")
    ax_scan.set_xlabel("Slowness (s/km)", fontsize=11, fontweight="bold")
    ax_scan.set_ylabel("Normalised peak amplitude", fontsize=11,
                       fontweight="bold")
    ax_scan.set_title("(a) Slowness scan -- peak beam power in P window",
                      fontsize=12, fontweight="bold")
    ax_scan.legend(fontsize=9)
    ax_scan.grid(True, alpha=0.2)

    # ── (b) Best beam vs single element for each filter ──
    for row, (filt_label, filt_kw) in enumerate(FILTERS):
        ax = fig.add_subplot(gs[1 + row // 2, row % 2 + (0 if row < 2 else 0)])
        # Actually, let's use a simpler layout: 3 panels in row 1 col 1,
        # and row 2 col 0, row 2 col 1
        # Re-do layout: (0,1), (1,0), (1,1), (2,0), (2,1)
        # positions: row=0 col=1 for scan, then rows 1-2 for 3 filters
        pass

    # Simpler: three panels below the scan
    # Reuse gs more carefully
    plt.close(fig)

    fig = plt.figure(figsize=(22, 20))
    gs = GridSpec(4, 2, figure=fig, hspace=0.4, wspace=0.3,
                  height_ratios=[1, 1, 1, 1])

    # (a) Slowness scan
    ax_scan = fig.add_subplot(gs[0, 0])
    ax_scan.plot(slownesses, np.array(peaks) / max(max(peaks), 1e-30),
                 "k-", lw=2)
    ax_scan.axvline(best_slow, color="red", ls="--", lw=1.5,
                    label=f"Best: {best_slow:.4f} s/km "
                          f"({1/best_slow:.1f} km/s)")
    ax_scan.axvline(1 / 7.9, color="blue", ls=":", lw=1.5, alpha=0.7,
                    label="Pn (7.9 km/s)")
    ax_scan.set_xlabel("Slowness (s/km)", fontsize=11, fontweight="bold")
    ax_scan.set_ylabel("Normalised peak amp.", fontsize=11, fontweight="bold")
    ax_scan.set_title("(a) Slowness scan -- peak beam power in P window",
                      fontsize=12, fontweight="bold")
    ax_scan.legend(fontsize=9)
    ax_scan.grid(True, alpha=0.2)

    # (b) Array layout with backazimuth
    ax_layout = fig.add_subplot(gs[0, 1])
    positions = compute_relative_positions(coords, ref_lat, ref_lon)
    for code, (dx, dy) in positions.items():
        color = "#d62728" if code == "MKAR" else "#1f77b4"
        marker = "s" if code == "MKAR" else "^"
        ax_layout.plot(dx, dy, marker, color=color, ms=10,
                       markeredgecolor="k", markeredgewidth=0.8, zorder=5)
        ax_layout.text(dx + 0.15, dy + 0.15, code, fontsize=7,
                       fontweight="bold", color="#333")
    # Backazimuth arrow
    arrow_len = 2.5
    ax_layout.annotate("",
                       xy=(arrow_len * np.sin(np.radians(baz)),
                           arrow_len * np.cos(np.radians(baz))),
                       xytext=(0, 0),
                       arrowprops=dict(arrowstyle="->", color="red", lw=2.5))
    ax_layout.text(0, -5.0,
                   f"To Lop Nor (BAZ={baz:.1f}\u00b0)",
                   fontsize=9, color="red", fontweight="bold", ha="center")
    ax_layout.set_xlabel("East-West (km)", fontsize=10)
    ax_layout.set_ylabel("North-South (km)", fontsize=10)
    ax_layout.set_title("(b) PS23 Array Layout", fontsize=12,
                        fontweight="bold")
    ax_layout.set_aspect("equal")
    ax_layout.grid(True, alpha=0.3)

    # (c-e) Best-slowness beam vs single element, one per filter
    panel_labels = ["(c)", "(d)", "(e)"]
    for idx, (filt_label, filt_kw) in enumerate(FILTERS):
        ax = fig.add_subplot(gs[1 + idx, :])

        # Filter and beam at best slowness
        st_f = stream_shz.copy()
        st_f.detrend("demean")
        st_f.taper(0.05)
        st_f.filter(**filt_kw)

        beam_f, btimes_f = delay_and_sum_beam(
            st_f, coords, ref_lat, ref_lon, baz, best_slow)
        btimes_f_origin = btimes_f + t0_off

        # Single element (MK01)
        tr_single = st_f.select(station="MK01")
        if not tr_single:
            tr_single = st_f[0:1]
        tr_s = tr_single[0]
        s_times = (np.arange(tr_s.stats.npts) / tr_s.stats.sampling_rate
                   + (tr_s.stats.starttime - EVENT_TIME))
        s_data = tr_s.data.astype(float)

        # Normalise both to single-element peak in P window
        s_pmask = (s_times > p_arrival - 2) & (s_times < p_arrival + 30)
        s_peak = np.max(np.abs(s_data[s_pmask])) if np.any(s_pmask) else 1.0
        if s_peak == 0:
            s_peak = 1.0

        ax.plot(s_times, s_data / s_peak, color="#bbbbbb", lw=0.5,
                label="Single element (MK01)")
        ax.plot(btimes_f_origin[:len(beam_f)], beam_f / s_peak,
                "k-", lw=1.2,
                label=f"Beam (9 elem, slow={best_slow:.4f} s/km)")

        # Beam envelope
        if len(beam_f) > 10:
            env = np.abs(hilbert(beam_f))
            ax.plot(btimes_f_origin[:len(env)], env / s_peak,
                    "r-", lw=0.8, alpha=0.5, label="Beam envelope")

        # Mark two pulses
        ax.axvline(p_arrival, color="red", ls="--", lw=1.5, alpha=0.7,
                   label="Expected Pn")
        ax.axvline(p_arrival + 12, color="orange", ls=":", lw=1.5,
                   alpha=0.7, label="Pn + 12 s")

        # SNR
        b_pmask = (btimes_f_origin[:len(beam_f)] > p_arrival - 2) & \
                  (btimes_f_origin[:len(beam_f)] < p_arrival + 20)
        b_nmask = (btimes_f_origin[:len(beam_f)] > 20) & \
                  (btimes_f_origin[:len(beam_f)] < p_arrival - 20)
        if np.any(b_pmask) and np.any(b_nmask):
            snr_beam = (np.sqrt(np.mean(beam_f[b_pmask]**2))
                        / max(np.sqrt(np.mean(beam_f[b_nmask]**2)), 1e-30))
        else:
            snr_beam = 0

        ax.text(0.01, 0.95,
                f"SNR (beam, P-window): {snr_beam:.1f}",
                transform=ax.transAxes, fontsize=10, fontweight="bold",
                va="top", bbox=dict(fc="white", ec="#ccc", alpha=0.9))

        ax.set_xlim(p_arrival - 15, p_arrival + 55)
        ax.set_xlabel("Time after origin (s)", fontsize=11)
        ax.set_ylabel("Normalised amplitude", fontsize=11)
        ax.set_title(f"{panel_labels[idx]} Minimum-lag beam stack -- "
                     f"{filt_label}", fontsize=12, fontweight="bold")
        ax.legend(fontsize=8, loc="upper right", ncol=3)
        ax.grid(True, alpha=0.2)
        ax.axhline(0, color="gray", lw=0.3)

    fig.suptitle(f"PS23 Array stack at minimum lag "
                 f"(BAZ={baz:.1f}\u00b0, best slow={best_slow:.4f} s/km, "
                 f"vel={1/best_slow:.1f} km/s)",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    p = os.path.join(outdir, "fig05_array_stack_minimum_lag.png")
    fig.savefig(p, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {p}")


# ═══════════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    os.makedirs(OUTDIR, exist_ok=True)

    print("=" * 72)
    print("  Lop Nor 2020-06-22 -- All Data Multi-Filter Plots")
    print("=" * 72)
    print(f"  Origin:   {EVENT_TIME}")
    print(f"  Location: {EVENT_LAT:.3f}\u00b0N  {EVENT_LON:.3f}\u00b0E")
    print(f"  mb:       ~{MB_OBSERVED}")
    print(f"  Filters:  {', '.join(f[0] for f in FILTERS)}")
    print(f"  Output:   {OUTDIR}")
    print()

    # ── Figure 1: PyGMT station map (no data needed) ──
    fig01_pygmt_station_map(OUTDIR)
    print()

    # ── Fetch broadband data ──
    print("Fetching broadband data from IRIS/FDSN ...")
    bb_results = fetch_broadband_data()
    print(f"  {len(bb_results)} station(s) retrieved\n")

    if bb_results:
        fig02_multifilter_record_section(OUTDIR, bb_results)
        print()
        fig03_two_pulse_timing(OUTDIR, bb_results)
        print()

    # ── Fetch PS23 array data ──
    print("Fetching PS23 Makanchi array data ...")
    st_shz, st_all, inv = fetch_ps23_array()
    print()

    if len(st_all) > 0:
        fig04_array_multifilter(OUTDIR, st_all)
        print()

    if len(st_shz) > 0:
        fig05_array_stack_minimum_lag(OUTDIR, st_shz, inv)
        print()

    print("=" * 72)
    print(f"  Done.  All figures saved to: {OUTDIR}")
    print("=" * 72)


if __name__ == "__main__":
    main()
