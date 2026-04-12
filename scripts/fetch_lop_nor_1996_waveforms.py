#!/usr/bin/env python3
"""
Fetch and display real seismic waveforms from the 1996-07-29 Lop Nor nuclear test.

Loads local KNDC archive SAC files (26 regional short-period stations from the
Kazakhstan/Kyrgyzstan network) and also attempts to download broadband data from
IRIS FDSN for any available stations that were operational in 1996.

The KNDC archive contains short-period (SKM/SKD) recordings from Soviet-era
Central Asian stations at distances of 761–1531 km.  These are supplemented by
any available broadband FDSN data from IMS, IRIS GSN, and Chinese networks.

Event:  1996-07-29  01:48:57.170 UTC, Lop Nor (~41.7163°N, 88.3748°E)
        45th and final confirmed Chinese underground nuclear test
        mb 4.90 (ISC), estimated yield ~5 kt (coupled in granite tunnel)

Usage:
    python scripts/fetch_lop_nor_1996_waveforms.py
"""

import os
import sys

from obspy import UTCDateTime, read
from obspy.geodetics import gps2dist_azimuth

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fsrm.data_fetching import fetch_raw_waveforms, process_trace
from fsrm.fetch_plotting import (
    fetch_plot_record_section,
    fetch_plot_zoomed,
    fetch_plot_distance_scaled,
    fetch_plot_long_period,
    fetch_plot_spectrograms,
    fetch_plot_closeup,
)

# ──────────────────────────────────────────────────────────────────────────────
# Event parameters
# ──────────────────────────────────────────────────────────────────────────────
EVENT_TIME = UTCDateTime("1996-07-29T01:48:57.170")
EVENT_LAT = 41.7163
EVENT_LON = 88.3748
EVENT_DEPTH_KM = 0.4

PRE_ORIGIN = 120
POST_ORIGIN = 600
FREQMIN = 0.5
FREQMAX = 8.0

OUTDIR = os.path.join(os.path.dirname(os.path.dirname(__file__)),
                      "figures", "lop_nor_1996")

EVENT_INFO = {
    "name": "Lop Nor Final Test 1996-07-29",
    "date": "1996-07-29",
    "time": "01:48:57 UTC",
    "lat": EVENT_LAT,
    "lon": EVENT_LON,
    "depth_km": EVENT_DEPTH_KM,
    "mb": "4.90",
    "yield_text": "~5 kt (coupled, from mb-yield)",
}

PHASE_VELOCITIES = {
    "Pn": 7.9,
    "Pg": 6.1,
    "Sn": 4.5,
    "Lg": 3.5,
}

# Key stations from the KNDC archive (nearest, best SNR)
KEY_STATIONS = ["DZHR", "UZB", "PRZ", "SATY", "KURM"]

# FDSN stations to attempt downloading (some may not have been operational in 1996)
FDSN_STATIONS = [
    ("II",  "KURK", "*", "BH",  "Kurchatov, Kazakhstan"),
    ("II",  "AAK",  "*", "BH",  "Ala Archa, Kyrgyzstan"),
    ("IU",  "MAKZ", "*", "BH",  "Makanchi, Kazakhstan (PS23)"),
    ("IC",  "WMQ",  "*", "BH",  "Urumqi, China"),
    ("KR",  "PRZ",  "*", "BH",  "Karakol, Kyrgyzstan"),
    ("KZ",  "KNDC", "*", "BH",  "KNDC Almaty, Kazakhstan"),
    ("IC",  "LSA",  "*", "BH",  "Lhasa, Tibet"),
    ("IC",  "XAN",  "*", "BH",  "Xi'an, China"),
    ("II",  "NIL",  "*", "BH",  "Nilore, Pakistan"),
    ("IU",  "ULN",  "*", "BH",  "Ulaanbaatar, Mongolia"),
]

_REPO_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _parse_site_coords():
    """Parse KNDC .site file to build station -> (lat, lon) lookup."""
    site_file = os.path.join(_REPO_DIR, "data", "historic_nuclear",
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


def load_kndc_archive():
    """
    Load all Z-component SAC files from the KNDC archive for this event.

    Returns list of (net, sta, dist_km, az, trace, inv, description) tuples
    compatible with fetch_raw_waveforms output (inv=None for KNDC data).
    """
    data_dir = os.path.join(
        _REPO_DIR,
        "data", "historic_nuclear", "03.Lopnor", "19960729.0148", "wf")

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

            # Extract station name
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

            # Get distance/azimuth
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

            results.append(("KZ", sta, dist_km, az, tr, None,
                           f"{sta} (KNDC, {dist_km:.0f} km)"))
        except Exception as exc:
            print(f"  WARNING: Could not read {fname}: {exc}")

    results.sort(key=lambda x: x[2])
    return results


def process_kndc_traces(raw_results, freqmin, freqmax):
    """
    Apply basic processing to KNDC traces (no instrument response removal
    since FAP response files require custom deconvolution, but apply detrend
    and bandpass).

    Returns processed list (net, sta, dist_km, az, tr_proc, inv, desc).
    """
    processed = []
    for net, sta, dist_km, az, tr, inv, desc in raw_results:
        try:
            tr_proc = tr.copy()
            tr_proc.detrend("demean")
            tr_proc.detrend("linear")
            tr_proc.taper(max_percentage=0.05, type="cosine")
            tr_proc.filter("bandpass", freqmin=freqmin, freqmax=freqmax,
                           corners=4, zerophase=True)
            processed.append((net, sta, dist_km, az, tr_proc, inv, desc))
            print(f"  \u2713 {sta} (KNDC, {dist_km:.0f} km)")
        except Exception as exc:
            print(f"  \u2717 {sta}: {exc}")
    return processed


def main():
    os.makedirs(OUTDIR, exist_ok=True)

    print("=" * 72)
    print("  Lop Nor 1996-07-29 — Waveform Gather (KNDC Archive + FDSN/IRIS)")
    print("  45th & Final Chinese Underground Nuclear Test")
    print("=" * 72)
    print(f"  Origin time:  {EVENT_TIME}")
    print(f"  Location:     {EVENT_LAT:.4f}\u00b0N  {EVENT_LON:.4f}\u00b0E")
    print(f"  Depth:        ~{EVENT_DEPTH_KM} km (horizontal tunnel)")
    print(f"  mb:           4.90 (ISC)")
    print(f"  Window:       {PRE_ORIGIN} s before \u2192 {POST_ORIGIN} s after origin")
    print(f"  Filter:       {FREQMIN}\u2013{FREQMAX} Hz bandpass")
    print()

    # ── 1. Load KNDC archive ──
    print("Loading KNDC archive SAC files ...")
    kndc_raw = load_kndc_archive()
    print(f"  Found {len(kndc_raw)} Z-component trace(s)")
    print()

    print("Processing KNDC traces (detrend + bandpass) ...")
    kndc_proc = process_kndc_traces(kndc_raw, FREQMIN, FREQMAX)
    print(f"  Processed {len(kndc_proc)} trace(s)")
    print()

    # ── 2. Download FDSN data ──
    print("Downloading FDSN waveforms (may fail for 1996 data) ...")
    fdsn_raw = fetch_raw_waveforms(
        EVENT_TIME, EVENT_LAT, EVENT_LON, FDSN_STATIONS,
        pre_origin=PRE_ORIGIN, post_origin=POST_ORIGIN)
    print(f"  Retrieved {len(fdsn_raw)} station(s) from FDSN")
    print()

    # Process FDSN data (with instrument response removal)
    fdsn_proc = []
    if fdsn_raw:
        print("Processing FDSN traces (response removal + bandpass) ...")
        for net, sta, dist_km, az, tr, inv, desc in fdsn_raw:
            try:
                tr_proc = process_trace(tr, inv, FREQMIN, FREQMAX)
                fdsn_proc.append((net, sta, dist_km, az, tr_proc, inv, desc))
                print(f"  \u2713 {net}.{sta}")
            except Exception as exc:
                print(f"  \u2717 {net}.{sta}: {exc}")
        print(f"  Processed {len(fdsn_proc)} FDSN trace(s)")
        print()

    # ── 3. Merge results ──
    fdsn_stas = {x[1] for x in fdsn_proc}
    all_results = list(fdsn_proc)
    for item in kndc_proc:
        if item[1] not in fdsn_stas:
            all_results.append(item)
    all_results.sort(key=lambda x: x[2])
    print(f"Total: {len(all_results)} unique station(s) (KNDC + FDSN)")
    print()

    if not all_results:
        print("ERROR: No waveforms available. Exiting.")
        sys.exit(1)

    # ── 4. Generate figures ──
    kw = dict(event_info=EVENT_INFO, phase_velocities=PHASE_VELOCITIES,
              freqmin=FREQMIN, fmax=FREQMAX, event_depth_km=EVENT_DEPTH_KM)

    print("Generating figures ...")

    print("  [1/6] Record section ...")
    fetch_plot_record_section(all_results, EVENT_TIME, OUTDIR, **kw)

    print("  [2/6] Zoomed record section ...")
    fetch_plot_zoomed(all_results, EVENT_TIME, OUTDIR, **kw,
                      tmin=-30, tmax=300)

    print("  [3/6] Distance-scaled record section ...")
    fetch_plot_distance_scaled(all_results, EVENT_TIME, OUTDIR, **kw,
                               tmin=-20, tmax=400)

    print("  [4/6] Long-period filter (surface waves) ...")
    fetch_plot_long_period(kndc_raw, EVENT_TIME, OUTDIR,
                           event_info=EVENT_INFO,
                           phase_velocities=PHASE_VELOCITIES,
                           event_depth_km=EVENT_DEPTH_KM)

    print("  [5/6] Spectrograms ...")
    fetch_plot_spectrograms(all_results, EVENT_TIME, OUTDIR, **kw)

    print("  [6/6] Close-up (key stations) ...")
    fetch_plot_closeup(all_results, EVENT_TIME, OUTDIR, **kw,
                       key_stations=KEY_STATIONS)

    print()
    print("=" * 72)
    print(f"  All figures saved to: {OUTDIR}")
    print("=" * 72)


if __name__ == "__main__":
    main()
