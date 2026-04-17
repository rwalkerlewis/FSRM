#!/usr/bin/env python3
"""
Fetch and display real seismic waveforms from the June 22, 2020 Lop Nor event.

Downloads broadband data from IRIS FDSN for all confirmed stations, removes
instrument response, bandpass-filters, and produces a distance-sorted record
section (waveform gather).

Event:  2020-06-22 ~09:18 UTC, Lop Nor test site (~41.735N, 88.730E)
        Detected at PS23 Makanchi (Kazakhstan), mb ~ 2.75

Usage:
    python scripts/fetch_lop_nor_2020_waveforms.py
"""

import os
import sys

from obspy import UTCDateTime

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
EVENT_TIME = UTCDateTime("2020-06-22T09:18:00")
EVENT_LAT = 41.735
EVENT_LON = 88.730
EVENT_DEPTH_KM = 0.3

PRE_ORIGIN = 120
POST_ORIGIN = 900
FREQMIN = 0.5
FREQMAX = 8.0

OUTDIR = os.path.join(os.path.dirname(os.path.dirname(__file__)),
                      "figures", "lop_nor_2020")

EVENT_INFO = {
    "name": "Lop Nor Event 2020-06-22",
    "date": "2020-06-22",
    "time": "~09:18 UTC",
    "lat": EVENT_LAT,
    "lon": EVENT_LON,
    "depth_km": EVENT_DEPTH_KM,
    "mb": "2.75",
    "yield_text": "unknown (possible ~1 kt)",
}

PHASE_VELOCITIES = {
    "Pn": 8.1,
    "Pg": 6.1,
    "Sn": 4.6,
    "Lg": 3.5,
}

KEY_STATIONS = ["MAKZ", "MKAR", "WMQ"]

STATIONS = [
    ("IC", "WMQ",  "*", "BH",  "Urumqi, China"),
    ("IU", "MAKZ", "*", "BH",  "Makanchi, Kazakhstan (PS23)"),
    ("KZ", "MKAR", "*", "BH",  "Makanchi Array, Kazakhstan"),
    ("G",  "WUS",  "*", "BH",  "Wushi, China"),
    ("KZ", "PDGK", "*", "BH",  "Podgonoye, Kazakhstan"),
    ("KR", "PRZ",  "*", "BH",  "Karakol, Kyrgyzstan"),
    ("KZ", "KNDC", "*", "BH",  "KNDC Almaty, Kazakhstan"),
    ("II", "AAK",  "*", "BH",  "Ala Archa, Kyrgyzstan"),
    ("II", "KURK", "*", "BH",  "Kurchatov, Kazakhstan"),
    ("IC", "LSA",  "*", "BH",  "Lhasa, Tibet"),
    ("II", "NIL",  "*", "BH",  "Nilore, Pakistan"),
    ("IC", "XAN",  "*", "BH",  "Xi'an, China"),
]


def main():
    os.makedirs(OUTDIR, exist_ok=True)

    print("=" * 72)
    print("  Lop Nor 2020-06-22 Event — Waveform Gather from FDSN/IRIS")
    print("=" * 72)
    print(f"  Origin time:  {EVENT_TIME}")
    print(f"  Location:     {EVENT_LAT:.3f}°N  {EVENT_LON:.3f}°E")
    print(f"  Depth:        ~{EVENT_DEPTH_KM} km")
    print(f"  Window:       {PRE_ORIGIN} s before → {POST_ORIGIN} s after origin")
    print(f"  Filter:       {FREQMIN}–{FREQMAX} Hz bandpass")
    print()

    # 1. Download
    print("Downloading waveforms ...")
    raw_results = fetch_raw_waveforms(
        EVENT_TIME, EVENT_LAT, EVENT_LON, STATIONS,
        pre_origin=PRE_ORIGIN, post_origin=POST_ORIGIN)
    print(f"\n  Retrieved {len(raw_results)} station(s)\n")

    if not raw_results:
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

    kw = dict(event_info=EVENT_INFO, phase_velocities=PHASE_VELOCITIES,
              freqmin=FREQMIN, fmax=FREQMAX, event_depth_km=EVENT_DEPTH_KM)

    print("Generating waveform gather (record section – full window) ...")
    fetch_plot_record_section(
        processed, EVENT_TIME, OUTDIR, **kw,
        pre_origin=PRE_ORIGIN, post_origin=POST_ORIGIN,
        title_extra="instrument-response-corrected velocity",
        fname="observed_waveform_gather")

    print("Generating waveform gather (record section – zoomed) ...")
    fetch_plot_zoomed(
        processed, EVENT_TIME, OUTDIR, **kw,
        fname="observed_waveform_gather_zoomed")

    print("Generating distance-scaled record section with phase arrivals ...")
    fetch_plot_distance_scaled(
        processed, EVENT_TIME, OUTDIR, **kw,
        fname="observed_waveform_gather_distance_scaled")

    print("Generating long-period filtered record section (0.02–0.10 Hz) ...")
    fetch_plot_long_period(
        raw_results, EVENT_TIME, OUTDIR,
        event_info=EVENT_INFO, phase_velocities=PHASE_VELOCITIES,
        event_depth_km=EVENT_DEPTH_KM,
        fname="observed_waveform_gather_long_period")

    print("Generating spectrograms ...")
    fetch_plot_spectrograms(
        processed, EVENT_TIME, OUTDIR, **kw,
        pre_origin=PRE_ORIGIN, post_origin=POST_ORIGIN,
        fname="observed_spectrograms")

    print("Generating key station close-up ...")
    fetch_plot_closeup(
        processed, EVENT_TIME, OUTDIR,
        event_info=EVENT_INFO, key_stations=KEY_STATIONS,
        phase_velocities=PHASE_VELOCITIES,
        freqmin=FREQMIN, fmax=FREQMAX, post_origin=POST_ORIGIN,
        event_depth_km=EVENT_DEPTH_KM,
        fname="observed_closeup_ps23")

    print()
    print("=" * 72)
    print("  Done.  Figures saved to:", OUTDIR)
    print("=" * 72)


if __name__ == "__main__":
    main()
