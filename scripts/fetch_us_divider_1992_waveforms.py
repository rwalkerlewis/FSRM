#!/usr/bin/env python3
"""
Fetch and display real seismic waveforms from the US Divider nuclear test
(1992-09-23), the last US underground nuclear test.

Downloads broadband data from IRIS FDSN for CI TERRAscope, BK, NN/LN, and
IU stations, removes instrument response, bandpass-filters, and produces a
distance-sorted record section (waveform gather).

Event:  1992-09-23 15:04:00 UTC, NTS Area 4 (~37.021N, -116.058W)
        mb ~4.0, depth ~800 m, yield <20 kt (classified)
        Last US underground nuclear test before moratorium

NOTE:   1992 data availability may be limited.  Several stations were not yet
        deployed or only had short-period instruments.  The script includes
        extra error handling for missing / partial data.

Usage:
    python scripts/fetch_us_divider_1992_waveforms.py
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
EVENT_TIME = UTCDateTime("1992-09-23T15:04:00")
EVENT_LAT = 37.021
EVENT_LON = -116.058
EVENT_DEPTH_KM = 0.8

PRE_ORIGIN = 120
POST_ORIGIN = 600
FREQMIN = 0.5
FREQMAX = 8.0

OUTDIR = os.path.join(os.path.dirname(os.path.dirname(__file__)),
                      "figures", "us_divider_1992")

EVENT_INFO = {
    "name": "US Divider Test 1992-09-23",
    "date": "1992-09-23",
    "time": "15:04:00 UTC",
    "lat": EVENT_LAT,
    "lon": EVENT_LON,
    "depth_km": EVENT_DEPTH_KM,
    "mb": "4.0",
    "yield_text": "<20 kt (classified)",
}

# Basin & Range / NTS regional velocities (thinner crust than East Asia)
PHASE_VELOCITIES = {
    "Pn": 7.8,
    "Pg": 5.0,
    "Sn": 4.4,
    "Lg": 3.3,
}

KEY_STATIONS = ["PAS", "GSC", "SAO"]

STATIONS = [
    ("CI", "PAS", "*", "BH", "Pasadena, CA"),
    ("CI", "GSC", "*", "BH", "Goldstone, CA"),
    ("CI", "ISA", "*", "BH", "Isabella, CA"),
    ("CI", "PFO", "*", "BH", "Piñon Flat, CA"),
    ("CI", "SBC", "*", "BH", "Santa Barbara, CA"),
    ("CI", "SVD", "*", "BH", "Seven Oaks Dam, CA"),
    ("BK", "SAO", "*", "BH", "San Andreas Obs, CA"),
    ("NN", "ELK", "*", "BH", "Elko, NV"),
    ("LN", "MNV", "*", "BH", "Mina, NV"),
    ("IU", "ANMO", "*", "BH", "Albuquerque, NM"),
]

# Channel fallbacks for 1992 — some stations may only have VH or LH
_CHAN_FALLBACKS = ["BH", "LH", "VH", "SH"]


def main():
    os.makedirs(OUTDIR, exist_ok=True)

    print("=" * 72)
    print("  US Divider 1992-09-23 — Waveform Gather from FDSN/IRIS")
    print("=" * 72)
    print(f"  Origin time:  {EVENT_TIME}")
    print(f"  Location:     {EVENT_LAT:.3f}°N  {EVENT_LON:.3f}°W")
    print(f"  Depth:        ~{EVENT_DEPTH_KM} km (~800 m)")
    print(f"  Window:       {PRE_ORIGIN} s before → {POST_ORIGIN} s after origin")
    print(f"  Filter:       {FREQMIN}–{FREQMAX} Hz bandpass")
    print()

    # 1. Download  (with channel fallbacks for 1992 data)
    print("Downloading waveforms (with channel fallbacks for 1992 data) ...")
    raw_results = fetch_raw_waveforms(
        EVENT_TIME, EVENT_LAT, EVENT_LON, STATIONS,
        pre_origin=PRE_ORIGIN, post_origin=POST_ORIGIN,
        timeout=90, chan_fallbacks=_CHAN_FALLBACKS)
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
        fname="observed_closeup_key_stations")

    print()
    print("=" * 72)
    print("  Done.  Figures saved to:", OUTDIR)
    print("=" * 72)


if __name__ == "__main__":
    main()
