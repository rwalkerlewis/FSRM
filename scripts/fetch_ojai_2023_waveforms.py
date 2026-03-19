#!/usr/bin/env python3
"""
Fetch and display real seismic waveforms from the 2023 Ojai M5.1 earthquake.

Downloads broadband data from IRIS FDSN for regional stations in southern
California, removes instrument response, bandpass-filters, and produces
distance-sorted record sections and spectrograms.

Event:  2023-08-20 22:41:51 UTC, M5.1
        34.443°N  119.482°W, depth ~5.7 km
        Reverse / thrust — Western Transverse Ranges

Usage:
    python scripts/fetch_ojai_2023_waveforms.py
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
EVENT_TIME = UTCDateTime("2023-08-20T22:41:51")
EVENT_LAT = 34.443
EVENT_LON = -119.482
EVENT_DEPTH_KM = 5.7

PRE_ORIGIN = 60
POST_ORIGIN = 300
FREQMIN = 0.5
FREQMAX = 8.0

OUTDIR = os.path.join(os.path.dirname(os.path.dirname(__file__)),
                      "figures", "ojai_2023")

EVENT_INFO = {
    "name": "2023-08-20 Ojai M5.1 Earthquake",
    "date": "2023-08-20",
    "time": "22:41:51 UTC",
    "lat": EVENT_LAT,
    "lon": EVENT_LON,
    "depth_km": EVENT_DEPTH_KM,
    "mb": "5.1",
    "yield_text": "N/A (tectonic earthquake)",
}

PHASE_VELOCITIES = {
    "Pn": 7.8,
    "Pg": 6.3,
    "Sn": 4.5,
    "Lg": 3.5,
}

KEY_STATIONS = ["SBC", "PASC", "ISA"]

STATIONS = [
    ("CI", "SBC",  "*", "BH", "Santa Barbara, CA"),
    ("CI", "SNCC", "*", "BH", "Santa Monica, CA"),
    ("CI", "PASC", "*", "BH", "Pasadena, CA"),
    ("CI", "PAS",  "*", "BH", "Pasadena (CIT), CA"),
    ("CI", "OSI",  "*", "BH", "Osito Audit, CA"),
    ("CI", "ISA",  "*", "BH", "Isabella, CA"),
    ("CI", "TPNP", "*", "BH", "Tejon Pass, CA"),
    ("CI", "SLA",  "*", "BH", "Searles, CA"),
    ("CI", "GSC",  "*", "BH", "Goldstone, CA"),
    ("II", "PFO",  "*", "BH", "Piñon Flat, CA"),
    ("IU", "TUC",  "*", "BH", "Tucson, AZ"),
]


def main():
    os.makedirs(OUTDIR, exist_ok=True)

    print("=" * 72)
    print("  Ojai 2023-08-20 M5.1 Earthquake — Waveform Gather from FDSN/IRIS")
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

    # Common kwargs
    kw = dict(event_info=EVENT_INFO, phase_velocities=PHASE_VELOCITIES,
              freqmin=FREQMIN, fmax=FREQMAX)

    # 3. Record section — full window
    print("Generating waveform gather (record section – full window) ...")
    fetch_plot_record_section(
        processed, EVENT_TIME, OUTDIR, **kw,
        pre_origin=PRE_ORIGIN, post_origin=POST_ORIGIN,
        title_extra="instrument-response-corrected velocity",
        fname="observed_waveform_gather")

    # 4. Zoomed on expected signal
    print("Generating waveform gather (record section – zoomed) ...")
    fetch_plot_zoomed(
        processed, EVENT_TIME, OUTDIR, **kw,
        fname="observed_waveform_gather_zoomed")

    # 5. Distance-scaled
    print("Generating distance-scaled record section with phase arrivals ...")
    fetch_plot_distance_scaled(
        processed, EVENT_TIME, OUTDIR, **kw,
        fname="observed_waveform_gather_distance_scaled")

    # 6. Long-period
    print("Generating long-period filtered record section (0.02–0.10 Hz) ...")
    fetch_plot_long_period(
        raw_results, EVENT_TIME, OUTDIR,
        event_info=EVENT_INFO, phase_velocities=PHASE_VELOCITIES,
        fname="observed_waveform_gather_long_period")

    # 7. Spectrograms
    print("Generating spectrograms ...")
    fetch_plot_spectrograms(
        processed, EVENT_TIME, OUTDIR, **kw,
        pre_origin=PRE_ORIGIN, post_origin=POST_ORIGIN,
        fname="observed_spectrograms")

    # 8. Close-up on key stations
    print("Generating key station close-up ...")
    fetch_plot_closeup(
        processed, EVENT_TIME, OUTDIR,
        event_info=EVENT_INFO, key_stations=KEY_STATIONS,
        phase_velocities=PHASE_VELOCITIES,
        freqmin=FREQMIN, fmax=FREQMAX, post_origin=POST_ORIGIN,
        fname="observed_closeup_key_stations")

    print()
    print("=" * 72)
    print("  Done.  Figures saved to:", OUTDIR)
    print("=" * 72)


if __name__ == "__main__":
    main()
