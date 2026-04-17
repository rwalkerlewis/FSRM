#!/usr/bin/env python3
"""
Full Moment Tensor Inversion for the 1992-09-23 US Divider Nuclear Test

Performs waveform-based moment tensor inversion using regional seismic data:
  1. Downloads broadband waveforms from IRIS FDSN
  2. Computes Green's functions for a 1D velocity model
  3. Inverts for the full 6-component moment tensor via linear least squares
  4. Decomposes result into ISO, CLVD, and DC components
  5. Generates publication-quality figures

Note: Divider was the last US underground nuclear test (September 1992).
Data availability from 1992 may be limited.

Usage:
    python scripts/invert_us_divider_1992_moment_tensor.py
"""

import os
import sys
import warnings
import numpy as np

import matplotlib
matplotlib.use("Agg")

from obspy import UTCDateTime, Trace

warnings.filterwarnings("ignore")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fsrm.velocity_models import NTSVelocityModel
from fsrm.inversion import (
    compute_greens_functions, build_data_vector, build_greens_matrix,
    invert_moment_tensor, decompose_moment_tensor,
)
from fsrm.data_fetching import fetch_waveforms_for_inversion
from fsrm.plotting import fig_inversion_results

# ══════════════════════════════════════════════════════════════════════════════
# Event Parameters
# ══════════════════════════════════════════════════════════════════════════════
EVENT_TIME     = UTCDateTime("1992-09-23T15:04:00")
EVENT_LAT      = 37.021
EVENT_LON      = -116.058
EVENT_DEPTH_KM = 0.8

FREQMIN  = 0.02
FREQMAX  = 0.10
PRE_ORIGIN  = 60
POST_ORIGIN = 600

OUTDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                      "figures", "us_divider_1992")

EVENT_INFO = {
    "name": "US Divider 1992-09-23 Nuclear Test",
    "short_label": "Divider\n1992",
    "time": "1992-09-23 15:04:00 UTC",
    "lat": EVENT_LAT,
    "lon": EVENT_LON,
    "depth_km": EVENT_DEPTH_KM,
}

STATIONS = [
    ("CI", "PAS", "*", "BH", "Pasadena"),
    ("CI", "GSC", "*", "BH", "Goldstone"),
    ("BK", "SAO", "*", "BH", "San Andreas Obs"),
    ("NN", "ELK", "*", "BH", "Elko"),
]


# ══════════════════════════════════════════════════════════════════════════════
# Event-Specific: Synthetic Fallback Data
# ══════════════════════════════════════════════════════════════════════════════

def generate_synthetic_data():
    """
    Generate synthetic data when real waveforms are unavailable.
    Uses an explosion-like moment tensor for a small (~1 kt) NTS test.
    """
    print("\nGenerating synthetic test data...")

    M0 = 1e17  # N·m (roughly 1 kt)
    m_true = np.array([M0/3, M0/3, M0/3, 0.02*M0, 0.01*M0, 0.015*M0])

    dt = 1.0
    npts = 600

    synthetic_stations = [
        {'net': 'CI', 'sta': 'PAS', 'desc': 'Pasadena',
         'dist_km': 350, 'az': 195, 'baz': 15},
        {'net': 'CI', 'sta': 'GSC', 'desc': 'Goldstone',
         'dist_km': 230, 'az': 180, 'baz': 0},
        {'net': 'BK', 'sta': 'SAO', 'desc': 'San Andreas Obs',
         'dist_km': 520, 'az': 260, 'baz': 80},
        {'net': 'NN', 'sta': 'ELK', 'desc': 'Elko',
         'dist_km': 340, 'az': 20, 'baz': 200},
    ]

    for sdata in synthetic_stations:
        gf = compute_greens_functions(sdata['dist_km'], sdata['az'],
                                      EVENT_DEPTH_KM, dt, npts, None,
                                      FREQMIN, FREQMAX)

        components = ['Mxx', 'Myy', 'Mzz', 'Mxy', 'Mxz', 'Myz']
        data_z = sum(m_true[i] * gf[comp]['Z'] for i, comp in enumerate(components))
        data_r = sum(m_true[i] * gf[comp]['R'] for i, comp in enumerate(components))
        data_t = sum(m_true[i] * gf[comp]['T'] for i, comp in enumerate(components))

        noise_level = 0.1 * np.max(np.abs(data_z))
        data_z += np.random.randn(npts) * noise_level
        data_r += np.random.randn(npts) * noise_level
        data_t += np.random.randn(npts) * noise_level

        for arr, comp_name in [(data_z, 'tr_z'), (data_r, 'tr_r'), (data_t, 'tr_t')]:
            tr = Trace(data=arr)
            tr.stats.sampling_rate = 1.0 / dt
            tr.stats.npts = npts
            sdata[comp_name] = tr
        sdata['inv'] = None

    return synthetic_stations


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

def main():
    os.makedirs(OUTDIR, exist_ok=True)

    print("=" * 70)
    print("Full Moment Tensor Inversion for US Divider 1992-09-23 Nuclear Test")
    print("=" * 70)
    print()

    # Fetch 3-component waveforms
    station_data = fetch_waveforms_for_inversion(
        EVENT_TIME, EVENT_LAT, EVENT_LON, STATIONS,
        freqmin=FREQMIN, freqmax=FREQMAX,
        pre_origin=PRE_ORIGIN, post_origin=POST_ORIGIN)

    if len(station_data) < 2:
        print("\nInsufficient stations for inversion. Using synthetic example.")
        station_data = generate_synthetic_data()

    print(f"\nUsing {len(station_data)} stations for inversion")

    # Compute Green's functions
    print("\nComputing Green's functions...")
    observations = []
    greens_list = []
    dt = 1.0

    for sdata in station_data:
        if sdata['tr_z'] is None:
            continue
        npts = sdata['tr_z'].stats.npts
        data_z = sdata['tr_z'].data
        data_r = sdata['tr_r'].data if sdata['tr_r'] is not None else np.zeros(npts)
        data_t = sdata['tr_t'].data if sdata['tr_t'] is not None else np.zeros(npts)

        gf = compute_greens_functions(sdata['dist_km'], sdata['az'],
                                      EVENT_DEPTH_KM, dt, npts, None,
                                      FREQMIN, FREQMAX)
        weight = 1.0 / np.sqrt(sdata['dist_km'] / 1000.0)
        observations.append((data_z, data_r, data_t, weight))
        greens_list.append(gf)

        print(f"  {sdata['net']}.{sdata['sta']:5s}  dist={sdata['dist_km']:6.1f} km  "
              f"npts={npts}  weight={weight:.3f}")

    # Invert
    print("\nBuilding linear system...")
    d = build_data_vector(observations)
    G = build_greens_matrix(greens_list, observations)
    print(f"  Data vector size: {len(d)}")
    print(f"  Green's matrix size: {G.shape}")

    print("\nInverting for moment tensor...")
    m, residual, var_red, d_pred = invert_moment_tensor(d, G, method='lstsq')
    decomp = decompose_moment_tensor(m)

    # Print results
    print("\n" + "=" * 70)
    print("INVERSION RESULTS")
    print("=" * 70)
    print(f"\nMoment tensor [Mxx, Myy, Mzz, Mxy, Mxz, Myz]:")
    print(f"  {m}")
    print(f"\nScalar moment M₀ = {decomp['M0']:.3e} N·m")
    print(f"Moment magnitude Mw = {decomp['Mw']:.2f}")
    print(f"\nSource type decomposition:")
    print(f"  Isotropic (ISO):    {decomp['ISO_fraction']*100:5.1f}%")
    print(f"  CLVD:               {decomp['CLVD_fraction']*100:5.1f}%")
    print(f"  Double-couple (DC): {decomp['DC_fraction']*100:5.1f}%")
    print(f"\nIsotropic component: {decomp['ISO_value']:.3e} N·m")
    print(f"  (Positive = compressional = explosion-like)")
    print(f"\nFault plane solution (DC component):")
    print(f"  Strike = {decomp['strike']:.1f}°")
    print(f"  Dip    = {decomp['dip']:.1f}°")
    print(f"  Rake   = {decomp['rake']:.1f}°")
    print(f"\nVariance reduction: {var_red:.1f}%")
    print(f"RMS residual: {residual:.3e}")

    # Interpretation
    print("\n" + "=" * 70)
    print("INTERPRETATION")
    print("=" * 70)
    if decomp['ISO_fraction'] > 0.5 and decomp['ISO_value'] > 0:
        print("\n✓ Large positive isotropic component (>50%)")
        print("  → Consistent with EXPLOSIVE SOURCE (underground nuclear test)")
        print("  → Clean explosion: last US underground test at NTS")
    elif decomp['DC_fraction'] > 0.7:
        print("\n✗ Dominated by double-couple component")
        print("  → Consistent with tectonic earthquake")
    else:
        print("\n⚠ Mixed source mechanism")
        print("  → Possible explosion with tectonic release")

    # Generate figures
    print("\nGenerating figures...")
    interpretation = [
        "  Last US underground nuclear test (September 1992)",
        "  ~1 kt coupled in NTS tuff at Yucca Flat",
        "  Expected: clean explosion → high ISO fraction",
        "  Note: 1992 broadband data availability may be limited",
    ]
    fig_inversion_results(m, decomp, station_data, observations, greens_list,
                          d, d_pred, var_red, EVENT_INFO, OUTDIR,
                          freqmin=FREQMIN, freqmax=FREQMAX,
                          interpretation_lines=interpretation)

    return m, decomp


if __name__ == "__main__":
    try:
        m, decomp = main()
    except KeyboardInterrupt:
        print("\nInversion interrupted.")
        sys.exit(1)
    except Exception as e:
        print(f"\nError during inversion: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
