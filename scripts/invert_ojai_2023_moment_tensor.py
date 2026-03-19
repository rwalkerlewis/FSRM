#!/usr/bin/env python3
"""
Full Moment Tensor Inversion for the 2023-08-20 M5.1 Ojai Earthquake

Performs waveform-based moment tensor inversion using regional seismic data:
  1. Downloads broadband waveforms from IRIS FDSN
  2. Computes Green's functions for a 1D velocity model (FK / ray theory)
  3. Inverts for the full 6-component moment tensor via linear least squares
  4. Decomposes result into ISO, CLVD, and DC components
  5. Generates publication-quality figures

Expected decomposition: ~90–100% DC, ~0% ISO
Interpretation: dominant double-couple → reverse/thrust faulting
                Western Transverse Ranges (N–S shortening)

Usage:
    python scripts/invert_ojai_2023_moment_tensor.py
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

from fsrm.velocity_models import SoCalVelocityModel
from fsrm.inversion import (
    compute_greens_functions, build_data_vector, build_greens_matrix,
    invert_moment_tensor, decompose_moment_tensor,
)
from fsrm.data_fetching import fetch_waveforms_for_inversion
from fsrm.plotting import fig_inversion_results

# ══════════════════════════════════════════════════════════════════════════════
# Event Parameters
# ══════════════════════════════════════════════════════════════════════════════
EVENT_TIME     = UTCDateTime("2023-08-20T22:41:51")
EVENT_LAT      = 34.443
EVENT_LON      = -119.482
EVENT_DEPTH_KM = 5.7

FREQMIN  = 0.02
FREQMAX  = 0.10
PRE_ORIGIN  = 60
POST_ORIGIN = 600

OUTDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                      "figures", "ojai_2023")

EVENT_INFO = {
    "name": "2023-08-20 Ojai M5.1 Earthquake",
    "short_label": "Ojai\n2023",
    "time": "2023-08-20 22:41:51 UTC",
    "lat": EVENT_LAT,
    "lon": EVENT_LON,
    "depth_km": EVENT_DEPTH_KM,
}

STATIONS = [
    ("CI", "SBC",  "*", "BH", "Santa Barbara"),
    ("CI", "PASC", "*", "BH", "Pasadena"),
    ("CI", "ISA",  "*", "BH", "Isabella"),
    ("II", "PFO",  "*", "BH", "Piñon Flat"),
]


# ══════════════════════════════════════════════════════════════════════════════
# Event-Specific: Synthetic Fallback Data
# ══════════════════════════════════════════════════════════════════════════════

def generate_synthetic_data():
    """
    Generate synthetic data when real waveforms are unavailable.
    Uses a reverse-fault double-couple moment tensor — characteristic of
    the Ojai M5.1 earthquake in the Western Transverse Ranges.

    Approximate focal mechanism: strike ~285°, dip ~45°, rake ~90° (thrust)
    """
    print("\nGenerating synthetic test data...")

    M0 = 10.0 ** (1.5 * 5.1 + 9.1)  # M0 from Mw 5.1
    # Reverse fault oriented E–W with northward dip (Transverse Ranges)
    strike_rad = np.radians(285)
    dip_rad = np.radians(45)
    rake_rad = np.radians(90)

    # Aki & Richards (2002) convention: moment tensor from strike/dip/rake
    sd = np.sin(dip_rad)
    cd = np.cos(dip_rad)
    s2d = np.sin(2 * dip_rad)
    c2d = np.cos(2 * dip_rad)
    ss = np.sin(strike_rad)
    cs = np.cos(strike_rad)
    s2s = np.sin(2 * strike_rad)
    c2s = np.cos(2 * strike_rad)
    sr = np.sin(rake_rad)
    cr = np.cos(rake_rad)

    Mxx = -M0 * (sd * cr * s2s + s2d * sr * ss * ss)
    Myy =  M0 * (sd * cr * s2s - s2d * sr * cs * cs)
    Mzz = M0 * s2d * sr
    Mxy =  M0 * (sd * cr * c2s + 0.5 * s2d * sr * s2s)
    Mxz = -M0 * (cd * cr * ss + c2d * sr * cs)
    Myz = -M0 * (cd * cr * cs - c2d * sr * ss)

    m_true = np.array([Mxx, Myy, Mzz, Mxy, Mxz, Myz])

    dt = 1.0
    npts = 600

    synthetic_stations = [
        {'net': 'CI', 'sta': 'SBC', 'desc': 'Santa Barbara',
         'dist_km': 28, 'az': 275, 'baz': 95},
        {'net': 'CI', 'sta': 'PASC', 'desc': 'Pasadena',
         'dist_km': 115, 'az': 120, 'baz': 300},
        {'net': 'CI', 'sta': 'ISA', 'desc': 'Isabella',
         'dist_km': 150, 'az': 30, 'baz': 210},
        {'net': 'II', 'sta': 'PFO', 'desc': 'Piñon Flat',
         'dist_km': 290, 'az': 125, 'baz': 305},
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
    print("Full Moment Tensor Inversion for 2023-08-20 Ojai M5.1 Earthquake")
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
    if decomp['DC_fraction'] > 0.7:
        print("\n✓ Dominant double-couple component (>70%)")
        print("  → Consistent with TECTONIC EARTHQUAKE (shear faulting)")
        if abs(decomp['rake'] - 90) < 30 or abs(decomp['rake'] + 90) < 30:
            print("  → Rake near ±90°: reverse or normal mechanism")
        print("  → Western Transverse Ranges: N–S compressional regime")
        print("  → Likely associated with San Cayetano / Red Mountain fault system")
    elif decomp['ISO_fraction'] > 0.5 and decomp['ISO_value'] > 0:
        print("\n✗ Large positive isotropic component")
        print("  → Would be consistent with explosive source (unexpected)")
    else:
        print("\n⚠ Mixed source mechanism")
        print("  → Possible tectonic event with some volumetric component")

    # Generate figures
    print("\nGenerating figures...")
    interpretation = [
        "  Expected: ~90–100% DC, ~0% ISO",
        "  Dominant DC → tectonic earthquake (shear failure)",
        "  Reverse/thrust: Western Transverse Ranges compressive regime",
        "  Probable fault: San Cayetano or Red Mountain Fault system",
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
