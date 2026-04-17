#!/usr/bin/env python3
"""
Cannikin 1971: Observed vs. Mueller-Murphy Synthetic Comparison

Loads digitised waveforms from the KNDC CSS3.0 archive
(data/historic_nuclear/10.AMCHITKA/19711106.2200/) and compares them
against Mueller-Murphy synthetics computed for the Cannikin event.

The observed data are SAC files from Soviet-era stations in Central Asia
(teleseismic distance ~60-70 degrees). Instrument response is provided
as KNDC FAP (Frequency-Amplitude-Phase) files.

Produces a publication-quality multi-panel figure:
  - Panel 1: Observed record section (Z-component, all stations)
  - Panel 2: Mueller-Murphy synthetics at matching distances
  - Panel 3: Station-by-station overlay comparison
  - Panel 4: Spectral comparison

Usage:
    python scripts/compare_cannikin_1971.py
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

from obspy import read, UTCDateTime
from obspy.geodetics import gps2dist_azimuth

warnings.filterwarnings("ignore")

# Add scripts directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fsrm.signal_processing import time_relative, spectral_amplitude
from fsrm.source_physics import (
    cavity_radius_empirical, corner_frequency_patton,
    scalar_moment_coupled, mb_from_yield, mueller_murphy_spectrum,
)

# ============================================================================
# Event parameters — Cannikin, 1971-11-06 22:00:00 UTC
# ============================================================================
EVENT_TIME = UTCDateTime("1971-11-06T22:00:00")
EVENT_LAT = 51.4719
EVENT_LON = 179.1069
EVENT_DEPTH_KM = 1.86
YIELD_KT = 5000.0

# ============================================================================
# Paths
# ============================================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.dirname(SCRIPT_DIR)
DATA_DIR = os.path.join(REPO_DIR, "data", "historic_nuclear", "10.AMCHITKA")
WF_DIR = os.path.join(DATA_DIR, "19711106.2200", "wf")
RESP_DIR = os.path.join(DATA_DIR, "resp_kndc")
OUTDIR = os.path.join(REPO_DIR, "figures", "cannikin_1971")
SIM_OUTDIR = os.path.join(REPO_DIR, "output", "historical_cannikin_1971")

# ============================================================================
# Station metadata (from .site file)
# Coordinates: lat, lon, elev_km, description, instrument_type, resp_file
# ============================================================================
STATIONS = {
    "FRU":  (42.8411,  74.6128, 0.836, "Bishkek, Kyrgyzstan"),
    "ANVS": (42.7853,  77.6669, 1.860, "Ananevo, Kyrgyzstan"),
    "EKS":  (42.6683,  73.7864, 1.180, "Erkin-Sai, Kyrgyzstan"),
    "ARK":  (41.7995,  71.9667, 1.280, "Arkit, Kyrgyzstan"),
    "ARSB": (41.3233,  72.9820, 1.510, "Arslanbob, Kyrgyzstan"),
    "BRVK": (53.0578,  70.2828, 0.315, "Borovoye, N. Kazakhstan"),
}

# Channel-to-response mapping (from .sensor + .instrument files)
# LP channels use SKD, SP channels use SKM
RESPONSE_MAP = {
    "LP": "skd-1.ts25_tg1.2.fap",   # SKD long-period (most stations)
    "LP_ARK": "skd-2.ts20_tg1.2.fap",  # ARK uses SKD type 2
    "SP": "skm-4.ts2.0_tg0.32.fap",  # SKM short-period (BRVK, FRU SP)
}

# Filter band for comparison
FREQMIN_LP = 0.02   # long-period band
FREQMAX_LP = 0.5
FREQMIN_SP = 0.2    # short-period band
FREQMAX_SP = 5.0


def compute_distance(sta_lat, sta_lon):
    """Epicentral distance in km from event to station."""
    d_m, az, baz = gps2dist_azimuth(EVENT_LAT, EVENT_LON, sta_lat, sta_lon)
    return d_m / 1000.0, az, baz


def load_fap_response(resp_file):
    """
    Load a KNDC FAP (Frequency-Amplitude-Phase) instrument response file.

    Returns (freq, amplitude, phase) arrays. Phase is in radians.
    """
    filepath = os.path.join(RESP_DIR, resp_file)
    freqs, amps, phases = [], [], []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("theoretical"):
                continue
            # Skip the line count
            parts = line.split()
            if len(parts) == 1:
                continue
            if len(parts) >= 3:
                try:
                    freq = float(parts[0])
                    amp = float(parts[1])
                    phase = float(parts[2])
                    freqs.append(freq)
                    amps.append(amp)
                    phases.append(phase)
                except ValueError:
                    continue
    return np.array(freqs), np.array(amps), np.array(phases)


def deconvolve_instrument(data, dt, resp_file, water_level=0.01):
    """
    Remove instrument response using FAP file via spectral division.

    Parameters:
        data:        waveform data array
        dt:          sample interval (s)
        resp_file:   FAP response filename
        water_level: water-level stabilisation (fraction of max response)

    Returns deconvolved data (displacement proportional).
    """
    resp_f, resp_a, resp_p = load_fap_response(resp_file)

    n = len(data)
    freqs = np.fft.rfftfreq(n, d=dt)
    spec = np.fft.rfft(data)

    # Interpolate FAP response to FFT frequencies
    resp_amp_interp = np.interp(freqs, resp_f, resp_a, left=resp_a[0], right=resp_a[-1])
    resp_pha_interp = np.interp(freqs, resp_f, resp_p, left=resp_p[0], right=resp_p[-1])

    # Build complex response
    resp_complex = resp_amp_interp * np.exp(1j * 2 * np.pi * resp_pha_interp)

    # Water-level deconvolution
    wl = water_level * np.max(np.abs(resp_complex))
    denom = resp_complex.copy()
    small = np.abs(denom) < wl
    denom[small] = wl * np.exp(1j * np.angle(denom[small]))

    spec_decon = spec / denom
    return np.fft.irfft(spec_decon, n=n)


def load_observed_waveforms():
    """
    Load all Z-component SAC files for the Cannikin event.

    Returns list of (station, dist_km, az, trace, band, resp_file) sorted
    by distance, using only the first time segment per station/band.
    """
    results = []

    for sta_name, (lat, lon, elev, desc) in STATIONS.items():
        dist_km, az, baz = compute_distance(lat, lon)

        # Find Z-component files for this station (first segment = *z1.sac)
        pattern_lp = f"{sta_name}71310z1.sac"
        pattern_sp = f"{sta_name.lower()}71310z1.sac"  # BRVK SP files are lowercase

        # Try LP first
        lp_path = os.path.join(WF_DIR, pattern_lp)
        if os.path.exists(lp_path):
            st = read(lp_path)
            tr = st[0]
            resp = RESPONSE_MAP.get(f"LP_{sta_name}", RESPONSE_MAP["LP"])
            results.append((sta_name, dist_km, az, tr, "LP", resp, desc))
            print(f"  {sta_name:6s}  LP  dist={dist_km:7.1f} km  "
                  f"az={az:5.1f}°  npts={tr.stats.npts}  "
                  f"sr={tr.stats.sampling_rate} Hz")

        # Try SP (lowercase filenames for BRVK)
        sp_path = os.path.join(WF_DIR, pattern_sp)
        if not os.path.exists(sp_path):
            # Also try uppercase
            sp_path = os.path.join(WF_DIR, f"{sta_name}71310z1.sac")

        # Only add SP if it's a different instrument (BRVK is SP-only)
        if sta_name == "BRVK" and os.path.exists(
                os.path.join(WF_DIR, "brvk71310z1.sac")):
            sp_path = os.path.join(WF_DIR, "brvk71310z1.sac")
            st = read(sp_path)
            tr = st[0]
            results.append((sta_name, dist_km, az, tr, "SP",
                            RESPONSE_MAP["SP"], desc))
            print(f"  {sta_name:6s}  SP  dist={dist_km:7.1f} km  "
                  f"az={az:5.1f}°  npts={tr.stats.npts}  "
                  f"sr={tr.stats.sampling_rate} Hz")

    results.sort(key=lambda x: x[1])
    return results


def teleseismic_travel_time(dist_km, phase="P"):
    """
    Approximate teleseismic travel time for P and S phases.

    Uses IASP91-like approximations for distances 50-100 degrees.
    """
    dist_deg = dist_km / 111.195
    if phase == "P":
        # P travel time: ~10 s/deg at teleseismic distances
        return 10.0 * dist_deg + 20.0  # rough P arrival
    elif phase == "S":
        return 18.0 * dist_deg + 40.0  # rough S arrival
    elif phase == "PP":
        return 12.0 * dist_deg + 60.0
    elif phase == "SS":
        return 22.0 * dist_deg + 80.0
    return 0.0


def generate_teleseismic_synthetic(dist_km, yield_kt, dt=0.025, duration=600.0,
                                   fmin=0.02, fmax=0.5):
    """
    Generate a simplified Mueller-Murphy teleseismic P-wave synthetic.

    At teleseismic distances, uses P-wave with depth phases (pP, sP)
    and mantle attenuation.
    """
    npts = int(duration / dt)
    t = np.arange(npts) * dt
    vel = np.zeros(npts)

    fc = corner_frequency_patton(yield_kt)
    M0 = scalar_moment_coupled(yield_kt)

    f = np.fft.rfftfreq(npts, d=dt)
    f[0] = 1e-10
    omega = 2 * np.pi * f

    # Bandpass taper
    bp = np.ones_like(f)
    bp[f < fmin] = 0
    bp[f > fmax] = 0
    taper_lo = (f >= fmin * 0.5) & (f <= fmin)
    bp[taper_lo] = 0.5 * (1 - np.cos(np.pi * (f[taper_lo] - fmin * 0.5) / (fmin * 0.5)))
    taper_hi = (f >= fmax) & (f <= fmax * 1.25)
    bp[taper_hi] = 0.5 * (1 + np.cos(np.pi * (f[taper_hi] - fmax) / (fmax * 0.25)))

    # Teleseismic P arrival
    t_P = teleseismic_travel_time(dist_km, "P")
    dist_deg = dist_km / 111.195

    # Geometric spreading for teleseismic
    spreading = 1.0 / (dist_km * 1000.0)  # 1/r in metres

    # Mantle Q ~ 400 for P
    Q_mantle = 400.0
    t_star = t_P / Q_mantle  # t* attenuation operator

    # Source spectrum
    src = mueller_murphy_spectrum(f, yield_kt, fc=fc)
    att = np.exp(-np.pi * f * t_star)
    spec = src * att * spreading

    # Velocity spectrum with bandpass
    spec_vel = spec * omega * bp

    # Direct P
    phase_shift = np.exp(-1j * omega * t_P)
    p_wavelet = np.fft.irfft(spec_vel * phase_shift, n=npts)

    # pP depth phase (surface reflection, polarity reversal)
    depth_m = EVENT_DEPTH_KM * 1000.0
    vp_source = 4500.0  # m/s from config
    t_pP_delay = 2.0 * depth_m / vp_source
    pP_shift = np.exp(-1j * omega * (t_P + t_pP_delay))
    pP_wavelet = np.fft.irfft(spec_vel * pP_shift * (-0.9), n=npts)

    # sP depth phase
    vs_source = 2600.0  # m/s from config
    t_sP_delay = depth_m / vs_source + depth_m / vp_source
    sP_shift = np.exp(-1j * omega * (t_P + t_sP_delay))
    sP_wavelet = np.fft.irfft(spec_vel * sP_shift * 0.5, n=npts)

    vel = p_wavelet + pP_wavelet + sP_wavelet

    return t, vel, t_P


def process_and_filter(data, dt, fmin, fmax):
    """Detrend, taper, and bandpass-filter a waveform."""
    from scipy.signal import butter, filtfilt

    # Detrend
    d = data - np.mean(data)
    d = d - np.polyval(np.polyfit(np.arange(len(d)), d, 1), np.arange(len(d)))

    # Cosine taper (5%)
    ntaper = int(0.05 * len(d))
    taper = np.ones(len(d))
    taper[:ntaper] = 0.5 * (1 - np.cos(np.pi * np.arange(ntaper) / ntaper))
    taper[-ntaper:] = 0.5 * (1 - np.cos(np.pi * np.arange(ntaper) / ntaper))[::-1]
    d *= taper

    # Butterworth bandpass
    nyq = 0.5 / dt
    low = fmin / nyq
    high = min(fmax / nyq, 0.99)
    if low >= high or low <= 0:
        return d
    b, a = butter(3, [low, high], btype='band')
    d = filtfilt(b, a, d)
    return d


def plot_comparison(observed_data):
    """
    Produce the multi-panel comparison figure.

    Parameters:
        observed_data: list of (sta, dist_km, az, trace, band, resp_file, desc)
    """
    # Separate LP and SP stations
    lp_data = [(s, d, az, tr, b, r, desc)
                for s, d, az, tr, b, r, desc in observed_data if b == "LP"]
    sp_data = [(s, d, az, tr, b, r, desc)
                for s, d, az, tr, b, r, desc in observed_data if b == "SP"]

    n_lp = len(lp_data)
    if n_lp == 0:
        print("ERROR: No LP waveforms found.")
        return

    fig = plt.figure(figsize=(24, max(14, 2.0 * n_lp)))
    gs = GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.25,
                  height_ratios=[3, 2])

    ax_obs = fig.add_subplot(gs[0, 0])
    ax_syn = fig.add_subplot(gs[0, 1])
    ax_overlay = fig.add_subplot(gs[1, 0])
    ax_spectra = fig.add_subplot(gs[1, 1])

    # Colours for stations
    sta_colors = plt.cm.Set1(np.linspace(0, 1, n_lp))

    obs_traces = []  # (label, dist_km, times, data_filtered)
    syn_traces = []  # (label, dist_km, times, data)

    print("\n  Processing LP waveforms:")
    for i, (sta, dist_km, az, tr, band, resp_file, desc) in enumerate(lp_data):
        dt = tr.stats.delta
        data = tr.data.astype(float)

        # Deconvolve instrument response
        data_decon = deconvolve_instrument(data, dt, resp_file)

        # Filter to LP band
        data_filt = process_and_filter(data_decon, dt, FREQMIN_LP, FREQMAX_LP)

        # Time relative to event origin
        t_rel = np.arange(len(data_filt)) * dt + (tr.stats.starttime - EVENT_TIME)

        obs_traces.append((f"{sta} ({dist_km:.0f} km)", dist_km, t_rel, data_filt))
        print(f"    {sta:6s}  dist={dist_km:.0f} km  "
              f"t_range=[{t_rel[0]:.1f}, {t_rel[-1]:.1f}] s")

        # Generate matching synthetic
        t_syn, v_syn, t_P = generate_teleseismic_synthetic(
            dist_km, YIELD_KT, dt=0.025, duration=max(1200, t_rel[-1] + 60),
            fmin=FREQMIN_LP, fmax=FREQMAX_LP)

        syn_traces.append((f"{sta} synth", dist_km, t_syn, v_syn))

    # ── Panel 1: Observed record section ──
    trace_scale = 0.38
    y_positions = np.arange(n_lp)

    for i, (label, dist_km, t, d) in enumerate(obs_traces):
        peak = np.max(np.abs(d))
        if peak > 0:
            d_norm = d / peak * trace_scale
        else:
            d_norm = d

        ax_obs.fill_between(t, d_norm + y_positions[i], y_positions[i],
                            where=(d_norm > 0), color="black", alpha=0.25)
        ax_obs.fill_between(t, d_norm + y_positions[i], y_positions[i],
                            where=(d_norm < 0), color="gray", alpha=0.12)
        ax_obs.plot(t, d_norm + y_positions[i], "k-", lw=0.5, alpha=0.85)

        # P arrival marker
        t_P = teleseismic_travel_time(dist_km, "P")
        ax_obs.axvline(t_P, color="#1f77b4", ls=":", lw=0.8, alpha=0.5)

        ax_obs.text(-10, y_positions[i], label,
                    fontsize=8, fontweight="bold", ha="right", va="center",
                    color="#222", clip_on=False)

    ax_obs.set_ylim(-0.8, n_lp - 0.2)
    t_min_obs = min(t[0] for _, _, t, _ in obs_traces)
    t_max_obs = max(t[-1] for _, _, t, _ in obs_traces)
    ax_obs.set_xlim(max(t_min_obs, 0), min(t_max_obs, 1800))
    ax_obs.set_yticks(y_positions)
    ax_obs.set_yticklabels([])
    ax_obs.set_xlabel("Time after origin (s)", fontsize=10, fontweight="bold")
    ax_obs.set_title("Observed (KNDC Archive, LP Z-component)",
                     fontsize=12, fontweight="bold")
    ax_obs.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax_obs.grid(True, which="major", axis="x", alpha=0.2)

    # ── Panel 2: Synthetic record section ──
    for i, (label, dist_km, t, d) in enumerate(syn_traces):
        peak = np.max(np.abs(d))
        if peak > 0:
            d_norm = d / peak * trace_scale
        else:
            d_norm = d

        ax_syn.fill_between(t, d_norm + y_positions[i], y_positions[i],
                            where=(d_norm > 0), color="#1f77b4", alpha=0.25)
        ax_syn.fill_between(t, d_norm + y_positions[i], y_positions[i],
                            where=(d_norm < 0), color="#aec7e8", alpha=0.15)
        ax_syn.plot(t, d_norm + y_positions[i], color="#1f77b4",
                    lw=0.5, alpha=0.85)

        # P arrival marker
        t_P = teleseismic_travel_time(dist_km, "P")
        ax_syn.axvline(t_P, color="#1f77b4", ls=":", lw=0.8, alpha=0.5)

        ax_syn.text(-10, y_positions[i], label,
                    fontsize=8, fontweight="bold", ha="right", va="center",
                    color="#333", clip_on=False)

    ax_syn.set_ylim(-0.8, n_lp - 0.2)
    ax_syn.set_xlim(ax_obs.get_xlim())
    ax_syn.set_yticks(y_positions)
    ax_syn.set_yticklabels([])
    ax_syn.set_xlabel("Time after origin (s)", fontsize=10, fontweight="bold")
    ax_syn.set_title("Mueller-Murphy Teleseismic Synthetic",
                     fontsize=12, fontweight="bold")
    ax_syn.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax_syn.grid(True, which="major", axis="x", alpha=0.2)

    # ── Panel 3: Overlay comparison (selected stations) ──
    # Show first 3 stations overlaid
    n_overlay = min(3, n_lp)
    for i in range(n_overlay):
        o_label, o_dist, o_t, o_d = obs_traces[i]
        s_label, s_dist, s_t, s_d = syn_traces[i]

        # Normalise both
        o_peak = np.max(np.abs(o_d)) if np.max(np.abs(o_d)) > 0 else 1.0
        s_peak = np.max(np.abs(s_d)) if np.max(np.abs(s_d)) > 0 else 1.0

        offset = i * 2.5
        ax_overlay.plot(o_t, o_d / o_peak + offset, "k-", lw=0.7,
                        alpha=0.85, label=f"Observed {o_label}" if i == 0 else "")
        ax_overlay.plot(s_t, s_d / s_peak + offset, color="#1f77b4", lw=0.7,
                        alpha=0.75, label="Synthetic" if i == 0 else "")

        ax_overlay.text(ax_obs.get_xlim()[0] - 10, offset, o_label.split(" (")[0],
                        fontsize=8, fontweight="bold", ha="right", va="center",
                        color="#222")

    ax_overlay.set_xlabel("Time after origin (s)", fontsize=10, fontweight="bold")
    ax_overlay.set_title("Overlay: Observed (black) vs. Synthetic (blue)",
                         fontsize=12, fontweight="bold")
    ax_overlay.set_xlim(ax_obs.get_xlim())
    ax_overlay.set_yticks([])
    ax_overlay.legend(loc="upper right", fontsize=9)
    ax_overlay.grid(True, which="major", axis="x", alpha=0.2)

    # ── Panel 4: Spectral comparison ──
    for i, (o_data, s_data) in enumerate(zip(obs_traces, syn_traces)):
        o_label, o_dist, o_t, o_d = o_data
        s_label, s_dist, s_t, s_d = s_data

        sta_name = o_label.split(" (")[0]
        color = sta_colors[i]

        # Observed spectrum
        dt_obs = o_t[1] - o_t[0] if len(o_t) > 1 else 0.025
        f_obs, a_obs = spectral_amplitude(o_d, dt_obs)
        a_obs_norm = a_obs / np.max(a_obs) if np.max(a_obs) > 0 else a_obs
        ax_spectra.loglog(f_obs, a_obs_norm, color=color, ls="-", lw=1.2,
                          alpha=0.8, label=f"{sta_name} obs")

        # Synthetic spectrum
        dt_syn = s_t[1] - s_t[0] if len(s_t) > 1 else 0.025
        f_syn, a_syn = spectral_amplitude(s_d, dt_syn)
        a_syn_norm = a_syn / np.max(a_syn) if np.max(a_syn) > 0 else a_syn
        ax_spectra.loglog(f_syn, a_syn_norm, color=color, ls="--", lw=1.0,
                          alpha=0.6, label=f"{sta_name} syn")

    # Mueller-Murphy source spectrum reference
    f_ref = np.logspace(-2, 1, 500)
    fc = corner_frequency_patton(YIELD_KT)
    mm_spec = mueller_murphy_spectrum(f_ref, YIELD_KT, fc=fc)
    mm_norm = mm_spec / np.max(mm_spec)
    ax_spectra.loglog(f_ref, mm_norm, "k-", lw=2.0, alpha=0.4,
                      label=f"M-M source (fc={fc:.3f} Hz)")
    ax_spectra.axvline(fc, color="green", ls="--", lw=1.0, alpha=0.5,
                       label=f"fc = {fc:.3f} Hz")

    ax_spectra.set_xlabel("Frequency (Hz)", fontsize=10, fontweight="bold")
    ax_spectra.set_ylabel("Normalised Amplitude", fontsize=10, fontweight="bold")
    ax_spectra.set_title("Spectral Comparison", fontsize=12, fontweight="bold")
    ax_spectra.set_xlim(0.005, 5.0)
    ax_spectra.legend(loc="lower left", fontsize=7, ncol=2)
    ax_spectra.grid(True, which="both", alpha=0.15)

    # ── Suptitle with event info ──
    Rc = cavity_radius_empirical(YIELD_KT, rho=2400.0)
    mb = mb_from_yield(YIELD_KT)
    info_text = (
        f"Yield: ~{YIELD_KT:.0f} kt (~5 Mt)  |  DOB: {EVENT_DEPTH_KM*1000:.0f} m  |  "
        f"mb(pred): {mb:.2f}  mb(obs): 6.97  |  "
        f"Rc: {Rc:.0f} m  |  fc: {fc:.3f} Hz  |  "
        f"Filter: {FREQMIN_LP}-{FREQMAX_LP} Hz (LP)"
    )

    fig.suptitle(
        "Cannikin (1971-11-06) — Observed vs. Mueller-Murphy Synthetic\n"
        "Largest US underground nuclear test, Amchitka Island, ~5 Mt",
        fontsize=14, fontweight="bold", y=0.99)
    fig.text(0.5, 0.01, info_text, ha="center", va="bottom",
             fontsize=9, fontfamily="monospace",
             bbox=dict(boxstyle="round,pad=0.5", fc="#fffff0",
                       ec="#cca", alpha=0.9))

    plt.tight_layout(rect=[0, 0.04, 1, 0.95])
    outpath = os.path.join(OUTDIR, "cannikin_observed_vs_synthetic.png")
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"\n  Saved: {outpath}")
    return outpath


def plot_simulation_comparison(observed_data):
    """
    Compare observed waveforms against FSRM simulation output (if available).

    Looks for seismogram files in the simulation output directory.
    """
    if not os.path.isdir(SIM_OUTDIR):
        print(f"\n  Simulation output not found at: {SIM_OUTDIR}")
        print("  Run the simulation first, then re-run this script.")
        return

    sim_files = [f for f in os.listdir(SIM_OUTDIR)
                 if f.endswith(('.mseed', '.sac', '.csv', '.h5', '.hdf5'))]

    if not sim_files:
        print(f"\n  No seismogram files found in: {SIM_OUTDIR}")
        print("  The simulation may not have produced output yet.")
        print("  Configure [OUTPUT.SEISMOGRAMS] in the .config file and re-run.")
        return

    print(f"\n  Found {len(sim_files)} simulation output file(s):")
    for f in sorted(sim_files):
        print(f"    {f}")

    # Load simulation seismograms
    fig, axes = plt.subplots(len(sim_files), 1, figsize=(16, 3 * len(sim_files)),
                             squeeze=False)

    for i, sim_file in enumerate(sorted(sim_files)):
        filepath = os.path.join(SIM_OUTDIR, sim_file)
        ax = axes[i, 0]

        if sim_file.endswith('.sac') or sim_file.endswith('.mseed'):
            st = read(filepath)
            tr = st[0]
            t = np.arange(tr.stats.npts) * tr.stats.delta
            ax.plot(t, tr.data, "b-", lw=0.5)
            ax.set_title(f"Simulation: {sim_file} ({tr.stats.station})",
                         fontsize=10)
        elif sim_file.endswith('.csv'):
            data = np.loadtxt(filepath, delimiter=',', skiprows=1)
            ax.plot(data[:, 0], data[:, 1], "b-", lw=0.5)
            ax.set_title(f"Simulation: {sim_file}", fontsize=10)

        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Amplitude")
        ax.grid(True, alpha=0.2)

    plt.tight_layout()
    outpath = os.path.join(OUTDIR, "cannikin_simulation_output.png")
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {outpath}")


def print_station_summary(observed_data):
    """Print a summary table of stations and distances."""
    print(f"\n  {'Station':<8s} {'Band':<4s} {'Dist (km)':>10s} "
          f"{'Dist (deg)':>10s} {'Az':>6s} {'Description'}")
    print("  " + "-" * 70)
    for sta, dist_km, az, tr, band, resp, desc in observed_data:
        dist_deg = dist_km / 111.195
        print(f"  {sta:<8s} {band:<4s} {dist_km:10.1f} {dist_deg:10.1f} "
              f"{az:6.1f} {desc}")


def main():
    os.makedirs(OUTDIR, exist_ok=True)

    print("=" * 72)
    print("  Cannikin 1971: Observed vs. Synthetic Comparison")
    print("  Event: 1971-11-06 22:00:00 UTC, Amchitka Island, ~5 Mt")
    print("=" * 72)

    # Step 1: Load observed waveforms
    print("\n[1/4] Loading observed waveforms from KNDC archive ...")
    observed_data = load_observed_waveforms()
    if not observed_data:
        print("ERROR: No waveforms loaded. Check data paths.")
        sys.exit(1)
    print(f"       Loaded {len(observed_data)} trace(s)")

    # Station summary
    print_station_summary(observed_data)

    # Step 2: Generate comparison figure
    print("\n[2/4] Generating observed vs. Mueller-Murphy comparison ...")
    plot_comparison(observed_data)

    # Step 3: Check for simulation output
    print("\n[3/4] Checking for FSRM simulation output ...")
    plot_simulation_comparison(observed_data)

    # Step 4: Summary
    print("\n[4/4] Summary")
    print(f"       Output directory: {OUTDIR}")
    fc = corner_frequency_patton(YIELD_KT)
    Rc = cavity_radius_empirical(YIELD_KT, rho=2400.0)
    print(f"       Yield: ~{YIELD_KT:.0f} kt  |  fc: {fc:.3f} Hz  |  Rc: {Rc:.0f} m")
    print(f"       Stations: {len(observed_data)} traces from "
          f"{len(set(s for s,_,_,_,_,_,_ in observed_data))} stations")

    print("\n" + "=" * 72)
    print("  Done.")
    print("=" * 72)


if __name__ == "__main__":
    main()
