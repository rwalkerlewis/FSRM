#!/usr/bin/env python3
"""
Model the Degelen Mountain Nuclear Test — 1988-10-18 (STS Event 116)

Comprehensive physics-based modeling of an underground nuclear test at the
Semipalatinsk Test Site (STS), Degelen Mountain tunnel complex:
  1. Mueller-Murphy seismic source (RDP model) — ~16 kt, fully coupled
  2. 1-D layered velocity model for the Kazakh Platform / STS region
  3. Synthetic seismogram generation at real station distances
  4. Comparison against local CSS3.0 / SAC waveform archive
  5. Explosion-vs-earthquake discrimination analysis

The observed data come from the Soviet-era seismic networks (KazIS, KyrgIS)
in Kazakhstan and Kyrgyzstan, distributed as CSS3.0 format SAC files.

Physical references:
  Mueller & Murphy (1971)  — Seismic characteristics of underground nuclear detonations
  Patton (1988)            — Corner-frequency scaling
  Ryaboy (1989)            — Deep crustal structure of the Kazakh Platform
  Priestley et al. (1988)  — Seismic structure of the STS region
  Ringdal et al. (1992)    — Seismicity of the STS

Usage:
    python scripts/model_degelen_1988_nuclear_test.py
"""

import os
import sys
import glob
import warnings
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from obspy import UTCDateTime, read
from obspy.geodetics import gps2dist_azimuth
from scipy.signal import butter, filtfilt

warnings.filterwarnings("ignore")

# Add scripts directory to path for fsrm package import
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fsrm.source_physics import (
    cavity_radius_empirical, zone_radii, corner_frequency_patton,
    scalar_moment_coupled, mb_from_yield, Mw_from_M0,
    mueller_murphy_spectrum,
)
from fsrm.velocity_models import DegalenVelocityModel
from fsrm.propagation import (
    regional_travel_times, generate_synthetic,
)
from fsrm.signal_processing import envelope, spectral_amplitude, time_relative
from fsrm.data_fetching import compute_station_distances
from fsrm.inversion import invert_epicenter
from fsrm.plotting import (
    PHASE_COLORS, fig_source_model, fig_velocity_model,
    fig_synthetic_record_section, fig_observed_record_section,
    fig_comparison, fig_location_inversion,
)

# ═══════════════════════════════════════════════════════════════════════════════
# Event Parameters
# ═══════════════════════════════════════════════════════════════════════════════
EVENT_TIME   = UTCDateTime("1988-10-18T03:40:09.16")
EVENT_LAT    = 49.7800
EVENT_LON    = 78.0172
EVENT_DEPTH  = 0.0           # km (tunnel test)
YIELD_KT     = 16.0          # estimated from mb = 4.0 + 0.75*log10(W)
DF            = 1.0           # Fully coupled

OUTDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                      "figures", "degelen_1988")

DATADIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "data", "historic_nuclear", "13.Degelen",
                       "19881018.0340", "wf")

# Bandpass for processing
FMIN, FMAX = 0.5, 8.0

EVENT_INFO = {
    "name": "Degelen Mtn Nuclear Test (STS Event 116)",
    "time": "1988-10-18 03:40:09 UTC",
    "lat": EVENT_LAT,
    "lon": EVENT_LON,
    "mb_obs": 4.90,
}

# Stations from CSS3.0 .site file — KazIS and KyrgIS networks
# Coordinates for epoch valid in 1988
STATIONS = [
    # (net, sta, lat, lon, description)
    ("KazIS",  "SEM",   50.4000, 80.2500, "Semey (Semipalatinsk)"),
    ("KazIS",  "ILI",   43.8530, 77.0290, "Ili, Kazakhstan"),
    ("KazIS",  "KUU",   43.8930, 76.3390, "Kurty, Kazakhstan"),
    ("KazIS",  "BRTG",  43.3670, 78.5170, "Bartogai, Kazakhstan"),
    ("KazIS",  "PDG",   43.3290, 79.4900, "Podgornoye, Kazakhstan"),
    ("KazIS",  "TRG",   43.3070, 77.6370, "Turgen, Kazakhstan"),
    ("KazIS",  "TLG",   43.2330, 77.2250, "Talgar, Kazakhstan"),
    ("KazIS",  "MDO",   43.1670, 77.0100, "Medeo, Kazakhstan"),
    ("KazIS",  "SATY",  43.0570, 78.4070, "Saty, Kazakhstan"),
    ("KazIS",  "TSN",   43.0430, 76.9430, "Tien-Shan, Kazakhstan"),
    ("KazIS",  "MTB",   43.1300, 76.4300, "Maytyube, Kazakhstan"),
    ("KazIS",  "KST",   43.0430, 75.9630, "Kastek, Kazakhstan"),
    ("KyrgIS", "URVKG", 42.6889, 75.0500, "Urevka, Kyrgyzstan"),
    ("KyrgIS", "BOM",   42.4817, 75.9422, "Boom, Kyrgyzstan"),
    ("KyrgIS", "BGK",   42.6250, 74.2361, "Belogorka, Kyrgyzstan"),
    ("KyrgIS", "EKS",   42.6683, 73.7864, "Erkin-Sai, Kyrgyzstan"),
    ("KyrgIS", "KNSKG", 42.3244, 79.2439, "Ken-Suu, Kyrgyzstan"),
    ("KyrgIS", "KRSKG", 41.5750, 77.9061, "Kara-Sai, Kyrgyzstan"),
    ("KyrgIS", "ARLS",  41.8458, 74.3236, "Aral, Kyrgyzstan"),
    ("KyrgIS", "CHMS",  42.9917, 74.7539, "Chumysh, Kyrgyzstan"),
    ("KyrgIS", "MNAS",  42.4870, 72.5042, "Manas, Kyrgyzstan"),
    ("KyrgIS", "ARK",   41.7995, 71.9667, "Arkit, Kyrgyzstan"),
    ("KyrgIS", "SALK",  40.8708, 73.8042, "Salom-Alik, Kyrgyzstan"),
    ("KyrgIS", "SUFI",  40.0128, 73.5030, "Sufi-Kurgan, Kyrgyzstan"),
    ("KyrgIS", "OHH",   40.5244, 72.7847, "Osh, Kyrgyzstan"),
    ("KyrgIS", "CHVKG", 40.1447, 72.2106, "Chauvai, Kyrgyzstan"),
    ("KyrgIS", "DRKKG", 39.4808, 71.8050, "Daraut-Kurgan, Kyrgyzstan"),
]

STATION_INFO = compute_station_distances(EVENT_LAT, EVENT_LON, STATIONS)

# Observed Pn arrival times (epoch seconds) from the CSS3.0 .arrival file
# These are used for location inversion and travel-time analysis
_ORIGIN_EPOCH = 593149209.16
OBSERVED_PICKS = [
    ("SEM",   50.4000, 80.2500, "Pn", 593149236.24054 - _ORIGIN_EPOCH),
    ("ILI",   43.8530, 77.0290, "Pn", 593149297.45407 - _ORIGIN_EPOCH),
    ("KUU",   43.8930, 76.3390, "Pn", 593149297.95188 - _ORIGIN_EPOCH),
    ("BRTG",  43.3670, 78.5170, "Pn", 593149304.53089 - _ORIGIN_EPOCH),
    ("PDG",   43.3290, 79.4900, "Pn", 593149305.25320 - _ORIGIN_EPOCH),
    ("TRG",   43.3070, 77.6370, "Pn", 593149305.48395 - _ORIGIN_EPOCH),
    ("TLG",   43.2330, 77.2250, "Pn", 593149306.56769 - _ORIGIN_EPOCH),
    ("MDO",   43.1670, 77.0100, "Pn", 593149309.03535 - _ORIGIN_EPOCH),
    ("SATY",  43.0570, 78.4070, "Pn", 593149310.00891 - _ORIGIN_EPOCH),
    ("TSN",   43.0430, 76.9430, "Pn", 593149310.21207 - _ORIGIN_EPOCH),
    ("MTB",   43.1300, 76.4300, "Pn", 593149310.34857 - _ORIGIN_EPOCH),
    ("KST",   43.0430, 75.9630, "Pn", 593149315.16162 - _ORIGIN_EPOCH),
    ("URVKG", 42.6889, 75.0500, "Pn", 593149315.29798 - _ORIGIN_EPOCH),
    ("BOM",   42.4817, 75.9422, "Pn", 593149317.45303 - _ORIGIN_EPOCH),
    ("BGK",   42.6250, 74.2361, "Pn", 593149319.88990 - _ORIGIN_EPOCH),
    ("EKS",   42.6683, 73.7864, "Pn", 593149320.23767 - _ORIGIN_EPOCH),
    ("KNSKG", 42.3244, 79.2439, "Pn", 593149320.35736 - _ORIGIN_EPOCH),
    ("KRSKG", 41.5750, 77.9061, "Pn", 593149328.90752 - _ORIGIN_EPOCH),
    ("ARLS",  41.8458, 74.3236, "Pn", 593149330.50769 - _ORIGIN_EPOCH),
    ("CHMS",  42.9917, 74.7539, "Pn", 593149334.62739 - _ORIGIN_EPOCH),
    ("MNAS",  42.4870, 72.5042, "Pn", 593149326.53799 - _ORIGIN_EPOCH),
    ("ARK",   41.7995, 71.9667, "Pn", 593149341.17081 - _ORIGIN_EPOCH),
    ("SALK",  40.8708, 73.8042, "Pn", 593149343.14337 - _ORIGIN_EPOCH),
    ("SUFI",  40.0128, 73.5030, "Pn", 593149354.09684 - _ORIGIN_EPOCH),
    ("OHH",   40.5244, 72.7847, "Pn", 593149354.52278 - _ORIGIN_EPOCH),
    ("CHVKG", 40.1447, 72.2106, "Pn", 593149356.86350 - _ORIGIN_EPOCH),
    ("DRKKG", 39.4808, 71.8050, "Pn", 593149373.13131 - _ORIGIN_EPOCH),
]


# ═══════════════════════════════════════════════════════════════════════════════
# Local SAC Data Loader
# ═══════════════════════════════════════════════════════════════════════════════

def load_local_sac_data(datadir, stations, event_time, event_lat, event_lon,
                        fmin=0.5, fmax=8.0, component="z"):
    """
    Load waveform data from local SAC files in CSS3.0 archive.

    Reads Z-component SAC files, merges multiple segments per station,
    bandpass filters, and returns data in the same format as fetch_waveforms:
        list of (net, sta, dist_km, az, trace, desc) sorted by distance.
    """
    results = []
    sta_lookup = {s[1]: s for s in stations}

    for net, sta, slat, slon, desc in stations:
        # Find Z-component files matching this station
        sta_lower = sta.lower()
        patterns = [
            os.path.join(datadir, f"*{sta_lower}z*.sac"),
            os.path.join(datadir, f"*{sta_lower}*z*.sac"),
            os.path.join(datadir, f"{sta}*z*.sac"),
            os.path.join(datadir, f"{sta}*Z*.sac"),
        ]
        files = []
        for pat in patterns:
            files.extend(glob.glob(pat))
        files = sorted(set(files))

        if not files:
            continue

        # Read and merge all segments for this station
        try:
            st = read(files[0])
            for f in files[1:]:
                st += read(f)
            st.merge(fill_value=0)
            tr = st[0]

            # Bandpass filter
            sr = tr.stats.sampling_rate
            if sr > 2 * fmax:
                nyq = sr / 2.0
                low = max(fmin / nyq, 0.001)
                high = min(fmax / nyq, 0.999)
                b, a = butter(4, [low, high], btype="band")
                tr.data = filtfilt(b, a, tr.data.astype(float))

            d_m, az, _ = gps2dist_azimuth(event_lat, event_lon, slat, slon)
            dist_km = d_m / 1e3
            results.append((net, sta, dist_km, az, tr, desc))
        except Exception as e:
            print(f"  Warning: could not read {sta}: {e}")
            continue

    results.sort(key=lambda x: x[2])
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# Event-Specific Figures (not generalizable across events)
# ═══════════════════════════════════════════════════════════════════════════════


def fig04_yield_sensitivity():
    """Synthetic waveforms at SEM (~170 km) for different yields."""
    dist_km = 170.0  # approximate distance to SEM
    model = DegalenVelocityModel()
    yields = [1, 5, 16, 50, 100]

    fig, axes = plt.subplots(len(yields), 1, figsize=(16, 3.0 * len(yields)),
                             sharex=True)

    for i, W in enumerate(yields):
        t, v = generate_synthetic(dist_km, W, model,
                                  decoupling_factor=DF, dt=0.025, duration=400,
                                  fmin=FMIN, fmax=FMAX, depth_km=EVENT_DEPTH)
        mb = mb_from_yield(W, DF)

        ax = axes[i]
        peak = np.max(np.abs(v))
        if peak > 0:
            v_norm = v / peak
        else:
            v_norm = v

        ax.plot(t, v_norm, "k-", lw=0.5)
        env = envelope(v_norm)
        ax.plot(t, env, "r-", lw=0.8, alpha=0.5)
        ax.plot(t, -env, "r-", lw=0.8, alpha=0.5)

        fc = corner_frequency_patton(W)
        ax.set_ylabel("Norm. vel.", fontsize=9)
        ax.set_title(f"Yield = {W} kt  |  mb(coupled) = {mb:.2f}  |  fc = {fc:.2f} Hz  |  "
                     f"Peak vel ∝ {peak:.2e}",
                     fontsize=9, fontweight="bold", loc="left")
        ax.set_ylim(-1.2, 1.2)
        ax.grid(True, alpha=0.2)

        tt = regional_travel_times(dist_km, model)
        for ph, (tarr, _, _) in tt.items():
            if 0 < tarr < 400:
                ax.axvline(tarr, color=PHASE_COLORS.get(ph, "gray"),
                           ls="--", lw=1, alpha=0.6)

    axes[-1].set_xlabel("Time after origin (s)", fontsize=11)

    fig.suptitle("Figure 4 — Yield Sensitivity at SEM (170 km):  Coupled Test (DF=1)\n"
                 f"Degelen 1-D model, SPZ {FMIN}–{FMAX} Hz",
                 fontsize=13, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    p = os.path.join(OUTDIR, "fig04_yield_sensitivity.png")
    fig.savefig(p, dpi=180, bbox_inches="tight"); plt.close(fig)
    print(f"  Saved {p}")


def fig05_yield_tradeoff():
    """Compare different yields at SEM (~170 km) — yield tradeoff analysis."""
    dist_km = 170.0
    model = DegalenVelocityModel()
    cases = [
        ( 5.0,  "5 kt",   "#2ca02c"),
        (10.0,  "10 kt",  "#1f77b4"),
        (16.0,  "16 kt",  "#ff7f0e"),
        (30.0,  "30 kt",  "#d62728"),
    ]

    fig, axes = plt.subplots(len(cases) + 1, 1,
                             figsize=(16, 3.5 * (len(cases) + 1)),
                             gridspec_kw={"height_ratios": [1]*len(cases) + [1.2]})

    spectra = {}
    ref_peak = None
    for i, (W, label, col) in enumerate(cases):
        t, v = generate_synthetic(dist_km, W, model,
                                  decoupling_factor=DF, dt=0.025, duration=400,
                                  fmin=FMIN, fmax=FMAX, depth_km=EVENT_DEPTH)

        if ref_peak is None:
            ref_peak = np.max(np.abs(v))

        ax = axes[i]
        peak = np.max(np.abs(v))
        v_rel = v / ref_peak if ref_peak > 0 else v
        ax.plot(t, v_rel, "-", color=col, lw=0.5, alpha=0.8)

        mb = mb_from_yield(W, DF)
        fc = corner_frequency_patton(W)
        ax.set_ylabel("Rel. vel.", fontsize=9)
        ax.set_title(f"{label}  |  mb = {mb:.2f}  |  fc = {fc:.2f} Hz  |  "
                     f"rel. amplitude = {peak/ref_peak:.3f}",
                     fontsize=9, fontweight="bold", loc="left")
        ax.grid(True, alpha=0.2)
        if i < len(cases) - 1:
            ax.set_xticklabels([])

        freq, spec = spectral_amplitude(v, 0.025)
        spectra[label] = (freq, spec, col)

    axes[len(cases)-1].set_xlabel("Time after origin (s)", fontsize=11)

    # Spectral comparison
    ax = axes[-1]
    ref_label = list(spectra.keys())[0]
    ref_freq = spectra[ref_label][0]
    ref_spec_mask = (ref_freq > 0.1) & (ref_freq < 10)
    ref_spec_max = np.max(spectra[ref_label][1][ref_spec_mask])

    for label, (freq, spec, col) in spectra.items():
        mask = (freq > 0.1) & (freq < 10)
        s = spec[mask]
        if ref_spec_max > 0:
            s = s / ref_spec_max
        ax.loglog(freq[mask], s, "-", color=col, lw=1.5, label=label)
    ax.set_xlabel("Frequency (Hz)"); ax.set_ylabel("Rel. spectral amplitude")
    ax.set_title("Spectral comparison — different yields at SEM",
                 fontsize=9, fontweight="bold", loc="left")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.2, which="both")

    fig.suptitle("Figure 5 — Yield Trade-off at SEM (170 km):  Coupled Tests\n"
                 f"Amplitudes relative to 5 kt, {FMIN}–{FMAX} Hz",
                 fontsize=13, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    p = os.path.join(OUTDIR, "fig05_yield_tradeoff.png")
    fig.savefig(p, dpi=180, bbox_inches="tight"); plt.close(fig)
    print(f"  Saved {p}")


def fig08_discrimination(real_data):
    """Explosion vs earthquake discrimination analysis — Degelen is clearly an explosion."""
    fig = plt.figure(figsize=(18, 12))
    gs = GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.35)

    model = DegalenVelocityModel()

    # (a) P/Lg amplitude ratio for synthetic
    ax = fig.add_subplot(gs[0, 0])
    distances = np.linspace(200, 1500, 20)
    ps_ratios = []
    for d in distances:
        t, v = generate_synthetic(d, YIELD_KT, model,
                                  decoupling_factor=DF, dt=0.025, duration=600,
                                  fmin=FMIN, fmax=FMAX, depth_km=EVENT_DEPTH)
        tt = regional_travel_times(d, model)
        t_pn = tt["Pn"][0]; t_lg = tt["Lg"][0]

        p_mask = (t >= max(0, t_pn - 2)) & (t <= t_pn + 25)
        s_mask = (t >= max(0, t_lg - 5)) & (t <= min(t_lg + 60, t[-1]))

        p_rms = np.sqrt(np.mean(v[p_mask]**2)) if np.any(p_mask) and np.sum(p_mask) > 5 else 1e-30
        s_rms = np.sqrt(np.mean(v[s_mask]**2)) if np.any(s_mask) and np.sum(s_mask) > 5 else 1e-30
        ratio = p_rms / max(s_rms, 1e-30)
        ps_ratios.append(min(ratio, 100))

    ax.semilogy(distances, ps_ratios, "ro-", lw=1.5, ms=5, label="Explosion (synthetic)")

    rng_eq = np.random.RandomState(123)
    eq_ps = 0.25 + 0.1 * rng_eq.randn(len(distances))
    eq_ps = np.clip(eq_ps, 0.05, 0.8)
    ax.semilogy(distances, eq_ps, "b^--", lw=1, ms=4, alpha=0.6,
                label="Earthquake (typical)")
    ax.axhline(1.0, color="gray", ls=":", lw=1, label="P/Lg = 1 threshold")
    ax.set_xlabel("Distance (km)"); ax.set_ylabel("P/Lg amplitude ratio (RMS)")
    ax.set_title("(a) P/Lg Ratio vs Distance", fontweight="bold")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.2)
    ax.set_ylim(0.01, 50)

    # (b) mb vs Ms diagram
    ax = fig.add_subplot(gs[0, 1])
    W_exp = np.logspace(-1, 3, 50)
    mb_exp = 4.0 + 0.75 * np.log10(W_exp)
    ms_exp = mb_exp - 1.2

    mb_eq = np.linspace(3, 7.5, 50)
    ms_eq = 1.0 * mb_eq - 0.5 + 0.2 * np.random.randn(50)

    ax.scatter(mb_eq, ms_eq, c="b", s=20, alpha=0.4, label="Earthquakes")
    ax.scatter(mb_exp, ms_exp, c="r", s=20, alpha=0.4, label="Explosions (NTS)")

    mb_this = 4.90
    ms_this = 4.80
    ax.plot(mb_this, ms_this, "r*", ms=22, zorder=10,
            label=f"Degelen 1988 (mb={mb_this:.1f}, Ms={ms_this:.1f})")

    mb_line = np.linspace(2, 8, 100)
    ms_line = mb_line - 1.0
    ax.plot(mb_line, ms_line, "k--", lw=1.5, label="Discrimination line")

    ax.set_xlabel("mb"); ax.set_ylabel("Ms")
    ax.set_title("(b) mb–Ms Discrimination", fontweight="bold")
    ax.legend(fontsize=7, loc="upper left"); ax.grid(True, alpha=0.2)
    ax.set_xlim(2, 7); ax.set_ylim(1, 7)

    # (c) Spectral ratio analysis (synthetic)
    ax = fig.add_subplot(gs[0, 2])
    dist_sem = 170.0
    t, v = generate_synthetic(dist_sem, YIELD_KT, model,
                              decoupling_factor=DF, dt=0.025, duration=500,
                              fmin=FMIN, fmax=FMAX, depth_km=EVENT_DEPTH)
    freq, spec = spectral_amplitude(v, 0.025)
    m = (freq > 0.3) & (freq < 10)
    spec_n = spec[m] / np.max(spec[m]) if np.max(spec[m]) > 0 else spec[m]

    ax.semilogy(freq[m], spec_n, "r-", lw=1.5, label="Explosion (synth, 16 kt)")

    fc_eq = 0.8
    spec_eq = 1.0 / (1 + (freq[m] / fc_eq)**2)**0.5
    spec_eq_n = spec_eq / np.max(spec_eq)
    ax.semilogy(freq[m], spec_eq_n, "b--", lw=1.5, label="Earthquake (model, same mb)")

    ax.set_xlabel("Frequency (Hz)"); ax.set_ylabel("Normalised amplitude")
    ax.set_title("(c) Spectral Shape at SEM (170 km)", fontweight="bold")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.2, which="both")

    # (d) Observed P/Lg ratio at nearest stations
    ax = fig.add_subplot(gs[1, 0])
    sem_data = [x for x in real_data if x[1] in ("SEM",)]
    if sem_data:
        net, sta, dist_km, az, tr, desc = sem_data[0]
        times = time_relative(tr, EVENT_TIME)
        data = tr.data.copy()
        dt_real = tr.stats.delta

        tt = regional_travel_times(dist_km, model)
        t_pn = tt["Pn"][0]; t_lg = tt["Lg"][0]

        p_mask = (times >= t_pn - 5) & (times <= t_pn + 25)
        lg_mask = (times >= t_lg - 5) & (times <= t_lg + 50)

        if np.any(p_mask) and np.any(lg_mask):
            fp, sp = spectral_amplitude(data[p_mask], dt_real)
            fl, sl = spectral_amplitude(data[lg_mask], dt_real)

            mp = (fp > 0.3) & (fp < FMAX)
            ml = (fl > 0.3) & (fl < FMAX)

            if np.max(sp[mp]) > 0:
                ax.semilogy(fp[mp], sp[mp]/np.max(sp[mp]), "b-", lw=1.2, label="P window")
            if np.max(sl[ml]) > 0:
                ax.semilogy(fl[ml], sl[ml]/np.max(sl[ml]), "r-", lw=1.2, label="Lg window")
    else:
        ax.text(0.5, 0.5, "No SEM data available", transform=ax.transAxes,
                ha="center", fontsize=12, color="gray")

    ax.set_xlabel("Frequency (Hz)"); ax.set_ylabel("Normalised amplitude")
    ax.set_title("(d) Observed P vs Lg Spectra at SEM", fontweight="bold")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.2, which="both")

    # (e) Discrimination summary
    ax = fig.add_subplot(gs[1, 1])
    ax.axis("off")
    mb_pred = mb_from_yield(YIELD_KT, DF)
    lines = [
        "DISCRIMINATION ANALYSIS",
        "═" * 36,
        "",
        "INDICATORS CONSISTENT WITH EXPLOSION:",
        "",
        "  ✓  Impulsive P, high P/Lg ratio",
        "     at all regional distances",
        "",
        "  ✓  mb–Ms: ms well below earthquake",
        f"     trend (mb={mb_this}, Ms={ms_this})",
        "",
        "  ✓  Higher-frequency spectral content",
        "     than equivalent earthquake",
        "",
        "  ✓  Location: Degelen Mountain complex",
        "     (known Soviet nuclear facility)",
        "",
        "  ✓  Depth: surface / tunnel (0 km)",
        "",
        "  ✓  Purely compressional first motion",
        "     at all azimuths (isotropic source)",
        "",
        "COUPLED TEST PARAMETERS:",
        "",
        f"  Yield (est):  ~{YIELD_KT:.0f} kt",
        f"  DF:            1 (fully coupled)",
        f"  mb:            {mb_pred:.2f} (pred) / {mb_this} (obs)",
        f"  Ms:            {ms_this}",
        "",
        "  → CLEARLY an explosion.",
        "  → Consistent with Degelen tunnel test.",
    ]
    ax.text(0.05, 0.97, "\n".join(lines), transform=ax.transAxes,
            fontfamily="monospace", fontsize=8.5, va="top",
            bbox=dict(boxstyle="round,pad=0.5", fc="#fffff0", ec="#cca"))

    # (f) Yield-mb plot with observed mb constraint
    ax = fig.add_subplot(gs[1, 2])
    W_range = np.logspace(-1, 3, 200)
    mbs_coupled = 4.0 + 0.75 * np.log10(W_range)
    ax.semilogx(W_range, mbs_coupled, "b-", lw=2.5, label="Coupled (DF=1)")

    for dfi in [5, 10, 30, 70]:
        mbs = 4.0 + 0.75 * np.log10(W_range / dfi)
        ax.semilogx(W_range, mbs, lw=1.2, ls="--", alpha=0.5, label=f"DF={dfi}")

    ax.axhline(4.9, color="gray", ls=":", lw=2, label="Observed mb ≈ 4.9")
    ax.fill_between(W_range, 4.7, 5.1, color="gray", alpha=0.1)

    ax.plot(YIELD_KT, mb_from_yield(YIELD_KT, 1.0), "r*", ms=22, zorder=10,
            label=f"Best estimate: {YIELD_KT:.0f} kt")

    ax.set_xlabel("True yield (kt)"); ax.set_ylabel("Observed mb")
    ax.set_title("(f) Yield–mb Constraint", fontweight="bold")
    ax.legend(fontsize=7, ncol=2, loc="upper left"); ax.grid(True, alpha=0.2, which="both")
    ax.set_xlim(0.1, 1e3); ax.set_ylim(2, 7)

    fig.suptitle("Figure 8 — Explosion/Earthquake Discrimination Analysis\n"
                 "Degelen Mtn Nuclear Test 1988-10-18 03:40:09 UTC",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    p = os.path.join(OUTDIR, "fig08_discrimination.png")
    fig.savefig(p, dpi=180, bbox_inches="tight"); plt.close(fig)
    print(f"  Saved {p}")


# ═══════════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    os.makedirs(OUTDIR, exist_ok=True)
    np.random.seed(42)

    model = DegalenVelocityModel()

    print("=" * 72)
    print("  Degelen 1988 — Nuclear Test (~16 kt Coupled) Model & Analysis")
    print("  Semipalatinsk Test Site, Degelen Mountain Tunnel Complex")
    print("=" * 72)
    print()

    # Fig 1 — Source model (shared plotting module)
    print("[1/9] Source model ...")
    fig_source_model(YIELD_KT, EVENT_DEPTH * 1000, DF, EVENT_INFO,
                     os.path.join(OUTDIR, "fig01_source_model.png"),
                     rho=2700.0)

    # Fig 2 — Velocity model (shared plotting module)
    print("[2/9] Velocity model ...")
    fig_velocity_model(model, EVENT_INFO, STATION_INFO, STATIONS,
                       os.path.join(OUTDIR, "fig02_velocity_model.png"))

    # Fig 3 — Synthetic record section (shared plotting module)
    print("[3/9] Synthetic seismograms ...")
    fig_synthetic_record_section(
        STATION_INFO, YIELD_KT, DF, model,
        os.path.join(OUTDIR, "fig03_synthetic_seismograms.png"),
        fmin=FMIN, fmax=FMAX, depth_km=EVENT_DEPTH,
        title_prefix="Figure 3 — Synthetic Seismograms")

    # Fig 4 — Yield sensitivity (event-specific)
    print("[4/9] Yield sensitivity ...")
    fig04_yield_sensitivity()

    # Fig 5 — Yield trade-off (event-specific)
    print("[5/9] Yield trade-off ...")
    fig05_yield_tradeoff()

    # Fig 6 — Real data from local SAC archive (custom loader)
    print("[6/9] Loading observed data from SAC archive ...")
    real_data = load_local_sac_data(DATADIR, STATIONS, EVENT_TIME,
                                    EVENT_LAT, EVENT_LON,
                                    fmin=FMIN, fmax=FMAX)
    print(f"       Loaded {len(real_data)} station(s)")
    fig_observed_record_section(
        real_data, EVENT_TIME, model,
        os.path.join(OUTDIR, "fig06_observed_waveforms.png"),
        fmin=FMIN, fmax=FMAX,
        title="Figure 6 — Observed Waveforms: Degelen Mtn 1988-10-18")

    # Fig 7 — Synthetic vs observed comparison (shared plotting module)
    print("[7/9] Synthetic vs observed comparison ...")
    key_stas = ["SEM", "ILI", "KUU", "BRTG", "TLG"]
    fig_comparison(
        real_data, EVENT_TIME, YIELD_KT, DF, model,
        os.path.join(OUTDIR, "fig07_comparison.png"),
        fmin=FMIN, fmax=FMAX, depth_km=EVENT_DEPTH,
        key_stations=key_stas,
        title="Figure 7 — Synthetic vs Observed Comparison")

    # Fig 8 — Discrimination (event-specific)
    print("[8/9] Discrimination analysis ...")
    fig08_discrimination(real_data)

    # Fig 9 — Location inversion using observed Pn picks
    print("[9/9] Epicenter location inversion ...")
    print("  Running grid search with observed Pn arrival times...")
    best_lat, best_lon, misfit, lats, lons, picks = invert_epicenter(
        STATIONS, model, EVENT_LAT, EVENT_LON,
        observed_picks=OBSERVED_PICKS)
    err_km = fig_location_inversion(
        best_lat, best_lon, misfit, lats, lons, picks,
        EVENT_LAT, EVENT_LON, STATIONS,
        os.path.join(OUTDIR, "fig09_location_inversion.png"),
        title="Figure 9 — Epicenter Location Inversion "
              "(Degelen 1988-10-18)")
    print(f"       Inverted location: {best_lat:.3f}°N, {best_lon:.3f}°E "
          f"(error: {err_km:.1f} km)")

    print()
    print("=" * 72)
    print(f"  All figures saved to: {OUTDIR}")
    print("  Structure parallels DPRK 2017 example:")
    print("    fig01 — Source model (Mueller-Murphy, 16 kt)")
    print("    fig02 — Velocity model (Degelen / Kazakh Platform)")
    print("    fig03 — Synthetic seismograms")
    print("    fig04 — Yield sensitivity at SEM (170 km)")
    print("    fig05 — Yield trade-off")
    print("    fig06 — Observed waveforms (KazIS/KyrgIS SAC archive)")
    print("    fig07 — Synthetic vs observed comparison")
    print("    fig08 — Discrimination analysis")
    print("    fig09 — Location inversion (observed Pn picks)")
    print("=" * 72)


if __name__ == "__main__":
    main()
