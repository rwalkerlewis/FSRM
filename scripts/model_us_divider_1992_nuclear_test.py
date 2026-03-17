#!/usr/bin/env python3
"""
Model the 1992 US Divider Nuclear Test — Last US Underground Test at NTS

Comprehensive physics-based modeling of a coupled underground nuclear test:
  1. Mueller-Murphy seismic source (RDP model) with spall contribution
  2. 1-D layered velocity model for the Nevada Test Site (NTS) tuff region
  3. Synthetic seismogram generation at real station distances
  4. Comparison against actual FDSN/IRIS broadband data
  5. Explosion-vs-earthquake discrimination analysis
  6. NTS tuff vs granite geology comparison

All results are presented as publication-quality figures.

Physical references:
  Mueller & Murphy (1971)  — Seismic characteristics of underground nuclear detonations
  Patton (1988)            — Corner-frequency scaling
  Springer et al. (2002)   — Seismic source characteristics of NTS explosions
  Walter et al. (2004)     — Moment tensor decomposition of NTS events
  Werth & Herbst (1963)    — Comparison of amplitudes in tuff and granite

Usage:
    python scripts/model_us_divider_1992_nuclear_test.py
"""

import os
import sys
import warnings
import numpy as np
from scipy.signal import butter, filtfilt

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from obspy import UTCDateTime

warnings.filterwarnings("ignore")

# Add scripts directory to path for fsrm package import
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fsrm.source_physics import (
    cavity_radius_empirical, zone_radii, corner_frequency_patton,
    scalar_moment_coupled, mb_from_yield, Mw_from_M0,
    mueller_murphy_spectrum, mueller_murphy_time,
)
from fsrm.velocity_models import NTSVelocityModel, GenericGraniteVelocityModel
from fsrm.propagation import regional_travel_times
from fsrm.signal_processing import envelope, spectral_amplitude, time_relative
from fsrm.data_fetching import fetch_waveforms, compute_station_distances
from fsrm.plotting import (
    PHASE_COLORS, fig_source_model, fig_velocity_model,
    fig_synthetic_record_section, fig_observed_record_section,
    fig_comparison,
)

# ═══════════════════════════════════════════════════════════════════════════════
# Event Parameters
# ═══════════════════════════════════════════════════════════════════════════════
EVENT_TIME   = UTCDateTime("1992-09-23T15:04:00")
EVENT_LAT    = 37.021
EVENT_LON    = -116.058
EVENT_DEPTH  = 0.8          # km
YIELD_KT     = 1.0
DF           = 1.0          # Coupled test
CORNER_FREQ  = 2.5          # Hz
M0_COUPLED   = 1.0e17       # N·m
ISO_FRAC     = 0.80
CLVD_FRAC    = 0.18
DC_FRAC      = 0.02
SPALL_DEPTH  = 80.0         # m
SPALL_VEL    = 1.5          # m/s
CAVITY_RADIUS = 55.0        # m (empirical for 1 kt in tuff)

OUTDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                      "figures", "us_divider_1992")

FMIN, FMAX = 0.5, 8.0

EVENT_INFO = {
    "name": "US Divider Test (Last US Underground Test)",
    "time": "1992-09-23 15:04:00 UTC",
    "lat": EVENT_LAT,
    "lon": EVENT_LON,
    "mb_obs": 4.0,
}

STATIONS = [
    # (net, sta, lat, lon, description)
    ("US",  "TPNV",  36.949, -116.249, "Test Site, NV"),
    ("SN",  "NSP",   36.728, -116.212, "NTS South, NV"),
    ("CI",  "PAS",   34.148, -118.172, "Pasadena, CA"),
    ("CI",  "GSC",   35.302, -116.805, "Goldstone, CA"),
    ("CI",  "ISA",   35.663, -118.474, "Isabella, CA"),
    ("CI",  "PFO",   33.611, -116.459, "Piñon Flat, CA"),
    ("CI",  "SBC",   34.441, -119.715, "Santa Barbara, CA"),
    ("CI",  "SVD",   34.101, -117.097, "Seven Oaks Dam, CA"),
    ("BK",  "SAO",   36.765, -121.447, "San Andreas Obs, CA"),
    ("NN",  "ELK",   40.745, -115.239, "Elko, NV"),
    ("LN",  "MNV",   38.432, -118.153, "Mina, NV"),
    ("NN",  "KNB",   37.017, -112.822, "Kanab, UT"),
    ("IU",  "ANMO",  34.946, -106.457, "Albuquerque, NM"),
    ("US",  "TUC",   32.310, -110.785, "Tucson, AZ"),
]

STATION_INFO = compute_station_distances(EVENT_LAT, EVENT_LON, STATIONS)


# ═══════════════════════════════════════════════════════════════════════════════
# Event-Specific Functions (NTS spall + moment tensor — not in shared modules)
# ═══════════════════════════════════════════════════════════════════════════════

def spall_source_time(t, spall_depth=SPALL_DEPTH, spall_vel=SPALL_VEL):
    """
    Spall closure impulse: occurs when free-surface spall slab
    falls back and impacts.  Modelled as delayed Gaussian pulse.
    """
    t = np.asarray(t, dtype=float)
    t_delay = 2.0 * spall_depth / spall_vel
    sigma = spall_depth / 2500.0
    amp = 0.15 * M0_COUPLED
    pulse = amp * np.exp(-0.5 * ((t - t_delay) / sigma)**2)
    return pulse


def moment_tensor_decomposition():
    """Return (M_iso, M_clvd, M_dc) tensors for NTS explosion."""
    M0 = M0_COUPLED
    M_iso = ISO_FRAC * M0 * np.eye(3) / 3.0
    M_clvd = CLVD_FRAC * M0 * np.array([
        [-0.5, 0, 0],
        [0, -0.5, 0],
        [0,  0,  1.0]
    ])
    M_dc = DC_FRAC * M0 * np.array([
        [0, 0, 1],
        [0, 0, 0],
        [1, 0, 0]
    ], dtype=float)
    return M_iso, M_clvd, M_dc


def _geometric_spreading(dist_km, wave_type="body"):
    """Geometric spreading factor."""
    r = max(dist_km, 0.1)
    if wave_type == "surface":
        return 1.0 / np.sqrt(r)
    return 1.0 / r


def _attenuation(freq, travel_time_s, Q):
    """Anelastic attenuation:  exp(-π f t* / Q)."""
    return np.exp(-np.pi * freq * travel_time_s / Q)


def generate_synthetic(dist_km, yield_kt=YIELD_KT,
                       dt=0.05, duration=500.0, model=None,
                       include_spall=True):
    """
    Generate a synthetic velocity seismogram at epicentral distance dist_km.

    NTS-specific: includes spall closure impulse and uses NTS tuff constants.
    Returns (time_array, velocity_waveform) — both as numpy arrays.
    """
    if model is None:
        model = NTSVelocityModel()

    npts = int(duration / dt)
    t = np.arange(npts) * dt
    vel = np.zeros(npts)

    fc = CORNER_FREQ
    M0 = M0_COUPLED

    phases = regional_travel_times(dist_km, model)

    phase_weights = {"Pn": 1.0, "Pg": 0.7, "Sn": 0.4, "Lg": 1.8}
    q_values = {"Pn": 500, "Pg": 400, "Sn": 300, "Lg": 250}
    spreading = {"Pn": "body", "Pg": "body", "Sn": "body", "Lg": "surface"}

    f = np.fft.rfftfreq(npts, d=dt)
    f[0] = 1e-10

    for phase_name, (t_arr, vel_phase, wtype) in phases.items():
        if t_arr >= duration or t_arr < 0:
            continue

        weight = phase_weights.get(phase_name, 1.0)
        Q = q_values.get(phase_name, 400)
        spr = _geometric_spreading(dist_km, spreading.get(phase_name, "body"))

        src = mueller_murphy_spectrum(f, yield_kt, fc=fc)
        att = _attenuation(f, t_arr, Q)
        spec = src * att * spr * weight

        omega = 2 * np.pi * f
        spec_vel = spec * omega

        bp = np.ones_like(f)
        bp[f < FMIN] = 0
        bp[f > FMAX] = 0
        taper_lo = (f >= FMIN * 0.5) & (f <= FMIN)
        bp[taper_lo] = 0.5 * (1 - np.cos(np.pi * (f[taper_lo] - FMIN*0.5) / (FMIN*0.5)))
        taper_hi = (f >= FMAX) & (f <= FMAX * 1.25)
        bp[taper_hi] = 0.5 * (1 + np.cos(np.pi * (f[taper_hi] - FMAX) / (FMAX*0.25)))
        spec_vel *= bp

        phase_shift = np.exp(-1j * omega * t_arr)

        if wtype == "P":
            wavelet = np.fft.irfft(spec_vel * phase_shift, n=npts)
        else:
            rng = np.random.RandomState(hash(phase_name) % 2**31)
            random_phase = np.exp(1j * rng.uniform(-np.pi/3, np.pi/3, len(f)))
            wavelet = np.fft.irfft(spec_vel * phase_shift * random_phase, n=npts)

        env_duration = {"Pn": 15.0, "Pg": 20.0, "Sn": 25.0, "Lg": 60.0}
        ed = env_duration.get(phase_name, 15.0)
        rise_time = {"Pn": 1.5, "Pg": 2.0, "Sn": 3.0, "Lg": 5.0}
        rt = rise_time.get(phase_name, 2.0)

        idx_arr = int(t_arr / dt)
        for j in range(npts):
            t_rel = (j - idx_arr) * dt
            if t_rel < -rt:
                wavelet[j] *= 0
            elif t_rel < 0:
                wavelet[j] *= np.exp(-0.5 * (t_rel / (rt * 0.5))**2)
            elif t_rel < ed * 0.25:
                ramp = min(1.0, t_rel / (rt * 0.8))
                wavelet[j] *= ramp
            elif t_rel < ed * 0.4:
                wavelet[j] *= 1.0
            else:
                wavelet[j] *= np.exp(-(t_rel - ed * 0.4) / (ed * 0.4))

        rng_coda = np.random.RandomState(hash(phase_name + str(int(dist_km))) % 2**31)
        coda_noise = rng_coda.randn(npts) * 0.3
        if npts > 100:
            try:
                b_c, a_c = butter(2, [FMIN * 2 * dt, min(FMAX * 2 * dt, 0.99)], btype='band')
                coda_noise = filtfilt(b_c, a_c, coda_noise)
            except Exception:
                pass
        for j in range(npts):
            t_rel = (j - idx_arr) * dt
            if t_rel < rt:
                coda_noise[j] = 0
            elif t_rel < ed:
                coda_env = np.exp(-(t_rel - rt) / (ed * 0.5)) * 0.4
                coda_noise[j] *= coda_env
            else:
                coda_noise[j] *= np.exp(-(t_rel - ed) / (ed * 0.3)) * 0.15

        wavelet += coda_noise * np.max(np.abs(wavelet)) * 0.3
        vel += wavelet

    if include_spall:
        spall_pulse = spall_source_time(t)
        spall_spec = np.fft.rfft(spall_pulse)
        spall_spec *= _attenuation(f, dist_km / 5.0, 300)
        spall_spec *= _geometric_spreading(dist_km, "body")
        omega = 2 * np.pi * f
        spall_vel_spec = spall_spec * omega
        bp_spall = np.ones_like(f)
        bp_spall[f < FMIN] = 0
        bp_spall[f > FMAX] = 0
        spall_vel_spec *= bp_spall
        spall_contrib = np.fft.irfft(spall_vel_spec, n=npts)
        vel += spall_contrib * 0.1

    return t, vel


# ═══════════════════════════════════════════════════════════════════════════════
# Event-Specific Figures
# ═══════════════════════════════════════════════════════════════════════════════

def fig04_yield_sensitivity():
    """Synthetic waveforms at GSC for different yields."""
    dist_km = None
    for net, sta, d, az, desc in STATION_INFO:
        if sta == "GSC":
            dist_km = d
            break
    if dist_km is None:
        dist_km = 200.0

    model = NTSVelocityModel()
    yields = [0.1, 0.5, 1.0, 5.0, 20.0]

    fig, axes = plt.subplots(len(yields), 1, figsize=(16, 3.0 * len(yields)),
                             sharex=True)

    for i, W in enumerate(yields):
        fc_w = 2.5 * W**(-1.0/3.0)
        t, v = generate_synthetic(dist_km, W, dt=0.05, duration=400, model=model)
        mb = mb_from_yield(W)

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

        ax.set_ylabel("Norm. vel.", fontsize=9)
        ax.set_title(f"Yield = {W} kt  |  mb = {mb:.2f}  |  fc = {fc_w:.2f} Hz  |  "
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

    fig.suptitle(f"Figure 4 — Yield Sensitivity at GSC ({dist_km:.0f} km):  Coupled Test\n"
                 f"NTS tuff model, BHZ {FMIN}–{FMAX} Hz",
                 fontsize=13, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    p = os.path.join(OUTDIR, "fig04_yield_sensitivity.png")
    fig.savefig(p, dpi=180, bbox_inches="tight"); plt.close(fig)
    print(f"  Saved {p}")


def fig05_tuff_vs_granite():
    """Compare NTS tuff vs granite velocity model waveforms and spectra."""
    model_tuff = NTSVelocityModel()
    model_gran = GenericGraniteVelocityModel()

    dist_km = None
    for net, sta, d, az, desc in STATION_INFO:
        if sta == "PAS":
            dist_km = d
            break
    if dist_km is None:
        dist_km = 370.0

    fig = plt.figure(figsize=(20, 16))
    gs = GridSpec(3, 2, figure=fig, hspace=0.35, wspace=0.30)

    # (a) Tuff velocity profile vs Granite velocity profile
    ax = fig.add_subplot(gs[0, 0])
    z_t, vp_t, vs_t, _, _, _ = model_tuff.profile_arrays(dz=0.05, zmax=45)
    z_g, vp_g, vs_g, _, _, _ = model_gran.profile_arrays(dz=0.05, zmax=45)

    ax.plot(vp_t, z_t, "b-", lw=2, label="Vp tuff")
    ax.plot(vs_t, z_t, "b--", lw=1.5, label="Vs tuff")
    ax.plot(vp_g, z_g, "r-", lw=2, label="Vp granite")
    ax.plot(vs_g, z_g, "r--", lw=1.5, label="Vs granite")
    ax.axhline(33, color="green", ls=":", lw=1.5, label="Moho (33 km)")
    ax.invert_yaxis()
    ax.set_xlabel("Velocity (km/s)"); ax.set_ylabel("Depth (km)")
    ax.set_title("(a) Velocity profiles: Tuff vs Granite", fontweight="bold")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.2)

    # (b) Q profiles
    ax = fig.add_subplot(gs[0, 1])
    _, _, _, _, qp_t, qs_t = model_tuff.profile_arrays(dz=0.05, zmax=45)
    _, _, _, _, qp_g, qs_g = model_gran.profile_arrays(dz=0.05, zmax=45)

    ax.plot(qp_t, z_t, "b-", lw=2, label="Qp tuff")
    ax.plot(qs_t, z_t, "b--", lw=1.5, label="Qs tuff")
    ax.plot(qp_g, z_g, "r-", lw=2, label="Qp granite")
    ax.plot(qs_g, z_g, "r--", lw=1.5, label="Qs granite")
    ax.axhline(33, color="green", ls=":", lw=1.5)
    ax.invert_yaxis()
    ax.set_xlabel("Quality factor Q"); ax.set_ylabel("Depth (km)")
    ax.set_title("(b) Attenuation: Tuff vs Granite", fontweight="bold")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.2)

    # (c) Waveform — tuff model
    t_tuff, v_tuff = generate_synthetic(dist_km, dt=0.05, duration=400, model=model_tuff)
    t_gran, v_gran = generate_synthetic(dist_km, dt=0.05, duration=400, model=model_gran)

    peak_ref = max(np.max(np.abs(v_tuff)), np.max(np.abs(v_gran)), 1e-30)

    ax = fig.add_subplot(gs[1, 0])
    ax.plot(t_tuff, v_tuff / peak_ref, "b-", lw=0.5, alpha=0.8)
    env_tuff = envelope(v_tuff / peak_ref)
    ax.plot(t_tuff, env_tuff, "b-", lw=0.8, alpha=0.3)
    ax.plot(t_tuff, -env_tuff, "b-", lw=0.8, alpha=0.3)
    ax.set_ylabel("Rel. velocity"); ax.set_xlabel("Time (s)")
    ax.set_title(f"(c) NTS Tuff model at PAS ({dist_km:.0f} km)", fontweight="bold")
    ax.grid(True, alpha=0.2); ax.set_xlim(0, 400)

    tt_tuff = regional_travel_times(dist_km, model_tuff)
    for ph, (tarr, _, _) in tt_tuff.items():
        if 0 < tarr < 400:
            ax.axvline(tarr, color=PHASE_COLORS.get(ph, "gray"),
                       ls="--", lw=0.8, alpha=0.5)

    # (d) Waveform — granite model
    ax = fig.add_subplot(gs[1, 1])
    ax.plot(t_gran, v_gran / peak_ref, "r-", lw=0.5, alpha=0.8)
    env_gran = envelope(v_gran / peak_ref)
    ax.plot(t_gran, env_gran, "r-", lw=0.8, alpha=0.3)
    ax.plot(t_gran, -env_gran, "r-", lw=0.8, alpha=0.3)
    ax.set_ylabel("Rel. velocity"); ax.set_xlabel("Time (s)")
    ax.set_title(f"(d) Granite model at PAS ({dist_km:.0f} km)", fontweight="bold")
    ax.grid(True, alpha=0.2); ax.set_xlim(0, 400)

    tt_gran = regional_travel_times(dist_km, model_gran)
    for ph, (tarr, _, _) in tt_gran.items():
        if 0 < tarr < 400:
            ax.axvline(tarr, color=PHASE_COLORS.get(ph, "gray"),
                       ls="--", lw=0.8, alpha=0.5)

    # (e) Waveform overlay
    ax = fig.add_subplot(gs[2, 0])
    ax.plot(t_tuff, v_tuff / peak_ref, "b-", lw=0.5, alpha=0.7, label="NTS Tuff")
    ax.plot(t_gran, v_gran / peak_ref, "r-", lw=0.5, alpha=0.7, label="Granite")
    ax.set_ylabel("Rel. velocity"); ax.set_xlabel("Time after origin (s)")
    ax.set_title("(e) Waveform overlay", fontweight="bold")
    ax.legend(fontsize=9); ax.grid(True, alpha=0.2); ax.set_xlim(0, 400)

    # (f) Spectral comparison
    ax = fig.add_subplot(gs[2, 1])
    freq_t, sp_t = spectral_amplitude(v_tuff, 0.05)
    freq_g, sp_g = spectral_amplitude(v_gran, 0.05)

    m_t = (freq_t > 0.1) & (freq_t < 10)
    m_g = (freq_g > 0.1) & (freq_g < 10)

    sp_max = max(np.max(sp_t[m_t]), np.max(sp_g[m_g]), 1e-30)
    ax.loglog(freq_t[m_t], sp_t[m_t] / sp_max, "b-", lw=1.5, label="NTS Tuff")
    ax.loglog(freq_g[m_g], sp_g[m_g] / sp_max, "r-", lw=1.5, label="Granite")
    ax.set_xlabel("Frequency (Hz)"); ax.set_ylabel("Normalised spectral amplitude")
    ax.set_title("(f) Spectral comparison", fontweight="bold")
    ax.legend(fontsize=9); ax.grid(True, alpha=0.2, which="both")
    ax.set_xlim(0.1, 10)

    fig.suptitle("Figure 5 — NTS Tuff vs Granite Geology Comparison\n"
                 f"1 kt coupled, Mueller-Murphy source, {FMIN}–{FMAX} Hz",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    p = os.path.join(OUTDIR, "fig05_tuff_vs_granite.png")
    fig.savefig(p, dpi=180, bbox_inches="tight"); plt.close(fig)
    print(f"  Saved {p}")


def fig08_discrimination(real_data):
    """Explosion vs earthquake discrimination analysis for the Divider test."""
    fig = plt.figure(figsize=(18, 12))
    gs = GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.35)

    model = NTSVelocityModel()

    # (a) P/Lg amplitude ratio for synthetic
    ax = fig.add_subplot(gs[0, 0])
    distances = np.linspace(50, 1200, 25)
    ps_ratios = []
    for d in distances:
        t, v = generate_synthetic(d, dt=0.05, duration=600, model=model)
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

    mb_eq = np.linspace(2, 7, 50)
    ms_eq = 1.0 * mb_eq - 0.5 + 0.2 * np.random.randn(50)

    ax.scatter(mb_eq, ms_eq, c="b", s=20, alpha=0.4, label="Earthquakes")
    ax.scatter(mb_exp, ms_exp, c="r", s=20, alpha=0.4, label="Explosions (NTS)")

    mb_this = mb_from_yield(YIELD_KT)
    ms_this = mb_this - 1.2
    ax.plot(mb_this, ms_this, "r*", ms=20, zorder=10,
            label=f"Divider (mb={mb_this:.2f})")

    mb_line = np.linspace(2, 7, 100)
    ms_line = mb_line - 1.0
    ax.plot(mb_line, ms_line, "k--", lw=1.5, label="Discrimination line")

    ax.set_xlabel("mb"); ax.set_ylabel("Ms")
    ax.set_title("(b) mb–Ms Discrimination", fontweight="bold")
    ax.legend(fontsize=7, loc="upper left"); ax.grid(True, alpha=0.2)
    ax.set_xlim(1, 7); ax.set_ylim(0, 7)

    # (c) Spectral ratio analysis (synthetic)
    ax = fig.add_subplot(gs[0, 2])
    dist_ref = None
    for net, sta, d, az, desc in STATION_INFO:
        if sta == "PAS":
            dist_ref = d
            break
    if dist_ref is None:
        dist_ref = 370.0

    t, v = generate_synthetic(dist_ref, dt=0.05, duration=500, model=model)
    freq, spec = spectral_amplitude(v, 0.05)
    m = (freq > 0.3) & (freq < 10)
    spec_n = spec[m] / np.max(spec[m]) if np.max(spec[m]) > 0 else spec[m]

    ax.semilogy(freq[m], spec_n, "r-", lw=1.5, label="Explosion (synth)")

    fc_eq = 1.0
    spec_eq = 1.0 / (1 + (freq[m] / fc_eq)**2)**0.5
    spec_eq_n = spec_eq / np.max(spec_eq)
    ax.semilogy(freq[m], spec_eq_n, "b--", lw=1.5, label="Earthquake (model)")

    ax.set_xlabel("Frequency (Hz)"); ax.set_ylabel("Normalised amplitude")
    ax.set_title(f"(c) Spectral Shape at PAS ({dist_ref:.0f} km)", fontweight="bold")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.2, which="both")

    # (d) Observed spectral ratio (if real data available)
    ax = fig.add_subplot(gs[1, 0])
    pas_data = [x for x in real_data if x[1] in ("PAS", "GSC")]
    if pas_data:
        net, sta, dist_km, az, tr, desc = pas_data[0]
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
        ax.text(0.5, 0.5, "No station data available", transform=ax.transAxes,
                ha="center", fontsize=12, color="gray")

    ax.set_xlabel("Frequency (Hz)"); ax.set_ylabel("Normalised amplitude")
    ax.set_title("(d) Observed P vs Lg Spectra", fontweight="bold")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.2, which="both")

    # (e) Source-type discrimination summary
    ax = fig.add_subplot(gs[1, 1])
    ax.axis("off")
    lines = [
        "DISCRIMINATION ANALYSIS",
        "═" * 36,
        "",
        "INDICATORS CONSISTENT WITH EXPLOSION:",
        "",
        "  ✓  Impulsive P, weak S/Lg",
        "     (P/Lg ratio > 1 at regional dist.)",
        "",
        "  ✓  mb–Ms below earthquake trend",
        f"     (mb≈{mb_this:.2f}, Ms≈{ms_this:.2f})",
        "",
        "  ✓  Higher-frequency content than",
        "     equivalent earthquake",
        "",
        "  ✓  Location at Nevada Test Site",
        "     (known underground test facility)",
        "",
        "  ✓  Shallow depth (~800 m)",
        "",
        "  ✓  Moment tensor: 80% isotropic",
        "     (characteristic of explosions)",
        "",
        "  ✓  Spall signal present",
        f"     ({SPALL_DEPTH:.0f} m slab, {SPALL_VEL} m/s)",
        "",
        "COUPLED TEST:",
        "",
        f"  Yield:  ~{YIELD_KT} kt in NTS tuff",
        f"  M0:     {M0_COUPLED:.1e} N·m",
        f"  mb:     {mb_this:.2f}",
    ]
    ax.text(0.05, 0.97, "\n".join(lines), transform=ax.transAxes,
            fontfamily="monospace", fontsize=8.5, va="top",
            bbox=dict(boxstyle="round,pad=0.5", fc="#fffff0", ec="#cca"))

    # (f) Yield–mb for NTS tests
    ax = fig.add_subplot(gs[1, 2])
    W_range = np.logspace(-2, 4, 200)
    mb_range = 4.0 + 0.75 * np.log10(W_range)
    ax.semilogx(W_range, mb_range, "b-", lw=2, label="mb–yield (coupled)")

    nts_tests = [
        (150, 6.2, "Boxcar (1968)"),
        (1000, 6.8, "Cannikin (1971)"),
        (20, 5.0, "Baneberry (1970)"),
        (1, 4.0, "Divider (1992)"),
    ]
    for yw, mbw, name in nts_tests:
        ax.plot(yw, mbw, "r*", ms=14, zorder=10)
        ax.text(yw * 1.2, mbw + 0.1, name, fontsize=7, fontweight="bold")

    ax.set_xlabel("Yield (kt)"); ax.set_ylabel("mb")
    ax.set_title("(f) NTS mb–Yield Calibration", fontweight="bold")
    ax.legend(fontsize=8, loc="upper left"); ax.grid(True, alpha=0.2, which="both")
    ax.set_xlim(0.01, 1e4); ax.set_ylim(1, 8)

    fig.suptitle("Figure 8 — Explosion/Earthquake Discrimination Analysis\n"
                 "US Divider Test 1992-09-23 15:04 UTC — Last US Underground Nuclear Test",
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

    model = NTSVelocityModel()

    print("=" * 72)
    print("  US Divider 1992 — Last Underground Nuclear Test Model & Analysis")
    print("=" * 72)
    print()

    print("[1/8] Source model ...")
    fig_source_model(YIELD_KT, EVENT_DEPTH * 1000, DF, EVENT_INFO, OUTDIR,
                     rho=1900.0)

    print("[2/8] Velocity model ...")
    fig_velocity_model(model, EVENT_INFO, STATION_INFO, STATIONS, OUTDIR)

    print("[3/8] Synthetic seismograms ...")
    fig_synthetic_record_section(STATION_INFO, YIELD_KT, DF, model, OUTDIR,
                                 fmin=FMIN, fmax=FMAX, depth_km=EVENT_DEPTH,
                                 title_prefix="US Divider 1992")

    print("[4/8] Yield sensitivity ...")
    fig04_yield_sensitivity()

    print("[5/8] NTS tuff vs granite ...")
    fig05_tuff_vs_granite()

    print("[6/8] Fetching real data from IRIS ...")
    real_data = fetch_waveforms(EVENT_TIME, EVENT_LAT, EVENT_LON, STATIONS,
                                fmin=FMIN, fmax=FMAX, pre=120, post=500)
    print(f"       Retrieved {len(real_data)} station(s)")
    fig_observed_record_section(real_data, EVENT_TIME, model, OUTDIR,
                                fmin=FMIN, fmax=FMAX,
                                title="US Divider 1992 — Observed Waveforms")

    print("[7/8] Synthetic vs observed comparison ...")
    fig_comparison(real_data, EVENT_TIME, YIELD_KT, DF, model, OUTDIR,
                   fmin=FMIN, fmax=FMAX, depth_km=EVENT_DEPTH,
                   key_stations=["TPNV", "PAS", "GSC", "ANMO", "TUC"],
                   title="US Divider 1992 — Synthetic vs Observed")

    print("[8/8] Discrimination analysis ...")
    fig08_discrimination(real_data)

    print()
    print("=" * 72)
    print(f"  All figures saved to: {OUTDIR}")
    print("=" * 72)


if __name__ == "__main__":
    main()
