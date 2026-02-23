#!/usr/bin/env python3
"""
Model the DPRK 6th Nuclear Test — 2017-09-03 (Mt. Mantap)

Comprehensive physics-based modeling of a fully coupled underground nuclear test:
  1. Mueller-Murphy seismic source (RDP model) — 250 kt, no decoupling
  2. 1-D layered velocity model for the Punggye-ri / Korean Peninsula region
  3. Synthetic seismogram generation at real station distances
  4. Comparison against actual FDSN/IRIS broadband data
  5. Explosion-vs-earthquake discrimination analysis

All results are presented as publication-quality figures.

Physical references:
  Mueller & Murphy (1971)  — Seismic characteristics of underground nuclear detonations
  Patton (1988)            — Corner-frequency scaling
  Kim & Richards (2007)    — Lg blockage at Korean Peninsula
  Zhao et al. (2017)       — Yield estimation for the 2017 DPRK test
  Wei (2017)               — Mt. Mantap collapse analysis
  NORSAR (2017)            — DPRK-6 detection bulletin

Usage:
    python scripts/model_dprk_2017_nuclear_test.py
"""

import os
import sys
import warnings
import numpy as np
from scipy.signal import hilbert, butter, filtfilt

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.ticker import AutoMinorLocator
from matplotlib.gridspec import GridSpec

from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.geodetics import gps2dist_azimuth

warnings.filterwarnings("ignore")

# ═══════════════════════════════════════════════════════════════════════════════
# Constants & Event Parameters
# ═══════════════════════════════════════════════════════════════════════════════
JOULES_PER_KT = 4.184e12

EVENT_TIME   = UTCDateTime("2017-09-03T03:30:01")
EVENT_LAT    = 41.300
EVENT_LON    = 129.076
EVENT_DEPTH  = 0.76          # km

OUTDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                      "figures", "dprk_2017")

# Bandpass for processing
FMIN, FMAX = 0.5, 8.0

# Stations for synthetic + real comparison
STATIONS = [
    # (net, sta, lat, lon, description)
    ("IC",  "MDJ",  44.616, 129.592, "Mudanjiang, China"),
    ("KS",  "KSRS", 36.539, 127.867, "KSRS Array, S. Korea"),
    ("KG",  "INCN", 37.479, 126.624, "Incheon, S. Korea"),
    ("IC",  "BJT",  40.018, 116.168, "Baijiatuan, China"),
    ("IC",  "HIA",  49.267, 119.742, "Hailar, China"),
    ("IC",  "SSE",  31.096, 121.187, "Shanghai, China"),
    ("IU",  "MAJO", 36.545, 138.204, "Matsushiro, Japan"),
    ("II",  "ERM",  42.015, 143.157, "Erimo, Japan"),
    ("IU",  "YSS",  46.959, 142.760, "Yuzhno, Russia"),
    ("IC",  "ENH",  30.275, 109.494, "Enshi, China"),
    ("IU",  "ULN",  47.865, 107.053, "Ulaanbaatar, Mongolia"),
    ("II",  "AAK",  42.638,  74.494, "Ala Archa, Kyrgyzstan"),
]

def _station_distances():
    """Return sorted list of (net, sta, dist_km, az, desc)."""
    out = []
    for net, sta, slat, slon, desc in STATIONS:
        d_m, az, _ = gps2dist_azimuth(EVENT_LAT, EVENT_LON, slat, slon)
        out.append((net, sta, d_m / 1e3, az, desc))
    out.sort(key=lambda x: x[2])
    return out

STATION_INFO = _station_distances()

# ═══════════════════════════════════════════════════════════════════════════════
# PART A — Source Physics
# ═══════════════════════════════════════════════════════════════════════════════

def cavity_radius_empirical(yield_kt, rho=2700.0):
    """NTS empirical:  R_c = 55 · W^0.295 · (ρ/2650)^(-1/3.4)  [m]
    For Korean granite, ρ ≈ 2700 kg/m³."""
    return 55.0 * yield_kt**0.295 * (rho / 2650.0)**(-1.0 / 3.4)

def zone_radii(Rc):
    return dict(cavity=Rc, crushed=2.5*Rc, fractured=5.0*Rc, damaged=10.0*Rc)

def corner_frequency_patton(yield_kt):
    """Patton (1988):  fc ≈ 2.5 · W^(-1/3)  Hz"""
    return 2.5 * yield_kt**(-1.0 / 3.0)

def scalar_moment_coupled(yield_kt):
    """Empirical:  log10(M0) ≈ 17.0 + log10(W_kt)  → N·m"""
    return 10.0**(17.0 + np.log10(max(yield_kt, 1e-6)))

def mb_from_yield(yield_kt, decoupling_factor=1.0):
    """mb ≈ 4.0 + 0.75·log10(W)  with decoupling reduction."""
    w_eff = yield_kt / decoupling_factor
    return 4.0 + 0.75 * np.log10(max(w_eff, 1e-6))

def Mw_from_M0(M0):
    return (np.log10(M0) - 9.1) / 1.5

def mueller_murphy_spectrum(f, yield_kt, fc=None, overshoot=1.1):
    """
    Mueller-Murphy source spectrum  Ψ(f) = Ψ_∞ · B / [1 + (f/fc)²]
    Returns displacement spectral amplitude (arbitrary units, relative).
    """
    if fc is None:
        fc = corner_frequency_patton(yield_kt)
    psi_inf = scalar_moment_coupled(yield_kt)
    return overshoot * psi_inf / (1.0 + (f / fc)**2)

def mueller_murphy_time(t, yield_kt, fc=None, overshoot=1.1):
    """
    Reduced Displacement Potential ψ(t) and moment-rate dψ/dt.
    Brune-style:  ψ(t) = Ψ_∞·B·[1 - (1+t/τ)·exp(-t/τ)]  for t ≥ 0.
    """
    if fc is None:
        fc = corner_frequency_patton(yield_kt)
    tau = 1.0 / (2.0 * np.pi * fc)
    psi_inf = scalar_moment_coupled(yield_kt)
    t = np.asarray(t, dtype=float)
    psi = np.zeros_like(t)
    dpsi = np.zeros_like(t)
    m = t > 0
    x = t[m] / tau
    psi[m]  = overshoot * psi_inf * (1.0 - (1.0 + x) * np.exp(-x))
    dpsi[m] = overshoot * psi_inf * x * np.exp(-x) / tau
    return psi, dpsi


# ═══════════════════════════════════════════════════════════════════════════════
# PART B — 1-D Velocity Model (Punggye-ri / Korean Peninsula)
# ═══════════════════════════════════════════════════════════════════════════════

class PunggyeRiVelocityModel:
    """
    Layered 1-D velocity model for the Punggye-ri test site, Korean Peninsula.

    Layers (top-of-layer depth, Vp, Vs, rho, Qp, Qs):
        0.00 km  Weathered surface       3.50  2.00  2100   80   40
        0.10 km  Upper granite/syenite   5.60  3.25  2700  350  175
        5.00 km  Upper crust             6.10  3.55  2800  500  250
       18.00 km  Middle crust            6.50  3.75  2900  600  300
       28.00 km  Lower crust             7.00  4.05  3100  700  350
       33.00 km  Upper mantle            8.05  4.55  3350 1000  500
    """
    _layers = np.array([
        # ztop(km)  Vp(km/s)  Vs(km/s)  rho(kg/m³)  Qp    Qs
        [  0.00,     3.50,     2.00,     2100,        80,   40],
        [  0.10,     5.60,     3.25,     2700,       350,  175],
        [  5.00,     6.10,     3.55,     2800,       500,  250],
        [ 18.00,     6.50,     3.75,     2900,       600,  300],
        [ 28.00,     7.00,     4.05,     3100,       700,  350],
        [ 33.00,     8.05,     4.55,     3350,      1000,  500],
    ])

    def __init__(self):
        self.layers = self._layers.copy()
        self.n = len(self.layers)

    def _layer_index(self, z_km):
        for i in range(self.n - 1, -1, -1):
            if z_km >= self.layers[i, 0]:
                return i
        return 0

    def get(self, z_km):
        """Return (Vp, Vs, rho, Qp, Qs) at depth z_km."""
        i = self._layer_index(z_km)
        return tuple(self.layers[i, 1:])

    @property
    def moho_depth(self):
        return 33.0   # km

    @property
    def pn_velocity(self):
        return 8.05    # km/s

    def taupy_npz_string(self):
        """Build an ObsPy TauPyModel-compatible .npz layer specification."""
        pass

    def profile_arrays(self, dz=0.1, zmax=50):
        """Return (z, Vp, Vs, rho, Qp, Qs) arrays for plotting."""
        z = np.arange(0, zmax, dz)
        vp = np.zeros_like(z); vs = np.zeros_like(z)
        rho = np.zeros_like(z); qp = np.zeros_like(z); qs = np.zeros_like(z)
        for j, zi in enumerate(z):
            vp[j], vs[j], rho[j], qp[j], qs[j] = self.get(zi)
        return z, vp, vs, rho, qp, qs


# ═══════════════════════════════════════════════════════════════════════════════
# PART C — Synthetic Seismogram Generation
# ═══════════════════════════════════════════════════════════════════════════════

def _travel_time(dist_km, velocity_km_s):
    return dist_km / velocity_km_s

def regional_travel_times(dist_km, model=None):
    """
    Approximate travel times for regional seismic phases.
    Returns dict  {phase_name: (arrival_time_s, velocity_km_s, wave_type)}.
    """
    if model is None:
        model = PunggyeRiVelocityModel()
    moho = model.moho_depth
    pn = model.pn_velocity

    vp_crust = 6.1    # rough average P for Korean Peninsula crust
    vs_crust = 3.55   # rough average S

    # Pn — refraction along Moho
    cos_ic_p = vp_crust / pn
    sin_ic_p = np.sqrt(1 - cos_ic_p**2) if cos_ic_p < 1 else 0.01
    t_pn = 2 * moho * cos_ic_p / vp_crust + dist_km / pn

    # Pg — direct crustal P
    t_pg = dist_km / vp_crust

    # Sn — S refraction along Moho
    sn_vel = 4.55
    cos_ic_s = vs_crust / sn_vel
    t_sn = 2 * moho * cos_ic_s / vs_crust + dist_km / sn_vel

    # Lg — crustal guided S/surface wave
    lg_vel = 3.5
    t_lg = dist_km / lg_vel

    return {
        "Pn": (t_pn, pn,       "P"),
        "Pg": (t_pg, vp_crust, "P"),
        "Sn": (t_sn, sn_vel,   "S"),
        "Lg": (t_lg, lg_vel,   "S"),
    }


def _geometric_spreading(dist_km, wave_type="body"):
    """Geometric spreading factor."""
    r = max(dist_km, 0.1)
    if wave_type == "surface":
        return 1.0 / np.sqrt(r)
    return 1.0 / r


def _attenuation(freq, travel_time_s, Q):
    """Anelastic attenuation:  exp(-π f t* / Q)."""
    return np.exp(-np.pi * freq * travel_time_s / Q)


def generate_synthetic(dist_km, yield_kt, decoupling_factor=1.0,
                       dt=0.05, duration=500.0, model=None):
    """
    Generate a synthetic velocity seismogram at epicentral distance dist_km.

    Sums Pn, Pg, Sn, Lg phases, each convolved with the Mueller-Murphy
    source spectrum and propagated with geometric spreading + Q attenuation.

    For the DPRK 2017 test, DF=1.0 (fully coupled).

    Returns (time_array, velocity_waveform) — both as numpy arrays.
    """
    if model is None:
        model = PunggyeRiVelocityModel()

    npts = int(duration / dt)
    t = np.arange(npts) * dt
    vel = np.zeros(npts)

    fc = corner_frequency_patton(yield_kt)
    M0 = scalar_moment_coupled(yield_kt) / decoupling_factor

    phases = regional_travel_times(dist_km, model)

    phase_weights = {"Pn": 1.0, "Pg": 0.7, "Sn": 0.4, "Lg": 1.8}
    q_values = {"Pn": 600, "Pg": 450, "Sn": 350, "Lg": 280}
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

        src = src * (M0 / scalar_moment_coupled(yield_kt))

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
        n_scatter = int(ed / (3 * dt))
        coda_noise = rng_coda.randn(npts) * 0.3
        if npts > 100:
            from scipy.signal import butter as _butter, filtfilt as _filtfilt
            try:
                b_c, a_c = _butter(2, [FMIN * 2 * dt, min(FMAX * 2 * dt, 0.99)], btype='band')
                coda_noise = _filtfilt(b_c, a_c, coda_noise)
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

    return t, vel


# ═══════════════════════════════════════════════════════════════════════════════
# PART D — Real Data Retrieval
# ═══════════════════════════════════════════════════════════════════════════════

def fetch_real_waveforms(stations_to_fetch=None, pre=120, post=500):
    """
    Download BHZ waveforms from IRIS, remove response, bandpass filter.
    The DPRK 2017 test (mb ~6.3) is very well recorded globally.
    Returns list of (net, sta, dist_km, az, trace, desc).
    """
    if stations_to_fetch is None:
        stations_to_fetch = STATIONS

    client = Client("IRIS", timeout=60)
    t1 = EVENT_TIME - pre
    t2 = EVENT_TIME + post
    results = []

    for net, sta, slat, slon, desc in stations_to_fetch:
        try:
            st = client.get_waveforms(net, sta, "*", "BHZ", t1, t2)
            if len(st) == 0:
                continue
            st.merge(fill_value="interpolate")
            tr = st[0]

            inv = client.get_stations(network=net, station=sta,
                                      starttime=EVENT_TIME, endtime=EVENT_TIME,
                                      level="response", channel="BHZ")

            d_m, az, _ = gps2dist_azimuth(EVENT_LAT, EVENT_LON, slat, slon)
            dist_km = d_m / 1e3

            tr = tr.copy()
            tr.detrend("demean"); tr.detrend("linear")
            tr.taper(0.05, type="cosine")
            try:
                tr.remove_response(inventory=inv, output="VEL",
                                   water_level=60,
                                   pre_filt=[FMIN*0.5, FMIN, FMAX, FMAX*1.25])
            except Exception:
                pass
            tr.filter("bandpass", freqmin=FMIN, freqmax=FMAX,
                      corners=4, zerophase=True)

            results.append((net, sta, dist_km, az, tr, desc))
        except Exception:
            pass

    results.sort(key=lambda x: x[2])
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# PART E — Waveform Comparison Utilities
# ═══════════════════════════════════════════════════════════════════════════════

def _time_relative(tr, origin):
    npts = tr.stats.npts
    dt = tr.stats.delta
    return np.arange(npts) * dt + (tr.stats.starttime - origin)

def envelope(data):
    return np.abs(hilbert(data))

def spectral_amplitude(data, dt):
    """Return (freq, amplitude) via FFT."""
    n = len(data)
    spec = np.abs(np.fft.rfft(data)) * dt
    freq = np.fft.rfftfreq(n, d=dt)
    return freq, spec


# ═══════════════════════════════════════════════════════════════════════════════
# PART F — Figure Generation
# ═══════════════════════════════════════════════════════════════════════════════

def fig01_source_model():
    """Source model summary: cavity, RDP, spectrum, mb-yield, damage, parameters."""
    fig = plt.figure(figsize=(18, 12))
    gs = GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.35)

    yield_kt = 250.0
    depth_m = 760.0
    df = 1.0
    fc = corner_frequency_patton(yield_kt)
    Rc = cavity_radius_empirical(yield_kt)
    zones = zone_radii(Rc)
    M0_c = scalar_moment_coupled(yield_kt)

    # (a) Cavity cross-section (760 m depth, large damage zones)
    ax = fig.add_subplot(gs[0, 0])
    colors = {"damaged": "#CCCC00", "fractured": "#FF8800",
              "crushed": "#CC3311", "cavity": "#111"}
    extent = 1.3 * zones["damaged"]
    for name in ["damaged", "fractured", "crushed", "cavity"]:
        ax.add_patch(Circle((0, -depth_m), zones[name],
                            fc=colors[name], ec="k", lw=1.2))
    ax.axhline(0, color="brown", lw=3)
    ax.fill_between([-extent, extent], 0, 80, color="#E8D8B8", alpha=0.4)
    ax.plot(0, -depth_m, "w*", ms=18, mec="k", mew=1.5, zorder=10)
    ax.plot(0, 0, "rv", ms=12, zorder=10)

    ax.annotate("Mt. Mantap surface", xy=(0, 0), xytext=(extent*0.3, 150),
                fontsize=7, color="brown", arrowprops=dict(arrowstyle="->", color="brown"))

    ax.set_xlim(-extent, extent); ax.set_ylim(-depth_m - 1.1*zones["damaged"], 200)
    ax.set_aspect("equal"); ax.grid(True, alpha=0.2)
    ax.set_xlabel("Distance (m)"); ax.set_ylabel("Depth (m)")
    ax.set_title("(a)  Cavity & Damage Zones (760 m depth)", fontweight="bold")
    from matplotlib.patches import Patch
    ax.legend([Patch(fc=c, ec="k") for c in colors.values()],
              [f"Damaged ({zones['damaged']:.0f} m)",
               f"Fractured ({zones['fractured']:.0f} m)",
               f"Crushed ({zones['crushed']:.0f} m)",
               f"Cavity ({zones['cavity']:.0f} m)"],
              fontsize=7, loc="upper right")

    # (b) RDP time function
    ax = fig.add_subplot(gs[0, 1])
    dur = 12.0 / fc
    tt = np.linspace(-0.05 * dur, dur, 600)
    psi, dpsi = mueller_murphy_time(tt, yield_kt, fc=fc)
    psi_n = psi / np.max(psi) if np.max(psi) > 0 else psi
    dpsi_n = dpsi / np.max(np.abs(dpsi)) if np.max(np.abs(dpsi)) > 0 else dpsi
    ax.plot(tt, psi_n, "b-", lw=2, label="ψ(t)")
    ax.plot(tt, dpsi_n, "r--", lw=1.5, label="dψ/dt")
    ax.axhline(0, color="gray", lw=0.5); ax.axvline(0, color="gray", lw=0.5)
    ax.set_xlabel("Time (s)"); ax.set_ylabel("Normalised amplitude")
    ax.set_title("(b)  Reduced Displacement Potential", fontweight="bold")
    ax.legend(fontsize=9); ax.grid(True, alpha=0.2)

    # (c) Source spectrum — coupled (no decoupled comparison for this event)
    ax = fig.add_subplot(gs[0, 2])
    f = np.logspace(-2, 1.5, 500)
    spec_c = mueller_murphy_spectrum(f, yield_kt, fc=fc)
    ax.loglog(f, spec_c / spec_c[1], "b-", lw=2.5, label=f"Coupled ({yield_kt:.0f} kt)")
    ax.axvline(fc, color="green", ls="--", lw=1.5, label=f"f_c = {fc:.2f} Hz")

    spec_1kt = mueller_murphy_spectrum(f, 1.0)
    ax.loglog(f, spec_1kt / spec_c[1], "gray", lw=1.2, ls="--", alpha=0.6,
              label="1 kt (reference)")

    ax.set_xlabel("Frequency (Hz)"); ax.set_ylabel("Relative amplitude")
    ax.set_title("(c)  Source Spectrum (Coupled)", fontweight="bold")
    ax.legend(fontsize=9); ax.grid(True, alpha=0.2, which="both")
    ax.set_xlim(0.01, 20)

    # (d) mb–yield diagram
    ax = fig.add_subplot(gs[1, 0])
    W = np.logspace(-1, 4, 200)
    mb_c = 4.0 + 0.75 * np.log10(W)
    ax.plot(W, mb_c, "b-", lw=2.2, label="Coupled (DF=1)")
    for dfi, col, ls in [(10, "orange", "--"), (70, "r", "--"), (200, "purple", ":")]:
        mb_d = 4.0 + 0.75 * np.log10(W / dfi)
        ax.plot(W, mb_d, color=col, ls=ls, lw=1.4, label=f"DF = {dfi}")

    mb_this = mb_from_yield(yield_kt, df)
    ax.plot(yield_kt, mb_this, "r*", ms=22, zorder=10,
            label=f"DPRK 2017 ({yield_kt:.0f} kt)\nmb_pred = {mb_this:.2f}")
    ax.axhline(6.3, color="gray", ls=":", lw=1.5, label="Observed mb ≈ 6.3")
    ax.set_xscale("log"); ax.set_xlabel("Yield (kt)"); ax.set_ylabel("mb")
    ax.set_title("(d)  mb–Yield Relationship", fontweight="bold")
    ax.legend(fontsize=7.5, loc="upper left"); ax.grid(True, alpha=0.2, which="both")
    ax.set_xlim(0.1, 1e4); ax.set_ylim(2, 8)

    # (e) Damage profile
    ax = fig.add_subplot(gs[1, 1])
    r = np.logspace(np.log10(zones["cavity"]*0.8), np.log10(zones["damaged"]*3), 300)
    D = np.zeros_like(r)
    for j, ri in enumerate(r):
        if ri < zones["cavity"]:
            D[j] = 1.0
        elif ri < zones["crushed"]:
            D[j] = 0.9 + 0.09*(1-(ri-zones["cavity"])/(zones["crushed"]-zones["cavity"]))
        elif ri < zones["fractured"]:
            D[j] = 0.4 + 0.5*(1-(ri-zones["crushed"])/(zones["fractured"]-zones["crushed"]))
        elif ri < zones["damaged"]:
            D[j] = 0.1 + 0.3*(1-(ri-zones["fractured"])/(zones["damaged"]-zones["fractured"]))
    ax.semilogx(r, D, "r-", lw=2)
    ax.fill_between(r, 0, D, color="red", alpha=0.15)
    for nm, rv, c in [("Cav", zones["cavity"], "k"), ("Cr", zones["crushed"], "#CC3311"),
                       ("Fr", zones["fractured"], "#FF8800"), ("Dm", zones["damaged"], "#CCCC00")]:
        ax.axvline(rv, color=c, ls="--", lw=1.2)
        ax.text(rv, 1.02, nm, ha="center", fontsize=8, color=c, fontweight="bold")
    ax.set_xlabel("Radius (m)"); ax.set_ylabel("Damage D")
    ax.set_title("(e)  Radial Damage Profile", fontweight="bold")
    ax.set_ylim(0, 1.1); ax.grid(True, alpha=0.2)

    mantap_note = (f"NOTE: Damage zone radius {zones['damaged']:.0f} m\n"
                   f"exceeds overburden depth {depth_m:.0f} m.\n"
                   "→ Mt. Mantap surface collapse observed.")
    ax.text(0.98, 0.55, mantap_note, transform=ax.transAxes, fontsize=7,
            ha="right", va="center", color="#CC3311", style="italic",
            bbox=dict(boxstyle="round,pad=0.3", fc="#fff8f0", ec="#CC3311", alpha=0.8))

    # (f) Parameter table
    ax = fig.add_subplot(gs[1, 2])
    ax.axis("off")
    lines = [
        "EVENT & SOURCE PARAMETERS",
        "═" * 38,
        f"Date/time:          2017-09-03  03:30:01 UTC",
        f"Location:           {EVENT_LAT:.3f}°N  {EVENT_LON:.3f}°E",
        f"Depth:              ~{depth_m:.0f} m  (under Mt. Mantap)",
        "",
        f"Yield (estimated):     {yield_kt:.0f} kt",
        f"Decoupling factor:     {df:.0f}  (fully coupled)",
        f"Cavity radius:         {Rc:.0f} m  (empirical)",
        "",
        f"Scalar moment (coupled):    {M0_c:.2e} N·m",
        f"Mw (coupled):    {Mw_from_M0(M0_c):.2f}",
        f"mb (predicted):  {mb_this:.2f}",
        f"mb (observed):   ~6.3",
        "",
        f"Corner frequency:  {fc:.2f} Hz  (Patton)",
        f"Source model:       Mueller-Murphy (1971)",
        "",
        "SPECIAL FEATURES:",
        "  • Mountain collapse at Mt. Mantap",
        "  • Landslide signals observed 8.5 s after",
        "  • Surface subsidence visible in InSAR",
    ]
    ax.text(0.05, 0.97, "\n".join(lines), transform=ax.transAxes,
            fontfamily="monospace", fontsize=8.5, va="top",
            bbox=dict(boxstyle="round,pad=0.4", fc="#f8f8f8", ec="#bbb"))

    fig.suptitle("Figure 1 — Source Model:  DPRK 6th Nuclear Test (250 kt, Coupled)",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    p = os.path.join(OUTDIR, "fig01_source_model.png")
    fig.savefig(p, dpi=180, bbox_inches="tight"); plt.close(fig)
    print(f"  Saved {p}")


def fig02_velocity_model():
    """Korean Peninsula 1-D velocity model + station map."""
    model = PunggyeRiVelocityModel()
    z, vp, vs, rho, qp, qs = model.profile_arrays(dz=0.05, zmax=50)

    fig, axes = plt.subplots(1, 4, figsize=(20, 8))

    # (a) Vp / Vs
    ax = axes[0]
    ax.plot(vp, z, "b-", lw=2, label="Vp")
    ax.plot(vs, z, "r-", lw=2, label="Vs")
    ax.axhline(33, color="green", ls="--", lw=1.5, label="Moho (33 km)")
    ax.invert_yaxis(); ax.set_xlabel("Velocity (km/s)"); ax.set_ylabel("Depth (km)")
    ax.set_title("(a) P- and S-wave velocity", fontweight="bold")
    ax.legend(fontsize=9); ax.grid(True, alpha=0.2)

    # (b) Density
    ax = axes[1]
    ax.plot(rho, z, "k-", lw=2)
    ax.axhline(33, color="green", ls="--", lw=1.5)
    ax.invert_yaxis(); ax.set_xlabel("Density (kg/m³)"); ax.set_ylabel("Depth (km)")
    ax.set_title("(b) Density profile", fontweight="bold")
    ax.grid(True, alpha=0.2)

    # (c) Q model
    ax = axes[2]
    ax.plot(qp, z, "b-", lw=2, label="Qp")
    ax.plot(qs, z, "r-", lw=2, label="Qs")
    ax.axhline(33, color="green", ls="--", lw=1.5)
    ax.invert_yaxis(); ax.set_xlabel("Quality factor Q"); ax.set_ylabel("Depth (km)")
    ax.set_title("(c) Attenuation model", fontweight="bold")
    ax.legend(fontsize=9); ax.grid(True, alpha=0.2)

    # (d) Station map
    ax = axes[3]
    ax.plot(EVENT_LON, EVENT_LAT, "r*", ms=18, zorder=10, label="Punggye-ri (source)")
    for net, sta, dist, az, desc in STATION_INFO:
        for n2, s2, slat, slon, d2 in STATIONS:
            if s2 == sta:
                ax.plot(slon, slat, "^", ms=8, color="#1f77b4", zorder=5)
                ax.text(slon + 0.4, slat + 0.3, f"{sta}\n({dist:.0f} km)",
                        fontsize=6, color="#333")
                break
    ax.set_xlabel("Longitude (°E)"); ax.set_ylabel("Latitude (°N)")
    ax.set_title("(d) Station map", fontweight="bold")
    ax.set_xlim(70, 150); ax.set_ylim(28, 55)
    ax.grid(True, alpha=0.2)
    ax.legend(fontsize=9, loc="lower left")
    ax.set_aspect(1.0 / np.cos(np.radians(42)))

    fig.suptitle("Figure 2 — Korean Peninsula: 1-D Velocity Model & Station Network",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    p = os.path.join(OUTDIR, "fig02_velocity_model.png")
    fig.savefig(p, dpi=180, bbox_inches="tight"); plt.close(fig)
    print(f"  Saved {p}")


def fig03_synthetic_seismograms():
    """Record section of synthetic seismograms at all 12 stations."""
    model = PunggyeRiVelocityModel()
    yield_kt = 250.0
    df = 1.0

    fig = plt.figure(figsize=(20, 14))
    gs = GridSpec(1, 24, figure=fig)
    ax = fig.add_subplot(gs[0, :20])
    ax_info = fig.add_subplot(gs[0, 20:]); ax_info.axis("off")

    n = len(STATION_INFO)
    y_pos = np.arange(n)

    phase_colors = {"Pn": "#1f77b4", "Pg": "#2ca02c", "Sn": "#d62728", "Lg": "#ff7f0e"}

    for i, (net, sta, dist_km, az, desc) in enumerate(STATION_INFO):
        t, v = generate_synthetic(dist_km, yield_kt, df, dt=0.05, duration=500, model=model)

        peak = np.max(np.abs(v))
        if peak > 0:
            v = v / peak * 0.38

        ax.fill_between(t, v + y_pos[i], y_pos[i],
                        where=(v > 0), color="#444", alpha=0.12)
        ax.fill_between(t, v + y_pos[i], y_pos[i],
                        where=(v < 0), color="#444", alpha=0.06)
        ax.plot(t, v + y_pos[i], "k-", lw=0.4, alpha=0.85)

        ax.text(-8, y_pos[i], f"{net}.{sta}  ({dist_km:.0f} km)",
                fontsize=8, fontweight="bold", ha="right", va="center",
                color="#222", clip_on=False)

        tt = regional_travel_times(dist_km, model)
        for ph, (tarr, _, _) in tt.items():
            if 0 < tarr < 500:
                ax.plot(tarr, y_pos[i], "d", color=phase_colors[ph], ms=5, alpha=0.6, zorder=3)

    for ph in phase_colors:
        ax.plot([], [], "d", color=phase_colors[ph], ms=6, label=f"{ph}")
    ax.legend(loc="lower right", fontsize=9, ncol=4, title="Phases")

    ax.set_xlabel("Time (s)", fontsize=12, fontweight="bold")
    ax.set_ylabel("Station (sorted by distance →)", fontsize=12, fontweight="bold")
    ax.set_ylim(-0.8, n - 0.2); ax.set_xlim(-10, 480)
    ax.set_yticks(y_pos); ax.set_yticklabels([])
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.grid(True, which="major", axis="x", alpha=0.2)
    for y in y_pos:
        ax.axhline(y, color="#ddd", lw=0.3)

    fc = corner_frequency_patton(yield_kt)
    info = [
        "Synthetic Parameters",
        "─" * 26,
        f"Yield:    {yield_kt:.0f} kt (coupled)",
        f"DF:       {df:.0f} (fully coupled)",
        f"fc:       {fc:.2f} Hz",
        f"mb_pred:  {mb_from_yield(yield_kt, df):.2f}",
        f"Filter:   {FMIN}–{FMAX} Hz",
        f"Model:    Mueller-Murphy",
        f"Velocity: Punggye-ri 1-D",
        f"M0:       {scalar_moment_coupled(yield_kt):.2e} N·m",
    ]
    ax_info.text(0.05, 0.98, "\n".join(info), transform=ax_info.transAxes,
                 fontfamily="monospace", fontsize=8, va="top",
                 bbox=dict(boxstyle="round,pad=0.4", fc="#f8f8f8", ec="#bbb"))

    fig.suptitle(f"Figure 3 — Synthetic Seismograms:  {yield_kt:.0f} kt Coupled Test (DF=1)\n"
                 f"Mueller-Murphy source, Punggye-ri 1-D model, BHZ {FMIN}–{FMAX} Hz",
                 fontsize=13, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    p = os.path.join(OUTDIR, "fig03_synthetic_seismograms.png")
    fig.savefig(p, dpi=180, bbox_inches="tight"); plt.close(fig)
    print(f"  Saved {p}")


def fig04_yield_sensitivity():
    """Synthetic waveforms at MDJ (~370 km) for different yields."""
    dist_km = 370.0
    model = PunggyeRiVelocityModel()
    yields = [10, 50, 100, 250, 500]
    df = 1.0

    fig, axes = plt.subplots(len(yields), 1, figsize=(16, 3.0 * len(yields)),
                             sharex=True)

    for i, W in enumerate(yields):
        t, v = generate_synthetic(dist_km, W, df, dt=0.05, duration=400, model=model)
        mb = mb_from_yield(W, df)

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
                ax.axvline(tarr, color={"Pn":"#1f77b4","Pg":"#2ca02c",
                                         "Sn":"#d62728","Lg":"#ff7f0e"}[ph],
                           ls="--", lw=1, alpha=0.6)

    axes[-1].set_xlabel("Time after origin (s)", fontsize=11)

    fig.suptitle("Figure 4 — Yield Sensitivity at MDJ (370 km):  Coupled Test (DF=1)\n"
                 f"Punggye-ri 1-D model, BHZ {FMIN}–{FMAX} Hz",
                 fontsize=13, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    p = os.path.join(OUTDIR, "fig04_yield_sensitivity.png")
    fig.savefig(p, dpi=180, bbox_inches="tight"); plt.close(fig)
    print(f"  Saved {p}")


def fig05_yield_tradeoff():
    """Compare 100 kt vs 200 kt vs 300 kt vs 400 kt at MDJ (yield tradeoff)."""
    dist_km = 370.0
    model = PunggyeRiVelocityModel()
    df = 1.0
    cases = [
        (100.0, "100 kt",  "#2ca02c"),
        (200.0, "200 kt",  "#1f77b4"),
        (300.0, "300 kt",  "#ff7f0e"),
        (400.0, "400 kt",  "#d62728"),
    ]

    fig, axes = plt.subplots(len(cases) + 1, 1,
                             figsize=(16, 3.5 * (len(cases) + 1)),
                             gridspec_kw={"height_ratios": [1]*len(cases) + [1.2]})

    spectra = {}
    ref_peak = None
    for i, (W, label, col) in enumerate(cases):
        t, v = generate_synthetic(dist_km, W, df, dt=0.05, duration=400, model=model)

        if ref_peak is None:
            ref_peak = np.max(np.abs(v))

        ax = axes[i]
        peak = np.max(np.abs(v))
        v_rel = v / ref_peak if ref_peak > 0 else v
        ax.plot(t, v_rel, "-", color=col, lw=0.5, alpha=0.8)

        mb = mb_from_yield(W, df)
        fc = corner_frequency_patton(W)
        ax.set_ylabel("Rel. vel.", fontsize=9)
        ax.set_title(f"{label}  |  mb = {mb:.2f}  |  fc = {fc:.2f} Hz  |  "
                     f"rel. amplitude = {peak/ref_peak:.3f}",
                     fontsize=9, fontweight="bold", loc="left")
        ax.grid(True, alpha=0.2)
        if i < len(cases) - 1:
            ax.set_xticklabels([])

        freq, spec = spectral_amplitude(v, 0.05)
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
    ax.set_title("Spectral comparison — different yields at MDJ",
                 fontsize=9, fontweight="bold", loc="left")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.2, which="both")

    fig.suptitle("Figure 5 — Yield Trade-off at MDJ (370 km):  Coupled Tests\n"
                 f"Amplitudes relative to 100 kt, {FMIN}–{FMAX} Hz",
                 fontsize=13, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    p = os.path.join(OUTDIR, "fig05_yield_tradeoff.png")
    fig.savefig(p, dpi=180, bbox_inches="tight"); plt.close(fig)
    print(f"  Saved {p}")


def fig06_real_data(real_data):
    """Display real observed seismograms from IRIS — DPRK 2017 is well-recorded."""
    n = len(real_data)
    if n == 0:
        print("  No real data — skipping fig06")
        return

    fig = plt.figure(figsize=(20, max(12, 1.5 * n)))
    gs = GridSpec(1, 24, figure=fig)
    ax = fig.add_subplot(gs[0, :20])
    ax_info = fig.add_subplot(gs[0, 20:]); ax_info.axis("off")

    y_pos = np.arange(n)
    TMIN, TMAX = -30, 420

    phase_colors = {"Pn": "#1f77b4", "Pg": "#2ca02c", "Sn": "#d62728", "Lg": "#ff7f0e"}

    for i, (net, sta, dist_km, az, tr, desc) in enumerate(real_data):
        times = _time_relative(tr, EVENT_TIME)
        data = tr.data.copy()

        mask = (times >= TMIN) & (times <= TMAX)
        t_w = times[mask]; d_w = data[mask]
        peak = np.max(np.abs(d_w)) if len(d_w) > 0 else 1.0
        if peak > 0:
            d_w = d_w / peak * 0.38

        ax.fill_between(t_w, d_w + y_pos[i], y_pos[i],
                        where=(d_w > 0), color="#444", alpha=0.12)
        ax.fill_between(t_w, d_w + y_pos[i], y_pos[i],
                        where=(d_w < 0), color="#444", alpha=0.06)
        ax.plot(t_w, d_w + y_pos[i], "k-", lw=0.4, alpha=0.85)

        ax.text(TMIN - 4, y_pos[i],
                f"{net}.{sta}  ({dist_km:.0f} km)",
                fontsize=8, fontweight="bold", ha="right", va="center",
                color="#222", clip_on=False)

        tt = regional_travel_times(dist_km)
        for ph, (tarr, _, _) in tt.items():
            if TMIN < tarr < TMAX:
                ax.plot(tarr, y_pos[i], "d", color=phase_colors[ph], ms=5, alpha=0.6, zorder=3)

    for ph in phase_colors:
        ax.plot([], [], "d", color=phase_colors[ph], ms=6, label=ph)
    ax.legend(loc="lower right", fontsize=9, ncol=4)

    ax.set_xlabel("Time after origin (s)", fontsize=12, fontweight="bold")
    ax.set_ylabel("Station (sorted by distance →)", fontsize=12, fontweight="bold")
    ax.set_ylim(-0.8, n - 0.2); ax.set_xlim(TMIN, TMAX)
    ax.set_yticks(y_pos); ax.set_yticklabels([])
    ax.grid(True, which="major", axis="x", alpha=0.2)
    for y in y_pos:
        ax.axhline(y, color="#ddd", lw=0.3)

    ax_info.text(0.05, 0.98, "Real Data\n" + "─"*24 + "\n"
                 f"Source: IRIS FDSN\n"
                 f"Channel: BHZ\n"
                 f"Filter: {FMIN}–{FMAX} Hz\n"
                 f"Resp: removed → vel\n"
                 f"Window: signal only\n"
                 f"Norm: per-trace peak\n"
                 f"\nEvent: DPRK-6\n"
                 f"mb ≈ 6.3 (very well\n"
                 f"recorded globally)",
                 transform=ax_info.transAxes, fontfamily="monospace",
                 fontsize=8, va="top",
                 bbox=dict(boxstyle="round,pad=0.4", fc="#f8f8f8", ec="#bbb"))

    fig.suptitle("Figure 6 — Observed Waveforms (zoomed signal window)\n"
                 f"DPRK 6th Nuclear Test 2017-09-03 03:30:01 UTC, BHZ {FMIN}–{FMAX} Hz",
                 fontsize=13, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    p = os.path.join(OUTDIR, "fig06_real_data.png")
    fig.savefig(p, dpi=180, bbox_inches="tight"); plt.close(fig)
    print(f"  Saved {p}")


def fig07_comparison(real_data):
    """Side-by-side synthetic vs real for key stations."""
    model = PunggyeRiVelocityModel()
    yield_kt = 250.0; df = 1.0

    key_stas = ["MDJ", "INCN", "KSRS", "MAJO", "BJT"]
    entries = [(n, s, d, a, tr, desc) for n, s, d, a, tr, desc in real_data
               if s in key_stas]

    if not entries:
        entries = real_data[:5]

    n = len(entries)
    if n == 0:
        print("  No data for comparison — skipping fig07")
        return

    fig, axes = plt.subplots(n, 2, figsize=(20, 3.5 * n), squeeze=False)

    for i, (net, sta, dist_km, az, tr, desc) in enumerate(entries):
        times_r = _time_relative(tr, EVENT_TIME)
        data_r = tr.data.copy()
        TMIN = -20; TMAX = 400

        mask = (times_r >= TMIN) & (times_r <= TMAX)
        t_r = times_r[mask]; d_r = data_r[mask]
        peak_r = np.max(np.abs(d_r)) if len(d_r) > 0 else 1.0
        if peak_r > 0: d_r = d_r / peak_r

        t_s, v_s = generate_synthetic(dist_km, yield_kt, df, dt=0.05, duration=420, model=model)
        mask_s = (t_s >= TMIN) & (t_s <= TMAX)
        t_s = t_s[mask_s]; v_s = v_s[mask_s]
        peak_s = np.max(np.abs(v_s)) if len(v_s) > 0 else 1.0
        if peak_s > 0: v_s = v_s / peak_s

        ax = axes[i][0]
        ax.plot(t_r, d_r, "k-", lw=0.4, alpha=0.7, label="Observed")
        ax.plot(t_s, v_s, "r-", lw=0.6, alpha=0.7, label="Synthetic")

        tt = regional_travel_times(dist_km, model)
        for ph, (tarr, _, _) in tt.items():
            if TMIN < tarr < TMAX:
                ax.axvline(tarr, color={"Pn":"#1f77b4","Pg":"#2ca02c",
                                         "Sn":"#d62728","Lg":"#ff7f0e"}[ph],
                           ls="--", lw=0.8, alpha=0.5)

        ax.set_xlim(TMIN, TMAX); ax.set_ylim(-1.2, 1.2)
        ax.set_title(f"{net}.{sta} ({dist_km:.0f} km) — Waveforms",
                     fontsize=9, fontweight="bold")
        ax.legend(fontsize=7, loc="upper right")
        ax.grid(True, alpha=0.2)
        if i == n - 1:
            ax.set_xlabel("Time after origin (s)")

        ax2 = axes[i][1]
        freq_r, sp_r = spectral_amplitude(d_r, tr.stats.delta)
        freq_s, sp_s = spectral_amplitude(v_s, 0.05)

        m_r = (freq_r > 0.1) & (freq_r < 10)
        m_s = (freq_s > 0.1) & (freq_s < 10)

        sp_r_n = sp_r[m_r] / np.max(sp_r[m_r]) if np.max(sp_r[m_r]) > 0 else sp_r[m_r]
        sp_s_n = sp_s[m_s] / np.max(sp_s[m_s]) if np.max(sp_s[m_s]) > 0 else sp_s[m_s]

        ax2.loglog(freq_r[m_r], sp_r_n, "k-", lw=1, alpha=0.7, label="Observed")
        ax2.loglog(freq_s[m_s], sp_s_n, "r-", lw=1, alpha=0.7, label="Synthetic")
        ax2.set_title(f"{net}.{sta} — Spectra", fontsize=9, fontweight="bold")
        ax2.legend(fontsize=7); ax2.grid(True, alpha=0.2, which="both")
        ax2.set_xlim(0.1, 10)
        if i == n - 1:
            ax2.set_xlabel("Frequency (Hz)")

    fig.suptitle("Figure 7 — Synthetic vs Observed Comparison\n"
                 f"250 kt coupled (DF=1), Mueller-Murphy source, {FMIN}–{FMAX} Hz",
                 fontsize=13, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    p = os.path.join(OUTDIR, "fig07_comparison.png")
    fig.savefig(p, dpi=180, bbox_inches="tight"); plt.close(fig)
    print(f"  Saved {p}")


def fig08_discrimination(real_data):
    """Explosion vs earthquake discrimination analysis — DPRK 2017 is clearly an explosion."""
    fig = plt.figure(figsize=(18, 12))
    gs = GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.35)

    model = PunggyeRiVelocityModel()
    yield_kt = 250.0; df = 1.0

    # (a) P/Lg amplitude ratio for synthetic
    ax = fig.add_subplot(gs[0, 0])
    distances = np.linspace(200, 2000, 20)
    ps_ratios = []
    for d in distances:
        t, v = generate_synthetic(d, yield_kt, df, dt=0.05, duration=600, model=model)
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

    mb_this = 6.3
    ms_this = 4.5
    ax.plot(mb_this, ms_this, "r*", ms=22, zorder=10,
            label=f"DPRK 2017 (mb={mb_this:.1f}, Ms={ms_this:.1f})")

    mb_line = np.linspace(2, 8, 100)
    ms_line = mb_line - 1.0
    ax.plot(mb_line, ms_line, "k--", lw=1.5, label="Discrimination line")

    ax.set_xlabel("mb"); ax.set_ylabel("Ms")
    ax.set_title("(b) mb–Ms Discrimination", fontweight="bold")
    ax.legend(fontsize=7, loc="upper left"); ax.grid(True, alpha=0.2)
    ax.set_xlim(2, 8); ax.set_ylim(1, 8)

    # (c) Spectral ratio analysis (synthetic)
    ax = fig.add_subplot(gs[0, 2])
    dist_mdj = 370.0
    t, v = generate_synthetic(dist_mdj, yield_kt, df, dt=0.05, duration=500, model=model)
    freq, spec = spectral_amplitude(v, 0.05)
    m = (freq > 0.3) & (freq < 10)
    spec_n = spec[m] / np.max(spec[m]) if np.max(spec[m]) > 0 else spec[m]

    ax.semilogy(freq[m], spec_n, "r-", lw=1.5, label="Explosion (synth, 250 kt)")

    fc_eq = 0.3
    spec_eq = 1.0 / (1 + (freq[m] / fc_eq)**2)**0.5
    spec_eq_n = spec_eq / np.max(spec_eq)
    ax.semilogy(freq[m], spec_eq_n, "b--", lw=1.5, label="Earthquake (model, same mb)")

    ax.set_xlabel("Frequency (Hz)"); ax.set_ylabel("Normalised amplitude")
    ax.set_title("(c) Spectral Shape at MDJ (370 km)", fontweight="bold")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.2, which="both")

    # (d) Observed spectral ratio at MDJ (if real data available)
    ax = fig.add_subplot(gs[1, 0])
    mdj_data = [x for x in real_data if x[1] in ("MDJ",)]
    if mdj_data:
        net, sta, dist_km, az, tr, desc = mdj_data[0]
        times = _time_relative(tr, EVENT_TIME)
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
        ax.text(0.5, 0.5, "No MDJ data available", transform=ax.transAxes,
                ha="center", fontsize=12, color="gray")

    ax.set_xlabel("Frequency (Hz)"); ax.set_ylabel("Normalised amplitude")
    ax.set_title("(d) Observed P vs Lg Spectra at MDJ", fontweight="bold")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.2, which="both")

    # (e) Discrimination summary
    ax = fig.add_subplot(gs[1, 1])
    ax.axis("off")
    mb_pred = mb_from_yield(yield_kt, df)
    lines = [
        "DISCRIMINATION ANALYSIS",
        "═" * 36,
        "",
        "INDICATORS CONSISTENT WITH EXPLOSION:",
        "",
        "  ✓  Impulsive P, high P/Lg ratio",
        "     at all regional distances",
        "",
        "  ✓  mb–Ms well below earthquake trend",
        f"     (mb≈{mb_this}, Ms≈{ms_this})",
        "",
        "  ✓  Higher-frequency spectral content",
        "     than equivalent earthquake",
        "",
        "  ✓  Location: Punggye-ri test site",
        "     (known DPRK nuclear facility)",
        "",
        "  ✓  Shallow depth (~760 m)",
        "",
        "  ✓  Mountain collapse signal observed",
        "     ~8.5 s after main event at",
        "     Mt. Mantap (secondary source)",
        "",
        "COUPLED TEST PARAMETERS:",
        "",
        f"  Yield:  ~{yield_kt:.0f} kt",
        f"  DF:     1 (fully coupled)",
        f"  mb:     {mb_pred:.2f} (pred) / {mb_this} (obs)",
        "",
        "  → CLEARLY an explosion.",
        "  → Largest DPRK test to date.",
    ]
    ax.text(0.05, 0.97, "\n".join(lines), transform=ax.transAxes,
            fontfamily="monospace", fontsize=8.5, va="top",
            bbox=dict(boxstyle="round,pad=0.5", fc="#fffff0", ec="#cca"))

    # (f) Yield-mb plot with observed mb constraint
    ax = fig.add_subplot(gs[1, 2])
    W_range = np.logspace(0, 4, 200)
    mbs_coupled = 4.0 + 0.75 * np.log10(W_range)
    ax.semilogx(W_range, mbs_coupled, "b-", lw=2.5, label="Coupled (DF=1)")

    for dfi in [5, 10, 30, 70]:
        mbs = 4.0 + 0.75 * np.log10(W_range / dfi)
        ax.semilogx(W_range, mbs, lw=1.2, ls="--", alpha=0.5, label=f"DF={dfi}")

    ax.axhline(6.3, color="gray", ls=":", lw=2, label="Observed mb ≈ 6.3")
    ax.fill_between(W_range, 6.1, 6.5, color="gray", alpha=0.1)

    ax.plot(250, mb_from_yield(250, 1.0), "r*", ms=22, zorder=10,
            label=f"Best estimate: 250 kt")

    ax.set_xlabel("True yield (kt)"); ax.set_ylabel("Observed mb")
    ax.set_title("(f) Yield–mb Constraint", fontweight="bold")
    ax.legend(fontsize=7, ncol=2, loc="upper left"); ax.grid(True, alpha=0.2, which="both")
    ax.set_xlim(1, 1e4); ax.set_ylim(3, 8)

    fig.suptitle("Figure 8 — Explosion/Earthquake Discrimination Analysis\n"
                 "DPRK 6th Nuclear Test 2017-09-03 03:30:01 UTC",
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

    print("=" * 72)
    print("  DPRK 2017 — 6th Nuclear Test (250 kt Coupled) Model & Analysis")
    print("=" * 72)
    print()

    print("[1/8] Source model ...")
    fig01_source_model()

    print("[2/8] Velocity model ...")
    fig02_velocity_model()

    print("[3/8] Synthetic seismograms ...")
    fig03_synthetic_seismograms()

    print("[4/8] Yield sensitivity ...")
    fig04_yield_sensitivity()

    print("[5/8] Yield trade-off ...")
    fig05_yield_tradeoff()

    print("[6/8] Fetching real data from IRIS ...")
    real_data = fetch_real_waveforms(pre=120, post=500)
    print(f"       Retrieved {len(real_data)} station(s)")
    fig06_real_data(real_data)

    print("[7/8] Synthetic vs observed comparison ...")
    fig07_comparison(real_data)

    print("[8/8] Discrimination analysis ...")
    fig08_discrimination(real_data)

    print()
    print("=" * 72)
    print(f"  All figures saved to: {OUTDIR}")
    print("=" * 72)


if __name__ == "__main__":
    main()
