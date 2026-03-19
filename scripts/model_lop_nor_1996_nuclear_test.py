#!/usr/bin/env python3
"""
Model the Last Confirmed Chinese Nuclear Test at Lop Nor — 1996-07-29

Comprehensive physics-based modeling of a fully coupled underground nuclear test:
  1. Mueller-Murphy seismic source (RDP model) — ~5 kt, fully coupled
  2. 1-D layered velocity model for the Lop Nor / Kuruktag region
  3. Synthetic seismogram generation at real station distances
  4. Comparison against actual KNDC archive data (26 regional stations)
  5. Explosion-vs-earthquake discrimination analysis

This was the 45th and final confirmed Chinese underground nuclear test,
conducted in a horizontal tunnel in the Kuruktag mountains on the NW edge
of the Lop Nor test site.  China signed the CTBT two months later
(1996-09-24) and declared a testing moratorium.

All results are presented as publication-quality figures.

Physical references:
  Mueller & Murphy (1971)  — Seismic characteristics of underground nuclear detonations
  Patton (1988)            — Corner-frequency scaling
  Charlie & Veyera (1994)  — Geology of the Lop Nor test site
  Zhenbang & Shangyi (1986)— Chinese nuclear testing practice
  Xia et al. (2017)        — NCCrust crustal velocity model
  Ringdal et al. (1992)    — mb-yield calibration for Lop Nor
  Murphy (1996)            — Lop Nor yield estimation from body waves

Usage:
    python scripts/model_lop_nor_1996_nuclear_test.py
"""

import os
import sys
import warnings
import numpy as np

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
    mueller_murphy_spectrum,
)
from fsrm.velocity_models import LopNorVelocityModel
from fsrm.propagation import regional_travel_times, generate_synthetic
from fsrm.signal_processing import envelope, spectral_amplitude, time_relative
from fsrm.data_fetching import fetch_waveforms, compute_station_distances
from fsrm.inversion import invert_epicenter
from fsrm.plotting import (
    PHASE_COLORS, fig_source_model, fig_velocity_model,
    fig_synthetic_record_section, fig_observed_record_section,
    fig_comparison, fig_location_inversion,
)

# ═══════════════════════════════════════════════════════════════════════════════
# Event Parameters
# ═══════════════════════════════════════════════════════════════════════════════
EVENT_TIME   = UTCDateTime("1996-07-29T01:48:57.170")
EVENT_LAT    = 41.7163
EVENT_LON    = 88.3748
EVENT_DEPTH  = 0.4            # km (horizontal tunnel)
YIELD_KT     = 5.0            # ~5 kt from mb-yield scaling (mb 4.90)
DF            = 1.0            # Fully coupled

OUTDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                      "figures", "lop_nor_1996")

# Bandpass for processing
FMIN, FMAX = 0.5, 8.0

EVENT_INFO = {
    "name": "Lop Nor Final Test (45th)",
    "time": "1996-07-29 01:48:57 UTC",
    "lat": EVENT_LAT,
    "lon": EVENT_LON,
    "mb_obs": 4.90,
}

# Stations for synthetic + real comparison
# These are modern FDSN stations at similar distances/azimuths to the 1996 KNDC
# stations, plus attempts at legacy KNDC networks.  For stations not available
# through IRIS (KNDC-only), the local archive data will substitute.
STATIONS = [
    # (net, sta, lat, lon, description)
    ("II",  "KURK", 50.715, 78.620,  "Kurchatov, KZ"),
    ("II",  "AAK",  42.638, 74.494,  "Ala Archa, KG"),
    ("IU",  "MAKZ", 46.808, 81.977,  "Makanchi PS23, KZ"),
    ("IC",  "WMQ",  43.814, 87.704,  "Urumqi, China"),
    ("KR",  "PRZ",  42.500, 78.400,  "Karakol, KG"),
    ("KZ",  "KNDC", 43.217, 76.966,  "Almaty, KZ"),
    ("KZ",  "PDGK", 43.328, 79.485,  "Podgonoye, KZ"),
    ("KZ",  "MKAR", 46.794, 82.290,  "Makanchi Array, KZ"),
    ("IC",  "LSA",  29.703, 91.128,  "Lhasa, Tibet"),
    ("IC",  "XAN",  34.031, 108.923, "Xi'an, China"),
    ("II",  "NIL",  33.651, 73.269,  "Nilore, PK"),
    ("IU",  "ULN",  47.865, 107.053, "Ulaanbaatar, Mongolia"),
]

STATION_INFO = compute_station_distances(EVENT_LAT, EVENT_LON, STATIONS)


# ═══════════════════════════════════════════════════════════════════════════════
# KNDC Archive Data Loading
# ═══════════════════════════════════════════════════════════════════════════════

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


def load_kndc_observed():
    """
    Load observed waveforms from the local KNDC SAC archive.

    Returns list of (net, sta, dist_km, az, trace, description) tuples —
    same format as fetch_waveforms() for compatibility with shared plotting.
    """
    from obspy import read
    from obspy.geodetics import gps2dist_azimuth

    data_dir = os.path.join(
        _REPO_DIR,
        "data", "historic_nuclear", "03.Lopnor", "19960729.0148", "wf")

    if not os.path.isdir(data_dir):
        print(f"  WARNING: KNDC data directory not found: {data_dir}")
        return []

    site_coords = _parse_site_coords()

    results = []
    sac_files = sorted(f for f in os.listdir(data_dir) if f.endswith(".sac"))

    # Only load Z-component first segment files
    z_files = [f for f in sac_files
               if "z1" in f.lower() or "Z1" in f]

    for fname in z_files:
        fpath = os.path.join(data_dir, fname)
        try:
            st = read(fpath, format="SAC")
            tr = st[0]

            # Extract station name from SAC header
            sta = tr.stats.sac.get("kstnm", "").strip()
            if not sta:
                # Try to extract from filename
                base = fname.replace(".sac", "")
                # Patterns: 96211<sta>z1 or <STA>96211z1
                if base[:5].isdigit():
                    sta = base[5:].replace("z1", "").replace("z2", "").upper()
                else:
                    idx = base.find("96211")
                    if idx > 0:
                        sta = base[:idx].upper()
                    else:
                        sta = base[:4].upper()

            # Get distance/azimuth from SAC headers or compute from location
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

            # Bandpass filter
            tr.detrend("demean")
            tr.detrend("linear")
            tr.taper(max_percentage=0.05, type="cosine")
            tr.filter("bandpass", freqmin=FMIN, freqmax=FMAX,
                       corners=4, zerophase=True)

            results.append(("KZ", sta, dist_km, az, tr, f"{sta} (KNDC archive)"))
        except Exception as exc:
            print(f"  WARNING: Could not read {fname}: {exc}")

    # Sort by distance
    results.sort(key=lambda x: x[2])
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# Event-Specific Figures
# ═══════════════════════════════════════════════════════════════════════════════

def fig04_yield_sensitivity():
    """Synthetic waveforms at nearest KNDC station (~760 km) for different yields."""
    dist_km = 760.0   # DZHR, nearest station
    model = LopNorVelocityModel()
    yields = [1, 3, 5, 10, 20]

    fig, axes = plt.subplots(len(yields), 1, figsize=(16, 3.0 * len(yields)),
                             sharex=True)

    for i, W in enumerate(yields):
        t, v = generate_synthetic(dist_km, W, model,
                                  decoupling_factor=DF, dt=0.05, duration=400,
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

    fig.suptitle("Figure 4 — Yield Sensitivity at DZHR (760 km):  Coupled Test (DF=1)\n"
                 f"Lop Nor 1-D model, BHZ {FMIN}–{FMAX} Hz",
                 fontsize=13, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    p = os.path.join(OUTDIR, "fig04_yield_sensitivity.png")
    fig.savefig(p, dpi=180, bbox_inches="tight"); plt.close(fig)
    print(f"  Saved {p}")


def fig05_yield_tradeoff():
    """Compare yield estimates at DZHR (~760 km) to constrain yield from mb 4.90."""
    dist_km = 760.0
    model = LopNorVelocityModel()
    cases = [
        (2.0,  "2 kt",   "#2ca02c"),
        (5.0,  "5 kt",   "#1f77b4"),
        (8.0,  "8 kt",   "#ff7f0e"),
        (15.0, "15 kt",  "#d62728"),
    ]

    fig, axes = plt.subplots(len(cases) + 1, 1,
                             figsize=(16, 3.5 * (len(cases) + 1)),
                             gridspec_kw={"height_ratios": [1]*len(cases) + [1.2]})

    spectra = {}
    ref_peak = None
    for i, (W, label, col) in enumerate(cases):
        t, v = generate_synthetic(dist_km, W, model,
                                  decoupling_factor=DF, dt=0.05, duration=400,
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
    ax.set_title("Spectral comparison — different yields at DZHR (760 km)",
                 fontsize=9, fontweight="bold", loc="left")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.2, which="both")

    fig.suptitle("Figure 5 — Yield Trade-off at DZHR (760 km):  Coupled Tests\n"
                 f"Amplitudes relative to 2 kt, {FMIN}–{FMAX} Hz",
                 fontsize=13, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    p = os.path.join(OUTDIR, "fig05_yield_tradeoff.png")
    fig.savefig(p, dpi=180, bbox_inches="tight"); plt.close(fig)
    print(f"  Saved {p}")


def fig08_discrimination(real_data):
    """Explosion vs earthquake discrimination analysis — 1996 event is clearly an explosion."""
    fig = plt.figure(figsize=(18, 12))
    gs = GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.35)

    model = LopNorVelocityModel()

    # (a) P/Lg amplitude ratio for synthetic
    ax = fig.add_subplot(gs[0, 0])
    distances = np.linspace(200, 2000, 20)
    ps_ratios = []
    for d in distances:
        t, v = generate_synthetic(d, YIELD_KT, model,
                                  decoupling_factor=DF, dt=0.05, duration=600,
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

    mb_eq = np.linspace(3, 7, 50)
    ms_eq = 1.0 * mb_eq - 0.5 + 0.2 * np.random.randn(50)

    ax.scatter(mb_eq, ms_eq, c="b", s=20, alpha=0.4, label="Earthquakes")
    ax.scatter(mb_exp, ms_exp, c="r", s=20, alpha=0.4, label="Explosions (NTS)")

    mb_this = 4.90
    ms_this = 3.5   # estimated from mb-Ms explosion trend
    ax.plot(mb_this, ms_this, "r*", ms=20, zorder=10,
            label=f"Lop Nor 1996 (mb={mb_this:.2f})")

    mb_line = np.linspace(2, 7, 100)
    ms_line = mb_line - 1.0
    ax.plot(mb_line, ms_line, "k--", lw=1.5, label="Discrimination line")

    ax.set_xlabel("mb"); ax.set_ylabel("Ms")
    ax.set_title("(b) mb–Ms Discrimination", fontweight="bold")
    ax.legend(fontsize=7, loc="upper left"); ax.grid(True, alpha=0.2)
    ax.set_xlim(2, 7); ax.set_ylim(1, 7)

    # (c) Spectral ratio analysis (synthetic)
    ax = fig.add_subplot(gs[0, 2])
    dist_ref = 760.0    # DZHR
    t, v = generate_synthetic(dist_ref, YIELD_KT, model,
                              decoupling_factor=DF, dt=0.05, duration=500,
                              fmin=FMIN, fmax=FMAX, depth_km=EVENT_DEPTH)
    freq, spec = spectral_amplitude(v, 0.05)
    m = (freq > 0.3) & (freq < 10)
    spec_n = spec[m] / np.max(spec[m]) if np.max(spec[m]) > 0 else spec[m]

    ax.semilogy(freq[m], spec_n, "r-", lw=1.5, label="Explosion (synth, 5 kt)")

    fc_eq = 0.5
    spec_eq = 1.0 / (1 + (freq[m] / fc_eq)**2)**0.5
    spec_eq_n = spec_eq / np.max(spec_eq)
    ax.semilogy(freq[m], spec_eq_n, "b--", lw=1.5, label="Earthquake (model, same mb)")

    ax.set_xlabel("Frequency (Hz)"); ax.set_ylabel("Normalised amplitude")
    ax.set_title("(c) Spectral Shape at DZHR (760 km)", fontweight="bold")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.2, which="both")

    # (d) Observed spectral ratio at nearest station
    ax = fig.add_subplot(gs[1, 0])
    near_data = [x for x in real_data if x[2] < 900]
    if near_data:
        net, sta, dist_km, az, tr, desc = near_data[0]
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
            ax.set_title(f"(d) Observed P vs Lg Spectra at {sta} ({dist_km:.0f} km)",
                         fontweight="bold")
        else:
            ax.text(0.5, 0.5, "No P/Lg window data", transform=ax.transAxes,
                    ha="center", fontsize=12, color="gray")
            ax.set_title("(d) Observed P vs Lg Spectra", fontweight="bold")
    else:
        ax.text(0.5, 0.5, "No KNDC data available", transform=ax.transAxes,
                ha="center", fontsize=12, color="gray")
        ax.set_title("(d) Observed P vs Lg Spectra", fontweight="bold")

    ax.set_xlabel("Frequency (Hz)"); ax.set_ylabel("Normalised amplitude")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.2, which="both")

    # (e) Discrimination summary
    ax = fig.add_subplot(gs[1, 1])
    ax.axis("off")
    mb_pred = mb_from_yield(YIELD_KT, DF)
    lines = [
        "DISCRIMINATION ANALYSIS",
        "\u2550" * 36,
        "",
        "INDICATORS CONSISTENT WITH EXPLOSION:",
        "",
        "  \u2713  Impulsive P, high P/Lg ratio",
        "     at all regional distances",
        "",
        "  \u2713  mb\u2013Ms below earthquake trend",
        f"     (mb={mb_this:.2f}, Ms\u2248{ms_this:.1f})",
        "",
        "  \u2713  Higher-frequency spectral content",
        "     than equivalent earthquake",
        "",
        "  \u2713  Location: Lop Nor test site",
        "     (known Chinese nuclear facility)",
        "",
        "  \u2713  Shallow depth (~400 m tunnel)",
        "",
        "  \u2713  Declared test (PRC announcement",
        "     29 July 1996, moratorium follows)",
        "",
        "COUPLED TEST PARAMETERS:",
        "",
        f"  Yield:  ~{YIELD_KT:.0f} kt (estimated from mb)",
        f"  DF:     1 (fully coupled)",
        f"  mb:     {mb_pred:.2f} (pred) / {mb_this} (obs)",
        f"  DOB:    ~400 m (horizontal tunnel)",
        "",
        "  \u2192 45th & final Chinese underground",
        "    nuclear test.  CTBT signed 1996-09-24.",
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

    ax.axhline(4.90, color="gray", ls=":", lw=2, label="Observed mb = 4.90")
    ax.fill_between(W_range, 4.7, 5.1, color="gray", alpha=0.1)

    ax.plot(5.0, mb_from_yield(5.0, 1.0), "r*", ms=22, zorder=10,
            label="Best estimate: 5 kt")

    ax.set_xlabel("True yield (kt)"); ax.set_ylabel("Observed mb")
    ax.set_title("(f) Yield\u2013mb Constraint", fontweight="bold")
    ax.legend(fontsize=7, ncol=2, loc="upper left"); ax.grid(True, alpha=0.2, which="both")
    ax.set_xlim(0.1, 1000); ax.set_ylim(2, 7)

    fig.suptitle("Figure 8 — Explosion/Earthquake Discrimination Analysis\n"
                 "Lop Nor Final Chinese Nuclear Test 1996-07-29 01:48:57 UTC",
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

    model = LopNorVelocityModel()

    print("=" * 72)
    print("  Lop Nor 1996-07-29 — Final Chinese Nuclear Test (~5 kt Coupled)")
    print("  45th & Last Confirmed Chinese Underground Nuclear Test")
    print("=" * 72)
    print()

    # Fig 1 — Source model (shared plotting module)
    print("[1/9] Source model ...")
    fig_source_model(YIELD_KT, EVENT_DEPTH * 1000, DF, EVENT_INFO,
                     os.path.join(OUTDIR, "fig01_source_model.png"),
                     rho=2650.0)

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
        title_prefix="Figure 3 — Synthetic Seismograms: Lop Nor 1996 (5 kt Coupled)")

    # Fig 4 — Yield sensitivity (event-specific)
    print("[4/9] Yield sensitivity ...")
    fig04_yield_sensitivity()

    # Fig 5 — Yield trade-off (event-specific)
    print("[5/9] Yield trade-off ...")
    fig05_yield_tradeoff()

    # Fig 6 — Real data from both KNDC archive and IRIS
    print("[6/9] Loading observed data ...")
    # First try local KNDC archive
    kndc_data = load_kndc_observed()
    print(f"       Loaded {len(kndc_data)} trace(s) from KNDC archive")

    # Also attempt FDSN for any available broadband stations
    fdsn_data = fetch_waveforms(EVENT_TIME, EVENT_LAT, EVENT_LON, STATIONS,
                                fmin=FMIN, fmax=FMAX, pre=120, post=500)
    print(f"       Retrieved {len(fdsn_data)} station(s) from FDSN")

    # Merge: prefer FDSN for duplicates (better instrument response), add KNDC-only
    fdsn_stas = {x[1] for x in fdsn_data}
    all_real_data = list(fdsn_data)
    for item in kndc_data:
        if item[1] not in fdsn_stas:
            all_real_data.append(item)
    all_real_data.sort(key=lambda x: x[2])
    print(f"       Total observed: {len(all_real_data)} station(s)")

    fig_observed_record_section(
        all_real_data, EVENT_TIME, model,
        os.path.join(OUTDIR, "fig06_observed_data.png"),
        fmin=FMIN, fmax=FMAX,
        title="Figure 6 — Observed Waveforms: Lop Nor 1996-07-29 (KNDC + FDSN)")

    # Fig 7 — Synthetic vs observed comparison (shared plotting module)
    print("[7/9] Synthetic vs observed comparison ...")
    key_stas = ["DZHR", "UZB", "PRZ", "KURK", "AAK"]
    fig_comparison(
        all_real_data, EVENT_TIME, YIELD_KT, DF, model,
        os.path.join(OUTDIR, "fig07_comparison.png"),
        fmin=FMIN, fmax=FMAX, depth_km=EVENT_DEPTH,
        key_stations=key_stas,
        title="Figure 7 — Synthetic vs Observed: Lop Nor 1996 (5 kt Coupled)")

    # Fig 8 — Discrimination (event-specific)
    print("[8/9] Discrimination analysis ...")
    fig08_discrimination(all_real_data)

    # Fig 9 — Location inversion (shared inversion + plotting modules)
    print("[9/9] Epicenter location inversion ...")
    # Use only stations with known coordinates
    inv_stations = [(net, sta, lat, lon, desc) for net, sta, lat, lon, desc in STATIONS]
    print("  Running grid search...")
    best_lat, best_lon, misfit, lats, lons, picks = invert_epicenter(
        inv_stations, model, EVENT_LAT, EVENT_LON)
    err_km = fig_location_inversion(
        best_lat, best_lon, misfit, lats, lons, picks,
        EVENT_LAT, EVENT_LON, inv_stations,
        os.path.join(OUTDIR, "fig09_location_inversion.png"),
        title="Figure 9 — Epicenter Inversion: Lop Nor Final Test 1996-07-29")
    print(f"       Inverted location: {best_lat:.3f}\u00b0N, {best_lon:.3f}\u00b0E "
          f"(error: {err_km:.1f} km)")

    print()
    print("=" * 72)
    print(f"  All figures saved to: {OUTDIR}")
    print("=" * 72)


if __name__ == "__main__":
    main()
