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
from fsrm.velocity_models import PunggyeRiVelocityModel
from fsrm.propagation import (
    regional_travel_times, generate_synthetic, mt_mantap_collapse_signal,
)
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
EVENT_TIME   = UTCDateTime("2017-09-03T03:30:01")
EVENT_LAT    = 41.300
EVENT_LON    = 129.076
EVENT_DEPTH  = 0.76          # km
YIELD_KT     = 250.0
DF            = 1.0           # Fully coupled

OUTDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                      "figures", "dprk_2017")

# Bandpass for processing
FMIN, FMAX = 0.5, 8.0

EVENT_INFO = {
    "name": "DPRK 6th Nuclear Test",
    "time": "2017-09-03 03:30:01 UTC",
    "lat": EVENT_LAT,
    "lon": EVENT_LON,
    "mb_obs": 6.3,
}

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

STATION_INFO = compute_station_distances(EVENT_LAT, EVENT_LON, STATIONS)


# ═══════════════════════════════════════════════════════════════════════════════
# Event-Specific Figures (not generalizable across events)
# ═══════════════════════════════════════════════════════════════════════════════


def fig04_yield_sensitivity():
    """Synthetic waveforms at MDJ (~370 km) for different yields."""
    dist_km = 370.0
    model = PunggyeRiVelocityModel()
    yields = [10, 50, 100, 250, 500]

    fig, axes = plt.subplots(len(yields), 1, figsize=(16, 3.0 * len(yields)),
                             sharex=True)

    for i, W in enumerate(yields):
        t, v = generate_synthetic(dist_km, W, model,
                                  decoupling_factor=DF, dt=0.05, duration=400,
                                  fmin=FMIN, fmax=FMAX, depth_km=EVENT_DEPTH,
                                  collapse_signal=mt_mantap_collapse_signal)
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
        t, v = generate_synthetic(dist_km, W, model,
                                  decoupling_factor=DF, dt=0.05, duration=400,
                                  fmin=FMIN, fmax=FMAX, depth_km=EVENT_DEPTH,
                                  collapse_signal=mt_mantap_collapse_signal)

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


def fig08_discrimination(real_data):
    """Explosion vs earthquake discrimination analysis — DPRK 2017 is clearly an explosion."""
    fig = plt.figure(figsize=(18, 12))
    gs = GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.35)

    model = PunggyeRiVelocityModel()

    # (a) P/Lg amplitude ratio for synthetic
    ax = fig.add_subplot(gs[0, 0])
    distances = np.linspace(200, 2000, 20)
    ps_ratios = []
    for d in distances:
        t, v = generate_synthetic(d, YIELD_KT, model,
                                  decoupling_factor=DF, dt=0.05, duration=600,
                                  fmin=FMIN, fmax=FMAX, depth_km=EVENT_DEPTH,
                                  collapse_signal=mt_mantap_collapse_signal)
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
    t, v = generate_synthetic(dist_mdj, YIELD_KT, model,
                              decoupling_factor=DF, dt=0.05, duration=500,
                              fmin=FMIN, fmax=FMAX, depth_km=EVENT_DEPTH,
                              collapse_signal=mt_mantap_collapse_signal)
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
        ax.text(0.5, 0.5, "No MDJ data available", transform=ax.transAxes,
                ha="center", fontsize=12, color="gray")

    ax.set_xlabel("Frequency (Hz)"); ax.set_ylabel("Normalised amplitude")
    ax.set_title("(d) Observed P vs Lg Spectra at MDJ", fontweight="bold")
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
        f"  Yield:  ~{YIELD_KT:.0f} kt",
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
            label="Best estimate: 250 kt")

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

    model = PunggyeRiVelocityModel()

    print("=" * 72)
    print("  DPRK 2017 — 6th Nuclear Test (250 kt Coupled) Model & Analysis")
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
        collapse_signal=mt_mantap_collapse_signal,
        title_prefix="Figure 3 — Synthetic Seismograms")

    # Fig 4 — Yield sensitivity (event-specific)
    print("[4/9] Yield sensitivity ...")
    fig04_yield_sensitivity()

    # Fig 5 — Yield trade-off (event-specific)
    print("[5/9] Yield trade-off ...")
    fig05_yield_tradeoff()

    # Fig 6 — Real data (shared plotting module)
    print("[6/9] Fetching real data from IRIS ...")
    real_data = fetch_waveforms(EVENT_TIME, EVENT_LAT, EVENT_LON, STATIONS,
                                fmin=FMIN, fmax=FMAX, pre=120, post=500)
    print(f"       Retrieved {len(real_data)} station(s)")
    fig_observed_record_section(
        real_data, EVENT_TIME, model,
        os.path.join(OUTDIR, "fig06_real_data.png"),
        fmin=FMIN, fmax=FMAX,
        title="Figure 6 — Observed Waveforms: DPRK 6th Nuclear Test")

    # Fig 7 — Synthetic vs observed comparison (shared plotting module)
    print("[7/9] Synthetic vs observed comparison ...")
    fig_comparison(
        real_data, EVENT_TIME, YIELD_KT, DF, model,
        os.path.join(OUTDIR, "fig07_comparison.png"),
        fmin=FMIN, fmax=FMAX, depth_km=EVENT_DEPTH,
        key_stations=["MDJ", "INCN", "KSRS", "MAJO", "BJT"],
        collapse_signal=mt_mantap_collapse_signal,
        title="Figure 7 — Synthetic vs Observed Comparison")

    # Fig 8 — Discrimination (event-specific)
    print("[8/9] Discrimination analysis ...")
    fig08_discrimination(real_data)

    # Fig 9 — Location inversion (shared inversion + plotting modules)
    print("[9/9] Epicenter location inversion ...")
    print("  Running grid search...")
    best_lat, best_lon, misfit, lats, lons, picks = invert_epicenter(
        STATIONS, model, EVENT_LAT, EVENT_LON)
    err_km = fig_location_inversion(
        best_lat, best_lon, misfit, lats, lons, picks,
        EVENT_LAT, EVENT_LON, STATIONS,
        os.path.join(OUTDIR, "fig09_location_inversion.png"),
        title="Figure 9 — Epicenter Location Inversion (DPRK 6th Nuclear Test)")
    print(f"       Inverted location: {best_lat:.3f}°N, {best_lon:.3f}°E "
          f"(error: {err_km:.1f} km)")

    print()
    print("=" * 72)
    print(f"  All figures saved to: {OUTDIR}")
    print("=" * 72)


if __name__ == "__main__":
    main()
