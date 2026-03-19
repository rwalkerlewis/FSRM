#!/usr/bin/env python3
"""
Mueller-Murphy Seismic Source Model: Self-Contained Demo

A standalone, zero-dependency (beyond numpy/matplotlib/scipy) implementation
of the Mueller-Murphy underground nuclear explosion source model. Produces six
publication-quality figures demonstrating:

  1. Cavity scaling and damage zone geometry
  2. Body-wave magnitude and mb/Ms discrimination
  3. Reduced Displacement Potential (RDP), amplitude spectrum, far-field
     P velocity, and corner frequency scaling
  4. Decoupling effects on mb and RDP
  5. DPRK nuclear test series progression with model predictions
  6. Config catalog global map of all historical test configurations

Physical references:
  Mueller & Murphy (1971)  -- Seismic characteristics of underground nuclear detonations
  Patton (1988)            -- Corner-frequency scaling
  Glasstone & Dolan (1977) -- Effects of Nuclear Weapons
  Ringdal et al. (1992)    -- IMS detection thresholds

Usage:
    python scripts/mueller_murphy_demo.py
"""

import os
import sys
import glob
import configparser

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# ============================================================================
# Physical constants
# ============================================================================
JOULES_PER_KT = 4.184e12


# ============================================================================
# Mueller-Murphy Source Model (self-contained)
# ============================================================================

class MuellerMurphySource:
    """
    Self-contained Mueller-Murphy explosion source model.

    Parameters:
        yield_kt:  yield in kilotons
        rho:       host rock density (kg/m3)
        vp:        host rock compressional wave velocity (m/s)
        depth_m:   depth of burial (m)
    """

    def __init__(self, yield_kt, rho=2650.0, vp=5600.0, depth_m=300.0):
        self.yield_kt = yield_kt
        self.rho = rho
        self.vp = vp
        self.depth_m = depth_m

    @property
    def cavity_radius(self):
        """NTS empirical cavity radius (m): R_c = 55 * W^0.295 * (rho/2650)^(-1/3.4)."""
        return 55.0 * self.yield_kt**0.295 * (self.rho / 2650.0)**(-1.0 / 3.4)

    @property
    def zone_radii(self):
        """Damage zone radii (m) as dict."""
        rc = self.cavity_radius
        return {
            "cavity": rc,
            "crushed": 2.5 * rc,
            "fractured": 5.0 * rc,
            "damaged": 10.0 * rc,
        }

    @property
    def corner_frequency(self):
        """Patton (1988): fc ~ 3.0 * W^(-1/3) Hz (hard rock calibration)."""
        return 3.0 * self.yield_kt**(-1.0 / 3.0)

    @property
    def scalar_moment(self):
        """Empirical fully-coupled scalar moment (N*m): log10(M0) ~ 17 + log10(W)."""
        return 10.0**(17.0 + np.log10(max(self.yield_kt, 1e-6)))

    @staticmethod
    def mb_from_yield(yield_kt, decoupling_factor=1.0):
        """Body-wave magnitude: mb ~ 4.45 + 0.75 * log10(W_eff) (IASPEI NTS hard rock)."""
        w_eff = yield_kt / max(decoupling_factor, 1e-12)
        return 4.45 + 0.75 * np.log10(max(w_eff, 1e-6))

    @staticmethod
    def Mw_from_M0(M0):
        """Moment magnitude from scalar moment (N*m)."""
        return (np.log10(M0) - 9.1) / 1.5

    def spectrum(self, f, overshoot=1.1):
        """
        Mueller-Murphy displacement spectrum:
            Psi(f) = Psi_inf * B / [1 + (f/fc)^2]
        """
        fc = self.corner_frequency
        psi_inf = self.scalar_moment
        return overshoot * psi_inf / (1.0 + (f / fc)**2)

    def rdp_time(self, t, overshoot=1.1):
        """
        Reduced Displacement Potential psi(t) and moment-rate dpsi/dt.

        Brune-style: psi(t) = Psi_inf * B * [1 - (1 + t/tau) * exp(-t/tau)]

        Returns (psi, dpsi) as numpy arrays.
        """
        fc = self.corner_frequency
        tau = 1.0 / (2.0 * np.pi * fc)
        psi_inf = self.scalar_moment
        t = np.asarray(t, dtype=float)
        psi = np.zeros_like(t)
        dpsi = np.zeros_like(t)
        m = t > 0
        x = t[m] / tau
        psi[m] = overshoot * psi_inf * (1.0 - (1.0 + x) * np.exp(-x))
        dpsi[m] = overshoot * psi_inf * x * np.exp(-x) / tau
        return psi, dpsi

    def far_field_velocity(self, t, dist_m):
        """
        Far-field P-wave velocity pulse (m/s) at distance dist_m.

        v(t) = (1 / 4 pi rho vp^3 r) * d^2(psi)/dt^2
        """
        fc = self.corner_frequency
        tau = 1.0 / (2.0 * np.pi * fc)
        psi_inf = self.scalar_moment * 1.1
        t = np.asarray(t, dtype=float)
        d2psi = np.zeros_like(t)
        m = t > 0
        x = t[m] / tau
        d2psi[m] = psi_inf * (1.0 - x) * np.exp(-x) / tau**2

        factor = 1.0 / (4.0 * np.pi * self.rho * self.vp**3 * dist_m)
        return factor * d2psi


# ============================================================================
# Figure 1: Cavity Scaling and Damage Zones
# ============================================================================

def fig01_cavity_scaling(outdir):
    """Cavity radius and damage zone radii vs. yield for different host rocks."""
    yields = np.logspace(-1, 4, 200)
    rocks = [
        ("Granite", 2650.0, "#1f77b4"),
        ("Salt", 2200.0, "#2ca02c"),
        ("Tuff", 2100.0, "#ff7f0e"),
        ("Shale", 2400.0, "#d62728"),
    ]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

    # Left: cavity radius vs yield for different rocks
    for name, rho, col in rocks:
        rc = np.array([MuellerMurphySource(W, rho=rho).cavity_radius for W in yields])
        ax1.loglog(yields, rc, "-", color=col, lw=2, label=f"{name} ({rho} kg/m3)")

    ax1.set_xlabel("Yield (kt)", fontsize=12, fontweight="bold")
    ax1.set_ylabel("Cavity Radius (m)", fontsize=12, fontweight="bold")
    ax1.set_title("(a) Cavity Radius vs. Yield", fontsize=13, fontweight="bold")
    ax1.legend(fontsize=10)
    ax1.grid(True, which="both", alpha=0.2)
    ax1.set_xlim(0.1, 1e4)

    # Mark known test points
    known_tests = [
        (250.0, 2700.0, "DPRK 2017\n(250 kt)", "o", "#1f77b4"),
        (1200.0, 2400.0, "Cannikin\n(~5 Mt)", "s", "#ff7f0e"),
        (104.0, 2200.0, "Sedan\n(104 kt)", "^", "#2ca02c"),
    ]
    for W, rho, label, marker, col in known_tests:
        rc = MuellerMurphySource(W, rho=rho).cavity_radius
        ax1.plot(W, rc, marker, color=col, ms=10, markeredgecolor="k",
                 markeredgewidth=1, zorder=5)
        ax1.annotate(label, (W, rc), textcoords="offset points", xytext=(10, 5),
                     fontsize=8, fontweight="bold")

    # Right: damage zones for 250 kt granite
    src = MuellerMurphySource(250.0, rho=2700.0)
    zones = src.zone_radii
    zone_names = ["cavity", "crushed", "fractured", "damaged"]
    zone_colors = ["#d62728", "#ff7f0e", "#2ca02c", "#1f77b4"]
    zone_labels = [
        f"Cavity ({zones['cavity']:.0f} m)",
        f"Crushed zone ({zones['crushed']:.0f} m)",
        f"Fractured zone ({zones['fractured']:.0f} m)",
        f"Damaged zone ({zones['damaged']:.0f} m)",
    ]

    for i in range(len(zone_names) - 1, -1, -1):
        r = zones[zone_names[i]]
        circle = plt.Circle((0, 0), r, color=zone_colors[i], alpha=0.3, label=zone_labels[i])
        ax2.add_patch(circle)
        ax2.plot(r, 0, "|", color=zone_colors[i], ms=15, mew=2)
        ax2.text(r + 20, 50 + i * 100, f"{r:.0f} m", fontsize=9,
                 color=zone_colors[i], fontweight="bold")

    max_r = zones["damaged"] * 1.15
    ax2.set_xlim(-max_r, max_r)
    ax2.set_ylim(-max_r, max_r)
    ax2.set_aspect("equal")
    ax2.set_xlabel("Distance (m)", fontsize=12, fontweight="bold")
    ax2.set_ylabel("Distance (m)", fontsize=12, fontweight="bold")
    ax2.set_title("(b) Damage Zones: 250 kt in Granite", fontsize=13, fontweight="bold")
    ax2.legend(loc="upper left", fontsize=9)
    ax2.grid(True, alpha=0.15)

    fig.suptitle("Figure 1: Cavity Scaling and Damage Zone Geometry\n"
                 "Mueller-Murphy/Glasstone-Dolan empirical models",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    p = os.path.join(outdir, "fig01_cavity_scaling.png")
    fig.savefig(p, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {p}")
    return zones


# ============================================================================
# Figure 2: Magnitude Discrimination (mb/Ms)
# ============================================================================

def fig02_magnitude_discrimination(outdir):
    """mb/Ms scaling, explosion-earthquake discrimination diagram."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

    # Left: mb vs yield
    yields = np.logspace(-1, 4, 200)
    mb_coupled = np.array([MuellerMurphySource.mb_from_yield(W) for W in yields])

    ax1.semilogx(yields, mb_coupled, "b-", lw=2.5, label="Coupled (DF=1)")

    for df, ls in [(5, "--"), (10, "-."), (70, ":")]:
        mb_dec = np.array([MuellerMurphySource.mb_from_yield(W, df) for W in yields])
        ax1.semilogx(yields, mb_dec, ls=ls, lw=1.5, alpha=0.7, label=f"Decoupled (DF={df})")

    ax1.axhline(3.5, color="gray", ls=":", lw=2, label="IMS threshold (~3.5)")
    ax1.fill_between(yields, 3.0, 4.0, color="gray", alpha=0.05)

    # Mark DPRK tests
    dprk_tests = [
        (0.7, 3.9, "2006"),
        (4.0, 4.5, "2009"),
        (10.0, 4.9, "2013"),
        (15.0, 5.1, "2016a"),
        (20.0, 5.3, "2016b"),
        (250.0, 6.3, "2017"),
    ]
    for W, mb, label in dprk_tests:
        ax1.plot(W, mb, "r*", ms=12, zorder=5)
        ax1.annotate(label, (W, mb), textcoords="offset points",
                     xytext=(8, 3), fontsize=8, fontweight="bold", color="r")

    ax1.set_xlabel("Yield (kt)", fontsize=12, fontweight="bold")
    ax1.set_ylabel("mb", fontsize=12, fontweight="bold")
    ax1.set_title("(a) Body-Wave Magnitude vs. Yield", fontsize=13, fontweight="bold")
    ax1.legend(fontsize=9, loc="upper left")
    ax1.grid(True, which="both", alpha=0.2)
    ax1.set_xlim(0.1, 1e4)
    ax1.set_ylim(1, 8)

    # Right: mb-Ms discrimination
    mb_exp = np.linspace(3, 7, 50)
    ms_exp = mb_exp - 1.2  # Explosion trend (Ms deficit)

    rng = np.random.RandomState(42)
    mb_eq = np.linspace(3, 7.5, 80)
    ms_eq = 1.0 * mb_eq - 0.5 + 0.3 * rng.randn(80)

    ax2.scatter(mb_eq, ms_eq, c="b", s=20, alpha=0.3, label="Earthquakes")
    ax2.plot(mb_exp, ms_exp, "r-", lw=2, alpha=0.6, label="Explosion trend")

    # Known explosion points
    ax2.plot(6.3, 4.5, "r*", ms=18, zorder=10,
             label="DPRK 2017 (mb=6.3, Ms~4.5)")
    ax2.plot(6.97, 5.2, "rs", ms=10, zorder=10, markeredgecolor="k",
             label="Cannikin (mb=6.97)")

    # Discrimination line
    mb_line = np.linspace(2, 8, 100)
    ms_line = mb_line - 1.0
    ax2.plot(mb_line, ms_line, "k--", lw=1.5, label="Discrimination line")
    ax2.fill_between(mb_line, ms_line - 3, ms_line, color="red", alpha=0.03)
    ax2.fill_between(mb_line, ms_line, ms_line + 3, color="blue", alpha=0.03)
    ax2.text(4.0, 2.0, "EXPLOSION", fontsize=11, color="red", alpha=0.5,
             fontweight="bold", fontstyle="italic")
    ax2.text(4.0, 5.0, "EARTHQUAKE", fontsize=11, color="blue", alpha=0.5,
             fontweight="bold", fontstyle="italic")

    ax2.set_xlabel("mb (body-wave magnitude)", fontsize=12, fontweight="bold")
    ax2.set_ylabel("Ms (surface-wave magnitude)", fontsize=12, fontweight="bold")
    ax2.set_title("(b) mb-Ms Discrimination Diagram", fontsize=13, fontweight="bold")
    ax2.legend(fontsize=8, loc="upper left")
    ax2.grid(True, alpha=0.2)
    ax2.set_xlim(2, 8)
    ax2.set_ylim(1, 8)

    fig.suptitle("Figure 2: Magnitude Scaling and Explosion/Earthquake Discrimination\n"
                 "mb = 4.0 + 0.75 log10(W)  |  Explosions show systematic Ms deficit",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    p = os.path.join(outdir, "fig02_magnitude_discrimination.png")
    fig.savefig(p, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {p}")

    # Return DPRK 2017 predicted mb for validation
    return MuellerMurphySource.mb_from_yield(250.0)


# ============================================================================
# Figure 3: RDP Waveforms, Spectrum, Far-field P, Corner Frequency
# ============================================================================

def fig03_rdp_waveforms(outdir):
    """Time-domain RDP, amplitude spectrum, far-field P velocity, corner frequency scaling."""
    fig = plt.figure(figsize=(18, 14))
    gs = GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.3)

    # (a) Time-domain RDP for different yields
    ax = fig.add_subplot(gs[0, 0])
    t = np.linspace(-0.5, 10, 4000)
    test_yields = [1, 10, 100, 250, 1000]
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(test_yields)))
    for W, col in zip(test_yields, colors):
        src = MuellerMurphySource(W, rho=2700.0)
        psi, _ = src.rdp_time(t)
        psi_norm = psi / src.scalar_moment if src.scalar_moment > 0 else psi
        ax.plot(t, psi_norm, "-", color=col, lw=1.5,
                label=f"{W} kt (fc={src.corner_frequency:.2f} Hz)")
    ax.set_xlabel("Time (s)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Normalised RDP (psi / M0)", fontsize=11, fontweight="bold")
    ax.set_title("(a) Reduced Displacement Potential", fontsize=12, fontweight="bold")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.2)
    ax.set_xlim(-0.5, 10)

    # (b) Amplitude spectrum
    ax = fig.add_subplot(gs[0, 1])
    f = np.logspace(-2, 2, 500)
    for W, col in zip(test_yields, colors):
        src = MuellerMurphySource(W, rho=2700.0)
        spec = src.spectrum(f)
        spec_norm = spec / src.scalar_moment if src.scalar_moment > 0 else spec
        ax.loglog(f, spec_norm, "-", color=col, lw=1.5,
                  label=f"{W} kt (fc={src.corner_frequency:.2f} Hz)")
        ax.axvline(src.corner_frequency, color=col, ls=":", lw=0.8, alpha=0.5)
    ax.set_xlabel("Frequency (Hz)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Normalised displacement spectrum", fontsize=11, fontweight="bold")
    ax.set_title("(b) Mueller-Murphy Source Spectrum", fontsize=12, fontweight="bold")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.2)

    # (c) Far-field P velocity at 1000 km for 250 kt
    ax = fig.add_subplot(gs[1, 0])
    src = MuellerMurphySource(250.0, rho=2700.0, vp=5600.0)
    t_ff = np.linspace(-0.5, 15, 5000)
    v_ff = src.far_field_velocity(t_ff, 1000e3)
    ax.plot(t_ff, v_ff * 1e9, "k-", lw=1)  # convert to nm/s
    ax.fill_between(t_ff, v_ff * 1e9, 0, where=(v_ff > 0), color="red", alpha=0.15)
    ax.fill_between(t_ff, v_ff * 1e9, 0, where=(v_ff < 0), color="blue", alpha=0.15)
    ax.set_xlabel("Time (s)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Velocity (nm/s)", fontsize=11, fontweight="bold")
    ax.set_title("(c) Far-field P Velocity at 1000 km (250 kt, granite)",
                 fontsize=12, fontweight="bold")
    ax.grid(True, alpha=0.2)
    ax.axhline(0, color="gray", lw=0.5)

    # (d) Corner frequency scaling
    ax = fig.add_subplot(gs[1, 1])
    yields_fc = np.logspace(-1, 4, 200)
    fc = np.array([MuellerMurphySource(W).corner_frequency for W in yields_fc])
    ax.loglog(yields_fc, fc, "b-", lw=2.5, label="Patton: fc = 2.5 W^(-1/3)")

    # Mark some specific values
    for W, marker, label in [(1, "o", "1 kt"), (10, "s", "10 kt"),
                              (100, "^", "100 kt"), (250, "*", "250 kt"),
                              (1000, "D", "1 Mt")]:
        fc_val = MuellerMurphySource(W).corner_frequency
        ax.plot(W, fc_val, marker, color="red", ms=10, markeredgecolor="k",
                markeredgewidth=0.8, zorder=5)
        ax.annotate(f"{label}\n{fc_val:.3f} Hz", (W, fc_val),
                    textcoords="offset points", xytext=(10, 5), fontsize=8)

    ax.set_xlabel("Yield (kt)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Corner frequency (Hz)", fontsize=11, fontweight="bold")
    ax.set_title("(d) Corner Frequency Scaling", fontsize=12, fontweight="bold")
    ax.legend(fontsize=10)
    ax.grid(True, which="both", alpha=0.2)
    ax.set_xlim(0.1, 1e4)

    fig.suptitle("Figure 3: Mueller-Murphy Source Model\n"
                 "RDP, Spectrum, Far-field P-wave, Corner Frequency",
                 fontsize=15, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    p = os.path.join(outdir, "fig03_rdp_waveforms.png")
    fig.savefig(p, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {p}")

    fc_250 = MuellerMurphySource(250.0).corner_frequency
    return fc_250


# ============================================================================
# Figure 4: Decoupling Effect
# ============================================================================

def fig04_decoupling(outdir):
    """Decoupling effect on mb, RDP comparison tamped vs. cavity-decoupled."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

    # (a) mb vs yield for different decoupling factors
    yields = np.logspace(-1, 4, 200)
    ax1.semilogx(yields, [MuellerMurphySource.mb_from_yield(W) for W in yields],
                 "b-", lw=2.5, label="Coupled (DF=1)")

    df_list = [5, 10, 70, 290]
    df_labels = ["DF=5", "DF=10", "DF=70 (hard rock cavity)", "DF=290 (salt cavity, theoretical max)"]
    styles = ["--", "-.", ":", (0, (3, 5, 1, 5))]
    for df, ls, dl in zip(df_list, styles, df_labels):
        mb_vals = [MuellerMurphySource.mb_from_yield(W, df) for W in yields]
        ax1.semilogx(yields, mb_vals, ls=ls, lw=1.5, alpha=0.7, label=dl)

    ax1.axhline(3.5, color="gray", ls=":", lw=2, label="IMS threshold (~3.5)")

    # Mark the 10 kt decoupled-in-salt case
    DF_SALT = 290  # Theoretical maximum for large cavity in salt
    mb_10kt_coupled = MuellerMurphySource.mb_from_yield(10.0)
    mb_10kt_dec = MuellerMurphySource.mb_from_yield(10.0, DF_SALT)
    ax1.plot(10.0, mb_10kt_coupled, "ro", ms=10, zorder=5)
    ax1.plot(10.0, mb_10kt_dec, "rv", ms=10, zorder=5)
    ax1.annotate(f"10 kt coupled: mb={mb_10kt_coupled:.2f}",
                 (10.0, mb_10kt_coupled), textcoords="offset points",
                 xytext=(15, 5), fontsize=9, color="r", fontweight="bold")
    ax1.annotate(f"10 kt dec (DF={DF_SALT}): mb={mb_10kt_dec:.2f}",
                 (10.0, mb_10kt_dec), textcoords="offset points",
                 xytext=(15, -15), fontsize=9, color="r", fontweight="bold")

    ax1.set_xlabel("Yield (kt)", fontsize=12, fontweight="bold")
    ax1.set_ylabel("mb", fontsize=12, fontweight="bold")
    ax1.set_title("(a) Decoupling Effect on mb", fontsize=13, fontweight="bold")
    ax1.legend(fontsize=9)
    ax1.grid(True, which="both", alpha=0.2)
    ax1.set_xlim(0.1, 1e4)
    ax1.set_ylim(0, 8)

    # (b) RDP comparison: 10 kt coupled vs. 10 kt decoupled (DF=290, salt)
    t = np.linspace(-0.2, 8, 3000)
    src_coupled = MuellerMurphySource(10.0, rho=2200.0)
    psi_c, _ = src_coupled.rdp_time(t)

    # Decoupled: effective yield = W/DF in terms of amplitude
    src_decoupled = MuellerMurphySource(10.0 / DF_SALT, rho=2200.0)
    psi_d, _ = src_decoupled.rdp_time(t)

    ax2.plot(t, psi_c / np.max(np.abs(psi_c)), "b-", lw=2, label="10 kt coupled")
    ax2.plot(t, psi_d / np.max(np.abs(psi_c)), "r--", lw=2,
             label=f"10 kt decoupled (DF={DF_SALT})")
    ax2.set_xlabel("Time (s)", fontsize=12, fontweight="bold")
    ax2.set_ylabel("Normalised RDP", fontsize=12, fontweight="bold")
    ax2.set_title("(b) RDP: Coupled vs. Cavity-Decoupled (10 kt, salt)",
                  fontsize=13, fontweight="bold")
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.2)

    mb_drop = mb_10kt_coupled - mb_10kt_dec
    fig.suptitle(f"Figure 4: Decoupling Effects on Seismic Observables\n"
                 f"Cavity decoupling reduces mb by up to {mb_drop:.2f} magnitude units (DF={DF_SALT})",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    p = os.path.join(outdir, "fig04_decoupling.png")
    fig.savefig(p, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {p}")

    # Compute equivalent coupled yield
    w_eff = 10.0 / DF_SALT
    return mb_10kt_dec, w_eff


# ============================================================================
# Figure 5: DPRK Test Series Progression
# ============================================================================

def fig05_dprk_series(outdir):
    """DPRK nuclear test series with Mueller-Murphy model predictions vs. observed."""
    # DPRK tests: (year, name, observed_mb, estimated_yield_kt)
    dprk_tests = [
        (2006, "DPRK-1", 4.1, 0.7),
        (2009, "DPRK-2", 4.52, 4.0),
        (2013, "DPRK-3", 4.9, 10.0),
        (2016.1, "DPRK-4", 4.85, 10.0),
        (2016.7, "DPRK-5", 5.04, 20.0),
        (2017, "DPRK-6", 6.3, 250.0),
    ]

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 7))

    years = [t[0] for t in dprk_tests]
    names = [t[1] for t in dprk_tests]
    obs_mb = np.array([t[2] for t in dprk_tests])
    est_yield = np.array([t[3] for t in dprk_tests])
    pred_mb = np.array([MuellerMurphySource.mb_from_yield(W) for W in est_yield])

    # (a) Observed vs predicted mb
    ax1.plot(years, obs_mb, "ro-", ms=10, lw=2, label="Observed mb", zorder=5)
    ax1.plot(years, pred_mb, "b^--", ms=10, lw=1.5, label="Predicted mb (M-M)")

    for i, name in enumerate(names):
        ax1.annotate(name, (years[i], obs_mb[i]), textcoords="offset points",
                     xytext=(0, 12), fontsize=8, fontweight="bold", ha="center")
        ax1.annotate(f"{est_yield[i]} kt", (years[i], obs_mb[i]),
                     textcoords="offset points", xytext=(0, -18), fontsize=7,
                     color="gray", ha="center")

    residuals = obs_mb - pred_mb
    ax1.set_xlabel("Year", fontsize=12, fontweight="bold")
    ax1.set_ylabel("mb", fontsize=12, fontweight="bold")
    ax1.set_title("(a) DPRK Test Series: mb Progression", fontsize=13, fontweight="bold")
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.2)

    # (b) Yield growth
    ax2.semilogy(years, est_yield, "rs-", ms=10, lw=2, label="Estimated yield")
    for i, name in enumerate(names):
        ax2.annotate(f"{est_yield[i]} kt", (years[i], est_yield[i]),
                     textcoords="offset points", xytext=(8, 5), fontsize=9,
                     fontweight="bold")

    ax2.set_xlabel("Year", fontsize=12, fontweight="bold")
    ax2.set_ylabel("Estimated Yield (kt)", fontsize=12, fontweight="bold")
    ax2.set_title("(b) Yield Growth", fontsize=13, fontweight="bold")
    ax2.legend(fontsize=10)
    ax2.grid(True, which="both", alpha=0.2)

    # (c) Spectra comparison for each test
    f = np.logspace(-2, 1, 300)
    colors = plt.cm.hot(np.linspace(0.2, 0.9, len(dprk_tests)))
    for i, (yr, name, mb_obs, W) in enumerate(dprk_tests):
        src = MuellerMurphySource(W, rho=2700.0)
        spec = src.spectrum(f) / MuellerMurphySource(250.0, rho=2700.0).scalar_moment
        ax3.loglog(f, spec, "-", color=colors[i], lw=1.5,
                   label=f"{name} ({W} kt)")
        ax3.axvline(src.corner_frequency, color=colors[i], ls=":", lw=0.6, alpha=0.5)

    ax3.set_xlabel("Frequency (Hz)", fontsize=12, fontweight="bold")
    ax3.set_ylabel("Relative spectrum (normalised to DPRK-6)", fontsize=12,
                   fontweight="bold")
    ax3.set_title("(c) Source Spectra", fontsize=13, fontweight="bold")
    ax3.legend(fontsize=8)
    ax3.grid(True, which="both", alpha=0.2)

    fig.suptitle("Figure 5: DPRK Nuclear Test Series (2006-2017)\n"
                 "Mueller-Murphy model predictions vs. observed mb",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    p = os.path.join(outdir, "fig05_dprk_series.png")
    fig.savefig(p, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {p}")

    return pred_mb[-1]  # predicted mb for DPRK 2017


# ============================================================================
# Figure 6: Global Config Catalog Map
# ============================================================================

def fig06_config_catalog(outdir, config_dir):
    """Plot all historical test config locations on a world map (matplotlib only)."""
    # Parse config files for locations
    configs = []
    config_files = sorted(glob.glob(os.path.join(config_dir, "*.config")))
    for cf in config_files:
        parser = configparser.ConfigParser(inline_comment_prefixes=("#",))
        try:
            parser.read(cf)
        except Exception:
            continue

        name = os.path.basename(cf).replace(".config", "")
        lat = lon = None

        # Try multiple places where location might be stored
        for section in parser.sections():
            for key in parser.options(section):
                val = parser.get(section, key).strip()
                if "local_origin_y" in key or key == "event_latitude":
                    try:
                        lat = float(val.split()[0])
                    except (ValueError, IndexError):
                        pass
                if "local_origin_x" in key or key == "event_longitude":
                    try:
                        lon = float(val.split()[0])
                    except (ValueError, IndexError):
                        pass

        # Also check [NUCLEAR_DEVICE] and [EXPLOSION_SOURCE]
        for section in ["NUCLEAR_DEVICE", "EXPLOSION_SOURCE", "GRID", "NUCLEAR_EMP"]:
            if parser.has_section(section):
                for key in ["location_x", "longitude", "detonation_longitude",
                            "burst_longitude", "epicenter_lon"]:
                    if parser.has_option(section, key):
                        try:
                            lon = float(parser.get(section, key).split()[0])
                        except (ValueError, IndexError):
                            pass
                for key in ["location_y", "latitude", "detonation_latitude",
                            "burst_latitude", "epicenter_lat"]:
                    if parser.has_option(section, key):
                        try:
                            lat = float(parser.get(section, key).split()[0])
                        except (ValueError, IndexError):
                            pass

        if lat is not None and lon is not None:
            is_historical = "historical" in name.lower()
            configs.append((name, lat, lon, is_historical))

    fig, ax = plt.subplots(figsize=(18, 9), subplot_kw={"projection": None})

    # Simple world map outline (coastlines approximation using rectangles)
    # Draw background
    ax.set_facecolor("#e8f0f8")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_aspect("equal")

    # Approximate land mass outlines as filled polygons (simplified)
    # North America
    land_patches = [
        ([-170, -170, -50, -50], [15, 72, 72, 15], "#d0e0c0"),  # N America
        ([-85, -85, -30, -30], [-60, 15, 15, -60], "#d0e0c0"),  # S America
        ([-20, -20, 55, 55], [-35, 40, 40, -35], "#d0e0c0"),    # Africa
        ([-15, -15, 45, 45], [35, 72, 72, 35], "#d0e0c0"),      # Europe
        ([25, 25, 180, 180], [-10, 72, 72, -10], "#d0e0c0"),    # Asia
        ([110, 110, 155, 155], [-45, -10, -10, -45], "#d0e0c0"),# Australia
    ]
    for xs, ys, col in land_patches:
        ax.fill(xs, ys, color=col, alpha=0.5, zorder=1)

    # Grid lines
    for lat_g in range(-60, 90, 30):
        ax.axhline(lat_g, color="#ccc", lw=0.3, zorder=0)
    for lon_g in range(-150, 180, 30):
        ax.axvline(lon_g, color="#ccc", lw=0.3, zorder=0)

    # Plot configs
    hist_lats = [c[1] for c in configs if c[3]]
    hist_lons = [c[2] for c in configs if c[3]]
    other_lats = [c[1] for c in configs if not c[3]]
    other_lons = [c[2] for c in configs if not c[3]]

    if other_lats:
        ax.scatter(other_lons, other_lats, c="#1f77b4", s=60, marker="o",
                   edgecolors="k", linewidths=0.5, zorder=5, alpha=0.7,
                   label=f"Simulation configs ({len(other_lats)})")
    if hist_lats:
        ax.scatter(hist_lons, hist_lats, c="#d62728", s=100, marker="*",
                   edgecolors="k", linewidths=0.5, zorder=6, alpha=0.9,
                   label=f"Historical test configs ({len(hist_lats)})")

    # Label historical configs
    for name, lat, lon, is_hist in configs:
        if is_hist:
            short = name.replace("historical_", "").replace("_", " ")
            if len(short) > 40:
                short = short[:37] + "..."
            ax.annotate(short, (lon, lat), textcoords="offset points",
                        xytext=(8, 5), fontsize=6.5, fontweight="bold",
                        color="#d62728", zorder=7)

    ax.set_xlabel("Longitude", fontsize=12, fontweight="bold")
    ax.set_ylabel("Latitude", fontsize=12, fontweight="bold")
    ax.set_title("Figure 6: FSRM Configuration Catalog\n"
                 f"All geolocated configs ({len(configs)} total)",
                 fontsize=14, fontweight="bold")
    ax.legend(fontsize=10, loc="lower left")

    plt.tight_layout()
    p = os.path.join(outdir, "fig06_config_catalog.png")
    fig.savefig(p, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {p}")
    return len(configs)


# ============================================================================
# Main
# ============================================================================

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_dir = os.path.dirname(script_dir)
    outdir = os.path.join(repo_dir, "figures", "mueller_murphy_demo")
    config_dir = os.path.join(repo_dir, "config")
    os.makedirs(outdir, exist_ok=True)

    print("=" * 72)
    print("  Mueller-Murphy Seismic Source Model: Self-Contained Demo")
    print("=" * 72)
    print()

    # Section 1: Cavity Scaling
    print("[1/6] Cavity scaling and damage zones ...")
    zones = fig01_cavity_scaling(outdir)
    rc_250 = zones["cavity"]
    print(f"       DPRK 2017 (250 kt, granite): cavity radius = {rc_250:.1f} m")
    print()

    # Section 2: Magnitude Discrimination
    print("[2/6] Magnitude discrimination ...")
    mb_pred = fig02_magnitude_discrimination(outdir)
    print(f"       DPRK 2017 predicted mb: {mb_pred:.2f} (observed: 6.3)")
    print()

    # Section 3: RDP Waveforms
    print("[3/6] RDP waveforms, spectrum, far-field P, corner frequency ...")
    fc_250 = fig03_rdp_waveforms(outdir)
    print(f"       Corner frequency (250 kt): {fc_250:.3f} Hz")
    print()

    # Section 4: Decoupling
    print("[4/6] Decoupling effects ...")
    mb_dec, w_eff = fig04_decoupling(outdir)
    print(f"       Decoupled 10 kt in salt: mb {mb_dec:.2f} "
          f"(below IMS threshold of 3.5)")
    print(f"       Equivalent coupled yield for decoupled 10 kt: {w_eff:.3f} kt")
    print()

    # Section 5: DPRK Series
    print("[5/6] DPRK nuclear test series ...")
    mb_dprk_pred = fig05_dprk_series(outdir)
    print(f"       DPRK 2017 predicted mb: {mb_dprk_pred:.2f}")
    print()

    # Section 6: Config Catalog
    print("[6/6] Config catalog map ...")
    n_configs = fig06_config_catalog(outdir, config_dir)
    print(f"       Mapped {n_configs} geolocated configs")
    print()

    print("=" * 72)
    print(f"  All figures saved to: {outdir}")
    print("=" * 72)


if __name__ == "__main__":
    main()
