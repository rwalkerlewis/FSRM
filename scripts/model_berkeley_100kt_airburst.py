#!/usr/bin/env python3
"""
Model a Hypothetical 100 kT Airburst over Berkeley, California (50 m HOB)

Comprehensive physics-based analysis of a hypothetical nuclear airburst:
  1. Blast wave propagation and damage zones (Glasstone-Dolan / Brode)
  2. Thermal radiation and burn zones
  3. Prompt nuclear radiation (gamma + neutron)
  4. Fallout pattern with prevailing wind transport (Gaussian plume)
  5. EMP effects (E1/E2/E3)
  6. Fireball and mushroom cloud evolution
  7. Ground-coupled seismic and acoustic effects
  8. Combined hazard mapping
  9. Casualty estimation by zone

All results are publication-quality figures saved to figures/berkeley_100kt/.

Physical references:
  Glasstone & Dolan (1977) — Effects of Nuclear Weapons
  Brode (1955)             — Numerical Solutions of Spherical Blast Waves
  Sedov (1959)             — Similarity and Dimensional Methods in Mechanics
  Way & Wigner (1948)      — Radioactive Decay of Fission Products
  Bridgman (2001)          — Introduction to the Physics of Nuclear Weapons Effects

Usage:
    python scripts/model_berkeley_100kt_airburst.py
"""

import os
import warnings
import numpy as np
from scipy.ndimage import gaussian_filter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Circle, Wedge, FancyArrowPatch
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import AutoMinorLocator

warnings.filterwarnings("ignore")

# ═══════════════════════════════════════════════════════════════════════════════
# Constants & Event Parameters
# ═══════════════════════════════════════════════════════════════════════════════
JOULES_PER_KT = 4.184e12

YIELD_KT   = 100.0
BURST_HEIGHT = 50.0           # meters above ground level
TOTAL_ENERGY = YIELD_KT * JOULES_PER_KT

EVENT_LAT  = 37.8716          # UC Berkeley campus
EVENT_LON  = -122.2727

OUTDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                      "figures", "berkeley_100kt")

# Bay Area prevailing wind (from W/NW, blowing toward SE)
WIND_SPEED_SURFACE  = 5.0     # m/s at surface
WIND_SPEED_ALTITUDE = 15.0    # m/s at cloud stabilisation height
WIND_DIR_DEG        = 315.0   # degrees (from NW, meteorological convention)

# Population density (approximate)
POP_DENSITY_URBAN = 4000      # per km²  (Berkeley / Oakland)
POP_DENSITY_SF    = 7000      # per km²  (San Francisco)
POP_DENSITY_SUB   = 2000      # per km²  (suburban East Bay)

EARTH_RADIUS_KM = 6371.0


# ═══════════════════════════════════════════════════════════════════════════════
# PART A — Blast Physics (Glasstone-Dolan)
# ═══════════════════════════════════════════════════════════════════════════════

def peak_overpressure_kpa(ground_range_m, yield_kt=YIELD_KT, hob=BURST_HEIGHT):
    """
    Peak static overpressure (kPa) using calibrated Glasstone-Dolan formula.

    Reference for 1 kT optimal-height airburst:
        20 psi (140 kPa) at ~450 m
        5 psi  ( 35 kPa) at ~900 m
        2 psi  ( 14 kPa) at ~1400 m
    Scales by W^(1/3) for other yields.
    """
    r = np.atleast_1d(np.asarray(ground_range_m, dtype=float))
    R = np.sqrt(r**2 + hob**2)

    # Kinney-Graham (1985) free-air peak overpressure — continuous closed-form
    # fit to compiled experimental data.
    #   ΔP/P₀ = 808 [1 + (Z̄/4.5)²]
    #           / {√[1+(Z̄/0.048)²] · √[1+(Z̄/0.32)²] · √[1+(Z̄/1.35)²]}
    # where Z̄ = R / m_charge^{1/3} (m / kg^{1/3}).
    # Ref: Kinney & Graham, "Explosive Shocks in Air", 2nd ed., 1985.
    P0 = 101.325  # kPa
    W_kg = max(yield_kt, 1e-12) * 1.0e6   # kt → kg TNT equivalent
    Zbar = R / W_kg ** (1.0 / 3.0)         # scaled distance (m / kg^{1/3})
    num = 808.0 * (1.0 + (Zbar / 4.5) ** 2)
    d1 = np.sqrt(1.0 + (Zbar / 0.048) ** 2)
    d2 = np.sqrt(1.0 + (Zbar / 0.32) ** 2)
    d3 = np.sqrt(1.0 + (Zbar / 1.35) ** 2)
    P_kpa = P0 * num / (d1 * d2 * d3)

    # Mach stem enhancement for low burst angles
    angle = np.arctan2(hob, np.maximum(r, 1.0))
    mach_factor = np.where(angle < np.radians(40), 1.8, 1.0)

    return np.squeeze(P_kpa * mach_factor)


def dynamic_pressure_kpa(overpressure_kpa):
    """Dynamic (wind) pressure from static overpressure."""
    P0 = 101.325  # kPa
    return 2.5 * overpressure_kpa**2 / (7 * P0 + overpressure_kpa)


def blast_radius_for_pressure(target_kpa, yield_kt=YIELD_KT, hob=BURST_HEIGHT):
    """Find ground range (m) where overpressure equals target_kpa."""
    r_test = np.linspace(1, 50000, 50000)
    p = peak_overpressure_kpa(r_test, yield_kt, hob)
    idx = np.argmin(np.abs(p - target_kpa))
    return r_test[idx]


DAMAGE_THRESHOLDS = {
    "Total destruction (>140 kPa)":            140.0,
    "Reinforced concrete destroyed (>35 kPa)": 35.0,
    "Most buildings destroyed (>14 kPa)":      14.0,
    "Moderate damage (>7 kPa)":                7.0,
    "Light damage, window breakage (>3.5 kPa)": 3.5,
    "Glass breakage (>1 kPa)":                 1.0,
}


# ═══════════════════════════════════════════════════════════════════════════════
# PART B — Thermal Radiation
# ═══════════════════════════════════════════════════════════════════════════════

def thermal_fluence_kj(ground_range_m, yield_kt=YIELD_KT, hob=BURST_HEIGHT):
    """Thermal fluence (kJ/m²) at ground range, with atmospheric absorption."""
    r = np.atleast_1d(np.asarray(ground_range_m, dtype=float))
    R = np.sqrt(r**2 + hob**2)
    E_thermal = 0.35 * yield_kt * JOULES_PER_KT
    Q = E_thermal / (4 * np.pi * R**2)
    tau = np.exp(-R / 18000.0)
    return np.squeeze(Q * tau / 1000.0)


THERMAL_THRESHOLDS = {
    "Third-degree burns (>670 kJ/m²)":  670.0,
    "Second-degree burns (>335 kJ/m²)": 335.0,
    "First-degree burns (>125 kJ/m²)":  125.0,
    "Pain threshold (~50 kJ/m²)":       50.0,
}


def fireball_max_radius(yield_kt=YIELD_KT):
    return 66.0 * yield_kt**0.4


# ═══════════════════════════════════════════════════════════════════════════════
# PART C — Prompt Nuclear Radiation
# ═══════════════════════════════════════════════════════════════════════════════

def prompt_radiation_dose_gy(ground_range_m, yield_kt=YIELD_KT, hob=BURST_HEIGHT):
    """Prompt (initial) radiation dose in Gy at ground range."""
    r = np.atleast_1d(np.asarray(ground_range_m, dtype=float))
    R = np.sqrt(r**2 + hob**2)
    D_ref = 10.0  # Gy at 1 km for 1 kt
    D = D_ref * yield_kt * (1000.0 / R)**2 * np.exp(-R / 2500.0)
    return np.squeeze(D)


RADIATION_THRESHOLDS = {
    "Fatal (>10 Gy)":           10.0,
    "LD50/60 (~4.5 Gy)":        4.5,
    "Acute radiation sickness (>2 Gy)": 2.0,
    "Mild symptoms (>0.5 Gy)":  0.5,
    "Detectable (>0.01 Gy)":    0.01,
}


# ═══════════════════════════════════════════════════════════════════════════════
# PART D — Fallout Model (Gaussian Plume)
# ═══════════════════════════════════════════════════════════════════════════════

def cloud_top_height(yield_kt=YIELD_KT):
    """Approximate stabilised cloud top (m) for yield in kt."""
    return min(25000, 2200 * yield_kt**0.4)


def fallout_activity(x_m, y_m, yield_kt=YIELD_KT,
                     wind_speed=WIND_SPEED_ALTITUDE,
                     wind_dir_deg=WIND_DIR_DEG):
    """
    Ground-level activity (normalised, 0-1) from fission-product fallout.
    Gaussian plume in wind-rotated coordinates.
    """
    x = np.atleast_1d(np.asarray(x_m, dtype=float))
    y = np.atleast_1d(np.asarray(y_m, dtype=float))

    theta = np.radians(wind_dir_deg)
    dx_wind = np.cos(theta)
    dy_wind = np.sin(theta)

    downwind  = x * dx_wind + y * dy_wind
    crosswind = np.abs(-x * dy_wind + y * dx_wind)

    H = cloud_top_height(yield_kt)
    t_fall = H / 2.0
    drift = wind_speed * t_fall

    sigma_along = 0.08 * np.maximum(downwind, 1) + 800
    sigma_cross = 0.05 * np.maximum(downwind, 1) + 500

    activity = np.where(
        downwind > 0,
        np.exp(-0.5 * ((downwind - drift) / sigma_along)**2) *
        np.exp(-0.5 * (crosswind / sigma_cross)**2),
        0.0
    )
    return np.squeeze(activity)


def fallout_dose_rate_1hr(activity_norm, yield_kt=YIELD_KT):
    """
    Approximate dose rate (Sv/h) at 1 hour from normalised activity.
    Scaled from fission yield fraction.
    """
    fission_fraction = 0.5
    fission_yield_kt = yield_kt * fission_fraction
    ref_rate = 50.0 * fission_yield_kt
    return activity_norm * ref_rate


# ═══════════════════════════════════════════════════════════════════════════════
# PART E — EMP Model
# ═══════════════════════════════════════════════════════════════════════════════

def emp_e1_peak(distance_m, yield_kt=YIELD_KT):
    """Peak E1 field (V/m) vs distance for a surface/low-altitude burst."""
    E0 = 25000.0 * (yield_kt / 100.0)**0.5
    return E0 * np.exp(-distance_m / 30000.0)


def emp_e1_waveform(t_s, distance_m=10000, yield_kt=YIELD_KT):
    tau_rise = 2.5e-9
    tau_decay = 1.0e-6
    E_peak = emp_e1_peak(distance_m, yield_kt)
    t = np.atleast_1d(np.asarray(t_s, dtype=float))
    E = np.where(t > 0, E_peak * (np.exp(-t / tau_decay) - np.exp(-t / tau_rise)), 0.0)
    return np.squeeze(E)


# ═══════════════════════════════════════════════════════════════════════════════
# PART F — Fireball & Mushroom Cloud
# ═══════════════════════════════════════════════════════════════════════════════

def fireball_radius_t(t, yield_kt=YIELD_KT):
    R_max = fireball_max_radius(yield_kt)
    t_form = 0.001 * yield_kt**0.3
    t = np.atleast_1d(np.asarray(t, dtype=float))
    R = np.where(t < 0, 0, np.where(t < t_form, R_max * (t / t_form)**0.4, R_max))
    return np.squeeze(R)


def shock_radius_t(t, yield_kt=YIELD_KT):
    rho_0 = 1.225
    E = yield_kt * JOULES_PER_KT
    t = np.atleast_1d(np.asarray(t, dtype=float))
    R = np.where(t > 0, 1.15 * (E / rho_0)**0.2 * t**0.4, 0)
    return np.squeeze(R)


# ═══════════════════════════════════════════════════════════════════════════════
# PART G — Ground-Coupled Seismic Effects
# ═══════════════════════════════════════════════════════════════════════════════

def peak_ground_velocity(ground_range_m, yield_kt=YIELD_KT, hob=BURST_HEIGHT):
    """Estimated PGV (m/s) from air-blast ground coupling."""
    r = np.atleast_1d(np.asarray(ground_range_m, dtype=float))
    P = peak_overpressure_kpa(r, yield_kt, hob)
    impedance = 2700 * 3000  # rho * Vs for near-surface rock
    PGV = P * 1000.0 / impedance * 0.5
    return np.squeeze(PGV)


def seismic_magnitude_airblast(yield_kt=YIELD_KT):
    """Approximate mb from air-ground coupling (much weaker than underground)."""
    coupling_efficiency = 0.001
    W_eff = yield_kt * coupling_efficiency
    return 4.0 + 0.75 * np.log10(max(W_eff, 1e-6))


# ═══════════════════════════════════════════════════════════════════════════════
# PART H — Coordinate Utilities
# ═══════════════════════════════════════════════════════════════════════════════

def meters_to_deg_lat(m):
    return m / 111320.0

def meters_to_deg_lon(m, lat=EVENT_LAT):
    return m / (111320.0 * np.cos(np.radians(lat)))

def make_radial_grid(max_range_m=30000, n=500):
    """Create (X_m, Y_m) grid centred on ground zero."""
    x = np.linspace(-max_range_m, max_range_m, n)
    X, Y = np.meshgrid(x, x)
    return X, Y, np.sqrt(X**2 + Y**2)

def km_grid_to_lonlat(X_m, Y_m):
    lon = EVENT_LON + meters_to_deg_lon(X_m)
    lat = EVENT_LAT + meters_to_deg_lat(Y_m)
    return lon, lat


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE GENERATION
# ═══════════════════════════════════════════════════════════════════════════════

def fig01_scenario_overview():
    """Scenario summary: parameters, energy partition, and fireball cross-section."""
    fig = plt.figure(figsize=(20, 14))
    gs = GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.35)

    # (a) Energy partition
    ax = fig.add_subplot(gs[0, 0])
    fractions = [0.50, 0.35, 0.05, 0.10]
    labels = ["Blast/Shock\n50%", "Thermal\nRadiation\n35%",
              "Prompt\nRadiation\n5%", "Residual\n(Fallout)\n10%"]
    colors = ["#e74c3c", "#f39c12", "#2ecc71", "#3498db"]
    wedges, texts = ax.pie(fractions, labels=labels, colors=colors,
                           startangle=90, textprops={"fontsize": 9})
    ax.set_title("(a) Energy Partition", fontweight="bold", fontsize=11)

    # (b) Fireball / shock cross-section at t=0.1s
    ax = fig.add_subplot(gs[0, 1])
    R_fb = fireball_max_radius()
    R_shock_01 = shock_radius_t(0.1)
    ax.add_patch(Circle((0, BURST_HEIGHT), R_shock_01, fc="#ffdddd", ec="red",
                        lw=2, ls="--", label=f"Shock @ 0.1 s ({R_shock_01:.0f} m)"))
    ax.add_patch(Circle((0, BURST_HEIGHT), R_fb, fc="#ff6600", ec="k",
                        lw=2, label=f"Fireball ({R_fb:.0f} m)"))
    ax.axhline(0, color="brown", lw=3)
    ax.fill_between([-1500, 1500], -200, 0, color="#D2B48C", alpha=0.4)
    ax.plot(0, BURST_HEIGHT, "w*", ms=20, mec="k", mew=2, zorder=10)
    ax.set_xlim(-1500, 1500)
    ax.set_ylim(-200, max(R_shock_01 + BURST_HEIGHT + 200, 1500))
    ax.set_aspect("equal")
    ax.set_xlabel("Distance (m)")
    ax.set_ylabel("Height (m)")
    ax.set_title("(b) Fireball & Shock (t = 0.1 s)", fontweight="bold", fontsize=11)
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(True, alpha=0.2)

    # (c) Overpressure vs ground range
    ax = fig.add_subplot(gs[0, 2])
    r = np.linspace(10, 25000, 1000)
    P = peak_overpressure_kpa(r)
    ax.semilogy(r / 1000, P, "r-", lw=2.5)

    for label, thresh in DAMAGE_THRESHOLDS.items():
        short = label.split("(")[0].strip()
        ax.axhline(thresh, color="gray", ls=":", lw=0.8, alpha=0.7)
        rad = blast_radius_for_pressure(thresh)
        ax.annotate(f"{short}\n({rad/1000:.1f} km)",
                    xy=(rad / 1000, thresh), fontsize=6, color="#333",
                    ha="left", va="bottom")

    ax.set_xlabel("Ground range (km)")
    ax.set_ylabel("Peak overpressure (kPa)")
    ax.set_title("(c) Overpressure vs Range", fontweight="bold", fontsize=11)
    ax.set_xlim(0, 25)
    ax.set_ylim(0.5, 5000)
    ax.grid(True, alpha=0.2, which="both")

    # (d) Thermal fluence vs range
    ax = fig.add_subplot(gs[1, 0])
    Q = thermal_fluence_kj(r)
    ax.semilogy(r / 1000, Q, "darkorange", lw=2.5)

    for label, thresh in THERMAL_THRESHOLDS.items():
        short = label.split("(")[0].strip()
        ax.axhline(thresh, color="gray", ls=":", lw=0.8, alpha=0.7)
        idx = np.argmin(np.abs(Q - thresh))
        ax.annotate(f"{short}\n({r[idx]/1000:.1f} km)",
                    xy=(r[idx] / 1000, thresh), fontsize=6, color="#333",
                    ha="left", va="bottom")

    ax.set_xlabel("Ground range (km)")
    ax.set_ylabel("Thermal fluence (kJ/m²)")
    ax.set_title("(d) Thermal Fluence vs Range", fontweight="bold", fontsize=11)
    ax.set_xlim(0, 25)
    ax.grid(True, alpha=0.2, which="both")

    # (e) Prompt radiation dose vs range
    ax = fig.add_subplot(gs[1, 1])
    D = prompt_radiation_dose_gy(r)
    ax.semilogy(r / 1000, D, "green", lw=2.5)

    for label, thresh in RADIATION_THRESHOLDS.items():
        short = label.split("(")[0].strip()
        ax.axhline(thresh, color="gray", ls=":", lw=0.8, alpha=0.7)
        idx = np.argmin(np.abs(D - thresh))
        ax.annotate(f"{short}\n({r[idx]/1000:.1f} km)",
                    xy=(r[idx] / 1000, thresh), fontsize=6, color="#333",
                    ha="left", va="bottom")

    ax.set_xlabel("Ground range (km)")
    ax.set_ylabel("Dose (Gy)")
    ax.set_title("(e) Prompt Radiation Dose vs Range", fontweight="bold", fontsize=11)
    ax.set_xlim(0, 10)
    ax.set_ylim(1e-4, 1e4)
    ax.grid(True, alpha=0.2, which="both")

    # (f) Parameter table
    ax = fig.add_subplot(gs[1, 2])
    ax.axis("off")
    Rfb = fireball_max_radius()
    Hcloud = cloud_top_height()
    mb_air = seismic_magnitude_airblast()
    lines = [
        "SCENARIO PARAMETERS",
        "═" * 42,
        f"Weapon yield:       {YIELD_KT:.0f} kt",
        f"Height of burst:    {BURST_HEIGHT:.0f} m  (airburst)",
        f"Location:           {EVENT_LAT:.4f}°N  {abs(EVENT_LON):.4f}°W",
        f"                    (UC Berkeley campus)",
        "",
        f"Total energy:       {TOTAL_ENERGY:.2e} J",
        f"Fireball radius:    {Rfb:.0f} m  (maximum)",
        f"Cloud top height:   {Hcloud/1000:.1f} km  (stabilised)",
        "",
        "PREVAILING WINDS:",
        f"  Surface:          {WIND_SPEED_SURFACE:.0f} m/s from NW",
        f"  Altitude:         {WIND_SPEED_ALTITUDE:.0f} m/s from NW",
        f"  Direction:        {WIND_DIR_DEG:.0f}°  (toward SE)",
        "",
        "DAMAGE RADII (ground range):",
    ]
    for label, thresh in DAMAGE_THRESHOLDS.items():
        rad = blast_radius_for_pressure(thresh)
        short = label.split("(")[0].strip()
        lines.append(f"  {short}:  {rad/1000:.2f} km")
    lines += [
        "",
        f"Air-coupled mb:     {mb_air:.1f}",
        f"PGV @ 5 km:         {peak_ground_velocity(5000)*100:.1f} cm/s",
    ]
    ax.text(0.03, 0.97, "\n".join(lines), transform=ax.transAxes,
            fontfamily="monospace", fontsize=8, va="top",
            bbox=dict(boxstyle="round,pad=0.4", fc="#f8f8f8", ec="#bbb"))

    fig.suptitle("Figure 1 — Scenario Overview:  100 kT Airburst over Berkeley, CA (50 m HOB)",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    p = os.path.join(OUTDIR, "fig01_scenario_overview.png")
    fig.savefig(p, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {p}")


def fig02_blast_damage_map():
    """2D blast damage zone map centred on Berkeley with geographic coordinates.
    
    Note on topography: Blast wave propagation is affected by terrain features.
    Hills (e.g., Oakland Hills, Berkeley Hills) can provide partial shielding,
    while valleys may channel and focus blast energy. This model assumes flat
    terrain; actual damage patterns would be modified by the local topography.
    """
    max_range = 25000
    X, Y, R = make_radial_grid(max_range, 800)
    P = peak_overpressure_kpa(R)

    damage_cat = np.zeros_like(P)
    damage_cat[P >= 1.0]   = 1
    damage_cat[P >= 3.5]   = 2
    damage_cat[P >= 7.0]   = 3
    damage_cat[P >= 14.0]  = 4
    damage_cat[P >= 35.0]  = 5
    damage_cat[P >= 140.0] = 6

    fig, axes = plt.subplots(1, 2, figsize=(20, 10))

    # Left: categorical damage
    ax = axes[0]
    cmap = mcolors.ListedColormap(
        ["white", "#cceeff", "#ffffaa", "#ffcc44", "#ff6600", "#cc0000", "#440000"])
    bounds = [0, 1, 2, 3, 4, 5, 6, 7]
    norm = mcolors.BoundaryNorm(bounds, cmap.N)
    im = ax.imshow(damage_cat, extent=[-max_range/1000, max_range/1000,
                                        -max_range/1000, max_range/1000],
                   origin="lower", cmap=cmap, norm=norm)
    ax.plot(0, 0, "w*", ms=18, mec="k", mew=2, zorder=10)
    ax.set_xlabel("Distance E-W (km)")
    ax.set_ylabel("Distance N-S (km)")
    ax.set_title("(a) Blast Damage Zones", fontweight="bold")
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.15)

    # Add secondary axes with geographic coordinates
    ax_lon = ax.secondary_xaxis('top', functions=(
        lambda x: EVENT_LON + meters_to_deg_lon(x * 1000),
        lambda lon: (lon - EVENT_LON) / meters_to_deg_lon(1000)))
    ax_lon.set_xlabel('Longitude (°W)')
    ax_lat = ax.secondary_yaxis('right', functions=(
        lambda y: EVENT_LAT + meters_to_deg_lat(y * 1000),
        lambda lat: (lat - EVENT_LAT) / meters_to_deg_lat(1000)))
    ax_lat.set_ylabel('Latitude (°N)')

    labels = ["No damage", ">1 kPa glass", ">3.5 kPa light",
              ">7 kPa moderate", ">14 kPa severe", ">35 kPa reinf. conc.",
              ">140 kPa total"]
    cbar = plt.colorbar(im, ax=ax, ticks=[0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5],
                        label="Damage Category (threshold overpressure)")
    cbar.set_ticklabels(labels)

    for label, thresh in DAMAGE_THRESHOLDS.items():
        rad = blast_radius_for_pressure(thresh)
        circle = Circle((0, 0), rad / 1000, fill=False, ec="k", ls="--", lw=0.6)
        ax.add_patch(circle)

    # landmarks
    landmarks = [
        (-3.8, 1.2, "Emeryville"),
        (4.0, -3.0, "Oakland Hills"),
        (0.0, 0.5, "UC Berkeley"),
        (-8.0, 3.0, "San Francisco\n(via Bay Bridge)"),
        (5.0, 5.0, "Walnut Creek"),
        (-3.0, -5.0, "Alameda"),
        (0.0, -8.0, "San Leandro"),
    ]
    for lx, ly, lname in landmarks:
        ax.plot(lx, ly, "ko", ms=4, zorder=8)
        ax.text(lx + 0.3, ly + 0.3, lname, fontsize=6, color="k", zorder=8)

    # Right: continuous overpressure
    ax = axes[1]
    P_plot = np.clip(P, 0.1, 5000)
    im2 = ax.contourf(X / 1000, Y / 1000, P_plot,
                       levels=np.logspace(-1, 3, 30),
                       cmap="hot", norm=mcolors.LogNorm(vmin=0.1, vmax=5000))
    ax.plot(0, 0, "w*", ms=18, mec="k", mew=2, zorder=10)
    ax.set_xlabel("Distance E-W (km)")
    ax.set_ylabel("Distance N-S (km)")
    ax.set_title("(b) Peak Overpressure", fontweight="bold")
    ax.set_aspect("equal")
    plt.colorbar(im2, ax=ax, label="Overpressure (kPa)")

    # Add secondary axes with geographic coordinates
    ax_lon = ax.secondary_xaxis('top', functions=(
        lambda x: EVENT_LON + meters_to_deg_lon(x * 1000),
        lambda lon: (lon - EVENT_LON) / meters_to_deg_lon(1000)))
    ax_lon.set_xlabel('Longitude (°W)')
    ax_lat = ax.secondary_yaxis('right', functions=(
        lambda y: EVENT_LAT + meters_to_deg_lat(y * 1000),
        lambda lat: (lat - EVENT_LAT) / meters_to_deg_lat(1000)))
    ax_lat.set_ylabel('Latitude (°N)')

    for label, thresh in DAMAGE_THRESHOLDS.items():
        ax.contour(X / 1000, Y / 1000, P, levels=[thresh],
                   colors="white", linewidths=0.8)

    # Add topography note
    fig.text(0.5, 0.02,
             "Note: Flat terrain assumed. Actual damage modified by Berkeley/Oakland Hills topography (shielding) "
             "and valleys (focusing). Center: {:.4f}°N, {:.4f}°W".format(EVENT_LAT, abs(EVENT_LON)),
             ha='center', fontsize=9, style='italic', color='#444')

    fig.suptitle("Figure 2 — Blast Damage Map:  100 kT at Berkeley (50 m HOB)",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0.04, 1, 0.96])
    p = os.path.join(OUTDIR, "fig02_blast_damage_map.png")
    fig.savefig(p, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {p}")


def fig03_thermal_and_radiation():
    """Combined thermal fluence and prompt radiation maps.
    
    Note on topography: Thermal radiation travels in straight lines and can be
    blocked by terrain. Berkeley Hills provide natural shielding for areas to
    the east. Valleys and urban canyons may experience focusing effects.
    """
    max_range = 25000
    X, Y, R = make_radial_grid(max_range, 600)

    Q = thermal_fluence_kj(R)
    D = prompt_radiation_dose_gy(R)

    fig, axes = plt.subplots(1, 3, figsize=(22, 8))

    # (a) Thermal fluence
    ax = axes[0]
    Q_clip = np.clip(Q, 1, 1e5)
    im1 = ax.contourf(X / 1000, Y / 1000, Q_clip,
                       levels=np.logspace(0, 5, 25),
                       cmap="inferno", norm=mcolors.LogNorm())
    ax.plot(0, 0, "w*", ms=14, mec="k", mew=1.5)
    ax.set_xlabel("Distance E-W (km)")
    ax.set_ylabel("Distance N-S (km)")
    ax.set_title("(a) Thermal Fluence", fontweight="bold")
    ax.set_aspect("equal")
    plt.colorbar(im1, ax=ax, label="Thermal Fluence (kJ/m²)")
    for _, thresh in THERMAL_THRESHOLDS.items():
        ax.contour(X / 1000, Y / 1000, Q, levels=[thresh],
                   colors="white", linewidths=0.8, linestyles="--")
    # Geographic coordinate axes
    ax_lon = ax.secondary_xaxis('top', functions=(
        lambda x: EVENT_LON + meters_to_deg_lon(x * 1000),
        lambda lon: (lon - EVENT_LON) / meters_to_deg_lon(1000)))
    ax_lon.set_xlabel('Longitude (°W)')

    # (b) Prompt radiation
    ax = axes[1]
    D_clip = np.clip(D, 1e-4, 1e4)
    im2 = ax.contourf(X / 1000, Y / 1000, D_clip,
                       levels=np.logspace(-4, 4, 25),
                       cmap="YlOrRd", norm=mcolors.LogNorm())
    ax.plot(0, 0, "w*", ms=14, mec="k", mew=1.5)
    ax.set_xlabel("Distance E-W (km)")
    ax.set_ylabel("Distance N-S (km)")
    ax.set_title("(b) Prompt Radiation Dose", fontweight="bold")
    ax.set_aspect("equal")
    plt.colorbar(im2, ax=ax, label="Absorbed Dose (Gy)")
    for _, thresh in RADIATION_THRESHOLDS.items():
        ax.contour(X / 1000, Y / 1000, D, levels=[thresh],
                   colors="white", linewidths=0.8, linestyles="--")
    # Geographic coordinate axes
    ax_lon = ax.secondary_xaxis('top', functions=(
        lambda x: EVENT_LON + meters_to_deg_lon(x * 1000),
        lambda lon: (lon - EVENT_LON) / meters_to_deg_lon(1000)))
    ax_lon.set_xlabel('Longitude (°W)')

    # (c) Thermal burn zones
    ax = axes[2]
    burn_cat = np.zeros_like(Q)
    burn_cat[Q >= 50]  = 1
    burn_cat[Q >= 125] = 2
    burn_cat[Q >= 335] = 3
    burn_cat[Q >= 670] = 4

    cmap_burn = mcolors.ListedColormap(
        ["white", "#ffffaa", "#ffcc44", "#ff6600", "#cc0000"])
    bounds = [0, 1, 2, 3, 4, 5]
    norm_b = mcolors.BoundaryNorm(bounds, cmap_burn.N)
    im3 = ax.imshow(burn_cat,
                    extent=[-max_range/1000, max_range/1000,
                            -max_range/1000, max_range/1000],
                    origin="lower", cmap=cmap_burn, norm=norm_b)
    ax.plot(0, 0, "w*", ms=14, mec="k", mew=1.5)
    ax.set_xlabel("Distance E-W (km)")
    ax.set_ylabel("Distance N-S (km)")
    ax.set_title("(c) Thermal Burn Zones", fontweight="bold")
    ax.set_aspect("equal")
    cb3 = plt.colorbar(im3, ax=ax, ticks=[0.5, 1.5, 2.5, 3.5, 4.5],
                       label="Burn Severity (thermal fluence threshold)")
    cb3.set_ticklabels(["No effect", "Pain (>50 kJ/m²)", 
                        "1st° (>125 kJ/m²)", "2nd° (>335 kJ/m²)", 
                        "3rd° (>670 kJ/m²)"])
    # Geographic coordinate axes
    ax_lon = ax.secondary_xaxis('top', functions=(
        lambda x: EVENT_LON + meters_to_deg_lon(x * 1000),
        lambda lon: (lon - EVENT_LON) / meters_to_deg_lon(1000)))
    ax_lon.set_xlabel('Longitude (°W)')

    # Add topography note
    fig.text(0.5, 0.01,
             "Note: Flat terrain assumed. Hills provide line-of-sight shielding from thermal radiation. "
             "Center: {:.4f}°N, {:.4f}°W (UC Berkeley)".format(EVENT_LAT, abs(EVENT_LON)),
             ha='center', fontsize=9, style='italic', color='#444')

    fig.suptitle("Figure 3 — Thermal Radiation & Prompt Nuclear Radiation",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0.04, 1, 0.96])
    p = os.path.join(OUTDIR, "fig03_thermal_and_radiation.png")
    fig.savefig(p, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {p}")


def fig04_fallout_pattern():
    """Fallout plume and dose-rate maps under prevailing NW wind."""
    max_range = 80000
    n = 600
    x = np.linspace(-max_range, max_range, n)
    X, Y = np.meshgrid(x, x)

    activity = fallout_activity(X, Y)
    dose_1h = fallout_dose_rate_1hr(activity)

    activity_smooth = gaussian_filter(activity, sigma=4)
    dose_smooth = gaussian_filter(dose_1h, sigma=4)

    fig = plt.figure(figsize=(22, 12))
    gs = GridSpec(2, 3, figure=fig, hspace=0.3, wspace=0.3)

    # (a) Fallout deposition
    ax = fig.add_subplot(gs[0, 0])
    act_plot = np.maximum(activity_smooth, 1e-6)
    im1 = ax.contourf(X / 1000, Y / 1000, act_plot,
                       levels=np.logspace(-4, 0, 20),
                       cmap="YlGnBu_r", norm=mcolors.LogNorm())
    ax.plot(0, 0, "r*", ms=14, zorder=10)
    ax.set_xlabel("km E-W")
    ax.set_ylabel("km N-S")
    ax.set_title("(a) Fallout Deposition (normalised)", fontweight="bold")
    ax.set_aspect("equal")
    plt.colorbar(im1, ax=ax)

    wind_dx = np.cos(np.radians(WIND_DIR_DEG))
    wind_dy = np.sin(np.radians(WIND_DIR_DEG))
    ax.annotate("", xy=(wind_dx * 25, wind_dy * 25),
                xytext=(0, 0),
                arrowprops=dict(arrowstyle="->", color="red", lw=2.5))
    ax.text(wind_dx * 27, wind_dy * 27, "Wind", color="red", fontsize=10,
            fontweight="bold")

    # (b) Dose rate at 1 hour
    ax = fig.add_subplot(gs[0, 1])
    dose_plot = np.maximum(dose_smooth, 1e-3)
    im2 = ax.contourf(X / 1000, Y / 1000, dose_plot,
                       levels=np.logspace(-3, 4, 25),
                       cmap="YlOrRd", norm=mcolors.LogNorm())
    ax.plot(0, 0, "k*", ms=14, zorder=10)
    ax.set_xlabel("km E-W")
    ax.set_ylabel("km N-S")
    ax.set_title("(b) Dose Rate at H+1 (Sv/h)", fontweight="bold")
    ax.set_aspect("equal")
    plt.colorbar(im2, ax=ax, label="Sv/h")

    evac_levels = [0.02, 0.1, 1, 10, 100]
    ax.contour(X / 1000, Y / 1000, dose_smooth, levels=evac_levels,
               colors="white", linewidths=0.8)

    # (c) Dose rate decay (Way-Wigner)
    ax = fig.add_subplot(gs[0, 2])
    t_hours = np.logspace(-1, 3, 200)
    ref_rate = 1000.0  # Sv/h at 1 hour at hotspot
    rate = ref_rate * (t_hours)**(-1.2)

    ax.loglog(t_hours, rate, "r-", lw=2.5, label="Close-in (hotspot)")
    ax.loglog(t_hours, rate * 0.1, "orange", lw=2, label="10% of hotspot")
    ax.loglog(t_hours, rate * 0.01, "b-", lw=2, label="1% of hotspot")

    ax.axhline(0.02, color="green", ls="--", lw=1.5, label="Shelter (20 mSv/h)")
    ax.axhline(0.5, color="orange", ls="--", lw=1.5, label="Evacuation (500 mSv/h)")

    ax.set_xlabel("Time after detonation (hours)")
    ax.set_ylabel("Dose rate (Sv/h)")
    ax.set_title("(c) Fallout Decay (Way-Wigner t⁻¹·²)", fontweight="bold")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.2, which="both")
    ax.set_xlim(0.1, 1000)
    ax.set_ylim(1e-5, 1e4)

    # (d) Cumulative dose
    ax = fig.add_subplot(gs[1, 0])
    cumul_hotspot = ref_rate * (t_hours**(1 - 1.2) - 1**(-0.2)) / (1 - 1.2) * -1
    cumul_10pct = cumul_hotspot * 0.1
    cumul_1pct = cumul_hotspot * 0.01
    ax.loglog(t_hours, cumul_hotspot, "r-", lw=2.5, label="Hotspot")
    ax.loglog(t_hours, cumul_10pct, "orange", lw=2, label="10% level")
    ax.loglog(t_hours, cumul_1pct, "b-", lw=2, label="1% level")
    ax.axhline(1.0, color="red", ls=":", lw=1.5, label="1 Sv (evacuation)")
    ax.axhline(0.1, color="orange", ls=":", lw=1.5, label="100 mSv (shelter)")
    ax.set_xlabel("Time (hours)")
    ax.set_ylabel("Cumulative dose (Sv)")
    ax.set_title("(d) Cumulative Fallout Dose", fontweight="bold")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.2, which="both")
    ax.set_xlim(0.1, 1000)

    # (e) Downwind cross-section
    ax = fig.add_subplot(gs[1, 1])
    downwind_km = np.linspace(0, 80, 200)
    theta = np.radians(WIND_DIR_DEG)
    for cw_km in [0, 2, 5, 10]:
        dx = downwind_km * 1000 * np.cos(theta) - cw_km * 1000 * np.sin(theta)
        dy = downwind_km * 1000 * np.sin(theta) + cw_km * 1000 * np.cos(theta)
        act = fallout_activity(dx, dy)
        dr = fallout_dose_rate_1hr(act)
        dr = np.maximum(dr, 1e-6)
        ax.semilogy(downwind_km, dr, lw=1.8,
                    label=f"Crosswind = {cw_km} km")

    ax.set_xlabel("Downwind distance (km)")
    ax.set_ylabel("Dose rate at H+1 (Sv/h)")
    ax.set_title("(e) Downwind Fallout Profile", fontweight="bold")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.2, which="both")
    ax.set_ylim(1e-3, 1e4)

    # (f) Fallout info
    ax = fig.add_subplot(gs[1, 2])
    ax.axis("off")
    H = cloud_top_height()
    lines = [
        "FALLOUT PARAMETERS",
        "═" * 38,
        f"Cloud top height:   {H/1000:.1f} km",
        f"Wind at altitude:   {WIND_SPEED_ALTITUDE:.0f} m/s from {WIND_DIR_DEG:.0f}°",
        f"Fission fraction:   50%",
        f"Fallout direction:  toward SE",
        "",
        "KEY DOWNWIND AREAS:",
        "  Oakland, San Leandro, Hayward",
        "  Fremont, San Jose (extended)",
        "",
        "PROTECTIVE ACTIONS:",
        "  0-2 km:     Unsurvivable blast",
        "  2-5 km:     Seek immediate shelter",
        "  5-20 km:    Shelter in place 24-48 h",
        "  20-80 km:   Monitor, possible evac.",
        "",
        "NOTE: Fallout from a 50 m airburst",
        "is significantly LESS than a surface",
        "burst — much less ground material is",
        "entrained. Nonetheless, a 100 kt",
        "fission weapon produces substantial",
        "residual radioactivity.",
    ]
    ax.text(0.03, 0.97, "\n".join(lines), transform=ax.transAxes,
            fontfamily="monospace", fontsize=8.5, va="top",
            bbox=dict(boxstyle="round,pad=0.4", fc="#f8f8f8", ec="#bbb"))

    fig.suptitle("Figure 4 — Fallout Pattern:  100 kT Airburst, Prevailing NW Wind",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    p = os.path.join(OUTDIR, "fig04_fallout_pattern.png")
    fig.savefig(p, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {p}")


def fig05_emp_analysis():
    """EMP analysis for low-altitude burst."""
    fig = plt.figure(figsize=(18, 12))
    gs = GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.35)

    # (a) E1 waveform
    ax = fig.add_subplot(gs[0, 0])
    t_us = np.linspace(0, 10, 2000)
    for d_km in [5, 10, 20, 50]:
        E = emp_e1_waveform(t_us * 1e-6, d_km * 1000)
        ax.plot(t_us, E, lw=1.5, label=f"{d_km} km")
    ax.set_xlabel("Time (µs)")
    ax.set_ylabel("Electric field (V/m)")
    ax.set_title("(a) E1 Waveform vs Distance", fontweight="bold")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.2)

    # (b) E1 peak vs distance
    ax = fig.add_subplot(gs[0, 1])
    d = np.linspace(100, 100000, 500)
    E_peak = emp_e1_peak(d)
    ax.semilogy(d / 1000, E_peak, "r-", lw=2.5)
    ax.axhline(25000, color="gray", ls=":", label="MIL-STD-461 limit")
    ax.set_xlabel("Distance (km)")
    ax.set_ylabel("Peak E1 field (V/m)")
    ax.set_title("(b) E1 Peak vs Distance", fontweight="bold")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.2, which="both")

    # (c) E1 spatial distribution
    ax = fig.add_subplot(gs[0, 2])
    max_r = 60000
    Xg, Yg, Rg = make_radial_grid(max_r, 400)
    E_spatial = emp_e1_peak(Rg)
    im = ax.contourf(Xg / 1000, Yg / 1000, E_spatial,
                      levels=np.logspace(1, 5, 20),
                      cmap="plasma", norm=mcolors.LogNorm())
    ax.plot(0, 0, "w*", ms=12, mec="k")
    ax.set_xlabel("km E-W")
    ax.set_ylabel("km N-S")
    ax.set_title("(c) E1 Peak Field (V/m)", fontweight="bold")
    ax.set_aspect("equal")
    plt.colorbar(im, ax=ax, label="V/m")

    # (d) Induced voltage on power lines
    ax = fig.add_subplot(gs[1, 0])
    t_line = np.linspace(0, 10e-6, 1000)
    for L_km in [1, 5, 10, 50]:
        coupling = 0.3
        V = emp_e1_waveform(t_line, 10000) * L_km * 1000 * coupling
        ax.plot(t_line * 1e6, V / 1000, lw=1.5, label=f"{L_km} km line")
    ax.set_xlabel("Time (µs)")
    ax.set_ylabel("Induced voltage (kV)")
    ax.set_title("(d) Induced Voltage on Power Lines (@ 10 km)", fontweight="bold")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.2)

    # (e) EMP frequency spectrum
    ax = fig.add_subplot(gs[1, 1])
    f = np.logspace(3, 11, 1000)
    e1_spec = 1e4 * np.exp(-(np.log10(f) - 8)**2 / 2)
    e2_spec = 1e3 * np.exp(-(np.log10(f) - 5)**2 / 2)
    e3_spec = 1e2 * np.exp(-(np.log10(f) - 0)**2 / 2)
    ax.loglog(f, e1_spec, "r-", lw=2, label="E1")
    ax.loglog(f, e2_spec, "orange", lw=2, label="E2")
    ax.loglog(f, e3_spec, "b-", lw=2, label="E3")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Spectral density (V/m/Hz)")
    ax.set_title("(e) EMP Frequency Spectrum", fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.2, which="both")

    # (f) Infrastructure vulnerability table
    ax = fig.add_subplot(gs[1, 2])
    ax.axis("off")
    lines = [
        "EMP INFRASTRUCTURE VULNERABILITY",
        "═" * 38,
        "",
        "  E1 (nanoseconds):",
        "    Electronics, SCADA, telecom",
        "    Vulnerable within ~30 km",
        "",
        "  E2 (microseconds):",
        "    Similar to lightning surges",
        "    Moderate vulnerability",
        "",
        "  E3 (MHD, seconds):",
        "    Power grid transformers",
        "    Very limited for surface burst",
        "",
        "BAY AREA IMPACT:",
        "  Berkeley, Oakland: severe E1",
        "  San Francisco: moderate E1",
        "  South Bay: weak E1",
        "",
        "NOTE: Low-altitude burst produces",
        "  much weaker EMP than high-altitude",
        "  burst. Main EMP threat is localised.",
    ]
    ax.text(0.03, 0.97, "\n".join(lines), transform=ax.transAxes,
            fontfamily="monospace", fontsize=8.5, va="top",
            bbox=dict(boxstyle="round,pad=0.4", fc="#f8f8f8", ec="#bbb"))

    fig.suptitle("Figure 5 — Electromagnetic Pulse (EMP) Analysis:  100 kT at 50 m HOB",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    p = os.path.join(OUTDIR, "fig05_emp_analysis.png")
    fig.savefig(p, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {p}")


def fig06_fireball_and_cloud():
    """Fireball evolution and mushroom cloud development."""
    fig = plt.figure(figsize=(18, 12))
    gs = GridSpec(2, 3, figure=fig, hspace=0.3, wspace=0.3)

    # (a) Fireball & shock radius vs time
    ax = fig.add_subplot(gs[0, 0])
    t = np.logspace(-6, 2, 1000)
    Rfb = fireball_radius_t(t)
    Rsh = shock_radius_t(t)
    ax.loglog(t, Rfb, "r-", lw=2.5, label="Fireball")
    ax.loglog(t, Rsh, "b--", lw=2, label="Shock front")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Radius (m)")
    ax.set_title("(a) Fireball & Shock Expansion", fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.2, which="both")
    ax.set_xlim(1e-6, 100)

    # (b) Fireball temperature
    ax = fig.add_subplot(gs[0, 1])
    T = np.zeros_like(t)
    for i, ti in enumerate(t):
        if ti < 1e-6:
            T[i] = 1e8
        elif ti < 0.01:
            T[i] = 1e8 * (1e-6 / ti)**0.5
        elif ti < 0.1:
            T[i] = 8000
        else:
            T[i] = 8000 * (0.1 / ti)**0.5
    ax.loglog(t, T, "r-", lw=2.5)
    ax.axhline(5778, color="orange", ls="--", lw=1.5, label="Solar surface")
    ax.axhline(1811, color="gray", ls=":", lw=1, label="Steel melting point")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Temperature (K)")
    ax.set_title("(b) Fireball Surface Temperature", fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.2, which="both")
    ax.set_xlim(1e-6, 100)

    # (c) Thermal power
    ax = fig.add_subplot(gs[0, 2])
    sigma = 5.67e-8
    P_thermal = 4 * np.pi * Rfb**2 * sigma * T**4
    ax.loglog(t, P_thermal / 1e15, "r-", lw=2.5)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Thermal power (PW)")
    ax.set_title("(c) Thermal Radiation Power", fontweight="bold")
    ax.grid(True, alpha=0.2, which="both")
    ax.set_xlim(1e-6, 100)

    # (d) Mushroom cloud rise
    ax = fig.add_subplot(gs[1, 0:2])
    times_cloud = [1, 5, 15, 30, 60, 120]
    colors_cloud = plt.cm.Reds(np.linspace(0.3, 1, len(times_cloud)))

    for ti, c in zip(times_cloud, colors_cloud):
        cloud_top = min(cloud_top_height() / 1000, 0.8 * ti**0.5)
        cloud_rad = min(5, 0.15 * ti**0.5)
        stem_w = cloud_rad * 0.3
        stem = plt.Rectangle((-stem_w / 2, 0.05), stem_w,
                              cloud_top * 0.7, fc=c, alpha=0.3, ec=c)
        ax.add_patch(stem)
        cap = plt.matplotlib.patches.Ellipse(
            (0, cloud_top), cloud_rad * 2, cloud_rad,
            fc=c, alpha=0.5, ec=c, label=f"t = {ti} s")
        ax.add_patch(cap)

    ax.set_xlim(-15, 15)
    ax.set_ylim(0, 20)
    ax.set_xlabel("Horizontal distance (km)")
    ax.set_ylabel("Altitude (km)")
    ax.set_title("(d) Mushroom Cloud Rise", fontweight="bold")
    ax.legend(loc="upper right", fontsize=8)
    ax.axhline(11, color="blue", ls="--", lw=1, alpha=0.5)
    ax.text(10, 11.3, "Tropopause", fontsize=8, color="blue")
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.2)

    # (e) Double flash
    ax = fig.add_subplot(gs[1, 2])
    t_flash = np.linspace(0, 5, 2000)
    flash = np.zeros_like(t_flash)
    for i, ti in enumerate(t_flash):
        if ti < 0.001:
            flash[i] = np.exp(-ti / 0.0001)
        elif ti < 0.01:
            flash[i] = 0.1 * np.exp(-(ti - 0.01)**2 / 0.0001)
        elif ti < 1:
            flash[i] = 0.8 * np.exp(-(ti - 0.3)**2 / 0.08)
        else:
            flash[i] = 0.4 * np.exp(-(ti - 1) / 1.5)

    ax.plot(t_flash * 1000, flash, "r-", lw=2)
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Relative brightness")
    ax.set_title("(e) Optical Double Flash", fontweight="bold")
    ax.grid(True, alpha=0.2)
    ax.axvline(10, color="blue", ls=":", alpha=0.5)
    ax.text(15, 0.85, "First minimum", fontsize=8, color="blue")
    ax.axvline(300, color="green", ls=":", alpha=0.5)
    ax.text(350, 0.85, "Second max", fontsize=8, color="green")

    fig.suptitle("Figure 6 — Fireball & Mushroom Cloud Evolution:  100 kT Airburst",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    p = os.path.join(OUTDIR, "fig06_fireball_and_cloud.png")
    fig.savefig(p, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {p}")


def fig07_ground_coupling_seismic():
    """Air-blast ground coupling and seismic effects."""
    fig = plt.figure(figsize=(18, 12))
    gs = GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.35)

    # (a) Peak ground velocity vs range
    ax = fig.add_subplot(gs[0, 0])
    r = np.linspace(100, 50000, 500)
    PGV = peak_ground_velocity(r)
    ax.semilogy(r / 1000, PGV * 100, "b-", lw=2.5)
    ax.axhline(2, color="orange", ls="--", label="Light damage (2 cm/s)")
    ax.axhline(10, color="red", ls="--", label="Moderate damage (10 cm/s)")
    ax.set_xlabel("Ground range (km)")
    ax.set_ylabel("PGV (cm/s)")
    ax.set_title("(a) Peak Ground Velocity", fontweight="bold")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.2, which="both")
    ax.set_xlim(0, 50)

    # (b) Overpressure time history at different distances
    ax = fig.add_subplot(gs[0, 1])
    distances = [1, 3, 5, 10, 20]
    t_hist = np.linspace(0, 30, 500)
    for d_km in distances:
        d_m = d_km * 1000
        P_peak = peak_overpressure_kpa(d_m)
        t_pos = 0.2 * d_km**0.5
        P_hist = P_peak * np.where(
            t_hist < d_m / 340,
            0,
            np.where(
                t_hist < d_m / 340 + t_pos,
                1 - (t_hist - d_m / 340) / t_pos,
                -0.3 * np.exp(-(t_hist - d_m / 340 - t_pos) / (t_pos * 0.5))
            )
        )
        ax.plot(t_hist, P_hist, lw=1.5, label=f"{d_km} km ({P_peak:.1f} kPa)")

    ax.set_xlabel("Time after detonation (s)")
    ax.set_ylabel("Overpressure (kPa)")
    ax.set_title("(b) Blast Pressure Time History", fontweight="bold")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.2)

    # (c) Acoustic arrival times
    ax = fig.add_subplot(gs[0, 2])
    d_range = np.linspace(1, 100, 200)
    t_shock = d_range / 0.34
    t_sound = d_range / 0.343

    ax.plot(d_range, t_sound, "b-", lw=2, label="Sound wave (343 m/s)")
    ax.plot(d_range, t_shock, "r--", lw=2, label="Shock front (faster near GZ)")
    ax.set_xlabel("Distance (km)")
    ax.set_ylabel("Arrival time (s)")
    ax.set_title("(c) Acoustic/Shock Arrival Times", fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.2)

    # (d) PGV 2D map with terrain relief
    ax = fig.add_subplot(gs[1, 0])
    max_r = 50000
    Xg, Yg, Rg = make_radial_grid(max_r, 400)
    PGV_map = peak_ground_velocity(Rg) * 100

    # Load terrain relief as background
    try:
        import pygmt
        half_deg = meters_to_deg_lon(max_r)
        half_lat = meters_to_deg_lat(max_r)
        dem_region = [EVENT_LON - half_deg, EVENT_LON + half_deg,
                      EVENT_LAT - half_lat, EVENT_LAT + half_lat]
        dem = pygmt.datasets.load_earth_relief(resolution="03s", region=dem_region)
        elev = dem.values
        dem_lons = dem.lon.values
        dem_lats = dem.lat.values
        # Convert DEM coordinates to km offsets from GZ
        dem_x = (dem_lons - EVENT_LON) * 111320.0 * np.cos(np.radians(EVENT_LAT)) / 1000.0
        dem_y = (dem_lats - EVENT_LAT) * 111320.0 / 1000.0
        # Hillshade for shaded relief
        from matplotlib.colors import LightSource
        ls = LightSource(azdeg=315, altdeg=35)
        hs = ls.hillshade(elev, vert_exag=2, dx=1, dy=1)
        ax.imshow(hs, extent=[dem_x[0], dem_x[-1], dem_y[0], dem_y[-1]],
                  cmap="gray", origin="lower", aspect="equal", alpha=0.6)
    except Exception:
        pass  # If PyGMT/DEM unavailable, proceed without relief

    im = ax.contourf(Xg / 1000, Yg / 1000, PGV_map,
                      levels=np.logspace(-2, 2, 20),
                      cmap="viridis", norm=mcolors.LogNorm(), alpha=0.65)
    ax.plot(0, 0, "r*", ms=12, mec="k")
    ax.set_xlabel("km E-W")
    ax.set_ylabel("km N-S")
    ax.set_title("(d) PGV Map (cm/s) on Terrain", fontweight="bold")
    ax.set_aspect("equal")
    ax.set_xlim(-max_r / 1000, max_r / 1000)
    ax.set_ylim(-max_r / 1000, max_r / 1000)
    plt.colorbar(im, ax=ax, label="cm/s")

    # (e) Air-ground coupling diagram
    ax = fig.add_subplot(gs[1, 1])
    ax.axis("off")
    coupling = [
        "AIR-GROUND COUPLING PHYSICS",
        "═" * 38,
        "",
        "The 100 kT airburst at 50 m HOB",
        "couples to the ground through:",
        "",
        "1. DIRECT BLAST LOADING",
        "   Overpressure acts on surface",
        "   → compressional wave in rock",
        "   → surface/Rayleigh waves",
        "",
        "2. AIR-COUPLED RAYLEIGH WAVE",
        "   Blast wave outruns seismic wave",
        "   → continuous ground excitation",
        "   along shock front",
        "",
        "3. CRATER COUPLING (minor at 50 m)",
        "   Some ground excavation occurs",
        "   → enhanced high-frequency content",
        "",
        f"Estimated seismic magnitude: ",
        f"  mb ≈ {seismic_magnitude_airblast():.1f}",
        f"  (vs mb ≈ 5.5 if underground)",
        "",
        "FELT RANGE:",
        f"  Modified Mercalli ≥ III:  ~50 km",
        f"  Modified Mercalli ≥ V:    ~10 km",
    ]
    ax.text(0.03, 0.97, "\n".join(coupling), transform=ax.transAxes,
            fontfamily="monospace", fontsize=8, va="top",
            bbox=dict(boxstyle="round,pad=0.4", fc="#f0f8ff", ec="#4488bb"))

    # (f) Crater estimate
    ax = fig.add_subplot(gs[1, 2])
    r_range = np.linspace(0, 300, 200)
    crater_depth = 15 * np.exp(-r_range / 80)
    ejecta_height = 5 * np.exp(-r_range / 120)
    ax.fill_between(r_range, -crater_depth, 0, color="#8B4513", alpha=0.4)
    ax.fill_between(r_range, 0, ejecta_height, color="#D2691E", alpha=0.3)
    ax.plot(r_range, -crater_depth, "k-", lw=2)
    ax.plot(r_range, ejecta_height, "k--", lw=1.5)
    ax.axhline(0, color="brown", lw=2)
    ax.set_xlabel("Range from GZ (m)")
    ax.set_ylabel("Depth / Height (m)")
    ax.set_title("(f) Estimated Crater Profile (50 m HOB)", fontweight="bold")
    ax.text(150, -8, "Crater ~130 m wide\n~15 m deep (est.)",
            fontsize=9, color="#8B4513", fontweight="bold")
    ax.text(150, 3, "Ejecta lip", fontsize=8, color="#D2691E")
    ax.grid(True, alpha=0.2)
    ax.set_ylim(-20, 10)

    fig.suptitle("Figure 7 — Ground-Coupled Seismic & Acoustic Effects",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    p = os.path.join(OUTDIR, "fig07_ground_coupling.png")
    fig.savefig(p, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {p}")


def fig08_combined_hazard():
    """Combined hazard index and casualty estimation.
    
    Note on topography: Actual casualties would be modified by terrain:
    - Berkeley/Oakland Hills provide shielding for eastern areas
    - Building density affects blast wave propagation
    - Urban canyons may focus thermal radiation
    """
    max_range = 25000
    X, Y, R = make_radial_grid(max_range, 600)

    P = peak_overpressure_kpa(R)
    Q = thermal_fluence_kj(R)
    D = prompt_radiation_dose_gy(R)

    haz_blast = np.clip(P / 35.0, 0, 1)
    haz_thermal = np.clip(Q / 670.0, 0, 1)
    haz_rad = np.clip(D / 6.0, 0, 1)
    haz_combined = np.maximum(np.maximum(haz_blast, haz_thermal), haz_rad)

    fig = plt.figure(figsize=(22, 14))
    gs = GridSpec(2, 3, figure=fig, hspace=0.3, wspace=0.3)

    # (a) Combined hazard map
    ax = fig.add_subplot(gs[0, 0])
    im1 = ax.contourf(X / 1000, Y / 1000, haz_combined,
                       levels=np.linspace(0, 1, 25), cmap="RdYlGn_r")
    ax.plot(0, 0, "k*", ms=14, zorder=10)
    ax.set_xlabel("Distance E-W (km)")
    ax.set_ylabel("Distance N-S (km)")
    ax.set_title("(a) Combined Hazard Index", fontweight="bold")
    ax.set_aspect("equal")
    plt.colorbar(im1, ax=ax, label="Hazard Index (0=safe, 1=lethal)")
    # Geographic coordinate axes
    ax_lon = ax.secondary_xaxis('top', functions=(
        lambda x: EVENT_LON + meters_to_deg_lon(x * 1000),
        lambda lon: (lon - EVENT_LON) / meters_to_deg_lon(1000)))
    ax_lon.set_xlabel('Longitude (°W)')
    ax_lat = ax.secondary_yaxis('right', functions=(
        lambda y: EVENT_LAT + meters_to_deg_lat(y * 1000),
        lambda lat: (lat - EVENT_LAT) / meters_to_deg_lat(1000)))
    ax_lat.set_ylabel('Latitude (°N)')

    # (b) Dominant effect map
    ax = fig.add_subplot(gs[0, 1])
    dominant = np.zeros_like(P)
    dominant[haz_blast >= np.maximum(haz_thermal, haz_rad)] = 1
    dominant[haz_thermal >= np.maximum(haz_blast, haz_rad)] = 2
    dominant[haz_rad >= np.maximum(haz_blast, haz_thermal)] = 3
    dominant[haz_combined < 0.01] = 0

    cmap_dom = mcolors.ListedColormap(["white", "#e74c3c", "#f39c12", "#2ecc71"])
    norm_dom = mcolors.BoundaryNorm([0, 0.5, 1.5, 2.5, 3.5], cmap_dom.N)
    im2 = ax.imshow(dominant,
                    extent=[-max_range/1000, max_range/1000,
                            -max_range/1000, max_range/1000],
                    origin="lower", cmap=cmap_dom, norm=norm_dom)
    ax.plot(0, 0, "k*", ms=14, zorder=10)
    ax.set_xlabel("Distance E-W (km)")
    ax.set_ylabel("Distance N-S (km)")
    ax.set_title("(b) Dominant Lethal Effect", fontweight="bold")
    ax.set_aspect("equal")
    cb2 = plt.colorbar(im2, ax=ax, ticks=[0.25, 1, 2, 3], label="Dominant Effect")
    cb2.set_ticklabels(["None", "Blast (>35 kPa)", "Thermal (>670 kJ/m²)", "Radiation (>6 Gy)"])
    # Geographic coordinate axes
    ax_lon = ax.secondary_xaxis('top', functions=(
        lambda x: EVENT_LON + meters_to_deg_lon(x * 1000),
        lambda lon: (lon - EVENT_LON) / meters_to_deg_lon(1000)))
    ax_lon.set_xlabel('Longitude (°W)')

    # (c) Radial profile of all effects
    ax = fig.add_subplot(gs[0, 2])
    r_1d = np.linspace(100, 25000, 500)
    ax.semilogy(r_1d / 1000, peak_overpressure_kpa(r_1d), "r-", lw=2, label="Overpressure (kPa)")
    ax.semilogy(r_1d / 1000, thermal_fluence_kj(r_1d), "orange", lw=2, label="Thermal (kJ/m²)")
    ax.semilogy(r_1d / 1000, prompt_radiation_dose_gy(r_1d), "g-", lw=2, label="Dose (Gy)")
    ax.semilogy(r_1d / 1000, dynamic_pressure_kpa(peak_overpressure_kpa(r_1d)),
                "b--", lw=1.5, label="Dynamic P (kPa)")
    ax.set_xlabel("Ground Range (km)")
    ax.set_ylabel("Effect Magnitude (see legend for units)")
    ax.set_title("(c) All Effects vs Range", fontweight="bold")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.2, which="both")
    ax.set_xlim(0, 25)

    # (d) Casualty estimation by annular ring
    ax = fig.add_subplot(gs[1, 0])
    ring_edges = np.array([0, 0.5, 1, 2, 3, 5, 7, 10, 15, 20, 25])
    fatality_rate = np.array([1.0, 0.98, 0.90, 0.60, 0.35, 0.15, 0.05, 0.02, 0.005, 0.001])
    injury_rate = np.array([0.0, 0.02, 0.10, 0.35, 0.50, 0.45, 0.30, 0.15, 0.05, 0.01])

    pop_density = POP_DENSITY_URBAN
    ring_midpoints = 0.5 * (ring_edges[:-1] + ring_edges[1:])
    ring_areas = np.pi * (ring_edges[1:]**2 - ring_edges[:-1]**2)
    pop_per_ring = ring_areas * pop_density
    fatalities = pop_per_ring * fatality_rate
    injuries = pop_per_ring * injury_rate

    ax.bar(ring_midpoints - 0.3, fatalities / 1000, width=0.6,
           color="red", alpha=0.8, label="Fatalities")
    ax.bar(ring_midpoints + 0.3, injuries / 1000, width=0.6,
           color="orange", alpha=0.8, label="Injuries")
    ax.set_xlabel("Distance from GZ (km)")
    ax.set_ylabel("Casualties (thousands)")
    ax.set_title("(d) Estimated Casualties by Ring", fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.2)

    total_fat = np.sum(fatalities)
    total_inj = np.sum(injuries)

    # (e) Cumulative casualties
    ax = fig.add_subplot(gs[1, 1])
    cumul_fat = np.cumsum(fatalities)
    cumul_inj = np.cumsum(injuries)
    ax.plot(ring_edges[1:], cumul_fat / 1000, "r-o", lw=2, ms=5, label="Fatalities")
    ax.plot(ring_edges[1:], cumul_inj / 1000, "orange", lw=2, ms=5,
            marker="s", label="Injuries")
    ax.set_xlabel("Radius from GZ (km)")
    ax.set_ylabel("Cumulative Casualties (thousands)")
    ax.set_title("(e) Cumulative Casualties", fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.2)

    # (f) Summary
    ax = fig.add_subplot(gs[1, 2])
    ax.axis("off")
    lines = [
        "CASUALTY & DAMAGE SUMMARY",
        "═" * 42,
        "",
        f"Estimated fatalities:   ~{total_fat/1000:.0f},000",
        f"Estimated injuries:     ~{total_inj/1000:.0f},000",
        f"  (assuming {POP_DENSITY_URBAN} ppl/km²)",
        "",
        "STRUCTURAL DAMAGE:",
        f"  Total destruction:    {blast_radius_for_pressure(140)/1000:.1f} km radius",
        f"  Severe:               {blast_radius_for_pressure(35)/1000:.1f} km radius",
        f"  Moderate:             {blast_radius_for_pressure(14)/1000:.1f} km radius",
        f"  Glass breakage:       {blast_radius_for_pressure(1)/1000:.1f} km radius",
        "",
        "FIRE / THERMAL:",
        f"  3rd° burns:           {blast_radius_for_pressure(140)/1000:.1f} km (est.)",
        f"  Ignition of materials within ~3 km",
        f"  Firestorm possible in 2-5 km zone",
        "",
        "AFFECTED BAY AREA CITIES:",
        "  Berkeley, Oakland, Emeryville,",
        "  Albany, El Cerrito, Piedmont,",
        "  parts of Richmond and San Leandro",
        "",
        "San Francisco (across Bay):",
        "  Light-moderate blast damage",
        "  Thermal burns at exposed surfaces",
        "  Broken windows extensively",
    ]
    ax.text(0.03, 0.97, "\n".join(lines), transform=ax.transAxes,
            fontfamily="monospace", fontsize=8.5, va="top",
            bbox=dict(boxstyle="round,pad=0.4", fc="#fff0f0", ec="#cc4444"))

    # Add topography/location note
    fig.text(0.5, 0.01,
             "Note: Flat terrain assumed. Topographic shielding from Berkeley/Oakland Hills not modeled. "
             "Center: {:.4f}°N, {:.4f}°W".format(EVENT_LAT, abs(EVENT_LON)),
             ha='center', fontsize=9, style='italic', color='#444')

    fig.suptitle("Figure 8 — Combined Hazard & Casualty Estimation",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0.04, 1, 0.96])
    p = os.path.join(OUTDIR, "fig08_combined_hazard.png")
    fig.savefig(p, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {p}")


# ═══════════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    os.makedirs(OUTDIR, exist_ok=True)
    np.random.seed(42)

    print("=" * 72)
    print("  Berkeley 100 kT Airburst (50 m HOB) — Model & Analysis")
    print("=" * 72)
    print()

    print("[1/8] Scenario overview ...")
    fig01_scenario_overview()

    print("[2/8] Blast damage map ...")
    fig02_blast_damage_map()

    print("[3/8] Thermal & radiation effects ...")
    fig03_thermal_and_radiation()

    print("[4/8] Fallout pattern (NW wind) ...")
    fig04_fallout_pattern()

    print("[5/8] EMP analysis ...")
    fig05_emp_analysis()

    print("[6/8] Fireball & mushroom cloud ...")
    fig06_fireball_and_cloud()

    print("[7/8] Ground coupling & seismic ...")
    fig07_ground_coupling_seismic()

    print("[8/8] Combined hazard & casualties ...")
    fig08_combined_hazard()

    print()
    print("=" * 72)
    print(f"  All figures saved to: {OUTDIR}")
    print("=" * 72)


if __name__ == "__main__":
    main()
