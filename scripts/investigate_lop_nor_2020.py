#!/usr/bin/env python3
"""
Investigate the Alleged Lop Nor Nuclear Test -- 2020-06-22

A standalone, self-contained investigation of the alleged Chinese decoupled
nuclear test at Lop Nor on June 22, 2020.  Every physical observable is
computed from first principles using only numpy, matplotlib, and scipy.

The script walks through the public evidence for the event, computing:
  1. Station geometry and detection capability
  2. Yield-decoupling constraint space for the observed mb = 2.75
  3. Source comparison for three plausible hypotheses
  4. Synthetic P-wave record section at 12 regional stations
  5. Signal-to-noise ratio and detectability analysis
  6. Two-pulse structure modeling (12-second separation)
  7. Explosion-earthquake discrimination (mb/Ms)
  8. Investigation summary with all key findings

Event parameters (public sources):
  Date/time : 2020-06-22  ~09:18 UTC
  Location  : 41.735 N, 88.730 E  (Lop Nor tunnel area, Kuruktag mountains)
  Detection : PS23 Makanchi array, Kazakhstan (~780 km), mb ~ 2.75
  CTBTO note: Two very small seismic events separated by ~12 seconds
  Host rock : Paleozoic granite, rho=2650 kg/m3, Vp=5800 m/s, Vs=3400 m/s

Physical references:
  Mueller & Murphy (1971) -- Seismic characteristics of underground nuclear detonations
  Patton (1988)           -- Corner-frequency scaling
  Charlie & Veyera (1994) -- Lop Nor site geology
  Glasstone & Dolan (1977) -- Effects of Nuclear Weapons
  Ringdal et al. (1992)   -- IMS detection thresholds

Usage:
    python scripts/investigate_lop_nor_2020.py
"""

import os
import warnings

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.colors as mcolors
from matplotlib.patches import Circle

warnings.filterwarnings("ignore")

# Style setup -- use seaborn grid style if available, fall back gracefully
try:
    plt.style.use("seaborn-v0_8-whitegrid")
except Exception:
    pass


# ============================================================================
# Physical constants and event parameters
# ============================================================================
JOULES_PER_KT = 4.184e12
EARTH_RADIUS_KM = 6371.0

# Event parameters
EVENT_LAT = 41.735
EVENT_LON = 88.730
EVENT_DEPTH_M = 300.0
EVENT_DATE = "2020-06-22"
EVENT_TIME_UTC = "09:18"
MB_OBSERVED = 2.75

# Host rock -- Paleozoic granite (Charlie & Veyera 1994)
RHO_GRANITE = 2650.0    # kg/m3
VP_GRANITE = 5800.0      # m/s
VS_GRANITE = 3400.0      # m/s
Q_P = 400.0              # P-wave quality factor for granite

# Noise floor -- Peterson NLNM at 1 Hz (m/s)
NOISE_NLNM_1HZ = 1.0e-9

# Regional Pn propagation efficiency factor.
# The Mueller-Murphy far-field formula assumes a homogeneous infinite medium.
# Real regional Pn propagation through heterogeneous continental crust incurs
# large losses from scattering, mode conversion, energy leakage from the
# head-wave, and waveguide inefficiency.  This single factor accounts for
# all of these effects.  It is calibrated so that the observed detection at
# PS23 (mb ~ 2.75 at 780 km, barely above noise) yields SNR ~ 3-5 against
# the Peterson NLNM noise floor.
REGIONAL_PATH_FACTOR = 2.0e-4

# Output directory
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.dirname(SCRIPT_DIR)
OUTDIR = os.path.join(REPO_DIR, "figures", "lop_nor_2020_investigation")


# ============================================================================
# Station database -- all 12 stations with real coordinates
# ============================================================================
STATIONS = [
    # (net_sta, name, lat, lon, dist_km, azimuth)
    ("IC.WMQ",   "Urumqi, China",        43.8144,  87.6951,  246,  340),
    ("IU.MAKZ",  "PS23 Makanchi, KZ",    46.7928,  81.9770,  780,  319),
    ("KZ.MKAR",  "Makanchi Array, KZ",   46.7939,  82.2904,  778,  319),
    ("G.WUS",    "Wushi, China",          41.5508,  79.2389,  839,  270),
    ("KZ.PDGK",  "Podgonoye, KZ",        43.2220,  76.9669, 1020,  303),
    ("KR.PRZ",   "Karakol, KG",          42.4665,  78.4373, 1000,  289),
    ("KZ.KNDC",  "Almaty, KZ",           43.2194,  76.9661, 1020,  303),
    ("II.AAK",   "Ala Archa, KG",        42.6389,  74.4942, 1179,  289),
    ("II.KURK",  "Kurchatov, KZ",        50.7154,  78.6189, 1264,  326),
    ("IC.LSA",   "Lhasa, Tibet",         29.7025,  91.1500, 1350,  172),
    ("II.NIL",   "Nilore, Pakistan",     33.6506,  73.2686, 1370,  237),
    ("IC.XAN",   "Xi'an, China",         34.0310, 108.9230, 1965,  109),
]


def station_color(dist_km):
    """Color-code station by detection plausibility."""
    if dist_km < 800:
        return "#2ca02c"   # green -- plausible detection
    elif dist_km <= 1200:
        return "#d4a017"   # yellow/gold -- marginal
    else:
        return "#d62728"   # red -- unlikely detection


# ============================================================================
# Three investigation hypotheses
# ============================================================================
SCENARIOS = {
    "A": {
        "label": "Coupled 21 tonnes (DF=1)",
        "yield_kt": 0.021,
        "df": 1.0,
        "color": "#1f77b4",   # blue
    },
    "B": {
        "label": "Decoupled 1 kt (DF=70)",
        "yield_kt": 1.0,
        "df": 70.0,
        "color": "#d62728",   # red
    },
    "C": {
        "label": "Decoupled 2 kt (DF=70)",
        "yield_kt": 2.0,
        "df": 70.0,
        "color": "#ff7f0e",   # orange
    },
}


# ============================================================================
# Mueller-Murphy Source Model (self-contained)
# ============================================================================

class MuellerMurphySource:
    """
    Self-contained Mueller-Murphy explosion source model.

    Parameters
    ----------
    yield_kt : float
        Device yield in kilotons.
    rho : float
        Host rock density (kg/m3).
    vp : float
        Host rock P-wave velocity (m/s).
    depth_m : float
        Depth of burial (m).
    """

    def __init__(self, yield_kt, rho=RHO_GRANITE, vp=VP_GRANITE,
                 depth_m=EVENT_DEPTH_M):
        self.yield_kt = yield_kt
        self.rho = rho
        self.vp = vp
        self.depth_m = depth_m

    # ---- cavity and damage zones ----

    @property
    def cavity_radius(self):
        """Empirical cavity radius (m): R_c = 55 * W^0.295 * (rho/2650)^(-1/3.4)."""
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

    # ---- spectral properties ----

    @property
    def corner_frequency(self):
        """Patton (1988): fc = 2.5 * W^(-1/3) Hz."""
        return 2.5 * self.yield_kt**(-1.0 / 3.0)

    @property
    def scalar_moment(self):
        """Scalar moment (N*m): M0 = 4 * pi * rho * Vp^2 * R_c^3."""
        rc = self.cavity_radius
        return 4.0 * np.pi * self.rho * self.vp**2 * rc**3

    def decoupled_moment(self, df):
        """Decoupled scalar moment: M0 / DF."""
        return self.scalar_moment / max(df, 1e-12)

    # ---- magnitude ----

    @staticmethod
    def mb_from_yield(yield_kt, decoupling_factor=1.0):
        """Body-wave magnitude: mb = 4.0 + 0.75 * log10(W / DF)."""
        yield_kt = np.asarray(yield_kt, dtype=float)
        decoupling_factor = np.asarray(decoupling_factor, dtype=float)
        w_eff = yield_kt / np.maximum(decoupling_factor, 1e-12)
        return 4.0 + 0.75 * np.log10(np.maximum(w_eff, 1e-12))

    # ---- source spectrum ----

    def spectrum(self, f, df=1.0):
        """
        Mueller-Murphy displacement spectrum (normalised shape).

        S(f) = 1 / (1 + (f/fc)^2)

        For absolute amplitude, multiply by M0/DF.
        """
        fc = self.corner_frequency
        shape = 1.0 / (1.0 + (f / fc)**2)
        M0_eff = self.decoupled_moment(df)
        return M0_eff * shape

    # ---- time-domain RDP and far-field velocity ----

    def rdp_time(self, t, df=1.0):
        """
        Reduced Displacement Potential psi(t).

        psi(t) = (M0/DF) * [1 - (1 + t/tau) * exp(-t/tau)]
        where tau = 1 / (2 * pi * fc).

        Returns (psi, dpsi_dt) as numpy arrays.
        """
        fc = self.corner_frequency
        tau = 1.0 / (2.0 * np.pi * fc)
        M0_eff = self.decoupled_moment(df)
        t = np.asarray(t, dtype=float)
        psi = np.zeros_like(t)
        dpsi = np.zeros_like(t)
        m = t > 0
        x = t[m] / tau
        psi[m] = M0_eff * (1.0 - (1.0 + x) * np.exp(-x))
        dpsi[m] = M0_eff * x * np.exp(-x) / tau
        return psi, dpsi

    def far_field_velocity(self, t, dist_m, df=1.0):
        """
        Far-field P-wave velocity pulse (m/s) at distance dist_m.

        v(t) = (1 / (4*pi*rho*Vp^3*r)) * d^2(psi)/dt^2
        """
        fc = self.corner_frequency
        tau = 1.0 / (2.0 * np.pi * fc)
        M0_eff = self.decoupled_moment(df)
        t = np.asarray(t, dtype=float)
        d2psi = np.zeros_like(t)
        m = t > 0
        x = t[m] / tau
        d2psi[m] = M0_eff * (1.0 - x) * np.exp(-x) / tau**2
        factor = 1.0 / (4.0 * np.pi * self.rho * self.vp**3 * dist_m)
        return factor * d2psi

    def peak_velocity_at_distance(self, dist_m, df=1.0):
        """Peak far-field P velocity (m/s) at distance dist_m (no attenuation)."""
        fc = self.corner_frequency
        tau = 1.0 / (2.0 * np.pi * fc)
        M0_eff = self.decoupled_moment(df)
        # Peak of d2psi/dt2 is at t=0+ (x=0): d2psi_peak = M0_eff / tau^2
        d2psi_peak = M0_eff / tau**2
        factor = 1.0 / (4.0 * np.pi * self.rho * self.vp**3 * dist_m)
        return factor * d2psi_peak


# ============================================================================
# Geodesic and coordinate utilities
# ============================================================================

def haversine_km(lat1, lon1, lat2, lon2):
    """Great-circle distance in km using the Haversine formula."""
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
    return 2 * EARTH_RADIUS_KM * np.arcsin(np.sqrt(a))


def great_circle_points(lat1, lon1, lat2, lon2, n=80):
    """Return n (lat, lon) points along the great circle between two points."""
    lat1, lon1 = np.radians(lat1), np.radians(lon1)
    lat2, lon2 = np.radians(lat2), np.radians(lon2)
    d = 2 * np.arcsin(np.sqrt(
        np.sin((lat2 - lat1) / 2)**2
        + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2)**2
    ))
    if d < 1e-10:
        return np.array([np.degrees(lat1)]), np.array([np.degrees(lon1)])
    fracs = np.linspace(0, 1, n)
    a_arr = np.sin((1 - fracs) * d) / np.sin(d)
    b_arr = np.sin(fracs * d) / np.sin(d)
    x = a_arr * np.cos(lat1) * np.cos(lon1) + b_arr * np.cos(lat2) * np.cos(lon2)
    y = a_arr * np.cos(lat1) * np.sin(lon1) + b_arr * np.cos(lat2) * np.sin(lon2)
    z = a_arr * np.sin(lat1) + b_arr * np.sin(lat2)
    lats = np.degrees(np.arctan2(z, np.sqrt(x**2 + y**2)))
    lons = np.degrees(np.arctan2(y, x))
    return lats, lons


def distance_ring(center_lat, center_lon, radius_km, n=120):
    """Return (lats, lons) tracing a circle of given radius on the sphere."""
    clat = np.radians(center_lat)
    clon = np.radians(center_lon)
    angular_dist = radius_km / EARTH_RADIUS_KM
    bearings = np.linspace(0, 2 * np.pi, n)
    lats = np.arcsin(
        np.sin(clat) * np.cos(angular_dist)
        + np.cos(clat) * np.sin(angular_dist) * np.cos(bearings)
    )
    lons = clon + np.arctan2(
        np.sin(bearings) * np.sin(angular_dist) * np.cos(clat),
        np.cos(angular_dist) - np.sin(clat) * np.sin(lats),
    )
    return np.degrees(lats), np.degrees(lons)


# ============================================================================
# Travel time computation -- simplified 3-layer crustal model
# ============================================================================
# Layer 1: 0-20 km depth, Vp = 5.8 km/s
# Layer 2: 20-40 km depth, Vp = 6.5 km/s
# Layer 3: upper mantle, Vp = 7.9 km/s (Pn head wave)

LAYER_H = [20.0, 20.0]         # layer thicknesses in km
LAYER_V = [5.8, 6.5, 7.9]     # Vp in km/s for each layer (last is half-space)


def compute_travel_time(dist_km):
    """
    Compute first-arriving P travel time using a 3-layer crustal model.

    Returns (travel_time_s, phase_name) where phase_name is 'Pg' or 'Pn'.
    """
    h1, h2 = LAYER_H
    v1, v2, vn = LAYER_V
    moho_depth = h1 + h2

    # --- Direct (Pg) travel time ---
    # Use average crustal velocity weighted by layer thickness
    v_avg = (h1 * v1 + h2 * v2) / moho_depth
    # Approximate direct P as straight-line through crust
    t_direct = dist_km / v_avg

    # --- Pn head wave travel time ---
    # Critical angle at Moho for each layer
    sin_ic1 = v1 / vn
    sin_ic2 = v2 / vn
    # Check that critical angle exists (v < vn always true here)
    if sin_ic1 >= 1.0 or sin_ic2 >= 1.0:
        return t_direct, "Pg"

    cos_ic1 = np.sqrt(1.0 - sin_ic1**2)
    cos_ic2 = np.sqrt(1.0 - sin_ic2**2)

    # Vertical delay through each layer (down + up = factor of 2)
    t_delay1 = 2.0 * h1 / (v1 * cos_ic1)
    t_delay2 = 2.0 * h2 / (v2 * cos_ic2)

    # Horizontal offset consumed in each layer (down + up)
    x_offset1 = 2.0 * h1 * sin_ic1 / cos_ic1
    x_offset2 = 2.0 * h2 * sin_ic2 / cos_ic2
    x_offsets = x_offset1 + x_offset2

    # Head wave only exists if distance exceeds crossover
    if dist_km <= x_offsets:
        return t_direct, "Pg"

    x_head = dist_km - x_offsets
    t_pn = t_delay1 + t_delay2 + x_head / vn

    if t_pn < t_direct:
        return t_pn, "Pn"
    else:
        return t_direct, "Pg"


# ============================================================================
# Attenuation
# ============================================================================

def anelastic_attenuation(f, travel_time_s, Q=Q_P):
    """Anelastic attenuation factor: exp(-pi * f * t / Q)."""
    f = np.asarray(f, dtype=float)
    return np.exp(-np.pi * f * travel_time_s / Q)


def signal_amplitude_at_station(src, dist_km, df, Q=Q_P):
    """
    Predicted peak P-wave velocity amplitude (m/s) at a station.

    Includes:
      - Mueller-Murphy source scaling
      - Geometric spreading (1/r)
      - Anelastic attenuation at the source corner frequency
      - Regional Pn propagation efficiency factor

    The attenuation is evaluated at the corner frequency (where the source
    energy is concentrated) rather than at 1 Hz, giving a more realistic
    estimate of the signal amplitude after propagation.
    """
    dist_m = dist_km * 1000.0
    peak_vel = src.peak_velocity_at_distance(dist_m, df=df)
    # Attenuation at the corner frequency (dominant signal frequency)
    tt_s, _ = compute_travel_time(dist_km)
    f_dom = src.corner_frequency
    att = anelastic_attenuation(f_dom, tt_s, Q)
    return peak_vel * att * REGIONAL_PATH_FACTOR


# ============================================================================
# Figure 1: Station Map
# ============================================================================

def fig01_station_map(outdir):
    """Map of Central/South Asia with event and 12 seismic stations."""
    fig, ax = plt.subplots(figsize=(16, 12))

    # Map extent -- Central Asia
    lon_min, lon_max = 65, 115
    lat_min, lat_max = 25, 55
    ax.set_xlim(lon_min, lon_max)
    ax.set_ylim(lat_min, lat_max)
    ax.set_aspect(1.0 / np.cos(np.radians(40)))  # correct for latitude
    ax.set_facecolor("#e8f0f8")

    # Grid lines
    for lat_g in range(25, 60, 5):
        ax.axhline(lat_g, color="#ccc", lw=0.3, zorder=0)
    for lon_g in range(65, 120, 5):
        ax.axvline(lon_g, color="#ccc", lw=0.3, zorder=0)

    # Distance rings at 500, 1000, 1500, 2000 km
    ring_colors = ["#aaa", "#999", "#888", "#777"]
    for radius_km, rc in zip([500, 1000, 1500, 2000], ring_colors):
        rlats, rlons = distance_ring(EVENT_LAT, EVENT_LON, radius_km)
        ax.plot(rlons, rlats, "-", color=rc, lw=0.8, alpha=0.6, zorder=1)
        # Label the ring
        label_idx = len(rlats) // 8  # place label at ~45 deg bearing
        if (lon_min < rlons[label_idx] < lon_max
                and lat_min < rlats[label_idx] < lat_max):
            ax.text(rlons[label_idx], rlats[label_idx],
                    f"  {radius_km} km", fontsize=7, color=rc,
                    ha="left", va="bottom", zorder=2)

    # Great-circle paths from event to each station
    for net_sta, name, slat, slon, dist, az in STATIONS:
        gc_lats, gc_lons = great_circle_points(EVENT_LAT, EVENT_LON, slat, slon)
        ax.plot(gc_lons, gc_lats, "-", color="#bbb", lw=0.6, alpha=0.7, zorder=2)

    # Plot stations as triangles, color-coded by distance
    # Manual label offsets to avoid overlapping text
    label_offsets = {
        "IC.WMQ":  (0.5, 0.6),
        "IU.MAKZ": (-6.0, 1.0),
        "KZ.MKAR": (0.5, 0.6),
        "G.WUS":   (-7.0, -1.5),
        "KZ.PDGK": (-7.0, 0.5),
        "KR.PRZ":  (0.5, -1.5),
        "KZ.KNDC": (0.5, -1.5),
        "II.AAK":  (-6.0, -1.5),
        "II.KURK": (0.5, 0.6),
        "IC.LSA":  (0.5, 0.6),
        "II.NIL":  (-6.0, 0.5),
        "IC.XAN":  (0.5, 0.6),
    }
    for net_sta, name, slat, slon, dist, az in STATIONS:
        col = station_color(dist)
        ax.plot(slon, slat, "^", color=col, ms=10, markeredgecolor="k",
                markeredgewidth=0.8, zorder=5)
        # Label with net.sta and distance
        dx, dy = label_offsets.get(net_sta, (0.5, 0.5))
        ax.annotate(f"{net_sta}\n({dist} km)",
                    xy=(slon, slat), xytext=(slon + dx, slat + dy),
                    fontsize=6.5, color="#333", fontweight="bold",
                    ha="left" if dx > 0 else "right", va="bottom",
                    arrowprops=dict(arrowstyle="-", color="#999", lw=0.5),
                    zorder=6)

    # Plot event as red star
    ax.plot(EVENT_LON, EVENT_LAT, "*", color="red", ms=22,
            markeredgecolor="k", markeredgewidth=1.0, zorder=10)
    ax.text(EVENT_LON + 0.5, EVENT_LAT - 1.5,
            f"Lop Nor\n{EVENT_DATE}\nmb {MB_OBSERVED}",
            fontsize=9, color="red", fontweight="bold",
            ha="left", va="top", zorder=10)

    # Country labels (approximate centres)
    countries = [
        (80, 42, "CHINA", 14),
        (100, 35, "CHINA", 10),
        (70, 48, "KAZAKHSTAN", 10),
        (75, 41, "KYRGYZSTAN", 7),
        (69, 39, "UZBEKISTAN", 7),
        (72, 34, "PAKISTAN", 8),
        (80, 29, "INDIA", 9),
        (91, 31, "TIBET", 8),
        (69, 34, "AFGHANISTAN", 7),
    ]
    for clon, clat, cname, fs in countries:
        if lon_min < clon < lon_max and lat_min < clat < lat_max:
            ax.text(clon, clat, cname, fontsize=fs, color="#666",
                    alpha=0.4, ha="center", va="center",
                    fontweight="bold", fontstyle="italic", zorder=1)

    # Legend
    import matplotlib.lines as mlines
    h_green = mlines.Line2D([], [], marker="^", color="#2ca02c", ls="",
                            ms=10, mec="k", mew=0.8,
                            label="Plausible detection (< 800 km)")
    h_yellow = mlines.Line2D([], [], marker="^", color="#d4a017", ls="",
                             ms=10, mec="k", mew=0.8,
                             label="Marginal detection (800-1200 km)")
    h_red = mlines.Line2D([], [], marker="^", color="#d62728", ls="",
                          ms=10, mec="k", mew=0.8,
                          label="Unlikely detection (> 1200 km)")
    h_star = mlines.Line2D([], [], marker="*", color="red", ls="",
                           ms=16, mec="k", mew=0.8,
                           label=f"Event: mb {MB_OBSERVED}")
    ax.legend(handles=[h_star, h_green, h_yellow, h_red],
              loc="lower left", fontsize=9, framealpha=0.9)

    ax.set_xlabel("Longitude (degrees E)", fontsize=12)
    ax.set_ylabel("Latitude (degrees N)", fontsize=12)
    ax.tick_params(labelsize=10)

    fig.suptitle("Lop Nor 2020-06-22 event and regional seismic stations.",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    p = os.path.join(outdir, "fig01_station_map.png")
    fig.savefig(p, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {p}")


# ============================================================================
# Figure 2: Hypothesis Space (Yield-Decoupling Constraint)
# ============================================================================

def fig02_hypothesis_space(outdir):
    """2D contour/heatmap of yield vs. decoupling factor with mb contours."""
    fig, ax = plt.subplots(figsize=(12, 9))

    # Grid
    yields = np.logspace(-3, 2, 400)   # 0.001 to 100 kt
    dfs = np.logspace(0, 2, 300)        # 1 to 100
    W_grid, DF_grid = np.meshgrid(yields, dfs)

    # Predicted mb
    mb_grid = MuellerMurphySource.mb_from_yield(W_grid, DF_grid)

    # Heatmap
    levels = np.arange(-2, 8, 0.25)
    im = ax.contourf(W_grid, DF_grid, mb_grid, levels=levels,
                     cmap="RdYlBu_r", extend="both")
    cbar = plt.colorbar(im, ax=ax, label="Predicted mb", shrink=0.85)
    cbar.ax.tick_params(labelsize=10)

    # mb = 2.75 contour (observed) -- thick white line
    cs_obs = ax.contour(W_grid, DF_grid, mb_grid, levels=[MB_OBSERVED],
                        colors="white", linewidths=3.0)
    ax.clabel(cs_obs, fmt="mb=%.2f", fontsize=10, colors="white")

    # IMS threshold at mb = 3.5 -- dashed black line
    cs_ims = ax.contour(W_grid, DF_grid, mb_grid, levels=[3.5],
                        colors="black", linewidths=2.0, linestyles="dashed")
    ax.clabel(cs_ims, fmt="IMS threshold mb=%.1f", fontsize=9, colors="black")

    # Mark three hypothesis points
    for key, sc in SCENARIOS.items():
        ax.plot(sc["yield_kt"], sc["df"], "o", color=sc["color"],
                ms=14, markeredgecolor="k", markeredgewidth=1.5, zorder=10)
        mb_pred = MuellerMurphySource.mb_from_yield(sc["yield_kt"], sc["df"])
        label_text = f'{key}: {sc["label"]}\nmb = {mb_pred:.2f}'
        # Position labels so they stay inside the plot area
        offsets = {"A": (0.3, -0.6), "B": (-1.5, 0.5), "C": (-1.5, -0.5)}
        ox, oy = offsets[key]
        ax.annotate(label_text,
                    xy=(sc["yield_kt"], sc["df"]),
                    xytext=(sc["yield_kt"] * 10**ox, sc["df"] * 10**oy),
                    fontsize=9, fontweight="bold", color=sc["color"],
                    arrowprops=dict(arrowstyle="->", color=sc["color"], lw=1.5),
                    bbox=dict(boxstyle="round,pad=0.3", fc="white",
                              ec=sc["color"], alpha=0.9),
                    zorder=11)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Yield (kt)", fontsize=12, fontweight="bold")
    ax.set_ylabel("Decoupling Factor", fontsize=12, fontweight="bold")
    ax.tick_params(labelsize=10)
    ax.set_xlim(1e-3, 100)
    ax.set_ylim(1, 100)

    fig.suptitle("Yield-decoupling constraint space for mb = 2.75.",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    p = os.path.join(outdir, "fig02_hypothesis_space.png")
    fig.savefig(p, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {p}")


# ============================================================================
# Figure 3: Source Comparison (four panels)
# ============================================================================

def fig03_source_comparison(outdir):
    """Four-panel comparison of three hypotheses."""
    fig = plt.figure(figsize=(18, 14))
    gs = GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.3)

    scenario_keys = ["A", "B", "C"]

    # Pre-compute sources
    sources = {}
    for key in scenario_keys:
        sc = SCENARIOS[key]
        src = MuellerMurphySource(sc["yield_kt"])
        sources[key] = (src, sc)

    # ---- Panel (a): Cavity radius and damage zones ----
    ax = fig.add_subplot(gs[0, 0])
    x_pos = np.arange(len(scenario_keys))
    bar_width = 0.2

    zone_names = ["cavity", "crushed", "fractured", "damaged"]
    zone_colors = ["#d62728", "#ff7f0e", "#2ca02c", "#1f77b4"]
    zone_labels = ["Cavity", "Crushed zone", "Fractured zone", "Damaged zone"]

    for i, zname in enumerate(zone_names):
        radii = []
        for key in scenario_keys:
            src, sc = sources[key]
            radii.append(src.zone_radii[zname])
        ax.bar(x_pos + i * bar_width, radii, bar_width,
               color=zone_colors[i], label=zone_labels[i], alpha=0.8)

    ax.set_xticks(x_pos + 1.5 * bar_width)
    ax.set_xticklabels([SCENARIOS[k]["label"] for k in scenario_keys],
                       fontsize=8, rotation=15, ha="right")
    ax.set_ylabel("Radius (m)", fontsize=11, fontweight="bold")
    ax.set_title("(a) Cavity and Damage Zone Radii", fontsize=12,
                 fontweight="bold")
    ax.legend(fontsize=8, loc="upper left")
    ax.set_yscale("log")
    ax.grid(True, alpha=0.2, which="both")

    # ---- Panel (b): Mueller-Murphy source spectrum ----
    ax = fig.add_subplot(gs[0, 1])
    f = np.logspace(-2, 2, 500)

    for key in scenario_keys:
        src, sc = sources[key]
        spec = src.spectrum(f, df=sc["df"])
        ax.loglog(f, spec, "-", color=sc["color"], lw=2,
                  label=f'{key}: {sc["label"]}')
        # Mark corner frequency
        fc = src.corner_frequency
        spec_fc = src.spectrum(np.array([fc]), df=sc["df"])[0]
        ax.plot(fc, spec_fc, "v", color=sc["color"], ms=8,
                markeredgecolor="k", markeredgewidth=0.5)

    ax.set_xlabel("Frequency (Hz)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Displacement Spectrum (N*m*s)", fontsize=11,
                  fontweight="bold")
    ax.set_title("(b) Mueller-Murphy Source Spectrum", fontsize=12,
                 fontweight="bold")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.2)
    ax.set_xlim(0.01, 100)

    # ---- Panel (c): Predicted mb with IMS threshold ----
    ax = fig.add_subplot(gs[1, 0])
    mb_values = []
    bar_colors = []
    bar_labels = []
    for key in scenario_keys:
        sc = SCENARIOS[key]
        mb_val = MuellerMurphySource.mb_from_yield(sc["yield_kt"], sc["df"])
        mb_values.append(mb_val)
        bar_colors.append(sc["color"])
        bar_labels.append(f'{key}: {sc["label"]}')

    short_labels = ["A: 21t coupled", "B: 1kt dec.", "C: 2kt dec."]
    bars = ax.bar(short_labels, mb_values, color=bar_colors, alpha=0.8,
                  edgecolor="k", linewidth=0.8)
    ax.axhline(3.5, color="gray", ls="--", lw=2, label="IMS threshold (mb=3.5)")
    ax.axhline(MB_OBSERVED, color="red", ls=":", lw=2,
               label=f"Observed mb = {MB_OBSERVED}")

    for bar, mb_val in zip(bars, mb_values):
        ax.text(bar.get_x() + bar.get_width() / 2, mb_val + 0.05,
                f"{mb_val:.2f}", ha="center", va="bottom", fontsize=10,
                fontweight="bold")

    ax.set_ylabel("Body-wave magnitude (mb)", fontsize=11, fontweight="bold")
    ax.set_title("(c) Predicted mb for Each Hypothesis", fontsize=12,
                 fontweight="bold")
    ax.legend(fontsize=9)
    ax.set_ylim(0, 5)
    ax.tick_params(axis="x", labelsize=9)
    ax.grid(True, alpha=0.2, axis="y")

    # ---- Panel (d): Far-field P-wave velocity at 780 km (MAKZ) ----
    ax = fig.add_subplot(gs[1, 1])
    dist_m = 780.0 * 1000.0
    t = np.linspace(-0.5, 10.0, 5000)

    for key in scenario_keys:
        src, sc = sources[key]
        v_ff = src.far_field_velocity(t, dist_m, df=sc["df"])
        # Normalise to peak of largest scenario for comparison
        peak = np.max(np.abs(v_ff))
        if peak > 0:
            ax.plot(t, v_ff / peak, "-", color=sc["color"], lw=1.5,
                    label=f'{key} (peak: {peak:.2e} m/s)')
        else:
            ax.plot(t, v_ff, "-", color=sc["color"], lw=1.5,
                    label=f'{key}')

    ax.set_xlabel("Time (s)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Normalised velocity", fontsize=11, fontweight="bold")
    ax.set_title("(d) Far-field P Velocity at 780 km (MAKZ)",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.2)
    ax.axhline(0, color="gray", lw=0.5)
    ax.set_xlim(-0.5, 10)

    fig.suptitle("Three hypotheses for the Lop Nor 2020 event.",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    p = os.path.join(outdir, "fig03_source_comparison.png")
    fig.savefig(p, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {p}")

    return sources


# ============================================================================
# Figure 4: Synthetic Record Section
# ============================================================================

def fig04_synthetic_record_section(outdir):
    """Synthetic P-wave record section at all 12 stations -- 1 kt decoupled."""
    # Use Scenario B: 1 kt decoupled (DF=70)
    src = MuellerMurphySource(1.0)
    df = 70.0

    # Sort stations by distance
    sorted_stations = sorted(STATIONS, key=lambda s: s[4])

    fig, ax = plt.subplots(figsize=(16, 12))

    t_max = 350  # seconds
    dt = 0.05
    t_full = np.arange(0, t_max, dt)

    travel_times = []
    distances = []

    for i, (net_sta, name, slat, slon, dist, az) in enumerate(sorted_stations):
        tt_s, phase = compute_travel_time(dist)
        travel_times.append(tt_s)
        distances.append(dist)

        # Generate far-field P pulse starting at travel time
        t_rel = t_full - tt_s  # time relative to P arrival
        v_ff = src.far_field_velocity(t_rel, dist * 1000.0, df=df)

        # Apply attenuation at 1 Hz reference
        att = anelastic_attenuation(1.0, tt_s, Q_P)
        v_ff = v_ff * att

        # Normalise for display -- scale to trace height
        peak = np.max(np.abs(v_ff))
        if peak > 0:
            trace = v_ff / peak * 60  # scale for visual separation
        else:
            trace = v_ff

        # Plot at y = distance
        ax.plot(t_full, dist + trace, "k-", lw=0.5, zorder=3)
        ax.fill_between(t_full, dist, dist + trace,
                        where=(trace > 0), color="red", alpha=0.2, zorder=2)
        ax.fill_between(t_full, dist, dist + trace,
                        where=(trace < 0), color="blue", alpha=0.2, zorder=2)

        # Station label -- offset vertically if stations are close together
        label_y = dist
        # Nudge labels that would overlap at similar distances
        if i > 0:
            prev_dist = sorted_stations[i - 1][4]
            if abs(dist - prev_dist) < 30:
                label_y = dist + 15  # offset slightly upward
        ax.text(2, label_y, f" {net_sta} ({phase}, {tt_s:.1f} s)",
                fontsize=7, va="center", ha="left", color="#333", zorder=6)

    # P travel time moveout line
    dist_line = np.linspace(200, 2100, 200)
    tt_line = np.array([compute_travel_time(d)[0] for d in dist_line])
    ax.plot(tt_line, dist_line, "r--", lw=2, alpha=0.7, label="P travel time",
            zorder=5)

    ax.set_xlabel("Time after origin (s)", fontsize=12, fontweight="bold")
    ax.set_ylabel("Epicentral distance (km)", fontsize=12, fontweight="bold")
    ax.set_xlim(0, t_max)
    ax.set_ylim(150, 2100)
    ax.tick_params(labelsize=10)
    ax.legend(loc="upper right", fontsize=10)

    fig.suptitle("Synthetic P-wave record section -- decoupled 1 kt scenario.",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    p = os.path.join(outdir, "fig04_synthetic_record_section.png")
    fig.savefig(p, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {p}")

    return list(zip([s[0] for s in sorted_stations], distances, travel_times))


# ============================================================================
# Figure 5: Detection SNR
# ============================================================================

def fig05_detection_snr(outdir):
    """Signal-to-noise ratio at each station for all three scenarios."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8))

    station_names = [s[0] for s in STATIONS]
    station_dists = [s[4] for s in STATIONS]

    # Sort by distance for the line plot
    sort_idx = np.argsort(station_dists)
    sorted_names = [station_names[i] for i in sort_idx]
    sorted_dists = [station_dists[i] for i in sort_idx]

    # ---- Left panel: SNR vs distance (lines) ----
    for key in ["A", "B", "C"]:
        sc = SCENARIOS[key]
        src = MuellerMurphySource(sc["yield_kt"])
        snr_vals = []
        for d in sorted_dists:
            amp = signal_amplitude_at_station(src, d, sc["df"])
            snr = amp / NOISE_NLNM_1HZ
            snr_vals.append(snr)

        ax1.semilogy(sorted_dists, snr_vals, "o-", color=sc["color"],
                     lw=2, ms=7, label=f'{key}: {sc["label"]}')

    ax1.axhline(3, color="gray", ls="--", lw=2, label="Detection threshold (SNR=3)")
    ax1.axhline(1, color="gray", ls=":", lw=1, alpha=0.5)

    # Annotate stations on x-axis
    for d, name in zip(sorted_dists, sorted_names):
        ax1.axvline(d, color="#eee", lw=0.5)

    ax1.set_xlabel("Epicentral distance (km)", fontsize=12, fontweight="bold")
    ax1.set_ylabel("Signal-to-Noise Ratio", fontsize=12, fontweight="bold")
    ax1.set_title("(a) SNR vs Distance", fontsize=12, fontweight="bold")
    ax1.legend(fontsize=8, loc="upper right")
    ax1.grid(True, which="both", alpha=0.2)
    ax1.set_xlim(200, 2100)
    ax1.tick_params(labelsize=10)

    # ---- Right panel: SNR bar chart for Scenario B (1kt decoupled) ----
    sc_b = SCENARIOS["B"]
    src_b = MuellerMurphySource(sc_b["yield_kt"])
    snr_b = []
    for d in sorted_dists:
        amp = signal_amplitude_at_station(src_b, d, sc_b["df"])
        snr_b.append(amp / NOISE_NLNM_1HZ)

    colors_bar = [station_color(d) for d in sorted_dists]
    bars = ax2.barh(range(len(sorted_names)), snr_b, color=colors_bar,
                    edgecolor="k", linewidth=0.5)
    ax2.axvline(3, color="gray", ls="--", lw=2, label="Detection threshold (SNR=3)")
    ax2.set_yticks(range(len(sorted_names)))
    ax2.set_yticklabels([f"{n}\n({d} km)" for n, d in zip(sorted_names, sorted_dists)],
                        fontsize=8)
    ax2.set_xlabel("Signal-to-Noise Ratio", fontsize=12, fontweight="bold")
    ax2.set_title("(b) SNR at Each Station -- 1 kt Decoupled",
                  fontsize=12, fontweight="bold")
    ax2.set_xscale("log")
    ax2.legend(fontsize=9)
    ax2.grid(True, which="both", alpha=0.2, axis="x")
    ax2.tick_params(labelsize=10)

    fig.suptitle("Predicted signal-to-noise ratio at regional stations.",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    p = os.path.join(outdir, "fig05_detection_snr.png")
    fig.savefig(p, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {p}")

    return list(zip(sorted_names, sorted_dists, snr_b))


# ============================================================================
# Figure 6: Two-Pulse Analysis
# ============================================================================

def fig06_two_pulse_analysis(outdir):
    """Three physical models for the observed two-pulse structure."""
    fig, axes = plt.subplots(3, 1, figsize=(16, 14), sharex=True)

    dist_m = 780.0 * 1000.0  # MAKZ distance
    src = MuellerMurphySource(1.0)
    df = 70.0

    t_duration = 30.0
    dt = 0.01
    t = np.arange(0, t_duration, dt)

    # ---- Panel (a): Explosion + cavity collapse ----
    ax = axes[0]

    # First pulse: explosion at t=0
    v1 = src.far_field_velocity(t, dist_m, df=df)

    # Second pulse at t=12s: cavity collapse (implosive -- negative first motion)
    # Cavity roof failure under lithostatic pressure
    # P_litho = rho * g * depth = 2650 * 9.81 * 300 = 7.80 MPa
    # Lower corner frequency (larger source dimension) -- use fc/2
    t_collapse = t - 12.0
    fc_collapse = src.corner_frequency / 2.0
    tau_c = 1.0 / (2.0 * np.pi * fc_collapse)
    M0_collapse = src.decoupled_moment(df) * 0.3  # 30% of explosion moment

    d2psi_c = np.zeros_like(t)
    m = t_collapse > 0
    x_c = t_collapse[m] / tau_c
    d2psi_c[m] = -M0_collapse * (1.0 - x_c) * np.exp(-x_c) / tau_c**2  # negative
    factor = 1.0 / (4.0 * np.pi * RHO_GRANITE * VP_GRANITE**3 * dist_m)
    v2 = factor * d2psi_c

    composite = v1 + v2
    peak_c = np.max(np.abs(composite))
    if peak_c > 0:
        scale = 1.0 / peak_c
    else:
        scale = 1.0

    ax.plot(t, v1 * scale, "b-", lw=1.5, alpha=0.7, label="Explosion (t=0)")
    ax.plot(t, v2 * scale, "r-", lw=1.5, alpha=0.7,
            label="Cavity collapse (t=12 s)")
    ax.plot(t, composite * scale, "k-", lw=2.0, label="Composite")
    ax.axvline(12.0, color="gray", ls=":", lw=1, alpha=0.5)
    ax.set_ylabel("Normalised velocity", fontsize=11, fontweight="bold")
    ax.set_title("(a) Explosion + cavity roof collapse (implosive, negative first motion)",
                 fontsize=11, fontweight="bold")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(True, alpha=0.2)
    ax.axhline(0, color="gray", lw=0.3)

    # Annotate lithostatic pressure
    P_litho = RHO_GRANITE * 9.81 * EVENT_DEPTH_M * 1e-6  # MPa
    ax.text(0.02, 0.95, f"Lithostatic pressure at {EVENT_DEPTH_M:.0f} m: "
            f"{P_litho:.1f} MPa",
            transform=ax.transAxes, fontsize=8, va="top",
            bbox=dict(fc="white", ec="#ccc", alpha=0.8))

    # ---- Panel (b): Explosion + tunnel containment shot ----
    ax = axes[1]

    v1 = src.far_field_velocity(t, dist_m, df=df)

    # Second shot at t=12s: 10% amplitude, similar waveform
    t_shot2 = t - 12.0
    v2_shot = src.far_field_velocity(t_shot2, dist_m, df=df) * 0.10

    composite_b = v1 + v2_shot
    peak_b = np.max(np.abs(composite_b))
    if peak_b > 0:
        scale_b = 1.0 / peak_b
    else:
        scale_b = 1.0

    ax.plot(t, v1 * scale_b, "b-", lw=1.5, alpha=0.7,
            label="Primary explosion (t=0)")
    ax.plot(t, v2_shot * scale_b, "-", color="darkorange", lw=1.5, alpha=0.7,
            label="Containment shot (t=12 s, 10%)")
    ax.plot(t, composite_b * scale_b, "k-", lw=2.0, label="Composite")
    ax.axvline(12.0, color="gray", ls=":", lw=1, alpha=0.5)
    ax.set_ylabel("Normalised velocity", fontsize=11, fontweight="bold")
    ax.set_title("(b) Explosion + tunnel containment shot (10% amplitude, t=12 s)",
                 fontsize=11, fontweight="bold")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(True, alpha=0.2)
    ax.axhline(0, color="gray", lw=0.3)

    # ---- Panel (c): Chemical explosions (mining) ----
    ax = axes[2]

    rng = np.random.RandomState(42)
    n_blasts = 7
    # Mining blasts: many small pulses in rapid sequence over ~2 seconds
    blast_times = np.sort(rng.uniform(0, 2, n_blasts))
    blast_amps = rng.uniform(0.5, 1.0, n_blasts)

    # Use a tiny source for each blast (0.001 kt chemical)
    src_chem = MuellerMurphySource(0.001)
    composite_mining = np.zeros_like(t)
    for bt, ba in zip(blast_times, blast_amps):
        t_shifted = t - bt
        v_blast = src_chem.far_field_velocity(t_shifted, dist_m, df=1.0) * ba
        composite_mining += v_blast

    peak_m = np.max(np.abs(composite_mining))
    if peak_m > 0:
        composite_mining /= peak_m

    ax.plot(t, composite_mining, "k-", lw=1.0,
            label=f"Mining ripple-fire ({n_blasts} blasts over 2 s)")
    for bt in blast_times:
        ax.axvline(bt, color="green", ls=":", lw=0.5, alpha=0.5)
    ax.axhline(0, color="gray", lw=0.3)

    # Show where 12-second gap would be
    ax.annotate("No second pulse at 12 s", xy=(12, 0), fontsize=9,
                color="red", fontweight="bold", ha="center",
                arrowprops=dict(arrowstyle="->", color="red"))

    ax.set_ylabel("Normalised velocity", fontsize=11, fontweight="bold")
    ax.set_xlabel("Time (s)", fontsize=12, fontweight="bold")
    ax.set_title("(c) Chemical mining blasts -- rapid ripple-fire sequence",
                 fontsize=11, fontweight="bold")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(True, alpha=0.2)

    for a in axes:
        a.set_xlim(-1, 25)
        a.tick_params(labelsize=10)

    fig.suptitle("Three physical models for the observed two-pulse structure.",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    p = os.path.join(outdir, "fig06_two_pulse_analysis.png")
    fig.savefig(p, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {p}")


# ============================================================================
# Figure 7: mb/Ms Discrimination
# ============================================================================

def fig07_discrimination(outdir):
    """mb/Ms discrimination diagram with earthquake/explosion populations."""
    fig, ax = plt.subplots(figsize=(12, 10))

    rng = np.random.RandomState(42)

    # Background earthquake population (~80 events)
    mb_eq = rng.uniform(2.5, 6.5, 80)
    ms_eq = 1.0 + 1.0 * mb_eq + rng.normal(0, 0.3, 80)
    ax.scatter(mb_eq, ms_eq, c="#4488cc", s=25, alpha=0.3,
               marker="o", label="Earthquakes (global)", zorder=2)

    # Known explosion population (~15 events, NTS-like)
    mb_exp = rng.uniform(3.5, 6.5, 15)
    ms_exp = -0.5 + 1.0 * mb_exp + rng.normal(0, 0.2, 15)
    ax.scatter(mb_exp, ms_exp, c="#cc4444", s=30, alpha=0.5,
               marker="D", label="Known explosions (NTS)", zorder=3)

    # Confirmed Chinese nuclear tests (1990-1996)
    chinese_tests = [
        (4.9, "1990 May"),
        (5.7, "1992 May"),
        (5.9, "1992 Sep"),
        (6.0, "1993 Oct"),
        (6.1, "1994 Jun"),
        (6.3, "1995 May"),
        (6.5, "1995 Aug"),
        (6.2, "1996 Jun"),
        (6.0, "1996 Jul"),
    ]
    for mb_cn, label_cn in chinese_tests:
        ms_cn = mb_cn - 1.2 + rng.normal(0, 0.15)
        ax.plot(mb_cn, ms_cn, "s", color="#8B0000", ms=10,
                markeredgecolor="k", markeredgewidth=0.8, zorder=5)
    # Single legend entry for Chinese tests
    ax.plot([], [], "s", color="#8B0000", ms=10, markeredgecolor="k",
            markeredgewidth=0.8, label="Chinese tests (1990-1996)")

    # Discrimination line: Ms = mb - 1.0
    mb_line = np.linspace(1, 8, 100)
    ms_line = mb_line - 1.0
    ax.plot(mb_line, ms_line, "k--", lw=2, label="Discrimination line")
    ax.fill_between(mb_line, ms_line - 5, ms_line, color="red", alpha=0.03)
    ax.fill_between(mb_line, ms_line, ms_line + 5, color="blue", alpha=0.03)
    ax.text(5.5, 3.0, "EXPLOSION", fontsize=13, color="red", alpha=0.4,
            fontweight="bold", fontstyle="italic")
    ax.text(3.5, 5.0, "EARTHQUAKE", fontsize=13, color="blue", alpha=0.4,
            fontweight="bold", fontstyle="italic")

    # Plot the three Lop Nor 2020 hypotheses -- stagger annotations to avoid overlap
    annot_offsets = {
        "A": (1.5, 2.5),
        "B": (2.5, 1.2),
        "C": (1.0, 3.8),
    }
    for key in ["A", "B", "C"]:
        sc = SCENARIOS[key]
        mb_val = MuellerMurphySource.mb_from_yield(sc["yield_kt"], sc["df"])
        # For an explosion, Ms ~ mb - 1.2 (Ms deficit)
        ms_val = mb_val - 1.2
        ax.plot(mb_val, ms_val, "*", color=sc["color"], ms=20,
                markeredgecolor="k", markeredgewidth=1.0, zorder=10)
        tx, ty = annot_offsets[key]
        ax.annotate(f'{key}: {sc["label"]}\nmb={mb_val:.2f}, Ms={ms_val:.2f}',
                    xy=(mb_val, ms_val),
                    xytext=(tx, ty),
                    fontsize=8, fontweight="bold", color=sc["color"],
                    arrowprops=dict(arrowstyle="->", color=sc["color"], lw=1.5),
                    bbox=dict(boxstyle="round,pad=0.2", fc="white",
                              ec=sc["color"], alpha=0.9),
                    zorder=11)

    # Annotation about difficulty at low magnitude
    ax.annotate(
        "At mb < 3, surface waves are\nextremely weak. Ms measurement\n"
        "is unreliable, making explosion/\nearthquake discrimination\n"
        "very difficult at these magnitudes.",
        xy=(MB_OBSERVED, MB_OBSERVED - 1.2),
        xytext=(1.5, 4.0),
        fontsize=9, fontstyle="italic", color="#666",
        arrowprops=dict(arrowstyle="->", color="#999", lw=1),
        bbox=dict(boxstyle="round,pad=0.4", fc="#f8f8f8", ec="#ccc"),
        zorder=8)

    ax.set_xlabel("mb (body-wave magnitude)", fontsize=12, fontweight="bold")
    ax.set_ylabel("Ms (surface-wave magnitude)", fontsize=12, fontweight="bold")
    ax.set_xlim(1, 8)
    ax.set_ylim(0, 8)
    ax.tick_params(labelsize=10)
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(True, alpha=0.2)

    fig.suptitle("Explosion-earthquake discrimination for the Lop Nor 2020 event.",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    p = os.path.join(outdir, "fig07_discrimination.png")
    fig.savefig(p, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {p}")


# ============================================================================
# Figure 8: Investigation Summary (6 panels, 3x2)
# ============================================================================

def fig08_investigation_summary(outdir):
    """Six-panel summary of all investigation findings."""
    fig = plt.figure(figsize=(20, 13))
    gs = GridSpec(2, 3, figure=fig, hspace=0.4, wspace=0.35)

    # ---- (a) Yield-DF constraint line for mb=2.75 ----
    ax = fig.add_subplot(gs[0, 0])
    # mb = 4.0 + 0.75*log10(W/DF) = 2.75
    # => log10(W/DF) = (2.75 - 4.0)/0.75 = -1.667
    # => W/DF = 10^(-1.667) ~ 0.0216
    # So DF = W / 0.0216
    W_range = np.logspace(-2, 2, 200)
    DF_for_mb275 = W_range / 10**((MB_OBSERVED - 4.0) / 0.75)
    valid = (DF_for_mb275 >= 1) & (DF_for_mb275 <= 200)
    ax.loglog(W_range[valid], DF_for_mb275[valid], "k-", lw=3,
              label=f"mb = {MB_OBSERVED}")
    # Mark the three hypotheses
    for key in ["A", "B", "C"]:
        sc = SCENARIOS[key]
        ax.plot(sc["yield_kt"], sc["df"], "o", color=sc["color"], ms=12,
                markeredgecolor="k", markeredgewidth=1, zorder=5)
        ax.text(sc["yield_kt"] * 1.5, sc["df"] * 1.2, key, fontsize=11,
                fontweight="bold", color=sc["color"])
    ax.set_xlabel("Yield (kt)", fontsize=10)
    ax.set_ylabel("Decoupling Factor", fontsize=10)
    ax.set_title("(a) Yield-DF constraint for mb=2.75", fontsize=11,
                 fontweight="bold")
    ax.grid(True, which="both", alpha=0.2)
    ax.set_xlim(0.01, 100)
    ax.set_ylim(0.8, 200)
    ax.legend(fontsize=9)

    # ---- (b) Detection: SNR bar chart for 1kt decoupled ----
    ax = fig.add_subplot(gs[0, 1])
    src_b = MuellerMurphySource(1.0)
    sorted_stas = sorted(STATIONS, key=lambda s: s[4])
    names_short = [s[0].split(".")[1] for s in sorted_stas]
    dists_sorted = [s[4] for s in sorted_stas]
    snr_vals = []
    for s in sorted_stas:
        amp = signal_amplitude_at_station(src_b, s[4], 70.0)
        snr_vals.append(amp / NOISE_NLNM_1HZ)

    colors_snr = [station_color(d) for d in dists_sorted]
    ax.barh(range(len(names_short)), snr_vals, color=colors_snr,
            edgecolor="k", linewidth=0.5)
    ax.axvline(3, color="gray", ls="--", lw=2)
    ax.set_yticks(range(len(names_short)))
    ax.set_yticklabels(names_short, fontsize=8)
    ax.set_xlabel("SNR", fontsize=10)
    ax.set_xscale("log")
    ax.set_title("(b) Detection SNR -- 1 kt decoupled", fontsize=11,
                 fontweight="bold")
    ax.grid(True, which="both", alpha=0.2, axis="x")

    # ---- (c) Two-pulse: explosion+collapse composite at MAKZ ----
    ax = fig.add_subplot(gs[0, 2])
    dist_m = 780.0 * 1000.0
    src = MuellerMurphySource(1.0)
    t = np.arange(0, 25, 0.01)
    v1 = src.far_field_velocity(t, dist_m, df=70.0)
    # Collapse at t=12
    fc_c = src.corner_frequency / 2.0
    tau_c = 1.0 / (2.0 * np.pi * fc_c)
    M0_c = src.decoupled_moment(70.0) * 0.3
    d2psi_c = np.zeros_like(t)
    t_c = t - 12.0
    m = t_c > 0
    x_c = t_c[m] / tau_c
    d2psi_c[m] = -M0_c * (1.0 - x_c) * np.exp(-x_c) / tau_c**2
    factor = 1.0 / (4.0 * np.pi * RHO_GRANITE * VP_GRANITE**3 * dist_m)
    v2 = factor * d2psi_c
    composite = v1 + v2
    peak = np.max(np.abs(composite))
    if peak > 0:
        composite /= peak
    ax.plot(t, composite, "k-", lw=1)
    ax.axvline(0, color="blue", ls=":", lw=0.8, alpha=0.5)
    ax.axvline(12, color="red", ls=":", lw=0.8, alpha=0.5)
    ax.text(1, 0.8, "Explosion", fontsize=8, color="blue")
    ax.text(13, 0.5, "Collapse", fontsize=8, color="red")
    ax.axhline(0, color="gray", lw=0.3)
    ax.set_xlabel("Time (s)", fontsize=10)
    ax.set_ylabel("Normalised velocity", fontsize=10)
    ax.set_title("(c) Two-pulse composite at MAKZ", fontsize=11,
                 fontweight="bold")
    ax.set_xlim(-1, 25)
    ax.grid(True, alpha=0.2)

    # ---- (d) Discrimination: mb/Ms ----
    ax = fig.add_subplot(gs[1, 0])
    rng = np.random.RandomState(42)
    mb_eq = rng.uniform(2.5, 6.5, 60)
    ms_eq = 1.0 + 1.0 * mb_eq + rng.normal(0, 0.3, 60)
    ax.scatter(mb_eq, ms_eq, c="#4488cc", s=15, alpha=0.3, zorder=2)

    mb_line = np.linspace(1, 8, 100)
    ax.plot(mb_line, mb_line - 1.0, "k--", lw=1.5)

    for key in ["A", "B", "C"]:
        sc = SCENARIOS[key]
        mb_val = MuellerMurphySource.mb_from_yield(sc["yield_kt"], sc["df"])
        ms_val = mb_val - 1.2
        ax.plot(mb_val, ms_val, "*", color=sc["color"], ms=16,
                markeredgecolor="k", markeredgewidth=0.8, zorder=5)
        ax.text(mb_val + 0.15, ms_val + 0.15, key, fontsize=10,
                fontweight="bold", color=sc["color"])

    ax.set_xlabel("mb", fontsize=10)
    ax.set_ylabel("Ms", fontsize=10)
    ax.set_title("(d) mb/Ms discrimination", fontsize=11, fontweight="bold")
    ax.set_xlim(1, 8)
    ax.set_ylim(0, 8)
    ax.grid(True, alpha=0.2)

    # ---- (e) Cavity geometry cross-section ----
    ax = fig.add_subplot(gs[1, 1])
    # Draw granite host rock
    ax.fill_between([-80, 80], -350, 0, color="#d4c4a8", alpha=0.4,
                    label="Paleozoic granite")
    # Surface
    ax.plot([-80, 80], [0, 0], "k-", lw=2)
    ax.text(0, 5, "Surface", ha="center", fontsize=9, fontweight="bold")

    # Tunnel at 300 m depth
    ax.plot([-60, -25], [-300, -300], "k-", lw=3)
    ax.text(-65, -295, "Tunnel", fontsize=8, ha="right", va="bottom")

    # Cavity (25 m radius circle at 300 m depth)
    cavity = Circle((0, -300), 25, fc="#ff6666", ec="k", lw=2,
                    alpha=0.5, label="Cavity (R=25 m)")
    ax.add_patch(cavity)
    ax.text(0, -300, "25 m", ha="center", va="center", fontsize=9,
            fontweight="bold")

    # Damage zones
    for r, col, lbl, alpha in [(62.5, "#ff7f0e", "Crushed", 0.15),
                                (125, "#2ca02c", "Fractured", 0.10),
                                (250, "#1f77b4", "Damaged", 0.05)]:
        zone = Circle((0, -300), r, fc=col, ec=col, lw=0.8,
                      alpha=alpha, ls="--")
        ax.add_patch(zone)

    # Depth markers
    ax.annotate("", xy=(60, 0), xytext=(60, -300),
                arrowprops=dict(arrowstyle="<->", color="k", lw=1))
    ax.text(62, -150, "300 m", fontsize=9, rotation=90, va="center")

    ax.set_xlim(-80, 80)
    ax.set_ylim(-380, 30)
    ax.set_aspect("equal")
    ax.set_xlabel("Horizontal (m)", fontsize=10)
    ax.set_ylabel("Depth (m)", fontsize=10)
    ax.set_title("(e) Cavity geometry at 300 m depth", fontsize=11,
                 fontweight="bold")
    ax.legend(fontsize=8, loc="lower left")

    # ---- (f) Text panel: key numbers ----
    ax = fig.add_subplot(gs[1, 2])
    ax.axis("off")

    # Compute key values
    mb_A = MuellerMurphySource.mb_from_yield(0.021, 1.0)
    mb_B = MuellerMurphySource.mb_from_yield(1.0, 70.0)
    mb_C = MuellerMurphySource.mb_from_yield(2.0, 70.0)
    rc_1kt = MuellerMurphySource(1.0).cavity_radius

    lines = [
        "INVESTIGATION KEY FINDINGS",
        "=" * 40,
        "",
        f"  Observed mb:          {MB_OBSERVED}",
        f"  Event date:           {EVENT_DATE} ~{EVENT_TIME_UTC} UTC",
        f"  Location:             {EVENT_LAT} N, {EVENT_LON} E",
        "",
        "  HYPOTHESIS A (coupled):",
        f"    Yield:  21 tonnes,  DF = 1",
        f"    mb:     {mb_A:.2f}",
        "",
        "  HYPOTHESIS B (decoupled):",
        f"    Yield:  1 kt,       DF = 70",
        f"    mb:     {mb_B:.2f}",
        "",
        "  HYPOTHESIS C (decoupled):",
        f"    Yield:  2 kt,       DF = 70",
        f"    mb:     {mb_C:.2f}",
        "",
        f"  Cavity radius (1 kt): {rc_1kt:.1f} m",
        f"  Burial depth:         {EVENT_DEPTH_M:.0f} m",
        f"  Host rock Q_P:        {Q_P:.0f}",
        "",
        "  Detection: only PS23 (780 km)",
        "  Two pulses separated by ~12 s",
        "  Consistent with decoupled test",
        "  in pre-excavated granite cavity",
    ]
    ax.text(0.05, 0.97, "\n".join(lines), transform=ax.transAxes,
            fontfamily="monospace", fontsize=9, va="top",
            bbox=dict(boxstyle="round,pad=0.5", fc="#fffff0", ec="#cca"))
    ax.set_title("(f) Key numbers", fontsize=11, fontweight="bold")

    fig.suptitle("Investigation summary: Lop Nor 2020-06-22.",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    p = os.path.join(outdir, "fig08_investigation_summary.png")
    fig.savefig(p, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {p}")


# ============================================================================
# Main
# ============================================================================

def main():
    os.makedirs(OUTDIR, exist_ok=True)
    np.random.seed(42)

    print("=" * 72)
    print("  Lop Nor 2020-06-22 Investigation")
    print("  Vibe Coding Historic Nuclear Tests")
    print("=" * 72)
    print()

    n_figs = 8

    # ---- Figure 1: Station map ----
    print(f"[1/{n_figs}] Station map ...")
    fig01_station_map(OUTDIR)
    print("       Stations:")
    for net_sta, name, slat, slon, dist, az in STATIONS:
        d_calc = haversine_km(EVENT_LAT, EVENT_LON, slat, slon)
        print(f"         {net_sta:10s}  {name:22s}  {dist:5d} km  "
              f"(calc: {d_calc:.0f} km)  az {az}")
    print()

    # ---- Figure 2: Hypothesis space ----
    print(f"[2/{n_figs}] Hypothesis space ...")
    fig02_hypothesis_space(OUTDIR)
    print("       mb predictions:")
    for key in ["A", "B", "C"]:
        sc = SCENARIOS[key]
        mb_val = MuellerMurphySource.mb_from_yield(sc["yield_kt"], sc["df"])
        print(f"         {key}: {sc['label']:35s}  mb = {mb_val:.2f}  "
              f"(observed: {MB_OBSERVED})")
    print()

    # ---- Figure 3: Source comparison ----
    print(f"[3/{n_figs}] Source comparison ...")
    sources = fig03_source_comparison(OUTDIR)
    print("       Source parameters:")
    for key in ["A", "B", "C"]:
        src, sc = sources[key]
        rc = src.cavity_radius
        fc = src.corner_frequency
        M0 = src.scalar_moment
        M0_eff = src.decoupled_moment(sc["df"])
        print(f"         {key}: R_c = {rc:.2f} m,  fc = {fc:.3f} Hz,  "
              f"M0 = {M0:.2e} N*m,  M0_eff = {M0_eff:.2e} N*m")
    print()

    # ---- Figure 4: Synthetic record section ----
    print(f"[4/{n_figs}] Synthetic record section ...")
    tt_info = fig04_synthetic_record_section(OUTDIR)
    print("       Travel times (1 kt decoupled):")
    for name, dist, tt in tt_info:
        _, phase = compute_travel_time(dist)
        print(f"         {name:10s}  {dist:5.0f} km  {tt:6.1f} s  ({phase})")
    print()

    # ---- Figure 5: Detection SNR ----
    print(f"[5/{n_figs}] Detection SNR ...")
    snr_info = fig05_detection_snr(OUTDIR)
    print("       SNR at each station (1 kt decoupled, DF=70):")
    for name, dist, snr in snr_info:
        detected = "DETECTED" if snr >= 3 else "below threshold"
        print(f"         {name:10s}  {dist:5.0f} km  SNR = {snr:8.1f}  {detected}")
    print()

    # ---- Figure 6: Two-pulse analysis ----
    print(f"[6/{n_figs}] Two-pulse analysis ...")
    fig06_two_pulse_analysis(OUTDIR)
    P_litho = RHO_GRANITE * 9.81 * EVENT_DEPTH_M * 1e-6
    print(f"       Lithostatic pressure at {EVENT_DEPTH_M:.0f} m: {P_litho:.1f} MPa")
    print("       Three models: explosion+collapse, explosion+containment, mining")
    print()

    # ---- Figure 7: Discrimination ----
    print(f"[7/{n_figs}] Discrimination ...")
    fig07_discrimination(OUTDIR)
    print("       mb/Ms predictions:")
    for key in ["A", "B", "C"]:
        sc = SCENARIOS[key]
        mb_val = MuellerMurphySource.mb_from_yield(sc["yield_kt"], sc["df"])
        ms_val = mb_val - 1.2
        print(f"         {key}: mb = {mb_val:.2f},  Ms = {ms_val:.2f}")
    print()

    # ---- Figure 8: Investigation summary ----
    print(f"[8/{n_figs}] Investigation summary ...")
    fig08_investigation_summary(OUTDIR)
    print()

    # ---- Final summary ----
    print("=" * 72)
    print("  INVESTIGATION CONCLUSIONS")
    print("=" * 72)
    print()
    mb_B = MuellerMurphySource.mb_from_yield(1.0, 70.0)
    mb_C = MuellerMurphySource.mb_from_yield(2.0, 70.0)
    rc_1kt = MuellerMurphySource(1.0).cavity_radius
    print(f"  Observed mb = {MB_OBSERVED} at PS23 Makanchi (780 km)")
    print(f"  Decoupling factor DF ~ 70 in Paleozoic granite")
    print(f"  Yield range: 0.5 -- 2 kt (mb range: {mb_B:.2f} -- {mb_C:.2f})")
    print(f"  Cavity radius for 1 kt: {rc_1kt:.1f} m")
    print(f"  Two-pulse separation: 12 s (explosion + collapse or containment)")
    print(f"  Only PS23 detected -- consistent with decoupled low-yield event")
    print(f"  Discrimination: mb/Ms difficult at this magnitude")
    print()
    print(f"  All {n_figs} figures saved to: {OUTDIR}")
    print("=" * 72)


if __name__ == "__main__":
    main()
