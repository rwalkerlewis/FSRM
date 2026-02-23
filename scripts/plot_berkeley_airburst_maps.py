#!/usr/bin/env python3
"""
PyGMT Maps for 100 kT Airburst over Berkeley, California (50 m HOB)

Generates publication-quality maps with topographic relief showing:
  1. Regional overview with Bay Area context
  2. Blast damage zones overlaid on terrain
  3. Thermal radiation contours on terrain
  4. Fallout plume with prevailing NW wind on terrain
  5. Combined effects / infrastructure map

All maps use @earth_relief background for terrain context.

Requirements:
    pip install pygmt numpy pandas

Usage:
    python scripts/plot_berkeley_airburst_maps.py
"""

import os
import tempfile
import numpy as np
import pygmt

# ═══════════════════════════════════════════════════════════════════════════════
# Configuration
# ═══════════════════════════════════════════════════════════════════════════════
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)
OUTDIR = os.path.join(PROJECT_DIR, "figures", "berkeley_100kt")
os.makedirs(OUTDIR, exist_ok=True)

EVENT_LON = -122.2727
EVENT_LAT = 37.8716
YIELD_KT = 100.0
BURST_HEIGHT = 50.0

WIND_DIR_DEG = 315.0
WIND_SPEED = 15.0

JOULES_PER_KT = 4.184e12


# ═══════════════════════════════════════════════════════════════════════════════
# Physics Functions (duplicated for standalone use)
# ═══════════════════════════════════════════════════════════════════════════════

def peak_overpressure_kpa(ground_range_m, yield_kt=YIELD_KT, hob=BURST_HEIGHT):
    """
    Peak static overpressure (kPa) using calibrated Glasstone-Dolan formula.
    
    Reference for 1 kT optimal-height airburst:
        20 psi (140 kPa) at ~450m
        5 psi (35 kPa) at ~900m  
        2 psi (14 kPa) at ~1400m
    Scales by W^(1/3) for other yields.
    """
    r = np.atleast_1d(np.asarray(ground_range_m, dtype=float))
    R = np.sqrt(r**2 + hob**2)
    
    # Scale factor for yield (cube root scaling)
    W_scale = yield_kt ** (1.0 / 3.0)
    
    # Scaled distance (meters per kt^1/3)
    Z = np.maximum(R / W_scale, 1.0)
    
    # Empirical fit calibrated to Glasstone-Dolan data for airburst
    # Overpressure in kPa as function of scaled distance
    P_kpa = np.where(
        Z < 50,
        # Very close range - extremely high overpressure
        5000.0 * (50.0 / Z) ** 2.5,
        np.where(
            Z < 200,
            # Close range (exponential decay)
            1000.0 * np.exp(-0.015 * (Z - 50)),
            np.where(
                Z < 600,
                # Mid range (power law)
                200.0 * (200.0 / Z) ** 2.0,
                np.where(
                    Z < 2500,
                    # Far range
                    8.0 * (600.0 / Z) ** 1.8,
                    # Very far range
                    0.3 * (2500.0 / Z) ** 1.5
                )
            )
        )
    )
    
    # Mach stem enhancement for low burst angles (increases overpressure)
    angle = np.arctan2(hob, np.maximum(r, 1.0))
    mach_factor = np.where(angle < np.radians(40), 1.8, 1.0)
    
    return np.squeeze(P_kpa * mach_factor)


def blast_radius_for_pressure(target_kpa):
    r = np.linspace(1, 60000, 60000)
    p = peak_overpressure_kpa(r)
    idx = np.argmin(np.abs(p - target_kpa))
    return r[idx]


def thermal_fluence_kj(ground_range_m, yield_kt=YIELD_KT, hob=BURST_HEIGHT):
    r = np.atleast_1d(np.asarray(ground_range_m, dtype=float))
    R = np.sqrt(r**2 + hob**2)
    E_thermal = 0.35 * yield_kt * JOULES_PER_KT
    Q = E_thermal / (4 * np.pi * R**2)
    tau = np.exp(-R / 18000.0)
    return np.squeeze(Q * tau / 1000.0)


def prompt_radiation_gy(ground_range_m, yield_kt=YIELD_KT, hob=BURST_HEIGHT):
    r = np.atleast_1d(np.asarray(ground_range_m, dtype=float))
    R = np.sqrt(r**2 + hob**2)
    D = 10.0 * yield_kt * (1000.0 / R)**2 * np.exp(-R / 2500.0)
    return np.squeeze(D)


def meters_to_deg_lat(m):
    return m / 111320.0

def meters_to_deg_lon(m, lat=EVENT_LAT):
    return m / (111320.0 * np.cos(np.radians(lat)))


def radius_to_circle_coords(radius_m, n=361, lon0=EVENT_LON, lat0=EVENT_LAT):
    """Convert radius (meters) to lon/lat circle coordinates."""
    angles = np.linspace(0, 2 * np.pi, n)
    lons = lon0 + meters_to_deg_lon(radius_m * np.cos(angles))
    lats = lat0 + meters_to_deg_lat(radius_m * np.sin(angles))
    return lons, lats


# ═══════════════════════════════════════════════════════════════════════════════
# Map 1: Regional Overview
# ═══════════════════════════════════════════════════════════════════════════════

def map01_regional_overview():
    """Bay Area regional overview with ground zero."""
    region = [-123.0, -121.5, 37.2, 38.2]
    projection = "M18c"

    fig = pygmt.Figure()

    fig.grdimage(
        grid="@earth_relief_03s",
        region=region,
        projection=projection,
        shading=True,
        cmap="geo",
    )

    fig.basemap(
        region=region,
        projection=projection,
        frame=["WSne+t100 kT Airburst over Berkeley, CA — Regional Overview",
               "xa0.5f0.25", "ya0.25f0.125"],
    )

    fig.coast(
        region=region,
        projection=projection,
        resolution="f",
        water="lightblue",
        borders="2/0.3p,gray50",
        shorelines="0.3p,gray30",
    )

    fig.plot(x=[EVENT_LON], y=[EVENT_LAT],
             style="a0.6c", fill="red", pen="1p,black")
    fig.text(x=EVENT_LON, y=EVENT_LAT - 0.04,
             text="Ground Zero (Berkeley)", font="8p,Helvetica-Bold,red",
             justify="TC")

    cities = [
        (-122.4194, 37.7749, "San Francisco"),
        (-122.2712, 37.8044, "Oakland"),
        (-122.3477, 37.5585, "San Mateo"),
        (-122.0322, 37.3230, "San Jose"),
        (-122.0264, 37.5630, "Fremont"),
        (-122.3106, 37.9477, "Richmond"),
        (-122.0574, 37.8853, "Walnut Creek"),
        (-121.9358, 37.6879, "Pleasanton"),
        (-122.0819, 37.3861, "Milpitas"),
    ]
    for clon, clat, cname in cities:
        fig.plot(x=[clon], y=[clat], style="c0.12c", fill="black", pen="0.3p,black")
        fig.text(x=clon + 0.02, y=clat + 0.02,
                 text=cname, font="6p,Helvetica,black", justify="BL")

    damage_levels = [
        (140.0, "red",    "Total destruction (>140 kPa)"),
        (35.0,  "red",    "Severe damage (>35 kPa)"),
        (14.0,  "orange", "Most buildings destroyed (>14 kPa)"),
        (7.0,   "yellow", "Moderate damage (>7 kPa)"),
        (3.5,   "#aaff00", "Light damage (>3.5 kPa)"),
        (1.0,   "cyan",   "Glass breakage (>1 kPa)"),
    ]
    for thresh, color, label in damage_levels:
        radius = blast_radius_for_pressure(thresh)
        lons, lats = radius_to_circle_coords(radius)
        fig.plot(x=lons, y=lats, pen=f"1.5p,{color}", label=label)

    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        f.write("L 8p,Helvetica-Bold,black L Blast Overpressure Zones:\n")
        for thresh, color, label in damage_levels:
            radius = blast_radius_for_pressure(thresh)
            f.write(f"S 0.2c - 0.4c - 1.5p,{color} 0.8c {label} (r={radius/1000:.1f} km)\n")
        legend_file = f.name

    fig.legend(spec=legend_file, position="JBL+w8c+o0.3c/0.3c",
               box="+p0.5p+gwhite@20+s")
    os.unlink(legend_file)

    fig.basemap(region=region, projection=projection,
                map_scale="g-121.7/37.3+c37.8+w20k+f+l")

    # Add terrain note
    fig.text(x=-122.25, y=37.25, text="Topography: SRTM 3-arcsec DEM. Hills may provide shielding.",
             font="7p,Helvetica-Oblique,gray40", justify="ML")

    with fig.inset(position="jTR+w4c+o0.3c", box="+p0.5p+gwhite"):
        fig.coast(region=[-130, -110, 30, 45], projection="M?",
                  resolution="l", land="gray80", water="lightblue",
                  shorelines="faint", borders="1/0.3p,gray50")
        fig.plot(x=[EVENT_LON], y=[EVENT_LAT], style="a0.3c",
                 fill="red", pen="0.3p,black")

    output = os.path.join(OUTDIR, "map01_regional_overview.png")
    fig.savefig(output, dpi=300)
    print(f"  Saved {output}")


# ═══════════════════════════════════════════════════════════════════════════════
# Map 2: Blast Damage Zones on Terrain
# ═══════════════════════════════════════════════════════════════════════════════

def map02_blast_damage():
    """Close-up blast damage zones overlaid on terrain."""
    half = 0.16
    region = [EVENT_LON - half, EVENT_LON + half,
              EVENT_LAT - half * 0.8, EVENT_LAT + half * 0.8]
    projection = "M18c"

    fig = pygmt.Figure()

    fig.grdimage(grid="@earth_relief_03s", region=region,
                 projection=projection, shading=True, cmap="geo")

    fig.basemap(region=region, projection=projection,
                frame=["WSne+tBlast Damage Zones — 100 kT at Berkeley (50 m HOB)",
                       "xa0.05f0.025", "ya0.05f0.025"])

    fig.coast(region=region, projection=projection, resolution="f",
              water="lightblue", shorelines="0.3p,gray30")

    damage_zones = [
        (140.0, "darkred",   "3p", "Total destruction (>140 kPa)"),
        (35.0,  "red",       "2.5p", "Reinf. concrete destroyed (>35 kPa)"),
        (14.0,  "orange",    "2p", "Most buildings destroyed (>14 kPa)"),
        (7.0,   "yellow",    "2p", "Moderate damage (>7 kPa)"),
        (3.5,   "green",     "1.5p", "Light damage (>3.5 kPa)"),
        (1.0,   "cyan",      "1p", "Glass breakage (>1 kPa)"),
    ]

    for thresh, color, pen_w, label in damage_zones:
        radius = blast_radius_for_pressure(thresh)
        lons, lats = radius_to_circle_coords(radius)
        fig.plot(x=lons, y=lats, pen=f"{pen_w},{color}")

    fig.plot(x=[EVENT_LON], y=[EVENT_LAT],
             style="a0.5c", fill="red", pen="1p,black")
    fig.text(x=EVENT_LON, y=EVENT_LAT + 0.008,
             text="GZ", font="8p,Helvetica-Bold,red", justify="BC")

    local_places = [
        (-122.260, 37.870, "UC Berkeley"),
        (-122.286, 37.869, "Berkeley\nMarina"),
        (-122.273, 37.858, "Downtown\nBerkeley"),
        (-122.250, 37.832, "Rockridge"),
        (-122.200, 37.876, "Tilden\nPark"),
        (-122.257, 37.843, "Temescal"),
        (-122.303, 37.832, "Emeryville"),
        (-122.270, 37.805, "Downtown\nOakland"),
        (-122.227, 37.893, "Kensington"),
        (-122.310, 37.890, "Albany"),
        (-122.350, 37.900, "El Cerrito"),
        (-122.345, 37.867, "Golden Gate\nFields"),
    ]
    for plon, plat, pname in local_places:
        if region[0] < plon < region[1] and region[2] < plat < region[3]:
            fig.plot(x=[plon], y=[plat], style="c0.08c", fill="white",
                     pen="0.3p,black")
            fig.text(x=plon + 0.003, y=plat + 0.003, text=pname,
                     font="5p,Helvetica,black", justify="BL")

    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        f.write("L 7p,Helvetica-Bold,black L Peak Overpressure Zones:\n")
        for thresh, color, pen_w, label in damage_zones:
            radius = blast_radius_for_pressure(thresh)
            f.write(f"S 0.2c - 0.4c - {pen_w},{color} 0.8c {label} (r={radius/1000:.1f} km)\n")
        f.write("L 6p,Helvetica-Oblique,gray40 L Terrain: hills may shield/focus blast\n")
        legend_file = f.name
    fig.legend(spec=legend_file, position="JBL+w8c+o0.3c/0.3c",
               box="+p0.5p+gwhite@20+s")
    os.unlink(legend_file)

    fig.basemap(region=region, projection=projection,
                map_scale=f"g{EVENT_LON+0.12}/{EVENT_LAT-0.1}+c{EVENT_LAT}+w5k+f+l")

    output = os.path.join(OUTDIR, "map02_blast_damage.png")
    fig.savefig(output, dpi=300)
    print(f"  Saved {output}")


# ═══════════════════════════════════════════════════════════════════════════════
# Map 3: Thermal and Radiation Zones
# ═══════════════════════════════════════════════════════════════════════════════

def map03_thermal_radiation():
    """Thermal and prompt radiation contours on terrain."""
    half = 0.15
    region = [EVENT_LON - half, EVENT_LON + half,
              EVENT_LAT - half * 0.85, EVENT_LAT + half * 0.85]
    projection = "M18c"

    fig = pygmt.Figure()
    fig.grdimage(grid="@earth_relief_03s", region=region,
                 projection=projection, shading=True, cmap="geo")
    fig.basemap(region=region, projection=projection,
                frame=["WSne+tThermal & Radiation Zones — 100 kT at Berkeley",
                       "xa0.05f0.025", "ya0.05f0.025"])
    fig.coast(region=region, projection=projection, resolution="f",
              water="lightblue", shorelines="0.3p,gray30")

    thermal_levels = [
        (670, "darkred",  "2.5p,darkred,-", "3rd° burns (>670 kJ/m²)"),
        (335, "red",      "2p,red,-",       "2nd° burns (>335 kJ/m²)"),
        (125, "orange",   "1.5p,orange,-",  "1st° burns (>125 kJ/m²)"),
        (50,  "yellow",   "1p,yellow,-",    "Pain threshold (~50 kJ/m²)"),
    ]

    for thresh_kj, color, pen, label in thermal_levels:
        r_test = np.linspace(1, 30000, 30000)
        Q = thermal_fluence_kj(r_test)
        idx = np.argmin(np.abs(Q - thresh_kj))
        radius = r_test[idx]
        lons, lats = radius_to_circle_coords(radius)
        fig.plot(x=lons, y=lats, pen=pen)

    rad_levels = [
        (10.0, "purple",  "2p,purple",  "Fatal (>10 Gy)"),
        (4.5,  "magenta", "1.5p,magenta", "LD50 (>4.5 Gy)"),
        (2.0,  "blue",    "1.5p,blue",  "Sickness (>2 Gy)"),
        (0.5,  "cyan",    "1p,cyan",    "Mild (>0.5 Gy)"),
    ]

    for thresh_gy, color, pen, label in rad_levels:
        r_test = np.linspace(1, 20000, 20000)
        D = prompt_radiation_gy(r_test)
        idx = np.argmin(np.abs(D - thresh_gy))
        radius = r_test[idx]
        lons, lats = radius_to_circle_coords(radius)
        fig.plot(x=lons, y=lats, pen=pen)

    fig.plot(x=[EVENT_LON], y=[EVENT_LAT],
             style="a0.5c", fill="red", pen="1p,black")

    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        f.write("L 7p,Helvetica-Bold,black L THERMAL FLUENCE:\n")
        for thresh_kj, color, pen, label in thermal_levels:
            r_test = np.linspace(1, 30000, 30000)
            Q = thermal_fluence_kj(r_test)
            idx = np.argmin(np.abs(Q - thresh_kj))
            r_km = r_test[idx] / 1000
            f.write(f"S 0.2c - 0.4c - {pen} 0.8c {label} (r={r_km:.1f} km)\n")
        f.write("L 7p,Helvetica-Bold,black L ABSORBED DOSE:\n")
        for thresh_gy, color, pen, label in rad_levels:
            r_test = np.linspace(1, 20000, 20000)
            D = prompt_radiation_gy(r_test)
            idx = np.argmin(np.abs(D - thresh_gy))
            r_km = r_test[idx] / 1000
            f.write(f"S 0.2c - 0.4c - {pen} 0.8c {label} (r={r_km:.1f} km)\n")
        f.write("L 6p,Helvetica-Oblique,gray40 L Hills block line-of-sight thermal\n")
        legend_file = f.name

    fig.legend(spec=legend_file, position="JBL+w8c+o0.3c/0.3c",
               box="+p0.5p+gwhite@20+s")
    os.unlink(legend_file)

    fig.basemap(region=region, projection=projection,
                map_scale=f"g{EVENT_LON+0.11}/{EVENT_LAT-0.1}+c{EVENT_LAT}+w5k+f+l")

    output = os.path.join(OUTDIR, "map03_thermal_radiation.png")
    fig.savefig(output, dpi=300)
    print(f"  Saved {output}")


# ═══════════════════════════════════════════════════════════════════════════════
# Map 4: Fallout Plume
# ═══════════════════════════════════════════════════════════════════════════════

def map04_fallout_plume():
    """Fallout plume under prevailing NW wind on relief."""
    region = [-122.8, -121.6, 37.2, 38.1]
    projection = "M18c"

    fig = pygmt.Figure()
    fig.grdimage(grid="@earth_relief_03s", region=region,
                 projection=projection, shading=True, cmap="geo")
    fig.basemap(region=region, projection=projection,
                frame=["WSne+tFallout Plume — 100 kT, Prevailing NW Wind",
                       "xa0.25f0.125", "ya0.25f0.125"])
    fig.coast(region=region, projection=projection, resolution="f",
              water="lightblue", shorelines="0.3p,gray30",
              borders="2/0.3p,gray50")

    n_contours = 6
    activity_levels = [0.5, 0.2, 0.1, 0.05, 0.01, 0.005]
    contour_colors = ["darkred", "red", "orange", "yellow", "green", "cyan"]
    contour_labels = [
        "Very high (>50% max)",
        "High (>20% max)",
        "Moderate (>10% max)",
        "Light (>5% max)",
        "Low-level (>1% max)",
        "Trace (>0.5% max)",
    ]

    theta = np.radians(WIND_DIR_DEG)
    dx_w = np.cos(theta)
    dy_w = np.sin(theta)

    for i, (level, color, label) in enumerate(zip(activity_levels,
                                                   contour_colors,
                                                   contour_labels)):
        scale = -2 * np.log(level)
        along_dist = WIND_SPEED * (10000 / 2.0) + 3000 * np.sqrt(scale)
        cross_dist = 2000 * np.sqrt(scale)

        angles = np.linspace(0, 2 * np.pi, 181)
        ellipse_x = along_dist * np.cos(angles)
        ellipse_y = cross_dist * np.sin(angles)

        rot_x = ellipse_x * dx_w - ellipse_y * dy_w
        rot_y = ellipse_x * dy_w + ellipse_y * dx_w

        center_offset_x = along_dist * 0.6 * dx_w
        center_offset_y = along_dist * 0.6 * dy_w

        lons = EVENT_LON + meters_to_deg_lon(rot_x + center_offset_x)
        lats = EVENT_LAT + meters_to_deg_lat(rot_y + center_offset_y)

        fig.plot(x=lons, y=lats, pen=f"1.5p,{color}")

    fig.plot(x=[EVENT_LON], y=[EVENT_LAT],
             style="a0.5c", fill="red", pen="1p,black")

    arrow_len = 0.08
    arrow_end_lon = EVENT_LON + arrow_len * dx_w / np.cos(np.radians(EVENT_LAT))
    arrow_end_lat = EVENT_LAT + arrow_len * dy_w
    fig.plot(x=[EVENT_LON, arrow_end_lon],
             y=[EVENT_LAT, arrow_end_lat],
             pen="3p,red", style="v0.3c+e+a40")
    fig.text(x=arrow_end_lon + 0.02, y=arrow_end_lat + 0.02,
             text="Wind (NW)", font="8p,Helvetica-Bold,red", justify="BL")

    cities = [
        (-122.4194, 37.7749, "San Francisco"),
        (-122.2712, 37.8044, "Oakland"),
        (-122.0264, 37.5630, "Fremont"),
        (-122.0322, 37.3230, "San Jose"),
        (-122.0574, 37.8853, "Walnut Creek"),
        (-122.3106, 37.9477, "Richmond"),
        (-121.9358, 37.6879, "Pleasanton"),
        (-122.0819, 37.3861, "Milpitas"),
        (-122.1600, 37.4450, "Santa Clara"),
        (-122.0580, 37.5510, "Union City"),
        (-122.2200, 37.5000, "Hayward"),
    ]
    for clon, clat, cname in cities:
        fig.plot(x=[clon], y=[clat], style="c0.1c", fill="white", pen="0.3p,black")
        fig.text(x=clon + 0.015, y=clat + 0.015,
                 text=cname, font="5p,Helvetica,black", justify="BL")

    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        f.write("L 7p,Helvetica-Bold,black L Normalized Fallout Activity:\n")
        for color, label in zip(contour_colors, contour_labels):
            f.write(f"S 0.2c - 0.4c - 1.5p,{color} 0.8c {label}\n")
        f.write("L 6p,Helvetica-Oblique,gray40 L Terrain affects deposition patterns\n")
        legend_file = f.name
    fig.legend(spec=legend_file, position="JBR+w6.5c+o0.3c/0.3c",
               box="+p0.5p+gwhite@20+s")
    os.unlink(legend_file)

    fig.basemap(region=region, projection=projection,
                map_scale=f"g-121.8/37.3+c37.7+w20k+f+l")

    output = os.path.join(OUTDIR, "map04_fallout_plume.png")
    fig.savefig(output, dpi=300)
    print(f"  Saved {output}")


# ═══════════════════════════════════════════════════════════════════════════════
# Map 5: Combined Effects with Infrastructure
# ═══════════════════════════════════════════════════════════════════════════════

def map05_combined_effects():
    """Combined effects map with key infrastructure."""
    half = 0.22
    region = [EVENT_LON - half, EVENT_LON + half,
              EVENT_LAT - half * 0.85, EVENT_LAT + half * 0.85]
    projection = "M18c"

    fig = pygmt.Figure()
    fig.grdimage(grid="@earth_relief_03s", region=region,
                 projection=projection, shading=True, cmap="geo")
    fig.basemap(region=region, projection=projection,
                frame=["WSne+tCombined Effects — 100 kT at Berkeley (50 m HOB)",
                       "xa0.1f0.05", "ya0.1f0.05"])
    fig.coast(region=region, projection=projection, resolution="f",
              water="lightblue", shorelines="0.3p,gray30")

    zones = [
        (140.0, "darkred@50",  "2p,darkred"),
        (35.0,  "red@50",      "2p,red"),
        (14.0,  "orange@60",   "1.5p,orange"),
        (7.0,   "yellow@70",   "1.5p,yellow"),
        (3.5,   "green@70",    "1p,green"),
    ]

    for thresh, fill_color, pen in zones:
        radius = blast_radius_for_pressure(thresh)
        lons, lats = radius_to_circle_coords(radius)
        fig.plot(x=lons, y=lats, pen=pen)

    fig.plot(x=[EVENT_LON], y=[EVENT_LAT],
             style="a0.6c", fill="red", pen="1p,black")

    infrastructure = [
        (-122.394, 37.795, "s0.2c", "blue", "Bay Bridge"),
        (-122.478, 37.819, "s0.2c", "blue", "Golden Gate Br."),
        (-122.353, 37.945, "s0.2c", "blue", "Richmond Br."),
        (-122.271, 37.876, "t0.2c", "darkgreen", "UC Berkeley"),
        (-122.258, 37.875, "t0.2c", "darkgreen", "Lawrence Lab"),
        (-122.257, 37.870, "d0.15c", "purple", "Berkeley BART"),
        (-122.267, 37.870, "d0.15c", "purple", "Downtown Berk BART"),
        (-122.268, 37.853, "d0.15c", "purple", "Ashby BART"),
        (-122.274, 37.803, "d0.15c", "purple", "MacArthur BART"),
        (-122.252, 37.828, "d0.15c", "purple", "Rockridge BART"),
        (-122.300, 37.849, "c0.15c", "brown", "Port of Oakland"),
        (-122.375, 37.619, "c0.15c", "brown", "SFO Airport"),
        (-122.206, 37.722, "c0.15c", "brown", "OAK Airport"),
    ]

    for ilon, ilat, style, color, name in infrastructure:
        if region[0] < ilon < region[1] and region[2] < ilat < region[3]:
            fig.plot(x=[ilon], y=[ilat], style=style, fill=color,
                     pen="0.3p,black")
            fig.text(x=ilon + 0.005, y=ilat + 0.005, text=name,
                     font="5p,Helvetica,black", justify="BL")

    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        f.write("L 7p,Helvetica-Bold,black L BLAST OVERPRESSURE ZONES:\n")
        labels = ["Total (>140 kPa)", "Severe (>35 kPa)",
                  "Destroyed (>14 kPa)", "Moderate (>7 kPa)", "Light (>3.5 kPa)"]
        pens = ["2p,darkred", "2p,red", "1.5p,orange", "1.5p,yellow", "1p,green"]
        for lbl, pen in zip(labels, pens):
            r_km = blast_radius_for_pressure(
                float(lbl.split(">")[1].split(" ")[0])) / 1000
            f.write(f"S 0.2c - 0.4c - {pen} 0.8c {lbl} (r={r_km:.1f} km)\n")
        f.write("L 7p,Helvetica-Bold,black L INFRASTRUCTURE:\n")
        f.write("S 0.2c t 0.15c darkgreen 0.3p,black 0.8c University\n")
        f.write("S 0.2c d 0.15c purple 0.3p,black 0.8c BART Station\n")
        f.write("S 0.2c s 0.15c blue 0.3p,black 0.8c Bridge\n")
        f.write("S 0.2c c 0.15c brown 0.3p,black 0.8c Port/Airport\n")
        f.write("L 6p,Helvetica-Oblique,gray40 L Terrain shielding not shown\n")
        legend_file = f.name

    fig.legend(spec=legend_file, position="JBL+w8c+o0.3c/0.3c",
               box="+p0.5p+gwhite@20+s")
    os.unlink(legend_file)

    fig.basemap(region=region, projection=projection,
                map_scale=f"g{EVENT_LON+0.17}/{EVENT_LAT-0.15}+c{EVENT_LAT}+w10k+f+l")

    output = os.path.join(OUTDIR, "map05_combined_effects.png")
    fig.savefig(output, dpi=300)
    print(f"  Saved {output}")


# ═══════════════════════════════════════════════════════════════════════════════
# Map 6: Prompt Radiation on Terrain
# ═══════════════════════════════════════════════════════════════════════════════

def map06_radiation_zones():
    """Prompt radiation dose contours on terrain."""
    half = 0.1
    region = [EVENT_LON - half, EVENT_LON + half,
              EVENT_LAT - half * 0.85, EVENT_LAT + half * 0.85]
    projection = "M18c"

    fig = pygmt.Figure()
    fig.grdimage(grid="@earth_relief_03s", region=region,
                 projection=projection, shading=True, cmap="geo")
    fig.basemap(region=region, projection=projection,
                frame=["WSne+tPrompt Radiation & Thermal Zones — 100 kT",
                       "xa0.025f0.0125", "ya0.025f0.0125"])
    fig.coast(region=region, projection=projection, resolution="f",
              water="lightblue", shorelines="0.3p,gray30")

    rad_contours = [
        (10.0, "2.5p,purple",  "Fatal (>10 Gy)"),
        (4.5,  "2p,magenta",   "LD50 (>4.5 Gy)"),
        (2.0,  "1.5p,blue",    "Rad. sickness (>2 Gy)"),
        (0.5,  "1p,cyan",      "Mild symptoms (>0.5 Gy)"),
        (0.01, "0.5p,green",   "Detectable (>10 mGy)"),
    ]

    for thresh, pen, label in rad_contours:
        r_test = np.linspace(1, 20000, 20000)
        D = prompt_radiation_gy(r_test)
        idx = np.argmin(np.abs(D - thresh))
        radius = r_test[idx]
        lons, lats = radius_to_circle_coords(radius)
        fig.plot(x=lons, y=lats, pen=pen)

    thermal_contours = [
        (670, "2.5p,red,-",    "3rd° burns (>670 kJ/m²)"),
        (335, "2p,darkorange,-", "2nd° burns (>335 kJ/m²)"),
        (125, "1.5p,orange,-", "1st° burns (>125 kJ/m²)"),
    ]

    for thresh, pen, label in thermal_contours:
        r_test = np.linspace(1, 30000, 30000)
        Q = thermal_fluence_kj(r_test)
        idx = np.argmin(np.abs(Q - thresh))
        radius = r_test[idx]
        lons, lats = radius_to_circle_coords(radius)
        fig.plot(x=lons, y=lats, pen=pen)

    fig.plot(x=[EVENT_LON], y=[EVENT_LAT],
             style="a0.5c", fill="red", pen="1p,black")

    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        f.write("L 7p,Helvetica-Bold,black L ABSORBED DOSE (PROMPT):\n")
        for thresh, pen, label in rad_contours:
            r_test = np.linspace(1, 20000, 20000)
            D = prompt_radiation_gy(r_test)
            idx = np.argmin(np.abs(D - thresh))
            r_km = r_test[idx] / 1000
            f.write(f"S 0.2c - 0.4c - {pen} 0.8c {label} (r={r_km:.1f} km)\n")
        f.write("L 7p,Helvetica-Bold,black L THERMAL FLUENCE:\n")
        for thresh, pen, label in thermal_contours:
            r_test = np.linspace(1, 30000, 30000)
            Q = thermal_fluence_kj(r_test)
            idx = np.argmin(np.abs(Q - thresh))
            r_km = r_test[idx] / 1000
            f.write(f"S 0.2c - 0.4c - {pen} 0.8c {label} (r={r_km:.1f} km)\n")
        f.write("L 6p,Helvetica-Oblique,gray40 L Terrain provides LOS shielding\n")
        legend_file = f.name

    fig.legend(spec=legend_file, position="JBL+w8c+o0.3c/0.3c",
               box="+p0.5p+gwhite@20+s")
    os.unlink(legend_file)

    fig.basemap(region=region, projection=projection,
                map_scale=f"g{EVENT_LON+0.07}/{EVENT_LAT-0.07}+c{EVENT_LAT}+w2k+f+l")

    output = os.path.join(OUTDIR, "map06_radiation_zones.png")
    fig.savefig(output, dpi=300)
    print(f"  Saved {output}")


# ═══════════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    print("=" * 72)
    print("  Berkeley 100 kT Airburst — PyGMT Maps with Terrain Relief")
    print("=" * 72)
    print()

    print("[1/6] Regional overview ...")
    map01_regional_overview()

    print("[2/6] Blast damage zones ...")
    map02_blast_damage()

    print("[3/6] Thermal & radiation zones ...")
    map03_thermal_radiation()

    print("[4/6] Fallout plume ...")
    map04_fallout_plume()

    print("[5/6] Combined effects ...")
    map05_combined_effects()

    print("[6/6] Radiation zones (close-up) ...")
    map06_radiation_zones()

    print()
    print("=" * 72)
    print(f"  All maps saved to: {OUTDIR}")
    print("=" * 72)


if __name__ == "__main__":
    main()
