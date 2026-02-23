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
import xarray as xr

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
# Terrain-Aware Damage Functions
# ═══════════════════════════════════════════════════════════════════════════════

def load_terrain_grid(region, resolution="03s"):
    """
    Load DEM data for a region using PyGMT's earth_relief.
    
    Args:
        region: [lon_min, lon_max, lat_min, lat_max]
        resolution: Grid resolution (e.g., "03s" for 3 arc-seconds)
    
    Returns:
        xarray DataArray with elevation in meters
    """
    grid = pygmt.datasets.load_earth_relief(resolution=resolution, region=region)
    return grid


def compute_line_of_sight_factor(dem_grid, burst_lon=EVENT_LON, burst_lat=EVENT_LAT,
                                  burst_height=BURST_HEIGHT):
    """
    Compute line-of-sight visibility from burst point to each grid cell.
    Fully vectorized — processes entire grid in one pass per ray sample.

    Returns:
        2D numpy array of visibility factors: 1 = visible, 0 = blocked.
    """
    lons = dem_grid.lon.values
    lats = dem_grid.lat.values
    elev = dem_grid.values
    ny, nx = elev.shape

    burst_elev = float(dem_grid.interp(lon=burst_lon, lat=burst_lat).values)
    burst_z = burst_elev + burst_height

    # Grid spacing
    dlon = lons[1] - lons[0] if nx > 1 else 1.0
    dlat = lats[1] - lats[0] if ny > 1 else 1.0

    # Burst in fractional grid indices
    burst_i = (burst_lon - lons[0]) / dlon
    burst_j = (burst_lat - lats[0]) / dlat

    # Target grid indices
    target_i, target_j = np.meshgrid(np.arange(nx, dtype=float),
                                      np.arange(ny, dtype=float))

    los_factor = np.ones((ny, nx), dtype=float)
    N_SAMPLES = 25

    for t in np.linspace(0.1, 0.9, N_SAMPLES):
        # Sample position along ray for every pixel simultaneously
        sample_i = burst_i + t * (target_i - burst_i)
        sample_j = burst_j + t * (target_j - burst_j)

        sample_i = np.clip(sample_i, 0, nx - 1.001)
        sample_j = np.clip(sample_j, 0, ny - 1.001)

        i0 = sample_i.astype(int)
        j0 = sample_j.astype(int)
        i1 = np.minimum(i0 + 1, nx - 1)
        j1 = np.minimum(j0 + 1, ny - 1)

        di = sample_i - i0
        dj = sample_j - j0

        # Bilinear interpolation of terrain
        sample_elev = (
            (1 - di) * (1 - dj) * elev[j0, i0] +
            di * (1 - dj) * elev[j0, i1] +
            (1 - di) * dj * elev[j1, i0] +
            di * dj * elev[j1, i1]
        )

        # Ray height at parameter t
        ray_height = burst_z + t * (elev - burst_z)

        # Block where terrain exceeds ray
        los_factor[sample_elev > ray_height] = 0.0

    return los_factor


def compute_blast_terrain_factor(dem_grid, burst_lon=EVENT_LON, burst_lat=EVENT_LAT,
                                  burst_height=BURST_HEIGHT):
    """
    Compute terrain modification factor for blast waves.
    Fully vectorized — no per-pixel loops.

    Returns 2D array of factors to multiply base overpressure by.
    """
    lons = dem_grid.lon.values
    lats = dem_grid.lat.values
    elev = dem_grid.values
    ny, nx = elev.shape

    burst_elev = float(dem_grid.interp(lon=burst_lon, lat=burst_lat).values)
    burst_z = burst_elev + burst_height

    dlat = lats[1] - lats[0] if ny > 1 else 0.0001
    dlon = lons[1] - lons[0] if nx > 1 else 0.0001
    cos_lat = np.cos(np.radians(burst_lat))

    grad_lon = np.gradient(elev, axis=1) / (dlon * 111320.0 * cos_lat)
    grad_lat = np.gradient(elev, axis=0) / (dlat * 111320.0)

    lon_grid, lat_grid = np.meshgrid(lons, lats)
    dx = (lon_grid - burst_lon) * 111320.0 * cos_lat
    dy = (lat_grid - burst_lat) * 111320.0
    dist = np.sqrt(dx**2 + dy**2)
    dist_safe = np.maximum(dist, 1.0)

    ux = dx / dist_safe
    uy = dy / dist_safe

    slope_toward_blast = -(grad_lon * ux + grad_lat * uy)
    height_diff = elev - burst_z

    blast_factor = np.ones((ny, nx), dtype=float)
    far = dist > 50

    # Forward-facing slopes — enhanced reflection
    m_fwd = far & (slope_toward_blast > 0.1)
    blast_factor = np.where(m_fwd,
                            1.0 + np.minimum(slope_toward_blast * 2, 1.0),
                            blast_factor)

    # Behind a ridge — significant shielding
    m_ridge = far & ~m_fwd & (height_diff > 100) & (slope_toward_blast < -0.1)
    blast_factor = np.where(m_ridge,
                            0.3 + 0.4 * np.exp(-height_diff / 300),
                            blast_factor)

    # Elevated terrain with mild shielding
    m_elev = far & ~m_fwd & ~m_ridge & (height_diff > 50)
    blast_factor = np.where(m_elev, 0.7, blast_factor)

    return blast_factor


def compute_terrain_aware_effects(region, resolution="03s"):
    """
    Compute terrain-aware damage grids for all effects.
    
    Returns dict with gridded damage values accounting for topography.
    """
    print("    Loading DEM data...")
    dem = load_terrain_grid(region, resolution)
    
    lons = dem.lon.values
    lats = dem.lat.values
    elev = dem.values
    ny, nx = elev.shape
    
    print(f"    Computing effects on {nx}x{ny} grid...")
    
    # Get burst ground elevation
    burst_elev = float(dem.interp(lon=EVENT_LON, lat=EVENT_LAT).values)
    
    # Compute distances and base effects
    lon_grid, lat_grid = np.meshgrid(lons, lats)
    dx = (lon_grid - EVENT_LON) * 111320.0 * np.cos(np.radians(EVENT_LAT))
    dy = (lat_grid - EVENT_LAT) * 111320.0
    ground_range = np.sqrt(dx**2 + dy**2)
    
    # Slant range accounting for terrain elevation relative to burst
    dz = elev - (burst_elev + BURST_HEIGHT)
    slant_range = np.sqrt(ground_range**2 + dz**2)
    
    print("    Computing line-of-sight visibility...")
    los_factor = compute_line_of_sight_factor(dem)
    
    print("    Computing blast terrain factors...")
    blast_terrain = compute_blast_terrain_factor(dem)
    
    # Base overpressure (flat terrain)
    base_overpressure = peak_overpressure_kpa(ground_range)
    
    # Terrain-modified overpressure
    terrain_overpressure = base_overpressure * blast_terrain
    
    # Thermal fluence with LOS blocking
    base_thermal = thermal_fluence_kj(slant_range)
    terrain_thermal = base_thermal * los_factor
    
    # Prompt radiation with LOS blocking
    base_radiation = prompt_radiation_gy(slant_range)
    terrain_radiation = base_radiation * los_factor
    
    return {
        'lons': lons,
        'lats': lats,
        'elevation': elev,
        'ground_range': ground_range,
        'slant_range': slant_range,
        'los_factor': los_factor,
        'blast_terrain_factor': blast_terrain,
        'overpressure_flat': base_overpressure,
        'overpressure_terrain': terrain_overpressure,
        'thermal_flat': base_thermal,
        'thermal_terrain': terrain_thermal,
        'radiation_flat': base_radiation,
        'radiation_terrain': terrain_radiation,
        'dem': dem,
    }


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
                map_scale="jBC+o0c/-1.2c+c37.8+w20k+f+l")

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
                map_scale=f"jBC+o0c/-1.2c+c{EVENT_LAT}+w5k+f+l")

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
                map_scale=f"jBC+o0c/-1.2c+c{EVENT_LAT}+w5k+f+l")

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
                map_scale="jBC+o0c/-1.2c+c37.7+w20k+f+l")

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
                map_scale=f"jBC+o0c/-1.2c+c{EVENT_LAT}+w10k+f+l")

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
                map_scale=f"jBC+o0c/-1.2c+c{EVENT_LAT}+w2k+f+l")

    output = os.path.join(OUTDIR, "map06_radiation_zones.png")
    fig.savefig(output, dpi=300)
    print(f"  Saved {output}")


# ═══════════════════════════════════════════════════════════════════════════════
# Map 7: Terrain-Aware Blast Damage
# ═══════════════════════════════════════════════════════════════════════════════

def map07_blast_terrain():
    """
    Blast damage zones with actual terrain interaction.
    
    Shows how Berkeley/Oakland Hills modify blast wave propagation:
    - Enhanced overpressure on forward-facing slopes
    - Reduced overpressure in terrain shadows
    - Valley channeling effects
    """
    half = 0.16
    region = [EVENT_LON - half, EVENT_LON + half,
              EVENT_LAT - half * 0.8, EVENT_LAT + half * 0.8]
    projection = "M18c"
    
    print("    Computing terrain-aware blast effects...")
    effects = compute_terrain_aware_effects(region, resolution="03s")
    
    # Create grid for PyGMT
    lons = effects['lons']
    lats = effects['lats']
    overpressure = effects['overpressure_terrain']
    terrain_factor = effects['blast_terrain_factor']
    
    fig = pygmt.Figure()
    
    # First panel: terrain-modified overpressure
    fig.grdimage(grid="@earth_relief_03s", region=region,
                 projection=projection, shading=True, cmap="geo")
    
    fig.basemap(region=region, projection=projection,
                frame=["WSne+tBlast Damage with Terrain Effects — 100 kT (50 m HOB)",
                       "xa0.05f0.025", "ya0.05f0.025"])
    
    fig.coast(region=region, projection=projection, resolution="f",
              water="lightblue", shorelines="0.3p,gray30")
    
    # Plot overpressure contours with terrain effects
    damage_levels = [
        (140.0, "darkred",   "2.5p", "Total destruction (>140 kPa)"),
        (35.0,  "red",       "2p", "Severe damage (>35 kPa)"),
        (14.0,  "orange",    "1.5p", "Buildings destroyed (>14 kPa)"),
        (7.0,   "yellow",    "1.5p", "Moderate damage (>7 kPa)"),
        (3.5,   "green",     "1p", "Light damage (>3.5 kPa)"),
    ]
    
    # Create a temporary grid file for contouring
    with tempfile.NamedTemporaryFile(suffix='.nc', delete=False) as f:
        grid_file = f.name
    
    ds = xr.DataArray(overpressure, coords=[lats, lons], dims=['lat', 'lon'])
    ds.to_netcdf(grid_file)
    
    for thresh, color, pen, label in damage_levels:
        try:
            fig.grdcontour(grid=grid_file, levels=thresh, limit=[thresh, thresh*1.001],
                          pen=f"{pen},{color}", region=region, projection=projection)
        except:
            pass  # Contour level may not exist
    
    os.unlink(grid_file)
    
    # Mark ground zero
    fig.plot(x=[EVENT_LON], y=[EVENT_LAT],
             style="a0.5c", fill="red", pen="1p,black")
    fig.text(x=EVENT_LON, y=EVENT_LAT + 0.008,
             text="GZ (50m HOB)", font="8p,Helvetica-Bold,red", justify="BC")
    
    # Mark terrain features
    terrain_labels = [
        (-122.200, 37.876, "Berkeley\nHills"),
        (-122.190, 37.830, "Oakland\nHills"),
        (-122.240, 37.860, "Claremont\nCanyon"),
    ]
    for tlon, tlat, tname in terrain_labels:
        if region[0] < tlon < region[1] and region[2] < tlat < region[3]:
            fig.text(x=tlon, y=tlat, text=tname,
                     font="6p,Helvetica-Oblique,brown", justify="MC")
    
    # Legend
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        f.write("L 7p,Helvetica-Bold,black L Terrain-Modified Overpressure:\n")
        for thresh, color, pen, label in damage_levels:
            f.write(f"S 0.2c - 0.4c - {pen},{color} 0.8c {label}\n")
        f.write("L 6p,Helvetica-Oblique,gray40 L\n")
        f.write("L 6p,Helvetica-Oblique,gray40 L Hills: shielding behind, enhancement on slopes\n")
        f.write("L 6p,Helvetica-Oblique,gray40 L Valleys: channeling may focus blast\n")
        legend_file = f.name
    
    fig.legend(spec=legend_file, position="JBL+w8.5c+o0.3c/0.3c",
               box="+p0.5p+gwhite@20+s")
    os.unlink(legend_file)
    
    fig.basemap(region=region, projection=projection,
                map_scale=f"jBC+o0c/-1.2c+c{EVENT_LAT}+w5k+f+l")
    
    output = os.path.join(OUTDIR, "map07_blast_terrain.png")
    fig.savefig(output, dpi=300)
    print(f"  Saved {output}")


# ═══════════════════════════════════════════════════════════════════════════════
# Map 8: Terrain-Aware Combined Effects
# ═══════════════════════════════════════════════════════════════════════════════

def map08_combined_terrain():
    """
    Combined effects map with terrain interaction.
    
    Shows:
    - Blast overpressure modified by terrain slope and shielding
    - Thermal radiation blocked by hills (line-of-sight)
    - Prompt radiation blocked by hills (line-of-sight)
    - Dominant effect at each location
    """
    half = 0.18
    region = [EVENT_LON - half, EVENT_LON + half,
              EVENT_LAT - half * 0.85, EVENT_LAT + half * 0.85]
    projection = "M18c"
    
    print("    Computing terrain-aware combined effects...")
    effects = compute_terrain_aware_effects(region, resolution="03s")
    
    lons = effects['lons']
    lats = effects['lats']
    P = effects['overpressure_terrain']
    Q = effects['thermal_terrain']
    D = effects['radiation_terrain']
    los = effects['los_factor']
    
    # Compute hazard indices
    haz_blast = np.clip(P / 35.0, 0, 1)
    haz_thermal = np.clip(Q / 670.0, 0, 1)
    haz_rad = np.clip(D / 6.0, 0, 1)
    haz_combined = np.maximum(np.maximum(haz_blast, haz_thermal), haz_rad)
    
    # Determine dominant effect
    dominant = np.zeros_like(P)
    dominant[haz_blast >= np.maximum(haz_thermal, haz_rad)] = 1
    dominant[haz_thermal >= np.maximum(haz_blast, haz_rad)] = 2
    dominant[haz_rad >= np.maximum(haz_blast, haz_thermal)] = 3
    dominant[haz_combined < 0.01] = 0
    
    fig = pygmt.Figure()
    
    # Background terrain
    fig.grdimage(grid="@earth_relief_03s", region=region,
                 projection=projection, shading=True, cmap="geo")
    
    fig.basemap(region=region, projection=projection,
                frame=["WSne+tCombined Effects with Terrain — 100 kT (50 m HOB)",
                       "xa0.05f0.025", "ya0.05f0.025"])
    
    fig.coast(region=region, projection=projection, resolution="f",
              water="lightblue", shorelines="0.3p,gray30")
    
    # Create grid files for plotting
    
    # Plot LOS shielded areas (where thermal/radiation is blocked)
    with tempfile.NamedTemporaryFile(suffix='.nc', delete=False) as f:
        los_file = f.name
    ds_los = xr.DataArray(1.0 - los, coords=[lats, lons], dims=['lat', 'lon'])
    ds_los.to_netcdf(los_file)
    
    # Areas with LOS shielding (>50% blocked)
    try:
        fig.grdcontour(grid=los_file, levels=0.5, limit=[0.5, 1.0],
                      pen="1p,purple,-", region=region, projection=projection,
                      annotation="0.5+f6p+gwhite")
    except:
        pass
    os.unlink(los_file)
    
    # Plot combined hazard contours
    with tempfile.NamedTemporaryFile(suffix='.nc', delete=False) as f:
        haz_file = f.name
    ds_haz = xr.DataArray(haz_combined, coords=[lats, lons], dims=['lat', 'lon'])
    ds_haz.to_netcdf(haz_file)
    
    hazard_levels = [
        (0.9, "darkred", "2p", "Very High (>90%)"),
        (0.5, "red",     "1.5p", "High (>50%)"),
        (0.25, "orange", "1p", "Moderate (>25%)"),
        (0.1, "yellow",  "0.5p", "Low (>10%)"),
    ]
    
    for thresh, color, pen, label in hazard_levels:
        try:
            fig.grdcontour(grid=haz_file, levels=thresh, limit=[thresh, thresh*1.001],
                          pen=f"{pen},{color}", region=region, projection=projection)
        except:
            pass
    os.unlink(haz_file)
    
    # Mark ground zero
    fig.plot(x=[EVENT_LON], y=[EVENT_LAT],
             style="a0.6c", fill="red", pen="1p,black")
    
    # Legend
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        f.write("L 7p,Helvetica-Bold,black L Combined Hazard Index:\n")
        for thresh, color, pen, label in hazard_levels:
            f.write(f"S 0.2c - 0.4c - {pen},{color} 0.8c {label}\n")
        f.write("L 7p,Helvetica-Bold,black L\n")
        f.write("L 7p,Helvetica-Bold,black L Terrain Effects:\n")
        f.write("S 0.2c - 0.4c - 1p,purple,- 0.8c LOS shielding (thermal/rad)\n")
        f.write("L 6p,Helvetica-Oblique,gray40 L Hills block thermal & radiation\n")
        f.write("L 6p,Helvetica-Oblique,gray40 L Blast diffracts but is weakened\n")
        legend_file = f.name
    
    fig.legend(spec=legend_file, position="JBL+w7.5c+o0.3c/0.3c",
               box="+p0.5p+gwhite@20+s")
    os.unlink(legend_file)
    
    fig.basemap(region=region, projection=projection,
                map_scale=f"jBC+o0c/-1.2c+c{EVENT_LAT}+w5k+f+l")
    
    output = os.path.join(OUTDIR, "map08_combined_terrain.png")
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

    print("[1/8] Regional overview ...")
    map01_regional_overview()

    print("[2/8] Blast damage zones (flat terrain) ...")
    map02_blast_damage()

    print("[3/8] Thermal & radiation zones ...")
    map03_thermal_radiation()

    print("[4/8] Fallout plume ...")
    map04_fallout_plume()

    print("[5/8] Combined effects (flat terrain) ...")
    map05_combined_effects()

    print("[6/8] Radiation zones (close-up) ...")
    map06_radiation_zones()

    print("[7/8] Blast damage with terrain interaction ...")
    map07_blast_terrain()

    print("[8/8] Combined effects with terrain interaction ...")
    map08_combined_terrain()

    print()
    print("=" * 72)
    print(f"  All maps saved to: {OUTDIR}")
    print("=" * 72)


if __name__ == "__main__":
    main()
