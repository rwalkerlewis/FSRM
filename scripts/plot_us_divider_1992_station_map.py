#!/usr/bin/env python3
"""
PyGMT Map of the US Divider Nuclear Test (1992-09-23) and Seismic Stations

Generates a publication-quality map showing the NTS Divider event location
and surrounding CI TERRAscope, BK, NN/LN, and IU seismic stations.

Event:
    1992-09-23 15:04:00 UTC
    Location: 37.021°N, 116.058°W (NTS Area 4, Nevada)
    mb ~4.0, depth ~800 m
    Last US underground nuclear test before the testing moratorium

Output:
    figures/us_divider_1992/station_map.png

Requirements:
    pip install pygmt

Usage:
    python scripts/plot_us_divider_1992_station_map.py
"""

import os
import tempfile
import pygmt
import pandas as pd
import numpy as np

# ==============================================================================
# Directories
# ==============================================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)
OUTDIR = os.path.join(PROJECT_DIR, "figures", "us_divider_1992")
os.makedirs(OUTDIR, exist_ok=True)

# ==============================================================================
# Event Parameters
# ==============================================================================
EVENT_LON = -116.058
EVENT_LAT = 37.021
EVENT_DATE = "1992-09-23"

# ==============================================================================
# Station Data
#
# 14 stations: CI TERRAscope (6), BK (1), NN (2), LN (1), IU (1),
# plus 3 more distant stations for regional coverage.
#
# Coordinates are approximate IRIS FDSN metadata values for the 1992 epoch.
# ==============================================================================
STATIONS = pd.DataFrame({
    "lon":     [-118.1711, -116.8056, -118.4744, -116.4553, -119.7147,
                -117.0978, -121.4472, -115.2390, -118.1542, -106.4572,
                -116.6722, -118.4039, -110.7847, -123.3046],
    "lat":     [34.1484, 35.3017, 35.6628, 33.6107, 34.4408,
                34.1086, 36.7639, 40.7448, 38.4328, 34.9459,
                32.6800, 33.7433, 32.3098, 44.5855],
    "net_sta": ["CI.PAS", "CI.GSC", "CI.ISA", "CI.PFO", "CI.SBC",
                "CI.SVD", "BK.SAO", "NN.ELK", "LN.MNV", "IU.ANMO",
                "CI.BAR", "CI.RPV", "US.TUC", "IU.COR"],
    "dist_km": [330, 197, 195, 385, 415,
                340, 510, 418, 205, 940,
                495, 385, 730, 1030],
    "desc":    ["Pasadena (TERRAscope)", "Goldstone (TERRAscope)",
                "Isabella (TERRAscope)", "Piñon Flat (TERRAscope)",
                "Santa Barbara (TERRAscope)", "Seven Oaks Dam (TERRAscope)",
                "San Andreas Obs", "Elko", "Mina",
                "Albuquerque", "Barrett", "Rancho Palos Verdes",
                "Tucson", "Corvallis"],
})


def plot_station_map():
    """Generate the station map using PyGMT."""

    region = [-125, -100, 30, 45]
    projection = "M18c"

    fig = pygmt.Figure()

    # ── Topographic relief ──
    fig.grdimage(
        grid="@earth_relief_30s",
        region=region,
        projection=projection,
        shading=True,
        cmap="geo",
    )

    # ── Basemap frame ──
    fig.basemap(
        region=region,
        projection=projection,
        frame=["WSne+tUS Divider Test 1992-09-23 and Seismic Stations",
               "xa5f1g5", "ya5f1g5"],
    )

    # ── Coastlines and borders ──
    fig.coast(
        region=region,
        projection=projection,
        resolution="l",
        water="lightblue",
        borders=["1/0.8p,gray30", "2/0.3p,gray50"],
        shorelines="thin,gray30",
        area_thresh=2000,
    )

    # ── Great-circle paths from event to stations ──
    for _, sta in STATIONS.iterrows():
        fig.plot(
            x=[EVENT_LON, sta["lon"]],
            y=[EVENT_LAT, sta["lat"]],
            pen="0.5p,gray70,--",
        )

    # ── Plot stations as triangles ──
    fig.plot(
        x=STATIONS["lon"],
        y=STATIONS["lat"],
        style="t0.35c",
        fill="blue",
        pen="0.5p,black",
    )

    # ── Plot event as star ──
    fig.plot(
        x=[EVENT_LON],
        y=[EVENT_LAT],
        style="a0.5c",
        fill="red",
        pen="0.5p,black",
    )

    # ── Station labels ──
    for _, sta in STATIONS.iterrows():
        fig.text(
            x=sta["lon"],
            y=sta["lat"] + 0.45,
            text=sta["net_sta"],
            font="6p,Helvetica,black",
            justify="BC",
        )

    # ── Event label ──
    fig.text(
        x=EVENT_LON,
        y=EVENT_LAT - 0.6,
        text="NTS Divider",
        font="8p,Helvetica-Bold,red",
        justify="TC",
    )

    # ── Event info box ──
    info_lines = [
        "Event Information:",
        f"Date: {EVENT_DATE}  15:04 UTC",
        f"Location: {EVENT_LAT:.3f}@.N, {abs(EVENT_LON):.3f}@.W",
        "Site: NTS Area 4",
        "Depth: @~\\176@~ 800 m",
        "mb @~\\176@~ 4.0",
        "Last US underground test",
        "Yield: <20 kt (classified)",
    ]

    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        for i, line in enumerate(info_lines):
            if i == 0:
                f.write(f"L 8p,Helvetica-Bold,black L {line}\n")
            else:
                f.write(f"L 7p,Helvetica,black L {line}\n")
        info_file = f.name

    fig.legend(
        spec=info_file,
        position="JTL+w6c+o0.5c/0.5c",
        box="+p0.5p+gwhite@30+s",
    )
    os.unlink(info_file)

    # ── Legend ──
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        f.write("S 0.2c a 0.4c red 0.5p,black 0.8c Event Origin (NTS)\n")
        f.write("S 0.2c t 0.3c blue 0.5p,black 0.8c Seismic Station\n")
        legend_file = f.name

    fig.legend(
        spec=legend_file,
        position="JBL+w5c+o0.5c/0.5c",
        box="+p0.5p+gwhite+s",
    )
    os.unlink(legend_file)

    # ── Scale bar ──
    fig.basemap(
        region=region,
        projection=projection,
        map_scale="g-103/31+c37+w500k+f+l",
    )

    # ── Inset globe showing western US context ──
    with fig.inset(position="jTR+w5c+o0.3c", box="+p0.5p+gwhite"):
        fig.coast(
            region="g",
            projection=f"G{EVENT_LON}/{EVENT_LAT}/?",
            resolution="c",
            land="gray80",
            water="lightblue",
            shorelines="faint",
            area_thresh=10000,
        )
        fig.plot(
            x=[EVENT_LON],
            y=[EVENT_LAT],
            style="a0.3c",
            fill="red",
            pen="0.3p,black",
        )

    # Save figure
    output_path = os.path.join(OUTDIR, "station_map.png")
    fig.savefig(output_path, dpi=300)
    print(f"Map saved to: {output_path}")

    return fig


if __name__ == "__main__":
    plot_station_map()
