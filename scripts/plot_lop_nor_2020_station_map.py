#!/usr/bin/env python3
"""
PyGMT Map of the 2020-06-22 Lop Nor Event and Seismic Stations

Generates a publication-quality map showing the Lop Nor 2020 event location
and surrounding IMS / IRIS seismic stations that recorded the event.

Event:
    2020-06-22 ~09:18 UTC
    Location: 41.735°N, 88.730°E (Lop Nor test site, Xinjiang, China)
    Detection: PS23 Makanchi array (Kazakhstan), mb ≈ 2.75

Output:
    figures/lop_nor_2020/station_map.png

Requirements:
    pip install pygmt

Usage:
    python scripts/plot_lop_nor_2020_station_map.py
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
OUTDIR = os.path.join(PROJECT_DIR, "figures", "lop_nor_2020")
os.makedirs(OUTDIR, exist_ok=True)

# ==============================================================================
# Event Parameters
# ==============================================================================
EVENT_LON = 88.730
EVENT_LAT = 41.735
EVENT_DATE = "2020-06-22"

# ==============================================================================
# Station Data (from IRIS FDSN metadata)
# ==============================================================================
STATIONS = pd.DataFrame({
    "lon":     [87.6951, 81.9770, 82.2904, 79.2389, 76.9669, 78.4373,
                76.9661, 74.4942, 78.6189, 91.1500, 73.2686, 108.9230],
    "lat":     [43.8144, 46.7928, 46.7939, 41.5508, 43.2220, 42.4665,
                43.2194, 42.6389, 50.7154, 29.7025, 33.6506, 34.0310],
    "net_sta": ["IC.WMQ", "IU.MAKZ", "KZ.MKAR", "G.WUS", "KZ.PDGK", "KR.PRZ",
                "KZ.KNDC", "II.AAK", "II.KURK", "IC.LSA", "II.NIL", "IC.XAN"],
    "dist_km": [246, 780, 778, 839, 1020, 1000, 1020, 1179, 1264, 1350, 1370, 1965],
    "desc":    ["Urumqi", "PS23 Makanchi", "Makanchi Array", "Wushi", "Podgonoye",
                "Karakol", "Almaty", "Ala Archa", "Kurchatov", "Lhasa", "Nilore", "Xi'an"],
})


def plot_station_map():
    """Generate the station map using PyGMT."""
    
    # Map region: Central Asia
    region = [65, 115, 25, 55]
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
        frame=["WSne+tLop Nor 2020-06-22 Event and Seismic Stations",
               "xa10f5g10", "ya10f5g10"],
    )
    
    # ── Coastlines and borders (overlay on relief) ──
    fig.coast(
        region=region,
        projection=projection,
        resolution="l",
        water="lightblue",
        borders="1/0.5p,gray40",
        shorelines="thin,gray30",
        area_thresh=5000,
    )
    
    # ── Plot great-circle paths from event to stations ──
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
            y=sta["lat"] + 0.6,
            text=sta["net_sta"],
            font="6p,Helvetica,black",
            justify="BC",
        )
    
    # ── Event label ──
    fig.text(
        x=EVENT_LON,
        y=EVENT_LAT - 0.8,
        text="Lop Nor",
        font="8p,Helvetica-Bold,red",
        justify="TC",
    )
    
    # ── Event info box ──
    # Create a text box with event details in the top-left corner
    info_lines = [
        "Event Information:",
        f"Date: {EVENT_DATE}  ~09:18 UTC",
        f"Location: {EVENT_LAT:.3f}@.N, {EVENT_LON:.3f}@.E",
        "Depth: ~300 m (tunnel)",
        "mb @~\\176@~ 2.75 (PS23 Makanchi)",
        "Decoupled underground test",
    ]
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        for i, line in enumerate(info_lines):
            # First line bold, rest normal
            if i == 0:
                f.write(f"L 8p,Helvetica-Bold,black L {line}\n")
            else:
                f.write(f"L 7p,Helvetica,black L {line}\n")
        info_file = f.name
    
    fig.legend(
        spec=info_file,
        position="JTL+w5.5c+o0.5c/0.5c",
        box="+p0.5p+gwhite@30+s",
    )
    os.unlink(info_file)
    
    # ── Legend ──
    # Create a temporary legend spec file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        f.write("S 0.2c a 0.4c red 0.5p,black 0.8c Event Origin\n")
        f.write("S 0.2c t 0.3c blue 0.5p,black 0.8c Seismic Station\n")
        legend_file = f.name
    
    fig.legend(
        spec=legend_file,
        position="JBL+w4c+o0.5c/0.5c",
        box="+p0.5p+gwhite+s",
    )
    os.unlink(legend_file)
    
    # ── Scale bar ──
    fig.basemap(
        region=region,
        projection=projection,
        map_scale="g72/27+c41+w500k+f+l",
    )
    
    # ── Inset map showing global location ──
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
