#!/usr/bin/env python3
"""
PyGMT Map of the 2017-09-03 DPRK Nuclear Test and Seismic Stations

Generates a publication-quality map showing the DPRK 2017 event location
and surrounding IMS / IRIS seismic stations that recorded the event.

Event:
    2017-09-03 03:30:01 UTC
    Location: 41.300°N, 129.076°E (Punggye-ri test site, DPRK)
    Depth: ~760 m, mb ~6.3, Fully coupled ~250 kt

Output:
    figures/dprk_2017/station_map.png

Requirements:
    pip install pygmt

Usage:
    python scripts/plot_dprk_2017_station_map.py
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
OUTDIR = os.path.join(PROJECT_DIR, "figures", "dprk_2017")
os.makedirs(OUTDIR, exist_ok=True)

# ==============================================================================
# Event Parameters
# ==============================================================================
EVENT_LON = 129.076
EVENT_LAT = 41.300
EVENT_DATE = "2017-09-03"

# ==============================================================================
# Station Data (from IRIS FDSN metadata)
# ==============================================================================
STATIONS = pd.DataFrame({
    "lon":     [129.592, 127.867, 126.624, 116.168, 119.742, 121.187,
                138.204, 143.157, 142.760, 109.494, 107.053,  74.494],
    "lat":     [ 44.616,  36.539,  37.479,  40.018,  49.267,  31.096,
                 36.545,  42.015,  46.959,  30.275,  47.865,  42.638],
    "net_sta": ["IC.MDJ", "KS.KSRS", "KG.INCN", "IC.BJT", "IC.HIA", "IC.SSE",
                "IU.MAJO", "II.ERM", "IU.YSS", "IC.ENH", "IU.ULN", "II.AAK"],
    "dist_km": [370, 520, 460, 1130, 930, 1250, 1050, 1300, 750, 2000, 1900, 4200],
    "desc":    ["Mudanjiang", "KSRS Array", "Incheon", "Baijiatuan", "Hailar", "Shanghai",
                "Matsushiro", "Erimo", "Yuzhno", "Enshi", "Ulaanbaatar", "Ala Archa"],
})


def plot_station_map():
    """Generate the station map using PyGMT."""

    # Map region: East Asia
    region = [100, 155, 25, 55]
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
        frame=["WSne+tDPRK 2017-09-03 Nuclear Test and Seismic Stations",
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
        text="Punggye-ri",
        font="8p,Helvetica-Bold,red",
        justify="TC",
    )

    # ── Event info box ──
    # Create a text box with event details in the top-left corner
    info_lines = [
        "Event Information:",
        f"Date: {EVENT_DATE}  03:30:01 UTC",
        f"Location: {EVENT_LAT:.3f}@.N, {EVENT_LON:.3f}@.E",
        "Depth: ~760 m",
        "mb @~\\176@~ 6.3",
        "Fully coupled ~250 kt",
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
        map_scale="g107/27+c41+w500k+f+l",
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
