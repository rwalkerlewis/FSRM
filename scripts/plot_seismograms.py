#!/usr/bin/env python3
"""
plot_seismograms.py -- Plot SAC seismograms produced by FSRM SeismometerNetwork.

Reads all .sac files from a directory and produces a multi-panel figure
with one trace per station-component pair.

Usage:
    python3 scripts/plot_seismograms.py <sac_directory> [--output <png_file>]

Examples:
    python3 scripts/plot_seismograms.py build/output/seismograms/
    python3 scripts/plot_seismograms.py build/output/punggye_ri/ --output figures/seismograms.png

Requires: obspy, matplotlib
"""

import argparse
import os
import sys
import glob

def main():
    parser = argparse.ArgumentParser(
        description="Plot SAC seismograms from FSRM output"
    )
    parser.add_argument("sac_dir", help="Directory containing .sac files")
    parser.add_argument(
        "--output", "-o",
        default=None,
        help="Output PNG file path (default: display or auto-name)"
    )
    parser.add_argument(
        "--title", "-t",
        default="FSRM Seismograms",
        help="Plot title"
    )
    args = parser.parse_args()

    if not os.path.isdir(args.sac_dir):
        print(f"Error: Directory not found: {args.sac_dir}")
        sys.exit(1)

    sac_files = sorted(glob.glob(os.path.join(args.sac_dir, "*.sac")))
    if not sac_files:
        print(f"Error: No .sac files found in {args.sac_dir}")
        sys.exit(1)

    try:
        from obspy import read as obspy_read
    except ImportError:
        print("Error: obspy is required. Install with: pip install obspy")
        sys.exit(1)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("Error: matplotlib is required. Install with: pip install matplotlib")
        sys.exit(1)

    # Read all SAC files
    traces = []
    for sac_file in sac_files:
        try:
            st = obspy_read(sac_file, format="SAC")
            for tr in st:
                traces.append(tr)
        except Exception as e:
            print(f"Warning: Could not read {sac_file}: {e}")

    if not traces:
        print("Error: No valid traces found")
        sys.exit(1)

    # Sort traces by station name, then component
    traces.sort(key=lambda tr: (tr.stats.station, tr.stats.channel))

    n_traces = len(traces)
    fig, axes = plt.subplots(n_traces, 1, figsize=(12, 2.5 * n_traces), sharex=True)
    if n_traces == 1:
        axes = [axes]

    for i, tr in enumerate(traces):
        t = tr.times()
        axes[i].plot(t, tr.data, "k-", linewidth=0.5)
        label = f"{tr.stats.station}.{tr.stats.channel}"
        axes[i].set_ylabel(label, fontsize=9)
        axes[i].ticklabel_format(axis="y", style="scientific", scilimits=(-2, 2))
        axes[i].grid(True, alpha=0.3)

    axes[-1].set_xlabel("Time (s)")
    axes[0].set_title(args.title)
    plt.tight_layout()

    if args.output:
        out_dir = os.path.dirname(args.output)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        plt.savefig(args.output, dpi=150, bbox_inches="tight")
        print(f"Saved: {args.output}")
    else:
        out_path = os.path.join(args.sac_dir, "seismograms.png")
        plt.savefig(out_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {out_path}")

    plt.close()


if __name__ == "__main__":
    main()
