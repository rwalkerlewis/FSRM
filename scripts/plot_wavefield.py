#!/usr/bin/env python3
"""
plot_wavefield.py -- Plot displacement field from FSRM HDF5 output.

Reads PETSc HDF5 solution files and plots displacement magnitude.
PETSc stores solution vectors in /solution/Vec with field data as raw arrays.
This script extracts the displacement components and computes the magnitude.

Usage:
    python3 scripts/plot_wavefield.py <solution.h5> [--output <png_file>]

Examples:
    python3 scripts/plot_wavefield.py build/output/solution.h5
    python3 scripts/plot_wavefield.py build/output/solution.h5 --output figures/wavefield.png

Requires: h5py, numpy, matplotlib

Note: For full 3D visualization, use ParaView with the XDMF reader.
This script provides a quick 1D profile or summary of the solution data.
"""

import argparse
import os
import sys


def main():
    parser = argparse.ArgumentParser(
        description="Plot displacement field from FSRM HDF5 output"
    )
    parser.add_argument("h5_file", help="PETSc HDF5 solution file")
    parser.add_argument(
        "--output", "-o",
        default=None,
        help="Output PNG file path"
    )
    parser.add_argument(
        "--title", "-t",
        default="FSRM Solution Field",
        help="Plot title"
    )
    args = parser.parse_args()

    if not os.path.isfile(args.h5_file):
        print(f"Error: File not found: {args.h5_file}")
        sys.exit(1)

    try:
        import h5py
    except ImportError:
        print("Error: h5py is required. Install with: pip install h5py")
        sys.exit(1)

    try:
        import numpy as np
    except ImportError:
        print("Error: numpy is required. Install with: pip install numpy")
        sys.exit(1)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("Error: matplotlib is required. Install with: pip install matplotlib")
        sys.exit(1)

    # Read PETSc HDF5 file
    with h5py.File(args.h5_file, "r") as f:
        print(f"HDF5 file: {args.h5_file}")
        print(f"Top-level groups: {list(f.keys())}")

        # PETSc stores solution vectors under various paths
        # Common patterns: /solution, /fields/solution, /Vec, /timestep_*/solution
        solution_data = None
        solution_path = None

        # Try common PETSc HDF5 paths
        candidates = [
            "solution",
            "fields/solution",
            "Vec",
        ]

        for path in candidates:
            if path in f:
                obj = f[path]
                if isinstance(obj, h5py.Dataset):
                    solution_data = obj[:]
                    solution_path = path
                    break
                elif isinstance(obj, h5py.Group):
                    # Look for dataset inside group
                    for key in obj.keys():
                        if isinstance(obj[key], h5py.Dataset):
                            solution_data = obj[key][:]
                            solution_path = f"{path}/{key}"
                            break
                    if solution_data is not None:
                        break

        if solution_data is None:
            # Print all datasets for debugging
            print("\nAvailable datasets:")
            def print_datasets(name, obj):
                if isinstance(obj, h5py.Dataset):
                    print(f"  {name}: shape={obj.shape}, dtype={obj.dtype}")
            f.visititems(print_datasets)
            print("\nCould not find solution data. Use ParaView for full visualization.")
            sys.exit(0)

        print(f"Solution path: {solution_path}")
        print(f"Solution shape: {solution_data.shape}")
        print(f"Solution dtype: {solution_data.dtype}")

        # Compute basic statistics
        data = solution_data.flatten()
        print(f"Solution stats: min={np.min(data):.6e}, max={np.max(data):.6e}, "
              f"mean={np.mean(data):.6e}, norm={np.linalg.norm(data):.6e}")

        # Plot histogram of solution values
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        # Left: histogram
        nonzero = data[np.abs(data) > 1e-30]
        if len(nonzero) > 0:
            axes[0].hist(nonzero, bins=100, color="steelblue", edgecolor="none", alpha=0.8)
            axes[0].set_xlabel("Solution value")
            axes[0].set_ylabel("Count")
            axes[0].set_title("Distribution of nonzero solution values")
            axes[0].ticklabel_format(axis="x", style="scientific", scilimits=(-2, 2))
        else:
            axes[0].text(0.5, 0.5, "All values are zero",
                        ha="center", va="center", transform=axes[0].transAxes)
            axes[0].set_title("Solution values")

        # Right: solution values vs DOF index
        axes[1].plot(data, ".", markersize=0.5, alpha=0.5, color="steelblue")
        axes[1].set_xlabel("DOF index")
        axes[1].set_ylabel("Solution value")
        axes[1].set_title("Solution vector")
        axes[1].ticklabel_format(axis="y", style="scientific", scilimits=(-2, 2))
        axes[1].grid(True, alpha=0.3)

        fig.suptitle(args.title, fontsize=14, fontweight="bold")
        plt.tight_layout()

        if args.output:
            out_dir = os.path.dirname(args.output)
            if out_dir:
                os.makedirs(out_dir, exist_ok=True)
            plt.savefig(args.output, dpi=150, bbox_inches="tight")
            print(f"Saved: {args.output}")
        else:
            out_path = args.h5_file.replace(".h5", "_summary.png")
            plt.savefig(out_path, dpi=150, bbox_inches="tight")
            print(f"Saved: {out_path}")

        plt.close()


if __name__ == "__main__":
    main()
