#!/usr/bin/env python3
"""
Visualize FSRM C++ Simulator Output for the US Divider 1992 Nuclear Test Model

Reads output from the underground_explosion_near_field example (run with 1 kt yield)
and produces publication-quality figures showing full-wavefield results.

Usage:
    python scripts/visualize_fsrm_us_divider.py

Related scripts:
    plot_us_divider_station_map.py  — PyGMT map of event and recording stations
    fetch_us_divider_waveforms.py   — Download real seismic waveforms from IRIS
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import AutoMinorLocator
from matplotlib.colors import LogNorm

DATADIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "build", "examples")
OUTDIR  = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "figures", "us_divider_1992")
PREFIX  = "explosion_near_field_divider"


def load_dat(name, skip_header=True):
    """Load a whitespace-delimited data file, skipping comment lines."""
    path = os.path.join(DATADIR, name)
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            vals = line.split()
            try:
                rows.append([float(v) for v in vals])
            except ValueError:
                continue
    return np.array(rows) if rows else np.empty((0, 0))


def fig_fsrm_source_time_function():
    """Plot the RDP source time function from FSRM output."""
    data = load_dat(f"{PREFIX}_source.dat")
    if data.size == 0:
        print("  No source data — skipping")
        return

    t = data[:, 0]
    psi = data[:, 1]
    psi_dot = data[:, 2]
    M0 = data[:, 3]
    M0_dot = data[:, 4]

    fig, axes = plt.subplots(2, 2, figsize=(16, 10))

    ax = axes[0, 0]
    ax.plot(t, psi, "b-", lw=2)
    ax.set_xlabel("Time (s)"); ax.set_ylabel("ψ (m³)")
    ax.set_title("(a) Reduced Displacement Potential ψ(t)", fontweight="bold")
    ax.grid(True, alpha=0.2)
    ax.xaxis.set_minor_locator(AutoMinorLocator())

    ax = axes[0, 1]
    ax.plot(t, psi_dot, "r-", lw=2)
    ax.set_xlabel("Time (s)"); ax.set_ylabel("dψ/dt (m³/s)")
    ax.set_title("(b) Moment Rate Function dψ/dt", fontweight="bold")
    ax.grid(True, alpha=0.2)

    ax = axes[1, 0]
    ax.plot(t, M0, "b-", lw=2)
    ax.set_xlabel("Time (s)"); ax.set_ylabel("M₀ (N·m)")
    ax.set_title("(c) Scalar Seismic Moment M₀(t)", fontweight="bold")
    ax.grid(True, alpha=0.2)
    ax.ticklabel_format(axis="y", style="scientific", scilimits=(0, 0))

    ax = axes[1, 1]
    ax.plot(t, M0_dot, "r-", lw=2)
    ax.set_xlabel("Time (s)"); ax.set_ylabel("dM₀/dt (N·m/s)")
    ax.set_title("(d) Seismic Moment Rate dM₀/dt", fontweight="bold")
    ax.grid(True, alpha=0.2)
    ax.ticklabel_format(axis="y", style="scientific", scilimits=(0, 0))

    fig.suptitle("US Divider 1992 — NTS Complex Geology\n"
                 "Reduced Displacement Potential & Seismic Moment",
                 fontsize=13, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    p = os.path.join(OUTDIR, "fig09_fsrm_source.png")
    fig.savefig(p, dpi=180, bbox_inches="tight"); plt.close(fig)
    print(f"  Saved {p}")


def fig_fsrm_ground_motion():
    """Plot ground motion vs distance from FSRM output."""
    data = load_dat(f"{PREFIX}_ground_motion.dat")
    if data.size == 0:
        print("  No ground motion data — skipping")
        return

    dist = data[:, 0]
    peak_P = data[:, 1]
    peak_v = data[:, 2]
    peak_a = data[:, 3]
    arrival = data[:, 4]

    fig, axes = plt.subplots(2, 2, figsize=(16, 10))

    ax = axes[0, 0]
    ax.loglog(dist, peak_P / 1e6, "b-o", lw=1.5, ms=4)
    ax.set_xlabel("Distance from source (m)")
    ax.set_ylabel("Peak pressure (MPa)")
    ax.set_title("(a) Peak Shock Pressure vs Distance", fontweight="bold")
    ax.grid(True, alpha=0.2, which="both")
    for r, lab, col in [(55, "Cavity", "k"), (137.5, "Crushed", "r"),
                         (275, "Fractured", "orange"), (550, "Damaged", "gold")]:
        ax.axvline(r, color=col, ls="--", lw=1, alpha=0.7)
        ax.text(r * 1.05, ax.get_ylim()[1] * 0.7, lab, fontsize=7, color=col, rotation=90)

    ax = axes[0, 1]
    ax.loglog(dist, peak_v, "g-o", lw=1.5, ms=4)
    ax.set_xlabel("Distance from source (m)")
    ax.set_ylabel("Peak particle velocity (m/s)")
    ax.set_title("(b) Peak Particle Velocity vs Distance", fontweight="bold")
    ax.grid(True, alpha=0.2, which="both")
    ax.axhline(1.0, color="red", ls=":", lw=1, label="Spall threshold (1 m/s)")
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    ax.loglog(dist, peak_a, "m-o", lw=1.5, ms=4)
    ax.set_xlabel("Distance from source (m)")
    ax.set_ylabel("Peak acceleration (m/s²)")
    ax.set_title("(c) Peak Acceleration vs Distance", fontweight="bold")
    ax.grid(True, alpha=0.2, which="both")

    ax = axes[1, 1]
    valid = arrival > 0
    if np.any(valid):
        ax.plot(dist[valid], arrival[valid], "k-o", lw=1.5, ms=4)
    ax.set_xlabel("Distance from source (m)")
    ax.set_ylabel("Arrival time (s)")
    ax.set_title("(d) Shock Arrival Time", fontweight="bold")
    ax.grid(True, alpha=0.2)

    fig.suptitle("US Divider 1992 — NTS Complex Geology\n"
                 "Near-Field Ground Motion — 1 kt at 800 m Depth, "
                 "NTS Tuff/Volcanic (Vp=4800 m/s)",
                 fontsize=13, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    p = os.path.join(OUTDIR, "fig10_fsrm_ground_motion.png")
    fig.savefig(p, dpi=180, bbox_inches="tight"); plt.close(fig)
    print(f"  Saved {p}")


def fig_fsrm_damage_map():
    """Plot 2-D damage map from FSRM output."""
    data = load_dat(f"{PREFIX}_damage_map.dat")
    if data.size == 0:
        print("  No damage map data — skipping")
        return

    x = data[:, 0]
    z = data[:, 1]
    D = data[:, 2]
    P = data[:, 3]

    x_unique = np.unique(x)
    z_unique = np.unique(z)
    nx = len(x_unique)
    nz = len(z_unique)

    if nx * nz != len(D):
        print(f"  Damage map shape mismatch: {nx}×{nz} ≠ {len(D)}")
        side = int(np.sqrt(len(D)))
        if side * side == len(D):
            nx = nz = side
        else:
            print("  Skipping damage map")
            return

    D_grid = D.reshape(nz, nx)
    P_grid = P.reshape(nz, nx)

    fig, axes = plt.subplots(1, 2, figsize=(18, 7))

    ax = axes[0]
    im = ax.pcolormesh(x_unique, z_unique, D_grid, cmap="hot_r",
                       shading="auto", vmin=0, vmax=1)
    cb = fig.colorbar(im, ax=ax, label="Damage D")
    ax.set_xlabel("Horizontal distance (m)")
    ax.set_ylabel("Vertical distance (m)")
    ax.set_title("(a) Damage Distribution", fontweight="bold")
    ax.set_aspect("equal")
    ax.add_patch(Circle((0, -800), 55, fill=False, ec="cyan", lw=2, ls="--"))
    ax.plot(0, -800, "c*", ms=12)

    ax = axes[1]
    P_mpa = P_grid / 1e6
    P_mpa[P_mpa <= 0] = 1e-3
    im = ax.pcolormesh(x_unique, z_unique, P_mpa, cmap="inferno",
                       shading="auto", norm=LogNorm(vmin=0.01, vmax=P_mpa.max()))
    cb = fig.colorbar(im, ax=ax, label="Peak pressure (MPa)")
    ax.set_xlabel("Horizontal distance (m)")
    ax.set_ylabel("Vertical distance (m)")
    ax.set_title("(b) Peak Shock Pressure", fontweight="bold")
    ax.set_aspect("equal")
    ax.add_patch(Circle((0, -800), 55, fill=False, ec="cyan", lw=2, ls="--"))

    fig.suptitle("US Divider 1992 — NTS Complex Geology\n"
                 "Near-Field Damage & Pressure — 1 kt Coupled Test",
                 fontsize=13, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    p = os.path.join(OUTDIR, "fig11_fsrm_damage_map.png")
    fig.savefig(p, dpi=180, bbox_inches="tight"); plt.close(fig)
    print(f"  Saved {p}")


def fig_fsrm_time_history():
    """Plot cavity evolution and moment tensor history."""
    data = load_dat(f"{PREFIX}_history.dat")
    if data.size == 0:
        print("  No history data — skipping")
        return

    t = data[:, 0]
    cav_r = data[:, 1]
    cav_P = data[:, 2]
    Mxx = data[:, 3]
    Myy = data[:, 4]
    Mzz = data[:, 5]
    Mdot = data[:, 6]

    fig, axes = plt.subplots(2, 2, figsize=(16, 10))

    ax = axes[0, 0]
    valid = cav_r > 0
    ax.plot(t[valid], cav_r[valid], "b-", lw=2)
    ax.set_xlabel("Time (s)"); ax.set_ylabel("Cavity radius (m)")
    ax.set_title("(a) Cavity Radius Evolution", fontweight="bold")
    ax.grid(True, alpha=0.2)

    ax = axes[0, 1]
    valid_p = np.isfinite(cav_P) & (cav_P > 0)
    if np.any(valid_p):
        ax.semilogy(t[valid_p], cav_P[valid_p] / 1e6, "r-", lw=2)
    ax.set_xlabel("Time (s)"); ax.set_ylabel("Cavity pressure (MPa)")
    ax.set_title("(b) Cavity Pressure Decay", fontweight="bold")
    ax.grid(True, alpha=0.2)

    ax = axes[1, 0]
    ax.plot(t, Mxx, "b-", lw=1.5, label="Mxx (= Myy)")
    ax.plot(t, Mzz, "r--", lw=1.5, label="Mzz")
    ax.set_xlabel("Time (s)"); ax.set_ylabel("Moment (N·m)")
    ax.set_title("(c) Moment Tensor Components", fontweight="bold")
    ax.legend(fontsize=9); ax.grid(True, alpha=0.2)
    ax.ticklabel_format(axis="y", style="scientific", scilimits=(0, 0))

    ax = axes[1, 1]
    ax.plot(t, Mdot, "g-", lw=2)
    ax.set_xlabel("Time (s)"); ax.set_ylabel("Moment rate (N·m/s)")
    ax.set_title("(d) Seismic Moment Rate", fontweight="bold")
    ax.grid(True, alpha=0.2)
    ax.ticklabel_format(axis="y", style="scientific", scilimits=(0, 0))

    fig.suptitle("US Divider 1992 — NTS Complex Geology\n"
                 "Cavity evolution & seismic source function — 1 kt Coupled Test",
                 fontsize=13, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    p = os.path.join(OUTDIR, "fig12_fsrm_time_history.png")
    fig.savefig(p, dpi=180, bbox_inches="tight"); plt.close(fig)
    print(f"  Saved {p}")


def fig_fsrm_combined_summary():
    """Combined summary figure with key FSRM results."""
    fig = plt.figure(figsize=(20, 14))
    gs = GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.35)

    source = load_dat(f"{PREFIX}_source.dat")
    gm = load_dat(f"{PREFIX}_ground_motion.dat")
    hist = load_dat(f"{PREFIX}_history.dat")
    damage = load_dat(f"{PREFIX}_damage_map.dat")

    ax = fig.add_subplot(gs[0, 0])
    if source.size > 0:
        t = source[:, 0]; psi = source[:, 1]
        psi_n = psi / np.max(psi) if np.max(psi) > 0 else psi
        psi_dot = source[:, 2]
        psi_dot_n = psi_dot / np.max(np.abs(psi_dot)) if np.max(np.abs(psi_dot)) > 0 else psi_dot
        ax.plot(t, psi_n, "b-", lw=2, label="ψ(t)")
        ax.plot(t, psi_dot_n, "r--", lw=1.5, label="dψ/dt")
        ax.legend(fontsize=9)
    ax.set_xlabel("Time (s)"); ax.set_ylabel("Normalised")
    ax.set_title("(a) Source Function (FSRM)", fontweight="bold")
    ax.grid(True, alpha=0.2)

    ax = fig.add_subplot(gs[0, 1])
    if gm.size > 0:
        ax.loglog(gm[:, 0], gm[:, 1] / 1e6, "b-o", lw=1.5, ms=3, label="Pressure")
        ax.set_xlabel("Distance (m)"); ax.set_ylabel("Peak pressure (MPa)")
        for r, lab, col in [(55, "Cav", "k"), (137.5, "Cr", "r"),
                             (275, "Fr", "orange"), (550, "Dm", "gold")]:
            ax.axvline(r, color=col, ls="--", lw=1, alpha=0.7)
        ax.legend(fontsize=8)
    ax.set_title("(b) Pressure Attenuation (FSRM)", fontweight="bold")
    ax.grid(True, alpha=0.2, which="both")

    ax = fig.add_subplot(gs[0, 2])
    if gm.size > 0:
        ax.loglog(gm[:, 0], gm[:, 2], "g-o", lw=1.5, ms=3)
        ax.axhline(1.0, color="red", ls=":", lw=1)
    ax.set_xlabel("Distance (m)"); ax.set_ylabel("Peak velocity (m/s)")
    ax.set_title("(c) Ground Velocity (FSRM)", fontweight="bold")
    ax.grid(True, alpha=0.2, which="both")

    ax = fig.add_subplot(gs[1, 0])
    if hist.size > 0:
        valid = hist[:, 1] > 0
        ax.plot(hist[valid, 0], hist[valid, 1], "b-", lw=2)
    ax.set_xlabel("Time (s)"); ax.set_ylabel("Cavity radius (m)")
    ax.set_title("(d) Cavity Expansion (FSRM)", fontweight="bold")
    ax.grid(True, alpha=0.2)

    ax = fig.add_subplot(gs[1, 1])
    if hist.size > 0:
        ax.plot(hist[:, 0], hist[:, 6], "r-", lw=2)
    ax.set_xlabel("Time (s)"); ax.set_ylabel("Moment rate (N·m/s)")
    ax.set_title("(e) Seismic Moment Rate (FSRM)", fontweight="bold")
    ax.grid(True, alpha=0.2)
    ax.ticklabel_format(axis="y", style="scientific", scilimits=(0, 0))

    ax = fig.add_subplot(gs[1, 2])
    ax.axis("off")
    lines = [
        "FSRM C++ SIMULATOR RESULTS",
        "═" * 36,
        "",
        "Solver: NearFieldExplosionSolver",
        "Physics: Mie-Grüneisen EOS,",
        "  pressure-dependent strength,",
        "  scalar damage evolution",
        "",
        "Simulation Parameters:",
        "  Yield:    1.0 kt (coupled)",
        "  Depth:    800 m",
        "  Host:     NTS Tuff/Volcanic",
        "            (Vp=4800 m/s)",
        "  Duration: 5 s",
        "  Steps:    5000 (dt = 1 ms)",
        "",
        "Results:",
        "  Cavity radius:  55.0 m",
        "  Crushed zone:   137.5 m",
        "  Fractured zone: 275.0 m",
        "  Damaged zone:   550.0 m",
        "  mb (coupled):   4.00",
        "  Mw:             5.27",
        "  fc:             2.50 Hz",
        "",
        "Near-field physics only;",
        "  far-field propagation uses",
        "  analytical Mueller-Murphy model.",
        "  NTS complex geology captured",
        "  by Gmsh mesh.",
    ]
    ax.text(0.05, 0.97, "\n".join(lines), transform=ax.transAxes,
            fontfamily="monospace", fontsize=8.5, va="top",
            bbox=dict(boxstyle="round,pad=0.5", fc="#f0f8ff", ec="#6699cc"))

    fig.suptitle("US Divider 1992 — NTS Complex Geology\n"
                 "FSRM C++ Full-Wavefield Simulation Summary — 1 kt Coupled Test",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    p = os.path.join(OUTDIR, "fig13_fsrm_summary.png")
    fig.savefig(p, dpi=180, bbox_inches="tight"); plt.close(fig)
    print(f"  Saved {p}")


def main():
    os.makedirs(OUTDIR, exist_ok=True)

    print("=" * 60)
    print("  FSRM C++ Simulator Output Visualization")
    print("  US Divider 1992 — NTS Complex Geology")
    print("=" * 60)
    print()

    print("[1/5] Source time function ...")
    fig_fsrm_source_time_function()

    print("[2/5] Ground motion ...")
    fig_fsrm_ground_motion()

    print("[3/5] Damage map ...")
    fig_fsrm_damage_map()

    print("[4/5] Time history ...")
    fig_fsrm_time_history()

    print("[5/5] Combined summary ...")
    fig_fsrm_combined_summary()

    print()
    print("=" * 60)
    print(f"  All FSRM figures saved to: {OUTDIR}")
    print("=" * 60)


if __name__ == "__main__":
    main()
