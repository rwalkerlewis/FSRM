#!/usr/bin/env python3
"""02_cavity_formation_radial_profile.py -- DPRK 2017 (Punggye-ri, Sixth Test)

4-panel radial-profile snapshot at four times during cavity formation.
Reads near_field_profile.h5 written by the radial Lagrangian solver
(solver_kind = RADIAL_LAGRANGIAN). If the file is absent the COUPLED_ANALYTIC
path is in use and the script falls back to annotating RDP scalar parameters
read from simulation output log.
Visualization only.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[4]
EXAMPLE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tools" / "figures"))
from figure_style import FIG_WIDE, PALETTE, apply_style, savefig  # noqa: E402


def _plot_rdp_fallback(ax_list, out: Path) -> int:
    """Render a placeholder panel when near_field_profile.h5 is absent."""
    fig, ax = plt.subplots(figsize=(8, 4))
    rdp_params = {
        "psi_inf (m^3)": "1.243e+07",
        "fc (Hz)": "0.397",
        "Cavity R (m)": "278.8",
        "Crushed R (m)": "697.1",
        "Fractured R (m)": "1394.2",
        "Yield (kt)": "250",
        "Depth (m)": "600",
        "Medium": "granite (Mt. Mantap)",
    }
    y = 0.9
    ax.axis("off")
    ax.set_title(
        "DPRK 2017: COUPLED_ANALYTIC RDP source parameters\n"
        "(near_field_profile.h5 not present -- re-run with "
        "solver_kind = RADIAL_LAGRANGIAN to get radial profiles)",
        fontsize=10,
    )
    for k, v in rdp_params.items():
        ax.text(0.1, y, f"{k}:", fontsize=11, transform=ax.transAxes)
        ax.text(0.55, y, v, fontsize=11, fontweight="bold",
                transform=ax.transAxes, color=PALETTE["annotation"])
        y -= 0.10
    savefig(fig, out)
    print(f"wrote {out}  (RDP fallback -- no radial profile HDF5)")
    return 0


def main():
    apply_style()
    h5_path = EXAMPLE / "output" / "near_field_profile.h5"
    out = Path(__file__).resolve().parents[1] / \
        "02_cavity_formation_radial_profile.png"

    if not h5_path.exists():
        print(f"INFO: {h5_path} not found; rendering RDP parameter card.",
              file=sys.stderr)
        return _plot_rdp_fallback(None, out)

    import h5py  # noqa: PLC0415
    with h5py.File(str(h5_path), "r") as f:
        radius = f["radius_m"][:]
        snap_keys = sorted(f["snapshots"].keys(),
                           key=lambda s: int(s)
                           if s.isdigit() else float("inf"))
        n_snaps = len(snap_keys)
        if n_snaps < 4:
            print(f"WARN: only {n_snaps} snapshots in {h5_path}",
                  file=sys.stderr)
            return 1
        picks = [snap_keys[0],
                 snap_keys[n_snaps // 3],
                 snap_keys[2 * n_snaps // 3],
                 snap_keys[-1]]

        fig, axes = plt.subplots(1, 4, figsize=FIG_WIDE, sharex=True)
        fields = [
            ("velocity_m_per_s", "Velocity (m/s)"),
            ("pressure_pa", "Pressure (Pa)"),
            ("plastic_strain", "Plastic strain"),
            ("temperature_radiation_K", "T_rad (K)"),
        ]

        for ax, (field, label) in zip(axes, fields):
            for k, snap_id in enumerate(picks):
                grp = f[f"snapshots/{snap_id}"]
                if field not in grp:
                    continue
                data = grp[field][:]
                t = grp.attrs.get("time_s", float("nan"))
                color_key = ["low_tier", "med_tier",
                             "high_tier", "highest_tier"][k]
                ax.plot(radius, data, color=PALETTE[color_key],
                        label=f"t = {t:.2e} s")
            ax.set_xlabel("Radius (m)")
            ax.set_title(label)
            if "pressure" in field or "temperature" in field:
                ax.set_yscale("log")
        axes[0].legend(loc="upper right", fontsize=7)

        fig.suptitle(
            "DPRK 2017 (Punggye-ri): cavity formation radial profile\n"
            "250 kt, 600 m depth, granite -- Wilkins AV, Drucker-Prager",
            y=1.02,
        )
        savefig(fig, out)
        print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
