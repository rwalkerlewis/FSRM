#!/usr/bin/env python3
"""02_cavity_formation_radial_profile.py -- Cannikin 1971

4-panel radial-profile snapshot at four times during cavity formation.
Reads near_field_profile.h5 written by the radial Lagrangian solver.
Visualization only.
"""

from __future__ import annotations

import sys
from pathlib import Path

import h5py
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[4]
EXAMPLE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tools" / "figures"))
from figure_style import FIG_WIDE, PALETTE, apply_style, savefig  # noqa: E402


def main():
    apply_style()
    h5_path = EXAMPLE / "output" / "near_field_profile.h5"
    if not h5_path.exists():
        print(f"WARN: {h5_path} not found; run the simulation first.",
              file=sys.stderr)
        return 1

    with h5py.File(str(h5_path), "r") as f:
        snap_keys = sorted(f["profiles"].keys(), key=lambda s: int(s))
        n_snaps = len(snap_keys)
        if n_snaps < 4:
            print(f"WARN: only {n_snaps} snapshots in {h5_path}",
                  file=sys.stderr)
            return 1
        # Pick four spread snapshots
        picks = [snap_keys[0],
                 snap_keys[n_snaps // 3],
                 snap_keys[2 * n_snaps // 3],
                 snap_keys[-1]]
        times = f["time"][:]

        # Read cell-centre radii from first snapshot (constant grid)
        r_cell = f[f"profiles/{picks[0]}/r_cell"][:]

        fig, axes = plt.subplots(1, 4, figsize=FIG_WIDE, sharex=True)
        # Map (h5_field, use_node_r, label, log_scale)
        fields = [
            ("v_r",    True,  "Velocity (m/s)",    False),
            ("p",      False, "Pressure (Pa)",      True),
            ("eps_p",  False, "Plastic strain (--)", False),
            ("damage", False, "Damage (--)",         False),
        ]

        for ax, (field, use_node, label, do_log) in zip(axes, fields):
            for k, snap_id in enumerate(picks):
                grp = f[f"profiles/{snap_id}"]
                if field not in grp:
                    continue
                raw = grp[field][:]
                if use_node:
                    # v_r lives on nodes (201 pts); average to cell centres
                    data = 0.5 * (raw[:-1] + raw[1:])
                    r = r_cell
                else:
                    data = raw
                    r = r_cell
                t_snap = times[int(snap_id)]
                color_key = ["low_tier", "med_tier",
                             "high_tier", "highest_tier"][k]
                ax.plot(r, data, color=PALETTE[color_key],
                        label=f"t = {t_snap:.2e} s")
            ax.set_xlabel("Radius (m)")
            ax.set_ylabel(label)
            ax.set_title(label)
            if do_log:
                ax.set_yscale("log")
        axes[0].legend(loc="upper right", fontsize=7)

        fig.suptitle("Cannikin 1971: cavity formation radial profile",
                     y=1.02)
        out = Path(__file__).resolve().parents[1] / \
            "02_cavity_formation_radial_profile.png"
        savefig(fig, out)
        print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
