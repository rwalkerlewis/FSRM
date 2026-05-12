#!/usr/bin/env python3
"""03_moment_tensor_history.py -- DPRK 2017 (Punggye-ri, Sixth Test)

M(t) and Mdot(t) for the six independent moment-tensor components.
Reads near_field_history.csv written by the COUPLED_ANALYTIC or
RADIAL_LAGRANGIAN solver. Also annotates the RDP corner frequency
fc = 0.397 Hz and the Mueller-Murphy psi_inf.
Visualization only.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[4]
EXAMPLE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tools" / "figures"))
from figure_style import FIG_WIDE, PALETTE, apply_style, savefig  # noqa: E402


def main():
    apply_style()
    csv_path = EXAMPLE / "output" / "near_field_history.csv"
    out = Path(__file__).resolve().parents[1] / "03_moment_tensor_history.png"

    if not csv_path.exists():
        # Render a placeholder with the analytic RDP prediction.
        fc = 0.397          # Hz, from COUPLED_ANALYTIC run output
        psi_inf = 1.243e7   # m^3
        # Mueller-Murphy (1971): g(t) = omega^2 * t * exp(-omega*t)
        omega = 2 * np.pi * fc
        t = np.linspace(0, 5, 2000)
        g = omega**2 * t * np.exp(-omega * t)
        # Mxx = Myy = Mzz = (lambda + 2mu/3) * psi_inf * g(t); approximate
        # isotropic: M_iso(t) = 3 * bulk_modulus * psi_inf * g(t)
        # Granite bulk modulus ~ 50 GPa
        K = 50e9
        M_iso = 3 * K * psi_inf * g

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=FIG_WIDE)
        for comp, lw, ls in [("xx", 2.0, "-"), ("yy", 1.5, "--"),
                               ("zz", 1.5, ":")]:
            ax1.plot(t, M_iso, lw=lw, ls=ls, color=PALETTE["synthetic"],
                     label=f"M_{comp} (isotropic approx)")
        ax1.set_xlabel("t (s)")
        ax1.set_ylabel("M (N·m)")
        ax1.set_title("Moment tensor M(t) -- RDP analytic estimate")
        ax1.legend(fontsize=7)

        Mdot = omega**2 * np.exp(-omega * t) * (1 - omega * t)
        Mdot_iso = 3 * K * psi_inf * Mdot
        ax2.plot(t, Mdot_iso, color=PALETTE["synthetic"],
                 label="Mdot_xx (isotropic approx)")
        ax2.axvline(1.0 / fc, color=PALETTE["annotation"],
                    linestyle="--", label=f"1/fc = {1/fc:.2f} s")
        ax2.set_xlabel("t (s)")
        ax2.set_ylabel("Mdot (N·m/s)")
        ax2.set_title("Moment-rate Mdot(t) -- RDP analytic estimate")
        ax2.legend(fontsize=7)

        fig.suptitle(
            "DPRK 2017 (Punggye-ri): moment tensor history\n"
            f"fc={fc} Hz, psi_inf={psi_inf:.3e} m^3 (COUPLED_ANALYTIC)\n"
            "(near_field_history.csv not found; showing analytic RDP preview)",
            y=1.02,
        )
        savefig(fig, out)
        print(f"wrote {out}  (analytic RDP preview -- no CSV found)")
        return 0

    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    if not rows:
        print(f"WARN: {csv_path} empty", file=sys.stderr)
        return 1

    t = np.array([float(r["time_s"]) for r in rows])
    components = ["xx", "yy", "zz", "xy", "xz", "yz"]
    M_keys = [f"M_{c}_Nm" for c in components]
    Mdot_keys = [f"Mdot_{c}_Nms" for c in components]

    M = {}
    Mdot = {}
    for c, mk, dk in zip(components, M_keys, Mdot_keys):
        if mk in rows[0]:
            M[c] = np.array([float(r[mk]) for r in rows])
        if dk in rows[0]:
            Mdot[c] = np.array([float(r[dk]) for r in rows])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=FIG_WIDE)
    color_keys = ["synthetic", "observed", "analytic",
                  "low_tier", "med_tier", "high_tier"]
    for c, key in zip(components, color_keys):
        if c in M:
            ax1.plot(t, M[c], color=PALETTE[key], label=f"M_{c}")
        if c in Mdot:
            ax2.plot(t, Mdot[c], color=PALETTE[key], label=f"Mdot_{c}")

    ax1.set_xlabel("t (s)")
    ax1.set_ylabel("M (N·m)")
    ax1.set_title("Moment tensor M(t)")
    ax1.legend(fontsize=7, ncol=2)

    ax2.axvline(1.0 / 0.397, color=PALETTE["annotation"],
                linestyle="--", linewidth=0.8, label="1/fc = 2.52 s")
    ax2.set_xlabel("t (s)")
    ax2.set_ylabel("Mdot (N·m/s)")
    ax2.set_title("Moment-rate Mdot(t)")
    ax2.legend(fontsize=7, ncol=2)

    fig.suptitle(
        "DPRK 2017 (Punggye-ri, 6th test): 6-component moment tensor history\n"
        "250 kt @ 600 m depth in granite, fc=0.397 Hz",
        y=1.02,
    )
    savefig(fig, out)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
