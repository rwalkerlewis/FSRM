#!/usr/bin/env python3
"""
Set Up, Validate, and Document FSRM Far-Field Wave Propagation — US Divider 1992

This script does NOT execute the FSRM C++ numerical simulation.  Instead it:
  A) Validates the FSRM config: CFL stability, memory, wall-clock estimate
  B) Prints the mpirun invocation command
  C) Provides post-processing functions for FSRM output (SAC / miniSEED)
  D) Generates comparison figures (numerical vs analytical) if output exists

The actual FSRM simulation requires the compiled C++ binary and HPC resources.

Usage:
    python scripts/run_fsrm_divider_1992_far_field.py
"""

import os
import sys
import warnings
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import AutoMinorLocator

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)

CONFIG_PATH = os.path.join(PROJECT_ROOT,
                           "config",
                           "us_divider_1992_far_field_seismograms.config")
OUTPUT_DIR = os.path.join(PROJECT_ROOT,
                          "output", "us_divider_1992_far_field", "seismograms")
OUTDIR = os.path.join(PROJECT_ROOT, "figures", "us_divider_1992")

# ---------------------------------------------------------------------------
# Grid / physics parameters for validation
# ---------------------------------------------------------------------------
NX, NY, NZ = 300, 300, 80
DX = 500.0            # m
VP_MAX = 7800.0       # m/s  (upper-mantle Pn beneath Basin & Range)
DOMAIN_KM = 150.0     # km  (lateral extent)
CFL_TARGET = 0.4
DT_SUGGESTED = 0.02   # s
NPROCS_MIN, NPROCS_MAX = 8, 32
BYTES_PER_DOF = 8     # double precision
DOF_PER_CELL = 9      # stress(6) + velocity(3)
DURATION = 250.0      # s  simulation time
THROUGHPUT = 1.0e7    # cells·steps / core / second  (rough estimate)

# ---------------------------------------------------------------------------
# Import analytical model for comparison
# ---------------------------------------------------------------------------
sys.path.insert(0, SCRIPT_DIR)
from model_us_divider_1992_nuclear_test import (
    NTSVelocityModel,
    generate_synthetic,
    mueller_murphy_spectrum,
    corner_frequency_patton,
    scalar_moment_coupled,
    spectral_amplitude,
    regional_travel_times,
    YIELD_KT,
    CORNER_FREQ,
    M0_COUPLED,
)


# ═══════════════════════════════════════════════════════════════════════════════
# PART A — Config Validation
# ═══════════════════════════════════════════════════════════════════════════════

def read_config(path):
    """Read a simple key=value config file, returning a dict of strings."""
    cfg = {}
    if not os.path.isfile(path):
        return cfg
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("["):
                continue
            if "=" in line:
                key, _, val = line.partition("=")
                key = key.strip()
                val = val.split("#")[0].strip()
                cfg[key] = val
    return cfg


def validate_config():
    """Compute and print CFL stability, memory, wall-clock, and resolution."""
    print("  ── Config Validation ──")
    print()

    cfg = read_config(CONFIG_PATH)
    if cfg:
        print(f"  Config file found : {CONFIG_PATH}")
        print(f"  Simulation name   : {cfg.get('name', '(unknown)')}")
    else:
        print(f"  Config file NOT found: {CONFIG_PATH}")
        print("  Using built-in parameter defaults for validation.")
    print()

    dt_max_cfl = CFL_TARGET * DX / VP_MAX
    print(f"  Grid dimensions   : {NX} x {NY} x {NZ}")
    print(f"  Grid spacing (dx) : {DX:.0f} m")
    print(f"  Domain extent     : {DOMAIN_KM:.0f} km  ({NX*DX/1e3:.0f} x {NY*DX/1e3:.0f} x {NZ*DX/1e3:.0f} km)")
    print(f"  Vp_max            : {VP_MAX:.0f} m/s")
    print()

    print(f"  CFL number        : {CFL_TARGET}")
    print(f"  dt_max (CFL)      : {dt_max_cfl:.6f} s")
    print(f"  dt (suggested)    : {DT_SUGGESTED} s")
    stable = DT_SUGGESTED <= dt_max_cfl
    tag = "STABLE" if stable else "UNSTABLE"
    print(f"  Stability check   : {tag}  (dt {'<=' if stable else '>'} dt_max)")
    print()

    total_cells = NX * NY * NZ
    mem_bytes = total_cells * DOF_PER_CELL * BYTES_PER_DOF
    mem_gb = mem_bytes / 1e9
    print(f"  Total grid cells  : {total_cells:,}")
    print(f"  DOF per cell      : {DOF_PER_CELL}  ({DOF_PER_CELL} x {BYTES_PER_DOF} B)")
    print(f"  Memory estimate   : {mem_gb:.2f} GB  (field arrays only)")
    print(f"  With ghost/buffers: ~{mem_gb * 2.5:.1f} GB  (×2.5 overhead)")
    print()

    n_steps = int(DURATION / DT_SUGGESTED)
    wall_serial = total_cells * n_steps / THROUGHPUT
    wall_min = wall_serial / NPROCS_MIN
    wall_max = wall_serial / NPROCS_MAX
    print(f"  Simulation time   : {DURATION:.0f} s  →  ~{n_steps:,} time-steps")
    print(f"  Throughput est.   : {THROUGHPUT:.0e} cell·steps/core/s")
    print(f"  Wall-clock ({NPROCS_MIN:2d} cores): ~{wall_min/3600:.1f} h")
    print(f"  Wall-clock ({NPROCS_MAX:2d} cores): ~{wall_max/3600:.1f} h")
    print()

    f_max_resolved = VP_MAX / (2.0 * 5 * DX)
    wavelength_min = VP_MAX / f_max_resolved
    print(f"  Frequency limit   : f_max ≈ {f_max_resolved:.2f} Hz  (5 pts/wavelength)")
    print(f"  Min wavelength    : λ_min ≈ {wavelength_min:.0f} m")
    print()


# ═══════════════════════════════════════════════════════════════════════════════
# PART B — Print MPI Invocation
# ═══════════════════════════════════════════════════════════════════════════════

def print_invocation(nprocs=16):
    """Print the mpirun command without executing it."""
    config_path = CONFIG_PATH
    cmd = f"mpirun -np {nprocs} fsrm -c {config_path}"
    print(f"\n  To run the FSRM numerical simulation:\n    {cmd}\n")
    print("  NOTE: This simulation is NOT executed by this script.")
    print("  It requires the FSRM C++ binary and sufficient compute resources.")
    print(f"  Recommended core count: {NPROCS_MIN}–{NPROCS_MAX}")
    print()


# ═══════════════════════════════════════════════════════════════════════════════
# PART C — Post-Processing Functions
# ═══════════════════════════════════════════════════════════════════════════════

def read_sac_seismograms(output_dir):
    """
    Read SAC-format seismograms from FSRM output directory.

    Returns a list of obspy.Trace objects, one per file found.
    """
    from obspy import read as obspy_read

    traces = []
    if not os.path.isdir(output_dir):
        print(f"  [read_sac] Directory does not exist: {output_dir}")
        return traces

    sac_files = sorted(f for f in os.listdir(output_dir) if f.endswith(".sac"))
    if not sac_files:
        print(f"  [read_sac] No .sac files in {output_dir}")
        return traces

    for fname in sac_files:
        path = os.path.join(output_dir, fname)
        try:
            st = obspy_read(path, format="SAC")
            for tr in st:
                tr.stats._fsrm_file = fname
                traces.append(tr)
        except Exception as exc:
            print(f"  [read_sac] Could not read {fname}: {exc}")
    print(f"  [read_sac] Loaded {len(traces)} trace(s) from {len(sac_files)} SAC file(s)")
    return traces


def read_mseed_seismograms(output_dir):
    """
    Read miniSEED-format seismograms from FSRM output directory.

    Returns a list of obspy.Trace objects.
    """
    from obspy import read as obspy_read

    traces = []
    if not os.path.isdir(output_dir):
        print(f"  [read_mseed] Directory does not exist: {output_dir}")
        return traces

    mseed_files = sorted(
        f for f in os.listdir(output_dir)
        if f.endswith(".mseed") or f.endswith(".miniseed")
    )
    if not mseed_files:
        print(f"  [read_mseed] No miniSEED files in {output_dir}")
        return traces

    for fname in mseed_files:
        path = os.path.join(output_dir, fname)
        try:
            st = obspy_read(path, format="MSEED")
            for tr in st:
                tr.stats._fsrm_file = fname
                traces.append(tr)
        except Exception as exc:
            print(f"  [read_mseed] Could not read {fname}: {exc}")
    print(f"  [read_mseed] Loaded {len(traces)} trace(s) from {len(mseed_files)} file(s)")
    return traces


def plot_numerical_record_section(traces, title_extra=""):
    """
    Plot a record section from FSRM numerical seismograms.

    Parameters
    ----------
    traces : list of obspy.Trace
        Seismogram traces with distance stored in ``tr.stats.sac.dist``
        or ``tr.stats.distance`` (km).

    Returns the matplotlib Figure.
    """
    if not traces:
        print("  [record_section] No traces to plot.")
        return None

    fig, ax = plt.subplots(figsize=(18, 12))

    dists = []
    for tr in traces:
        d = getattr(tr.stats.sac, "dist", None) if hasattr(tr.stats, "sac") else None
        if d is None:
            d = getattr(tr.stats, "distance", 0.0)
        dists.append(d)

    order = np.argsort(dists)
    y_positions = np.arange(len(traces))

    for yi, idx in enumerate(order):
        tr = traces[idx]
        dist_km = dists[idx]
        t_arr = np.arange(tr.stats.npts) * tr.stats.delta
        data = tr.data.copy().astype(float)
        peak = np.max(np.abs(data))
        if peak > 0:
            data /= peak
        data *= 0.4

        ax.plot(t_arr, data + yi, "k-", lw=0.35, alpha=0.85)
        ax.fill_between(t_arr, data + yi, yi, where=(data > 0),
                        color="#444", alpha=0.10)
        sta_label = getattr(tr.stats, "station", "?")
        ax.text(-2, yi, f"{sta_label} ({dist_km:.0f} km)",
                fontsize=7, ha="right", va="center", fontweight="bold",
                clip_on=False)

    ax.set_xlabel("Time (s)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Station (sorted by distance)", fontsize=11, fontweight="bold")
    ax.set_yticks(y_positions)
    ax.set_yticklabels([])
    ax.set_ylim(-0.8, len(traces) - 0.2)
    ax.grid(True, axis="x", alpha=0.2)
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.set_title(f"FSRM Numerical Record Section — US Divider 1992 Far-Field{title_extra}",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    return fig


def compare_numerical_vs_analytical(numerical_traces, analytical_model,
                                    yield_kt=YIELD_KT, dt_synth=0.05,
                                    duration=500.0):
    """
    Compare FSRM PDE numerical seismograms with analytical Mueller-Murphy
    synthetics at matching distances.

    Returns a matplotlib Figure with waveform + spectral comparison panels.
    """
    if not numerical_traces:
        print("  [compare] No numerical traces.")
        return None

    num_dists = []
    for tr in numerical_traces:
        d = getattr(tr.stats.sac, "dist", None) if hasattr(tr.stats, "sac") else None
        if d is None:
            d = getattr(tr.stats, "distance", 0.0)
        num_dists.append(d)

    n = min(len(numerical_traces), 8)
    order = np.argsort(num_dists)[:n]

    fig, axes = plt.subplots(n, 2, figsize=(20, 3.5 * n), squeeze=False)

    model = analytical_model

    for row, idx in enumerate(order):
        tr = numerical_traces[idx]
        dist_km = num_dists[idx]

        t_num = np.arange(tr.stats.npts) * tr.stats.delta
        d_num = tr.data.copy().astype(float)
        peak_num = np.max(np.abs(d_num))
        if peak_num > 0:
            d_num /= peak_num

        t_ana, v_ana = generate_synthetic(dist_km, yield_kt,
                                          dt=dt_synth, duration=duration,
                                          model=model)
        peak_ana = np.max(np.abs(v_ana))
        if peak_ana > 0:
            v_ana /= peak_ana

        ax = axes[row, 0]
        ax.plot(t_num, d_num, "k-", lw=0.5, alpha=0.7, label="FSRM numerical")
        ax.plot(t_ana, v_ana, "r-", lw=0.6, alpha=0.6, label="Mueller-Murphy analytic")
        ax.set_xlim(0, min(t_num[-1], duration))
        ax.set_ylim(-1.3, 1.3)
        sta = getattr(tr.stats, "station", "?")
        ax.set_title(f"{sta} ({dist_km:.0f} km) — Waveforms",
                     fontsize=9, fontweight="bold")
        ax.legend(fontsize=7, loc="upper right")
        ax.grid(True, alpha=0.2)
        if row == n - 1:
            ax.set_xlabel("Time (s)")

        ax2 = axes[row, 1]
        freq_n, sp_n = spectral_amplitude(d_num, tr.stats.delta)
        freq_a, sp_a = spectral_amplitude(v_ana, dt_synth)
        m_n = (freq_n > 0.1) & (freq_n < 10)
        m_a = (freq_a > 0.1) & (freq_a < 10)
        if np.any(m_n) and np.max(sp_n[m_n]) > 0:
            ax2.loglog(freq_n[m_n], sp_n[m_n] / np.max(sp_n[m_n]),
                       "k-", lw=1, alpha=0.7, label="Numerical")
        if np.any(m_a) and np.max(sp_a[m_a]) > 0:
            ax2.loglog(freq_a[m_a], sp_a[m_a] / np.max(sp_a[m_a]),
                       "r-", lw=1, alpha=0.7, label="Analytical")
        ax2.set_title(f"{sta} — Spectra", fontsize=9, fontweight="bold")
        ax2.legend(fontsize=7)
        ax2.grid(True, alpha=0.2, which="both")
        ax2.set_xlim(0.1, 10)
        if row == n - 1:
            ax2.set_xlabel("Frequency (Hz)")

    fig.suptitle("FSRM Numerical vs Analytical (Mueller-Murphy) — US Divider 1992\n"
                 f"~{YIELD_KT} kt coupled, NTS tuff velocity model",
                 fontsize=13, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    return fig


def extrapolate_to_far_field(traces, true_distances_km, domain_distances_km,
                             Q_avg=400.0, f_ref=1.0):
    """
    Apply Q-corrected geometric spreading to extend domain-edge waveforms
    to true station distances.

    For each trace recorded at the domain boundary (``domain_distances_km``),
    scale amplitude to the actual station distance (``true_distances_km``)
    using:
        A_far = A_dom × (r_dom / r_true) × exp(-π f_ref Δt* / Q)
    where Δt* is the additional travel time beyond the domain edge.

    Parameters
    ----------
    traces : list of obspy.Trace
    true_distances_km : array-like, true epicentral distance per trace (km)
    domain_distances_km : array-like, distance at domain edge per trace (km)
    Q_avg : float, average path Q for attenuation correction
    f_ref : float, reference frequency (Hz) for t* correction

    Returns list of corrected obspy.Trace copies.
    """
    corrected = []
    Vp_avg = 5.5  # km/s average crustal P for Basin & Range

    for i, tr in enumerate(traces):
        r_dom = max(domain_distances_km[i], 1.0)
        r_true = max(true_distances_km[i], r_dom + 0.1)

        geom_factor = r_dom / r_true

        delta_t = (r_true - r_dom) / Vp_avg
        atten_factor = np.exp(-np.pi * f_ref * delta_t / Q_avg)

        scale = geom_factor * atten_factor

        tr_out = tr.copy()
        tr_out.data = tr_out.data.astype(float) * scale
        if hasattr(tr_out.stats, "sac"):
            tr_out.stats.sac.dist = r_true
        tr_out.stats.distance = r_true
        corrected.append(tr_out)
        print(f"  [extrapolate] {getattr(tr.stats, 'station', '?'):>6s}: "
              f"{r_dom:.0f} → {r_true:.0f} km  "
              f"(geom={geom_factor:.4f}, atten={atten_factor:.4f}, "
              f"scale={scale:.5f})")

    return corrected


# ═══════════════════════════════════════════════════════════════════════════════
# PART D — Comparison Figure Generation
# ═══════════════════════════════════════════════════════════════════════════════

def _get_trace_distance(tr):
    d = getattr(tr.stats.sac, "dist", None) if hasattr(tr.stats, "sac") else None
    if d is None:
        d = getattr(tr.stats, "distance", 0.0)
    return d


def fig_record_section(traces):
    """Record section of FSRM numerical seismograms."""
    fig = plot_numerical_record_section(traces,
                                        title_extra=f" (NTS ~{YIELD_KT} kt)")
    if fig is None:
        return
    p = os.path.join(OUTDIR, "fig_fsrm_far_field_record_section.png")
    fig.savefig(p, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {p}")


def fig_numerical_vs_analytical(traces):
    """Numerical vs analytical waveform comparison at matching distances."""
    model = NTSVelocityModel()
    fig = compare_numerical_vs_analytical(traces, model,
                                          yield_kt=YIELD_KT,
                                          dt_synth=0.05,
                                          duration=500.0)
    if fig is None:
        return
    p = os.path.join(OUTDIR, "fig_fsrm_far_field_num_vs_analytical.png")
    fig.savefig(p, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {p}")


def fig_spectral_comparison(traces):
    """Spectral comparison between numerical and analytical for key distances."""
    model = NTSVelocityModel()
    yield_kt = YIELD_KT
    fc = CORNER_FREQ

    dists = [_get_trace_distance(tr) for tr in traces]
    order = np.argsort(dists)
    n = min(len(traces), 6)
    indices = order[:n]

    fig, axes = plt.subplots(2, 3, figsize=(20, 10))
    axes = axes.ravel()

    f_source = np.logspace(-2, 1.5, 500)
    spec_source = mueller_murphy_spectrum(f_source, yield_kt, fc=fc)
    spec_source_norm = spec_source / np.max(spec_source)

    for panel, idx in enumerate(indices):
        if panel >= 6:
            break
        tr = traces[idx]
        dist_km = dists[idx]
        ax = axes[panel]

        freq_n, sp_n = spectral_amplitude(tr.data.astype(float), tr.stats.delta)
        m_n = (freq_n > 0.05) & (freq_n < 12)
        if np.any(m_n) and np.max(sp_n[m_n]) > 0:
            ax.loglog(freq_n[m_n], sp_n[m_n] / np.max(sp_n[m_n]),
                      "k-", lw=1.2, label="FSRM numerical")

        t_ana, v_ana = generate_synthetic(dist_km, yield_kt, dt=0.05,
                                          duration=500.0, model=model)
        freq_a, sp_a = spectral_amplitude(v_ana, 0.05)
        m_a = (freq_a > 0.05) & (freq_a < 12)
        if np.any(m_a) and np.max(sp_a[m_a]) > 0:
            ax.loglog(freq_a[m_a], sp_a[m_a] / np.max(sp_a[m_a]),
                      "r--", lw=1, alpha=0.7, label="Analytical")

        ax.loglog(f_source, spec_source_norm, "b:", lw=0.8, alpha=0.4,
                  label="Source spectrum")
        ax.axvline(fc, color="green", ls=":", lw=0.8, alpha=0.5)

        sta = getattr(tr.stats, "station", "?")
        ax.set_title(f"{sta} ({dist_km:.0f} km)", fontsize=9, fontweight="bold")
        ax.set_xlim(0.05, 12)
        ax.legend(fontsize=6)
        ax.grid(True, alpha=0.2, which="both")
        if panel >= 3:
            ax.set_xlabel("Frequency (Hz)")
        if panel % 3 == 0:
            ax.set_ylabel("Normalised amplitude")

    for j in range(panel + 1, 6):
        axes[j].axis("off")

    fig.suptitle("Spectral Comparison — FSRM Numerical vs Analytical\n"
                 f"US Divider 1992, ~{YIELD_KT} kt coupled, NTS tuff model",
                 fontsize=13, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    p = os.path.join(OUTDIR, "fig_fsrm_far_field_spectral_comparison.png")
    fig.savefig(p, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {p}")


def fig_amplitude_decay(traces):
    """Amplitude decay curve: numerical peak amplitudes vs 1/r analytical."""
    dists = np.array([_get_trace_distance(tr) for tr in traces])
    peaks = np.array([np.max(np.abs(tr.data.astype(float))) for tr in traces])

    valid = (dists > 0) & (peaks > 0)
    if not np.any(valid):
        print("  [amplitude_decay] No valid data for plot.")
        return

    dists_v = dists[valid]
    peaks_v = peaks[valid]

    fig, ax = plt.subplots(figsize=(12, 7))

    ax.loglog(dists_v, peaks_v, "ko", ms=7, label="FSRM numerical peak amplitude")

    r_ana = np.linspace(dists_v.min() * 0.8, dists_v.max() * 1.2, 200)
    ref_idx = np.argmin(dists_v)
    A0 = peaks_v[ref_idx]
    r0 = dists_v[ref_idx]
    ana_1r = A0 * (r0 / r_ana)
    ax.loglog(r_ana, ana_1r, "r--", lw=2, label="1/r geometric spreading")

    ana_sqrt = A0 * np.sqrt(r0 / r_ana)
    ax.loglog(r_ana, ana_sqrt, "b:", lw=1.5, alpha=0.6,
              label=r"$1/\sqrt{r}$ (surface wave)")

    ax.set_xlabel("Distance (km)", fontsize=11)
    ax.set_ylabel("Peak amplitude (arb. units)", fontsize=11)
    ax.set_title("Amplitude Decay — FSRM Numerical vs Analytical Spreading\n"
                 f"US Divider 1992, ~{YIELD_KT} kt coupled",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.2, which="both")
    plt.tight_layout()

    p = os.path.join(OUTDIR, "fig_fsrm_far_field_amplitude_decay.png")
    fig.savefig(p, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {p}")


# ═══════════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    os.makedirs(OUTDIR, exist_ok=True)

    print("=" * 72)
    print(f"  FSRM Far-Field Wave Propagation — US Divider 1992 (~{YIELD_KT} kt)")
    print("  Setup, Validation & Post-Processing")
    print("=" * 72)
    print()

    # Step 1: config validation (always)
    print("[1/3] Validating FSRM configuration ...")
    print()
    validate_config()

    # Step 2: print invocation (always)
    print("[2/3] FSRM invocation command:")
    print_invocation(nprocs=16)

    # Step 3: post-processing (only if output exists)
    print("[3/3] Post-processing ...")
    print()

    if os.path.isdir(OUTPUT_DIR):
        print(f"  Output directory found: {OUTPUT_DIR}")
        print("  Attempting to read FSRM seismograms ...")
        print()

        traces = read_sac_seismograms(OUTPUT_DIR)
        if not traces:
            traces = read_mseed_seismograms(OUTPUT_DIR)

        if traces:
            print()
            print(f"  Loaded {len(traces)} numerical trace(s). Generating figures ...")
            print()
            fig_record_section(traces)
            fig_numerical_vs_analytical(traces)
            fig_spectral_comparison(traces)
            fig_amplitude_decay(traces)
            print()
            print(f"  Figures saved to: {OUTDIR}")
        else:
            print("  No SAC or miniSEED traces found in output directory.")
            print("  Ensure the FSRM simulation has completed and written seismograms.")
    else:
        print(f"  Output directory NOT found: {OUTPUT_DIR}")
        print()
        print("  ┌─────────────────────────────────────────────────────────────┐")
        print("  │  Run the FSRM simulation first to generate seismograms.    │")
        print("  │  Then re-run this script for post-processing & figures.    │")
        print("  └─────────────────────────────────────────────────────────────┘")

    print()
    print("=" * 72)
    print("  Done.")
    print("=" * 72)


if __name__ == "__main__":
    main()
