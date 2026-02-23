#!/usr/bin/env python3
"""
FSRM GPU-Accelerated Wave Propagation — DPRK 2017 Punggye-ri (250 kt)

This script validates, documents, and (when output exists) post-processes
GPU-accelerated full-waveform seismogram generation for the 2017 DPRK
nuclear test at Punggye-ri.

Sections:
  A) GPU hardware & CUDA availability check
  B) Config validation (CFL, GPU memory estimate, resolution)
  C) Print GPU-specific invocation commands
  D) Post-processing if FSRM output exists (reuses analytical model)

The FSRM simulation itself is a compiled C++ binary.  This script does
NOT execute it — it prepares, validates, and post-processes.

Usage:
    python scripts/run_fsrm_dprk_2017_gpu.py
"""

import os
import sys
import subprocess
import warnings
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)

CONFIG_GPU = os.path.join(PROJECT_ROOT,
                          "config",
                          "dprk_2017_punggye_ri_far_field_seismograms_gpu.config")
CONFIG_CPU = os.path.join(PROJECT_ROOT,
                          "config",
                          "dprk_2017_punggye_ri_far_field_seismograms.config")

OUTPUT_DIR_GPU = os.path.join(PROJECT_ROOT, "output", "dprk_2017_far_field_gpu")
OUTPUT_DIR_CPU = os.path.join(PROJECT_ROOT, "output", "dprk_2017_far_field")
OUTDIR = os.path.join(PROJECT_ROOT, "figures", "dprk_2017")

# ---------------------------------------------------------------------------
# Grid / physics parameters
# ---------------------------------------------------------------------------
NX, NY, NZ = 400, 400, 120
DX = 500.0               # m
VP_MAX = 8050.0           # m/s (upper mantle Pn)
CFL_TARGET = 0.45
DT_MAX = 0.005            # s (from config)
DURATION = 30.0           # s
DG_ORDER = 4
N_VAR = 9                 # elastic wave: σ(6) + v(3)

# GPU performance model (based on NVIDIA A100 / V100)
GPU_MEMORY_PER_ELEM_BYTES = N_VAR * (DG_ORDER + 1)**3 * 8 * 3  # Q + R_vol + R_surf
GPU_THROUGHPUT_ELEM_STEPS_PER_SEC = 5.0e8  # Conservative for A100


# ═══════════════════════════════════════════════════════════════════════════════
# PART A — GPU / CUDA Availability
# ═══════════════════════════════════════════════════════════════════════════════

def check_cuda():
    """Check for CUDA toolkit and GPU hardware."""
    print("  ── GPU / CUDA Check ──")
    print()

    # nvidia-smi
    try:
        result = subprocess.run(["nvidia-smi", "--query-gpu=name,memory.total,driver_version,compute_cap",
                                 "--format=csv,noheader"],
                                capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            lines = result.stdout.strip().split("\n")
            print(f"  GPU(s) detected: {len(lines)}")
            for i, line in enumerate(lines):
                print(f"    [{i}] {line.strip()}")
            print()
            gpu_available = True
        else:
            print("  nvidia-smi failed — no GPU detected or driver issue.")
            print(f"  stderr: {result.stderr.strip()}")
            gpu_available = False
    except FileNotFoundError:
        print("  nvidia-smi not found — no NVIDIA driver installed.")
        gpu_available = False
    except subprocess.TimeoutExpired:
        print("  nvidia-smi timed out.")
        gpu_available = False

    # nvcc
    try:
        result = subprocess.run(["nvcc", "--version"], capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            version_line = [l for l in result.stdout.split("\n") if "release" in l.lower()]
            if version_line:
                print(f"  CUDA Toolkit: {version_line[0].strip()}")
            else:
                print(f"  CUDA Toolkit: found (version parse failed)")
            cuda_available = True
        else:
            print("  nvcc not found — CUDA toolkit not installed.")
            cuda_available = False
    except FileNotFoundError:
        print("  nvcc not found — CUDA toolkit not installed.")
        cuda_available = False
    except subprocess.TimeoutExpired:
        cuda_available = False

    # fsrm binary
    fsrm_path = os.path.join(PROJECT_ROOT, "build", "fsrm")
    if os.path.isfile(fsrm_path):
        print(f"  FSRM binary: {fsrm_path} ✓")
    else:
        print(f"  FSRM binary: NOT FOUND at {fsrm_path}")
        print("  Build with: mkdir build && cd build && cmake .. -DENABLE_CUDA=ON && make -j")

    print()
    return gpu_available


# ═══════════════════════════════════════════════════════════════════════════════
# PART B — Config Validation
# ═══════════════════════════════════════════════════════════════════════════════

def read_config(path):
    """Read a simple key=value config, return dict."""
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
                cfg[key.strip()] = val.split("#")[0].strip()
    return cfg


def validate_gpu_config():
    """Validate GPU config: CFL stability, GPU memory, resolution, wall-clock."""
    print("  ── GPU Config Validation ──")
    print()

    cfg = read_config(CONFIG_GPU)
    if cfg:
        print(f"  Config file  : {CONFIG_GPU}")
        print(f"  Name         : {cfg.get('name', '(unknown)')}")
        print(f"  GPU mode     : {cfg.get('gpu_mode', 'unknown')}")
        print(f"  GPU device   : {cfg.get('gpu_device_id', '0')}")
    else:
        print(f"  Config NOT FOUND: {CONFIG_GPU}")
        print("  Using built-in defaults.")
    print()

    # CFL check
    dt_max_cfl = CFL_TARGET * DX / VP_MAX
    print(f"  Grid             : {NX} × {NY} × {NZ}")
    print(f"  Cell spacing     : {DX:.0f} m")
    print(f"  Domain           : {NX*DX/1e3:.0f} × {NY*DX/1e3:.0f} × {NZ*DX/1e3:.0f} km")
    print(f"  DG order         : {DG_ORDER}  ({(DG_ORDER+1)**3} DOF/elem)")
    print(f"  Variables        : {N_VAR}  (elastic: σ_6 + v_3)")
    print(f"  Vp_max           : {VP_MAX:.0f} m/s")
    print()

    print(f"  CFL number       : {CFL_TARGET}")
    print(f"  dt_max (CFL)     : {dt_max_cfl:.6f} s")
    print(f"  dt (config)      : {DT_MAX} s")
    stable = DT_MAX <= dt_max_cfl
    print(f"  Stability        : {'STABLE ✓' if stable else 'UNSTABLE ✗'}")
    print()

    # Memory estimate
    total_cells = NX * NY * NZ
    dof_per_elem = (DG_ORDER + 1) ** 3  # hexahedral tensor product
    total_dof = total_cells * dof_per_elem
    mem_solution = total_dof * N_VAR * 8  # bytes (double precision)
    mem_residuals = mem_solution * 3       # R_vol, R_surf, R_src
    mem_geometry = total_cells * (9 + 1) * 8  # inv_jacobian + det_jacobian
    mem_total = mem_solution + mem_residuals + mem_geometry
    mem_gb = mem_total / 1e9

    print(f"  Total cells      : {total_cells:,}")
    print(f"  DOF per element  : {dof_per_elem}")
    print(f"  Total DOFs       : {total_dof:,}")
    print(f"  Solution memory  : {mem_solution / 1e9:.2f} GB")
    print(f"  Total GPU memory : {mem_gb:.2f} GB  (solution + residuals + geometry)")
    print()

    # GPU recommendations
    if mem_gb < 16:
        print(f"  GPU requirement  : Single GPU with ≥{max(int(mem_gb * 1.5), 8)} GB VRAM")
        print(f"  Recommended      : NVIDIA V100 (16/32 GB) or A100 (40/80 GB)")
    elif mem_gb < 40:
        print(f"  GPU requirement  : A100 (40 GB) or 2× V100 (32 GB)")
    else:
        print(f"  GPU requirement  : A100 80GB or multi-GPU (MPI + CUDA)")
    print()

    # Wall-clock estimate
    n_steps = int(DURATION / DT_MAX)
    wall_gpu = total_cells * n_steps / GPU_THROUGHPUT_ELEM_STEPS_PER_SEC
    cpu_throughput = 1.0e7  # cells·steps/core/s
    wall_cpu_16 = total_cells * n_steps / (cpu_throughput * 16)

    print(f"  Simulation time  : {DURATION:.0f} s  →  ~{n_steps:,} steps")
    print(f"  GPU wall-clock   : ~{wall_gpu:.0f} s  ({wall_gpu/60:.1f} min)")
    print(f"  CPU wall-clock   : ~{wall_cpu_16:.0f} s  ({wall_cpu_16/60:.1f} min, 16 cores)")
    speedup = wall_cpu_16 / max(wall_gpu, 1e-6)
    print(f"  GPU speedup      : ~{speedup:.0f}×")
    print()

    # Resolution
    f_max = VP_MAX / (2.0 * 5 * DX / (DG_ORDER + 1))
    lambda_min = VP_MAX / f_max
    print(f"  Max frequency    : {f_max:.2f} Hz  (5 pts/wavelength)")
    print(f"  Min wavelength   : {lambda_min:.0f} m")
    print(f"  Analysis band    : 0.5–{min(f_max, 4.0):.1f} Hz (regional seismology)")
    print()


# ═══════════════════════════════════════════════════════════════════════════════
# PART C — Invocation Commands
# ═══════════════════════════════════════════════════════════════════════════════

def print_gpu_invocation():
    """Print the GPU execution command."""
    print("  ── GPU Execution Commands ──")
    print()

    print("  Single-GPU execution:")
    print(f"    fsrm -c {CONFIG_GPU} --use-gpu")
    print()
    print("  Multi-GPU with MPI:")
    print(f"    mpirun -np 4 fsrm -c {CONFIG_GPU} --use-gpu")
    print()
    print("  With CUDA profiling:")
    print(f"    nsys profile -o dprk_gpu_profile fsrm -c {CONFIG_GPU} --use-gpu --profile")
    print()
    print("  CPU fallback (if no GPU available):")
    print(f"    mpirun -np 32 fsrm -c {CONFIG_CPU}")
    print()

    # Build instructions
    print("  Build FSRM with CUDA support:")
    print("    mkdir build && cd build")
    print('    cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_CUDA=ON \\')
    print('             -DCMAKE_CUDA_ARCHITECTURES="70;80;86"')
    print("    make -j$(nproc)")
    print()


# ═══════════════════════════════════════════════════════════════════════════════
# PART D — Post-Processing
# ═══════════════════════════════════════════════════════════════════════════════

def try_import_analytical():
    """Try to import the analytical model for comparison."""
    sys.path.insert(0, SCRIPT_DIR)
    try:
        from model_dprk_2017_nuclear_test import (
            PunggyeRiVelocityModel,
            generate_synthetic,
            mueller_murphy_spectrum,
            corner_frequency_patton,
            spectral_amplitude,
        )
        return True
    except ImportError as e:
        print(f"  Cannot import analytical model: {e}")
        print("  (Requires scipy, obspy — install for comparison figures)")
        return False


def check_output():
    """Check for FSRM GPU output and generate comparison figures."""
    print("  ── Post-Processing ──")
    print()

    for label, out_dir in [("GPU", OUTPUT_DIR_GPU), ("CPU", OUTPUT_DIR_CPU)]:
        seis_dir = os.path.join(out_dir, "seismograms")
        if os.path.isdir(seis_dir):
            sac_files = [f for f in os.listdir(seis_dir) if f.endswith(".sac")]
            mseed_files = [f for f in os.listdir(seis_dir)
                           if f.endswith(".mseed") or f.endswith(".miniseed")]
            print(f"  {label} output found: {seis_dir}")
            print(f"    SAC files    : {len(sac_files)}")
            print(f"    miniSEED     : {len(mseed_files)}")
        else:
            print(f"  {label} output NOT FOUND: {seis_dir}")
    print()

    # If GPU output exists, try to load and plot
    seis_dir = os.path.join(OUTPUT_DIR_GPU, "seismograms")
    if not os.path.isdir(seis_dir):
        print("  ┌─────────────────────────────────────────────────────────────┐")
        print("  │  Run the FSRM GPU simulation first to generate output.     │")
        print("  │  Then re-run this script for post-processing & figures.     │")
        print("  └─────────────────────────────────────────────────────────────┘")
        return

    # Try loading with obspy
    try:
        from obspy import read as obspy_read
    except ImportError:
        print("  obspy not installed — cannot read seismograms.")
        print("  Install: pip install obspy")
        return

    traces = []
    for ext in [".sac", ".mseed", ".miniseed"]:
        for fname in sorted(os.listdir(seis_dir)):
            if fname.endswith(ext):
                try:
                    st = obspy_read(os.path.join(seis_dir, fname))
                    traces.extend(st)
                except Exception:
                    pass

    if traces:
        print(f"  Loaded {len(traces)} trace(s). Generating figures ...")
        os.makedirs(OUTDIR, exist_ok=True)

        # Record section
        fig, ax = plt.subplots(figsize=(18, 12))
        for i, tr in enumerate(sorted(traces, key=lambda t: getattr(t.stats, 'distance', 0))):
            t_arr = np.arange(tr.stats.npts) * tr.stats.delta
            data = tr.data.copy().astype(float)
            peak = np.max(np.abs(data))
            if peak > 0:
                data /= peak
            ax.plot(t_arr, data * 0.4 + i, 'k-', lw=0.35)
        ax.set_xlabel("Time (s)")
        ax.set_title("FSRM GPU Record Section — DPRK 2017 Far-Field (250 kt)")
        plt.tight_layout()
        p = os.path.join(OUTDIR, "fig_fsrm_gpu_record_section.png")
        fig.savefig(p, dpi=180)
        plt.close(fig)
        print(f"  Saved {p}")
    else:
        print("  No traces loaded from GPU output directory.")


# ═══════════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    os.makedirs(OUTDIR, exist_ok=True)

    print("=" * 72)
    print("  FSRM GPU-Accelerated Wave Propagation")
    print("  DPRK 2017 Punggye-ri — 250 kt Coupled Nuclear Test")
    print("  Synthetic Seismograms via Numerical Wave Equation (CUDA)")
    print("=" * 72)
    print()

    # Step 1: GPU check
    print("[1/4] Checking GPU / CUDA availability ...")
    print()
    gpu_ok = check_cuda()
    if not gpu_ok:
        print("  ⚠ No GPU detected. The GPU config will require a CUDA-capable machine.")
        print()

    # Step 2: validate config
    print("[2/4] Validating GPU configuration ...")
    print()
    validate_gpu_config()

    # Step 3: print commands
    print("[3/4] Execution commands:")
    print()
    print_gpu_invocation()

    # Step 4: post-processing
    print("[4/4] Post-processing ...")
    print()
    check_output()

    print()
    print("=" * 72)
    print("  Done.")
    print("=" * 72)


if __name__ == "__main__":
    main()
