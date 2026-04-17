#!/usr/bin/env python3
"""
Set Up, Validate, and Document FSRM Far-Field Wave Propagation — DPRK 2017

This script does NOT execute the FSRM C++ numerical simulation.  Instead it:
  A) Validates the FSRM config: CFL stability, memory, wall-clock estimate
  B) Prints the mpirun invocation command
  C) Provides post-processing functions for FSRM output (SAC / miniSEED)
  D) Generates comparison figures (numerical vs analytical) if output exists

The actual FSRM simulation requires the compiled C++ binary and HPC resources.

Usage:
    python scripts/run_fsrm_dprk_2017_far_field.py
"""

import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)

sys.path.insert(0, SCRIPT_DIR)

from fsrm.run_utils import (
    validate_config, print_invocation, postprocess_far_field,
)
from fsrm.source_physics import (
    mueller_murphy_spectrum, corner_frequency_patton,
)
from fsrm.signal_processing import spectral_amplitude
from fsrm.propagation import generate_synthetic
from fsrm.velocity_models import PunggyeRiVelocityModel

# ── Event / grid parameters ─────────────────────────────────────────────────
CONFIG_PATH = os.path.join(PROJECT_ROOT, "config",
                           "dprk_2017_punggye_ri_far_field_seismograms.config")
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "output",
                          "dprk_2017_far_field", "seismograms")
OUTDIR = os.path.join(PROJECT_ROOT, "figures", "dprk_2017")

YIELD_KT = 250.0

GRID = dict(
    nx=400, ny=400, nz=120,
    dx=500.0,
    vp_max=8050.0,
    domain_km=200.0,
    cfl_target=0.4,
    dt_suggested=0.02,
    nprocs_min=16, nprocs_max=64,
    duration=300.0,
)

EVENT_LABEL = "DPRK 2017, 250 kt coupled, Punggye-ri model"


def main():
    os.makedirs(OUTDIR, exist_ok=True)

    print("=" * 72)
    print("  FSRM Far-Field Wave Propagation — DPRK 2017 Punggye-ri (250 kt)")
    print("  Setup, Validation & Post-Processing")
    print("=" * 72)
    print()

    print("[1/3] Validating FSRM configuration ...")
    print()
    validate_config(CONFIG_PATH, **GRID)

    print("[2/3] FSRM invocation command:")
    print_invocation(CONFIG_PATH, nprocs=32,
                     nprocs_min=GRID["nprocs_min"],
                     nprocs_max=GRID["nprocs_max"])

    print("[3/3] Post-processing ...")
    print()

    model = PunggyeRiVelocityModel()
    fc = corner_frequency_patton(YIELD_KT)

    postprocess_far_field(
        OUTPUT_DIR, OUTDIR,
        model=model,
        generate_synth_fn=generate_synthetic,
        spectral_fn=spectral_amplitude,
        mm_spectrum_fn=mueller_murphy_spectrum,
        yield_kt=YIELD_KT,
        fc=fc,
        event_label=EVENT_LABEL,
    )

    print()
    print("=" * 72)
    print("  Done.")
    print("=" * 72)


if __name__ == "__main__":
    main()
