#!/bin/bash
# Pass-12 followup 2 (V&V hardening): smoke runner for the
# examples_runtime CTest gate.
#
# Invokes a single examples/<N>_<event>/run.sh under a smoke
# budget (FSRM_FINAL_TIME_OVERRIDE), then asserts that the
# example produced at least one non-empty file under its
# output/ directory. Exits 0 on success, non-zero on failure
# with a diagnostic message naming the example, the MPI rank
# count, and the assertion that failed.
#
# Usage:
#   tests/integration/run_example_smoke.sh <example_dir> <mpi_ranks> [<final_time>]
#
# Example:
#   tests/integration/run_example_smoke.sh examples/01_uniaxial_compression 4 0.01
#
# The wrapper cleans the example's output/ directory before
# invoking run.sh so a stale artifact from a previous run cannot
# satisfy the post-run assertions.

set -uo pipefail

if [ "$#" -lt 2 ]; then
    echo "[smoke] usage: $0 <example_dir> <mpi_ranks> [<final_time>]" >&2
    exit 64
fi

EXAMPLE_DIR="$1"
MPI_RANKS="$2"
FINAL_TIME="${3:-0.01}"

if [ ! -d "$EXAMPLE_DIR" ]; then
    echo "[smoke] FAIL: $EXAMPLE_DIR is not a directory" >&2
    exit 65
fi

RUN_SH="$EXAMPLE_DIR/run.sh"
if [ ! -x "$RUN_SH" ] && [ ! -f "$RUN_SH" ]; then
    echo "[smoke] FAIL: $RUN_SH is missing" >&2
    exit 66
fi

OUT_DIR="$EXAMPLE_DIR/output"
rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR"

# Prefer a config_ci.config when present. The smoke wrapper points
# the example's run.sh at the CI variant via FSRM_CONFIG_OVERRIDE,
# which run.sh consults when set. Falls back to the default
# config.config so examples without a CI variant still smoke-test.
if [ -f "$EXAMPLE_DIR/config_ci.config" ]; then
    export FSRM_CONFIG_OVERRIDE="$EXAMPLE_DIR/config_ci.config"
    echo "[smoke] using config_ci.config"
fi

echo "[smoke] $EXAMPLE_DIR ranks=$MPI_RANKS final_time=$FINAL_TIME"

export MPI_RANKS
export FSRM_FINAL_TIME_OVERRIDE="$FINAL_TIME"

# Run the example. Capture stdout and stderr. Any non-zero exit
# fails the gate; do not silence it.
bash "$RUN_SH"
RC=$?

if [ "$RC" -ne 0 ]; then
    echo "[smoke] FAIL: $EXAMPLE_DIR run.sh exited $RC at MPI_RANKS=$MPI_RANKS" >&2
    exit "$RC"
fi

# Post-run: assert at least one non-empty regular file under the
# example's output/ directory. Some examples write to a nested
# subdirectory; descend recursively.
NON_EMPTY_FILES=$(find "$OUT_DIR" -type f -size +0c 2>/dev/null | wc -l)
if [ "$NON_EMPTY_FILES" -le 0 ]; then
    echo "[smoke] FAIL: $EXAMPLE_DIR produced no non-empty output" >&2
    echo "[smoke]       Listing $OUT_DIR for diagnostic:" >&2
    ls -laR "$OUT_DIR" >&2 || true
    exit 67
fi

echo "[smoke] OK: $EXAMPLE_DIR produced $NON_EMPTY_FILES file(s)"
exit 0
