#!/bin/bash
# tests/diagnostics/capture_parallel_ksp_baseline.sh
#
# Pass-12 followup 3 (parallel KSP / PC fix): diagnostic capture
# helper for the EXAMPLE_SMOKE_MPI4_KNOWN_BROKEN "SNES diverges step 0"
# examples.
#
# Runs one examples/<N>_<event>/run.sh at MPI_RANKS=4 under a smoke
# budget, with a verbose PETSc monitor string, and writes the combined
# stdout / stderr to
#   tests/diagnostics/parallel_ksp_baseline/<N>_<event>.mpi4.log
#
# The script always exits 0 as long as the log file was produced and is
# non-empty. The example itself may converge (post-fix) or diverge
# (baseline); either way the failure signature or the successful
# convergence is captured verbatim.
#
# Usage:
#   tests/diagnostics/capture_parallel_ksp_baseline.sh <example_dir> [<final_time>]
#
# The monitor string is fixed so the captured logs are comparable and
# stay small enough to commit as fixtures:
#   -snes_monitor                show the SNES residual trajectory
#   -ksp_monitor_true_residual   show the KSP residual trajectory (the
#                                preconditioned + true residual pair
#                                makes the PREONLY-vs-GMRES difference
#                                obvious)
#   -snes_converged_reason       print the terminal SNES reason code
#   -ksp_converged_reason        print the terminal KSP reason code
#   -snes_linesearch_monitor     show the line-search step / gnorm
#   -snes_view                   dump the full SNES/KSP/PC/Mat tree once
#                                after the (first) SNESSolve: this is
#                                where the KSP type (preonly vs gmres),
#                                the PC composition (bjacobi + per-rank
#                                LU), and the Jacobian sparsity
#                                (rows/cols/nonzeros) are recorded
#   -snes_max_it 12              cap the Newton loop so a diverging run
#                                terminates within the test timeout
#   -ksp_max_it 200              cap the Krylov loop similarly
#   -ts_max_steps 1              run only the first TS step: enough to
#                                show step-0 convergence (post-fix) or
#                                the step-0 divergence (baseline)
#                                without dumping the whole transient
#   -options_left                flag any unused options
#
# FSRM_MPI_PETSC_OPTS is intentionally NOT set here: the run inherits
# whatever scripts/run_with_mpi.sh injects (the production default), so
# the baseline log shows the production parallel PC path and the
# post-fix log shows the fixed one.

set -uo pipefail

REPO_DIR="$(cd "$(dirname "$0")/../.." && pwd)"

if [ "$#" -lt 1 ]; then
    echo "[ksp-baseline] usage: $0 <example_dir> [<final_time>]" >&2
    exit 64
fi

EXAMPLE_DIR="$1"
# Make the path absolute so we can chdir freely.
case "$EXAMPLE_DIR" in
    /*) : ;;
    *)  EXAMPLE_DIR="$(cd "$EXAMPLE_DIR" && pwd)" ;;
esac
FINAL_TIME="${2:-0.05}"

if [ ! -d "$EXAMPLE_DIR" ]; then
    echo "[ksp-baseline] FAIL: $EXAMPLE_DIR is not a directory" >&2
    exit 65
fi
if [ ! -f "$EXAMPLE_DIR/run.sh" ]; then
    echo "[ksp-baseline] FAIL: $EXAMPLE_DIR/run.sh is missing" >&2
    exit 66
fi

EX_NAME="$(basename "$EXAMPLE_DIR")"
LOG_DIR="$REPO_DIR/tests/diagnostics/parallel_ksp_baseline"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/${EX_NAME}.mpi4.log"

rm -rf "$EXAMPLE_DIR/output"
mkdir -p "$EXAMPLE_DIR/output"

export MPI_RANKS=4
export FSRM_FINAL_TIME_OVERRIDE="$FINAL_TIME"
# Honour root-in-Docker OpenMPI just like scripts/run_with_mpi.sh does.
export OMPI_ALLOW_RUN_AS_ROOT="${OMPI_ALLOW_RUN_AS_ROOT:-1}"
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="${OMPI_ALLOW_RUN_AS_ROOT_CONFIRM:-1}"
export PETSC_OPTIONS="-snes_monitor -ksp_monitor_true_residual -snes_converged_reason -ksp_converged_reason -snes_linesearch_monitor -snes_view -snes_max_it 12 -ksp_max_it 200 -ts_max_steps 1 -options_left"

{
    echo "=================================================================="
    echo "FSRM parallel-KSP diagnostic capture"
    echo "  example      : $EX_NAME"
    echo "  MPI ranks    : $MPI_RANKS"
    echo "  final time   : $FINAL_TIME (FSRM_FINAL_TIME_OVERRIDE)"
    echo "  PETSC_OPTIONS: $PETSC_OPTIONS"
    echo "  date         : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "=================================================================="
    echo
    ( cd "$EXAMPLE_DIR" && bash run.sh )
    rc=$?
    echo
    echo "=== run.sh exit code: $rc ==="
} > "$LOG_FILE" 2>&1

if [ ! -s "$LOG_FILE" ]; then
    echo "[ksp-baseline] FAIL: no log produced at $LOG_FILE" >&2
    exit 67
fi

echo "[ksp-baseline] OK: captured $(wc -l < "$LOG_FILE") lines to $LOG_FILE"
exit 0
