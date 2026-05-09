#!/bin/bash
# scripts/run_with_mpi.sh
#
# Centralized MPI launcher for FSRM example run.sh scripts.
#
# Detects whether the linked MPI is OpenMPI or MPICH and applies the
# right flags: --allow-run-as-root when running as uid 0 in Docker,
# --bind-to core for predictable per-rank performance, --oversubscribe
# if MPI_RANKS exceeds the detected core count. Falls back to direct
# invocation when MPI_RANKS=1.
#
# Honors OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 from the environment if set.
#
# Usage:
#   source scripts/run_with_mpi.sh
#   run_with_mpi /path/to/fsrm -c config.config [args...]
#
# Configuration via env vars:
#   MPI_RANKS  - number of ranks (default: 4)
#   MPI_DEBUG  - set non-empty to print the chosen mpirun command
#
# Source script (intended to be sourced by example run.sh files).

run_with_mpi() {
    local ranks="${MPI_RANKS:-4}"

    if [ "$ranks" -le 1 ]; then
        [ -n "$MPI_DEBUG" ] && echo "[run_with_mpi] direct: $*" >&2
        "$@"
        return $?
    fi

    if ! command -v mpirun >/dev/null 2>&1; then
        echo "[run_with_mpi] WARNING: mpirun not found; running serially" >&2
        "$@"
        return $?
    fi

    local mpi_args=()
    local impl="unknown"

    local mpi_version
    mpi_version="$(mpirun --version 2>&1 | head -3)"

    if echo "$mpi_version" | grep -qiE 'open[- ]?mpi'; then
        impl="openmpi"
        if [ "$(id -u)" -eq 0 ]; then
            mpi_args+=("--allow-run-as-root")
        fi
        mpi_args+=("--bind-to" "core")
        local cores
        cores="$(nproc 2>/dev/null || echo 1)"
        if [ "$ranks" -gt "$cores" ]; then
            mpi_args+=("--oversubscribe")
        fi
    elif echo "$mpi_version" | grep -qiE 'mpich|hydra'; then
        impl="mpich"
    fi

    [ -n "$MPI_DEBUG" ] && \
        echo "[run_with_mpi] $impl, ranks=$ranks: mpirun -n $ranks ${mpi_args[*]} $*" >&2

    mpirun -n "$ranks" "${mpi_args[@]}" "$@"
}

# Allow direct invocation: bash scripts/run_with_mpi.sh /path/to/fsrm ...
if [ "${BASH_SOURCE[0]}" = "$0" ]; then
    run_with_mpi "$@"
fi
