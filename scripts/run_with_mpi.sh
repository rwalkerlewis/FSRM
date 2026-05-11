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
# Also injects parallel-friendly PETSc KSP/PC options when MPI_RANKS>1,
# because the supported PETSc 3.25 build does not include MUMPS or
# SuperLU_DIST and so the default -pc_type lu (serial KLU) fails on
# MPIAIJ matrices. The override is
#   -ksp_type gmres -pc_type bjacobi -sub_pc_type lu
# i.e. block-Jacobi (per-rank LU) wrapped in a GMRES Krylov loop.
#
# Pass-12 followup 3: the previous default omitted -ksp_type. Several
# physics paths (Simulator::setupSolvers, elastodynamics branch) set
# KSPPREONLY + PCLU, so a serial run does an exact direct solve and
# Newton converges quadratically. In parallel, -pc_type bjacobi
# overrides PCLU but the KSP type stays PREONLY, which means each
# "Newton step" is a single block-Jacobi application rather than a
# converged linear solve. That is a fixed-point sweep, not a Newton
# step, and it stalls (SNES DIVERGED_MAX_IT at step 0). Forcing
# -ksp_type gmres restores a real Krylov solve around the block-Jacobi
# preconditioner. See docs/PARALLEL_KSP.md for the diagnosis, the PC
# sweep, and the per-example evidence.
#
# Override or extend by setting FSRM_MPI_PETSC_OPTS in the environment,
# or pass `-ksp_type ... -pc_type ...` directly on the command line
# (the helper appends its defaults; PETSc CLI lets later flags win).
#
# Honors OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 from the environment if set.
#
# Usage:
#   source scripts/run_with_mpi.sh
#   run_with_mpi /path/to/fsrm -c config.config [args...]
#
# Configuration via env vars:
#   MPI_RANKS               - number of ranks (default: 4)
#   FSRM_MPI_PETSC_OPTS     - PETSc options for parallel runs.
#                             Default:
#                               "-ksp_type gmres -pc_type bjacobi -sub_pc_type lu"
#                             Set to empty to disable injection.
#   MPI_DEBUG               - set non-empty to print the chosen mpirun command
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

    # Parallel-friendly PETSc options. The PETSc 3.25 build linked
    # into the FSRM image does not include MUMPS or SuperLU_DIST, so
    # the default -pc_type lu fails on MPIAIJ matrices. Inject a GMRES
    # Krylov loop around block-Jacobi + per-rank LU as the default.
    # -ksp_type gmres is load-bearing: it overrides the KSPPREONLY that
    # Simulator::setupSolvers sets on the elastodynamics path, which
    # would otherwise turn each Newton step into a single block-Jacobi
    # sweep (see the header comment and docs/PARALLEL_KSP.md). Disable
    # the injection entirely by setting FSRM_MPI_PETSC_OPTS to "" in
    # the environment.
    local default_petsc_opts="-ksp_type gmres -pc_type bjacobi -sub_pc_type lu"
    local petsc_opts="${FSRM_MPI_PETSC_OPTS-$default_petsc_opts}"
    local petsc_args=()
    if [ -n "$petsc_opts" ]; then
        # shellcheck disable=SC2206
        petsc_args=($petsc_opts)
    fi

    [ -n "$MPI_DEBUG" ] && \
        echo "[run_with_mpi] $impl, ranks=$ranks: mpirun -n $ranks ${mpi_args[*]} $* ${petsc_args[*]}" >&2

    mpirun -n "$ranks" "${mpi_args[@]}" "$@" "${petsc_args[@]}"
}

# Allow direct invocation: bash scripts/run_with_mpi.sh /path/to/fsrm ...
if [ "${BASH_SOURCE[0]}" = "$0" ]; then
    run_with_mpi "$@"
fi
