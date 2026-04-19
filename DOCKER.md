# FSRM Docker Images

All four Dockerfiles build PETSc from gitlab.com/petsc/petsc `main` pinned to
commit `267b8824abdb82fedbefdec8ea8d0cb3c505ddfd` (dated 2026-04-17). They
share the same configure flags except where noted; only the surrounding base
image, optimization level, and ancillary tooling differ.

| File | Purpose | Base image | When to use |
|------|---------|------------|-------------|
| `Dockerfile` | Production multi-stage build (compiles FSRM, ships runtime + visualization Python) | `ubuntu:22.04` | Distributing a runnable image of FSRM with `numpy`, `matplotlib`, `obspy`, `pygmt`, etc. |
| `Dockerfile.ci` | Build environment only, FSRM source mounted at runtime | `ubuntu:22.04` | CI; matches `CLAUDE.md` build/test instructions |
| `Dockerfile.cuda` | CI build env + PETSc CUDA backend (`--with-cuda`) | `nvidia/cuda:12.2.0-devel-ubuntu22.04` | GPU runs (`-vec_type cuda -mat_type aijcusparse`) |
| `Dockerfile.dev` | Full developer environment (gmsh, gdb, valgrind, editors, oh-my-zsh) | `nvidia/cuda:12.6.3-devel-ubuntu24.04` | VS Code Dev Containers; interactive development |

## Bumping PETSc

Update the `PETSC_SHA` ARG in all four Dockerfiles and in this document. The
SHA is pinned (not a branch) so the build is reproducible.

## Why not the apt PETSc package

`Dockerfile.dev` previously used `libpetsc-real-dev` from Ubuntu apt. That
shipped a different PETSc version than the CI/production images, which masked
API drift between containers. All four images now build the same source.
