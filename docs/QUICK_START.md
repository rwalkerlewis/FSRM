# FSRM Quick Start

Get FSRM building, running, and producing seismograms in five minutes.

## Prerequisites

The supported build path is Docker. The CI image bundles a verified PETSc
3.25.0 build with MPI, HDF5, ctetgen, and Gmsh.

- Docker (or any OCI runtime)
- ~10 GB free disk for the build image and intermediate output

A native build is documented in [DEVELOPMENT.md](DEVELOPMENT.md) but is
not the recommended starting path.

## 1. Clone and build

```bash
git clone https://github.com/rwalkerlewis/FSRM.git
cd FSRM

docker build -f Dockerfile.ci -t fsrm-ci:local .
docker run --rm -v "$(pwd)":/workspace -w /workspace fsrm-ci:local bash -c \
  'mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release \
   -DENABLE_TESTING=ON -DENABLE_CUDA=OFF && make -j$(nproc)'
```

The `fsrm` executable lands in `build/fsrm`.

## 2. Run the test suite

```bash
docker run --rm -v "$(pwd)":/workspace -w /workspace/build fsrm-ci:local \
  ctest -j$(nproc) --output-on-failure
```

Expected: 110 of 116 tests pass. The six failures are pre-existing and
documented in [SOLVER_STATE.md](SOLVER_STATE.md).

## 3. Run an example end-to-end

Each `examples/N_<event>/` directory is self-contained: config, runner,
and output landing site live together. To run Salmon 1964 (the canonical
V&V anchor on the Marshak-tier source physics ladder):

```bash
docker run --rm -v "$(pwd)":/workspace \
  -w /workspace/examples/20_salmon_1964 fsrm-ci:local \
  ./run.sh
```

This produces:

- `output/seismograms/*.sac`: SAC-format synthetic seismograms.
- `output/near_field_history.csv`: 6-component moment-rate tensor and
  cavity-radius history from the radial Lagrangian solver.
- `output/near_field_profile.h5` + `.xdmf`: per-snapshot radial state of
  the source-ball solver (open in ParaView via the XDMF wrapper).
- `output/solution.h5` + `.xmf`: full 3-D wavefield (large; opt out via
  `[OUTPUT] hdf5_enabled = false`).

Run any historic-nuclear example by changing the directory. The
`examples/` index in [USER_GUIDE.md](USER_GUIDE.md) lists what each one
exercises.

## 4. Run in parallel

The `run.sh` scripts launch under MPI by default. Override the rank
count with the `MPI_RANKS` env var:

```bash
MPI_RANKS=8 ./run.sh
```

The default is `MPI_RANKS=4`; showcase scripts default to
`MPI_RANKS=8`. The 1-D radial source-ball solver remains serial; only
the FEM far-field is parallelized. See the "Running in parallel"
section of [USER_GUIDE.md](USER_GUIDE.md).

## 5. Visualize the output

The visualization scripts ship in `scripts/` and read the simulator
output without modification:

```bash
pip install matplotlib obspy h5py numpy
python3 scripts/plot_seismograms.py examples/20_salmon_1964/output/seismograms/
python3 scripts/plot_wavefield.py examples/20_salmon_1964/output/solution.h5
```

For the showcase events (Sedan 1962, Salmon 1964, Punggye-ri 2017,
Cannikin 1971, Sterling 1966) a `figures/regenerate.sh` script
regenerates a six-figure presentation pack from `output/` using the
shared style infrastructure in `tools/figures/`.

## Next steps

- [USER_GUIDE.md](USER_GUIDE.md): end-to-end manual covering config blocks,
  fidelity-ladder selection, and output catalogs.
- [FIDELITY_LADDER_GUIDE.md](FIDELITY_LADDER_GUIDE.md): one-page guide to
  picking LOW / MED / HIGH / HIGHEST tiers.
- [CONFIGURATION.md](CONFIGURATION.md): full config-key reference.
- [HISTORIC_NUCLEAR_FIDELITY.md](HISTORIC_NUCLEAR_FIDELITY.md): per-pass
  detail of what the source physics shipped pass-by-pass.
- [AXIS_1A_FIDELITY_REPORT.md](AXIS_1A_FIDELITY_REPORT.md): the
  canonical cross-pass V&V result.

## Troubleshooting

### Build fails on PETSc detection (native path)

The supported build is Docker. If you have a native PETSc, set
`PETSC_DIR` and `PETSC_ARCH` to a 3.25.0 build (older PETSc will fail
because the cohesive-cell BdResidual API and the `DMSetAuxiliaryVec`
signature change between minor versions).

### `mpirun` permission errors when running as root in Docker

The wrapper `scripts/run_with_mpi.sh` (sourced by every example
`run.sh`) detects the OpenMPI vs MPICH ABI and applies
`--allow-run-as-root` automatically. If you bypass the wrapper, add
`--allow-run-as-root` and `--bind-to core` manually under OpenMPI.

### `ctest` shows fault-test failures

The six listed in [SOLVER_STATE.md](SOLVER_STATE.md) are honest
pre-existing failures behind a PETSc 3.25 BdResidual limitation on
cohesive geometry. They do not block any historic-nuclear or
source-physics run.
