# FSRM Production Docker Image (multi-stage)
# Stage 1 builds PETSc main (gitlab.com/petsc/petsc) at the pinned SHA below
# and compiles FSRM. Stage 2 ships only runtime artifacts plus visualization
# packages (numpy, matplotlib, h5py, scipy, pandas, obspy, pygmt, gmt).
#
# Build: docker build -t fsrm:local .
# Run:   docker run --rm -it fsrm:local
FROM ubuntu:22.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get clean && \
    apt-get update --fix-missing && \
    apt-get install -y --fix-missing \
    build-essential \
    cmake \
    git \
    wget \
    pkg-config \
    libopenmpi-dev \
    openmpi-bin \
    libhdf5-openmpi-dev \
    liblapack-dev \
    libblas-dev \
    python3 \
    python3-pip \
    gfortran \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# PETSc main pinned to commit 267b8824abdb82fedbefdec8ea8d0cb3c505ddfd
# (2026-04-17, "Merge remote-tracking branch 'origin/release'").
ARG PETSC_SHA=267b8824abdb82fedbefdec8ea8d0cb3c505ddfd
RUN git clone https://gitlab.com/petsc/petsc.git /opt/petsc-main \
    && cd /opt/petsc-main \
    && git checkout ${PETSC_SHA} \
    && git rev-parse HEAD > /opt/petsc-main/COMMIT_SHA

WORKDIR /opt/petsc-main
RUN ./configure \
        --with-cc=mpicc \
        --with-cxx=mpicxx \
        --with-fc=mpif90 \
        --with-mpi=1 \
        --download-fblaslapack \
        --download-ctetgen \
        --with-hdf5-include=/usr/include/hdf5/openmpi \
        --with-hdf5-lib="-L/usr/lib/x86_64-linux-gnu/hdf5/openmpi -lhdf5" \
        --with-debugging=0 \
        COPTFLAGS='-O3 -march=native' \
        CXXOPTFLAGS='-O3 -march=native' \
        FOPTFLAGS='-O3 -march=native' \
    && make PETSC_DIR=/opt/petsc-main PETSC_ARCH=arch-linux-c-opt all

ENV PETSC_DIR=/opt/petsc-main
ENV PETSC_ARCH=arch-linux-c-opt
ENV PKG_CONFIG_PATH=/opt/petsc-main/arch-linux-c-opt/lib/pkgconfig:${PKG_CONFIG_PATH}

WORKDIR /build
COPY . .

RUN mkdir -p build && cd build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON && \
    make -j$(nproc)

# Stage 2: Runtime environment
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get clean && \
    apt-get update --fix-missing && \
    apt-get install -y --fix-missing \
    libopenmpi3 \
    openmpi-bin \
    libhdf5-openmpi-103 \
    liblapack3 \
    libblas3 \
    libgomp1 \
    gnuplot \
    python3 \
    python3-pip \
    gmt \
    gmt-dcw \
    gmt-gshhg \
    ghostscript \
    && rm -rf /var/lib/apt/lists/*

RUN if [ -f /usr/lib/x86_64-linux-gnu/libgmt.so.6 ] && [ ! -f /usr/lib/x86_64-linux-gnu/libgmt.so ]; then \
        ln -sf /usr/lib/x86_64-linux-gnu/libgmt.so.6 /usr/lib/x86_64-linux-gnu/libgmt.so && ldconfig; \
    fi

RUN pip3 install --break-system-packages \
    numpy \
    matplotlib \
    h5py \
    scipy \
    pandas \
    obspy \
    pygmt

COPY --from=builder /opt/petsc-main/arch-linux-c-opt/lib /usr/local/lib/petsc
COPY --from=builder /opt/petsc-main/COMMIT_SHA /usr/local/lib/petsc/COMMIT_SHA
ENV LD_LIBRARY_PATH=/usr/local/lib/petsc:$LD_LIBRARY_PATH

COPY --from=builder /build/build /app
COPY --from=builder /build/config /app/config
COPY --from=builder /build/docs /app/docs
COPY --from=builder /build/examples /app/examples

RUN mkdir -p /app/output /data

WORKDIR /app

RUN useradd -m -s /bin/bash simuser && \
    chown -R simuser:simuser /app /data
USER simuser

CMD ["/bin/bash"]
