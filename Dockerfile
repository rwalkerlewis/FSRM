# Multi-stage Docker build for FSRM
# Stage 1: Build environment
FROM ubuntu:22.04 AS builder

# Avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install build dependencies
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
    && rm -rf /var/lib/apt/lists/*

# Install PETSc
WORKDIR /opt
RUN wget https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.20.0.tar.gz && \
    tar -xzf petsc-3.20.0.tar.gz && \
    rm petsc-3.20.0.tar.gz

WORKDIR /opt/petsc-3.20.0
RUN ./configure \
    --with-cc=mpicc \
    --with-cxx=mpicxx \
    --with-fc=mpif90 \
    --download-fblaslapack \
    --with-debugging=0 \
    COPTFLAGS='-O3 -march=native' \
    CXXOPTFLAGS='-O3 -march=native' \
    FOPTFLAGS='-O3 -march=native' && \
    make PETSC_DIR=/opt/petsc-3.20.0 PETSC_ARCH=arch-linux-c-opt all

ENV PETSC_DIR=/opt/petsc-3.20.0
ENV PETSC_ARCH=arch-linux-c-opt
ENV PKG_CONFIG_PATH=/opt/petsc-3.20.0/arch-linux-c-opt/lib/pkgconfig:$PKG_CONFIG_PATH

# Copy source code
WORKDIR /build
COPY . .

# Build FSRM
RUN mkdir -p build && cd build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON && \
    make -j$(nproc)

# Stage 2: Runtime environment
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Install runtime dependencies
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
    && rm -rf /var/lib/apt/lists/*

# Install Python packages for visualization
RUN pip3 install numpy matplotlib h5py

# Copy PETSc libraries
COPY --from=builder /opt/petsc-3.20.0/arch-linux-c-opt/lib /usr/local/lib/petsc
ENV LD_LIBRARY_PATH=/usr/local/lib/petsc:$LD_LIBRARY_PATH

# Copy built application
COPY --from=builder /build/build /app
COPY --from=builder /build/config /app/config
COPY --from=builder /build/docs /app/docs
COPY --from=builder /build/examples /app/examples

# Create output directory
RUN mkdir -p /app/output /data

WORKDIR /app

# Set up user (non-root for security)
RUN useradd -m -s /bin/bash simuser && \
    chown -R simuser:simuser /app /data
USER simuser

# Default command
CMD ["/bin/bash"]
