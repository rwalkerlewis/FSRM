#!/bin/bash
# User data script for EC2 instance initialization

set -e

# Update system
apt-get update
apt-get upgrade -y

# Install dependencies
DEBIAN_FRONTEND=noninteractive apt-get install -y \
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
    awscli \
    htop \
    iotop \
    sysstat \
    gnuplot \
    docker.io

# Install Python packages
pip3 install numpy matplotlib h5py boto3

# Install PETSc
cd /opt
wget https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.20.0.tar.gz
tar -xzf petsc-3.20.0.tar.gz
rm petsc-3.20.0.tar.gz

cd petsc-3.20.0
./configure \
    --with-cc=mpicc \
    --with-cxx=mpicxx \
    --with-fc=mpif90 \
    --download-fblaslapack \
    --with-debugging=0 \
    COPTFLAGS='-O3 -march=native' \
    CXXOPTFLAGS='-O3 -march=native' \
    FOPTFLAGS='-O3 -march=native'

make PETSC_DIR=/opt/petsc-3.20.0 PETSC_ARCH=arch-linux-c-opt all

# Set up environment variables
cat >> /etc/environment << EOF
PETSC_DIR=/opt/petsc-3.20.0
PETSC_ARCH=arch-linux-c-opt
PKG_CONFIG_PATH=/opt/petsc-3.20.0/arch-linux-c-opt/lib/pkgconfig:\$PKG_CONFIG_PATH
EOF

# Set up for ubuntu user
cat >> /home/ubuntu/.bashrc << EOF

# PETSc environment
export PETSC_DIR=/opt/petsc-3.20.0
export PETSC_ARCH=arch-linux-c-opt
export PKG_CONFIG_PATH=/opt/petsc-3.20.0/arch-linux-c-opt/lib/pkgconfig:\$PKG_CONFIG_PATH
export PATH=\$PATH:/home/ubuntu/ReservoirSim/build

# AWS configuration
export S3_BUCKET=${s3_bucket}
EOF

# Clone and build ReservoirSim
cd /home/ubuntu
git clone https://github.com/yourusername/ReservoirSim.git || \
    echo "Repository not yet public - will need to upload manually"

# Create directories
mkdir -p /home/ubuntu/simulations/{input,output}
mkdir -p /data

# Set up S3 sync script
cat > /home/ubuntu/sync_to_s3.sh << 'EOFSCRIPT'
#!/bin/bash
# Sync simulation outputs to S3
aws s3 sync /home/ubuntu/simulations/output/ s3://${s3_bucket}/outputs/ --exclude "*" --include "*.vtu" --include "*.h5" --include "*.log"
EOFSCRIPT

chmod +x /home/ubuntu/sync_to_s3.sh

# Set up automatic S3 sync on shutdown
cat > /etc/systemd/system/s3-sync-shutdown.service << EOFSERVICE
[Unit]
Description=Sync simulation data to S3 on shutdown
DefaultDependencies=no
Before=shutdown.target

[Service]
Type=oneshot
ExecStart=/home/ubuntu/sync_to_s3.sh
User=ubuntu

[Install]
WantedBy=shutdown.target
EOFSERVICE

systemctl enable s3-sync-shutdown.service

# Fix permissions
chown -R ubuntu:ubuntu /home/ubuntu

# Log completion
echo "Instance initialization completed at $(date)" > /var/log/user-data.log
