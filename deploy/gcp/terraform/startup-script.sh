#!/bin/bash
# Startup script for GCE instance initialization

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
    htop \
    iotop \
    sysstat \
    gnuplot \
    docker.io

# Install Google Cloud SDK (if not already installed)
if ! command -v gcloud &> /dev/null; then
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | \
        tee -a /etc/apt/sources.list.d/google-cloud-sdk.list
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | \
        apt-key --keyring /usr/share/keyrings/cloud.google.gpg add -
    apt-get update && apt-get install -y google-cloud-sdk
fi

# Install Python packages
pip3 install numpy matplotlib h5py google-cloud-storage

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
if ! id -u ubuntu &>/dev/null; then
    useradd -m -s /bin/bash ubuntu
fi

cat >> /home/ubuntu/.bashrc << EOF

# PETSc environment
export PETSC_DIR=/opt/petsc-3.20.0
export PETSC_ARCH=arch-linux-c-opt
export PKG_CONFIG_PATH=/opt/petsc-3.20.0/arch-linux-c-opt/lib/pkgconfig:\$PKG_CONFIG_PATH
export PATH=\$PATH:/home/ubuntu/FSRM/build

# GCP configuration
export GCS_BUCKET=${bucket_name}
EOF

# Clone and build FSRM
cd /home/ubuntu
git clone https://github.com/yourusername/FSRM.git || \
    echo "Repository not yet public - will need to upload manually"

# Create directories
mkdir -p /home/ubuntu/simulations/{input,output}
mkdir -p /data

# Mount local SSDs if available
if [ -d /dev/disk/by-id ]; then
    for ssd in /dev/disk/by-id/google-local-ssd-*; do
        if [ -e "$ssd" ]; then
            mkfs.ext4 -F "$ssd"
            mkdir -p /mnt/localssd
            mount "$ssd" /mnt/localssd
            chmod 777 /mnt/localssd
        fi
    done
fi

# Set up GCS sync script
cat > /home/ubuntu/sync_to_gcs.sh << 'EOFSCRIPT'
#!/bin/bash
# Sync simulation outputs to Google Cloud Storage
gsutil -m rsync -r /home/ubuntu/simulations/output/ gs://${bucket_name}/outputs/
EOFSCRIPT

chmod +x /home/ubuntu/sync_to_gcs.sh

# Set up automatic GCS sync on shutdown
cat > /etc/systemd/system/gcs-sync-shutdown.service << EOFSERVICE
[Unit]
Description=Sync simulation data to GCS on shutdown
DefaultDependencies=no
Before=shutdown.target

[Service]
Type=oneshot
ExecStart=/home/ubuntu/sync_to_gcs.sh
User=ubuntu

[Install]
WantedBy=shutdown.target
EOFSERVICE

systemctl enable gcs-sync-shutdown.service

# Fix permissions
chown -R ubuntu:ubuntu /home/ubuntu

# Log completion
echo "Instance initialization completed at $(date)" > /var/log/startup-script.log
