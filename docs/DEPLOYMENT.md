# FSRM Deployment Guide

Comprehensive guide for deploying FSRM on cloud and HPC platforms.

---

## Table of Contents

1. [Local Installation](#local-installation)
2. [Docker Deployment](#docker-deployment)
3. [AWS Deployment](#aws-deployment)
4. [Google Cloud Deployment](#google-cloud-deployment)
5. [HPC Clusters](#hpc-clusters)
6. [GPU Configuration](#gpu-configuration)
7. [Scaling Guidelines](#scaling-guidelines)

---

## Local Installation

### Prerequisites

- C++17 compiler (GCC 9+, Clang 10+)
- CMake 3.16+
- MPI implementation (OpenMPI or MPICH)
- PETSc 3.18+ with MPI

### Build Steps

```bash
# Clone repository
git clone https://github.com/your-org/fsrm.git
cd fsrm

# Set PETSc environment
export PETSC_DIR=/path/to/petsc
export PETSC_ARCH=arch-linux-c-opt

# Build
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)

# Test
ctest --output-on-failure

# Install (optional)
sudo make install
```

### Quick Test

```bash
mpirun -np 4 ./fsrm -c ../config/default.config
```

---

## Docker Deployment

### Using Pre-built Image

```bash
# Pull image
docker pull ghcr.io/your-org/fsrm:latest

# Run simulation
docker run -v $(pwd)/data:/data fsrm \
  mpirun -np 4 fsrm -c /data/config.config

# Interactive shell
docker run -it -v $(pwd)/data:/data fsrm bash
```

### Building Custom Image

```dockerfile
# Dockerfile
FROM ubuntu:22.04

RUN apt-get update && apt-get install -y \
    build-essential cmake git \
    libopenmpi-dev openmpi-bin \
    libpetsc-real-dev \
    && rm -rf /var/lib/apt/lists/*

COPY . /fsrm
WORKDIR /fsrm/build

RUN cmake .. -DCMAKE_BUILD_TYPE=Release && make -j$(nproc)

ENTRYPOINT ["./fsrm"]
```

```bash
docker build -t fsrm:custom .
docker run -v $(pwd)/data:/data fsrm:custom -c /data/config.config
```

### Docker Compose (Multi-node)

```yaml
# docker-compose.yml
version: '3.8'
services:
  head:
    image: fsrm:latest
    volumes:
      - ./data:/data
      - ./results:/results
    command: >
      mpirun --hostfile /data/hostfile -np 16
      fsrm -c /data/config.config
    networks:
      - mpi-network

  worker:
    image: fsrm:latest
    deploy:
      replicas: 3
    networks:
      - mpi-network

networks:
  mpi-network:
    driver: overlay
```

---

## AWS Deployment

### Option 1: EC2 with Terraform

**Infrastructure setup:**

```hcl
# terraform/main.tf
provider "aws" {
  region = var.region
}

resource "aws_instance" "fsrm_head" {
  ami           = data.aws_ami.ubuntu.id
  instance_type = var.instance_type
  key_name      = var.key_name

  vpc_security_group_ids = [aws_security_group.fsrm.id]
  
  user_data = file("user_data.sh")

  root_block_device {
    volume_size = 100
    volume_type = "gp3"
  }

  tags = {
    Name = "fsrm-head"
  }
}

resource "aws_instance" "fsrm_worker" {
  count         = var.worker_count
  ami           = data.aws_ami.ubuntu.id
  instance_type = var.instance_type
  key_name      = var.key_name

  vpc_security_group_ids = [aws_security_group.fsrm.id]
  
  user_data = file("user_data.sh")

  tags = {
    Name = "fsrm-worker-${count.index}"
  }
}
```

**Deployment:**

```bash
cd deploy/terraform
terraform init
terraform plan -var="worker_count=4"
terraform apply
```

### Option 2: AWS ParallelCluster

**Cluster configuration:**

```yaml
# parallelcluster-config.yaml
Region: us-east-1
Image:
  Os: ubuntu2204
HeadNode:
  InstanceType: c6i.xlarge
  Networking:
    SubnetId: subnet-xxxxx
  Ssh:
    KeyName: your-key
Scheduling:
  Scheduler: slurm
  SlurmQueues:
    - Name: compute
      ComputeResources:
        - Name: c6i-large
          InstanceType: c6i.4xlarge
          MinCount: 0
          MaxCount: 100
      Networking:
        SubnetIds:
          - subnet-xxxxx
```

**Commands:**

```bash
# Create cluster
pcluster create-cluster -n fsrm-cluster -c parallelcluster-config.yaml

# Connect
pcluster ssh -n fsrm-cluster

# Submit job
sbatch --nodes=4 --ntasks-per-node=16 run_fsrm.sh
```

### Option 3: AWS Batch

```bash
# Submit batch job
aws batch submit-job \
  --job-name fsrm-simulation \
  --job-queue fsrm-queue \
  --job-definition fsrm-definition \
  --container-overrides '{
    "command": ["mpirun", "-np", "64", "fsrm", "-c", "/data/config.config"],
    "resourceRequirements": [
      {"type": "VCPU", "value": "64"},
      {"type": "MEMORY", "value": "128000"}
    ]
  }'
```

### Recommended EC2 Instances

| Workload | Instance | vCPU | Memory | Notes |
|----------|----------|------|--------|-------|
| Small | c6i.2xlarge | 8 | 16 GB | Dev/testing |
| Medium | c6i.8xlarge | 32 | 64 GB | Production |
| Large | c6i.24xlarge | 96 | 192 GB | Large models |
| GPU | p4d.24xlarge | 96 | 1152 GB | 8x A100 GPUs |
| HPC | hpc6a.48xlarge | 96 | 384 GB | EFA networking |

---

## Google Cloud Deployment

### Option 1: Compute Engine with Terraform

```hcl
# terraform/main.tf
provider "google" {
  project = var.project_id
  region  = var.region
}

resource "google_compute_instance" "fsrm" {
  count        = var.node_count
  name         = "fsrm-node-${count.index}"
  machine_type = var.machine_type
  zone         = var.zone

  boot_disk {
    initialize_params {
      image = "ubuntu-2204-lts"
      size  = 100
    }
  }

  network_interface {
    network = "default"
    access_config {}
  }

  metadata_startup_script = file("startup-script.sh")
}
```

### Option 2: GKE with MPI Operator

```yaml
# mpi-job.yaml
apiVersion: kubeflow.org/v2beta1
kind: MPIJob
metadata:
  name: fsrm-simulation
spec:
  slotsPerWorker: 4
  runPolicy:
    cleanPodPolicy: Running
  mpiReplicaSpecs:
    Launcher:
      replicas: 1
      template:
        spec:
          containers:
          - image: gcr.io/your-project/fsrm:latest
            name: launcher
            command:
            - mpirun
            - -np
            - "16"
            - fsrm
            - -c
            - /data/config.config
    Worker:
      replicas: 4
      template:
        spec:
          containers:
          - image: gcr.io/your-project/fsrm:latest
            name: worker
            resources:
              limits:
                cpu: "4"
                memory: 16Gi
```

### Recommended GCP Instances

| Workload | Machine Type | vCPU | Memory |
|----------|--------------|------|--------|
| Small | c2-standard-8 | 8 | 32 GB |
| Medium | c2-standard-30 | 30 | 120 GB |
| Large | c2-standard-60 | 60 | 240 GB |
| GPU | a2-highgpu-8g | 96 | 680 GB |

---

## HPC Clusters

### SLURM Job Script

```bash
#!/bin/bash
#SBATCH --job-name=fsrm
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=32
#SBATCH --time=24:00:00
#SBATCH --partition=compute

module load gcc/11.2.0
module load openmpi/4.1.1
module load petsc/3.18.0

cd $SLURM_SUBMIT_DIR

mpirun -np $SLURM_NTASKS ./fsrm -c simulation.config \
  -log_view :performance.log
```

### PBS/Torque Job Script

```bash
#!/bin/bash
#PBS -N fsrm
#PBS -l nodes=4:ppn=32
#PBS -l walltime=24:00:00
#PBS -q batch

module load gcc openmpi petsc

cd $PBS_O_WORKDIR

mpirun -np 128 ./fsrm -c simulation.config
```

### LSF Job Script

```bash
#!/bin/bash
#BSUB -J fsrm
#BSUB -n 128
#BSUB -R "span[ptile=32]"
#BSUB -W 24:00

module load gcc openmpi petsc

mpirun -np 128 ./fsrm -c simulation.config
```

---

## GPU Configuration

> **Architecture Details:** See [PHYSICS_AND_GPU_ARCHITECTURE.md](PHYSICS_AND_GPU_ARCHITECTURE.md) for comprehensive GPU physics kernel documentation.

### Build with CUDA (NVIDIA GPUs)

```bash
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DENABLE_CUDA=ON \
  -DCMAKE_CUDA_ARCHITECTURES="70;80;86;90"

make -j$(nproc)
```

**Supported CUDA Architectures:**

| Architecture | GPUs | Compute Capability |
|--------------|------|-------------------|
| `70` | V100 | 7.0 |
| `80` | A100 | 8.0 |
| `86` | RTX 30xx | 8.6 |
| `89` | RTX 40xx | 8.9 |
| `90` | H100 | 9.0 |

### Build with ROCm/HIP (AMD GPUs)

```bash
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DENABLE_HIP=ON \
  -DCMAKE_HIP_ARCHITECTURES="gfx908;gfx90a;gfx942"

make -j$(nproc)
```

**Supported AMD Architectures:**

| Architecture | GPUs |
|--------------|------|
| `gfx908` | MI100 |
| `gfx90a` | MI210/MI250 |
| `gfx942` | MI300 |

### GPU Execution Modes

FSRM supports four GPU execution modes:

| Mode | Description | Use Case |
|------|-------------|----------|
| `CPU_ONLY` | Run entirely on CPU | Testing, small problems |
| `GPU_ONLY` | Run entirely on GPU | Large problems, dedicated GPU |
| `CPU_FALLBACK` | Try GPU, fallback to CPU | Production (default) |
| `HYBRID` | Load balance CPU+GPU | Multi-GPU + CPU systems |

### Configuration for GPU

```ini
[SIMULATION]
# Enable GPU acceleration
use_gpu = true
gpu_mode = CPU_FALLBACK       # AUTO, GPU_ONLY, CPU_FALLBACK, HYBRID
gpu_device_id = 0             # GPU device to use (0 = first GPU)
gpu_memory_fraction = 0.8     # Fraction of GPU memory to use
gpu_verbose = false           # Print GPU info at startup

# GPU solver options
use_gpu_preconditioner = true
use_gpu_matrix_assembly = true
pin_host_memory = true        # Faster CPU-GPU transfers

# Physics-specific GPU options
enable_elastodynamics = true
enable_poroelastodynamics = false
```

### GPU-Accelerated Physics Kernels

| Kernel | Status | Speedup | Min Elements |
|--------|--------|---------|--------------|
| Single-Phase Flow | ✓ Available | 2-5× | 1,000 |
| Elastodynamics | ✓ Available | 5-20× | 1,000 |
| Poroelastodynamics | ✓ Available | 5-50× | 1,000 |
| Black Oil | Planned | 3-8× | — |
| Compositional | Planned | 5-15× | — |
| Thermal | Planned | 2-5× | — |

### Multi-GPU Setup

```bash
# One MPI rank per GPU (recommended)
mpirun -np 4 ./fsrm -c config.config

# Query GPU mapping
nvidia-smi -L  # List available GPUs
```

**MPI-GPU Binding:**
```bash
# Explicit GPU assignment per rank
mpirun -np 4 --bind-to numa \
  --map-by ppr:1:numa \
  ./fsrm -c config.config
```

### GPU Memory Estimation

```
GPU Memory (bytes) ≈ n_cells × bytes_per_cell × factor

Bytes per cell by physics:
- Single-phase flow:    ~100 bytes
- Elastodynamics:       ~300 bytes (displacement, velocity, strain, stress)
- Poroelastodynamics:   ~500 bytes (solid + fluid arrays)
- Black oil:            ~400 bytes (three phases + properties)

Factor: 1.5-2.0 for working memory

Example: 1M cells with poroelastodynamics
  1,000,000 × 500 × 2.0 = 1 GB GPU memory required
```

### GPU Performance Tips

1. **Large grids benefit most**: Use grids > 100K cells for significant speedup
2. **Keep data on GPU**: Minimize CPU-GPU transfers between timesteps
3. **Batch operations**: GPU kernels process all elements simultaneously
4. **Use appropriate precision**: Double precision for physics, single for visualization
5. **Profile execution**: Use `--gpu-verbose` flag to monitor GPU utilization

### Verifying GPU Usage

```bash
# Check GPU support at runtime
./fsrm --version          # Shows GPU support status
./fsrm --gpu-info         # Lists available GPUs

# Monitor GPU during simulation
nvidia-smi -l 1           # Update every 1 second (NVIDIA)
rocm-smi -l               # AMD equivalent

# Enable verbose GPU logging
mpirun -np 4 ./fsrm -c config.config -gpu_verbose
```

### Troubleshooting GPU Issues

**"GPU not available" warning:**
```bash
# Check CUDA installation
nvidia-smi
nvcc --version

# Verify library paths
echo $LD_LIBRARY_PATH | grep -i cuda
```

**Out of GPU memory:**
```ini
# Reduce memory usage
gpu_memory_fraction = 0.5   # Use only 50% of GPU memory

# Or use CPU fallback for large problems
gpu_mode = CPU_FALLBACK
```

**Poor GPU performance:**
- Check problem size (too small may be slower on GPU)
- Verify MPI-GPU binding with `nvidia-smi`
- Profile with `nsys` (NVIDIA) or `rocprof` (AMD)

---

## Scaling Guidelines

### Problem Sizing

| Grid Size | Cells | Recommended Cores | Memory/Core |
|-----------|-------|-------------------|-------------|
| Small | < 100K | 4-16 | 2 GB |
| Medium | 100K-1M | 16-64 | 4 GB |
| Large | 1M-10M | 64-256 | 4-8 GB |
| Very Large | > 10M | 256+ | 8+ GB |

### Optimal Cells per Core

- **Flow only**: 5,000-20,000 cells/core
- **Flow + Geomechanics**: 2,000-10,000 cells/core
- **Full coupling**: 1,000-5,000 cells/core

### Scaling Example

```bash
# 1M cell model
# 50,000 cells/core → 20 cores minimum
# 10,000 cells/core → 100 cores optimal
# 5,000 cells/core → 200 cores good scaling

mpirun -np 100 ./fsrm -c million_cell.config
```

### Communication Overhead

- **Strong scaling limit**: ~1,000 cells/core
- **Weak scaling**: Maintain cells/core constant
- **Network**: Use InfiniBand for > 64 cores

### Memory Estimation

```
Memory ≈ cells × DOF × 8 bytes × factor

where:
  DOF = degrees of freedom (1 for flow, 4 for coupled)
  factor = 10-50 (matrix storage, solvers)

Example: 1M cells, coupled (4 DOF)
  1e6 × 4 × 8 × 30 ≈ 1 GB per core
```

---

## Monitoring and Debugging

### PETSc Logging

```bash
mpirun -np 16 ./fsrm -c config.config \
  -log_view :log.txt \
  -log_view_memory \
  -malloc_view
```

### Performance Profiling

```bash
# With Score-P
scorep mpirun -np 16 ./fsrm -c config.config

# With TAU
tau_exec mpirun -np 16 ./fsrm -c config.config
```

### Common Issues

**Slow convergence:**
```ini
[SIMULATION]
rtol = 1e-4              # Relax tolerance
max_nonlinear_iterations = 100
preconditioner = hypre   # Use AMG
```

**Out of memory:**
```bash
# Increase per-process memory
mpirun -np 8 ./fsrm -c config.config  # Fewer processes
```

**Load imbalance:**
- Use `-log_view` to check time per rank
- Consider graph partitioning for unstructured grids

---

## Cost Optimization

### Spot/Preemptible Instances

**AWS Spot:**
```bash
aws ec2 request-spot-instances \
  --instance-count 4 \
  --type "one-time" \
  --launch-specification file://spot-spec.json
```

**GCP Preemptible:**
```hcl
resource "google_compute_instance" "fsrm" {
  scheduling {
    preemptible = true
    automatic_restart = false
  }
}
```

### Checkpointing for Fault Tolerance

```ini
[OUTPUT]
checkpoint_frequency = 100   # Every 100 steps
checkpoint_path = checkpoints/
```

```bash
# Resume from checkpoint
mpirun -np 16 ./fsrm -c config.config -restart checkpoints/step_1000.h5
```

### Cost Estimates (AWS, 2024)

| Configuration | Instance | $/hour | Typical Job Cost |
|---------------|----------|--------|------------------|
| Small (4 cores) | c6i.xlarge | $0.17 | $1-5 |
| Medium (32 cores) | c6i.8xlarge | $1.36 | $10-50 |
| Large (96 cores) | c6i.24xlarge | $4.08 | $50-200 |
| GPU (8x A100) | p4d.24xlarge | $32.77 | $100-500 |

*Use spot instances for 60-80% savings*
