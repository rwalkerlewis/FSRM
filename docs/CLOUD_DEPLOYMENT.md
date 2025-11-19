# ReservoirSim Cloud Deployment Guide

This guide provides detailed instructions for deploying ReservoirSim on AWS and Google Cloud Platform for high-performance computing workloads.

## Table of Contents

1. [Overview](#overview)
2. [Prerequisites](#prerequisites)
3. [AWS Deployment](#aws-deployment)
4. [Google Cloud Deployment](#google-cloud-deployment)
5. [Docker Deployment](#docker-deployment)
6. [Cost Optimization](#cost-optimization)
7. [Monitoring and Management](#monitoring-and-management)
8. [Troubleshooting](#troubleshooting)

---

## Overview

ReservoirSim is a computationally intensive application that benefits from cloud deployment for:

- **Scalability**: Run large simulations on high-core-count instances
- **Flexibility**: Scale resources up/down based on simulation needs
- **Parallel Computing**: Leverage MPI across multiple nodes for massive simulations
- **Cost Efficiency**: Pay only for compute resources when needed
- **Data Storage**: Persistent storage in S3/GCS for simulation results

### Deployment Options

| Option | Best For | Cost | Complexity |
|--------|----------|------|------------|
| Single EC2/GCE Instance | Testing, small-to-medium simulations | Low | Low |
| AWS ParallelCluster | Large-scale parallel simulations | Medium-High | Medium |
| Docker on Cloud | Containerized, reproducible deployments | Low-Medium | Low-Medium |
| Kubernetes (GKE/EKS) | Production, multi-tenant, auto-scaling | Medium-High | High |

---

## Prerequisites

### Local Requirements

Before deploying, ensure you have the following installed on your local machine:

1. **Cloud CLI Tools**:
   - AWS: [AWS CLI](https://aws.amazon.com/cli/)
   - GCP: [Google Cloud SDK](https://cloud.google.com/sdk/docs/install)

2. **Infrastructure Tools**:
   - [Terraform](https://www.terraform.io/downloads) (>= 1.0)
   - [Docker](https://docs.docker.com/get-docker/) (for containerized deployments)

3. **Authentication**:
   - AWS: Configure credentials with `aws configure`
   - GCP: Authenticate with `gcloud auth login`

### Cloud Account Setup

#### AWS
```bash
# Configure AWS credentials
aws configure
# Enter: Access Key ID, Secret Access Key, Region, Output format

# Create an SSH key pair in the AWS console or via CLI
aws ec2 create-key-pair --key-name reservoirsim-key --query 'KeyMaterial' --output text > reservoirsim-key.pem
chmod 400 reservoirsim-key.pem
```

#### GCP
```bash
# Login and set project
gcloud auth login
gcloud config set project YOUR_PROJECT_ID

# Create SSH key (optional - gcloud manages this automatically)
ssh-keygen -t rsa -f ~/.ssh/reservoirsim-key -C ubuntu
```

---

## AWS Deployment

### Method 1: Automated Setup Script (Recommended)

The easiest way to deploy on AWS:

```bash
cd deploy/aws
./setup.sh
```

This interactive script will:
1. Check prerequisites
2. Prompt for configuration (region, instance type, etc.)
3. Deploy infrastructure using Terraform
4. Provide connection details

**Example Session:**
```bash
$ ./setup.sh
==========================================
ReservoirSim AWS Deployment Setup
==========================================

[INFO] Checking prerequisites...
[INFO] Prerequisites check passed!

Select deployment type:
1) Single EC2 instance (recommended for testing, small simulations)
2) AWS ParallelCluster (for large-scale parallel simulations)
3) Docker on EC2 (containerized deployment)
Enter choice [1-3]: 1

Enter AWS region (default: us-east-1): us-west-2
Enter SSH key pair name: reservoirsim-key

Select instance type:
1) c5.2xlarge  (8 vCPUs, 16 GB RAM) - Small simulations
2) c5.4xlarge  (16 vCPUs, 32 GB RAM) - Medium simulations [Recommended]
3) c5.9xlarge  (36 vCPUs, 72 GB RAM) - Large simulations
4) c5.18xlarge (72 vCPUs, 144 GB RAM) - Very large simulations
5) Custom instance type
Enter choice [1-5]: 2

[INFO] Using instance type: c5.4xlarge
[INFO] Deploying single EC2 instance with Terraform...
```

### Method 2: Manual Terraform Deployment

For more control over the deployment:

```bash
cd deploy/aws/terraform

# Create terraform.tfvars file
cat > terraform.tfvars << EOF
aws_region = "us-east-1"
key_name = "reservoirsim-key"
instance_type = "c5.4xlarge"
root_volume_size = 100
allowed_ssh_cidr = ["YOUR_IP/32"]  # Replace with your IP
EOF

# Initialize Terraform
terraform init

# Review the plan
terraform plan

# Apply the configuration
terraform apply
```

### Method 3: AWS ParallelCluster (Multi-Node HPC)

For large-scale simulations requiring multiple nodes:

```bash
# Install AWS ParallelCluster
pip3 install aws-parallelcluster

# Edit the configuration file
cd deploy/aws
nano parallelcluster-config.yaml

# Update the following fields:
# - SubnetId: Your VPC subnet ID
# - KeyName: Your SSH key name

# Create the cluster
pcluster create-cluster \
  --cluster-name reservoirsim-cluster \
  --cluster-configuration parallelcluster-config.yaml

# Check cluster status
pcluster describe-cluster --cluster-name reservoirsim-cluster

# Connect to head node
pcluster ssh --cluster-name reservoirsim-cluster
```

### Connecting to Your Instance

After deployment, connect via SSH:

```bash
# Get connection details
cd deploy/aws/terraform
terraform output

# SSH to instance
ssh -i ~/reservoirsim-key.pem ubuntu@<INSTANCE_PUBLIC_IP>
```

### Building ReservoirSim on AWS

Once connected to your instance:

```bash
# The startup script automatically installs dependencies
# Check installation status
sudo journalctl -u cloud-final -f

# Wait for installation to complete (5-10 minutes)

# Clone your repository (if not done automatically)
cd ~
git clone https://github.com/yourusername/ReservoirSim.git
cd ReservoirSim

# Build
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)

# Verify installation
./reservoirsim -help
```

### Running Simulations on AWS

```bash
cd ~/ReservoirSim/build

# Run a single-node simulation
mpirun -np 16 ./reservoirsim -c ../config/shale_reservoir.config

# For AWS ParallelCluster (multi-node):
# Create a SLURM job script
cat > run_simulation.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=reservoirsim
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00
#SBATCH --output=simulation_%j.out

module load openmpi

mpirun ./reservoirsim -c ../config/shale_reservoir.config
EOF

# Submit job
sbatch run_simulation.sh

# Check job status
squeue
```

### Data Management on AWS

The deployment automatically configures S3 integration:

```bash
# Environment variable is set automatically
echo $S3_BUCKET

# Upload input data to S3
aws s3 cp input_data.dat s3://$S3_BUCKET/inputs/

# Download from S3
aws s3 cp s3://$S3_BUCKET/inputs/input_data.dat ./

# Sync outputs to S3 (done automatically on shutdown)
~/sync_to_s3.sh

# Manual sync
aws s3 sync ~/simulations/output/ s3://$S3_BUCKET/outputs/
```

---

## Google Cloud Deployment

### Method 1: Automated Setup Script (Recommended)

```bash
cd deploy/gcp
./setup.sh
```

This interactive script will:
1. Check prerequisites and authenticate
2. Enable required GCP APIs
3. Prompt for configuration
4. Deploy infrastructure using Terraform

**Example Session:**
```bash
$ ./setup.sh
==========================================
ReservoirSim GCP Deployment Setup
==========================================

[INFO] Checking prerequisites...
[INFO] Prerequisites check passed!
[INFO] Using GCP project: my-project-id

Enter GCP region (default: us-central1): us-central1
Enter GCP zone (default: us-central1-a): us-central1-a

Select machine type:
1) c2-standard-8   (8 vCPUs, 32 GB RAM) - Small simulations
2) c2-standard-16  (16 vCPUs, 64 GB RAM) - Medium simulations [Recommended]
3) c2-standard-30  (30 vCPUs, 120 GB RAM) - Large simulations
4) c2-standard-60  (60 vCPUs, 240 GB RAM) - Very large simulations
5) n2-highcpu-32   (32 vCPUs, 32 GB RAM) - High CPU/memory ratio
6) Custom machine type
Enter choice [1-6]: 2

Use preemptible (spot) instances for cost savings? (yes/no): no
Number of local SSDs for high-performance temporary storage (0-8, default 0): 1
```

### Method 2: Manual Terraform Deployment

```bash
cd deploy/gcp/terraform

# Create terraform.tfvars
cat > terraform.tfvars << EOF
project_id = "your-gcp-project-id"
region = "us-central1"
zone = "us-central1-a"
machine_type = "c2-standard-16"
boot_disk_size = 100
local_ssd_count = 1
use_preemptible = false
EOF

# Initialize and apply
terraform init
terraform plan
terraform apply
```

### Connecting to Your GCE Instance

```bash
# Get connection details
cd deploy/gcp/terraform
terraform output

# SSH using gcloud (recommended)
gcloud compute ssh ubuntu@reservoirsim-compute --zone=us-central1-a

# Or use standard SSH
ssh ubuntu@<INSTANCE_EXTERNAL_IP>
```

### Building ReservoirSim on GCP

```bash
# Check startup script progress
sudo journalctl -u google-startup-scripts.service -f

# Wait for installation to complete (5-10 minutes)

# Clone repository (if not done automatically)
cd ~
git clone https://github.com/yourusername/ReservoirSim.git
cd ReservoirSim

# Build
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)

# Verify
./reservoirsim -help
```

### Running Simulations on GCP

```bash
cd ~/ReservoirSim/build

# Run simulation with all available cores
NCORES=$(nproc)
mpirun -np $NCORES ./reservoirsim -c ../config/geothermal.config

# Use local SSD for temporary storage (if configured)
mpirun -np $NCORES ./reservoirsim -c ../config/geothermal.config -o /mnt/localssd/output
```

### Data Management on GCP

```bash
# Environment variable is set automatically
echo $GCS_BUCKET

# Upload to Cloud Storage
gsutil cp input_data.dat gs://$GCS_BUCKET/inputs/

# Download from Cloud Storage
gsutil cp gs://$GCS_BUCKET/inputs/input_data.dat ./

# Sync outputs (done automatically on shutdown)
~/sync_to_gcs.sh

# Manual sync
gsutil -m rsync -r ~/simulations/output/ gs://$GCS_BUCKET/outputs/
```

---

## Docker Deployment

### Building the Docker Image

```bash
# Build locally
docker build -t reservoirsim:latest .

# Test locally
docker run -it reservoirsim:latest /bin/bash
```

### Running with Docker

```bash
# Run interactively
docker run -it -v $(pwd)/data:/data reservoirsim:latest

# Inside container
cd /app
mpirun -np 4 ./reservoirsim -c config/default.config

# Run simulation directly
docker run -v $(pwd)/config:/config -v $(pwd)/output:/output \
  reservoirsim:latest \
  mpirun -np 4 /app/reservoirsim -c /config/my_simulation.config -o /output/results
```

### Deploying Docker to AWS

```bash
# Build and save image
docker build -t reservoirsim:latest .
docker save reservoirsim:latest | gzip > reservoirsim-docker.tar.gz

# Upload to EC2 instance
scp -i reservoirsim-key.pem reservoirsim-docker.tar.gz ubuntu@<INSTANCE_IP>:~/

# On EC2 instance
ssh -i reservoirsim-key.pem ubuntu@<INSTANCE_IP>
docker load < reservoirsim-docker.tar.gz
docker run -it reservoirsim:latest
```

### Deploying Docker to GCP

```bash
# Tag for Google Container Registry
docker tag reservoirsim:latest gcr.io/YOUR_PROJECT_ID/reservoirsim:latest

# Push to GCR
docker push gcr.io/YOUR_PROJECT_ID/reservoirsim:latest

# On GCE instance
gcloud compute ssh ubuntu@reservoirsim-compute --zone=us-central1-a
docker pull gcr.io/YOUR_PROJECT_ID/reservoirsim:latest
docker run -it gcr.io/YOUR_PROJECT_ID/reservoirsim:latest
```

---

## Cost Optimization

### Instance Type Selection

Choose the right instance type for your simulation:

#### AWS Instance Types

| Instance Type | vCPUs | RAM (GB) | Cost/hour* | Best For |
|--------------|-------|----------|------------|----------|
| c5.2xlarge | 8 | 16 | $0.34 | Small simulations, testing |
| c5.4xlarge | 16 | 32 | $0.68 | Medium simulations |
| c5.9xlarge | 36 | 72 | $1.53 | Large simulations |
| c5.18xlarge | 72 | 144 | $3.06 | Very large simulations |
| c5n.18xlarge | 72 | 192 | $3.89 | Network-intensive parallel runs |

*Prices are approximate for us-east-1 on-demand. Check current pricing.

#### GCP Machine Types

| Machine Type | vCPUs | RAM (GB) | Cost/hour* | Best For |
|--------------|-------|----------|------------|----------|
| c2-standard-8 | 8 | 32 | $0.36 | Small-medium simulations |
| c2-standard-16 | 16 | 64 | $0.71 | Medium-large simulations |
| c2-standard-30 | 30 | 120 | $1.33 | Large simulations |
| c2-standard-60 | 60 | 240 | $2.67 | Very large simulations |
| n2-highcpu-32 | 32 | 32 | $0.95 | CPU-intensive workloads |

*Prices are approximate for us-central1. Check current pricing.

### Cost-Saving Strategies

1. **Use Spot/Preemptible Instances**
   - AWS Spot: Up to 90% savings
   - GCP Preemptible: Up to 80% savings
   - Best for fault-tolerant simulations with checkpointing

```bash
# AWS: Enable spot instances in terraform.tfvars
instance_market_options {
  market_type = "spot"
}

# GCP: Enable in terraform.tfvars
use_preemptible = true
```

2. **Auto-Stop Instances**

Create a CloudWatch/Cloud Scheduler rule to stop instances during off-hours:

```bash
# AWS - Stop instance at 6 PM
aws ec2 stop-instances --instance-ids i-1234567890abcdef0

# GCP - Stop instance
gcloud compute instances stop reservoirsim-compute --zone=us-central1-a
```

3. **Right-Size Your Instances**

Monitor CPU and memory usage:

```bash
# Install monitoring
htop
iotop

# AWS CloudWatch
aws cloudwatch get-metric-statistics ...

# GCP Monitoring
gcloud monitoring time-series list ...
```

4. **Data Storage Optimization**

```bash
# Use lifecycle policies for old data
# AWS S3
aws s3api put-bucket-lifecycle-configuration ...

# GCP Storage
gsutil lifecycle set lifecycle.json gs://your-bucket
```

5. **Reserved Instances**

For long-term use, purchase reserved instances (up to 75% savings):
- AWS: 1-year or 3-year reserved instances
- GCP: Committed use contracts

---

## Monitoring and Management

### Performance Monitoring

```bash
# CPU and memory usage
htop

# Disk I/O
iotop

# Network usage
nethogs

# MPI profiling
mpirun -np 16 --mca btl tcp,self --report-bindings ./reservoirsim ...
```

### Application Monitoring

```bash
# PETSc performance monitoring
mpirun -np 16 ./reservoirsim -c config/simulation.config \
  -log_view \
  -snes_monitor \
  -ksp_monitor \
  -ts_monitor

# Output to file
mpirun -np 16 ./reservoirsim -c config/simulation.config \
  -log_view :performance.log:ascii_info_detail
```

### Cloud-Native Monitoring

#### AWS CloudWatch

```bash
# Install CloudWatch agent
wget https://s3.amazonaws.com/amazoncloudwatch-agent/ubuntu/amd64/latest/amazon-cloudwatch-agent.deb
sudo dpkg -i amazon-cloudwatch-agent.deb

# Configure and start
sudo /opt/aws/amazon-cloudwatch-agent/bin/amazon-cloudwatch-agent-ctl -a fetch-config ...
```

#### GCP Cloud Monitoring

```bash
# Install monitoring agent
curl -sSO https://dl.google.com/cloudagents/add-google-cloud-ops-agent-repo.sh
sudo bash add-google-cloud-ops-agent-repo.sh --also-install
```

### Checkpointing

Enable checkpointing for long-running simulations:

```bash
# Add to your config file
[SIMULATION]
enable_checkpointing = true
checkpoint_frequency = 3600  # Every hour

# Restart from checkpoint
mpirun -np 16 ./reservoirsim -c config/simulation.config -restart checkpoint_0001.h5
```

---

## Troubleshooting

### Common Issues

#### 1. PETSc Not Found

```bash
# Check environment variables
echo $PETSC_DIR
echo $PETSC_ARCH
echo $PKG_CONFIG_PATH

# Reinstall if necessary
cd /opt/petsc-3.20.0
./configure ...
make all
```

#### 2. MPI Errors

```bash
# Check MPI installation
which mpirun
mpirun --version

# Test MPI
mpirun -np 4 hostname

# Use correct MPI launcher
mpirun -np 4 --mca btl tcp,self ./reservoirsim ...
```

#### 3. Out of Memory

```bash
# Check available memory
free -h

# Monitor during run
watch -n 1 free -h

# Solutions:
# - Use larger instance type
# - Reduce grid resolution
# - Enable out-of-core solvers in PETSc
```

#### 4. Slow Performance

```bash
# Check CPU frequency
lscpu | grep MHz

# Enable performance mode
sudo cpupower frequency-set --governor performance

# Check NUMA topology
numactl --hardware

# Optimize MPI binding
mpirun --bind-to core --map-by socket ...
```

#### 5. Network Issues (Multi-Node)

```bash
# Check connectivity
ping other-node

# Check security groups/firewall
# AWS: Ensure security group allows all traffic within the group
# GCP: Ensure firewall rules allow internal traffic

# Use appropriate MPI transport
mpirun --mca btl tcp,self ...  # For Ethernet
mpirun --mca btl openib,self ...  # For InfiniBand
```

### Getting Help

1. **Check Logs**
```bash
# Startup script logs
# AWS
sudo cat /var/log/cloud-init-output.log

# GCP
sudo journalctl -u google-startup-scripts.service
```

2. **Verify Installation**
```bash
# Run tests
cd ~/ReservoirSim/build
ctest -V
```

3. **Enable Debug Mode**
```bash
# Rebuild with debug symbols
cd ~/ReservoirSim/build
cmake .. -DCMAKE_BUILD_TYPE=Debug
make

# Run with debugger
mpirun -np 4 gdb --args ./reservoirsim -c config/simulation.config
```

---

## Example Workflows

### End-to-End AWS Workflow

```bash
# 1. Deploy infrastructure
cd deploy/aws
./setup.sh
# Select option 1, configure as needed

# 2. Wait for instance to initialize
# (5-10 minutes)

# 3. Connect
ssh -i reservoirsim-key.pem ubuntu@<INSTANCE_IP>

# 4. Verify installation
cd ~/ReservoirSim/build
./reservoirsim -help

# 5. Run simulation
mpirun -np 16 ./reservoirsim -c ../config/shale_reservoir.config -o ~/simulations/output/run1

# 6. Monitor progress
tail -f ~/simulations/output/run1.log

# 7. Sync results to S3
~/sync_to_s3.sh

# 8. Download results locally
exit
aws s3 sync s3://<BUCKET_NAME>/outputs/ ./local_results/

# 9. Cleanup (when done)
cd deploy/aws/terraform
terraform destroy
```

### End-to-End GCP Workflow

```bash
# 1. Deploy infrastructure
cd deploy/gcp
./setup.sh
# Select option 1, configure as needed

# 2. Connect
gcloud compute ssh ubuntu@reservoirsim-compute --zone=us-central1-a

# 3. Verify installation
cd ~/ReservoirSim/build
./reservoirsim -help

# 4. Run simulation
mpirun -np 16 ./reservoirsim -c ../config/geothermal.config

# 5. Sync to Cloud Storage
~/sync_to_gcs.sh

# 6. Download locally
exit
gsutil -m rsync -r gs://<BUCKET_NAME>/outputs/ ./local_results/

# 7. Cleanup
cd deploy/gcp/terraform
terraform destroy
```

---

## Best Practices

1. **Always use version control** - Keep your configuration files in Git
2. **Tag your cloud resources** - For cost tracking and organization
3. **Enable checkpointing** - For long-running simulations
4. **Monitor costs** - Set up billing alerts
5. **Use spot/preemptible instances** - For cost savings on fault-tolerant workloads
6. **Automate shutdown** - Stop instances when not in use
7. **Backup results** - Sync to S3/GCS regularly
8. **Test locally first** - Use Docker to test configurations before deploying
9. **Document your workflows** - Keep notes on what works for your simulations
10. **Security** - Restrict SSH access, use IAM roles, encrypt data

---

## Next Steps

- Review [Wave Physics Documentation](./Wave_Physics_Theory.md)
- Check [LEFM Theory](./LEFM_Theory.md) for fracture modeling
- Explore [example configurations](../config/)
- Set up automated testing pipelines
- Scale to multi-node clusters for massive simulations

---

## Support and Resources

- GitHub Issues: [Report problems or ask questions]
- Documentation: See `docs/` directory
- Examples: See `examples/` directory
- AWS Documentation: https://docs.aws.amazon.com/
- GCP Documentation: https://cloud.google.com/docs
- PETSc Documentation: https://petsc.org/release/docs/
