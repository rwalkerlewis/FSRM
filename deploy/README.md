# FSRM Deployment Tools

This directory contains deployment scripts and configurations for running FSRM on cloud platforms.

## Directory Structure

```
deploy/
├── aws/                          # AWS deployment files
│   ├── terraform/                # Terraform infrastructure as code
│   │   ├── main.tf              # Main Terraform configuration
│   │   ├── variables.tf         # Input variables
│   │   └── user_data.sh         # Instance initialization script
│   ├── parallelcluster-config.yaml  # AWS ParallelCluster config
│   └── setup.sh                 # Automated setup script
│
├── gcp/                          # Google Cloud deployment files
│   ├── terraform/                # Terraform infrastructure as code
│   │   ├── main.tf              # Main Terraform configuration
│   │   ├── variables.tf         # Input variables
│   │   └── startup-script.sh    # Instance initialization script
│   └── setup.sh                 # Automated setup script
│
└── scripts/                      # Utility scripts
    ├── install_petsc.sh         # PETSc installation script
    ├── build_fsrm.sh    # Build script
    └── run_example.sh           # Example runner

```

## Quick Start

### AWS Deployment

```bash
cd deploy/aws
./setup.sh
```

### Google Cloud Deployment

```bash
cd deploy/gcp
./setup.sh
```

## Documentation

- **Quick Start Guide**: `../docs/QUICK_START_CLOUD.md`
- **Comprehensive Guide**: `../docs/CLOUD_DEPLOYMENT.md`
- **Main README**: `../README.md`

## Deployment Options

### 1. Single Instance Deployment

**Best for**: Testing, development, small-to-medium simulations

**Cost**: $0.34-$3.00/hour depending on instance size

**Setup time**: 5-10 minutes

```bash
# AWS
cd deploy/aws && ./setup.sh
# Select option 1

# GCP
cd deploy/gcp && ./setup.sh
# Select option 1
```

### 2. AWS ParallelCluster (Multi-Node)

**Best for**: Large-scale parallel simulations, production workloads

**Cost**: Variable, scales with number of nodes

**Setup time**: 15-30 minutes

```bash
cd deploy/aws
# Edit parallelcluster-config.yaml with your settings
pcluster create-cluster --cluster-name fsrm-cluster \
  --cluster-configuration parallelcluster-config.yaml
```

### 3. Docker Deployment

**Best for**: Reproducibility, development, CI/CD

**Cost**: Same as underlying instance

**Setup time**: 5 minutes

```bash
# Build locally
docker build -t fsrm:latest .

# Or use docker-compose
docker-compose up -d fsrm

# Run simulation
docker exec -it fsrm mpirun -np 4 /app/fsrm -c /config/default.config
```

## Prerequisites

### Local Machine

1. **Cloud CLI**:
   - AWS: [AWS CLI](https://aws.amazon.com/cli/)
   - GCP: [gcloud SDK](https://cloud.google.com/sdk)

2. **Infrastructure Tools**:
   - [Terraform](https://www.terraform.io/downloads)
   - [Docker](https://docs.docker.com/get-docker/) (optional)

3. **Authentication**:
   ```bash
   # AWS
   aws configure
   
   # GCP
   gcloud auth login
   gcloud config set project YOUR_PROJECT_ID
   ```

## Terraform Deployments

### AWS

```bash
cd deploy/aws/terraform

# Create configuration
cat > terraform.tfvars << EOF
aws_region = "us-east-1"
key_name = "your-key-name"
instance_type = "c5.4xlarge"
EOF

# Deploy
terraform init
terraform plan
terraform apply

# Get outputs
terraform output

# Cleanup
terraform destroy
```

### GCP

```bash
cd deploy/gcp/terraform

# Create configuration
cat > terraform.tfvars << EOF
project_id = "your-project-id"
region = "us-central1"
zone = "us-central1-a"
machine_type = "c2-standard-16"
EOF

# Deploy
terraform init
terraform plan
terraform apply

# Get outputs
terraform output

# Cleanup
terraform destroy
```

## Utility Scripts

### Install PETSc

Installs PETSc on any Linux system:

```bash
cd deploy/scripts
./install_petsc.sh
```

Customize installation:

```bash
PETSC_VERSION=3.20.0 INSTALL_DIR=/opt ./install_petsc.sh
```

### Build FSRM

Builds FSRM from source:

```bash
cd deploy/scripts
./build_fsrm.sh
```

Customize build:

```bash
BUILD_TYPE=Debug ./build_fsrm.sh
```

### Run Examples

Interactive example runner:

```bash
cd deploy/scripts
./run_example.sh
```

Or specify parameters:

```bash
NPROCS=16 OUTPUT_DIR=./results ./run_example.sh
```

## Instance Type Selection

### AWS Recommendations

| Use Case | Instance Type | vCPUs | RAM | Cost/hr* |
|----------|---------------|-------|-----|----------|
| Testing | t3.xlarge | 4 | 16 GB | $0.17 |
| Small sims | c5.2xlarge | 8 | 16 GB | $0.34 |
| **Recommended** | **c5.4xlarge** | **16** | **32 GB** | **$0.68** |
| Large sims | c5.9xlarge | 36 | 72 GB | $1.53 |
| Very large | c5.18xlarge | 72 | 144 GB | $3.06 |
| HPC | hpc6a.48xlarge | 96 | 384 GB | $2.88 |

### GCP Recommendations

| Use Case | Machine Type | vCPUs | RAM | Cost/hr* |
|----------|--------------|-------|-----|----------|
| Testing | n2-standard-4 | 4 | 16 GB | $0.19 |
| Small sims | c2-standard-8 | 8 | 32 GB | $0.36 |
| **Recommended** | **c2-standard-16** | **16** | **64 GB** | **$0.71** |
| Large sims | c2-standard-30 | 30 | 120 GB | $1.33 |
| Very large | c2-standard-60 | 60 | 240 GB | $2.67 |
| High CPU | n2-highcpu-32 | 32 | 32 GB | $0.95 |

*Prices are approximate on-demand rates. Check current pricing.

## Cost Optimization

### Use Spot/Preemptible Instances

Save 70-80% with interruptible instances:

```bash
# AWS Terraform
instance_market_options {
  market_type = "spot"
}

# GCP Terraform
use_preemptible = true
```

### Auto-Shutdown

Set up automatic shutdown during off-hours:

```bash
# AWS - Stop at 6 PM, start at 8 AM
aws ec2 stop-instances --instance-ids i-xxx
aws ec2 start-instances --instance-ids i-xxx

# GCP
gcloud compute instances stop INSTANCE_NAME
gcloud compute instances start INSTANCE_NAME
```

### Right-Size Instances

Monitor usage and adjust:

```bash
# Monitor CPU/memory
htop

# Check if you're using all resources
# If CPU < 80%, consider smaller instance
```

## Data Management

### AWS S3

```bash
# Upload inputs
aws s3 cp input.dat s3://$S3_BUCKET/inputs/

# Download results
aws s3 sync s3://$S3_BUCKET/outputs/ ./local_results/

# Automatic sync on shutdown (configured by deployment)
~/sync_to_s3.sh
```

### GCP Cloud Storage

```bash
# Upload inputs
gsutil cp input.dat gs://$GCS_BUCKET/inputs/

# Download results
gsutil -m rsync -r gs://$GCS_BUCKET/outputs/ ./local_results/

# Automatic sync on shutdown (configured by deployment)
~/sync_to_gcs.sh
```

## Monitoring

### System Monitoring

```bash
# CPU and memory
htop

# Disk I/O
iotop

# Network
nethogs

# Disk usage
df -h
```

### Application Monitoring

```bash
# PETSc performance logging
mpirun -np 16 ./fsrm -c config.file -log_view

# Monitor convergence
mpirun -np 16 ./fsrm -c config.file -snes_monitor -ksp_monitor -ts_monitor
```

### Cloud Monitoring

```bash
# AWS CloudWatch
aws cloudwatch get-metric-statistics ...

# GCP Cloud Monitoring
gcloud monitoring time-series list ...
```

## Troubleshooting

### Instance won't connect

```bash
# Wait for initialization (5-10 minutes)
# Check security groups allow SSH

# AWS
aws ec2 describe-instance-status --instance-ids i-xxx

# GCP
gcloud compute instances describe INSTANCE_NAME
```

### PETSc not found

```bash
# Check environment
echo $PETSC_DIR
source /etc/profile.d/petsc.sh

# Reinstall
cd deploy/scripts && ./install_petsc.sh
```

### Build failures

```bash
# Check dependencies
cmake --version
mpicc --version
pkg-config --modversion PETSc

# Clean rebuild
rm -rf build && mkdir build && cd build
cmake .. && make clean && make -j$(nproc)
```

### Out of memory

Solutions:
1. Use larger instance type
2. Reduce grid resolution
3. Enable out-of-core solvers
4. Add swap space

## Getting Help

- **Documentation**: `../docs/CLOUD_DEPLOYMENT.md`
- **Quick Start**: `../docs/QUICK_START_CLOUD.md`
- **Issues**: GitHub Issues
- **Examples**: `../examples/`

## Contributing

Improvements to deployment scripts are welcome! Please test thoroughly before submitting PRs.

## License

This project is licensed under the [MIT License](../LICENSE).
