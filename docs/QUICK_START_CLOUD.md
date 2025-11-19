# ReservoirSim Cloud Quick Start Guide

Get ReservoirSim running on AWS or Google Cloud in 15 minutes.

## Prerequisites

- AWS or GCP account
- AWS CLI or Google Cloud SDK installed locally
- Basic familiarity with cloud computing

## AWS Quick Start

### Step 1: Install Prerequisites

```bash
# Install AWS CLI (if not already installed)
# macOS
brew install awscli

# Linux
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip
sudo ./aws/install

# Install Terraform
brew install terraform  # macOS
# OR
wget https://releases.hashicorp.com/terraform/1.6.0/terraform_1.6.0_linux_amd64.zip
unzip terraform_*.zip
sudo mv terraform /usr/local/bin/
```

### Step 2: Configure AWS

```bash
# Set up credentials
aws configure
# Enter your Access Key ID
# Enter your Secret Access Key
# Default region: us-east-1
# Output format: json

# Create SSH key pair
aws ec2 create-key-pair --key-name reservoirsim-key \
    --query 'KeyMaterial' --output text > ~/.ssh/reservoirsim-key.pem
chmod 400 ~/.ssh/reservoirsim-key.pem
```

### Step 3: Deploy

```bash
# Clone the repository
git clone https://github.com/yourusername/ReservoirSim.git
cd ReservoirSim

# Run the automated setup
cd deploy/aws
./setup.sh
```

Follow the prompts:
- Select deployment type: `1` (Single EC2 instance)
- AWS region: Press Enter for default (us-east-1)
- SSH key name: `reservoirsim-key`
- Instance type: `2` (c5.4xlarge - recommended)

Wait 2-3 minutes for Terraform to provision infrastructure.

### Step 4: Connect and Test

```bash
# Get the instance IP from output
export INSTANCE_IP=$(cd terraform && terraform output -raw instance_public_ip)

# SSH to instance (wait 5-10 minutes for initialization to complete)
ssh -i ~/.ssh/reservoirsim-key.pem ubuntu@$INSTANCE_IP

# Check installation status
sudo tail -f /var/log/cloud-init-output.log
# Wait until you see "Instance initialization completed"
# Press Ctrl+C when done

# Build ReservoirSim
cd ReservoirSim
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)

# Run a test
./reservoirsim -c ../config/default.config
```

### Step 5: Run Your First Simulation

```bash
# Run a shale reservoir simulation with 16 cores
mpirun -np 16 ./reservoirsim -c ../config/shale_reservoir.config

# Monitor progress
tail -f output.log
```

### Step 6: Cleanup (When Done)

```bash
# Exit from instance
exit

# Destroy infrastructure (to stop charges)
cd deploy/aws/terraform
terraform destroy
```

**Estimated Cost**: ~$0.68/hour for c5.4xlarge instance

---

## Google Cloud Quick Start

### Step 1: Install Prerequisites

```bash
# Install Google Cloud SDK
# macOS
brew install --cask google-cloud-sdk

# Linux
curl https://sdk.cloud.google.com | bash
exec -l $SHELL

# Install Terraform
brew install terraform  # macOS
# OR
wget https://releases.hashicorp.com/terraform/1.6.0/terraform_1.6.0_linux_amd64.zip
unzip terraform_*.zip
sudo mv terraform /usr/local/bin/
```

### Step 2: Configure GCP

```bash
# Login to Google Cloud
gcloud auth login

# Set your project
gcloud config set project YOUR_PROJECT_ID

# Enable billing (if not already enabled)
gcloud beta billing projects link YOUR_PROJECT_ID \
    --billing-account=YOUR_BILLING_ACCOUNT_ID
```

### Step 3: Deploy

```bash
# Clone the repository
git clone https://github.com/yourusername/ReservoirSim.git
cd ReservoirSim

# Run the automated setup
cd deploy/gcp
./setup.sh
```

Follow the prompts:
- GCP region: Press Enter for default (us-central1)
- GCP zone: Press Enter for default (us-central1-a)
- Machine type: `2` (c2-standard-16 - recommended)
- Use preemptible: `no` (for first run)
- Local SSDs: `0` (for first run)

Wait 2-3 minutes for deployment.

### Step 4: Connect and Test

```bash
# Connect via gcloud (easiest)
gcloud compute ssh ubuntu@reservoirsim-compute --zone=us-central1-a

# Check installation status
sudo journalctl -u google-startup-scripts.service -f
# Wait until initialization completes
# Press Ctrl+C when done

# Build ReservoirSim
cd ReservoirSim
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)

# Run a test
./reservoirsim -help
```

### Step 5: Run Your First Simulation

```bash
# Run a geothermal simulation
mpirun -np 16 ./reservoirsim -c ../config/geothermal.config

# Monitor progress
watch -n 1 'tail -20 output.log'
```

### Step 6: Cleanup (When Done)

```bash
# Exit from instance
exit

# Destroy infrastructure
cd deploy/gcp/terraform
terraform destroy
```

**Estimated Cost**: ~$0.71/hour for c2-standard-16 instance

---

## Docker Quick Start (Local or Cloud)

### Run Locally with Docker

```bash
# Build the Docker image
docker build -t reservoirsim:latest .

# Run interactively
docker run -it reservoirsim:latest /bin/bash

# Inside container
cd /app
mpirun -np 4 ./reservoirsim -c config/default.config
```

### Deploy Docker to AWS/GCP

Already deployed an instance? Add Docker:

```bash
# On your cloud instance
docker pull yourusername/reservoirsim:latest
docker run -it -v ~/data:/data reservoirsim:latest
```

---

## Common First Simulations

### 1. Default Test (2 minutes)

Quick test to verify everything works:

```bash
mpirun -np 4 ./reservoirsim -c ../config/default.config
```

### 2. Shale Reservoir (10 minutes)

Hydraulic fracturing simulation:

```bash
mpirun -np 16 ./reservoirsim -c ../config/shale_reservoir.config
```

### 3. Geothermal System (15 minutes)

Enhanced geothermal system with thermal effects:

```bash
mpirun -np 16 ./reservoirsim -c ../config/geothermal.config
```

### 4. CO2 Storage (20 minutes)

Carbon sequestration simulation:

```bash
mpirun -np 16 ./reservoirsim -c ../config/co2_storage.config
```

---

## Tips for First-Time Users

### Monitoring Your Simulation

```bash
# Watch CPU usage
htop

# Monitor output
tail -f output.log

# Check disk space
df -h
```

### Transferring Results

#### AWS
```bash
# Sync to S3 (automatically configured)
~/sync_to_s3.sh

# Download locally
aws s3 sync s3://YOUR-BUCKET/outputs/ ./local_results/
```

#### GCP
```bash
# Sync to Cloud Storage
~/sync_to_gcs.sh

# Download locally
gsutil -m rsync -r gs://YOUR-BUCKET/outputs/ ./local_results/
```

### Visualizing Results

```bash
# VTK files - use ParaView
paraview output/*.vtu

# HDF5 files - use Python
python3
>>> import h5py
>>> f = h5py.File('output.h5', 'r')
>>> print(list(f.keys()))
```

### Saving Costs

1. **Stop instance when not in use**
```bash
# AWS
aws ec2 stop-instances --instance-ids i-XXXXX

# GCP
gcloud compute instances stop reservoirsim-compute --zone=us-central1-a
```

2. **Use spot/preemptible instances** (up to 80% savings)
   - Re-run setup script and select preemptible option
   - Enable checkpointing in your simulation config

3. **Choose right instance size**
   - Small test: c5.2xlarge (AWS) or c2-standard-8 (GCP)
   - Production: c5.4xlarge (AWS) or c2-standard-16 (GCP)

---

## Troubleshooting

### "PETSc not found" error

```bash
# Check if installation completed
echo $PETSC_DIR
source /etc/profile.d/petsc.sh

# Or reinstall
cd /workspace/deploy/scripts
./install_petsc.sh
```

### "Out of memory" error

- Use larger instance type
- Reduce grid resolution in config file
- Enable disk-based solvers

### Slow performance

```bash
# Check if using all cores
htop

# Verify MPI is working
mpirun -np 4 hostname

# Should see 4 different responses
```

### Can't connect to instance

- Wait 5-10 minutes for initialization
- Check security group allows SSH (port 22)
- Verify key permissions: `chmod 400 key.pem`

---

## Next Steps

1. **Explore configurations**: Check `config/` directory for more examples
2. **Read documentation**: See `docs/CLOUD_DEPLOYMENT.md` for advanced usage
3. **Customize simulations**: Edit config files for your specific case
4. **Scale up**: Try AWS ParallelCluster for multi-node simulations
5. **Optimize costs**: Set up automated shutdown and monitoring

---

## Getting Help

- **Documentation**: `docs/CLOUD_DEPLOYMENT.md`
- **Examples**: `examples/` directory
- **Issues**: GitHub Issues page
- **Community**: [Your community forum/chat]

---

## Approximate Costs

### AWS (us-east-1)
- c5.2xlarge: $0.34/hour (~$8/day)
- c5.4xlarge: $0.68/hour (~$16/day)
- c5.9xlarge: $1.53/hour (~$37/day)
- Storage (S3): ~$0.023/GB/month

### GCP (us-central1)
- c2-standard-8: $0.36/hour (~$9/day)
- c2-standard-16: $0.71/hour (~$17/day)
- c2-standard-30: $1.33/hour (~$32/day)
- Storage (GCS): ~$0.020/GB/month

**Note**: Prices are approximate. Check current pricing on cloud provider websites.

Use spot/preemptible instances for 70-80% savings!

---

**You're all set!** Start running reservoir simulations in the cloud. 🚀
