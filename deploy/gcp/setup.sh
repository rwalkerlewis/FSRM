#!/bin/bash
# Google Cloud Platform deployment setup script for ReservoirSim
# This script automates the deployment of ReservoirSim on GCP

set -e

echo "=========================================="
echo "ReservoirSim GCP Deployment Setup"
echo "=========================================="
echo ""

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored output
print_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Check prerequisites
print_info "Checking prerequisites..."

if ! command -v gcloud &> /dev/null; then
    print_error "Google Cloud SDK not found. Please install it first."
    echo "Visit: https://cloud.google.com/sdk/docs/install"
    exit 1
fi

if ! command -v terraform &> /dev/null; then
    print_error "Terraform not found. Please install it first."
    echo "Visit: https://www.terraform.io/downloads"
    exit 1
fi

print_info "Prerequisites check passed!"
echo ""

# Get GCP project ID
print_info "Getting GCP project information..."
PROJECT_ID=$(gcloud config get-value project 2>/dev/null)

if [ -z "$PROJECT_ID" ]; then
    print_error "No GCP project configured. Run: gcloud config set project YOUR_PROJECT_ID"
    exit 1
fi

print_info "Using GCP project: $PROJECT_ID"
echo ""

# Prompt for region and zone
read -p "Enter GCP region (default: us-central1): " REGION
REGION=${REGION:-us-central1}

read -p "Enter GCP zone (default: ${REGION}-a): " ZONE
ZONE=${ZONE:-${REGION}-a}

# Prompt for machine type
echo ""
echo "Select machine type:"
echo "1) c2-standard-8   (8 vCPUs, 32 GB RAM) - Small simulations"
echo "2) c2-standard-16  (16 vCPUs, 64 GB RAM) - Medium simulations [Recommended]"
echo "3) c2-standard-30  (30 vCPUs, 120 GB RAM) - Large simulations"
echo "4) c2-standard-60  (60 vCPUs, 240 GB RAM) - Very large simulations"
echo "5) n2-highcpu-32   (32 vCPUs, 32 GB RAM) - High CPU/memory ratio"
echo "6) Custom machine type"
read -p "Enter choice [1-6]: " MACHINE_CHOICE

case $MACHINE_CHOICE in
    1) MACHINE_TYPE="c2-standard-8";;
    2) MACHINE_TYPE="c2-standard-16";;
    3) MACHINE_TYPE="c2-standard-30";;
    4) MACHINE_TYPE="c2-standard-60";;
    5) MACHINE_TYPE="n2-highcpu-32";;
    6) read -p "Enter custom machine type: " MACHINE_TYPE;;
    *) MACHINE_TYPE="c2-standard-16";;
esac

print_info "Using machine type: $MACHINE_TYPE"
echo ""

# Ask about preemptible instances
read -p "Use preemptible (spot) instances for cost savings? (yes/no): " USE_PREEMPTIBLE
if [ "$USE_PREEMPTIBLE" = "yes" ]; then
    PREEMPTIBLE="true"
    print_warning "Preemptible instances can be terminated at any time. Make sure to checkpoint your simulations!"
else
    PREEMPTIBLE="false"
fi

# Ask about local SSDs
read -p "Number of local SSDs for high-performance temporary storage (0-8, default 0): " SSD_COUNT
SSD_COUNT=${SSD_COUNT:-0}

# Enable required APIs
print_info "Enabling required GCP APIs..."
gcloud services enable compute.googleapis.com
gcloud services enable storage-api.googleapis.com
gcloud services enable logging.googleapis.com
gcloud services enable monitoring.googleapis.com

print_info "APIs enabled successfully!"
echo ""

# Prompt for deployment type
echo "Select deployment type:"
echo "1) Single VM instance (recommended)"
echo "2) Docker on GCE"
echo "3) GKE cluster (for containerized deployment)"
read -p "Enter choice [1-3]: " DEPLOY_TYPE

case $DEPLOY_TYPE in
    1)
        print_info "Deploying single GCE instance with Terraform..."
        cd terraform
        
        # Create terraform.tfvars
        cat > terraform.tfvars << EOF
project_id = "$PROJECT_ID"
region = "$REGION"
zone = "$ZONE"
machine_type = "$MACHINE_TYPE"
use_preemptible = $PREEMPTIBLE
local_ssd_count = $SSD_COUNT
use_static_ip = false
EOF
        
        terraform init
        terraform plan
        
        read -p "Proceed with deployment? (yes/no): " CONFIRM
        if [ "$CONFIRM" = "yes" ]; then
            terraform apply -auto-approve
            
            print_info "Deployment complete!"
            print_info "Getting instance information..."
            
            INSTANCE_NAME=$(terraform output -raw instance_name)
            INSTANCE_IP=$(terraform output -raw instance_ip)
            BUCKET_NAME=$(terraform output -raw bucket_name)
            
            echo ""
            echo "=========================================="
            echo "Deployment Information"
            echo "=========================================="
            echo "Instance Name: $INSTANCE_NAME"
            echo "Instance IP: $INSTANCE_IP"
            echo "Cloud Storage Bucket: $BUCKET_NAME"
            echo "SSH Command: gcloud compute ssh ubuntu@$INSTANCE_NAME --zone=$ZONE"
            echo ""
            print_warning "Wait 5-10 minutes for instance initialization to complete"
            echo ""
            
            print_info "You can check startup script progress with:"
            echo "gcloud compute ssh ubuntu@$INSTANCE_NAME --zone=$ZONE --command='sudo journalctl -u google-startup-scripts.service'"
        fi
        ;;
        
    2)
        print_info "Deploying Docker-based GCE instance..."
        cd terraform
        
        cat > terraform.tfvars << EOF
project_id = "$PROJECT_ID"
region = "$REGION"
zone = "$ZONE"
machine_type = "$MACHINE_TYPE"
use_preemptible = $PREEMPTIBLE
local_ssd_count = $SSD_COUNT
EOF
        
        terraform init
        terraform apply -auto-approve
        
        INSTANCE_NAME=$(terraform output -raw instance_name)
        
        print_info "Building Docker image..."
        cd ../..
        docker build -t reservoirsim:latest .
        
        print_info "Pushing to Google Container Registry..."
        docker tag reservoirsim:latest gcr.io/$PROJECT_ID/reservoirsim:latest
        docker push gcr.io/$PROJECT_ID/reservoirsim:latest
        
        print_info "Pulling image on GCE instance..."
        gcloud compute ssh ubuntu@$INSTANCE_NAME --zone=$ZONE --command="docker pull gcr.io/$PROJECT_ID/reservoirsim:latest"
        
        print_info "Docker deployment complete!"
        echo "Connect with: gcloud compute ssh ubuntu@$INSTANCE_NAME --zone=$ZONE"
        ;;
        
    3)
        print_info "GKE cluster deployment..."
        print_warning "GKE deployment requires additional configuration."
        print_warning "Please refer to deploy/gcp/kubernetes/ directory for GKE deployment manifests."
        ;;
        
    *)
        print_error "Invalid choice"
        exit 1
        ;;
esac

print_info "Setup complete!"
