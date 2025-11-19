#!/bin/bash
# AWS deployment setup script for ReservoirSim
# This script automates the deployment of ReservoirSim on AWS

set -e

echo "=========================================="
echo "ReservoirSim AWS Deployment Setup"
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

if ! command -v aws &> /dev/null; then
    print_error "AWS CLI not found. Please install it first."
    echo "Visit: https://aws.amazon.com/cli/"
    exit 1
fi

if ! command -v terraform &> /dev/null; then
    print_error "Terraform not found. Please install it first."
    echo "Visit: https://www.terraform.io/downloads"
    exit 1
fi

print_info "Prerequisites check passed!"
echo ""

# Prompt for deployment type
echo "Select deployment type:"
echo "1) Single EC2 instance (recommended for testing, small simulations)"
echo "2) AWS ParallelCluster (for large-scale parallel simulations)"
echo "3) Docker on EC2 (containerized deployment)"
read -p "Enter choice [1-3]: " DEPLOY_TYPE

# Prompt for AWS region
read -p "Enter AWS region (default: us-east-1): " AWS_REGION
AWS_REGION=${AWS_REGION:-us-east-1}

# Prompt for SSH key name
print_warning "Make sure you have created an SSH key pair in AWS EC2!"
read -p "Enter SSH key pair name: " KEY_NAME

if [ -z "$KEY_NAME" ]; then
    print_error "SSH key name is required!"
    exit 1
fi

# Prompt for instance type
echo ""
echo "Select instance type:"
echo "1) c5.2xlarge  (8 vCPUs, 16 GB RAM) - Small simulations"
echo "2) c5.4xlarge  (16 vCPUs, 32 GB RAM) - Medium simulations [Recommended]"
echo "3) c5.9xlarge  (36 vCPUs, 72 GB RAM) - Large simulations"
echo "4) c5.18xlarge (72 vCPUs, 144 GB RAM) - Very large simulations"
echo "5) Custom instance type"
read -p "Enter choice [1-5]: " INSTANCE_CHOICE

case $INSTANCE_CHOICE in
    1) INSTANCE_TYPE="c5.2xlarge";;
    2) INSTANCE_TYPE="c5.4xlarge";;
    3) INSTANCE_TYPE="c5.9xlarge";;
    4) INSTANCE_TYPE="c5.18xlarge";;
    5) read -p "Enter custom instance type: " INSTANCE_TYPE;;
    *) INSTANCE_TYPE="c5.4xlarge";;
esac

print_info "Using instance type: $INSTANCE_TYPE"
echo ""

# Deploy based on selection
case $DEPLOY_TYPE in
    1)
        print_info "Deploying single EC2 instance with Terraform..."
        cd terraform
        
        # Create terraform.tfvars
        cat > terraform.tfvars << EOF
aws_region = "$AWS_REGION"
key_name = "$KEY_NAME"
instance_type = "$INSTANCE_TYPE"
use_parallel_cluster = false
EOF
        
        terraform init
        terraform plan
        
        read -p "Proceed with deployment? (yes/no): " CONFIRM
        if [ "$CONFIRM" = "yes" ]; then
            terraform apply -auto-approve
            
            print_info "Deployment complete!"
            print_info "Getting instance information..."
            
            INSTANCE_IP=$(terraform output -raw instance_public_ip)
            S3_BUCKET=$(terraform output -raw s3_bucket_name)
            
            echo ""
            echo "=========================================="
            echo "Deployment Information"
            echo "=========================================="
            echo "Instance IP: $INSTANCE_IP"
            echo "S3 Bucket: $S3_BUCKET"
            echo "SSH Command: ssh -i $KEY_NAME.pem ubuntu@$INSTANCE_IP"
            echo ""
            print_warning "Wait 5-10 minutes for instance initialization to complete"
            echo ""
        fi
        ;;
        
    2)
        print_info "Setting up AWS ParallelCluster..."
        
        if ! command -v pcluster &> /dev/null; then
            print_info "Installing AWS ParallelCluster CLI..."
            pip3 install aws-parallelcluster
        fi
        
        print_warning "Please edit deploy/aws/parallelcluster-config.yaml with your subnet IDs"
        print_warning "Then run: pcluster create-cluster --cluster-name reservoirsim-cluster --cluster-configuration parallelcluster-config.yaml"
        ;;
        
    3)
        print_info "Deploying Docker-based EC2 instance..."
        cd terraform
        
        cat > terraform.tfvars << EOF
aws_region = "$AWS_REGION"
key_name = "$KEY_NAME"
instance_type = "$INSTANCE_TYPE"
use_parallel_cluster = false
EOF
        
        terraform init
        terraform plan
        terraform apply -auto-approve
        
        INSTANCE_IP=$(terraform output -raw instance_public_ip)
        
        print_info "Building and deploying Docker container..."
        
        # Build Docker image locally
        cd ../..
        docker build -t reservoirsim:latest .
        
        # Save and upload to instance
        docker save reservoirsim:latest | gzip > reservoirsim-docker.tar.gz
        
        print_info "Uploading Docker image to EC2 instance..."
        scp -i "$KEY_NAME.pem" reservoirsim-docker.tar.gz ubuntu@$INSTANCE_IP:/home/ubuntu/
        
        print_info "Loading Docker image on EC2 instance..."
        ssh -i "$KEY_NAME.pem" ubuntu@$INSTANCE_IP "docker load < reservoirsim-docker.tar.gz"
        
        print_info "Docker deployment complete!"
        echo "Connect with: ssh -i $KEY_NAME.pem ubuntu@$INSTANCE_IP"
        ;;
        
    *)
        print_error "Invalid choice"
        exit 1
        ;;
esac

print_info "Setup complete!"
