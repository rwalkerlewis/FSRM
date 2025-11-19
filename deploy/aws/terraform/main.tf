# Terraform configuration for ReservoirSim on AWS
# This creates an HPC-optimized EC2 instance or cluster for running simulations

terraform {
  required_version = ">= 1.0"
  required_providers {
    aws = {
      source  = "hashicorp/aws"
      version = "~> 5.0"
    }
  }
}

provider "aws" {
  region = var.aws_region
}

# VPC Configuration
resource "aws_vpc" "reservoirsim_vpc" {
  cidr_block           = "10.0.0.0/16"
  enable_dns_hostnames = true
  enable_dns_support   = true

  tags = {
    Name = "reservoirsim-vpc"
  }
}

resource "aws_subnet" "reservoirsim_subnet" {
  vpc_id                  = aws_vpc.reservoirsim_vpc.id
  cidr_block              = "10.0.1.0/24"
  availability_zone       = data.aws_availability_zones.available.names[0]
  map_public_ip_on_launch = true

  tags = {
    Name = "reservoirsim-subnet"
  }
}

resource "aws_internet_gateway" "reservoirsim_igw" {
  vpc_id = aws_vpc.reservoirsim_vpc.id

  tags = {
    Name = "reservoirsim-igw"
  }
}

resource "aws_route_table" "reservoirsim_rt" {
  vpc_id = aws_vpc.reservoirsim_vpc.id

  route {
    cidr_block = "0.0.0.0/0"
    gateway_id = aws_internet_gateway.reservoirsim_igw.id
  }

  tags = {
    Name = "reservoirsim-rt"
  }
}

resource "aws_route_table_association" "reservoirsim_rta" {
  subnet_id      = aws_subnet.reservoirsim_subnet.id
  route_table_id = aws_route_table.reservoirsim_rt.id
}

# Security Group
resource "aws_security_group" "reservoirsim_sg" {
  name        = "reservoirsim-sg"
  description = "Security group for ReservoirSim instances"
  vpc_id      = aws_vpc.reservoirsim_vpc.id

  # SSH access
  ingress {
    from_port   = 22
    to_port     = 22
    protocol    = "tcp"
    cidr_blocks = var.allowed_ssh_cidr
  }

  # MPI communication
  ingress {
    from_port = 0
    to_port   = 65535
    protocol  = "tcp"
    self      = true
  }

  ingress {
    from_port = 0
    to_port   = 65535
    protocol  = "udp"
    self      = true
  }

  # Outbound internet access
  egress {
    from_port   = 0
    to_port     = 0
    protocol    = "-1"
    cidr_blocks = ["0.0.0.0/0"]
  }

  tags = {
    Name = "reservoirsim-sg"
  }
}

# IAM Role for EC2 instances
resource "aws_iam_role" "reservoirsim_role" {
  name = "reservoirsim-ec2-role"

  assume_role_policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Action = "sts:AssumeRole"
        Effect = "Allow"
        Principal = {
          Service = "ec2.amazonaws.com"
        }
      }
    ]
  })
}

resource "aws_iam_role_policy_attachment" "ssm_policy" {
  role       = aws_iam_role.reservoirsim_role.name
  policy_arn = "arn:aws:iam::aws:policy/AmazonSSMManagedInstanceCore"
}

resource "aws_iam_role_policy" "s3_policy" {
  name = "reservoirsim-s3-policy"
  role = aws_iam_role.reservoirsim_role.id

  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Effect = "Allow"
        Action = [
          "s3:GetObject",
          "s3:PutObject",
          "s3:ListBucket"
        ]
        Resource = [
          aws_s3_bucket.reservoirsim_data.arn,
          "${aws_s3_bucket.reservoirsim_data.arn}/*"
        ]
      }
    ]
  })
}

resource "aws_iam_instance_profile" "reservoirsim_profile" {
  name = "reservoirsim-instance-profile"
  role = aws_iam_role.reservoirsim_role.name
}

# S3 Bucket for data storage
resource "aws_s3_bucket" "reservoirsim_data" {
  bucket = "${var.project_name}-data-${random_string.suffix.result}"

  tags = {
    Name = "reservoirsim-data"
  }
}

resource "aws_s3_bucket_versioning" "reservoirsim_data_versioning" {
  bucket = aws_s3_bucket.reservoirsim_data.id

  versioning_configuration {
    status = "Enabled"
  }
}

resource "random_string" "suffix" {
  length  = 8
  special = false
  upper   = false
}

# EC2 Instance for single-node deployment
resource "aws_instance" "reservoirsim_instance" {
  count = var.use_parallel_cluster ? 0 : 1

  ami                    = data.aws_ami.ubuntu.id
  instance_type          = var.instance_type
  key_name               = var.key_name
  subnet_id              = aws_subnet.reservoirsim_subnet.id
  vpc_security_group_ids = [aws_security_group.reservoirsim_sg.id]
  iam_instance_profile   = aws_iam_instance_profile.reservoirsim_profile.name

  root_block_device {
    volume_size = var.root_volume_size
    volume_type = "gp3"
    iops        = 3000
    throughput  = 125
  }

  user_data = templatefile("${path.module}/user_data.sh", {
    s3_bucket = aws_s3_bucket.reservoirsim_data.id
  })

  tags = {
    Name = "reservoirsim-compute"
  }
}

# Data source for Ubuntu AMI
data "aws_ami" "ubuntu" {
  most_recent = true
  owners      = ["099720109477"] # Canonical

  filter {
    name   = "name"
    values = ["ubuntu/images/hvm-ssd/ubuntu-jammy-22.04-amd64-server-*"]
  }

  filter {
    name   = "virtualization-type"
    values = ["hvm"]
  }
}

data "aws_availability_zones" "available" {
  state = "available"
}

# Elastic IP
resource "aws_eip" "reservoirsim_eip" {
  count    = var.use_parallel_cluster ? 0 : 1
  instance = aws_instance.reservoirsim_instance[0].id
  domain   = "vpc"

  tags = {
    Name = "reservoirsim-eip"
  }
}

# Outputs
output "instance_public_ip" {
  value       = var.use_parallel_cluster ? null : aws_eip.reservoirsim_eip[0].public_ip
  description = "Public IP address of the ReservoirSim instance"
}

output "s3_bucket_name" {
  value       = aws_s3_bucket.reservoirsim_data.id
  description = "Name of the S3 bucket for data storage"
}

output "ssh_command" {
  value       = var.use_parallel_cluster ? null : "ssh -i ${var.key_name}.pem ubuntu@${aws_eip.reservoirsim_eip[0].public_ip}"
  description = "SSH command to connect to the instance"
}
