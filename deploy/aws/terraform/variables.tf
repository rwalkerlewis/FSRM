# Variables for AWS deployment

variable "aws_region" {
  description = "AWS region for deployment"
  type        = string
  default     = "us-east-1"
}

variable "project_name" {
  description = "Project name for resource naming"
  type        = string
  default     = "fsrm"
}

variable "instance_type" {
  description = "EC2 instance type (use compute-optimized for HPC)"
  type        = string
  default     = "c5.4xlarge" # 16 vCPUs, 32 GB RAM
  # Options:
  # c5.2xlarge  - 8 vCPUs, 16 GB RAM
  # c5.4xlarge  - 16 vCPUs, 32 GB RAM
  # c5.9xlarge  - 36 vCPUs, 72 GB RAM
  # c5.18xlarge - 72 vCPUs, 144 GB RAM
  # c5n.18xlarge - 72 vCPUs, 192 GB RAM, 100 Gbps network
  # hpc6a.48xlarge - 96 vCPUs, 384 GB RAM (HPC optimized)
}

variable "key_name" {
  description = "Name of the SSH key pair (must exist in AWS)"
  type        = string
}

variable "allowed_ssh_cidr" {
  description = "CIDR blocks allowed to SSH to the instance"
  type        = list(string)
  default     = ["0.0.0.0/0"] # WARNING: Restrict this in production!
}

variable "root_volume_size" {
  description = "Size of root volume in GB"
  type        = number
  default     = 100
}

variable "use_parallel_cluster" {
  description = "Whether to use AWS ParallelCluster for multi-node deployment"
  type        = bool
  default     = false
}
