# Variables for GCP deployment

variable "project_id" {
  description = "GCP Project ID"
  type        = string
}

variable "region" {
  description = "GCP region for deployment"
  type        = string
  default     = "us-central1"
}

variable "zone" {
  description = "GCP zone for deployment"
  type        = string
  default     = "us-central1-a"
}

variable "machine_type" {
  description = "GCE machine type (use compute-optimized for HPC)"
  type        = string
  default     = "c2-standard-16" # 16 vCPUs, 64 GB RAM
  # Options:
  # c2-standard-4   - 4 vCPUs, 16 GB RAM
  # c2-standard-8   - 8 vCPUs, 32 GB RAM
  # c2-standard-16  - 16 vCPUs, 64 GB RAM
  # c2-standard-30  - 30 vCPUs, 120 GB RAM
  # c2-standard-60  - 60 vCPUs, 240 GB RAM
  # c2d-standard-112 - 112 vCPUs, 448 GB RAM (3rd gen)
  # n2-highcpu-32   - 32 vCPUs, 32 GB RAM (high CPU)
}

variable "boot_disk_size" {
  description = "Size of boot disk in GB"
  type        = number
  default     = 100
}

variable "local_ssd_count" {
  description = "Number of local SSDs (375 GB each, for temporary high-performance storage)"
  type        = number
  default     = 0
}

variable "ssh_key" {
  description = "SSH public key for instance access"
  type        = string
  default     = ""
}

variable "allowed_ssh_ranges" {
  description = "CIDR ranges allowed to SSH to the instance"
  type        = list(string)
  default     = ["0.0.0.0/0"] # WARNING: Restrict this in production!
}

variable "use_preemptible" {
  description = "Use preemptible (spot) instances for cost savings"
  type        = bool
  default     = false
}

variable "use_static_ip" {
  description = "Use static IP address"
  type        = bool
  default     = false
}
