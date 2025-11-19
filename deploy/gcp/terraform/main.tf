# Terraform configuration for ReservoirSim on Google Cloud Platform
# This creates a compute-optimized VM instance for running simulations

terraform {
  required_version = ">= 1.0"
  required_providers {
    google = {
      source  = "hashicorp/google"
      version = "~> 5.0"
    }
  }
}

provider "google" {
  project = var.project_id
  region  = var.region
  zone    = var.zone
}

# VPC Network
resource "google_compute_network" "reservoirsim_network" {
  name                    = "reservoirsim-network"
  auto_create_subnetworks = false
}

resource "google_compute_subnetwork" "reservoirsim_subnet" {
  name          = "reservoirsim-subnet"
  ip_cidr_range = "10.0.1.0/24"
  region        = var.region
  network       = google_compute_network.reservoirsim_network.id
}

# Firewall rules
resource "google_compute_firewall" "reservoirsim_ssh" {
  name    = "reservoirsim-allow-ssh"
  network = google_compute_network.reservoirsim_network.name

  allow {
    protocol = "tcp"
    ports    = ["22"]
  }

  source_ranges = var.allowed_ssh_ranges
  target_tags   = ["reservoirsim"]
}

resource "google_compute_firewall" "reservoirsim_internal" {
  name    = "reservoirsim-allow-internal"
  network = google_compute_network.reservoirsim_network.name

  allow {
    protocol = "tcp"
    ports    = ["0-65535"]
  }

  allow {
    protocol = "udp"
    ports    = ["0-65535"]
  }

  allow {
    protocol = "icmp"
  }

  source_tags = ["reservoirsim"]
  target_tags = ["reservoirsim"]
}

# Service Account
resource "google_service_account" "reservoirsim_sa" {
  account_id   = "reservoirsim-compute"
  display_name = "ReservoirSim Compute Service Account"
}

resource "google_project_iam_member" "reservoirsim_sa_storage" {
  project = var.project_id
  role    = "roles/storage.objectAdmin"
  member  = "serviceAccount:${google_service_account.reservoirsim_sa.email}"
}

resource "google_project_iam_member" "reservoirsim_sa_logging" {
  project = var.project_id
  role    = "roles/logging.logWriter"
  member  = "serviceAccount:${google_service_account.reservoirsim_sa.email}"
}

resource "google_project_iam_member" "reservoirsim_sa_monitoring" {
  project = var.project_id
  role    = "roles/monitoring.metricWriter"
  member  = "serviceAccount:${google_service_account.reservoirsim_sa.email}"
}

# Cloud Storage bucket for data
resource "google_storage_bucket" "reservoirsim_data" {
  name     = "${var.project_id}-reservoirsim-data"
  location = var.region

  uniform_bucket_level_access = true

  versioning {
    enabled = true
  }

  lifecycle_rule {
    action {
      type = "Delete"
    }
    condition {
      age = 365
    }
  }
}

# Compute Instance
resource "google_compute_instance" "reservoirsim_instance" {
  name         = "reservoirsim-compute"
  machine_type = var.machine_type
  zone         = var.zone

  tags = ["reservoirsim"]

  boot_disk {
    initialize_params {
      image = "ubuntu-os-cloud/ubuntu-2204-lts"
      size  = var.boot_disk_size
      type  = "pd-ssd"
    }
  }

  # Optional: Add local SSD for temporary high-performance storage
  dynamic "scratch_disk" {
    for_each = var.local_ssd_count > 0 ? range(var.local_ssd_count) : []
    content {
      interface = "NVME"
    }
  }

  network_interface {
    subnetwork = google_compute_subnetwork.reservoirsim_subnet.id

    access_config {
      # Ephemeral public IP
    }
  }

  metadata = {
    enable-oslogin = "FALSE"
    ssh-keys       = var.ssh_key != "" ? "ubuntu:${var.ssh_key}" : ""
  }

  metadata_startup_script = templatefile("${path.module}/startup-script.sh", {
    bucket_name = google_storage_bucket.reservoirsim_data.name
  })

  service_account {
    email  = google_service_account.reservoirsim_sa.email
    scopes = ["cloud-platform"]
  }

  scheduling {
    # Use preemptible instances for cost savings (optional)
    preemptible       = var.use_preemptible
    automatic_restart = !var.use_preemptible
  }
}

# Static IP
resource "google_compute_address" "reservoirsim_ip" {
  name   = "reservoirsim-ip"
  region = var.region
}

resource "google_compute_instance" "reservoirsim_instance_static_ip" {
  count        = var.use_static_ip ? 1 : 0
  name         = "reservoirsim-compute-static"
  machine_type = var.machine_type
  zone         = var.zone

  tags = ["reservoirsim"]

  boot_disk {
    initialize_params {
      image = "ubuntu-os-cloud/ubuntu-2204-lts"
      size  = var.boot_disk_size
      type  = "pd-ssd"
    }
  }

  network_interface {
    subnetwork = google_compute_subnetwork.reservoirsim_subnet.id

    access_config {
      nat_ip = google_compute_address.reservoirsim_ip.address
    }
  }

  metadata_startup_script = templatefile("${path.module}/startup-script.sh", {
    bucket_name = google_storage_bucket.reservoirsim_data.name
  })

  service_account {
    email  = google_service_account.reservoirsim_sa.email
    scopes = ["cloud-platform"]
  }
}

# Outputs
output "instance_name" {
  value       = var.use_static_ip ? google_compute_instance.reservoirsim_instance_static_ip[0].name : google_compute_instance.reservoirsim_instance.name
  description = "Name of the compute instance"
}

output "instance_ip" {
  value       = var.use_static_ip ? google_compute_address.reservoirsim_ip.address : google_compute_instance.reservoirsim_instance.network_interface[0].access_config[0].nat_ip
  description = "Public IP address of the instance"
}

output "bucket_name" {
  value       = google_storage_bucket.reservoirsim_data.name
  description = "Name of the Cloud Storage bucket"
}

output "ssh_command" {
  value       = "gcloud compute ssh ubuntu@${var.use_static_ip ? google_compute_instance.reservoirsim_instance_static_ip[0].name : google_compute_instance.reservoirsim_instance.name} --zone=${var.zone}"
  description = "Command to SSH into the instance"
}

output "zone" {
  value       = var.zone
  description = "GCP zone where instance is deployed"
}
