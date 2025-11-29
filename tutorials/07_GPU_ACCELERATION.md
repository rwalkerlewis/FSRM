# Tutorial 07: GPU Acceleration

Leverage GPU computing for large-scale simulations in FSRM.

## Table of Contents

1. [GPU Support Overview](#gpu-support-overview)
2. [Hardware Requirements](#hardware-requirements)
3. [Building with GPU Support](#building-with-gpu-support)
4. [Configuration](#configuration)
5. [GPU-Accelerated Physics](#gpu-accelerated-physics)
6. [Memory Management](#memory-management)
7. [Performance Optimization](#performance-optimization)
8. [Troubleshooting](#troubleshooting)

---

## GPU Support Overview

FSRM supports GPU acceleration through CUDA (NVIDIA) and HIP (AMD):

```
┌──────────────────────────────────────────────────────────────┐
│                  GPU Acceleration Stack                      │
├──────────────────────────────────────────────────────────────┤
│                                                               │
│  ┌─────────────┐    ┌─────────────┐    ┌─────────────┐       │
│  │   CUDA      │    │    HIP      │    │   OpenCL    │       │
│  │  (NVIDIA)   │    │   (AMD)     │    │  (Future)   │       │
│  └──────┬──────┘    └──────┬──────┘    └─────────────┘       │
│         │                  │                                  │
│         └────────┬─────────┘                                  │
│                  ▼                                            │
│         ┌────────────────┐                                   │
│         │  GPU Manager   │  ← Device selection, memory       │
│         └────────┬───────┘                                   │
│                  ▼                                            │
│         ┌────────────────┐                                   │
│         │ Physics Kernels│  ← GPU-optimized computations     │
│         └────────┬───────┘                                   │
│                  ▼                                            │
│  ┌───────────────────────────────────────────────────┐       │
│  │  cuSPARSE  │  cuBLAS  │  AmgX  │  Thrust  │ CUB  │       │
│  └───────────────────────────────────────────────────┘       │
│                                                               │
└──────────────────────────────────────────────────────────────┘
```

### Accelerated Components

| Component | GPU Speedup | Memory Requirement |
|-----------|-------------|-------------------|
| Linear algebra (SpMV) | 5-20× | High |
| Physics residuals | 10-50× | Medium |
| Jacobian assembly | 5-15× | High |
| Time integration | 10-30× | Medium |
| Post-processing | 3-10× | Low |

---

## Hardware Requirements

### Supported GPUs

**NVIDIA (CUDA):**
- Compute Capability ≥ 6.0 (Pascal or newer)
- Recommended: V100, A100, H100 for HPC
- Consumer: RTX 2000/3000/4000 series

**AMD (HIP):**
- GCN 3.0 or newer
- Recommended: MI100, MI200 series
- Consumer: RX 6000/7000 series

### Memory Requirements

| Problem Size | Minimum VRAM | Recommended |
|-------------|--------------|-------------|
| < 100K cells | 4 GB | 8 GB |
| 100K - 1M cells | 8 GB | 16 GB |
| 1M - 10M cells | 16 GB | 32 GB |
| > 10M cells | 32 GB+ | Multi-GPU |

### System Requirements

```bash
# Check NVIDIA GPU
nvidia-smi

# Check CUDA version
nvcc --version

# Check available memory
nvidia-smi --query-gpu=memory.free --format=csv
```

---

## Building with GPU Support

### CUDA Build

```bash
# Prerequisites
sudo apt install -y nvidia-cuda-toolkit

# Configure with CUDA
cd build
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DENABLE_CUDA=ON \
  -DCMAKE_CUDA_ARCHITECTURES="70;80;86" \
  -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda

# Build
make -j$(nproc)
```

### HIP Build (AMD)

```bash
# Prerequisites (ROCm)
sudo apt install -y rocm-dev hip-dev

# Configure with HIP
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DENABLE_HIP=ON \
  -DHIP_PATH=/opt/rocm/hip

make -j$(nproc)
```

### Verify GPU Support

```bash
# Check if GPU support is enabled
./fsrm -gpu_info

# Expected output:
# GPU Support: Enabled
# Backend: CUDA
# Device Count: 2
# Device 0: NVIDIA A100-SXM4-40GB
#   Compute Capability: 8.0
#   Memory: 40536 MB
#   Multiprocessors: 108
# Device 1: NVIDIA A100-SXM4-40GB
#   ...
```

---

## Configuration

### Basic GPU Configuration

```ini
[SIMULATION]
# Enable GPU
use_gpu = true

# GPU mode
gpu_mode = CPU_FALLBACK        # GPU_ONLY, CPU_FALLBACK, AUTO

# Device selection
gpu_device_id = 0              # Which GPU to use (0-indexed)

# Memory management
gpu_memory_fraction = 0.8      # Use up to 80% of VRAM
```

### GPU Mode Options

| Mode | Description | Use Case |
|------|-------------|----------|
| `GPU_ONLY` | All computation on GPU | Large problems, sufficient VRAM |
| `CPU_FALLBACK` | GPU primary, CPU backup | Typical use |
| `AUTO` | Automatic selection | Development/testing |

### Multi-GPU Configuration

```ini
[SIMULATION]
use_gpu = true
num_gpus = 2                   # Use multiple GPUs

# Distribution strategy
gpu_distribution = DOMAIN      # DOMAIN, PHYSICS, HYBRID

# Device assignment
gpu_devices = 0, 1             # Specific devices to use
```

### Memory Configuration

```ini
[GPU]
# Memory pool
enable_memory_pool = true
pool_initial_size = 1073741824  # 1 GB
pool_growth_factor = 2.0

# Unified memory (simpler but slower)
use_unified_memory = false

# Out-of-core for very large problems
enable_out_of_core = false
out_of_core_buffer_size = 4294967296  # 4 GB
```

---

## GPU-Accelerated Physics

### Supported Physics on GPU

| Physics | GPU Support | Notes |
|---------|-------------|-------|
| Single-phase flow | ✅ Full | Best speedup |
| Black oil | ✅ Full | |
| Compositional | ⚠️ Partial | Flash calc on CPU |
| Linear elasticity | ✅ Full | |
| Poroelasticity | ✅ Full | |
| Elastodynamics | ✅ Full | Excellent for waves |
| Plasticity | ⚠️ Partial | Return mapping on CPU |
| Fractures | ⚠️ Partial | DFN on GPU, propagation CPU |

### Enabling GPU for Specific Physics

```ini
[SIMULATION]
use_gpu = true

# Per-physics GPU settings
gpu_flow = true                # Flow equations on GPU
gpu_mechanics = true           # Mechanics on GPU
gpu_thermal = true             # Thermal on GPU
gpu_assembly = true            # Matrix assembly on GPU
gpu_solver = true              # Linear solver on GPU
```

### GPU Solver Options

```ini
[SIMULATION]
# GPU-accelerated linear solver
solver_type = amgx             # Use AmgX (NVIDIA)

# AmgX configuration
amgx_config_file = amgx_config.json

# Or use PETSc with GPU
# solver_type = gmres
# preconditioner = hypre_boomeramg
# use_petsc_gpu = true
```

### AmgX Configuration

Create `amgx_config.json`:

```json
{
  "config_version": 2,
  "solver": {
    "solver": "AMG",
    "presweeps": 1,
    "postsweeps": 1,
    "max_levels": 10,
    "coarsest_sweeps": 2,
    "cycle": "V",
    "selector": "PMIS",
    "coarse_solver": "DENSE_LU_SOLVER",
    "min_coarse_rows": 128,
    "interpolator": "D2",
    "smoother": {
      "solver": "MULTICOLOR_DILU",
      "coloring_level": 1
    }
  }
}
```

---

## Memory Management

### GPU Arrays

FSRM provides RAII wrappers for GPU memory:

```cpp
// In C++ code (for developers)
#include "GPUManager.hpp"

// Allocate GPU array
FSRM::GPUArray<double> pressure(num_cells);

// Copy from host
pressure.copyFromHost(host_pressure.data());

// Use in kernel
kernel<<<blocks, threads>>>(pressure.data(), ...);

// Copy back
pressure.copyToHost(host_pressure.data());
// Automatic cleanup on destruction
```

### Memory Monitoring

```ini
[OUTPUT]
# Monitor GPU memory usage
gpu_memory_stats = true
gpu_memory_log_frequency = 100  # Every 100 timesteps
```

Output format:
```
Step  100: GPU Memory: Used 2.4 GB / 8.0 GB (30%), Peak 3.1 GB
Step  200: GPU Memory: Used 2.4 GB / 8.0 GB (30%), Peak 3.1 GB
```

### Handling Out-of-Memory

```ini
[GPU]
# Automatic handling of OOM
on_oom = CPU_FALLBACK          # CPU_FALLBACK, ERROR, REDUCE

# Reduce problem size on GPU OOM
reduce_strategy = DOMAIN_DECOMP  # Split domain across CPU/GPU
```

---

## Performance Optimization

### Occupancy Optimization

```ini
[GPU]
# Thread block configuration
block_size_x = 256             # Threads per block
block_size_y = 1
block_size_z = 1

# Grid size (auto-calculated if not specified)
# grid_size_x = ...
```

### Memory Access Patterns

```ini
[GPU]
# Optimize memory access
enable_coalescing = true       # Coalesced memory access
enable_texture_cache = true    # Use texture cache for read-only data
enable_shared_memory = true    # Use shared memory for reductions
```

### Asynchronous Execution

```ini
[GPU]
# Overlap computation and communication
enable_async = true
num_streams = 4                # CUDA streams

# Overlap CPU-GPU transfers
enable_pinned_memory = true
```

### Kernel Fusion

```ini
[GPU]
# Fuse multiple small kernels
enable_kernel_fusion = true

# Reduces kernel launch overhead
# Example: residual computation + boundary conditions → single kernel
```

### Performance Monitoring

```bash
# Run with profiling
mpirun -np 4 ./fsrm -c config.config -gpu_profile

# Or use NVIDIA profiler
nsys profile mpirun -np 4 ./fsrm -c config.config
ncu --target-processes all mpirun -np 4 ./fsrm -c config.config
```

---

## Multi-GPU and Multi-Node

### Single Node, Multiple GPUs

```ini
[SIMULATION]
use_gpu = true
num_gpus = 4

# MPI rank to GPU mapping
gpu_mapping = ROUND_ROBIN      # ROUND_ROBIN, SEQUENTIAL, CUSTOM

# For CUSTOM mapping:
# gpu_map = 0:0, 1:1, 2:2, 3:3  # rank:device
```

### Multi-Node GPU Clusters

```bash
# Launch with GPU-aware MPI
mpirun -np 8 --map-by ppr:2:node ./fsrm -c config.config

# Environment for GPU-aware MPI
export OMPI_MCA_btl_smcuda_use_cuda_ipc=1
export OMPI_MCA_mpi_cuda_support=1
```

### NCCL for Multi-GPU Communication

```ini
[GPU]
# Use NCCL for GPU-GPU communication
enable_nccl = true
nccl_rings = 2                 # Number of rings

# P2P access
enable_p2p = true              # Direct GPU-GPU memory access
```

---

## Troubleshooting

### Common Issues

#### CUDA Out of Memory

```
Error: CUDA out of memory. Tried to allocate 2.00 GB
```

**Solutions:**
1. Reduce grid size or use domain decomposition
2. Increase `num_gpus` with MPI
3. Set `gpu_memory_fraction = 0.6`
4. Enable `use_unified_memory = true` (slower)

#### Wrong GPU Selected

```bash
# Check which GPU MPI ranks use
mpirun -np 4 ./fsrm -c config.config -gpu_debug

# Set CUDA visible devices
export CUDA_VISIBLE_DEVICES=0,1
```

#### Driver/Toolkit Mismatch

```
Error: CUDA driver version is insufficient for CUDA runtime version
```

**Solution:**
```bash
# Check versions
nvidia-smi  # Driver version
nvcc --version  # Toolkit version

# Update driver or rebuild with matching toolkit
```

#### Poor Performance

1. **Check GPU utilization:**
   ```bash
   nvidia-smi dmon -s u
   ```

2. **Profile for bottlenecks:**
   ```bash
   nsys profile ./fsrm -c config.config
   ```

3. **Common causes:**
   - Problem too small (GPU underutilized)
   - Too much CPU-GPU transfer
   - Suboptimal solver settings

### Debugging

```ini
[GPU]
# Enable debug output
debug = true
debug_level = 2                # 0=off, 1=errors, 2=warnings, 3=all

# Synchronous execution for debugging
enable_async = false

# Check all CUDA errors
check_errors = true
```

```bash
# Run with CUDA debugging
CUDA_LAUNCH_BLOCKING=1 ./fsrm -c config.config

# Memory debugging
cuda-memcheck ./fsrm -c config.config
```

---

## Example: Large-Scale GPU Simulation

### Configuration for 10M Cell Simulation

```ini
[SIMULATION]
name = large_scale_gpu
fluid_model = SINGLE_COMPONENT
solid_model = ELASTIC
enable_geomechanics = true

# GPU acceleration
use_gpu = true
gpu_mode = GPU_ONLY
gpu_device_id = 0
gpu_memory_fraction = 0.9

# Solver for GPU
solver_type = amgx
amgx_config_file = amgx_pcg.json

end_time = 86400.0

[GRID]
nx = 500
ny = 500
nz = 40                        # ~10M cells
Lx = 10000.0
Ly = 10000.0
Lz = 1000.0

[ROCK]
porosity = 0.2
permeability_x = 100.0
youngs_modulus = 20.0e9
poisson_ratio = 0.25

[FLUID]
density = 1000.0
viscosity = 0.001
compressibility = 4.5e-10

[OUTPUT]
format = HDF5                  # Efficient for large data
frequency = 100
compression = true
```

### Performance Comparison

| Configuration | 1M Cells | 10M Cells | 100M Cells |
|--------------|----------|-----------|------------|
| CPU (8 cores) | 45 min | 12 hr | - |
| 1× V100 | 5 min | 45 min | 8 hr |
| 4× V100 | 2 min | 15 min | 2 hr |
| 8× A100 | 1 min | 5 min | 45 min |

---

## GPU Best Practices

1. **Right-size the problem**: GPUs need enough work to be efficient (>100K cells)

2. **Minimize transfers**: Keep data on GPU, avoid frequent CPU-GPU copies

3. **Use appropriate precision**: `float` is 2× faster than `double` on many GPUs

4. **Profile before optimizing**: Use `nsys` and `ncu` to find bottlenecks

5. **Consider memory**: VRAM is limited; plan for larger problems

6. **Test CPU fallback**: Ensure code works without GPU for debugging

7. **Match MPI ranks to GPUs**: One rank per GPU typically optimal

---

**Previous**: [← Fault Mechanics](06_FAULT_MECHANICS.md) | **Next**: [Adaptive Mesh Refinement →](08_ADAPTIVE_MESH_REFINEMENT.md)
