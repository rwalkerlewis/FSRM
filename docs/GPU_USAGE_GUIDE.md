# GPU Usage Guide for FSRM

This guide provides detailed instructions for using GPU acceleration in FSRM.

## Table of Contents
1. [System Requirements](#system-requirements)
2. [Installation](#installation)
3. [Configuration](#configuration)
4. [Performance Optimization](#performance-optimization)
5. [Troubleshooting](#troubleshooting)
6. [Best Practices](#best-practices)

## System Requirements

### NVIDIA GPUs (CUDA)
- **Minimum**: GTX 1060 (Pascal, SM 6.0), 6 GB VRAM
- **Recommended**: RTX 3080/3090, A100, 24+ GB VRAM
- **CUDA**: Version 11.0 or later (12.0+ recommended)
- **Driver**: Latest NVIDIA driver (525.x or later)

### AMD GPUs (ROCm/HIP)
- **Minimum**: MI50, 16 GB VRAM
- **Recommended**: MI100/MI200 series, 32+ GB VRAM
- **ROCm**: Version 5.0 or later
- **Driver**: Latest AMD ROCm driver

### Memory Requirements
For a simulation with `N` cells:
- **GPU Memory**: ~150 bytes/cell for basic physics
- **GPU Memory**: ~250 bytes/cell for poroelastodynamics with permeability tracking
- **CPU Memory**: ~100 bytes/cell (always required for host data)

Example memory requirements:
- 1M cells: ~250 MB GPU + 100 MB CPU
- 10M cells: ~2.5 GB GPU + 1 GB CPU
- 100M cells: ~25 GB GPU + 10 GB CPU

## Installation

### 1. Install CUDA Toolkit (for NVIDIA GPUs)

```bash
# Ubuntu/Debian
wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/cuda-keyring_1.1-1_all.deb
sudo dpkg -i cuda-keyring_1.1-1_all.deb
sudo apt-get update
sudo apt-get install cuda-toolkit-12-0

# Add to PATH
export PATH=/usr/local/cuda-12.0/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda-12.0/lib64:$LD_LIBRARY_PATH

# Verify installation
nvcc --version
nvidia-smi
```

### 2. Install ROCm (for AMD GPUs)

```bash
# Ubuntu/Debian
wget https://repo.radeon.com/rocm/rocm.gpg.key -O - | gpg --dearmor | sudo tee /etc/apt/keyrings/rocm.gpg
echo 'deb [arch=amd64 signed-by=/etc/apt/keyrings/rocm.gpg] https://repo.radeon.com/rocm/apt/debian jammy main' | sudo tee /etc/apt/sources.list.d/rocm.list
sudo apt-get update
sudo apt-get install rocm-hip-sdk

# Add to PATH
export PATH=/opt/rocm/bin:$PATH
export LD_LIBRARY_PATH=/opt/rocm/lib:$LD_LIBRARY_PATH

# Verify installation
rocminfo
```

### 3. Build FSRM with GPU Support

```bash
# Clone repository
git clone https://github.com/your-org/fsrm.git
cd fsrm

# Create build directory
mkdir build && cd build

# Configure with CUDA (NVIDIA)
cmake .. -DCMAKE_BUILD_TYPE=Release \
         -DENABLE_CUDA=ON \
         -DCMAKE_CUDA_ARCHITECTURES="80;86"  # Adjust for your GPU

# Or configure with HIP (AMD)
cmake .. -DCMAKE_BUILD_TYPE=Release \
         -DENABLE_HIP=ON

# Build
make -j$(nproc)

# Verify GPU support
./fsrm --version
./fsrm --gpu-info
```

## Configuration

### Basic GPU Configuration

Add to your `.config` file:

```ini
[SIMULATION]
use_gpu = true
gpu_mode = CPU_FALLBACK
gpu_device_id = 0
gpu_verbose = false
```

### GPU Execution Modes

1. **CPU_ONLY**: No GPU acceleration (default)
```ini
gpu_mode = CPU_ONLY
```

2. **GPU_ONLY**: Force GPU execution, fail if unavailable
```ini
gpu_mode = GPU_ONLY
```

3. **HYBRID**: Automatic work partitioning between CPU and GPU
```ini
gpu_mode = HYBRID
gpu_memory_fraction = 0.9
hybrid_cpu_fraction = 0.3
```

4. **CPU_FALLBACK**: Try GPU first, fallback to CPU (recommended)
```ini
gpu_mode = CPU_FALLBACK
```

### Advanced GPU Options

```ini
[SIMULATION]
# Memory management
gpu_memory_fraction = 0.85      # Use 85% of GPU memory
use_unified_memory = false      # Use unified memory (slower but larger problems)
pin_host_memory = true          # Pin CPU memory for faster transfers

# Solver options
use_gpu_preconditioner = true   # GPU-accelerated preconditioner
use_gpu_matrix_assembly = true  # Assemble matrices on GPU
gpu_block_size = 256            # CUDA block size (advanced)

# Hybrid execution
hybrid_partition_strategy = AUTOMATIC
hybrid_rebalance_frequency = 100

# Performance monitoring
enable_profiling = true
print_gpu_utilization = true
```

### Multi-GPU Configuration

For multiple GPUs, use MPI:

```bash
# Run on 4 GPUs
mpirun -np 4 fsrm -c config.ini --use-gpu

# Each MPI rank automatically uses a different GPU
# GPU assignment: rank 0 -> GPU 0, rank 1 -> GPU 1, etc.
```

Manual GPU assignment:
```bash
# Assign specific GPUs
mpirun -np 2 \
  -x CUDA_VISIBLE_DEVICES=0,2 \
  fsrm -c config.ini --use-gpu
```

## Performance Optimization

### 1. Choosing the Right Problem Size

GPU acceleration is most effective for:
- **Cell count**: >10,000 cells
- **Timesteps**: >100 timesteps
- **Physics**: Elastodynamics, poroelastodynamics
- **Grid type**: Structured grids

Less effective for:
- Small problems (<1000 cells)
- Single timestep
- I/O-bound simulations
- Highly irregular meshes

### 2. Memory Optimization

```ini
# For large problems that don't fit in GPU memory

# Option 1: Reduce memory usage
gpu_memory_fraction = 0.7
output_frequency = 1000         # Less frequent output

# Option 2: Use hybrid mode
gpu_mode = HYBRID
hybrid_cpu_fraction = 0.4

# Option 3: Use unified memory (slower)
use_unified_memory = true
gpu_memory_fraction = 1.5       # Can exceed GPU memory
```

### 3. Minimizing Data Transfers

```ini
# Keep data on GPU as much as possible
output_frequency = 100          # Don't output every timestep
checkpointing = false           # Disable frequent checkpointing
enable_compression = true       # Compress output data
```

### 4. Optimal GPU Block Sizes

For different problem types:

```ini
# Single-phase flow (1D problems)
gpu_block_size = 256

# Elastodynamics (3D problems)
gpu_block_size_3d = 8,8,8

# Poroelastodynamics (coupled)
gpu_block_size = 128
```

### 5. Benchmarking

Compare CPU vs GPU performance:

```bash
# Benchmark GPU performance
fsrm -c config.ini --benchmark-gpu

# Profile GPU kernels
fsrm -c config.ini --use-gpu --profile

# Generate performance report
fsrm -c config.ini --use-gpu --performance-report
```

## Troubleshooting

### GPU Not Detected

**Problem**: "No GPU devices found"

**Solutions**:
1. Check CUDA/ROCm installation:
```bash
nvidia-smi  # NVIDIA
rocm-smi    # AMD
```

2. Verify driver version:
```bash
nvidia-smi  # Should show CUDA version
```

3. Check build configuration:
```bash
./fsrm --version
# Should show "CUDA support: Enabled"
```

4. Set CUDA device explicitly:
```bash
export CUDA_VISIBLE_DEVICES=0
./fsrm -c config.ini --use-gpu
```

### Out of Memory Errors

**Problem**: "CUDA out of memory"

**Solutions**:

1. Reduce problem size:
```ini
gpu_memory_fraction = 0.6
```

2. Enable hybrid mode:
```ini
gpu_mode = HYBRID
```

3. Use unified memory:
```ini
use_unified_memory = true
```

4. Reduce output frequency:
```ini
output_frequency = 1000
```

5. Check memory usage:
```bash
nvidia-smi  # Monitor GPU memory
fsrm -c config.ini --use-gpu --print-memory-usage
```

### Poor Performance

**Problem**: GPU slower than CPU

**Possible causes**:

1. **Problem too small**:
   - Need >10k cells for GPU benefit
   - Try larger grid or longer simulation

2. **Data transfer overhead**:
   - Reduce output frequency
   - Keep data on GPU longer

3. **CPU bottleneck**:
   - Check if I/O is limiting performance
   - Use `--profile` to identify bottlenecks

4. **Wrong GPU mode**:
   - Try different gpu_mode settings
   - Benchmark each mode

### Multi-GPU Issues

**Problem**: Not utilizing all GPUs

**Solutions**:

1. Check MPI configuration:
```bash
# Should show all GPUs
mpirun -np 4 nvidia-smi
```

2. Verify GPU assignment:
```bash
fsrm -c config.ini --use-gpu --gpu-verbose
```

3. Enable GPU-aware MPI:
```bash
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_RDMA_ENABLED_CUDA=1
```

## Best Practices

### 1. Start Small, Scale Up
- Test with small problem first
- Verify correctness against CPU results
- Gradually increase problem size
- Benchmark at each scale

### 2. Choose Appropriate GPU Mode
- **Small problems (<10k cells)**: CPU_ONLY
- **Medium problems (10k-1M cells)**: GPU_ONLY or CPU_FALLBACK
- **Large problems (>1M cells)**: HYBRID
- **Very large (>100M cells)**: HYBRID with unified memory

### 3. Monitor Performance
```ini
enable_profiling = true
print_gpu_utilization = true
print_memory_usage = true
```

### 4. Optimize Data Layout
- Use structured grids when possible
- Minimize heterogeneity in material properties
- Regular mesh preferred over irregular

### 5. Balance Computation and I/O
```ini
# Don't output every timestep
output_frequency = 100

# Use compressed HDF5 for large data
output_format = HDF5
enable_compression = true

# Checkpoint less frequently on GPU
checkpoint_frequency = 1000
```

### 6. Test Correctness
```bash
# Compare CPU and GPU results
fsrm -c config.ini --use-gpu --compare-cpu-gpu
```

## Performance Guidelines

### Expected Speedups

| Problem Type | Problem Size | GPU Model | Speedup |
|--------------|--------------|-----------|---------|
| Single-phase flow | 1M cells | RTX 3080 | 10-15x |
| Elastodynamics | 1M cells | A100 | 30-50x |
| Poroelastodynamics | 1M cells | A100 | 20-40x |
| Black oil | 10M cells | A100 | 15-25x |

### Scaling Guidelines

| Grid Size | Cells | CPU Cores | GPU Model | Mode |
|-----------|-------|-----------|-----------|------|
| 50x50x50 | 125k | 4 | RTX 3070 | GPU_ONLY |
| 100x100x100 | 1M | 8 | RTX 3080 | GPU_ONLY |
| 200x200x100 | 4M | 16 | RTX 3090 | GPU_ONLY |
| 300x300x150 | 13M | 32 | A100 (40GB) | GPU_ONLY |
| 500x500x200 | 50M | 64 | A100 (80GB) | HYBRID |

## Support

For GPU-related issues:
- Check the troubleshooting section above
- Review example configurations in `config/gpu_*.config`
- Run diagnostic: `fsrm --gpu-info --gpu-verbose`
- Report issues on GitHub with output from `--gpu-info`
