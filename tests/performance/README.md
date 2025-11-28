# FSRM Performance Benchmarks

This directory contains comprehensive performance benchmarks for the FSRM reservoir simulator.

## Overview

The benchmark suite tests performance across multiple dimensions:

### 1. **Kernel Benchmarks** (`test_benchmarks.cpp`)
- Single-phase flow kernel performance
- Geomechanics kernel performance
- Mathematical operations (matrix-vector multiply, Lame parameters)
- Wave speed calculations
- Memory access patterns
- DOF counting

### 2. **Physics Benchmarks** (`test_physics_benchmarks.cpp`)
- **Poroelasticity**: Biot coefficient, coupled flow-mechanics
- **Fracture Mechanics**: Stress intensity factors, crack propagation
- **Wave Propagation**: Elastic waves, poroelastic waves
- **Two-Phase Flow**: Relative permeability, capillary pressure
- **Thermal**: Heat diffusion
- **Grid Size Scalability**: Performance vs problem size

### 3. **GPU Benchmarks** (`test_gpu_benchmarks.cpp`)
- Memory bandwidth (host-device, device-host)
- Kernel performance (vector operations, physics kernels)
- CPU vs GPU speedup comparisons
- GPU strong scaling
- Single-phase and poroelastic kernels on GPU

### 4. **Memory & I/O Benchmarks** (`test_memory_io_benchmarks.cpp`)
- Memory allocation/deallocation
- Vector reallocation strategies
- Cache performance (stride patterns, matrix access)
- Binary file I/O
- HDF5 I/O performance
- Memory bandwidth (copy, STREAM triad)

### 5. **Parallel Scaling** (`test_scaling.cpp`)
- MPI operations (broadcast, reduce, allreduce)
- Domain decomposition
- Load balancing
- Communication overhead
- Global reduction performance

### 6. **Scenario Benchmarks** (`test_scenario_benchmarks.cpp`)
Real-world simulation scenarios:
- **Hydraulic Fracturing**: Small (4k cells) and Medium (50k cells)
- **Geothermal Systems**: THM coupling (27k cells)
- **CO2 Storage**: Two-phase flow (32k cells)
- **Wave Propagation**: Elastodynamics (125k cells)
- **Parallel Scaling**: Strong scaling tests
- **Problem Size Scaling**: Weak scaling tests

## SPE Benchmarks

Industry-standard SPE (Society of Petroleum Engineers) benchmarks:

### SPE1 (`../examples/spe1.cpp`)
- **Problem**: 3-phase black oil with gas dissolution
- **Grid**: 10×10×3 (300 cells)
- **Purpose**: Validation of black oil formulation
- **Reference**: Odeh (1981)

### SPE3 (`../examples/spe3.cpp`)
- **Problem**: Gas cycling of retrograde condensate
- **Grid**: 9×9×4 (324 cells)
- **Purpose**: Compositional modeling with 4 components
- **Reference**: Kenyon & Behie (1987)

### SPE9 (`../examples/spe9.cpp`)
- **Problem**: North Sea reservoir model
- **Grid**: 24×25×15 (9,000 cells)
- **Purpose**: Complex heterogeneous black oil simulation
- **Reference**: Killough (1995)

### SPE10 (`../examples/spe10.cpp`)
- **Problem**: Large-scale heterogeneous reservoir
- **Grid**: 60×220×85 (1.1 million cells)
- **Purpose**: Scalability with extreme permeability variations
- **Reference**: Christie & Blunt (2001)

## Running Benchmarks

### Quick Start

```bash
# Run all benchmarks with default settings (4 processes)
cd tests
./run_benchmarks.sh

# Run specific benchmark categories
./run_benchmarks.sh --kernel      # Kernel benchmarks only
./run_benchmarks.sh --physics     # Physics benchmarks only
./run_benchmarks.sh --gpu         # GPU benchmarks only
./run_benchmarks.sh --memory      # Memory/IO benchmarks only
./run_benchmarks.sh --scaling     # Parallel scaling tests
./run_benchmarks.sh --scenario    # Scenario benchmarks (long)
./run_benchmarks.sh --spe         # SPE benchmarks (very long)

# Specify number of MPI processes
./run_benchmarks.sh --physics -n 8

# Verbose output
./run_benchmarks.sh --kernel -v

# Custom output directory
./run_benchmarks.sh -o my_results
```

### Using CTest

```bash
cd build

# Run all performance tests
ctest -L performance

# Run specific test categories
ctest -R "Performance.Benchmarks"
ctest -R "Performance.PhysicsBenchmarks"
ctest -R "Performance.GPUBenchmarks"
ctest -R "Performance.MemoryIO"
ctest -R "Performance.Scenarios"

# Run with verbose output
ctest -L performance -V

# Run only GPU benchmarks
ctest -L gpu
```

### Direct Execution

```bash
cd build/tests

# Run all performance tests
mpirun -np 4 ./run_performance_tests

# Run specific benchmark groups
mpirun -np 4 ./run_performance_tests --gtest_filter=BenchmarkTest.*
mpirun -np 4 ./run_performance_tests --gtest_filter=PhysicsBenchmark.*
mpirun -np 4 ./run_performance_tests --gtest_filter=GPUBenchmark.*
mpirun -np 4 ./run_performance_tests --gtest_filter=MemoryIOBenchmark.*
mpirun -np 4 ./run_performance_tests --gtest_filter=ScalingTest.*
mpirun -np 4 ./run_performance_tests --gtest_filter=ScenarioBenchmark.*
```

### Running SPE Benchmarks

```bash
cd build/examples

# SPE1 - Small benchmark (quick)
mpirun -np 4 ./spe1 -c config/spe1_benchmark.config

# SPE3 - Compositional model
mpirun -np 4 ./spe3 -c config/spe3_benchmark.config

# SPE9 - Medium-scale benchmark
mpirun -np 8 ./spe9 -c config/spe9_benchmark.config

# SPE10 - Large-scale benchmark (requires >= 16 processes)
mpirun -np 32 ./spe10 -c config/spe10_benchmark.config
```

## Expected Performance

### Kernel Performance (typical workstation)
| Kernel | Time per eval | Throughput |
|--------|--------------|------------|
| Single-phase flow | ~10-50 μs | 20k-100k eval/s |
| Geomechanics | ~20-80 μs | 12k-50k eval/s |
| Poroelasticity | ~50-200 μs | 5k-20k eval/s |

### GPU Speedup (NVIDIA A100)
| Physics Model | CPU Time | GPU Time | Speedup |
|--------------|----------|----------|---------|
| Single-phase | 100 ms | 5 ms | 20x |
| Elastodynamics | 200 ms | 5 ms | 40x |
| Poroelasticity | 300 ms | 15 ms | 20x |

### Memory Bandwidth
| Operation | Expected Performance |
|-----------|---------------------|
| Memory copy | 10-50 GB/s |
| STREAM triad | 10-40 GB/s |
| GPU memory (H2D) | 10-15 GB/s |
| GPU memory (D2H) | 10-15 GB/s |

### Parallel Scaling Efficiency
| Processes | Expected Efficiency |
|-----------|-------------------|
| 1-4 | > 90% |
| 4-16 | > 80% |
| 16-64 | > 70% |
| 64-256 | > 60% |

## Interpreting Results

### Kernel Benchmarks
- **Time per eval**: Lower is better
- **Typical values**: 10-200 μs depending on complexity
- **Watch for**: Significant increases may indicate optimization opportunities

### GPU Benchmarks
- **Speedup**: GPU time / CPU time
- **Good speedup**: > 10x for large problems
- **Poor speedup**: < 2x may indicate:
  - Problem too small
  - Memory bandwidth limited
  - Excessive data transfer

### Memory Benchmarks
- **Cache performance**: Row-major should be faster than column-major
- **Memory bandwidth**: Compare to theoretical peak (check specs)
- **I/O performance**: 
  - HDF5: 100-500 MB/s typical
  - Binary: 200-1000 MB/s typical

### Scaling Benchmarks
- **Parallel efficiency** = (Speedup / # processes)
- **Good efficiency**: > 80%
- **Poor efficiency**: < 50% indicates:
  - Communication overhead
  - Load imbalance
  - Serial bottlenecks

## Benchmark Configuration

Benchmarks can be configured through:

1. **Command-line options** (via run_benchmarks.sh)
2. **Environment variables**:
   ```bash
   export OMP_NUM_THREADS=8       # OpenMP threads
   export CUDA_VISIBLE_DEVICES=0  # GPU device selection
   ```
3. **Config files** (for scenario benchmarks)

## Output Files

Benchmark results are saved to `benchmark_results/` directory:

```
benchmark_results/
├── kernel_benchmarks_YYYYMMDD_HHMMSS.log
├── physics_benchmarks_YYYYMMDD_HHMMSS.log
├── gpu_benchmarks_YYYYMMDD_HHMMSS.log
├── memory_io_benchmarks_YYYYMMDD_HHMMSS.log
├── scaling_tests_YYYYMMDD_HHMMSS.log
├── scenario_benchmarks_YYYYMMDD_HHMMSS.log
├── spe1_YYYYMMDD_HHMMSS.log
├── spe3_YYYYMMDD_HHMMSS.log
├── spe9_YYYYMMDD_HHMMSS.log
├── spe10_YYYYMMDD_HHMMSS.log
└── summary_YYYYMMDD_HHMMSS.txt
```

## Troubleshooting

### GPU benchmarks fail
- **Check**: `nvidia-smi` shows GPU
- **Solution**: Rebuild with `-DENABLE_CUDA=ON`

### MPI errors
- **Check**: Proper MPI installation
- **Solution**: Verify with `mpirun --version`

### Out of memory
- **For SPE10**: Use >= 32 GB RAM
- **Solution**: Reduce problem size or use more processes

### Poor performance
1. Check CPU frequency scaling
2. Disable turbo boost for consistent results
3. Run on dedicated node (no other jobs)
4. Check for thermal throttling

## Performance Regression Testing

To track performance over time:

```bash
# Run benchmarks and save results
./run_benchmarks.sh -o results_v1.0
./run_benchmarks.sh -o results_v1.1

# Compare results
diff results_v1.0/summary_*.txt results_v1.1/summary_*.txt
```

## Contributing

When adding new benchmarks:

1. Add test to appropriate file (or create new file)
2. Update CMakeLists.txt
3. Update this README
4. Verify with `./run_benchmarks.sh`
5. Document expected performance

## References

- **SPE1**: Odeh, A.S., JPT (1981)
- **SPE3**: Kenyon & Behie, JPT (1987)
- **SPE9**: Killough, SPE 29110 (1995)
- **SPE10**: Christie & Blunt, SPE 66599 (2001)
- **STREAM**: McCalpin, https://www.cs.virginia.edu/stream/

## Support

For issues or questions about benchmarks:
- Check existing documentation
- Review log files in benchmark_results/
- Open GitHub issue with:
  - System specifications
  - Benchmark command used
  - Complete log file
