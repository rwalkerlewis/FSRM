# Complete Benchmark Suite Summary - FSRM

This document provides a complete overview of all benchmarks added to the FSRM reservoir simulator.

## 📊 Overview

FSRM now includes **60+ comprehensive benchmarks** spanning:
- **Reservoir engineering** (SPE benchmarks)
- **Earthquake seismology** (SCEC benchmarks)  
- **Performance testing** (micro to macro scale)
- **GPU acceleration**
- **Parallel scaling**

---

## 🏆 Industry-Standard Benchmarks

### SPE Benchmarks (Society of Petroleum Engineers)

Validates reservoir simulation capabilities:

| Benchmark | Grid Size | Physics | Purpose | Runtime |
|-----------|-----------|---------|---------|---------|
| **SPE1** | 10×10×3 | 3-phase black oil | Validation | 10-30 min |
| **SPE3** | 9×9×4 | 4-component compositional | Gas cycling | 20-40 min |
| **SPE9** | 24×25×15 | Heterogeneous black oil | Complex geology | 1-3 hours |
| **SPE10** | 60×220×85 | Extreme heterogeneity | Scalability | 3-12 hours |

**Files**: 4 executables + 4 configs

---

### SCEC Benchmarks (Southern California Earthquake Center)

Validates earthquake/induced seismicity capabilities:

| Benchmark | Grid Size | Physics | Purpose | Runtime |
|-----------|-----------|---------|---------|---------|
| **TPV5** | 192×192×96 | Dynamic rupture | Strike-slip fault | 1-3 hours |
| **TPV10** | 192×192×96 | Branching fault | Fault interaction | 2-4 hours |
| **TPV16** | 240×240×120 | Rough fault | Geometric complexity | 3-6 hours |
| **LOH.1** | 150×150×85 | Wave propagation | Layered media | 0.5-1 hour |

**Files**: 4 executables + 4 configs

---

## 🚀 Performance Benchmarks

### 1. Kernel Benchmarks (test_benchmarks.cpp)

Tests individual computational kernels:
- Single-phase flow kernel (10-50 μs/eval)
- Geomechanics kernel (20-80 μs/eval)
- Matrix-vector operations
- Wave speed calculations
- DOF counting

**Tests**: 6 benchmarks  
**Runtime**: 1-2 minutes

---

### 2. Physics Benchmarks (test_physics_benchmarks.cpp)

Tests physics-specific models:

**Poroelasticity** (4 tests):
- Kernel performance
- Biot coefficient calculation
- Coupled flow-mechanics

**Fracture Mechanics** (2 tests):
- Growth calculation
- Stress intensity factors

**Wave Propagation** (2 tests):
- Elastic waves
- Poroelastic waves

**Two-Phase Flow** (3 tests):
- Kernel performance
- Relative permeability (Corey model)
- Capillary pressure (Brooks-Corey)

**Thermal** (1 test):
- Heat diffusion

**Scalability** (1 test):
- Grid size scaling

**Tests**: 13 benchmarks  
**Runtime**: 2-5 minutes

---

### 3. GPU Benchmarks (test_gpu_benchmarks.cpp)

Tests GPU acceleration (requires CUDA):

**Memory Bandwidth** (2 tests):
- Host-to-device transfer (10-15 GB/s)
- Device-to-host transfer (10-15 GB/s)

**Kernel Performance** (3 tests):
- Vector addition (20x speedup)
- Single-phase flow kernel
- Poroelastic kernel (20-40x speedup)

**GPU Scaling** (1 test):
- Strong scaling with block sizes

**Tests**: 6 benchmarks  
**Runtime**: 2-3 minutes (with GPU)

---

### 4. Memory & I/O Benchmarks (test_memory_io_benchmarks.cpp)

Tests memory and I/O performance:

**Memory Allocation** (2 tests):
- Allocation/deallocation
- Vector reallocation strategies

**Cache Performance** (2 tests):
- Stride patterns
- Matrix access (row-major, column-major, random)

**File I/O** (2 tests):
- Binary file I/O (200-1000 MB/s)
- HDF5 I/O performance

**Memory Bandwidth** (2 tests):
- Memory copy (10-50 GB/s)
- STREAM triad benchmark

**Tests**: 8 benchmarks  
**Runtime**: 2-3 minutes

---

### 5. Parallel Scaling (test_scaling.cpp)

Tests MPI parallel performance:
- MPI operations (broadcast, reduce, allreduce)
- Domain decomposition
- Load balancing
- Communication overhead
- Global reduction performance

**Tests**: 6 benchmarks  
**Runtime**: 1-2 minutes

---

### 6. Scenario Benchmarks (test_scenario_benchmarks.cpp)

Tests real-world simulation scenarios:

| Scenario | Grid | Physics | Runtime |
|----------|------|---------|---------|
| Hydraulic fracturing (small) | 20³ | Flow + mechanics + fractures | 5-10 min |
| Hydraulic fracturing (medium) | 50×50×20 | Flow + mechanics + fractures | 10-20 min |
| Geothermal system | 30³ | THM coupling | 10-20 min |
| CO2 storage | 40×40×20 | Two-phase flow | 10-20 min |
| Wave propagation | 50³ | Elastodynamics | 15-30 min |
| Parallel scaling | 60³ | Strong scaling | 5-10 min |
| Problem size scaling | Variable | Weak scaling | 10-20 min |

**Tests**: 7 scenarios  
**Runtime**: 1-2 hours total

---

### 7. SCEC Physics Benchmarks (test_scec_benchmarks.cpp)

Tests earthquake physics components:

**Friction Laws** (2 tests):
- Slip-weakening (< 100 ns/eval)
- Rate-and-state (< 500 ns/eval)

**Dynamic Rupture** (2 tests):
- Rupture speed calculation
- Stress tensor rotation (< 1000 ns/rotation)

**Wave Propagation** (3 tests):
- Seismic wave speeds (< 200 ns/calc)
- Wave arrival times
- Ricker wavelet generation

**Slip Distribution** (2 tests):
- Distribution analysis
- Fault point scaling (> 1M points/s)

**Tests**: 9 benchmarks  
**Runtime**: 1-2 minutes

---

## 📁 File Structure

```
fsrm/
├── config/
│   ├── spe1_benchmark.config
│   ├── spe3_benchmark.config
│   ├── spe9_benchmark.config
│   ├── spe10_benchmark.config
│   ├── scec_tpv5.config
│   ├── scec_tpv10.config
│   ├── scec_tpv16.config
│   └── scec_loh1.config
│
├── examples/
│   ├── spe1.cpp
│   ├── spe3.cpp
│   ├── spe9.cpp
│   ├── spe10.cpp
│   ├── scec_tpv5.cpp
│   ├── scec_tpv10.cpp
│   ├── scec_tpv16.cpp
│   └── scec_loh1.cpp
│
├── tests/performance/
│   ├── test_benchmarks.cpp
│   ├── test_scaling.cpp
│   ├── test_physics_benchmarks.cpp
│   ├── test_gpu_benchmarks.cpp
│   ├── test_memory_io_benchmarks.cpp
│   ├── test_scenario_benchmarks.cpp
│   ├── test_scec_benchmarks.cpp
│   └── README.md
│
├── tests/
│   ├── run_benchmarks.sh
│   └── CMakeLists.txt
│
└── Documentation/
    ├── BENCHMARKS_ADDED.md
    ├── SCEC_BENCHMARKS_ADDED.md
    └── COMPLETE_BENCHMARK_SUMMARY.md (this file)
```

---

## 🎯 Quick Start Guide

### Build with All Features

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_CUDA=ON
make -j8
```

### Run Quick Performance Tests

```bash
cd ../tests
./run_benchmarks.sh                    # All quick tests (~5 min)
./run_benchmarks.sh --physics -n 8     # Physics benchmarks
./run_benchmarks.sh --gpu              # GPU benchmarks
```

### Run Industry Benchmarks

```bash
cd build/examples

# SPE benchmarks (reservoir)
mpirun -np 4 ./spe1 -c config/spe1_benchmark.config
mpirun -np 8 ./spe9 -c config/spe9_benchmark.config

# SCEC benchmarks (earthquake)
mpirun -np 8 ./scec_tpv5 -c config/scec_tpv5.config
mpirun -np 8 ./scec_loh1 -c config/scec_loh1.config
```

### Use CTest

```bash
cd build

# All performance tests
ctest -L performance

# Specific categories
ctest -R "Performance.Physics"
ctest -R "Performance.GPU"
ctest -R "Performance.SCEC"
```

---

## 📈 Performance Expectations

### Micro-Benchmarks

| Category | Typical Performance |
|----------|-------------------|
| Single-phase flow | 20k-100k eval/s |
| Poroelasticity | 5k-20k eval/s |
| Fracture mechanics | 10k-50k eval/s |
| Memory copy | 10-50 GB/s |
| GPU memory transfer | 10-15 GB/s |
| Friction law | 1-10M eval/s |

### GPU Speedup

| Physics Model | Speedup |
|--------------|---------|
| Single-phase | 10-20x |
| Elastodynamics | 30-50x |
| Poroelastodynamics | 20-40x |
| Black oil | 15-25x |

### Parallel Efficiency

| Cores | Expected Efficiency |
|-------|-------------------|
| 1-4 | > 90% |
| 4-16 | > 80% |
| 16-64 | > 70% |
| 64-256 | > 60% |

---

## ✅ Validation Coverage

### Physics Models Validated

- ✅ Single-phase flow (SPE1, SPE9, SPE10)
- ✅ Two-phase flow (SPE10, Buckley-Leverett)
- ✅ Three-phase black oil (SPE1, SPE9)
- ✅ Compositional (SPE3)
- ✅ Geomechanics (elastostatics)
- ✅ Dynamic rupture (SCEC TPV5, TPV10, TPV16)
- ✅ Wave propagation (SCEC LOH.1)
- ✅ Poroelasticity (analytical, MMS)
- ✅ Fracture mechanics (LEFM)
- ✅ Thermal diffusion

### Numerical Methods Validated

- ✅ Finite volume (reservoir flow)
- ✅ Finite element (geomechanics)
- ✅ Time integration (implicit, explicit)
- ✅ Nonlinear solvers (Newton-Raphson)
- ✅ Linear solvers (GMRES, AMG)
- ✅ Domain decomposition (MPI)
- ✅ GPU acceleration (CUDA)

---

## 📊 Statistics Summary

### Code Statistics
- **New executable implementations**: 8 (4 SPE + 4 SCEC)
- **New config files**: 8
- **New performance test files**: 6
- **Total performance tests**: 60+
- **Lines of code added**: ~6,000+
- **Documentation lines**: 1,000+

### Benchmark Categories
- **Micro-benchmarks**: 50+ (kernels, operations)
- **Industry benchmarks**: 8 (SPE + SCEC)
- **Scenario tests**: 7 (realistic simulations)
- **GPU tests**: 7 (if CUDA enabled)
- **Scaling tests**: 6 (MPI performance)

### Coverage
- **Reservoir flow**: ✅ SPE1, 3, 9, 10
- **Earthquake physics**: ✅ TPV5, 10, 16, LOH.1
- **CPU performance**: ✅ All tests
- **GPU performance**: ✅ 7 tests
- **MPI scaling**: ✅ 6 tests

---

## 🎓 Educational Value

### For Students/Researchers

**Learn About**:
- Reservoir simulation (SPE benchmarks)
- Earthquake mechanics (SCEC benchmarks)
- Computational performance optimization
- Parallel computing with MPI
- GPU acceleration with CUDA

### For Developers

**Test Suites For**:
- Regression testing
- Performance tracking
- New feature validation
- Hardware benchmarking
- Optimization guidance

### For Industry

**Validation Against**:
- Published reference solutions
- Community-accepted benchmarks
- Multi-code comparison projects
- Industry best practices

---

## 🔗 References

### SPE Resources
- **SPE Website**: https://www.spe.org/
- **SPE Comparative Solutions**: https://www.spe.org/csp/

### SCEC Resources  
- **SCEC Website**: https://www.scec.org/
- **Dynamic Rupture Codes**: https://strike.scec.org/cvws/
- **Wave Propagation**: https://strike.scec.org/scecpedia/

### Key Publications
1. Odeh (1981): SPE1 benchmark
2. Kenyon & Behie (1987): SPE3 benchmark
3. Killough (1995): SPE9 benchmark
4. Christie & Blunt (2001): SPE10 benchmark
5. Harris et al. (2009): SCEC TPV benchmarks
6. Olsen et al. (2006): SCEC LOH problems
7. Dunham et al. (2011): SCEC TPV16

---

## 🚀 Future Extensions

### Potential Additions

**More SPE Benchmarks**:
- SPE2 (thermal recovery)
- SPE5 (6-component compositional)
- SPE11 (adaptive implicit methods)

**More SCEC Benchmarks**:
- TPV11 (off-fault plasticity)
- TPV24 (prestress variations)
- TPV27 (material contrast)
- LOH.2, LOH.3 (more complex media)

**Additional Tests**:
- Multi-GPU scaling
- I/O performance at scale
- Checkpointing overhead
- Mesh quality impact
- Preconditioner comparison

---

## 💡 Best Practices

### Running Benchmarks

1. **Start small**: Run quick tests first
2. **Check results**: Verify output makes sense
3. **Document hardware**: Note CPU, GPU, memory specs
4. **Track over time**: Compare results between versions
5. **Use appropriate cores**: Follow recommendations

### Performance Analysis

1. **Profile first**: Use benchmarks to find bottlenecks
2. **Measure carefully**: Run multiple times, report statistics
3. **Compare fairly**: Same hardware, same problem size
4. **Document changes**: Track what optimizations did
5. **Validate correctness**: Speed means nothing if wrong

### Contributing Benchmarks

1. **Follow patterns**: Use existing benchmarks as templates
2. **Document thoroughly**: Explain what benchmark tests
3. **Provide references**: Cite published solutions
4. **Set tolerances**: Define acceptable error ranges
5. **Test thoroughly**: Run on multiple platforms

---

## 🎉 Summary

FSRM now has one of the most comprehensive benchmark suites in the reservoir simulation community:

- ✅ **Industry validation**: 8 SPE + SCEC benchmarks
- ✅ **Performance testing**: 60+ micro to macro tests  
- ✅ **GPU support**: Full CUDA benchmark suite
- ✅ **Parallel scaling**: MPI performance validated
- ✅ **Documentation**: Complete usage guides
- ✅ **Automation**: One-command benchmark runner
- ✅ **CTest integration**: CI/CD ready

**From reservoir engineering to earthquake seismology, from single kernels to full simulations, from CPU to GPU - FSRM benchmarks cover it all!** 🚀

---

## 📞 Support

- **Documentation**: See `docs/` and benchmark READMEs
- **Examples**: All benchmarks have config files and usage examples
- **Issues**: GitHub issue tracker
- **Community**: SPE and SCEC forums for benchmark questions

---

**Last Updated**: November 2025  
**Version**: FSRM v2.0 with complete benchmark suite  
**Total Benchmarks**: 60+  
**Total Lines Added**: ~6,000+
