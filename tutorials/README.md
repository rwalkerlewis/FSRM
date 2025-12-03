# FSRM Tutorials

Welcome to the comprehensive tutorial series for **FSRM (Full-Scale Reservoir Modeling)**!

These tutorials will guide you from basic installation to advanced multi-physics simulations.

---

## 📚 Tutorial Index

### Getting Started

| Tutorial | Description | Time |
|----------|-------------|------|
| [01. Getting Started](01_GETTING_STARTED.md) | Installation, first simulation, project structure | 30 min |
| [02. Configuration System](02_CONFIGURATION.md) | Config file format, sections, parameters | 45 min |

### Core Physics

| Tutorial | Description | Time |
|----------|-------------|------|
| [03. Physics Models](03_PHYSICS_MODELS.md) | Flow, mechanics, thermal models | 60 min |
| [04. Well Modeling](04_WELL_MODELING.md) | Producers, injectors, advanced wells | 45 min |
| [05. Fracture Modeling](05_FRACTURE_MODELING.md) | DFN, hydraulic fracturing, proppant | 60 min |
| [06. Fault Mechanics](06_FAULT_MECHANICS.md) | Friction laws, seismicity, IMEX | 60 min |

### Advanced Features

| Tutorial | Description | Time |
|----------|-------------|------|
| [07. GPU Acceleration](07_GPU_ACCELERATION.md) | CUDA setup, optimization, multi-GPU | 45 min |
| [08. Adaptive Mesh Refinement](08_ADAPTIVE_MESH_REFINEMENT.md) | AMR criteria, strategies, load balancing | 45 min |
| [09. Wave Propagation](09_WAVE_PROPAGATION.md) | Elastodynamics, sources, benchmarks | 60 min |
| [10. Advanced Simulations](10_ADVANCED_SIMULATIONS.md) | CO2 storage, geothermal, EOR, workflows | 90 min |

---

## 🚀 Quick Start

### Minimum Requirements

- CMake ≥ 3.15
- C++17 compiler (GCC 7+, Clang 6+)
- PETSc ≥ 3.15 with MPI

### Installation

```bash
# Install dependencies (Ubuntu/Debian)
sudo apt install -y cmake g++ libpetsc-dev libhdf5-mpi-dev

# Build FSRM
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

### First Simulation

```bash
# Run with default configuration
./fsrm -c ../config/default.config

# View results
paraview output/*.vtu
```

---

## 📋 Learning Paths

### Path 1: Reservoir Engineer

1. [Getting Started](01_GETTING_STARTED.md)
2. [Configuration System](02_CONFIGURATION.md)
3. [Physics Models](03_PHYSICS_MODELS.md) (focus on flow)
4. [Well Modeling](04_WELL_MODELING.md)
5. [Advanced Simulations](10_ADVANCED_SIMULATIONS.md) (EOR section)

### Path 2: Geomechanics Specialist

1. [Getting Started](01_GETTING_STARTED.md)
2. [Configuration System](02_CONFIGURATION.md)
3. [Physics Models](03_PHYSICS_MODELS.md) (focus on mechanics)
4. [Fracture Modeling](05_FRACTURE_MODELING.md)
5. [Fault Mechanics](06_FAULT_MECHANICS.md)

### Path 3: Seismologist

1. [Getting Started](01_GETTING_STARTED.md)
2. [Configuration System](02_CONFIGURATION.md)
3. [Fault Mechanics](06_FAULT_MECHANICS.md)
4. [Wave Propagation](09_WAVE_PROPAGATION.md)
5. [Advanced Simulations](10_ADVANCED_SIMULATIONS.md) (seismicity sections)

### Path 4: HPC Developer

1. [Getting Started](01_GETTING_STARTED.md)
2. [GPU Acceleration](07_GPU_ACCELERATION.md)
3. [Adaptive Mesh Refinement](08_ADAPTIVE_MESH_REFINEMENT.md)
4. [Advanced Simulations](10_ADVANCED_SIMULATIONS.md) (workflow automation)

---

## 🔧 Example Configurations

FSRM includes many example configurations in the `config/` directory:

| Category | Configuration Files |
|----------|-------------------|
| **Basic Flow** | `default.config`, `single_phase.config` |
| **Reservoir** | `spe1_benchmark.config`, `spe10_benchmark.config` |
| **Shale** | `shale_reservoir.config`, `hydraulic_fracturing.config` |
| **Geothermal** | `geothermal.config` |
| **CO2 Storage** | `co2_storage.config` |
| **Seismicity** | `induced_seismicity.config`, `induced_seismicity_imex.config` |
| **Earthquakes** | `scec_tpv3.config` through `scec_tpv37.config` |
| **Waves** | `elastodynamic_waves.config`, `poroelastodynamic_waves.config` |
| **Complete** | `complete_template.config` (all options) |

---

## 📖 Additional Documentation

- [User Guide](../docs/USER_GUIDE.md) - Complete reference documentation
- [Configuration Reference](../docs/CONFIGURATION.md) - All configuration options
- [Physics Models](../docs/PHYSICS_MODELS.md) - Mathematical background
- [API Reference](../docs/API_REFERENCE.md) - For developers

---

## 🆘 Getting Help

### Troubleshooting

Common issues and solutions are covered in each tutorial. See also:

- [Getting Started - Troubleshooting](01_GETTING_STARTED.md#troubleshooting)
- [GPU - Troubleshooting](07_GPU_ACCELERATION.md#troubleshooting)

### Debug Mode

```bash
# Build with debug symbols
cmake .. -DCMAKE_BUILD_TYPE=Debug

# Run with verbose output
./fsrm -c config.config -snes_monitor -ksp_monitor -log_view
```

### Community

- GitHub Issues: Report bugs and request features
- Discussions: Ask questions and share experiences

---

## 📊 Tutorial Prerequisites

| Tutorial | Prerequisites |
|----------|---------------|
| 01 | None |
| 02 | Tutorial 01 |
| 03 | Tutorials 01-02 |
| 04 | Tutorials 01-03 |
| 05 | Tutorials 01-03 |
| 06 | Tutorials 01-03, 05 |
| 07 | Tutorials 01-02, CUDA toolkit |
| 08 | Tutorials 01-03 |
| 09 | Tutorials 01-03, 06 |
| 10 | All previous tutorials |

---

## 🎯 Practice Exercises

Each tutorial includes hands-on exercises. Here's a summary:

### Beginner

1. Run `default.config` and visualize results
2. Modify grid resolution and observe changes
3. Add a production well to an existing config

### Intermediate

4. Create a five-spot well pattern simulation
5. Add hydraulic fractures to a horizontal well
6. Configure adaptive mesh refinement for a waterflood

### Advanced

7. Set up an induced seismicity simulation with IMEX
8. Create a multi-GPU simulation for a large grid
9. Implement a parameter sensitivity study
10. Reproduce a SCEC benchmark and validate results

---

## 📝 Changelog

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2024-01 | Initial tutorial series |

---

## 📜 License

FSRM and these tutorials are provided under the [MIT License](../LICENSE).

---

**Happy Simulating! 🚀**

Start your journey: [Tutorial 01: Getting Started →](01_GETTING_STARTED.md)
