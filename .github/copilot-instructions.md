# Copilot Instructions for FSRM

This file provides context and guidelines for GitHub Copilot when working in this repository.

## Repository Overview

FSRM (Full Service Reservoir Model) is a C++17/PETSc/MPI coupled multiphysics simulator for nuclear explosion monitoring, seismic wave propagation, dynamic fault rupture, and coupled THM poroelasticity. MIT licensed. Config-driven: a single executable reads `.config` files.

### Key Features (Integration-Tested)
- Elastostatics and elastodynamics with absorbing BCs
- Poroelasticity (Biot coupling, Terzaghi consolidation)
- Mueller-Murphy underground explosion source with moment tensor injection
- COUPLED_ANALYTIC NearField 1D solver to 3D FEM coupling
- Layered heterogeneous material via auxiliary fields
- Cohesive fault cells (locked, prescribed slip, time-dependent slip)
- Drucker-Prager elastoplasticity
- Gmsh mesh import with per-region materials
- SAC seismometer output
- HDF5/VTK output with checkpoint/restart
- 5 historic nuclear test digital twins (Gasbuggy, Gnome, Sedan, Degelen, Pahute Mesa)

## Technology Stack

- **Language**: C++17
- **Build System**: CMake >= 3.15
- **Parallel Computing**: PETSc 3.22.2, MPI
- **I/O**: HDF5
- **GPU**: PETSc CUDA backend (no FSRM code changes needed)
- **Testing**: Google Test (108 tests)
- **Build Environment**: Docker (Dockerfile.ci)

## Code Style Guidelines

FSRM follows the [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html) with modifications:

### Naming Conventions
- **Classes**: `PascalCase` (e.g., `FluidModel`, `FaultModel`)
- **Functions**: `camelCase` (e.g., `computeFlux`, `updateState`)
- **Variables**: `snake_case` (e.g., `permeability_x`, `total_stress`)
- **Constants**: `UPPER_SNAKE_CASE` (e.g., `MAX_ITERATIONS`, `PI`)
- **File names**: `PascalCase.cpp` / `PascalCase.hpp`
- **Namespace**: `FSRM`

### Formatting
- Indentation: 2 spaces (no tabs)
- Line length: 100 characters maximum
- Braces: Allman style (opening brace on new line)
- Comments: Use `//` for single line, `/* */` for multi-line

### Code Example
```cpp
namespace FSRM
{

class FluidModel
{
public:
  FluidModel(double density, double viscosity)
    : density_(density),
      viscosity_(viscosity)
  {
  }

  double computeVelocity(double pressure_gradient)
  {
    const double permeability = 1.0e-13;  // m^2
    return -permeability / viscosity_ * pressure_gradient;
  }

private:
  double density_;
  double viscosity_;
};

} // namespace FSRM
```

### Documentation
Use Doxygen-style comments for classes and functions:

```cpp
/**
 * @brief Compute Darcy flux from pressure gradient
 * 
 * @param pressure_gradient Pressure gradient in Pa/m
 * @param permeability Permeability in m^2
 * @param viscosity Dynamic viscosity in Pa·s
 * @return Darcy flux in m/s
 */
double computeDarcyFlux(double pressure_gradient, double permeability, double viscosity);
```

## Project Structure

```
FSRM/
├── include/          # Header files (.hpp)
├── src/              # Source files (.cpp, .cu for CUDA)
├── tests/            # Test files
│   ├── unit/         # Unit tests
│   ├── functional/   # Functional tests
│   ├── integration/  # Integration tests
│   ├── physics/      # Physics/MMS tests
│   └── performance/  # Performance tests
├── config/           # Configuration file templates
├── docs/             # Documentation
├── examples/         # Example configurations and executables
└── external/         # External dependencies
```

## Building and Testing

### Build Commands
```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON
make -j$(nproc)
```

### Running Tests
```bash
cd build
ctest                          # Run all tests
ctest -L "unit"                # Run unit tests only
ctest -R test_name -V          # Run specific test with verbose output
```

### Test Labels
- `unit`: Fast unit tests
- `functional`: Functional tests
- `physics`: Physics/MMS validation tests
- `integration`: Integration tests
- `performance`: Performance benchmarks

## Configuration Files

FSRM uses INI-style configuration files. When creating or modifying configurations:

- Use `[SECTION]` headers to organize parameters
- Include comments explaining parameters with units
- Specify units explicitly (e.g., `permeability = 150.0  # milliDarcy`)

Example:
```ini
[SIMULATION]
start_time = 0.0
end_time = 5.0            # seconds
fluid_model = NONE
solid_model = ELASTIC
enable_elastodynamics = true

[ROCK]
density = 2650.0          # kg/m^3
youngs_modulus = 60.0e9   # Pa
poissons_ratio = 0.25
```

## Git Commit Messages

Follow [Conventional Commits](https://www.conventionalcommits.org/):

```
<type>(<scope>): <subject>
```

Types: `feat`, `fix`, `docs`, `style`, `refactor`, `test`, `perf`, `build`, `ci`

## Important Notes

1. **PETSc 3.22.2**: Verify API signatures before calling any PETSc function. They change between versions.
2. **MPI**: Code must be MPI-safe for parallel execution.
3. **Error Handling**: Use `CHKERRQ(ierr)` for PETSc error propagation (not `PetscCall` -- this codebase uses the older pattern).
4. **DS/BC Ordering**: NEVER change the ordering in setupFields(): DMCreateDS -> setupBoundaryConditions -> DMSetLocalSection(nullptr) -> DMSetUp -> DMGetDS.
5. **Do NOT modify**: callback math in PetscFEElasticity.cpp, PetscFEPoroelasticity.cpp, PetscFEFluidFlow.cpp, FaultMeshManager::splitMeshAlongFault, or CohesiveFaultKernel::registerWithDS.
6. **Units**: Internal calculations use SI units; conversions at I/O boundaries.
7. **Build/test in Docker**: `docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:local bash -c 'make -j$(nproc) && ctest --output-on-failure'`
8. **No Python in simulator**: Python is ONLY for post-processing visualization scripts.
9. **Dead code**: Archived in `archive/`. Do not reference or use.
10. **Executable**: `./fsrm -c ../config/examples/example.config` (not `fsrm_simulator`).

## Common Patterns

### PETSc Error Handling
```cpp
PetscErrorCode MyFunction(PetscInt n)
{
  PetscFunctionBeginUser;
  
  // Use PetscCheck for input validation and assertions
  PetscCheck(n > 0, PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE, "n must be positive");
  
  Vec x;
  // Use PetscCall for wrapping PETSc function calls
  PetscCall(VecCreate(PETSC_COMM_WORLD, &x));
  PetscCall(VecSetSizes(x, PETSC_DECIDE, n));
  // ... operations ...
  PetscCall(VecDestroy(&x));
  PetscFunctionReturn(PETSC_SUCCESS);
}
```

- **PetscCall**: Wrap all PETSc function calls to propagate errors
- **PetscCheck**: Use for input validation and runtime assertions

### Configuration Reading
```cpp
ConfigReader config;
config.loadFile("simulation.config");
double permeability = config.getDouble("ROCK", "permeability_x", 100.0);
```

## Documentation References

- `docs/DEVELOPMENT.md` - Development guide
- `docs/CONFIGURATION.md` - Configuration reference
- `docs/USER_GUIDE.md` - User manual
- `docs/API_REFERENCE.md` - API documentation
