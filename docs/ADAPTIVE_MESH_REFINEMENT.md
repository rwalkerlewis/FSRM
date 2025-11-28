# Adaptive Mesh Refinement (AMR) for Unstructured Grids

## Overview

The FSRM framework provides adaptive mesh refinement (AMR) capabilities for unstructured grids using PETSc's DMPlex and DMForest components. AMR enables efficient simulations by:

- **Locally refining** regions with high solution gradients (e.g., saturation fronts, near wells)
- **Coarsening** regions with smooth solutions to reduce computational cost
- **Dynamically adapting** the mesh as the simulation progresses
- **Preserving features** like wells, faults, and boundaries

## PETSc Foundation

The AMR implementation leverages several PETSc components:

| Component | Purpose |
|-----------|---------|
| **DMPlex** | Unstructured mesh management and topology |
| **DMForest** | Octree/quadtree-based hierarchical refinement |
| **DMPlexAdapt** | Metric-based mesh adaptation |
| **VecScatter** | Solution transfer between meshes |
| **ParMETIS/PTScotch** | Parallel mesh partitioning for load balancing |

## Error Estimators

### Gradient-Based (Kelly Estimator)

The Kelly error estimator measures the jump in normal gradient across cell faces:

$$\eta_K^2 = h_K \sum_{e \in \partial K} \|[\nabla u \cdot n]\|^2_{L^2(e)}$$

Best for: Elliptic problems, pressure fields

```cpp
#include "AdaptiveMeshRefinement.hpp"

FSRM::GradientErrorEstimator estimator;
estimator.setFieldIndex(0);      // Pressure field
estimator.setHPower(0.5);        // h-scaling exponent

std::vector<FSRM::CellErrorIndicator> errors;
estimator.estimate(dm, solution, errors);
```

### Hessian-Based

Uses recovered Hessian for smooth solutions:

$$\eta_K = h_K^2 \|H(u)\|_F$$

Best for: Smooth elliptic problems without shocks

```cpp
FSRM::HessianErrorEstimator estimator;
estimator.setRecoveryMethod(FSRM::HessianErrorEstimator::RecoveryMethod::SPR);
```

### Jump-Based

Tracks discontinuities in the solution (saturation fronts, phase boundaries):

```cpp
FSRM::JumpErrorEstimator estimator;
estimator.setFieldIndex(1);           // Saturation field
estimator.setJumpThreshold(0.1);      // Mark if jump > 0.1
```

### Feature-Based

Refines around geometric features (wells, faults, specific regions):

```cpp
FSRM::FeatureErrorEstimator estimator;

// Refine around wells
estimator.addWellLocation(100.0, 100.0, 0.0,  // x, y, z
                          20.0,                 // radius
                          3);                   // refinement level

// Refine along fault trace
std::vector<std::array<PetscReal, 3>> fault_trace = {
    {0.0, 0.0, 0.0},
    {500.0, 500.0, 0.0},
    {1000.0, 0.0, 0.0}
};
estimator.addFaultTrace(fault_trace, 10.0, 2);

// Refine in a box region
estimator.addRefinementBox(0.0, 100.0,    // x range
                           0.0, 100.0,    // y range
                           0.0, 50.0,     // z range
                           2);            // level
```

### Combined Estimator

Weighted combination of multiple criteria:

```cpp
auto grad_est = std::make_shared<FSRM::GradientErrorEstimator>();
auto jump_est = std::make_shared<FSRM::JumpErrorEstimator>();
auto feature_est = std::make_shared<FSRM::FeatureErrorEstimator>();

FSRM::CombinedErrorEstimator combined;
combined.addEstimator(grad_est, 0.3);     // 30% weight for pressure gradient
combined.addEstimator(jump_est, 1.0);     // 100% weight for saturation front
combined.addEstimator(feature_est, 0.5);  // 50% weight for features
```

## Refinement Strategies

### Fixed Fraction

Refine/coarsen a fixed fraction of cells based on error ranking:

```cpp
FSRM::RefinementParameters params;
params.strategy = FSRM::RefinementStrategy::FIXED_FRACTION;
params.refine_fraction = 0.20;   // Refine top 20%
params.coarsen_fraction = 0.05;  // Coarsen bottom 5%
```

### Threshold

Refine/coarsen based on error thresholds:

```cpp
params.strategy = FSRM::RefinementStrategy::THRESHOLD;
params.refine_threshold = 0.5;   // Refine if error/max_error > 0.5
params.coarsen_threshold = 0.05; // Coarsen if error/max_error < 0.05
```

### Equilibration

Target equal error distribution across all cells:

```cpp
params.strategy = FSRM::RefinementStrategy::EQUILIBRATION;
// Cells with error > 2x mean are refined
// Cells with error < 0.25x mean are coarsened
```

## Mesh Constraints

```cpp
FSRM::RefinementParameters params;

// Level limits
params.max_level = 5;           // Maximum refinement depth
params.min_level = 0;           // Minimum level (coarsest)

// Cell count limits
params.max_cells = 1000000;     // Prevent excessive refinement
params.min_cells = 100;         // Ensure minimum resolution

// Cell size limits
params.min_cell_size = 1.0;     // Minimum cell dimension [m]
params.max_cell_size = 100.0;   // Maximum cell dimension [m]

// Quality constraints
params.min_quality = 0.1;       // Minimum cell quality (0-1)
params.max_aspect_ratio = 10.0; // Maximum aspect ratio

// Feature preservation
params.preserve_boundaries = PETSC_TRUE;
params.preserve_wells = PETSC_TRUE;
params.preserve_faults = PETSC_TRUE;

// Buffer zones
params.buffer_layers = 2;       // Extra refinement layers around marked cells
```

## Adaptation Methods

### DMPlex Uniform Refinement

Simple bisection refinement:

```cpp
amr.setMethod(FSRM::AdaptationMethod::PLEX_REFINE);
```

### DMPlexAdapt (Metric-Based)

Uses a metric tensor field to guide anisotropic adaptation:

```cpp
amr.setMethod(FSRM::AdaptationMethod::PLEX_ADAPT);
// Requires metric tensor at vertices
```

### DMForest (p4est)

Octree/quadtree-based hierarchical refinement:

```cpp
amr.setMethod(FSRM::AdaptationMethod::FOREST_P4EST);
// Requires DMForest-compatible mesh
```

### External Adapters

Integration with external mesh adaptation libraries:

```cpp
// PRAgMaTIc (anisotropic)
amr.setMethod(FSRM::AdaptationMethod::PRAGMATIC);

// MMG (remeshing)
amr.setMethod(FSRM::AdaptationMethod::MMG);
```

## Complete Usage Example

```cpp
#include "AdaptiveMeshRefinement.hpp"

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    
    // Load or create mesh
    DM dm;
    DMPlexCreateBoxMesh(PETSC_COMM_WORLD, 2, PETSC_TRUE,
                        (PetscInt[]){10, 10},
                        (PetscReal[]){0.0, 0.0},
                        (PetscReal[]){1000.0, 1000.0},
                        nullptr, PETSC_TRUE, &dm);
    
    // Initialize AMR
    FSRM::AdaptiveMeshRefinement amr;
    amr.initialize(dm);
    
    // Configure parameters
    FSRM::RefinementParameters params;
    params.strategy = FSRM::RefinementStrategy::FIXED_FRACTION;
    params.refine_fraction = 0.25;
    params.coarsen_fraction = 0.10;
    params.max_level = 4;
    params.max_cells = 50000;
    params.buffer_layers = 1;
    amr.setParameters(params);
    
    // Set up combined error estimator
    auto grad_est = std::make_shared<FSRM::GradientErrorEstimator>();
    grad_est->setFieldIndex(0);  // Pressure
    
    auto jump_est = std::make_shared<FSRM::JumpErrorEstimator>();
    jump_est->setFieldIndex(1);  // Saturation
    jump_est->setJumpThreshold(0.1);
    
    auto combined = std::make_shared<FSRM::CombinedErrorEstimator>();
    combined->addEstimator(grad_est, 0.3);
    combined->addEstimator(jump_est, 1.0);
    amr.setErrorEstimator(combined);
    
    // Add well refinement
    amr.addWellRefinement(50.0, 50.0, 0.0, 25.0, 3);    // Producer
    amr.addWellRefinement(950.0, 950.0, 0.0, 30.0, 3);  // Injector
    
    // Set callbacks (optional)
    amr.setPreAdaptCallback([](DM dm, Vec sol) {
        PetscPrintf(PETSC_COMM_WORLD, "Pre-adaptation checkpoint\n");
        return PETSC_SUCCESS;
    });
    
    amr.setPostAdaptCallback([](DM dm_old, DM dm_new, Vec sol_old, Vec sol_new) {
        PetscPrintf(PETSC_COMM_WORLD, "Post-adaptation: solution transferred\n");
        return PETSC_SUCCESS;
    });
    
    // Time stepping loop
    Vec solution;
    DMGetGlobalVector(dm, &solution);
    // ... initialize solution ...
    
    for (int step = 0; step < 100; step++) {
        // ... advance physics ...
        
        // Check if adaptation needed
        PetscBool need_adapt;
        amr.checkAdaptation(solution, &need_adapt);
        
        if (need_adapt || step % 10 == 0) {
            DM dm_new;
            Vec solution_new;
            
            // Perform adaptation
            amr.adapt(solution, &dm_new, &solution_new);
            
            // Print statistics
            amr.printStatistics();
            
            // Update references
            DMRestoreGlobalVector(dm, &solution);
            DMDestroy(&dm);
            dm = dm_new;
            solution = solution_new;
        }
        
        // Output
        if (step % 10 == 0) {
            char filename[256];
            snprintf(filename, sizeof(filename), "mesh_%04d.vtk", step);
            amr.writeAdaptedMesh(filename);
        }
    }
    
    DMRestoreGlobalVector(dm, &solution);
    DMDestroy(&dm);
    PetscFinalize();
    
    return 0;
}
```

## Configuration File

```ini
[amr]
enabled = true
method = "plex_adapt"

[amr.criterion]
type = "combined"
weight_pressure = 0.3
weight_saturation = 1.0

[amr.strategy]
type = "fixed_fraction"
refine_fraction = 0.20
coarsen_fraction = 0.05

[amr.limits]
max_refinement_level = 4
max_cells = 100000
min_cell_size = 1.0

[amr.quality]
min_quality = 0.2
max_aspect_ratio = 5.0

[amr.features]
buffer_layers = 2

[[amr.features.wells]]
name = "PROD1"
x = 50.0
y = 50.0
z = 0.0
radius = 20.0
level = 3
```

## Solution Transfer

The `SolutionTransfer` class handles moving solution fields between meshes:

```cpp
FSRM::SolutionTransfer transfer;
transfer.setMethod(FSRM::SolutionTransfer::TransferMethod::INTERPOLATION);
transfer.setConservative(PETSC_TRUE);  // Mass-conserving transfer

Vec solution_new;
transfer.transfer(dm_old, dm_new, solution_old, solution_new);

// Transfer multiple fields
std::vector<Vec> fields_old = {pressure, saturation, temperature};
std::vector<Vec> fields_new;
transfer.transferFields(dm_old, dm_new, fields_old, fields_new);
```

### Transfer Methods

| Method | Description | Conservation |
|--------|-------------|--------------|
| `INTERPOLATION` | Finite element interpolation | Approximate |
| `PROJECTION` | L2 projection | Exact (weak) |
| `INJECTION` | Direct value copy (coarsening) | Approximate |
| `CONSERVATIVE` | Mass-weighted averaging | Exact |

## Mesh Quality

The `MeshQuality` class evaluates mesh health:

```cpp
FSRM::MeshQuality quality;
quality.setQualityMeasure(FSRM::MeshQuality::QualityMeasure::SCALED_JACOBIAN);
quality.computeQuality(dm);

PetscReal min_q = quality.getMinQuality();      // 0 to 1
PetscReal mean_q = quality.getMeanQuality();
PetscReal max_ar = quality.getMaxAspectRatio();
PetscInt bad_cells = quality.getNumBadCells();

// Check if acceptable
if (!quality.isAcceptable(params)) {
    // Mesh needs improvement
}
```

### Quality Measures

- **Scaled Jacobian**: Ratio of minimum to maximum Jacobian (ideal = 1)
- **Shape**: Deviation from ideal element shape
- **Aspect Ratio**: Ratio of longest to shortest edge
- **Skewness**: Angular deviation from ideal
- **Condition Number**: Jacobian matrix condition

## Parallel Considerations

### Load Balancing

After adaptation, the mesh may become unbalanced. Use `rebalance()`:

```cpp
// After adaptation
amr.adapt(solution, &dm_new, &solution_new);

// Rebalance across processors
amr.rebalance(dm_new);
```

### PETSc Options

```bash
# Command-line options for AMR
-amr_refine_fraction 0.2
-amr_coarsen_fraction 0.05
-amr_max_level 4
-amr_max_cells 100000
-amr_min_cell_size 1.0

# Partitioning options
-dm_plex_partition_type parmetis
-dm_plex_partition_balance

# Debugging
-dm_plex_check_all
-dm_view
```

## Statistics and Output

```cpp
const FSRM::AMRStatistics& stats = amr.getStatistics();

PetscPrintf(PETSC_COMM_WORLD,
    "AMR Step: %d cells -> %d cells\n"
    "  Refined: %d, Coarsened: %d\n"
    "  Error: min=%.2e, max=%.2e, mean=%.2e\n"
    "  Time: adapt=%.3fs, transfer=%.3fs, balance=%.3fs\n",
    stats.num_cells_before, stats.num_cells_after,
    stats.num_refined, stats.num_coarsened,
    stats.min_error, stats.max_error, stats.mean_error,
    stats.adaptation_time, stats.transfer_time, stats.balance_time);

// Write mesh and error field
amr.writeAdaptedMesh("adapted_mesh.vtk");
amr.writeErrorField("error_field.vtk");
```

## Best Practices

1. **Start Coarse**: Begin with a coarse mesh and let AMR add resolution where needed
2. **Limit Levels**: Keep `max_level` reasonable (4-6) to avoid too small cells
3. **Use Buffers**: Set `buffer_layers ≥ 1` to prevent refinement boundaries from affecting solution
4. **Conservative Transfer**: Enable for mass-sensitive problems (reservoir simulation)
5. **Rebalance Regularly**: Especially after major refinement changes
6. **Monitor Quality**: Check mesh quality after adaptation; poor quality can cause solver issues
7. **Adapt Strategically**: Don't adapt every time step; use `checkAdaptation()` or fixed intervals

## Performance Tips

- Use `FOREST_P4EST` for geometrically simple domains (faster adaptation)
- Use `PLEX_ADAPT` with external libraries (PRAgMaTIc, MMG) for complex geometries
- Set appropriate `max_cells` to control memory usage
- Use parallel partitioners (ParMETIS, PTScotch) for large-scale problems
- Consider amortizing adaptation cost over multiple time steps

## References

1. PETSc DMPlex Documentation: https://petsc.org/release/docs/manual/dmplex/
2. p4est Library: http://www.p4est.org/
3. Kelly, D.W., et al. "A posteriori error analysis and adaptive processes in the finite element method" (1983)
4. Zienkiewicz, O.C., Zhu, J.Z. "The superconvergent patch recovery and a posteriori error estimates" (1992)
