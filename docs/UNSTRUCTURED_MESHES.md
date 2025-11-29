# Unstructured Mesh Support with DMPlex

FSRM uses PETSc's DMPlex for all spatial domain representation. This provides consistent unstructured mesh handling across all simulation types, enabling complex geometries, adaptive mesh refinement, and advanced mesh operations.

## Overview

All problems in FSRM use DMPlex-based unstructured grids for the spatial domain:
- **Cartesian grids** are created as DMPlex box meshes (hexahedral or tetrahedral cells)
- **External meshes** (Gmsh, Exodus) are loaded directly into DMPlex
- **Adaptive mesh refinement** is supported through DMPlex's native capabilities

Gmsh is an open-source finite element mesh generator that produces high-quality tetrahedral, hexahedral, and mixed-element meshes. FSRM can directly load Gmsh `.msh` files and convert them to PETSc DMPlex for simulation.

## Supported Features

- **Mesh Formats**: Gmsh MSH v2.2 and v4.x (ASCII)
- **Element Types**: Triangles, quads, tetrahedra, hexahedra, prisms, pyramids
- **Higher Order**: Second-order elements supported
- **Physical Groups**: Named regions for material properties and boundary conditions
- **Parallel Distribution**: Automatic mesh partitioning via PETSc

## Configuration

### Basic Usage

```ini
[GRID]
mesh_type = GMSH
mesh_file = meshes/my_reservoir.msh
use_unstructured = true
```

### With Physical Groups

Physical groups defined in Gmsh can be referenced for boundary conditions:

```ini
[GRID]
mesh_type = GMSH
mesh_file = meshes/reservoir.msh

# Physical group names must match those in the Gmsh file
gmsh_physical_volume = reservoir_volume
gmsh_boundaries = inlet, outlet, walls, top, bottom

[BC1]
type = DIRICHLET
field = PRESSURE
location = inlet           # Reference physical group name
value = 15.0e6

[BC2]
type = NEUMANN
field = PRESSURE
location = walls
value = 0.0               # No-flow
```

### With Coordinate Transformation

Gmsh meshes can be combined with coordinate system transformations:

```ini
[GRID]
mesh_type = GMSH
mesh_file = meshes/reservoir.msh

# Input coordinates in WGS84
input_crs = EPSG:4326

# Model in UTM Zone 10N
model_crs = EPSG:32610

# Center on project site
use_local_coordinates = true
local_origin_x = -122.4194    # Longitude
local_origin_y = 37.7749      # Latitude
```

## Creating Gmsh Meshes

### Using Gmsh GUI

1. Open Gmsh
2. Create geometry (CAD import or direct modeling)
3. Define Physical Groups for regions and boundaries
4. Generate mesh
5. Save as `.msh` file

### Using Gmsh Scripting

Example `.geo` file for a simple reservoir:

```geo
// Simple reservoir geometry
SetFactory("OpenCASCADE");

// Create box (reservoir domain)
Box(1) = {0, 0, -100, 1000, 500, 100};

// Define physical groups
Physical Volume("reservoir") = {1};
Physical Surface("inlet") = {1};      // West face
Physical Surface("outlet") = {2};     // East face
Physical Surface("walls") = {3, 4};   // North/South
Physical Surface("top") = {5};
Physical Surface("bottom") = {6};

// Mesh settings
Mesh.CharacteristicLengthMin = 20;
Mesh.CharacteristicLengthMax = 50;
Mesh.Algorithm3D = 1;  // Delaunay
```

Generate mesh:
```bash
gmsh -3 reservoir.geo -o reservoir.msh
```

### Complex Geometries

For complex geometries with faults or fractures:

```geo
// Reservoir with fault
SetFactory("OpenCASCADE");

// Main reservoir block
Box(1) = {0, 0, -100, 1000, 500, 100};

// Fault plane (thin layer)
Box(2) = {400, 0, -100, 5, 500, 100};

// Boolean fragment to connect meshes
BooleanFragments{ Volume{1}; Delete; }{ Volume{2}; Delete; }

// Physical groups
Physical Volume("matrix") = {1};
Physical Volume("fault") = {2};
Physical Surface("boundary") = {1:6};

// Refined mesh near fault
Field[1] = Distance;
Field[1].NNodesByEdge = 100;
Field[1].FacesList = {7:12};  // Fault faces

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = 5;    // Fine near fault
Field[2].LcMax = 50;   // Coarse far away
Field[2].DistMin = 10;
Field[2].DistMax = 200;

Background Field = 2;
```

## API Usage

### C++ Example

```cpp
#include "GmshIO.hpp"
#include "Simulator.hpp"

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    
    MPI_Comm comm = PETSC_COMM_WORLD;
    
    // Load mesh directly
    FSRM::GmshIO gmsh;
    if (!gmsh.readMshFile("meshes/reservoir.msh")) {
        std::cerr << "Failed to load mesh\n";
        return 1;
    }
    
    // Print statistics
    gmsh.printStatistics();
    
    // Create DMPlex for PETSc
    DM dm;
    gmsh.createDMPlex(comm, &dm);
    
    // Use in simulation...
    
    // Cleanup
    DMDestroy(&dm);
    PetscFinalize();
    return 0;
}
```

### Accessing Mesh Data

```cpp
FSRM::GmshIO gmsh;
gmsh.readMshFile("mesh.msh");

auto mesh = gmsh.getMesh();

// Access nodes
for (const auto& node : mesh->nodes) {
    std::cout << "Node " << node.tag 
              << ": (" << node.x << ", " << node.y << ", " << node.z << ")\n";
}

// Access elements
for (const auto& elem : mesh->elements) {
    std::cout << "Element " << elem.tag 
              << " type: " << static_cast<int>(elem.type)
              << " physical: " << elem.physical_tag << "\n";
}

// Get elements by physical group
std::vector<int> inlet_faces = mesh->getElementsByPhysicalTag(1);
```

## Best Practices

1. **Physical Groups**: Always define physical groups in Gmsh for boundaries and regions
2. **Mesh Quality**: Use Gmsh's mesh optimization tools for better element quality
3. **Refinement**: Refine near wells, fractures, and regions of interest
4. **Partitioning**: Let PETSc handle mesh partitioning for parallel runs
5. **Validation**: Use `gmsh.validate()` to check mesh consistency

## Troubleshooting

### Common Issues

1. **"Failed to read Gmsh mesh file"**
   - Check file path
   - Verify MSH format version (v2.2 or v4.x)
   - Ensure file is ASCII (binary not yet supported)

2. **"No physical groups found"**
   - Define Physical Volume and Physical Surface in Gmsh
   - These are required for material and boundary assignment

3. **Poor solver convergence**
   - Check mesh quality with `gmsh -check mesh.msh`
   - Look for degenerate elements or high aspect ratios
   - Refine mesh in problem areas

## See Also

- [Gmsh Documentation](https://gmsh.info/doc/texinfo/gmsh.html)
- [Coordinate Systems](COORDINATE_SYSTEMS.md)
- [Configuration Reference](CONFIGURATION.md)
