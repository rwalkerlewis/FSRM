# Unstructured Mesh Support with Gmsh

FSRM supports unstructured meshes using the Gmsh format, enabling complex reservoir geometries that cannot be represented with structured grids.

## Overview

Gmsh is an open-source finite element mesh generator that produces high-quality tetrahedral, hexahedral, and mixed-element meshes. FSRM can directly load Gmsh `.msh` files and convert them to PETSc DMPlex for simulation.

> **For detailed instructions on creating and configuring Gmsh meshes, see the [Complete Gmsh Mesh Guide](GMSH_MESH_GUIDE.md).**

## Supported Features

| Feature | Status |
|---------|--------|
| MSH v2.2 (ASCII) | ✅ Full support |
| MSH v4.x (ASCII) | ✅ Full support |
| MSH Binary | ⏳ Not yet supported |
| Triangles/Quads (2D) | ✅ Full support |
| Tetrahedra/Hexahedra (3D) | ✅ Full support |
| Prisms/Pyramids | ✅ Full support |
| Second-order elements | ✅ Full support |
| Physical Groups | ✅ Full support |
| **Material Domains** | ✅ Full support |
| **Fault Surfaces** | ✅ Full support |
| **Boundary Surfaces** | ✅ Full support |
| Parallel Distribution | ✅ Via PETSc |
| Mesh Quality Validation | ✅ Built-in |

## Quick Start

### 1. Create a Gmsh Mesh

Create a `.geo` file defining your geometry with physical groups:

```geo
SetFactory("OpenCASCADE");

// Create reservoir geometry
Box(1) = {0, 0, -100, 1000, 500, 100};

// Define physical groups
Physical Volume("reservoir") = {1};
Physical Surface("inlet") = {1};
Physical Surface("outlet") = {2};
Physical Surface("walls") = {3, 4};
Physical Surface("top") = {5};
Physical Surface("bottom") = {6};

// Mesh settings
Mesh.CharacteristicLengthMin = 20;
Mesh.CharacteristicLengthMax = 50;
```

Generate the mesh:
```bash
gmsh -3 reservoir.geo -o reservoir.msh
```

### 2. Configure FSRM

```ini
[GRID]
mesh_type = GMSH
mesh_file = meshes/reservoir.msh
use_unstructured = true

# Map physical groups to materials
gmsh_material_mapping = reservoir:ROCK1

# Specify boundaries
gmsh_boundaries = inlet, outlet, walls, top, bottom

[ROCK1]
porosity = 0.20
permeability_x = 100 mD
permeability_y = 100 mD
permeability_z = 10 mD
youngs_modulus = 15 GPa

[BC1]
type = DIRICHLET
field = PRESSURE
location = inlet
value = 15.0e6

[BC2]
type = NEUMANN
field = PRESSURE
location = walls
value = 0.0
```

### 3. Run Simulation

```bash
mpirun -np 4 fsrm -c config/my_simulation.config
```

## Material Domain Mapping

FSRM allows mapping physical volumes in your Gmsh mesh to different material properties:

```ini
[GRID]
mesh_type = GMSH
mesh_file = meshes/layered_reservoir.msh

# Format: physical_group:MATERIAL_SECTION, ...
gmsh_material_mapping = reservoir:ROCK1, caprock:ROCK2, aquifer:ROCK3

[ROCK1]
name = Reservoir_Sandstone
porosity = 0.20
permeability_x = 100 mD

[ROCK2]
name = Caprock_Shale
porosity = 0.05
permeability_x = 0.001 mD

[ROCK3]
name = Aquifer_Limestone
porosity = 0.15
permeability_x = 50 mD
```

## Fault Surface Mapping

Map physical surfaces to fault models for seismicity simulation:

```ini
[GRID]
# Format: physical_group:FAULT_SECTION[:split]
# Add :split to enable split nodes for discontinuous displacement
gmsh_fault_mapping = main_fault:FAULT1:split, secondary_fault:FAULT2

[FAULT1]
name = Main_Fault
friction_law = RATE_STATE_AGING
static_friction = 0.6
dynamic_friction = 0.4
cohesion = 1 MPa
use_split_nodes = true
split_node_method = PENALTY
```

## Coordinate System Support

Gmsh meshes can use geographic coordinates:

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
local_origin_x = -122.4194
local_origin_y = 37.7749
```

## C++ API

### Basic Usage

```cpp
#include "GmshIO.hpp"

FSRM::GmshIO gmsh;
if (gmsh.readMshFile("mesh.msh")) {
    // Print statistics
    gmsh.printStatistics();
    
    // Set up domain mappings
    gmsh.addMaterialMapping("reservoir", "ROCK1", 0);
    gmsh.addFaultMapping("main_fault", "FAULT1", true);
    gmsh.setBoundaryNames({"inlet", "outlet"});
    
    // Process domains
    gmsh.processDomains();
    gmsh.printDomainInfo();
    
    // Create DMPlex
    DM dm;
    gmsh.createDMPlexFull(PETSC_COMM_WORLD, &dm);
}
```

### Accessing Mesh Data

```cpp
auto mesh = gmsh.getMesh();

// Get material domain info
for (const auto& name : mesh->getMaterialDomainNames()) {
    auto* domain = mesh->getMaterialDomain(name);
    std::cout << "Domain: " << name 
              << ", Volume: " << domain->volume << " m³\n";
}

// Get fault surface info
for (const auto& name : mesh->getFaultNames()) {
    auto* fault = mesh->getFaultSurface(name);
    std::cout << "Fault: " << name 
              << ", Area: " << fault->area << " m²"
              << ", Strike: " << fault->strike * 180/M_PI << "°\n";
}

// Get boundary info
for (const auto& name : mesh->getBoundaryNames()) {
    auto nodes = mesh->getBoundaryNodes(name);
    std::cout << "Boundary: " << name 
              << ", Nodes: " << nodes.size() << "\n";
}
```

### Mesh Quality Validation

```cpp
// Get quality statistics
double min_q, max_q, avg_q;
gmsh.getQualityStats(min_q, max_q, avg_q);

// Validate quality
if (!gmsh.validateQuality(0.1, 10.0)) {
    std::cerr << "Warning: Some elements have poor quality\n";
}

// Write VTK with material IDs for visualization
gmsh.writeVtkFileWithMaterials("mesh_with_materials.vtu");
```

## Creating Complex Meshes

### Reservoir with Fault

```geo
SetFactory("OpenCASCADE");

// Reservoir with fault zone
Box(1) = {0, 0, -100, 1000, 500, 100};
Box(2) = {400, 0, -100, 5, 500, 100};  // Fault zone

BooleanFragments{ Volume{1}; Delete; }{ Volume{2}; Delete; }

Physical Volume("matrix") = {1, 3};
Physical Volume("fault_zone") = {2};
Physical Surface("fault_interface") = {7, 8};

// Refined mesh near fault
Field[1] = Distance;
Field[1].FacesList = {7, 8};

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = 5;
Field[2].LcMax = 50;
Field[2].DistMin = 10;
Field[2].DistMax = 200;

Background Field = 2;
```

### Multi-Layer System

```geo
SetFactory("OpenCASCADE");

// Layers
Box(1) = {0, 0, -200, 2000, 2000, 50};   // Caprock
Box(2) = {0, 0, -250, 2000, 2000, 50};   // Reservoir
Box(3) = {0, 0, -400, 2000, 2000, 150};  // Aquifer

BooleanFragments{ Volume{1,2,3}; Delete; }{}

Physical Volume("caprock") = {1};
Physical Volume("reservoir") = {2};
Physical Volume("aquifer") = {3};
Physical Surface("top") = {...};
Physical Surface("bottom") = {...};
Physical Surface("sides") = {...};
```

## Best Practices

1. **Always define physical groups** - Required for material and boundary assignment
2. **Use mesh refinement near features** - Wells, faults, and regions of interest
3. **Validate mesh quality** - Use `gmsh -check mesh.msh`
4. **Optimize the mesh** - Add `Mesh.Optimize = 1;` to .geo files
5. **Use ASCII format** - Binary not yet supported

## Troubleshooting

| Problem | Solution |
|---------|----------|
| "Failed to read mesh file" | Check file path; ensure ASCII format |
| "No physical groups found" | Add Physical Volume/Surface to .geo file |
| "Physical group not found" | Check spelling matches exactly |
| Poor solver convergence | Improve mesh quality; reduce aspect ratio |
| Memory issues | Use coarser mesh or parallel execution |

## See Also

- [Complete Gmsh Mesh Guide](GMSH_MESH_GUIDE.md) - Detailed mesh creation tutorial
- [Configuration Reference](CONFIGURATION.md) - All configuration options
- [Coordinate Systems](COORDINATE_SYSTEMS.md) - Geographic coordinate handling
- [Fault Model Documentation](PHYSICS_MODELS.md#fault-mechanics)
- [Gmsh Official Documentation](https://gmsh.info/doc/texinfo/gmsh.html)
