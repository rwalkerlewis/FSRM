# Complete Guide to Gmsh Mesh Input for FSRM

This guide provides detailed instructions on creating, configuring, and loading Gmsh meshes into FSRM simulations. It covers the complete workflow from geometry creation to simulation setup.

## Table of Contents

1. [Overview](#overview)
2. [Prerequisites](#prerequisites)
3. [Creating a Gmsh Mesh](#creating-a-gmsh-mesh)
4. [Physical Groups](#physical-groups)
5. [Material Domains](#material-domains)
6. [Fault Surfaces](#fault-surfaces)
7. [Boundary Surfaces](#boundary-surfaces)
8. [Configuration File Setup](#configuration-file-setup)
9. [C++ API Usage](#c-api-usage)
10. [Best Practices](#best-practices)
11. [Troubleshooting](#troubleshooting)
12. [Complete Examples](#complete-examples)

---

## Overview

FSRM supports unstructured meshes in Gmsh MSH format (versions 2.2 and 4.x). Using Gmsh meshes enables:

- **Complex geometries**: Model realistic reservoir shapes, unconventional boundaries
- **Multiple material domains**: Different rock types in different regions
- **Fault representation**: Explicit fault surfaces for seismicity modeling
- **Local refinement**: Higher resolution near wells, faults, and regions of interest
- **CAD integration**: Import complex geometries from engineering software

### Supported Features

| Feature | Support Level |
|---------|---------------|
| MSH v2.2 ASCII | Full |
| MSH v4.x ASCII | Full |
| MSH Binary | Not yet |
| Triangles (2D) | Full |
| Quadrilaterals (2D) | Full |
| Tetrahedra (3D) | Full |
| Hexahedra (3D) | Full |
| Prisms, Pyramids | Full |
| Second-order elements | Full |
| Physical Groups | Full |
| Parallel Distribution | Via PETSc |

---

## Prerequisites

### Installing Gmsh

**Linux (Ubuntu/Debian):**
```bash
sudo apt-get install gmsh
```

**macOS:**
```bash
brew install gmsh
```

**From Source:**
```bash
git clone https://gitlab.onelab.info/gmsh/gmsh.git
cd gmsh
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=/usr/local ..
make -j$(nproc)
sudo make install
```

### Verify Installation
```bash
gmsh --version
```

---

## Creating a Gmsh Mesh

### Method 1: Gmsh GUI

1. Launch Gmsh: `gmsh`
2. Create geometry using **Modules → Geometry**
3. Define Physical Groups (see below)
4. Generate mesh: **Mesh → 3D**
5. Save: **File → Save As → my_mesh.msh**

### Method 2: Gmsh Scripting (.geo files)

Create a `.geo` script file and generate the mesh:

```bash
gmsh -3 my_model.geo -o my_mesh.msh
```

Options:
- `-3`: Generate 3D mesh
- `-2`: Generate 2D mesh
- `-order 2`: Generate second-order elements
- `-format msh2`: Use MSH 2.2 format
- `-format msh4`: Use MSH 4.x format (default)

### Basic .geo Syntax

```geo
// Set geometry kernel (OpenCASCADE recommended)
SetFactory("OpenCASCADE");

// Create primitive shapes
Box(1) = {x0, y0, z0, dx, dy, dz};
Cylinder(2) = {x, y, z, ax, ay, az, r};
Sphere(3) = {x, y, z, r};

// Boolean operations
BooleanUnion{ Volume{1}; Delete; }{ Volume{2}; Delete; }
BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Delete; }
BooleanIntersection{ Volume{1}; Delete; }{ Volume{2}; Delete; }
BooleanFragments{ Volume{1}; Delete; }{ Volume{2}; Delete; }

// Mesh size control
Mesh.CharacteristicLengthMin = 10;
Mesh.CharacteristicLengthMax = 100;
Mesh.Algorithm3D = 1;  // 1=Delaunay, 4=Frontal, 7=MMG3D
```

---

## Physical Groups

Physical groups are **essential** for FSRM. They define:
- Material regions (rock types)
- Fault surfaces
- Boundary conditions

### Defining Physical Groups in .geo Files

```geo
// Physical Volumes (3D material domains)
Physical Volume("reservoir") = {1};
Physical Volume("caprock") = {2};
Physical Volume("aquifer") = {3};

// Physical Surfaces (boundaries and faults)
Physical Surface("inlet") = {10};
Physical Surface("outlet") = {11};
Physical Surface("top") = {12};
Physical Surface("bottom") = {13};
Physical Surface("main_fault") = {20, 21};  // Multiple surfaces
Physical Surface("secondary_fault") = {22};

// Physical Lines (for 2D problems)
Physical Curve("well_trajectory") = {100};
```

### Naming Conventions

Use clear, descriptive names:
- Material domains: `reservoir`, `caprock`, `basement`, `salt_dome`
- Boundaries: `inlet`, `outlet`, `top`, `bottom`, `north`, `south`
- Faults: `main_fault`, `fault_1`, `san_andreas_fault`

### Verifying Physical Groups

After creating a mesh, verify physical groups:
```bash
gmsh my_mesh.msh -info
```

Or in the GUI: **Tools → Visibility → Physical groups**

---

## Material Domains

Material domains allow different material properties in different mesh regions.

### Step 1: Define Physical Volumes in Gmsh

```geo
SetFactory("OpenCASCADE");

// Create a layered reservoir model
Box(1) = {0, 0, -1000, 5000, 5000, 200};  // Reservoir
Box(2) = {0, 0, -800, 5000, 5000, 800};   // Caprock
Box(3) = {0, 0, -1200, 5000, 5000, 200};  // Aquifer

// Fragment to ensure conforming mesh
BooleanFragments{ Volume{1,2,3}; Delete; }{}

// Assign physical groups
Physical Volume("reservoir") = {1};
Physical Volume("caprock") = {2};
Physical Volume("aquifer") = {3};

// Mesh settings
Mesh.CharacteristicLengthMin = 50;
Mesh.CharacteristicLengthMax = 200;
```

### Step 2: Configure Material Mapping in FSRM

In your `.config` file:

```ini
[GRID]
mesh_type = GMSH
mesh_file = meshes/layered_reservoir.msh
use_unstructured = true

# Map physical groups to material sections
# Format: physical_group:MATERIAL_SECTION, ...
gmsh_material_mapping = reservoir:ROCK1, caprock:ROCK2, aquifer:ROCK3

# Define boundaries from physical surfaces
gmsh_boundaries = inlet, outlet, top, bottom

[ROCK1]
name = Reservoir_Sandstone
porosity = 0.20
permeability_x = 100 mD
permeability_y = 100 mD
permeability_z = 10 mD
youngs_modulus = 15 GPa
poisson_ratio = 0.25
density = 2500 kg/m3

[ROCK2]
name = Caprock_Shale
porosity = 0.05
permeability_x = 0.001 mD
permeability_y = 0.001 mD
permeability_z = 0.0001 mD
youngs_modulus = 25 GPa
poisson_ratio = 0.30
density = 2650 kg/m3

[ROCK3]
name = Aquifer_Limestone
porosity = 0.15
permeability_x = 50 mD
permeability_y = 50 mD
permeability_z = 5 mD
youngs_modulus = 40 GPa
poisson_ratio = 0.28
density = 2700 kg/m3
```

---

## Fault Surfaces

Faults are defined as physical surfaces in Gmsh and can be configured for:
- Friction-based slip behavior
- Rate-and-state friction laws
- Seismicity modeling with split nodes

### Creating Fault Geometry in Gmsh

#### Method 1: Thin Layer (Recommended for Simple Faults)

```geo
SetFactory("OpenCASCADE");

// Main domain
Box(1) = {0, 0, -2000, 10000, 5000, 2000};

// Fault as thin layer
// Strike: N-S, Dip: 60° to the east
Box(2) = {4500, 0, -2000, 100, 5000, 2000};

// Fragment for conforming mesh
BooleanFragments{ Volume{1}; Delete; }{ Volume{2}; Delete; }

// Identify fault surfaces (internal faces of fault layer)
Physical Volume("matrix") = {1, 3};  // Material on both sides
Physical Volume("fault_zone") = {2}; // Fault layer
Physical Surface("fault_surface") = {7, 8};  // Internal interfaces
```

#### Method 2: Embedded Surface (For Discontinuous Faults)

```geo
SetFactory("OpenCASCADE");

// Main domain
Box(1) = {0, 0, -3000, 10000, 10000, 3000};

// Fault plane as rectangle
Point(100) = {3000, 0, -3000};
Point(101) = {7000, 0, -3000};
Point(102) = {7000, 10000, -1000};
Point(103) = {3000, 10000, -1000};

Line(100) = {100, 101};
Line(101) = {101, 102};
Line(102) = {102, 103};
Line(103) = {103, 100};

Curve Loop(100) = {100, 101, 102, 103};
Plane Surface(100) = {100};

// Embed fault in volume
BooleanFragments{ Volume{1}; Delete; }{ Surface{100}; Delete; }

// Mark fault surface
Physical Surface("main_fault") = {100};
Physical Volume("rock") = {1};
```

### Configuring Fault Properties

```ini
[GRID]
mesh_type = GMSH
mesh_file = meshes/faulted_domain.msh
use_unstructured = true

# Map fault surfaces to fault configurations
# Format: physical_group:FAULT_SECTION[:split]
# Add :split to enable split nodes for discontinuous displacement
gmsh_fault_mapping = main_fault:FAULT1:split, secondary_fault:FAULT2

gmsh_material_mapping = rock:ROCK1, fault_zone:ROCK_FAULT
gmsh_boundaries = top, bottom, north, south, east, west

[FAULT1]
name = Main_Fault
friction_law = RATE_STATE_AGING
static_friction = 0.6
dynamic_friction = 0.4
cohesion = 1 MPa

# Rate-state parameters
rate_state_a = 0.010
rate_state_b = 0.015    # b > a = velocity weakening (unstable)
rate_state_dc = 0.1 mm
rate_state_v0 = 1e-6 m/s
rate_state_f0 = 0.6

# Split node configuration
use_split_nodes = true
split_node_method = PENALTY
penalty_normal = 1e12
penalty_tangent = 1e10

[FAULT2]
name = Secondary_Fault
friction_law = COULOMB
static_friction = 0.65
dynamic_friction = 0.45
cohesion = 500 kPa
use_split_nodes = false
```

---

## Boundary Surfaces

### Defining Boundary Conditions

```geo
SetFactory("OpenCASCADE");

Box(1) = {0, 0, 0, 1000, 500, 100};

// Get boundary surfaces
Physical Surface("west") = {1};    // x = 0
Physical Surface("east") = {2};    // x = 1000
Physical Surface("south") = {3};   // y = 0
Physical Surface("north") = {4};   // y = 500
Physical Surface("bottom") = {5};  // z = 0
Physical Surface("top") = {6};     // z = 100

Physical Volume("domain") = {1};
```

### Applying Boundary Conditions in Config

```ini
[GRID]
mesh_type = GMSH
mesh_file = meshes/reservoir.msh
gmsh_boundaries = west, east, south, north, bottom, top
gmsh_material_mapping = domain:ROCK1

# Boundary conditions reference physical group names
[BC1]
type = DIRICHLET
field = PRESSURE
location = west              # Reference physical group name
value = 20 MPa

[BC2]
type = DIRICHLET
field = PRESSURE
location = east
value = 15 MPa

[BC3]
type = NEUMANN
field = PRESSURE
location = north, south      # Multiple boundaries
value = 0                    # No-flow

[BC4]
type = DIRICHLET
field = DISPLACEMENT
location = bottom
value = 0, 0, 0              # Fixed base
```

---

## Configuration File Setup

### Complete Configuration Example

```ini
# =============================================================================
# FSRM Configuration with Gmsh Mesh
# =============================================================================

[SIMULATION]
name = faulted_reservoir
start_time = 0.0
end_time = 365 days
dt_initial = 1 hour
dt_min = 1 second
dt_max = 1 day
max_timesteps = 100000

fluid_model = SINGLE_COMPONENT
solid_model = POROELASTIC

enable_geomechanics = true
enable_faults = true
enable_thermal = false

rtol = 1.0e-6
atol = 1.0e-8

[GRID]
# Unstructured mesh from Gmsh
mesh_type = GMSH
mesh_file = meshes/complex_reservoir.msh
use_unstructured = true

# Material domain mappings
gmsh_material_mapping = reservoir:ROCK1, caprock:ROCK2, basement:ROCK3

# Fault surface mappings
gmsh_fault_mapping = main_fault:FAULT1:split, branch_fault:FAULT2

# Boundary surface names
gmsh_boundaries = inlet, outlet, top, bottom, sides

# Coordinate system (optional)
input_crs = EPSG:4326
model_crs = EPSG:32610
use_local_coordinates = true
local_origin_x = -122.4
local_origin_y = 37.8

[ROCK1]
name = Reservoir
porosity = 0.20
permeability_x = 100 mD
permeability_y = 100 mD
permeability_z = 10 mD
youngs_modulus = 15 GPa
poisson_ratio = 0.25
density = 2500 kg/m3
biot_coefficient = 0.8

[ROCK2]
name = Caprock
porosity = 0.05
permeability_x = 0.001 mD
permeability_y = 0.001 mD
permeability_z = 0.0001 mD
youngs_modulus = 25 GPa
poisson_ratio = 0.30
density = 2650 kg/m3

[ROCK3]
name = Basement
porosity = 0.02
permeability_x = 0.0001 mD
permeability_y = 0.0001 mD
permeability_z = 0.00001 mD
youngs_modulus = 60 GPa
poisson_ratio = 0.22
density = 2800 kg/m3

[FLUID]
type = SINGLE_PHASE
density = 1000 kg/m3
viscosity = 1 cP
compressibility = 4.5e-10 1/Pa

[FAULT1]
name = Main_Fault
friction_law = RATE_STATE_AGING
static_friction = 0.6
dynamic_friction = 0.4
cohesion = 1 MPa
rate_state_a = 0.010
rate_state_b = 0.015
rate_state_dc = 0.0001 m
use_split_nodes = true
split_node_method = PENALTY

[FAULT2]
name = Branch_Fault
friction_law = COULOMB
static_friction = 0.6
dynamic_friction = 0.4
cohesion = 500 kPa

[WELL1]
name = INJECTOR
type = INJECTOR
x = 2000
y = 2500
z = -1500
control_mode = RATE
target_value = 0.01 m3/s

[BC1]
type = DIRICHLET
field = PRESSURE
location = inlet
value = 25 MPa

[BC2]
type = DIRICHLET
field = PRESSURE
location = outlet
value = 15 MPa

[BC3]
type = NEUMANN
field = PRESSURE
location = sides, top, bottom
value = 0

[IC1]
field = PRESSURE
distribution = GRADIENT
value = 20 MPa
gradient = 0, 0, 10000

[OUTPUT]
format = VTK
path = output/faulted_reservoir
frequency = 10
pressure = true
displacement = true
stress = true
fault_slip = true
```

---

## C++ API Usage

### Loading and Processing a Gmsh Mesh

```cpp
#include "GmshIO.hpp"
#include "Simulator.hpp"

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    
    // Load Gmsh mesh
    FSRM::GmshIO gmsh;
    if (!gmsh.readMshFile("meshes/reservoir.msh")) {
        std::cerr << "Failed to load mesh\n";
        return 1;
    }
    
    // Print mesh statistics
    gmsh.printStatistics();
    
    // Set up material domain mappings
    gmsh.addMaterialMapping("reservoir", "ROCK1", 0);
    gmsh.addMaterialMapping("caprock", "ROCK2", 1);
    gmsh.addMaterialMapping("basement", "ROCK3", 2);
    
    // Set up fault mappings
    gmsh.addFaultMapping("main_fault", "FAULT1", true);  // use split nodes
    gmsh.addFaultMapping("branch_fault", "FAULT2", false);
    
    // Set boundary names
    gmsh.setBoundaryNames({"inlet", "outlet", "top", "bottom", "sides"});
    
    // Process all domain mappings
    if (!gmsh.processDomains()) {
        std::cerr << "Failed to process domains\n";
        return 1;
    }
    
    // Print domain information
    gmsh.printDomainInfo();
    
    // Validate mesh quality
    double min_q, max_q, avg_q;
    gmsh.getQualityStats(min_q, max_q, avg_q);
    std::cout << "Mesh quality: min=" << min_q 
              << ", max=" << max_q 
              << ", avg=" << avg_q << "\n";
    
    // Create PETSc DMPlex with all labels
    DM dm;
    gmsh.createDMPlexFull(PETSC_COMM_WORLD, &dm);
    
    // Access mesh data
    auto mesh = gmsh.getMesh();
    
    // Iterate over material domains
    for (const auto& name : mesh->getMaterialDomainNames()) {
        auto* domain = mesh->getMaterialDomain(name);
        std::cout << "Domain: " << name << "\n";
        std::cout << "  Volume: " << domain->volume << " m³\n";
        std::cout << "  Elements: " << domain->num_elements << "\n";
    }
    
    // Access fault surfaces
    for (const auto& name : mesh->getFaultNames()) {
        auto* fault = mesh->getFaultSurface(name);
        std::cout << "Fault: " << name << "\n";
        std::cout << "  Area: " << fault->area << " m²\n";
        std::cout << "  Strike: " << fault->strike * 180/M_PI << "°\n";
        std::cout << "  Dip: " << fault->dip * 180/M_PI << "°\n";
    }
    
    // Write VTK output with material IDs
    gmsh.writeVtkFileWithMaterials("output/mesh_with_materials.vtu");
    
    // Cleanup
    DMDestroy(&dm);
    PetscFinalize();
    return 0;
}
```

### Integrating with Simulator

```cpp
#include "Simulator.hpp"
#include "ConfigReader.hpp"

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    
    FSRM::Simulator sim(PETSC_COMM_WORLD);
    
    // Initialize from config file (handles Gmsh mesh loading automatically)
    PetscErrorCode ierr = sim.initializeFromConfigFile("config/my_simulation.config");
    CHKERRQ(ierr);
    
    // Run simulation
    ierr = sim.run();
    CHKERRQ(ierr);
    
    PetscFinalize();
    return 0;
}
```

---

## Best Practices

### Mesh Quality

1. **Check element quality** before simulation:
   ```bash
   gmsh my_mesh.msh -check
   ```

2. **Optimize mesh** after generation:
   ```geo
   Mesh.Optimize = 1;
   Mesh.OptimizeNetgen = 1;
   ```

3. **Target quality metrics**:
   - Minimum element quality > 0.1
   - Average quality > 0.5
   - Maximum aspect ratio < 10

### Mesh Refinement

1. **Use fields for local refinement**:
   ```geo
   // Distance-based refinement near fault
   Field[1] = Distance;
   Field[1].SurfacesList = {fault_surface};
   
   Field[2] = Threshold;
   Field[2].IField = 1;
   Field[2].LcMin = 10;      // Fine mesh size
   Field[2].LcMax = 100;     // Coarse mesh size
   Field[2].DistMin = 50;    // Distance for fine mesh
   Field[2].DistMax = 500;   // Distance for coarse mesh
   
   Background Field = 2;
   ```

2. **Refine near wells**:
   ```geo
   Field[3] = Ball;
   Field[3].XCenter = 500;
   Field[3].YCenter = 500;
   Field[3].ZCenter = -1000;
   Field[3].Radius = 100;
   Field[3].VIn = 5;         // Mesh size inside ball
   Field[3].VOut = 50;       // Mesh size outside
   ```

### Physical Groups

1. **Always define physical groups** - elements without physical groups may be ignored
2. **Use descriptive names** - easier to maintain configuration files
3. **Verify group assignments** - ensure all elements are assigned

### File Organization

```
project/
├── meshes/
│   ├── reservoir.geo        # Gmsh geometry script
│   ├── reservoir.msh        # Generated mesh
│   └── README.md            # Mesh documentation
├── config/
│   └── simulation.config    # FSRM configuration
└── output/
    └── ...                  # Simulation results
```

---

## Troubleshooting

### "Failed to read Gmsh mesh file"

**Causes & Solutions:**
- Wrong file path → Check path in config file
- Binary format → Convert to ASCII: `gmsh mesh.msh -o mesh_ascii.msh -format msh2`
- Corrupted file → Regenerate mesh

### "No physical groups found"

**Solutions:**
1. Add physical groups to your .geo file
2. Ensure you're using `Physical Volume` (3D) or `Physical Surface` (2D)
3. Regenerate the mesh after adding physical groups

### "Physical group 'name' not found"

**Solutions:**
1. Check spelling in config file matches .geo file exactly
2. Physical group names are case-sensitive
3. Run `gmsh mesh.msh -info` to list available groups

### "Poor solver convergence"

**Mesh-related causes:**
1. Poor element quality → Run `gmsh mesh.msh -check`
2. Very small elements → Increase minimum mesh size
3. High aspect ratio elements → Use mesh optimization

### "Memory issues with large meshes"

**Solutions:**
1. Use parallel meshing: `gmsh -nt 4 mesh.geo -3`
2. Reduce mesh resolution where possible
3. Use coarser mesh with AMR (adaptive mesh refinement)

---

## Complete Examples

### Example 1: Simple Reservoir with Fault

**File: `simple_faulted.geo`**
```geo
SetFactory("OpenCASCADE");

// Parameters
Lx = 5000;  // Domain length (m)
Ly = 2000;  // Domain width (m)
Lz = 1000;  // Domain height (m)
fault_x = 2500;  // Fault x-position

// Create main domain
Box(1) = {0, 0, -Lz, fault_x, Ly, Lz};
Box(2) = {fault_x, 0, -Lz, Lx-fault_x, Ly, Lz};

// Fragment for conforming mesh at fault
BooleanFragments{ Volume{1,2}; Delete; }{}

// Physical groups
Physical Volume("west_block") = {1};
Physical Volume("east_block") = {2};

// Get fault surface (interface between volumes)
Physical Surface("fault") = {7};  // Check numbering with GUI

// Boundaries
Physical Surface("west") = {1};
Physical Surface("east") = {6};
Physical Surface("south") = {3, 9};
Physical Surface("north") = {4, 10};
Physical Surface("bottom") = {2, 8};
Physical Surface("top") = {5, 11};

// Mesh settings
Mesh.CharacteristicLengthMin = 50;
Mesh.CharacteristicLengthMax = 200;

// Refine near fault
Field[1] = Distance;
Field[1].SurfacesList = {7};

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = 20;
Field[2].LcMax = 200;
Field[2].DistMin = 100;
Field[2].DistMax = 500;

Background Field = 2;
```

**File: `simple_faulted.config`**
```ini
[SIMULATION]
name = simple_faulted_reservoir
end_time = 30 days
dt_initial = 1 hour
fluid_model = SINGLE_COMPONENT
solid_model = POROELASTIC
enable_geomechanics = true
enable_faults = true

[GRID]
mesh_type = GMSH
mesh_file = meshes/simple_faulted.msh
use_unstructured = true
gmsh_material_mapping = west_block:ROCK1, east_block:ROCK1
gmsh_fault_mapping = fault:FAULT1:split
gmsh_boundaries = west, east, south, north, bottom, top

[ROCK1]
porosity = 0.20
permeability_x = 100 mD
youngs_modulus = 15 GPa
poisson_ratio = 0.25

[FAULT1]
friction_law = COULOMB
static_friction = 0.6
dynamic_friction = 0.4
use_split_nodes = true

[FLUID]
density = 1000 kg/m3
viscosity = 1 cP

[BC1]
type = DIRICHLET
field = PRESSURE
location = west
value = 25 MPa

[BC2]
type = DIRICHLET
field = PRESSURE
location = east
value = 20 MPa

[OUTPUT]
format = VTK
path = output/simple_faulted
pressure = true
displacement = true
fault_slip = true
```

**Generate and run:**
```bash
gmsh -3 meshes/simple_faulted.geo -o meshes/simple_faulted.msh
mpirun -np 4 fsrm -c config/simple_faulted.config
```

### Example 2: Multi-Layer Reservoir

See `config/unstructured_gmsh_example.config` for a complete multi-layer example.

---

## See Also

- [Unstructured Meshes](UNSTRUCTURED_MESHES.md) - Overview of unstructured mesh support
- [Configuration Reference](CONFIGURATION.md) - All configuration options
- [Fault Model](PHYSICS_MODELS.md#fault-mechanics) - Fault physics documentation
- [Coordinate Systems](COORDINATE_SYSTEMS.md) - Geographic coordinate handling
- [Gmsh Documentation](https://gmsh.info/doc/texinfo/gmsh.html) - Official Gmsh manual
