// ============================================================================
// Binary Gmsh mesh example (for PETSc's Gmsh reader in FSRM)
// ============================================================================
// Generate a *binary* .msh (recommended):
//   gmsh -3 meshes/binary_box/binary_box.geo -format msh4 -bin 1 -o meshes/binary_box/binary_box_3d_bin.msh
//
// Or generate ASCII:
//   gmsh -3 meshes/binary_box/binary_box.geo -format msh4 -bin 0 -o meshes/binary_box/binary_box_3d_ascii.msh
// ============================================================================

SetFactory("OpenCASCADE");

Lx = 1000;
Ly =  500;
Lz =  100;

lc = 50;

Box(1) = {0, 0, 0, Lx, Ly, Lz};

// Physical groups used by config/gmsh_binary_box_example.config
Physical Volume("reservoir", 1) = {1};

// Robust boundary selection by bounding boxes (stable across OpenCASCADE tag changes)
eps = 1e-6;
xmin[] = Surface In BoundingBox{-eps, -eps, -eps, eps, Ly + eps, Lz + eps};
xmax[] = Surface In BoundingBox{Lx - eps, -eps, -eps, Lx + eps, Ly + eps, Lz + eps};
ymin[] = Surface In BoundingBox{-eps, -eps, -eps, Lx + eps, eps, Lz + eps};
ymax[] = Surface In BoundingBox{-eps, Ly - eps, -eps, Lx + eps, Ly + eps, Lz + eps};
zmin[] = Surface In BoundingBox{-eps, -eps, -eps, Lx + eps, Ly + eps, eps};
zmax[] = Surface In BoundingBox{-eps, -eps, Lz - eps, Lx + eps, Ly + eps, Lz + eps};

Physical Surface("xmin") = {xmin[]};
Physical Surface("xmax") = {xmax[]};
Physical Surface("ymin") = {ymin[]};
Physical Surface("ymax") = {ymax[]};
Physical Surface("zmin") = {zmin[]};
Physical Surface("zmax") = {zmax[]};

Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;
Mesh.MshFileVersion = 4.1;
