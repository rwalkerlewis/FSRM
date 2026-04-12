// Minimal 10x10x10 box for Gmsh import test
// Generate: gmsh -3 meshes/test_box_gmsh.geo -format msh2 -o meshes/test_box_gmsh.msh

SetFactory("OpenCASCADE");

Lx = 10;
Ly = 10;
Lz = 10;

lc = 5;

Box(1) = {0, 0, 0, Lx, Ly, Lz};

Physical Volume("rock", 1) = {1};

// Boundary selection by bounding boxes
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
