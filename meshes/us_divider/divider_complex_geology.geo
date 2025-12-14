// ============================================================================
// Divider (1992) - Complex Geological Model + Detailed Gmsh Mesh (FSRM example)
// ============================================================================
// Purpose
//   - Provide a *rich* Gmsh model for a “Divider-style” contained underground shot.
//   - Demonstrate: multi-layer stratigraphy, a “divider” fault (hydraulic/mechanical
//     barrier), a secondary fault, and an intrusive/salt-like body.
//   - Provide stable Physical Groups (names) for FSRM material + fault mappings.
//
// Generate:
//   gmsh -3 meshes/us_divider/divider_complex_geology.geo \
//        -format msh4 -bin 0 -o meshes/us_divider/divider_complex_geology_3d.msh
//
// Notes
//   - Keep the mesh ASCII (FSRM does not support binary .msh).
//   - The physical group names are used verbatim in the config file.
// ============================================================================

SetFactory("OpenCASCADE");

// ----------------------------
// Geometry: domain + layers
// ----------------------------
Lx = 12000;   // m
Ly =  8000;   // m
zTop =    0;  // m (surface)
zA   = -200;  // alluvium base
zT   = -800;  // tuff base
zWT  = -1500; // welded tuff base
zC   = -2500; // carbonate base
zB   = -4000; // basement base (domain bottom)

// Main “divider” fault plane is vertical at x = xDiv.
// This gives a robust selection strategy (slab by X) for Physical Surface.
xDiv = 6000;

// Secondary fault is also vertical (parallel) for robust Physical Surface selection.
xF2 = 9000;

// Characteristic lengths (global)
lc_far   = 400;  // coarse
lc_mid   = 200;  // medium
lc_near  = 80;   // near-source + near-fault

// Layer boxes (stacked). Boxes are defined with {x,y,z, dx,dy,dz}.
// Note: dz must be positive.
Box(1) = {0, 0, zA,  Lx, Ly, zTop - zA};      // alluvium: 0 .. -200
Box(2) = {0, 0, zT,  Lx, Ly, zA   - zT};      // tuff:     -200 .. -800
Box(3) = {0, 0, zWT, Lx, Ly, zT   - zWT};     // welded:   -800 .. -1500
Box(4) = {0, 0, zC,  Lx, Ly, zWT  - zC};      // carbonate:-1500 .. -2500
Box(5) = {0, 0, zB,  Lx, Ly, zC   - zB};      // basement: -2500 .. -4000

// ----------------------------
// Intrusive / salt-like body
// ----------------------------
// An ellipsoidal “dome/stock” that cuts multiple layers.
// We keep it off-center and mostly on the downstream side to enforce asymmetry.
dome_x = 8500;
dome_y = 4200;
dome_z = -1800;
dome_r = 700;

Sphere(50) = {dome_x, dome_y, dome_z, dome_r};
// Turn sphere into an ellipsoid (wider in X, slightly compressed in Z)
Dilate {{dome_x, dome_y, dome_z}, 1.25, 1.00, 0.80} { Volume{50}; }

// ----------------------------
// Fault planes (as surfaces)
// ----------------------------
// Divider fault: x = xDiv, spanning full y and z range.
Point(1000) = {xDiv, 0,  zB, lc_mid};
Point(1001) = {xDiv, Ly, zB, lc_mid};
Point(1002) = {xDiv, Ly, zTop, lc_mid};
Point(1003) = {xDiv, 0,  zTop, lc_mid};

Line(1000) = {1000, 1001};
Line(1001) = {1001, 1002};
Line(1002) = {1002, 1003};
Line(1003) = {1003, 1000};
Curve Loop(1000) = {1000, 1001, 1002, 1003};
Plane Surface(1000) = {1000};

// Secondary fault: x = xF2
Point(1100) = {xF2, 0,  zB, lc_mid};
Point(1101) = {xF2, Ly, zB, lc_mid};
Point(1102) = {xF2, Ly, zTop, lc_mid};
Point(1103) = {xF2, 0,  zTop, lc_mid};

Line(1100) = {1100, 1101};
Line(1101) = {1101, 1102};
Line(1102) = {1102, 1103};
Line(1103) = {1103, 1100};
Curve Loop(1100) = {1100, 1101, 1102, 1103};
Plane Surface(1100) = {1100};

// ----------------------------
// Boolean fragmentation
// ----------------------------
// Split layers by: (1) divider fault, (2) secondary fault, (3) dome volume.
// This produces many conforming sub-volumes while retaining clean interfaces.
BooleanFragments{ Volume{1,2,3,4,5}; Delete; }{ Volume{50}; Surface{1000,1100}; }{};

// ----------------------------
// Physical groups: volumes
// ----------------------------
// Strategy: robust selection by bounding boxes.
// We produce upstream (x < xDiv) and downstream (x > xDiv) material domains
// for each layer, plus the intrusive dome volume.

bb_eps = 0.5; // meters

// Upstream side (west of divider fault)
Physical Volume("alluvium_upstream") = Volume In BoundingBox{0, 0, zA,   xDiv - bb_eps, Ly, zTop};
Physical Volume("tuff_upstream")     = Volume In BoundingBox{0, 0, zT,   xDiv - bb_eps, Ly, zA};
Physical Volume("welded_tuff_upstream") = Volume In BoundingBox{0, 0, zWT, xDiv - bb_eps, Ly, zT};
Physical Volume("carbonate_upstream")   = Volume In BoundingBox{0, 0, zC,  xDiv - bb_eps, Ly, zWT};
Physical Volume("basement_upstream")    = Volume In BoundingBox{0, 0, zB,  xDiv - bb_eps, Ly, zC};

// Downstream side (east of divider fault)
Physical Volume("alluvium_downstream") = Volume In BoundingBox{xDiv + bb_eps, 0, zA,   Lx, Ly, zTop};
Physical Volume("tuff_downstream")     = Volume In BoundingBox{xDiv + bb_eps, 0, zT,   Lx, Ly, zA};
Physical Volume("welded_tuff_downstream") = Volume In BoundingBox{xDiv + bb_eps, 0, zWT, Lx, Ly, zT};
Physical Volume("carbonate_downstream")   = Volume In BoundingBox{xDiv + bb_eps, 0, zC,  Lx, Ly, zWT};
Physical Volume("basement_downstream")    = Volume In BoundingBox{xDiv + bb_eps, 0, zB,  Lx, Ly, zC};

// Intrusive dome volume (tight bounding box around the ellipsoid)
Physical Volume("intrusive_dome") = Volume In BoundingBox{
  dome_x - 1.35*dome_r, dome_y - 1.10*dome_r, dome_z - 0.95*dome_r,
  dome_x + 1.35*dome_r, dome_y + 1.10*dome_r, dome_z + 0.95*dome_r
};

// ----------------------------
// Physical groups: boundaries
// ----------------------------
// Name boundaries for easy BC application in FSRM.
// “upstream” and “downstream” are the west/east outer faces.

Physical Surface("surface")  = Surface In BoundingBox{0, 0, zTop - bb_eps, Lx, Ly, zTop + bb_eps};
Physical Surface("bottom")   = Surface In BoundingBox{0, 0, zB   - bb_eps, Lx, Ly, zB   + bb_eps};
Physical Surface("upstream") = Surface In BoundingBox{0 - bb_eps, 0, zB, 0 + bb_eps, Ly, zTop};
Physical Surface("downstream") = Surface In BoundingBox{Lx - bb_eps, 0, zB, Lx + bb_eps, Ly, zTop};
Physical Surface("south")    = Surface In BoundingBox{0, 0 - bb_eps, zB, Lx, 0 + bb_eps, zTop};
Physical Surface("north")    = Surface In BoundingBox{0, Ly - bb_eps, zB, Lx, Ly + bb_eps, zTop};

// ----------------------------
// Physical groups: faults
// ----------------------------
// Internal fault surfaces (use a tight slab in X so the groups are stable).
Physical Surface("divider_fault") = Surface In BoundingBox{xDiv - bb_eps, 0, zB, xDiv + bb_eps, Ly, zTop};
Physical Surface("secondary_fault") = Surface In BoundingBox{xF2 - bb_eps, 0, zB, xF2 + bb_eps, Ly, zTop};

// ----------------------------
// Mesh sizing fields
// ----------------------------
// Near-source refinement around a representative ground-zero location.
// (Used by the config’s [EXPLOSION_SOURCE] location.)
GZx = 5600;
GZy = 4000;
GZz = -800;

// Refine near the explosion point
Field[1] = Ball;
Field[1].XCenter = GZx;
Field[1].YCenter = GZy;
Field[1].ZCenter = GZz;
Field[1].Radius  = 900;
Field[1].VIn     = lc_near;
Field[1].VOut    = lc_mid;

// Refine near the divider fault plane
Field[2] = Distance;
Field[2].SurfacesList = {1000};

Field[3] = Threshold;
Field[3].IField = 2;
Field[3].LcMin = lc_near;
Field[3].LcMax = lc_far;
Field[3].DistMin = 150;
Field[3].DistMax = 1200;

// Refine near the dome (optional additional heterogeneity focus)
Field[4] = Ball;
Field[4].XCenter = dome_x;
Field[4].YCenter = dome_y;
Field[4].ZCenter = dome_z;
Field[4].Radius  = 1200;
Field[4].VIn     = lc_mid;
Field[4].VOut    = lc_far;

Field[10] = Min;
Field[10].FieldsList = {1, 3, 4};
Background Field = 10;

// Global meshing options
Mesh.CharacteristicLengthMin = lc_near;
Mesh.CharacteristicLengthMax = lc_far;
Mesh.Algorithm3D = 1;     // Delaunay
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;

// Ensure ASCII MSH v4 output by default
Mesh.MshFileVersion = 4.1;
Mesh.Binary = 0;
