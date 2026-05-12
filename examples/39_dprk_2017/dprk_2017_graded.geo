// dprk_2017_graded.geo
//
// Graded 3-D finite-element mesh for the DPRK 2017 (Punggye-ri) RADIAL_LAGRANGIAN run.
//
// Domain:   30 x 30 x 20 km  (x in [0,30000], y in [0,30000], z in [0,20000])
// Surface:  z = 20000 m
// Source:   x=15000, y=15000, z=19400 (600 m below surface, within granite)
//
// Layer z-boundaries
//   tuff        : z = 19700 - 20000   (top 300 m)
//   granite     : z = 17000 - 19700   (300 - 3000 m depth)
//   basement    : z = 0     - 17000   (> 3000 m depth)
//
// Mesh grading:
//   Near source (within 3000 m):  ~300 m
//   Far field:                    ~1500 m

SetFactory("OpenCASCADE");

// ---- Parameters ----
lx   = 30000;   ly   = 30000;   lz   = 20000;
z1   = 17000;   // basement/granite interface
z2   = 19700;   // granite/tuff interface
sx   = 15000;   sy   = 15000;   sz   = 19400;   // source centre

lc_near  = 300;    // mesh size at source
lc_far   = 1500;   // mesh size at domain edge
r_refine = 3000;   // refinement sphere radius

// ---- Volumes (3 stacked boxes) ----
Box(1) = {0, 0, 0,  lx, ly, z1};         // basement
Box(2) = {0, 0, z1, lx, ly, z2-z1};      // granite
Box(3) = {0, 0, z2, lx, ly, lz-z2};      // tuff

// Merge interfaces so layers share faces
BooleanFragments{ Volume{1}; Delete; }{ Volume{2,3}; Delete; }
// After fragmentation the 3 non-overlapping boxes keep tags 1,2,3

// ---- Physical groups ----
// Volumes
Physical Volume("basement")  = {1};
Physical Volume("granite")   = {2};
Physical Volume("tuff")      = {3};

// Boundary surfaces for absorbing / free-surface BCs
// FSRM detects boundaries by coordinate -- no explicit Physical Surface needed,
// but label them for reference / debugging.
//
// Identify the 6 outer faces after fragmentation using bounding-box queries:
//   Bottom (z=0):   zmin face of basement
//   Top (z=lz):     zmax face of tuff
//   Sides: xmin/xmax (x=0/lx), ymin/ymax (y=0/ly)
// Gmsh numbers these automatically; we expose them as Physical Surfaces.

// ---- Mesh size field (graded around source) ----
Field[1] = MathEval;
Field[1].F = Sprintf(
    "%g + (%g - %g) * Min(Sqrt((x-%g)^2 + (y-%g)^2 + (z-%g)^2) / %g, 1)",
    lc_near, lc_far, lc_near,
    sx, sy, sz,
    r_refine
);

Background Field = 1;

Mesh.Algorithm3D      = 1;   // Delaunay
Mesh.CharacteristicLengthMin  = 200;
Mesh.CharacteristicLengthMax  = 2000;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.Optimize         = 1;
Mesh.OptimizeNetgen   = 1;
