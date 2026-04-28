// =============================================================================
// Sedan 1962 near-field mesh
// =============================================================================
//
// Generates a graded tetrahedral mesh for the Sedan (Project Plowshare,
// 1962-07-06, 104 kt, depth of burial 194 m) near-field PPV-versus-range
// validation in Milestone M1.4 of the FSRM Nuclear Explosion Monitoring
// roadmap.
//
// Domain
// ------
//   - Box, 8 km x 8 km laterally (4 km half-width about the surface
//     ground zero) and 2 km deep.
//   - Coordinate origin at the surface directly above the shot point.
//   - Free surface at z = 0.  The shot point is at (0, 0, -194).
//
// Element grading
// ---------------
//   - 50 m target element edge length within a 500 m radius sphere
//     around the shot point.
//   - Graded out to a target size of ~300 m near the lateral and
//     bottom outer boundaries.  The roadmap M1.1 spec calls out 500 m
//     at the boundaries; the actual setting is tightened to 300 m so
//     that the longest tet edges produced by the 3D Delaunay mesher
//     stay below the 750 m sanity bound asserted in
//     tests/integration/test_sedan_1962_mesh.cpp.  This is a
//     conservative choice -- it costs about a factor of 2 in element
//     count but does not change the high-resolution sphere or the
//     frequency band targeted in M1.4.
//   - Tetrahedral cells (algorithm 1 / Delaunay 3D).
//
// Frequency-band justification
// ----------------------------
//   Target faithful representation up to 5 Hz.  The slowest wave in
//   the model (alluvium S-wave at Vs = 700 m/s, see models/sedan_1962.vel)
//   has a wavelength of 140 m at 5 Hz.  With 50 m elements that is 2.8
//   nodes per wavelength which is below the standard 6-10 nodes-per-
//   wavelength rule for second-order finite elements.  Therefore the
//   high-resolution sphere is sized to keep alluvium S-waves accurate
//   to roughly 3 Hz, and the M1.4 validation is restricted to the band
//   relevant to peak particle velocity (PPV is dominated by the source
//   corner frequency, which for Sedan is approximately 1 - 2 Hz).
//   If the validation in M1.4 fails because of insufficient resolution,
//   drop the high-resolution element size to 25 m as a remediation
//   step before changing source physics.
//
// Physical groups (Gmsh -> FSRM mapping)
// --------------------------------------
//   Volume "rock"   (id 1) -> single material region.
//   Surface "zmax"  (id 16) -> free surface (top, z = 0).
//   Surfaces "xmin"/"xmax"/"ymin"/"ymax"/"zmin" -> absorbing-BC faces.
//   The xmin/xmax/ymin/ymax/zmin/zmax naming is the convention used by
//   tests/integration/test_gmsh_import.cpp and meshes/test_box_gmsh.geo.
//
// Build:
//   bash scripts/meshing/build_sedan_1962_mesh.sh
//   (writes meshes/historical/sedan_1962_nearfield.msh in MSH 2.2 format)
// =============================================================================

SetFactory("OpenCASCADE");

// ----- Domain extents (metres) ------------------------------------------------
Lx = 8000.0;    // total x extent
Ly = 8000.0;    // total y extent
Lz = 2000.0;    // depth (z extent)

// Source at surface ground zero, depth of burial 194 m.
src_x = 0.0;
src_y = 0.0;
src_z = -194.0;

// High-resolution sphere radius around the source.
r_fine = 500.0;

// Target element sizes.
lc_fine   =  50.0;   // inside r_fine of the source
lc_coarse = 300.0;   // at the outer boundaries

// ----- Geometry ---------------------------------------------------------------
// Box() takes the lower-left corner and the side lengths.
Box(1) = {-Lx/2, -Ly/2, -Lz, Lx, Ly, Lz};

Physical Volume("rock", 1) = {1};

// Boundary surfaces by bounding-box selection (robust against face renumbering).
eps = 1e-6;
xmin_s[] = Surface In BoundingBox{-Lx/2 - eps, -Ly/2 - eps, -Lz - eps,
                                  -Lx/2 + eps,  Ly/2 + eps,        eps};
xmax_s[] = Surface In BoundingBox{ Lx/2 - eps, -Ly/2 - eps, -Lz - eps,
                                   Lx/2 + eps,  Ly/2 + eps,        eps};
ymin_s[] = Surface In BoundingBox{-Lx/2 - eps, -Ly/2 - eps, -Lz - eps,
                                   Lx/2 + eps, -Ly/2 + eps,        eps};
ymax_s[] = Surface In BoundingBox{-Lx/2 - eps,  Ly/2 - eps, -Lz - eps,
                                   Lx/2 + eps,  Ly/2 + eps,        eps};
zmin_s[] = Surface In BoundingBox{-Lx/2 - eps, -Ly/2 - eps, -Lz - eps,
                                   Lx/2 + eps,  Ly/2 + eps, -Lz + eps};
zmax_s[] = Surface In BoundingBox{-Lx/2 - eps, -Ly/2 - eps,       -eps,
                                   Lx/2 + eps,  Ly/2 + eps,        eps};

Physical Surface("xmin", 11) = {xmin_s[]};
Physical Surface("xmax", 12) = {xmax_s[]};
Physical Surface("ymin", 13) = {ymin_s[]};
Physical Surface("ymax", 14) = {ymax_s[]};
Physical Surface("zmin", 15) = {zmin_s[]};
Physical Surface("zmax", 16) = {zmax_s[]};   // free surface

// ----- Size field: graded refinement around the shot point --------------------
// A single auxiliary point at the source location is the distance origin.
// The Threshold field then maps distance -> element size, linearly from
// lc_fine inside r_fine to lc_coarse at distance ~ 3 * r_fine.  Setting the
// transition over ~1.5 km keeps the grading wavelength longer than the local
// element size.
src_pt = newp;
Point(src_pt) = {src_x, src_y, src_z, lc_fine};
// Embed the source point inside the volume so the mesher anchors a node
// at the shot location and respects the local fine size there.
Point{src_pt} In Volume{1};

Field[1] = Distance;
Field[1].PointsList = {src_pt};

Field[2] = Threshold;
Field[2].InField  = 1;
Field[2].SizeMin  = lc_fine;
Field[2].SizeMax  = lc_coarse;
Field[2].DistMin  = r_fine;
Field[2].DistMax  = 3.0 * r_fine;
Field[2].StopAtDistMax = 0;

Background Field = 2;

// Disable the default mesh-size-from-points / extend-from-boundary behaviour
// so the size field controls the mesh density everywhere.
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.MeshSizeExtendFromBoundary = 0;

Mesh.CharacteristicLengthMin = lc_fine;
Mesh.CharacteristicLengthMax = lc_coarse;

// 3D Delaunay produces high-quality tets.
Mesh.Algorithm   = 5;     // Delaunay (2D)
Mesh.Algorithm3D = 1;     // Delaunay (3D)

// Linear (first-order) tetrahedra to match the rest of the FSRM
// elastodynamic test suite (see tests/integration/test_explosion_seismogram.cpp).
Mesh.ElementOrder = 1;
