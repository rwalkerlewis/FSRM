// ============================================================================
// Divider (1992) - US Underground Nuclear Test Style
// Complex Geological Model + Detailed Gmsh Mesh (FSRM example)
// ============================================================================
// This model is intentionally *not simplified*.
// It includes:
//   - Many stratigraphic units (9 background layers)
//   - Primary “divider” fault surface + explicit fault-core gouge volume
//   - Dipping secondary fault surface + explicit fault-core corridor volume
//   - Intrusive/salt-like ellipsoidal dome (cuts multiple units)
//   - Channel/lens body within tuffs (facies heterogeneity)
//   - Chimney/breccia pipe volume above the source
//   - Nested near-source damage halos (crushed / fractured / damaged)
//   - Strong local mesh refinement near faults, source, halos, and heterogeneities
//
// Output can be ASCII or binary (FSRM uses PETSc's Gmsh reader):
//   gmsh -3 meshes/us_divider/divider_complex_geology.geo \
//        -format msh4 -bin 0 -o meshes/us_divider/divider_complex_geology_3d_ascii.msh
//   gmsh -3 meshes/us_divider/divider_complex_geology.geo \
//        -format msh4 -bin 1 -o meshes/us_divider/divider_complex_geology_3d_bin.msh
// ============================================================================

SetFactory("OpenCASCADE");

// ----------------------------
// Domain + stratigraphy (meters)
// ----------------------------
Lx = 12000;
Ly =  8000;
zTop  =    0;
zSoil =  -50;
zA    = -200;
zT1   = -550;
zT2   = -900;
zWT   = -1500;
zV    = -2050;
zC    = -2700;
zM    = -3300;
zB    = -4200;

// Divider “fault” is represented as an explicit gouge/core *volume* (a barrier),
// not as an embedded internal fault surface. This avoids PLC self-intersection
// failures in some gmsh builds while still creating a strong hydraulic/mechanical divider.
xDiv = 6000;

// Secondary fault corridor (also represented as a volume barrier)
xF2_top = 8900;
xF2_bot = 9650;

// Mesh characteristic lengths
// (Keep geology complex; keep mesh size practical for a repo-shipped example.)
lc_far   = 500;
lc_mid   = 250;
lc_near  = 120;
lc_ultra =  80;

// Background stratigraphic volumes
Box(1) = {0, 0, zSoil, Lx, Ly, zTop  - zSoil};  // soil/colluvium
Box(2) = {0, 0, zA,    Lx, Ly, zSoil - zA};     // alluvium
Box(3) = {0, 0, zT1,   Lx, Ly, zA    - zT1};    // tuff-1
Box(4) = {0, 0, zT2,   Lx, Ly, zT1   - zT2};    // tuff-2
Box(5) = {0, 0, zWT,   Lx, Ly, zT2   - zWT};    // welded tuff
Box(6) = {0, 0, zV,    Lx, Ly, zWT   - zV};     // vitric tuff
Box(7) = {0, 0, zC,    Lx, Ly, zV    - zC};     // carbonate
Box(8) = {0, 0, zM,    Lx, Ly, zC    - zM};     // metasediments/argillite
Box(9) = {0, 0, zB,    Lx, Ly, zM    - zB};     // basement

// ----------------------------
// Heterogeneities
// ----------------------------
// Intrusive / salt-like dome (ellipsoid)
dome_x = 8500;
dome_y = 4200;
dome_z = -1800;
dome_r = 700;
Sphere(50) = {dome_x, dome_y, dome_z, dome_r};

// Channel/lens body within tuffs (ellipsoid)
lens_x = 7600;
lens_y = 2400;
lens_z = -700;
lens_rx = 1200;
lens_ry = 450;
lens_rz = 220;
// Note: This build’s .geo parser rejects Dilate/Scale transforms; use a prismatic
// lens body (still a strong heterogeneity + meshing target).
Box(60) = {lens_x - lens_rx, lens_y - lens_ry, lens_z - lens_rz, 2*lens_rx, 2*lens_ry, 2*lens_rz};

// Chimney / breccia pipe above source
pipe_x = 5650;
pipe_y = 3950;
chimney_r = 180;
// Avoid exact coincidence with stratigraphic interfaces (helps PLC robustness)
pipe_z0 = zWT + 5;
pipe_z1 = zSoil - 5;
// Use a prismatic chimney to avoid cylinder PLC corner-cases in some gmsh builds.
Box(70) = {pipe_x - chimney_r, pipe_y - chimney_r, pipe_z0, 2*chimney_r, 2*chimney_r, pipe_z1 - pipe_z0};

// Source point (used for refinement and for FSRM [EXPLOSION_SOURCE] location)
GZx = 5600;
GZy = 4000;
GZz = -800;

// Nested damage halos
R_crushed   = 260;
R_fractured = 650;
R_damaged   = 1300;
Sphere(80) = {GZx, GZy, GZz, R_crushed};
// Create *disjoint* shells to avoid nested-overlap PLC issues in boolean fragmenting:
// - crushed core: volume 80
// - fractured shell: 81 \ 80
// - damaged shell:   82 \ 83  (83 is a duplicate fractured sphere)
Sphere(81) = {GZx, GZy, GZz, R_fractured};
Sphere(82) = {GZx, GZy, GZz, R_damaged};
v_frac_shell() = BooleanDifference{ Volume{81}; Delete; }{ Volume{80}; };
Sphere(83) = {GZx, GZy, GZz, R_fractured};
v_dmg_shell() = BooleanDifference{ Volume{82}; Delete; }{ Volume{83}; Delete; };

// ----------------------------
// Explicit fault-core gouge volumes
// ----------------------------
bb_eps = 0.5;
core_t = 80;
xDivL = xDiv - core_t/2;
xDivR = xDiv + core_t/2;
Box(90) = {xDiv - core_t/2, 0, zB, core_t, Ly, zTop - zB};

core2_t = 120;
// Note: for this specific geometry we choose xF2_top < xF2_bot, so min/max are stable.
xF2_min = xF2_top;
xF2_max = xF2_bot;
Box(91) = {xF2_min - core2_t/2, 0, zB, (xF2_max - xF2_min) + core2_t, Ly, zTop - zB};

// ----------------------------
// Boolean fragmentation: make everything conforming
// ----------------------------
// Fragment stratigraphy by all internal heterogeneities and barrier corridors.
v_all() = BooleanFragments{ Volume{1:9}; Delete; }{ Volume{50,60,70,90,91,80, v_frac_shell(), v_dmg_shell()}; Delete; };

// ----------------------------
// Robust physical groups (disjoint)
// ----------------------------
// We collect special volumes first and subtract them from background layer selections.

// Special volume selections
v_dome()  = Volume In BoundingBox{dome_x - 1.35*dome_r, dome_y - 1.10*dome_r, dome_z - 0.95*dome_r,
                                 dome_x + 1.35*dome_r, dome_y + 1.10*dome_r, dome_z + 0.95*dome_r};

v_lens()  = Volume In BoundingBox{lens_x - 1.05*lens_rx, lens_y - 1.05*lens_ry, lens_z - 1.05*lens_rz,
                                 lens_x + 1.05*lens_rx, lens_y + 1.05*lens_ry, lens_z + 1.05*lens_rz};

v_pipe()  = Volume In BoundingBox{pipe_x - 1.10*chimney_r, pipe_y - 1.10*chimney_r, pipe_z0 - 10,
                                 pipe_x + 1.10*chimney_r, pipe_y + 1.10*chimney_r, pipe_z1 + 10};

v_crush_all() = Volume In BoundingBox{GZx - R_crushed,   GZy - R_crushed,   GZz - R_crushed,
                                     GZx + R_crushed,   GZy + R_crushed,   GZz + R_crushed};

v_frac_all()  = Volume In BoundingBox{GZx - R_fractured, GZy - R_fractured, GZz - R_fractured,
                                     GZx + R_fractured, GZy + R_fractured, GZz + R_fractured};

v_dmg_all()   = Volume In BoundingBox{GZx - R_damaged,   GZy - R_damaged,   GZz - R_damaged,
                                     GZx + R_damaged,   GZy + R_damaged,   GZz + R_damaged};

v_divcore() = Volume In BoundingBox{xDiv - core_t/2 - bb_eps, 0, zB, xDiv + core_t/2 + bb_eps, Ly, zTop};

v_f2core()  = Volume In BoundingBox{xF2_min - core2_t/2 - bb_eps, 0, zB, xF2_max + core2_t/2 + bb_eps, Ly, zTop};

// Make damage halos disjoint (ring-style)
v_crush() = v_crush_all();

v_frac() = v_frac_all();
v_frac() -= v_crush();

v_dmg() = v_dmg_all();
v_dmg() -= v_frac_all();

// Keep chimney and fault cores separate from damage halos
v_crush() -= v_pipe();
v_crush() -= v_divcore();
v_crush() -= v_f2core();

v_frac()  -= v_pipe();
v_frac()  -= v_divcore();
v_frac()  -= v_f2core();

v_dmg()   -= v_pipe();
v_dmg()   -= v_divcore();
v_dmg()   -= v_f2core();

// Upstream (x < xDiv) background volumes
v_soil_u()  = Volume In BoundingBox{0, 0, zSoil, xDivL - bb_eps, Ly, zTop};
v_alluv_u() = Volume In BoundingBox{0, 0, zA,    xDivL - bb_eps, Ly, zSoil};
v_tuff1_u() = Volume In BoundingBox{0, 0, zT1,   xDivL - bb_eps, Ly, zA};
v_tuff2_u() = Volume In BoundingBox{0, 0, zT2,   xDivL - bb_eps, Ly, zT1};
v_weld_u()  = Volume In BoundingBox{0, 0, zWT,   xDivL - bb_eps, Ly, zT2};
v_vit_u()   = Volume In BoundingBox{0, 0, zV,    xDivL - bb_eps, Ly, zWT};
v_carb_u()  = Volume In BoundingBox{0, 0, zC,    xDivL - bb_eps, Ly, zV};
v_meta_u()  = Volume In BoundingBox{0, 0, zM,    xDivL - bb_eps, Ly, zC};
v_base_u()  = Volume In BoundingBox{0, 0, zB,    xDivL - bb_eps, Ly, zM};

v_soil_u()  -= v_dome(); v_soil_u()  -= v_lens(); v_soil_u()  -= v_pipe(); v_soil_u()  -= v_crush_all(); v_soil_u()  -= v_frac_all(); v_soil_u()  -= v_dmg_all(); v_soil_u()  -= v_divcore(); v_soil_u()  -= v_f2core();
v_alluv_u() -= v_dome(); v_alluv_u() -= v_lens(); v_alluv_u() -= v_pipe(); v_alluv_u() -= v_crush_all(); v_alluv_u() -= v_frac_all(); v_alluv_u() -= v_dmg_all(); v_alluv_u() -= v_divcore(); v_alluv_u() -= v_f2core();
v_tuff1_u() -= v_dome(); v_tuff1_u() -= v_lens(); v_tuff1_u() -= v_pipe(); v_tuff1_u() -= v_crush_all(); v_tuff1_u() -= v_frac_all(); v_tuff1_u() -= v_dmg_all(); v_tuff1_u() -= v_divcore(); v_tuff1_u() -= v_f2core();
v_tuff2_u() -= v_dome(); v_tuff2_u() -= v_lens(); v_tuff2_u() -= v_pipe(); v_tuff2_u() -= v_crush_all(); v_tuff2_u() -= v_frac_all(); v_tuff2_u() -= v_dmg_all(); v_tuff2_u() -= v_divcore(); v_tuff2_u() -= v_f2core();
v_weld_u()  -= v_dome(); v_weld_u()  -= v_lens(); v_weld_u()  -= v_pipe(); v_weld_u()  -= v_crush_all(); v_weld_u()  -= v_frac_all(); v_weld_u()  -= v_dmg_all(); v_weld_u()  -= v_divcore(); v_weld_u()  -= v_f2core();
v_vit_u()   -= v_dome(); v_vit_u()   -= v_lens(); v_vit_u()   -= v_pipe(); v_vit_u()   -= v_crush_all(); v_vit_u()   -= v_frac_all(); v_vit_u()   -= v_dmg_all(); v_vit_u()   -= v_divcore(); v_vit_u()   -= v_f2core();
v_carb_u()  -= v_dome(); v_carb_u()  -= v_lens(); v_carb_u()  -= v_pipe(); v_carb_u()  -= v_crush_all(); v_carb_u()  -= v_frac_all(); v_carb_u()  -= v_dmg_all(); v_carb_u()  -= v_divcore(); v_carb_u()  -= v_f2core();
v_meta_u()  -= v_dome(); v_meta_u()  -= v_lens(); v_meta_u()  -= v_pipe(); v_meta_u()  -= v_crush_all(); v_meta_u()  -= v_frac_all(); v_meta_u()  -= v_dmg_all(); v_meta_u()  -= v_divcore(); v_meta_u()  -= v_f2core();
v_base_u()  -= v_dome(); v_base_u()  -= v_lens(); v_base_u()  -= v_pipe(); v_base_u()  -= v_crush_all(); v_base_u()  -= v_frac_all(); v_base_u()  -= v_dmg_all(); v_base_u()  -= v_divcore(); v_base_u()  -= v_f2core();

// Downstream (x > xDiv) background volumes
v_soil_d()  = Volume In BoundingBox{xDivR + bb_eps, 0, zSoil, Lx, Ly, zTop};
v_alluv_d() = Volume In BoundingBox{xDivR + bb_eps, 0, zA,    Lx, Ly, zSoil};
v_tuff1_d() = Volume In BoundingBox{xDivR + bb_eps, 0, zT1,   Lx, Ly, zA};
v_tuff2_d() = Volume In BoundingBox{xDivR + bb_eps, 0, zT2,   Lx, Ly, zT1};
v_weld_d()  = Volume In BoundingBox{xDivR + bb_eps, 0, zWT,   Lx, Ly, zT2};
v_vit_d()   = Volume In BoundingBox{xDivR + bb_eps, 0, zV,    Lx, Ly, zWT};
v_carb_d()  = Volume In BoundingBox{xDivR + bb_eps, 0, zC,    Lx, Ly, zV};
v_meta_d()  = Volume In BoundingBox{xDivR + bb_eps, 0, zM,    Lx, Ly, zC};
v_base_d()  = Volume In BoundingBox{xDivR + bb_eps, 0, zB,    Lx, Ly, zM};

v_soil_d()  -= v_dome(); v_soil_d()  -= v_lens(); v_soil_d()  -= v_pipe(); v_soil_d()  -= v_crush_all(); v_soil_d()  -= v_frac_all(); v_soil_d()  -= v_dmg_all(); v_soil_d()  -= v_divcore(); v_soil_d()  -= v_f2core();
v_alluv_d() -= v_dome(); v_alluv_d() -= v_lens(); v_alluv_d() -= v_pipe(); v_alluv_d() -= v_crush_all(); v_alluv_d() -= v_frac_all(); v_alluv_d() -= v_dmg_all(); v_alluv_d() -= v_divcore(); v_alluv_d() -= v_f2core();
v_tuff1_d() -= v_dome(); v_tuff1_d() -= v_lens(); v_tuff1_d() -= v_pipe(); v_tuff1_d() -= v_crush_all(); v_tuff1_d() -= v_frac_all(); v_tuff1_d() -= v_dmg_all(); v_tuff1_d() -= v_divcore(); v_tuff1_d() -= v_f2core();
v_tuff2_d() -= v_dome(); v_tuff2_d() -= v_lens(); v_tuff2_d() -= v_pipe(); v_tuff2_d() -= v_crush_all(); v_tuff2_d() -= v_frac_all(); v_tuff2_d() -= v_dmg_all(); v_tuff2_d() -= v_divcore(); v_tuff2_d() -= v_f2core();
v_weld_d()  -= v_dome(); v_weld_d()  -= v_lens(); v_weld_d()  -= v_pipe(); v_weld_d()  -= v_crush_all(); v_weld_d()  -= v_frac_all(); v_weld_d()  -= v_dmg_all(); v_weld_d()  -= v_divcore(); v_weld_d()  -= v_f2core();
v_vit_d()   -= v_dome(); v_vit_d()   -= v_lens(); v_vit_d()   -= v_pipe(); v_vit_d()   -= v_crush_all(); v_vit_d()   -= v_frac_all(); v_vit_d()   -= v_dmg_all(); v_vit_d()   -= v_divcore(); v_vit_d()   -= v_f2core();
v_carb_d()  -= v_dome(); v_carb_d()  -= v_lens(); v_carb_d()  -= v_pipe(); v_carb_d()  -= v_crush_all(); v_carb_d()  -= v_frac_all(); v_carb_d()  -= v_dmg_all(); v_carb_d()  -= v_divcore(); v_carb_d()  -= v_f2core();
v_meta_d()  -= v_dome(); v_meta_d()  -= v_lens(); v_meta_d()  -= v_pipe(); v_meta_d()  -= v_crush_all(); v_meta_d()  -= v_frac_all(); v_meta_d()  -= v_dmg_all(); v_meta_d()  -= v_divcore(); v_meta_d()  -= v_f2core();
v_base_d()  -= v_dome(); v_base_d()  -= v_lens(); v_base_d()  -= v_pipe(); v_base_d()  -= v_crush_all(); v_base_d()  -= v_frac_all(); v_base_d()  -= v_dmg_all(); v_base_d()  -= v_divcore(); v_base_d()  -= v_f2core();

// Physical volumes: background (upstream/downstream)
Physical Volume("soil_upstream") = {v_soil_u()};
Physical Volume("alluvium_upstream") = {v_alluv_u()};
Physical Volume("tuff1_upstream") = {v_tuff1_u()};
Physical Volume("tuff2_upstream") = {v_tuff2_u()};
Physical Volume("welded_tuff_upstream") = {v_weld_u()};
Physical Volume("vitric_tuff_upstream") = {v_vit_u()};
Physical Volume("carbonate_upstream") = {v_carb_u()};
Physical Volume("metasediment_upstream") = {v_meta_u()};
Physical Volume("basement_upstream") = {v_base_u()};

Physical Volume("soil_downstream") = {v_soil_d()};
Physical Volume("alluvium_downstream") = {v_alluv_d()};
Physical Volume("tuff1_downstream") = {v_tuff1_d()};
Physical Volume("tuff2_downstream") = {v_tuff2_d()};
Physical Volume("welded_tuff_downstream") = {v_weld_d()};
Physical Volume("vitric_tuff_downstream") = {v_vit_d()};
Physical Volume("carbonate_downstream") = {v_carb_d()};
Physical Volume("metasediment_downstream") = {v_meta_d()};
Physical Volume("basement_downstream") = {v_base_d()};

// Physical volumes: special domains (disjoint)
Physical Volume("intrusive_dome") = {v_dome()};
Physical Volume("tuff_channel_lens") = {v_lens()};
Physical Volume("chimney_pipe") = {v_pipe()};
Physical Volume("damage_crushed") = {v_crush()};
Physical Volume("damage_fractured") = {v_frac()};
Physical Volume("damage_damaged") = {v_dmg()};
Physical Volume("divider_fault_core") = {v_divcore()};
Physical Volume("secondary_fault_core") = {v_f2core()};

// ----------------------------
// Physical surfaces: boundaries
// ----------------------------
Physical Surface("surface")  = Surface In BoundingBox{0, 0, zTop - bb_eps, Lx, Ly, zTop + bb_eps};
Physical Surface("bottom")   = Surface In BoundingBox{0, 0, zB   - bb_eps, Lx, Ly, zB   + bb_eps};
Physical Surface("upstream") = Surface In BoundingBox{0 - bb_eps, 0, zB, 0 + bb_eps, Ly, zTop};
Physical Surface("downstream") = Surface In BoundingBox{Lx - bb_eps, 0, zB, Lx + bb_eps, Ly, zTop};
Physical Surface("south")    = Surface In BoundingBox{0, 0 - bb_eps, zB, Lx, 0 + bb_eps, zTop};
Physical Surface("north")    = Surface In BoundingBox{0, Ly - bb_eps, zB, Lx, Ly + bb_eps, zTop};

// ----------------------------
// Mesh sizing fields
// ----------------------------
// Explosion-centric refinement
Field[1] = Ball;
Field[1].XCenter = GZx;
Field[1].YCenter = GZy;
Field[1].ZCenter = GZz;
Field[1].Radius  = R_damaged + 500;
Field[1].VIn     = lc_near;
Field[1].VOut    = lc_far;

Field[2] = Ball;
Field[2].XCenter = GZx;
Field[2].YCenter = GZy;
Field[2].ZCenter = GZz;
Field[2].Radius  = R_fractured + 250;
Field[2].VIn     = lc_ultra;
Field[2].VOut    = lc_near;

// Divider corridor refinement (ball around divider fault-core volume)
Field[3] = Ball;
Field[3].XCenter = xDiv;
Field[3].YCenter = Ly/2;
Field[3].ZCenter = (zTop + zB)/2;
Field[3].Radius  = 1200;
Field[3].VIn     = lc_ultra;
Field[3].VOut    = lc_mid;

// Secondary corridor refinement (ball around secondary fault-core corridor)
Field[4] = Ball;
Field[4].XCenter = (xF2_top + xF2_bot)/2;
Field[4].YCenter = Ly/2;
Field[4].ZCenter = (zTop + zB)/2;
Field[4].Radius  = 1400;
Field[4].VIn     = lc_ultra;
Field[4].VOut    = lc_mid;

// Intrusive dome refinement
Field[7] = Ball;
Field[7].XCenter = dome_x;
Field[7].YCenter = dome_y;
Field[7].ZCenter = dome_z;
Field[7].Radius  = 1700;
Field[7].VIn     = lc_mid;
Field[7].VOut    = lc_far;

// Lens refinement
Field[8] = Ball;
Field[8].XCenter = lens_x;
Field[8].YCenter = lens_y;
Field[8].ZCenter = lens_z;
Field[8].Radius  = 1600;
Field[8].VIn     = lc_near;
Field[8].VOut    = lc_far;

// Chimney refinement
Field[9] = Ball;
Field[9].XCenter = pipe_x;
Field[9].YCenter = pipe_y;
Field[9].ZCenter = (pipe_z0 + pipe_z1)/2;
Field[9].Radius  = 750;
Field[9].VIn     = lc_ultra;
Field[9].VOut    = lc_mid;

Field[20] = Min;
Field[20].FieldsList = {1,2,3,4,7,8,9};
Background Field = 20;

// Global meshing options
Mesh.CharacteristicLengthMin = lc_ultra;
Mesh.CharacteristicLengthMax = lc_far;
// Use classic Delaunay tetrahedralization (more forgiving for complex PLCs here).
Mesh.Algorithm3D = 1;
Mesh.Optimize = 1;
// Avoid “overlapping facets” false-positives for nearly-coplanar fragment faces.
Mesh.AngleToleranceFacetOverlap = 0.05;

// Output MSH v4.1 (ASCII or binary; set Mesh.Binary accordingly)
Mesh.MshFileVersion = 4.1;
Mesh.Binary = 0;
