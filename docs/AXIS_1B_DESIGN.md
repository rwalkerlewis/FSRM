# Axis-1b Design: 3D Source Ball

This document is the canonical design for the axis-1b 3D source-ball
solver. Pass-11 shipped the scaffold (header, dispatch, factory that
threw on construction). Pass-13a (foundation slice) ships the leaf
data structures: TetGen mesh generation tool, DMPlex mesh load and
distribute, and a `Source3DBallImpl` class whose `step()` still throws
but whose `initialize()` is real. Pass-13b implements the 3D
constitutive update, asymmetric overburden initial state, and grey
radiation diffusion. Pass-13c implements the surface-integral
moment-tensor extraction, host-side delegation from
`RadialLagrangianSolver`, end-to-end MPI=4 Salmon validation, and
the spherical-symmetry regression-equivalence headline gate.

The original pass-12 plan was to land the entire 3D solver in one
pass. Pass-12 became housekeeping in practice (see
`docs/HISTORIC_NUCLEAR_FIDELITY.md` 4k); pass-13 is therefore split
into three foundation / physics / validation slices.

## Motivation

Pass-7 through pass-10 axis-1a uses a 1D radial Lagrangian
elastoplastic shock solver in spherical symmetry. Two observed
residuals are limited by spherical symmetry:

- **Salmon CavityRadius** (factor 3 envelope vs spec target 5%).
- **Chagan CavityRadius** (factor 5 envelope vs spec target 3).

The actual underground cavity is asymmetric under overburden: the
upward-facing hemisphere experiences less confining pressure than
the downward hemisphere, so the cavity grows preferentially upward
(Patton 1991 LANL technical report). This asymmetry generates CLVD
moment-tensor content that the spherically-symmetric 1D solver
cannot represent (Stevens & Day 1985).

Closing the cavity-radius gates to spec requires a 3D source-ball
solver that resolves:

1. Asymmetric overburden initial stress.
2. 3D Drucker-Prager constitutive on a 3D unstructured mesh.
3. 3D radiation block (matching the multigroup tabulated-opacity
   path of pass-10).
4. Surface-integral moment-tensor extraction at the elastic radius
   for handoff to the far-field FEM (axis-1d).

## Pass-13 slicing

| Slice | Deliverable | Status |
|-------|-------------|--------|
| pass-13a foundation | TetGen pre-process tool, `Source3DBallMesh` (DMPlex create + distribute + per-vertex marker label), `Source3DBallImpl` skeleton (initialize loads mesh; step throws), unit tests, factory throw replaced. | LANDED (PR #129) |
| pass-13b physics | 3D Drucker-Prager radial-return constitutive (Simo & Hughes 1998), asymmetric overburden initial state, cell-centred FV grey radiation diffusion (PETSc Mat + KSP + matter Newton outer loop), host-side `RadialLagrangianSolver -> Source3DBallImpl` delegation, ConfigReader plumbing for `[NEAR_FIELD_SOURCE]` 3D sub-keys, headline cell-level equivalence gate vs 1D scalar reduction. | LANDED (this branch) |
| pass-13c validation | Surface-integral moment-tensor extraction (Day & McLaughlin 1991), [OUTPUT] config block + wavefield + source-ball 3D HDF5/XDMF output, ParaView state files, two new academic verification anchors (Lamb's problem, layered halfspace explosion), CLVD-content-with-overburden best-effort gate, cavity aspect ratio measurement (no gate; ~1.00 in pass-13c). The end-to-end MPI=4 Salmon-with-overburden integration test running to simulation completion is deferred to a follow-up PR (the missing piece is the Salmon-specific TetGen mesh; the host delegation + surface integral + wavefield output infrastructure are in place). The literature-range CLVD ratio (0.05-0.30) and the cavity aspect ratio gate are pass-14 deliverables once 3D Lagrangian advection lands. | LANDED (this branch) |

Pass-14+ continues with multigroup-3D radiation, tabulated EOS /
opacity in 3D, RK3-SSP / BDF2 in 3D, and axis-1d 3D far-field FEM
coupling.

## Interface

`include/domain/explosion/Source3DBall.hpp` defines the contract:

```cpp
class Source3DBall {
public:
    virtual ~Source3DBall() = default;
    virtual void initialize(const Source3DBallConfig& cfg) = 0;
    virtual void step(double dt) = 0;
    virtual void getMomentTensor(std::array<double, 6>& M) const = 0;
    virtual void getMomentRateTensor(std::array<double, 6>& Mdot) const = 0;
    virtual void getState(Source3DBallState& state) const = 0;
    virtual const char* name() const = 0;
};
```

Pass-13a foundation adds:

- `cavity_radius_m`, `mesh_path`, `overburden_K0` fields to
  `Source3DBallConfig`.
- `Source3DBallImpl` final class implementing the interface (mesh
  load works; `step` and friends are pass-13b/c work and throw or
  return zeros).

Pass-13b lands the host wiring: `RadialLagrangianSolver::setConfig`
constructs a `Source3DBallImpl` when `cavity_geometry =
THREE_DIMENSIONAL` and calls `initialize()` with the full sub-config
populated from the parsed `[NEAR_FIELD_SOURCE]` keys. `step()`
forwards `dt` to the impl; `getMomentTensor` and
`getMomentRateTensor` forward as well (returning zeros from the impl
in pass-13b; surface-integral extraction is pass-13c). A new
`name()` accessor reports
`RadialLagrangianSolver+Source3DBallImpl_v1_pass13b_physics` when
the delegation is active.

## Mesh strategy

Pass-13a foundation uses **TetGen as a build-host pre-process CLI**,
not as a runtime library. The Python tool
`tools/mesh_generation/build_source_ball_mesh.py` generates a `.poly`
PSLG describing two concentric icospheres (cavity surface, elastic
radius outer surface) with a hole marker at the origin, runs
`tetgen -pq1.4Ya`, and emits the resulting `.node` / `.ele` pair.
The C++ runtime parses the plain-text TetGen output via
`Source3DBallMesh` (`include/domain/explosion/Source3DBallMesh.hpp`).

This decision keeps `fsrm-ci:local` free of TetGen, keeps the meshes
deterministic across CI runs (cached under
`cache/source_ball_meshes/`), and removes a build-time dependency.
Pass-13b will wire the runtime to consume the cached meshes.

Pass-13 baseline:

- Cavity surface: sphere at the initial cavity radius `R_c`
  (from the pass-7 PHYSICS_BASED Newton solve), refined with cell
  size `~ 0.2 R_c`.
- Outer surface: sphere at the elastic radius
  `r_elastic = factor * R_c` (factor in [3, 5] from existing pass-9
  / pass-10 configuration).
- Bulk: graded-size unstructured tetrahedra, cell size growing
  smoothly from `~ 0.2 R_c` near the cavity to `~ R_c` at the
  elastic radius.

Resolution argument: pass-13 baseline `~ 50000 tetrahedra` (an
elastic-radius sphere of `~ 100 R_c^3` volume divided by `~ 0.2 R_c`
graded cells). For Salmon at `R_c ~ 17 m` and `r_elastic ~ 70 m`,
this is `~ 50000 cells` with `~ 5 m` average edge length.
Computational cost scales linearly with cell count; explicit time-
stepping at CFL `~ 0.4` with `vp ~ 5500 m/s` gives `dt ~ 9e-4 s`,
comparable to the pass-10 1D radial CFL. Production runs in CI
budget.

## DMPlex pattern (pass-13a foundation)

`Source3DBallMesh::loadFromTetGen` follows the existing
`src/io/GmshIO.cpp` pattern in this repository:

1. Rank 0 parses the `.node` and `.ele` files.
2. Rank 0 calls `DMPlexCreateFromCellListPetsc` with the cell list
   and vertex coordinates, `interpolate=PETSC_TRUE` to build the
   full Hasse diagram. Other ranks call the same with empty arrays.
3. Per-vertex markers from the `.node` file (1 = outer-elastic,
   2 = inner-cavity, 0 = interior) are written into a DMLabel named
   `SourceBallVertexMarker` on the rank-0 DM and an empty label on
   other ranks. `DMPlexDistribute` migrates the label.
4. `DMPlexDistribute(*dm, 0, NULL, &dm_dist)` partitions across
   ranks with overlap=0.
5. The mesh exposes `numLocalCells()`, `numGlobalCells()`,
   `getDM()`, the marker label, and a few sanity-check accessors
   (`localMinVertexRadius`, `localMaxVertexRadius`,
   `numCavityMarkedVertices`, `numElasticMarkedVertices`).

The pass-12 parallel KSP fix (`bjacobi+sub_lu` for MPI > 1) does not
yet engage in pass-13a because nothing solves a system on the 3D
mesh; pass-13b's radiation diffusion will hook into the same
parallel KSP convention as the rest of the project.

## 3D Drucker-Prager constitutive (pass-13b)

Replace the 1D radial-return formula with a 3D tensor solve. Pass-7
/ pass-10 uses a closed-form radial-return on the deviatoric stress
with the spherical-symmetry assumption `s_tt = -s_rr / 2`. Pass-13b
generalises to:

- Per-cell stress tensor `sigma_ij` (6-component Voigt vector).
- Per-cell elastic strain tensor.
- Drucker-Prager yield surface:
  `sqrt(J_2(s)) - alpha * I_1(sigma) / 3 - k = 0` with the existing
  pass-7 `(alpha, k)` parameters.
- Trial elastic predictor: `sigma_trial = sigma + C : delta_eps`.
- Radial-return projection: explicit per-step plastic increment
  `delta_lambda` chosen such that `f(sigma_new) ~ 0` to first order
  in the deviatoric direction.
- Plastic strain accumulation in the standard 6-component form.

The 3D radial return is well-known (Simo & Hughes 1998 sec 3.6,
"Computational Inelasticity"). Pass-13b uses an explicit form (no
inner Newton on the plastic multiplier) to maintain CI tractability.
Pass-15+ may switch to an implicit return if accuracy demands.

## Asymmetric overburden initial state (pass-13b)

Depth-dependent gravity-loaded initial stress:

- At `t = 0`, set `sigma_zz(z) = -rho * g * z` (negative =
  compressive) with `z` measured downward from the host rock free
  surface or a configured datum.
- Set `sigma_xx(z) = sigma_yy(z) = K_0 * sigma_zz` with `K_0` the
  configured at-rest Earth coefficient (existing pass-9 default
  `K_0 = 0.5`; user-configurable via the new
  `Source3DBallConfig::overburden_K0` field, present in pass-13a).
- Apply the same depth-dependent stress to every cell of the source
  ball, keyed to its centroid depth.

The asymmetric overburden creates a non-spherical free-stress
condition at the cavity wall, which seeds the asymmetric cavity
growth.

## 3D radiation block (pass-13b)

Pass-13b ships **cell-centred FV grey-diffusion** as the only
working radiation discretization. Reasons (unchanged from pass-12
plan):

1. The pass-10 multigroup solver is the calibrated source of truth;
   reusing its kernel maintains byte-identicalness for the radiation
   physics on the spherical-symmetric reference path.
2. The 3D mesh is unstructured tetrahedral; cell-centred FV needs
   only the tet adjacency and the per-face area / normal, both of
   which DMPlex provides directly.
3. FEM-nodal radiation requires a Petsc DM with a new auxiliary
   field; the build-out has more moving parts.

`Source3DBallConfig::RadiationDiscretization::FEM_NODAL` remains a
named scaffold. Pass-13a's `Source3DBallImpl::initialize` throws a
clear "pass-15+" message when this option is selected.

Pass-14 extends to multigroup 3D (G coupled diffusion solves per
Newton iteration) and the tabulated EOS / opacity patches.

## Handoff to far-field FEM (axis-1d, pass-13c)

Pass-13c implements the surface-integral moment-tensor extraction at
the elastic radius. Following Day & McLaughlin 1991:

- Identify the set of surface-integral facets via DMLabel: the tet
  faces whose centroids lie on the elastic-radius surface.
- Sum
  `M_ij = sum_face (n_i * sigma_jk * n_k * A_face)` over those
  facets.
- Asymmetric content of `M` survives the integration when the cavity
  is asymmetric, which is the whole point of axis-1b.

The far-field FEM coupling itself is axis-1d and is independent of
axis-1b: pass-13c produces a non-isotropic moment-rate tensor;
axis-1d (a future pass) consumes that tensor in a 3D coupled-far-
field-FEM solve.

## Validation strategy (pass-13c)

Pass-13c gates (analogous to the pass-10 multigroup gates):

- **Source3DBall.SymmetricEqualsRadial** (headline): with isotropic
  IC, no overburden, and spherical mesh refinement, the 3D `M(t)`
  for Salmon matches pass-10 1D `M(t)` within factor 1.1 at all
  time samples. This is the "the 3D mesh and 3D constitutive don't
  break the physics" gate.
- **Source3DBall.AsymmetricOverburdenProducesCLVD**: with
  `K_0 = 0.5` overburden, `M(t)` has nonzero CLVD content. The
  CLVD axis aligns with vertical (within 5 degrees) and the sign is
  compressive-vertical-dominant.
- **Source3DBall.SalmonCavityRadius_3D**: Salmon 1964 cavity radius
  closes from the pass-10 factor-3 envelope to factor 1.5 (or
  better) under the 3D solver.
- **Source3DBall.MomentTensorCLVDContent**: under asymmetric
  overburden, the CLVD content of the moment tensor is non-zero and
  consistent with Day & McLaughlin 1991 sec 4.

## What pass-13 does NOT do

Pass-13 (in any of its slices):

- 3D far-field FEM coupling (axis-1d, future pass).
- Topography (axis-2, future pass).
- Layered-medium fidelity (axis-3, future pass).
- ANEOS Tillotson refit (axis-4, future pass).
- 3D multigroup radiation (pass-14).
- Tabulated EOS / opacity patches in 3D (pass-14).
- Higher-order time integration in 3D (pass-14).
- Anchors beyond Salmon (Chagan, DPRK 2017, etc.; pass-15+).

## References

- Day, S. M. and McLaughlin, K. L. (1991), "Effects of plasticity on
  the seismic source of underground explosions", J. Geophys. Res.
  96(B2), pp 1955-1975.
- Patton, H. J. (1991), "Decoupling and topological factors at NTS:
  Cavity asymmetry under overburden", LANL technical report.
- Stevens, J. L. and Day, S. M. (1985), "The physical basis of mb:Ms
  and variable frequency magnitude methods for earthquake/explosion
  discrimination", J. Geophys. Res. 90(B4), pp 3009-3020.
- Simo, J. C. and Hughes, T. J. R. (1998), "Computational
  Inelasticity", Springer, sec 3.6 (3D Drucker-Prager radial return).
- Aki, K. and Richards, P. G. (2002), "Quantitative Seismology" 2nd
  ed, ch 4 (surface-integral moment-tensor extraction).
- Si, H. (2015), "TetGen, a Delaunay-based quality tetrahedral mesh
  generator", ACM Trans. Math. Software 41(2), Article 11.
- Liu, A. and Joe, B. (1995), "Quality local refinement of
  tetrahedral meshes based on bisection", J. Sci. Comput. 16(6).
- Hoek, E. and Brown, E. T. (1980), "Underground excavations in
  rock", Institution of Mining and Metallurgy (K_0 default).
