# Axis-1b Design Stub: 3D Source Ball

This document is the design stub for the axis-1b 3D source-ball
solver. Pass-11 ships only the scaffold (header, dispatch, factory
that throws); pass-12 implements the concrete subclass into the
frame defined here.

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

The host RadialLagrangianSolver (or a sibling Source3DBallSolver
class) owns the Source3DBall instance under
`cavity_geometry = THREE_DIMENSIONAL`.

## Mesh strategy

**Unstructured tetrahedra conforming to a spherical outer boundary
at the elastic radius.** The pass-12 baseline mesh:

- Cavity surface: sphere at the initial cavity radius R_c (from the
  pass-7 PHYSICS_BASED Newton solve), refined with cell size
  ~ 0.2 R_c.
- Outer surface: sphere at the elastic radius r_elastic = factor *
  R_c (factor in [3, 5] from existing pass-9/10 configuration).
- Bulk: graded-size unstructured tetrahedra, cell size growing
  smoothly from ~ 0.2 R_c near the cavity to ~ R_c at the elastic
  radius.
- Mesh generator: TetGen (BSD-licensed). The pass-12 build does not
  depend on a new external runtime library; TetGen's mesh-emitting
  CLI runs as a one-shot pre-process step and the C++ runtime
  consumes the resulting `.node` / `.ele` files.

**Resolution argument.** Pass-12 baseline ~ 50000 tetrahedra (an
elastic-radius sphere of ~ 100 R_c^3 volume divided by ~ 0.2 R_c
graded cells). For Salmon at R_c ~ 17 m and r_elastic ~ 70 m, this
is ~ 50000 cells with ~ 5 m average edge length. Computational cost
scales linearly with cell count; explicit time-stepping at CFL ~
0.4 with vp ~ 5500 m/s gives dt ~ 9e-4 s, comparable to the pass-10
1D radial CFL. Production runs in CI budget.

## 3D Drucker-Prager constitutive

**Replace the 1D radial-return formula with a 3D tensor solve.**
Pass-7 / pass-10 uses a closed-form radial-return on the deviatoric
stress with the spherical-symmetry assumption s_tt = -s_rr/2. Pass-12
generalises to:

- Per-cell stress tensor sigma_ij (6-component Voigt vector).
- Per-cell elastic strain tensor.
- Drucker-Prager yield surface: sqrt(J_2) - alpha * I_1 / 3 - k = 0
  with the existing pass-7 (alpha, k) parameters.
- Radial-return projection: explicit per-step plastic increment
  proportional to the deviatoric stress excess above yield.
- Plastic strain accumulation in the standard 6-component form.

The 3D radial return is well-known (Simo & Hughes 1998 sec 3.6,
"Computational Inelasticity"). Pass-12 uses an explicit form (no
inner Newton on the plastic multiplier) to maintain CI tractability.

## Asymmetric overburden initial state

**Depth-dependent gravity-loaded initial stress.** Pass-12:

- At t = 0, set sigma_zz(z) = -rho * g * z (negative = compressive)
  with z measured downward from the host rock free surface (or a
  configured datum).
- Set sigma_xx(z) = sigma_yy(z) = K_0 * sigma_zz with K_0 the
  configured at-rest Earth coefficient (existing pass-9 K_0 = 0.5
  default; user-configurable).
- Apply the same depth-dependent stress to every cell of the source
  ball, keyed to its centroid depth.

The asymmetric overburden creates a non-spherical free-stress
condition at the cavity wall, which seeds the asymmetric cavity
growth.

## 3D radiation block

**Design choice deferred between FEM (nodal) and cell-centred FV.**
Pass-12 must pick:

- **FEM (nodal)**: matches the rest of the host's PetscDS pointwise
  callback infrastructure (`include/numerics/PetscFEElasticity.hpp`
  pattern). Nodal radiation field couples to the nodal displacement
  via Petsc's auxiliary field mechanism. Implementation cost: high
  (build new PetscFE callbacks); long-term maintenance: low (uniform
  pattern).
- **Cell-centred FV**: matches the pass-10
  MultigroupRadiationDiffusion (cell-centred state, harmonic-mean
  face diffusion). Implementation cost: medium (port pass-10
  multigroup to a 3D unstructured cell layout); long-term maintenance:
  medium (new infrastructure).

**Recommendation**: cell-centred FV for pass-12. Reasons:
1. The pass-10 multigroup solver is the calibrated source of truth;
   reusing its kernel maintains byte-identicalness for the radiation
   physics.
2. The 3D mesh is unstructured tetrahedral; cell-centred FV needs
   only the tet adjacency and the per-face area/normal, which TetGen
   emits directly.
3. FEM-nodal radiation requires a Petsc DM with a new auxiliary
   field; the build-out has more moving parts.

Pass-12 implementation should preserve the option to switch (the
`Source3DBallConfig::radiation_discretization` enum exists for this).

## Handoff to far-field FEM (axis-1d)

**Surface-integral moment-tensor extraction at the elastic radius.**
The pass-7 / pass-10 1D solver computes M_ij = -A * sigma_rr * n_i n_j
integrated over a fixed Eulerian sphere at the elastic radius
(Day & McLaughlin 1991). The 3D solver:

- Identifies the set of surface-integral facets: the tet faces whose
  centroids fall on the spherical surface at the elastic radius.
- Sums M_ij = sum_face (n_i * sigma_jk * n_k * A_face) over those
  facets. Asymmetric content of M survives the integration when the
  cavity is asymmetric, which is the whole point of axis-1b.

The far-field FEM coupling itself is axis-1d and is independent of
axis-1b: pass-12 produces a non-isotropic moment-rate tensor; pass-X
(axis-1d) consumes that tensor in a 3D coupled-far-field-FEM solve.

## Validation strategy

Pass-12 ships with the following gates (analogous to the pass-10
multigroup gates):

- **Source3DBall.SymmetricEqualsRadial**: when the source ball is
  driven with no asymmetric overburden, the moment tensor is
  isotropic and matches the 1D radial solver result within 5%.
  This is the byte-identical regression for the new path.
- **Source3DBall.SalmonCavityRadius_3D**: Salmon 1964 cavity radius
  closes from the pass-10 factor-3 envelope to factor 1.5 (or
  better) under the 3D solver. This is the headline validation.
- **Source3DBall.ChaganCavityRadius_3D**: Chagan 1965 cavity radius
  closes from factor 5 to factor 2.
- **Source3DBall.MomentTensorCLVDContent**: under asymmetric
  overburden, the CLVD content of the moment tensor is non-zero and
  consistent with Day & McLaughlin 1991 sec 4.

## What pass-12 does NOT do

- 3D far-field FEM coupling (axis-1d).
- Topography (axis-2).
- Layered-medium fidelity (axis-3).
- ANEOS Tillotson refit (axis-4).

These are tracked separately in `docs/AXIS_1A_FIDELITY_REPORT.md`.

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
