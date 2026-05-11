/**
 * @file Source3DBall.hpp
 * @brief Pass-11 (axis 1b) scaffolded interface for the 3D source-ball
 *        solver. Pass-11 ships only the header, type definitions, and
 *        dispatch wiring; the concrete implementation is named pass-12
 *        work.
 *
 * Constraint and motivation
 * -------------------------
 * The pass-7-through-pass-10 axis-1a path uses a 1D radial Lagrangian
 * elastoplastic shock solver that assumes spherical symmetry. The
 * Salmon and Chagan cavity-radius gates remain at factor 3-5 envelope
 * (vs spec target of 5% relative) because the actual underground
 * cavity is asymmetric under overburden: the upward-facing
 * hemisphere experiences less confining pressure than the downward
 * hemisphere, so the cavity grows preferentially upward.
 *
 * Axis-1b removes the spherical-symmetry assumption by replacing the
 * 1D radial Lagrangian with a 3D unstructured-tetrahedral source-ball
 * solver. The 3D mesh conforms to a spherical outer boundary at the
 * elastic radius, where the surface-integral moment-tensor extraction
 * hands off to the far-field FEM (axis-1d separately resolves the
 * 3D far-field coupling).
 *
 * Pass-11 ships:
 *  - cavity_geometry config knob in RadialLagrangianSolver::Config
 *    selecting SPHERICAL (pass-10 default) or THREE_DIMENSIONAL
 *    (scaffold; throws on selection).
 *  - Source3DBall abstract interface (this file).
 *  - Throw-on-selection dispatch in RadialLagrangianSolver::setConfig
 *    so callers cannot silently opt into a non-existent code path.
 *  - docs/AXIS_1B_DESIGN.md design stub describing the prospective
 *    3D mesh strategy, 3D Drucker-Prager extension, asymmetric
 *    overburden BC, and 3D radiation diffusion choice.
 *
 * Pass-12 will replace the throw with a concrete Source3DBallImpl
 * subclass that carries this interface to a working unstructured
 * tetrahedral solver.
 *
 * References
 * ----------
 *  - Day & McLaughlin (1991), "Effects of plasticity on the seismic
 *    source of underground explosions", J. Geophys. Res. 96, pp 1955-1975.
 *  - Patton, H. J. (1991), "Decoupling and topological factors at NTS:
 *    Cavity asymmetry under overburden", LANL technical report.
 *  - Stevens, J. L. and Day, S. M. (1985), "The physical basis of mb:Ms
 *    and variable frequency magnitude methods for earthquake/explosion
 *    discrimination", J. Geophys. Res. 90, pp 3009-3020 (3D source
 *    region effects on radiated wavefield).
 *  - Aki, K. and Richards, P. G. (2002), "Quantitative Seismology"
 *    2nd ed, ch 4 (surface-integral moment-tensor extraction).
 */

#ifndef NEAR_FIELD_SOURCE_3D_BALL_HPP
#define NEAR_FIELD_SOURCE_3D_BALL_HPP

#include <array>
#include <memory>
#include <string>

namespace FSRM {

/// Configuration sub-block for the axis-1b 3D source ball.
/// Pass-11 scaffold: factory threw on construction.
/// Pass-13a foundation: factory returns a real Source3DBallImpl whose
/// initialize() loads a TetGen-generated mesh from disk via
/// Source3DBallMesh and constructs a distributed DMPlex. step() and
/// getMomentTensor() remain pass-13b/c work and throw with a clear
/// message naming the follow-on pass.
struct Source3DBallConfig
{
    /// Outer radius of the 3D source-ball domain in metres. Defaults
    /// to a multiple of the elastic radius set by the host solver.
    double outer_radius_m = -1.0;

    /// Cavity radius of the 3D source-ball domain in metres. Pass-13a
    /// foundation uses this only as a sanity check against the loaded
    /// mesh's minimum vertex radius. Pass-13b uses it to set the
    /// inner free-stress boundary.
    double cavity_radius_m = -1.0;

    /// Approximate cell size for the unstructured-tet mesh in metres.
    /// Smaller -> more cells, finer resolution. The pass-13 baseline
    /// uses ~0.5 * cavity-radius cell size near the cavity and graded
    /// outward.
    double mesh_cell_size_m = -1.0;

    /// Path *without* extension to the TetGen-generated mesh files.
    /// Source3DBallImpl::initialize() reads `mesh_path + ".node"` and
    /// `mesh_path + ".ele"`. Empty disables mesh loading (foundation
    /// no-op mode used by the factory-instantiation regression test).
    std::string mesh_path;

    /// When true, the source ball is loaded with a depth-dependent
    /// initial stress matching the overburden profile of the host
    /// rock column. Pass-13a foundation honours the field but does
    /// not yet apply the IC; pass-13b implements the asymmetric
    /// overburden initial state.
    bool asymmetric_overburden = true;

    /// At-rest Earth coefficient for the asymmetric overburden initial
    /// stress: sigma_xx = sigma_yy = K_0 * sigma_zz with sigma_zz the
    /// lithostatic component. Default 0.5 matches typical crustal-rock
    /// values (Hoek & Brown 1980); only consulted when
    /// asymmetric_overburden is true. Pass-13a stores the value;
    /// pass-13b applies it.
    double overburden_K0 = 0.5;

    /// Pass-13b material parameters used by the per-cell constitutive
    /// update. Defaults reproduce a granite-class host rock; the host
    /// (RadialLagrangianSolver) overrides these from the configured
    /// medium during cavity_geometry = THREE_DIMENSIONAL delegation.
    double host_density_kg_per_m3 = 2700.0;
    double bulk_modulus_K_pa = 30.0e9;
    double shear_modulus_G_pa = 18.0e9;
    double source_depth_m = 500.0;
    /// Per-cell heat capacity at constant volume [J/(kg*K)]. Used by
    /// the matter-coupling update inside the radiation diffusion
    /// solver.
    double cv_J_per_kg_K = 1000.0;
    /// Initial (and ambient) matter temperature [K]. Used to seed the
    /// per-cell T_m field and as the radiation-diffusion boundary
    /// reservoir.
    double T_ambient_K = 300.0;
    /// Drucker-Prager 3D yield-surface parameters. Defaults to granite
    /// (DruckerPrager3DSets::granite()). Set by the host's medium
    /// dispatch during delegation.
    double dp3d_alpha = 0.30;
    double dp3d_k_pa = 70.0e6;
    /// Medium label, propagated to the diagnostics. "GRANITE" by
    /// default. Used by the host to select dp3d / opacity parameter
    /// sets.
    std::string medium_label = "GRANITE";

    /// 3D constitutive model selector. Pass-11 enumerated only the
    /// names; pass-13b will populate the radial-return implementation.
    enum class Constitutive
    {
        DRUCKER_PRAGER_3D,         ///< Pass-13b target
        VON_MISES_3D,              ///< Test only
        ELASTIC_LINEAR_3D          ///< For pure-elastic V&V
    };
    Constitutive constitutive = Constitutive::DRUCKER_PRAGER_3D;

    /// Radiation block selector inside the 3D source ball. Pass-13a
    /// foundation only stores the choice; pass-13b implements the
    /// cell-centred FV grey-diffusion path. NODAL_FEM remains scaffold
    /// and Source3DBallImpl::initialize() throws a clear pass-15+
    /// message when it is selected.
    enum class RadiationDiscretization
    {
        FEM_NODAL,
        FV_CELL_CENTRED
    };
    RadiationDiscretization radiation_discretization =
        RadiationDiscretization::FV_CELL_CENTRED;

    /// Pass-13c performance knob: run the radiation diffusion sub-step
    /// only every Nth call to `Source3DBallImpl::step()`. The 1D host
    /// runs at ~10 us substeps for CFL on the radial mesh; the
    /// radiation diffusion is implicit (unconditionally stable) and
    /// physically only needs ~1 ms cadence to resolve the radiation-
    /// transport timescale at Salmon-class yields. Default 1 (every
    /// step) preserves the pass-13b numerical answer; set to 100 in
    /// large-mesh end-to-end runs.
    int radiation_substep_cadence = 1;

    // ----------------------------------------------------------------
    // Pass-14a (axis-1b source-forcing slice) -- volumetric energy
    // deposition that drives the matter pressure (Tillotson EOS) inside
    // the inner-cavity cells and radiates an acoustic pressure pulse
    // outward to the elastic-radius extraction surface. Before pass-14a
    // the 3D path had no internal source term: the cavity wall sat
    // motionless and the surface-integral moment tensor stayed zero.
    // ----------------------------------------------------------------

    /// Source-time-function shapes for the deposited mechanical energy.
    /// All shapes are normalised so the time integral of the rate is
    /// unity; the rate is multiplied by the total deposited energy
    /// (source_yield_kt * 4.184e12 J/kt * source_deposition_efficiency)
    /// to give J/s. MUELLER_MURPHY / BRUNE use the critically-damped
    /// reduced-displacement-potential pulse shape (Mueller & Murphy
    /// 1971 BSSA 61(6); Brune 1970 JGR 75). RAMP is a boxcar (linear
    /// cumulative energy over the deposition window). DELTA dumps the
    /// whole energy in the first step.
    enum class SourceTimeFunction
    {
        MUELLER_MURPHY,
        BRUNE,
        RAMP,
        DELTA
    };

    /// Master switch for the source forcing. The host (RadialLagrangian
    /// delegation) defaults this to true when cavity_geometry =
    /// THREE_DIMENSIONAL and false otherwise. With source_forcing_enabled
    /// = false the 3D path reproduces pass-13c byte-for-byte (the cavity
    /// wall stays motionless, the moment tensor stays zero) so the
    /// pass-13c "pipeline completes" regression guard still applies.
    bool source_forcing_enabled = false;

    /// Source-time-function shape selector.
    SourceTimeFunction source_time_function = SourceTimeFunction::MUELLER_MURPHY;

    /// Yield in kilotons of TNT. The deposited mechanical energy is
    /// source_yield_kt * 4.184e12 J/kt * source_deposition_efficiency.
    /// Negative -> the host fills this from the configured explosion
    /// yield during delegation; if still non-positive at initialize()
    /// the source forcing is disabled with a one-time warning.
    double source_yield_kt = -1.0;

    /// Full timescale of the source-time function [s]. Mueller-Murphy
    /// default 1 ms; this sets the deposition rate, not the eventual
    /// seismic corner frequency.
    double source_deposition_duration_s = 1.0e-3;

    /// Fraction of the yield deposited as mechanical (cavity) energy.
    /// 1.0 deposits the whole nominal yield; smaller values model the
    /// radiative / vaporization losses. Clamped to (0, 1] at initialize.
    double source_deposition_efficiency = 1.0;
};

/// Snapshot of the 3D source-ball state at a single time. Used by the
/// host RadialLagrangianSolver-equivalent to record the moment-tensor
/// volume integral and spatial profile snapshots. Pass-12 populates
/// the fields with actual data.
struct Source3DBallState
{
    double time = 0.0;
    /// 6-component Cartesian moment tensor [Mxx, Myy, Mzz, Mxy, Mxz,
    /// Myz] in N*m. Pass-12 distinguishes asymmetric content
    /// (CLVD, double-couple) from the isotropic baseline.
    std::array<double, 6> moment_tensor_N_m = {0, 0, 0, 0, 0, 0};
    std::array<double, 6> moment_rate_tensor_N_m_per_s =
        {0, 0, 0, 0, 0, 0};
    /// Maximum cavity-face displacement [m]. Pass-12 reports per-
    /// hemisphere extrema for the asymmetry diagnostic.
    double cavity_radius_max_m = 0.0;
    double cavity_radius_min_m = 0.0;
};

/**
 * @brief Pass-11 scaffolded interface for the axis-1b 3D source ball.
 *
 * Pass-12 implementation will provide:
 *   - Unstructured-tet mesh generation conforming to a spherical
 *     outer boundary at the elastic radius (TetGen / Gmsh API).
 *   - 3D Drucker-Prager constitutive with explicit radial-return.
 *   - Depth-dependent overburden stress as the initial state.
 *   - 3D radiation block (FEM or cell-centred FV; design choice
 *     deferred to docs/AXIS_1B_DESIGN.md).
 *   - Surface-integral moment-tensor extraction at the elastic radius
 *     for handoff to the far-field FEM (axis-1d).
 *
 * This interface is the contract between the host driver and the
 * eventual pass-12 concrete subclass. Pass-11 keeps the interface
 * minimal so pass-12 can extend it without breaking callers.
 */
class Source3DBall
{
public:
    virtual ~Source3DBall() = default;

    /// Configure the 3D source-ball solver. Throws std::runtime_error
    /// in pass-11 because no concrete implementation exists.
    virtual void initialize(const Source3DBallConfig& cfg) = 0;

    /// Advance the 3D state by dt. The host's outer loop calls this
    /// with the same cadence as the pass-10 1D radial step().
    virtual void step(double dt) = 0;

    /// Volume-integrated 6-component Cartesian moment tensor at the
    /// current time. Used by the far-field FEM coupling at the
    /// elastic radius.
    virtual void getMomentTensor(std::array<double, 6>& M) const = 0;

    /// Volume-integrated 6-component Cartesian moment-rate tensor.
    virtual void getMomentRateTensor(std::array<double, 6>& Mdot) const = 0;

    /// Snapshot of the 3D state at the current time. Used by the
    /// pass-12 HDF5 / XDMF spatial-profile writer.
    virtual void getState(Source3DBallState& state) const = 0;

    /// Diagnostic name. Pass-11 returns the pass-12 placeholder
    /// "axis_1b_3d_source_ball_pass12_pending".
    virtual const char* name() const = 0;
};

/// Pass-11 throw-on-construct factory placeholder. Pass-12 will
/// replace the implementation with a real 3D-mesh-and-FEM subclass.
std::unique_ptr<Source3DBall> makeSource3DBall(
    const Source3DBallConfig& cfg);

}  // namespace FSRM

#endif  // NEAR_FIELD_SOURCE_3D_BALL_HPP
