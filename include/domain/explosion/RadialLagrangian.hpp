/**
 * @file RadialLagrangian.hpp
 * @brief 1D radial Lagrangian finite-volume elastoplastic shock solver.
 *
 * Pass-6 (axis-1a, see docs/HISTORIC_NUCLEAR_ROADMAP.md) replaces the
 * closed-form exponential cavity-expansion kernel inside
 * NearFieldExplosionSolver with a real shock-physics solve in spherical
 * symmetry: cell-centred conserved state, face-centred velocity (staggered),
 * Wilkins linear+quadratic artificial viscosity for shock capture, explicit
 * predictor for deviatoric stress, Drucker-Prager radial-return per cell,
 * Mie-Gruneisen EOS evaluation per cell, energy update including plastic
 * dissipation, damage evolution, non-reflecting characteristic outer BC,
 * and surface-integral moment-tensor extraction at a fixed Eulerian
 * elastic radius.
 *
 * The class is intentionally decoupled from NearFieldExplosionSolver: the
 * Simulator's pass-5 setup loop chooses between the existing closed-form
 * kernel (rdp_source.momentRateTensor) and this new solver based on the
 * solver_kind sub-key under [NEAR_FIELD_SOURCE]. CLOSED_FORM preserves
 * the pass-5 published-test behaviour exactly; RADIAL_LAGRANGIAN runs
 * this solver and substitutes its Mdot at the recording site.
 *
 * State convention. sigma > 0 is tension, sigma < 0 is compression. The
 * scalar pressure p > 0 in compression follows p = -trace(sigma)/3. In
 * spherical symmetry sigma is diagonal in spherical coordinates with
 * sigma_rr and sigma_tt = sigma_phiphi; the deviator obeys
 * s_rr + 2*s_tt = 0 so we track s_rr and derive s_tt = -s_rr/2.
 *
 * References.
 *  - Wilkins (1980), Computer Simulation of Dynamic Phenomena, Springer
 *    (artificial viscosity, hourglassing, radial return formulation).
 *  - Trangenstein (2009), Numerical Solution of Hyperbolic Partial
 *    Differential Equations, ch. 6 (1D Lagrangian schemes).
 *  - Day & McLaughlin (1991), Effects of plasticity on the seismic source
 *    of underground explosions (moment-tensor extraction at the elastic
 *    radius for spherically symmetric explosion sources).
 *  - Sedov (1959), Similarity and Dimensional Methods (point-source blast
 *    self-similar solution used as the no-strength validation target).
 */

#ifndef NEAR_FIELD_RADIAL_LAGRANGIAN_HPP
#define NEAR_FIELD_RADIAL_LAGRANGIAN_HPP

#include <array>
#include <cstddef>
#include <string>
#include <vector>

#include "domain/explosion/NearFieldExplosion.hpp"
#include "domain/explosion/TillotsonEOS.hpp"

namespace FSRM {

/**
 * @brief 1D radial Lagrangian finite-volume elastoplastic shock solver.
 *
 * Public API is small and stable: configure, initialize, step, and a
 * handful of accessors that the Simulator records into the pass-5
 * history vectors. Internal state is exposed through getRadialProfile
 * for HDF5 / XDMF spatial-profile writing during DYNAMIC_PLASTIC runs.
 */
class RadialLagrangianSolver
{
public:
    /// Configuration sub-keys parsed under [NEAR_FIELD_SOURCE]. Defaults
    /// reproduce a reasonable cavity-formation transient for ~100 kt
    /// nuclear-yield shots in alluvium-class media.
    ///
    /// Pass-7 default switches:
    ///  - Wilkins AV coefficients land at the Wilkins (1980) prescription
    ///    c_l = 0.06, c_q = 1.5 (literature values for production shock
    ///    capture). Pass-6 used 0.5 / 2.0 for early-development
    ///    robustness which over-dissipated the leading shock.
    ///  - cavity_eos defaults to TILLOTSON: the inner cavity is
    ///    fully-vaporized rock plasma at temperatures four to six
    ///    orders of magnitude above chemical-detonation conditions,
    ///    NOT chemical detonation gas at gamma = 1.4. The Tillotson
    ///    parameter set is selected from the medium type unless
    ///    explicitly overridden.
    ///  - cavity_initialization defaults to PHYSICS_BASED: the
    ///    initial cavity radius and (rho, e) state are computed by
    ///    a Newton solve on a closed energy-partition equation
    ///    rather than chosen by hand. See solveCavityInitialState
    ///    in RadialLagrangian.cpp.
    enum class CavityEOS { IDEAL_GAS, TILLOTSON };
    enum class CavityInitialization { PHYSICS_BASED, MANUAL };

    struct Config
    {
        int radial_cells = 200;             ///< Number of FV cells along radius.
        double radial_outer_factor = 3.0;   ///< Outer radius / elastic radius.
        double cfl = 0.4;                   ///< CFL number on max(c_p + |v|).
        double art_visc_linear = 0.06;      ///< Wilkins linear AV coefficient (pass-7 literature default).
        double art_visc_quadratic = 1.5;    ///< Wilkins quadratic AV coefficient (pass-7 literature default).
        double gas_eos_gamma = 1.4;         ///< Adiabatic index when cavity_eos = IDEAL_GAS (pass-6 placeholder).
        int profile_output_cadence_us = 100;///< Spatial-profile snapshot cadence.

        /// Pass-7 cavity inner-state controls.
        CavityEOS cavity_eos = CavityEOS::TILLOTSON;
        TillotsonParameters tillotson_params =
            TillotsonParameterSets::granite();
        CavityInitialization cavity_initialization =
            CavityInitialization::PHYSICS_BASED;
        /// MANUAL only: explicit initial cavity radius [m]. Ignored
        /// under PHYSICS_BASED.
        double initial_cavity_radius_m = 0.0;
        /// Radiation-to-hydrodynamic transition time [s]. The cavity
        /// initial state under PHYSICS_BASED is set at this point.
        /// Default derived from yield via Zel'dovich-Raizer scaling
        /// in solveCavityInitialState() when negative or zero.
        double radiation_transition_time_s = -1.0;

        /// When true, plasticity / strength terms are zeroed: the solver
        /// becomes a pure-elastic shock solver. Used by the
        /// PureElasticSphericalWave / OutgoingBC / Sedov physics tests.
        bool disable_plasticity = false;
    };

    /// Snapshot of the radial state at a single time. Layout matches the
    /// HDF5 group structure (/profiles/<i>/{r,v_r,rho,p,sigma_rr,sigma_tt,
    /// eps_p,damage,yield_indicator}) the Simulator writes during pass-6
    /// DYNAMIC_PLASTIC runs.
    struct RadialProfile
    {
        double time = 0.0;
        std::vector<double> r_face;        ///< (N+1) face radii [m].
        std::vector<double> r_cell;        ///< (N) cell-centred radii [m].
        std::vector<double> v_r;           ///< (N+1) face velocities [m/s].
        std::vector<double> rho;           ///< (N) density [kg/m^3].
        std::vector<double> p;             ///< (N) pressure [Pa], compression positive.
        std::vector<double> sigma_rr;      ///< (N) total radial stress [Pa].
        std::vector<double> sigma_tt;      ///< (N) total hoop stress [Pa].
        std::vector<double> eps_p;         ///< (N) equivalent plastic strain.
        std::vector<double> damage;        ///< (N) scalar damage [0,1].
        std::vector<double> yield_indicator;///< (N) 1 if yielded this step, else 0.
    };

    RadialLagrangianSolver();
    ~RadialLagrangianSolver() = default;

    void setSource(const UndergroundExplosionSource& src);
    void setEOS(const MieGruneisenEOS& eos);
    void setStrength(const PressureDependentStrength& strength);
    void setDamage(const DamageEvolutionModel& damage);
    void setConfig(const Config& cfg);

    /// Build the spherical mesh, deposit yield-derived energy in the
    /// initial cavity, set hydrostatic stress in solid cells, zero
    /// velocities. Idempotent.
    void initialize();

    /// Advance the radial state by dt_target. The CFL bound may force
    /// internal sub-stepping; the public time advances by exactly
    /// dt_target on return so the Simulator's outer loop can keep its
    /// own cadence.
    void step(double dt_target);

    double getCurrentTime() const { return current_time_; }

    /// Cavity radius at the current time. Tracks the inner face position
    /// of the gas / solid interface.
    double getCavityRadius() const;

    /// Elastic radius (the fixed Eulerian sphere where moment-tensor
    /// surface integrals are evaluated). Set at initialize() from the
    /// source-derived NTS analytic.
    double getElasticRadius() const { return r_elastic_; }

    /// Plastic radius at the current time: outermost cell that has
    /// accumulated nonzero equivalent plastic strain.
    double getPlasticRadius() const;

    /// 6-component Cartesian moment-rate tensor (Voigt: xx,yy,zz,xy,xz,yz)
    /// from the spherical-symmetry surface integral at the elastic
    /// radius. By construction the solver produces a purely diagonal
    /// isotropic tensor; CLVD content requires axis-1b (3D subdomain).
    void getMomentRateTensor(std::array<double, 6>& Mdot) const;

    /// 6-component Cartesian moment tensor (time-integrated Mdot).
    void getMomentTensor(std::array<double, 6>& M) const;

    /// Diagnostic: kinetic, internal, plastic-dissipated, and radiated
    /// energy crossing the outer boundary. Each is in joules and the
    /// sum (plus the residual elastic strain energy) should match the
    /// initially-deposited yield to within the Wilkins-AV dissipation
    /// budget.
    double getKineticEnergy() const { return kinetic_energy_; }
    double getInternalEnergy() const { return internal_energy_; }
    double getPlasticDissipation() const { return plastic_dissipation_; }
    double getRadiatedEnergyOut() const { return radiated_energy_out_; }
    double getInitialDepositedEnergy() const { return initial_energy_; }

    /// Snapshot of the spherically-symmetric radial state. Used by the
    /// pass-6 HDF5 / XDMF writer in Simulator.cpp and by the
    /// MeshRefinementConvergence test to compare resolutions.
    void getRadialProfile(RadialProfile& profile) const;

    /// Number of cells. Useful in tests for refinement studies.
    int getNumCells() const { return N_; }

    /// Pass-7 diagnostic accessors for the physics-based cavity
    /// initialization (PhysicsBasedCavityEnergyConservation /
    /// PhysicsBasedCavityRadius validation gates):
    ///   getInitialCavityRadius - the R_v solved by the energy
    ///     partition Newton iteration (or the user-supplied value
    ///     under MANUAL initialization).
    ///   getInitialCavityVaporSpecificEnergy - the e_v at t = t_rh
    ///     under PHYSICS_BASED. Zero under MANUAL.
    ///   getRadiationTransitionTime - the t_rh used. Either the
    ///     user-supplied value or the Zel'dovich-Raizer-scaled
    ///     default.
    double getInitialCavityRadius() const { return Rc_init_; }
    double getInitialCavityVaporSpecificEnergy() const { return e_v_init_; }
    double getRadiationTransitionTime() const { return t_rh_used_; }

private:
    void allocate(int N);
    void cflLimit(double& dt) const;
    void substep(double dt);

    void advanceFaces(double dt);
    void updateDensity();
    void computeArtificialViscosity(double dt);
    void momentumUpdate(double dt);
    void deviatoricElasticPredictor(double dt);
    void radialReturnPlasticity(double dt);
    void updateInternalEnergy(double dt);
    void updateEOS();
    void updateDamage(double dt);
    void absorbingOuterBC();
    void recordMomentExtraction();

    /// Pass-7: solve for the (R_v, rho_v, e_v) inner-cavity state at
    /// t = t_rh using a first-principles energy partition. Sets
    /// Rc_init_, rho_v_init_, e_v_init_, and t_rh_used_ in place.
    /// Implementation in RadialLagrangian.cpp; see the comment block
    /// there for the algorithm and citations (Zel'dovich-Raizer 1967,
    /// Melosh 1989).
    void solveCavityInitialState();

    /// Pass-7: cavity-cell pressure evaluation. Branches on
    /// config_.cavity_eos: IDEAL_GAS reproduces the pass-6 placeholder
    /// (p = (gamma - 1) rho e); TILLOTSON evaluates the configured
    /// parameter set against (rho, e) of the cavity cell.
    double cavityPressure(double rho, double e) const;

    static double cellVolumeSpherical(double r_lo, double r_hi);
    static double faceAreaSpherical(double r);

    UndergroundExplosionSource src_;
    MieGruneisenEOS eos_;
    PressureDependentStrength strength_;
    DamageEvolutionModel damage_model_;
    Config config_;

    int N_ = 0;
    int gas_cells_ = 1;       ///< Number of inner cells assigned the gas EOS.
    double r_elastic_ = 0.0;  ///< Fixed Eulerian extraction sphere [m].
    double r_outer_ = 0.0;    ///< Initial outer face position [m].
    double Rc_init_ = 0.0;    ///< Initial inner cavity radius [m].

    /// Pass-7 cavity-init bookkeeping. rho_v_init_ and e_v_init_ are
    /// the (rho, e) state imposed on the inner cavity cells under
    /// PHYSICS_BASED initialization (zero / unset otherwise);
    /// t_rh_used_ is the radiation-to-hydrodynamic transition time
    /// (either user-supplied or Zel'dovich-Raizer scaling default).
    double rho_v_init_ = 0.0;
    double e_v_init_ = 0.0;
    double t_rh_used_ = 0.0;
    /// Cached Tillotson EOS for the cavity. setConfig copies
    /// config_.tillotson_params into this evaluator so the inner-cell
    /// pressure update does not allocate a new EOS each step.
    TillotsonEOS cavity_tillotson_;

    // Material parameters (cached at initialize for speed).
    double mu_solid_ = 0.0;    ///< Shear modulus G in solid cells [Pa].
    double K_solid_ = 0.0;     ///< Bulk modulus K in solid cells [Pa].
    double cp_ref_ = 0.0;      ///< Reference P-wave speed for outer BC [m/s].

    // Cell-centred conserved quantities, length N.
    std::vector<double> r_cell_;
    std::vector<double> mass_;
    std::vector<double> rho_;
    std::vector<double> e_int_;     ///< Specific internal energy [J/kg].
    std::vector<double> p_;         ///< Pressure (compression positive) [Pa].
    std::vector<double> s_rr_;      ///< Deviatoric radial stress [Pa].
    std::vector<double> q_visc_;    ///< Wilkins artificial viscosity [Pa].
    std::vector<double> eps_p_;     ///< Equivalent plastic strain.
    std::vector<double> damage_;    ///< Scalar damage [0,1].
    std::vector<int> yielded_;      ///< 1 if cell yielded in latest step.
    std::vector<int> is_gas_;       ///< 1 for inner gas cell, 0 for solid.

    // Face-centred state, length N+1.
    std::vector<double> r_face_;
    std::vector<double> v_face_;

    double current_time_ = 0.0;

    // Pass-5 byte-identical: the recorded Mdot is the FULL Cartesian
    // 6-vector. From the 1D radial solver the deviatoric components are
    // identically zero (CLVD / DC requires asymmetric 3D physics, axis-1b).
    std::array<double, 6> M_iso_ = {0, 0, 0, 0, 0, 0};
    std::array<double, 6> Mdot_iso_ = {0, 0, 0, 0, 0, 0};
    double sigma_rr_extract_prev_ = 0.0;
    double sigma_rr_extract_initial_ = 0.0;
    bool extract_initialized_ = false;

    // Energy accounting.
    double initial_energy_ = 0.0;
    double kinetic_energy_ = 0.0;
    double internal_energy_ = 0.0;
    double plastic_dissipation_ = 0.0;
    double radiated_energy_out_ = 0.0;

    bool initialized_ = false;
};

} // namespace FSRM

#endif // NEAR_FIELD_RADIAL_LAGRANGIAN_HPP
