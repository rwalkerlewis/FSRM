/**
 * @file Source3DBallImpl.hpp
 * @brief Pass-13b (axis-1b physics) concrete implementation of the
 *        Source3DBall interface. Pass-13a foundation (PR #129) loaded
 *        the mesh; pass-13b adds the matter and radiation update.
 *
 * What pass-13b delivers (the physics slice):
 *  - initialize(cfg): loads the mesh, builds per-cell state, applies the
 *    asymmetric overburden initial state via SourceBallOverburden, and
 *    builds the cell / face POD lists for the 3D radiation diffusion
 *    solver.
 *  - step(dt): advances each owned cell through the 3D Drucker-Prager
 *    radial return (DruckerPrager3D) using the host-supplied strain
 *    rate, then advances the matter / radiation state through the
 *    SourceBallRadiation3DSolver. Replaces the pass-13a throw.
 *  - setCellStrainRate(): hook used by the host (RadialLagrangianSolver
 *    delegation in pass-13b; pass-13c will couple the rate through the
 *    surface-integral moment-tensor extraction). Defaults to zero per
 *    cell when not explicitly set.
 *  - getMomentTensor / getMomentRateTensor: pass-13b returned the zero
 *    tensor; pass-13c integrates the cell-centred stress drop over the
 *    elastic-radius surface (Day & McLaughlin 1991; Aki & Richards
 *    2002 ch 4), producing real M(t) and Mdot(t) signals.
 *  - getState: returns a snapshot with mesh stats and energy diagnostics
 *    but the 6-component moment tensor remains zero.
 *  - getCellState (new in pass-13b): per-cell snapshot used by the
 *    physics-validation gates.
 *  - name(): "Source3DBallImpl_v1_pass13b_physics".
 *
 * What pass-13c adds on top of pass-13b:
 *  - Surface-integral moment-tensor extraction. Each rank sums
 *    M_ij = integral_S (sigma_jk(t) - sigma_jk(0)) * n_k * x_i dA over
 *    its owned elastic-radius boundary faces; MPI_Allreduce produces
 *    the global tensor. Citation: Day & McLaughlin 1991 sec 4 and
 *    Aki & Richards 2002 ch 4.
 *  - getMomentTensor / getMomentRateTensor return real signals.
 *  - HDF5 / XDMF spatial-profile output (writeSpatialSnapshot).
 *  - Cavity vertical / horizontal radius measurement (no gate; reported
 *    in PR body and HDF5 metadata).
 *  - name() bumps to "Source3DBallImpl_v1_pass13c_validation".
 *
 * What pass-13c explicitly does not implement (named for follow-on
 * passes; throws or no-ops with clear messages):
 *  - 3D Lagrangian face advection / mass-conservation hydro (pass-14+;
 *    the pass-13c driver still expects strain rates from the host, the
 *    same pattern the pass-10 1D radial solver uses internally).
 *  - 3D far-field FEM coupling (axis-1d, future pass).
 *
 * Backward compatibility: cavity_geometry = SPHERICAL remains the
 * default and entirely bypasses this code path. The 32 historic-event
 * integration tests do not exercise Source3DBallImpl. Pass-13b also
 * lands the host-side delegation (RadialLagrangianSolver constructs
 * Source3DBallImpl when cavity_geometry = THREE_DIMENSIONAL).
 */

#ifndef NEAR_FIELD_SOURCE_3D_BALL_IMPL_HPP
#define NEAR_FIELD_SOURCE_3D_BALL_IMPL_HPP

#include "domain/explosion/Source3DBall.hpp"
#include "domain/explosion/DruckerPrager3D.hpp"
#include "domain/explosion/SourceBallRadiation3D.hpp"
#include "domain/explosion/SourceBallOverburden.hpp"

#include <array>
#include <memory>
#include <vector>
#include <mpi.h>

namespace FSRM {

class Source3DBallMesh;

/// Per-cell state snapshot exposed by Source3DBallImpl::getCellStates.
/// Used by the physics-validation gates. Length = number of locally
/// owned cells.
struct Source3DCellState
{
    std::array<double, 3> centroid;          ///< Local frame [m].
    std::array<double, 6> sigma;              ///< Voigt stress [Pa].
    double eps_p_eq = 0.0;                    ///< Equivalent plastic strain.
    double rho = 0.0;                         ///< Density [kg/m^3].
    double e_int = 0.0;                       ///< Specific internal energy [J/kg].
    double T_m = 0.0;                          ///< Matter temperature [K].
    double E_r = 0.0;                          ///< Radiation energy density [J/m^3].
};

class Source3DBallImpl final : public Source3DBall
{
public:
    Source3DBallImpl();
    ~Source3DBallImpl() override;

    Source3DBallImpl(const Source3DBallImpl&) = delete;
    Source3DBallImpl& operator=(const Source3DBallImpl&) = delete;

    /// Override the MPI communicator the underlying DMPlex lives on.
    /// Defaults to PETSC_COMM_WORLD when initialize() is called
    /// without first calling this. Pass-13b uses this from the host's
    /// existing communicator selection logic.
    void setComm(MPI_Comm comm);

    // Source3DBall interface --------------------------------------

    void initialize(const Source3DBallConfig& cfg) override;
    void step(double dt) override;
    void getMomentTensor(std::array<double, 6>& M) const override;
    void getMomentRateTensor(std::array<double, 6>& Mdot) const override;
    void getState(Source3DBallState& state) const override;
    const char* name() const override;

    // Pass-13a foundation diagnostics -----------------------------

    /// True after initialize() has loaded a mesh successfully.
    bool isInitialized() const { return mesh_loaded_; }

    /// Pointer to the underlying mesh. Null before initialize() is
    /// called or after initialize() returned without a mesh_path.
    const Source3DBallMesh* mesh() const { return mesh_.get(); }

    // Pass-13b physics accessors ----------------------------------

    /// Per-cell strain rate (Voigt 6 with engineering shear in 3..5).
    /// Used by the host to drive the constitutive update from external
    /// kinematics (the pass-13b host delegation passes the spherical-
    /// symmetric rates from the 1D solver in this configuration). The
    /// vector must have one entry per local cell; defaults to zero.
    void setCellStrainRate(const std::vector<std::array<double, 6>>& rates);

    /// Set per-cell density (e.g. constant host-rock density). Length
    /// must equal numLocalCells(). Defaults to cfg.host_density at
    /// initialize().
    void setCellDensity(const std::vector<double>& rho);

    /// Test-only: directly write the per-cell Voigt stress field. Used
    /// by the pass-13c surface-integral unit gates to drive the
    /// moment-tensor computation from a synthetic (analytically
    /// tractable) stress state. Length must equal numLocalCells().
    void setCellStressForTest(const std::vector<std::array<double, 6>>& sig);

    /// Test-only: reset the cached sigma_initial_ snapshot to the
    /// supplied per-cell stress field (or zero by default). Used so the
    /// surface-integral gates can dial the reference state independent
    /// of the overburden IC.
    void setCellStressInitialForTest(
        const std::vector<std::array<double, 6>>& sig0);

    /// Number of locally owned cells (after initialize). 0 in foundation
    /// no-op mode.
    int numLocalCells() const { return n_local_cells_; }

    /// Per-cell state snapshot. Returns an empty vector before
    /// initialize() or in foundation no-op mode.
    void getCellStates(std::vector<Source3DCellState>& out) const;

    /// Direct access to the per-cell stress, plastic-strain, density,
    /// internal-energy, matter-temperature, and radiation-energy
    /// vectors used internally. Used by pass-13b unit gates that
    /// inspect post-step state without copying through getCellStates.
    const std::vector<std::array<double, 6>>& sigmaCells() const { return sigma_; }
    const std::vector<double>& eps_p_eq_cells() const { return eps_p_eq_; }
    const std::vector<double>& rhoCells() const { return rho_; }
    const std::vector<double>& eIntCells() const { return e_int_; }
    const std::vector<double>& TmCells() const { return T_m_; }
    const std::vector<double>& ErCells() const { return E_r_; }

    /// Diagnostic: per-step Newton iterations of the radiation solver.
    int lastRadiationNewtonIters() const { return last_rad_newton_iters_; }

    /// Diagnostic: number of cells that yielded in the most recent step.
    int lastYieldedCellCount() const { return last_yielded_count_; }

    /// Configured cavity radius (from `Source3DBallConfig::cavity_radius_m`
    /// at initialize time). Used by the host's getCavityRadius shim
    /// which now early-outs for the 3D path.
    double cfg_cavity_radius_m() const { return cfg_.cavity_radius_m; }

    /// Number of locally owned elastic-radius surface faces (pass-13c).
    int numLocalElasticSurfaceFaces() const
    { return static_cast<int>(elastic_surface_faces_.size()); }

    /// Pass-13c cavity radius diagnostic. Returns the maximum and
    /// minimum vertex distance from origin among vertices marked as
    /// cavity surface, MPI-reduced over the communicator. Both equal
    /// the initial cavity radius before any displacement is added; in
    /// pass-13c the cavity is rigid in the no-advection architecture
    /// (cavity_max == cavity_min). Pass-14 will populate the asymmetry
    /// once 3D Lagrangian advection lands.
    void cavityRadiusExtremes(double& r_min_out, double& r_max_out) const;

    /// Force a recomputation of the surface-integral moment tensor from
    /// the current per-cell sigma_ state. Called by step() after the
    /// constitutive + radiation update. Exposed so tests can drive the
    /// integral from a synthetic stress field.
    void recomputeMomentTensorFromSurface(double dt);

private:
    void buildCellAndFaceLists();
    void applyOverburdenIC();
    void cacheElasticSurfaceGeometry();

    MPI_Comm comm_ = MPI_COMM_NULL;
    bool comm_set_ = false;
    bool mesh_loaded_ = false;
    Source3DBallConfig cfg_{};
    std::unique_ptr<Source3DBallMesh> mesh_;
    Source3DBallState state_snapshot_{};

    // Pass-13b per-cell state vectors. Length = n_local_cells_.
    int n_local_cells_ = 0;
    std::vector<std::array<double, 3>> centroid_;
    std::vector<double> volume_;
    std::vector<std::array<double, 6>> sigma_;
    std::vector<std::array<double, 6>> strain_rate_;
    std::vector<double> eps_p_eq_;
    std::vector<double> rho_;
    std::vector<double> e_int_;
    std::vector<double> T_m_;
    std::vector<double> E_r_;

    // Mesh topology lists for the radiation solver (built from DMPlex
    // once during initialize).
    std::vector<CellGeometry3D> rad_cells_;
    std::vector<InternalFace3D> rad_internal_faces_;
    std::vector<BoundaryFace3D> rad_boundary_faces_;
    std::unique_ptr<SourceBallRadiation3DSolver> rad_solver_;

    DruckerPrager3DParameters dp_params_;

    int last_rad_newton_iters_ = 0;
    int last_yielded_count_ = 0;
    double current_time_ = 0.0;
    int rad_substep_counter_ = 0;

    // Pass-13c surface-integral moment-tensor extraction.
    /// Per-cell stress at the start of the run (before any step() call).
    /// Used as the reference state for the stress-drop integral
    /// M_ij = integral_S [sigma_jk(t) - sigma_jk(0)] * n_k * x_i dA
    /// (Day & McLaughlin 1991; Aki & Richards 2002 ch 4).
    std::vector<std::array<double, 6>> sigma_initial_;

    /// Per-elastic-surface-face geometry cache. cell_local indexes into
    /// sigma_ / sigma_initial_; outward_normal points away from the
    /// owning cell's centroid; centroid is the face centroid in the
    /// local mesh frame (origin = source center).
    struct ElasticSurfaceFace
    {
        int cell_local = -1;
        double area = 0.0;
        std::array<double, 3> centroid = {0.0, 0.0, 0.0};
        std::array<double, 3> outward_normal = {0.0, 0.0, 0.0};
    };
    std::vector<ElasticSurfaceFace> elastic_surface_faces_;

    /// Global cumulative moment tensor [Mxx, Myy, Mzz, Mxy, Mxz, Myz]
    /// in N*m. Updated by recomputeMomentTensorFromSurface() each step.
    std::array<double, 6> M_global_ = {0, 0, 0, 0, 0, 0};
    /// Global instantaneous moment-rate tensor in N*m/s. Computed by
    /// finite difference of M_global_ at successive step() calls.
    std::array<double, 6> Mdot_global_ = {0, 0, 0, 0, 0, 0};
    std::array<double, 6> M_global_prev_ = {0, 0, 0, 0, 0, 0};

    // Material parameters used by the constitutive update. Supplied
    // through Source3DBallConfig in pass-13b (defaults match a granite-
    // class host rock).
    double bulk_modulus_K_ = 30.0e9;
    double shear_modulus_G_ = 18.0e9;
};

}  // namespace FSRM

#endif  // NEAR_FIELD_SOURCE_3D_BALL_IMPL_HPP
