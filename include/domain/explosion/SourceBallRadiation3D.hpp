/**
 * @file SourceBallRadiation3D.hpp
 * @brief Pass-13b (axis-1b physics) cell-centred FV grey radiation
 *        diffusion on the unstructured 3D source-ball mesh.
 *
 * Pass-13b ships GREY only. Pass-14 multigroup port is named-only here.
 *
 * Physics (matching the pass-8 1D MARSHAK_GREY problem statement, just
 * on 3D unstructured cells instead of a 1D radial layout):
 *
 *   dE_r/dt = div( D grad E_r ) + c kappa_P rho ( a T_m^4 - E_r )
 *   rho cv dT_m/dt = c kappa_P rho ( E_r - a T_m^4 )
 *
 *   D = c / (3 kappa_R rho) (Rosseland diffusion coefficient)
 *   a = 4 sigma_SB / c (radiation constant)
 *
 * Discretisation: cell-centred FV, two-point flux approximation with
 * harmonic mean of the adjacent face densities, backward Euler in time,
 * Newton outer iteration on the matter T^4 coupling.
 *
 * Mesh access: this solver is decoupled from PETSc DMPlex at the
 * interface. The host (Source3DBallImpl) traverses the DMPlex once
 * during initialise() and presents the cells, internal faces, and
 * boundary faces to the solver as POD lists. This keeps the unit gates
 * simple (a hand-built two-cell mesh suffices).
 *
 * Boundary conditions:
 *   CAVITY: prescribed flux (positive outward into the source ball)
 *           supplied per-cell each step. Pass-13b uses zero by default;
 *           callers (the host's coupled hydro substep) can plug in
 *           the yield-deposition source.
 *   ELASTIC: vacuum Marshak boundary E_r = 0. The diffusion-flux
 *            contribution is computed using the centroid-to-face
 *            distance dist_to_face.
 *
 * Parallelism: the linear solve uses PETSc KSP with bjacobi+sub_lu (the
 * pass-12 parallel-KSP convention). At MPI=1 this collapses to plain LU
 * which the unit gates exercise.
 *
 * References:
 *   - Mihalas, D. and Mihalas, B. W. (1984). Foundations of Radiation
 *     Hydrodynamics. Oxford. Section 81 (Marshak boundary), Section 82
 *     (Rosseland diffusion).
 *   - LeVeque, R. J. (2002). Finite Volume Methods for Hyperbolic
 *     Problems. Cambridge. Section 9 (FV on unstructured meshes).
 *   - PETSc 3.25 manual, "KSP" chapter (assembly + solve patterns).
 */

#ifndef NEAR_FIELD_SOURCE_BALL_RADIATION_3D_HPP
#define NEAR_FIELD_SOURCE_BALL_RADIATION_3D_HPP

#include <petscksp.h>
#include <petscmat.h>
#include <petscvec.h>

#include <array>
#include <vector>

#include "domain/explosion/OpacityModel.hpp"

namespace FSRM {

/// Per-cell geometry needed by the FV diffusion assembly.
struct CellGeometry3D
{
    double volume;
    std::array<double, 3> centroid;
};

/// Internal-face geometry (between two owned cells, or between an owned
/// and a ghost cell). The cell indices are local to this rank; ghost
/// indices use a host-supplied global numbering carried in
/// `cell_right_global` so the FV assembly can place off-diagonal
/// entries into the parallel matrix correctly.
struct InternalFace3D
{
    /// Local cell index on this rank.
    int cell_left_local = -1;
    /// Local cell index on this rank (-1 if cell_right is ghost).
    int cell_right_local = -1;
    /// Global cell index for the right cell. Equal to global numbering
    /// of cell_left when cell_right_local >= 0; the "global" form is
    /// what MatSetValues uses.
    PetscInt cell_left_global = -1;
    PetscInt cell_right_global = -1;
    /// Face area [m^2].
    double area = 0.0;
    /// Distance between the two cell centroids [m].
    double dist_centroid = 1.0;
};

enum class BoundaryFaceType3D
{
    CAVITY,    ///< Inner cavity: prescribed flux (default zero).
    ELASTIC    ///< Outer elastic radius: vacuum Marshak (E_r = 0).
};

struct BoundaryFace3D
{
    int cell_local = -1;
    PetscInt cell_global = -1;
    double area = 0.0;
    double dist_to_face = 1.0;
    BoundaryFaceType3D type = BoundaryFaceType3D::ELASTIC;
};

/// Per-step diagnostic.
struct SourceBallRadiation3DStepResult
{
    bool converged = true;
    int newton_iters = 0;
    double residual_inf_norm = 0.0;
    /// Sum of E_r * volume over locally owned cells [J] before the
    /// MPI_Allreduce. Host can reduce for the global energy.
    double local_radiation_energy_J = 0.0;
    /// Local matter-energy change [J]. Sign positive when radiation
    /// transfers energy into the matter.
    double local_matter_energy_change_J = 0.0;
};

/// Pass-13b grey radiation diffusion solver on a 3D cell-centred FV mesh.
class SourceBallRadiation3DSolver
{
public:
    struct Config
    {
        int max_newton_iter = 12;
        double newton_tolerance = 1.0e-6;
        /// KSP relative tolerance for the inner backward-Euler solve.
        double ksp_rtol = 1.0e-9;
        /// KSP absolute tolerance.
        double ksp_atol = 1.0e-50;
        /// Maximum KSP iterations.
        PetscInt ksp_max_it = 1000;
        /// Per-cell heat capacity at constant volume [J/(kg K)].
        double cv_J_per_kg_K = 1000.0;
        /// Optional uniform constant opacity floor [m^2/kg]. When
        /// `kappa_constant_m2_per_kg` is positive it overrides the
        /// power-law model, matching the pass-8 1D MARSHAK_GREY
        /// CONSTANT-opacity convention.
        double kappa_constant_m2_per_kg = 0.0;
        /// Power-law opacity parameters (used when
        /// kappa_constant_m2_per_kg <= 0).
        PowerLawOpacityParameters opacity_params =
            PowerLawOpacitySets::granite();
    };

    SourceBallRadiation3DSolver();
    ~SourceBallRadiation3DSolver();

    SourceBallRadiation3DSolver(const SourceBallRadiation3DSolver&) = delete;
    SourceBallRadiation3DSolver& operator=(
        const SourceBallRadiation3DSolver&) = delete;

    void setConfig(const Config& cfg) { cfg_ = cfg; }
    const Config& getConfig() const { return cfg_; }

    /// Build PETSc Mat / Vec / KSP from the supplied mesh topology.
    /// `local_offset_global` is this rank's first global cell index in
    /// the natural numbering used by `cell_left_global` /
    /// `cell_right_global` on the supplied face lists. Idempotent: a
    /// second call rebuilds the operator state.
    void initialize(MPI_Comm comm,
                    int n_local_cells,
                    PetscInt local_offset_global,
                    PetscInt n_global_cells,
                    const std::vector<CellGeometry3D>& cells,
                    const std::vector<InternalFace3D>& internal_faces,
                    const std::vector<BoundaryFace3D>& boundary_faces);

    /// Advance the coupled (E_r, T_m, e_int) system by dt. The vectors
    /// are length `n_local_cells`. On return:
    ///   E_r mutated to the post-step radiation energy density.
    ///   T_m mutated to the post-step matter temperature.
    ///   e_int incremented by the per-cell matter-energy change [J/kg].
    /// rho is read-only (held fixed across the radiation substep; the
    /// hydro update happens elsewhere in the host's substep loop).
    SourceBallRadiation3DStepResult step(double dt,
                                         const std::vector<double>& rho,
                                         std::vector<double>& T_m,
                                         std::vector<double>& E_r,
                                         std::vector<double>& e_int);

    /// Diagnostic accessors used by the unit gates.
    int numLocalCells() const { return n_local_; }
    PetscInt numGlobalCells() const { return n_global_; }
    bool isInitialized() const { return initialized_; }
    MPI_Comm getComm() const { return comm_; }

    /// Power-law opacity at (rho, T). Public so the unit gates can
    /// exercise the same kappa(rho, T) the solver uses.
    double rosselandKappa(double rho, double T) const;
    double planckKappa(double rho, double T) const;

private:
    void destroyPetsc();
    /// Assemble the implicit operator at the current Newton iterate.
    /// On entry T_iter holds the current matter temperature iterate.
    /// On exit Mat_ and rhs_ hold the linearised system for E_r^{n+1}.
    void assembleOperator(double dt,
                          const std::vector<double>& rho,
                          const std::vector<double>& T_iter,
                          const std::vector<double>& T_m_old,
                          const std::vector<double>& E_r_old);

    void updateMatterTemperature(double dt,
                                 const std::vector<double>& rho,
                                 const std::vector<double>& T_iter,
                                 const std::vector<double>& E_r_new,
                                 std::vector<double>& T_m,
                                 std::vector<double>& e_int);

    Config cfg_;
    MPI_Comm comm_ = MPI_COMM_NULL;
    bool initialized_ = false;

    int n_local_ = 0;
    PetscInt n_global_ = 0;
    PetscInt local_offset_global_ = 0;

    std::vector<CellGeometry3D> cells_;
    std::vector<InternalFace3D> internal_faces_;
    std::vector<BoundaryFace3D> boundary_faces_;

    Mat A_ = nullptr;
    Vec rhs_ = nullptr;
    Vec sol_ = nullptr;
    KSP ksp_ = nullptr;
    PowerLawOpacity opacity_;
};

}  // namespace FSRM

#endif  // NEAR_FIELD_SOURCE_BALL_RADIATION_3D_HPP
