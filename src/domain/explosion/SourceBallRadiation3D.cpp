/**
 * @file SourceBallRadiation3D.cpp
 * @brief Pass-13b grey-radiation diffusion implementation. See
 *        SourceBallRadiation3D.hpp for the physics statement and the
 *        public contract.
 *
 * Algorithm (per substep):
 *   Newton outer iteration on T_iter:
 *     1. Update kappa_R(rho, T_iter), kappa_P(rho, T_iter) per cell.
 *     2. Linearise the matter source about T_iter using
 *        d(T^4)/dT = 4 T_iter^3.
 *     3. Assemble the cell-centred FV backward-Euler operator on
 *        E_r^{n+1}. For each internal face:
 *          flux_face = D_face * area / dist_centroid * (E_R - E_L)
 *        contributes -tau to the L diagonal, +tau to the L row's R
 *        column, and symmetrically for the R row, with
 *          tau = D_face * area / dist_centroid.
 *        Volume term V/dt on the diagonal, RHS V/dt * E_r_old.
 *        Coupling term: V * c kappa_P rho (a T_m^4 - E_r) linearised:
 *          source = c kappa_P rho * a * T_iter^4
 *                  + c kappa_P rho * a * 4 T_iter^3 (T_m^new - T_iter)
 *                  - c kappa_P rho * E_r^{n+1}
 *        With the inner T_m update solved per cell (small ODE) we
 *        linearise the coefficient and let the per-cell matter update
 *        absorb T_m^{n+1}. The diagonal contribution is
 *          + V * c kappa_P rho (1 + d_factor)
 *        with d_factor a small term encoding the T-dependence; we keep
 *        it implicit on E_r and explicit on T_m, then iterate.
 *     4. KSP solve for E_r^{n+1}.
 *     5. Per-cell update of T_m^{n+1} from energy balance.
 *     6. Convergence check on max |T_iter - T_m^{n+1}|.
 *
 * The ELASTIC boundary face contributes a vacuum Marshak BC: half-flux
 * sink at the outer boundary modeled as a one-sided diffusion to a
 * fictitious cell at E_r = 0 located dist_to_face away. The CAVITY
 * boundary face supplies a prescribed (zero in pass-13b default)
 * inflow, rolled into the RHS as a flux source.
 */

#include "domain/explosion/SourceBallRadiation3D.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace FSRM {

namespace
{

constexpr double C_LIGHT = RadiationConstants::SPEED_OF_LIGHT_M_PER_S;
constexpr double A_RAD = RadiationConstants::RADIATION_CONSTANT_A_J_PER_M3_K4;

inline double harmonicMean(double a, double b)
{
    if (a <= 0.0 || b <= 0.0) return 0.0;
    return 2.0 * a * b / (a + b);
}

}  // namespace

SourceBallRadiation3DSolver::SourceBallRadiation3DSolver() = default;

SourceBallRadiation3DSolver::~SourceBallRadiation3DSolver()
{
    destroyPetsc();
}

void SourceBallRadiation3DSolver::destroyPetsc()
{
    if (ksp_) { KSPDestroy(&ksp_); ksp_ = nullptr; }
    if (A_)   { MatDestroy(&A_);   A_ = nullptr; }
    if (rhs_) { VecDestroy(&rhs_); rhs_ = nullptr; }
    if (sol_) { VecDestroy(&sol_); sol_ = nullptr; }
}

double SourceBallRadiation3DSolver::rosselandKappa(double rho, double T) const
{
    if (cfg_.kappa_constant_m2_per_kg > 0.0) {
        return cfg_.kappa_constant_m2_per_kg;
    }
    return opacity_.rosseland(rho, T);
}

double SourceBallRadiation3DSolver::planckKappa(double rho, double T) const
{
    if (cfg_.kappa_constant_m2_per_kg > 0.0) {
        return cfg_.kappa_constant_m2_per_kg;
    }
    return opacity_.planck(rho, T);
}

void SourceBallRadiation3DSolver::initialize(MPI_Comm comm,
                                             int n_local_cells,
                                             PetscInt local_offset_global,
                                             PetscInt n_global_cells,
                                             const std::vector<CellGeometry3D>& cells,
                                             const std::vector<InternalFace3D>& internal_faces,
                                             const std::vector<BoundaryFace3D>& boundary_faces)
{
    if (static_cast<int>(cells.size()) != n_local_cells) {
        throw std::runtime_error(
            "SourceBallRadiation3DSolver::initialize: cells.size()="
            + std::to_string(cells.size())
            + " does not match n_local_cells=" + std::to_string(n_local_cells));
    }
    destroyPetsc();
    comm_ = comm;
    n_local_ = n_local_cells;
    n_global_ = n_global_cells;
    local_offset_global_ = local_offset_global;
    cells_ = cells;
    internal_faces_ = internal_faces;
    boundary_faces_ = boundary_faces;
    opacity_.setParameters(cfg_.opacity_params);

    // Build a sparse AIJ matrix sized [n_local x n_global]. Each cell
    // has at most 1 (diagonal) + degree(cell) off-diagonals. A loose
    // upper bound on the per-row nonzero count uses 8 (a tet has 4
    // faces; surface cells have fewer; +1 for the diagonal).
    PetscInt d_nz = 8;
    PetscInt o_nz = 4;
    MatCreate(comm_, &A_);
    MatSetSizes(A_, n_local_, n_local_, n_global_, n_global_);
    MatSetType(A_, MATAIJ);
    MatMPIAIJSetPreallocation(A_, d_nz, nullptr, o_nz, nullptr);
    MatSeqAIJSetPreallocation(A_, d_nz, nullptr);
    MatSetOption(A_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

    VecCreate(comm_, &rhs_);
    VecSetSizes(rhs_, n_local_, n_global_);
    VecSetType(rhs_, VECSTANDARD);
    VecDuplicate(rhs_, &sol_);

    KSPCreate(comm_, &ksp_);
    KSPSetOperators(ksp_, A_, A_);
    KSPSetType(ksp_, KSPGMRES);
    PC pc;
    KSPGetPC(ksp_, &pc);
    PetscMPIInt size = 1;
    MPI_Comm_size(comm_, &size);
    if (size > 1) {
        // Pass-12 parallel KSP convention.
        PCSetType(pc, PCBJACOBI);
    } else {
        PCSetType(pc, PCLU);
    }
    KSPSetTolerances(ksp_, cfg_.ksp_rtol, cfg_.ksp_atol, PETSC_DEFAULT,
                     cfg_.ksp_max_it);
    KSPSetFromOptions(ksp_);

    initialized_ = true;
}

void SourceBallRadiation3DSolver::assembleOperator(
    double dt,
    const std::vector<double>& rho,
    const std::vector<double>& T_iter,
    const std::vector<double>& T_m_old,
    const std::vector<double>& E_r_old)
{
    MatZeroEntries(A_);
    VecSet(rhs_, 0.0);

    // Volume term V/dt * E^{n+1} = V/dt * E^n.
    for (int i = 0; i < n_local_; ++i) {
        const double V = cells_[i].volume;
        const PetscInt row = local_offset_global_ + i;
        const double diag = V / dt;
        const double rhs_val = (V / dt) * E_r_old[i];
        MatSetValue(A_, row, row, diag, ADD_VALUES);
        VecSetValue(rhs_, row, rhs_val, ADD_VALUES);
    }

    // Internal face flux contributions (two-point flux approximation).
    for (const auto& f : internal_faces_) {
        const int iL = f.cell_left_local;
        // Compute D_face = c / (3 * kappa_R * rho_face). rho_face
        // harmonic-mean of left/right cell densities; same for kappa_R.
        double rhoL = (iL >= 0 && iL < n_local_) ? rho[iL] : 0.0;
        double TL   = (iL >= 0 && iL < n_local_) ? T_iter[iL] : 300.0;
        double rhoR = rhoL;  // ghost: replicate (one-sided)
        double TR   = TL;
        const int iR = f.cell_right_local;
        if (iR >= 0 && iR < n_local_) {
            rhoR = rho[iR];
            TR   = T_iter[iR];
        }
        const double rho_face = harmonicMean(rhoL, rhoR);
        if (rho_face <= 0.0) continue;
        const double kR_L = rosselandKappa(rhoL, TL);
        const double kR_R = rosselandKappa(rhoR, TR);
        const double kR_face = harmonicMean(kR_L, kR_R);
        if (kR_face <= 0.0) continue;
        const double D_face = C_LIGHT / (3.0 * kR_face * rho_face);
        const double tau = D_face * f.area / std::max(1.0e-12, f.dist_centroid);

        // Contribution to L row: A[L,L] += tau, A[L,R] -= tau.
        if (iL >= 0 && iL < n_local_) {
            MatSetValue(A_, f.cell_left_global, f.cell_left_global,  tau, ADD_VALUES);
            MatSetValue(A_, f.cell_left_global, f.cell_right_global, -tau, ADD_VALUES);
        }
        // Contribution to R row: A[R,R] += tau, A[R,L] -= tau.
        // Only stamp for owned R cells (otherwise the off-rank owner of
        // R will stamp the symmetric entry from its own face traversal).
        if (iR >= 0 && iR < n_local_) {
            MatSetValue(A_, f.cell_right_global, f.cell_right_global,  tau, ADD_VALUES);
            MatSetValue(A_, f.cell_right_global, f.cell_left_global,  -tau, ADD_VALUES);
        }
    }

    // Boundary faces.
    for (const auto& bf : boundary_faces_) {
        const int i = bf.cell_local;
        if (i < 0 || i >= n_local_) continue;
        const double rho_i = rho[i];
        if (rho_i <= 0.0) continue;
        const double kR = rosselandKappa(rho_i, T_iter[i]);
        if (kR <= 0.0) continue;
        const double D_face = C_LIGHT / (3.0 * kR * rho_i);
        const double tau = D_face * bf.area / std::max(1.0e-12, bf.dist_to_face);

        if (bf.type == BoundaryFaceType3D::ELASTIC) {
            // Vacuum Marshak: E_face = 0. One-sided flux equivalent to
            // a Dirichlet at the ghost: A[i,i] += tau, RHS unchanged
            // (target value = 0).
            MatSetValue(A_, bf.cell_global, bf.cell_global, tau, ADD_VALUES);
        } else {
            // CAVITY: prescribed-flux (default zero). A non-zero
            // prescribed flux would add to the RHS as +area * flux.
            // Pass-13b ships the no-flux baseline; pass-13c hooks up
            // the yield-deposition source.
        }
    }

    // Matter coupling. Per cell: source = V c kappa_P rho (a T_m^4 - E_r).
    // Linearise about T_iter:
    //   a T_m^4 ~ a T_iter^4 + 4 a T_iter^3 (T_m - T_iter)
    // We treat E_r implicitly and T_m_new explicitly through one outer
    // iteration. Diagonal contribution: + V c kappa_P rho.
    // RHS contribution:                   + V c kappa_P rho a T_iter^4
    //                                       (linearised at T_iter).
    // T_m feedback enters via the next Newton iterate.
    for (int i = 0; i < n_local_; ++i) {
        const double V = cells_[i].volume;
        const double rho_i = rho[i];
        if (rho_i <= 0.0) continue;
        const double kP = planckKappa(rho_i, T_iter[i]);
        const double coef = V * C_LIGHT * kP * rho_i;
        const PetscInt row = local_offset_global_ + i;
        MatSetValue(A_, row, row, coef, ADD_VALUES);
        VecSetValue(rhs_, row, coef * A_RAD * std::pow(T_iter[i], 4.0),
                    ADD_VALUES);
        // Suppress unused warnings.
        (void)T_m_old;
    }

    MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(rhs_);
    VecAssemblyEnd(rhs_);
}

void SourceBallRadiation3DSolver::updateMatterTemperature(
    double dt,
    const std::vector<double>& rho,
    const std::vector<double>& T_iter,
    const std::vector<double>& E_r_new,
    std::vector<double>& T_m,
    std::vector<double>& e_int)
{
    // Per-cell ODE: rho cv dT/dt = c kappa_P rho (E_r - a T^4).
    // Linearised about T_iter:
    //   T_new = T_iter + dt c kappa_P (E_r - a T_iter^4) / cv
    //   subject to keeping the explicit step stable; if the step would
    //   reduce T below ambient too quickly we clamp to a positive
    //   non-zero floor (pass-13b uses 1 K).
    constexpr double T_floor = 1.0;
    for (int i = 0; i < n_local_; ++i) {
        const double rho_i = rho[i];
        if (rho_i <= 0.0) continue;
        const double kP = planckKappa(rho_i, T_iter[i]);
        const double T_old = T_m[i];
        const double dE = (E_r_new[i] - A_RAD * std::pow(T_iter[i], 4.0));
        const double dT = dt * C_LIGHT * kP * dE / cfg_.cv_J_per_kg_K;
        double T_new = T_old + dT;
        if (T_new < T_floor) T_new = T_floor;
        T_m[i] = T_new;
        // e_int change [J/kg]: cv * (T_new - T_old).
        e_int[i] += cfg_.cv_J_per_kg_K * (T_new - T_old);
    }
}

SourceBallRadiation3DStepResult SourceBallRadiation3DSolver::step(
    double dt,
    const std::vector<double>& rho,
    std::vector<double>& T_m,
    std::vector<double>& E_r,
    std::vector<double>& e_int)
{
    SourceBallRadiation3DStepResult res;
    if (!initialized_) {
        throw std::runtime_error(
            "SourceBallRadiation3DSolver::step called before initialize().");
    }

    std::vector<double> T_iter = T_m;
    std::vector<double> T_m_old = T_m;
    std::vector<double> E_r_old = E_r;
    std::vector<double> e_int_in = e_int;
    std::vector<double> E_r_new(n_local_, 0.0);

    int it = 0;
    double resid_inf = 0.0;
    for (it = 0; it < cfg_.max_newton_iter; ++it) {
        assembleOperator(dt, rho, T_iter, T_m_old, E_r_old);

        KSPSolve(ksp_, rhs_, sol_);
        KSPConvergedReason reason;
        KSPGetConvergedReason(ksp_, &reason);
        if (reason < 0) {
            // KSP failed; record and bail to caller.
            res.converged = false;
            res.newton_iters = it + 1;
            return res;
        }

        // Pull local solution into E_r_new.
        const PetscScalar* sa = nullptr;
        VecGetArrayRead(sol_, &sa);
        for (int i = 0; i < n_local_; ++i) {
            E_r_new[i] = static_cast<double>(sa[i]);
        }
        VecRestoreArrayRead(sol_, &sa);

        // Update T_m / e_int from per-cell energy balance.
        T_m = T_m_old;
        e_int = e_int_in;
        updateMatterTemperature(dt, rho, T_iter, E_r_new, T_m, e_int);

        // Convergence: max |T_iter - T_m|.
        double dmax = 0.0;
        for (int i = 0; i < n_local_; ++i) {
            const double dT = std::abs(T_m[i] - T_iter[i]);
            if (dT > dmax) dmax = dT;
        }
        resid_inf = dmax;
        T_iter = T_m;
        if (dmax < cfg_.newton_tolerance) {
            ++it;
            break;
        }
    }

    E_r = E_r_new;

    // Diagnostics.
    double localE = 0.0;
    for (int i = 0; i < n_local_; ++i) {
        localE += E_r[i] * cells_[i].volume;
    }
    double localdMatter = 0.0;
    for (int i = 0; i < n_local_; ++i) {
        const double dE = e_int[i] - e_int_in[i];
        localdMatter += dE * rho[i] * cells_[i].volume;
    }
    res.converged = (resid_inf < cfg_.newton_tolerance);
    res.newton_iters = it;
    res.residual_inf_norm = resid_inf;
    res.local_radiation_energy_J = localE;
    res.local_matter_energy_change_J = localdMatter;
    return res;
}

}  // namespace FSRM
