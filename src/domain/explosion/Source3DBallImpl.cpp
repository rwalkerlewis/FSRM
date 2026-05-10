/**
 * @file Source3DBallImpl.cpp
 * @brief Pass-13b (axis-1b physics) Source3DBallImpl implementation.
 *        See Source3DBallImpl.hpp for the contract and the explicit
 *        pass-13a / pass-13b / pass-13c slice boundaries.
 */

#include "domain/explosion/Source3DBallImpl.hpp"

#include "domain/explosion/Source3DBallMesh.hpp"

#include <petscdmplex.h>
#include <petscdmlabel.h>
#include <petscsys.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace FSRM {

namespace
{

constexpr const char* kPhysicsName = "Source3DBallImpl_v1_pass13b_physics";

/// True if all vertices in a face's transitive closure carry the
/// requested marker value. Used to classify boundary faces as
/// CAVITY (marker = 2) vs ELASTIC (marker = 1).
bool faceHasUniformVertexMarker(DM dm, DMLabel marker_label, PetscInt face,
                                int target_value)
{
    if (!marker_label) return false;
    PetscInt closure_size = 0;
    PetscInt* closure = nullptr;
    DMPlexGetTransitiveClosure(dm, face, PETSC_TRUE, &closure_size, &closure);
    if (!closure) return false;
    PetscInt vStart = 0;
    PetscInt vEnd = 0;
    DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd);
    bool all_match = true;
    bool any_vertex = false;
    for (PetscInt i = 0; i < closure_size; ++i) {
        const PetscInt p = closure[2 * i];
        if (p < vStart || p >= vEnd) continue;
        any_vertex = true;
        PetscInt v = 0;
        DMLabelGetValue(marker_label, p, &v);
        if (v != target_value) { all_match = false; break; }
    }
    DMPlexRestoreTransitiveClosure(dm, face, PETSC_TRUE, &closure_size, &closure);
    return any_vertex && all_match;
}

}  // namespace

Source3DBallImpl::Source3DBallImpl() = default;
Source3DBallImpl::~Source3DBallImpl() = default;

void Source3DBallImpl::setComm(MPI_Comm comm)
{
    comm_ = comm;
    comm_set_ = true;
}

void Source3DBallImpl::initialize(const Source3DBallConfig& cfg)
{
    cfg_ = cfg;

    if (cfg.radiation_discretization
        == Source3DBallConfig::RadiationDiscretization::FEM_NODAL) {
        throw std::runtime_error(
            "Source3DBallImpl: radiation_discretization=FEM_NODAL is the "
            "axis-1b nodal-FEM radiation block: pass-15+ work, not yet "
            "implemented. Use FV_CELL_CENTRED for the pass-13b grey "
            "diffusion path.");
    }

    if (!comm_set_) {
        comm_ = PETSC_COMM_WORLD;
    }

    // Cache material constants used by the constitutive + radiation
    // updates. All come from cfg with safe granite-class defaults.
    bulk_modulus_K_ = cfg.bulk_modulus_K_pa;
    shear_modulus_G_ = cfg.shear_modulus_G_pa;
    dp_params_.alpha_dp = cfg.dp3d_alpha;
    dp_params_.k_dp = cfg.dp3d_k_pa;
    dp_params_.medium_label = cfg.medium_label;

    if (cfg.mesh_path.empty()) {
        // Foundation no-op mode: caller wants to verify the factory
        // returns a valid object without touching the disk.
        mesh_loaded_ = false;
        n_local_cells_ = 0;
        sigma_.clear();
        strain_rate_.clear();
        eps_p_eq_.clear();
        rho_.clear();
        e_int_.clear();
        T_m_.clear();
        E_r_.clear();
        rad_cells_.clear();
        rad_internal_faces_.clear();
        rad_boundary_faces_.clear();
        state_snapshot_ = Source3DBallState{};
        return;
    }

    mesh_ = std::make_unique<Source3DBallMesh>();
    mesh_->loadFromTetGen(comm_, cfg.mesh_path);
    mesh_loaded_ = true;

    // Pass-13a sanity-check the loaded mesh against the configured
    // cavity / elastic radii.
    double local_min = mesh_->localMinVertexRadius();
    double local_max = mesh_->localMaxVertexRadius();
    if (mesh_->numLocalVertices() == 0) {
        local_min = std::numeric_limits<double>::infinity();
        local_max = 0.0;
    }
    double global_min = local_min;
    double global_max = local_max;
    MPI_Allreduce(&local_min, &global_min, 1, MPI_DOUBLE, MPI_MIN, comm_);
    MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, comm_);

    if (cfg.cavity_radius_m > 0.0) {
        const double rel = std::abs(global_min - cfg.cavity_radius_m)
                           / cfg.cavity_radius_m;
        if (rel > 0.20) {
            throw std::runtime_error(
                "Source3DBallImpl: configured cavity_radius_m="
                + std::to_string(cfg.cavity_radius_m)
                + " m disagrees with mesh min radius="
                + std::to_string(global_min)
                + " m by >20 percent. Verify cfg.mesh_path.");
        }
    }
    if (cfg.outer_radius_m > 0.0) {
        const double rel = std::abs(global_max - cfg.outer_radius_m)
                           / cfg.outer_radius_m;
        if (rel > 0.20) {
            throw std::runtime_error(
                "Source3DBallImpl: configured outer_radius_m="
                + std::to_string(cfg.outer_radius_m)
                + " m disagrees with mesh max radius="
                + std::to_string(global_max)
                + " m by >20 percent. Verify cfg.mesh_path.");
        }
    }

    buildCellAndFaceLists();

    // Allocate per-cell state vectors.
    sigma_.assign(n_local_cells_, {0, 0, 0, 0, 0, 0});
    strain_rate_.assign(n_local_cells_, {0, 0, 0, 0, 0, 0});
    eps_p_eq_.assign(n_local_cells_, 0.0);
    rho_.assign(n_local_cells_, cfg.host_density_kg_per_m3);
    e_int_.assign(n_local_cells_, 0.0);
    T_m_.assign(n_local_cells_, cfg.T_ambient_K);
    E_r_.assign(n_local_cells_, 0.0);

    // Apply asymmetric overburden initial state (when enabled).
    if (cfg.asymmetric_overburden) {
        applyOverburdenIC();
    }

    // Build the radiation solver on the same DM topology.
    rad_solver_ = std::make_unique<SourceBallRadiation3DSolver>();
    SourceBallRadiation3DSolver::Config rcfg;
    rcfg.cv_J_per_kg_K = cfg.cv_J_per_kg_K;
    rcfg.kappa_constant_m2_per_kg = 0.0;  // power-law default
    rcfg.opacity_params = PowerLawOpacitySets::byName(cfg.medium_label);
    rad_solver_->setConfig(rcfg);

    // Compute global cell offset: the natural-numbering convention is
    // that ranks own contiguous global ids. Use MPI_Scan on
    // n_local_cells_ to find this rank's first global id.
    PetscInt local_count = n_local_cells_;
    PetscInt scan_count = 0;
    MPI_Scan(&local_count, &scan_count, 1, MPIU_INT, MPI_SUM, comm_);
    PetscInt local_offset = scan_count - local_count;
    PetscInt n_global = mesh_->numGlobalCells();

    // Patch the rad face/cell lists with global ids (the caller-side
    // builder used local-only ids). For pass-13b at MPI=1 this is a
    // simple offset; multi-rank consistency comes via the MPI_Scan
    // above. Faces with off-rank neighbours retain
    // cell_right_local = -1 and have to skip the off-rank stamp; the
    // off-rank owner would emit the symmetric entry. Pass-13b ships
    // the single-rank Mat assembly path (the multi-rank face exchange
    // is pass-13c work; the pass-12 parallel KSP convention still
    // applies to the rank-local Mat).
    for (auto& f : rad_internal_faces_) {
        f.cell_left_global = local_offset + f.cell_left_local;
        f.cell_right_global = (f.cell_right_local >= 0)
                              ? local_offset + f.cell_right_local
                              : -1;
    }
    for (auto& bf : rad_boundary_faces_) {
        bf.cell_global = local_offset + bf.cell_local;
    }

    rad_solver_->initialize(comm_, n_local_cells_, local_offset, n_global,
                            rad_cells_, rad_internal_faces_, rad_boundary_faces_);

    state_snapshot_.time = 0.0;
    state_snapshot_.cavity_radius_max_m = global_min;
    state_snapshot_.cavity_radius_min_m = global_min;
    state_snapshot_.moment_tensor_N_m = {0, 0, 0, 0, 0, 0};
    state_snapshot_.moment_rate_tensor_N_m_per_s = {0, 0, 0, 0, 0, 0};
}

void Source3DBallImpl::buildCellAndFaceLists()
{
    DM dm = mesh_->getDM();
    if (!dm) {
        n_local_cells_ = 0;
        return;
    }

    PetscInt cStart = 0, cEnd = 0;
    DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd);
    PetscInt fStart = 0, fEnd = 0;
    DMPlexGetHeightStratum(dm, 1, &fStart, &fEnd);
    n_local_cells_ = static_cast<int>(cEnd - cStart);

    // Cell geometry: volume + centroid via DMPlexComputeCellGeometryFVM.
    centroid_.assign(n_local_cells_, {0, 0, 0});
    volume_.assign(n_local_cells_, 0.0);
    rad_cells_.clear();
    rad_cells_.reserve(n_local_cells_);
    for (PetscInt c = cStart; c < cEnd; ++c) {
        PetscReal vol = 0.0;
        PetscReal centroid[3] = {0, 0, 0};
        DMPlexComputeCellGeometryFVM(dm, c, &vol, centroid, nullptr);
        const int li = static_cast<int>(c - cStart);
        volume_[li] = vol;
        centroid_[li] = {centroid[0], centroid[1], centroid[2]};
        CellGeometry3D cg;
        cg.volume = vol;
        cg.centroid = {centroid[0], centroid[1], centroid[2]};
        rad_cells_.push_back(cg);
    }

    // Face traversal.
    rad_internal_faces_.clear();
    rad_boundary_faces_.clear();
    DMLabel marker = mesh_->getVertexMarkerLabel();

    for (PetscInt f = fStart; f < fEnd; ++f) {
        PetscReal area = 0.0;
        PetscReal fcentroid[3] = {0, 0, 0};
        PetscReal fnormal[3] = {0, 0, 0};
        DMPlexComputeCellGeometryFVM(dm, f, &area, fcentroid, fnormal);
        if (area <= 0.0) continue;

        PetscInt support_size = 0;
        const PetscInt* support = nullptr;
        DMPlexGetSupportSize(dm, f, &support_size);
        DMPlexGetSupport(dm, f, &support);
        if (support_size == 0) continue;

        if (support_size == 2) {
            const int iL = static_cast<int>(support[0] - cStart);
            const int iR = static_cast<int>(support[1] - cStart);
            if (iL < 0 || iL >= n_local_cells_) continue;
            if (iR < 0 || iR >= n_local_cells_) continue;
            const auto& pL = centroid_[iL];
            const auto& pR = centroid_[iR];
            const double dx = pR[0] - pL[0];
            const double dy = pR[1] - pL[1];
            const double dz = pR[2] - pL[2];
            const double dist = std::sqrt(dx * dx + dy * dy + dz * dz);
            InternalFace3D ifc;
            ifc.cell_left_local = iL;
            ifc.cell_right_local = iR;
            ifc.cell_left_global = iL;   // patched after MPI_Scan
            ifc.cell_right_global = iR;
            ifc.area = area;
            ifc.dist_centroid = (dist > 1.0e-12) ? dist : 1.0e-12;
            rad_internal_faces_.push_back(ifc);
        } else if (support_size == 1) {
            const int iL = static_cast<int>(support[0] - cStart);
            if (iL < 0 || iL >= n_local_cells_) continue;
            const auto& pL = centroid_[iL];
            const double dx = fcentroid[0] - pL[0];
            const double dy = fcentroid[1] - pL[1];
            const double dz = fcentroid[2] - pL[2];
            const double dist = std::sqrt(dx * dx + dy * dy + dz * dz);
            BoundaryFace3D bf;
            bf.cell_local = iL;
            bf.cell_global = iL;  // patched after MPI_Scan
            bf.area = area;
            bf.dist_to_face = (dist > 1.0e-12) ? dist : 1.0e-12;
            // Classify boundary face by vertex marker.
            const bool is_cavity = faceHasUniformVertexMarker(
                dm, marker, f, SourceBallVertexMarkerValues::INNER_CAVITY_SURFACE);
            const bool is_elastic = faceHasUniformVertexMarker(
                dm, marker, f, SourceBallVertexMarkerValues::OUTER_ELASTIC_SURFACE);
            if (is_cavity) {
                bf.type = BoundaryFaceType3D::CAVITY;
            } else if (is_elastic) {
                bf.type = BoundaryFaceType3D::ELASTIC;
            } else {
                // Faces without uniform marker (e.g. mixed-marker faces
                // arising from the icosphere tessellation) default to
                // ELASTIC. Pass-13b conservative default; pass-13c
                // refines this with a proper face-marker label.
                bf.type = BoundaryFaceType3D::ELASTIC;
            }
            rad_boundary_faces_.push_back(bf);
        }
    }
}

void Source3DBallImpl::applyOverburdenIC()
{
    SourceBallOverburdenConfig oc;
    oc.source_depth_m = cfg_.source_depth_m;
    oc.rho_solid = cfg_.host_density_kg_per_m3;
    oc.g = 9.81;
    oc.K_0 = cfg_.overburden_K0;
    for (int i = 0; i < n_local_cells_; ++i) {
        sigma_[i] = overburdenStressAtCell(centroid_[i], oc);
    }
}

void Source3DBallImpl::setCellStrainRate(
    const std::vector<std::array<double, 6>>& rates)
{
    if (static_cast<int>(rates.size()) != n_local_cells_) {
        throw std::runtime_error(
            "Source3DBallImpl::setCellStrainRate size mismatch: got "
            + std::to_string(rates.size()) + " entries for "
            + std::to_string(n_local_cells_) + " cells.");
    }
    strain_rate_ = rates;
}

void Source3DBallImpl::setCellDensity(const std::vector<double>& rho)
{
    if (static_cast<int>(rho.size()) != n_local_cells_) {
        throw std::runtime_error(
            "Source3DBallImpl::setCellDensity size mismatch: got "
            + std::to_string(rho.size()) + " entries for "
            + std::to_string(n_local_cells_) + " cells.");
    }
    rho_ = rho;
}

void Source3DBallImpl::step(double dt)
{
    if (!mesh_loaded_) {
        throw std::runtime_error(
            "Source3DBallImpl::step: foundation no-op mode (mesh_path "
            "was empty at initialize). Set cfg.mesh_path to a TetGen "
            "mesh basename before calling step.");
    }

    last_yielded_count_ = 0;
    double total_plastic_dissipation = 0.0;
    for (int i = 0; i < n_local_cells_; ++i) {
        std::array<double, 6> deps;
        for (int k = 0; k < 6; ++k) deps[k] = strain_rate_[i][k] * dt;
        std::array<double, 6> sigma_new{};
        std::array<double, 6> deps_p{};
        auto res = druckerPrager3DRadialReturn(sigma_[i], deps,
                                                bulk_modulus_K_, shear_modulus_G_,
                                                dp_params_, sigma_new, deps_p);
        sigma_[i] = sigma_new;
        if (res.yielded) {
            ++last_yielded_count_;
            eps_p_eq_[i] += res.delta_eps_p_eq;
            // Plastic work per unit mass = sqrt(J2) * delta_eps_p_eq /
            // rho. Use the Y-equivalent at the projection point as a
            // conservative estimate. Adds to internal energy.
            const double Y_proj = dp_params_.alpha_dp
                * (-(sigma_new[0] + sigma_new[1] + sigma_new[2]) / 3.0)
                + dp_params_.k_dp;
            const double dW = Y_proj * res.delta_eps_p_eq;
            const double rho_safe = (rho_[i] > 1.0e-9) ? rho_[i] : 1.0e-9;
            e_int_[i] += dW / rho_safe;
            total_plastic_dissipation += dW * volume_[i];
        }
    }

    if (rad_solver_ && rad_solver_->isInitialized()) {
        auto rres = rad_solver_->step(dt, rho_, T_m_, E_r_, e_int_);
        last_rad_newton_iters_ = rres.newton_iters;
    }

    current_time_ += dt;
    state_snapshot_.time = current_time_;
}

void Source3DBallImpl::getMomentTensor(std::array<double, 6>& M) const
{
    // Pass-13b: surface-integral moment-tensor extraction is pass-13c
    // work. The 3D solver runs the constitutive / radiation update on
    // every owned cell but does not yet integrate stresses across the
    // elastic-radius surface. M remains identically zero so the host
    // gets a consistent (zero) signal that 3D mode is selected but
    // moment-tensor extraction is not yet wired.
    M = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
}

void Source3DBallImpl::getMomentRateTensor(std::array<double, 6>& Mdot) const
{
    Mdot = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
}

void Source3DBallImpl::getState(Source3DBallState& state) const
{
    state = state_snapshot_;
}

const char* Source3DBallImpl::name() const
{
    return kPhysicsName;
}

void Source3DBallImpl::getCellStates(std::vector<Source3DCellState>& out) const
{
    out.clear();
    out.reserve(n_local_cells_);
    for (int i = 0; i < n_local_cells_; ++i) {
        Source3DCellState s;
        s.centroid = centroid_[i];
        s.sigma = sigma_[i];
        s.eps_p_eq = eps_p_eq_[i];
        s.rho = rho_[i];
        s.e_int = e_int_[i];
        s.T_m = T_m_[i];
        s.E_r = E_r_[i];
        out.push_back(s);
    }
}

}  // namespace FSRM
