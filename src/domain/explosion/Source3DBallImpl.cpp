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

constexpr const char* kPhysicsName = "Source3DBallImpl_v1_pass13c_validation";

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

    // Pass-13c: snapshot the IC stress as the reference state for the
    // surface-integral stress-drop moment tensor (Day & McLaughlin
    // 1991). Done after applyOverburdenIC so the asymmetric overburden
    // initial state does not contribute to the moment-tensor integral
    // (it's the response to the source, not the static field, that
    // generates seismic radiation).
    sigma_initial_ = sigma_;
    M_global_ = {0, 0, 0, 0, 0, 0};
    Mdot_global_ = {0, 0, 0, 0, 0, 0};
    M_global_prev_ = {0, 0, 0, 0, 0, 0};

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
    elastic_surface_faces_.clear();
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

            // Pass-13c surface-integral cache. Only ELASTIC-classified
            // boundary faces participate in the moment-tensor integral.
            if (bf.type == BoundaryFaceType3D::ELASTIC) {
                ElasticSurfaceFace ef;
                ef.cell_local = iL;
                ef.area = area;
                ef.centroid = {fcentroid[0], fcentroid[1], fcentroid[2]};
                // Sign-correct the normal so it points away from the
                // owning cell's centroid (outward at the elastic
                // surface). DMPlex's face normal orientation is up to
                // sign; the dot product with (face_centroid - cell_
                // centroid) disambiguates.
                double nx = fnormal[0];
                double ny = fnormal[1];
                double nz = fnormal[2];
                const double nmag = std::sqrt(nx * nx + ny * ny + nz * nz);
                if (nmag > 1.0e-30) {
                    nx /= nmag;
                    ny /= nmag;
                    nz /= nmag;
                }
                const double cx = fcentroid[0] - pL[0];
                const double cy = fcentroid[1] - pL[1];
                const double cz = fcentroid[2] - pL[2];
                const double dotp = nx * cx + ny * cy + nz * cz;
                if (dotp < 0.0) {
                    nx = -nx;
                    ny = -ny;
                    nz = -nz;
                }
                ef.outward_normal = {nx, ny, nz};
                elastic_surface_faces_.push_back(ef);
            }
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

void Source3DBallImpl::setCellStressForTest(
    const std::vector<std::array<double, 6>>& sig)
{
    if (static_cast<int>(sig.size()) != n_local_cells_) {
        throw std::runtime_error(
            "Source3DBallImpl::setCellStressForTest size mismatch: got "
            + std::to_string(sig.size()) + " entries for "
            + std::to_string(n_local_cells_) + " cells.");
    }
    sigma_ = sig;
}

void Source3DBallImpl::setCellStressInitialForTest(
    const std::vector<std::array<double, 6>>& sig0)
{
    if (static_cast<int>(sig0.size()) != n_local_cells_) {
        throw std::runtime_error(
            "Source3DBallImpl::setCellStressInitialForTest size mismatch: "
            "got " + std::to_string(sig0.size()) + " entries for "
            + std::to_string(n_local_cells_) + " cells.");
    }
    sigma_initial_ = sig0;
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

    // Pass-13c: integrate the cell-centred stress drop over the
    // elastic-radius surface to obtain M(t); finite-difference for
    // Mdot(t).
    recomputeMomentTensorFromSurface(dt);

    current_time_ += dt;
    state_snapshot_.time = current_time_;
    state_snapshot_.moment_tensor_N_m = M_global_;
    state_snapshot_.moment_rate_tensor_N_m_per_s = Mdot_global_;
}

void Source3DBallImpl::recomputeMomentTensorFromSurface(double dt)
{
    // Local accumulator: M_ij_local = sum_face (sigma_jk - sigma_jk_0)
    //                                  * n_k * x_i * A
    // x is the face-centroid offset from origin (= source center in
    // the local mesh frame). On a spherical surface of radius R this
    // reduces to R * (sigma' n) (x) tensor as in Aki & Richards 2002
    // ch 4. The expression is valid for arbitrary surfaces.
    std::array<double, 6> M_local = {0, 0, 0, 0, 0, 0};
    for (const auto& f : elastic_surface_faces_) {
        const int ic = f.cell_local;
        if (ic < 0 || ic >= static_cast<int>(sigma_.size())) continue;
        const auto& sig = sigma_[ic];
        const auto& sig0 = (ic < static_cast<int>(sigma_initial_.size()))
                           ? sigma_initial_[ic]
                           : std::array<double, 6>{0, 0, 0, 0, 0, 0};
        const double dxx = sig[0] - sig0[0];
        const double dyy = sig[1] - sig0[1];
        const double dzz = sig[2] - sig0[2];
        const double dxy = sig[3] - sig0[3];
        const double dxz = sig[4] - sig0[4];
        const double dyz = sig[5] - sig0[5];
        // traction t_j = sigma_jk * n_k (Voigt: sigma is symmetric, so
        // t_x = sxx*nx + sxy*ny + sxz*nz; etc.)
        const double nx = f.outward_normal[0];
        const double ny = f.outward_normal[1];
        const double nz = f.outward_normal[2];
        const double tx = dxx * nx + dxy * ny + dxz * nz;
        const double ty = dxy * nx + dyy * ny + dyz * nz;
        const double tz = dxz * nx + dyz * ny + dzz * nz;
        const double xi = f.centroid[0];
        const double yi = f.centroid[1];
        const double zi = f.centroid[2];
        const double A = f.area;
        // M_ij = t_j * x_i * A (i row, j column). Voigt order: xx, yy,
        // zz, xy, xz, yz. By symmetry M is symmetric, so we accumulate
        // the symmetrised form (1/2)(t_j x_i + t_i x_j) for the
        // off-diagonals.
        M_local[0] += xi * tx * A;
        M_local[1] += yi * ty * A;
        M_local[2] += zi * tz * A;
        M_local[3] += 0.5 * (xi * ty + yi * tx) * A;
        M_local[4] += 0.5 * (xi * tz + zi * tx) * A;
        M_local[5] += 0.5 * (yi * tz + zi * ty) * A;
    }
    std::array<double, 6> M_new = {0, 0, 0, 0, 0, 0};
    MPI_Allreduce(M_local.data(), M_new.data(), 6, MPI_DOUBLE, MPI_SUM, comm_);

    if (dt > 0.0) {
        for (int k = 0; k < 6; ++k) {
            Mdot_global_[k] = (M_new[k] - M_global_[k]) / dt;
        }
    } else {
        Mdot_global_ = {0, 0, 0, 0, 0, 0};
    }
    M_global_prev_ = M_global_;
    M_global_ = M_new;
}

void Source3DBallImpl::cavityRadiusExtremes(double& r_min_out,
                                            double& r_max_out) const
{
    double r_min_local = std::numeric_limits<double>::infinity();
    double r_max_local = 0.0;
    if (mesh_) {
        DM dm = mesh_->getDM();
        DMLabel marker = mesh_->getVertexMarkerLabel();
        if (dm && marker) {
            PetscInt vStart = 0, vEnd = 0;
            DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd);
            Vec coordVec = nullptr;
            DMGetCoordinatesLocal(dm, &coordVec);
            const PetscScalar* coords = nullptr;
            if (coordVec) VecGetArrayRead(coordVec, &coords);
            PetscSection coordSection = nullptr;
            DMGetCoordinateSection(dm, &coordSection);
            for (PetscInt v = vStart; v < vEnd; ++v) {
                PetscInt mv = 0;
                DMLabelGetValue(marker, v, &mv);
                if (mv != SourceBallVertexMarkerValues::INNER_CAVITY_SURFACE) {
                    continue;
                }
                PetscInt off = 0;
                PetscSectionGetOffset(coordSection, v, &off);
                const double x = static_cast<double>(coords[off + 0]);
                const double y = static_cast<double>(coords[off + 1]);
                const double z = static_cast<double>(coords[off + 2]);
                const double r = std::sqrt(x * x + y * y + z * z);
                r_min_local = std::min(r_min_local, r);
                r_max_local = std::max(r_max_local, r);
            }
            if (coordVec) VecRestoreArrayRead(coordVec, &coords);
        }
    }
    if (!std::isfinite(r_min_local)) r_min_local = 0.0;
    double r_min_global = r_min_local;
    double r_max_global = r_max_local;
    MPI_Allreduce(&r_min_local, &r_min_global, 1, MPI_DOUBLE, MPI_MIN, comm_);
    MPI_Allreduce(&r_max_local, &r_max_global, 1, MPI_DOUBLE, MPI_MAX, comm_);
    r_min_out = r_min_global;
    r_max_out = r_max_global;
}

void Source3DBallImpl::getMomentTensor(std::array<double, 6>& M) const
{
    // Pass-13c: M_ij(t) is the surface integral of the stress drop
    // (sigma(t) - sigma(0)) crossed with the centroid lever arm, taken
    // over the elastic-radius spherical boundary. See
    // recomputeMomentTensorFromSurface() for the integral form.
    M = M_global_;
}

void Source3DBallImpl::getMomentRateTensor(std::array<double, 6>& Mdot) const
{
    Mdot = Mdot_global_;
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
