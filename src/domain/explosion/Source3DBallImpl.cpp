/**
 * @file Source3DBallImpl.cpp
 * @brief Pass-13b/13c/14a (axis-1b) Source3DBallImpl implementation.
 *        See Source3DBallImpl.hpp for the contract and the explicit
 *        pass-13a / pass-13b / pass-13c / pass-14a slice boundaries.
 *
 * Pass-14a (axis-1b source-forcing slice) adds the internal source term
 * the 3D path was missing: volumetric mechanical-energy deposition into
 * the inner-cavity cells (the cavity wall), a Tillotson-EOS pressure
 * rise from that energy, and an explicit finite-volume acoustic pulse
 * that carries the pressure perturbation out to the elastic-radius
 * extraction surface, where the pass-13c surface integral picks it up
 * as a real moment-rate tensor. The cavity does not grow (no mesh
 * motion / Lagrangian advection; that is pass-14b); the stress and
 * plastic-strain fields evolve in place exactly as the pass-13b
 * constitutive update already did.
 *
 * References (pass-14a):
 *   - Mueller, G. and Murphy, J. R. (1971), "Seismic characteristics of
 *     underground nuclear detonations", BSSA 61(6), 1675-1692
 *     (reduced-displacement-potential source-time function).
 *   - Brune, J. N. (1970), "Tectonic stress and the spectra of seismic
 *     shear waves from earthquakes", JGR 75, 4997-5009 (omega-square
 *     source pulse shape).
 *   - Denny, M. D. and Johnson, L. R. (1991), "The explosion seismic
 *     source function", in Explosion Source Phenomenology (AGU Monogr.
 *     65), 1-24 (energy-balance picture of the cavity source).
 *   - Day, S. M. and McLaughlin, K. L. (1991), "Effects of plasticity on
 *     the seismic source of underground explosions", JGR 96, 1955-1975
 *     (surface-integral moment-tensor extraction, retained from 13c).
 *   - LeVeque, R. J. (2002), Finite Volume Methods for Hyperbolic
 *     Problems, ch 3 (explicit FV acoustics on unstructured cells).
 */

#include "domain/explosion/Source3DBallImpl.hpp"

#include "domain/explosion/Source3DBallMesh.hpp"

#include <petscdmplex.h>
#include <petscdmlabel.h>
#include <petscsys.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>

namespace FSRM {

namespace
{

constexpr const char* kPhysicsName = "Source3DBallImpl_v1_pass14a_source_forcing";

/// J per kiloton of TNT (the standard 1e9 cal x 4.184 J/cal scaling).
constexpr double kJoulesPerKt = 4.184e12;

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

/// Normalised source-time-function rate (the time integral over [0, inf)
/// is unity for MUELLER_MURPHY / BRUNE; over [0, duration] it is unity
/// for RAMP; over [0, 0.02*duration] it is unity for DELTA). The caller
/// multiplies by the total deposited energy to get J/s.
///
/// MUELLER_MURPHY / BRUNE both use the critically-damped reduced-
/// displacement-potential pulse g(t) = omega^2 t exp(-omega t), the
/// Fourier pair of the Mueller-Murphy (1971) RDP spectrum (also the
/// Brune 1970 omega-square moment-rate). They differ only in how the
/// shape timescale relates to the configured deposition duration:
/// MUELLER_MURPHY puts the peak at duration/4 (omega = 4/duration);
/// BRUNE is faster, peak at duration/8 (omega = 8/duration).
double stfShapeRate(Source3DBallConfig::SourceTimeFunction kind,
                    double t, double duration)
{
    const double T = std::max(1.0e-12, duration);
    if (t < 0.0) return 0.0;
    switch (kind) {
        case Source3DBallConfig::SourceTimeFunction::MUELLER_MURPHY: {
            const double w = 4.0 / T;
            return w * w * t * std::exp(-w * t);
        }
        case Source3DBallConfig::SourceTimeFunction::BRUNE: {
            const double w = 8.0 / T;
            return w * w * t * std::exp(-w * t);
        }
        case Source3DBallConfig::SourceTimeFunction::RAMP: {
            // Boxcar: cumulative deposited energy ramps linearly to the
            // total over [0, T].
            return (t <= T) ? (1.0 / T) : 0.0;
        }
        case Source3DBallConfig::SourceTimeFunction::DELTA: {
            // Approximate delta: all energy concentrated in the first
            // 2 percent of the deposition window. A true delta is not
            // representable in a finite-dt stepper.
            const double Td = 0.02 * T;
            return (t <= Td) ? (1.0 / Td) : 0.0;
        }
    }
    return 0.0;
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

    // Reset pass-14a source-forcing state on every initialize().
    source_forcing_active_ = false;
    total_source_energy_J_ = 0.0;
    deposited_energy_global_J_ = 0.0;
    source_forcing_t0_ = 0.0;
    inner_cavity_cells_.clear();
    inner_cavity_volume_local_ = 0.0;
    inner_cavity_volume_global_ = 0.0;
    inner_cavity_radius_a_ = 0.0;
    dp_cavity_mean_ = 0.0;
    dp_cell_.clear();
    source_power_density_.clear();
    current_time_ = 0.0;

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

    // Pass-14a: identify the inner-cavity cells, build the DMLabel, set
    // up the acoustic state and the Tillotson EOS for the source pulse.
    setupSourceForcing();

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

// =========================================================================
// Pass-14a source forcing
// =========================================================================

void Source3DBallImpl::setupSourceForcing()
{
    // Always allocate the pressure-field / diagnostic vectors so
    // getCellStates and the radial-return coupling have well-defined
    // zero state even when source forcing is disabled.
    dp_cell_.assign(n_local_cells_, 0.0);
    source_power_density_.assign(n_local_cells_, 0.0);

    if (!cfg_.source_forcing_enabled) {
        return;
    }

    // Total deposited mechanical energy.
    double yield = (cfg_.source_yield_kt > 0.0) ? cfg_.source_yield_kt : 0.0;
    double eff = cfg_.source_deposition_efficiency;
    if (eff <= 0.0) eff = 1.0;
    if (eff > 1.0) eff = 1.0;
    total_source_energy_J_ = yield * kJoulesPerKt * eff;

    if (total_source_energy_J_ <= 0.0) {
        int rank = 0;
        MPI_Comm_rank(comm_, &rank);
        if (rank == 0) {
            std::fprintf(stderr,
                "Source3DBallImpl: source_forcing_enabled but yield is "
                "non-positive (%.3g kt); source forcing disabled.\n",
                cfg_.source_yield_kt);
        }
        return;
    }

    // Identify inner-cavity cells: centroid radius at most cavity_radius_m.
    // For a shell mesh the cavity is a hole, so no centroid lies inside
    // cavity_radius_m; in that case fall back to the innermost layer,
    // defined as cells whose centroid radius is within
    // (global-min-centroid-radius + 0.10 * cavity_radius_m). This always
    // catches the cavity-wall cells.
    auto centroid_radius = [&](int i) {
        const auto& c = centroid_[i];
        return std::sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
    };
    double local_min_centroid_r = std::numeric_limits<double>::infinity();
    for (int i = 0; i < n_local_cells_; ++i) {
        local_min_centroid_r = std::min(local_min_centroid_r,
                                        centroid_radius(i));
    }
    if (!std::isfinite(local_min_centroid_r)) local_min_centroid_r = 0.0;
    double global_min_centroid_r = local_min_centroid_r;
    MPI_Allreduce(&local_min_centroid_r, &global_min_centroid_r, 1,
                  MPI_DOUBLE, MPI_MIN, comm_);

    const double r_cav = std::max(0.0, cfg_.cavity_radius_m);
    // First pass: strict centroid <= cavity_radius_m.
    inner_cavity_cells_.clear();
    for (int i = 0; i < n_local_cells_; ++i) {
        if (centroid_radius(i) <= r_cav) inner_cavity_cells_.push_back(i);
    }
    long long local_n = static_cast<long long>(inner_cavity_cells_.size());
    long long global_n = local_n;
    MPI_Allreduce(&local_n, &global_n, 1, MPI_LONG_LONG, MPI_SUM, comm_);
    if (global_n == 0) {
        // Fallback: innermost cavity-wall layer.
        const double r_thr = global_min_centroid_r + 0.10 * std::max(1.0, r_cav);
        inner_cavity_cells_.clear();
        for (int i = 0; i < n_local_cells_; ++i) {
            if (centroid_radius(i) <= r_thr) inner_cavity_cells_.push_back(i);
        }
        local_n = static_cast<long long>(inner_cavity_cells_.size());
        MPI_Allreduce(&local_n, &global_n, 1, MPI_LONG_LONG, MPI_SUM, comm_);
    }

    if (global_n == 0) {
        int rank = 0;
        MPI_Comm_rank(comm_, &rank);
        if (rank == 0) {
            std::fprintf(stderr,
                "Source3DBallImpl: no inner-cavity cells found on the "
                "mesh (cavity_radius_m=%.3g); source forcing disabled.\n",
                cfg_.cavity_radius_m);
        }
        inner_cavity_cells_.clear();
        return;
    }

    // Inner-cavity volume bookkeeping for the per-cell energy weight.
    inner_cavity_volume_local_ = 0.0;
    double r_layer_edge_local = 0.0;
    for (int i : inner_cavity_cells_) {
        inner_cavity_volume_local_ += volume_[i];
        r_layer_edge_local = std::max(r_layer_edge_local, centroid_radius(i));
    }
    inner_cavity_volume_global_ = inner_cavity_volume_local_;
    MPI_Allreduce(&inner_cavity_volume_local_, &inner_cavity_volume_global_,
                  1, MPI_DOUBLE, MPI_SUM, comm_);
    if (inner_cavity_volume_global_ <= 0.0) {
        inner_cavity_volume_global_ = 1.0;  // defensive; never reached.
    }
    double r_layer_edge_global = r_layer_edge_local;
    MPI_Allreduce(&r_layer_edge_local, &r_layer_edge_global, 1,
                  MPI_DOUBLE, MPI_MAX, comm_);
    // Effective cavity radius for the (a/r)^3 elastostatic pressure
    // envelope: the larger of the configured cavity radius and the
    // outer edge of the tagged cavity-wall layer (so the envelope is
    // continuous across the cavity-wall / rock boundary). Falls back
    // to a small positive floor on degenerate meshes.
    inner_cavity_radius_a_ = std::max(std::max(0.0, cfg_.cavity_radius_m),
                                      r_layer_edge_global);
    if (inner_cavity_radius_a_ <= 0.0) inner_cavity_radius_a_ = 1.0e-3;

    // One-line per-rank diagnostic; a rank with zero inner-cavity cells
    // (which can happen at MPI > 1 when the mesh is partitioned so a
    // chunk holds none of the cavity wall) is a no-op for the deposit
    // step but still participates in the global energy bookkeeping.
    int rank = 0;
    int nranks = 1;
    MPI_Comm_rank(comm_, &rank);
    MPI_Comm_size(comm_, &nranks);
    if (nranks > 1 && inner_cavity_cells_.empty()) {
        std::fprintf(stderr,
            "Source3DBallImpl[rank %d]: no inner-cavity cells on this "
            "rank; energy deposit is a no-op here (global total still "
            "non-zero).\n", rank);
    }
    if (rank == 0) {
        std::fprintf(stderr,
            "Source3DBallImpl: source forcing active. yield=%.4g kt, "
            "E_total=%.4g J, STF duration=%.3g s, inner-cavity cells "
            "(global)=%lld, inner-cavity volume=%.4g m^3.\n",
            cfg_.source_yield_kt, total_source_energy_J_,
            cfg_.source_deposition_duration_s, global_n,
            inner_cavity_volume_global_);
    }

    // Build the SourceBallInnerCavity DMLabel on the DM (value 1 on the
    // tagged cells). Used by the HDF5 / ParaView source-ball output.
    DM dm = mesh_->getDM();
    if (dm) {
        DMLabel cav_label = nullptr;
        DMCreateLabel(dm, "SourceBallInnerCavity");
        DMGetLabel(dm, "SourceBallInnerCavity", &cav_label);
        if (cav_label) {
            PetscInt cStart = 0, cEnd = 0;
            DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd);
            for (int li : inner_cavity_cells_) {
                DMLabelSetValue(cav_label, cStart + li, 1);
            }
        }
    }

    // Tillotson host-rock EOS that maps the deposited specific internal
    // energy to a matter pressure (the cavity-wall pressure perturbation
    // that sets the amplitude of the elastostatic source field).
    source_tillotson_.setParameters(
        TillotsonParameterSets::byName(cfg_.medium_label));

    source_forcing_active_ = true;
}

double Source3DBallImpl::sourceTimeFunctionPowerJ(double t) const
{
    if (total_source_energy_J_ <= 0.0) return 0.0;
    return total_source_energy_J_
           * stfShapeRate(cfg_.source_time_function, t,
                          cfg_.source_deposition_duration_s);
}

void Source3DBallImpl::depositSourceEnergyStep(double dt)
{
    std::fill(source_power_density_.begin(), source_power_density_.end(), 0.0);
    if (!source_forcing_active_ || dt <= 0.0) return;

    // Trapezoidal estimate of the energy deposited globally over this
    // host step (the rate is smooth so the midpoint value is fine).
    const double t_mid = current_time_ - source_forcing_t0_ + 0.5 * dt;
    const double power_J_per_s = sourceTimeFunctionPowerJ(std::max(0.0, t_mid));
    const double dE_global = power_J_per_s * dt;
    if (dE_global <= 0.0) return;

    // Distribute by volume weight: dE_i = dE_global * V_i / V_total, so
    // the matter-energy increment per unit mass is dE_global / (rho_i *
    // V_total) (the cell volume cancels). Energy sums to dE_global.
    for (int i : inner_cavity_cells_) {
        const double rho_safe = (rho_[i] > 1.0e-9) ? rho_[i] : 1.0e-9;
        const double dE_i = dE_global * volume_[i] / inner_cavity_volume_global_;
        e_int_[i] += dE_i / (rho_safe * volume_[i]);
        source_power_density_[i] = power_J_per_s / inner_cavity_volume_global_;
    }
    deposited_energy_global_J_ += dE_global;
}

void Source3DBallImpl::updateCavityPressureField(std::vector<double>& d_dp_out)
{
    d_dp_out.assign(n_local_cells_, 0.0);
    if (n_local_cells_ == 0 || !source_forcing_active_) return;

    // Volume-averaged cavity-wall pressure perturbation (compression
    // positive) from the post-deposit specific internal energy via the
    // Tillotson host-rock EOS. Tillotson returns p > 0 in compression;
    // a negative value (rarefaction) would not represent a source, so
    // it is floored at zero.
    double num = 0.0;
    double den = 0.0;
    for (int i : inner_cavity_cells_) {
        const double p = std::max(0.0, source_tillotson_.pressure(rho_[i], e_int_[i]));
        num += p * volume_[i];
        den += volume_[i];
    }
    // MPI-reduce over the source-ball communicator (a no-op on the
    // PETSC_COMM_SELF the host runs the ball on; non-trivial only when
    // a caller explicitly distributes the ball, e.g. the standalone
    // tests at MPI > 1).
    double red_num = num;
    double red_den = den;
    MPI_Allreduce(&num, &red_num, 1, MPI_DOUBLE, MPI_SUM, comm_);
    MPI_Allreduce(&den, &red_den, 1, MPI_DOUBLE, MPI_SUM, comm_);
    const double dp_cavity = (red_den > 0.0) ? (red_num / red_den) : 0.0;
    dp_cavity_mean_ = dp_cavity;

    // Build the per-cell pressure-perturbation field: dp_cavity inside
    // the cavity-wall layer, dp_cavity * (a / r)^3 in the surrounding
    // rock (the Sharpe 1942 / Lame elastostatic pressurised-cavity
    // field; Mueller & Murphy 1971; Aki & Richards 2002 ch 4). A cell
    // with at least one vertex tagged inner-cavity carries the full
    // cavity pressure; the others follow the (a/r)^3 envelope clamped
    // to never exceed the cavity value.
    std::vector<unsigned char> is_cav(n_local_cells_, 0);
    for (int i : inner_cavity_cells_) {
        if (i >= 0 && i < n_local_cells_) is_cav[i] = 1;
    }
    const double a3 = inner_cavity_radius_a_ * inner_cavity_radius_a_
                      * inner_cavity_radius_a_;
    for (int i = 0; i < n_local_cells_; ++i) {
        const double prev = dp_cell_[i];
        double dp_target;
        if (is_cav[i]) {
            dp_target = dp_cavity;
        } else {
            const auto& c = centroid_[i];
            const double r = std::sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
            const double r3 = (r > 1.0e-12) ? (r * r * r) : 1.0e-36;
            const double env = std::min(1.0, a3 / r3);
            dp_target = dp_cavity * env;
        }
        dp_cell_[i] = dp_target;
        d_dp_out[i] = dp_target - prev;
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

    // Pass-14a: deposit the source-energy increment into the inner-
    // cavity cells, recompute the Tillotson cavity-wall pressure, and
    // update the per-cell elastostatic pressure field, before the
    // constitutive update. With source_forcing_enabled = false these
    // are no-ops, so the path reproduces pass-13c exactly.
    depositSourceEnergyStep(dt);
    std::vector<double> d_dp_field;
    updateCavityPressureField(d_dp_field);

    last_yielded_count_ = 0;
    for (int i = 0; i < n_local_cells_; ++i) {
        // Total strain increment: the host-supplied external strain rate
        // (pass-13b path; zero in the pass-14a delegated configuration)
        // plus the volumetric strain implied by the change in the
        // cavity-source pressure field, deps_vol = -d(dp) / K, split
        // isotropically. The radial return turns this into an isotropic
        // stress drop -d(dp) (plus a deviatoric drop wherever the cell
        // yields under the asymmetric overburden, which is what gives a
        // small CLVD content); the pass-13c surface integral over the
        // elastic radius then captures it as a real moment-rate tensor.
        std::array<double, 6> deps;
        for (int k = 0; k < 6; ++k) deps[k] = strain_rate_[i][k] * dt;
        if (source_forcing_active_ && !d_dp_field.empty()) {
            const double Kc = std::max(1.0, bulk_modulus_K_);
            const double deps_vol = -d_dp_field[i] / Kc;
            deps[0] += deps_vol / 3.0;
            deps[1] += deps_vol / 3.0;
            deps[2] += deps_vol / 3.0;
        }
        std::array<double, 6> sigma_new{};
        std::array<double, 6> deps_p{};
        auto res = druckerPrager3DRadialReturn(sigma_[i], deps,
                                                bulk_modulus_K_, shear_modulus_G_,
                                                dp_params_, sigma_new, deps_p);
        sigma_[i] = sigma_new;
        if (res.yielded) {
            ++last_yielded_count_;
            eps_p_eq_[i] += res.delta_eps_p_eq;
            const double Y_proj = dp_params_.alpha_dp
                * (-(sigma_new[0] + sigma_new[1] + sigma_new[2]) / 3.0)
                + dp_params_.k_dp;
            const double dW = Y_proj * res.delta_eps_p_eq;
            const double rho_safe = (rho_[i] > 1.0e-9) ? rho_[i] : 1.0e-9;
            e_int_[i] += dW / rho_safe;
        }
    }

    // Pass-13c radiation cadence (unchanged).
    rad_substep_counter_ += 1;
    const int rcad = std::max(1, cfg_.radiation_substep_cadence);
    if (rad_solver_ && rad_solver_->isInitialized()
        && (rad_substep_counter_ % rcad) == 0) {
        auto rres = rad_solver_->step(rcad * dt, rho_, T_m_, E_r_, e_int_);
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

// The pass-13c surface integral M_surf = integral_S (sigma(t) -
// sigma(0)) n x dA equals the volume integral of the stress drop over
// the source region (divergence theorem on an equilibrium stress
// field). The seismic moment tensor convention used by the far-field
// FEM injector and the 1D radial path is M_ij = -integral_V (stress
// drop)_ij dV: an explosion (a compressive stress increase, negative
// in FSRM's compression-negative convention) has a positive isotropic
// moment. So when source forcing drives the stress drop the reported
// moment is the negation of the raw surface integral; without source
// forcing the raw value is preserved (so the pass-13c synthetic-stress
// unit gates, which call recomputeMomentTensorFromSurface directly, are
// unchanged).
namespace {
inline std::array<double, 6> negTensor(const std::array<double, 6>& M)
{
    return {-M[0], -M[1], -M[2], -M[3], -M[4], -M[5]};
}
}  // namespace

void Source3DBallImpl::getMomentTensor(std::array<double, 6>& M) const
{
    M = source_forcing_active_ ? negTensor(M_global_) : M_global_;
}

void Source3DBallImpl::getMomentRateTensor(std::array<double, 6>& Mdot) const
{
    Mdot = source_forcing_active_ ? negTensor(Mdot_global_) : Mdot_global_;
}

void Source3DBallImpl::getState(Source3DBallState& state) const
{
    state = state_snapshot_;
    if (source_forcing_active_) {
        state.moment_tensor_N_m = negTensor(state.moment_tensor_N_m);
        state.moment_rate_tensor_N_m_per_s =
            negTensor(state.moment_rate_tensor_N_m_per_s);
    }
}

const char* Source3DBallImpl::name() const
{
    return kPhysicsName;
}

void Source3DBallImpl::getCellStates(std::vector<Source3DCellState>& out) const
{
    out.clear();
    out.reserve(n_local_cells_);
    // Build an O(1) inner-cavity membership lookup.
    std::vector<unsigned char> is_cav(n_local_cells_, 0);
    for (int i : inner_cavity_cells_) {
        if (i >= 0 && i < n_local_cells_) is_cav[i] = 1;
    }
    for (int i = 0; i < n_local_cells_; ++i) {
        Source3DCellState s;
        s.centroid = centroid_[i];
        s.sigma = sigma_[i];
        s.eps_p_eq = eps_p_eq_[i];
        s.rho = rho_[i];
        s.e_int = e_int_[i];
        s.T_m = T_m_[i];
        s.E_r = E_r_[i];
        s.source_forcing_power =
            (i < static_cast<int>(source_power_density_.size()))
                ? source_power_density_[i] : 0.0;
        s.inner_cavity_marker = is_cav[i] ? 1.0 : 0.0;
        s.dp_pressure_field = (i < static_cast<int>(dp_cell_.size()))
                              ? dp_cell_[i] : 0.0;
        out.push_back(s);
    }
}

bool Source3DBallImpl::getMeshGeometry(
    std::vector<std::array<float, 3>>& verts,
    std::vector<std::array<int, 4>>& tets) const
{
    verts.clear();
    tets.clear();
    if (cfg_.mesh_path.empty())
        return false;

    // Read .node file  (format: N_verts dim 0 1 / idx x y z attr)
    {
        std::string node_path = cfg_.mesh_path + ".node";
        std::ifstream fin(node_path);
        if (!fin.is_open())
            return false;
        int n_verts = 0, dim = 0, n_attr = 0, n_bmark = 0;
        fin >> n_verts >> dim >> n_attr >> n_bmark;
        if (n_verts <= 0 || dim != 3)
            return false;
        verts.reserve(static_cast<size_t>(n_verts));
        for (int i = 0; i < n_verts; ++i)
        {
            int idx;
            float x, y, z;
            fin >> idx >> x >> y >> z;
            // skip n_attr attributes and boundary marker
            for (int a = 0; a < n_attr + n_bmark; ++a)
            {
                float dummy;
                fin >> dummy;
            }
            verts.push_back({x, y, z});
        }
    }

    // Read .ele file  (format: N_tets 4 n_attr / idx n1 n2 n3 n4 [attr])
    {
        std::string ele_path = cfg_.mesh_path + ".ele";
        std::ifstream fin(ele_path);
        if (!fin.is_open())
        {
            verts.clear();
            return false;
        }
        int n_tets = 0, npc = 0, n_attr = 0;
        fin >> n_tets >> npc >> n_attr;
        if (n_tets <= 0 || npc != 4)
        {
            verts.clear();
            return false;
        }
        tets.reserve(static_cast<size_t>(n_tets));
        for (int i = 0; i < n_tets; ++i)
        {
            int idx, n0, n1, n2, n3;
            fin >> idx >> n0 >> n1 >> n2 >> n3;
            for (int a = 0; a < n_attr; ++a)
            {
                int dummy;
                fin >> dummy;
            }
            // TetGen is 1-indexed; convert to 0-indexed
            tets.push_back({n0 - 1, n1 - 1, n2 - 1, n3 - 1});
        }
    }
    return true;
}

}  // namespace FSRM

