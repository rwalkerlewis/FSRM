/**
 * @file Source3DBallImpl.cpp
 * @brief Pass-13a (axis-1b foundation) Source3DBallImpl implementation.
 *        See Source3DBallImpl.hpp for the contract and the explicit
 *        pass-13a / pass-13b / pass-13c slice boundaries.
 */

#include "domain/explosion/Source3DBallImpl.hpp"

#include "domain/explosion/Source3DBallMesh.hpp"

#include <petscsys.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace FSRM {

namespace
{

constexpr const char* kFoundationName = "Source3DBallImpl_pass-13a_foundation_skeleton";

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

    // Pass-13a: NODAL_FEM is named-only on the radiation ladder; the
    // cell-centred FV path is the only viable choice for pass-13b
    // grey diffusion. Reject NODAL_FEM with a clear pass-15+ message
    // so callers do not silently opt into a non-existent code path.
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

    if (cfg.mesh_path.empty()) {
        // Foundation no-op mode: caller wants to verify the factory
        // returns a valid object without touching the disk. Pass-13b
        // will require a mesh path; pass-13a permits the empty-path
        // construction so the BackwardCompat factory-instantiation
        // regression test does not require a fixture mesh.
        mesh_loaded_ = false;
        state_snapshot_ = Source3DBallState{};
        return;
    }

    mesh_ = std::make_unique<Source3DBallMesh>();
    mesh_->loadFromTetGen(comm_, cfg.mesh_path);
    mesh_loaded_ = true;

    // Pass-13a sanity-check the loaded mesh against the configured
    // cavity / elastic radii. The mesh's local extrema are gathered
    // across ranks for a global comparison.
    double local_min = mesh_->localMinVertexRadius();
    double local_max = mesh_->localMaxVertexRadius();
    double global_min = local_min;
    double global_max = local_max;
    if (mesh_->numLocalVertices() == 0) {
        // Rank with no local vertices contributes a neutral value to
        // the reduction. Use +inf for min and 0 for max.
        local_min = std::numeric_limits<double>::infinity();
        local_max = 0.0;
    }
    MPI_Allreduce(&local_min, &global_min, 1, MPI_DOUBLE, MPI_MIN, comm_);
    MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, comm_);

    if (cfg.cavity_radius_m > 0.0) {
        const double rel = std::abs(global_min - cfg.cavity_radius_m)
                           / cfg.cavity_radius_m;
        if (rel > 0.20) {
            // Pass-13a tolerates a 20 percent envelope between the
            // configured cavity radius and the loaded mesh's minimum
            // vertex radius because the icosphere tessellation is not
            // a perfect sphere. A larger drift signals a config /
            // mesh-file mismatch and the foundation slice should
            // refuse to proceed silently.
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

    state_snapshot_.time = 0.0;
    state_snapshot_.cavity_radius_max_m = global_min;
    state_snapshot_.cavity_radius_min_m = global_min;
    state_snapshot_.moment_tensor_N_m = {0, 0, 0, 0, 0, 0};
    state_snapshot_.moment_rate_tensor_N_m_per_s = {0, 0, 0, 0, 0, 0};
}

void Source3DBallImpl::step(double dt)
{
    (void)dt;
    throw std::runtime_error(
        "Source3DBallImpl::step is pass-13b work, not yet implemented. "
        "Pass-13a foundation only loads the DMPlex mesh and validates "
        "geometry; the 3D Drucker-Prager radial return, asymmetric "
        "overburden initial state, and cell-centred FV grey radiation "
        "diffusion ship in pass-13b. See docs/AXIS_1B_DESIGN.md.");
}

void Source3DBallImpl::getMomentTensor(std::array<double, 6>& M) const
{
    // Pass-13a: no constitutive update has run, so the moment tensor
    // is identically zero. Pass-13c surface-integral extraction
    // populates this from the boundary stress at the elastic radius.
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
    return kFoundationName;
}

}  // namespace FSRM
