/**
 * @file Source3DBallImpl.hpp
 * @brief Pass-13a (axis-1b foundation) concrete implementation of the
 *        Source3DBall interface. Replaces the pass-11 throw-on-construct
 *        scaffold with a real instantiation that loads a TetGen-generated
 *        mesh into a distributed DMPlex via Source3DBallMesh.
 *
 * What pass-13a delivers (the foundation slice):
 *  - initialize(cfg): loads the mesh referenced by cfg.mesh_path via
 *    Source3DBallMesh (DMPlex create + DMPlexDistribute), validates
 *    the geometry against cfg.cavity_radius_m / cfg.outer_radius_m,
 *    captures vertex-marker counts.
 *  - getMomentTensor / getMomentRateTensor: return the zero tensor.
 *    The 3D constitutive update and surface-integral extraction are
 *    pass-13b/c work; until then there is no moment tensor to report.
 *  - getState: returns a snapshot with mesh stats but zero kinematics.
 *  - step(dt): throws std::runtime_error referencing pass-13b. The 3D
 *    hydro substep is the pass-13b deliverable.
 *  - name(): "Source3DBallImpl_v1_pass13a_skeleton". Diagnostic only.
 *
 * What pass-13a explicitly does not implement (named for follow-on
 * passes; throws or no-ops with clear messages):
 *  - 3D Drucker-Prager radial return (pass-13b)
 *  - Asymmetric overburden initial stress (pass-13b)
 *  - Cell-centred FV grey radiation diffusion (pass-13b)
 *  - Surface-integral moment-tensor extraction (pass-13c)
 *  - End-to-end Salmon at MPI=4 (pass-13c gate)
 *  - HDF5 / XDMF output (pass-13c)
 *
 * Backward compatibility: cavity_geometry = SPHERICAL remains the
 * default and entirely bypasses this code path. The 32 historic-event
 * integration tests do not exercise Source3DBallImpl. Pass-13a
 * intentionally keeps the RadialLagrangianSolver throw on
 * cavity_geometry = THREE_DIMENSIONAL because the host-side
 * delegation from RadialLagrangianSolver into Source3DBallImpl is
 * pass-13b/c work; the foundation slice provides only the leaf
 * Source3DBallImpl that the host will call into later.
 */

#ifndef NEAR_FIELD_SOURCE_3D_BALL_IMPL_HPP
#define NEAR_FIELD_SOURCE_3D_BALL_IMPL_HPP

#include "domain/explosion/Source3DBall.hpp"

#include <memory>
#include <mpi.h>

namespace FSRM {

class Source3DBallMesh;

class Source3DBallImpl final : public Source3DBall
{
public:
    Source3DBallImpl();
    ~Source3DBallImpl() override;

    Source3DBallImpl(const Source3DBallImpl&) = delete;
    Source3DBallImpl& operator=(const Source3DBallImpl&) = delete;

    /// Override the MPI communicator the underlying DMPlex lives on.
    /// Defaults to PETSC_COMM_WORLD when initialize() is called
    /// without first calling this. Pass-13b will wire this from the
    /// host's existing communicator selection logic.
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

private:
    MPI_Comm comm_ = MPI_COMM_NULL;
    bool comm_set_ = false;
    bool mesh_loaded_ = false;
    Source3DBallConfig cfg_{};
    std::unique_ptr<Source3DBallMesh> mesh_;
    Source3DBallState state_snapshot_{};
};

}  // namespace FSRM

#endif  // NEAR_FIELD_SOURCE_3D_BALL_IMPL_HPP
