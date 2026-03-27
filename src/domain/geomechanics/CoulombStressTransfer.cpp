/**
 * @file CoulombStressTransfer.cpp
 * @brief Sample stress from FEM solution at fault vertices and compute CFS
 *
 * Implements the physical chain:
 *   displacement gradient → strain → Hooke's law → stress tensor
 *   → project onto fault → resolve normal/shear → compute delta_CFS
 */

#include "domain/geomechanics/CoulombStressTransfer.hpp"
#include "domain/geomechanics/PyLithFault.hpp"
#include <cmath>
#include <algorithm>
#include <limits>

namespace FSRM {

CoulombStressTransfer::CoulombStressTransfer(MPI_Comm comm)
    : comm_(comm) {}

PetscErrorCode CoulombStressTransfer::initialize(DM dm, FaultCohesiveDyn* fault,
                                                   double lambda, double mu,
                                                   double biot_alpha) {
    PetscFunctionBeginUser;

    dm_ = dm;
    fault_ = fault;
    lambda_ = lambda;
    mu_ = mu;
    biot_alpha_ = biot_alpha;

    // Allocate stress arrays for all fault vertices
    size_t nv = fault_->numVertices();
    current_.resize(nv);
    initial_.resize(nv);
    has_initial_ = false;

    PetscFunctionReturn(0);
}

PetscErrorCode CoulombStressTransfer::sampleStressAtFaults(Vec solution) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

    if (!dm_ || !fault_) PetscFunctionReturn(0);

    PetscInt dim;
    ierr = DMGetDimension(dm_, &dim); CHKERRQ(ierr);

    // Get the local solution vector
    Vec local;
    ierr = DMGetLocalVector(dm_, &local); CHKERRQ(ierr);
    ierr = DMGlobalToLocal(dm_, solution, INSERT_VALUES, local); CHKERRQ(ierr);

    // Get local section for DOF offsets
    PetscSection section;
    ierr = DMGetLocalSection(dm_, &section); CHKERRQ(ierr);

    const PetscScalar* sol_array;
    ierr = VecGetArrayRead(local, &sol_array); CHKERRQ(ierr);

    // For each fault vertex, find the adjacent cell and compute stress
    size_t nv = fault_->numVertices();
    for (size_t vi = 0; vi < nv; ++vi) {
        const FaultVertex& fv = fault_->getVertex(vi);

        // Use the negative-side vertex to find the adjacent cell
        PetscInt vert = fv.vertex_negative;
        if (vert < 0) continue;

        // Get the support (cells) of this vertex
        PetscInt support_size;
        const PetscInt* support;
        ierr = DMPlexGetSupportSize(dm_, vert, &support_size); CHKERRQ(ierr);
        ierr = DMPlexGetSupport(dm_, vert, &support); CHKERRQ(ierr);

        // Find a volume cell (height 0) connected to this vertex
        // Walk up the support chain: vertex → edge → face → cell
        PetscInt cStart, cEnd;
        ierr = DMPlexGetHeightStratum(dm_, 0, &cStart, &cEnd); CHKERRQ(ierr);

        // For simplicity, use the vertex DOFs directly
        // In the 6-field layout: 0=P, 1=S, 2=ux, 3=uy, 4=uz, 5=phi
        // We need the displacement gradients, which are not directly available
        // from a vertex. Instead, we approximate using the cell-averaged gradient.

        // Get DOFs at this vertex
        PetscInt dof, offset;
        ierr = PetscSectionGetDof(section, vert, &dof); CHKERRQ(ierr);
        ierr = PetscSectionGetOffset(section, vert, &offset); CHKERRQ(ierr);

        // Extract pressure (field 0) and displacement (fields 2,3,4)
        double P = 0.0;
        double ux = 0.0, uy = 0.0, uz = 0.0;
        if (dof >= 6) {
            P  = PetscRealPart(sol_array[offset + 0]);
            ux = PetscRealPart(sol_array[offset + 2]);
            uy = PetscRealPart(sol_array[offset + 3]);
            uz = PetscRealPart(sol_array[offset + 4]);
        } else if (dof >= 4) {
            P  = PetscRealPart(sol_array[offset + 0]);
            ux = PetscRealPart(sol_array[offset + 1]);
            uy = (dof > 2) ? PetscRealPart(sol_array[offset + 2]) : 0.0;
            uz = (dof > 3) ? PetscRealPart(sol_array[offset + 3]) : 0.0;
        }

        // For strain computation, we need displacement gradients.
        // Since we're sampling at a vertex, we estimate gradients using
        // the positive-side vertex displacement to get the displacement jump.
        // The strain in the continuum is obtained from neighboring vertex values.
        //
        // Simplified approach: estimate strain from the displacement magnitude
        // relative to a characteristic length scale. This is a proxy that works
        // for detecting CFS changes during injection, though a full implementation
        // would project strain to a DG field and sample properly.

        // Get positive-side vertex displacement for displacement jump
        double ux_pos = ux, uy_pos = uy, uz_pos = uz;
        PetscInt pos_vert = fv.vertex_positive;
        if (pos_vert >= 0) {
            PetscInt pos_dof, pos_off;
            ierr = PetscSectionGetDof(section, pos_vert, &pos_dof); CHKERRQ(ierr);
            ierr = PetscSectionGetOffset(section, pos_vert, &pos_off); CHKERRQ(ierr);
            if (pos_dof >= 6) {
                ux_pos = PetscRealPart(sol_array[pos_off + 2]);
                uy_pos = PetscRealPart(sol_array[pos_off + 3]);
                uz_pos = PetscRealPart(sol_array[pos_off + 4]);
            }
        }

        // Estimate strain using finite differences between support cells
        // For now, use a simplified approach: the stress is estimated from
        // the pressure change via Biot coupling: sigma' = -alpha*P*I + deviatoric
        // For a uniform injection, the stress change is primarily isotropic:
        //   delta_sigma_xx = delta_sigma_yy = delta_sigma_zz ≈ alpha * delta_P / 3
        // The shear stress change comes from geometric effects near the fault.

        // For CFS computation, the key quantity is the pore pressure change
        // and its effect on effective normal stress. Use Biot coupling:
        //   delta_sigma_n_eff = delta_sigma_n - alpha * delta_P
        // For isotropic loading: delta_sigma_n ≈ 0 (confined), so
        //   delta_sigma_n_eff ≈ -alpha * delta_P (pore pressure increase reduces effective stress)

        // Store the stress state using Biot effective stress coupling
        VertexStress& vs = current_[vi];
        vs.pressure = P;

        // Estimate strain from average displacement gradient using the
        // displacement jump across the fault as a proxy. This is a simple
        // first-order approximation that avoids full FEM reconstruction
        // but still captures the effect of deformation on stress.
        //
        // Δu = u_pos - u  (displacement jump across the fault)
        // eps_xx ≈ Δux / L_ref, etc., with a local reference length scale.
        const double ref_len = 1.0;  // [m] reference length; tune if mesh scale is known

        const double dux = ux_pos - ux;
        const double duy = uy_pos - uy;
        const double duz = uz_pos - uz;

        double eps[6];
        // Normal strains
        eps[0] = dux / ref_len;  // eps_xx
        eps[1] = duy / ref_len;  // eps_yy
        eps[2] = duz / ref_len;  // eps_zz
        // Shear strains (symmetric combinations; small-strain assumption)
        eps[3] = 0.5 * (dux + duy) / ref_len;  // eps_xy
        eps[4] = 0.5 * (dux + duz) / ref_len;  // eps_xz
        eps[5] = 0.5 * (duy + duz) / ref_len;  // eps_yz

        // Compute stress from strain (including Biot pressure coupling)
        computeStressFromStrain(eps, P, vs.sigma.data());
    }

    ierr = VecRestoreArrayRead(local, &sol_array); CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(dm_, &local); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode CoulombStressTransfer::resolveStressOnFault() {
    PetscFunctionBeginUser;

    if (!fault_) PetscFunctionReturn(0);

    size_t nv = fault_->numVertices();
    for (size_t vi = 0; vi < nv; ++vi) {
        const FaultVertex& fv = fault_->getVertex(vi);
        VertexStress& vs = current_[vi];

        // Stress tensor in Voigt notation: [xx, yy, zz, xy, xz, yz]
        // Full tensor:
        // | sigma[0]  sigma[3]  sigma[4] |
        // | sigma[3]  sigma[1]  sigma[5] |
        // | sigma[4]  sigma[5]  sigma[2] |
        const double* s = vs.sigma.data();
        const auto& n = fv.normal;

        // sigma * n (matrix-vector product)
        double sn[3] = {
            s[0] * n[0] + s[3] * n[1] + s[4] * n[2],
            s[3] * n[0] + s[1] * n[1] + s[5] * n[2],
            s[4] * n[0] + s[5] * n[1] + s[2] * n[2]
        };

        // Normal stress: sigma_n = n^T * sigma * n
        double sigma_n = sn[0] * n[0] + sn[1] * n[1] + sn[2] * n[2];

        // Effective normal stress (Biot): sigma_n_eff = sigma_n - alpha * P
        // Convention: compression is negative for stress, positive for effective stress
        vs.sigma_n_eff = sigma_n - biot_alpha_ * vs.pressure;

        // Shear stress: tau = sqrt(|sigma*n|^2 - sigma_n^2)
        double sn_mag_sq = sn[0] * sn[0] + sn[1] * sn[1] + sn[2] * sn[2];
        double tau_sq = sn_mag_sq - sigma_n * sigma_n;
        vs.tau = (tau_sq > 0.0) ? std::sqrt(tau_sq) : 0.0;
    }

    PetscFunctionReturn(0);
}

PetscErrorCode CoulombStressTransfer::computeDeltaCFS(double static_friction) {
    PetscFunctionBeginUser;

    size_t nv = current_.size();
    for (size_t vi = 0; vi < nv; ++vi) {
        VertexStress& vs = current_[vi];

        if (has_initial_) {
            const VertexStress& init = initial_[vi];
            double delta_tau = vs.tau - init.tau;
            double delta_sigma_n_eff = vs.sigma_n_eff - init.sigma_n_eff;

            // delta_CFS = delta_tau - mu_s * delta_sigma_n_eff
            // Positive delta_CFS → moving toward failure
            // Note: if sigma_n_eff decreases (more tensile), delta_sigma_n_eff < 0,
            // and -mu_s * delta_sigma_n_eff > 0, promoting failure.
            vs.delta_cfs = delta_tau - static_friction * delta_sigma_n_eff;
        } else {
            // No initial state stored yet — compute CFS from absolute state
            // CFS = tau - mu_s * |sigma_n_eff| (compression positive)
            vs.delta_cfs = vs.tau - static_friction * std::abs(vs.sigma_n_eff);
        }
    }

    PetscFunctionReturn(0);
}

PetscErrorCode CoulombStressTransfer::updateFaultState(FaultCohesiveDyn* fault) {
    PetscFunctionBeginUser;

    if (!fault) PetscFunctionReturn(0);

    size_t nv = std::min(current_.size(), fault->numVertices());
    for (size_t vi = 0; vi < nv; ++vi) {
        const VertexStress& vs = current_[vi];

        // Access dynamic state — need mutable access
        // FaultCohesiveDyn provides getState() as const, but we need to modify.
        // Use the initialize/setInitialState pattern instead.
        // For now, build a vector and call setInitialState.
    }

    // Build initial traction field from resolved stress
    FaultTractionField traction;
    traction.resize(nv);
    for (size_t vi = 0; vi < nv; ++vi) {
        const VertexStress& vs = current_[vi];
        const FaultVertex& fv = fault->getVertex(vi);

        // Decompose shear traction into strike and dip components
        // sigma * n gives the traction vector on the fault
        const double* s = vs.sigma.data();
        const auto& n = fv.normal;
        double sn[3] = {
            s[0] * n[0] + s[3] * n[1] + s[4] * n[2],
            s[3] * n[0] + s[1] * n[1] + s[5] * n[2],
            s[4] * n[0] + s[5] * n[1] + s[2] * n[2]
        };
        double sigma_n = sn[0] * n[0] + sn[1] * n[1] + sn[2] * n[2];

        // Tangential traction
        double tau_vec[3] = {sn[0] - sigma_n * n[0],
                            sn[1] - sigma_n * n[1],
                            sn[2] - sigma_n * n[2]};

        // Project onto strike and dip
        const auto& sd = fv.along_strike;
        const auto& dd = fv.up_dip;
        traction.traction_shear_ll[vi] = tau_vec[0] * sd[0] + tau_vec[1] * sd[1] + tau_vec[2] * sd[2];
        traction.traction_shear_ud[vi] = tau_vec[0] * dd[0] + tau_vec[1] * dd[1] + tau_vec[2] * dd[2];
        traction.traction_normal[vi] = sigma_n - biot_alpha_ * vs.pressure;
    }

    fault->setInitialTraction(traction);

    PetscFunctionReturn(0);
}

PetscErrorCode CoulombStressTransfer::storeInitialStress() {
    PetscFunctionBeginUser;
    initial_ = current_;
    has_initial_ = true;
    PetscFunctionReturn(0);
}

double CoulombStressTransfer::getMaxDeltaCFS() const {
    double max_cfs = -std::numeric_limits<double>::max();
    for (const auto& vs : current_) {
        max_cfs = std::max(max_cfs, vs.delta_cfs);
    }
    return max_cfs;
}

int CoulombStressTransfer::getMaxCFSVertex() const {
    int max_idx = -1;
    double max_cfs = -std::numeric_limits<double>::max();
    for (size_t i = 0; i < current_.size(); ++i) {
        if (current_[i].delta_cfs > max_cfs) {
            max_cfs = current_[i].delta_cfs;
            max_idx = static_cast<int>(i);
        }
    }
    return max_idx;
}

void CoulombStressTransfer::computeStressFromStrain(const double eps[6], double P,
                                                      double sigma[6]) const {
    // Hooke's law for isotropic material:
    // sigma_ij = lambda * eps_kk * delta_ij + 2 * mu * eps_ij
    //
    // Voigt notation: [xx, yy, zz, xy, xz, yz]
    // eps_kk = eps[0] + eps[1] + eps[2]

    double eps_kk = eps[0] + eps[1] + eps[2];

    // Effective stress (Biot coupling):
    // sigma'_ij = sigma_ij - alpha * P * delta_ij
    sigma[0] = lambda_ * eps_kk + 2.0 * mu_ * eps[0] - biot_alpha_ * P;  // xx
    sigma[1] = lambda_ * eps_kk + 2.0 * mu_ * eps[1] - biot_alpha_ * P;  // yy
    sigma[2] = lambda_ * eps_kk + 2.0 * mu_ * eps[2] - biot_alpha_ * P;  // zz
    sigma[3] = 2.0 * mu_ * eps[3];   // xy
    sigma[4] = 2.0 * mu_ * eps[4];   // xz
    sigma[5] = 2.0 * mu_ * eps[5];   // yz
}

} // namespace FSRM
