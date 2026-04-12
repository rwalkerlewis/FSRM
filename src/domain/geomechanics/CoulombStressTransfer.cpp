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
#include <petscfe.h>
#include <petscdt.h>
#include <cmath>
#include <algorithm>
#include <limits>
#include <vector>

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

    // Get cell range for height-0 (volume cells)
    PetscInt cStart, cEnd;
    ierr = DMPlexGetHeightStratum(dm_, 0, &cStart, &cEnd); CHKERRQ(ierr);

    // Get the displacement FE from the DM.
    // Field ordering: displacement is field 0 for elastostatics,
    // field 1 for poroelastic (field 0 = pressure).
    // Determine displacement field index by checking field component count.
    PetscInt nFields;
    ierr = DMGetNumFields(dm_, &nFields); CHKERRQ(ierr);

    PetscInt disp_field = 0;
    PetscFE fe_disp = nullptr;
    for (PetscInt f = 0; f < nFields; ++f) {
      PetscObject obj;
      ierr = DMGetField(dm_, f, nullptr, &obj); CHKERRQ(ierr);
      PetscFE fe = (PetscFE)obj;
      PetscInt nc;
      ierr = PetscFEGetNumComponents(fe, &nc); CHKERRQ(ierr);
      if (nc == dim) {
        disp_field = f;
        fe_disp = fe;
        break;
      }
    }

    // Tabulate basis function derivatives at the cell centroid (reference origin).
    // For simplex affine elements, the gradients are constant over the cell,
    // so the reference point does not matter; we use the origin {0,0,0}.
    PetscTabulation tab = nullptr;
    if (fe_disp) {
      std::vector<PetscReal> refPoint(dim, 0.0);
      // K=1: tabulate values (T[0]) and first derivatives (T[1])
      ierr = PetscFECreateTabulation(fe_disp, 1, 1, refPoint.data(), 1, &tab); CHKERRQ(ierr);
    }

    PetscInt Nb = 0;  // number of basis functions
    if (tab) {
      Nb = tab->Nb;
    }

    // For each fault vertex, find an adjacent volume cell and compute
    // the displacement gradient using FEM basis function derivatives
    size_t nv = fault_->numVertices();
    for (size_t vi = 0; vi < nv; ++vi) {
      const FaultVertex& fv = fault_->getVertex(vi);

      // Use the negative-side vertex to find the adjacent cell
      PetscInt vert = fv.vertex_negative;
      if (vert < 0) continue;

      // Walk the DAG upward from the vertex to find a volume cell.
      // DMPlexGetTransitiveClosure with useCone=PETSC_FALSE gives the
      // star (support closure) of the point: all mesh entities that
      // contain this vertex, including edges, faces, and cells.
      PetscInt *closure = nullptr;
      PetscInt closureSize = 0;
      ierr = DMPlexGetTransitiveClosure(dm_, vert, PETSC_FALSE, &closureSize, &closure);
      CHKERRQ(ierr);

      PetscInt cell = -1;
      for (PetscInt ci = 0; ci < closureSize; ++ci) {
        PetscInt p = closure[2 * ci];
        if (p >= cStart && p < cEnd) {
          cell = p;
          break;
        }
      }
      ierr = DMPlexRestoreTransitiveClosure(dm_, vert, PETSC_FALSE, &closureSize, &closure);
      CHKERRQ(ierr);

      if (cell < 0) continue;

      // Get the inverse Jacobian for this cell (affine geometry)
      // v0[dim]: cell vertex, J[dim*dim]: Jacobian, invJ[dim*dim]: inverse, detJ
      std::vector<PetscReal> v0(dim), J(dim * dim), invJ(dim * dim);
      PetscReal detJ;
      ierr = DMPlexComputeCellGeometryAffineFEM(dm_, cell, v0.data(), J.data(),
                                                 invJ.data(), &detJ);
      CHKERRQ(ierr);

      // Get the displacement DOFs for this cell via closure
      PetscScalar *cellDofs = nullptr;
      PetscInt cellDofSize;
      ierr = DMPlexVecGetClosure(dm_, section, local, cell, &cellDofSize, &cellDofs);
      CHKERRQ(ierr);

      // Extract displacement DOFs.
      // The closure contains DOFs for ALL fields, ordered as:
      //   [field_0: Nb0 values] [field_1: Nb1 values] ...
      // where Nb includes all components (e.g., Nb=12 for P1 vector in 3D).
      PetscInt disp_offset = 0;
      for (PetscInt f = 0; f < disp_field; ++f) {
        PetscObject obj;
        ierr = DMGetField(dm_, f, nullptr, &obj); CHKERRQ(ierr);
        PetscFE fe_f = (PetscFE)obj;
        PetscInt nb_f;
        ierr = PetscFEGetDimension(fe_f, &nb_f); CHKERRQ(ierr);
        disp_offset += nb_f;
      }

      // Sample pressure from the solution at this vertex (for Biot coupling)
      double P = 0.0;
      {
        const PetscScalar* sol_array;
        ierr = VecGetArrayRead(local, &sol_array); CHKERRQ(ierr);
        PetscInt vdof, voff;
        ierr = PetscSectionGetDof(section, vert, &vdof); CHKERRQ(ierr);
        ierr = PetscSectionGetOffset(section, vert, &voff); CHKERRQ(ierr);
        // Pressure is field 0 if there are more fields than just displacement.
        // Check if the displacement field is > 0 (meaning pressure is field 0).
        if (disp_field > 0 && vdof > 0) {
          P = PetscRealPart(sol_array[voff]);
        }
        ierr = VecRestoreArrayRead(local, &sol_array); CHKERRQ(ierr);
      }

      // Compute displacement gradient: du_c/dx_j
      // For PETSc vector Lagrange FE, each basis function b controls a single
      // component c = b % Nc. The tabulation D[b][c][k] = dN_b^c / dX_k is
      // nonzero only when c == b % Nc.
      //
      // Physical gradient:
      //   du_c/dx_j = sum_b closureDofs[disp_offset + b] *
      //               sum_k D[b*Nc*cdim + c*cdim + k] * invJ[k*dim + j]
      // where c = b % dim.
      double grad_u[3][3] = {{0}};

      if (tab && Nb > 0) {
        const PetscReal *D = tab->T[1];

        for (PetscInt b = 0; b < Nb; ++b) {
          PetscInt c = b % dim;  // component this basis function controls
          double u_val = PetscRealPart(cellDofs[disp_offset + b]);

          for (PetscInt j = 0; j < dim; ++j) {
            double phys_grad = 0.0;
            for (PetscInt k = 0; k < dim; ++k) {
              phys_grad += D[b * dim * dim + c * dim + k] * invJ[k * dim + j];
            }
            grad_u[c][j] += u_val * phys_grad;
          }
        }
      }

      ierr = DMPlexVecRestoreClosure(dm_, section, local, cell, &cellDofSize, &cellDofs);
      CHKERRQ(ierr);

      // Symmetrize to get strain: eps_ij = 0.5 * (du_i/dx_j + du_j/dx_i)
      double eps[6];
      eps[0] = grad_u[0][0];                                // eps_xx
      eps[1] = grad_u[1][1];                                // eps_yy
      eps[2] = (dim > 2) ? grad_u[2][2] : 0.0;             // eps_zz
      eps[3] = 0.5 * (grad_u[0][1] + grad_u[1][0]);        // eps_xy
      eps[4] = (dim > 2) ? 0.5 * (grad_u[0][2] + grad_u[2][0]) : 0.0; // eps_xz
      eps[5] = (dim > 2) ? 0.5 * (grad_u[1][2] + grad_u[2][1]) : 0.0; // eps_yz

      // Compute stress from strain (including Biot pressure coupling)
      VertexStress& vs = current_[vi];
      vs.pressure = P;
      computeStressFromStrain(eps, P, vs.sigma.data());
    }

    if (tab) {
      ierr = PetscTabulationDestroy(&tab); CHKERRQ(ierr);
    }

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
