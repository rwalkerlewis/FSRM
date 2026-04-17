/**
 * @file DerivedFieldComputer.cpp
 * @brief Compute cell-centered stress, strain, and CFS from the FEM solution
 *
 * Uses PetscFECreateTabulation to compute basis function derivatives at cell
 * centroids, then maps displacement gradients to physical space via the
 * inverse Jacobian. This is the same FEM strain recovery approach used by
 * CoulombStressTransfer::sampleStressAtFaults().
 */

#include "core/DerivedFieldComputer.hpp"
#include <petscviewerhdf5.h>
#include <cmath>
#include <algorithm>

namespace FSRM
{

void DerivedFieldComputer::setMaterialParameters(double lambda, double mu, double biot_alpha)
{
  lambda_ = lambda;
  mu_ = mu;
  biot_alpha_ = biot_alpha;
}

void DerivedFieldComputer::setReceiverOrientation(double strike_deg, double dip_deg,
                                                   double friction)
{
  strike_deg_ = strike_deg;
  dip_deg_ = dip_deg;
  friction_ = friction;
  computeReceiverNormal();
}

void DerivedFieldComputer::computeReceiverNormal()
{
  // Convert strike/dip to a unit normal vector.
  // Strike is clockwise from north (y-axis), dip is from horizontal.
  // For a vertical fault (dip=90) striking north (strike=0), the normal
  // points east (+x direction).
  const double deg2rad = M_PI / 180.0;
  double s = strike_deg_ * deg2rad;
  double d = dip_deg_ * deg2rad;

  // Normal to the fault plane:
  //   n = (sin(strike)*cos(dip), -cos(strike)*cos(dip), sin(dip))
  // Wait, convention: fault normal points into the footwall.
  // For a plane with strike s and dip d:
  //   n_x =  sin(s) * cos(d)   (but actually depends on convention)
  // Simplified: normal = cross(strike_dir, dip_dir)
  // strike_dir = (sin(s), cos(s), 0)
  // dip_dir = (-cos(s)*cos(d), sin(s)*cos(d), -sin(d))
  // Simpler: for strike=0 (north), dip=90 (vertical), normal = (1, 0, 0)
  receiver_normal_[0] =  std::sin(s) * std::cos(d);
  receiver_normal_[1] = -std::cos(s) * std::cos(d);
  receiver_normal_[2] =  std::sin(d);

  // Normalize (should already be unit, but be safe)
  double mag = std::sqrt(receiver_normal_[0] * receiver_normal_[0] +
                         receiver_normal_[1] * receiver_normal_[1] +
                         receiver_normal_[2] * receiver_normal_[2]);
  if (mag > 0.0)
  {
    receiver_normal_[0] /= mag;
    receiver_normal_[1] /= mag;
    receiver_normal_[2] /= mag;
  }
}

void DerivedFieldComputer::computeStressFromStrain(const double eps[6], double P,
                                                    double sigma[6]) const
{
  double eps_kk = eps[0] + eps[1] + eps[2];
  sigma[0] = lambda_ * eps_kk + 2.0 * mu_ * eps[0] - biot_alpha_ * P;
  sigma[1] = lambda_ * eps_kk + 2.0 * mu_ * eps[1] - biot_alpha_ * P;
  sigma[2] = lambda_ * eps_kk + 2.0 * mu_ * eps[2] - biot_alpha_ * P;
  sigma[3] = 2.0 * mu_ * eps[3];
  sigma[4] = 2.0 * mu_ * eps[4];
  sigma[5] = 2.0 * mu_ * eps[5];
}

double DerivedFieldComputer::computeCFS(const double sigma[6]) const
{
  // Resolve stress onto receiver fault orientation
  // Traction on plane: t = sigma * n
  const double* s = sigma;
  const double* n = receiver_normal_.data();

  double sn[3] = {
    s[0] * n[0] + s[3] * n[1] + s[4] * n[2],
    s[3] * n[0] + s[1] * n[1] + s[5] * n[2],
    s[4] * n[0] + s[5] * n[1] + s[2] * n[2]
  };

  // Normal stress
  double sigma_n = sn[0] * n[0] + sn[1] * n[1] + sn[2] * n[2];

  // Shear stress magnitude
  double sn_mag_sq = sn[0] * sn[0] + sn[1] * sn[1] + sn[2] * sn[2];
  double tau_sq = sn_mag_sq - sigma_n * sigma_n;
  double tau = (tau_sq > 0.0) ? std::sqrt(tau_sq) : 0.0;

  // CFS = |tau| - friction * |sigma_n| (compression negative convention)
  return tau - friction_ * std::abs(sigma_n);
}

PetscErrorCode DerivedFieldComputer::compute(DM dm, Vec solution)
{
  PetscFunctionBeginUser;
  PetscErrorCode ierr;

  PetscInt dim;
  ierr = DMGetDimension(dm, &dim); CHKERRQ(ierr);

  // Get local solution
  Vec local;
  ierr = DMGetLocalVector(dm, &local); CHKERRQ(ierr);
  ierr = DMGlobalToLocal(dm, solution, INSERT_VALUES, local); CHKERRQ(ierr);
  // Insert Dirichlet BC values into constrained DOFs of the local vector.
  // DMGlobalToLocal fills constrained DOFs with zero; we need the actual BC values
  // for correct gradient computation at boundary cells.
  ierr = DMPlexInsertBoundaryValues(dm, PETSC_TRUE, local, 0.0,
                                     nullptr, nullptr, nullptr); CHKERRQ(ierr);

  PetscSection section;
  ierr = DMGetLocalSection(dm, &section); CHKERRQ(ierr);

  // Cell range
  PetscInt cStart, cEnd;
  ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd); CHKERRQ(ierr);
  num_cells_ = cEnd - cStart;

  // Identify displacement field (field with dim components)
  PetscInt nFields;
  ierr = DMGetNumFields(dm, &nFields); CHKERRQ(ierr);

  PetscInt disp_field = 0;
  PetscFE fe_disp = nullptr;
  for (PetscInt f = 0; f < nFields; ++f)
  {
    PetscObject obj;
    ierr = DMGetField(dm, f, nullptr, &obj); CHKERRQ(ierr);
    PetscFE fe = (PetscFE)obj;
    PetscInt nc;
    ierr = PetscFEGetNumComponents(fe, &nc); CHKERRQ(ierr);
    if (nc == dim)
    {
      disp_field = f;
      fe_disp = fe;
      break;
    }
  }

  if (!fe_disp)
  {
    ierr = DMRestoreLocalVector(dm, &local); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  // Tabulate basis derivatives at reference centroid
  PetscTabulation tab = nullptr;
  std::vector<PetscReal> refPoint(dim, 0.0);
  ierr = PetscFECreateTabulation(fe_disp, 1, 1, refPoint.data(), 1, &tab); CHKERRQ(ierr);

  PetscInt Nb = tab->Nb;

  // Compute displacement DOF offset (sum of all field dimensions before disp_field)
  PetscInt disp_offset = 0;
  for (PetscInt f = 0; f < disp_field; ++f)
  {
    PetscObject obj;
    ierr = DMGetField(dm, f, nullptr, &obj); CHKERRQ(ierr);
    PetscFE fe_f = (PetscFE)obj;
    PetscInt nb_f;
    ierr = PetscFEGetDimension(fe_f, &nb_f); CHKERRQ(ierr);
    disp_offset += nb_f;
  }

  // Get coordinate data for self-consistent Jacobian computation.
  // DMPlexComputeCellGeometryAffineFEM uses an internally-corrected Jacobian,
  // but the tabulation D assumes the FE reference vertex ordering. The closure
  // DOFs follow the mesh DAG vertex ordering. These can differ on hex cells,
  // causing sign errors in the gradient. By computing J from the coordinate
  // closure using the same tabulation D and same vertex ordering as the DOF
  // closure, any permutation mismatch cancels out.
  DM cdm;
  ierr = DMGetCoordinateDM(dm, &cdm); CHKERRQ(ierr);
  Vec coordsLocal;
  ierr = DMGetCoordinatesLocal(dm, &coordsLocal); CHKERRQ(ierr);
  PetscSection coordSection;
  ierr = DMGetLocalSection(cdm, &coordSection); CHKERRQ(ierr);

  // Allocate output arrays
  stress_.resize(6 * num_cells_);
  strain_.resize(6 * num_cells_);
  cfs_.resize(num_cells_);

  const PetscReal *D = tab->T[1];

  for (PetscInt c = cStart; c < cEnd; ++c)
  {
    PetscInt ci = c - cStart;

    // Get coordinate closure (same vertex ordering as DOF closure)
    PetscScalar *coordDofs = nullptr;
    PetscInt coordDofSize;
    ierr = DMPlexVecGetClosure(cdm, coordSection, coordsLocal, c,
                                &coordDofSize, &coordDofs);
    CHKERRQ(ierr);

    // Compute Jacobian from coordinate closure using tabulation D
    // J[i][j] = dx_i/dxi_j = sum_b coord_i(b) * dN_b/dxi_j
    // Using same (D, vertex ordering) as for displacement ensures consistency.
    double Jm[3][3] = {{0}};
    for (PetscInt b = 0; b < Nb; ++b)
    {
      PetscInt comp = b % dim;
      double x_val = PetscRealPart(coordDofs[b]);
      for (PetscInt j = 0; j < dim; ++j)
      {
        Jm[comp][j] += x_val * D[b * dim * dim + comp * dim + j];
      }
    }

    ierr = DMPlexVecRestoreClosure(cdm, coordSection, coordsLocal, c,
                                    &coordDofSize, &coordDofs);
    CHKERRQ(ierr);

    // Invert 3x3 Jacobian
    double det = Jm[0][0] * (Jm[1][1] * Jm[2][2] - Jm[1][2] * Jm[2][1])
               - Jm[0][1] * (Jm[1][0] * Jm[2][2] - Jm[1][2] * Jm[2][0])
               + Jm[0][2] * (Jm[1][0] * Jm[2][1] - Jm[1][1] * Jm[2][0]);
    double invDet = 1.0 / det;

    double invJm[3][3];
    invJm[0][0] =  (Jm[1][1] * Jm[2][2] - Jm[1][2] * Jm[2][1]) * invDet;
    invJm[0][1] = -(Jm[0][1] * Jm[2][2] - Jm[0][2] * Jm[2][1]) * invDet;
    invJm[0][2] =  (Jm[0][1] * Jm[1][2] - Jm[0][2] * Jm[1][1]) * invDet;
    invJm[1][0] = -(Jm[1][0] * Jm[2][2] - Jm[1][2] * Jm[2][0]) * invDet;
    invJm[1][1] =  (Jm[0][0] * Jm[2][2] - Jm[0][2] * Jm[2][0]) * invDet;
    invJm[1][2] = -(Jm[0][0] * Jm[1][2] - Jm[0][2] * Jm[1][0]) * invDet;
    invJm[2][0] =  (Jm[1][0] * Jm[2][1] - Jm[1][1] * Jm[2][0]) * invDet;
    invJm[2][1] = -(Jm[0][0] * Jm[2][1] - Jm[0][1] * Jm[2][0]) * invDet;
    invJm[2][2] =  (Jm[0][0] * Jm[1][1] - Jm[0][1] * Jm[1][0]) * invDet;

    // Get cell DOFs
    PetscScalar *cellDofs = nullptr;
    PetscInt cellDofSize;
    ierr = DMPlexVecGetClosure(dm, section, local, c, &cellDofSize, &cellDofs);
    CHKERRQ(ierr);

    if (ci == 0 || ci == 12)
    {
      PetscPrintf(PETSC_COMM_SELF, "Cell %d: Jm diag=%.4f,%.4f,%.4f det=%.6e coordSize=%d dofSize=%d\n",
        (int)c, Jm[0][0], Jm[1][1], Jm[2][2], det, (int)coordDofSize, (int)cellDofSize);
      PetscPrintf(PETSC_COMM_SELF, "  invJm diag=%.4f,%.4f,%.4f\n",
        invJm[0][0], invJm[1][1], invJm[2][2]);
    }

    // Compute displacement reference gradient: du_comp/dxi_j
    double grad_u_ref[3][3] = {{0}};
    for (PetscInt b = 0; b < Nb; ++b)
    {
      PetscInt comp = b % dim;
      double u_val = PetscRealPart(cellDofs[disp_offset + b]);
      for (PetscInt j = 0; j < dim; ++j)
      {
        grad_u_ref[comp][j] += u_val * D[b * dim * dim + comp * dim + j];
      }
    }

    ierr = DMPlexVecRestoreClosure(dm, section, local, c, &cellDofSize, &cellDofs);
    CHKERRQ(ierr);

    // Physical gradient: du_i/dx_j = sum_k (du_i/dxi_k) * invJ[k][j]
    double grad_u[3][3] = {{0}};
    for (PetscInt i = 0; i < dim; ++i)
    {
      for (PetscInt j = 0; j < dim; ++j)
      {
        for (PetscInt k = 0; k < dim; ++k)
        {
          grad_u[i][j] += grad_u_ref[i][k] * invJm[k][j];
        }
      }
    }

    // Symmetrize to strain (Voigt)
    double eps[6];
    eps[0] = grad_u[0][0];
    eps[1] = grad_u[1][1];
    eps[2] = (dim > 2) ? grad_u[2][2] : 0.0;
    eps[3] = 0.5 * (grad_u[0][1] + grad_u[1][0]);
    eps[4] = (dim > 2) ? 0.5 * (grad_u[0][2] + grad_u[2][0]) : 0.0;
    eps[5] = (dim > 2) ? 0.5 * (grad_u[1][2] + grad_u[2][1]) : 0.0;

    // Stress (no pressure coupling for pure elastostatics; P=0)
    double sigma[6];
    computeStressFromStrain(eps, 0.0, sigma);

    // Store
    for (int i = 0; i < 6; ++i)
    {
      strain_[6 * ci + i] = eps[i];
      stress_[6 * ci + i] = sigma[i];
    }

    // CFS
    cfs_[ci] = computeCFS(sigma);
  }

  ierr = PetscTabulationDestroy(&tab); CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dm, &local); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode DerivedFieldComputer::writeStressHDF5(MPI_Comm comm, const char* filename,
                                                      PetscInt step, PetscReal time)
{
  PetscFunctionBeginUser;
  PetscErrorCode ierr;

  if (num_cells_ == 0) PetscFunctionReturn(0);

  // Create a simple sequential Vec for cell-centered stress (6 components)
  Vec stress_vec;
  ierr = VecCreateSeq(PETSC_COMM_SELF, 6 * num_cells_, &stress_vec); CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)stress_vec, "stress"); CHKERRQ(ierr);

  PetscScalar *arr;
  ierr = VecGetArray(stress_vec, &arr); CHKERRQ(ierr);
  for (PetscInt i = 0; i < 6 * num_cells_; ++i)
  {
    arr[i] = stress_[i];
  }
  ierr = VecRestoreArray(stress_vec, &arr); CHKERRQ(ierr);

  // Write to HDF5
  PetscViewer viewer;
  ierr = PetscViewerHDF5Open(comm, filename, FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
  ierr = VecView(stress_vec, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  ierr = VecDestroy(&stress_vec); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode DerivedFieldComputer::writeCfsHDF5(MPI_Comm comm, const char* filename,
                                                    PetscInt step, PetscReal time)
{
  PetscFunctionBeginUser;
  PetscErrorCode ierr;

  if (num_cells_ == 0) PetscFunctionReturn(0);

  Vec cfs_vec;
  ierr = VecCreateSeq(PETSC_COMM_SELF, num_cells_, &cfs_vec); CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)cfs_vec, "cfs"); CHKERRQ(ierr);

  PetscScalar *arr;
  ierr = VecGetArray(cfs_vec, &arr); CHKERRQ(ierr);
  for (PetscInt i = 0; i < num_cells_; ++i)
  {
    arr[i] = cfs_[i];
  }
  ierr = VecRestoreArray(cfs_vec, &arr); CHKERRQ(ierr);

  PetscViewer viewer;
  ierr = PetscViewerHDF5Open(comm, filename, FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
  ierr = VecView(cfs_vec, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  ierr = VecDestroy(&cfs_vec); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

} // namespace FSRM
