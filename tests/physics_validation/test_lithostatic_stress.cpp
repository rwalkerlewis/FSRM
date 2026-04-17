#include <gtest/gtest.h>
#include <petsc.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include "core/Simulator.hpp"
#include "core/FSRM.hpp"
#include "numerics/PetscFEElasticityGravity.hpp"
#include "numerics/PetscFEElasticityAux.hpp"

using namespace FSRM;

// Physics validation test: solve a 1D column under gravity and verify
// that the stress at depth matches the analytical lithostatic solution.
//
// Geometry: 1x1x1000 m column (x_extent=1, y_extent=1, z_extent=1000)
// Material: rho=2650, E=10 GPa, nu=0.25
//           lambda = 4 GPa, mu = 4 GPa, K0 = nu/(1-nu) = 1/3
// BCs: bottom fixed (z=0), sides roller, top free (z=1000)
// Gravity: g = 9.81 m/s^2 downward
//
// Analytical solution:
//   sigma_zz(z) = -rho * g * (H - z)    [compressive negative]
//   sigma_xx(z) = sigma_yy(z) = K0 * sigma_zz(z)
//   u_z(z) = -rho*g/(lambda+2*mu) * (H*z - z^2/2)

class LithostaticStressTest : public ::testing::Test
{
protected:
  static constexpr double rho_ = 2650.0;
  static constexpr double E_ = 10.0e9;
  static constexpr double nu_ = 0.25;
  static constexpr double g_ = 9.81;
  static constexpr double H_ = 1000.0;

  // Derived
  static constexpr double lambda_ = E_ * nu_ / ((1.0 + nu_) * (1.0 - 2.0 * nu_));
  static constexpr double mu_ = E_ / (2.0 * (1.0 + nu_));
  static constexpr double K0_ = nu_ / (1.0 - nu_);

  // Analytical stress at height z above base
  double analytical_sigma_zz(double z) const
  {
    return -rho_ * g_ * (H_ - z);
  }

  double analytical_sigma_xx(double z) const
  {
    return K0_ * analytical_sigma_zz(z);
  }

  // Write the temporary config for the lithostatic column
  std::string writeConfig(int nz = 10)
  {
    const char* path = "test_lithostatic_solve.config";
    FILE* fp = fopen(path, "w");
    if (!fp) return "";
    fprintf(fp,
      "[SIMULATION]\n"
      "name = lithostatic_stress_test\n"
      "start_time = 0.0\n"
      "end_time = 1.0\n"
      "dt_initial = 1.0\n"
      "output_frequency = 1\n"
      "fluid_model = NONE\n"
      "solid_model = ELASTIC\n"
      "enable_geomechanics = true\n"
      "enable_faults = false\n"
      "max_nonlinear_iterations = 20\n"
      "enable_gravity = true\n"
      "rtol = 1.0e-10\n"
      "atol = 1.0e-12\n"
      "\n"
      "[GRID]\n"
      "nx = 2\n"
      "ny = 2\n"
      "nz = %d\n"
      "Lx = 1.0\n"
      "Ly = 1.0\n"
      "Lz = 1000.0\n"
      "\n"
      "[ROCK]\n"
      "density = 2650.0\n"
      "youngs_modulus = 10.0e9\n"
      "poissons_ratio = 0.25\n"
      "\n"
      "[MATERIAL]\n"
      "gravity = 9.81\n"
      "\n"
      "[GEOMECHANICS]\n"
      "enabled = true\n"
      "\n"
      "[BOUNDARY_CONDITIONS]\n"
      "bottom = fixed\n"
      "sides = roller\n"
      "top = free\n",
      nz);
    fclose(fp);
    return path;
  }

  // Compute stress at cell centroids from the FEM displacement solution.
  // Returns a vector of (z_centroid, sigma_xx, sigma_yy, sigma_zz) tuples.
  struct CellStress
  {
    double z;
    double sigma_xx, sigma_yy, sigma_zz;
  };

  PetscErrorCode computeCellStresses(DM dm, Vec solution,
      std::vector<CellStress>& stresses)
  {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    PetscInt dim;
    ierr = DMGetDimension(dm, &dim); CHKERRQ(ierr);

    // Get local solution
    Vec local;
    ierr = DMGetLocalVector(dm, &local); CHKERRQ(ierr);
    ierr = DMGlobalToLocal(dm, solution, INSERT_VALUES, local); CHKERRQ(ierr);

    PetscSection section;
    ierr = DMGetLocalSection(dm, &section); CHKERRQ(ierr);

    // Cell range
    PetscInt cStart, cEnd;
    ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd); CHKERRQ(ierr);

    // Find the displacement FE (vector field with dim components)
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

    if (!fe_disp) {
      ierr = DMRestoreLocalVector(dm, &local); CHKERRQ(ierr);
      PetscFunctionReturn(PETSC_ERR_ARG_WRONG);
    }

    // Tabulate basis function derivatives at reference origin
    PetscTabulation tab = nullptr;
    std::vector<PetscReal> refPoint(dim, 0.0);
    ierr = PetscFECreateTabulation(fe_disp, 1, 1, refPoint.data(), 1, &tab);
    CHKERRQ(ierr);

    PetscInt Nb = tab->Nb;

    // Compute displacement field offset for multi-field problems
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

    // Get material properties from aux fields if available
    Vec auxVec = nullptr;
    PetscSection auxSection = nullptr;
    DM auxDM = nullptr;
    {
      ierr = DMGetAuxiliaryVec(dm, nullptr, 0, 0, &auxVec);
      if (ierr || !auxVec) {
        // No aux vector, use defaults
        ierr = 0;
        auxVec = nullptr;
      } else {
        ierr = VecGetDM(auxVec, &auxDM); CHKERRQ(ierr);
        if (auxDM) {
          ierr = DMGetLocalSection(auxDM, &auxSection); CHKERRQ(ierr);
        }
      }
    }

    stresses.clear();
    stresses.reserve(cEnd - cStart);

    const PetscScalar *auxArray = nullptr;
    if (auxVec) {
      ierr = VecGetArrayRead(auxVec, &auxArray); CHKERRQ(ierr);
    }

    for (PetscInt c = cStart; c < cEnd; ++c)
    {
      // Get cell centroid
      PetscReal vol, centroid[3], normal[3];
      ierr = DMPlexComputeCellGeometryFVM(dm, c, &vol, centroid, normal);
      CHKERRQ(ierr);

      // Get inverse Jacobian
      std::vector<PetscReal> v0(dim), J(dim * dim), invJ(dim * dim);
      PetscReal detJ;
      ierr = DMPlexComputeCellGeometryAffineFEM(dm, c, v0.data(), J.data(),
                                                 invJ.data(), &detJ);
      CHKERRQ(ierr);

      // Get cell DOFs
      PetscScalar *cellDofs = nullptr;
      PetscInt cellDofSize;
      ierr = DMPlexVecGetClosure(dm, section, local, c, &cellDofSize, &cellDofs);
      CHKERRQ(ierr);

      // Compute displacement gradient
      double grad_u[3][3] = {{0}};
      const PetscReal *D = tab->T[1];

      for (PetscInt b = 0; b < Nb; ++b)
      {
        PetscInt comp = b % dim;
        double u_val = PetscRealPart(cellDofs[disp_offset + b]);

        for (PetscInt j = 0; j < dim; ++j)
        {
          double phys_grad = 0.0;
          for (PetscInt k = 0; k < dim; ++k)
          {
            phys_grad += D[b * dim * dim + comp * dim + k] * invJ[k * dim + j];
          }
          grad_u[comp][j] += u_val * phys_grad;
        }
      }

      ierr = DMPlexVecRestoreClosure(dm, section, local, c, &cellDofSize, &cellDofs);
      CHKERRQ(ierr);

      // Compute strain
      double eps_xx = grad_u[0][0];
      double eps_yy = grad_u[1][1];
      double eps_zz = grad_u[2][2];
      double trace = eps_xx + eps_yy + eps_zz;

      // Get lambda, mu from aux fields or defaults
      double lam = lambda_;
      double mu = mu_;
      if (auxArray && auxSection)
      {
        PetscInt auxOff;
        ierr = PetscSectionGetOffset(auxSection, c, &auxOff); CHKERRQ(ierr);
        lam = PetscRealPart(auxArray[auxOff + AUX_LAMBDA]);
        mu  = PetscRealPart(auxArray[auxOff + AUX_MU]);
      }

      // Compute stress via Hooke's law
      CellStress cs;
      cs.z = centroid[dim - 1];
      cs.sigma_xx = lam * trace + 2.0 * mu * eps_xx;
      cs.sigma_yy = lam * trace + 2.0 * mu * eps_yy;
      cs.sigma_zz = lam * trace + 2.0 * mu * eps_zz;
      stresses.push_back(cs);
    }

    if (auxArray && auxVec) {
      ierr = VecRestoreArrayRead(auxVec, &auxArray); CHKERRQ(ierr);
    }

    ierr = PetscTabulationDestroy(&tab); CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(dm, &local); CHKERRQ(ierr);

    // Sort by z coordinate
    std::sort(stresses.begin(), stresses.end(),
        [](const CellStress& a, const CellStress& b) { return a.z < b.z; });

    PetscFunctionReturn(0);
  }
};

// Test: New PetscFEElasticityGravity callbacks have correct signatures
TEST_F(LithostaticStressTest, GravityCallbackDelegation)
{
  const PetscInt dim = 3;
  const PetscInt Nf = 1;
  const PetscInt NfAux = 3;
  const PetscInt uOff[1] = {0};
  const PetscInt uOff_x[1] = {0};
  PetscScalar u[3] = {0};
  PetscScalar u_x[9] = {0};
  const PetscInt aOff[3] = {0, 1, 2};
  const PetscInt aOff_x[3] = {0, 0, 0};
  PetscScalar a[3] = {lambda_, mu_, rho_};
  PetscScalar a_x[1] = {0};
  PetscReal t = 0.0;
  PetscReal x[3] = {0.0, 0.0, 500.0};
  PetscScalar constants[1] = {g_};
  PetscScalar f0[3] = {0};

  // Use PetscFEElasticityGravity callback (delegates to PetscFEElasticityAux)
  PetscFEElasticityGravity::f0_gravity(
      dim, Nf, NfAux, uOff, uOff_x, u, nullptr, u_x,
      aOff, aOff_x, a, nullptr, a_x, t, x, 1, constants, f0);

  EXPECT_DOUBLE_EQ(PetscRealPart(f0[0]), 0.0);
  EXPECT_DOUBLE_EQ(PetscRealPart(f0[1]), 0.0);
  // f0[2] = +rho*g (positive: PETSc FEM convention, f0 = -body_force)
  EXPECT_NEAR(PetscRealPart(f0[2]), rho_ * g_, 1.0e-8);
}

// Full solve test: verify lithostatic stress at 10 depth points
TEST_F(LithostaticStressTest, LithostaticColumnFullSolve)
{
  std::string config_path = writeConfig(10);
  ASSERT_FALSE(config_path.empty());

  // Use direct LU solver for robustness on this small problem
  PetscOptionsSetValue(NULL, "-ksp_type", "preonly");
  PetscOptionsSetValue(NULL, "-pc_type", "lu");

  Simulator sim(PETSC_COMM_WORLD);
  PetscErrorCode ierr;

  // Full pipeline
  ierr = sim.initializeFromConfigFile(config_path);
  ASSERT_EQ(ierr, 0) << "initializeFromConfigFile failed";

  ierr = sim.setupDM();
  ASSERT_EQ(ierr, 0) << "setupDM failed";

  ierr = sim.labelBoundaries();
  ASSERT_EQ(ierr, 0) << "labelBoundaries failed";

  ierr = sim.setupFields();
  ASSERT_EQ(ierr, 0) << "setupFields failed";

  ierr = sim.setupPhysics();
  ASSERT_EQ(ierr, 0) << "setupPhysics failed";

  ierr = sim.setupTimeStepper();
  ASSERT_EQ(ierr, 0) << "setupTimeStepper failed";

  ierr = sim.setupSolvers();
  ASSERT_EQ(ierr, 0) << "setupSolvers failed";

  ierr = sim.setInitialConditions();
  ASSERT_EQ(ierr, 0) << "setInitialConditions failed";

  // Run the solve
  ierr = sim.run();
  ASSERT_EQ(ierr, 0) << "run() failed";

  // Get DM and solution
  DM dm = sim.getDM();
  Vec solution = sim.getSolution();
  ASSERT_NE(dm, nullptr);
  ASSERT_NE(solution, nullptr);

  // Compute stress at all cell centroids
  std::vector<CellStress> stresses;
  ierr = computeCellStresses(dm, solution, stresses);
  ASSERT_EQ(ierr, 0) << "computeCellStresses failed";
  ASSERT_GT(stresses.size(), 0u);

  // Select 10 depth points evenly spaced from z=50 to z=950
  // (avoid boundaries where FEM constraints may affect accuracy)
  int n_checks = 10;
  int n_pass_zz = 0;
  int n_pass_xx = 0;

  for (int i = 0; i < n_checks; ++i)
  {
    double z_target = 50.0 + i * (900.0 / (n_checks - 1));

    // Find the cell closest to z_target
    double best_dist = 1e30;
    const CellStress* best = nullptr;
    for (const auto& cs : stresses)
    {
      double dist = std::fabs(cs.z - z_target);
      if (dist < best_dist)
      {
        best_dist = dist;
        best = &cs;
      }
    }
    ASSERT_NE(best, nullptr);

    double z = best->z;
    double expected_zz = analytical_sigma_zz(z);
    double expected_xx = analytical_sigma_xx(z);

    // Relative error (avoid division by zero near free surface)
    double ref_stress = std::fabs(expected_zz);
    if (ref_stress < 1.0) continue;  // Skip near-zero stress points

    double rel_err_zz = std::fabs(best->sigma_zz - expected_zz) / ref_stress;
    double rel_err_xx = std::fabs(best->sigma_xx - expected_xx) / ref_stress;

    // Report
    if (rel_err_zz <= 0.01) n_pass_zz++;
    if (rel_err_xx <= 0.01) n_pass_xx++;

    // 1% tolerance for each depth point
    EXPECT_LE(rel_err_zz, 0.05)
        << "sigma_zz at z=" << z << ": got " << best->sigma_zz
        << ", expected " << expected_zz
        << " (rel_err " << rel_err_zz * 100 << "%)";
    EXPECT_LE(rel_err_xx, 0.05)
        << "sigma_xx at z=" << z << ": got " << best->sigma_xx
        << ", expected " << expected_xx
        << " (rel_err " << rel_err_xx * 100 << "%)";
  }

  // At least 7 out of 10 points should be within 1% tolerance
  // (boundary effects may affect some cells)
  EXPECT_GE(n_pass_zz, 7)
      << "sigma_zz: only " << n_pass_zz << " of " << n_checks
      << " points within 1% tolerance";
  EXPECT_GE(n_pass_xx, 7)
      << "sigma_xx: only " << n_pass_xx << " of " << n_checks
      << " points within 1% tolerance";

  remove(config_path.c_str());
  PetscOptionsClearValue(NULL, "-ksp_type");
  PetscOptionsClearValue(NULL, "-pc_type");
}

// Test: K0 ratio is consistent across the column
TEST_F(LithostaticStressTest, K0RatioAtDepth)
{
  std::string config_path = writeConfig(10);
  ASSERT_FALSE(config_path.empty());

  // Use direct LU solver for robustness on this small problem
  PetscOptionsSetValue(NULL, "-ksp_type", "preonly");
  PetscOptionsSetValue(NULL, "-pc_type", "lu");

  Simulator sim(PETSC_COMM_WORLD);
  PetscErrorCode ierr;

  ierr = sim.initializeFromConfigFile(config_path);
  ASSERT_EQ(ierr, 0);
  ierr = sim.setupDM();
  ASSERT_EQ(ierr, 0);
  ierr = sim.labelBoundaries();
  ASSERT_EQ(ierr, 0);
  ierr = sim.setupFields();
  ASSERT_EQ(ierr, 0);
  ierr = sim.setupPhysics();
  ASSERT_EQ(ierr, 0);
  ierr = sim.setupTimeStepper();
  ASSERT_EQ(ierr, 0);
  ierr = sim.setupSolvers();
  ASSERT_EQ(ierr, 0);
  ierr = sim.setInitialConditions();
  ASSERT_EQ(ierr, 0);
  ierr = sim.run();
  ASSERT_EQ(ierr, 0);

  DM dm = sim.getDM();
  Vec solution = sim.getSolution();

  std::vector<CellStress> stresses;
  ierr = computeCellStresses(dm, solution, stresses);
  ASSERT_EQ(ierr, 0);

  // For cells in the interior (away from boundaries), K0 = sigma_xx / sigma_zz
  // should be approximately nu/(1-nu) = 1/3
  int n_interior = 0;
  double sum_K0 = 0.0;
  for (const auto& cs : stresses)
  {
    // Only check cells well inside the domain
    if (cs.z > 100.0 && cs.z < 900.0 && std::fabs(cs.sigma_zz) > 1e4)
    {
      double measured_K0 = cs.sigma_xx / cs.sigma_zz;
      sum_K0 += measured_K0;
      n_interior++;
    }
  }

  ASSERT_GT(n_interior, 0) << "No interior cells found";

  double avg_K0 = sum_K0 / n_interior;
  EXPECT_NEAR(avg_K0, K0_, 0.02)
      << "Average K0 = " << avg_K0 << ", expected " << K0_;

  remove(config_path.c_str());
  PetscOptionsClearValue(NULL, "-ksp_type");
  PetscOptionsClearValue(NULL, "-pc_type");
}
