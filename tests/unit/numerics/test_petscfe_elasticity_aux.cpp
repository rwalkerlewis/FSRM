/**
 * @file test_petscfe_elasticity_aux.cpp
 * @brief Unit tests for PetscFEElasticityAux callbacks
 *
 * Verifies that auxiliary-field-aware elasticity callbacks produce correct
 * stress and stiffness tensor values for known inputs.
 */

#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include <petsc.h>

#include "numerics/PetscFEElasticityAux.hpp"

using namespace FSRM;

class PetscFEElasticityAuxTest : public ::testing::Test
{
protected:
  static constexpr PetscInt DIM = 3;

  // Aux field values: granite
  static constexpr PetscScalar TEST_LAMBDA = 30.0e9;  // Pa
  static constexpr PetscScalar TEST_MU     = 25.0e9;  // Pa
  static constexpr PetscScalar TEST_RHO    = 2650.0;   // kg/m^3
};

TEST_F(PetscFEElasticityAuxTest, F1StressFromUniaxialStrain)
{
  // Apply uniaxial strain: eps_zz = 0.001, all other eps = 0
  // Expected: sigma_zz = (lambda + 2*mu) * 0.001 = (30e9 + 50e9)*0.001 = 80 MPa
  // sigma_xx = sigma_yy = lambda * 0.001 = 30 MPa

  // Set up aux fields: [lambda, mu, rho]
  PetscScalar aux_vals[3] = {TEST_LAMBDA, TEST_MU, TEST_RHO};
  PetscInt aOff[3] = {0, 1, 2};
  PetscInt aOff_x[3] = {0, 0, 0};

  // Set up solution gradient: du_z/dz = 0.001, all others = 0
  // u_x layout: [du_x/dx, du_x/dy, du_x/dz, du_y/dx, du_y/dy, du_y/dz, du_z/dx, du_z/dy, du_z/dz]
  PetscScalar u_x[9] = {0};
  u_x[8] = 0.001;  // du_z/dz = eps_zz

  PetscInt uOff[1] = {0};
  PetscInt uOff_x[1] = {0};
  PetscScalar u[3] = {0};

  PetscScalar f1[9] = {0};
  PetscReal x[3] = {0};

  PetscFEElasticityAux::f1_elastostatics_aux(
      DIM, 1, NUM_AUX_FIELDS,
      uOff, uOff_x, u, nullptr, u_x,
      aOff, aOff_x, aux_vals, nullptr, nullptr,
      0.0, x, 0, nullptr, f1);

  // Check stress values
  // sigma_xx = lambda * eps_zz = 30e9 * 0.001 = 30e6
  EXPECT_NEAR(PetscRealPart(f1[0]), 30.0e6, 1.0e3);  // sigma_xx
  // sigma_yy = lambda * eps_zz = 30e6
  EXPECT_NEAR(PetscRealPart(f1[4]), 30.0e6, 1.0e3);  // sigma_yy
  // sigma_zz = (lambda + 2*mu) * eps_zz = 80e9 * 0.001 = 80e6
  EXPECT_NEAR(PetscRealPart(f1[8]), 80.0e6, 1.0e3);  // sigma_zz

  // Off-diagonal should be zero
  EXPECT_NEAR(PetscRealPart(f1[1]), 0.0, 1.0e-6);  // sigma_xy
  EXPECT_NEAR(PetscRealPart(f1[2]), 0.0, 1.0e-6);  // sigma_xz
  EXPECT_NEAR(PetscRealPart(f1[3]), 0.0, 1.0e-6);  // sigma_yx
  EXPECT_NEAR(PetscRealPart(f1[5]), 0.0, 1.0e-6);  // sigma_yz
  EXPECT_NEAR(PetscRealPart(f1[6]), 0.0, 1.0e-6);  // sigma_zx
  EXPECT_NEAR(PetscRealPart(f1[7]), 0.0, 1.0e-6);  // sigma_zy
}

TEST_F(PetscFEElasticityAuxTest, F1StressFromShearStrain)
{
  // Apply pure shear: eps_xz = 0.001 (from du_x/dz = 0.001, du_z/dx = 0.001)
  // Expected: sigma_xz = 2*mu*eps_xz = 2*25e9*0.001 = 50 MPa

  PetscScalar aux_vals[3] = {TEST_LAMBDA, TEST_MU, TEST_RHO};
  PetscInt aOff[3] = {0, 1, 2};
  PetscInt aOff_x[3] = {0, 0, 0};

  PetscScalar u_x[9] = {0};
  u_x[2] = 0.001;  // du_x/dz
  u_x[6] = 0.001;  // du_z/dx

  PetscInt uOff[1] = {0};
  PetscInt uOff_x[1] = {0};
  PetscScalar u[3] = {0};

  PetscScalar f1[9] = {0};
  PetscReal x[3] = {0};

  PetscFEElasticityAux::f1_elastostatics_aux(
      DIM, 1, NUM_AUX_FIELDS,
      uOff, uOff_x, u, nullptr, u_x,
      aOff, aOff_x, aux_vals, nullptr, nullptr,
      0.0, x, 0, nullptr, f1);

  // sigma_xz = 2*mu*eps_xz = 2*25e9*0.001 = 50e6
  EXPECT_NEAR(PetscRealPart(f1[2]), 50.0e6, 1.0e3);  // sigma_xz
  EXPECT_NEAR(PetscRealPart(f1[6]), 50.0e6, 1.0e3);  // sigma_zx (symmetric)

  // Diagonal stresses should be zero (no volumetric strain)
  EXPECT_NEAR(PetscRealPart(f1[0]), 0.0, 1.0e-6);  // sigma_xx
  EXPECT_NEAR(PetscRealPart(f1[4]), 0.0, 1.0e-6);  // sigma_yy
  EXPECT_NEAR(PetscRealPart(f1[8]), 0.0, 1.0e-6);  // sigma_zz
}

TEST_F(PetscFEElasticityAuxTest, G3ElasticityTensorNonzeroEntries)
{
  PetscScalar aux_vals[3] = {TEST_LAMBDA, TEST_MU, TEST_RHO};
  PetscInt aOff[3] = {0, 1, 2};
  PetscInt aOff_x[3] = {0, 0, 0};

  PetscScalar g3[81] = {0};  // 3^4 = 81
  PetscScalar u_x[9] = {0};
  PetscInt uOff[1] = {0};
  PetscInt uOff_x[1] = {0};
  PetscScalar u[3] = {0};
  PetscReal x[3] = {0};

  PetscFEElasticityAux::g3_elastostatics_aux(
      DIM, 1, NUM_AUX_FIELDS,
      uOff, uOff_x, u, nullptr, u_x,
      aOff, aOff_x, aux_vals, nullptr, nullptr,
      0.0, 0.0, x, 0, nullptr, g3);

  // C_0000 = lambda + 2*mu = 80e9
  EXPECT_NEAR(PetscRealPart(g3[0]), 80.0e9, 1.0e3);
  // C_1111 = lambda + 2*mu = 80e9
  EXPECT_NEAR(PetscRealPart(g3[40]), 80.0e9, 1.0e3);
  // C_2222 = lambda + 2*mu = 80e9
  EXPECT_NEAR(PetscRealPart(g3[80]), 80.0e9, 1.0e3);

  // C_0011 = lambda = 30e9
  EXPECT_NEAR(PetscRealPart(g3[4]), TEST_LAMBDA, 1.0e3);
  // C_0022 = lambda = 30e9
  EXPECT_NEAR(PetscRealPart(g3[8]), TEST_LAMBDA, 1.0e3);
  // C_1100 = lambda = 30e9
  EXPECT_NEAR(PetscRealPart(g3[36]), TEST_LAMBDA, 1.0e3);

  // C_0101 = mu = 25e9 (index: ((0*3+1)*3+0)*3+1 = 10)
  EXPECT_NEAR(PetscRealPart(g3[10]), TEST_MU, 1.0e3);
  // C_0110 = mu = 25e9 (index: ((0*3+1)*3+1)*3+0 = 12)
  EXPECT_NEAR(PetscRealPart(g3[12]), TEST_MU, 1.0e3);
}

TEST_F(PetscFEElasticityAuxTest, F0GravityBodyForce)
{
  PetscScalar aux_vals[3] = {TEST_LAMBDA, TEST_MU, TEST_RHO};
  PetscInt aOff[3] = {0, 1, 2};
  PetscInt aOff_x[3] = {0, 0, 0};

  PetscScalar f0[3] = {0};
  PetscScalar u_x[9] = {0};
  PetscInt uOff[1] = {0};
  PetscInt uOff_x[1] = {0};
  PetscScalar u[3] = {0};
  PetscReal x[3] = {0};

  // Test with gravity = 9.81
  PetscScalar constants[1] = {9.81};

  PetscFEElasticityAux::f0_elastostatics_aux(
      DIM, 1, NUM_AUX_FIELDS,
      uOff, uOff_x, u, nullptr, u_x,
      aOff, aOff_x, aux_vals, nullptr, nullptr,
      0.0, x, 1, constants, f0);

  // f0[0] = f0[1] = 0
  EXPECT_NEAR(PetscRealPart(f0[0]), 0.0, 1.0e-6);
  EXPECT_NEAR(PetscRealPart(f0[1]), 0.0, 1.0e-6);
  // f0[2] = +rho * g = +2650 * 9.81 = +25996.5
  // (positive: PETSc FEM convention where f0 = -body_force)
  EXPECT_NEAR(PetscRealPart(f0[2]), 25996.5, 1.0);
}

TEST_F(PetscFEElasticityAuxTest, F0NoGravity)
{
  PetscScalar aux_vals[3] = {TEST_LAMBDA, TEST_MU, TEST_RHO};
  PetscInt aOff[3] = {0, 1, 2};
  PetscInt aOff_x[3] = {0, 0, 0};

  PetscScalar f0[3] = {0};
  PetscScalar u_x[9] = {0};
  PetscInt uOff[1] = {0};
  PetscInt uOff_x[1] = {0};
  PetscScalar u[3] = {0};
  PetscReal x[3] = {0};

  // Test with no constants (gravity disabled)
  PetscFEElasticityAux::f0_elastostatics_aux(
      DIM, 1, NUM_AUX_FIELDS,
      uOff, uOff_x, u, nullptr, u_x,
      aOff, aOff_x, aux_vals, nullptr, nullptr,
      0.0, x, 0, nullptr, f0);

  EXPECT_NEAR(PetscRealPart(f0[0]), 0.0, 1.0e-6);
  EXPECT_NEAR(PetscRealPart(f0[1]), 0.0, 1.0e-6);
  EXPECT_NEAR(PetscRealPart(f0[2]), 0.0, 1.0e-6);
}

TEST_F(PetscFEElasticityAuxTest, DifferentMaterialProperties)
{
  // Test with alluvium properties (soft material)
  PetscScalar lambda_alluvium = 3.0e9;
  PetscScalar mu_alluvium = 1.5e9;
  PetscScalar rho_alluvium = 1800.0;

  PetscScalar aux_vals[3] = {lambda_alluvium, mu_alluvium, rho_alluvium};
  PetscInt aOff[3] = {0, 1, 2};
  PetscInt aOff_x[3] = {0, 0, 0};

  // Apply same uniaxial strain as first test
  PetscScalar u_x[9] = {0};
  u_x[8] = 0.001;  // du_z/dz = eps_zz

  PetscInt uOff[1] = {0};
  PetscInt uOff_x[1] = {0};
  PetscScalar u[3] = {0};
  PetscScalar f1[9] = {0};
  PetscReal x[3] = {0};

  PetscFEElasticityAux::f1_elastostatics_aux(
      DIM, 1, NUM_AUX_FIELDS,
      uOff, uOff_x, u, nullptr, u_x,
      aOff, aOff_x, aux_vals, nullptr, nullptr,
      0.0, x, 0, nullptr, f1);

  // sigma_zz = (3e9 + 2*1.5e9)*0.001 = 6e9*0.001 = 6e6
  EXPECT_NEAR(PetscRealPart(f1[8]), 6.0e6, 1.0e3);
  // sigma_xx = lambda * eps_zz = 3e9*0.001 = 3e6
  EXPECT_NEAR(PetscRealPart(f1[0]), 3.0e6, 1.0e3);
}
