#include <gtest/gtest.h>
#include <petsc.h>
#include <cmath>
#include "core/Simulator.hpp"
#include "core/FSRM.hpp"
#include "numerics/PetscFEElasticityAux.hpp"

using namespace FSRM;

// Test fixture for gravity callback tests
class GravityLithostaticTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
  }
  void TearDown() override
  {
  }
};

// Test that the aux callback correctly computes gravity body force
TEST_F(GravityLithostaticTest, GravityBodyForceCallback)
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
  PetscScalar a[3] = {30.0e9, 25.0e9, 2650.0};  // lambda, mu, rho
  PetscScalar a_x[1] = {0};
  PetscReal t = 0.0;
  PetscReal x[3] = {0.0, 0.0, 0.0};
  PetscScalar constants[1] = {9.81};  // gravity magnitude
  PetscScalar f0[3] = {0};

  PetscFEElasticityAux::f0_elastostatics_aux(
      dim, Nf, NfAux, uOff, uOff_x, u, nullptr, u_x,
      aOff, aOff_x, a, nullptr, a_x, t, x, 1, constants, f0);

  // f0[0] and f0[1] should be zero (no lateral gravity)
  EXPECT_DOUBLE_EQ(PetscRealPart(f0[0]), 0.0);
  EXPECT_DOUBLE_EQ(PetscRealPart(f0[1]), 0.0);
  // f0[2] = rho * g = 2650 * 9.81 = 25996.5
  // (positive: in PETSc FEM convention, f0 = -body_force = rho*g*ez)
  double expected_f0z = 2650.0 * 9.81;
  EXPECT_NEAR(PetscRealPart(f0[2]), expected_f0z, 1.0e-8);
}

// Test that zero gravity produces zero body force
TEST_F(GravityLithostaticTest, ZeroGravityZeroBodyForce)
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
  PetscScalar a[3] = {30.0e9, 25.0e9, 2650.0};
  PetscScalar a_x[1] = {0};
  PetscReal t = 0.0;
  PetscReal x[3] = {0.0, 0.0, 0.0};
  PetscScalar constants[1] = {0.0};  // no gravity
  PetscScalar f0[3] = {0};

  PetscFEElasticityAux::f0_elastostatics_aux(
      dim, Nf, NfAux, uOff, uOff_x, u, nullptr, u_x,
      aOff, aOff_x, a, nullptr, a_x, t, x, 1, constants, f0);

  EXPECT_DOUBLE_EQ(PetscRealPart(f0[0]), 0.0);
  EXPECT_DOUBLE_EQ(PetscRealPart(f0[1]), 0.0);
  EXPECT_DOUBLE_EQ(PetscRealPart(f0[2]), 0.0);
}

// Test that gravity body force scales with density
TEST_F(GravityLithostaticTest, GravityScalesWithDensity)
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
  PetscScalar a_x[1] = {0};
  PetscReal t = 0.0;
  PetscReal x[3] = {0.0, 0.0, 0.0};
  PetscScalar constants[1] = {9.81};

  // Alluvium: rho = 1800
  PetscScalar a1[3] = {3.0e9, 1.5e9, 1800.0};
  PetscScalar f0_1[3] = {0};
  PetscFEElasticityAux::f0_elastostatics_aux(
      dim, Nf, NfAux, uOff, uOff_x, u, nullptr, u_x,
      aOff, aOff_x, a1, nullptr, a_x, t, x, 1, constants, f0_1);

  // Granite: rho = 2650
  PetscScalar a2[3] = {30.0e9, 25.0e9, 2650.0};
  PetscScalar f0_2[3] = {0};
  PetscFEElasticityAux::f0_elastostatics_aux(
      dim, Nf, NfAux, uOff, uOff_x, u, nullptr, u_x,
      aOff, aOff_x, a2, nullptr, a_x, t, x, 1, constants, f0_2);

  // Ratio of body forces should equal ratio of densities
  double ratio = PetscRealPart(f0_2[2]) / PetscRealPart(f0_1[2]);
  EXPECT_NEAR(ratio, 2650.0 / 1800.0, 1.0e-10);
}

// Test fixture for Simulator-level lithostatic tests
class LithostaticSimulatorTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
  }
  void TearDown() override
  {
  }
};

// Test that the simulator setup runs without crash when gravity is enabled
// with heterogeneous material, using the full pipeline
TEST_F(LithostaticSimulatorTest, GravitySetupDoesNotCrash)
{
  // Write a temporary config file
  const char* config_path = "test_gravity_setup.config";
  FILE* fp = fopen(config_path, "w");
  ASSERT_NE(fp, nullptr);
  fprintf(fp,
    "[SIMULATION]\n"
    "name = gravity_test\n"
    "dt = 1.0\n"
    "total_time = 1.0\n"
    "max_time_steps = 1\n"
    "max_nonlinear_iterations = 5\n"
    "\n"
    "[MESH]\n"
    "type = box\n"
    "x_extent = 1000.0\n"
    "y_extent = 1000.0\n"
    "z_extent = 1000.0\n"
    "nx = 3\n"
    "ny = 3\n"
    "nz = 3\n"
    "\n"
    "[ROCK]\n"
    "density = 2650.0\n"
    "youngs_modulus = 80.0e9\n"
    "poissons_ratio = 0.25\n"
    "\n"
    "[MATERIAL]\n"
    "heterogeneous = true\n"
    "gravity = 9.81\n"
    "\n"
    "[LAYER_1]\n"
    "z_top = 0.0\n"
    "z_bottom = -1000.0\n"
    "lambda = 30.0e9\n"
    "mu = 25.0e9\n"
    "rho = 2650.0\n"
    "\n"
    "[GEOMECHANICS]\n"
    "enabled = true\n"
    "\n"
    "[BOUNDARY_CONDITIONS]\n"
    "bottom = fixed\n"
    "sides = roller\n"
    "top = free\n"
  );
  fclose(fp);

  Simulator sim(PETSC_COMM_WORLD);
  PetscErrorCode ierr;

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

  remove(config_path);
}

// Test that the simulator setup runs without crash when gravity is enabled
// with homogeneous material (uses aux callbacks with uniform field values)
TEST_F(LithostaticSimulatorTest, HomogeneousGravitySetupDoesNotCrash)
{
  const char* config_path = "test_homogeneous_gravity.config";
  FILE* fp = fopen(config_path, "w");
  ASSERT_NE(fp, nullptr);
  fprintf(fp,
    "[SIMULATION]\n"
    "name = homogeneous_gravity_test\n"
    "dt = 1.0\n"
    "total_time = 1.0\n"
    "max_time_steps = 1\n"
    "max_nonlinear_iterations = 5\n"
    "\n"
    "[MESH]\n"
    "type = box\n"
    "x_extent = 1000.0\n"
    "y_extent = 1000.0\n"
    "z_extent = 1000.0\n"
    "nx = 3\n"
    "ny = 3\n"
    "nz = 3\n"
    "\n"
    "[ROCK]\n"
    "density = 2650.0\n"
    "youngs_modulus = 80.0e9\n"
    "poissons_ratio = 0.25\n"
    "\n"
    "[MATERIAL]\n"
    "heterogeneous = false\n"
    "gravity = 9.81\n"
    "\n"
    "[GEOMECHANICS]\n"
    "enabled = true\n"
    "\n"
    "[BOUNDARY_CONDITIONS]\n"
    "bottom = fixed\n"
    "sides = roller\n"
    "top = free\n"
  );
  fclose(fp);

  Simulator sim(PETSC_COMM_WORLD);
  PetscErrorCode ierr;

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

  remove(config_path);
}
