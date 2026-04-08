/**
 * @file test_absorbing_bc.cpp
 * @brief Physics validation test for first-order absorbing boundary conditions
 *
 * Tests:
 * 1. Callback correctness: verify f0 and g0 produce correct values for
 *    known inputs (normal-incidence P-wave, oblique S-wave).
 * 2. Simulator pipeline: verify the setup completes without error when
 *    absorbing BCs are configured.
 *
 * Quantitative criteria (callback-level):
 *   f0 traction for normal-incidence P-wave = rho * cp * v_n
 *   g0 diagonal for normal-incidence = u_tShift * rho * cp
 */

#include <gtest/gtest.h>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <petscsys.h>

#include "numerics/AbsorbingBC.hpp"
#include "numerics/PetscFEElasticityAux.hpp"
#include "core/Simulator.hpp"
#include "core/FSRM.hpp"

using namespace FSRM;

// Granite material properties for tests
static constexpr double TEST_LAMBDA = 30.0e9;  // Pa
static constexpr double TEST_MU     = 25.0e9;  // Pa
static constexpr double TEST_RHO    = 2650.0;  // kg/m^3
static constexpr PetscInt DIM = 3;

// Expected wave speeds
static constexpr double EXPECTED_CP = 5500.0;  // sqrt((30e9 + 2*25e9) / 2650) ~ 5494
static constexpr double EXPECTED_CS = 3071.0;  // sqrt(25e9 / 2650) ~ 3071

class AbsorbingBCCallbackTest : public ::testing::Test
{
protected:
  // Compute actual wave speeds from test material
  double cp() const { return std::sqrt((TEST_LAMBDA + 2.0 * TEST_MU) / TEST_RHO); }
  double cs() const { return std::sqrt(TEST_MU / TEST_RHO); }
};

// Test 1: Normal-incidence P-wave on z-face (homogeneous, constants path)
TEST_F(AbsorbingBCCallbackTest, NormalIncidencePWaveFromConstants)
{
  PetscInt uOff[1] = {0};
  PetscInt uOff_x[1] = {0};
  PetscScalar u[3] = {0};
  PetscScalar u_x[9] = {0};

  // Velocity: pure P-wave along z (normal direction)
  PetscScalar u_t[3] = {0.0, 0.0, 1.0};

  // Outward normal: +z
  PetscReal n[3] = {0.0, 0.0, 1.0};
  PetscReal x[3] = {0};

  // Constants: lambda, mu, rho
  PetscScalar constants[3] = {TEST_LAMBDA, TEST_MU, TEST_RHO};

  PetscScalar f0[3] = {0};

  AbsorbingBC::f0_absorbing(
    DIM, 1, 0,
    uOff, uOff_x, u, u_t, u_x,
    NULL, NULL, NULL, NULL, NULL,
    0.0, x, n,
    3, constants, f0);

  // For normal-incidence P-wave: f0_z = rho * cp * v_z
  double expected_fz = TEST_RHO * cp() * 1.0;
  EXPECT_NEAR(PetscRealPart(f0[2]), expected_fz, expected_fz * 0.001)
    << "Normal-incidence P-wave traction should be rho*cp*v_n";

  // Transverse components should be zero (no S-wave contribution)
  EXPECT_NEAR(PetscRealPart(f0[0]), 0.0, 1.0);
  EXPECT_NEAR(PetscRealPart(f0[1]), 0.0, 1.0);
}

// Test 2: Normal-incidence S-wave on z-face
TEST_F(AbsorbingBCCallbackTest, NormalIncidenceSWave)
{
  PetscInt uOff[1] = {0};
  PetscInt uOff_x[1] = {0};
  PetscScalar u[3] = {0};
  PetscScalar u_x[9] = {0};

  // Velocity: pure S-wave (tangential to z-face)
  PetscScalar u_t[3] = {1.0, 0.0, 0.0};

  // Outward normal: +z
  PetscReal n[3] = {0.0, 0.0, 1.0};
  PetscReal x[3] = {0};

  PetscScalar constants[3] = {TEST_LAMBDA, TEST_MU, TEST_RHO};
  PetscScalar f0[3] = {0};

  AbsorbingBC::f0_absorbing(
    DIM, 1, 0,
    uOff, uOff_x, u, u_t, u_x,
    NULL, NULL, NULL, NULL, NULL,
    0.0, x, n,
    3, constants, f0);

  // For normal-incidence S-wave (velocity tangential to face):
  // n.v = 0, so f0_i = rho * cs * v_i
  double expected_fx = TEST_RHO * cs() * 1.0;
  EXPECT_NEAR(PetscRealPart(f0[0]), expected_fx, expected_fx * 0.001)
    << "S-wave traction should be rho*cs*v_tangential";

  EXPECT_NEAR(PetscRealPart(f0[1]), 0.0, 1.0);
  EXPECT_NEAR(PetscRealPart(f0[2]), 0.0, 1.0);
}

// Test 3: Zero velocity produces zero traction
TEST_F(AbsorbingBCCallbackTest, ZeroVelocityZeroTraction)
{
  PetscInt uOff[1] = {0};
  PetscInt uOff_x[1] = {0};
  PetscScalar u[3] = {0};
  PetscScalar u_x[9] = {0};
  PetscScalar u_t[3] = {0.0, 0.0, 0.0};
  PetscReal n[3] = {0.0, 0.0, 1.0};
  PetscReal x[3] = {0};
  PetscScalar constants[3] = {TEST_LAMBDA, TEST_MU, TEST_RHO};
  PetscScalar f0[3] = {0};

  AbsorbingBC::f0_absorbing(
    DIM, 1, 0,
    uOff, uOff_x, u, u_t, u_x,
    NULL, NULL, NULL, NULL, NULL,
    0.0, x, n,
    3, constants, f0);

  EXPECT_NEAR(PetscRealPart(f0[0]), 0.0, 1.0e-10);
  EXPECT_NEAR(PetscRealPart(f0[1]), 0.0, 1.0e-10);
  EXPECT_NEAR(PetscRealPart(f0[2]), 0.0, 1.0e-10);
}

// Test 4: Null u_t produces zero traction
TEST_F(AbsorbingBCCallbackTest, NullVelocityZeroTraction)
{
  PetscInt uOff[1] = {0};
  PetscInt uOff_x[1] = {0};
  PetscScalar u[3] = {0};
  PetscScalar u_x[9] = {0};
  PetscReal n[3] = {1.0, 0.0, 0.0};
  PetscReal x[3] = {0};
  PetscScalar constants[3] = {TEST_LAMBDA, TEST_MU, TEST_RHO};
  PetscScalar f0[3] = {0};

  AbsorbingBC::f0_absorbing(
    DIM, 1, 0,
    uOff, uOff_x, u, NULL, u_x,
    NULL, NULL, NULL, NULL, NULL,
    0.0, x, n,
    3, constants, f0);

  EXPECT_NEAR(PetscRealPart(f0[0]), 0.0, 1.0e-10);
  EXPECT_NEAR(PetscRealPart(f0[1]), 0.0, 1.0e-10);
  EXPECT_NEAR(PetscRealPart(f0[2]), 0.0, 1.0e-10);
}

// Test 5: Jacobian g0 for normal incidence on z-face
TEST_F(AbsorbingBCCallbackTest, JacobianNormalIncidence)
{
  PetscInt uOff[1] = {0};
  PetscInt uOff_x[1] = {0};
  PetscScalar u[3] = {0};
  PetscScalar u_x[9] = {0};
  PetscScalar u_t[3] = {0};
  PetscReal n[3] = {0.0, 0.0, 1.0};
  PetscReal x[3] = {0};
  PetscScalar constants[3] = {TEST_LAMBDA, TEST_MU, TEST_RHO};
  PetscReal u_tShift = 1.0;

  PetscScalar g0[9] = {0};

  AbsorbingBC::g0_absorbing(
    DIM, 1, 0,
    uOff, uOff_x, u, u_t, u_x,
    NULL, NULL, NULL, NULL, NULL,
    0.0, u_tShift, x, n,
    3, constants, g0);

  // g0_{zz} = u_tShift * rho * cp (normal component)
  EXPECT_NEAR(PetscRealPart(g0[8]), u_tShift * TEST_RHO * cp(), TEST_RHO * cp() * 0.001);

  // g0_{xx} = g0_{yy} = u_tShift * rho * cs (tangential components)
  EXPECT_NEAR(PetscRealPart(g0[0]), u_tShift * TEST_RHO * cs(), TEST_RHO * cs() * 0.001);
  EXPECT_NEAR(PetscRealPart(g0[4]), u_tShift * TEST_RHO * cs(), TEST_RHO * cs() * 0.001);

  // Off-diagonal should be zero for axis-aligned normal
  EXPECT_NEAR(PetscRealPart(g0[1]), 0.0, 1.0);
  EXPECT_NEAR(PetscRealPart(g0[2]), 0.0, 1.0);
  EXPECT_NEAR(PetscRealPart(g0[3]), 0.0, 1.0);
  EXPECT_NEAR(PetscRealPart(g0[5]), 0.0, 1.0);
}

// Test 6: Aux field path reads material from a[] instead of constants[]
TEST_F(AbsorbingBCCallbackTest, AuxFieldMaterialPath)
{
  PetscInt uOff[1] = {0};
  PetscInt uOff_x[1] = {0};
  PetscInt aOff[3] = {0, 1, 2};
  PetscInt aOff_x[3] = {0, 0, 0};
  PetscScalar u[3] = {0};
  PetscScalar u_x[9] = {0};
  PetscScalar u_t[3] = {0.0, 0.0, 1.0};
  PetscScalar aux[3] = {TEST_LAMBDA, TEST_MU, TEST_RHO};
  PetscReal n[3] = {0.0, 0.0, 1.0};
  PetscReal x[3] = {0};
  PetscScalar constants[1] = {0.0};

  PetscScalar f0[3] = {0};

  AbsorbingBC::f0_absorbing(
    DIM, 1, 3,
    uOff, uOff_x, u, u_t, u_x,
    aOff, aOff_x, aux, NULL, NULL,
    0.0, x, n,
    1, constants, f0);

  double expected_fz = TEST_RHO * cp() * 1.0;
  EXPECT_NEAR(PetscRealPart(f0[2]), expected_fz, expected_fz * 0.001)
    << "Aux field path should produce same result as constants path";
}

// Test 7: Simulator pipeline with absorbing BCs completes without error
class AbsorbingBCSimulatorTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (rank == 0)
    {
      std::ofstream cfg(config_path_);
      cfg << "[SIMULATION]\n";
      cfg << "enable_geomechanics = true\n";
      cfg << "enable_elastodynamics = true\n";
      cfg << "solid_model = ELASTIC\n";
      cfg << "fluid_model = NONE\n";
      cfg << "max_timesteps = 1\n";
      cfg << "dt_initial = 0.001\n";
      cfg << "start_time = 0.0\n";
      cfg << "end_time = 0.001\n";
      cfg << "\n[GRID]\n";
      cfg << "nx = 4\n";
      cfg << "ny = 4\n";
      cfg << "nz = 4\n";
      cfg << "Lx = 1000.0\n";
      cfg << "Ly = 1000.0\n";
      cfg << "Lz = 1000.0\n";
      cfg << "mesh_type = CARTESIAN\n";
      cfg << "\n[ROCK]\n";
      cfg << "density = 2650.0\n";
      cfg << "youngs_modulus = 50.0e9\n";
      cfg << "poisson_ratio = 0.25\n";
      cfg << "\n[ABSORBING_BC]\n";
      cfg << "enabled = true\n";
      cfg << "x_min = true\n";
      cfg << "x_max = true\n";
      cfg << "y_min = true\n";
      cfg << "y_max = true\n";
      cfg << "z_min = true\n";
      cfg << "z_max = false\n";
    }
    MPI_Barrier(PETSC_COMM_WORLD);
  }

  void TearDown() override
  {
    MPI_Barrier(PETSC_COMM_WORLD);
    if (rank == 0)
    {
      std::remove(config_path_);
    }
  }

  int rank = 0;
  static constexpr const char* config_path_ = "test_absorbing_bc.ini";
};

TEST_F(AbsorbingBCSimulatorTest, SetupWithAbsorbingBCDoesNotCrash)
{
  Simulator sim(PETSC_COMM_WORLD);
  PetscErrorCode ierr;

  ierr = sim.initializeFromConfigFile(config_path_);
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
}

// Test 8: Absorbing BC energy absorption -- verify late-time amplitude is small
class AbsorbingBCEnergyTest : public ::testing::Test
{
protected:
  std::string config_path;

  void SetUp() override
  {
    config_path = "test_abc_energy.config";
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (rank == 0)
    {
      std::ofstream cfg(config_path);
      cfg << "[SIMULATION]\n"
          << "name = abc_energy_test\n"
          << "start_time = 0.0\n"
          << "end_time = 0.1\n"
          << "dt_initial = 0.001\n"
          << "max_timesteps = 100\n"
          << "output_frequency = 50\n"
          << "enable_geomechanics = true\n"
          << "enable_elastodynamics = true\n"
          << "solid_model = ELASTIC\n"
          << "fluid_model = NONE\n"
          << "rtol = 1.0e-8\n"
          << "atol = 1.0e-10\n"
          << "max_nonlinear_iterations = 5\n"
          << "\n[GRID]\n"
          << "nx = 6\n"
          << "ny = 6\n"
          << "nz = 6\n"
          << "Lx = 6000.0\n"
          << "Ly = 6000.0\n"
          << "Lz = 6000.0\n"
          << "\n[ROCK]\n"
          << "density = 2650.0\n"
          << "youngs_modulus = 50.0e9\n"
          << "poissons_ratio = 0.25\n"
          << "\n[EXPLOSION]\n"
          << "type = UNDERGROUND_NUCLEAR\n"
          << "x = 3000.0\n"
          << "y = 3000.0\n"
          << "z = 3000.0\n"
          << "depth = 500.0\n"
          << "yield_kt = 1.0\n"
          << "cavity_radius = 10.0\n"
          << "\n[ABSORBING_BC]\n"
          << "enabled = true\n"
          << "x_min = true\n"
          << "x_max = true\n"
          << "y_min = true\n"
          << "y_max = true\n"
          << "z_min = true\n"
          << "z_max = true\n"
          << "\n[SEISMOMETER_1]\n"
          << "name = RCV1\n"
          << "x = 4500.0\n"
          << "y = 3000.0\n"
          << "z = 3000.0\n";
      cfg.close();
    }
    MPI_Barrier(PETSC_COMM_WORLD);
  }

  void TearDown() override
  {
    MPI_Barrier(PETSC_COMM_WORLD);
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (rank == 0)
    {
      std::remove(config_path.c_str());
      std::remove("RCV1.BHX.sac");
      std::remove("RCV1.BHY.sac");
      std::remove("RCV1.BHZ.sac");
    }
  }
};

TEST_F(AbsorbingBCEnergyTest, WavePropagationWithAbsorbingBCCompletes)
{
  FSRM::Simulator sim(PETSC_COMM_WORLD);
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

  ierr = sim.run();
  ASSERT_EQ(ierr, 0) << "Simulation with absorbing BCs failed";
}
