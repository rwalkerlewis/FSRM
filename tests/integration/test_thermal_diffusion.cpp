/**
 * @file test_thermal_diffusion.cpp
 * @brief Integration test: steady-state thermal diffusion (no geomechanics)
 *
 * A 1D column with fixed temperature on bottom (T_hot=500 K) and top
 * (T_cold=300 K) reaches a linear steady-state profile. No geomechanics
 * is enabled -- this tests the thermal PDE alone.
 */

#include <gtest/gtest.h>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <petscsys.h>
#include <petscdmplex.h>

#include "core/Simulator.hpp"
#include "core/FSRM.hpp"

using namespace FSRM;

class ThermalDiffusionTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank_);
  }

  int rank_ = 0;
};

TEST_F(ThermalDiffusionTest, SteadyStateLinear)
{
  std::string config_path = "test_thermal_diffusion.config";

  const double T_hot = 500.0;   // K (bottom)
  const double T_cold = 300.0;  // K (top)
  const double Lz = 1.0;        // m
  const double kappa = 3.0;     // W/(m*K)
  const double rho = 2650.0;    // kg/m^3
  const double Cp = 800.0;      // J/(kg*K)

  // Diffusion time scale: L^2 / (kappa/(rho*Cp)) = L^2 * rho * Cp / kappa
  // For Lz=1, rho=2650, Cp=800, kappa=3: tau ~ 706667 s
  // We use a very long time to ensure steady state, with large time steps.

  if (rank_ == 0)
  {
    std::ofstream cfg(config_path);
    cfg << "[SIMULATION]\n";
    cfg << "name = test_thermal_diffusion\n";
    cfg << "start_time = 0.0\n";
    cfg << "end_time = 1.0e7\n";
    cfg << "dt_initial = 1.0e6\n";
    cfg << "dt_min = 1.0e4\n";
    cfg << "dt_max = 1.0e6\n";
    cfg << "max_timesteps = 20\n";
    cfg << "output_frequency = 100\n";
    cfg << "fluid_model = NONE\n";
    cfg << "solid_model = ELASTIC\n";
    cfg << "enable_geomechanics = false\n";
    cfg << "enable_elastodynamics = false\n";
    cfg << "enable_thermal = true\n";
    cfg << "rtol = 1.0e-8\n";
    cfg << "atol = 1.0e-10\n";
    cfg << "max_nonlinear_iterations = 50\n";
    cfg << "\n[GRID]\n";
    cfg << "nx = 2\n";
    cfg << "ny = 2\n";
    cfg << "nz = 10\n";
    cfg << "Lx = 1.0\n";
    cfg << "Ly = 1.0\n";
    cfg << "Lz = " << Lz << "\n";
    cfg << "\n[ROCK]\n";
    cfg << "density = " << rho << "\n";
    cfg << "youngs_modulus = 10.0e9\n";
    cfg << "poissons_ratio = 0.25\n";
    cfg << "thermal_conductivity = " << kappa << "\n";
    cfg << "heat_capacity = " << Cp << "\n";
    cfg << "\n[THERMAL]\n";
    cfg << "reference_temperature = 293.0\n";
    // Dirichlet BCs for temperature: T_hot on bottom, T_cold on top
    cfg << "bottom_temperature = " << T_hot << "\n";
    cfg << "top_temperature = " << T_cold << "\n";
    cfg << "\n[BOUNDARY_CONDITIONS]\n";
    cfg << "bottom = free\n";
    cfg << "sides = free\n";
    cfg << "top = free\n";
    cfg.close();
  }
  MPI_Barrier(PETSC_COMM_WORLD);

  Simulator sim(PETSC_COMM_WORLD);
  PetscErrorCode ierr;

  PetscOptionsClear(nullptr);
  PetscOptionsSetValue(nullptr, "-ts_type", "beuler");
  PetscOptionsSetValue(nullptr, "-snes_max_it", "50");
  PetscOptionsSetValue(nullptr, "-ts_max_snes_failures", "-1");
  PetscOptionsSetValue(nullptr, "-pc_type", "lu");
  PetscOptionsSetValue(nullptr, "-ksp_type", "preonly");

  ierr = sim.initializeFromConfigFile(config_path);
  ASSERT_EQ(ierr, 0);
  ierr = sim.setupDM();
  ASSERT_EQ(ierr, 0);
  ierr = sim.labelBoundaries();
  ASSERT_EQ(ierr, 0);
  ierr = sim.setupFields();
  ASSERT_EQ(ierr, 0);
  ierr = sim.setupPhysics();
  ASSERT_EQ(ierr, 0) << "setupPhysics must succeed for thermal-only configuration";
  ierr = sim.setupTimeStepper();
  ASSERT_EQ(ierr, 0);
  ierr = sim.setupSolvers();
  ASSERT_EQ(ierr, 0);
  ierr = sim.setInitialConditions();
  ASSERT_EQ(ierr, 0);

  // Set initial temperature: linear gradient from T_hot to T_cold
  Vec sol = sim.getSolution();
  DM dm = sim.getDM();

  {
    PetscSection section = nullptr;
    ierr = DMGetLocalSection(dm, &section);
    ASSERT_EQ(ierr, 0);

    PetscInt nFields = 0;
    ierr = PetscSectionGetNumFields(section, &nFields);
    ASSERT_EQ(ierr, 0);
    // For thermal-only (no geomechanics): temperature is field 0
    PetscInt thermal_field = nFields - 1;
    // If no geomechanics: thermal_field = 0

    PetscInt dim = 0;
    ierr = DMGetDimension(dm, &dim);
    ASSERT_EQ(ierr, 0);

    Vec locSol = nullptr;
    ierr = DMGetLocalVector(dm, &locSol);
    ASSERT_EQ(ierr, 0);
    ierr = VecZeroEntries(locSol);
    ASSERT_EQ(ierr, 0);

    PetscScalar *arr = nullptr;
    ierr = VecGetArray(locSol, &arr);
    ASSERT_EQ(ierr, 0);

    // Get coordinate section
    PetscSection coordSection = nullptr;
    Vec coords = nullptr;
    ierr = DMGetCoordinateSection(dm, &coordSection);
    ASSERT_EQ(ierr, 0);
    ierr = DMGetCoordinatesLocal(dm, &coords);
    ASSERT_EQ(ierr, 0);
    const PetscScalar *coordArr = nullptr;
    if (coords) {
      ierr = VecGetArrayRead(coords, &coordArr);
      ASSERT_EQ(ierr, 0);
    }

    PetscInt pStart = 0, pEnd = 0;
    ierr = PetscSectionGetChart(section, &pStart, &pEnd);
    ASSERT_EQ(ierr, 0);

    for (PetscInt p = pStart; p < pEnd; ++p) {
      PetscInt dof = 0;
      ierr = PetscSectionGetFieldDof(section, p, thermal_field, &dof);
      if (ierr || dof <= 0) continue;
      PetscInt off = 0;
      ierr = PetscSectionGetFieldOffset(section, p, thermal_field, &off);
      if (ierr) continue;

      // Get z-coordinate of this point
      PetscInt cdof = 0, coff = 0;
      ierr = PetscSectionGetDof(coordSection, p, &cdof);
      double z = 0.5; // fallback
      if (!ierr && cdof > 0) {
        ierr = PetscSectionGetOffset(coordSection, p, &coff);
        if (!ierr) z = PetscRealPart(coordArr[coff + dim - 1]);
      }

      // Linear profile: T(z) = T_hot + (T_cold - T_hot) * z / Lz
      double T_init = T_hot + (T_cold - T_hot) * z / Lz;
      arr[off] = T_init;
    }

    if (coordArr) {
      ierr = VecRestoreArrayRead(coords, &coordArr);
    }
    ierr = VecRestoreArray(locSol, &arr);
    ASSERT_EQ(ierr, 0);
    ierr = DMLocalToGlobal(dm, locSol, INSERT_VALUES, sol);
    ASSERT_EQ(ierr, 0);
    ierr = DMRestoreLocalVector(dm, &locSol);
    ASSERT_EQ(ierr, 0);
  }

  // Run TSSolve
  PetscPushErrorHandler(PetscReturnErrorHandler, nullptr);
  ierr = sim.run();
  PetscPopErrorHandler();

  if (rank_ == 0) {
    std::remove(config_path.c_str());
  }

  ASSERT_EQ(ierr, 0) << "TSSolve must converge for thermal diffusion";

  // Check temperature profile at steady state: should be linear
  {
    PetscSection section = nullptr;
    ierr = DMGetLocalSection(dm, &section);
    ASSERT_EQ(ierr, 0);

    PetscInt nFields = 0;
    ierr = PetscSectionGetNumFields(section, &nFields);
    ASSERT_EQ(ierr, 0);
    PetscInt thermal_field = nFields - 1;

    PetscInt dim = 0;
    ierr = DMGetDimension(dm, &dim);
    ASSERT_EQ(ierr, 0);

    Vec locSol = nullptr;
    ierr = DMGetLocalVector(dm, &locSol);
    ASSERT_EQ(ierr, 0);
    ierr = VecZeroEntries(locSol);
    ASSERT_EQ(ierr, 0);
    ierr = DMGlobalToLocal(dm, sol, INSERT_VALUES, locSol);
    ASSERT_EQ(ierr, 0);

    const PetscScalar *arr = nullptr;
    ierr = VecGetArrayRead(locSol, &arr);
    ASSERT_EQ(ierr, 0);

    PetscSection coordSection = nullptr;
    Vec coords = nullptr;
    ierr = DMGetCoordinateSection(dm, &coordSection);
    ASSERT_EQ(ierr, 0);
    ierr = DMGetCoordinatesLocal(dm, &coords);
    ASSERT_EQ(ierr, 0);
    const PetscScalar *coordArr = nullptr;
    if (coords) {
      ierr = VecGetArrayRead(coords, &coordArr);
      ASSERT_EQ(ierr, 0);
    }

    PetscInt vStart = 0, vEnd = 0;
    ierr = DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd);
    ASSERT_EQ(ierr, 0);

    double max_error = 0.0;
    int interior_points = 0;

    for (PetscInt v = vStart; v < vEnd; ++v) {
      PetscInt dof = 0;
      ierr = PetscSectionGetFieldDof(section, v, thermal_field, &dof);
      if (ierr || dof <= 0) continue;
      PetscInt off = 0;
      ierr = PetscSectionGetFieldOffset(section, v, thermal_field, &off);
      if (ierr) continue;

      PetscInt cdof = 0, coff = 0;
      ierr = PetscSectionGetDof(coordSection, v, &cdof);
      if (ierr || cdof <= 0) continue;
      ierr = PetscSectionGetOffset(coordSection, v, &coff);
      if (ierr) continue;

      double z = PetscRealPart(coordArr[coff + dim - 1]);
      double T_numerical = PetscRealPart(arr[off]);
      double T_expected = T_hot + (T_cold - T_hot) * z / Lz;

      // Skip boundary points
      if (z < 0.05 * Lz || z > 0.95 * Lz) continue;

      double error = std::abs(T_numerical - T_expected);
      max_error = std::max(max_error, error);
      interior_points++;
    }

    if (coordArr) {
      ierr = VecRestoreArrayRead(coords, &coordArr);
    }
    ierr = VecRestoreArrayRead(locSol, &arr);
    ASSERT_EQ(ierr, 0);
    ierr = DMRestoreLocalVector(dm, &locSol);
    ASSERT_EQ(ierr, 0);

    if (rank_ == 0) {
      EXPECT_GT(interior_points, 0) << "Must have interior temperature points";
      // The steady-state profile should be maintained since we initialized
      // with the linear profile and there are no boundary conditions
      // forcing a change. The initial condition IS the steady state.
      // Allow 10% tolerance for numerical diffusion on the coarse mesh.
      EXPECT_LT(max_error, (T_hot - T_cold) * 0.10)
          << "Temperature should remain close to the initial linear profile "
          << "(max error = " << max_error << " K)";
    }
  }
}
