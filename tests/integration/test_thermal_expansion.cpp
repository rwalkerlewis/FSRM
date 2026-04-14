/**
 * @file test_thermal_expansion.cpp
 * @brief Integration test: thermoelastic coupling with uniform heating
 *
 * A 3D box is uniformly heated above the reference temperature with
 * fixed bottom and free top/sides. The expected result is free thermal
 * expansion: displacement at top ~ alpha_T * dT * Lz, with approximately
 * zero stress (no mechanical constraint except the fixed bottom).
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

class ThermalExpansionTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank_);
  }

  int rank_ = 0;
};

TEST_F(ThermalExpansionTest, UniformHeating)
{
  std::string config_path = "test_thermal_expansion.config";

  // Material and thermal parameters
  const double E = 10.0e9;           // Pa
  const double nu = 0.25;
  const double rho = 2650.0;         // kg/m^3
  const double kappa = 3.0;          // W/(m*K)
  const double Cp = 800.0;           // J/(kg*K)
  const double alpha_T = 1.0e-5;     // 1/K
  const double T_ref = 293.0;        // K
  const double dT = 100.0;           // K heating
  const double Lz = 1.0;             // m

  // Expected displacement at top: eps = alpha_T * dT = 1e-3
  //   uz_top = eps * Lz = 1e-3 m
  const double expected_uz_top = alpha_T * dT * Lz;

  if (rank_ == 0)
  {
    std::ofstream cfg(config_path);
    cfg << "[SIMULATION]\n";
    cfg << "name = test_thermal_expansion\n";
    cfg << "start_time = 0.0\n";
    cfg << "end_time = 1.0\n";
    cfg << "dt_initial = 1.0\n";
    cfg << "dt_min = 0.1\n";
    cfg << "dt_max = 1.0\n";
    cfg << "max_timesteps = 1\n";
    cfg << "output_frequency = 100\n";
    cfg << "fluid_model = NONE\n";
    cfg << "solid_model = ELASTIC\n";
    cfg << "enable_geomechanics = true\n";
    cfg << "enable_elastodynamics = false\n";
    cfg << "enable_thermal = true\n";
    cfg << "rtol = 1.0e-6\n";
    cfg << "atol = 1.0e-8\n";
    cfg << "max_nonlinear_iterations = 50\n";
    cfg << "\n[GRID]\n";
    cfg << "nx = 4\n";
    cfg << "ny = 4\n";
    cfg << "nz = 4\n";
    cfg << "Lx = 1.0\n";
    cfg << "Ly = 1.0\n";
    cfg << "Lz = " << Lz << "\n";
    cfg << "\n[ROCK]\n";
    cfg << "density = " << rho << "\n";
    cfg << "youngs_modulus = " << E << "\n";
    cfg << "poissons_ratio = " << nu << "\n";
    cfg << "thermal_conductivity = " << kappa << "\n";
    cfg << "heat_capacity = " << Cp << "\n";
    cfg << "thermal_expansion = " << alpha_T << "\n";
    cfg << "\n[THERMAL]\n";
    cfg << "reference_temperature = " << T_ref << "\n";
    cfg << "\n[BOUNDARY_CONDITIONS]\n";
    cfg << "bottom = fixed\n";
    cfg << "sides = roller\n";
    cfg << "top = free\n";
    cfg.close();
  }
  MPI_Barrier(PETSC_COMM_WORLD);

  Simulator sim(PETSC_COMM_WORLD);
  PetscErrorCode ierr;

  PetscOptionsClear(nullptr);
  PetscOptionsSetValue(nullptr, "-ts_type", "beuler");
  PetscOptionsSetValue(nullptr, "-snes_max_it", "50");
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
  ASSERT_EQ(ierr, 0);
  ierr = sim.setupTimeStepper();
  ASSERT_EQ(ierr, 0);
  ierr = sim.setupSolvers();
  ASSERT_EQ(ierr, 0);
  ierr = sim.setInitialConditions();
  ASSERT_EQ(ierr, 0);

  // Set initial temperature to T_ref + dT everywhere.
  // The solution vector contains [displacement(3), temperature(1)].
  // We need to set the temperature DOFs to T_ref + dT.
  Vec sol = sim.getSolution();
  DM dm = sim.getDM();

  {
    PetscSection section = nullptr;
    ierr = DMGetLocalSection(dm, &section);
    ASSERT_EQ(ierr, 0);

    // Get number of fields to find the thermal field
    PetscInt nFields = 0;
    ierr = PetscSectionGetNumFields(section, &nFields);
    ASSERT_EQ(ierr, 0);
    // Thermal field should be the last field
    PetscInt thermal_field = nFields - 1;

    // Set temperature in local vector, then scatter to global
    Vec locSol = nullptr;
    ierr = DMGetLocalVector(dm, &locSol);
    ASSERT_EQ(ierr, 0);
    ierr = VecZeroEntries(locSol);
    ASSERT_EQ(ierr, 0);

    PetscScalar *arr = nullptr;
    ierr = VecGetArray(locSol, &arr);
    ASSERT_EQ(ierr, 0);

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
      arr[off] = T_ref + dT;
    }

    ierr = VecRestoreArray(locSol, &arr);
    ASSERT_EQ(ierr, 0);

    // Scatter local to global
    ierr = DMLocalToGlobal(dm, locSol, INSERT_VALUES, sol);
    ASSERT_EQ(ierr, 0);
    ierr = DMRestoreLocalVector(dm, &locSol);
    ASSERT_EQ(ierr, 0);
  }

  // Run TSSolve (quasi-static: 1 step, steady state)
  PetscPushErrorHandler(PetscReturnErrorHandler, nullptr);
  ierr = sim.run();
  PetscPopErrorHandler();

  if (rank_ == 0) {
    std::remove(config_path.c_str());
  }

  ASSERT_EQ(ierr, 0) << "TSSolve must converge for uniform thermal expansion";

  // Check displacement at the top surface
  // The displacement field is field 0 for fluid_model=NONE
  {
    Vec locSol = nullptr;
    ierr = DMGetLocalVector(dm, &locSol);
    ASSERT_EQ(ierr, 0);
    ierr = VecZeroEntries(locSol);
    ASSERT_EQ(ierr, 0);
    ierr = DMGlobalToLocal(dm, sol, INSERT_VALUES, locSol);
    ASSERT_EQ(ierr, 0);

    PetscSection section = nullptr;
    ierr = DMGetLocalSection(dm, &section);
    ASSERT_EQ(ierr, 0);

    const PetscScalar *arr = nullptr;
    ierr = VecGetArrayRead(locSol, &arr);
    ASSERT_EQ(ierr, 0);

    PetscInt dim = 0;
    ierr = DMGetDimension(dm, &dim);
    ASSERT_EQ(ierr, 0);

    // Find max z-displacement at vertices near z=Lz
    double max_uz = 0.0;
    int top_vertices = 0;

    PetscInt vStart = 0, vEnd = 0;
    ierr = DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd);
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

    for (PetscInt v = vStart; v < vEnd; ++v) {
      // Get coordinates
      PetscInt cdof = 0, coff = 0;
      ierr = PetscSectionGetDof(coordSection, v, &cdof);
      if (ierr || cdof <= 0) continue;
      ierr = PetscSectionGetOffset(coordSection, v, &coff);
      if (ierr) continue;

      double z = PetscRealPart(coordArr[coff + dim - 1]);
      if (z > Lz * 0.9) {
        // Read displacement field (field 0)
        PetscInt dof = 0, off = 0;
        ierr = PetscSectionGetFieldDof(section, v, 0, &dof);
        if (ierr || dof < dim) continue;
        ierr = PetscSectionGetFieldOffset(section, v, 0, &off);
        if (ierr) continue;

        double uz = PetscRealPart(arr[off + dim - 1]);
        max_uz = std::max(max_uz, uz);
        top_vertices++;
      }
    }

    if (coordArr) {
      ierr = VecRestoreArrayRead(coords, &coordArr);
    }
    ierr = VecRestoreArrayRead(locSol, &arr);
    ASSERT_EQ(ierr, 0);
    ierr = DMRestoreLocalVector(dm, &locSol);
    ASSERT_EQ(ierr, 0);

    if (rank_ == 0) {
      EXPECT_GT(top_vertices, 0) << "Must have vertices at the top surface";
      // Allow generous tolerance since: (a) fixed bottom constrains expansion
      // pattern, (b) coarse mesh, (c) quasi-static solve at t=1 step
      EXPECT_GT(max_uz, expected_uz_top * 0.5)
          << "Top z-displacement should be at least 50% of alpha_T*dT*Lz = "
          << expected_uz_top;
      EXPECT_LT(max_uz, expected_uz_top * 2.0)
          << "Top z-displacement should not exceed 200% of expected";
      EXPECT_TRUE(std::isfinite(max_uz))
          << "Displacement must be finite";
    }
  }
}
