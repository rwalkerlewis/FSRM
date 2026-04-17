/**
 * @file test_explosion_damage_zone.cpp
 * @brief Physics validation for FEM-coupled explosion damage-zone material degradation
 */

#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <petscdmplex.h>
#include <petscsys.h>
#include <petscvec.h>

#include "core/Simulator.hpp"
#include "core/FSRM.hpp"
#include "numerics/PetscFEElasticityAux.hpp"

using namespace FSRM;

namespace
{

struct AuxProbe
{
  double radius = 0.0;
  double lambda = 0.0;
};

AuxProbe findClosestCellLambda(Simulator& sim, double tx, double ty, double tz)
{
  AuxProbe best;
  best.radius = std::numeric_limits<double>::max();

  DM dm = sim.getDM();
  DM aux_dm = sim.getAuxDM();
  Vec aux_vec = sim.getAuxVector();

  PetscSection aux_section = nullptr;
  PetscErrorCode ierr = DMGetLocalSection(aux_dm, &aux_section);
  EXPECT_EQ(ierr, 0);

  const PetscScalar* aux_array = nullptr;
  ierr = VecGetArrayRead(aux_vec, &aux_array);
  EXPECT_EQ(ierr, 0);

  PetscInt cStart = 0, cEnd = 0;
  ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd);
  EXPECT_EQ(ierr, 0);

  for (PetscInt c = cStart; c < cEnd; ++c)
  {
    PetscReal volume = 0.0, centroid[3] = {0.0, 0.0, 0.0}, normal[3] = {0.0, 0.0, 0.0};
    ierr = DMPlexComputeCellGeometryFVM(dm, c, &volume, centroid, normal);
    EXPECT_EQ(ierr, 0);

    const double dx = centroid[0] - tx;
    const double dy = centroid[1] - ty;
    const double dz = centroid[2] - tz;
    const double radius = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (radius >= best.radius) continue;

    PetscInt offset = 0;
    ierr = PetscSectionGetOffset(aux_section, c, &offset);
    EXPECT_EQ(ierr, 0);
    best.radius = radius;
    best.lambda = PetscRealPart(aux_array[offset + AUX_LAMBDA]);
  }

  ierr = VecRestoreArrayRead(aux_vec, &aux_array);
  EXPECT_EQ(ierr, 0);
  return best;
}

PetscErrorCode runDamagePipeline(Simulator& sim, const char* config_path)
{
  PetscErrorCode ierr = sim.initializeFromConfigFile(config_path);
  if (ierr) return ierr;
  ierr = sim.setupDM();
  if (ierr) return ierr;
  ierr = sim.labelBoundaries();
  if (ierr) return ierr;
  ierr = sim.setupFields();
  if (ierr) return ierr;
  ierr = sim.setupPhysics();
  if (ierr) return ierr;
  ierr = sim.setupTimeStepper();
  if (ierr) return ierr;
  ierr = sim.setupSolvers();
  if (ierr) return ierr;
  return sim.setInitialConditions();
}

} // namespace

class ExplosionDamageZoneTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank_);
    if (rank_ == 0)
    {
      writeConfig(damaged_config_path_, true);
      writeConfig(undamaged_config_path_, false);
    }
    MPI_Barrier(PETSC_COMM_WORLD);
  }

  void TearDown() override
  {
    MPI_Barrier(PETSC_COMM_WORLD);
    if (rank_ == 0)
    {
      std::remove(damaged_config_path_);
      std::remove(undamaged_config_path_);
    }
  }

  void writeConfig(const char* path, bool apply_damage_zone)
  {
    std::ofstream cfg(path);
    cfg << "[SIMULATION]\n";
    cfg << "name = damage_zone_test\n";
    cfg << "start_time = 0.0\n";
    cfg << "end_time = 0.05\n";
    cfg << "dt_initial = 0.0025\n";
    cfg << "dt_min = 0.001\n";
    cfg << "dt_max = 0.005\n";
    cfg << "max_timesteps = 20\n";
    cfg << "output_frequency = 100\n";
    cfg << "fluid_model = NONE\n";
    cfg << "solid_model = ELASTIC\n";
    cfg << "enable_geomechanics = true\n";
    cfg << "enable_elastodynamics = true\n";
    cfg << "enable_faults = false\n";
    cfg << "rtol = 1.0e-6\n";
    cfg << "atol = 1.0e-8\n";
    cfg << "\n[GRID]\n";
    cfg << "nx = 8\n";
    cfg << "ny = 8\n";
    cfg << "nz = 8\n";
    cfg << "Lx = 2000.0\n";
    cfg << "Ly = 2000.0\n";
    cfg << "Lz = 2000.0\n";
    cfg << "\n[ROCK]\n";
    cfg << "density = 2650.0\n";
    cfg << "youngs_modulus = 60.0e9\n";
    cfg << "poisson_ratio = 0.25\n";
    cfg << "\n[MATERIAL]\n";
    cfg << "heterogeneous = true\n";
    cfg << "\n[EXPLOSION_SOURCE]\n";
    cfg << "type = UNDERGROUND_NUCLEAR\n";
    cfg << "yield_kt = 10.0\n";
    cfg << "depth_of_burial = 1000.0\n";
    cfg << "location_x = 1000.0\n";
    cfg << "location_y = 1000.0\n";
    cfg << "location_z = 1000.0\n";
    cfg << "onset_time = 0.0\n";
    cfg << "rise_time = 0.01\n";
    cfg << "cavity_overpressure = 1.0e10\n";
    cfg << "apply_damage_zone = " << (apply_damage_zone ? "true" : "false") << "\n";
    cfg << "\n[BOUNDARY_CONDITIONS]\n";
    cfg << "bottom = free\n";
    cfg << "sides = free\n";
    cfg << "top = free\n";
    cfg << "\n[ABSORBING_BC]\n";
    cfg << "enabled = true\n";
    cfg << "x_min = true\n";
    cfg << "x_max = true\n";
    cfg << "y_min = true\n";
    cfg << "y_max = true\n";
    cfg << "z_min = true\n";
    cfg << "z_max = false\n";
  }

  int rank_ = 0;
  static constexpr const char* damaged_config_path_ = "test_explosion_damage_zone_on.ini";
  static constexpr const char* undamaged_config_path_ = "test_explosion_damage_zone_off.ini";
};

TEST_F(ExplosionDamageZoneTest, AuxFieldsAreDegradedNearSourceAndResponseChanges)
{
  Simulator damaged_sim(PETSC_COMM_WORLD);
  PetscErrorCode ierr = runDamagePipeline(damaged_sim, damaged_config_path_);
  ASSERT_EQ(ierr, 0);

  const AuxProbe near_damaged = findClosestCellLambda(damaged_sim, 1000.0, 1000.0, 1000.0);
  const AuxProbe far_damaged = findClosestCellLambda(damaged_sim, 0.0, 0.0, 0.0);

  EXPECT_LT(near_damaged.lambda, far_damaged.lambda)
      << "The aux lambda near the explosion should be reduced by the damage-zone model";
    EXPECT_LT(near_damaged.lambda / far_damaged.lambda, 0.95)
      << "The nearest resolved cell should show measurable stiffness degradation";

  ierr = damaged_sim.run();
  ASSERT_EQ(ierr, 0);
  PetscReal damaged_norm = 0.0;
  ierr = VecNorm(damaged_sim.getSolution(), NORM_2, &damaged_norm);
  ASSERT_EQ(ierr, 0);

  Simulator undamaged_sim(PETSC_COMM_WORLD);
  ierr = runDamagePipeline(undamaged_sim, undamaged_config_path_);
  ASSERT_EQ(ierr, 0);

  const AuxProbe near_undamaged = findClosestCellLambda(undamaged_sim, 1000.0, 1000.0, 1000.0);
  const AuxProbe far_undamaged = findClosestCellLambda(undamaged_sim, 0.0, 0.0, 0.0);

  EXPECT_NEAR(near_undamaged.lambda, far_undamaged.lambda, far_undamaged.lambda * 1.0e-12)
      << "Without damage enabled, aux lambda should remain homogeneous";

  ierr = undamaged_sim.run();
  ASSERT_EQ(ierr, 0);
  PetscReal undamaged_norm = 0.0;
  ierr = VecNorm(undamaged_sim.getSolution(), NORM_2, &undamaged_norm);
  ASSERT_EQ(ierr, 0);

  EXPECT_GT(std::abs(damaged_norm - undamaged_norm), 1.0e-12)
      << "Damage-zone degradation should change the dynamic response norm";
}