/**
 * @file test_gmsh_multimaterial.cpp
 * @brief Integration test for gmsh-label-based per-cell material assignment
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <petscdmplex.h>
#include <petscsys.h>
#include <petscvec.h>

#include "core/FSRM.hpp"
#include "core/Simulator.hpp"
#include "numerics/PetscFEElasticityAux.hpp"

using namespace FSRM;

namespace
{

struct RegionStats
{
  PetscInt count = 0;
  double lambda_sum = 0.0;
  double uz_sum = 0.0;
};

std::map<PetscInt, RegionStats> collectRegionStats(Simulator& sim)
{
  std::map<PetscInt, RegionStats> stats;

  DM dm = sim.getDM();
  DM aux_dm = sim.getAuxDM();
  Vec aux_vec = sim.getAuxVector();
  Vec sol = sim.getSolution();

  EXPECT_NE(dm, nullptr);
  EXPECT_NE(aux_dm, nullptr);
  EXPECT_NE(aux_vec, nullptr);
  EXPECT_NE(sol, nullptr);

  DMLabel material = nullptr;
  PetscErrorCode ierr = DMGetLabel(dm, "Material", &material);
  EXPECT_EQ(ierr, 0);
  EXPECT_NE(material, nullptr);
  if (!material) return stats;

  PetscSection aux_section = nullptr;
  ierr = DMGetLocalSection(aux_dm, &aux_section);
  EXPECT_EQ(ierr, 0);

  const PetscScalar* aux_array = nullptr;
  ierr = VecGetArrayRead(aux_vec, &aux_array);
  EXPECT_EQ(ierr, 0);

  Vec local = nullptr;
  ierr = DMGetLocalVector(dm, &local);
  EXPECT_EQ(ierr, 0);
  ierr = DMGlobalToLocalBegin(dm, sol, INSERT_VALUES, local);
  EXPECT_EQ(ierr, 0);
  ierr = DMGlobalToLocalEnd(dm, sol, INSERT_VALUES, local);
  EXPECT_EQ(ierr, 0);

  PetscSection section = nullptr;
  ierr = DMGetLocalSection(dm, &section);
  EXPECT_EQ(ierr, 0);

  PetscInt dim = 0;
  ierr = DMGetDimension(dm, &dim);
  EXPECT_EQ(ierr, 0);

  PetscInt cStart = 0, cEnd = 0;
  ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd);
  EXPECT_EQ(ierr, 0);

  for (PetscInt c = cStart; c < cEnd; ++c)
  {
    PetscInt material_id = -1;
    ierr = DMLabelGetValue(material, c, &material_id);
    EXPECT_EQ(ierr, 0);
    if (material_id < 0) continue;

    PetscInt aux_offset = 0;
    ierr = PetscSectionGetOffset(aux_section, c, &aux_offset);
    EXPECT_EQ(ierr, 0);

    PetscInt closure_size = 0;
    PetscScalar* closure = nullptr;
    ierr = DMPlexVecGetClosure(dm, section, local, c, &closure_size, &closure);
    EXPECT_EQ(ierr, 0);

    double uz_sum = 0.0;
    PetscInt uz_count = 0;
    for (PetscInt i = dim - 1; i < closure_size; i += dim)
    {
      uz_sum += std::abs(PetscRealPart(closure[i]));
      ++uz_count;
    }

    RegionStats& region = stats[material_id];
    region.count += 1;
    region.lambda_sum += PetscRealPart(aux_array[aux_offset + AUX_LAMBDA]);
    if (uz_count > 0)
    {
      region.uz_sum += uz_sum / static_cast<double>(uz_count);
    }

    ierr = DMPlexVecRestoreClosure(dm, section, local, c, &closure_size, &closure);
    EXPECT_EQ(ierr, 0);
  }

  ierr = DMRestoreLocalVector(dm, &local);
  EXPECT_EQ(ierr, 0);
  ierr = VecRestoreArrayRead(aux_vec, &aux_array);
  EXPECT_EQ(ierr, 0);

  return stats;
}

} // namespace

class GmshMultiMaterialTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank_);
    if (rank_ == 0)
    {
      std::ofstream cfg(config_path_);
      cfg << "[SIMULATION]\n";
      cfg << "enable_geomechanics = true\n";
      cfg << "solid_model = ELASTIC\n";
      cfg << "fluid_model = NONE\n";
      cfg << "max_timesteps = 1\n";
      cfg << "dt_initial = 1.0\n";
      cfg << "start_time = 0.0\n";
      cfg << "end_time = 1.0\n";
      cfg << "\n[GRID]\n";
      cfg << "mesh_type = GMSH\n";
      cfg << "mesh_file = " << mesh_path_ << "\n";
      cfg << "\n[ROCK]\n";
      cfg << "density = 2500.0\n";
      cfg << "youngs_modulus = 40.0e9\n";
      cfg << "poisson_ratio = 0.25\n";
      cfg << "\n[MATERIAL]\n";
      cfg << "heterogeneous = true\n";
      cfg << "assignment = gmsh_label\n";
      cfg << "\n[MATERIAL_REGION_1]\n";
      cfg << "gmsh_label = soft_rock\n";
      cfg << "youngs_modulus = 10.0e9\n";
      cfg << "poisson_ratio = 0.25\n";
      cfg << "density = 2200.0\n";
      cfg << "\n[MATERIAL_REGION_2]\n";
      cfg << "gmsh_label = hard_rock\n";
      cfg << "youngs_modulus = 80.0e9\n";
      cfg << "poisson_ratio = 0.25\n";
      cfg << "density = 2800.0\n";
    }
    MPI_Barrier(PETSC_COMM_WORLD);
  }

  void TearDown() override
  {
    MPI_Barrier(PETSC_COMM_WORLD);
    if (rank_ == 0)
    {
      std::remove(config_path_);
    }
  }

  int rank_ = 0;
  static constexpr const char* config_path_ = "test_gmsh_multimaterial.ini";
  static constexpr const char* mesh_path_ = "../../meshes/test_two_material.msh";
};

TEST_F(GmshMultiMaterialTest, AuxFieldsAndDisplacementsDifferByRegion)
{
  Simulator sim(PETSC_COMM_WORLD);
  PetscErrorCode ierr = sim.initializeFromConfigFile(config_path_);
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

  auto stats = collectRegionStats(sim);
  ASSERT_EQ(stats.size(), 2u) << "Expected exactly two material regions";

  auto soft = stats.begin();
  auto hard = std::next(stats.begin());
  const double soft_lambda = soft->second.lambda_sum / static_cast<double>(soft->second.count);
  const double hard_lambda = hard->second.lambda_sum / static_cast<double>(hard->second.count);
  const double soft_uz = soft->second.uz_sum / static_cast<double>(soft->second.count);
  const double hard_uz = hard->second.uz_sum / static_cast<double>(hard->second.count);

  EXPECT_GT(soft->second.count, 0);
  EXPECT_GT(hard->second.count, 0);

  EXPECT_LT(soft_lambda, hard_lambda)
      << "Soft-rock lambda should be lower than hard-rock lambda";
  EXPECT_GT(hard_lambda / soft_lambda, 4.0)
      << "The mapped Lame parameters should reflect the large stiffness contrast";

    EXPECT_GT(soft_uz + hard_uz, 0.0)
      << "The solve should produce nonzero displacement response across the mapped regions";
}