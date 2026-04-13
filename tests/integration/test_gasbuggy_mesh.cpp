/**
 * @file test_gasbuggy_mesh.cpp
 * @brief Integration test for the historical Gasbuggy gmsh mesh
 */

#include <gtest/gtest.h>

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

std::map<PetscInt, double> collectAverageLambdaByMaterial(Simulator& sim)
{
  std::map<PetscInt, double> sums;
  std::map<PetscInt, PetscInt> counts;

  DM dm = sim.getDM();
  DM aux_dm = sim.getAuxDM();
  Vec aux_vec = sim.getAuxVector();

  DMLabel material = nullptr;
  PetscErrorCode ierr = DMGetLabel(dm, "Material", &material);
  EXPECT_EQ(ierr, 0);
  EXPECT_NE(material, nullptr);
  if (!material) return sums;

  PetscSection aux_section = nullptr;
  ierr = DMGetLocalSection(aux_dm, &aux_section);
  EXPECT_EQ(ierr, 0);

  const PetscScalar* aux_array = nullptr;
  ierr = VecGetArrayRead(aux_vec, &aux_array);
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

    PetscInt offset = 0;
    ierr = PetscSectionGetOffset(aux_section, c, &offset);
    EXPECT_EQ(ierr, 0);

    sums[material_id] += PetscRealPart(aux_array[offset + AUX_LAMBDA]);
    counts[material_id] += 1;
  }

  ierr = VecRestoreArrayRead(aux_vec, &aux_array);
  EXPECT_EQ(ierr, 0);

  for (auto& entry : sums)
  {
    entry.second /= static_cast<double>(counts[entry.first]);
  }
  return sums;
}

} // namespace

class GasbuggyMeshTest : public ::testing::Test
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
      cfg << "gmsh_label = overburden\n";
      cfg << "youngs_modulus = 20.0e9\n";
      cfg << "poisson_ratio = 0.30\n";
      cfg << "density = 2200.0\n";
      cfg << "\n[MATERIAL_REGION_2]\n";
      cfg << "gmsh_label = target\n";
      cfg << "youngs_modulus = 40.0e9\n";
      cfg << "poisson_ratio = 0.25\n";
      cfg << "density = 2500.0\n";
      cfg << "\n[MATERIAL_REGION_3]\n";
      cfg << "gmsh_label = basement\n";
      cfg << "youngs_modulus = 70.0e9\n";
      cfg << "poisson_ratio = 0.20\n";
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
  static constexpr const char* config_path_ = "test_gasbuggy_mesh.ini";
  static constexpr const char* mesh_path_ = "../../meshes/historical/gasbuggy_layered_3d.msh";
};

TEST_F(GasbuggyMeshTest, HistoricalMeshMapsThreeMaterialRegions)
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

  auto lambdas = collectAverageLambdaByMaterial(sim);
  ASSERT_EQ(lambdas.size(), 3u) << "Expected three mapped material regions on Gasbuggy mesh";

  double min_lambda = lambdas.begin()->second;
  double max_lambda = lambdas.begin()->second;
  for (const auto& entry : lambdas)
  {
    min_lambda = std::min(min_lambda, entry.second);
    max_lambda = std::max(max_lambda, entry.second);
  }
  EXPECT_GT(max_lambda, min_lambda);
  EXPECT_GT(max_lambda / std::max(min_lambda, 1.0e-12), 2.0)
      << "Expected a clear stiffness contrast across historical layers";

  ierr = sim.setupTimeStepper();
  ASSERT_EQ(ierr, 0);
  ierr = sim.setupSolvers();
  ASSERT_EQ(ierr, 0);
  ierr = sim.setInitialConditions();
  ASSERT_EQ(ierr, 0);
  ierr = sim.run();
  ASSERT_EQ(ierr, 0);

  PetscReal norm = 0.0;
  ierr = VecNorm(sim.getSolution(), NORM_2, &norm);
  ASSERT_EQ(ierr, 0);
  EXPECT_GT(norm, 0.0) << "Historical Gasbuggy solve should produce a nonzero displacement field";
}