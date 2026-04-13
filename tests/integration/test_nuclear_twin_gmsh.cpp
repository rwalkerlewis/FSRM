/**
 * @file test_nuclear_twin_gmsh.cpp
 * @brief Integration test for a compact Gmsh-based nuclear twin workflow
 */

#include <gtest/gtest.h>

#include <cstdio>
#include <filesystem>
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

std::map<PetscInt, double> averageLambdaByMaterial(Simulator& sim)
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
    PetscInt mat_id = -1;
    ierr = DMLabelGetValue(material, c, &mat_id);
    EXPECT_EQ(ierr, 0);
    if (mat_id < 0) continue;

    PetscInt offset = 0;
    ierr = PetscSectionGetOffset(aux_section, c, &offset);
    EXPECT_EQ(ierr, 0);

    sums[mat_id] += PetscRealPart(aux_array[offset + AUX_LAMBDA]);
    counts[mat_id] += 1;
  }

  ierr = VecRestoreArrayRead(aux_vec, &aux_array);
  EXPECT_EQ(ierr, 0);

  for (auto& kv : sums)
  {
    kv.second /= static_cast<double>(counts[kv.first]);
  }

  return sums;
}

} // namespace

class NuclearTwinGmshTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank_);
    output_dir_ = "test_nuclear_twin_gmsh_hdf5";

    if (rank_ == 0)
    {
      std::filesystem::remove_all(output_dir_);

      std::string mesh_path = "../meshes/test_two_material.msh";
      if (!std::filesystem::exists(mesh_path))
      {
        mesh_path = "../../meshes/test_two_material.msh";
      }

      std::ofstream cfg(config_path_);
      cfg << "[SIMULATION]\n";
      cfg << "name = test_nuclear_twin_gmsh\n";
      cfg << "start_time = 0.0\n";
      cfg << "end_time = 0.2\n";
      cfg << "dt_initial = 0.005\n";
      cfg << "dt_min = 0.001\n";
      cfg << "dt_max = 0.01\n";
      cfg << "max_timesteps = 100\n";
      cfg << "output_frequency = 50\n";
      cfg << "output_format = HDF5\n";
      cfg << "fluid_model = NONE\n";
      cfg << "solid_model = ELASTIC\n";
      cfg << "enable_geomechanics = true\n";
      cfg << "enable_elastodynamics = true\n";
      cfg << "enable_faults = false\n";
      cfg << "rtol = 1.0e-6\n";
      cfg << "atol = 1.0e-8\n";
      cfg << "\n[GRID]\n";
      cfg << "mesh_type = GMSH\n";
      cfg << "mesh_file = " << mesh_path << "\n";
      cfg << "\n[ROCK]\n";
      cfg << "density = 2500.0\n";
      cfg << "youngs_modulus = 40.0e9\n";
      cfg << "poisson_ratio = 0.25\n";
      cfg << "\n[MATERIAL]\n";
      cfg << "heterogeneous = true\n";
      cfg << "assignment = gmsh_label\n";
      cfg << "\n[MATERIAL_REGION_1]\n";
      cfg << "gmsh_label = soft_rock\n";
      cfg << "youngs_modulus = 12.0e9\n";
      cfg << "poissons_ratio = 0.25\n";
      cfg << "density = 2200.0\n";
      cfg << "\n[MATERIAL_REGION_2]\n";
      cfg << "gmsh_label = hard_rock\n";
      cfg << "youngs_modulus = 75.0e9\n";
      cfg << "poissons_ratio = 0.25\n";
      cfg << "density = 2800.0\n";
      cfg << "\n[EXPLOSION_SOURCE]\n";
      cfg << "type = UNDERGROUND_NUCLEAR\n";
      cfg << "yield_kt = 5.0\n";
      cfg << "depth_of_burial = 500.0\n";
      cfg << "location_x = 500.0\n";
      cfg << "location_y = 500.0\n";
      cfg << "location_z = -500.0\n";
      cfg << "onset_time = 0.0\n";
      cfg << "rise_time = 0.01\n";
      cfg << "cavity_overpressure = 1.0e10\n";
      cfg << "apply_damage_zone = true\n";
      cfg << "\n[BOUNDARY_CONDITIONS]\n";
      cfg << "bottom = free\n";
      cfg << "sides = free\n";
      cfg << "top = free\n";
      cfg << "\n[OUTPUT]\n";
      cfg << "output_directory = " << output_dir_ << "\n";
    }
    MPI_Barrier(PETSC_COMM_WORLD);
  }

  void TearDown() override
  {
    MPI_Barrier(PETSC_COMM_WORLD);
    if (rank_ == 0)
    {
      std::remove(config_path_);
      std::filesystem::remove_all(output_dir_);
    }
  }

  int rank_ = 0;
  std::string output_dir_;
  static constexpr const char* config_path_ = "test_nuclear_twin_gmsh.ini";
};

TEST_F(NuclearTwinGmshTest, GmshMaterialExplosionPipelineRuns)
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

  auto lambdas = averageLambdaByMaterial(sim);
  ASSERT_EQ(lambdas.size(), 2u);
  double min_lambda = lambdas.begin()->second;
  double max_lambda = lambdas.begin()->second;
  for (const auto& kv : lambdas)
  {
    min_lambda = std::min(min_lambda, kv.second);
    max_lambda = std::max(max_lambda, kv.second);
  }
  EXPECT_GT(max_lambda / std::max(min_lambda, 1.0e-12), 3.0);

  ierr = sim.run();
  ASSERT_EQ(ierr, 0);

  PetscReal norm = 0.0;
  ierr = VecNorm(sim.getSolution(), NORM_2, &norm);
  ASSERT_EQ(ierr, 0);
  EXPECT_GT(norm, 0.0);

  if (rank_ != 0) return;

  const std::string h5_path = output_dir_ + "/solution.h5";
  ASSERT_TRUE(std::filesystem::exists(h5_path));
  EXPECT_GT(std::filesystem::file_size(h5_path), 0u);
}
