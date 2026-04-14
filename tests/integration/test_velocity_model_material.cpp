/**
 * @file test_velocity_model_material.cpp
 * @brief Integration test: per-cell material properties from binary velocity model
 *
 * Generates a small synthetic two-layer velocity model, writes it to a
 * temporary binary file, configures the Simulator to read it, and verifies
 * that the auxiliary field values (lambda, mu, rho) differ between the top
 * and bottom halves of the mesh according to the velocity model.
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
#include "io/VelocityModelReader.hpp"
#include "numerics/PetscFEElasticityAux.hpp"

using namespace FSRM;

class VelocityModelMaterialTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank_);
  }

  int rank_ = 0;
};

TEST_F(VelocityModelMaterialTest, TwoLayerAuxFields)
{
  // -------------------------------------------------------------------
  // 1. Generate a two-layer velocity model binary file
  // -------------------------------------------------------------------
  std::string vel_path = "test_velocity_model.bin";
  std::string config_path = "test_velocity_model.config";

  // Layer parameters:
  // Top half  (z > 2500): Vp=3000, Vs=1732, rho=2200 (soft sediment)
  // Bottom half (z <= 2500): Vp=6000, Vs=3464, rho=2700 (hard basement)
  const float vp_top = 3000.0f, vs_top = 1732.0f, rho_top = 2200.0f;
  const float vp_bot = 6000.0f, vs_bot = 3464.0f, rho_bot = 2700.0f;

  if (rank_ == 0)
  {
    VelocityModel model;
    model.nx = 4;
    model.ny = 4;
    model.nz = 4;
    model.x_min = 0.0;
    model.x_max = 5000.0;
    model.y_min = 0.0;
    model.y_max = 5000.0;
    model.z_min = 0.0;
    model.z_max = 5000.0;

    size_t n = static_cast<size_t>(model.nx * model.ny * model.nz);
    model.vp.resize(n);
    model.vs.resize(n);
    model.rho.resize(n);

    for (int ix = 0; ix < model.nx; ++ix) {
      for (int iy = 0; iy < model.ny; ++iy) {
        for (int iz = 0; iz < model.nz; ++iz) {
          size_t idx = static_cast<size_t>(ix) * static_cast<size_t>(model.ny * model.nz)
                     + static_cast<size_t>(iy) * static_cast<size_t>(model.nz)
                     + static_cast<size_t>(iz);
          double z = model.z_min
                   + (model.z_max - model.z_min) * iz / (model.nz - 1);
          if (z > 2500.0) {
            model.vp[idx]  = vp_top;
            model.vs[idx]  = vs_top;
            model.rho[idx] = rho_top;
          } else {
            model.vp[idx]  = vp_bot;
            model.vs[idx]  = vs_bot;
            model.rho[idx] = rho_bot;
          }
        }
      }
    }

    int err = writeVelocityModel(vel_path, model);
    ASSERT_EQ(err, 0) << "Failed to write velocity model binary";

    // -------------------------------------------------------------------
    // 2. Write config file
    // -------------------------------------------------------------------
    std::ofstream cfg(config_path);
    cfg << "[SIMULATION]\n";
    cfg << "name = test_velocity_model\n";
    cfg << "start_time = 0.0\n";
    cfg << "end_time = 0.001\n";
    cfg << "dt_initial = 0.001\n";
    cfg << "dt_min = 0.0001\n";
    cfg << "dt_max = 0.001\n";
    cfg << "max_timesteps = 1\n";
    cfg << "output_frequency = 100\n";
    cfg << "fluid_model = NONE\n";
    cfg << "solid_model = ELASTIC\n";
    cfg << "enable_geomechanics = true\n";
    cfg << "enable_elastodynamics = false\n";
    cfg << "\n[GRID]\n";
    cfg << "nx = 4\n";
    cfg << "ny = 4\n";
    cfg << "nz = 4\n";
    cfg << "Lx = 5000.0\n";
    cfg << "Ly = 5000.0\n";
    cfg << "Lz = 5000.0\n";
    cfg << "\n[ROCK]\n";
    cfg << "density = 2650.0\n";
    cfg << "youngs_modulus = 10.0e9\n";
    cfg << "poissons_ratio = 0.25\n";
    cfg << "\n[MATERIAL]\n";
    cfg << "heterogeneous = true\n";
    cfg << "assignment = velocity_model\n";
    cfg << "velocity_model_file = " << vel_path << "\n";
    cfg << "\n[BOUNDARY_CONDITIONS]\n";
    cfg << "bottom = fixed\n";
    cfg << "sides = roller\n";
    cfg << "top = free\n";
    cfg.close();
  }
  MPI_Barrier(PETSC_COMM_WORLD);

  // -------------------------------------------------------------------
  // 3. Run the Simulator pipeline through setupPhysics (to populate aux)
  // -------------------------------------------------------------------
  Simulator sim(PETSC_COMM_WORLD);
  PetscErrorCode ierr;

  PetscOptionsClear(nullptr);
  PetscOptionsSetValue(nullptr, "-ts_type", "beuler");
  PetscOptionsSetValue(nullptr, "-pc_type", "lu");
  PetscOptionsSetValue(nullptr, "-ksp_type", "preonly");

  ierr = sim.initializeFromConfigFile(config_path);
  ASSERT_EQ(ierr, 0) << "initializeFromConfigFile must succeed";
  ierr = sim.setupDM();
  ASSERT_EQ(ierr, 0);
  ierr = sim.labelBoundaries();
  ASSERT_EQ(ierr, 0);
  ierr = sim.setupFields();
  ASSERT_EQ(ierr, 0);
  ierr = sim.setupPhysics();
  ASSERT_EQ(ierr, 0) << "setupPhysics must succeed with velocity_model assignment";

  // -------------------------------------------------------------------
  // 4. Verify auxiliary field values
  // -------------------------------------------------------------------
  DM dm = sim.getDM();
  DM auxDM = sim.getAuxDM();
  Vec auxVec = sim.getAuxVector();
  ASSERT_NE(auxDM, nullptr) << "Auxiliary DM must be created for velocity_model assignment";
  ASSERT_NE(auxVec, nullptr) << "Auxiliary vector must be created";

  PetscInt dim = 0;
  ierr = DMGetDimension(dm, &dim);
  ASSERT_EQ(ierr, 0);

  PetscSection section = nullptr;
  ierr = DMGetLocalSection(auxDM, &section);
  ASSERT_EQ(ierr, 0);

  PetscInt cStart = 0, cEnd = 0;
  ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd);
  ASSERT_EQ(ierr, 0);
  ASSERT_GT(cEnd - cStart, 0) << "Mesh must have cells";

  const PetscScalar *a = nullptr;
  ierr = VecGetArrayRead(auxVec, &a);
  ASSERT_EQ(ierr, 0);

  // Expected Lame parameters from velocity model:
  // Top:    mu = 2200 * 1732^2 = 6.601e9,  lambda = 2200*(3000^2 - 2*1732^2) = 6.599e9
  // Bottom: mu = 2700 * 3464^2 = 3.240e10, lambda = 2700*(6000^2 - 2*3464^2) = 3.239e10
  double expected_mu_top = static_cast<double>(rho_top) * static_cast<double>(vs_top) * static_cast<double>(vs_top);
  double expected_mu_bot = static_cast<double>(rho_bot) * static_cast<double>(vs_bot) * static_cast<double>(vs_bot);

  int cells_top = 0, cells_bot = 0;
  double mu_sum_top = 0.0, mu_sum_bot = 0.0;
  double rho_sum_top = 0.0, rho_sum_bot = 0.0;

  for (PetscInt c = cStart; c < cEnd; ++c) {
    // Skip cohesive cells
    DMPolytopeType ct;
    ierr = DMPlexGetCellType(dm, c, &ct);
    if (ct == DM_POLYTOPE_SEG_PRISM_TENSOR ||
        ct == DM_POLYTOPE_TRI_PRISM_TENSOR ||
        ct == DM_POLYTOPE_QUAD_PRISM_TENSOR) continue;

    PetscReal vol, centroid[3], normal[3];
    ierr = DMPlexComputeCellGeometryFVM(dm, c, &vol, centroid, normal);
    if (ierr) continue;

    PetscInt offset = 0;
    ierr = PetscSectionGetOffset(section, c, &offset);
    if (ierr) continue;

    double mu_val = PetscRealPart(a[offset + AUX_MU]);
    double rho_val = PetscRealPart(a[offset + AUX_RHO]);

    // z > 2500 is the soft top layer
    if (centroid[dim - 1] > 2500.0) {
      mu_sum_top += mu_val;
      rho_sum_top += rho_val;
      cells_top++;
    } else {
      mu_sum_bot += mu_val;
      rho_sum_bot += rho_val;
      cells_bot++;
    }
  }

  ierr = VecRestoreArrayRead(auxVec, &a);
  ASSERT_EQ(ierr, 0);

  if (rank_ == 0) {
    std::remove(vel_path.c_str());
    std::remove(config_path.c_str());
  }

  ASSERT_GT(cells_top, 0) << "Must have cells in the top layer";
  ASSERT_GT(cells_bot, 0) << "Must have cells in the bottom layer";

  double avg_mu_top = mu_sum_top / cells_top;
  double avg_mu_bot = mu_sum_bot / cells_bot;
  double avg_rho_top = rho_sum_top / cells_top;
  double avg_rho_bot = rho_sum_bot / cells_bot;

  // Top layer: soft sediment (lower mu, lower rho)
  // Bottom layer: hard basement (higher mu, higher rho)
  EXPECT_GT(avg_mu_bot, avg_mu_top * 2.0)
      << "Bottom layer mu should be significantly larger than top layer mu";
  EXPECT_GT(avg_rho_bot, avg_rho_top)
      << "Bottom layer rho should be larger than top layer rho";

  // Check that values are in the right ballpark (within 30% of expected,
  // accounting for trilinear interpolation across the interface)
  EXPECT_NEAR(avg_mu_top, expected_mu_top, expected_mu_top * 0.3)
      << "Top layer mu should be near expected value from Vs=1732, rho=2200";
  EXPECT_NEAR(avg_mu_bot, expected_mu_bot, expected_mu_bot * 0.3)
      << "Bottom layer mu should be near expected value from Vs=3464, rho=2700";
  EXPECT_NEAR(avg_rho_top, static_cast<double>(rho_top), 300.0)
      << "Top layer rho should be near 2200 kg/m^3";
  EXPECT_NEAR(avg_rho_bot, static_cast<double>(rho_bot), 300.0)
      << "Bottom layer rho should be near 2700 kg/m^3";
}
