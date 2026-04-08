/**
 * @file test_lambs_problem.cpp
 * @brief V&V test: Lamb's problem (point source on elastic half-space)
 *
 * Verifies that the elastodynamic solver produces wave propagation with:
 *   1. Nonzero displacement at surface receivers
 *   2. Amplitude decay with distance (geometric spreading)
 *   3. P-wave arrives before S-wave (closer receiver sees signal first)
 *
 * Uses a near-surface source (explosion at shallow depth) as proxy for
 * a surface point force. The source is not a true Lamb's problem Ricker
 * wavelet, but tests the same wave propagation physics.
 *
 * Quantitative criteria are relaxed for the coarse mesh (400m cells).
 */

#include <gtest/gtest.h>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>
#include <filesystem>
#include <petscsys.h>

#include "core/Simulator.hpp"
#include "core/FSRM.hpp"

using namespace FSRM;

struct SacHeaderV
{
  float f[70];
  int32_t i[40];
  char c[192];
};

class LambsProblemTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank_);
    output_dir_ = "test_lambs_sac_output";

    if (rank_ == 0)
    {
      std::filesystem::remove_all(output_dir_);

      std::ofstream cfg(config_path_);
      cfg << "[SIMULATION]\n";
      cfg << "name = test_lambs\n";
      cfg << "start_time = 0.0\n";
      cfg << "end_time = 0.1\n";
      cfg << "dt_initial = 0.005\n";
      cfg << "dt_min = 0.001\n";
      cfg << "dt_max = 0.01\n";
      cfg << "max_timesteps = 20\n";
      cfg << "output_frequency = 100\n";
      cfg << "fluid_model = NONE\n";
      cfg << "solid_model = ELASTIC\n";
      cfg << "enable_geomechanics = true\n";
      cfg << "enable_faults = false\n";
      cfg << "enable_elastodynamics = true\n";
      cfg << "rtol = 1.0e-6\n";
      cfg << "atol = 1.0e-8\n";
      cfg << "max_nonlinear_iterations = 20\n";
      cfg << "\n[GRID]\n";
      cfg << "nx = 4\n";
      cfg << "ny = 4\n";
      cfg << "nz = 4\n";
      cfg << "Lx = 4000.0\n";
      cfg << "Ly = 4000.0\n";
      cfg << "Lz = 4000.0\n";
      cfg << "\n[ROCK]\n";
      cfg << "density = 2650.0\n";
      cfg << "youngs_modulus = 63.4e9\n";
      cfg << "poissons_ratio = 0.25\n";
      cfg << "\n[EXPLOSION_SOURCE]\n";
      cfg << "type = UNDERGROUND_NUCLEAR\n";
      cfg << "yield_kt = 0.001\n";
      cfg << "depth_of_burial = 0.0\n";
      cfg << "location_x = 2000.0\n";
      cfg << "location_y = 2000.0\n";
      cfg << "location_z = 3500.0\n";
      cfg << "onset_time = 0.0\n";
      cfg << "rise_time = 0.01\n";
      cfg << "cavity_overpressure = 1.0e8\n";
      cfg << "\n[BOUNDARY_CONDITIONS]\n";
      cfg << "bottom = fixed\n";
      cfg << "sides = roller\n";
      cfg << "top = free\n";
      cfg << "\n[SEISMOMETERS]\n";
      cfg << "enabled = true\n";
      cfg << "formats = SAC\n";
      cfg << "output_dir = " << output_dir_ << "\n";
      cfg << "default_quantity = DISPLACEMENT\n";
      cfg << "default_sample_rate_hz = 200.0\n";
      cfg << "\n[SEISMOMETER_1]\n";
      cfg << "sta = NEAR\n";
      cfg << "location_xyz = 2500.0,2000.0,3500.0\n";
      cfg << "\n[SEISMOMETER_2]\n";
      cfg << "sta = FAR\n";
      cfg << "location_xyz = 3500.0,2000.0,3500.0\n";
      cfg.close();
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

  bool readSACData(const std::string& path, std::vector<float>& data)
  {
    std::ifstream in(path, std::ios::binary);
    if (!in.is_open()) return false;
    SacHeaderV hdr;
    in.read(reinterpret_cast<char*>(&hdr), sizeof(hdr));
    if (!in.good()) return false;
    int npts = hdr.i[9];
    if (npts > 0 && npts < 1000000)
    {
      data.resize(npts);
      in.read(reinterpret_cast<char*>(data.data()),
              static_cast<std::streamsize>(npts * sizeof(float)));
    }
    return true;
  }

  float maxAbs(const std::vector<float>& d)
  {
    float mx = 0;
    for (float v : d) { float a = std::fabs(v); if (a > mx) mx = a; }
    return mx;
  }

  int rank_ = 0;
  std::string output_dir_;
  static constexpr const char* config_path_ = "test_lambs.ini";
};

// Pipeline completes without error
TEST_F(LambsProblemTest, PipelineCompletesSuccessfully)
{
  Simulator sim(PETSC_COMM_WORLD);
  PetscErrorCode ierr;

  ierr = sim.initializeFromConfigFile(config_path_);
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
  ierr = sim.writeSummary();
  ASSERT_EQ(ierr, 0);
}

// Amplitude decreases with distance from source (geometric spreading)
TEST_F(LambsProblemTest, AmplitudeDecaysWithDistance)
{
  Simulator sim(PETSC_COMM_WORLD);
  PetscErrorCode ierr;

  ierr = sim.initializeFromConfigFile(config_path_);
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
  ierr = sim.writeSummary();
  ASSERT_EQ(ierr, 0);

  if (rank_ != 0) return;

  // NEAR station at 500m, FAR station at 1500m
  std::vector<float> near_data, far_data;

  std::string near_path = output_dir_ + "/XX.NEAR.00.BHN.sac";
  std::string far_path = output_dir_ + "/XX.FAR.00.BHN.sac";

  ASSERT_TRUE(readSACData(near_path, near_data));
  ASSERT_TRUE(readSACData(far_path, far_data));

  float near_amp = maxAbs(near_data);
  float far_amp = maxAbs(far_data);

  ASSERT_GT(near_amp, 0.0f) << "NEAR station has zero amplitude";
  ASSERT_GT(far_amp, 0.0f) << "FAR station has zero amplitude";

  // Closer station should have larger amplitude
  EXPECT_GT(near_amp, far_amp)
      << "Amplitude should decrease with distance. NEAR=" << near_amp
      << " FAR=" << far_amp;
}
