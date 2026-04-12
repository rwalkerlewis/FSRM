/**
 * @file test_explosion_seismogram.cpp
 * @brief Integration test: explosion source end-to-end with seismogram output
 *
 * Verifies the full pipeline: explosion moment tensor injection, TSALPHA2
 * second-order time integration, DMInterpolation sampling at seismometer
 * stations, and SAC file output.
 *
 * Quantitative checks:
 *   1. SAC files are produced for each station and each component
 *   2. SAC data contains nonzero samples (explosion produces displacement)
 *   3. Closer station records larger amplitudes (wave amplitude decay)
 */

#include <gtest/gtest.h>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>
#include <filesystem>
#include <petscsys.h>

#include "core/Simulator.hpp"
#include "core/FSRM.hpp"

using namespace FSRM;

// Minimal SAC header for reading back
struct SacHeaderRead
{
  float f[70];
  int32_t i[40];
  char c[192];
};

class ExplosionSeismogramTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank_);

    // Use a unique output directory to avoid conflicts
    output_dir_ = "test_explosion_sac_output";

    if (rank_ == 0)
    {
      // Clean up any previous output
      std::filesystem::remove_all(output_dir_);

      std::ofstream cfg(config_path_);
      cfg << "[SIMULATION]\n";
      cfg << "name = test_explosion\n";
      cfg << "start_time = 0.0\n";
      cfg << "end_time = 0.05\n";
      cfg << "dt_initial = 0.005\n";
      cfg << "dt_min = 0.001\n";
      cfg << "dt_max = 0.01\n";
      cfg << "max_timesteps = 10\n";
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
      cfg << "youngs_modulus = 60.0e9\n";
      cfg << "poissons_ratio = 0.25\n";
      cfg << "\n[EXPLOSION_SOURCE]\n";
      cfg << "type = UNDERGROUND_NUCLEAR\n";
      cfg << "yield_kt = 10.0\n";
      cfg << "depth_of_burial = 500.0\n";
      cfg << "location_x = 2000.0\n";
      cfg << "location_y = 2000.0\n";
      cfg << "location_z = 2000.0\n";
      cfg << "onset_time = 0.0\n";
      cfg << "rise_time = 0.01\n";
      cfg << "cavity_overpressure = 1.0e10\n";
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
      cfg << "location_xyz = 2500.0,2000.0,2000.0\n";
      cfg << "\n[SEISMOMETER_2]\n";
      cfg << "sta = FAR\n";
      cfg << "location_xyz = 3500.0,2000.0,2000.0\n";
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

  // Read a SAC file and return the data samples and npts
  bool readSACFile(const std::string& path, int& npts, float& depmax,
                   std::vector<float>& data)
  {
    std::ifstream in(path, std::ios::binary);
    if (!in.is_open()) return false;

    SacHeaderRead hdr;
    in.read(reinterpret_cast<char*>(&hdr), sizeof(hdr));
    if (!in.good()) return false;

    npts = hdr.i[9];  // NPTS field
    depmax = hdr.f[2]; // DEPMAX field

    if (npts > 0 && npts < 1000000)
    {
      data.resize(npts);
      in.read(reinterpret_cast<char*>(data.data()),
              static_cast<std::streamsize>(npts * sizeof(float)));
    }
    return in.good() || in.eof();
  }

  float maxAbsAmplitude(const std::vector<float>& data)
  {
    float mx = 0.0f;
    for (float v : data)
    {
      float a = std::fabs(v);
      if (a > mx) mx = a;
    }
    return mx;
  }

  int rank_ = 0;
  std::string output_dir_;
  static constexpr const char* config_path_ = "test_explosion_seis.ini";
};

// Test 1: Full explosion pipeline completes without error
TEST_F(ExplosionSeismogramTest, PipelineCompletesWithoutError)
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

  ierr = sim.run();
  ASSERT_EQ(ierr, 0) << "run() failed";

  ierr = sim.writeSummary();
  ASSERT_EQ(ierr, 0) << "writeSummary failed";
}

// Test 2: SAC files are produced with nonzero data
TEST_F(ExplosionSeismogramTest, SACFilesProducedWithNonzeroData)
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

  if (rank_ != 0) return;  // SAC files written by rank 0 only

  // Check that SAC files exist for both stations, all 3 components
  const std::vector<std::string> components = {"BHN", "BHE", "BHZ"};
  const std::vector<std::string> stations = {"NEAR", "FAR"};

  for (const auto& sta : stations)
  {
    for (const auto& comp : components)
    {
      std::string sac_path = output_dir_ + "/XX." + sta + ".00." + comp + ".sac";
      ASSERT_TRUE(std::filesystem::exists(sac_path))
          << "SAC file missing: " << sac_path;

      int npts = 0;
      float depmax = 0.0f;
      std::vector<float> data;
      ASSERT_TRUE(readSACFile(sac_path, npts, depmax, data))
          << "Failed to read SAC file: " << sac_path;

      // Must have at least 2 data points (need 2+ samples for SAC write)
      EXPECT_GT(npts, 1) << "SAC file has too few samples: " << sac_path;

      // Data must be nonzero (explosion creates displacement)
      float max_amp = maxAbsAmplitude(data);
      EXPECT_GT(max_amp, 0.0f)
          << "SAC data is all zeros for " << sac_path;
    }
  }
}

// Test 3: Closer station records larger amplitude (wave propagation check)
TEST_F(ExplosionSeismogramTest, AmplitudeDecreasesWithDistance)
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

  // Read the radial component (BHN = x-direction) for both stations
  // NEAR is at 500m from source, FAR is at 1500m
  std::string near_path = output_dir_ + "/XX.NEAR.00.BHN.sac";
  std::string far_path = output_dir_ + "/XX.FAR.00.BHN.sac";

  int near_npts = 0, far_npts = 0;
  float near_depmax = 0, far_depmax = 0;
  std::vector<float> near_data, far_data;

  ASSERT_TRUE(readSACFile(near_path, near_npts, near_depmax, near_data));
  ASSERT_TRUE(readSACFile(far_path, far_npts, far_depmax, far_data));

  float near_amp = maxAbsAmplitude(near_data);
  float far_amp = maxAbsAmplitude(far_data);

  // Both must be nonzero
  ASSERT_GT(near_amp, 0.0f) << "NEAR station has zero amplitude";
  ASSERT_GT(far_amp, 0.0f) << "FAR station has zero amplitude";

  // Closer station should have larger amplitude
  // For body waves in 3D, amplitude decays as 1/r, so at 3x distance
  // amplitude should be ~1/3. With coarse mesh, just check > relation.
  EXPECT_GT(near_amp, far_amp)
      << "NEAR amplitude (" << near_amp << ") should exceed FAR amplitude ("
      << far_amp << ") due to geometric spreading";
}
