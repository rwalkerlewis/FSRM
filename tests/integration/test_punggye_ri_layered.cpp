/**
 * @file test_punggye_ri_layered.cpp
 * @brief Compact integration test for the Punggye-ri layered explosion workflow
 */

#include <gtest/gtest.h>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include <petscsys.h>

#include "core/FSRM.hpp"
#include "core/Simulator.hpp"

using namespace FSRM;

namespace
{

struct SacHeaderRead
{
  float f[70];
  int32_t i[40];
  char c[192];
};

bool readSACFile(const std::string& path, int& npts, std::vector<float>& data)
{
  std::ifstream in(path, std::ios::binary);
  if (!in.is_open()) return false;

  SacHeaderRead hdr;
  in.read(reinterpret_cast<char*>(&hdr), sizeof(hdr));
  if (!in.good()) return false;

  npts = hdr.i[9];
  if (npts <= 0 || npts >= 1000000) return false;

  data.resize(npts);
  in.read(reinterpret_cast<char*>(data.data()), static_cast<std::streamsize>(npts * sizeof(float)));
  return in.good() || in.eof();
}

float maxAbsAmplitude(const std::vector<float>& data)
{
  float mx = 0.0f;
  for (float value : data)
  {
    mx = std::max(mx, std::fabs(value));
  }
  return mx;
}

} // namespace

class PunggyeRiLayeredTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank_);
    sac_output_dir_ = "test_punggye_ri_sac";
    hdf5_output_dir_ = "test_punggye_ri_hdf5";

    if (rank_ == 0)
    {
      std::filesystem::remove_all(sac_output_dir_);
      std::filesystem::remove_all(hdf5_output_dir_);

      std::ofstream cfg(config_path_);
      cfg << "[SIMULATION]\n";
      cfg << "name = test_punggye_ri_layered\n";
      cfg << "start_time = 0.0\n";
      cfg << "end_time = 1.0\n";
      cfg << "dt_initial = 0.005\n";
      cfg << "dt_min = 0.001\n";
      cfg << "dt_max = 0.01\n";
      cfg << "max_timesteps = 200\n";
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
      cfg << "nx = 8\n";
      cfg << "ny = 8\n";
      cfg << "nz = 6\n";
      cfg << "Lx = 10000.0\n";
      cfg << "Ly = 10000.0\n";
      cfg << "Lz = 5000.0\n";
      cfg << "\n[ROCK]\n";
      cfg << "density = 2650.0\n";
      cfg << "youngs_modulus = 70.0e9\n";
      cfg << "poisson_ratio = 0.25\n";
      cfg << "\n[MATERIAL]\n";
      cfg << "heterogeneous = true\n";
      cfg << "gravity = 0.0\n";
      cfg << "\n[LAYER_1]\n";
      cfg << "z_top = 5000.0\n";
      cfg << "z_bottom = 4800.0\n";
      cfg << "lambda = 9.00e9\n";
      cfg << "mu = 8.98e9\n";
      cfg << "rho = 2200.0\n";
      cfg << "\n[LAYER_2]\n";
      cfg << "z_top = 4800.0\n";
      cfg << "z_bottom = 4000.0\n";
      cfg << "lambda = 1.69e10\n";
      cfg << "mu = 1.69e10\n";
      cfg << "rho = 2500.0\n";
      cfg << "\n[LAYER_3]\n";
      cfg << "z_top = 4000.0\n";
      cfg << "z_bottom = 0.0\n";
      cfg << "lambda = 2.97e10\n";
      cfg << "mu = 2.97e10\n";
      cfg << "rho = 2650.0\n";
      cfg << "\n[EXPLOSION_SOURCE]\n";
      cfg << "type = UNDERGROUND_NUCLEAR\n";
      cfg << "yield_kt = 250.0\n";
      cfg << "depth_of_burial = 800.0\n";
      cfg << "location_x = 5000.0\n";
      cfg << "location_y = 5000.0\n";
      cfg << "location_z = 4200.0\n";
      cfg << "onset_time = 0.0\n";
      cfg << "rise_time = 0.05\n";
      cfg << "cavity_overpressure = 1.0e10\n";
      cfg << "apply_damage_zone = true\n";
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
      cfg << "\n[OUTPUT]\n";
      cfg << "output_directory = " << hdf5_output_dir_ << "\n";
      cfg << "\n[SEISMOMETERS]\n";
      cfg << "enabled = true\n";
      cfg << "formats = SAC\n";
      cfg << "output_dir = " << sac_output_dir_ << "\n";
      cfg << "default_quantity = VELOCITY\n";
      cfg << "default_sample_rate_hz = 100.0\n";
      cfg << "\n[SEISMOMETER_1]\n";
      cfg << "sta = ABOVE\n";
      cfg << "location_xyz = 5000.0,5000.0,5000.0\n";
      cfg << "\n[SEISMOMETER_2]\n";
      cfg << "sta = EAST2\n";
      cfg << "location_xyz = 7000.0,5000.0,5000.0\n";
      cfg << "\n[SEISMOMETER_3]\n";
      cfg << "sta = EAST4\n";
      cfg << "location_xyz = 9000.0,5000.0,5000.0\n";
    }
    MPI_Barrier(PETSC_COMM_WORLD);
  }

  void TearDown() override
  {
    MPI_Barrier(PETSC_COMM_WORLD);
    if (rank_ == 0)
    {
      std::remove(config_path_);
      std::filesystem::remove_all(sac_output_dir_);
      std::filesystem::remove_all(hdf5_output_dir_);
    }
  }

  int rank_ = 0;
  std::string sac_output_dir_;
  std::string hdf5_output_dir_;
  static constexpr const char* config_path_ = "test_punggye_ri_layered.ini";
};

TEST_F(PunggyeRiLayeredTest, RunProducesSACAndHDF5Outputs)
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

  ierr = sim.writeSummary();
  ASSERT_EQ(ierr, 0);

  PetscReal norm = 0.0;
  ierr = VecNorm(sim.getSolution(), NORM_2, &norm);
  ASSERT_EQ(ierr, 0);
  EXPECT_GT(norm, 0.0);

  if (rank_ != 0) return;

  const std::string h5_path = hdf5_output_dir_ + "/solution.h5";
  ASSERT_TRUE(std::filesystem::exists(h5_path));
  EXPECT_GT(std::filesystem::file_size(h5_path), 0u);

  const std::vector<std::string> stations = {"ABOVE", "EAST2", "EAST4"};
  const std::vector<std::string> components = {"BHN", "BHE", "BHZ"};

  for (const auto& sta : stations)
  {
    for (const auto& comp : components)
    {
      const std::string sac_path = sac_output_dir_ + "/XX." + sta + ".00." + comp + ".sac";
      ASSERT_TRUE(std::filesystem::exists(sac_path)) << sac_path;

      int npts = 0;
      std::vector<float> data;
      ASSERT_TRUE(readSACFile(sac_path, npts, data));
      EXPECT_GT(npts, 1);
      EXPECT_GT(maxAbsAmplitude(data), 0.0f);
    }
  }

  int above_npts = 0, east2_npts = 0, east4_npts = 0;
  std::vector<float> above_data, east2_data, east4_data;
  ASSERT_TRUE(readSACFile(sac_output_dir_ + "/XX.ABOVE.00.BHZ.sac", above_npts, above_data));
  ASSERT_TRUE(readSACFile(sac_output_dir_ + "/XX.EAST2.00.BHE.sac", east2_npts, east2_data));
  ASSERT_TRUE(readSACFile(sac_output_dir_ + "/XX.EAST4.00.BHE.sac", east4_npts, east4_data));

  const float above_amp = maxAbsAmplitude(above_data);
  const float east2_amp = maxAbsAmplitude(east2_data);
  const float east4_amp = maxAbsAmplitude(east4_data);

  EXPECT_GT(above_amp, 0.0f);
  EXPECT_GT(east2_amp, 0.0f);
  EXPECT_GT(east4_amp, 0.0f);
  EXPECT_GT(east2_amp, east4_amp)
      << "Amplitude should decay between the 2 km and 4 km east stations";
}