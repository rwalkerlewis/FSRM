/**
 * @file test_nearfield_coupled.cpp
 * @brief Integration test: COUPLED_ANALYTIC NearField-FEM source coupling
 *
 * Verifies that the COUPLED_ANALYTIC explosion mode (RDP-based source time
 * function from the 1D NearFieldExplosionSolver) produces valid seismograms
 * that differ from the default PROXY (Mueller-Murphy) mode.
 *
 * Quantitative checks:
 *   1. Both PROXY and COUPLED_ANALYTIC complete successfully
 *   2. Both produce nonzero seismograms
 *   3. Seismograms differ (different source time functions)
 *   4. COUPLED_ANALYTIC near-station amplitude differs from PROXY by >10%
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

struct SacHeaderNFC
{
  float f[70];
  int32_t i[40];
  char c[192];
};

class NearFieldCoupledTest : public ::testing::Test
{
protected:
  void writeConfig(const std::string& path, const std::string& mode,
                   const std::string& out_dir)
  {
    std::ofstream cfg(path);
    cfg << "[SIMULATION]\n";
    cfg << "name = test_nf_" << mode << "\n";
    cfg << "start_time = 0.0\n";
    cfg << "end_time = 0.5\n";
    cfg << "dt_initial = 0.0025\n";
    cfg << "dt_min = 0.001\n";
    cfg << "dt_max = 0.005\n";
    cfg << "max_timesteps = 200\n";
    cfg << "output_frequency = 100\n";
    cfg << "fluid_model = NONE\n";
    cfg << "solid_model = ELASTIC\n";
    cfg << "enable_geomechanics = true\n";
    cfg << "enable_faults = false\n";
    cfg << "enable_elastodynamics = true\n";
    cfg << "enable_near_field_damage = true\n";
    cfg << "explosion_solve_mode = " << mode << "\n";
    cfg << "rtol = 1.0e-6\n";
    cfg << "atol = 1.0e-8\n";
    cfg << "max_nonlinear_iterations = 20\n";
    cfg << "\n[GRID]\n";
    cfg << "nx = 4\nny = 4\nnz = 4\n";
    cfg << "Lx = 4000.0\nLy = 4000.0\nLz = 4000.0\n";
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
    cfg << "apply_damage_zone = true\n";
    cfg << "\n[BOUNDARY_CONDITIONS]\n";
    cfg << "bottom = free\nsides = free\ntop = free\n";
    cfg << "\n[ABSORBING_BC]\n";
    cfg << "enabled = true\n";
    cfg << "x_min = true\nx_max = true\n";
    cfg << "y_min = true\ny_max = true\n";
    cfg << "z_min = true\nz_max = false\n";
    cfg << "\n[SEISMOMETERS]\n";
    cfg << "enabled = true\n";
    cfg << "formats = SAC\n";
    cfg << "output_dir = " << out_dir << "\n";
    cfg << "default_quantity = DISPLACEMENT\n";
    cfg << "default_sample_rate_hz = 200.0\n";
    cfg << "\n[SEISMOMETER_1]\n";
    cfg << "sta = NEAR\n";
    cfg << "location_xyz = 2100.0,2000.0,2000.0\n";
    cfg << "\n[SEISMOMETER_2]\n";
    cfg << "sta = FAR\n";
    cfg << "location_xyz = 2500.0,2000.0,2000.0\n";
    cfg.close();
  }

  PetscErrorCode runSimulation(const std::string& config_path)
  {
    Simulator sim(PETSC_COMM_WORLD);
    PetscErrorCode ierr;

    ierr = sim.initializeFromConfigFile(config_path); CHKERRQ(ierr);
    ierr = sim.setupDM(); CHKERRQ(ierr);
    ierr = sim.labelBoundaries(); CHKERRQ(ierr);
    ierr = sim.setupFields(); CHKERRQ(ierr);
    ierr = sim.setupPhysics(); CHKERRQ(ierr);
    ierr = sim.setupTimeStepper(); CHKERRQ(ierr);
    ierr = sim.setupSolvers(); CHKERRQ(ierr);
    ierr = sim.setInitialConditions(); CHKERRQ(ierr);
    ierr = sim.run(); CHKERRQ(ierr);
    ierr = sim.writeSummary(); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  bool readSACFile(const std::string& path, int& npts,
                   std::vector<float>& data)
  {
    std::ifstream in(path, std::ios::binary);
    if (!in.is_open()) return false;
    SacHeaderNFC hdr;
    in.read(reinterpret_cast<char*>(&hdr), sizeof(hdr));
    if (!in.good()) return false;
    npts = hdr.i[9];
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

  void SetUp() override
  {
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank_);
    proxy_dir_ = "test_nf_proxy_output";
    coupled_dir_ = "test_nf_coupled_output";
    if (rank_ == 0)
    {
      std::filesystem::remove_all(proxy_dir_);
      std::filesystem::remove_all(coupled_dir_);
      writeConfig(proxy_config_, "PROXY", proxy_dir_);
      writeConfig(coupled_config_, "COUPLED_ANALYTIC", coupled_dir_);
    }
    MPI_Barrier(PETSC_COMM_WORLD);
  }

  void TearDown() override
  {
    MPI_Barrier(PETSC_COMM_WORLD);
    if (rank_ == 0)
    {
      std::remove(proxy_config_);
      std::remove(coupled_config_);
      std::filesystem::remove_all(proxy_dir_);
      std::filesystem::remove_all(coupled_dir_);
    }
  }

  int rank_ = 0;
  std::string proxy_dir_;
  std::string coupled_dir_;
  static constexpr const char* proxy_config_ = "test_nf_proxy.ini";
  static constexpr const char* coupled_config_ = "test_nf_coupled.ini";
};

// Both PROXY and COUPLED_ANALYTIC complete and produce different seismograms
TEST_F(NearFieldCoupledTest, CoupledAnalyticProducesDifferentSource)
{
  // Run PROXY mode
  PetscErrorCode ierr = runSimulation(proxy_config_);
  ASSERT_EQ(ierr, 0) << "PROXY simulation failed";

  // Run COUPLED_ANALYTIC mode
  ierr = runSimulation(coupled_config_);
  ASSERT_EQ(ierr, 0) << "COUPLED_ANALYTIC simulation failed";

  if (rank_ != 0) return;

  // Read SAC files from both runs
  std::string proxy_near = proxy_dir_ + "/XX.NEAR.00.BHN.sac";
  std::string coupled_near = coupled_dir_ + "/XX.NEAR.00.BHN.sac";

  ASSERT_TRUE(std::filesystem::exists(proxy_near))
      << "PROXY SAC file missing: " << proxy_near;
  ASSERT_TRUE(std::filesystem::exists(coupled_near))
      << "COUPLED_ANALYTIC SAC file missing: " << coupled_near;

  int proxy_npts = 0, coupled_npts = 0;
  std::vector<float> proxy_data, coupled_data;

  ASSERT_TRUE(readSACFile(proxy_near, proxy_npts, proxy_data));
  ASSERT_TRUE(readSACFile(coupled_near, coupled_npts, coupled_data));

  float proxy_amp = maxAbsAmplitude(proxy_data);
  float coupled_amp = maxAbsAmplitude(coupled_data);

  // Both must produce nonzero seismograms
  EXPECT_GT(proxy_amp, 0.0f) << "PROXY produced zero-amplitude seismogram";
  EXPECT_GT(coupled_amp, 0.0f) << "COUPLED_ANALYTIC produced zero-amplitude seismogram";

  // Peak amplitudes must differ by more than 10% (different source functions)
  float relative_diff = std::fabs(proxy_amp - coupled_amp)
                        / std::max(proxy_amp, coupled_amp);
  EXPECT_GT(relative_diff, 0.10f)
      << "PROXY amp=" << proxy_amp << ", COUPLED amp=" << coupled_amp
      << ", diff=" << relative_diff * 100.0f << "% (expected >10%)";
}
