/**
 * @file test_historic_nuclear.cpp
 * @brief Integration tests for 5 historic nuclear test configurations
 *
 * Each test verifies that the full Simulator pipeline completes for a
 * specific historic nuclear test with layered geology, Mueller-Murphy
 * explosion source, absorbing BCs, and SAC seismometer output.
 *
 * Configs are written inline with CI-friendly parameters (4x4x4 grid,
 * 0.1s simulation time) to keep each test under 30 seconds.
 *
 * Quantitative checks:
 *   1. Pipeline completes without error
 *   2. Solution norm is nonzero and finite
 *   3. SAC output files are produced with nonzero data
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

struct LayerDef
{
  double z_top;
  double z_bottom;
  double lambda;
  double mu;
  double rho;
};

class HistoricNuclearTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank_);
  }

  void TearDown() override
  {
    MPI_Barrier(PETSC_COMM_WORLD);
    if (rank_ == 0)
    {
      if (!config_path_.empty()) std::remove(config_path_.c_str());
      if (!output_dir_.empty()) std::filesystem::remove_all(output_dir_);
    }
  }

  // Write a CI-quick config for a historic nuclear test
  void writeConfig(const std::string& name, double yield_kt, double depth_m,
                   double domain_z, const std::vector<LayerDef>& layers)
  {
    config_path_ = "test_historic_" + name + ".config";
    output_dir_ = "test_historic_" + name + "_output";

    if (rank_ != 0) return;

    std::filesystem::remove_all(output_dir_);

    double loc_z = domain_z - depth_m;
    double half_xy = domain_z * 2.0;

    std::ofstream cfg(config_path_);
    cfg << "[SIMULATION]\n";
    cfg << "name = " << name << "\n";
    cfg << "start_time = 0.0\n";
    cfg << "end_time = 0.1\n";
    cfg << "dt_initial = 0.001\n";
    cfg << "dt_min = 0.0001\n";
    cfg << "dt_max = 0.005\n";
    cfg << "max_timesteps = 100\n";
    cfg << "output_frequency = 50\n";
    cfg << "output_format = HDF5\n";
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
    cfg << "Lx = " << half_xy << "\n";
    cfg << "Ly = " << half_xy << "\n";
    cfg << "Lz = " << domain_z << "\n";
    cfg << "\n[ROCK]\n";
    // Use deepest layer as fallback
    cfg << "density = " << layers.back().rho << "\n";
    cfg << "lambda = " << layers.back().lambda << "\n";
    cfg << "shear_modulus = " << layers.back().mu << "\n";
    cfg << "\n[MATERIAL]\n";
    cfg << "heterogeneous = true\n";
    cfg << "num_layers = " << layers.size() << "\n";
    for (size_t i = 0; i < layers.size(); ++i)
    {
      cfg << "\n[LAYER_" << (i + 1) << "]\n";
      cfg << "z_top = " << layers[i].z_top << "\n";
      cfg << "z_bottom = " << layers[i].z_bottom << "\n";
      cfg << "lambda = " << layers[i].lambda << "\n";
      cfg << "mu = " << layers[i].mu << "\n";
      cfg << "rho = " << layers[i].rho << "\n";
    }
    cfg << "\n[EXPLOSION_SOURCE]\n";
    cfg << "type = UNDERGROUND_NUCLEAR\n";
    cfg << "yield_kt = " << yield_kt << "\n";
    cfg << "depth_of_burial = " << depth_m << "\n";
    cfg << "location_x = " << half_xy / 2.0 << "\n";
    cfg << "location_y = " << half_xy / 2.0 << "\n";
    cfg << "location_z = " << loc_z << "\n";
    cfg << "onset_time = 0.0\n";
    cfg << "rise_time = 0.01\n";
    cfg << "cavity_overpressure = 1.0e10\n";
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
    cfg << "\n[SEISMOMETERS]\n";
    cfg << "enabled = true\n";
    cfg << "formats = SAC\n";
    cfg << "output_dir = " << output_dir_ << "\n";
    cfg << "default_quantity = DISPLACEMENT\n";
    cfg << "default_sample_rate_hz = 200.0\n";
    cfg << "\n[SEISMOMETER_1]\n";
    cfg << "sta = SPALL\n";
    cfg << "location_xyz = " << half_xy / 2.0 << ","
        << half_xy / 2.0 << "," << domain_z << "\n";
    cfg.close();
  }

  // Run the full pipeline for the current config
  PetscErrorCode runPipeline(PetscReal& sol_norm)
  {
    MPI_Barrier(PETSC_COMM_WORLD);

    Simulator sim(PETSC_COMM_WORLD);
    PetscErrorCode ierr;

    PetscOptionsClear(nullptr);

    ierr = sim.initializeFromConfigFile(config_path_);
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
    ierr = sim.setInitialConditions();
    if (ierr) return ierr;

    PetscPushErrorHandler(PetscReturnErrorHandler, nullptr);
    ierr = sim.run();
    PetscPopErrorHandler();
    if (ierr) return ierr;

    ierr = sim.writeSummary();
    if (ierr) return ierr;

    Vec sol = sim.getSolution();
    if (sol)
    {
      VecNorm(sol, NORM_2, &sol_norm);
    }
    else
    {
      sol_norm = -1.0;
    }
    return PETSC_SUCCESS;
  }

  // Check that at least one SAC file was produced with nonzero data
  bool checkSACOutput()
  {
    if (rank_ != 0) return true;

    std::vector<std::string> components = {"BHN", "BHE", "BHZ"};
    for (const auto& comp : components)
    {
      std::string sac_path = output_dir_ + "/XX.SPALL.00." + comp + ".sac";
      std::ifstream f(sac_path, std::ios::binary);
      if (f.good())
      {
        f.seekg(0, std::ios::end);
        std::streamsize size = f.tellg();
        if (size > 632)
        {
          return true;
        }
      }
    }
    return false;
  }

  int rank_ = 0;
  std::string config_path_;
  std::string output_dir_;
};

// Gasbuggy, NM (1967): 29 kt, 1280m, Lewis Shale
// Elastic moduli from: lambda = rho*(vp^2 - 2*vs^2), mu = rho*vs^2
TEST_F(HistoricNuclearTest, Gasbuggy1967)
{
  std::vector<LayerDef> layers = {
    {5000.0, 4800.0, 7.07e9,  3.02e9,  2100.0},  // Alluvium/sandstone
    {4800.0, 4200.0, 1.33e10, 1.30e10, 2450.0},  // Mesaverde Ss
    {4200.0, 3400.0, 1.12e10, 1.02e10, 2550.0},  // Lewis Shale
    {3400.0,    0.0, 1.68e10, 1.69e10, 2500.0},   // Pictured Cliffs Ss
  };
  writeConfig("gasbuggy_1967", 29.0, 1280.0, 5000.0, layers);

  PetscReal sol_norm = 0.0;
  PetscErrorCode ierr = runPipeline(sol_norm);

  ASSERT_EQ(ierr, 0) << "Gasbuggy 1967 (29 kt, Lewis Shale) pipeline must complete";
  EXPECT_GT(sol_norm, 0.0) << "Solution must be nonzero";
  EXPECT_TRUE(std::isfinite(sol_norm)) << "Solution norm must be finite";
  if (rank_ == 0)
  {
    EXPECT_TRUE(checkSACOutput()) << "SAC output files must be produced";
  }
}

// Gnome, NM (1961): 3.1 kt, 361m, Salado Salt
TEST_F(HistoricNuclearTest, Gnome1961)
{
  std::vector<LayerDef> layers = {
    {2000.0, 1900.0, 3.80e9,  2.25e9,  1900.0},  // Alluvium
    {1900.0, 1500.0, 2.56e10, 1.57e10, 2750.0},  // Rustler anhydrite
    {1500.0,  800.0, 1.82e10, 1.16e10, 2200.0},  // Salado Salt
    { 800.0,    0.0, 2.56e10, 1.57e10, 2750.0},   // Castile anhydrite
  };
  writeConfig("gnome_1961", 3.1, 361.0, 2000.0, layers);

  PetscReal sol_norm = 0.0;
  PetscErrorCode ierr = runPipeline(sol_norm);

  ASSERT_EQ(ierr, 0) << "Gnome 1961 (3.1 kt, Salado Salt) pipeline must complete";
  EXPECT_GT(sol_norm, 0.0) << "Solution must be nonzero";
  EXPECT_TRUE(std::isfinite(sol_norm)) << "Solution norm must be finite";
  if (rank_ == 0)
  {
    EXPECT_TRUE(checkSACOutput()) << "SAC output files must be produced";
  }
}

// Sedan, NV (1962): 104 kt, 194m, alluvium (shallow cratering)
TEST_F(HistoricNuclearTest, Sedan1962)
{
  std::vector<LayerDef> layers = {
    {2000.0, 1700.0, 4.36e9,  2.69e9,  1800.0},  // Alluvium
    {1700.0, 1000.0, 1.02e10, 7.50e9,  2300.0},  // Welded tuff
    {1000.0,    0.0, 1.87e10, 1.34e10, 2650.0},   // Paleozoic carbonate
  };
  writeConfig("sedan_1962", 104.0, 194.0, 2000.0, layers);

  PetscReal sol_norm = 0.0;
  PetscErrorCode ierr = runPipeline(sol_norm);

  ASSERT_EQ(ierr, 0) << "Sedan 1962 (104 kt, alluvium) pipeline must complete";
  EXPECT_GT(sol_norm, 0.0) << "Solution must be nonzero";
  EXPECT_TRUE(std::isfinite(sol_norm)) << "Solution norm must be finite";
  if (rank_ == 0)
  {
    EXPECT_TRUE(checkSACOutput()) << "SAC output files must be produced";
  }
}

// Degelen Mountain, Kazakhstan: 50 kt, 300m, granite
TEST_F(HistoricNuclearTest, DegelenMountain)
{
  std::vector<LayerDef> layers = {
    {3000.0, 2900.0, 8.40e9,  6.00e9,  2400.0},  // Weathered granite
    {2900.0, 2200.0, 1.82e10, 2.00e10, 2650.0},  // Competent granite
    {2200.0,    0.0, 2.40e10, 2.50e10, 2700.0},   // Deep granite
  };
  writeConfig("degelen_mountain", 50.0, 300.0, 3000.0, layers);

  PetscReal sol_norm = 0.0;
  PetscErrorCode ierr = runPipeline(sol_norm);

  ASSERT_EQ(ierr, 0) << "Degelen Mountain (50 kt, granite) pipeline must complete";
  EXPECT_GT(sol_norm, 0.0) << "Solution must be nonzero";
  EXPECT_TRUE(std::isfinite(sol_norm)) << "Solution norm must be finite";
  if (rank_ == 0)
  {
    EXPECT_TRUE(checkSACOutput()) << "SAC output files must be produced";
  }
}

// NTS Pahute Mesa, NV: 150 kt, 600m, tuff
TEST_F(HistoricNuclearTest, NtsPahuteMesa)
{
  std::vector<LayerDef> layers = {
    {3000.0, 2800.0, 4.36e9,  2.69e9,  1800.0},  // Alluvium/soil
    {2800.0, 2200.0, 5.10e9,  3.70e9,  2050.0},  // Nonwelded tuff
    {2200.0, 1400.0, 1.02e10, 7.50e9,  2300.0},  // Welded tuff
    {1400.0,    0.0, 1.87e10, 1.34e10, 2650.0},   // Paleozoic basement
  };
  writeConfig("nts_pahute_mesa", 150.0, 600.0, 3000.0, layers);

  PetscReal sol_norm = 0.0;
  PetscErrorCode ierr = runPipeline(sol_norm);

  ASSERT_EQ(ierr, 0) << "NTS Pahute Mesa (150 kt, tuff) pipeline must complete";
  EXPECT_GT(sol_norm, 0.0) << "Solution must be nonzero";
  EXPECT_TRUE(std::isfinite(sol_norm)) << "Solution norm must be finite";
  if (rank_ == 0)
  {
    EXPECT_TRUE(checkSACOutput()) << "SAC output files must be produced";
  }
}
