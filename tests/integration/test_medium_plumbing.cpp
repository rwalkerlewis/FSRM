/**
 * @file test_medium_plumbing.cpp
 * @brief Integration test that EXPLOSION_SOURCE.medium_type plumbs end
 *        to end through the FEM residual.
 *
 * The pass-2 medium-type plumbing route is:
 *   config -> ExplosionCoupling::medium_type
 *          -> parseMediumType helper
 *          -> MuellerMurphySource::setMedium
 *          -> NuclearSourceParameters::scalar_moment(medium)
 *          -> moment-rate residual in addExplosionSourceToResidual.
 *
 * The default explosion_solve_mode is COUPLED_ANALYTIC, which uses the
 * RDP-source branch rather than the Mueller-Murphy proxy. To exercise
 * the medium-aware path these tests force explosion_solve_mode = PROXY.
 *
 * SAC-peak verification approach (chosen over a Simulator getter):
 *   - Two runs at the same yield (1 kt), the same depth, the same mesh
 *     and the same uniform medium_properties, differing only in the
 *     EXPLOSION_SOURCE.medium_type string.
 *   - The Mueller-Murphy scalar moment is M0 = 4*pi*rho*vp^2*Rc^3 with
 *     Rc = C(medium) * cbrt(W) * cbrt(rho_ref / rho). Thus
 *     M0(SALT)/M0(GRANITE) = (16/11)^3 ~ 3.073 at fixed density.
 *   - In the linear-response regime the FEM peak vertical displacement
 *     scales linearly with M0 (the moment tensor enters the residual
 *     linearly, the elastic operator is linear, no nonlinearity is
 *     introduced at these CI yields/depths). So the SAC peak ratio
 *     must match the cube-of-coefficient ratio within the mesh
 *     tolerance.
 *   - Tolerance: factor 0.7-1.3 of (16/11)^3 ~ 3.07 (i.e. a peak ratio
 *     in [2.15, 4.0]) absorbs CI mesh-cell quadrature noise.
 *
 * The unknown-medium fallback test runs a third configuration with
 * medium_type = "MARS_REGOLITH" and asserts (a) the run completes
 * successfully (no crash on unknown string), (b) peak amplitude is
 * within 0.1% of the GENERIC-explicit run (parseMediumType warns once
 * and returns GENERIC).
 */

#include <gtest/gtest.h>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <filesystem>
#include <string>
#include <vector>
#include <petscsys.h>

#include "core/Simulator.hpp"
#include "core/FSRM.hpp"
#include "sac_test_reader.hpp"

using namespace FSRM;

class MediumPlumbingTest : public ::testing::Test
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
      for (auto& cfg : created_configs_) std::remove(cfg.c_str());
      for (auto& dir : created_dirs_) std::filesystem::remove_all(dir);
    }
  }

  // Write a minimal CI config: 4x4x4 mesh, 1 kt, 300 m depth, single
  // bulk material (rho = 2650, vp = 5500, vs = 3100 -- same parameters
  // for every medium so the only difference between runs is the
  // medium_type string used by MuellerMurphySource::setMedium).
  std::string writeConfig(const std::string& tag,
                          const std::string& medium_type)
  {
    const std::string config_path = "test_medium_" + tag + ".config";
    const std::string output_dir = "test_medium_" + tag + "_output";

    if (rank_ != 0) {
      created_configs_.push_back(config_path);
      created_dirs_.push_back(output_dir);
      return config_path;
    }

    std::filesystem::remove_all(output_dir);

    const double half_xy = 4000.0;
    const double domain_z = 1000.0;
    const double yield_kt = 1.0;
    const double depth_m = 300.0;
    const double end_time = 0.05;

    std::ofstream cfg(config_path);
    cfg << "[SIMULATION]\n";
    cfg << "name = medium_" << tag << "\n";
    cfg << "start_time = 0.0\n";
    cfg << "end_time = " << end_time << "\n";
    cfg << "dt_initial = 0.001\n";
    cfg << "dt_min = 0.0001\n";
    cfg << "dt_max = 0.005\n";
    cfg << "max_timesteps = 200\n";
    cfg << "output_frequency = 100\n";
    cfg << "output_format = HDF5\n";
    cfg << "fluid_model = NONE\n";
    cfg << "solid_model = ELASTIC\n";
    cfg << "enable_geomechanics = true\n";
    cfg << "enable_faults = false\n";
    cfg << "enable_elastodynamics = true\n";
    // Force PROXY so the residual path actually invokes
    // MuellerMurphySource::setMedium. Default COUPLED_ANALYTIC bypasses
    // the medium-plumbing branch entirely (rdp_source is used instead).
    cfg << "explosion_solve_mode = PROXY\n";
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
    cfg << "density = 2650.0\n";
    cfg << "lambda = 2.5e10\n";
    cfg << "shear_modulus = 2.5e10\n";
    cfg << "\n[EXPLOSION_SOURCE]\n";
    cfg << "type = UNDERGROUND_NUCLEAR\n";
    cfg << "yield_kt = " << yield_kt << "\n";
    cfg << "depth_of_burial = " << depth_m << "\n";
    cfg << "location_x = " << half_xy / 2.0 << "\n";
    cfg << "location_y = " << half_xy / 2.0 << "\n";
    cfg << "location_z = " << (domain_z - depth_m) << "\n";
    cfg << "onset_time = 0.0\n";
    cfg << "rise_time = 0.01\n";
    cfg << "cavity_overpressure = 1.0e10\n";
    cfg << "medium_type = " << medium_type << "\n";
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
    cfg << "output_dir = " << output_dir << "\n";
    cfg << "default_quantity = DISPLACEMENT\n";
    cfg << "default_sample_rate_hz = 200.0\n";
    cfg << "\n[SEISMOMETER_1]\n";
    cfg << "sta = SPALL\n";
    cfg << "location_xyz = " << half_xy / 2.0 << ","
        << half_xy / 2.0 << "," << domain_z << "\n";
    cfg.close();

    created_configs_.push_back(config_path);
    created_dirs_.push_back(output_dir);
    return config_path;
  }

  // Run the full Simulator pipeline for the current config and return
  // the absolute peak BHZ amplitude at the spall station (rank 0 only).
  double runAndPeak(const std::string& config_path,
                    const std::string& tag)
  {
    MPI_Barrier(PETSC_COMM_WORLD);

    Simulator sim(PETSC_COMM_WORLD);
    PetscOptionsClear(nullptr);

    PetscErrorCode ierr = sim.initializeFromConfigFile(config_path);
    EXPECT_EQ(ierr, 0) << tag << ": initializeFromConfigFile failed";
    if (ierr) return 0.0;

    ierr = sim.setupDM();           EXPECT_EQ(ierr, 0); if (ierr) return 0.0;
    ierr = sim.labelBoundaries();   EXPECT_EQ(ierr, 0); if (ierr) return 0.0;
    ierr = sim.setupFields();       EXPECT_EQ(ierr, 0); if (ierr) return 0.0;
    ierr = sim.setupPhysics();      EXPECT_EQ(ierr, 0); if (ierr) return 0.0;
    ierr = sim.setupTimeStepper();  EXPECT_EQ(ierr, 0); if (ierr) return 0.0;
    ierr = sim.setupSolvers();      EXPECT_EQ(ierr, 0); if (ierr) return 0.0;
    ierr = sim.setInitialConditions(); EXPECT_EQ(ierr, 0); if (ierr) return 0.0;

    PetscPushErrorHandler(PetscReturnErrorHandler, nullptr);
    ierr = sim.run();
    PetscPopErrorHandler();
    EXPECT_EQ(ierr, 0) << tag << ": run() failed";
    if (ierr) return 0.0;

    ierr = sim.writeSummary();
    EXPECT_EQ(ierr, 0);
    if (ierr) return 0.0;

    if (rank_ != 0) return 0.0;

    const std::string sac_path =
        "test_medium_" + tag + "_output/XX.SPALL.00.BHZ.sac";
    auto trace = FSRM::test_helpers::readSAC(sac_path);
    EXPECT_TRUE(trace.valid)
        << tag << ": SAC trace must parse (" << sac_path << ")";
    if (!trace.valid) return 0.0;

    auto peak = FSRM::test_helpers::tracePeak(trace);
    return peak.peak_abs;
  }

  int rank_ = 0;
  std::vector<std::string> created_configs_;
  std::vector<std::string> created_dirs_;
};

// SaltExceedsGranite: same yield, same depth, same mesh, same bulk
// material; only EXPLOSION_SOURCE.medium_type differs. Salt cavity
// coefficient C = 16, granite C = 11. Cube-of-ratio (16/11)^3 ~ 3.07
// is the analytic M0 ratio. The SAC peak ratio must land in [2.15,
// 4.0] -- factor 0.7-1.3 of the analytic ratio, which absorbs CI
// mesh-cell quadrature noise.
TEST_F(MediumPlumbingTest, SaltExceedsGranite)
{
  const std::string cfg_g = writeConfig("granite", "GRANITE");
  const std::string cfg_s = writeConfig("salt", "SALT");

  const double peak_g = runAndPeak(cfg_g, "granite");
  const double peak_s = runAndPeak(cfg_s, "salt");

  if (rank_ != 0) return;

  ASSERT_GT(peak_g, 0.0)
      << "GRANITE peak BHZ amplitude must be positive";
  ASSERT_GT(peak_s, 0.0)
      << "SALT peak BHZ amplitude must be positive";

  const double ratio = peak_s / peak_g;
  const double expected = std::pow(16.0 / 11.0, 3.0);  // ~3.07

  EXPECT_GT(ratio, 0.7 * expected)
      << "SALT/GRANITE peak ratio " << ratio
      << " more than 30% below the analytic (16/11)^3 = " << expected;
  EXPECT_LT(ratio, 1.3 * expected)
      << "SALT/GRANITE peak ratio " << ratio
      << " more than 30% above the analytic (16/11)^3 = " << expected;
}

// GenericFallbackOnUnknownString: an unrecognised medium_type string
// must not crash the simulator. parseMediumType warns once and
// substitutes GENERIC. The peak amplitude must therefore equal the
// GENERIC-explicit run within tight tolerance.
TEST_F(MediumPlumbingTest, GenericFallbackOnUnknownString)
{
  const std::string cfg_explicit = writeConfig("generic", "GENERIC");
  const std::string cfg_unknown  = writeConfig("unknown", "MARS_REGOLITH");

  const double peak_g = runAndPeak(cfg_explicit, "generic");
  const double peak_u = runAndPeak(cfg_unknown,  "unknown");

  if (rank_ != 0) return;

  ASSERT_GT(peak_g, 0.0) << "Explicit GENERIC peak must be positive";
  ASSERT_GT(peak_u, 0.0)
      << "Unknown-medium peak must be positive (fallback to GENERIC)";

  // 0.1% tolerance: the only source of drift is the rank-0 stderr
  // warning print, which has no effect on the residual. Float-cmp on
  // the BHZ peak should agree to many digits.
  EXPECT_NEAR(peak_u, peak_g, 0.001 * peak_g)
      << "Unknown-medium peak " << peak_u
      << " differs from GENERIC-explicit " << peak_g
      << " by more than 0.1%";
}
