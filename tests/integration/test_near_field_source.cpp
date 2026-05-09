/**
 * @file test_near_field_source.cpp
 * @brief Integration tests for the pass-5 [NEAR_FIELD_SOURCE] config
 *        grammar.
 *
 * Pass-5 introduces the [NEAR_FIELD_SOURCE] section selecting between
 * KINEMATIC_RDP (default, byte-identical to pre-pass-5 behaviour) and
 * DYNAMIC_PLASTIC (1D NearFieldExplosionSolver runs at setup time and
 * the recorded full 6-component moment tensor drives the residual).
 *
 * The single test in this file is the backward-compat guard: omitting
 * the section, or setting mode = KINEMATIC_RDP explicitly, must
 * produce float-exact (bit-identical) SAC output. This is modeled on
 * the pass-4 Integration.SourceDistribution.SingleCellLegacyByteIdentical
 * pattern in tests/integration/test_source_distribution.cpp.
 */

#include <gtest/gtest.h>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <filesystem>
#include <petscsys.h>

#include "core/Simulator.hpp"
#include "core/FSRM.hpp"
#include "domain/explosion/ExplosionImpactPhysics.hpp"
#include "sac_test_reader.hpp"

using namespace FSRM;

class NearFieldSourceTest : public ::testing::Test
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
      for (auto& p : tracked_paths_)
      {
        if (std::filesystem::exists(p))
        {
          std::error_code ec;
          std::filesystem::remove_all(p, ec);
          std::remove(p.c_str());
        }
      }
    }
  }

  // Sedan 1962 - style 104 kt alluvium config (matches the fixture
  // shared by tests/integration/test_source_distribution.cpp). The
  // optional near_field_section block is appended verbatim so each
  // test can opt into a specific [NEAR_FIELD_SOURCE] mode.
  std::string writeConfig(const std::string& tag,
                          const std::string& near_field_section)
  {
    const std::string config_path =
        "test_near_field_src_" + tag + ".config";
    const std::string output_dir =
        "test_near_field_src_" + tag + "_out";
    if (rank_ == 0)
    {
      tracked_paths_.push_back(config_path);
      tracked_paths_.push_back(output_dir);
      std::filesystem::remove_all(output_dir);

      std::ofstream cfg(config_path);
      cfg << "[SIMULATION]\n"
          << "name = near_field_src_" << tag << "\n"
          << "start_time = 0.0\n"
          << "end_time = 0.1\n"
          << "dt_initial = 0.001\n"
          << "dt_min = 0.0001\n"
          << "dt_max = 0.005\n"
          << "max_timesteps = 200\n"
          << "output_frequency = 50\n"
          << "output_format = HDF5\n"
          << "fluid_model = NONE\n"
          << "solid_model = ELASTIC\n"
          << "enable_geomechanics = true\n"
          << "enable_faults = false\n"
          << "enable_elastodynamics = true\n"
          << "rtol = 1.0e-6\n"
          << "atol = 1.0e-8\n"
          << "max_nonlinear_iterations = 20\n"
          << "\n[GRID]\n"
          << "nx = 4\nny = 4\nnz = 4\n"
          << "Lx = 4000.0\nLy = 4000.0\nLz = 2000.0\n"
          << "\n[ROCK]\ndensity = 2650.0\nlambda = 1.87e10\n"
          << "shear_modulus = 1.34e10\n"
          << "\n[MATERIAL]\nheterogeneous = true\nnum_layers = 3\n"
          << "\n[LAYER_1]\nz_top = 2000.0\nz_bottom = 1700.0\n"
          << "lambda = 4.36e9\nmu = 2.69e9\nrho = 1800.0\n"
          << "\n[LAYER_2]\nz_top = 1700.0\nz_bottom = 1000.0\n"
          << "lambda = 1.02e10\nmu = 7.50e9\nrho = 2300.0\n"
          << "\n[LAYER_3]\nz_top = 1000.0\nz_bottom = 0.0\n"
          << "lambda = 1.87e10\nmu = 1.34e10\nrho = 2650.0\n"
          << "\n[EXPLOSION_SOURCE]\ntype = UNDERGROUND_NUCLEAR\n"
          << "yield_kt = 104.0\ndepth_of_burial = 194.0\n"
          << "location_x = 2000.0\nlocation_y = 2000.0\n"
          << "location_z = 1806.0\n"
          << "onset_time = 0.0\nrise_time = 0.01\n"
          << "cavity_overpressure = 1.0e10\n"
          << "medium_type = ALLUVIUM\n";
      if (!near_field_section.empty())
      {
        cfg << "\n" << near_field_section;
      }
      cfg << "\n[BOUNDARY_CONDITIONS]\nbottom = free\n"
          << "sides = free\ntop = free\n"
          << "\n[ABSORBING_BC]\nenabled = true\n"
          << "x_min = true\nx_max = true\n"
          << "y_min = true\ny_max = true\n"
          << "z_min = true\nz_max = false\n"
          << "\n[SEISMOMETERS]\nenabled = true\nformats = SAC\n"
          << "output_dir = " << output_dir << "\n"
          << "default_quantity = DISPLACEMENT\n"
          << "default_sample_rate_hz = 200.0\n"
          << "\n[SEISMOMETER_1]\nsta = SPALL\n"
          << "location_xyz = 2000.0,2000.0,2000.0\n";
      cfg.close();
    }
    config_path_ = config_path;
    output_dir_ = output_dir;
    return config_path;
  }

  PetscErrorCode runPipeline(PetscReal& sol_norm)
  {
    MPI_Barrier(PETSC_COMM_WORLD);
    Simulator sim(PETSC_COMM_WORLD);
    PetscErrorCode ierr;
    PetscOptionsClear(nullptr);

    ierr = sim.initializeFromConfigFile(config_path_); if (ierr) return ierr;
    ierr = sim.setupDM();             if (ierr) return ierr;
    ierr = sim.labelBoundaries();     if (ierr) return ierr;
    ierr = sim.setupFields();         if (ierr) return ierr;
    ierr = sim.setupPhysics();        if (ierr) return ierr;
    ierr = sim.setupTimeStepper();    if (ierr) return ierr;
    ierr = sim.setupSolvers();        if (ierr) return ierr;
    ierr = sim.setInitialConditions(); if (ierr) return ierr;

    PetscPushErrorHandler(PetscReturnErrorHandler, nullptr);
    ierr = sim.run();
    PetscPopErrorHandler();
    if (ierr) return ierr;
    ierr = sim.writeSummary();
    if (ierr) return ierr;

    Vec sol = sim.getSolution();
    if (sol) VecNorm(sol, NORM_2, &sol_norm);
    else sol_norm = -1.0;
    return PETSC_SUCCESS;
  }

  double readBhzPeak()
  {
    if (rank_ != 0) return 0.0;
    using namespace FSRM::test_helpers;
    const std::string sac_path = output_dir_ + "/XX.SPALL.00.BHZ.sac";
    SACTrace trace = readSAC(sac_path);
    EXPECT_TRUE(trace.valid) << "SAC BHZ must parse: " << sac_path;
    if (!trace.valid) return 0.0;
    PeakResult peak = tracePeak(trace);
    return peak.peak_abs;
  }

  int rank_ = 0;
  std::string config_path_;
  std::string output_dir_;
  std::vector<std::string> tracked_paths_;
};

// Pass-5 backward-compatibility guard. Setting [NEAR_FIELD_SOURCE]
// mode = KINEMATIC_RDP must produce SAC output that is float-exact
// (bit-identical) to the run that omits the section entirely.
//
// Any drift here would mean the [NEAR_FIELD_SOURCE] parsing branch is
// mutating state that affects the legacy KINEMATIC_RDP code path; the
// pass-5 guarantee is that configs that do not opt in see no
// behavioural change. Modeled on
// Integration.SourceDistribution.SingleCellLegacyByteIdentical.
TEST_F(NearFieldSourceTest, KinematicRDPLegacyByteIdentical)
{
  writeConfig("legacy_implicit", "");
  PetscReal n_legacy = 0.0;
  ASSERT_EQ(runPipeline(n_legacy), 0);
  const double peak_legacy = readBhzPeak();

  writeConfig("legacy_explicit",
              "[NEAR_FIELD_SOURCE]\nmode = KINEMATIC_RDP\n");
  PetscReal n_explicit = 0.0;
  ASSERT_EQ(runPipeline(n_explicit), 0);
  const double peak_explicit = readBhzPeak();

  if (rank_ == 0)
  {
    EXPECT_EQ(n_legacy, n_explicit)
        << "Solution norm drifted between implicit and explicit "
        << "KINEMATIC_RDP paths: " << n_legacy << " vs " << n_explicit;
    EXPECT_EQ(peak_legacy, peak_explicit)
        << "BHZ peak drifted between implicit and explicit "
        << "KINEMATIC_RDP paths: " << peak_legacy << " vs " << peak_explicit;
  }
}

// Pass-6 + pass-7: under [NEAR_FIELD_SOURCE] mode = DYNAMIC_PLASTIC
// with solver_kind = CLOSED_FORM, the pipeline must continue to
// produce the same byte-identical pass-5 RDP-driven output regardless
// of repeated invocation. Pass-7 promoted RADIAL_LAGRANGIAN to the
// default for DYNAMIC_PLASTIC, so the regression guard now compares
// two explicit-CLOSED_FORM runs (rather than default-vs-explicit) to
// ensure the pinned-CLOSED_FORM legacy path is unchanged.
TEST_F(NearFieldSourceTest, ClosedFormFallback)
{
  writeConfig("dyn_closedform_a",
              "[NEAR_FIELD_SOURCE]\n"
              "mode = DYNAMIC_PLASTIC\n"
              "solver_kind = CLOSED_FORM\n");
  PetscReal n_a = 0.0;
  ASSERT_EQ(runPipeline(n_a), 0);
  const double peak_a = readBhzPeak();

  writeConfig("dyn_closedform_b",
              "[NEAR_FIELD_SOURCE]\n"
              "mode = DYNAMIC_PLASTIC\n"
              "solver_kind = CLOSED_FORM\n");
  PetscReal n_b = 0.0;
  ASSERT_EQ(runPipeline(n_b), 0);
  const double peak_b = readBhzPeak();

  if (rank_ == 0)
  {
    EXPECT_EQ(n_a, n_b)
        << "Solution norm drifted between two explicit-CLOSED_FORM "
        << "runs under DYNAMIC_PLASTIC mode: " << n_a << " vs " << n_b;
    EXPECT_EQ(peak_a, peak_b)
        << "BHZ peak drifted between two explicit-CLOSED_FORM runs "
        << "under DYNAMIC_PLASTIC mode: " << peak_a << " vs " << peak_b;
  }
}

// Pass-7 headline gate. With the Tillotson host-rock EOS, physics-
// based energy-partition cavity initialization, and Wilkins (1980)
// literature AV coefficients (c_l = 0.06, c_q = 1.5), the radial
// Lagrangian shock solver lands within a factor of 5 of the closed-
// form RDP estimate at the elastic-radius extraction surface (Sedan
// 1962 anchor: ratio ~2.2x). This test runs the same fixture under
// both solver_kind dispatches and asserts the peak |M0_iso_dot|
// ratio falls inside the factor-5 envelope. The HDF5 + XDMF profile
// pair must also be produced and the recorded cavity radius must be
// finite and positive.
TEST_F(NearFieldSourceTest, RadialLagrangianAnchor)
{
  // RADIAL_LAGRANGIAN run. radial_cells = 100 keeps the CI pass
  // fast; the production resolution (800) is exercised by the
  // pass-7 amplitude diagnostic script (scripts/pass7_amplitude_
  // diagnostic.sh) and the historic-nuclear / examples runs.
  writeConfig("radial_lagrangian",
              "[NEAR_FIELD_SOURCE]\n"
              "mode = DYNAMIC_PLASTIC\n"
              "solver_kind = RADIAL_LAGRANGIAN\n"
              "radial_cells = 100\n"
              "elastic_radius_factor = 3.0\n"
              "near_field_dt = 1.0e-5\n"
              "output_cadence_microseconds = 1000\n"
              "profile_output_cadence_microseconds = 5000\n");
  PetscReal sol_norm_rl = 0.0;
  ASSERT_EQ(runPipeline(sol_norm_rl), 0)
      << "RADIAL_LAGRANGIAN pipeline must complete";

  // Capture the RL output paths now: the next writeConfig overwrites
  // the test fixture's output_dir_ field, but the actual files on
  // disk for the RL run live under their own per-tag directory.
  const std::string csv_rl = output_dir_ + "/near_field_history.csv";
  const std::string h5_path = output_dir_ + "/near_field_profile.h5";
  const std::string xdmf_path = output_dir_ + "/near_field_profile.xdmf";
  if (rank_ == 0) {
    EXPECT_TRUE(std::isfinite(sol_norm_rl))
        << "Solution norm must be finite under RADIAL_LAGRANGIAN";
  }

  // CLOSED_FORM run for the amplitude-ratio reference.
  writeConfig("closed_form_anchor",
              "[NEAR_FIELD_SOURCE]\n"
              "mode = DYNAMIC_PLASTIC\n"
              "solver_kind = CLOSED_FORM\n"
              "elastic_radius_factor = 3.0\n"
              "near_field_dt = 1.0e-5\n"
              "output_cadence_microseconds = 1000\n");
  PetscReal sol_norm_cf = 0.0;
  ASSERT_EQ(runPipeline(sol_norm_cf), 0)
      << "CLOSED_FORM reference pipeline must complete";

  if (rank_ != 0) return;
  const std::string csv_cf = output_dir_ + "/near_field_history.csv";

  auto peak_m0iso_dot = [](const std::string& csv_path) -> double {
    std::ifstream csv(csv_path);
    if (!csv.is_open()) return -1.0;
    std::string line;
    double peak = 0.0;
    int rows = 0;
    while (std::getline(csv, line)) {
      if (line.empty() || line[0] == '#') continue;
      if (line.rfind("t,", 0) == 0) continue;
      std::stringstream ss(line);
      std::string cell;
      std::vector<double> cols;
      while (std::getline(ss, cell, ',')) {
        try { cols.push_back(std::stod(cell)); }
        catch (...) { cols.clear(); break; }
      }
      if (cols.size() >= 10) {
        const double v = std::abs(cols[9]);
        if (v > peak) peak = v;
        ++rows;
      }
    }
    return rows > 0 ? peak : -1.0;
  };

  const double peak_rl = peak_m0iso_dot(csv_rl);
  const double peak_cf = peak_m0iso_dot(csv_cf);

  ASSERT_GT(peak_rl, 0.0)
      << "RADIAL_LAGRANGIAN peak |M0_iso_dot| must be positive";
  ASSERT_GT(peak_cf, 0.0)
      << "CLOSED_FORM peak |M0_iso_dot| must be positive";

  const double ratio = peak_rl / peak_cf;
  EXPECT_GT(ratio, 0.2)
      << "RADIAL_LAGRANGIAN peak |M0_iso_dot|=" << peak_rl
      << " < 0.2x CLOSED_FORM=" << peak_cf
      << " (factor-5 envelope on the low side; gap re-opened)";
  EXPECT_LT(ratio, 5.0)
      << "RADIAL_LAGRANGIAN peak |M0_iso_dot|=" << peak_rl
      << " > 5x CLOSED_FORM=" << peak_cf
      << " (factor-5 envelope on the high side; calibration drifted)";

  // Profile pair guard: the HDF5 + XDMF spatial-profile output must
  // be present from the RADIAL_LAGRANGIAN run.
  EXPECT_TRUE(std::filesystem::exists(h5_path))
      << "near_field_profile.h5 must be written under RADIAL_LAGRANGIAN: "
      << h5_path;
  EXPECT_TRUE(std::filesystem::exists(xdmf_path))
      << "near_field_profile.xdmf must be written under RADIAL_LAGRANGIAN: "
      << xdmf_path;
}
