/**
 * @file test_sedan_1962_run.cpp
 * @brief Integration test: run the Sedan 1962 validated configuration end-to-end.
 *
 * This is the M1.3 (run + receivers) test produced by Session 2 of the
 * nuclear-explosion-monitoring (NEM) roadmap.  It is deliberately a
 * "runs end-to-end" test, not a validation test.  PPV-vs-range
 * comparison against published Sedan data lives in M1.4 (Session 3).
 *
 * Required behavior:
 *   1. The `config/examples/sedan_1962_validated.config` configuration
 *      drives Simulator through the full pipeline without error.
 *   2. SAC files are produced for all 8 stations, 3 components each,
 *      with at least 100 non-zero samples per file.  At 200 Hz that is
 *      0.5 s, comfortably longer than the slowest direct P arrival at
 *      the 2 km surface station (~570 ms via basement vp = 3500 m/s).
 *   3. On the free surface, peak |BHZ| amplitude monotonically
 *      decreases with radial distance (sanity check against a
 *      sign-flipped source or backwards source-time function).
 *   4. The first-arrival time for BHZ at SUR_R2000 is between 0.4 s
 *      and 0.8 s (theoretical ~570 ms; ±200 ms slack covers the layered
 *      ray path and the 5%-of-peak detection threshold).
 *
 * Pattern: clones the structure of test_explosion_seismogram.cpp and
 * test_punggye_ri_layered.cpp but uses the committed Gmsh mesh and
 * shipped configuration file under config/examples/.
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
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

bool readSACFile(const std::string& path, int& npts, float& delta,
                 std::vector<float>& data)
{
  std::ifstream in(path, std::ios::binary);
  if (!in.is_open()) return false;

  SacHeaderRead hdr;
  in.read(reinterpret_cast<char*>(&hdr), sizeof(hdr));
  if (!in.good()) return false;

  delta = hdr.f[0];   // DELTA field
  npts = hdr.i[9];    // NPTS field
  if (npts <= 0 || npts >= 10000000) return false;

  data.resize(static_cast<size_t>(npts));
  in.read(reinterpret_cast<char*>(data.data()),
          static_cast<std::streamsize>(npts * sizeof(float)));
  return in.good() || in.eof();
}

float maxAbsAmplitude(const std::vector<float>& data)
{
  float mx = 0.0f;
  for (float v : data)
  {
    mx = std::max(mx, std::fabs(v));
  }
  return mx;
}

int countNonzeroSamples(const std::vector<float>& data)
{
  int count = 0;
  for (float v : data)
  {
    if (v != 0.0f) ++count;
  }
  return count;
}

// Run the full Simulator pipeline on the validated config file.
PetscErrorCode runSimulationOnce(const char* config_path)
{
  PetscFunctionBeginUser;
  Simulator sim(PETSC_COMM_WORLD);
  PetscErrorCode ierr;
  ierr = sim.initializeFromConfigFile(config_path); CHKERRQ(ierr);
  ierr = sim.setupDM();                              CHKERRQ(ierr);
  ierr = sim.labelBoundaries();                      CHKERRQ(ierr);
  ierr = sim.setupFields();                          CHKERRQ(ierr);
  ierr = sim.setupPhysics();                         CHKERRQ(ierr);
  ierr = sim.setupTimeStepper();                     CHKERRQ(ierr);
  ierr = sim.setupSolvers();                         CHKERRQ(ierr);
  ierr = sim.setInitialConditions();                 CHKERRQ(ierr);
  ierr = sim.run();                                  CHKERRQ(ierr);
  ierr = sim.writeSummary();                         CHKERRQ(ierr);
  PetscFunctionReturn(PETSC_SUCCESS);
}

// Resolve paths relative to the build directory used by ctest, matching
// the convention in test_sedan_1962_mesh.cpp and other Gmsh-based
// integration tests (run from <repo>/build).
constexpr const char* kConfigPath =
    "../config/examples/sedan_1962_validated.config";
constexpr const char* kSacDir =
    "output/sedan_1962_validated";

const std::vector<std::string> kSurfaceStations = {
    "SUR_R0100", "SUR_R0500", "SUR_R1000", "SUR_R2000"};
const std::vector<std::string> kDepthStations = {
    "DEP_R0100", "DEP_R0500", "DEP_R1000", "DEP_R2000"};
const std::vector<std::string> kComponents = {"BHN", "BHE", "BHZ"};

std::string sacPath(const std::string& sta, const std::string& comp)
{
  return std::string(kSacDir) + "/XX." + sta + ".00." + comp + ".sac";
}

}  // namespace

class Sedan1962RunTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank_);
    if (rank_ == 0)
    {
      // Clear any stale SAC output from a previous run so the
      // "non-zero samples" assertions are not satisfied by leftovers.
      std::filesystem::remove_all(kSacDir);
    }
    MPI_Barrier(PETSC_COMM_WORLD);
  }

  void TearDown() override
  {
    MPI_Barrier(PETSC_COMM_WORLD);
    if (rank_ == 0)
    {
      std::filesystem::remove_all(kSacDir);
    }
  }

  int rank_ = 0;
};

// 1. The full pipeline runs to completion without error.
TEST_F(Sedan1962RunTest, SimulationCompletes)
{
  ASSERT_TRUE(std::filesystem::exists(kConfigPath))
      << "Config file not found at " << kConfigPath
      << " (cwd should be the build directory)";

  PetscErrorCode ierr = runSimulationOnce(kConfigPath);
  ASSERT_EQ(ierr, PETSC_SUCCESS) << "Sedan 1962 validated run failed";
}

// 2. SAC files are produced for every station / component and contain
//    at least 100 non-zero samples (0.5 s of signal at 200 Hz).
TEST_F(Sedan1962RunTest, AllSACFilesProducedAndNonzero)
{
  ASSERT_EQ(runSimulationOnce(kConfigPath), PETSC_SUCCESS);
  if (rank_ != 0) return;

  std::vector<std::string> all_stations = kSurfaceStations;
  all_stations.insert(all_stations.end(),
                      kDepthStations.begin(), kDepthStations.end());

  for (const auto& sta : all_stations)
  {
    for (const auto& comp : kComponents)
    {
      const std::string p = sacPath(sta, comp);
      ASSERT_TRUE(std::filesystem::exists(p)) << "SAC file missing: " << p;

      int npts = 0;
      float delta = 0.0f;
      std::vector<float> data;
      ASSERT_TRUE(readSACFile(p, npts, delta, data))
          << "Failed to read SAC file: " << p;
      EXPECT_GT(npts, 100)
          << "Expected > 100 samples in " << p << ", got " << npts;
      EXPECT_GT(countNonzeroSamples(data), 100)
          << "Expected > 100 non-zero samples in " << p;
    }
  }
}

// 3. Surface BHZ peak amplitude monotonically decays with range.
//    Sanity check against a sign-flipped source-time function or a
//    swapped channel ordering.  Geometric-spreading slope check is M1.4.
TEST_F(Sedan1962RunTest, AmplitudeMonotonicallyDecaysOnSurface)
{
  ASSERT_EQ(runSimulationOnce(kConfigPath), PETSC_SUCCESS);
  if (rank_ != 0) return;

  std::array<float, 4> peaks{};
  for (size_t i = 0; i < kSurfaceStations.size(); ++i)
  {
    int npts = 0;
    float delta = 0.0f;
    std::vector<float> data;
    const std::string p = sacPath(kSurfaceStations[i], "BHZ");
    ASSERT_TRUE(readSACFile(p, npts, delta, data)) << p;
    peaks[i] = maxAbsAmplitude(data);
    EXPECT_GT(peaks[i], 0.0f) << "Zero peak at " << kSurfaceStations[i];
  }

  EXPECT_GT(peaks[0], peaks[1])
      << "BHZ peak at 100 m should exceed 500 m (got "
      << peaks[0] << " vs " << peaks[1] << ")";
  EXPECT_GT(peaks[1], peaks[2])
      << "BHZ peak at 500 m should exceed 1000 m (got "
      << peaks[1] << " vs " << peaks[2] << ")";
  EXPECT_GT(peaks[2], peaks[3])
      << "BHZ peak at 1000 m should exceed 2000 m (got "
      << peaks[2] << " vs " << peaks[3] << ")";
}

// 4. First-arrival time for BHZ at SUR_R2000 lies in [0.4, 0.8] s.
//    Direct-P theoretical arrival ~ 2000 / 3500 = 0.571 s; we allow
//    +-200 ms for the actual layered ray path and the 5%-of-peak
//    threshold detection.
TEST_F(Sedan1962RunTest, PFirstArrivalTimeReasonable)
{
  ASSERT_EQ(runSimulationOnce(kConfigPath), PETSC_SUCCESS);
  if (rank_ != 0) return;

  int npts = 0;
  float delta = 0.0f;
  std::vector<float> data;
  const std::string p = sacPath("SUR_R2000", "BHZ");
  ASSERT_TRUE(readSACFile(p, npts, delta, data)) << p;
  ASSERT_GT(npts, 1);
  ASSERT_GT(delta, 0.0f);

  const float peak = maxAbsAmplitude(data);
  ASSERT_GT(peak, 0.0f) << "SUR_R2000 BHZ peak is zero";

  const float threshold = 0.05f * peak;  // 5% of peak for first-arrival detection
  int first_idx = -1;
  for (int i = 0; i < npts; ++i)
  {
    if (std::fabs(data[i]) > threshold)
    {
      first_idx = i;
      break;
    }
  }
  ASSERT_GE(first_idx, 0) << "No sample above 5% threshold at SUR_R2000 BHZ";

  const float t_arrival = static_cast<float>(first_idx) * delta;
  EXPECT_GE(t_arrival, 0.4f)
      << "First-arrival time " << t_arrival
      << " s is suspiciously early (theoretical ~0.57 s)";
  EXPECT_LE(t_arrival, 0.8f)
      << "First-arrival time " << t_arrival
      << " s is suspiciously late (theoretical ~0.57 s)";
}
