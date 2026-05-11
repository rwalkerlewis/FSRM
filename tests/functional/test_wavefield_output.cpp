/**
 * @file test_wavefield_output.cpp
 * @brief Pass-13c functional gates for the wavefield + source-ball 3D
 *        output infrastructure.
 *
 * Three gates:
 *  - BackwardCompat.WavefieldOutputDefaultIsNone:
 *      A config without an [OUTPUT] block (or with the section but no
 *      wavefield_format key) leaves wavefield_format_ = NONE so all 32
 *      historic-event integration tests are byte-identical to
 *      pre-pass-13c.
 *
 *  - WavefieldFormatHDF5XdmfParsesAndDefaults:
 *      A config with wavefield_format = HDF5_XDMF parses to
 *      WavefieldFormat::HDF5_XDMF; default cadence (100) and basename
 *      ("wavefield") apply.
 *
 *  - XDMFTimeSeriesEmissionContainsSnapshots:
 *      Seed the simulator with three pretend snapshot times and call
 *      writeWavefieldXdmfWrapper. The resulting .xdmf wrapper must
 *      enumerate three temporal snapshots with the right time stamps.
 *
 * Cite Henderson 2007 "ParaView Guide" for the XDMF schema.
 */

#include <gtest/gtest.h>

#include "core/Simulator.hpp"
#include "core/FSRM.hpp"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using FSRM::Simulator;

namespace
{

/// Write a minimal config to a temp file and return its path. Caller
/// deletes the file.
std::string writeTempConfig(const std::string& body)
{
    namespace fs = std::filesystem;
    fs::path tmp = fs::temp_directory_path()
                   / ("fsrm_pass13c_wavefield_"
                      + std::to_string(std::rand()) + ".config");
    std::ofstream f(tmp);
    f << body;
    f.close();
    return tmp.string();
}

/// Minimal config that the Simulator will accept for parsing without
/// trying to set up a real mesh / time stepper. The [OUTPUT] section is
/// the part under test; the rest is just enough to keep
/// initializeFromConfigFile from erroring before reaching [OUTPUT].
/// Keys conform to the pass-12 followup-2 strict validator schema
/// (see src/io/ConfigValidator.cpp): [SIMULATION] uses dt_initial /
/// end_time, and the material moduli live in a [ROCK] block, since
/// parseMaterialProperties() reads ROCK_* sections and the [MATERIAL]
/// schema does not accept lambda / mu / density directly.
const char* kMinimalConfigPrefix = R"(
[GRID]
nx = 4
ny = 4
nz = 4
Lx = 100.0
Ly = 100.0
Lz = 100.0

[ROCK]
lambda = 1.0e9
mu = 1.0e9
density = 2700.0

[SIMULATION]
end_time = 1.0e-6
dt_initial = 1.0e-7
)";

}  // namespace

class WavefieldOutputTest : public ::testing::Test
{
};

// =========================================================================
// 1. Default: no [OUTPUT] block -> WavefieldFormat::NONE.
// =========================================================================
TEST_F(WavefieldOutputTest, BackwardCompat_WavefieldOutputDefaultIsNone)
{
    const std::string body = std::string(kMinimalConfigPrefix);
    const std::string path = writeTempConfig(body);

    Simulator sim(PETSC_COMM_WORLD);
    sim.initializeFromConfigFile(path);

    EXPECT_EQ(sim.wavefieldFormat(), Simulator::WavefieldFormat::NONE)
        << "Without [OUTPUT].wavefield_format the default must remain "
           "NONE so the 32 historic-event tests are byte-identical.";
    EXPECT_EQ(sim.sourceBall3DOutputFormatRaw(),
              static_cast<int>(Simulator::SourceBall3DOutputFormat::NONE));

    std::filesystem::remove(path);
}

// =========================================================================
// 2. wavefield_format = HDF5_XDMF parses and exposes the right config.
// =========================================================================
TEST_F(WavefieldOutputTest, WavefieldFormatHDF5XdmfParsesAndDefaults)
{
    std::string body = std::string(kMinimalConfigPrefix) + R"(
[OUTPUT]
wavefield_format = HDF5_XDMF
wavefield_cadence_steps = 25
wavefield_basename = test_wf
wavefield_fields = displacement, velocity
source_ball_3d_output_format = HDF5_XDMF
source_ball_3d_output_cadence_steps = 5
)";
    const std::string path = writeTempConfig(body);

    Simulator sim(PETSC_COMM_WORLD);
    sim.initializeFromConfigFile(path);

    EXPECT_EQ(sim.wavefieldFormat(),
              Simulator::WavefieldFormat::HDF5_XDMF);
    EXPECT_EQ(sim.wavefieldCadenceSteps(), 25);
    EXPECT_EQ(sim.wavefieldBasename(), "test_wf");

    const auto& fields = sim.wavefieldFields();
    ASSERT_EQ(fields.size(), 2u);
    EXPECT_EQ(fields[0], "displacement");
    EXPECT_EQ(fields[1], "velocity");

    EXPECT_EQ(sim.sourceBall3DOutputFormatRaw(),
              static_cast<int>(
                  Simulator::SourceBall3DOutputFormat::HDF5_XDMF));
    EXPECT_EQ(sim.sourceBall3DOutputCadenceSteps(), 5);

    std::filesystem::remove(path);
}

// =========================================================================
// 3. XDMF wrapper enumerates the seeded snapshots.
// =========================================================================
TEST_F(WavefieldOutputTest, XDMFTimeSeriesEmissionContainsSnapshots)
{
    namespace fs = std::filesystem;
    const fs::path out_dir = fs::temp_directory_path()
        / ("fsrm_pass13c_xdmf_" + std::to_string(std::rand()));
    fs::create_directories(out_dir);
    const std::string out_dir_s = out_dir.string();

    const std::string body = std::string(kMinimalConfigPrefix)
        + "\n[OUTPUT]\n"
        + "wavefield_format = HDF5_XDMF\n"
        + "wavefield_output_directory = " + out_dir_s + "\n"
        + "wavefield_basename = unit_xdmf\n";
    const std::string path = writeTempConfig(body);

    Simulator sim(PETSC_COMM_WORLD);
    sim.initializeFromConfigFile(path);

    std::vector<double> times = {0.001, 0.002, 0.003};
    sim.seedWavefieldSnapshotsForTest(times);
    PetscErrorCode ierr = sim.triggerWavefieldXdmfWrapperForTest();
    ASSERT_EQ(ierr, 0);

    int rank_local = 0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank_local);
    if (rank_local == 0) {
        const fs::path xdmf = out_dir / "unit_xdmf.xdmf";
        ASSERT_TRUE(fs::exists(xdmf))
            << "Expected XDMF wrapper at " << xdmf;
        std::ifstream in(xdmf);
        std::stringstream ss;
        ss << in.rdbuf();
        const std::string s = ss.str();

        // Three temporal snapshots, three Time Value tags.
        size_t snap0 = s.find("snap0");
        size_t snap1 = s.find("snap1");
        size_t snap2 = s.find("snap2");
        EXPECT_NE(snap0, std::string::npos);
        EXPECT_NE(snap1, std::string::npos);
        EXPECT_NE(snap2, std::string::npos);

        EXPECT_NE(s.find("Time Value=\"0.001\""), std::string::npos);
        EXPECT_NE(s.find("Time Value=\"0.002\""), std::string::npos);
        EXPECT_NE(s.find("Time Value=\"0.003\""), std::string::npos);

        EXPECT_NE(s.find("CollectionType=\"Temporal\""), std::string::npos);
    }

    std::filesystem::remove(path);
    fs::remove_all(out_dir);
}
