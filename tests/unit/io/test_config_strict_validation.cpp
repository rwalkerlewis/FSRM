// Pass-12 followup 2 (V&V hardening): strict configuration
// validation. Each test asserts on a single rejection class so a
// regression points at exactly one rule, and the error messages
// are inspected for the named field per the test-suite convention.

#include <gtest/gtest.h>

#include "core/ConfigReader.hpp"
#include "io/ConfigValidator.hpp"

#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <mpi.h>
#include <string>

using FSRM::ConfigReader;
using FSRM::ConfigValidator;

namespace {

// Concatenate a vector of strings for substring inspection.
std::string concat(const std::vector<std::string>& v) {
    std::string s;
    for (const auto& e : v) { s += e; s += "\n"; }
    return s;
}

class ConfigStrictValidationTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
        // Make sure the env-var bypass is not set; rank-0 alone
        // sees this filesystem so the unset is global enough.
        unsetenv("FSRM_DISABLE_STRICT_VALIDATION");
        cfg_path_ = "test_config_strict_validation_" +
                    std::to_string(rank_) + ".cfg";
    }

    void TearDown() override {
        std::error_code ec;
        std::filesystem::remove(cfg_path_, ec);
    }

    void writeConfig(const std::string& body) {
        std::ofstream f(cfg_path_);
        f << body;
    }

    int rank_ = 0;
    std::string cfg_path_;
};

const char* kMinimalValidExplosion = R"CFG(
[SIMULATION]
name = unit_test_explosion
start_time = 0.0
end_time = 0.001
dt_initial = 0.0001
max_timesteps = 10
output_format = HDF5
fluid_model = NONE
solid_model = ELASTIC
enable_geomechanics = true
enable_elastodynamics = true

[GRID]
nx = 4
ny = 4
nz = 4
Lx = 1000.0
Ly = 1000.0
Lz = 500.0

[ROCK]
density = 2650.0
youngs_modulus = 30.0e9
poisson_ratio = 0.25

[EXPLOSION_SOURCE]
type = UNDERGROUND_NUCLEAR
yield_kt = 1.0
depth_of_burial = 250.0
location_x = 500.0
location_y = 500.0
location_z = 250.0
onset_time = 0.0
rise_time = 0.005
cavity_overpressure = 1.0e10

[BOUNDARY_CONDITIONS]
bottom = free
sides = free
top = free
)CFG";

}  // namespace

// ---------------------------------------------------------------------------
// Sanity: the canonical fixture passes strict validation.
// ---------------------------------------------------------------------------
TEST_F(ConfigStrictValidationTest, MinimalExplosionConfigValidates) {
    writeConfig(kMinimalValidExplosion);
    ConfigReader reader;
    ASSERT_TRUE(reader.loadFile(cfg_path_));
    auto r = ConfigValidator::validate(reader, cfg_path_);
    EXPECT_TRUE(r.valid) << "errors=" << concat(r.errors);
    EXPECT_TRUE(r.errors.empty());
}

// ---------------------------------------------------------------------------
// Rejection class 1: unknown top-level section.
// ---------------------------------------------------------------------------
TEST_F(ConfigStrictValidationTest, UnknownSectionRejected) {
    writeConfig(std::string(kMinimalValidExplosion) +
                "\n[NOT_A_SECTION]\nfoo = bar\n");
    ConfigReader reader;
    ASSERT_TRUE(reader.loadFile(cfg_path_));
    auto r = ConfigValidator::validate(reader, cfg_path_);
    EXPECT_FALSE(r.valid);
    const std::string all = concat(r.errors);
    EXPECT_NE(all.find("NOT_A_SECTION"), std::string::npos)
        << "Error must name the offending section: " << all;
}

// ---------------------------------------------------------------------------
// Rejection class 2: unknown key inside a known section.
// ---------------------------------------------------------------------------
TEST_F(ConfigStrictValidationTest, UnknownKeyInKnownSectionRejected) {
    writeConfig(std::string(kMinimalValidExplosion) +
                "\n[SIMULATION]\nbogus_key_does_not_exist = 1\n");
    ConfigReader reader;
    ASSERT_TRUE(reader.loadFile(cfg_path_));
    auto r = ConfigValidator::validate(reader, cfg_path_);
    EXPECT_FALSE(r.valid);
    const std::string all = concat(r.errors);
    EXPECT_NE(all.find("bogus_key_does_not_exist"), std::string::npos)
        << "Error must name the offending key: " << all;
}

// ---------------------------------------------------------------------------
// Rejection class 3: deprecated SEISMOMETERS station_<n> form
// (the original PR #128 finding).
// ---------------------------------------------------------------------------
TEST_F(ConfigStrictValidationTest, DeprecatedStationKeyInSeismometersRejected) {
    writeConfig(std::string(kMinimalValidExplosion) + R"CFG(
[SEISMOMETERS]
enabled = true
formats = SAC
output_dir = output
station_1 = 1000.0, 1000.0, 0.0
station_1_name = STA01
)CFG");
    ConfigReader reader;
    ASSERT_TRUE(reader.loadFile(cfg_path_));
    auto r = ConfigValidator::validate(reader, cfg_path_);
    EXPECT_FALSE(r.valid);
    const std::string all = concat(r.errors);
    EXPECT_NE(all.find("station_1"), std::string::npos);
    EXPECT_NE(all.find("[SEISMOMETER_"), std::string::npos)
        << "Error must direct the user at the per-station block: "
        << all;
}

// ---------------------------------------------------------------------------
// Rejection class 3 (variant): orphan [EXPLOSION] section (the
// other PR #128 finding -- block name without _SOURCE suffix).
// ---------------------------------------------------------------------------
TEST_F(ConfigStrictValidationTest, DeprecatedOrphanExplosionSectionRejected) {
    // Build a minimal config that swaps EXPLOSION_SOURCE for
    // EXPLOSION (which the simulator parser silently ignores). The
    // schema rejects it via the EXPLOSION-* deprecated form.
    const std::string body = R"CFG(
[SIMULATION]
name = orphan_explosion
start_time = 0.0
end_time = 0.001
fluid_model = NONE
solid_model = ELASTIC
enable_geomechanics = true

[GRID]
nx = 4
ny = 4
nz = 4
Lx = 100.0
Ly = 100.0
Lz = 100.0

[ROCK]
density = 2650.0
youngs_modulus = 30.0e9
poisson_ratio = 0.25

[EXPLOSION]
enabled = true
yield_kt = 1.0
depth = 250.0

[BOUNDARY_CONDITIONS]
bottom = free
sides = free
top = free
)CFG";
    writeConfig(body);
    ConfigReader reader;
    ASSERT_TRUE(reader.loadFile(cfg_path_));
    auto r = ConfigValidator::validate(reader, cfg_path_);
    EXPECT_FALSE(r.valid);
    const std::string all = concat(r.errors);
    EXPECT_NE(all.find("EXPLOSION_SOURCE"), std::string::npos)
        << "Error must point at the canonical replacement: " << all;
}

// ---------------------------------------------------------------------------
// Rejection class 4: required-section omission.
// EXPLOSION_SOURCE present, BOUNDARY_CONDITIONS missing.
// ---------------------------------------------------------------------------
TEST_F(ConfigStrictValidationTest, ExplosionWithoutBoundaryConditionsRejected) {
    // Minimal explosion config WITHOUT the BC block. Mirrors the
    // PR #127 bug where a missing [BOUNDARY_CONDITIONS] silently
    // suppressed the explosion source coupling.
    const std::string body = R"CFG(
[SIMULATION]
name = no_bc_explosion
start_time = 0.0
end_time = 0.001
fluid_model = NONE
solid_model = ELASTIC
enable_geomechanics = true
enable_elastodynamics = true

[GRID]
nx = 4
ny = 4
nz = 4
Lx = 1000.0
Ly = 1000.0
Lz = 500.0

[ROCK]
density = 2650.0
youngs_modulus = 30.0e9
poisson_ratio = 0.25

[EXPLOSION_SOURCE]
type = UNDERGROUND_NUCLEAR
yield_kt = 1.0
depth_of_burial = 250.0
location_x = 500.0
location_y = 500.0
location_z = 250.0
onset_time = 0.0
rise_time = 0.005
cavity_overpressure = 1.0e10
)CFG";
    writeConfig(body);
    ConfigReader reader;
    ASSERT_TRUE(reader.loadFile(cfg_path_));
    auto r = ConfigValidator::validate(reader, cfg_path_);
    EXPECT_FALSE(r.valid);
    const std::string all = concat(r.errors);
    EXPECT_NE(all.find("BOUNDARY_CONDITIONS"), std::string::npos)
        << "Error must name the missing BC block: " << all;
    EXPECT_NE(all.find("EXPLOSION_SOURCE"), std::string::npos)
        << "Error must name the section that triggered the rule: "
        << all;
}

// ---------------------------------------------------------------------------
// Opt-out: [META] strict_validation = false bypasses validation
// with a warning, and Result::valid stays true.
// ---------------------------------------------------------------------------
TEST_F(ConfigStrictValidationTest, MetaOptOutEmitsWarningAndPasses) {
    writeConfig(R"CFG(
[META]
strict_validation = false

[SIMULATION]
name = opt_out_test
end_time = 1.0

[NOT_A_SECTION]
this_should_be_invalid = true
)CFG");
    ConfigReader reader;
    ASSERT_TRUE(reader.loadFile(cfg_path_));
    auto r = ConfigValidator::validate(reader, cfg_path_);
    EXPECT_TRUE(r.valid);
    EXPECT_TRUE(r.errors.empty());
    EXPECT_FALSE(r.warnings.empty())
        << "Opt-out must emit a stderr warning so it never goes "
           "silent.";
    const std::string warns = concat(r.warnings);
    EXPECT_NE(warns.find("strict_validation"), std::string::npos);
}

// ---------------------------------------------------------------------------
// Env-var bypass: FSRM_DISABLE_STRICT_VALIDATION=1 is the
// emergency override used by legacy callers and the examples-runtime
// CTest gate when the underlying example pre-dates strict mode.
// ---------------------------------------------------------------------------
TEST_F(ConfigStrictValidationTest, EnvVarBypassSkipsValidation) {
    setenv("FSRM_DISABLE_STRICT_VALIDATION", "1", 1);
    writeConfig(std::string(kMinimalValidExplosion) +
                "\n[NOT_A_SECTION]\nfoo = bar\n");
    ConfigReader reader;
    ASSERT_TRUE(reader.loadFile(cfg_path_));
    auto r = ConfigValidator::validate(reader, cfg_path_);
    EXPECT_TRUE(r.valid);
    EXPECT_TRUE(r.errors.empty());
    unsetenv("FSRM_DISABLE_STRICT_VALIDATION");
}

// ---------------------------------------------------------------------------
// Dynamic-prefix sections: LAYER_<n>, SEISMOMETER_<n>, etc. are
// accepted with a numeric suffix and validated against the prefix
// schema.
// ---------------------------------------------------------------------------
TEST_F(ConfigStrictValidationTest, DynamicPrefixSectionsAccepted) {
    writeConfig(std::string(kMinimalValidExplosion) + R"CFG(
[LAYER_1]
z_top = 500.0
z_bottom = 250.0
lambda = 1.0e10
mu = 5.0e9
rho = 2650.0

[LAYER_42]
z_top = 250.0
z_bottom = 0.0
lambda = 2.0e10
mu = 1.0e10
rho = 2700.0

[SEISMOMETERS]
enabled = true
formats = SAC
output_dir = out

[SEISMOMETER_3]
sta = STA03
location_xyz = 100.0,100.0,500.0
)CFG");
    ConfigReader reader;
    ASSERT_TRUE(reader.loadFile(cfg_path_));
    auto r = ConfigValidator::validate(reader, cfg_path_);
    EXPECT_TRUE(r.valid) << "errors=" << concat(r.errors);
}

// ---------------------------------------------------------------------------
// Coverage: walk every config in examples/ and assert it validates
// under strict mode. This is the "every example you can run today
// passes the schema" claim.
// ---------------------------------------------------------------------------
TEST_F(ConfigStrictValidationTest, AllShippedExampleConfigsValidate) {
    if (rank_ != 0) GTEST_SKIP() << "rank-0 only";
    namespace fs = std::filesystem;
    // CMake injects FSRM_REPO_ROOT_DIR via target_compile_definitions
    // so this test is independent of the ctest working directory.
#ifndef FSRM_REPO_ROOT_DIR
#error "FSRM_REPO_ROOT_DIR must be defined by CMake (see tests/CMakeLists.txt)."
#endif
    const fs::path base = fs::path(FSRM_REPO_ROOT_DIR) / "examples";
    ASSERT_TRUE(fs::exists(base) && fs::is_directory(base))
        << "examples/ not found under FSRM_REPO_ROOT_DIR="
        << FSRM_REPO_ROOT_DIR;

    int n_validated = 0;
    int n_failed = 0;
    std::vector<std::string> failure_messages;
    for (const auto& entry : fs::recursive_directory_iterator(base)) {
        if (!entry.is_regular_file()) continue;
        const std::string p = entry.path().string();
        if (p.size() < 8) continue;
        if (p.substr(p.size() - 7) != ".config") continue;

        ConfigReader reader;
        ASSERT_TRUE(reader.loadFile(p)) << "loadFile failed: " << p;
        auto r = ConfigValidator::validate(reader, p);
        ++n_validated;
        if (!r.valid) {
            ++n_failed;
            failure_messages.push_back("FAILED: " + p + "\n" +
                                       concat(r.errors));
        }
    }
    EXPECT_GT(n_validated, 0)
        << "Expected at least one example config under " << base;
    EXPECT_EQ(n_failed, 0)
        << "Some example configs do not validate under strict mode "
           "(repair the configs or extend the schema in "
           "src/io/ConfigValidator.cpp):\n"
        << concat(failure_messages);
}
