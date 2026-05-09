/**
 * @file test_paraview_state_load.cpp
 * @brief Pass-7 functional test: verify the hand-authored ParaView
 *        state file loads under pvpython without errors.
 *
 * Closes the pass-6 verification debt for
 * examples/11_sedan_1962/paraview/near_field_cavity.pvsm. Pass-6
 * shipped the .pvsm with documented ParaView 5.10+ proxy types but
 * never verified the rendered output. Pass-7 invokes pvpython on the
 * file via scripts/verify_pvsm.sh and asserts it loads cleanly.
 *
 * The fsrm-ci:local Docker image does not include ParaView. Under
 * that environment the script returns exit code 2 (pvpython missing)
 * and the test records GTEST_SKIP with an explicit message. Users
 * who want the verification on their workstation must install
 * ParaView 5.10+ and run the script standalone:
 *
 *     scripts/verify_pvsm.sh
 *
 * In a development environment that has ParaView installed the test
 * gates the .pvsm against silent rendering bugs (proxy-type drift,
 * missing data files, schema mismatches) that the hand-authored XML
 * is sensitive to.
 */

#include <gtest/gtest.h>

#include <cstdlib>
#include <filesystem>
#include <string>

class ParaViewStateLoadTest : public ::testing::Test {};

TEST_F(ParaViewStateLoadTest, NearFieldCavityStateLoads)
{
    // Locate the script. ctest runs each test in its own working
    // directory; we walk up from the CWD looking for the project
    // root layout (a sibling "scripts" directory containing the
    // verification script).
    namespace fs = std::filesystem;
    fs::path script;
    fs::path cur = fs::current_path();
    for (int up = 0; up < 6; ++up) {
        const fs::path candidate = cur / "scripts" / "verify_pvsm.sh";
        if (fs::exists(candidate)) {
            script = candidate;
            break;
        }
        if (!cur.has_parent_path()) break;
        cur = cur.parent_path();
    }
    ASSERT_FALSE(script.empty())
        << "scripts/verify_pvsm.sh not found by walking parent paths "
        << "from CWD=" << fs::current_path();

    std::string cmd = "bash \"" + script.string() + "\" > /dev/null 2>&1";
    const int rc = std::system(cmd.c_str());
    const int exit_code = WEXITSTATUS(rc);

    if (exit_code == 2) {
        GTEST_SKIP()
            << "pvpython not available in the test environment "
            << "(scripts/verify_pvsm.sh exited 2). Install ParaView "
            << "5.10+ to enable this verification.";
    }

    ASSERT_EQ(exit_code, 0)
        << "scripts/verify_pvsm.sh failed with exit code "
        << exit_code << "; the .pvsm did not load cleanly. Run the "
        << "script manually for full diagnostics: " << script;
}
