/**
 * @file test_pass14a_source_forcing.cpp
 * @brief Pass-14a (axis-1b source-forcing slice) physics-validation
 *        gates for the 3D source ball.
 *
 * The pass-14a source forcing deposits the explosion yield as
 * mechanical energy into the inner-cavity cells (Mueller & Murphy 1971
 * source-time function), turns it into a Tillotson matter-pressure rise,
 * spreads it through the rock as the Sharpe (1942) / Lame elastostatic
 * pressurised-cavity field, and the pass-13c surface integral over the
 * elastic radius reads it back as a real moment-rate tensor (Day &
 * McLaughlin 1991; Aki & Richards 2002 ch 4).
 *
 * Gates:
 *   SourceForcingProducesNonZeroMoment (HEADLINE)
 *     -- with Salmon's yield (5.3 kt) the moment-tensor trace at the
 *        elastic radius exceeds 1e-2 * (yield * 4.184e12 J) by the end
 *        of the deposition window. (Pass-13c had M = 0 here; this is
 *        the headline non-regression-to-non-trivial flip.)
 *   SphericalSymmetrySourceForcingGivesIsotropicMoment
 *     -- with an isotropic IC (asymmetric_overburden = false) the
 *        moment is diagonal and isotropic: the deviator is below 10
 *        percent of the isotropic part, confirming the source forcing
 *        does not introduce spurious asymmetry. (The literature-range
 *        CLVD content with K_0 < 1 needs the pass-14b cavity-asymmetry
 *        physics; the pass-14a stress-only contribution under pure
 *        compression stays purely isotropic, which this gate documents.)
 */

#include <gtest/gtest.h>

#include <petscsys.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <vector>

#include "domain/explosion/Source3DBall.hpp"
#include "domain/explosion/Source3DBallImpl.hpp"

using FSRM::Source3DBallConfig;
using FSRM::Source3DBallImpl;

namespace
{

constexpr double kJoulesPerKt = 4.184e12;

std::string fixtureBasename(const std::string& fixture_name)
{
    const char* env = std::getenv("FSRM_TEST_DATA_DIR");
    std::string root;
    if (env && *env) {
        root = env;
    } else {
        const std::vector<std::string> candidates = {
            "../tests/data/source_ball",
            "../../tests/data/source_ball",
            "tests/data/source_ball",
            "./tests/data/source_ball",
        };
        for (const auto& c : candidates) {
            if (std::filesystem::exists(std::filesystem::path(c))) {
                root = c;
                break;
            }
        }
        if (root.empty()) root = "../tests/data/source_ball";
    }
    return root + "/" + fixture_name + "/source_ball";
}

bool salmonMeshAvailable()
{
    const std::string base = fixtureBasename("salmon_default");
    return std::filesystem::exists(base + ".node")
           && std::filesystem::exists(base + ".ele");
}

Source3DBallConfig salmonCfg()
{
    Source3DBallConfig cfg;
    cfg.mesh_path = fixtureBasename("salmon_default");
    cfg.cavity_radius_m = 17.0;
    cfg.outer_radius_m = 70.0;
    cfg.host_density_kg_per_m3 = 2160.0;     // Tatum salt
    cfg.bulk_modulus_K_pa = 30.0e9;
    cfg.shear_modulus_G_pa = 18.0e9;
    cfg.dp3d_alpha = 0.30;
    cfg.dp3d_k_pa = 70.0e6;
    cfg.cv_J_per_kg_K = 850.0;
    cfg.T_ambient_K = 300.0;
    cfg.source_depth_m = 828.0;
    cfg.medium_label = "SALT";
    cfg.radiation_substep_cadence = 100;     // implicit, unconditionally stable
    cfg.source_forcing_enabled = true;
    cfg.source_yield_kt = 5.3;
    cfg.source_time_function = Source3DBallConfig::SourceTimeFunction::MUELLER_MURPHY;
    cfg.source_deposition_duration_s = 1.0e-3;
    cfg.source_deposition_efficiency = 1.0;
    return cfg;
}

// Step the ball past the end of the deposition window.
void runDepositionWindow(Source3DBallImpl& impl, double duration_s)
{
    const double dt = 1.0e-4;                 // 10 steps per ms duration
    const int n = static_cast<int>(std::ceil(3.0 * duration_s / dt)) + 5;
    for (int k = 0; k < n; ++k) impl.step(dt);
}

double traceVoigt(const std::array<double, 6>& M) { return M[0] + M[1] + M[2]; }

double deviatorFrobeniusVoigt(const std::array<double, 6>& M)
{
    const double tr3 = (M[0] + M[1] + M[2]) / 3.0;
    const double dxx = M[0] - tr3;
    const double dyy = M[1] - tr3;
    const double dzz = M[2] - tr3;
    return std::sqrt(dxx * dxx + dyy * dyy + dzz * dzz
                     + 2.0 * (M[3] * M[3] + M[4] * M[4] + M[5] * M[5]));
}

}  // namespace

class Pass14aSourceForcing : public ::testing::Test {};

// =========================================================================
// HEADLINE: source forcing produces a non-zero moment tensor.
// =========================================================================
TEST_F(Pass14aSourceForcing, SourceForcingProducesNonZeroMoment)
{
    if (!salmonMeshAvailable()) {
        GTEST_SKIP() << "salmon_default source-ball mesh fixture not found "
                        "under tests/data/source_ball/salmon_default/";
    }
    Source3DBallConfig cfg = salmonCfg();
    cfg.asymmetric_overburden = true;
    cfg.overburden_K0 = 0.5;
    Source3DBallImpl impl;
    impl.initialize(cfg);
    ASSERT_TRUE(impl.sourceForcingActive());
    ASSERT_GT(impl.numInnerCavityCells(), 0);
    ASSERT_GT(impl.innerCavityVolumeGlobal(), 0.0);

    runDepositionWindow(impl, cfg.source_deposition_duration_s);

    const double E_total = cfg.source_yield_kt * kJoulesPerKt;
    EXPECT_NEAR(impl.depositedEnergyGlobalJ(), E_total, 1.0e-2 * E_total)
        << "after the deposition window the cumulative deposited energy "
           "should be (essentially) the full yield energy";

    std::array<double, 6> M{};
    impl.getMomentTensor(M);
    const double tr = traceVoigt(M);
    EXPECT_GT(std::abs(tr), 1.0e-2 * E_total)
        << "moment-tensor trace |" << tr << "| must exceed 1e-2 * E_total ("
        << 1.0e-2 * E_total << " N*m); pass-13c had M = 0";
    for (double m : M) EXPECT_TRUE(std::isfinite(m));
}

// =========================================================================
// Spherical symmetry: an isotropic IC -> a (near) isotropic moment.
// =========================================================================
TEST_F(Pass14aSourceForcing, SphericalSymmetrySourceForcingGivesIsotropicMoment)
{
    if (!salmonMeshAvailable()) {
        GTEST_SKIP() << "salmon_default source-ball mesh fixture not found "
                        "under tests/data/source_ball/salmon_default/";
    }
    Source3DBallConfig cfg = salmonCfg();
    cfg.asymmetric_overburden = false;        // zero IC -> no IC deviator
    cfg.overburden_K0 = 1.0;
    Source3DBallImpl impl;
    impl.initialize(cfg);
    ASSERT_TRUE(impl.sourceForcingActive());

    runDepositionWindow(impl, cfg.source_deposition_duration_s);

    std::array<double, 6> M{};
    impl.getMomentTensor(M);
    const double tr = traceVoigt(M);
    ASSERT_GT(std::abs(tr), 0.0);
    const double dev = deviatorFrobeniusVoigt(M);
    // The deviator should be a small fraction of the isotropic part.
    // The residual comes only from the icosphere tessellation not being
    // perfectly spherical (the per-cell stress drop is purely isotropic
    // under pure compression), so 10 percent of |trace| is comfortable.
    EXPECT_LT(dev, 0.10 * std::abs(tr))
        << "deviator " << dev << " should be below 10 percent of |trace| "
        << std::abs(tr) << " for the isotropic-IC source forcing";
}

// =========================================================================
// Backward compat: source_forcing_enabled = false reproduces the
// pass-13c "no source signal" behaviour under cavity_geometry =
// THREE_DIMENSIONAL (the regression guard for the new feature).
// =========================================================================
TEST_F(Pass14aSourceForcing, ThreeDimensionalWithoutSourceForcingGivesZeroMoment)
{
    if (!salmonMeshAvailable()) {
        GTEST_SKIP() << "salmon_default source-ball mesh fixture not found";
    }
    Source3DBallConfig cfg = salmonCfg();
    cfg.asymmetric_overburden = true;
    cfg.overburden_K0 = 0.5;
    cfg.source_forcing_enabled = false;       // the regression-guard mode
    Source3DBallImpl impl;
    impl.initialize(cfg);
    EXPECT_FALSE(impl.sourceForcingActive());
    EXPECT_EQ(impl.numInnerCavityCells(), 0);

    for (int k = 0; k < 10; ++k) impl.step(1.0e-4);

    EXPECT_DOUBLE_EQ(impl.depositedEnergyGlobalJ(), 0.0);
    std::array<double, 6> M{};
    impl.getMomentTensor(M);
    for (double m : M) EXPECT_DOUBLE_EQ(m, 0.0)
        << "no source forcing -> the surface-integral moment must stay "
           "exactly zero (pass-13c behaviour)";
}
