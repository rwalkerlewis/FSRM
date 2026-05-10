/**
 * @file test_drucker_prager_3d.cpp
 * @brief Pass-13b (axis-1b physics) unit tests for the 3D Drucker-Prager
 *        explicit radial-return projection used inside the
 *        Source3DBallImpl per-cell update.
 *
 * Gates:
 *   1. PureElasticPredictorNoPlasticStrain: trial stress below yield
 *      passes through unchanged.
 *   2. YieldHitProjectsToSurface: trial stress above yield projects to
 *      the surface within machine precision.
 *   3. SphericalSymmetryReducesTo1D: with isotropic IC plus a radial
 *      deviator the 3D return reproduces the 1D scalar return.
 *   4. MediumParameterRoundTrip: each per-medium preset loads and
 *      produces the documented yield-surface scale.
 */

#include <gtest/gtest.h>

#include <array>
#include <cmath>

#include "domain/explosion/DruckerPrager3D.hpp"

using FSRM::DruckerPrager3DParameters;
using FSRM::DruckerPrager3DReturnResult;
using FSRM::druckerPrager3DRadialReturn;
using FSRM::Voigt6::deviator;
using FSRM::Voigt6::sqrtJ2;
using FSRM::Voigt6::trace;

namespace
{

constexpr double K_TEST = 30.0e9;  // 30 GPa bulk
constexpr double G_TEST = 18.0e9;  // 18 GPa shear

}  // namespace

TEST(DruckerPrager3D, PureElasticPredictorNoPlasticStrain)
{
    // Initial state: small isotropic compression.
    std::array<double, 6> sigma = {-1.0e6, -1.0e6, -1.0e6, 0, 0, 0};
    // Tiny strain increment that keeps trial well inside yield surface.
    std::array<double, 6> deps = {-1.0e-6, -1.0e-6, -1.0e-6, 0, 0, 0};

    DruckerPrager3DParameters params;
    params.alpha_dp = 0.3;
    params.k_dp = 50.0e6;

    std::array<double, 6> sigma_out{};
    std::array<double, 6> deps_p_out{};
    auto res = druckerPrager3DRadialReturn(sigma, deps, K_TEST, G_TEST,
                                           params, sigma_out, deps_p_out);

    EXPECT_FALSE(res.yielded);
    EXPECT_DOUBLE_EQ(res.delta_lambda, 0.0);
    EXPECT_DOUBLE_EQ(res.delta_eps_p_eq, 0.0);
    for (int i = 0; i < 6; ++i) {
        EXPECT_DOUBLE_EQ(deps_p_out[i], 0.0);
    }
    // Stress is the elastic-predicted stress.
    const double tr = -1.0e-6 - 1.0e-6 - 1.0e-6;
    const double expected_diag = -1.0e6 + (K_TEST - 2.0/3.0 * G_TEST) * tr
                                 + 2.0 * G_TEST * (-1.0e-6);
    EXPECT_NEAR(sigma_out[0], expected_diag, 1.0e-3);
    EXPECT_NEAR(sigma_out[1], expected_diag, 1.0e-3);
    EXPECT_NEAR(sigma_out[2], expected_diag, 1.0e-3);
    EXPECT_DOUBLE_EQ(sigma_out[3], 0.0);
    EXPECT_DOUBLE_EQ(sigma_out[4], 0.0);
    EXPECT_DOUBLE_EQ(sigma_out[5], 0.0);
}

TEST(DruckerPrager3D, YieldHitProjectsToSurface)
{
    // Initial state at zero stress.
    std::array<double, 6> sigma = {0, 0, 0, 0, 0, 0};
    // Apply pure shear strain large enough to push past yield.
    std::array<double, 6> deps = {0, 0, 0, 0.01, 0, 0};

    DruckerPrager3DParameters params;
    params.alpha_dp = 0.3;
    params.k_dp = 50.0e6;

    std::array<double, 6> sigma_out{};
    std::array<double, 6> deps_p_out{};
    auto res = druckerPrager3DRadialReturn(sigma, deps, K_TEST, G_TEST,
                                           params, sigma_out, deps_p_out);

    EXPECT_TRUE(res.yielded);
    EXPECT_GT(res.delta_lambda, 0.0);
    EXPECT_GT(res.delta_eps_p_eq, 0.0);

    // Projected stress lies on the yield surface.
    const double I1 = trace(sigma_out);
    const auto s = deviator(sigma_out);
    const double sJ2 = sqrtJ2(s);
    const double f_after = sJ2 - params.alpha_dp * I1 / 3.0 - params.k_dp;
    EXPECT_NEAR(f_after, 0.0, 1.0e-3);
}

TEST(DruckerPrager3D, SphericalSymmetryReducesTo1D)
{
    // With pure-Mises (alpha=0) and an isotropic IC plus a radial deviator
    // s_rr = +S, s_tt = s_phiphi = -S/2, the 3D return should produce the
    // same projected stress as the 1D scalar return on s_rr.
    //
    // Drive with a strain rate that mimics spherical compression.
    DruckerPrager3DParameters params;
    params.alpha_dp = 0.0;
    params.k_dp = 30.0e6;

    // Isotropic initial compressive stress + radial deviator.
    const double p_iso = -50.0e6;
    const double S = 25.0e6;  // s_rr
    std::array<double, 6> sigma = {
        p_iso + S,
        p_iso - 0.5 * S,
        p_iso - 0.5 * S,
        0, 0, 0
    };

    // Radial-symmetric strain increment: eps_rr = -2e-3, eps_tt = -1e-3.
    // (Larger compression along radial than tangential, drives more deviator.)
    std::array<double, 6> deps = {-2.0e-3, -1.0e-3, -1.0e-3, 0, 0, 0};

    std::array<double, 6> sigma_out{};
    std::array<double, 6> deps_p_out{};
    auto res = druckerPrager3DRadialReturn(sigma, deps, K_TEST, G_TEST,
                                           params, sigma_out, deps_p_out);

    EXPECT_TRUE(res.yielded);

    // Compute the 1D-equivalent: the trial deviator s_rr_trial =
    //   s_rr_in + (4G/3)(eps_rr - eps_tt) dt
    // Then sigma_eq = 1.5 * |s_rr_trial|; if > Y, scale by Y/sigma_eq.
    const double ds_rr = (4.0 / 3.0) * G_TEST * (-2.0e-3 - (-1.0e-3));
    const double s_rr_trial = S + ds_rr;
    const double sigma_eq = 1.5 * std::abs(s_rr_trial);
    const double Y = params.k_dp;
    ASSERT_GT(sigma_eq, Y) << "Test driver chose subyield trial; tighten "
                              "strain increment.";
    const double scale = Y / sigma_eq;
    const double s_rr_1d = s_rr_trial * scale;
    const double s_tt_1d = -0.5 * s_rr_1d;

    // Recover the deviator out of the 3D result.
    const auto s_out = deviator(sigma_out);
    const double s_rr_3d = s_out[0];
    const double s_tt_3d = s_out[1];

    EXPECT_NEAR(s_rr_3d, s_rr_1d, 1.0e-3 * std::abs(s_rr_1d) + 1.0)
        << "3D s_rr=" << s_rr_3d << " vs 1D s_rr=" << s_rr_1d;
    EXPECT_NEAR(s_tt_3d, s_tt_1d, 1.0e-3 * std::abs(s_tt_1d) + 1.0)
        << "3D s_tt=" << s_tt_3d << " vs 1D s_tt=" << s_tt_1d;
}

TEST(DruckerPrager3D, MediumParameterRoundTrip)
{
    using namespace FSRM::DruckerPrager3DSets;
    // Each preset loads with positive alpha and k. The granite preset
    // is the default fallback for unrecognised names.
    const auto p_g = granite();
    const auto p_t = tuff();
    const auto p_s = salt();
    const auto p_a = alluvium();
    EXPECT_GT(p_g.alpha_dp, 0.0);  EXPECT_GT(p_g.k_dp, 0.0);
    EXPECT_GT(p_t.alpha_dp, 0.0);  EXPECT_GT(p_t.k_dp, 0.0);
    EXPECT_GT(p_s.alpha_dp, 0.0);  EXPECT_GT(p_s.k_dp, 0.0);
    EXPECT_GT(p_a.alpha_dp, 0.0);  EXPECT_GT(p_a.k_dp, 0.0);
    EXPECT_EQ(p_g.medium_label, "GRANITE");
    EXPECT_EQ(p_t.medium_label, "TUFF");
    EXPECT_EQ(p_s.medium_label, "SALT");
    EXPECT_EQ(p_a.medium_label, "ALLUVIUM");

    // byName dispatch.
    EXPECT_EQ(byName("GRANITE").medium_label, "GRANITE");
    EXPECT_EQ(byName("TUFF").medium_label, "TUFF");
    EXPECT_EQ(byName("SALT").medium_label, "SALT");
    EXPECT_EQ(byName("ALLUVIUM").medium_label, "ALLUVIUM");
    EXPECT_EQ(byName("UNKNOWN_MEDIUM").medium_label, "GRANITE");

    // Each preset gives a distinct yield-surface scale: granite > tuff >
    // salt at the same confining pressure (k decreases from granite to
    // salt; alpha is comparable).
    const double P_ref = 100.0e6;
    auto Y_at = [&](const DruckerPrager3DParameters& p) {
        return p.alpha_dp * P_ref + p.k_dp;
    };
    EXPECT_GT(Y_at(p_g), Y_at(p_t));
    EXPECT_GT(Y_at(p_t), Y_at(p_a));
    // Salt has lowest both alpha and k.
    EXPECT_LT(Y_at(p_s), Y_at(p_g));
}
