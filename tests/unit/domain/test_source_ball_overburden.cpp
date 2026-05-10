/**
 * @file test_source_ball_overburden.cpp
 * @brief Pass-13b unit tests for the asymmetric overburden initial state.
 */

#include <gtest/gtest.h>

#include <array>
#include <cmath>

#include "domain/explosion/SourceBallOverburden.hpp"

using FSRM::SourceBallOverburdenConfig;
using FSRM::overburdenStressAtCell;

TEST(Overburden, K0EqualsOneIsotropic)
{
    SourceBallOverburdenConfig cfg;
    cfg.source_depth_m = 500.0;
    cfg.rho_solid = 2700.0;
    cfg.g = 9.81;
    cfg.K_0 = 1.0;

    const std::array<std::array<double, 3>, 3> samples = {{
        {{0.0, 0.0, 0.0}},
        {{5.0, -3.0, 7.0}},
        {{-2.0, 1.5, -4.0}}
    }};
    for (const auto& xyz : samples) {
        auto s = overburdenStressAtCell(xyz, cfg);
        EXPECT_NEAR(s[0], s[2], 1.0e-9 * std::abs(s[2]) + 1.0)
            << "K_0=1 must give isotropic stress: "
            << "sxx=" << s[0] << " szz=" << s[2];
        EXPECT_NEAR(s[1], s[2], 1.0e-9 * std::abs(s[2]) + 1.0);
        EXPECT_DOUBLE_EQ(s[3], 0.0);
        EXPECT_DOUBLE_EQ(s[4], 0.0);
        EXPECT_DOUBLE_EQ(s[5], 0.0);
    }
}

TEST(Overburden, K0EqualsHalfAnisotropic)
{
    SourceBallOverburdenConfig cfg;
    cfg.source_depth_m = 500.0;
    cfg.rho_solid = 2700.0;
    cfg.g = 9.81;
    cfg.K_0 = 0.5;

    const std::array<double, 3> xyz = {2.0, -1.0, 4.0};
    auto s = overburdenStressAtCell(xyz, cfg);
    ASSERT_LT(s[2], 0.0) << "Compression must be negative.";
    EXPECT_NEAR(s[0] / s[2], 0.5, 1.0e-9);
    EXPECT_NEAR(s[1] / s[2], 0.5, 1.0e-9);
}

TEST(Overburden, DepthGradientFromGravity)
{
    SourceBallOverburdenConfig cfg;
    cfg.source_depth_m = 500.0;
    cfg.rho_solid = 2700.0;
    cfg.g = 9.81;
    cfg.K_0 = 1.0;

    // Two cells at the same horizontal location, different z. Their
    // sigma_zz difference equals -rho * g * (d_high - d_low) = -rho * g
    // * (-(z_high - z_low)) = +rho * g * (z_high - z_low) but with sign:
    // higher cell (larger z) has shallower depth so smaller compressive
    // stress in magnitude. Gradient d sigma_zz / dz = +rho * g.
    const double dz = 1.0;
    auto s_low  = overburdenStressAtCell({0.0, 0.0, 0.0}, cfg);
    auto s_high = overburdenStressAtCell({0.0, 0.0, dz},  cfg);
    const double dszz = s_high[2] - s_low[2];
    const double expected = cfg.rho_solid * cfg.g * dz;
    EXPECT_NEAR(dszz, expected, 1.0e-3 * expected);
}
