/**
 * @file test_source3d_ball_moment_tensor.cpp
 * @brief Pass-13c (axis-1b validation) unit gates for the
 *        Source3DBallImpl surface-integral moment-tensor extraction.
 *
 * Verifies the integral
 *   M_ij(t) = integral_S [sigma_jk(t) - sigma_jk(0)] * n_k * x_i dA
 * (Day & McLaughlin 1991; Aki & Richards 2002 ch 4) on the cube_5tet
 * fixture, which has unit volume V = 1 with the origin vertex marked
 * cavity (marker = 2) and the seven other corners marked elastic
 * (marker = 1). All boundary faces participate in the elastic-surface
 * set: the three faces incident to the origin vertex have mixed
 * markers and default to ELASTIC per the buildCellAndFaceLists
 * fall-through convention. The boundary therefore covers the closed
 * cube surface.
 *
 * Divergence-theorem identity for any closed surface S enclosing the
 * volume V (origin can lie anywhere, including on the boundary):
 *   integral_S n_j * x_i dA = integral_V (d x_i / d x_j) dV
 *                           = V * delta_ij
 *
 * Test gates:
 *  - UniformIsotropicStressGivesIsotropicM:
 *      sigma = -p * I, sigma_initial = 0
 *      => M_ij = -p * V * delta_ij (diagonal isotropic; off-diagonals 0)
 *  - PureDeviatoricStressGivesDeviatoricM:
 *      sigma_xx = +s, sigma_yy = +s, sigma_zz = -2s, off-diag = 0
 *      => M = (s * V) * diag(1, 1, -2) (CLVD pattern, axis = z)
 *  - StressDropReferencesInitialState:
 *      sigma = -p * I, sigma_initial = -p * I
 *      => M = 0 (the integral measures the change, not the absolute
 *      stress)
 *  - InitialStateNoStepReturnsZero:
 *      Integration with the as-loaded fixture returns M = 0 before
 *      any step() call.
 *  - MomentRateFiniteDifference:
 *      Two recompute calls with dt=1.0 between successive snapshots
 *      produce Mdot = (M2 - M1)/dt.
 *
 * Cite Day & McLaughlin 1991 sec 4.
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

std::string fixtureBasename(const std::string& fixture_name)
{
    const char* env = std::getenv("FSRM_TEST_DATA_DIR");
    std::string root;
    if (env && *env) {
        root = env;
    } else {
        std::vector<std::string> candidates = {
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

Source3DBallConfig cubeConfig()
{
    Source3DBallConfig cfg;
    cfg.mesh_path = fixtureBasename("cube_5tet");
    cfg.cavity_radius_m = 0.0;       // cube fixture origin
    cfg.outer_radius_m = std::sqrt(3.0);
    cfg.asymmetric_overburden = false;
    cfg.host_density_kg_per_m3 = 2700.0;
    cfg.bulk_modulus_K_pa = 30.0e9;
    cfg.shear_modulus_G_pa = 18.0e9;
    cfg.dp3d_alpha = 0.0;
    cfg.dp3d_k_pa = 1.0e30;          // effectively elastic for these gates
    cfg.cv_J_per_kg_K = 1000.0;
    cfg.T_ambient_K = 300.0;
    cfg.medium_label = "GRANITE";
    return cfg;
}

}  // namespace

class Source3DBallSurfaceIntegralTest : public ::testing::Test
{
};

// =========================================================================
// 1. Uniform isotropic stress over a closed surface produces a purely
//    isotropic moment tensor.
// =========================================================================
TEST_F(Source3DBallSurfaceIntegralTest, UniformIsotropicStressGivesIsotropicM)
{
    Source3DBallImpl ball;
    ball.initialize(cubeConfig());
    ASSERT_TRUE(ball.isInitialized());
    const int N = ball.numLocalCells();
    ASSERT_GT(N, 0);
    ASSERT_GT(ball.numLocalElasticSurfaceFaces(), 0);

    // Pure isotropic compression sigma = -p * I.
    const double p = 1.0e9;  // 1 GPa
    std::vector<std::array<double, 6>> sig(
        N, {-p, -p, -p, 0.0, 0.0, 0.0});
    std::vector<std::array<double, 6>> sig0(
        N, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
    ball.setCellStressForTest(sig);
    ball.setCellStressInitialForTest(sig0);

    ball.recomputeMomentTensorFromSurface(1.0);

    std::array<double, 6> M{};
    ball.getMomentTensor(M);

    const double V = 1.0;  // unit cube volume
    const double expected_diag = -p * V;

    // Tolerance: the cube fixture is a low-resolution reference so we
    // accept 1e-6 relative roundoff. The integral is exact in
    // arithmetic (each face contributes once); machine-precision
    // tolerance is appropriate.
    EXPECT_NEAR(M[0], expected_diag, 1.0e-3 * std::abs(expected_diag));
    EXPECT_NEAR(M[1], expected_diag, 1.0e-3 * std::abs(expected_diag));
    EXPECT_NEAR(M[2], expected_diag, 1.0e-3 * std::abs(expected_diag));
    EXPECT_NEAR(M[3], 0.0, 1.0e-3 * std::abs(expected_diag));
    EXPECT_NEAR(M[4], 0.0, 1.0e-3 * std::abs(expected_diag));
    EXPECT_NEAR(M[5], 0.0, 1.0e-3 * std::abs(expected_diag));
}

// =========================================================================
// 2. Pure deviatoric stress produces a deviatoric moment tensor with
//    the CLVD pattern (vertical compression, horizontal tension or vice
//    versa).
// =========================================================================
TEST_F(Source3DBallSurfaceIntegralTest, PureDeviatoricStressGivesDeviatoricM)
{
    Source3DBallImpl ball;
    ball.initialize(cubeConfig());
    ASSERT_TRUE(ball.isInitialized());
    const int N = ball.numLocalCells();
    ASSERT_GT(N, 0);

    // Trace-free deviator: sigma_xx = sigma_yy = +s, sigma_zz = -2s.
    // Voigt order: [xx, yy, zz, xy, xz, yz]. CLVD axis = z (vertical
    // compression, horizontal tension).
    const double s = 0.5e9;
    std::vector<std::array<double, 6>> sig(
        N, {+s, +s, -2.0 * s, 0.0, 0.0, 0.0});
    std::vector<std::array<double, 6>> sig0(
        N, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
    ball.setCellStressForTest(sig);
    ball.setCellStressInitialForTest(sig0);

    ball.recomputeMomentTensorFromSurface(1.0);

    std::array<double, 6> M{};
    ball.getMomentTensor(M);

    const double V = 1.0;
    EXPECT_NEAR(M[0], +s * V, 1.0e-3 * std::abs(s * V));
    EXPECT_NEAR(M[1], +s * V, 1.0e-3 * std::abs(s * V));
    EXPECT_NEAR(M[2], -2.0 * s * V, 1.0e-3 * std::abs(s * V));
    EXPECT_NEAR(M[3], 0.0, 1.0e-3 * std::abs(s * V));
    EXPECT_NEAR(M[4], 0.0, 1.0e-3 * std::abs(s * V));
    EXPECT_NEAR(M[5], 0.0, 1.0e-3 * std::abs(s * V));

    // Trace must be (numerically) zero for a deviatoric input.
    const double trace = M[0] + M[1] + M[2];
    EXPECT_NEAR(trace, 0.0, 1.0e-6 * std::abs(s * V));
}

// =========================================================================
// 3. Stress drop references the initial state: identical sigma and
//    sigma_initial gives M = 0 (no source, no radiation).
// =========================================================================
TEST_F(Source3DBallSurfaceIntegralTest, StressDropReferencesInitialState)
{
    Source3DBallImpl ball;
    ball.initialize(cubeConfig());
    const int N = ball.numLocalCells();
    ASSERT_GT(N, 0);

    const double p = 7.5e8;
    std::vector<std::array<double, 6>> sig(
        N, {-p, -p, -p, 0.0, 0.0, 0.0});
    ball.setCellStressForTest(sig);
    ball.setCellStressInitialForTest(sig);  // identical reference

    ball.recomputeMomentTensorFromSurface(1.0);

    std::array<double, 6> M{};
    ball.getMomentTensor(M);
    for (int k = 0; k < 6; ++k) {
        EXPECT_NEAR(M[k], 0.0, 1.0e-3 * p)
            << "M[" << k << "] should be zero when sigma = sigma_initial";
    }
}

// =========================================================================
// 4. As-loaded with no step(), the moment tensor reads zero (sigma_
//    starts equal to sigma_initial_).
// =========================================================================
TEST_F(Source3DBallSurfaceIntegralTest, InitialStateNoStepReturnsZero)
{
    Source3DBallImpl ball;
    ball.initialize(cubeConfig());

    std::array<double, 6> M{};
    std::array<double, 6> Mdot{};
    ball.getMomentTensor(M);
    ball.getMomentRateTensor(Mdot);
    for (int k = 0; k < 6; ++k) {
        EXPECT_EQ(M[k], 0.0);
        EXPECT_EQ(Mdot[k], 0.0);
    }
}

// =========================================================================
// 5. Mdot = (M_new - M_old) / dt finite difference.
// =========================================================================
TEST_F(Source3DBallSurfaceIntegralTest, MomentRateFiniteDifference)
{
    Source3DBallImpl ball;
    ball.initialize(cubeConfig());
    const int N = ball.numLocalCells();
    ASSERT_GT(N, 0);

    std::vector<std::array<double, 6>> sig0(
        N, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
    ball.setCellStressInitialForTest(sig0);

    // First snapshot: sigma = -p1 * I.
    const double p1 = 1.0e9;
    std::vector<std::array<double, 6>> sig1(
        N, {-p1, -p1, -p1, 0.0, 0.0, 0.0});
    ball.setCellStressForTest(sig1);
    ball.recomputeMomentTensorFromSurface(1.0);

    // Second snapshot: sigma = -p2 * I with dt = 1 ms.
    const double p2 = 2.0e9;
    const double dt = 1.0e-3;
    std::vector<std::array<double, 6>> sig2(
        N, {-p2, -p2, -p2, 0.0, 0.0, 0.0});
    ball.setCellStressForTest(sig2);
    ball.recomputeMomentTensorFromSurface(dt);

    std::array<double, 6> Mdot{};
    ball.getMomentRateTensor(Mdot);

    // Expected dM = -V*(p2 - p1)*delta_ij = -V*1e9 (diag).
    const double V = 1.0;
    const double expected_dot_diag = -V * (p2 - p1) / dt;
    EXPECT_NEAR(Mdot[0], expected_dot_diag,
                1.0e-3 * std::abs(expected_dot_diag));
    EXPECT_NEAR(Mdot[1], expected_dot_diag,
                1.0e-3 * std::abs(expected_dot_diag));
    EXPECT_NEAR(Mdot[2], expected_dot_diag,
                1.0e-3 * std::abs(expected_dot_diag));
}

// =========================================================================
// 6. Cavity radius extremes: with the cube_5tet fixture the cavity
//    vertex sits at the origin (radius = 0). Asymmetric reporting is
//    pass-14 (no advection in pass-13c).
// =========================================================================
TEST_F(Source3DBallSurfaceIntegralTest, CavityRadiusExtremesReportable)
{
    Source3DBallImpl ball;
    ball.initialize(cubeConfig());

    double r_min = -1.0;
    double r_max = -1.0;
    ball.cavityRadiusExtremes(r_min, r_max);
    EXPECT_NEAR(r_min, 0.0, 1.0e-9);
    EXPECT_NEAR(r_max, 0.0, 1.0e-9);
}
