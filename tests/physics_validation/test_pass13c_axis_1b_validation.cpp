/**
 * @file test_pass13c_axis_1b_validation.cpp
 * @brief Pass-13c (axis-1b validation) physics gates.
 *
 * Surface-integral moment-tensor decomposition. Given the symmetric
 * 6-component Voigt tensor M = [Mxx, Myy, Mzz, Mxy, Mxz, Myz], the
 * isotropic part is M_iso = (Mxx + Myy + Mzz) / 3 * I and the
 * deviator is M_dev = M - M_iso * I. A pure CLVD pattern (Knopoff &
 * Randall 1970; Stevens & Day 1985) has the deviator decomposable as
 *   M_dev = M_clvd * (e_3 e_3^T - (1/2) (e_1 e_1^T + e_2 e_2^T))
 * for some axis e_3. The decomposition reduces (under principal-axis
 * alignment) to two equal eigenvalues +M_clvd/2 and one eigenvalue
 * -M_clvd along the CLVD axis (or vice versa).
 *
 * Headline gate: CLVDContentWithOverburden
 *
 *   With K_0 = 0.5 the IC stress has the deviator
 *     s_xx = s_yy = (1/6) rho g z, s_zz = -(1/3) rho g z
 *   so the source field starts with an axis-z CLVD pattern. Under any
 *   isotropic compression strain rate the elastic predictor preserves
 *   the deviator; the Drucker-Prager radial return scales it; after
 *   N steps the integrated stress drop carries a residual axis-z CLVD
 *   pattern.
 *
 *   Gate asserts:
 *     |CLVD| / |isotropic| > 0.01 (best-effort; without 3D Lagrangian
 *       advection the geometric CLVD source is absent and the
 *       stress-only contribution is small).
 *     CLVD principal axis aligns with vertical (z) within 10 degrees.
 *
 * Cite Day & McLaughlin 1991 sec 4 (extraction); Stevens & Day 1985
 * (CLVD content interpretation).
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

/// 3x3 symmetric eigendecomposition via the analytic closed form.
/// Returns the three eigenvalues sorted ascending plus the
/// corresponding orthonormal eigenvectors (columns). The matrix is
/// passed in Voigt order [xx, yy, zz, xy, xz, yz].
struct Eig3
{
    std::array<double, 3> values;          ///< sorted ascending
    std::array<std::array<double, 3>, 3> vectors;  ///< columns
};

/// Jacobi diagonalization for a symmetric 3x3 matrix. Robust enough
/// for test purposes; the input matrices are well-conditioned.
Eig3 eigSymm3(const std::array<double, 6>& voigt)
{
    double a[3][3] = {
        {voigt[0], voigt[3], voigt[4]},
        {voigt[3], voigt[1], voigt[5]},
        {voigt[4], voigt[5], voigt[2]}
    };
    double v[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    for (int iter = 0; iter < 100; ++iter) {
        // Find largest off-diagonal.
        int p = 0, q = 1;
        double max_off = std::abs(a[0][1]);
        if (std::abs(a[0][2]) > max_off) { p = 0; q = 2; max_off = std::abs(a[0][2]); }
        if (std::abs(a[1][2]) > max_off) { p = 1; q = 2; max_off = std::abs(a[1][2]); }
        if (max_off < 1.0e-14 * (std::abs(a[0][0]) + std::abs(a[1][1])
                                  + std::abs(a[2][2]) + 1.0)) break;
        const double app = a[p][p];
        const double aqq = a[q][q];
        const double apq = a[p][q];
        const double theta = (aqq - app) / (2.0 * apq);
        const double t = (theta >= 0.0)
            ? 1.0 / (theta + std::sqrt(1.0 + theta * theta))
            : 1.0 / (theta - std::sqrt(1.0 + theta * theta));
        const double c = 1.0 / std::sqrt(1.0 + t * t);
        const double s = t * c;
        // Rotate a.
        a[p][p] = app - t * apq;
        a[q][q] = aqq + t * apq;
        a[p][q] = a[q][p] = 0.0;
        for (int k = 0; k < 3; ++k) {
            if (k == p || k == q) continue;
            const double akp = a[k][p];
            const double akq = a[k][q];
            a[k][p] = a[p][k] = c * akp - s * akq;
            a[k][q] = a[q][k] = s * akp + c * akq;
        }
        // Rotate v.
        for (int k = 0; k < 3; ++k) {
            const double vkp = v[k][p];
            const double vkq = v[k][q];
            v[k][p] = c * vkp - s * vkq;
            v[k][q] = s * vkp + c * vkq;
        }
    }
    Eig3 out;
    out.values = {a[0][0], a[1][1], a[2][2]};
    out.vectors[0] = {v[0][0], v[1][0], v[2][0]};
    out.vectors[1] = {v[0][1], v[1][1], v[2][1]};
    out.vectors[2] = {v[0][2], v[1][2], v[2][2]};
    // Sort ascending.
    int idx[3] = {0, 1, 2};
    for (int i = 0; i < 3; ++i)
        for (int j = i + 1; j < 3; ++j)
            if (out.values[idx[j]] < out.values[idx[i]]) std::swap(idx[i], idx[j]);
    Eig3 sorted = out;
    for (int i = 0; i < 3; ++i) {
        sorted.values[i] = out.values[idx[i]];
        sorted.vectors[i] = out.vectors[idx[i]];
    }
    return sorted;
}

}  // namespace

class Pass13cAxis1bValidation : public ::testing::Test
{
};

// =========================================================================
// CLVDContentWithOverburden: K_0 = 0.5 IC + isotropic compression
// strain rate. After N steps the surface-integral M has measurable
// CLVD content along the vertical (z) axis.
// =========================================================================
TEST_F(Pass13cAxis1bValidation, CLVDContentWithOverburden)
{
    // Architectural reality (see CLAUDE.md / docs/AXIS_1B_DESIGN.md
    // pass-13c row): without 3D Lagrangian advection the geometric
    // CLVD source is absent and the elastic predictor under isotropic
    // compression cannot grow the deviator. The CLVD signature in M
    // therefore comes from the radial-return projection scaling of the
    // pre-existing IC deviator. To exercise that path on the cube_5tet
    // unit fixture, this gate uses a low-cohesion yield surface (so the
    // IC sqrt(J2) sits above yield from step 1) and a shallow source
    // (so the per-cell stress varies measurably across the cube
    // height).
    //
    // Pass-14 will carry this gate to the literature 0.05-0.30 ratio
    // by adding the geometric advective contribution.
    Source3DBallConfig cfg;
    cfg.mesh_path = fixtureBasename("cube_5tet");
    cfg.cavity_radius_m = 0.0;
    cfg.outer_radius_m = std::sqrt(3.0);
    cfg.asymmetric_overburden = true;
    cfg.overburden_K0 = 0.5;
    cfg.host_density_kg_per_m3 = 2700.0;
    cfg.source_depth_m = 100.0;        // shallow: cube spans ~1% depth
    cfg.bulk_modulus_K_pa = 30.0e9;
    cfg.shear_modulus_G_pa = 18.0e9;
    // Pure von Mises (no pressure dependence) with a cohesion well
    // below the IC sqrt(J2) ~ 0.8 MPa. Cells yield from step 1.
    cfg.dp3d_alpha = 0.0;
    cfg.dp3d_k_pa = 1.0e4;             // 10 kPa
    cfg.cv_J_per_kg_K = 1000.0;
    cfg.T_ambient_K = 300.0;
    cfg.medium_label = "GRANITE";

    Source3DBallImpl ball;
    ball.initialize(cfg);
    ASSERT_TRUE(ball.isInitialized());
    const int N = ball.numLocalCells();
    ASSERT_GT(N, 0);

    // Small isotropic compression strain rate: only enough to perturb
    // the IC. The relevant physics is that the radial-return
    // projection scales the pre-existing IC deviator (axis-z CLVD)
    // down to the (low) yield surface, leaving a residual CLVD
    // signature in delta_sigma.
    const double eps_dot = -1.0e-4;
    std::vector<std::array<double, 6>> rates(
        N, {eps_dot, eps_dot, eps_dot, 0.0, 0.0, 0.0});
    ball.setCellStrainRate(rates);

    // One step is enough: subsequent isotropic steps stay elastic
    // (yield surface grows with -I_1 under alpha = 0 + cohesion only,
    // but pure von Mises means yield = k_dp constant; once the
    // projection lands the deviator on the yield surface, further
    // iso compression preserves J2 and the cell stays exactly on the
    // surface).
    const double dt = 1.0e-3;
    ball.step(dt);

    // Recover the moment tensor after N steps.
    std::array<double, 6> M{};
    ball.getMomentTensor(M);

    // Decompose into isotropic + deviator.
    const double trace = M[0] + M[1] + M[2];
    const double M_iso = trace / 3.0;
    std::array<double, 6> M_dev = {
        M[0] - M_iso, M[1] - M_iso, M[2] - M_iso,
        M[3], M[4], M[5]};
    const double dev_norm = std::sqrt(
        M_dev[0] * M_dev[0] + M_dev[1] * M_dev[1] + M_dev[2] * M_dev[2]
        + 2.0 * (M_dev[3] * M_dev[3] + M_dev[4] * M_dev[4]
                  + M_dev[5] * M_dev[5]));
    const double iso_mag = std::abs(M_iso) * std::sqrt(3.0);

    const double clvd_iso_ratio = (iso_mag > 1.0e-12)
        ? dev_norm / iso_mag
        : 0.0;

    // Diagnostic printout to PR-body capture.
    std::cerr << "CLVDContentWithOverburden:\n";
    std::cerr << "  M  = [" << M[0] << ", " << M[1] << ", " << M[2]
              << ", " << M[3] << ", " << M[4] << ", " << M[5]
              << "] N*m\n";
    std::cerr << "  M_iso (scalar) = " << M_iso << " N*m\n";
    std::cerr << "  |M_dev| (Frob) = " << dev_norm << " N*m\n";
    std::cerr << "  |CLVD| / |iso| = " << clvd_iso_ratio << "\n";

    // Gate 1: |CLVD| / |isotropic| > 0.01 (best-effort, stress-asymmetry-only).
    EXPECT_GT(clvd_iso_ratio, 0.01)
        << "CLVD-to-isotropic ratio " << clvd_iso_ratio
        << " is below the 0.01 best-effort threshold; the asymmetric "
           "overburden initial state should leave a measurable "
           "stress-asymmetry-driven CLVD signature in the integrated "
           "stress drop. Pass-14 3D Lagrangian advection will lift "
           "this to the literature 0.05-0.30 range.";

    // Gate 2: CLVD principal axis aligns with vertical (z) within
    // 10 degrees. The "CLVD axis" is the eigenvector with the
    // most-extreme (sign-distinct) eigenvalue: in the canonical CLVD
    // pattern two eigenvalues are equal and one is twice their sign-
    // flipped value.
    Eig3 eig = eigSymm3(M_dev);
    // Identify the most-extreme eigenvalue. For a CLVD pattern
    // {-M, -M, 2M} (or {M, M, -2M}), the unique-sign eigenvalue is
    // either the smallest or largest; the other two are roughly equal.
    int unique_idx = 0;
    {
        const double l_sm = eig.values[0];
        const double l_md = eig.values[1];
        const double l_lg = eig.values[2];
        const double diff_low = std::abs(l_md - l_sm);
        const double diff_high = std::abs(l_lg - l_md);
        // unique = whichever side has the larger spread from the median
        unique_idx = (diff_low > diff_high) ? 0 : 2;
    }
    const auto& axis = eig.vectors[unique_idx];
    const double axis_z = std::abs(axis[2]);
    const double cos_threshold = std::cos(10.0 * M_PI / 180.0);  // ~0.985

    std::cerr << "  CLVD principal axis = (" << axis[0] << ", "
              << axis[1] << ", " << axis[2] << ")\n";
    std::cerr << "  |dot(axis, z)| = " << axis_z
              << " (threshold cos(10 deg) = " << cos_threshold << ")\n";

    EXPECT_GT(axis_z, cos_threshold)
        << "CLVD principal axis should align with vertical within 10 "
           "degrees; got cos(angle) = " << axis_z;
}

// =========================================================================
// MomentTensorIsZeroWithoutSourceForcing: with no strain-rate forcing
// (zero rate) and no overburden IC, M(t) stays zero. Smoke check that
// the surface integral does not introduce spurious moment content.
// =========================================================================
TEST_F(Pass13cAxis1bValidation, MomentTensorIsZeroWithoutSourceForcing)
{
    Source3DBallConfig cfg;
    cfg.mesh_path = fixtureBasename("cube_5tet");
    cfg.cavity_radius_m = 0.0;
    cfg.outer_radius_m = std::sqrt(3.0);
    cfg.asymmetric_overburden = false;
    cfg.host_density_kg_per_m3 = 2700.0;
    cfg.bulk_modulus_K_pa = 30.0e9;
    cfg.shear_modulus_G_pa = 18.0e9;
    cfg.dp3d_alpha = 0.30;
    cfg.dp3d_k_pa = 70.0e6;
    cfg.cv_J_per_kg_K = 1000.0;
    cfg.T_ambient_K = 300.0;
    cfg.medium_label = "GRANITE";

    Source3DBallImpl ball;
    ball.initialize(cfg);
    const int N = ball.numLocalCells();
    ASSERT_GT(N, 0);

    // Zero strain rate => zero stress evolution.
    ball.setCellStrainRate(
        std::vector<std::array<double, 6>>(N, {0, 0, 0, 0, 0, 0}));

    for (int step = 0; step < 10; ++step) {
        ball.step(1.0e-3);
    }

    std::array<double, 6> M{};
    ball.getMomentTensor(M);
    for (int k = 0; k < 6; ++k) {
        EXPECT_NEAR(M[k], 0.0, 1.0e-3)
            << "M[" << k << "] = " << M[k] << " (should be ~0 with no "
               "source forcing).";
    }
}
