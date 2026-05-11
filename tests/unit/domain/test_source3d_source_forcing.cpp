/**
 * @file test_source3d_source_forcing.cpp
 * @brief Pass-14a (axis-1b source-forcing slice) unit gates for the
 *        volumetric energy deposition that drives the 3D source ball.
 *
 * Gates:
 *   InnerCavityCellsIdentifiedCorrectly       -- the cavity-wall cells
 *       are tagged (DMLabel SourceBallInnerCavity + the impl's list),
 *       everything outside the cavity-wall layer is not.
 *   SourceTimeFunctionDeliversTotalYield      -- integrating the STF
 *       over the pulse produces yield_kt * 4.184e12 * efficiency to
 *       0.1 percent, for all four STF shapes.
 *   VolumetricEnergyDistributionConservation  -- one host step in the
 *       deposition window deposits exactly rate(t_mid) * dt globally,
 *       distributed uniformly per unit mass over the inner-cavity cells.
 *
 * References: Mueller & Murphy 1971 BSSA 61(6); Brune 1970 JGR 75;
 * Denny & Johnson 1991 (AGU Monogr. 65).
 */

#include <gtest/gtest.h>

#include <petscsys.h>
#include <petscdmplex.h>
#include <petscdmlabel.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <vector>

#include "domain/explosion/Source3DBall.hpp"
#include "domain/explosion/Source3DBallImpl.hpp"
#include "domain/explosion/Source3DBallMesh.hpp"

using FSRM::Source3DBallConfig;
using FSRM::Source3DBallImpl;
using FSRM::Source3DCellState;

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

Source3DBallConfig baseCfgCube()
{
    Source3DBallConfig cfg;
    cfg.mesh_path = fixtureBasename("cube_5tet");
    cfg.cavity_radius_m = 0.0;       // cube fixture passes through origin
    cfg.outer_radius_m = std::sqrt(3.0);
    cfg.asymmetric_overburden = false;  // zero initial stress
    cfg.host_density_kg_per_m3 = 2700.0;
    cfg.bulk_modulus_K_pa = 30.0e9;
    cfg.shear_modulus_G_pa = 18.0e9;
    cfg.dp3d_alpha = 0.30;
    cfg.dp3d_k_pa = 70.0e6;
    cfg.cv_J_per_kg_K = 1000.0;
    cfg.T_ambient_K = 300.0;
    cfg.medium_label = "GRANITE";
    cfg.source_forcing_enabled = true;
    cfg.source_yield_kt = 5.3;
    cfg.source_deposition_duration_s = 1.0e-3;
    cfg.source_deposition_efficiency = 1.0;
    return cfg;
}

double centroidRadius(const std::array<double, 3>& c)
{
    return std::sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
}

}  // namespace

class Source3DSourceForcing : public ::testing::Test {};

// =========================================================================
// InnerCavityCellsIdentifiedCorrectly
// =========================================================================
TEST_F(Source3DSourceForcing, InnerCavityCellsIdentifiedCorrectly)
{
    Source3DBallConfig cfg = baseCfgCube();
    Source3DBallImpl impl;
    impl.initialize(cfg);

    ASSERT_TRUE(impl.sourceForcingActive())
        << "source forcing should activate (positive yield, mesh has a "
           "cavity-wall layer)";
    const int n = impl.numLocalCells();
    ASSERT_GT(n, 0);

    std::vector<Source3DCellState> states;
    impl.getCellStates(states);
    ASSERT_EQ(static_cast<int>(states.size()), n);

    const auto& inner = impl.innerCavityCellIndices();
    ASSERT_FALSE(inner.empty()) << "at least the cavity-wall layer must "
                                   "be tagged";

    // The strict criterion (centroid <= cavity_radius_m = 0) tags no
    // cell, so the impl falls back to the innermost layer: the largest
    // centroid radius among tagged cells must be strictly below the
    // smallest centroid radius among untagged cells.
    std::vector<unsigned char> tagged(n, 0);
    double max_tagged_r = 0.0;
    for (int i : inner) {
        ASSERT_GE(i, 0);
        ASSERT_LT(i, n);
        tagged[i] = 1;
        max_tagged_r = std::max(max_tagged_r, centroidRadius(states[i].centroid));
        EXPECT_DOUBLE_EQ(states[i].inner_cavity_marker, 1.0);
    }
    double min_untagged_r = 1.0e300;
    int n_untagged = 0;
    for (int i = 0; i < n; ++i) {
        if (tagged[i]) continue;
        ++n_untagged;
        min_untagged_r = std::min(min_untagged_r, centroidRadius(states[i].centroid));
        EXPECT_DOUBLE_EQ(states[i].inner_cavity_marker, 0.0);
    }
    EXPECT_GT(n_untagged, 0)
        << "cube_5tet has one corner-far cell that should stay untagged";
    EXPECT_LT(max_tagged_r, min_untagged_r)
        << "tagged (inner-cavity) cells must be radially inside the "
           "untagged ones";

    // The DMLabel must agree with the impl's list.
    DM dm = impl.mesh() ? impl.mesh()->getDM() : nullptr;
    ASSERT_NE(dm, nullptr);
    DMLabel cav_label = nullptr;
    DMGetLabel(dm, "SourceBallInnerCavity", &cav_label);
    ASSERT_NE(cav_label, nullptr) << "SourceBallInnerCavity label must "
                                     "be created at initialize()";
    PetscInt label_size = 0;
    DMLabelGetStratumSize(cav_label, 1, &label_size);
    EXPECT_EQ(static_cast<int>(label_size), impl.numInnerCavityCells());
}

// =========================================================================
// SourceTimeFunctionDeliversTotalYield
// =========================================================================
TEST_F(Source3DSourceForcing, SourceTimeFunctionDeliversTotalYield)
{
    const std::array<Source3DBallConfig::SourceTimeFunction, 4> shapes = {
        Source3DBallConfig::SourceTimeFunction::MUELLER_MURPHY,
        Source3DBallConfig::SourceTimeFunction::BRUNE,
        Source3DBallConfig::SourceTimeFunction::RAMP,
        Source3DBallConfig::SourceTimeFunction::DELTA,
    };
    for (auto shape : shapes) {
        Source3DBallConfig cfg = baseCfgCube();
        cfg.source_time_function = shape;
        Source3DBallImpl impl;
        impl.initialize(cfg);
        ASSERT_TRUE(impl.sourceForcingActive());

        const double E_total = 5.3 * kJoulesPerKt;  // efficiency 1.0
        EXPECT_NEAR(impl.totalSourceEnergyJ(), E_total, 1.0e-6 * E_total);

        // Midpoint-rule integral of the STF rate over a long window
        // (50 deposition durations -- effectively infinity for the
        // exponentially-decaying shapes; exact for the boxcar RAMP /
        // DELTA shapes since the midpoint rule integrates piecewise-
        // constant functions exactly).
        const double T = cfg.source_deposition_duration_s;
        const double t_end = 50.0 * T;
        const int N = 400000;
        const double dt = t_end / N;
        double integral = 0.0;
        for (int k = 0; k < N; ++k) {
            integral += impl.sourceTimeFunctionPowerJ((k + 0.5) * dt) * dt;
        }
        EXPECT_NEAR(integral, E_total, 1.0e-3 * E_total)
            << "STF shape " << static_cast<int>(shape)
            << " must integrate to the total deposited energy";
    }
}

// =========================================================================
// VolumetricEnergyDistributionConservation
// =========================================================================
TEST_F(Source3DSourceForcing, VolumetricEnergyDistributionConservation)
{
    Source3DBallConfig cfg = baseCfgCube();
    cfg.source_time_function = Source3DBallConfig::SourceTimeFunction::RAMP;
    Source3DBallImpl impl;
    impl.initialize(cfg);
    ASSERT_TRUE(impl.sourceForcingActive());

    // One host step squarely inside the deposition window.
    const double dt = 0.10 * cfg.source_deposition_duration_s;
    impl.step(dt);

    // The deposited energy must equal rate(t_mid) * dt (the impl uses
    // the midpoint rate of the step).
    const double rate_mid = impl.sourceTimeFunctionPowerJ(0.5 * dt);
    const double expected_dE = rate_mid * dt;
    EXPECT_GT(expected_dE, 0.0);
    EXPECT_NEAR(impl.depositedEnergyGlobalJ(), expected_dE,
                1.0e-9 * expected_dE)
        << "global deposited energy must match rate(t_mid) * dt";

    // Each inner-cavity cell carries the same instantaneous deposited
    // power density rate / V_inner_global (uniform per unit volume).
    const double V_inner = impl.innerCavityVolumeGlobal();
    ASSERT_GT(V_inner, 0.0);
    const double expected_power_density = rate_mid / V_inner;

    std::vector<Source3DCellState> states;
    impl.getCellStates(states);
    const auto& inner = impl.innerCavityCellIndices();
    ASSERT_FALSE(inner.empty());
    double power_sum_vol_weighted = 0.0;  // sum power_density * V approx
    for (int i : inner) {
        EXPECT_NEAR(states[i].source_forcing_power, expected_power_density,
                    1.0e-6 * expected_power_density + 1.0e-30)
            << "inner-cavity cell " << i
            << " power density must be rate / V_inner";
    }
    (void)power_sum_vol_weighted;

    // Cells outside the cavity wall carry zero deposited power.
    std::vector<unsigned char> tagged(impl.numLocalCells(), 0);
    for (int i : inner) tagged[i] = 1;
    for (int i = 0; i < impl.numLocalCells(); ++i) {
        if (tagged[i]) continue;
        EXPECT_DOUBLE_EQ(states[i].source_forcing_power, 0.0);
    }
}
