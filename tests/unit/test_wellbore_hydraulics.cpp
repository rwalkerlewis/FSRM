/**
 * @file test_wellbore_hydraulics.cpp
 * @brief Unit tests for wellbore hydraulics (THP/BHP conversion and pressure drop components)
 */
 
#include <gtest/gtest.h>
#include "WellboreHydraulics.hpp"
#include "WellModel.hpp"
#include <cmath>
 
using namespace FSRM;
 
TEST(WellboreHydraulicsTest, HydrostaticOnlyVertical) {
    WellboreHydraulicsCalculator calc;
 
    std::vector<double> md = {0.0, 1000.0};
    std::vector<double> inc = {0.0, 0.0};
    std::vector<double> azi = {0.0, 0.0};
    std::vector<double> diam = {0.1, 0.1};
    calc.buildFromSurvey(md, inc, azi, diam, 1e-5);
 
    const double thp = 1.0e6;
    const double q = 0.0;
    const double rho = 1000.0;
    const double mu = 0.001;
 
    const double bhp = calc.calculateBHP(thp, q, rho, mu);
    const double expected = thp + rho * 9.80665 * 1000.0;
    EXPECT_NEAR(bhp, expected, 1e-3 * expected);
}
 
TEST(WellboreHydraulicsTest, FrictionAndDragReduction) {
    WellboreHydraulicsCalculator calc;
 
    std::vector<double> md = {0.0, 1000.0};
    std::vector<double> inc = {0.0, 0.0};
    std::vector<double> azi = {0.0, 0.0};
    std::vector<double> diam = {0.1, 0.1};
    calc.buildFromSurvey(md, inc, azi, diam, 1e-5);
 
    const double thp = 1.0e6;
    const double q = 0.01;  // m3/s
    const double rho = 1000.0;
    const double mu = 0.001;
 
    (void)calc.calculateBHP(thp, q, rho, mu);
    auto base = calc.getPressureBreakdown();
    EXPECT_GT(base["friction"], 0.0);
 
    calc.setDragReduction(0.5);
    (void)calc.calculateBHP(thp, q, rho, mu);
    auto dr = calc.getPressureBreakdown();
    EXPECT_LT(dr["friction"], base["friction"]);
}
 
TEST(WellboreHydraulicsTest, WellModelTHPtoBHPConsistency) {
    ProductionWell well("PROD1");
 
    WellModel::WellboreHydraulicsOptions opts;
    opts.measured_depth = 1000.0;
    opts.tubing_id = 0.1;
    opts.roughness = 1e-5;
    opts.drag_reduction_factor = 0.0;
    opts.fluid_density = 1000.0;
    opts.fluid_viscosity = 0.001;
    opts.friction_correlation = FrictionModel::Correlation::HAALAND;
    well.configureWellboreHydraulics(opts);
 
    const double thp = 1.0e6;
    const double rate = 0.01;
    const double bhp = well.computeBHPFromTHP(thp, rate);
 
    // Basic sanity: BHP must exceed THP by at least hydrostatic head.
    EXPECT_GT(bhp, thp + opts.fluid_density * 9.80665 * 900.0);
 
    const double thp_back = well.computeTHPFromBHP(bhp, rate);
    EXPECT_NEAR(thp_back, thp, 1e-3 * thp);
}

