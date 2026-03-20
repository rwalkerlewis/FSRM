/**
 * @file test_unit_system.cpp
 * @brief Smoke tests for UnitSystem conversions
 */

#include "util/UnitSystem.hpp"
#include <cmath>
#include <gtest/gtest.h>

using namespace FSRM;

class UnitSystemTest : public ::testing::Test {
protected:
    UnitSystem units;
};

TEST_F(UnitSystemTest, PressurePaPsiRoundTrip) {
    const double Pa = 101325.0;
    double psi = units.convert(Pa, "Pa", "psi");
    double back = units.convert(psi, "psi", "Pa");
    EXPECT_NEAR(back, Pa, 1.0);
}

TEST_F(UnitSystemTest, LengthMetersFeetRoundTrip) {
    const double m = 10.0;
    double ft = units.convert(m, "m", "ft");
    double m2 = units.convert(ft, "ft", "m");
    EXPECT_NEAR(m2, m, 1e-9);
}

TEST_F(UnitSystemTest, PermeabilityM2MillidarcyRoundTrip) {
    const double k_si = 1e-15;
    double md = units.convert(k_si, "m2", "mD");
    double k2 = units.convert(md, "mD", "m2");
    EXPECT_NEAR(k2, k_si, 1e-30);
}

TEST_F(UnitSystemTest, TemperatureK_C_F_RoundTrip) {
    const double T_K = 373.15;
    double C = units.convert(T_K, "K", "degC");
    double K_from_C = units.convert(C, "degC", "K");
    EXPECT_NEAR(K_from_C, T_K, 1e-6);

    double F = units.convert(T_K, "K", "degF");
    double K_from_F = units.convert(F, "degF", "K");
    EXPECT_NEAR(K_from_F, T_K, 1e-4);
}

TEST_F(UnitSystemTest, CompoundStressFromPsi) {
    double sigma = units.convert(1000.0, "psi", "Pa");
    EXPECT_GT(sigma, 0.0);
    double back = units.convert(sigma, "Pa", "psi");
    EXPECT_NEAR(back, 1000.0, 1e-6);
}

TEST_F(UnitSystemTest, SIToFieldDensityStyle) {
    const double rho_si = 1000.0;
    double lb_ft3 = units.convert(rho_si, "kg/m3", "lbm/ft3");
    EXPECT_GT(lb_ft3, 0.0);
    double rho2 = units.convert(lb_ft3, "lbm/ft3", "kg/m3");
    EXPECT_NEAR(rho2, rho_si, 0.01);
}
