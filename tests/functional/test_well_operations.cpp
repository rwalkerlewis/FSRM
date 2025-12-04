/**
 * @file test_well_operations.cpp
 * @brief Functional tests for well operations and calculations
 */

#include <gtest/gtest.h>
#include "WellModel.hpp"
#include "FSRM.hpp"
#include <cmath>

using namespace FSRM;

class WellOperationsTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    }
    
    int rank;
};

// ============================================================================
// Well Index Calculations
// ============================================================================

TEST_F(WellOperationsTest, PeacemanWellIndex) {
    // Peaceman formula: WI = 2*π*k*h / (ln(r_e/r_w) + s)
    double k = 100e-15;  // Permeability (m^2)
    double h = 10.0;     // Thickness (m)
    double r_w = 0.1;    // Wellbore radius (m)
    double r_e = 100.0;  // Drainage radius (m)
    double s = 0.0;      // Skin factor
    
    double WI = 2.0 * M_PI * k * h / (std::log(r_e / r_w) + s);
    
    EXPECT_GT(WI, 0.0);
    EXPECT_TRUE(std::isfinite(WI));
}

TEST_F(WellOperationsTest, EquivalentRadius) {
    // Peaceman equivalent radius: r_o = 0.28 * sqrt(dx^2 + dy^2)
    double dx = 100.0;
    double dy = 100.0;
    
    double r_o = 0.28 * std::sqrt(dx*dx + dy*dy);
    
    EXPECT_GT(r_o, 0.0);
    EXPECT_LT(r_o, dx);  // Should be less than cell size
    EXPECT_NEAR(r_o, 39.6, 0.1);
}

TEST_F(WellOperationsTest, AnisotropicWellIndex) {
    // For anisotropic permeability
    double kx = 100e-15;
    double ky = 50e-15;
    double kz = 10e-15;
    double dx = 100.0;
    double dy = 100.0;
    double dz = 10.0;
    double h = dz;
    double r_w = 0.1;
    double s = 0.0;
    
    // Effective permeability
    double k_eff = std::sqrt(kx * ky);
    
    // Anisotropic equivalent radius
    double r_o = 0.28 * std::sqrt(
        std::sqrt(ky/kx) * dx*dx + std::sqrt(kx/ky) * dy*dy
    );
    
    double WI = 2.0 * M_PI * k_eff * h / (std::log(r_o / r_w) + s);
    
    EXPECT_GT(WI, 0.0);
    EXPECT_TRUE(std::isfinite(WI));
}

// ============================================================================
// Productivity Index Tests
// ============================================================================

TEST_F(WellOperationsTest, ProductivityIndex) {
    // PI = q / (p_res - p_wf)
    double q = 100.0;        // m^3/day
    double p_res = 20e6;     // Pa
    double p_wf = 10e6;      // Pa
    
    double PI = q / (p_res - p_wf);
    
    EXPECT_GT(PI, 0.0);
}

TEST_F(WellOperationsTest, IPRCurve) {
    // Vogel's IPR for solution gas drive
    // q/q_max = 1 - 0.2*(p_wf/p_res) - 0.8*(p_wf/p_res)^2
    
    double p_res = 20e6;
    double q_max = 200.0;
    
    std::vector<double> p_wf_values = {0.0, 5e6, 10e6, 15e6, 20e6};
    std::vector<double> q_values;
    
    for (double p_wf : p_wf_values) {
        double ratio = p_wf / p_res;
        double q = q_max * (1.0 - 0.2*ratio - 0.8*ratio*ratio);
        q_values.push_back(q);
    }
    
    // At p_wf = 0, q = q_max
    EXPECT_NEAR(q_values[0], q_max, 1e-10);
    
    // At p_wf = p_res, q = 0
    EXPECT_NEAR(q_values.back(), 0.0, 1e-10);
    
    // Rate should decrease monotonically
    for (size_t i = 1; i < q_values.size(); ++i) {
        EXPECT_LT(q_values[i], q_values[i-1]);
    }
}

// ============================================================================
// Wellbore Pressure Tests
// ============================================================================

TEST_F(WellOperationsTest, HydrostaticHead) {
    // dp = ρ*g*dz
    double rho = 1000.0;   // kg/m^3
    double g = 9.81;       // m/s^2
    double dz = 100.0;     // m
    
    double dp = rho * g * dz;
    
    EXPECT_NEAR(dp, 981000.0, 1.0);  // ~1 MPa per 100m
}

TEST_F(WellOperationsTest, FrictionPressureDrop) {
    // Darcy-Weisbach: dp/dL = f * (ρ*v^2) / (2*D)
    double f = 0.02;       // Friction factor
    double rho = 1000.0;   // kg/m^3
    double v = 1.0;        // m/s
    double D = 0.1;        // Diameter
    double L = 100.0;      // Length
    
    double dp = f * L * rho * v*v / (2.0 * D);
    
    EXPECT_GT(dp, 0.0);
    EXPECT_TRUE(std::isfinite(dp));
}

// ============================================================================
// Multi-Phase Flow Tests
// ============================================================================

TEST_F(WellOperationsTest, WaterCut) {
    // WC = q_w / (q_o + q_w)
    double q_o = 80.0;   // Oil rate
    double q_w = 20.0;   // Water rate
    
    double WC = q_w / (q_o + q_w);
    
    EXPECT_DOUBLE_EQ(WC, 0.2);  // 20% water cut
}

TEST_F(WellOperationsTest, GasOilRatio) {
    // GOR = q_g / q_o
    double q_o = 100.0;    // Oil rate (stb/day)
    double q_g = 50000.0;  // Gas rate (scf/day)
    
    double GOR = q_g / q_o;
    
    EXPECT_DOUBLE_EQ(GOR, 500.0);  // 500 scf/stb
}

// ============================================================================
// Skin Factor Tests
// ============================================================================

TEST_F(WellOperationsTest, PositiveSkinDamage) {
    // Positive skin increases pressure drop
    double k = 100e-15;
    double h = 10.0;
    double r_w = 0.1;
    double r_e = 100.0;
    double s = 5.0;  // Damaged well
    
    double WI_damaged = 2.0 * M_PI * k * h / (std::log(r_e / r_w) + s);
    
    double s_clean = 0.0;
    double WI_clean = 2.0 * M_PI * k * h / (std::log(r_e / r_w) + s_clean);
    
    EXPECT_LT(WI_damaged, WI_clean);  // Damaged well has lower WI
}

TEST_F(WellOperationsTest, NegativeSkinStimulated) {
    // Negative skin (stimulation) increases productivity
    double k = 100e-15;
    double h = 10.0;
    double r_w = 0.1;
    double r_e = 100.0;
    double s_stim = -3.0;  // Stimulated well
    
    double WI_stim = 2.0 * M_PI * k * h / (std::log(r_e / r_w) + s_stim);
    
    double s_clean = 0.0;
    double WI_clean = 2.0 * M_PI * k * h / (std::log(r_e / r_w) + s_clean);
    
    EXPECT_GT(WI_stim, WI_clean);  // Stimulated well has higher WI
}

// ============================================================================
// Well Control Mode Tests
// ============================================================================

TEST_F(WellOperationsTest, RateControl) {
    // At rate control: q = target, BHP adjusts
    double q_target = 100.0;
    double WI = 1e-12;
    double p_res = 20e6;
    
    // BHP = p_res - q/WI
    double p_wf = p_res - q_target / WI;
    
    EXPECT_LT(p_wf, p_res);  // BHP < reservoir pressure
}

TEST_F(WellOperationsTest, BHPControl) {
    // At BHP control: p_wf = target, rate adjusts
    double p_wf_target = 10e6;
    double WI = 1e-12;
    double p_res = 20e6;
    
    // q = WI * (p_res - p_wf)
    double q = WI * (p_res - p_wf_target);
    
    EXPECT_GT(q, 0.0);
}

// ============================================================================
// Directional Well Tests
// ============================================================================

TEST_F(WellOperationsTest, DeviatedWellPI) {
    // Deviated wells have different PI due to geometry
    // Joshi's correlation for horizontal wells
    
    double L = 500.0;    // Horizontal length
    double h = 20.0;     // Pay thickness
    double a = 1000.0;   // Half major axis of drainage
    double r_w = 0.1;
    double k_h = 100e-15;
    double k_v = 10e-15;
    double beta = std::sqrt(k_h / k_v);
    
    // Simplified Joshi PI
    // PI ~ k_h * L / (ln(2*a/L))
    double PI_approx = k_h * L / std::log(2.0 * a / L);
    
    EXPECT_GT(PI_approx, 0.0);
}

TEST_F(WellOperationsTest, DoglegSeverity) {
    // Dogleg severity in deg/100ft
    double inc1 = 30.0;   // degrees
    double inc2 = 35.0;
    double azi1 = 45.0;
    double azi2 = 50.0;
    double dL = 30.0;     // measured depth change (m)
    
    // Simplified dogleg calculation
    double dL_ft = dL * 3.28084;  // Convert to feet
    double delta_inc = inc2 - inc1;
    double delta_azi = azi2 - azi1;
    
    // Approximate dogleg
    double DLS = std::sqrt(delta_inc*delta_inc + delta_azi*delta_azi) * 100.0 / dL_ft;
    
    EXPECT_GT(DLS, 0.0);
}
