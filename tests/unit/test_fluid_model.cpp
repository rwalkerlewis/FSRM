/**
 * @file test_fluid_model.cpp
 * @brief Unit tests for FluidModel classes
 */

#include <gtest/gtest.h>
#include "FluidModel.hpp"
#include <cmath>

using namespace FSRM;

class FluidModelTest : public ::testing::Test {
protected:
    void SetUp() override {}
};

// ============================================================================
// FluidType Parsing Tests
// ============================================================================

TEST_F(FluidModelTest, ParseFluidTypeSinglePhase) {
    EXPECT_EQ(parseFluidType("SINGLE_PHASE"), FluidType::SINGLE_PHASE);
    EXPECT_EQ(parseFluidType("single_phase"), FluidType::SINGLE_PHASE);
}

TEST_F(FluidModelTest, ParseFluidTypeBlackOil) {
    EXPECT_EQ(parseFluidType("BLACK_OIL"), FluidType::BLACK_OIL);
    EXPECT_EQ(parseFluidType("black_oil"), FluidType::BLACK_OIL);
}

TEST_F(FluidModelTest, ParseFluidTypeCompositional) {
    EXPECT_EQ(parseFluidType("COMPOSITIONAL"), FluidType::COMPOSITIONAL);
}

// ============================================================================
// SinglePhaseFluid Tests
// ============================================================================

TEST_F(FluidModelTest, SinglePhaseFluidCreate) {
    SinglePhaseFluid fluid;
    EXPECT_EQ(fluid.getType(), FluidType::SINGLE_PHASE);
}

TEST_F(FluidModelTest, SinglePhaseFluidDensity) {
    SinglePhaseFluid fluid;
    std::map<std::string, std::string> config;
    config["density"] = "1000.0";
    config["compressibility"] = "4.5e-10";
    config["reference_pressure"] = "1e5";
    fluid.configure(config);
    
    // At reference pressure, should return reference density
    double rho = fluid.getDensity(1e5, 300.0);
    EXPECT_NEAR(rho, 1000.0, 0.1);
    
    // At higher pressure, density should be higher (compressible)
    double rho_high = fluid.getDensity(20e6, 300.0);
    EXPECT_GT(rho_high, rho);
}

TEST_F(FluidModelTest, SinglePhaseFluidViscosity) {
    SinglePhaseFluid fluid;
    std::map<std::string, std::string> config;
    config["viscosity"] = "0.001";
    fluid.configure(config);
    
    double mu = fluid.getViscosity(10e6, 300.0);
    EXPECT_NEAR(mu, 0.001, 1e-6);
}

TEST_F(FluidModelTest, SinglePhaseFluidCompressibility) {
    SinglePhaseFluid fluid;
    std::map<std::string, std::string> config;
    config["compressibility"] = "4.5e-10";
    fluid.configure(config);
    
    double c = fluid.getCompressibility(10e6, 300.0);
    EXPECT_NEAR(c, 4.5e-10, 1e-12);
}

// ============================================================================
// BlackOilFluid Tests
// ============================================================================

TEST_F(FluidModelTest, BlackOilFluidCreate) {
    BlackOilFluid fluid;
    EXPECT_EQ(fluid.getType(), FluidType::BLACK_OIL);
}

TEST_F(FluidModelTest, BlackOilFluidDensity) {
    BlackOilFluid fluid;
    std::map<std::string, std::string> config;
    config["oil_density_std"] = "850.0";
    config["solution_gor"] = "100.0";
    config["pvt_correlation"] = "STANDING";
    fluid.configure(config);
    
    double rho = fluid.getDensity(20e6, 350.0);
    EXPECT_GT(rho, 0.0);
    EXPECT_TRUE(std::isfinite(rho));
}

TEST_F(FluidModelTest, BlackOilFluidViscosity) {
    BlackOilFluid fluid;
    std::map<std::string, std::string> config;
    config["oil_viscosity_std"] = "0.002";
    fluid.configure(config);
    
    double mu = fluid.getViscosity(20e6, 350.0);
    EXPECT_GT(mu, 0.0);
    EXPECT_TRUE(std::isfinite(mu));
}

TEST_F(FluidModelTest, BlackOilFluidSolutionGOR) {
    BlackOilFluid fluid;
    std::map<std::string, std::string> config;
    config["solution_gor"] = "100.0";
    config["bubble_point"] = "15e6";
    fluid.configure(config);
    
    // Below bubble point, GOR should be less than max
    double Rs_below = fluid.getSolutionGOR(10e6);
    double Rs_at = fluid.getSolutionGOR(15e6);
    
    EXPECT_LE(Rs_below, Rs_at);
    EXPECT_GT(Rs_at, 0.0);
}

// ============================================================================
// CompositionalFluid Tests
// ============================================================================

TEST_F(FluidModelTest, CompositionalFluidCreate) {
    CompositionalFluid fluid;
    EXPECT_EQ(fluid.getType(), FluidType::COMPOSITIONAL);
}

TEST_F(FluidModelTest, CompositionalFluidEOS) {
    CompositionalFluid fluid;
    std::map<std::string, std::string> config;
    config["num_components"] = "2";
    config["eos_type"] = "PENG_ROBINSON";
    fluid.configure(config);
    
    // Verify EOS is set
    double rho = fluid.getDensity(10e6, 350.0);
    EXPECT_TRUE(std::isfinite(rho));
}

// ============================================================================
// BrineFluid Tests
// ============================================================================

TEST_F(FluidModelTest, BrineFluidCreate) {
    BrineFluid fluid;
    EXPECT_EQ(fluid.getType(), FluidType::BRINE);
}

TEST_F(FluidModelTest, BrineFluidSalinityEffect) {
    BrineFluid fluid;
    std::map<std::string, std::string> config;
    config["salinity"] = "100000";  // 100,000 ppm
    fluid.configure(config);
    
    // Brine should be denser than freshwater
    double rho_brine = fluid.getDensity(10e6, 300.0);
    EXPECT_GT(rho_brine, 1000.0);
}

// ============================================================================
// CO2Fluid Tests
// ============================================================================

TEST_F(FluidModelTest, CO2FluidCreate) {
    CO2Fluid fluid;
    EXPECT_EQ(fluid.getType(), FluidType::CO2);
}

TEST_F(FluidModelTest, CO2FluidSupercritical) {
    CO2Fluid fluid;
    
    // Above critical point (73.8 bar, 31.1°C)
    double rho = fluid.getDensity(10e6, 320.0);  // ~100 bar, 47°C
    EXPECT_GT(rho, 0.0);
    EXPECT_TRUE(std::isfinite(rho));
}

// ============================================================================
// Factory Function Tests
// ============================================================================

TEST_F(FluidModelTest, CreateFluidModelSinglePhase) {
    std::map<std::string, std::string> config;
    auto fluid = createFluidModel("SINGLE_PHASE", config);
    EXPECT_NE(fluid, nullptr);
    EXPECT_EQ(fluid->getType(), FluidType::SINGLE_PHASE);
}

TEST_F(FluidModelTest, CreateFluidModelBlackOil) {
    std::map<std::string, std::string> config;
    auto fluid = createFluidModel("BLACK_OIL", config);
    EXPECT_NE(fluid, nullptr);
    EXPECT_EQ(fluid->getType(), FluidType::BLACK_OIL);
}

TEST_F(FluidModelTest, CreateFluidModelWithConfig) {
    std::map<std::string, std::string> config;
    config["density"] = "850.0";
    config["viscosity"] = "0.002";
    
    auto fluid = createFluidModel("SINGLE_PHASE", config);
    EXPECT_NE(fluid, nullptr);
    
    double rho = fluid->getDensity(1e5, 300.0);
    EXPECT_NEAR(rho, 850.0, 10.0);
}

// ============================================================================
// Physical Validity Tests
// ============================================================================

TEST_F(FluidModelTest, DensityPositive) {
    SinglePhaseFluid fluid;
    std::map<std::string, std::string> config;
    config["density"] = "1000.0";
    fluid.configure(config);
    
    // Density must always be positive
    for (double p = 1e5; p <= 50e6; p += 5e6) {
        for (double T = 280.0; T <= 400.0; T += 20.0) {
            double rho = fluid.getDensity(p, T);
            EXPECT_GT(rho, 0.0) << "P=" << p << " T=" << T;
        }
    }
}

TEST_F(FluidModelTest, ViscosityPositive) {
    SinglePhaseFluid fluid;
    std::map<std::string, std::string> config;
    config["viscosity"] = "0.001";
    fluid.configure(config);
    
    // Viscosity must always be positive
    for (double p = 1e5; p <= 50e6; p += 5e6) {
        for (double T = 280.0; T <= 400.0; T += 20.0) {
            double mu = fluid.getViscosity(p, T);
            EXPECT_GT(mu, 0.0) << "P=" << p << " T=" << T;
        }
    }
}

TEST_F(FluidModelTest, CompressibilityPositive) {
    SinglePhaseFluid fluid;
    std::map<std::string, std::string> config;
    config["compressibility"] = "4.5e-10";
    fluid.configure(config);
    
    // Compressibility must be positive
    double c = fluid.getCompressibility(10e6, 300.0);
    EXPECT_GT(c, 0.0);
}
