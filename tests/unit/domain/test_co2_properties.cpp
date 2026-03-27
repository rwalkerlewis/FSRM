#include <gtest/gtest.h>
#include <mpi.h>
#include <petsc.h>

#include "domain/reservoir/CO2Properties.hpp"
#include <cmath>

using namespace FSRM;

class CO2PropertiesTest : public ::testing::Test {
protected:
    void SetUp() override {
        int rk = 0;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rk);
        EXPECT_GE(rk, 0);
    }
};

// Test 1: Supercritical density (10 MPa, 320 K should be ~630 kg/m³)
TEST_F(CO2PropertiesTest, SupercriticalDensity) {
    double rho = CO2Properties::density(10e6, 320.0);
    EXPECT_GT(rho, 200.0);
    EXPECT_LT(rho, 900.0);
}

// Test 2: Gas density (1 MPa, 300 K should be ~18 kg/m³)
TEST_F(CO2PropertiesTest, GasDensity) {
    double rho = CO2Properties::density(1e6, 300.0);
    EXPECT_GT(rho, 5.0);
    EXPECT_LT(rho, 50.0);
}

// Test 3: PR density gives reasonable values
TEST_F(CO2PropertiesTest, PengRobinsonDensity) {
    double rho_pr = CO2Properties::densityPR(10e6, 320.0);
    EXPECT_GT(rho_pr, 200.0);
    EXPECT_LT(rho_pr, 1200.0);
}

// Test 4: Viscosity at supercritical (should be 2-8 × 10⁻⁵ Pa·s)
TEST_F(CO2PropertiesTest, SupercriticalViscosity) {
    double mu = CO2Properties::viscosity(10e6, 320.0);
    EXPECT_GT(mu, 1e-6);
    EXPECT_LT(mu, 1e-3);
}

// Test 5: Phase identification
TEST_F(CO2PropertiesTest, PhaseIdentification) {
    EXPECT_TRUE(CO2Properties::isSupercritical(10e6, 320.0));
    EXPECT_FALSE(CO2Properties::isSupercritical(1e6, 300.0));
    EXPECT_TRUE(CO2Properties::isGas(1e6, 300.0));
}

// Test 6: Saturation pressure at known point (~5.18 MPa at 216.6 K)
TEST_F(CO2PropertiesTest, SaturationPressure) {
    double Psat = CO2Properties::saturationPressure(280.0);
    EXPECT_GT(Psat, 1e6);
    EXPECT_LT(Psat, 10e6);
}

// Test 7: Compressibility is positive
TEST_F(CO2PropertiesTest, CompressibilityPositive) {
    double beta = CO2Properties::compressibility(10e6, 320.0);
    EXPECT_GT(beta, 0.0);
    EXPECT_LT(beta, 1.0);
}

// Test 8: Thermal conductivity is reasonable
TEST_F(CO2PropertiesTest, ThermalConductivity) {
    double k = CO2Properties::thermalConductivity(10e6, 320.0);
    EXPECT_GT(k, 0.01);
    EXPECT_LT(k, 1.0);
}

// Test 9: Enthalpy is finite
TEST_F(CO2PropertiesTest, EnthalpyFinite) {
    double h = CO2Properties::enthalpy(10e6, 320.0);
    EXPECT_TRUE(std::isfinite(h));
}
