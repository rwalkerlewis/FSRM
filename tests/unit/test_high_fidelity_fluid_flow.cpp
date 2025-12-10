/**
 * @file test_high_fidelity_fluid_flow.cpp
 * @brief Unit tests for HighFidelityFluidFlow module
 * 
 * Tests cover:
 * - Non-Darcy flow (Forchheimer, Klinkenberg)
 * - Dual/Triple porosity-permeability models
 * - Non-isothermal multiphase flow
 * - Dynamic relative permeability
 * - Miscible flow mixing models
 */

#include <gtest/gtest.h>
#include "HighFidelityFluidFlow.hpp"
#include <cmath>

using namespace FSRM::HighFidelity;

// =============================================================================
// Non-Darcy Flow Tests
// =============================================================================

class NonDarcyFlowTest : public ::testing::Test {
protected:
    NonDarcyFlow ndf;
    
    void SetUp() override {
        NonDarcyFlow::Parameters params;
        params.model = NonDarcyModel::FORCHHEIMER;
        params.forchheimer_beta = 1.0e8;  // Typical value for sandstone
        params.klinkenberg_b = 1.0e5;     // Pa
        params.ergun_a = 150.0;
        params.ergun_b = 1.75;
        ndf.setParameters(params);
    }
};

TEST_F(NonDarcyFlowTest, ForchheimerCorrection) {
    // Darcy velocity: 0.01 m/s (high for porous media)
    double v_darcy = 0.01;
    double rho = 1000.0;   // kg/m³ (water)
    double mu = 0.001;     // Pa·s (water)
    double K = 1e-12;      // m² (permeability)
    
    double correction = ndf.calculateCorrection(v_darcy, rho, mu, K);
    
    // Forchheimer correction reduces effective permeability
    EXPECT_GT(correction, 0.0);
    EXPECT_LT(correction, 1.0);
}

TEST_F(NonDarcyFlowTest, ForchheimerLowVelocity) {
    // At low velocity, correction should be close to 1 (Darcy regime)
    double v_darcy = 1e-6;  // Very low velocity
    double rho = 1000.0;
    double mu = 0.001;
    double K = 1e-12;
    
    double correction = ndf.calculateCorrection(v_darcy, rho, mu, K);
    
    // Should be close to 1.0 at low velocity
    EXPECT_NEAR(correction, 1.0, 0.01);
}

TEST_F(NonDarcyFlowTest, KlinkenbergCorrection) {
    NonDarcyFlow::Parameters params;
    params.model = NonDarcyModel::KLINKENBERG;
    params.klinkenberg_b = 1.0e5;  // Pa
    ndf.setParameters(params);
    
    double K_abs = 1e-15;  // m² (tight formation)
    double pressure = 1e6;  // Pa (1 MPa)
    
    double K_app = ndf.klinkenbergPermeability(K_abs, pressure);
    
    // Apparent permeability should be higher than absolute
    EXPECT_GT(K_app, K_abs);
}

TEST_F(NonDarcyFlowTest, KlinkenbergHighPressure) {
    NonDarcyFlow::Parameters params;
    params.model = NonDarcyModel::KLINKENBERG;
    params.klinkenberg_b = 1.0e5;  // Pa
    ndf.setParameters(params);
    
    double K_abs = 1e-15;
    double pressure = 1e8;  // High pressure (100 MPa)
    
    double K_app = ndf.klinkenbergPermeability(K_abs, pressure);
    
    // At high pressure, Klinkenberg effect diminishes: K_app → K_abs
    EXPECT_NEAR(K_app / K_abs, 1.0, 0.01);
}

TEST_F(NonDarcyFlowTest, ReynoldsNumber) {
    double Re = ndf.reynoldsNumber(0.01, 1000.0, 0.001, 1e-12);
    
    // Re = ρ·v·√K / μ = 1000 · 0.01 · 1e-6 / 0.001 = 0.01
    // For porous media, typically Re < 1 is Darcy regime
    EXPECT_NEAR(Re, 0.01, 0.001);
}

TEST_F(NonDarcyFlowTest, EffectiveVelocity) {
    double grad_p = 1e6;  // Pa/m
    double rho = 1000.0;
    double mu = 0.001;
    double K = 1e-12;
    
    double v_eff = ndf.effectiveVelocity(grad_p, rho, mu, K);
    
    // Should be positive (flow in pressure gradient direction)
    EXPECT_GT(v_eff, 0.0);
    
    // Should be less than Darcy velocity (Forchheimer reduces it)
    double v_darcy = K * grad_p / mu;
    EXPECT_LT(v_eff, v_darcy);
}

TEST_F(NonDarcyFlowTest, Configure) {
    std::map<std::string, std::string> config;
    config["non_darcy_model"] = "ergun";
    config["non_darcy_ergun_a"] = "200.0";
    config["non_darcy_ergun_b"] = "2.0";
    
    NonDarcyFlow ndf2;
    ndf2.configure(config);
    
    // Verify it doesn't crash
    double correction = ndf2.calculateCorrection(0.01, 1000.0, 0.001, 1e-12);
    EXPECT_GT(correction, 0.0);
}


// =============================================================================
// Dual Porosity-Permeability Tests
// =============================================================================

class DualPorosityTest : public ::testing::Test {
protected:
    DualPorosityPermeability dpp;
    
    void SetUp() override {
        DualPorosityPermeability::Parameters params;
        params.model = DualPorosityModel::WARREN_ROOT;
        params.matrix_porosity = 0.1;
        params.fracture_porosity = 0.01;
        params.matrix_permeability = 1e-15;  // m²
        params.fracture_permeability = 1e-12;  // m²
        params.shape_factor = 1.0;  // m⁻²
        params.transfer_coeff = 1e-8;  // s⁻¹
        dpp.setParameters(params);
    }
};

TEST_F(DualPorosityTest, TransferRate) {
    double p_matrix = 10e6;    // Pa
    double p_fracture = 5e6;   // Pa (lower pressure in fracture)
    
    double rate = dpp.calculateTransferRate(p_matrix, p_fracture);
    
    // Flow from matrix to fracture (positive rate)
    EXPECT_GT(rate, 0.0);
}

TEST_F(DualPorosityTest, TransferRateEquilibrium) {
    double p_matrix = 10e6;
    double p_fracture = 10e6;  // Same pressure
    
    double rate = dpp.calculateTransferRate(p_matrix, p_fracture);
    
    // No transfer at equilibrium
    EXPECT_NEAR(rate, 0.0, 1e-20);
}

TEST_F(DualPorosityTest, EffectivePorosity) {
    double phi_eff = dpp.effectivePorosity();
    
    // φ_eff = φ_m + φ_f = 0.1 + 0.01 = 0.11
    EXPECT_NEAR(phi_eff, 0.11, 1e-10);
}

TEST_F(DualPorosityTest, EffectivePermeability) {
    double K_eff = dpp.effectivePermeability();
    
    // For fractured system, K_eff dominated by fracture permeability
    EXPECT_GT(K_eff, 1e-15);  // Higher than matrix
    EXPECT_LE(K_eff, 1e-12);  // Not higher than fracture
}

TEST_F(DualPorosityTest, StorageCapacity) {
    double c_matrix = 1e-9;   // Pa⁻¹
    double c_fracture = 1e-10;  // Pa⁻¹
    
    double omega = dpp.storativityRatio(c_matrix, c_fracture);
    
    // ω = (φ·c)_f / [(φ·c)_m + (φ·c)_f]
    // Should be between 0 and 1
    EXPECT_GT(omega, 0.0);
    EXPECT_LT(omega, 1.0);
}

TEST_F(DualPorosityTest, InterporosityFlowCoeff) {
    double lambda = dpp.interporosityFlowCoefficient(100.0);  // L = 100 m characteristic length
    
    // λ = σ · K_m · L² / K_f
    // Should be positive
    EXPECT_GT(lambda, 0.0);
}

TEST_F(DualPorosityTest, TriplePorosityModel) {
    DualPorosityPermeability::Parameters params;
    params.model = DualPorosityModel::TRIPLE_POROSITY;
    params.matrix_porosity = 0.15;
    params.fracture_porosity = 0.02;
    params.vug_porosity = 0.05;
    params.matrix_permeability = 1e-16;
    params.fracture_permeability = 1e-13;
    params.vug_permeability = 1e-11;
    params.matrix_fracture_transfer = 1e-8;
    params.fracture_vug_transfer = 1e-6;
    dpp.setParameters(params);
    
    double phi_eff = dpp.effectivePorosity();
    
    // φ_eff = φ_m + φ_f + φ_v = 0.15 + 0.02 + 0.05 = 0.22
    EXPECT_NEAR(phi_eff, 0.22, 1e-10);
}

TEST_F(DualPorosityTest, Configure) {
    std::map<std::string, std::string> config;
    config["dual_porosity_model"] = "kazemi";
    config["matrix_porosity"] = "0.2";
    config["fracture_porosity"] = "0.05";
    
    DualPorosityPermeability dpp2;
    dpp2.configure(config);
    
    SUCCEED();
}


// =============================================================================
// Non-Isothermal Multiphase Flow Tests
// =============================================================================

class NonIsothermalFlowTest : public ::testing::Test {
protected:
    NonIsothermalMultiphaseFlow nimf;
    
    void SetUp() override {
        NonIsothermalMultiphaseFlow::Parameters params;
        params.reference_temperature = 293.15;  // K (20°C)
        params.thermal_conductivity_solid = 2.5;  // W/(m·K)
        params.thermal_conductivity_water = 0.6;
        params.thermal_conductivity_oil = 0.15;
        params.thermal_conductivity_gas = 0.03;
        params.heat_capacity_solid = 800.0;     // J/(kg·K)
        params.heat_capacity_water = 4186.0;
        params.heat_capacity_oil = 2000.0;
        params.heat_capacity_gas = 1000.0;
        params.joule_thomson_water = -2.2e-7;   // K/Pa
        params.joule_thomson_oil = -1e-6;
        params.joule_thomson_gas = 3e-6;
        params.enable_viscous_heating = true;
        params.enable_joule_thomson = true;
        nimf.setParameters(params);
    }
};

TEST_F(NonIsothermalFlowTest, EffectiveThermalConductivity) {
    double porosity = 0.2;
    double S_w = 0.5;  // Water saturation
    double S_o = 0.3;  // Oil saturation
    double S_g = 0.2;  // Gas saturation
    
    double k_eff = nimf.effectiveThermalConductivity(porosity, S_w, S_o, S_g);
    
    // Should be between min and max component values
    EXPECT_GT(k_eff, 0.03);   // > gas conductivity
    EXPECT_LT(k_eff, 2.5);    // < solid conductivity
}

TEST_F(NonIsothermalFlowTest, EffectiveHeatCapacity) {
    double porosity = 0.2;
    double S_w = 0.5;
    double S_o = 0.3;
    double S_g = 0.2;
    double rho_s = 2650.0;  // kg/m³ (sandstone)
    double rho_w = 1000.0;
    double rho_o = 800.0;
    double rho_g = 100.0;
    
    double cp_eff = nimf.effectiveHeatCapacity(porosity, S_w, S_o, S_g,
                                               rho_s, rho_w, rho_o, rho_g);
    
    // Should be positive
    EXPECT_GT(cp_eff, 0.0);
}

TEST_F(NonIsothermalFlowTest, ViscosityTemperatureDependence) {
    double mu_ref = 0.001;  // Pa·s at 20°C
    double T_ref = 293.15;  // K
    double T = 373.15;      // K (100°C)
    double activation_energy = 15000.0;  // J/mol (water)
    
    double mu_T = nimf.viscosityAtTemperature(mu_ref, T_ref, T, activation_energy);
    
    // Viscosity decreases with temperature
    EXPECT_LT(mu_T, mu_ref);
}

TEST_F(NonIsothermalFlowTest, ViscosityLowTemperature) {
    double mu_ref = 0.001;
    double T_ref = 293.15;
    double T = 273.15;  // K (0°C)
    double activation_energy = 15000.0;
    
    double mu_T = nimf.viscosityAtTemperature(mu_ref, T_ref, T, activation_energy);
    
    // Viscosity increases at lower temperature
    EXPECT_GT(mu_T, mu_ref);
}

TEST_F(NonIsothermalFlowTest, ViscousHeating) {
    double mu = 0.001;      // Pa·s
    double grad_p[3] = {1e5, 0.0, 0.0};  // Pa/m
    double K = 1e-12;       // m²
    
    double Q_vis = nimf.viscousHeatingRate(mu, grad_p, K);
    
    // Should be positive (dissipation)
    EXPECT_GT(Q_vis, 0.0);
}

TEST_F(NonIsothermalFlowTest, JouleThomsonCooling) {
    double JT_coeff = 3e-6;  // K/Pa (positive for gas, cooling on expansion)
    double dp = -1e6;        // Pa (pressure drop)
    
    double dT = nimf.jouleThomsonTemperatureChange(JT_coeff, dp);
    
    // Positive JT + negative dp = negative dT (cooling)
    EXPECT_LT(dT, 0.0);
}

TEST_F(NonIsothermalFlowTest, Configure) {
    std::map<std::string, std::string> config;
    config["thermal_cond_solid"] = "3.0";
    config["thermal_cond_water"] = "0.65";
    config["enable_viscous_heating"] = "false";
    config["enable_joule_thomson"] = "true";
    
    NonIsothermalMultiphaseFlow nimf2;
    nimf2.configure(config);
    
    SUCCEED();
}


// =============================================================================
// Dynamic Relative Permeability Tests
// =============================================================================

class DynamicRelPermTest : public ::testing::Test {
protected:
    DynamicRelativePermeability drp;
    
    void SetUp() override {
        DynamicRelativePermeability::Parameters params;
        params.model = DynamicRelPermModel::CAPILLARY_NUMBER;
        params.reference_capillary_number = 1e-5;
        params.desaturation_exponent = 0.5;
        params.kr_endpoint_residual = 0.2;
        params.kr_endpoint_max = 0.95;
        drp.setParameters(params);
    }
};

TEST_F(DynamicRelPermTest, CapillaryNumber) {
    double mu = 0.001;   // Pa·s
    double v = 1e-5;     // m/s
    double sigma = 0.03; // N/m (interfacial tension)
    
    double Nc = drp.capillaryNumber(mu, v, sigma);
    
    // Nc = μ·v / σ = 0.001 · 1e-5 / 0.03 ≈ 3.3e-7
    EXPECT_NEAR(Nc, 3.33e-7, 1e-8);
}

TEST_F(DynamicRelPermTest, DesaturationCoefficient) {
    double Nc_high = 1e-3;   // High capillary number
    double Nc_low = 1e-7;    // Low capillary number
    
    double D_high = drp.desaturationCoefficient(Nc_high);
    double D_low = drp.desaturationCoefficient(Nc_low);
    
    // Higher capillary number -> more desaturation
    EXPECT_GT(D_high, D_low);
    
    // Coefficients should be bounded [0, 1]
    EXPECT_GE(D_high, 0.0);
    EXPECT_LE(D_high, 1.0);
    EXPECT_GE(D_low, 0.0);
    EXPECT_LE(D_low, 1.0);
}

TEST_F(DynamicRelPermTest, DynamicResidualSaturation) {
    double Sor_static = 0.2;  // Static residual oil saturation
    double D = 0.5;           // Desaturation coefficient
    
    double Sor_dynamic = drp.dynamicResidualSaturation(Sor_static, D);
    
    // Sor_dyn = Sor_static · (1 - D) = 0.2 · 0.5 = 0.1
    EXPECT_NEAR(Sor_dynamic, 0.1, 1e-10);
}

TEST_F(DynamicRelPermTest, DynamicEndpoint) {
    double kr_max_static = 0.8;
    double D = 0.3;
    
    double kr_max_dynamic = drp.dynamicEndpoint(kr_max_static, D);
    
    // Dynamic endpoint should be higher
    EXPECT_GT(kr_max_dynamic, kr_max_static);
    EXPECT_LE(kr_max_dynamic, 1.0);
}

TEST_F(DynamicRelPermTest, VelocityDependentRelPerm) {
    double S_w = 0.6;
    double mu = 0.001;
    double v = 1e-4;
    double sigma = 0.03;
    
    double kr_static = 0.5;  // Static kr at this saturation
    double kr_dynamic = drp.modifyRelativePermeability(kr_static, S_w, mu, v, sigma);
    
    // Dynamic kr should be different from static
    // At higher velocities, dynamic effects become significant
    EXPECT_NE(kr_dynamic, kr_static);
}

TEST_F(DynamicRelPermTest, MobilityNumber) {
    double mu_o = 0.01;   // Pa·s (oil)
    double mu_w = 0.001;  // Pa·s (water)
    double kr_o = 0.3;
    double kr_w = 0.5;
    
    double M = drp.mobilityRatio(kr_w, mu_w, kr_o, mu_o);
    
    // M = (kr_w/μ_w) / (kr_o/μ_o) = (0.5/0.001) / (0.3/0.01) = 500/30 ≈ 16.7
    EXPECT_NEAR(M, 16.67, 0.1);
}

TEST_F(DynamicRelPermTest, FingeringModel) {
    DynamicRelativePermeability::Parameters params;
    params.model = DynamicRelPermModel::FINGERING;
    params.fingering_exponent = 0.5;
    drp.setParameters(params);
    
    double mu_o = 0.01;
    double mu_w = 0.001;
    double kr_w = 0.5;
    double kr_o = 0.3;
    
    double fractional_flow = drp.fractionalFlowWithFingering(kr_w, mu_w, kr_o, mu_o);
    
    // Should be between 0 and 1
    EXPECT_GE(fractional_flow, 0.0);
    EXPECT_LE(fractional_flow, 1.0);
}

TEST_F(DynamicRelPermTest, Configure) {
    std::map<std::string, std::string> config;
    config["dynamic_relperm_model"] = "fingering";
    config["fingering_exponent"] = "0.6";
    
    DynamicRelativePermeability drp2;
    drp2.configure(config);
    
    SUCCEED();
}


// =============================================================================
// Miscible Flow Tests
// =============================================================================

class MiscibleFlowTest : public ::testing::Test {
protected:
    MiscibleFlow mf;
    
    void SetUp() override {
        MiscibleFlow::Parameters params;
        params.model = MiscibleFlowModel::TODD_LONGSTAFF;
        params.mixing_parameter = 0.67;  // ω
        params.mmp = 15e6;               // Pa (minimum miscibility pressure)
        params.dispersion_longitudinal = 0.1;  // m
        params.dispersion_transverse = 0.01;   // m
        mf.setParameters(params);
    }
};

TEST_F(MiscibleFlowTest, MiscibilityFraction) {
    double P_above_MMP = 20e6;  // Above MMP
    double P_below_MMP = 10e6;  // Below MMP
    double P_at_MMP = 15e6;     // At MMP
    
    double f_above = mf.miscibilityFraction(P_above_MMP);
    double f_below = mf.miscibilityFraction(P_below_MMP);
    double f_at = mf.miscibilityFraction(P_at_MMP);
    
    // Above MMP: fully miscible
    EXPECT_NEAR(f_above, 1.0, 0.01);
    
    // Below MMP: immiscible
    EXPECT_NEAR(f_below, 0.0, 0.01);
    
    // At MMP: transition
    EXPECT_NEAR(f_at, 0.5, 0.1);
}

TEST_F(MiscibleFlowTest, EffectiveViscosityToddLongstaff) {
    double mu_o = 0.01;   // Pa·s
    double mu_s = 0.0001; // Pa·s (solvent)
    double S_s = 0.3;     // Solvent saturation
    double f_misc = 1.0;  // Fully miscible
    
    double mu_eff = mf.effectiveViscosity(mu_o, mu_s, S_s, f_misc);
    
    // Should be between oil and solvent viscosity
    EXPECT_LT(mu_eff, mu_o);
    EXPECT_GT(mu_eff, mu_s);
}

TEST_F(MiscibleFlowTest, EffectiveViscosityImmiscible) {
    double mu_o = 0.01;
    double mu_s = 0.0001;
    double S_s = 0.3;
    double f_misc = 0.0;  // Immiscible
    
    // In immiscible case, behaves like two-phase flow
    double mu_eff = mf.effectiveViscosity(mu_o, mu_s, S_s, f_misc);
    
    // Should be closer to oil viscosity (no mixing)
    EXPECT_GT(mu_eff, 0.005);  // Not as low as miscible case
}

TEST_F(MiscibleFlowTest, EffectiveDensity) {
    double rho_o = 800.0;   // kg/m³
    double rho_s = 300.0;   // kg/m³
    double S_s = 0.4;
    double f_misc = 0.8;
    
    double rho_eff = mf.effectiveDensity(rho_o, rho_s, S_s, f_misc);
    
    // Should be between oil and solvent density
    EXPECT_LT(rho_eff, rho_o);
    EXPECT_GT(rho_eff, rho_s);
}

TEST_F(MiscibleFlowTest, DispersionCoefficient) {
    double v[3] = {1e-5, 0.0, 0.0};  // m/s velocity in x-direction
    double D_mol = 1e-9;             // m²/s molecular diffusion
    
    double D_eff[3][3];
    mf.dispersionTensor(v, D_mol, D_eff);
    
    // Longitudinal dispersion (x-x) should be largest
    EXPECT_GT(D_eff[0][0], D_mol);
    
    // Transverse dispersion (y-y, z-z) should be smaller
    EXPECT_LT(D_eff[1][1], D_eff[0][0]);
    EXPECT_LT(D_eff[2][2], D_eff[0][0]);
}

TEST_F(MiscibleFlowTest, FirstContactMiscibility) {
    MiscibleFlow::Parameters params;
    params.model = MiscibleFlowModel::FIRST_CONTACT;
    params.mmp = 15e6;
    mf.setParameters(params);
    
    double f_below = mf.miscibilityFraction(14e6);
    double f_above = mf.miscibilityFraction(16e6);
    
    // Sharp transition at MMP
    EXPECT_NEAR(f_below, 0.0, 0.01);
    EXPECT_NEAR(f_above, 1.0, 0.01);
}

TEST_F(MiscibleFlowTest, MultiContactMiscibility) {
    MiscibleFlow::Parameters params;
    params.model = MiscibleFlowModel::MULTI_CONTACT;
    params.mmp = 15e6;
    params.fcm_pressure = 20e6;  // First contact miscibility pressure
    mf.setParameters(params);
    
    double f = mf.miscibilityFraction(17e6);  // Between MMP and FCM
    
    // Gradual transition
    EXPECT_GT(f, 0.0);
    EXPECT_LT(f, 1.0);
}

TEST_F(MiscibleFlowTest, Configure) {
    std::map<std::string, std::string> config;
    config["miscible_model"] = "first_contact";
    config["mmp"] = "12e6";
    config["mixing_parameter"] = "0.5";
    
    MiscibleFlow mf2;
    mf2.configure(config);
    
    SUCCEED();
}


// =============================================================================
// Configuration Validation Tests
// =============================================================================

class FluidFlowConfigTest : public ::testing::Test {
protected:
    HighFidelityFluidFlowConfig config;
};

TEST_F(FluidFlowConfigTest, ValidConfiguration) {
    config.enable_non_darcy = true;
    config.non_darcy_params.model = NonDarcyModel::FORCHHEIMER;
    config.non_darcy_params.forchheimer_beta = 1e8;
    
    config.enable_dual_porosity = true;
    config.dual_porosity_params.matrix_porosity = 0.1;
    config.dual_porosity_params.fracture_porosity = 0.01;
    
    std::string error_msg;
    bool valid = config.validate(error_msg);
    
    EXPECT_TRUE(valid);
    EXPECT_TRUE(error_msg.empty());
}

TEST_F(FluidFlowConfigTest, InvalidPorosity) {
    config.enable_dual_porosity = true;
    config.dual_porosity_params.matrix_porosity = 1.5;  // Invalid: > 1
    
    std::string error_msg;
    bool valid = config.validate(error_msg);
    
    EXPECT_FALSE(valid);
    EXPECT_FALSE(error_msg.empty());
}

TEST_F(FluidFlowConfigTest, InvalidNegativePermeability) {
    config.enable_dual_porosity = true;
    config.dual_porosity_params.matrix_permeability = -1e-15;  // Invalid: negative
    
    std::string error_msg;
    bool valid = config.validate(error_msg);
    
    EXPECT_FALSE(valid);
    EXPECT_FALSE(error_msg.empty());
}

TEST_F(FluidFlowConfigTest, ParseConfig) {
    std::map<std::string, std::string> cfg;
    cfg["enable_non_darcy"] = "true";
    cfg["enable_dual_porosity"] = "true";
    cfg["enable_non_isothermal"] = "false";
    cfg["enable_dynamic_relperm"] = "true";
    cfg["enable_miscible"] = "false";
    
    config.parseConfig(cfg);
    
    EXPECT_TRUE(config.enable_non_darcy);
    EXPECT_TRUE(config.enable_dual_porosity);
    EXPECT_FALSE(config.enable_non_isothermal);
    EXPECT_TRUE(config.enable_dynamic_relperm);
    EXPECT_FALSE(config.enable_miscible);
}


// =============================================================================
// Factory Function Tests
// =============================================================================

TEST(FluidFlowFactoryTest, CreateNonDarcyFlow) {
    std::map<std::string, std::string> config;
    config["non_darcy_model"] = "forchheimer";
    config["forchheimer_beta"] = "1e9";
    
    auto ndf = createNonDarcyFlow(config);
    
    ASSERT_NE(ndf, nullptr);
    
    double correction = ndf->calculateCorrection(0.01, 1000.0, 0.001, 1e-12);
    EXPECT_GT(correction, 0.0);
}

TEST(FluidFlowFactoryTest, CreateDualPorosity) {
    std::map<std::string, std::string> config;
    config["dual_porosity_model"] = "warren_root";
    
    auto dpp = createDualPorosityModel(config);
    
    ASSERT_NE(dpp, nullptr);
}

TEST(FluidFlowFactoryTest, CreateMiscibleFlow) {
    std::map<std::string, std::string> config;
    config["miscible_model"] = "todd_longstaff";
    config["mixing_parameter"] = "0.67";
    
    auto mf = createMiscibleFlow(config);
    
    ASSERT_NE(mf, nullptr);
}
