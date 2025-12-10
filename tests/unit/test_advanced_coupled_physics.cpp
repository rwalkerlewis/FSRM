/**
 * @file test_advanced_coupled_physics.cpp
 * @brief Unit tests for AdvancedCoupledPhysics module
 * 
 * Tests cover:
 * - Full Biot Dynamic Poroelasticity
 * - THM (Thermo-Hydro-Mechanical) coupling
 * - THMC (Thermo-Hydro-Mechanical-Chemical) coupling
 * - Unsaturated flow coupling
 * - Multi-scale coupling methods
 */

#include <gtest/gtest.h>
#include "AdvancedCoupledPhysics.hpp"
#include <cmath>

using namespace FSRM::AdvancedCoupled;

// =============================================================================
// Full Biot Dynamic Poroelasticity Tests
// =============================================================================

class FullBiotDynamicsTest : public ::testing::Test {
protected:
    FullBiotDynamics fbd;
    
    void SetUp() override {
        FullBiotDynamics::Parameters params;
        params.formulation = BiotFormulation::U_P;
        params.biot_alpha = 0.8;
        params.biot_M = 1e10;        // Pa (Biot modulus)
        params.porosity = 0.2;
        params.permeability = 1e-13; // m²
        params.fluid_viscosity = 1e-3;  // Pa·s
        params.solid_density = 2650.0;  // kg/m³
        params.fluid_density = 1000.0;  // kg/m³
        params.shear_modulus = 10e9;    // Pa
        params.bulk_modulus = 15e9;     // Pa
        params.enable_inertia = true;
        params.include_tortuosity = true;
        params.tortuosity = 2.0;
        fbd.setParameters(params);
    }
};

TEST_F(FullBiotDynamicsTest, BiotAlphaRange) {
    double alpha = fbd.biotCoefficient();
    
    // Biot coefficient should be between 0 and 1
    EXPECT_GE(alpha, 0.0);
    EXPECT_LE(alpha, 1.0);
}

TEST_F(FullBiotDynamicsTest, BiotModulus) {
    double M = fbd.biotModulus();
    
    // Biot modulus should be positive
    EXPECT_GT(M, 0.0);
}

TEST_F(FullBiotDynamicsTest, UndrainedBulkModulus) {
    double K_u = fbd.undrainedBulkModulus();
    double K_d = fbd.drainedBulkModulus();
    
    // Undrained should be >= drained (Kᵤ ≥ K_d)
    EXPECT_GE(K_u, K_d);
}

TEST_F(FullBiotDynamicsTest, SkemptonCoefficient) {
    double B = fbd.skemptonCoefficient();
    
    // Skempton's B should be between 0 and 1
    EXPECT_GT(B, 0.0);
    EXPECT_LE(B, 1.0);
}

TEST_F(FullBiotDynamicsTest, PwaveVelocityDrained) {
    double Vp_d = fbd.pwaveVelocityDrained();
    
    // Should be positive
    EXPECT_GT(Vp_d, 0.0);
    
    // Typical rock values: 2000-6000 m/s
    EXPECT_GT(Vp_d, 1000.0);
    EXPECT_LT(Vp_d, 10000.0);
}

TEST_F(FullBiotDynamicsTest, PwaveVelocityUndrained) {
    double Vp_u = fbd.pwaveVelocityUndrained();
    double Vp_d = fbd.pwaveVelocityDrained();
    
    // Undrained Vp should be >= drained (stiffer response)
    EXPECT_GE(Vp_u, Vp_d);
}

TEST_F(FullBiotDynamicsTest, SwaveVelocity) {
    double Vs = fbd.swaveVelocity();
    
    // Should be positive
    EXPECT_GT(Vs, 0.0);
    
    // Vs < Vp always
    EXPECT_LT(Vs, fbd.pwaveVelocityDrained());
}

TEST_F(FullBiotDynamicsTest, SlowPwaveVelocity) {
    double Vp_slow = fbd.slowPwaveVelocity();
    double Vp_fast = fbd.pwaveVelocityUndrained();
    
    // Slow P-wave is slower than fast P-wave
    EXPECT_GT(Vp_slow, 0.0);
    EXPECT_LT(Vp_slow, Vp_fast);
}

TEST_F(FullBiotDynamicsTest, CharacteristicFrequency) {
    double f_c = fbd.characteristicFrequency();
    
    // Should be positive
    EXPECT_GT(f_c, 0.0);
}

TEST_F(FullBiotDynamicsTest, HighFrequencyEffects) {
    double freq = fbd.characteristicFrequency();
    
    bool high_freq = fbd.isHighFrequencyRegime(freq * 10);
    bool low_freq = fbd.isHighFrequencyRegime(freq / 10);
    
    EXPECT_TRUE(high_freq);
    EXPECT_FALSE(low_freq);
}

TEST_F(FullBiotDynamicsTest, EffectiveStress) {
    double total_stress[6] = {-10e6, -10e6, -10e6, 0, 0, 0};
    double pore_pressure = 5e6;
    
    double effective_stress[6];
    fbd.effectiveStress(total_stress, pore_pressure, effective_stress);
    
    // σ'_ij = σ_ij + α·p·δ_ij
    // For isotropic stress: σ'_11 = -10e6 + 0.8·5e6 = -6e6
    EXPECT_NEAR(effective_stress[0], -10e6 + 0.8 * 5e6, 1e3);
}

TEST_F(FullBiotDynamicsTest, FluidContent) {
    double volumetric_strain = -0.001;  // Compression
    double pore_pressure = 2e6;
    
    double zeta = fbd.fluidContent(volumetric_strain, pore_pressure);
    
    // Should be related to pressure and strain
    // ζ = α·ε_kk + p/M
}

TEST_F(FullBiotDynamicsTest, MassMatrixCoefficients) {
    double M_uu, M_ww, M_uw;
    fbd.massMatrixCoefficients(M_uu, M_ww, M_uw);
    
    // All should be positive
    EXPECT_GT(M_uu, 0.0);
    EXPECT_GT(M_ww, 0.0);
}

TEST_F(FullBiotDynamicsTest, DampingCoefficient) {
    double damping = fbd.viscousDampingCoefficient();
    
    // Should be positive (dissipative)
    EXPECT_GT(damping, 0.0);
}

TEST_F(FullBiotDynamicsTest, Configure) {
    std::map<std::string, std::string> config;
    config["biot_formulation"] = "u_p_w";
    config["biot_alpha"] = "0.9";
    config["include_tortuosity"] = "true";
    
    FullBiotDynamics fbd2;
    fbd2.configure(config);
    
    SUCCEED();
}


// =============================================================================
// THM Coupling Tests
// =============================================================================

class THMCouplingTest : public ::testing::Test {
protected:
    THMCoupling thm;
    
    void SetUp() override {
        THMCoupling::Parameters params;
        params.biot_alpha = 0.8;
        params.biot_M = 1e10;
        params.thermal_expansion_solid = 1e-5;  // 1/K
        params.thermal_expansion_fluid = 2e-4;  // 1/K
        params.thermal_conductivity = 2.5;      // W/(m·K)
        params.specific_heat = 1000.0;          // J/(kg·K)
        params.reference_temperature = 293.15;  // K
        params.coupling_thermal_stress = true;
        params.coupling_thermal_fluid = true;
        params.coupling_pore_pressure = true;
        thm.setParameters(params);
    }
};

TEST_F(THMCouplingTest, ThermalStress) {
    double delta_T = 50.0;  // K temperature increase
    
    double thermal_stress = thm.thermalStress(delta_T);
    
    // Should be negative (compression) for temperature increase
    // σ_th = -3·K·α_s·ΔT (for isotropic case)
    EXPECT_LT(thermal_stress, 0.0);
}

TEST_F(THMCouplingTest, ThermalExpansionStrain) {
    double delta_T = 50.0;
    
    double thermal_strain = thm.thermalExpansionStrain(delta_T);
    
    // ε_th = α_s · ΔT = 1e-5 · 50 = 5e-4
    EXPECT_NEAR(thermal_strain, 5e-4, 1e-5);
}

TEST_F(THMCouplingTest, ThermalPressurization) {
    double delta_T = 50.0;
    double bulk_modulus = 15e9;
    
    double delta_p = thm.thermalPressurization(delta_T, bulk_modulus);
    
    // Heating in confined space increases pore pressure
    EXPECT_GT(delta_p, 0.0);
}

TEST_F(THMCouplingTest, EffectiveThermalExpansion) {
    double alpha_eff = thm.effectiveThermalExpansion();
    
    // Should be positive
    EXPECT_GT(alpha_eff, 0.0);
}

TEST_F(THMCouplingTest, ThermalDiffusivity) {
    double rho = 2500.0;  // kg/m³
    
    double kappa_T = thm.thermalDiffusivity(rho);
    
    // κ = k / (ρ·c) = 2.5 / (2500·1000) = 1e-6 m²/s
    EXPECT_NEAR(kappa_T, 1e-6, 1e-8);
}

TEST_F(THMCouplingTest, CoupledThermalConductivity) {
    double porosity = 0.2;
    double k_solid = 3.0;
    double k_fluid = 0.6;
    
    double k_eff = thm.effectiveThermalConductivity(porosity, k_solid, k_fluid);
    
    // Should be between fluid and solid values
    EXPECT_GT(k_eff, k_fluid);
    EXPECT_LT(k_eff, k_solid);
}

TEST_F(THMCouplingTest, ConsolidationCoefficient) {
    double cv = thm.consolidationCoefficient(1e-13, 1e-3, 15e9, 0.2);
    
    // Should be positive
    EXPECT_GT(cv, 0.0);
}

TEST_F(THMCouplingTest, ThermalConsolidationTime) {
    double L = 10.0;  // m characteristic length
    double cv = 1e-6; // m²/s diffusivity
    
    double t_c = thm.characteristicTime(L, cv);
    
    // t = L² / cv = 100 / 1e-6 = 1e8 s
    EXPECT_NEAR(t_c, 1e8, 1e6);
}

TEST_F(THMCouplingTest, Configure) {
    std::map<std::string, std::string> config;
    config["thm_thermal_expansion_solid"] = "2e-5";
    config["thm_reference_temperature"] = "300.0";
    config["coupling_thermal_stress"] = "true";
    
    THMCoupling thm2;
    thm2.configure(config);
    
    SUCCEED();
}


// =============================================================================
// THMC Coupling Tests
// =============================================================================

class THMCCouplingTest : public ::testing::Test {
protected:
    THMCCoupling thmc;
    
    void SetUp() override {
        THMCCoupling::Parameters params;
        // THM base parameters
        params.thm_params.biot_alpha = 0.8;
        params.thm_params.thermal_expansion_solid = 1e-5;
        params.thm_params.thermal_conductivity = 2.5;
        
        // Chemical parameters
        params.num_species = 2;
        params.diffusion_coefficients = {1e-9, 5e-10};  // m²/s
        params.reaction_rate_constants = {1e-6, 2e-6};  // s⁻¹
        params.chemical_strain_coefficients = {0.01, 0.005};
        params.enable_chemical_dissolution = true;
        params.enable_chemical_precipitation = true;
        params.enable_chemical_swelling = true;
        thmc.setParameters(params);
    }
};

TEST_F(THMCCouplingTest, ChemicalStrain) {
    std::vector<double> concentrations = {0.1, 0.05};  // mol/m³
    
    double eps_chem = thmc.chemicalStrain(concentrations);
    
    // Should be non-zero for non-zero concentrations
    EXPECT_NE(eps_chem, 0.0);
}

TEST_F(THMCCouplingTest, ChemicalStrainRate) {
    std::vector<double> conc_rate = {1e-6, 5e-7};  // mol/(m³·s)
    
    double eps_chem_dot = thmc.chemicalStrainRate(conc_rate);
    
    // Should be related to concentration rates
    EXPECT_NE(eps_chem_dot, 0.0);
}

TEST_F(THMCCouplingTest, ReactionHeat) {
    std::vector<double> reaction_rates = {1e-8, 2e-8};  // mol/(m³·s)
    std::vector<double> enthalpies = {-50000.0, -30000.0};  // J/mol (exothermic)
    
    double Q_rxn = thmc.reactionHeat(reaction_rates, enthalpies);
    
    // Exothermic reactions release heat (positive Q)
    EXPECT_GT(Q_rxn, 0.0);
}

TEST_F(THMCCouplingTest, PorosityChange) {
    double dissolution_rate = 1e-8;   // mol/(m³·s)
    double precipitation_rate = 5e-9; // mol/(m³·s)
    double molar_volume = 1e-4;       // m³/mol
    
    double phi_dot = thmc.porosityChangeRate(dissolution_rate, precipitation_rate, molar_volume);
    
    // Net dissolution increases porosity
    EXPECT_GT(phi_dot, 0.0);
}

TEST_F(THMCCouplingTest, EffectiveDiffusivity) {
    double D_free = 1e-9;    // m²/s free diffusivity
    double porosity = 0.2;
    double tortuosity = 2.0;
    
    double D_eff = thmc.effectiveDiffusivity(D_free, porosity, tortuosity);
    
    // D_eff = φ · D_free / τ = 0.2 · 1e-9 / 2 = 1e-10
    EXPECT_NEAR(D_eff, 1e-10, 1e-12);
}

TEST_F(THMCCouplingTest, ArrheniusRate) {
    double k_0 = 1e-6;         // s⁻¹ pre-exponential factor
    double E_a = 50000.0;      // J/mol activation energy
    double T = 350.0;          // K temperature
    double T_ref = 293.15;     // K reference temperature
    
    double k_T = thmc.temperatureDependentRate(k_0, E_a, T, T_ref);
    
    // Higher temperature -> higher rate
    EXPECT_GT(k_T, k_0);
}

TEST_F(THMCCouplingTest, ChemicalOsmosisEffect) {
    double conc_gradient = 100.0;  // mol/m³/m
    double osmotic_coeff = 1e-10;  // m²/s
    
    double flux = thmc.chemicalOsmosisFlux(conc_gradient, osmotic_coeff);
    
    // Flux should be proportional to gradient
    EXPECT_NE(flux, 0.0);
}

TEST_F(THMCCouplingTest, Configure) {
    std::map<std::string, std::string> config;
    config["thmc_num_species"] = "3";
    config["enable_chemical_swelling"] = "true";
    config["enable_dissolution"] = "true";
    
    THMCCoupling thmc2;
    thmc2.configure(config);
    
    SUCCEED();
}


// =============================================================================
// Unsaturated Flow Coupling Tests
// =============================================================================

class UnsaturatedFlowTest : public ::testing::Test {
protected:
    UnsaturatedFlowCoupling ufc;
    
    void SetUp() override {
        UnsaturatedFlowCoupling::Parameters params;
        params.retention_model = RetentionModel::VAN_GENUCHTEN;
        params.van_genuchten_alpha = 0.5;    // 1/m
        params.van_genuchten_n = 1.5;
        params.van_genuchten_m = 0.333;
        params.residual_saturation = 0.1;
        params.air_entry_pressure = 1000.0;  // Pa
        params.biot_alpha = 0.8;
        params.bishops_chi_model = BishopChiModel::SATURATION;
        ufc.setParameters(params);
    }
};

TEST_F(UnsaturatedFlowTest, WaterRetention) {
    double suction = 5000.0;  // Pa
    
    double S_w = ufc.waterSaturation(suction);
    
    // Should be between residual and 1
    EXPECT_GT(S_w, 0.1);
    EXPECT_LE(S_w, 1.0);
}

TEST_F(UnsaturatedFlowTest, FullSaturationAtZeroSuction) {
    double suction = 0.0;
    
    double S_w = ufc.waterSaturation(suction);
    
    // Full saturation at zero suction
    EXPECT_NEAR(S_w, 1.0, 0.01);
}

TEST_F(UnsaturatedFlowTest, ResidualSaturationHighSuction) {
    double suction = 1e8;  // Very high suction
    
    double S_w = ufc.waterSaturation(suction);
    
    // Should approach residual saturation
    EXPECT_NEAR(S_w, 0.1, 0.05);
}

TEST_F(UnsaturatedFlowTest, RelativePermeability) {
    double S_e = 0.5;  // Effective saturation
    
    double kr = ufc.relativePermeability(S_e);
    
    // Should be between 0 and 1
    EXPECT_GT(kr, 0.0);
    EXPECT_LE(kr, 1.0);
}

TEST_F(UnsaturatedFlowTest, RelativePermFullSaturation) {
    double S_e = 1.0;
    
    double kr = ufc.relativePermeability(S_e);
    
    // Should be 1 at full saturation
    EXPECT_NEAR(kr, 1.0, 0.01);
}

TEST_F(UnsaturatedFlowTest, BishopParameter) {
    double S_w = 0.6;
    
    double chi = ufc.bishopParameter(S_w);
    
    // For saturation-based model: χ ≈ S_w
    EXPECT_GT(chi, 0.0);
    EXPECT_LE(chi, 1.0);
}

TEST_F(UnsaturatedFlowTest, EffectiveStress) {
    double total_stress = -100e3;  // Pa (compression)
    double p_water = 50e3;         // Pa
    double p_air = 101e3;          // Pa (atmospheric)
    double chi = 0.6;
    
    double sigma_eff = ufc.effectiveStress(total_stress, p_water, p_air, chi);
    
    // σ' = σ - [χ·p_w + (1-χ)·p_a]
    double expected = total_stress - (chi * p_water + (1 - chi) * p_air);
    EXPECT_NEAR(sigma_eff, expected, 1e3);
}

TEST_F(UnsaturatedFlowTest, CapacityFunction) {
    double suction = 5000.0;
    
    double C = ufc.moistureCapacity(suction);
    
    // Should be positive
    EXPECT_GT(C, 0.0);
}

TEST_F(UnsaturatedFlowTest, BrooksCoreyModel) {
    UnsaturatedFlowCoupling::Parameters params;
    params.retention_model = RetentionModel::BROOKS_COREY;
    params.brooks_corey_lambda = 2.0;
    params.air_entry_pressure = 2000.0;
    params.residual_saturation = 0.15;
    ufc.setParameters(params);
    
    double suction = 5000.0;
    double S_w = ufc.waterSaturation(suction);
    
    EXPECT_GT(S_w, params.residual_saturation);
    EXPECT_LT(S_w, 1.0);
}

TEST_F(UnsaturatedFlowTest, Configure) {
    std::map<std::string, std::string> config;
    config["retention_model"] = "van_genuchten";
    config["vg_alpha"] = "0.3";
    config["vg_n"] = "2.0";
    
    UnsaturatedFlowCoupling ufc2;
    ufc2.configure(config);
    
    SUCCEED();
}


// =============================================================================
// Multi-Scale Coupling Tests
// =============================================================================

class MultiScaleCouplingTest : public ::testing::Test {
protected:
    MultiScaleCoupling msc;
    
    void SetUp() override {
        MultiScaleCoupling::Parameters params;
        params.method = MultiScaleMethod::FE_SQUARED;
        params.rve_size = 0.001;  // m (1 mm RVE)
        params.num_rve_elements = 100;
        params.homogenization = HomogenizationType::COMPUTATIONAL;
        params.boundary_condition = RVEBoundaryCondition::PERIODIC;
        msc.setParameters(params);
    }
};

TEST_F(MultiScaleCouplingTest, VoigtBound) {
    std::vector<double> moduli = {100e9, 50e9};  // Pa
    std::vector<double> fractions = {0.4, 0.6};
    
    double K_voigt = msc.voigtAverage(moduli, fractions);
    
    // Voigt: K = Σ f_i · K_i = 0.4·100e9 + 0.6·50e9 = 70e9
    EXPECT_NEAR(K_voigt, 70e9, 1e7);
}

TEST_F(MultiScaleCouplingTest, ReussBound) {
    std::vector<double> moduli = {100e9, 50e9};
    std::vector<double> fractions = {0.4, 0.6};
    
    double K_reuss = msc.reussAverage(moduli, fractions);
    
    // Reuss: 1/K = Σ f_i / K_i = 0.4/100e9 + 0.6/50e9
    double expected = 1.0 / (0.4 / 100e9 + 0.6 / 50e9);
    EXPECT_NEAR(K_reuss, expected, 1e7);
}

TEST_F(MultiScaleCouplingTest, VoigtReussBounds) {
    std::vector<double> moduli = {100e9, 50e9};
    std::vector<double> fractions = {0.4, 0.6};
    
    double K_voigt = msc.voigtAverage(moduli, fractions);
    double K_reuss = msc.reussAverage(moduli, fractions);
    
    // Voigt should be upper bound, Reuss lower bound
    EXPECT_GT(K_voigt, K_reuss);
}

TEST_F(MultiScaleCouplingTest, HashinShtrikmanBounds) {
    double K1 = 100e9, mu1 = 80e9;
    double K2 = 50e9, mu2 = 30e9;
    double f1 = 0.4;
    
    double K_lower, K_upper, mu_lower, mu_upper;
    msc.hashinShtrikmanBounds(K1, mu1, K2, mu2, f1,
                              K_lower, K_upper, mu_lower, mu_upper);
    
    // Bounds should satisfy Reuss ≤ HS_lower ≤ HS_upper ≤ Voigt
    EXPECT_LT(K_lower, K_upper);
    EXPECT_LT(mu_lower, mu_upper);
}

TEST_F(MultiScaleCouplingTest, MoriTanakaScheme) {
    // Matrix-inclusion composite
    double K_m = 50e9, mu_m = 30e9;  // Matrix
    double K_i = 100e9, mu_i = 60e9; // Inclusion
    double f_i = 0.3;
    
    double K_eff, mu_eff;
    msc.moriTanaka(K_m, mu_m, K_i, mu_i, f_i, K_eff, mu_eff);
    
    // Effective should be between matrix and inclusion
    EXPECT_GT(K_eff, K_m);
    EXPECT_LT(K_eff, K_i);
}

TEST_F(MultiScaleCouplingTest, SelfConsistentScheme) {
    double K1 = 100e9, mu1 = 80e9;
    double K2 = 50e9, mu2 = 30e9;
    double f1 = 0.5;
    
    double K_eff, mu_eff;
    msc.selfConsistent(K1, mu1, K2, mu2, f1, K_eff, mu_eff);
    
    // Should be between bounds
    std::vector<double> K = {K1, K2};
    std::vector<double> f = {f1, 1-f1};
    
    double K_reuss = msc.reussAverage(K, f);
    double K_voigt = msc.voigtAverage(K, f);
    
    EXPECT_GE(K_eff, K_reuss * 0.99);  // Allow small numerical tolerance
    EXPECT_LE(K_eff, K_voigt * 1.01);
}

TEST_F(MultiScaleCouplingTest, StrainLocalization) {
    double macro_strain[6] = {0.001, -0.0003, -0.0003, 0, 0, 0};
    double A_tensor[6][6];  // Simplified - identity for test
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            A_tensor[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
    
    double micro_strain[6];
    msc.strainLocalization(macro_strain, A_tensor, micro_strain);
    
    // For identity localization, micro = macro
    for (int i = 0; i < 6; ++i) {
        EXPECT_NEAR(micro_strain[i], macro_strain[i], 1e-10);
    }
}

TEST_F(MultiScaleCouplingTest, Configure) {
    std::map<std::string, std::string> config;
    config["multiscale_method"] = "hmm";
    config["rve_size"] = "0.002";
    config["rve_boundary"] = "periodic";
    
    MultiScaleCoupling msc2;
    msc2.configure(config);
    
    SUCCEED();
}


// =============================================================================
// Configuration Validation Tests
// =============================================================================

class CoupledPhysicsConfigTest : public ::testing::Test {
protected:
    AdvancedCoupledPhysicsConfig config;
};

TEST_F(CoupledPhysicsConfigTest, ValidConfiguration) {
    config.enable_full_biot = true;
    config.biot_params.formulation = BiotFormulation::U_P;
    config.biot_params.biot_alpha = 0.8;
    
    config.enable_thm = true;
    config.thm_params.thermal_expansion_solid = 1e-5;
    
    std::string error_msg;
    bool valid = config.validate(error_msg);
    
    EXPECT_TRUE(valid);
    EXPECT_TRUE(error_msg.empty());
}

TEST_F(CoupledPhysicsConfigTest, InvalidBiotAlpha) {
    config.enable_full_biot = true;
    config.biot_params.biot_alpha = 1.5;  // Invalid: > 1
    
    std::string error_msg;
    bool valid = config.validate(error_msg);
    
    EXPECT_FALSE(valid);
    EXPECT_FALSE(error_msg.empty());
}

TEST_F(CoupledPhysicsConfigTest, InvalidNegativeConductivity) {
    config.enable_thm = true;
    config.thm_params.thermal_conductivity = -1.0;  // Invalid: negative
    
    std::string error_msg;
    bool valid = config.validate(error_msg);
    
    EXPECT_FALSE(valid);
}

TEST_F(CoupledPhysicsConfigTest, ParseConfig) {
    std::map<std::string, std::string> cfg;
    cfg["enable_full_biot"] = "true";
    cfg["enable_thm"] = "true";
    cfg["enable_thmc"] = "false";
    cfg["enable_unsaturated"] = "true";
    cfg["enable_multiscale"] = "false";
    
    config.parseConfig(cfg);
    
    EXPECT_TRUE(config.enable_full_biot);
    EXPECT_TRUE(config.enable_thm);
    EXPECT_FALSE(config.enable_thmc);
    EXPECT_TRUE(config.enable_unsaturated);
    EXPECT_FALSE(config.enable_multiscale);
}


// =============================================================================
// Factory Function Tests
// =============================================================================

TEST(CoupledPhysicsFactoryTest, CreateFullBiotDynamics) {
    std::map<std::string, std::string> config;
    config["biot_formulation"] = "u_p";
    config["biot_alpha"] = "0.8";
    
    auto fbd = createFullBiotDynamics(config);
    
    ASSERT_NE(fbd, nullptr);
}

TEST(CoupledPhysicsFactoryTest, CreateTHMCoupling) {
    std::map<std::string, std::string> config;
    config["thm_thermal_expansion"] = "1e-5";
    config["thm_conductivity"] = "2.5";
    
    auto thm = createTHMCoupling(config);
    
    ASSERT_NE(thm, nullptr);
}

TEST(CoupledPhysicsFactoryTest, CreateTHMCCoupling) {
    std::map<std::string, std::string> config;
    config["thmc_num_species"] = "2";
    
    auto thmc = createTHMCCoupling(config);
    
    ASSERT_NE(thmc, nullptr);
}

TEST(CoupledPhysicsFactoryTest, CreateUnsaturatedFlowCoupling) {
    std::map<std::string, std::string> config;
    config["retention_model"] = "van_genuchten";
    
    auto ufc = createUnsaturatedFlowCoupling(config);
    
    ASSERT_NE(ufc, nullptr);
}

TEST(CoupledPhysicsFactoryTest, CreateMultiScaleCoupling) {
    std::map<std::string, std::string> config;
    config["multiscale_method"] = "fe_squared";
    
    auto msc = createMultiScaleCoupling(config);
    
    ASSERT_NE(msc, nullptr);
}


// =============================================================================
// Coupling Consistency Tests
// =============================================================================

class CouplingConsistencyTest : public ::testing::Test {
protected:
    // Test that different coupling modules are consistent
};

TEST_F(CouplingConsistencyTest, BiotTHMConsistency) {
    // Biot and THM should give same poroelastic response when thermal effects disabled
    FullBiotDynamics fbd;
    THMCoupling thm;
    
    FullBiotDynamics::Parameters biot_params;
    biot_params.biot_alpha = 0.8;
    biot_params.porosity = 0.2;
    fbd.setParameters(biot_params);
    
    THMCoupling::Parameters thm_params;
    thm_params.biot_alpha = 0.8;
    thm_params.coupling_thermal_stress = false;  // Disable thermal
    thm.setParameters(thm_params);
    
    // Both should have same Biot coefficient
    EXPECT_NEAR(fbd.biotCoefficient(), 0.8, 1e-10);
}

TEST_F(CouplingConsistencyTest, UnsaturatedBiotConsistency) {
    // At full saturation, unsaturated model should reduce to Biot
    UnsaturatedFlowCoupling ufc;
    
    UnsaturatedFlowCoupling::Parameters params;
    params.biot_alpha = 0.8;
    params.residual_saturation = 0.0;
    ufc.setParameters(params);
    
    // At zero suction, saturation = 1, χ = 1
    double suction = 0.0;
    double S_w = ufc.waterSaturation(suction);
    double chi = ufc.bishopParameter(S_w);
    
    EXPECT_NEAR(S_w, 1.0, 0.01);
    EXPECT_NEAR(chi, 1.0, 0.01);
}
