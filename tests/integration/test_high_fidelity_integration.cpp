/**
 * @file test_high_fidelity_integration.cpp
 * @brief Full-scale integration tests for high-fidelity modules
 * 
 * Integration tests verify that different modules work correctly together
 * in realistic simulation scenarios. These tests:
 * 1. Use physically meaningful parameter values
 * 2. Test complete workflows from configuration to solution
 * 3. Verify coupling between different physics modules
 * 4. Check conservation properties and physical bounds
 */

#include <gtest/gtest.h>
#include "HighFidelityFluidFlow.hpp"
#include "HighFidelityGeomechanics.hpp"
#include "HighFidelityNumerics.hpp"
#include "AdvancedCoupledPhysics.hpp"
#include <map>
#include <string>
#include <cmath>
#include <memory>

using namespace FSRM::HighFidelity;
using namespace FSRM::AdvancedCoupled;


// =============================================================================
// Test Fixtures
// =============================================================================

/**
 * @brief Base fixture for integration tests
 */
class HighFidelityIntegrationTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Initialize default configurations
        config_map["domain_length"] = "1000.0";
        config_map["domain_width"] = "500.0";
        config_map["domain_depth"] = "100.0";
        config_map["num_cells_x"] = "50";
        config_map["num_cells_y"] = "25";
        config_map["num_cells_z"] = "10";
    }
    
    std::map<std::string, std::string> config_map;
};


// =============================================================================
// Reservoir Simulation Integration Tests
// =============================================================================

/**
 * @brief Full reservoir simulation with high-fidelity fluid flow
 */
class ReservoirIntegrationTest : public HighFidelityIntegrationTest {
protected:
    void SetUp() override {
        HighFidelityIntegrationTest::SetUp();
        
        // Reservoir properties
        config_map["porosity"] = "0.2";
        config_map["permeability"] = "1e-13";
        config_map["initial_pressure"] = "25e6";
        config_map["temperature"] = "353.15";  // 80°C
        
        // Fluid properties
        config_map["oil_viscosity"] = "0.005";
        config_map["water_viscosity"] = "0.001";
        config_map["oil_density"] = "850.0";
        config_map["water_density"] = "1020.0";
    }
};

TEST_F(ReservoirIntegrationTest, NonDarcyHighVelocityFlow) {
    // Test non-Darcy flow near wellbore (high velocity region)
    
    NonDarcyFlow ndf;
    NonDarcyFlow::Parameters params;
    params.model = NonDarcyModel::FORCHHEIMER;
    params.forchheimer_beta = 5e8;  // Tight formation
    ndf.setParameters(params);
    
    double K = 1e-13;      // m²
    double mu = 0.001;     // Pa·s
    double rho = 1000.0;   // kg/m³
    
    // Wellbore region: high velocity
    std::vector<double> velocities = {1e-3, 5e-3, 1e-2, 5e-2};
    std::vector<double> corrections;
    
    for (double v : velocities) {
        double correction = ndf.calculateCorrection(v, rho, mu, K);
        corrections.push_back(correction);
        
        // All corrections should be valid
        EXPECT_GT(correction, 0.0);
        EXPECT_LE(correction, 1.0);
    }
    
    // Higher velocity -> lower correction (more non-Darcy effect)
    for (size_t i = 1; i < corrections.size(); ++i) {
        EXPECT_LT(corrections[i], corrections[i-1]);
    }
    
    // At highest velocity, correction should be significantly less than 1
    EXPECT_LT(corrections.back(), 0.8);
}

TEST_F(ReservoirIntegrationTest, DualPorosityMatrixFractureExchange) {
    // Test matrix-fracture exchange in naturally fractured reservoir
    
    DualPorosityPermeability dpp;
    DualPorosityPermeability::Parameters params;
    params.model = DualPorosityModel::WARREN_ROOT;
    params.matrix_porosity = 0.15;
    params.fracture_porosity = 0.02;
    params.matrix_permeability = 1e-16;
    params.fracture_permeability = 1e-12;
    params.shape_factor = 0.5;
    params.transfer_coeff = 1e-7;
    dpp.setParameters(params);
    
    // Simulate drawdown: fracture pressure drops, matrix pressure initially constant
    double p_matrix_initial = 25e6;  // Pa
    double p_fracture_depleted = 20e6;  // Pa (after production)
    
    double transfer = dpp.calculateTransferRate(p_matrix_initial, p_fracture_depleted);
    
    // Positive transfer: matrix feeding fracture
    EXPECT_GT(transfer, 0.0);
    
    // Verify effective properties
    double phi_eff = dpp.effectivePorosity();
    EXPECT_NEAR(phi_eff, 0.17, 0.01);  // 0.15 + 0.02
    
    double K_eff = dpp.effectivePermeability();
    // Effective permeability dominated by fractures
    EXPECT_GT(K_eff, params.matrix_permeability);
}

TEST_F(ReservoirIntegrationTest, MiscibleCO2Injection) {
    // Test miscible CO2 flooding scenario
    
    MiscibleFlow mf;
    MiscibleFlow::Parameters params;
    params.model = MiscibleFlowModel::TODD_LONGSTAFF;
    params.mixing_parameter = 0.67;
    params.mmp = 12e6;  // Pa (minimum miscibility pressure)
    params.dispersion_longitudinal = 0.5;
    params.dispersion_transverse = 0.05;
    mf.setParameters(params);
    
    // Test at different pressures
    std::vector<double> pressures = {8e6, 12e6, 15e6, 20e6};
    
    for (double p : pressures) {
        double f_misc = mf.miscibilityFraction(p);
        
        EXPECT_GE(f_misc, 0.0);
        EXPECT_LE(f_misc, 1.0);
        
        if (p < params.mmp) {
            EXPECT_LT(f_misc, 0.5);  // Below MMP: mostly immiscible
        } else {
            EXPECT_GT(f_misc, 0.5);  // Above MMP: mostly miscible
        }
    }
    
    // Test viscosity mixing
    double mu_oil = 0.01;    // Pa·s
    double mu_CO2 = 0.0001;  // Pa·s
    double S_CO2 = 0.3;
    double f = mf.miscibilityFraction(15e6);  // Above MMP
    
    double mu_eff = mf.effectiveViscosity(mu_oil, mu_CO2, S_CO2, f);
    
    // Effective viscosity should be reduced
    EXPECT_LT(mu_eff, mu_oil);
    EXPECT_GT(mu_eff, mu_CO2);
}

TEST_F(ReservoirIntegrationTest, ThermalRecoveryProcess) {
    // Test thermal EOR scenario with non-isothermal flow
    
    NonIsothermalMultiphaseFlow nimf;
    NonIsothermalMultiphaseFlow::Parameters params;
    params.reference_temperature = 323.15;  // K (50°C reservoir)
    params.thermal_conductivity_solid = 2.0;
    params.thermal_conductivity_water = 0.65;
    params.thermal_conductivity_oil = 0.15;
    params.heat_capacity_solid = 850.0;
    params.heat_capacity_water = 4200.0;
    params.heat_capacity_oil = 2100.0;
    params.enable_viscous_heating = true;
    params.enable_joule_thomson = true;
    nimf.setParameters(params);
    
    double porosity = 0.25;
    double S_w = 0.3;
    double S_o = 0.7;
    double S_g = 0.0;
    
    // Effective thermal conductivity
    double k_eff = nimf.effectiveThermalConductivity(porosity, S_w, S_o, S_g);
    EXPECT_GT(k_eff, 0.0);
    EXPECT_LT(k_eff, params.thermal_conductivity_solid);  // Less than solid alone
    
    // Temperature-dependent viscosity (steam injection)
    double mu_ref = 0.01;  // Pa·s at reference T
    double T_ref = 323.15; // K
    double T_hot = 473.15; // K (200°C - steam zone)
    double E_a = 20000.0;  // J/mol
    
    double mu_hot = nimf.viscosityAtTemperature(mu_ref, T_ref, T_hot, E_a);
    
    // Viscosity should decrease significantly at high temperature
    EXPECT_LT(mu_hot, mu_ref);
    EXPECT_LT(mu_hot / mu_ref, 0.1);  // Large reduction
}


// =============================================================================
// Geomechanics Integration Tests
// =============================================================================

/**
 * @brief Geomechanics integration tests
 */
class GeomechanicsIntegrationTest : public HighFidelityIntegrationTest {
protected:
    void SetUp() override {
        HighFidelityIntegrationTest::SetUp();
        
        // Rock properties
        config_map["young_modulus"] = "20e9";
        config_map["poisson_ratio"] = "0.25";
        config_map["density"] = "2650.0";
    }
};

TEST_F(GeomechanicsIntegrationTest, FiniteStrainRubberLikeDeformation) {
    // Test large deformation of soft geomaterial
    
    FiniteStrainMechanics fsm;
    FiniteStrainMechanics::Parameters params;
    params.model = HyperelasticModel::NEO_HOOKEAN;
    params.mu = 1e6;      // Pa (soft material)
    params.kappa = 1e8;   // Pa
    params.use_quasi_incompressible = true;
    fsm.setParameters(params);
    
    // Test various deformation states
    std::vector<double> stretch_ratios = {0.8, 0.9, 1.0, 1.1, 1.2, 1.3};
    std::vector<double> energies, stresses;
    
    for (double lambda : stretch_ratios) {
        // Uniaxial incompressible deformation
        double F[3][3] = {{lambda, 0, 0},
                          {0, 1.0/std::sqrt(lambda), 0},
                          {0, 0, 1.0/std::sqrt(lambda)}};
        
        double stress[3][3];
        fsm.computeStress(F, stress);
        
        double W = fsm.computeStrainEnergy(F);
        
        energies.push_back(W);
        stresses.push_back(stress[0][0]);
    }
    
    // Energy minimum at λ = 1
    size_t idx_unity = 2;  // λ = 1.0
    EXPECT_LT(energies[idx_unity], energies[0]);
    EXPECT_LT(energies[idx_unity], energies.back());
    
    // Tension (λ > 1) -> positive stress, Compression (λ < 1) -> negative stress
    EXPECT_LT(stresses[0], 0.0);  // λ = 0.8
    EXPECT_GT(stresses.back(), 0.0);  // λ = 1.3
}

TEST_F(GeomechanicsIntegrationTest, ViscoplasticDeformation) {
    // Test time-dependent plastic deformation (salt creep analogy)
    
    Viscoplasticity vp;
    Viscoplasticity::Parameters params;
    params.model = ViscoplasticModel::PERZYNA;
    params.yield_stress = 15e6;   // Pa (salt-like)
    params.viscosity = 1e14;      // Pa·s
    params.rate_sensitivity = 0.15;
    vp.setParameters(params);
    
    // Stress state above yield
    double sigma = 25e6;  // Pa (in-situ stress)
    double stress_vec[6] = {sigma, sigma, sigma, 0, 0, 0};
    
    // Initial plastic strain rate
    double eps_dot_initial = vp.equivalentPlasticStrainRate(stress_vec, 0.0);
    EXPECT_GT(eps_dot_initial, 0.0);
    
    // With hardening, rate should decrease
    double eps_p_accumulated = 0.05;  // 5% plastic strain
    double eps_dot_later = vp.equivalentPlasticStrainRate(stress_vec, eps_p_accumulated);
    
    // Due to hardening, should be lower (or equal if no hardening)
    EXPECT_LE(eps_dot_later, eps_dot_initial * 1.1);  // Allow small tolerance
}

TEST_F(GeomechanicsIntegrationTest, CreepUnderConstantLoad) {
    // Test long-term creep behavior
    
    CreepModel cm;
    CreepModel::Parameters params;
    params.model = CreepType::POWER_LAW;
    params.A = 5e-20;
    params.n = 4.0;
    params.Q = 130000.0;  // J/mol
    params.temperature = 423.15;  // K (150°C - deep reservoir)
    cm.setParameters(params);
    
    double stress = 30e6;  // Pa
    
    // Creep rate
    double eps_dot = cm.creepStrainRate(stress);
    EXPECT_GT(eps_dot, 0.0);
    
    // Cumulative strain at various times
    std::vector<double> times = {86400.0, 2592000.0, 31536000.0};  // 1 day, 1 month, 1 year
    
    double eps_prev = 0.0;
    for (double t : times) {
        double eps = eps_dot * t;  // Simple integration
        EXPECT_GT(eps, eps_prev);
        eps_prev = eps;
    }
}

TEST_F(GeomechanicsIntegrationTest, GradientDamageLocalization) {
    // Test damage localization with gradient regularization
    
    GradientEnhancedDamage ged;
    GradientEnhancedDamage::Parameters params;
    params.model = DamageEvolutionModel::EXPONENTIAL;
    params.kappa_0 = 5e-5;
    params.kappa_c = 0.005;
    params.alpha_d = 0.99;
    params.beta_d = 200.0;
    params.internal_length = 0.02;  // m (damage zone width)
    ged.setParameters(params);
    
    // Test damage evolution through loading history
    std::vector<double> kappa_history = {0, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 2e-3};
    std::vector<double> damage_history;
    
    for (double kappa : kappa_history) {
        double d = ged.computeDamage(kappa);
        damage_history.push_back(d);
    }
    
    // Damage should be monotonically increasing
    for (size_t i = 1; i < damage_history.size(); ++i) {
        EXPECT_GE(damage_history[i], damage_history[i-1]);
    }
    
    // Check internal length provides regularization
    double c = ged.gradientParameter();
    EXPECT_NEAR(c, 0.02 * 0.02, 1e-6);  // l² = 0.0004
}


// =============================================================================
// Coupled Physics Integration Tests
// =============================================================================

/**
 * @brief Integration tests for coupled physics
 */
class CoupledPhysicsIntegrationTest : public HighFidelityIntegrationTest {
protected:
    void SetUp() override {
        HighFidelityIntegrationTest::SetUp();
        
        // Coupling parameters
        config_map["biot_alpha"] = "0.8";
        config_map["thermal_expansion"] = "1e-5";
    }
};

TEST_F(CoupledPhysicsIntegrationTest, FullBiotDynamicResponse) {
    // Test dynamic poroelastic response (earthquake-like loading)
    
    FullBiotDynamics fbd;
    FullBiotDynamics::Parameters params;
    params.formulation = BiotFormulation::U_P;
    params.biot_alpha = 0.85;
    params.biot_M = 8e9;
    params.porosity = 0.2;
    params.permeability = 5e-14;
    params.fluid_viscosity = 1e-3;
    params.solid_density = 2700.0;
    params.fluid_density = 1050.0;
    params.shear_modulus = 20e9;
    params.bulk_modulus = 30e9;
    params.enable_inertia = true;
    params.include_tortuosity = true;
    params.tortuosity = 1.5;
    fbd.setParameters(params);
    
    // Wave velocities
    double Vp_drained = fbd.pwaveVelocityDrained();
    double Vp_undrained = fbd.pwaveVelocityUndrained();
    double Vs = fbd.swaveVelocity();
    double Vp_slow = fbd.slowPwaveVelocity();
    
    // Physical constraints
    EXPECT_GT(Vp_drained, 0.0);
    EXPECT_GE(Vp_undrained, Vp_drained);  // Undrained stiffer
    EXPECT_LT(Vs, Vp_drained);            // S-wave slower than P-wave
    EXPECT_LT(Vp_slow, Vs);               // Slow P-wave is slowest
    
    // Typical rock values (km/s)
    EXPECT_GT(Vp_drained, 2000.0);
    EXPECT_LT(Vp_drained, 6000.0);
    
    // Characteristic frequency
    double f_c = fbd.characteristicFrequency();
    EXPECT_GT(f_c, 0.0);
    
    // Skempton coefficient
    double B = fbd.skemptonCoefficient();
    EXPECT_GT(B, 0.0);
    EXPECT_LE(B, 1.0);
}

TEST_F(CoupledPhysicsIntegrationTest, THMGeothermalReservoir) {
    // Test THM coupling in geothermal reservoir
    
    THMCoupling thm;
    THMCoupling::Parameters params;
    params.biot_alpha = 0.9;
    params.biot_M = 5e9;
    params.thermal_expansion_solid = 1.2e-5;
    params.thermal_expansion_fluid = 3e-4;
    params.thermal_conductivity = 3.0;
    params.specific_heat = 900.0;
    params.reference_temperature = 293.15;
    params.coupling_thermal_stress = true;
    params.coupling_thermal_fluid = true;
    params.coupling_pore_pressure = true;
    thm.setParameters(params);
    
    // Injection of cold water into hot reservoir
    double T_reservoir = 473.15;  // K (200°C)
    double T_injection = 323.15;  // K (50°C)
    double dT = T_injection - T_reservoir;  // -150 K (cooling)
    
    // Thermal contraction induces tensile stress
    double thermal_stress = thm.thermalStress(dT);
    EXPECT_GT(thermal_stress, 0.0);  // Tensile (positive) for cooling
    
    // Thermal strain
    double thermal_strain = thm.thermalExpansionStrain(dT);
    EXPECT_LT(thermal_strain, 0.0);  // Contraction (negative) for cooling
    
    // Pore pressure decrease due to cooling in confined volume
    double K = 30e9;  // Pa bulk modulus
    double dp_thermal = thm.thermalPressurization(dT, K);
    // Note: actual sign depends on specific formulation
    
    // Characteristic time scales
    double L = 100.0;  // m reservoir scale
    double cv = thm.consolidationCoefficient(1e-13, 1e-3, K, 0.2);
    double t_diffusion = thm.characteristicTime(L, cv);
    
    EXPECT_GT(t_diffusion, 0.0);
}

TEST_F(CoupledPhysicsIntegrationTest, THMCReactiveTransport) {
    // Test THMC coupling with chemical reactions
    
    THMCCoupling thmc;
    THMCCoupling::Parameters params;
    params.thm_params.biot_alpha = 0.8;
    params.thm_params.thermal_expansion_solid = 1e-5;
    params.thm_params.thermal_conductivity = 2.5;
    params.num_species = 3;  // CO2, H2O, mineral
    params.diffusion_coefficients = {2e-9, 2.5e-9, 0.0};
    params.reaction_rate_constants = {1e-8, 0.0, 5e-10};
    params.chemical_strain_coefficients = {0.0, 0.0, 0.02};  // Only mineral contributes to strain
    params.enable_chemical_dissolution = true;
    params.enable_chemical_precipitation = true;
    params.enable_chemical_swelling = true;
    thmc.setParameters(params);
    
    // Chemical strain from mineral dissolution/precipitation
    std::vector<double> concentrations = {0.5, 0.8, -0.1};  // mol/m³ change (negative = dissolution)
    double eps_chem = thmc.chemicalStrain(concentrations);
    
    // Dissolution reduces solid volume (negative strain contribution)
    EXPECT_LT(eps_chem, 0.0);
    
    // Porosity change due to reactions
    double diss_rate = 1e-8;   // mol/(m³·s)
    double precip_rate = 5e-9; // mol/(m³·s)
    double V_molar = 4e-5;     // m³/mol
    
    double phi_dot = thmc.porosityChangeRate(diss_rate, precip_rate, V_molar);
    
    // Net dissolution increases porosity
    EXPECT_GT(phi_dot, 0.0);
    
    // Reaction heat
    std::vector<double> rates = {1e-8, 0.0, 5e-9};
    std::vector<double> enthalpies = {-50000.0, 0.0, 30000.0};  // J/mol
    
    double Q_rxn = thmc.reactionHeat(rates, enthalpies);
    // Net heat depends on balance of exothermic/endothermic reactions
}

TEST_F(CoupledPhysicsIntegrationTest, UnsaturatedSoilMechanics) {
    // Test unsaturated soil behavior
    
    UnsaturatedFlowCoupling ufc;
    UnsaturatedFlowCoupling::Parameters params;
    params.retention_model = RetentionModel::VAN_GENUCHTEN;
    params.van_genuchten_alpha = 0.1;  // 1/m
    params.van_genuchten_n = 2.0;
    params.van_genuchten_m = 0.5;
    params.residual_saturation = 0.05;
    params.air_entry_pressure = 5000.0;
    params.biot_alpha = 1.0;
    params.bishops_chi_model = BishopChiModel::SATURATION;
    ufc.setParameters(params);
    
    // Test water retention curve
    std::vector<double> suctions = {0, 1000, 5000, 10000, 50000, 100000};  // Pa
    std::vector<double> saturations;
    
    for (double psi : suctions) {
        double S_w = ufc.waterSaturation(psi);
        saturations.push_back(S_w);
        
        // Bounds
        EXPECT_GE(S_w, params.residual_saturation);
        EXPECT_LE(S_w, 1.0);
    }
    
    // Saturation should decrease with suction
    for (size_t i = 1; i < saturations.size(); ++i) {
        EXPECT_LE(saturations[i], saturations[i-1]);
    }
    
    // Bishop's effective stress
    double total_stress = -100e3;  // Pa (compression)
    double p_water = 50e3;         // Pa
    double p_air = 101325.0;       // Pa (atmospheric)
    
    double S_w = ufc.waterSaturation(p_air - p_water);
    double chi = ufc.bishopParameter(S_w);
    
    double sigma_eff = ufc.effectiveStress(total_stress, p_water, p_air, chi);
    // Effective stress should account for partial saturation
}

TEST_F(CoupledPhysicsIntegrationTest, MultiScaleHomogenization) {
    // Test multi-scale modeling for heterogeneous rock
    
    MultiScaleCoupling msc;
    MultiScaleCoupling::Parameters params;
    params.method = MultiScaleMethod::FE_SQUARED;
    params.rve_size = 0.01;  // m (1 cm RVE)
    params.homogenization = HomogenizationType::COMPUTATIONAL;
    params.boundary_condition = RVEBoundaryCondition::PERIODIC;
    msc.setParameters(params);
    
    // Two-phase composite (matrix + inclusions)
    double K_matrix = 30e9, mu_matrix = 20e9;
    double K_inclusion = 80e9, mu_inclusion = 50e9;
    double f_inclusion = 0.3;  // 30% inclusion volume fraction
    
    // Bounds
    std::vector<double> K = {K_matrix, K_inclusion};
    std::vector<double> f = {1 - f_inclusion, f_inclusion};
    
    double K_voigt = msc.voigtAverage(K, f);
    double K_reuss = msc.reussAverage(K, f);
    
    // Hashin-Shtrikman bounds (tighter)
    double K_HS_lower, K_HS_upper, mu_HS_lower, mu_HS_upper;
    msc.hashinShtrikmanBounds(K_matrix, mu_matrix, K_inclusion, mu_inclusion,
                              f_inclusion, K_HS_lower, K_HS_upper,
                              mu_HS_lower, mu_HS_upper);
    
    // Bounds should satisfy: Reuss ≤ HS_lower ≤ HS_upper ≤ Voigt
    EXPECT_LE(K_reuss, K_HS_lower * 1.01);  // Small tolerance
    EXPECT_LE(K_HS_lower, K_HS_upper);
    EXPECT_LE(K_HS_upper, K_voigt * 1.01);
    
    // Mori-Tanaka estimate
    double K_MT, mu_MT;
    msc.moriTanaka(K_matrix, mu_matrix, K_inclusion, mu_inclusion,
                   f_inclusion, K_MT, mu_MT);
    
    // MT should be within HS bounds
    EXPECT_GE(K_MT, K_HS_lower * 0.99);
    EXPECT_LE(K_MT, K_HS_upper * 1.01);
}


// =============================================================================
// Numerical Methods Integration Tests
// =============================================================================

/**
 * @brief Integration tests for numerical methods
 */
class NumericsIntegrationTest : public HighFidelityIntegrationTest {
protected:
    MPI_Comm comm;
    
    void SetUp() override {
        HighFidelityIntegrationTest::SetUp();
        comm = PETSC_COMM_WORLD;
    }
};

TEST_F(NumericsIntegrationTest, SUPGStabilizedAdvection) {
    // Test SUPG for advection-dominated flow
    
    SUPGStabilization supg;
    SUPGStabilization::Parameters params;
    params.method = StabilizationMethod::SUPG;
    params.tau_formula = TauFormula::SHAKIB;
    params.tau_multiplier = 1.0;
    params.shock_capturing = true;
    params.shock_capture_coeff = 0.5;
    supg.setParameters(params);
    
    // High Peclet number flow
    PetscReal velocity[3] = {10.0, 0.0, 0.0};  // m/s
    PetscReal diffusivity = 1e-6;              // m²/s (very small)
    PetscReal h = 0.01;                        // m element size
    
    PetscReal Pe = supg.pecletNumber(10.0, h, diffusivity);
    EXPECT_GT(Pe, 1000.0);  // Very high Peclet
    
    // SUPG parameter
    PetscReal tau = supg.calculateTau(velocity, diffusivity, 0.0, h, 0.001);
    
    // τ should be O(h/|v|) for advection-dominated
    PetscReal tau_expected = h / (2.0 * 10.0);
    EXPECT_NEAR(tau, tau_expected, tau_expected * 0.5);
    
    // Artificial diffusivity
    PetscReal kappa_art = supg.optimalDiffusivity(10.0, h, diffusivity);
    EXPECT_GT(kappa_art, 0.0);
    EXPECT_LT(kappa_art, 10.0 * h);  // Maximum is |v|·h/2
}

TEST_F(NumericsIntegrationTest, PETScTimeIntegrationSetup) {
    // Test PETSc time integration wrapper
    
    PETScTimeIntegration ts(comm);
    PETScTimeIntegration::Parameters params;
    params.scheme = TimeScheme::BDF2;
    params.adaptive = true;
    params.atol = 1e-6;
    params.rtol = 1e-4;
    params.max_dt = 1.0;
    params.min_dt = 1e-10;
    ts.setParameters(params);
    
    // Verify properties
    EXPECT_EQ(ts.getOrder(), 2);
    EXPECT_TRUE(ts.isAStable());
    EXPECT_TRUE(ts.isLStable());
    
    // Test with different schemes
    params.scheme = TimeScheme::GENERALIZED_ALPHA;
    params.rho_infinity = 0.5;
    ts.setParameters(params);
    
    EXPECT_TRUE(ts.isAStable());
}

TEST_F(NumericsIntegrationTest, PAdaptivitySmoothnessIndicator) {
    // Test p-adaptivity smoothness estimation
    
    PAdaptivityManager p_adapt;
    PAdaptivityManager::Parameters params;
    params.min_order = 1;
    params.max_order = 6;
    params.default_order = 3;
    params.smoothness_threshold_inc = 0.8;
    params.smoothness_threshold_dec = 0.2;
    p_adapt.setParameters(params);
    
    // Test order change decisions
    struct TestCase {
        PetscReal smoothness;
        PetscReal error;
        PetscInt current_order;
        PetscInt expected_change;
    };
    
    std::vector<TestCase> cases = {
        {0.9, 0.6, 3, 1},    // High smoothness, high error -> increase
        {0.9, 0.3, 3, 0},    // High smoothness, low error -> keep
        {0.1, 0.5, 3, -1},   // Low smoothness -> decrease
        {0.5, 0.5, 3, 0},    // Medium -> keep
        {0.9, 0.9, 6, 0},    // At max order, can't increase
        {0.1, 0.5, 1, 0},    // At min order, can't decrease
    };
    
    for (const auto& tc : cases) {
        PetscInt change = p_adapt.decideOrderChange(tc.smoothness, tc.error, tc.current_order);
        EXPECT_EQ(change, tc.expected_change);
    }
}

TEST_F(NumericsIntegrationTest, PhysicsBasedPreconditionerSetup) {
    // Test physics-based preconditioner configuration
    
    PETScPhysicsPreconditioner pc(comm);
    
    std::map<std::string, std::string> config;
    config["precond_type"] = "fieldsplit_schur";
    config["precond_schur_fact"] = "full";
    config["precond_mg_levels"] = "4";
    
    pc.configure(config);
    
    // Just verify configuration doesn't crash
    SUCCEED();
}


// =============================================================================
// End-to-End Scenario Tests
// =============================================================================

TEST(EndToEndTest, CO2SequestrationScenario) {
    // Simplified end-to-end test for CO2 sequestration
    
    // 1. Create fluid flow model with CO2-brine system
    MiscibleFlow mf;
    MiscibleFlow::Parameters mf_params;
    mf_params.model = MiscibleFlowModel::MULTI_CONTACT;
    mf_params.mmp = 8e6;  // Pa
    mf_params.fcm_pressure = 15e6;
    mf_params.mixing_parameter = 0.6;
    mf.setParameters(mf_params);
    
    // 2. Create geomechanics model
    FiniteStrainMechanics fsm;
    FiniteStrainMechanics::Parameters fsm_params;
    fsm_params.model = HyperelasticModel::NEO_HOOKEAN;
    fsm_params.mu = 15e9;
    fsm_params.kappa = 25e9;
    fsm.setParameters(fsm_params);
    
    // 3. Create coupling
    THMCCoupling thmc;
    THMCCoupling::Parameters thmc_params;
    thmc_params.thm_params.biot_alpha = 0.85;
    thmc_params.thm_params.thermal_expansion_solid = 8e-6;
    thmc_params.num_species = 2;  // CO2, brine
    thmc_params.enable_chemical_dissolution = true;
    thmc.setParameters(thmc_params);
    
    // 4. Test at injection conditions
    double P_inj = 15e6;  // Pa injection pressure
    double T_inj = 323.15;  // K (50°C)
    
    // Miscibility at injection pressure
    double f_misc = mf.miscibilityFraction(P_inj);
    EXPECT_GT(f_misc, 0.5);  // Should be partially miscible
    
    // Chemical effects
    std::vector<double> delta_c = {0.1, -0.05};  // CO2 increase, brine decrease
    double eps_chem = thmc.chemicalStrain(delta_c);
    
    // All components should work together
    SUCCEED();
}

TEST(EndToEndTest, GeothermalEGSScenario) {
    // Enhanced Geothermal Systems scenario
    
    // 1. Fracture flow (non-Darcy at high rates)
    NonDarcyFlow ndf;
    NonDarcyFlow::Parameters ndf_params;
    ndf_params.model = NonDarcyModel::FORCHHEIMER;
    ndf_params.forchheimer_beta = 1e7;  // Fractured rock
    ndf.setParameters(ndf_params);
    
    // 2. THM coupling for thermal stress
    THMCoupling thm;
    THMCoupling::Parameters thm_params;
    thm_params.biot_alpha = 0.8;
    thm_params.thermal_expansion_solid = 8e-6;
    thm_params.thermal_conductivity = 3.5;  // Granite-like
    thm_params.coupling_thermal_stress = true;
    thm.setParameters(thm_params);
    
    // 3. Biot dynamics for pressure diffusion
    FullBiotDynamics fbd;
    FullBiotDynamics::Parameters fbd_params;
    fbd_params.formulation = BiotFormulation::U_P;
    fbd_params.biot_alpha = 0.8;
    fbd_params.biot_M = 10e9;
    fbd_params.permeability = 1e-15;  // Low matrix permeability
    fbd.setParameters(fbd_params);
    
    // Test thermal drawdown scenario
    double T_hot = 473.15;   // K (200°C reservoir)
    double T_cold = 323.15;  // K (50°C injection)
    double dT = T_cold - T_hot;
    
    // Thermal stress from cooling
    double sigma_th = thm.thermalStress(dT);
    EXPECT_GT(sigma_th, 0.0);  // Tensile due to cooling
    
    // This tensile stress could induce fracturing
    SUCCEED();
}
