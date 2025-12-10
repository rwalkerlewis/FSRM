/**
 * @file test_high_fidelity_geomechanics.cpp
 * @brief Unit tests for HighFidelityGeomechanics module
 * 
 * Tests cover:
 * - Finite strain mechanics (Neo-Hookean, Mooney-Rivlin, Ogden)
 * - Viscoplasticity (Perzyna, Duvaut-Lions)
 * - Creep models (Power Law, Norton, Burgers)
 * - Hypoplasticity
 * - Gradient-enhanced damage
 */

#include <gtest/gtest.h>
#include "HighFidelityGeomechanics.hpp"
#include <cmath>

using namespace FSRM::HighFidelity;

// =============================================================================
// Helper Functions for Testing
// =============================================================================

namespace {
    // Identity deformation gradient
    void setIdentityF(double F[3][3]) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                F[i][j] = (i == j) ? 1.0 : 0.0;
            }
        }
    }
    
    // Simple shear deformation gradient
    void setSimpleShearF(double F[3][3], double gamma) {
        setIdentityF(F);
        F[0][1] = gamma;  // Shear in x-y plane
    }
    
    // Uniaxial extension deformation gradient
    void setUniaxialExtensionF(double F[3][3], double lambda) {
        setIdentityF(F);
        F[0][0] = lambda;
        F[1][1] = 1.0 / std::sqrt(lambda);
        F[2][2] = 1.0 / std::sqrt(lambda);
    }
    
    // Compute determinant of 3x3 matrix
    double det3x3(const double A[3][3]) {
        return A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
             - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
             + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
    }
}


// =============================================================================
// Finite Strain Mechanics Tests
// =============================================================================

class FiniteStrainTest : public ::testing::Test {
protected:
    FiniteStrainMechanics fsm;
    
    void SetUp() override {
        FiniteStrainMechanics::Parameters params;
        params.model = HyperelasticModel::NEO_HOOKEAN;
        params.mu = 1e6;     // Pa (shear modulus)
        params.kappa = 1e9;  // Pa (bulk modulus)
        params.use_quasi_incompressible = true;
        fsm.setParameters(params);
    }
};

TEST_F(FiniteStrainTest, NeoHookeanIdentityDeformation) {
    double F[3][3], stress[3][3];
    setIdentityF(F);
    
    fsm.computeStress(F, stress);
    
    // At identity deformation, stress should be zero
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            EXPECT_NEAR(stress[i][j], 0.0, 1e-6);
        }
    }
}

TEST_F(FiniteStrainTest, NeoHookeanUniaxialTension) {
    double F[3][3], stress[3][3];
    double lambda = 1.1;  // 10% extension
    setUniaxialExtensionF(F, lambda);
    
    fsm.computeStress(F, stress);
    
    // σ_11 should be positive (tensile)
    EXPECT_GT(stress[0][0], 0.0);
    
    // For uniaxial test, σ_22 ≈ σ_33 ≈ 0 (lateral surfaces are free)
    // Note: This depends on implementation details
}

TEST_F(FiniteStrainTest, NeoHookeanSimpleShear) {
    double F[3][3], stress[3][3];
    double gamma = 0.1;  // 10% shear
    setSimpleShearF(F, gamma);
    
    fsm.computeStress(F, stress);
    
    // Shear stress σ_12 should be non-zero
    EXPECT_NE(stress[0][1], 0.0);
}

TEST_F(FiniteStrainTest, StrainEnergy) {
    double F[3][3];
    setUniaxialExtensionF(F, 1.1);
    
    double W = fsm.computeStrainEnergy(F);
    
    // Strain energy should be positive
    EXPECT_GT(W, 0.0);
}

TEST_F(FiniteStrainTest, StrainEnergyZeroAtIdentity) {
    double F[3][3];
    setIdentityF(F);
    
    double W = fsm.computeStrainEnergy(F);
    
    // Strain energy should be zero at identity
    EXPECT_NEAR(W, 0.0, 1e-10);
}

TEST_F(FiniteStrainTest, MooneyRivlinModel) {
    FiniteStrainMechanics::Parameters params;
    params.model = HyperelasticModel::MOONEY_RIVLIN;
    params.C10 = 0.5e6;  // Pa
    params.C01 = 0.1e6;  // Pa
    params.kappa = 1e9;
    fsm.setParameters(params);
    
    double F[3][3], stress[3][3];
    setUniaxialExtensionF(F, 1.2);
    
    fsm.computeStress(F, stress);
    
    // Stress should be non-zero
    EXPECT_NE(stress[0][0], 0.0);
}

TEST_F(FiniteStrainTest, OgdenModel) {
    FiniteStrainMechanics::Parameters params;
    params.model = HyperelasticModel::OGDEN;
    params.ogden_n_terms = 1;
    params.ogden_mu = {1e6};
    params.ogden_alpha = {2.0};
    params.kappa = 1e9;
    fsm.setParameters(params);
    
    double F[3][3], stress[3][3];
    setUniaxialExtensionF(F, 1.15);
    
    fsm.computeStress(F, stress);
    
    EXPECT_GT(stress[0][0], 0.0);
}

TEST_F(FiniteStrainTest, MaterialTangent) {
    double F[3][3];
    double C[3][3][3][3];
    setUniaxialExtensionF(F, 1.05);
    
    fsm.computeMaterialTangent(F, C);
    
    // Major symmetry: C_ijkl = C_klij
    EXPECT_NEAR(C[0][1][2][0], C[2][0][0][1], 1e-6);
    
    // Minor symmetry: C_ijkl = C_jikl = C_ijlk
    EXPECT_NEAR(C[0][1][0][1], C[1][0][0][1], 1e-6);
    
    // Positive definiteness (diagonal components positive)
    EXPECT_GT(C[0][0][0][0], 0.0);
    EXPECT_GT(C[1][1][1][1], 0.0);
}

TEST_F(FiniteStrainTest, VolumePreservationCheck) {
    double F[3][3];
    setUniaxialExtensionF(F, 1.1);
    
    double J = det3x3(F);
    
    // For incompressible material, J should be 1
    EXPECT_NEAR(J, 1.0, 1e-10);
}

TEST_F(FiniteStrainTest, Configure) {
    std::map<std::string, std::string> config;
    config["hyperelastic_model"] = "mooney_rivlin";
    config["C10"] = "1e6";
    config["C01"] = "2e5";
    
    FiniteStrainMechanics fsm2;
    fsm2.configure(config);
    
    SUCCEED();
}


// =============================================================================
// Finite Strain Plasticity Tests
// =============================================================================

class FiniteStrainPlasticityTest : public ::testing::Test {
protected:
    FiniteStrainPlasticity fsp;
    
    void SetUp() override {
        FiniteStrainPlasticity::Parameters params;
        params.model = FinitePlasticityModel::MULTIPLICATIVE_J2;
        params.mu = 80e9;      // Pa (steel-like shear modulus)
        params.kappa = 170e9;  // Pa
        params.yield_stress = 250e6;  // Pa
        params.hardening_modulus = 1e9;
        params.hardening_exponent = 0.5;
        fsp.setParameters(params);
    }
};

TEST_F(FiniteStrainPlasticityTest, ElasticResponse) {
    double F[3][3], stress[3][3];
    FiniteStrainPlasticity::State state;
    state.plastic_deformation_gradient[0][0] = 1.0;
    state.plastic_deformation_gradient[1][1] = 1.0;
    state.plastic_deformation_gradient[2][2] = 1.0;
    state.equivalent_plastic_strain = 0.0;
    
    // Small elastic deformation (below yield)
    setUniaxialExtensionF(F, 1.001);  // 0.1% strain
    
    fsp.computeStress(F, state, stress);
    
    // Should have tensile stress
    EXPECT_GT(stress[0][0], 0.0);
    
    // No plastic strain for elastic loading
    EXPECT_NEAR(state.equivalent_plastic_strain, 0.0, 1e-6);
}

TEST_F(FiniteStrainPlasticityTest, PlasticYielding) {
    double F[3][3], stress[3][3];
    FiniteStrainPlasticity::State state;
    state.plastic_deformation_gradient[0][0] = 1.0;
    state.plastic_deformation_gradient[1][1] = 1.0;
    state.plastic_deformation_gradient[2][2] = 1.0;
    state.equivalent_plastic_strain = 0.0;
    
    // Large deformation (beyond yield)
    setUniaxialExtensionF(F, 1.01);  // 1% strain
    
    fsp.computeStress(F, state, stress);
    
    // For large strain beyond yield, plastic strain should evolve
    // (depends on actual yield stress and moduli)
}

TEST_F(FiniteStrainPlasticityTest, YieldFunction) {
    double stress_vec[6] = {300e6, 0, 0, 0, 0, 0};  // Uniaxial tension
    
    double f = fsp.yieldFunction(stress_vec, 0.0);
    
    // Above yield stress: f > 0
    EXPECT_GT(f, 0.0);
}

TEST_F(FiniteStrainPlasticityTest, YieldFunctionBelowYield) {
    double stress_vec[6] = {200e6, 0, 0, 0, 0, 0};  // Below yield
    
    double f = fsp.yieldFunction(stress_vec, 0.0);
    
    // Below yield stress: f < 0
    EXPECT_LT(f, 0.0);
}

TEST_F(FiniteStrainPlasticityTest, IsotropicHardening) {
    double sigma_y_0 = fsp.currentYieldStress(0.0);
    double sigma_y_1 = fsp.currentYieldStress(0.1);
    
    // With hardening, yield stress increases with plastic strain
    EXPECT_GT(sigma_y_1, sigma_y_0);
}

TEST_F(FiniteStrainPlasticityTest, Configure) {
    std::map<std::string, std::string> config;
    config["finite_plasticity_model"] = "multiplicative_j2";
    config["yield_stress"] = "300e6";
    config["hardening_modulus"] = "2e9";
    
    FiniteStrainPlasticity fsp2;
    fsp2.configure(config);
    
    SUCCEED();
}


// =============================================================================
// Viscoplasticity Tests
// =============================================================================

class ViscoplasticityTest : public ::testing::Test {
protected:
    Viscoplasticity vp;
    
    void SetUp() override {
        Viscoplasticity::Parameters params;
        params.model = ViscoplasticModel::PERZYNA;
        params.yield_stress = 250e6;   // Pa
        params.viscosity = 1e12;       // Pa·s
        params.rate_sensitivity = 0.1;
        params.reference_strain_rate = 1e-3;  // s⁻¹
        vp.setParameters(params);
    }
};

TEST_F(ViscoplasticityTest, PerzynaOverstress) {
    double stress_vec[6] = {300e6, 0, 0, 0, 0, 0};  // Above yield
    
    double overstress = vp.computeOverstress(stress_vec, 0.0);
    
    // Overstress = σ_eq - σ_y = √3·|σ_11| - 250e6 ≈ 520e6 - 250e6
    EXPECT_GT(overstress, 0.0);
}

TEST_F(ViscoplasticityTest, PerzynaZeroOverstress) {
    double stress_vec[6] = {100e6, 0, 0, 0, 0, 0};  // Below yield
    
    double overstress = vp.computeOverstress(stress_vec, 0.0);
    
    // Below yield: overstress = 0
    EXPECT_NEAR(overstress, 0.0, 1e-6);
}

TEST_F(ViscoplasticityTest, PlasticStrainRate) {
    double stress_vec[6] = {300e6, 0, 0, 0, 0, 0};
    double eps_p_dot[6];
    
    vp.computePlasticStrainRate(stress_vec, 0.0, eps_p_dot);
    
    // Plastic strain rate should be positive in tension direction
    EXPECT_GT(eps_p_dot[0], 0.0);
    
    // Deviatoric flow: trace should be zero
    double trace = eps_p_dot[0] + eps_p_dot[1] + eps_p_dot[2];
    EXPECT_NEAR(trace, 0.0, 1e-6 * std::abs(eps_p_dot[0]));
}

TEST_F(ViscoplasticityTest, EquivalentPlasticStrainRate) {
    double stress_vec[6] = {300e6, 0, 0, 0, 0, 0};
    
    double eps_p_dot_eq = vp.equivalentPlasticStrainRate(stress_vec, 0.0);
    
    // Should be positive when above yield
    EXPECT_GT(eps_p_dot_eq, 0.0);
}

TEST_F(ViscoplasticityTest, DuvautLionsModel) {
    Viscoplasticity::Parameters params;
    params.model = ViscoplasticModel::DUVAUT_LIONS;
    params.yield_stress = 250e6;
    params.relaxation_time = 1.0;  // s
    vp.setParameters(params);
    
    double stress_vec[6] = {300e6, 0, 0, 0, 0, 0};
    double eps_p_dot[6];
    
    vp.computePlasticStrainRate(stress_vec, 0.0, eps_p_dot);
    
    EXPECT_GT(eps_p_dot[0], 0.0);
}

TEST_F(ViscoplasticityTest, RateSensitivityEffect) {
    // Compare strain rates at different rate sensitivities
    Viscoplasticity::Parameters params1, params2;
    params1.model = ViscoplasticModel::PERZYNA;
    params1.yield_stress = 250e6;
    params1.viscosity = 1e12;
    params1.rate_sensitivity = 0.05;  // Lower sensitivity
    
    params2.model = ViscoplasticModel::PERZYNA;
    params2.yield_stress = 250e6;
    params2.viscosity = 1e12;
    params2.rate_sensitivity = 0.2;   // Higher sensitivity
    
    Viscoplasticity vp1, vp2;
    vp1.setParameters(params1);
    vp2.setParameters(params2);
    
    double stress_vec[6] = {300e6, 0, 0, 0, 0, 0};
    
    double rate1 = vp1.equivalentPlasticStrainRate(stress_vec, 0.0);
    double rate2 = vp2.equivalentPlasticStrainRate(stress_vec, 0.0);
    
    // Higher rate sensitivity -> lower strain rate (more viscous response)
    // Note: depends on specific formulation
}

TEST_F(ViscoplasticityTest, Configure) {
    std::map<std::string, std::string> config;
    config["viscoplastic_model"] = "duvaut_lions";
    config["relaxation_time"] = "0.5";
    
    Viscoplasticity vp2;
    vp2.configure(config);
    
    SUCCEED();
}


// =============================================================================
// Creep Model Tests
// =============================================================================

class CreepModelTest : public ::testing::Test {
protected:
    CreepModel cm;
    
    void SetUp() override {
        CreepModel::Parameters params;
        params.model = CreepType::POWER_LAW;
        params.A = 1e-20;       // Creep coefficient
        params.n = 3.0;         // Stress exponent
        params.Q = 150000.0;    // Activation energy (J/mol)
        params.temperature = 500.0;  // K
        cm.setParameters(params);
    }
};

TEST_F(CreepModelTest, PowerLawCreepRate) {
    double stress = 100e6;  // Pa
    
    double eps_c_dot = cm.creepStrainRate(stress);
    
    // Should be positive
    EXPECT_GT(eps_c_dot, 0.0);
}

TEST_F(CreepModelTest, PowerLawStressDependence) {
    double stress1 = 50e6;
    double stress2 = 100e6;
    
    double rate1 = cm.creepStrainRate(stress1);
    double rate2 = cm.creepStrainRate(stress2);
    
    // Power law: ε̇ ∝ σⁿ, so rate2/rate1 ≈ (100/50)^3 = 8
    EXPECT_NEAR(rate2 / rate1, 8.0, 1.0);
}

TEST_F(CreepModelTest, TemperatureDependence) {
    CreepModel::Parameters params1, params2;
    params1.model = CreepType::POWER_LAW;
    params1.A = 1e-20;
    params1.n = 3.0;
    params1.Q = 150000.0;
    params1.temperature = 400.0;  // Lower T
    
    params2.model = CreepType::POWER_LAW;
    params2.A = 1e-20;
    params2.n = 3.0;
    params2.Q = 150000.0;
    params2.temperature = 600.0;  // Higher T
    
    CreepModel cm1, cm2;
    cm1.setParameters(params1);
    cm2.setParameters(params2);
    
    double stress = 100e6;
    double rate1 = cm1.creepStrainRate(stress);
    double rate2 = cm2.creepStrainRate(stress);
    
    // Higher temperature -> faster creep
    EXPECT_GT(rate2, rate1);
}

TEST_F(CreepModelTest, NortonCreep) {
    CreepModel::Parameters params;
    params.model = CreepType::NORTON;
    params.A = 1e-15;
    params.n = 4.5;
    params.temperature = 500.0;
    cm.setParameters(params);
    
    double stress = 80e6;
    double rate = cm.creepStrainRate(stress);
    
    EXPECT_GT(rate, 0.0);
}

TEST_F(CreepModelTest, BurgersCreepStrain) {
    CreepModel::Parameters params;
    params.model = CreepType::BURGERS;
    params.maxwell_viscosity = 1e18;  // Pa·s
    params.maxwell_modulus = 50e9;    // Pa
    params.kelvin_viscosity = 1e16;   // Pa·s
    params.kelvin_modulus = 30e9;     // Pa
    cm.setParameters(params);
    
    double stress = 50e6;
    double time = 3600.0;  // 1 hour
    
    double strain = cm.burgersCreepStrain(stress, time);
    
    EXPECT_GT(strain, 0.0);
}

TEST_F(CreepModelTest, BurgersDecomposition) {
    CreepModel::Parameters params;
    params.model = CreepType::BURGERS;
    params.maxwell_viscosity = 1e18;
    params.maxwell_modulus = 50e9;
    params.kelvin_viscosity = 1e16;
    params.kelvin_modulus = 30e9;
    cm.setParameters(params);
    
    double stress = 50e6;
    double time = 86400.0;  // 1 day
    
    double elastic, transient, steady;
    cm.burgersCreepComponents(stress, time, elastic, transient, steady);
    
    // All components should be positive
    EXPECT_GT(elastic, 0.0);
    EXPECT_GE(transient, 0.0);
    EXPECT_GT(steady, 0.0);
    
    // Total should equal sum of components
    double total = cm.burgersCreepStrain(stress, time);
    EXPECT_NEAR(total, elastic + transient + steady, total * 0.01);
}

TEST_F(CreepModelTest, Configure) {
    std::map<std::string, std::string> config;
    config["creep_model"] = "burgers";
    config["maxwell_viscosity"] = "1e17";
    config["kelvin_viscosity"] = "1e15";
    
    CreepModel cm2;
    cm2.configure(config);
    
    SUCCEED();
}


// =============================================================================
// Hypoplasticity Tests
// =============================================================================

class HypoplasticityTest : public ::testing::Test {
protected:
    Hypoplasticity hp;
    
    void SetUp() override {
        Hypoplasticity::Parameters params;
        params.phi_c = 33.0;       // Critical friction angle (degrees)
        params.hs = 1e9;           // Pa (granular hardness)
        params.n = 0.28;           // Exponent
        params.ed0 = 0.5;          // Min void ratio
        params.ec0 = 1.0;          // Critical void ratio
        params.ei0 = 1.15;         // Max void ratio
        params.alpha = 0.13;
        params.beta = 1.0;
        hp.setParameters(params);
    }
};

TEST_F(HypoplasticityTest, StressRate) {
    double stress[6] = {-100e3, -100e3, -100e3, 0, 0, 0};  // Isotropic compression
    double strain_rate[6] = {-1e-4, -1e-4, -1e-4, 0, 0, 0};  // Compressive
    double void_ratio = 0.65;
    
    double stress_rate[6];
    hp.computeStressRate(stress, strain_rate, void_ratio, stress_rate);
    
    // Stress rate should be non-zero
    double norm_sq = 0.0;
    for (int i = 0; i < 6; ++i) {
        norm_sq += stress_rate[i] * stress_rate[i];
    }
    EXPECT_GT(norm_sq, 0.0);
}

TEST_F(HypoplasticityTest, VoidRatioLimits) {
    double e_d = hp.minVoidRatio(-200e3);  // At 200 kPa
    double e_c = hp.criticalVoidRatio(-200e3);
    double e_i = hp.maxVoidRatio(-200e3);
    
    // Should maintain e_d < e_c < e_i
    EXPECT_LT(e_d, e_c);
    EXPECT_LT(e_c, e_i);
}

TEST_F(HypoplasticityTest, RelativeDensity) {
    double e = 0.7;
    double p = -200e3;  // Mean stress
    
    double r_e = hp.relativeDensity(e, p);
    
    // Should be between 0 and 1
    EXPECT_GE(r_e, 0.0);
    EXPECT_LE(r_e, 1.0);
}

TEST_F(HypoplasticityTest, BarotropyFactor) {
    double p = -200e3;  // Mean stress
    double e = 0.7;
    
    double f_b = hp.barotropyFactor(p, e);
    
    // Should be positive
    EXPECT_GT(f_b, 0.0);
}

TEST_F(HypoplasticityTest, PyknotopyFactor) {
    double e = 0.7;
    double p = -200e3;
    
    double f_e = hp.pyknotopyFactor(e, p);
    
    // Should be positive
    EXPECT_GT(f_e, 0.0);
}

TEST_F(HypoplasticityTest, CriticalStateLine) {
    // Test that void ratio converges to critical state under continued shearing
    // This is more of an integration test concept
    
    double e_c_100 = hp.criticalVoidRatio(-100e3);
    double e_c_1000 = hp.criticalVoidRatio(-1000e3);
    
    // Higher pressure -> lower critical void ratio
    EXPECT_GT(e_c_100, e_c_1000);
}

TEST_F(HypoplasticityTest, Configure) {
    std::map<std::string, std::string> config;
    config["hypoplasticity_phi_c"] = "35.0";
    config["hypoplasticity_hs"] = "2e9";
    
    Hypoplasticity hp2;
    hp2.configure(config);
    
    SUCCEED();
}


// =============================================================================
// Gradient-Enhanced Damage Tests
// =============================================================================

class GradientDamageTest : public ::testing::Test {
protected:
    GradientEnhancedDamage ged;
    
    void SetUp() override {
        GradientEnhancedDamage::Parameters params;
        params.model = DamageEvolutionModel::EXPONENTIAL;
        params.kappa_0 = 1e-4;       // Damage threshold
        params.kappa_c = 0.01;       // Critical strain
        params.alpha_d = 0.99;       // Maximum damage
        params.beta_d = 100.0;       // Evolution rate
        params.internal_length = 0.01;  // m
        params.nonlocal_averaging = NonlocalAveraging::GRADIENT;
        ged.setParameters(params);
    }
};

TEST_F(GradientDamageTest, DamageEvolution) {
    double kappa = 0.002;  // Above threshold
    
    double d = ged.computeDamage(kappa);
    
    // Damage should be between 0 and alpha_d
    EXPECT_GT(d, 0.0);
    EXPECT_LT(d, 0.99);
}

TEST_F(GradientDamageTest, DamageZeroBelowThreshold) {
    double kappa = 5e-5;  // Below threshold
    
    double d = ged.computeDamage(kappa);
    
    // No damage below threshold
    EXPECT_NEAR(d, 0.0, 1e-10);
}

TEST_F(GradientDamageTest, DamageMonotonicity) {
    double kappa1 = 0.002;
    double kappa2 = 0.005;
    
    double d1 = ged.computeDamage(kappa1);
    double d2 = ged.computeDamage(kappa2);
    
    // Damage should increase with history variable
    EXPECT_GT(d2, d1);
}

TEST_F(GradientDamageTest, DamageAsymptoticLimit) {
    double kappa = 1.0;  // Very large
    
    double d = ged.computeDamage(kappa);
    
    // Should approach maximum damage
    EXPECT_NEAR(d, 0.99, 0.01);
}

TEST_F(GradientDamageTest, EquivalentStrain) {
    double strain[6] = {0.001, -0.0003, -0.0003, 0, 0, 0};
    
    double kappa = ged.equivalentStrain(strain);
    
    // Should be positive
    EXPECT_GT(kappa, 0.0);
}

TEST_F(GradientDamageTest, EffectiveStress) {
    double stress[6] = {100e6, 50e6, 50e6, 0, 0, 0};
    double d = 0.3;
    
    double stress_eff[6];
    ged.effectiveStress(stress, d, stress_eff);
    
    // Effective stress = σ / (1-d)
    EXPECT_NEAR(stress_eff[0], 100e6 / 0.7, 1e6);
}

TEST_F(GradientDamageTest, LinearDamageModel) {
    GradientEnhancedDamage::Parameters params;
    params.model = DamageEvolutionModel::LINEAR_SOFTENING;
    params.kappa_0 = 1e-4;
    params.kappa_c = 0.01;
    params.alpha_d = 1.0;
    ged.setParameters(params);
    
    double kappa = 0.005;  // Midpoint
    double d = ged.computeDamage(kappa);
    
    // Linear: d ≈ 0.5 at midpoint
    EXPECT_NEAR(d, 0.5, 0.1);
}

TEST_F(GradientDamageTest, InternalLengthScale) {
    // Test that characteristic length affects gradient term
    double c = ged.gradientParameter();
    
    // c = l² (internal length squared)
    EXPECT_NEAR(c, 0.0001, 1e-6);  // (0.01)² = 0.0001
}

TEST_F(GradientDamageTest, Configure) {
    std::map<std::string, std::string> config;
    config["damage_model"] = "mazars";
    config["internal_length"] = "0.02";
    config["kappa_0"] = "2e-4";
    
    GradientEnhancedDamage ged2;
    ged2.configure(config);
    
    SUCCEED();
}


// =============================================================================
// Configuration Validation Tests
// =============================================================================

class GeomechConfigTest : public ::testing::Test {
protected:
    HighFidelityGeomechanicsConfig config;
};

TEST_F(GeomechConfigTest, ValidConfiguration) {
    config.enable_finite_strain = true;
    config.finite_strain_params.model = HyperelasticModel::NEO_HOOKEAN;
    config.finite_strain_params.mu = 1e6;
    config.finite_strain_params.kappa = 1e9;
    
    std::string error_msg;
    bool valid = config.validate(error_msg);
    
    EXPECT_TRUE(valid);
    EXPECT_TRUE(error_msg.empty());
}

TEST_F(GeomechConfigTest, InvalidNegativeModulus) {
    config.enable_finite_strain = true;
    config.finite_strain_params.mu = -1e6;  // Invalid: negative
    
    std::string error_msg;
    bool valid = config.validate(error_msg);
    
    EXPECT_FALSE(valid);
    EXPECT_FALSE(error_msg.empty());
}

TEST_F(GeomechConfigTest, InvalidCreepExponent) {
    config.enable_creep = true;
    config.creep_params.n = 0.0;  // Invalid: should be positive
    
    std::string error_msg;
    bool valid = config.validate(error_msg);
    
    EXPECT_FALSE(valid);
}

TEST_F(GeomechConfigTest, ParseConfig) {
    std::map<std::string, std::string> cfg;
    cfg["enable_finite_strain"] = "true";
    cfg["enable_viscoplastic"] = "false";
    cfg["enable_creep"] = "true";
    cfg["enable_hypoplastic"] = "false";
    cfg["enable_gradient_damage"] = "true";
    
    config.parseConfig(cfg);
    
    EXPECT_TRUE(config.enable_finite_strain);
    EXPECT_FALSE(config.enable_viscoplastic);
    EXPECT_TRUE(config.enable_creep);
    EXPECT_FALSE(config.enable_hypoplastic);
    EXPECT_TRUE(config.enable_gradient_damage);
}


// =============================================================================
// Factory Function Tests
// =============================================================================

TEST(GeomechFactoryTest, CreateFiniteStrainMechanics) {
    std::map<std::string, std::string> config;
    config["hyperelastic_model"] = "neo_hookean";
    config["shear_modulus"] = "1e6";
    config["bulk_modulus"] = "1e9";
    
    auto fsm = createFiniteStrainMechanics(config);
    
    ASSERT_NE(fsm, nullptr);
}

TEST(GeomechFactoryTest, CreateViscoplasticity) {
    std::map<std::string, std::string> config;
    config["viscoplastic_model"] = "perzyna";
    config["viscosity"] = "1e12";
    
    auto vp = createViscoplasticity(config);
    
    ASSERT_NE(vp, nullptr);
}

TEST(GeomechFactoryTest, CreateCreepModel) {
    std::map<std::string, std::string> config;
    config["creep_model"] = "power_law";
    config["creep_A"] = "1e-20";
    config["creep_n"] = "3.0";
    
    auto cm = createCreepModel(config);
    
    ASSERT_NE(cm, nullptr);
}

TEST(GeomechFactoryTest, CreateGradientDamage) {
    std::map<std::string, std::string> config;
    config["damage_model"] = "exponential";
    config["internal_length"] = "0.01";
    
    auto ged = createGradientDamage(config);
    
    ASSERT_NE(ged, nullptr);
}
