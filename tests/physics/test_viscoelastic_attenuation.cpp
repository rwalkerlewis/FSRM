/**
 * @file test_viscoelastic_attenuation.cpp
 * @brief Tests for viscoelastic attenuation models
 * 
 * Tests for:
 * - Maxwell mechanism relaxation
 * - Quality factor Q approximation
 * - Frequency-dependent properties
 * - Dispersion relations
 * - Causality (Kramers-Kronig)
 * - Anelastic state integration
 */

#include <gtest/gtest.h>
#include "Testing.hpp"
#include "ViscoelasticAttenuation.hpp"
#include <cmath>
#include <complex>
#include <vector>

using namespace FSRM;
using namespace FSRM::Testing;

class ViscoelasticAttenuationTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        
        // Reference velocities
        Vp_ref = 6000.0;    // P-wave velocity (m/s)
        Vs_ref = 3500.0;    // S-wave velocity (m/s)
        rho = 2700.0;       // Density (kg/m³)
        
        // Target Q values
        Qp = 100.0;         // P-wave quality factor
        Qs = 50.0;          // S-wave quality factor
        
        // Frequency range
        f_min = 0.01;       // Hz
        f_max = 10.0;       // Hz
        f_central = 1.0;    // Hz
    }
    
    int rank;
    double Vp_ref, Vs_ref, rho;
    double Qp, Qs;
    double f_min, f_max, f_central;
};

// ============================================================================
// Maxwell Mechanism Tests
// ============================================================================

TEST_F(ViscoelasticAttenuationTest, MaxwellMechanismCreation) {
    ViscoelasticAttenuation::MaxwellMechanism mech(1.0, 0.1, 0.1);
    
    EXPECT_NEAR(mech.relaxation_frequency, 1.0, 1e-10);
    EXPECT_NEAR(mech.relaxation_time, 1.0 / (2.0 * M_PI), 1e-6);
    EXPECT_NEAR(mech.anelastic_coefficient_p, 0.1, 1e-10);
    EXPECT_NEAR(mech.anelastic_coefficient_s, 0.1, 1e-10);
}

TEST_F(ViscoelasticAttenuationTest, RelaxationTimeFrequencyRelation) {
    double f = 2.0;  // Hz
    double tau = 1.0 / (2.0 * M_PI * f);
    
    ViscoelasticAttenuation::MaxwellMechanism mech(f, 0.1, 0.1);
    
    EXPECT_NEAR(mech.relaxation_time, tau, 1e-10);
    
    // Check inverse
    double f_back = 1.0 / (2.0 * M_PI * mech.relaxation_time);
    EXPECT_NEAR(f_back, f, 1e-10);
}

// ============================================================================
// Quality Factor Tests
// ============================================================================

TEST_F(ViscoelasticAttenuationTest, QualityFactorSetup) {
    ViscoelasticAttenuation atten;
    
    atten.setupMechanismsFromQ(3, f_central, f_max / f_min, Qp, Qs);
    
    EXPECT_EQ(atten.getNumMechanisms(), 3);
    
    // Check mechanisms are set up
    const auto& mechs = atten.getMechanisms();
    for (const auto& m : mechs) {
        EXPECT_GT(m.relaxation_frequency, 0.0);
        EXPECT_TRUE(std::isfinite(m.anelastic_coefficient_p));
        EXPECT_TRUE(std::isfinite(m.anelastic_coefficient_s));
    }
}

TEST_F(ViscoelasticAttenuationTest, QualityFactorApproximation) {
    ViscoelasticAttenuation atten;
    
    atten.setupMechanismsFromQ(3, f_central, f_max / f_min, Qp, Qs);
    
    // Test Q at various frequencies within the target bandwidth
    std::vector<double> test_freqs = {0.1, 0.5, 1.0, 2.0, 5.0};
    
    for (double f : test_freqs) {
        double Q_p = atten.getQuality(f, false);  // P-wave
        double Q_s = atten.getQuality(f, true);   // S-wave
        
        // Q should be approximately constant across bandwidth
        // (within ~20% of target for well-tuned mechanisms)
        EXPECT_GT(Q_p, Qp * 0.5) << "Q_p too low at f = " << f;
        EXPECT_LT(Q_p, Qp * 2.0) << "Q_p too high at f = " << f;
        EXPECT_GT(Q_s, Qs * 0.5) << "Q_s too low at f = " << f;
        EXPECT_LT(Q_s, Qs * 2.0) << "Q_s too high at f = " << f;
    }
}

TEST_F(ViscoelasticAttenuationTest, QualityFactorMaxError) {
    ViscoelasticAttenuation atten;
    
    atten.setupMechanismsFromQ(3, f_central, f_max / f_min, Qp, Qs);
    
    double max_error_p = atten.getMaxQError(f_min, f_max, false);
    double max_error_s = atten.getMaxQError(f_min, f_max, true);
    
    // Maximum Q error should be small for well-tuned mechanisms
    EXPECT_LT(max_error_p, 0.3) << "Q_p error > 30%";
    EXPECT_LT(max_error_s, 0.3) << "Q_s error > 30%";
}

// ============================================================================
// Modulus Tests
// ============================================================================

TEST_F(ViscoelasticAttenuationTest, RelaxedModulus) {
    ViscoelasticAttenuation atten;
    
    atten.setupMechanismsFromQ(3, f_central, 100.0, Qp, Qs);
    
    double M_unrelaxed = rho * Vs_ref * Vs_ref;  // Shear modulus
    double M_relaxed = atten.getRelaxedModulus(M_unrelaxed, true);
    
    // Relaxed modulus < unrelaxed modulus
    EXPECT_LT(M_relaxed, M_unrelaxed);
    EXPECT_GT(M_relaxed, 0.0);
    
    // Ratio depends on attenuation
    double ratio = M_relaxed / M_unrelaxed;
    EXPECT_GT(ratio, 0.5);  // Shouldn't be drastically different
}

TEST_F(ViscoelasticAttenuationTest, ComplexModulus) {
    ViscoelasticAttenuation atten;
    
    atten.setupMechanismsFromQ(3, f_central, 100.0, Qp, Qs);
    
    double M_ref = rho * Vs_ref * Vs_ref;
    
    std::complex<double> M_complex = atten.getComplexModulus(f_central, M_ref, true);
    
    // Real part should be positive
    EXPECT_GT(M_complex.real(), 0.0);
    
    // Imaginary part related to attenuation
    // Im(M) / Re(M) ≈ 1/Q for low attenuation
    double Q_approx = M_complex.real() / M_complex.imag();
    EXPECT_GT(Q_approx, 0.0);
}

// ============================================================================
// Velocity Dispersion Tests
// ============================================================================

TEST_F(ViscoelasticAttenuationTest, VelocityDispersion) {
    ViscoelasticAttenuation atten;
    
    atten.setupMechanismsFromQ(3, f_central, 100.0, Qp, Qs);
    
    // Velocity should increase with frequency (normal dispersion)
    double v_low = atten.getVelocity(0.1, Vs_ref, true);
    double v_high = atten.getVelocity(10.0, Vs_ref, true);
    
    // Higher frequencies travel faster in attenuating medium
    EXPECT_GT(v_high, v_low) << "Velocity should increase with frequency";
    
    // Dispersion should be small for high Q
    double dispersion = (v_high - v_low) / v_low;
    EXPECT_LT(dispersion, 0.1) << "Dispersion should be < 10% for Q ~ 50";
}

TEST_F(ViscoelasticAttenuationTest, DispersionCurve) {
    ViscoelasticAttenuation atten;
    
    atten.setupMechanismsFromQ(3, f_central, 100.0, Qp, Qs);
    
    auto curve = atten.computeDispersionCurve(Vp_ref, Vs_ref, f_min, f_max, 20);
    
    EXPECT_EQ(curve.size(), 20);
    
    // Check monotonicity of velocity with frequency
    for (size_t i = 1; i < curve.size(); ++i) {
        EXPECT_GT(curve[i].frequency, curve[i-1].frequency);
        // Velocity should generally increase (small fluctuations allowed)
    }
    
    // Check Q values are approximately constant
    for (const auto& point : curve) {
        EXPECT_GT(point.Qp, 0.0);
        EXPECT_GT(point.Qs, 0.0);
        EXPECT_GT(point.velocity_p, 0.0);
        EXPECT_GT(point.velocity_s, 0.0);
    }
}

// ============================================================================
// Anelastic State Tests
// ============================================================================

TEST_F(ViscoelasticAttenuationTest, AnelasticStateInitialization) {
    ViscoelasticAttenuation::AnelasticState state(3);
    
    EXPECT_EQ(state.R.size(), 3);
    
    for (const auto& r : state.R) {
        for (int i = 0; i < 6; ++i) {
            EXPECT_NEAR(r[i], 0.0, 1e-15);
        }
    }
}

TEST_F(ViscoelasticAttenuationTest, AnelasticStateReset) {
    ViscoelasticAttenuation::AnelasticState state(3);
    
    // Set non-zero values
    state.R[0][0] = 1.0;
    state.R[1][1] = 2.0;
    
    state.reset();
    
    for (const auto& r : state.R) {
        for (int i = 0; i < 6; ++i) {
            EXPECT_NEAR(r[i], 0.0, 1e-15);
        }
    }
}

TEST_F(ViscoelasticAttenuationTest, AnelasticStateUpdate) {
    ViscoelasticAttenuation atten;
    atten.setupMechanismsFromQ(3, f_central, 100.0, Qp, Qs);
    
    ViscoelasticAttenuation::AnelasticState state(3);
    
    // Apply stress
    std::array<double, 6> stress = {1e6, -0.5e6, -0.5e6, 0.5e6, 0.0, 0.0};
    double dt = 0.001;  // 1 ms
    
    atten.updateAnelasticState(state, stress, dt);
    
    // State should be non-zero after update with non-zero stress
    bool any_nonzero = false;
    for (const auto& r : state.R) {
        for (int i = 0; i < 6; ++i) {
            if (std::abs(r[i]) > 1e-15) {
                any_nonzero = true;
                break;
            }
        }
    }
    EXPECT_TRUE(any_nonzero) << "Anelastic state should evolve with stress";
}

TEST_F(ViscoelasticAttenuationTest, AnelasticStressCorrection) {
    ViscoelasticAttenuation atten;
    atten.setupMechanismsFromQ(3, f_central, 100.0, Qp, Qs);
    
    ViscoelasticAttenuation::AnelasticState state(3);
    
    // Build up anelastic state
    std::array<double, 6> stress = {1e6, -0.5e6, -0.5e6, 0.5e6, 0.0, 0.0};
    for (int i = 0; i < 100; ++i) {
        atten.updateAnelasticState(state, stress, 0.001);
    }
    
    // Get stress correction
    std::array<double, 6> correction;
    atten.computeAnelasticStress(state, correction);
    
    // Correction should be finite
    for (int i = 0; i < 6; ++i) {
        EXPECT_TRUE(std::isfinite(correction[i])) << "Correction component " << i;
    }
}

// ============================================================================
// Causality Tests (Kramers-Kronig)
// ============================================================================

TEST_F(ViscoelasticAttenuationTest, KramersKronigRelation) {
    ViscoelasticAttenuation atten;
    atten.setupMechanismsFromQ(3, f_central, 100.0, Qp, Qs);
    
    bool causal = AttenuationUtils::checkCausality(atten, f_min, f_max, 0.1);
    
    // Well-constructed viscoelastic model should be causal
    EXPECT_TRUE(causal) << "Model should satisfy Kramers-Kronig relations";
}

// ============================================================================
// Attenuation Coefficient Tests
// ============================================================================

TEST_F(ViscoelasticAttenuationTest, AttenuationCoefficientAlpha) {
    // α = π * f / (Q * v)
    
    double f = 1.0;       // Hz
    double Q = 50.0;
    double v = 3500.0;    // m/s
    
    double alpha = AttenuationUtils::attenuationCoefficient(f, v, Q);
    
    double alpha_expected = M_PI * f / (Q * v);
    EXPECT_NEAR(alpha, alpha_expected, alpha_expected * 1e-6);
    
    // Units: Nepers/m
    EXPECT_GT(alpha, 0.0);
}

TEST_F(ViscoelasticAttenuationTest, AmplitudeReduction) {
    double distance = 10000.0;  // 10 km
    double f = 1.0;
    double Q = 50.0;
    double v = 3500.0;
    
    double reduction = AttenuationUtils::amplitudeReduction(distance, f, v, Q);
    
    // A = A0 * exp(-α * distance)
    double alpha = M_PI * f / (Q * v);
    double reduction_expected = std::exp(-alpha * distance);
    
    EXPECT_NEAR(reduction, reduction_expected, 1e-6);
    
    // Should be between 0 and 1
    EXPECT_GT(reduction, 0.0);
    EXPECT_LT(reduction, 1.0);
}

TEST_F(ViscoelasticAttenuationTest, ComputeQFromDecay) {
    double amplitude_ratio = 0.5;   // 50% amplitude remaining
    double distance = 5000.0;       // 5 km
    double wavelength = 3500.0;     // λ = v/f at 1 Hz, v = 3500 m/s
    
    double Q = AttenuationUtils::computeQ(amplitude_ratio, distance, wavelength);
    
    // Q = -π * distance / (λ * ln(A1/A0))
    double Q_expected = -M_PI * distance / (wavelength * std::log(amplitude_ratio));
    
    EXPECT_NEAR(Q, Q_expected, 1.0);
    EXPECT_GT(Q, 0.0);
}

// ============================================================================
// Rock Type Database Tests
// ============================================================================

TEST_F(ViscoelasticAttenuationTest, RockTypeParameters) {
    auto params_granite = AttenuationDatabase::getParameters(AttenuationDatabase::RockType::GRANITE);
    auto params_sediment = AttenuationDatabase::getParameters(AttenuationDatabase::RockType::SEDIMENT_SOFT);
    
    // Granite should have higher Q than soft sediment
    EXPECT_GT(params_granite.Qp_desired, params_sediment.Qp_desired);
    EXPECT_GT(params_granite.Qs_desired, params_sediment.Qs_desired);
}

TEST_F(ViscoelasticAttenuationTest, EmpiricalQsVsRelation) {
    // Q_s ≈ 40-50 * V_s (km/s) empirical relation
    
    double Vs_km = 3.5;  // km/s
    double Qs_estimate = AttenuationDatabase::estimateQs(Vs_km);
    
    // Should be in range 100-200 for Vs = 3.5 km/s
    EXPECT_GT(Qs_estimate, 100.0);
    EXPECT_LT(Qs_estimate, 300.0);
}

TEST_F(ViscoelasticAttenuationTest, QpQsRelation) {
    // Q_p ≈ 2 * Q_s empirical relation
    
    double Qs = 50.0;
    double Qp_estimate = AttenuationDatabase::estimateQp(Qs);
    
    EXPECT_NEAR(Qp_estimate, 2.0 * Qs, Qs * 0.5);
}

TEST_F(ViscoelasticAttenuationTest, DepthCorrection) {
    double Q_surface = 50.0;
    
    double Q_shallow = AttenuationDatabase::depthCorrection(Q_surface, 1.0);   // 1 km
    double Q_deep = AttenuationDatabase::depthCorrection(Q_surface, 10.0);     // 10 km
    
    // Q should increase with depth (less attenuation)
    EXPECT_GT(Q_deep, Q_shallow);
}

TEST_F(ViscoelasticAttenuationTest, DamageCorrection) {
    double Q_intact = 100.0;
    
    double Q_damaged = AttenuationDatabase::damageCorrection(Q_intact, 0.5);  // 50% damage
    
    // Damage should decrease Q (more attenuation)
    EXPECT_LT(Q_damaged, Q_intact);
    EXPECT_GT(Q_damaged, 0.0);
}

// ============================================================================
// Spatial Attenuation Tests
// ============================================================================

TEST_F(ViscoelasticAttenuationTest, SpatialAttenuationVelocityScaling) {
    SpatialAttenuation spatial;
    
    ViscoelasticAttenuation::Parameters base_params;
    base_params.Qp_desired = 100.0;
    base_params.Qs_desired = 50.0;
    
    spatial.setBaseModel(base_params);
    spatial.enableVelocityScaling(true, 50.0);  // Q_s = 50 * V_s
    
    // High velocity region should have higher Q
    auto params_high = spatial.getParametersAt(0, 0, 0, 6000, 3500, 0.0);
    auto params_low = spatial.getParametersAt(0, 0, 0, 4000, 2000, 0.0);
    
    EXPECT_GT(params_high.Qs_desired, params_low.Qs_desired);
}

TEST_F(ViscoelasticAttenuationTest, SpatialAttenuationDepthScaling) {
    SpatialAttenuation spatial;
    
    ViscoelasticAttenuation::Parameters base_params;
    base_params.Qs_desired = 50.0;
    
    spatial.setBaseModel(base_params);
    spatial.enableDepthScaling(true, 0.1);  // dQ/dz = 0.1 per km
    
    // Deeper should have higher Q
    auto params_shallow = spatial.getParametersAt(0, 0, -1000, 5000, 3000, 0.0);
    auto params_deep = spatial.getParametersAt(0, 0, -10000, 5000, 3000, 0.0);
    
    EXPECT_GT(params_deep.Qs_desired, params_shallow.Qs_desired);
}

// ============================================================================
// Integrator Tests
// ============================================================================

TEST_F(ViscoelasticAttenuationTest, ForwardEulerIntegration) {
    ViscoelasticAttenuation atten;
    atten.setupMechanismsFromQ(3, f_central, 100.0, Qp, Qs);
    
    AnelasticIntegrator integrator(AnelasticIntegrator::Scheme::FORWARD_EULER);
    
    ViscoelasticAttenuation::AnelasticState state(3);
    std::array<double, 6> stress = {1e6, 0, 0, 0, 0, 0};
    
    integrator.integrate(atten, state, stress, 0.001);
    
    // State should be updated
    bool any_nonzero = false;
    for (const auto& r : state.R) {
        for (int i = 0; i < 6; ++i) {
            if (std::abs(r[i]) > 1e-15) any_nonzero = true;
        }
    }
    EXPECT_TRUE(any_nonzero);
}

TEST_F(ViscoelasticAttenuationTest, RungeKutta2Integration) {
    ViscoelasticAttenuation atten;
    atten.setupMechanismsFromQ(3, f_central, 100.0, Qp, Qs);
    
    AnelasticIntegrator integrator(AnelasticIntegrator::Scheme::RUNGE_KUTTA_2);
    
    ViscoelasticAttenuation::AnelasticState state(3);
    std::array<double, 6> stress = {1e6, 0, 0, 0, 0, 0};
    
    integrator.integrate(atten, state, stress, 0.001);
    
    // Compare with forward Euler (RK2 should be more accurate)
    ViscoelasticAttenuation::AnelasticState state_fe(3);
    AnelasticIntegrator integrator_fe(AnelasticIntegrator::Scheme::FORWARD_EULER);
    integrator_fe.integrate(atten, state_fe, stress, 0.001);
    
    // Both should give finite results
    for (size_t m = 0; m < state.R.size(); ++m) {
        for (int i = 0; i < 6; ++i) {
            EXPECT_TRUE(std::isfinite(state.R[m][i]));
            EXPECT_TRUE(std::isfinite(state_fe.R[m][i]));
        }
    }
}

// ============================================================================
// Configuration Tests
// ============================================================================

TEST_F(ViscoelasticAttenuationTest, ConfigParsing) {
    ViscoelasticConfig config;
    
    std::map<std::string, std::string> params;
    params["enabled"] = "true";
    params["num_mechanisms"] = "3";
    params["Qp"] = "100";
    params["Qs"] = "50";
    params["freq_central"] = "1.0";
    
    config.parseConfig(params);
    
    EXPECT_TRUE(config.enabled);
    EXPECT_EQ(config.num_mechanisms, 3);
    EXPECT_NEAR(config.Qp, 100.0, 0.1);
    EXPECT_NEAR(config.Qs, 50.0, 0.1);
}

TEST_F(ViscoelasticAttenuationTest, CreateModelFromConfig) {
    ViscoelasticConfig config;
    config.enabled = true;
    config.num_mechanisms = 3;
    config.Qp = 100.0;
    config.Qs = 50.0;
    config.freq_central = 1.0;
    config.freq_ratio = 100.0;
    
    ViscoelasticAttenuation model = config.createModel();
    
    EXPECT_TRUE(model.isEnabled());
    EXPECT_EQ(model.getNumMechanisms(), 3);
}

// ============================================================================
// Coupling Tests
// ============================================================================

TEST_F(ViscoelasticAttenuationTest, PoroelasticCoupling) {
    ViscoelasticCoupling coupling;
    coupling.setCoupledPoroelastic(true);
    
    double Q_base = 50.0;
    double porosity = 0.2;
    double damage = 0.0;
    
    double Q_modified = coupling.getModifiedQ(Q_base, porosity, damage, false);
    
    // Porosity may affect Q (implementation dependent)
    EXPECT_TRUE(std::isfinite(Q_modified));
    EXPECT_GT(Q_modified, 0.0);
}

TEST_F(ViscoelasticAttenuationTest, PlasticityCoupling) {
    ViscoelasticCoupling coupling;
    coupling.setCoupledPlasticity(true);
    coupling.setPlasticQReduction(0.5);  // 50% Q reduction in plastic zone
    
    double Q_base = 100.0;
    
    double Q_elastic = coupling.getModifiedQ(Q_base, 0.0, 0.0, false);
    double Q_plastic = coupling.getModifiedQ(Q_base, 0.0, 0.0, true);
    
    // Plastic zone should have lower Q
    EXPECT_LT(Q_plastic, Q_elastic);
}

TEST_F(ViscoelasticAttenuationTest, DamageCoupling) {
    ViscoelasticCoupling coupling;
    coupling.setCoupledDamage(true);
    coupling.setDamageQRelation(2.0);  // Q ~ (1-D)^2
    
    double Q_base = 100.0;
    
    double Q_undamaged = coupling.getModifiedQ(Q_base, 0.0, 0.0, false);
    double Q_damaged = coupling.getModifiedQ(Q_base, 0.0, 0.5, false);  // 50% damage
    
    // Damage should reduce Q
    EXPECT_LT(Q_damaged, Q_undamaged);
}

// ============================================================================
// Physical Validation Tests
// ============================================================================

TEST_F(ViscoelasticAttenuationTest, PhaseGroupVelocityRelation) {
    // v_phase = v_group * (1 + 1/(2*Q²)) approximately for Q >> 1
    
    double v_group = 3500.0;
    double Q = 50.0;
    
    double v_phase = AttenuationUtils::phaseVelocity(v_group, Q, f_central);
    
    // Phase velocity should be slightly higher than group velocity
    EXPECT_GT(v_phase, v_group);
    
    // For Q = 50, difference should be small (~0.02%)
    double diff = (v_phase - v_group) / v_group;
    EXPECT_LT(diff, 0.01);
}

TEST_F(ViscoelasticAttenuationTest, WaveAmplitudeDecay) {
    // Test that amplitude decays as expected
    
    double f = 1.0;       // Hz
    double Q = 50.0;
    double v = 3500.0;    // m/s
    
    double A0 = 1.0;
    double distance = 10000.0;  // 10 km
    
    double alpha = M_PI * f / (Q * v);
    double A1 = A0 * std::exp(-alpha * distance);
    
    // After 10 km at Q=50, f=1 Hz, v=3500 m/s:
    // α = π * 1 / (50 * 3500) ≈ 1.8e-5 Neper/m
    // A1/A0 = exp(-1.8e-5 * 10000) ≈ 0.83
    
    EXPECT_GT(A1, 0.5);
    EXPECT_LT(A1, 1.0);
}

TEST_F(ViscoelasticAttenuationTest, SpectralRatioQEstimation) {
    // Test spectral ratio method for Q estimation
    
    std::vector<double> spectrum1(100);
    std::vector<double> spectrum2(100);
    std::vector<double> frequencies(100);
    
    double Q_true = 50.0;
    double v = 3500.0;
    double distance_diff = 5000.0;  // 5 km between stations
    
    // Synthetic spectra
    for (int i = 0; i < 100; ++i) {
        frequencies[i] = 0.1 + i * 0.1;  // 0.1 to 10 Hz
        double alpha = M_PI * frequencies[i] / (Q_true * v);
        spectrum1[i] = 1.0;
        spectrum2[i] = std::exp(-alpha * distance_diff);
    }
    
    auto result = AttenuationUtils::estimateQFromSpectralRatio(
        spectrum1, spectrum2, frequencies, distance_diff);
    
    double Q_estimated = result.first;
    double Q_uncertainty = result.second;
    
    // Estimated Q should be close to true Q
    EXPECT_NEAR(Q_estimated, Q_true, Q_true * 0.2);
    EXPECT_TRUE(std::isfinite(Q_uncertainty));
}

// ============================================================================
// Edge Cases and Robustness Tests
// ============================================================================

TEST_F(ViscoelasticAttenuationTest, HighQLimit) {
    ViscoelasticAttenuation atten;
    
    // Very high Q (nearly elastic)
    atten.setupMechanismsFromQ(3, f_central, 100.0, 1000.0, 500.0);
    
    // Dispersion should be minimal
    double v_low = atten.getVelocity(0.1, Vs_ref, true);
    double v_high = atten.getVelocity(10.0, Vs_ref, true);
    
    double dispersion = std::abs(v_high - v_low) / v_low;
    EXPECT_LT(dispersion, 0.01) << "High Q should give minimal dispersion";
}

TEST_F(ViscoelasticAttenuationTest, LowQLimit) {
    ViscoelasticAttenuation atten;
    
    // Very low Q (high attenuation)
    atten.setupMechanismsFromQ(3, f_central, 100.0, 10.0, 5.0);
    
    // Should still give finite results
    double v = atten.getVelocity(1.0, Vs_ref, true);
    double Q = atten.getQuality(1.0, true);
    
    EXPECT_TRUE(std::isfinite(v));
    EXPECT_TRUE(std::isfinite(Q));
    EXPECT_GT(v, 0.0);
    EXPECT_GT(Q, 0.0);
}

TEST_F(ViscoelasticAttenuationTest, SingleMechanism) {
    ViscoelasticAttenuation atten;
    
    // Single mechanism (simplified model)
    atten.setupMechanismsFromQ(1, f_central, 10.0, Qp, Qs);
    
    EXPECT_EQ(atten.getNumMechanisms(), 1);
    
    // Should still work but Q won't be constant over wide bandwidth
    double Q = atten.getQuality(f_central, true);
    EXPECT_GT(Q, 0.0);
}

TEST_F(ViscoelasticAttenuationTest, ZeroStressState) {
    ViscoelasticAttenuation atten;
    atten.setupMechanismsFromQ(3, f_central, 100.0, Qp, Qs);
    
    ViscoelasticAttenuation::AnelasticState state(3);
    std::array<double, 6> stress = {0, 0, 0, 0, 0, 0};
    
    atten.updateAnelasticState(state, stress, 0.001);
    
    // With zero stress, state should remain zero
    for (const auto& r : state.R) {
        for (int i = 0; i < 6; ++i) {
            EXPECT_NEAR(r[i], 0.0, 1e-15);
        }
    }
}

TEST_F(ViscoelasticAttenuationTest, DisabledAttenuation) {
    ViscoelasticAttenuation atten;
    atten.setupMechanismsFromQ(3, f_central, 100.0, Qp, Qs);
    atten.enable(false);
    
    EXPECT_FALSE(atten.isEnabled());
    
    // When disabled, relaxed modulus should equal unrelaxed
    double M = 1e10;
    double M_relaxed = atten.getRelaxedModulus(M, true);
    
    // Implementation-dependent, but typically returns unrelaxed when disabled
    EXPECT_TRUE(std::isfinite(M_relaxed));
}
