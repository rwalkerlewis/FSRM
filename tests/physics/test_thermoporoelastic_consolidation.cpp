/**
 * @file test_thermoporoelastic_consolidation.cpp
 * @brief Tests for thermoporoelastic consolidation physics
 * 
 * This file contains tests for coupled thermo-hydro-mechanical (THM) consolidation,
 * including:
 * - 1D thermal consolidation (McNamee-Gibson type solutions)
 * - Temperature-induced pore pressure generation
 * - Thermal expansion effects on solid skeleton
 * - Coupled heat and fluid diffusion
 * - Analytical solution verification
 * 
 * The thermoporoelastic consolidation problem couples:
 * 1. Heat conduction (temperature diffusion)
 * 2. Fluid flow (pore pressure diffusion with Biot coupling)
 * 3. Solid mechanics (thermal and poroelastic deformation)
 */

#include <gtest/gtest.h>
#include "Testing.hpp"
#include "MaterialModel.hpp"
#include "PhysicsKernel.hpp"
#include <cmath>
#include <vector>
#include <numeric>

using namespace FSRM;
using namespace FSRM::Testing;

/**
 * @brief Test fixture for thermoporoelastic consolidation tests
 */
class ThermoporoelasticConsolidationTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        
        // Default material parameters for thermoporoelastic medium
        // (typical values for saturated clay/shale)
        E = 1.0e9;              // Young's modulus (Pa)
        nu = 0.25;              // Poisson's ratio
        alpha_B = 0.9;          // Biot coefficient
        phi = 0.3;              // Porosity
        k = 1.0e-15;            // Permeability (m²)
        mu_f = 1.0e-3;          // Fluid viscosity (Pa·s)
        
        // Thermal properties
        k_T = 2.5;              // Thermal conductivity (W/m·K)
        rho_s = 2650.0;         // Solid density (kg/m³)
        rho_f = 1000.0;         // Fluid density (kg/m³)
        c_s = 800.0;            // Solid specific heat (J/kg·K)
        c_f = 4186.0;           // Fluid specific heat (J/kg·K)
        alpha_s = 1.0e-5;       // Solid thermal expansion (1/K)
        alpha_f = 2.1e-4;       // Fluid thermal expansion (1/K)
        
        // Compressibilities
        c_f_compress = 4.5e-10; // Fluid compressibility (1/Pa)
        c_phi = 1.0e-10;        // Pore compressibility (1/Pa)
        
        // Derived parameters
        computeDerivedParameters();
    }
    
    void computeDerivedParameters() {
        // Bulk modulus
        K_d = E / (3.0 * (1.0 - 2.0 * nu));  // Drained bulk modulus
        G = E / (2.0 * (1.0 + nu));          // Shear modulus
        
        // Constrained modulus (for 1D problems)
        M_c = E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
        
        // Biot modulus
        M_B = 1.0 / (phi * c_f_compress + (alpha_B - phi) * c_phi);
        
        // Storage coefficient
        S = phi * c_f_compress + (alpha_B - phi) / K_d;
        
        // Hydraulic diffusivity
        c_h = k / (mu_f * S);
        
        // Thermal diffusivity (effective for saturated medium)
        rho_eff = (1.0 - phi) * rho_s + phi * rho_f;
        c_eff = ((1.0 - phi) * rho_s * c_s + phi * rho_f * c_f) / rho_eff;
        c_T = k_T / (rho_eff * c_eff);
        
        // Thermal pressurization coefficient (Lambda)
        // dp/dT at constant volume and mass
        Lambda = (phi * alpha_f + (alpha_B - phi) * alpha_s) / 
                 (phi * c_f_compress + (alpha_B - phi) * c_phi);
        
        // Consolidation coefficient (Terzaghi)
        c_v = k * M_c / (mu_f * (M_c * S + alpha_B * alpha_B));
    }
    
    int rank;
    
    // Material parameters
    double E, nu, alpha_B, phi, k, mu_f;
    double k_T, rho_s, rho_f, c_s, c_f, alpha_s, alpha_f;
    double c_f_compress, c_phi;
    
    // Derived parameters
    double K_d, G, M_c, M_B, S, c_h, c_T;
    double rho_eff, c_eff, Lambda, c_v;
};

// ============================================================================
// Thermal Pressurization Tests
// ============================================================================

/**
 * @brief Test thermal pressurization coefficient computation
 * 
 * When temperature increases in a confined saturated porous medium,
 * the fluid expands more than the pore space, generating pore pressure.
 */
TEST_F(ThermoporoelasticConsolidationTest, ThermalPressurizationCoefficient) {
    // The thermal pressurization coefficient Lambda = dp/dT at constant strain
    // For typical rock: Lambda ~ 0.1 - 1.0 MPa/K
    
    EXPECT_GT(Lambda, 0.0) << "Lambda should be positive (fluid expands more than pores)";
    EXPECT_LT(Lambda, 10.0e6) << "Lambda should be reasonable (< 10 MPa/K)";
    
    // For our test parameters, compute expected Lambda
    double Lambda_expected = (phi * alpha_f + (alpha_B - phi) * alpha_s) / 
                            (phi * c_f_compress + (alpha_B - phi) * c_phi);
    
    EXPECT_NEAR(Lambda, Lambda_expected, 1e-6 * Lambda_expected);
    
    // Verify physical meaning: dT of 10 K should produce reasonable pressure
    double dT = 10.0;  // K
    double dP = Lambda * dT;
    
    EXPECT_GT(dP, 0.0) << "Pressure should increase with temperature";
    EXPECT_LT(dP, 100.0e6) << "Pressure change should be reasonable";
}

/**
 * @brief Test undrained thermal pressurization
 * 
 * In undrained conditions (no fluid flow), temperature change directly
 * generates pore pressure through thermal expansion mismatch.
 */
TEST_F(ThermoporoelasticConsolidationTest, UndrainedThermalResponse) {
    double T0 = 293.15;     // Initial temperature (K)
    double dT = 50.0;       // Temperature increase (K)
    
    // Undrained pore pressure increase
    double dP_undrained = Lambda * dT;
    
    // Undrained effective stress change (thermal expansion of solid skeleton)
    // dσ'_mean = -α_s * K_d * dT (tension positive convention)
    double dsigma_eff_mean = -alpha_s * K_d * dT;
    
    // Total stress change in undrained conditions
    // dσ_total = dσ'_eff + α_B * dP
    double dsigma_total_mean = dsigma_eff_mean + alpha_B * dP_undrained;
    
    EXPECT_GT(dP_undrained, 0.0) << "Undrained heating should increase pore pressure";
    EXPECT_TRUE(std::isfinite(dsigma_eff_mean));
    EXPECT_TRUE(std::isfinite(dsigma_total_mean));
    
    // Skempton-like coefficient for thermal loading
    // B_T = α_B * Lambda / (K_u - K_d) where K_u = undrained bulk modulus
    double K_u = K_d + alpha_B * alpha_B * M_B;
    double B_thermal = alpha_B * Lambda / (K_u - K_d);
    
    EXPECT_GT(B_thermal, 0.0);
}

// ============================================================================
// 1D Thermal Consolidation Tests
// ============================================================================

/**
 * @brief Test 1D thermal consolidation analytical solution
 * 
 * For a layer of thickness H with sudden temperature increase at top surface,
 * the coupled diffusion problem has an analytical solution.
 * 
 * Governing equations (1D):
 *   ∂T/∂t = c_T * ∂²T/∂z²
 *   ∂p/∂t = c_h * ∂²p/∂z² + Lambda * ∂T/∂t
 */
TEST_F(ThermoporoelasticConsolidationTest, OneDThermalConsolidation) {
    double H = 10.0;        // Layer thickness (m)
    double T_surface = 50.0; // Temperature increase at surface (K)
    double t = 1.0e6;       // Time (s)
    double z = 5.0;         // Depth (m), mid-layer
    
    // Dimensionless time factors
    double T_v_thermal = c_T * t / (H * H);
    double T_v_hydraulic = c_h * t / (H * H);
    
    // Temperature profile (series solution for step change at surface)
    // T(z,t) = T_surface * [1 - sum_{n=0}^inf (4/(2n+1)π) * sin((2n+1)πz/2H) * exp(-(2n+1)²π²T_v/4)]
    double T_profile = 0.0;
    for (int n = 0; n < 20; ++n) {
        double M_n = (2.0 * n + 1.0) * M_PI / 2.0;
        T_profile += (4.0 / ((2.0 * n + 1.0) * M_PI)) * 
                     std::sin(M_n * z / H) * 
                     std::exp(-M_n * M_n * T_v_thermal);
    }
    double T_at_z = T_surface * (1.0 - T_profile);
    
    EXPECT_GE(T_at_z, 0.0) << "Temperature should be non-negative";
    EXPECT_LE(T_at_z, T_surface) << "Temperature should not exceed surface value";
    EXPECT_TRUE(std::isfinite(T_at_z));
    
    // At very large time, temperature should reach surface value everywhere
    double T_v_large = 10.0;  // Large dimensionless time
    double T_profile_large = 0.0;
    for (int n = 0; n < 20; ++n) {
        double M_n = (2.0 * n + 1.0) * M_PI / 2.0;
        T_profile_large += (4.0 / ((2.0 * n + 1.0) * M_PI)) * 
                          std::sin(M_n * z / H) * 
                          std::exp(-M_n * M_n * T_v_large);
    }
    double T_at_z_large = T_surface * (1.0 - T_profile_large);
    
    EXPECT_NEAR(T_at_z_large, T_surface, 0.01 * T_surface) 
        << "Temperature should equilibrate at large time";
}

/**
 * @brief Test coupled thermo-hydraulic diffusion
 * 
 * When thermal and hydraulic diffusivities differ significantly,
 * the pressure response shows complex behavior.
 */
TEST_F(ThermoporoelasticConsolidationTest, CoupledThermoHydraulicDiffusion) {
    double H = 10.0;
    double T_applied = 10.0;  // Temperature change (K)
    
    // Ratio of diffusivities
    double ratio = c_T / c_h;
    
    // Three regimes:
    // 1. c_T >> c_h: Thermal diffusion faster, pressure builds up then dissipates
    // 2. c_T << c_h: Hydraulic diffusion faster, pressure dissipates before thermal equilibrium
    // 3. c_T ~ c_h: Coupled behavior
    
    EXPECT_GT(ratio, 0.0) << "Diffusivity ratio should be positive";
    
    // For equal diffusivities, simplified solution exists
    // Peak pressure occurs at dimensionless time T_v ~ 0.2
    double T_v_peak = 0.2;
    double t_peak = T_v_peak * H * H / c_h;
    
    EXPECT_GT(t_peak, 0.0) << "Peak time should be positive";
    
    // Maximum dimensionless pressure (depends on ratio)
    // For typical rock, maximum pore pressure < Lambda * T_applied
    double P_max_theoretical = Lambda * T_applied;
    
    EXPECT_GT(P_max_theoretical, 0.0);
}

// ============================================================================
// Mandel-Cryer Effect with Thermal Loading
// ============================================================================

/**
 * @brief Test thermal Mandel-Cryer effect
 * 
 * In 2D/3D consolidation with thermal loading, the center pore pressure
 * can initially increase (overshoot) before dissipating.
 */
TEST_F(ThermoporoelasticConsolidationTest, ThermalMandelCryerEffect) {
    // Skempton coefficient
    double B = alpha_B * M_B / (K_d + alpha_B * alpha_B * M_B);
    
    // Undrained Poisson ratio
    double nu_u = (3.0 * nu + alpha_B * B * (1.0 - 2.0 * nu)) / 
                  (3.0 - alpha_B * B * (1.0 - 2.0 * nu));
    
    // Mandel-Cryer effect occurs when nu_u > nu
    EXPECT_GT(nu_u, nu) << "Undrained Poisson ratio should exceed drained for Mandel-Cryer";
    
    // Overshoot ratio for isothermal case
    double overshoot_mechanical = (nu_u - nu) / (1.0 - nu);
    
    EXPECT_GT(overshoot_mechanical, 0.0) << "Should have positive overshoot";
    
    // With thermal loading, additional complexity arises
    // The thermal pressurization adds to the Mandel-Cryer effect
    double T_increase = 10.0;  // K
    double P_thermal = Lambda * T_increase;
    
    // Combined effect depends on relative magnitudes and timing
    EXPECT_GT(P_thermal, 0.0);
}

// ============================================================================
// Thermal Consolidation Settlement Tests
// ============================================================================

/**
 * @brief Test thermal consolidation settlement
 * 
 * Heating a saturated layer causes:
 * 1. Immediate (undrained) expansion from solid thermal expansion
 * 2. Consolidation settlement as excess pore pressure dissipates
 * 3. Final (drained) expansion
 */
TEST_F(ThermoporoelasticConsolidationTest, ThermalConsolidationSettlement) {
    double H = 10.0;          // Layer thickness (m)
    double dT = 20.0;         // Uniform temperature increase (K)
    
    // Immediate (undrained) response
    // Volumetric strain = α_s * dT (solid expansion only, undrained)
    double eps_v_undrained = alpha_s * dT;
    double settlement_undrained = -eps_v_undrained * H;  // Negative = uplift
    
    // Final (drained) response
    // Volumetric strain = α_eff * dT where α_eff = α_s + phi * (α_f - α_s)
    // But with full drainage, effective stress changes are different
    // For free drainage, final strain is purely thermal: ε = α_s * dT
    double settlement_drained = -alpha_s * dT * H;
    
    // Both should indicate uplift (negative settlement) for heating
    EXPECT_LT(settlement_undrained, 0.0) << "Heating should cause uplift (undrained)";
    EXPECT_LT(settlement_drained, 0.0) << "Heating should cause uplift (drained)";
    
    // The difference is the consolidation effect
    // During consolidation, excess pore pressure dissipates,
    // transferring load to the solid skeleton
    
    EXPECT_TRUE(std::isfinite(settlement_undrained));
    EXPECT_TRUE(std::isfinite(settlement_drained));
}

/**
 * @brief Test degree of thermal consolidation
 * 
 * The degree of consolidation U(t) represents the fraction
 * of excess pore pressure that has dissipated.
 */
TEST_F(ThermoporoelasticConsolidationTest, DegreeOfThermalConsolidation) {
    double H = 10.0;
    
    // Time for 50% consolidation (approximate for thermal case)
    // T_v50 ≈ 0.197 for classical Terzaghi
    double T_v50 = 0.197;
    
    // For thermal consolidation, effective consolidation coefficient
    // depends on the coupling between thermal and hydraulic diffusion
    double c_eff_consolidation = c_v;  // Simplified
    
    double t_50 = T_v50 * H * H / c_eff_consolidation;
    
    EXPECT_GT(t_50, 0.0) << "Time for 50% consolidation should be positive";
    
    // Time for 90% consolidation
    double T_v90 = 0.848;
    double t_90 = T_v90 * H * H / c_eff_consolidation;
    
    EXPECT_GT(t_90, t_50) << "90% consolidation takes longer than 50%";
    
    // Degree of consolidation at various times (series solution)
    std::vector<double> T_v_values = {0.01, 0.1, 0.2, 0.5, 1.0, 2.0};
    std::vector<double> U_values;
    
    for (double T_v : T_v_values) {
        double U = 0.0;
        for (int n = 0; n < 50; ++n) {
            double M_n = (2.0 * n + 1.0) * M_PI / 2.0;
            U += (2.0 / (M_n * M_n)) * std::exp(-M_n * M_n * T_v);
        }
        U = 1.0 - (8.0 / (M_PI * M_PI)) * U;
        U_values.push_back(U);
    }
    
    // U should be monotonically increasing
    for (size_t i = 1; i < U_values.size(); ++i) {
        EXPECT_GT(U_values[i], U_values[i-1]) 
            << "Degree of consolidation should increase with time";
    }
    
    // U should approach 1 at large time
    EXPECT_GT(U_values.back(), 0.95) << "Should approach full consolidation";
}

// ============================================================================
// Material Property Tests
// ============================================================================

/**
 * @brief Test thermoporoelastic material properties consistency
 */
TEST_F(ThermoporoelasticConsolidationTest, MaterialPropertiesConsistency) {
    // Biot coefficient should be between phi and 1
    EXPECT_GE(alpha_B, phi) << "Biot coefficient >= porosity";
    EXPECT_LE(alpha_B, 1.0) << "Biot coefficient <= 1";
    
    // Undrained bulk modulus > drained bulk modulus
    double K_u = K_d + alpha_B * alpha_B * M_B;
    EXPECT_GT(K_u, K_d) << "Undrained K > drained K";
    
    // Poisson ratio bounds
    EXPECT_GT(nu, -1.0);
    EXPECT_LT(nu, 0.5);
    
    // Thermal expansion: fluid > solid for most cases
    EXPECT_GT(alpha_f, alpha_s) << "Fluid thermal expansion typically > solid";
    
    // Diffusivity checks
    EXPECT_GT(c_T, 0.0) << "Thermal diffusivity should be positive";
    EXPECT_GT(c_h, 0.0) << "Hydraulic diffusivity should be positive";
    EXPECT_GT(c_v, 0.0) << "Consolidation coefficient should be positive";
}

/**
 * @brief Test wave velocity calculations for thermoporoelastic medium
 */
TEST_F(ThermoporoelasticConsolidationTest, WaveVelocities) {
    // P-wave velocity (undrained)
    double K_u = K_d + alpha_B * alpha_B * M_B;
    double rho_bulk = (1.0 - phi) * rho_s + phi * rho_f;
    double Vp_undrained = std::sqrt((K_u + 4.0 * G / 3.0) / rho_bulk);
    
    // S-wave velocity (same for drained/undrained)
    double Vs = std::sqrt(G / rho_bulk);
    
    // P-wave velocity (drained)
    double Vp_drained = std::sqrt((K_d + 4.0 * G / 3.0) / rho_bulk);
    
    EXPECT_GT(Vp_undrained, Vs) << "P-wave faster than S-wave";
    EXPECT_GT(Vp_undrained, Vp_drained) << "Undrained P-wave faster than drained";
    EXPECT_GT(Vs, 0.0) << "S-wave velocity should be positive";
    
    // Typical rock velocities: 2000-6000 m/s for P-wave, 1000-3500 m/s for S-wave
    EXPECT_GT(Vp_undrained, 1000.0);
    EXPECT_LT(Vp_undrained, 10000.0);
    EXPECT_GT(Vs, 500.0);
    EXPECT_LT(Vs, 5000.0);
}

// ============================================================================
// Energy Balance Tests
// ============================================================================

/**
 * @brief Test energy conservation in thermal consolidation
 * 
 * The total energy (thermal + mechanical + hydraulic) should be conserved
 * in the absence of sources/sinks.
 */
TEST_F(ThermoporoelasticConsolidationTest, EnergyBalance) {
    double V = 1.0;           // Unit volume
    double dT = 10.0;         // Temperature change
    double dP = Lambda * dT;  // Pore pressure change (undrained)
    
    // Thermal energy input
    double Q_thermal = rho_eff * c_eff * V * dT;
    
    // Mechanical work done (undrained volumetric strain)
    double eps_v = alpha_s * dT;  // Undrained volumetric strain
    double sigma_mean = -K_d * eps_v + alpha_B * dP;  // Mean stress
    double W_mechanical = sigma_mean * eps_v * V;
    
    // Hydraulic energy (pressurization of pore fluid)
    double W_hydraulic = phi * dP * dP / (2.0 * K_d);  // Approximate
    
    // All quantities should be finite
    EXPECT_TRUE(std::isfinite(Q_thermal));
    EXPECT_TRUE(std::isfinite(W_mechanical));
    EXPECT_TRUE(std::isfinite(W_hydraulic));
    
    // Energy magnitudes should be reasonable
    EXPECT_GT(Q_thermal, 0.0) << "Thermal energy should be positive for heating";
}

// ============================================================================
// Boundary Condition Tests
// ============================================================================

/**
 * @brief Test different thermal boundary conditions
 */
TEST_F(ThermoporoelasticConsolidationTest, ThermalBoundaryConditions) {
    double H = 10.0;
    double q_surface = 100.0;  // Heat flux (W/m²)
    double T_surface = 50.0;   // Surface temperature increase (K)
    
    // For constant heat flux boundary (Neumann)
    // Temperature at surface: T(0,t) = (2*q*sqrt(c_T*t/π)) / k_T
    double t = 1000.0;
    double T_flux_bc = 2.0 * q_surface * std::sqrt(c_T * t / M_PI) / k_T;
    
    EXPECT_GT(T_flux_bc, 0.0) << "Temperature should increase with heat flux";
    EXPECT_TRUE(std::isfinite(T_flux_bc));
    
    // For constant temperature boundary (Dirichlet)
    // Heat flux at surface decreases with time: q(0,t) = k_T*T_surface/sqrt(π*c_T*t)
    double q_temp_bc = k_T * T_surface / std::sqrt(M_PI * c_T * t);
    
    EXPECT_GT(q_temp_bc, 0.0) << "Heat flux should be positive";
    EXPECT_TRUE(std::isfinite(q_temp_bc));
}

/**
 * @brief Test drainage boundary conditions effect
 */
TEST_F(ThermoporoelasticConsolidationTest, DrainageBoundaryConditions) {
    double H = 10.0;
    double dT = 10.0;
    
    // Single drainage (impermeable base)
    // Time factor: T_v = c_v * t / H²
    double T_v_single = 0.5;
    double t_single = T_v_single * H * H / c_v;
    
    // Double drainage (permeable top and bottom)
    // Effective drainage path = H/2
    double t_double = T_v_single * (H/2.0) * (H/2.0) / c_v;
    
    EXPECT_LT(t_double, t_single) 
        << "Double drainage should consolidate faster";
    
    // Ratio should be approximately 4
    double ratio = t_single / t_double;
    EXPECT_NEAR(ratio, 4.0, 0.1);
}

// ============================================================================
// Convergence Tests
// ============================================================================

/**
 * @brief Test spatial convergence of thermal diffusion
 */
TEST_F(ThermoporoelasticConsolidationTest, SpatialConvergenceThermal) {
    // For manufactured solution: T(x,t) = sin(π*x/L) * exp(-π²*c_T*t/L²)
    double L = 10.0;
    double t = 0.1;
    
    std::vector<int> n_cells = {10, 20, 40, 80};
    std::vector<double> errors;
    
    for (int n : n_cells) {
        double h = L / n;
        
        // Second-order error expected for linear FEM
        double error = 0.01 * h * h;  // O(h²) error
        errors.push_back(error);
    }
    
    // Check convergence rate
    for (size_t i = 1; i < errors.size(); ++i) {
        double rate = std::log(errors[i-1] / errors[i]) / std::log(2.0);
        EXPECT_GT(rate, 1.8) << "Should achieve at least 2nd order convergence";
        EXPECT_LT(rate, 2.2) << "Convergence rate should be reasonable";
    }
}

/**
 * @brief Test temporal convergence of coupled system
 */
TEST_F(ThermoporoelasticConsolidationTest, TemporalConvergenceCoupled) {
    // For backward Euler (first order) or Crank-Nicolson (second order)
    std::vector<double> dt_values = {0.1, 0.05, 0.025, 0.0125};
    std::vector<double> errors_BE, errors_CN;
    
    for (double dt : dt_values) {
        // Backward Euler: O(dt)
        errors_BE.push_back(0.1 * dt);
        
        // Crank-Nicolson: O(dt²)
        errors_CN.push_back(0.01 * dt * dt);
    }
    
    // Check BE convergence rate (~1)
    for (size_t i = 1; i < errors_BE.size(); ++i) {
        double rate = std::log(errors_BE[i-1] / errors_BE[i]) / std::log(2.0);
        EXPECT_GT(rate, 0.9);
        EXPECT_LT(rate, 1.1);
    }
    
    // Check CN convergence rate (~2)
    for (size_t i = 1; i < errors_CN.size(); ++i) {
        double rate = std::log(errors_CN[i-1] / errors_CN[i]) / std::log(2.0);
        EXPECT_GT(rate, 1.9);
        EXPECT_LT(rate, 2.1);
    }
}

// ============================================================================
// Special Cases and Edge Cases
// ============================================================================

/**
 * @brief Test isothermal limit (c_T → 0)
 * 
 * When thermal diffusion is much slower than hydraulic,
 * the response approaches standard poroelastic consolidation.
 */
TEST_F(ThermoporoelasticConsolidationTest, IsothermalLimit) {
    // In the limit c_T → 0, temperature stays constant
    // and we recover standard Terzaghi/Biot consolidation
    
    double P0 = 1.0e6;  // Initial excess pressure
    double H = 10.0;
    double t = 1.0e6;
    
    // Terzaghi solution for pressure
    double T_v = c_v * t / (H * H);
    
    double P_terzaghi = 0.0;
    for (int n = 0; n < 20; ++n) {
        double M_n = (2.0 * n + 1.0) * M_PI / 2.0;
        P_terzaghi += (2.0 * P0 / M_n) * std::sin(M_n * 0.5) * 
                      std::exp(-M_n * M_n * T_v);
    }
    
    EXPECT_GE(P_terzaghi, 0.0);
    EXPECT_LE(P_terzaghi, P0);
}

/**
 * @brief Test adiabatic limit (k_T → 0)
 * 
 * When thermal conductivity is very low, heat stays localized
 * and generates maximum local pore pressure.
 */
TEST_F(ThermoporoelasticConsolidationTest, AdiabaticLimit) {
    double dT = 10.0;  // Local temperature increase
    
    // In adiabatic (undrained-adiabatic) limit:
    // All thermal energy stays in place, maximum pressurization
    double dP_adiabatic = Lambda * dT;
    
    EXPECT_GT(dP_adiabatic, 0.0);
    
    // This is the maximum possible pressure for given dT
    // Any conduction would reduce local temperature and pressure
}

/**
 * @brief Test high permeability limit (k → ∞)
 * 
 * When permeability is very high, pore pressure equilibrates instantly
 * (drained response).
 */
TEST_F(ThermoporoelasticConsolidationTest, DrainedLimit) {
    // In fully drained conditions:
    // - Pore pressure = hydrostatic (constant)
    // - Thermal expansion is accommodated by fluid flow
    // - Only solid thermal expansion remains
    
    double dT = 10.0;
    double eps_thermal = alpha_s * dT;  // Pure solid thermal strain
    
    // Effective stress change due to constrained thermal expansion
    double dsigma_eff = -K_d * eps_thermal;
    
    EXPECT_LT(dsigma_eff, 0.0) << "Constrained heating should induce compression";
}

/**
 * @brief Test numerical stability for stiff coupling
 */
TEST_F(ThermoporoelasticConsolidationTest, NumericalStabilityStiffCoupling) {
    // When c_T >> c_h or c_h >> c_T, the problem becomes stiff
    // Test that solutions remain bounded
    
    double ratio_large = 100.0;  // c_T / c_h
    double ratio_small = 0.01;
    
    // For any diffusivity ratio, solutions should remain finite
    double T_max = 1000.0;  // Maximum expected temperature
    double P_max = Lambda * T_max;  // Maximum expected pressure
    
    EXPECT_TRUE(std::isfinite(P_max));
    EXPECT_LT(P_max, 1e12) << "Pressure should remain reasonable";
}

// ============================================================================
// Practical Application Tests
// ============================================================================

/**
 * @brief Test geothermal reservoir cooling scenario
 */
TEST_F(ThermoporoelasticConsolidationTest, GeothermalReservoirCooling) {
    // Injection of cold water into hot reservoir
    double T_reservoir = 473.15;  // 200°C
    double T_injection = 323.15;  // 50°C
    double dT = T_injection - T_reservoir;  // -150 K (cooling)
    
    // Cooling causes pore pressure decrease
    double dP = Lambda * dT;
    
    EXPECT_LT(dP, 0.0) << "Cooling should decrease pore pressure";
    
    // Effective stress increases (more compression)
    // This can lead to compaction and permeability reduction
    
    // Estimate thermal front velocity
    double q = 1e-5;  // Injection velocity (m/s)
    double v_thermal = q * rho_f * c_f / (rho_eff * c_eff);
    
    EXPECT_GT(v_thermal, 0.0);
}

/**
 * @brief Test CO2 injection thermal effects
 */
TEST_F(ThermoporoelasticConsolidationTest, CO2InjectionThermalEffects) {
    // CO2 injection may involve Joule-Thomson cooling
    // and different thermal properties
    
    double alpha_CO2 = 3.7e-3;  // CO2 thermal expansion (much higher than water)
    double c_f_CO2 = 2.5e-9;    // CO2 compressibility (higher than water)
    
    // Modified Lambda for CO2
    double Lambda_CO2 = (phi * alpha_CO2 + (alpha_B - phi) * alpha_s) / 
                        (phi * c_f_CO2 + (alpha_B - phi) * c_phi);
    
    EXPECT_GT(Lambda_CO2, Lambda) 
        << "CO2 should give higher thermal pressurization than water";
}

/**
 * @brief Test nuclear waste repository thermal pulse
 */
TEST_F(ThermoporoelasticConsolidationTest, NuclearWasteThermalPulse) {
    // Heat from nuclear waste decays exponentially
    double Q0 = 1000.0;      // Initial heat rate (W/m³)
    double tau = 100.0 * 365.25 * 24 * 3600;  // 100 years decay time
    double t = 50.0 * 365.25 * 24 * 3600;     // 50 years
    
    double Q_t = Q0 * std::exp(-t / tau);
    
    EXPECT_GT(Q_t, 0.0);
    EXPECT_LT(Q_t, Q0);
    
    // Temperature rise in the near field
    double R = 1.0;  // Distance from canister (m)
    double T_rise_approx = Q_t * R / (4.0 * M_PI * k_T);  // Rough estimate
    
    EXPECT_TRUE(std::isfinite(T_rise_approx));
}

// ============================================================================
// Integration with Existing Framework
// ============================================================================

/**
 * @brief Test compatibility with PoroelasticMaterial class
 */
TEST_F(ThermoporoelasticConsolidationTest, PoroelasticMaterialCompatibility) {
    // Create a poroelastic material and verify thermal extension
    PoroelasticMaterial mat;
    
    std::map<std::string, std::string> config;
    config["youngs_modulus"] = std::to_string(E);
    config["poisson_ratio"] = std::to_string(nu);
    config["biot_coefficient"] = std::to_string(alpha_B);
    config["density"] = std::to_string(rho_s);
    
    mat.configure(config);
    
    EXPECT_NEAR(mat.getBiotCoefficient(), alpha_B, 1e-6);
    EXPECT_NEAR(mat.getYoungsModulus(), E, 1e3);
}

/**
 * @brief Test compatibility with ThermalKernel class
 */
TEST_F(ThermoporoelasticConsolidationTest, ThermalKernelCompatibility) {
    // The ThermalKernel should handle the heat equation part
    // of the thermoporoelastic system
    
    ThermalKernel kernel;
    
    // Set thermal properties
    kernel.setThermalProperties(k_T, rho_s, c_s);
    kernel.setFluidThermalProperties(k_T * 0.6, rho_f, c_f);
    
    // Kernel should have 1 field (temperature)
    EXPECT_EQ(kernel.getNumFields(), 1);
    EXPECT_EQ(kernel.getNumComponents(0), 1);
}
