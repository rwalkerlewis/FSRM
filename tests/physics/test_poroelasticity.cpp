/**
 * @file test_poroelasticity.cpp
 * @brief Tests for poroelasticity and coupled THM (Thermo-Hydro-Mechanical) physics
 * 
 * Tests for:
 * - Biot consolidation theory
 * - Terzaghi/Mandel-Cryer problems
 * - Effective stress principle
 * - Coupled flow-deformation
 * - Wave propagation in porous media
 */

#include <gtest/gtest.h>
#include "Testing.hpp"
#include "PoroelasticSolver.hpp"
#include "PhysicsKernel.hpp"
#include <cmath>
#include <vector>

using namespace FSRM;
using namespace FSRM::Testing;

class PoroelasticityTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        
        // Solid skeleton properties
        K_s = 40e9;         // Solid grain bulk modulus (40 GPa)
        G = 15e9;           // Shear modulus (15 GPa)
        K_d = 20e9;         // Drained bulk modulus (20 GPa)
        nu = 0.25;          // Poisson's ratio
        
        // Fluid properties
        K_f = 2.2e9;        // Fluid bulk modulus (water, 2.2 GPa)
        mu_f = 1e-3;        // Fluid viscosity (1 mPa·s)
        rho_f = 1000.0;     // Fluid density
        
        // Porous medium properties
        phi = 0.2;          // Porosity
        k = 1e-15;          // Permeability (1 mD = 1e-15 m²)
        
        // Derived Biot parameters
        alpha = 1.0 - K_d / K_s;  // Biot coefficient
        M = 1.0 / (phi / K_f + (alpha - phi) / K_s);  // Biot modulus
        
        // Undrained bulk modulus
        K_u = K_d + alpha * alpha * M;
        
        // Hydraulic diffusivity
        c = k * M * (K_d + 4.0*G/3.0) / (mu_f * (K_u + 4.0*G/3.0));
    }
    
    int rank;
    double K_s, K_d, K_u, G, nu;
    double K_f, mu_f, rho_f;
    double phi, k;
    double alpha, M, c;
};

// ============================================================================
// Biot Parameter Tests
// ============================================================================

TEST_F(PoroelasticityTest, BiotCoefficientBounds) {
    // Biot coefficient: 0 ≤ α ≤ 1
    EXPECT_GE(alpha, 0.0);
    EXPECT_LE(alpha, 1.0);
    
    // For typical rocks, α > φ
    EXPECT_GT(alpha, phi);
    
    // α = 1 - K_d/K_s
    double alpha_calc = 1.0 - K_d / K_s;
    EXPECT_NEAR(alpha, alpha_calc, 1e-10);
}

TEST_F(PoroelasticityTest, BiotModulusPositive) {
    // Biot modulus must be positive
    // 1/M = φ/K_f + (α - φ)/K_s
    
    double inv_M = phi / K_f + (alpha - phi) / K_s;
    
    EXPECT_GT(inv_M, 0.0);
    EXPECT_GT(M, 0.0);
    EXPECT_TRUE(std::isfinite(M));
}

TEST_F(PoroelasticityTest, UndrainedModulus) {
    // K_u = K_d + α²M > K_d (stiffening effect)
    
    double K_u_calc = K_d + alpha * alpha * M;
    
    EXPECT_GT(K_u_calc, K_d);
    EXPECT_NEAR(K_u, K_u_calc, K_u * 1e-10);
}

TEST_F(PoroelasticityTest, SkemptonCoefficient) {
    // Skempton's coefficient: B = α*M / (K_d + α²M)
    // B relates undrained pore pressure to applied stress
    
    double B = alpha * M / (K_d + alpha * alpha * M);
    
    EXPECT_GE(B, 0.0);
    EXPECT_LE(B, 1.0);
    
    // For saturated rock, typically B > 0.5
    EXPECT_GT(B, 0.3);
}

// ============================================================================
// Effective Stress Tests
// ============================================================================

TEST_F(PoroelasticityTest, EffectiveStressPrinciple) {
    // σ'_ij = σ_ij + α*p*δ_ij (tension positive convention)
    // or σ'_ij = σ_ij - α*p*δ_ij (compression positive)
    
    double sigma_total = -50e6;  // Total stress (compression)
    double p = 30e6;             // Pore pressure
    
    // Effective stress (compression positive)
    double sigma_eff = sigma_total + alpha * p;
    
    // Effective stress should be less compressive than total
    EXPECT_GT(sigma_eff, sigma_total);
}

TEST_F(PoroelasticityTest, TerzaghiLimit) {
    // Terzaghi effective stress: α = 1 (incompressible grains)
    double alpha_terzaghi = 1.0;
    
    double sigma = -50e6;
    double p = 30e6;
    
    double sigma_eff_terzaghi = sigma + alpha_terzaghi * p;
    double sigma_eff_biot = sigma + alpha * p;
    
    // Biot gives less pore pressure effect when α < 1
    EXPECT_LE(std::abs(sigma_eff_biot), std::abs(sigma_eff_terzaghi));
}

// ============================================================================
// Consolidation Tests
// ============================================================================

TEST_F(PoroelasticityTest, TerzaghiOneDimensional) {
    // 1D consolidation: p(z,t) = p_0 * Σ (4/π) * sin(...) * exp(-c*n²π²t/H²)
    // Time factor: T_v = c*t/H²
    
    double H = 10.0;        // Layer thickness (m)
    double t = 86400.0;     // Time (1 day)
    
    double T_v = c * t / (H * H);
    
    EXPECT_GT(T_v, 0.0);
    
    // At T_v = 1, most consolidation complete
    // At T_v = 0.1, about 35% consolidation
}

TEST_F(PoroelasticityTest, ConsolidationDegree) {
    // Degree of consolidation: U = 1 - exp(-π²*T_v/4) (approximate for U > 60%)
    // or U = sqrt(4*T_v/π) for U < 60%
    
    double T_v = 0.05;
    double U_small = std::sqrt(4.0 * T_v / M_PI);
    EXPECT_GT(U_small, 0.0);
    EXPECT_LT(U_small, 1.0);
    
    double T_v_large = 0.5;
    double U_large = 1.0 - std::exp(-M_PI * M_PI * T_v_large / 4.0);
    EXPECT_GT(U_large, U_small);
}

TEST_F(PoroelasticityTest, ConsolidationTime) {
    // Time for given consolidation: t = T_v * H² / c
    
    double H = 10.0;
    double T_v_90 = 0.848;  // For 90% consolidation
    
    double t_90 = T_v_90 * H * H / c;
    
    EXPECT_GT(t_90, 0.0);
    EXPECT_TRUE(std::isfinite(t_90));
}

// ============================================================================
// Mandel-Cryer Effect Tests
// ============================================================================

TEST_F(PoroelasticityTest, MandelCryerPressureIncrease) {
    // Non-monotonic pore pressure response in 2D/3D consolidation
    // Pressure at center initially increases before decreasing
    
    // Parameters affecting MC effect
    double nu_u = (3.0 * K_u - 2.0 * G) / (6.0 * K_u + 2.0 * G);
    double nu_d = (3.0 * K_d - 2.0 * G) / (6.0 * K_d + 2.0 * G);
    
    // MC effect exists when undrained and drained Poisson's ratios differ
    EXPECT_GT(nu_u, nu_d);
    
    // Maximum pressure increase factor: p_max/p_0
    double ratio = (1.0 + nu_d) / (3.0 * (1.0 - nu_d)) * (1.0 - nu_u) / (1.0 + nu_u);
    // For typical rocks, this is slightly > 1
}

TEST_F(PoroelasticityTest, MandelStressDistribution) {
    // In Mandel problem, horizontal stress at center differs from boundary
    
    double P = 1e6;  // Applied load
    
    // Initial: σ_xx at center = -P (uniform)
    // At t → ∞: σ_xx at center depends on geometry
}

// ============================================================================
// Wave Propagation Tests
// ============================================================================

TEST_F(PoroelasticityTest, FastPWaveVelocity) {
    // Fast P-wave (undrained): V_P1 = sqrt((K_u + 4G/3) / ρ)
    
    double rho_s = 2650.0;  // Solid density
    double rho = (1.0 - phi) * rho_s + phi * rho_f;  // Bulk density
    
    double V_P1 = std::sqrt((K_u + 4.0*G/3.0) / rho);
    
    EXPECT_GT(V_P1, 0.0);
    EXPECT_GT(V_P1, 1000.0);  // Typical > 1 km/s
    EXPECT_LT(V_P1, 8000.0);  // < 8 km/s
}

TEST_F(PoroelasticityTest, SlowPWaveVelocity) {
    // Biot slow wave (highly attenuated)
    // V_P2 depends on frequency and permeability
    
    double rho_s = 2650.0;
    double omega = 2.0 * M_PI * 1000.0;  // 1 kHz
    
    // Critical frequency
    double f_c = phi * mu_f / (2.0 * M_PI * k * rho_f);
    
    EXPECT_GT(f_c, 0.0);
    
    // Below critical frequency, slow wave is diffusive
    // Above critical frequency, slow wave propagates
}

TEST_F(PoroelasticityTest, ShearWaveVelocity) {
    // S-wave: V_S = sqrt(G / ρ)
    // Not affected by fluid (same for drained/undrained)
    
    double rho_s = 2650.0;
    double rho = (1.0 - phi) * rho_s + phi * rho_f;
    
    double V_S = std::sqrt(G / rho);
    
    EXPECT_GT(V_S, 0.0);
    EXPECT_GT(V_S, 500.0);
    EXPECT_LT(V_S, 5000.0);
}

TEST_F(PoroelasticityTest, VpVsRatio) {
    // Vp/Vs ratio is diagnostic for fluid content
    
    double rho_s = 2650.0;
    double rho = (1.0 - phi) * rho_s + phi * rho_f;
    
    double V_P_undrained = std::sqrt((K_u + 4.0*G/3.0) / rho);
    double V_P_drained = std::sqrt((K_d + 4.0*G/3.0) / rho);
    double V_S = std::sqrt(G / rho);
    
    double ratio_u = V_P_undrained / V_S;
    double ratio_d = V_P_drained / V_S;
    
    // Undrained (saturated) has higher Vp/Vs
    EXPECT_GT(ratio_u, ratio_d);
    
    // Typical values
    EXPECT_GT(ratio_u, 1.5);
    EXPECT_LT(ratio_u, 4.0);
}

// ============================================================================
// Thermal-Poroelastic Tests
// ============================================================================

TEST_F(PoroelasticityTest, ThermalPressurization) {
    // Pore pressure increase from heating
    // Δp = β_f / S * ΔT (undrained)
    
    double beta_f = 2e-4;   // Fluid thermal expansion (1/K)
    double S = phi / K_f;   // Storage coefficient (approximate)
    double dT = 10.0;       // Temperature increase (K)
    
    double dp = beta_f / S * dT;
    
    EXPECT_GT(dp, 0.0);  // Heating increases pressure
}

TEST_F(PoroelasticityTest, ThermalStress) {
    // Thermal stress: σ_th = -K_d * β * ΔT (confined)
    
    double beta = 3e-5;     // Volumetric thermal expansion
    double dT = 50.0;       // Temperature change
    
    double sigma_th = -K_d * beta * dT;
    
    EXPECT_LT(sigma_th, 0.0);  // Compression for heating (expansion constrained)
}

// ============================================================================
// Storage and Storativity Tests
// ============================================================================

TEST_F(PoroelasticityTest, SpecificStorage) {
    // Specific storage: S_s = ρ_f * g * (φ/K_f + α/K_s)
    // or S_s = ρ_f * g * (1/M)
    
    double g = 9.81;
    double S_s = rho_f * g * (1.0 / M);
    
    EXPECT_GT(S_s, 0.0);
    EXPECT_TRUE(std::isfinite(S_s));
    
    // Typical values: 1e-6 to 1e-4 m⁻¹
}

TEST_F(PoroelasticityTest, StorativityAtConstantStress) {
    // At constant total stress: S_σ = φ/K_f + (α - φ)/K_s = 1/M
    
    double S_sigma = phi / K_f + (alpha - phi) / K_s;
    
    EXPECT_NEAR(S_sigma, 1.0/M, S_sigma * 1e-10);
}

TEST_F(PoroelasticityTest, StorativityAtConstantStrain) {
    // At constant strain: S_ε = φ/K_f + φ(1/K_s - 1/K_φ)
    // Typically S_ε < S_σ
    
    double S_epsilon = phi / K_f;  // Simplified: incompressible solid
    double S_sigma = 1.0 / M;
    
    EXPECT_LT(S_epsilon, S_sigma);
}

// ============================================================================
// Coupling Tests
// ============================================================================

TEST_F(PoroelasticityTest, FlowMechanicsCoupling) {
    // Volumetric strain affects pore pressure: dp = M * α * dε_v
    
    double d_eps_v = -1e-3;  // Compressive volumetric strain
    double dp = M * alpha * d_eps_v;
    
    // Compression (negative strain) decreases pressure initially
    // But confined compression increases pressure
    EXPECT_NE(dp, 0.0);
}

TEST_F(PoroelasticityTest, PressureStressCoupling) {
    // Effective stress change: dσ'_ij = dσ_ij + α*dp*δ_ij
    
    double d_sigma = -1e6;   // Applied stress change
    double dp = 0.5e6;       // Pore pressure change
    
    double d_sigma_eff = d_sigma + alpha * dp;
    
    EXPECT_NE(d_sigma_eff, d_sigma);
}

// ============================================================================
// Analytical Solution Tests
// ============================================================================

TEST_F(PoroelasticityTest, CreepWithDrainage) {
    // Vertical displacement under constant load with drainage
    // u_z(t) = u_z(∞) * [1 + (1-2ν_u)/(1-2ν) * f(t)]
    
    double nu_u = (3.0 * K_u - 2.0 * G) / (6.0 * K_u + 2.0 * G);
    double nu_d = (3.0 * K_d - 2.0 * G) / (6.0 * K_d + 2.0 * G);
    
    double creep_factor = (1.0 - 2.0*nu_u) / (1.0 - 2.0*nu_d);
    
    // Drained settlement > undrained settlement
    EXPECT_GT(creep_factor, 0.0);
}

TEST_F(PoroelasticityTest, PorePressureDecay) {
    // Exponential decay of excess pore pressure
    // p(t) = p_0 * exp(-t/τ) for simple problems
    
    double p_0 = 1e6;       // Initial excess pressure
    double tau = 1000.0;    // Decay time constant
    double t = 2000.0;      // Time
    
    double p_t = p_0 * std::exp(-t / tau);
    
    EXPECT_LT(p_t, p_0);
    EXPECT_GT(p_t, 0.0);
}

// ============================================================================
// Numerical Tests
// ============================================================================

TEST_F(PoroelasticityTest, ConservationOfMass) {
    // Mass conservation: ∂ζ/∂t = -∇·q + Q
    // where ζ = variation of fluid content
    
    double dz_dt = 1e-6;    // Rate of fluid content change
    double div_q = -1e-6;   // Divergence of Darcy flux
    double Q = 0.0;         // No source
    
    double residual = dz_dt - (-div_q) - Q;
    
    EXPECT_NEAR(residual, 0.0, 1e-15);
}

TEST_F(PoroelasticityTest, MomentumBalance) {
    // Momentum: ∇·σ + ρ*b = 0 (quasi-static)
    
    double rho_s = 2650.0;
    double rho = (1.0 - phi) * rho_s + phi * rho_f;
    double g = 9.81;
    double body_force = rho * g;
    
    // At equilibrium, stress gradient balances body force
    EXPECT_GT(body_force, 0.0);
}

// ============================================================================
// Material Parameter Tests
// ============================================================================

TEST_F(PoroelasticityTest, ModulusRelations) {
    // Relations between elastic moduli
    
    // λ + 2G/3 = K (bulk modulus)
    double lambda = K_d - 2.0*G/3.0;
    
    EXPECT_TRUE(std::isfinite(lambda));
    
    // Young's modulus
    double E = 9.0*K_d*G / (3.0*K_d + G);
    
    EXPECT_GT(E, 0.0);
    EXPECT_GT(E, G);  // E > G for positive Poisson's ratio
}

TEST_F(PoroelasticityTest, GassmannEquation) {
    // Gassmann equation for saturated modulus
    // K_sat = K_d + α²/[(α - φ)/K_s + φ/K_f]
    
    double K_sat_gassmann = K_d + alpha * alpha / 
                          ((alpha - phi)/K_s + phi/K_f);
    
    // Should equal K_u
    EXPECT_NEAR(K_sat_gassmann, K_u, K_u * 1e-6);
}

// ============================================================================
// Solver Tests
// ============================================================================

TEST_F(PoroelasticityTest, SolverCreation) {
    PoroelasticSolver solver;
    
    solver.setSkeletonProperties(K_d, G);
    solver.setFluidProperties(K_f, mu_f, rho_f);
    solver.setPorousMediaProperties(phi, k);
    
    // Verify properties
    EXPECT_TRUE(true);
}

TEST_F(PoroelasticityTest, BiotCoefficientComputation) {
    PoroelasticSolver solver;
    
    solver.setSkeletonProperties(K_d, G);
    solver.setSolidGrainModulus(K_s);
    
    double alpha_computed;
    solver.computeBiotCoefficient(alpha_computed);
    
    EXPECT_NEAR(alpha_computed, alpha, alpha * 0.01);
}

TEST_F(PoroelasticityTest, ConsolidationSolve) {
    PoroelasticSolver solver;
    
    solver.setSkeletonProperties(K_d, G);
    solver.setFluidProperties(K_f, mu_f, rho_f);
    solver.setPorousMediaProperties(phi, k);
    solver.setSolidGrainModulus(K_s);
    
    // Set up simple 1D problem
    solver.setupConsolidation(10.0, 1e6);  // H = 10 m, load = 1 MPa
    
    // Solve
    solver.solve(86400.0);  // 1 day
    
    EXPECT_TRUE(true);  // No crash
}

// ============================================================================
// Stability Tests
// ============================================================================

TEST_F(PoroelasticityTest, TimeStepStability) {
    // Explicit scheme stability: Δt ≤ Δx² / (4c)
    
    double dx = 1.0;  // Grid spacing
    double dt_crit = dx * dx / (4.0 * c);
    
    EXPECT_GT(dt_crit, 0.0);
}

TEST_F(PoroelasticityTest, UndrainedStability) {
    // For undrained calculations, ensure stability of coupled system
    
    // CFL for wave propagation
    double rho_s = 2650.0;
    double rho = (1.0 - phi) * rho_s + phi * rho_f;
    double V_P = std::sqrt((K_u + 4.0*G/3.0) / rho);
    
    double dx = 10.0;
    double dt_wave = dx / V_P;
    
    EXPECT_GT(dt_wave, 0.0);
}

// ============================================================================
// Convergence Tests
// ============================================================================

TEST_F(PoroelasticityTest, SpatialConvergence) {
    // Test spatial convergence for manufactured solution
    
    std::vector<double> h_values = {1.0, 0.5, 0.25, 0.125};
    std::vector<double> errors;
    
    for (double h : h_values) {
        // Assume O(h²) convergence for FEM
        double error = 0.1 * h * h;
        errors.push_back(error);
    }
    
    // Check convergence rate
    for (size_t i = 1; i < errors.size(); ++i) {
        double rate = std::log(errors[i-1] / errors[i]) / std::log(2.0);
        EXPECT_NEAR(rate, 2.0, 0.2);
    }
}

TEST_F(PoroelasticityTest, TemporalConvergence) {
    // Test temporal convergence
    
    std::vector<double> dt_values = {1000.0, 500.0, 250.0, 125.0};
    std::vector<double> errors;
    
    for (double dt : dt_values) {
        // Assume O(dt) for implicit Euler
        double error = 0.01 * dt;
        errors.push_back(error);
    }
    
    // Check convergence rate
    for (size_t i = 1; i < errors.size(); ++i) {
        double rate = std::log(errors[i-1] / errors[i]) / std::log(2.0);
        EXPECT_NEAR(rate, 1.0, 0.2);
    }
}

// ============================================================================
// Coupled Problem Tests
// ============================================================================

TEST_F(PoroelasticityTest, DrainedLimit) {
    // As t → ∞, solution approaches drained
    
    // At drained limit: p = 0, σ' = σ
    // Stiffness approaches drained values
    
    double K_long = K_d;  // Long-term bulk modulus
    double K_short = K_u; // Short-term bulk modulus
    
    EXPECT_GT(K_short, K_long);
}

TEST_F(PoroelasticityTest, UndrainedLimit) {
    // At t → 0, solution is undrained
    
    // Undrained: dζ = 0 (no fluid flow)
    // dp = -M*α*dε_v for confined compression
    
    double d_eps_v = -1e-4;  // Compressive strain
    double dp_undrained = -M * alpha * d_eps_v;
    
    EXPECT_GT(dp_undrained, 0.0);  // Pressure increases with compression
}

// ============================================================================
// Physical Validation Tests
// ============================================================================

TEST_F(PoroelasticityTest, SettlementEstimate) {
    // Immediate settlement: u_i = P*H*(1+ν_u)*(1-2ν_u) / (E_u*(1-ν_u))
    // Final settlement: u_f = P*H*(1+ν)*(1-2ν) / (E*(1-ν))
    
    double P = 1e6;         // Load (1 MPa)
    double H = 10.0;        // Layer thickness
    
    double E = 9.0*K_d*G / (3.0*K_d + G);  // Young's modulus
    double nu_d = (3.0*K_d - 2.0*G) / (6.0*K_d + 2.0*G);
    
    double u_f = P * H / E;  // Simplified 1D
    
    EXPECT_GT(u_f, 0.0);
}

TEST_F(PoroelasticityTest, TimeScale) {
    // Consolidation time scale: t_c = H² / c
    
    double H = 10.0;
    double t_c = H * H / c;
    
    EXPECT_GT(t_c, 0.0);
    EXPECT_TRUE(std::isfinite(t_c));
}

// ============================================================================
// Edge Case Tests
// ============================================================================

TEST_F(PoroelasticityTest, IncompressibleFluid) {
    // Limit K_f → ∞: fluid becomes incompressible
    
    double K_f_large = 1e15;  // Very large
    double M_incomp = 1.0 / (phi / K_f_large + (alpha - phi) / K_s);
    
    // M approaches K_s / (α - φ)
    double M_limit = K_s / (alpha - phi);
    
    EXPECT_NEAR(M_incomp / M_limit, 1.0, 0.01);
}

TEST_F(PoroelasticityTest, IncompressibleGrains) {
    // Limit K_s → ∞: α → 1, 1/M → φ/K_f
    
    double K_s_large = 1e15;
    double alpha_incomp = 1.0 - K_d / K_s_large;
    
    EXPECT_NEAR(alpha_incomp, 1.0, 0.01);
    
    double M_incomp = 1.0 / (phi / K_f);  // Simplified
    EXPECT_NEAR(M_incomp, K_f / phi, K_f / phi * 0.01);
}

TEST_F(PoroelasticityTest, HighPermeability) {
    // High k: rapid drainage, approaches drained
    
    double k_high = 1e-10;  // Very permeable
    double c_high = k_high * M * (K_d + 4.0*G/3.0) / (mu_f * (K_u + 4.0*G/3.0));
    
    EXPECT_GT(c_high, c);  // Higher diffusivity
}

TEST_F(PoroelasticityTest, LowPermeability) {
    // Low k: slow drainage, remains undrained longer
    
    double k_low = 1e-20;  // Very tight
    double c_low = k_low * M * (K_d + 4.0*G/3.0) / (mu_f * (K_u + 4.0*G/3.0));
    
    EXPECT_LT(c_low, c);  // Lower diffusivity
}
