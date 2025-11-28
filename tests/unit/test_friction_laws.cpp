/**
 * @file test_friction_laws.cpp
 * @brief Unit tests for friction laws
 * 
 * Tests cover:
 * - Coulomb friction
 * - Rate-and-state friction (aging/slip laws)
 * - Flash heating
 * - Thermal pressurization
 * - Strong velocity weakening
 */

#include "FaultModel.hpp"
#include "Testing.hpp"
#include <cmath>
#include <iostream>
#include <vector>

namespace FSRM {
namespace Testing {

// =============================================================================
// Coulomb Friction Tests
// =============================================================================

/**
 * @test Test static Coulomb friction
 */
bool test_coulomb_static() {
    std::cout << "Testing static Coulomb friction..." << std::endl;
    
    double tol = 1e-10;
    
    CoulombFriction cf;
    cf.setStaticFriction(0.6);
    cf.setDynamicFriction(0.4);
    cf.setCriticalSlipDistance(0.4);  // Dc for transition
    
    double sigma_n = 100e6;  // 100 MPa normal stress
    
    // At zero slip, should be static friction
    double mu = cf.getFriction(sigma_n, 0.0, 0.0);
    
    if (std::abs(mu - 0.6) > tol) {
        std::cerr << "  FAIL: Static friction = " << mu << ", expected 0.6" << std::endl;
        return false;
    }
    
    // Shear strength
    double tau = mu * sigma_n;
    double expected_tau = 0.6 * 100e6;
    
    if (std::abs(tau - expected_tau) > tol) {
        std::cerr << "  FAIL: Shear strength incorrect" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Static Coulomb friction correct" << std::endl;
    return true;
}

/**
 * @test Test slip-weakening Coulomb friction
 */
bool test_coulomb_slip_weakening() {
    std::cout << "Testing slip-weakening friction..." << std::endl;
    
    double tol = 1e-6;
    
    CoulombFriction cf;
    cf.setStaticFriction(0.7);
    cf.setDynamicFriction(0.3);
    cf.setCriticalSlipDistance(0.5);  // 0.5m Dc
    
    double sigma_n = 50e6;
    
    // At half Dc, friction should be interpolated
    double mu_half = cf.getFriction(sigma_n, 0.25, 0.0);
    double expected_half = 0.7 - 0.5 * (0.7 - 0.3);
    
    if (std::abs(mu_half - expected_half) > tol) {
        std::cerr << "  FAIL: Half-way friction = " << mu_half 
                  << ", expected " << expected_half << std::endl;
        return false;
    }
    
    // Beyond Dc, should be dynamic
    double mu_full = cf.getFriction(sigma_n, 1.0, 0.0);
    
    if (std::abs(mu_full - 0.3) > tol) {
        std::cerr << "  FAIL: Full slip friction = " << mu_full 
                  << ", expected 0.3" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Slip-weakening friction correct" << std::endl;
    return true;
}

// =============================================================================
// Rate-and-State Tests
// =============================================================================

/**
 * @test Test rate-and-state steady-state friction
 */
bool test_rate_state_steady_state() {
    std::cout << "Testing rate-and-state steady-state..." << std::endl;
    
    double tol = 1e-6;
    
    RateStateFriction rsf;
    rsf.setReferenceFriction(0.6);
    rsf.setDirectEffect_a(0.01);
    rsf.setEvolutionEffect_b(0.015);
    rsf.setCriticalSlipDistance(0.01);  // 10 mm
    rsf.setReferenceVelocity(1e-6);
    
    double sigma_n = 50e6;
    double V = 1e-3;  // 1 mm/s
    
    // At steady-state, θ = Dc/V
    double Dc = 0.01;
    double theta_ss = Dc / V;
    
    // Steady-state friction
    // μ = μ0 + a*ln(V/V0) + b*ln(θV0/Dc)
    // At steady-state: θ = Dc/V, so b*ln(θV0/Dc) = b*ln(V0/V)
    // μ_ss = μ0 + (a-b)*ln(V/V0)
    
    double a = 0.01;
    double b = 0.015;
    double mu0 = 0.6;
    double V0 = 1e-6;
    
    double mu_ss_expected = mu0 + (a - b) * std::log(V / V0);
    double mu_computed = rsf.getSteadyStateFriction(sigma_n, V);
    
    if (std::abs(mu_computed - mu_ss_expected) > tol) {
        std::cerr << "  FAIL: Steady-state μ = " << mu_computed 
                  << ", expected " << mu_ss_expected << std::endl;
        return false;
    }
    
    // For a < b (velocity weakening), μ should decrease with V
    if (a < b) {
        double mu_fast = rsf.getSteadyStateFriction(sigma_n, 1.0);
        double mu_slow = rsf.getSteadyStateFriction(sigma_n, 1e-6);
        
        if (mu_fast >= mu_slow) {
            std::cerr << "  FAIL: Velocity weakening not observed" << std::endl;
            return false;
        }
    }
    
    std::cout << "  PASS: Rate-and-state steady-state correct" << std::endl;
    return true;
}

/**
 * @test Test rate-and-state aging law
 */
bool test_rate_state_aging_law() {
    std::cout << "Testing rate-and-state aging law..." << std::endl;
    
    double tol = 1e-10;
    
    RateStateFriction rsf(RateStateLaw::AGING);
    rsf.setCriticalSlipDistance(0.01);
    
    // Aging law: dθ/dt = 1 - θV/Dc
    // At V=0: dθ/dt = 1 (linear healing)
    // At steady-state: dθ/dt = 0 => θ = Dc/V
    
    double Dc = 0.01;
    double V = 0;
    double theta = 0.1;
    
    double dtheta_dt = rsf.getStateEvolutionRate(theta, V);
    
    // At V=0, should be 1
    if (std::abs(dtheta_dt - 1.0) > tol) {
        std::cerr << "  FAIL: Aging rate at V=0 = " << dtheta_dt 
                  << ", expected 1.0" << std::endl;
        return false;
    }
    
    // At steady-state
    V = 1e-3;
    theta = Dc / V;  // 10 s
    dtheta_dt = rsf.getStateEvolutionRate(theta, V);
    
    if (std::abs(dtheta_dt) > tol) {
        std::cerr << "  FAIL: Not at steady-state, dθ/dt = " << dtheta_dt << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Rate-and-state aging law correct" << std::endl;
    return true;
}

/**
 * @test Test rate-and-state slip law
 */
bool test_rate_state_slip_law() {
    std::cout << "Testing rate-and-state slip law..." << std::endl;
    
    double tol = 1e-10;
    
    RateStateFriction rsf(RateStateLaw::SLIP);
    rsf.setCriticalSlipDistance(0.01);
    
    // Slip law: dθ/dt = -θV/Dc * ln(θV/Dc)
    // At steady-state: ln(θV/Dc) = 0 => θV/Dc = 1 => θ = Dc/V
    
    double Dc = 0.01;
    double V = 1e-3;
    double theta_ss = Dc / V;
    
    double dtheta_dt = rsf.getStateEvolutionRate(theta_ss, V);
    
    if (std::abs(dtheta_dt) > tol) {
        std::cerr << "  FAIL: Not at steady-state, dθ/dt = " << dtheta_dt << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Rate-and-state slip law correct" << std::endl;
    return true;
}

// =============================================================================
// Flash Heating Tests
// =============================================================================

/**
 * @test Test flash heating weakening velocity
 */
bool test_flash_heating_weakening_velocity() {
    std::cout << "Testing flash heating weakening velocity..." << std::endl;
    
    FlashHeatingFriction fh;
    fh.setLowVelocityFriction(0.7);
    fh.setContactDiameter(10e-6);  // 10 μm
    fh.setWeakeningTemperature(900.0);  // 900°C
    fh.setThermalDiffusivity(1e-6);
    fh.setHeatCapacity(1000.0);
    fh.setDensity(2700.0);
    fh.setAmbientTemperature(100.0);
    
    // Compute weakening velocity
    // Vw ≈ π * α * ρ * c * (Tw - T0) / (μ * σn * D)
    double Vw = fh.computeWeakeningVelocity(100e6);  // 100 MPa
    
    // Should be on order of 0.1 m/s for typical rock
    if (Vw < 0.01 || Vw > 10.0) {
        std::cerr << "  FAIL: Unreasonable weakening velocity Vw = " << Vw << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Flash heating weakening velocity reasonable" << std::endl;
    return true;
}

/**
 * @test Test flash heating friction reduction
 */
bool test_flash_heating_friction() {
    std::cout << "Testing flash heating friction reduction..." << std::endl;
    
    FlashHeatingFriction fh;
    fh.setLowVelocityFriction(0.7);
    fh.setContactDiameter(10e-6);
    fh.setWeakeningTemperature(900.0);
    fh.setThermalDiffusivity(1e-6);
    fh.setHeatCapacity(1000.0);
    fh.setDensity(2700.0);
    fh.setAmbientTemperature(100.0);
    
    double sigma_n = 50e6;
    
    // Low velocity - should be near static friction
    double mu_slow = fh.getFriction(sigma_n, 1e-6, 0.0);
    if (std::abs(mu_slow - 0.7) > 0.05) {
        std::cerr << "  FAIL: Low-velocity friction not near static" << std::endl;
        return false;
    }
    
    // High velocity - should be reduced
    double mu_fast = fh.getFriction(sigma_n, 5.0, 0.0);
    if (mu_fast >= 0.7) {
        std::cerr << "  FAIL: High-velocity friction not reduced" << std::endl;
        return false;
    }
    
    // Should always be positive
    if (mu_fast <= 0) {
        std::cerr << "  FAIL: Negative friction at high velocity" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Flash heating friction reduction correct" << std::endl;
    return true;
}

// =============================================================================
// Thermal Pressurization Tests
// =============================================================================

/**
 * @test Test thermal pressurization pressure evolution
 */
bool test_thermal_pressurization_pressure() {
    std::cout << "Testing thermal pressurization pressure evolution..." << std::endl;
    
    ThermalPressurizationFriction tp;
    tp.setBaseFriction(0.6);
    tp.setHydraulicDiffusivity(1e-4);
    tp.setThermalDiffusivity(1e-6);
    tp.setSlipZoneWidth(0.001);  // 1 mm
    tp.setPressurization(1e6);   // 1 MPa/K
    tp.setHeatCapacity(1000.0);
    tp.setDensity(2700.0);
    
    double sigma_n = 100e6;
    double initial_pore_pressure = 20e6;
    
    // Evolve during slip
    double delta_p = 0.0;
    double delta_T = 0.0;
    double dt = 0.001;
    
    for (int i = 0; i < 100; ++i) {
        double tau = (sigma_n - initial_pore_pressure - delta_p) * 0.6;
        double V = 1.0;  // 1 m/s slip rate
        
        tp.evolveTPState(tau, V, dt, delta_p, delta_T);
    }
    
    // Pore pressure should have increased
    if (delta_p <= 0) {
        std::cerr << "  FAIL: Pore pressure didn't increase during slip" << std::endl;
        return false;
    }
    
    // Temperature should have increased
    if (delta_T <= 0) {
        std::cerr << "  FAIL: Temperature didn't increase during slip" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Thermal pressurization evolution correct" << std::endl;
    return true;
}

/**
 * @test Test thermal pressurization friction reduction
 */
bool test_thermal_pressurization_friction() {
    std::cout << "Testing thermal pressurization friction reduction..." << std::endl;
    
    ThermalPressurizationFriction tp;
    tp.setBaseFriction(0.6);
    tp.setHydraulicDiffusivity(1e-4);
    tp.setThermalDiffusivity(1e-6);
    tp.setSlipZoneWidth(0.001);
    tp.setPressurization(1e6);
    
    double sigma_n = 100e6;
    double pore_pressure_0 = 20e6;
    
    // Initial effective normal stress and strength
    double sigma_eff_0 = sigma_n - pore_pressure_0;  // 80 MPa
    double tau_0 = tp.getFriction(sigma_n, 0.0, 0.0) * sigma_eff_0;
    
    // After pore pressure increase
    double delta_p = 30e6;  // 30 MPa increase
    double sigma_eff_1 = sigma_n - pore_pressure_0 - delta_p;  // 50 MPa
    
    // Strength should be reduced
    double tau_1 = tp.getFriction(sigma_n, 0.0, 0.0) * sigma_eff_1;
    
    if (tau_1 >= tau_0) {
        std::cerr << "  FAIL: Strength not reduced by TP" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Thermal pressurization friction reduction correct" << std::endl;
    return true;
}

/**
 * @test Test TP characteristic slip distance
 */
bool test_tp_characteristic_slip() {
    std::cout << "Testing TP characteristic slip distance..." << std::endl;
    
    ThermalPressurizationFriction tp;
    tp.setHydraulicDiffusivity(1e-4);
    tp.setThermalDiffusivity(1e-6);
    tp.setSlipZoneWidth(0.001);
    
    double L_star = tp.getCharacteristicSlip();
    
    // Should be positive and reasonable (mm to cm scale)
    if (L_star <= 0 || L_star > 1.0) {
        std::cerr << "  FAIL: Unreasonable L* = " << L_star << std::endl;
        return false;
    }
    
    std::cout << "  PASS: TP characteristic slip distance reasonable" << std::endl;
    return true;
}

// =============================================================================
// Strong Velocity Weakening Tests
// =============================================================================

/**
 * @test Test SVW friction at different velocities
 */
bool test_svw_friction() {
    std::cout << "Testing strong velocity weakening friction..." << std::endl;
    
    StrongVelocityWeakeningFriction svw;
    svw.setStaticFriction(0.7);
    svw.setDynamicFriction(0.2);
    svw.setCriticalVelocity(0.1);  // 0.1 m/s
    svw.setWeakeningRate(5.0);
    
    double sigma_n = 50e6;
    
    // At very low velocity, should be near static
    double mu_slow = svw.getFriction(sigma_n, 1e-6, 0.0);
    if (std::abs(mu_slow - 0.7) > 0.05) {
        std::cerr << "  FAIL: Low-velocity friction incorrect" << std::endl;
        return false;
    }
    
    // At high velocity, should approach dynamic
    double mu_fast = svw.getFriction(sigma_n, 10.0, 0.0);
    if (mu_fast > 0.3) {  // Should be close to 0.2
        std::cerr << "  FAIL: High-velocity friction too high" << std::endl;
        return false;
    }
    
    // Should weaken monotonically with velocity
    double mu_mid = svw.getFriction(sigma_n, 0.5, 0.0);
    if (mu_mid >= mu_slow || mu_mid <= mu_fast) {
        std::cerr << "  FAIL: Friction not monotonically weakening" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Strong velocity weakening friction correct" << std::endl;
    return true;
}

/**
 * @test Test SVW slip-weakening distance
 */
bool test_svw_slip_weakening_distance() {
    std::cout << "Testing SVW slip-weakening distance..." << std::endl;
    
    StrongVelocityWeakeningFriction svw;
    svw.setStaticFriction(0.7);
    svw.setDynamicFriction(0.2);
    svw.setCriticalVelocity(0.1);
    
    double Dc = svw.getSlipWeakeningDistance(50e6);
    
    // Should be positive and reasonable
    if (Dc <= 0 || Dc > 10.0) {
        std::cerr << "  FAIL: Unreasonable Dc = " << Dc << std::endl;
        return false;
    }
    
    std::cout << "  PASS: SVW slip-weakening distance reasonable" << std::endl;
    return true;
}

/**
 * @test Test SVW fracture energy
 */
bool test_svw_fracture_energy() {
    std::cout << "Testing SVW fracture energy..." << std::endl;
    
    StrongVelocityWeakeningFriction svw;
    svw.setStaticFriction(0.7);
    svw.setDynamicFriction(0.2);
    svw.setCriticalVelocity(0.1);
    
    double sigma_n = 50e6;
    double Gc = svw.getFractureEnergy(sigma_n);
    
    // Should be positive
    if (Gc <= 0) {
        std::cerr << "  FAIL: Non-positive fracture energy" << std::endl;
        return false;
    }
    
    // Should be on order of MJ/m² for typical earthquakes
    if (Gc < 1e3 || Gc > 1e8) {
        std::cerr << "  FAIL: Unreasonable fracture energy Gc = " << Gc << std::endl;
        return false;
    }
    
    std::cout << "  PASS: SVW fracture energy reasonable" << std::endl;
    return true;
}

// =============================================================================
// Friction Model Factory Tests
// =============================================================================

/**
 * @test Test friction model factory
 */
bool test_friction_factory() {
    std::cout << "Testing friction model factory..." << std::endl;
    
    // Create different friction models
    auto coulomb = createFrictionModel(FrictionLaw::COULOMB);
    if (!coulomb) {
        std::cerr << "  FAIL: Could not create Coulomb model" << std::endl;
        return false;
    }
    
    auto rs_aging = createFrictionModel(FrictionLaw::RATE_STATE_AGING);
    if (!rs_aging) {
        std::cerr << "  FAIL: Could not create RS aging model" << std::endl;
        return false;
    }
    
    auto rs_slip = createFrictionModel(FrictionLaw::RATE_STATE_SLIP);
    if (!rs_slip) {
        std::cerr << "  FAIL: Could not create RS slip model" << std::endl;
        return false;
    }
    
    auto flash = createFrictionModel(FrictionLaw::FLASH_HEATING);
    if (!flash) {
        std::cerr << "  FAIL: Could not create flash heating model" << std::endl;
        return false;
    }
    
    auto tp = createFrictionModel(FrictionLaw::THERMAL_PRESSURIZATION);
    if (!tp) {
        std::cerr << "  FAIL: Could not create TP model" << std::endl;
        return false;
    }
    
    auto svw = createFrictionModel(FrictionLaw::STRONG_VELOCITY_WEAKENING);
    if (!svw) {
        std::cerr << "  FAIL: Could not create SVW model" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Friction model factory works" << std::endl;
    return true;
}

/**
 * @test Test friction law parsing
 */
bool test_friction_parsing() {
    std::cout << "Testing friction law parsing..." << std::endl;
    
    // Test parsing from strings
    auto law1 = parseFrictionLaw("coulomb");
    if (law1 != FrictionLaw::COULOMB) {
        std::cerr << "  FAIL: Could not parse 'coulomb'" << std::endl;
        return false;
    }
    
    auto law2 = parseFrictionLaw("rate_state_aging");
    if (law2 != FrictionLaw::RATE_STATE_AGING) {
        std::cerr << "  FAIL: Could not parse 'rate_state_aging'" << std::endl;
        return false;
    }
    
    auto law3 = parseFrictionLaw("FH");  // Alias for flash heating
    if (law3 != FrictionLaw::FLASH_HEATING) {
        std::cerr << "  FAIL: Could not parse 'FH'" << std::endl;
        return false;
    }
    
    auto law4 = parseFrictionLaw("SVW");  // Alias for SVW
    if (law4 != FrictionLaw::STRONG_VELOCITY_WEAKENING) {
        std::cerr << "  FAIL: Could not parse 'SVW'" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Friction law parsing works" << std::endl;
    return true;
}

// =============================================================================
// Test Runner
// =============================================================================

int runFrictionLawTests() {
    std::cout << "\n=== Friction Law Unit Tests ===" << std::endl;
    
    int passed = 0;
    int failed = 0;
    
    // Coulomb tests
    if (test_coulomb_static()) ++passed; else ++failed;
    if (test_coulomb_slip_weakening()) ++passed; else ++failed;
    
    // Rate-and-state tests
    if (test_rate_state_steady_state()) ++passed; else ++failed;
    if (test_rate_state_aging_law()) ++passed; else ++failed;
    if (test_rate_state_slip_law()) ++passed; else ++failed;
    
    // Flash heating tests
    if (test_flash_heating_weakening_velocity()) ++passed; else ++failed;
    if (test_flash_heating_friction()) ++passed; else ++failed;
    
    // Thermal pressurization tests
    if (test_thermal_pressurization_pressure()) ++passed; else ++failed;
    if (test_thermal_pressurization_friction()) ++passed; else ++failed;
    if (test_tp_characteristic_slip()) ++passed; else ++failed;
    
    // SVW tests
    if (test_svw_friction()) ++passed; else ++failed;
    if (test_svw_slip_weakening_distance()) ++passed; else ++failed;
    if (test_svw_fracture_energy()) ++passed; else ++failed;
    
    // Factory tests
    if (test_friction_factory()) ++passed; else ++failed;
    if (test_friction_parsing()) ++passed; else ++failed;
    
    std::cout << "\n=== Friction Law Test Summary ===" << std::endl;
    std::cout << "Passed: " << passed << std::endl;
    std::cout << "Failed: " << failed << std::endl;
    std::cout << "Total:  " << (passed + failed) << std::endl;
    
    return failed;
}

} // namespace Testing
} // namespace FSRM

#ifndef FSRM_TEST_NO_MAIN
int main() {
    return FSRM::Testing::runFrictionLawTests();
}
#endif
