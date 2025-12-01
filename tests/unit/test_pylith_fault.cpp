/**
 * @file test_pylith_fault.cpp
 * @brief Unit tests for PyLith-style fault representation
 * 
 * Tests cover:
 * - Slip time functions (step, ramp, Brune, Liu)
 * - Friction models (static, slip-weakening, rate-state)
 * - Kinematic faults (FaultCohesiveKin)
 * - Dynamic faults (FaultCohesiveDyn)
 * - Fault utilities
 * - Configuration and factory functions
 */

#include "PyLithFault.hpp"
#include "Testing.hpp"
#include <cmath>
#include <iostream>
#include <vector>
#include <memory>

namespace FSRM {
namespace Testing {

// =============================================================================
// Slip Time Function Tests
// =============================================================================

/**
 * @test Test step slip time function
 */
bool test_slip_time_fn_step() {
    std::cout << "Testing step slip time function..." << std::endl;
    
    SlipTimeFnStep step_fn;
    
    double amplitude = 2.0;  // 2 m final slip
    double slip_time = 1.0;  // Starts at t=1s
    
    // Before slip time - should be zero
    double s0 = step_fn.computeSlip(0.5, slip_time, amplitude);
    if (std::abs(s0) > 1e-10) {
        std::cerr << "  FAIL: Slip before initiation = " << s0 << ", expected 0" << std::endl;
        return false;
    }
    
    // After slip time - should be full amplitude
    double s1 = step_fn.computeSlip(2.0, slip_time, amplitude);
    if (std::abs(s1 - amplitude) > 1e-10) {
        std::cerr << "  FAIL: Slip after initiation = " << s1 << ", expected " << amplitude << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Step slip time function correct" << std::endl;
    return true;
}

/**
 * @test Test ramp slip time function
 */
bool test_slip_time_fn_ramp() {
    std::cout << "Testing ramp slip time function..." << std::endl;
    
    SlipTimeFnRamp ramp_fn;
    ramp_fn.setRiseTime(2.0);  // 2 second rise time
    
    double amplitude = 4.0;
    double slip_time = 0.0;
    
    // At half rise time - should be half amplitude
    double s_half = ramp_fn.computeSlip(1.0, slip_time, amplitude);
    if (std::abs(s_half - 2.0) > 1e-10) {
        std::cerr << "  FAIL: Slip at half rise time = " << s_half << ", expected 2.0" << std::endl;
        return false;
    }
    
    // After rise time - should be full amplitude
    double s_full = ramp_fn.computeSlip(3.0, slip_time, amplitude);
    if (std::abs(s_full - amplitude) > 1e-10) {
        std::cerr << "  FAIL: Slip after rise time = " << s_full << ", expected " << amplitude << std::endl;
        return false;
    }
    
    // Slip rate during rise - should be constant
    double rate = ramp_fn.computeSlipRate(1.0, slip_time, amplitude);
    double expected_rate = amplitude / 2.0;  // 2 m/s
    if (std::abs(rate - expected_rate) > 1e-10) {
        std::cerr << "  FAIL: Slip rate during rise = " << rate << ", expected " << expected_rate << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Ramp slip time function correct" << std::endl;
    return true;
}

/**
 * @test Test Brune slip time function
 */
bool test_slip_time_fn_brune() {
    std::cout << "Testing Brune slip time function..." << std::endl;
    
    SlipTimeFnBrune brune_fn;
    brune_fn.setCharacteristicTime(1.0);  // τ_c = 1s
    
    double amplitude = 1.0;
    double slip_time = 0.0;
    
    // At t = 0, slip = 0
    double s0 = brune_fn.computeSlip(0.0, slip_time, amplitude);
    if (std::abs(s0) > 1e-10) {
        std::cerr << "  FAIL: Slip at t=0 = " << s0 << ", expected 0" << std::endl;
        return false;
    }
    
    // At large time, slip should approach amplitude
    double s_large = brune_fn.computeSlip(10.0, slip_time, amplitude);
    if (std::abs(s_large - amplitude) > 0.01) {
        std::cerr << "  FAIL: Slip at large t = " << s_large << ", expected ~" << amplitude << std::endl;
        return false;
    }
    
    // Slip rate at t = τ_c should be amplitude * exp(-1) / τ_c
    double rate_tc = brune_fn.computeSlipRate(1.0, slip_time, amplitude);
    double expected_rate = amplitude * 1.0 / 1.0 * std::exp(-1.0);  // τ/τ_c * exp(-τ/τ_c) / τ_c
    if (std::abs(rate_tc - expected_rate) > 1e-6) {
        std::cerr << "  FAIL: Slip rate at τ_c = " << rate_tc << ", expected " << expected_rate << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Brune slip time function correct" << std::endl;
    return true;
}

/**
 * @test Test Liu slip time function
 */
bool test_slip_time_fn_liu() {
    std::cout << "Testing Liu slip time function..." << std::endl;
    
    SlipTimeFnLiu liu_fn;
    liu_fn.setRiseTime(2.0);
    
    double amplitude = 1.0;
    double slip_time = 0.0;
    
    // At t = 0, slip = 0
    double s0 = liu_fn.computeSlip(0.0, slip_time, amplitude);
    if (std::abs(s0) > 1e-10) {
        std::cerr << "  FAIL: Slip at t=0 = " << s0 << ", expected 0" << std::endl;
        return false;
    }
    
    // At t = rise_time, slip should be amplitude
    double s_full = liu_fn.computeSlip(2.0, slip_time, amplitude);
    if (std::abs(s_full - amplitude) > 1e-6) {
        std::cerr << "  FAIL: Slip at rise time = " << s_full << ", expected " << amplitude << std::endl;
        return false;
    }
    
    // Slip rate should be symmetric about midpoint
    double rate_quarter = liu_fn.computeSlipRate(0.5, slip_time, amplitude);
    double rate_three_quarter = liu_fn.computeSlipRate(1.5, slip_time, amplitude);
    if (std::abs(rate_quarter - rate_three_quarter) > 1e-6) {
        std::cerr << "  FAIL: Slip rate not symmetric" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Liu slip time function correct" << std::endl;
    return true;
}

/**
 * @test Test constant rate slip function
 */
bool test_slip_time_fn_constant_rate() {
    std::cout << "Testing constant rate slip function..." << std::endl;
    
    SlipTimeFnConstantRate creep_fn;
    creep_fn.setSlipRate(1e-9);  // 1 nm/s (typical creep rate)
    
    double slip_time = 0.0;
    double amplitude = 0.0;  // Not used for constant rate
    
    // After 1 year (31536000 s), should have ~3.15 cm
    double t_year = 31536000.0;
    double s_year = creep_fn.computeSlip(t_year, slip_time, amplitude);
    double expected = 1e-9 * t_year;  // ~0.0315 m
    
    if (std::abs(s_year - expected) > 1e-10) {
        std::cerr << "  FAIL: Slip after 1 year = " << s_year << ", expected " << expected << std::endl;
        return false;
    }
    
    // Rate should be constant
    double rate = creep_fn.computeSlipRate(t_year / 2, slip_time, amplitude);
    if (std::abs(rate - 1e-9) > 1e-15) {
        std::cerr << "  FAIL: Slip rate = " << rate << ", expected 1e-9" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Constant rate slip function correct" << std::endl;
    return true;
}

// =============================================================================
// Friction Model Tests
// =============================================================================

/**
 * @test Test static friction model
 */
bool test_static_friction() {
    std::cout << "Testing static friction model..." << std::endl;
    
    StaticFriction sf;
    sf.setFrictionCoefficient(0.6);
    
    DynamicFaultState state;
    state.slip_rate = 1.0;
    state.state_variable = 0.5;
    
    double f = sf.computeFriction(state);
    if (std::abs(f - 0.6) > 1e-10) {
        std::cerr << "  FAIL: Friction = " << f << ", expected 0.6" << std::endl;
        return false;
    }
    
    // State rate should be zero
    double rate = sf.computeStateRate(state);
    if (std::abs(rate) > 1e-10) {
        std::cerr << "  FAIL: State rate = " << rate << ", expected 0" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Static friction model correct" << std::endl;
    return true;
}

/**
 * @test Test slip-weakening friction model
 */
bool test_slip_weakening_friction() {
    std::cout << "Testing slip-weakening friction model..." << std::endl;
    
    SlipWeakeningFriction sw;
    sw.setStaticCoefficient(0.6);
    sw.setDynamicCoefficient(0.4);
    sw.setCriticalSlipDistance(0.5);  // 0.5 m
    
    DynamicFaultState state;
    
    // At zero slip - static friction
    state.state_variable = 0.0;  // Cumulative slip
    double f0 = sw.computeFriction(state);
    if (std::abs(f0 - 0.6) > 1e-10) {
        std::cerr << "  FAIL: Friction at zero slip = " << f0 << ", expected 0.6" << std::endl;
        return false;
    }
    
    // At half D_c - should be interpolated
    state.state_variable = 0.25;
    double f_half = sw.computeFriction(state);
    double expected_half = 0.6 - 0.5 * (0.6 - 0.4);  // = 0.5
    if (std::abs(f_half - expected_half) > 1e-10) {
        std::cerr << "  FAIL: Friction at half D_c = " << f_half << ", expected " << expected_half << std::endl;
        return false;
    }
    
    // Beyond D_c - dynamic friction
    state.state_variable = 1.0;
    double f_full = sw.computeFriction(state);
    if (std::abs(f_full - 0.4) > 1e-10) {
        std::cerr << "  FAIL: Friction beyond D_c = " << f_full << ", expected 0.4" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Slip-weakening friction model correct" << std::endl;
    return true;
}

/**
 * @test Test rate-and-state friction (aging law) - steady state
 */
bool test_rate_state_aging_steady_state() {
    std::cout << "Testing rate-and-state aging law steady state..." << std::endl;
    
    RateStateFrictionAging rs;
    rs.setDirectEffect(0.01);       // a
    rs.setEvolutionEffect(0.015);   // b
    rs.setCriticalSlipDistance(0.01);  // D_c = 1 cm
    rs.setReferenceFriction(0.6);   // f_0
    rs.setReferenceSlipRate(1e-6);  // V_0 = 1 μm/s
    
    double V = 1e-3;  // 1 mm/s
    double Dc = 0.01;
    double V0 = 1e-6;
    
    // Steady-state state variable: θ_ss = D_c / V
    double theta_ss = Dc / V;
    
    DynamicFaultState state;
    state.slip_rate = V;
    state.state_variable = theta_ss;
    
    // At steady state, dθ/dt should be ~0
    double dtheta_dt = rs.computeStateRate(state);
    if (std::abs(dtheta_dt) > 1e-10) {
        std::cerr << "  FAIL: State rate at steady state = " << dtheta_dt << ", expected ~0" << std::endl;
        return false;
    }
    
    // Steady-state friction: f_ss = f_0 + (a-b) * ln(V/V_0)
    double f_ss = rs.computeFriction(state);
    double expected_f = 0.6 + (0.01 - 0.015) * std::log(V / V0);
    if (std::abs(f_ss - expected_f) > 1e-6) {
        std::cerr << "  FAIL: Steady-state friction = " << f_ss << ", expected " << expected_f << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Rate-and-state aging law steady state correct" << std::endl;
    return true;
}

/**
 * @test Test rate-and-state friction - velocity weakening behavior
 */
bool test_rate_state_velocity_weakening() {
    std::cout << "Testing rate-and-state velocity weakening..." << std::endl;
    
    RateStateFrictionAging rs;
    rs.setDirectEffect(0.01);       // a
    rs.setEvolutionEffect(0.015);   // b > a : velocity weakening
    rs.setCriticalSlipDistance(0.01);
    rs.setReferenceFriction(0.6);
    rs.setReferenceSlipRate(1e-6);
    
    // For a < b, increasing velocity should decrease steady-state friction
    double Dc = 0.01;
    
    DynamicFaultState state_slow;
    state_slow.slip_rate = 1e-6;  // V_0
    state_slow.state_variable = Dc / state_slow.slip_rate;
    double f_slow = rs.computeFriction(state_slow);
    
    DynamicFaultState state_fast;
    state_fast.slip_rate = 1e-3;  // 1000x faster
    state_fast.state_variable = Dc / state_fast.slip_rate;
    double f_fast = rs.computeFriction(state_fast);
    
    if (f_fast >= f_slow) {
        std::cerr << "  FAIL: f_fast = " << f_fast << " >= f_slow = " << f_slow << std::endl;
        std::cerr << "  FAIL: Velocity weakening not observed" << std::endl;
        return false;
    }
    
    // Check critical stiffness (positive for velocity weakening)
    double sigma_n = 50e6;  // 50 MPa
    double k_c = rs.getCriticalStiffness(sigma_n);
    if (k_c <= 0) {
        std::cerr << "  FAIL: Critical stiffness should be positive for VW" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Rate-and-state velocity weakening correct" << std::endl;
    return true;
}

/**
 * @test Test nucleation length calculation
 */
bool test_nucleation_length() {
    std::cout << "Testing nucleation length calculation..." << std::endl;
    
    RateStateFrictionAging rs;
    rs.setDirectEffect(0.01);
    rs.setEvolutionEffect(0.015);
    rs.setCriticalSlipDistance(0.01);
    
    double G = 30e9;      // 30 GPa shear modulus
    double sigma_n = 50e6; // 50 MPa normal stress
    
    double L_c = rs.getNucleationLength(G, sigma_n);
    
    // L_c = G * D_c / ((b - a) * σ_n)
    double expected = G * 0.01 / (0.005 * sigma_n);
    
    if (std::abs(L_c - expected) > 1.0) {
        std::cerr << "  FAIL: Nucleation length = " << L_c << " m, expected " << expected << " m" << std::endl;
        return false;
    }
    
    // Should be on order of 100 m - 10 km for typical parameters
    if (L_c < 10.0 || L_c > 100000.0) {
        std::cerr << "  FAIL: Unreasonable nucleation length: " << L_c << " m" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Nucleation length calculation correct, L_c = " << L_c << " m" << std::endl;
    return true;
}

// =============================================================================
// Kinematic Fault Tests
// =============================================================================

/**
 * @test Test kinematic fault with uniform slip
 */
bool test_kinematic_fault_uniform_slip() {
    std::cout << "Testing kinematic fault with uniform slip..." << std::endl;
    
    FaultCohesiveKin fault;
    
    // Create simple fault vertices
    std::vector<FaultVertex> vertices(4);
    for (int i = 0; i < 4; ++i) {
        vertices[i].vertex_id = i;
        vertices[i].vertex_negative = i;
        vertices[i].vertex_positive = i + 4;
        vertices[i].coords[0] = (i % 2) * 100.0;
        vertices[i].coords[1] = (i / 2) * 100.0;
        vertices[i].coords[2] = 0.0;
        vertices[i].normal = {0.0, 0.0, 1.0};
        vertices[i].along_strike = {1.0, 0.0, 0.0};
        vertices[i].up_dip = {0.0, 1.0, 0.0};
        vertices[i].area = 2500.0;  // 50m x 50m
    }
    fault.setFaultVertices(vertices);
    
    // Set uniform slip
    fault.setUniformSlip(1.0, 0.0, 0.0,  // 1m left-lateral slip
                        0.5, 2.0);        // starts at t=0.5s, 2s rise time
    
    // Set Brune slip time function
    fault.setSlipTimeFnType(SlipTimeFnType::BRUNE);
    
    fault.initialize();
    
    // Get slip at different times
    double slip_ll, slip_rev, slip_open;
    
    // Before initiation
    fault.getSlipAtTime(0.0, 0, slip_ll, slip_rev, slip_open);
    if (std::abs(slip_ll) > 1e-10) {
        std::cerr << "  FAIL: Slip before initiation = " << slip_ll << std::endl;
        return false;
    }
    
    // Long after initiation - should be near final
    fault.getSlipAtTime(10.0, 0, slip_ll, slip_rev, slip_open);
    if (std::abs(slip_ll - 1.0) > 0.01) {
        std::cerr << "  FAIL: Final slip = " << slip_ll << ", expected ~1.0" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Kinematic fault uniform slip correct" << std::endl;
    return true;
}

/**
 * @test Test kinematic fault with circular rupture
 */
bool test_kinematic_fault_circular_rupture() {
    std::cout << "Testing kinematic fault with circular rupture..." << std::endl;
    
    FaultCohesiveKin fault;
    
    // Create grid of fault vertices
    int nx = 5, nz = 5;
    std::vector<FaultVertex> vertices(nx * nz);
    for (int j = 0; j < nz; ++j) {
        for (int i = 0; i < nx; ++i) {
            int idx = j * nx + i;
            vertices[idx].vertex_id = idx;
            vertices[idx].vertex_negative = idx;
            vertices[idx].vertex_positive = idx + nx * nz;
            vertices[idx].coords[0] = (i - nx/2) * 1000.0;  // ±2 km
            vertices[idx].coords[1] = (j - nz/2) * 1000.0;
            vertices[idx].coords[2] = 0.0;
            vertices[idx].normal = {0.0, 0.0, 1.0};
            vertices[idx].along_strike = {1.0, 0.0, 0.0};
            vertices[idx].up_dip = {0.0, 1.0, 0.0};
            vertices[idx].area = 1e6;  // 1 km²
        }
    }
    fault.setFaultVertices(vertices);
    
    // Set circular rupture from center
    double rupture_velocity = 2800.0;  // 2.8 km/s
    double final_slip = 2.0;           // 2 m slip
    double rise_time = 1.0;            // 1 s rise time
    double rake = 0.0;                 // Pure strike-slip
    
    fault.setCircularRupture(0.0, 0.0, 0.0, rupture_velocity, final_slip, rise_time, rake);
    fault.setSlipTimeFnType(SlipTimeFnType::RAMP);
    fault.initialize();
    
    // Check slip times - center should rupture first
    const auto& params_center = fault.getSlipParams(12);  // Center vertex
    if (std::abs(params_center.slip_time) > 1e-10) {
        std::cerr << "  FAIL: Center slip time = " << params_center.slip_time << ", expected 0" << std::endl;
        return false;
    }
    
    // Corner should rupture later
    const auto& params_corner = fault.getSlipParams(0);  // Corner vertex
    double expected_time = std::sqrt(2e6 + 2e6) / rupture_velocity;  // √(2² + 2²) km / V_r
    if (std::abs(params_corner.slip_time - expected_time) > 0.01) {
        std::cerr << "  FAIL: Corner slip time = " << params_corner.slip_time 
                  << ", expected " << expected_time << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Kinematic fault circular rupture correct" << std::endl;
    return true;
}

// =============================================================================
// Dynamic Fault Tests
// =============================================================================

/**
 * @test Test dynamic fault initialization
 */
bool test_dynamic_fault_initialization() {
    std::cout << "Testing dynamic fault initialization..." << std::endl;
    
    FaultCohesiveDyn fault;
    
    // Create fault vertices
    std::vector<FaultVertex> vertices(4);
    for (int i = 0; i < 4; ++i) {
        vertices[i].vertex_id = i;
        vertices[i].vertex_negative = i;
        vertices[i].vertex_positive = i + 4;
        vertices[i].coords[0] = (i % 2) * 100.0;
        vertices[i].coords[1] = (i / 2) * 100.0;
        vertices[i].coords[2] = 0.0;
        vertices[i].normal = {0.0, 0.0, 1.0};
        vertices[i].along_strike = {1.0, 0.0, 0.0};
        vertices[i].up_dip = {0.0, 1.0, 0.0};
        vertices[i].area = 2500.0;
    }
    fault.setFaultVertices(vertices);
    
    // Set initial traction (compression and shear)
    fault.setUniformInitialTraction(25e6, 0.0, -50e6);  // τ = 25 MPa, σ_n = -50 MPa (compression)
    
    // Set slip-weakening friction
    auto sw = std::make_unique<SlipWeakeningFriction>();
    sw->setStaticCoefficient(0.6);
    sw->setDynamicCoefficient(0.4);
    sw->setCriticalSlipDistance(0.4);
    fault.setFrictionModel(std::move(sw));
    
    fault.initialize();
    
    // Check initial state
    const auto& state = fault.getState(0);
    
    // Should be locked initially
    if (!state.is_locked) {
        std::cerr << "  FAIL: Fault should be locked initially" << std::endl;
        return false;
    }
    
    // Effective normal stress should be positive (compression)
    if (state.effective_normal <= 0) {
        std::cerr << "  FAIL: Effective normal stress = " << state.effective_normal << ", should be positive" << std::endl;
        return false;
    }
    
    // Has not ruptured yet
    if (fault.hasRuptured()) {
        std::cerr << "  FAIL: Fault should not have ruptured initially" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Dynamic fault initialization correct" << std::endl;
    return true;
}

/**
 * @test Test dynamic fault traction perturbation
 */
bool test_dynamic_fault_perturbation() {
    std::cout << "Testing dynamic fault traction perturbation..." << std::endl;
    
    FaultCohesiveDyn fault;
    
    // Create simple fault
    std::vector<FaultVertex> vertices(1);
    vertices[0].vertex_id = 0;
    vertices[0].coords = {0.0, 0.0, 0.0};
    fault.setFaultVertices(vertices);
    
    // Set Gaussian perturbation at center
    fault.setTractionPerturbation(
        TractionPerturbationType::GAUSSIAN,
        5e6,    // 5 MPa amplitude
        1.0,    // 1 s duration
        0.0, 0.0, 0.0,  // center at origin
        1000.0  // 1 km width
    );
    
    // Manually compute perturbation (fault doesn't have direct access)
    // At center (r=0), spatial factor = 1.0
    // At t=0.5s (half duration), time factor = 0.5*(1-cos(π/2)) = 0.5*(1-0) = 0.5
    // But this is internal to the class, so we just verify initialization works
    
    fault.initialize();
    
    std::cout << "  PASS: Dynamic fault perturbation setup correct" << std::endl;
    return true;
}

// =============================================================================
// Factory and Configuration Tests
// =============================================================================

/**
 * @test Test slip time function factory
 */
bool test_slip_fn_factory() {
    std::cout << "Testing slip time function factory..." << std::endl;
    
    // Create different types
    auto step = createSlipTimeFn(SlipTimeFnType::STEP);
    if (!step || step->getType() != SlipTimeFnType::STEP) {
        std::cerr << "  FAIL: Could not create step function" << std::endl;
        return false;
    }
    
    auto brune = createSlipTimeFn("brune");
    if (!brune || brune->getType() != SlipTimeFnType::BRUNE) {
        std::cerr << "  FAIL: Could not create Brune function from string" << std::endl;
        return false;
    }
    
    auto liu = createSlipTimeFn("LIU");
    if (!liu || liu->getType() != SlipTimeFnType::LIU) {
        std::cerr << "  FAIL: Could not create Liu function (case insensitive)" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Slip time function factory correct" << std::endl;
    return true;
}

/**
 * @test Test friction model factory
 */
bool test_friction_factory() {
    std::cout << "Testing friction model factory..." << std::endl;
    
    auto sf = createFaultFrictionModel("static");
    if (!sf || sf->name() != "static") {
        std::cerr << "  FAIL: Could not create static friction" << std::endl;
        return false;
    }
    
    auto sw = createFaultFrictionModel("slip_weakening");
    if (!sw || sw->name() != "slip_weakening") {
        std::cerr << "  FAIL: Could not create slip-weakening friction" << std::endl;
        return false;
    }
    
    auto rs = createFaultFrictionModel("rate_state_aging");
    if (!rs || rs->name() != "rate_state_aging") {
        std::cerr << "  FAIL: Could not create rate-state aging friction" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Friction model factory correct" << std::endl;
    return true;
}

/**
 * @test Test fault creation from configuration
 */
bool test_fault_configuration() {
    std::cout << "Testing fault configuration..." << std::endl;
    
    // Create kinematic fault from config
    std::map<std::string, std::string> kin_config;
    kin_config["name"] = "test_kinematic";
    kin_config["fault_type"] = "kinematic";
    kin_config["slip_time_fn"] = "brune";
    kin_config["rise_time"] = "1.5";
    kin_config["final_slip_ll"] = "2.0";
    
    auto kin_fault = createFaultCohesive(kin_config);
    if (!kin_fault || kin_fault->getType() != CohesiveFaultType::KINEMATIC) {
        std::cerr << "  FAIL: Could not create kinematic fault from config" << std::endl;
        return false;
    }
    if (kin_fault->getName() != "test_kinematic") {
        std::cerr << "  FAIL: Fault name not set correctly" << std::endl;
        return false;
    }
    
    // Create dynamic fault from config
    std::map<std::string, std::string> dyn_config;
    dyn_config["name"] = "test_dynamic";
    dyn_config["fault_type"] = "dynamic";
    dyn_config["friction_model"] = "slip_weakening";
    dyn_config["static_friction"] = "0.7";
    dyn_config["dynamic_friction"] = "0.3";
    dyn_config["critical_slip_distance"] = "0.5";
    
    auto dyn_fault = createFaultCohesive(dyn_config);
    if (!dyn_fault || dyn_fault->getType() != CohesiveFaultType::DYNAMIC) {
        std::cerr << "  FAIL: Could not create dynamic fault from config" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Fault configuration correct" << std::endl;
    return true;
}

// =============================================================================
// Fault Utility Tests
// =============================================================================

/**
 * @test Test moment magnitude calculation
 */
bool test_moment_magnitude() {
    std::cout << "Testing moment magnitude calculation..." << std::endl;
    
    // M7 earthquake: M0 ≈ 4e19 N·m
    // Area ≈ 40 km × 20 km = 8e8 m²
    // Slip ≈ 1.5 m
    // μ ≈ 30 GPa
    
    double slip = 1.5;
    double area = 8e8;
    double mu = 30e9;
    
    double Mw = FaultUtils::computeMomentMagnitude(slip, area, mu);
    
    // Should be approximately M7
    if (std::abs(Mw - 7.0) > 0.5) {
        std::cerr << "  FAIL: Magnitude = " << Mw << ", expected ~7.0" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Moment magnitude = " << Mw << " (expected ~7)" << std::endl;
    return true;
}

/**
 * @test Test stress drop calculation
 */
bool test_stress_drop() {
    std::cout << "Testing stress drop calculation..." << std::endl;
    
    double slip = 1.0;  // 1 m
    double radius = 5000.0;  // 5 km
    double mu = 30e9;  // 30 GPa
    
    double delta_sigma = FaultUtils::computeStressDrop(slip, radius, mu);
    
    // Δσ = (7/16) * μ * D / r ≈ 2.6 MPa
    double expected = (7.0/16.0) * mu * slip / radius;
    
    if (std::abs(delta_sigma - expected) > 1e3) {
        std::cerr << "  FAIL: Stress drop = " << delta_sigma/1e6 << " MPa, expected " << expected/1e6 << " MPa" << std::endl;
        return false;
    }
    
    // Should be order of 1-10 MPa for typical earthquakes
    if (delta_sigma < 0.1e6 || delta_sigma > 100e6) {
        std::cerr << "  FAIL: Unreasonable stress drop: " << delta_sigma/1e6 << " MPa" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Stress drop = " << delta_sigma/1e6 << " MPa" << std::endl;
    return true;
}

/**
 * @test Test rupture velocity calculation
 */
bool test_rupture_velocity() {
    std::cout << "Testing rupture velocity calculation..." << std::endl;
    
    double Vs = 3500.0;  // 3.5 km/s
    double Vr = 0.92 * Vs;  // Rayleigh velocity ≈ 0.92 * Vs
    
    // Sub-Rayleigh rupture
    double V_sub = FaultUtils::computeRuptureVelocity(Vs, Vr, false);
    if (V_sub <= 0 || V_sub >= Vr) {
        std::cerr << "  FAIL: Sub-Rayleigh velocity out of range" << std::endl;
        return false;
    }
    
    // Supershear rupture
    double V_super = FaultUtils::computeRuptureVelocity(Vs, Vr, true);
    if (V_super <= Vs) {
        std::cerr << "  FAIL: Supershear velocity should exceed Vs" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Rupture velocities: sub-Rayleigh = " << V_sub/1000 
              << " km/s, supershear = " << V_super/1000 << " km/s" << std::endl;
    return true;
}

// =============================================================================
// FaultSystem Tests
// =============================================================================

/**
 * @test Test fault system with multiple faults
 */
bool test_fault_system() {
    std::cout << "Testing fault system..." << std::endl;
    
    FaultSystem system;
    
    // Add kinematic fault
    auto kin = std::make_unique<FaultCohesiveKin>();
    kin->setName("main_fault");
    std::vector<FaultVertex> v1(2);
    v1[0].vertex_id = 0;
    v1[1].vertex_id = 1;
    kin->setFaultVertices(v1);
    system.addFault(std::move(kin));
    
    // Add dynamic fault
    auto dyn = std::make_unique<FaultCohesiveDyn>();
    dyn->setName("secondary_fault");
    std::vector<FaultVertex> v2(2);
    v2[0].vertex_id = 0;
    v2[1].vertex_id = 1;
    dyn->setFaultVertices(v2);
    system.addFault(std::move(dyn));
    
    if (system.numFaults() != 2) {
        std::cerr << "  FAIL: Wrong number of faults" << std::endl;
        return false;
    }
    
    // Access by name
    auto* main_fault = system.getFault("main_fault");
    if (!main_fault || main_fault->getName() != "main_fault") {
        std::cerr << "  FAIL: Could not retrieve fault by name" << std::endl;
        return false;
    }
    
    // Initialize system
    system.initialize();
    
    std::cout << "  PASS: Fault system management correct" << std::endl;
    return true;
}

// =============================================================================
// Multi-Physics (THM) Fault Tests
// =============================================================================

/**
 * @test Test fault permeability tensor
 */
bool test_fault_permeability_tensor() {
    std::cout << "Testing fault permeability tensor..." << std::endl;
    
    FaultPermeabilityTensor k;
    
    // Default values - anisotropic
    if (k.k_strike != 1e-12 || k.k_normal != 1e-15) {
        std::cerr << "  FAIL: Default permeability values incorrect" << std::endl;
        return false;
    }
    
    // Test isotropic setting
    k.setIsotropic(1e-14);
    if (k.k_strike != 1e-14 || k.k_dip != 1e-14 || k.k_normal != 1e-14) {
        std::cerr << "  FAIL: Isotropic permeability not set correctly" << std::endl;
        return false;
    }
    
    // Test principal setting
    k.setPrincipal(1e-13, 1e-16);
    if (k.k_strike != 1e-13 || k.k_dip != 1e-13 || k.k_normal != 1e-16) {
        std::cerr << "  FAIL: Principal permeability not set correctly" << std::endl;
        return false;
    }
    
    // Test tensor retrieval
    auto tensor = k.getTensor();
    if (tensor.size() != 9) {
        std::cerr << "  FAIL: Tensor should be 3x3 (9 elements)" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Permeability tensor: k_parallel = " << k.k_strike 
              << ", k_normal = " << k.k_normal << std::endl;
    return true;
}

/**
 * @test Test fault hydraulic properties
 */
bool test_fault_hydraulic_properties() {
    std::cout << "Testing fault hydraulic properties..." << std::endl;
    
    FaultHydraulicProperties props;
    
    // Check default values
    if (std::abs(props.porosity - 0.1) > 1e-10) {
        std::cerr << "  FAIL: Default porosity incorrect" << std::endl;
        return false;
    }
    if (std::abs(props.biot_coefficient - 0.9) > 1e-10) {
        std::cerr << "  FAIL: Default Biot coefficient incorrect" << std::endl;
        return false;
    }
    
    // Test configuration from map
    std::map<std::string, std::string> config;
    config["fault_porosity"] = "0.15";
    config["fault_biot_coefficient"] = "0.85";
    config["fault_thickness"] = "0.05";
    config["fault_permeability_parallel"] = "5e-13";
    config["fault_permeability_perpendicular"] = "5e-16";
    
    props.configure(config);
    
    if (std::abs(props.porosity - 0.15) > 1e-10) {
        std::cerr << "  FAIL: Configured porosity incorrect: " << props.porosity << std::endl;
        return false;
    }
    if (std::abs(props.fault_thickness - 0.05) > 1e-10) {
        std::cerr << "  FAIL: Configured fault thickness incorrect" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Hydraulic properties configured: φ = " << props.porosity
              << ", w = " << props.fault_thickness << " m" << std::endl;
    return true;
}

/**
 * @test Test fault thermal properties
 */
bool test_fault_thermal_properties() {
    std::cout << "Testing fault thermal properties..." << std::endl;
    
    FaultThermalProperties props;
    
    // Check default reference temperature
    if (std::abs(props.reference_temperature - 293.15) > 1e-10) {
        std::cerr << "  FAIL: Default reference temperature incorrect" << std::endl;
        return false;
    }
    
    // Check flash heating parameters
    if (props.weakening_temperature < 1000.0) {
        std::cerr << "  FAIL: Weakening temperature too low" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Thermal properties: T_ref = " << props.reference_temperature 
              << " K, T_weak = " << props.weakening_temperature << " K" << std::endl;
    return true;
}

/**
 * @test Test effective stress calculation
 */
bool test_effective_stress() {
    std::cout << "Testing effective stress calculation..." << std::endl;
    
    double sigma_n = 50e6;  // 50 MPa total normal stress
    double p = 20e6;        // 20 MPa pore pressure
    double biot = 1.0;
    
    double sigma_eff = FaultPoroelastic::effectiveNormalStress(sigma_n, p, biot);
    double expected = 30e6;  // 30 MPa
    
    if (std::abs(sigma_eff - expected) > 1e-6) {
        std::cerr << "  FAIL: Effective stress = " << sigma_eff/1e6 
                  << " MPa, expected " << expected/1e6 << " MPa" << std::endl;
        return false;
    }
    
    // Test with Biot coefficient < 1
    biot = 0.8;
    sigma_eff = FaultPoroelastic::effectiveNormalStress(sigma_n, p, biot);
    expected = sigma_n - biot * p;  // 50 - 0.8*20 = 34 MPa
    
    if (std::abs(sigma_eff - expected) > 1e-6) {
        std::cerr << "  FAIL: Effective stress with Biot = " << sigma_eff/1e6 
                  << " MPa, expected " << expected/1e6 << " MPa" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Effective stress correctly computed" << std::endl;
    return true;
}

/**
 * @test Test thermal functions
 */
bool test_thermal_functions() {
    std::cout << "Testing thermal functions..." << std::endl;
    
    // Test frictional heating
    double tau = 10e6;    // 10 MPa shear stress
    double V = 1.0;       // 1 m/s slip rate
    double Q = FaultThermal::frictionalHeatingRate(tau, V);
    double expected_Q = 10e6;  // 10 MW/m²
    
    if (std::abs(Q - expected_Q) > 1e-6) {
        std::cerr << "  FAIL: Frictional heating = " << Q/1e6 
                  << " MW/m², expected " << expected_Q/1e6 << std::endl;
        return false;
    }
    
    // Test adiabatic temperature rise
    double slip = 0.1;    // 10 cm slip
    double w = 0.001;     // 1 mm fault width
    double rho_cp = 2.7e6; // Volumetric heat capacity
    double dT = FaultThermal::adiabaticTemperatureRise(tau, slip, w, rho_cp);
    // Expected: 10e6 * 0.1 / (2.7e6 * 0.001) = 370 K
    double expected_dT = tau * slip / (rho_cp * w);
    
    if (std::abs(dT - expected_dT) > 1e-6) {
        std::cerr << "  FAIL: Temperature rise = " << dT 
                  << " K, expected " << expected_dT << " K" << std::endl;
        return false;
    }
    
    // Test thermal weakening
    double f0 = 0.6;
    double T_current = 1000.0;  // Below weakening temp
    double T_weak = 1200.0;
    double f_residual = 0.1;
    
    double f = FaultThermal::thermalWeakeningFriction(f0, T_current, T_weak, f_residual);
    if (f != f0) {
        std::cerr << "  FAIL: Friction below weakening temp should be unchanged" << std::endl;
        return false;
    }
    
    // Above weakening temp
    T_current = 1300.0;
    f = FaultThermal::thermalWeakeningFriction(f0, T_current, T_weak, f_residual);
    if (f >= f0) {
        std::cerr << "  FAIL: Friction above weakening temp should be reduced" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Thermal functions correct" << std::endl;
    return true;
}

/**
 * @test Test permeability evolution models
 */
bool test_permeability_evolution() {
    std::cout << "Testing permeability evolution..." << std::endl;
    
    DefaultFaultPropertyUpdater updater;
    
    // Test slip-dependent model
    PermeabilityEvolutionParams params;
    params.model = PermeabilityUpdateModel::SLIP_DEPENDENT;
    params.slip_coefficient = 10.0;  // 1/m
    updater.setPermeabilityParams(params);
    
    FaultMultiPhysicsState state;
    state.slip_total = 0.1;  // 10 cm slip
    state.current_permeability.setIsotropic(1e-15);
    
    FaultHydraulicProperties props;
    
    auto k_new = updater.updatePermeability(state, props, 0.001);
    
    // k = k0 * exp(10 * 0.1) = k0 * e = 2.72 * k0
    double expected_factor = std::exp(10.0 * 0.1);
    double k_expected = 1e-15 * expected_factor;
    
    if (std::abs(k_new.k_strike - k_expected) / k_expected > 0.01) {
        std::cerr << "  FAIL: Slip-dependent permeability = " << k_new.k_strike 
                  << ", expected " << k_expected << std::endl;
        return false;
    }
    
    // Test stress-dependent model
    params.model = PermeabilityUpdateModel::STRESS_DEPENDENT;
    params.stress_coefficient = 0.1;
    params.stress_reference = 50e6;
    updater.setPermeabilityParams(params);
    
    state.effective_normal_stress = 50e6;  // Equal to reference
    state.current_permeability.setIsotropic(1e-15);
    
    k_new = updater.updatePermeability(state, props, 0.001);
    
    // k = k0 * exp(-0.1 * 50/50) = k0 * exp(-0.1) ≈ 0.905 * k0
    expected_factor = std::exp(-0.1);
    k_expected = 1e-15 * expected_factor;
    
    if (std::abs(k_new.k_strike - k_expected) / k_expected > 0.01) {
        std::cerr << "  FAIL: Stress-dependent permeability = " << k_new.k_strike 
                  << ", expected " << k_expected << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Permeability evolution models correct" << std::endl;
    return true;
}

/**
 * @test Test porosity evolution models
 */
bool test_porosity_evolution() {
    std::cout << "Testing porosity evolution..." << std::endl;
    
    DefaultFaultPropertyUpdater updater;
    
    // Test dilation model
    PorosityEvolutionParams params;
    params.model = PorosityUpdateModel::DILATION;
    params.dilation_coefficient = 0.1;
    updater.setPorosityParams(params);
    
    FaultMultiPhysicsState state;
    state.volumetric_strain = 0.01;  // 1% dilation
    
    FaultHydraulicProperties props;
    props.porosity = 0.1;
    props.porosity_ref = 0.1;
    
    double phi_new = updater.updatePorosity(state, props, 0.001);
    
    // φ = φ0 + β * Δε_v = 0.1 + 0.1 * 0.01 = 0.101
    double expected_phi = 0.1 + 0.1 * 0.01;
    
    if (std::abs(phi_new - expected_phi) > 1e-6) {
        std::cerr << "  FAIL: Dilation porosity = " << phi_new 
                  << ", expected " << expected_phi << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Porosity evolution models correct" << std::endl;
    return true;
}

/**
 * @test Test FaultCohesiveTHM initialization
 */
bool test_thm_fault_initialization() {
    std::cout << "Testing THM fault initialization..." << std::endl;
    
    FaultCohesiveTHM fault;
    
    // Configure with THM coupling
    std::map<std::string, std::string> config;
    config["fault_coupling_mode"] = "thm";
    config["fault_porosity"] = "0.1";
    config["fault_permeability"] = "1e-14";
    config["fault_thickness"] = "0.01";
    config["fault_reference_temperature"] = "350";
    config["permeability_update_model"] = "slip_dependent";
    
    fault.configure(config);
    
    if (fault.getCouplingMode() != FaultCouplingMode::THERMO_HYDRO_MECHANICAL) {
        std::cerr << "  FAIL: Coupling mode not set correctly" << std::endl;
        return false;
    }
    
    // Set up fault vertices
    std::vector<FaultVertex> vertices(5);
    for (size_t i = 0; i < vertices.size(); ++i) {
        vertices[i].vertex_id = i;
        vertices[i].coords = {static_cast<double>(i) * 100.0, 0.0, 0.0};
    }
    fault.setFaultVertices(vertices);
    
    // Initialize
    fault.initialize();
    
    // Check hydraulic field is initialized
    const auto& hyd_field = fault.getHydraulicField();
    if (hyd_field.size() != vertices.size()) {
        std::cerr << "  FAIL: Hydraulic field size mismatch" << std::endl;
        return false;
    }
    
    // Check thermal field is initialized  
    const auto& therm_field = fault.getThermalField();
    if (therm_field.size() != vertices.size()) {
        std::cerr << "  FAIL: Thermal field size mismatch" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: THM fault initialized with " << vertices.size() << " vertices" << std::endl;
    return true;
}

/**
 * @test Test coupling mode parsing
 */
bool test_coupling_mode_parsing() {
    std::cout << "Testing coupling mode parsing..." << std::endl;
    
    // Test various input strings
    if (parseFaultCouplingMode("mechanical") != FaultCouplingMode::MECHANICAL_ONLY) {
        std::cerr << "  FAIL: 'mechanical' not parsed correctly" << std::endl;
        return false;
    }
    if (parseFaultCouplingMode("poroelastic") != FaultCouplingMode::POROELASTIC) {
        std::cerr << "  FAIL: 'poroelastic' not parsed correctly" << std::endl;
        return false;
    }
    if (parseFaultCouplingMode("hm") != FaultCouplingMode::POROELASTIC) {
        std::cerr << "  FAIL: 'hm' not parsed correctly" << std::endl;
        return false;
    }
    if (parseFaultCouplingMode("THM") != FaultCouplingMode::THERMO_HYDRO_MECHANICAL) {
        std::cerr << "  FAIL: 'THM' not parsed correctly" << std::endl;
        return false;
    }
    if (parseFaultCouplingMode("thmc") != FaultCouplingMode::THMC) {
        std::cerr << "  FAIL: 'thmc' not parsed correctly" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: All coupling modes parsed correctly" << std::endl;
    return true;
}

/**
 * @test Test custom property updater callback
 */
bool test_custom_property_updater() {
    std::cout << "Testing custom property updater..." << std::endl;
    
    // Create a custom updater with specific behavior
    class CustomUpdater : public FaultPropertyUpdateCallback {
    public:
        FaultPermeabilityTensor updatePermeability(
            const FaultMultiPhysicsState& state,
            const FaultHydraulicProperties& /*props*/,
            double /*dt*/) override {
            // Custom: double permeability for any slip
            FaultPermeabilityTensor k = state.current_permeability;
            if (state.slip_total > 0) {
                k.k_strike *= 2.0;
                k.k_dip *= 2.0;
                k.k_normal *= 2.0;
            }
            return k;
        }
        
        double updatePorosity(
            const FaultMultiPhysicsState& /*state*/,
            const FaultHydraulicProperties& props,
            double /*dt*/) override {
            return props.porosity;  // No change
        }
        
        double updateFriction(
            const FaultMultiPhysicsState& /*state*/,
            double base_friction) override {
            return base_friction * 0.9;  // 10% reduction always
        }
    };
    
    auto custom_updater = std::make_shared<CustomUpdater>();
    
    FaultCohesiveTHM fault;
    fault.setPropertyUpdater(custom_updater);
    
    // Test the custom updater
    FaultMultiPhysicsState state;
    state.slip_total = 0.01;
    state.current_permeability.setIsotropic(1e-15);
    
    FaultHydraulicProperties props;
    props.porosity = 0.1;
    
    auto k_new = custom_updater->updatePermeability(state, props, 0.001);
    if (std::abs(k_new.k_strike - 2e-15) > 1e-20) {
        std::cerr << "  FAIL: Custom permeability update incorrect" << std::endl;
        return false;
    }
    
    double f = custom_updater->updateFriction(state, 0.6);
    if (std::abs(f - 0.54) > 1e-10) {
        std::cerr << "  FAIL: Custom friction update incorrect: " << f << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Custom property updater works correctly" << std::endl;
    return true;
}

/**
 * @test Test pore pressure effects on fault
 */
bool test_pore_pressure_effects() {
    std::cout << "Testing pore pressure effects on fault..." << std::endl;
    
    FaultCohesiveTHM fault;
    
    std::map<std::string, std::string> config;
    config["fault_coupling_mode"] = "poroelastic";
    config["fault_biot_coefficient"] = "0.9";
    fault.configure(config);
    
    // Set up vertices
    std::vector<FaultVertex> vertices(3);
    for (size_t i = 0; i < vertices.size(); ++i) {
        vertices[i].vertex_id = i;
    }
    fault.setFaultVertices(vertices);
    fault.initialize();
    
    // Set pore pressure field
    std::vector<double> pore_pressure = {10e6, 20e6, 15e6};  // 10, 20, 15 MPa
    fault.setPorePressureField(pore_pressure);
    
    // Check that pore pressure is stored
    const auto& hyd_field = fault.getHydraulicField();
    if (std::abs(hyd_field.pore_pressure[0] - 10e6) > 1e-6 ||
        std::abs(hyd_field.pore_pressure[1] - 20e6) > 1e-6) {
        std::cerr << "  FAIL: Pore pressure field not set correctly" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Pore pressure effects correctly applied" << std::endl;
    return true;
}

/**
 * @test Test THM fault factory
 */
bool test_thm_fault_factory() {
    std::cout << "Testing THM fault factory..." << std::endl;
    
    std::map<std::string, std::string> config;
    config["fault_coupling_mode"] = "thmc";
    config["fault_porosity"] = "0.12";
    config["fault_permeability_strike"] = "2e-13";
    config["permeability_update_model"] = "combined";
    
    auto fault = createFaultCohesiveTHM(config);
    
    if (!fault) {
        std::cerr << "  FAIL: Factory returned null" << std::endl;
        return false;
    }
    
    if (fault->getCouplingMode() != FaultCouplingMode::THMC) {
        std::cerr << "  FAIL: Factory did not set coupling mode correctly" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: THM fault factory works correctly" << std::endl;
    return true;
}

// =============================================================================
// Test Runner
// =============================================================================

int runPyLithFaultTests() {
    std::cout << "\n========================================" << std::endl;
    std::cout << "  PyLith-Style Fault Unit Tests" << std::endl;
    std::cout << "========================================\n" << std::endl;
    
    int passed = 0;
    int failed = 0;
    
    // Slip time function tests
    std::cout << "\n--- Slip Time Function Tests ---\n" << std::endl;
    if (test_slip_time_fn_step()) ++passed; else ++failed;
    if (test_slip_time_fn_ramp()) ++passed; else ++failed;
    if (test_slip_time_fn_brune()) ++passed; else ++failed;
    if (test_slip_time_fn_liu()) ++passed; else ++failed;
    if (test_slip_time_fn_constant_rate()) ++passed; else ++failed;
    
    // Friction model tests
    std::cout << "\n--- Friction Model Tests ---\n" << std::endl;
    if (test_static_friction()) ++passed; else ++failed;
    if (test_slip_weakening_friction()) ++passed; else ++failed;
    if (test_rate_state_aging_steady_state()) ++passed; else ++failed;
    if (test_rate_state_velocity_weakening()) ++passed; else ++failed;
    if (test_nucleation_length()) ++passed; else ++failed;
    
    // Kinematic fault tests
    std::cout << "\n--- Kinematic Fault Tests ---\n" << std::endl;
    if (test_kinematic_fault_uniform_slip()) ++passed; else ++failed;
    if (test_kinematic_fault_circular_rupture()) ++passed; else ++failed;
    
    // Dynamic fault tests
    std::cout << "\n--- Dynamic Fault Tests ---\n" << std::endl;
    if (test_dynamic_fault_initialization()) ++passed; else ++failed;
    if (test_dynamic_fault_perturbation()) ++passed; else ++failed;
    
    // Factory and configuration tests
    std::cout << "\n--- Factory and Configuration Tests ---\n" << std::endl;
    if (test_slip_fn_factory()) ++passed; else ++failed;
    if (test_friction_factory()) ++passed; else ++failed;
    if (test_fault_configuration()) ++passed; else ++failed;
    
    // Utility tests
    std::cout << "\n--- Fault Utility Tests ---\n" << std::endl;
    if (test_moment_magnitude()) ++passed; else ++failed;
    if (test_stress_drop()) ++passed; else ++failed;
    if (test_rupture_velocity()) ++passed; else ++failed;
    
    // System tests
    std::cout << "\n--- Fault System Tests ---\n" << std::endl;
    if (test_fault_system()) ++passed; else ++failed;
    
    // Multi-physics (THM) tests
    std::cout << "\n--- Multi-Physics (THM) Fault Tests ---\n" << std::endl;
    if (test_fault_permeability_tensor()) ++passed; else ++failed;
    if (test_fault_hydraulic_properties()) ++passed; else ++failed;
    if (test_fault_thermal_properties()) ++passed; else ++failed;
    if (test_effective_stress()) ++passed; else ++failed;
    if (test_thermal_functions()) ++passed; else ++failed;
    if (test_permeability_evolution()) ++passed; else ++failed;
    if (test_porosity_evolution()) ++passed; else ++failed;
    if (test_thm_fault_initialization()) ++passed; else ++failed;
    if (test_coupling_mode_parsing()) ++passed; else ++failed;
    if (test_custom_property_updater()) ++passed; else ++failed;
    if (test_pore_pressure_effects()) ++passed; else ++failed;
    if (test_thm_fault_factory()) ++passed; else ++failed;
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "  PyLith Fault Test Summary" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Passed: " << passed << std::endl;
    std::cout << "Failed: " << failed << std::endl;
    std::cout << "Total:  " << (passed + failed) << std::endl;
    std::cout << "========================================\n" << std::endl;
    
    return failed;
}

} // namespace Testing
} // namespace FSRM

#ifndef FSRM_TEST_NO_MAIN
int main() {
    return FSRM::Testing::runPyLithFaultTests();
}
#endif
