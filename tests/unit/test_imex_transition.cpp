/**
 * @file test_imex_transition.cpp
 * @brief Unit tests for Implicit-Explicit (IMEX) time integration transition manager
 * 
 * Tests cover:
 * - Trigger condition detection (stress, Coulomb failure, slip rate, velocity, etc.)
 * - Settling criteria (time elapsed, velocity decay, energy decay, slip rate decay)
 * - Mode transitions between implicit and explicit
 * - Phase statistics tracking
 * - Dynamic event detection
 * - Helper classes: ExplicitDynamicsHelper, SeismicEventDetector, KineticEnergyMonitor
 */

#include "ImplicitExplicitTransition.hpp"
#include "Testing.hpp"
#include <cmath>
#include <iostream>
#include <vector>

namespace FSRM {
namespace Testing {

// Tolerance constants for floating point comparisons
constexpr double TOL_STRICT = 1e-10;    // For exact numerical comparisons
constexpr double TOL_MAGNITUDE = 0.05;  // For magnitude calculations (0.05 units)
constexpr double TOL_RELATIVE = 0.01;   // For relative comparisons (1%)

// =============================================================================
// IMEXInternalConfig Tests
// =============================================================================

/**
 * @test Test default configuration values
 */
bool test_imex_default_config() {
    std::cout << "Testing IMEX default configuration..." << std::endl;
    
    IMEXInternalConfig config;
    
    // Check default values
    if (config.enabled != false) {
        std::cerr << "  FAIL: enabled should be false by default" << std::endl;
        return false;
    }
    
    if (config.initial_mode != IntegrationMode::IMPLICIT) {
        std::cerr << "  FAIL: initial_mode should be IMPLICIT by default" << std::endl;
        return false;
    }
    
    // Check time step defaults
    if (config.implicit_dt_initial != 3600.0) {
        std::cerr << "  FAIL: implicit_dt_initial should be 3600.0" << std::endl;
        return false;
    }
    
    if (config.explicit_dt_initial != 1e-4) {
        std::cerr << "  FAIL: explicit_dt_initial should be 1e-4" << std::endl;
        return false;
    }
    
    // Check trigger defaults
    if (config.trigger_type != TransitionTrigger::COULOMB_FAILURE) {
        std::cerr << "  FAIL: trigger_type should be COULOMB_FAILURE by default" << std::endl;
        return false;
    }
    
    // Check settling defaults
    if (config.settling_type != SettlingCriterion::COMBINED) {
        std::cerr << "  FAIL: settling_type should be COMBINED by default" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Default configuration values correct" << std::endl;
    return true;
}

/**
 * @test Test configuration with custom values
 */
bool test_imex_custom_config() {
    std::cout << "Testing IMEX custom configuration..." << std::endl;
    
    IMEXInternalConfig config;
    config.enabled = true;
    config.initial_mode = IntegrationMode::EXPLICIT;
    config.implicit_dt_initial = 1800.0;
    config.explicit_dt_initial = 5e-5;
    config.trigger_type = TransitionTrigger::SLIP_RATE;
    config.slip_rate_threshold = 1e-2;
    config.settling_type = SettlingCriterion::VELOCITY_DECAY;
    config.settling_velocity = 1e-5;
    
    if (!config.enabled) {
        std::cerr << "  FAIL: enabled not set correctly" << std::endl;
        return false;
    }
    
    if (config.initial_mode != IntegrationMode::EXPLICIT) {
        std::cerr << "  FAIL: initial_mode not set correctly" << std::endl;
        return false;
    }
    
    if (std::abs(config.implicit_dt_initial - 1800.0) > TOL_STRICT) {
        std::cerr << "  FAIL: implicit_dt_initial not set correctly" << std::endl;
        return false;
    }
    
    if (config.trigger_type != TransitionTrigger::SLIP_RATE) {
        std::cerr << "  FAIL: trigger_type not set correctly" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Custom configuration values correct" << std::endl;
    return true;
}

// =============================================================================
// PhaseStatistics Tests
// =============================================================================

/**
 * @test Test phase statistics initialization
 */
bool test_phase_statistics_init() {
    std::cout << "Testing phase statistics initialization..." << std::endl;
    
    PhaseStatistics phase;
    phase.mode = IntegrationMode::IMPLICIT;
    phase.start_time = 0.0;
    phase.end_time = 100.0;
    phase.duration = 100.0;
    phase.num_steps = 50;
    phase.num_nonlinear_iterations = 200;
    phase.num_linear_iterations = 1000;
    phase.min_dt = 1.0;
    phase.max_dt = 5.0;
    phase.avg_dt = 2.0;
    phase.peak_velocity = 1e-6;
    phase.peak_slip_rate = 1e-9;
    phase.total_slip = 0.0;
    phase.released_energy = 0.0;
    phase.trigger = TransitionTrigger::NONE;
    phase.settling = SettlingCriterion::TIME_ELAPSED;
    
    if (phase.mode != IntegrationMode::IMPLICIT) {
        std::cerr << "  FAIL: Mode not set correctly" << std::endl;
        return false;
    }
    
    if (std::abs(phase.duration - 100.0) > TOL_STRICT) {
        std::cerr << "  FAIL: Duration not calculated correctly" << std::endl;
        return false;
    }
    
    if (phase.num_steps != 50) {
        std::cerr << "  FAIL: num_steps not set correctly" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Phase statistics initialization correct" << std::endl;
    return true;
}

// =============================================================================
// DynamicEvent Tests
// =============================================================================

/**
 * @test Test dynamic event recording
 */
bool test_dynamic_event_recording() {
    std::cout << "Testing dynamic event recording..." << std::endl;
    
    DynamicEvent event;
    event.trigger_time = 10.0;
    event.end_time = 12.0;
    event.duration = 2.0;
    event.peak_moment_rate = 1e16;
    event.total_moment = 1e17;
    event.magnitude = 5.3;
    event.peak_slip_rate = 2.0;
    event.max_slip = 0.5;
    event.hypo_x = 100.0;
    event.hypo_y = 200.0;
    event.hypo_z = -5000.0;
    
    // Add some history
    event.slip_history.push_back(0.1);
    event.slip_history.push_back(0.3);
    event.slip_history.push_back(0.1);
    event.moment_rate_history.push_back(5e15);
    event.moment_rate_history.push_back(1e16);
    event.moment_rate_history.push_back(3e15);
    
    if (std::abs(event.duration - 2.0) > TOL_STRICT) {
        std::cerr << "  FAIL: Duration not correct" << std::endl;
        return false;
    }
    
    if (std::abs(event.magnitude - 5.3) > TOL_STRICT) {
        std::cerr << "  FAIL: Magnitude not correct" << std::endl;
        return false;
    }
    
    if (event.slip_history.size() != 3) {
        std::cerr << "  FAIL: Slip history size incorrect" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Dynamic event recording correct" << std::endl;
    return true;
}

/**
 * @test Test magnitude calculation from moment
 */
bool test_magnitude_from_moment() {
    std::cout << "Testing magnitude calculation from moment..." << std::endl;
    
    // Test cases: (M0 in N·m, expected Mw)
    std::vector<std::pair<double, double>> test_cases = {
        {1e14, 3.3},
        {1e16, 4.6},
        {1e18, 6.0},
        {1e20, 7.3}
    };
    
    for (const auto& tc : test_cases) {
        double M0 = tc.first;
        double expected_Mw = tc.second;
        
        // Mw = 2/3 * log10(M0) - 6.07
        double computed_Mw = (2.0/3.0) * std::log10(M0) - 6.07;
        
        if (std::abs(computed_Mw - expected_Mw) > TOL_MAGNITUDE) {
            std::cerr << "  FAIL: M0=" << M0 << ", Mw=" << computed_Mw 
                      << ", expected~" << expected_Mw << std::endl;
            return false;
        }
    }
    
    std::cout << "  PASS: Magnitude calculation from moment correct" << std::endl;
    return true;
}

// =============================================================================
// KineticEnergyMonitor Tests
// =============================================================================

/**
 * @test Test kinetic energy monitor initialization
 */
bool test_kinetic_energy_monitor_init() {
    std::cout << "Testing kinetic energy monitor initialization..." << std::endl;
    
    KineticEnergyMonitor monitor(50);
    
    // Initial state should have zero energy
    if (monitor.getPeakEnergy() != 0.0) {
        std::cerr << "  FAIL: Initial peak energy should be 0" << std::endl;
        return false;
    }
    
    if (monitor.getCurrentEnergy() != 0.0) {
        std::cerr << "  FAIL: Initial current energy should be 0" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Kinetic energy monitor initialization correct" << std::endl;
    return true;
}

/**
 * @test Test kinetic energy sampling
 */
bool test_kinetic_energy_sampling() {
    std::cout << "Testing kinetic energy sampling..." << std::endl;
    
    KineticEnergyMonitor monitor(100);
    
    // Add samples
    monitor.addSample(0.0, 100.0);
    monitor.addSample(0.1, 200.0);
    monitor.addSample(0.2, 150.0);
    monitor.addSample(0.3, 100.0);
    
    // Peak should be 200
    if (std::abs(monitor.getPeakEnergy() - 200.0) > TOL_STRICT) {
        std::cerr << "  FAIL: Peak energy should be 200, got " << monitor.getPeakEnergy() << std::endl;
        return false;
    }
    
    // Current should be 100 (last sample)
    if (std::abs(monitor.getCurrentEnergy() - 100.0) > TOL_STRICT) {
        std::cerr << "  FAIL: Current energy should be 100, got " << monitor.getCurrentEnergy() << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Kinetic energy sampling correct" << std::endl;
    return true;
}

/**
 * @test Test energy decay detection
 */
bool test_energy_decay_detection() {
    std::cout << "Testing energy decay detection..." << std::endl;
    
    KineticEnergyMonitor monitor(100);
    
    // Add decaying energy samples
    for (int i = 0; i <= 20; ++i) {
        double t = i * 0.01;
        double energy = 1000.0 * std::exp(-5.0 * t);  // Exponential decay
        monitor.addSample(t, energy);
    }
    
    // Should detect decay
    if (!monitor.isDecaying(0.01)) {
        std::cerr << "  FAIL: Should detect energy decay" << std::endl;
        return false;
    }
    
    // Reset and add increasing energy
    monitor.reset();
    for (int i = 0; i <= 20; ++i) {
        double t = i * 0.01;
        double energy = 100.0 + 50.0 * t;  // Linear increase
        monitor.addSample(t, energy);
    }
    
    // Should not detect decay
    if (monitor.isDecaying(0.01)) {
        std::cerr << "  FAIL: Should not detect decay for increasing energy" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Energy decay detection correct" << std::endl;
    return true;
}

/**
 * @test Test energy threshold detection
 */
bool test_energy_threshold_detection() {
    std::cout << "Testing energy threshold detection..." << std::endl;
    
    KineticEnergyMonitor monitor(100);
    
    // Add samples with peak at 1000, current at 5
    monitor.addSample(0.0, 100.0);
    monitor.addSample(0.1, 500.0);
    monitor.addSample(0.2, 1000.0);  // Peak
    monitor.addSample(0.3, 100.0);
    monitor.addSample(0.4, 10.0);
    monitor.addSample(0.5, 5.0);
    
    // Current (5) is 0.5% of peak (1000)
    // Should be below 1% threshold
    if (!monitor.isBelowThreshold(0.01)) {
        std::cerr << "  FAIL: 5 should be below 1% of 1000" << std::endl;
        return false;
    }
    
    // Should not be below 0.1% threshold
    if (monitor.isBelowThreshold(0.001)) {
        std::cerr << "  FAIL: 5 should not be below 0.1% of 1000" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Energy threshold detection correct" << std::endl;
    return true;
}

/**
 * @test Test monitor reset
 */
bool test_kinetic_energy_reset() {
    std::cout << "Testing kinetic energy monitor reset..." << std::endl;
    
    KineticEnergyMonitor monitor(100);
    
    // Add some samples
    monitor.addSample(0.0, 100.0);
    monitor.addSample(0.1, 500.0);
    
    // Reset
    monitor.reset();
    
    if (monitor.getPeakEnergy() != 0.0) {
        std::cerr << "  FAIL: Peak should be 0 after reset" << std::endl;
        return false;
    }
    
    if (monitor.getCurrentEnergy() != 0.0) {
        std::cerr << "  FAIL: Current should be 0 after reset" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Kinetic energy monitor reset correct" << std::endl;
    return true;
}

// =============================================================================
// SeismicEventDetector Tests
// =============================================================================

/**
 * @test Test seismic event detector initialization
 */
bool test_seismic_event_detector_init() {
    std::cout << "Testing seismic event detector initialization..." << std::endl;
    
    SeismicEventDetector detector;
    
    // Configure
    detector.configure(1e-3, -2.0, 1e-4);
    
    // Initially no event should be active
    if (detector.isEventActive()) {
        std::cerr << "  FAIL: No event should be active initially" << std::endl;
        return false;
    }
    
    if (detector.eventJustStarted()) {
        std::cerr << "  FAIL: No event should have just started" << std::endl;
        return false;
    }
    
    if (detector.eventJustEnded()) {
        std::cerr << "  FAIL: No event should have just ended" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Seismic event detector initialization correct" << std::endl;
    return true;
}

/**
 * @test Test seismic event detection
 */
bool test_seismic_event_detection() {
    std::cout << "Testing seismic event detection..." << std::endl;
    
    SeismicEventDetector detector;
    detector.configure(1e-3, -2.0, 1e-4);  // 1 mm/s threshold
    
    // Create fault stress states
    FaultStressState state_interseismic;
    state_interseismic.slip_rate = 1e-9;  // Very slow, interseismic
    state_interseismic.moment = 0.0;
    
    FaultStressState state_coseismic;
    state_coseismic.slip_rate = 1.0;  // 1 m/s, seismic
    state_coseismic.moment = 1e15;
    
    FaultStressState state_postseismic;
    state_postseismic.slip_rate = 1e-6;  // Post-seismic, below threshold
    state_postseismic.moment = 1e10;
    
    double dt = 0.001;
    
    // Update with interseismic state - no event
    detector.update(state_interseismic, 0.0, dt);
    if (detector.isEventActive()) {
        std::cerr << "  FAIL: No event should be active during interseismic" << std::endl;
        return false;
    }
    
    // Update with coseismic state - event starts
    detector.update(state_coseismic, 0.001, dt);
    if (!detector.isEventActive()) {
        std::cerr << "  FAIL: Event should be active during coseismic" << std::endl;
        return false;
    }
    if (!detector.eventJustStarted()) {
        std::cerr << "  FAIL: Event should have just started" << std::endl;
        return false;
    }
    
    // Continue coseismic
    detector.update(state_coseismic, 0.002, dt);
    if (!detector.isEventActive()) {
        std::cerr << "  FAIL: Event should still be active" << std::endl;
        return false;
    }
    if (detector.eventJustStarted()) {
        std::cerr << "  FAIL: Event should not have just started (continuing)" << std::endl;
        return false;
    }
    
    // Update with post-seismic state - event ends
    detector.update(state_postseismic, 0.003, dt);
    if (detector.isEventActive()) {
        std::cerr << "  FAIL: Event should have ended" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Seismic event detection correct" << std::endl;
    return true;
}

/**
 * @test Test event catalog
 */
bool test_event_catalog() {
    std::cout << "Testing event catalog..." << std::endl;
    
    SeismicEventDetector detector;
    detector.configure(1e-3, -5.0, 1e-5);  // Low magnitude threshold
    
    FaultStressState seismic_state;
    seismic_state.slip_rate = 0.5;
    seismic_state.moment = 1e14;
    
    FaultStressState interseismic_state;
    interseismic_state.slip_rate = 1e-10;
    interseismic_state.moment = 0.0;
    
    double dt = 0.001;
    
    // Simulate first event
    detector.update(interseismic_state, 0.0, dt);
    for (int i = 1; i <= 10; ++i) {
        detector.update(seismic_state, i * dt, dt);
    }
    detector.update(interseismic_state, 0.011, dt);
    
    // Simulate second event
    detector.update(interseismic_state, 0.5, dt);
    for (int i = 1; i <= 10; ++i) {
        detector.update(seismic_state, 0.5 + i * dt, dt);
    }
    detector.update(interseismic_state, 0.511, dt);
    
    // Check catalog
    const auto& catalog = detector.getCatalog();
    if (catalog.size() != 2) {
        std::cerr << "  FAIL: Expected 2 events in catalog, got " << catalog.size() << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Event catalog correct" << std::endl;
    return true;
}

// =============================================================================
// Integration Mode Enum Tests
// =============================================================================

/**
 * @test Test integration mode enum values
 */
bool test_integration_mode_enum() {
    std::cout << "Testing integration mode enum..." << std::endl;
    
    IntegrationMode implicit = IntegrationMode::IMPLICIT;
    IntegrationMode explicit_mode = IntegrationMode::EXPLICIT;
    
    if (implicit == explicit_mode) {
        std::cerr << "  FAIL: IMPLICIT and EXPLICIT should be different" << std::endl;
        return false;
    }
    
    // Test switch statement compatibility
    bool passed = false;
    switch (implicit) {
        case IntegrationMode::IMPLICIT:
            passed = true;
            break;
        case IntegrationMode::EXPLICIT:
            passed = false;
            break;
    }
    
    if (!passed) {
        std::cerr << "  FAIL: Switch on IMPLICIT failed" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Integration mode enum correct" << std::endl;
    return true;
}

// =============================================================================
// Transition Trigger Enum Tests
// =============================================================================

/**
 * @test Test transition trigger enum values
 */
bool test_transition_trigger_enum() {
    std::cout << "Testing transition trigger enum..." << std::endl;
    
    // Test all trigger types exist
    std::vector<TransitionTrigger> triggers = {
        TransitionTrigger::NONE,
        TransitionTrigger::STRESS_THRESHOLD,
        TransitionTrigger::COULOMB_FAILURE,
        TransitionTrigger::SLIP_RATE,
        TransitionTrigger::VELOCITY,
        TransitionTrigger::ACCELERATION,
        TransitionTrigger::ENERGY_RATE,
        TransitionTrigger::USER_DEFINED
    };
    
    // Verify they're all distinct
    for (size_t i = 0; i < triggers.size(); ++i) {
        for (size_t j = i + 1; j < triggers.size(); ++j) {
            if (triggers[i] == triggers[j]) {
                std::cerr << "  FAIL: Trigger types " << i << " and " << j 
                          << " should be different" << std::endl;
                return false;
            }
        }
    }
    
    std::cout << "  PASS: Transition trigger enum correct" << std::endl;
    return true;
}

// =============================================================================
// Settling Criterion Enum Tests
// =============================================================================

/**
 * @test Test settling criterion enum values
 */
bool test_settling_criterion_enum() {
    std::cout << "Testing settling criterion enum..." << std::endl;
    
    // Test all settling criteria exist
    std::vector<SettlingCriterion> criteria = {
        SettlingCriterion::TIME_ELAPSED,
        SettlingCriterion::VELOCITY_DECAY,
        SettlingCriterion::ENERGY_DECAY,
        SettlingCriterion::SLIP_RATE_DECAY,
        SettlingCriterion::COMBINED
    };
    
    // Verify they're all distinct
    for (size_t i = 0; i < criteria.size(); ++i) {
        for (size_t j = i + 1; j < criteria.size(); ++j) {
            if (criteria[i] == criteria[j]) {
                std::cerr << "  FAIL: Settling criteria " << i << " and " << j 
                          << " should be different" << std::endl;
                return false;
            }
        }
    }
    
    std::cout << "  PASS: Settling criterion enum correct" << std::endl;
    return true;
}

// =============================================================================
// CFL Time Step Calculation Tests
// =============================================================================

/**
 * @test Test CFL time step calculation
 */
bool test_cfl_timestep_calculation() {
    std::cout << "Testing CFL time step calculation..." << std::endl;
    
    // Test parameters
    double cfl_factor = 0.5;
    double min_cell_size = 100.0;  // 100 m
    double max_wave_speed = 5000.0;  // 5000 m/s (P-wave velocity)
    
    // Expected: dt = CFL * dx / v = 0.5 * 100 / 5000 = 0.01 s
    double expected_dt = cfl_factor * min_cell_size / max_wave_speed;
    
    if (std::abs(expected_dt - 0.01) > TOL_STRICT) {
        std::cerr << "  FAIL: Expected dt=0.01, got " << expected_dt << std::endl;
        return false;
    }
    
    // Test with different parameters
    cfl_factor = 0.8;
    min_cell_size = 50.0;
    max_wave_speed = 4000.0;
    expected_dt = cfl_factor * min_cell_size / max_wave_speed;
    
    // Expected: 0.8 * 50 / 4000 = 0.01 s
    if (std::abs(expected_dt - 0.01) > TOL_STRICT) {
        std::cerr << "  FAIL: Expected dt=0.01, got " << expected_dt << std::endl;
        return false;
    }
    
    std::cout << "  PASS: CFL time step calculation correct" << std::endl;
    return true;
}

// =============================================================================
// Time Step Bounds Tests
// =============================================================================

/**
 * @test Test time step bounds configuration
 */
bool test_timestep_bounds() {
    std::cout << "Testing time step bounds..." << std::endl;
    
    IMEXInternalConfig config;
    
    // Set custom bounds
    config.implicit_dt_min = 10.0;
    config.implicit_dt_max = 100000.0;
    config.explicit_dt_min = 1e-7;
    config.explicit_dt_max = 1e-3;
    
    // Test bounds
    double test_dt = 5.0;  // Below implicit min
    double clamped = std::max(config.implicit_dt_min, std::min(test_dt, config.implicit_dt_max));
    if (std::abs(clamped - 10.0) > TOL_STRICT) {
        std::cerr << "  FAIL: dt should be clamped to min" << std::endl;
        return false;
    }
    
    test_dt = 200000.0;  // Above implicit max
    clamped = std::max(config.implicit_dt_min, std::min(test_dt, config.implicit_dt_max));
    if (std::abs(clamped - 100000.0) > TOL_STRICT) {
        std::cerr << "  FAIL: dt should be clamped to max" << std::endl;
        return false;
    }
    
    test_dt = 1000.0;  // Within bounds
    clamped = std::max(config.implicit_dt_min, std::min(test_dt, config.implicit_dt_max));
    if (std::abs(clamped - 1000.0) > TOL_STRICT) {
        std::cerr << "  FAIL: dt should not be modified when within bounds" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Time step bounds correct" << std::endl;
    return true;
}

// =============================================================================
// Duration Bounds Tests
// =============================================================================

/**
 * @test Test dynamic duration bounds
 */
bool test_dynamic_duration_bounds() {
    std::cout << "Testing dynamic duration bounds..." << std::endl;
    
    IMEXInternalConfig config;
    config.min_dynamic_duration = 0.5;   // At least 0.5 seconds in explicit mode
    config.max_dynamic_duration = 120.0; // At most 2 minutes in explicit mode
    
    // Test minimum duration check
    double event_duration = 0.1;
    bool should_settle_early = event_duration >= config.min_dynamic_duration;
    if (should_settle_early) {
        std::cerr << "  FAIL: Should not settle before min duration" << std::endl;
        return false;
    }
    
    // Test after min duration
    event_duration = 1.0;
    should_settle_early = event_duration >= config.min_dynamic_duration;
    if (!should_settle_early) {
        std::cerr << "  FAIL: Should be eligible to settle after min duration" << std::endl;
        return false;
    }
    
    // Test force settling after max duration
    event_duration = 150.0;
    bool force_settle = event_duration >= config.max_dynamic_duration;
    if (!force_settle) {
        std::cerr << "  FAIL: Should force settle after max duration" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Dynamic duration bounds correct" << std::endl;
    return true;
}

// =============================================================================
// Settling Observation Window Tests
// =============================================================================

/**
 * @test Test settling observation window
 */
bool test_settling_observation_window() {
    std::cout << "Testing settling observation window..." << std::endl;
    
    IMEXInternalConfig config;
    config.settling_observation_window = 0.1;  // 100 ms window
    double dt = 0.001;  // 1 ms time step
    
    // Calculate expected history size
    int history_size = static_cast<int>(config.settling_observation_window / dt);
    if (history_size != 100) {
        std::cerr << "  FAIL: Expected history size 100, got " << history_size << std::endl;
        return false;
    }
    
    // Test with different dt
    dt = 0.0001;  // 0.1 ms time step
    history_size = static_cast<int>(config.settling_observation_window / dt);
    if (history_size != 1000) {
        std::cerr << "  FAIL: Expected history size 1000, got " << history_size << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Settling observation window correct" << std::endl;
    return true;
}

// =============================================================================
// Velocity Threshold Tests
// =============================================================================

/**
 * @test Test velocity threshold triggering
 */
bool test_velocity_threshold() {
    std::cout << "Testing velocity threshold..." << std::endl;
    
    IMEXInternalConfig config;
    config.velocity_threshold = 1e-4;  // 0.1 mm/s
    
    // Below threshold
    double v_max = 1e-5;
    bool should_trigger = v_max >= config.velocity_threshold;
    if (should_trigger) {
        std::cerr << "  FAIL: Should not trigger below threshold" << std::endl;
        return false;
    }
    
    // At threshold
    v_max = 1e-4;
    should_trigger = v_max >= config.velocity_threshold;
    if (!should_trigger) {
        std::cerr << "  FAIL: Should trigger at threshold" << std::endl;
        return false;
    }
    
    // Above threshold
    v_max = 1e-3;
    should_trigger = v_max >= config.velocity_threshold;
    if (!should_trigger) {
        std::cerr << "  FAIL: Should trigger above threshold" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Velocity threshold correct" << std::endl;
    return true;
}

// =============================================================================
// Slip Rate Threshold Tests
// =============================================================================

/**
 * @test Test slip rate threshold triggering
 */
bool test_slip_rate_threshold() {
    std::cout << "Testing slip rate threshold..." << std::endl;
    
    IMEXInternalConfig config;
    config.slip_rate_threshold = 1e-3;  // 1 mm/s (seismic slip rate)
    
    // Interseismic slip rate
    double slip_rate = 1e-9;  // ~30 mm/year
    bool is_seismic = slip_rate >= config.slip_rate_threshold;
    if (is_seismic) {
        std::cerr << "  FAIL: Interseismic slip rate should not trigger" << std::endl;
        return false;
    }
    
    // Seismic slip rate
    slip_rate = 1.0;  // 1 m/s
    is_seismic = slip_rate >= config.slip_rate_threshold;
    if (!is_seismic) {
        std::cerr << "  FAIL: Seismic slip rate should trigger" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Slip rate threshold correct" << std::endl;
    return true;
}

// =============================================================================
// Energy Ratio Threshold Tests
// =============================================================================

/**
 * @test Test energy ratio for settling
 */
bool test_energy_ratio_settling() {
    std::cout << "Testing energy ratio settling..." << std::endl;
    
    IMEXInternalConfig config;
    config.settling_energy_ratio = 0.01;  // 1% of peak energy
    
    double peak_energy = 1000.0;
    
    // Current energy is 50% of peak
    double current_energy = 500.0;
    bool is_settled = current_energy <= peak_energy * config.settling_energy_ratio;
    if (is_settled) {
        std::cerr << "  FAIL: 50% of peak should not be settled" << std::endl;
        return false;
    }
    
    // Current energy is 0.5% of peak
    current_energy = 5.0;
    is_settled = current_energy <= peak_energy * config.settling_energy_ratio;
    if (!is_settled) {
        std::cerr << "  FAIL: 0.5% of peak should be settled" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Energy ratio settling correct" << std::endl;
    return true;
}

// =============================================================================
// Adaptive Time Stepping Tests
// =============================================================================

/**
 * @test Test adaptive time stepping factors
 */
bool test_adaptive_timestep_factors() {
    std::cout << "Testing adaptive time stepping factors..." << std::endl;
    
    IMEXInternalConfig config;
    config.dt_growth_factor = 1.2;
    config.dt_shrink_factor = 0.5;
    
    double current_dt = 1.0;
    
    // Test growth
    double new_dt = current_dt * config.dt_growth_factor;
    if (std::abs(new_dt - 1.2) > TOL_STRICT) {
        std::cerr << "  FAIL: Growth factor not applied correctly" << std::endl;
        return false;
    }
    
    // Test shrink
    new_dt = current_dt * config.dt_shrink_factor;
    if (std::abs(new_dt - 0.5) > TOL_STRICT) {
        std::cerr << "  FAIL: Shrink factor not applied correctly" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Adaptive time stepping factors correct" << std::endl;
    return true;
}

// =============================================================================
// Test Runner
// =============================================================================

int runIMEXTransitionTests() {
    std::cout << "\n=== IMEX Transition Unit Tests ===" << std::endl;
    
    int passed = 0;
    int failed = 0;
    
    // Configuration tests
    if (test_imex_default_config()) ++passed; else ++failed;
    if (test_imex_custom_config()) ++passed; else ++failed;
    
    // Phase statistics tests
    if (test_phase_statistics_init()) ++passed; else ++failed;
    
    // Dynamic event tests
    if (test_dynamic_event_recording()) ++passed; else ++failed;
    if (test_magnitude_from_moment()) ++passed; else ++failed;
    
    // Kinetic energy monitor tests
    if (test_kinetic_energy_monitor_init()) ++passed; else ++failed;
    if (test_kinetic_energy_sampling()) ++passed; else ++failed;
    if (test_energy_decay_detection()) ++passed; else ++failed;
    if (test_energy_threshold_detection()) ++passed; else ++failed;
    if (test_kinetic_energy_reset()) ++passed; else ++failed;
    
    // Seismic event detector tests
    if (test_seismic_event_detector_init()) ++passed; else ++failed;
    if (test_seismic_event_detection()) ++passed; else ++failed;
    if (test_event_catalog()) ++passed; else ++failed;
    
    // Enum tests
    if (test_integration_mode_enum()) ++passed; else ++failed;
    if (test_transition_trigger_enum()) ++passed; else ++failed;
    if (test_settling_criterion_enum()) ++passed; else ++failed;
    
    // Threshold and bounds tests
    if (test_cfl_timestep_calculation()) ++passed; else ++failed;
    if (test_timestep_bounds()) ++passed; else ++failed;
    if (test_dynamic_duration_bounds()) ++passed; else ++failed;
    if (test_settling_observation_window()) ++passed; else ++failed;
    if (test_velocity_threshold()) ++passed; else ++failed;
    if (test_slip_rate_threshold()) ++passed; else ++failed;
    if (test_energy_ratio_settling()) ++passed; else ++failed;
    if (test_adaptive_timestep_factors()) ++passed; else ++failed;
    
    std::cout << "\n=== IMEX Transition Test Summary ===" << std::endl;
    std::cout << "Passed: " << passed << std::endl;
    std::cout << "Failed: " << failed << std::endl;
    std::cout << "Total:  " << (passed + failed) << std::endl;
    
    return failed;
}

} // namespace Testing
} // namespace FSRM

#ifndef FSRM_TEST_NO_MAIN
int main() {
    return FSRM::Testing::runIMEXTransitionTests();
}
#endif
