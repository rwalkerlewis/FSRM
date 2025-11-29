/*
 * Implicit-Explicit Time Integration Transition Manager
 * 
 * Implementation of automatic transitions between implicit (quasi-static)
 * and explicit (dynamic) time integration for induced seismicity modeling.
 * 
 * Physical Motivation:
 * - Pore pressure diffusion and slow fault loading: IMPLICIT (hours-days)
 * - Dynamic rupture and wave propagation: EXPLICIT (milliseconds)
 * - The transition captures the full physics of induced seismicity
 */

#include "ImplicitExplicitTransition.hpp"
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace FSRM {

// ============================================================================
// ImplicitExplicitTransitionManager Implementation
// ============================================================================

ImplicitExplicitTransitionManager::ImplicitExplicitTransitionManager(MPI_Comm comm_in)
    : comm(comm_in), ts(nullptr), dm(nullptr), solution(nullptr),
      velocity_old(nullptr), acceleration(nullptr), mass_matrix(nullptr),
      mass_lumped(nullptr), current_mode(IntegrationMode::IMPLICIT),
      current_dt(1.0), current_time(0.0), current_step(0),
      in_dynamic_event(false), event_start_time(0.0), peak_kinetic_energy(0.0),
      fault_model(nullptr), fault_network(nullptr) {
    
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
}

ImplicitExplicitTransitionManager::~ImplicitExplicitTransitionManager() {
    if (velocity_old) VecDestroy(&velocity_old);
    if (acceleration) VecDestroy(&acceleration);
    if (mass_lumped) MatDestroy(&mass_lumped);
}

PetscErrorCode ImplicitExplicitTransitionManager::initialize(const IMEXInternalConfig& config_in) {
    PetscFunctionBeginUser;
    
    config = config_in;
    current_mode = config.initial_mode;
    current_dt = (current_mode == IntegrationMode::IMPLICIT) 
                 ? config.implicit_dt_initial 
                 : config.explicit_dt_initial;
    
    // Initialize phase tracking
    current_phase.mode = current_mode;
    current_phase.start_time = 0.0;
    current_phase.num_steps = 0;
    current_phase.num_nonlinear_iterations = 0;
    current_phase.num_linear_iterations = 0;
    current_phase.min_dt = current_dt;
    current_phase.max_dt = current_dt;
    current_phase.peak_velocity = 0.0;
    current_phase.peak_slip_rate = 0.0;
    current_phase.total_slip = 0.0;
    current_phase.released_energy = 0.0;
    
    if (rank == 0 && config.log_transitions) {
        PetscPrintf(comm, "IMEX Manager: Initialized in %s mode\n",
                   current_mode == IntegrationMode::IMPLICIT ? "IMPLICIT" : "EXPLICIT");
        PetscPrintf(comm, "  Trigger type: %d, Threshold: %g\n",
                   static_cast<int>(config.trigger_type), config.stress_threshold);
        PetscPrintf(comm, "  Settling type: %d, Duration bounds: [%g, %g] s\n",
                   static_cast<int>(config.settling_type),
                   config.min_dynamic_duration, config.max_dynamic_duration);
    }
    
    PetscFunctionReturn(0);
}

PetscErrorCode ImplicitExplicitTransitionManager::setup(TS ts_in, DM dm_in, 
                                                        Vec solution_in, Mat mass_in) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    ts = ts_in;
    dm = dm_in;
    solution = solution_in;
    mass_matrix = mass_in;
    
    // Create work vectors
    ierr = VecDuplicate(solution, &velocity_old); CHKERRQ(ierr);
    ierr = VecDuplicate(solution, &acceleration); CHKERRQ(ierr);
    ierr = VecSet(velocity_old, 0.0); CHKERRQ(ierr);
    ierr = VecSet(acceleration, 0.0); CHKERRQ(ierr);
    
    // Configure initial mode
    if (current_mode == IntegrationMode::IMPLICIT) {
        ierr = configureImplicitMode(); CHKERRQ(ierr);
    } else {
        ierr = configureExplicitMode(); CHKERRQ(ierr);
    }
    
    // Register monitor for automatic transition checking
    ierr = TSMonitorSet(ts, IMEXMonitor, this, nullptr); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

void ImplicitExplicitTransitionManager::setFaultModel(SeismicFaultModel* fault) {
    fault_model = fault;
}

void ImplicitExplicitTransitionManager::setFaultNetwork(FaultNetwork* network) {
    fault_network = network;
}

void ImplicitExplicitTransitionManager::setCustomTrigger(std::function<bool(Vec, double)> func) {
    custom_trigger = func;
}

void ImplicitExplicitTransitionManager::setCustomSettling(std::function<bool(Vec, double, double)> func) {
    custom_settling = func;
}

PetscErrorCode ImplicitExplicitTransitionManager::configureImplicitMode() {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    // Set time stepping method for implicit integration
    if (config.implicit_method == "BEULER") {
        ierr = TSSetType(ts, TSBEULER); CHKERRQ(ierr);
    } else if (config.implicit_method == "BDF") {
        ierr = TSSetType(ts, TSBDF); CHKERRQ(ierr);
    } else if (config.implicit_method == "ALPHA2") {
        ierr = TSSetType(ts, TSALPHA2); CHKERRQ(ierr);
    } else {
        // Default to backward Euler
        ierr = TSSetType(ts, TSBEULER); CHKERRQ(ierr);
    }
    
    // Set time step bounds
    ierr = TSSetTimeStep(ts, config.implicit_dt_initial); CHKERRQ(ierr);
    
    // Get SNES for solver settings
    SNES snes;
    ierr = TSGetSNES(ts, &snes); CHKERRQ(ierr);
    ierr = SNESSetTolerances(snes, config.implicit_atol, config.implicit_rtol,
                            PETSC_DEFAULT, config.implicit_max_iterations,
                            PETSC_DEFAULT); CHKERRQ(ierr);
    
    current_dt = config.implicit_dt_initial;
    
    if (rank == 0) {
        PetscPrintf(comm, "IMEX: Configured IMPLICIT mode (dt=%g s, method=%s)\n",
                   current_dt, config.implicit_method.c_str());
    }
    
    PetscFunctionReturn(0);
}

PetscErrorCode ImplicitExplicitTransitionManager::configureExplicitMode() {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    // Set time stepping method for explicit integration
    if (config.explicit_method == "EULER") {
        ierr = TSSetType(ts, TSEULER); CHKERRQ(ierr);
    } else if (config.explicit_method == "RK") {
        ierr = TSSetType(ts, TSRK); CHKERRQ(ierr);
        ierr = TSRKSetType(ts, TSRK4); CHKERRQ(ierr);
    } else if (config.explicit_method == "NEWMARK") {
        // Generalized-alpha with explicit parameters
        ierr = TSSetType(ts, TSALPHA2); CHKERRQ(ierr);
        // Set spectral radius to 1 for explicit Newmark
    } else {
        // Default to forward Euler
        ierr = TSSetType(ts, TSEULER); CHKERRQ(ierr);
    }
    
    // Set time step (should be CFL-limited)
    ierr = TSSetTimeStep(ts, config.explicit_dt_initial); CHKERRQ(ierr);
    
    // For explicit methods, we typically don't need SNES iterations
    // But keep settings reasonable in case of implicit parts
    SNES snes;
    ierr = TSGetSNES(ts, &snes); CHKERRQ(ierr);
    ierr = SNESSetTolerances(snes, config.explicit_atol, config.explicit_rtol,
                            PETSC_DEFAULT, config.explicit_max_iterations,
                            PETSC_DEFAULT); CHKERRQ(ierr);
    
    current_dt = config.explicit_dt_initial;
    
    if (rank == 0) {
        PetscPrintf(comm, "IMEX: Configured EXPLICIT mode (dt=%g s, method=%s)\n",
                   current_dt, config.explicit_method.c_str());
    }
    
    PetscFunctionReturn(0);
}

PetscErrorCode ImplicitExplicitTransitionManager::checkAndTransition(
    Vec solution, Vec velocity, double time) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    if (!config.enabled) {
        PetscFunctionReturn(0);
    }
    
    current_time = time;
    current_step++;
    
    // Update phase statistics
    current_phase.num_steps++;
    double v_max = computeMaxVelocity(velocity);
    current_phase.peak_velocity = std::max(current_phase.peak_velocity, v_max);
    current_phase.min_dt = std::min(current_phase.min_dt, current_dt);
    current_phase.max_dt = std::max(current_phase.max_dt, current_dt);
    
    // Track slip rate if fault model available
    if (fault_model || fault_network) {
        double slip_rate = computeMaxSlipRate();
        current_phase.peak_slip_rate = std::max(current_phase.peak_slip_rate, slip_rate);
    }
    
    // Check for mode transition
    if (current_mode == IntegrationMode::IMPLICIT) {
        // Check if we should transition to EXPLICIT
        if (checkTriggerCondition(solution, velocity, time)) {
            // Log transition
            ierr = logTransition(IntegrationMode::IMPLICIT, IntegrationMode::EXPLICIT,
                               time, "Trigger condition met"); CHKERRQ(ierr);
            
            // End current phase
            endCurrentPhase(time, config.trigger_type, SettlingCriterion::TIME_ELAPSED);
            
            // Start dynamic event tracking
            // TODO: Extract hypocenter from fault model
            startDynamicEvent(time, 0.0, 0.0, 0.0);
            
            // Configure explicit mode
            ierr = configureExplicitMode(); CHKERRQ(ierr);
            current_mode = IntegrationMode::EXPLICIT;
            in_dynamic_event = true;
            event_start_time = time;
            peak_kinetic_energy = computeKineticEnergy(velocity);
            
            // Start new phase
            startNewPhase(IntegrationMode::EXPLICIT, time);
            
            // Optionally smooth transition with ramped damping
            if (config.smooth_transition) {
                ierr = adjustDamping(0.0); CHKERRQ(ierr);  // Start with zero extra damping
            }
        }
    } else {
        // EXPLICIT mode - check if we should return to IMPLICIT
        double event_duration = time - event_start_time;
        
        // Update dynamic event
        double slip_rate = computeMaxSlipRate();
        double kinetic_energy = computeKineticEnergy(velocity);
        double moment_rate = 0.0;  // TODO: compute from fault model
        updateDynamicEvent(slip_rate, moment_rate);
        
        // Track peak kinetic energy
        peak_kinetic_energy = std::max(peak_kinetic_energy, kinetic_energy);
        
        // Add to history buffers for settling detection
        velocity_history.push_back(v_max);
        energy_history.push_back(kinetic_energy);
        slip_rate_history.push_back(slip_rate);
        
        // Keep history buffers bounded
        int history_size = static_cast<int>(config.settling_observation_window / current_dt);
        while (velocity_history.size() > static_cast<size_t>(history_size)) {
            velocity_history.pop_front();
            energy_history.pop_front();
            slip_rate_history.pop_front();
        }
        
        // Check settling condition (with minimum duration)
        bool should_settle = false;
        if (event_duration >= config.min_dynamic_duration) {
            should_settle = checkSettlingCondition(solution, velocity, time);
        }
        
        // Force settling after max duration
        if (event_duration >= config.max_dynamic_duration) {
            should_settle = true;
            if (rank == 0) {
                PetscPrintf(comm, "IMEX: Forcing return to IMPLICIT after max duration %g s\n",
                           config.max_dynamic_duration);
            }
        }
        
        if (should_settle) {
            // Log transition
            ierr = logTransition(IntegrationMode::EXPLICIT, IntegrationMode::IMPLICIT,
                               time, "Settling condition met"); CHKERRQ(ierr);
            
            // End dynamic event
            endDynamicEvent(time);
            
            // End current phase
            endCurrentPhase(time, TransitionTrigger::NONE, config.settling_type);
            
            // Configure implicit mode
            ierr = configureImplicitMode(); CHKERRQ(ierr);
            current_mode = IntegrationMode::IMPLICIT;
            in_dynamic_event = false;
            
            // Start new phase
            startNewPhase(IntegrationMode::IMPLICIT, time);
            
            // Clear history buffers
            velocity_history.clear();
            energy_history.clear();
            slip_rate_history.clear();
            
            // Ramp damping if smooth transition enabled
            if (config.smooth_transition) {
                ierr = adjustDamping(1.0); CHKERRQ(ierr);
            }
        }
    }
    
    // Store velocity for next step
    ierr = VecCopy(velocity, velocity_old); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

bool ImplicitExplicitTransitionManager::checkTriggerCondition(
    Vec solution, Vec velocity, double time) {
    
    switch (config.trigger_type) {
        case TransitionTrigger::STRESS_THRESHOLD: {
            double vm_stress = computeVonMisesStress(solution);
            return vm_stress >= config.stress_threshold;
        }
        
        case TransitionTrigger::COULOMB_FAILURE: {
            double cff = computeCoulombFailure(solution);
            return cff >= config.coulomb_threshold;
        }
        
        case TransitionTrigger::SLIP_RATE: {
            double slip_rate = computeMaxSlipRate();
            return slip_rate >= config.slip_rate_threshold;
        }
        
        case TransitionTrigger::VELOCITY: {
            double v_max = computeMaxVelocity(velocity);
            return v_max >= config.velocity_threshold;
        }
        
        case TransitionTrigger::ACCELERATION: {
            double a_max = computeMaxAcceleration(velocity);
            return a_max >= config.acceleration_threshold;
        }
        
        case TransitionTrigger::ENERGY_RATE: {
            // Compute rate of change of kinetic energy
            double ke_current = computeKineticEnergy(velocity);
            static double ke_prev = 0.0;
            static double time_prev = 0.0;
            double rate = (time > time_prev) ? (ke_current - ke_prev) / (time - time_prev) : 0.0;
            ke_prev = ke_current;
            time_prev = time;
            return rate >= config.energy_rate_threshold;
        }
        
        case TransitionTrigger::USER_DEFINED: {
            if (custom_trigger) {
                return custom_trigger(solution, time);
            }
            return false;
        }
        
        case TransitionTrigger::NONE:
        default:
            return false;
    }
}

bool ImplicitExplicitTransitionManager::checkSettlingCondition(
    Vec solution, Vec velocity, double time) {
    
    double event_duration = time - event_start_time;
    
    switch (config.settling_type) {
        case SettlingCriterion::TIME_ELAPSED:
            return event_duration >= config.min_dynamic_duration;
        
        case SettlingCriterion::VELOCITY_DECAY: {
            double v_max = computeMaxVelocity(velocity);
            return v_max <= config.settling_velocity;
        }
        
        case SettlingCriterion::ENERGY_DECAY: {
            double ke_current = computeKineticEnergy(velocity);
            return ke_current <= peak_kinetic_energy * config.settling_energy_ratio;
        }
        
        case SettlingCriterion::SLIP_RATE_DECAY: {
            double slip_rate = computeMaxSlipRate();
            return slip_rate <= config.settling_slip_rate;
        }
        
        case SettlingCriterion::COMBINED: {
            // All criteria must be met
            double v_max = computeMaxVelocity(velocity);
            double ke_current = computeKineticEnergy(velocity);
            double slip_rate = computeMaxSlipRate();
            
            bool velocity_settled = v_max <= config.settling_velocity;
            bool energy_settled = ke_current <= peak_kinetic_energy * config.settling_energy_ratio;
            bool slip_settled = slip_rate <= config.settling_slip_rate;
            
            return velocity_settled && energy_settled && slip_settled;
        }
        
        default:
            return false;
    }
}

PetscErrorCode ImplicitExplicitTransitionManager::forceTransition(IntegrationMode target_mode) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    if (target_mode == current_mode) {
        PetscFunctionReturn(0);
    }
    
    ierr = logTransition(current_mode, target_mode, current_time, "Forced transition"); 
    CHKERRQ(ierr);
    
    if (target_mode == IntegrationMode::EXPLICIT) {
        ierr = configureExplicitMode(); CHKERRQ(ierr);
        in_dynamic_event = true;
        event_start_time = current_time;
    } else {
        ierr = configureImplicitMode(); CHKERRQ(ierr);
        in_dynamic_event = false;
    }
    
    current_mode = target_mode;
    
    PetscFunctionReturn(0);
}

double ImplicitExplicitTransitionManager::computeMaxVelocity(Vec velocity) {
    PetscReal v_max;
    VecNorm(velocity, NORM_INFINITY, &v_max);
    return static_cast<double>(v_max);
}

double ImplicitExplicitTransitionManager::computeMaxAcceleration(Vec velocity) {
    // Compute acceleration as (v_new - v_old) / dt
    PetscErrorCode ierr;
    Vec diff;
    VecDuplicate(velocity, &diff);
    VecCopy(velocity, diff);
    VecAXPY(diff, -1.0, velocity_old);  // diff = v - v_old
    if (current_dt > 0.0) {
        VecScale(diff, 1.0 / current_dt);
    }
    
    PetscReal a_max;
    VecNorm(diff, NORM_INFINITY, &a_max);
    VecDestroy(&diff);
    
    return static_cast<double>(a_max);
}

double ImplicitExplicitTransitionManager::computeKineticEnergy(Vec velocity) {
    // KE = 0.5 * v^T * M * v
    // For simplicity, assume uniform mass and compute 0.5 * ||v||^2
    // TODO: Use actual mass matrix
    PetscReal v_norm;
    VecNorm(velocity, NORM_2, &v_norm);
    return 0.5 * v_norm * v_norm;
}

double ImplicitExplicitTransitionManager::computeMaxSlipRate() {
    if (fault_model) {
        // Get slip rate from fault model
        // TODO: Access actual slip rate from fault state
        return 0.0;
    }
    
    if (fault_network) {
        // Get maximum slip rate across all faults
        double max_rate = 0.0;
        for (size_t i = 0; i < fault_network->numFaults(); ++i) {
            auto* fault = fault_network->getFault(i);
            // TODO: Access actual slip rate
        }
        return max_rate;
    }
    
    return 0.0;
}

double ImplicitExplicitTransitionManager::computeCoulombFailure(Vec solution) {
    // Compute Coulomb Failure Function from stress state
    // CFF = τ - μ(σ_n - p) - c
    // TODO: Extract stress from solution and compute CFF
    
    if (fault_model) {
        // Get CFF from fault model
        // TODO: Implement
    }
    
    return -1e10;  // Return large negative value (stable)
}

double ImplicitExplicitTransitionManager::computeVonMisesStress(Vec solution) {
    // Compute von Mises stress from solution
    // σ_vm = sqrt(0.5 * ((σ1-σ2)² + (σ2-σ3)² + (σ3-σ1)²))
    // TODO: Extract stress from solution
    return 0.0;
}

double ImplicitExplicitTransitionManager::computeCFLTimeStep(
    double min_cell_size, double max_wave_speed) {
    
    return config.cfl_factor * min_cell_size / max_wave_speed;
}

double ImplicitExplicitTransitionManager::getRecommendedDt(Vec solution, Vec velocity) {
    if (current_mode == IntegrationMode::IMPLICIT) {
        // For implicit, use error-based adaptation
        // TODO: Compute error estimate and adapt
        return current_dt;
    } else {
        // For explicit, use CFL condition
        // TODO: Get actual cell size and wave speed from DM
        double min_cell_size = 10.0;  // meters (placeholder)
        double max_wave_speed = 5000.0;  // m/s (placeholder)
        double cfl_dt = computeCFLTimeStep(min_cell_size, max_wave_speed);
        return std::min(cfl_dt, config.explicit_dt_max);
    }
}

PetscErrorCode ImplicitExplicitTransitionManager::interpolateSolution(
    double old_dt, double new_dt) {
    PetscFunctionBeginUser;
    // TODO: Implement solution interpolation for smooth transition
    PetscFunctionReturn(0);
}

PetscErrorCode ImplicitExplicitTransitionManager::adjustDamping(double ramp_factor) {
    PetscFunctionBeginUser;
    // TODO: Adjust Rayleigh damping parameters in physics kernels
    PetscFunctionReturn(0);
}

void ImplicitExplicitTransitionManager::startNewPhase(IntegrationMode mode, double time) {
    current_phase = PhaseStatistics();
    current_phase.mode = mode;
    current_phase.start_time = time;
    current_phase.num_steps = 0;
    current_phase.num_nonlinear_iterations = 0;
    current_phase.num_linear_iterations = 0;
    current_phase.min_dt = (mode == IntegrationMode::IMPLICIT) 
                          ? config.implicit_dt_max : config.explicit_dt_max;
    current_phase.max_dt = 0.0;
    current_phase.peak_velocity = 0.0;
    current_phase.peak_slip_rate = 0.0;
    current_phase.total_slip = 0.0;
    current_phase.released_energy = 0.0;
}

void ImplicitExplicitTransitionManager::endCurrentPhase(
    double time, TransitionTrigger trigger, SettlingCriterion settling) {
    
    current_phase.end_time = time;
    current_phase.duration = time - current_phase.start_time;
    current_phase.avg_dt = current_phase.duration / std::max(1, current_phase.num_steps);
    current_phase.trigger = trigger;
    current_phase.settling = settling;
    
    phase_history.push_back(current_phase);
}

void ImplicitExplicitTransitionManager::startDynamicEvent(
    double time, double x, double y, double z) {
    
    current_event = DynamicEvent();
    current_event.trigger_time = time;
    current_event.hypo_x = x;
    current_event.hypo_y = y;
    current_event.hypo_z = z;
    current_event.peak_moment_rate = 0.0;
    current_event.total_moment = 0.0;
    current_event.peak_slip_rate = 0.0;
    current_event.max_slip = 0.0;
}

void ImplicitExplicitTransitionManager::updateDynamicEvent(
    double slip_rate, double moment_rate) {
    
    current_event.peak_slip_rate = std::max(current_event.peak_slip_rate, slip_rate);
    current_event.peak_moment_rate = std::max(current_event.peak_moment_rate, moment_rate);
    current_event.total_moment += moment_rate * current_dt;
    
    current_event.slip_history.push_back(slip_rate * current_dt);
    current_event.moment_rate_history.push_back(moment_rate);
}

void ImplicitExplicitTransitionManager::endDynamicEvent(double time) {
    current_event.end_time = time;
    current_event.duration = time - current_event.trigger_time;
    
    // Compute cumulative slip
    current_event.max_slip = 0.0;
    for (double slip : current_event.slip_history) {
        current_event.max_slip += slip;
    }
    
    // Compute moment magnitude: Mw = 2/3 * log10(M0) - 6.07
    if (current_event.total_moment > 0.0) {
        current_event.magnitude = (2.0/3.0) * std::log10(current_event.total_moment) - 6.07;
    } else {
        current_event.magnitude = -99.0;  // No significant moment
    }
    
    dynamic_events.push_back(current_event);
    
    if (rank == 0 && config.log_transitions) {
        PetscPrintf(comm, "IMEX: Dynamic event ended (duration=%g s, Mw=%g, slip=%g m)\n",
                   current_event.duration, current_event.magnitude, current_event.max_slip);
    }
}

PetscErrorCode ImplicitExplicitTransitionManager::logTransition(
    IntegrationMode from, IntegrationMode to, double time, const std::string& reason) {
    PetscFunctionBeginUser;
    
    if (rank == 0) {
        const char* from_str = (from == IntegrationMode::IMPLICIT) ? "IMPLICIT" : "EXPLICIT";
        const char* to_str = (to == IntegrationMode::IMPLICIT) ? "IMPLICIT" : "EXPLICIT";
        
        PetscPrintf(comm, "\n============================================================\n");
        PetscPrintf(comm, "IMEX TRANSITION: %s -> %s at t = %g s\n", from_str, to_str, time);
        PetscPrintf(comm, "Reason: %s\n", reason.c_str());
        PetscPrintf(comm, "============================================================\n\n");
        
        // Write to log file if enabled
        if (config.log_transitions && !config.transition_log_file.empty()) {
            std::ofstream log(config.transition_log_file, std::ios::app);
            if (log.is_open()) {
                log << std::setprecision(12);
                log << time << " " << from_str << " -> " << to_str << " : " << reason << "\n";
                log.close();
            }
        }
    }
    
    PetscFunctionReturn(0);
}

PetscErrorCode ImplicitExplicitTransitionManager::writeTransitionLog() {
    PetscFunctionBeginUser;
    
    if (rank != 0) {
        PetscFunctionReturn(0);
    }
    
    std::ofstream log(config.transition_log_file);
    if (!log.is_open()) {
        PetscFunctionReturn(0);
    }
    
    log << "# IMEX Transition Log\n";
    log << "# Time (s), From Mode, To Mode, Trigger, Settling\n";
    log << "\n";
    
    log << "# Phase Summary:\n";
    for (size_t i = 0; i < phase_history.size(); ++i) {
        const auto& phase = phase_history[i];
        log << "Phase " << i << ": " 
            << (phase.mode == IntegrationMode::IMPLICIT ? "IMPLICIT" : "EXPLICIT")
            << " t=[" << phase.start_time << ", " << phase.end_time << "] "
            << "steps=" << phase.num_steps << " "
            << "dt=[" << phase.min_dt << ", " << phase.max_dt << "]\n";
    }
    
    log << "\n# Dynamic Events:\n";
    for (size_t i = 0; i < dynamic_events.size(); ++i) {
        const auto& event = dynamic_events[i];
        log << "Event " << i << ": "
            << "t=" << event.trigger_time << " "
            << "duration=" << event.duration << " "
            << "Mw=" << event.magnitude << " "
            << "slip=" << event.max_slip << "\n";
    }
    
    log.close();
    
    PetscFunctionReturn(0);
}

void ImplicitExplicitTransitionManager::getSummaryStatistics(
    int& num_implicit_phases, int& num_explicit_phases,
    double& total_implicit_time, double& total_explicit_time,
    int& num_events, double& total_moment) {
    
    num_implicit_phases = 0;
    num_explicit_phases = 0;
    total_implicit_time = 0.0;
    total_explicit_time = 0.0;
    
    for (const auto& phase : phase_history) {
        if (phase.mode == IntegrationMode::IMPLICIT) {
            num_implicit_phases++;
            total_implicit_time += phase.duration;
        } else {
            num_explicit_phases++;
            total_explicit_time += phase.duration;
        }
    }
    
    num_events = static_cast<int>(dynamic_events.size());
    total_moment = 0.0;
    for (const auto& event : dynamic_events) {
        total_moment += event.total_moment;
    }
}

// ============================================================================
// PETSc Callback Functions
// ============================================================================

PetscErrorCode IMEXMonitor(TS ts, PetscInt step, PetscReal time, Vec U, void* ctx) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    auto* manager = static_cast<ImplicitExplicitTransitionManager*>(ctx);
    
    // Get velocity from TS (for second-order systems)
    Vec velocity;
    ierr = VecDuplicate(U, &velocity); CHKERRQ(ierr);
    
    // Try to get velocity from TS context or compute from solution
    // For now, use a simple finite difference approximation
    static Vec U_old = nullptr;
    static double time_old = 0.0;
    
    if (U_old == nullptr) {
        ierr = VecDuplicate(U, &U_old); CHKERRQ(ierr);
        ierr = VecCopy(U, U_old); CHKERRQ(ierr);
    }
    
    PetscReal dt;
    ierr = TSGetTimeStep(ts, &dt); CHKERRQ(ierr);
    
    // velocity ≈ (U - U_old) / dt
    ierr = VecCopy(U, velocity); CHKERRQ(ierr);
    ierr = VecAXPY(velocity, -1.0, U_old); CHKERRQ(ierr);
    if (dt > 0.0) {
        ierr = VecScale(velocity, 1.0 / dt); CHKERRQ(ierr);
    }
    
    // Check for transition
    ierr = manager->checkAndTransition(U, velocity, time); CHKERRQ(ierr);
    
    // Store current solution for next step
    ierr = VecCopy(U, U_old); CHKERRQ(ierr);
    time_old = time;
    
    ierr = VecDestroy(&velocity); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

PetscErrorCode IMEXPostStep(TS ts) {
    PetscFunctionBeginUser;
    // Post-step processing if needed
    PetscFunctionReturn(0);
}

PetscErrorCode IMEXAdapt(TS ts, PetscReal* dt, PetscBool* accept, void* ctx) {
    PetscFunctionBeginUser;
    
    auto* manager = static_cast<ImplicitExplicitTransitionManager*>(ctx);
    
    // Get current solution for adaptation
    Vec U;
    TSGetSolution(ts, &U);
    
    // Compute recommended dt based on current state
    // TODO: Get velocity vector properly
    Vec velocity;
    VecDuplicate(U, &velocity);
    VecSet(velocity, 0.0);
    
    double recommended_dt = manager->getRecommendedDt(U, velocity);
    
    // Clamp to mode-appropriate bounds
    if (manager->getCurrentMode() == IntegrationMode::IMPLICIT) {
        *dt = std::min(std::max(recommended_dt, manager->config.implicit_dt_min),
                      manager->config.implicit_dt_max);
    } else {
        *dt = std::min(std::max(recommended_dt, manager->config.explicit_dt_min),
                      manager->config.explicit_dt_max);
    }
    
    *accept = PETSC_TRUE;
    
    VecDestroy(&velocity);
    
    PetscFunctionReturn(0);
}

// ============================================================================
// ExplicitDynamicsHelper Implementation
// ============================================================================

ExplicitDynamicsHelper::ExplicitDynamicsHelper(MPI_Comm comm_in)
    : comm(comm_in), dm(nullptr), mass_inv(nullptr), work1(nullptr), work2(nullptr),
      density(2500.0), lambda(0.0), mu(0.0) {}

ExplicitDynamicsHelper::~ExplicitDynamicsHelper() {
    if (mass_inv) VecDestroy(&mass_inv);
    if (work1) VecDestroy(&work1);
    if (work2) VecDestroy(&work2);
}

PetscErrorCode ExplicitDynamicsHelper::setup(DM dm_in, double rho, 
                                             double E, double nu) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    dm = dm_in;
    density = rho;
    
    // Compute Lame parameters
    lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    mu = E / (2.0 * (1.0 + nu));
    
    // Create work vectors
    Vec local;
    ierr = DMGetLocalVector(dm, &local); CHKERRQ(ierr);
    ierr = DMGetGlobalVector(dm, &work1); CHKERRQ(ierr);
    ierr = VecDuplicate(work1, &work2); CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(dm, &local); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

PetscErrorCode ExplicitDynamicsHelper::createLumpedMass(Mat& M_lumped) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    // Create diagonal matrix for lumped mass
    // M_lumped = diag(ρ * V_node) where V_node is nodal volume
    
    // For simplicity, create a vector to store inverse diagonal
    ierr = VecDuplicate(work1, &mass_inv); CHKERRQ(ierr);
    ierr = VecSet(mass_inv, 1.0 / density); CHKERRQ(ierr);  // Simplified
    
    // TODO: Proper lumped mass computation based on element volumes
    
    PetscFunctionReturn(0);
}

PetscErrorCode ExplicitDynamicsHelper::computeInternalForce(Vec u, Vec f_int) {
    PetscFunctionBeginUser;
    // TODO: Implement internal force computation from displacement
    // f_int = K * u (in matrix form) or element-by-element assembly
    VecSet(f_int, 0.0);
    PetscFunctionReturn(0);
}

PetscErrorCode ExplicitDynamicsHelper::centralDifferenceStep(
    Vec u_new, Vec u_curr, Vec u_old, Vec f_ext, double dt) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    // Central difference: u_{n+1} = 2*u_n - u_{n-1} + dt² * M^{-1} * f
    
    // Compute internal force
    ierr = computeInternalForce(u_curr, work1); CHKERRQ(ierr);
    
    // f = f_ext - f_int
    ierr = VecWAXPY(work2, -1.0, work1, f_ext); CHKERRQ(ierr);
    
    // M^{-1} * f (using lumped mass)
    ierr = VecPointwiseMult(work1, mass_inv, work2); CHKERRQ(ierr);
    
    // u_{n+1} = 2*u_n - u_{n-1} + dt² * M^{-1} * f
    ierr = VecAXPBYPCZ(u_new, 2.0, -1.0, 0.0, u_curr, u_old); CHKERRQ(ierr);
    ierr = VecAXPY(u_new, dt * dt, work1); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

PetscErrorCode ExplicitDynamicsHelper::explicitNewmarkStep(
    Vec u_new, Vec u, Vec v, Vec a, Vec f_ext, double dt) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    // Explicit Newmark (β=0, γ=1/2):
    // u_{n+1} = u_n + dt*v_n + 0.5*dt²*a_n
    // a_{n+1} = M^{-1} * (f_ext - f_int(u_{n+1}))
    // v_{n+1} = v_n + 0.5*dt*(a_n + a_{n+1})
    
    // Predict u_{n+1}
    ierr = VecCopy(u, u_new); CHKERRQ(ierr);
    ierr = VecAXPY(u_new, dt, v); CHKERRQ(ierr);
    ierr = VecAXPY(u_new, 0.5 * dt * dt, a); CHKERRQ(ierr);
    
    // Compute new acceleration
    ierr = computeInternalForce(u_new, work1); CHKERRQ(ierr);
    ierr = VecWAXPY(work2, -1.0, work1, f_ext); CHKERRQ(ierr);
    ierr = VecPointwiseMult(a, mass_inv, work2); CHKERRQ(ierr);  // a_{n+1}
    
    // Update velocity (v_{n+1} = v_n + dt*a_{n+1} for γ=1/2 with central diff assumption)
    ierr = VecAXPY(v, dt, a); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

PetscErrorCode ExplicitDynamicsHelper::computeVelocity(Vec v, Vec u_new, Vec u_old, double dt) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    // v = (u_new - u_old) / (2*dt) for central difference
    ierr = VecWAXPY(v, -1.0, u_old, u_new); CHKERRQ(ierr);
    ierr = VecScale(v, 0.5 / dt); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

PetscErrorCode ExplicitDynamicsHelper::computeAcceleration(Vec a, Vec f_ext, Vec f_int, Vec f_damp) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    // a = M^{-1} * (f_ext - f_int - f_damp)
    ierr = VecCopy(f_ext, work1); CHKERRQ(ierr);
    ierr = VecAXPY(work1, -1.0, f_int); CHKERRQ(ierr);
    ierr = VecAXPY(work1, -1.0, f_damp); CHKERRQ(ierr);
    ierr = VecPointwiseMult(a, mass_inv, work1); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

PetscErrorCode ExplicitDynamicsHelper::applyDamping(Vec f_damp, Vec v, double alpha, double beta) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    // Rayleigh damping: f_damp = (α*M + β*K) * v
    // For lumped mass: f_damp ≈ α*M*v (mass-proportional part only)
    
    // Mass-proportional damping
    ierr = VecCopy(v, f_damp); CHKERRQ(ierr);
    ierr = VecScale(f_damp, alpha * density); CHKERRQ(ierr);  // Simplified
    
    // Stiffness-proportional damping would require K*v, which is expensive
    // Often omitted for explicit methods
    
    PetscFunctionReturn(0);
}

// ============================================================================
// SeismicEventDetector Implementation
// ============================================================================

SeismicEventDetector::SeismicEventDetector()
    : slip_rate_threshold(1e-3), min_magnitude(-2.0), min_duration(1e-4),
      event_active(false), just_started(false), just_ended(false),
      total_slip(0.0), total_moment(0.0), peak_slip_rate(0.0), event_start_time(0.0) {}

void SeismicEventDetector::configure(double slip_threshold, double min_mag, double min_dur) {
    slip_rate_threshold = slip_threshold;
    min_magnitude = min_mag;
    min_duration = min_dur;
}

void SeismicEventDetector::update(const FaultStressState& state, double time, double dt) {
    just_started = false;
    just_ended = false;
    
    bool is_seismic = state.slip_rate >= slip_rate_threshold;
    
    if (!event_active && is_seismic) {
        // Event just started
        event_active = true;
        just_started = true;
        event_start_time = time;
        total_slip = 0.0;
        total_moment = 0.0;
        peak_slip_rate = 0.0;
        
        current_event = SeismicEvent();
        current_event.time = time;
        current_event.hypo_x = 0.0;  // TODO: Get from fault
        current_event.hypo_y = 0.0;
        current_event.hypo_z = 0.0;
    }
    
    if (event_active) {
        // Accumulate during event
        double slip_increment = state.slip_rate * dt;
        total_slip += slip_increment;
        total_moment += state.moment;  // Should be moment rate * dt
        peak_slip_rate = std::max(peak_slip_rate, state.slip_rate);
        
        if (!is_seismic) {
            // Event just ended
            double duration = time - event_start_time;
            
            if (duration >= min_duration) {
                // Finalize event
                current_event.duration = duration;
                current_event.slip = total_slip;
                current_event.moment = total_moment;
                
                // Compute magnitude
                if (total_moment > 0.0) {
                    current_event.magnitude = (2.0/3.0) * std::log10(total_moment) - 6.07;
                } else {
                    current_event.magnitude = -99.0;
                }
                
                if (current_event.magnitude >= min_magnitude) {
                    catalog.push_back(current_event);
                    just_ended = true;
                }
            }
            
            event_active = false;
        }
    }
}

// ============================================================================
// KineticEnergyMonitor Implementation
// ============================================================================

KineticEnergyMonitor::KineticEnergyMonitor(int size)
    : history_size(size), peak_energy(0.0) {}

void KineticEnergyMonitor::addSample(double time, double energy) {
    history.push_back({time, energy});
    peak_energy = std::max(peak_energy, energy);
    
    while (static_cast<int>(history.size()) > history_size) {
        history.pop_front();
    }
}

bool KineticEnergyMonitor::isDecaying(double decay_rate_threshold) const {
    if (history.size() < 3) return false;
    
    // Fit linear trend to recent history
    int n = static_cast<int>(history.size());
    double sum_t = 0.0, sum_e = 0.0, sum_te = 0.0, sum_t2 = 0.0;
    
    for (const auto& sample : history) {
        sum_t += sample.first;
        sum_e += sample.second;
        sum_te += sample.first * sample.second;
        sum_t2 += sample.first * sample.first;
    }
    
    double slope = (n * sum_te - sum_t * sum_e) / (n * sum_t2 - sum_t * sum_t);
    
    return slope < -decay_rate_threshold;
}

bool KineticEnergyMonitor::isBelowThreshold(double ratio_threshold) const {
    if (peak_energy <= 0.0) return true;
    double current = getCurrentEnergy();
    return current <= peak_energy * ratio_threshold;
}

double KineticEnergyMonitor::getCurrentEnergy() const {
    if (history.empty()) return 0.0;
    return history.back().second;
}

void KineticEnergyMonitor::reset() {
    history.clear();
    peak_energy = 0.0;
}

} // namespace FSRM
