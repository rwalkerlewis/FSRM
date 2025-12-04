#ifndef IMPLICIT_EXPLICIT_TRANSITION_HPP
#define IMPLICIT_EXPLICIT_TRANSITION_HPP

#include "FSRM.hpp"
#include "FaultModel.hpp"
#include <petscts.h>
#include <petscvec.h>
#include <petscmat.h>
#include <vector>
#include <memory>
#include <functional>
#include <deque>

namespace FSRM {

/**
 * @brief Time integration mode
 * 
 * IMPLICIT: Backward Euler, BDF, or Generalized-alpha with large time steps
 *           Used for quasi-static problems (slow deformation, diffusion)
 * 
 * EXPLICIT: Central difference, explicit Newmark, or Runge-Kutta
 *           Used for wave propagation and dynamic rupture
 */
enum class IntegrationMode {
    IMPLICIT,       // Quasi-static (large dt)
    EXPLICIT        // Dynamic (small dt for CFL)
};

/**
 * @brief Trigger condition for mode transition
 */
enum class TransitionTrigger {
    NONE,                   // No trigger active
    STRESS_THRESHOLD,       // Von Mises or max shear stress exceeded
    COULOMB_FAILURE,        // Coulomb failure function >= 0
    SLIP_RATE,              // Fault slip rate exceeds seismic threshold
    VELOCITY,               // Particle velocity exceeds threshold
    ACCELERATION,           // Acceleration exceeds threshold
    ENERGY_RATE,            // Kinetic energy rate of change
    USER_DEFINED            // Custom trigger function
};

/**
 * @brief Criterion for returning to implicit mode
 */
enum class SettlingCriterion {
    TIME_ELAPSED,           // Fixed duration after trigger
    VELOCITY_DECAY,         // Velocity below threshold
    ENERGY_DECAY,           // Kinetic energy below threshold
    SLIP_RATE_DECAY,        // Slip rate back to interseismic values
    COMBINED                // Multiple criteria must be met
};

/**
 * @brief Explicit time integration scheme type for velocity computation
 */
enum class ExplicitSchemeType {
    CENTRAL_DIFFERENCE,     // v = (u_new - u_old) / (2*dt)
    NEWMARK,                // v_{n+1} = v_n + 0.5*dt*(a_n + a_{n+1})
    FORWARD_EULER           // v = (u_new - u_old) / dt
};

/**
 * @brief Internal configuration for implicit-explicit transition
 * 
 * This is the internal representation used by the IMEX manager.
 * For user-facing configuration, see ConfigReader::IMEXConfig.
 */
struct IMEXInternalConfig {
    // Enable/disable
    bool enabled = false;
    IntegrationMode initial_mode = IntegrationMode::IMPLICIT;
    
    // Implicit mode settings
    double implicit_dt_initial = 3600.0;    // Initial dt for implicit (s)
    double implicit_dt_min = 1.0;           // Minimum dt for implicit (s)
    double implicit_dt_max = 86400.0;       // Maximum dt for implicit (s)
    std::string implicit_method = "BEULER"; // TS method: BEULER, BDF, ALPHA2
    
    // Explicit mode settings
    double explicit_dt_initial = 1e-4;      // Initial dt for explicit (s)
    double explicit_dt_min = 1e-8;          // Minimum dt for explicit (s)
    double explicit_dt_max = 1e-2;          // Maximum dt for explicit (s)
    double cfl_factor = 0.5;                // CFL stability factor
    std::string explicit_method = "EULER";  // TS method: EULER, RK, NEWMARK
    
    // Trigger settings (implicit -> explicit)
    TransitionTrigger trigger_type = TransitionTrigger::COULOMB_FAILURE;
    double stress_threshold = 10.0e6;       // Pa (for STRESS_THRESHOLD)
    double coulomb_threshold = 0.0;         // Pa (for COULOMB_FAILURE, typically 0)
    double slip_rate_threshold = 1e-3;      // m/s (seismic slip rate)
    double velocity_threshold = 1e-4;       // m/s (particle velocity)
    double acceleration_threshold = 1e-2;   // m/s² (particle acceleration)
    double energy_rate_threshold = 1e3;     // J/s (kinetic energy rate)
    
    // Settling criteria (explicit -> implicit)
    SettlingCriterion settling_type = SettlingCriterion::COMBINED;
    double min_dynamic_duration = 1.0;      // Minimum time in explicit mode (s)
    double max_dynamic_duration = 60.0;     // Maximum time in explicit mode (s)
    double settling_velocity = 1e-6;        // Velocity threshold for settling (m/s)
    double settling_energy_ratio = 0.01;    // Energy ratio to peak for settling
    double settling_slip_rate = 1e-9;       // Slip rate threshold for settling (m/s)
    double settling_observation_window = 0.1; // Time to observe settling (s)
    
    // Adaptive time stepping within modes
    bool adaptive_timestep = true;
    double dt_growth_factor = 1.2;          // Max dt increase per step
    double dt_shrink_factor = 0.5;          // dt decrease on failure
    
    // Solver settings
    int implicit_max_iterations = 50;
    int explicit_max_iterations = 10;
    double implicit_rtol = 1e-6;
    double explicit_rtol = 1e-8;
    double implicit_atol = 1e-8;
    double explicit_atol = 1e-10;
    
    // Mass matrix handling
    bool use_lumped_mass = true;            // Lumped mass for explicit
    bool use_consistent_mass = false;       // Consistent mass for implicit
    
    // Output control
    int implicit_output_frequency = 10;     // Output every N implicit steps
    int explicit_output_frequency = 100;    // Output every N explicit steps
    bool log_transitions = true;            // Log mode transitions
    std::string transition_log_file = "imex_transitions.log";
    
    // Advanced settings
    bool smooth_transition = true;          // Ramp damping during transition
    double transition_ramp_time = 0.01;     // Damping ramp duration (s)
    bool preserve_solution = true;          // Interpolate solution at transition
    bool recompute_jacobian = true;         // Recompute Jacobian after transition
};

/**
 * @brief Statistics for a single integration phase
 */
struct PhaseStatistics {
    IntegrationMode mode;
    double start_time;
    double end_time;
    double duration;
    int num_steps;
    int num_nonlinear_iterations;
    int num_linear_iterations;
    double min_dt;
    double max_dt;
    double avg_dt;
    double peak_velocity;
    double peak_slip_rate;
    double total_slip;
    double released_energy;
    TransitionTrigger trigger;
    SettlingCriterion settling;
};

/**
 * @brief Record of a seismic event detected during explicit phase
 */
struct DynamicEvent {
    double trigger_time;
    double end_time;
    double duration;
    double peak_moment_rate;
    double total_moment;
    double magnitude;
    double peak_slip_rate;
    double max_slip;
    double hypo_x, hypo_y, hypo_z;
    std::vector<double> slip_history;
    std::vector<double> moment_rate_history;
};

/**
 * @brief Implicit-Explicit Time Integration Transition Manager
 * 
 * This class manages automatic transitions between implicit (quasi-static)
 * and explicit (dynamic) time integration schemes. This is essential for
 * induced seismicity modeling where:
 * 
 * 1. Pore pressure changes occur slowly (days-months) -> IMPLICIT
 * 2. Fault reaches criticality -> TRIGGER
 * 3. Dynamic rupture and wave propagation (seconds) -> EXPLICIT  
 * 4. Waves dissipate, system settles -> IMPLICIT
 * 
 * The manager monitors physical quantities (stress, slip rate, velocity)
 * to detect when transitions should occur, and reconfigures the PETSc
 * time stepper (TS) accordingly.
 */
class ImplicitExplicitTransitionManager {
public:
    ImplicitExplicitTransitionManager(MPI_Comm comm);
    ~ImplicitExplicitTransitionManager();
    
    /**
     * @brief Initialize with configuration
     */
    PetscErrorCode initialize(const IMEXInternalConfig& config);
    
    /**
     * @brief Setup with PETSc TS and solution vectors
     */
    PetscErrorCode setup(TS ts, DM dm, Vec solution, Mat mass_matrix = nullptr);
    
    /**
     * @brief Set fault model for slip rate monitoring
     */
    void setFaultModel(SeismicFaultModel* fault);
    void setFaultNetwork(FaultNetwork* network);
    
    /**
     * @brief Set custom trigger function
     * 
     * @param func Function returning true when transition should occur
     */
    void setCustomTrigger(std::function<bool(Vec, double)> func);
    
    /**
     * @brief Set custom settling function
     * 
     * @param func Function returning true when system has settled
     */
    void setCustomSettling(std::function<bool(Vec, double, double)> func);
    
    /**
     * @brief Main check called at each time step
     * 
     * Evaluates transition criteria and performs mode switch if needed.
     * Should be called from the TS monitor or post-step function.
     */
    PetscErrorCode checkAndTransition(Vec solution, Vec velocity, double time);
    
    /**
     * @brief Force transition to specific mode
     */
    PetscErrorCode forceTransition(IntegrationMode target_mode);
    
    /**
     * @brief Get current integration mode
     */
    IntegrationMode getCurrentMode() const { return current_mode; }
    
    /**
     * @brief Check if system is in dynamic event
     */
    bool isInDynamicEvent() const { return in_dynamic_event; }
    
    /**
     * @brief Get current time step (mode-appropriate)
     */
    double getCurrentDt() const { return current_dt; }
    
    /**
     * @brief Get recommended time step based on current state
     */
    double getRecommendedDt(Vec solution, Vec velocity);
    
    /**
     * @brief Compute CFL-limited time step for explicit mode
     */
    double computeCFLTimeStep(double min_cell_size, double max_wave_speed);
    
    /**
     * @brief Get the IMEX configuration (read-only access for callbacks)
     */
    const IMEXInternalConfig& getConfig() const { return config; }
    
    /**
     * @brief Get phase statistics
     */
    const std::vector<PhaseStatistics>& getPhaseHistory() const { return phase_history; }
    
    /**
     * @brief Get detected dynamic events
     */
    const std::vector<DynamicEvent>& getDynamicEvents() const { return dynamic_events; }
    
    /**
     * @brief Write transition log
     */
    PetscErrorCode writeTransitionLog();
    
    /**
     * @brief Get summary statistics
     */
    void getSummaryStatistics(int& num_implicit_phases, int& num_explicit_phases,
                             double& total_implicit_time, double& total_explicit_time,
                             int& num_events, double& total_moment);
    
    /**
     * @brief Get previous solution for velocity computation
     * 
     * Used by IMEXMonitor to compute velocity via finite differences.
     * Returns nullptr if no previous solution has been stored.
     */
    Vec getPreviousSolution() const { return solution_prev; }
    
    /**
     * @brief Store current solution as previous for next step
     * 
     * @param U Current solution to copy
     */
    PetscErrorCode storePreviousSolution(Vec U);
    
    // ========================================================================
    // Internal helper methods
    // ========================================================================
    
private:
    MPI_Comm comm;
    int rank, size;
    
    // Configuration
    IMEXInternalConfig config;
    
    // PETSc objects
    TS ts;
    DM dm;
    Vec solution;
    Vec solution_prev;      // Previous solution for velocity computation
    Vec velocity_old;       // For velocity-based criteria
    Vec acceleration;       // For acceleration-based criteria  
    Mat mass_matrix;
    Mat mass_lumped;        // Lumped mass for explicit
    
    // Current state
    IntegrationMode current_mode;
    double current_dt;
    double current_time;
    int current_step;
    bool in_dynamic_event;
    double event_start_time;
    double peak_kinetic_energy;
    
    // Fault models (for slip rate monitoring)
    SeismicFaultModel* fault_model;
    FaultNetwork* fault_network;
    
    // Custom functions
    std::function<bool(Vec, double)> custom_trigger;
    std::function<bool(Vec, double, double)> custom_settling;
    
    // History tracking
    std::vector<PhaseStatistics> phase_history;
    PhaseStatistics current_phase;
    std::vector<DynamicEvent> dynamic_events;
    DynamicEvent current_event;
    
    // Settling detection buffers
    std::deque<double> velocity_history;
    std::deque<double> energy_history;
    std::deque<double> slip_rate_history;
    
    // Internal methods
    PetscErrorCode configureImplicitMode();
    PetscErrorCode configureExplicitMode();
    
    bool checkTriggerCondition(Vec solution, Vec velocity, double time);
    bool checkSettlingCondition(Vec solution, Vec velocity, double time);
    
    double computeMaxVelocity(Vec velocity);
    double computeMaxAcceleration(Vec velocity);
    double computeKineticEnergy(Vec velocity);
    double computeMaxSlipRate();
    double computeCoulombFailure(Vec solution);
    double computeVonMisesStress(Vec solution);
    
    void getGridAndMaterialParameters(double& min_cell_size, double& max_wave_speed);
    
    PetscErrorCode interpolateSolution(double old_dt, double new_dt);
    PetscErrorCode adjustDamping(double ramp_factor);
    
    void startNewPhase(IntegrationMode mode, double time);
    void endCurrentPhase(double time, TransitionTrigger trigger, SettlingCriterion settling);
    
    void startDynamicEvent(double time, double x, double y, double z);
    void updateDynamicEvent(double slip_rate, double moment_rate);
    void endDynamicEvent(double time);
    
    PetscErrorCode logTransition(IntegrationMode from, IntegrationMode to, 
                                 double time, const std::string& reason);
};

/**
 * @brief PETSc TS monitor function for IMEX transitions
 * 
 * This function should be registered with TSMonitorSet to enable
 * automatic transition checking at each time step.
 */
PetscErrorCode IMEXMonitor(TS ts, PetscInt step, PetscReal time, Vec U, void* ctx);

/**
 * @brief PETSc TS post-step function for IMEX
 */
PetscErrorCode IMEXPostStep(TS ts);

/**
 * @brief PETSc TS adapt function that respects IMEX constraints
 */
PetscErrorCode IMEXAdapt(TS ts, PetscReal* dt, PetscBool* accept, void* ctx);

/**
 * @brief Helper class for explicit time integration with lumped mass
 * 
 * Provides efficient matrix-free operations for explicit dynamics
 * using lumped mass matrices.
 */
class ExplicitDynamicsHelper {
public:
    ExplicitDynamicsHelper(MPI_Comm comm);
    ~ExplicitDynamicsHelper();
    
    /**
     * @brief Setup with DM and material properties
     */
    PetscErrorCode setup(DM dm, double density, double youngs_modulus, 
                        double poisson_ratio);
    
    /**
     * @brief Create lumped mass matrix
     */
    PetscErrorCode createLumpedMass(Mat& M_lumped);
    
    /**
     * @brief Compute internal forces (stiffness * displacement)
     */
    PetscErrorCode computeInternalForce(Vec u, Vec f_int);
    
    /**
     * @brief Central difference time step
     * 
     * u_{n+1} = 2*u_n - u_{n-1} + dt² * M^{-1} * (f_ext - f_int - C*v)
     */
    PetscErrorCode centralDifferenceStep(Vec u_new, Vec u_curr, Vec u_old,
                                         Vec f_ext, double dt);
    
    /**
     * @brief Explicit Newmark step (beta=0, gamma=0.5)
     */
    PetscErrorCode explicitNewmarkStep(Vec u_new, Vec u, Vec v, Vec a,
                                       Vec f_ext, double dt);
    
    /**
     * @brief Compute velocity from displacement history
     * 
     * @param v Output velocity vector
     * @param u_new Current displacement
     * @param u_old Previous displacement
     * @param v_old Previous velocity (for Newmark scheme, can be nullptr)
     * @param a_old Previous acceleration (for Newmark scheme, can be nullptr)
     * @param a_new Current acceleration (for Newmark scheme, can be nullptr)
     * @param dt Time step
     * @param scheme Integration scheme type
     */
    PetscErrorCode computeVelocity(Vec v, Vec u_new, Vec u_old, Vec v_old, 
                                   Vec a_old, Vec a_new, double dt,
                                   ExplicitSchemeType scheme = ExplicitSchemeType::CENTRAL_DIFFERENCE);
    
    /**
     * @brief Compute acceleration from force balance
     */
    PetscErrorCode computeAcceleration(Vec a, Vec f_ext, Vec f_int, Vec f_damp);
    
    /**
     * @brief Apply damping forces (mass-proportional Rayleigh damping)
     * 
     * Note: Stiffness-proportional damping (β*K*v) is intentionally omitted for
     * explicit methods as it is expensive to compute and reduces the stable time step.
     * 
     * @param f_damp Output damping force vector
     * @param v Velocity vector
     * @param alpha Mass-proportional damping coefficient
     */
    PetscErrorCode applyDamping(Vec f_damp, Vec v, double alpha);
    
private:
    MPI_Comm comm;
    DM dm;
    Vec mass_inv;       // Inverse lumped mass (diagonal)
    Vec work1, work2;   // Work vectors
    double density;
    double lambda, mu;  // Lame parameters
};

/**
 * @brief Seismic event detector for IMEX transitions
 * 
 * Monitors fault slip and determines when a seismic event has occurred,
 * computing event properties (magnitude, slip, moment).
 */
class SeismicEventDetector {
public:
    SeismicEventDetector();
    
    /**
     * @brief Configure detection parameters
     */
    void configure(double slip_rate_threshold, double min_magnitude,
                  double min_duration);
    
    /**
     * @brief Update with current fault state
     */
    void update(const FaultStressState& state, double time, double dt);
    
    /**
     * @brief Check if event is currently active
     */
    bool isEventActive() const { return event_active; }
    
    /**
     * @brief Check if event just started
     */
    bool eventJustStarted() const { return just_started; }
    
    /**
     * @brief Check if event just ended
     */
    bool eventJustEnded() const { return just_ended; }
    
    /**
     * @brief Get current event (if active)
     */
    const SeismicEvent& getCurrentEvent() const { return current_event; }
    
    /**
     * @brief Get event catalog
     */
    const std::vector<SeismicEvent>& getCatalog() const { return catalog; }
    
private:
    double slip_rate_threshold;
    double min_magnitude;
    double min_duration;
    
    bool event_active;
    bool just_started;
    bool just_ended;
    
    SeismicEvent current_event;
    std::vector<SeismicEvent> catalog;
    
    // Accumulation during event
    double total_slip;
    double total_moment;
    double peak_slip_rate;
    double event_start_time;
};

/**
 * @brief Kinetic energy monitor for settling detection
 */
class KineticEnergyMonitor {
public:
    KineticEnergyMonitor(int history_size = 100);
    
    /**
     * @brief Add new energy sample
     */
    void addSample(double time, double energy);
    
    /**
     * @brief Check if energy is decaying
     */
    bool isDecaying(double decay_rate_threshold = 0.1) const;
    
    /**
     * @brief Check if energy is below threshold relative to peak
     */
    bool isBelowThreshold(double ratio_threshold = 0.01) const;
    
    /**
     * @brief Get peak energy since last reset
     */
    double getPeakEnergy() const { return peak_energy; }
    
    /**
     * @brief Get current energy
     */
    double getCurrentEnergy() const;
    
    /**
     * @brief Reset (call when entering explicit mode)
     */
    void reset();
    
private:
    int history_size;
    std::deque<std::pair<double, double>> history;  // (time, energy)
    double peak_energy;
};

} // namespace FSRM

#endif // IMPLICIT_EXPLICIT_TRANSITION_HPP
