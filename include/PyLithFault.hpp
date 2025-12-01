#ifndef PYLITH_FAULT_HPP
#define PYLITH_FAULT_HPP

/**
 * @file PyLithFault.hpp
 * @brief PyLith-style fault representation using cohesive cells
 * 
 * This module implements fault modeling following the PyLith approach
 * (https://github.com/geodynamics/pylith), which uses:
 * 
 * - Cohesive cells: Zero-volume cells inserted along fault surfaces
 * - Split nodes: Nodes duplicated on either side of the fault
 * - Lagrange multipliers: For enforcing traction/slip conditions
 * - Kinematic faults: Prescribed slip history
 * - Dynamic faults: Spontaneous rupture with friction laws
 * 
 * Key classes:
 * - FaultCohesive: Base class for all cohesive fault types
 * - FaultCohesiveKin: Kinematic fault with prescribed slip
 * - FaultCohesiveDyn: Dynamic fault with spontaneous rupture
 * - SlipTimeFn: Base class for slip time functions
 * - FrictionModel: Friction interface for dynamic rupture
 * 
 * References:
 * - Aagaard et al. (2013): A domain decomposition approach to implementing
 *   fault slip in finite-element models of quasi-static and dynamic crustal
 *   deformation
 * - PyLith manual: https://geodynamics.org/resources/pylith
 */

#include <string>
#include <vector>
#include <map>
#include <memory>
#include <cmath>
#include <array>
#include <functional>

namespace FSRM {

// Forward declarations
class FaultCohesive;
class FaultCohesiveKin;
class FaultCohesiveDyn;

// =============================================================================
// Enumerations
// =============================================================================

/**
 * @brief Types of cohesive fault implementations
 */
enum class CohesiveFaultType {
    KINEMATIC,          // Prescribed slip history
    DYNAMIC,            // Spontaneous rupture with friction
    IMPULSE,            // Impulsive sources for Green's functions
    ABSORBING           // Absorbing boundary conditions on fault
};

/**
 * @brief Types of slip time functions for kinematic faults
 */
enum class SlipTimeFnType {
    STEP,               // Heaviside step function
    RAMP,               // Linear ramp
    BRUNE,              // Brune-style exponential
    LIU,                // Liu et al. regularized slip rate
    TIME_HISTORY,       // User-supplied time history
    CONSTANT_RATE       // Constant slip rate (creeping)
};

/**
 * @brief Traction perturbation types for dynamic faults
 */
enum class TractionPerturbationType {
    UNIFORM,            // Spatially uniform perturbation
    GAUSSIAN,           // Gaussian spatial distribution
    ASPERITY,           // Asperity-based heterogeneity
    SPATIALLY_VARIABLE  // Fully heterogeneous from data
};

/**
 * @brief Output fields available for faults
 */
enum class FaultOutputField {
    SLIP,               // Slip vector (m)
    SLIP_RATE,          // Slip rate vector (m/s)
    TRACTION,           // Traction vector (Pa)
    TRACTION_CHANGE,    // Change in traction (Pa)
    STATE_VARIABLE,     // Friction state variable
    RUPTURE_TIME,       // Rupture initiation time (s)
    RISE_TIME,          // Slip duration (s)
    VERTEX_LAGRANGE_MULTIPLIER  // Lagrange multiplier at vertices
};

// =============================================================================
// Slip Time Function Classes (Kinematic Sources)
// =============================================================================

/**
 * @brief Base class for slip time functions
 * 
 * Slip time functions define how slip evolves in time for kinematic faults.
 * The slip at a point is: slip(t) = amplitude * f(t - slip_time)
 * where f is the normalized slip time function.
 */
class SlipTimeFn {
public:
    SlipTimeFn(SlipTimeFnType type) : fn_type(type) {}
    virtual ~SlipTimeFn() = default;
    
    /**
     * @brief Compute slip at given time
     * @param t Current time (s)
     * @param slip_time Initiation time for slip at this point (s)
     * @param amplitude Final slip amplitude (m)
     * @return Current slip value (m)
     */
    virtual double computeSlip(double t, double slip_time, double amplitude) const = 0;
    
    /**
     * @brief Compute slip rate at given time
     * @param t Current time (s)
     * @param slip_time Initiation time for slip (s)
     * @param amplitude Final slip amplitude (m)
     * @return Current slip rate (m/s)
     */
    virtual double computeSlipRate(double t, double slip_time, double amplitude) const = 0;
    
    /**
     * @brief Configure from key-value pairs
     */
    virtual void configure(const std::map<std::string, std::string>& config) = 0;
    
    SlipTimeFnType getType() const { return fn_type; }
    
protected:
    SlipTimeFnType fn_type;
};

/**
 * @brief Step slip time function (Heaviside)
 * 
 * Instantaneous slip at rupture time.
 * f(τ) = H(τ) where τ = t - slip_time
 */
class SlipTimeFnStep : public SlipTimeFn {
public:
    SlipTimeFnStep() : SlipTimeFn(SlipTimeFnType::STEP) {}
    
    double computeSlip(double t, double slip_time, double amplitude) const override {
        return (t >= slip_time) ? amplitude : 0.0;
    }
    
    double computeSlipRate(double t, double slip_time, double amplitude) const override {
        // Delta function - return large value at slip_time
        double dt = 1e-6;  // Small time window
        if (std::abs(t - slip_time) < dt) {
            return amplitude / dt;
        }
        return 0.0;
    }
    
    void configure(const std::map<std::string, std::string>&) override {}
};

/**
 * @brief Linear ramp slip time function
 * 
 * Linear increase from 0 to amplitude over rise_time.
 * f(τ) = min(τ/rise_time, 1) for τ >= 0
 */
class SlipTimeFnRamp : public SlipTimeFn {
public:
    SlipTimeFnRamp() : SlipTimeFn(SlipTimeFnType::RAMP), rise_time(1.0) {}
    
    void setRiseTime(double tau) { rise_time = tau; }
    double getRiseTime() const { return rise_time; }
    
    double computeSlip(double t, double slip_time, double amplitude) const override {
        double tau = t - slip_time;
        if (tau <= 0.0) return 0.0;
        if (tau >= rise_time) return amplitude;
        return amplitude * tau / rise_time;
    }
    
    double computeSlipRate(double t, double slip_time, double amplitude) const override {
        double tau = t - slip_time;
        if (tau <= 0.0 || tau >= rise_time) return 0.0;
        return amplitude / rise_time;
    }
    
    void configure(const std::map<std::string, std::string>& config) override;
    
private:
    double rise_time;  // Duration of slip (s)
};

/**
 * @brief Brune-style exponential slip time function
 * 
 * f(τ) = 1 - (1 + τ/τ_c) * exp(-τ/τ_c) for τ >= 0
 * 
 * This models an exponentially decaying slip rate, commonly used
 * in seismological source models.
 */
class SlipTimeFnBrune : public SlipTimeFn {
public:
    SlipTimeFnBrune() : SlipTimeFn(SlipTimeFnType::BRUNE), tau_c(1.0) {}
    
    void setCharacteristicTime(double tc) { tau_c = tc; }
    double getCharacteristicTime() const { return tau_c; }
    
    double computeSlip(double t, double slip_time, double amplitude) const override {
        double tau = t - slip_time;
        if (tau <= 0.0) return 0.0;
        double x = tau / tau_c;
        return amplitude * (1.0 - (1.0 + x) * std::exp(-x));
    }
    
    double computeSlipRate(double t, double slip_time, double amplitude) const override {
        double tau = t - slip_time;
        if (tau <= 0.0) return 0.0;
        double x = tau / tau_c;
        return amplitude * x / tau_c * std::exp(-x);
    }
    
    void configure(const std::map<std::string, std::string>& config) override;
    
private:
    double tau_c;  // Characteristic time (s)
};

/**
 * @brief Liu et al. (2006) regularized slip rate function
 * 
 * Smooth, regularized slip rate with controlled rise time.
 * More physical than Brune for large earthquakes.
 */
class SlipTimeFnLiu : public SlipTimeFn {
public:
    SlipTimeFnLiu() : SlipTimeFn(SlipTimeFnType::LIU), 
                     rise_time(1.0), cnorm(1.0) {
        computeNormalization();
    }
    
    void setRiseTime(double tr) { 
        rise_time = tr; 
        computeNormalization();
    }
    double getRiseTime() const { return rise_time; }
    
    double computeSlip(double t, double slip_time, double amplitude) const override;
    double computeSlipRate(double t, double slip_time, double amplitude) const override;
    
    void configure(const std::map<std::string, std::string>& config) override;
    
private:
    double rise_time;
    double cnorm;  // Normalization constant
    
    void computeNormalization();
};

/**
 * @brief Constant slip rate function (creeping fault)
 * 
 * f(τ) = slip_rate * τ for τ >= 0
 * 
 * Used for modeling steady-state aseismic creep.
 */
class SlipTimeFnConstantRate : public SlipTimeFn {
public:
    SlipTimeFnConstantRate() : SlipTimeFn(SlipTimeFnType::CONSTANT_RATE), 
                               slip_rate(1e-9) {}
    
    void setSlipRate(double v) { slip_rate = v; }
    double getSlipRate() const { return slip_rate; }
    
    double computeSlip(double t, double slip_time, double amplitude) const override {
        (void)amplitude;  // Not used - slip grows indefinitely
        double tau = t - slip_time;
        if (tau <= 0.0) return 0.0;
        return slip_rate * tau;
    }
    
    double computeSlipRate(double t, double slip_time, double amplitude) const override {
        (void)amplitude;
        return (t >= slip_time) ? slip_rate : 0.0;
    }
    
    void configure(const std::map<std::string, std::string>& config) override;
    
private:
    double slip_rate;  // Constant rate (m/s)
};

/**
 * @brief User-supplied time history
 */
class SlipTimeFnTimeHistory : public SlipTimeFn {
public:
    SlipTimeFnTimeHistory() : SlipTimeFn(SlipTimeFnType::TIME_HISTORY) {}
    
    void setTimeHistory(const std::vector<double>& times,
                       const std::vector<double>& values) {
        time_points = times;
        slip_values = values;
    }
    
    double computeSlip(double t, double slip_time, double amplitude) const override;
    double computeSlipRate(double t, double slip_time, double amplitude) const override;
    
    void configure(const std::map<std::string, std::string>& config) override;
    
private:
    std::vector<double> time_points;
    std::vector<double> slip_values;  // Normalized (0 to 1)
    
    double interpolate(double t, double slip_time) const;
};

// =============================================================================
// Cohesive Cell Data Structures
// =============================================================================

/**
 * @brief Vertex on fault surface
 */
struct FaultVertex {
    int vertex_id;              // Global vertex index
    int vertex_negative;        // Vertex on negative side of fault
    int vertex_positive;        // Vertex on positive side of fault
    int vertex_lagrange;        // Lagrange multiplier vertex
    
    std::array<double, 3> coords;       // Original coordinates
    std::array<double, 3> normal;       // Fault normal at vertex
    std::array<double, 3> along_strike; // Strike direction
    std::array<double, 3> up_dip;       // Up-dip direction
    
    double area;                // Area associated with vertex
    
    FaultVertex() : vertex_id(-1), vertex_negative(-1), vertex_positive(-1),
                   vertex_lagrange(-1), coords({0,0,0}), normal({0,0,1}),
                   along_strike({1,0,0}), up_dip({0,1,0}), area(0.0) {}
};

/**
 * @brief Cohesive cell (zero-volume element) along fault
 */
struct CohesiveCell {
    int cell_id;                // Global cell index
    std::vector<int> vertices_negative;  // Vertices on - side
    std::vector<int> vertices_positive;  // Vertices on + side
    std::vector<int> vertices_lagrange;  // Lagrange multiplier vertices
    
    double area;                // Cell area
    std::array<double, 3> centroid;     // Cell center
    std::array<double, 3> normal;       // Average normal
    
    CohesiveCell() : cell_id(-1), area(0.0), centroid({0,0,0}), normal({0,0,1}) {}
};

/**
 * @brief Slip field on fault surface
 */
struct FaultSlipField {
    std::vector<double> slip_left_lateral;   // Left-lateral slip (m)
    std::vector<double> slip_reverse;        // Reverse/thrust slip (m)
    std::vector<double> slip_opening;        // Opening/normal slip (m)
    
    std::vector<double> slip_rate_ll;        // Left-lateral rate (m/s)
    std::vector<double> slip_rate_reverse;   // Reverse slip rate (m/s)
    std::vector<double> slip_rate_opening;   // Opening rate (m/s)
    
    void resize(size_t n) {
        slip_left_lateral.resize(n, 0.0);
        slip_reverse.resize(n, 0.0);
        slip_opening.resize(n, 0.0);
        slip_rate_ll.resize(n, 0.0);
        slip_rate_reverse.resize(n, 0.0);
        slip_rate_opening.resize(n, 0.0);
    }
    
    size_t size() const { return slip_left_lateral.size(); }
};

/**
 * @brief Traction field on fault surface
 */
struct FaultTractionField {
    std::vector<double> traction_shear_ll;   // Left-lateral shear traction (Pa)
    std::vector<double> traction_shear_ud;   // Up-dip shear traction (Pa)
    std::vector<double> traction_normal;     // Normal traction (Pa)
    
    void resize(size_t n) {
        traction_shear_ll.resize(n, 0.0);
        traction_shear_ud.resize(n, 0.0);
        traction_normal.resize(n, 0.0);
    }
    
    size_t size() const { return traction_shear_ll.size(); }
};

/**
 * @brief Kinematic slip parameters at a vertex
 */
struct KinematicSlipParams {
    double slip_time;           // Rupture initiation time (s)
    double rise_time;           // Slip duration (s)
    double final_slip_ll;       // Final left-lateral slip (m)
    double final_slip_reverse;  // Final reverse slip (m)
    double final_slip_opening;  // Final opening (m)
    
    KinematicSlipParams() : slip_time(0.0), rise_time(1.0),
                           final_slip_ll(0.0), final_slip_reverse(0.0),
                           final_slip_opening(0.0) {}
};

/**
 * @brief Dynamic fault state at a vertex
 */
struct DynamicFaultState {
    double state_variable;      // Friction state variable (θ)
    double slip_rate;           // Current slip rate magnitude (m/s)
    double effective_normal;    // Effective normal stress (Pa)
    double shear_stress;        // Shear stress magnitude (Pa)
    double friction;            // Current friction coefficient
    bool is_locked;             // Currently locked (stick)?
    bool has_ruptured;          // Has rupture initiated?
    double rupture_time;        // Time of rupture initiation (s)
    
    DynamicFaultState() : state_variable(1e6), slip_rate(0.0),
                         effective_normal(0.0), shear_stress(0.0),
                         friction(0.6), is_locked(true), has_ruptured(false),
                         rupture_time(-1.0) {}
};

// =============================================================================
// Friction Models for Dynamic Faults
// =============================================================================

/**
 * @brief Base friction model interface for dynamic faults
 */
class FaultFrictionModel {
public:
    virtual ~FaultFrictionModel() = default;
    
    /**
     * @brief Compute friction coefficient
     */
    virtual double computeFriction(const DynamicFaultState& state) const = 0;
    
    /**
     * @brief Compute state variable evolution rate
     */
    virtual double computeStateRate(const DynamicFaultState& state) const = 0;
    
    /**
     * @brief Update state for given slip increment
     */
    virtual void updateState(DynamicFaultState& state, double dt) const = 0;
    
    /**
     * @brief Configure from parameters
     */
    virtual void configure(const std::map<std::string, std::string>& config) = 0;
    
    /**
     * @brief Get friction model name
     */
    virtual std::string name() const = 0;
};

/**
 * @brief Static friction model
 */
class StaticFriction : public FaultFrictionModel {
public:
    StaticFriction() : coef(0.6) {}
    
    double computeFriction(const DynamicFaultState&) const override {
        return coef;
    }
    
    double computeStateRate(const DynamicFaultState&) const override { return 0.0; }
    void updateState(DynamicFaultState&, double) const override {}
    
    void configure(const std::map<std::string, std::string>& config) override;
    std::string name() const override { return "static"; }
    
    void setFrictionCoefficient(double f) { coef = f; }
    
private:
    double coef;
};

/**
 * @brief Slip-weakening friction model
 */
class SlipWeakeningFriction : public FaultFrictionModel {
public:
    SlipWeakeningFriction() : static_coef(0.6), dynamic_coef(0.4), 
                              d0(0.4), cohesion(0.0) {}
    
    double computeFriction(const DynamicFaultState& state) const override;
    double computeStateRate(const DynamicFaultState& state) const override;
    void updateState(DynamicFaultState& state, double dt) const override;
    
    void configure(const std::map<std::string, std::string>& config) override;
    std::string name() const override { return "slip_weakening"; }
    
    void setStaticCoefficient(double fs) { static_coef = fs; }
    void setDynamicCoefficient(double fd) { dynamic_coef = fd; }
    void setCriticalSlipDistance(double d) { d0 = d; }
    void setCohesion(double c) { cohesion = c; }
    
    double getStaticCoefficient() const { return static_coef; }
    double getDynamicCoefficient() const { return dynamic_coef; }
    double getCriticalSlipDistance() const { return d0; }
    double getCohesion() const { return cohesion; }
    
private:
    double static_coef;
    double dynamic_coef;
    double d0;          // Critical slip distance (m)
    double cohesion;    // Cohesive strength (Pa)
};

/**
 * @brief Rate-and-state friction with aging law
 */
class RateStateFrictionAging : public FaultFrictionModel {
public:
    RateStateFrictionAging() : a(0.01), b(0.015), dc(0.01), f0(0.6),
                               v0(1e-6), linear_threshold(1e-12) {}
    
    double computeFriction(const DynamicFaultState& state) const override;
    double computeStateRate(const DynamicFaultState& state) const override;
    void updateState(DynamicFaultState& state, double dt) const override;
    
    void configure(const std::map<std::string, std::string>& config) override;
    std::string name() const override { return "rate_state_aging"; }
    
    // Parameter setters
    void setDirectEffect(double a_) { a = a_; }
    void setEvolutionEffect(double b_) { b = b_; }
    void setCriticalSlipDistance(double dc_) { dc = dc_; }
    void setReferenceFriction(double f0_) { f0 = f0_; }
    void setReferenceSlipRate(double v0_) { v0 = v0_; }
    
    // Stability analysis
    bool isVelocityWeakening() const { return b > a; }
    double getCriticalStiffness(double sigma_n) const;
    double getNucleationLength(double G, double sigma_n) const;
    
private:
    double a;           // Direct effect parameter
    double b;           // Evolution effect parameter
    double dc;          // Critical slip distance (m)
    double f0;          // Reference friction
    double v0;          // Reference slip rate (m/s)
    double linear_threshold;  // Linearization threshold
};

/**
 * @brief Rate-and-state friction with slip law
 */
class RateStateFrictionSlip : public FaultFrictionModel {
public:
    RateStateFrictionSlip() : a(0.01), b(0.015), dc(0.01), f0(0.6),
                              v0(1e-6), linear_threshold(1e-12) {}
    
    double computeFriction(const DynamicFaultState& state) const override;
    double computeStateRate(const DynamicFaultState& state) const override;
    void updateState(DynamicFaultState& state, double dt) const override;
    
    void configure(const std::map<std::string, std::string>& config) override;
    std::string name() const override { return "rate_state_slip"; }
    
    void setDirectEffect(double a_) { a = a_; }
    void setEvolutionEffect(double b_) { b = b_; }
    void setCriticalSlipDistance(double dc_) { dc = dc_; }
    void setReferenceFriction(double f0_) { f0 = f0_; }
    void setReferenceSlipRate(double v0_) { v0 = v0_; }
    
private:
    double a, b, dc, f0, v0;
    double linear_threshold;
};

// =============================================================================
// Base Cohesive Fault Class
// =============================================================================

/**
 * @brief Base class for cohesive fault implementations
 * 
 * Implements the common interface and functionality for all fault types.
 * Handles mesh topology, vertex/cell data, and output.
 */
class FaultCohesive {
public:
    FaultCohesive(CohesiveFaultType type);
    virtual ~FaultCohesive() = default;
    
    // Configuration
    virtual void configure(const std::map<std::string, std::string>& config);
    void setName(const std::string& n) { name = n; }
    std::string getName() const { return name; }
    
    // Mesh topology
    void setFaultVertices(const std::vector<FaultVertex>& verts);
    void setCohesiveCells(const std::vector<CohesiveCell>& cells);
    
    /**
     * @brief Create cohesive cells from fault surface definition
     * @param node_coords Node coordinates [x0,y0,z0, x1,y1,z1, ...]
     * @param fault_node_indices Indices of nodes on fault surface
     * @param fault_orientation Strike and dip angles (radians)
     */
    void createCohesiveCells(const std::vector<double>& node_coords,
                            const std::vector<int>& fault_node_indices,
                            double strike, double dip);
    
    // Integration with solver
    virtual void initialize() = 0;
    virtual void prestep(double t, double dt) = 0;
    virtual void computeResidual(const double* solution, double* residual) = 0;
    virtual void computeJacobian(const double* solution, double* jacobian) = 0;
    virtual void poststep(double t, double dt) = 0;
    
    // Access to fields
    const FaultSlipField& getSlipField() const { return slip_field; }
    const FaultTractionField& getTractionField() const { return traction_field; }
    
    // Query functions
    size_t numVertices() const { return vertices.size(); }
    size_t numCells() const { return cells.size(); }
    const FaultVertex& getVertex(size_t i) const { return vertices[i]; }
    const CohesiveCell& getCell(size_t i) const { return cells[i]; }
    
    // Output
    void setOutputFields(const std::vector<FaultOutputField>& fields) {
        output_fields = fields;
    }
    virtual void writeOutput(const std::string& filename, double time);
    
    CohesiveFaultType getType() const { return fault_type; }
    
protected:
    std::string name;
    CohesiveFaultType fault_type;
    
    // Mesh data
    std::vector<FaultVertex> vertices;
    std::vector<CohesiveCell> cells;
    
    // Field data
    FaultSlipField slip_field;
    FaultTractionField traction_field;
    
    // Output configuration
    std::vector<FaultOutputField> output_fields;
    
    // Material properties for traction computation
    double shear_modulus;
    double bulk_modulus;
    
    // Helper functions
    void computeFaultOrientation(const FaultVertex& v, 
                                std::array<double, 3>& normal,
                                std::array<double, 3>& strike_dir,
                                std::array<double, 3>& dip_dir) const;
    
    void transformToFaultCoords(const double* global_vector,
                               const FaultVertex& v,
                               double& ll_component,
                               double& ud_component,
                               double& normal_component) const;
    
    void transformToGlobalCoords(double ll, double ud, double normal,
                                const FaultVertex& v,
                                double* global_vector) const;
};

// =============================================================================
// Kinematic Fault (FaultCohesiveKin)
// =============================================================================

/**
 * @brief Kinematic fault with prescribed slip history
 * 
 * Implements PyLith's FaultCohesiveKin:
 * - Slip is prescribed as a function of space and time
 * - Lagrange multipliers enforce the kinematic constraint
 * - Multiple slip time functions supported
 * 
 * Usage:
 * 1. Set fault geometry via setFaultVertices() or createCohesiveCells()
 * 2. Set slip time function via setSlipTimeFn()
 * 3. Set slip parameters for each vertex via setSlipParams()
 * 4. Call initialize() before time stepping
 * 5. Call prestep()/computeResidual()/computeJacobian()/poststep() each step
 */
class FaultCohesiveKin : public FaultCohesive {
public:
    FaultCohesiveKin();
    
    void configure(const std::map<std::string, std::string>& config) override;
    
    // Set slip time function
    void setSlipTimeFn(std::unique_ptr<SlipTimeFn> fn);
    void setSlipTimeFnType(SlipTimeFnType type);
    
    // Set slip parameters at each vertex
    void setSlipParams(const std::vector<KinematicSlipParams>& params);
    void setUniformSlip(double slip_ll, double slip_reverse, double slip_opening,
                       double slip_time, double rise_time);
    
    // Set slip from rupture propagation (circular crack)
    void setCircularRupture(double hypocenter_x, double hypocenter_y, double hypocenter_z,
                           double rupture_velocity, double final_slip,
                           double rise_time, double rake_angle);
    
    // Time stepping interface
    void initialize() override;
    void prestep(double t, double dt) override;
    void computeResidual(const double* solution, double* residual) override;
    void computeJacobian(const double* solution, double* jacobian) override;
    void poststep(double t, double dt) override;
    
    // Get prescribed slip at given time
    void getSlipAtTime(double t, size_t vertex_idx, 
                      double& slip_ll, double& slip_reverse, double& slip_opening) const;
    void getSlipRateAtTime(double t, size_t vertex_idx,
                          double& rate_ll, double& rate_reverse, double& rate_opening) const;
    
    // Access kinematic parameters
    const KinematicSlipParams& getSlipParams(size_t i) const { return kin_params[i]; }
    
private:
    std::unique_ptr<SlipTimeFn> slip_time_fn;
    std::vector<KinematicSlipParams> kin_params;
    
    // Current time for slip computation
    double current_time;
};

// =============================================================================
// Dynamic Fault (FaultCohesiveDyn)
// =============================================================================

/**
 * @brief Dynamic fault with spontaneous rupture
 * 
 * Implements PyLith's FaultCohesiveDyn:
 * - Rupture propagates spontaneously based on friction law
 * - Supports static, slip-weakening, and rate-state friction
 * - Traction perturbations can be applied to trigger rupture
 * 
 * The fault is governed by the fault constitutive relation:
 * τ = f(V, θ) × σ_n^eff
 * 
 * where τ is shear traction, f is the friction coefficient,
 * V is slip rate, θ is the state variable, and σ_n^eff is
 * effective normal stress.
 */
class FaultCohesiveDyn : public FaultCohesive {
public:
    FaultCohesiveDyn();
    
    void configure(const std::map<std::string, std::string>& config) override;
    
    // Set friction model
    void setFrictionModel(std::unique_ptr<FaultFrictionModel> model);
    void setFrictionModelType(const std::string& type);
    
    // Initial stress state
    void setInitialTraction(const FaultTractionField& initial_traction);
    void setUniformInitialTraction(double tau_ll, double tau_ud, double sigma_n);
    void setInitialState(const std::vector<DynamicFaultState>& states);
    
    // Traction perturbation to trigger rupture
    void setTractionPerturbation(TractionPerturbationType type,
                                double amplitude, double duration,
                                double center_x, double center_y, double center_z,
                                double width);
    void setTractionPerturbationFunction(
        std::function<void(double t, double x, double y, double z, double& dtau)> fn);
    
    // Time stepping interface
    void initialize() override;
    void prestep(double t, double dt) override;
    void computeResidual(const double* solution, double* residual) override;
    void computeJacobian(const double* solution, double* jacobian) override;
    void poststep(double t, double dt) override;
    
    // Access dynamic state
    const DynamicFaultState& getState(size_t i) const { return dynamic_states[i]; }
    
    // Rupture analysis
    bool hasRuptured() const;
    double getRuptureTime() const;
    double getMaxSlipRate() const;
    double getMaxSlip() const;
    
    // Energy computations
    double getFrictionalWork() const;
    double getFractureEnergy() const;
    double getRadiatedEnergy() const;
    
private:
    std::unique_ptr<FaultFrictionModel> friction_model;
    std::vector<DynamicFaultState> dynamic_states;
    
    // Initial traction
    FaultTractionField initial_traction;
    
    // Traction perturbation
    TractionPerturbationType pert_type;
    double pert_amplitude;
    double pert_duration;
    double pert_center[3];
    double pert_width;
    std::function<void(double, double, double, double, double&)> pert_function;
    
    // Current time
    double current_time;
    
    // Energy tracking
    mutable double cumulative_frictional_work;
    
    // Helper functions
    void updateFrictionState(double dt);
    double computeTractionPerturbation(double t, double x, double y, double z) const;
    void solveFrictionEquation(size_t vertex_idx, double sigma_n_eff,
                              double trial_tau, double& tau, double& slip_rate);
};

// =============================================================================
// Factory Functions
// =============================================================================

/**
 * @brief Create slip time function from type
 */
std::unique_ptr<SlipTimeFn> createSlipTimeFn(SlipTimeFnType type);

/**
 * @brief Create slip time function from string
 */
std::unique_ptr<SlipTimeFn> createSlipTimeFn(const std::string& type_str);

/**
 * @brief Create friction model from string
 */
std::unique_ptr<FaultFrictionModel> createFaultFrictionModel(const std::string& type_str);

/**
 * @brief Create cohesive fault from configuration
 */
std::unique_ptr<FaultCohesive> createFaultCohesive(
    const std::map<std::string, std::string>& config);

/**
 * @brief Parse slip time function type from string
 */
SlipTimeFnType parseSlipTimeFnType(const std::string& str);

/**
 * @brief Parse cohesive fault type from string
 */
CohesiveFaultType parseCohesiveFaultType(const std::string& str);

// =============================================================================
// Fault System (Multiple Faults)
// =============================================================================

/**
 * @brief System of multiple interacting faults
 * 
 * Manages a collection of faults (kinematic or dynamic) and handles
 * their interaction through stress transfer.
 */
class FaultSystem {
public:
    FaultSystem();
    
    void configure(const std::vector<std::map<std::string, std::string>>& configs);
    
    // Add faults
    void addFault(std::unique_ptr<FaultCohesive> fault);
    
    // Time stepping (calls all faults)
    void initialize();
    void prestep(double t, double dt);
    void computeResidual(const double* solution, double* residual);
    void computeJacobian(const double* solution, double* jacobian);
    void poststep(double t, double dt);
    
    // Access
    size_t numFaults() const { return faults.size(); }
    FaultCohesive* getFault(size_t i) { return faults[i].get(); }
    FaultCohesive* getFault(const std::string& name);
    
    // Output all faults
    void writeOutput(const std::string& base_filename, double time);
    
    // Enable/disable stress interaction
    void enableStressInteraction(bool enable) { compute_interaction = enable; }
    
private:
    std::vector<std::unique_ptr<FaultCohesive>> faults;
    bool compute_interaction;
};

// =============================================================================
// Utility Functions
// =============================================================================

namespace FaultUtils {
    
/**
 * @brief Compute moment magnitude from slip and area
 */
double computeMomentMagnitude(double slip, double area, double shear_modulus);

/**
 * @brief Compute seismic moment
 */
double computeSeismicMoment(double slip, double area, double shear_modulus);

/**
 * @brief Compute stress drop for circular rupture
 */
double computeStressDrop(double slip, double radius, double shear_modulus);

/**
 * @brief Compute corner frequency from rupture area and stress drop
 */
double computeCornerFrequency(double moment, double stress_drop, 
                             double shear_velocity, double density);

/**
 * @brief Compute nucleation length for rate-state friction
 */
double computeNucleationLength(double a, double b, double dc, 
                              double sigma_n, double shear_modulus);

/**
 * @brief Compute critical stiffness for spring-slider stability
 */
double computeCriticalStiffness(double b_minus_a, double sigma_n, double dc);

/**
 * @brief Estimate rise time from slip and slip rate
 */
double estimateRiseTime(double slip, double slip_rate);

/**
 * @brief Compute rupture velocity for Mode II crack
 */
double computeRuptureVelocity(double shear_wave_velocity, 
                             double rayleigh_velocity,
                             bool supershear = false);

} // namespace FaultUtils

} // namespace FSRM

#endif // PYLITH_FAULT_HPP
