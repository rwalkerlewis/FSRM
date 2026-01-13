# FSRM API Reference

## Overview

This document provides a comprehensive API reference for the FSRM (Fully-coupled Seismic Reservoir Model) simulator. FSRM combines reservoir engineering with earthquake simulation capabilities.

## Table of Contents

1. [Core Modules](#core-modules)
2. [Discontinuous Galerkin Solver](#discontinuous-galerkin-solver)
3. [Seismic Sources and Receivers](#seismic-sources-and-receivers)
4. [Boundary Conditions](#boundary-conditions)
5. [Friction Laws](#friction-laws)
6. [Plasticity Models](#plasticity-models)
7. [Viscoelastic Attenuation](#viscoelastic-attenuation)
8. [Performance Optimization](#performance-optimization)
9. [Configuration System](#configuration-system)

---

## Core Modules

### Simulator

The main simulation controller.

```cpp
#include "Simulator.hpp"

class Simulator {
public:
    // Initialization
    void initialize(const std::string& config_file);
    void initialize(const SimulationConfig& config);
    
    // Time stepping
    void step(double dt);
    void run();
    
    // Queries
    double getCurrentTime() const;
    double getTimeStep() const;
    int getStep() const;
    
    // Physics modules
    void setFaultModel(FaultModel& fault);
    void setMaterialModel(MaterialModel& material);
    void setBoundaryConditions(BoundaryConditionManager& bc_mgr);
    void setSourceReceiverManager(SourceReceiverManager& srm);
    
    // Output
    void writeOutput(const std::string& filename);
    void checkpoint(const std::string& filename);
    void restart(const std::string& filename);
};
```

### SimulationConfig

Configuration structure for simulations.

```cpp
struct SimulationConfig {
    // Domain
    double domain_xmin, domain_xmax;
    double domain_ymin, domain_ymax;
    double domain_zmin, domain_zmax;
    
    // Discretization
    double element_size;
    int dg_order;           // Polynomial order (1-10)
    double cfl;             // CFL number
    
    // Time
    double t_start;
    double t_end;
    double dt_max;
    
    // Physics
    bool enable_plasticity;
    bool enable_attenuation;
    bool use_local_time_stepping;
    
    // Output
    std::string output_directory;
    int output_interval;
    OutputFormat output_format;
};
```

---

## Discontinuous Galerkin Solver

### DGSolver

High-order Discontinuous Galerkin solver with ADER time integration.

```cpp
#include "DiscontinuousGalerkin.hpp"

class DGSolver {
public:
    // Setup
    void initialize(const SimulationConfig& config, const MaterialModel& material);
    void setMesh(const Mesh& mesh);
    void setBoundaryConditions(BoundaryConditionManager& bc_mgr);
    
    // Time stepping
    double computeTimeStep();
    void step(double dt);
    
    // ADER predictor
    void computeADERPredictor(int order);
    
    // Local time stepping
    void enableLTS(int num_clusters);
    void disableLTS();
    
    // Queries
    void getSolution(int element, double* u) const;
    void setSolution(int element, const double* u);
    
    // Flux methods
    void setRiemannSolver(FluxMethod method);
};

// Riemann solver options
enum class FluxMethod {
    RUSANOV,    // Rusanov (Local Lax-Friedrichs)
    ROE,        // Roe linearization
    HLL,        // Harten-Lax-van Leer
    HLLC,       // HLL with contact wave
    GODUNOV     // Exact Godunov
};
```

### QuadratureRule

Numerical quadrature rules.

```cpp
struct QuadratureRule {
    std::vector<std::array<double, 3>> points;  // Quadrature points
    std::vector<double> weights;                 // Weights
    
    static QuadratureRule gaussLegendre(int n, int dim);
    static QuadratureRule gaussLobatto(int n, int dim);
    static QuadratureRule dunavant(int order, ElementType type);
    static QuadratureRule grundmannMoeller(int order, ElementType type);
};
```

### BasisFunctions

Polynomial basis functions for DG.

```cpp
class BasisFunctions {
public:
    BasisFunctions(BasisType type, int order, ElementType element);
    
    void evaluate(double x, double y, double z, std::vector<double>& values) const;
    void evaluateGradients(double x, double y, double z, std::vector<double>& gradients) const;
    
    int getNumDOFs() const;
    int getOrder() const;
};

enum class BasisType {
    LAGRANGE,   // Lagrange interpolating polynomials
    LEGENDRE,   // Legendre orthogonal polynomials
    DUBINER,    // Dubiner basis (triangles/tetrahedra)
    MODAL,      // Modal (hierarchical) basis
    NODAL       // Nodal (interpolating) basis
};
```

### RiemannSolver

Numerical flux computation at element interfaces.

```cpp
class RiemannSolver {
public:
    RiemannSolver(FluxMethod method);
    
    void solve(const double* u_L, const double* u_R,
               const double* normal,
               double rho, double lambda, double mu,
               double* flux);
    
    void setMaxWaveSpeed(double c_max);
};
```

---

## Seismic Sources and Receivers

### MomentTensor

Seismic moment tensor representation.

```cpp
#include "SeismicSource.hpp"

struct MomentTensor {
    double Mxx, Myy, Mzz;    // Diagonal components
    double Mxy, Mxz, Myz;    // Off-diagonal components
    
    // Constructors
    MomentTensor();
    MomentTensor(double mxx, double myy, double mzz,
                 double mxy, double mxz, double myz);
    
    // Factory methods
    static MomentTensor doubleCouple(double strike, double dip, double rake, double M0);
    static MomentTensor explosion(double M0);
    static MomentTensor clvd(double M0, const std::array<double, 3>& axis);
    
    // Properties
    double scalarMoment() const;
    double magnitude() const;     // Moment magnitude Mw
    void decompose(MomentTensor& iso, MomentTensor& dc, MomentTensor& clvd) const;
    
    // Operations
    MomentTensor operator+(const MomentTensor& other) const;
    MomentTensor operator*(double scale) const;
    void toArray(double* arr) const;
    void toMatrix(double mat[3][3]) const;
};
```

### SourceTimeFunction

Source time function evaluation.

```cpp
enum class SourceTimeFunction {
    GAUSSIAN,       // Gaussian pulse
    RICKER,         // Ricker wavelet
    STEP,           // Step function
    RAMP,           // Linear ramp
    TRIANGLE,       // Triangular pulse
    SMOOTHED_RAMP,  // Smoothed ramp (erf-based)
    YOFFE,          // Yoffe regularized slip rate
    BRUNE,          // Brune's omega-squared spectrum
    CUSTOM          // User-defined function
};

class SourceTimeFunctionEvaluator {
public:
    SourceTimeFunctionEvaluator(SourceTimeFunction type);
    
    // Parameters
    void setDuration(double T);
    void setPeakTime(double t_peak);
    void setOnsetTime(double t_onset);
    void setFrequency(double f);      // For Ricker
    void setRiseTime(double t_rise);  // For Yoffe
    void setCustomFunction(std::function<double(double)> func);
    
    // Evaluation
    double momentRate(double t) const;  // M'(t) - moment rate
    double moment(double t) const;      // M(t) - cumulative moment
};
```

### PointSource

Point seismic source.

```cpp
class PointSource {
public:
    // Location
    void setLocation(double x, double y, double z);
    void getLocation(double& x, double& y, double& z) const;
    std::array<double, 3> getLocation() const;
    
    // Source mechanism
    void setMomentTensor(const MomentTensor& M);
    void setFromFault(double strike, double dip, double rake, double M0,
                      double x, double y, double z);
    void setSourceTimeFunction(const SourceTimeFunctionEvaluator& stf);
    
    // Spatial smoothing
    void setSpatialSmoothing(double sigma);
    double getSpatialWeight(double x, double y, double z) const;
    
    // Queries
    double getMomentRate(double t) const;
    void getSourceTerm(double x, double y, double z, double t, double* f) const;
};
```

### KinematicSource

Finite fault kinematic source.

```cpp
struct KinematicSubfault {
    double x, y, z;           // Center position
    double area;              // Subfault area
    double strike, dip;       // Orientation
    double slip;              // Total slip
    double slip_angle;        // Rake angle
    double rupture_time;      // Rupture onset time
    double rise_time;         // Slip duration
};

class KinematicSource {
public:
    // Fault geometry
    void createRectangularFault(double cx, double cy, double cz,
                                double strike, double dip,
                                double length, double width,
                                int n_strike, int n_dip);
    void loadFromSRF(const std::string& filename);
    void loadFromFSP(const std::string& filename);
    
    // Slip distribution
    void setUniformSlip(double slip, double rake);
    void setSlipFromFile(const std::string& filename);
    void setEllipticalAsperity(double cx, double cy, double cz,
                               double major, double minor, double max_slip);
    
    // Rupture propagation
    void setCircularRupture(double hx, double hy, double hz, double Vr);
    void setEllipticalRupture(double hx, double hy, double hz,
                              double Vr_strike, double Vr_dip);
    void setRuptureFrontFromFile(const std::string& filename);
    
    // Properties
    void setShearModulus(double mu);
    double getTotalMoment() const;
    double getMagnitude() const;
    
    // Access
    const std::vector<KinematicSubfault>& getSubfaults() const;
    std::vector<KinematicSubfault>& getSubfaults();
};
```

### SeismicReceiver

Surface or volume receiver.

```cpp
class SeismicReceiver {
public:
    // Setup
    void setLocation(double x, double y, double z);
    void setName(const std::string& name);
    void setSamplingRate(double dt);
    
    // Recording
    void record(double time, const double* velocity);
    void recordStress(double time, const double* stress);
    
    // Access
    const std::vector<double>& getTimes() const;
    const std::vector<std::array<double, 3>>& getData() const;
    
    // Analysis
    double getPeakGroundVelocity() const;
    double getPeakGroundAcceleration() const;
    std::array<double, 3> getPGVComponents() const;
    
    // Output
    void writeASCII(const std::string& filename) const;
    void writeSAC(const std::string& filename) const;
    void writeHDF5(const std::string& filename) const;
};
```

### FaultReceiver

On-fault receiver.

```cpp
class FaultReceiver {
public:
    void setLocation(double x, double y, double z);
    void setName(const std::string& name);
    void setSamplingRate(double dt);
    
    // Recording
    void record(double time, double slip, double slip_rate,
                double sigma_n, double tau, double pore_pressure,
                double state_variable);
    
    // Queries
    double getPeakSlipRate() const;
    double getFinalSlip() const;
    double getRuptureTime(double threshold) const;
    
    // Output
    void writeASCII(const std::string& filename) const;
};
```

### SourceReceiverManager

Manages multiple sources and receivers.

```cpp
class SourceReceiverManager {
public:
    // Sources
    void addPointSource(std::unique_ptr<PointSource> src);
    void addKinematicSource(std::unique_ptr<KinematicSource> src);
    void addSingleForceSource(std::unique_ptr<SingleForceSource> src);
    
    // Receivers
    void addReceiver(std::unique_ptr<SeismicReceiver> rec);
    void addFaultReceiver(std::unique_ptr<FaultReceiver> rec);
    void addReceiverLine(double x0, double y0, double z0,
                         double x1, double y1, double z1,
                         int n_receivers, const std::string& prefix);
    void addReceiverGrid(double xmin, double xmax, double ymin, double ymax,
                         double z, int nx, int ny, const std::string& prefix);
    
    // Evaluation
    void getSourceTerm(double x, double y, double z, double t, double* f) const;
    void recordAll(double t, const FieldData& solution);
    
    // Output
    void writeAllReceivers(const std::string& directory) const;
    
    // Access
    const std::vector<std::unique_ptr<PointSource>>& getPointSources() const;
    const std::vector<std::unique_ptr<SeismicReceiver>>& getReceivers() const;
};
```

---

## Boundary Conditions

### BoundaryConditionBase

Base class for boundary conditions.

```cpp
#include "BoundaryConditions.hpp"

struct BoundaryFace {
    int element_id;
    int local_face_id;
    std::array<double, 3> centroid;
    std::array<double, 3> normal;
    double area;
    double rho, vp, vs;
    double impedance_p, impedance_s;
};

class BoundaryConditionBase {
public:
    virtual ~BoundaryConditionBase() = default;
    
    virtual BoundaryType getType() const = 0;
    virtual void getBoundaryState(const BoundaryFace& face,
                                  const double* u_int,
                                  double* u_ext,
                                  double time) const = 0;
    virtual void computeBoundaryFlux(const BoundaryFace& face,
                                     const double* u_int,
                                     double* flux,
                                     double time) const = 0;
};

enum class BoundaryType {
    FREE_SURFACE,
    PML,
    CLAYTON_ENGQUIST,
    LYSMER,
    DIRICHLET,
    NEUMANN,
    PERIODIC,
    SYMMETRY
};
```

### FreeSurfaceBC

Traction-free surface boundary.

```cpp
class FreeSurfaceBC : public BoundaryConditionBase {
public:
    BoundaryType getType() const override { return BoundaryType::FREE_SURFACE; }
    
    void getBoundaryState(const BoundaryFace& face,
                          const double* u_int,
                          double* u_ext,
                          double time) const override;
    
    // Mirror method for stress
    void mirrorStress(const double* stress_int, const double* normal,
                      double* stress_ext) const;
    
    // Compute surface traction
    void getSurfaceTraction(const BoundaryFace& face,
                           const double* stress,
                           double* traction) const;
};
```

### PMLBoundary

Perfectly Matched Layer absorbing boundary.

```cpp
class PMLBoundary : public BoundaryConditionBase {
public:
    BoundaryType getType() const override { return BoundaryType::PML; }
    
    // Setup
    void setThickness(double d);
    void setOrder(int n);  // Polynomial order for damping profile
    void setReflectionCoefficient(double R);
    void setDomainBounds(double xmin, double xmax, double ymin, double ymax,
                         double zmin, double zmax);
    void setMaterialProperties(double rho, double vp, double vs);
    
    // Queries
    bool isInPML(double x, double y, double z) const;
    int getPMLDirection(double x, double y, double z) const;
    void getDampingCoefficients(double x, double y, double z,
                                double& d_x, double& d_y, double& d_z) const;
    
    // PML-specific operations
    void initializeAuxiliaryVariables(int n_elements);
    void evolveAuxiliaryVariables(double dt);
    void computePMLFlux(const BoundaryFace& face,
                        const double* u, double* flux) const;
};
```

### ClaytonEngquistBC

First/second order paraxial absorbing BC.

```cpp
class ClaytonEngquistBC : public BoundaryConditionBase {
public:
    BoundaryType getType() const override { return BoundaryType::CLAYTON_ENGQUIST; }
    
    void setOrder(int order);  // 1 or 2
    void setMaterialProperties(double rho, double vp, double vs);
    
    // Impedance matrix
    void computeImpedanceMatrix(const double* normal, double* Z) const;
    
    // Apply absorbing condition
    void applyFirstOrder(const BoundaryFace& face,
                         const double* velocity, double* traction) const;
    void applySecondOrder(const BoundaryFace& face,
                          const double* velocity, const double* acceleration,
                          double* traction) const;
};
```

### BoundaryConditionManager

Manages boundary conditions on a mesh.

```cpp
class BoundaryConditionManager {
public:
    // Setup
    void addBoundaryCondition(int boundary_tag, std::unique_ptr<BoundaryConditionBase> bc);
    void setDefaultBC(std::unique_ptr<BoundaryConditionBase> bc);
    void classifyBoundaryFaces(const Mesh& mesh);
    
    // PML setup
    void setupPML(double thickness, double reflection_coeff);
    void setPMLRegion(double xmin, double xmax, double ymin, double ymax,
                      double zmin, double zmax);
    
    // Application
    void applyBoundaryConditions(double time);
    void computeBoundaryFluxes(double time);
    void evolvePMLAuxiliaryVariables(double dt);
    
    // Queries
    BoundaryConditionBase* getBCForFace(int face_id);
};
```

---

## Friction Laws

### FrictionModel

Base class for friction laws.

```cpp
#include "FaultModel.hpp"

enum class FrictionLaw {
    COULOMB,
    RATE_STATE_AGING,
    RATE_STATE_SLIP,
    SLIP_WEAKENING,
    FLASH_HEATING,
    THERMAL_PRESSURIZATION,
    STRONG_VELOCITY_WEAKENING
};

class FrictionModel {
public:
    virtual ~FrictionModel() = default;
    
    virtual FrictionLaw getType() const = 0;
    virtual double getFriction(double sigma_n, double slip, double slip_rate) const = 0;
    virtual double getStateEvolutionRate(double state, double slip_rate) const = 0;
    virtual void configure(const std::map<std::string, double>& params) = 0;
};

// Factory function
std::unique_ptr<FrictionModel> createFrictionModel(FrictionLaw type);

// Parsing
FrictionLaw parseFrictionLaw(const std::string& name);
```

### RateStateFriction

Rate-and-state friction with aging or slip law.

```cpp
class RateStateFriction : public FrictionModel {
public:
    RateStateFriction(RateStateLaw law = RateStateLaw::AGING);
    
    // Parameters
    void setReferenceFriction(double f0);
    void setDirectEffect_a(double a);
    void setEvolutionEffect_b(double b);
    void setCriticalSlipDistance(double Dc);
    void setReferenceVelocity(double V0);
    
    // Implementation
    double getFriction(double sigma_n, double slip, double slip_rate) const override;
    double getStateEvolutionRate(double state, double slip_rate) const override;
    
    // Steady-state properties
    double getSteadyStateFriction(double sigma_n, double V) const;
    double getSteadyStateTheta(double V) const;
    
    // Stability analysis
    bool isVelocityWeakening() const;
    double getNucleationLength(double sigma_n) const;
};

enum class RateStateLaw {
    AGING,    // dθ/dt = 1 - θV/Dc
    SLIP      // dθ/dt = -θV/Dc * ln(θV/Dc)
};
```

### FlashHeatingFriction

Flash heating at high slip rates.

```cpp
class FlashHeatingFriction : public FrictionModel {
public:
    // Parameters
    void setLowVelocityFriction(double f_lv);
    void setContactDiameter(double D);
    void setWeakeningTemperature(double T_w);
    void setThermalDiffusivity(double alpha);
    void setHeatCapacity(double c);
    void setDensity(double rho);
    void setAmbientTemperature(double T0);
    
    // Derived quantities
    double computeWeakeningVelocity(double sigma_n) const;
    double getFlashTemperature(double sigma_n, double V) const;
    
    // Implementation
    double getFriction(double sigma_n, double slip, double slip_rate) const override;
};
```

### ThermalPressurizationFriction

Thermal pressurization model.

```cpp
class ThermalPressurizationFriction : public FrictionModel {
public:
    // Parameters
    void setBaseFriction(double f0);
    void setHydraulicDiffusivity(double alpha_hy);
    void setThermalDiffusivity(double alpha_th);
    void setSlipZoneWidth(double w);
    void setPressurization(double Lambda);  // Δp/ΔT
    void setHeatCapacity(double c);
    void setDensity(double rho);
    
    // State evolution
    void evolveTPState(double tau, double V, double dt,
                       double& delta_p, double& delta_T);
    double getPorePressureChange() const;
    double getTemperatureRise() const;
    void resetTPState();
    
    // Properties
    double getCharacteristicSlip() const;
    
    // Implementation
    double getFriction(double sigma_n, double slip, double slip_rate) const override;
};
```

### StrongVelocityWeakeningFriction

Strong velocity weakening model.

```cpp
class StrongVelocityWeakeningFriction : public FrictionModel {
public:
    // Parameters
    void setStaticFriction(double f_s);
    void setDynamicFriction(double f_d);
    void setCriticalVelocity(double V_w);
    void setWeakeningRate(double L);
    
    // Implementation
    double getFriction(double sigma_n, double slip, double slip_rate) const override;
    
    // Properties
    double getSlipWeakeningDistance(double sigma_n) const;
    double getFractureEnergy(double sigma_n) const;
    double getNucleationLength(double sigma_n, double mu) const;
};
```

---

## Plasticity Models

### PlasticityModel

Base class for plasticity models.

```cpp
#include "PlasticityModel.hpp"

class PlasticityModel {
public:
    virtual ~PlasticityModel() = default;
    
    // Yield function
    virtual double yieldFunction(const double* sigma) const = 0;
    
    // Return mapping
    virtual bool returnMap(const double* sigma_trial,
                          double* sigma_return,
                          double& plastic_multiplier) const = 0;
    
    // Consistent tangent
    virtual void computeConsistentTangent(const double* sigma,
                                          const double* eps_p,
                                          double* C_tan) const = 0;
};
```

### DruckerPragerModel

Pressure-dependent plasticity.

```cpp
class DruckerPragerModel : public PlasticityModel {
public:
    void setFrictionAngle(double phi);
    void setCohesion(double c);
    void setDilatancyAngle(double psi);
    void setElasticModuli(double E, double nu);
    
    // Yield function: f = sqrt(J2) + α*I1 - k
    double yieldFunction(const double* sigma) const override;
    
    // Return mapping with full algorithmic treatment
    bool returnMap(const double* sigma_trial,
                  double* sigma_return,
                  double& plastic_multiplier) const override;
};
```

### VonMisesModel

Isochoric (pressure-independent) plasticity.

```cpp
class VonMisesModel : public PlasticityModel {
public:
    void setYieldStress(double sigma_y);
    void setIsotropicHardening(double H);
    
    // Yield function: f = sqrt(3*J2) - σ_y
    double yieldFunction(const double* sigma) const override;
    
    // Radial return mapping
    bool returnMap(const double* sigma_trial,
                  double* sigma_return,
                  double& plastic_multiplier) const override;
    
    // Hardening
    void updateHardening(double delta_eps_p);
    double getCurrentYieldStress() const;
};
```

### MohrCoulombModel

Mohr-Coulomb plasticity with corners.

```cpp
class MohrCoulombModel : public PlasticityModel {
public:
    void setFrictionAngle(double phi);
    void setCohesion(double c);
    void setDilatancyAngle(double psi);
    
    double yieldFunction(const double* sigma) const override;
    
    // Return mapping handles singular apex and edges
    bool returnMap(const double* sigma_trial,
                  double* sigma_return,
                  double& plastic_multiplier) const override;
};
```

### CapModel

Cap model with shear yield and compaction cap.

```cpp
class CapModel : public PlasticityModel {
public:
    void setFrictionAngle(double phi);
    void setCohesion(double c);
    void setCapHardening(double H_cap);
    void setInitialCapPosition(double X_0);
    
    // Shear + cap yield
    double yieldFunction(const double* sigma) const override;
    
    // Return mapping for both surfaces
    bool returnMap(const double* sigma_trial,
                  double* sigma_return,
                  double* eps_plastic,
                  double& plastic_multiplier) const override;
    
    double getCurrentCapPosition() const;
};
```

### DamageModel

Continuum damage mechanics.

```cpp
class DamageModel {
public:
    void setCriticalStrain(double eps_c);
    void setDamageRate(double k);
    
    void updateDamage(double& damage, double plastic_work, double strain_eq);
    void applyDamage(const double* C_undamaged, double* C_damaged, double damage);
};
```

---

## Viscoelastic Attenuation

### ViscoelasticAttenuation

Multi-mechanism Q-factor attenuation.

```cpp
#include "ViscoelasticAttenuation.hpp"

class ViscoelasticAttenuation {
public:
    // Setup
    void setReferenceModulus(double M);
    void setReferenceDensity(double rho);
    void setTargetQ(double Q, double f_min, double f_max);
    void setFrequencyDependentQ(double Q0, double f0, double alpha);
    void setNumMechanisms(int n);
    void setFrequencyRange(double f_min, double f_max);
    
    // Compute coefficients via optimization
    void computeCoefficients();
    
    // State evolution
    void evolveAnelasticState(const double* strain_rate,
                              double* anelastic_state, double dt);
    void computeAnelasticStress(const double* anelastic_state,
                                double* sigma_anelastic) const;
    
    // Queries
    int getNumMechanisms() const;
    double getQualityFactor(double omega) const;
    double getVelocity(double omega) const;
    double getRelaxedModulus() const;
    
    // Validation
    bool checkCausality() const;
    double getQError(double f_min, double f_max) const;
};
```

### AnelasticIntegrator

Time integration for anelastic state.

```cpp
enum class IntegrationScheme {
    FORWARD_EULER,
    RK2,
    RK4,
    ADER
};

class AnelasticIntegrator {
public:
    AnelasticIntegrator(IntegrationScheme scheme);
    
    void initialize(const ViscoelasticAttenuation& model);
    void setADEROrder(int order);
    
    void step(double* state, const double* strain_rate, double dt);
    void stepADER(double* state,
                  std::function<void(double, double*)> strain_rate_func,
                  double t, double dt);
};
```

### AttenuationDatabase

Pre-configured attenuation values.

```cpp
class AttenuationDatabase {
public:
    double getQForRockType(const std::string& type) const;
    double getQAtDepth(double z) const;
    
    // Empirical relations
    double getQFromVelocity(double vs) const;
    double getQFromDamage(double Q0, double damage) const;
};
```

---

## Performance Optimization

### SIMD Utilities

```cpp
#include "PerformanceOptimizations.hpp"

// Aligned memory allocation
template<typename T>
using AlignedVector = std::vector<T, AlignedAllocator<T>>;

// SIMD vector class (4 doubles)
class SimdVec4d {
public:
    static SimdVec4d load(const double* ptr);
    static SimdVec4d loadu(const double* ptr);
    void store(double* ptr) const;
    
    SimdVec4d operator+(const SimdVec4d& b) const;
    SimdVec4d operator-(const SimdVec4d& b) const;
    SimdVec4d operator*(const SimdVec4d& b) const;
    SimdVec4d sqrt() const;
    double sum() const;
    
    static SimdVec4d fmadd(const SimdVec4d& a, const SimdVec4d& b, const SimdVec4d& c);
    static SimdVec4d max(const SimdVec4d& a, const SimdVec4d& b);
};
```

### Vectorized Kernels

```cpp
// Dot product
double dot_simd(const double* a, const double* b, size_t n);

// AXPY: y = a*x + y
void axpy_simd(double a, const double* x, double* y, size_t n);

// Matrix-vector product
void matvec_simd(const double* A, const double* x, double* y, int m, int n);

// Stress update
void stress_update_simd(const double* strain_rate, double* stress,
                        double lambda, double mu, double dt, int n_points);
```

### ThreadPool

```cpp
class ThreadPool {
public:
    ThreadPool(size_t num_threads = std::thread::hardware_concurrency());
    
    template<typename F>
    void enqueue(F&& f);
    
    void parallelFor(size_t start, size_t end,
                     std::function<void(size_t, size_t)> body);
    
    size_t numThreads() const;
};
```

### Profiler

```cpp
class Profiler {
public:
    void start(const std::string& name);
    void stop(const std::string& name);
    void report(std::ostream& os = std::cout) const;
    void reset();
};

extern Profiler g_profiler;
```

---

## Configuration System

### ConfigReader

```cpp
#include "ConfigReader.hpp"

class ConfigReader {
public:
    bool load(const std::string& filename);
    
    // Typed access
    template<typename T>
    T get(const std::string& section, const std::string& key) const;
    
    template<typename T>
    T get(const std::string& section, const std::string& key, T default_value) const;
    
    bool hasSection(const std::string& section) const;
    bool hasKey(const std::string& section, const std::string& key) const;
    
    std::vector<std::string> getSections() const;
    std::vector<std::string> getKeys(const std::string& section) const;
};
```

### Configuration Parsing

```cpp
// Parse complete simulation config
SimulationConfig parseSimulationConfig(const ConfigReader& reader);

// Parse individual components
MaterialModel parseMaterialConfig(const ConfigReader& reader, const std::string& section);
FaultModel parseFaultConfig(const ConfigReader& reader, const std::string& section);
WellModel parseWellConfig(const ConfigReader& reader, const std::string& section);
FractureModel parseFractureConfig(const ConfigReader& reader, const std::string& section);
BoundaryConditionConfig parseBCConfig(const ConfigReader& reader, const std::string& section);
```

---

## Error Handling

FSRM uses exceptions for error reporting:

```cpp
namespace FSRM {
    class FSRMException : public std::runtime_error {
    public:
        FSRMException(const std::string& msg);
    };
    
    class ConfigException : public FSRMException { /* ... */ };
    class MeshException : public FSRMException { /* ... */ };
    class SolverException : public FSRMException { /* ... */ };
    class IOException : public FSRMException { /* ... */ };
}
```

---

## Version Information

```cpp
namespace FSRM {
    constexpr int VERSION_MAJOR = 2;
    constexpr int VERSION_MINOR = 0;
    constexpr int VERSION_PATCH = 0;
    
    std::string getVersionString();
    bool hasGPUSupport();
    bool hasCUDASupport();
    bool hasROCmSupport();
    bool hasPROJSupport();
}
```
