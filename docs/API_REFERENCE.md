# FSRM C++ API Reference

Reference documentation for the FSRM C++ programming interface.

---

## Table of Contents

1. [Core Classes](#core-classes)
2. [Physics Kernels](#physics-kernels)
3. [Fluid Models](#fluid-models)
4. [Material Models](#material-models)
5. [Fault Models](#fault-models)
6. [Configuration](#configuration)
7. [I/O Classes](#io-classes)
8. [Utilities](#utilities)

---

## Core Classes

### `FSRM::Simulator`

Main simulation controller.

```cpp
#include "Simulator.hpp"

namespace FSRM {

class Simulator {
public:
    // Construction
    Simulator();
    ~Simulator();
    
    // Initialization
    PetscErrorCode initializeFromConfigFile(const std::string& filename);
    PetscErrorCode initialize(const SimulationConfig& config,
                              const GridConfig& grid);
    
    // Setup
    PetscErrorCode setupDM();
    PetscErrorCode setupFields();
    PetscErrorCode setupPhysics();
    PetscErrorCode setupTS();
    
    // Configuration
    void setMaterialProperties(const MaterialProperties& props);
    void setFluidProperties(const FluidProperties& props);
    void addWell(const WellConfig& well);
    void addFracture(const FractureConfig& frac);
    void addFault(const FaultConfig& fault);
    void setBoundaryCondition(const BCConfig& bc);
    void setInitialCondition(const ICConfig& ic);
    
    // Execution
    PetscErrorCode run();
    PetscErrorCode step();
    PetscErrorCode finalize();
    
    // Accessors
    DM getDM() const;
    Vec getSolution() const;
    TS getTS() const;
    double getCurrentTime() const;
    int getCurrentStep() const;
    
    // Output
    PetscErrorCode writeOutput(const std::string& filename);
    PetscErrorCode writeCheckpoint(const std::string& filename);
    PetscErrorCode loadCheckpoint(const std::string& filename);
};

}
```

**Example:**

```cpp
FSRM::Simulator sim;
sim.initializeFromConfigFile("simulation.config");
sim.run();
sim.finalize();
```

---

### Configuration Structures

```cpp
#include "ReservoirSim.hpp"

namespace FSRM {

// Simulation settings
struct SimulationConfig {
    std::string name;
    double startTime, endTime;
    double dtInitial, dtMin, dtMax;
    FluidModelType fluidModel;
    SolidModelType solidModel;
    bool enableGeomechanics;
    bool enableThermal;
    bool enableFractures;
    bool enableFaults;
    bool enableElastodynamics;
    double rtol, atol;
    int maxNonlinearIter;
    bool useGPU;
    GPUMode gpuMode;
};

// Grid configuration
struct GridConfig {
    int nx, ny, nz;
    double Lx, Ly, Lz;
    double originX, originY, originZ;
    GridType gridType;
    int refinementLevel;
};

// Well configuration
struct WellConfig {
    std::string name;
    WellType type;
    int i, j, k;                    // Grid location
    ControlMode controlMode;
    double targetValue;
    double maxRate, minBHP, maxBHP;
    double diameter;
    double skin;
};

// Fracture configuration
struct FractureConfig {
    FractureType type;
    std::vector<double> location;   // x,y,z,nx,ny,nz
    double aperture;
    double permeability;
    bool enablePropagation;
    double toughness;
};

// Fault configuration
struct FaultConfig {
    std::string name;
    double x, y, z;                 // Center
    double strike, dip, rake;
    double length, width;
    FrictionLaw frictionLaw;
    double staticFriction;
    double dynamicFriction;
    // Rate-state parameters
    double a, b, Dc, V0, f0;
};

// Boundary condition
struct BCConfig {
    BCType type;
    FieldType field;
    BCLocation location;
    double value;
};

// Initial condition
struct ICConfig {
    FieldType field;
    ICDistribution distribution;
    double value;
    std::vector<double> gradient;
};

}
```

---

## Physics Kernels

### `FSRM::PhysicsKernel`

Abstract base class for physics implementations.

```cpp
#include "PhysicsKernel.hpp"

namespace FSRM {

class PhysicsKernel {
public:
    virtual ~PhysicsKernel() = default;
    
    // Core methods (pure virtual)
    virtual PetscErrorCode residual(
        DM dm, Vec X, Vec R, void* ctx) = 0;
    
    virtual PetscErrorCode jacobian(
        DM dm, Vec X, Mat J, Mat P, void* ctx) = 0;
    
    // Physics type
    virtual PhysicsType getType() const = 0;
    
    // Configuration
    void setMaterialProperties(const MaterialProperties& props);
    void setFluidProperties(const FluidProperties& props);
    void setTimeStep(double dt);
};

// Factory function
std::unique_ptr<PhysicsKernel> createPhysicsKernel(PhysicsType type);

}
```

### Available Kernels

| Class | Physics |
|-------|---------|
| `SinglePhaseFlowKernel` | Single-phase compressible flow |
| `BlackOilKernel` | Three-phase black oil |
| `CompositionalKernel` | Multi-component with EOS |
| `GeomechanicsKernel` | Linear/nonlinear mechanics |
| `ThermalKernel` | Heat conduction/convection |
| `CoupledFlowGeomechKernel` | Poroelastic coupling |
| `FractureKernel` | DFN and hydraulic fracturing |
| `ElastodynamicsKernel` | Elastic wave propagation |
| `PoroelastodynamicsKernel` | Biot wave equations |

---

## Fluid Models

### `FSRM::FluidModelBase`

```cpp
#include "FluidModel.hpp"

namespace FSRM {

class FluidModelBase {
public:
    virtual ~FluidModelBase() = default;
    
    // Core property methods
    virtual double getDensity(double pressure, double temperature) const = 0;
    virtual double getViscosity(double pressure, double temperature) const = 0;
    virtual double getCompressibility(double pressure, double temperature) const = 0;
    
    // Configuration
    virtual void configure(const std::map<std::string, std::string>& config) = 0;
    
    // Type identification
    virtual FluidType getType() const = 0;
};

// Concrete implementations
class SinglePhaseFluid : public FluidModelBase { /* ... */ };
class BlackOilFluid : public FluidModelBase { /* ... */ };
class CompositionalFluid : public FluidModelBase { /* ... */ };
class BrineFluid : public FluidModelBase { /* ... */ };
class CO2Fluid : public FluidModelBase { /* ... */ };

// Factory function
std::unique_ptr<FluidModelBase> createFluidModel(
    const std::string& type,
    const std::map<std::string, std::string>& config);

}
```

### `FSRM::BlackOilFluid`

```cpp
class BlackOilFluid : public FluidModelBase {
public:
    // Standard API
    double getDensity(double p, double T) const override;
    double getViscosity(double p, double T) const override;
    double getCompressibility(double p, double T) const override;
    
    // Black-oil specific
    double getSolutionGOR(double pressure) const;
    double getOilFVF(double pressure, double Rs) const;
    double getGasFVF(double pressure, double temperature) const;
    double getBubblePoint() const;
    
    // Phase properties
    double getOilDensity(double pressure) const;
    double getGasDensity(double pressure, double temperature) const;
    double getWaterDensity(double pressure) const;
    
    // Configuration
    void setPVTCorrelation(PVTCorrelation corr);
    void setStockTankProperties(double oilAPI, double gasGravity);
};
```

### `FSRM::CompositionalFluid`

```cpp
class CompositionalFluid : public FluidModelBase {
public:
    // Component setup
    void setComponents(const std::vector<std::string>& names,
                       const std::vector<double>& Tc,
                       const std::vector<double>& Pc,
                       const std::vector<double>& omega,
                       const std::vector<double>& MW);
    
    // EOS
    void setEOS(EOSType eos);
    
    // Flash calculation
    FlashResult flash(double pressure, double temperature,
                      const std::vector<double>& z) const;
    
    // Fugacity
    std::vector<double> getFugacity(double P, double T,
                                    const std::vector<double>& x,
                                    Phase phase) const;
};

struct FlashResult {
    double V;                       // Vapor fraction
    std::vector<double> x;          // Liquid composition
    std::vector<double> y;          // Vapor composition
    double liquidDensity;
    double vaporDensity;
};
```

---

## Material Models

### `FSRM::MaterialModelBase`

```cpp
#include "MaterialModel.hpp"

namespace FSRM {

class MaterialModelBase {
public:
    virtual ~MaterialModelBase() = default;
    
    // Stress computation
    virtual StressTensor computeStress(
        const StrainTensor& strain,
        const StressTensor* prevStress = nullptr,
        double dt = 0.0) const = 0;
    
    // Stiffness matrix
    virtual std::array<double, 36> getStiffnessMatrix() const = 0;
    
    // Properties
    virtual double getDensity() const = 0;
    virtual double getPorosity() const = 0;
    virtual double getYoungsModulus() const = 0;
    virtual double getPoissonsRatio() const = 0;
    
    // Configuration
    virtual void configure(const std::map<std::string, std::string>& config) = 0;
};

// Implementations
class LinearElasticMaterial : public MaterialModelBase { /* ... */ };
class ViscoelasticMaterial : public MaterialModelBase { /* ... */ };
class PoroelasticMaterial : public MaterialModelBase { /* ... */ };
class ElastoplasticMaterial : public MaterialModelBase { /* ... */ };
class AnisotropicMaterial : public MaterialModelBase { /* ... */ };

// Factory
std::unique_ptr<MaterialModelBase> createMaterialModel(
    const std::string& type,
    const std::map<std::string, std::string>& config);

}
```

### Stress and Strain Tensors

```cpp
struct StressTensor {
    double xx, yy, zz;      // Normal components
    double xy, xz, yz;      // Shear components
    
    double meanStress() const;
    double vonMises() const;
    double maxPrincipal() const;
    double minPrincipal() const;
    std::array<double, 3> principals() const;
};

struct StrainTensor {
    double xx, yy, zz;
    double xy, xz, yz;
    
    double volumetric() const;
    double deviatoric() const;
};
```

### `FSRM::PoroelasticMaterial`

```cpp
class PoroelasticMaterial : public MaterialModelBase {
public:
    // Poroelastic-specific
    double getBiotCoefficient() const;
    double getBiotModulus() const;
    double getUndrainedPoissonsRatio() const;
    double getSkemptonCoefficient() const;
    
    // Effective stress
    StressTensor getEffectiveStress(
        const StressTensor& total, double porePressure) const;
    
    // Storage coefficient
    double getStorageCoefficient() const;
};
```

### `FSRM::RockProperties`

Aggregates all rock property models.

```cpp
class RockProperties {
public:
    // Access sub-models
    MaterialModelBase* getMechanicalModel();
    PermeabilityModel* getPermeabilityModel();
    ThermalProperties* getThermalProperties();
    FractureProperties* getFractureProperties();
    
    // Quick access
    double getPorosity() const;
    double getPermeability(double pressure, double strain) const;
    double getThermalConductivity() const;
};
```

---

## Fault Models

### `FSRM::FrictionModelBase`

```cpp
#include "FaultModel.hpp"

namespace FSRM {

class FrictionModelBase {
public:
    virtual ~FrictionModelBase() = default;
    
    virtual double getFriction(
        double slipVelocity,
        double stateVariable,
        double normalStress) const = 0;
    
    virtual double getStateEvolutionRate(
        double slipVelocity,
        double stateVariable) const = 0;
    
    virtual void configure(const std::map<std::string, std::string>& config) = 0;
};

class CoulombFriction : public FrictionModelBase { /* ... */ };
class RateStateFriction : public FrictionModelBase { /* ... */ };

// Factory
std::unique_ptr<FrictionModelBase> createFrictionModel(
    FrictionLaw law,
    const std::map<std::string, std::string>& config);

}
```

### `FSRM::FaultModel`

```cpp
class FaultModel {
public:
    FaultModel(const FaultGeometry& geom,
               std::unique_ptr<FrictionModelBase> friction);
    
    // Geometry
    const FaultGeometry& getGeometry() const;
    
    // Stress analysis
    FaultStressState resolveStress(const StressTensor& stress) const;
    double getCoulombStress(const FaultStressState& state,
                            double porePressure) const;
    
    // Slip update
    void updateSlip(const FaultStressState& state,
                    double porePressure, double dt);
    
    // Current state
    double getSlip() const;
    double getSlipRate() const;
    double getStateVariable() const;
    
    // Event detection
    bool isSlipping() const;
    std::optional<SeismicEvent> checkForEvent();
    const std::vector<SeismicEvent>& getEventHistory() const;
};
```

### `FSRM::FaultNetwork`

```cpp
class FaultNetwork {
public:
    void addFault(std::unique_ptr<FaultModel> fault);
    void addFault(const FaultConfig& config);
    
    // Update all faults
    void update(const StressTensor& regionalStress,
                double porePressure, double dt);
    
    // Stress transfer
    void computeStressTransfer();
    
    // Access
    size_t size() const;
    FaultModel* getFault(size_t index);
    FaultModel* getFaultByName(const std::string& name);
    
    // Events
    std::vector<SeismicEvent> getAllEvents() const;
    void writeSeismicCatalog(const std::string& filename) const;
};
```

### Seismicity Analysis

```cpp
namespace SeismicityAnalysis {
    // Magnitude calculations
    double momentToMagnitude(double moment);
    double magnitudeToMoment(double Mw);
    double computeMoment(double shearModulus, double area, double slip);
    
    // Statistics
    double estimateBValue(const std::vector<SeismicEvent>& events);
    double gutenbergRichterRate(double a, double b, double M);
    
    // Aftershocks
    double omoriRate(double K, double c, double p, double time);
}
```

---

## Configuration

### `FSRM::ConfigReader`

```cpp
#include "ConfigReader.hpp"

namespace FSRM {

class ConfigReader {
public:
    ConfigReader();
    
    // Loading
    bool loadFile(const std::string& filename);
    void mergeFile(const std::string& filename);
    
    // Parsing main structures
    SimulationConfig parseSimulation();
    GridConfig parseGrid();
    MaterialProperties parseMaterial();
    FluidProperties parseFluid();
    
    // Parsing collections
    std::vector<WellConfig> parseWells();
    std::vector<FractureConfig> parseFractures();
    std::vector<FaultConfig> parseFaults();
    std::vector<BCConfig> parseBCs();
    std::vector<ICConfig> parseICs();
    
    // New generic model parsing
    std::unique_ptr<FluidModelBase> parseFluidModel();
    std::vector<RockProperties> parseRockProperties();
    FaultNetwork parseFaultNetwork();
    
    // Validation
    bool validate() const;
    std::vector<std::string> getErrors() const;
    
    // Raw access
    std::string getValue(const std::string& section,
                         const std::string& key,
                         const std::string& defaultVal = "") const;
    
    // Template generation
    static void generateCompleteTemplate(const std::string& filename);
};

}
```

---

## I/O Classes

### `FSRM::EclipseIO`

```cpp
#include "EclipseIO.hpp"

namespace FSRM {

class EclipseIO {
public:
    // Reading
    bool readDATAFile(const std::string& filename);
    GridConfig getGrid() const;
    std::vector<WellConfig> getWells() const;
    
    // Writing
    bool writeRestartFile(const std::string& filename,
                          Vec solution, double time);
    bool writeUnrstFile(const std::string& filename,
                        const std::vector<Vec>& solutions,
                        const std::vector<double>& times);
};

}
```

### `FSRM::Visualization`

```cpp
#include "Visualization.hpp"

namespace FSRM {

class Visualization {
public:
    // VTK output
    PetscErrorCode writeVTK(DM dm, Vec solution,
                            const std::string& filename);
    
    // Field selection
    void enableField(const std::string& fieldName);
    void disableField(const std::string& fieldName);
    
    // Time series
    void writeTimeSeries(DM dm, const std::vector<Vec>& solutions,
                         const std::vector<double>& times,
                         const std::string& basename);
};

}
```

---

## Utilities

### GPU Management

```cpp
#include "GPUManager.hpp"

namespace FSRM {

class GPUManager {
public:
    static GPUManager& instance();
    
    bool initialize(int deviceId = 0);
    bool isAvailable() const;
    void setMemoryFraction(double fraction);
    
    // Memory management
    void* allocate(size_t bytes);
    void free(void* ptr);
    void copyToDevice(void* dst, const void* src, size_t bytes);
    void copyToHost(void* dst, const void* src, size_t bytes);
};

}
```

### Testing Framework

```cpp
#include "Testing.hpp"

namespace FSRM {

class TestFramework {
public:
    // Run tests
    static int runAllTests(int argc, char** argv);
    static int runTest(const std::string& testName);
    
    // Verification
    static bool verifyConvergence(Vec numerical, Vec analytical,
                                  double tolerance);
    static double computeL2Error(Vec numerical, Vec analytical);
};

}
```

---

## Error Handling

FSRM uses PETSc error codes. Check returns with:

```cpp
PetscErrorCode ierr;
ierr = sim.initialize(config, grid); CHKERRQ(ierr);
ierr = sim.run(); CHKERRQ(ierr);
```

Common error codes:
- `0` = success
- Non-zero = error (use PETSc error macros)

---

## Thread Safety

- `Simulator` instances are NOT thread-safe
- `FluidModelBase` implementations are thread-safe for const methods
- `MaterialModelBase` implementations are thread-safe for const methods
- Use MPI for parallelism, not threads

---

## Memory Management

- Use smart pointers (`std::unique_ptr`, `std::shared_ptr`)
- PETSc objects (Vec, Mat, DM) managed via PETSc destroy functions
- Factory functions return `std::unique_ptr` for ownership transfer

```cpp
auto fluid = createFluidModel("BLACK_OIL", config);  // unique_ptr
// fluid automatically destroyed when out of scope
```
