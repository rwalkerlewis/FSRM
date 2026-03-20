/**
 * @file PhysicsModuleInterface.hpp
 * @brief Abstract interface for all physics modules in FSRM
 *
 * Every physics domain (reservoir flow, geomechanics, seismic, volcanic, etc.)
 * implements this interface to participate in the simulation pipeline.
 *
 * Modules are discovered and instantiated through the PhysicsModuleRegistry.
 * The Simulator calls lifecycle methods in order:
 *   1. configure()     — read parameters from config
 *   2. validateConfig() — check parameter consistency
 *   3. createKernels()  — return PhysicsKernel instances for FE assembly
 *   4. initialize()     — allocate PETSc objects, set ICs
 *   5. preStep() / postStep() — per-timestep hooks
 *   6. finalize()       — cleanup
 */

#ifndef PHYSICS_MODULE_INTERFACE_HPP
#define PHYSICS_MODULE_INTERFACE_HPP

#include "core/FSRM.hpp"
#include "physics/PhysicsKernel.hpp"
#include <memory>
#include <string>
#include <vector>

namespace FSRM {

// Forward declarations
class Simulator;
class ConfigReader;

/**
 * @brief Interface for all physics modules in FSRM
 */
class PhysicsModuleInterface {
public:
    virtual ~PhysicsModuleInterface() = default;

    /// Unique string identifier for this module (e.g., "reservoir.fluid_flow")
    virtual std::string name() const = 0;

    /// Human-readable description
    virtual std::string description() const = 0;

    /// Physics types this module provides
    virtual std::vector<PhysicsType> providedPhysics() const = 0;

    /// Physics types this module requires from other modules
    virtual std::vector<PhysicsType> requiredPhysics() const = 0;

    /// Configure from config reader. Return 0 on success.
    virtual PetscErrorCode configure(const ConfigReader& config) = 0;

    /// Validate configuration consistency. Return 0 on success.
    virtual PetscErrorCode validateConfig() const = 0;

    /// Return PhysicsKernel instances for FE assembly
    virtual std::vector<std::shared_ptr<PhysicsKernel>> createKernels() = 0;

    /// Initialize PETSc objects, set initial conditions on the DM
    virtual PetscErrorCode initialize(DM dm, Vec solution) = 0;

    /// Called before each timestep
    virtual PetscErrorCode preStep(PetscReal time, PetscReal dt) {
        (void)time; (void)dt;
        return 0;
    }

    /// Called after each timestep
    virtual PetscErrorCode postStep(PetscReal time, Vec solution) {
        (void)time; (void)solution;
        return 0;
    }

    /// Cleanup
    virtual PetscErrorCode finalize() { return 0; }

    /// Whether this module is enabled for the current configuration
    virtual bool isEnabled(const SimulationConfig& config) const = 0;
};

} // namespace FSRM

#endif // PHYSICS_MODULE_INTERFACE_HPP
