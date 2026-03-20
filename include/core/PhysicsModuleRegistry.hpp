/**
 * @file PhysicsModuleRegistry.hpp
 * @brief Registry for physics modules with self-registration support
 *
 * Modules register themselves at static initialization time using the
 * FSRM_REGISTER_PHYSICS_MODULE macro. The Simulator queries the registry
 * to discover available physics. This replaces the boolean enable_ flags
 * and the switch cascade in Simulator::setupPhysics().
 */

#ifndef PHYSICS_MODULE_REGISTRY_HPP
#define PHYSICS_MODULE_REGISTRY_HPP

#include "core/PhysicsModuleInterface.hpp"
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace FSRM {

class PhysicsModuleRegistry {
public:
    using FactoryFn = std::function<std::unique_ptr<PhysicsModuleInterface>()>;

    /// Get the singleton instance (Meyer's singleton — safe for static init)
    static PhysicsModuleRegistry& instance();

    /// Register a module factory. Called at static init time.
    void registerModule(const std::string& name, FactoryFn factory);

    /// Get all registered module names
    std::vector<std::string> availableModules() const;

    /// Check if a module is registered
    bool hasModule(const std::string& name) const;

    /// Create a module instance by name. Returns nullptr if not found.
    std::unique_ptr<PhysicsModuleInterface> createModule(const std::string& name) const;

    /// Create all modules that are enabled for the given config
    std::vector<std::unique_ptr<PhysicsModuleInterface>>
    createEnabledModules(const SimulationConfig& config) const;

    /// Create modules from explicit list of module names
    std::vector<std::unique_ptr<PhysicsModuleInterface>>
    createModulesByName(const std::vector<std::string>& names) const;

    /// Check dependency graph: are all required physics provided?
    /// Returns 0 if all dependencies are satisfied, non-zero otherwise.
    PetscErrorCode validateDependencies(
        const std::vector<PhysicsModuleInterface*>& modules) const;

    /// Clear all registrations (for testing only)
    void clear();

private:
    PhysicsModuleRegistry() = default;
    PhysicsModuleRegistry(const PhysicsModuleRegistry&) = delete;
    PhysicsModuleRegistry& operator=(const PhysicsModuleRegistry&) = delete;

    std::map<std::string, FactoryFn> factories_;
};

/// Macro for module self-registration at static initialization time
#define FSRM_REGISTER_PHYSICS_MODULE(ModuleClass)                             \
    namespace {                                                               \
    static const bool ModuleClass##_registered_ = []() {                      \
        ::FSRM::PhysicsModuleRegistry::instance().registerModule(             \
            ModuleClass().name(),                                             \
            []() -> std::unique_ptr<::FSRM::PhysicsModuleInterface> {         \
                return std::make_unique<ModuleClass>();                        \
            });                                                               \
        return true;                                                          \
    }();                                                                      \
    } // anonymous namespace

} // namespace FSRM

#endif // PHYSICS_MODULE_REGISTRY_HPP
