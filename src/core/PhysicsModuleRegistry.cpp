/**
 * @file PhysicsModuleRegistry.cpp
 * @brief Implementation of the PhysicsModuleRegistry singleton
 */

#include "core/PhysicsModuleRegistry.hpp"
#include <algorithm>
#include <set>
#include <sstream>
#include <iostream>

namespace FSRM {

PhysicsModuleRegistry& PhysicsModuleRegistry::instance() {
    static PhysicsModuleRegistry registry;
    return registry;
}

void PhysicsModuleRegistry::registerModule(const std::string& name, FactoryFn factory) {
    if (factories_.count(name)) {
        std::cerr << "Warning: PhysicsModule '" << name
                  << "' is already registered. Skipping duplicate." << std::endl;
        return;
    }
    factories_[name] = std::move(factory);
}

std::vector<std::string> PhysicsModuleRegistry::availableModules() const {
    std::vector<std::string> names;
    names.reserve(factories_.size());
    for (const auto& kv : factories_) {
        names.push_back(kv.first);
    }
    return names;
}

bool PhysicsModuleRegistry::hasModule(const std::string& name) const {
    return factories_.count(name) > 0;
}

std::unique_ptr<PhysicsModuleInterface>
PhysicsModuleRegistry::createModule(const std::string& name) const {
    auto it = factories_.find(name);
    if (it == factories_.end()) {
        return nullptr;
    }
    return it->second();
}

std::vector<std::unique_ptr<PhysicsModuleInterface>>
PhysicsModuleRegistry::createEnabledModules(const SimulationConfig& config) const {
    std::vector<std::unique_ptr<PhysicsModuleInterface>> result;
    for (const auto& kv : factories_) {
        auto module = kv.second();
        if (module && module->isEnabled(config)) {
            result.push_back(std::move(module));
        }
    }
    return result;
}

std::vector<std::unique_ptr<PhysicsModuleInterface>>
PhysicsModuleRegistry::createModulesByName(const std::vector<std::string>& names) const {
    std::vector<std::unique_ptr<PhysicsModuleInterface>> result;
    for (const auto& name : names) {
        auto module = createModule(name);
        if (module) {
            result.push_back(std::move(module));
        } else {
            std::cerr << "Error: Unknown physics module '" << name << "'" << std::endl;
        }
    }
    return result;
}

PetscErrorCode PhysicsModuleRegistry::validateDependencies(
    const std::vector<PhysicsModuleInterface*>& modules) const {
    PetscFunctionBeginUser;

    // Collect all provided physics types
    std::set<PhysicsType> provided;
    for (const auto* mod : modules) {
        for (auto pt : mod->providedPhysics()) {
            provided.insert(pt);
        }
    }

    // Check that every required physics type is provided
    std::vector<std::string> missing;
    for (const auto* mod : modules) {
        for (auto pt : mod->requiredPhysics()) {
            if (provided.find(pt) == provided.end()) {
                std::ostringstream oss;
                oss << "Module '" << mod->name()
                    << "' requires physics type " << static_cast<int>(pt)
                    << " which is not provided by any enabled module.";
                missing.push_back(oss.str());
            }
        }
    }

    if (!missing.empty()) {
        std::cerr << "Dependency validation failed:" << std::endl;
        for (const auto& msg : missing) {
            std::cerr << "  " << msg << std::endl;
        }
        PetscFunctionReturn(PETSC_ERR_ARG_INCOMP);
    }

    PetscFunctionReturn(0);
}

void PhysicsModuleRegistry::clear() {
    factories_.clear();
}

} // namespace FSRM
