/**
 * @file OceanPhysicsModule.hpp
 * @brief Physics module wrapper for ocean circulation and coastal hydrodynamics
 */

#ifndef OCEAN_PHYSICS_MODULE_HPP
#define OCEAN_PHYSICS_MODULE_HPP

#include "core/PhysicsModuleInterface.hpp"
#include "core/PhysicsModuleRegistry.hpp"

namespace FSRM {

class OceanPhysicsModule : public PhysicsModuleInterface {
public:
    std::string name() const override { return "ocean_physics"; }
    std::string description() const override {
        return "Ocean circulation, coastal hydrodynamics, and marine physics";
    }

    std::vector<PhysicsType> providedPhysics() const override {
        return {
            PhysicsType::OCEAN_CIRCULATION,
            PhysicsType::COASTAL_HYDRODYNAMICS,
            PhysicsType::STORM_SURGE
        };
    }
    std::vector<PhysicsType> requiredPhysics() const override { return {}; }

    PetscErrorCode configure(const ConfigReader& config) override {
        (void)config;
        return 0;
    }
    PetscErrorCode validateConfig() const override { return 0; }

    std::vector<std::shared_ptr<PhysicsKernel>> createKernels() override {
        return {};
    }

    PetscErrorCode initialize(DM dm, Vec solution) override {
        (void)dm;
        (void)solution;
        return 0;
    }

    bool isEnabled(const SimulationConfig& config) const override {
        return config.enable_ocean_circulation || config.enable_coastal_hydrodynamics
            || config.enable_storm_surge;
    }
};

} // namespace FSRM

#endif // OCEAN_PHYSICS_MODULE_HPP
