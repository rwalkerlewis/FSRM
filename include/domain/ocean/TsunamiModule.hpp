/**
 * @file TsunamiModule.hpp
 * @brief Physics module wrapper for shallow water tsunami propagation
 */

#ifndef TSUNAMI_MODULE_HPP
#define TSUNAMI_MODULE_HPP

#include "core/PhysicsModuleInterface.hpp"
#include "core/PhysicsModuleRegistry.hpp"

namespace FSRM {

class TsunamiModule : public PhysicsModuleInterface {
public:
    std::string name() const override { return "tsunami"; }
    std::string description() const override {
        return "Shallow water tsunami wave propagation";
    }

    std::vector<PhysicsType> providedPhysics() const override {
        return {PhysicsType::TSUNAMI};
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
        return config.enable_tsunami;
    }
};

} // namespace FSRM

#endif // TSUNAMI_MODULE_HPP
