/**
 * @file ThermalModule.hpp
 * @brief Physics module for thermal conduction and convection
 */

#ifndef THERMAL_MODULE_HPP
#define THERMAL_MODULE_HPP

#include "core/PhysicsModuleInterface.hpp"
#include "core/PhysicsModuleRegistry.hpp"

namespace FSRM {

class ThermalModule : public PhysicsModuleInterface {
public:
    std::string name() const override { return "thermal"; }
    std::string description() const override {
        return "Heat conduction and convection";
    }

    std::vector<PhysicsType> providedPhysics() const override {
        return {PhysicsType::THERMAL};
    }
    std::vector<PhysicsType> requiredPhysics() const override { return {}; }

    PetscErrorCode configure(const ConfigReader& config) override {
        (void)config;
        return 0;
    }
    PetscErrorCode validateConfig() const override { return 0; }

    std::vector<std::shared_ptr<PhysicsKernel>> createKernels() override {
        return {std::make_shared<ThermalKernel>()};
    }

    PetscErrorCode initialize(DM dm, Vec solution) override {
        (void)dm; (void)solution;
        return 0;
    }

    bool isEnabled(const SimulationConfig& config) const override {
        return config.enable_thermal;
    }
};

} // namespace FSRM

#endif // THERMAL_MODULE_HPP
