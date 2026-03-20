/**
 * @file FaultMechanicsModule.hpp
 * @brief Physics module wrapper for fault mechanics (service module, no PDE kernels)
 */

#ifndef FAULT_MECHANICS_MODULE_HPP
#define FAULT_MECHANICS_MODULE_HPP

#include "core/PhysicsModuleInterface.hpp"
#include "core/PhysicsModuleRegistry.hpp"

namespace FSRM {

class FaultMechanicsModule : public PhysicsModuleInterface {
public:
    std::string name() const override { return "fault_mechanics"; }
    std::string description() const override {
        return "Fault mechanics and slip modeling";
    }

    std::vector<PhysicsType> providedPhysics() const override { return {}; }
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
        return config.enable_faults;
    }
};

} // namespace FSRM

#endif // FAULT_MECHANICS_MODULE_HPP
