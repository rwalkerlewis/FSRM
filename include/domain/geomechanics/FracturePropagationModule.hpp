/**
 * @file FracturePropagationModule.hpp
 * @brief Physics module for fracture propagation (cohesive zone / LEFM)
 */

#ifndef FRACTURE_PROPAGATION_MODULE_HPP
#define FRACTURE_PROPAGATION_MODULE_HPP

#include "core/PhysicsModuleInterface.hpp"
#include "core/PhysicsModuleRegistry.hpp"

namespace FSRM {

class FracturePropagationModule : public PhysicsModuleInterface {
public:
    std::string name() const override { return "fracture_propagation"; }
    std::string description() const override {
        return "Cohesive zone / LEFM fracture propagation";
    }

    std::vector<PhysicsType> providedPhysics() const override {
        return {PhysicsType::FRACTURE_PROPAGATION};
    }
    std::vector<PhysicsType> requiredPhysics() const override { return {}; }

    PetscErrorCode configure(const ConfigReader& config) override {
        (void)config;
        return 0;
    }
    PetscErrorCode validateConfig() const override { return 0; }

    std::vector<std::shared_ptr<PhysicsKernel>> createKernels() override {
        return {std::make_shared<FracturePropagationKernel>()};
    }

    PetscErrorCode initialize(DM dm, Vec solution) override {
        (void)dm; (void)solution;
        return 0;
    }

    bool isEnabled(const SimulationConfig& config) const override {
        return config.enable_fractures;
    }
};

} // namespace FSRM

#endif // FRACTURE_PROPAGATION_MODULE_HPP
