/**
 * @file FluidFlowModule.hpp
 * @brief Physics module for fluid flow (single-phase, black oil, compositional)
 */

#ifndef FLUID_FLOW_MODULE_HPP
#define FLUID_FLOW_MODULE_HPP

#include "core/PhysicsModuleInterface.hpp"
#include "core/PhysicsModuleRegistry.hpp"

namespace FSRM {

class FluidFlowModule : public PhysicsModuleInterface {
public:
    std::string name() const override { return "fluid_flow"; }
    std::string description() const override {
        return "Darcy fluid flow (single-phase, black oil, compositional)";
    }

    std::vector<PhysicsType> providedPhysics() const override {
        return {PhysicsType::FLUID_FLOW};
    }
    std::vector<PhysicsType> requiredPhysics() const override { return {}; }

    PetscErrorCode configure(const ConfigReader& config) override {
        (void)config;
        return 0;
    }
    PetscErrorCode validateConfig() const override { return 0; }

    std::vector<std::shared_ptr<PhysicsKernel>> createKernels() override;

    PetscErrorCode initialize(DM dm, Vec solution) override {
        (void)dm; (void)solution;
        return 0;
    }

    bool isEnabled(const SimulationConfig& config) const override {
        return config.fluid_model != FluidModelType::NONE;
    }

private:
    FluidModelType model_type_ = FluidModelType::SINGLE_COMPONENT;
};

} // namespace FSRM

#endif // FLUID_FLOW_MODULE_HPP
