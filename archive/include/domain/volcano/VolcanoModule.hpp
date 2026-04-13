/**
 * @file VolcanoModule.hpp
 * @brief Physics module wrapper for volcanic physics models
 */

#ifndef VOLCANO_MODULE_HPP
#define VOLCANO_MODULE_HPP

#include "core/PhysicsModuleInterface.hpp"
#include "core/PhysicsModuleRegistry.hpp"

namespace FSRM {

class VolcanoModule : public PhysicsModuleInterface {
public:
    std::string name() const override { return "volcano"; }
    std::string description() const override {
        return "Volcanic physics (magma chamber, conduit flow, eruption column)";
    }

    std::vector<PhysicsType> providedPhysics() const override {
        return {
            PhysicsType::MAGMA_CHAMBER,
            PhysicsType::CONDUIT_FLOW,
            PhysicsType::ERUPTION_COLUMN,
            PhysicsType::PYROCLASTIC_FLOW,
            PhysicsType::LAVA_FLOW
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
        return config.enable_magma_chamber || config.enable_conduit_flow
            || config.enable_eruption_column || config.enable_pyroclastic_flow
            || config.enable_lava_flow;
    }
};

} // namespace FSRM

#endif // VOLCANO_MODULE_HPP
