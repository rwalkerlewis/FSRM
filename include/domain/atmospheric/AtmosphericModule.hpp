/**
 * @file AtmosphericModule.hpp
 * @brief Physics module wrapper for atmospheric blast, acoustic, and infrasound
 */

#ifndef ATMOSPHERIC_MODULE_HPP
#define ATMOSPHERIC_MODULE_HPP

#include "core/PhysicsModuleInterface.hpp"
#include "core/PhysicsModuleRegistry.hpp"

namespace FSRM {

class AtmosphericModule : public PhysicsModuleInterface {
public:
    std::string name() const override { return "atmospheric"; }
    std::string description() const override {
        return "Atmospheric blast, acoustic, and infrasound propagation";
    }

    std::vector<PhysicsType> providedPhysics() const override {
        return {
            PhysicsType::ATMOSPHERIC_BLAST,
            PhysicsType::ATMOSPHERIC_ACOUSTIC,
            PhysicsType::INFRASOUND
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
        return config.enable_atmospheric_blast || config.enable_atmospheric_acoustic
            || config.enable_infrasound;
    }
};

} // namespace FSRM

#endif // ATMOSPHERIC_MODULE_HPP
