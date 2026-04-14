/**
 * @file ElastodynamicsModule.hpp
 * @brief Physics module for elastic wave propagation with inertia
 */

#ifndef ELASTODYNAMICS_MODULE_HPP
#define ELASTODYNAMICS_MODULE_HPP

#include "core/PhysicsModuleInterface.hpp"
#include "core/PhysicsModuleRegistry.hpp"

namespace FSRM {

class ElastodynamicsModule : public PhysicsModuleInterface {
public:
    std::string name() const override { return "elastodynamics"; }
    std::string description() const override {
        return "Elastic wave propagation with inertia";
    }

    std::vector<PhysicsType> providedPhysics() const override {
        return {PhysicsType::ELASTODYNAMICS};
    }
    std::vector<PhysicsType> requiredPhysics() const override { return {}; }

    PetscErrorCode configure(const ConfigReader& config) override {
        (void)config;
        return 0;
    }
    PetscErrorCode validateConfig() const override { return 0; }

    std::vector<std::shared_ptr<PhysicsKernel>> createKernels() override {
        return {std::make_shared<ElastodynamicsKernel>()};
    }

    PetscErrorCode initialize(DM dm, Vec solution) override {
        (void)dm; (void)solution;
        return 0;
    }

    bool isEnabled(const SimulationConfig& config) const override {
        return config.enable_elastodynamics ||
               (config.enable_geomechanics &&
                config.solid_model == SolidModelType::ELASTODYNAMIC);
    }
};

} // namespace FSRM

#endif // ELASTODYNAMICS_MODULE_HPP
