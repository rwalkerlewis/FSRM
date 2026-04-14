/**
 * @file PoroelastodynamicsModule.hpp
 * @brief Physics module for coupled fluid-solid wave propagation (Biot)
 */

#ifndef POROELASTODYNAMICS_MODULE_HPP
#define POROELASTODYNAMICS_MODULE_HPP

#include "core/PhysicsModuleInterface.hpp"
#include "core/PhysicsModuleRegistry.hpp"

namespace FSRM {

class PoroelastodynamicsModule : public PhysicsModuleInterface {
public:
    std::string name() const override { return "poroelastodynamics"; }
    std::string description() const override {
        return "Coupled fluid-solid wave propagation (Biot equations)";
    }

    std::vector<PhysicsType> providedPhysics() const override {
        return {PhysicsType::POROELASTODYNAMICS};
    }
    std::vector<PhysicsType> requiredPhysics() const override { return {}; }

    PetscErrorCode configure(const ConfigReader& config) override {
        (void)config;
        return 0;
    }
    PetscErrorCode validateConfig() const override { return 0; }

    std::vector<std::shared_ptr<PhysicsKernel>> createKernels() override {
        return {std::make_shared<PoroelastodynamicsKernel>()};
    }

    PetscErrorCode initialize(DM dm, Vec solution) override {
        (void)dm; (void)solution;
        return 0;
    }

    bool isEnabled(const SimulationConfig& config) const override {
        return config.enable_poroelastodynamics ||
               (config.enable_geomechanics &&
                config.solid_model == SolidModelType::POROELASTODYNAMIC);
    }
};

} // namespace FSRM

#endif // POROELASTODYNAMICS_MODULE_HPP
