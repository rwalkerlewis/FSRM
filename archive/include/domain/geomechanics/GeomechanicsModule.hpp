/**
 * @file GeomechanicsModule.hpp
 * @brief Physics module for geomechanics (elastic, viscoelastic, poroelastic)
 */

#ifndef GEOMECHANICS_MODULE_HPP
#define GEOMECHANICS_MODULE_HPP

#include "core/PhysicsModuleInterface.hpp"
#include "core/PhysicsModuleRegistry.hpp"

namespace FSRM {

class GeomechanicsModule : public PhysicsModuleInterface {
public:
    std::string name() const override { return "geomechanics"; }
    std::string description() const override {
        return "Static geomechanics (elastic, viscoelastic, poroelastic)";
    }

    std::vector<PhysicsType> providedPhysics() const override {
        return {PhysicsType::GEOMECHANICS};
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
        return config.enable_geomechanics &&
               config.solid_model != SolidModelType::ELASTODYNAMIC &&
               config.solid_model != SolidModelType::POROELASTODYNAMIC;
    }

private:
    SolidModelType model_type_ = SolidModelType::ELASTIC;
};

} // namespace FSRM

#endif // GEOMECHANICS_MODULE_HPP
