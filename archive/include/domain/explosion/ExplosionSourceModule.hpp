/**
 * @file ExplosionSourceModule.hpp
 * @brief Physics module for explosion source models and near-field damage
 */

#ifndef EXPLOSION_SOURCE_MODULE_HPP
#define EXPLOSION_SOURCE_MODULE_HPP

#include "core/PhysicsModuleInterface.hpp"
#include "core/PhysicsModuleRegistry.hpp"
#include "physics/ExplosionDamageKernels.hpp"

namespace FSRM {

class ExplosionSourceModule : public PhysicsModuleInterface {
public:
    std::string name() const override { return "explosion_source"; }
    std::string description() const override {
        return "Nuclear/chemical explosion source, near-field damage, crater formation";
    }

    std::vector<PhysicsType> providedPhysics() const override {
        return {
            PhysicsType::EXPLOSION_SOURCE,
            PhysicsType::NEAR_FIELD_DAMAGE,
            PhysicsType::HYDRODYNAMIC,
            PhysicsType::CRATER_FORMATION
        };
    }
    std::vector<PhysicsType> requiredPhysics() const override { return {}; }

    PetscErrorCode configure(const ConfigReader& config) override {
        (void)config;
        return 0;
    }
    PetscErrorCode validateConfig() const override { return 0; }

    std::vector<std::shared_ptr<PhysicsKernel>> createKernels() override {
        std::vector<std::shared_ptr<PhysicsKernel>> kernels;
        kernels.push_back(std::make_shared<ExplosionSourceKernel>());
        if (enable_near_field_) {
            kernels.push_back(std::make_shared<NearFieldDamageKernel>());
        }
        if (enable_crater_) {
            kernels.push_back(std::make_shared<CraterFormationKernel>());
        }
        return kernels;
    }

    PetscErrorCode initialize(DM dm, Vec solution) override {
        (void)dm; (void)solution;
        return 0;
    }

    bool isEnabled(const SimulationConfig& config) const override {
        return config.enable_explosion_source;
    }

private:
    bool enable_near_field_ = true;
    bool enable_crater_ = false;
};

} // namespace FSRM

#endif // EXPLOSION_SOURCE_MODULE_HPP
