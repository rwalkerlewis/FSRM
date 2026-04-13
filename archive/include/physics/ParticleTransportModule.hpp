/**
 * @file ParticleTransportModule.hpp
 * @brief Physics module for particle/proppant/tracer transport
 */

#ifndef PARTICLE_TRANSPORT_MODULE_HPP
#define PARTICLE_TRANSPORT_MODULE_HPP

#include "core/PhysicsModuleInterface.hpp"
#include "core/PhysicsModuleRegistry.hpp"

namespace FSRM {

class ParticleTransportModule : public PhysicsModuleInterface {
public:
    std::string name() const override { return "particle_transport"; }
    std::string description() const override {
        return "Proppant, tracer, and solids transport";
    }

    std::vector<PhysicsType> providedPhysics() const override {
        return {PhysicsType::PARTICLE_TRANSPORT};
    }
    std::vector<PhysicsType> requiredPhysics() const override { return {}; }

    PetscErrorCode configure(const ConfigReader& config) override {
        (void)config;
        return 0;
    }
    PetscErrorCode validateConfig() const override { return 0; }

    std::vector<std::shared_ptr<PhysicsKernel>> createKernels() override {
        return {std::make_shared<ParticleTransportKernel>()};
    }

    PetscErrorCode initialize(DM dm, Vec solution) override {
        (void)dm; (void)solution;
        return 0;
    }

    bool isEnabled(const SimulationConfig& config) const override {
        return config.enable_particle_transport;
    }
};

} // namespace FSRM

#endif // PARTICLE_TRANSPORT_MODULE_HPP
