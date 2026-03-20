/**
 * @file GeomechanicsModule.cpp
 * @brief Implementation of GeomechanicsModule
 */

#include "domain/geomechanics/GeomechanicsModule.hpp"

namespace FSRM {

std::vector<std::shared_ptr<PhysicsKernel>> GeomechanicsModule::createKernels() {
    std::vector<std::shared_ptr<PhysicsKernel>> kernels;
    kernels.push_back(std::make_shared<GeomechanicsKernel>(model_type_));
    return kernels;
}

FSRM_REGISTER_PHYSICS_MODULE(GeomechanicsModule)

} // namespace FSRM
