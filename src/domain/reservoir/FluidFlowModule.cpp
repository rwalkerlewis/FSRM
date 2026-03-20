/**
 * @file FluidFlowModule.cpp
 * @brief Implementation of FluidFlowModule
 */

#include "domain/reservoir/FluidFlowModule.hpp"

namespace FSRM {

std::vector<std::shared_ptr<PhysicsKernel>> FluidFlowModule::createKernels() {
    std::vector<std::shared_ptr<PhysicsKernel>> kernels;
    switch (model_type_) {
        case FluidModelType::SINGLE_COMPONENT:
            kernels.push_back(std::make_shared<SinglePhaseFlowKernel>());
            break;
        case FluidModelType::BLACK_OIL:
            kernels.push_back(std::make_shared<BlackOilKernel>());
            break;
        case FluidModelType::COMPOSITIONAL:
            kernels.push_back(std::make_shared<CompositionalKernel>(3));
            break;
        default:
            break;
    }
    return kernels;
}

FSRM_REGISTER_PHYSICS_MODULE(FluidFlowModule)

} // namespace FSRM
