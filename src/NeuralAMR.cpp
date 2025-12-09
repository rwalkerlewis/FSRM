#include "NeuralAMR.hpp"
#include <iostream>

namespace FSRM {
namespace ML {

NeuralAMRController::NeuralAMRController(const AMRConfig& config)
    : config_(config) {
    std::cout << "[NeuralAMR] Initialized\n";
}

bool NeuralAMRController::shouldAdapt(int /* step */, double /* time */) {
    return false;
}

Tensor NeuralAMRController::computeRefinementIndicator(DM /* dm */, Vec /* solution */, double /* time */) {
    return Tensor({1}, {0.5});
}

void NeuralAMRController::getRefinementFlags(DM /* dm */, Vec /* solution */, double /* time */,
                                             std::vector<PetscInt>& /* refine_cells */,
                                             std::vector<PetscInt>& /* coarsen_cells */) {
}

PetscErrorCode NeuralAMRController::adapt(DM& /* dm */, Vec& /* solution */, double /* time */) {
    return PETSC_SUCCESS;
}

void NeuralAMRController::train(const std::vector<std::vector<std::pair<DM, Vec>>>& /* trajectories */,
                                const std::vector<std::vector<double>>& /* times */,
                                const std::vector<std::vector<double>>& /* errors */) {
}

void NeuralAMRController::save(const std::string& /* path */) const {
}

void NeuralAMRController::load(const std::string& /* path */) {
}

} // namespace ML
} // namespace FSRM
