/**
 * @file test_neural_stubs.cpp
 * @brief Smoke tests for experimental ML / AMR headers
 */

#include "core/FSRM.hpp"
#include "experimental/MultiFidelityLearning.hpp"
#include "experimental/NeuralAMR.hpp"
#include <gtest/gtest.h>

using namespace FSRM::ML;

class ExperimentalNeuralTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    }
    int rank = 0;
};

TEST_F(ExperimentalNeuralTest, MultiFidelityControllerConstructs) {
    MultiFidelityController::MFConfig cfg{};
    MultiFidelityController ctrl(cfg);
    Tensor in({1}, {1.0});
    Tensor out = ctrl.predict(in);
    EXPECT_EQ(out.data.size(), 1u);
    (void)rank;
}

TEST_F(ExperimentalNeuralTest, NeuralAMRControllerConstructs) {
    NeuralAMRController::AMRConfig cfg{};
    NeuralAMRController amr(cfg);
    EXPECT_FALSE(amr.shouldAdapt(0, 0.0));
}

TEST_F(ExperimentalNeuralTest, FidelityLevelDefaultConstructible) {
    FidelityLevel L{};
    EXPECT_EQ(L.level, 0);
}
