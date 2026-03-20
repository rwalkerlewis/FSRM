/**
 * @file test_memory.cpp
 * @brief Repeated kernel evaluations should not grow memory without bound
 */

#include "core/FSRM.hpp"
#include "physics/PhysicsKernel.hpp"
#include <gtest/gtest.h>
#include <petscsys.h>

using namespace FSRM;

class MemoryStabilityTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    }

    int rank = 0;
    void* scratch_ = nullptr;
};

TEST_F(MemoryStabilityTest, RepeatedResidualEvaluationsStayBounded) {
    PetscLogDouble mem_before = 0.0;
    PetscLogDouble mem_after = 0.0;
    PetscMemoryGetCurrentUsage(&mem_before);

    SinglePhaseFlowKernel kernel;
    kernel.setProperties(0.2, 1.0e-13, 1.0e-9, 1.0e-3, 1000.0);

    PetscScalar u[1] = {0.0};
    PetscScalar u_t[1] = {0.0};
    PetscScalar u_x[3] = {0.0};
    PetscScalar a[1] = {0.0};
    PetscReal x[3] = {0.0, 0.0, 0.0};
    PetscScalar f[1] = {0.0};

    for (int round = 0; round < 20; ++round) {
        for (int i = 0; i < 500; ++i) {
            u[0] = PetscScalar(static_cast<double>(i + round));
            kernel.residual(u, u_t, u_x, a, x, f);
        }
        PetscMalloc(4096, &scratch_);
        PetscFree(scratch_);
        scratch_ = nullptr;
    }

    PetscMemoryGetCurrentUsage(&mem_after);
    const PetscLogDouble delta = mem_after - mem_before;
    if (rank == 0) {
        RecordProperty("petsc_mem_delta_bytes", static_cast<double>(delta));
    }
    EXPECT_LT(delta, 50.0 * 1024.0 * 1024.0);
}
