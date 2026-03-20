/**
 * @file test_physics_kernel_base.cpp
 * @brief Contract tests for PhysicsKernel base behavior
 */

#include <gtest/gtest.h>
#include "physics/PhysicsKernel.hpp"
#include <cmath>
#include <cstring>
#include <mpi.h>

using namespace FSRM;

namespace {

class ConcreteKernel : public PhysicsKernel {
public:
    explicit ConcreteKernel(PhysicsType t) : PhysicsKernel(t) {}

    PetscErrorCode setup(DM, PetscFE) override { return 0; }

    void residual(const PetscScalar u[], const PetscScalar[], const PetscScalar[],
                  const PetscScalar[], const PetscReal[], PetscScalar f[]) override {
        f[0] = u[0];
    }

    void jacobian(const PetscScalar[], const PetscScalar[], const PetscScalar[],
                  const PetscScalar[], const PetscReal[], PetscScalar J[]) override {
        J[0] = 1.0;
    }

    int getNumFields() const override { return 1; }
    int getNumComponents(int field) const override {
        (void)field;
        return 1;
    }
};

} // namespace

class PhysicsKernelInterfaceTest : public ::testing::Test {
protected:
    void SetUp() override { MPI_Comm_rank(PETSC_COMM_WORLD, &rank); }
    int rank = 0;
};

TEST_F(PhysicsKernelInterfaceTest, GetType) {
    ConcreteKernel k(PhysicsType::THERMAL);
    EXPECT_EQ(k.getType(), PhysicsType::THERMAL);
}

TEST_F(PhysicsKernelInterfaceTest, GetTypeNameNonNullForCommonTypes) {
    const PhysicsType types[] = {PhysicsType::FLUID_FLOW,   PhysicsType::GEOMECHANICS,
                                 PhysicsType::THERMAL,       PhysicsType::ELASTODYNAMICS,
                                 PhysicsType::POROELASTODYNAMICS, PhysicsType::PARTICLE_TRANSPORT};
    for (PhysicsType t : types) {
        ConcreteKernel k(t);
        const char* nm = k.getTypeName();
        ASSERT_NE(nm, nullptr);
        EXPECT_GT(std::strlen(nm), 0u);
    }
}

TEST_F(PhysicsKernelInterfaceTest, GetTotalDOF) {
    ConcreteKernel k(PhysicsType::FLUID_FLOW);
    EXPECT_EQ(k.getTotalDOF(), 1);
}

TEST_F(PhysicsKernelInterfaceTest, DefaultExecutionTraitsCpuOnly) {
    ConcreteKernel k(PhysicsType::FLUID_FLOW);
    const ExecutionTraits& tr = k.getExecutionTraits();
    EXPECT_FALSE(tr.supports_gpu);
    EXPECT_FALSE(tr.requires_gpu);
}

TEST_F(PhysicsKernelInterfaceTest, DefaultCapabilitiesAllCpu) {
    ConcreteKernel k(PhysicsType::FLUID_FLOW);
    EXPECT_TRUE(k.hasCapability(KernelCapability::RESIDUAL));
    EXPECT_TRUE(k.hasCapability(KernelCapability::JACOBIAN));
    EXPECT_FALSE(k.supportsGPU());
}

TEST_F(PhysicsKernelInterfaceTest, SelectBackendCpuForCpuOnlyKernel) {
    ConcreteKernel k(PhysicsType::FLUID_FLOW);
    EXPECT_EQ(k.selectBackend(1000000), ExecutionBackend::CPU);
}

TEST_F(PhysicsKernelInterfaceTest, ResidualAndJacobianIdentity) {
    ConcreteKernel k(PhysicsType::FLUID_FLOW);
    PetscScalar u[1] = {2.5};
    PetscScalar u_t[1] = {};
    PetscScalar u_x[3] = {};
    PetscScalar a[1] = {1.0};
    PetscReal x[3] = {0.5, 0.5, 0.5};
    PetscScalar f[1] = {};
    PetscScalar J[1] = {};
    k.residual(u, u_t, u_x, a, x, f);
    k.jacobian(u, u_t, u_x, a, x, J);
    EXPECT_DOUBLE_EQ(PetscRealPart(f[0]), 2.5);
    EXPECT_DOUBLE_EQ(PetscRealPart(J[0]), 1.0);
}
