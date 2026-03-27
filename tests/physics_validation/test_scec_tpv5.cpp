/**
 * @file test_scec_tpv5.cpp
 * @brief Placeholder for SCEC TPV5 (vertical strike-slip, rate-state) benchmark validation
 *
 * When implemented, this test should run a full dynamic rupture simulation and compare
 * time histories of slip rate at SCEC-specified recording points against the community
 * reference solution, requiring agreement within 5% (peak slip rate and/or RMS) over the
 * rupture window. Currently skipped because it requires the full FSRM/PETSc runtime,
 * verified mesh, absorbing boundaries, and reference datasets.
 */

#include <gtest/gtest.h>
#include <petscsys.h>

class SCECTPV5Test : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    }

    int rank = 0;
};

TEST_F(SCECTPV5Test, BenchmarkSlipRateWithinFivePercent) {
    GTEST_SKIP() << "Full SCEC TPV5 simulation not wired in CI; future test will compare slip "
                    "rate at recording stations to the SCEC reference within 5%";
}

TEST_F(SCECTPV5Test, DocumentedVerificationTargets) {
    GTEST_SKIP()
        << "TPV5 validation targets: (1) spontaneous rupture on a vertical strike-slip fault "
           "with rate-and-state friction, (2) slip-rate time series at prescribed fault "
           "receivers vs. published benchmark, (3) pass criterion max|ΔV|/V_ref < 5% on peaks "
           "or equivalent integrated metric — requires simulation driver + reference data";
}
