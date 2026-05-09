// tests/unit/test_mpi_correctness.cpp
//
// Pass-12 sanity check that the linked MPI library produces correct
// reduction results. This is a guard against MPI-stack misconfiguration
// (e.g. a Docker image that links a partial OpenMPI build) which
// otherwise surfaces as confusing simulation failures further into the
// stack.

#include <gtest/gtest.h>
#include <mpi.h>
#include <petsc.h>

#include <vector>

class MPICorrectnessTest : public ::testing::Test
{
protected:
  int rank_ = 0;
  int size_ = 1;

  void SetUp() override
  {
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank_);
    MPI_Comm_size(PETSC_COMM_WORLD, &size_);
  }
};

TEST_F(MPICorrectnessTest, AllreduceSumIsCorrect)
{
  // Each rank contributes its rank id; sum across N ranks is N(N-1)/2.
  const int local = rank_;
  int global = -1;
  ASSERT_EQ(MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_SUM,
                          PETSC_COMM_WORLD),
            MPI_SUCCESS)
      << "MPI_Allreduce returned error";
  const int expected = (size_ * (size_ - 1)) / 2;
  EXPECT_EQ(global, expected)
      << "MPI_Allreduce sum incorrect: expected " << expected
      << " got " << global << " on " << size_ << " ranks";
}

TEST_F(MPICorrectnessTest, AllreduceMaxIsCorrect)
{
  const double local = static_cast<double>(rank_);
  double global = -1.0;
  ASSERT_EQ(MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX,
                          PETSC_COMM_WORLD),
            MPI_SUCCESS);
  EXPECT_DOUBLE_EQ(global, static_cast<double>(size_ - 1))
      << "MPI_Allreduce max incorrect on " << size_ << " ranks";
}

TEST_F(MPICorrectnessTest, AllreduceVectorSumIsCorrect)
{
  // 4-element vector: each rank contributes (rank, rank+1, rank+2, rank+3).
  // Element-wise sum is sum_i (rank+i) = 4*size + offset_sum;
  // we just check element 0 against the integer formula.
  const std::vector<int> local = {rank_, rank_ + 1, rank_ + 2, rank_ + 3};
  std::vector<int> global(4, -1);
  ASSERT_EQ(MPI_Allreduce(local.data(), global.data(), 4, MPI_INT,
                          MPI_SUM, PETSC_COMM_WORLD),
            MPI_SUCCESS);
  const int expected_0 = (size_ * (size_ - 1)) / 2;       // sum of ranks
  EXPECT_EQ(global[0], expected_0);
  EXPECT_EQ(global[3], expected_0 + 3 * size_);
}
