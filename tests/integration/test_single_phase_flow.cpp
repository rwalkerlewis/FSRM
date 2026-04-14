/**
 * @file test_single_phase_flow.cpp
 * @brief Integration test: single-phase pressure diffusion through TSSolve
 *
 * Tests the single-phase fluid flow equation:
 *   phi * ct * dP/dt - div( (k/mu) * grad(P) ) = 0
 *
 * Setup: 1D pressure diffusion in a 100m x 1m x 1m porous medium.
 * Dirichlet BCs: P(x=0) = 20 MPa, P(x=100) = 10 MPa.
 * Other boundaries: no-flow (natural Neumann BC).
 *
 * At steady state the analytical solution is linear:
 *   P(x) = 20e6 - (10e6 / 100) * x = 20e6 - 1e5 * x
 *
 * The test verifies:
 *   1. TSSolve completes without error
 *   2. Solution is nonzero
 *   3. Pressure monotonically decreases from x_min to x_max
 *   4. Center pressure approaches 15 MPa (within 3 MPa tolerance)
 */

#include <gtest/gtest.h>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <petscsys.h>
#include <petscdmplex.h>

#include "core/Simulator.hpp"
#include "core/FSRM.hpp"

using namespace FSRM;

class SinglePhaseFlowTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank_);
    }

    void writeConfig(const std::string& path)
    {
        if (rank_ != 0) return;
        std::ofstream cfg(path);
        cfg << "[SIMULATION]\n";
        cfg << "name = test_single_phase_flow\n";
        cfg << "start_time = 0.0\n";
        cfg << "end_time = 100.0\n";
        cfg << "dt_initial = 1.0\n";
        cfg << "dt_min = 0.1\n";
        cfg << "dt_max = 10.0\n";
        cfg << "max_timesteps = 100\n";
        cfg << "output_frequency = 1000\n";
        cfg << "fluid_model = SINGLE_COMPONENT\n";
        cfg << "solid_model = ELASTIC\n";
        cfg << "enable_geomechanics = false\n";
        cfg << "enable_elastodynamics = false\n";
        cfg << "enable_faults = false\n";
        cfg << "rtol = 1.0e-6\n";
        cfg << "atol = 1.0e-8\n";
        cfg << "max_nonlinear_iterations = 20\n";
        cfg << "\n[GRID]\n";
        cfg << "nx = 20\n";
        cfg << "ny = 1\n";
        cfg << "nz = 1\n";
        cfg << "Lx = 100.0\n";
        cfg << "Ly = 1.0\n";
        cfg << "Lz = 1.0\n";
        cfg << "\n[ROCK]\n";
        cfg << "density = 2500.0\n";
        cfg << "youngs_modulus = 10.0e9\n";
        cfg << "poissons_ratio = 0.25\n";
        cfg << "porosity = 0.20\n";
        cfg << "permeability_x = 100.0\n";  // mD (converted to m^2 by ConfigReader)
        cfg << "permeability_y = 100.0\n";
        cfg << "permeability_z = 100.0\n";
        cfg << "\n[FLUID]\n";
        cfg << "water_compressibility = 1.0e-9\n";
        cfg << "water_viscosity = 0.001\n";
        cfg << "\n[BOUNDARY_CONDITIONS]\n";
        cfg << "x_min = dirichlet_pressure\n";
        cfg << "x_min_pressure = 20.0e6\n";
        cfg << "x_max = dirichlet_pressure\n";
        cfg << "x_max_pressure = 10.0e6\n";
        cfg << "\n[ABSORBING_BC]\n";
        cfg << "enabled = false\n";
        cfg.close();
    }

    int rank_ = 0;
};

TEST_F(SinglePhaseFlowTest, PressureDiffusion1D)
{
    std::string config_path = "test_single_phase_flow.config";
    writeConfig(config_path);
    MPI_Barrier(PETSC_COMM_WORLD);

    Simulator sim(PETSC_COMM_WORLD);

    PetscOptionsClear(nullptr);
    PetscOptionsSetValue(nullptr, "-ts_type", "beuler");
    PetscOptionsSetValue(nullptr, "-pc_type", "lu");
    PetscOptionsSetValue(nullptr, "-ksp_type", "preonly");
    PetscOptionsSetValue(nullptr, "-snes_rtol", "1e-8");
    PetscOptionsSetValue(nullptr, "-snes_atol", "1e-10");
    PetscOptionsSetValue(nullptr, "-snes_max_it", "50");

    PetscErrorCode ierr;
    ierr = sim.initializeFromConfigFile(config_path);
    ASSERT_EQ(ierr, 0) << "initializeFromConfigFile failed";
    ierr = sim.setupDM();
    ASSERT_EQ(ierr, 0) << "setupDM failed";
    ierr = sim.labelBoundaries();
    ASSERT_EQ(ierr, 0) << "labelBoundaries failed";
    ierr = sim.setupFields();
    ASSERT_EQ(ierr, 0) << "setupFields failed";
    ierr = sim.setupPhysics();
    ASSERT_EQ(ierr, 0) << "setupPhysics failed";
    ierr = sim.setupTimeStepper();
    ASSERT_EQ(ierr, 0) << "setupTimeStepper failed";
    ierr = sim.setupSolvers();
    ASSERT_EQ(ierr, 0) << "setupSolvers failed";
    ierr = sim.setInitialConditions();
    ASSERT_EQ(ierr, 0) << "setInitialConditions failed";

    PetscPushErrorHandler(PetscReturnErrorHandler, nullptr);
    ierr = sim.run();
    PetscPopErrorHandler();
    ASSERT_EQ(ierr, 0) << "TSSolve failed for single-phase pressure diffusion";

    // Verify solution
    Vec sol = sim.getSolution();
    ASSERT_TRUE(sol != nullptr) << "Solution vector is null";

    PetscReal sol_norm = 0.0;
    VecNorm(sol, NORM_2, &sol_norm);
    EXPECT_GT(sol_norm, 0.0) << "Solution must be nonzero";
    EXPECT_TRUE(std::isfinite(sol_norm)) << "Solution norm must be finite";

    // Sample pressure at three x-locations: near x_min, center, near x_max
    DM dm = sim.getDM();
    PetscInt dim = 0;
    DMGetDimension(dm, &dim);

    Vec local_sol;
    DMGetLocalVector(dm, &local_sol);
    VecZeroEntries(local_sol);
    DMGlobalToLocal(dm, sol, INSERT_VALUES, local_sol);

    PetscSection section;
    DMGetLocalSection(dm, &section);

    // Iterate over cells to sample pressure at different x-locations
    PetscInt cStart, cEnd;
    DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd);

    double p_left = 0.0, p_center = 0.0, p_right = 0.0;
    double x_left = 1e30, x_center_dist = 1e30, x_right = -1e30;
    const double x_mid = 50.0;

    const PetscScalar *sol_array;
    VecGetArrayRead(local_sol, &sol_array);

    for (PetscInt c = cStart; c < cEnd; ++c) {
        // Get cell centroid
        PetscReal vol, centroid[3], normal[3];
        DMPlexComputeCellGeometryFVM(dm, c, &vol, centroid, normal);
        if (vol <= 0.0) continue;

        double cx = centroid[0];

        // Get pressure at this cell (field 0, component 0)
        PetscInt dof, off;
        PetscSectionGetDof(section, c, &dof);
        if (dof <= 0) continue;
        PetscSectionGetOffset(section, c, &off);
        double p = PetscRealPart(sol_array[off]);

        // Track leftmost, center, rightmost pressure
        if (cx < x_left) { x_left = cx; p_left = p; }
        if (cx > x_right) { x_right = cx; p_right = p; }
        if (std::abs(cx - x_mid) < x_center_dist) {
            x_center_dist = std::abs(cx - x_mid);
            p_center = p;
        }
    }

    VecRestoreArrayRead(local_sol, &sol_array);
    DMRestoreLocalVector(dm, &local_sol);

    if (rank_ == 0) {
        PetscPrintf(PETSC_COMM_WORLD,
            "Pressure samples: P(x=%.1f)=%.4e, P(x=%.1f)=%.4e, P(x=%.1f)=%.4e\n",
            x_left, p_left, x_mid, p_center, x_right, p_right);
    }

    // Assertions with quantitative tolerances
    // 1. Monotonicity: P(x_min) > P(center) > P(x_max)
    EXPECT_GT(p_left, p_center) << "Pressure at x_min must exceed pressure at center";
    EXPECT_GT(p_center, p_right) << "Pressure at center must exceed pressure at x_max";

    // 2. Center pressure should approach 15 MPa (within 3 MPa tolerance, given finite time)
    EXPECT_NEAR(p_center, 15.0e6, 3.0e6)
        << "Center pressure should approach 15 MPa at steady state";

    // 3. Left pressure should be near 20 MPa (Dirichlet BC)
    EXPECT_NEAR(p_left, 20.0e6, 2.0e6)
        << "Left-side pressure should be near the 20 MPa Dirichlet BC";

    // 4. Right pressure should be near 10 MPa (Dirichlet BC)
    EXPECT_NEAR(p_right, 10.0e6, 2.0e6)
        << "Right-side pressure should be near the 10 MPa Dirichlet BC";

    if (rank_ == 0) std::remove(config_path.c_str());
}
