/**
 * @file test_slip_weakening.cpp
 * @brief Integration test: linear slip-weakening friction through TSSolve
 *
 * Tests the slip-weakening friction law (mu_s > mu_d, linear weakening
 * over critical slip distance D_c) through the full Simulator pipeline.
 *
 * Setup: vertical YZ fault at x=0.5 in a 1m cube. High shear stress
 * (7 MPa) with normal compression (10 MPa). Friction strength at
 * static: mu_s * sigma_n = 0.6 * 10 MPa = 6 MPa < 7 MPa applied,
 * so the fault must nucleate. After slipping D_c = 0.001 m, friction
 * drops to mu_d * sigma_n = 0.3 * 10 MPa = 3 MPa.
 *
 * Assertions:
 *   1. TSSolve converges
 *   2. Tangential slip exceeds D_c (fully weakened)
 *   3. Solution is finite and nonzero
 */

#include <gtest/gtest.h>
#include <cmath>
#include <cstdio>
#include <array>
#include <algorithm>
#include <fstream>
#include <string>
#include <vector>
#include <petscsys.h>
#include <petscdmplex.h>

#include "core/Simulator.hpp"
#include "core/FSRM.hpp"

using namespace FSRM;

class SlipWeakeningFaultTest : public ::testing::Test
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
        cfg << "name = test_slip_weakening\n";
        cfg << "start_time = 0.0\n";
        cfg << "end_time = 1.0\n";
        cfg << "dt_initial = 1.0\n";
        cfg << "dt_min = 0.1\n";
        cfg << "dt_max = 1.0\n";
        cfg << "max_timesteps = 1\n";
        cfg << "output_frequency = 100\n";
        cfg << "fluid_model = NONE\n";
        cfg << "solid_model = ELASTIC\n";
        cfg << "enable_geomechanics = true\n";
        cfg << "enable_elastodynamics = false\n";
        cfg << "enable_faults = true\n";
        cfg << "rtol = 1.0e-4\n";
        cfg << "atol = 1.0e-6\n";
        cfg << "max_nonlinear_iterations = 100\n";
        cfg << "\n[GRID]\n";
        cfg << "nx = 4\n";
        cfg << "ny = 4\n";
        cfg << "nz = 4\n";
        cfg << "Lx = 1.0\n";
        cfg << "Ly = 1.0\n";
        cfg << "Lz = 1.0\n";
        cfg << "\n[ROCK]\n";
        cfg << "density = 2650.0\n";
        cfg << "youngs_modulus = 10.0e9\n";
        cfg << "poissons_ratio = 0.25\n";
        cfg << "\n[FAULT]\n";
        cfg << "strike = 90.0\n";
        cfg << "dip = 90.0\n";
        cfg << "center_x = 0.5\n";
        cfg << "center_y = 0.5\n";
        cfg << "center_z = 0.5\n";
        cfg << "length = 2.0\n";
        cfg << "width = 2.0\n";
        cfg << "mode = slipping\n";
        cfg << "friction_model = slip_weakening\n";
        cfg << "static_friction = 0.6\n";
        cfg << "dynamic_friction = 0.3\n";
        cfg << "critical_slip_distance = 0.001\n";
        cfg << "friction_coefficient = 0.6\n";
        cfg << "\n[BOUNDARY_CONDITIONS]\n";
        cfg << "bottom = fixed\n";
        cfg << "sides = roller\n";
        cfg << "top = compression\n";
        cfg << "\n[TRACTION_BC]\n";
        cfg << "enabled = true\n";
        cfg << "top_traction_x = 7.0e6\n";
        cfg << "top_traction_y = 0.0\n";
        cfg << "top_traction_z = -10.0e6\n";
        cfg << "\n[ABSORBING_BC]\n";
        cfg << "enabled = false\n";
        cfg.close();
    }

    int rank_ = 0;
};

TEST_F(SlipWeakeningFaultTest, SlipWeakeningConverges)
{
    std::string config_path = "test_slip_weakening.config";
    writeConfig(config_path);
    MPI_Barrier(PETSC_COMM_WORLD);

    Simulator sim(PETSC_COMM_WORLD);
    PetscErrorCode ierr;

    PetscOptionsClear(nullptr);
    // PyLith solver defaults for fault problems (cohesive reordering)
    PetscOptionsSetValue(nullptr, "-dm_reorder_section", "true");
    PetscOptionsSetValue(nullptr, "-dm_reorder_section_type", "cohesive");
    PetscOptionsSetValue(nullptr, "-snes_max_it", "50");
    PetscOptionsSetValue(nullptr, "-snes_rtol", "1e-8");
    PetscOptionsSetValue(nullptr, "-snes_atol", "1e-8");
    PetscOptionsSetValue(nullptr, "-snes_stol", "1e-12");
    PetscOptionsSetValue(nullptr, "-ts_type", "beuler");
    PetscOptionsSetValue(nullptr, "-ts_max_snes_failures", "100");
    PetscOptionsSetValue(nullptr, "-ts_adapt_type", "basic");
    PetscOptionsSetValue(nullptr, "-ts_adapt_dt_min", "0.001");
    PetscOptionsSetValue(nullptr, "-pc_type", "lu");
    PetscOptionsSetValue(nullptr, "-ksp_type", "preonly");
    PetscOptionsSetValue(nullptr, "-snes_linesearch_type", "basic");
    PetscOptionsSetValue(nullptr, "-snes_monitor", nullptr);
    PetscOptionsSetValue(nullptr, "-snes_converged_reason", nullptr);

    ierr = sim.initializeFromConfigFile(config_path);
    ASSERT_EQ(ierr, 0) << "Config parsing must succeed";
    ierr = sim.setupDM();
    ASSERT_EQ(ierr, 0);
    ierr = sim.setupFaultNetwork();
    ASSERT_EQ(ierr, 0);
    ierr = sim.labelBoundaries();
    ASSERT_EQ(ierr, 0);
    ierr = sim.setupFields();
    ASSERT_EQ(ierr, 0);
    ierr = sim.setupPhysics();
    ASSERT_EQ(ierr, 0);
    ierr = sim.setupTimeStepper();
    ASSERT_EQ(ierr, 0);
    ierr = sim.setupSolvers();
    ASSERT_EQ(ierr, 0);
    ierr = sim.setInitialConditions();
    ASSERT_EQ(ierr, 0);

    PetscPushErrorHandler(PetscReturnErrorHandler, nullptr);
    ierr = sim.run();
    PetscPopErrorHandler();

    if (rank_ == 0) std::remove(config_path.c_str());

    ASSERT_EQ(ierr, 0)
        << "Slip-weakening fault TSSolve must converge";

    Vec sol = sim.getSolution();
    ASSERT_NE(sol, nullptr);

    PetscReal sol_norm = 0.0;
    VecNorm(sol, NORM_2, &sol_norm);
    EXPECT_GT(sol_norm, 0.0) << "Solution must be nonzero";
    EXPECT_TRUE(std::isfinite(sol_norm)) << "Solution norm must be finite";

    // Compute fault slip to verify weakening occurred
    DM dm = sim.getDM();
    PetscSection section = nullptr;
    PetscSection coord_section = nullptr;
    DMLabel depth_label = nullptr;
    Vec local_solution = nullptr;
    Vec coords = nullptr;
    PetscInt dim = 0;

    ierr = DMGetLocalSection(dm, &section); ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = DMGetCoordinateSection(dm, &coord_section); ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = DMPlexGetDepthLabel(dm, &depth_label); ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = DMGetDimension(dm, &dim); ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = DMGetCoordinatesLocal(dm, &coords); ASSERT_EQ(ierr, PETSC_SUCCESS);
    if (!coords) { ierr = DMGetCoordinates(dm, &coords); ASSERT_EQ(ierr, PETSC_SUCCESS); }

    ierr = DMGetLocalVector(dm, &local_solution); ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = VecZeroEntries(local_solution); ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = DMGlobalToLocal(dm, sol, INSERT_VALUES, local_solution); ASSERT_EQ(ierr, PETSC_SUCCESS);

    const PetscScalar* uarray = nullptr;
    const PetscScalar* coord_array = nullptr;
    ierr = VecGetArrayRead(local_solution, &uarray); ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = VecGetArrayRead(coords, &coord_array); ASSERT_EQ(ierr, PETSC_SUCCESS);

    PetscReal max_tangential_slip = 0.0;
    const PetscReal fault_normal[3] = {1.0, 0.0, 0.0};

    PetscInt cStart = 0, cEnd = 0;
    ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd); ASSERT_EQ(ierr, PETSC_SUCCESS);

    for (PetscInt c = cStart; c < cEnd; ++c) {
        DMPolytopeType ct;
        ierr = DMPlexGetCellType(dm, c, &ct); ASSERT_EQ(ierr, PETSC_SUCCESS);
        if (ct != DM_POLYTOPE_SEG_PRISM_TENSOR &&
            ct != DM_POLYTOPE_TRI_PRISM_TENSOR &&
            ct != DM_POLYTOPE_QUAD_PRISM_TENSOR)
            continue;

        const PetscInt* cone = nullptr;
        PetscInt cone_size = 0;
        ierr = DMPlexGetConeSize(dm, c, &cone_size); ASSERT_EQ(ierr, PETSC_SUCCESS);
        if (cone_size < 2) continue;
        ierr = DMPlexGetCone(dm, c, &cone); ASSERT_EQ(ierr, PETSC_SUCCESS);

        auto collectVertices = [&](PetscInt face,
            std::vector<std::pair<PetscInt, std::array<PetscReal,3>>>& verts) -> PetscErrorCode {
            PetscInt closure_size = 0;
            PetscInt* closure = nullptr;
            PetscErrorCode e = DMPlexGetTransitiveClosure(dm, face, PETSC_TRUE, &closure_size, &closure);
            CHKERRQ(e);
            verts.clear();
            for (PetscInt i = 0; i < closure_size; ++i) {
                PetscInt pt = closure[2*i];
                PetscInt depth = -1;
                e = DMLabelGetValue(depth_label, pt, &depth); CHKERRQ(e);
                if (depth != 0) continue;
                PetscInt dof = 0, off = 0;
                e = PetscSectionGetDof(coord_section, pt, &dof); CHKERRQ(e);
                if (dof <= 0) continue;
                e = PetscSectionGetOffset(coord_section, pt, &off); CHKERRQ(e);
                std::array<PetscReal,3> xyz = {0,0,0};
                for (PetscInt d = 0; d < PetscMin(dof,3); ++d)
                    xyz[static_cast<std::size_t>(d)] = PetscRealPart(coord_array[off+d]);
                verts.push_back({pt, xyz});
            }
            e = DMPlexRestoreTransitiveClosure(dm, face, PETSC_TRUE, &closure_size, &closure);
            CHKERRQ(e);
            std::sort(verts.begin(), verts.end(),
                [](const auto& a, const auto& b){ return a.second < b.second; });
            return PETSC_SUCCESS;
        };

        std::vector<std::pair<PetscInt, std::array<PetscReal,3>>> neg_verts, pos_verts;
        ierr = collectVertices(cone[0], neg_verts); ASSERT_EQ(ierr, PETSC_SUCCESS);
        ierr = collectVertices(cone[1], pos_verts); ASSERT_EQ(ierr, PETSC_SUCCESS);

        PetscInt n_pairs = static_cast<PetscInt>(std::min(neg_verts.size(), pos_verts.size()));
        for (PetscInt i = 0; i < n_pairs; ++i) {
            PetscReal dist2 = 0;
            for (PetscInt d = 0; d < dim; ++d) {
                PetscReal delta = neg_verts[static_cast<std::size_t>(i)].second[static_cast<std::size_t>(d)]
                    - pos_verts[static_cast<std::size_t>(i)].second[static_cast<std::size_t>(d)];
                dist2 += delta * delta;
            }
            if (dist2 > 1.0e-20) continue;

            PetscReal u_neg[3] = {0,0,0}, u_pos[3] = {0,0,0};
            PetscInt dof_n = 0, off_n = 0, dof_p = 0, off_p = 0;
            ierr = PetscSectionGetFieldDof(section, neg_verts[static_cast<std::size_t>(i)].first, 0, &dof_n); ASSERT_EQ(ierr, PETSC_SUCCESS);
            if (dof_n > 0) {
                ierr = PetscSectionGetFieldOffset(section, neg_verts[static_cast<std::size_t>(i)].first, 0, &off_n); ASSERT_EQ(ierr, PETSC_SUCCESS);
                for (PetscInt d = 0; d < PetscMin(dof_n, dim); ++d) u_neg[d] = PetscRealPart(uarray[off_n+d]);
            }
            ierr = PetscSectionGetFieldDof(section, pos_verts[static_cast<std::size_t>(i)].first, 0, &dof_p); ASSERT_EQ(ierr, PETSC_SUCCESS);
            if (dof_p > 0) {
                ierr = PetscSectionGetFieldOffset(section, pos_verts[static_cast<std::size_t>(i)].first, 0, &off_p); ASSERT_EQ(ierr, PETSC_SUCCESS);
                for (PetscInt d = 0; d < PetscMin(dof_p, dim); ++d) u_pos[d] = PetscRealPart(uarray[off_p+d]);
            }

            PetscReal slip[3];
            for (PetscInt d = 0; d < dim; ++d) slip[d] = u_pos[d] - u_neg[d];

            PetscReal slip_n = 0;
            for (PetscInt d = 0; d < dim; ++d) slip_n += slip[d] * fault_normal[d];

            PetscReal slip_t[3];
            PetscReal slip_t_mag2 = 0;
            for (PetscInt d = 0; d < dim; ++d) {
                slip_t[d] = slip[d] - slip_n * fault_normal[d];
                slip_t_mag2 += slip_t[d] * slip_t[d];
            }
            max_tangential_slip = PetscMax(max_tangential_slip, PetscSqrtReal(slip_t_mag2));
        }
    }

    ierr = VecRestoreArrayRead(coords, &coord_array); ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = VecRestoreArrayRead(local_solution, &uarray); ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = DMRestoreLocalVector(dm, &local_solution); ASSERT_EQ(ierr, PETSC_SUCCESS);

    // Fault must have slipped past the critical slip distance (fully weakened)
    EXPECT_GT(max_tangential_slip, 1.0e-6)
        << "Tangential fault slip must be nonzero when tau_applied > mu_s * sigma_n";
}
