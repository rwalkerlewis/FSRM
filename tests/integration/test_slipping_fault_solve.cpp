/**
 * @file test_slipping_fault_solve.cpp
 * @brief Integration test: slipping fault with Coulomb friction through TSSolve
 *
 * Tests the slipping (friction-governed) fault mode through the full Simulator
 * pipeline including TSSolve. The test applies shear loading that exceeds the
 * Coulomb friction strength (tau_applied > mu_f * sigma_n).
 *
 * Uses the semi-smooth Newton tangent operator for Coulomb friction in
 * addCohesivePenaltyToJacobian, which provides the correct derivatives:
 *   - d(constraint)/d(u): tangential slip direction derivative
 *   - d(constraint)/d(lambda): friction strength dependency on normal traction
 *
 * The test verifies:
 *   1. Setup pipeline completes (mesh split, cohesive cells, BCs)
 *   2. TSSolve converges
 *   3. Tangential slip is nonzero (fault slides under applied shear)
 *   4. Normal opening is near zero (compressive loading)
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

class SlippingFaultSolveTest : public ::testing::Test
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
    // Quasi-static solve: 1 timestep of 1 s
    // Grid: 1m x 1m x 1m box, 4x4x4 hexes
    // Fault: vertical YZ plane at x=0.5 (strike=90, dip=90)
    // BCs: bottom fixed, sides roller, top with shear traction
    // Friction: mu_f = 0.3
    // Loading: compression -10 MPa normal (z), shear 5 MPa (x on top)
    // Friction strength: tau_f = 0.3 * 10 MPa = 3 MPa
    // Applied shear: 5 MPa > 3 MPa => fault must slip in x

    cfg << "[SIMULATION]\n";
    cfg << "name = test_slipping_fault\n";
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
    cfg << "friction_coefficient = 0.3\n";
    cfg << "\n[BOUNDARY_CONDITIONS]\n";
    cfg << "bottom = fixed\n";
    cfg << "sides = roller\n";
    cfg << "top = compression\n";
    cfg << "\n[TRACTION_BC]\n";
    cfg << "enabled = true\n";
    cfg << "top_traction_x = 5.0e6\n";
    cfg << "top_traction_y = 0.0\n";
    cfg << "top_traction_z = -10.0e6\n";
    cfg << "\n[ABSORBING_BC]\n";
    cfg << "enabled = false\n";
    cfg.close();
  }

  PetscErrorCode runSlippingFault(const std::string& config_path,
                                   PetscReal& sol_norm,
                                   PetscReal& max_tangential_slip,
                                   PetscReal& max_normal_opening)
  {
    Simulator sim(PETSC_COMM_WORLD);
    PetscErrorCode ierr;

    PetscOptionsClear(nullptr);
    // Use the analytical (approximate) Jacobian with augmented Lagrangian
    // regularization -- addCohesivePenaltyToJacobian now contributes for
    // slipping mode with reduced penalty factor.
    PetscOptionsSetValue(nullptr, "-snes_max_it", "500");
    PetscOptionsSetValue(nullptr, "-snes_rtol", "1e-4");
    PetscOptionsSetValue(nullptr, "-snes_atol", "1e-2");
    PetscOptionsSetValue(nullptr, "-snes_stol", "1e-12");
    PetscOptionsSetValue(nullptr, "-ts_type", "beuler");
    PetscOptionsSetValue(nullptr, "-ts_max_snes_failures", "50");
    PetscOptionsSetValue(nullptr, "-ts_adapt_type", "basic");
    PetscOptionsSetValue(nullptr, "-pc_type", "lu");
    PetscOptionsSetValue(nullptr, "-ksp_type", "preonly");
    PetscOptionsSetValue(nullptr, "-snes_linesearch_type", "bt");
    PetscOptionsSetValue(nullptr, "-snes_linesearch_maxstep", "1e10");
    PetscOptionsSetValue(nullptr, "-snes_monitor", nullptr);
    PetscOptionsSetValue(nullptr, "-snes_converged_reason", nullptr);

    ierr = sim.initializeFromConfigFile(config_path);
    if (ierr) return ierr;
    ierr = sim.setupDM();
    if (ierr) return ierr;
    ierr = sim.setupFaultNetwork();
    if (ierr) return ierr;
    ierr = sim.labelBoundaries();
    if (ierr) return ierr;
    ierr = sim.setupFields();
    if (ierr) return ierr;
    ierr = sim.setupPhysics();
    if (ierr) return ierr;
    ierr = sim.setupTimeStepper();
    if (ierr) return ierr;
    ierr = sim.setupSolvers();
    if (ierr) return ierr;
    ierr = sim.setInitialConditions();
    if (ierr) return ierr;

    PetscPushErrorHandler(PetscReturnErrorHandler, nullptr);
    ierr = sim.run();
    PetscPopErrorHandler();
    if (ierr) return ierr;

    Vec sol = sim.getSolution();
    if (!sol)
    {
      sol_norm = -1.0;
      max_tangential_slip = -1.0;
      max_normal_opening = -1.0;
      return 0;
    }

    VecNorm(sol, NORM_2, &sol_norm);

    // Compute fault slip: tangential and normal components
    DM dm = sim.getDM();
    PetscSection section = nullptr;
    PetscSection coord_section = nullptr;
    DMLabel depth_label = nullptr;
    Vec local_solution = nullptr;
    Vec coords = nullptr;
    PetscInt dim = 0;

    ierr = DMGetLocalSection(dm, &section); CHKERRQ(ierr);
    ierr = DMGetCoordinateSection(dm, &coord_section); CHKERRQ(ierr);
    ierr = DMPlexGetDepthLabel(dm, &depth_label); CHKERRQ(ierr);
    ierr = DMGetDimension(dm, &dim); CHKERRQ(ierr);
    ierr = DMGetCoordinatesLocal(dm, &coords); CHKERRQ(ierr);
    if (!coords) { ierr = DMGetCoordinates(dm, &coords); CHKERRQ(ierr); }

    ierr = DMGetLocalVector(dm, &local_solution); CHKERRQ(ierr);
    ierr = VecZeroEntries(local_solution); CHKERRQ(ierr);
    ierr = DMGlobalToLocal(dm, sol, INSERT_VALUES, local_solution); CHKERRQ(ierr);

    const PetscScalar* uarray = nullptr;
    const PetscScalar* coord_array = nullptr;
    ierr = VecGetArrayRead(local_solution, &uarray); CHKERRQ(ierr);
    ierr = VecGetArrayRead(coords, &coord_array); CHKERRQ(ierr);

    max_tangential_slip = 0.0;
    max_normal_opening = 0.0;

    // Fault normal: x-direction (strike=90, dip=90 means fault is YZ plane at x=0.5)
    const PetscReal fault_normal[3] = {1.0, 0.0, 0.0};

    PetscInt cStart = 0, cEnd = 0;
    ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd); CHKERRQ(ierr);

    for (PetscInt c = cStart; c < cEnd; ++c)
    {
      DMPolytopeType ct;
      ierr = DMPlexGetCellType(dm, c, &ct); CHKERRQ(ierr);
      if (ct != DM_POLYTOPE_SEG_PRISM_TENSOR &&
          ct != DM_POLYTOPE_TRI_PRISM_TENSOR &&
          ct != DM_POLYTOPE_QUAD_PRISM_TENSOR)
      {
        continue;
      }

      const PetscInt* cone = nullptr;
      PetscInt cone_size = 0;
      ierr = DMPlexGetConeSize(dm, c, &cone_size); CHKERRQ(ierr);
      if (cone_size < 2) continue;
      ierr = DMPlexGetCone(dm, c, &cone); CHKERRQ(ierr);

      // Collect vertices from negative and positive faces
      auto collectVertices = [&](PetscInt face, std::vector<std::pair<PetscInt, std::array<PetscReal,3>>>& verts)
          -> PetscErrorCode
      {
        PetscInt closure_size = 0;
        PetscInt* closure = nullptr;
        PetscErrorCode e = DMPlexGetTransitiveClosure(dm, face, PETSC_TRUE, &closure_size, &closure);
        CHKERRQ(e);
        verts.clear();
        for (PetscInt i = 0; i < closure_size; ++i)
        {
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
      ierr = collectVertices(cone[0], neg_verts); CHKERRQ(ierr);
      ierr = collectVertices(cone[1], pos_verts); CHKERRQ(ierr);

      PetscInt n_pairs = static_cast<PetscInt>(std::min(neg_verts.size(), pos_verts.size()));
      for (PetscInt i = 0; i < n_pairs; ++i)
      {
        // Check vertices are at the same position
        PetscReal dist2 = 0;
        for (PetscInt d = 0; d < dim; ++d)
        {
          PetscReal delta = neg_verts[static_cast<std::size_t>(i)].second[static_cast<std::size_t>(d)]
                          - pos_verts[static_cast<std::size_t>(i)].second[static_cast<std::size_t>(d)];
          dist2 += delta * delta;
        }
        if (dist2 > 1.0e-20) continue;

        // Read displacements
        PetscReal u_neg[3] = {0,0,0}, u_pos[3] = {0,0,0};
        PetscInt dof_n = 0, off_n = 0, dof_p = 0, off_p = 0;
        ierr = PetscSectionGetFieldDof(section, neg_verts[static_cast<std::size_t>(i)].first, 0, &dof_n); CHKERRQ(ierr);
        if (dof_n > 0) {
          ierr = PetscSectionGetFieldOffset(section, neg_verts[static_cast<std::size_t>(i)].first, 0, &off_n); CHKERRQ(ierr);
          for (PetscInt d = 0; d < PetscMin(dof_n, dim); ++d) u_neg[d] = PetscRealPart(uarray[off_n+d]);
        }
        ierr = PetscSectionGetFieldDof(section, pos_verts[static_cast<std::size_t>(i)].first, 0, &dof_p); CHKERRQ(ierr);
        if (dof_p > 0) {
          ierr = PetscSectionGetFieldOffset(section, pos_verts[static_cast<std::size_t>(i)].first, 0, &off_p); CHKERRQ(ierr);
          for (PetscInt d = 0; d < PetscMin(dof_p, dim); ++d) u_pos[d] = PetscRealPart(uarray[off_p+d]);
        }

        // Compute slip = u+ - u-
        PetscReal slip[3];
        for (PetscInt d = 0; d < dim; ++d) slip[d] = u_pos[d] - u_neg[d];

        // Normal component
        PetscReal slip_n = 0;
        for (PetscInt d = 0; d < dim; ++d) slip_n += slip[d] * fault_normal[d];

        // Tangential component
        PetscReal slip_t[3];
        PetscReal slip_t_mag2 = 0;
        for (PetscInt d = 0; d < dim; ++d)
        {
          slip_t[d] = slip[d] - slip_n * fault_normal[d];
          slip_t_mag2 += slip_t[d] * slip_t[d];
        }

        max_tangential_slip = PetscMax(max_tangential_slip, PetscSqrtReal(slip_t_mag2));
        max_normal_opening = PetscMax(max_normal_opening, PetscAbsReal(slip_n));
      }
    }

    ierr = VecRestoreArrayRead(coords, &coord_array); CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(local_solution, &uarray); CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(dm, &local_solution); CHKERRQ(ierr);
    return PETSC_SUCCESS;
  }

  int rank_ = 0;
};

// Slipping fault setup test: verifies the full pipeline can be configured
// and the solve is attempted. Currently SNES diverges due to missing
// tangent operator for the nonlinear Coulomb friction constraint.
TEST_F(SlippingFaultSolveTest, SlippingFaultQuasiStatic)
{
  std::string config_path = "test_slipping_fault_solve.config";
  writeConfig(config_path);
  MPI_Barrier(PETSC_COMM_WORLD);

  PetscReal sol_norm = 0.0;
  PetscReal max_tangential_slip = 0.0;
  PetscReal max_normal_opening = 0.0;

  PetscErrorCode ierr = runSlippingFault(config_path, sol_norm,
                                          max_tangential_slip, max_normal_opening);

  if (rank_ == 0) std::remove(config_path.c_str());

  ASSERT_EQ(ierr, 0) << "Slipping fault TSSolve must converge with semi-smooth Newton Jacobian";

  // If the solver converges in the future, verify the solution
  EXPECT_GT(sol_norm, 0.0) << "Solution must be nonzero";
  EXPECT_TRUE(std::isfinite(sol_norm)) << "Solution norm must be finite";

  EXPECT_GT(max_tangential_slip, 1.0e-8)
      << "Tangential fault slip must be nonzero when tau_applied > mu_f * sigma_n";

  EXPECT_LT(max_normal_opening, 1.0e-4)
      << "Normal opening should be near zero under compressive loading";
}
