/**
 * @file test_source_ball_radiation_3d.cpp
 * @brief Pass-13b unit tests for the cell-centred FV grey radiation
 *        diffusion solver on the 3D source ball.
 *
 * Hand-built minimal mesh: a chain of N cells of unit volume separated
 * by unit-area faces of unit centroid spacing. The chain models a 1D
 * slice of the 3D radial layout; that is enough to exercise the FV
 * assembly, KSP solve, and Newton outer loop without depending on a
 * DMPlex mesh fixture.
 */

#include <gtest/gtest.h>

#include <petscsys.h>

#include <array>
#include <cmath>
#include <vector>

#include "domain/explosion/SourceBallRadiation3D.hpp"
#include "domain/explosion/OpacityModel.hpp"

using FSRM::BoundaryFace3D;
using FSRM::BoundaryFaceType3D;
using FSRM::CellGeometry3D;
using FSRM::InternalFace3D;
using FSRM::SourceBallRadiation3DSolver;
using FSRM::PowerLawOpacitySets;

namespace
{

struct ChainMesh {
    std::vector<CellGeometry3D> cells;
    std::vector<InternalFace3D> internal_faces;
    std::vector<BoundaryFace3D> boundary_faces;
};

ChainMesh makeChain(int N, double cell_volume = 1.0,
                    double face_area = 1.0, double face_dist = 1.0)
{
    ChainMesh m;
    m.cells.resize(N);
    for (int i = 0; i < N; ++i) {
        m.cells[i].volume = cell_volume;
        m.cells[i].centroid = {static_cast<double>(i) * face_dist, 0.0, 0.0};
    }
    for (int i = 0; i < N - 1; ++i) {
        InternalFace3D f;
        f.cell_left_local = i;
        f.cell_right_local = i + 1;
        f.cell_left_global = i;
        f.cell_right_global = i + 1;
        f.area = face_area;
        f.dist_centroid = face_dist;
        m.internal_faces.push_back(f);
    }
    // Two boundary faces (vacuum) at the ends. For the
    // EnergyConservationNoCoupling gate we make these reflective by
    // omitting them; the gate requires no flux across boundaries and
    // no matter coupling.
    return m;
}

}  // namespace

TEST(SourceBallRadiation3D, UniformFieldStable)
{
    ChainMesh m = makeChain(8);

    SourceBallRadiation3DSolver::Config cfg;
    // Tiny constant opacity -> diffusion coefficient large but
    // matter-coupling source c*kP*rho is bounded (kP * rho ~ 1e-15).
    // The diffusion terms cancel exactly for a uniform initial field
    // (zero face flux), so the solution stays uniform regardless of D.
    cfg.kappa_constant_m2_per_kg = 1.0e-3;
    cfg.cv_J_per_kg_K = 1000.0;
    cfg.max_newton_iter = 5;
    cfg.newton_tolerance = 1.0e-8;

    SourceBallRadiation3DSolver solver;
    solver.setConfig(cfg);
    solver.initialize(PETSC_COMM_SELF, 8, 0, 8, m.cells, m.internal_faces,
                      m.boundary_faces);

    // Uniform IC, no boundary faces. With matter at radiative
    // equilibrium (T_m chosen so a T_m^4 = E_r), the source term
    // vanishes. Combined with zero face flux -> the field is exact
    // steady state.
    const double E0 = 100.0;
    const double T_eq = std::pow(E0 / FSRM::RadiationConstants::
                                     RADIATION_CONSTANT_A_J_PER_M3_K4,
                                 0.25);
    std::vector<double> rho(8, 1.0);
    std::vector<double> T_m(8, T_eq);
    std::vector<double> E_r(8, E0);
    std::vector<double> e_int(8, 0.0);

    auto res = solver.step(1.0e-9, rho, T_m, E_r, e_int);
    EXPECT_TRUE(res.converged) << "iters=" << res.newton_iters
                                << " resid=" << res.residual_inf_norm;

    for (int i = 0; i < 8; ++i) {
        EXPECT_NEAR(E_r[i], E0, 1.0e-6 * E0)
            << "Uniform E_r should remain uniform; cell " << i
            << " drifted to " << E_r[i];
    }
}

TEST(SourceBallRadiation3D, EnergyConservationNoCoupling)
{
    ChainMesh m = makeChain(16);

    SourceBallRadiation3DSolver::Config cfg;
    cfg.kappa_constant_m2_per_kg = 1.0e-3;  // small constant kappa
    cfg.cv_J_per_kg_K = 1000.0;
    cfg.max_newton_iter = 5;
    cfg.newton_tolerance = 1.0e-10;

    SourceBallRadiation3DSolver solver;
    solver.setConfig(cfg);
    solver.initialize(PETSC_COMM_SELF, 16, 0, 16, m.cells, m.internal_faces,
                      m.boundary_faces);

    // Near-zero matter (kappa_P * rho is non-zero through constant
    // kappa, but a T^4 ~ 0 so the source is just absorption of E_r).
    // To suppress coupling completely set kP * rho * V * dt = 0 by
    // setting rho = 0. The diffusion coefficient D = c / (3 kappa rho)
    // diverges in that limit; instead use a tiny rho and a tiny dt so
    // the matter coupling is negligible.
    std::vector<double> rho(16, 1.0e-3);
    std::vector<double> T_m(16, 1.0);
    std::vector<double> E_r(16, 0.0);
    // Localized initial pulse at cell 7 / 8.
    E_r[7] = 100.0;
    E_r[8] = 100.0;
    std::vector<double> e_int(16, 0.0);

    double E_initial = 0.0;
    for (int i = 0; i < 16; ++i) E_initial += E_r[i] * m.cells[i].volume;

    const double dt = 1.0e-9;  // very small to suppress matter coupling.
    for (int step = 0; step < 50; ++step) {
        auto res = solver.step(dt, rho, T_m, E_r, e_int);
        ASSERT_TRUE(res.converged) << "step " << step;
    }

    double E_final = 0.0;
    for (int i = 0; i < 16; ++i) E_final += E_r[i] * m.cells[i].volume;

    // No boundary faces -> no outflow. Matter coupling at this dt is
    // negligible. Total energy is conserved within Newton tolerance.
    const double rel = std::abs(E_final - E_initial) / E_initial;
    EXPECT_LT(rel, 1.0e-3)
        << "E_initial=" << E_initial << " E_final=" << E_final
        << " rel diff=" << rel;
}

TEST(SourceBallRadiation3D, EquilibrationToMatter)
{
    // Single cell, no boundary faces, no internal faces. Only the
    // matter-coupling term acts. Should drive E_r toward a T^4 (LTE).
    ChainMesh m = makeChain(1);

    SourceBallRadiation3DSolver::Config cfg;
    cfg.kappa_constant_m2_per_kg = 1.0;  // strong coupling
    cfg.cv_J_per_kg_K = 1.0e6;            // large heat capacity:
                                          // matter T won't drift much
    cfg.max_newton_iter = 5;
    cfg.newton_tolerance = 1.0e-8;

    SourceBallRadiation3DSolver solver;
    solver.setConfig(cfg);
    solver.initialize(PETSC_COMM_SELF, 1, 0, 1, m.cells, m.internal_faces,
                      m.boundary_faces);

    // Hot matter (T_m=1e4 K), cold radiation (E_r=0).
    std::vector<double> rho(1, 1.0);
    std::vector<double> T_m(1, 1.0e4);
    std::vector<double> E_r(1, 0.0);
    std::vector<double> e_int(1, 0.0);
    const double T_target = T_m[0];
    const double aT4_target =
        FSRM::RadiationConstants::RADIATION_CONSTANT_A_J_PER_M3_K4
        * std::pow(T_target, 4.0);

    const double dt = 1.0e-6;
    for (int step = 0; step < 200; ++step) {
        auto res = solver.step(dt, rho, T_m, E_r, e_int);
        ASSERT_TRUE(res.converged);
    }
    // E_r approaches a T_m^4 within 1% (with cv chosen large enough that
    // T_m barely drifts).
    const double rel = std::abs(E_r[0] - aT4_target) / aT4_target;
    EXPECT_LT(rel, 0.01)
        << "E_r=" << E_r[0] << " target a T_m^4=" << aT4_target
        << " rel diff=" << rel;
}
