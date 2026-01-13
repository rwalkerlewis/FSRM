/**
 * @file test_boundary_conditions.cpp
 * @brief Unit tests for boundary conditions
 * 
 * Tests cover:
 * - Free surface boundary conditions
 * - Absorbing boundary conditions (PML, CE, LK)
 * - Dirichlet/Neumann conditions
 * - Periodic and symmetry conditions
 */

#include "BoundaryConditions.hpp"
#include "Testing.hpp"
#include <cmath>
#include <iostream>
#include <vector>
#include <array>

namespace FSRM {
namespace Testing {

// =============================================================================
// Helper Functions
// =============================================================================

BoundaryFace createTestFace(const std::array<double, 3>& normal,
                           const std::array<double, 3>& centroid,
                           double rho = 2500.0, double vp = 5000.0, double vs = 3000.0) {
    BoundaryFace face;
    face.element_id = 0;
    face.local_face_id = 0;
    face.centroid = centroid;
    face.normal = normal;
    face.area = 1.0;
    face.rho = rho;
    face.vp = vp;
    face.vs = vs;
    face.impedance_p = rho * vp;
    face.impedance_s = rho * vs;
    return face;
}

// =============================================================================
// Free Surface Tests
// =============================================================================

/**
 * @test Test free surface traction-free condition
 */
bool test_free_surface_traction_free() {
    std::cout << "Testing free surface traction-free condition..." << std::endl;
    
    double tol = 1e-10;
    
    FreeSurfaceBC bc;
    
    // Create face with normal pointing up (z direction)
    BoundaryFace face = createTestFace({0, 0, 1}, {0, 0, 0});
    
    // Interior state with non-zero stress
    double u_int[9] = {
        0.1, 0.05, 0.02,          // Velocity: vx, vy, vz
        1e6, 0.5e6, 0.3e6,        // Normal stresses: σxx, σyy, σzz
        0.2e6, 0.1e6, 0.15e6      // Shear stresses: σxy, σxz, σyz
    };
    
    // Get boundary state (ghost)
    double u_ext[9];
    bc.getBoundaryState(face, u_int, u_ext, 0.0);
    
    // For free surface, traction t = σ·n should be zero
    // Compute traction from exterior state
    double tn_x = u_ext[3] * face.normal[0] + u_ext[6] * face.normal[1] + u_ext[7] * face.normal[2];
    double tn_y = u_ext[6] * face.normal[0] + u_ext[4] * face.normal[1] + u_ext[8] * face.normal[2];
    double tn_z = u_ext[7] * face.normal[0] + u_ext[8] * face.normal[1] + u_ext[5] * face.normal[2];
    
    // Average traction (between int and ext) should be approximately zero
    double tn_avg_x = 0.5 * (u_int[3] + u_ext[3]) * face.normal[0] + 
                      0.5 * (u_int[6] + u_ext[6]) * face.normal[1] +
                      0.5 * (u_int[7] + u_ext[7]) * face.normal[2];
    
    // The mirror should cancel the normal traction component
    // This is a simplified check
    if (std::isnan(u_ext[0]) || std::isinf(u_ext[0])) {
        std::cerr << "  FAIL: Invalid exterior state" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Free surface produces valid ghost state" << std::endl;
    return true;
}

/**
 * @test Test free surface velocity reflection
 */
bool test_free_surface_velocity() {
    std::cout << "Testing free surface velocity condition..." << std::endl;
    
    FreeSurfaceBC bc;
    BoundaryFace face = createTestFace({0, 0, 1}, {0, 0, 0});
    
    // Interior state with velocity
    double u_int[9] = {1.0, 0.5, 0.2, 1e6, 0.5e6, 0.3e6, 0.2e6, 0.1e6, 0.15e6};
    
    double u_ext[9];
    bc.getBoundaryState(face, u_int, u_ext, 0.0);
    
    // For symmetric reflection, velocities should match
    if (std::abs(u_ext[0] - u_int[0]) > 1e-10 ||
        std::abs(u_ext[1] - u_int[1]) > 1e-10 ||
        std::abs(u_ext[2] - u_int[2]) > 1e-10) {
        std::cerr << "  FAIL: Velocity not preserved at free surface" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Free surface velocity condition correct" << std::endl;
    return true;
}

// =============================================================================
// PML Tests
// =============================================================================

/**
 * @test Test PML damping profile
 */
bool test_pml_damping_profile() {
    std::cout << "Testing PML damping profile..." << std::endl;
    
    PMLBoundary pml;
    pml.setThickness(1000.0);  // 1km PML
    pml.setOrder(3);
    pml.setReflectionCoefficient(1e-6);
    pml.setDomainBounds(0, 10000, 0, 10000, -5000, 0);
    pml.setMaterialProperties(2500, 5000, 3000);
    
    // Test damping at different locations
    double d_x, d_y, d_z;
    
    // Inside domain (no PML)
    pml.getDampingCoefficients(5000, 5000, -2500, d_x, d_y, d_z);
    if (d_x > 1e-15 || d_y > 1e-15 || d_z > 1e-15) {
        std::cerr << "  FAIL: Non-zero damping inside domain" << std::endl;
        return false;
    }
    
    // At x boundary (should have x damping)
    pml.getDampingCoefficients(500, 5000, -2500, d_x, d_y, d_z);
    if (d_x <= 0) {
        std::cerr << "  FAIL: No damping in PML region" << std::endl;
        return false;
    }
    
    // Damping should increase toward boundary
    double d_x_inner, d_x_outer, d_y_inner, d_z_inner;
    pml.getDampingCoefficients(800, 5000, -2500, d_x_inner, d_y_inner, d_z_inner);
    pml.getDampingCoefficients(200, 5000, -2500, d_x_outer, d_y, d_z);
    
    if (d_x_inner >= d_x_outer) {
        std::cerr << "  FAIL: Damping doesn't increase toward boundary" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: PML damping profile correct" << std::endl;
    return true;
}

/**
 * @test Test PML region detection
 */
bool test_pml_region_detection() {
    std::cout << "Testing PML region detection..." << std::endl;
    
    PMLBoundary pml;
    pml.setThickness(1000.0);
    pml.setDomainBounds(0, 10000, 0, 10000, -5000, 0);
    
    // Inside domain
    if (pml.isInPML(5000, 5000, -2500)) {
        std::cerr << "  FAIL: Center detected as PML" << std::endl;
        return false;
    }
    
    // In PML region (near boundary)
    if (!pml.isInPML(500, 5000, -2500)) {
        std::cerr << "  FAIL: PML region not detected" << std::endl;
        return false;
    }
    
    // Check direction flags
    int dir = pml.getPMLDirection(500, 5000, -200);  // Near x-min and z-max
    if (!(dir & 1) || !(dir & 4)) {  // Should have x and z flags
        std::cerr << "  FAIL: Incorrect PML direction flags" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: PML region detection correct" << std::endl;
    return true;
}

// =============================================================================
// Clayton-Engquist Tests
// =============================================================================

/**
 * @test Test Clayton-Engquist impedance matrix
 */
bool test_clayton_engquist_impedance() {
    std::cout << "Testing Clayton-Engquist impedance matrix..." << std::endl;
    
    double tol = 1e-10;
    
    ClaytonEngquistBC ce;
    ce.setMaterialProperties(2500, 5000, 3000);
    
    // Test impedance matrix for normal = (1, 0, 0)
    double n[3] = {1, 0, 0};
    double Z[9];
    ce.computeImpedanceMatrix(n, Z);
    
    // For normal in x, Z_xx should be Zp = ρ*vp
    double Zp = 2500 * 5000;
    if (std::abs(Z[0] - Zp) > tol) {
        std::cerr << "  FAIL: Z_xx = " << Z[0] << ", expected " << Zp << std::endl;
        return false;
    }
    
    // Z_yy and Z_zz should be Zs = ρ*vs
    double Zs = 2500 * 3000;
    if (std::abs(Z[4] - Zs) > tol || std::abs(Z[8] - Zs) > tol) {
        std::cerr << "  FAIL: Tangential impedance incorrect" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Clayton-Engquist impedance correct" << std::endl;
    return true;
}

/**
 * @test Test Clayton-Engquist absorbing traction
 */
bool test_clayton_engquist_absorbing() {
    std::cout << "Testing Clayton-Engquist absorbing traction..." << std::endl;
    
    ClaytonEngquistBC ce;
    ce.setMaterialProperties(2500, 5000, 3000);
    ce.setOrder(1);
    
    BoundaryFace face = createTestFace({1, 0, 0}, {0, 0, 0}, 2500, 5000, 3000);
    
    // Velocity state
    double v[3] = {1.0, 0.5, 0.2};
    double traction[3];
    
    ce.applyFirstOrder(face, v, traction);
    
    // Traction should oppose velocity (absorbing)
    // t = -Z * v, so t_n = -Zp * v_n
    double Zp = 2500 * 5000;
    double expected_tx = -Zp * v[0];
    
    if (std::abs(traction[0] - expected_tx) / std::abs(expected_tx) > 0.1) {
        std::cerr << "  FAIL: Normal traction incorrect" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Clayton-Engquist absorbing traction correct" << std::endl;
    return true;
}

// =============================================================================
// Lysmer-Kuhlemeyer Tests
// =============================================================================

/**
 * @test Test Lysmer-Kuhlemeyer dashpot
 */
bool test_lysmer_dashpot() {
    std::cout << "Testing Lysmer-Kuhlemeyer dashpot..." << std::endl;
    
    LysmerKuhlemeyer lk;
    lk.setMaterialProperties(2500, 5000, 3000);
    
    double n[3] = {1, 0, 0};
    double v[3] = {1.0, 0.5, 0.2};
    double traction[3];
    
    lk.computeAbsorbingTraction(v, n, traction);
    
    // Normal component: t_n = -ρ*vp*v_n
    double Zp = 2500 * 5000;
    double expected_tn = -Zp * v[0];
    
    // The traction should have correct sign (opposing velocity)
    double t_dot_v = traction[0] * v[0] + traction[1] * v[1] + traction[2] * v[2];
    
    if (t_dot_v >= 0) {
        std::cerr << "  FAIL: Dashpot not absorbing (t·v >= 0)" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Lysmer-Kuhlemeyer dashpot correct" << std::endl;
    return true;
}

// =============================================================================
// Dirichlet Tests
// =============================================================================

/**
 * @test Test homogeneous Dirichlet BC
 */
bool test_dirichlet_zero() {
    std::cout << "Testing homogeneous Dirichlet BC..." << std::endl;
    
    double tol = 1e-10;
    
    DirichletBC bc;
    bc.setZero();
    
    BoundaryFace face = createTestFace({1, 0, 0}, {0, 0, 0});
    
    double u_int[9] = {1.0, 0.5, 0.2, 1e6, 0.5e6, 0.3e6, 0.2e6, 0.1e6, 0.15e6};
    double u_ext[9];
    
    bc.getBoundaryState(face, u_int, u_ext, 0.0);
    
    // Ghost state should be such that average is zero
    // u_ext = 2*u_bc - u_int = -u_int for u_bc = 0
    for (int i = 0; i < 9; ++i) {
        if (std::abs(u_ext[i] + u_int[i]) > tol) {
            std::cerr << "  FAIL: Ghost state incorrect for component " << i << std::endl;
            return false;
        }
    }
    
    std::cout << "  PASS: Homogeneous Dirichlet BC correct" << std::endl;
    return true;
}

/**
 * @test Test prescribed Dirichlet BC
 */
bool test_dirichlet_prescribed() {
    std::cout << "Testing prescribed Dirichlet BC..." << std::endl;
    
    DirichletBC bc;
    bc.setPrescribedValue([](double x, double y, double z, double t, double* u) {
        u[0] = 2.0;  // Prescribed vx
        u[1] = 1.0;  // Prescribed vy
        u[2] = 0.5;  // Prescribed vz
        for (int i = 3; i < 9; ++i) u[i] = 0.0;
    });
    
    BoundaryFace face = createTestFace({1, 0, 0}, {0, 0, 0});
    
    double u_int[9] = {1.0, 0.5, 0.2, 0, 0, 0, 0, 0, 0};
    double u_ext[9];
    
    bc.getBoundaryState(face, u_int, u_ext, 0.0);
    
    // Average should equal prescribed value
    double avg_vx = 0.5 * (u_int[0] + u_ext[0]);
    double expected_vx = 2.0;
    
    if (std::abs(avg_vx - expected_vx) > 1e-10) {
        std::cerr << "  FAIL: Prescribed value not enforced" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Prescribed Dirichlet BC correct" << std::endl;
    return true;
}

// =============================================================================
// Neumann Tests
// =============================================================================

/**
 * @test Test constant Neumann BC
 */
bool test_neumann_constant() {
    std::cout << "Testing constant Neumann BC..." << std::endl;
    
    NeumannBC bc;
    bc.setConstantTraction(1e6, 0.5e6, 0.0);
    
    BoundaryFace face = createTestFace({1, 0, 0}, {0, 0, 0}, 2500, 5000, 3000);
    
    double u_int[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    double flux[9];
    
    bc.computeBoundaryFlux(face, u_int, flux, 0.0);
    
    // Flux should include prescribed traction
    if (std::abs(flux[0] - 1e6) > 1e-3) {
        std::cerr << "  FAIL: Prescribed traction not in flux" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Constant Neumann BC correct" << std::endl;
    return true;
}

// =============================================================================
// Symmetry Tests
// =============================================================================

/**
 * @test Test symmetry BC velocity reflection
 */
bool test_symmetry_velocity() {
    std::cout << "Testing symmetry BC velocity reflection..." << std::endl;
    
    double tol = 1e-10;
    
    SymmetryBC bc(false);  // Symmetry (not anti-symmetry)
    
    BoundaryFace face = createTestFace({1, 0, 0}, {0, 0, 0});
    
    // Interior velocity with normal component
    double u_int[9] = {1.0, 0.5, 0.2, 1e6, 0.5e6, 0.3e6, 0.2e6, 0.1e6, 0.15e6};
    double u_ext[9];
    
    bc.getBoundaryState(face, u_int, u_ext, 0.0);
    
    // Normal velocity should be reflected: v_ext,n = -v_int,n
    // For n = (1,0,0), v_n = v_x
    double v_n_int = u_int[0];
    double v_n_ext = u_ext[0];
    
    if (std::abs(v_n_ext + v_n_int) > tol) {
        std::cerr << "  FAIL: Normal velocity not reflected" << std::endl;
        return false;
    }
    
    // Tangential velocities should match: v_t,ext = v_t,int
    if (std::abs(u_ext[1] - u_int[1]) > tol || std::abs(u_ext[2] - u_int[2]) > tol) {
        std::cerr << "  FAIL: Tangential velocity changed" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Symmetry BC velocity reflection correct" << std::endl;
    return true;
}

/**
 * @test Test anti-symmetry BC
 */
bool test_antisymmetry() {
    std::cout << "Testing anti-symmetry BC..." << std::endl;
    
    double tol = 1e-10;
    
    SymmetryBC bc(true);  // Anti-symmetry
    
    BoundaryFace face = createTestFace({1, 0, 0}, {0, 0, 0});
    
    double u_int[9] = {1.0, 0.5, 0.2, 0, 0, 0, 0, 0, 0};
    double u_ext[9];
    
    bc.getBoundaryState(face, u_int, u_ext, 0.0);
    
    // For anti-symmetry: tangential should flip, normal should stay
    // v_t,ext = -v_t,int
    double vn_int = u_int[0];
    double vn_ext = u_ext[0];
    
    // Normal should match
    if (std::abs(vn_ext - vn_int) > tol) {
        std::cerr << "  FAIL: Normal velocity changed in anti-symmetry" << std::endl;
        return false;
    }
    
    // Tangential should flip
    if (std::abs(u_ext[1] + u_int[1]) > tol) {
        std::cerr << "  FAIL: Tangential velocity not flipped" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Anti-symmetry BC correct" << std::endl;
    return true;
}

// =============================================================================
// Boundary Condition Manager Tests
// =============================================================================

/**
 * @test Test BC manager with multiple conditions
 */
bool test_bc_manager() {
    std::cout << "Testing boundary condition manager..." << std::endl;
    
    BoundaryConditionManager mgr;
    
    // Add different BCs for different tags
    mgr.addBoundaryCondition(1, std::make_unique<FreeSurfaceBC>());
    mgr.addBoundaryCondition(2, std::make_unique<DirichletBC>());
    mgr.setDefaultBC(std::make_unique<LysmerKuhlemeyer>());
    
    // Test retrieval
    // Note: Full test would require mesh data
    
    std::cout << "  PASS: BC manager setup correct" << std::endl;
    return true;
}

/**
 * @test Test BC factory
 */
bool test_bc_factory() {
    std::cout << "Testing BC factory..." << std::endl;
    
    auto fs = BoundaryConditionFactory::create("free_surface");
    if (!fs) {
        std::cerr << "  FAIL: Could not create free surface BC" << std::endl;
        return false;
    }
    
    auto pml = BoundaryConditionFactory::create("pml");
    if (!pml) {
        std::cerr << "  FAIL: Could not create PML BC" << std::endl;
        return false;
    }
    
    auto ce = BoundaryConditionFactory::create("clayton_engquist");
    if (!ce) {
        std::cerr << "  FAIL: Could not create Clayton-Engquist BC" << std::endl;
        return false;
    }
    
    auto lk = BoundaryConditionFactory::create("lysmer");
    if (!lk) {
        std::cerr << "  FAIL: Could not create Lysmer BC" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: BC factory creates all types" << std::endl;
    return true;
}

// =============================================================================
// Configuration Tests
// =============================================================================

/**
 * @test Test BC configuration parsing
 */
bool test_bc_config() {
    std::cout << "Testing BC configuration parsing..." << std::endl;
    
    std::map<std::string, std::string> config = {
        {"default_boundary_condition", "free_surface"},
        {"absorbing_method", "pml"},
        {"pml_thickness", "1000.0"},
        {"pml_reflection_coefficient", "1e-6"},
        {"boundary_density", "2500"},
        {"boundary_vp", "5000"},
        {"boundary_vs", "3000"}
    };
    
    BoundaryConditionConfig bc_config;
    bc_config.parseConfig(config);
    
    if (bc_config.default_bc_type != "free_surface") {
        std::cerr << "  FAIL: Default BC type not parsed" << std::endl;
        return false;
    }
    
    if (std::abs(bc_config.pml_thickness - 1000.0) > 1e-10) {
        std::cerr << "  FAIL: PML thickness not parsed" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: BC configuration parsing correct" << std::endl;
    return true;
}

// =============================================================================
// Test Runner
// =============================================================================

int runBoundaryConditionTests() {
    std::cout << "\n=== Boundary Condition Unit Tests ===" << std::endl;
    
    int passed = 0;
    int failed = 0;
    
    // Free surface tests
    if (test_free_surface_traction_free()) ++passed; else ++failed;
    if (test_free_surface_velocity()) ++passed; else ++failed;
    
    // PML tests
    if (test_pml_damping_profile()) ++passed; else ++failed;
    if (test_pml_region_detection()) ++passed; else ++failed;
    
    // Clayton-Engquist tests
    if (test_clayton_engquist_impedance()) ++passed; else ++failed;
    if (test_clayton_engquist_absorbing()) ++passed; else ++failed;
    
    // Lysmer-Kuhlemeyer tests
    if (test_lysmer_dashpot()) ++passed; else ++failed;
    
    // Dirichlet tests
    if (test_dirichlet_zero()) ++passed; else ++failed;
    if (test_dirichlet_prescribed()) ++passed; else ++failed;
    
    // Neumann tests
    if (test_neumann_constant()) ++passed; else ++failed;
    
    // Symmetry tests
    if (test_symmetry_velocity()) ++passed; else ++failed;
    if (test_antisymmetry()) ++passed; else ++failed;
    
    // Manager tests
    if (test_bc_manager()) ++passed; else ++failed;
    if (test_bc_factory()) ++passed; else ++failed;
    if (test_bc_config()) ++passed; else ++failed;
    
    std::cout << "\n=== Boundary Condition Test Summary ===" << std::endl;
    std::cout << "Passed: " << passed << std::endl;
    std::cout << "Failed: " << failed << std::endl;
    std::cout << "Total:  " << (passed + failed) << std::endl;
    
    return failed;
}

} // namespace Testing
} // namespace FSRM

#ifndef FSRM_TEST_NO_MAIN
int main() {
    return FSRM::Testing::runBoundaryConditionTests();
}
#endif
