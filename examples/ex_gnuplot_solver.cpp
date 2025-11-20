#include "PoroelasticSolver.hpp"
#include "GnuplotViz.hpp"
#include <mpi.h>
#include <iostream>

/**
 * @brief Simple example using PoroelasticSolver library
 * 
 * This demonstrates how to run a poroelastic simulation in ~20 lines
 * instead of manually implementing 1000+ lines of PETSc setup, physics,
 * analytical Jacobian computation, and visualization.
 * 
 * The library handles:
 * - DMPlex mesh creation (unstructured grid support)
 * - PetscFE field setup (5 DOF: P, Sw, ux, uz, phi)
 * - PETSc TS implicit time integration
 * - Analytical Jacobian computation (hand-coded, no AD)
 * - Result extraction and visualization
 */

using namespace ResSim;

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    
    // ========================================================================
    // Setup (10 lines vs 200+ manual)
    // ========================================================================
    
    // Domain: 1000m x 500m with 50x25 cells
    Vec3 domain = {1000.0, 0.0, 500.0};
    Vec3i cells = {50, 0, 25};
    
    PoroelasticSolver solver(PETSC_COMM_WORLD);
    solver.setDomain(domain, cells);
    
    // Physics parameters
    PoroelasticSolver::PhysicsParams params;
    params.porosity0 = 0.25;
    params.permeability0 = 200e-15;  // 200 mD
    params.youngs_modulus = 15e9;    // 15 GPa
    params.poisson_ratio = 0.2;
    params.biot_coefficient = 0.8;
    params.initial_pressure = 25e6;  // 25 MPa
    params.initial_saturation = 0.3;
    solver.setPhysicsParams(params);
    
    // Add injection well
    PoroelasticSolver::WellData well;
    well.name = "INJ-1";
    well.position = {500.0, 0.0, 250.0};
    well.is_injector = true;
    well.target_rate = 100.0;  // 100 m³/day
    well.target_bhp = 40e6;    // 40 MPa
    solver.addWell(well);
    
    solver.initialize();
    
    // ========================================================================
    // Solve (1 line vs 300+ manual PETSc TS setup + analytical Jacobian)
    // ========================================================================
    
    double sim_time = 365.0 * 86400.0;  // 1 year in seconds
    int num_steps = 100;
    
    std::cout << "Running poroelastic simulation..." << std::endl;
    std::cout << "  Domain: " << domain.x << "m x " << domain.z << "m" << std::endl;
    std::cout << "  Grid: " << cells.x << " x " << cells.z << std::endl;
    std::cout << "  Time: " << sim_time / 86400.0 << " days" << std::endl;
    std::cout << "  Steps: " << num_steps << std::endl;
    
    solver.solve(sim_time, num_steps);
    
    // ========================================================================
    // Visualization (5 lines vs 500+ manual plotting)
    // ========================================================================
    
    std::vector<std::vector<double>> P, Sw, ux, uz, phi;
    solver.getPressure(P);
    solver.getSaturation(Sw);
    solver.getDisplacement(ux, uz);
    solver.getPorosity(phi);
    
    GnuplotViz viz("output_poroelastic");
    
    // Individual field plots
    viz.plot2DField(P, domain.x, domain.z, "Pressure Distribution", 
                   "X (m)", "Z (m)", "viridis", "pressure", 0, 0);
    viz.plot2DField(Sw, domain.x, domain.z, "Water Saturation", 
                   "X (m)", "Z (m)", "plasma", "saturation", 0, 1);
    viz.plot2DField(ux, domain.x, domain.z, "X-Displacement", 
                   "X (m)", "Z (m)", "coolwarm", "disp_x", 0, 0);
    viz.plot2DField(uz, domain.x, domain.z, "Z-Displacement", 
                   "X (m)", "Z (m)", "coolwarm", "disp_z", 0, 0);
    viz.plot2DField(phi, domain.x, domain.z, "Porosity", 
                   "X (m)", "Z (m)", "jet", "porosity", 0, 0);
    
    // Multi-panel summary
    std::vector<std::vector<std::vector<double>>> fields = {P, Sw, ux, uz};
    std::vector<std::string> titles = {"Pressure (Pa)", "Saturation", "Ux (m)", "Uz (m)"};
    std::vector<std::pair<double,double>> ranges = {{0,0}, {0,1}, {0,0}, {0,0}};
    viz.plotMultiPanel(fields, titles, domain.x, domain.z, "viridis", "multi_panel", ranges);
    
    // Pressure with well locations
    std::vector<PoroelasticSolver::WellData> wells = solver.getWells();
    std::vector<GnuplotViz::WellMarker> well_markers;
    for (const auto& w : wells) {
        well_markers.push_back({w.position.x, w.position.z, w.name, w.is_injector});
    }
    viz.addWells(well_markers);
    viz.plot2DField(P, domain.x, domain.z, "Pressure with Wells", 
                   "X (m)", "Z (m)", "viridis", "pressure_wells", 0, 0);
    
    // Well performance summary
    std::cout << "\n=== Well Performance ===" << std::endl;
    for (const auto& w : wells) {
        std::cout << w.name << ":" << std::endl;
        std::cout << "  BHP: " << w.current_bhp / 1e6 << " MPa" << std::endl;
        std::cout << "  Oil rate: " << w.oil_rate << " m³/day" << std::endl;
        std::cout << "  Water rate: " << w.water_rate << " m³/day" << std::endl;
        std::cout << "  Cum. oil: " << w.cumulative_oil << " m³" << std::endl;
        std::cout << "  Cum. water: " << w.cumulative_water << " m³" << std::endl;
    }
    
    std::cout << "\nResults saved to PNG files:" << std::endl;
    std::cout << "  output_poroelastic/pressure.png" << std::endl;
    std::cout << "  output_poroelastic/saturation.png" << std::endl;
    std::cout << "  output_poroelastic/disp_x.png, disp_z.png" << std::endl;
    std::cout << "  output_poroelastic/porosity.png" << std::endl;
    std::cout << "  output_poroelastic/multi_panel.png (2x2 summary)" << std::endl;
    std::cout << "  output_poroelastic/pressure_wells.png" << std::endl;
    
    PetscFinalize();
    return 0;
}
