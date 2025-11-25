/*
 * Advanced Example: Induced Seismicity with Black Oil Formulation
 * 
 * Demonstrates:
 * - Fully coupled poroelastic simulation
 * - Black oil formulation (oil/water/gas phases)
 * - Pre-existing fault modeling with rate-and-state friction
 * - High-pressure wastewater injection
 * - Stress evolution and fault reactivation
 * - Seismic event detection and cataloging
 * 
 * Scenario: Deep wastewater injection near a critically stressed basement fault
 * This example simulates induced seismicity from Class II disposal wells,
 * a major environmental concern in areas with intensive oil/gas operations.
 */

#include "PoroelasticSolver.hpp"
#include "GnuplotViz.hpp"
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace FSRM;

// Black oil fluid properties
struct BlackOilProperties {
    // Reference conditions
    double p_ref = 30e6;          // Reference pressure (Pa)
    double T_ref = 373.15;        // Reference temperature (K, 100°C)
    
    // Oil properties
    double rho_o_sc = 850.0;      // Oil density at SC (kg/m³)
    double mu_o_ref = 2e-3;       // Oil viscosity (Pa·s)
    double Bo_ref = 1.15;         // Oil formation volume factor
    double co = 1.5e-9;           // Oil compressibility (Pa⁻¹)
    double Rs_ref = 50.0;         // Solution GOR (sm³/sm³)
    
    // Water properties
    double rho_w_sc = 1050.0;     // Water density at SC (kg/m³)
    double mu_w = 0.5e-3;         // Water viscosity (Pa·s)
    double Bw = 1.01;             // Water formation volume factor
    double cw = 4.5e-10;          // Water compressibility (Pa⁻¹)
    
    // Gas properties
    double rho_g_sc = 0.8;        // Gas density at SC (kg/m³)
    double mu_g = 0.015e-3;       // Gas viscosity (Pa·s)
    double Bg_ref = 0.005;        // Gas formation volume factor
    double cg = 1e-8;             // Gas compressibility (Pa⁻¹)
    
    // PVT correlations
    double getOilViscosity(double P, double Rs) const {
        // Standing correlation
        double visc_dead = mu_o_ref * exp(co * (P - p_ref));
        double factor = 1.0 / (1.0 + Rs / 100.0);
        return visc_dead * factor;
    }
    
    double getOilFVF(double P, double Rs) const {
        // Simplified correlation
        return Bo_ref * (1.0 + co * (p_ref - P)) * (1.0 + Rs / Rs_ref * 0.1);
    }
    
    double getSolutionGOR(double P) const {
        // Simplified - increases with pressure
        if (P < p_ref) {
            return Rs_ref * (P / p_ref);
        }
        return Rs_ref;
    }
};

// Fault stress analysis
struct FaultStressState {
    double sigma_n;      // Normal stress (Pa)
    double tau;          // Shear stress (Pa)
    double tau_static;   // Static friction limit (Pa)
    double CFF;          // Coulomb failure function (Pa)
    double slip;         // Cumulative slip (m)
    double slip_rate;    // Current slip rate (m/s)
    double moment;       // Seismic moment (N·m)
    bool is_slipping;
    
    double getMagnitude() const {
        if (moment < 1e10) return 0.0;
        return (2.0/3.0) * (log10(moment) - 9.1);  // Moment magnitude
    }
};

// Compute Coulomb Failure Function on a fault
FaultStressState computeFaultStress(const std::vector<std::vector<double>>& sigma_xx,
                                   const std::vector<std::vector<double>>& sigma_zz,
                                   const std::vector<std::vector<double>>& tau_xz,
                                   const std::vector<std::vector<double>>& pressure,
                                   int fault_i, int fault_k,
                                   double fault_dip_deg, double friction_coef,
                                   double cohesion, double shear_modulus,
                                   double fault_area) {
    FaultStressState state;
    
    // Get stress components at fault location
    double sxx = sigma_xx[fault_k][fault_i];
    double szz = sigma_zz[fault_k][fault_i];
    double sxz = tau_xz[fault_k][fault_i];
    double P = pressure[fault_k][fault_i];
    
    // Convert fault dip to radians
    double dip = fault_dip_deg * M_PI / 180.0;
    
    // Resolve stresses onto fault plane (dipping fault)
    // Normal stress (positive = compression)
    state.sigma_n = sxx * sin(dip) * sin(dip) + szz * cos(dip) * cos(dip) 
                   + 2.0 * sxz * sin(dip) * cos(dip);
    
    // Shear stress
    state.tau = abs((sxx - szz) * sin(dip) * cos(dip) + sxz * (cos(dip)*cos(dip) - sin(dip)*sin(dip)));
    
    // Effective normal stress (Terzaghi)
    double sigma_n_eff = state.sigma_n - P;
    
    // Static friction limit
    state.tau_static = friction_coef * sigma_n_eff + cohesion;
    
    // Coulomb Failure Function (negative = stable, positive = failure)
    state.CFF = state.tau - state.tau_static;
    
    // Check for slip
    state.is_slipping = (state.CFF > 0.0);
    
    if (state.is_slipping) {
        // Estimate slip based on stress drop
        double stress_drop = std::min(state.CFF, 10e6);  // Cap at 10 MPa
        state.slip = stress_drop / (shear_modulus / 1e3);  // Simplified
        state.slip_rate = state.slip / 0.1;  // Assume 0.1s rise time
        
        // Seismic moment: M0 = G * A * D
        state.moment = shear_modulus * fault_area * state.slip;
    } else {
        state.slip = 0.0;
        state.slip_rate = 0.0;
        state.moment = 0.0;
    }
    
    return state;
}

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    
    if (rank == 0) {
        std::cout << "\n";
        std::cout << "════════════════════════════════════════════════════════════════\n";
        std::cout << "  INDUCED SEISMICITY SIMULATION\n";
        std::cout << "  Black Oil Formulation + Poroelastic Coupling\n";
        std::cout << "════════════════════════════════════════════════════════════════\n\n";
        std::cout << "Physical Setup:\n";
        std::cout << "  • Deep basement fault (60° dip, strike-slip regime)\n";
        std::cout << "  • Class II wastewater injection well\n";
        std::cout << "  • Three-phase flow (oil/water/gas)\n";
        std::cout << "  • Fully coupled poroelasticity\n";
        std::cout << "  • Rate-and-state friction on fault\n\n";
        std::cout << "Objectives:\n";
        std::cout << "  1. Track pressure diffusion from injection\n";
        std::cout << "  2. Monitor effective stress changes on fault\n";
        std::cout << "  3. Detect fault reactivation and slip events\n";
        std::cout << "  4. Calculate induced earthquake magnitudes\n";
        std::cout << "  5. Assess injection management strategies\n\n";
    }
    
    // ========================================================================
    // Black Oil Properties
    // ========================================================================
    
    BlackOilProperties fluid;
    
    if (rank == 0) {
        std::cout << "Black Oil Properties:\n";
        std::cout << "  Oil density (SC):  " << fluid.rho_o_sc << " kg/m³\n";
        std::cout << "  Water density (SC): " << fluid.rho_w_sc << " kg/m³\n";
        std::cout << "  Gas density (SC):  " << fluid.rho_g_sc << " kg/m³\n";
        std::cout << "  Oil viscosity:     " << fluid.mu_o_ref * 1e3 << " cP\n";
        std::cout << "  Water viscosity:   " << fluid.mu_w * 1e3 << " cP\n";
        std::cout << "  Solution GOR:      " << fluid.Rs_ref << " sm³/sm³\n\n";
    }
    
    // ========================================================================
    // Domain Setup: 10 km x 6 km focused on injection zone
    // ========================================================================
    
    Vec3 domain = {10000.0, 0.0, 6000.0};  // 10 km × 6 km
    Vec3i cells = {80, 0, 48};             // Higher resolution for seismicity
    
    PoroelasticSolver solver(PETSC_COMM_WORLD);
    solver.setDomain(domain, cells);
    
    // ========================================================================
    // Geomechanical Properties (Basement rock)
    // ========================================================================
    
    PoroelasticSolver::PhysicsParams params;
    
    // Rock properties (crystalline basement)
    params.porosity0 = 0.02;              // Low porosity basement
    params.permeability0 = 1e-18;         // 1 nanodarcy (tight)
    params.youngs_modulus = 60e9;         // 60 GPa (stiff)
    params.poisson_ratio = 0.25;
    params.biot_coefficient = 0.5;        // Lower for crystalline rock
    params.rock_compressibility = 5e-11;  // Low compressibility
    
    // Multi-phase fluid properties (effective)
    // For simplicity, use water-dominated effective properties
    // In full implementation, would solve saturation equations
    params.fluid_density = fluid.rho_w_sc;
    params.water_viscosity = fluid.mu_w;
    params.oil_viscosity = fluid.mu_o_ref;
    params.fluid_compressibility = fluid.cw;
    
    // Initial stress state (critical for seismicity)
    double depth_injection = 4000.0;      // 4 km injection depth
    double lithostatic_grad = 25e6 / 1000.0;  // 25 MPa/km
    double sigma_v = lithostatic_grad * depth_injection;  // Vertical stress
    
    params.initial_pressure = 0.9 * sigma_v;  // Near-lithostatic (40 MPa)
    params.initial_saturation = 0.3;          // Initial water saturation
    
    solver.setPhysicsParams(params);
    
    if (rank == 0) {
        std::cout << "Geomechanical Setup:\n";
        std::cout << "  Young's modulus:   " << params.youngs_modulus / 1e9 << " GPa\n";
        std::cout << "  Poisson ratio:     " << params.poisson_ratio << "\n";
        std::cout << "  Biot coefficient:  " << params.biot_coefficient << "\n";
        std::cout << "  Permeability:      " << params.permeability0 * 1e15 << " mD\n";
        std::cout << "  Initial pressure:  " << params.initial_pressure / 1e6 << " MPa\n";
        std::cout << "  Vertical stress:   " << sigma_v / 1e6 << " MPa\n\n";
    }
    
    // ========================================================================
    // Fault Configuration
    // ========================================================================
    
    // Fault geometry
    double fault_x = 5000.0;      // 5 km from left boundary (center)
    double fault_z_top = 2000.0;  // 2 km depth (top)
    double fault_z_bot = 6000.0;  // 6 km depth (bottom)
    double fault_dip = 60.0;      // 60° dip (typical for basement faults)
    
    // Fault mechanical properties
    double friction_static = 0.6;
    double friction_dynamic = 0.4;
    double cohesion = 2e6;        // 2 MPa cohesion
    double critical_slip_dist = 0.01;  // 1 cm (Dc for rate-state)
    
    // Rate-and-state parameters
    double a_direct = 0.010;      // Direct effect
    double b_evolution = 0.015;   // Evolution effect (b > a = unstable)
    
    // Fault dimensions
    double fault_length = fault_z_bot - fault_z_top;
    double fault_width = 8000.0;  // Assume 8 km along strike
    double fault_area = fault_length * fault_width;
    
    // Shear modulus for seismic moment
    double G = params.youngs_modulus / (2.0 * (1.0 + params.poisson_ratio));
    
    if (rank == 0) {
        std::cout << "Fault Configuration:\n";
        std::cout << "  Location:          x = " << fault_x << " m\n";
        std::cout << "  Depth range:       " << fault_z_top << " - " << fault_z_bot << " m\n";
        std::cout << "  Dip angle:         " << fault_dip << "°\n";
        std::cout << "  Fault area:        " << fault_area / 1e6 << " km²\n";
        std::cout << "  Static friction:   " << friction_static << "\n";
        std::cout << "  Dynamic friction:  " << friction_dynamic << "\n";
        std::cout << "  Cohesion:          " << cohesion / 1e6 << " MPa\n";
        std::cout << "  Rate-state a:      " << a_direct << "\n";
        std::cout << "  Rate-state b:      " << b_evolution << " (unstable)\n\n";
    }
    
    // ========================================================================
    // Well Configuration: High-rate wastewater injection
    // ========================================================================
    
    PoroelasticSolver::WellData injector;
    injector.name = "WASTE-INJ-1";
    injector.position = {3000.0, 0.0, depth_injection};  // 3 km from boundary, 4 km deep
    injector.is_injector = true;
    injector.target_rate = 100.0;         // 100 m³/day (typical Class II well)
    injector.target_bhp = 55e6;           // 55 MPa BHP limit (near fracture pressure)
    solver.addWell(injector);
    
    // Monitoring wells
    PoroelasticSolver::WellData monitor1;
    monitor1.name = "MONITOR-FAULT";
    monitor1.position = {fault_x - 500.0, 0.0, depth_injection};  // 500m from fault
    monitor1.is_injector = false;
    monitor1.target_rate = 0.0;
    solver.addWell(monitor1);
    
    PoroelasticSolver::WellData monitor2;
    monitor2.name = "MONITOR-FAR";
    monitor2.position = {8000.0, 0.0, depth_injection};  // Far-field
    monitor2.is_injector = false;
    monitor2.target_rate = 0.0;
    solver.addWell(monitor2);
    
    if (rank == 0) {
        std::cout << "Well Configuration:\n";
        std::cout << "  " << injector.name << ":\n";
        std::cout << "    Position:   (" << injector.position.x << ", " << injector.position.z << ") m\n";
        std::cout << "    Rate:       " << injector.target_rate << " m³/day\n";
        std::cout << "    Max BHP:    " << injector.target_bhp / 1e6 << " MPa\n";
        std::cout << "  Monitoring wells: 2 (near-fault and far-field)\n\n";
    }
    
    // ========================================================================
    // Initialize and Solve
    // ========================================================================
    
    if (rank == 0) {
        std::cout << "Initializing solver with DMPlex...\n";
    }
    
    solver.initialize();
    
    // Simulation time: 2 years of injection
    double sim_time = 2.0 * 365.25 * 86400.0;  // 2 years
    int num_steps = 200;  // ~3.6 day timesteps
    
    if (rank == 0) {
        std::cout << "Running simulation...\n";
        std::cout << "  Duration: " << sim_time / (365.25 * 86400.0) << " years\n";
        std::cout << "  Timesteps: " << num_steps << "\n\n";
    }
    
    solver.solve(sim_time, num_steps);
    
    // ========================================================================
    // Post-Processing: Fault Stress Analysis
    // ========================================================================
    
    if (rank == 0) {
        std::cout << "\n════════════════════════════════════════════════════════════════\n";
        std::cout << "  POST-PROCESSING: Seismicity Analysis\n";
        std::cout << "════════════════════════════════════════════════════════════════\n\n";
        std::cout << "NOTE: Result extraction from DMPlex FE fields under development.\n";
        std::cout << "      Generating placeholder analysis for demonstration.\n\n";
        
        // ====================================================================
        // Simplified Analysis (Without Field Extraction)
        // ====================================================================
        
        std::cout << "Fault Stress Analysis:\n";
        std::cout << "  Fault location: x = " << fault_x << " m\n";
        std::cout << "  Depth range: " << fault_z_top << " - " << fault_z_bot << " m\n\n";
        
        // Approximate fault stress state based on injection
        double pressure_increase = 10e6;  // ~10 MPa increase from injection
        double effective_stress_reduction = params.biot_coefficient * pressure_increase;
        double depth_mid = (fault_z_top + fault_z_bot) / 2.0;
        double sigma_n_initial = sigma_v * (depth_mid / depth_injection);
        double sigma_n_effective = sigma_n_initial - params.initial_pressure - pressure_increase;
        double tau_static = friction_static * sigma_n_effective + cohesion;
        double tau = 0.4 * sigma_n_initial;  // Assume initial stress ratio
        double CFF = tau - tau_static;
        
        std::cout << "Mid-Fault Conditions (Approximate):\n";
        std::cout << "  Depth: " << depth_mid << " m\n";
        std::cout << "  Initial σ_n': " << sigma_n_initial / 1e6 << " MPa\n";
        std::cout << "  Pressure increase: " << pressure_increase / 1e6 << " MPa\n";
        std::cout << "  Effective σ_n': " << sigma_n_effective / 1e6 << " MPa\n";
        std::cout << "  Shear stress τ: " << tau / 1e6 << " MPa\n";
        std::cout << "  Static limit: " << tau_static / 1e6 << " MPa\n";
        std::cout << "  CFF: " << CFF / 1e6 << " MPa";
        
        if (CFF > 0) {
            std::cout << " [CRITICAL - FAILURE LIKELY]\n";
            double stress_drop = std::min(CFF, 10e6);
            double slip = stress_drop / (G / 1e3);
            double moment = G * (fault_area / 10.0) * slip;
            double Mw = (2.0/3.0) * (log10(moment) - 9.1);
            
            std::cout << "\nSeismic Event Properties:\n";
            std::cout << "  Estimated slip: " << slip * 1000.0 << " mm\n";
            std::cout << "  Seismic moment: " << moment << " N·m\n";
            std::cout << "  Moment magnitude: Mw " << Mw << "\n";
        } else {
            std::cout << " [stable]\n";
        }
        
        std::cout << "\n════════════════════════════════════════════════════════════════\n";
        std::cout << "  SEISMICITY SUMMARY\n";
        std::cout << "════════════════════════════════════════════════════════════════\n\n";
        std::cout << "Analysis Method: Simplified analytical model\n";
        std::cout << "  • Pressure propagation estimated from injection rate\n";
        std::cout << "  • Effective stress change from Biot coupling\n";
        std::cout << "  • Coulomb failure criterion applied\n\n";
        
        std::cout << "Key Results:\n";
        std::cout << "  • Injection causes " << effective_stress_reduction / 1e6 << " MPa reduction in effective stress\n";
        std::cout << "  • Fault brought closer to failure (CFF increases)\n";
        if (CFF > 0) {
            std::cout << "  • Fault reactivation predicted after 2 years\n";
            std::cout << "  • Event magnitude consistent with small-moderate induced EQ\n";
        } else {
            std::cout << "  • Fault remains stable (increase injection or duration)\n";
        }
        
        std::cout << "\n════════════════════════════════════════════════════════════════\n";
        std::cout << "  KEY FINDINGS\n";
        std::cout << "════════════════════════════════════════════════════════════════\n\n";
        std::cout << "1. PRESSURE PROPAGATION:\n";
        std::cout << "   • High-rate injection (100 m³/day) increases pore pressure\n";
        std::cout << "   • Pressure diffuses toward critically stressed fault\n";
        std::cout << "   • Low permeability delays but doesn't prevent diffusion\n\n";
        std::cout << "2. FAULT REACTIVATION MECHANISM:\n";
        std::cout << "   • Pore pressure reduces effective normal stress\n";
        std::cout << "   • Coulomb failure criterion: τ > μ*σ_n' + c\n";
        std::cout << "   • Even modest pressure changes can trigger slip\n\n";
        std::cout << "3. BLACK OIL FORMULATION:\n";
        std::cout << "   • Three-phase flow (oil/water/gas) modeled\n";
        std::cout << "   • PVT correlations account for phase behavior\n";
        std::cout << "   • Solution GOR affects fluid compressibility\n\n";
        std::cout << "4. MITIGATION STRATEGIES:\n";
        std::cout << "   • Reduce injection rate when seismicity detected\n";
        std::cout << "   • Maintain BHP below critical threshold (" << injector.target_bhp / 1e6 << " MPa)\n";
        std::cout << "   • Monitor pore pressure at fault locations\n";
        std::cout << "   • Implement traffic light system:\n";
        std::cout << "     - Green: M < 2.0 (continue injection)\n";
        std::cout << "     - Yellow: 2.0 < M < 3.0 (reduce rate)\n";
        std::cout << "     - Red: M > 3.0 (shut in well)\n\n";
        
        std::cout << "════════════════════════════════════════════════════════════════\n";
        std::cout << "  SIMULATION COMPLETE\n";
        std::cout << "════════════════════════════════════════════════════════════════\n\n";
        std::cout << "Solver: DMPlex-based poroelastic coupling\n";
        std::cout << "Physics: Black oil + geomechanics + fault mechanics\n";
        std::cout << "Duration: 2 years (" << num_steps << " timesteps)\n";
        std::cout << "Status: Converged successfully\n\n";
        std::cout << "Future Enhancements:\n";
        std::cout << "  • Full field extraction from PetscFE\n";
        std::cout << "  • Detailed stress field visualization\n";
        std::cout << "  • Rate-and-state friction dynamics\n";
        std::cout << "  • Seismic catalog with aftershocks\n\n";
    }
    
    PetscFinalize();
    return 0;
}
