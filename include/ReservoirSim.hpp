#ifndef RESERVOIR_SIM_HPP
#define RESERVOIR_SIM_HPP

#include <petsc.h>
#include <petscts.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <petscdmplex.h>
#include <petscfe.h>

#include <string>
#include <vector>
#include <memory>
#include <map>
#include <functional>

namespace ResSim {

// Forward declarations
class Simulator;
class PhysicsKernel;
class FluidModel;
class SolidModel;
class FractureModel;
class WellModel;
class EclipseIO;

// Enumerations
enum class FluidModelType {
    SINGLE_COMPONENT,
    BLACK_OIL,
    COMPOSITIONAL,
    THERMAL_COMPOSITIONAL
};

enum class SolidModelType {
    ELASTIC,
    VISCOELASTIC,
    POROELASTIC,
    THERMOELASTIC
};

enum class FractureType {
    NATURAL,
    INDUCED_HYDRAULIC,
    INDUCED_THERMAL
};

enum class WellType {
    PRODUCER,
    INJECTOR,
    OBSERVATION
};

enum class PhysicsType {
    FLUID_FLOW,
    GEOMECHANICS,
    THERMAL,
    PARTICLE_TRANSPORT,
    FRACTURE_PROPAGATION,
    TIDAL_FORCES,
    CHEMICAL_REACTION
};

// Configuration structures
struct SimulationConfig {
    std::string input_file;
    std::string output_file;
    std::string output_format = "ECLIPSE"; // ECLIPSE, VTK, HDF5
    
    double start_time = 0.0;
    double end_time = 1.0;
    double dt_initial = 0.01;
    double dt_min = 1e-10;
    double dt_max = 1.0;
    
    int max_timesteps = 10000;
    int output_frequency = 10;
    
    bool enable_adaptive_timestepping = true;
    bool enable_checkpointing = true;
    int checkpoint_frequency = 100;
    
    // Solver options
    double rtol = 1e-6;
    double atol = 1e-8;
    int max_nonlinear_iterations = 50;
    int max_linear_iterations = 1000;
    
    // Physics flags
    bool enable_geomechanics = false;
    bool enable_thermal = false;
    bool enable_fractures = false;
    bool enable_particle_transport = false;
    bool enable_faults = false;
    bool enable_tidal_forces = false;
    
    FluidModelType fluid_model = FluidModelType::SINGLE_COMPONENT;
    SolidModelType solid_model = SolidModelType::ELASTIC;
};

struct GridConfig {
    int nx = 10, ny = 10, nz = 10;
    double Lx = 1000.0, Ly = 1000.0, Lz = 100.0; // meters
    bool use_unstructured = false;
    std::string mesh_file;
};

struct MaterialProperties {
    // Rock properties
    double porosity = 0.2;
    double permeability_x = 100.0; // mD
    double permeability_y = 100.0;
    double permeability_z = 10.0;
    double compressibility = 1e-9; // 1/Pa
    
    // Mechanical properties
    double youngs_modulus = 10e9; // Pa
    double poisson_ratio = 0.25;
    double density = 2500.0; // kg/m^3
    double biot_coefficient = 1.0;
    
    // Viscoelastic properties (if enabled)
    double relaxation_time = 1e6; // seconds
    double viscosity = 1e19; // Pa·s
    
    // Thermal properties
    double thermal_conductivity = 2.5; // W/(m·K)
    double heat_capacity = 900.0; // J/(kg·K)
    double thermal_expansion = 1e-5; // 1/K
    
    // Fracture properties
    double fracture_toughness = 1e6; // Pa·m^0.5
    double fracture_energy = 100.0; // J/m^2
};

struct FluidProperties {
    // Single phase/component
    double density = 1000.0; // kg/m^3
    double viscosity = 0.001; // Pa·s
    double compressibility = 1e-9; // 1/Pa
    
    // Black oil properties
    double oil_density_std = 800.0;
    double gas_density_std = 1.0;
    double water_density_std = 1000.0;
    double oil_viscosity = 0.005;
    double gas_viscosity = 0.00001;
    double water_viscosity = 0.001;
    double solution_GOR = 100.0; // scf/stb
    
    // Compositional properties
    std::vector<double> component_mw; // Molecular weights
    std::vector<double> component_Tc; // Critical temperatures
    std::vector<double> component_Pc; // Critical pressures
    std::vector<double> component_omega; // Acentric factors
};

} // namespace ResSim

#endif // RESERVOIR_SIM_HPP
