#ifndef FSRM_HPP
#define FSRM_HPP

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

namespace FSRM {

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
    THERMOELASTIC,
    ELASTODYNAMIC,           // Full dynamic elasticity with inertia
    POROELASTODYNAMIC        // Biot's dynamic poroelasticity
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
    CHEMICAL_REACTION,
    ELASTODYNAMICS,          // Elastic wave propagation
    POROELASTODYNAMICS       // Coupled fluid-solid wave propagation
};

// GPU execution mode
enum class GPUExecutionMode {
    CPU_ONLY,           // Run entirely on CPU
    GPU_ONLY,           // Run entirely on GPU
    HYBRID,             // Automatic load balancing between CPU and GPU
    CPU_FALLBACK        // Try GPU first, fallback to CPU if not available
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
    
    // GPU configuration
    bool use_gpu = false;
    GPUExecutionMode gpu_mode = GPUExecutionMode::CPU_FALLBACK;
    int gpu_device_id = 0;
    bool gpu_verbose = false;
    double gpu_memory_fraction = 0.8;  // Fraction of GPU memory to use
    
    // Solver options
    double rtol = 1e-6;
    double atol = 1e-8;
    int max_nonlinear_iterations = 50;
    int max_linear_iterations = 1000;
    
    // GPU solver options (when use_gpu = true)
    bool use_gpu_preconditioner = true;
    bool use_gpu_matrix_assembly = true;
    bool pin_host_memory = true;  // Pin CPU memory for faster GPU transfers
    
    // Physics flags
    bool enable_geomechanics = false;
    bool enable_thermal = false;
    bool enable_fractures = false;
    bool enable_particle_transport = false;
    bool enable_faults = false;
    bool enable_tidal_forces = false;
    bool enable_elastodynamics = false;           // Enable wave propagation
    bool enable_poroelastodynamics = false;       // Enable coupled fluid-solid waves
    
    // Dynamic simulation mode
    bool use_dynamic_mode = false;                 // Full dynamic vs quasi-static
    bool use_static_triggering = false;            // Static stress triggers dynamic event
    double dynamic_trigger_threshold = 1.0e6;      // Stress threshold for triggering (Pa)
    double dynamic_event_duration = 10.0;          // Duration of dynamic event (s)
    
    // Permeability change from dynamic waves
    bool enable_dynamic_permeability_change = false;
    double permeability_sensitivity = 1.0;         // Sensitivity to strain/stress
    double permeability_recovery_time = 100.0;     // Recovery time constant (s)
    
    FluidModelType fluid_model = FluidModelType::SINGLE_COMPONENT;
    SolidModelType solid_model = SolidModelType::ELASTIC;
};

/**
 * @brief Mesh type enumeration
 */
enum class MeshType {
    CARTESIAN,         ///< Structured Cartesian grid
    CORNER_POINT,      ///< Corner-point grid (Eclipse style)
    GMSH,              ///< Gmsh unstructured mesh
    EXODUS,            ///< Exodus II mesh format
    CUSTOM             ///< Custom mesh loader
};

struct GridConfig {
    // Grid dimensions (for structured grids)
    int nx = 10, ny = 10, nz = 10;
    double Lx = 1000.0, Ly = 1000.0, Lz = 100.0; // meters
    
    // Domain origin (can be geo-referenced)
    double origin_x = 0.0;
    double origin_y = 0.0;
    double origin_z = 0.0;
    
    // Mesh type
    MeshType mesh_type = MeshType::CARTESIAN;
    bool use_unstructured = false;
    std::string mesh_file;
    
    // Gmsh-specific options
    std::string gmsh_physical_volume;       ///< Physical group for main volume
    std::vector<std::string> gmsh_boundaries;  ///< Physical groups for boundaries
    int gmsh_refinement_level = 0;          ///< Additional mesh refinement
    
    // Coordinate Reference System (CRS)
    std::string input_crs;                  ///< Input coordinates CRS (e.g., "EPSG:4326")
    std::string model_crs;                  ///< Model/simulation CRS (e.g., "EPSG:32610")
    bool use_local_coordinates = true;      ///< Shift to local origin
    double local_origin_x = 0.0;            ///< Local coordinate origin X
    double local_origin_y = 0.0;            ///< Local coordinate origin Y
    double local_origin_z = 0.0;            ///< Local coordinate origin Z
    bool auto_detect_utm = false;           ///< Auto-detect UTM zone from coordinates
    
    // Grid quality settings
    double min_cell_volume = 1e-10;         ///< Minimum cell volume (m³)
    double max_aspect_ratio = 100.0;        ///< Maximum cell aspect ratio
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
    
    // Wave propagation properties
    double p_wave_velocity = 5000.0;     // m/s (typical for rock)
    double s_wave_velocity = 3000.0;     // m/s
    double damping_alpha = 0.01;         // Rayleigh damping (mass proportional)
    double damping_beta = 0.001;         // Rayleigh damping (stiffness proportional)
    double quality_factor = 100.0;       // Seismic quality factor Q
    
    // Permeability dynamics
    double permeability_strain_coeff = 1e-7;     // dK/dε (1/strain)
    double permeability_stress_coeff = 1e-15;    // dK/dσ (1/Pa)
    double min_permeability = 0.001;             // mD (minimum allowable)
    double max_permeability = 10000.0;           // mD (maximum allowable)
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

} // namespace FSRM

#endif // FSRM_HPP
