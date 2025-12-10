#ifndef FSRM_HPP
#define FSRM_HPP

#include <petsc.h>
#include <petscts.h>
#include <petscdm.h>
#include <petscdmplex.h>
#include <petscfe.h>

#include <string>
#include <vector>
#include <memory>
#include <map>
#include <functional>
#include <tuple>
#include <utility>

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

/**
 * @brief Physics types supported by FSRM kernels
 * 
 * Each physics type represents a distinct set of governing equations
 * that can be solved independently or coupled with other physics.
 * All calculations are performed in SI base units (m, kg, s).
 */
enum class PhysicsType {
    // Core reservoir physics
    FLUID_FLOW,              ///< Darcy flow (single/multi-phase)
    GEOMECHANICS,            ///< Static stress/strain (elastic, poroelastic)
    THERMAL,                 ///< Heat conduction and convection
    PARTICLE_TRANSPORT,      ///< Proppant, tracers, solids transport
    FRACTURE_PROPAGATION,    ///< Cohesive zone / LEFM fracture growth
    TIDAL_FORCES,            ///< Earth tides and loading effects
    CHEMICAL_REACTION,       ///< Reactive transport / geochemistry
    
    // Dynamic wave propagation
    ELASTODYNAMICS,          ///< Elastic wave propagation with inertia
    POROELASTODYNAMICS,      ///< Coupled fluid-solid wave (Biot waves)
    
    // Explosion and impact physics
    EXPLOSION_SOURCE,        ///< Nuclear/chemical explosion source models
    NEAR_FIELD_DAMAGE,       ///< Cavity, crushing, fracturing zones
    HYDRODYNAMIC,            ///< High-pressure hydrodynamic flow (shock)
    CRATER_FORMATION,        ///< Impact crater excavation and modification
    
    // Atmospheric physics
    ATMOSPHERIC_BLAST,       ///< Atmospheric blast wave propagation
    ATMOSPHERIC_ACOUSTIC,    ///< Linear acoustic wave in atmosphere
    INFRASOUND,              ///< Low-frequency atmospheric sound (< 20 Hz)
    THERMAL_RADIATION,       ///< Fireball thermal pulse
    EMP,                     ///< Electromagnetic pulse
    FALLOUT,                 ///< Radioactive debris transport
    
    // Tsunami and surface water
    TSUNAMI,                 ///< Shallow water wave propagation
    SURFACE_DEFORMATION,     ///< Ground surface uplift/subsidence
    
    // Hydrodynamics and Ocean Physics
    OCEAN_CIRCULATION,       ///< Large-scale ocean current modeling
    COASTAL_HYDRODYNAMICS,   ///< Nearshore and coastal flows
    ESTUARINE_DYNAMICS,      ///< Estuary and river-ocean interaction
    STORM_SURGE,             ///< Storm-driven coastal flooding
    OCEAN_ACOUSTICS,         ///< Underwater sound propagation
    SEDIMENT_TRANSPORT,      ///< Coastal/marine sediment dynamics
    WAVE_DYNAMICS,           ///< Surface gravity wave modeling
    INTERNAL_WAVES,          ///< Internal/baroclinic wave propagation
    OCEAN_MIXING,            ///< Turbulent mixing and diffusion
    THERMOHALINE,            ///< Temperature/salinity dynamics
    
    // Volcanic physics
    MAGMA_CHAMBER,           ///< Magma chamber thermodynamics
    CONDUIT_FLOW,            ///< Magma conduit ascent/fragmentation
    ERUPTION_COLUMN,         ///< Buoyant volcanic plume
    PYROCLASTIC_FLOW,        ///< Pyroclastic density currents
    LAVA_FLOW,               ///< Lava flow propagation
    LAHAR,                   ///< Volcanic debris flows
    VOLCANIC_DEFORMATION,    ///< Volcanic ground deformation (Mogi, etc.)
    VOLCANIC_SEISMICITY,     ///< LP, VT, VLP, tremor sources
    VOLCANIC_GAS,            ///< SO2, CO2 emission and dispersal
    TEPHRA_DISPERSAL         ///< Ash fall modeling
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
    
    // =========================================================================
    // Physics Flags - Core Reservoir
    // =========================================================================
    bool enable_geomechanics = false;
    bool enable_thermal = false;
    bool enable_fractures = false;
    bool enable_particle_transport = false;
    bool enable_faults = false;
    bool enable_tidal_forces = false;
    bool enable_elastodynamics = false;           // Elastic wave propagation
    bool enable_poroelastodynamics = false;       // Coupled fluid-solid waves
    
    // =========================================================================
    // Physics Flags - Explosion and Impact
    // =========================================================================
    bool enable_explosion_source = false;          // Explosion source models
    bool enable_near_field_damage = false;         // Cavity/damage zone evolution
    bool enable_hydrodynamic = false;              // High-pressure shock flow
    bool enable_crater_formation = false;          // Impact cratering
    
    // =========================================================================
    // Physics Flags - Atmospheric
    // =========================================================================
    bool enable_atmospheric_blast = false;         // Atmospheric blast wave
    bool enable_atmospheric_acoustic = false;      // Linear acoustic propagation
    bool enable_infrasound = false;                // Low-frequency sound (< 20 Hz)
    bool enable_thermal_radiation = false;         // Fireball thermal pulse
    bool enable_emp = false;                       // Electromagnetic pulse
    bool enable_fallout = false;                   // Radioactive debris transport
    
    // =========================================================================
    // Physics Flags - Surface and Water
    // =========================================================================
    bool enable_tsunami = false;                   // Shallow water waves
    bool enable_surface_deformation = false;       // Ground surface changes
    
    // =========================================================================
    // Physics Flags - Hydrodynamics and Ocean Physics
    // =========================================================================
    bool enable_ocean_circulation = false;         // Large-scale ocean currents
    bool enable_coastal_hydrodynamics = false;     // Nearshore/coastal flows
    bool enable_estuarine_dynamics = false;        // Estuary/river-ocean
    bool enable_storm_surge = false;               // Storm-driven flooding
    bool enable_ocean_acoustics = false;           // Underwater acoustics
    bool enable_sediment_transport = false;        // Sediment dynamics
    bool enable_wave_dynamics = false;             // Surface gravity waves
    bool enable_internal_waves = false;            // Internal waves
    bool enable_ocean_mixing = false;              // Turbulent mixing
    bool enable_thermohaline = false;              // Temperature/salinity
    
    // =========================================================================
    // Physics Flags - Volcanic Systems
    // =========================================================================
    bool enable_magma_chamber = false;             // Magma chamber dynamics
    bool enable_conduit_flow = false;              // Conduit ascent/fragmentation
    bool enable_eruption_column = false;           // Buoyant eruption column
    bool enable_pyroclastic_flow = false;          // PDC modeling
    bool enable_lava_flow = false;                 // Lava flow propagation
    bool enable_lahar = false;                     // Volcanic debris flows
    bool enable_volcanic_deformation = false;      // Volcanic deformation (Mogi)
    bool enable_volcanic_seismicity = false;       // LP, VT, tremor sources
    bool enable_volcanic_gas = false;              // SO2, CO2 dispersal
    bool enable_tephra_dispersal = false;          // Ash fall modeling
    
    // =========================================================================
    // Dynamic Simulation Mode
    // =========================================================================
    bool use_dynamic_mode = false;                 // Full dynamic vs quasi-static
    bool use_static_triggering = false;            // Static stress triggers dynamic
    double dynamic_trigger_threshold = 1.0e6;      // Stress threshold (Pa)
    double dynamic_event_duration = 10.0;          // Dynamic event duration (s)
    
    // =========================================================================
    // Permeability Dynamics
    // =========================================================================
    bool enable_dynamic_permeability_change = false;
    double permeability_sensitivity = 1.0;         // Sensitivity to strain/stress
    double permeability_recovery_time = 100.0;     // Recovery time constant (s)
    
    // =========================================================================
    // IMEX Time Integration
    // =========================================================================
    bool enable_imex_transition = false;           // Enable IMEX mode switching
    std::string imex_trigger_type = "COULOMB_FAILURE";  // Trigger condition
    std::string imex_settling_type = "COMBINED";   // Settling condition
    double imex_min_dynamic_duration = 1.0;        // Min explicit duration (s)
    double imex_max_dynamic_duration = 60.0;       // Max explicit duration (s)
    
    // =========================================================================
    // Infrasound Modeling Configuration
    // =========================================================================
    bool infrasound_include_topography = true;     // Include terrain effects
    bool infrasound_include_wind = true;           // Include wind effects
    bool infrasound_include_attenuation = true;    // Include atmospheric absorption
    double infrasound_frequency_min = 0.01;        // Minimum frequency (Hz)
    double infrasound_frequency_max = 20.0;        // Maximum frequency (Hz)
    int infrasound_num_frequencies = 100;          // Number of frequency bins
    std::string infrasound_source_type = "POINT";  // POINT, EXPLOSION, VOLCANIC
    
    // =========================================================================
    // Unit System Configuration
    // =========================================================================
    // All internal calculations use SI base units (m, kg, s).
    // Input/output units can be configured per quantity.
    std::string input_unit_system = "SI";          // SI, FIELD, METRIC
    std::string output_unit_system = "SI";         // SI, FIELD, METRIC
    
    // Per-quantity output units (override system default)
    std::string output_pressure_unit = "";         // e.g., "psi", "MPa"
    std::string output_length_unit = "";           // e.g., "ft", "m"
    std::string output_time_unit = "";             // e.g., "day", "s"
    std::string output_permeability_unit = "";     // e.g., "mD", "m2"
    std::string output_temperature_unit = "";      // e.g., "degC", "K"
    
    // =========================================================================
    // Model Types
    // =========================================================================
    FluidModelType fluid_model = FluidModelType::SINGLE_COMPONENT;
    SolidModelType solid_model = SolidModelType::ELASTIC;
    
    // =========================================================================
    // HIGH-FIDELITY OPTIONS (Optional Advanced Features)
    // =========================================================================
    // These are optional extensions that provide higher physical fidelity
    // at the cost of additional computational expense. Enable only when
    // the additional accuracy is required.
    
    // --- High-Fidelity Fluid Flow ---
    bool enable_high_fidelity_flow = false;      ///< Master switch for advanced flow
    bool enable_non_darcy_flow = false;          ///< Forchheimer non-Darcy correction
    bool enable_klinkenberg = false;             ///< Gas slippage in tight formations
    bool enable_dual_porosity = false;           ///< Matrix-fracture continuum
    bool enable_dual_permeability = false;       ///< Matrix flow in dual-porosity
    bool enable_triple_porosity = false;         ///< Matrix-fracture-vug system
    bool enable_non_isothermal_flow = false;     ///< Temperature-dependent flow
    bool enable_dynamic_relperm = false;         ///< Rate-dependent rel perm (CDC)
    bool enable_miscible_flow = false;           ///< Miscible/near-miscible EOR
    
    // Non-Darcy parameters
    double forchheimer_beta = 1e5;               ///< Non-Darcy coefficient [1/m]
    std::string beta_correlation = "geertsma";   ///< constant, cooke, geertsma, etc.
    
    // Dual-porosity parameters
    double matrix_porosity = 0.2;                ///< Matrix continuum porosity
    double matrix_permeability = 1e-18;          ///< Matrix permeability [m²]
    double matrix_block_size = 1.0;              ///< Matrix block dimension [m]
    double fracture_porosity = 0.01;             ///< Fracture continuum porosity
    double fracture_permeability = 1e-12;        ///< Fracture permeability [m²]
    std::string shape_factor_model = "kazemi";   ///< kazemi, coats, gilman, lemonnier
    
    // Miscible flow parameters
    double minimum_miscibility_pressure = 15e6;  ///< MMP [Pa]
    double todd_longstaff_omega = 0.67;          ///< Mixing parameter
    double longitudinal_dispersivity = 0.1;      ///< α_L [m]
    double transverse_dispersivity = 0.01;       ///< α_T [m]
    
    // --- High-Fidelity Geomechanics ---
    bool enable_high_fidelity_geomech = false;   ///< Master switch for advanced geomech
    bool enable_finite_strain = false;           ///< Large deformation analysis
    bool enable_viscoplasticity = false;         ///< Rate-dependent plasticity
    bool enable_creep = false;                   ///< Time-dependent deformation
    bool enable_hypoplasticity = false;          ///< Granular material model
    bool enable_gradient_damage = false;         ///< Non-local damage regularization
    
    // Finite strain model
    std::string hyperelastic_model = "neo_hookean"; ///< neo_hookean, mooney_rivlin, ogden
    std::string plasticity_formulation = "multiplicative"; ///< multiplicative, hypoelastic
    
    // Viscoplasticity parameters
    std::string viscoplastic_model = "perzyna";  ///< perzyna, duvaut_lions
    double viscoplastic_fluidity = 1e-6;         ///< Fluidity parameter [1/(Pa·s)]
    double viscoplastic_exponent = 3.0;          ///< Rate exponent
    
    // Creep parameters
    std::string creep_model = "power_law";       ///< power_law, norton, lemaitre
    double creep_coefficient_A = 1e-30;          ///< Creep coefficient
    double creep_stress_exponent = 4.0;          ///< Stress exponent n
    double creep_activation_energy = 54e3;       ///< Q [J/mol]
    
    // --- Advanced Coupled Physics ---
    bool enable_advanced_coupling = false;       ///< Master switch for advanced coupling
    bool enable_full_biot_dynamics = false;      ///< Full Biot with slow P-wave
    bool enable_thm_coupling = false;            ///< Thermo-Hydro-Mechanical
    bool enable_thmc_coupling = false;           ///< Add Chemical to THM
    bool enable_unsaturated_flow = false;        ///< Partially saturated mechanics
    bool enable_multiscale_coupling = false;     ///< FE² / homogenization
    
    // Full Biot parameters
    std::string biot_formulation = "u_p";        ///< u_p, u_w, u_p_w
    bool biot_high_frequency = true;             ///< Include inertial coupling
    double biot_tortuosity = 2.0;                ///< Porous media tortuosity
    
    // THM coupling
    std::string thm_coupling_type = "iterative"; ///< sequential, iterative, monolithic
    double thermal_expansion_solid = 1e-5;       ///< α_T [1/K]
    double thermal_expansion_fluid = 2.1e-4;     ///< β_T [1/K]
    double thermal_pressurization_coeff = 0.1e6; ///< Λ [Pa/K]
    
    // Unsaturated flow
    std::string swrc_model = "van_genuchten";    ///< van_genuchten, brooks_corey
    double vg_alpha = 0.01;                      ///< Van Genuchten α [1/Pa]
    double vg_n = 1.5;                           ///< Van Genuchten n
    double residual_saturation = 0.1;            ///< S_r
    
    // --- High-Fidelity Numerics ---
    bool enable_high_fidelity_numerics = false;  ///< Master switch for advanced numerics
    bool enable_dg_stabilization = false;        ///< Discontinuous Galerkin
    bool enable_supg_stabilization = false;      ///< SUPG/GLS for advection
    bool enable_advanced_time_integration = false; ///< Higher-order time schemes
    bool enable_amr = false;                     ///< Adaptive mesh refinement
    bool enable_p_adaptivity = false;            ///< Polynomial adaptivity
    bool enable_physics_preconditioner = false;  ///< Physics-based preconditioning
    
    // DG parameters
    std::string dg_flux_type = "upwind";         ///< upwind, lax_friedrichs, roe
    std::string dg_penalty_type = "sipg";        ///< sipg, nipg, bassi_rebay
    double dg_penalty_coefficient = 10.0;        ///< Penalty σ
    
    // Time integration
    std::string time_scheme = "backward_euler";  ///< backward_euler, newmark, generalized_alpha
    double newmark_beta = 0.25;                  ///< Newmark β
    double newmark_gamma = 0.5;                  ///< Newmark γ
    double generalized_alpha_rho = 0.5;          ///< ρ_∞ for gen-α
    
    // AMR parameters
    std::string refinement_indicator = "gradient"; ///< gradient, residual, recovery
    double refine_fraction = 0.3;                ///< Top fraction to refine
    double coarsen_fraction = 0.1;               ///< Bottom fraction to coarsen
    int max_refinement_level = 5;                ///< Maximum refinement depth
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
    
    // Gmsh material domain mappings (physical_group_name -> material_section)
    std::vector<std::pair<std::string, std::string>> gmsh_material_domains;
    
    // Gmsh fault surface mappings (physical_group_name -> fault_section, use_split_nodes)
    std::vector<std::tuple<std::string, std::string, bool>> gmsh_fault_surfaces;
    
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
    
    // Residual saturations for relative permeability calculations
    double residual_saturation = 0.2;         // Oil residual saturation (Sor)
    double water_residual_saturation = 0.2;   // Water residual/connate saturation (Swc/Swr)
    double gas_residual_saturation = 0.05;    // Critical gas saturation (Sgc)
    double water_saturation_max = 0.8;        // Maximum water saturation (1 - Sor)
    
    // Corey relative permeability model parameters
    double corey_exponent_oil = 2.0;          // Corey exponent for oil (no)
    double corey_exponent_water = 2.0;        // Corey exponent for water (nw)
    double corey_exponent_gas = 2.0;          // Corey exponent for gas (ng)
    double kr_max_oil = 1.0;                  // Maximum oil relative permeability
    double kr_max_water = 0.5;                // Maximum water relative permeability (endpoint)
    double kr_max_gas = 0.8;                  // Maximum gas relative permeability (endpoint)
    
    // Phase compressibilities
    double oil_compressibility = 1.5e-9;      // Oil compressibility [1/Pa]
    double water_compressibility = 4.5e-10;   // Water compressibility [1/Pa]
    double gas_compressibility = 1e-8;        // Gas compressibility [1/Pa]
    
    // Temperature for PVT calculations
    double reservoir_temperature = 350.0;     // Reservoir temperature [K]
    
    // Gas properties
    double gas_molecular_weight = 0.016;      // Molecular weight [kg/mol] (methane default)
    
    // Transport properties
    double longitudinal_dispersivity = 0.1;   // Longitudinal dispersivity [m]
    double transverse_dispersivity = 0.01;    // Transverse dispersivity [m]
    
    // Compositional properties
    std::vector<double> component_mw; // Molecular weights
    std::vector<double> component_Tc; // Critical temperatures
    std::vector<double> component_Pc; // Critical pressures
    std::vector<double> component_omega; // Acentric factors
};

} // namespace FSRM

#endif // FSRM_HPP
