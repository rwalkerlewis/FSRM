#include "Simulator.hpp"
#include "Visualization.hpp"
#include "ConfigReader.hpp"
#include "ImplicitExplicitTransition.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

namespace FSRM {

Simulator::Simulator(MPI_Comm comm_in) 
    : comm(comm_in), dm(nullptr), ts(nullptr), 
      solution(nullptr), solution_old(nullptr),
      jacobian(nullptr), snes(nullptr), ksp(nullptr), pc(nullptr),
      prob(nullptr), current_time(0.0), dt(0.01), timestep(0), output_counter(0) {
    
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    
    eclipse_io = std::make_unique<EclipseIO>();
    gmsh_io = std::make_unique<GmshIO>();
    coord_manager = std::make_unique<CoordinateSystemManager>();
}

Simulator::~Simulator() {
    if (solution) VecDestroy(&solution);
    if (solution_old) VecDestroy(&solution_old);
    if (solution_prev_) VecDestroy(&solution_prev_);
    if (velocity_est_) VecDestroy(&velocity_est_);
    if (jacobian) MatDestroy(&jacobian);
    if (ts) TSDestroy(&ts);
    if (dm) DMDestroy(&dm);
}

PetscErrorCode Simulator::initialize(const SimulationConfig& config_in) {
    PetscFunctionBeginUser;
    
    config = config_in;
    current_time = config.start_time;
    dt = config.dt_initial;
    
    if (rank == 0) {
        PetscPrintf(comm, "Configuration:\n");
        PetscPrintf(comm, "  Start time: %g\n", config.start_time);
        PetscPrintf(comm, "  End time: %g\n", config.end_time);
        PetscPrintf(comm, "  Initial dt: %g\n", config.dt_initial);
        PetscPrintf(comm, "  Fluid model: %d\n", static_cast<int>(config.fluid_model));
        PetscPrintf(comm, "  Solid model: %d\n", static_cast<int>(config.solid_model));
    }
    
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::initializeFromConfigFile(const std::string& config_file) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    if (rank == 0) {
        PetscPrintf(comm, "Loading configuration from: %s\n", config_file.c_str());
    }
    
    ConfigReader reader;
    if (!reader.loadFile(config_file)) {
        SETERRQ(comm, PETSC_ERR_FILE_OPEN, "Failed to load configuration file");
    }
    
    // Parse simulation configuration
    reader.parseSimulationConfig(config);
    
    // Parse grid configuration
    reader.parseGridConfig(grid_config);
    
    // Parse material properties
    std::vector<MaterialProperties> props;
    if (reader.parseMaterialProperties(props)) {
        material_props = props;
    }
    
    // Parse fluid properties
    FluidProperties fluid;
    if (reader.parseFluidProperties(fluid)) {
        fluid_props.clear();
        fluid_props.push_back(fluid);
    }

    // Parse fault network (for induced seismicity / IMEX)
    if (config.enable_faults) {
        fault_network = reader.parseFaultNetwork();
    }

    // Parse IMEX configuration (if present)
    ConfigReader::IMEXConfig imex_cfg;
    if (reader.parseIMEXConfig(imex_cfg) && imex_cfg.enabled) {
        // IMEX manager is configured later after TS/DM exist; we only store for now.
        imex_config = imex_cfg;
    }

    // Parse seismicity config and catalog path (optional)
    ConfigReader::SeismicityConfig seis_cfg;
    if (reader.parseSeismicityConfig(seis_cfg)) {
        seismicity_config_ = seis_cfg;
        seismic_catalog_file_ = reader.getString("SEISMICITY", "catalog_file", "");
    }

    // Optional synthetic stress pulse to emulate an underground nuclear test trigger
    if (reader.hasSection("NUCLEAR_TRIGGER")) {
        nuclear_trigger_enabled_ = reader.getBool("NUCLEAR_TRIGGER", "enabled", true);
        nuclear_trigger_time0_ = reader.getDouble("NUCLEAR_TRIGGER", "time0", 10.0);
        nuclear_trigger_rise_ = reader.getDouble("NUCLEAR_TRIGGER", "rise_time", 0.5);
        nuclear_trigger_decay_ = reader.getDouble("NUCLEAR_TRIGGER", "decay_time", 60.0);
        nuclear_trigger_dsigma_ = reader.getDouble("NUCLEAR_TRIGGER", "delta_sigma_n", 0.0);
        nuclear_trigger_dtau_ = reader.getDouble("NUCLEAR_TRIGGER", "delta_tau", 2.0e6);
        nuclear_trigger_dpore_ = reader.getDouble("NUCLEAR_TRIGGER", "delta_pore_pressure", 0.0);
    }
    
    // Parse and configure wells
    auto well_cfgs = reader.parseWells();
    for (const auto& wc : well_cfgs) {
        WellType wt = WellType::PRODUCER;
        if (wc.type == "INJECTOR") wt = WellType::INJECTOR;
        if (wc.type == "PRODUCER") wt = WellType::PRODUCER;
        
        ierr = addWell(wc.name, wt); CHKERRQ(ierr);
        
        // Find well object and configure it
        for (auto& w : wells) {
            if (w->getName() != wc.name) continue;
            
            // Basic completion (single-cell completion)
            WellCompletion comp;
            comp.i = wc.i;
            comp.j = wc.j;
            comp.k = wc.k;
            comp.diameter = wc.diameter;
            comp.skin_factor = wc.skin;
            comp.is_open = true;
            comp.well_index = 1e-12; // default; can be overridden later
            w->addCompletion(comp);
            
            // Control mode
            if (wc.control_mode == "BHP") {
                w->setControl(WellControlMode::BHP_CONTROL, wc.target_value);
            } else if (wc.control_mode == "THP") {
                w->setControl(WellControlMode::THP_CONTROL, wc.target_value);
            } else if (wc.control_mode == "RESERVOIR_VOIDAGE") {
                w->setControl(WellControlMode::RESERVOIR_VOIDAGE, wc.target_value);
            } else {
                w->setControl(WellControlMode::RATE_CONTROL, wc.target_value);
            }
            
            // Wellbore model (legacy toggles)
            w->enableWellboreModel(wc.enable_wellbore_model);
            w->setWellboreRoughness(wc.tubing_roughness);
            w->setWellboreGeometry(wc.tubing_id * 0.5, wc.tubing_id * 0.5);
            
            // High-fidelity hydraulics configuration for THP/BHP conversion
            if (wc.enable_wellbore_model) {
                WellModel::WellboreHydraulicsOptions opts;
                opts.measured_depth = wc.wellbore_depth;
                opts.tubing_id = wc.tubing_id;
                opts.roughness = wc.tubing_roughness;
                opts.drag_reduction_factor = wc.drag_reduction_factor;
                opts.fluid_density = wc.wellbore_fluid_density;
                opts.fluid_viscosity = wc.wellbore_fluid_viscosity;
                opts.use_power_law = wc.use_power_law;
                opts.k_prime = wc.k_prime;
                opts.n_prime = wc.n_prime;
                
                std::string corr = wc.friction_correlation;
                std::transform(corr.begin(), corr.end(), corr.begin(), ::toupper);
                if (corr == "COLEBROOK_WHITE") {
                    opts.friction_correlation = FrictionModel::Correlation::COLEBROOK_WHITE;
                } else if (corr == "SWAMEE_JAIN") {
                    opts.friction_correlation = FrictionModel::Correlation::SWAMEE_JAIN;
                } else if (corr == "CHEN") {
                    opts.friction_correlation = FrictionModel::Correlation::CHEN;
                } else {
                    opts.friction_correlation = FrictionModel::Correlation::HAALAND;
                }
                
                w->configureWellboreHydraulics(opts);
            }
            
            break;
        }
    }
    
    // Initialize with loaded config
    ierr = initialize(config); CHKERRQ(ierr);
    
    if (rank == 0) {
        PetscPrintf(comm, "Configuration loaded successfully\n");
    }
    
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::loadMaterialPropertiesFromFile(const std::string& filename) {
    PetscFunctionBeginUser;
    
    ConfigReader reader;
    if (!reader.loadFile(filename)) {
        SETERRQ(comm, PETSC_ERR_FILE_OPEN, "Failed to load material properties file");
    }
    
    std::vector<MaterialProperties> props;
    if (reader.parseMaterialProperties(props)) {
        material_props = props;
    }
    
    FluidProperties fluid;
    if (reader.parseFluidProperties(fluid)) {
        fluid_props.clear();
        fluid_props.push_back(fluid);
    }
    
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::loadEclipseInput(const std::string& filename) {
    PetscFunctionBeginUser;
    
    if (rank == 0) {
        bool success = eclipse_io->readDeckFile(filename);
        if (!success) {
            SETERRQ(comm, PETSC_ERR_FILE_OPEN, "Failed to read Eclipse input file");
        }
        
        // Parse sections
        eclipse_io->parseRunspec();
        eclipse_io->parseGrid();
        eclipse_io->parseProps();
        eclipse_io->parseRegions();
        eclipse_io->parseSolution();
        eclipse_io->parseSchedule();
        eclipse_io->parseSummary();
    }
    
    // Get configuration from Eclipse file
    grid_config = eclipse_io->getGridConfig();
    
    if (rank == 0) {
        PetscPrintf(comm, "Grid: %d x %d x %d cells\n", 
                   grid_config.nx, grid_config.ny, grid_config.nz);
        PetscPrintf(comm, "Domain: %.1f x %.1f x %.1f m\n",
                   grid_config.Lx, grid_config.Ly, grid_config.Lz);
    }
    
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::setupDM() {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    // Setup coordinate system first
    ierr = setupCoordinateSystem(); CHKERRQ(ierr);
    
    // Choose grid setup based on mesh type
    switch (grid_config.mesh_type) {
        case MeshType::GMSH:
            ierr = setupGmshGrid(); CHKERRQ(ierr);
            break;
        case MeshType::CORNER_POINT:
        case MeshType::EXODUS:
        case MeshType::CUSTOM:
            if (grid_config.use_unstructured) {
                ierr = setupUnstructuredGrid(); CHKERRQ(ierr);
            } else {
                ierr = setupStructuredGrid(); CHKERRQ(ierr);
            }
            break;
        case MeshType::CARTESIAN:
        default:
            if (grid_config.use_unstructured) {
                ierr = setupUnstructuredGrid(); CHKERRQ(ierr);
            } else {
                ierr = setupStructuredGrid(); CHKERRQ(ierr);
            }
            break;
    }
    
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::setupCoordinateSystem() {
    PetscFunctionBeginUser;
    
    // Setup coordinate reference system if specified
    if (!grid_config.input_crs.empty()) {
        coord_manager->setInputCRS(grid_config.input_crs);
        
        if (!grid_config.model_crs.empty()) {
            coord_manager->setModelCRS(grid_config.model_crs);
        } else if (grid_config.auto_detect_utm) {
            // Auto-detect UTM zone based on local origin
            GeoPoint origin(grid_config.local_origin_x, grid_config.local_origin_y, grid_config.local_origin_z);
            coord_manager->setLocalOrigin(origin);
            coord_manager->useAutoLocalCRS();
        }
        
        if (grid_config.use_local_coordinates) {
            GeoPoint origin(grid_config.local_origin_x, grid_config.local_origin_y, grid_config.local_origin_z);
            coord_manager->setLocalOrigin(origin);
        }
        
        if (rank == 0 && coord_manager->isConfigured()) {
            PetscPrintf(comm, "Coordinate System:\n");
            PetscPrintf(comm, "  Input CRS: %s\n", grid_config.input_crs.c_str());
            PetscPrintf(comm, "  Model CRS: %s\n", coord_manager->getModelCRS().epsg_code.c_str());
            if (grid_config.use_local_coordinates) {
                PetscPrintf(comm, "  Local origin: (%.2f, %.2f, %.2f)\n",
                           grid_config.local_origin_x, grid_config.local_origin_y, grid_config.local_origin_z);
            }
        }
    }
    
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::setupGmshGrid() {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    if (grid_config.mesh_file.empty()) {
        SETERRQ(comm, PETSC_ERR_ARG_NULL, "Gmsh mesh file not specified");
    }
    
    if (rank == 0) {
        PetscPrintf(comm, "Loading Gmsh mesh: %s\n", grid_config.mesh_file.c_str());
    }
    
    // Load the Gmsh file
    if (!gmsh_io->readMshFile(grid_config.mesh_file)) {
        SETERRQ(comm, PETSC_ERR_FILE_OPEN, "Failed to read Gmsh mesh file");
    }
    
    if (rank == 0) {
        gmsh_io->printStatistics();
    }
    
    // Validate mesh
    if (!gmsh_io->validate()) {
        SETERRQ(comm, PETSC_ERR_ARG_WRONG, "Gmsh mesh validation failed");
    }

    // Configure domain mappings (materials / faults / boundaries) if provided.
    // Note: mappings are parsed into grid_config by ConfigReader::parseGridConfig().
    bool have_mappings = !grid_config.gmsh_material_domains.empty() || !grid_config.gmsh_fault_surfaces.empty();
    if (have_mappings) {
        for (const auto& [physical_group, material_section] : grid_config.gmsh_material_domains) {
            gmsh_io->addMaterialMapping(physical_group, material_section);
        }
        for (const auto& [physical_group, fault_section, use_split] : grid_config.gmsh_fault_surfaces) {
            gmsh_io->addFaultMapping(physical_group, fault_section, use_split);
        }
    }
    if (!grid_config.gmsh_boundaries.empty()) {
        gmsh_io->setBoundaryNames(grid_config.gmsh_boundaries);
        have_mappings = true;
    }

    if (have_mappings) {
        if (!gmsh_io->processDomains()) {
            SETERRQ(comm, PETSC_ERR_ARG_WRONG, "Failed to process Gmsh domain mappings");
        }
        ierr = gmsh_io->createDMPlexFull(comm, &dm); CHKERRQ(ierr);
    } else {
        ierr = gmsh_io->createDMPlex(comm, &dm); CHKERRQ(ierr);
    }
    
    // Update grid config with actual mesh dimensions
    auto mesh = gmsh_io->getMesh();
    if (mesh) {
        mesh->computeBoundingBox();
        grid_config.Lx = mesh->getLx();
        grid_config.Ly = mesh->getLy();
        grid_config.Lz = mesh->getLz();
        
        // Apply coordinate transformation if configured
        if (coord_manager->isConfigured()) {
            // Transform the mesh coordinates
            // Note: This would require modifying the DM coordinates which is done below
        }
    }
    
    ierr = DMSetFromOptions(dm); CHKERRQ(ierr);
    ierr = DMViewFromOptions(dm, nullptr, "-dm_view"); CHKERRQ(ierr);
    
    if (rank == 0) {
        PetscPrintf(comm, "Gmsh mesh loaded successfully\n");
        PetscPrintf(comm, "  Domain size: %.1f x %.1f x %.1f m\n",
                   grid_config.Lx, grid_config.Ly, grid_config.Lz);
    }
    
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::setupStructuredGrid() {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    // Create DMPlex for structured grid
    PetscInt faces[3] = {grid_config.nx, grid_config.ny, grid_config.nz};
    PetscReal lower[3] = {0.0, 0.0, 0.0};
    PetscReal upper[3] = {grid_config.Lx, grid_config.Ly, grid_config.Lz};
    ierr = DMPlexCreateBoxMesh(comm, 3, PETSC_FALSE, 
                              faces,
                              lower, upper, nullptr, PETSC_TRUE, &dm); CHKERRQ(ierr);
    
    ierr = DMSetFromOptions(dm); CHKERRQ(ierr);
    ierr = DMSetUp(dm); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::setupUnstructuredGrid() {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    // Create DMPlex for unstructured grid
    PetscInt faces[3] = {grid_config.nx, grid_config.ny, grid_config.nz};
    ierr = DMPlexCreateBoxMesh(comm, 3, PETSC_FALSE, 
                              faces,
                              nullptr, nullptr, nullptr, PETSC_TRUE, &dm); CHKERRQ(ierr);
    
    ierr = DMSetFromOptions(dm); CHKERRQ(ierr);
    ierr = DMViewFromOptions(dm, nullptr, "-dm_view"); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::setupFields() {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    ierr = createFieldsFromConfig(); CHKERRQ(ierr);
    
    // Create solution vectors
    ierr = DMCreateGlobalVector(dm, &solution); CHKERRQ(ierr);
    ierr = VecDuplicate(solution, &solution_old); CHKERRQ(ierr);
    
    ierr = PetscObjectSetName((PetscObject)solution, "solution"); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::createFieldsFromConfig() {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    // Create discretization based on fluid model
    PetscFE fe;
    
    switch (config.fluid_model) {
        case FluidModelType::SINGLE_COMPONENT:
            // Single pressure field
            ierr = PetscFECreateDefault(comm, 3, 1, PETSC_FALSE, "pressure_", -1, &fe); CHKERRQ(ierr);
            ierr = DMAddField(dm, nullptr, (PetscObject)fe); CHKERRQ(ierr);
            fe_fields.push_back(fe);
            break;
            
        case FluidModelType::BLACK_OIL:
            // Pressure + saturations
            ierr = PetscFECreateDefault(comm, 3, 1, PETSC_FALSE, "pressure_", -1, &fe); CHKERRQ(ierr);
            ierr = DMAddField(dm, nullptr, (PetscObject)fe); CHKERRQ(ierr);
            fe_fields.push_back(fe);
            
            ierr = PetscFECreateDefault(comm, 3, 1, PETSC_FALSE, "saturation_w_", -1, &fe); CHKERRQ(ierr);
            ierr = DMAddField(dm, nullptr, (PetscObject)fe); CHKERRQ(ierr);
            fe_fields.push_back(fe);
            
            ierr = PetscFECreateDefault(comm, 3, 1, PETSC_FALSE, "saturation_g_", -1, &fe); CHKERRQ(ierr);
            ierr = DMAddField(dm, nullptr, (PetscObject)fe); CHKERRQ(ierr);
            fe_fields.push_back(fe);
            break;
            
        case FluidModelType::COMPOSITIONAL:
            // Pressure + compositions
            ierr = PetscFECreateDefault(comm, 3, 1, PETSC_FALSE, "pressure_", -1, &fe); CHKERRQ(ierr);
            ierr = DMAddField(dm, nullptr, (PetscObject)fe); CHKERRQ(ierr);
            fe_fields.push_back(fe);
            
            // Add composition fields (number depends on components)
            for (int i = 0; i < 3; ++i) {  // Example: 3 components
                char name[256];
                snprintf(name, sizeof(name), "composition_%d_", i);
                ierr = PetscFECreateDefault(comm, 3, 1, PETSC_FALSE, name, -1, &fe); CHKERRQ(ierr);
                ierr = DMAddField(dm, nullptr, (PetscObject)fe); CHKERRQ(ierr);
                fe_fields.push_back(fe);
            }
            break;
            
        default:
            SETERRQ(comm, PETSC_ERR_ARG_WRONG, "Unknown fluid model type");
    }
    
    // Add geomechanics field if enabled
    if (config.enable_geomechanics) {
        ierr = PetscFECreateDefault(comm, 3, 3, PETSC_FALSE, "displacement_", -1, &fe); CHKERRQ(ierr);
        ierr = DMAddField(dm, nullptr, (PetscObject)fe); CHKERRQ(ierr);
        fe_fields.push_back(fe);
    }
    
    // Add thermal field if enabled
    if (config.enable_thermal) {
        ierr = PetscFECreateDefault(comm, 3, 1, PETSC_FALSE, "temperature_", -1, &fe); CHKERRQ(ierr);
        ierr = DMAddField(dm, nullptr, (PetscObject)fe); CHKERRQ(ierr);
        fe_fields.push_back(fe);
    }
    
    // Create DS
    ierr = DMCreateDS(dm); CHKERRQ(ierr);
    ierr = DMGetDS(dm, &prob); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::setupPhysics() {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    // Add physics kernels based on configuration
    switch (config.fluid_model) {
        case FluidModelType::SINGLE_COMPONENT: {
            auto kernel = std::make_shared<SinglePhaseFlowKernel>();
            ierr = addPhysicsKernel(kernel); CHKERRQ(ierr);
            break;
        }
        case FluidModelType::BLACK_OIL: {
            auto kernel = std::make_shared<BlackOilKernel>();
            ierr = addPhysicsKernel(kernel); CHKERRQ(ierr);
            break;
        }
        case FluidModelType::COMPOSITIONAL: {
            auto kernel = std::make_shared<CompositionalKernel>(3); // 3 components
            ierr = addPhysicsKernel(kernel); CHKERRQ(ierr);
            break;
        }
        default:
            break;
    }
    
    // Add geomechanics if enabled
    if (config.enable_geomechanics) {
        // Choose between static and dynamic geomechanics
        if (config.solid_model == SolidModelType::ELASTODYNAMIC) {
            auto kernel = std::make_shared<ElastodynamicsKernel>();
            if (!material_props.empty()) {
                const MaterialProperties& mat = material_props[0];
                kernel->setMaterialProperties(mat.youngs_modulus, mat.poisson_ratio, mat.density);
                kernel->setWaveProperties(mat.p_wave_velocity, mat.s_wave_velocity, mat.quality_factor);
                kernel->setDamping(mat.damping_alpha, mat.damping_beta);
                
                if (config.use_static_triggering) {
                    kernel->setStaticTriggeringMode(true, config.dynamic_trigger_threshold, 
                                                   config.dynamic_event_duration);
                }
            }
            ierr = addPhysicsKernel(kernel); CHKERRQ(ierr);
        } else if (config.solid_model == SolidModelType::POROELASTODYNAMIC) {
            auto kernel = std::make_shared<PoroelastodynamicsKernel>();
            if (!material_props.empty()) {
                const MaterialProperties& mat = material_props[0];
                const FluidProperties& fluid = fluid_props.empty() ? FluidProperties() : fluid_props[0];
                
                kernel->setMaterialProperties(mat.youngs_modulus, mat.poisson_ratio, 
                                             mat.density, mat.porosity);
                kernel->setFluidProperties(fluid.density, fluid.viscosity, 2.2e9);
                kernel->setBiotParameters(mat.biot_coefficient, 1e10);
                kernel->setWaveProperties(mat.p_wave_velocity, mat.s_wave_velocity, 1500.0);
                kernel->setDamping(mat.damping_alpha, mat.damping_beta);
                
                if (config.use_static_triggering) {
                    kernel->setStaticTriggeringMode(true, config.dynamic_trigger_threshold,
                                                   config.dynamic_event_duration);
                }
                
                if (config.enable_dynamic_permeability_change) {
                    kernel->enableDynamicPermeabilityChange(true, 
                                                           mat.permeability_strain_coeff,
                                                           mat.permeability_stress_coeff,
                                                           config.permeability_recovery_time);
                }
            }
            ierr = addPhysicsKernel(kernel); CHKERRQ(ierr);
        } else {
            auto kernel = std::make_shared<GeomechanicsKernel>(config.solid_model);
            ierr = addPhysicsKernel(kernel); CHKERRQ(ierr);
        }
    }
    
    // Add standalone elastodynamics if enabled
    if (config.enable_elastodynamics) {
        auto kernel = std::make_shared<ElastodynamicsKernel>();
        if (!material_props.empty()) {
            const MaterialProperties& mat = material_props[0];
            kernel->setMaterialProperties(mat.youngs_modulus, mat.poisson_ratio, mat.density);
            kernel->setWaveProperties(mat.p_wave_velocity, mat.s_wave_velocity, mat.quality_factor);
            kernel->setDamping(mat.damping_alpha, mat.damping_beta);
            
            if (config.use_static_triggering) {
                kernel->setStaticTriggeringMode(true, config.dynamic_trigger_threshold,
                                               config.dynamic_event_duration);
            }
        }
        ierr = addPhysicsKernel(kernel); CHKERRQ(ierr);
    }
    
    // Add standalone poroelastodynamics if enabled
    if (config.enable_poroelastodynamics) {
        auto kernel = std::make_shared<PoroelastodynamicsKernel>();
        if (!material_props.empty()) {
            const MaterialProperties& mat = material_props[0];
            const FluidProperties& fluid = fluid_props.empty() ? FluidProperties() : fluid_props[0];
            
            kernel->setMaterialProperties(mat.youngs_modulus, mat.poisson_ratio,
                                         mat.density, mat.porosity);
            kernel->setFluidProperties(fluid.density, fluid.viscosity, 2.2e9);
            kernel->setBiotParameters(mat.biot_coefficient, 1e10);
            kernel->setWaveProperties(mat.p_wave_velocity, mat.s_wave_velocity, 1500.0);
            kernel->setDamping(mat.damping_alpha, mat.damping_beta);
            
            if (config.use_static_triggering) {
                kernel->setStaticTriggeringMode(true, config.dynamic_trigger_threshold,
                                               config.dynamic_event_duration);
            }
            
            if (config.enable_dynamic_permeability_change) {
                kernel->enableDynamicPermeabilityChange(true,
                                                       mat.permeability_strain_coeff,
                                                       mat.permeability_stress_coeff,
                                                       config.permeability_recovery_time);
            }
        }
        ierr = addPhysicsKernel(kernel); CHKERRQ(ierr);
    }
    
    // Add thermal if enabled
    if (config.enable_thermal) {
        auto kernel = std::make_shared<ThermalKernel>();
        ierr = addPhysicsKernel(kernel); CHKERRQ(ierr);
    }
    
    // Add particle transport if enabled
    if (config.enable_particle_transport) {
        auto kernel = std::make_shared<ParticleTransportKernel>();
        ierr = addPhysicsKernel(kernel); CHKERRQ(ierr);
    }
    
    // Setup residual and Jacobian functions
    ierr = PetscDSSetResidual(prob, 0, f0_SinglePhase, f1_SinglePhase); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::addPhysicsKernel(std::shared_ptr<PhysicsKernel> kernel) {
    PetscFunctionBeginUser;
    
    kernels.push_back(kernel);
    kernel_map[kernel->getType()] = kernels.size() - 1;
    
    // Setup kernel with DM and FE
    if (fe_fields.size() > 0) {
        kernel->setup(dm, fe_fields[0]);
    }
    
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::setupTimeStepper() {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    // Create time stepper
    ierr = TSCreate(comm, &ts); CHKERRQ(ierr);
    ierr = TSSetDM(ts, dm); CHKERRQ(ierr);
    ierr = TSSetProblemType(ts, TS_NONLINEAR); CHKERRQ(ierr);
    
    // Set time stepping method
    // Use implicit methods for quasi-static, explicit or implicit for dynamic
    if (config.use_dynamic_mode || config.enable_elastodynamics || config.enable_poroelastodynamics) {
        // For wave propagation, use second-order time integrator
        ierr = TSSetType(ts, TSALPHA2); CHKERRQ(ierr);  // Generalized-alpha for dynamics
        
        // Set parameters for optimal stability and accuracy
        PetscReal alpha_m = 0.5;  // Spectral radius
        PetscReal alpha_f = 0.5;
        // PetscReal gamma = 0.5 + alpha_m - alpha_f;  // Not used in current implementation
        // These can be set via TSAlpha2SetParams if needed
    } else {
        // Quasi-static: use backward Euler
        ierr = TSSetType(ts, TSBEULER); CHKERRQ(ierr);
    }
    
    // Set time parameters
    ierr = TSSetTime(ts, config.start_time); CHKERRQ(ierr);
    ierr = TSSetMaxTime(ts, config.end_time); CHKERRQ(ierr);
    ierr = TSSetTimeStep(ts, config.dt_initial); CHKERRQ(ierr);
    ierr = TSSetMaxSteps(ts, config.max_timesteps); CHKERRQ(ierr);
    ierr = TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP); CHKERRQ(ierr);
    
    // Set residual and Jacobian functions
    ierr = TSSetIFunction(ts, nullptr, FormFunction, this); CHKERRQ(ierr);
    ierr = TSSetIJacobian(ts, nullptr, nullptr, FormJacobian, this); CHKERRQ(ierr);
    
    // Set monitor (may be overridden by IMEX setup, see below)
    ierr = TSMonitorSet(ts, MonitorFunction, this, nullptr); CHKERRQ(ierr);

    // If IMEX was configured during config parsing, initialize it now (TS/DM exist).
    if (imex_config.enabled) {
        ierr = setupIMEXTransition(imex_config); CHKERRQ(ierr);
        // IMEX setup registers its own monitor; restore our monitor so we can:
        // - write outputs
        // - update fault network
        // - call imex_manager->checkAndTransition() ourselves
        ierr = TSMonitorSet(ts, MonitorFunction, this, nullptr); CHKERRQ(ierr);
    }
    
    // Set from options
    ierr = TSSetFromOptions(ts); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::setupSolvers() {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    // Get SNES from TS
    ierr = TSGetSNES(ts, &snes); CHKERRQ(ierr);
    
    // Set tolerances
    ierr = SNESSetTolerances(snes, config.atol, config.rtol, 
                            PETSC_DEFAULT, config.max_nonlinear_iterations, 
                            PETSC_DEFAULT); CHKERRQ(ierr);
    
    // Get KSP from SNES
    ierr = SNESGetKSP(snes, &ksp); CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp, config.rtol, config.atol, 
                           PETSC_DEFAULT, config.max_linear_iterations); CHKERRQ(ierr);
    
    // Get preconditioner
    ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
    
    // Set solver options from command line
    ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
    ierr = PCSetFromOptions(pc); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::setupIMEXTransition(const ConfigReader::IMEXConfig& imex_cfg) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    if (!imex_cfg.enabled) {
        if (rank == 0) {
            PetscPrintf(comm, "IMEX transition: disabled\n");
        }
        PetscFunctionReturn(0);
    }
    
    if (rank == 0) {
        PetscPrintf(comm, "Setting up IMEX (Implicit-Explicit) time integration transition...\n");
        PetscPrintf(comm, "  Initial mode: %s\n", imex_cfg.initial_mode.c_str());
        PetscPrintf(comm, "  Trigger type: %s\n", imex_cfg.trigger_type.c_str());
        PetscPrintf(comm, "  Settling type: %s\n", imex_cfg.settling_type.c_str());
    }
    
    // Store config
    imex_config = imex_cfg;
    
    // Create IMEX manager
    imex_manager = std::make_unique<ImplicitExplicitTransitionManager>(comm);
    
    // Convert string config to internal IMEX config structure
    IMEXInternalConfig internal_config;
    internal_config.enabled = imex_cfg.enabled;
    internal_config.initial_mode = (imex_cfg.initial_mode == "EXPLICIT") 
                                   ? IntegrationMode::EXPLICIT 
                                   : IntegrationMode::IMPLICIT;
    
    // Time step settings
    internal_config.implicit_dt_initial = imex_cfg.implicit_dt_initial;
    internal_config.implicit_dt_min = imex_cfg.implicit_dt_min;
    internal_config.implicit_dt_max = imex_cfg.implicit_dt_max;
    internal_config.implicit_method = imex_cfg.implicit_method;
    
    internal_config.explicit_dt_initial = imex_cfg.explicit_dt_initial;
    internal_config.explicit_dt_min = imex_cfg.explicit_dt_min;
    internal_config.explicit_dt_max = imex_cfg.explicit_dt_max;
    internal_config.explicit_method = imex_cfg.explicit_method;
    internal_config.cfl_factor = imex_cfg.cfl_factor;
    
    // Trigger settings
    if (imex_cfg.trigger_type == "STRESS_THRESHOLD" || imex_cfg.trigger_type == "STRESS") {
        internal_config.trigger_type = TransitionTrigger::STRESS_THRESHOLD;
    } else if (imex_cfg.trigger_type == "COULOMB_FAILURE" || imex_cfg.trigger_type == "COULOMB") {
        internal_config.trigger_type = TransitionTrigger::COULOMB_FAILURE;
    } else if (imex_cfg.trigger_type == "SLIP_RATE") {
        internal_config.trigger_type = TransitionTrigger::SLIP_RATE;
    } else if (imex_cfg.trigger_type == "VELOCITY") {
        internal_config.trigger_type = TransitionTrigger::VELOCITY;
    } else if (imex_cfg.trigger_type == "ACCELERATION") {
        internal_config.trigger_type = TransitionTrigger::ACCELERATION;
    } else if (imex_cfg.trigger_type == "ENERGY_RATE") {
        internal_config.trigger_type = TransitionTrigger::ENERGY_RATE;
    } else {
        internal_config.trigger_type = TransitionTrigger::COULOMB_FAILURE;
    }
    
    internal_config.stress_threshold = imex_cfg.stress_threshold;
    internal_config.coulomb_threshold = imex_cfg.coulomb_threshold;
    internal_config.slip_rate_threshold = imex_cfg.slip_rate_threshold;
    internal_config.velocity_threshold = imex_cfg.velocity_threshold;
    internal_config.acceleration_threshold = imex_cfg.acceleration_threshold;
    internal_config.energy_rate_threshold = imex_cfg.energy_rate_threshold;
    
    // Settling settings
    if (imex_cfg.settling_type == "TIME_ELAPSED" || imex_cfg.settling_type == "TIME") {
        internal_config.settling_type = SettlingCriterion::TIME_ELAPSED;
    } else if (imex_cfg.settling_type == "VELOCITY_DECAY" || imex_cfg.settling_type == "VELOCITY") {
        internal_config.settling_type = SettlingCriterion::VELOCITY_DECAY;
    } else if (imex_cfg.settling_type == "ENERGY_DECAY" || imex_cfg.settling_type == "ENERGY") {
        internal_config.settling_type = SettlingCriterion::ENERGY_DECAY;
    } else if (imex_cfg.settling_type == "SLIP_RATE_DECAY" || imex_cfg.settling_type == "SLIP_RATE") {
        internal_config.settling_type = SettlingCriterion::SLIP_RATE_DECAY;
    } else if (imex_cfg.settling_type == "COMBINED") {
        internal_config.settling_type = SettlingCriterion::COMBINED;
    } else {
        internal_config.settling_type = SettlingCriterion::COMBINED;
    }
    
    internal_config.min_dynamic_duration = imex_cfg.min_dynamic_duration;
    internal_config.max_dynamic_duration = imex_cfg.max_dynamic_duration;
    internal_config.settling_velocity = imex_cfg.settling_velocity;
    internal_config.settling_energy_ratio = imex_cfg.settling_energy_ratio;
    internal_config.settling_slip_rate = imex_cfg.settling_slip_rate;
    internal_config.settling_observation_window = imex_cfg.settling_observation_window;
    
    // Solver settings
    internal_config.adaptive_timestep = imex_cfg.adaptive_timestep;
    internal_config.dt_growth_factor = imex_cfg.dt_growth_factor;
    internal_config.dt_shrink_factor = imex_cfg.dt_shrink_factor;
    internal_config.implicit_max_iterations = imex_cfg.implicit_max_iterations;
    internal_config.explicit_max_iterations = imex_cfg.explicit_max_iterations;
    internal_config.implicit_rtol = imex_cfg.implicit_rtol;
    internal_config.explicit_rtol = imex_cfg.explicit_rtol;
    internal_config.implicit_atol = imex_cfg.implicit_atol;
    internal_config.explicit_atol = imex_cfg.explicit_atol;
    
    // Output settings
    internal_config.use_lumped_mass = imex_cfg.use_lumped_mass;
    internal_config.implicit_output_frequency = imex_cfg.implicit_output_frequency;
    internal_config.explicit_output_frequency = imex_cfg.explicit_output_frequency;
    internal_config.log_transitions = imex_cfg.log_transitions;
    internal_config.transition_log_file = imex_cfg.transition_log_file;
    internal_config.smooth_transition = imex_cfg.smooth_transition;
    internal_config.transition_ramp_time = imex_cfg.transition_ramp_time;
    
    // Initialize IMEX manager
    ierr = imex_manager->initialize(internal_config); CHKERRQ(ierr);
    
    // Setup with TS and solution
    ierr = imex_manager->setup(ts, dm, solution, nullptr); CHKERRQ(ierr);
    
    // Connect fault network if available
    if (fault_network) {
        imex_manager->setFaultNetwork(fault_network.get());
    }
    
    if (rank == 0) {
        PetscPrintf(comm, "IMEX transition manager initialized successfully\n");
    }
    
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::setupFaultNetwork() {
    PetscFunctionBeginUser;
    
    if (!config.enable_faults) {
        PetscFunctionReturn(0);
    }
    
    if (rank == 0) {
        PetscPrintf(comm, "Setting up fault network for induced seismicity...\n");
    }
    
    fault_network = std::make_unique<FaultNetwork>();
    
    // Fault network will be populated from config file parsing
    // This is called after config parsing to finalize setup
    
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::setInitialConditions() {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    // Set initial conditions from Eclipse data
    Vec local;
    ierr = DMGetLocalVector(dm, &local); CHKERRQ(ierr);
    
    // Set uniform initial condition (for now)
    ierr = VecSet(solution, 1.0e7); CHKERRQ(ierr); // 1 bar initial pressure
    
    ierr = DMRestoreLocalVector(dm, &local); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::setMaterialProperties() {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    // If material_props already loaded from config file, use those
    if (!material_props.empty()) {
        ierr = applyMaterialPropertiesToKernels(); CHKERRQ(ierr);
        PetscFunctionReturn(0);
    }
    
    // Otherwise, get material properties from Eclipse input for each cell
    int ncells = grid_config.nx * grid_config.ny * grid_config.nz;
    material_props.resize(ncells);
    
    for (int k = 0; k < grid_config.nz; ++k) {
        for (int j = 0; j < grid_config.ny; ++j) {
            for (int i = 0; i < grid_config.nx; ++i) {
                int idx = i + j * grid_config.nx + k * grid_config.nx * grid_config.ny;
                material_props[idx] = eclipse_io->getMaterialProperties(i, j, k);
            }
        }
    }
    
    ierr = applyMaterialPropertiesToKernels(); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::applyMaterialPropertiesToKernels() {
    PetscFunctionBeginUser;
    
    if (material_props.empty()) {
        PetscFunctionReturn(0);
    }
    
    // Use first material properties (can be extended for heterogeneous)
    const MaterialProperties& mat = material_props[0];
    const FluidProperties& fluid = fluid_props.empty() ? FluidProperties() : fluid_props[0];
    
    // Apply to physics kernels
    for (auto& kernel : kernels) {
        switch (kernel->getType()) {
            case PhysicsType::FLUID_FLOW: {
                if (auto sp_kernel = std::dynamic_pointer_cast<SinglePhaseFlowKernel>(kernel)) {
                    sp_kernel->setProperties(mat.porosity, mat.permeability_x, 
                                            mat.compressibility, fluid.viscosity, fluid.density);
                } else if (auto bo_kernel = std::dynamic_pointer_cast<BlackOilKernel>(kernel)) {
                    bo_kernel->setRockProperties(mat.porosity, mat.permeability_x,
                                                mat.permeability_y, mat.permeability_z);
                    bo_kernel->setFluidProperties(fluid);
                }
                break;
            }
            case PhysicsType::GEOMECHANICS: {
                if (auto geo_kernel = std::dynamic_pointer_cast<GeomechanicsKernel>(kernel)) {
                    geo_kernel->setMaterialProperties(mat.youngs_modulus, mat.poisson_ratio, mat.density);
                    if (config.solid_model == SolidModelType::VISCOELASTIC) {
                        geo_kernel->setViscoelasticProperties(mat.relaxation_time, mat.viscosity);
                    }
                    if (config.solid_model == SolidModelType::POROELASTIC) {
                        geo_kernel->setPoroelasticCoupling(mat.biot_coefficient, 1e10);
                    }
                }
                break;
            }
            case PhysicsType::THERMAL: {
                if (auto therm_kernel = std::dynamic_pointer_cast<ThermalKernel>(kernel)) {
                    therm_kernel->setThermalProperties(mat.thermal_conductivity, 
                                                      mat.density, mat.heat_capacity);
                }
                break;
            }
            default:
                break;
        }
    }
    
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::run() {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    // Set initial time
    ierr = TSSetSolution(ts, solution); CHKERRQ(ierr);
    
    // Solve
    ierr = TSSolve(ts, solution); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

// Static callback functions
PetscErrorCode Simulator::FormFunction(TS ts, PetscReal t, Vec U, Vec U_t, Vec F, void *ctx) {
    (void)ts;  // Part of PETSc callback interface - accessed via ctx
    
    PetscFunctionBeginUser;
    Simulator *sim = static_cast<Simulator*>(ctx);
    PetscErrorCode ierr;
    
    // Zero out residual
    ierr = VecSet(F, 0.0); CHKERRQ(ierr);
    
    // Compute residual from DMPlex
    ierr = DMPlexTSComputeIFunctionFEM(sim->dm, t, U, U_t, F, ctx); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::FormJacobian(TS ts, PetscReal t, Vec U, Vec U_t, 
                                       PetscReal a, Mat J, Mat P, void *ctx) {
    (void)ts;  // Part of PETSc callback interface - accessed via ctx
    
    PetscFunctionBeginUser;
    Simulator *sim = static_cast<Simulator*>(ctx);
    PetscErrorCode ierr;
    
    // Compute Jacobian from DMPlex
    ierr = DMPlexTSComputeIJacobianFEM(sim->dm, t, U, U_t, a, J, P, ctx); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::MonitorFunction(TS ts, PetscInt step, PetscReal t, 
                                          Vec U, void *ctx) {
    (void)ts;  // Part of PETSc callback interface - accessed via ctx
    (void)U;   // Solution available through ctx->solution
    
    PetscFunctionBeginUser;
    Simulator *sim = static_cast<Simulator*>(ctx);
    PetscErrorCode ierr;
    
    // Print progress
    if (sim->rank == 0 && step % 10 == 0) {
        PetscPrintf(sim->comm, "Step %d, Time = %g\n", (int)step, (double)t);
    }
    
    // Write output if necessary
    if (step % sim->config.output_frequency == 0) {
        ierr = sim->writeOutput(step); CHKERRQ(ierr);
    }

    // Estimate velocity from finite difference (needed by IMEX manager)
    if (sim->imex_manager) {
        if (!sim->solution_prev_) {
            ierr = VecDuplicate(U, &sim->solution_prev_); CHKERRQ(ierr);
            ierr = VecCopy(U, sim->solution_prev_); CHKERRQ(ierr);
            sim->prev_time_ = static_cast<double>(t);
        } else {
            double dt = static_cast<double>(t) - sim->prev_time_;
            if (dt <= 0.0) dt = sim->config.dt_initial;
            if (!sim->velocity_est_) {
                ierr = VecDuplicate(U, &sim->velocity_est_); CHKERRQ(ierr);
            }
            // velocity_est = (U - U_prev) / dt
            ierr = VecCopy(U, sim->velocity_est_); CHKERRQ(ierr);
            ierr = VecAXPY(sim->velocity_est_, -1.0, sim->solution_prev_); CHKERRQ(ierr);
            ierr = VecScale(sim->velocity_est_, 1.0 / dt); CHKERRQ(ierr);
            ierr = VecCopy(U, sim->solution_prev_); CHKERRQ(ierr);
            sim->prev_time_ = static_cast<double>(t);
        }
    }

    // Update fault network with a simplified background stress + optional nuclear trigger pulse.
    if (sim->fault_network && sim->seismicity_config_.enable_seismicity) {
        double dt = sim->config.dt_initial;
        if (sim->last_fault_update_time_ >= 0.0) {
            dt = static_cast<double>(t) - sim->last_fault_update_time_;
            if (dt <= 0.0) dt = sim->config.dt_initial;
        }
        sim->last_fault_update_time_ = static_cast<double>(t);

        // Background stress model: simple lithostatic + horizontal ratio.
        const double rho = 2700.0;
        const double g = 9.81;
        // Use median fault depth if available
        double depth = 3000.0;
        if (sim->fault_network->numFaults() > 0) {
            auto* f0 = sim->fault_network->getFault(0);
            if (f0) depth = std::max(1.0, f0->getGeometry().z);
        }
        const double szz = rho * g * depth;          // Pa
        const double k0 = 0.7;                       // Shmin/Shv ratio (typical)
        const double sxx = k0 * szz;
        const double syy = k0 * szz;

        double sxy = 0.0, sxz = 0.0, syz = 0.0;
        double pore = 0.0;

        // Optional nuclear trigger pulse: apply a transient shear/normal stress change.
        if (sim->nuclear_trigger_enabled_) {
            const double tt = static_cast<double>(t);
            const double t0 = sim->nuclear_trigger_time0_;
            if (tt >= t0) {
                const double tr = std::max(1e-6, sim->nuclear_trigger_rise_);
                const double td = std::max(1e-6, sim->nuclear_trigger_decay_);
                // Smooth rise then exponential decay
                double rise = 1.0 - std::exp(-(tt - t0) / tr);
                double decay = std::exp(-(tt - t0) / td);
                double amp = rise * decay;
                // Apply as shear stress on xz component (illustrative)
                sxz += amp * sim->nuclear_trigger_dtau_;
                // Apply as normal stress change on vertical (illustrative)
                // Positive = more compression; negative could unclamp
                // Here we use the provided delta_sigma_n directly.
                // (Applied to szz as simplest proxy.)
                (void)szz;
                pore += amp * sim->nuclear_trigger_dpore_;
                // Note: delta_sigma_n is not directly applied to full tensor here;
                // it is available for future extensions.
            }
        }

        const double shear_modulus = 30.0e9;  // Pa (order-of-magnitude)
        if (sim->seismicity_config_.enable_stress_transfer) {
            sim->fault_network->enableStressTransfer(true);
            sim->fault_network->updateWithCascade(sxx, syy, szz, sxy, sxz, syz, pore, shear_modulus, dt);
        } else {
            sim->fault_network->update(sxx, syy, szz, sxy, sxz, syz, pore, shear_modulus, dt);
        }
    }

    // Let IMEX manager check for transitions (requires a velocity vector)
    if (sim->imex_manager && sim->velocity_est_) {
        ierr = sim->imex_manager->checkAndTransition(U, sim->velocity_est_, static_cast<double>(t)); CHKERRQ(ierr);
    }
    
    PetscFunctionReturn(0);
}

// Pointwise residual functions
void Simulator::f0_SinglePhase(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                               const PetscInt uOff[], const PetscInt uOff_x[],
                               const PetscScalar u[], const PetscScalar u_t[],
                               const PetscScalar u_x[], const PetscInt aOff[],
                               const PetscInt aOff_x[], const PetscScalar a[],
                               const PetscScalar a_x[], const PetscScalar a_t[],
                               PetscReal t,
                               const PetscReal x[], PetscInt numConstants,
                               const PetscScalar constants[], PetscScalar f0[]) {
    // Suppress unused parameter warnings - part of PETSc pointwise function interface
    (void)dim; (void)Nf; (void)NfAux;
    (void)uOff; (void)uOff_x;
    (void)u; (void)u_x;
    (void)aOff; (void)aOff_x;
    (void)a; (void)a_x; (void)a_t;
    (void)t; (void)x;
    (void)numConstants; (void)constants;
    
    // Accumulation term: phi * ct * dP/dt
    const PetscReal phi = 0.2;  // porosity
    const PetscReal ct = 1.0e-9; // compressibility
    
    f0[0] = phi * ct * u_t[0];
}

void Simulator::f1_SinglePhase(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                               const PetscInt uOff[], const PetscInt uOff_x[],
                               const PetscScalar u[], const PetscScalar u_t[],
                               const PetscScalar u_x[], const PetscInt aOff[],
                               const PetscInt aOff_x[], const PetscScalar a[],
                               const PetscScalar a_x[], const PetscScalar a_t[],
                               PetscReal t,
                               const PetscReal x[], PetscInt numConstants,
                               const PetscScalar constants[], PetscScalar f1[]) {
    // Suppress unused parameter warnings - part of PETSc pointwise function interface
    (void)Nf; (void)NfAux;
    (void)uOff; (void)uOff_x;
    (void)u; (void)u_t;
    (void)aOff; (void)aOff_x;
    (void)a; (void)a_x; (void)a_t;
    (void)t; (void)x;
    (void)numConstants; (void)constants;
    
    // Flux term: -k/mu * grad(P)
    const PetscReal k = 100.0e-15;  // permeability (m^2)
    const PetscReal mu = 0.001;      // viscosity (Pa·s)
    const PetscReal mobility = k / mu;
    
    for (int d = 0; d < dim; ++d) {
        f1[d] = mobility * u_x[d];
    }
}

PetscErrorCode Simulator::writeOutput(int step) {
    PetscFunctionBeginUser;
    // PetscErrorCode ierr;  // Not used in current implementation
    
    // Disable VTK output for now to avoid format incompatibility
    // TODO: Properly implement output based on DM type (DMPlex)
    if (rank == 0 && step % 10 == 0) {
        PetscPrintf(comm, "Output at step %d would be written here\n", step);
    }
    
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::writeCheckpoint(int step) {
    PetscFunctionBeginUser;
    // Implementation for checkpoint writing
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::writeSummary() {
    PetscFunctionBeginUser;
    // PetscErrorCode ierr;  // Not used in current implementation
    
    if (rank == 0) {
        std::ofstream summary(config.output_file + "_SUMMARY.txt");
        summary << "Simulation Summary\n";
        summary << "==================\n\n";
        summary << "Total timesteps: " << timestep << "\n";
        summary << "Final time: " << current_time << "\n";
        summary.close();
    }

    // Write seismic catalog if requested and events exist
    if (rank == 0 && fault_network && !seismic_catalog_file_.empty()) {
        std::ofstream cat(seismic_catalog_file_);
        cat << "time_s,magnitude,moment_Nm,slip_m,stress_drop_Pa,duration_s,area_m2,x_m,y_m,z_m,type\n";
        for (const auto& ev : fault_network->getAllEvents()) {
            cat << ev.time << ","
                << ev.magnitude << ","
                << ev.moment << ","
                << ev.slip << ","
                << ev.stress_drop << ","
                << ev.duration << ","
                << ev.area << ","
                << ev.x << ","
                << ev.y << ","
                << ev.z << ","
                << ev.type << "\n";
        }
    }
    
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::writeVisualization(int step) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    ierr = writeOutput(step); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::computePerformanceMetrics() {
    PetscFunctionBeginUser;
    // Compute and store performance metrics
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::generatePlots() {
    PetscFunctionBeginUser;
    // Generate visualization plots
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::writeStatistics() {
    PetscFunctionBeginUser;
    // Write statistics file
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::addWell(const std::string& name, WellType type) {
    PetscFunctionBeginUser;
    
    if (type == WellType::PRODUCER) {
        wells.push_back(std::make_shared<ProductionWell>(name));
    } else if (type == WellType::INJECTOR) {
        wells.push_back(std::make_shared<InjectionWell>(name));
    }
    
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::setWellControl(const std::string& name, double target) {
    PetscFunctionBeginUser;
    
    for (auto& well : wells) {
        if (well->getName() == name) {
            well->setControl(WellControlMode::RATE_CONTROL, target);
            break;
        }
    }
    
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::addFracture(FractureType type, const std::vector<double>& coords) {
    PetscFunctionBeginUser;
    
    if (type == FractureType::NATURAL) {
        auto frac = std::make_shared<NaturalFractureNetwork>();
        frac->setGeometry(coords);
        fractures.push_back(frac);
    } else if (type == FractureType::INDUCED_HYDRAULIC) {
        auto frac = std::make_shared<HydraulicFractureModel>();
        frac->setGeometry(coords);
        fractures.push_back(frac);
    }
    
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::updateFractureNetwork() {
    PetscFunctionBeginUser;
    
    for (auto& frac : fractures) {
        frac->updateGeometry(dt);
    }
    
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::enableCoupling(PhysicsType type1, PhysicsType type2) {
    PetscFunctionBeginUser;
    // Enable coupling between physics
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::adaptTimeStep(double& dt_new) {
    PetscFunctionBeginUser;
    // Adaptive time stepping logic
    dt_new = dt;
    PetscFunctionReturn(0);
}

bool Simulator::checkConvergence() {
    return true;
}

void Simulator::transformToModelCoordinates(double& x, double& y, double& z) const {
    if (coord_manager && coord_manager->isConfigured()) {
        GeoPoint pt(x, y, z);
        GeoPoint result = coord_manager->toModelCoords(pt);
        x = result.x;
        y = result.y;
        z = result.z;
    }
}

void Simulator::transformToInputCoordinates(double& x, double& y, double& z) const {
    if (coord_manager && coord_manager->isConfigured()) {
        GeoPoint pt(x, y, z);
        GeoPoint result = coord_manager->toInputCoords(pt);
        x = result.x;
        y = result.y;
        z = result.z;
    }
}

PetscErrorCode Simulator::setupBoundaryConditions() {
    PetscFunctionBeginUser;
    // Setup boundary conditions
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::interpolateMaterialProperties() {
    PetscFunctionBeginUser;
    // Interpolate material properties to quadrature points
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::solve() {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    ierr = run(); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

} // namespace FSRM
