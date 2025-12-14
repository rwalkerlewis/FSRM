#include "Simulator.hpp"
#include "Visualization.hpp"
#include "ConfigReader.hpp"
#include "ImplicitExplicitTransition.hpp"
#include "ExplosionImpactPhysics.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <limits>

namespace {

struct GmshPhysicalName {
    int dim = -1;
    int tag = -1;
    std::string name;
};

// Parse only the $PhysicalNames section from a .msh file.
// This works for both ASCII and binary Gmsh meshes because the header sections are ASCII.
static std::unordered_map<std::string, GmshPhysicalName>
parseGmshPhysicalNames(const std::string& filename) {
    std::unordered_map<std::string, GmshPhysicalName> out;

    std::ifstream file(filename, std::ios::in | std::ios::binary);
    if (!file.is_open()) return out;

    std::string line;
    while (std::getline(file, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line != "$PhysicalNames") continue;

        // Next line: number of physical names
        if (!std::getline(file, line)) break;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        int n = 0;
        {
            std::istringstream iss(line);
            iss >> n;
        }

        for (int i = 0; i < n; ++i) {
            if (!std::getline(file, line)) break;
            if (!line.empty() && line.back() == '\r') line.pop_back();

            std::istringstream iss(line);
            int dim = -1, tag = -1;
            iss >> dim >> tag;

            // The rest of the line is the quoted name; allow spaces.
            std::string rest;
            std::getline(iss, rest);
            auto q0 = rest.find('"');
            auto q1 = rest.rfind('"');
            std::string name;
            if (q0 != std::string::npos && q1 != std::string::npos && q1 > q0) {
                name = rest.substr(q0 + 1, q1 - q0 - 1);
            } else {
                // Fallback: third token as name (may still be quoted)
                std::istringstream iss2(line);
                std::string nameTok;
                iss2 >> dim >> tag >> nameTok;
                if (!nameTok.empty() && nameTok.front() == '"') nameTok.erase(nameTok.begin());
                if (!nameTok.empty() && nameTok.back() == '"') nameTok.pop_back();
                name = nameTok;
            }

            if (!name.empty()) out[name] = GmshPhysicalName{dim, tag, name};
        }

        // Consume until end marker, then stop scanning to avoid hitting binary blocks.
        while (std::getline(file, line)) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (line == "$EndPhysicalNames") break;
        }
        break;
    }

    return out;
}

static PetscErrorCode computeDMBoundingBox(MPI_Comm comm, DM dm, std::array<double, 3>& bmin, std::array<double, 3>& bmax) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

    bmin = {std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max()};
    bmax = {std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest()};

    PetscInt cdim = 0;
    ierr = DMGetCoordinateDim(dm, &cdim); CHKERRQ(ierr);
    if (cdim <= 0) PetscFunctionReturn(0);

    Vec coords = nullptr;
    ierr = DMGetCoordinatesLocal(dm, &coords); CHKERRQ(ierr);
    if (!coords) {
        ierr = DMGetCoordinates(dm, &coords); CHKERRQ(ierr);
    }
    if (!coords) PetscFunctionReturn(0);

    PetscSection csec = nullptr;
    ierr = DMGetCoordinateSection(dm, &csec); CHKERRQ(ierr);
    if (!csec) PetscFunctionReturn(0);

    PetscInt vStart = 0, vEnd = 0;
    ierr = DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd); CHKERRQ(ierr);

    const PetscScalar* a = nullptr;
    ierr = VecGetArrayRead(coords, &a); CHKERRQ(ierr);
    for (PetscInt v = vStart; v < vEnd; ++v) {
        PetscInt dof = 0, off = 0;
        ierr = PetscSectionGetDof(csec, v, &dof); CHKERRQ(ierr);
        if (dof <= 0) continue;
        ierr = PetscSectionGetOffset(csec, v, &off); CHKERRQ(ierr);
        for (PetscInt d = 0; d < std::min<PetscInt>(cdim, 3); ++d) {
            double x = static_cast<double>(a[off + d]);
            bmin[d] = std::min(bmin[d], x);
            bmax[d] = std::max(bmax[d], x);
        }
        if (cdim < 3) {
            for (PetscInt d = cdim; d < 3; ++d) {
                bmin[d] = std::min(bmin[d], 0.0);
                bmax[d] = std::max(bmax[d], 0.0);
            }
        }
    }
    ierr = VecRestoreArrayRead(coords, &a); CHKERRQ(ierr);

    // Reduce across ranks
    double gmin[3] = {bmin[0], bmin[1], bmin[2]};
    double gmax[3] = {bmax[0], bmax[1], bmax[2]};
    MPI_Allreduce(MPI_IN_PLACE, gmin, 3, MPI_DOUBLE, MPI_MIN, comm);
    MPI_Allreduce(MPI_IN_PLACE, gmax, 3, MPI_DOUBLE, MPI_MAX, comm);
    bmin = {gmin[0], gmin[1], gmin[2]};
    bmax = {gmax[0], gmax[1], gmax[2]};

    PetscFunctionReturn(0);
}

static PetscErrorCode applyGmshNameMappingsToDM(
    MPI_Comm comm,
    DM dm,
    const std::unordered_map<std::string, GmshPhysicalName>& physicalNames,
    const std::vector<std::pair<std::string, std::string>>& materialMappings,
    const std::vector<std::tuple<std::string, std::string, bool>>& faultMappings,
    const std::vector<std::string>& boundaryNames
) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

    struct NamedLabel {
        std::string name;
        DMLabel label = nullptr;
    };

    // Collect DM labels (prefer PETSc's "* Sets" labels; avoid generic labels like "depth").
    PetscInt nlabels = 0;
    ierr = DMGetNumLabels(dm, &nlabels); CHKERRQ(ierr);
    std::vector<NamedLabel> labels;
    labels.reserve(static_cast<size_t>(std::max<PetscInt>(0, nlabels)));
    for (PetscInt i = 0; i < nlabels; ++i) {
        const char* lname = nullptr;
        DMLabel lbl = nullptr;
        ierr = DMGetLabelName(dm, i, &lname); CHKERRQ(ierr);
        if (!lname) continue;
        ierr = DMGetLabel(dm, lname, &lbl); CHKERRQ(ierr);
        if (!lbl) continue;

        std::string n(lname);
        // Skip labels that are not physical-group markers.
        if (n == "depth" || n == "celltype" || n == "Material" || n == "Boundary" || n == "Fault") continue;
        labels.push_back(NamedLabel{n, lbl});
    }

    auto findPhysicalTag = [&](const std::string& group) -> int {
        auto it = physicalNames.find(group);
        if (it == physicalNames.end()) return -1;
        return it->second.tag;
    };

    // Build mapping tables keyed by Gmsh physical tag.
    std::unordered_map<int, int> physTagToMaterialId;
    std::unordered_set<int> materialPhysTags;
    {
        int nextId = 0;
        for (const auto& [group, materialSection] : materialMappings) {
            (void)materialSection; // material section is stored elsewhere; here we only label numeric IDs
            int tag = findPhysicalTag(group);
            if (tag < 0) {
                PetscPrintf(comm, "Warning: material physical group '%s' not found in $PhysicalNames\n", group.c_str());
                continue;
            }
            if (!physTagToMaterialId.count(tag)) {
                physTagToMaterialId[tag] = nextId++;
            }
            materialPhysTags.insert(tag);
        }
    }

    std::unordered_map<int, int> physTagToFaultId;
    std::unordered_set<int> faultPhysTags;
    {
        int nextId = 0;
        for (const auto& tup : faultMappings) {
            const auto& group = std::get<0>(tup);
            int tag = findPhysicalTag(group);
            if (tag < 0) {
                PetscPrintf(comm, "Warning: fault physical group '%s' not found in $PhysicalNames\n", group.c_str());
                continue;
            }
            if (!physTagToFaultId.count(tag)) {
                physTagToFaultId[tag] = nextId++;
            }
            faultPhysTags.insert(tag);
        }
    }

    std::unordered_map<int, int> physTagToBoundaryId;
    std::unordered_set<int> boundaryPhysTags;
    {
        int nextId = 0;
        for (const auto& group : boundaryNames) {
            int tag = findPhysicalTag(group);
            if (tag < 0) {
                PetscPrintf(comm, "Warning: boundary physical group '%s' not found in $PhysicalNames\n", group.c_str());
                continue;
            }
            if (!physTagToBoundaryId.count(tag)) {
                physTagToBoundaryId[tag] = nextId++;
            }
            boundaryPhysTags.insert(tag);
        }
    }

    auto matchTagForPoint = [&](PetscInt point,
                                const std::unordered_set<int>& wanted,
                                const std::vector<NamedLabel>& sources) -> int {
        if (wanted.empty()) return -1;
        for (const auto& nl : sources) {
            PetscInt v = -1;
            DMLabelGetValue(nl.label, point, &v);
            if (v >= 0 && wanted.count(static_cast<int>(v))) return static_cast<int>(v);
        }
        return -1;
    };

    // Material label on cells (height 0)
    if (!materialPhysTags.empty()) {
        DMLabel matLbl = nullptr;
        ierr = DMCreateLabel(dm, "Material"); CHKERRQ(ierr);
        ierr = DMGetLabel(dm, "Material", &matLbl); CHKERRQ(ierr);

        // Prefer PETSc's Gmsh "Cell Sets" label if present.
        std::vector<NamedLabel> cellSources;
        for (const auto& nl : labels) {
            if (nl.name == "Cell Sets") cellSources.push_back(nl);
        }
        const auto& sources = cellSources.empty() ? labels : cellSources;

        PetscInt cStart = 0, cEnd = 0;
        ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd); CHKERRQ(ierr);
        for (PetscInt c = cStart; c < cEnd; ++c) {
            int tag = matchTagForPoint(c, materialPhysTags, sources);
            if (tag < 0) continue;
            auto it = physTagToMaterialId.find(tag);
            if (it == physTagToMaterialId.end()) continue;
            ierr = DMLabelSetValue(matLbl, c, it->second); CHKERRQ(ierr);
        }
    }

    // Fault/Boundary labels on faces (height 1)
    if (!faultPhysTags.empty() || !boundaryPhysTags.empty()) {
        PetscInt fStart = 0, fEnd = 0;
        ierr = DMPlexGetHeightStratum(dm, 1, &fStart, &fEnd); CHKERRQ(ierr);

        // Prefer PETSc's Gmsh "Face Sets" label if present.
        std::vector<NamedLabel> faceSources;
        for (const auto& nl : labels) {
            if (nl.name == "Face Sets") faceSources.push_back(nl);
        }
        const auto& sources = faceSources.empty() ? labels : faceSources;

        DMLabel faultLbl = nullptr;
        if (!faultPhysTags.empty()) {
            ierr = DMCreateLabel(dm, "Fault"); CHKERRQ(ierr);
            ierr = DMGetLabel(dm, "Fault", &faultLbl); CHKERRQ(ierr);
        }
        DMLabel bndLbl = nullptr;
        if (!boundaryPhysTags.empty()) {
            ierr = DMCreateLabel(dm, "Boundary"); CHKERRQ(ierr);
            ierr = DMGetLabel(dm, "Boundary", &bndLbl); CHKERRQ(ierr);
        }

        for (PetscInt f = fStart; f < fEnd; ++f) {
            if (faultLbl) {
                int tag = matchTagForPoint(f, faultPhysTags, sources);
                if (tag >= 0) {
                    auto it = physTagToFaultId.find(tag);
                    if (it != physTagToFaultId.end()) {
                        ierr = DMLabelSetValue(faultLbl, f, it->second); CHKERRQ(ierr);
                    }
                }
            }
            if (bndLbl) {
                int tag = matchTagForPoint(f, boundaryPhysTags, sources);
                if (tag >= 0) {
                    auto it = physTagToBoundaryId.find(tag);
                    if (it != physTagToBoundaryId.end()) {
                        ierr = DMLabelSetValue(bndLbl, f, it->second); CHKERRQ(ierr);
                    }
                }
            }
        }
    }

    PetscFunctionReturn(0);
}

} // namespace

namespace FSRM {

// =============================================================================
// ExplosionCoupling (PIMPL)
// =============================================================================

struct Simulator::ExplosionCoupling {
    // Source timing
    double t0 = 0.0;          // detonation time (s)

    // Source location
    double sx = 0.0, sy = 0.0, sz = 0.0;

    // Medium properties for spherical cavity model
    double rho = 2700.0;
    double vp = 5500.0;
    double vs = 3200.0;

    // Cavity model (spherically symmetric proxy)
    SphericalCavitySource cavity;

    // Configure from basic nuclear underground parameters
    void configureUndergroundNuclear(double yield_kt, double depth_m,
                                     double x, double y, double z,
                                     double rho_in, double vp_in, double vs_in,
                                     double rise_time_s, double overpressure_pa,
                                     double detonation_time_s) {
        t0 = detonation_time_s;
        sx = x; sy = y; sz = z;
        rho = rho_in; vp = vp_in; vs = vs_in;

        // Compute elastic moduli from wave speeds (isotropic)
        double G = rho * vs * vs;
        double K = rho * (vp * vp - 4.0/3.0 * vs * vs);
        if (K < 1e8) K = 1e8;

        // Empirical cavity radius (from NuclearSourceParameters scaling)
        NuclearSourceParameters params;
        params.yield_kt = yield_kt;
        params.depth_of_burial = depth_m;
        double Rc = params.cavity_radius(rho);

        cavity.setCavityRadius(Rc);
        cavity.setMediumProperties(rho, K, G);
        cavity.setOverpressure(overpressure_pa);
        cavity.setRiseTime(rise_time_s);
    }

    // Compute spherically symmetric stress tensor at a point due to cavity source.
    // Returns full tensor components in global coordinates (x,y,z).
    void stressTensorAt(double x, double y, double z, double time,
                        double& sxx, double& syy, double& szz,
                        double& sxy, double& sxz, double& syz) const {
        sxx = syy = szz = sxy = sxz = syz = 0.0;

        double t = time - t0;
        if (t < 0.0) return;

        // Radial distance from source
        double dx = x - sx;
        double dy = y - sy;
        double dz = z - sz;
        double r = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (r < 1.0) r = 1.0;

        double n1 = dx / r;
        double n2 = dy / r;
        double n3 = dz / r;

        // Spherical cavity stresses (radial/tangential)
        double sigma_rr = 0.0, sigma_tt = 0.0;
        cavity.stress(r, t, sigma_rr, sigma_tt);

        // Construct tensor: σ = σ_tt I + (σ_rr - σ_tt) n⊗n
        double d = sigma_rr - sigma_tt;
        sxx = sigma_tt + d * n1 * n1;
        syy = sigma_tt + d * n2 * n2;
        szz = sigma_tt + d * n3 * n3;
        sxy = d * n1 * n2;
        sxz = d * n1 * n3;
        syz = d * n2 * n3;
    }
};

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

    // Explosion PDE solve mode selection (explicit knob)
    if (!config.explosion_solve_mode.empty()) {
        if (config.explosion_solve_mode != "PROXY" &&
            config.explosion_solve_mode != "COUPLED_ANALYTIC" &&
            config.explosion_solve_mode != "FULL_PDE") {
            SETERRQ(comm, PETSC_ERR_ARG_WRONG,
                    ("Invalid SIMULATION.explosion_solve_mode='" + config.explosion_solve_mode +
                     "'. Valid: PROXY, COUPLED_ANALYTIC, FULL_PDE").c_str());
        }
        if (config.explosion_solve_mode == "FULL_PDE") {
            SETERRQ(comm, PETSC_ERR_SUP,
                    "SIMULATION.explosion_solve_mode=FULL_PDE requested, but a full coupled blast+shock+wave PDE solver is not implemented in this build. "
                    "Use SIMULATION.explosion_solve_mode=COUPLED_ANALYTIC (integrated spherical-cavity stress coupling) or PROXY.");
        }
    }
    
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

    // Explosion coupling (integrated into fault updates and IMEX triggering)
    if (reader.hasSection("EXPLOSION_SOURCE")) {
        std::string ex_type = reader.getString("EXPLOSION_SOURCE", "type", "");
        // Currently only support underground nuclear coupling via spherical cavity proxy.
        if (!ex_type.empty() && (ex_type == "NUCLEAR_UNDERGROUND" || ex_type == "UNDERGROUND_CONTAINED")) {
            double yield_kt = reader.getDouble("EXPLOSION_SOURCE", "yield_kt", 1.0);
            double depth_m = reader.getDouble("EXPLOSION_SOURCE", "depth_of_burial", 1000.0);
            double x = reader.getDouble("EXPLOSION_SOURCE", "location_x",
                                        reader.getDouble("EXPLOSION_SOURCE", "source_x", 0.0));
            double y = reader.getDouble("EXPLOSION_SOURCE", "location_y",
                                        reader.getDouble("EXPLOSION_SOURCE", "source_y", 0.0));
            double z = reader.getDouble("EXPLOSION_SOURCE", "location_z",
                                        reader.getDouble("EXPLOSION_SOURCE", "source_z", -depth_m));
            double t0 = reader.getDouble("EXPLOSION_SOURCE", "onset_time",
                                         reader.getDouble("EXPLOSION_SOURCE", "time0", 20.0));

            // Medium properties: use first material by default
            double rho = 2700.0, vp = 5500.0, vs = 3200.0;
            if (!material_props.empty()) {
                rho = material_props[0].density;
                vp = material_props[0].p_wave_velocity;
                vs = material_props[0].s_wave_velocity;
            }

            // Rise time: allow override, else use a yield-scaled heuristic (seconds)
            double rise = reader.getDouble("EXPLOSION_SOURCE", "rise_time", 0.05);
            // Overpressure at cavity wall: allow override, else use a rough scaling
            double P0 = reader.getDouble("EXPLOSION_SOURCE", "cavity_overpressure", 1.0e11 * std::pow(std::max(0.1, yield_kt), 0.25));

            explosion_ = std::make_unique<ExplosionCoupling>();
            explosion_->configureUndergroundNuclear(yield_kt, depth_m, x, y, z, rho, vp, vs, rise, P0, t0);
        }
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
    
    // Use PETSc's built-in Gmsh reader (supports ASCII and binary .msh).
    ierr = DMPlexCreateGmshFromFile(comm, grid_config.mesh_file.c_str(), PETSC_TRUE, &dm); CHKERRQ(ierr);

    // Configure domain mappings (materials / faults / boundaries) if provided.
    // We map by *physical group name* using $PhysicalNames (ASCII even for binary meshes),
    // and apply labels onto the PETSc DM based on whatever labels the reader created.
    const bool have_mappings =
        !grid_config.gmsh_material_domains.empty() ||
        !grid_config.gmsh_fault_surfaces.empty() ||
        !grid_config.gmsh_boundaries.empty();

    if (have_mappings) {
        auto physicalNames = parseGmshPhysicalNames(grid_config.mesh_file);
        if (physicalNames.empty()) {
            PetscPrintf(comm,
                        "Warning: no $PhysicalNames section found in '%s'. "
                        "Name-based mappings (materials/faults/boundaries) may not be applied.\n",
                        grid_config.mesh_file.c_str());
        } else {
            ierr = applyGmshNameMappingsToDM(
                comm, dm, physicalNames,
                grid_config.gmsh_material_domains,
                grid_config.gmsh_fault_surfaces,
                grid_config.gmsh_boundaries
            ); CHKERRQ(ierr);
        }
    }

    // Update grid config with actual mesh dimensions from DM coordinates
    std::array<double, 3> bmin, bmax;
    ierr = computeDMBoundingBox(comm, dm, bmin, bmax); CHKERRQ(ierr);
    grid_config.Lx = bmax[0] - bmin[0];
    grid_config.Ly = bmax[1] - bmin[1];
    grid_config.Lz = bmax[2] - bmin[2];
    
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
    PetscReal lower[3] = {grid_config.origin_x, grid_config.origin_y, grid_config.origin_z};
    PetscReal upper[3] = {grid_config.origin_x + grid_config.Lx,
                          grid_config.origin_y + grid_config.Ly,
                          grid_config.origin_z + grid_config.Lz};
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
        case FluidModelType::NONE:
            // No fluid fields
            break;
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
        case FluidModelType::NONE:
            break;
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

    // Update fault network with a simplified background stress + integrated explosion coupling.
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

        const double shear_modulus = 30.0e9;  // Pa (order-of-magnitude)

        // Per-fault update so explosion-induced stresses can vary spatially.
        for (size_t i = 0; i < sim->fault_network->numFaults(); ++i) {
            auto* fault = sim->fault_network->getFault(i);
            if (!fault) continue;

            double ex_sxx = 0.0, ex_syy = 0.0, ex_szz = 0.0, ex_sxy = 0.0, ex_sxz = 0.0, ex_syz = 0.0;
            if (sim->explosion_) {
                const auto& geom = fault->getGeometry();
                sim->explosion_->stressTensorAt(geom.x, geom.y, geom.z, static_cast<double>(t),
                                                ex_sxx, ex_syy, ex_szz, ex_sxy, ex_sxz, ex_syz);
            } else if (sim->nuclear_trigger_enabled_) {
                // Backward-compatible synthetic trigger if no integrated explosion is configured.
                const double tt = static_cast<double>(t);
                const double t0 = sim->nuclear_trigger_time0_;
                if (tt >= t0) {
                    const double tr = std::max(1e-6, sim->nuclear_trigger_rise_);
                    const double td = std::max(1e-6, sim->nuclear_trigger_decay_);
                    double rise = 1.0 - std::exp(-(tt - t0) / tr);
                    double decay = std::exp(-(tt - t0) / td);
                    double amp = rise * decay;
                    ex_sxz += amp * sim->nuclear_trigger_dtau_;
                }
            }

            // Combine background + explosion-induced stresses
            double Sxx = sxx + ex_sxx;
            double Syy = syy + ex_syy;
            double Szz = szz + ex_szz;
            double Sxy = ex_sxy;
            double Sxz = ex_sxz;
            double Syz = ex_syz;

            double pore = 0.0;
            FaultStressState state = fault->computeStressState(Sxx, Syy, Szz, Sxy, Sxz, Syz, pore);
            FaultStressState new_state = fault->updateSlip(state, shear_modulus, dt);

            // Record event
            if (new_state.slip - state.slip > 1e-6 && new_state.is_seismic) {
                SeismicEvent event = fault->computeEvent(state, new_state, shear_modulus, static_cast<double>(t));
                fault->addEvent(event);
            }
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
