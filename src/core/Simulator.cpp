#include "core/Simulator.hpp"
#include "core/PhysicsModuleRegistry.hpp"
#include "numerics/PetscFEFluidFlow.hpp"
#include "numerics/PetscFEElasticity.hpp"
#include "numerics/PetscFEElasticityAux.hpp"
#include "numerics/PetscFEElastoplasticity.hpp"
#include "numerics/PetscFEViscoelastic.hpp"
#include "numerics/PetscFEPoroelasticity.hpp"
#include "numerics/PetscFEHydrofrac.hpp"
#include "numerics/PetscFEThermal.hpp"
#include "numerics/AbsorbingBC.hpp"
#include "io/Visualization.hpp"
#include "core/ConfigReader.hpp"
#include "numerics/ImplicitExplicitTransition.hpp"
#include "domain/explosion/ExplosionImpactPhysics.hpp"
#include "domain/explosion/NearFieldExplosion.hpp"
#include "domain/seismic/SeismometerNetwork.hpp"
#include "io/VelocityModelReader.hpp"

// GPU kernel support — conditionally include GPU headers when built with CUDA
#ifdef USE_CUDA
#include "gpu/GPUManager.hpp"
#include "physics/PhysicsKernel_GPU.hpp"
#include "physics/PhysicsKernel_GPU_Extended.hpp"
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <cctype>
#include <cmath>
#include <algorithm>
#include <vector>
#include <filesystem>
#include <unordered_map>
#include <unordered_set>
#include <limits>
#include <sys/stat.h>
#include <petscviewerhdf5.h>
#include <petscfe.h>

namespace {

PetscErrorCode ensureDirectoryExists(MPI_Comm comm, int rank, const std::string& path) {
    PetscFunctionBeginUser;

    if (rank == 0 && !path.empty()) {
        std::error_code ec;
        std::filesystem::create_directories(path, ec);
        PetscCheck(!ec, comm, PETSC_ERR_FILE_OPEN,
                   "Failed to create output directory %s: %s",
                   path.c_str(), ec.message().c_str());
    }
    MPI_Barrier(comm);

    PetscFunctionReturn(PETSC_SUCCESS);
}

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

    // Stored source parameters for Mueller-Murphy moment rate
    double yield_kt = 0.0;
    double depth_of_burial = 0.0;

    // Cavity model (spherically symmetric proxy)
    SphericalCavitySource cavity;

    // COUPLED_ANALYTIC: 1D solver results
    bool use_nearfield_coupling = false;
    RDPSeismicSource rdp_source;
    NearFieldExplosionSolver nf_solver;
    double nf_cavity_radius = 0.0;
    double nf_crushed_radius = 0.0;
    double nf_fractured_radius = 0.0;

    // Configure from basic nuclear underground parameters
    void configureUndergroundNuclear(double yield_kt, double depth_m,
                                     double x, double y, double z,
                                     double rho_in, double vp_in, double vs_in,
                                     double rise_time_s, double overpressure_pa,
                                     double detonation_time_s) {
        t0 = detonation_time_s;
        sx = x; sy = y; sz = z;
        rho = rho_in; vp = vp_in; vs = vs_in;
        this->yield_kt = yield_kt;
        depth_of_burial = depth_m;

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
    if (auxVec_) VecDestroy(&auxVec_);
    if (auxDM_) DMDestroy(&auxDM_);
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

    // Parse injection source (optional)
    if (reader.hasSection("INJECTION")) {
        injection_enabled_ = reader.getBool("INJECTION", "enabled", true);
        injection_x_ = reader.getDouble("INJECTION", "x", 0.0);
        injection_y_ = reader.getDouble("INJECTION", "y", 0.0);
        injection_z_ = reader.getDouble("INJECTION", "z", 0.0);
        injection_rate_ = reader.getDouble("INJECTION", "rate", 1e-4);
        injection_start_ = reader.getDouble("INJECTION", "start_time", 0.0);
        injection_end_ = reader.getDouble("INJECTION", "end_time", 1e30);
        if (rank == 0 && injection_enabled_) {
            PetscPrintf(comm, "Injection source: (%.1f, %.1f, %.1f), rate=%.2e m^3/s\n",
                        injection_x_, injection_y_, injection_z_, injection_rate_);
        }
    }

    // Parse hydraulic fracture model (optional)
    if (reader.hasSection("HYDRAULIC_FRACTURE")) {
        hydrofrac_ = std::make_unique<HydraulicFractureModel>();
        hydrofrac_->setFormationProperties(
            reader.getDouble("HYDRAULIC_FRACTURE", "youngs_modulus", 10e9),
            reader.getDouble("HYDRAULIC_FRACTURE", "poissons_ratio", 0.25));
        hydrofrac_->setFractureToughness(
            reader.getDouble("HYDRAULIC_FRACTURE", "toughness", 1.5e6));
        hydrofrac_->setStressState(
            reader.getDouble("HYDRAULIC_FRACTURE", "min_horizontal_stress", 20e6),
            reader.getDouble("HYDRAULIC_FRACTURE", "max_horizontal_stress", 25e6),
            reader.getDouble("HYDRAULIC_FRACTURE", "vertical_stress", 30e6));
        hydrofrac_->setHeight(
            reader.getDouble("HYDRAULIC_FRACTURE", "height", 50.0));
        hydrofrac_->setFluidProperties(
            reader.getDouble("HYDRAULIC_FRACTURE", "fluid_density", 1000.0),
            reader.getDouble("HYDRAULIC_FRACTURE", "fluid_viscosity", 0.001));
        std::string model = reader.getString("HYDRAULIC_FRACTURE", "model", "PKN");
        if (model == "KGD") hydrofrac_->setModel("KGD");
        else if (model == "P3D") hydrofrac_->setModel("P3D");
        else hydrofrac_->setModel("PKN");

        hydrofrac_fem_pressurized_mode_ =
            reader.getBool("HYDRAULIC_FRACTURE", "enable_fem_pressurized", false);
        hydrofrac_uniform_pressure_pa_ =
            reader.getDouble("HYDRAULIC_FRACTURE", "uniform_pressure_pa", 0.0);

        if (rank == 0) {
            PetscPrintf(comm, "Hydraulic fracture model: %s, K_Ic=%.2e\n",
                        model.c_str(), hydrofrac_->getFractureToughness());
            if (hydrofrac_fem_pressurized_mode_) {
                PetscPrintf(comm, "Hydrofrac FEM pressurized mode enabled, p_f=%.3e Pa\n",
                            hydrofrac_uniform_pressure_pa_);
            }
        }
    }

    // Parse fracture plane for cohesive fracture propagation (optional)
    if (reader.hasSection("FRACTURE_PLANE")) {
        fracture_plane_enabled_ = reader.getBool("FRACTURE_PLANE", "enabled", true);
        if (fracture_plane_enabled_) {
            double strike_deg = reader.getDouble("FRACTURE_PLANE", "strike", 0.0);
            double dip_deg = reader.getDouble("FRACTURE_PLANE", "dip", 90.0);
            fracture_plane_strike_ = strike_deg * M_PI / 180.0;
            fracture_plane_dip_ = dip_deg * M_PI / 180.0;
            fracture_plane_center_[0] = reader.getDouble("FRACTURE_PLANE", "center_x",
                                                          grid_config.Lx / 2.0);
            fracture_plane_center_[1] = reader.getDouble("FRACTURE_PLANE", "center_y",
                                                          grid_config.Ly / 2.0);
            fracture_plane_center_[2] = reader.getDouble("FRACTURE_PLANE", "center_z",
                                                          grid_config.Lz / 2.0);
            fracture_plane_length_ = reader.getDouble("FRACTURE_PLANE", "length", 200.0);
            fracture_plane_width_ = reader.getDouble("FRACTURE_PLANE", "width", 100.0);
            fracture_plane_tensile_strength_ = reader.getDouble("FRACTURE_PLANE",
                "tensile_strength", 5.0e6);
            // Enable faults so setupFaultNetwork() runs
            config.enable_faults = true;
            if (rank == 0) {
                PetscPrintf(comm, "Fracture plane: strike=%.1f deg, dip=%.1f deg, "
                    "T_s=%.2e Pa, L=%.1f m, W=%.1f m\n",
                    strike_deg, dip_deg, fracture_plane_tensile_strength_,
                    fracture_plane_length_, fracture_plane_width_);
            }
        }
    }

    // Parse fault mode and prescribed slip (optional)
    if (reader.hasSection("FAULT")) {
        fault_mode_ = reader.getString("FAULT", "mode", "locked");
        fault_slip_strike_ = reader.getDouble("FAULT", "slip_strike", 0.0);
        fault_slip_dip_ = reader.getDouble("FAULT", "slip_dip", 0.0);
        fault_slip_opening_ = reader.getDouble("FAULT", "slip_opening", 0.0);
        fault_friction_coefficient_ = reader.getDouble("FAULT", "friction_coefficient", 0.6);
        // Read slip-weakening friction model parameters
        if (reader.hasKey("FAULT", "friction_model")) {
            std::string fm = reader.getString("FAULT", "friction_model", "constant");
            if (fm == "slip_weakening") {
                use_slip_weakening_ = true;
                mu_static_ = reader.getDouble("FAULT", "static_friction", 0.677);
                mu_dynamic_ = reader.getDouble("FAULT", "dynamic_friction", 0.525);
                critical_slip_distance_ = reader.getDouble("FAULT", "critical_slip_distance", 0.40);
                if (rank == 0) {
                    PetscPrintf(comm,
                        "Slip-weakening friction: mu_s=%.3f, mu_d=%.3f, Dc=%.4f m\n",
                        mu_static_, mu_dynamic_, critical_slip_distance_);
                }
            }
        }
        // Read optional fault geometry (center, orientation, dimensions)
        if (reader.hasKey("FAULT", "center_x") || reader.hasKey("FAULT", "center_y") ||
            reader.hasKey("FAULT", "center_z") || reader.hasKey("FAULT", "strike") ||
            reader.hasKey("FAULT", "dip") || reader.hasKey("FAULT", "length") ||
            reader.hasKey("FAULT", "width")) {
            fault_geometry_from_config_ = true;
            double strike_deg = reader.getDouble("FAULT", "strike", 0.0);
            double dip_deg = reader.getDouble("FAULT", "dip", 90.0);
            fault_strike_ = strike_deg * M_PI / 180.0;
            fault_dip_ = dip_deg * M_PI / 180.0;
            // Use -1 sentinel to indicate "read from config" vs "use default"
            fault_center_[0] = reader.getDouble("FAULT", "center_x", -1.0);
            fault_center_[1] = reader.getDouble("FAULT", "center_y", -1.0);
            fault_center_[2] = reader.getDouble("FAULT", "center_z", -1.0);
            fault_length_ = reader.getDouble("FAULT", "length", -1.0);
            fault_width_ = reader.getDouble("FAULT", "width", -1.0);
        }
        // Time-dependent prescribed slip parameters
        fault_slip_onset_time_ = reader.getDouble("FAULT", "slip_onset_time", 0.0);
        fault_slip_rise_time_ = reader.getDouble("FAULT", "slip_rise_time", 0.0);
        // Initial fault stress (for SCEC TPV-type benchmarks)
        if (reader.hasKey("FAULT", "initial_normal_stress")) {
            fault_has_initial_stress_ = true;
            fault_initial_normal_stress_ = reader.getDouble("FAULT", "initial_normal_stress", -120.0e6);
            fault_initial_shear_stress_ = reader.getDouble("FAULT", "initial_shear_stress", 70.0e6);
            fault_nucleation_center_x_ = reader.getDouble("FAULT", "nucleation_center_x", 0.0);
            fault_nucleation_center_z_ = reader.getDouble("FAULT", "nucleation_center_z", -7500.0);
            fault_nucleation_radius_ = reader.getDouble("FAULT", "nucleation_radius", 1500.0);
            fault_nucleation_shear_stress_ = reader.getDouble("FAULT", "nucleation_shear_stress", 81.6e6);
            if (rank == 0) {
                PetscPrintf(comm,
                    "Initial fault stress: sigma_n=%.3e Pa, tau=%.3e Pa\n"
                    "  Nucleation patch: center=(%.1f, %.1f), radius=%.1f m, tau=%.3e Pa\n",
                    fault_initial_normal_stress_, fault_initial_shear_stress_,
                    fault_nucleation_center_x_, fault_nucleation_center_z_,
                    fault_nucleation_radius_, fault_nucleation_shear_stress_);
            }
        }
        config.enable_faults = true;
    }

    // Parse elastoplasticity configuration (optional)
    if (reader.hasSection("PLASTICITY")) {
        config.enable_elastoplasticity = reader.getBool("PLASTICITY", "enabled", false);
        config.ep_cohesion = reader.getDouble("PLASTICITY", "cohesion", 1.0e6);
        double friction_deg = reader.getDouble("PLASTICITY", "friction_angle", 30.0);
        double dilation_deg = reader.getDouble("PLASTICITY", "dilation_angle", 0.0);
        config.ep_friction_angle = friction_deg * M_PI / 180.0;
        config.ep_dilation_angle = dilation_deg * M_PI / 180.0;
        config.ep_hardening_modulus = reader.getDouble("PLASTICITY", "hardening_modulus", 0.0);
    }

    // Parse viscoelastic attenuation configuration (optional)
    if (reader.hasSection("VISCOELASTIC")) {
        config.enable_viscoelastic = reader.getBool("VISCOELASTIC", "enabled", false);
        config.visco_num_mechanisms = static_cast<int>(reader.getDouble("VISCOELASTIC", "num_mechanisms", 3.0));
        config.visco_q_p = reader.getDouble("VISCOELASTIC", "q_p", 200.0);
        config.visco_q_s = reader.getDouble("VISCOELASTIC", "q_s", 100.0);
        config.visco_f_min = reader.getDouble("VISCOELASTIC", "f_min", 0.1);
        config.visco_f_max = reader.getDouble("VISCOELASTIC", "f_max", 10.0);
        if (config.visco_num_mechanisms < 1) config.visco_num_mechanisms = 1;
        if (config.visco_num_mechanisms > 5) config.visco_num_mechanisms = 5;
        if (rank == 0 && config.enable_viscoelastic) {
            PetscPrintf(comm, "Viscoelastic attenuation: N=%d, Qp=%.1f, Qs=%.1f, f=[%.2f, %.2f] Hz\n",
                config.visco_num_mechanisms, config.visco_q_p, config.visco_q_s,
                config.visco_f_min, config.visco_f_max);
        }
    }

    // Parse heterogeneous material configuration (optional)
    if (reader.hasSection("MATERIAL")) {
        config.use_heterogeneous_material = reader.getBool("MATERIAL", "heterogeneous", false);
        config.material_assignment = reader.getString("MATERIAL", "assignment", "depth");
        std::transform(config.material_assignment.begin(), config.material_assignment.end(),
                       config.material_assignment.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        config.gravity = reader.getDouble("MATERIAL", "gravity", 0.0);
        config.velocity_model_file = reader.getString("MATERIAL", "velocity_model_file", "");
        // Velocity model assignment implies heterogeneous material
        if (config.material_assignment == "velocity_model") {
            config.use_heterogeneous_material = true;
        }
    }
    // Also accept gravity from [SIMULATION] section
    if (reader.hasSection("SIMULATION")) {
        config.enable_gravity = reader.getBool("SIMULATION", "enable_gravity", config.gravity > 0.0);
        if (config.enable_gravity && config.gravity == 0.0) {
            config.gravity = reader.getDouble("SIMULATION", "gravity", 9.81);
        }
        config.K0 = reader.getDouble("SIMULATION", "K0", 0.5);
    }
    if (config.use_heterogeneous_material) {
        // Look for LAYER_N sections
        config.material_layers.clear();
        for (int i = 1; i <= 100; ++i) {
            std::string section_name = "LAYER_" + std::to_string(i);
            if (!reader.hasSection(section_name)) break;
            MaterialLayer layer;
            layer.z_top    = reader.getDouble(section_name, "z_top", 0.0);
            layer.z_bottom = reader.getDouble(section_name, "z_bottom", 0.0);
            layer.lambda   = reader.getDouble(section_name, "lambda", 30.0e9);
            layer.mu       = reader.getDouble(section_name, "mu", 25.0e9);
            layer.rho      = reader.getDouble(section_name, "rho", 2650.0);
            config.material_layers.push_back(layer);
        }
        if (rank == 0) {
            PetscPrintf(comm, "Heterogeneous material: %d layers defined\n",
                        (int)config.material_layers.size());
        }
    }

    config.material_regions.clear();
    for (int i = 1; i <= 100; ++i) {
        std::string section_name = "MATERIAL_REGION_" + std::to_string(i);
        if (!reader.hasSection(section_name)) break;

        MaterialRegion region;
        region.gmsh_label_name = reader.getString(section_name, "gmsh_label", "");
        if (region.gmsh_label_name.empty()) continue;

        const double youngs_modulus = reader.getDouble(section_name, "youngs_modulus", 30.0e9);
        const double poissons_ratio = reader.getDouble(section_name, "poissons_ratio", 0.25);
        region.lambda = youngs_modulus * poissons_ratio /
                        ((1.0 + poissons_ratio) * (1.0 - 2.0 * poissons_ratio));
        region.mu = youngs_modulus / (2.0 * (1.0 + poissons_ratio));
        region.rho = reader.getDouble(section_name, "density", 2650.0);
        config.material_regions.push_back(region);
    }
    if (!config.material_regions.empty()) {
        config.use_heterogeneous_material = true;
        if (config.material_assignment.empty() || config.material_assignment == "depth") {
            config.material_assignment = "gmsh_label";
        }
        if (rank == 0) {
            PetscPrintf(comm, "Material-region mapping: %d gmsh regions defined\n",
                        (int)config.material_regions.size());
        }
    }

    // Parse absorbing boundary conditions (optional)
    if (reader.hasSection("ABSORBING_BC")) {
        config.absorbing_bc_enabled = reader.getBool("ABSORBING_BC", "enabled", false);
        config.absorbing_bc_x_min = reader.getBool("ABSORBING_BC", "x_min", true);
        config.absorbing_bc_x_max = reader.getBool("ABSORBING_BC", "x_max", true);
        config.absorbing_bc_y_min = reader.getBool("ABSORBING_BC", "y_min", true);
        config.absorbing_bc_y_max = reader.getBool("ABSORBING_BC", "y_max", true);
        config.absorbing_bc_z_min = reader.getBool("ABSORBING_BC", "z_min", true);
        config.absorbing_bc_z_max = reader.getBool("ABSORBING_BC", "z_max", false);
        if (config.absorbing_bc_enabled && rank == 0) {
            PetscPrintf(comm, "Absorbing BCs: x=[%d,%d] y=[%d,%d] z=[%d,%d]\n",
                (int)config.absorbing_bc_x_min, (int)config.absorbing_bc_x_max,
                (int)config.absorbing_bc_y_min, (int)config.absorbing_bc_y_max,
                (int)config.absorbing_bc_z_min, (int)config.absorbing_bc_z_max);
        }
    }

    // Parse thermal configuration (optional)
    if (reader.hasSection("THERMAL")) {
        config.enable_thermal = true;
        // Override thermal properties from [THERMAL] section if specified
        if (reader.hasKey("THERMAL", "thermal_conductivity") && !material_props.empty()) {
            material_props[0].thermal_conductivity =
                reader.getDouble("THERMAL", "thermal_conductivity",
                                 material_props[0].thermal_conductivity);
        }
        if (reader.hasKey("THERMAL", "specific_heat") && !material_props.empty()) {
            material_props[0].heat_capacity =
                reader.getDouble("THERMAL", "specific_heat",
                                 material_props[0].heat_capacity);
        }
        if (reader.hasKey("THERMAL", "thermal_expansion") && !material_props.empty()) {
            material_props[0].thermal_expansion =
                reader.getDouble("THERMAL", "thermal_expansion",
                                 material_props[0].thermal_expansion);
        }
        config.reference_temperature =
            reader.getDouble("THERMAL", "reference_temperature", 293.0);
        if (rank == 0) {
            PetscPrintf(comm, "Thermal coupling enabled: T_ref = %.1f K\n",
                        config.reference_temperature);
        }
    }

    // Parse boundary condition configuration (optional)
    if (reader.hasSection("BOUNDARY_CONDITIONS")) {
        config.bottom_bc = reader.getString("BOUNDARY_CONDITIONS", "bottom", "fixed");
        config.side_bc = reader.getString("BOUNDARY_CONDITIONS", "sides", "free");
        config.top_bc = reader.getString("BOUNDARY_CONDITIONS", "top", "compression");

        // Parse per-face boundary conditions (overrides bottom/sides/top)
        // Face ordering: 0=x_min, 1=x_max, 2=y_min, 3=y_max, 4=z_min, 5=z_max
        static const char *face_keys[6] = {
            "x_min", "x_max", "y_min", "y_max", "z_min", "z_max"
        };
        for (int fi = 0; fi < 6; ++fi) {
            std::string ftype = reader.getString("BOUNDARY_CONDITIONS", face_keys[fi], "");
            if (!ftype.empty()) {
                config.face_bc[fi].configured = true;
                config.face_bc[fi].type = ftype;
                if (ftype == "traction") {
                    // Read traction vector: x_min_traction_x, x_min_traction_y, x_min_traction_z
                    // or compact form: x_min_traction = tx, ty, tz
                    std::string tkey = std::string(face_keys[fi]) + "_traction";
                    if (reader.hasKey("BOUNDARY_CONDITIONS", tkey)) {
                        auto tv = reader.getDoubleArray("BOUNDARY_CONDITIONS", tkey);
                        for (int d = 0; d < 3 && d < static_cast<int>(tv.size()); ++d) {
                            config.face_bc[fi].traction[d] = tv[static_cast<std::size_t>(d)];
                        }
                    } else {
                        config.face_bc[fi].traction[0] = reader.getDouble("BOUNDARY_CONDITIONS",
                            std::string(face_keys[fi]) + "_traction_x", 0.0);
                        config.face_bc[fi].traction[1] = reader.getDouble("BOUNDARY_CONDITIONS",
                            std::string(face_keys[fi]) + "_traction_y", 0.0);
                        config.face_bc[fi].traction[2] = reader.getDouble("BOUNDARY_CONDITIONS",
                            std::string(face_keys[fi]) + "_traction_z", 0.0);
                    }
                } else if (ftype == "dirichlet") {
                    // Read constrained components and values
                    std::string ckey = std::string(face_keys[fi]) + "_components";
                    std::string vkey = std::string(face_keys[fi]) + "_values";
                    if (reader.hasKey("BOUNDARY_CONDITIONS", ckey)) {
                        auto cv = reader.getDoubleArray("BOUNDARY_CONDITIONS", ckey);
                        for (auto c : cv) {
                            config.face_bc[fi].components.push_back(static_cast<int>(c));
                        }
                    } else {
                        config.face_bc[fi].components = {0, 1, 2};
                    }
                    if (reader.hasKey("BOUNDARY_CONDITIONS", vkey)) {
                        config.face_bc[fi].values = reader.getDoubleArray("BOUNDARY_CONDITIONS", vkey);
                    }
                    // Fill to match component count
                    while (config.face_bc[fi].values.size() < config.face_bc[fi].components.size()) {
                        config.face_bc[fi].values.push_back(0.0);
                    }
                } else if (ftype == "dirichlet_pressure") {
                    // Read pressure value for Dirichlet pressure BC on fluid field
                    std::string pkey = std::string(face_keys[fi]) + "_pressure";
                    config.face_bc[fi].pressure = reader.getDouble("BOUNDARY_CONDITIONS", pkey, 0.0);
                }
            }
        }

        if (rank == 0) {
            PetscPrintf(comm, "Boundary conditions: bottom=%s, sides=%s, top=%s\n",
                        config.bottom_bc.c_str(), config.side_bc.c_str(), config.top_bc.c_str());
            for (int fi = 0; fi < 6; ++fi) {
                if (config.face_bc[fi].configured) {
                    PetscPrintf(comm, "  %s = %s", face_keys[fi], config.face_bc[fi].type.c_str());
                    if (config.face_bc[fi].type == "traction") {
                        PetscPrintf(comm, " (%.2e, %.2e, %.2e)",
                                    config.face_bc[fi].traction[0],
                                    config.face_bc[fi].traction[1],
                                    config.face_bc[fi].traction[2]);
                    }
                    PetscPrintf(comm, "\n");
                }
            }
        }
    }

    // Parse initial conditions (optional)
    if (reader.hasSection("INITIAL_CONDITIONS")) {
        config.initial_condition_type = reader.getString("INITIAL_CONDITIONS", "type", "zero");
        config.initial_condition_path = reader.getString("INITIAL_CONDITIONS", "path", "");
        if (rank == 0) {
            PetscPrintf(comm, "Initial conditions: type=%s", config.initial_condition_type.c_str());
            if (!config.initial_condition_path.empty()) {
                PetscPrintf(comm, ", path=%s", config.initial_condition_path.c_str());
            }
            PetscPrintf(comm, "\n");
        }
    }

    // Parse solution save path and output directory (optional)
    if (reader.hasSection("OUTPUT")) {
        config.save_solution_path = reader.getString("OUTPUT", "save_solution", "");
        output_directory_ = reader.getString("OUTPUT", "output_directory", "output");
        std::string fmt = reader.getString("OUTPUT", "format", "");
        if (!fmt.empty()) {
            config.output_format = fmt;
        }
        config.output_stress = reader.getBool("OUTPUT", "output_stress", false);
        config.output_cfs = reader.getBool("OUTPUT", "output_cfs", false);
        config.cfs_receiver_strike = reader.getDouble("OUTPUT", "cfs_receiver_strike", 0.0);
        config.cfs_receiver_dip = reader.getDouble("OUTPUT", "cfs_receiver_dip", 90.0);
        config.cfs_friction = reader.getDouble("OUTPUT", "cfs_friction", 0.4);
    }

    // Parse seismometers (optional)
    seismo_specs_.clear();
    seismo_out_cfg_ = SeismometerOutputConfig{};
    if (reader.hasSection("SEISMOMETERS")) {
        auto strtoupper_local = [](std::string s) {
            std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
            return s;
        };
        auto split_local = [](const std::string& s, char delim) {
            std::vector<std::string> out;
            std::stringstream ss(s);
            std::string item;
            while (std::getline(ss, item, delim)) {
                // trim
                auto is_space = [](unsigned char c) { return std::isspace(c) != 0; };
                size_t b = 0;
                while (b < item.size() && is_space(static_cast<unsigned char>(item[b]))) ++b;
                size_t e = item.size();
                while (e > b && is_space(static_cast<unsigned char>(item[e - 1]))) --e;
                std::string t = item.substr(b, e - b);
                if (!t.empty()) out.push_back(t);
            }
            return out;
        };

        seismo_out_cfg_.enabled = reader.getBool("SEISMOMETERS", "enabled", true);
        seismo_out_cfg_.output_dir = reader.getString("SEISMOMETERS", "output_dir", seismo_out_cfg_.output_dir);
        seismo_out_cfg_.start_time_utc = reader.getString("SEISMOMETERS", "start_time_utc", seismo_out_cfg_.start_time_utc);

        std::string formats = reader.getString("SEISMOMETERS", "formats", "SAC,MSEED");
        formats = strtoupper_local(formats);
        seismo_out_cfg_.write_sac = (formats.find("SAC") != std::string::npos);
        seismo_out_cfg_.write_mseed = (formats.find("MSEED") != std::string::npos) || (formats.find("MINISEED") != std::string::npos);

        double default_sr = reader.getDouble("SEISMOMETERS", "default_sample_rate_hz", 100.0);
        std::string default_quantity = strtoupper_local(reader.getString("SEISMOMETERS", "default_quantity", "VELOCITY"));

        auto parseQuantity = [&](const std::string& q) {
            std::string uq = strtoupper_local(q);
            if (uq == "DISPLACEMENT") return SeismoQuantity::DISPLACEMENT;
            if (uq == "ACCELERATION") return SeismoQuantity::ACCELERATION;
            return SeismoQuantity::VELOCITY;
        };

        // Station sections: [SEISMOMETER_*]
        for (const auto& sec : reader.getSections()) {
            if (sec.find("SEISMOMETER_") != 0) continue;

            SeismometerSpec s;
            s.network = reader.getString(sec, "net", reader.getString(sec, "network", s.network));
            s.station = reader.getString(sec, "sta", reader.getString(sec, "station", s.station));
            s.location = reader.getString(sec, "loc", reader.getString(sec, "location", s.location));
            s.sample_rate_hz = reader.getDouble(sec, "sample_rate_hz", default_sr);
            s.quantity = parseQuantity(reader.getString(sec, "quantity", default_quantity));

            // Channels: "BHN,BHE,BHZ"
            std::string ch = reader.getString(sec, "channels", "");
            if (!ch.empty()) {
                auto parts = split_local(ch, ',');
                if (parts.size() >= 3) {
                    s.channels = {parts[0], parts[1], parts[2]};
                }
            }

            // Location by model coords
            if (reader.hasKey(sec, "location_xyz")) {
                auto v = reader.getDoubleArray(sec, "location_xyz");
                if (v.size() >= 3) {
                    s.coord_type = SeismoCoordinateType::MODEL_XYZ;
                    s.x = v[0]; s.y = v[1]; s.z = v[2];
                }
            } else if (reader.hasKey(sec, "location")) {
                // Back-compat: location = x,y,z
                auto v = reader.getDoubleArray(sec, "location");
                if (v.size() >= 3) {
                    s.coord_type = SeismoCoordinateType::MODEL_XYZ;
                    s.x = v[0]; s.y = v[1]; s.z = v[2];
                }
            }

            // Location by grid indices: grid = i,j,k
            if (reader.hasKey(sec, "grid")) {
                auto v = reader.getDoubleArray(sec, "grid");
                if (v.size() >= 3) {
                    s.coord_type = SeismoCoordinateType::GRID_INDEX;
                    s.i = static_cast<int>(std::llround(v[0]));
                    s.j = static_cast<int>(std::llround(v[1]));
                    s.k = static_cast<int>(std::llround(v[2]));
                    std::string gm = strtoupper_local(reader.getString(sec, "grid_mode", "CELL_CENTER"));
                    s.grid_mode = (gm == "NODE") ? GridIndexMode::NODE : GridIndexMode::CELL_CENTER;
                }
            }

            // Location by geographic coords: geo = lon,lat,elev
            if (reader.hasKey(sec, "geo")) {
                auto v = reader.getDoubleArray(sec, "geo");
                if (v.size() >= 2) {
                    s.coord_type = SeismoCoordinateType::GEOGRAPHIC;
                    s.lon = v[0];
                    s.lat = v[1];
                    s.elev = (v.size() >= 3) ? v[2] : 0.0;
                }
            }

            // Instrument performance knobs
            s.instrument.highpass_corner_hz = reader.getDouble(sec, "highpass_corner_hz", 0.0);
            s.instrument.lowpass_corner_hz = reader.getDouble(sec, "lowpass_corner_hz", 0.0);
            s.instrument.noise_std = reader.getDouble(sec, "noise_std", 0.0);
            s.instrument.gain = reader.getDouble(sec, "gain", 1.0);
            s.instrument.adc_bits = reader.getInt(sec, "adc_bits", 0);
            s.instrument.full_scale = reader.getDouble(sec, "full_scale", 0.0);
            s.instrument.clip = reader.getDouble(sec, "clip", 0.0);

            seismo_specs_.push_back(std::move(s));
        }
    }
    
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
        config.enable_near_field_damage =
            config.enable_near_field_damage || reader.getBool("EXPLOSION_SOURCE", "apply_damage_zone", false);
        std::string ex_type = reader.getString("EXPLOSION_SOURCE", "type", "");
        // Currently only support underground nuclear coupling via spherical cavity proxy.
        if (!ex_type.empty() && (ex_type == "NUCLEAR_UNDERGROUND" || ex_type == "UNDERGROUND_NUCLEAR" || ex_type == "UNDERGROUND_CONTAINED")) {
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

            // COUPLED_ANALYTIC: run 1D NearField solver and store RDP source
            if (config.explosion_solve_mode == "COUPLED_ANALYTIC") {
                UndergroundExplosionSource nf_src;
                nf_src.yield_kt = yield_kt;
                nf_src.depth = depth_m;
                nf_src.location = {x, y, z};
                nf_src.host_density = rho;
                nf_src.host_vp = vp;
                nf_src.host_vs = vs;
                nf_src.rise_time = rise;
                nf_src.overburden_stress = rho * 9.81 * depth_m;

                explosion_->nf_solver.setSource(nf_src);
                explosion_->nf_solver.initialize();

                // Run the 1D solver to completion (fast, sub-second)
                double nf_dt = 0.0001;
                double nf_end = std::max(1.0, 10.0 * rise);
                for (double nf_t = 0.0; nf_t < nf_end; nf_t += nf_dt) {
                    explosion_->nf_solver.prestep(nf_t, nf_dt);
                    explosion_->nf_solver.step(nf_dt);
                    explosion_->nf_solver.poststep(nf_t + nf_dt, nf_dt);
                }

                // Store zone radii from the 1D solver
                explosion_->nf_cavity_radius = explosion_->nf_solver.getCavityRadius(nf_end);
                explosion_->nf_crushed_radius = explosion_->nf_solver.getCrushedZoneRadius();
                explosion_->nf_fractured_radius = explosion_->nf_solver.getFracturedZoneRadius();

                // Initialize RDP source from the same explosion parameters
                explosion_->rdp_source.initialize(nf_src);
                explosion_->use_nearfield_coupling = true;

                if (rank == 0) {
                    PetscPrintf(comm, "COUPLED_ANALYTIC: NearField 1D solver completed\n");
                    PetscPrintf(comm, "  RDP psi_inf = %.3e m^3, fc = %.3f Hz\n",
                                explosion_->rdp_source.psi_infinity,
                                explosion_->rdp_source.corner_frequency);
                    PetscPrintf(comm, "  Cavity R = %.1f m, Crushed R = %.1f m, Fractured R = %.1f m\n",
                                explosion_->nf_cavity_radius,
                                explosion_->nf_crushed_radius,
                                explosion_->nf_fractured_radius);
                }
            }
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
    
    // =========================================================================
    // GPU Initialization
    // =========================================================================
    // When use_gpu is requested, attempt to initialize the GPU backend.
    // Depending on gpu_mode, either hard-fail (GPU_ONLY) or fall back to CPU.
#ifdef USE_CUDA
    if (config.use_gpu) {
        GPUManager& gpu = GPUManager::getInstance();
        if (gpu.initialize()) {
            if (config.gpu_device_id >= 0) {
                gpu.setDevice(config.gpu_device_id);
            }
            if (rank == 0) {
                PetscPrintf(comm, "\n");
                PetscPrintf(comm, "============================================================\n");
                PetscPrintf(comm, "  GPU Acceleration Enabled\n");
                PetscPrintf(comm, "============================================================\n");
                gpu.printDeviceInfo();
                PetscPrintf(comm, "\n");
            }
        } else {
            if (config.gpu_mode == GPUExecutionMode::GPU_ONLY) {
                SETERRQ(comm, PETSC_ERR_SUP,
                        "GPU_ONLY mode requested but no CUDA-capable GPU is available. "
                        "Set gpu_mode = CPU_FALLBACK to allow CPU execution.");
            }
            config.use_gpu = false;
            if (rank == 0) {
                PetscPrintf(comm, "Warning: GPU requested but not available. Falling back to CPU.\n");
            }
        }
    }
#else
    if (config.use_gpu) {
        if (rank == 0) {
            PetscPrintf(comm, "Warning: use_gpu = true in config, but FSRM was built without "
                              "CUDA support (ENABLE_CUDA=OFF). All physics will run on CPU.\n");
        }
        config.use_gpu = false;
    }
#endif

    if (rank == 0) {
        PetscPrintf(comm, "Configuration loaded successfully\n");
        if (config.use_gpu) {
            PetscPrintf(comm, "GPU execution: ENABLED (mode=%s)\n",
                        config.gpu_mode == GPUExecutionMode::GPU_ONLY ? "GPU_ONLY" :
                        config.gpu_mode == GPUExecutionMode::HYBRID ? "HYBRID" :
                        config.gpu_mode == GPUExecutionMode::CPU_FALLBACK ? "CPU_FALLBACK" :
                        "CPU_ONLY");
        } else {
            PetscPrintf(comm, "GPU execution: DISABLED (CPU-only mode)\n");
        }
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
    std::vector<std::pair<std::string, std::string>> materialMappings = grid_config.gmsh_material_domains;
    if (materialMappings.empty() && !config.material_regions.empty()) {
        for (const auto& region : config.material_regions) {
            materialMappings.emplace_back(region.gmsh_label_name, region.gmsh_label_name);
        }
    }

    const bool have_mappings =
        !materialMappings.empty() ||
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
            std::unordered_map<int, int> physTagToMaterialId;
            int nextMaterialId = 0;
            for (const auto& mapping : materialMappings) {
                auto pit = physicalNames.find(mapping.first);
                if (pit == physicalNames.end()) continue;
                if (!physTagToMaterialId.count(pit->second.tag)) {
                    physTagToMaterialId[pit->second.tag] = nextMaterialId++;
                }
            }
            for (auto& region : config.material_regions) {
                auto pit = physicalNames.find(region.gmsh_label_name);
                if (pit == physicalNames.end()) {
                    region.label_id = -1;
                    continue;
                }
                auto mit = physTagToMaterialId.find(pit->second.tag);
                region.label_id = (mit != physTagToMaterialId.end()) ? mit->second : -1;
            }

            ierr = applyGmshNameMappingsToDM(
                comm, dm, physicalNames,
                materialMappings,
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

    // NOTE: PETSc 3.20 cohesive cell insertion works reliably only on simplex meshes.
    // Force simplex (tetrahedral) mesh when faults are enabled, even for CARTESIAN configs.
    PetscBool use_simplex = config.enable_faults ? PETSC_TRUE : PETSC_FALSE;

    ierr = DMPlexCreateBoxMesh(comm, 3, use_simplex,
                              faces,
                              lower, upper, nullptr, PETSC_TRUE, 0, PETSC_FALSE, &dm); CHKERRQ(ierr);

    ierr = DMSetFromOptions(dm); CHKERRQ(ierr);
    ierr = DMSetUp(dm); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::setupUnstructuredGrid() {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

    // Create DMPlex for unstructured grid
    PetscInt faces[3] = {grid_config.nx, grid_config.ny, grid_config.nz};

    // NOTE: PETSc 3.20 cohesive cell insertion works reliably only on simplex meshes.
    // Force simplex (tetrahedral) mesh when faults are enabled.
    PetscBool use_simplex = config.enable_faults ? PETSC_TRUE : PETSC_FALSE;

    ierr = DMPlexCreateBoxMesh(comm, 3, use_simplex,
                              faces,
                              nullptr, nullptr, nullptr, PETSC_TRUE, 0, PETSC_FALSE, &dm); CHKERRQ(ierr);

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

    // Initialize seismometers after fields exist (so we can subselect displacement_)
    if (seismo_out_cfg_.enabled && !seismo_specs_.empty()) {
        seismometers_ = std::make_unique<SeismometerNetwork>(comm);
        seismometers_->setOutputConfig(seismo_out_cfg_);
        seismometers_->setStations(seismo_specs_);
        ierr = seismometers_->initialize(dm, grid_config, coord_manager.get()); CHKERRQ(ierr);
    }
    
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::createFieldsFromConfig() {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    // Create discretization based on fluid model
    // Use PetscFECreateLagrange with degree 1 (Q1/P1) for continuous Galerkin FEM.
    // This places DOFs on vertices, enabling proper Dirichlet BC enforcement.
    PetscFE fe;

    // Detect if the mesh is simplex (tets) or tensor-product (hexes).
    // When faults are enabled the mesh is forced to simplex in setupDM().
    PetscBool isSimplex = PETSC_FALSE;
    {
        PetscInt cStart, cEnd;
        ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd); CHKERRQ(ierr);
        if (cStart < cEnd) {
            DMPolytopeType ct;
            ierr = DMPlexGetCellType(dm, cStart, &ct); CHKERRQ(ierr);
            isSimplex = (ct == DM_POLYTOPE_TETRAHEDRON || ct == DM_POLYTOPE_TRIANGLE)
                        ? PETSC_TRUE : PETSC_FALSE;
        }
    }
    
    switch (config.fluid_model) {
        case FluidModelType::NONE:
            // No fluid fields
            break;
        case FluidModelType::SINGLE_COMPONENT:
            // Single pressure field
            ierr = PetscFECreateLagrange(comm, 3, 1, isSimplex, 1, -1, &fe); CHKERRQ(ierr);
            ierr = PetscObjectSetName((PetscObject)fe, "pressure_"); CHKERRQ(ierr);
            ierr = DMAddField(dm, nullptr, (PetscObject)fe); CHKERRQ(ierr);
            fe_fields.push_back(fe);
            break;
            
        case FluidModelType::BLACK_OIL:
            // Pressure + saturations
            ierr = PetscFECreateLagrange(comm, 3, 1, isSimplex, 1, -1, &fe); CHKERRQ(ierr);
            ierr = PetscObjectSetName((PetscObject)fe, "pressure_"); CHKERRQ(ierr);
            ierr = DMAddField(dm, nullptr, (PetscObject)fe); CHKERRQ(ierr);
            fe_fields.push_back(fe);
            
            ierr = PetscFECreateLagrange(comm, 3, 1, isSimplex, 1, -1, &fe); CHKERRQ(ierr);
            ierr = PetscObjectSetName((PetscObject)fe, "saturation_w_"); CHKERRQ(ierr);
            ierr = DMAddField(dm, nullptr, (PetscObject)fe); CHKERRQ(ierr);
            fe_fields.push_back(fe);
            
            ierr = PetscFECreateLagrange(comm, 3, 1, isSimplex, 1, -1, &fe); CHKERRQ(ierr);
            ierr = PetscObjectSetName((PetscObject)fe, "saturation_g_"); CHKERRQ(ierr);
            ierr = DMAddField(dm, nullptr, (PetscObject)fe); CHKERRQ(ierr);
            fe_fields.push_back(fe);
            break;
            
        case FluidModelType::COMPOSITIONAL:
            // Pressure + compositions
            ierr = PetscFECreateLagrange(comm, 3, 1, isSimplex, 1, -1, &fe); CHKERRQ(ierr);
            ierr = PetscObjectSetName((PetscObject)fe, "pressure_"); CHKERRQ(ierr);
            ierr = DMAddField(dm, nullptr, (PetscObject)fe); CHKERRQ(ierr);
            fe_fields.push_back(fe);
            
            // Add composition fields (number depends on components)
            for (int i = 0; i < 3; ++i) {  // Example: 3 components
                char name[256];
                snprintf(name, sizeof(name), "composition_%d_", i);
                ierr = PetscFECreateLagrange(comm, 3, 1, isSimplex, 1, -1, &fe); CHKERRQ(ierr);
                ierr = PetscObjectSetName((PetscObject)fe, name); CHKERRQ(ierr);
                ierr = DMAddField(dm, nullptr, (PetscObject)fe); CHKERRQ(ierr);
                fe_fields.push_back(fe);
            }
            break;
            
        default:
            SETERRQ(comm, PETSC_ERR_ARG_WRONG, "Unknown fluid model type");
    }
    
    // Add displacement field if geomechanics or elastodynamics is enabled
    if (config.enable_geomechanics || config.enable_elastodynamics) {
        ierr = PetscFECreateLagrange(comm, 3, 3, isSimplex, 1, -1, &fe); CHKERRQ(ierr);
        ierr = PetscObjectSetName((PetscObject)fe, "displacement_"); CHKERRQ(ierr);
        ierr = DMAddField(dm, nullptr, (PetscObject)fe); CHKERRQ(ierr);
        fe_fields.push_back(fe);
    }

    // Add Lagrange multiplier field for cohesive faults if enabled.
    // The field is on ALL cells (nullptr label). PETSc 3.25 region DS with
    // label-restricted fields does not support volume residual assembly via
    // DMPlexTSComputeIFunctionFEM. Instead, a weak regularization (epsilon=1e-8)
    // keeps interior Lagrange DOFs non-singular for LU while being negligible
    // compared to the BdResidual constraint on cohesive faces.
    PetscInt lagrange_field_idx = -1;
    if (config.enable_faults && cohesive_kernel_) {
        PetscFE fe_lagrange;
        // Lagrange multiplier: 3-component vector (traction) on cohesive cells
        ierr = PetscFECreateLagrange(comm, 3, 3, isSimplex, 1, -1, &fe_lagrange); CHKERRQ(ierr);
        ierr = PetscObjectSetName((PetscObject)fe_lagrange, "lagrange_"); CHKERRQ(ierr);
        lagrange_field_idx = static_cast<PetscInt>(fe_fields.size());
        ierr = DMAddField(dm, nullptr, (PetscObject)fe_lagrange); CHKERRQ(ierr);
        fe_fields.push_back(fe_lagrange);
        if (rank == 0) {
            PetscPrintf(comm, "Added Lagrange multiplier field for fault traction\n");
        }
    }

    // Add thermal field if enabled
    if (config.enable_thermal) {
        ierr = PetscFECreateLagrange(comm, 3, 1, isSimplex, 1, -1, &fe); CHKERRQ(ierr);
        ierr = PetscObjectSetName((PetscObject)fe, "temperature_"); CHKERRQ(ierr);
        ierr = DMAddField(dm, nullptr, (PetscObject)fe); CHKERRQ(ierr);
        fe_fields.push_back(fe);
    }

    // Create DS (required before adding boundaries in PETSc 3.25)
    ierr = DMCreateDS(dm); CHKERRQ(ierr);

    // Add boundary conditions to the DS
    ierr = setupBoundaryConditions(); CHKERRQ(ierr);

    // Force PETSc to rebuild the local section with BC constraints
    ierr = DMSetLocalSection(dm, nullptr); CHKERRQ(ierr);
    ierr = DMSetUp(dm); CHKERRQ(ierr);

    // Get the DS for later use
    ierr = DMGetDS(dm, &prob); CHKERRQ(ierr);

    if (rank == 0 && config.enable_faults && cohesive_kernel_) {
        PetscPrintf(comm, "Fields setup complete with Lagrange multiplier field for fault traction\n");
    }

    PetscFunctionReturn(0);
}

// ============================================================================
// Weak Lagrange volume regularization callbacks
// ============================================================================

// Weak regularization: f = epsilon * lambda, g = epsilon * I.
// epsilon = 1e-4 is small enough that the regularization force
// (epsilon * lambda ~ 1e-4 * stress) is negligible compared to the
// BdResidual constraint (O(stress * area / n_nodes)) on ALL meshes.
// The old epsilon=1 competed with BdResidual on coarse meshes.
//
// The small epsilon Jacobian diagonal (1e-4) at INTERIOR Lagrange DOFs
// is complemented by a penalty-scaled diagonal added manually at COHESIVE
// vertices in addCohesivePenaltyToJacobian. This gives:
//   interior: diag = 1e-4 (small but nonzero for LU)
//   cohesive: diag = 1e-4 + penalty*coeff (dominated by penalty ~ E/h)
// Condition number: stiffness / epsilon = O(1e9 / 1e-4) = O(1e13),
// within double precision (16 digits).
static void f0_weak_lagrange(PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[],
    const PetscScalar u_x[], const PetscInt aOff[],
    const PetscInt aOff_x[], const PetscScalar a[],
    const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, const PetscReal x[],
    PetscInt numConstants, const PetscScalar constants[],
    PetscScalar f[])
{
    (void)Nf; (void)NfAux; (void)uOff_x; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)x; (void)numConstants; (void)constants;
    const PetscScalar epsilon = 1.0e-4;
    const PetscInt lagr_off = uOff[Nf - 1];
    for (PetscInt d = 0; d < dim; ++d) {
        f[d] = epsilon * u[lagr_off + d];
    }
}

static void g0_weak_lagrange(PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[],
    const PetscScalar u_x[], const PetscInt aOff[],
    const PetscInt aOff_x[], const PetscScalar a[],
    const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, PetscReal u_tShift, const PetscReal x[],
    PetscInt numConstants, const PetscScalar constants[],
    PetscScalar g0[])
{
    (void)Nf; (void)NfAux; (void)uOff; (void)uOff_x; (void)u; (void)u_t;
    (void)u_x; (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)u_tShift; (void)x; (void)numConstants; (void)constants;
    const PetscScalar epsilon = 1.0e-4;
    for (PetscInt d = 0; d < dim; ++d) {
        g0[d * dim + d] = epsilon;
    }
}

// ============================================================================
// Boundary condition callback functions
// ============================================================================

// BC callback: all components = 0 (for fixed boundaries)
static PetscErrorCode bc_zero(PetscInt dim, PetscReal time, const PetscReal x[],
                               PetscInt Nc, PetscScalar *u, void *ctx) {
    (void)dim; (void)time; (void)x; (void)ctx;  // Suppress unused warnings
    for (PetscInt c = 0; c < Nc; c++) {
        u[c] = 0.0;
    }
    return PETSC_SUCCESS;
}

// BC callback: applied compression on top (u_z = -0.001, u_x = u_y = 0)
static PetscErrorCode bc_compression(PetscInt dim, PetscReal time, const PetscReal x[],
                                      PetscInt Nc, PetscScalar *u, void *ctx) {
    (void)dim; (void)time; (void)x; (void)ctx;  // Suppress unused warnings
    if (Nc >= 3) {
        u[0] = 0.0;   // u_x = 0
        u[1] = 0.0;   // u_y = 0
        u[2] = -0.001; // u_z = -1mm (compression)
    }
    return PETSC_SUCCESS;
}

// BC callback: drained boundary (pressure = 0)
static PetscErrorCode bc_drained(PetscInt dim, PetscReal time, const PetscReal x[],
                                  PetscInt Nc, PetscScalar *u, void *ctx) {
    (void)dim; (void)time; (void)x; (void)Nc; (void)ctx;  // Suppress unused warnings
    u[0] = 0.0; // pressure = 0
    return PETSC_SUCCESS;
}

// BC callback: prescribed Dirichlet values read from context
// ctx points to a static double[3] array with the prescribed displacement values
static PetscErrorCode bc_dirichlet_values(PetscInt dim, PetscReal time, const PetscReal x[],
                                           PetscInt Nc, PetscScalar *u, void *ctx) {
    (void)dim; (void)time; (void)x;
    const double *vals = static_cast<const double *>(ctx);
    for (PetscInt c = 0; c < Nc; c++) {
        u[c] = vals ? vals[c] : 0.0;
    }
    return PETSC_SUCCESS;
}

// BC callback: Dirichlet pressure (single scalar field).
// The ctx pointer points to a double holding the pressure value.
static PetscErrorCode bc_pressure_dirichlet(PetscInt dim, PetscReal time, const PetscReal x[],
                                             PetscInt Nc, PetscScalar *u, void *ctx) {
    (void)dim; (void)time; (void)x; (void)Nc;
    const double *pval = static_cast<const double *>(ctx);
    u[0] = pval ? *pval : 0.0;
    return PETSC_SUCCESS;
}

// Traction BC boundary residual callback (f0)
// Reads traction vector from constants array using the face normal to determine face index.
// Face ordering: 0=x_min, 1=x_max, 2=y_min, 3=y_max, 4=z_min, 5=z_max
// Traction constants start at slot TractionBC::TRACTION_CONST_BASE (36).
static void f0_traction_bc(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                            const PetscInt uOff[], const PetscInt uOff_x[],
                            const PetscScalar u[], const PetscScalar u_t[],
                            const PetscScalar u_x[],
                            const PetscInt aOff[], const PetscInt aOff_x[],
                            const PetscScalar a[], const PetscScalar a_t[],
                            const PetscScalar a_x[],
                            PetscReal t, const PetscReal x[],
                            const PetscReal n[],
                            PetscInt numConstants,
                            const PetscScalar constants[],
                            PetscScalar f[]) {
    (void)Nf; (void)NfAux; (void)uOff; (void)uOff_x;
    (void)u; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)x;
    if (numConstants < FSRM::TractionBC::TRACTION_CONST_COUNT) {
        for (PetscInt d = 0; d < dim; ++d) f[d] = 0.0;
        return;
    }
    // Determine face index from outward normal direction
    PetscInt face_idx = 0;
    PetscReal max_abs = 0.0;
    for (PetscInt d = 0; d < dim; ++d) {
        PetscReal absn = PetscAbsReal(PetscRealPart(n[d]));
        if (absn > max_abs) {
            max_abs = absn;
            face_idx = d * 2 + (PetscRealPart(n[d]) > 0.0 ? 1 : 0);
        }
    }
    PetscInt base = FSRM::TractionBC::TRACTION_CONST_BASE + face_idx * 3;
    for (PetscInt d = 0; d < dim; ++d) {
        f[d] = PetscRealPart(constants[base + d]);
    }
}

// IC callback: Terzaghi consolidation initial conditions
// Nc=1 for pressure field (set to 1 MPa undrained response), Nc=3 for displacement (zero)
static PetscErrorCode terzaghi_ic(PetscInt dim, PetscReal time, const PetscReal x[],
                                   PetscInt Nc, PetscScalar *u, void *ctx) {
    (void)dim; (void)time; (void)x; (void)ctx;
    if (Nc == 1) {
        u[0] = 1.0e6;  // 1 MPa initial pore pressure
    } else {
        for (PetscInt c = 0; c < Nc; c++) u[c] = 0.0;
    }
    return PETSC_SUCCESS;
}

PetscErrorCode Simulator::setupPhysics() {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    // =========================================================================
    // Helper lambda: create a GPU kernel if use_gpu is enabled and the GPU
    // implementation exists, otherwise fall back to the CPU kernel.
    // =========================================================================
    auto tryGPUKernel = [&](PhysicsType ptype) -> std::shared_ptr<PhysicsKernel> {
#ifdef USE_CUDA
        if (config.use_gpu) {
            auto gpu_kernel = createGPUKernelExtended(ptype);
            if (gpu_kernel) {
                if (rank == 0) {
                    PetscPrintf(comm, "  GPU kernel created for %s\n",
                                gpu_kernel->getTypeName());
                }
                return gpu_kernel;
            }
            if (rank == 0 && config.gpu_verbose) {
                PetscPrintf(comm, "  No GPU kernel for physics type %d, using CPU\n",
                            static_cast<int>(ptype));
            }
        }
#else
        (void)ptype;  // suppress unused warning
#endif
        return nullptr;
    };

    // Add physics kernels based on configuration
    switch (config.fluid_model) {
        case FluidModelType::NONE:
            break;
        case FluidModelType::SINGLE_COMPONENT: {
            auto kernel = tryGPUKernel(PhysicsType::FLUID_FLOW);
            if (!kernel) kernel = std::make_shared<SinglePhaseFlowKernel>();
            ierr = addPhysicsKernel(kernel); CHKERRQ(ierr);
            break;
        }
        case FluidModelType::BLACK_OIL: {
            auto kernel = tryGPUKernel(PhysicsType::FLUID_FLOW);
            if (!kernel) kernel = std::make_shared<BlackOilKernel>();
            ierr = addPhysicsKernel(kernel); CHKERRQ(ierr);
            break;
        }
        case FluidModelType::COMPOSITIONAL: {
            // No GPU kernel for compositional yet; use CPU
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
            auto gpu_kernel = tryGPUKernel(PhysicsType::ELASTODYNAMICS);
            auto kernel = gpu_kernel ? std::dynamic_pointer_cast<ElastodynamicsKernel>(gpu_kernel)
                                     : std::make_shared<ElastodynamicsKernel>();
            if (!kernel) kernel = std::make_shared<ElastodynamicsKernel>();
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
            auto gpu_kernel = tryGPUKernel(PhysicsType::POROELASTODYNAMICS);
            auto kernel = gpu_kernel ? std::dynamic_pointer_cast<PoroelastodynamicsKernel>(gpu_kernel)
                                     : std::make_shared<PoroelastodynamicsKernel>();
            if (!kernel) kernel = std::make_shared<PoroelastodynamicsKernel>();
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
            auto gpu_kernel = tryGPUKernel(PhysicsType::GEOMECHANICS);
            if (gpu_kernel) {
                ierr = addPhysicsKernel(gpu_kernel); CHKERRQ(ierr);
            } else {
                auto kernel = std::make_shared<GeomechanicsKernel>(config.solid_model);
                ierr = addPhysicsKernel(kernel); CHKERRQ(ierr);
            }
        }
    }
    
    // Add standalone elastodynamics if enabled (wave equation solver)
    if (config.enable_elastodynamics) {
        auto gpu_kernel = tryGPUKernel(PhysicsType::ELASTODYNAMICS);
        auto kernel = gpu_kernel ? std::dynamic_pointer_cast<ElastodynamicsKernel>(gpu_kernel)
                                 : std::make_shared<ElastodynamicsKernel>();
        if (!kernel) kernel = std::make_shared<ElastodynamicsKernel>();
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
    
    // Add standalone poroelastodynamics if enabled (coupled fluid-solid waves)
    if (config.enable_poroelastodynamics) {
        auto gpu_kernel = tryGPUKernel(PhysicsType::POROELASTODYNAMICS);
        auto kernel = gpu_kernel ? std::dynamic_pointer_cast<PoroelastodynamicsKernel>(gpu_kernel)
                                 : std::make_shared<PoroelastodynamicsKernel>();
        if (!kernel) kernel = std::make_shared<PoroelastodynamicsKernel>();
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
        auto gpu_kernel = tryGPUKernel(PhysicsType::THERMAL);
        auto kernel = gpu_kernel ? gpu_kernel : std::make_shared<ThermalKernel>();
        ierr = addPhysicsKernel(kernel); CHKERRQ(ierr);
    }
    
    // Add particle transport if enabled
    if (config.enable_particle_transport) {
        auto gpu_kernel = tryGPUKernel(PhysicsType::PARTICLE_TRANSPORT);
        auto kernel = gpu_kernel ? gpu_kernel : std::make_shared<ParticleTransportKernel>();
        ierr = addPhysicsKernel(kernel); CHKERRQ(ierr);
    }
    
    // Setup residuals for PETSc FEM time integration.
    // NOTE: At present we only have a consistent PETSc FEM residual for the
    // single-phase pressure equation. Other physics kernels are not wired into
    // DMPlexTSComputeIFunctionFEM yet, so we must not enable the FEM path or
    // PETSc will crash during assembly.
    use_fem_time_residual_ = false;

    // Unified PetscDS constants array (shared by all physics callbacks)
    // Layout:
    //   [0]  = lambda (first Lame parameter)
    //   [1]  = mu (shear modulus)
    //   [2]  = rho_solid
    //   [3]  = porosity
    //   [4]  = kx (permeability x, m^2)
    //   [5]  = ky (permeability y, m^2)
    //   [6]  = kz (permeability z, m^2)
    //   [7]  = cw (water compressibility, 1/Pa)
    //   [8]  = co (oil compressibility, 1/Pa)
    //   [9]  = cg (gas compressibility, 1/Pa)
    //   [10] = mu_w (water viscosity, Pa*s)
    //   [11] = mu_o (oil viscosity, Pa*s)
    //   [12] = mu_g (gas viscosity, Pa*s)
    //   [13] = Swr (water residual saturation)
    //   [14] = Sor (oil residual saturation)
    //   [15] = Sgr (gas residual saturation)
    //   [16] = nw (Corey exponent water)
    //   [17] = no (Corey exponent oil)
    //   [18] = ng (Corey exponent gas)
    //   [19] = krw0 (max relative perm water)
    //   [20] = kro0 (max relative perm oil)
    //   [21] = krg0 (max relative perm gas)
    //   [22] = biot_coefficient (alpha)
    //   [23] = biot_modulus_inv (1/M)
    //   [24] = rho_fluid
    //   [36-53] = traction BC values (6 faces x 3 components)
    // Maximum constant array size: must accommodate all physics modules.
    // TractionBC uses up to slot 53, Viscoelastic uses up to slot 69.
    static constexpr PetscInt MAX_UNIFIED_CONSTANTS = 80;
    PetscScalar unified_constants[MAX_UNIFIED_CONSTANTS] = {0};

    // Always set material properties if geomechanics is enabled
    if (config.enable_geomechanics || config.enable_elastodynamics) {
        const MaterialProperties& mat = material_props.empty() ? MaterialProperties() : material_props[0];
        double E = mat.youngs_modulus;
        double nu = mat.poisson_ratio;
        double rho = mat.density;

        // Convert Young's modulus and Poisson's ratio to Lame parameters
        double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        double mu = E / (2.0 * (1.0 + nu));

        unified_constants[0] = lambda;
        unified_constants[1] = mu;
        unified_constants[2] = rho;
        unified_constants[22] = mat.biot_coefficient;
    }

    // Set fluid properties if fluid flow is enabled
    if (config.fluid_model != FluidModelType::NONE) {
        const MaterialProperties mat = material_props.empty() ? MaterialProperties() : material_props[0];
        const FluidProperties flu = fluid_props.empty() ? FluidProperties() : fluid_props[0];

        // ConfigReader::getDoubleWithUnit already converts permeability from
        // mD to m^2 via UnitSystem::toBase("mD"). Do NOT multiply by mD_to_m2
        // again. The MaterialProperties default values (100.0) are expressed
        // in mD and must be converted when no config file is loaded.
        constexpr double mD_to_m2 = 9.869233e-16;
        bool from_config = !material_props.empty();

        unified_constants[3]  = mat.porosity;
        unified_constants[4]  = from_config ? mat.permeability_x : mat.permeability_x * mD_to_m2;
        unified_constants[5]  = from_config ? mat.permeability_y : mat.permeability_y * mD_to_m2;
        unified_constants[6]  = from_config ? mat.permeability_z : mat.permeability_z * mD_to_m2;
        unified_constants[7]  = flu.water_compressibility;
        unified_constants[8]  = flu.oil_compressibility;
        unified_constants[9]  = flu.gas_compressibility;
        unified_constants[10] = flu.water_viscosity;
        unified_constants[11] = flu.oil_viscosity;
        unified_constants[12] = flu.gas_viscosity;
        unified_constants[13] = flu.water_residual_saturation;
        unified_constants[14] = flu.residual_saturation;       // Sor
        unified_constants[15] = flu.gas_residual_saturation;
        unified_constants[16] = flu.corey_exponent_water;
        unified_constants[17] = flu.corey_exponent_oil;
        unified_constants[18] = flu.corey_exponent_gas;
        unified_constants[19] = flu.kr_max_water;
        unified_constants[20] = flu.kr_max_oil;
        unified_constants[21] = flu.kr_max_gas;
        unified_constants[24] = flu.density;

        // Biot modulus: 1/M = phi * cf (simplified, assumes incompressible grains)
        if (config.enable_geomechanics) {
            unified_constants[23] = mat.porosity * flu.water_compressibility;
        }
    }

    // When using heterogeneous material or gravity with homogeneous material,
    // material properties (lambda, mu, rho) come from auxiliary fields instead of
    // constants. Repurpose constants[0] as gravity magnitude.
    bool use_aux_callbacks = config.use_heterogeneous_material ||
                             config.enable_near_field_damage ||
                             config.enable_elastoplasticity ||
                             config.enable_viscoelastic ||
                             (config.gravity > 0.0);
    if (use_aux_callbacks) {
        // When absorbing BCs are enabled with gravity, the constants[0-2] slots
        // are repurposed for gravity and cannot hold lambda/mu/rho. The AbsorbingBC
        // callback must read material properties from aux fields.
        // Aux fields are always populated when use_aux_callbacks is true
        // (setupAuxiliaryDM + populateAuxFieldsByDepth), so this is safe.
        if (config.absorbing_bc_enabled) {
            PetscCheck(!material_props.empty(), PETSC_COMM_WORLD, PETSC_ERR_USER,
                "Absorbing BCs with gravity require material properties to populate auxiliary fields");
        }
        unified_constants[0] = config.gravity;
        unified_constants[1] = 0.0;
        unified_constants[2] = 0.0;
    }

    if (hydrofrac_fem_pressurized_mode_) {
        unified_constants[PetscFEHydrofrac::HYDROFRAC_CONST_PRESSURE] = hydrofrac_uniform_pressure_pa_;
    }

    // Set elastoplasticity constants if enabled
    if (config.enable_elastoplasticity) {
        unified_constants[PetscFEElastoplasticity::EP_CONST_COHESION] = config.ep_cohesion;
        unified_constants[PetscFEElastoplasticity::EP_CONST_FRICTION_ANGLE] = config.ep_friction_angle;
        unified_constants[PetscFEElastoplasticity::EP_CONST_DILATION_ANGLE] = config.ep_dilation_angle;
        unified_constants[PetscFEElastoplasticity::EP_CONST_HARDENING_MODULUS] = config.ep_hardening_modulus;
    }

    // Set viscoelastic constants if enabled
    if (config.enable_viscoelastic) {
        int N = config.visco_num_mechanisms;
        double tau_arr[5] = {0}, dmu_arr[5] = {0};
        PetscFEViscoelastic::computeMechanismWeights(N,
            config.visco_f_min, config.visco_f_max, config.visco_q_s,
            tau_arr, dmu_arr);
        unified_constants[PetscFEViscoelastic::VISCO_CONST_N] = static_cast<PetscScalar>(N);
        for (int m = 0; m < N && m < 5; ++m) {
            unified_constants[PetscFEViscoelastic::VISCO_CONST_TAU_BASE + m] = tau_arr[m];
            unified_constants[PetscFEViscoelastic::VISCO_CONST_DMU_BASE + m] = dmu_arr[m];
            // Use same weights for bulk modulus (simplified: Q_p approx 2*Q_s)
            unified_constants[PetscFEViscoelastic::VISCO_CONST_DK_BASE + m] = dmu_arr[m];
        }
    }

    // Set cohesive fault constants if enabled
    if (config.enable_faults && cohesive_kernel_) {
        double fault_mode_val = 0.0;  // locked
        if (fault_mode_ == "slipping" || fault_mode_ == "dynamic") {
            fault_mode_val = 1.0;
        }
        unified_constants[CohesiveFaultKernel::COHESIVE_CONST_MODE] = fault_mode_val;
        unified_constants[CohesiveFaultKernel::COHESIVE_CONST_MU_F] = fault_friction_coefficient_;
        if (use_slip_weakening_) {
            unified_constants[CohesiveFaultKernel::COHESIVE_CONST_FRICTION_MODEL] = 1.0;
            unified_constants[CohesiveFaultKernel::COHESIVE_CONST_MU_S] = mu_static_;
            unified_constants[CohesiveFaultKernel::COHESIVE_CONST_MU_D] = mu_dynamic_;
            unified_constants[CohesiveFaultKernel::COHESIVE_CONST_DC] = critical_slip_distance_;
        }
    }

    // Set traction BC constants (6 faces x 3 components, starting at slot 36)
    bool has_traction_bc = false;
    for (int fi = 0; fi < 6; ++fi) {
        if (config.face_bc[fi].configured && config.face_bc[fi].type == "traction") {
            has_traction_bc = true;
            PetscInt base = TractionBC::TRACTION_CONST_BASE + fi * 3;
            for (int d = 0; d < 3; ++d) {
                unified_constants[base + d] = config.face_bc[fi].traction[d];
            }
        }
    }

    // Set constants once for all physics (eliminates collision bug)
    // [0-24]=fluid/elasticity/poroelasticity, [25-30]=cohesive fault, [31]=hydrofrac pressure,
    // [32-35]=elastoplasticity, [36-53]=traction BC
    PetscInt nconst_unified = hydrofrac_fem_pressurized_mode_
                              ? PetscFEHydrofrac::HYDROFRAC_CONST_COUNT
                              : 27;
    if (config.enable_faults && cohesive_kernel_) {
        nconst_unified = std::max(nconst_unified,
            static_cast<PetscInt>(CohesiveFaultKernel::COHESIVE_CONST_COUNT));
    }
    if (config.enable_elastoplasticity) {
        nconst_unified = std::max(nconst_unified,
            static_cast<PetscInt>(PetscFEElastoplasticity::EP_CONST_COUNT));
    }
    if (has_traction_bc) {
        nconst_unified = std::max(nconst_unified,
            static_cast<PetscInt>(TractionBC::TRACTION_CONST_COUNT));
    }
    if (config.enable_viscoelastic) {
        nconst_unified = std::max(nconst_unified,
            static_cast<PetscInt>(PetscFEViscoelastic::VISCO_CONST_COUNT));
    }
    has_traction_bc_ = has_traction_bc;
    ierr = PetscDSSetConstants(prob, nconst_unified, unified_constants); CHKERRQ(ierr);

    // Register PetscFEPoroelasticity callbacks for fully coupled Biot poroelasticity
    if (config.solid_model == SolidModelType::POROELASTIC &&
        config.fluid_model == FluidModelType::SINGLE_COMPONENT) {
        // Two-field system: field 0 = pressure (1 component), field 1 = displacement (3 components)

        // Pressure equation: (1/M)*dp/dt + alpha*div(du/dt) + div((k/mu)*grad(p)) = 0
        ierr = PetscDSSetResidual(prob, 0,
                                  PetscFEPoroelasticity::f0_pressure,
                                  PetscFEPoroelasticity::f1_pressure); CHKERRQ(ierr);

        // Displacement equation: -div(sigma - alpha*p*I) = 0
        ierr = PetscDSSetResidual(prob, 1,
                                  PetscFEPoroelasticity::f0_displacement,
                                  PetscFEPoroelasticity::f1_displacement); CHKERRQ(ierr);

        // Jacobian block: pressure-pressure (diagonal)
        ierr = PetscDSSetJacobian(prob, 0, 0,
                                  PetscFEPoroelasticity::g0_pp, nullptr, nullptr,
                                  PetscFEPoroelasticity::g3_pp); CHKERRQ(ierr);

        // Jacobian block: pressure-displacement (off-diagonal coupling)
        // The Biot coupling -alpha*u_t appears in f1_pressure (test gradient),
        // so its Jacobian goes in the g2 slot (test_gradient, trial_basis)
        ierr = PetscDSSetJacobian(prob, 0, 1,
                                  nullptr, nullptr,
                                  PetscFEPoroelasticity::g2_pu, nullptr); CHKERRQ(ierr);

        // Jacobian block: displacement-pressure (off-diagonal coupling)
        ierr = PetscDSSetJacobian(prob, 1, 0,
                                  nullptr, nullptr,
                                  PetscFEPoroelasticity::g2_up, nullptr); CHKERRQ(ierr);

        // Jacobian block: displacement-displacement (diagonal)
        ierr = PetscDSSetJacobian(prob, 1, 1,
                                  nullptr, nullptr, nullptr,
                                  PetscFEPoroelasticity::g3_uu); CHKERRQ(ierr);

        use_fem_time_residual_ = true;
    }
    // Register fluid flow residuals and jacobians (uncoupled or non-poroelastic)
    else if (config.fluid_model == FluidModelType::SINGLE_COMPONENT) {
        ierr = PetscDSSetResidual(prob, 0, PetscFEFluidFlow::f0_SinglePhase, PetscFEFluidFlow::f1_SinglePhase); CHKERRQ(ierr);
        // Dedicated single-phase Jacobian: g0 = phi*ct*shift, g3 = (k/mu)*I
        // Uses only single-field data (no Sw/Sg access like the BlackOil callbacks).
        ierr = PetscDSSetJacobian(prob, 0, 0, PetscFEFluidFlow::g0_SinglePhasePressure, nullptr, nullptr, PetscFEFluidFlow::g3_SinglePhasePressure); CHKERRQ(ierr);
        use_fem_time_residual_ = true;
    } else if (config.fluid_model == FluidModelType::BLACK_OIL) {
        // Pressure equation
        ierr = PetscDSSetResidual(prob, 0, PetscFEFluidFlow::f0_BlackOilPressure, PetscFEFluidFlow::f1_BlackOilPressure); CHKERRQ(ierr);
        ierr = PetscDSSetJacobian(prob, 0, 0, PetscFEFluidFlow::g0_BlackOilPressurePressure, nullptr, nullptr, PetscFEFluidFlow::g3_BlackOilPressurePressure); CHKERRQ(ierr);

        // Safe placeholders for saturation fields (keep constant)
        ierr = PetscDSSetResidual(prob, 1, PetscFEFluidFlow::f0_BlackOilSw, PetscFEFluidFlow::f1_Zero); CHKERRQ(ierr);
        ierr = PetscDSSetJacobian(prob, 1, 1, PetscFEFluidFlow::g0_Identity, nullptr, nullptr, nullptr); CHKERRQ(ierr);

        ierr = PetscDSSetResidual(prob, 2, PetscFEFluidFlow::f0_BlackOilSg, PetscFEFluidFlow::f1_Zero); CHKERRQ(ierr);
        ierr = PetscDSSetJacobian(prob, 2, 2, PetscFEFluidFlow::g0_Identity, nullptr, nullptr, nullptr); CHKERRQ(ierr);

        use_fem_time_residual_ = true;
    }

    // Register PETSc FEM callbacks for elasticity if enabled (skip if poroelastic, already handled above)
    if ((config.enable_geomechanics || config.enable_elastodynamics) &&
        !(config.solid_model == SolidModelType::POROELASTIC &&
          config.fluid_model == FluidModelType::SINGLE_COMPONENT)) {
        // Determine field index for displacement
        // Field ordering: [fluid fields...] [displacement] [thermal] [...]
        PetscInt displacement_field_idx = 0;
        if (config.fluid_model == FluidModelType::SINGLE_COMPONENT) {
            displacement_field_idx = 1;
        } else if (config.fluid_model == FluidModelType::BLACK_OIL) {
            displacement_field_idx = 3;
        } else if (config.fluid_model == FluidModelType::COMPOSITIONAL) {
            displacement_field_idx = 4;  // P + 3 compositions
        }
        // else NONE: displacement_field_idx = 0

        // Register callbacks
        if (config.enable_elastodynamics) {
            // Full dynamics with inertia
            // For TSALPHA2, we use the standard first-order form callbacks
            // The time stepper handles the second-order time integration internally
            if (use_aux_callbacks) {
                // Use auxiliary-field-aware callbacks (heterogeneous material)
                // f1: viscoelastic > elastoplastic > elastic (priority order)
                // f1: elastoplastic if enabled, otherwise standard elastic aux.
                // Viscoelastic uses the same elastic f1 (memory variables contribute
                // zero at initial relaxed state; future: separate memory variable update).
                auto f1_cb = config.enable_elastoplasticity
                    ? PetscFEElastoplasticity::f1_elastoplastic_aux
                    : PetscFEElasticityAux::f1_elastostatics_aux;
                auto g3_cb = config.enable_elastoplasticity
                    ? PetscFEElastoplasticity::g3_elastoplastic_aux
                    : PetscFEElasticityAux::g3_elastostatics_aux;
                ierr = PetscDSSetResidual(prob, displacement_field_idx,
                                          PetscFEElasticityAux::f0_elastodynamics_aux,
                                          f1_cb); CHKERRQ(ierr);
                ierr = PetscDSSetJacobian(prob, displacement_field_idx, displacement_field_idx,
                                          PetscFEElasticityAux::g0_elastodynamics_aux, nullptr, nullptr,
                                          g3_cb); CHKERRQ(ierr);
            } else {
                ierr = PetscDSSetResidual(prob, displacement_field_idx,
                                          PetscFEElasticity::f0_elastodynamics,
                                          PetscFEElasticity::f1_elastodynamics); CHKERRQ(ierr);
                ierr = PetscDSSetJacobian(prob, displacement_field_idx, displacement_field_idx,
                                          PetscFEElasticity::g0_elastodynamics, nullptr, nullptr,
                                          PetscFEElasticity::g3_elastodynamics); CHKERRQ(ierr);
            }
        } else {
            // Quasi-static elasticity (no inertia)
            if (use_aux_callbacks) {
                // Use auxiliary-field-aware callbacks (heterogeneous material or gravity)
                // f1: viscoelastic > elastoplastic > elastic (priority order)
                // f1: same as above (no custom viscoelastic callback for now)
                auto f1_cb = config.enable_elastoplasticity
                    ? PetscFEElastoplasticity::f1_elastoplastic_aux
                    : PetscFEElasticityAux::f1_elastostatics_aux;
                auto g3_cb = config.enable_elastoplasticity
                        ? PetscFEElastoplasticity::g3_elastoplastic_aux
                        : PetscFEElasticityAux::g3_elastostatics_aux;
                ierr = PetscDSSetResidual(prob, displacement_field_idx,
                                          PetscFEElasticityAux::f0_elastostatics_aux,
                                          f1_cb); CHKERRQ(ierr);
                ierr = PetscDSSetJacobian(prob, displacement_field_idx, displacement_field_idx,
                                          nullptr, nullptr, nullptr,
                                          g3_cb); CHKERRQ(ierr);
            } else {
                ierr = PetscDSSetResidual(prob, displacement_field_idx,
                                          PetscFEElasticity::f0_elastostatics,
                                          PetscFEElasticity::f1_elastostatics); CHKERRQ(ierr);
                ierr = PetscDSSetJacobian(prob, displacement_field_idx, displacement_field_idx,
                                          nullptr, nullptr, nullptr,
                                          PetscFEElasticity::g3_elastostatics); CHKERRQ(ierr);
            }
        }

        use_fem_time_residual_ = true;
    }

    // Register thermal PetscDS callbacks if enabled
    if (config.enable_thermal) {
        // Compute thermal field index: it is the last field added to the DM.
        // Field ordering: [fluid fields] [displacement] [lagrange] [temperature]
        PetscInt thermal_field_idx = 0;
        if (config.fluid_model == FluidModelType::SINGLE_COMPONENT) {
            thermal_field_idx += 1;
        } else if (config.fluid_model == FluidModelType::BLACK_OIL) {
            thermal_field_idx += 3;
        } else if (config.fluid_model == FluidModelType::COMPOSITIONAL) {
            thermal_field_idx += 4;
        }
        if (config.enable_geomechanics || config.enable_elastodynamics) {
            thermal_field_idx += 1; // displacement field
        }
        if (config.enable_faults && cohesive_kernel_) {
            thermal_field_idx += 1; // lagrange field
        }

        // Store thermal field index and material properties in constants
        const MaterialProperties& mat = material_props.empty() ? MaterialProperties() : material_props[0];
        unified_constants[PetscFEThermal::THERMAL_CONST_KAPPA] = mat.thermal_conductivity;
        unified_constants[PetscFEThermal::THERMAL_CONST_CP] = mat.heat_capacity;
        unified_constants[PetscFEThermal::THERMAL_CONST_ALPHA_T] = mat.thermal_expansion;
        unified_constants[PetscFEThermal::THERMAL_CONST_T_REF] = config.reference_temperature;
        unified_constants[PetscFEThermal::THERMAL_CONST_FIELD_IDX] = static_cast<PetscScalar>(thermal_field_idx);

        // Register thermal diffusion callbacks
        ierr = PetscDSSetResidual(prob, thermal_field_idx,
                                  PetscFEThermal::f0_thermal,
                                  PetscFEThermal::f1_thermal); CHKERRQ(ierr);
        ierr = PetscDSSetJacobian(prob, thermal_field_idx, thermal_field_idx,
                                  PetscFEThermal::g0_thermal_thermal, nullptr, nullptr,
                                  PetscFEThermal::g3_thermal_thermal); CHKERRQ(ierr);

        // If geomechanics/elastodynamics is also enabled, override the elastic
        // stress callback with the thermoelastic version that subtracts thermal strain.
        // The elastic Jacobian (g3) is unchanged since thermal strain is independent
        // of displacement gradient.
        if (config.enable_geomechanics || config.enable_elastodynamics) {
            PetscInt disp_idx = 0;
            if (config.fluid_model == FluidModelType::SINGLE_COMPONENT) disp_idx = 1;
            else if (config.fluid_model == FluidModelType::BLACK_OIL) disp_idx = 3;
            else if (config.fluid_model == FluidModelType::COMPOSITIONAL) disp_idx = 4;

            if (use_aux_callbacks) {
                // Override f1 with thermoelastic aux version
                // f0 remains the same (elastostatics/elastodynamics)
                // g3 remains the same (elastic tangent does not change)
                if (config.enable_elastodynamics) {
                    ierr = PetscDSSetResidual(prob, disp_idx,
                                              PetscFEElasticityAux::f0_elastodynamics_aux,
                                              PetscFEThermal::f1_thermoelastic_aux); CHKERRQ(ierr);
                } else {
                    ierr = PetscDSSetResidual(prob, disp_idx,
                                              PetscFEElasticityAux::f0_elastostatics_aux,
                                              PetscFEThermal::f1_thermoelastic_aux); CHKERRQ(ierr);
                }
            } else {
                if (config.enable_elastodynamics) {
                    ierr = PetscDSSetResidual(prob, disp_idx,
                                              PetscFEElasticity::f0_elastodynamics,
                                              PetscFEThermal::f1_thermoelastic); CHKERRQ(ierr);
                } else {
                    ierr = PetscDSSetResidual(prob, disp_idx,
                                              PetscFEElasticity::f0_elastostatics,
                                              PetscFEThermal::f1_thermoelastic); CHKERRQ(ierr);
                }
            }
        }

        use_fem_time_residual_ = true;
        if (rank == 0) {
            PetscPrintf(comm, "Registered thermal PetscDS callbacks: field %d, "
                        "kappa=%.2f, Cp=%.1f, alpha_T=%.2e, T_ref=%.1f\n",
                        (int)thermal_field_idx,
                        mat.thermal_conductivity, mat.heat_capacity,
                        mat.thermal_expansion, config.reference_temperature);
        }
    }

    // Determine field indices used by cohesive and absorbing boundary callbacks.
    PetscInt cohesive_disp_field = -1;
    PetscInt cohesive_lagrange_field = -1;
    PetscInt absorbing_disp_field = -1;
    if (config.fluid_model == FluidModelType::SINGLE_COMPONENT) {
        cohesive_disp_field = 1;
        absorbing_disp_field = 1;
    } else if (config.fluid_model == FluidModelType::BLACK_OIL) {
        cohesive_disp_field = 3;
        absorbing_disp_field = 3;
    } else if (config.fluid_model == FluidModelType::COMPOSITIONAL) {
        cohesive_disp_field = 4;
        absorbing_disp_field = 4;
    } else {
        cohesive_disp_field = 0;
        absorbing_disp_field = 0;
    }
    cohesive_lagrange_field = cohesive_disp_field + 1;

    const bool use_region_ds_split =
        config.enable_faults && cohesive_kernel_ && config.absorbing_bc_enabled &&
        (config.enable_geomechanics || config.enable_elastodynamics);

    // Register cohesive fault callbacks if enabled and region split is not needed.
    if (config.enable_faults && cohesive_kernel_ && !use_region_ds_split) {
        // Determine field indices
        // Field ordering: [fluid fields...] [displacement] [lagrange] [thermal] [...]
        PetscInt disp_field = cohesive_disp_field;
        PetscInt lagrange_field = cohesive_lagrange_field;

        // Lagrange field volume callback: zero residual.
        // Weak regularization: f = epsilon * lambda, g = epsilon * I with
        // epsilon = 1e-4. The small residual does not overwhelm BdResidual
        // on coarse meshes (unlike the old epsilon=1). The small Jacobian
        // diagonal is sufficient for interior Lagrange LU non-singularity.
        // On cohesive vertices, addCohesivePenaltyToJacobian adds a
        // penalty-scaled diagonal that dominates epsilon.
        // BdJacobian does NOT work in PETSc 3.25, so the Jacobian is
        // assembled manually in addCohesivePenaltyToJacobian.
        ierr = PetscDSSetResidual(prob, lagrange_field,
            f0_weak_lagrange, nullptr); CHKERRQ(ierr);
        ierr = PetscDSSetJacobian(prob, lagrange_field, lagrange_field,
            g0_weak_lagrange, nullptr,
            nullptr, nullptr); CHKERRQ(ierr);

        // Register BdResidual on cohesive cells. This works in PETSc 3.25+
        // (verified by Physics.CohesiveBdResidual test). DMPlexTSComputeIFunctionFEM
        // will evaluate these callbacks on the hybrid prism/tensor cells.
        if (hydrofrac_fem_pressurized_mode_) {
            ierr = PetscDSSetBdResidual(prob, lagrange_field,
                PetscFEHydrofrac::f0_lagrange_pressure_balance, nullptr); CHKERRQ(ierr);
        } else if (fault_mode_ == "prescribed_slip") {
            ierr = PetscDSSetBdResidual(prob, lagrange_field,
                CohesiveFaultKernel::f0_prescribed_slip, nullptr); CHKERRQ(ierr);
        } else {
            ierr = PetscDSSetBdResidual(prob, lagrange_field,
                CohesiveFaultKernel::f0_lagrange_constraint, nullptr); CHKERRQ(ierr);
        }
        ierr = PetscDSSetBdResidual(prob, disp_field,
            CohesiveFaultKernel::f0_displacement_cohesive, nullptr); CHKERRQ(ierr);

        if (rank == 0) {
            PetscPrintf(comm,
                        "Registered BdResidual cohesive callbacks on fields %d (disp), %d (lagrange)\n",
                        disp_field, lagrange_field);
        }
    }

    // Region-specific DS fix for cohesive + absorbing coexistence.
    if (use_region_ds_split) {
        DMLabel fault_label = nullptr;
        ierr = DMGetLabel(dm, "fault", &fault_label); CHKERRQ(ierr);
        if (!fault_label) {
            ierr = DMGetLabel(dm, "Fault", &fault_label); CHKERRQ(ierr);
        }
        PetscCheck(fault_label, comm, PETSC_ERR_USER,
                   "Fault + absorbing coexistence requested but fault label was not found on DM");

        DMLabel absorbing_label = nullptr;
        ierr = DMCreateLabel(dm, "absorbing_boundary"); CHKERRQ(ierr);
        ierr = DMGetLabel(dm, "absorbing_boundary", &absorbing_label); CHKERRQ(ierr);
        PetscCheck(absorbing_label, comm, PETSC_ERR_USER,
                   "Failed to create absorbing boundary label");

        auto addBoundaryLabelPoints = [&](const char *label_name, bool enabled) -> PetscErrorCode {
            if (!enabled) {
                return PETSC_SUCCESS;
            }
            DMLabel src = nullptr;
            PetscErrorCode ierr_local = DMGetLabel(dm, label_name, &src);CHKERRQ(ierr_local);
            if (!src) {
                return PETSC_SUCCESS;
            }
            IS stratum = nullptr;
            ierr_local = DMLabelGetStratumIS(src, 1, &stratum);CHKERRQ(ierr_local);
            if (!stratum) {
                return PETSC_SUCCESS;
            }
            const PetscInt *pts = nullptr;
            PetscInt npts = 0;
            ierr_local = ISGetLocalSize(stratum, &npts);CHKERRQ(ierr_local);
            ierr_local = ISGetIndices(stratum, &pts);CHKERRQ(ierr_local);
            for (PetscInt i = 0; i < npts; ++i) {
                ierr_local = DMLabelSetValue(absorbing_label, pts[i], 1);CHKERRQ(ierr_local);
            }
            ierr_local = ISRestoreIndices(stratum, &pts);CHKERRQ(ierr_local);
            ierr_local = ISDestroy(&stratum);CHKERRQ(ierr_local);
            return PETSC_SUCCESS;
        };

        ierr = addBoundaryLabelPoints("boundary_x_min", config.absorbing_bc_x_min); CHKERRQ(ierr);
        ierr = addBoundaryLabelPoints("boundary_x_max", config.absorbing_bc_x_max); CHKERRQ(ierr);
        ierr = addBoundaryLabelPoints("boundary_y_min", config.absorbing_bc_y_min); CHKERRQ(ierr);
        ierr = addBoundaryLabelPoints("boundary_y_max", config.absorbing_bc_y_max); CHKERRQ(ierr);
        ierr = addBoundaryLabelPoints("boundary_z_min", config.absorbing_bc_z_min); CHKERRQ(ierr);
        ierr = addBoundaryLabelPoints("boundary_z_max", config.absorbing_bc_z_max); CHKERRQ(ierr);

        PetscInt nfields = 0;
        ierr = PetscDSGetNumFields(prob, &nfields); CHKERRQ(ierr);

        PetscDS fault_ds = nullptr;
        PetscDS absorbing_ds = nullptr;
        ierr = PetscDSCreate(comm, &fault_ds); CHKERRQ(ierr);
        ierr = PetscDSCreate(comm, &absorbing_ds); CHKERRQ(ierr);
        ierr = PetscDSCopy(prob, 0, nfields, dm, fault_ds); CHKERRQ(ierr);
        ierr = PetscDSCopy(prob, 0, nfields, dm, absorbing_ds); CHKERRQ(ierr);

        if (hydrofrac_fem_pressurized_mode_) {
            // Hydrofrac mode: BdResidual still called on fault_ds because
            // addFaultPressureToResidual handles the real pressure contribution.
            // The BdResidual here is for the Lagrange multiplier coupling.
            ierr = PetscDSSetBdResidual(fault_ds, cohesive_disp_field,
                                        CohesiveFaultKernel::f0_displacement_cohesive,
                                        nullptr); CHKERRQ(ierr);
            ierr = PetscDSSetBdResidual(fault_ds, cohesive_lagrange_field,
                                        PetscFEHydrofrac::f0_lagrange_pressure_balance,
                                        nullptr); CHKERRQ(ierr);
            ierr = PetscDSSetBdJacobian(fault_ds, cohesive_disp_field, cohesive_lagrange_field,
                                        CohesiveFaultKernel::g0_displacement_lagrange, nullptr,
                                        nullptr, nullptr); CHKERRQ(ierr);
            ierr = PetscDSSetBdJacobian(fault_ds, cohesive_lagrange_field, cohesive_disp_field,
                                        CohesiveFaultKernel::g0_lagrange_displacement, nullptr,
                                        nullptr, nullptr); CHKERRQ(ierr);
            ierr = PetscDSSetBdJacobian(fault_ds, cohesive_lagrange_field, cohesive_lagrange_field,
                                        CohesiveFaultKernel::g0_lagrange_lagrange, nullptr,
                                        nullptr, nullptr); CHKERRQ(ierr);
        } else {
            // Non-hydrofrac fault modes: register BdResidual on fault_ds.
            // PETSc 3.25 supports BdResidual on cohesive cells.
            if (fault_mode_ == "prescribed_slip") {
                ierr = PetscDSSetBdResidual(fault_ds, cohesive_lagrange_field,
                                            CohesiveFaultKernel::f0_prescribed_slip,
                                            nullptr); CHKERRQ(ierr);
            } else {
                ierr = PetscDSSetBdResidual(fault_ds, cohesive_lagrange_field,
                                            CohesiveFaultKernel::f0_lagrange_constraint,
                                            nullptr); CHKERRQ(ierr);
            }
            ierr = PetscDSSetBdResidual(fault_ds, cohesive_disp_field,
                                        CohesiveFaultKernel::f0_displacement_cohesive,
                                        nullptr); CHKERRQ(ierr);
        }

        ierr = PetscDSSetBdResidual(absorbing_ds, absorbing_disp_field,
                                    AbsorbingBC::f0_absorbing, nullptr); CHKERRQ(ierr);
        ierr = PetscDSSetBdJacobian(absorbing_ds, absorbing_disp_field, absorbing_disp_field,
                                    AbsorbingBC::g0_absorbing, nullptr, nullptr, nullptr); CHKERRQ(ierr);

        ierr = DMSetRegionDS(dm, fault_label, nullptr, fault_ds, nullptr); CHKERRQ(ierr);
        ierr = DMSetRegionDS(dm, absorbing_label, nullptr, absorbing_ds, nullptr); CHKERRQ(ierr);

        ierr = PetscDSDestroy(&fault_ds); CHKERRQ(ierr);
        ierr = PetscDSDestroy(&absorbing_ds); CHKERRQ(ierr);

        if (rank == 0) {
            PetscPrintf(comm,
                        "Registered region DS split: cohesive callbacks on fault label, absorbing callbacks on external boundaries\n");
        }
    }

    if (config.absorbing_bc_enabled &&
        (config.enable_geomechanics || config.enable_elastodynamics) && !use_region_ds_split) {
        PetscInt displacement_field_idx = absorbing_disp_field;

        ierr = PetscDSSetBdResidual(prob, displacement_field_idx,
            AbsorbingBC::f0_absorbing, NULL); CHKERRQ(ierr);
        ierr = PetscDSSetBdJacobian(prob, displacement_field_idx, displacement_field_idx,
            AbsorbingBC::g0_absorbing, NULL, NULL, NULL); CHKERRQ(ierr);

        if (rank == 0) {
            PetscPrintf(comm, "Registered absorbing BC callbacks on field %d\n",
                        displacement_field_idx);
        }
    }

    // Register traction boundary residual callbacks
    // Must use PetscWeakFormAddBdResidual with the exact label and value that
    // DMAddBoundary used, because PETSc weak form lookup is a strict hash match
    // on (label, value, field, part) with no fallback to NULL label.
    if (has_traction_bc) {
        PetscInt disp_field = absorbing_disp_field;
        if (disp_field < 0) {
            disp_field = 0;  // Fallback: displacement is first field
        }

        PetscWeakForm wf = nullptr;
        ierr = PetscDSGetWeakForm(prob, &wf); CHKERRQ(ierr);

        DMLabel face_labels_phys[6];
        ierr = DMGetLabel(dm, "boundary_x_min", &face_labels_phys[0]); CHKERRQ(ierr);
        ierr = DMGetLabel(dm, "boundary_x_max", &face_labels_phys[1]); CHKERRQ(ierr);
        ierr = DMGetLabel(dm, "boundary_y_min", &face_labels_phys[2]); CHKERRQ(ierr);
        ierr = DMGetLabel(dm, "boundary_y_max", &face_labels_phys[3]); CHKERRQ(ierr);
        ierr = DMGetLabel(dm, "boundary_z_min", &face_labels_phys[4]); CHKERRQ(ierr);
        ierr = DMGetLabel(dm, "boundary_z_max", &face_labels_phys[5]); CHKERRQ(ierr);

        for (int fi = 0; fi < 6; ++fi) {
            if (!config.face_bc[fi].configured || config.face_bc[fi].type != "traction") continue;
            if (!face_labels_phys[fi]) continue;
            PetscInt label_value = 1;
            ierr = PetscWeakFormAddBdResidual(wf, face_labels_phys[fi], label_value,
                disp_field, 0, f0_traction_bc, nullptr); CHKERRQ(ierr);
            if (rank == 0) {
                static const char *face_names[6] = {
                    "x_min", "x_max", "y_min", "y_max", "z_min", "z_max"
                };
                PetscPrintf(comm, "  Registered traction BdResidual on %s, field %d\n",
                            face_names[fi], disp_field);
            }
        }
    }

    // Setup auxiliary DM for heterogeneous material properties or gravity
    if (use_aux_callbacks) {
        ierr = setupAuxiliaryDM(); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::setupAuxiliaryDM()
{
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    PetscInt dim;
    ierr = DMGetDimension(dm, &dim); CHKERRQ(ierr);

    // Determine if mesh is simplex
    PetscBool isSimplex = PETSC_FALSE;
    DMPolytopeType ct;
    PetscInt cStart, cEnd;
    ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd); CHKERRQ(ierr);
    if (cEnd > cStart) {
        ierr = DMPlexGetCellType(dm, cStart, &ct); CHKERRQ(ierr);
        if (ct == DM_POLYTOPE_TETRAHEDRON || ct == DM_POLYTOPE_TRIANGLE) {
            isSimplex = PETSC_TRUE;
        }
    }

    // Clone the primary DM topology
    ierr = DMClone(dm, &auxDM_); CHKERRQ(ierr);

    // Get quadrature from the primary FE to match tabulation points
    PetscQuadrature quad = NULL;
    if (!fe_fields.empty()) {
        ierr = PetscFEGetQuadrature(fe_fields[0], &quad); CHKERRQ(ierr);
    }

    // Create scalar FE fields for lambda, mu, rho on the auxiliary DM
    // Use degree-0 (constant per cell) for material properties
    for (PetscInt f = 0; f < NUM_AUX_FIELDS; ++f) {
        PetscFE fe_aux;
        ierr = PetscFECreateLagrange(PETSC_COMM_SELF, dim, 1, isSimplex, 0, PETSC_DETERMINE, &fe_aux); CHKERRQ(ierr);
        // Match quadrature to primary FE so tabulation point counts agree
        if (quad) {
            ierr = PetscFESetQuadrature(fe_aux, quad); CHKERRQ(ierr);
        }
        ierr = DMSetField(auxDM_, f, NULL, (PetscObject)fe_aux); CHKERRQ(ierr);
        ierr = PetscFEDestroy(&fe_aux); CHKERRQ(ierr);
    }
    // Memory variables for viscoelastic attenuation are NOT stored as aux DM
    // fields. The aux DM keeps only the 3 material property fields (lambda,
    // mu, rho). Memory variables (which start at zero for the relaxed initial
    // state) are updated externally by the TSPostStep callback using a
    // separate storage mechanism. The f1_viscoelastic_aux callback reads only
    // the relaxed moduli and adds zero for the memory variables at startup.
    //
    // This avoids PETSc DMPlex FEM tabulation issues with large numbers of
    // aux fields that can cause SEGV during pointwise callback evaluation.

    ierr = DMCreateDS(auxDM_); CHKERRQ(ierr);

    // Create local vector and zero-initialize (important for viscoelastic
    // memory variable fields that start at relaxed state = zero)
    ierr = DMCreateLocalVector(auxDM_, &auxVec_); CHKERRQ(ierr);
    ierr = VecZeroEntries(auxVec_); CHKERRQ(ierr);
    if (config.material_assignment == "gmsh_label") {
        ierr = populateAuxFieldsByMaterialLabel(); CHKERRQ(ierr);
    } else if (config.material_assignment == "velocity_model") {
        ierr = populateAuxFieldsByVelocityModel(); CHKERRQ(ierr);
    } else {
        ierr = populateAuxFieldsByDepth(); CHKERRQ(ierr);
    }
    ierr = applyExplosionDamageToAuxFields(); CHKERRQ(ierr);

    // Attach auxiliary data to primary DM
    ierr = DMSetAuxiliaryVec(dm, NULL, 0, 0, auxVec_); CHKERRQ(ierr);

    if (rank == 0) {
        PetscPrintf(comm, "Auxiliary DM setup complete: %d layers, %d aux fields\n",
                    (int)config.material_layers.size(), (int)NUM_AUX_FIELDS);
    }

    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::populateAuxFieldsByDepth()
{
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    PetscInt dim;
    ierr = DMGetDimension(dm, &dim); CHKERRQ(ierr);

    // Get cell range
    PetscInt cStart, cEnd;
    ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd); CHKERRQ(ierr);

    // Get the section for the auxiliary DM
    PetscSection section;
    ierr = DMGetLocalSection(auxDM_, &section); CHKERRQ(ierr);

    // Default material properties (from first material_props if available)
    double default_lambda = 30.0e9, default_mu = 25.0e9, default_rho = 2650.0;
    if (!material_props.empty()) {
        double E = material_props[0].youngs_modulus;
        double nu = material_props[0].poisson_ratio;
        default_lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        default_mu = E / (2.0 * (1.0 + nu));
        default_rho = material_props[0].density;
    }

    PetscScalar *a;
    ierr = VecGetArray(auxVec_, &a); CHKERRQ(ierr);

    for (PetscInt c = cStart; c < cEnd; ++c) {
        // Compute cell centroid via FVM geometry
        PetscReal vol, centroid[3], normal[3];
        ierr = DMPlexComputeCellGeometryFVM(dm, c, &vol, centroid, normal); CHKERRQ(ierr);

        // Determine z-coordinate (last dim)
        PetscReal z = centroid[dim - 1];

        // Look up material properties by depth
        PetscScalar lambda = default_lambda;
        PetscScalar mu = default_mu;
        PetscScalar rho = default_rho;

        for (const auto& layer : config.material_layers) {
            if (z <= layer.z_top && z > layer.z_bottom) {
                lambda = layer.lambda;
                mu = layer.mu;
                rho = layer.rho;
                break;
            }
        }

        // Get offset for this cell in the auxiliary vector
        PetscInt offset;
        ierr = PetscSectionGetOffset(section, c, &offset); CHKERRQ(ierr);

        a[offset + AUX_LAMBDA] = lambda;
        a[offset + AUX_MU]     = mu;
        a[offset + AUX_RHO]    = rho;
    }

    ierr = VecRestoreArray(auxVec_, &a); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::populateAuxFieldsByMaterialLabel()
{
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

    DMLabel matLabel = nullptr;
    ierr = DMGetLabel(dm, "Material", &matLabel); CHKERRQ(ierr);
    PetscCheck(matLabel, comm, PETSC_ERR_ARG_WRONG,
        "material assignment=gmsh_label but no Material label exists on the DM");

    PetscInt cStart = 0, cEnd = 0;
    ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd); CHKERRQ(ierr);

    PetscSection section;
    ierr = DMGetLocalSection(auxDM_, &section); CHKERRQ(ierr);

    double default_lambda = 30.0e9, default_mu = 25.0e9, default_rho = 2650.0;
    if (!material_props.empty()) {
        const double youngs_modulus = material_props[0].youngs_modulus;
        const double poissons_ratio = material_props[0].poisson_ratio;
        default_lambda = youngs_modulus * poissons_ratio /
                         ((1.0 + poissons_ratio) * (1.0 - 2.0 * poissons_ratio));
        default_mu = youngs_modulus / (2.0 * (1.0 + poissons_ratio));
        default_rho = material_props[0].density;
    }

    PetscScalar* a = nullptr;
    ierr = VecGetArray(auxVec_, &a); CHKERRQ(ierr);

    for (PetscInt c = cStart; c < cEnd; ++c) {
        PetscScalar lambda = default_lambda;
        PetscScalar mu = default_mu;
        PetscScalar rho = default_rho;

        PetscInt mat_id = -1;
        ierr = DMLabelGetValue(matLabel, c, &mat_id); CHKERRQ(ierr);
        if (mat_id >= 0) {
            for (const auto& region : config.material_regions) {
                if (region.label_id == mat_id) {
                    lambda = region.lambda;
                    mu = region.mu;
                    rho = region.rho;
                    break;
                }
            }
        }

        PetscInt offset = 0;
        ierr = PetscSectionGetOffset(section, c, &offset); CHKERRQ(ierr);
        a[offset + AUX_LAMBDA] = lambda;
        a[offset + AUX_MU] = mu;
        a[offset + AUX_RHO] = rho;
    }

    ierr = VecRestoreArray(auxVec_, &a); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::populateAuxFieldsByVelocityModel()
{
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

    VelocityModel model;
    int read_err = readVelocityModel(config.velocity_model_file, model);
    PetscCheck(read_err == 0, comm, PETSC_ERR_FILE_OPEN,
        "Failed to read velocity model file: %s (error code %d)",
        config.velocity_model_file.c_str(), read_err);

    PetscInt dim = 0;
    ierr = DMGetDimension(dm, &dim); CHKERRQ(ierr);

    PetscInt cStart = 0, cEnd = 0;
    ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd); CHKERRQ(ierr);

    PetscSection section;
    ierr = DMGetLocalSection(auxDM_, &section); CHKERRQ(ierr);

    PetscScalar *a = nullptr;
    ierr = VecGetArray(auxVec_, &a); CHKERRQ(ierr);

    for (PetscInt c = cStart; c < cEnd; ++c) {
        PetscReal vol, centroid[3], normal[3];
        ierr = DMPlexComputeCellGeometryFVM(dm, c, &vol, centroid, normal); CHKERRQ(ierr);

        double vp_val, vs_val, rho_val;
        model.interpolate(centroid[0], centroid[1], centroid[2],
                          vp_val, vs_val, rho_val);

        double mu_val = rho_val * vs_val * vs_val;
        double lambda_val = rho_val * (vp_val * vp_val - 2.0 * vs_val * vs_val);

        PetscInt offset = 0;
        ierr = PetscSectionGetOffset(section, c, &offset); CHKERRQ(ierr);

        a[offset + AUX_LAMBDA] = lambda_val;
        a[offset + AUX_MU]     = mu_val;
        a[offset + AUX_RHO]    = rho_val;
    }

    ierr = VecRestoreArray(auxVec_, &a); CHKERRQ(ierr);

    if (rank == 0) {
        PetscPrintf(comm, "Populated auxiliary fields from velocity model: %s "
                    "(%d x %d x %d grid)\n",
                    config.velocity_model_file.c_str(),
                    (int)model.nx, (int)model.ny, (int)model.nz);
    }

    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::applyExplosionDamageToAuxFields()
{
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

    if (!config.enable_near_field_damage || !auxVec_ || !explosion_) {
        PetscFunctionReturn(0);
    }

    PetscSection section;
    ierr = DMGetLocalSection(auxDM_, &section); CHKERRQ(ierr);

    PetscInt dim = 0;
    ierr = DMGetDimension(dm, &dim); CHKERRQ(ierr);

    NuclearSourceParameters source_params;
    source_params.yield_kt = explosion_->yield_kt;
    source_params.depth_of_burial = explosion_->depth_of_burial;

    // Use 1D solver zone radii when COUPLED_ANALYTIC, else empirical scaling
    double cavity_radius, crushed_radius, fractured_radius;
    if (explosion_->use_nearfield_coupling) {
        cavity_radius = explosion_->nf_cavity_radius;
        crushed_radius = explosion_->nf_crushed_radius;
        fractured_radius = explosion_->nf_fractured_radius;
    } else {
        cavity_radius = source_params.cavity_radius(explosion_->rho);
        crushed_radius = source_params.crushed_zone_radius();
        fractured_radius = source_params.fractured_zone_radius();
    }

    PetscInt cStart = 0, cEnd = 0;
    ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd); CHKERRQ(ierr);

    PetscScalar* a = nullptr;
    ierr = VecGetArray(auxVec_, &a); CHKERRQ(ierr);

    for (PetscInt c = cStart; c < cEnd; ++c) {
        PetscReal volume = 0.0, centroid[3] = {0.0, 0.0, 0.0}, normal[3] = {0.0, 0.0, 0.0};
        ierr = DMPlexComputeCellGeometryFVM(dm, c, &volume, centroid, normal); CHKERRQ(ierr);

        const double dx = centroid[0] - explosion_->sx;
        const double dy = centroid[1] - explosion_->sy;
        const double dz = ((dim > 2) ? centroid[2] : 0.0) - explosion_->sz;
        const double radius = std::sqrt(dx * dx + dy * dy + dz * dz);

        PetscInt offset = 0;
        ierr = PetscSectionGetOffset(section, c, &offset); CHKERRQ(ierr);

        PetscScalar& lambda = a[offset + AUX_LAMBDA];
        PetscScalar& mu = a[offset + AUX_MU];
        PetscScalar& rho = a[offset + AUX_RHO];

        const PetscScalar base_lambda = lambda;
        const PetscScalar base_mu = mu;

        if (explosion_->use_nearfield_coupling) {
            // COUPLED_ANALYTIC: use 1D solver continuous damage profile
            std::array<double, 3> pt = {centroid[0], centroid[1],
                                        (dim > 2) ? centroid[2] : 0.0};
            double D = explosion_->nf_solver.getDamage(pt);
            if (D > 0.9) {
                // Cavity/crushed zone
                lambda = 1.0e6;
                mu = 1.0e6;
                rho = 100.0;
            } else if (D > 0.01) {
                double factor = 1.0 - D;
                factor = std::max(0.1, factor);
                lambda = factor * base_lambda;
                mu = factor * base_mu;
            }
        } else {
            // PROXY: simple 3-zone model
            if (radius < cavity_radius) {
                lambda = 1.0e6;
                mu = 1.0e6;
                rho = 100.0;
            } else if (radius < crushed_radius) {
                lambda = 0.1 * base_lambda;
                mu = 0.1 * base_mu;
            } else if (radius < fractured_radius) {
                const double frac = (radius - crushed_radius) /
                                    std::max(1.0e-12, fractured_radius - crushed_radius);
                const double factor = 0.5 + 0.5 * std::clamp(frac, 0.0, 1.0);
                lambda = factor * base_lambda;
                mu = factor * base_mu;
            }
        }
    }

    ierr = VecRestoreArray(auxVec_, &a); CHKERRQ(ierr);
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

    // Create Jacobian matrix (required by TSSetIJacobian callbacks)
    if (!jacobian) {
        ierr = DMCreateMatrix(dm, &jacobian); CHKERRQ(ierr);
        ierr = MatZeroEntries(jacobian); CHKERRQ(ierr);
    }
    
    // Set time stepping method
    // Use TSALPHA2 (generalized-alpha) for elastodynamics: second-order time integrator
    // TSALPHA2 handles the second-order ODE M*u_tt + K*u = f using generalized-alpha method.
    // With IFunction form, u_t represents velocity and TSALPHA2 internally tracks acceleration.
    // Use TSBEULER for other cases (single-phase flow, poroelasticity): first-order, more stable
    if (config.enable_elastodynamics) {
        ierr = TSSetType(ts, TSALPHA2); CHKERRQ(ierr);
        // Set spectral radius for numerical dissipation (1.0 = no dissipation, 0.0 = maximum)
        ierr = TSAlpha2SetRadius(ts, 1.0); CHKERRQ(ierr);
        if (rank == 0) {
            PetscPrintf(comm, "Using TSALPHA2 for elastodynamics (generalized-alpha, second-order)\n");
        }
    } else {
        ierr = TSSetType(ts, TSBEULER); CHKERRQ(ierr);
    }

    // Set time parameters
    ierr = TSSetTime(ts, config.start_time); CHKERRQ(ierr);
    ierr = TSSetMaxTime(ts, config.end_time); CHKERRQ(ierr);
    ierr = TSSetTimeStep(ts, config.dt_initial); CHKERRQ(ierr);
    ierr = TSSetMaxSteps(ts, config.max_timesteps); CHKERRQ(ierr);
    ierr = TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP); CHKERRQ(ierr);

    // Set residual and Jacobian functions
    // Use standard IFunction/IJacobian form for all cases
    // TSALPHA2 handles second-order time integration internally using first-order form
    ierr = TSSetIFunction(ts, nullptr, FormFunction, this); CHKERRQ(ierr);
    ierr = TSSetIJacobian(ts, jacobian, jacobian, FormJacobian, this); CHKERRQ(ierr);
    
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
    
    // Register viscoelastic memory variable update as a post-step callback.
    // Currently disabled: the TSPostStep is a no-op placeholder. Memory variable
    // updates will be re-enabled when a dedicated storage mechanism is implemented.
    if (config.enable_viscoelastic && rank == 0) {
        PetscPrintf(comm, "Viscoelastic attenuation enabled (relaxed moduli mode, no memory variable update)\n");
    }

    // Set from options
    ierr = TSSetFromOptions(ts); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

// =============================================================================
// ViscoelasticPostStep
//
// Update memory variables R_i after each converged time step using the
// exponentially-fitted update formula:
//   R_i^{n+1} = R_i^n * exp(-dt/tau_i) + 2 * delta_mu_i * mu * deps * phi_i
// where phi_i = tau_i/dt * (1 - exp(-dt/tau_i)) and deps is the strain increment.
//
// This callback reads the current displacement gradient to compute the strain
// increment relative to the previous step (stored in solution_old).
// =============================================================================
PetscErrorCode Simulator::ViscoelasticPostStep(TS ts)
{
    PetscFunctionBeginUser;
    void *ctx;
    PetscErrorCode ierr;
    ierr = TSGetApplicationContext(ts, &ctx); CHKERRQ(ierr);
    Simulator *sim = static_cast<Simulator*>(ctx);

    if (!sim->config.enable_viscoelastic || !sim->auxDM_ || !sim->auxVec_)
        PetscFunctionReturn(PETSC_SUCCESS);

    PetscReal dt;
    ierr = TSGetTimeStep(ts, &dt); CHKERRQ(ierr);
    if (dt <= 0.0) PetscFunctionReturn(PETSC_SUCCESS);

    const int N = sim->config.visco_num_mechanisms;
    if (N <= 0) PetscFunctionReturn(PETSC_SUCCESS);

    // Read constants for relaxation times and weights
    PetscDS prob;
    ierr = DMGetDS(sim->dm, &prob); CHKERRQ(ierr);
    PetscInt nconst = 0;
    const PetscScalar *constants = nullptr;
    ierr = PetscDSGetConstants(prob, &nconst, &constants); CHKERRQ(ierr);

    double tau_arr[5] = {0}, dmu_arr[5] = {0};
    for (int m = 0; m < N && m < 5; ++m) {
        tau_arr[m] = (nconst > PetscFEViscoelastic::VISCO_CONST_TAU_BASE + m)
            ? PetscRealPart(constants[PetscFEViscoelastic::VISCO_CONST_TAU_BASE + m]) : 1.0;
        dmu_arr[m] = (nconst > PetscFEViscoelastic::VISCO_CONST_DMU_BASE + m)
            ? PetscRealPart(constants[PetscFEViscoelastic::VISCO_CONST_DMU_BASE + m]) : 0.0;
    }

    // Memory variables are not stored in the aux DM to avoid tabulation
    // issues. The relaxed-modulus-only approach currently used by the
    // f1_viscoelastic_aux callback does not require post-step updates.
    //
    // Future: when a separate memory variable storage mechanism is
    // implemented (e.g., a dedicated Vec), the exponential decay update
    //   R_i^{n+1} = R_i^n * exp(-dt/tau_i) + delta_mu_i * deps * phi_i
    // would be applied here.
    (void)tau_arr;
    (void)dmu_arr;

    PetscFunctionReturn(PETSC_SUCCESS);
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

    // For elastodynamics with TSALPHA2, the system matrix is K + shift*M.
    // Use direct LU solver for robustness. For larger problems, can be overridden
    // with command-line options: -ksp_type cg -pc_type gamg
    if (config.enable_elastodynamics) {
        ierr = KSPSetType(ksp, KSPPREONLY); CHKERRQ(ierr);
        ierr = PCSetType(pc, PCLU); CHKERRQ(ierr);
    }
    
    // Set solver options from command line (override programmatic defaults)
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
    if (!config.enable_faults) PetscFunctionReturn(0);
    PetscErrorCode ierr;

    if (rank == 0) PetscPrintf(comm, "Setting up fault network...\n");

    // 1. Create fault mesh manager
    fault_mesh_manager_ = std::make_unique<FaultMeshManager>(comm);

    // 2. Read fault geometry from config
    //    ConfigReader should have parsed [FAULT] section with:
    //    strike, dip, center_x/y/z, length, width, tolerance
    //    For now, use defaults if not in config
    double strike = 0.0;        // radians
    double dip = M_PI / 2.0;    // vertical fault
    double center[3] = {0.0, 0.0, 0.0};
    double length = 1e10;       // large enough to cut the whole mesh
    double width = 1e10;
    // Simplex mesh face centroids are not at exact grid coords, use smaller tolerance
    double tol = 0.05;

    if (fracture_plane_enabled_) {
        // Use fracture plane geometry from [FRACTURE_PLANE] config
        strike = fracture_plane_strike_;
        dip = fracture_plane_dip_;
        center[0] = fracture_plane_center_[0];
        center[1] = fracture_plane_center_[1];
        center[2] = fracture_plane_center_[2];
        length = fracture_plane_length_;
        width = fracture_plane_width_;
    } else if (fault_geometry_from_config_) {
        // Use fault geometry from [FAULT] config section
        strike = fault_strike_;
        dip    = fault_dip_;
        // Replace any -1 sentinel with grid-center defaults
        center[0] = (fault_center_[0] >= 0.0) ? fault_center_[0] : grid_config.Lx / 2.0;
        center[1] = (fault_center_[1] >= 0.0) ? fault_center_[1] : grid_config.Ly / 2.0;
        center[2] = (fault_center_[2] >= 0.0) ? fault_center_[2] : grid_config.Lz / 2.0;
        double max_dim = std::max({grid_config.Lx, grid_config.Ly, grid_config.Lz});
        length = (fault_length_ > 0.0) ? fault_length_ : 2.0 * max_dim;
        width  = (fault_width_  > 0.0) ? fault_width_  : 2.0 * max_dim;
    } else {
        // Default: vertical fault at mesh center cutting entire domain
        center[0] = grid_config.Lx / 2.0;
        center[1] = grid_config.Ly / 2.0;
        center[2] = grid_config.Lz / 2.0;
        length = 2.0 * std::max({grid_config.Lx, grid_config.Ly, grid_config.Lz});
        width = length;
    }

    // 3. Label fault faces on the DM
    DMLabel fault_label = nullptr;
    ierr = fault_mesh_manager_->createPlanarFaultLabel(
        dm, &fault_label, strike, dip, center, length, width, tol);
    if (ierr != 0) {
        if (rank == 0) PetscPrintf(comm, "Warning: fault labeling failed, skipping fault setup\n");
        fault_mesh_manager_.reset();
        PetscFunctionReturn(0);
    }

    // 4. Split mesh along fault (inserts cohesive cells)
    //    THIS CHANGES THE DM TOPOLOGY
    ierr = fault_mesh_manager_->splitMeshAlongFault(&dm, "fault");
    if (ierr != 0) {
        if (rank == 0) PetscPrintf(comm, "Warning: mesh splitting failed, skipping fault setup\n");
        fault_mesh_manager_.reset();
        PetscFunctionReturn(0);
    }

    // 5. Create FaultCohesiveDyn and extract topology from cohesive cells
    cohesive_fault_ = std::make_unique<FaultCohesiveDyn>();
    ierr = fault_mesh_manager_->extractCohesiveTopology(dm, cohesive_fault_.get());
    CHKERRQ(ierr);

    // 6. Set friction model
    //    TODO: Read from config. Default: slip-weakening
    auto friction = std::make_unique<SlipWeakeningFriction>();
    friction->setStaticCoefficient(0.6);
    friction->setDynamicCoefficient(0.4);
    friction->setCriticalSlipDistance(0.4);
    cohesive_fault_->setFrictionModel(std::move(friction));

    // 7. Set initial traction (pre-stress on fault)
    //    For locked fault test: doesn't matter, but set something physical
    cohesive_fault_->setUniformInitialTraction(25e6, 0.0, -50e6);

    // 8. Initialize fault state
    cohesive_fault_->initialize();

    // 9. Create cohesive kernel
    cohesive_kernel_ = std::make_unique<CohesiveFaultKernel>();
    
    if (fault_mode_ == "prescribed_slip") {
        // Convert from fault-local (strike, dip, opening) to Cartesian (x, y, z)
        double fault_strike = fault_geometry_from_config_
            ? fault_strike_
            : (fracture_plane_enabled_ ? fracture_plane_strike_ : 0.0);
        double fault_dip = fault_geometry_from_config_
            ? fault_dip_
            : (fracture_plane_enabled_ ? fracture_plane_dip_ : M_PI / 2.0);
        
        double cs = std::cos(fault_strike);
        double ss = std::sin(fault_strike);
        double cd = std::cos(fault_dip);
        double sd = std::sin(fault_dip);
        
        // Strike direction (horizontal, along-fault)
        double strike_dir[3] = {cs, ss, 0.0};
        // Normal direction (into hanging wall)
        double normal_dir[3] = {-ss * sd, cs * sd, -cd};
        // Dip direction = normal x strike (downdip)
        double dip_dir[3] = {
            normal_dir[1] * strike_dir[2] - normal_dir[2] * strike_dir[1],
            normal_dir[2] * strike_dir[0] - normal_dir[0] * strike_dir[2],
            normal_dir[0] * strike_dir[1] - normal_dir[1] * strike_dir[0]
        };
        
        double slip_x = fault_slip_strike_ * strike_dir[0] + fault_slip_dip_ * dip_dir[0] + fault_slip_opening_ * normal_dir[0];
        double slip_y = fault_slip_strike_ * strike_dir[1] + fault_slip_dip_ * dip_dir[1] + fault_slip_opening_ * normal_dir[1];
        double slip_z = fault_slip_strike_ * strike_dir[2] + fault_slip_dip_ * dip_dir[2] + fault_slip_opening_ * normal_dir[2];
        
        cohesive_kernel_->setPrescribedSlip(slip_x, slip_y, slip_z);
        cohesive_kernel_->setMode(true);
        
        if (rank == 0) {
            PetscPrintf(comm, "Prescribed slip mode: strike=%.3f, dip=%.3f, opening=%.3f\n",
                        fault_slip_strike_, fault_slip_dip_, fault_slip_opening_);
            PetscPrintf(comm, "Prescribed slip (Cartesian): (%.6f, %.6f, %.6f)\n",
                        slip_x, slip_y, slip_z);
        }
    } else if (fault_mode_ == "slipping") {
        cohesive_kernel_->setMode(false);
        cohesive_kernel_->setFrictionCoefficient(fault_friction_coefficient_);
    } else {
        cohesive_kernel_->setMode(true);  // locked
        cohesive_kernel_->setFrictionCoefficient(fault_friction_coefficient_);
    }
    
    if (fracture_plane_enabled_) {
        cohesive_kernel_->setTensileStrength(fracture_plane_tensile_strength_);
    }

    if (rank == 0) {
        PetscPrintf(comm, "Fault network setup complete: %zu fault vertices\n",
                    cohesive_fault_->numVertices());
    }

    // createCohesiveCellLabel creates the cohesive_cells label for diagnostics.
    // Temporarily disabled to verify it does not cause regressions.
    // ierr = createCohesiveCellLabel(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

// =============================================================================
// createCohesiveCellLabel
//
// Walk height-0 cells and label every cohesive prism cell (SEG_PRISM_TENSOR,
// TRI_PRISM_TENSOR, QUAD_PRISM_TENSOR) with value 1 in a new DMLabel named
// "cohesive_cells". After marking the cells, DMPlexLabelComplete is called to
// extend the label to include all closure points (faces, edges, vertices) so
// that PETSc correctly assigns DOFs when the label is passed to DMAddField().
// =============================================================================
PetscErrorCode Simulator::createCohesiveCellLabel()
{
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

    if (!config.enable_faults) PetscFunctionReturn(PETSC_SUCCESS);

    ierr = DMCreateLabel(dm, "cohesive_cells"); CHKERRQ(ierr);
    DMLabel cohesive_label = nullptr;
    ierr = DMGetLabel(dm, "cohesive_cells", &cohesive_label); CHKERRQ(ierr);

    PetscInt cStart, cEnd;
    ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd); CHKERRQ(ierr);

    PetscInt n_cohesive = 0;
    for (PetscInt c = cStart; c < cEnd; ++c) {
        DMPolytopeType ct;
        ierr = DMPlexGetCellType(dm, c, &ct); CHKERRQ(ierr);
        if (ct == DM_POLYTOPE_SEG_PRISM_TENSOR ||
            ct == DM_POLYTOPE_TRI_PRISM_TENSOR ||
            ct == DM_POLYTOPE_QUAD_PRISM_TENSOR) {
            ierr = DMLabelSetValue(cohesive_label, c, 1); CHKERRQ(ierr);
            ++n_cohesive;
        }
    }

    // DMPlexLabelComplete walks the closure of every labeled point and
    // adds faces, edges, and vertices to the label. Without this call,
    // PETSc's DS builder cannot allocate DOFs consistently because the
    // FE basis has DOFs at vertices but the label only marks cells.
    // This is the step PR #105 missed.
    // Pattern: PyLith libsrc/pylith/topology/MeshOps.cc line ~247
    ierr = DMPlexLabelComplete(dm, cohesive_label); CHKERRQ(ierr);

    // Diagnostic: print cells vs. total closure points to verify completion
    if (rank == 0) {
        PetscInt n_total = 0;
        IS stratum_is = nullptr;
        ierr = DMLabelGetStratumIS(cohesive_label, 1, &stratum_is); CHKERRQ(ierr);
        if (stratum_is) {
            ierr = ISGetLocalSize(stratum_is, &n_total); CHKERRQ(ierr);
            ierr = ISDestroy(&stratum_is); CHKERRQ(ierr);
        }
        PetscPrintf(comm,
            "cohesive_cells label: %d cohesive cells -> %d total closure points (DMPlexLabelComplete)\n",
            (int)n_cohesive, (int)n_total);
    }

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode Simulator::getOrCreateInterfacesLabel(DMLabel *interfacesLabel) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

    const char *labelName = "cohesive interface";
    PetscBool hasLabel = PETSC_FALSE;
    ierr = DMHasLabel(dm, labelName, &hasLabel); CHKERRQ(ierr);

    if (hasLabel) {
        ierr = DMGetLabel(dm, labelName, interfacesLabel); CHKERRQ(ierr);
        PetscFunctionReturn(PETSC_SUCCESS);
    }

    PetscInt dim = 0;
    ierr = DMGetDimension(dm, &dim); CHKERRQ(ierr);
    ierr = DMCreateLabel(dm, labelName); CHKERRQ(ierr);
    ierr = DMGetLabel(dm, labelName, interfacesLabel); CHKERRQ(ierr);

    for (PetscInt iDim = 0; iDim <= dim; ++iDim) {
        PetscInt pStart = 0, pEnd = 0, pMax = 0;
        ierr = DMPlexGetHeightStratum(dm, iDim, &pStart, &pEnd); CHKERRQ(ierr);
        ierr = DMPlexGetSimplexOrBoxCells(dm, iDim, NULL, &pMax); CHKERRQ(ierr);
        for (PetscInt p = pMax; p < pEnd; ++p) {
            ierr = DMLabelSetValue(*interfacesLabel, p, 1); CHKERRQ(ierr);
        }
    }

    if (rank == 0) {
        PetscInt n_total = 0;
        IS stratumIS = NULL;
        ierr = DMLabelGetStratumIS(*interfacesLabel, 1, &stratumIS); CHKERRQ(ierr);
        if (stratumIS) {
            ierr = ISGetLocalSize(stratumIS, &n_total); CHKERRQ(ierr);
            ierr = ISDestroy(&stratumIS); CHKERRQ(ierr);
        }
        PetscPrintf(comm, "'cohesive interface' label created: %d hybrid points\n", (int)n_total);
    }

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode Simulator::locateInjectionCell() {
    PetscFunctionBeginUser;
    if (!injection_enabled_) PetscFunctionReturn(0);
    PetscErrorCode ierr;

    Vec pointVec;
    ierr = VecCreateSeq(PETSC_COMM_SELF, 3, &pointVec); CHKERRQ(ierr);
    ierr = VecSetBlockSize(pointVec, 3); CHKERRQ(ierr);
    PetscScalar *pv;
    ierr = VecGetArray(pointVec, &pv); CHKERRQ(ierr);
    pv[0] = injection_x_; pv[1] = injection_y_; pv[2] = injection_z_;
    ierr = VecRestoreArray(pointVec, &pv); CHKERRQ(ierr);

    PetscSF cellSF = nullptr;
    ierr = DMLocatePoints(dm, pointVec, DM_POINTLOCATION_NONE, &cellSF); CHKERRQ(ierr);

    // Extract the cell index from the PetscSF
    const PetscSFNode *remotePoints;
    PetscInt nFound;
    ierr = PetscSFGetGraph(cellSF, nullptr, &nFound, nullptr, &remotePoints); CHKERRQ(ierr);
    if (nFound > 0 && remotePoints[0].index >= 0) {
        injection_cell_ = remotePoints[0].index;
    } else {
        injection_cell_ = -1;
    }

    ierr = PetscSFDestroy(&cellSF); CHKERRQ(ierr);
    ierr = VecDestroy(&pointVec); CHKERRQ(ierr);

    if (rank == 0) {
        PetscPrintf(comm, "Injection point located in cell %d\n", (int)injection_cell_);
    }
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::addInjectionToResidual(PetscReal t, Vec locF) {
    PetscFunctionBeginUser;
    if (injection_cell_ < 0) PetscFunctionReturn(0);
    PetscErrorCode ierr;

    // Compute cell volume
    PetscReal vol;
    ierr = DMPlexComputeCellGeometryFVM(dm, injection_cell_, &vol, nullptr, nullptr);
    CHKERRQ(ierr);

    // Get the local section to find pressure DOFs
    PetscSection section;
    ierr = DMGetLocalSection(dm, &section); CHKERRQ(ierr);

    // Pressure field index: 0 for poroelastic and single-phase
    PetscInt pressure_field = 0;

    // Get the number of pressure DOFs for this cell
    PetscInt nPressureDofs;
    ierr = PetscSectionGetFieldDof(section, injection_cell_, pressure_field, &nPressureDofs);
    CHKERRQ(ierr);

    if (nPressureDofs <= 0) PetscFunctionReturn(0);

    // Get the closure for this cell
    PetscScalar *closure = nullptr;
    PetscInt closureSize;
    ierr = DMPlexVecGetClosure(dm, section, locF, injection_cell_, &closureSize, &closure);
    CHKERRQ(ierr);

    // Source term: Q / V_cell distributed over pressure DOFs
    // Negative because IFunction form: F(t,U,U_t) = 0, source goes on RHS
    PetscScalar source = -injection_rate_ / vol;

    // Pressure DOFs come first in the closure for poroelastic (field 0)
    for (PetscInt i = 0; i < nPressureDofs; i++) {
        closure[i] += source / nPressureDofs;
    }

    ierr = DMPlexVecSetClosure(dm, section, locF, injection_cell_, closure, INSERT_VALUES);
    CHKERRQ(ierr);
    ierr = DMPlexVecRestoreClosure(dm, section, locF, injection_cell_, &closureSize, &closure);
    CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::addTractionToResidual(PetscReal t, Vec locF)
{
    (void)t;
    PetscFunctionBeginUser;
    if (!has_traction_bc_) PetscFunctionReturn(0);

    PetscErrorCode ierr;
    PetscSection section;
    ierr = DMGetLocalSection(dm, &section); CHKERRQ(ierr);

    PetscInt dim;
    ierr = DMGetDimension(dm, &dim); CHKERRQ(ierr);

    DMLabel depth_label;
    ierr = DMPlexGetDepthLabel(dm, &depth_label); CHKERRQ(ierr);

    PetscInt disp_field = 0;
    if (config.fluid_model == FluidModelType::SINGLE_COMPONENT) {
        disp_field = 1;
    } else if (config.fluid_model == FluidModelType::BLACK_OIL) {
        disp_field = 3;
    } else if (config.fluid_model == FluidModelType::COMPOSITIONAL) {
        disp_field = 4;
    }

    static const char *label_names[6] = {
        "boundary_x_min", "boundary_x_max",
        "boundary_y_min", "boundary_y_max",
        "boundary_z_min", "boundary_z_max"
    };

    PetscScalar *farray = nullptr;
    ierr = VecGetArray(locF, &farray); CHKERRQ(ierr);

    for (int fi = 0; fi < 6; ++fi) {
        if (!config.face_bc[fi].configured || config.face_bc[fi].type != "traction") continue;

        DMLabel face_label = nullptr;
        ierr = DMGetLabel(dm, label_names[fi], &face_label); CHKERRQ(ierr);
        if (!face_label) continue;

        IS stratum_is = nullptr;
        ierr = DMLabelGetStratumIS(face_label, 1, &stratum_is); CHKERRQ(ierr);
        if (!stratum_is) continue;

        const PetscInt *pts = nullptr;
        PetscInt npts = 0;
        ierr = ISGetLocalSize(stratum_is, &npts); CHKERRQ(ierr);
        ierr = ISGetIndices(stratum_is, &pts); CHKERRQ(ierr);

        const double *traction = config.face_bc[fi].traction;

        for (PetscInt i = 0; i < npts; ++i) {
            const PetscInt face = pts[i];

            // Only process face entities (depth dim-1)
            PetscInt depth = -1;
            ierr = DMLabelGetValue(depth_label, face, &depth); CHKERRQ(ierr);
            if (depth != dim - 1) continue;

            // Must be a boundary face (support size == 1)
            PetscInt support_size = 0;
            ierr = DMPlexGetSupportSize(dm, face, &support_size); CHKERRQ(ierr);
            if (support_size != 1) continue;

            // Get face geometry
            PetscReal face_area = 0.0;
            PetscReal centroid[3] = {0.0, 0.0, 0.0};
            PetscReal normal[3] = {0.0, 0.0, 0.0};
            ierr = DMPlexComputeCellGeometryFVM(dm, face, &face_area, centroid, normal); CHKERRQ(ierr);
            if (face_area <= 0.0) continue;

            // Get the support cell
            const PetscInt *support = nullptr;
            ierr = DMPlexGetSupport(dm, face, &support); CHKERRQ(ierr);
            const PetscInt cell = support[0];

            // Get closure of the face to find vertices on the face
            PetscInt nclosure_face = 0;
            PetscInt *closure_face = nullptr;
            ierr = DMPlexGetTransitiveClosure(dm, face, PETSC_TRUE, &nclosure_face, &closure_face); CHKERRQ(ierr);

            // Count nodes on the face that are in the cell and have displacement DOFs
            PetscInt nface_nodes = 0;
            for (PetscInt c = 0; c < 2 * nclosure_face; c += 2) {
                PetscInt point = closure_face[c];
                PetscInt fdof = 0;
                ierr = PetscSectionGetFieldDof(section, point, disp_field, &fdof); CHKERRQ(ierr);
                if (fdof > 0) nface_nodes++;
            }

            if (nface_nodes <= 0) {
                ierr = DMPlexRestoreTransitiveClosure(dm, face, PETSC_TRUE, &nclosure_face, &closure_face); CHKERRQ(ierr);
                continue;
            }

            // Distribute traction force equally among face nodes
            // Total force = traction * face_area
            // Per-node force = traction * face_area / nface_nodes
            // Sign: traction enters the residual as -t (weak form: F_int - F_ext = 0)
            for (PetscInt c = 0; c < 2 * nclosure_face; c += 2) {
                PetscInt point = closure_face[c];
                PetscInt fdof = 0;
                ierr = PetscSectionGetFieldDof(section, point, disp_field, &fdof); CHKERRQ(ierr);
                if (fdof <= 0) continue;

                PetscInt field_off = 0;
                ierr = PetscSectionGetFieldOffset(section, point, disp_field, &field_off); CHKERRQ(ierr);

                for (PetscInt d = 0; d < PetscMin(fdof, dim); ++d) {
                    farray[field_off + d] -= traction[d] * face_area / (PetscReal)nface_nodes;
                }
            }

            ierr = DMPlexRestoreTransitiveClosure(dm, face, PETSC_TRUE, &nclosure_face, &closure_face); CHKERRQ(ierr);
        }

        ierr = ISRestoreIndices(stratum_is, &pts); CHKERRQ(ierr);
        ierr = ISDestroy(&stratum_is); CHKERRQ(ierr);
    }

    ierr = VecRestoreArray(locF, &farray); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::addFaultPressureToResidual(PetscReal t, Vec locF) {
    (void)t;
    PetscFunctionBeginUser;
    if (!hydrofrac_fem_pressurized_mode_ || hydrofrac_uniform_pressure_pa_ <= 0.0 || !config.enable_faults) {
        PetscFunctionReturn(0);
    }

    PetscErrorCode ierr;
    DMLabel fault_label = nullptr;
    ierr = DMGetLabel(dm, "fault", &fault_label); CHKERRQ(ierr);
    if (!fault_label) {
        ierr = DMGetLabel(dm, "Fault", &fault_label); CHKERRQ(ierr);
    }
    if (!fault_label) {
        PetscFunctionReturn(0);
    }

    PetscInt disp_field = 0;
    if (config.fluid_model == FluidModelType::SINGLE_COMPONENT) {
        disp_field = 1;
    } else if (config.fluid_model == FluidModelType::BLACK_OIL) {
        disp_field = 3;
    } else if (config.fluid_model == FluidModelType::COMPOSITIONAL) {
        disp_field = 4;
    }

    PetscSection section;
    ierr = DMGetLocalSection(dm, &section); CHKERRQ(ierr);

    DMLabel depth_label;
    ierr = DMPlexGetDepthLabel(dm, &depth_label); CHKERRQ(ierr);
    PetscInt dim;
    ierr = DMGetDimension(dm, &dim); CHKERRQ(ierr);

    IS stratum_is = nullptr;
    ierr = DMLabelGetStratumIS(fault_label, 1, &stratum_is); CHKERRQ(ierr);
    if (!stratum_is) {
        PetscFunctionReturn(0);
    }

    const PetscInt *pts = nullptr;
    PetscInt npts = 0;
    ierr = ISGetLocalSize(stratum_is, &npts); CHKERRQ(ierr);
    ierr = ISGetIndices(stratum_is, &pts); CHKERRQ(ierr);

    PetscScalar *farray = nullptr;
    ierr = VecGetArray(locF, &farray); CHKERRQ(ierr);

    for (PetscInt i = 0; i < npts; ++i) {
        const PetscInt face = pts[i];
        PetscInt depth = -1;
        ierr = DMLabelGetValue(depth_label, face, &depth); CHKERRQ(ierr);
        if (depth != dim - 1) {
            continue;
        }

        PetscInt support_size = 0;
        ierr = DMPlexGetSupportSize(dm, face, &support_size); CHKERRQ(ierr);
        if (support_size < 2) {
            continue;
        }

        const PetscInt *support = nullptr;
        ierr = DMPlexGetSupport(dm, face, &support); CHKERRQ(ierr);

        PetscReal face_measure = 0.0;
        PetscReal centroid[3] = {0.0, 0.0, 0.0};
        PetscReal normal[3] = {0.0, 0.0, 0.0};
        ierr = DMPlexComputeCellGeometryFVM(dm, face, &face_measure, centroid, normal); CHKERRQ(ierr);
        if (face_measure <= 0.0) {
            continue;
        }

        PetscReal nmag = 0.0;
        for (PetscInt d = 0; d < dim; ++d) {
            nmag += normal[d] * normal[d];
        }
        nmag = PetscSqrtReal(nmag);
        if (nmag <= PETSC_SMALL) {
            continue;
        }

        PetscReal nunit[3] = {0.0, 0.0, 0.0};
        for (PetscInt d = 0; d < dim; ++d) {
            nunit[d] = normal[d] / nmag;
        }

        for (PetscInt side = 0; side < 2; ++side) {
            const PetscInt cell = support[side];
            const PetscReal side_sign = (side == 0) ? 1.0 : -1.0;

            PetscInt nclosure = 0;
            PetscInt *closure = nullptr;
            ierr = DMPlexGetTransitiveClosure(dm, cell, PETSC_TRUE, &nclosure, &closure); CHKERRQ(ierr);

            PetscInt total_disp_dof = 0;
            for (PetscInt c = 0; c < 2 * nclosure; c += 2) {
                PetscInt point = closure[c];
                PetscInt fdof = 0;
                ierr = PetscSectionGetFieldDof(section, point, disp_field, &fdof); CHKERRQ(ierr);
                total_disp_dof += fdof;
            }

            if (total_disp_dof > 0) {
                const PetscInt nshape = PetscMax(1, total_disp_dof / dim);
                for (PetscInt c = 0; c < 2 * nclosure; c += 2) {
                    PetscInt point = closure[c];
                    PetscInt field_off = 0;
                    PetscInt fdof = 0;
                    ierr = PetscSectionGetFieldDof(section, point, disp_field, &fdof); CHKERRQ(ierr);
                    if (fdof <= 0) {
                        continue;
                    }
                    ierr = PetscSectionGetFieldOffset(section, point, disp_field, &field_off); CHKERRQ(ierr);
                    if (field_off < 0) {
                        continue;
                    }

                    for (PetscInt d = 0; d < fdof; ++d) {
                        PetscInt comp = d % dim;
                        farray[field_off + d] +=
                            side_sign * hydrofrac_uniform_pressure_pa_ * face_measure * nunit[comp] /
                            static_cast<PetscReal>(nshape);
                    }
                }
            }

            ierr = DMPlexRestoreTransitiveClosure(dm, cell, PETSC_TRUE, &nclosure, &closure); CHKERRQ(ierr);
        }
    }

    ierr = VecRestoreArray(locF, &farray); CHKERRQ(ierr);
    ierr = ISRestoreIndices(stratum_is, &pts); CHKERRQ(ierr);
    ierr = ISDestroy(&stratum_is); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

// =============================================================================
// subtractLagrangeRegularizationOnCohesive
//
// After DMPlexTSComputeIFunctionFEM assembles the residual (which includes
// f=lambda volume regularization on ALL cells plus BdResidual on cohesive
// faces), this function subtracts the f=lambda contribution from vertices
// in the closure of cohesive cells. The net effect:
//   - Interior vertices: Lagrange residual = lambda (regularization, keeps LU happy)
//   - Cohesive vertices: Lagrange residual = BdResidual only (constraint equation)
//
// This prevents the regularization from competing with BdResidual on coarse
// meshes where f=lambda overwhelms the constraint and drives lambda to zero.
// =============================================================================
PetscErrorCode Simulator::subtractLagrangeRegularizationOnCohesive(Vec locF, Vec locU)
{
    PetscFunctionBeginUser;
    if (!config.enable_faults || !cohesive_kernel_) PetscFunctionReturn(PETSC_SUCCESS);
    PetscErrorCode ierr;

    PetscSection section = nullptr;
    ierr = DMGetLocalSection(dm, &section); CHKERRQ(ierr);

    PetscInt disp_field = 0;
    if (config.fluid_model == FluidModelType::SINGLE_COMPONENT) disp_field = 1;
    else if (config.fluid_model == FluidModelType::BLACK_OIL) disp_field = 3;
    else if (config.fluid_model == FluidModelType::COMPOSITIONAL) disp_field = 4;
    const PetscInt lagrange_field = disp_field + 1;

    PetscInt dim = 0;
    ierr = DMGetDimension(dm, &dim); CHKERRQ(ierr);

    PetscScalar *farray = nullptr;
    const PetscScalar *uarray = nullptr;
    ierr = VecGetArray(locF, &farray); CHKERRQ(ierr);
    ierr = VecGetArrayRead(locU, &uarray); CHKERRQ(ierr);

    DMLabel depth_label = nullptr;
    ierr = DMPlexGetDepthLabel(dm, &depth_label); CHKERRQ(ierr);

    // Collect all vertices in the closure of cohesive cells
    std::set<PetscInt> cohesive_vertices;
    PetscInt cStart = 0, cEnd = 0;
    ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd); CHKERRQ(ierr);

    for (PetscInt c = cStart; c < cEnd; ++c) {
        DMPolytopeType ct;
        ierr = DMPlexGetCellType(dm, c, &ct); CHKERRQ(ierr);
        if (ct != DM_POLYTOPE_SEG_PRISM_TENSOR &&
            ct != DM_POLYTOPE_TRI_PRISM_TENSOR &&
            ct != DM_POLYTOPE_QUAD_PRISM_TENSOR) {
            continue;
        }
        PetscInt closure_size = 0;
        PetscInt *closure = nullptr;
        ierr = DMPlexGetTransitiveClosure(dm, c, PETSC_TRUE,
            &closure_size, &closure); CHKERRQ(ierr);
        for (PetscInt i = 0; i < closure_size; ++i) {
            const PetscInt point = closure[2 * i];
            PetscInt depth = -1;
            ierr = DMLabelGetValue(depth_label, point, &depth); CHKERRQ(ierr);
            if (depth == 0) cohesive_vertices.insert(point);
        }
        ierr = DMPlexRestoreTransitiveClosure(dm, c, PETSC_TRUE,
            &closure_size, &closure); CHKERRQ(ierr);
    }

    // For each cohesive vertex, compute the lumped mass contribution from
    // ALL cells in its support (both interior and cohesive) and subtract
    // the f=lambda regularization from the Lagrange residual.
    // For P1 on tets: M_lumped(v) = sum_cells_c { V_c / (dim+1) }
    // The PetscDS assembles f=lambda using the consistent mass form, but
    // lumped mass is a good approximation for the correction.
    for (PetscInt v : cohesive_vertices) {
        PetscInt lag_dof = 0, lag_off = 0;
        ierr = PetscSectionGetFieldDof(section, v, lagrange_field, &lag_dof); CHKERRQ(ierr);
        if (lag_dof <= 0) continue;
        ierr = PetscSectionGetFieldOffset(section, v, lagrange_field, &lag_off); CHKERRQ(ierr);

        // Compute lumped mass at this vertex from its support cells
        PetscReal mass_lumped = 0.0;
        PetscInt support_size = 0;
        const PetscInt *support = nullptr;
        ierr = DMPlexGetTransitiveClosure(dm, v, PETSC_FALSE,
            &support_size, const_cast<PetscInt**>(&support)); CHKERRQ(ierr);
        // Actually, DMPlexGetTransitiveClosure with useCone=FALSE gives the star.
        // We need cells touching this vertex. Walk the star and pick height-0 points.
        for (PetscInt i = 0; i < support_size; ++i) {
            const PetscInt pt = support[2 * i];
            PetscInt pt_depth = -1;
            ierr = DMLabelGetValue(depth_label, pt, &pt_depth); CHKERRQ(ierr);
            if (pt_depth != dim) continue; // only cells (depth == dim)
            PetscReal vol = 0.0;
            ierr = DMPlexComputeCellGeometryFVM(dm, pt, &vol, nullptr, nullptr); CHKERRQ(ierr);
            mass_lumped += PetscAbsReal(vol) / static_cast<PetscReal>(dim + 1);
        }
        ierr = DMPlexRestoreTransitiveClosure(dm, v, PETSC_FALSE,
            &support_size, const_cast<PetscInt**>(&support)); CHKERRQ(ierr);

        // Subtract: R_lag -= mass_lumped * lambda
        // This cancels the f=lambda volume integral at this vertex
        for (PetscInt d = 0; d < PetscMin(lag_dof, dim); ++d) {
            farray[lag_off + d] -= mass_lumped * uarray[lag_off + d];
        }
    }

    ierr = VecRestoreArrayRead(locU, &uarray); CHKERRQ(ierr);
    ierr = VecRestoreArray(locF, &farray); CHKERRQ(ierr);

    PetscFunctionReturn(PETSC_SUCCESS);
}

// =============================================================================
// zeroLagrangeDiagonalOnCohesive
//
// Walk cohesive cells and zero the Lagrange-Lagrange diagonal block at each
// vertex in the closure. This removes the PetscDS g=I volume regularization
// that competes with BdResidual on cohesive vertices.
// Must be called between MatAssemblyBegin/End(MAT_FLUSH_ASSEMBLY) calls when
// switching between ADD_VALUES and INSERT_VALUES.
// =============================================================================
PetscErrorCode Simulator::zeroLagrangeDiagonalOnCohesive(Mat J)
{
    PetscFunctionBeginUser;
    if (!config.enable_faults || !cohesive_kernel_) PetscFunctionReturn(PETSC_SUCCESS);
    PetscErrorCode ierr;

    PetscSection gsection = nullptr;
    ierr = DMGetGlobalSection(dm, &gsection); CHKERRQ(ierr);
    PetscInt dim = 0;
    ierr = DMGetDimension(dm, &dim); CHKERRQ(ierr);

    PetscInt disp_field = 0;
    if (config.fluid_model == FluidModelType::SINGLE_COMPONENT) disp_field = 1;
    else if (config.fluid_model == FluidModelType::BLACK_OIL) disp_field = 3;
    else if (config.fluid_model == FluidModelType::COMPOSITIONAL) disp_field = 4;
    const PetscInt lagrange_field = disp_field + 1;

    DMLabel depth_label = nullptr;
    ierr = DMPlexGetDepthLabel(dm, &depth_label); CHKERRQ(ierr);

    PetscInt cStart = 0, cEnd = 0;
    ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd); CHKERRQ(ierr);

    for (PetscInt c = cStart; c < cEnd; ++c) {
        DMPolytopeType ct;
        ierr = DMPlexGetCellType(dm, c, &ct); CHKERRQ(ierr);
        if (ct != DM_POLYTOPE_SEG_PRISM_TENSOR &&
            ct != DM_POLYTOPE_TRI_PRISM_TENSOR &&
            ct != DM_POLYTOPE_QUAD_PRISM_TENSOR) {
            continue;
        }

        // Walk closure to find vertices
        PetscInt closure_size = 0;
        PetscInt *closure = nullptr;
        ierr = DMPlexGetTransitiveClosure(dm, c, PETSC_TRUE,
            &closure_size, &closure); CHKERRQ(ierr);

        for (PetscInt i = 0; i < closure_size; ++i) {
            const PetscInt point = closure[2 * i];
            PetscInt depth = -1;
            ierr = DMLabelGetValue(depth_label, point, &depth); CHKERRQ(ierr);
            if (depth != 0) continue;

            PetscInt lag_dof = 0, lag_off = -1;
            ierr = PetscSectionGetFieldDof(gsection, point, lagrange_field, &lag_dof); CHKERRQ(ierr);
            if (lag_dof <= 0) continue;
            ierr = PetscSectionGetFieldOffset(gsection, point, lagrange_field, &lag_off); CHKERRQ(ierr);
            if (lag_off < 0) continue;

            // Zero the full Lagrange block (d x d) at this vertex
            for (PetscInt d = 0; d < PetscMin(lag_dof, dim); ++d) {
                for (PetscInt e = 0; e < PetscMin(lag_dof, dim); ++e) {
                    ierr = MatSetValue(J, lag_off + d, lag_off + e,
                        0.0, INSERT_VALUES); CHKERRQ(ierr);
                }
            }
        }

        ierr = DMPlexRestoreTransitiveClosure(dm, c, PETSC_TRUE,
            &closure_size, &closure); CHKERRQ(ierr);
    }

    PetscFunctionReturn(PETSC_SUCCESS);
}

// =============================================================================
// addCohesivePenaltyToJacobian
//
// Analytical interface Jacobian for all cohesive fault modes:
//   - locked fault: diagonal identity coupling + penalty
//   - prescribed slip: same as locked
//   - slipping (Coulomb friction): full semi-smooth Newton tangent with
//     off-diagonal slip direction derivatives and friction-normal coupling,
//     matching CohesiveFaultKernel g0 kernels exactly.
// =============================================================================
PetscErrorCode Simulator::addCohesivePenaltyToJacobian(Mat J, Vec locU)
{
    PetscFunctionBeginUser;
    if (!config.enable_faults || !cohesive_kernel_) PetscFunctionReturn(PETSC_SUCCESS);

    // No penalty augmentation for slipping mode -- the Jacobian coupling
    // (traction +/- lambda and constraint identity) is sufficient.
    // Full penalty for locked/prescribed modes.
    const bool is_slipping = (fault_mode_ == "slipping");
    const PetscReal penalty_scale = is_slipping ? 0.0 : 1.0;

    PetscErrorCode ierr;
    PetscSection gsection = nullptr;
    PetscSection coord_section = nullptr;
    Vec coords = nullptr;
    const PetscScalar* coord_array = nullptr;
    DMLabel depth_label = nullptr;
    PetscInt dim = 0;

    ierr = DMGetGlobalSection(dm, &gsection); CHKERRQ(ierr);
    ierr = DMGetCoordinateSection(dm, &coord_section); CHKERRQ(ierr);
    ierr = DMGetCoordinatesLocal(dm, &coords); CHKERRQ(ierr);
    if (!coords)
    {
        ierr = DMGetCoordinates(dm, &coords); CHKERRQ(ierr);
    }
    PetscCheck(coords && coord_section, comm, PETSC_ERR_ARG_WRONGSTATE,
                         "Cohesive Jacobian assembly requires mesh coordinates");
    ierr = VecGetArrayRead(coords, &coord_array); CHKERRQ(ierr);
    ierr = DMPlexGetDepthLabel(dm, &depth_label); CHKERRQ(ierr);
    ierr = DMGetDimension(dm, &dim); CHKERRQ(ierr);

    PetscInt disp_field = 0;
    if (config.fluid_model == FluidModelType::SINGLE_COMPONENT)
    {
        disp_field = 1;
    }
    else if (config.fluid_model == FluidModelType::BLACK_OIL)
    {
        disp_field = 3;
    }
    else if (config.fluid_model == FluidModelType::COMPOSITIONAL)
    {
        disp_field = 4;
    }
    const PetscInt lagrange_field = disp_field + 1;
    PetscReal youngs_modulus = 10.0e9;
    if (!material_props.empty())
    {
        youngs_modulus = material_props[0].youngs_modulus;
    }

    struct FaultVertex
    {
        PetscInt point = -1;
        std::array<PetscReal, 3> xyz = {0.0, 0.0, 0.0};
    };

    auto loadCoordinates = [&](PetscInt point, std::array<PetscReal, 3>& xyz) -> PetscErrorCode
    {
        PetscInt dof = 0;
        PetscInt off = 0;
        xyz = {0.0, 0.0, 0.0};
        PetscErrorCode ierr_local = PetscSectionGetDof(coord_section, point, &dof);CHKERRQ(ierr_local);
        if (dof <= 0)
        {
            return PETSC_SUCCESS;
        }
        ierr_local = PetscSectionGetOffset(coord_section, point, &off);CHKERRQ(ierr_local);
        for (PetscInt d = 0; d < PetscMin(dof, 3); ++d)
        {
            xyz[static_cast<std::size_t>(d)] = PetscRealPart(coord_array[off + d]);
        }
        return PETSC_SUCCESS;
    };

    auto collectFaceVertices = [&](PetscInt face, std::vector<FaultVertex>& vertices) -> PetscErrorCode
    {
        PetscInt closure_size = 0;
        PetscInt* closure = nullptr;
        PetscErrorCode ierr_local = DMPlexGetTransitiveClosure(dm, face, PETSC_TRUE,
                                                                                                                     &closure_size, &closure);CHKERRQ(ierr_local);
        vertices.clear();
        for (PetscInt i = 0; i < closure_size; ++i)
        {
            const PetscInt point = closure[2 * i];
            PetscInt depth = -1;
            ierr_local = DMLabelGetValue(depth_label, point, &depth);CHKERRQ(ierr_local);
            if (depth != 0)
            {
                continue;
            }
            FaultVertex vertex;
            vertex.point = point;
            ierr_local = loadCoordinates(point, vertex.xyz);CHKERRQ(ierr_local);
            vertices.push_back(vertex);
        }
        ierr_local = DMPlexRestoreTransitiveClosure(dm, face, PETSC_TRUE,
                                                                                                &closure_size, &closure);CHKERRQ(ierr_local);
        std::sort(vertices.begin(), vertices.end(),
                            [](const FaultVertex& lhs, const FaultVertex& rhs)
                            {
                                return lhs.xyz < rhs.xyz;
                            });
        return PETSC_SUCCESS;
    };

    PetscInt cStart = 0;
    PetscInt cEnd = 0;
    ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd); CHKERRQ(ierr);
    for (PetscInt c = cStart; c < cEnd; ++c)
    {
        DMPolytopeType ct;
        ierr = DMPlexGetCellType(dm, c, &ct); CHKERRQ(ierr);
        const bool is_cohesive = (ct == DM_POLYTOPE_SEG_PRISM_TENSOR ||
                                                            ct == DM_POLYTOPE_TRI_PRISM_TENSOR ||
                                                            ct == DM_POLYTOPE_QUAD_PRISM_TENSOR);
        if (!is_cohesive)
        {
            continue;
        }

        const PetscInt* cone = nullptr;
        PetscInt cone_size = 0;
        ierr = DMPlexGetConeSize(dm, c, &cone_size); CHKERRQ(ierr);
        if (cone_size < 2)
        {
            continue;
        }
        ierr = DMPlexGetCone(dm, c, &cone); CHKERRQ(ierr);

        std::vector<FaultVertex> neg_vertices;
        std::vector<FaultVertex> pos_vertices;
        ierr = collectFaceVertices(cone[0], neg_vertices); CHKERRQ(ierr);
        ierr = collectFaceVertices(cone[1], pos_vertices); CHKERRQ(ierr);

        const PetscInt n_pairs = static_cast<PetscInt>(std::min(neg_vertices.size(), pos_vertices.size()));
        if (n_pairs <= 0)
        {
            continue;
        }

        PetscReal face_area = 0.0;
        PetscReal centroid[3] = {0.0, 0.0, 0.0};
        PetscReal normal[3] = {0.0, 0.0, 0.0};
        ierr = DMPlexComputeCellGeometryFVM(dm, cone[0], &face_area, centroid, normal); CHKERRQ(ierr);
        if (face_area <= 0.0)
        {
            continue;
        }

        // Normalize the face normal (DMPlexComputeCellGeometryFVM may return unnormalized)
        PetscReal normal_mag_j = 0.0;
        for (PetscInt d = 0; d < dim; ++d) normal_mag_j += normal[d] * normal[d];
        normal_mag_j = PetscSqrtReal(normal_mag_j);
        if (normal_mag_j <= PETSC_SMALL) continue;
        for (PetscInt d = 0; d < dim; ++d) normal[d] /= normal_mag_j;

        const PetscScalar coeff = static_cast<PetscScalar>(face_area / static_cast<PetscReal>(n_pairs));
        const PetscReal h_char_jac = PetscMax(PetscSqrtReal(face_area), PETSC_SMALL);
        const PetscScalar penalty = static_cast<PetscScalar>(penalty_scale * 10.0 * youngs_modulus / h_char_jac);

        // For slipping mode, read solution data to compute friction derivatives
        PetscSection lsection = nullptr;
        const PetscScalar *locU_array = nullptr;
        if (is_slipping && locU) {
            ierr = DMGetLocalSection(dm, &lsection); CHKERRQ(ierr);
            ierr = VecGetArrayRead(locU, &locU_array); CHKERRQ(ierr);
        }

        for (PetscInt i = 0; i < n_pairs; ++i)
        {
            const FaultVertex& neg_vertex = neg_vertices[static_cast<std::size_t>(i)];
            const FaultVertex& pos_vertex = pos_vertices[static_cast<std::size_t>(i)];

            PetscReal pair_distance2 = 0.0;
            for (PetscInt d = 0; d < dim; ++d)
            {
                const PetscReal delta = neg_vertex.xyz[static_cast<std::size_t>(d)] -
                                                                pos_vertex.xyz[static_cast<std::size_t>(d)];
                pair_distance2 += delta * delta;
            }
            if (pair_distance2 > 1.0e-80)
            {
                continue;
            }

            PetscInt neg_disp_dof = 0;
            PetscInt pos_disp_dof = 0;
            PetscInt neg_lag_dof = 0;
            PetscInt neg_disp_off = 0;
            PetscInt pos_disp_off = 0;
            PetscInt neg_lag_off = 0;
            ierr = PetscSectionGetFieldDof(gsection, neg_vertex.point, disp_field, &neg_disp_dof); CHKERRQ(ierr);
            ierr = PetscSectionGetFieldDof(gsection, pos_vertex.point, disp_field, &pos_disp_dof); CHKERRQ(ierr);
            ierr = PetscSectionGetFieldDof(gsection, neg_vertex.point, lagrange_field, &neg_lag_dof); CHKERRQ(ierr);
            ierr = PetscSectionGetFieldOffset(gsection, neg_vertex.point, disp_field, &neg_disp_off); CHKERRQ(ierr);
            ierr = PetscSectionGetFieldOffset(gsection, pos_vertex.point, disp_field, &pos_disp_off); CHKERRQ(ierr);
            ierr = PetscSectionGetFieldOffset(gsection, neg_vertex.point, lagrange_field, &neg_lag_off); CHKERRQ(ierr);

            if (neg_disp_dof <= 0 || pos_disp_dof <= 0 || neg_lag_dof <= 0)
            {
                continue;
            }
            if (neg_disp_off < 0 || pos_disp_off < 0 || neg_lag_off < 0)
            {
                continue;
            }

            const PetscInt ndof = PetscMin(PetscMin(neg_disp_dof, pos_disp_dof), PetscMin(neg_lag_dof, dim));

            if (is_slipping && locU_array && lsection) {
                // =============================================================
                // Semi-smooth Newton Jacobian for Coulomb friction
                // =============================================================
                // Read current displacement and Lagrange multiplier from solution
                PetscReal u_neg_v[3] = {0,0,0}, u_pos_v[3] = {0,0,0}, lam[3] = {0,0,0};
                {
                    PetscInt dof_l = 0, off_l = 0;
                    ierr = PetscSectionGetFieldDof(lsection, neg_vertex.point, disp_field, &dof_l); CHKERRQ(ierr);
                    if (dof_l > 0) {
                        ierr = PetscSectionGetFieldOffset(lsection, neg_vertex.point, disp_field, &off_l); CHKERRQ(ierr);
                        for (PetscInt d = 0; d < PetscMin(dof_l, dim); ++d)
                            u_neg_v[d] = PetscRealPart(locU_array[off_l + d]);
                    }
                    ierr = PetscSectionGetFieldDof(lsection, pos_vertex.point, disp_field, &dof_l); CHKERRQ(ierr);
                    if (dof_l > 0) {
                        ierr = PetscSectionGetFieldOffset(lsection, pos_vertex.point, disp_field, &off_l); CHKERRQ(ierr);
                        for (PetscInt d = 0; d < PetscMin(dof_l, dim); ++d)
                            u_pos_v[d] = PetscRealPart(locU_array[off_l + d]);
                    }
                    ierr = PetscSectionGetFieldDof(lsection, neg_vertex.point, lagrange_field, &dof_l); CHKERRQ(ierr);
                    if (dof_l > 0) {
                        ierr = PetscSectionGetFieldOffset(lsection, neg_vertex.point, lagrange_field, &off_l); CHKERRQ(ierr);
                        for (PetscInt d = 0; d < PetscMin(dof_l, dim); ++d)
                            lam[d] = PetscRealPart(locU_array[off_l + d]);
                    }
                }

                // Compute slip = u_pos - u_neg
                PetscReal slip[3] = {0,0,0};
                for (PetscInt d = 0; d < dim; ++d) slip[d] = u_pos_v[d] - u_neg_v[d];

                // Normal/tangential decomposition
                PetscReal slip_n = 0.0, lam_n = 0.0;
                for (PetscInt d = 0; d < dim; ++d) {
                    slip_n += slip[d] * normal[d];
                    lam_n += lam[d] * normal[d];
                }

                PetscReal slip_t[3] = {0,0,0}, lam_t[3] = {0,0,0};
                PetscReal slip_t_mag2 = 0.0;
                for (PetscInt d = 0; d < dim; ++d) {
                    slip_t[d] = slip[d] - slip_n * normal[d];
                    lam_t[d] = lam[d] - lam_n * normal[d];
                    slip_t_mag2 += slip_t[d] * slip_t[d];
                }
                PetscReal slip_t_mag = PetscSqrtReal(slip_t_mag2);

                // Friction coefficient: constant or slip-weakening
                PetscReal mu_f;
                if (use_slip_weakening_) {
                    mu_f = mu_static_ - (mu_static_ - mu_dynamic_) *
                           PetscMin(slip_t_mag, critical_slip_distance_) / critical_slip_distance_;
                } else {
                    mu_f = fault_friction_coefficient_;
                }
                // Friction strength: tau_f = mu_f * |compressive normal traction|
                // Convention: lam_n > 0 is tensile, compressive = -lam_n
                const PetscReal sigma_n_comp = PetscMax(0.0, -lam_n);
                const PetscReal tau_f = mu_f * sigma_n_comp;
                const PetscReal slip_eps = 1.0e-14;

                // Compute slip direction (s_hat) if slipping
                PetscReal s_hat[3] = {0,0,0};
                if (slip_t_mag > slip_eps) {
                    for (PetscInt d = 0; d < dim; ++d) s_hat[d] = slip_t[d] / slip_t_mag;
                }

                // Augmented Lagrangian penalty: moderate penalty for the semi-smooth
                // Newton method. Too large causes overshoot, too small gives slow convergence.
                const PetscReal eff_penalty = 0.1 * 10.0 * youngs_modulus / h_char_jac;

                // ---------------------------------------------------------
                // Assemble full Jacobian blocks for all (d, e) pairs.
                // The constraint is:
                //   C_d = lam_t_d - tau_f * s_hat_d  (tangential)
                //       + max(-slip_n, 0) * n_d      (normal non-penetration)
                // Derivatives w.r.t. u (via slip) and lambda are assembled
                // in the (d,e) double loop below, matching the pointwise
                // kernels g0_lagrange_displacement and g0_lagrange_lagrange
                // in CohesiveFaultKernel.cpp exactly.
                // ---------------------------------------------------------

                // Full semi-smooth Newton Jacobian for Coulomb friction.
                // Matches the pointwise kernels in CohesiveFaultKernel.cpp
                // (g0_lagrange_displacement, g0_lagrange_lagrange,
                //  g0_displacement_lagrange) exactly, with penalty augmentation.
                for (PetscInt d = 0; d < ndof; ++d)
                {
                    const PetscInt row_neg_disp = neg_disp_off + d;
                    const PetscInt row_pos_disp = pos_disp_off + d;
                    const PetscInt row_lag = neg_lag_off + d;

                    for (PetscInt e = 0; e < ndof; ++e)
                    {
                        const PetscInt col_neg_disp = neg_disp_off + e;
                        const PetscInt col_pos_disp = pos_disp_off + e;
                        const PetscInt col_lag = neg_lag_off + e;

                        const PetscReal P_de = ((d == e) ? 1.0 : 0.0) - normal[d] * normal[e];

                        // --- J_lam_lam: d(C_d)/d(lam_e) ---
                        // From g0_lagrange_lagrange:
                        //   d(lam_t_d)/d(lam_e) = P_de
                        //   d(-tau_f * s_hat_d)/d(lam_e) = mu_f * s_hat_d * n_e
                        PetscReal J_ll = P_de;
                        if (sigma_n_comp > 0.0 && slip_t_mag > slip_eps) {
                            J_ll += mu_f * s_hat[d] * normal[e];
                        }
                        ierr = MatSetValue(J, row_lag, col_lag,
                            static_cast<PetscScalar>(J_ll) * coeff, ADD_VALUES); CHKERRQ(ierr);

                        // --- J_lam_u: d(C_d)/d(u_e) ---
                        // From g0_lagrange_displacement:
                        //   d(s_hat_d)/d(u_e) = (P_de - s_hat_d * s_hat_e) / |slip_t|
                        //   d(max(-slip_n,0)*n_d)/d(u_e) = -H(-slip_n) * n_d * n_e
                        PetscReal dS_de = 0.0;
                        if (slip_t_mag > slip_eps) {
                            dS_de = (P_de - s_hat[d] * s_hat[e]) / slip_t_mag;
                        }
                        PetscReal J_lu = -tau_f * dS_de;
                        // Slip-weakening: d(mu_f)/d(|slip_t|) * sigma_n * s_hat_d * s_hat_e
                        if (use_slip_weakening_ && slip_t_mag < critical_slip_distance_
                            && slip_t_mag > slip_eps) {
                            PetscReal dmu_dslip = -(mu_static_ - mu_dynamic_) / critical_slip_distance_;
                            J_lu += -dmu_dslip * sigma_n_comp * s_hat[d] * s_hat[e];
                        }
                        if (slip_n < 0.0) {
                            J_lu += -normal[d] * normal[e];
                        }
                        // slip = u_pos - u_neg: d(slip)/d(u_pos)=+1, d(slip)/d(u_neg)=-1
                        ierr = MatSetValue(J, row_lag, col_pos_disp,
                            static_cast<PetscScalar>(J_lu) * coeff, ADD_VALUES); CHKERRQ(ierr);
                        ierr = MatSetValue(J, row_lag, col_neg_disp,
                            static_cast<PetscScalar>(-J_lu) * coeff, ADD_VALUES); CHKERRQ(ierr);

                        // --- J_u_lam: d(R_u)/d(lambda) = +/- delta_de ---
                        // From g0_displacement_lagrange (identity)
                        PetscReal J_ul = (d == e) ? 1.0 : 0.0;
                        ierr = MatSetValue(J, row_neg_disp, col_lag,
                            static_cast<PetscScalar>(-J_ul) * coeff, ADD_VALUES); CHKERRQ(ierr);
                        ierr = MatSetValue(J, row_pos_disp, col_lag,
                            static_cast<PetscScalar>(J_ul) * coeff, ADD_VALUES); CHKERRQ(ierr);

                        // --- J_u_u: penalty augmentation from constraint ---
                        PetscReal pen_uu = eff_penalty * J_lu;
                        ierr = MatSetValue(J, row_neg_disp, col_neg_disp,
                            static_cast<PetscScalar>(pen_uu) * coeff, ADD_VALUES); CHKERRQ(ierr);
                        ierr = MatSetValue(J, row_neg_disp, col_pos_disp,
                            static_cast<PetscScalar>(-pen_uu) * coeff, ADD_VALUES); CHKERRQ(ierr);
                        ierr = MatSetValue(J, row_pos_disp, col_neg_disp,
                            static_cast<PetscScalar>(-pen_uu) * coeff, ADD_VALUES); CHKERRQ(ierr);
                        ierr = MatSetValue(J, row_pos_disp, col_pos_disp,
                            static_cast<PetscScalar>(pen_uu) * coeff, ADD_VALUES); CHKERRQ(ierr);
                    }
                }
            } else {
                // =============================================================
                // Locked / prescribed slip: diagonal coupling + penalty
                // =============================================================
                for (PetscInt d = 0; d < ndof; ++d)
                {
                    const PetscInt row_neg_disp = neg_disp_off + d;
                    const PetscInt row_pos_disp = pos_disp_off + d;
                    const PetscInt row_lag = neg_lag_off + d;
                    const PetscInt col_neg_disp = neg_disp_off + d;
                    const PetscInt col_pos_disp = pos_disp_off + d;
                    const PetscInt col_lag = neg_lag_off + d;

                    ierr = MatSetValue(J, row_neg_disp, col_lag, -coeff, ADD_VALUES); CHKERRQ(ierr);
                    ierr = MatSetValue(J, row_pos_disp, col_lag, coeff, ADD_VALUES); CHKERRQ(ierr);
                    ierr = MatSetValue(J, row_lag, col_neg_disp, -coeff, ADD_VALUES); CHKERRQ(ierr);
                    ierr = MatSetValue(J, row_lag, col_pos_disp, coeff, ADD_VALUES); CHKERRQ(ierr);
                    // Penalty-scaled Lagrange diagonal at cohesive vertices.
                    // Makes the Lagrange diagonal O(penalty) at fault vertices,
                    // matching displacement stiffness scale for LU conditioning.
                    // The penalty acts as augmented Lagrangian; it slows
                    // convergence for the Lagrange multiplier but does not
                    // prevent convergence (the BdResidual drives lambda to
                    // the correct traction).
                    ierr = MatSetValue(J, row_lag, col_lag, penalty * coeff, ADD_VALUES); CHKERRQ(ierr);
                    ierr = MatSetValue(J, row_neg_disp, col_neg_disp, penalty * coeff, ADD_VALUES); CHKERRQ(ierr);
                    ierr = MatSetValue(J, row_neg_disp, col_pos_disp, -penalty * coeff, ADD_VALUES); CHKERRQ(ierr);
                    ierr = MatSetValue(J, row_pos_disp, col_neg_disp, -penalty * coeff, ADD_VALUES); CHKERRQ(ierr);
                    ierr = MatSetValue(J, row_pos_disp, col_pos_disp, penalty * coeff, ADD_VALUES); CHKERRQ(ierr);
                }
            }
        }

        if (is_slipping && locU_array) {
            ierr = VecRestoreArrayRead(locU, &locU_array); CHKERRQ(ierr);
        }
    }

    ierr = VecRestoreArrayRead(coords, &coord_array); CHKERRQ(ierr);
    ierr = MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode Simulator::locateExplosionCell() {
    PetscFunctionBeginUser;
    if (!explosion_) PetscFunctionReturn(0);
    PetscErrorCode ierr;

    Vec pointVec;
    ierr = VecCreateSeq(PETSC_COMM_SELF, 3, &pointVec); CHKERRQ(ierr);
    ierr = VecSetBlockSize(pointVec, 3); CHKERRQ(ierr);
    PetscScalar *pv;
    ierr = VecGetArray(pointVec, &pv); CHKERRQ(ierr);
    pv[0] = explosion_->sx; pv[1] = explosion_->sy; pv[2] = explosion_->sz;
    ierr = VecRestoreArray(pointVec, &pv); CHKERRQ(ierr);

    PetscSF cellSF = nullptr;
    ierr = DMLocatePoints(dm, pointVec, DM_POINTLOCATION_NONE, &cellSF); CHKERRQ(ierr);

    const PetscSFNode *remotePoints;
    PetscInt nFound;
    ierr = PetscSFGetGraph(cellSF, nullptr, &nFound, nullptr, &remotePoints); CHKERRQ(ierr);
    if (nFound > 0 && remotePoints[0].index >= 0) {
        explosion_cell_ = remotePoints[0].index;
    } else {
        explosion_cell_ = -1;
    }

    ierr = PetscSFDestroy(&cellSF); CHKERRQ(ierr);
    ierr = VecDestroy(&pointVec); CHKERRQ(ierr);

    if (rank == 0) {
        PetscPrintf(comm, "Explosion source located in cell %d\n", (int)explosion_cell_);
    }
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::addExplosionSourceToResidual(PetscReal t, Vec locF) {
    PetscFunctionBeginUser;
    if (explosion_cell_ < 0 || !explosion_) PetscFunctionReturn(0);
    PetscErrorCode ierr;

    // Get moment tensor at current time
    double M[6];  // Mxx, Myy, Mzz, Mxy, Mxz, Myz
    double mr = 0.0;
    double elapsed = static_cast<double>(t) - explosion_->t0;
    if (elapsed >= 0.0) {
        if (explosion_->use_nearfield_coupling) {
            // COUPLED_ANALYTIC: moment rate from RDP derivative
            // M_dot(t) = 4 * pi * rho * vp^2 * psi_dot(t)
            double psi_dot = explosion_->rdp_source.psiDot(elapsed);
            double K = explosion_->rho * explosion_->vp * explosion_->vp;
            mr = 4.0 * M_PI * K * psi_dot;
        } else {
            // PROXY: Mueller-Murphy moment rate
            NuclearSourceParameters nsp;
            nsp.yield_kt = explosion_->yield_kt;
            nsp.depth_of_burial = explosion_->depth_of_burial;
            MuellerMurphySource mm;
            mm.setMediumProperties(explosion_->rho, explosion_->vp, explosion_->vs);
            mm.setParameters(nsp);
            mr = mm.momentRate(elapsed);
        }
    }
    double scale = mr / 3.0;
    M[0] = M[1] = M[2] = scale;
    M[3] = M[4] = M[5] = 0.0;

    // Get displacement field index (0 for elastostatics, 1 for poroelastic)
    PetscInt disp_field = 0;
    if (config.solid_model == SolidModelType::POROELASTIC &&
        config.fluid_model == FluidModelType::SINGLE_COMPONENT) {
        disp_field = 1;
    }

    PetscSection section;
    ierr = DMGetLocalSection(dm, &section); CHKERRQ(ierr);
    DMLabel depth_label = nullptr;
    ierr = DMPlexGetDepthLabel(dm, &depth_label); CHKERRQ(ierr);
    PetscScalar *farray = nullptr;
    ierr = VecGetArray(locF, &farray); CHKERRQ(ierr);

    // Proper moment tensor equivalent nodal forces (BUG 5 fix)
    // For a point moment tensor M_ij * delta(x - x_s):
    //   F_i^a = -(sum_j dPhi_a/dx_j * M_ij)
    // where dPhi_a is the scalar basis function of node a.
    PetscInt dim = 3;

    // Build full symmetric moment tensor [3][3]
    double Mmat[3][3] = {
        {M[0], M[3], M[4]},
        {M[3], M[1], M[5]},
        {M[4], M[5], M[2]}
    };

    // Get the FE for the displacement field
    PetscFE fe;
    ierr = DMGetField(dm, disp_field, NULL, (PetscObject*)&fe); CHKERRQ(ierr);

    // Get total number of basis functions and number of components
    PetscInt Nb;
    ierr = PetscFEGetDimension(fe, &Nb); CHKERRQ(ierr);
    PetscInt nNodes = Nb / dim;

    // Determine reference centroid based on cell type
    DMPolytopeType cellType;
    ierr = DMPlexGetCellType(dm, explosion_cell_, &cellType); CHKERRQ(ierr);

    PetscReal refPt[3];
    if (cellType == DM_POLYTOPE_TETRAHEDRON) {
        refPt[0] = 0.25; refPt[1] = 0.25; refPt[2] = 0.25;
    } else {
        refPt[0] = 0.0; refPt[1] = 0.0; refPt[2] = 0.0;
    }

    // Tabulate basis derivatives at cell centroid (K=1 for first derivatives)
    PetscTabulation T;
    ierr = PetscFECreateTabulation(fe, 1, 1, refPt, 1, &T); CHKERRQ(ierr);

    // Get inverse Jacobian for reference-to-physical coordinate mapping
    PetscReal v0[3], J[9], invJ_g[9], detJ_val;
    ierr = DMPlexComputeCellGeometryFEM(dm, explosion_cell_, NULL, v0, J, invJ_g, &detJ_val);
    CHKERRQ(ierr);

    // T->T[1] layout: [b * Nc * cdim + c * cdim + d]
    // For vector FE with dim components:
    //   basis b = node_a * dim + comp_k
    //   dPhi_a/dxi_d = T->T[1][(a*dim + k) * dim * dim + k * dim + d] (for any k with c=k)
    // Using k=0, c=0: dPhi_a/dxi_d = T->T[1][a * dim^3 + d]
    // Physical gradient: dPhi_a/dx_j = sum_d dPhi_a/dxi_d * invJ[d*dim+j]

    PetscInt closureSize = 0;
    PetscInt *cellClosure = nullptr;
    ierr = DMPlexGetTransitiveClosure(dm, explosion_cell_, PETSC_TRUE, &closureSize, &cellClosure);
    CHKERRQ(ierr);

    std::vector<PetscInt> cell_vertices;
    cell_vertices.reserve(static_cast<std::size_t>(nNodes));
    for (PetscInt c = 0; c < closureSize; ++c) {
        const PetscInt point = cellClosure[2 * c];
        PetscInt depth = -1;
        PetscInt disp_dof = 0;
        ierr = DMLabelGetValue(depth_label, point, &depth); CHKERRQ(ierr);
        if (depth != 0) {
            continue;
        }
        ierr = PetscSectionGetFieldDof(section, point, disp_field, &disp_dof); CHKERRQ(ierr);
        if (disp_dof >= dim) {
            cell_vertices.push_back(point);
        }
    }

    PetscCheck(static_cast<PetscInt>(cell_vertices.size()) >= nNodes, comm, PETSC_ERR_PLIB,
               "Explosion source cell %d has %d displacement vertices, expected at least %d",
               static_cast<int>(explosion_cell_), static_cast<int>(cell_vertices.size()), static_cast<int>(nNodes));

    for (PetscInt a = 0; a < nNodes; ++a) {
        // Compute physical gradient of scalar basis function Phi_a
        double dPhi_dx[3] = {0.0, 0.0, 0.0};
        for (PetscInt j = 0; j < dim; ++j) {
            for (PetscInt d = 0; d < dim; ++d) {
                PetscInt tabIdx = a * dim * dim * dim + d;
                dPhi_dx[j] += T->T[1][tabIdx] * invJ_g[d * dim + j];
            }
        }

        // F_i^a = -(sum_j dPhi_a/dx_j * M_ij)
        PetscInt field_off = 0;
        ierr = PetscSectionGetFieldOffset(section, cell_vertices[static_cast<std::size_t>(a)],
                                          disp_field, &field_off); CHKERRQ(ierr);
        if (field_off < 0) {
            continue;
        }
        for (PetscInt i = 0; i < dim; ++i) {
            double force = 0.0;
            for (PetscInt j = 0; j < dim; ++j) {
                force += dPhi_dx[j] * Mmat[i][j];
            }
            farray[field_off + i] += -force;
        }
    }

    ierr = DMPlexRestoreTransitiveClosure(dm, explosion_cell_, PETSC_TRUE, &closureSize, &cellClosure);
    CHKERRQ(ierr);
    ierr = PetscTabulationDestroy(&T); CHKERRQ(ierr);
    ierr = VecRestoreArray(locF, &farray); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::setInitialConditions() {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

    // Locate injection cell if injection is enabled (needs DM fully set up)
    ierr = locateInjectionCell(); CHKERRQ(ierr);

    // Locate explosion source cell
    ierr = locateExplosionCell(); CHKERRQ(ierr);

    // Load initial conditions from file if requested
    if (config.initial_condition_type == "from_file" && !config.initial_condition_path.empty()) {
        PetscViewer viewer;
        ierr = PetscViewerBinaryOpen(comm, config.initial_condition_path.c_str(),
                                     FILE_MODE_READ, &viewer); CHKERRQ(ierr);
        ierr = VecLoad(solution, viewer); CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
        if (rank == 0) {
            PetscPrintf(comm, "  Loaded initial conditions from %s\n",
                        config.initial_condition_path.c_str());
        }
        PetscFunctionReturn(0);
    }

    // Start from zero (displacement at rest, zero pressure perturbation)
    ierr = VecZeroEntries(solution); CHKERRQ(ierr);

    // For poroelasticity: set initial pore pressure (undrained response to sudden load)
    if (config.solid_model == SolidModelType::POROELASTIC &&
        config.fluid_model == FluidModelType::SINGLE_COMPONENT) {
        PetscInt nfields;
        ierr = DMGetNumFields(dm, &nfields); CHKERRQ(ierr);
        if (nfields >= 3 && config.enable_faults) {
            // Hybrid mesh with cohesive cells: DMProjectFunction fails on cohesive
            // cells. Manually set pressure DOFs to 1 MPa via the PetscSection.
            PetscSection section;
            ierr = DMGetLocalSection(dm, &section); CHKERRQ(ierr);
            Vec localVec_ic;
            ierr = DMGetLocalVector(dm, &localVec_ic); CHKERRQ(ierr);
            ierr = VecZeroEntries(localVec_ic); CHKERRQ(ierr);
            PetscInt pStart, pEnd;
            ierr = PetscSectionGetChart(section, &pStart, &pEnd); CHKERRQ(ierr);
            PetscScalar *larray;
            ierr = VecGetArray(localVec_ic, &larray); CHKERRQ(ierr);
            for (PetscInt p = pStart; p < pEnd; ++p) {
                PetscInt dof, off, fdof, foff;
                ierr = PetscSectionGetDof(section, p, &dof); CHKERRQ(ierr);
                if (dof <= 0) continue;
                ierr = PetscSectionGetOffset(section, p, &off); CHKERRQ(ierr);
                // Field 0 = pressure (1 component)
                ierr = PetscSectionGetFieldDof(section, p, 0, &fdof); CHKERRQ(ierr);
                if (fdof > 0) {
                    ierr = PetscSectionGetFieldOffset(section, p, 0, &foff); CHKERRQ(ierr);
                    for (PetscInt i = 0; i < fdof; ++i) {
                        larray[off + foff + i] = 1.0e6;  // 1 MPa
                    }
                }
            }
            ierr = VecRestoreArray(localVec_ic, &larray); CHKERRQ(ierr);
            ierr = DMLocalToGlobal(dm, localVec_ic, INSERT_VALUES, solution); CHKERRQ(ierr);
            ierr = DMRestoreLocalVector(dm, &localVec_ic); CHKERRQ(ierr);
        } else {
            // Standard mesh: use DMProjectFunction
            // terzaghi_ic returns 1 MPa for pressure (Nc=1) and 0 for displacement (Nc=3)
            PetscErrorCode (*funcs[2])(PetscInt, PetscReal, const PetscReal[], PetscInt, PetscScalar*, void*) = {terzaghi_ic, terzaghi_ic};
            void *ctxs[2] = {nullptr, nullptr};
            ierr = DMProjectFunction(dm, 0.0, funcs, ctxs, INSERT_VALUES, solution); CHKERRQ(ierr);
        }
        if (rank == 0) {
            PetscPrintf(comm, "  Set Terzaghi initial conditions: p0 = 1 MPa\n");
        }
    }

    // Ensure DM section is created (required for DMPlexInsertBoundaryValues to work)
    ierr = DMSetUp(dm); CHKERRQ(ierr);

    // Check how many boundaries are registered on the DS
    PetscDS ds;
    PetscInt numBd;
    ierr = DMGetDS(dm, &ds); CHKERRQ(ierr);
    ierr = PetscDSGetNumBoundary(ds, &numBd); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "  Number of boundaries registered: %d\n", numBd); CHKERRQ(ierr);

    // Insert boundary condition values into solution vector at boundary DOFs
    // This creates a non-equilibrium initial state for SNES to solve
    // DMPlexInsertBoundaryValues requires a local vector
    Vec localVec;
    ierr = DMGetLocalVector(dm, &localVec); CHKERRQ(ierr);
    ierr = DMGlobalToLocal(dm, solution, INSERT_VALUES, localVec); CHKERRQ(ierr);
    ierr = DMPlexInsertBoundaryValues(dm, PETSC_TRUE, localVec, 0.0, NULL, NULL, NULL); CHKERRQ(ierr);
    ierr = DMLocalToGlobal(dm, localVec, INSERT_VALUES, solution); CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(dm, &localVec); CHKERRQ(ierr);

    PetscReal norm;
    ierr = VecNorm(solution, NORM_2, &norm); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "  Solution norm after BC insertion: %.6e\n", norm); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

// =============================================================================
// applyInitialFaultStress
//
// Set Lagrange multiplier DOFs on cohesive cells to prescribed initial traction.
// Used for SCEC TPV-type benchmarks where the fault starts with known stress
// state. Walks cohesive cells and sets lambda = tau * strike_dir + sigma_n * normal
// at each fault vertex, with an optional overstressed nucleation patch.
// =============================================================================
PetscErrorCode Simulator::applyInitialFaultStress()
{
    PetscFunctionBeginUser;
    if (!fault_has_initial_stress_ || !config.enable_faults || !cohesive_kernel_)
        PetscFunctionReturn(PETSC_SUCCESS);

    PetscErrorCode ierr;
    PetscSection section = nullptr;
    PetscSection coord_section = nullptr;
    PetscInt dim = 0;
    Vec coords = nullptr;
    DMLabel depth_label = nullptr;

    ierr = DMGetLocalSection(dm, &section); CHKERRQ(ierr);
    ierr = DMGetCoordinateSection(dm, &coord_section); CHKERRQ(ierr);
    ierr = DMGetDimension(dm, &dim); CHKERRQ(ierr);
    ierr = DMGetCoordinatesLocal(dm, &coords); CHKERRQ(ierr);
    if (!coords) { ierr = DMGetCoordinates(dm, &coords); CHKERRQ(ierr); }
    ierr = DMPlexGetDepthLabel(dm, &depth_label); CHKERRQ(ierr);

    PetscInt disp_field = 0;
    if (config.fluid_model == FluidModelType::SINGLE_COMPONENT) disp_field = 1;
    else if (config.fluid_model == FluidModelType::BLACK_OIL) disp_field = 3;
    else if (config.fluid_model == FluidModelType::COMPOSITIONAL) disp_field = 4;
    const PetscInt lagrange_field = disp_field + 1;

    // Compute fault-plane directions from strike/dip
    const PetscReal strike = fault_geometry_from_config_ ? fault_strike_ : 0.0;
    const PetscReal dip = fault_geometry_from_config_ ? fault_dip_ : (PetscReal)(M_PI / 2.0);
    const PetscReal cs = PetscCosReal(strike);
    const PetscReal ss = PetscSinReal(strike);
    const PetscReal cd = PetscCosReal(dip);
    const PetscReal sd = PetscSinReal(dip);
    const PetscReal strike_dir[3] = {cs, ss, 0.0};
    const PetscReal normal_dir[3] = {-ss * sd, cs * sd, -cd};

    // Work on local vector
    Vec localVec;
    ierr = DMGetLocalVector(dm, &localVec); CHKERRQ(ierr);
    ierr = DMGlobalToLocal(dm, solution, INSERT_VALUES, localVec); CHKERRQ(ierr);

    PetscScalar *larray = nullptr;
    const PetscScalar *coord_array = nullptr;
    ierr = VecGetArray(localVec, &larray); CHKERRQ(ierr);
    ierr = VecGetArrayRead(coords, &coord_array); CHKERRQ(ierr);

    struct FaultVertex {
        PetscInt point = -1;
        std::array<PetscReal, 3> xyz = {0.0, 0.0, 0.0};
    };

    auto loadCoordinates = [&](PetscInt point, std::array<PetscReal,3>& xyz) -> PetscErrorCode {
        PetscInt dof = 0, off = 0;
        xyz = {0.0, 0.0, 0.0};
        PetscErrorCode e = PetscSectionGetDof(coord_section, point, &dof); CHKERRQ(e);
        if (dof <= 0) return PETSC_SUCCESS;
        e = PetscSectionGetOffset(coord_section, point, &off); CHKERRQ(e);
        for (PetscInt d = 0; d < PetscMin(dof, 3); ++d)
            xyz[static_cast<std::size_t>(d)] = PetscRealPart(coord_array[off + d]);
        return PETSC_SUCCESS;
    };

    auto collectFaceVertices = [&](PetscInt face, std::vector<FaultVertex>& vertices) -> PetscErrorCode {
        PetscInt closure_size = 0;
        PetscInt* closure = nullptr;
        PetscErrorCode e = DMPlexGetTransitiveClosure(dm, face, PETSC_TRUE, &closure_size, &closure); CHKERRQ(e);
        vertices.clear();
        for (PetscInt i = 0; i < closure_size; ++i) {
            PetscInt point = closure[2*i];
            PetscInt depth = -1;
            e = DMLabelGetValue(depth_label, point, &depth); CHKERRQ(e);
            if (depth != 0) continue;
            FaultVertex v;
            v.point = point;
            e = loadCoordinates(point, v.xyz); CHKERRQ(e);
            vertices.push_back(v);
        }
        e = DMPlexRestoreTransitiveClosure(dm, face, PETSC_TRUE, &closure_size, &closure); CHKERRQ(e);
        return PETSC_SUCCESS;
    };

    PetscInt cStart = 0, cEnd = 0;
    ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd); CHKERRQ(ierr);
    PetscInt nset = 0;

    for (PetscInt c = cStart; c < cEnd; ++c) {
        DMPolytopeType ct;
        ierr = DMPlexGetCellType(dm, c, &ct); CHKERRQ(ierr);
        if (ct != DM_POLYTOPE_SEG_PRISM_TENSOR &&
            ct != DM_POLYTOPE_TRI_PRISM_TENSOR &&
            ct != DM_POLYTOPE_QUAD_PRISM_TENSOR)
            continue;

        const PetscInt* cone = nullptr;
        PetscInt cone_size = 0;
        ierr = DMPlexGetConeSize(dm, c, &cone_size); CHKERRQ(ierr);
        if (cone_size < 2) continue;
        ierr = DMPlexGetCone(dm, c, &cone); CHKERRQ(ierr);

        std::vector<FaultVertex> neg_vertices;
        ierr = collectFaceVertices(cone[0], neg_vertices); CHKERRQ(ierr);

        for (const auto& vertex : neg_vertices) {
            PetscInt lag_dof = 0, lag_off = 0;
            ierr = PetscSectionGetFieldDof(section, vertex.point, lagrange_field, &lag_dof); CHKERRQ(ierr);
            if (lag_dof <= 0) continue;
            ierr = PetscSectionGetFieldOffset(section, vertex.point, lagrange_field, &lag_off); CHKERRQ(ierr);

            // Determine shear stress: background or nucleation patch
            PetscReal tau = fault_initial_shear_stress_;
            if (fault_nucleation_radius_ > 0.0) {
                // Distance from nucleation center in the fault plane
                // nucleation_center_x is along-strike distance from domain center
                // nucleation_center_z is depth (vertical)
                PetscReal dx = vertex.xyz[0] - (fault_center_[0] >= 0 ? fault_center_[0] : 0.0)
                             - fault_nucleation_center_x_ * strike_dir[0];
                PetscReal dy = vertex.xyz[1] - (fault_center_[1] >= 0 ? fault_center_[1] : 0.0)
                             - fault_nucleation_center_x_ * strike_dir[1];
                PetscReal dz = vertex.xyz[2] - fault_nucleation_center_z_;
                PetscReal dist = PetscSqrtReal(dx*dx + dy*dy + dz*dz);
                if (dist <= fault_nucleation_radius_) {
                    tau = fault_nucleation_shear_stress_;
                }
            }

            // Set lambda = tau * strike_dir + sigma_n * normal_dir
            // sigma_n is compressive negative, lambda convention: compression = negative
            for (PetscInt d = 0; d < PetscMin(lag_dof, dim); ++d) {
                larray[lag_off + d] = tau * strike_dir[d] + fault_initial_normal_stress_ * normal_dir[d];
            }
            nset++;
        }
    }

    ierr = VecRestoreArrayRead(coords, &coord_array); CHKERRQ(ierr);
    ierr = VecRestoreArray(localVec, &larray); CHKERRQ(ierr);
    ierr = DMLocalToGlobal(dm, localVec, INSERT_VALUES, solution); CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(dm, &localVec); CHKERRQ(ierr);

    if (rank == 0) {
        PetscPrintf(comm, "Applied initial fault stress to %d Lagrange DOF sets\n", (int)nset);
    }

    PetscFunctionReturn(PETSC_SUCCESS);
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
    
    ierr = ensureDirectoryExists(comm, rank, output_directory_); CHKERRQ(ierr);
    
    // Set initial conditions for the time stepper
    if (config.enable_elastodynamics) {
        // TSALPHA2 requires both displacement (U) and velocity (V) initial conditions
        Vec velocity;
        ierr = VecDuplicate(solution, &velocity); CHKERRQ(ierr);
        ierr = VecZeroEntries(velocity); CHKERRQ(ierr);
        ierr = TS2SetSolution(ts, solution, velocity); CHKERRQ(ierr);
        ierr = VecDestroy(&velocity); CHKERRQ(ierr);
    } else {
        ierr = TSSetSolution(ts, solution); CHKERRQ(ierr);
    }
    
    // Solve
    ierr = TSSolve(ts, solution); CHKERRQ(ierr);

    // Compute and write derived fields (stress, strain, CFS) if requested
    if (config.output_stress || config.output_cfs)
    {
      // Get material parameters from the first material
      double lambda = 0.0, mu_mat = 0.0, biot = 0.0;
      if (!material_props.empty())
      {
        const auto& mat = material_props[0];
        double E = mat.youngs_modulus;
        double nu = mat.poisson_ratio;
        lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        mu_mat = E / (2.0 * (1.0 + nu));
        biot = mat.biot_coefficient;
      }
      derived_fields_.setMaterialParameters(lambda, mu_mat, biot);

      if (config.output_cfs)
      {
        derived_fields_.setReceiverOrientation(
          config.cfs_receiver_strike, config.cfs_receiver_dip, config.cfs_friction);
      }

      ierr = derived_fields_.compute(dm, solution); CHKERRQ(ierr);

      if (config.output_stress)
      {
        std::string stress_path = output_directory_ + "/stress.h5";
        ierr = derived_fields_.writeStressHDF5(comm, stress_path.c_str(), timestep, current_time);
        CHKERRQ(ierr);
        if (rank == 0)
        {
          PetscPrintf(comm, "Stress field written to %s (%d cells, 6 components)\n",
                      stress_path.c_str(), (int)derived_fields_.numCells());
        }
      }

      if (config.output_cfs)
      {
        std::string cfs_path = output_directory_ + "/cfs.h5";
        ierr = derived_fields_.writeCfsHDF5(comm, cfs_path.c_str(), timestep, current_time);
        CHKERRQ(ierr);
        if (rank == 0)
        {
          PetscPrintf(comm, "CFS field written to %s (%d cells)\n",
                      cfs_path.c_str(), (int)derived_fields_.numCells());
        }
      }
    }

    // Save solution to binary file if requested
    if (!config.save_solution_path.empty()) {
        PetscViewer viewer;
        ierr = PetscViewerBinaryOpen(comm, config.save_solution_path.c_str(),
                                     FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
        ierr = VecView(solution, viewer); CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
        if (rank == 0) {
            PetscPrintf(comm, "Saved solution to %s\n", config.save_solution_path.c_str());
        }
    }
    
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

    if (sim->use_fem_time_residual_) {
        DM dm = sim->dm;

        // DMPlexTSComputeIFunctionFEM expects LOCAL vectors.
        // The TS callback provides GLOBAL vectors, so we must convert.
        Vec locU, locU_t = NULL, locF;
        ierr = DMGetLocalVector(dm, &locU); CHKERRQ(ierr);
        ierr = DMGetLocalVector(dm, &locF); CHKERRQ(ierr);
        ierr = VecZeroEntries(locU); CHKERRQ(ierr);
        ierr = VecZeroEntries(locF); CHKERRQ(ierr);

        // Scatter global solution to local and insert boundary values
        ierr = DMGlobalToLocal(dm, U, INSERT_VALUES, locU); CHKERRQ(ierr);
        ierr = DMPlexInsertBoundaryValues(dm, PETSC_TRUE, locU, t, NULL, NULL, NULL); CHKERRQ(ierr);

        if (U_t) {
            ierr = DMGetLocalVector(dm, &locU_t); CHKERRQ(ierr);
            ierr = VecZeroEntries(locU_t); CHKERRQ(ierr);
            ierr = DMGlobalToLocal(dm, U_t, INSERT_VALUES, locU_t); CHKERRQ(ierr);
        }

        // Compute FEM residual on local vectors
        ierr = DMPlexTSComputeIFunctionFEM(dm, t, locU, locU_t, locF, ctx); CHKERRQ(ierr);

        // Add injection source to the local residual (modifies pressure DOFs)
        if (sim->injection_enabled_ && sim->injection_cell_ >= 0) {
            PetscReal t_local = static_cast<double>(t);
            if (t_local >= sim->injection_start_ && t_local <= sim->injection_end_) {
                ierr = sim->addInjectionToResidual(t, locF); CHKERRQ(ierr);
            }
        }

        // Add explosion moment tensor source to displacement DOFs
        if (sim->explosion_ && sim->explosion_cell_ >= 0 && sim->use_fem_time_residual_) {
            ierr = sim->addExplosionSourceToResidual(t, locF); CHKERRQ(ierr);
        }

        // Apply traction boundary conditions via manual assembly
        if (sim->has_traction_bc_) {
            ierr = sim->addTractionToResidual(t, locF); CHKERRQ(ierr);
        }

        // Apply uniform pressurized-fracture traction on cohesive faces.
        if (sim->hydrofrac_fem_pressurized_mode_ && sim->config.enable_faults) {
            ierr = sim->addFaultPressureToResidual(t, locF); CHKERRQ(ierr);
        }

        // Cohesive fault constraint residual is now handled by PetscDS BdResidual
        // callbacks registered in setupPhysics(). DMPlexTSComputeIFunctionFEM
        // evaluates them on cohesive cells automatically (PETSc 3.25+).
        // The manual Jacobian (addCohesivePenaltyToJacobian) remains in
        // FormJacobian because BdJacobian is not functional in PETSc 3.25.

        // Scatter local residual back to global (ADD_VALUES to accumulate)
        ierr = DMLocalToGlobal(dm, locF, ADD_VALUES, F); CHKERRQ(ierr);

        if (U_t) {
            ierr = DMRestoreLocalVector(dm, &locU_t); CHKERRQ(ierr);
        }
        ierr = DMRestoreLocalVector(dm, &locF); CHKERRQ(ierr);
        ierr = DMRestoreLocalVector(dm, &locU); CHKERRQ(ierr);
    } else {
        // Safe fallback: enforce U_t = 0 (no time evolution) to keep simulations
        // runnable even when residuals aren't configured for the selected fields.
        // F = U_t
        ierr = VecCopy(U_t, F); CHKERRQ(ierr);
    }
    
    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::FormJacobian(TS ts, PetscReal t, Vec U, Vec U_t, 
                                       PetscReal a, Mat J, Mat P, void *ctx) {
    (void)ts;  // Part of PETSc callback interface - accessed via ctx
    
    PetscFunctionBeginUser;
    Simulator *sim = static_cast<Simulator*>(ctx);
    PetscErrorCode ierr;
    
    if (sim->use_fem_time_residual_) {
        DM dm = sim->dm;

        // DMPlexTSComputeIJacobianFEM expects LOCAL vectors.
        Vec locU, locU_t = NULL;
        ierr = DMGetLocalVector(dm, &locU); CHKERRQ(ierr);
        ierr = VecZeroEntries(locU); CHKERRQ(ierr);
        ierr = DMGlobalToLocal(dm, U, INSERT_VALUES, locU); CHKERRQ(ierr);
        ierr = DMPlexInsertBoundaryValues(dm, PETSC_TRUE, locU, t, NULL, NULL, NULL); CHKERRQ(ierr);

        if (U_t) {
            ierr = DMGetLocalVector(dm, &locU_t); CHKERRQ(ierr);
            ierr = VecZeroEntries(locU_t); CHKERRQ(ierr);
            ierr = DMGlobalToLocal(dm, U_t, INSERT_VALUES, locU_t); CHKERRQ(ierr);
        }

        ierr = DMPlexTSComputeIJacobianFEM(dm, t, locU, locU_t, a, J, P, ctx);
        CHKERRQ(ierr);

        // The cohesive interface currently relies on finite-difference
        // coloring for the hybrid-cell coupling terms. Keep the placeholder
        // hook so the explicit interface Jacobian can be added later without
        // changing the call pattern here.
        if (sim->config.enable_faults && !sim->hydrofrac_fem_pressurized_mode_) {
            ierr = MatSetOption(P, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE); CHKERRQ(ierr);
            ierr = sim->addCohesivePenaltyToJacobian(P, locU);
            CHKERRQ(ierr);
            if (J != P) {
                ierr = MatSetOption(J, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE); CHKERRQ(ierr);
                ierr = sim->addCohesivePenaltyToJacobian(J, locU); CHKERRQ(ierr);
            }
        }

        if (U_t) {
            ierr = DMRestoreLocalVector(dm, &locU_t); CHKERRQ(ierr);
        }
        ierr = DMRestoreLocalVector(dm, &locU); CHKERRQ(ierr);
    } else {
        // Jacobian for F = U_t is simply shift * I in PETSc's implicit form.
        ierr = MatZeroEntries(P); CHKERRQ(ierr);
        ierr = MatShift(P, a); CHKERRQ(ierr);
        if (J != P) {
            ierr = MatZeroEntries(J); CHKERRQ(ierr);
            ierr = MatShift(J, a); CHKERRQ(ierr);
        }
    }
    
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

    // Sample seismometers (if enabled)
    if (sim->seismometers_) {
        ierr = sim->seismometers_->sample(static_cast<double>(t), U); CHKERRQ(ierr);
    }

    // Update hydraulic fracture model (if enabled)
    if (sim->hydrofrac_ && sim->injection_enabled_) {
        // Estimate wellbore pressure from injection rate and formation transmissibility
        // For a proper estimate we would extract pressure from the FEM solution at the
        // injection point. Use a simplified analytical estimate instead:
        // p = Q * mu / (4*pi*k*H) * ln(re/rw) + sigma_min
        double k = 10.0 * 9.869233e-16;  // 10 mD in m^2
        double H = sim->hydrofrac_->getHeight();
        double mu_f = 0.001;
        double rw = 0.1;
        double re = 100.0;
        double sigma_min = sim->hydrofrac_->getMinHorizontalStress();
        double wellbore_pressure = sigma_min +
            sim->injection_rate_ / (4.0 * M_PI * k * H / mu_f) * std::log(re / rw);

        // Compute stress intensity factor
        double half_len = sim->hydrofrac_->getLength() / 2.0;
        if (half_len < 1.0) half_len = 1.0;
        double net_pressure = wellbore_pressure - sigma_min;
        double K_I = net_pressure * std::sqrt(M_PI * half_len);

        if (K_I > sim->hydrofrac_->getFractureToughness()) {
            if (!sim->hydrofrac_initiated_) {
                sim->hydrofrac_initiated_ = true;
                if (sim->rank == 0) {
                    PetscPrintf(sim->comm,
                        "FRACTURE INITIATED at t=%.2f, K_I=%.2e > K_Ic=%.2e\n",
                        (double)t, K_I, sim->hydrofrac_->getFractureToughness());
                }
            }
            sim->hydrofrac_->propagate(wellbore_pressure, sim->dt);

            double frac_perm = sim->hydrofrac_->getWidth() * sim->hydrofrac_->getWidth() / 12.0;
            if (sim->rank == 0 && ((int)step % 10 == 0)) {
                PetscPrintf(sim->comm,
                    "FRACTURE: L=%.2f m, W=%.4e m, H=%.1f m, k_f=%.2e m^2\n",
                    sim->hydrofrac_->getLength(), sim->hydrofrac_->getWidth(),
                    sim->hydrofrac_->getHeight(), frac_perm);
            }
        }
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


PetscErrorCode Simulator::writeOutput(int step) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

    ierr = ensureDirectoryExists(comm, rank, output_directory_); CHKERRQ(ierr);

    // Determine output format: HDF5 (default) or VTK
    std::string fmt = config.output_format;
    std::transform(fmt.begin(), fmt.end(), fmt.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

    if (fmt == "HDF5" || fmt == "XDMF") {
        // Build HDF5 output file path
        std::string h5_path = output_directory_ + "/solution.h5";

        PetscViewer viewer;
        if (!output_topology_written_) {
            // First write: create the file and write the DM topology
            ierr = PetscViewerHDF5Open(comm, h5_path.c_str(), FILE_MODE_WRITE, &viewer);
            CHKERRQ(ierr);

            // Write mesh topology once
            ierr = DMView(dm, viewer); CHKERRQ(ierr);

            output_topology_written_ = true;
        } else {
            // Subsequent writes: append to the existing file
            ierr = PetscViewerHDF5Open(comm, h5_path.c_str(), FILE_MODE_APPEND, &viewer);
            CHKERRQ(ierr);
        }

        // Use DMSetOutputSequenceNumber to associate a timestep index and time
        // with this output. DMPlex VecView reads this to label HDF5 datasets.
        ierr = DMSetOutputSequenceNumber(dm, static_cast<PetscInt>(step),
                                         static_cast<PetscReal>(current_time));
        CHKERRQ(ierr);

        // Name the solution vector so HDF5 datasets have meaningful labels
        ierr = PetscObjectSetName(reinterpret_cast<PetscObject>(solution), "solution");
        CHKERRQ(ierr);

        // Write the solution vector
        ierr = VecView(solution, viewer); CHKERRQ(ierr);

        ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

        if (rank == 0 && step % 10 == 0) {
            PetscPrintf(comm, "HDF5 output written at step %d to %s\n", step, h5_path.c_str());
        }
    } else if (fmt == "VTK") {
        // Write VTK output (.vtu file per step)
        char vtu_path[PETSC_MAX_PATH_LEN];
        ierr = PetscSNPrintf(vtu_path, sizeof(vtu_path), "%s/solution_%06d.vtu",
                             output_directory_.c_str(), step); CHKERRQ(ierr);

        PetscViewer viewer;
        ierr = PetscViewerVTKOpen(comm, vtu_path, FILE_MODE_WRITE, &viewer);
        CHKERRQ(ierr);

        ierr = PetscObjectSetName(reinterpret_cast<PetscObject>(solution), "solution");
        CHKERRQ(ierr);
        ierr = VecView(solution, viewer); CHKERRQ(ierr);

        ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

        if (rank == 0 && step % 10 == 0) {
            PetscPrintf(comm, "VTK output written at step %d to %s\n", step, vtu_path);
        }
    } else {
        // Unknown format or ECLIPSE -- print stub message
        if (rank == 0 && step % 10 == 0) {
            PetscPrintf(comm, "Output at step %d (format=%s not implemented)\n",
                        step, config.output_format.c_str());
        }
    }

    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::writeCheckpoint(int step) {
    (void)step;
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

    // Write seismometer outputs (rank 0)
    if (seismometers_) {
        PetscErrorCode ierr = seismometers_->finalizeAndWrite(); CHKERRQ(ierr);
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
    (void)type1;
    (void)type2;
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

PetscErrorCode Simulator::labelBoundaries() {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

    if (rank == 0) {
        PetscPrintf(comm, "Labeling boundary faces...\n");
    }

    // Get mesh dimension
    PetscInt dim;
    ierr = DMGetDimension(dm, &dim); CHKERRQ(ierr);

    // Get bounding box
    Vec coords;
    ierr = DMGetCoordinatesLocal(dm, &coords); CHKERRQ(ierr);
    const PetscScalar *coordArray;
    ierr = VecGetArrayRead(coords, &coordArray); CHKERRQ(ierr);

    PetscInt coordSize;
    ierr = VecGetLocalSize(coords, &coordSize); CHKERRQ(ierr);
    PetscInt nCoords = coordSize / dim;

    // Compute bounding box
    double xmin = 1e30, xmax = -1e30;
    double ymin = 1e30, ymax = -1e30;
    double zmin = 1e30, zmax = -1e30;

    for (PetscInt i = 0; i < nCoords; i++) {
        double x = PetscRealPart(coordArray[i * dim + 0]);
        double y = (dim >= 2) ? PetscRealPart(coordArray[i * dim + 1]) : 0.0;
        double z = (dim >= 3) ? PetscRealPart(coordArray[i * dim + 2]) : 0.0;

        if (x < xmin) xmin = x;
        if (x > xmax) xmax = x;
        if (y < ymin) ymin = y;
        if (y > ymax) ymax = y;
        if (z < zmin) zmin = z;
        if (z > zmax) zmax = z;
    }

    ierr = VecRestoreArrayRead(coords, &coordArray); CHKERRQ(ierr);

    // Reduce across all ranks
    double local_min[3] = {xmin, ymin, zmin};
    double local_max[3] = {xmax, ymax, zmax};
    double global_min[3], global_max[3];
    ierr = MPI_Allreduce(local_min, global_min, 3, MPI_DOUBLE, MPI_MIN, comm); CHKERRQ(ierr);
    ierr = MPI_Allreduce(local_max, global_max, 3, MPI_DOUBLE, MPI_MAX, comm); CHKERRQ(ierr);

    xmin = global_min[0]; ymin = global_min[1]; zmin = global_min[2];
    xmax = global_max[0]; ymax = global_max[1]; zmax = global_max[2];

    // Tolerance for boundary detection (1% of domain size)
    double eps_x = 0.01 * (xmax - xmin);
    double eps_y = 0.01 * (ymax - ymin);
    double eps_z = 0.01 * (zmax - zmin);

    if (rank == 0) {
        PetscPrintf(comm, "  Bounding box: [%.3f, %.3f] x [%.3f, %.3f] x [%.3f, %.3f]\n",
                    xmin, xmax, ymin, ymax, zmin, zmax);
    }

    // Create DMLabels for each boundary face
    const char *label_names[6] = {
        "boundary_x_min", "boundary_x_max",
        "boundary_y_min", "boundary_y_max",
        "boundary_z_min", "boundary_z_max"
    };

    DMLabel labels[6];
    for (int i = 0; i < 6; i++) {
        ierr = DMCreateLabel(dm, label_names[i]); CHKERRQ(ierr);
        ierr = DMGetLabel(dm, label_names[i], &labels[i]); CHKERRQ(ierr);
    }

    // Get face stratum (height 1)
    PetscInt fStart, fEnd;
    ierr = DMPlexGetHeightStratum(dm, 1, &fStart, &fEnd); CHKERRQ(ierr);

    // Iterate over all faces
    PetscInt nLabeled[6] = {0, 0, 0, 0, 0, 0};
    for (PetscInt face = fStart; face < fEnd; face++) {
        // Check if this is a boundary face (support size == 1)
        PetscInt supportSize;
        ierr = DMPlexGetSupportSize(dm, face, &supportSize); CHKERRQ(ierr);
        if (supportSize != 1) continue; // Interior face

        // Get face centroid
        PetscReal centroid[3] = {0.0, 0.0, 0.0};
        PetscReal normal[3] = {0.0, 0.0, 0.0};
        ierr = DMPlexComputeCellGeometryFVM(dm, face, nullptr, centroid, normal); CHKERRQ(ierr);

        // Assign to appropriate label based on centroid
        if (dim >= 1 && PetscAbsReal(centroid[0] - xmin) < eps_x) {
            ierr = DMLabelSetValue(labels[0], face, 1); CHKERRQ(ierr);
            nLabeled[0]++;
        } else if (dim >= 1 && PetscAbsReal(centroid[0] - xmax) < eps_x) {
            ierr = DMLabelSetValue(labels[1], face, 1); CHKERRQ(ierr);
            nLabeled[1]++;
        }

        if (dim >= 2 && PetscAbsReal(centroid[1] - ymin) < eps_y) {
            ierr = DMLabelSetValue(labels[2], face, 1); CHKERRQ(ierr);
            nLabeled[2]++;
        } else if (dim >= 2 && PetscAbsReal(centroid[1] - ymax) < eps_y) {
            ierr = DMLabelSetValue(labels[3], face, 1); CHKERRQ(ierr);
            nLabeled[3]++;
        }

        if (dim >= 3 && PetscAbsReal(centroid[2] - zmin) < eps_z) {
            ierr = DMLabelSetValue(labels[4], face, 1); CHKERRQ(ierr);
            nLabeled[4]++;
        } else if (dim >= 3 && PetscAbsReal(centroid[2] - zmax) < eps_z) {
            ierr = DMLabelSetValue(labels[5], face, 1); CHKERRQ(ierr);
            nLabeled[5]++;
        }
    }

    // Call DMPlexLabelComplete on each label to add vertices and edges in closure
    for (int i = 0; i < 6; i++) {
        ierr = DMPlexLabelComplete(dm, labels[i]); CHKERRQ(ierr);
    }

    // Report labeled faces
    if (rank == 0) {
        for (int i = 0; i < 6; i++) {
            PetscInt global_count;
            ierr = MPI_Reduce(&nLabeled[i], &global_count, 1, MPIU_INT, MPI_SUM, 0, comm); CHKERRQ(ierr);
            PetscPrintf(comm, "  %s: %d faces\n", label_names[i], global_count);
        }
    } else {
        for (int i = 0; i < 6; i++) {
            ierr = MPI_Reduce(&nLabeled[i], nullptr, 1, MPIU_INT, MPI_SUM, 0, comm); CHKERRQ(ierr);
        }
    }

    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::setupBoundaryConditions() {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

    if (rank == 0) {
        PetscPrintf(comm, "Setting up boundary conditions...\n");
    }

    // Determine field indices based on configuration
    // NOTE: This is called BEFORE DMCreateDS, so we cannot use DMGetDS yet
    PetscInt pressure_field = -1;
    PetscInt displacement_field = -1;

    if (config.solid_model == SolidModelType::POROELASTIC &&
        config.fluid_model == FluidModelType::SINGLE_COMPONENT) {
        // Poroelasticity: field 0 = pressure, field 1 = displacement
        pressure_field = 0;
        displacement_field = 1;
    } else if (config.fluid_model == FluidModelType::SINGLE_COMPONENT) {
        // Fluid flow only: field 0 = pressure
        pressure_field = 0;
    } else if (config.enable_geomechanics || config.enable_elastodynamics) {
        // Geomechanics only: displacement field is 0 for NONE fluid model, else 1
        displacement_field = 0;
        if (config.fluid_model == FluidModelType::SINGLE_COMPONENT) {
            displacement_field = 1;
        } else if (config.fluid_model == FluidModelType::BLACK_OIL) {
            displacement_field = 3;
        } else if (config.fluid_model == FluidModelType::COMPOSITIONAL) {
            displacement_field = 4;
        }
    }

    // Apply boundary conditions based on physics
    DMLabel label_xmin, label_xmax, label_ymin, label_ymax, label_zmin, label_zmax;
    ierr = DMGetLabel(dm, "boundary_x_min", &label_xmin); CHKERRQ(ierr);
    ierr = DMGetLabel(dm, "boundary_x_max", &label_xmax); CHKERRQ(ierr);
    ierr = DMGetLabel(dm, "boundary_y_min", &label_ymin); CHKERRQ(ierr);
    ierr = DMGetLabel(dm, "boundary_y_max", &label_ymax); CHKERRQ(ierr);
    ierr = DMGetLabel(dm, "boundary_z_min", &label_zmin); CHKERRQ(ierr);
    ierr = DMGetLabel(dm, "boundary_z_max", &label_zmax); CHKERRQ(ierr);

    if (config.solid_model == SolidModelType::POROELASTIC &&
        config.fluid_model == FluidModelType::SINGLE_COMPONENT) {
        // POROELASTIC case: field 0 = pressure, field 1 = displacement
        PetscInt pfield = 0;
        PetscInt ufield = 1;
        PetscInt label_value = 1;

        // Bottom: fix all displacement components (impermeable = natural BC on pressure)
        if (label_zmin) {
            PetscInt comps[3] = {0, 1, 2};
            ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "fixed_bottom",
                label_zmin, 1, &label_value, ufield, 3, comps,
                (void (*)(void))bc_zero, nullptr, nullptr, nullptr); CHKERRQ(ierr);
            if (rank == 0) {
                PetscPrintf(comm, "  Applied BC: fixed bottom (z_min), field %d, all components = 0\n", ufield);
            }
        }
        // Top: drained (p = 0)
        if (label_zmax) {
            PetscInt comp = 0;
            ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "drained_top",
                label_zmax, 1, &label_value, pfield, 1, &comp,
                (void (*)(void))bc_zero, nullptr, nullptr, nullptr); CHKERRQ(ierr);
            if (rank == 0) {
                PetscPrintf(comm, "  Applied BC: drained top (z_max), field %d, pressure = 0\n", pfield);
            }
        }
        // Top: applied compression (u_z = -0.001)
        if (label_zmax) {
            PetscInt comps[3] = {0, 1, 2};
            ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "compression_top",
                label_zmax, 1, &label_value, ufield, 3, comps,
                (void (*)(void))bc_compression, nullptr, nullptr, nullptr); CHKERRQ(ierr);
            if (rank == 0) {
                PetscPrintf(comm, "  Applied BC: compression top (z_max), field %d, u_z = -0.001\n", ufield);
            }
        }
        // Lateral rollers
        if (label_xmin) {
            PetscInt comp = 0;
            ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "roller_xmin",
                label_xmin, 1, &label_value, ufield, 1, &comp,
                (void (*)(void))bc_zero, nullptr, nullptr, nullptr); CHKERRQ(ierr);
            if (rank == 0) {
                PetscPrintf(comm, "  Applied BC: roller x_min, field %d, u_x = 0\n", ufield);
            }
        }
        if (label_xmax) {
            PetscInt comp = 0;
            ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "roller_xmax",
                label_xmax, 1, &label_value, ufield, 1, &comp,
                (void (*)(void))bc_zero, nullptr, nullptr, nullptr); CHKERRQ(ierr);
            if (rank == 0) {
                PetscPrintf(comm, "  Applied BC: roller x_max, field %d, u_x = 0\n", ufield);
            }
        }
        if (label_ymin) {
            PetscInt comp = 1;
            ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "roller_ymin",
                label_ymin, 1, &label_value, ufield, 1, &comp,
                (void (*)(void))bc_zero, nullptr, nullptr, nullptr); CHKERRQ(ierr);
            if (rank == 0) {
                PetscPrintf(comm, "  Applied BC: roller y_min, field %d, u_y = 0\n", ufield);
            }
        }
        if (label_ymax) {
            PetscInt comp = 1;
            ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "roller_ymax",
                label_ymax, 1, &label_value, ufield, 1, &comp,
                (void (*)(void))bc_zero, nullptr, nullptr, nullptr); CHKERRQ(ierr);
            if (rank == 0) {
                PetscPrintf(comm, "  Applied BC: roller y_max, field %d, u_y = 0\n", ufield);
            }
        }
    } else if (displacement_field >= 0) {
        // Pure elastostatics: field 0 = displacement
        PetscInt field = displacement_field;
        PetscInt comps_all[3] = {0, 1, 2};
        PetscInt label_value = 1;

        // Check whether per-face BCs are configured
        bool use_per_face = false;
        for (int fi = 0; fi < 6; ++fi) {
            if (config.face_bc[fi].configured) {
                use_per_face = true;
                break;
            }
        }

        if (use_per_face) {
            // Per-face boundary conditions (new system)
            DMLabel face_labels[6] = {label_xmin, label_xmax, label_ymin, label_ymax, label_zmin, label_zmax};
            static const char *face_names[6] = {
                "x_min", "x_max", "y_min", "y_max", "z_min", "z_max"
            };
            // Normal component index for roller BCs: x faces constrain comp 0, y faces comp 1, z faces comp 2
            PetscInt roller_comp[6] = {0, 0, 1, 1, 2, 2};

            for (int fi = 0; fi < 6; ++fi) {
                if (!config.face_bc[fi].configured || !face_labels[fi]) continue;
                const auto &fbc = config.face_bc[fi];
                std::string bc_name = std::string("bc_") + face_names[fi];

                if (fbc.type == "fixed") {
                    ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, bc_name.c_str(),
                        face_labels[fi], 1, &label_value, field, 3, comps_all,
                        (void (*)(void))bc_zero, nullptr, nullptr, nullptr); CHKERRQ(ierr);
                    if (rank == 0) {
                        PetscPrintf(comm, "  Applied BC: fixed %s, field %d, all components = 0\n",
                                    face_names[fi], field);
                    }
                } else if (fbc.type == "roller") {
                    PetscInt comp = roller_comp[fi];
                    ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, bc_name.c_str(),
                        face_labels[fi], 1, &label_value, field, 1, &comp,
                        (void (*)(void))bc_zero, nullptr, nullptr, nullptr); CHKERRQ(ierr);
                    if (rank == 0) {
                        PetscPrintf(comm, "  Applied BC: roller %s, field %d, u_%d = 0\n",
                                    face_names[fi], field, (int)comp);
                    }
                } else if (fbc.type == "roller_y") {
                    // Plane-strain roller: constrain u_y = 0
                    PetscInt comp = 1;
                    ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, bc_name.c_str(),
                        face_labels[fi], 1, &label_value, field, 1, &comp,
                        (void (*)(void))bc_zero, nullptr, nullptr, nullptr); CHKERRQ(ierr);
                    if (rank == 0) {
                        PetscPrintf(comm, "  Applied BC: roller_y %s, field %d, u_y = 0\n",
                                    face_names[fi], field);
                    }
                } else if (fbc.type == "compression") {
                    ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, bc_name.c_str(),
                        face_labels[fi], 1, &label_value, field, 3, comps_all,
                        (void (*)(void))bc_compression, nullptr, nullptr, nullptr); CHKERRQ(ierr);
                    if (rank == 0) {
                        PetscPrintf(comm, "  Applied BC: compression %s, field %d\n",
                                    face_names[fi], field);
                    }
                } else if (fbc.type == "dirichlet") {
                    // General Dirichlet: constrain specified components to specified values
                    PetscInt ncomps = static_cast<PetscInt>(fbc.components.size());
                    std::vector<PetscInt> pcomps(fbc.components.begin(), fbc.components.end());
                    bool all_zero = true;
                    for (const auto &v : fbc.values) {
                        if (std::abs(v) > 1.0e-40) { all_zero = false; break; }
                    }
                    if (all_zero) {
                        ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, bc_name.c_str(),
                            face_labels[fi], 1, &label_value, field, ncomps, pcomps.data(),
                            (void (*)(void))bc_zero, nullptr, nullptr, nullptr); CHKERRQ(ierr);
                    } else {
                        // Store values in static storage so callback can read them
                        // Use face index to avoid collisions
                        static double dirichlet_vals[6][3] = {};
                        for (int d = 0; d < 3; ++d) {
                            dirichlet_vals[fi][d] = (d < static_cast<int>(fbc.values.size()))
                                                     ? fbc.values[static_cast<std::size_t>(d)] : 0.0;
                        }
                        ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, bc_name.c_str(),
                            face_labels[fi], 1, &label_value, field, ncomps, pcomps.data(),
                            (void (*)(void))bc_dirichlet_values, nullptr,
                            dirichlet_vals[fi], nullptr); CHKERRQ(ierr);
                    }
                    if (rank == 0) {
                        PetscPrintf(comm, "  Applied BC: dirichlet %s, field %d, %d components\n",
                                    face_names[fi], field, (int)ncomps);
                    }
                } else if (fbc.type == "traction") {
                    // Natural BC: register DM_BC_NATURAL so PETSc does not constrain DOFs.
                    // The traction residual contribution is registered later via PetscWeakForm.
                    ierr = DMAddBoundary(dm, DM_BC_NATURAL, bc_name.c_str(),
                        face_labels[fi], 1, &label_value, field, 0, NULL,
                        NULL, NULL, NULL, NULL); CHKERRQ(ierr);
                    if (rank == 0) {
                        PetscPrintf(comm, "  Applied BC: traction %s, field %d, t=(%.2e, %.2e, %.2e)\n",
                                    face_names[fi], field,
                                    fbc.traction[0], fbc.traction[1], fbc.traction[2]);
                    }
                }
                // "free" type: no Dirichlet constraint, natural zero-traction (do nothing)
            }
        } else {
            // Legacy boundary conditions (bottom/sides/top strings)
            // Bottom BC
            if (label_zmin && config.bottom_bc == "fixed") {
                ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "fixed_bottom", label_zmin, 1, &label_value,
                                     field, 3, comps_all, (void (*)(void))bc_zero, nullptr, nullptr, nullptr);
                CHKERRQ(ierr);
                if (rank == 0) {
                    PetscPrintf(comm, "  Applied BC: fixed bottom (z_min), field %d, all components = 0\n", field);
                }
            }

            // Top BC
            if (label_zmax && config.top_bc == "compression") {
                ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "compression_top", label_zmax, 1, &label_value,
                                     field, 3, comps_all, (void (*)(void))bc_compression, nullptr, nullptr, nullptr);
                CHKERRQ(ierr);
                if (rank == 0) {
                    PetscPrintf(comm, "  Applied BC: compression top (z_max), field %d, u_z = -0.001\n", field);
                }
            }
            // top_bc == "free" means no Dirichlet BC on top (natural zero-traction)

            // Side BCs: roller constrains normal displacement to zero
            if (config.side_bc == "roller") {
                if (label_xmin) {
                    PetscInt comp = 0;
                    ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "roller_xmin", label_xmin, 1, &label_value,
                                         field, 1, &comp, (void (*)(void))bc_zero, nullptr, nullptr, nullptr);
                    CHKERRQ(ierr);
                    if (rank == 0) {
                        PetscPrintf(comm, "  Applied BC: roller x_min, field %d, u_x = 0\n", field);
                    }
                }
                if (label_xmax) {
                    PetscInt comp = 0;
                    ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "roller_xmax", label_xmax, 1, &label_value,
                                         field, 1, &comp, (void (*)(void))bc_zero, nullptr, nullptr, nullptr);
                    CHKERRQ(ierr);
                    if (rank == 0) {
                        PetscPrintf(comm, "  Applied BC: roller x_max, field %d, u_x = 0\n", field);
                    }
                }
                if (label_ymin) {
                    PetscInt comp = 1;
                    ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "roller_ymin", label_ymin, 1, &label_value,
                                         field, 1, &comp, (void (*)(void))bc_zero, nullptr, nullptr, nullptr);
                    CHKERRQ(ierr);
                    if (rank == 0) {
                        PetscPrintf(comm, "  Applied BC: roller y_min, field %d, u_y = 0\n", field);
                    }
                }
                if (label_ymax) {
                    PetscInt comp = 1;
                    ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "roller_ymax", label_ymax, 1, &label_value,
                                         field, 1, &comp, (void (*)(void))bc_zero, nullptr, nullptr, nullptr);
                    CHKERRQ(ierr);
                    if (rank == 0) {
                        PetscPrintf(comm, "  Applied BC: roller y_max, field %d, u_y = 0\n", field);
                    }
                }
            }

            // Side BCs: fixed constrains all displacement components to zero
            if (config.side_bc == "fixed") {
                if (label_xmin) {
                    ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "fixed_xmin", label_xmin, 1, &label_value,
                                         field, 3, comps_all, (void (*)(void))bc_zero, nullptr, nullptr, nullptr);
                    CHKERRQ(ierr);
                    if (rank == 0) {
                        PetscPrintf(comm, "  Applied BC: fixed x_min, field %d, all components = 0\n", field);
                    }
                }
                if (label_xmax) {
                    ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "fixed_xmax", label_xmax, 1, &label_value,
                                         field, 3, comps_all, (void (*)(void))bc_zero, nullptr, nullptr, nullptr);
                    CHKERRQ(ierr);
                    if (rank == 0) {
                        PetscPrintf(comm, "  Applied BC: fixed x_max, field %d, all components = 0\n", field);
                    }
                }
                if (label_ymin) {
                    ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "fixed_ymin", label_ymin, 1, &label_value,
                                         field, 3, comps_all, (void (*)(void))bc_zero, nullptr, nullptr, nullptr);
                    CHKERRQ(ierr);
                    if (rank == 0) {
                        PetscPrintf(comm, "  Applied BC: fixed y_min, field %d, all components = 0\n", field);
                    }
                }
                if (label_ymax) {
                    ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "fixed_ymax", label_ymax, 1, &label_value,
                                         field, 3, comps_all, (void (*)(void))bc_zero, nullptr, nullptr, nullptr);
                    CHKERRQ(ierr);
                    if (rank == 0) {
                        PetscPrintf(comm, "  Applied BC: fixed y_max, field %d, all components = 0\n", field);
                    }
                }
            }

            // Top BC: fixed constrains all displacement components to zero
            if (config.top_bc == "fixed") {
                if (label_zmax) {
                    ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "fixed_top", label_zmax, 1, &label_value,
                                         field, 3, comps_all, (void (*)(void))bc_zero, nullptr, nullptr, nullptr);
                    CHKERRQ(ierr);
                    if (rank == 0) {
                        PetscPrintf(comm, "  Applied BC: fixed top (z_max), field %d, all components = 0\n", field);
                    }
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // Standalone fluid flow: Dirichlet pressure BCs on the pressure field
    // -------------------------------------------------------------------------
    if (pressure_field >= 0 && displacement_field < 0) {
        // Pure fluid flow (no geomechanics) -- apply pressure Dirichlet BCs
        PetscInt pfield = pressure_field;
        PetscInt label_value = 1;
        DMLabel face_labels_p[6] = {label_xmin, label_xmax, label_ymin, label_ymax, label_zmin, label_zmax};
        static const char *face_names_p[6] = {
            "x_min", "x_max", "y_min", "y_max", "z_min", "z_max"
        };

        // Static storage for pressure values indexed by face (up to 6 faces).
        // The BC callback reads from this array via the context pointer.
        static double pressure_bc_vals[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        for (int fi = 0; fi < 6; ++fi) {
            if (!config.face_bc[fi].configured || !face_labels_p[fi]) continue;
            const auto &fbc = config.face_bc[fi];

            if (fbc.type == "dirichlet_pressure") {
                pressure_bc_vals[fi] = fbc.pressure;

                std::string bc_name = std::string("pressure_") + face_names_p[fi];
                PetscInt comp = 0;
                ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, bc_name.c_str(),
                    face_labels_p[fi], 1, &label_value, pfield, 1, &comp,
                    (void (*)(void))bc_pressure_dirichlet, nullptr,
                    &pressure_bc_vals[fi], nullptr); CHKERRQ(ierr);
                if (rank == 0) {
                    PetscPrintf(comm, "  Applied BC: dirichlet_pressure %s, field %d, P = %.4e Pa\n",
                                face_names_p[fi], pfield, fbc.pressure);
                }
            }
        }
    }

    // Register natural (Neumann) BCs on fault faces for cohesive traction assembly
    if (config.enable_faults && displacement_field >= 0) {
        DMLabel fault_label = nullptr;
        ierr = DMGetLabel(dm, "fault", &fault_label); CHKERRQ(ierr);
        if (fault_label) {
            PetscInt label_value = 1;
            PetscInt lagrange_field_idx = displacement_field + 1;

            ierr = DMAddBoundary(dm, DM_BC_NATURAL, "fault_traction",
                fault_label, 1, &label_value, displacement_field, 0, NULL,
                NULL, NULL, NULL, NULL); CHKERRQ(ierr);

            ierr = DMAddBoundary(dm, DM_BC_NATURAL, "fault_constraint",
                fault_label, 1, &label_value, lagrange_field_idx, 0, NULL,
                NULL, NULL, NULL, NULL); CHKERRQ(ierr);

            if (rank == 0) {
                // Diagnostic: dump the fault label stratification
                IS       stratum_is;
                PetscInt stratum_size;
                PetscInt depth_val;
                DMLabel  depth_label;
                ierr = DMPlexGetDepthLabel(dm, &depth_label); CHKERRQ(ierr);
                ierr = DMLabelGetStratumIS(fault_label, label_value, &stratum_is); CHKERRQ(ierr);
                if (stratum_is) {
                    ierr = ISGetSize(stratum_is, &stratum_size); CHKERRQ(ierr);
                    PetscPrintf(comm, "  Fault label (value=%d): %d points total\n",
                        (int)label_value, (int)stratum_size);

                    const PetscInt *indices;
                    ierr = ISGetIndices(stratum_is, &indices); CHKERRQ(ierr);
                    PetscInt n_cells = 0, n_faces = 0, n_edges = 0, n_vertices = 0;
                    PetscInt dim_mesh;
                    ierr = DMGetDimension(dm, &dim_mesh); CHKERRQ(ierr);
                    for (PetscInt i = 0; i < stratum_size; ++i)
                    {
                        PetscInt pt_depth;
                        ierr = DMLabelGetValue(depth_label, indices[i], &pt_depth); CHKERRQ(ierr);
                        if (pt_depth == dim_mesh) n_cells++;
                        else if (pt_depth == dim_mesh - 1) n_faces++;
                        else if (pt_depth == 1) n_edges++;
                        else if (pt_depth == 0) n_vertices++;
                    }
                    PetscPrintf(comm, "    cells=%d, faces=%d, edges=%d, vertices=%d\n",
                        (int)n_cells, (int)n_faces, (int)n_edges, (int)n_vertices);
                    ierr = ISRestoreIndices(stratum_is, &indices); CHKERRQ(ierr);
                    ierr = ISDestroy(&stratum_is); CHKERRQ(ierr);
                } else {
                    PetscPrintf(comm, "  Fault label (value=%d): no points!\n",
                        (int)label_value);
                }

                // Also check other label values
                PetscInt num_strata;
                ierr = DMLabelGetNumValues(fault_label, &num_strata); CHKERRQ(ierr);
                PetscPrintf(comm, "  Fault label has %d strata total\n", (int)num_strata);
                IS val_is;
                ierr = DMLabelGetValueIS(fault_label, &val_is); CHKERRQ(ierr);
                const PetscInt *vals;
                ierr = ISGetIndices(val_is, &vals); CHKERRQ(ierr);
                for (PetscInt s = 0; s < num_strata; ++s)
                {
                    PetscInt ssize;
                    ierr = DMLabelGetStratumSize(fault_label, vals[s], &ssize); CHKERRQ(ierr);
                    PetscPrintf(comm, "    stratum value=%d: %d points\n",
                        (int)vals[s], (int)ssize);
                }
                ierr = ISRestoreIndices(val_is, &vals); CHKERRQ(ierr);
                ierr = ISDestroy(&val_is); CHKERRQ(ierr);
            }

            if (rank == 0) {
                PetscPrintf(comm,
                    "  Applied BC: fault natural BCs for traction (field %d) and constraint (field %d)\n",
                    displacement_field, lagrange_field_idx);
            }
        }
    }

    // For fluid flow only (no geomechanics):
    // - Fixed pressure on x_min: p = 0 (reference)
    // Skip this legacy BC if per-face dirichlet_pressure BCs are configured.
    {
        bool has_dirichlet_pressure_bc = false;
        for (int fi = 0; fi < 6; ++fi) {
            if (config.face_bc[fi].configured && config.face_bc[fi].type == "dirichlet_pressure") {
                has_dirichlet_pressure_bc = true;
                break;
            }
        }
        if (pressure_field >= 0 && !config.enable_geomechanics &&
            config.solid_model != SolidModelType::POROELASTIC &&
            !has_dirichlet_pressure_bc) {
            PetscInt field = pressure_field;
            PetscInt comp = 0;
            PetscInt label_value = 1;

            if (label_xmin) {
                ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "fixed_pressure", label_xmin, 1, &label_value,
                                     field, 1, &comp, (void (*)(void))bc_drained, nullptr, nullptr, nullptr);
                CHKERRQ(ierr);
                if (rank == 0) {
                    PetscPrintf(comm, "  Applied BC: fixed pressure (x_min), field %d, pressure = 0\n", field);
                }
            }
        }
    }

    // Register absorbing boundary conditions (DM_BC_NATURAL)
    if (config.absorbing_bc_enabled && displacement_field >= 0) {
        PetscInt field = displacement_field;
        PetscInt label_value = 1;

        struct AbsorbingFace {
            bool enabled;
            DMLabel label;
            const char* name;
        };
        AbsorbingFace faces[6] = {
            {config.absorbing_bc_x_min, label_xmin, "absorbing_x_min"},
            {config.absorbing_bc_x_max, label_xmax, "absorbing_x_max"},
            {config.absorbing_bc_y_min, label_ymin, "absorbing_y_min"},
            {config.absorbing_bc_y_max, label_ymax, "absorbing_y_max"},
            {config.absorbing_bc_z_min, label_zmin, "absorbing_z_min"},
            {config.absorbing_bc_z_max, label_zmax, "absorbing_z_max"},
        };

        for (int i = 0; i < 6; ++i) {
            if (faces[i].enabled && faces[i].label) {
                ierr = DMAddBoundary(dm, DM_BC_NATURAL, faces[i].name,
                    faces[i].label, 1, &label_value, field, 0, NULL,
                    NULL, NULL, NULL, NULL); CHKERRQ(ierr);
                if (rank == 0) {
                    PetscPrintf(comm, "  Applied BC: %s (absorbing), field %d\n",
                        faces[i].name, field);
                }
            }
        }
    }

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
