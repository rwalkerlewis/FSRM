#include "core/Simulator.hpp"
#include "core/PhysicsModuleRegistry.hpp"
#include "numerics/PetscFEFluidFlow.hpp"
#include "numerics/PetscFEElasticity.hpp"
#include "numerics/PetscFEPoroelasticity.hpp"
#include "io/Visualization.hpp"
#include "core/ConfigReader.hpp"
#include "numerics/ImplicitExplicitTransition.hpp"
#include "domain/explosion/ExplosionImpactPhysics.hpp"
#include "domain/seismic/SeismometerNetwork.hpp"

// GPU kernel support — conditionally include GPU headers when built with CUDA
#ifdef USE_CUDA
#include "gpu/GPUManager.hpp"
#include "physics/PhysicsKernel_GPU.hpp"
#include "physics/PhysicsKernel_GPU_Extended.hpp"
#endif

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
    
    switch (config.fluid_model) {
        case FluidModelType::NONE:
            // No fluid fields
            break;
        case FluidModelType::SINGLE_COMPONENT:
            // Single pressure field
            ierr = PetscFECreateLagrange(comm, 3, 1, PETSC_FALSE, 1, -1, &fe); CHKERRQ(ierr);
            ierr = PetscObjectSetName((PetscObject)fe, "pressure_"); CHKERRQ(ierr);
            ierr = DMAddField(dm, nullptr, (PetscObject)fe); CHKERRQ(ierr);
            fe_fields.push_back(fe);
            break;
            
        case FluidModelType::BLACK_OIL:
            // Pressure + saturations
            ierr = PetscFECreateLagrange(comm, 3, 1, PETSC_FALSE, 1, -1, &fe); CHKERRQ(ierr);
            ierr = PetscObjectSetName((PetscObject)fe, "pressure_"); CHKERRQ(ierr);
            ierr = DMAddField(dm, nullptr, (PetscObject)fe); CHKERRQ(ierr);
            fe_fields.push_back(fe);
            
            ierr = PetscFECreateLagrange(comm, 3, 1, PETSC_FALSE, 1, -1, &fe); CHKERRQ(ierr);
            ierr = PetscObjectSetName((PetscObject)fe, "saturation_w_"); CHKERRQ(ierr);
            ierr = DMAddField(dm, nullptr, (PetscObject)fe); CHKERRQ(ierr);
            fe_fields.push_back(fe);
            
            ierr = PetscFECreateLagrange(comm, 3, 1, PETSC_FALSE, 1, -1, &fe); CHKERRQ(ierr);
            ierr = PetscObjectSetName((PetscObject)fe, "saturation_g_"); CHKERRQ(ierr);
            ierr = DMAddField(dm, nullptr, (PetscObject)fe); CHKERRQ(ierr);
            fe_fields.push_back(fe);
            break;
            
        case FluidModelType::COMPOSITIONAL:
            // Pressure + compositions
            ierr = PetscFECreateLagrange(comm, 3, 1, PETSC_FALSE, 1, -1, &fe); CHKERRQ(ierr);
            ierr = PetscObjectSetName((PetscObject)fe, "pressure_"); CHKERRQ(ierr);
            ierr = DMAddField(dm, nullptr, (PetscObject)fe); CHKERRQ(ierr);
            fe_fields.push_back(fe);
            
            // Add composition fields (number depends on components)
            for (int i = 0; i < 3; ++i) {  // Example: 3 components
                char name[256];
                snprintf(name, sizeof(name), "composition_%d_", i);
                ierr = PetscFECreateLagrange(comm, 3, 1, PETSC_FALSE, 1, -1, &fe); CHKERRQ(ierr);
                ierr = PetscObjectSetName((PetscObject)fe, name); CHKERRQ(ierr);
                ierr = DMAddField(dm, nullptr, (PetscObject)fe); CHKERRQ(ierr);
                fe_fields.push_back(fe);
            }
            break;
            
        default:
            SETERRQ(comm, PETSC_ERR_ARG_WRONG, "Unknown fluid model type");
    }
    
    // Add geomechanics field if enabled
    if (config.enable_geomechanics) {
        ierr = PetscFECreateLagrange(comm, 3, 3, PETSC_FALSE, 1, -1, &fe); CHKERRQ(ierr);
        ierr = PetscObjectSetName((PetscObject)fe, "displacement_"); CHKERRQ(ierr);
        ierr = DMAddField(dm, nullptr, (PetscObject)fe); CHKERRQ(ierr);
        fe_fields.push_back(fe);
    }

    // Add Lagrange multiplier field for cohesive faults if enabled
    if (config.enable_faults && cohesive_kernel_) {
        PetscFE fe_lagrange;
        // Lagrange multiplier: 3-component vector (traction) on cohesive cells
        ierr = PetscFECreateLagrange(comm, 3, 3, PETSC_FALSE, 1, -1, &fe_lagrange); CHKERRQ(ierr);
        ierr = PetscObjectSetName((PetscObject)fe_lagrange, "lagrange_"); CHKERRQ(ierr);
        ierr = DMAddField(dm, nullptr, (PetscObject)fe_lagrange); CHKERRQ(ierr);
        fe_fields.push_back(fe_lagrange);
        if (rank == 0) {
            PetscPrintf(comm, "Added Lagrange multiplier field for fault traction\n");
        }
    }

    // Add thermal field if enabled
    if (config.enable_thermal) {
        ierr = PetscFECreateLagrange(comm, 3, 1, PETSC_FALSE, 1, -1, &fe); CHKERRQ(ierr);
        ierr = PetscObjectSetName((PetscObject)fe, "temperature_"); CHKERRQ(ierr);
        ierr = DMAddField(dm, nullptr, (PetscObject)fe); CHKERRQ(ierr);
        fe_fields.push_back(fe);
    }

    // Create DS (required before adding boundaries in PETSc 3.22.2)
    ierr = DMCreateDS(dm); CHKERRQ(ierr);

    // Add boundary conditions to the DS
    // This must be done after DMCreateDS but before the section is finalized
    ierr = setupBoundaryConditions(); CHKERRQ(ierr);

    // Get the DS for later use
    ierr = DMGetDS(dm, &prob); CHKERRQ(ierr);

    PetscFunctionReturn(0);
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
    PetscScalar unified_constants[27] = {0};  // [0-24]=fluid/elasticity/poroelasticity, [25-26]=cohesive fault

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

        // Convert permeability from mD -> m^2
        constexpr double mD_to_m2 = 9.869233e-16;

        unified_constants[3]  = mat.porosity;
        unified_constants[4]  = mat.permeability_x * mD_to_m2;
        unified_constants[5]  = mat.permeability_y * mD_to_m2;
        unified_constants[6]  = mat.permeability_z * mD_to_m2;
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

    // Set constants ONCE for all physics (eliminates collision bug)
    // [0-24]=fluid/elasticity/poroelasticity, [25-26]=cohesive fault (if enabled)
    ierr = PetscDSSetConstants(prob, 27, unified_constants); CHKERRQ(ierr);

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
        ierr = PetscDSSetJacobian(prob, 0, 1,
                                  nullptr, PetscFEPoroelasticity::g1_pu,
                                  nullptr, nullptr); CHKERRQ(ierr);

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
        // Provide a basic Jacobian for pressure diffusion
        // g0 = phi*ct*shift, g3 = (k/mu)*I
        // For now, reuse the multiphase pressure Jacobian with default constants.
        ierr = PetscDSSetJacobian(prob, 0, 0, PetscFEFluidFlow::g0_BlackOilPressurePressure, nullptr, nullptr, PetscFEFluidFlow::g3_BlackOilPressurePressure); CHKERRQ(ierr);
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
            ierr = PetscDSSetResidual(prob, displacement_field_idx,
                                      PetscFEElasticity::f0_elastodynamics,
                                      PetscFEElasticity::f1_elastodynamics); CHKERRQ(ierr);
            ierr = PetscDSSetJacobian(prob, displacement_field_idx, displacement_field_idx,
                                      PetscFEElasticity::g0_elastodynamics, nullptr, nullptr,
                                      PetscFEElasticity::g3_elastodynamics); CHKERRQ(ierr);
        } else {
            // Quasi-static elasticity (no inertia)
            ierr = PetscDSSetResidual(prob, displacement_field_idx,
                                      PetscFEElasticity::f0_elastostatics,
                                      PetscFEElasticity::f1_elastostatics); CHKERRQ(ierr);
            ierr = PetscDSSetJacobian(prob, displacement_field_idx, displacement_field_idx,
                                      nullptr, nullptr, nullptr,
                                      PetscFEElasticity::g3_elastostatics); CHKERRQ(ierr);
        }

        use_fem_time_residual_ = true;
    }

    // Register cohesive fault callbacks if enabled
    if (config.enable_faults && cohesive_kernel_) {
        // Determine field indices
        // Field ordering: [fluid fields...] [displacement] [lagrange] [thermal] [...]
        PetscInt disp_field = 0;
        if (config.fluid_model == FluidModelType::SINGLE_COMPONENT) {
            disp_field = 1;
        } else if (config.fluid_model == FluidModelType::BLACK_OIL) {
            disp_field = 3;
        } else if (config.fluid_model == FluidModelType::COMPOSITIONAL) {
            disp_field = 4;  // P + 3 compositions
        }
        // else NONE: disp_field = 0

        // Lagrange multiplier field is next after displacement
        PetscInt lagrange_field = disp_field + 1;

        // Register cohesive callbacks via CohesiveFaultKernel::registerWithDS
        // This expands the constants array and registers PetscDSSetBdResidual/BdJacobian
        ierr = cohesive_kernel_->registerWithDS(prob, disp_field, lagrange_field);
        CHKERRQ(ierr);

        if (rank == 0) {
            PetscPrintf(comm, "Registered cohesive fault callbacks on fields %d (disp), %d (lagrange)\n",
                        disp_field, lagrange_field);
        }
    }

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
    // Use TSCN (Crank-Nicolson) for elastodynamics: second-order accurate, implicit
    // Use TSBEULER for other cases (single-phase flow, poroelasticity): first-order, more stable
    // Note: TSALPHA2 with I2Function/I2Jacobian is not available in PETSc 3.20 DMPlex FEM
    if (config.enable_elastodynamics) {
        ierr = TSSetType(ts, TSCN); CHKERRQ(ierr);
        if (rank == 0) {
            PetscPrintf(comm, "Using TSCN for elastodynamics (implicit, second-order)\n");
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

    // TODO: Read from config once FaultConfig parsing is implemented
    // For now, hardcode a vertical fault at x=Lx/2
    center[0] = grid_config.Lx / 2.0;
    center[1] = grid_config.Ly / 2.0;
    center[2] = grid_config.Lz / 2.0;
    length = 2.0 * std::max({grid_config.Lx, grid_config.Ly, grid_config.Lz});
    width = length;

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
    cohesive_kernel_->setMode(true);  // locked
    cohesive_kernel_->setFrictionCoefficient(0.6);

    if (rank == 0) {
        PetscPrintf(comm, "Fault network setup complete: %zu fault vertices\n",
                    cohesive_fault_->numVertices());
    }

    PetscFunctionReturn(0);
}

PetscErrorCode Simulator::setInitialConditions() {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

    // Start from zero (displacement at rest, zero pressure perturbation)
    ierr = VecZeroEntries(solution); CHKERRQ(ierr);

    // Ensure DM section is created (required for DMPlexInsertBoundaryValues to work)
    ierr = DMSetUp(dm); CHKERRQ(ierr);

    // Debug: Check if we have boundary IDs registered
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

        ierr = DMPlexTSComputeIJacobianFEM(dm, t, locU, locU_t, a, J, P, ctx); CHKERRQ(ierr);

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
    // PetscErrorCode ierr;  // Not used in current implementation
    
    // Disable VTK output for now to avoid format incompatibility
    // TODO: Properly implement output based on DM type (DMPlex)
    if (rank == 0 && step % 10 == 0) {
        PetscPrintf(comm, "Output at step %d would be written here\n", step);
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

    // For elastostatics/elastodynamics:
    // - Fixed bottom (z_min): all displacement components = 0
    // - Applied compression on top (z_max): u_z = -0.001
    if (displacement_field >= 0) {
        PetscInt field = displacement_field;
        PetscInt comps_all[3] = {0, 1, 2}; // All displacement components
        PetscInt label_value = 1; // We labeled boundary faces with value 1

        // Fixed bottom (z_min): all components = 0
        if (label_zmin) {
            ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "fixed_bottom", label_zmin, 1, &label_value,
                                 field, 3, comps_all, (void (*)(void))bc_zero, nullptr, nullptr, nullptr);
            CHKERRQ(ierr);
            if (rank == 0) {
                PetscPrintf(comm, "  Applied BC: fixed bottom (z_min), field %d, all components = 0\n", field);
            }
        }

        // Applied compression on top (z_max): all components specified
        if (label_zmax) {
            ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "compression_top", label_zmax, 1, &label_value,
                                 field, 3, comps_all, (void (*)(void))bc_compression, nullptr, nullptr, nullptr);
            CHKERRQ(ierr);
            if (rank == 0) {
                PetscPrintf(comm, "  Applied BC: compression top (z_max), field %d, u_z = -0.001\n", field);
            }
        }
    }

    // For poroelasticity: drained top (z_max, pressure = 0)
    if (pressure_field >= 0 && config.solid_model == SolidModelType::POROELASTIC) {
        PetscInt field = pressure_field;
        PetscInt comp = 0; // Pressure is scalar (1 component)
        PetscInt label_value = 1;

        if (label_zmax) {
            ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "drained_top", label_zmax, 1, &label_value,
                                 field, 1, &comp, (void (*)(void))bc_drained, nullptr, nullptr, nullptr);
            CHKERRQ(ierr);
            if (rank == 0) {
                PetscPrintf(comm, "  Applied BC: drained top (z_max), field %d, pressure = 0\n", field);
            }
        }
    }

    // For fluid flow only (no geomechanics):
    // - Fixed pressure on x_min: p = 0 (reference)
    if (pressure_field >= 0 && !config.enable_geomechanics &&
        config.solid_model != SolidModelType::POROELASTIC) {
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
