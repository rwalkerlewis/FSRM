#include "io/ConfigValidator.hpp"
#include "core/ConfigReader.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace FSRM {

namespace {

// ---------------------------------------------------------------------------
// Static schema. One canonical place for every section and key the
// FSRM parser actually consumes. New sections or keys must land
// here before they are accepted by strict validation.
//
// Sources:
//   * src/core/ConfigReader.cpp dispatch sites
//   * src/core/Simulator.cpp hand-rolled reader.getXxx("SECTION", ...)
//     calls (NEAR_FIELD_SOURCE, EXPLOSION_SOURCE, BOUNDARY_CONDITIONS,
//     ABSORBING_BC, SEISMOMETERS, FAULT, FRACTURE_PLANE,
//     HYDRAULIC_FRACTURE, INJECTION, INITIAL_CONDITIONS, MATERIAL,
//     MESH_REFINEMENT, NUCLEAR_TRIGGER, OUTPUT, PLASTICITY,
//     SOURCE_DISTRIBUTION, SEISMICITY, THERMAL, VISCOELASTIC,
//     WAVEFORM_VV, etc.)
//   * Per-instance keys for LAYER_<n>, SEISMOMETER_<n>,
//     MATERIAL_REGION_<n>, ROCK_<name>, FAULT_<name>.
// ---------------------------------------------------------------------------
using KeySet = std::unordered_set<std::string>;
using SectionMap = std::unordered_map<std::string, KeySet>;

const SectionMap& knownSections() {
    static const SectionMap kSchema = {
        {"META", {
            "strict_validation",
        }},
        {"SIMULATION", {
            "name", "start_time", "end_time", "dt_initial", "dt_min",
            "dt_max", "max_timesteps", "output_frequency",
            "adaptive_timestepping", "enable_checkpointing",
            "checkpoint_frequency", "rtol", "atol",
            "max_nonlinear_iterations", "max_linear_iterations",
            "enable_geomechanics", "enable_thermal", "enable_fractures",
            "enable_particle_transport", "enable_faults",
            "enable_tidal_forces", "enable_elastodynamics",
            "enable_poroelastodynamics", "enable_explosion_source",
            "enable_near_field_damage", "enable_hydrodynamic",
            "enable_crater_formation", "enable_atmospheric_blast",
            "enable_atmospheric_acoustic", "enable_infrasound",
            "enable_thermal_radiation", "enable_emp", "enable_fallout",
            "explosion_solve_mode", "fluid_model", "solid_model",
            "output_format", "enable_gravity", "gravity", "K0",
        }},
        {"MODULES", { "enabled" }},
        {"GRID", {
            "nx", "ny", "nz", "Lx", "Ly", "Lz",
            "origin_x", "origin_y", "origin_z",
            "mesh_type", "use_unstructured", "mesh_file",
            "gmsh_physical_volume", "gmsh_boundaries",
            "gmsh_refinement_level", "gmsh_material_mapping",
            "gmsh_fault_mapping",
            "input_crs", "model_crs", "use_local_coordinates",
            "local_origin_x", "local_origin_y", "local_origin_z",
            "auto_detect_utm", "min_cell_volume", "max_aspect_ratio",
            // Bounding-box style aliases used by some examples /
            // unit tests (xmin/xmax/ymin/ymax/zmin/zmax). The
            // simulator parser computes Lx/Ly/Lz from these.
            "xmin", "xmax", "ymin", "ymax", "zmin", "zmax",
        }},
        {"MATERIAL", {
            "assignment", "gravity", "heterogeneous",
            "velocity_model_file", "default_density",
            "default_youngs_modulus", "default_poisson_ratio",
            // Informational count hint frequently written by
            // historic-event fixtures and templates. The parser
            // discovers layers from getSectionsMatching("LAYER_").
            "num_layers",
        }},
        {"FLUID", {
            "type", "density", "viscosity", "compressibility",
            "oil_density_std", "gas_density_std", "water_density_std",
            "oil_viscosity", "gas_viscosity", "water_viscosity",
            "solution_GOR", "component_mw", "component_Tc",
            "component_Pc", "component_omega",
            "water_compressibility",
            "reference_pressure",
        }},
        {"ROCK", {
            "name", "density", "porosity", "youngs_modulus",
            "poisson_ratio", "poissons_ratio",
            "permeability_x", "permeability_y", "permeability_z",
            "thermal_conductivity", "thermal_expansion",
            "heat_capacity",
            "constitutive_model",
            // Aliases / fundamental moduli the test fixtures often
            // write directly. The parser computes mu / lambda from
            // youngs+poisson by default; these forms are accepted
            // for backward compatibility.
            "lambda", "mu", "shear_modulus",
            "biot_coefficient", "reference_pressure",
            "compressibility",
        }},
        {"PARTICLE", {
            "diameter", "density", "concentration", "diffusivity",
            "enable_settling", "enable_bridging",
        }},
        {"DYNAMICS", {
            "enable", "static_triggering", "trigger_threshold",
            "event_duration", "damping_alpha", "damping_beta",
            "quality_factor", "dynamic_permeability",
            "permeability_sensitivity", "permeability_recovery_time",
        }},
        {"SEISMICITY", {
            "enable", "friction_law", "nucleation_size",
            "seismic_slip_rate", "b_value", "aftershocks",
            "stress_transfer", "catalog_file",
        }},
        {"IMEX", {
            "enabled", "initial_mode", "implicit_dt_initial",
            "implicit_dt_min", "implicit_dt_max", "implicit_method",
            "explicit_dt_initial", "explicit_dt_min", "explicit_dt_max",
            "cfl_factor", "explicit_method", "trigger_type",
            "stress_threshold", "coulomb_threshold",
            "slip_rate_threshold", "velocity_threshold",
            "acceleration_threshold", "energy_rate_threshold",
            "settling_type", "min_dynamic_duration",
            "max_dynamic_duration", "settling_velocity",
            "settling_energy_ratio", "settling_slip_rate",
            "settling_observation_window", "adaptive_timestep",
            "dt_growth_factor", "dt_shrink_factor",
            "implicit_max_iterations", "explicit_max_iterations",
            "implicit_rtol", "explicit_rtol", "implicit_atol",
            "explicit_atol", "use_lumped_mass",
            "implicit_output_frequency", "explicit_output_frequency",
            "log_transitions", "transition_log_file",
            "smooth_transition", "transition_ramp_time",
        }},
        {"OUTPUT", {
            "format", "path", "frequency", "directory",
            "output_directory", "save_solution", "output_stress",
            "output_cfs", "cfs_friction", "cfs_receiver_strike",
            "cfs_receiver_dip",
            "pressure", "displacement", "stress", "strain",
            "velocity", "permeability", "temperature", "saturation",
            "fault_slip", "seismic_catalog",
            "pressure_unit", "displacement_unit", "stress_unit",
            "permeability_unit", "temperature_unit", "density_unit",
            "viscosity_unit", "length_unit", "time_unit",
            // Pass-13c wavefield + source-ball-3D output sub-keys.
            "wavefield_format", "wavefield_cadence_steps",
            "wavefield_fields", "wavefield_output_directory",
            "wavefield_basename",
            "source_ball_3d_output_format",
            "source_ball_3d_output_cadence_steps",
        }},
        {"BOUNDARY_CONDITIONS", {
            // Conventional per-face shorthand
            "top", "bottom", "sides",
            // Per-face explicit forms used by the simulator parser
            "x_min", "x_max", "y_min", "y_max", "z_min", "z_max",
            "x_min_traction", "x_max_traction",
            "y_min_traction", "y_max_traction",
            "z_min_traction", "z_max_traction",
            "x_min_pressure", "x_max_pressure",
            "y_min_pressure", "y_max_pressure",
            "z_min_pressure", "z_max_pressure",
            "x_min_displacement", "x_max_displacement",
            "y_min_displacement", "y_max_displacement",
            "z_min_displacement", "z_max_displacement",
            // Component-specific traction values used by tests.
            "top_traction_x", "top_traction_y", "top_traction_z",
        }},
        {"INITIAL_CONDITIONS", {
            "type", "path", "field", "distribution", "value",
            "gradient", "file",
        }},
        {"ABSORBING_BC", {
            "enabled", "x_min", "x_max", "y_min", "y_max",
            "z_min", "z_max",
        }},
        {"PLASTICITY", {
            "enabled", "cohesion", "friction_angle", "dilation_angle",
            "hardening_modulus",
        }},
        {"VISCOELASTIC", {
            "enabled", "num_mechanisms", "f_min", "f_max",
            "q_p", "q_s",
        }},
        {"THERMAL", {
            "thermal_conductivity", "specific_heat",
            "thermal_expansion", "reference_temperature",
            "bottom_temperature", "top_temperature",
        }},
        {"INJECTION", {
            "enabled", "x", "y", "z", "rate", "start_time", "end_time",
        }},
        {"FAULT", {
            "mode", "center_x", "center_y", "center_z",
            "strike", "dip", "length", "width",
            "friction_coefficient", "friction_model",
            "static_friction", "dynamic_friction",
            "critical_slip_distance", "cohesion",
            "initial_normal_stress", "initial_shear_stress",
            "initial_normal_traction", "initial_shear_traction",
            "nucleation_center_x", "nucleation_center_z",
            "nucleation_radius", "nucleation_shear_stress",
            "slip_strike", "slip_dip", "slip_opening",
            "slip_onset_time", "slip_rise_time",
        }},
        {"FRACTURE_PLANE", {
            "enabled", "center_x", "center_y", "center_z",
            "strike", "dip", "length", "width",
        }},
        {"HYDRAULIC_FRACTURE", {
            "model", "enable_fem_pressurized", "uniform_pressure_pa",
            "fluid_density", "fluid_viscosity", "youngs_modulus",
            "poissons_ratio", "toughness", "height",
            "max_horizontal_stress", "min_horizontal_stress",
            "vertical_stress",
        }},
        {"MESH_REFINEMENT", {
            "enabled", "refinement_levels", "source_radius_factor",
        }},
        {"NUCLEAR_TRIGGER", {
            "enabled", "time0", "rise_time", "decay_time",
            "delta_pore_pressure", "delta_sigma_n", "delta_tau",
        }},
        {"EXPLOSION_SOURCE", {
            "type", "yield_kt", "depth_of_burial",
            "location_x", "location_y", "location_z",
            "source_x", "source_y", "source_z",
            "onset_time", "time0", "detonation_time",
            "rise_time", "cavity_overpressure", "medium_type",
            "apply_damage_zone", "explosion_solve_mode", "mode",
        }},
        {"NEAR_FIELD_SOURCE", {
            "mode", "damage_model", "elastic_radius_factor",
            "near_field_dt", "output_cadence_microseconds",
            "profile_output_cadence_microseconds",
            // Pass-6 / pass-7
            "solver_kind", "radial_cells", "radial_outer_factor",
            "cfl", "art_visc_linear", "art_visc_quadratic",
            "gas_eos_gamma", "cavity_eos", "tillotson_parameter_set",
            "cavity_initialization", "initial_cavity_radius_m",
            "radiation_transition_time_s",
            // Pass-8 / pass-9
            "radiation_phase", "opacity_model", "operator_splitting",
            "tabulated_eos_table_path",
            "tabulated_opacity_rosseland_path",
            "tabulated_opacity_planck_path",
            "tabulated_eos_blend_lower_pa",
            "tabulated_eos_blend_upper_pa",
            "tabulated_opacity_blend_lower_k",
            "tabulated_opacity_blend_upper_k",
            "kappa_constant_m2_per_kg",
            "tillotson_extrapolation_warning_threshold_pa",
            "radial_outer_radius_m",
            // Pass-10
            "time_integrator", "radiation_n_groups",
            "radiation_freq_min_hz", "radiation_freq_max_hz",
            "radiation_simpson_points", "output_per_group_radiation",
            "operator_splitting_convergence_diagnostic",
            "radiation_max_newton_iter", "radiation_newton_tolerance",
            "radiation_handoff_debounce_steps",
            // Pass-11 / pass-13b
            "time_integrator_diffusion", "sponge_layer_enabled",
            "sponge_layer_thickness_factor",
            "cavity_geometry", "cavity_radius_m", "outer_radius_m",
            "mesh_path", "overburden_K0",
            "source_ball_radiation_discretization",
            // Pass-13c
            "source_ball_radiation_substep_cadence",
            // Pass-14a source forcing
            "source_forcing_enabled", "source_time_function",
            "source_yield_kt", "source_deposition_duration_s",
            "source_deposition_efficiency",
        }},
        {"SOURCE_DISTRIBUTION", {
            "mode", "support_radius_factor", "gaussian_sigma_factor",
            "min_cells",
        }},
        {"SEISMOMETERS", {
            "enabled", "output_dir", "start_time_utc",
            "formats", "default_sample_rate_hz", "default_quantity",
        }},
        {"WAVEFORM_VV", {
            "enabled", "cache_root", "events",
            "metric_freq_band_hz_min", "metric_freq_band_hz_max",
            "metric_lag_window_s",
        }},
        {"TRACTION_BC", {
            "enabled",
            "top_traction_x", "top_traction_y", "top_traction_z",
            "bottom_traction_x", "bottom_traction_y", "bottom_traction_z",
        }},
        // Optional sections kept for forward compatibility with
        // historical templates; the parser tolerates them but they
        // do not feed any active code path.
        {"GEOMECHANICS", { "*" }},
        {"MESH",         { "*" }},
        {"BOUNDARY",     { "*" }},
    };
    return kSchema;
}

// Dynamic-prefix sections. The validator accepts any section
// whose name starts with one of these prefixes followed by an
// integer or a name token; the keys inside that section are
// validated against the prefix's allowed-key set.
struct DynamicSection {
    std::string prefix;            // e.g. "LAYER_", "SEISMOMETER_"
    bool numeric_suffix_only;      // require [PREFIX][0-9]+
    KeySet keys;
};

const std::vector<DynamicSection>& dynamicSections() {
    static const std::vector<DynamicSection> kPrefixes = {
        {"LAYER_", true, {
            "name", "z_top", "z_bottom", "z_min", "z_max",
            "lambda", "mu", "rho", "density", "vp", "vs",
            "shear_modulus",
            // Per-layer attenuation overrides used by the layered-Q
            // integration tests.
            "q_p", "q_s",
        }},
        {"SEISMOMETER_", true, {
            "sta", "location_xyz", "loc",
            "quantity", "sample_rate_hz", "format",
        }},
        {"MATERIAL_REGION_", true, {
            "gmsh_label", "density", "youngs_modulus",
            "poissons_ratio", "poisson_ratio",
            "lambda", "mu", "rho", "name",
        }},
        // ROCK_<name> sections (parsed via getSectionsMatching("ROCK"))
        // share the same keys as the primary [ROCK] block.
        {"ROCK_", false, {
            "name", "density", "porosity", "youngs_modulus",
            "poisson_ratio", "poissons_ratio",
            "permeability_x", "permeability_y", "permeability_z",
            "thermal_conductivity", "thermal_expansion",
            "heat_capacity",
            "constitutive_model",
        }},
        // FAULT_<name> sections (parsed via getSectionsMatching("FAULT"))
        // accept the per-fault extended schema documented in
        // ConfigReader::parseFaultsExtended.
        {"FAULT_", false, {
            "name", "strike", "dip", "length", "width",
            "static_friction", "dynamic_friction", "cohesion",
            "use_rate_state", "a_parameter", "b_parameter",
            "Dc_parameter", "fault_mode", "slip_strike", "slip_dip",
            "slip_opening", "use_split_nodes", "split_node_method",
            "traction_type", "penalty_normal", "penalty_tangent",
            "prescribed_traction_normal", "prescribed_traction_strike",
            "prescribed_traction_dip", "cohesive_strength",
            "critical_opening", "critical_slip", "allow_separation",
            // Numeric-suffix variants like FAULT_1
            "mode", "center_x", "center_y", "center_z",
            "friction_coefficient", "friction_model",
            "critical_slip_distance",
            "initial_normal_stress", "initial_shear_stress",
            "nucleation_center_x", "nucleation_center_z",
            "nucleation_radius", "nucleation_shear_stress",
            "slip_onset_time", "slip_rise_time",
        }},
    };
    return kPrefixes;
}

// Sections under archived / aspirational templates that the
// validator silently ignores. These never reach a runnable example;
// keeping them here lets `config/complete_template.config` parse
// without error while strict validation still catches typos in the
// example tree.
const std::unordered_set<std::string>& archivedTemplateSections() {
    static const std::unordered_set<std::string> kArchived = {
        "WELLS", "FRACTURES", "FAULTS", "PARTICLES",
        "BOUNDARY_CONDITIONS_LIST", "INITIAL_CONDITIONS_LIST",
    };
    return kArchived;
}

// AMR sub-section family ("amr", "amr.criterion", ...). Permissive
// because the AMR config parser uses dot-paths and is a separate
// schema universe.
bool isAMRSection(const std::string& s) {
    return s == "amr" || (s.size() >= 4 && s.substr(0, 4) == "amr.");
}

// Deprecated key forms with custom error messages. Each entry
// matches a (section-pattern, key-pattern) pair from the pre-PR
// #128 historic-event configs.
struct DeprecatedForm {
    std::string section_pattern;   // exact match or empty (any)
    std::string key_pattern;       // regex-lite: "*" matches anything
    std::string message;
};

// Match a key with a glob-lite pattern. Supports "*" suffix only:
// "station_*" matches "station_1", "station_42", etc.
bool matchKey(const std::string& pattern, const std::string& key) {
    if (pattern.empty()) return true;
    auto star = pattern.find('*');
    if (star == std::string::npos) return pattern == key;
    if (key.size() < star) return false;
    return key.compare(0, star, pattern, 0, star) == 0;
}

const std::vector<DeprecatedForm>& deprecatedForms() {
    static const std::vector<DeprecatedForm> kDeprecated = {
        // Pre-PR #128 inline-station block (config grammar drift)
        {"SEISMOMETERS", "station_*",
         "Inline station_<n> keys in [SEISMOMETERS] are no longer "
         "consumed. Use a [SEISMOMETER_<n>] block per station with "
         "keys `sta = <name>` and `location_xyz = x,y,z`. See "
         "docs/CONFIGURATION_VALIDATION.md \"deprecated forms\"."},
        {"SEISMOMETERS", "station_count",
         "station_count is no longer consumed. The seismometer count "
         "is inferred from the number of [SEISMOMETER_<n>] blocks. "
         "Remove this key."},
        {"SEISMOMETERS", "sac_sampling_rate_hz",
         "sac_sampling_rate_hz was renamed to default_sample_rate_hz "
         "after PR #128. Update the key name."},
        {"SEISMOMETERS", "hdf5_enabled",
         "hdf5_enabled in [SEISMOMETERS] was retired in PR #128. "
         "Set [SEISMOMETERS] formats = HDF5 (or comma-separated "
         "list) instead."},
        {"SEISMOMETERS", "sac_enabled",
         "sac_enabled in [SEISMOMETERS] was retired in PR #128. "
         "Set [SEISMOMETERS] formats = SAC instead."},
        // PR #127 shipped some configs with an [EXPLOSION] block
        // (no _SOURCE) that the parser silently ignores.
        {"EXPLOSION", "*",
         "[EXPLOSION] is not a recognised section. Rename to "
         "[EXPLOSION_SOURCE]. Required keys: type, yield_kt, "
         "depth_of_burial, location_x, location_y, location_z, "
         "onset_time, rise_time, cavity_overpressure."},
    };
    return kDeprecated;
}

// Sections that, when present, force [BOUNDARY_CONDITIONS] to
// also be present. Catches the PR #127 bug: an explosion-source
// example with no BCs falls back to "compression top, fixed
// bottom" Dirichlet defaults that suppress the explosion coupling.
const std::vector<std::string>& sectionsRequiringBCs() {
    static const std::vector<std::string> kSections = {
        "EXPLOSION_SOURCE",
    };
    return kSections;
}

// True if `s` matches one of the configured dynamic prefixes.
const DynamicSection* findDynamicSection(const std::string& s) {
    for (const auto& d : dynamicSections()) {
        if (s.size() <= d.prefix.size()) continue;
        if (s.compare(0, d.prefix.size(), d.prefix) != 0) continue;
        const std::string suffix = s.substr(d.prefix.size());
        if (d.numeric_suffix_only) {
            bool all_digits = !suffix.empty() &&
                std::all_of(suffix.begin(), suffix.end(),
                            [](char c) { return std::isdigit(static_cast<unsigned char>(c)); });
            if (!all_digits) continue;
        }
        return &d;
    }
    return nullptr;
}

std::string formatLocation(const std::string& filename,
                           const std::string& section,
                           const std::string& key) {
    std::ostringstream os;
    if (!filename.empty()) os << filename << ": ";
    os << "[" << section << "]";
    if (!key.empty()) os << "." << key;
    return os.str();
}

// Build a small "did you mean" hint for the closest known key
// inside a given section. Naive Levenshtein over the small
// per-section key vocabulary; cheap and focused.
int editDistance(const std::string& a, const std::string& b) {
    const size_t m = a.size(), n = b.size();
    std::vector<std::vector<int>> d(m + 1, std::vector<int>(n + 1, 0));
    for (size_t i = 0; i <= m; ++i) d[i][0] = static_cast<int>(i);
    for (size_t j = 0; j <= n; ++j) d[0][j] = static_cast<int>(j);
    for (size_t i = 1; i <= m; ++i) {
        for (size_t j = 1; j <= n; ++j) {
            int cost = (a[i - 1] == b[j - 1]) ? 0 : 1;
            d[i][j] = std::min({d[i - 1][j] + 1, d[i][j - 1] + 1,
                                d[i - 1][j - 1] + cost});
        }
    }
    return d[m][n];
}

std::string suggestKey(const std::string& key, const KeySet& known) {
    if (known.empty()) return "";
    int best = std::numeric_limits<int>::max();
    std::string best_match;
    for (const auto& k : known) {
        const int d = editDistance(key, k);
        if (d < best) { best = d; best_match = k; }
    }
    if (best <= 3 && best_match != key) return best_match;
    return "";
}

} // namespace

ConfigValidator::Result ConfigValidator::validate(
    const ConfigReader& reader, const std::string& source_filename) {
    Result out;

    // Env-var emergency bypass. Used by the examples-runtime CTest
    // gate to selectively re-enable strict validation, and by
    // legacy callers that need to skip validation outright. Set
    // FSRM_DISABLE_STRICT_VALIDATION=1 to skip validation entirely.
    if (const char* env =
            std::getenv("FSRM_DISABLE_STRICT_VALIDATION")) {
        if (env[0] == '1' || env[0] == 't' || env[0] == 'T' ||
            env[0] == 'y' || env[0] == 'Y') {
            out.valid = true;
            return out;
        }
    }

    // Opt-out flag. Honour it but always warn.
    const bool strict_opt_out =
        reader.hasSection("META") &&
        reader.getBool("META", "strict_validation", true) == false;
    if (strict_opt_out) {
        std::ostringstream w;
        if (!source_filename.empty()) w << source_filename << ": ";
        w << "[META] strict_validation = false. Strict config "
             "validation has been disabled for this file. This is "
             "intended only for short-lived debugging; remove the "
             "opt-out before committing.";
        out.warnings.push_back(w.str());
        out.valid = true;
        return out;
    }

    const auto& schema = knownSections();
    const auto& archived = archivedTemplateSections();
    const auto& deprecated = deprecatedForms();

    // 1. Walk every section in the parsed config.
    const auto sections = reader.getSections();
    for (const auto& section : sections) {
        // Permissive paths first.
        if (isAMRSection(section)) continue;
        if (archived.count(section)) continue;

        // Known fixed section?
        auto schema_it = schema.find(section);
        const KeySet* allowed = nullptr;
        const DynamicSection* dyn = nullptr;
        if (schema_it != schema.end()) {
            allowed = &schema_it->second;
        } else {
            dyn = findDynamicSection(section);
            if (dyn != nullptr) allowed = &dyn->keys;
        }

        if (!allowed) {
            std::ostringstream e;
            if (!source_filename.empty()) e << source_filename << ": ";
            e << "Unknown section [" << section << "]. ";
            e << "If this is a new section, add it to the registry "
                 "in src/io/ConfigValidator.cpp (knownSections "
                 "or dynamicSections). If this is a typo, the "
                 "closest known sections include [SIMULATION], "
                 "[GRID], [EXPLOSION_SOURCE], [SEISMOMETERS], "
                 "[BOUNDARY_CONDITIONS].";
            out.errors.push_back(e.str());
            out.valid = false;
            continue;
        }

        // 2. Walk every key inside the section.
        const auto keys = reader.getKeys(section);
        // Wildcard-permissive sections (e.g. archived templates).
        const bool wildcard = allowed->count("*") > 0;
        for (const auto& key : keys) {
            // 2a. Deprecated form?
            bool was_deprecated = false;
            for (const auto& dep : deprecated) {
                if (dep.section_pattern != section &&
                    dep.section_pattern != "") {
                    continue;
                }
                if (!matchKey(dep.key_pattern, key)) continue;
                std::ostringstream e;
                e << "Deprecated form at "
                  << formatLocation(source_filename, section, key)
                  << ": " << dep.message;
                out.errors.push_back(e.str());
                out.valid = false;
                was_deprecated = true;
                break;
            }
            if (was_deprecated) continue;
            if (wildcard) continue;

            // 2b. Unknown key?
            if (!allowed->count(key)) {
                std::ostringstream e;
                e << "Unknown key at "
                  << formatLocation(source_filename, section, key);
                const std::string hint = suggestKey(key, *allowed);
                if (!hint.empty()) {
                    e << ". Did you mean `" << hint << "`?";
                } else {
                    e << ". This key is not in the schema for ["
                      << section << "]. If it is a new key, add it "
                      << "to src/io/ConfigValidator.cpp.";
                }
                out.errors.push_back(e.str());
                out.valid = false;
            }
        }
    }

    // 3. Cross-section requirement: explosion-source configs must
    //    declare [BOUNDARY_CONDITIONS] explicitly so the default
    //    "compression top, fixed bottom" Dirichlet does not
    //    silently suppress the explosion coupling.
    for (const auto& trigger : sectionsRequiringBCs()) {
        if (!reader.hasSection(trigger)) continue;
        if (reader.hasSection("BOUNDARY_CONDITIONS")) continue;
        std::ostringstream e;
        if (!source_filename.empty()) e << source_filename << ": ";
        e << "[" << trigger << "] is present but "
             "[BOUNDARY_CONDITIONS] is missing. Without an explicit "
             "BC block the simulator falls back to a Dirichlet "
             "compression-top / fixed-bottom default that suppresses "
             "the explosion source coupling. Add a [BOUNDARY_CONDITIONS] "
             "block (e.g. bottom = free, sides = free, top = free for "
             "a standard underground-test config). See "
             "docs/CONFIGURATION_VALIDATION.md.";
        out.errors.push_back(e.str());
        out.valid = false;
    }

    return out;
}

bool ConfigValidator::validateAndReport(const ConfigReader& reader,
                                        const std::string& source_filename) {
    const Result r = validate(reader, source_filename);
    for (const auto& w : r.warnings) {
        std::cerr << "[ConfigValidator][warn] " << w << std::endl;
    }
    for (const auto& e : r.errors) {
        std::cerr << "[ConfigValidator][error] " << e << std::endl;
    }
    return r.valid;
}

} // namespace FSRM
