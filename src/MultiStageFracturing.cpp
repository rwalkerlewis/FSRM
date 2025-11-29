#include "MultiStageFracturing.hpp"
#include <algorithm>
#include <numeric>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace FSRM {

// ============================================================================
// PerforationCluster Implementation
// ============================================================================

PerforationCluster::PerforationCluster()
    : cluster_id(0), measured_depth(0.0), true_vertical_depth(0.0),
      num_perforations(6), perforation_diameter(0.01), perforation_length(0.3),
      phasing(60.0), shot_density(13.0),
      current_diameter(0.01), discharge_coefficient(0.85), friction_multiplier(1.0),
      is_open(true), total_flow_area(0.0),
      min_horizontal_stress(30e6), max_horizontal_stress(35e6), vertical_stress(45e6),
      pore_pressure(25e6), stress_shadow_delta(0.0),
      allocated_rate(0.0), cumulative_volume(0.0), proppant_placed(0.0)
{
    total_flow_area = num_perforations * M_PI * current_diameter * current_diameter / 4.0;
}

double PerforationCluster::calculatePerfFriction(double rate, double fluid_density) const {
    if (!is_open || num_perforations <= 0 || current_diameter <= 0) {
        return 1e20; // Infinite friction if closed
    }
    
    // Orifice equation: dP = rho * Q^2 / (2 * Cd^2 * A^2)
    double area = num_perforations * M_PI * current_diameter * current_diameter / 4.0;
    double dp = fluid_density * rate * rate / 
                (2.0 * discharge_coefficient * discharge_coefficient * area * area);
    
    return dp * friction_multiplier;
}

void PerforationCluster::updateErosion(double rate, double proppant_conc, double dt) {
    if (rate <= 0 || !is_open) return;
    
    // Calculate velocity through perforations
    double area = num_perforations * M_PI * current_diameter * current_diameter / 4.0;
    double velocity = rate / area;
    
    // Erosion rate proportional to velocity^2 * proppant concentration
    // Based on API RP 19B erosion model
    double erosion_coeff = 1e-11; // Calibration coefficient
    double hardness_factor = 1.0; // Sand = 1.0, ceramic = 0.5
    
    double erosion_rate = erosion_coeff * velocity * velocity * proppant_conc * hardness_factor;
    
    // Update diameter
    current_diameter += erosion_rate * dt;
    
    // Update flow area
    total_flow_area = num_perforations * M_PI * current_diameter * current_diameter / 4.0;
    
    // Update discharge coefficient (increases with erosion)
    discharge_coefficient = std::min(0.95, 0.7 + 0.25 * (current_diameter / perforation_diameter - 1.0));
}

double PerforationCluster::calculateBreakdownPressure() const {
    // Hubbert-Willis breakdown pressure for vertical fracture:
    // Pb = 3*Shmin - SHmax + T - Pp
    // where T is tensile strength (~5-10 MPa for most rocks)
    double tensile_strength = 5e6; // Pa
    
    double breakdown = 3.0 * min_horizontal_stress - max_horizontal_stress 
                     + tensile_strength - pore_pressure + stress_shadow_delta;
    
    return breakdown;
}

// ============================================================================
// TreatmentStage Implementation
// ============================================================================

TreatmentStage::TreatmentStage()
    : stage_id(0), stage_name("Stage 1"), start_md(0.0), end_md(0.0),
      isolation_method(IsolationMethod::PLUG), isolation_depth(0.0),
      isolation_confirmed(false), state(StageState::PENDING),
      start_time(0.0), end_time(0.0), total_volume(0.0), total_proppant(0.0),
      max_rate(0.0), avg_pressure(0.0), max_pressure(0.0),
      isip(0.0), closure_pressure(0.0),
      fracture_half_length(0.0), fracture_height(0.0),
      propped_half_length(0.0), propped_height(0.0),
      average_conductivity(0.0), cluster_efficiency(0.0), srw_volume(0.0)
{}

// ============================================================================
// PumpScheduleStep Implementation
// ============================================================================

PumpScheduleStep::PumpScheduleStep()
    : fluid_name("Slickwater"), rate(0.05), volume(100.0), duration(0.0),
      proppant_concentration(0.0), proppant_type("100 mesh"),
      is_pad(false), is_flush(false),
      fluid_viscosity(0.001), fluid_density(1000.0)
{
    if (rate > 0) {
        duration = volume / rate;
    }
}

// ============================================================================
// FracFluid Implementation
// ============================================================================

FracFluid::FracFluid()
    : name("Slickwater"), type(FluidType::SLICKWATER),
      base_viscosity(0.001), density(1000.0),
      n_prime(1.0), k_prime(0.001),
      viscosity_temp_coeff(-0.02), breakdown_temp(393.15), breakdown_time(3600.0),
      leakoff_coefficient(5e-5), spurt_loss(0.0), wall_building_coeff(0.0),
      friction_reduction(0.7)
{}

double FracFluid::getViscosity(double temperature, double shear_rate, double time) const {
    // Power-law viscosity with temperature dependence
    double mu_base = k_prime * std::pow(shear_rate, n_prime - 1.0);
    
    // Temperature correction (Arrhenius-type)
    double T_ref = 298.15; // 25°C
    double temp_factor = std::exp(viscosity_temp_coeff * (temperature - T_ref));
    
    // Time-dependent breakdown for crosslinked gels
    double breakdown_factor = 1.0;
    if (type == FluidType::CROSSLINKED_GEL && temperature > breakdown_temp) {
        breakdown_factor = std::exp(-time / breakdown_time);
    }
    
    return mu_base * temp_factor * breakdown_factor;
}

double FracFluid::getLeakoff(double time, double pressure_diff) const {
    // Carter leakoff model: vL = CL / sqrt(t)
    // Total leakoff = CL * sqrt(t) + Sp
    if (time <= 0) return spurt_loss;
    
    double leakoff = leakoff_coefficient * std::sqrt(time) + spurt_loss;
    
    // Pressure-dependent correction
    double pressure_factor = std::sqrt(pressure_diff / 1e6); // Normalized to 1 MPa
    
    return leakoff * pressure_factor;
}

// ============================================================================
// Proppant Implementation
// ============================================================================

Proppant::Proppant()
    : name("100 Mesh Sand"), type(ProppantType::SAND),
      diameter(0.00015), specific_gravity(2.65), bulk_density(1600.0),
      absolute_density(2650.0), porosity(0.35), sphericity(0.7), roundness(0.6),
      crush_strength(35e6), fines_generated(0.05),
      beta_factor(1e8)
{}

double Proppant::getSettlingVelocity(double fluid_viscosity, double fluid_density) const {
    // Stokes settling for spherical particles
    double g = 9.81;
    double delta_rho = absolute_density - fluid_density;
    
    // Stokes velocity
    double v_stokes = delta_rho * g * diameter * diameter / (18.0 * fluid_viscosity);
    
    // Sphericity correction
    double shape_factor = 1.0 / (1.0 + 0.15 * std::pow(1.0 - sphericity, 0.5));
    
    // Hindered settling correction for concentrated suspensions
    // Will be applied externally based on local concentration
    
    return v_stokes * shape_factor;
}

double Proppant::getPermeability(double closure_stress) const {
    // Interpolate from stress-permeability table
    if (stress_points.empty() || permeability.empty()) {
        // Use correlation if no table
        double k0 = 1e-10; // Base permeability (m²)
        double stress_factor = std::exp(-2e-8 * closure_stress);
        return k0 * stress_factor;
    }
    
    // Linear interpolation
    for (size_t i = 0; i < stress_points.size() - 1; ++i) {
        if (closure_stress >= stress_points[i] && closure_stress <= stress_points[i+1]) {
            double t = (closure_stress - stress_points[i]) / 
                      (stress_points[i+1] - stress_points[i]);
            return permeability[i] + t * (permeability[i+1] - permeability[i]);
        }
    }
    
    return permeability.back();
}

double Proppant::getConductivity(double closure_stress, double width) const {
    double k = getPermeability(closure_stress);
    return k * width; // m² * m = m³, convert to mD·ft in output
}

double Proppant::getPackPorosity() const {
    // Correlation for random packing
    return 0.259 + 0.178 * sphericity;
}

// ============================================================================
// LimitedEntryDesign Implementation
// ============================================================================

LimitedEntryDesign::LimitedEntryDesign()
    : min_perf_friction_(3.5e6), max_perf_friction_(10e6), discharge_coeff_(0.85)
{}

std::vector<int> LimitedEntryDesign::designPerfCount(
    const std::vector<PerforationCluster>& clusters,
    double target_rate,
    double fluid_density,
    double target_perf_friction)
{
    std::vector<int> perf_counts(clusters.size());
    
    // Calculate required flow area for target friction
    // dP = rho * Q^2 / (2 * Cd^2 * A^2)
    // A = sqrt(rho * Q^2 / (2 * Cd^2 * dP))
    
    double total_area_required = std::sqrt(fluid_density * target_rate * target_rate /
                                           (2.0 * discharge_coeff_ * discharge_coeff_ * target_perf_friction));
    
    // Distribute among clusters based on stress (more perfs where stress is higher)
    std::vector<double> stress_factors(clusters.size());
    double max_stress = 0.0;
    for (const auto& c : clusters) {
        max_stress = std::max(max_stress, c.min_horizontal_stress + c.stress_shadow_delta);
    }
    
    double total_factor = 0.0;
    for (size_t i = 0; i < clusters.size(); ++i) {
        // More perfs at higher stress locations
        stress_factors[i] = (clusters[i].min_horizontal_stress + clusters[i].stress_shadow_delta) / max_stress;
        total_factor += stress_factors[i];
    }
    
    // Calculate perfs per cluster
    for (size_t i = 0; i < clusters.size(); ++i) {
        double area_per_cluster = total_area_required * stress_factors[i] / total_factor;
        double perf_area = M_PI * clusters[i].perforation_diameter * clusters[i].perforation_diameter / 4.0;
        perf_counts[i] = std::max(1, static_cast<int>(std::ceil(area_per_cluster / perf_area)));
    }
    
    return perf_counts;
}

std::pair<double, std::vector<double>> LimitedEntryDesign::predictClusterEfficiency(
    const std::vector<PerforationCluster>& clusters,
    double rate,
    const FracFluid& fluid_props)
{
    std::vector<double> flow_dist = calculateFlowDistribution(clusters, rate, fluid_props.density);
    
    // Calculate efficiency (fraction taking > 10% of fair share)
    double fair_share = rate / clusters.size();
    int active_clusters = 0;
    
    for (double flow : flow_dist) {
        if (flow > 0.1 * fair_share) {
            ++active_clusters;
        }
    }
    
    double efficiency = static_cast<double>(active_clusters) / clusters.size();
    
    return {efficiency, flow_dist};
}

std::vector<double> LimitedEntryDesign::calculateFlowDistribution(
    const std::vector<PerforationCluster>& clusters,
    double total_rate,
    double fluid_density)
{
    // Iterative flow distribution based on pressure equilibrium
    // All clusters see same treating pressure, but different friction + closure
    
    int n = clusters.size();
    std::vector<double> rates(n, total_rate / n); // Initial guess
    
    const int max_iter = 100;
    const double tol = 1e-6;
    
    for (int iter = 0; iter < max_iter; ++iter) {
        // Calculate pressure drop for each cluster at current rate
        std::vector<double> dp(n);
        for (int i = 0; i < n; ++i) {
            dp[i] = clusters[i].calculatePerfFriction(rates[i], fluid_density);
            dp[i] += clusters[i].min_horizontal_stress + clusters[i].stress_shadow_delta - clusters[i].pore_pressure;
        }
        
        // Find average pressure
        double avg_dp = std::accumulate(dp.begin(), dp.end(), 0.0) / n;
        
        // Adjust rates to equalize pressure
        double sum_rates = 0.0;
        for (int i = 0; i < n; ++i) {
            // Rate proportional to sqrt(target_dp - closure)
            double effective_dp = std::max(0.0, avg_dp - (clusters[i].min_horizontal_stress + 
                                            clusters[i].stress_shadow_delta - clusters[i].pore_pressure));
            rates[i] = std::sqrt(2.0 * effective_dp / fluid_density) * 
                      clusters[i].total_flow_area * clusters[i].discharge_coefficient;
            sum_rates += rates[i];
        }
        
        // Normalize to total rate
        for (int i = 0; i < n; ++i) {
            rates[i] *= total_rate / sum_rates;
        }
        
        // Check convergence
        double max_change = 0.0;
        for (int i = 0; i < n; ++i) {
            max_change = std::max(max_change, std::abs(rates[i] - total_rate / n) / (total_rate / n));
        }
        if (max_change < tol) break;
    }
    
    return rates;
}

// ============================================================================
// MultiStageFracManager Implementation
// ============================================================================

MultiStageFracManager::MultiStageFracManager()
    : casing_id_(0.1), tubing_id_(0.073), roughness_(1e-5)
{}

void MultiStageFracManager::setWellGeometry(const std::vector<double>& md,
                                             const std::vector<double>& tvd,
                                             const std::vector<double>& inclination,
                                             const std::vector<double>& azimuth) {
    well_md_ = md;
    well_tvd_ = tvd;
    well_inc_ = inclination;
    well_azi_ = azimuth;
}

void MultiStageFracManager::setWellboreProperties(double casing_id, double tubing_id,
                                                   double roughness) {
    casing_id_ = casing_id;
    tubing_id_ = tubing_id;
    roughness_ = roughness;
}

void MultiStageFracManager::setStressProfile(const std::vector<double>& depths,
                                              const std::vector<double>& shmin,
                                              const std::vector<double>& shmax,
                                              const std::vector<double>& sv,
                                              const std::vector<double>& pp) {
    stress_depths_ = depths;
    shmin_ = shmin;
    shmax_ = shmax;
    sv_ = sv;
    pp_ = pp;
}

void MultiStageFracManager::setRockProfile(const std::vector<double>& depths,
                                            const std::vector<double>& youngs_modulus,
                                            const std::vector<double>& poisson_ratio,
                                            const std::vector<double>& toughness) {
    rock_depths_ = depths;
    youngs_modulus_ = youngs_modulus;
    poisson_ratio_ = poisson_ratio;
    toughness_ = toughness;
}

void MultiStageFracManager::addCluster(const PerforationCluster& cluster) {
    clusters_.push_back(cluster);
    
    // Interpolate stress at cluster depth
    clusters_.back().min_horizontal_stress = interpolateStress(cluster.measured_depth, 
                                                                stress_depths_, shmin_);
    clusters_.back().max_horizontal_stress = interpolateStress(cluster.measured_depth,
                                                                stress_depths_, shmax_);
    clusters_.back().vertical_stress = interpolateStress(cluster.measured_depth,
                                                          stress_depths_, sv_);
    clusters_.back().pore_pressure = interpolateStress(cluster.measured_depth,
                                                        stress_depths_, pp_);
}

void MultiStageFracManager::addStage(const TreatmentStage& stage) {
    stages_.push_back(stage);
}

void MultiStageFracManager::addFluidSystem(const FracFluid& fluid) {
    fluids_[fluid.name] = fluid;
}

void MultiStageFracManager::addProppant(const Proppant& proppant) {
    proppants_[proppant.name] = proppant;
}

void MultiStageFracManager::setPumpSchedule(int stage_id, 
                                             const std::vector<PumpScheduleStep>& schedule) {
    pump_schedules_[stage_id] = schedule;
}

MultiStageFracManager::StageSimResult 
MultiStageFracManager::simulateStage(int stage_id, double dt) {
    StageSimResult result;
    
    // Find stage
    TreatmentStage* stage = getStage(stage_id);
    if (!stage) {
        throw std::runtime_error("Stage not found: " + std::to_string(stage_id));
    }
    
    // Get clusters for this stage
    std::vector<PerforationCluster*> stage_clusters;
    for (int cid : stage->cluster_ids) {
        PerforationCluster* cluster = getCluster(cid);
        if (cluster) {
            stage_clusters.push_back(cluster);
        }
    }
    
    if (stage_clusters.empty()) {
        throw std::runtime_error("No clusters for stage " + std::to_string(stage_id));
    }
    
    // Get pump schedule
    auto it = pump_schedules_.find(stage_id);
    if (it == pump_schedules_.end()) {
        throw std::runtime_error("No pump schedule for stage " + std::to_string(stage_id));
    }
    const auto& schedule = it->second;
    
    // Initialize result vectors
    result.cluster_volumes.resize(stage_clusters.size(), 0.0);
    result.cluster_proppant.resize(stage_clusters.size(), 0.0);
    result.cluster_efficiency.resize(stage_clusters.size(), 0.0);
    
    // Simulation state
    double time = 0.0;
    double total_volume = 0.0;
    double total_proppant = 0.0;
    double half_length = 0.0;
    double height = 50.0; // Initial assumed height (m)
    double width = 0.001; // Initial width (m)
    
    // Rock properties at stage depth
    double E = interpolateRock(stage_clusters[0]->measured_depth, rock_depths_, youngs_modulus_);
    double nu = interpolateRock(stage_clusters[0]->measured_depth, rock_depths_, poisson_ratio_);
    double E_prime = E / (1.0 - nu * nu);
    
    double closure = stage_clusters[0]->min_horizontal_stress;
    double leakoff_coeff = 5e-5; // Default Carter leakoff
    
    // Simulate each schedule step
    for (const auto& step : schedule) {
        FracFluid fluid;
        if (fluids_.count(step.fluid_name)) {
            fluid = fluids_[step.fluid_name];
            leakoff_coeff = fluid.leakoff_coefficient;
        }
        
        double step_time = 0.0;
        while (step_time < step.duration) {
            // Calculate flow distribution
            auto [efficiency, rates] = limited_entry_.predictClusterEfficiency(
                std::vector<PerforationCluster>(stage_clusters.begin(), stage_clusters.end()),
                step.rate, fluid);
            
            // Update each cluster
            for (size_t i = 0; i < stage_clusters.size(); ++i) {
                double cluster_rate = rates[i];
                result.cluster_volumes[i] += cluster_rate * dt;
                result.cluster_proppant[i] += cluster_rate * step.proppant_concentration * dt;
                
                stage_clusters[i]->cumulative_volume += cluster_rate * dt;
                stage_clusters[i]->proppant_placed += cluster_rate * step.proppant_concentration * dt;
            }
            
            // Update fracture geometry (simplified P3D)
            total_volume += step.rate * dt;
            total_proppant += step.rate * step.proppant_concentration * dt;
            
            // Leakoff volume
            double leakoff_volume = 2.0 * leakoff_coeff * std::sqrt(time) * 2.0 * half_length * height;
            double frac_volume = total_volume - leakoff_volume;
            
            // Estimate geometry from volume and pressure
            // Simplified: V = w * L * H
            double net_pressure = 2e6; // Simplified net pressure assumption
            width = 4.0 * net_pressure * height / E_prime; // PKN width
            
            if (frac_volume > 0 && width > 0) {
                half_length = frac_volume / (2.0 * width * height);
            }
            
            // Record history
            result.time_history.push_back(time);
            result.rate_history.push_back(step.rate);
            result.pressure_history.push_back(closure + net_pressure);
            
            time += dt;
            step_time += dt;
        }
    }
    
    // Calculate final results
    result.final_half_length = half_length;
    result.final_height = height;
    
    // Propped geometry (assume 80% of created)
    result.propped_length = 0.8 * half_length;
    result.propped_height = 0.9 * height;
    
    // Leakoff
    result.total_leakoff = leakoff_coeff * std::sqrt(time) * 2.0 * half_length * height;
    
    // Cluster efficiency
    double total_cluster_volume = std::accumulate(result.cluster_volumes.begin(),
                                                   result.cluster_volumes.end(), 0.0);
    for (size_t i = 0; i < result.cluster_efficiency.size(); ++i) {
        result.cluster_efficiency[i] = result.cluster_volumes[i] / 
                                       (total_cluster_volume / result.cluster_efficiency.size());
    }
    
    // Pressures
    if (!result.pressure_history.empty()) {
        result.isip = result.pressure_history.back();
        result.closure_pressure = closure;
        result.net_pressure = result.isip - closure;
    }
    
    return result;
}

std::vector<MultiStageFracManager::StageSimResult>
MultiStageFracManager::simulateFullTreatment(double dt) {
    std::vector<StageSimResult> results;
    
    // Sort stages by ID (toe to heel)
    std::vector<int> stage_order;
    for (const auto& stage : stages_) {
        stage_order.push_back(stage.stage_id);
    }
    std::sort(stage_order.begin(), stage_order.end());
    
    for (int stage_id : stage_order) {
        // Update stress shadowing before each stage
        updateStressShadowing();
        
        // Simulate stage
        auto result = simulateStage(stage_id, dt);
        results.push_back(result);
    }
    
    return results;
}

void MultiStageFracManager::updateStressShadowing() {
    // Simplified stress shadowing calculation
    // For each cluster, add stress from all existing fractures
    
    for (auto& cluster : clusters_) {
        double delta_stress = 0.0;
        
        // Calculate stress shadow from all other clusters that have been treated
        for (const auto& other : clusters_) {
            if (&other == &cluster) continue;
            if (other.cumulative_volume < 1.0) continue; // Not yet treated
            
            // Distance between clusters
            double distance = std::abs(cluster.measured_depth - other.measured_depth);
            
            // Simplified stress shadow model
            // Delta_Shmin ~ P_net * (w/d)
            double net_pressure = 2e6; // Assumed
            double frac_width = 0.005; // Assumed 5mm
            
            if (distance > 0.1) { // Avoid division by zero
                delta_stress += net_pressure * frac_width / distance;
            }
        }
        
        cluster.stress_shadow_delta = delta_stress;
    }
}

TreatmentStage* MultiStageFracManager::getStage(int stage_id) {
    for (auto& stage : stages_) {
        if (stage.stage_id == stage_id) {
            return &stage;
        }
    }
    return nullptr;
}

PerforationCluster* MultiStageFracManager::getCluster(int cluster_id) {
    for (auto& cluster : clusters_) {
        if (cluster.cluster_id == cluster_id) {
            return &cluster;
        }
    }
    return nullptr;
}

double MultiStageFracManager::interpolateStress(double depth, 
                                                 const std::vector<double>& depths,
                                                 const std::vector<double>& values) const {
    if (depths.empty() || values.empty()) return 30e6; // Default
    if (depth <= depths.front()) return values.front();
    if (depth >= depths.back()) return values.back();
    
    for (size_t i = 0; i < depths.size() - 1; ++i) {
        if (depth >= depths[i] && depth <= depths[i+1]) {
            double t = (depth - depths[i]) / (depths[i+1] - depths[i]);
            return values[i] + t * (values[i+1] - values[i]);
        }
    }
    
    return values.back();
}

double MultiStageFracManager::interpolateRock(double depth,
                                               const std::vector<double>& depths,
                                               const std::vector<double>& values) const {
    return interpolateStress(depth, depths, values); // Same interpolation logic
}

void MultiStageFracManager::writeStageSummary(std::ostream& os, int stage_id) const {
    const TreatmentStage* stage = nullptr;
    for (const auto& s : stages_) {
        if (s.stage_id == stage_id) {
            stage = &s;
            break;
        }
    }
    
    if (!stage) {
        os << "Stage " << stage_id << " not found\n";
        return;
    }
    
    os << "========================================\n";
    os << "Stage " << stage_id << " Summary\n";
    os << "========================================\n";
    os << "Number of clusters: " << stage->cluster_ids.size() << "\n";
    os << "Total volume: " << stage->total_volume << " m³\n";
    os << "Total proppant: " << stage->total_proppant << " kg\n";
    os << "ISIP: " << stage->isip / 1e6 << " MPa\n";
    os << "Closure pressure: " << stage->closure_pressure / 1e6 << " MPa\n";
    os << "Fracture half-length: " << stage->fracture_half_length << " m\n";
    os << "Fracture height: " << stage->fracture_height << " m\n";
    os << "Cluster efficiency: " << stage->cluster_efficiency * 100 << " %\n";
}

// ============================================================================
// PerforationErosion Implementation
// ============================================================================

PerforationErosion::PerforationErosion()
    : erosion_coeff_(1e-11), hardness_exp_(1.0), velocity_exp_(2.0)
{}

double PerforationErosion::calculateErosionRate(double velocity, double proppant_conc,
                                                  double proppant_hardness, double time) const {
    // API RP 19B erosion model
    // Erosion rate = C * V^n * conc * hardness^m
    double rate = erosion_coeff_ * std::pow(velocity, velocity_exp_) * 
                 proppant_conc * std::pow(proppant_hardness, hardness_exp_);
    
    return rate;
}

void PerforationErosion::updatePerforation(PerforationCluster& cluster, double rate,
                                            double proppant_conc, double dt) {
    cluster.updateErosion(rate, proppant_conc, dt);
}

// ============================================================================
// IsolationSimulator Implementation
// ============================================================================

IsolationSimulator::IsolationSimulator()
    : ball_density_(1200.0), ball_cd_(0.44)
{}

std::pair<double, double> IsolationSimulator::simulateBallDrop(
    double ball_diameter,
    double seat_diameter,
    double flow_rate,
    double fluid_viscosity,
    const std::vector<double>& well_md,
    const std::vector<double>& well_inc)
{
    if (well_md.empty() || well_inc.empty()) {
        return {0.0, 0.0};
    }
    
    double position = 0.0;
    double time = 0.0;
    double dt = 0.1; // Time step
    
    // Simulate ball falling
    while (position < well_md.back()) {
        // Find current inclination
        size_t idx = 0;
        for (size_t i = 0; i < well_md.size() - 1; ++i) {
            if (position >= well_md[i] && position <= well_md[i+1]) {
                idx = i;
                break;
            }
        }
        
        double inc_rad = well_inc[idx] * M_PI / 180.0;
        double fluid_density = 1000.0; // Assumed
        
        // Ball velocity
        double v_ball = calculateBallVelocity(ball_diameter, fluid_density,
                                               fluid_viscosity, inc_rad);
        
        // Account for fluid flow
        double pipe_area = M_PI * 0.1 * 0.1 / 4.0; // 10cm pipe
        double fluid_velocity = flow_rate / pipe_area;
        
        double net_velocity = v_ball + fluid_velocity * std::cos(inc_rad);
        
        position += net_velocity * dt;
        time += dt;
        
        // Check if ball reached seat
        // (This is simplified - would need actual seat depths)
    }
    
    return {time, position};
}

double IsolationSimulator::calculateBallVelocity(double diameter, double fluid_density,
                                                   double fluid_viscosity, double inclination) const {
    // Terminal velocity considering buoyancy and drag
    double g = 9.81;
    double delta_rho = ball_density_ - fluid_density;
    
    // Stokes settling adjusted for inclination
    double v_terminal = std::sqrt(4.0 * g * diameter * delta_rho / (3.0 * ball_cd_ * fluid_density));
    
    // Only vertical component contributes to downward motion
    return v_terminal * std::cos(inclination);
}

std::pair<bool, double> IsolationSimulator::simulatePlugSet(
    double plug_depth,
    double setting_pressure,
    double wellbore_pressure)
{
    // Check if pressure sufficient to set plug
    if (wellbore_pressure >= setting_pressure) {
        return {true, plug_depth};
    }
    
    return {false, 0.0};
}

} // namespace FSRM
