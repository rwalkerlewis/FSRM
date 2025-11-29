#ifndef FIBER_OPTIC_INTEGRATION_HPP
#define FIBER_OPTIC_INTEGRATION_HPP

/**
 * @file FiberOpticIntegration.hpp
 * @brief Fiber optic diagnostics integration (DAS/DTS)
 * 
 * ResFrac-equivalent capabilities:
 * - Distributed Acoustic Sensing (DAS) analysis
 * - Distributed Temperature Sensing (DTS) analysis
 * - Cluster efficiency from DAS
 * - Fracture height from DTS warm-back
 * - Flow allocation monitoring
 * - Cross-well strain monitoring
 * - Microseismic integration
 */

#include <vector>
#include <memory>
#include <string>
#include <map>
#include <array>
#include <functional>
#include <complex>
#include <cmath>

namespace FSRM {

/**
 * @brief DAS data point
 */
struct DASDataPoint {
    double time;                        ///< Time (s)
    double depth;                       ///< Measured depth (m)
    double strain_rate;                 ///< Strain rate (nε/s)
    double amplitude;                   ///< Raw amplitude
    double phase;                       ///< Phase (radians)
    
    DASDataPoint();
};

/**
 * @brief DTS data point
 */
struct DTSDataPoint {
    double time;                        ///< Time (s)
    double depth;                       ///< Measured depth (m)
    double temperature;                 ///< Temperature (K)
    double stokes;                      ///< Stokes intensity
    double anti_stokes;                 ///< Anti-Stokes intensity
    
    DTSDataPoint();
};

/**
 * @brief Fiber specification
 */
struct FiberSpecification {
    std::string fiber_type;             ///< "SINGLE_MODE", "MULTI_MODE"
    double gauge_length;                ///< DAS gauge length (m)
    double spatial_sampling;            ///< Spatial sampling interval (m)
    double temporal_sampling;           ///< Temporal sampling rate (Hz)
    double depth_offset;                ///< Offset from wellhead (m)
    
    // DAS properties
    double sensitivity;                 ///< Strain sensitivity (nε/s)
    double noise_floor;                 ///< Noise floor (nε/s)
    
    // DTS properties
    double temperature_resolution;      ///< Temperature resolution (K)
    double spatial_resolution;          ///< Spatial resolution (m)
    
    FiberSpecification();
};

/**
 * @brief Cluster identification from DAS
 */
struct DASClusterResult {
    int cluster_id;                     ///< Cluster identifier
    double depth_min;                   ///< Minimum depth (m)
    double depth_max;                   ///< Maximum depth (m)
    double depth_center;                ///< Center depth (m)
    double relative_activity;           ///< Relative acoustic activity (0-1)
    double efficiency;                  ///< Estimated cluster efficiency (0-1)
    double dominant_frequency;          ///< Dominant frequency (Hz)
    double duration;                    ///< Activity duration (s)
    
    DASClusterResult();
};

/**
 * @brief Temperature anomaly from DTS
 */
struct DTSAnomaly {
    double depth_top;                   ///< Top of anomaly (m)
    double depth_bottom;                ///< Bottom of anomaly (m)
    double temperature_change;          ///< Temperature change from baseline (K)
    double anomaly_type;                ///< "COOLING" or "HEATING"
    double confidence;                  ///< Detection confidence (0-1)
    
    DTSAnomaly();
};

/**
 * @brief DAS waterfall analyzer
 */
class DASWaterfallAnalyzer {
public:
    DASWaterfallAnalyzer();
    
    /**
     * @brief Set fiber specification
     */
    void setFiberSpec(const FiberSpecification& spec);
    
    /**
     * @brief Load DAS data
     * 
     * @param data 2D array [time_idx][depth_idx] of strain rate
     * @param times Time values (s)
     * @param depths Depth values (m)
     */
    void loadData(const std::vector<std::vector<double>>& data,
                  const std::vector<double>& times,
                  const std::vector<double>& depths);
    
    /**
     * @brief Set cluster locations (from completion)
     * 
     * @param cluster_depths Measured depths of clusters (m)
     */
    void setClusterLocations(const std::vector<double>& cluster_depths);
    
    /**
     * @brief Identify active clusters
     * 
     * @param time_window Time window for analysis (s)
     * @return Vector of cluster results
     */
    std::vector<DASClusterResult> identifyActiveClusters(double time_window);
    
    /**
     * @brief Calculate cluster efficiency
     * 
     * Based on relative acoustic activity at each cluster.
     * 
     * @return Vector of efficiencies (0-1)
     */
    std::vector<double> calculateClusterEfficiency();
    
    /**
     * @brief Detect fracture hits
     * 
     * Sudden acoustic activity at unexpected depths.
     * 
     * @return Depths of suspected frac hits (m)
     */
    std::vector<double> detectFracHits();
    
    /**
     * @brief Calculate frequency spectrum
     * 
     * @param depth Depth for analysis (m)
     * @param start_time Start time (s)
     * @param end_time End time (s)
     * @return {frequencies, power spectral density}
     */
    std::pair<std::vector<double>, std::vector<double>> frequencySpectrum(
        double depth, double start_time, double end_time);
    
    /**
     * @brief Generate time-depth waterfall plot data
     * 
     * @param freq_min Minimum frequency for bandpass (Hz)
     * @param freq_max Maximum frequency for bandpass (Hz)
     * @return Filtered waterfall data
     */
    std::vector<std::vector<double>> getFilteredWaterfall(
        double freq_min, double freq_max);
    
    /**
     * @brief Calculate RMS energy over time
     * 
     * @param depth_min Minimum depth (m)
     * @param depth_max Maximum depth (m)
     * @return {times, RMS energy}
     */
    std::pair<std::vector<double>, std::vector<double>> calculateRMSEnergy(
        double depth_min, double depth_max);
    
private:
    FiberSpecification fiber_spec_;
    std::vector<std::vector<double>> data_;
    std::vector<double> times_;
    std::vector<double> depths_;
    std::vector<double> cluster_depths_;
    
    // Signal processing
    void bandpassFilter(std::vector<double>& signal, 
                        double freq_min, double freq_max);
    std::vector<std::complex<double>> fft(const std::vector<double>& signal);
    double calculateRMS(const std::vector<double>& signal);
};

/**
 * @brief DTS analyzer
 */
class DTSAnalyzer {
public:
    DTSAnalyzer();
    
    /**
     * @brief Set fiber specification
     */
    void setFiberSpec(const FiberSpecification& spec);
    
    /**
     * @brief Load DTS data
     * 
     * @param data 2D array [time_idx][depth_idx] of temperatures
     * @param times Time values (s)
     * @param depths Depth values (m)
     */
    void loadData(const std::vector<std::vector<double>>& data,
                  const std::vector<double>& times,
                  const std::vector<double>& depths);
    
    /**
     * @brief Set baseline (pre-treatment) temperature
     * 
     * @param baseline Temperature profile before treatment (K)
     */
    void setBaseline(const std::vector<double>& baseline);
    
    /**
     * @brief Detect temperature anomalies
     * 
     * @param threshold_temperature Temperature change threshold (K)
     * @return Vector of anomalies
     */
    std::vector<DTSAnomaly> detectAnomalies(double threshold_temperature = 1.0);
    
    /**
     * @brief Analyze warm-back after shut-in
     * 
     * Used to estimate fracture height from temperature recovery.
     * 
     * @param shut_in_time Time of shut-in (s)
     * @param analysis_duration Duration to analyze (s)
     * @return Estimated fracture height per cluster (m)
     */
    std::vector<double> analyzeWarmBack(double shut_in_time,
                                         double analysis_duration);
    
    /**
     * @brief Estimate injection depth from cooling
     * 
     * @param injection_start Start of injection (s)
     * @return Depth range of injection (m)
     */
    std::pair<double, double> estimateInjectionDepth(double injection_start);
    
    /**
     * @brief Calculate flow allocation from temperature
     * 
     * Using cooling rate at each cluster.
     * 
     * @param cluster_depths Cluster measured depths (m)
     * @return Flow allocation fraction for each cluster
     */
    std::vector<double> calculateFlowAllocation(
        const std::vector<double>& cluster_depths);
    
    /**
     * @brief Get temperature profile at time
     * 
     * @param time Query time (s)
     * @return Temperature profile (K)
     */
    std::vector<double> getProfileAtTime(double time);
    
    /**
     * @brief Get temperature difference from baseline
     * 
     * @param time Query time (s)
     * @return Temperature difference profile (K)
     */
    std::vector<double> getTemperatureDelta(double time);
    
private:
    FiberSpecification fiber_spec_;
    std::vector<std::vector<double>> data_;
    std::vector<double> times_;
    std::vector<double> depths_;
    std::vector<double> baseline_;
    
    // Interpolation
    double interpolateTemperature(double time, double depth);
    
    // Warm-back analysis
    double fitWarmBackExponential(const std::vector<double>& times,
                                   const std::vector<double>& temps);
};

/**
 * @brief Production DAS analysis
 */
class ProductionDASAnalyzer {
public:
    ProductionDASAnalyzer();
    
    /**
     * @brief Analyze flow contribution from DAS
     * 
     * During production, DAS can indicate flow from different depths.
     * 
     * @param data DAS data during production
     * @param times Time values (s)
     * @param depths Depth values (m)
     * @return Flow contribution per depth interval
     */
    std::vector<double> analyzeFlowContribution(
        const std::vector<std::vector<double>>& data,
        const std::vector<double>& times,
        const std::vector<double>& depths);
    
    /**
     * @brief Detect water breakthrough
     * 
     * Changes in acoustic signature may indicate water production.
     * 
     * @param current_data Current DAS data
     * @param baseline_data Baseline (dry gas/oil) DAS data
     * @return Depths with suspected water
     */
    std::vector<double> detectWaterBreakthrough(
        const std::vector<std::vector<double>>& current_data,
        const std::vector<std::vector<double>>& baseline_data);
    
    /**
     * @brief Calculate production allocation
     * 
     * @param cluster_depths Cluster depths (m)
     * @return Production allocation per cluster (0-1)
     */
    std::vector<double> calculateProductionAllocation(
        const std::vector<double>& cluster_depths);
    
private:
    FiberSpecification fiber_spec_;
};

/**
 * @brief Cross-well DAS analysis
 */
class CrossWellDASAnalyzer {
public:
    CrossWellDASAnalyzer();
    
    /**
     * @brief Set monitoring well fiber
     */
    void setMonitoringWell(const std::string& well_name,
                           const FiberSpecification& spec);
    
    /**
     * @brief Load monitoring well DAS data
     */
    void loadMonitoringData(const std::vector<std::vector<double>>& data,
                            const std::vector<double>& times,
                            const std::vector<double>& depths);
    
    /**
     * @brief Set treatment well location
     */
    void setTreatmentWell(double heel_x, double heel_y, double heel_z,
                          double toe_x, double toe_y, double toe_z);
    
    /**
     * @brief Detect strain from treatment well
     * 
     * @return Strain events detected in monitoring well
     */
    struct StrainEvent {
        double time;                    ///< Detection time (s)
        double depth;                   ///< Detection depth in monitor well (m)
        double strain_magnitude;        ///< Strain magnitude (nε)
        double estimated_distance;      ///< Estimated distance to source (m)
        double azimuth;                 ///< Estimated azimuth to source (degrees)
    };
    
    std::vector<StrainEvent> detectStrainEvents(double threshold);
    
    /**
     * @brief Estimate SRV from strain
     * 
     * Using extent of strain signals in monitoring well.
     * 
     * @return Estimated SRV dimensions {length, width, height}
     */
    std::array<double, 3> estimateSRV();
    
    /**
     * @brief Detect frac hit timing
     * 
     * @return Time and depth of suspected frac hit
     */
    std::pair<double, double> detectFracHit();
    
private:
    FiberSpecification monitor_spec_;
    std::vector<std::vector<double>> monitor_data_;
    std::vector<double> monitor_times_;
    std::vector<double> monitor_depths_;
    
    double treatment_heel_[3];
    double treatment_toe_[3];
    
    // Distance calculation
    double calculateDistance(double monitor_depth, double treatment_md);
};

/**
 * @brief Microseismic event
 */
struct MicroseismicEvent {
    double time;                        ///< Origin time (s)
    double x, y, z;                     ///< Location (m)
    double magnitude;                   ///< Moment magnitude
    double azimuth;                     ///< Strike azimuth (degrees)
    double dip;                         ///< Dip angle (degrees)
    std::string mechanism;              ///< Focal mechanism ("strike-slip", "normal", etc.)
    double location_uncertainty;        ///< Location uncertainty (m)
    
    MicroseismicEvent();
};

/**
 * @brief Microseismic analysis
 */
class MicroseismicAnalyzer {
public:
    MicroseismicAnalyzer();
    
    /**
     * @brief Load microseismic events
     */
    void loadEvents(const std::vector<MicroseismicEvent>& events);
    
    /**
     * @brief Set well trajectory
     */
    void setWellTrajectory(const std::vector<double>& md,
                           const std::vector<double>& x,
                           const std::vector<double>& y,
                           const std::vector<double>& z);
    
    /**
     * @brief Estimate SRV from events
     * 
     * @param hull_method "CONVEX_HULL", "ELLIPSOID", "RECTANGULAR"
     * @return SRV volume (m³)
     */
    double estimateSRV(const std::string& hull_method = "CONVEX_HULL");
    
    /**
     * @brief Get SRV dimensions
     * 
     * @return {length, width, height} in meters
     */
    std::array<double, 3> getSRVDimensions();
    
    /**
     * @brief Calculate b-value (Gutenberg-Richter)
     * 
     * @return b-value
     */
    double calculateBValue();
    
    /**
     * @brief Identify event clusters
     * 
     * @param min_cluster_size Minimum events per cluster
     * @param eps DBSCAN epsilon parameter (m)
     * @return Cluster assignments for each event
     */
    std::vector<int> clusterEvents(int min_cluster_size = 5,
                                    double eps = 50.0);
    
    /**
     * @brief Calculate fracture azimuth from events
     * 
     * @return Dominant azimuth (degrees)
     */
    double calculateDominantAzimuth();
    
    /**
     * @brief Get event density map
     * 
     * @param cell_size Grid cell size (m)
     * @return 3D density map
     */
    std::vector<std::vector<std::vector<int>>> getEventDensityMap(
        double cell_size);
    
    /**
     * @brief Correlate events with treatment
     * 
     * @param treatment_times Treatment time points (s)
     * @param treatment_pressures Pressures (Pa)
     * @param treatment_rates Rates (m³/s)
     * @return Correlation results
     */
    struct CorrelationResult {
        double pressure_event_correlation;
        double rate_event_correlation;
        double lag_time;                ///< Time lag for max correlation
    };
    
    CorrelationResult correlateWithTreatment(
        const std::vector<double>& treatment_times,
        const std::vector<double>& treatment_pressures,
        const std::vector<double>& treatment_rates);
    
private:
    std::vector<MicroseismicEvent> events_;
    std::vector<double> well_md_, well_x_, well_y_, well_z_;
    
    // Convex hull calculation
    std::vector<std::array<double, 3>> calculateConvexHull();
    
    // Ellipsoid fitting
    void fitEllipsoid(double& a, double& b, double& c);
};

/**
 * @brief Integrated fiber optic manager
 */
class FiberOpticManager {
public:
    FiberOpticManager();
    
    // Setup
    void setDASAnalyzer(std::shared_ptr<DASWaterfallAnalyzer> analyzer);
    void setDTSAnalyzer(std::shared_ptr<DTSAnalyzer> analyzer);
    void setMicroseismicAnalyzer(std::shared_ptr<MicroseismicAnalyzer> analyzer);
    
    /**
     * @brief Set cluster locations
     */
    void setClusterLocations(const std::vector<double>& cluster_depths);
    
    /**
     * @brief Generate integrated completion diagnostics
     */
    struct CompletionDiagnostics {
        std::vector<double> cluster_efficiency_das;
        std::vector<double> cluster_efficiency_dts;
        std::vector<double> flow_allocation;
        std::vector<double> fracture_height;
        double estimated_srv;
        std::vector<double> frac_hit_depths;
        std::string overall_assessment;
    };
    
    CompletionDiagnostics generateDiagnostics();
    
    /**
     * @brief Export integrated report
     */
    void exportReport(const std::string& filename);
    
private:
    std::shared_ptr<DASWaterfallAnalyzer> das_analyzer_;
    std::shared_ptr<DTSAnalyzer> dts_analyzer_;
    std::shared_ptr<MicroseismicAnalyzer> microseismic_analyzer_;
    
    std::vector<double> cluster_depths_;
};

} // namespace FSRM

#endif // FIBER_OPTIC_INTEGRATION_HPP
