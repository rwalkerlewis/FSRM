/**
 * @file SeismicSource.hpp
 * @brief Seismic source models for earthquake simulation
 * 
 * Implements various seismic source representations used in SeisSol:
 * - Point sources (moment tensor, single force)
 * - Kinematic sources (finite faults with prescribed slip)
 * - Dynamic rupture sources (coupled with friction)
 * - Moment rate functions (Gaussian, regularized Yoffe, etc.)
 * 
 * Also includes receiver models for recording waveforms at specified locations.
 */

#ifndef SEISMIC_SOURCE_HPP
#define SEISMIC_SOURCE_HPP

#include <vector>
#include <array>
#include <string>
#include <map>
#include <memory>
#include <functional>
#include <complex>

namespace FSRM {

// Forward declarations
class DiscontinuousGalerkin;

/**
 * @brief Moment tensor representation
 * 
 * Full 3x3 symmetric moment tensor:
 * | M_xx  M_xy  M_xz |
 * | M_xy  M_yy  M_yz |
 * | M_xz  M_yz  M_zz |
 */
struct MomentTensor {
    double Mxx, Myy, Mzz;  // Diagonal components
    double Mxy, Mxz, Myz;  // Off-diagonal components
    
    MomentTensor() : Mxx(0), Myy(0), Mzz(0), Mxy(0), Mxz(0), Myz(0) {}
    
    MomentTensor(double mxx, double myy, double mzz, 
                 double mxy, double mxz, double myz)
        : Mxx(mxx), Myy(myy), Mzz(mzz), Mxy(mxy), Mxz(mxz), Myz(myz) {}
    
    // Create from strike, dip, rake and scalar moment
    static MomentTensor fromFaultGeometry(double strike, double dip, double rake, double M0);
    
    // Create double-couple source
    static MomentTensor doubleCouple(double strike, double dip, double rake, double M0);
    
    // Create explosion/implosion source
    static MomentTensor explosion(double M0);
    
    // Create CLVD (Compensated Linear Vector Dipole)
    static MomentTensor clvd(double M0, double strike, double dip);
    
    // Scalar moment (√(Σ M_ij² / 2))
    double scalarMoment() const;
    
    // Moment magnitude Mw
    double magnitude() const;
    
    // Decompose into isotropic, CLVD, and double-couple components
    void decompose(double& M0_iso, double& M0_clvd, double& M0_dc,
                  double& strike, double& dip, double& rake) const;
    
    // Get as array
    void toArray(double* M) const;  // Length 6
    void toMatrix(double* M) const; // 3x3 matrix
    
    // Operators
    MomentTensor operator+(const MomentTensor& other) const;
    MomentTensor operator*(double scale) const;
};

/**
 * @brief Source time function types
 */
enum class SourceTimeFunction {
    GAUSSIAN,           // Gaussian pulse
    RICKER,             // Ricker (Mexican hat) wavelet
    STEP,               // Heaviside step function
    RAMP,               // Linear ramp
    TRIANGLE,           // Triangular pulse
    SMOOTHED_RAMP,      // Smoothed ramp (regularized)
    YOFFE,              // Regularized Yoffe function
    BRUNE,              // Brune ω² model
    CUSTOM              // User-defined function
};

/**
 * @brief Source time function evaluator
 * 
 * Provides moment rate M'(t) for various source time functions.
 * The moment at time t is M(t) = ∫₀ᵗ M'(τ) dτ
 */
class SourceTimeFunctionEvaluator {
public:
    SourceTimeFunctionEvaluator(SourceTimeFunction type = SourceTimeFunction::GAUSSIAN);
    
    // Set parameters
    void setType(SourceTimeFunction type);
    void setDuration(double T);      // Characteristic duration
    void setRiseTime(double tau);    // Rise time for ramp functions
    void setOnsetTime(double t0);    // Onset time
    void setPeakTime(double tp);     // Time of peak moment rate
    void setFrequency(double f0);    // Dominant frequency (for Ricker)
    
    // Custom function
    void setCustomFunction(std::function<double(double)> func);
    
    // Evaluate moment rate M'(t)
    double momentRate(double t) const;
    
    // Evaluate moment M(t) = ∫₀ᵗ M'(τ) dτ
    double moment(double t) const;
    
    // Normalized moment rate (peak = 1)
    double normalizedRate(double t) const;
    
    // Get frequency content (dominant frequency)
    double getDominantFrequency() const;
    
    // Get corner frequency (for Brune model)
    double getCornerFrequency() const;
    
private:
    SourceTimeFunction type;
    double duration;       // T or σ for Gaussian
    double rise_time;      // τ
    double onset_time;     // t₀
    double peak_time;      // t_p
    double frequency;      // f₀
    
    std::function<double(double)> custom_func;
    
    // Specific implementations
    double gaussianRate(double t) const;
    double rickerRate(double t) const;
    double stepRate(double t) const;
    double rampRate(double t) const;
    double triangleRate(double t) const;
    double smoothedRampRate(double t) const;
    double yoffeRate(double t) const;
    double bruneRate(double t) const;
};

/**
 * @brief Point source (moment tensor source at a single location)
 * 
 * Represents a point source with:
 * - Location (x, y, z)
 * - Moment tensor M_ij
 * - Source time function s(t)
 * 
 * The source term in the equations is:
 * f_i = -∂j[M_ij * s(t) * δ(x - x_s)]
 */
class PointSource {
public:
    PointSource();
    
    // Configuration
    void setLocation(double x, double y, double z);
    void setMomentTensor(const MomentTensor& M);
    void setSourceTimeFunction(const SourceTimeFunctionEvaluator& stf);
    void setScalarMoment(double M0);
    
    // Configure from fault parameters
    void setFromFault(double strike, double dip, double rake, double M0,
                     double x, double y, double z);
    
    // Get location
    void getLocation(double& x, double& y, double& z) const;
    std::array<double, 3> getLocation() const;
    
    // Evaluate source term at time t
    MomentTensor getMomentTensor(double t) const;
    double getScalarMoment(double t) const;
    double getMomentRate(double t) const;
    
    // Source term for wave equation (9 components for elastic)
    void getSourceTerm(double t, double* f) const;
    
    // Finite-extent source approximation (Gaussian smoothing)
    void setSpatialSmoothing(double sigma);
    double getSpatialWeight(double x, double y, double z) const;
    
    // Is source active at time t?
    bool isActive(double t) const;
    double getOnsetTime() const;
    double getEndTime() const;
    
    std::string name;
    
private:
    double loc_x, loc_y, loc_z;
    MomentTensor moment_tensor;
    SourceTimeFunctionEvaluator stf;
    double scalar_moment;
    double spatial_sigma;
    bool use_spatial_smoothing;
};

/**
 * @brief Subfault for kinematic source model
 * 
 * Represents a small patch of a finite fault with:
 * - Location and orientation
 * - Slip amplitude, rake, and rupture time
 * - Rise time and slip rate function
 */
struct KinematicSubfault {
    // Location (center of subfault)
    double x, y, z;
    
    // Dimensions
    double length, width;
    double area;
    
    // Orientation
    double strike, dip;
    
    // Slip parameters
    double slip;           // Final slip amplitude (m)
    double rake;           // Slip direction (radians)
    double rupture_time;   // Rupture arrival time (s)
    double rise_time;      // Slip rise time (s)
    
    // Optional: slip rate function parameters
    double peak_slip_rate; // Peak slip rate (m/s)
    int slip_function;     // 0=triangle, 1=regularized Yoffe, 2=Brune
    
    // Derived quantities
    double seismic_moment; // M0 = μ * A * D
    MomentTensor moment_tensor;
    
    KinematicSubfault() :
        x(0), y(0), z(0), length(100), width(100), area(10000),
        strike(0), dip(M_PI/2), slip(1.0), rake(0), rupture_time(0),
        rise_time(1.0), peak_slip_rate(0), slip_function(0), seismic_moment(0) {}
    
    // Compute moment tensor from slip
    void computeMomentTensor(double shear_modulus);
    
    // Get slip at time t
    double getSlip(double t) const;
    double getSlipRate(double t) const;
};

/**
 * @brief Kinematic source model (finite fault)
 * 
 * Represents a finite fault discretized into subfaults.
 * Each subfault has prescribed slip as a function of time.
 * Used for:
 * - Earthquake source studies
 * - Ground motion simulation
 * - SCEC TPV benchmarks
 */
class KinematicSource {
public:
    KinematicSource();
    
    // Create fault geometry
    void createRectangularFault(double center_x, double center_y, double center_z,
                               double strike, double dip,
                               double length, double width,
                               int n_along_strike, int n_down_dip);
    
    // Load from file (SRF format or similar)
    void loadFromSRF(const std::string& filename);
    void loadFromFSP(const std::string& filename);
    
    // Set rupture parameters
    void setUniformSlip(double slip, double rake);
    void setCircularRupture(double hypo_x, double hypo_y, double hypo_z,
                           double rupture_velocity);
    void setSlipFromFile(const std::string& filename);
    
    // Set rise time
    void setUniformRiseTime(double rise_time);
    void setVariableRiseTime(std::function<double(double, double, double)> func);
    
    // Set shear modulus (for moment computation)
    void setShearModulus(double mu);
    
    // Get total seismic moment
    double getTotalMoment() const;
    double getMagnitude() const;
    
    // Get moment rate at time t
    double getMomentRate(double t) const;
    
    // Get source contribution at a point
    void getSourceTerm(double x, double y, double z, double t, double* f) const;
    
    // Get subfaults
    const std::vector<KinematicSubfault>& getSubfaults() const { return subfaults; }
    size_t getNumSubfaults() const { return subfaults.size(); }
    
    // Time range
    double getStartTime() const;
    double getEndTime() const;
    
    std::string name;
    
private:
    std::vector<KinematicSubfault> subfaults;
    double shear_modulus;
    
    // Hypocenter
    double hypo_x, hypo_y, hypo_z;
    
    // Spatial smoothing for distributed sources
    double spatial_sigma;
};

/**
 * @brief Single force source (for comparison with analytic solutions)
 */
class SingleForceSource {
public:
    SingleForceSource();
    
    void setLocation(double x, double y, double z);
    void setForce(double fx, double fy, double fz);
    void setSourceTimeFunction(const SourceTimeFunctionEvaluator& stf);
    void setAmplitude(double F0);
    
    void getForce(double t, double& fx, double& fy, double& fz) const;
    
private:
    double loc_x, loc_y, loc_z;
    double force_x, force_y, force_z;
    double amplitude;
    SourceTimeFunctionEvaluator stf;
};

/**
 * @brief Receiver for recording seismic waveforms
 * 
 * Records displacement, velocity, or stress at specified locations.
 * Supports output in various formats (ASCII, SAC, etc.)
 */
class SeismicReceiver {
public:
    enum class RecordType {
        DISPLACEMENT,
        VELOCITY,
        ACCELERATION,
        STRESS,
        STRAIN,
        ROTATION
    };
    
    SeismicReceiver();
    
    // Configuration
    void setLocation(double x, double y, double z);
    void setName(const std::string& name);
    void setRecordType(RecordType type);
    void setSamplingRate(double dt);
    
    // Record data
    void record(double t, const double* u);
    void recordVelocity(double t, const double* v);
    void recordStress(double t, const double* sigma);
    
    // Get recorded data
    const std::vector<double>& getTimes() const { return times; }
    const std::vector<std::array<double, 3>>& getData() const { return data; }
    
    // Output
    void writeASCII(const std::string& filename) const;
    void writeSAC(const std::string& filename) const;
    void writeHDF5(const std::string& filename) const;
    
    // Get peak values
    double getPeakValue(int component = -1) const;
    double getPeakGroundVelocity() const;
    double getPeakGroundAcceleration() const;
    
    // Location
    void getLocation(double& x, double& y, double& z) const;
    
    std::string name;
    
private:
    double loc_x, loc_y, loc_z;
    RecordType record_type;
    double sampling_dt;
    double last_record_time;
    
    std::vector<double> times;
    std::vector<std::array<double, 3>> data;
    std::vector<std::array<double, 6>> stress_data;
};

/**
 * @brief Fault output receiver for dynamic rupture
 * 
 * Records fault quantities at specified fault points:
 * - Slip and slip rate
 * - Traction components
 * - State variable (for rate-state)
 * - Rupture time
 */
class FaultReceiver {
public:
    FaultReceiver();
    
    // Configuration
    void setLocation(double x, double y, double z);
    void setName(const std::string& name);
    void setSamplingRate(double dt);
    
    // Record fault data
    void record(double t, double slip, double slip_rate,
               double traction_normal, double traction_strike, double traction_dip,
               double state_variable = 0.0);
    
    // Get recorded data
    const std::vector<double>& getTimes() const { return times; }
    const std::vector<double>& getSlip() const { return slip; }
    const std::vector<double>& getSlipRate() const { return slip_rate; }
    const std::vector<double>& getTractionNormal() const { return traction_n; }
    const std::vector<double>& getTractionStrike() const { return traction_s; }
    const std::vector<double>& getTractionDip() const { return traction_d; }
    const std::vector<double>& getStateVariable() const { return state_var; }
    
    // Get rupture time (when slip rate first exceeds threshold)
    double getRuptureTime(double threshold = 1e-3) const;
    
    // Get peak slip rate
    double getPeakSlipRate() const;
    
    // Output
    void writeASCII(const std::string& filename) const;
    
    std::string name;
    
private:
    double loc_x, loc_y, loc_z;
    double sampling_dt;
    double last_record_time;
    
    std::vector<double> times;
    std::vector<double> slip;
    std::vector<double> slip_rate;
    std::vector<double> traction_n;
    std::vector<double> traction_s;
    std::vector<double> traction_d;
    std::vector<double> state_var;
};

/**
 * @brief Manager for multiple sources and receivers
 */
class SourceReceiverManager {
public:
    SourceReceiverManager();
    
    // Add sources
    void addPointSource(std::unique_ptr<PointSource> source);
    void addKinematicSource(std::unique_ptr<KinematicSource> source);
    void addSingleForceSource(std::unique_ptr<SingleForceSource> source);
    
    // Add receivers
    void addReceiver(std::unique_ptr<SeismicReceiver> receiver);
    void addFaultReceiver(std::unique_ptr<FaultReceiver> receiver);
    
    // Add receivers on a line
    void addReceiverLine(double x1, double y1, double z1,
                        double x2, double y2, double z2,
                        int num_receivers, const std::string& prefix);
    
    // Add receivers on a grid
    void addReceiverGrid(double x_min, double x_max, double y_min, double y_max,
                        double z, int nx, int ny, const std::string& prefix);
    
    // Evaluate total source term at a point
    void getSourceTerm(double x, double y, double z, double t, double* f) const;
    
    // Record at all receivers
    void recordAll(double t, 
                  std::function<void(double, double, double, double*)> get_solution);
    
    // Output all receivers
    void writeAllReceivers(const std::string& output_dir) const;
    
    // Get sources
    size_t getNumPointSources() const { return point_sources.size(); }
    size_t getNumKinematicSources() const { return kinematic_sources.size(); }
    
    // Get receivers
    size_t getNumReceivers() const { return receivers.size(); }
    size_t getNumFaultReceivers() const { return fault_receivers.size(); }
    const std::vector<std::unique_ptr<SeismicReceiver>>& getReceivers() const { return receivers; }
    
    // Time range
    double getSourceStartTime() const;
    double getSourceEndTime() const;
    
private:
    std::vector<std::unique_ptr<PointSource>> point_sources;
    std::vector<std::unique_ptr<KinematicSource>> kinematic_sources;
    std::vector<std::unique_ptr<SingleForceSource>> force_sources;
    std::vector<std::unique_ptr<SeismicReceiver>> receivers;
    std::vector<std::unique_ptr<FaultReceiver>> fault_receivers;
};

/**
 * @brief Configuration helper for sources
 */
struct SourceConfig {
    // Point source parameters
    bool use_point_source = false;
    std::string source_time_function = "gaussian";
    double source_duration = 1.0;
    double source_onset = 0.0;
    double source_x = 0.0, source_y = 0.0, source_z = 0.0;
    double strike = 0.0, dip = 90.0, rake = 0.0;
    double scalar_moment = 1e18;
    
    // Kinematic source parameters
    bool use_kinematic_source = false;
    std::string kinematic_file = "";
    double rupture_velocity = 3000.0;
    double rise_time = 1.0;
    
    // Parse from config map
    void parseConfig(const std::map<std::string, std::string>& config);
    
    // Create source from config
    std::unique_ptr<PointSource> createPointSource() const;
    std::unique_ptr<KinematicSource> createKinematicSource() const;
};

/**
 * @brief Configuration helper for receivers
 */
struct ReceiverConfig {
    std::string output_format = "ASCII";
    double sampling_rate = 100.0;
    
    // Receiver list file
    std::string receiver_file = "";
    
    // Inline receiver specification
    std::vector<std::array<double, 3>> receiver_locations;
    std::vector<std::string> receiver_names;
    
    // Parse from config map
    void parseConfig(const std::map<std::string, std::string>& config);
    
    // Create receivers from config
    std::vector<std::unique_ptr<SeismicReceiver>> createReceivers() const;
};

} // namespace FSRM

#endif // SEISMIC_SOURCE_HPP
