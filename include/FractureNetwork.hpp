#ifndef FRACTURE_NETWORK_HPP
#define FRACTURE_NETWORK_HPP

/**
 * @file FractureNetwork.hpp
 * @brief Advanced discrete fracture network (DFN) algorithms
 * 
 * Comprehensive fracture network modeling capabilities:
 * - Stochastic DFN generation with multiple distribution types
 * - Fracture set characterization with orientation distributions
 * - Network connectivity analysis using graph theory
 * - Percolation threshold calculations
 * - Fracture-matrix upscaling (MINC, EDFM)
 * - Complex fracture interaction models
 * - Fracture clustering and spatial analysis
 */

#include "FSRM.hpp"
#include <vector>
#include <memory>
#include <string>
#include <map>
#include <set>
#include <array>
#include <functional>
#include <random>
#include <cmath>

namespace FSRM {

// Forward declarations
class FractureNetwork;
class FractureSet;
class FractureConnectivityGraph;

/**
 * @brief Single discrete fracture representation
 */
struct DiscreteFracture {
    int id;                              ///< Unique fracture identifier
    int set_id;                          ///< Parent fracture set ID
    
    // Geometry (3D planar disc or polygon)
    std::array<double, 3> center;        ///< Center point (x, y, z)
    std::array<double, 3> normal;        ///< Unit normal vector
    std::array<double, 3> strike_dir;    ///< Strike direction vector
    std::array<double, 3> dip_dir;       ///< Dip direction vector
    double radius;                       ///< Equivalent radius (for disc)
    double length;                       ///< Length along strike
    double height;                       ///< Height along dip
    
    // Orientation (degrees)
    double strike;                       ///< Strike angle (0-360°)
    double dip;                          ///< Dip angle (0-90°)
    
    // Hydraulic properties
    double aperture;                     ///< Hydraulic aperture (m)
    double permeability;                 ///< Fracture permeability (m²)
    double transmissivity;               ///< T = k * b (m³)
    double storativity;                  ///< Fracture storativity
    
    // Mechanical properties
    double normal_stiffness;             ///< Normal stiffness (Pa/m)
    double shear_stiffness;              ///< Shear stiffness (Pa/m)
    double friction_angle;               ///< Friction angle (degrees)
    double cohesion;                     ///< Cohesion (Pa)
    double dilation_angle;               ///< Dilation angle (degrees)
    
    // State
    bool is_conductive;                  ///< Whether fracture conducts flow
    double current_aperture;             ///< Stress-dependent aperture
    double cumulative_slip;              ///< Accumulated shear slip
    
    // Connectivity
    std::vector<int> connected_fractures; ///< IDs of intersecting fractures
    std::vector<int> matrix_cells;        ///< IDs of matrix cells in contact
    
    DiscreteFracture();
    
    // Compute derived properties
    void computeVectors();               ///< Compute strike/dip vectors from angles
    void computePermeability();          ///< k = b²/12 (cubic law)
    double computeArea() const;          ///< Fracture area
    
    // Stress-dependent aperture
    double computeAperture(double normal_stress, double reference_aperture,
                           double reference_stress) const;
    
    // Check if point lies on fracture plane
    bool containsPoint(const std::array<double, 3>& point, double tolerance = 1e-6) const;
    
    // Get vertices of fracture polygon
    std::vector<std::array<double, 3>> getVertices(int n_sides = 8) const;
};

/**
 * @brief Orientation distribution types
 */
enum class OrientationDistribution {
    UNIFORM,                 ///< Uniform random orientation
    FISHER,                  ///< Fisher (von Mises-Fisher) distribution
    BINGHAM,                 ///< Bingham distribution (girdle/cluster)
    KENT,                    ///< Kent (Fisher-Bingham) distribution
    BOOTSTRAP                ///< Bootstrap from measured data
};

/**
 * @brief Size distribution types
 */
enum class SizeDistribution {
    CONSTANT,                ///< Fixed size
    UNIFORM,                 ///< Uniform distribution
    NORMAL,                  ///< Normal/Gaussian distribution
    LOGNORMAL,               ///< Log-normal distribution
    POWER_LAW,               ///< Power-law (scale-free)
    EXPONENTIAL,             ///< Exponential distribution
    PARETO                   ///< Pareto distribution
};

/**
 * @brief Spatial distribution types
 */
enum class SpatialDistribution {
    UNIFORM_POISSON,         ///< Uniform Poisson point process
    CLUSTERED_POISSON,       ///< Non-homogeneous Poisson (clustered)
    FRACTAL,                 ///< Fractal/box-counting distribution
    PARENT_DAUGHTER,         ///< Parent-daughter clustering
    LEVY_FLIGHT,             ///< Lévy flight spatial distribution
    STRATIGRAPHIC            ///< Layer-controlled distribution
};

/**
 * @brief Fracture intensity measures
 */
struct FractureIntensity {
    double P10;              ///< Fracture count per unit length (1/m)
    double P20;              ///< Fracture count per unit area (1/m²)
    double P21;              ///< Fracture trace length per unit area (m/m²)
    double P30;              ///< Fracture count per unit volume (1/m³)
    double P32;              ///< Fracture area per unit volume (m²/m³)
    double P33;              ///< Fracture volume per unit volume (dimensionless)
    
    FractureIntensity();
    void computeFromP32(double p32, double mean_radius);
    void computeFromP21(double p21, double mean_trace_length);
};

/**
 * @brief Parameters for Fisher distribution
 */
struct FisherParameters {
    double mean_strike;      ///< Mean strike angle (degrees)
    double mean_dip;         ///< Mean dip angle (degrees)
    double kappa;            ///< Concentration parameter (0 = uniform, inf = deterministic)
    
    FisherParameters(double strike = 0, double dip = 90, double k = 10)
        : mean_strike(strike), mean_dip(dip), kappa(k) {}
};

/**
 * @brief Parameters for Bingham distribution
 */
struct BinghamParameters {
    std::array<double, 3> principal_axes[3];  ///< Principal axes
    double kappa1;                             ///< Concentration along axis 1
    double kappa2;                             ///< Concentration along axis 2
    // kappa3 is constrained: kappa1 + kappa2 + kappa3 = 0
    
    BinghamParameters();
};

/**
 * @brief Fracture set definition
 * 
 * A fracture set represents a statistically homogeneous group
 * of fractures with shared orientation, size, and property distributions.
 */
class FractureSet {
public:
    FractureSet(int id, const std::string& name);
    
    // Set identification
    int getId() const { return id_; }
    const std::string& getName() const { return name_; }
    
    // Orientation distribution
    void setOrientationDistribution(OrientationDistribution type);
    void setFisherParameters(const FisherParameters& params);
    void setBinghamParameters(const BinghamParameters& params);
    void setBootstrapData(const std::vector<std::pair<double, double>>& strike_dip_pairs);
    
    // Size distribution
    void setSizeDistribution(SizeDistribution type);
    void setSizeParameters(double param1, double param2 = 0.0, double param3 = 0.0);
    void setMinMaxSize(double min_size, double max_size);
    
    // Spatial distribution
    void setSpatialDistribution(SpatialDistribution type);
    void setIntensity(double p32);  ///< Set P32 intensity
    void setClusteringParameters(double cluster_radius, double cluster_intensity);
    
    // Aperture distribution
    void setApertureDistribution(SizeDistribution type);
    void setApertureParameters(double param1, double param2 = 0.0);
    void setApertureCorrelation(double exponent);  ///< b ~ L^exponent
    
    // Mechanical properties
    void setMechanicalProperties(double Kn, double Ks, double phi, double c);
    
    // Termination rules
    enum class TerminationType {
        NONE,               ///< No termination
        OLDER_SET,          ///< Terminate against older sets
        ALL_INTERSECTING,   ///< Terminate at all intersections
        PROBABILISTIC       ///< Probabilistic termination
    };
    void setTerminationRule(TerminationType type, double probability = 1.0);
    void addTerminatingSet(int set_id);
    
    // Generate fractures
    std::vector<DiscreteFracture> generateFractures(
        const std::array<double, 3>& domain_min,
        const std::array<double, 3>& domain_max,
        std::mt19937& rng) const;
    
    // Sample single orientation
    std::pair<double, double> sampleOrientation(std::mt19937& rng) const;
    
    // Sample single size
    double sampleSize(std::mt19937& rng) const;
    
    // Sample aperture
    double sampleAperture(double fracture_size, std::mt19937& rng) const;
    
private:
    int id_;
    std::string name_;
    
    // Orientation distribution
    OrientationDistribution orientation_type_;
    FisherParameters fisher_params_;
    BinghamParameters bingham_params_;
    std::vector<std::pair<double, double>> bootstrap_data_;
    
    // Size distribution
    SizeDistribution size_type_;
    double size_param1_, size_param2_, size_param3_;
    double min_size_, max_size_;
    
    // Spatial distribution
    SpatialDistribution spatial_type_;
    double p32_;
    double cluster_radius_;
    double cluster_intensity_;
    
    // Aperture
    SizeDistribution aperture_type_;
    double aperture_param1_, aperture_param2_;
    double aperture_correlation_exp_;
    
    // Mechanical properties
    double normal_stiffness_;
    double shear_stiffness_;
    double friction_angle_;
    double cohesion_;
    
    // Termination
    TerminationType termination_type_;
    double termination_probability_;
    std::vector<int> terminating_sets_;
    
    // Helper functions
    std::pair<double, double> sampleFisherOrientation(std::mt19937& rng) const;
    std::pair<double, double> sampleBinghamOrientation(std::mt19937& rng) const;
    double samplePowerLaw(double x_min, double x_max, double alpha, std::mt19937& rng) const;
    int computeFractureCount(double volume) const;
};

/**
 * @brief Fracture intersection representation
 */
struct FractureIntersection {
    int fracture1_id;
    int fracture2_id;
    std::array<double, 3> point1;        ///< Start of intersection line
    std::array<double, 3> point2;        ///< End of intersection line
    double length;                        ///< Intersection length
    double transmissivity;                ///< Intersection transmissivity
    
    FractureIntersection();
};

/**
 * @brief Graph-based connectivity analysis
 */
class FractureConnectivityGraph {
public:
    FractureConnectivityGraph();
    
    // Build graph from fractures
    void buildFromFractures(const std::vector<DiscreteFracture>& fractures);
    
    // Add nodes and edges
    void addNode(int fracture_id);
    void addEdge(int frac1_id, int frac2_id, double weight = 1.0);
    
    // Connectivity queries
    bool isConnected(int frac1_id, int frac2_id) const;
    std::vector<int> getConnectedComponent(int fracture_id) const;
    std::vector<std::vector<int>> getAllComponents() const;
    int getLargestComponentSize() const;
    
    // Path finding
    std::vector<int> findShortestPath(int from_id, int to_id) const;
    double getPathLength(int from_id, int to_id) const;
    std::vector<std::vector<int>> findAllPaths(int from_id, int to_id, int max_depth = 100) const;
    
    // Network properties
    double getConnectivity() const;      ///< Fraction in largest component
    double getAverageClusteringCoefficient() const;
    double getAveragePathLength() const;
    int getNumComponents() const;
    
    // Percolation analysis
    struct PercolationResult {
        bool percolates_x;
        bool percolates_y;
        bool percolates_z;
        double percolation_probability;
        std::vector<int> backbone_fractures;  ///< Fractures in percolating cluster
    };
    
    PercolationResult analyzePercolation(
        const std::vector<DiscreteFracture>& fractures,
        const std::array<double, 3>& domain_min,
        const std::array<double, 3>& domain_max) const;
    
    // Backbone extraction
    std::vector<int> extractBackbone(int inlet_id, int outlet_id) const;
    
private:
    std::map<int, std::set<std::pair<int, double>>> adjacency_;  ///< id -> {(neighbor_id, weight)}
    std::map<int, int> component_id_;    ///< Fracture ID to component ID
    int num_components_;
    
    void computeComponents();
    void dfs(int node, int component, std::set<int>& visited);
};

/**
 * @brief Main fracture network class
 */
class FractureNetwork {
public:
    FractureNetwork();
    ~FractureNetwork() = default;
    
    // Domain setup
    void setDomain(const std::array<double, 3>& min_corner,
                   const std::array<double, 3>& max_corner);
    void setDomainFromMesh(DM dm);
    
    // Fracture set management
    void addFractureSet(const FractureSet& set);
    FractureSet* getFractureSet(int set_id);
    const std::vector<FractureSet>& getFractureSets() const { return sets_; }
    
    // DFN generation
    void setRandomSeed(unsigned int seed);
    void generate();
    void generateWithIntensityControl(double target_p32, double tolerance = 0.05);
    
    // Import/export
    void importFromFile(const std::string& filename, const std::string& format = "auto");
    void exportToFile(const std::string& filename, const std::string& format = "vtk") const;
    void exportToGmsh(const std::string& filename) const;
    
    // Access fractures
    const std::vector<DiscreteFracture>& getFractures() const { return fractures_; }
    std::vector<DiscreteFracture>& getFractures() { return fractures_; }
    DiscreteFracture* getFracture(int id);
    const DiscreteFracture* getFracture(int id) const;
    size_t getNumFractures() const { return fractures_.size(); }
    
    // Intersection detection
    void computeIntersections();
    const std::vector<FractureIntersection>& getIntersections() const { return intersections_; }
    
    // Connectivity analysis
    void buildConnectivityGraph();
    const FractureConnectivityGraph& getConnectivityGraph() const { return connectivity_; }
    FractureConnectivityGraph& getConnectivityGraph() { return connectivity_; }
    
    // Network properties
    FractureIntensity computeIntensity() const;
    double computePercolationThreshold(int n_realizations = 100);
    
    // Upscaling methods
    
    /**
     * @brief MINC (Multiple INteracting Continua) upscaling
     * 
     * Creates nested continua for dual-porosity/permeability modeling.
     * 
     * @param n_continua Number of interacting continua
     * @param volume_fractions Volume fraction for each continuum
     * @return Shape factors and transmissibilities
     */
    struct MINCResult {
        std::vector<double> volume_fractions;
        std::vector<double> shape_factors;
        std::vector<std::vector<double>> transmissibilities;  ///< Between continua
    };
    MINCResult computeMINC(int n_continua) const;
    
    /**
     * @brief EDFM (Embedded Discrete Fracture Model) upscaling
     * 
     * Computes fracture-matrix and fracture-fracture transmissibilities
     * for embedding DFN in continuum grid.
     * 
     * @param dm PETSc DM for matrix grid
     * @return Connection list with transmissibilities
     */
    struct EDFMConnection {
        enum Type { FRACTURE_MATRIX, FRACTURE_FRACTURE };
        Type type;
        int id1, id2;                     ///< Fracture/cell IDs
        double transmissibility;          ///< Connection transmissibility
        double distance;                  ///< Average distance
        double area;                      ///< Connection area
    };
    std::vector<EDFMConnection> computeEDFM(DM dm) const;
    
    /**
     * @brief Oda's crack tensor upscaling
     * 
     * Computes equivalent permeability tensor from DFN.
     * 
     * @return 3x3 permeability tensor (m²)
     */
    std::array<std::array<double, 3>, 3> computeOdaTensor() const;
    
    /**
     * @brief Snow's model for parallel plate permeability
     * 
     * @return Directional permeabilities (kx, ky, kz)
     */
    std::array<double, 3> computeSnowPermeability() const;
    
    // Stress-dependent properties
    void updateApertures(const std::vector<double>& stress_field);
    void updatePermeabilities();
    
    // Visualization
    void writeVTK(const std::string& filename) const;
    void writeFractureTraces(const std::string& filename, double z_level) const;
    
private:
    // Domain
    std::array<double, 3> domain_min_;
    std::array<double, 3> domain_max_;
    double domain_volume_;
    
    // Fracture sets and fractures
    std::vector<FractureSet> sets_;
    std::vector<DiscreteFracture> fractures_;
    std::vector<FractureIntersection> intersections_;
    
    // Connectivity
    FractureConnectivityGraph connectivity_;
    
    // Random number generator
    std::mt19937 rng_;
    
    // Helper functions
    bool detectIntersection(const DiscreteFracture& f1, const DiscreteFracture& f2,
                           FractureIntersection& intersection) const;
    void applyTerminationRules();
    void clipToDomain();
};

/**
 * @brief Fracture mechanics propagation model
 * 
 * Extends DFN with propagation capabilities using LEFM.
 */
class FracturePropagationModel {
public:
    FracturePropagationModel();
    
    // Material properties
    void setRockProperties(double E, double nu, double KIC);
    void setFluidProperties(double viscosity, double compressibility);
    
    // Stress state
    void setFarFieldStress(double Sxx, double Syy, double Szz,
                           double Sxy, double Sxz, double Syz);
    
    // Propagation analysis
    struct TipState {
        int fracture_id;
        std::array<double, 3> tip_position;
        std::array<double, 3> propagation_direction;
        double K_I, K_II, K_III;         ///< Stress intensity factors
        double energy_release_rate;       ///< G = (K_I² + K_II²)/E' + K_III²/(2μ)
        bool will_propagate;
    };
    
    std::vector<TipState> computeTipStates(const FractureNetwork& network) const;
    
    // Propagate fractures
    void propagateNetwork(FractureNetwork& network, double dt);
    
    // Interaction effects
    void computeStressShadowing(FractureNetwork& network) const;
    
    // Coalescence detection
    std::vector<std::pair<int, int>> detectCoalescence(
        const FractureNetwork& network, double threshold_distance) const;
    
private:
    double youngs_modulus_;
    double poisson_ratio_;
    double fracture_toughness_;
    double fluid_viscosity_;
    double fluid_compressibility_;
    std::array<std::array<double, 3>, 3> far_field_stress_;
    
    double computeKI(const DiscreteFracture& frac, double net_pressure) const;
    std::array<double, 3> computePropagationDirection(const TipState& state) const;
};

/**
 * @brief Fracture flow simulator (standalone, doesn't require full reservoir)
 */
class FractureFlowSimulator {
public:
    FractureFlowSimulator();
    
    // Setup
    void setNetwork(FractureNetwork* network);
    void setFluidProperties(double viscosity, double compressibility);
    
    // Boundary conditions
    void setInletPressure(int fracture_id, double pressure);
    void setOutletPressure(int fracture_id, double pressure);
    void setInletFlowRate(int fracture_id, double rate);
    
    // Solve steady-state flow
    struct FlowSolution {
        std::vector<double> pressures;       ///< Pressure at each fracture
        std::vector<double> flow_rates;      ///< Flow rate through each fracture
        double total_inflow;
        double total_outflow;
        double effective_transmissivity;
    };
    
    FlowSolution solveSteadyState();
    
    // Compute equivalent permeability
    std::array<double, 3> computeEquivalentPermeability();
    
private:
    FractureNetwork* network_;
    double viscosity_;
    double compressibility_;
    
    std::map<int, double> inlet_pressures_;
    std::map<int, double> outlet_pressures_;
    std::map<int, double> inlet_rates_;
};

/**
 * @brief Statistical analysis utilities for DFN
 */
class DFNStatistics {
public:
    // Orientation analysis
    static std::vector<FisherParameters> clusterOrientations(
        const std::vector<DiscreteFracture>& fractures,
        int n_clusters);
    
    static double computeFisherKappa(
        const std::vector<std::pair<double, double>>& orientations);
    
    // Size distribution fitting
    static std::pair<SizeDistribution, std::vector<double>> fitSizeDistribution(
        const std::vector<double>& sizes);
    
    // Spatial analysis
    static double computeFractalDimension(
        const std::vector<DiscreteFracture>& fractures,
        const std::array<double, 3>& domain_size);
    
    static std::vector<double> computeRipleyK(
        const std::vector<DiscreteFracture>& fractures,
        const std::vector<double>& distances);
    
    // Intensity computation
    static FractureIntensity computeIntensityFromScanline(
        const std::vector<DiscreteFracture>& fractures,
        const std::array<double, 3>& line_start,
        const std::array<double, 3>& line_end);
    
    static FractureIntensity computeIntensityFromWindow(
        const std::vector<DiscreteFracture>& fractures,
        const std::array<double, 3>& window_center,
        const std::array<double, 3>& window_size);
};

/**
 * @brief Parse DFN configuration from file
 */
void parseDFNConfig(const std::string& filename, FractureNetwork& network);

} // namespace FSRM

#endif // FRACTURE_NETWORK_HPP
