#ifndef AQUIFER_MODEL_HPP
#define AQUIFER_MODEL_HPP

/**
 * @file AquiferModel.hpp
 * @brief ECLIPSE-compatible aquifer models
 * 
 * Implements analytical and numerical aquifer models used in ECLIPSE:
 * - Carter-Tracy infinite-acting aquifer
 * - Fetkovich aquifer (pseudo-steady state)
 * - Numerical aquifer (gridded)
 * - Constant pressure/flux aquifers
 * 
 * ECLIPSE Keywords: AQUANCON, AQUCT, AQUFETP, AQUNUM
 */

#include <vector>
#include <memory>
#include <string>
#include <map>
#include <cmath>
#include <functional>

namespace FSRM {

/**
 * @brief Aquifer type enumeration (ECLIPSE compatible)
 */
enum class AquiferType {
    CARTER_TRACY,       ///< Analytical infinite-acting (AQUCT)
    FETKOVICH,          ///< Pseudo-steady state (AQUFETP)
    NUMERICAL,          ///< Gridded aquifer (AQUNUM)
    CONSTANT_PRESSURE,  ///< Constant pressure boundary
    CONSTANT_FLUX       ///< Constant influx rate
};

/**
 * @brief Aquifer geometry for Carter-Tracy model
 */
enum class AquiferGeometry {
    RADIAL,             ///< Radial flow (most common)
    LINEAR,             ///< Linear flow
    BOTTOM,             ///< Bottom water drive
    EDGE                ///< Edge water drive
};

/**
 * @brief Connection face between aquifer and reservoir
 */
enum class AquiferFace {
    I_MINUS,    ///< -X face
    I_PLUS,     ///< +X face
    J_MINUS,    ///< -Y face
    J_PLUS,     ///< +Y face
    K_MINUS,    ///< -Z face (bottom)
    K_PLUS      ///< +Z face (top)
};

/**
 * @brief Aquifer-reservoir connection specification (AQUANCON)
 */
struct AquiferConnection {
    int aquifer_id;             ///< Aquifer ID
    int i1, i2, j1, j2, k1, k2; ///< Grid block range
    AquiferFace face;           ///< Connection face
    double trans_mult;          ///< Transmissibility multiplier
    double area;                ///< Connection area (m²)
    bool allow_crossflow;       ///< Allow crossflow to aquifer
    
    AquiferConnection() : aquifer_id(1), i1(1), i2(1), j1(1), j2(1), k1(1), k2(1),
                          face(AquiferFace::I_MINUS), trans_mult(1.0), 
                          area(0.0), allow_crossflow(true) {}
};

/**
 * @brief Dimensionless time and pressure functions for Carter-Tracy
 */
class DimensionlessInflux {
public:
    /**
     * @brief Van Everdingen-Hurst dimensionless pressure function
     * @param tD Dimensionless time
     * @param reD Dimensionless outer radius (re/rw)
     * @return Dimensionless cumulative influx
     */
    static double pD_infinite(double tD);
    
    /**
     * @brief Finite aquifer dimensionless pressure
     * @param tD Dimensionless time
     * @param reD Dimensionless outer radius
     * @return Dimensionless cumulative influx
     */
    static double pD_finite(double tD, double reD);
    
    /**
     * @brief Derivative dpD/dtD for superposition
     * @param tD Dimensionless time
     * @return Dimensionless pressure derivative
     */
    static double dpD_dtD(double tD);
    
    /**
     * @brief Linear aquifer dimensionless influx
     * @param tD Dimensionless time
     * @return Dimensionless influx
     */
    static double pD_linear(double tD);
    
private:
    // Polynomial coefficients for van Everdingen-Hurst approximations
    static constexpr double a1 = 1.5;
    static constexpr double a2 = 4.29881;
    static constexpr double a3 = 2.02566;
};

/**
 * @brief Base class for aquifer models
 */
class AquiferModelBase {
public:
    AquiferModelBase(int id, AquiferType type);
    virtual ~AquiferModelBase() = default;
    
    /**
     * @brief Calculate water influx rate
     * @param reservoir_pressure Current average reservoir pressure (Pa)
     * @param dt Time step (s)
     * @return Water influx rate (m³/s)
     */
    virtual double calculateInflux(double reservoir_pressure, double dt) = 0;
    
    /**
     * @brief Update aquifer state after time step
     * @param influx_rate Calculated influx rate (m³/s)
     * @param dt Time step (s)
     */
    virtual void updateState(double influx_rate, double dt) = 0;
    
    /**
     * @brief Reset aquifer to initial conditions
     */
    virtual void reset() = 0;
    
    /**
     * @brief Configure from key-value pairs
     */
    virtual void configure(const std::map<std::string, std::string>& config) = 0;
    
    // Common accessors
    int getId() const { return aquifer_id; }
    AquiferType getType() const { return aquifer_type; }
    double getInitialPressure() const { return initial_pressure; }
    double getCurrentPressure() const { return current_pressure; }
    double getCumulativeInflux() const { return cumulative_influx; }
    double getTime() const { return time; }
    
    // Add connections
    void addConnection(const AquiferConnection& conn);
    const std::vector<AquiferConnection>& getConnections() const { return connections; }
    
protected:
    int aquifer_id;
    AquiferType aquifer_type;
    double initial_pressure;        ///< Initial aquifer pressure (Pa)
    double current_pressure;        ///< Current aquifer pressure (Pa)
    double cumulative_influx;       ///< Cumulative water influx (m³)
    double time;                    ///< Current time (s)
    std::vector<AquiferConnection> connections;
};

/**
 * @brief Carter-Tracy aquifer model (AQUCT)
 * 
 * Infinite-acting or finite radial aquifer using superposition.
 * Most accurate for unsteady-state aquifer response.
 * 
 * ECLIPSE keyword: AQUCT
 */
class CarterTracyAquifer : public AquiferModelBase {
public:
    CarterTracyAquifer(int id = 1);
    
    double calculateInflux(double reservoir_pressure, double dt) override;
    void updateState(double influx_rate, double dt) override;
    void reset() override;
    void configure(const std::map<std::string, std::string>& config) override;
    
    // Setters for aquifer properties
    void setAquiferProperties(double perm_md, double porosity, 
                              double compressibility, double viscosity,
                              double thickness, double encroachment_angle);
    
    void setAquiferRadius(double inner_radius, double outer_radius);
    void setInitialPressure(double P0);
    void setGeometry(AquiferGeometry geom);
    
    // Getters
    double getInfluxConstant() const { return influx_constant; }
    double getTimeConstant() const { return time_constant; }
    double getDimensionlessRadius() const { return reD; }
    
private:
    // Aquifer properties
    double permeability;            ///< Aquifer permeability (m²)
    double porosity;                ///< Aquifer porosity
    double total_compressibility;   ///< Total compressibility (1/Pa)
    double water_viscosity;         ///< Water viscosity (Pa·s)
    double thickness;               ///< Aquifer thickness (m)
    double encroachment_angle;      ///< Fraction of circle (0-1)
    
    // Geometry
    AquiferGeometry geometry;
    double inner_radius;            ///< Inner (reservoir) radius (m)
    double outer_radius;            ///< Outer (aquifer) radius (m)
    double reD;                     ///< Dimensionless outer radius
    
    // Calculated constants
    double influx_constant;         ///< B (m³/Pa)
    double time_constant;           ///< (φ μ ct ro²) / k
    
    // Superposition variables
    std::vector<double> pressure_history;
    std::vector<double> time_history;
    int num_steps;
    
    // Calculate influx using superposition
    double superpositionInflux(double current_pressure);
    void computeConstants();
};

/**
 * @brief Fetkovich aquifer model (AQUFETP)
 * 
 * Pseudo-steady state aquifer, simpler than Carter-Tracy.
 * Good for depleted aquifers or strong water drive.
 * 
 * ECLIPSE keyword: AQUFETP
 */
class FetkovichAquifer : public AquiferModelBase {
public:
    FetkovichAquifer(int id = 1);
    
    double calculateInflux(double reservoir_pressure, double dt) override;
    void updateState(double influx_rate, double dt) override;
    void reset() override;
    void configure(const std::map<std::string, std::string>& config) override;
    
    // Setters
    void setProductivityIndex(double J);  ///< m³/s/Pa
    void setInitialVolume(double V0);     ///< Initial water volume (m³)
    void setCompressibility(double ct);   ///< Total compressibility (1/Pa)
    void setInitialPressure(double P0);
    
    // Getters
    double getProductivityIndex() const { return productivity_index; }
    double getInitialVolume() const { return initial_volume; }
    double getRemainingVolume() const { return remaining_volume; }
    
private:
    double productivity_index;      ///< Aquifer PI (m³/s/Pa)
    double initial_volume;          ///< Initial encroachable volume (m³)
    double remaining_volume;        ///< Current encroachable volume (m³)
    double total_compressibility;   ///< Total compressibility (1/Pa)
};

/**
 * @brief Numerical aquifer model (AQUNUM)
 * 
 * Gridded aquifer cells for complex geometries.
 * More flexible but computationally expensive.
 * 
 * ECLIPSE keyword: AQUNUM
 */
class NumericalAquifer : public AquiferModelBase {
public:
    NumericalAquifer(int id = 1);
    
    double calculateInflux(double reservoir_pressure, double dt) override;
    void updateState(double influx_rate, double dt) override;
    void reset() override;
    void configure(const std::map<std::string, std::string>& config) override;
    
    // Add aquifer cells
    struct AquiferCell {
        int i, j, k;                ///< Grid indices
        double volume;              ///< Cell volume (m³)
        double porosity;            ///< Cell porosity
        double permeability;        ///< Cell permeability (m²)
        double pressure;            ///< Cell pressure (Pa)
        double depth;               ///< Cell depth (m)
    };
    
    void addCell(const AquiferCell& cell);
    void setWaterProperties(double density, double viscosity, double compressibility);
    
    // Getters
    int getNumCells() const { return static_cast<int>(cells.size()); }
    double getTotalPoreVolume() const;
    
private:
    std::vector<AquiferCell> cells;
    double water_density;
    double water_viscosity;
    double water_compressibility;
    
    // Calculate transmissibility between cells
    double calculateTransmissibility(const AquiferCell& c1, const AquiferCell& c2,
                                    double dx, double dy, double dz) const;
};

/**
 * @brief Constant pressure aquifer (simple boundary condition)
 */
class ConstantPressureAquifer : public AquiferModelBase {
public:
    ConstantPressureAquifer(int id = 1);
    
    double calculateInflux(double reservoir_pressure, double dt) override;
    void updateState(double influx_rate, double dt) override;
    void reset() override;
    void configure(const std::map<std::string, std::string>& config) override;
    
    void setTransmissibility(double T);  ///< m³/s/Pa
    void setConstantPressure(double P);
    
private:
    double transmissibility;
    double constant_pressure;
};

/**
 * @brief Constant flux aquifer (simple boundary condition)
 */
class ConstantFluxAquifer : public AquiferModelBase {
public:
    ConstantFluxAquifer(int id = 1);
    
    double calculateInflux(double reservoir_pressure, double dt) override;
    void updateState(double influx_rate, double dt) override;
    void reset() override;
    void configure(const std::map<std::string, std::string>& config) override;
    
    void setFluxRate(double Q);  ///< m³/s
    
private:
    double flux_rate;
};

/**
 * @brief Aquifer manager handling multiple aquifers
 */
class AquiferManager {
public:
    AquiferManager();
    
    // Add aquifers
    void addCarterTracyAquifer(int id, const std::map<std::string, std::string>& config);
    void addFetkovichAquifer(int id, const std::map<std::string, std::string>& config);
    void addNumericalAquifer(int id, const std::map<std::string, std::string>& config);
    void addConstantPressureAquifer(int id, double pressure, double transmissibility);
    void addConstantFluxAquifer(int id, double flux_rate);
    
    // Add connections (AQUANCON)
    void addConnection(const AquiferConnection& conn);
    
    // Calculate total influx from all aquifers
    double calculateTotalInflux(double reservoir_pressure, double dt);
    
    // Calculate influx to specific cell
    double calculateCellInflux(int i, int j, int k, 
                               double cell_pressure, double dt);
    
    // Update all aquifers
    void updateAll(double dt);
    
    // Reset all aquifers
    void resetAll();
    
    // Accessors
    AquiferModelBase* getAquifer(int id);
    const std::vector<std::unique_ptr<AquiferModelBase>>& getAquifers() const { 
        return aquifers; 
    }
    
    // Output
    double getTotalCumulativeInflux() const;
    void printSummary() const;
    
private:
    std::vector<std::unique_ptr<AquiferModelBase>> aquifers;
    std::map<int, size_t> aquifer_index;  // ID to index mapping
    
    // Cell connections
    std::map<std::tuple<int,int,int>, std::vector<int>> cell_aquifer_map;
};

/**
 * @brief Factory function to create aquifer from ECLIPSE keywords
 */
std::unique_ptr<AquiferModelBase> createAquifer(
    AquiferType type, int id,
    const std::map<std::string, std::string>& config);

/**
 * @brief Parse aquifer type from string
 */
AquiferType parseAquiferType(const std::string& type_str);

/**
 * @brief Parse aquifer geometry from string
 */
AquiferGeometry parseAquiferGeometry(const std::string& geom_str);

/**
 * @brief Parse aquifer face from string
 */
AquiferFace parseAquiferFace(const std::string& face_str);

} // namespace FSRM

#endif // AQUIFER_MODEL_HPP
