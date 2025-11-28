#ifndef UNIT_SYSTEM_HPP
#define UNIT_SYSTEM_HPP

#include <string>
#include <map>
#include <vector>
#include <memory>
#include <stdexcept>
#include <cmath>

namespace FSRM {

/**
 * @brief Unit dimension in terms of Length, Mass, Time (L M T)
 * 
 * Every physical quantity can be expressed as a combination of
 * fundamental dimensions: Length^a * Mass^b * Time^c
 */
struct Dimension {
    double L;  // Length exponent
    double M;  // Mass exponent
    double T;  // Time exponent
    
    Dimension(double length = 0, double mass = 0, double time = 0)
        : L(length), M(mass), T(time) {}
    
    bool operator==(const Dimension& other) const {
        return (std::abs(L - other.L) < 1e-10 &&
                std::abs(M - other.M) < 1e-10 &&
                std::abs(T - other.T) < 1e-10);
    }
    
    bool operator!=(const Dimension& other) const {
        return !(*this == other);
    }
    
    // Get human-readable dimension string
    std::string toString() const;
};

/**
 * @brief Unit definition with conversion factor to base SI units
 * 
 * Base SI units are:
 * - Length: meter (m)
 * - Mass: kilogram (kg)
 * - Time: second (s)
 */
struct Unit {
    std::string name;           // Full name (e.g., "meter")
    std::string symbol;         // Short symbol (e.g., "m")
    Dimension dimension;        // Dimensional formula
    double to_base;            // Conversion factor to base SI units
    double offset;             // Offset for affine conversions (e.g., temperature)
    std::string category;      // Category for organization
    std::vector<std::string> aliases;  // Alternative names/symbols
    
    Unit() : to_base(1.0), offset(0.0) {}
    
    Unit(const std::string& n, const std::string& s, 
         const Dimension& d, double factor, const std::string& cat = "")
        : name(n), symbol(s), dimension(d), to_base(factor), offset(0.0), category(cat) {}
    
    // Convert value from this unit to base SI
    double convertToBase(double value) const {
        return (value + offset) * to_base;
    }
    
    // Convert value from base SI to this unit
    double convertFromBase(double value) const {
        return value / to_base - offset;
    }
};

/**
 * @brief Comprehensive unit system with database and conversion utilities
 * 
 * This class provides:
 * - Database of common units across all physical quantities
 * - Conversion between any compatible units
 * - Dimensional analysis and validation
 * - Parsing of unit strings (e.g., "100 psi", "50 mD", "5 cP")
 */
class UnitSystem {
public:
    UnitSystem();
    ~UnitSystem() = default;
    
    // =========================================================================
    // Database Access
    // =========================================================================
    
    /**
     * @brief Get unit by name or symbol
     * @param name_or_symbol Unit name or symbol (case-insensitive)
     * @return Pointer to Unit, or nullptr if not found
     */
    const Unit* getUnit(const std::string& name_or_symbol) const;
    
    /**
     * @brief Check if unit exists in database
     */
    bool hasUnit(const std::string& name_or_symbol) const;
    
    /**
     * @brief Get all units in a category
     * @param category Category name (e.g., "length", "pressure")
     * @return Vector of unit pointers
     */
    std::vector<const Unit*> getUnitsInCategory(const std::string& category) const;
    
    /**
     * @brief Get all available categories
     */
    std::vector<std::string> getCategories() const;
    
    /**
     * @brief Get dimension for a unit
     */
    Dimension getDimension(const std::string& unit_name) const;
    
    // =========================================================================
    // Conversion Functions
    // =========================================================================
    
    /**
     * @brief Convert value between two units
     * @param value Input value
     * @param from_unit Source unit
     * @param to_unit Destination unit
     * @return Converted value
     * @throws std::runtime_error if units are incompatible
     */
    double convert(double value, const std::string& from_unit, 
                   const std::string& to_unit) const;
    
    /**
     * @brief Convert value to base SI units (m, kg, s)
     * @param value Input value
     * @param from_unit Source unit
     * @return Value in base SI units
     */
    double toBase(double value, const std::string& from_unit) const;
    
    /**
     * @brief Convert value from base SI units to specified unit
     * @param value Value in base SI units
     * @param to_unit Destination unit
     * @return Converted value
     */
    double fromBase(double value, const std::string& to_unit) const;
    
    /**
     * @brief Convert array of values
     */
    std::vector<double> convert(const std::vector<double>& values,
                               const std::string& from_unit,
                               const std::string& to_unit) const;
    
    // =========================================================================
    // Parsing Functions
    // =========================================================================
    
    /**
     * @brief Parse value with unit string (e.g., "100 psi", "50.5 mD")
     * @param value_with_unit String containing value and unit
     * @param[out] value Parsed value in base SI units
     * @param[out] unit Unit name
     * @return true if parsing successful
     */
    bool parseValueWithUnit(const std::string& value_with_unit,
                           double& value, std::string& unit) const;
    
    /**
     * @brief Parse value with unit and convert to base SI
     * @param value_with_unit String containing value and unit
     * @return Value in base SI units
     */
    double parseAndConvertToBase(const std::string& value_with_unit) const;
    
    // =========================================================================
    // Dimensional Analysis
    // =========================================================================
    
    /**
     * @brief Check if two units have compatible dimensions
     */
    bool areCompatible(const std::string& unit1, const std::string& unit2) const;
    
    /**
     * @brief Get base SI unit for a given dimension
     */
    std::string getBaseUnit(const Dimension& dim) const;
    
    /**
     * @brief Get suggested unit for display (commonly used in field)
     */
    std::string getSuggestedDisplayUnit(const std::string& quantity) const;
    
    // =========================================================================
    // Custom Unit Registration
    // =========================================================================
    
    /**
     * @brief Add custom unit to database
     */
    void addUnit(const Unit& unit);
    
    /**
     * @brief Add unit alias
     */
    void addAlias(const std::string& unit_name, const std::string& alias);
    
    // =========================================================================
    // Utility Functions
    // =========================================================================
    
    /**
     * @brief Format value with unit for display
     */
    std::string formatValue(double value, const std::string& unit, 
                           int precision = 6) const;
    
    /**
     * @brief Print unit database to stream (for documentation)
     */
    void printDatabase(std::ostream& os) const;
    
    /**
     * @brief Generate markdown documentation of all units
     */
    std::string generateDocumentation() const;
    
private:
    // Unit database: maps name/symbol -> Unit
    std::map<std::string, Unit> units_;
    
    // Category index: category -> list of unit names
    std::map<std::string, std::vector<std::string>> categories_;
    
    // Suggested display units for common quantities
    std::map<std::string, std::string> display_units_;
    
    // Initialize the comprehensive unit database
    void initializeDatabase();
    
    // Add units for each category
    void addLengthUnits();
    void addMassUnits();
    void addTimeUnits();
    void addPressureUnits();
    void addPermeabilityUnits();
    void addViscosityUnits();
    void addDensityUnits();
    void addVelocityUnits();
    void addAccelerationUnits();
    void addForceUnits();
    void addEnergyUnits();
    void addPowerUnits();
    void addTemperatureUnits();
    void addAngleUnits();
    void addAreaUnits();
    void addVolumeUnits();
    void addVolumetricRateUnits();
    void addMassRateUnits();
    void addThermalUnits();
    void addStressUnits();
    void addModulusUnits();
    void addFractureToughnessUnits();
    void addCompressibilityUnits();
    void addProductivityUnits();
    void addTransmissibilityUnits();
    
    // Helper to add unit with all variations
    void registerUnit(const Unit& unit);
    
    // String utilities
    std::string toLowerCase(const std::string& str) const;
    std::string trim(const std::string& str) const;
    std::vector<std::string> split(const std::string& str, char delim) const;
};

/**
 * @brief Global unit system instance (singleton pattern)
 */
class UnitSystemManager {
public:
    static UnitSystem& getInstance() {
        static UnitSystem instance;
        return instance;
    }
    
private:
    UnitSystemManager() = default;
};

// Convenience functions for quick access
inline double convertUnits(double value, const std::string& from, const std::string& to) {
    return UnitSystemManager::getInstance().convert(value, from, to);
}

inline double toSI(double value, const std::string& unit) {
    return UnitSystemManager::getInstance().toBase(value, unit);
}

inline double fromSI(double value, const std::string& unit) {
    return UnitSystemManager::getInstance().fromBase(value, unit);
}

} // namespace FSRM

#endif // UNIT_SYSTEM_HPP
