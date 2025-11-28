#include "UnitSystem.hpp"
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cctype>

namespace FSRM {

// =============================================================================
// Dimension Implementation
// =============================================================================

std::string Dimension::toString() const {
    std::stringstream ss;
    bool first = true;
    
    if (std::abs(L) > 1e-10) {
        ss << "L";
        if (std::abs(L - 1.0) > 1e-10) ss << "^" << L;
        first = false;
    }
    
    if (std::abs(M) > 1e-10) {
        if (!first) ss << " ";
        ss << "M";
        if (std::abs(M - 1.0) > 1e-10) ss << "^" << M;
        first = false;
    }
    
    if (std::abs(T) > 1e-10) {
        if (!first) ss << " ";
        ss << "T";
        if (std::abs(T - 1.0) > 1e-10) ss << "^" << T;
    }
    
    return ss.str().empty() ? "dimensionless" : ss.str();
}

// =============================================================================
// UnitSystem Implementation
// =============================================================================

UnitSystem::UnitSystem() {
    initializeDatabase();
}

void UnitSystem::initializeDatabase() {
    // Initialize all unit categories
    addLengthUnits();
    addMassUnits();
    addTimeUnits();
    addAreaUnits();
    addVolumeUnits();
    addAngleUnits();
    addVelocityUnits();
    addAccelerationUnits();
    addForceUnits();
    addPressureUnits();
    addEnergyUnits();
    addPowerUnits();
    addDensityUnits();
    addViscosityUnits();
    addPermeabilityUnits();
    addVolumetricRateUnits();
    addMassRateUnits();
    addTemperatureUnits();
    addThermalUnits();
    addStressUnits();
    addModulusUnits();
    addFractureToughnessUnits();
    addCompressibilityUnits();
    addProductivityUnits();
    addTransmissibilityUnits();
    
    // Set suggested display units for common quantities
    display_units_["length"] = "m";
    display_units_["depth"] = "m";
    display_units_["pressure"] = "MPa";
    display_units_["permeability"] = "mD";
    display_units_["viscosity"] = "cP";
    display_units_["density"] = "kg/m3";
    display_units_["time"] = "s";
    display_units_["temperature"] = "degC";
    display_units_["porosity"] = "fraction";
    display_units_["modulus"] = "GPa";
    display_units_["stress"] = "MPa";
}

// =============================================================================
// Length Units
// =============================================================================

void UnitSystem::addLengthUnits() {
    Dimension length(1, 0, 0);
    
    // Metric
    registerUnit(Unit("meter", "m", length, 1.0, "length"));
    registerUnit(Unit("centimeter", "cm", length, 0.01, "length"));
    registerUnit(Unit("millimeter", "mm", length, 0.001, "length"));
    registerUnit(Unit("kilometer", "km", length, 1000.0, "length"));
    registerUnit(Unit("micrometer", "um", length, 1e-6, "length"));
    registerUnit(Unit("nanometer", "nm", length, 1e-9, "length"));
    
    // Imperial/US
    registerUnit(Unit("foot", "ft", length, 0.3048, "length"));
    registerUnit(Unit("inch", "in", length, 0.0254, "length"));
    registerUnit(Unit("yard", "yd", length, 0.9144, "length"));
    registerUnit(Unit("mile", "mi", length, 1609.344, "length"));
    
    // Oil field
    registerUnit(Unit("angstrom", "angstrom", length, 1e-10, "length"));
}

// =============================================================================
// Mass Units
// =============================================================================

void UnitSystem::addMassUnits() {
    Dimension mass(0, 1, 0);
    
    // Metric
    registerUnit(Unit("kilogram", "kg", mass, 1.0, "mass"));
    registerUnit(Unit("gram", "g", mass, 0.001, "mass"));
    registerUnit(Unit("milligram", "mg", mass, 1e-6, "mass"));
    registerUnit(Unit("tonne", "tonne", mass, 1000.0, "mass"));
    registerUnit(Unit("metric ton", "t", mass, 1000.0, "mass"));
    
    // Imperial/US
    registerUnit(Unit("pound mass", "lbm", mass, 0.45359237, "mass"));
    registerUnit(Unit("ounce", "oz", mass, 0.028349523125, "mass"));
    registerUnit(Unit("slug", "slug", mass, 14.5939029, "mass"));
    registerUnit(Unit("short ton", "ton", mass, 907.18474, "mass"));
}

// =============================================================================
// Time Units
// =============================================================================

void UnitSystem::addTimeUnits() {
    Dimension time(0, 0, 1);
    
    registerUnit(Unit("second", "s", time, 1.0, "time"));
    registerUnit(Unit("millisecond", "ms", time, 0.001, "time"));
    registerUnit(Unit("microsecond", "us", time, 1e-6, "time"));
    registerUnit(Unit("minute", "min", time, 60.0, "time"));
    registerUnit(Unit("hour", "hr", time, 3600.0, "time"));
    registerUnit(Unit("day", "day", time, 86400.0, "time"));
    registerUnit(Unit("week", "week", time, 604800.0, "time"));
    registerUnit(Unit("year", "year", time, 31536000.0, "time"));
}

// =============================================================================
// Area Units
// =============================================================================

void UnitSystem::addAreaUnits() {
    Dimension area(2, 0, 0);
    
    // Metric
    registerUnit(Unit("square meter", "m2", area, 1.0, "area"));
    registerUnit(Unit("square centimeter", "cm2", area, 1e-4, "area"));
    registerUnit(Unit("square millimeter", "mm2", area, 1e-6, "area"));
    registerUnit(Unit("square kilometer", "km2", area, 1e6, "area"));
    registerUnit(Unit("hectare", "ha", area, 1e4, "area"));
    
    // Imperial/US
    registerUnit(Unit("square foot", "ft2", area, 0.09290304, "area"));
    registerUnit(Unit("square inch", "in2", area, 0.00064516, "area"));
    registerUnit(Unit("square yard", "yd2", area, 0.83612736, "area"));
    registerUnit(Unit("square mile", "mi2", area, 2589988.110336, "area"));
    registerUnit(Unit("acre", "acre", area, 4046.8564224, "area"));
}

// =============================================================================
// Volume Units
// =============================================================================

void UnitSystem::addVolumeUnits() {
    Dimension volume(3, 0, 0);
    
    // Metric
    registerUnit(Unit("cubic meter", "m3", volume, 1.0, "volume"));
    registerUnit(Unit("cubic centimeter", "cm3", volume, 1e-6, "volume"));
    registerUnit(Unit("cubic millimeter", "mm3", volume, 1e-9, "volume"));
    registerUnit(Unit("liter", "L", volume, 0.001, "volume"));
    registerUnit(Unit("milliliter", "mL", volume, 1e-6, "volume"));
    
    // Imperial/US
    registerUnit(Unit("cubic foot", "ft3", volume, 0.028316846592, "volume"));
    registerUnit(Unit("cubic inch", "in3", volume, 1.6387064e-5, "volume"));
    registerUnit(Unit("gallon US", "gal", volume, 0.003785411784, "volume"));
    registerUnit(Unit("barrel", "bbl", volume, 0.158987294928, "volume"));
    
    // Oil field
    registerUnit(Unit("stock tank barrel", "stb", volume, 0.158987294928, "volume"));
    registerUnit(Unit("reservoir barrel", "rb", volume, 0.158987294928, "volume"));
    registerUnit(Unit("standard cubic meter", "sm3", volume, 1.0, "volume"));
    registerUnit(Unit("standard cubic foot", "scf", volume, 0.028316846592, "volume"));
    registerUnit(Unit("thousand cubic feet", "Mcf", volume, 28.316846592, "volume"));
    registerUnit(Unit("million cubic feet", "MMcf", volume, 28316.846592, "volume"));
}

// =============================================================================
// Angle Units
// =============================================================================

void UnitSystem::addAngleUnits() {
    Dimension angle(0, 0, 0);  // Dimensionless
    
    registerUnit(Unit("radian", "rad", angle, 1.0, "angle"));
    registerUnit(Unit("degree", "deg", angle, M_PI/180.0, "angle"));
    registerUnit(Unit("gradian", "grad", angle, M_PI/200.0, "angle"));
}

// =============================================================================
// Velocity Units
// =============================================================================

void UnitSystem::addVelocityUnits() {
    Dimension velocity(1, 0, -1);
    
    // Metric
    registerUnit(Unit("meter per second", "m/s", velocity, 1.0, "velocity"));
    registerUnit(Unit("centimeter per second", "cm/s", velocity, 0.01, "velocity"));
    registerUnit(Unit("kilometer per hour", "km/h", velocity, 1.0/3.6, "velocity"));
    
    // Imperial/US
    registerUnit(Unit("foot per second", "ft/s", velocity, 0.3048, "velocity"));
    registerUnit(Unit("foot per minute", "ft/min", velocity, 0.00508, "velocity"));
    registerUnit(Unit("mile per hour", "mph", velocity, 0.44704, "velocity"));
    
    // Oil field
    registerUnit(Unit("meter per day", "m/day", velocity, 1.0/86400.0, "velocity"));
    registerUnit(Unit("foot per day", "ft/day", velocity, 0.3048/86400.0, "velocity"));
}

// =============================================================================
// Acceleration Units
// =============================================================================

void UnitSystem::addAccelerationUnits() {
    Dimension acceleration(1, 0, -2);
    
    registerUnit(Unit("meter per second squared", "m/s2", acceleration, 1.0, "acceleration"));
    registerUnit(Unit("foot per second squared", "ft/s2", acceleration, 0.3048, "acceleration"));
    registerUnit(Unit("gal", "gal", acceleration, 0.01, "acceleration"));  // cm/s²
    registerUnit(Unit("standard gravity", "g", acceleration, 9.80665, "acceleration"));
}

// =============================================================================
// Force Units
// =============================================================================

void UnitSystem::addForceUnits() {
    Dimension force(1, 1, -2);
    
    registerUnit(Unit("newton", "N", force, 1.0, "force"));
    registerUnit(Unit("kilonewton", "kN", force, 1000.0, "force"));
    registerUnit(Unit("dyne", "dyne", force, 1e-5, "force"));
    registerUnit(Unit("pound force", "lbf", force, 4.4482216152605, "force"));
    registerUnit(Unit("kip", "kip", force, 4448.2216152605, "force"));
}

// =============================================================================
// Pressure Units
// =============================================================================

void UnitSystem::addPressureUnits() {
    Dimension pressure(−1, 1, -2);
    
    // SI
    registerUnit(Unit("pascal", "Pa", pressure, 1.0, "pressure"));
    registerUnit(Unit("kilopascal", "kPa", pressure, 1000.0, "pressure"));
    registerUnit(Unit("megapascal", "MPa", pressure, 1e6, "pressure"));
    registerUnit(Unit("gigapascal", "GPa", pressure, 1e9, "pressure"));
    
    // Other metric
    registerUnit(Unit("bar", "bar", pressure, 1e5, "pressure"));
    registerUnit(Unit("millibar", "mbar", pressure, 100.0, "pressure"));
    registerUnit(Unit("atmosphere", "atm", pressure, 101325.0, "pressure"));
    
    // Imperial/US
    registerUnit(Unit("pounds per square inch", "psi", pressure, 6894.757293168, "pressure"));
    registerUnit(Unit("pounds per square foot", "psf", pressure, 47.88025898033584, "pressure"));
    registerUnit(Unit("ksi", "ksi", pressure, 6894757.293168, "pressure"));
    
    // Other
    registerUnit(Unit("torr", "torr", pressure, 133.322368421, "pressure"));
    registerUnit(Unit("millimeter mercury", "mmHg", pressure, 133.322368421, "pressure"));
    registerUnit(Unit("inch mercury", "inHg", pressure, 3386.389, "pressure"));
}

// =============================================================================
// Energy Units
// =============================================================================

void UnitSystem::addEnergyUnits() {
    Dimension energy(2, 1, -2);
    
    // SI
    registerUnit(Unit("joule", "J", energy, 1.0, "energy"));
    registerUnit(Unit("kilojoule", "kJ", energy, 1000.0, "energy"));
    registerUnit(Unit("megajoule", "MJ", energy, 1e6, "energy"));
    
    // Other metric
    registerUnit(Unit("erg", "erg", energy, 1e-7, "energy"));
    registerUnit(Unit("calorie", "cal", energy, 4.184, "energy"));
    registerUnit(Unit("kilocalorie", "kcal", energy, 4184.0, "energy"));
    
    // Imperial/US
    registerUnit(Unit("british thermal unit", "BTU", energy, 1055.05585262, "energy"));
    registerUnit(Unit("foot pound force", "ft-lbf", energy, 1.3558179483314, "energy"));
    
    // Electrical
    registerUnit(Unit("watt hour", "Wh", energy, 3600.0, "energy"));
    registerUnit(Unit("kilowatt hour", "kWh", energy, 3.6e6, "energy"));
    registerUnit(Unit("electron volt", "eV", energy, 1.602176634e-19, "energy"));
}

// =============================================================================
// Power Units
// =============================================================================

void UnitSystem::addPowerUnits() {
    Dimension power(2, 1, -3);
    
    // SI
    registerUnit(Unit("watt", "W", power, 1.0, "power"));
    registerUnit(Unit("kilowatt", "kW", power, 1000.0, "power"));
    registerUnit(Unit("megawatt", "MW", power, 1e6, "power"));
    
    // Other
    registerUnit(Unit("horsepower", "hp", power, 745.69987158227, "power"));
    registerUnit(Unit("metric horsepower", "PS", power, 735.49875, "power"));
    registerUnit(Unit("BTU per hour", "BTU/hr", power, 0.29307107017, "power"));
}

// =============================================================================
// Density Units
// =============================================================================

void UnitSystem::addDensityUnits() {
    Dimension density(−3, 1, 0);
    
    // SI
    registerUnit(Unit("kilogram per cubic meter", "kg/m3", density, 1.0, "density"));
    registerUnit(Unit("gram per cubic centimeter", "g/cm3", density, 1000.0, "density"));
    
    // Imperial/US
    registerUnit(Unit("pound mass per cubic foot", "lbm/ft3", density, 16.018463373960142, "density"));
    registerUnit(Unit("pound mass per gallon", "lbm/gal", density, 119.82642730074315, "density"));
    
    // Oil field - API gravity is special (needs different handling)
    registerUnit(Unit("API gravity", "API", density, 1.0, "density"));  // Special conversion
}

// =============================================================================
// Viscosity Units
// =============================================================================

void UnitSystem::addViscosityUnits() {
    Dimension viscosity(−1, 1, -1);
    
    // Dynamic viscosity
    registerUnit(Unit("pascal second", "Pa-s", viscosity, 1.0, "viscosity"));
    registerUnit(Unit("centipoise", "cP", viscosity, 0.001, "viscosity"));
    registerUnit(Unit("poise", "P", viscosity, 0.1, "viscosity"));
    registerUnit(Unit("millipascal second", "mPa-s", viscosity, 0.001, "viscosity"));
    registerUnit(Unit("pound per foot second", "lb/ft-s", viscosity, 1.4881639435695542, "viscosity"));
    
    // Kinematic viscosity (area/time, L^2 T^-1)
    Dimension kin_viscosity(2, 0, -1);
    registerUnit(Unit("square meter per second", "m2/s", kin_viscosity, 1.0, "kinematic_viscosity"));
    registerUnit(Unit("centistokes", "cSt", kin_viscosity, 1e-6, "kinematic_viscosity"));
    registerUnit(Unit("stokes", "St", kin_viscosity, 1e-4, "kinematic_viscosity"));
}

// =============================================================================
// Permeability Units
// =============================================================================

void UnitSystem::addPermeabilityUnits() {
    Dimension permeability(2, 0, 0);
    
    // SI
    registerUnit(Unit("square meter", "m2", permeability, 1.0, "permeability"));
    registerUnit(Unit("square centimeter", "cm2", permeability, 1e-4, "permeability"));
    
    // Darcy units (most common in reservoir engineering)
    registerUnit(Unit("darcy", "D", permeability, 9.869233e-13, "permeability"));
    registerUnit(Unit("millidarcy", "mD", permeability, 9.869233e-16, "permeability"));
    registerUnit(Unit("microdarcy", "uD", permeability, 9.869233e-19, "permeability"));
}

// =============================================================================
// Volumetric Rate Units
// =============================================================================

void UnitSystem::addVolumetricRateUnits() {
    Dimension rate(3, 0, -1);
    
    // SI
    registerUnit(Unit("cubic meter per second", "m3/s", rate, 1.0, "volumetric_rate"));
    registerUnit(Unit("cubic meter per day", "m3/day", rate, 1.0/86400.0, "volumetric_rate"));
    registerUnit(Unit("liter per second", "L/s", rate, 0.001, "volumetric_rate"));
    
    // Imperial/US
    registerUnit(Unit("cubic foot per second", "ft3/s", rate, 0.028316846592, "volumetric_rate"));
    registerUnit(Unit("cubic foot per day", "ft3/day", rate, 0.028316846592/86400.0, "volumetric_rate"));
    registerUnit(Unit("gallon per minute", "gpm", rate, 0.003785411784/60.0, "volumetric_rate"));
    
    // Oil field
    registerUnit(Unit("barrel per day", "bbl/day", rate, 0.158987294928/86400.0, "volumetric_rate"));
    registerUnit(Unit("stock tank barrel per day", "stb/day", rate, 0.158987294928/86400.0, "volumetric_rate"));
    registerUnit(Unit("thousand cubic feet per day", "Mcf/day", rate, 28.316846592/86400.0, "volumetric_rate"));
    registerUnit(Unit("million cubic feet per day", "MMcf/day", rate, 28316.846592/86400.0, "volumetric_rate"));
}

// =============================================================================
// Mass Rate Units
// =============================================================================

void UnitSystem::addMassRateUnits() {
    Dimension rate(0, 1, -1);
    
    registerUnit(Unit("kilogram per second", "kg/s", rate, 1.0, "mass_rate"));
    registerUnit(Unit("kilogram per day", "kg/day", rate, 1.0/86400.0, "mass_rate"));
    registerUnit(Unit("tonne per day", "tonne/day", rate, 1000.0/86400.0, "mass_rate"));
    registerUnit(Unit("pound mass per second", "lbm/s", rate, 0.45359237, "mass_rate"));
    registerUnit(Unit("pound mass per day", "lbm/day", rate, 0.45359237/86400.0, "mass_rate"));
}

// =============================================================================
// Temperature Units
// =============================================================================

void UnitSystem::addTemperatureUnits() {
    // Temperature has special handling due to offset (Celsius, Fahrenheit)
    Dimension temperature(0, 0, 0);  // Treated as base dimension
    
    // Absolute temperatures
    registerUnit(Unit("kelvin", "K", temperature, 1.0, "temperature"));
    registerUnit(Unit("rankine", "R", temperature, 5.0/9.0, "temperature"));
    
    // Relative temperatures (with offset)
    Unit celsius("celsius", "degC", temperature, 1.0, "temperature");
    celsius.offset = -273.15;  // 0°C = 273.15 K
    registerUnit(celsius);
    
    Unit fahrenheit("fahrenheit", "degF", temperature, 5.0/9.0, "temperature");
    fahrenheit.offset = -459.67 * 5.0/9.0;  // 0°F = 459.67°R = 255.372 K
    registerUnit(fahrenheit);
}

// =============================================================================
// Thermal Units
// =============================================================================

void UnitSystem::addThermalUnits() {
    // Thermal conductivity: W/(m·K) = kg·m/(s³·K)
    Dimension thermal_cond(1, 1, -3);  // Simplified: L M T^-3
    registerUnit(Unit("watt per meter kelvin", "W/(m-K)", thermal_cond, 1.0, "thermal_conductivity"));
    registerUnit(Unit("BTU per hour foot fahrenheit", "BTU/(hr-ft-degF)", 
                     thermal_cond, 1.7307346563862, "thermal_conductivity"));
    
    // Heat capacity: J/(kg·K) = m²/(s²·K)
    Dimension heat_cap(2, 0, -2);  // Simplified
    registerUnit(Unit("joule per kilogram kelvin", "J/(kg-K)", heat_cap, 1.0, "heat_capacity"));
    registerUnit(Unit("BTU per pound fahrenheit", "BTU/(lbm-degF)", heat_cap, 4186.8, "heat_capacity"));
    
    // Thermal expansion: 1/K = 1/T (dimensionless per temperature)
    Dimension thermal_exp(0, 0, 0);
    registerUnit(Unit("per kelvin", "1/K", thermal_exp, 1.0, "thermal_expansion"));
    registerUnit(Unit("per celsius", "1/degC", thermal_exp, 1.0, "thermal_expansion"));
    registerUnit(Unit("per fahrenheit", "1/degF", thermal_exp, 1.8, "thermal_expansion"));
}

// =============================================================================
// Stress/Modulus Units
// =============================================================================

void UnitSystem::addStressUnits() {
    // Stress has same dimension as pressure
    Dimension stress(−1, 1, -2);
    
    registerUnit(Unit("pascal stress", "Pa", stress, 1.0, "stress"));
    registerUnit(Unit("kilopascal stress", "kPa", stress, 1000.0, "stress"));
    registerUnit(Unit("megapascal stress", "MPa", stress, 1e6, "stress"));
    registerUnit(Unit("gigapascal stress", "GPa", stress, 1e9, "stress"));
    registerUnit(Unit("psi stress", "psi", stress, 6894.757293168, "stress"));
    registerUnit(Unit("ksi stress", "ksi", stress, 6894757.293168, "stress"));
}

void UnitSystem::addModulusUnits() {
    // Modulus has same dimension as pressure/stress
    Dimension modulus(−1, 1, -2);
    
    registerUnit(Unit("pascal modulus", "Pa", modulus, 1.0, "modulus"));
    registerUnit(Unit("kilopascal modulus", "kPa", modulus, 1000.0, "modulus"));
    registerUnit(Unit("megapascal modulus", "MPa", modulus, 1e6, "modulus"));
    registerUnit(Unit("gigapascal modulus", "GPa", modulus, 1e9, "modulus"));
    registerUnit(Unit("psi modulus", "psi", modulus, 6894.757293168, "modulus"));
}

// =============================================================================
// Fracture Toughness Units
// =============================================================================

void UnitSystem::addFractureToughnessUnits() {
    // Fracture toughness: Pa·m^0.5 = kg/(s²·m^0.5) dimension: L^0.5 M T^-2
    Dimension toughness(0.5, 1, -2);
    
    registerUnit(Unit("pascal square root meter", "Pa-m0.5", toughness, 1.0, "fracture_toughness"));
    registerUnit(Unit("megapascal square root meter", "MPa-m0.5", toughness, 1e6, "fracture_toughness"));
    registerUnit(Unit("psi square root inch", "psi-in0.5", toughness, 1099.685, "fracture_toughness"));
    registerUnit(Unit("ksi square root inch", "ksi-in0.5", toughness, 1.099685e6, "fracture_toughness"));
}

// =============================================================================
// Compressibility Units
// =============================================================================

void UnitSystem::addCompressibilityUnits() {
    // Compressibility: 1/Pa dimension: L T² M^-1
    Dimension compressibility(1, -1, 2);
    
    registerUnit(Unit("per pascal", "1/Pa", compressibility, 1.0, "compressibility"));
    registerUnit(Unit("per kilopascal", "1/kPa", compressibility, 0.001, "compressibility"));
    registerUnit(Unit("per megapascal", "1/MPa", compressibility, 1e-6, "compressibility"));
    registerUnit(Unit("per psi", "1/psi", compressibility, 1.0/6894.757293168, "compressibility"));
    registerUnit(Unit("per bar", "1/bar", compressibility, 1e-5, "compressibility"));
}

// =============================================================================
// Productivity Index Units
// =============================================================================

void UnitSystem::addProductivityUnits() {
    // Productivity index: (m³/s)/Pa = m³/(Pa·s) dimension: L^4 T M^-1
    Dimension productivity(4, -1, 1);
    
    registerUnit(Unit("cubic meter per second per pascal", "m3/(s-Pa)", 
                     productivity, 1.0, "productivity_index"));
    registerUnit(Unit("cubic meter per day per bar", "m3/(day-bar)", 
                     productivity, 1.0/(86400.0 * 1e5), "productivity_index"));
    registerUnit(Unit("barrel per day per psi", "bbl/(day-psi)", 
                     productivity, 0.158987294928/(86400.0 * 6894.757293168), "productivity_index"));
}

// =============================================================================
// Transmissibility Units
// =============================================================================

void UnitSystem::addTransmissibilityUnits() {
    // Transmissibility: k*A/μ*L dimension: L^3 (permeability*area/(viscosity*length))
    Dimension transmissibility(3, 0, 0);
    
    registerUnit(Unit("cubic meter transmissibility", "m3", transmissibility, 1.0, "transmissibility"));
    registerUnit(Unit("darcy foot squared per centipoise foot", "D-ft2/(cP-ft)", 
                     transmissibility, 9.869233e-13 * 0.09290304 / (0.001 * 0.3048), 
                     "transmissibility"));
}

// =============================================================================
// Helper Functions
// =============================================================================

void UnitSystem::registerUnit(const Unit& unit) {
    // Store by name (lowercase)
    std::string key = toLowerCase(unit.name);
    units_[key] = unit;
    
    // Store by symbol (case-sensitive primary, lowercase secondary)
    if (!unit.symbol.empty()) {
        units_[unit.symbol] = unit;
        units_[toLowerCase(unit.symbol)] = unit;
    }
    
    // Add aliases
    for (const auto& alias : unit.aliases) {
        units_[toLowerCase(alias)] = unit;
    }
    
    // Add to category index
    if (!unit.category.empty()) {
        categories_[unit.category].push_back(key);
    }
}

std::string UnitSystem::toLowerCase(const std::string& str) const {
    std::string result = str;
    std::transform(result.begin(), result.end(), result.begin(), ::tolower);
    return result;
}

std::string UnitSystem::trim(const std::string& str) const {
    size_t first = str.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    size_t last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, last - first + 1);
}

std::vector<std::string> UnitSystem::split(const std::string& str, char delim) const {
    std::vector<std::string> result;
    std::stringstream ss(str);
    std::string item;
    while (std::getline(ss, item, delim)) {
        item = trim(item);
        if (!item.empty()) {
            result.push_back(item);
        }
    }
    return result;
}

// =============================================================================
// Database Access
// =============================================================================

const Unit* UnitSystem::getUnit(const std::string& name_or_symbol) const {
    // Try exact match first
    auto it = units_.find(name_or_symbol);
    if (it != units_.end()) {
        return &(it->second);
    }
    
    // Try lowercase
    it = units_.find(toLowerCase(name_or_symbol));
    if (it != units_.end()) {
        return &(it->second);
    }
    
    return nullptr;
}

bool UnitSystem::hasUnit(const std::string& name_or_symbol) const {
    return getUnit(name_or_symbol) != nullptr;
}

std::vector<const Unit*> UnitSystem::getUnitsInCategory(const std::string& category) const {
    std::vector<const Unit*> result;
    auto it = categories_.find(category);
    if (it != categories_.end()) {
        for (const auto& unit_name : it->second) {
            auto unit_it = units_.find(unit_name);
            if (unit_it != units_.end()) {
                result.push_back(&(unit_it->second));
            }
        }
    }
    return result;
}

std::vector<std::string> UnitSystem::getCategories() const {
    std::vector<std::string> result;
    for (const auto& pair : categories_) {
        result.push_back(pair.first);
    }
    return result;
}

Dimension UnitSystem::getDimension(const std::string& unit_name) const {
    const Unit* unit = getUnit(unit_name);
    if (unit) {
        return unit->dimension;
    }
    throw std::runtime_error("Unit not found: " + unit_name);
}

// =============================================================================
// Conversion Functions
// =============================================================================

double UnitSystem::convert(double value, const std::string& from_unit, 
                           const std::string& to_unit) const {
    const Unit* from = getUnit(from_unit);
    const Unit* to = getUnit(to_unit);
    
    if (!from) {
        throw std::runtime_error("Unknown source unit: " + from_unit);
    }
    if (!to) {
        throw std::runtime_error("Unknown destination unit: " + to_unit);
    }
    
    if (from->dimension != to->dimension) {
        throw std::runtime_error("Incompatible dimensions: " + 
                                from->dimension.toString() + " vs " + 
                                to->dimension.toString());
    }
    
    // Convert: from_unit -> base -> to_unit
    double base_value = from->convertToBase(value);
    return to->convertFromBase(base_value);
}

double UnitSystem::toBase(double value, const std::string& from_unit) const {
    const Unit* unit = getUnit(from_unit);
    if (!unit) {
        throw std::runtime_error("Unknown unit: " + from_unit);
    }
    return unit->convertToBase(value);
}

double UnitSystem::fromBase(double value, const std::string& to_unit) const {
    const Unit* unit = getUnit(to_unit);
    if (!unit) {
        throw std::runtime_error("Unknown unit: " + to_unit);
    }
    return unit->convertFromBase(value);
}

std::vector<double> UnitSystem::convert(const std::vector<double>& values,
                                       const std::string& from_unit,
                                       const std::string& to_unit) const {
    std::vector<double> result;
    result.reserve(values.size());
    for (double value : values) {
        result.push_back(convert(value, from_unit, to_unit));
    }
    return result;
}

// =============================================================================
// Parsing Functions
// =============================================================================

bool UnitSystem::parseValueWithUnit(const std::string& value_with_unit,
                                    double& value, std::string& unit) const {
    std::string trimmed = trim(value_with_unit);
    if (trimmed.empty()) return false;
    
    // Find where the number ends and unit begins
    size_t i = 0;
    
    // Skip sign
    if (trimmed[i] == '+' || trimmed[i] == '-') i++;
    
    // Skip digits and decimal point
    bool has_digits = false;
    bool has_decimal = false;
    while (i < trimmed.length()) {
        if (std::isdigit(trimmed[i])) {
            has_digits = true;
            i++;
        } else if (trimmed[i] == '.' && !has_decimal) {
            has_decimal = true;
            i++;
        } else if (trimmed[i] == 'e' || trimmed[i] == 'E') {
            // Scientific notation
            i++;
            if (i < trimmed.length() && (trimmed[i] == '+' || trimmed[i] == '-')) {
                i++;
            }
        } else {
            break;
        }
    }
    
    if (!has_digits) return false;
    
    // Extract number and unit parts
    std::string num_str = trim(trimmed.substr(0, i));
    std::string unit_str = trim(trimmed.substr(i));
    
    try {
        value = std::stod(num_str);
        unit = unit_str;
        return true;
    } catch (...) {
        return false;
    }
}

double UnitSystem::parseAndConvertToBase(const std::string& value_with_unit) const {
    double value;
    std::string unit;
    
    if (!parseValueWithUnit(value_with_unit, value, unit)) {
        throw std::runtime_error("Failed to parse: " + value_with_unit);
    }
    
    if (unit.empty()) {
        // No unit specified, assume base SI
        return value;
    }
    
    return toBase(value, unit);
}

// =============================================================================
// Dimensional Analysis
// =============================================================================

bool UnitSystem::areCompatible(const std::string& unit1, const std::string& unit2) const {
    const Unit* u1 = getUnit(unit1);
    const Unit* u2 = getUnit(unit2);
    
    if (!u1 || !u2) return false;
    return u1->dimension == u2->dimension;
}

std::string UnitSystem::getBaseUnit(const Dimension& dim) const {
    // Construct base SI unit string from dimension
    std::stringstream ss;
    bool first = true;
    
    if (std::abs(dim.L) > 1e-10) {
        ss << "m";
        if (std::abs(dim.L - 1.0) > 1e-10) ss << "^" << dim.L;
        first = false;
    }
    
    if (std::abs(dim.M) > 1e-10) {
        if (!first) ss << " ";
        ss << "kg";
        if (std::abs(dim.M - 1.0) > 1e-10) ss << "^" << dim.M;
        first = false;
    }
    
    if (std::abs(dim.T) > 1e-10) {
        if (!first) ss << " ";
        ss << "s";
        if (std::abs(dim.T - 1.0) > 1e-10) ss << "^" << dim.T;
    }
    
    return ss.str().empty() ? "dimensionless" : ss.str();
}

std::string UnitSystem::getSuggestedDisplayUnit(const std::string& quantity) const {
    auto it = display_units_.find(quantity);
    if (it != display_units_.end()) {
        return it->second;
    }
    return "";
}

// =============================================================================
// Custom Units
// =============================================================================

void UnitSystem::addUnit(const Unit& unit) {
    registerUnit(unit);
}

void UnitSystem::addAlias(const std::string& unit_name, const std::string& alias) {
    const Unit* unit = getUnit(unit_name);
    if (unit) {
        Unit modified = *unit;
        modified.aliases.push_back(alias);
        registerUnit(modified);
    }
}

// =============================================================================
// Utility Functions
// =============================================================================

std::string UnitSystem::formatValue(double value, const std::string& unit, 
                                   int precision) const {
    std::stringstream ss;
    ss << std::setprecision(precision) << value << " " << unit;
    return ss.str();
}

void UnitSystem::printDatabase(std::ostream& os) const {
    os << "Unit System Database\n";
    os << "====================\n\n";
    
    for (const auto& cat_pair : categories_) {
        os << "Category: " << cat_pair.first << "\n";
        os << std::string(40, '-') << "\n";
        
        for (const auto& unit_name : cat_pair.second) {
            auto it = units_.find(unit_name);
            if (it != units_.end()) {
                const Unit& u = it->second;
                os << std::setw(25) << std::left << u.name 
                   << " [" << std::setw(10) << u.symbol << "] "
                   << " = " << u.to_base << " * base SI"
                   << " (" << u.dimension.toString() << ")\n";
            }
        }
        os << "\n";
    }
}

std::string UnitSystem::generateDocumentation() const {
    std::stringstream ss;
    
    ss << "# FSRM Unit System Documentation\n\n";
    ss << "This document lists all supported units in the FSRM simulation framework.\n";
    ss << "All calculations are performed internally in SI base units (m, kg, s).\n\n";
    
    for (const auto& cat_pair : categories_) {
        ss << "## " << cat_pair.first << "\n\n";
        ss << "| Name | Symbol | Conversion to SI | Dimension |\n";
        ss << "|------|--------|------------------|------------|\n";
        
        for (const auto& unit_name : cat_pair.second) {
            auto it = units_.find(unit_name);
            if (it != units_.end()) {
                const Unit& u = it->second;
                ss << "| " << u.name 
                   << " | " << u.symbol 
                   << " | " << u.to_base << " × base"
                   << " | " << u.dimension.toString() << " |\n";
            }
        }
        ss << "\n";
    }
    
    return ss.str();
}

} // namespace FSRM
