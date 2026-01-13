/**
 * @file unit_converter.cpp
 * @brief Simple command-line unit converter utility
 * 
 * This utility demonstrates the FSRM unit system and provides
 * a convenient tool for quick unit conversions.
 * 
 * Usage:
 *   ./unit_converter <value> <from_unit> <to_unit>
 *   ./unit_converter --list [category]
 *   ./unit_converter --help
 * 
 * Examples:
 *   ./unit_converter 5000 psi MPa
 *   ./unit_converter 150 mD m2
 *   ./unit_converter 100 degC degF
 *   ./unit_converter --list pressure
 */

#include "UnitSystem.hpp"
#include <iostream>
#include <iomanip>
#include <string>
#include <cstring>

using namespace FSRM;

void printHelp() {
    std::cout << "\n";
    std::cout << "FSRM Unit Converter\n";
    std::cout << "===================\n\n";
    std::cout << "Usage:\n";
    std::cout << "  unit_converter <value> <from_unit> <to_unit>\n";
    std::cout << "  unit_converter --list [category]\n";
    std::cout << "  unit_converter --help\n\n";
    std::cout << "Examples:\n";
    std::cout << "  unit_converter 5000 psi MPa\n";
    std::cout << "  unit_converter 150 mD m2\n";
    std::cout << "  unit_converter 100 degC degF\n";
    std::cout << "  unit_converter 30 day s\n";
    std::cout << "  unit_converter --list\n";
    std::cout << "  unit_converter --list pressure\n\n";
    std::cout << "Common Units:\n";
    std::cout << "  Pressure: Pa, kPa, MPa, GPa, psi, bar, atm\n";
    std::cout << "  Permeability: m2, D, mD, uD\n";
    std::cout << "  Viscosity: Pa-s, cP, P\n";
    std::cout << "  Length: m, cm, mm, km, ft, in\n";
    std::cout << "  Volume: m3, L, bbl, ft3, gal, scf\n";
    std::cout << "  Rate: m3/s, bbl/day, stb/day, Mcf/day\n";
    std::cout << "  Time: s, min, hr, day, year\n";
    std::cout << "  Temperature: K, degC, degF\n";
    std::cout << "  Density: kg/m3, g/cm3, lbm/ft3\n\n";
}

void listUnits(UnitSystem& units, const std::string& category = "") {
    std::cout << "\n";
    
    if (category.empty()) {
        // List all categories
        std::cout << "Available Unit Categories:\n";
        std::cout << "==========================\n\n";
        
        auto categories = units.getCategories();
        for (const auto& cat : categories) {
            auto cat_units = units.getUnitsInCategory(cat);
            std::cout << std::setw(25) << std::left << cat 
                     << " (" << cat_units.size() << " units)\n";
        }
        std::cout << "\nUse: unit_converter --list <category> to see units in a category\n\n";
    } else {
        // List units in specific category
        auto cat_units = units.getUnitsInCategory(category);
        
        if (cat_units.empty()) {
            std::cout << "Category '" << category << "' not found.\n";
            std::cout << "Use: unit_converter --list to see available categories\n\n";
            return;
        }
        
        std::cout << "Units in category: " << category << "\n";
        std::cout << std::string(50, '=') << "\n\n";
        std::cout << std::setw(25) << std::left << "Name" 
                 << std::setw(12) << "Symbol" 
                 << "To SI Base\n";
        std::cout << std::string(50, '-') << "\n";
        
        for (const auto* unit : cat_units) {
            std::cout << std::setw(25) << std::left << unit->name
                     << std::setw(12) << unit->symbol
                     << std::scientific << std::setprecision(6) << unit->to_base;
            if (unit->offset != 0.0) {
                std::cout << " (offset: " << unit->offset << ")";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }
}

void performConversion(UnitSystem& units, double value, 
                      const std::string& from_unit, 
                      const std::string& to_unit) {
    try {
        // Check if units exist
        if (!units.hasUnit(from_unit)) {
            std::cerr << "Error: Unknown source unit '" << from_unit << "'\n";
            std::cerr << "Use --list to see available units\n";
            return;
        }
        
        if (!units.hasUnit(to_unit)) {
            std::cerr << "Error: Unknown destination unit '" << to_unit << "'\n";
            std::cerr << "Use --list to see available units\n";
            return;
        }
        
        // Check compatibility
        if (!units.areCompatible(from_unit, to_unit)) {
            auto from_dim = units.getDimension(from_unit);
            auto to_dim = units.getDimension(to_unit);
            std::cerr << "Error: Incompatible units\n";
            std::cerr << "  " << from_unit << " has dimension: " << from_dim.toString() << "\n";
            std::cerr << "  " << to_unit << " has dimension: " << to_dim.toString() << "\n";
            return;
        }
        
        // Perform conversion
        double result = units.convert(value, from_unit, to_unit);
        
        // Also show SI base value
        double si_value = units.toBase(value, from_unit);
        std::string si_unit = units.getBaseUnit(units.getDimension(from_unit));
        
        // Print results
        std::cout << "\n";
        std::cout << "Conversion Result:\n";
        std::cout << "==================\n\n";
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "  Input:   " << value << " " << from_unit << "\n";
        std::cout << "  Output:  " << result << " " << to_unit << "\n";
        std::cout << "\n";
        std::cout << std::scientific << std::setprecision(6);
        std::cout << "  SI Base: " << si_value << " " << si_unit << "\n";
        std::cout << "\n";
        
        // Show conversion factor
        double factor = result / value;
        std::cout << "Conversion Factor: 1 " << from_unit << " = " 
                 << std::scientific << std::setprecision(9) << factor 
                 << " " << to_unit << "\n\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
    }
}

int main(int argc, char* argv[]) {
    UnitSystem& units = UnitSystemManager::getInstance();
    
    // Parse command line arguments
    if (argc == 1 || (argc == 2 && strcmp(argv[1], "--help") == 0)) {
        printHelp();
        return 0;
    }
    
    if (argc >= 2 && strcmp(argv[1], "--list") == 0) {
        if (argc == 2) {
            listUnits(units);
        } else {
            listUnits(units, argv[2]);
        }
        return 0;
    }
    
    if (argc != 4) {
        std::cerr << "Error: Invalid number of arguments\n";
        printHelp();
        return 1;
    }
    
    // Parse conversion arguments
    try {
        double value = std::stod(argv[1]);
        std::string from_unit = argv[2];
        std::string to_unit = argv[3];
        
        performConversion(units, value, from_unit, to_unit);
        
    } catch (const std::invalid_argument& e) {
        std::cerr << "Error: Invalid value '" << argv[1] << "'\n";
        return 1;
    } catch (const std::out_of_range& e) {
        std::cerr << "Error: Value out of range\n";
        return 1;
    }
    
    return 0;
}
