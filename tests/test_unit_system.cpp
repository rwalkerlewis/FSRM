#include "UnitSystem.hpp"
#include "Testing.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace FSRM;

// Helper function for approximate equality
bool approxEqual(double a, double b, double tol = 1e-9) {
    return std::abs(a - b) < tol;
}

void testLengthConversions() {
    std::cout << "Testing Length Conversions..." << std::endl;
    UnitSystem units;
    
    // m to ft
    double ft = units.convert(1.0, "m", "ft");
    assert(approxEqual(ft, 3.28084, 1e-5));
    
    // km to m
    double m = units.convert(1.0, "km", "m");
    assert(approxEqual(m, 1000.0));
    
    // in to cm
    double cm = units.convert(1.0, "in", "cm");
    assert(approxEqual(cm, 2.54));
    
    // mi to km
    double km = units.convert(1.0, "mi", "km");
    assert(approxEqual(km, 1.609344));
    
    std::cout << "  ✓ Length conversions passed" << std::endl;
}

void testPressureConversions() {
    std::cout << "Testing Pressure Conversions..." << std::endl;
    UnitSystem units;
    
    // psi to Pa
    double pa = units.convert(1.0, "psi", "Pa");
    assert(approxEqual(pa, 6894.757293168, 1e-6));
    
    // bar to Pa
    pa = units.convert(1.0, "bar", "Pa");
    assert(approxEqual(pa, 100000.0));
    
    // MPa to psi
    double psi = units.convert(1.0, "MPa", "psi");
    assert(approxEqual(psi, 145.03774, 1e-4));
    
    // atm to Pa
    pa = units.convert(1.0, "atm", "Pa");
    assert(approxEqual(pa, 101325.0));
    
    std::cout << "  ✓ Pressure conversions passed" << std::endl;
}

void testPermeabilityConversions() {
    std::cout << "Testing Permeability Conversions..." << std::endl;
    UnitSystem units;
    
    // mD to m²
    double m2 = units.convert(1.0, "mD", "m2");
    assert(approxEqual(m2, 9.869233e-16, 1e-22));
    
    // D to mD
    double md = units.convert(1.0, "D", "mD");
    assert(approxEqual(md, 1000.0));
    
    // mD to D
    double d = units.convert(100.0, "mD", "D");
    assert(approxEqual(d, 0.1));
    
    std::cout << "  ✓ Permeability conversions passed" << std::endl;
}

void testViscosityConversions() {
    std::cout << "Testing Viscosity Conversions..." << std::endl;
    UnitSystem units;
    
    // cP to Pa·s
    double pas = units.convert(1.0, "cP", "Pa-s");
    assert(approxEqual(pas, 0.001));
    
    // P to Pa·s
    pas = units.convert(1.0, "P", "Pa-s");
    assert(approxEqual(pas, 0.1));
    
    // cP to P
    double p = units.convert(100.0, "cP", "P");
    assert(approxEqual(p, 1.0));
    
    std::cout << "  ✓ Viscosity conversions passed" << std::endl;
}

void testVolumeConversions() {
    std::cout << "Testing Volume Conversions..." << std::endl;
    UnitSystem units;
    
    // bbl to m³
    double m3 = units.convert(1.0, "bbl", "m3");
    assert(approxEqual(m3, 0.158987294928, 1e-9));
    
    // gal to L
    double l = units.convert(1.0, "gal", "L");
    assert(approxEqual(l, 3.785411784, 1e-8));
    
    // ft³ to m³
    m3 = units.convert(1.0, "ft3", "m3");
    assert(approxEqual(m3, 0.028316846592, 1e-9));
    
    std::cout << "  ✓ Volume conversions passed" << std::endl;
}

void testTemperatureConversions() {
    std::cout << "Testing Temperature Conversions..." << std::endl;
    UnitSystem units;
    
    // Celsius to Kelvin
    double k = units.convert(0.0, "degC", "K");
    assert(approxEqual(k, 273.15, 1e-6));
    
    k = units.convert(100.0, "degC", "K");
    assert(approxEqual(k, 373.15, 1e-6));
    
    // Fahrenheit to Kelvin
    k = units.convert(32.0, "degF", "K");
    assert(approxEqual(k, 273.15, 1e-6));
    
    k = units.convert(212.0, "degF", "K");
    assert(approxEqual(k, 373.15, 1e-6));
    
    // Celsius to Fahrenheit
    double f = units.convert(0.0, "degC", "degF");
    assert(approxEqual(f, 32.0, 1e-6));
    
    f = units.convert(100.0, "degC", "degF");
    assert(approxEqual(f, 212.0, 1e-6));
    
    std::cout << "  ✓ Temperature conversions passed" << std::endl;
}

void testTimeConversions() {
    std::cout << "Testing Time Conversions..." << std::endl;
    UnitSystem units;
    
    // day to s
    double s = units.convert(1.0, "day", "s");
    assert(approxEqual(s, 86400.0));
    
    // hr to s
    s = units.convert(1.0, "hr", "s");
    assert(approxEqual(s, 3600.0));
    
    // year to day
    double day = units.convert(1.0, "year", "day");
    assert(approxEqual(day, 365.0, 1.0));  // Approximate
    
    std::cout << "  ✓ Time conversions passed" << std::endl;
}

void testRateConversions() {
    std::cout << "Testing Rate Conversions..." << std::endl;
    UnitSystem units;
    
    // bbl/day to m³/s
    double m3s = units.convert(1.0, "bbl/day", "m3/s");
    assert(approxEqual(m3s, 1.8401307e-6, 1e-12));
    
    // stb/day to m³/s
    m3s = units.convert(1000.0, "stb/day", "m3/s");
    assert(approxEqual(m3s, 1.8401307e-3, 1e-9));
    
    std::cout << "  ✓ Rate conversions passed" << std::endl;
}

void testDensityConversions() {
    std::cout << "Testing Density Conversions..." << std::endl;
    UnitSystem units;
    
    // g/cm³ to kg/m³
    double kgm3 = units.convert(1.0, "g/cm3", "kg/m3");
    assert(approxEqual(kgm3, 1000.0));
    
    // lbm/ft³ to kg/m³
    kgm3 = units.convert(1.0, "lbm/ft3", "kg/m3");
    assert(approxEqual(kgm3, 16.018463, 1e-5));
    
    std::cout << "  ✓ Density conversions passed" << std::endl;
}

void testAngleConversions() {
    std::cout << "Testing Angle Conversions..." << std::endl;
    UnitSystem units;
    
    // deg to rad
    double rad = units.convert(180.0, "deg", "rad");
    assert(approxEqual(rad, M_PI, 1e-10));
    
    // deg to rad (90°)
    rad = units.convert(90.0, "deg", "rad");
    assert(approxEqual(rad, M_PI/2.0, 1e-10));
    
    std::cout << "  ✓ Angle conversions passed" << std::endl;
}

void testParsingWithUnits() {
    std::cout << "Testing Parsing with Units..." << std::endl;
    UnitSystem units;
    
    // Parse "5000 psi"
    double value;
    std::string unit;
    bool success = units.parseValueWithUnit("5000 psi", value, unit);
    assert(success);
    assert(approxEqual(value, 5000.0));
    assert(unit == "psi");
    
    // Parse "100.5 mD"
    success = units.parseValueWithUnit("100.5 mD", value, unit);
    assert(success);
    assert(approxEqual(value, 100.5));
    assert(unit == "mD");
    
    // Parse with scientific notation
    success = units.parseValueWithUnit("1.5e-6 m2", value, unit);
    assert(success);
    assert(approxEqual(value, 1.5e-6));
    assert(unit == "m2");
    
    // Parse and convert to base
    double base = units.parseAndConvertToBase("5000 psi");
    assert(approxEqual(base, 5000.0 * 6894.757293168, 1e-3));
    
    std::cout << "  ✓ Parsing tests passed" << std::endl;
}

void testDimensionalAnalysis() {
    std::cout << "Testing Dimensional Analysis..." << std::endl;
    UnitSystem units;
    
    // Compatible units
    assert(units.areCompatible("psi", "Pa"));
    assert(units.areCompatible("mD", "m2"));
    assert(units.areCompatible("cP", "Pa-s"));
    assert(units.areCompatible("ft", "m"));
    
    // Incompatible units
    assert(!units.areCompatible("psi", "mD"));
    assert(!units.areCompatible("m", "s"));
    assert(!units.areCompatible("Pa", "cP"));
    
    std::cout << "  ✓ Dimensional analysis passed" << std::endl;
}

void testIncompatibleConversion() {
    std::cout << "Testing Incompatible Conversion Error..." << std::endl;
    UnitSystem units;
    
    bool caught_error = false;
    try {
        // This should throw an exception
        double bad = units.convert(100.0, "mD", "psi");
    } catch (const std::runtime_error& e) {
        caught_error = true;
        std::cout << "  Expected error caught: " << e.what() << std::endl;
    }
    
    assert(caught_error);
    std::cout << "  ✓ Error handling passed" << std::endl;
}

void testToBaseAndFromBase() {
    std::cout << "Testing toBase and fromBase..." << std::endl;
    UnitSystem units;
    
    // Convert to base (SI)
    double base = units.toBase(5000.0, "psi");
    assert(approxEqual(base, 5000.0 * 6894.757293168, 1e-3));
    
    // Convert back from base
    double psi = units.fromBase(base, "psi");
    assert(approxEqual(psi, 5000.0, 1e-6));
    
    // Test with permeability
    base = units.toBase(100.0, "mD");
    assert(approxEqual(base, 100.0 * 9.869233e-16, 1e-22));
    
    double md = units.fromBase(base, "mD");
    assert(approxEqual(md, 100.0, 1e-6));
    
    std::cout << "  ✓ toBase/fromBase tests passed" << std::endl;
}

void testFieldUnitsWorkflow() {
    std::cout << "\nTesting Realistic Field Units Workflow..." << std::endl;
    UnitSystem units;
    
    // Typical reservoir simulation inputs
    std::cout << "  Converting field units to SI:" << std::endl;
    
    double perm_md = 150.0;  // millidarcy
    double perm_si = units.toBase(perm_md, "mD");
    std::cout << "    Permeability: " << perm_md << " mD = " 
              << std::scientific << perm_si << " m²" << std::endl;
    
    double visc_cp = 5.0;  // centipoise
    double visc_si = units.toBase(visc_cp, "cP");
    std::cout << "    Viscosity: " << visc_cp << " cP = " 
              << visc_si << " Pa·s" << std::endl;
    
    double pres_psi = 5000.0;  // psi
    double pres_si = units.toBase(pres_psi, "psi");
    std::cout << "    Pressure: " << pres_psi << " psi = " 
              << std::fixed << std::setprecision(0) << pres_si << " Pa" << std::endl;
    
    double rate_bbl = 2000.0;  // bbl/day
    double rate_si = units.toBase(rate_bbl, "bbl/day");
    std::cout << "    Rate: " << rate_bbl << " bbl/day = " 
              << std::scientific << rate_si << " m³/s" << std::endl;
    
    // Convert back for display
    std::cout << "\n  Converting SI back to field units for display:" << std::endl;
    double perm_display = units.fromBase(perm_si, "mD");
    std::cout << "    Permeability: " << perm_display << " mD (recovered)" << std::endl;
    
    double pres_mpa = units.fromBase(pres_si, "MPa");
    std::cout << "    Pressure: " << std::fixed << std::setprecision(2) 
              << pres_mpa << " MPa" << std::endl;
    
    std::cout << "  ✓ Field units workflow passed" << std::endl;
}

void testUnitCategories() {
    std::cout << "\nTesting Unit Categories..." << std::endl;
    UnitSystem units;
    
    auto categories = units.getCategories();
    std::cout << "  Available categories (" << categories.size() << "):" << std::endl;
    
    for (const auto& cat : categories) {
        auto cat_units = units.getUnitsInCategory(cat);
        std::cout << "    " << cat << ": " << cat_units.size() << " units" << std::endl;
    }
    
    // Check some specific categories exist
    assert(units.getUnitsInCategory("pressure").size() > 10);
    assert(units.getUnitsInCategory("permeability").size() >= 4);
    assert(units.getUnitsInCategory("viscosity").size() >= 5);
    
    std::cout << "  ✓ Category tests passed" << std::endl;
}

int main() {
    std::cout << "\n";
    std::cout << "========================================" << std::endl;
    std::cout << "  FSRM Unit System Test Suite" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "\n";
    
    try {
        testLengthConversions();
        testPressureConversions();
        testPermeabilityConversions();
        testViscosityConversions();
        testVolumeConversions();
        testTemperatureConversions();
        testTimeConversions();
        testRateConversions();
        testDensityConversions();
        testAngleConversions();
        testParsingWithUnits();
        testDimensionalAnalysis();
        testIncompatibleConversion();
        testToBaseAndFromBase();
        testFieldUnitsWorkflow();
        testUnitCategories();
        
        std::cout << "\n";
        std::cout << "========================================" << std::endl;
        std::cout << "  ✓ ALL TESTS PASSED!" << std::endl;
        std::cout << "========================================" << std::endl;
        std::cout << "\n";
        
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "\n";
        std::cerr << "========================================" << std::endl;
        std::cerr << "  ✗ TEST FAILED!" << std::endl;
        std::cerr << "  Error: " << e.what() << std::endl;
        std::cerr << "========================================" << std::endl;
        std::cerr << "\n";
        return 1;
    }
}
