# FSRM Unit System Documentation

## Overview

The FSRM (Fault-Slip Reservoir Modeling) framework includes a comprehensive unit system that allows users to specify input parameters and output results in any supported unit system. 

**All calculations are performed internally using SI base units:**
- **Length**: meters (m)
- **Mass**: kilograms (kg)  
- **Time**: seconds (s)

This design ensures:
- Numerical consistency across all physics models
- No unit conversion errors during computation
- Easy verification against analytical solutions
- Compatibility with PETSc and other numerical libraries

## Key Features

- **Automatic Conversion**: Values are automatically converted to SI base units for calculations
- **Flexible Input**: Specify values with units directly in config files (e.g., `5000 psi`, `100 mD`)
- **Comprehensive Database**: Support for 150+ units across all physical quantities
- **Dimensional Analysis**: Automatic validation of unit compatibility
- **User-Friendly Output**: Display results in preferred units (field units, metric, etc.)
- **No Code Changes Required**: Unit conversion happens transparently

## SI Base Units (LMT System)

All internal calculations use the Length-Mass-Time (LMT) base unit system:

| Dimension | SI Base Unit | Symbol |
|-----------|--------------|--------|
| Length    | meter        | m      |
| Mass      | kilogram     | kg     |
| Time      | second       | s      |

All other units are derived from these base units or converted to them.

## Usage in Configuration Files

### Basic Syntax

```ini
parameter = value unit
```

### Examples

```ini
[ROCK]
permeability_x = 150 mD              # Millidarcy → m²
youngs_modulus = 15 GPa              # Gigapascal → Pa
density = 2.55 g/cm3                 # g/cm³ → kg/m³

[FLUID]
viscosity = 5 cP                     # Centipoise → Pa·s
density = 850 kg/m3                  # Already in SI

[WELL1]
target_value = 5000 bbl/day          # Barrels/day → m³/s
min_bhp = 2000 psi                   # psi → Pa

[SIMULATION]
end_time = 30 day                    # days → seconds
dt_initial = 1 hr                    # hours → seconds

[GRID]
Lx = 2 km                            # kilometers → meters
Ly = 5000 ft                         # feet → meters
```

### Default Units

If no unit is specified, the system uses default units (usually SI). You can specify default units programmatically:

```cpp
double perm = config.getDoubleWithUnit("ROCK", "permeability", 100e-15, "mD");
// If config has "permeability = 150", it's interpreted as 150 mD
// If config has "permeability = 150 D", it's interpreted as 150 Darcy
```

## Supported Units by Category

### Length Units

| Unit | Symbol | Conversion to m |
|------|--------|----------------|
| meter | m | 1.0 |
| centimeter | cm | 0.01 |
| millimeter | mm | 0.001 |
| kilometer | km | 1000.0 |
| micrometer | um | 1×10⁻⁶ |
| nanometer | nm | 1×10⁻⁹ |
| foot | ft | 0.3048 |
| inch | in | 0.0254 |
| yard | yd | 0.9144 |
| mile | mi | 1609.344 |
| angstrom | angstrom | 1×10⁻¹⁰ |

### Mass Units

| Unit | Symbol | Conversion to kg |
|------|--------|-----------------|
| kilogram | kg | 1.0 |
| gram | g | 0.001 |
| milligram | mg | 1×10⁻⁶ |
| tonne | tonne | 1000.0 |
| metric ton | t | 1000.0 |
| pound mass | lbm | 0.45359237 |
| ounce | oz | 0.028349523 |
| slug | slug | 14.5939029 |
| short ton | ton | 907.18474 |

### Time Units

| Unit | Symbol | Conversion to s |
|------|--------|----------------|
| second | s | 1.0 |
| millisecond | ms | 0.001 |
| microsecond | us | 1×10⁻⁶ |
| minute | min | 60.0 |
| hour | hr | 3600.0 |
| day | day | 86400.0 |
| week | week | 604800.0 |
| year | year | 3.1536×10⁷ |

### Pressure Units

| Unit | Symbol | Conversion to Pa |
|------|--------|-----------------|
| pascal | Pa | 1.0 |
| kilopascal | kPa | 1000.0 |
| megapascal | MPa | 1×10⁶ |
| gigapascal | GPa | 1×10⁹ |
| bar | bar | 1×10⁵ |
| millibar | mbar | 100.0 |
| atmosphere | atm | 101325.0 |
| pounds per square inch | psi | 6894.757 |
| pounds per square foot | psf | 47.88026 |
| ksi | ksi | 6.894757×10⁶ |
| torr | torr | 133.3224 |
| millimeter mercury | mmHg | 133.3224 |
| inch mercury | inHg | 3386.389 |

### Permeability Units

| Unit | Symbol | Conversion to m² |
|------|--------|-----------------|
| square meter | m2 | 1.0 |
| darcy | D | 9.869233×10⁻¹³ |
| millidarcy | mD | 9.869233×10⁻¹⁶ |
| microdarcy | uD | 9.869233×10⁻¹⁹ |

### Viscosity Units (Dynamic)

| Unit | Symbol | Conversion to Pa·s |
|------|--------|-------------------|
| pascal second | Pa-s | 1.0 |
| centipoise | cP | 0.001 |
| poise | P | 0.1 |
| millipascal second | mPa-s | 0.001 |
| pound per foot second | lb/ft-s | 1.488164 |

### Viscosity Units (Kinematic)

| Unit | Symbol | Conversion to m²/s |
|------|--------|--------------------|
| square meter per second | m2/s | 1.0 |
| centistokes | cSt | 1×10⁻⁶ |
| stokes | St | 1×10⁻⁴ |

### Volume Units

| Unit | Symbol | Conversion to m³ |
|------|--------|-----------------|
| cubic meter | m3 | 1.0 |
| cubic centimeter | cm3 | 1×10⁻⁶ |
| liter | L | 0.001 |
| milliliter | mL | 1×10⁻⁶ |
| cubic foot | ft3 | 0.028316847 |
| cubic inch | in3 | 1.6387064×10⁻⁵ |
| gallon US | gal | 0.003785412 |
| barrel | bbl | 0.158987295 |
| stock tank barrel | stb | 0.158987295 |
| standard cubic foot | scf | 0.028316847 |
| thousand cubic feet | Mcf | 28.316847 |
| million cubic feet | MMcf | 28316.847 |

### Volumetric Rate Units

| Unit | Symbol | Conversion to m³/s |
|------|--------|--------------------|
| cubic meter per second | m3/s | 1.0 |
| cubic meter per day | m3/day | 1.1574×10⁻⁵ |
| liter per second | L/s | 0.001 |
| cubic foot per second | ft3/s | 0.028316847 |
| gallon per minute | gpm | 6.3090197×10⁻⁵ |
| barrel per day | bbl/day | 1.8401307×10⁻⁶ |
| stock tank barrel per day | stb/day | 1.8401307×10⁻⁶ |
| thousand cubic feet per day | Mcf/day | 3.2774128×10⁻⁴ |
| million cubic feet per day | MMcf/day | 0.32774128 |

### Temperature Units

| Unit | Symbol | Notes |
|------|--------|-------|
| kelvin | K | Absolute temperature |
| celsius | degC | 0°C = 273.15 K |
| fahrenheit | degF | 32°F = 273.15 K |
| rankine | R | Absolute (°F + 459.67) |

**Note**: Temperature conversions handle both absolute values and differences correctly.

### Energy Units

| Unit | Symbol | Conversion to J |
|------|--------|----------------|
| joule | J | 1.0 |
| kilojoule | kJ | 1000.0 |
| megajoule | MJ | 1×10⁶ |
| erg | erg | 1×10⁻⁷ |
| calorie | cal | 4.184 |
| kilocalorie | kcal | 4184.0 |
| british thermal unit | BTU | 1055.056 |
| foot pound force | ft-lbf | 1.355818 |
| watt hour | Wh | 3600.0 |
| kilowatt hour | kWh | 3.6×10⁶ |

### Density Units

| Unit | Symbol | Conversion to kg/m³ |
|------|--------|---------------------|
| kilogram per cubic meter | kg/m3 | 1.0 |
| gram per cubic centimeter | g/cm3 | 1000.0 |
| pound mass per cubic foot | lbm/ft3 | 16.018463 |
| pound mass per gallon | lbm/gal | 119.82643 |

### Angle Units

| Unit | Symbol | Conversion to rad |
|------|--------|------------------|
| radian | rad | 1.0 |
| degree | deg | π/180 |
| gradian | grad | π/200 |

### Stress/Modulus Units

Same as pressure units (Pa, kPa, MPa, GPa, psi, ksi, etc.)

### Fracture Toughness Units

| Unit | Symbol | Conversion to Pa·m^0.5 |
|------|--------|------------------------|
| pascal square root meter | Pa-m0.5 | 1.0 |
| megapascal square root meter | MPa-m0.5 | 1×10⁶ |
| psi square root inch | psi-in0.5 | 1099.685 |
| ksi square root inch | ksi-in0.5 | 1.099685×10⁶ |

### Compressibility Units

| Unit | Symbol | Conversion to 1/Pa |
|------|--------|--------------------|
| per pascal | 1/Pa | 1.0 |
| per kilopascal | 1/kPa | 0.001 |
| per megapascal | 1/MPa | 1×10⁻⁶ |
| per psi | 1/psi | 1.4503774×10⁻⁴ |
| per bar | 1/bar | 1×10⁻⁵ |

## Output Unit Configuration

You can specify preferred units for output/visualization:

```ini
[OUTPUT]
format = VTK
frequency = 10

# Specify output units
pressure_unit = psi                  # Display pressure in psi
displacement_unit = mm               # Display displacement in mm
stress_unit = MPa                    # Display stress in MPa
permeability_unit = mD               # Display permeability in mD
temperature_unit = degC              # Display temperature in Celsius
viscosity_unit = cP                  # Display viscosity in cP
density_unit = kg/m3                 # Display density in kg/m³
time_unit = day                      # Display time in days
```

## Programmatic Usage

### C++ API

```cpp
#include "UnitSystem.hpp"

using namespace FSRM;

// Get the global unit system instance
UnitSystem& units = UnitSystemManager::getInstance();

// Convert between units
double pressure_pa = units.convert(5000.0, "psi", "Pa");
double perm_m2 = units.convert(150.0, "mD", "m2");

// Convert to/from SI base units
double value_si = units.toBase(100.0, "mD");
double value_field = units.fromBase(1e-15, "mD");

// Parse value with unit string
double parsed = units.parseAndConvertToBase("5000 psi");

// Check unit compatibility
bool compatible = units.areCompatible("psi", "MPa");  // true
bool incompatible = units.areCompatible("psi", "mD"); // false

// Format output
std::string output = units.formatValue(6.895e6, "MPa", 2);  // "6.90 MPa"

// Convenience functions
double p = toSI(5000, "psi");          // Convert to SI
double v = fromSI(1e-15, "mD");        // Convert from SI
double c = convertUnits(100, "mD", "D"); // Convert between units
```

### Using with ConfigReader

```cpp
#include "ConfigReader.hpp"

ConfigReader config;
config.loadFile("simulation.config");

// Automatic unit conversion from config file
double perm = config.getDoubleWithUnit("ROCK", "permeability_x", 100e-15, "mD");
double visc = config.getDoubleWithUnit("FLUID", "viscosity", 0.001, "cP");
double time = config.getDoubleWithUnit("SIMULATION", "end_time", 86400.0, "day");

// Arrays with units
std::vector<double> coords = config.getDoubleArrayWithUnit("FAULT", "location", "m");
```

## Best Practices

### 1. Always Specify Units in Config Files

✅ **Good**:
```ini
permeability_x = 150 mD
pressure = 5000 psi
time = 30 day
```

❌ **Avoid**:
```ini
permeability_x = 1.48e-13    # What unit is this?
pressure = 5000              # psi? Pa? bar?
```

### 2. Use Field-Standard Units

For petroleum engineering, use commonly recognized units:
- Pressure: `psi`, `bar`
- Permeability: `mD`, `D`
- Viscosity: `cP`
- Rate: `bbl/day`, `Mcf/day`
- Length/Depth: `ft`, `m`

### 3. Use SI for Mechanics

For geomechanics, prefer SI units:
- Stress/Modulus: `MPa`, `GPa`
- Displacement: `m`, `mm`
- Force: `N`, `kN`

### 4. Be Consistent Within Sections

Within each config section, try to use consistent unit systems (all field units or all SI) for readability.

### 5. Document Non-Standard Choices

If using unusual units, add comments:
```ini
thermal_expansion = 1.4e-5 1/degF    # Using Fahrenheit (field standard)
```

## Dimensional Analysis

The unit system performs automatic dimensional analysis to prevent errors:

```cpp
// This works - compatible dimensions
double p1 = units.convert(5000, "psi", "MPa");

// This throws an exception - incompatible dimensions
try {
    double bad = units.convert(100, "mD", "psi");  // Error!
} catch (const std::runtime_error& e) {
    std::cout << "Error: " << e.what() << std::endl;
    // Output: "Incompatible dimensions: L^2 vs L^-1 M T^-2"
}
```

## Common Conversions

### Pressure

| From | To | Factor |
|------|----|----|
| psi | Pa | 6894.757 |
| psi | MPa | 0.006894757 |
| bar | Pa | 100000 |
| atm | Pa | 101325 |

### Permeability

| From | To | Factor |
|------|----|----|
| mD | m² | 9.869233×10⁻¹⁶ |
| D | m² | 9.869233×10⁻¹³ |
| D | mD | 1000 |

### Viscosity

| From | To | Factor |
|------|----|----|
| cP | Pa·s | 0.001 |
| P | Pa·s | 0.1 |
| cP | P | 0.01 |

### Volume

| From | To | Factor |
|------|----|----|
| bbl | m³ | 0.158987295 |
| scf | m³ | 0.028316847 |
| gal | m³ | 0.003785412 |

### Temperature

- 0°C = 32°F = 273.15 K = 491.67°R
- ΔT(°C) = ΔT(K)
- ΔT(°F) = ΔT(°R)
- ΔT(°C) = ΔT(°F) × 5/9

## Extending the Unit System

### Adding Custom Units

```cpp
UnitSystem& units = UnitSystemManager::getInstance();

// Define new unit
Unit my_unit("my_length", "myl", Dimension(1,0,0), 1.5, "length");
units.addUnit(my_unit);

// Add alias
units.addAlias("meter", "metre");  // British spelling

// Use custom unit
double len = units.convert(10, "myl", "m");  // 15 m
```

### Creating Compound Units

For compound units not in the database, use the base units:

```ini
# Thermal diffusivity: m²/s
thermal_diffusivity = 1.5e-6 m2/s

# Stress rate: Pa/s
stress_rate = 1000 Pa/s
```

## Troubleshooting

### Unit Not Found

**Error**: `Unknown unit: millidarcy`

**Solution**: Check spelling. Unit names and symbols are case-insensitive for most units:
- Use `mD`, `milidarcy`, or `millidarcy`
- Use `psi`, `PSI`, or `pounds per square inch`

### Incompatible Dimensions

**Error**: `Incompatible dimensions: L^2 vs L^-1 M T^-2`

**Solution**: You're trying to convert between physically incompatible quantities (e.g., permeability to pressure). Check your config file for typos.

### Parsing Failed

**Error**: `Failed to parse: 5000psi`

**Solution**: Add space between number and unit: `5000 psi`

## Unit Database Reference

To generate a complete list of all available units:

```cpp
UnitSystem& units = UnitSystemManager::getInstance();
units.printDatabase(std::cout);  // Print to console

std::string docs = units.generateDocumentation();  // Generate markdown
```

Or use the command-line tool:

```bash
./fsrm --list-units              # List all units
./fsrm --list-units pressure     # List pressure units only
./fsrm --convert 5000 psi MPa    # Convert values
```

## Performance Considerations

- Unit conversions are lightweight (simple multiplications)
- Conversions happen once during input parsing
- No performance penalty during simulation
- All internal calculations use consistent SI base units

## References

1. International System of Units (SI): [BIPM SI Brochure](https://www.bipm.org/en/publications/si-brochure/)
2. Oil Field Units: SPE Petroleum Engineering Handbook
3. Conversion Factors: NIST Special Publication 811

## Summary

The FSRM unit system provides:

- ✅ Automatic conversion to SI base units (L M T)
- ✅ 150+ supported units across all physical quantities
- ✅ Flexible input in any supported unit
- ✅ User-configurable output units
- ✅ Dimensional analysis and error checking
- ✅ Zero performance overhead during simulation
- ✅ Transparent integration with configuration system

All you need to do is specify units in your config file, and the system handles the rest!
