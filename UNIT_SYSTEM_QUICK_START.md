# Unit System Quick Start Guide

## 5-Minute Quick Start

### What is it?

The FSRM unit system lets you use **any units** in your config files. The system automatically converts everything to SI base units (meters, kilograms, seconds) for calculations.

### Basic Usage

Just add units after numbers in your config file:

```ini
[ROCK]
permeability_x = 150 mD              # Instead of: 1.48e-13
youngs_modulus = 15 GPa              # Instead of: 15e9
density = 2.55 g/cm3                 # Instead of: 2550

[FLUID]
viscosity = 5 cP                     # Instead of: 0.005
density = 850 kg/m3                  # Already SI, but clearer

[WELL1]
target_value = 5000 bbl/day          # Instead of: 9.2e-3
min_bhp = 2000 psi                   # Instead of: 1.38e7

[SIMULATION]
end_time = 30 day                    # Instead of: 2592000
dt_initial = 1 hr                    # Instead of: 3600
```

### Supported Units

**Most Common:**
- Pressure: `Pa`, `kPa`, `MPa`, `GPa`, `psi`, `bar`, `atm`
- Permeability: `m2`, `D`, `mD`, `uD`
- Viscosity: `Pa-s`, `cP`, `P`
- Length: `m`, `cm`, `mm`, `km`, `ft`, `in`
- Volume: `m3`, `L`, `bbl`, `ft3`, `gal`, `scf`, `Mcf`
- Rate: `m3/s`, `m3/day`, `bbl/day`, `stb/day`, `Mcf/day`
- Time: `s`, `min`, `hr`, `day`, `year`
- Temperature: `K`, `degC`, `degF`, `R`
- Density: `kg/m3`, `g/cm3`, `lbm/ft3`

**See full list:** `docs/UNIT_SYSTEM.md`

### Examples

#### Petroleum Engineering
```ini
[ROCK]
permeability_x = 150 mD
porosity = 0.22

[FLUID]
oil_viscosity = 5 cP
water_density = 1020 kg/m3

[WELL1]
type = PRODUCER
target_value = 2000 psi
max_rate = 5000 bbl/day
```

#### Geomechanics
```ini
[ROCK]
youngs_modulus = 15 GPa
density = 2.55 g/cm3

[FAULT1]
cohesion = 2 MPa
strike = 45 deg
dip = 70 deg
```

#### Mixed Units
```ini
[GRID]
Lx = 1 km
Ly = 3280 ft
Lz = 100 m

[SIMULATION]
end_time = 30 day
dt_initial = 1 hr
```

### Output Units

Control how results are displayed:

```ini
[OUTPUT]
format = VTK
frequency = 10

# Display units
pressure_unit = psi
displacement_unit = mm
stress_unit = MPa
permeability_unit = mD
temperature_unit = degC
```

### Testing

Try the example config:
```bash
# Copy the example
cp config/with_units_example.config config/my_simulation.config

# Edit it with your values
nano config/my_simulation.config

# Run simulation
./fsrm config/my_simulation.config
```

### Rules

1. ✅ **DO** add space between number and unit: `5000 psi` ✓
2. ❌ **DON'T** forget the space: `5000psi` ✗
3. ✅ **DO** use standard abbreviations: `mD`, `cP`, `GPa`
4. ❌ **DON'T** make up units: `milliD`, `centiP`
5. ✅ **DO** specify units when possible
6. ❌ **DON'T** mix values with and without units

### Common Conversions

| From | To | Multiply by |
|------|----|----|
| psi | Pa | 6894.757 |
| mD | m² | 9.869×10⁻¹⁶ |
| cP | Pa·s | 0.001 |
| bbl | m³ | 0.15899 |
| ft | m | 0.3048 |
| bbl/day | m³/s | 1.840×10⁻⁶ |

*Don't worry about these - the system does it automatically!*

### Need Help?

- **Full docs**: `docs/UNIT_SYSTEM.md`
- **Example config**: `config/with_units_example.config`
- **Test program**: `./test_unit_system`
- **Implementation details**: `UNIT_SYSTEM_IMPLEMENTATION.md`

### That's it!

Just add units to your config files and everything works automatically. 🎉

---

**Pro tip**: Start with the example config (`with_units_example.config`) and modify it for your case.
