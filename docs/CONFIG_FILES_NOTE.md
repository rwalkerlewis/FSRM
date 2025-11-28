# Configuration Files for New Benchmarks

## Required Configuration Files

The following configuration files are required for the newly added benchmarks. These follow the same format as existing config files in the `config/` directory.

### New SPE Benchmark Configs (3)

1. **config/spe5_benchmark.config**
   - For: SPE5 - Volatile Oil/Gas Compositional
   - Grid: 7×7×3 (147 cells)
   - Key settings:
     - Compositional fluid model (4 components: C1, C3, C6, C10)
     - 1 injector (center), 4 producers (corners)
     - 1,500 days simulation
     - Timestep: 1-10 days

2. **config/spe11_benchmark.config**
   - For: SPE11 - CO2 Storage CSP
   - Grid: Variable (case A: 840, B: 8,400, C: 168,000 cells)
   - Key settings:
     - Two-phase (water + supercritical CO2)
     - Domain: 2.8×1.2×1.2 km
     - 50 years (25 injection + 25 post-injection)
     - Injection rate: Variable based on case
     - Timestep: 10-100 days

3. **config/spe13_benchmark.config**
   - For: SPE13 - Well Control and Constraints
   - Grid: 24×25×15 (9,000 cells)
   - Key settings:
     - Three-phase black oil
     - 26 wells (mixed producers/injectors)
     - Complex control switching (rate/BHP)
     - 3,000 days simulation
     - Timestep: 1-10 days

### New SCEC Benchmark Configs (5)

4. **config/scec_tpv11.config**
   - For: SCEC TPV11 - Supershear Rupture
   - Grid: 192×192×96 (1.8M cells)
   - Domain: 48×48×24 km
   - Key settings:
     - Vertical strike-slip fault
     - Low strength (enables supershear)
     - Duration: 12 seconds
     - Timestep: 1 ms
     - Material: Homogeneous elastic (Vp=6000 m/s, Vs=3464 m/s)

5. **config/scec_tpv14.config**
   - For: SCEC TPV14 - Bimaterial Fault
   - Grid: 240×192×96 (2.2M cells)
   - Domain: 60×48×24 km
   - Key settings:
     - Vertical strike-slip fault
     - Material contrast across fault
     - Duration: 15 seconds
     - Timestep: 1 ms
     - Side 1: Vp=6000 m/s, Vs=3464 m/s, ρ=2670 kg/m³
     - Side 2: Vp=5196 m/s, Vs=3000 m/s, ρ=2670 kg/m³

6. **config/scec_tpv24.config**
   - For: SCEC TPV24 - Dynamic Triggering
   - Grid: 192×192×96 (1.8M cells)
   - Domain: 48×48×24 km
   - Key settings:
     - Two parallel vertical strike-slip faults
     - Fault 1: Spontaneous rupture
     - Fault 2: Initially stable, triggered by Fault 1
     - Duration: 20 seconds
     - Timestep: 1 ms
     - Fault separation: 15 km

7. **config/scec_loh2.config**
   - For: SCEC LOH.2 - Basin Edge Effects
   - Grid: 200×200×100 (2.0M cells)
   - Domain: 40×40×20 km
   - Key settings:
     - Sedimentary basin with sharp edge
     - Point explosion source at depth
     - Duration: 20 seconds
     - Timestep: 2 ms
     - Surface receivers for monitoring
     - Basin: Vp=2000 m/s, Vs=1000 m/s
     - Halfspace: Vp=4000 m/s, Vs=2000 m/s

8. **config/scec_loh3.config**
   - For: SCEC LOH.3 - Layered Medium
   - Grid: 200×200×100 (2.0M cells)
   - Domain: 40×40×20 km
   - Key settings:
     - 5 horizontal layers with varying properties
     - Point explosion at 10 km depth
     - Duration: 20 seconds
     - Timestep: 2 ms
     - Vertical receiver array
     - Layer properties:
       * 0-2 km: Vp=6000, Vs=3464
       * 2-5 km: Vp=5500, Vs=3175
       * 5-10 km: Vp=5000, Vs=2887
       * 10-15 km: Vp=4500, Vs=2598
       * 15-20 km: Vp=4000, Vs=2309

## Configuration File Format

All config files follow the standard FSRM INI format:

```ini
[SIMULATION]
start_time = 0.0
end_time = <duration>
dt_initial = <timestep>
fluid_model = <BLACK_OIL|COMPOSITIONAL|SINGLE_PHASE>
enable_geomechanics = <true|false>
enable_dynamics = <true|false>

[GRID]
nx = <cells_x>
ny = <cells_y>
nz = <cells_z>
Lx = <domain_x>
Ly = <domain_y>
Lz = <domain_z>

[ROCK]
porosity = <value>
permeability_x = <value>
youngs_modulus = <value>
poissons_ratio = <value>
density = <value>

[FLUID]
# Specific to fluid model
...

[WELL1]
# Well definitions
...

[FAULT1]  # For SCEC benchmarks
# Fault properties
...
```

## How to Create Config Files

### Option 1: Generate Template
```bash
fsrm -generate_config my_benchmark.config
```

### Option 2: Copy and Modify
```bash
cp config/spe1_benchmark.config config/spe5_benchmark.config
# Edit the new file with appropriate parameters
```

### Option 3: Use Existing Complete Template
```bash
cp config/complete_template.config config/spe5_benchmark.config
# Uncomment and set relevant sections
```

## Notes

- All executables will fall back to reasonable defaults if config file is missing
- Config files enable full customization without recompilation
- See `docs/CONFIGURATION.md` for complete parameter reference
- See existing config files in `config/` directory for examples

## Status

✅ Executables implemented and will load these config files  
⏳ Config files can be created by user based on requirements  
✅ All parameters documented in existing configs and docs  
✅ Template generation tool available
