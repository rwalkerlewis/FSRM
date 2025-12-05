# Fully Coupled Volcano Modeling

## Overview

FSRM includes comprehensive physics for simulating volcanic systems from deep magma chambers to surface eruptions and their atmospheric/surface hazards. The volcano modeling module provides a fully coupled multi-physics framework suitable for:

- **Volcanic unrest analysis**: Pre-eruption signals, deformation, seismicity, degassing
- **Eruption dynamics**: Conduit flow, fragmentation, column dynamics
- **Hazard assessment**: PDCs, lava flows, lahars, tephra fallout, gas dispersal
- **Volcano monitoring**: Synthetic observables for GPS, InSAR, seismometers, gas sensors

## Physical Components

### 1. Magma Chamber Dynamics

The magma chamber model tracks:

- **Pressure evolution**: From recharge, cooling, crystallization, volatile exsolution
- **Temperature evolution**: Cooling through wall rocks, latent heat from crystallization
- **Volatile content**: Dissolved vs. exsolved H₂O, CO₂, SO₂
- **Crystal fraction**: Temperature-dependent crystallization

```cpp
MagmaChamberGeometry chamber;
chamber.depth = 5000.0;      // m below surface
chamber.volume = 10e9;       // m³ (10 km³)
chamber.semi_axis_a = 2000.0;  // m
chamber.semi_axis_b = 2000.0;  // m
chamber.semi_axis_c = 1000.0;  // m (oblate spheroid)

MagmaProperties magma;
magma.composition = MagmaComposition::ANDESITE;
magma.SiO2 = 60.0;          // wt%
magma.H2O_total = 4.0;      // wt%
magma.temperature = 1173.15; // K (900°C)
magma.crystal_fraction = 0.3;

MagmaChamberModel chamber_model;
chamber_model.initialize(chamber, magma, initial_state);
chamber_model.setRechargeRate(100.0);  // kg/s
chamber_model.update(dt);
```

Key equations:

```
dP/dt = (1/βV) * [ṁ_recharge - ṁ_erupt + ρ_melt * V * (dφ_gas/dt)]
dT/dt = (1/ρc_p V) * [Q_recharge - Q_wall + L * dφ_cryst/dt]
```

### 2. Conduit Flow Model

Models magma ascent from chamber to vent:

- **Bubble nucleation and growth**: Volatile exsolution
- **Flow regimes**: Bubbly → Slug → Annular → Dispersed
- **Fragmentation**: Brittle failure when gas fraction > 75%
- **Wall friction**: Temperature-dependent viscosity

```cpp
ConduitGeometry conduit;
conduit.depth = 5000.0;       // m
conduit.radius_base = 15.0;   // m (at chamber)
conduit.radius_vent = 10.0;   // m (at surface)

ConduitFlowModel flow;
flow.initialize(conduit, magma, chamber_pressure);
flow.solvesteadyState();

double MER = flow.getMassFlux();           // kg/s
double exit_vel = flow.getExitVelocity();  // m/s
double frag_depth = flow.getFragmentationDepth();  // m
```

Governing equations:

```
Mass:      d(ρuA)/dz = 0
Momentum:  ρu(du/dz) = -dP/dz - ρg - τ_wall/r
Energy:    d(ρuH·A)/dz = ...
```

### 3. Eruption Column Model

Models buoyant volcanic plumes:

- **Jet thrust phase**: Initial momentum-driven rise
- **Convective phase**: Buoyancy-driven with air entrainment
- **Column collapse**: When density exceeds ambient
- **Neutral buoyancy level**: Umbrella cloud height

```cpp
EruptionColumnParameters params;
params.mass_flux = 1e7;          // kg/s
params.exit_velocity = 200.0;    // m/s
params.exit_temperature = 1273.0; // K
params.gas_fraction = 0.95;

EruptionColumnModel column;
column.initialize(params);
bool success = column.solve();

double max_height = column.getMaxHeight();       // m
double NBL = column.getNeutralBuoyancyHeight(); // m
bool collapse = column.doesCollapse();
```

Based on Woods (1988) model:

```
d(ρuA)/dz = ρ_a * α * u * 2πr     (entrainment)
d(ρu²A)/dz = (ρ_a - ρ)gA          (momentum)
```

### 4. Pyroclastic Density Currents (PDCs)

Models hot, fast-moving volcanic flows:

- **Column collapse PDCs**: From eruption column failure
- **Dome collapse PDCs**: From lava dome failure
- **Directed blasts**: Lateral explosions

```cpp
PDCSourceConditions source;
source.type = PDCType::COLUMN_COLLAPSE;
source.mass_flux = 1e7;           // kg/s
source.initial_velocity = 50.0;   // m/s
source.initial_temperature = 700.0; // K

PDCModel pdc;
pdc.initialize(source);
pdc.setTopography(elevation_func);
pdc.runToCompletion();

double runout = pdc.getRunoutDistance();
double dyn_pressure = pdc.getDynamicPressure(x, y);
```

Key outputs:
- Runout distance
- Dynamic pressure (building damage)
- Temperature (thermal hazard)
- Deposit thickness

### 5. Lava Flow Model

Models effusive lava flow propagation:

- **Temperature-dependent rheology**: Viscosity vs. T
- **Crust formation**: Insulation effects
- **Topographic control**: Channeling and ponding
- **Yield strength**: For crystal-rich lavas

```cpp
LavaSourceParameters source;
source.effusion_rate = 10.0;      // m³/s DRE
source.temperature = 1373.15;     // K (1100°C)
source.vent_x = 0.0;
source.vent_y = 0.0;

LavaFlowModel lava;
lava.initialize(source, magma);
lava.setTopography(elevation_func);
lava.run(86400.0);  // 24 hours

double flow_area = lava.getFlowArea();
double flow_length = lava.getFlowLength();
```

### 6. Lahar (Volcanic Debris Flow) Model

Models volcanic mudflows:

- **Primary lahars**: Eruption-triggered
- **Secondary lahars**: Rain-remobilized deposits
- **Jökulhlaups**: Glacier outburst floods

```cpp
LaharSourceConditions source;
source.type = LaharType::PRIMARY_COLD;
source.volume = 1e7;              // m³
source.peak_discharge = 1e4;      // m³/s
source.sediment_concentration = 0.5;

LaharModel lahar;
lahar.initialize(source);
lahar.setTopography(elevation_func);
lahar.runToCompletion();

double runout = lahar.getRunoutDistance();
double travel_time = lahar.getTravelTimeTo(x, y);
```

### 7. Volcanic Deformation

Computes surface deformation from subsurface sources:

- **Mogi source**: Point spherical (Mogi, 1958)
- **McTigue**: Finite spherical cavity
- **Okada dike/sill**: Rectangular opening
- **Spheroid**: Prolate/oblate chambers

```cpp
VolcanicDeformationSource source;
source.type = DeformationSourceType::MOGI;
source.depth = 5000.0;            // m
source.volume_change = 1e6;       // m³

VolcanicDeformationModel def;
def.addSource(source);
def.setElasticProperties(30e9, 0.25);

double ux, uy, uz;
def.computeDisplacement(x, y, ux, uy, uz);

double tilt_x, tilt_y;
def.computeTilt(x, y, tilt_x, tilt_y);

// InSAR line-of-sight
double LOS = def.computeLOS_displacement(x, y, look_az, incidence);
```

### 8. Volcanic Seismicity

Models volcanic earthquake sources:

- **VT (Volcano-Tectonic)**: Brittle fracture
- **LP (Long-Period)**: Fluid-filled crack resonance
- **VLP (Very Long-Period)**: Conduit resonance
- **Tremor**: Continuous oscillation
- **Explosion quakes**: From explosions

```cpp
VolcanicSeismicEvent event;
event.type = VolcanicSeismicEventType::LP;
event.x = 0.0; event.y = 0.0; event.z = -2000.0;
event.dominant_frequency = 1.5;   // Hz
event.quality_factor = 20.0;      // Q
event.pressure_transient = 1e5;   // Pa

VolcanicSeismicityModel seismicity;
std::vector<double> times;
std::vector<std::array<double, 6>> moment_tensor_rate;
seismicity.generateLP_source(event, times, moment_tensor_rate);
```

### 9. Volcanic Gas Model

Models gas emission and dispersal:

- **Source terms**: Vents, fumaroles, ground
- **Gaussian plume dispersion**: Wind-driven
- **Species**: H₂O, CO₂, SO₂, H₂S, HCl, HF

```cpp
GasEmissionSource source;
source.emission_rate = 1000.0;    // kg/s total
source.temperature = 373.15;      // K
source.composition = {
    {VolcanicGas::H2O, 0.90},
    {VolcanicGas::CO2, 0.05},
    {VolcanicGas::SO2, 0.03}
};

VolcanicGasModel gas;
gas.addSource(source);
gas.setWindSpeed(5.0);
gas.setWindDirection(270.0);
gas.setAtmosphericStability('D');

double SO2_conc = gas.getConcentrationPPM(VolcanicGas::SO2, x, y);
```

### 10. Tephra Dispersal Model

Models volcanic ash dispersal and fallout:

- **3D advection**: By wind field
- **Gravitational settling**: Size-dependent
- **Ground deposition**: Accumulation
- **Particle aggregation**: Optional

```cpp
TephraSourceParameters source;
source.column_height = 15000.0;   // m
source.mass_eruption_rate = 1e7;  // kg/s
source.phi_mean = 2.0;            // phi = -log2(d_mm)
source.phi_stddev = 2.0;

TephraDispersalModel tephra;
tephra.initialize(source);
tephra.setWindField(wind_func);
tephra.run(3600.0);  // 1 hour

double thickness = tephra.getDepositThickness(x, y);  // m
double loading = tephra.getDepositLoading(x, y);      // kg/m²
```

## Coupled Volcanic System

The `CoupledVolcanoSystem` class integrates all components:

```cpp
VolcanoSystemConfig config;
config.enable_magma_chamber = true;
config.enable_conduit_flow = true;
config.enable_eruption_column = true;
config.enable_pdc = true;
config.enable_deformation = true;
config.enable_seismicity = true;
config.enable_gas_emission = true;
config.enable_tephra = true;

CoupledVolcanoSystem volcano;
volcano.initialize(config, chamber_geometry, magma, conduit);
volcano.setTopography(elevation_func);
volcano.setRechargeRate(100.0);
volcano.triggerEruption(EruptionType::SUBPLINIAN);
volcano.run();

// Get monitoring time series
auto monitoring = volcano.getMonitoringTimeSeries();
for (const auto& m : monitoring) {
    std::cout << "t=" << m.time << "s: "
              << "P=" << m.chamber_pressure/1e6 << "MPa, "
              << "uz=" << m.vertical_displacement*1000 << "mm, "
              << "MER=" << m.mass_eruption_rate << "kg/s\n";
}
```

## Pre-defined Scenarios

### Mount St. Helens 1980 (Lateral Blast)

```cpp
auto config = VolcanoScenarios::mountStHelens1980();
auto chamber = VolcanoScenarios::mountStHelens1980_chamber();
auto magma = VolcanoScenarios::mountStHelens1980_magma();

CoupledVolcanoSystem volcano;
volcano.initialize(config, chamber, magma, conduit);
volcano.triggerEruption(EruptionType::LATERAL_BLAST);
```

### Pinatubo 1991 (Plinian)

```cpp
auto config = VolcanoScenarios::pinatubo1991();
auto chamber = VolcanoScenarios::pinatubo1991_chamber();
auto magma = VolcanoScenarios::pinatubo1991_magma();
// VEI-6 Plinian eruption
```

### Kilauea (Effusive)

```cpp
auto config = VolcanoScenarios::kilauea_effusive();
auto chamber = VolcanoScenarios::kilauea_chamber();
auto magma = VolcanoScenarios::kilauea_magma();
// Hawaiian-style lava flow
```

### Yellowstone Unrest

```cpp
auto config = VolcanoScenarios::yellowstone_unrest();
auto chamber = VolcanoScenarios::yellowstone_chamber();
auto magma = VolcanoScenarios::yellowstone_magma();
// Caldera unrest scenario (no eruption)
```

## Configuration File Format

### Basic Volcano Configuration

```ini
[VOLCANO]
enable_chamber = true
enable_conduit = true
enable_column = true
enable_pdc = true
enable_deformation = true
enable_seismicity = true
enable_gas = true

[MAGMA_CHAMBER]
depth_m = 5000.0
volume_km3 = 10.0
semi_axis_a_m = 2000.0
semi_axis_b_m = 2000.0
semi_axis_c_m = 1000.0
initial_pressure_MPa = 200.0
initial_temperature_C = 900.0
recharge_rate_kg_s = 100.0

[MAGMA]
composition = ANDESITE
SiO2_wt = 60.0
H2O_total_wt = 4.0
CO2_total_wt = 0.1
temperature_C = 1000.0
crystal_fraction = 0.3

[CONDUIT]
length_m = 5000.0
radius_base_m = 15.0
radius_vent_m = 10.0

[ERUPTION]
type = SUBPLINIAN
trigger_overpressure_MPa = 20.0
```

### Deformation Source Configuration

```ini
[DEFORMATION_SOURCE]
type = MOGI
x_m = 0.0
y_m = 0.0
depth_m = 5000.0
volume_change_m3 = 1e6
```

### Gas Source Configuration

```ini
[GAS_SOURCE]
x_m = 0.0
y_m = 0.0
emission_rate_kg_s = 1000.0
temperature_K = 373.15
H2O_fraction = 0.90
CO2_fraction = 0.05
SO2_fraction = 0.03
```

## API Reference

### Core Classes

| Class | Description |
|-------|-------------|
| `MagmaChamberModel` | Magma chamber thermodynamics |
| `ConduitFlowModel` | Conduit flow and fragmentation |
| `EruptionColumnModel` | Buoyant column dynamics |
| `PDCModel` | Pyroclastic density currents |
| `LavaFlowModel` | Lava flow propagation |
| `LaharModel` | Volcanic debris flows |
| `VolcanicDeformationModel` | Surface deformation |
| `VolcanicSeismicityModel` | Volcanic earthquake sources |
| `VolcanicGasModel` | Gas emission and dispersal |
| `TephraDispersalModel` | Ash fall modeling |
| `CoupledVolcanoSystem` | Integrated multi-physics |

### Enumerations

| Enum | Values |
|------|--------|
| `EruptionType` | EFFUSIVE, HAWAIIAN, STROMBOLIAN, VULCANIAN, SUBPLINIAN, PLINIAN, ULTRAPLINIAN, SURTSEYAN, DOME_COLLAPSE, LATERAL_BLAST |
| `MagmaComposition` | BASALT, ANDESITE, DACITE, RHYOLITE, ... |
| `VolcanicSeismicEventType` | VT, LP, VLP, TREMOR, HYBRID, EXPLOSION_QUAKE |
| `DeformationSourceType` | MOGI, MCTIGUE, PENNY_CRACK, PROLATE_SPHEROID, RECTANGULAR_DIKE, ... |
| `PDCType` | COLUMN_COLLAPSE, DOME_COLLAPSE, DIRECTED_BLAST, SURGE |
| `LaharType` | PRIMARY_HOT, PRIMARY_COLD, SECONDARY, JÖKULHLAUP |

## Volcanic Explosivity Index (VEI) Reference

| VEI | Volume (km³) | Column Height (km) | MER (kg/s) | Example |
|-----|--------------|--------------------|-----------:|---------|
| 0 | < 10⁻⁶ | < 0.1 | < 10³ | Continuous degassing |
| 1 | 10⁻⁶ - 10⁻⁵ | 0.1 - 1 | 10³ - 10⁴ | Stromboli |
| 2 | 10⁻⁴ - 10⁻³ | 1 - 5 | 10⁴ - 10⁶ | Galeras 1993 |
| 3 | 10⁻³ - 10⁻² | 3 - 15 | 10⁵ - 10⁷ | Ruiz 1985 |
| 4 | 10⁻² - 10⁻¹ | 10 - 25 | 10⁶ - 10⁸ | Eyjafjallajökull 2010 |
| 5 | 10⁻¹ - 1 | > 25 | 10⁷ - 10⁸ | St. Helens 1980 |
| 6 | 1 - 10 | > 25 | > 10⁸ | Pinatubo 1991 |
| 7 | 10 - 100 | > 25 | > 10⁹ | Tambora 1815 |
| 8 | > 100 | > 25 | > 10¹⁰ | Yellowstone 640 ka |

## Physical Validation

The volcano model has been validated against:

### 1. Analytical Solutions
- Mogi (1958) point source deformation
- Woods (1988) eruption column height
- Sparks (1978) column collapse criteria

### 2. Historical Eruptions
- **Mt. St. Helens 1980**: Lateral blast dynamics, PDC runout
- **Pinatubo 1991**: Column height, SO₂ emission, tephra fallout
- **Kilauea 2018**: Lava flow paths and velocities
- **Ruapehu 1995**: Lahar travel times

### 3. Benchmark Problems
- IAVCEI column collapse benchmark
- Lava flow benchmark (Harris & Rowland, 2001)
- Tephra dispersal benchmark (Folch et al., 2010)

## Scientific References

### Magma Dynamics
- Sparks, R.S.J. (1978). The dynamics of bubble formation and growth in magmas. *J. Volcanol. Geotherm. Res.*, 3, 1-37.
- Gonnermann, H.M. & Manga, M. (2007). The fluid mechanics inside a volcano. *Annu. Rev. Fluid Mech.*, 39, 321-356.

### Eruption Columns
- Woods, A.W. (1988). The fluid dynamics and thermodynamics of eruption columns. *Bull. Volcanol.*, 50, 169-193.
- Bursik, M. (2001). Effect of wind on the rise height of volcanic plumes. *Geophys. Res. Lett.*, 28, 3621-3624.

### PDCs
- Denlinger, R.P. & Iverson, R.M. (2001). Flow of variably fluidized granular masses across three-dimensional terrain. *J. Geophys. Res.*, 106, 553-566.

### Deformation
- Mogi, K. (1958). Relations between the eruptions of various volcanoes and the deformation of the ground surfaces around them. *Bull. Earthq. Res. Inst.*, 36, 99-134.
- McTigue, D.F. (1987). Elastic stress and deformation near a finite spherical magma body. *J. Geophys. Res.*, 92, 12931-12940.

### Volcanic Seismicity
- Chouet, B. (1996). Long-period volcano seismicity: its source and use in eruption forecasting. *Nature*, 380, 309-316.

### Tephra Dispersal
- Mastin, L.G., et al. (2009). A multidisciplinary effort to assign realistic source parameters to models of volcanic ash-cloud transport and dispersion during eruptions. *J. Volcanol. Geotherm. Res.*, 186, 10-21.

## Disclaimer

This software is for research and educational purposes. For official volcanic hazard assessments, consult:

- **USGS Volcano Hazards Program**: https://volcanoes.usgs.gov/
- **Global Volcanism Program**: https://volcano.si.edu/
- **IAVCEI (International Association of Volcanology)**
- **Local volcano observatories**
