# Tsunami Modeling and Ocean Coupling

## Overview

FSRM includes comprehensive tsunami modeling capabilities that couple solid-earth earthquake rupture with oceanic wave propagation. This enables end-to-end simulation of earthquake-generated tsunamis from source to coastal inundation.

## Key Features

### 1. Earthquake-Tsunami Coupling

- **Coseismic Seafloor Deformation**: Okada (1985) elastic dislocation model
- **Kinematic Rupture**: Prescribed slip with variable rise time
- **Dynamic Rupture**: Coupling with `SeismicFaultModel` for spontaneous rupture
- **Kajiura Filter**: Proper treatment of water column response

### 2. Tsunami Propagation

- **Nonlinear Shallow Water Equations (NSWE)**:
  ```
  ∂h/∂t + ∇·(hu) = 0
  ∂(hu)/∂t + ∇·(huu + ½gh²I) = -gh∇b - τ_b/ρ
  ```
- **Numerical Solvers**:
  - F-wave method (GeoClaw-style)
  - HLLC Riemann solver
  - Roe approximate Riemann solver
  - Godunov exact Riemann solver
- **High-Order Reconstruction**: 2nd order with minmod limiter
- **Positivity Preserving**: Ensures non-negative water depth

### 3. Coastal Inundation

- **Wetting/Drying**: Momentum-conserving algorithm
- **Bottom Friction**: Manning's n roughness model
- **Runup Calculation**: Automatic tracking of maximum runup elevation
- **Land Cover Effects**: Variable friction by terrain type

### 4. Bathymetry Handling

- **GEBCO**: Global 15 arc-second bathymetry
- **ETOPO**: 1 arc-minute global relief
- **Custom Grids**: NetCDF, ASCII, GeoTIFF formats
- **Nested Grids**: Multi-resolution for coastal detail

## Cascadia Subduction Zone Example

### Background

The Cascadia Subduction Zone extends ~1,100 km from Cape Mendocino, California to Vancouver Island, British Columbia. The last great earthquake occurred on January 26, 1700 CE, producing:

- Estimated magnitude: Mw 9.0
- Rupture length: ~1,000 km
- Average slip: ~17 m (peak ~35 m)
- Trans-Pacific tsunami observed in Japan

### Running the Example

```bash
# Full margin M9.0 scenario (4 hours simulation)
mpirun -np 8 ./cascadia_tsunami -scenario full_margin -end_time 14400

# Southern segment rupture (M8.0-8.5)
mpirun -np 4 ./cascadia_tsunami -scenario southern

# Worst-case planning scenario
mpirun -np 16 ./cascadia_tsunami -scenario worst_case -end_time 28800
```

### Available Scenarios

| Scenario | Magnitude | Extent | Average Slip |
|----------|-----------|--------|--------------|
| `full_margin` | Mw 9.0+ | Cape Mendocino to Vancouver Island | 17 m |
| `southern` | Mw 8.0-8.5 | Cape Mendocino to central Oregon | 8 m |
| `northern` | Mw 8.0-8.5 | Washington to Vancouver Island | 10 m |
| `central` | Mw 8.0 | Central Oregon | 12 m |
| `worst_case` | Mw 9.2 | Full margin with maximum slip | 22 m |

### Fault Model Parameters

The Cascadia fault model is based on:

- **USGS SLAB2.0**: Subduction geometry
- **Wang et al. (2013)**: Heterogeneous slip models
- **Goldfinger et al. (2012)**: Turbidite paleoseismology
- **Witter et al. (2013)**: Paleotsunami constraints

```cpp
CascadiaFaultModel model;
model.north_latitude = 50.0;      // Vancouver Island
model.south_latitude = 40.5;      // Cape Mendocino
model.rupture_width = 120.0;      // km
model.average_slip = 17.0;        // m
model.peak_slip = 35.0;           // m
model.rupture_velocity = 2.8;     // km/s
```

### Expected Tsunami Arrivals

| Location | First Arrival | Max Height |
|----------|--------------|------------|
| DART Buoys | 15-25 min | 0.3-0.5 m |
| Crescent City | 20-30 min | 5-10 m |
| Oregon Coast | 20-35 min | 5-10 m |
| Westport, WA | 25-35 min | 6-10 m |
| Seattle | 2-3 hours | 1-3 m |
| San Francisco | 1.5-2 hours | 1-3 m |
| Los Angeles | 3-4 hours | 0.5-2 m |

## Configuration File Format

### Basic Tsunami Configuration

```ini
[TSUNAMI]
solver = FWAVE
spatial_order = 2
time_integration = SSP_RK3
cfl_number = 0.8
adaptive_timestep = true
source_type = KINEMATIC
source_rise_time = 60.0
boundary_west = OPEN
boundary_east = REFLECTIVE

[INUNDATION]
enabled = true
method = MOMENTUM_CONSERVING
dry_tolerance_m = 0.01
manning_n_ocean = 0.025
manning_n_land = 0.03
```

### Fault Model Configuration

```ini
[FAULT_MODEL]
north_latitude = 50.0
south_latitude = 40.5
trench_longitude = -125.0
average_strike = 350.0
shallow_dip = 8.0
deep_dip = 18.0
rupture_width_km = 120.0
num_subfaults_along_strike = 60
num_subfaults_down_dip = 12

[SLIP_MODEL]
type = HETEROGENEOUS
average_slip = 17.0
peak_slip = 35.0
rake = 90.0
```

### Observation Stations

```ini
[GAUGE_CRESCENT_CITY]
name = Crescent_City
type = TIDE_GAUGE
station_id = 9419750
longitude = -124.183
latitude = 41.745

[GAUGE_DART_46404]
name = DART_46404
type = DART
longitude = -128.894
latitude = 45.857
depth = 2785
```

## API Reference

### Core Classes

#### `CascadiaFaultModel`
```cpp
CascadiaFaultModel model;
model.generateSubfaults();     // Returns vector<TsunamiSubfault>
model.getTotalMoment();        // Returns seismic moment (N·m)
model.getMagnitude();          // Returns Mw
```

#### `OkadaModel`
```cpp
OkadaModel okada;
okada.computeDisplacement(lon, lat, subfault, ux, uy, uz);
okada.computeDisplacementField(grid, subfaults, uz_field);
okada.applyKajiuraFilter(grid, uz_field);
```

#### `ShallowWaterSolver`
```cpp
ShallowWaterSolver swe;
swe.initialize(bathymetry, config);
swe.setInitialCondition(eta_initial);
swe.setSeafloorMotion(motion_func);
swe.run(end_time);
swe.getMaximumValues(max_eta, max_speed, max_mflux);
```

#### `CoupledTsunamiEarthquake`
```cpp
CoupledTsunamiEarthquake sim;
sim.initializeCascadia(fault_model, bathymetry, config);
sim.run();
auto arrivals = sim.getArrivalTimes();
auto amplitudes = sim.getMaxAmplitudes();
sim.writeOutput(output_dir);
```

### Utility Functions

```cpp
// Tsunami physics
double speed = TsunamiUtils::waveSpeed(depth);
double arrival = TsunamiUtils::estimateArrivalTime(lon1, lat1, lon2, lat2, depth);
double amplified = TsunamiUtils::greensLaw(eta_deep, h_deep, h_shallow);
double runup = TsunamiUtils::estimateInundation(wave_height, slope);

// Bathymetry
BathymetryGrid bathy = TsunamiUtils::generateTestBathymetry(
    lon_min, lon_max, lat_min, lat_max, resolution);
```

### Pre-defined Gauge Locations

```cpp
// Get all West Coast stations
auto gauges = WestCoastGaugeNetwork::getAllStations();

// Specific locations
auto cc = WestCoastGaugeNetwork::crescent_city();
auto dart = WestCoastGaugeNetwork::dart_46404();
```

## Output Files

| File | Description |
|------|-------------|
| `*_timeseries.dat` | Sea surface elevation vs time at gauges |
| `maximum_values.dat` | Max η, velocity, momentum flux grid |
| `inundation_extent.dat` | Areas that were flooded |
| `arrival_time.dat` | First arrival time at each grid point |
| `seafloor_displacement.dat` | Coseismic vertical deformation |

## Physical Validation

The tsunami model has been validated against:

1. **Analytical Solutions**
   - Carrier & Greenspan (1958) runup on slopes
   - Synolakis (1987) solitary wave runup

2. **Laboratory Experiments**
   - Solitary wave over a conical island (Briggs et al., 1995)
   - Tsunami runup on a plane beach (Liu et al., 1991)

3. **Historical Tsunamis**
   - 1700 Cascadia (Japanese historical records)
   - 1964 Alaska (tide gauge records)
   - 2011 Tohoku (DART and tide gauge records)

## References

### Cascadia Seismotectonics
- Goldfinger, C., et al. (2012). Turbidite event history. USGS Professional Paper 1661-F.
- Wang, K., et al. (2013). Earthquake Science, 26(1), 1-7.
- Witter, R.C., et al. (2013). Natural Hazards and Earth System Sciences, 13, 1847-1869.

### Tsunami Modeling
- Okada, Y. (1985). Surface deformation due to shear and tensile faults. Bull. Seismol. Soc. Am., 75, 1135-1154.
- Kajiura, K. (1963). The leading wave of a tsunami. Bull. Earthquake Res. Inst., 41, 535-571.
- LeVeque, R.J., et al. (2011). GeoClaw tsunami modeling. Adv. Water Resources, 34, 1195-1206.

### Coastal Hazards
- Synolakis, C.E. (1987). The runup of solitary waves. J. Fluid Mech., 185, 523-545.
- NOAA Center for Tsunami Research: https://nctr.pmel.noaa.gov/

## Disclaimer

This software is for research and educational purposes. For official tsunami hazard assessments and evacuation planning, consult:

- **NOAA National Tsunami Warning Center**
- **FEMA Region X (Cascadia Rising Exercise)**
- **State/Provincial Emergency Management Agencies**
