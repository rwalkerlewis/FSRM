# Coordinate System Support

FSRM supports geographic and projected coordinate systems using EPSG codes and the PROJ library, enabling simulations with real-world coordinates.

## Overview

Reservoir simulations often use data from GIS systems, well surveys, and seismic surveys that reference specific coordinate systems. FSRM can:

- Accept input in any EPSG-defined coordinate reference system (CRS)
- Transform coordinates to a local model CRS for simulation
- Apply local origin offsets to reduce numerical precision issues
- Auto-detect appropriate UTM zones from geographic coordinates

## Prerequisites

For full coordinate transformation support, build FSRM with PROJ:

```bash
# Install PROJ library
sudo apt install libproj-dev   # Debian/Ubuntu
brew install proj              # macOS

# Build with PROJ support
cmake .. -DENABLE_PROJ=ON
make
```

Without PROJ, only local coordinate offsets are available.

## Configuration

### Basic CRS Transformation

Transform from WGS84 (GPS coordinates) to UTM for simulation:

```ini
[GRID]
# Input data is in WGS84 (lat/lon)
input_crs = EPSG:4326

# Simulate in UTM Zone 10N (meters)
model_crs = EPSG:32610
```

### Local Coordinates

Center the model at a local origin to improve numerical precision:

```ini
[GRID]
input_crs = EPSG:4326
model_crs = EPSG:32610

use_local_coordinates = true
local_origin_x = -122.4194    # Project center longitude
local_origin_y = 37.7749      # Project center latitude
local_origin_z = -2000.0      # Depth reference (e.g., mean sea level)
```

### Auto-detect UTM Zone

Let FSRM automatically select the appropriate UTM zone:

```ini
[GRID]
input_crs = EPSG:4326
auto_detect_utm = true

use_local_coordinates = true
local_origin_x = -122.4194
local_origin_y = 37.7749
```

## Common CRS Codes

### Geographic (Latitude/Longitude)

| EPSG Code | Name | Usage |
|-----------|------|-------|
| 4326 | WGS 84 | GPS, global standard |
| 4269 | NAD83 | North America modern |
| 4267 | NAD27 | North America legacy |
| 4258 | ETRS89 | Europe |

### UTM Zones (WGS84)

Northern Hemisphere: EPSG:326XX (XX = zone number)
Southern Hemisphere: EPSG:327XX

| Region | UTM Zone | EPSG Code |
|--------|----------|-----------|
| California | 10N | 32610 |
| Texas (West) | 13N | 32613 |
| Texas (Central) | 14N | 32614 |
| Louisiana | 15N | 32615 |
| North Sea | 31N | 32631 |
| Australia (East) | 56S | 32756 |

### State Plane (NAD83)

| EPSG Code | Name |
|-----------|------|
| 2227 | California Zone 3 |
| 2277 | Texas Central |
| 2868 | Louisiana South |

### Web Mapping

| EPSG Code | Name | Usage |
|-----------|------|-------|
| 3857 | Web Mercator | Google Maps, web |

## API Usage

### Basic Transformation

```cpp
#include "CoordinateSystem.hpp"

using namespace FSRM;

// Create transformer
CoordinateTransformer transformer;
transformer.setSourceCRS("EPSG:4326");  // WGS84
transformer.setTargetCRS("EPSG:32610"); // UTM 10N
transformer.initialize();

// Transform a point
GeoPoint wgs84(-122.4194, 37.7749, 0);  // San Francisco
GeoPoint utm = transformer.transform(wgs84);

std::cout << "UTM: E=" << utm.x << " N=" << utm.y << "\n";
```

### Using CoordinateSystemManager

```cpp
CoordinateSystemManager manager;
manager.setInputCRS("EPSG:4326");
manager.setModelCRS("EPSG:32610");

// Set local origin
GeoPoint origin(-122.4194, 37.7749, 0);
manager.setLocalOrigin(origin);

// Transform well location
GeoPoint well_input(-122.42, 37.78, -1500);
GeoPoint well_model = manager.toModelCoords(well_input);

// Now well_model has coordinates relative to local origin in meters
```

### Geodetic Calculations

```cpp
using namespace FSRM::Geodetic;

// Distance between two points (Vincenty method - ellipsoid)
double dist = vincentyDistance(37.7749, -122.4194,   // San Francisco
                               34.0522, -118.2437);   // Los Angeles
std::cout << "Distance: " << dist / 1000 << " km\n";

// Bearing calculation
double bearing = calculateBearing(37.7749, -122.4194,
                                  34.0522, -118.2437);
std::cout << "Bearing: " << bearing << " degrees\n";

// Destination point
GeoPoint dest = destinationPoint(37.7749, -122.4194, 
                                 180.0,      // South
                                 100000.0);  // 100 km
```

## Well Coordinates in Config

Wells can be specified in input CRS coordinates:

```ini
[GRID]
input_crs = EPSG:4326
model_crs = EPSG:32610
use_local_coordinates = true
local_origin_x = -122.4194
local_origin_y = 37.7749

[WELL1]
name = INJECTOR
type = INJECTOR
# Coordinates in WGS84 (will be transformed)
x = -122.4190
y = 37.7750
z = -500.0        # Depth below local origin datum
control_mode = RATE
target_value = 0.001
```

## Best Practices

1. **Use Appropriate CRS**: Choose a projected CRS (meters) for simulation accuracy
2. **Local Origin**: Always use a local origin for large domains to avoid precision loss
3. **Consistent Units**: Model CRS should use meters for compatibility with FSRM
4. **Document CRS**: Include CRS information in output files for reproducibility
5. **Validate**: Test coordinate transformations before running large simulations

## Limitations

- Binary Gmsh files with embedded coordinates require PROJ for transformation
- Some exotic CRS may not be supported by PROJ
- Real-time transformation of solution fields not yet implemented

## References

- [EPSG Registry](https://epsg.io/)
- [PROJ Documentation](https://proj.org/)
- [UTM Zone Calculator](https://mangomap.com/robertyoung/maps/69585/what-utm-zone-am-i-in-)

## See Also

- [Unstructured Meshes](UNSTRUCTURED_MESHES.md)
- [Configuration Reference](CONFIGURATION.md)
