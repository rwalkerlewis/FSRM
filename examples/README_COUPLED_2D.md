# 2D Coupled Reservoir Simulation Example

## Overview

This example demonstrates a fully integrated reservoir simulation combining:

- **Fluid Flow**: Two-phase (oil-water) flow with Darcy's law
- **Well Operations**: 1 water injector + 2 oil producers
- **Geomechanics**: Pressure-dependent subsidence and stress-strain coupling
- **Property Coupling**: Dynamic porosity and permeability from geomechanics

## Config-Driven Approach

All parameters are now specified in the configuration file:

```bash
./simulator -c config/coupled_reservoir_2d.config
```

See `config/coupled_reservoir_2d.config` for the complete parameter set.

## Features

### Advanced Visualization

The example generates **publication-quality 2×2 subplot visualizations** with:

- **Pressure field** (blue → cyan → green → yellow → red color map)
- **Water saturation** (Jet color map)
- **Surface subsidence** (white → yellow → orange → red)
- **Permeability changes** (Viridis color map)

Each plot includes:
- Well locations (green circles = injectors, red circles = producers)
- Borders and clean formatting
- Native C++ PPM generation → PNG conversion via ImageMagick

### Physics Models

#### Fluid Flow
- Darcy velocity: `q = -(k/μ) ∇P`
- Compressible single-phase equivalent
- Two-phase fractional flow with Corey relative permeabilities

#### Geomechanics
- Uniaxial compaction model
- Subsidence from pressure depletion: `ε_v = α ΔP / K`
- Lateral displacement from Poisson effects

#### Coupling
- Porosity update: `φ = φ₀(1 + ε_v)`
- Kozeny-Carman permeability: `k = k₀ (φ/φ₀)³ [(1-φ₀)/(1-φ)]²`

## Running the Example

```bash
cd build/examples
./ex_coupled_reservoir_2d
```

### Output

Generated in `output_coupled_2d/`:

- `coupled_sim_XXXX.png` - 1600×1200 pixel visualizations (21 timesteps)
- `info_XXXX.txt` - Detailed statistics and metadata

## Simulation Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| Domain | 500 m × 200 m | Reservoir dimensions |
| Grid | 150 × 60 cells | Spatial discretization |
| Duration | 1 year | Simulation time |
| Timesteps | 100 | Temporal discretization |
| Porosity | 0.2 (20%) | Initial porosity |
| Permeability | 100 mD | Initial permeability |
| Young's Modulus | 10 GPa | Rock stiffness |
| Poisson's Ratio | 0.25 | Lateral strain ratio |
| Initial Pressure | 30 MPa | Reservoir pressure |
| Biot Coefficient | 0.7 | Poroelastic coupling |

## Well Configuration

### INJ-1 (Injector)
- Location: (100, 100) m
- Function: Water injection
- Pressure: 35 MPa (5 MPa overpressure)

### PROD-1 (Producer)
- Location: (400, 50) m
- Function: Oil production (north)
- Pressure: 20 MPa (10 MPa drawdown)

### PROD-2 (Producer)
- Location: (400, 150) m
- Function: Oil production (south)
- Pressure: 20 MPa (10 MPa drawdown)

## Visualization Details

### Color Maps

The example implements several scientific color maps:

1. **PressureColorMap**: Blue → Cyan → Green → Yellow → Red
   - Ideal for pressure gradients
   - High contrast for flow patterns

2. **JetColorMap**: Blue → Cyan → Yellow → Red
   - Classic rainbow for saturation
   - Widely recognized in petroleum industry

3. **SubsidenceColorMap**: White → Yellow → Orange → Red → Dark Red
   - Highlights deformation severity
   - White indicates no subsidence

4. **ViridisColorMap**: Purple → Blue → Teal → Green → Yellow
   - Perceptually uniform
   - Colorblind-friendly

### Plot Layout

```
┌─────────────────┬─────────────────┐
│  Pressure (MPa) │  Water Sat (-)  │
│                 │                 │
│  [Pressure map] │  [Saturation]   │
│                 │                 │
├─────────────────┼─────────────────┤
│ Subsidence (m)  │ Permeability    │
│                 │     (mD)        │
│  [Subsidence]   │  [Perm change]  │
│                 │                 │
└─────────────────┴─────────────────┘
```

Each subplot: 700×500 pixels  
Total image: 1600×1200 pixels

## Technical Implementation

### Image Generation

Uses **native C++ PPM format** (Portable Pixmap):
- No external dependencies for image creation
- P6 binary format for efficiency
- Post-processing with ImageMagick for PNG conversion

### Color Class

```cpp
struct Color {
    unsigned char r, g, b;
    static Color lerp(const Color& c1, const Color& c2, double t);
};
```

Provides smooth color interpolation for gradients.

### PlotGenerator Class

High-level interface for creating complex visualizations:

```cpp
PlotGenerator plot(width, height, margins...);
plot.drawField(data, colormap, vmin, vmax);
plot.drawWell(x, y, is_injector);
plot.drawColorbar(colormap, vmin, vmax);
plot.writePPM(filename);
plot.convertToPNG(ppm_file, png_file);
```

## Expected Results

After 1 year of waterflooding:

- **Pressure**: 
  - Elevated near injector (~35 MPa)
  - Depleted near producers (~20 MPa)
  - Pressure communication across reservoir

- **Saturation**:
  - Water bank propagating from injector
  - Oil drainage toward producers
  - Mixing front with ~50% saturation

- **Subsidence**:
  - Maximum near producers (pressure depletion)
  - ~1-2 mm surface subsidence
  - Bowls of subsidence around each producer

- **Permeability**:
  - Slight reduction near producers (compaction)
  - Slight increase near injector (dilation)
  - Changes on order of ±0.5%

## Code Structure

```
ex_coupled_reservoir_2d.cpp
├── Color Maps (100 lines)
│   ├── JetColorMap
│   ├── ViridisColorMap
│   ├── PressureColorMap
│   └── SubsidenceColorMap
│
├── PlotGenerator Class (200 lines)
│   ├── Field rendering
│   ├── Well drawing
│   ├── Colorbar generation
│   └── PPM/PNG export
│
├── Physics Models (300 lines)
│   ├── CoupledFSRMulator
│   ├── solveFlow() - Darcy equation
│   ├── solveGeomechanics() - Subsidence
│   └── updatePoroPermeFromGeomechanics()
│
├── Visualization (200 lines)
│   └── generateCombinedPlot() - 2×2 layout
│
└── Main (100 lines)
    ├── Setup and initialization
    ├── Time stepping loop
    └── Output generation
```

Total: ~900 lines of well-documented C++

## Customization

### Modify Well Locations

```cpp
sim.addWell(Well(x, y, is_injector, rate, "NAME"));
```

### Change Domain Size

```cpp
const double Lx = 500.0;  // Length (m)
const double Ly = 200.0;  // Width (m)
```

### Adjust Time Duration

```cpp
const double total_time = 3600.0 * 24 * 365;  // 1 year in seconds
const int num_steps = 100;  // Number of timesteps
```

### Modify Rock Properties

```cpp
sim.setRockProperties(
    phi,    // Porosity (0-1)
    k,      // Permeability (m²)
    E,      // Young's modulus (Pa)
    nu      // Poisson's ratio (0-0.5)
);
```

## Limitations

This is a **demonstration example** with simplified physics:

- ⚠️ Explicit finite difference (can be unstable)
- ⚠️ Simple boundary conditions
- ⚠️ No adaptive timestepping
- ⚠️ Uniaxial compaction only (not full 3D stress)
- ⚠️ Simplified two-phase flow

For production simulations, use the full PETSc-based solvers with:
- Implicit time integration
- Proper boundary conditions
- Flux limiters
- Full tensor geomechanics

## Dependencies

### Required
- PETSc (for initialization, though not heavily used in this example)
- MPI (for parallel support)
- Standard C++17 compiler

### Optional
- ImageMagick (`convert` command) for PNG generation
- Otherwise outputs PPM files (viewable with most image viewers)

## Viewing Results

### Linux
```bash
eog output_coupled_2d/coupled_sim_0100.png
# or
xdg-open output_coupled_2d/coupled_sim_0100.png
```

### Convert all to animated GIF
```bash
convert -delay 20 -loop 0 output_coupled_2d/coupled_sim_*.png animation.gif
```

### Create video
```bash
ffmpeg -framerate 5 -pattern_type glob -i 'output_coupled_2d/coupled_sim_*.png' \
       -c:v libx264 -pix_fmt yuv420p reservoir_simulation.mp4
```

## References

- Darcy's Law for porous media flow
- Terzaghi's consolidation theory
- Biot's poroelasticity
- Corey relative permeability correlations
- Kozeny-Carman equation for permeability

## Author Notes

This example showcases:
- ✅ Beautiful scientific visualization in pure C++
- ✅ Multi-physics coupling
- ✅ Realistic reservoir geometry with wells
- ✅ Publication-ready plots without Python
- ✅ Self-contained (no external plotting libraries)

Perfect for demonstrating FSRM capabilities!
